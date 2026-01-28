"""
ERSPC34-RC Calculator 
TUM, December 2025
"""

import csv
import os
import json
import sys  # Required for command line arguments
import requests
from abc import ABC, abstractmethod
from typing import Literal, Optional, List

# Requirement: pip install pydantic requests
from pydantic import BaseModel, field_validator, confloat

# --- 1. Base Interfaces ---

class BaseRequestService(ABC):
    @abstractmethod
    def get_app_url(self) -> str: pass
    @abstractmethod
    def build_payload(self, data: dict) -> dict: pass
    @abstractmethod
    def parse_output(self, response_data: dict) -> Optional[float]: pass

class BaseClient(ABC):
    @abstractmethod
    def evaluate(self, inputs: List[dict]) -> List[Optional[float]]: pass

# --- 2. Model Logic (Validation) ---

class ERSPC34_RC_Model(BaseModel):
    age: confloat(ge=50, le=75)
    psa: confloat(ge=0.4, le=50)
    volume: confloat(ge=10, le=110)
    dre: int = 0 
    pirads: Literal[0, 1, 2, 3, 4] = 0
    priornegbiopsy: int = 1 

    model_config = {"validate_default": True}

    @field_validator("dre", mode="before")
    @classmethod
    def normalize_dre(cls, v):
        if isinstance(v, str):
            mapping = {"normal": 0, "abnormal": 1}
            return mapping.get(v.lower().strip(), 0)
        return v

    @field_validator("pirads", mode="before")
    @classmethod
    def normalize_pirads(cls, v):
        try: return int(v) - 1
        except: return 0

    @field_validator("priornegbiopsy", mode="before")
    @classmethod
    def normalize_priornegbiopsy(cls, v):
        if isinstance(v, str):
            mapping = {"yes": 0, "no": 1}
            return mapping.get(v.lower().strip(), 1)
        return v

# --- 3. Service Implementation ---

class ERSPC34_RC_Service(BaseRequestService):
    def get_app_url(self) -> str:
        return "https://calculatorallrpcrc.prostatecancer-riskcalculator.com/result.php"

    def build_payload(self, data: dict) -> dict:
        try:
            validated = ERSPC34_RC_Model(**data)
            return {
                "age_mri": str(validated.age),
                "dre_result": validated.dre,
                "mri": 0, "mri_dre": 0, "phi": 1,
                "pirads": validated.pirads,
                "previous_biopsy": validated.priornegbiopsy,
                "prostate_volume": str(validated.volume),
                "prostate_volume_range": 0,
                "psa": str(validated.psa),
            }
        except Exception:
            return {}

    def parse_output(self, response_data: dict) -> Optional[float]:
        try:
            if "result" in response_data:
                risk = response_data["result"].get("significant_cancer_risk")
                if risk is not None:
                    return round(float(risk) * 100, 2)
        except: return None
        return None

# --- 4. Client Implementation ---

class RequestClient(BaseClient):
    def __init__(self, service: BaseRequestService):
        self.service = service

    def evaluate(self, inputs: List[dict]) -> List[Optional[float]]:
        url = self.service.get_app_url()
        headers = {
            "User-Agent": (
                "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/143.0.0.0 Safari/537.36"
            ),
            "Accept": "*/*",
            "Accept-Language": "en-US,en;q=0.9",
            "Referer": "https://calculatorallrpcrc.prostatecancer-riskcalculator.com/?calc_3_default_lang=en",
        }
        
        values = []
        for item in inputs:
            payload = self.service.build_payload(item)
            if not payload:
                values.append(None)
                continue
            
            try:
                response = requests.post(url, data=payload, headers=headers, timeout=15)
                if response.status_code == 200:
                    data = response.json()
                    prob = self.service.parse_output(data)
                    values.append(prob)
                else:
                    values.append(None)
            except Exception:
                values.append(None)
        return values

# --- 5. Main Execution Logic ---

def main(input_csv: str, output_csv: str):
    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found.")
        return

    service = ERSPC34_RC_Service()
    client = RequestClient(service)

    with open(input_csv, mode='r', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    
    if not rows:
        print("No data found in the input file.")
        return

    results = client.evaluate(rows)

    for row, risk in zip(rows, results):
        row["ERSPC34-RC"] = f"{risk}" if risk is not None else "NA"

    with open(output_csv, mode='w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Success! Processed {len(rows)} rows into {output_csv}")

if __name__ == "__main__":
    # Use command line arguments if provided, else use defaults
    input_file = sys.argv[1] if len(sys.argv) > 1 else "data_input.csv"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "erspc_output.csv"
    main(input_file, output_file)
