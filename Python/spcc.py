"""
SPCC  
TUM, December 2025
"""

import csv
import sys
import os
import math
from typing import Optional

def calculate_spcc(
    age: float,
    psa: float,
    volume: float,
    history: str,
    race: str,
    pirads: int
) -> Optional[float]:
    """
    Python implementation of the Stanford Prostate Cancer Risk Model.
    """
    # 1. Handle Missing Data or PIRADS out of range (Model only valid for 3-5)
    if volume == 0 or pirads < 3 or pirads > 5:
        return None

    # 2. Variable Normalization
    history_clean = str(history).strip()
    race_clean = str(race).strip()
    
    # 3. Parameter Mapping
    history_coeff_map = {
        "Biopsy Naive": 0.0,
        "Active Surveillance": 0.139,
        "Prior Negative Biopsy": -0.6112
    }
    
    pirads_coeff_map = {
        3: 0.0,
        4: 1.051,
        5: 1.717
    }
    
    race_coeff_map = {
        "Asian": -0.1757,
        "Black": -0.145,
        "Hispanic": -0.288,
        "White": 0.4736,
        "Other or Unknown": 0.0
    }

    # Fixed Constants
    intercept = -4.326
    psa_d_coeff = 7.245
    age_coeff = 3.099

    # 4. Extract specific coefficients
    # Defaulting to 0.0 (Naive/Other) if mapping fails to prevent crashes
    ahx = history_coeff_map.get(history_clean, 0.0)
    brace = race_coeff_map.get(race_clean, 0.4736) # Defaults to White if not found
    b_pirads = pirads_coeff_map.get(int(pirads))

    if b_pirads is None:
        return None

    # 5. Logic Calculations
    psa_d = psa / volume
    
    # 6. Calculate Log-Odds
    # Note: Age is divided by 100 as per R code (const_bAge * age / 100)
    log_odds = (
        intercept +
        ahx +
        (psa_d * psa_d_coeff) +
        brace +
        (age_coeff * age / 100) +
        b_pirads
    )

    # 7. Convert Log-Odds to Probability
    risk = math.exp(log_odds) / (1 + math.exp(log_odds))
    return round(risk * 100, 2)

def main(input_file: str, output_file: str):
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    results = []
    with open(input_file, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames + ["SPCC"]
        
        for row in reader:
            try:
                # Handle Missing Race: Assume "White" if column missing or empty
                row_race = row.get('race', 'White')
                if not row_race or str(row_race).lower() == 'nan':
                    row_race = 'White'
                # Map priornegbiopsy to history
                prior_neg = str(row.get('priornegbiopsy', '')).lower().strip()
                history_val = "Prior Negative Biopsy" if prior_neg == "yes" else "Biopsy Naive"

                risk = calculate_spcc(
                    age=float(row['age']),
                    psa=float(row['psa']),
                    volume=float(row['volume']),
                    history=history_val,
                    race=row_race,
                    pirads=int(float(row['pirads'])) # float wrap handles '4.0' strings
                )
                
                row["SPCC"] = f"{risk}" if risk is not None else "NA"
                
            except:
                row["SPCC"] = "NA"
            
            results.append(row)

    with open(output_file, mode='w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Success! SPCC processed. Results in {output_file}")

if __name__ == "__main__":
    input_csv = sys.argv[1] if len(sys.argv) > 1 else "spcc_input.csv"
    output_csv = sys.argv[2] if len(sys.argv) > 2 else "spcc_output.csv"
    main(input_csv, output_csv)