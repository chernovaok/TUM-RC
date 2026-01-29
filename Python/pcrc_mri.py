"""
PCRC-MRI  
TUM, December 2025
"""

import csv
import sys
import os
import math
from typing import Optional

def calculate_pcrcmri_risk(
    age: float,
    psa: float,
    dre: str,
    priorbiopsy: str,
    race: str,
    volume: Optional[float],
    pirads: int
) -> Optional[float]:
    # 1. Handle Missing Data
    if volume is None or not priorbiopsy or str(priorbiopsy).lower() == 'nan':
        return None

    # 2. Variable Normalization
    dre_abnormal = 1 if str(dre).lower().strip() == "abnormal" else 0
    priorbiopsy_yes = 1 if str(priorbiopsy).lower().strip() == "yes" else 0
    
    # 3. Parameter Mapping
    race_coeff_map = {
        "African American": 0.23780,
        "Asian": -0.77190,
        "Caucasian": 0.0,
        "White": 0.0,
        "Other": 0.18840,
        "Unknown": -0.08560
    }
    
    pirads_coeff_map = {
        2: 0.0,
        3: 0.51020,
        4: 1.37510,
        5: 2.76270
    }

    # Coefficients
    intercept = -4.62990
    coef_age = 0.05470
    coef_psa = 0.01920
    coef_dre = 0.86460
    coef_prev_biopsy = -0.61630
    coef_prostate_volume = -0.01790
    coef_psa_density_binary = 0.88020

    # 4. Logic Calculations
    psad_value = psa / volume
    binary_psad = 1 if psad_value >= 0.15 else 0

    # Map race (Defaults to White/Caucasian)
    race_coef = race_coeff_map.get(race.strip(), race_coeff_map["White"])
    pirads_coef = pirads_coeff_map.get(int(pirads), 0.0)

    # 5. Calculate Log-Odds
    log_odds = (
        intercept +
        (coef_age * age) +
        race_coef +
        (coef_psa * psa) +
        (coef_dre * dre_abnormal) +
        (coef_prev_biopsy * priorbiopsy_yes) +
        (coef_prostate_volume * volume) +
        (coef_psa_density_binary * binary_psad) +
        pirads_coef
    )

    # 6. Convert Log-Odds to Probability
    risk = math.exp(log_odds) / (1 + math.exp(log_odds))
    return round(risk * 100, 2)

def main(input_file: str, output_file: str):
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    results = []
    with open(input_file, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames + ["PCRC_MRI"]
        
        for row in reader:
            try:
                # Handle volume safely
                vol_raw = row.get('volume')
                vol_val = float(vol_raw) if vol_raw and vol_raw.lower() != 'na' else None
                
                # Check for race column, otherwise default to "White"
                row_race = row.get('race', 'White')
                if not row_race: 
                    row_race = 'White'

                risk = calculate_pcrcmri_risk(
                    age=float(row['age']),
                    psa=float(row['psa']),
                    dre=row['dre'],
                    priorbiopsy=row['priornegbiopsy'],
                    race=row_race,
                    volume=vol_val,
                    pirads=int(row['pirads'])
                )
                
                row["PCRC_MRI"] = f"{risk}" if risk is not None else "NA"
                
            except:
                row["PCRC_MRI"] = "NA"
            
            results.append(row)

    with open(output_file, mode='w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Process complete. Results saved to {output_file}")

if __name__ == "__main__":
    input_csv = sys.argv[1] if len(sys.argv) > 1 else "pcrcmri_input.csv"
    output_csv = sys.argv[2] if len(sys.argv) > 2 else "pcrcmri_output.csv"
    main(input_csv, output_csv)
