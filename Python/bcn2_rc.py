"""
BCN2RC Simple Web Scraper.

A simple, sequential (single-browser, no parallelism) scraper for the BCN2RC
calculator at https://mripcaprediction.shinyapps.io/MRIPCaPrediction/

Uses generous wait times to avoid Shiny timing issues.
Extracts probability from result images via OCR (pytesseract).

Required columns:
    - age or Age: Patient age (numeric)
    - psa or PSA: PSA level (ng/mL)
    - volume or PV or prostate_volume: Prostate volume (cc)
    - pirads or PIRADS: PI-RADS score (1-5)
    - dre or DRE: Digital rectal exam ("Normal" or "Suspicious")
    - fhx or family_history: Family history ("Yes" or "No")
    - priornegbiopsy or prior_neg_bx: Prior negative biopsy ("Yes" or "No")

Outputs:
    - bcn2rc_cspca_risk: Probability of clinically significant cancer (0-1)

Usage:
    python scraper_simple.py input.csv output.csv
"""

import sys
import os
import time
import re
import base64
import io
import pandas as pd
from PIL import Image
import pytesseract

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


# ---------------------------------------------------------------------------
# Inline utilities (so this script is fully self-contained, no extra files)
# ---------------------------------------------------------------------------

class InvalidInputError(ValueError):
    """Raised when input validation fails."""
    pass


def parse_yes_no(val, default="No"):
    """Parse a yes/no value to 'Yes' or 'No'."""
    if val is None or str(val).strip() == "":
        return default
    val_lower = str(val).strip().lower()
    if val_lower in ("yes", "1", "true", "y"):
        return "Yes"
    return "No"


def parse_dre(val, default="Normal"):
    """Parse a DRE value to 'Normal' or 'Abnormal'."""
    if val is None or str(val).strip() == "":
        return default
    val_lower = str(val).strip().lower()
    if val_lower in ("abnormal", "suspicious", "yes", "1", "true"):
        return "Abnormal"
    return "Normal"


# Configuration
BCN2RC_URL = "https://mripcaprediction.shinyapps.io/MRIPCaPrediction/"
WAIT_AFTER_INPUT = 3      # seconds to wait after setting inputs
RETRY_WAIT = 5            # seconds to wait before retry on OCR failure
POLL_INTERVAL = 1         # seconds between polls when waiting for image
POLL_TIMEOUT = 15         # max seconds to poll for a valid image


def make_driver():
    """Create a headless Chrome driver (no performance logs needed)."""
    options = webdriver.ChromeOptions()
    options.add_argument("--headless=new")
    options.add_argument("--disable-gpu")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--window-size=1920,1080")
    return webdriver.Chrome(options=options)


def validate_row(row, idx):
    """Validate input row data. Raises InvalidInputError on bad data."""
    errors = []

    age = row.get("age", row.get("Age", None))
    psa = row.get("psa", row.get("PSA", None))
    volume = row.get("volume", row.get("PV", row.get("prostate_volume", None)))
    pirads = row.get("pirads", row.get("PIRADS", None))

    if age is None:
        errors.append("missing age")
    else:
        try:
            float(age)
        except (ValueError, TypeError):
            errors.append(f"age={age} (not numeric)")

    if psa is None:
        errors.append("missing psa")
    else:
        try:
            float(psa)
        except (ValueError, TypeError):
            errors.append(f"psa={psa} (not numeric)")

    if volume is None:
        errors.append("missing volume")
    else:
        try:
            float(volume)
        except (ValueError, TypeError):
            errors.append(f"volume={volume} (not numeric)")

    if pirads is None:
        errors.append("missing pirads")
    else:
        try:
            p = int(pirads)
            if p not in [1, 2, 3, 4, 5]:
                errors.append(f"pirads={p} (expected 1-5)")
        except (ValueError, TypeError):
            errors.append(f"pirads={pirads} (invalid)")

    if errors:
        raise InvalidInputError(f"Row {idx+1}: Invalid values - {', '.join(errors)}")


def _extract_row_inputs(row):
    """
    Extract and normalize input values from a row dict.

    Returns a dict with keys: age, psa, volume, dre, fhx, prior, pirads
    """
    age = float(row.get("age", row.get("Age", 68)))
    psa = float(row.get("psa", row.get("PSA", 6)))
    volume = float(row.get("volume", row.get("PV", row.get("prostate_volume", 50))))
    pirads = str(int(row.get("pirads", row.get("PIRADS", 3))))

    dre_raw = row.get("dre", row.get("DRE", "Normal"))
    dre = "Suspicious" if parse_dre(dre_raw) != "Normal" else "Normal"

    fhx_raw = row.get("fhx", row.get("family_history", row.get("Family_history", "No")))
    fhx = parse_yes_no(fhx_raw, "No")

    prior_raw = row.get("prior_neg_bx", row.get("NBiopsy", row.get("priornegbiopsy", "First")))
    prior_str = str(prior_raw).strip().lower()
    if prior_str in ("first", "no", "none", "no prior biopsy", "0", "false"):
        prior = "First"
    elif prior_str in ("negative", "previous negative", "yes", "1", "true"):
        prior = "Negative"
    else:
        prior = "First"

    return {
        "age": age,
        "psa": psa,
        "volume": volume,
        "dre": dre,
        "fhx": fhx,
        "prior": prior,
        "pirads": pirads,
    }


def get_base64_image_src(driver):
    """Find the displayed base64 PNG image on the page (the BCN2RC result plot)."""
    try:
        images = driver.find_elements(By.TAG_NAME, "img")
        for img in images:
            if not img.is_displayed():
                continue
            src = img.get_attribute("src")
            if src and src.startswith("data:image/png;base64,"):
                return src
    except Exception:
        pass
    return None


def parse_probability_from_ocr_text(text):
    """
    Parse CsPCa probability from OCR text.

    Looks for patterns like:
        'Probability of CsPCa= 66.14 %'
        'CsPCa=7.5%'
        '42.3 %'

    Returns the numeric value (e.g. 66.14) or None.
    """
    if not text:
        return None

    # Primary pattern: "Probability of CsPCa= XX.XX %"
    match = re.search(
        r'Probability\s+of\s+CsPCa\s*=\s*(\d+\.?\d*)\s*%', text, re.IGNORECASE
    )
    if match:
        return float(match.group(1))

    # Fallback: "CsPCa= XX.XX %"
    match = re.search(r'CsPCa\s*=\s*(\d+\.?\d*)\s*%', text, re.IGNORECASE)
    if match:
        return float(match.group(1))

    # Last resort: any "XX.XX %"
    match = re.search(r'(\d+\.?\d*)\s*%', text)
    if match:
        return float(match.group(1))

    return None


def extract_probability_from_base64(base64_src):
    """Extract CsPCa probability from a base64-encoded PNG image via OCR."""
    if not base64_src or not base64_src.startswith("data:image/png;base64,"):
        return None

    try:
        base64_data = base64_src.split(",")[1]
        img_data = base64.b64decode(base64_data)
        image = Image.open(io.BytesIO(img_data))

        # Crop top ~15% where the probability text is
        width, height = image.size
        top_region = image.crop((0, 0, width, int(height * 0.15)))

        text = pytesseract.image_to_string(top_region)
        prob = parse_probability_from_ocr_text(text)
        if prob is not None:
            return prob

        # If cropped region failed, try full image
        full_text = pytesseract.image_to_string(image)
        return parse_probability_from_ocr_text(full_text)

    except Exception:
        return None


def set_inputs(driver, values):
    """
    Set BCN2RC form inputs and wait for Shiny to recompute.

    Args:
        driver: Selenium WebDriver
        values: dict from _extract_row_inputs()
    """
    # Set numeric fields via clear + send_keys
    for input_id, value in [("Age", values["age"]), ("PSA", values["psa"]), ("PV", values["volume"])]:
        element = driver.find_element(By.ID, input_id)
        element.click()
        element.clear()
        element.send_keys(str(value))

    # Set selectize dropdowns via JavaScript
    driver.execute_script(f"""
        function setSelectize(id, value) {{
            var $sel = $('#' + id);
            if ($sel.length && $sel[0].selectize) {{
                $sel[0].selectize.setValue(value, false);
            }}
        }}
        setSelectize('PIRADS', '{values["pirads"]}');
        setSelectize('DRE', '{values["dre"]}');
        setSelectize('Family_history', '{values["fhx"]}');
        setSelectize('NBiopsy', '{values["prior"]}');
    """)

    # Click body to blur inputs and trigger Shiny reactive update
    driver.find_element(By.TAG_NAME, "body").click()

    # Wait for Shiny to recompute and render the new image
    time.sleep(WAIT_AFTER_INPUT)


def extract_result(driver):
    """
    Poll for a valid probability from the displayed image.

    Returns:
        (probability_0_to_1, None) on success
        (None, failure_reason) on failure
    """
    start = time.time()

    while time.time() - start < POLL_TIMEOUT:
        base64_src = get_base64_image_src(driver)
        if base64_src:
            prob = extract_probability_from_base64(base64_src)
            if prob is not None:
                return round(prob / 100, 3), None
        time.sleep(POLL_INTERVAL)

    # Final attempt
    base64_src = get_base64_image_src(driver)
    if not base64_src:
        return None, "NO_IMAGE"
    prob = extract_probability_from_base64(base64_src)
    if prob is not None:
        return round(prob / 100, 3), None
    return None, "OCR_FAILED"


def process_row(driver, row, idx):
    """
    Process a single patient row.

    Returns:
        (probability_0_to_1, failure_reason_or_None)
    """
    values = _extract_row_inputs(row)
    set_inputs(driver, values)

    prob, reason = extract_result(driver)
    if prob is not None:
        return prob, None

    # Retry once: wait longer and try again
    print(f"    Retry row {idx+1} after {RETRY_WAIT}s...")
    time.sleep(RETRY_WAIT)
    prob, reason = extract_result(driver)
    return prob, reason


def run_batch(input_csv, output_csv, num_workers=1):
    """
    Process a CSV file sequentially with a single browser.

    Args:
        input_csv: Path to input CSV file.
        output_csv: Path to output CSV file.
        num_workers: Ignored (always uses 1 browser). Kept for interface compatibility.

    Returns:
        DataFrame with bcn2rc_cspca_risk column added.
    """
    df = pd.read_csv(input_csv)
    num_rows = len(df)

    print(f"BCN2RC (simple): Processing {num_rows} rows sequentially...")
    print(f"  URL: {BCN2RC_URL}")
    print(f"  Wait per row: ~{WAIT_AFTER_INPUT}s + polling up to {POLL_TIMEOUT}s\n")

    driver = make_driver()
    results = [None] * num_rows
    start_time = time.time()

    try:
        # Load page
        print("  Loading page...")
        driver.get(BCN2RC_URL)

        # Wait for BCN2RC tab to be clickable
        WebDriverWait(driver, 30).until(
            EC.element_to_be_clickable((By.XPATH, "//a[contains(text(), 'BCN2RC')]"))
        )

        # Click BCN2RC tab
        bcn2rc_tab = driver.find_element(By.XPATH, "//a[contains(text(), 'BCN2RC')]")
        bcn2rc_tab.click()
        time.sleep(3)

        # Wait for initial image to render
        print("  Waiting for initial result image...")
        for _ in range(30):
            src = get_base64_image_src(driver)
            if src and extract_probability_from_base64(src) is not None:
                break
            time.sleep(0.5)
        print("  Browser ready.\n")

        # Process each row sequentially
        for i, row in df.iterrows():
            try:
                validate_row(row, i)
                prob, reason = process_row(driver, row, i)

                if prob is not None:
                    print(f"  Row {i+1}/{num_rows}: CsPCa={prob:.3f}")
                    results[i] = (prob, None)
                else:
                    print(f"  Row {i+1}/{num_rows}: FAILED - {reason}")
                    results[i] = (None, reason)

            except InvalidInputError as e:
                print(f"  Row {i+1}/{num_rows}: INVALID - {e}")
                results[i] = (None, str(e))
            except Exception as e:
                print(f"  Row {i+1}/{num_rows}: ERROR - {e}")
                results[i] = (None, str(e))

        elapsed = time.time() - start_time

        # Write output
        df["bcn2rc_cspca_risk"] = [r[0] if r else None for r in results]
        df.to_csv(output_csv, index=False)

        # Summary
        errors = sum(1 for r in results if r is None or r[0] is None)
        print(f"\n{'='*50}")
        print(f"Results saved to {output_csv}")
        print(f"Successful: {num_rows - errors}/{num_rows}")
        print(f"Time: {elapsed:.1f}s ({elapsed/num_rows:.2f}s per row)")
        if errors > 0:
            print(f"Errors: {errors}")

    finally:
        print("Closing browser...")
        driver.quit()

    return df


if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else "bcn2rc_test_input.csv"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "bcn2rc_output.csv"
    run_batch(input_file, output_file)
