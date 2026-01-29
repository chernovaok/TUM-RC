"""
MSP-RC Calculator - Parallel Browser
TUM, December 2025
"""

import time
import re
import sys
import json
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue, Empty
from threading import Lock

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.common.exceptions import WebDriverException

MSP_RC_URL = "https://darasriskcalcs.shinyapps.io/MSP-RC/"
NUM_WORKERS = 5

# Valid values for validation
VALID_DRE = {"Normal", "Abnormal"}
VALID_YES_NO = {"Yes", "No"}
VALID_PIRADS = {"1", "2", "3", "4", "5"}
VALID_PREV_BX = {"Negative", "ASAP/HGPIN", "GG1 cancer", "No prior biopsy"}


class InvalidInputError(ValueError):
    pass


def normalize_value(val, valid_set):
    """Normalize a value to match the valid set (case-insensitive)."""
    if val is None:
        return None
    val_lower = str(val).strip().lower()
    for valid in valid_set:
        if valid.lower() == val_lower:
            return valid
    return str(val)  # Return original if no match


def validate_row(row, row_idx):
    """Validate row values, raise error if invalid."""
    errors = []
    
    dre = row.get("dre", row.get("DRE", None))
    if dre is not None:
        dre_normalized = normalize_value(dre, VALID_DRE)
        if dre_normalized not in VALID_DRE:
            errors.append(f"dre='{dre}' (expected: {VALID_DRE})")
    
    # aa defaults to "No" if missing
    aa = row.get("aa", row.get("AA", None))
    if aa is not None:
        aa_normalized = normalize_value(aa, VALID_YES_NO)
        if aa_normalized not in VALID_YES_NO:
            errors.append(f"aa='{aa}' (expected: {VALID_YES_NO})")
    
    # fhx defaults to "No" if missing
    fhx = row.get("fhx", row.get("family_history", None))
    if fhx is not None:
        fhx_normalized = normalize_value(fhx, VALID_YES_NO)
        if fhx_normalized not in VALID_YES_NO:
            errors.append(f"fhx='{fhx}' (expected: {VALID_YES_NO})")
    
    pirads = str(row.get("pirads", row.get("PIRADS", "")))
    if pirads and pirads not in VALID_PIRADS:
        errors.append(f"pirads='{pirads}' (expected: {VALID_PIRADS})")
    
    # prev_bx: also accept priornegbiopsy column
    prev_bx = row.get("prev_bx", row.get("previous_biopsy", None))
    priornegbiopsy = row.get("priornegbiopsy", None)
    
    # Only validate prev_bx if it's provided (not priornegbiopsy)
    if prev_bx is not None and priornegbiopsy is None:
        prev_bx_normalized = normalize_value(prev_bx, VALID_PREV_BX)
        if prev_bx_normalized not in VALID_PREV_BX:
            errors.append(f"prev_bx='{prev_bx}' (expected: {VALID_PREV_BX})")
    
    if errors:
        raise InvalidInputError(f"Row {row_idx + 1}: Invalid values - " + "; ".join(errors))


def make_driver():
    """Create Chrome driver with performance logging for WebSocket monitoring."""
    options = webdriver.ChromeOptions()
    options.add_argument("--headless=new")
    options.add_argument("--disable-gpu")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.set_capability("goog:loggingPrefs", {"performance": "ALL"})
    return webdriver.Chrome(options=options)


def get_result_text(driver):
    """Get the current result percentages as text tuple."""
    text = driver.find_element(By.TAG_NAME, "body").text
    pos = re.search(r"risk of a positive biopsy is:\s*[\n\r]*\s*(\d+(?:\.\d+)?)\s*%", text, re.I)
    cs = re.search(r"risk of clinically significant.*?:\s*[\n\r]*\s*(\d+(?:\.\d+)?)\s*%", text, re.I)
    return (pos.group(1) if pos else None, cs.group(1) if cs else None)


def clear_logs(driver):
    """Clear existing performance logs."""
    try:
        driver.get_log("performance")
    except:
        pass


def check_for_ack(driver):
    """Check if ACK received in performance logs."""
    try:
        logs = driver.get_log("performance")
        for log in logs:
            try:
                msg = json.loads(log["message"])
                if msg.get("message", {}).get("method") == "Network.webSocketFrameReceived":
                    payload = msg["message"]["params"]["response"]["payloadData"]
                    if "ACK" in payload:
                        return True
            except:
                continue
    except:
        pass
    return False


def wait_for_result_or_ack(driver, old_result, timeout=10):
    """
    Wait until result text changes from old_result.
    Uses ACK as a signal that processing happened, but always waits for DOM update.
    Returns the new result tuple.
    """
    start = time.time()
    ack_received = False
    ack_time = None
    
    while time.time() - start < timeout:
        # Check if result changed first (primary condition)
        new_result = get_result_text(driver)
        if new_result != old_result and new_result[0] is not None:
            return new_result
        
        # Check for ACK (secondary signal that server processed)
        if not ack_received and check_for_ack(driver):
            ack_received = True
            ack_time = time.time()
        
        # If ACK received, wait up to 500ms more for DOM to update
        if ack_received and (time.time() - ack_time) > 0.5:
            # ACK received but result still hasn't changed after 500ms
            # This means the calculation produced the same result (rare but possible)
            return get_result_text(driver)
        
        time.sleep(0.02)  # 20ms poll interval
    
    # Timeout - return whatever we have
    return get_result_text(driver)


def set_inputs(driver, age, dre, aa, fhx, psa, volume, pirads, prev_bx):
    """Set all inputs atomically via JavaScript."""
    driver.execute_script(f"""
        Shiny.setInputValue('age', {age});
        Shiny.setInputValue('dre', '{dre}');
        Shiny.setInputValue('aarace', '{aa}');
        Shiny.setInputValue('familyhx', '{fhx}');
        Shiny.setInputValue('psa', {psa});
        Shiny.setInputValue('mri_vol', {volume});
        Shiny.setInputValue('mri_pirads', '{pirads}');
        Shiny.setInputValue('px_bx', '{prev_bx}');
    """)


def process_row(driver, row, row_idx):
    """Process a single row and return (pos_risk, cs_risk)."""
    # Get current result before changing inputs
    old_result = get_result_text(driver)
    
    # Clear logs before setting inputs
    clear_logs(driver)
    
    # Convert values
    age = row.get("age", row.get("patient_age", 55))
    
    dre_value = str(row.get("dre", row.get("DRE", "Normal"))).strip().lower()
    dre = "1" if dre_value == "abnormal" else "0"
    
    # aa defaults to "No" (baseline/lowest risk) if missing
    aa_value = str(row.get("aa", row.get("AA", "No"))).strip().lower()
    aa = "1" if aa_value == "yes" else "0"
    
    # fhx defaults to "No" (baseline/lowest risk) if missing
    fhx_value = str(row.get("fhx", row.get("family_history", "No"))).strip().lower()
    fhx = "1" if fhx_value == "yes" else "0"
    
    psa = row.get("psa", row.get("PSA", 3))
    volume = row.get("volume", row.get("volume_cc", 45))
    pirads = str(row.get("pirads", row.get("PIRADS", 3)))
    
    # Handle prev_bx: accept both prev_bx column and priornegbiopsy column
    prev_bx_value = row.get("prev_bx", row.get("previous_biopsy", None))
    priornegbiopsy = row.get("priornegbiopsy", None)
    
    if priornegbiopsy is not None:
        # Map priornegbiopsy: "yes" = had negative biopsy, "no" = no prior biopsy
        priorneg_str = str(priornegbiopsy).strip().lower()
        if priorneg_str in ("yes", "1", "true"):
            prev_bx_value = "Negative"
        else:
            prev_bx_value = "No prior biopsy"
    elif prev_bx_value is None:
        # Default to "No prior biopsy" if neither column provided
        prev_bx_value = "No prior biopsy"
    
    prev_bx_map = {"Negative": "0", "ASAP/HGPIN": "1", "GG1 cancer": "2", "No prior biopsy": "3"}
    prev_bx = prev_bx_map.get(prev_bx_value, "3")
    
    # Set inputs
    set_inputs(driver, age, dre, aa, fhx, psa, volume, pirads, prev_bx)
    
    # Wait for result change OR ACK
    new_result = wait_for_result_or_ack(driver, old_result)
    
    # Convert to floats and round to 3 decimal places
    pos = round(float(new_result[0]) / 100, 3) if new_result[0] else None
    cs = round(float(new_result[1]) / 100, 3) if new_result[1] else None
    
    return pos, cs


class BrowserPool:
    """Pool of reusable browser instances for parallel processing."""
    
    def __init__(self, num_browsers):
        self.browsers = Queue()
        self.num_browsers = num_browsers
        self.lock = Lock()
    
    def initialize(self):
        """Initialize all browsers and load the page."""
        print(f"Initializing {self.num_browsers} browsers...")
        for i in range(self.num_browsers):
            driver = make_driver()
            driver.get(MSP_RC_URL)
            # Wait for initial load
            WebDriverWait(driver, 20).until(
                lambda d: "risk of a positive biopsy" in d.find_element(By.TAG_NAME, "body").text
            )
            # Wait for initial result to be calculated (not None)
            for _ in range(100):  # up to 10 seconds
                result = get_result_text(driver)
                if result[0] is not None:
                    break
                time.sleep(0.1)
            self.browsers.put(driver)
            print(f"  Browser {i+1}/{self.num_browsers} ready")
        print("All browsers initialized!\n")
    
    def get_browser(self, timeout=120):
        """Get an available browser from the pool."""
        try:
            return self.browsers.get(timeout=timeout)
        except Empty:
            raise RuntimeError(f"Timed out waiting for browser after {timeout}s")
    
    def return_browser(self, driver):
        """Return a browser to the pool."""
        self.browsers.put(driver)
    
    def close_all(self):
        """Close all browsers."""
        while not self.browsers.empty():
            try:
                driver = self.browsers.get_nowait()
                driver.quit()
            except:
                pass


def process_row_with_pool(args):
    """Process a single row using the browser pool."""
    pool, idx, row, num_rows = args
    driver = None
    
    try:
        validate_row(row, idx)
        driver = pool.get_browser()
        pos, cs = process_row(driver, row, idx)
        
        if pos is not None and cs is not None:
            print(f"  Row {idx+1}/{num_rows}: Pos={pos:.3f}, CS={cs:.3f}")
            return idx, pos, cs, None
        else:
            print(f"  Row {idx+1}/{num_rows}: FAILED - Could not extract results")
            return idx, None, None, "Could not extract results"
            
    except InvalidInputError as e:
        print(f"  Row {idx+1}/{num_rows}: INVALID - {e}")
        return idx, None, None, str(e)
    except Exception as e:
        print(f"  Row {idx+1}/{num_rows}: ERROR - {e}")
        return idx, None, None, str(e)
    finally:
        if driver is not None:
            pool.return_browser(driver)


def run_batch(input_csv, output_csv, num_workers=NUM_WORKERS):
    """Process CSV file with parallel browsers, no reloads."""
    df = pd.read_csv(input_csv)
    num_rows = len(df)
    
    print(f"Processing {num_rows} rows with {num_workers} parallel browsers...")
    
    pool = BrowserPool(num_workers)
    pool.initialize()
    
    results = [None] * num_rows
    start_time = time.time()
    
    try:
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            tasks = [(pool, i, row, num_rows) for i, row in df.iterrows()]
            futures = {executor.submit(process_row_with_pool, task): task[1] for task in tasks}
            
            for future in as_completed(futures):
                idx, pos, cs, error = future.result()
                results[idx] = (pos, cs, error)
        
        elapsed = time.time() - start_time
        
        # Add results to dataframe
        df["positive_biopsy_risk"] = [r[0] for r in results]
        df["cspca_risk"] = [r[1] for r in results]
        
        df.to_csv(output_csv, index=False)
        
        # Summary
        errors = sum(1 for r in results if r[0] is None)
        print(f"\n{'='*50}")
        print(f"Results saved to {output_csv}")
        print(f"Successful: {num_rows - errors}/{num_rows}")
        print(f"Time: {elapsed:.1f}s ({elapsed/num_rows:.2f}s per row)")
        if errors > 0:
            print(f"Errors: {errors}")
        
    finally:
        print("Closing browsers...")
        pool.close_all()
    
    return df


if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else "msp_rc_input.csv"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "msp_rc_output.csv"
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else NUM_WORKERS
    run_batch(input_file, output_file, workers)
