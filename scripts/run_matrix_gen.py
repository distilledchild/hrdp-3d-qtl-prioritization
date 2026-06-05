import os
import shutil
import subprocess
import time
from datetime import datetime


samples = ["592BB", "607", "74AA", "A2DB", "D765A", "DA08A", "DA21A", "DA68A", "DBA9A", "DE8BA"]
dropbox_input_template = "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/hic/2023A/hic30_w_sb_options/{sample}/inter_30.hic"
dropbox_output_dir = "/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/hic2/phase-h3"
local_temp_dir = "/Users/pete/Desktop/playground/hic2_temp"

resolutions = ["5000", "10000", "25000", "50000", "100000", "250000", "500000", "1000000", "2500000"]
resolution_arg = ",".join(resolutions)
run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

os.makedirs(dropbox_output_dir, exist_ok=True)
os.makedirs(local_temp_dir, exist_ok=True)

master_log_path = os.path.join(dropbox_output_dir, f"phase_h3_matrix_gen_{run_id}.log")
summary_path = os.path.join(dropbox_output_dir, f"phase_h3_matrix_gen_summary_{run_id}.tsv")


def log(message, sample_log_path=None):
    print(message, flush=True)
    with open(master_log_path, "a") as handle:
        handle.write(message + "\n")
    if sample_log_path is not None:
        with open(sample_log_path, "a") as handle:
            handle.write(message + "\n")


def run_cmd(cmd, sample_log_path=None):
    log(f"Running command: {' '.join(cmd)}", sample_log_path)
    start_time = time.time()
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in iter(process.stdout.readline, ""):
        log(f"  [LOG] {line.rstrip()}", sample_log_path)
    process.wait()
    elapsed = time.time() - start_time
    if process.returncode != 0:
        log(f"Error: Command failed with exit code {process.returncode} (took {elapsed:.2f}s)", sample_log_path)
        return False
    log(f"Command succeeded (took {elapsed:.2f}s)", sample_log_path)
    return True


def is_valid_mcool(path, sample_log_path=None):
    if not os.path.exists(path):
        return False
    cmd = ["uv", "run", "--with", "cooler", "cooler", "ls", path]
    process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if process.returncode != 0:
        log(f"Validation failed for {path}: {process.stdout.strip()}", sample_log_path)
        return False
    found = set()
    for line in process.stdout.splitlines():
        marker = "::/resolutions/"
        if marker in line:
            found.add(line.rsplit(marker, 1)[1].strip())
    missing = [res for res in resolutions if res not in found]
    if missing:
        log(f"Validation failed for {path}: missing resolutions {','.join(missing)}", sample_log_path)
        return False
    return True


def copy_and_validate_mcool(src, dest, sample_log_path=None):
    tmp_dest = f"{dest}.tmp_{run_id}"
    if os.path.exists(tmp_dest):
        os.remove(tmp_dest)
    log(f"Copying final .mcool: {src} -> {dest}", sample_log_path)
    start_time = time.time()
    shutil.copy2(src, tmp_dest)
    log(f"Copied to temporary Dropbox file in {time.time() - start_time:.2f}s", sample_log_path)
    if not is_valid_mcool(tmp_dest, sample_log_path):
        log(f"Error: copied .mcool did not pass validation: {tmp_dest}", sample_log_path)
        return False
    os.replace(tmp_dest, dest)
    log(f"Final .mcool validated and moved into place: {dest}", sample_log_path)
    return True


def remove_file(path, sample_log_path=None):
    if os.path.exists(path):
        try:
            os.remove(path)
            log(f"Removed temporary file: {path}", sample_log_path)
        except Exception as exc:
            log(f"Warning: failed to remove temporary file {path}: {exc}", sample_log_path)


def write_summary(rows):
    with open(summary_path, "w") as handle:
        handle.write("sample\tstatus\tinput_hic\toutput_mcool\tnote\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


log("Starting Fast Local-Caching Phase H3: Matrix Generation Pipeline...")
log(f"Run ID: {run_id}")
log(f"Master log: {master_log_path}")
log(f"Summary file: {summary_path}")
log(f"Dropbox output directory: {dropbox_output_dir}")
log(f"Local temporary directory: {local_temp_dir}")

summary_rows = []

for i, sample in enumerate(samples, 1):
    sample_log_path = os.path.join(dropbox_output_dir, f"{sample}.matrix_gen_{run_id}.log")

    log("\n=========================================", sample_log_path)
    log(f"[{i}/{len(samples)}] Processing Sample: {sample}", sample_log_path)
    log("=========================================", sample_log_path)

    db_mcool = os.path.join(dropbox_output_dir, f"{sample}.mcool")
    db_input_hic = dropbox_input_template.format(sample=sample)
    local_hic = os.path.join(local_temp_dir, f"{sample}_inter_30.hic")
    local_temp_cool = os.path.join(local_temp_dir, f"{sample}_temp_5k.cool")
    local_mcool = os.path.join(local_temp_dir, f"{sample}.mcool")

    if is_valid_mcool(db_mcool, sample_log_path):
        log(f"Output already exists and is valid. Skipping: {db_mcool}", sample_log_path)
        summary_rows.append((sample, "skipped_existing_valid", db_input_hic, db_mcool, "valid output already present"))
        write_summary(summary_rows)
        continue

    if os.path.exists(db_mcool):
        invalid_path = f"{db_mcool}.invalid_{run_id}"
        log(f"Existing output is invalid or incomplete. Renaming to {invalid_path}", sample_log_path)
        os.replace(db_mcool, invalid_path)

    if is_valid_mcool(local_mcool, sample_log_path):
        log(f"Found valid local .mcool from an interrupted run: {local_mcool}", sample_log_path)
        if copy_and_validate_mcool(local_mcool, db_mcool, sample_log_path):
            for path in (local_hic, local_temp_cool, local_mcool):
                remove_file(path, sample_log_path)
            summary_rows.append((sample, "completed_from_existing_local_mcool", db_input_hic, db_mcool, "copied valid local output from interrupted run"))
            write_summary(summary_rows)
            continue
        summary_rows.append((sample, "failed_copy_existing_local_mcool", db_input_hic, db_mcool, "valid local output existed but copy/validation failed"))
        write_summary(summary_rows)
        continue

    for path in (local_hic, local_temp_cool, local_mcool):
        remove_file(path, sample_log_path)

    if not os.path.exists(db_input_hic):
        log(f"Error: Dropbox input .hic file not found: {db_input_hic}", sample_log_path)
        summary_rows.append((sample, "missing_input", db_input_hic, db_mcool, "input .hic missing"))
        write_summary(summary_rows)
        continue

    log(f"--- [Step 0] Copying .hic to local SSD for {sample} ---", sample_log_path)
    copy_start = time.time()
    try:
        shutil.copy2(db_input_hic, local_hic)
        log(f"Copied .hic to local SSD in {time.time() - copy_start:.2f}s", sample_log_path)
    except Exception as exc:
        log(f"Error: Failed to copy {db_input_hic} to {local_hic}: {exc}", sample_log_path)
        summary_rows.append((sample, "failed_copy_input", db_input_hic, db_mcool, str(exc)))
        write_summary(summary_rows)
        continue

    log(f"--- [Step 1] Converting .hic to 5kb base .cool for {sample} ---", sample_log_path)
    convert_cmd = [
        "uv", "run", "--with", "hic2cool",
        "hic2cool", "convert",
        local_hic, local_temp_cool,
        "-r", "5000",
    ]
    if not run_cmd(convert_cmd, sample_log_path):
        summary_rows.append((sample, "failed_hic2cool_convert", db_input_hic, db_mcool, "see log"))
        write_summary(summary_rows)
        continue

    log(f"--- [Step 2] Zoomifying and balancing for {sample} ---", sample_log_path)
    zoomify_cmd = [
        "uv", "run", "--with", "cooler",
        "cooler", "zoomify",
        "--balance",
        "-p", "4",
        "--resolutions", resolution_arg,
        "-o", local_mcool,
        local_temp_cool,
    ]
    if not run_cmd(zoomify_cmd, sample_log_path):
        summary_rows.append((sample, "failed_cooler_zoomify", db_input_hic, db_mcool, "see log"))
        write_summary(summary_rows)
        continue

    if not is_valid_mcool(local_mcool, sample_log_path):
        summary_rows.append((sample, "failed_local_mcool_validation", db_input_hic, db_mcool, "local .mcool missing expected resolutions"))
        write_summary(summary_rows)
        continue

    log(f"--- [Step 3] Copying final .mcool to Dropbox for {sample} ---", sample_log_path)
    if not copy_and_validate_mcool(local_mcool, db_mcool, sample_log_path):
        summary_rows.append((sample, "failed_copy_output", db_input_hic, db_mcool, "local .mcool left in temp for resume"))
        write_summary(summary_rows)
        continue

    log(f"--- [Step 4] Cleaning up local temporary files for {sample} ---", sample_log_path)
    for path in (local_hic, local_temp_cool, local_mcool):
        remove_file(path, sample_log_path)

    summary_rows.append((sample, "completed", db_input_hic, db_mcool, "ok"))
    write_summary(summary_rows)

try:
    os.rmdir(local_temp_dir)
    log(f"Removed empty local temporary directory: {local_temp_dir}")
except OSError:
    log(f"Local temporary directory retained because it is not empty: {local_temp_dir}")

completed = sum(1 for row in summary_rows if row[1].startswith("completed") or row[1] == "skipped_existing_valid")
failed = len(summary_rows) - completed
log(f"Phase H3 matrix generation finished. completed_or_skipped={completed}, failed={failed}")
log(f"Summary written: {summary_path}")

if failed:
    log("One or more samples failed. Inspect the master log and per-sample .log files before rerunning.")
else:
    log("All samples completed or were already valid.")
