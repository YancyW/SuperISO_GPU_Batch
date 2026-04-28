#!/usr/bin/env python3

import argparse
import concurrent.futures
import csv
import os
import re
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Check np_data scan rows with SuperIso by generating an SLHA/LHA file via 2HDMC "
            "and then running slha_chi2.x or slha.x."
        )
    )
    parser.add_argument("input_file", help="Path to scan_after_*.dat file")
    parser.add_argument(
        "--rows",
        default="",
        help=(
            "Rows to test: '1', '1,5,9', or '1-10'. Line numbers are 1-based and exclude the header. "
            "When omitted, non-table mode tests row 1 and append-table mode processes the full file."
        ),
    )
    parser.add_argument(
        "--mode",
        choices=["chi2", "obs"],
        default="chi2",
        help="Use slha_chi2.x for a global chi2 check or slha.x for observable output.",
    )
    parser.add_argument(
        "--superiso-root",
        default=str(Path(__file__).resolve().parents[1]),
        help="SuperIso root directory. Defaults to the parent of np_data/.",
    )
    parser.add_argument(
        "--twohdmc-bin",
        default="/home/yancy/Software/Bin/HiggsRelated/2HDMC-1.8.0/install/bin",
        help="Real 2HDMC bin directory containing CalcPhys and CalcPhys.bin.",
    )
    parser.add_argument(
        "--twohdmc-batch-bin",
        default="",
        help=(
            "Optional path to a CalcPhysBatch executable. When available, append-table mode uses it "
            "to generate LHA files for many rows in one 2HDMC process."
        ),
    )
    parser.add_argument(
        "--dataset-type-to-superiso",
        default="auto",
        help=(
            "Map the 'type' column before passing it to 2HDMC/SuperIso. "
            "Use 'auto' to keep the dataset value, or explicit mappings like '4:3,3:3'."
        ),
    )
    parser.add_argument(
        "--append-table",
        action="store_true",
        help=(
            "Write a new table with SuperIso results appended after the existing columns. "
            "In this mode, rows not selected by --rows are also processed unless --rows is explicitly set."
        ),
    )
    parser.add_argument(
        "--output-file",
        default="",
        help="Output path for the appended table. Defaults to <input>_superiso.dat when --append-table is used.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel workers to use in append-table mode. Default: 1.",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=50,
        help="Write a partial appended table every N processed rows in append-table mode. Default: 50.",
    )
    return parser.parse_args()


def parse_row_selection(spec):
    selected = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" in chunk:
            start_str, end_str = chunk.split("-", 1)
            start = int(start_str)
            end = int(end_str)
            if end < start:
                raise ValueError(f"Invalid row range: {chunk}")
            selected.extend(range(start, end + 1))
        else:
            selected.append(int(chunk))
    if not selected:
        raise ValueError("No rows selected")
    return sorted(set(selected))


def parse_type_map(spec):
    if spec == "auto":
        return None
    mapping = {}
    for part in spec.split(","):
        left, right = part.split(":", 1)
        mapping[int(left.strip())] = int(right.strip())
    return mapping


def load_rows(input_path, wanted_rows):
    rows = []
    with open(input_path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        for line_number, row in enumerate(reader, start=1):
            if line_number in wanted_rows:
                cleaned_row = {key.lstrip("\ufeff"): value for key, value in row.items()}
                rows.append((line_number, cleaned_row))
    return rows


def load_all_rows(input_path):
    rows = []
    with open(input_path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        for line_number, row in enumerate(reader, start=1):
            cleaned_row = {key.lstrip("\ufeff"): value for key, value in row.items()}
            rows.append((line_number, cleaned_row))
    return rows


def run_command(command, cwd=None, env=None):
    return subprocess.run(command, cwd=cwd, env=env, text=True, capture_output=True)


def log_message(message=""):
    if message:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] {message}", flush=True)
    else:
        print("", flush=True)


def detect_runtime_lib():
    candidates = []
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.append(Path(conda_prefix) / "lib")

    candidates.extend(
        [
            Path("/home/yancy/Package_Management/anaconda3/envs/root63604py311/lib"),
            Path("/usr/lib64"),
            Path("/usr/lib/x86_64-linux-gnu"),
        ]
    )

    for candidate in candidates:
        if (candidate / "libgsl.so").exists() or (candidate / "libgsl.so.25").exists():
            return str(candidate)
    return ""


def extract_chi2(text):
    match = re.search(r"chi2\s*=\s*([0-9.+\-eE]+)", text)
    return float(match.group(1)) if match else None


def extract_nobs(text):
    match = re.search(r"n_obs\s*=\s*(\d+)", text)
    return int(match.group(1)) if match else None


def load_myobs(superiso_root):
    obs_path = Path(superiso_root) / "chi2_input" / "myobs.in"
    return [line.strip() for line in obs_path.read_text().splitlines() if line.strip()]


def sanitize_observable_name(name):
    return (
        name.replace("/", "_over_")
        .replace("-", "minus")
        .replace("'", "prime")
        .replace("(", "_")
        .replace(")", "_")
        .replace(",", "_")
        .replace(" ", "_")
    )


def parse_slha_chi2_output(text, observable_names):
    result = {name: "" for name in observable_names}
    obs_lines = [
        line for line in text.splitlines() if line.startswith("\t") and "=" in line and "(exp:" in line
    ]

    for idx, line in enumerate(obs_lines):
        if idx >= len(observable_names):
            break
        match = re.match(r"^\t(.+?)\s*=\s*\t?([^\s]+)", line)
        if match:
            result[observable_names[idx]] = match.group(2)
    return result


def determine_status(text, lha_exists, calcphys_returncode):
    if not lha_exists:
        return "2hdmc_failed" if calcphys_returncode != 0 else "lha_missing"
    if "Invalid point" in text:
        return "invalid_point"
    if "Invalid SLHA file" in text:
        return "invalid_slha"
    if "Model not yet implemented" in text:
        return "model_not_implemented"
    return "ok"


def build_lha(row, line_number, temp_dir, calcphys_bin, runtime_lib):
    lha_path = Path(temp_dir) / f"row_{line_number}.lha"
    mh = float(row["mh"])
    mH = float(row["mH"])
    mA = float(row["mA"])
    mC = float(row["mC"])
    sba = float(row["sba"])
    tb = float(row["tb"])
    l6 = float(row["l6"])
    l7 = float(row["l7"])
    m12_2 = float(row["m12_2"])
    model_type = int(float(row["type"]))

    env = os.environ.copy()
    if runtime_lib:
        env["LD_LIBRARY_PATH"] = f"{runtime_lib}:{env.get('LD_LIBRARY_PATH', '')}".rstrip(":")

    result = run_command(
        [
            str(calcphys_bin),
            str(mh),
            str(mH),
            str(mA),
            str(mC),
            str(sba),
            str(l6),
            str(l7),
            str(m12_2),
            str(tb),
            str(model_type),
            str(lha_path),
        ],
        env=env,
    )
    return lha_path, result


def build_lhas_batch(prepared_rows, temp_dir, calcphys_batch_bin, runtime_lib):
    batch_input_path = Path(temp_dir) / "calcphys_batch_input.tsv"
    lha_output_dir = Path(temp_dir) / "lha_batch"
    lha_output_dir.mkdir(parents=True, exist_ok=True)

    with open(batch_input_path, "w", encoding="utf-8") as handle:
        for line_number, _original_row, working_row, _dataset_type, mapped_type in prepared_rows:
            handle.write(
                f"{line_number} {working_row['mh']} {working_row['mH']} {working_row['mA']} {working_row['mC']} "
                f"{working_row['sba']} {working_row['l6']} {working_row['l7']} {working_row['m12_2']} {working_row['tb']} {mapped_type}\n"
            )

    env = os.environ.copy()
    if runtime_lib:
        env["LD_LIBRARY_PATH"] = f"{runtime_lib}:{env.get('LD_LIBRARY_PATH', '')}".rstrip(":")

    batch_run = run_command([str(calcphys_batch_bin), str(batch_input_path), str(lha_output_dir)], env=env)
    results = {}

    for line in (batch_run.stdout or "").splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        row_id, status, output_file, message = parts[0], parts[1], parts[2], parts[3]
        try:
            line_number = int(row_id)
        except ValueError:
            continue

        returncode = 0 if status == "ok" else 1
        calcphys = subprocess.CompletedProcess(
            args=[str(calcphys_batch_bin), row_id],
            returncode=returncode,
            stdout=line,
            stderr="",
        )
        results[line_number] = {
            "line_number": line_number,
            "lha_path": Path(output_file),
            "calcphys": calcphys,
            "message": message,
            "batch_returncode": batch_run.returncode,
        }

    for line_number, *_rest in prepared_rows:
        if line_number not in results:
            missing_lha = lha_output_dir / f"{line_number}.lha"
            results[line_number] = {
                "line_number": line_number,
                "lha_path": missing_lha,
                "calcphys": subprocess.CompletedProcess(
                    args=[str(calcphys_batch_bin), str(line_number)],
                    returncode=1,
                    stdout="",
                    stderr=batch_run.stderr or "",
                ),
                "message": "missing_batch_result",
                "batch_returncode": batch_run.returncode,
            }

    return results, batch_run


def evaluate_existing_lha(batch_result, checker_bin, observable_names, cwd):
    lha_path = batch_result["lha_path"]
    combined = ""
    check = None

    if lha_path.exists():
        check = run_command([str(checker_bin), str(lha_path)], cwd=str(cwd))
        combined = (check.stdout or "") + (check.stderr or "")

    status = determine_status(combined, lha_path.exists(), batch_result["calcphys"].returncode)
    parsed_obs = parse_slha_chi2_output(combined, observable_names) if checker_bin.name == "slha_chi2.x" else {}

    return {
        "line_number": batch_result["line_number"],
        "lha_path": lha_path,
        "calcphys": batch_result["calcphys"],
        "check": check,
        "combined": combined,
        "status": status,
        "chi2": extract_chi2(combined),
        "n_obs": extract_nobs(combined),
        "observables": parsed_obs,
    }


def evaluate_row(row, line_number, temp_dir, calcphys_bin, runtime_lib, checker_bin, observable_names, cwd):
    lha_path, calcphys = build_lha(row, line_number, temp_dir, calcphys_bin, runtime_lib)
    combined = ""
    check = None
    if lha_path.exists():
        check = run_command([str(checker_bin), str(lha_path)], cwd=str(cwd))
        combined = (check.stdout or "") + (check.stderr or "")

    status = determine_status(combined, lha_path.exists(), calcphys.returncode)
    parsed_obs = parse_slha_chi2_output(combined, observable_names) if checker_bin.name == "slha_chi2.x" else {}

    return {
        "line_number": line_number,
        "lha_path": lha_path,
        "calcphys": calcphys,
        "check": check,
        "combined": combined,
        "status": status,
        "chi2": extract_chi2(combined),
        "n_obs": extract_nobs(combined),
        "observables": parsed_obs,
    }


def make_appended_row(original_row, dataset_type, mapped_type, evaluation, observable_names):
    out_row = dict(original_row)
    out_row["superiso_dataset_type"] = str(dataset_type)
    out_row["superiso_type_used"] = str(mapped_type)
    out_row["superiso_status"] = evaluation["status"]
    out_row["superiso_calcphys_rc"] = str(evaluation["calcphys"].returncode)
    out_row["superiso_rc"] = str(evaluation["check"].returncode if evaluation["check"] else "")
    out_row["superiso_chi2"] = "" if evaluation["chi2"] is None else f"{evaluation['chi2']:.8g}"
    out_row["superiso_n_obs"] = "" if evaluation["n_obs"] is None else str(evaluation["n_obs"])
    for observable_name in observable_names:
        out_row[f"superiso_{sanitize_observable_name(observable_name)}"] = evaluation["observables"].get(observable_name, "")
    return out_row


def default_output_path(input_path):
    if input_path.suffix:
        return input_path.with_name(f"{input_path.stem}_superiso{input_path.suffix}")
    return input_path.with_name(f"{input_path.name}_superiso.dat")


def write_appended_table(output_path, original_fieldnames, rows_out):
    extra_fieldnames = []
    if rows_out:
        for key in rows_out[0].keys():
            if key not in original_fieldnames:
                extra_fieldnames.append(key)
    fieldnames = original_fieldnames + extra_fieldnames

    temp_output_path = output_path.with_name(f".{output_path.name}.tmp")

    with open(temp_output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=" ", extrasaction="ignore")
        writer.writeheader()
        for row in rows_out:
            writer.writerow(row)

    temp_output_path.replace(output_path)


def build_partial_rows(prepared_rows, result_map):
    partial_rows = []
    for line_number, *_rest in prepared_rows:
        if line_number in result_map:
            partial_rows.append(result_map[line_number])
    return partial_rows


def main():
    args = parse_args()
    input_path = Path(args.input_file).resolve()
    superiso_root = Path(args.superiso_root).resolve()
    twohdmc_bin = Path(args.twohdmc_bin).resolve()
    calcphys_batch_bin = (
        Path(args.twohdmc_batch_bin).resolve() if args.twohdmc_batch_bin else (twohdmc_bin / "CalcPhysBatch")
    )

    calcphys_bin = twohdmc_bin / "CalcPhys.bin"
    slha_chi2_bin = superiso_root / "slha_chi2.x"
    slha_obs_bin = superiso_root / "slha.x"

    if not input_path.exists():
        print(f"Input file not found: {input_path}", file=sys.stderr, flush=True)
        return 1
    if not calcphys_bin.exists():
        print(f"CalcPhys.bin not found: {calcphys_bin}", file=sys.stderr, flush=True)
        return 1

    checker_bin = slha_chi2_bin if args.mode == "chi2" else slha_obs_bin
    if not checker_bin.exists():
        print(f"SuperIso binary not found: {checker_bin}", file=sys.stderr, flush=True)
        return 1

    type_map = parse_type_map(args.dataset_type_to_superiso)
    runtime_lib = detect_runtime_lib()
    observable_names = load_myobs(superiso_root) if args.mode == "chi2" else []

    if args.append_table:
        rows = load_all_rows(input_path)
        if args.rows:
            wanted_rows = set(parse_row_selection(args.rows))
            rows = [item for item in rows if item[0] in wanted_rows]
    else:
        wanted_rows = parse_row_selection(args.rows) if args.rows else [1]
        rows = load_rows(input_path, wanted_rows)

    if not rows:
        print("No selected rows were found in the input file.", file=sys.stderr, flush=True)
        return 1

    log_message(f"Input file: {input_path}")
    log_message(f"SuperIso root: {superiso_root}")
    log_message(f"2HDMC bin: {twohdmc_bin}")
    if args.append_table and calcphys_batch_bin.exists():
        log_message(f"2HDMC batch bin: {calcphys_batch_bin}")
    log_message(f"Check mode: {args.mode}")
    if args.append_table:
        log_message("Append-table mode: enabled")
        log_message(f"Parallel jobs: {max(1, args.jobs)}")
        log_message(f"Checkpoint every: {max(1, args.checkpoint_every)}")
    log_message()

    with tempfile.TemporaryDirectory(prefix="superiso_npdata_") as temp_dir:
        appended_rows = []
        output_path = None
        original_fieldnames = None
        prepared_rows = []
        for line_number, row in rows:
            original_row = dict(row)
            dataset_type = int(float(original_row["type"]))
            mapped_type = type_map.get(dataset_type, dataset_type) if type_map else dataset_type
            working_row = dict(original_row)
            working_row["type"] = str(mapped_type)
            prepared_rows.append((line_number, original_row, working_row, dataset_type, mapped_type))

        if args.append_table:
            all_rows = load_all_rows(input_path)
            original_fieldnames = list(all_rows[0][1].keys())
            output_path = Path(args.output_file).resolve() if args.output_file else default_output_path(input_path)

        batch_results = None
        if args.append_table and calcphys_batch_bin.exists():
            batch_results, batch_run = build_lhas_batch(prepared_rows, temp_dir, calcphys_batch_bin, runtime_lib)
            log_message(
                f"2HDMC batch generation finished with rc={batch_run.returncode} for {len(prepared_rows)} requested rows"
            )

        if args.append_table and max(1, args.jobs) > 1:
            def task(payload):
                line_number, original_row, working_row, dataset_type, mapped_type = payload
                if batch_results is not None:
                    evaluation = evaluate_existing_lha(
                        batch_results[line_number],
                        checker_bin,
                        observable_names,
                        superiso_root,
                    )
                else:
                    evaluation = evaluate_row(
                        working_row,
                        line_number,
                        temp_dir,
                        calcphys_bin,
                        runtime_lib,
                        checker_bin,
                        observable_names,
                        superiso_root,
                    )
                appended_row = make_appended_row(original_row, dataset_type, mapped_type, evaluation, observable_names)
                return line_number, appended_row

            result_map = {}
            with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as executor:
                future_map = {executor.submit(task, payload): payload[0] for payload in prepared_rows}
                completed = 0
                for future in concurrent.futures.as_completed(future_map):
                    line_number, appended_row = future.result()
                    result_map[line_number] = appended_row
                    completed += 1
                    if completed % 10 == 0 or completed == len(prepared_rows):
                        log_message(f"Processed {completed}/{len(prepared_rows)} rows...")
                    if completed % max(1, args.checkpoint_every) == 0 or completed == len(prepared_rows):
                        partial_rows = build_partial_rows(prepared_rows, result_map)
                        write_appended_table(output_path, original_fieldnames, partial_rows)
                        log_message(f"Checkpoint written: {len(partial_rows)}/{len(prepared_rows)} rows -> {output_path}")

            appended_rows = [result_map[line_number] for line_number, *_ in prepared_rows]
        else:
            for line_number, original_row, working_row, dataset_type, mapped_type in prepared_rows:
                if batch_results is not None:
                    evaluation = evaluate_existing_lha(
                        batch_results[line_number],
                        checker_bin,
                        observable_names,
                        superiso_root,
                    )
                else:
                    evaluation = evaluate_row(
                        working_row,
                        line_number,
                        temp_dir,
                        calcphys_bin,
                        runtime_lib,
                        checker_bin,
                        observable_names,
                        superiso_root,
                    )

                if args.append_table:
                    appended_rows.append(
                        make_appended_row(original_row, dataset_type, mapped_type, evaluation, observable_names)
                    )
                    if len(appended_rows) % 10 == 0 or len(appended_rows) == len(prepared_rows):
                        log_message(f"Processed {len(appended_rows)}/{len(prepared_rows)} rows...")
                    if len(appended_rows) % max(1, args.checkpoint_every) == 0 or len(appended_rows) == len(prepared_rows):
                        write_appended_table(output_path, original_fieldnames, appended_rows)
                        log_message(f"Checkpoint written: {len(appended_rows)}/{len(prepared_rows)} rows -> {output_path}")
                    continue

                print(f"=== Row {line_number} ===")
                print(
                    "params: "
                    f"mh={working_row['mh']} mH={working_row['mH']} mA={working_row['mA']} mC={working_row['mC']} "
                    f"sba={working_row['sba']} tb={working_row['tb']} m12_2={working_row['m12_2']} "
                    f"type(dataset={dataset_type} -> superiso={mapped_type})"
                )
                print(f"2HDMC return code: {evaluation['calcphys'].returncode}")
                stdout = evaluation["calcphys"].stdout.strip()
                stderr = evaluation["calcphys"].stderr.strip()
                if stdout:
                    print("2HDMC stdout:")
                    print(stdout)
                if stderr:
                    print("2HDMC stderr:")
                    print(stderr)

                if evaluation["status"] == "2hdmc_failed" or evaluation["status"] == "lha_missing":
                    print("SuperIso status: skipped, because 2HDMC did not produce an LHA file")
                    print()
                    continue

                print(f"SuperIso return code: {evaluation['check'].returncode}")
                print(f"SuperIso status: {evaluation['status']}")

                if evaluation["chi2"] is not None:
                    print(f"SuperIso chi2: {evaluation['chi2']:.4f}")

                key_lines = []
                for line in evaluation["combined"].splitlines():
                    if (
                        "chi2" in line
                        or "Invalid point" in line
                        or "BR(b->s gamma)" in line
                        or "BR(Bs->mu mu)" in line
                        or "a_muon" in line
                    ):
                        key_lines.append(line)

                if key_lines:
                    print("SuperIso excerpt:")
                    for line in key_lines[:12]:
                        print(line)
                print()

        if args.append_table:
            write_appended_table(output_path, original_fieldnames, appended_rows)
            log_message(f"Appended table written to: {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())