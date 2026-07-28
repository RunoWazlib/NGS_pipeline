import pipeline.shape_analysis.mutation_profiler

import argparse, time, sys, subprocess
from pathlib import Path

# Start up shape analysis driver
def main():
    print("[*] Starting SHAPE analysis...")
    # Track evaluation time
    start = time.perf_counter()

    # Setup argument parser
    parser = argparse.ArgumentParser(description="SHAPE-analysis of processed sequencing data REQUIRES mpileup analysis of data")
    parser.add_argument("--modified", "-m", type=str, help="Path to the modified sample pipeline output directory")
    parser.add_argument("--unmodified", "-u", type=str, help="Path to the unmodified sample pipeline output directory")
    args = parser.parse_args()

    # If no args, print help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    # Determine valid paths
    # Modified sample
    if not Path(args.modified).is_dir():
        raise NotADirectoryError(f"[!] Modified sample path is not a directory.")
    else:
        modified_dir = Path(args.modified)
    # Unmodified sample
    if not Path(args.unmodified).is_dir():
        raise NotADirectoryError(f"[!] Unmodified sample path is not a directory.")
    else:
        unmodified_dir = Path(args.unmodified)

    # Find mpileup from output
    # Only use the first file in the list since there should only be one of these
    try:
        modified_analysis = list(modified_dir.rglob("Full_analysis.txt"))[0]
        unmodified_analysis = list(unmodified_dir.rglob("Full_analysis.txt"))[0]
    # rglob did not find file
    except IndexError:
        raise IndexError(f"[!] Dir(s) missing mpileup analyses - run 'ngs-pipeline' with config 'do-mpileup-fullanalysis': true")

    # Make a new output directory on modified_dir path for SHAPE analysis results
    shape_output_dir = f"{modified_dir.parent}/SHAPE_analysis"
    subprocess.run(["mkdir", "-p", shape_output_dir])
    # Pass this analysis on to profiler
    shape_data = pipeline.shape_analysis.mutation_profiler.mutation_profiler(modified_analysis, unmodified_analysis, shape_output_dir)

    # End evaluation time
    end = time.perf_counter()
    print(f"[-] SHAPE analysis completed in {end - start:.2f} seconds.")