#!/usr/bin/env python3
import src.lib.fastq_benchmarks, src.lib.aligner, src.lib.analysis, src.lib.association_analysis
import json, argparse, time, subprocess

def load_config(config_file_path):
    with open(config_file_path, 'r') as f:
        config = json.load(f)
    return config

def main():
    print("[*] Starting next-gen sequencing analysis pipeline...")
    # Keep track of evalutation time for benchmarking
    start = time.perf_counter()

    # Set up argument parser
    parser = argparse.ArgumentParser(description='Next-gen sequencing analysis driver.')
    parser.add_argument('--config', type=str, help='Path to the configuration file.')
    args = parser.parse_args()

    # Load configuration
    config = load_config(args.config)

    # Create output directory if it doesn't exist
    subprocess.run(["mkdir", "-p", config["output-directory"]], check=True)

    # Initialize and run scripts based on the configuration
    if config["analysis-parameters"]["do-benchmarks"]:
        print("[*] Starting fastqc benchmarks...")
        # Print fastqc version
        command = f"fastqc --version"
        subprocess.run(command, shell=True, check=True)

        # Pass benchmarking parameters to the fastq_benchmarks script
        src.lib.fastq_benchmarks.main(config["mode"], config[config["mode"]], config["output-directory"])
    else:
        print("[*] Skipping fastqc benchmarks as per configuration.")

    if config["analysis-parameters"]["do-alignment"]:
        print("[*] Starting alignment...")
        # Print bowtie2 version
        command = f"bowtie2 --version"
        subprocess.run(command, shell=True, check=True)
        # Pass alignment parameters to the aligner script
        src.lib.aligner.main(config["mode"], config["reference-fasta"], config[config["mode"]], config["output-directory"])
    else:
        print("[*] Skipping alignment as per configuration.")

    if config["analysis-parameters"]["do-analysis"]:
        print("[*] Starting analysis...")
        src.lib.analysis.main(config["analysis-parameters"], config["reference-fasta"], config["output-directory"])
        src.lib.association_analysis.main(config["analysis-parameters"], config["reference-fasta"], config["output-directory"])
    else:
        print("[*] Skipping analysis as per configuration.")

    # Print total evaluation time
    end = time.perf_counter()
    print(f"[-] Pipeline completed in {end - start:.2f} seconds.")

# Run the main function when the script is executed
if __name__ == '__main__':
    main()