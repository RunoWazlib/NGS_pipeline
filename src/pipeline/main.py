#!/usr/bin/env python3
import pipeline.processing.fastqc_processing, pipeline.processing.q_trimmer, pipeline.alignment.aligner, pipeline.analysis.analysis_controller
from pipeline.config_validator import check_config_options
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

    # Validate configuration options
    try:
        check_config_options(config)
    except ValueError as e:
        print(f"[!] Invalid config! Hint - use 'config-builder' to generate configuration files")
        print(f"Error: {e}")
        return None

    # Create output directory if it doesn't exist
    subprocess.run(["mkdir", "-p", config["output-directory"]], check=True)

    # Initialize and run scripts based on the configuration
    ### FastQC processing ###
    if config["core-parameters"]["do-benchmarks"] == True:
        print("[*] Starting fastqc benchmarks...")
        # Pass parameters to the fastq benchmarks script
        pipeline.processing.fastqc_processing.main(config["mode"], config[config["mode"]], config["output-directory"])
    else:
        print("[*] Skipping fastqc benchmarks as per configuration.")

    ### Sequence Pre-processing ###
    if config["core-parameters"]["do-processing"] == True:
        print("[*] Starting sequence pre-processing...")
        if config["processing-parameters"]["do-qtrimming"] == True:
            print("[*] Starting quality trimming...")
            # Pass parameters to the quality trimming script
            pipeline.processing.q_trimmer.main(config["processing-parameters"], config["mode"], config[config["mode"]], config["output-directory"])
        # TODO - other pre-processing steps should be implemented here
    else:
        print("[*] Skipping processing as per configuration.")

    ### Bowtie2 Alignment ###
    if config["core-parameters"]["do-alignment"] == True:
        print("[*] Starting alignment...")
        # Pass alignment parameters to the aligner script
        pipeline.alignment.aligner.main(config["mode"], config["reference-fasta"], config[config["mode"]], config["output-directory"])
    else:
        print("[*] Skipping alignment as per configuration.")

    ### Analysis of Aligned Reads ###
    if config["core-parameters"]["do-analysis"] == True:
        print("[*] Starting analysis...")
        # Pass analysis parameters to the analysis script
        pipeline.analysis.analysis_controller.main(config["analysis-parameters"], config["reference-fasta"], config["output-directory"])
    else:
        print("[*] Skipping analysis as per configuration.")

    # Print total evaluation time
    end = time.perf_counter()
    print(f"[-] Pipeline completed in {end - start:.2f} seconds.")

# Run the main function when the script is executed
if __name__ == '__main__':
    main()