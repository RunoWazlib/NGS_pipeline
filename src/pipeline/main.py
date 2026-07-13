#!/usr/bin/env python3
import pipeline.processing.fastqc_processing, pipeline.processing.q_trimmer, pipeline.alignment.aligner, pipeline.analysis.analysis, pipeline.analysis.association_analysis
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
        print(f"[!] Error: {e}")
        return None

    # Create output directory if it doesn't exist
    subprocess.run(["mkdir", "-p", config["output-directory"]], check=True)

    # Initialize and run scripts based on the configuration
    ### FastQC processing ###
    if config["analysis-parameters"]["do-benchmarks"] == True:
        print("[*] Starting fastqc benchmarks...")
        # Attempt to call on fastqc
        try:
            # Get fastqc version
            command = f"fastqc --version"
            fastqc_ver = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
            with open(f"{config['output-directory']}/dep_versions.log", "a+") as f:
                f.write("-"*10+" FastQC "+"-"*10+"\n"+f"{fastqc_ver.stdout}")
            
            # Pass parameters to the fastq benchmarks script
            pipeline.processing.fastqc_processing.main(config["mode"], config[config["mode"]], config["output-directory"])

        # If initial command failed...
        except subprocess.CalledProcessError:
            print(f"[!] fastqc failed!")
            print(f"{fastqc_ver.stderr}")

    else:
        print("[*] Skipping fastqc benchmarks as per configuration.")

    ### Sequence Pre-processing ###
    if config["analysis-parameters"]["do-processing"] == True:
        print("[*] Starting sequence pre-processing...")
        if config["processing-parameters"]["do-qtrimming"] == True:
            print("[*] Starting quality trimming...")
            try:
                # Pass parameters to the quality trimming script
                pipeline.processing.q_trimmer.main(config["processing-parameters"], config["mode"], config[config["mode"]], config["output-directory"])
            except Exception as e:
                print(f"[!] Quality trimming failed!")
                print(f"{e}")
    else:
        print("[*] Skipping processing as per configuration.")

    ### Bowtie2 Alignment ###
    if config["analysis-parameters"]["do-alignment"] == True:
        print("[*] Starting alignment...")
        try:
            # Print bowtie2 version
            command = f"bowtie2 --version"
            bowtie2_ver = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
            with open(f"{config['output-directory']}/dep_versions.log", "a+") as f:
                    f.write("-"*10+" Bowtie2 "+"-"*10+"\n"+f"{bowtie2_ver.stdout}")

            # Pass alignment parameters to the aligner script
            pipeline.alignment.aligner.main(config["mode"], config["reference-fasta"], config[config["mode"]], config["output-directory"])

        except subprocess.CalledProcessError:
            print(f"[!] bowtie2 failed!")
            print(f"{bowtie2_ver.stderr}")
    else:
        print("[*] Skipping alignment as per configuration.")

    ### Analysis of Aligned Reads ###
    if config["analysis-parameters"]["do-analysis"] == True:
        print("[*] Starting analysis...")
        try:
            # Print samtools version
            command = f"samtools --version"
            samtools_ver = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
            with open(f"{config['output-directory']}/dep_versions.log","a+") as f:
                f.write("-"*10+" samtools "+"-"*10+"\n"+f"{samtools_ver.stdout}")
            
            # Pass analysis parameters to the analysis script
            pipeline.analysis.analysis.main(config["analysis-parameters"], config["reference-fasta"], config["output-directory"])
            pipeline.analysis.association_analysis.main(config["analysis-parameters"], config["reference-fasta"], config["output-directory"])
        except subprocess.CalledProcessError:
            print(f"[!] samtools failed!")
            print(f"{samtools_ver.stderr}")

    else:
        print("[*] Skipping analysis as per configuration.")

    # Print total evaluation time
    end = time.perf_counter()
    print(f"[-] Pipeline completed in {end - start:.2f} seconds.")

# Run the main function when the script is executed
if __name__ == '__main__':
    main()