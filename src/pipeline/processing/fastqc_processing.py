import subprocess, time
from pathlib import Path

def fastq_benchmark_data_paired(reads_path, output_directory):
    try:
        # Evaluate raw read files using fastqc
        command = f"fastqc -o {output_directory}/fastqc_results/ {reads_path['R1']} {reads_path['R2']}"
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
    
        # Unzip data analysis for access later...
        unzip_analysis(f"{output_directory}/fastqc_results/")

    except subprocess.CalledProcessError:
        print(f"[!] FastQC failed! Error: {result.stderr}")

def fastq_benchmark_data_merged(reads_path, output_directory):
    try:
        # Evaluate raw read file using fastqc
        command = f"fastqc -o {output_directory}/fastqc_results/ {reads_path['R1']}"
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)

        # Unzip data analysis for access later...
        unzip_analysis(f"{output_directory}/fastqc_results/")
    
    except subprocess.CalledProcessError:
        print(f"[!] FastQC Failed! Error: {result.stderr}")

def unzip_analysis(directory):
    # Unzip the data analysis yielded from FastQC
    for file in Path(directory).glob("*.zip"):
        command = f"unzip {file}"
        subprocess.run(command, shell=True, check=True)
        
def main(alignment_mode, reads_path, output_directory):
    """Main function to run FastQC program on sequencing read data

    Args:
        alignment_mode (str): The alignment mode being run - included in config.json. This is used to determine what read files to expect (i.e. 2 read files or 1 read file)
        reads_path (str): Where are the read files - included in config.json
        output_directory (_type_): Where output should be written - included in config.json
    """
    # Keep track of evalutation time for benchmarking
    start = time.perf_counter()
    
    # Create output directory for fastqc results if it doesn't exist
    subprocess.run(["mkdir", "-p", f"{output_directory}/fastqc_results"], check=True)

    # Evaluate raw read files using fastqc
    if alignment_mode == "paired-end-mode":
        fastq_benchmark_data_paired(reads_path, output_directory)
    else:
        fastq_benchmark_data_merged(reads_path, output_directory)
    
    # Unzip fastqc analysis to grab in a summary document downstream


    # Print total evaluation time
    end = time.perf_counter()
    print(f"[*] Benchmarking completed in {end - start:.2f} seconds.")