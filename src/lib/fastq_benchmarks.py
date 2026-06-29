import subprocess, time

def fastq_benchmark_data_paired(reads_path, output_directory):
    # Evaluate raw read files using fastqc
    command = f"fastqc -o {output_directory}/fastqc_results/ {reads_path['R1']} {reads_path['R2']}"
    subprocess.run(command, shell=True, check=True)

def fastq_benchmark_data_merged(reads_path, output_directory):
    # Evaluate raw read file using fastqc
    command = f"fastqc -o {output_directory}/fastqc_results/ {reads_path['R1']}"
    subprocess.run(command, shell=True, check=True)

def main(alignment_mode, reads_path, output_directory):
    # Keep track of evalutation time for benchmarking
    start = time.perf_counter()
    
    # Create output directory for fastqc results if it doesn't exist
    subprocess.run(["mkdir", "-p", f"{output_directory}/fastqc_results"], check=True)

    # Evaluate raw read files using fastqc
    if alignment_mode == "paired-end-mode":
        fastq_benchmark_data_paired(reads_path, output_directory)
    else:
        fastq_benchmark_data_merged(reads_path, output_directory)
    
    # Print total evaluation time
    end = time.perf_counter()
    print(f"[*] Benchmarking completed in {end - start:.2f} seconds.")