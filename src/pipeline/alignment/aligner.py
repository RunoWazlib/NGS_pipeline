import subprocess, time

def generate_ref_library(reference_fasta_path, output_directory=None):
    # Generate reference library using bowtie2-build
    try:
        subprocess.run(["mkdir", "-p", f"{output_directory}/ref_lib"], check=True)
        subprocess.run(["bowtie2-build", reference_fasta_path, f"{output_directory}/ref_lib/{reference_fasta_path.split('/')[-1]}"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[!] Error building reference index!")
        print(f"{e}")
        return None
    
    return f"{output_directory}/ref_lib/{reference_fasta_path.split('/')[-1]}"

def align_reads_paired(reference_fasta_path, R1_path, R2_path, output_bam_path):
    # Align reads to the reference library using bowtie2
    try:
        command = f"bowtie2 -x {reference_fasta_path} -1 {R1_path} -2 {R2_path} --local --sensitive-local --ignore-quals --maxins 500 2> {output_bam_path}.log | samtools view -bS - | samtools sort -o {output_bam_path} && samtools index {output_bam_path}"
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print("[!] Error aligning reads")
        print(f"{e}")
        return None

def align_reads_merged(reference_fasta_path, merged_reads_path, output_bam_path):
    # Align merged reads to the reference library using bowtie2
    try:
        command = f"bowtie2 -x {reference_fasta_path} -U {merged_reads_path} --local --sensitive-local --ignore-quals --maxins 500 2> {output_bam_path}.log | samtools view -bS - | samtools sort -o {output_bam_path} && samtools index {output_bam_path}"
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print("[!] Error aligning reads")
        print(f"{e}")
        return None

def main(alignment_mode, reference_fasta_path, reads_path, output_directory):
    # Track the time taken for alignment
    start_time = time.perf_counter()
    # Determine the paired reads path based on the alignment mode
    if alignment_mode == "paired-end-mode":
        R1_path = reads_path["R1"]
        R2_path = reads_path["R2"]
    elif alignment_mode == "merged-mode":
        merged_path = reads_path["R1"]

    # Check if the reference file exists as index for bowtie2
    try:
        # Does the index exist in the output directory?
        with open(f"{output_directory}/ref_lib/{reference_fasta_path.split('/')[-1]}.1.bt2", 'r') as f:
            print(f"[*] Reference file index found in {output_directory}, skipping index building.")
            ref_library = f"{output_directory}/ref_lib/{reference_fasta_path.split('/')[-1]}"

    except FileNotFoundError:
        print(f"[!] Reference file index not found - Building index for '{reference_fasta_path}'...")
        # Build a bowtie2 index for the reference file if it doesn't exist
        ref_library = generate_ref_library(reference_fasta_path, output_directory)
    
    # TODO: move else to main() config validation
    # Align reads to the reference library
    print(f"[*] Aligning reads to reference library...")
    if alignment_mode == "paired-end-mode":
        align_reads_paired(ref_library, R1_path, R2_path, output_directory + "/aligned_reads.bam")
    elif alignment_mode == "merged-mode":
        align_reads_merged(ref_library, merged_path, output_directory + "/aligned_reads.bam")
    else:
        print(f"[!] Invalid alignment mode specified: {alignment_mode}")

    # Track the time taken for alignment
    end_time = time.perf_counter()
    print(f"[*] Alignment completed in {end_time - start_time:.2f} seconds.")