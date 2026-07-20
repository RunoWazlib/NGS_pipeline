import subprocess
import matplotlib.pyplot as plt

def generate_alignment_stats(aligned_bam_path, output_directory):
    # Generate alignment statistics using samtools flagstat and samtools stat
    command = f"samtools flagstat {aligned_bam_path} > {output_directory}/alignment_stats_summary.txt"
    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)

    # Generate in-depth alignment statistics using samtools stat
    stats_file = f"{output_directory}/full_alignment_stats.txt"
    command = f"samtools stat {aligned_bam_path} > {stats_file}"
    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)

    # Define a mapping of prefixes to output filenames
    prefix_mapping = {
        "SN": "full_alignment_stats_summary.txt",
        "FFQ": "full_alignment_stats_R1_qualities.tsv",
        "LFQ": "full_alignment_stats_R2_qualities.tsv",
        "IS": "full_alignment_stats_insert_sizes.tsv",
        "RL": "full_alignment_stats_all_read_lengths.tsv",
        "FRL": "full_alignment_stats_R1_read_lengths.tsv",
        "LRL": "full_alignment_stats_R2_read_lengths.tsv",
        "MAPQ": "full_alignment_stats_mapping_quality.tsv",
        "ID": "full_alignment_stats_indel_dist.tsv",
        "IC": "full_alignment_stats_indels_per_cycle.tsv",
    }

    # Open all target files at once
    out_files = {}
    for prefix, filename in prefix_mapping.items():
        out_files[prefix] = open(f"{output_directory}/{filename}", "w")

    # Iter over the stats file and write to the appropriate output files based on the prefix
    with open(stats_file, "r") as f:
        for line in f:
            for prefix, out_f in out_files.items():
                if line.startswith(f"{prefix}\t"):
                    # Split off the prefix
                    out_f.write(line.split("\t", 1)[1])
    
    # Remove the full alignment stats file after processing
    subprocess.run(["rm", stats_file], check=True)

    # Close all files cleanly
    for out_f in out_files.values():
        out_f.close()

def generate_alignment_score_plot(aligned_bam_path, output_directory):
    # Generate alignment score plot using samtools view to extract alignment scores and matplotlib to create a histogram
    command = f"samtools view {aligned_bam_path} | awk '{{print $5}}' > {output_directory}/alignment_scores.txt"
    subprocess.run(command, shell=True, check=True)
    
    # Read alignment scores and generate a histogram
    with open(f"{output_directory}/alignment_scores.txt", 'r') as f:
        scores = [int(line.strip()) for line in f]
        
        # Generate and save the alignment score histogram
        plt.hist(scores, bins=50, color='cyan', edgecolor='black')
        plt.title('Alignment Score Distribution')
        plt.xlabel('Alignment Score')
        plt.ylabel('Frequency')
        plt.xlim(0, round(max(scores)))
        plt.tight_layout()
        plt.savefig(f"{output_directory}/alignment_score_plot.png")
        plt.clf()