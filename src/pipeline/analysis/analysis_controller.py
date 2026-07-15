from . import alignment_stats, alignment_analysis, mpileup_analysis, association_analysis
import time, subprocess

def main(analysis_params, reference_fasta_path, output_directory):
    # Track the time taken for analysis
    start_time = time.perf_counter()

    # Perform alignment statistics analysis if enabled
    # Alignment statistics analysis includes generating summary statistics and detailed statistics for the aligned reads
    if analysis_params.get("do-alignment-stats", True):
        print("[*] Performing alignment statistics analysis...")
        # Make sub-output directory - no error/overwrite if it already exists
        subprocess.run(["mkdir", "-p", f"{output_directory}/alignment_stats"], check=True)
        alignment_stats.generate_alignment_stats(f"{output_directory}/aligned_reads.bam", f"{output_directory}/alignment_stats")
    
    # Alignment score plot is a histogram of the alignment scores for the aligned reads
    if analysis_params.get("do-alignment-score-plot", True):
        print("[*] Generating alignment score plot...")
        # Make sub-output directory - no error/overwrite if it already exists
        subprocess.run(["mkdir", "-p", f"{output_directory}/alignment_stats"], check=True)
        alignment_stats.generate_alignment_score_plot(f"{output_directory}/aligned_reads.bam", f"{output_directory}/alignment_stats")

    # Alignment visualization is a global alignment-esque visual of all sequences
    if analysis_params.get("do-alignment-visualization", True):
        print("[*] Generating alignment visualization...")
        # Make sub-output directory - no error/overwrite if it already exists
        subprocess.run(["mkdir","-p",f"{output_directory}/alignment"])
        alignment_analysis.generate_alignment_visualization(f"{output_directory}/aligned_reads.bam", reference_fasta_path, f"{output_directory}/alignment")
        print("[*] Alignment visualization generated.")

    # mpileup analysis counts the base calls in each read and generates a table
    if analysis_params.get("do-mpileup", True):
        print("[*] Processing mpileup...")
        # Make sub-output directory - no error/overwrite if it already exists
        subprocess.run(["mkdir","-p",f"{output_directory}/alignment"])
        mpileup_file = mpileup_analysis.generate_mpileup(f"{output_directory}/aligned_reads.bam", reference_fasta_path, f"{output_directory}/alignment")
        
        # Further analysis of mpileup output
        if analysis_params.get("do-mpileup-fullanalysis", True):
            mpileup_analysis.generate_mpileup_full_analysis(mpileup_file, f"{output_directory}/alignment/")
        if analysis_params.get("do-mpileup-simpleanalysis", True):
            mpileup_analysis.generate_mpileup_simple_analysis(mpileup_file, f"{output_directory}/alignment/")
        if analysis_params.get("do-mpileup-visualization", True):
            mpileup_analysis.generate_mpileup_visualization(mpileup_file, f"{output_directory}/alignment")
    
    # Association analysis compares all mutation rates for full read sequences for correlation - indicates tertiary interactions
    if analysis_params.get("do-association-analysis", True):
        print("[*] Generating association analysis...")
        # Make sub-output directory - no error/overwrite if it already exists
        subprocess.run(["mkdir","-p",f"{output_directory}/associations"])
        association_analysis.main(reference_fasta_path, f"{output_directory}/aligned_reads.bam", f"{output_directory}/associations")

    # Track the time taken for analysis
    end_time = time.perf_counter()
    print(f"[*] Analysis completed in {end_time - start_time:.2f} seconds.")
