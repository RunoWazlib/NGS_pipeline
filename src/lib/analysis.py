import subprocess, time, re
import matplotlib.pyplot as plt

def generate_alignment_stats(aligned_bam_path, output_directory):
    # Make sub-output directory
    subprocess.run(["mkdir", "-p", f"{output_directory}/alignment_stats"], check=True)
    
    # Generate alignment statistics using samtools flagstat
    command = f"samtools flagstat {aligned_bam_path} > {output_directory}/alignment_stats/alignment_stats.txt"
    subprocess.run(command, shell=True, check=True)
    
    # Generate in-depth alignment statistics using samtools stat
    command = f"samtools stat {aligned_bam_path} > {output_directory}/alignment_stats/full_alignment_stats.txt"
    subprocess.run(command, shell=True, check=True)
    # Parse out each component of the larger alignment stats file
    command = f"grep ^SN {output_directory}/alignment_stats/full_alignment_stats.txt | cut -f 2- > {output_directory}/alignment_stats/full_alignment_stats_summary.txt"
    subprocess.run(command,shell=True, check=True)
    command = f"grep ^FFQ {output_directory}/alignment_stats/full_alignment_stats.txt | cut -f 2- > {output_directory}/alignment_stats/full_alignment_stats_R1_qualities.tsv"
    subprocess.run(command,shell=True, check=True)
    command = f"grep ^LFQ {output_directory}/alignment_stats/full_alignment_stats.txt | cut -f 2- > {output_directory}/alignment_stats/full_alignment_stats_R2_qualities.tsv"
    subprocess.run(command,shell=True, check=True)
    command = f"grep ^IS {output_directory}/alignment_stats/full_alignment_stats.txt | cut -f 2- > {output_directory}/alignment_stats/full_alignment_stats_insert_sizes.tsv"
    subprocess.run(command,shell=True, check=True)
    command = f"grep ^RL {output_directory}/alignment_stats/full_alignment_stats.txt | cut -f 2- > {output_directory}/alignment_stats/full_alignment_stats_all_read_lengths.tsv"
    subprocess.run(command,shell=True, check=True)
    command = f"grep ^FRL {output_directory}/alignment_stats/full_alignment_stats.txt | cut -f 2- > {output_directory}/alignment_stats/full_alignment_stats_R1_read_lengths.tsv"
    subprocess.run(command,shell=True, check=True)
    command = f"grep ^LRL {output_directory}/alignment_stats/full_alignment_stats.txt | cut -f 2- > {output_directory}/alignment_stats/full_alignment_stats_R2_read_lengths.tsv"
    subprocess.run(command,shell=True, check=True)
    command = f"grep ^MAPQ {output_directory}/alignment_stats/full_alignment_stats.txt | cut -f 2- > {output_directory}/alignment_stats/full_alignment_stats_mapping_quality.tsv"
    subprocess.run(command,shell=True, check=True)
    command = f"grep ^ID {output_directory}/alignment_stats/full_alignment_stats.txt | cut -f 2- > {output_directory}/alignment_stats/full_alignment_stats_indel_dist.tsv"
    subprocess.run(command,shell=True, check=True)
    command = f"grep ^IC {output_directory}/alignment_stats/full_alignment_stats.txt | cut -f 2- > {output_directory}/alignment_stats/full_alignment_stats_indel_dist.tsv"
    subprocess.run(command,shell=True, check=True)

def generate_alignment_visualization(aligned_bam_path, reference_file, output_directory):
    # # Generate alignment visualization using samtools tview
    # command = f"samtools tview -d T -w 250 {aligned_bam_path} > {output_directory}/alignment_visualization.txt"
    # subprocess.run(command, shell=True, check=True)

    # Get python-parsable file of all sequences
    command = f"samtools view {aligned_bam_path} > {output_directory}/all_reads.txt"
    subprocess.run(command, shell=True, check=True)

    with open(f"{output_directory}/global_alignment_matrix.txt", "w") as output_file:
        total_reads = 0
        aligned_reads = 0
        data_out = {}
        with open(f"{output_directory}/all_reads.txt","r") as source_file:
            for line in source_file:
                # Count total number of reads (including unaligned)
                total_reads += 1
                # Get entry data
                fields = line.strip().split("\t")
                read_name = fields[0]
                start_ref = int(fields[3])
                map_qual = int(fields[4])
                cigar_str = fields[5]
                read_seq = fields[9]
                # If read is unaligned, move on
                if cigar_str == "*":
                    data_out[total_reads] = [read_name, start_ref, cigar_str, read_seq, map_qual]
                else:
                    aligned_reads += 1
                    data_out[total_reads] = [read_name, start_ref, cigar_str, read_seq, map_qual]
            
        # Output file header
        output_file.write(f"# Number of input reads: {total_reads}\n")
        output_file.write(f"# Number of aligned/mapped reads: {aligned_reads} - {aligned_reads*100/total_reads:.2f}%\n")
        output_file.write("# Read Name, Aligned Read\n")
        
        # Get reference sequence
        ref_sequence = {}
        header = None
        with open(reference_file,"r") as f:
            for line in f:
                line = line.strip()
                
                # Skip empty lines
                if not line:
                    continue

                # Save reference name
                if line.startswith(">"):
                    # Start new header + sequence
                    header = line[1:]
                    seq_lines = []
                # Otherwise its sequence data
                else:
                    seq_lines.append(line)
                # Save header + sequence
                if header:
                    ref_sequence[header] = "".join(seq_lines)
        for header, sequence in ref_sequence.items():
            output_file.write(f"[REF] {header}\t{sequence}\n")

        # Decode Cigar String
        for read in data_out.keys():
            # Get data we pulled out from file
            read_name = data_out[read][0]
            ref_start = data_out[read][1]
            cigar = data_out[read][2]
            read_seq = data_out[read][3]
            # map_qual = data_out[read][4]

            # Bin for read alignment and initial read_position (default ref offset is +1)
            aligned_read = []
            read_position = 0
            ref_position = ref_start - 1 # SAM is 1-based

            # Handle unaligned sequences
            if cigar == "*":
                # No alignment
                aligned_read.append("(No Alignment Found)")
            else:
                # Parse CIGAR
                parsed_cigar = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
                
                # Add initial offset to read_position (+1 or more)
                aligned_read.append(" "*ref_position)

                # Loop over parsed_cigar to get global alignment
                for length, operator in parsed_cigar:
                    if operator in ["M", "=", "X"]:
                        # Match / Mismatch
                        aligned_read.append(read_seq[read_position:read_position+int(length)])
                        read_position += int(length)
                        ref_position += int(length)
                    elif operator in ["N", "D"]:
                        # Skipped or Deleted regions
                        aligned_read.append("-" * int(length))
                        ref_position += int(length)
                    elif operator == "I":
                        # Inserted Regions
                        aligned_read.append(read_seq[read_position:read_position+int(length)])
                        read_position += int(length)
                    elif operator == "S":
                        # Soft clipped regions (not aligned)
                        read_position += int(length)
                    else:
                        # Skip anything else, like hard clipping or unmapped reads
                        continue
            output_file.write(f"{read_name}\t{''.join(aligned_read)}\n")

    # # Delete source file
    # command = f"rm {output_directory}/all_reads.txt"
    # subprocess.run(command, shell=True, check=True)

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

def generate_mpileup(aligned_bam_path, reference_fasta_path, output_directory):
    # Generate new index for the reference fasta file for samtools mpileup if it doesn't exist
    try:
        with open(f"{output_directory}/ref_lib/{reference_fasta_path.split('/')[-1]}.fai", 'r') as f:
            print(f"[*] Reference fasta index found for samtools mpileup, skipping index generation.")
            pass
    except FileNotFoundError:
        print(f"[!] Reference fasta index not found for samtools mpileup, generating index...")
        command = f"samtools faidx {reference_fasta_path} && mv {reference_fasta_path}.fai {output_directory}/ref_lib/"
        subprocess.run(command, shell=True, check=True)
    
    # Generate mpileup file using samtools mpileup
    command = f"samtools mpileup -f {reference_fasta_path} -d 500000 {aligned_bam_path} > {output_directory}/{aligned_bam_path.split('/')[-1]}_mpileup.txt"
    subprocess.run(command, shell=True, check=True)
    return f"{output_directory}/{aligned_bam_path.split('/')[-1]}_mpileup.txt"
    
def mpileup_cleaner(read_bases):
    # This function processes the mpileup file to remove read starts
    cleaned_read = ""
    i = 0
    n = len(read_bases)
    while i < n:
        char = read_bases[i]

        # Skip read starts (indicated by '^' followed by a quality score)
        if char == '^':
            i += 2
            continue
        # Skip read ends (indicated by '$')
        if char == '$':
            i += 1
            continue
        
        # Skip Insertions and Deletions (indicated by '+' or '-' followed by the length of the indel and the inserted/deleted bases)
        if char == '+':
            i += 1
            length_str = ""

            # Read the length of the indel
            while i < n and read_bases[i].isdigit():
                length_str += read_bases[i]
                i += 1
            
            if length_str:
                indel_length = int(length_str)
                cleaned_read += char # Keep the indel symbol but not data in the cleaned read for counting later
                i += indel_length  # Skip the inserted bases
            continue
        # Deletions
        if char == '-':
            i += 1
            length_str = ""

            # Read the length of the indel
            while i < n and read_bases[i].isdigit():
                length_str += read_bases[i]
                i += 1
            
            if length_str:
                indel_length = int(length_str)
                cleaned_read += char # Keep the indel symbol but not data in the cleaned read for counting later
                i += indel_length  # Skip the deleted bases
            continue

        # Deletion placeholders should be removed from the read bases
        if char == '*':
            i += 1
            continue

        # Otherwise, add the character to the cleaned read bases
        cleaned_read += char
        i += 1
    
    return cleaned_read

def generate_mpileup_full_analysis(mpileup_file_path, output_directory):
    # Write a header to output file for the full mpileup analysis
    with open(f"{output_directory}/mpileup_full_analysis.txt", 'w') as out_f:
        out_f.write("Position\tPercent_Identical\tReference_Base\tPercent_A\tPercent_C\tPercent_G\tPercent_T\tPercent_Insertions\tPercent_Deletions\tDepth\n")

        # Generate a full analysis of the mpileup file, including percent identical to reference, base composition, and indel frequencies at each position
        with open(mpileup_file_path, 'r') as f:
                for line in f:
                    # Parse the mpileup line to extract relevant information
                    fields = line.strip().split('\t')
                    pos = fields[1]
                    ref_base = fields[2]
                    depth = fields[3]
                    read_bases = fields[4]
                    
                    # Clean the read bases to remove read starts, ends, and indel data beyond their presence
                    cleaned_read_bases = mpileup_cleaner(read_bases)
                    # Convert cleaned read bases exact base calls for counting
                    substituted_bases = cleaned_read_bases.replace('.', ref_base).replace(',', ref_base.lower())
                    
                    # Count the occurrences of each base and indel symbol in the cleaned read bases
                    base_counts = {
                        'A': substituted_bases.count('A') + substituted_bases.count('a'),
                        'C': substituted_bases.count('C') + substituted_bases.count('c'),
                        'G': substituted_bases.count('G') + substituted_bases.count('g'),
                        'T': substituted_bases.count('T') + substituted_bases.count('t'),
                        '+': substituted_bases.count('+'),  # Count insertions
                        '-': substituted_bases.count('-')   # Count deletions
                    }
                    # Convert counts to percentages
                    nt_count = sum(base_counts[base] for base in ['A', 'C', 'G', 'T'])
                    base_percentages = {base: (count / nt_count) * 100 for base, count in base_counts.items()}
                    # Calculate percent identical to reference
                    percent_identical = base_counts[ref_base.upper()] / nt_count * 100 if ref_base.upper() in base_counts else 0
                    # Output the analysis results for this position
                    out_f.write(f"{pos}\t{percent_identical:.2f}%\t{ref_base}\t{base_percentages['A']:.2f}%\t{base_percentages['C']:.2f}%\t{base_percentages['G']:.2f}%\t{base_percentages['T']:.2f}%\t{base_percentages['+']:.2f}%\t{base_percentages['-']:.2f}%\t{depth}\n")

def generate_mpileup_simple_analysis(mpileup_file_path, output_directory):
    # Write a header to output file for the simple mpileup analysis
    with open(f"{output_directory}/mpileup_simple_analysis.txt", 'w') as out_f:
        out_f.write("Position\tPercent_Identical\tReference_Base\tDepth\n")

        # Generate a simple analysis of the mpileup file, including only percent identical to reference and depth at each position
        with open(mpileup_file_path, 'r') as f:
            for line in f:
                # Parse the mpileup line to extract relevant information
                fields = line.strip().split('\t')
                pos = fields[1]
                ref_base = fields[2]
                depth = fields[3]
                read_bases = fields[4]
                
                # Clean the read bases to remove read starts, ends, and indel data beyond their presence
                cleaned_read_bases = mpileup_cleaner(read_bases)
                # Convert cleaned read bases exact base calls for counting
                substituted_bases = cleaned_read_bases.replace('.', ref_base).replace(',', ref_base.lower())
                
                # Count the occurrences of each base and indel symbol in the cleaned read bases
                base_counts = {
                    'A': substituted_bases.count('A') + substituted_bases.count('a'),
                    'C': substituted_bases.count('C') + substituted_bases.count('c'),
                    'G': substituted_bases.count('G') + substituted_bases.count('g'),
                    'T': substituted_bases.count('T') + substituted_bases.count('t'),
                    '+': substituted_bases.count('+'),  # Count insertions
                    '-': substituted_bases.count('-')   # Count deletions
                }
                # Convert counts to percentages
                nt_count = sum(base_counts[base] for base in ['A', 'C', 'G', 'T'])
                # Calculate percent identical to reference
                percent_identical = base_counts[ref_base.upper()] / nt_count * 100 if ref_base.upper() in base_counts else 0
                # Output the analysis results for this position
                out_f.write(f"{pos}\t{percent_identical:.2f}%\t{ref_base}\t{depth}\n")

def generate_mpileup_visualization(mpileup_file_path, output_directory):
    # Generate plots of the mpileup file - percent identical, percent mutation, indel frequencies across positions, and sequencing depth across positions
    output_data = []
    with open(mpileup_file_path, 'r') as f:
        for line in f:
            # Parse the mpileup line to extract relevant information
            fields = line.strip().split('\t')
            pos = fields[1]
            ref_base = fields[2]
            depth = fields[3]
            read_bases = fields[4]
            
            # Clean the read bases to remove read starts, ends, and indel data beyond their presence
            cleaned_read_bases = mpileup_cleaner(read_bases)
            # Convert cleaned read bases exact base calls for counting
            substituted_bases = cleaned_read_bases.replace('.', ref_base).replace(',', ref_base.lower())
            
            # Count the occurrences of each base and indel symbol in the cleaned read bases
            base_counts = {
                'A': substituted_bases.count('A') + substituted_bases.count('a'),
                'C': substituted_bases.count('C') + substituted_bases.count('c'),
                'G': substituted_bases.count('G') + substituted_bases.count('g'),
                'T': substituted_bases.count('T') + substituted_bases.count('t'),
                '+': substituted_bases.count('+'),  # Count insertions
                '-': substituted_bases.count('-')   # Count deletions
            }
            # Convert counts to percentages
            nt_count = sum(base_counts[base] for base in ['A', 'C', 'G', 'T'])
            # Calculate percent identical to reference
            percent_identical = base_counts[ref_base.upper()] / nt_count * 100 if ref_base.upper() in base_counts else 0
            output_data.append((pos, percent_identical, base_counts['-'], base_counts['+'], depth))
    
    # Separate the output data into lists for plotting
    pos = [int(data[0]) for data in output_data]
    percent_identical = [data[1] for data in output_data]
    percent_mutation = [100 - data[1] for data in output_data]
    indel_deletions = [data[2] for data in output_data]
    indel_insertions = [data[3] for data in output_data]
    depth = [int(data[4]) for data in output_data]

    # Generate and save the percent identical plot
    plt.plot(pos, percent_identical, color='blue')
    plt.title('Percent Identical to Reference Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Percent Identical (%)')
    plt.xlim(0, max(pos))
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(f"{output_directory}/mpileup_percent_identical_plot.png")
    plt.clf()
    # Generate and save the percent mutation plot
    plt.plot(pos, percent_mutation, color='red')
    plt.title('Percent Mutation from Reference Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Percent Mutation (%)')
    plt.xlim(0, max(pos))
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(f"{output_directory}/mpileup_percent_mutation_plot.png")
    plt.clf()
    # Generate and save the indel frequency plot
    plt.plot(pos, indel_insertions, color='green', label='Insertions')
    plt.plot(pos, indel_deletions, color='orange', label='Deletions')
    plt.title('Indel Frequencies Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Indel Count')
    plt.xlim(0, max(pos))
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_directory}/mpileup_indel_frequencies_plot.png")
    plt.clf()
    # Generate and save the sequencing depth plot
    plt.plot(pos, depth, color='purple')
    plt.title('Sequencing Depth Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.xlim(0, max(pos))
    plt.tight_layout()
    plt.savefig(f"{output_directory}/mpileup_depth_plot.png")
    plt.clf()

def main(analysis_params, reference_fasta_path, output_directory):
    # Track the time taken for analysis
    start_time = time.perf_counter()

    # Print samtools version
    command = f"samtools --version"
    subprocess.run(command, shell=True, check=True)

    # Run the analysis steps based on the specified parameters
    if analysis_params.get("do-alignment-stats", True):
        print("[*] Generating alignment statistics...")
        generate_alignment_stats(f"{output_directory}/aligned_reads.bam", output_directory)
        print("[*] Alignment statistics generated.")
    
    if analysis_params.get("do-alignment-visualization", True):
        print("[*] Generating alignment visualization...")
        try:
            generate_alignment_visualization(f"{output_directory}/aligned_reads.bam", reference_fasta_path, output_directory)
            print("[*] Alignment visualization generated.")
        except Exception as e:
            print(f"[!] Error generating alignment visualization: {e}")

    if analysis_params.get("do-alignment-score-plot", True):
        print("[*] Generating alignment score plot...")
        try:
            generate_alignment_score_plot(f"{output_directory}/aligned_reads.bam", output_directory)
            print("[*] Alignment score plot generated.")
        except Exception as e:
            print(f"[!] Error generating alignment score plot: {e}")

    if analysis_params.get("do-mpileup", True):
        print("[*] Generating mpileup file...")
        try:
            mpileup_file = generate_mpileup(f"{output_directory}/aligned_reads.bam", reference_fasta_path, output_directory)
            print("[*] Mpileup file generated.")
        except Exception as e:
            print(f"[!] Error generating mpileup file: {e}")

    if analysis_params.get("do-mpileup-fullanalysis", True):
        print("[*] Generating full mpileup analysis...")
        try:
            generate_mpileup_full_analysis(mpileup_file, output_directory)
            print("[*] Full mpileup analysis generated.")
        except Exception as e:
            print(f"[!] Error generating full mpileup analysis: {e}")

    if analysis_params.get("do-mpileup-simpleanalysis", True):
        print("[*] Generating simple mpileup analysis...")
        try:
            generate_mpileup_simple_analysis(mpileup_file, output_directory)
            print("[*] Simple mpileup analysis generated.")
        except Exception as e:
            print(f"[!] Error generating simple mpileup analysis: {e}")

    if analysis_params.get("do-mpileup-visualization", True):
        print("[*] Generating mpileup visualizations...")
        try:
            generate_mpileup_visualization(mpileup_file, output_directory)
            print("[*] Mpileup visualizations generated.")
        except Exception as e:
            print(f"[!] Error generating mpileup visualizations: {e}")

    # Track the time taken for analysis
    end_time = time.perf_counter()
    print(f"[*] Analysis completed in {end_time - start_time:.2f} seconds.")
