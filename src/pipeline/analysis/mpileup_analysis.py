import subprocess, shutil
from pathlib import Path
import matplotlib.pyplot as plt

def generate_mpileup(aligned_bam_path, reference_fasta_path, output_directory):
    # Generate new index for the reference fasta file for samtools mpileup if it doesn't exist
    try:
        open(f"{Path(output_directory).parent}/ref_lib/{reference_fasta_path.split('/')[-1]}.fai", 'r')
        print(f"[*] Reference fasta index found for samtools mpileup, skipping index generation.")
        pass
    except FileNotFoundError:
        print("[!] Reference fasta index not found for samtools mpileup, generating index...")
        command = f"samtools faidx {reference_fasta_path}"
        subprocess.run(command, shell=True, check=True)
        
    # Generate mpileup file using samtools mpileup
    command = f"samtools mpileup -f {reference_fasta_path} -d 500000 {aligned_bam_path} > {output_directory}/{aligned_bam_path.split('/')[-1]}_mpileup.txt"
    subprocess.run(command, shell=True)

    # Move this new index into the "output/ref_lib/"
    # TODO - diagnose why we can't move this file around before mpileup is done; maybe move the unindexed file in here too and point to ref_lib during mpileup
    shutil.move(Path(f"{reference_fasta_path}.fai"), Path(f"{Path(output_directory).parent}/ref_lib/ref.fai"))
    
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
    with open(f"{output_directory}/Full_analysis.txt", 'w') as out_f:
        out_f.write("Position\tPercent_Identical\tReference_Base\tPercent_A\tPercent_C\tPercent_G\tPercent_T\tPercent_Insertions\tPercent_Deletions\tDepth\n")

        # Generate a full analysis of the mpileup file, including percent identical to reference, base composition, and indel frequencies at each position
        with open(mpileup_file_path, 'r') as f:
                for line in f:
                    # Parse the mpileup line to extract relevant information
                    fields = line.strip().split('\t')
                    print(fields)
                    pos = fields[1]
                    ref_base = fields[2]
                    depth = fields[3]
                    read_bases = fields[4]
                    
                    # Clean the read bases to remove read starts, ends, and indel data beyond their presence
                    cleaned_read_bases = mpileup_cleaner(read_bases)
                    # Convert cleaned read bases exact base calls for counting
                    substituted_bases = cleaned_read_bases.replace('.', ref_base).replace(',', ref_base.lower())
                    print(substituted_bases)
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
                    try:
                        base_percentages = {base: (count / nt_count) * 100 for base, count in base_counts.items()}
                        # Calculate percent identical to reference
                        percent_identical = base_counts[ref_base.upper()] / nt_count * 100 if ref_base.upper() in base_counts else 0
                    except ZeroDivisionError:
                        base_percentages = {"A":0,"C":0,"G":0,"T":0,"+":0,"-":0}
                        percent_identical = 0.00

                    # Output the analysis results for this position
                    out_f.write(f"{pos}\t{percent_identical:.2f}%\t{ref_base}\t{base_percentages['A']:.2f}%\t{base_percentages['C']:.2f}%\t{base_percentages['G']:.2f}%\t{base_percentages['T']:.2f}%\t{base_percentages['+']:.2f}%\t{base_percentages['-']:.2f}%\t{depth}\n")

def generate_mpileup_simple_analysis(mpileup_file_path, output_directory):
    # Write a header to output file for the simple mpileup analysis
    with open(f"{output_directory}/Simple_analysis.txt", 'w') as out_f:
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
                try:
                    # Calculate percent identical to reference
                    percent_identical = base_counts[ref_base.upper()] / nt_count * 100 if ref_base.upper() in base_counts else 0
                except ZeroDivisionError:
                    percent_identical = 0.00
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
            try:
                # Calculate percent identical to reference
                percent_identical = base_counts[ref_base.upper()] / nt_count * 100 if ref_base.upper() in base_counts else 0
            except ZeroDivisionError:
                percent_identical = 0.00
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
    plt.savefig(f"{output_directory}/Percent_identical_plot.png")
    plt.clf()
    # Generate and save the percent mutation plot
    plt.plot(pos, percent_mutation, color='red')
    plt.title('Percent Mutation from Reference Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Percent Mutation (%)')
    plt.xlim(0, max(pos))
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(f"{output_directory}/Percent_mutation_plot.png")
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
    plt.savefig(f"{output_directory}/Indel_frequencies_plot.png")
    plt.clf()
    # Generate and save the sequencing depth plot
    plt.plot(pos, depth, color='purple')
    plt.title('Sequencing Depth Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.xlim(0, max(pos))
    plt.tight_layout()
    plt.savefig(f"{output_directory}/Depth_plot.png")
    plt.clf()
