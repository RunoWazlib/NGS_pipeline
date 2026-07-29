import subprocess, shutil
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
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
        
        # Skip Insertions (indicated by '+' followed by the length of the indel and the inserted/deleted bases)
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
        out_f.write("Position\tIdentical_Calls\tReference_Base\tCount_A\tCount_C\tCount_G\tCount_T\tCount_Insertion\tCount_Deletion\tMutation_Count\tDepth\n")

        # Generate a full analysis of the mpileup file, including percent identical to reference, base composition, and indel frequencies at each position
        with open(mpileup_file_path, 'r') as f:
                for line in f:
                    # Parse the mpileup line to extract relevant information
                    fields = line.strip().split('\t')
                    pos = fields[1]
                    ref_base = fields[2]
                    depth = int(fields[3])
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
                    # Normalize by depth: read coverage / position
                    try:
                        base_rates = {base: (count) for base, count in base_counts.items()}
                        
                        # count calls identical to reference
                        count_identical = (base_counts[ref_base.upper()]) if ref_base.upper() in base_counts else 0
                        
                        # Sum mismatch + insertion + deletion counts (anything not reference base), then normalize by depth
                        # Since these mutations are not mutually exclusive in a read, this is not a mutation frequency (proportion of mutant reads / coverage) - this is a per position event rate (mutant event / coverage)
                        # This mutation rate must be <= 2 since each position can at maximum have a mismatch and an indel (either insertion or deletion) per read
                        # in essence this is averaging the mutation count of every read at every position by the coverage of that position
                        mutation_rate = (sum(base_counts[base] for base in base_counts.keys() if base != ref_base.upper()) / depth)
                    # If there's no depth, we don't know anything!
                    except ZeroDivisionError:
                        base_rates = {"A":0,"C":0,"G":0,"T":0,"+":0,"-":0}
                        count_identical = 0
                        mutation_rate = 0

                    # Output the analysis results for this position
                    # TODO - Depth normalization, add in :.2f 
                    out_f.write(f"{pos}\t{count_identical}\t{ref_base}\t{base_rates['A']}\t{base_rates['C']}\t{base_rates['G']}\t{base_rates['T']}\t{base_rates['+']}\t{base_rates['-']}\t{mutation_rate}\t{depth}\n")

def generate_mpileup_simple_analysis(mpileup_file_path, output_directory):
    # Write a header to output file for the simple mpileup analysis
    with open(f"{output_directory}/Simple_analysis.txt", 'w') as out_f:
        out_f.write("Position\tRate_Identical\tReference_Base\tDepth\n")

        # Generate a simple analysis of the mpileup file, including only percent identical to reference and depth at each position
        with open(mpileup_file_path, 'r') as f:
            for line in f:
                # Parse the mpileup line to extract relevant information
                fields = line.strip().split('\t')
                pos = fields[1]
                ref_base = fields[2]
                depth = int(fields[3])
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
                try:
                    # Calculate percent identical to reference
                    percent_identical = base_counts[ref_base.upper()] / depth if ref_base.upper() in base_counts else 0
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
            pos = int(fields[1])
            ref_base = fields[2]
            depth = int(fields[3])
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
            # Normalize by depth: read coverage / position ("*" do contribute to read coverage as it was an opportunity to call a base, but was deleted according to alignment)
            try:
                base_rates = {base: (count / depth) for base, count in base_counts.items()}
                
                # count calls identical to reference
                rate_identical = (base_counts[ref_base.upper()] / depth) if ref_base.upper() in base_counts else 0
                
                # Sum mismatch + insertion + deletion counts (anything not reference base), then normalize by depth
                # Since these mutations are not mutually exclusive in a read, this is not a mutation frequency (proportion of mutant reads / coverage) - this is a per position event rate (mutant event / coverage)
                # This mutation rate must be <= 2 since each position can at maximum have a mismatch and an indel (either insertion or deletion) per read
                # in essence this is averaging the mutation count of every read at every position by the coverage of that position
                mutation_rate = (sum(base_counts[base] for base in base_counts.keys() if base != ref_base.upper()) / depth)
            # If there's no depth, we don't know anything!
            except ZeroDivisionError:
                base_rates = {"A":0,"C":0,"G":0,"T":0,"+":0,"-":0}
                rate_identical = 0
                mutation_rate = 0
            # Append output to list
            output_data.append((pos, rate_identical, mutation_rate, base_rates['-'], base_rates['+'], depth))
    
    # Separate the output data into lists for plotting
    pos = [data[0] for data in output_data]
    rate_identical = [data[1] for data in output_data]
    mutation_rate = [data[2] for data in output_data]
    indel_deletions = [data[3] for data in output_data]
    indel_insertions = [data[4] for data in output_data]
    depth = [data[5] for data in output_data]

    # Generate and save the percent identical plot
    plt.plot(pos, rate_identical, color='blue')
    plt.title('Percent Identical to Reference Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Rate Identical (identical base per read)')
    plt.xlim(0, max(pos))
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(f"{output_directory}/Rate_identical_plot.png")
    plt.clf()
    # Generate and save the percent mutation plot
    plt.plot(pos, mutation_rate, color='red')
    plt.title('Percent Mutation from Reference Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Mutation rate (mutations per read)')
    plt.xlim(0, max(pos))
    plt.ylim(0, 2)
    plt.tight_layout()
    plt.savefig(f"{output_directory}/Rate_mutation_plot.png")
    plt.clf()
    # Generate and save the indel frequency plot
    plt.plot(pos, indel_insertions, color='green', label='Insertions')
    plt.plot(pos, indel_deletions, color='orange', label='Deletions')
    plt.title('Indel Rates Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Indel Rate (indel start per read)')
    plt.xlim(0, max(pos))
    plt.legend()
    # Plot tick formatting
    plt.gca().yaxis.set_major_locator(mticker.MaxNLocator(5)) # 5 major ticks maximum
    plt.tight_layout()
    plt.savefig(f"{output_directory}/Indel_rates_plot.png")
    plt.clf()
    # Generate and save the sequencing depth plot
    plt.plot(pos, depth, color='purple')
    plt.title('Sequencing Coverage Across Positions')
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.xlim(0, max(pos))
    plt.tight_layout()
    plt.savefig(f"{output_directory}/Depth_plot.png")
    plt.clf()

def generate_indel_analysis(mpileup_file_path, output_directory):
    with open(f"{output_directory}/Indel_analysis.txt", "w") as out_f:
        # Read in mpileup
        with open(mpileup_file_path, "r") as in_f:
            for line in in_f:
                # Parse mpileup for indels
                fields = line.strip().split('\t')
                pos = fields[1]
                ref_base = fields[2]
                read_bases = fields[4]

                # Catch indels and mismatches
                i = 0
                n = len(read_bases)
                out_f.write(f"Position: {pos}\n")
                while i < n:
                    char = read_bases[i]
                    # Found Insertion (indicated by '+' followed by the length of the indel and the inserted/deleted bases)
                    if char == '+':
                        i += 1
                        length_str = ""

                        # Read the length of the indel
                        while i < n and read_bases[i].isdigit():
                            length_str += read_bases[i]
                            i += 1
                        if length_str:
                            indel_length = int(length_str)
                            indel = ""
                            while i < n and i < i + indel_length:
                                indel += read_bases[i]
                                i += 1

                            # Write indel to output file
                            out_f.write(f"[+] Insertion\t+{indel_length}\t{indel}\n")
                            # Next base
                            i += 1
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
                            indel = ""
                            while i < n and i < i + indel_length:
                                indel += read_bases[i]
                                i += 1

                            # Write indel to output file
                            out_f.write(f"[-] Deletion\t+{indel_length}\t{indel}\n")
                            # Next base
                            i += 1
                        continue
                    
                    # Found mismatch
                    if char.upper() in "ACGT" and char.upper() != ref_base.upper():
                        # Write mismatch to output file
                        out_f.write(f"[!] Mismatch\t{char} != {ref_base}\n")
                        # Next base
                        i += 1
                        continue

                    # If you didn't catch anything, keep going!
                    i += 1
