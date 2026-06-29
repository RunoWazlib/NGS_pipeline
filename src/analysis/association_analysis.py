import subprocess, re, time
class ReadAlignment:
    def __init__(self, aligned_bam_path, reference_file, output_directory):
        self.aligned_bam_path = aligned_bam_path
        self.reference_file = reference_file
        self.output_directory = output_directory
        self.aligned_reads = {}
        self.mutation_map = {}
        self.contingency_tables = {}
        self.chi_squares = {}
        self.total_reads_count = 0
        self.aligned_reads_count = 0

    def get_nucleotide_sequences(self):
        # Check if the "all_reads.txt" already exists
        try:
            with open(f"{self.output_directory}/all_reads.txt", 'r') as f:
                print(f"[*] All reads file found in {self.output_directory}, skipping sequence extraction.")
        except FileNotFoundError:
            print(f"[!] All reads file not found - Extracting sequences from '{self.aligned_bam_path}'...")
            # Get python-parsable file of all sequences
            command = f"samtools view {self.aligned_bam_path} > {self.output_directory}/all_reads.txt"
            subprocess.run(command, shell=True, check=True)

        # Start parsing "all_reads.txt" to get aligned reads and reference sequence
        total_read_count = 0
        aligned_read_count = 0
        data_out = {}
        with open(f"{self.output_directory}/all_reads.txt","r") as source_file:
            for line in source_file:
                # Count total number of reads (including unaligned)
                total_read_count += 1
                # Get entry data
                fields = line.strip().split("\t")
                read_name = fields[0]
                start_ref = int(fields[3])
                map_qual = int(fields[4])
                cigar_str = fields[5]
                read_seq = fields[9]
                # If read is unaligned, move on
                if cigar_str == "*":
                    data_out[total_read_count] = [read_name, start_ref, cigar_str, read_seq, map_qual]
                else:
                    aligned_read_count += 1
                    data_out[total_read_count] = [read_name, start_ref, cigar_str, read_seq, map_qual]
            
        # Define read_counts
        self.total_reads_count = total_read_count
        self.aligned_reads_count = aligned_read_count

        # Get reference sequence
        ref_sequence = {}
        header = None
        with open(self.reference_file,"r") as f:
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

                if header:
                    ref_sequence[header] = "".join(seq_lines)
            # Save last reference sequence to aligned reads dict
            self.aligned_reads["REF"] = ref_sequence[header]

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
                pass
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
            # Save aligned read
            self.aligned_reads[read_name] = "".join(aligned_read)
    
    def generate_mutation_map(self):
        ref_seq = self.aligned_reads["REF"]
        # Build a mutation map for each position in a read compared to the reference sequence (1 if mutated, 0 if not)
        for read_name in self.aligned_reads.keys():
            if read_name == "REF":
                # Skip reference sequence
                continue
            read_seq = self.aligned_reads[read_name]
            if len(read_seq) != len(ref_seq):
                # Skip reads that don't have the same length as the reference sequence (i.e. unaligned reads / truncated reads)
                continue
            # Determine if read is mutated or not at each position
            self.mutation_map[read_name] = []
            for i in range(len(ref_seq)):
                ref_nucleotide = ref_seq[i]
                read_nucleotide = read_seq[i]
                if read_nucleotide == ref_nucleotide:
                    # Read matches reference at this position
                    mutation_status = 0
                else:
                    mutation_status = 1
                self.mutation_map[read_name].append(mutation_status)
        return self.mutation_map

    def statistical_association_analysis(self):
        # Perform statistical association analysis on aligned reads as proposed by Homan et.al. 2014
        # Get mutation map for each read with respect to reference sequence
        self.generate_mutation_map()
        # Generate contingency table for each position in the reference sequence with respect to another position
        ref_seq = self.aligned_reads["REF"]
        for i in range(len(ref_seq)):
            for j in range(i+1, len(ref_seq)):
                # Get counts for contingency table
                count_00 = 0 # Both positions not mutated, i.e. A
                count_01 = 0 # Position i not mutated, position j mutated i.e. C
                count_10 = 0 # Position i mutated, position j not mutated i.e. B
                count_11 = 0 # Both positions mutated i.e. D
                for read_name in self.mutation_map.keys():
                    pos_i_mutated = self.mutation_map[read_name][i]
                    pos_j_mutated = self.mutation_map[read_name][j]
                    if pos_i_mutated == 0 and pos_j_mutated == 0:
                        count_00 += 1
                    elif pos_i_mutated == 0 and pos_j_mutated == 1:
                        count_01 += 1
                    elif pos_i_mutated == 1 and pos_j_mutated == 0:
                        count_10 += 1
                    elif pos_i_mutated == 1 and pos_j_mutated == 1:
                        count_11 += 1
                # [[A,C], [B,D]]
                self.contingency_tables[(i,j)] = [[count_00, count_01], [count_10, count_11]]

        # Perform chi-square analysis
        for positions, table in self.contingency_tables.items():
            count_00, count_01 = table[0]
            count_10, count_11 = table[1]
            # Calculate chi-square statistic with Yates correction
            numerator = (abs(count_00 * count_11 - count_01 * count_10)-0.5)**2
            # # Calculate chi-square statistic under normal, Pearson chi-square
            # numerator = (count_00 * count_11 - count_01 * count_10)**2
            denominator = (count_00 + count_01) * (count_10 + count_11) * (count_00 + count_10) * (count_01 + count_11)
            if denominator == 0:
                print(f"[!] Warning: Zero denominator for positions {positions}.")
                chi_square = 0
            else:
                chi_square = numerator / denominator
            self.chi_squares[positions] = chi_square
        return self.chi_squares
    
def main(analysis_params, reference_file, output_directory):
    # Track the time taken for association analysis
    start_time = time.perf_counter()
    # Does config include association analysis?
    if analysis_params.get("do-association-analysis", True):
        print("[*] Generating association analysis...")
        # Create an instance of the ReadAlignment class
        read_alignment = ReadAlignment(f"{output_directory}/aligned_reads.bam", reference_file, output_directory)
        # Get nucleotide sequences from aligned BAM file and reference file
        read_alignment.get_nucleotide_sequences()
        # Perform statistical association analysis on aligned reads
        chi_square_results = read_alignment.statistical_association_analysis()
        # Output results to file
        with open(f"{output_directory}/association_analysis_summary.txt", "w") as f:
            f.write(f"Total Reads: {read_alignment.total_reads_count}\n")
            f.write(f"Aligned Reads: {read_alignment.aligned_reads_count}\n")
            f.write(f"Full Length Reads: {len(read_alignment.mutation_map)}\n")
            significant_associations = 0
            for positions, chi_square in chi_square_results.items():
                if chi_square > 3.841:  # Chi-square threshold for p < 0.05, 1 degree of freedom
                    f.write(f"Positions {positions}: Chi-Square = {chi_square:.4f}\n")
                    significant_associations += 1
            if significant_associations == 0:
                f.write("No significant associations found (p < 0.05, 1 df).")

        # Output full chi-square results to a TSV file   
        with open(f"{output_directory}/association_analysis_results.tsv", "w") as f:
            f.write("Position1\tPosition2\tChi-Square\n")
            for positions, chi_square in chi_square_results.items():
                f.write(f"{positions[0]}\t{positions[1]}\t{chi_square:.4f}\n")
    # Print total time taken for association analysis
    end_time = time.perf_counter()
    print(f"[-] Association analysis completed in {end_time - start_time:.2f} seconds.")