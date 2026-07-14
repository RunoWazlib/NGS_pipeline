import subprocess, re

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