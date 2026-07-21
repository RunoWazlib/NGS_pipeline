import time, gzip
from Bio import SeqIO
from Bio.Seq import Seq

def rolling_trim(sequence_file, window_size=5, quality_threshold=20):
    # Trim sequences based on a rolling window quality score
    # Q threshold of 20 is ~ 99% base call accuracy (-10*log10(0.01))
    # Window size of 5 is standard, 
    #   shorter windows increase sensitivity but may be overly aggressive, 
    #   whereas longer windows are generally more permissive but may miss low-quality regions
    trimmed_sequences = []
    for record in SeqIO.parse(sequence_file, "fastq"):
        # Get the quality scores for the sequence
        q_scores = record.letter_annotations["phred_quality"]
        # For all bases/scores in the sequence, check if the quality score is below the threshold
        for i in range(len(q_scores) - window_size + 1):
            # Get the quality scores for the current window
            window = q_scores[i:i + window_size]
            # Check if the average quality score in the window is below the threshold
            if sum(window) / window_size < quality_threshold:
                # Replace low-quality bases with 'N'
                # record.seq is immutable, convert to list, modify, and reassign
                seq_list = list(str(record.seq))
                for j in range(i, i + window_size):
                    seq_list[j] = "N"
                record.seq = Seq(''.join(seq_list))
        # Append the trimmed record to the list of trimmed sequences
        trimmed_sequences.append(record)
    return trimmed_sequences

def simple_trim(sequence_file, quality_threshold=20):
    # Trim sequences based on a simple quality threshold
    # Q threshold of 20 is ~ 99% base call accuracy (-10*log10(0.01))
    trimmed_sequences = []
    # For all sequences...
    for record in SeqIO.parse(sequence_file, "fastq"):
        # Get the quality scores for the sequence
        q_scores = record.letter_annotations["phred_quality"]
        # For all bases/scores in the sequence, check if the quality score is below the threshold
        # record.seq is immutable, convert to list once and modify
        seq_list = list(str(record.seq))
        for i, score in enumerate(q_scores):
            if score < quality_threshold:
                # Replace low-quality bases with 'N'
                seq_list[i] = "N"
        record.seq = Seq(''.join(seq_list))
        # Append the trimmed record to the list of trimmed sequences
        trimmed_sequences.append(record)
    return trimmed_sequences

def main(processing_parameters, mode, mode_config, output_directory):
    # Start timer for benchmarking
    start_time = time.perf_counter()
    output_files = {}
    # Perform quality trimming on the input sequences based on the specified parameters
    if processing_parameters["qtrimming-method"] == "rolling-trim":
        print("[*] Using rolling-trim method for quality trimming...")
        for key, sequence_file in mode_config[mode].items():
            trimmed_sequences = rolling_trim(gzip.open(sequence_file, "rt", encoding="utf-8"), processing_parameters["trimming-window-size"], processing_parameters["trimming-quality-threshold"])
            
            # Write the trimmed sequences to a new FASTQ file in the output directory
            output_file = f"{output_directory}/trimmed_{sequence_file.split('/')[-1]}"
            with gzip.open(output_file, "wt", encoding="utf-8") as out_f:
                SeqIO.write(trimmed_sequences, out_f, "fastq")
            # Add to out file list
            output_files[key] = output_file
            print(f"[*] Quality trimming complete. Trimmed sequences written to {output_file}.")
    
    elif processing_parameters["qtrimming-method"] == "simple-trim":
        print("[*] Using simple-trim method for quality trimming...")
        for key, sequence_file in mode_config[mode].items():
            trimmed_sequences = simple_trim(gzip.open(sequence_file, "rt",encoding="utf-8"), processing_parameters["trimming-quality-threshold"])
            
            # Write the trimmed sequences to a new FASTQ file in the output directory
            output_file = f"{output_directory}/trimmed_{sequence_file.split('/')[-1]}"
            with gzip.open(output_file, "wt", encoding="utf-8") as out_f:
                SeqIO.write(trimmed_sequences, out_f, "fastq")
            # Add to out file list
            output_files[key] = output_file
            print(f"[*] Quality trimming complete. Trimmed sequences written to {output_file}.")

    # End timer for benchmarking
    end_time = time.perf_counter()
    print(f"[*] Quality trimming completed in {end_time - start_time:.2f} seconds.")
    return output_files