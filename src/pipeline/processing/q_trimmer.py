import time
from Bio import SeqIO

def rolling_trim(sequence_file, window_size=5, quality_threshold=20):
    # Trim sequences based on a rolling window quality score
    # Q threshold of 20 is ~ 99% base call accuracy (-10*log10(0.01))
    # Window size of 5 is standard, 
    #   shorter windows increase sensitivity but may be overly aggressive, 
    #   whereas longer windows are generally more permissive but may miss low-quality regions
    trimmed_sequences = []
    for record in SeqIO.parse(sequence_file, "fastq"):
        # Get the quality scores for the sequence
        q_scores = record.letter_annotations["phred-quality"]
        # For all bases/scores in the sequence, check if the quality score is below the threshold
        for i in range(len(q_scores) - window_size + 1):
            # Get the quality scores for the current window
            window = q_scores[i:i + window_size]
            # Check if the average quality score in the window is below the threshold
            if sum(window) / window_size < quality_threshold:
                # Replace low-quality bases with 'N'
                for j in range(i, i + window_size):
                    record.seq[j] = "N"
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
        for i, score in enumerate(q_scores):
            if score < quality_threshold:
                # Replace low-quality bases with 'N'
                record.seq[i] = "N"
        # Append the trimmed record to the list of trimmed sequences
        trimmed_sequences.append(record)
    return trimmed_sequences

def main(processing_parameters, mode, mode_config, output_directory):
    # Start timer for benchmarking
    start_time = time.perf_counter()

    # Perform quality trimming on the input sequences based on the specified parameters
    if processing_parameters["qtrimming-method"] == "rolling-trim":
        print("[*] Using rolling-trim method for quality trimming...")
        for sequence_file in mode_config[mode].values():
            trimmed_sequences = rolling_trim(sequence_file, processing_parameters["trimming-window-size"], processing_parameters["trimming-quality-threshold"])
    elif processing_parameters["qtrimming-method"] == "simple-trim":
        print("[*] Using simple-trim method for quality trimming...")
        for sequence_file in mode_config[mode].values():
            trimmed_sequences = simple_trim(sequence_file, processing_parameters["trimming-quality-threshold"])

    # Write the trimmed sequences to a new FASTQ file in the output directory
    output_file = f"{output_directory}/trimmed_R1.fastq"
    SeqIO.write(trimmed_sequences, output_file, "fastq")
    print(f"[*] Quality trimming complete. Trimmed sequences written to {output_file}.")

    # End timer for benchmarking
    end_time = time.perf_counter()
    print(f"[*] Quality trimming completed in {end_time - start_time:.2f} seconds.")