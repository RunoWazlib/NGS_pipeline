#!/bin/bash
#SBATCH --job-name=ngs_pipeline # Job name
#SBATCH --output=ngs_pipeline_%j.log # Standard output file, %j will be replaced with the job ID
#SBATCH --time=02:00:00 # Adjust as needed hh:mm:ss
#SBATCH --cpus-per-task=1 
#SBATCH --mem=16g # Adjust as needed
#SBATCH --partition=mb # Adjust as needed for your desired cluster

# This is a template for submitting the NGS pipeline to the cluster
module load conda
conda activate seq

# Add the path to your `ngs_pipeline` directory to $PATH
export PATH="$WORK/ngs_pipeline/":$PATH

# Run pipeline over all config (.json) files in the current directory
for config_file in *.json; do
    echo "Processing $config_file"
    ngs_pipeline.py --config "$config_file"
done

# Review ngs_pipeline_%j.log for all typically printed outputs, like total evaluation time