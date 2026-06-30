from pipeline.main import load_config
from pipeline.config_builder import generate_json_data_paired, generate_json_data_merged

import os, json, subprocess, pytest, shutil
from pathlib import Path

@pytest.fixture
def successful_config(tmp_path):
    """
    Fixture to create a successful configuration file for testing purposes.
    """
    # Copy necessary test files to the temporary directory
    source_dir = f"{Path(__file__).parent}/test_data"
    
    fileCopys = {
        f"{source_dir}/x98_query.fasta": f"{tmp_path}/reference.fasta",
        f"{source_dir}/X98-w-Mg_R1_001.fastq.gz": f"{tmp_path}/sample1_R1.fastq.gz",
        f"{source_dir}/X98-w-Mg_R2_001.fastq.gz": f"{tmp_path}/sample1_R2.fastq.gz"
    }
    
    for source, target in fileCopys.items():
        shutil.copy(source, target)

    # generate config data for paired-end mode
    config_data = {
        "mode": "paired-end-mode",
        "paired-end-mode":{
            "R1": f"{tmp_path}/sample1_R1.fastq.gz",
            "R2": f"{tmp_path}/sample1_R2.fastq.gz"
        },
        "reference-fasta": f"{tmp_path}/reference.fasta",
        "output-directory": str(tmp_path / "output"),
        "analysis-parameters": {
            "do-benchmarks": True,
            "do-alignment": True,
            "do-analysis": True
        }
    }
    config_file = f"{tmp_path}/config.json"
    with open(config_file, 'w') as f:
        json.dump(config_data, f)
    return config_file

class TestBasicInitialization:
    def test_paired_config_generation(self, tmp_path):
        """This test checks out "config_builder.py" script to ensure that sequencing configs are generated correctly for both paired and merged data

        Args:
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the raw fastq.gz files
        """
        # Change to the temporary directory for testing
        os.chdir(tmp_path)

        # Test paired-end config generation
        reference_fasta_path = "reference.fasta"
        output_directory_name = "output"
        sample_r1 = f"{tmp_path}/sample1_R1.fastq.gz"
        sample_r2 = f"{tmp_path}/sample1_R2.fastq.gz"
        with open(sample_r1, 'w') as f:
            f.write("dummy R1 content")
        with open(sample_r2, 'w') as f:
            f.write("dummy R2 content")
        
        paired_config = generate_json_data_paired(reference_fasta_path, output_directory_name)
        
        assert isinstance(paired_config, dict) # ensure output is a dict object for JSON serialization
        assert sample_r1 in paired_config["sample1"]["paired-end-mode"]["R1"] # ensure that the R1 file is correctly included in the paired-end config
        assert sample_r2 in paired_config["sample1"]["paired-end-mode"]["R2"] # ensure that the R2 file is correctly included in the paired-end config
        assert paired_config["sample1"]["mode"] == "paired-end-mode" # ensure that the mode is correctly set to "paired-end-mode" for the sample
        assert paired_config["sample1"]["reference-fasta"] == reference_fasta_path # ensure that the reference fasta path is correctly set in the config
        assert paired_config["sample1"]["output-directory"] == f"sample1_{output_directory_name}" # ensure that the output directory is correctly set in the config

    def test_merged_config_generation(self, tmp_path):
        """This test checks out "config_builder.py" script to ensure that sequencing configs are generated correctly for merged data
        Args:
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the raw fastq.gz files
        """
        
        # Change to the temporary directory for testing
        os.chdir(tmp_path)

        # Test merged config generation
        reference_fasta_path = "reference.fasta"
        output_directory_name = "output"
        sample_r1 = f"{tmp_path}/sample1_R1.fastq.gz"
        sample_r2 = f"{tmp_path}/sample1_R2.fastq.gz"
        with open(sample_r1, 'w') as f:
            f.write("dummy R1 content")
        with open(sample_r2, 'w') as f:
            f.write("dummy R2 content")

        # Test merged config generation
        merged_config = generate_json_data_merged(reference_fasta_path, output_directory_name)
        
        assert isinstance(merged_config, dict) # ensure output is a dict object for JSON serialization
        assert sample_r1 in merged_config["sample1"]["merged-mode"]["R1"] # ensure that the R1 file is correctly included in the merged config
        assert "R2" not in merged_config["sample1"]["merged-mode"].keys() # ensure that the R2 file is not included in the merged config
        assert merged_config["sample1"]["mode"] == "merged-mode" # ensure that the mode is correctly set to "merged-mode" for the sample
        assert merged_config["sample1"]["reference-fasta"] == reference_fasta_path # ensure that the reference fasta path is correctly set in the config
        assert merged_config["sample1"]["output-directory"] == f"sample1_{output_directory_name}" # ensure that the output directory is correctly set in the config

    def test_load_config(self, tmp_path):
        """This test confirms that the "load_config" method correctly loads intended config information - essentially verifying the JSON lib dump + load methods

        Args:
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the config file
        """
        # Create a temporary config file
        config_data = {
            "mode": "paired-end-mode",
            "reference-fasta": "reference.fasta",
            "output-directory": "output",
            "analysis-parameters": {
                "do-benchmarks": True,
                "do-alignment": True,
                "do-analysis": True
            }
        }
        
        config_file = f"{tmp_path}/config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        loaded_config = load_config(str(config_file))
        
        assert loaded_config == config_data

    def test_config_builder(self, tmp_path):
        """This test checks out "config_builder" link in python build to ensure that script is called correctly

        Args:
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the raw fastq.gz files
        """
        # Change to the temporary directory for testing
        os.chdir(tmp_path)

        # Test paired-end config generation
        reference_fasta_path = "reference.fasta"
        output_directory_name = "output"
        sample_r1 = f"{tmp_path}/sample1_R1.fastq.gz"
        sample_r2 = f"{tmp_path}/sample1_R2.fastq.gz"
        with open(sample_r1, 'w') as f:
            f.write("dummy R1 content")
        with open(sample_r2, 'w') as f:
            f.write("dummy R2 content")
        command = f"config-builder"
        subprocess.run(command, shell=True, check=True)

        assert Path(f"{tmp_path}/generated_config_sample1.json").exists()

    def test_output_dir_generation(self, tmp_path):
        """This test verifies that a. ngs_driver can read in a config file and b. that output directory is correctly generated

        Args:
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the config file
        """
        
        # Create a temporary config file with an output directory
        config_data = {
            "mode": "paired-end-mode",
            "reference-fasta": "reference.fasta",
            "output-directory": str(tmp_path / "output"),
            "analysis-parameters": {
                "do-benchmarks": False,
                "do-alignment": False,
                "do-analysis": False
            }
        }
        
        config_file = f"{tmp_path}/config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        # Run ngs_driver to create the output directory
        command = f"ngs-pipeline --config {config_file}"
        subprocess.run(command, shell=True, check=True)
        
        assert os.path.exists(tmp_path / "output") # ensure that the output directory is created
    
    def test_config_flags_disabled(self, tmp_path):
        """This test verifies that the flags in the config file are correctly interpreted by ngs_driver

        Args:
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the config file
        """
        
        # Create a temporary config file with all flags set to False
        config_data = {
            "mode": "paired-end-mode",
            "reference-fasta": "reference.fasta",
            "output-directory": str(tmp_path / "output"),
            "analysis-parameters": {
                "do-benchmarks": False,
                "do-alignment": False,
                "do-analysis": False
            }
        }
        
        config_file = tmp_path / "config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        # Run ngs_driver to check that no analysis is performed
        command = "ngs-pipeline --config " + str(config_file)
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        assert "[*] Skipping fastqc benchmarks as per configuration." in result.stdout
        assert "[*] Skipping alignment as per configuration." in result.stdout
        assert "[*] Skipping analysis as per configuration." in result.stdout

@pytest.mark.skip(reason="incomplete test")
class TestBasicBenchmarks:
    def test_benchmark_execution(self, successful_config):
        """This test verifies that the benchmarking step is executed when the corresponding flag is set to True in the config file

        Args:
            successful_config (_type_): pytest fixture for a successful configuration file
        """
        
        # Run ngs_driver to check that benchmarking is performed
        command = f"ngs-pipeline --config {successful_config}"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        assert "[*] Starting fastqc benchmarks..." in result.stdout

@pytest.mark.skip(reason="incomplete test")
class TestBasicAlignment:
    def test_alignment_execution(self, successful_config):
        """This test verifies that the alignment step is executed when the corresponding flag is set to True in the config file

        Args:
            successful_config (_type_): pytest fixture for a successful configuration file
        """
        
        # Run ngs_driver to check that alignment is performed
        command = f"ngs-pipeline --config {successful_config}"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        assert "[*] Starting alignment..." in result.stdout

@pytest.mark.skip(reason="incomplete test")
class TestBasicAnalysis:
    def test_analysis_execution(self, successful_config):
        """This test verifies that the analysis step is executed when the corresponding flag is set to True in the config file

        Args:
            successful_config (_type_): pytest fixture for a successful configuration file
        """
        
        # Run ngs_driver to check that analysis is performed
        command = f"ngs-pipeline --config {successful_config}"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        assert "[*] Starting analysis..." in result.stdout
        
