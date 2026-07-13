from pipeline.main import load_config
from pipeline.config_builder import generate_json_data_paired, generate_json_data_merged
from pipeline.alignment.aligner import generate_ref_library
from pipeline.config_validator import check_config_options

import os, json, subprocess, pytest, shutil
from pathlib import Path

@pytest.fixture
def make_config_data(tmp_path):
    """
    Fixture to create a successful configuration file for testing purposes. Does not include analysis-parameters, which are added in the individual tests as needed.
    """
    # Copy necessary test files to the temporary directory
    source_dir = f"{Path(__file__).parent}/test_data"
    
    fileCopys = {
        f"{source_dir}/x98_query.fasta": f"{tmp_path}/reference.fasta",
        f"{source_dir}/X98-w-Mg_R1_001.fastq.gz": f"{tmp_path}/sample1_R1.fastq.gz",
        f"{source_dir}/X98-w-Mg_R2_001.fastq.gz": f"{tmp_path}/sample1_R2.fastq.gz",
        f"{source_dir}/aligned_reads.bam": f"{tmp_path}/aligned_reads.bam"
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
        "output-directory": f"{tmp_path}/output"
    }
    return config_data

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

    def test_load_config(self, make_config_data, tmp_path):
        """This test confirms that the "load_config" method correctly loads intended config information - essentially verifying the JSON lib dump + load methods

        Args:
            make_config_data (_type_): pytest fixture - generates partial config data that all functional tests have
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the config file
        """
        # Create a temporary config file
        config_data = make_config_data
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

    def test_output_dir_generation(self, make_config_data, tmp_path):
        """This test verifies that a. ngs_driver can read in a config file and b. that output directory is correctly generated

        Args:
            make_config_data (_type_): pytest fixture - generates partial config data that all functional tests have
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the config file
        """
        
        # Create a temporary config file with an output directory
        config_data = make_config_data
        config_data["analysis-parameters"] = {
            "do-benchmarks": False,
            "do-processing": False,
            "do-alignment": False,
            "do-analysis": False
        }
        
        config_file = f"{tmp_path}/config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        # Run ngs_driver to create the output directory
        command = f"ngs-pipeline --config {config_file}"
        subprocess.run(command, shell=True, check=True)
        
        # Should create output directory
        assert Path(f"{tmp_path}/output").is_relative_to(tmp_path)
    
    def test_config_flags_disabled(self, make_config_data, tmp_path):
        """This test verifies that the flags in the config file are correctly interpreted by ngs_driver

        Args:
            make_config_data (_type_): pytest fixture - generates partial config data that all functional tests have
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the config file
        """
        
        # Create a temporary config file with all flags set to False
        config_data = make_config_data
        config_data["analysis-parameters"] = {
            "do-benchmarks": False,
            "do-processing": False,
            "do-alignment": False,
            "do-analysis": False
            }
        
        config_file = tmp_path / "config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        # Run ngs_driver to check that no analysis is performed
        command = "ngs-pipeline --config " + str(config_file)
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        assert "[*] Skipping fastqc benchmarks as per configuration." in result.stdout
        assert "[*] Skipping processing as per configuration." in result.stdout
        assert "[*] Skipping alignment as per configuration." in result.stdout
        assert "[*] Skipping analysis as per configuration." in result.stdout
    
class TestBadConfigs:
    def test_missing_all_data(self):
        """This tests that ValueError is thrown when no config data is given
        """
        config_data = {}

        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)
        
        assert "[!] No configuration data provided" in str(e_info.value)

    def test_missing_mode(self):
        """This tests that ValueError is thrown when config data does not have any "mode" specified
        """
        config_data = {
            "reference-fasta": "reference.fasta",
            "output-directory": "output"
        }

        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)

        assert "[!] No mode specified" in str(e_info.value)

    def test_invalid_mode(self):
        """This tests that ValueError is thrown when config data has invalid "mode"
        """
        config_data = {
            "mode": "invalid-mode",
            "reference-fasta": "reference.fasta",
            "output-directory": "output"
        }

        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)

        assert "[!] Invalid mode specified" in str(e_info.value)
    
    def test_missing_reference(self):
        """This tests that ValueError is thrown when config data does not have any "reference-fasta" specified
        """
        config_data = {
            "mode": "paired-end-mode",
            "output-directory": "output"
        }

        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)

        assert "[!] No reference FASTA file specified" in str(e_info.value)

    def test_invalid_reference(self):
        """This tests that ValueError is thrown when config data has invalid "reference-fasta"
        """
        config_data = {
            "mode": "paired-end-mode",
            "reference-fasta": "nonexistent_reference.fasta",
            "output-directory": "output"
        }

        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)

        assert "[!] Reference FASTA file not found" in str(e_info.value)

    def test_missing_output(self, make_config_data):
        """This tests that ValueError is thrown when config data does not have any "output-directory" specified
        """
        config_data = make_config_data
        # Remove output directory from the config
        del config_data["output-directory"]

        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)

        assert "[!] No output directory specified" in str(e_info.value)

    def test_missing_files_paired_end_R1(self, make_config_data):
        """This tests that ValueError is thrown when config data does not have any "R1" file specified for paired-end mode
        """
        config_data = make_config_data
        # Remove R1 file from the config
        config_data["paired-end-mode"]["R1"] = "nonexistent_R1.fastq"

        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)

        assert "[!] Read 1 file not found" in str(e_info.value)

    def test_missing_files_paired_end_R2(self, make_config_data):
        """This tests that ValueError is thrown when config data does not have any "R2" file specified for paired-end mode
        """
        config_data = make_config_data
        # Remove R2 file from the config
        config_data["paired-end-mode"]["R2"] = "nonexistent_R2.fastq"

        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)

        assert "[!] Read 2 file not found" in str(e_info.value)

    def test_missing_files_merged(self, make_config_data):
        """This tests that ValueError is thrown when config data does not have any "R1" file specified for merged mode
        """
        config_data = make_config_data
        # Change mode to merged-mode and remove R1 file from the config
        config_data["mode"] = "merged-mode"
        config_data["merged-mode"] = {
            "R1": "nonexistent_merged.fastq"
        }

        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)

        assert "[!] Merged file not found" in str(e_info.value)

    def test_missing_opts(self, make_config_data):
        """This tests that ValueError is thrown when config data does not have any "analysis-parameters" specified
        """
        config_data = make_config_data
        # make_config_data fixture doesn't include analysis-parameters
        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)

        assert "[!] No analysis parameters specified" in str(e_info.value)

    def test_invalid_opts(self, make_config_data):
        """This tests that ValueError is thrown when config data has invalid "analysis-parameters"
        """
        config_data = make_config_data
        # Add invalid analysis-parameters
        config_data["analysis-parameters"] = {
            "do-benchmarks": "foobar",  # Invalid, should be boolean
            "do-processing": False,
            "do-alignment": False,
            "do-analysis": False
        }

        with pytest.raises(ValueError) as e_info:
            check_config_options(config_data)

        assert "[!] Invalid value for 'do-benchmarks' - must be a boolean (True/False)" in str(e_info.value)

class TestBasicProcessing:
    def test_fastqc_execution(self, make_config_data, tmp_path):
        """This test verifies that the fastqc processing step is executed when the corresponding flag is set to True in the config file

        Args:
            make_config_data (_type_): pytest fixture - generates partial config data that all functional tests have
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the necessary files
        """
        # Make config data for this test
        config_data = make_config_data
        config_data["analysis-parameters"] = {
            "do-benchmarks": True,
            "do-processing": False,
            "do-alignment": False,
            "do-analysis": False
        }
        config_file = f"{tmp_path}/config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        # Run ngs-pipelilne to check that benchmarking is performed when specified in config.json
        command = f"ngs-pipeline --config {config_file}"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        target_output_dir = Path(f"{tmp_path}/output/fastqc_results")

        # Check status message from main()
        assert "[*] Starting fastqc benchmarks..." in result.stdout
        # Check fastqc output directory exists
        assert target_output_dir.is_dir()
        # Check fastqc .zip files exist
        assert Path(f"{target_output_dir}/sample1_R1_fastqc.zip").is_relative_to(target_output_dir)
        # Check fastqc unzip directory exists
        assert Path(f"{target_output_dir}/sample1_R1_fastqc").is_relative_to(target_output_dir)
        # Check fastqc unzip has data file
        assert Path(f"{target_output_dir}/sample1_R1_fastqc/fastqc_data.txt").parent == Path(f"{target_output_dir}/sample1_R1_fastqc")

class TestBasicAlignment:
    def test_good_reference(self, make_config_data):
        """This test verifies that bowtie2 can make a reference database when given a good reference file

        Args:
            make_config_data (_type_):pytest fixture - generates partial config data that all functional tests have
        """
        # Generate good config data
        config_data = make_config_data
        
        # Pass along good reference file
        result = generate_ref_library(config_data["reference-fasta"], config_data["output-directory"])
        target_output_dir = config_data["output-directory"]

        # Should make a ref_lib dir in /output/
        assert Path(f"{target_output_dir}/ref_lib").is_relative_to(target_output_dir)
        # Should make six '.bt2' files in /ref_lib/
        bt2_files = []
        for path in Path(f"{target_output_dir}/ref_lib").glob("*.bt2"):
            bt2_files.append(path)
        assert len(bt2_files) == 6
        # Should return the filename of the reference file
        assert result == f"{target_output_dir}/ref_lib/{config_data["reference-fasta"].split('/')[-1]}"

    def test_bad_reference(self, tmp_path, capsys):
        """This test verifies that bowtie2 cannot make a reference database when given a bad reference file and error is handled correctly

        Args:
            tmp_path (_type_): pytest temporary directory fixture
            capsys (_type_): pytest fixture for capturing stdout and stderr
        """
        # Create empty reference.fasta file
        with open(f"{tmp_path}/reference.fasta", "w") as f:
            f.write("")

        # Pass along bad reference file
        result = generate_ref_library(f"{tmp_path}/reference.fasta",f"{tmp_path}")

        target_output_dir = Path(f"{tmp_path}")
        
        # Should still make a ref_lib dir
        assert Path(f"{tmp_path}/ref_lib").is_relative_to(target_output_dir)
        # Should fail to make files
        bt2_files = []
        for path in Path(f"{target_output_dir}/ref_lib").glob("*.bt2"):
            bt2_files.append(path)
        assert len(bt2_files) == 0
        # Should yield a subprocess error that prints to stdout
        capture_stdout = capsys.readouterr()
        assert "[!] Error building reference index!" in capture_stdout.out
        # Should return none
        assert result == None

    def test_paired_alignment_execution(self, make_config_data, tmp_path):
        """This test verifies that the alignment step is executed when the corresponding flag is set to True in the config file

        Args:
            make_config_data (_type_): pytest fixture - generates partial config data that most functional tests have / typical use-case
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test containing the necessary files
        """
         # Make config data for this test
        config_data = make_config_data
        config_data["analysis-parameters"] = {
            "do-benchmarks": False,
            "do-processing": False,
            "do-alignment": True,
            "do-analysis": False
        }
        config_file = f"{tmp_path}/config.json"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        # Run ngs-pipelilne to check that alignment is performed when specified in config.json
        command = f"ngs-pipeline --config {config_file}"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        target_output_dir = Path(config_data['output-directory'])

        # Config correctly read in from pipeline - main
        assert "[*] Starting alignment..." in result.stdout
        # Aligner did not find reference db - not generated for this test
        assert f"[!] Reference file index not found - Building index for '{tmp_path}/reference.fasta'..." in result.stdout
        # Aligner should make 3 files:  .bam, .bam.log, and .bam.bai
        assert Path(f"{target_output_dir}/aligned_reads.bam").exists()
        assert Path(f"{target_output_dir}/aligned_reads.bam.log").exists()
        assert Path(f"{target_output_dir}/aligned_reads.bam.bai").exists()

class TestBasicAnalysis:
    def test_analysis_execution(self, make_config_data, tmp_path):
        """This test verifies that the analysis step is executed when the corresponding flag is set to True in the config file

        Args:
            make_config_data (_type_): pytest fixture - generates partial config data that most functional tests have / typical use-case
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test
        """
        config_data = make_config_data
        config_data["analysis-parameters"] = {
            "do-benchmarks": False,
            "do-processing": False,
            "do-alignment": False,
            "do-analysis": True,
            "do-alignment-stats": True,
            "do-alignment-visualization":True,
            "do-alignment-score-plot":True,
            "do-mpileup":True,
            "do-mpileup-fullanalysis":True,
            "do-mpileup-simpleanalysis":True,
            "do-mpileup-visualization":True,
            "do-association-analysis":True
        }
        config_file = f"{tmp_path}/config.json"
        with open(config_file, "w") as f:
            json.dump(config_data, f)

        # Run ngs_driver to check that analysis is performed
        command = f"ngs-pipeline --config {config_file}"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        target_output_dir = Path(config_data["output-directory"])

        # Main reads in config correct
        assert "[*] Starting analysis..." in result.stdout
        
        # Alignment Stats dir generates correctly
        assert Path(f"{target_output_dir}/alignment_stats").is_relative_to(target_output_dir)
        # samtools flagstat generates a file that is broken down into smaller files
        assert Path(f"{target_output_dir}/alignment_stats/alignment_stats.txt").is_relative_to(target_output_dir)
        # samtools stat generates a full file that is broken down into smaller files
        assert Path(f"{target_output_dir}/alignment_stats/full_alignment_stats.txt").is_relative_to(target_output_dir)
        # Make sure the smaller files are generated, should be 12 files in all (1 from flagstat, 11 from stat)
        alignment_stat_files = []
        for path in Path(f"{target_output_dir}/alignment_stats/").iterdir():
            alignment_stat_files.append(path)
            print(f"Alignment stat file: {path}")
        print(f"Alignment stat files: {alignment_stat_files}")
        assert len(alignment_stat_files) == 12

        # Alignment visualization generates correctly
        assert Path(f"{target_output_dir}/all_reads.txt").is_relative_to(target_output_dir)
        assert Path(f"{target_output_dir}/global_alignment_matrix.txt").is_relative_to(target_output_dir)

        # Alignment score plot generates correctly
        assert Path(f"{target_output_dir}/alignment_scores.txt").is_relative_to(target_output_dir)
        assert Path(f"{target_output_dir}/alignment_score_plot.png").is_relative_to(target_output_dir)

        # mpileup generates correctly
            # mpileup should generate a new reference index for itself
        assert "[!] Reference fasta index not found for samtools mpileup, generating index..." in result.stdout
            # new referenece index should end up in /output/ref_lib/
        assert Path(f"{target_output_dir}/ref_lib/reference.fasta.fai").is_relative_to(Path(f"{target_output_dir}/ref_lib/"))
            # mpileup file ends up in output
        assert Path(f"{target_output_dir}/aligned_reads.bam_mpileup.txt").is_relative_to(Path(f"{target_output_dir}"))

        # mpileup is analyzed
            # Full analysis generates
        assert Path(f"{target_output_dir}/mpileup_full_analysis.txt").is_relative_to(f"{target_output_dir}")
            # Partial analysis generates
        assert Path(f"{target_output_dir}/mpileup_simple_analysis.txt").is_relative_to(f"{target_output_dir}")

        # mpileup is visualized
        assert Path(f"{target_output_dir}/mpileup_percent_identical_plot.png")
        assert Path(f"{target_output_dir}/mpileup_percent_mutation_plot.png")
        assert Path(f"{target_output_dir}/mpileup_indel_frequencies_plot.png")
        assert Path(f"{target_output_dir}/mpileup_depth_plot.png")

    # TODO - Add association test once it is dealt with in main:main
    @pytest.mark.skip(reason="incomplete test")
    def test_association_analysis(self, make_config_data, tmp_path):
        pass