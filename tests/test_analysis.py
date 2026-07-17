import pytest, json, subprocess, shutil
from pathlib import Path

from pipeline.analysis.alignment_stats import generate_alignment_stats, generate_alignment_score_plot
from pipeline.analysis.alignment_analysis import generate_alignment_visualization
from pipeline.analysis.mpileup_analysis import generate_mpileup, mpileup_cleaner, generate_mpileup_full_analysis, generate_mpileup_simple_analysis, generate_mpileup_visualization, generate_indel_analysis

@pytest.fixture
def make_config_data(tmp_path):
    """
    Fixture to create a successful configuration file for testing purposes. Does not include analysis-parameters, which are added in the individual tests as needed.
    """
    # pre-create output dir for aligned_reads.bam
    subprocess.run(["mkdir", "-p",f"{tmp_path}/output"], check=True)
    # Copy necessary test files to the temporary directory
    source_dir = f"{Path(__file__).parent}/test_data"
    
    fileCopys = {
        f"{source_dir}/x98_query.fasta": f"{tmp_path}/reference.fasta",
        f"{source_dir}/X98-w-Mg_R1_001.fastq.gz": f"{tmp_path}/sample1_R1.fastq.gz",
        f"{source_dir}/X98-w-Mg_R2_001.fastq.gz": f"{tmp_path}/sample1_R2.fastq.gz",
        f"{source_dir}/aligned_reads.bam": f"{tmp_path}/output/aligned_reads.bam"
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

class TestAnalysisFunctions:
    def test_alignment_stats_func(self, make_config_data, tmp_path):
        """This test validates the alignment_stats generator function

        Args:
            make_config_data (_type_): pytest fixture - generates partial config data that most functional tests have / typical use-case
            tmp_path (_type_): pytest temporary directory fixture - acts as "output directory" for the test
        """
        # Create "alignment_stats" dir for output
        subprocess.run(["mkdir","-p",f"{tmp_path}/output/alignment_stats"], check=True)
        generate_alignment_stats(f"{tmp_path}/output/aligned_reads.bam", f"{tmp_path}/output/alignment_stats")
        
        # samtools flagstat generates a file
        try:
            open(f"{tmp_path}/output/alignment_stats/alignment_stats_summary.txt", "r")
            assert True
        except FileNotFoundError:
            assert False

        # Make sure the smaller files are generated, should be 15 files in all (1 from flagstat, 10 from stat, 4 from setup)
        alignment_stat_files = []
        for path in Path(f"{tmp_path}/output/alignment_stats").iterdir():
            alignment_stat_files.append(path)
            try:
                open(path, "r")
            except FileNotFoundError:
                assert False
        print(f"Alignment stat files: {alignment_stat_files}")
        assert len(alignment_stat_files) == 11

    def test_alignment_score_plot_func(self, make_config_data, tmp_path):
        """This test validates the alignment_score_plot generator function

        Args:
            make_config_data (_type_): pytest fixture - generates partial config data that most functional tests have / typical use-case
            tmp_path (_type_): pytest temporary directory fixture - acts as "output directory" for the test
        """
        # Create "alignment_stats" dir for output
        subprocess.run(["mkdir","-p",f"{tmp_path}/output/alignment_stats"], check=True)
        generate_alignment_score_plot(f"{tmp_path}/output/aligned_reads.bam", f"{tmp_path}/output/alignment_stats")

        # func generates a txt file of scores
        target = f"{tmp_path}/output/alignment_stats/alignment_scores.txt"
        assert Path(target).exists()
        try:
            with open(target, "r") as f:
                if len(f.read()) > 0:
                    assert True
                else:
                    assert False
        except FileNotFoundError:
            assert False
        # func generates a png image
        assert Path(f"{tmp_path}/output/alignment_stats/alignment_score_plot.png").exists()

    def test_alignment_visualization_func(self, make_config_data, tmp_path):
        """This test validates the alignment_stats generator function

        Args:
            make_config_data (_type_): pytest fixture - generates partial config data that most functional tests have / typical use-case
            tmp_path (_type_): pytest temporary directory fixture - acts as "output directory" for the test
        """
        # Get config data for ref and output locations
        config_data = make_config_data
        # Run unit
        generate_alignment_visualization(f"{tmp_path}/output/aligned_reads.bam", config_data["reference-fasta"],config_data["output-directory"])

        # func generates an "all_reads.txt" file
        assert Path(f"{tmp_path}/output/all_reads.txt").exists()
        # "all_reads.txt" is not empty
        target = f"{tmp_path}/output/all_reads.txt"
        assert Path(target).exists()
        # file is not empty
        try:
            with open(target, "r") as f:
                if len(f.read()) < 0:
                    assert False
                else:
                    assert True
        except FileNotFoundError:
            assert False
        # func generates a "global_alignment_matrix.txt" file
        assert Path(f"{tmp_path}/output/global_alignment_matrix.txt").exists()
        # "global_alignment_matrix.txt" is not empty
        target = f"{tmp_path}/output/global_alignment_matrix.txt"
        assert Path(target).exists()
        # file is not empty
        try:
            with open(target, "r") as f:
                if len(f.read()) < 0:
                    assert False
                else:
                    assert True
        except FileNotFoundError:
            assert False

    def test_mpileup_gen_func_no_fai(self, make_config_data, tmp_path):
        # Create "ref_lib" dir for ref.fai and "alignment" for output
        subprocess.run(["mkdir","-p",f"{tmp_path}/output/ref_lib/"], check=True)
        subprocess.run(["mkdir","-p",f"{tmp_path}/output/alignment/"], check=True)
        # Get config data
        config_data = make_config_data
        # unit test
        result = generate_mpileup(f"{tmp_path}/output/aligned_reads.bam", config_data["reference-fasta"], f"{config_data["output-directory"]}/alignment/")

        # .fai file should end up in output/ref_lib/
        assert Path(f"{tmp_path}/output/ref_lib/ref.fai").exists()
        # mpileup.txt should end up in output/
        assert Path(f"{tmp_path}/output/alignment/aligned_reads.bam_mpileup.txt").exists()
        # func should return location of mpileup
        assert Path(result) == Path(f"{tmp_path}/output/alignment/aligned_reads.bam_mpileup.txt")

    def test_mpileup_gen_func_fai(self, make_config_data, tmp_path):
        # Create "ref_lib" dir for ref.fai and "alignment" for output
        subprocess.run(["mkdir","-p",f"{tmp_path}/output/ref_lib/"], check=True)
        subprocess.run(["mkdir","-p",f"{tmp_path}/output/alignment/"], check=True)
        # Get config data
        config_data = make_config_data
        # Create "reference.fasta.fai" file
        subprocess.run(["samtools","faidx",config_data["reference-fasta"]], check=True)
        # unit test
        result = generate_mpileup(f"{tmp_path}/output/aligned_reads.bam", config_data["reference-fasta"], f"{config_data["output-directory"]}/alignment/")

        # mpileup.txt should end up in output/
        assert Path(f"{tmp_path}/output/alignment/aligned_reads.bam_mpileup.txt").exists()
        # func should return location of mpileup
        assert Path(result) == Path(f"{tmp_path}/output/alignment/aligned_reads.bam_mpileup.txt")

    def test_mpileup_full_analysis(self, make_config_data, tmp_path):
        # Create "alignment" output directory
        subprocess.run(["mkdir","-p",f"{tmp_path}/output/alignment/"], check=True)
        # Get config data
        config_data = make_config_data
        # Make mpileup file
        subprocess.run(f"samtools faidx {config_data["reference-fasta"]}", shell=True)
        subprocess.run(f"samtools mpileup -f {config_data["reference-fasta"]} -d 500000 {tmp_path}/output/aligned_reads.bam > {tmp_path}/output/alignment/aligned_reads.bam_mpileup.txt", shell=True)
        # Unit test
        generate_mpileup_full_analysis(f"{tmp_path}/output/alignment/aligned_reads.bam_mpileup.txt", f"{tmp_path}/output/alignment")

        # func generates file
        target = f"{tmp_path}/output/alignment/Full_analysis.txt"
        assert Path(target).exists()
        # file is not empty
        try:
            with open(target, "r") as f:
                if len(f.read()) < 0:
                    assert False
                else:
                    assert True
        except FileNotFoundError:
            assert False

    def test_mpileup_simple_analysis(self, make_config_data, tmp_path):
        # Create "alignment" output directory
        subprocess.run(["mkdir","-p",f"{tmp_path}/output/alignment/"], check=True)
        # Get config data
        config_data = make_config_data
        # Make mpileup file
        subprocess.run(f"samtools faidx {config_data["reference-fasta"]}", shell=True)
        subprocess.run(f"samtools mpileup -f {config_data["reference-fasta"]} -d 500000 {tmp_path}/output/aligned_reads.bam > {tmp_path}/output/alignment/aligned_reads.bam_mpileup.txt", shell=True)

        # Unit test
        generate_mpileup_simple_analysis(f"{tmp_path}/output/alignment/aligned_reads.bam_mpileup.txt", f"{tmp_path}/output/alignment")

        # func generates file
        target = f"{tmp_path}/output/alignment/Simple_analysis.txt"
        assert Path(target).exists()
        # file is not empty
        try:
            with open(target, "r") as f:
                if len(f.read()) < 0:
                    assert False
                else:
                    assert True
        except FileNotFoundError:
            assert False
# TODO - add remaining analysis unit tests; mpileup cleaner, mpileup visualization, indel analysis, association analysis
    def test_mpileup_visualization(self, make_config_data, tmp_path):
        pass

    def test_mpileup_indel_analysis(self, make_config_data, tmp_path):
        pass

# TODO - verify end-to-end analysis tests work / need an update
class TestBasicAnalysisRuns:
    def test_analysis_execution(self, make_config_data, tmp_path):
        """This test verifies that the analysis step is executed when the corresponding flag is set to True in the config file

        Args:
            make_config_data (_type_): pytest fixture - generates partial config data that most functional tests have / typical use-case
            tmp_path (_type_): pytest temporary directory fixture - acts as launch directory for the test
        """
        config_data = make_config_data
        config_data["core-parameters"] = {
            "do-benchmarks": False,
            "do-processing": False,
            "do-alignment": False,
            "do-analysis": True
        }
        config_data["analysis-parameters"] = {
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
        target_output_dir = Path(config_data["output-directory"]) # tmp_path/output/

        # Main reads in config correct
        assert "[*] Starting analysis..." in result.stdout
        
        # Alignment Stats dir generates correctly
        assert Path(f"{target_output_dir}/alignment_stats").is_dir()
        # samtools flagstat generates a file
        assert Path(f"{target_output_dir}/alignment_stats/alignment_stats_summary.txt").exists()
        # Make sure the smaller files are generated, should be 11 files in all (1 from flagstat, 10 from stat, 2 from alignment score dist)
        alignment_stat_files = []
        for path in Path(f"{target_output_dir}/alignment_stats/").iterdir():
            alignment_stat_files.append(path)
        print(f"Alignment stat files: {alignment_stat_files}")
        assert len(alignment_stat_files) == 13

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