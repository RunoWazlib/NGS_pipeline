import subprocess
from pathlib import Path

target_dir = Path(__file__).resolve().parent.parent.parent

def run_setup():
    subprocess.run(f"bash {target_dir}/setup.sh", shell=True, check=True)

def run_dep_update():
    subprocess.run(f"bash {target_dir}/update_env_doc.sh", shell=True, check=True)