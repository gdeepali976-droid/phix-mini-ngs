"""Run the whole pipeline in order. Takes under a minute."""
import subprocess
import sys

for script in ["00_fetch_reference.py", "01_simulate_reads.py",
               "02_read_qc.py", "03_assemble.py", "04_evaluate.py",
               "05_parameter_sweep.py", "test_known_values.py"]:
    print(f"\n=== {script} ===")
    subprocess.run([sys.executable, script], check=True)
