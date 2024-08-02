import os
import apbs
import pdb2pqr
import subprocess

def run_apbs(apbs_input_path):
    apbs_input_path = os.path.abspath(apbs_input_path)

    print("Running APBS with:")
    print("Input file path:", apbs_input_path)

    apbs_cmd = [
        '/usr/bin/apbs',
        apbs_input_path
    ]

    # Running the command and capturing output and errors
    process = subprocess.run(apbs_cmd)

apbs_input_path = '/home/max/stayahead/snellius2/outputs/apbs/test/ACE2.in'
run_apbs(apbs_input_path)