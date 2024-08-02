from pathlib import Path
import os
import csv
import argparse
import apbs
import numpy as np
import pdb2pqr
import subprocess
import re

# python /home/max/stayahead/snellius2/scripts/analysis/calc_epc.py -i /home/max/stayahead/analysis/datasets/tmp/omicron_variants/af/pdb/ba1_100 -o /home/max/stayahead/analysis/datasets/tmp/omicron_variants/af/outputs/epc -t alphafold -b True
# python /home/max/stayahead/snellius2/scripts/analysis/calc_epc.py -i /home/max/stayahead/analysis/datasets/tmp/omicron_variants/esm/ba1_100 -o /home/max/stayahead/analysis/datasets/tmp/omicron_variants/esm/ep/ba1 -t esmfold

def run_pdb2pqr(pdb_path, pqr_path, apbs_input_path, forcefield='PARSE'):
    """
    Convert a PDB file to a PQR file using PDB2PQR.
    :param pdb_path: Path to the input PDB file.
    :param pqr_path: Path to the output PQR file.
    :param forcefield: Forcefield to use, default is 'parse'.
    """
    directory, input_file = os.path.split(pqr_path)

    pdb2pqr_cmd = [
        'pdb2pqr',
        '--ff={}'.format(forcefield),
        "--apbs-input={}".format(apbs_input_path),
        pdb_path,
        input_file
    ]
    subprocess.run(pdb2pqr_cmd, cwd=directory, text=True)

        # "--nodebump",
        # "--noopt",

def run_apbs(apbs_input_path):
    directory, input_file = os.path.split(apbs_input_path)

    print("Running APBS with:")
    print("Input file path:", apbs_input_path)

    apbs_cmd = ['apbs', input_file]

    # Running the command and capturing output and errors
    process = subprocess.run(apbs_cmd, cwd=directory, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Print the complete output (optional)
    # print(process.stdout)

    # Print any errors
    if process.stderr:
        print("Errors encountered:", process.stderr)

    # Parse output to find energy values
    energy_values = []
    # extra_values = []
    for line in process.stdout.split('\n'):
        if "Local net energy (PE 0)" in line:
            # Extract the energy value using regex
            match = re.search(r"=\s*([\d.E+-]+)\s*kJ/mol", line)
            if match:
                energy_values.append(float(match.group(1)))
        elif "Global net ELEC energy" in line:
            match = re.search(r"=\s*([\d.E+-]+)\s*kJ/mol", line)
            if match:
                energy_values.append(float(match.group(1)))
        # else:
        #     extra_values.append(line)

    # Print the extracted values (optional)
    print('Local net energy', energy_values[0])
    print('Global net ELEC energy', energy_values[1])

    # You can now store these values in a file or return them
    return energy_values


def process_directory(input_dir,output_dir, type, batch):
    results = []
    if batch == True:
        for subdir, dirs, files in os.walk(input_dir):
            for file in files:
                if file.endswith("ranked_0.pdb") and type == 'alphafold':
                    pdb_path = os.path.join(subdir, file)
                    key_name = Path(subdir).stem
                    pqr_dir = os.path.join(output_dir, key_name)
                    os.makedirs(pqr_dir, exist_ok=True)
                    pqr_path = os.path.abspath(f'{pqr_dir}/{key_name}.pqr')
                    apbs_input_path = os.path.abspath(f'{pqr_dir}/{key_name}.in')
                    run_pdb2pqr(pdb_path, pqr_path, apbs_input_path)
                    energy_values = run_apbs(apbs_input_path)
                    results.append([key_name, energy_values[0], energy_values[1]])
                elif file.endswith(".pdb") and type == 'esmfold':
                    pdb_path = os.path.join(subdir, file)
                    key_name = file.split('.')[0]
                    pqr_dir = os.path.join(output_dir, key_name)
                    os.makedirs(pqr_dir, exist_ok=True)
                    pqr_path = os.path.abspath(f'{pqr_dir}/{key_name}.pqr')
                    apbs_input_path = os.path.abspath(f'{pqr_dir}/{key_name}.in')
                    run_pdb2pqr(pdb_path, pqr_path, apbs_input_path)
                    energy_values = run_apbs(apbs_input_path)
                    results.append([key_name, energy_values[0], energy_values[1]])
    
    else:
        for subdir, dirs, files in os.walk(input_dir):
            for file in files:
                    if file.endswith(".pdb") and type == 'alphafold':
                        pdb_path = os.path.join(input_dir, file)
                        print(pdb_path)
                        key_name = file.split('.')[0]
                        name = key_name.split('_')[0]
                        pqr_dir = os.path.join(output_dir, key_name)
                        os.makedirs(pqr_dir, exist_ok=True)
                        pqr_path = os.path.abspath(f'{pqr_dir}/{key_name}.pqr')
                        apbs_input_path = os.path.abspath(f'{pqr_dir}/{key_name}.in')
                        run_pdb2pqr(pdb_path, pqr_path, apbs_input_path)
                        energy_values = run_apbs(apbs_input_path)
                        results.append([name, energy_values[0], energy_values[1]])
                    elif file.endswith(".pdb") and type == 'esmfold':
                        pdb_path = os.path.join(input_dir, file)
                        key_name = file.split('.')[0]
                        name = key_name.split('_')[1]
                        pqr_dir = os.path.join(output_dir, key_name)
                        os.makedirs(pqr_dir, exist_ok=True)
                        pqr_path = os.path.abspath(f'{pqr_dir}/{key_name}.pqr')
                        apbs_input_path = os.path.abspath(f'{pqr_dir}/{key_name}.in')
                        run_pdb2pqr(pdb_path, pqr_path, apbs_input_path)
                        energy_values = run_apbs(apbs_input_path)
                        results.append([name, energy_values[0], energy_values[1]])


    # Save results to CSV file
    output_path = os.path.join(output_dir, 'AF_BA1_100_EP.csv')
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['sequence name', 'local_net_energy', 'global_net_energy'])
        writer.writerows(results)
    # txt_path = os.path.join(output_dir, 'AF_epc_parse.txt')
    # with open(txt_path, 'w') as f:
    #     for item in extra_values:
    #         f.write("%s\n" % item)

def create_parser():
    parser = argparse.ArgumentParser(description="Compute Electrostatic Potential Calculations.")
    parser.add_argument(
        "-i", "--input_dir", help="Path to input directory", type=Path, required=True,
    )
    parser.add_argument(
        "-o", "--out_dir", help="Path to save directory for sequences", type=Path, required=True,
    )
    parser.add_argument(
        "-t", "--type", help="Type of algorithm used", type=str, required=True,
    )
    parser.add_argument(
        "-b", "--batch", help="Is batch dir", type=bool, default=False,
    )
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    process_directory(args.input_dir, args.out_dir, args.type, args.batch)


if __name__ == "__main__":
    main()