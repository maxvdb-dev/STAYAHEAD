from Bio import PDB
from Bio.PDB import Superimposer
from Bio.PDB import DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from pathlib import Path
import os
import csv
import argparse
# from scipy.spatial import distance
import numpy as np
# import matplotlib.pyplot as plt
import math
from tmtools.io import get_structure, get_residue_data
from tmtools import tm_align
import json
import freesasa
from Bio.PDB import PDBParser, Polypeptide
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd
from Bio.SeqUtils.ProtParamData import hw

# python /home/max/stayahead/snellius2/scripts/analysis/pred_analysis.py -i /home/max/stayahead/analysis/datasets/tmp/omicron_variants/af/pdb/ba1_100 -o /home/max/stayahead/analysis/datasets/tmp/omicron_variants/af/outputs -r /home/max/stayahead/analysis/datasets/references/alphafold/pdb/SARS-CoV-ACE2-AF_1.pdb -t alphafold -m 1step -b True
# python /home/max/stayahead/snellius2/scripts/analysis/pred_analysis.py -i /home/max/stayahead/out_tmp/esmfold/1step -o /home/max/stayahead/analysis/datasets/predictions/pred_ref -r /home/max/stayahead/snellius2/outputs/references/SARS-CoV-ACE2-ESM.pdb -t esmfold -b True -m 1step
# python /home/max/stayahead/snellius2/scripts/analysis/pred_analysis.py -i /home/max/stayahead/snellius2/outputs/references/alphafold/pdb -o /home/max/stayahead/snellius2/outputs/references/alphafold -r /home/max/stayahead/snellius2/outputs/references/alphafold/SARS-CoV-ACE2-AF_1.pdb -t alphafold -m 1step


class ComputeMetrics():
    def __init__(self, input_path, ref_path):
        self.input = input_path
        self.ref = ref_path
        # q = query, r = reference, n = name
        # self.q, self.r, self.n = self.load_structures()

    
    def load_structures(self):
        p = PDB.PDBParser(QUIET=True)
        struct_path = Path(self.input)
        ref_path = Path(self.ref)
        parent_dir = struct_path.parent
        key_name = parent_dir.name
        fixed_struct = p.get_structure('fixed', ref_path)
        moving_struct = p.get_structure('moving', struct_path)

        return fixed_struct, moving_struct, key_name

    def calculate_rmsd(self, fixed_atoms, moving_atoms, structure):
        super_imposer = Superimposer()
        super_imposer.set_atoms(fixed_atoms, moving_atoms)
        super_imposer.apply(structure.get_atoms())
        return super_imposer.rms
    
    def rmsd(self):
        # fixed_model = self.r[0]
        # p = PDB.PDBParser(QUIET=True)
        # struct_path = Path(self.input)
        # ref_path = Path(self.ref)
        # parent_dir = struct_path.parent
        # key_name = str(parent_dir.name)
        # fixed_struct = p.get_structure('fixed', ref_path)
        # moving_struct = p.get_structure('moving', struct_path)
        fixed_struct, moving_struct, key_name = self.load_structures()
        fixed_model = fixed_struct[0]
        for chain in fixed_model:
            # if chain.id == 'E':
            all_atoms_fixed = [atom for residue in chain.get_residues() for atom in residue]
            all_atoms_moving = [atom for residue in moving_struct.get_residues() for atom in residue]
            if len(all_atoms_fixed) > len(all_atoms_moving):
                dif = len(all_atoms_fixed) - len(all_atoms_moving)
                all_atoms_fixed = all_atoms_fixed[:len(all_atoms_fixed) - dif]
            else:
                dif = len(all_atoms_moving) - len(all_atoms_fixed)
                all_atoms_moving = all_atoms_moving[:len(all_atoms_moving) - dif]
                # If dif is positive, you might need to address how to handle this case,
                # such as choosing which atoms to remove from fixed_crys_atoms or reconsidering alignment.
            # print(len(all_atoms_fixed), len(all_atoms_moving))
            assert len(all_atoms_fixed) == len(all_atoms_moving)
            rmsd = self.calculate_rmsd(all_atoms_fixed, all_atoms_moving, chain)
                # print(rmsd)
        return rmsd

    def tm_score(self):
        ref_struct, query_struct, key_name = self.load_structures()
        chain = next(ref_struct.get_chains())
        ref_coords, ref_seq = get_residue_data(chain)
                    # If dif is positive, you might need to address how to handle this case,
                    # such as choosing which atoms to remove from fixed_crys_atoms or reconsidering alignment.
        chain = next(query_struct.get_chains())
        coords, seq = get_residue_data(chain)
        res = tm_align(ref_coords, coords, ref_seq, seq)
        tm_score = res.tm_norm_chain1
        # print(tm_score)

        return tm_score
    

    
    def sasa(self):
        structure = freesasa.Structure(str(self.input))
        result = freesasa.calc(structure)
        total_sasa = result.totalArea()
        # res_dict = {}
        # for i, residue in enumerate(structure.residueAreas()):
        #     residue_sasa = sum(residue.total)
        #     res_dict[i] = residue_sasa
        return total_sasa
    
    def b_factor(self, al_type):
        fixed_struct, moving_struct, key_name = self.load_structures()
        len_prot = [atom for residue in moving_struct.get_residues() for atom in residue]
        bfactor = 0
        for model in moving_struct:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        bfactor += atom.get_bfactor()
        b_factor = bfactor / len(len_prot)
        if al_type == 'esmfold':
            b_factor = b_factor * 100
        else:
            return b_factor
        return b_factor

    def extract_sequence_from_pdb(self):
        fixed_struct, moving_struct, key_name = self.load_structures()
        model = moving_struct[0]
        ppb = Polypeptide.PPBuilder()
        for pp in ppb.build_peptides(model):
            return str(pp.get_sequence())
        
    def calculate_average_hydrophobicity(self, sequence):
        protein_analysis = ProteinAnalysis(sequence)
        hydrophobicity = protein_analysis.protein_scale(hw, window=7, edge=1.0)
        if len(hydrophobicity) == 0:
            return 'N/A'
        average_hydrophobicity = sum(hydrophobicity) / len(hydrophobicity)
        return average_hydrophobicity
    
    def calculate_total_hydrophobicity(self, sequence):
        total_hydrophobicity = 0
        for aa in sequence:
            if aa in kd:
                total_hydrophobicity += kd[aa]
        return total_hydrophobicity / len(sequence)

def process_control_directory(input_dir, control_dir, output_dir, al_type):
    results = []
    for subdir, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".pdb"):
                pdb_path = os.path.join(input_dir, file)
                control_path = os.path.join(control_dir, file)
                # Check if control_path exists
                if not os.path.exists(control_path):
                    continue
                print(pdb_path, control_path)
                key_name = file.split('.')[0]
                metrics = ComputeMetrics(pdb_path, control_path)
                rmsd = metrics.rmsd()
                sasa = metrics.sasa()
                tm_score = metrics.tm_score()
                plddt = metrics.b_factor(al_type)
                sequence = metrics.extract_sequence_from_pdb()
                avg_hydrophobicity = metrics.calculate_average_hydrophobicity(sequence)
                control_metrics = ComputeMetrics(control_path, control_path)
                c_sasa = control_metrics.sasa()
                c_b_factor = control_metrics.b_factor(al_type)
                c_sequence = control_metrics.extract_sequence_from_pdb()
                c_avg_hydrophobicity = control_metrics.calculate_average_hydrophobicity(c_sequence)
                results.append([key_name, rmsd, tm_score, sasa, avg_hydrophobicity, plddt, c_sasa, c_avg_hydrophobicity, c_b_factor])
    
    output_path = os.path.join(output_dir, 'controls_metrics.csv')
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['seq_id', 'rmsd', 'tm_score', 'sasa', 'avg_hydro', 'plddt', 'c_sasa', 'c_avg_hydro', 'c_b_factor'])
        writer.writerows(results)


def process_directory(input_dir, ref_file, output_dir, al_type, batch, mut):
    metrics = ComputeMetrics(input_dir, ref_file)
    results = []

    if batch == True:
        for subdir, dirs, files in os.walk(input_dir):
            for file in files:
                if file.endswith("ranked_0.pdb") and al_type == 'alphafold':
                    pdb_path = os.path.join(subdir, file)
                    print(pdb_path)
                    key_name = Path(subdir).stem
                    name = key_name.split('_')[1]
                    split_name = list(name)
                    or_aa = split_name[0]
                    pos = int(''.join(split_name[1:-1]))
                    mut_aa = split_name[-1]
                    metrics = ComputeMetrics(pdb_path, ref_file)
                    rmsd = metrics.rmsd()
                    sasa = metrics.sasa()
                    tm_score = metrics.tm_score()
                    plddt = metrics.b_factor(al_type)
                    sequence = metrics.extract_sequence_from_pdb()
                    avg_hydrophobicity = metrics.calculate_average_hydrophobicity(sequence)
                    results.append([name, or_aa, pos, mut_aa, rmsd, tm_score, sasa, avg_hydrophobicity, plddt])
                elif file.endswith(".pdb") and al_type == 'esmfold':
                    pdb_path = os.path.join(subdir, file)
                    key_name = file.split('.')[0]
                    name = key_name.split('-')[1]
                    print(name)
                    split_name = list(name)
                    or_aa = split_name[0]
                    pos = int(''.join(split_name[1:-1]))
                    mut_aa = split_name[-1]
                    metrics = ComputeMetrics(pdb_path, ref_file)
                    rmsd = metrics.rmsd()
                    sasa = metrics.sasa()
                    tm_score = metrics.tm_score()
                    plddt = metrics.b_factor(al_type)
                    sequence = metrics.extract_sequence_from_pdb()
                    avg_hydrophobicity = metrics.calculate_average_hydrophobicity(sequence)
                    total_hydrophobicity = metrics.calculate_total_hydrophobicity(sequence)
                    results.append([name, or_aa, pos, mut_aa, rmsd, tm_score, sasa, avg_hydrophobicity, total_hydrophobicity, plddt])

    else:
        for subdir, dirs, files in os.walk(input_dir):
            for file in files:
                    if file.endswith(".pdb") and al_type == 'alphafold':
                        pdb_path = os.path.join(input_dir, file)
                        print(pdb_path)
                        key_name = file.split('.')[0]
                        if mut == '2step':
                            name_1 = key_name.split('_')[1]
                            name_2 = key_name.split('_')[2]
                            name = name_1 + '-' + name_2
                            print(name)
                        else:
                            name = key_name.split('-')[1]
                            split_name = list(name)
                            or_aa = split_name[0]
                            pos = int(''.join(split_name[1:-1]))
                            mut_aa = split_name[-1]
                        metrics = ComputeMetrics(pdb_path, ref_file)
                        rmsd = metrics.rmsd()
                        sasa = metrics.sasa()
                        tm_score = metrics.tm_score()
                        plddt = metrics.b_factor(al_type)
                        sequence = metrics.extract_sequence_from_pdb()
                        print(sequence)
                        avg_hydrophobicity = metrics.calculate_average_hydrophobicity(sequence)
                        total_hydrophobicity = metrics.calculate_total_hydrophobicity(sequence)
                        results.append([name, or_aa, pos, mut_aa, rmsd, tm_score, sasa, avg_hydrophobicity, total_hydrophobicity, plddt])
                        # results.append([name, rmsd, tm_score, sasa, avg_hydrophobicity, total_hydrophobicity, plddt])
                    elif file.endswith(".pdb") and al_type == 'esmfold':
                        pdb_path = os.path.join(input_dir, file)
                        key_name = file.split('.')[0]
                        if mut == '2step':
                            name_1 = key_name.split('_')[1]
                            name_2 = key_name.split('_')[2]
                            name = name_1 + '-' + name_2
                            split_name_1 = list(name_1)
                            split_name_2 = list(name_2)
                            or_aa_1 = split_name_1[0]
                            or_aa_2 = split_name_2[0]
                            pos_1 = int(''.join(split_name_1[1:-1]))
                            pos_2 = int(''.join(split_name_2[1:-1]))
                            mut_aa_1 = split_name_1[-1]
                            mut_aa_2 = split_name_2[-1]
                        else:
                            name = key_name.split('-')[1]
                            print(name)
                            split_name = list(name)
                            or_aa = split_name[0]
                            pos = int(''.join(split_name[1:-1]))
                            mut_aa = split_name[-1]
                        metrics = ComputeMetrics(pdb_path, ref_file)
                        rmsd = metrics.rmsd()
                        sasa = metrics.sasa()
                        tm_score = metrics.tm_score()
                        b_factor = metrics.b_factor()
                        sequence = metrics.extract_sequence_from_pdb()
                        avg_hydrophobicity = metrics.calculate_average_hydrophobicity(sequence)
                        total_hydrophobicity = metrics.calculate_total_hydrophobicity(sequence)
                        if mut == '2step':
                            results.append([name, or_aa_1, pos_1, mut_aa_1, or_aa_2, pos_2, mut_aa_2, rmsd, tm_score, sasa, b_factor])
                        else:
                            results.append([name, or_aa, pos, mut_aa, rmsd, tm_score, sasa, avg_hydrophobicity, total_hydrophobicity, b_factor])

    # Save results to CSV file
    if al_type == 'alphafold':
        output_path = os.path.join(output_dir, 'AF_BA1_100.csv')
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            if mut == '2step':
                writer.writerow(['seq_id', 'wildtype_1', 'site_1', 'mutation_1', 'wildtype_2', 'site_2', 'mutation_2', 'rmsd', 'tm_score', 'sasa', 'plddt'])
            else:
                writer.writerow(['seq_id', 'wildtype', 'site', 'mutation', 'rmsd', 'tm_score', 'sasa', 'avg_hydro', 'plddt'])
            writer.writerows(results)
    elif al_type == 'esmfold':
        output_path = os.path.join(output_dir, 'ESM_1step_pred_14-06.csv')
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            if mut == '2step':
                writer.writerow(['seq_id', 'wildtype_1', 'site_1', 'mutation_1', 'wildtype_2', 'site_2', 'mutation_2', 'rmsd', 'tm_score', 'sasa', 'b_factor'])
            else:
                writer.writerow(['seq_id', 'wildtype', 'site', 'mutation', 'rmsd', 'tm_score', 'sasa', 'avg_hydro', 'total_hydro', 'b_factor'])
            writer.writerows(results)

def create_parser():
    parser = argparse.ArgumentParser(description="Compute structural metrics.")
    parser.add_argument(
        "-i", "--input_dir", help="Path to input directory", type=Path, required=True,
    )
    parser.add_argument(
        "-o", "--out_dir", help="Path to save directory for sequences", type=Path, required=True,
    )
    parser.add_argument(
        "-r", "--ref", help="Path to reference structure or directory", type=Path, required=True,
    )
    parser.add_argument(
        "-t", "--type", help="Type of algorithm used", type=str, required=True,
    )
    parser.add_argument(
        "-b", "--batch", help="Is batch dir", type=bool, default=False,
    )
    parser.add_argument(
        "-m", "--mutation", help="Mutation type", type=str, required=True,
    )
    parser.add_argument(
        "-c", "--control", help="Input dataset consists of control sequences", type=bool, default=False,
    )
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    if args.control == False:
        process_directory(args.input_dir, args.ref, args.out_dir, args.type, args.batch, args.mutation)
    else:
        process_control_directory(args.input_dir, args.ref, args.out_dir, args.type)


if __name__ == "__main__":
    main()


