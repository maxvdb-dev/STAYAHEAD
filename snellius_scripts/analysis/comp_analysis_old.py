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

# /home/max/stayahead/snellius2/scripts/analysis/comp_analysis.py -i 

class ComputeMetrics():
    def __init__(self, args):
        self.args = args
        self.input = args.input
        self.output = args.out_dir
        self.ref = args.ref
        # q = query, r = reference, n = name
        self.q, self.r, self.n = self.load_structures()

    
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
    
    def rmsd(self, dict):
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
            if chain.id == 'E':
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
                print(len(all_atoms_fixed), len(all_atoms_moving))
                assert len(all_atoms_fixed) == len(all_atoms_moving)
                rmsd = self.calculate_rmsd(all_atoms_fixed, all_atoms_moving, chain)
                print(rmsd)
                dict[key_name] = rmsd

        return dict

    def tm_score(self, dict):
        # fixed_model = self.r[0]
        # p = PDB.PDBParser(QUIET=True)
        # struct_path = Path(self.input)
        # ref_path = Path(self.ref)
        # parent_dir = struct_path.parent
        # key_name = str(parent_dir.name)
        # fixed_struct = p.get_structure('fixed', ref_path)
        # moving_struct = p.get_structure('moving', struct_path)
        ref_struct, query_struct, key_name = self.load_structures()
        chain = next(ref_struct.get_chains())
        ref_coords, ref_seq = get_residue_data(chain)
                    # If dif is positive, you might need to address how to handle this case,
                    # such as choosing which atoms to remove from fixed_crys_atoms or reconsidering alignment.
        chain = next(query_struct.get_chains())
        coords, seq = get_residue_data(chain)
        res = tm_align(ref_coords, coords, ref_seq, seq)
        tm_score = res.tm_norm_chain1
        print(tm_score)
        dict[key_name] = tm_score

        return dict
    
    def sasa(self, dict):
        structure = freesasa.Structure(str(self.input))
        result = freesasa.calc(structure)
        total_sasa = result.totalArea()
        dict[self.n] = total_sasa
        # res_dict = {}
        # for i, residue in enumerate(structure.residueAreas()):
        #     residue_sasa = sum(residue.total)
        #     res_dict[i] = residue_sasa
        return dict

def create_parser():
    parser = argparse.ArgumentParser(description="Compute structural metrics.")
    parser.add_argument(
        "-i", "--input", help="Path to input sequence", type=Path, required=True,
    )
    parser.add_argument(
        "-o", "--out_dir", help="Path to save directory for sequences", type=Path, required=True,
    )
    parser.add_argument(
        "-r", "--ref", help="Path to reference structure", type=Path, required=True,
    )
    return parser


def calculate_metrics(args):
    out_dict = {}
    metrics = ComputeMetrics(args)
    rmsd_dict = {}
    tm_score_dict = {}
    sasa_dict = {}
    sasa_dict = metrics.sasa(sasa_dict)
    rmsd_dict = metrics.rmsd(rmsd_dict)
    tm_score_dict = metrics.tm_score(tm_score_dict)
    out_dict['rmsd'] = rmsd_dict
    out_dict['tm_score'] = tm_score_dict
    out_dict['sasa'] = sasa_dict

    return out_dict, metrics.n

def main():
    parser = create_parser()
    args = parser.parse_args()
    out_dict, name = calculate_metrics(args)
    out_dir = os.path.join(args.out_dir, f'metrics_{name}.json')
    with open(out_dir, 'w') as file:
        json.dump(out_dict, file, indent=4)


if __name__ == "__main__":
    main()






























    # def rmsd(self):
    #     p = PDB.PDBParser(QUIET=True)
    #     struct_path = Path(self.input)
    #     ref_path = Path(self.ref)
    #     parent_dir = struct_path.parent
    #     key_name = parent_dir.name
    #     fixed_struct = p.get_structure('fixed', ref_path)
    #     moving_struct = p.get_structure('moving', struct_path)
    #     fixed_model = fixed_struct[0]
    #     for chain in fixed_model:
    #         all_atoms_fixed = [atom for residue in fixed_struct.get_residues() for atom in residue]
    #         all_atoms_moving = [atom for residue in moving_struct.get_residues() for atom in residue]
    #         if len(all_atoms_fixed) > len(all_atoms_moving):
    #             dif = len(all_atoms_fixed) - len(all_atoms_moving)
    #             all_atoms_fixed = all_atoms_fixed[:len(all_atoms_fixed) - dif]
    #         else:
    #             dif = len(all_atoms_moving) - len(all_atoms_fixed)
    #             all_atoms_moving = all_atoms_moving[:len(all_atoms_moving) - dif]
    #             # If dif is positive, you might need to address how to handle this case,
    #             # such as choosing which atoms to remove from fixed_crys_atoms or reconsidering alignment.
    #             pass
    #     assert len(all_atoms_fixed) == len(all_atoms_moving)
    #     rmsd = self.calculate_rmsd(all_atoms_fixed, all_atoms_moving, chain)
    #     self.dict[key_name] = rmsd