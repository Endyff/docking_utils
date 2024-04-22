import argparse
import pandas as pd
from pathlib import Path
from rdkit import Chem
from spyrmsd.rmsd import symmrmsd
import numpy as np
from collections import Counter
from posebusters import PoseBusters


argparser = argparse.ArgumentParser(description='Evaluate docking results')
argparser.add_argument('--csv', default='data/simple_val_pred.csv', help='CSV file with docking results')
argparser.add_argument('-o',  '--output', default='out.csv', help='Output CSV file')
argparser.add_argument('--pb_check', default=False, action='store_true', help='Check if the pose is correct using PoseBusters')

def load_molecule_from_sdf(path):
    mol = Chem.SDMolSupplier(path, sanitize=False)[0]
    mol = Chem.RemoveHs(mol)
    return mol

def get_spyrmsd_format(mol):
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    positions = mol.GetConformer().GetPositions()
    adj = Chem.GetAdjacencyMatrix(mol)
    return elements, positions, adj
     

def calculate_symmrmsd(mol_pred, mol_true):
    elem_pred, coords_pred, adj_pred = get_spyrmsd_format(mol_pred)
    elem_true, coords_true, adj_true = get_spyrmsd_format(mol_true)
    c1 = Counter(elem_pred)
    c2 = Counter(elem_true)
    if c1 != c2:
        print(f'Elements are different: {c1} vs {c2}')
        return np.nan
    if coords_pred.shape != coords_true.shape or adj_pred.shape != adj_true.shape:
        raise ValueError(f'Different number of atoms or bonds.')

    rmsd = symmrmsd(coords_pred, coords_true, elem_pred, elem_true, adj_pred, adj_true)[0]
    return rmsd

def calculate_geometric_center_dist(mol1, mol2):
    coords1 = mol1.GetConformer().GetPositions()
    coords2 = mol2.GetConformer().GetPositions()
    return np.linalg.norm(coords1.mean(axis=0) - coords2.mean(axis=0))


def posebusters_redocking_check(path_preds, path_conds, mode='dock'):

    if mode not in ['dock', 'redock', 'mol', 'gen']:
        raise ValueError(f'Invalid mode: {mode}')
     
    pb = PoseBusters(config=mode)

    pb_df = pd.DataFrame({'mol_cond': path_conds, 'mol_pred': path_preds})
    pb_results = pb.bust_table(pb_df).reset_index()
    # print(pb_results)
    return pb_results



if __name__ == '__main__':
    args = argparser.parse_args()

    data_df = pd.read_csv(args.csv)
    out_df = data_df.copy()
    rmsds = []
    for i, row in data_df.iterrows():

        path_true, path_pred, path_cond = row['path_true'], row['path_pred'], row['path_cond']
        mol_true = load_molecule_from_sdf(row['path_true'])
        mol_pred = load_molecule_from_sdf(row['path_pred'])
        rmsd = calculate_symmrmsd(mol_pred, mol_true)
        rmsds.append(rmsd)
        mean_rmsd = calculate_geometric_center_dist(mol_pred, mol_true)

    if args.pb_check:
        pb = PoseBusters()
        pb_df = data_df.copy()
        pb_df.columns = ['mol_cond', 'mol_true', 'mol_pred']
        # print(pb_df)
        pb_df.to_csv('pb_df.csv', index=False)
        pb_results = pb.bust_table(pb_df).reset_index()
        print(pb_results)
        pb_results.to_csv('pb_results.csv', index=False)
        print(set(pb_results['file']) == set(out_df['path_pred']))

        out_df = out_df.merge(pb_results, left_on='path_pred', right_on='file')
        print(out_df)
    out_df['rmsd'] = rmsds
    out_df.to_csv(args.output, index=False)
    # data_df.to_csv('posebusters_benchmark_set_diffdock_rmsd.csv', index=False)
