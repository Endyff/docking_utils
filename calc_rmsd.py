import argparse
import pandas as pd
from rdkit import Chem
import numpy as np
from spyrmsd.rmsd import symmrmsd

parser = argparse.ArgumentParser(description='Calculate RMSD from CSV containing paths, columns are \'path_pred\' and \'path_true\'')
parser.add_argument('-i', '--input_file', type=str, help='Input CSV file')
parser.add_argument('-o', '--output_file', type=str, help='Output CSV file')
parser.add_argument('--rmsd_type', type=str, choices=['all_atom', 'mean'], default='all_atom', help='RMSD calculation type')
args = parser.parse_args()

data_df = pd.read_csv(args.input_file)

# Load the reference protein and ligand
def calculate_symmrmsd(p1, p2, mean=False):
        
    m1 = Chem.MolFromMolFile(str(p1))
    m2 = Chem.MolFromMolFile(str(p2))

    amref = Chem.rdmolops.GetAdjacencyMatrix(m1)
    am = Chem.rdmolops.GetAdjacencyMatrix(m2)

    for c in m1.GetConformers():
        coordsref = np.asarray(c.GetPositions())

    for c in m2.GetConformers():
        coords = np.asarray(c.GetPositions())

    apropsref = np.asarray([x.GetSymbol() for x in m1.GetAtoms()])
    aprops = np.asarray([x.GetSymbol() for x in m2.GetAtoms()])
    if mean:
        return np.linalg.norm(coordsref.mean(axis=0) - coords.mean(axis=0))
    else:
        return symmrmsd(
            coordsref=coordsref,
            coords=coords,
            apropsref=apropsref,
            aprops=aprops,
            amref=amref,
            am=am,
            )

rmsds = []
mean_rmsds = []
for i, row in data_df.iterrows():
    print(i)
    # try:
    mol_pred = Chem.MolFromMolFile(row['path_pred'])
    mol_true = Chem.MolFromMolFile(row['path_true'])

    
    cur_rmsd, _ = calculate_symmrmsd(row['path_pred'], row['path_true'], False)
    code = str(row['path_pred']).split('/')[-2]
    rmsds.append(cur_rmsd)

    # except Exception as e:
    #     print(f'Error calculating RMSD for {row["path_pred"]} and {row["path_true"]}: {e}')
    #     rmsds.append(np.nan)
    
    # try:
    #     mean_rmsd, _ = calculate_symmrmsd(row['path_pred'], row['path_true'], True)
    #     mean_rmsds.append(mean_rmsd)
    # except Exception as e:
    #     mean_rmsds.append(np.nan)



# if args.rmsd_type == 'mean':
data_df['mean_rmsd'] = rmsds
# elif args.rmsd_type == 'all_atom':
data_df['all_atom_rmsd'] = mean_rmsds

data_df.to_csv(args.output_file, index=False)