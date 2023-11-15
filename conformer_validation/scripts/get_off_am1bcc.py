import os
import glob
import pandas as pd
import numpy as np
import sys
import pickle

from openff.toolkit.topology import Molecule, Topology

from rdkit import Chem

from multiprocessing import Pool

class MolConfig:
    def __init__(self, file=None, ik=None, output=None):
        self.file = file
        self.ik = ik
        self.output = output

    def createOpenffMolFromMolFile(self):
        rdmol = Chem.MolFromMolFile(self.file, removeHs=False)
        off_mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
        off_mol.compute_partial_charges_am1bcc()
        pickle.dump(off_mol, open(self.output+"/"+self.ik+".pickle", 'wb'))


def run(mol_file, ik, output_folder):
    mol = MolConfig(file=mol_file, ik=ik, output=output_folder)
    try:
        mol.createOpenffMolFromMolFile()
    except:
        pass

def main():
    input_folder, output_folder = sys.argv[1:3]
    pool = Pool(processes=4)
    # We actually used the lowest energy conformer for each molecule
    all_configs = glob.glob(input_folder+"/*.mol")
    jobs = []
    for f in all_configs:
        ik = os.path.basename(f)
        ik = ik.split(".")[0]
        jobs.append(pool.apply_async(run, args=(f,ik,output_folder)))

    for job in jobs:
        _ = job.get()

if __name__=="__main__":
    main()
