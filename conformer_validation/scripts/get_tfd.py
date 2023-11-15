from multiprocessing import Pool
import os
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import TorsionFingerprints


class MolConfig:
    def __init__(self, f1, f2):
        self.f1 = f1
        self.f2 = f2
    '''
    def calcTFD(self):
        mol1 = Molecule.from_file(self.f1)
        mol2 = Molecule.from_file(self.f2)
        tors_list, ring_tors_list = TorsionFingerprints.CalculateTorsionLists(mol1.rdmol)
        weights = TorsionFingerprints.CalculateTorsionWeights(mol1.rdmol)
        torsions1 = CalculateTorsionAngles(mol1.rdmol, tors_list, ring_tors_list)
        torsions2 = CalculateTorsionAngles(mol2.rdmol, tors_list, ring_tors_list)
        return CalculateTFD(torsions1, torsions2, weights=weights)
    '''
    def calcTFD(self):
        que_mol = Chem.MolFromPDBFile(self.f2, removeHs=False)
        ref_mol = Chem.MolFromMolFile(self.f1, removeHs=False)
        ref_mol.AddConformer(que_mol.GetConformer(), assignId=True)
        try:
            tfd = TorsionFingerprints.GetTFDBetweenConformers(ref_mol, confIds1=[0], confIds2=[1])[0]
        # triggered for molecules such as urea
        except:
            tfd = TorsionFingerprints.GetTFDBetweenConformers(ref_mol, confIds1=[0], confIds2=[1], useWeights=False)[0]
        return tfd

def run(f1, f2, ik, id, ene, rmsd):
    mol = MolConfig(f1, f2)
    tfd = mol.calcTFD()
    return ik, id, ene, rmsd, tfd

def main():
    fin, input_folder, output_folder = sys.argv[1:4]
    
    iks = []
    ids = []
    enes = []
    rmsds = []
    tfds = []
    
    jobs = []
    pool = Pool(processes=36)
    all_configs = pd.read_csv(fin)
    for row in all_configs.iterrows():
        ik = str(row[1]["ik"])
        id = str(row[1]["id"])
        #ene = row[1]["off-ene"]
        ene = row[1]['openff-ene']
        rmsd = row[1]['openff-rmsd']
        f1 = os.path.join(input_folder, ik+"_"+id+".mol")
        f2 = os.path.join(output_folder, ik+"_"+id+".pdb")
        jobs.append(pool.apply_async(run, args=(f1, f2, ik, id, ene, rmsd)))

    for job in jobs:
        ik, id, ene, rmsd, tfd = job.get()
        iks.append(ik)
        ids.append(id)
        enes.append(ene)
        rmsds.append(rmsd)
        tfds.append(tfd)
    data = pd.DataFrame({"ik":iks, "id":ids, "openff-ene":enes, "openff-rmsd":rmsds, "openff-tfd":tfds})
    data.to_csv("acs2016-openff_unconstrained-paper-result-tfd.csv")

if __name__=="__main__":
    main()
