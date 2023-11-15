import os
import glob
import pandas as pd
import numpy as np
import sys
import pickle

import simtk.openmm as mm
from simtk.openmm import app
from simtk.openmm import *
from simtk.openmm.app import *

import parmed
from parmed import unit as u

from rdkit import Chem
from rdkit.Chem import TorsionFingerprints, AllChem

from multiprocessing import Pool

from math_utils import bat_value, calcrms

from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField

class MolConfig:
    def __init__(self, file=None, restraints=None, output=None, ik=None, char_folder=None):
        self.file = file
        self.output = output
        self.restraints = restraints
        self.ik = ik
        self.char_folder = char_folder

    def createOpenffMolFromMolFile(self, mol_file):

        rdmol = Chem.MolFromMolFile(mol_file, removeHs=False)
        off_mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
        conf = rdmol.GetConformer(0)
        coords = conf.GetPositions() / 10.0
        return off_mol, coords

    def getCoordsFromFile(self):
        rdmol = Chem.MolFromMolFile(self.file, removeHs=False)
        conf = rdmol.GetConformer(0)
        coords = conf.GetPositions() / 10.0
        return coords

    def calcTFD(self, opt_file):
        que_mol = Chem.MolFromPDBFile(opt_file, removeHs=False)
        ref_mol = Chem.MolFromMolFile(self.file, removeHs=False)
        ref_mol.AddConformer(que_mol.GetConformer(), assignId=True)
        try:
            tfd = TorsionFingerprints.GetTFDBetweenConformers(ref_mol, confIds1=[0], confIds2=[1])[0]
        # triggered for molecules such as urea
        except:
            tfd = np.nan
        return tfd

    def min_ffxml(self, ffxml):
        try:
            self.offmol, self.coords = self.createOpenffMolFromMolFile(self.file)

        #try:
            # load in force field
            ff = ForceField(ffxml)

            # Load charge mol
            char_mol = pickle.load(open(os.path.join(self.char_folder, self.ik+".pickle"), 'rb'))

            # create components for OpenMM system
            topology = Topology.from_molecules(molecules=[self.offmol])

            # create openmm system ==> prone to triggering Exception
            #system = ff.create_openmm_system(topology, charge_from_molecules=[off_mol])
            system = ff.create_openmm_system(topology, charge_from_molecules=[char_mol], allow_nonintegral_charges=True)

            # minimize structure with ffxml
            energy = self.run_openmm(topology, system, self.coords, self.restraints)
            # Calculate structure RMS
            rdmol1 = Chem.MolFromMolFile(self.file)
            basename = os.path.basename(self.file).replace(".mol","")
            output_file = os.path.join(self.output, basename+".pdb")
            rdmol2 = Chem.MolFromPDBFile(output_file)
            rms = calcrms(rdmol1, rdmol2)
        
        except Exception:
            print('Unable to min with openff')
            return None, None

        return energy, rms

    def min_xff(self):

        amber_home=os.environ.get('AMBERHOME')
        if amber_home is not None:
            tleap = os.path.join(amber_home, 'bin/tleap')
        else:
            print("Environment variable AMBERHOME cannot be found!")
            exit()

        # Generate frcmod file
        fname = os.path.basename(self.file).split(".")[0]
        fbname = fname.split('_')[0]
        frcmod = os.path.join(self.char_folder, fbname + ".frcmod")
        mol2 = os.path.join(self.char_folder, fbname + ".mol2")
        prmtop = fname + ".prmtop"
        inpcrd = fname + ".inpcrd"
        leaprc = os.path.join(self.char_folder, fbname+"_leaprc.dat")

        # construct inpcrd & prmtop
        tleap_tmpl_with_leaprc = """
        source {}
        MOL = loadmol2 {}
        loadamberparams {}
        saveamberparm MOL {} {}
        quit

        """

        try:
            tleap_cmd = tleap_tmpl_with_leaprc.format(leaprc, mol2, frcmod, prmtop, inpcrd)
            leap_in_file = "{}.leap.in".format(fname)
            leap_log = "{}.log".format(fname)
            with open(leap_in_file, "w") as fw:
                fw.write(tleap_cmd)
            os.system("{} -s -f {} > {}".format(tleap, leap_in_file, leap_log))

            # create OpenMM system from inpcrd & prmtop
            prmtop_a = AmberPrmtopFile(prmtop)
            inpcrd_a = AmberInpcrdFile(inpcrd)
            system = prmtop_a.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1*u.nanometer, constraints=None)
            positions = self.getCoordsFromFile()

            energy = energy = self.run_openmm(prmtop_a.topology, system, positions, self.restraints, openff=False)
            # tfd = self.calcTFD(self.output+"/"+fname+".pdb")

            # Cleaning
            os.system("rm "+leap_in_file+" "+inpcrd+" "+prmtop+" "+leap_log)

            # Calculate structure RMS
            rdmol1 = Chem.MolFromMolFile(self.file)
            basename = os.path.basename(self.file).replace(".mol","")
            output_file = os.path.join(self.output, basename+".pdb")
            rdmol2 = Chem.MolFromPDBFile(output_file)
            rms = calcrms(rdmol1, rdmol2)

        except:
            print("Unable to min with XFF!")
            # Cleaning
            os.system("rm "+leap_in_file+" "+inpcrd+" "+prmtop+" "+leap_log)
            return None, None

        return energy, rms

    def min_gaff(self, gaff2=False):

        if gaff2:
            opt_ff = "gaff2"
        else:
            opt_ff = "gaff"

        amber_home=os.environ.get('AMBERHOME')
        if amber_home is not None:
            tleap = os.path.join(amber_home, 'bin/tleap')

        # Generate frcmod file
        fname = os.path.basename(self.file).split(".")[0]
        fbname = fname.split('_')[0]
        frcmod = os.path.join(self.char_folder, fbname + ".frcmod")
        mol2 = os.path.join(self.char_folder, fbname + ".mol2")
        prmtop = fname + ".prmtop"
        inpcrd = fname + ".inpcrd"
        leaprc = os.path.join(self.char_folder, fbname+"_leaprc.dat")

        # construct inpcrd & prmtop
        tleap_tmpl_with_leaprc = """
        source leaprc.{}
        MOL = loadmol2 {}
        check MOL
        loadamberparams {}
        saveamberparm MOL {} {}
        quit

        """

        try:
            tleap_cmd = tleap_tmpl_with_leaprc.format(opt_ff, mol2, frcmod, prmtop, inpcrd)
            leap_log = "{}.log".format(fname)
            leap_in_file = "{}.leap.in".format(fname)
            with open(leap_in_file, "w") as fw:
                fw.write(tleap_cmd)
            os.system("{} -s -f {} > {}".format(tleap, leap_in_file, leap_log))

            # create OpenMM system from inpcrd & prmtop
            prmtop_a = AmberPrmtopFile(prmtop)
            inpcrd_a = AmberInpcrdFile(inpcrd)
            system = prmtop_a.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1*u.nanometer, constraints=None)
            positions = self.getCoordsFromFile()

            energy = self.run_openmm(prmtop_a.topology, system, positions, self.restraints, openff=False)
            # tfd = self.calcTFD(self.output+"/"+fname+".pdb")

            # Cleaning
            os.system("rm "+leap_in_file+" "+inpcrd+" "+prmtop+" "+leap_log)

            # Calculate structure RMS
            rdmol1 = Chem.MolFromMolFile(self.file)
            basename = os.path.basename(self.file).replace(".mol","")
            output_file = os.path.join(self.output, basename+".pdb")
            rdmol2 = Chem.MolFromPDBFile(output_file)
            rms = calcrms(rdmol1, rdmol2)

        except:
            print("Unable to min with {}!".format(opt_ff))
            # Cleaning
            os.system("rm "+leap_in_file+" "+inpcrd+" "+prmtop+" "+leap_log)
            return None, None

        return energy, rms

    def run_openmm(self, topology, system, positions, restraints, openff=True):

        if restraints:
            # Add dihedral restraints
            current_torsion_restrain = CustomTorsionForce("k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = 3.1415926535")
            current_torsion_restrain.addGlobalParameter("k", 40000.0)
            dihe = bat_value(positions[restraints])
            current_torsion_restrain.addPerTorsionParameter('theta0')
            current_torsion_restrain.addTorsion(restraints[0], restraints[1], restraints[2], restraints[3], [dihe])
            sysForce = system.addForce(current_torsion_restrain)

        # need to create integrator but don't think it's used
        integrator = mm.LangevinIntegrator(
                300.0 * u.kelvin,
                1.0 / u.picosecond,
                2.0 * u.femtosecond)

        # create simulation object then minimize
        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        simulation.minimizeEnergy(tolerance=1.0E-8, maxIterations=3000)
        if restraints: system.removeForce(sysForce)

        # get minimized positions
        positions = simulation.context.getState(getPositions=True).getPositions()
        basename = os.path.basename(self.file).replace(".mol","")
        output_file = open(os.path.join(self.output, basename+".pdb"), 'w')
        if openff:
            omm_topology = simulation.topology.to_openmm()
        else:
            omm_topology = simulation.topology
        PDBFile.writeFile(omm_topology, positions, output_file)

        # get minimized energy
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        energy = energy.value_in_unit(u.kilocalories_per_mole)

        return energy

def run(mol_file, ff, restraints, output_folder, ik, id, char_folder, ffn):
    try:
        mol = MolConfig(file=mol_file, restraints=restraints, output=output_folder, ik=ik, char_folder=char_folder)
    except:
        return None, None, None, None
    if ff == "openff":
        ene, rms = mol.min_ffxml(ffn)
    elif ff == "gaff":
        ene, rms = mol.min_gaff(gaff2=False)
    elif ff == "gaff2":
        ene, rms = mol.min_gaff(gaff2=True)
    elif ff == "xff":
        ene, rms = mol.min_xff()
    else:
        print("Unsupported force field!")
        exit() 
    if ene:
        return ene, ik, id, rms
    else:
        return None, None, None, None

def main():
    from multiprocessing import Pool
    import argparse

    # Arg parse
    desc = '''Scripts to run OpenFF validation.
'''

    help_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                            formatter_class=help_formatter)
    #parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument("-i", help='path to the input file')
    required.add_argument("-f", help='force field, choose among openff, xff, gaff, and gaff2. Note if xff is chosen, mol2/frcmod/leaprc file \
                                      should be prepared in advance and put into the charge folder.')
    required.add_argument("-s", help='folder containing input structures')
    required.add_argument("-o", help='folder for output structures')
    required.add_argument("-c", help='folder containing charge files')
    required.add_argument("-r", help='number of cores')

    optional.add_argument("-n", help='If use openff, version of openff')

    opt = parser.parse_args()

    fin = opt.i
    ff = opt.f
    fin_folder = opt.s
    fout = opt.o
    char_folder = opt.c
    cores = int(opt.r)
    ffn = opt.n

    # Input validity check
    if not os.path.exists(fin): 
        print("Input file {} does not exist!".format(fin))
        exit()
    if not ff in ("openff", "xff", "gaff", "gaff2"):
        print("Force field {} is not supported!".format(ff))
        exit()
    if not os.path.exists(fin_folder):
        print("Input structure folder {} does not exist!".format(fin_folder))
        exit()
    if not os.path.exists(char_folder):
        print("Charge file folder {} does not exist!".format(char_folder))
        exit()
    if not os.path.exists(fout): os.system("mkdir -p "+fout)

    data = pd.read_csv(fin)
    jobs = []
    pool = Pool(processes=cores)
    iks = []
    ids = []
    enes = []
    rmsds = []
    for row in data.iterrows():
        ik = str(row[1]["ik"])
        id = str(row[1]["id"])
        dih = row[1]["dih"]
        restr = [int(x) for x in dih.split("_")]
        mol_file = os.path.join(fin_folder, ik + "_" + id + ".mol")
        print(ik+"\t"+id)
    #    ene, ik, id, rms = run(mol_file, ff, restr, fout, ik, id, char_folder, ffn)
    #    if ene:
    #        iks.append(ik)
    #        ids.append(id)
    #        enes.append(ene)
    #        rmsds.append(rms)
        jobs.append(pool.apply_async(run,(mol_file, ff, restr, fout, ik, id, char_folder, ffn)))

    for job in jobs:
        re = job.get()
        ene, ik, id, rms = re
        if ene:
            iks.append(ik)
            ids.append(id)
            enes.append(ene)
            rmsds.append(rms)

    data_out = pd.DataFrame({"ik":iks, "id":ids, "{}-ene".format(ff):enes, "rmsd":rmsds})
    data_out.to_csv("{}-result.csv".format(ff))

def test():
    mol = MolConfig(file="acs2016/AAALMKOPEYUEBW-YAQRNVERNA-N_4556759.mol", restraints=None,output="./", ik="AAALMKOPEYUEBW-YAQRNVERNA-N", char_folder="acs2016-offchar")
    print(mol.min_ffxml("openff_unconstrained-2.0.0.offxml"))

if __name__=="__main__":
    main()
    #test()
