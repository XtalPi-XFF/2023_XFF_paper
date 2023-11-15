# R group substitution set
FEP input structures and results for the R group substitution set. The R group substitution set contains 8 test cases (Bace, CDK2, Jnk1, MCL1, p38, PTP1B, Thrombin, Tyk2) originally from DOI [10.1021/ja512751q](https://pubs.acs.org/doi/10.1021/ja512751q).

## Data
Each test case in each data set has the following folders and files:

**ligands:** contains XFF parameter files for ligand molecules in AMBER format (mol2, frcmod and leaprc file)

**proteins:** contains protein PDB files

**map:** contains atom mapping file for each pair of transformation. The atom mapping file is a two-column file. Each row corresponds to the indices (1-based) of a pair of mapped atoms for ligand1 (left) and ligand2 (left).

**XXX\_edges.csv:** FEP ddG results for each pair. All free energy values are in unit of kcal/mol.
Contents of columns:
- Pairs: transformation pair name
- EXP: experimental ddG value
- FEP: raw RBFE ddG value
- FEP\_err: raw RBFE ddG error
- CCC: RBFE ddG value after cycle-closure correction
- CCC\_err: cycle-closure correction error

**XXX\_dg.csv:** FEP dG results for each ligand. All free energy values are in unit of kcal/mol.
Contents of columns:
- lig: ligand name
- EXP: experimental dG value
- Pred dG: predicted dG value
- Pred dG err: predicted dG error

The starting structure for each ligand can be found in the mol2 file that has the same name as the ligand.

