# Merck KGaA set
FEP input structures and results for the Merck KGaA set. The Merck KGaA set contains 8 test cases (CDK8, Eg5, HIF2a, c-Met, SYK, TNKS2, PFKFB3, SHP-2) originally from DOI [10.1021/acs.jcim.0c00900](https://pubs.acs.org/doi/10.1021/acs.jcim.0c00900).

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

## Note
For PFKFB3 test case, pfkfb3\_XXX is the results with the pyrophosphate and citrate anion cofactors. PFKFB3\_no\_cofac\_XXX is the results without cofactors.

For TNKS2 test case, the following files are present in the folder:

**TNKS2\_4UI5\_edges.csv:** FEP ddG results for each pair using structure 4UI5

**TNKS2\_4UI5\_dg.csv:** FEP dG results for each ligand using structure 4UI5

**TNKS2\_3MHK\_raw\_edges.csv:** raw FEP ddG results for each pair using structure 3MHK

**TNKS2\_3MHK\_edges.csv:** FEP ddG results after pKa correction for each pair using structure 3MHK. ddG values are corrected by cycle-closure correction before pKa correction

**TNKS2\_3MHK\_dg.csv:** FEP dG results from pKa corrected ddG values using structure 3MHK


