# Charge change set
FEP input structures and results for the charge change set. The charge change set contains 11 test cases (CDK2, DLK, EGFR, EPHX2, IRAK4\_S2, IRAK4\_S3, ITK, JAK1, JNK1, PTP1B, Thrombin, TYK2) originally from DOI [10.1021/acs.jctc.8b00825](https://pubs.acs.org/doi/10.1021/acs.jctc.8b00825) and DOI [10.1021/acs.jctc.1c00302](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00302).

## Data
Each test case in each data set has the following folders and files:

**ligands:** contains XFF parameter files for ligand molecules in AMBER format (mol2, frcmod and leaprc file)

**proteins:** contains protein PDB files

**map:** contains atom mapping file for each pair of transformation. The atom mapping file is a two-column file. Each row corresponds to the indices (1-based) of a pair of mapped atoms for ligand1 (left) and ligand2 (left).

**XXX\_raw\_edge.csv:** raw FEP ddG results for each pair. All free energy values are in unit of kcal/mol.
Contents of columns:
- Pairs: transformation pair name
- EXP: experimental ddG value
- FEP: raw RBFE ddG value
- FEP\_err: raw RBFE ddG error
- CCC: RBFE ddG value after cycle-closure correction
- CCC\_err: cycle-closure correction error

**XXX\_pka_edge.csv:** pKa corrected FEP ddG results for each pair. ddG values are corrected by the cycle-closure correction method before pKa corrrection. All free energy values are in unit of kcal/mol.
Contents of columns:
- Pairs: transformation pair name
- EXP: experimental ddG value
- FEP: RBFE ddG value after cycle-closure and pKa correction
- FEP\_err: raw RBFE ddG error: RBFE ddG error

**XXX\_pka_dg.csv:** FEP dG results from pka corrected ddG values. All free energy values are in unit of kcal/mol.
Contents of columns:
- lig: ligand name
- EXP: experimental dG value
- Pred dG: predicted dG value
- Pred dG err: predicted dG error

**XXX\_ion\_raw\_edge.csv:** raw FEP ddG results with additional ions for each pair

**XXX\_ion\_pka\_edge.csv:** pKa corrected FEP ddG results with additional ions for each pair. ddG values are corrected by the cycle-closure correction method before pKa corrrection.

**XXX\_ion\_pka\_dg.csv:** FEP dG results with additional ions from pka corrected ddG values

For Thrombin, since no pair needs pKa correction, only the following files exist:

**Thrombin\_raw\_edge.csv:** raw FEP ddG results for each pair

**Thrombin\_raw\_dg.csv:** FEP dG results

**Thrombin\_ion\_raw\_edge.csv:** raw FEP ddG results with additional ions for each pair

**Thrombin\_ion\_raw\_dg.csv:** FEP dG results with additional ions

The starting structure for each ligand can be found in the mol2 file that has the same name as the ligand.

