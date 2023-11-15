# Conformer energy and geometry validation
Molecule and conformer IDs for the Fragment Set and their XFF results

## Data
**Fragment-Set-result.csv:** contains information for molecules that are in the public domain in the Fragment Set and their XFF energies, RMSDs, and TFDs.
Contents of columns:
- mol\_id: molecule id
- conf\_id: conformer id
- dih: the indices of four atoms (0-based) for the torsion profile this conformer belongs to
- deg: the torsion angle (defined by the four atoms) value of this conformer, unit degree
- qm-ene: QM energy, unit kcal/mol
- xff-ene: XFF energy, unit kcal/mol
- xff-rmsd: XFF RMSD, unit Angstrom
- xff-tfd: XFF TFD

**input_structures.zip:** contains input structures for conformers of these molecules in the Fragment Set

**xff-params.zip:** XFF parameters for these molecules
