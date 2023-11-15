# Example Amber input files
Example Amber .in files used in FEP simulations

## Content
**complex\_equil.in:** complex equilibration in solvent

**ligand\_equil.in:** ligand equilibration in solvent

**complex\_ti\_min.in:** complex TI minimization

**complex\_ti\_heat.in:** complex TI heat equilibration

**complex\_remd.in:** complex REMD production

**solvated\_ti\_min.in:** solvated TI minimization

**solvated\_ti\_heat.in:** solvated TI heat equilibration

**solvated\_remd.in:** solvated REMD production

Please note: The exact input files for each tranformation pair are not shown here since the input files alone cannot completely reflect our FEP workflow. For example, our enhanced sampling algorithm requires setting "gti\_add\_sc" to 3. However, this is not the whole algorithm and we used additional codes to handle the scaling of torsion parameters in the soft-core region. The purpose of showing the inpput files is to display some common parameters we used in the Amber simulations.
