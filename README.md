# scripts

Collection of scripts written to help with my research

## confint.py

Program to calculate estimates and confidence intervals of portions
(uses normal approximation, Clopper-Pearson, Aggresti-Coull, Wilson score with and without continouity correction)

Arguments: number_of_hits sample_size conf_int_percentage

the last argument is optional, 95 is assumed if missing

Example:

``python3 confint.py 32 48 99``

This means 32 positive cases in 48 data points, 99% conf. interval requested.

See the code (or a textbook) for comments on the various methods.


## rasx2xy.py

Extracts the diffraction profile from a Rigaku {basename}.rasx file to {basename}.xy file. Only tries to extract a single (i.e. the first) profile.

Type ``python3 rasx2xy.py -h`` to see command line options.

## splitmol2.pl

Splits a mol2 file into separate files for each residue (requires perl5).

``./splitmol2.pl name.mol2`` will create name-1.mol2, name-2.mol2 ..., containing
residues 1, 2, ... of the original molecule, respectively. It is desinged to work with single-molecule mol2 files.
