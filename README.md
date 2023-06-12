# Anti-Fouling Peptide Prediction

*Author: Clyde Overby*


## Installation

Install using 
```sh
$ pip install -e .
```
## Usage

After installation, calling
```sh
$ AFPepPro
```
Will initate the script, upon which it will as for an input file and a destination file.

## Overview

This is a program intended to screen for semi-randomly constructed peptides' ability to resist self-interaction based on their charge sequence.  Evaluation is performed between discrete members of the total group of possible peptides from a single charge sequence.  Input parameters are the allowable types (+,-, neutral) of amino acids, the type's abundance, and the probability of specific amino acids for each type.  Testing is performed via the PASTA algorithm, with the assumption that peptides will be relatively unstructured while anchored on a soft surface.

Example input file:

#number_AAs (+/-/0), probability of charges (+/-/0), (+/-/0) [Amino acids, probability of amino acid incorporation], sequence length, # of charge motifs to generate, # of Monte Carlo iterations to perform per charge motif 
1, 2, 2, .4, .4, .2, K, 1, D, E, .5, .5, S, T, .5, .5, 20, 5000, 200


```sh
###Index:
(integer)
###V0:
(integer or float)
###constant:
(integer or float)
###bs_coefficients:
(integer)
###basis set:
(char) either l for Legendre or f for Fourier
###domain:
(float, float) (-1, 1) for nearly all cases
```

Output file for peptide-peptide mode is in the format of charge motif, minimum interaction energy average of tested peptides for the charge motif, and standard deviation of interaction energy.

Output file format for peptide-protein interactions is a matrix of peptide-protein interaction energy minimums.

(c) 2023
