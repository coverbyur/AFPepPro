# Anti-Fouling Peptide Prediction

*Author: Clyde Overby*


## Installation

Install using 
```sh
$ pip install setup.py
```
## Usage

After installation, calling
```sh
$ AFPep
```
Will initate the script, upon which it will as for an input file and a destination file.

## Overview

This is a program intended to screen for semi-randomly constructed peptides' ability to resist self-interaction based on their charge sequence.  Evaluation is performed between discrete members of the total group of possible peptides from a single charge sequence.  Input parameters are the allowable types (+,-, neutral) of amino acids, the type's abundance, and the probability of specific amino acids for each type.  Testing is performed via the PASTA algorithm, with the assumption that peptides will be relatively unstructured while anchored on a soft surface.

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

Output file is in the form of index: [coefficients] where the coefficients are those that satisfy the Schrodinger equation for the basis set chosen.

(c) 2019