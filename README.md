This package, `SL2Reps`, provides methods for constructing and testing matrix
presentations of the representations of SL(2,Z) whose kernels are congruence subgroups of SL(2,Z).

Irreducible representations of prime-power level are constructed individually by using the Weil representations of quadratic modules, and from these a list of all representations of a given degree or level can be produced. The format is designed for the study of modular tensor categories in particular, providing
symmetric matrix presentations of each representation.

## Installation

To install `SL2Reps`, first download it from `https://snw-0.github.io/sl2-reps/`, then place it in the `pkg`
subdirectory of your GAP installation (or in the `pkg` subdirectory of any other GAP
root directory, for example one added with the `-l` argument).

`SL2Reps` is then loaded with the GAP command

`gap> LoadPackage( "SL2Reps" );`

## Contact

Siu-Hung Ng, Louisiana State University (`rng@math.lsu.edu`)
<br> Yilong Wang[^1], Louisiana State University (`wyl@bimsa.cn`)
<br> Samuel Wilson, Louisiana State University (`swil311@lsu.edu`)

[^1]: Currently at Beijing Institute of Mathematical Science and Applications (BIMSA).

## Acknowledgements

This project is partially supported by NSF grant DMS 1664418.
