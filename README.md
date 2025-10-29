# eigen_derivatives

This Matlab package is an implementation of the pseudocode presented in chapter 3 of the PhD thesis

*Uncertainty Quantification of Elliptic Eigenvalue Problems with Degenerate Eigenvalues*. 2026.
**David Christoph Ebert**.
University of Bonn.
doi.org/TBA .

No additional Matlab toolboxes are required.

## Description

This package calculates derivatives of eigenpairs (eigenvalues and eigenfunctions) given the derivatives of matrices. 
For eigenpairs to degenerate eigenvalues derivatives with respect to the eigenspace and in the traditional sense can be calculated.
In the degenerate case the derivatives of the eigenpairs in the traditional sense are calculated via the polarization matrix and its derivatives.
Eigenvalue problems (EVPs) can be given as standard matrix EVPs or as generalized EVPs, i.e., as a "stiffness" and "mass" matrix. 
The implementation supports dense of sparse matrices as input. The input matrices of the EVP are assumed to be symmetric.

## Contents

To include this package use

`addpath *PATH*/eigen_derivatives`

### Functions

| File              | Details
| :-----------------| :------ |
| `group_eig.m`     | groups eigenvalues into degenerate sets
| `eig_der.m`       | derivatives of eigenpair (with respect to eigenspace if multiplicity m>1)
| `pol.m`           | (initial) polarization matrix given derivatives of eigenvalues with respect to eigenspace, throws warning if initial polarization cannot be fully determined by input
| `pol_der.m`       | derivatives of polarization matrix
| `du_pol.m`        | polarized eigenfunction derivatives given derivatives with respect to eigenspace and polarization (initial and derivatives)
| *auxiliary* |
| `multiindexsum.m` | compiles list of multiindices of certain length that sum to a give order
| `multinom.m`      | multinomial coefficients for multiindex given output of multiindexsum.m

### Examples

To access examples, please navigate to folder "examples" and launch scripts "example*.m".
The script "demo.m" runs the above functions for each example and visualizes a Taylor expansion of the eigenvalues.

| File                       | Details
| :------------------------- | :------ |
| `example_nondegenerate.m`  | non-degenerate eigenvalue
| `example_cross_friswell.m` | degenerate eigenvalues (multiplicity 2) that form cones, such that the trajectories cross, first derivatives of eigenvalues differ
| `example_deflect.m`        | degenerate eigenvalues (multiplicity 2), such that the trajectories deflect, first derivatives of eigenvalues coincide, second derivatives differ
| `example_cross_cross.m`    | degenerate eigenvalues (multiplicity 3), such that all trajectories cross pairwise
| `example_cross_deflect.m`  | degenerate eigenvalues (multiplicity 3), such that two trajectories deflect while another crosses both
| `example_pair_cross.m`     | degenerate eigenvalues (multiplicity 4), such that the eigenvalues form two deflecting pairs, which as pairs cross each other

