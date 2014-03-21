StructLDLT: Approximating matrices with multiple symmetries
=========

This git repo contains the slides and Matlab code explaining StructLDLT.

Abstract

In many matrix applications it is possible to utilize structure for improved computational efficiency and obtain a low-rank approximation which preserves the original structure. We seek to translate efficient, low-rank, structure preserving properties to tensors, or high dimensional matrices. After unfolding a tensor into a matrix, tensor computations turn into matrix computations. When the original tensor is structured, the matrix unfolding is also structured. Motivated by a problem in computational quantum chemistry which involves a four level nested summation over a low-rank symmetric 4-tensor, we obtain a structured matrix unfolding and apply StructLDLT: a rank-revealing lazy-evaluation symmetric-pivoting LDL^T factorization algorithm which preserves block symmetry, symmetric blocks, and perfect shuffle permutation symmetry. StructLDLT is implemented in Matlab and allows us to compute the quantum chemistry summation in just O(r^2n^2) where r is the Kronecker rank of the symmetric tensor. This is joint work with Charles Van Loan. 