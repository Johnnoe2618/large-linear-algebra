# Changelog
This is high-level changelog for each release of this extension.

## Version 1.0.0
Initial Release

### Features:
- Sparse and Dense Matrix representations with double precision arithmetic.
- Efficient Matrix-Matrix and Matrix-Vector multiplication for Dense and Sparse Matrices.
- Heterogenous arithmetic between different Matrix representations and Vectors, including Multiplication, Addition, and Subtraction.
- Routines for the conversion between Dense and Sparse Matrices, and Godot built-in PackedFloat64Array types.
- Both In-place and Not-in-place variants of a number of arithmetical operations, for efficiency. Some operations are in-place only, and Matrices should be explicitly cloned if that is the desired behaviour.
- Direct solving using Gaussian Elimination with Partial Pivoting for Dense Matrices
- Matrix inverses calculated using Gaussian Elimination with Partial Pivoting for Dense Matrices
- Iterative solvers for Dense and Sparse Matrices, using the Conjugate Gradient Method.

## Version 1.0.1
- Bug Fix (DenseMatrix.solve/DenseMatrix.inverse): pivot selection error when a column contains no positive numbers.
- Bug Fix: Conversion from DenseMatrix to SparseMatrix now correctly works with negatives.