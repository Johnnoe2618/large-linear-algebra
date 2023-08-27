# large linear algebra
 Large Linear Algebra (GDExtension)

This is an extension written for Godot 4.1+ in C++.
The extension provides efficient representations for the algebraic manipulation of **arbitrarily sized** Dense Matrices, Sparse Matrices, and Vectors with double precision floating point arithmetic. The interface is accessible directly from GDScript, using the GDExtension defined types *DenseMatrix*, *SparseMatrix*, and *VectorN*.

## Features:
- Sparse and Dense Matrix representations with double precision arithmetic.
- Efficient Matrix-Matrix and Matrix-Vector multiplication for Dense and Sparse Matrices.
- Heterogenous arithmetic between different Matrix representations and Vectors, including Multiplication, Addition, and Subtraction.
- Routines for the conversion between Dense and Sparse Matrices, and Godot built-in PackedFloat64Array types.
- Both In-place and Not-in-place variants of a number of arithmetical operations, for efficiency. Some operations are in-place only, and Matrices should be explicitly cloned if that is the desired behaviour.
- Direct solving using Gaussian Elimination with Partial Pivoting for Dense Matrices
- Matrix inverses calculated using Gaussian Elimination with Partial Pivoting for Dense Matrices
- Iterative solvers for Dense and Sparse Matrices, using the Conjugate Gradient Method.

## Documentation
Documentation can be found on the [GitHub Wiki](https://github.com/Johnnoe2618/large-linear-algebra/wiki)

## Compilation of Dynamically Linked Libraries from Source
If you are using this extension on a platform where precompiled binaries have not been provided (only provided currently for windows), you must compile them yourself from this repository.
The SConstruct file for building the extension is within addons/largelinearalgebra. It is a derivative of the template provided from the [Godot Documentation](https://docs.godotengine.org/en/stable/tutorials/scripting/gdextension/gdextension_cpp_example.html).

To compile the binaries, go to the extension directory in the command prompt, and run ``scons platform=[platform] target=[template_debug|template_release]``. Then, for the resulting binary in the ``bin`` directory, add the platform to the ``largelinearalgebra.gdextension`` file.