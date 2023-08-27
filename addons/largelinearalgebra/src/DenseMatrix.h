#ifndef DENSE_MATRIX_H
#define DENSE_MATRIX_H

#include "SparseMatrix.h"
#include "VectorN.h"
#include <godot_cpp/classes/ref.hpp>
#include <vector>

namespace godot {

// Note: processed for row major
class DenseMatrix : public RefCounted {
    GDCLASS(DenseMatrix, RefCounted)

private:
    int row_number = 0;
	int col_number = 0;

	std::vector<double> values;
    
    // helper functions for matrix internal manipulations
    _FORCE_INLINE_ void blit_rect_unsafe(const DenseMatrix& source, Rect2i src_rect, Vector2i dst); // x is the row, y is the column
    _FORCE_INLINE_ void blit_rect_unsafe_ref(Ref<DenseMatrix> source, Rect2i src_rect, Vector2i dst);
    _FORCE_INLINE_ void swap_rows(int row_a, int row_b);
    _FORCE_INLINE_ void multiply_row_after_col(int row, double coefficient, int first_col);
    _FORCE_INLINE_ void add_row_multiple_after_col(int row_dest, int row_src, double multiple, int first_col);
    _FORCE_INLINE_ double& get_cell_unsafe(int row, int col) { return values[col_number * row + col]; };

    
protected:
    static void _bind_methods();

public:
    friend class SparseMatrix;
    friend class VectorN;

    DenseMatrix() {};
    ~DenseMatrix() {};

    // static methods
    static Ref<DenseMatrix> identity(int size);
    static Ref<DenseMatrix> zero(int size);
    static Ref<DenseMatrix> from_packed_array(const PackedFloat64Array& from, int rows, int columns);

    // data manipulation methods
    void set_dimensions(int rows, int cols); 
    Vector2i get_dimensions() const { return Vector2i(row_number, col_number); }
    void set_element(int row, int col, double value);
    double get_element(int row, int col) const;
    Ref<DenseMatrix> transposed() const;
    Ref<DenseMatrix> clone() const;
    PackedFloat64Array to_packed_array() const;

    Ref<SparseMatrix> to_sparse(double zero_threshold) const;

    // Algebra
    Ref<DenseMatrix> multiply_dense(Ref<DenseMatrix> other) const;
    Ref<DenseMatrix> multiply_sparse(Ref<SparseMatrix> other) const;
    Ref<VectorN> multiply_vector(Ref<VectorN> other) const;

    Ref<DenseMatrix> add_sparse(Ref<SparseMatrix> other) const;
    Ref<DenseMatrix> subtract_sparse(Ref<SparseMatrix> other) const;
    Ref<DenseMatrix> add_dense(Ref<DenseMatrix> other) const;
    Ref<DenseMatrix> subtract_dense(Ref<DenseMatrix> other) const;

    void multiply_scaler_in_place(double scaler);
    void add_sparse_in_place(Ref<SparseMatrix> other);
    void subtract_sparse_in_place(Ref<SparseMatrix> other);
    void add_dense_in_place(Ref<DenseMatrix> other);
    void subtract_dense_in_place(Ref<DenseMatrix> other);

    // Solves: AX=B, where B is the parameter, and solves for X (for any matrix X, including column vectors, or rectangular), using Gaussian Elimination
    Ref<DenseMatrix> solve(Ref<DenseMatrix> B) const; 

    // Solve iteratively with the conjugate gradient method. This is guaranteed to converge for symmetric positive-definite matrices (but can often converge otherwise aswell). 
    // Note: it is relatively easy for end-users to create their own iterative solvers with the given functions, as the GDScript overhead is dominated by the matrix algebra.
    Ref<VectorN> solve_iterative_cg(Ref<VectorN> B, Ref<VectorN> initial_guess, int max_iterations) const;

    Ref<DenseMatrix> inverse() const;

    //data
    double norm_squared() const;
    double norm() const;
};

}

#endif