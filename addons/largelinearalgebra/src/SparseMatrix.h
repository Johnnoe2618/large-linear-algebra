#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "DenseMatrix.h"
#include "VectorN.h"
#include <godot_cpp/classes/ref.hpp>
#include <vector>

namespace godot {

template<class T> class move_only : public T{
public:
   move_only(){}
   move_only(const move_only&) = delete;
   move_only(move_only&&) noexcept {};
   ~move_only() noexcept {};

   using T::T;   
}; /* Enforces the moving of std::vector when reallocating values, rather than copying */

struct Value {
    double val;
    int column;
};

class SparseMatrixCoordinate : public RefCounted {
    GDCLASS(SparseMatrixCoordinate, RefCounted)

private:
    int row = 0;
    int column = 0;
    double value = 0.0; 

protected:
    static void _bind_methods();

public:
    friend class SparseMatrix;

    SparseMatrixCoordinate() {};
    ~SparseMatrixCoordinate() {};

    void set_row(int new_row) { row = new_row; }
    int get_row() { return row; }

    void set_column(int new_column) { column = new_column; }
    int get_column() { return column; }

    void set_value(double new_value) { value = new_value; }
    double get_value() { return value; }

    Ref<SparseMatrixCoordinate> clone() {
        Ref<SparseMatrixCoordinate> cloned;
        cloned.instantiate();

        cloned->row = row;
        cloned->column = column;
        cloned->value = value;

        return cloned;
    }

};

// Note: processed for row major
class SparseMatrix : public RefCounted {
    GDCLASS(SparseMatrix, RefCounted)

private:
    int column_number = 0;
    std::vector<move_only<std::vector<Value>>> values; // values[row] => { value, column }
    
protected:
    static void _bind_methods();

public:
    friend class DenseMatrix;
    friend class VectorN;

    SparseMatrix() {};
    ~SparseMatrix() {};

    // static methods
    static Ref<SparseMatrix> identity(int size);
    static Ref<SparseMatrix> zero(int size);

    // data manipulation methods
    void reserve_per_row(int number_columns_per_row);
    void set_dimensions(int rows, int cols); 
    Vector2i get_dimensions() const { return Vector2i(values.size(), column_number); }
    void set_element(int row, int col, double value);
    double get_element(int row, int col) const;
    Ref<SparseMatrix> transposed() const;
    Ref<SparseMatrix> clone() const;
    TypedArray<SparseMatrixCoordinate> to_coordinate_array() const;

    Ref<DenseMatrix> to_dense() const;

    // algebra methods
    Ref<SparseMatrix> multiply_sparse(Ref<SparseMatrix> other) const;
    Ref<DenseMatrix> multiply_sparse_to_dense(Ref<SparseMatrix> other) const;
    Ref<DenseMatrix> multiply_dense(Ref<DenseMatrix> other) const;
    Ref<VectorN> multiply_vector(Ref<VectorN> other) const;

    Ref<SparseMatrix> add_sparse(Ref<SparseMatrix> other) const;
    Ref<SparseMatrix> subtract_sparse(Ref<SparseMatrix> other) const;
    Ref<DenseMatrix> add_dense(Ref<DenseMatrix> other) const;
    Ref<DenseMatrix> subtract_dense(Ref<DenseMatrix> other) const;

    void multiply_scaler_in_place(double scaler);

    Ref<VectorN> solve_iterative_cg(Ref<VectorN> B, Ref<VectorN> initial_guess, int max_iterations) const;

    //data
    double norm_squared() const;
    double norm() const;
};



}

#endif