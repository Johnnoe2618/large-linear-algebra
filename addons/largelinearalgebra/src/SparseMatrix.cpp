#include "SparseMatrix.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/error_macros.hpp>

using namespace godot;

void SparseMatrixCoordinate::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_row", "row"), &SparseMatrixCoordinate::set_row);
    ClassDB::bind_method(D_METHOD("set_column", "column"), &SparseMatrixCoordinate::set_column);
    ClassDB::bind_method(D_METHOD("set_value", "value"), &SparseMatrixCoordinate::set_value);

    ClassDB::bind_method(D_METHOD("get_row"), &SparseMatrixCoordinate::get_row);
    ClassDB::bind_method(D_METHOD("get_column"), &SparseMatrixCoordinate::get_column);
    ClassDB::bind_method(D_METHOD("get_value"), &SparseMatrixCoordinate::get_value);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "row"), "set_row", "get_row");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "column"), "set_column", "get_column");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "value"), "set_value", "get_value");

    ClassDB::bind_method(D_METHOD("clone"), &SparseMatrixCoordinate::clone);
}

void SparseMatrix::_bind_methods() {
    ClassDB::bind_static_method("SparseMatrix", D_METHOD("identity", "size"), &SparseMatrix::identity);
	ClassDB::bind_static_method("SparseMatrix", D_METHOD("zero", "size"), &SparseMatrix::zero);

    ClassDB::bind_method(D_METHOD("reserve_per_row", "number_of_columns_reserved"), &SparseMatrix::reserve_per_row);
    ClassDB::bind_method(D_METHOD("set_dimensions", "rows", "columns"), &SparseMatrix::set_dimensions);
    ClassDB::bind_method(D_METHOD("get_dimensions"), &SparseMatrix::get_dimensions);
    ClassDB::bind_method(D_METHOD("set_element", "row", "column", "value"), &SparseMatrix::set_element);
    ClassDB::bind_method(D_METHOD("get_element", "row", "column"), &SparseMatrix::get_element);
    ClassDB::bind_method(D_METHOD("transposed"), &SparseMatrix::transposed);
    ClassDB::bind_method(D_METHOD("clone"), &SparseMatrix::clone);
    ClassDB::bind_method(D_METHOD("to_dense"), &SparseMatrix::to_dense);
    ClassDB::bind_method(D_METHOD("to_coordinate_array"), &SparseMatrix::to_coordinate_array);

    ClassDB::bind_method(D_METHOD("multiply_sparse", "other"), &SparseMatrix::multiply_sparse);
    ClassDB::bind_method(D_METHOD("multiply_sparse_to_dense", "other"), &SparseMatrix::multiply_sparse_to_dense);
    ClassDB::bind_method(D_METHOD("multiply_dense", "other"), &SparseMatrix::multiply_dense);
    ClassDB::bind_method(D_METHOD("multiply_vector", "vector"), &SparseMatrix::multiply_vector);

    ClassDB::bind_method(D_METHOD("add_sparse", "other"), &SparseMatrix::add_sparse);
    ClassDB::bind_method(D_METHOD("subtract_sparse", "other"), &SparseMatrix::subtract_sparse);
    ClassDB::bind_method(D_METHOD("add_dense", "other"), &SparseMatrix::add_dense);
    ClassDB::bind_method(D_METHOD("subtract_dense", "other"), &SparseMatrix::subtract_dense);

    ClassDB::bind_method(D_METHOD("multiply_scaler_in_place", "scalar"), &SparseMatrix::multiply_scaler_in_place);
    ClassDB::bind_method(D_METHOD("solve_iterative_cg", "B", "initial_guess", "max_iterations"), &SparseMatrix::solve_iterative_cg);

    ClassDB::bind_method(D_METHOD("norm_squared"), &SparseMatrix::norm_squared);
    ClassDB::bind_method(D_METHOD("norm"), &SparseMatrix::norm);
    
}

Ref<SparseMatrix> SparseMatrix::identity(int size) {
	if (unlikely(size < 0)) {
		ERR_PRINT_ED("ERROR: cannot create matrix of negative size.");
		Ref<SparseMatrix> null_ref;
		return null_ref;
	}
    Ref<SparseMatrix> identity;
    identity.instantiate();
    identity->set_dimensions(size, size);
    for (int i = 0; i < identity->values.size(); i++) {
        Value new_value;
        new_value.column = i;
        new_value.val = 1.0;
        identity->values[i].push_back(new_value);
    }

    return identity;
}

Ref<SparseMatrix> SparseMatrix::zero(int size) {
	if (unlikely(size < 0)) {
		ERR_PRINT_ED("ERROR: cannot create matrix of negative size.");
		Ref<SparseMatrix> null_ref;
		return null_ref;
	}
    Ref<SparseMatrix> zero_matrix;
    zero_matrix.instantiate();
    zero_matrix->set_dimensions(size, size);

    return zero_matrix;
}

void SparseMatrix::set_dimensions(int rows, int cols) {
	if (unlikely(rows < 0 || cols < 0)) {
		ERR_PRINT_ED("ERROR: cannot set dimensions to negative.");
		return;
	}
    // truncate/expand the rows
    values.resize(rows);

    if (cols < column_number) {
        for (int i = 0; i < values.size(); i++) {
            int end_iterator = values[i].size() - 1;

            for (;;) {
                if (end_iterator == 0 || values[i][end_iterator].column < cols) break;
            }

            values[i].resize(end_iterator); // truncate to the last column that is within the dimensions
        }
    } // else, nothing is required,

    column_number = cols;
}

void SparseMatrix::reserve_per_row(int number_columns_per_row) {
    for (auto& item : values) item.reserve(number_columns_per_row);
}

void SparseMatrix::set_element(int row, int col, double value) {
    if (unlikely(row >= values.size() || col >= column_number || row < 0 || col < 0)) {
        ERR_PRINT_ED("ERROR: cannot set cell outside dimension range.");
        return;
    }

    // set cell code here: iterator from the last element to the first element (reverse order), and place the record so that the list remains sorted (TODO: may see performance improvement with binary search - test this later). 
    int end_iterator = values[row].size() - 1;
    for (;;) {
        if (end_iterator == -1) {
            if (value == 0.0) break;
            // if == -1, then the column is less than everything in the list, and must be inserted at the start
            Value new_value;
            new_value.val = value;
            new_value.column = col;
            values[row].insert(values[row].begin(), new_value);
            break;
        } else if (values[row][end_iterator].column < col) {
            if (value == 0.0) break;
            // if column < col, then we are currently looking at the element that is right behind the place we want it to be. Insert in front.
            Value new_value;
            new_value.val = value;
            new_value.column = col;
            values[row].insert(values[row].begin() + (end_iterator + 1), new_value);
            break;
        } else if (values[row][end_iterator].column == col) {
            if (value == 0.0) {
                // delete the value from the array
                values[row].erase(values[row].begin() + end_iterator);
                break;
            }
            // if column == col, then we are currently looking at the element we want to edit, so edit it.
            values[row][end_iterator].val = value;
            break;
        }
    }
}

double SparseMatrix::get_element(int row, int col) const {
    if (unlikely(row >= values.size() || col >= column_number || row < 0 || col < 0)) {
        ERR_PRINT_ED("ERROR: cannot get cell outside dimension range.");
        return 0.0;
    }
    // find the cell, use binary search or linear search (TODO: test the options for performance, right now it is linear search)
    for (auto it = values[row].begin(); it != values[row].end(); ++it) {
        if (it->column == col) return it->val;
    }
    return 0.0;
}

Ref<SparseMatrix> SparseMatrix::clone() const {
    Ref<SparseMatrix> copy;
    copy.instantiate();
    // set dimensions
    copy->column_number = column_number;
    copy->values.resize(values.size());

    for (int i = 0; i < values.size(); ++i) {
        copy->values[i].resize(values[i].size());
        const Value* value_pointer = values[i].data();
        Value* destination_pointer = copy->values[i].data();
        memcpy(destination_pointer, value_pointer, values[i].size() * sizeof(Value));
    }
    return copy;
}

TypedArray<SparseMatrixCoordinate> SparseMatrix::to_coordinate_array() const {
    TypedArray<SparseMatrixCoordinate> result;

    for (int row = 0; row < values.size(); row++) {
        for (auto& column : values[row]) {
            Ref<SparseMatrixCoordinate> new_coordinate;
            new_coordinate.instantiate();
            new_coordinate->row = row;
            new_coordinate->column = column.column;
            new_coordinate->value = column.val;

            result.append(new_coordinate);
        }
    }

    return result;
}

Ref<SparseMatrix> SparseMatrix::transposed() const {
    Ref<SparseMatrix> transposed;
    transposed.instantiate();
    
    transposed->column_number = values.size();
    transposed->values.resize(column_number);
    transposed->reserve_per_row(8); // default to 8, so that reallocations are minimised. (though, more reallocations may occur, this reduces the number of them.)

    for (int row = 0; row < values.size(); row++) {
        for (const auto& column : values[row]) {
            Value new_value;
            new_value.val = column.val;
            new_value.column = row;
            transposed->values[column.column].push_back(new_value);
        }
    }

    return transposed;
}

Ref<SparseMatrix> SparseMatrix::multiply_sparse(Ref<SparseMatrix> other) const {
    if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
        Ref<SparseMatrix> null_ref;
        return null_ref;
	}
    if (unlikely(column_number != other->values.size())) {
        ERR_PRINT_ED("ERROR: cannot multiply matrices with incompatible dimensions. Left-hand row length must equal Right-hand column length.");
        Ref<SparseMatrix> null_ref;
        return null_ref;
    }
    // Algorithm: take 'other' and transpose it. Then, based on the new transposed matrix, do row-dot-row (which mimicks row-dot-column, due to transposing).
    // The complexity is MxNxl time, where M and N are the dimensions, and l is the number of numbers in each row (assuming that l is on average the same value).
    Ref<SparseMatrix> transposed = other->transposed();
    Ref<SparseMatrix> result;
    result.instantiate();
    result->set_dimensions(values.size(), other->column_number);
    result->reserve_per_row(8); // reserve 8 to minimise reallocations, though reallocations will still occur

    for (int row = 0; row < result->values.size(); row++) {
        for (int col = 0; col < result->column_number; col++) {
            // Then do the dot product of the relevant rows.
            auto a = values[row].begin();
            auto b = transposed->values[col].begin();

            auto end_a = values[row].end();
            auto end_b = transposed->values[col].end();

            double dot_product = 0.0;

            while (a != end_a && b != end_b) {
                if (a->column == b->column) {
                    dot_product += a->val * b->val;
                    a++; b++; // progress forward
                } else if(a->column < b->column) {
                    a++; // if a is less then b, then iterate on a to see if it can 'catch up'
                } else { // b->column < a->column
                    b++;
                }
            }

            if (dot_product != 0.0) {
                Value new_value;
                new_value.column = col;
                new_value.val = dot_product;
                result->values[row].push_back(new_value);
            }
        }
    }

    return result;
}

Ref<DenseMatrix> SparseMatrix::multiply_sparse_to_dense(Ref<SparseMatrix> other) const {
    if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
        Ref<DenseMatrix> null_ref;
        return null_ref;
	}
    if (unlikely(column_number != other->values.size())) {
        ERR_PRINT_ED("ERROR: cannot multiply matrices with incompatible dimensions. Left-hand row length must equal Right-hand column length.");
        Ref<DenseMatrix> null_ref;
        return null_ref;
    }

    Ref<DenseMatrix> result;
    result.instantiate();
    result->set_dimensions(values.size(), other->column_number);

    Ref<SparseMatrix> transposed = other->transposed();

    for (int row = 0; row < result->row_number; row++) {
        for (int col = 0; col < result->col_number; col++) {
            // Then do the dot product of the relevant rows.
            auto a = values[row].begin();
            auto b = transposed->values[col].begin();

            auto end_a = values[row].end();
            auto end_b = transposed->values[col].end();

            double dot_product = 0.0;

            while (a != end_a && b != end_b) {
                if (a->column == b->column) {
                    dot_product += a->val * b->val;
                    a++; b++; // progress forward
                } else if(a->column < b->column) {
                    a++; // if a is less then b, then iterate on a to see if it can 'catch up'
                } else { // b->column < a->column
                    b++;
                }
            }

            result->get_cell_unsafe(row, col) = dot_product;
        }
    }

    return result;
}

void SparseMatrix::multiply_scaler_in_place(double scaler) {
    for (int row = 0; row < values.size(); row++) {
        for (auto& value : values[row]) {
            value.val *= scaler;
        }
    }
}

Ref<DenseMatrix> SparseMatrix::multiply_dense(Ref<DenseMatrix> other) const {
    if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
        Ref<DenseMatrix> null_ref;
        return null_ref;
	}
    if (unlikely(column_number != other->row_number)) {
        ERR_PRINT_ED("ERROR: cannot multiply matrices with incompatible dimensions. Left-hand row length must equal Right-hand column length.");
        Ref<DenseMatrix> null_ref;
        return null_ref;
    }

    Ref<DenseMatrix> result;
	result.instantiate();
	result->set_dimensions(values.size(), other->col_number);

    for (int row = 0; row < values.size(); row++) {
		for (int col = 0; col < other->col_number; col++) {
			double total = 0.0;
			
			for (auto& it : values[row]) {
				total += it.val * other->values[other->col_number * it.column + col];
			}

			result->values[result->col_number * row + col] = total;
		}
	}

    return result;
}

Ref<SparseMatrix> SparseMatrix::add_sparse(Ref<SparseMatrix> other) const {
    if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
        Ref<SparseMatrix> null_ref;
        return null_ref;
	}	
    if (unlikely(values.size() != other->values.size() || column_number != other->column_number)) {
        ERR_PRINT_ED("ERROR: cannot add matrices with dimensions that are not equal.");
        Ref<SparseMatrix> null_ref;
        return null_ref;
    }

    Ref<SparseMatrix> result;
    result.instantiate();
    result->set_dimensions(values.size(), column_number);

    for (int i = 0; i < values.size(); i++) {
        auto left_iterator = values[i].begin();
        auto right_iterator = other->values[i].begin();

        auto left_end = values[i].end();
        auto right_end = values[i].end();

        result->values[i].reserve(values[i].size() + other->values[i].size()); // reserve to minimise re-allocations

        for (;;) {
            if (left_iterator == left_end || right_iterator == right_end) {
                // one or both iterators have finished traversal.
                break; // go into the next loops to complete the left or right traversal remainder
            }
            if (left_iterator->column < right_iterator->column) {
                // if left is behind, then add it to the result, and move forward
                Value result_val;
                result_val.column = left_iterator->column;
                result_val.val = left_iterator->val;
                result->values[i].push_back(result_val);
                left_iterator++;
            } else if (left_iterator->column == right_iterator->column) {
                // if they are equal, add both together
                Value result_val;
                result_val.column = left_iterator->column;
                result_val.val = left_iterator->val + right_iterator->val;
                result->values[i].push_back(result_val);
                left_iterator++;
                right_iterator++;
            } else {
                // left_iterator->column > right_iterator->column: right is behind. Add it, and move forward.
                Value result_val;
                result_val.column = right_iterator->column;
                result_val.val = right_iterator->val;
                result->values[i].push_back(result_val);
                right_iterator++;
            }
        }
        if (left_iterator != left_end) { // complete the left traversal
            while (left_iterator != left_end) {
                Value result_val;
                result_val.column = left_iterator->column;
                result_val.val = left_iterator->val;
                result->values[i].push_back(result_val);
                left_iterator++;
            }
        }
        if (right_iterator != right_end) { // complete the right traversal
            while (right_iterator != left_end) {
                Value result_val;
                result_val.column = right_iterator->column;
                result_val.val = right_iterator->val;
                result->values[i].push_back(result_val);
                right_iterator++;
            }
        }
    }

    return result;
}
Ref<SparseMatrix> SparseMatrix::subtract_sparse(Ref<SparseMatrix> other) const {
    if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
        Ref<SparseMatrix> null_ref;
        return null_ref;
	}	
    if (unlikely(values.size() != other->values.size() || column_number != other->column_number)) {
        ERR_PRINT_ED("ERROR: cannot subtract matrices with dimensions that are not equal.");
        Ref<SparseMatrix> null_ref;
        return null_ref;
    }

    Ref<SparseMatrix> result;
    result.instantiate();
    result->set_dimensions(values.size(), column_number);

    for (int i = 0; i < values.size(); i++) {
        auto left_iterator = values[i].begin();
        auto right_iterator = other->values[i].begin();

        auto left_end = values[i].end();
        auto right_end = values[i].end();

        result->values[i].reserve(values[i].size() + other->values[i].size()); // reserve to minimise re-allocations

        for (;;) {
            if (left_iterator == left_end || right_iterator == right_end) {
                // one or both iterators have finished traversal.
                break; // go into the next loops to complete the left or right traversal remainder
            }
            if (left_iterator->column < right_iterator->column) {
                // if left is behind, then add it to the result, and move forward
                Value result_val;
                result_val.column = left_iterator->column;
                result_val.val = left_iterator->val;
                result->values[i].push_back(result_val);
                left_iterator++;
            } else if (left_iterator->column == right_iterator->column) {
                // if they are equal, add both together
                Value result_val;
                result_val.column = left_iterator->column;
                result_val.val = left_iterator->val - right_iterator->val;
                result->values[i].push_back(result_val);
                left_iterator++;
                right_iterator++;
            } else {
                // left_iterator->column > right_iterator->column: right is behind. Add it, and move forward.
                Value result_val;
                result_val.column = right_iterator->column;
                result_val.val = -right_iterator->val;
                result->values[i].push_back(result_val);
                right_iterator++;
            }
        }
        if (left_iterator != left_end) { // complete the left traversal
            while (left_iterator != left_end) {
                Value result_val;
                result_val.column = left_iterator->column;
                result_val.val = left_iterator->val;
                result->values[i].push_back(result_val);
                left_iterator++;
            }
        }
        if (right_iterator != right_end) { // complete the right traversal
            while (right_iterator != left_end) {
                Value result_val;
                result_val.column = right_iterator->column;
                result_val.val = -right_iterator->val;
                result->values[i].push_back(result_val);
                right_iterator++;
            }
        }
    }

    return result;
}

Ref<DenseMatrix> SparseMatrix::add_dense(Ref<DenseMatrix> other) const {
    if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
        Ref<DenseMatrix> null_ref;
        return null_ref;
	}
	if (unlikely(other->row_number != values.size() || other->col_number != column_number)) {
        ERR_PRINT_ED("ERROR: cannot add matrices with dimensions that are not equal.");
        Ref<DenseMatrix> null_ref;
        return null_ref;
    }

	Ref<DenseMatrix> result = other->clone(); // clone values (many of them won't be touched. Only a few will be editted by the sparse matrix)
	for (int row = 0; row < values.size(); row++) {
		for (auto& col : values[row]) {
			result->values[column_number * row + col.column] += col.val;
		}
	}

	return result;
}
Ref<DenseMatrix> SparseMatrix::subtract_dense(Ref<DenseMatrix> other) const {
    if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
        Ref<DenseMatrix> null_ref;
        return null_ref;
	}
	if (unlikely(other->row_number != values.size() || other->col_number != column_number)) {
        ERR_PRINT_ED("ERROR: cannot subtract matrices with dimensions that are not equal.");
        Ref<DenseMatrix> null_ref;
        return null_ref;
    }

	Ref<DenseMatrix> result = other->clone();
    result->multiply_scaler_in_place(-1.0); // set them all negative, as other is on the RHS of the subtraction

	for (int row = 0; row < other->row_number; row++) {
		for (auto& col : values[row]) {
			result->values[column_number * row + col.column] += col.val;
		}
	}

	return result;
}

double SparseMatrix::norm_squared() const {
    double total = 0.0;
    for (int i = 0; i < values.size(); i++) {
        for (auto& col : values[i]) {
            total += col.val * col.val;
        }
    }

    return total;
}

double SparseMatrix::norm() const {
    return sqrt(norm_squared());
}

Ref<VectorN> SparseMatrix::multiply_vector(Ref<VectorN> other) const {
    if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
        Ref<VectorN> null_ref;
        return null_ref;
	}
	if (unlikely(column_number != other->values.size())) {
        ERR_PRINT_ED("ERROR: cannot multiply matrices with vectors of a length not equal to the matrix column number.");
        Ref<VectorN> null_ref;
        return null_ref;
    }

    Ref<VectorN> result;
	result.instantiate();
	result->values.resize(values.size());

    for (int row = 0; row < values.size(); row++) {
        double total = 0.0;

        for (auto& col : values[row]) {
            total += col.val * other->values[col.column];
        }

        result->values[row] = total;
    }

    return result;
}

Ref<DenseMatrix> SparseMatrix::to_dense() const {
    Ref<DenseMatrix> result;
    result.instantiate();
    result->set_dimensions(values.size(), column_number);

    for (int row = 0; row < values.size(); row++) {
        for (auto& col : values[row]) {
            result->values[column_number * row + col.column] = col.val;
        }
    }

    return result;
}

// This is identical to the DenseMatrix version algorithm
Ref<VectorN> SparseMatrix::solve_iterative_cg(Ref<VectorN> B, Ref<VectorN> initial_guess, int max_iterations) const {
	if (unlikely(B.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
        Ref<VectorN> null_ref;
        return null_ref;
	}
	if (unlikely(values.size() != column_number)) {
		ERR_PRINT_ED("ERROR: matrix must be square.");
        Ref<VectorN> null_ref;
        return null_ref;
	}

	if (initial_guess.is_null()) {
		initial_guess = VectorN::filled(0.0, B->values.size());
	}

	if (max_iterations > column_number) max_iterations = column_number; // the conjugate method is theoretically perfect after N steps
    if (max_iterations < 0) max_iterations = column_number;

	// algorithm was stolen from Wikipedia: (https://en.wikipedia.org/wiki/Conjugate_gradient_method)

	Ref<VectorN> X = initial_guess->clone();
	Ref<VectorN> R = B->subtract(multiply_vector(X));
	Ref<VectorN> P = R->clone();

	for (int i = 0; i < max_iterations; i++) {
		if (R->is_approximately_zero(std::numeric_limits<double>::epsilon())) break;

		Ref<VectorN> Ap = multiply_vector(P);
		double RR = R->dot(R);
		double a = RR / P->dot(Ap);
		X->add_multiple_in_place(P, a);
		R->add_multiple_in_place(Ap, -a);

		double b = R->dot(R) / RR;
		P->multiply_scalar_in_place(b);
		P->add_in_place(R);
	}

	// return the guess
	return X;
}