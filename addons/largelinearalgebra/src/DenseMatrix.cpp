#include "DenseMatrix.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/error_macros.hpp>
#include <string>

using namespace godot;

void DenseMatrix::_bind_methods() {
	ClassDB::bind_static_method("DenseMatrix", D_METHOD("identity", "size"), &DenseMatrix::identity);
	ClassDB::bind_static_method("DenseMatrix", D_METHOD("zero", "size"), &DenseMatrix::zero);
	ClassDB::bind_static_method("DenseMatrix", D_METHOD("from_packed_array", "array", "rows", "columns"), &DenseMatrix::from_packed_array);

	ClassDB::bind_method(D_METHOD("set_dimensions", "rows", "columns"), &DenseMatrix::set_dimensions);
	ClassDB::bind_method(D_METHOD("get_dimensions"), &DenseMatrix::get_dimensions);
	ClassDB::bind_method(D_METHOD("set_element", "row", "column", "value"), &DenseMatrix::set_element);
	ClassDB::bind_method(D_METHOD("get_element", "row", "column"), &DenseMatrix::get_element);
	ClassDB::bind_method(D_METHOD("transposed"), &DenseMatrix::transposed);
	ClassDB::bind_method(D_METHOD("clone"), &DenseMatrix::clone);
	ClassDB::bind_method(D_METHOD("to_sparse", "zero_threshold"), &DenseMatrix::to_sparse);
	ClassDB::bind_method(D_METHOD("to_packed_array"), &DenseMatrix::to_packed_array);	

	ClassDB::bind_method(D_METHOD("multiply_dense", "other"), &DenseMatrix::multiply_dense);
	ClassDB::bind_method(D_METHOD("multiply_sparse", "other"), &DenseMatrix::multiply_sparse);
	ClassDB::bind_method(D_METHOD("multiply_vector", "vector"), &DenseMatrix::multiply_vector);

	ClassDB::bind_method(D_METHOD("add_sparse", "other"), &DenseMatrix::add_sparse);
	ClassDB::bind_method(D_METHOD("subtract_sparse", "other"), &DenseMatrix::subtract_sparse);
	ClassDB::bind_method(D_METHOD("add_dense", "other"), &DenseMatrix::add_dense);
	ClassDB::bind_method(D_METHOD("subtract_dense", "other"), &DenseMatrix::subtract_dense);

	ClassDB::bind_method(D_METHOD("multiply_scaler_in_place", "scalar"), &DenseMatrix::multiply_scaler_in_place);
	ClassDB::bind_method(D_METHOD("add_sparse_in_place", "other"), &DenseMatrix::add_sparse_in_place);
	ClassDB::bind_method(D_METHOD("subtract_sparse_in_place", "other"), &DenseMatrix::subtract_sparse_in_place);
	ClassDB::bind_method(D_METHOD("add_dense_in_place", "other"), &DenseMatrix::add_dense_in_place);
	ClassDB::bind_method(D_METHOD("subtract_dense_in_place", "other"), &DenseMatrix::subtract_dense_in_place);

	ClassDB::bind_method(D_METHOD("solve", "B"), &DenseMatrix::solve);
	ClassDB::bind_method(D_METHOD("solve_iterative_cg", "B", "initial_guess", "max_iterations"), &DenseMatrix::solve_iterative_cg);
	ClassDB::bind_method(D_METHOD("inverse"), &DenseMatrix::inverse);

	ClassDB::bind_method(D_METHOD("norm_squared"), &DenseMatrix::norm_squared);
	ClassDB::bind_method(D_METHOD("norm"), &DenseMatrix::norm);
}

/* Helper Functions */

void DenseMatrix::blit_rect_unsafe(const DenseMatrix& source, Rect2i src_rect, Vector2i dst) {
	int copy_length = src_rect.size.y;
	for (int row = 0; row < src_rect.size.x; row++) {
		int source_row = src_rect.position.x + row;
		int dest_row = dst.x + row; // destination 
		int source_col = src_rect.position.y;
		int dest_col = dst.y;
		memcpy(&values[col_number * dest_row + dest_col], &source.values[source.col_number * source_row + source_col], sizeof(double) * copy_length);
	}
}

void DenseMatrix::blit_rect_unsafe_ref(Ref<DenseMatrix> source, Rect2i src_rect, Vector2i dst) {
	int copy_length = src_rect.size.y;
	for (int row = 0; row < src_rect.size.x; row++) {
		int source_row = src_rect.position.x + row;
		int dest_row = dst.x + row; // destination 
		int source_col = src_rect.position.y;
		int dest_col = dst.y;
		memcpy(&values[col_number * dest_row + dest_col], &source->values[source->col_number * source_row + source_col], sizeof(double) * copy_length);
	}
}

void DenseMatrix::swap_rows(int row_a, int row_b) {
	for (int i = 0; i < col_number; i++) {
		double a = values[col_number * row_a + i];
		values[col_number * row_a + i] = values[col_number * row_b + i];
		values[col_number * row_b + i] = a;
	}
}

void DenseMatrix::multiply_row_after_col(int row, double coefficient, int first_col) {
	for (int i = first_col; i < col_number; i++) {
		values[col_number * row + i] *= coefficient;
	}
}

void DenseMatrix::add_row_multiple_after_col(int row_dest, int row_src, double multiple, int first_col) {
	for (int i = first_col; i < col_number; i++) {
		values[col_number * row_dest + i] += values[col_number * row_src + i] * multiple;
	}
}

/* Public Functions */

Ref<DenseMatrix> DenseMatrix::from_packed_array(const PackedFloat64Array& from, int rows, int columns) {
	if (unlikely(rows < 0 || columns < 0 || from.size() != rows * columns)) {
		ERR_PRINT_ED("ERROR: invalid dimensions.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}

	Ref<DenseMatrix> result;
	result.instantiate();
	result->values.resize(from.size());

	memcpy(result->values.data(), from.ptr(), sizeof(double) * from.size());

	result->row_number = rows;
	result->col_number = columns;
	return result;
}

PackedFloat64Array DenseMatrix::to_packed_array() const {
	PackedFloat64Array result;
	result.resize(values.size());
	for (int i = 0; i < values.size(); i++) {
		result[i] = values[i];
	}

	return result;
}

Ref<DenseMatrix> DenseMatrix::identity(int size) {
	if (unlikely(size < 0)) {
		ERR_PRINT_ED("ERROR: cannot create matrix of negative size.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	Ref<DenseMatrix> identity;
	identity.instantiate();

	identity->row_number = size;
	identity->col_number = size;
	identity->values.resize(size * size, 0.0);
	for (int i = 0; i < size; i++) {
		identity->values[size * i + i] = 1.0;
	}

	return identity;
}

Ref<DenseMatrix> DenseMatrix::zero(int size) {
	if (unlikely(size < 0)) {
		ERR_PRINT_ED("ERROR: cannot create matrix of negative size.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	Ref<DenseMatrix> zero;
	zero.instantiate();

	zero->row_number = size;
	zero->col_number = size;
	zero->values.resize(size * size, 0.0);

	return zero;
}

void DenseMatrix::set_dimensions(int rows, int cols) {
	if (unlikely(rows < 0 || cols < 0)) {
		ERR_PRINT_ED("ERROR: cannot set dimensions to negative.");
		return;
	}
	std::vector<double> new_values;
	new_values.resize(rows * cols, 0.0);
	
	int row_copy_count = 0;
	if (rows < row_number) {
		row_copy_count = rows;
	} else {
		row_copy_count = row_number;
	}

	// copy data, to keep data after truncation
	if (cols < col_number) {
		for (int copy_row = 0; copy_row < row_copy_count; copy_row++) {
			memcpy(&new_values[cols * copy_row], &values[col_number * copy_row], sizeof(double) * cols);
		}
	} else { // col_number <= cols, therefore, expand instead of truncate
		for (int copy_row = 0; copy_row < row_copy_count; copy_row++) {
			memcpy(&new_values[cols * copy_row], &values[col_number * copy_row], sizeof(double) * col_number);
		}
	}
	values = new_values;
	row_number = rows;
	col_number = cols;
}

void DenseMatrix::set_element(int row, int col, double value) {
	if (unlikely(row >= row_number || col >= col_number || row < 0 || col < 0)) {
        ERR_PRINT_ED("ERROR: cannot set cell outside dimension range.");
		return;
    }
	values[col_number * row + col] = value;
}

double DenseMatrix::get_element(int row, int col) const {
	if (unlikely(row >= values.size() || col >= col_number || row < 0 || col < 0)) {
        ERR_PRINT_ED("ERROR: cannot get cell outside dimension range.");
		return 0.0;
    }
	return values[col_number * row + col];
}

Ref<DenseMatrix> DenseMatrix::transposed() const {
	Ref<DenseMatrix> transposed;
	transposed.instantiate();
	transposed->row_number = col_number;
	transposed->col_number = row_number;
	transposed->values.resize(row_number * col_number);
	// do not use set_dimensions, so that value initialisation does not occur, and is skipped.

	for (int new_row = 0; new_row < col_number; new_row++) {
		for (int new_col = 0; new_col < row_number; new_col++) {
			transposed->values[row_number * new_row + new_col] = values[col_number * new_col + new_row];
		}
	}

	return transposed;
}

Ref<DenseMatrix> DenseMatrix::clone() const {
	Ref<DenseMatrix> cloned;
	cloned.instantiate();

	cloned->row_number = row_number;
	cloned->col_number = col_number;
	cloned->values = values;

	return cloned;
}


void DenseMatrix::multiply_scaler_in_place(double scaler) {
	for (int row = 0; row < row_number; row++) {
		for (int col = 0; col < col_number; col++) {
			values[col_number * row + col] *= scaler;
		}
	}
}

Ref<DenseMatrix> DenseMatrix::multiply_dense(Ref<DenseMatrix> other) const {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	if (unlikely(col_number != other->row_number)) {
        ERR_PRINT_ED("ERROR: cannot multiply matrices with incompatible dimensions. Left-hand row length must equal Right-hand column length.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
    }

	Ref<DenseMatrix> result;
	result.instantiate();
	result->set_dimensions(row_number, other->col_number);

	for (int row = 0; row < row_number; row++) {
		for (int col = 0; col < other->col_number; col++) {
			double total = 0.0;

			for (int i = 0; i < col_number; i++) {
				total += values[col_number * row + i] * other->values[other->col_number * i + col];
			}

			result->values[result->col_number * row + col] = total;
		}
	}

	return result;
}


Ref<DenseMatrix> DenseMatrix::multiply_sparse(Ref<SparseMatrix> other) const {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	if (unlikely(col_number != other->values.size())) {
        ERR_PRINT_ED("ERROR: cannot multiply matrices with incompatible dimensions. Left-hand row length must equal Right-hand column length.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
    }

	Ref<SparseMatrix> transposed = other->transposed();

	Ref<DenseMatrix> result;
	result.instantiate();
	result->set_dimensions(row_number, other->column_number);

	for (int row = 0; row < row_number; row++) {
		for (int col = 0; col < other->column_number; col++) {
			double total = 0.0;
			
			for (auto& it : transposed->values[col]) {
				total += it.val * values[col_number * row + it.column];
			}

			result->values[result->col_number * row + col] = total;
		}
	}

	return result;
}

Ref<DenseMatrix> DenseMatrix::add_sparse(Ref<SparseMatrix> other) const {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	if (unlikely(row_number != other->values.size() || col_number != other->column_number)) {
        ERR_PRINT_ED("ERROR: cannot add matrices with dimensions that are not equal.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
    }

	Ref<DenseMatrix> result = clone(); // clone values (many of them won't be touched. Only a few will be editted by the sparse matrix)
	for (int row = 0; row < row_number; row++) {
		for (auto& col : other->values[row]) {
			result->values[col_number * row + col.column] += col.val;
		}
	}

	return result;
}

Ref<DenseMatrix> DenseMatrix::subtract_sparse(Ref<SparseMatrix> other) const {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	if (unlikely(row_number != other->values.size() || col_number != other->column_number)) {
        ERR_PRINT_ED("ERROR: cannot subtract matrices with dimensions that are not equal.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
    }

	Ref<DenseMatrix> result = clone(); // clone values (many of them won't be touched. Only a few will be editted by the sparse matrix)
	for (int row = 0; row < row_number; row++) {
		for (auto& col : other->values[row]) {
			result->values[col_number * row + col.column] -= col.val;
		}
	}

	return result;
}

Ref<DenseMatrix> DenseMatrix::add_dense(Ref<DenseMatrix> other) const {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	if (unlikely(row_number != other->row_number || col_number != other->col_number)) {
        ERR_PRINT_ED("ERROR: cannot add matrices with dimensions that are not equal.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
    }

	Ref<DenseMatrix> result;
	result.instantiate();
	result->set_dimensions(row_number, col_number);

	for (int i = 0; i < row_number * col_number; i++) {
		result->values[i] = values[i] + other->values[i];
	}

	return result;
}

Ref<DenseMatrix> DenseMatrix::subtract_dense(Ref<DenseMatrix> other) const {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	if (unlikely(row_number != other->row_number || col_number != other->col_number)) {
        ERR_PRINT_ED("ERROR: cannot subtract matrices with dimensions that are not equal.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
    }

	Ref<DenseMatrix> result;
	result.instantiate();
	result->set_dimensions(row_number, col_number);

	for (int i = 0; i < row_number * col_number; i++) {
		result->values[i] = values[i] - other->values[i];
	}

	return result;
}

void DenseMatrix::add_sparse_in_place(Ref<SparseMatrix> other) {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		return;
	}
	if (unlikely(row_number != other->values.size() || col_number != other->column_number)) {
        ERR_PRINT_ED("ERROR: cannot add matrices with dimensions that are not equal.");
		return;
    }

	for (int row = 0; row < row_number; row++) {
		for (auto& col : other->values[row]) {
			values[col_number * row + col.column] += col.val;
		}
	}
}

void DenseMatrix::subtract_sparse_in_place(Ref<SparseMatrix> other) {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		return;
	}
	if (unlikely(row_number != other->values.size() || col_number != other->column_number)) {
        ERR_PRINT_ED("ERROR: cannot subtract matrices with dimensions that are not equal.");
		return;
    }

	for (int row = 0; row < row_number; row++) {
		for (auto& col : other->values[row]) {
			values[col_number * row + col.column] -= col.val;
		}
	}
}

void DenseMatrix::add_dense_in_place(Ref<DenseMatrix> other) {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		return;
	}
	if (unlikely(row_number != other->row_number || col_number != other->col_number)) {
        ERR_PRINT_ED("ERROR: cannot add matrices with dimensions that are not equal.");
		return;
    }

	for (int i = 0; i < row_number * col_number; i++) {
		values[i] += other->values[i];
	}
}

void DenseMatrix::subtract_dense_in_place(Ref<DenseMatrix> other) {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		return;
	}
	if (unlikely(row_number != other->row_number || col_number != other->col_number)) {
        ERR_PRINT_ED("ERROR: cannot subtract matrices with dimensions that are not equal.");
		return;
    }

	for (int i = 0; i < row_number * col_number; i++) {
		values[i] -= other->values[i];
	}
}

double DenseMatrix::norm_squared() const {
	double total = 0.0;

	for (int i = 0; i < values.size(); i++) {
		total += values[i] * values[i];
	}

	return total;
}

double DenseMatrix::norm() const {
	return sqrt(norm_squared());
}

Ref<VectorN> DenseMatrix::multiply_vector(Ref<VectorN> other) const {
	if (unlikely(other.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		Ref<VectorN> null_ref;
		return null_ref;
	}
	if (unlikely(col_number != other->values.size())) {
        ERR_PRINT_ED("ERROR: cannot multiply matrices with vectors of a length not equal to the matrix column number.");
		Ref<VectorN> null_ref;
		return null_ref;
    }

	Ref<VectorN> result;
	result.instantiate();
	result->values.resize(row_number);

	for (int row = 0; row < row_number; row++) {
		double total = 0.0;

		for (int col = 0; col < col_number; col++) {
			total += values[col_number * row + col] * other->values[col];
		}

		result->values[row] = total;
	}

	return result;
}

Ref<SparseMatrix> DenseMatrix::to_sparse(double zero_threshold) const {
	Ref<SparseMatrix> result;
	result.instantiate();
	result->set_dimensions(row_number, col_number);

	for (int row = 0; row < row_number; row++) {
		for (int col = 0; col < col_number; col++) {
			double value = values[col_number * row + col];

			if (value > zero_threshold) {
				Value newvalue;
				newvalue.val = value;
				newvalue.column = col;
				result->values[row].push_back(newvalue);
			}
		}
	}

	return result;
}

Ref<DenseMatrix> DenseMatrix::solve(Ref<DenseMatrix> B) const {
	if (unlikely(B.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	if (unlikely(row_number != col_number)) {
		ERR_PRINT_ED("ERROR: cannot solve with a non-square coefficient matrix.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	if (unlikely(col_number != B->row_number)) {
		ERR_PRINT_ED("ERROR: B row number must be equal to the coefficient matrix width.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	// solves using guassian elimination - then, out of row-echelon form, and then applies back-substitution/(this step is isomorphic to back-elimination)

	Ref<DenseMatrix> result;
	result.instantiate();
	result->set_dimensions(row_number, col_number + B->col_number); // create it as an augmented matrix
	result->blit_rect_unsafe(*this, Rect2i(0, 0, row_number, row_number), Vector2i(0, 0)); // add the current matrix into the LHS
	result->blit_rect_unsafe(*B.ptr(), Rect2i(0, 0, B->row_number, B->col_number), Vector2i(0, col_number)); // add the solution matrix into the RIGHT of the augmented

	std::string return_string = "";
	auto matrix_values = result->to_packed_array();
	

	for (int row = 0; row < row_number; row++) {
		// step 1: search for the pivot element (greatest element), in rows including this, and above
		int pivot = row;
		double pivot_value = result->get_cell_unsafe(row, row);
		for (int search_row = row + 1; search_row < row_number; search_row++) {
			double search_value = result->get_cell_unsafe(search_row, row);
			if (search_value > pivot_value) {
				pivot_value = search_value;
				pivot = search_row;
			}
		}

		if (pivot != row) { // then switch the selected pivot row to the top
			result->swap_rows(row, pivot); // "after col" is inclusive
		}

		// Now that a good pivot has been selected at the top of the current working row, then start.
		if (std::abs(pivot_value) < std::numeric_limits<double>::epsilon()) { // check if the pivot equals zero. A value of atleast 1, or higher, is extremely likely. So, this epsilon factor is used.
			Ref<DenseMatrix> empty_ref; // null reference (can be compared in Godot using == null). This is the failure case.
			WARN_PRINT_ED("ERROR: this matrix is non-invertible/degenerate. Returned null.");
			return empty_ref;
		}
		// Divide by the pivot value to reduce the row.
		result->get_cell_unsafe(row, row) = 1.0;
		result->multiply_row_after_col(row, 1.0/pivot_value, row + 1);

		// Now with the reduced row, eliminate on the below rows
		for (int eliminate_row = row + 1; eliminate_row < row_number; eliminate_row++) {
			auto& root_cell = result->get_cell_unsafe(eliminate_row, row);
			if (std::abs(root_cell) < std::numeric_limits<double>::epsilon() * 1.0) { // if it is zero, there is no need to subtract
				continue;
			}
			double multiple = -root_cell;
			root_cell = 0.0;
			result->add_row_multiple_after_col(eliminate_row, row, multiple, row+1);
		}
		// Now the matrix, up to "row", is in row echelon form. Repeat for every row until the end.
	}

	// Now that it is in row echelon form, reduce backwards from the bottom.
	for (int row = row_number - 1; row >= 0; row--) {
		for (int eliminate_row = 0; eliminate_row < row; eliminate_row++) {
			double multiple = -result->get_cell_unsafe(eliminate_row, row);
			if (std::abs(multiple) < std::numeric_limits<double>::epsilon() * 1.0) { // if it is zero, there is no effect
				continue;
			}
			result->add_row_multiple_after_col(eliminate_row, row, multiple, col_number); // note: do not need to record the reductions on the LSH matrix.
		}
	}

	// then, crop result to the RHS of the augmented part
	Ref<DenseMatrix> result_cropped;
	result_cropped.instantiate();
	result_cropped->row_number = row_number;
	result_cropped->col_number = B->col_number;
	result_cropped->values.resize(row_number * B->col_number);
	result_cropped->blit_rect_unsafe_ref(result, Rect2i(0, col_number, row_number, B->col_number), Vector2i(0, 0));

	return result_cropped;
}

Ref<DenseMatrix> DenseMatrix::inverse() const {
	if (unlikely(row_number != col_number)) {
		ERR_PRINT_ED("ERROR: cannot invert non-square matrices.");
		Ref<DenseMatrix> null_ref;
		return null_ref;
	}
	
	return solve(identity(row_number));
}

Ref<VectorN> DenseMatrix::solve_iterative_cg(Ref<VectorN> B, Ref<VectorN> initial_guess, int max_iterations) const {
	if (unlikely(B.is_null())) {
		ERR_PRINT_ED("ERROR: parameter is null.");
		Ref<VectorN> null_ref;
		return null_ref;
	}
	if (unlikely(row_number != col_number)) {
		ERR_PRINT_ED("ERROR: matrix must be square.");
		Ref<VectorN> null_ref;
		return null_ref;
	}

	if (initial_guess.is_null()) {
		initial_guess = VectorN::filled(0.0, B->values.size());
	}

	if (max_iterations > row_number) max_iterations = row_number; // the conjugate method is theoretically perfect after N steps
	if (max_iterations < 0) max_iterations = row_number;

	// algorithm was stolen from Wikipedia: (https://en.wikipedia.org/wiki/Conjugate_gradient_method)

	Ref<VectorN> X = initial_guess->clone();
	Ref<VectorN> R = B->subtract(multiply_vector(X));
	Ref<VectorN> P = R->clone();

	for (int i = 0; i < max_iterations; i++) {
		if (R->is_approximately_zero(std::numeric_limits<double>::epsilon())) {
			break;
		}

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