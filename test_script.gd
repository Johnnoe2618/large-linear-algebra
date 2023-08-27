extends Node2D

func _ready():
	## Creating the Matrices
	assert(Array(DenseMatrix.identity(3).to_packed_array()).hash() == [1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0].hash())
	assert(Array(DenseMatrix.zero(2).to_packed_array()).hash() == [0.0,0.0,0.0,0.0].hash())
	assert(Array(DenseMatrix.from_packed_array([5.4, 3.7, 0.0, 1.2, 9.8, -8.5], 2, 3).to_packed_array()).hash() == [5.4, 3.7, 0.0, 1.2, 9.8, -8.5].hash())
	assert(DenseMatrix.from_packed_array([5.4, 3.7, 0.0, 1.2, 9.8, -8.5], 2, 3).get_dimensions() == Vector2i(2,3))
	
	var original = DenseMatrix.zero(2)
	original.set_element(0,1, 5.8)
	var new_clone = original.clone()
	new_clone.set_element(0,0, 7.8)
	
	assert(Array(new_clone.to_packed_array()).hash() == [7.8, 5.8, 0.0, 0.0].hash())
	assert(Array(original.to_packed_array()).hash() == [0.0, 5.8, 0.0, 0.0].hash())
	assert(Array(new_clone.to_packed_array()).hash() != [0.0, 5.8, 0.0, 0.0].hash())
	
	new_clone.set_element(1, 0, 0.1)
	assert(Array(new_clone.multiply_vector(VectorN.from_packed_array([1.0, 1.0])).to_packed_array()).hash() == [13.6, 0.1].hash())
	assert(Array(new_clone.to_sparse(0.0).multiply_vector(VectorN.from_packed_array([1.0, 1.0])).to_packed_array()).hash() == [13.6, 0.1].hash())
	
	var sparse = new_clone.to_sparse(0.0)
	
	assert(Array(sparse.add_dense(new_clone).to_packed_array()).hash() == Array(new_clone.add_sparse(sparse).to_packed_array()).hash())
	
	assert(Array(sparse.add_dense(new_clone).to_packed_array()).hash() == [7.8*2, 5.8*2, 0.1*2, 0.0].hash())
	assert(Array(sparse.subtract_dense(new_clone).to_packed_array()).hash() == [0.0, 0.0, 0.0, 0.0].hash())
	assert(Array(sparse.add_sparse(sparse).to_dense().to_packed_array()).hash() == [7.8*2, 5.8*2, 0.1*2, 0.0].hash())
	assert(Array(sparse.subtract_sparse(sparse).to_dense().to_packed_array()).hash() == [0.0, 0.0, 0.0, 0.0].hash())
	
	assert(Array(new_clone.add_dense(new_clone).to_packed_array()).hash() == [7.8*2, 5.8*2, 0.1*2, 0.0].hash())
	assert(Array(new_clone.add_sparse(sparse).to_packed_array()).hash() == [7.8*2, 5.8*2, 0.1*2, 0.0].hash())
	assert(Array(new_clone.subtract_dense(new_clone).to_packed_array()).hash() == [0.0, 0.0, 0.0, 0.0].hash())
	assert(Array(new_clone.subtract_sparse(sparse).to_packed_array()).hash() == [0.0, 0.0, 0.0, 0.0].hash())
	
	var in_place_test = new_clone.clone()
	in_place_test.add_dense_in_place(new_clone)
	assert(Array(in_place_test.to_packed_array()).hash() == [7.8*2, 5.8*2, 0.1*2, 0.0].hash())
	
	in_place_test = new_clone.clone()
	in_place_test.add_sparse_in_place(sparse)
	assert(Array(in_place_test.to_packed_array()).hash() == [7.8*2, 5.8*2, 0.1*2, 0.0].hash())
	
	in_place_test = new_clone.clone()
	in_place_test.subtract_dense_in_place(new_clone)
	assert(Array(in_place_test.to_packed_array()).hash() == [0.0, 0.0, 0.0, 0.0].hash())
	
	in_place_test = new_clone.clone()
	in_place_test.subtract_sparse_in_place(sparse)
	assert(Array(in_place_test.to_packed_array()).hash() == [0.0, 0.0, 0.0, 0.0].hash())
	
	assert(Array(new_clone.solve(VectorN.from_packed_array([1.0, 1.0]).column_vector()).to_packed_array()).hash() == [5.8/0.58, -7.7/0.58].hash())
	print((DenseMatrix.from_packed_array([3.0, 1.0, 1.0, 3.0], 2, 2)).solve_iterative_cg(VectorN.from_packed_array([3.0, 0.25]), VectorN.from_packed_array([3.093275, -2.28125]), 2).to_packed_array(), " should be approximately ", [1.09375, -0.28125])
	print((DenseMatrix.from_packed_array([3.0, 1.0, 1.0, 3.0], 2, 2).to_sparse(0.0)).solve_iterative_cg(VectorN.from_packed_array([3.0, 0.25]), VectorN.from_packed_array([3.093275, -2.28125]), 2).to_packed_array(), " should be approximately ", [1.09375, -0.28125])
	
	assert(DenseMatrix.from_packed_array([3.0, 1.0, 2.0, -1.0], 2, 2).norm_squared() == 15.0)
	assert(DenseMatrix.from_packed_array([3.0, 1.0, 2.0, -1.0], 2, 2).norm() == sqrt(15.0))
	
	assert(Array(SparseMatrix.identity(2).to_dense().to_packed_array()).hash() == [1.0, 0.0, 0.0, 1.0].hash())
	assert(Array(SparseMatrix.zero(3).to_dense().to_packed_array()).hash() == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0].hash())
	SparseMatrix.new().reserve_per_row(5) # check if it crashes
	SparseMatrix.identity(3).reserve_per_row(5) # check if it crashes
	
	var new_sparse_matrix = SparseMatrix.new()
	new_sparse_matrix.set_dimensions(2,3)
	assert(new_sparse_matrix.get_element(1,1) == 0.0)
	new_sparse_matrix.set_element(0,0,0.5)
	new_sparse_matrix.set_element(1,2,2.3)
	
	assert(Array(new_sparse_matrix.to_dense().to_packed_array()).hash() == [0.5, 0.0, 0.0, 0.0 ,0.0, 2.3].hash())
	
	var clone = new_sparse_matrix.clone()
	clone.set_element(0,1, 1.3)
	assert(Array(new_sparse_matrix.to_dense().to_packed_array()).hash() == [0.5, 0.0, 0.0, 0.0 ,0.0, 2.3].hash())
	assert(Array(clone.to_dense().to_packed_array()).hash() == [0.5, 1.3, 0.0, 0.0 ,0.0, 2.3].hash())
	
	assert(new_sparse_matrix.get_element(1,2) == 2.3)
	assert(new_sparse_matrix.get_element(0,0) == 0.5)
	assert(clone.get_element(0,1) == 1.3)
	assert(clone.get_element(1,1) == 0.0)
	
	assert(Array(clone.transposed().to_dense().to_packed_array()).hash() == [0.5, 0.0, 1.3, 0.0, 0.0, 2.3].hash())
	assert(clone.to_coordinate_array()[2].value == 2.3)
	assert(clone.to_coordinate_array()[2].row == 1)
	assert(clone.to_coordinate_array()[2].column == 2)
	
	assert(Array(clone.multiply_vector(VectorN.from_packed_array([1.0, 1.0, 1.0])).to_packed_array()).hash() == [1.8, 2.3].hash())
	assert(Array(clone.to_dense().multiply_vector(VectorN.from_packed_array([1.0, 1.0, 1.0])).to_packed_array()).hash() == [1.8, 2.3].hash())
	assert(Array(clone.multiply_dense(VectorN.from_packed_array([1.0, 1.0, 1.0]).column_vector()).to_packed_array()).hash() == [1.8, 2.3].hash())
	assert(Array(clone.multiply_vector(VectorN.from_packed_array([1.0, 1.0, 1.0])).to_packed_array()).hash() == [1.8, 2.3].hash())
	assert(Array(clone.to_dense().multiply_dense(VectorN.from_packed_array([1.0, 1.0, 1.0]).column_vector()).to_packed_array()).hash() == [1.8, 2.3].hash())
	assert(Array(clone.multiply_sparse(VectorN.from_packed_array([1.0, 1.0, 1.0]).column_vector().to_sparse(0.0)).to_dense().to_packed_array()).hash() == [1.8, 2.3].hash())
	assert(Array(clone.to_dense().multiply_sparse(VectorN.from_packed_array([1.0, 1.0, 1.0]).column_vector().to_sparse(0.0)).to_packed_array()).hash() == [1.8, 2.3].hash())
	assert(Array(clone.multiply_sparse_to_dense(VectorN.from_packed_array([1.0, 1.0, 1.0]).column_vector().to_sparse(0.0)).to_packed_array()).hash() == [1.8, 2.3].hash())
	assert(Array(clone.multiply_sparse_to_dense(VectorN.from_packed_array([1.0, 0.0, 1.0]).column_vector().to_sparse(0.0)).to_packed_array()).hash() == [0.5, 2.3].hash())
	
	var new_matrix_again =  SparseMatrix.new()
	new_matrix_again.set_dimensions(2,3)
	new_matrix_again.set_element(0,0,0.5)
	new_matrix_again.set_element(1,2,2.5)
	new_matrix_again.multiply_scaler_in_place(3)
	assert(Array(new_matrix_again.to_dense().to_packed_array()).hash() == [1.5, 0.0, 0.0, 0.0, 0.0, 7.5].hash())
	
	
