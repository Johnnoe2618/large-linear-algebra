#ifndef VECTOR_N_H
#define VECTOR_N_H

#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include <godot_cpp/classes/ref.hpp>
#include <vector>

namespace godot {

class VectorN : public RefCounted {
    GDCLASS(VectorN, RefCounted)

private:
    
protected:
    static void _bind_methods();

	std::vector<double> values;
    bool is_approximately_zero(double total_threshold) const;
public:
    friend class DenseMatrix;
	friend class SparseMatrix;

    VectorN() {};
    ~VectorN() {};

    static Ref<VectorN> filled(double value, int size);
    static Ref<VectorN> from_packed_array(const PackedFloat64Array& from);

    Ref<VectorN> clone() const;
    void resize(int size);
    int get_size() const { return values.size(); };

    void set_element(int idx, double value);
    double get_element(int idx) const;

    double dot(Ref<VectorN> other) const;
    void normalise();

    Ref<VectorN> add(Ref<VectorN> other) const;
    Ref<VectorN> subtract(Ref<VectorN> other) const;
    void add_in_place(Ref<VectorN> other);
    void subtract_in_place(Ref<VectorN> other);
    void add_multiple_in_place(Ref<VectorN> other, double multiple);
    void multiply_scalar_in_place(double scalar);

    double length_squared() const;
    double length() const;
    Ref<DenseMatrix> column_vector() const;
    Ref<DenseMatrix> row_vector() const;

    PackedFloat64Array to_packed_array() const;
};

}

#endif