/* godot-cpp integration testing project.
 *
 * This is free and unencumbered software released into the public domain.
 */

#ifndef EXAMPLE_CLASS_H
#define EXAMPLE_CLASS_H

// We don't need windows.h in this example plugin but many others do, and it can
// lead to annoying situations due to the ton of macros it defines.
// So we include it and make sure CI warns us if we use something that conflicts
// with a Windows define.
#ifdef WIN32
#include <windows.h>
#endif

#include <godot_cpp/classes/control.hpp>
#include <godot_cpp/classes/global_constants.hpp>
#include <godot_cpp/classes/image.hpp>
#include <godot_cpp/classes/input_event_key.hpp>
#include <godot_cpp/classes/viewport.hpp>

#include <godot_cpp/core/binder_common.hpp>

using namespace godot;

class ExampleRef : public RefCounted {
	GDCLASS(ExampleRef, RefCounted);

private:
	static int instance_count;
	static int last_id;

	int id;

protected:
	static void _bind_methods();

public:
	ExampleRef();
	~ExampleRef();

	void set_id(int p_id);
	int get_id() const;
};

class ExampleMin : public Control {
	GDCLASS(ExampleMin, Control);

protected:
	static void _bind_methods(){};
};

class Example : public Control {
	GDCLASS(Example, Control);

protected:
	static void _bind_methods();

	void _notification(int p_what);
	bool _set(const StringName &p_name, const Variant &p_value);
	bool _get(const StringName &p_name, Variant &r_ret) const;
	void _get_property_list(List<PropertyInfo> *p_list) const;
	bool _property_can_revert(const StringName &p_name) const;
	bool _property_get_revert(const StringName &p_name, Variant &r_property) const;

	String _to_string() const;

private:
	Vector2 custom_position;
	Vector3 property_from_list;
	Vector2 dprop[3];
	int last_rpc_arg = 0;

public:
	// Constants.
	enum Constants {
		FIRST,
		ANSWER_TO_EVERYTHING = 42,
	};

	enum Flags {
		FLAG_ONE = 1,
		FLAG_TWO = 2,
	};

	enum {
		CONSTANT_WITHOUT_ENUM = 314,
	};

	Example();
	~Example();

	// Functions.
	void simple_func();
	void simple_const_func() const;
	int custom_ref_func(Ref<ExampleRef> p_ref);
	int custom_const_ref_func(const Ref<ExampleRef> &p_ref);
	String image_ref_func(Ref<Image> p_image);
	String image_const_ref_func(const Ref<Image> &p_image);
	String return_something(const String &base);
	Viewport *return_something_const() const;
	Ref<ExampleRef> return_ref() const;
	Ref<ExampleRef> return_empty_ref() const;
	ExampleRef *return_extended_ref() const;
	Ref<ExampleRef> extended_ref_checks(Ref<ExampleRef> p_ref) const;
	Variant varargs_func(const Variant **args, GDExtensionInt arg_count, GDExtensionCallError &error);
	int varargs_func_nv(const Variant **args, GDExtensionInt arg_count, GDExtensionCallError &error);
	void varargs_func_void(const Variant **args, GDExtensionInt arg_count, GDExtensionCallError &error);
	void emit_custom_signal(const String &name, int value);
	int def_args(int p_a = 100, int p_b = 200);

	Array test_array() const;
	int test_tarray_arg(const TypedArray<int64_t> &p_array);
	TypedArray<Vector2> test_tarray() const;
	Dictionary test_dictionary() const;
	Example *test_node_argument(Example *p_node) const;
	String test_string_ops() const;
	String test_str_utility() const;
	bool test_string_is_fourty_two(const String &p_str) const;
	int test_vector_ops() const;

	BitField<Flags> test_bitfield(BitField<Flags> flags);

	// RPC
	void test_rpc(int p_value);
	void test_send_rpc(int p_value);
	int return_last_rpc_arg();

	// Property.
	void set_custom_position(const Vector2 &pos);
	Vector2 get_custom_position() const;
	Vector4 get_v4() const;

	// Static method.
	static int test_static(int p_a, int p_b);
	static void test_static2();

	// Virtual function override (no need to bind manually).
	virtual bool _has_point(const Vector2 &point) const override;
	virtual void _input(const Ref<InputEvent> &event) override;
};

VARIANT_ENUM_CAST(Example::Constants);
VARIANT_BITFIELD_CAST(Example::Flags);

enum EnumWithoutClass {
	OUTSIDE_OF_CLASS = 512
};
VARIANT_ENUM_CAST(EnumWithoutClass);

class ExampleVirtual : public Object {
	GDCLASS(ExampleVirtual, Object);

protected:
	static void _bind_methods() {}
};

class ExampleAbstract : public Object {
	GDCLASS(ExampleAbstract, Object);

protected:
	static void _bind_methods() {}
};

#endif // EXAMPLE_CLASS_H
