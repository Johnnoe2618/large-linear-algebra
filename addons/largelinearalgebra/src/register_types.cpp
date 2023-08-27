#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "VectorN.h"

#include <gdextension_interface.h>
#include <godot_cpp/core/defs.hpp>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/godot.hpp>

using namespace godot;

void initialize_largelinearalgebra_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) return;

    ClassDB::register_class<SparseMatrixCoordinate>();
    ClassDB::register_class<VectorN>();
    ClassDB::register_class<SparseMatrix>();
    ClassDB::register_class<DenseMatrix>();
}

void uninitialize_largelinearalgebra_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) return;

}

extern "C" {
// Initialization.
GDExtensionBool GDE_EXPORT largelinearalgebra_library_init(GDExtensionInterfaceGetProcAddress p_get_proc_address, const GDExtensionClassLibraryPtr p_library, GDExtensionInitialization *r_initialization) {
    godot::GDExtensionBinding::InitObject init_obj(p_get_proc_address, p_library, r_initialization);

    init_obj.register_initializer(initialize_largelinearalgebra_module);
    init_obj.register_terminator(uninitialize_largelinearalgebra_module);
    init_obj.set_minimum_library_initialization_level(MODULE_INITIALIZATION_LEVEL_CORE); // maybe at MODULE_INITIALIZATION_LEVEL_CORE, MODULE_INITIALIZATION_LEVEL_SCENE, or other

    return init_obj.init();
}
}