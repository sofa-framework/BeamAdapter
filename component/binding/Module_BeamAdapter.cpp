#include <SofaPython/PythonMacros.h>
#include <SofaPython/PythonFactory.h>
#include "Binding_BeamInterpolation.h"
namespace sofa
{

static PyObject * BeamAdapter_toto(PyObject* /*self*/, PyObject *args)
{
    std::cout << "Hello world" << std::endl;
    Py_RETURN_NONE;
}

/// Methods of the module
SP_MODULE_METHODS_BEGIN(BeamAdapter)
SP_MODULE_METHOD(BeamAdapter, toto)
SP_MODULE_METHODS_END

/// shortcut for SP_ADD_CLASS with fixed sofa module
/// #define SP_ADD_CLASS_IN_SOFAMODULE(C) SP_ADD_CLASS( PythonFactory::s_sofaPythonModule, C )

/// This is the macro to call to bind a new type inherited from sofa::core::objectmodel::Base
/// #define SP_ADD_CLASS_IN_FACTORY( PYTHONNAME, CPPCLASSNAME ) {\
///    SP_ADD_CLASS_IN_SOFAMODULE( PYTHONNAME )  \
///    PythonFactory::add<CPPCLASSNAME>( &SP_SOFAPYTYPEOBJECT(PYTHONNAME) ); \
///    }


void initBinding()
{
    static std::string docstring=R"()";

    if( PythonFactory::s_sofaPythonModule ) // add the module only if the Sofa module exists (SofaPython is loaded)
    {
        simulation::PythonEnvironment::gil lock(__func__);
        static PyObject *s_beamAdapterPythonModule = SP_INIT_MODULE(BeamAdapter,
                                                                     docstring.c_str());

        addBeamInterpolation(s_beamAdapterPythonModule);
    }
}

}
