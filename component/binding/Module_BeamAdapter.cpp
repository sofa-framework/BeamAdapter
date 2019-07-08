#include <SofaPython/PythonMacros.h>
#include <SofaPython/PythonFactory.h>
#include "Binding_BeamInterpolation.h"
namespace sofa
{

static PyObject * BeamAdapter_getVersion(PyObject* /*self*/, PyObject *args)
{
    PyObject* tuple=PyTuple_New(2);
    PyTuple_SetItem(tuple, 0, PyLong_FromLong(1));
    PyTuple_SetItem(tuple, 1, PyLong_FromLong(0));
    return tuple;
}

/// Methods of the module
SP_MODULE_METHODS_BEGIN(BeamAdapter)
SP_MODULE_METHOD(BeamAdapter, getVersion)
SP_MODULE_METHODS_END

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
