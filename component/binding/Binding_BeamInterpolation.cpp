/*********************************************************************
Copyright 2019, Inria, CNRS, University of Lille

This file is part of BeamAdapter

BeamAdapter is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BeamAdapter is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with sofaqtquick. If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/
/********************************************************************
 Contributors:
    - damien.marchal@univ-lille.fr
********************************************************************/
#include <sofa/defaulttype/RigidTypes.h>
using sofa::defaulttype::Rigid3Types;

#include <SofaPython/PythonToSofa.inl>
#include <SofaPython/PythonFactory.h>
#include "../../component/BeamInterpolation.h"
using sofa::component::fem::BeamInterpolation;

#include "Binding_BeamInterpolation.h"

namespace sofa
{

typedef BeamInterpolation<Rigid3Types> BeamInterpolationRigid3;

static BeamInterpolation<Rigid3Types>* get_beaminterpolation(PyObject* self)
{
    return sofa::py::unwrap<BeamInterpolationRigid3>(self);
}


static PyObject* BeamInterpolation_getValuesAt(PyObject *self, PyObject *args)
{
    /// Extract the Sofa type from the self argument.
    BeamInterpolationRigid3* beam = get_beaminterpolation(self);

    /// Extracts the argument from python (type is a string while the other is an array)
    PyObject* arg {nullptr};
    PyObject* stype {nullptr};
    const char* ctype {"current"};

    /// Parse two arguments, the second one is optional.
    if (!PyArg_ParseTuple(args, "O|O", &arg, &stype))
        return nullptr;

    PyObject *iter = PyObject_GetIter(arg);
    if (!iter)
    {
        PyErr_SetString(PyExc_TypeError, "Expecting an iterable object");
        return nullptr;
    }

    if ( stype != nullptr )
    {
        if (PyString_Check(stype))
        {
            PyErr_SetString(PyExc_TypeError, "Expecting a string as second parameter");
            return nullptr;
        }
        ctype = PyString_AsString(stype);
    }

    sofa::component::fem::_beaminterpolation_::CoordType type;
    if(!strcmp(ctype, "current"))
        type = sofa::component::fem::_beaminterpolation_::CoordType::CurrentPosition;
    else if(!strcmp(ctype, "free"))
        type = sofa::component::fem::_beaminterpolation_::CoordType::FreePosition;
    else if(!strcmp(ctype, "rest"))
        type = sofa::component::fem::_beaminterpolation_::CoordType::RestPosition;
    else
    {
        PyErr_SetString(PyExc_RuntimeError, "Invalid 'field' value. Expecting 'current, free or rest'");
        return nullptr;
    }

    BeamInterpolationRigid3::LocalCoordinates coords;
    while(PyObject* item = PyIter_Next(iter))
    {
        BeamInterpolationRigid3::LocalCoordinate coord {PyLong_AsUnsignedLong(PyTuple_GetItem(item, 0)),
                                                        PyFloat_AsDouble(PyTuple_GetItem(item, 1))} ;
        if(coord.indice>=beam->getNumBeams())
        {
            std::stringstream s;
            s << "Out of bound indice. '" << coord.indice << "' should be < " << beam->getNumBeams();
            PyErr_SetString(PyExc_RuntimeError, s.str().c_str());
            Py_DECREF(item);
            Py_DECREF(iter);
            return nullptr;
        }

        coords.push_back(coord);
        Py_DECREF(item);
    }
    Py_DECREF(iter);

    BeamInterpolationRigid3::InterpolatedValues values;
    beam->getValuesAt(type, coords, values);

    PyObject *list = PyList_New(coords.size());
    unsigned int idest=0;
    for(const auto& coord : coords)
    {
        auto& pos = values[idest].position;
        auto& vel = values[idest].velocity;

        PyObject *v = PyList_New(coords.size());
        for(size_t j=0;j<pos.size();++j)
            PyList_SetItem(v, j, PyFloat_FromDouble(pos[j]));

        PyList_SetItem(list, idest, v);
        idest++;
    }

    return list;
}

static PyObject* BeamInterpolation_getLength(PyObject * self, PyObject *args)
{
    BeamInterpolationRigid3* beam = get_beaminterpolation(self);
    if(!beam){
        PyErr_SetString(PyExc_TypeError, "Expecting a Sofa.BeamInterpolation");
        return nullptr;
    }

    size_t indice;
    if(!PyArg_ParseTuple(args, "i", &indice))
    {
        return nullptr;
    }

    if(indice>=beam->getNumBeams())
    {
        PyErr_SetString(PyExc_RuntimeError, "The 'index' parameter must be smaller than the number of beams.");
        return nullptr;
    }

    return PyFloat_FromDouble(beam->getLength(indice));
}

static PyObject* BeamInterpolation_getRestTotalLength(PyObject * self, PyObject *args)
{

    BeamInterpolationRigid3* beam = get_beaminterpolation(self);
    if(!beam){
        PyErr_SetString(PyExc_TypeError, "Expecting a Sofa.BeamInterpolation");
        return nullptr;
    }

    return PyFloat_FromDouble(beam->getRestTotalLength());
}

static PyObject* BeamInterpolation_getNumBeams(PyObject * self, PyObject *args)
{
    BeamInterpolationRigid3* beam = get_beaminterpolation(self);
    if(!beam){
        PyErr_SetString(PyExc_TypeError, "Expecting a Sofa.BeamInterpolation");
        return nullptr;
    }

    return PyLong_FromLong(beam->getNumBeams());
}


SP_CLASS_METHODS_BEGIN(BeamInterpolation)
SP_CLASS_METHOD(BeamInterpolation, getValuesAt)
SP_CLASS_METHOD_DOC(BeamInterpolation, getLength,
                    "Returns the length of a beam element at a given index.\n"
                    "\n"
                    "Example: c.getLength(2)")
SP_CLASS_METHOD_DOC(BeamInterpolation, getRestTotalLength,
                "Returns the total length of the beam in its rest position.\n"
                "\n"
                "Example: c.getRestTotalLength()")
SP_CLASS_METHOD_DOC(BeamInterpolation, getNumBeams,
                    "Returns the number of beam elements.\n"
                    "\n"
                    "Example: c.getNumBeams()")
SP_CLASS_METHODS_END

SP_CLASS_TYPE_SPTR(BeamInterpolation,
                   BeamInterpolationRigid3, BaseObject);

void addBeamInterpolation(PyObject *module)
{
    SP_ADD_CLASS(module, BeamInterpolation);
    PythonFactory::add<BeamInterpolationRigid3>( &SP_SOFAPYTYPEOBJECT(BeamInterpolation) );
}

}
