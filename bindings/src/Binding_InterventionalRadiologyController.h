/*********************************************************************
Copyright 2019, CNRS, University of Lille, INRIA

This file is part of sofaPython3

sofaPython3 is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

sofaPython3 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with sofaqtquick. If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <sofa/core/objectmodel/Base.h>
#include <sofa/core/objectmodel/BaseData.h>
#include <SofaPython3/DataHelper.h>
#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/behavior/BaseController.h>

#include "../../component/controller/InterventionalRadiologyController.h"

using sofa::defaulttype::Rigid3Types;
using sofa::defaulttype::Vec3dTypes;

template class pybind11::class_<sofa::component::controller::InterventionalRadiologyController<Rigid3Types>,
                                sofa::core::behavior::BaseController,
                                sofa::core::sptr<sofa::component::controller::InterventionalRadiologyController<Rigid3Types> >>;


namespace sofapython3 {
    
    using sofa::component::controller::InterventionalRadiologyController;

    void moduleAddInterventionalRadiologyController(pybind11::module& m);

} // namespace sofapython3
