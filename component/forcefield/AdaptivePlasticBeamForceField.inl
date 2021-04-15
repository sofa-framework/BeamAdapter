/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
 //
 // C++ Implementation : AdaptivePlasticBeamForceField
 //
 // Description: Implementation of plasticity over the AdaptiveBeamForceFieldAndMass
 // force field interface.
 //
 //
 // Author: Camille Krewcun, INRIA
 //
 // Copyright: See COPYING file that comes with this distribution
 //
 //
#pragma once

#include "AdaptivePlasticBeamForceField.h"

namespace sofa::plugin::beamadapter::component::forcefield
{

namespace _adaptiveplasticbeamforcefield_
{

template <class DataTypes>
AdaptivePlasticBeamForceField<DataTypes>::AdaptivePlasticBeamForceField()
{
}

template <class DataTypes>
AdaptivePlasticBeamForceField<DataTypes>::~AdaptivePlasticBeamForceField()
{
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::init()
{

}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::reinit()
{

}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::draw(const VisualParams* vparams)
{

}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::addForce(const MechanicalParams* mparams,
                                                        DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v)
{

}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::addDForce(const MechanicalParams* mparams,
                                                         DataVecDeriv& datadF, const DataVecDeriv& datadX)
{

}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::addKToMatrix(const MechanicalParams* mparams,
                                                            const MultiMatrixAccessor* matrix)
{

}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::computeStiffness(int beam, BeamLocalMatrices& beamLocalMatrices)
{

}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::applyStiffnessLarge(VecDeriv& df, const VecDeriv& dx, int beam, Index nd0Id, Index nd1Id, const double& factor)
{

}

} // namespace _adaptiveplasticbeamforcefield_
} // namespace sofa::plugin::beamadapter::component::forcefield