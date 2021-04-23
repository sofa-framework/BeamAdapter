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

using sofa::core::objectmodel::BaseContext;

template <class DataTypes>
AdaptivePlasticBeamForceField<DataTypes>::AdaptivePlasticBeamForceField()
{
}

template <class DataTypes>
AdaptivePlasticBeamForceField<DataTypes>::~AdaptivePlasticBeamForceField()
{
    //TO DO : should m_gaussPoints and m_integrationIntervals be deallocated explicitly ?
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::init()
{
    if (!l_interpolation)
        l_interpolation.set(dynamic_cast<BaseContext*>(this->getContext())->get<BPInterpolation>(BaseContext::Local));

    if (!l_interpolation)
        msg_error() << "No Beam Interpolation found, the component can not work.";

    ForceField<DataTypes>::init();

    const vector<BeamGeometry> beamGeometryParams = l_interpolation->getBeamGeometryParameters();
    unsigned int numBeams = l_interpolation->getNumBeams();

    m_gaussPoints.clear();
    m_gaussPoints.resize(numBeams);
    m_integrationIntervals.clear();
    m_integrationIntervals.resize(numBeams);

    for (unsigned int b = 0; b < numBeams; b++)
    {
        initialiseInterval(b, m_integrationIntervals, beamGeometryParams[b]);
        initialiseGaussPoints(b, m_gaussPoints, m_integrationIntervals[b]);
    }

}


template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::initialiseGaussPoints(int beam, vector<beamGaussPoints>& gaussPoints,
                                                                     const Interval3& integrationInterval)
{
    //Gaussian nodes coordinates and weights for a 1D integration on [-1,1]
    const double sqrt3_5 = helper::rsqrt( 3/5 );
    Vec3 canonical3NodesCoordinates = {-sqrt3_5, 0, sqrt3_5};
    Vec3 canonical3NodesWeights = {5/9, 8/9, 5/9};

    //Compute actual Gauss points coordinates and weights, with a 3D integration
    //NB: 3 loops because integration is in 3D, 3 iterations per loop because it's a 3 point integration
    unsigned int gaussPointIt = 0;
    for (unsigned int i = 0; i < 3; i++)
    {
        double x = canonical3NodesCoordinates[i];
        double w1 = canonical3NodesWeights[i];
        // Changing first coordinate and weight to adapt to the integration interval
        double a1 = integrationInterval.m_a1;
        double b1 = integrationInterval.m_b1;
        double xChanged = changeCoordinate(x, a1, b1);
        double w1Changed = changeWeight(w1, a1, b1);

        for (unsigned int j = 0; j < 3; j++)
        {
            double y = canonical3NodesCoordinates[j];
            double w2 = canonical3NodesWeights[j];
            // Changing second coordinate and weight to adapt to the integration interval
            double a2 = integrationInterval.m_a2;
            double b2 = integrationInterval.m_b2;
            double yChanged = changeCoordinate(y, a2, b2);
            double w2Changed = changeWeight(w2, a2, b2);

            for (unsigned int k = 0; k < 3; k++)
            {
                double z = canonical3NodesCoordinates[k];
                double w3 = canonical3NodesWeights[k];
                // Changing third coordinate and weight to adapt to the integration interval
                double a3 = integrationInterval.m_a3;
                double b3 = integrationInterval.m_b3;
                double zChanged = changeCoordinate(z, a3, b3);
                double w3Changed = changeWeight(w3, a3, b3);

                GaussPoint3 newGaussPoint = GaussPoint3(xChanged, yChanged, zChanged, w1Changed, w2Changed, w3Changed);
                gaussPoints[beam][gaussPointIt] = newGaussPoint;
                gaussPointIt++;
            }
        }
    }
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::initialiseInterval(int beam, vector<Interval3>& integrationIntervals,
                                                                  const BeamGeometry& beamGeometry)
{
    if (beamGeometry.sectionShape == "rectangular")
    {
        Real L = beamGeometry._L;
        Real Ly = beamGeometry._Ly;
        Real Lz = beamGeometry._Lz;

        // Integration interval definition for a local frame at the centre of the beam
        integrationIntervals.push_back(Interval3(-L/2, L/2, -Ly/2, Ly/2, -Lz/2, -Lz/2));
    }
    else if (beamGeometry.sectionShape == "elliptic")
    {
        // TO DO: implement quadrature method for elliptic cross section
        msg_error() << "Quadrature method for " << beamGeometry.sectionShape
            << " shape cross section has not been implemented yet. Methods for rectangular cross sections are available";
    }
    else if (beamGeometry.sectionShape == "circular")
    {
        //TO DO: implement quadrature method for a disc and a hollow-disc cross section
        msg_error() << "Quadrature method for " << beamGeometry.sectionShape
            << " shape cross section has not been implemented yet. Methods for rectangular cross sections are available";
    }
    else
    {
        msg_error() << "Quadrature method for " << beamGeometry.sectionShape
            << " shape cross section has not been implemented yet. Methods for rectangular cross sections are available";
    }
}


template <class DataTypes>
template <typename LambdaType>
void AdaptivePlasticBeamForceField<DataTypes>::integrate(const beamGaussPoints& gaussPoints, LambdaType integrationFun)
{
    //Apply a generic (lambda) integration function to each Gauss point of a beam element
    for (unsigned int gp = 0; gp < gaussPoints.size(); gp++)
    {
        integrationFun(gaussPoints[gp].m_coordinates, gaussPoints[gp].m_weights);
    }
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