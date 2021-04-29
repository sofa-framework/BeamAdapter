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
        double a1 = integrationInterval.geta1();
        double b1 = integrationInterval.getb1();
        double xChanged = changeCoordinate(x, a1, b1);
        double w1Changed = changeWeight(w1, a1, b1);

        for (unsigned int j = 0; j < 3; j++)
        {
            double y = canonical3NodesCoordinates[j];
            double w2 = canonical3NodesWeights[j];
            // Changing second coordinate and weight to adapt to the integration interval
            double a2 = integrationInterval.geta2();
            double b2 = integrationInterval.getb2();
            double yChanged = changeCoordinate(y, a2, b2);
            double w2Changed = changeWeight(w2, a2, b2);

            for (unsigned int k = 0; k < 3; k++)
            {
                double z = canonical3NodesCoordinates[k];
                double w3 = canonical3NodesWeights[k];
                // Changing third coordinate and weight to adapt to the integration interval
                double a3 = integrationInterval.geta3();
                double b3 = integrationInterval.getb3();
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
                                                        DataVecDeriv& dataf, const DataVecCoord& datax, const DataVecDeriv& v)
{
}

/// Computes current displacement, expressed in the beam central frame (local), with respect to
/// the beam rest position.
template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::computeLocalDisplacement(const VecCoord& x, Vec6& U0Local, Vec6& U1Local,
                                                                        unsigned int beamIndex,
                                                                        unsigned int node0Idx, unsigned int node1Idx,
                                                                        bool updateBeamMatrices)
{
    ///////////// new : Calcul du repère local de la beam & des transformations adequates///////////////
    Transform global_H_local0, global_H_local1;

    /// 1. get the current transform of the beam:
    dmsg_info() << "in addForce";
    l_interpolation->computeTransform2(beamIndex, global_H_local0, global_H_local1, x);

    /// 2. Computes the frame of the beam based on the spline interpolation:
    Transform global_H_local;
    Real baryX = 0.5;
    Real L = l_interpolation->getLength(beamIndex);

    l_interpolation->InterpolateTransformUsingSpline(global_H_local, baryX, global_H_local0, global_H_local1, L);

    Transform local_H_local0 = global_H_local.inversed() * global_H_local0;
    Transform local_H_local1 = global_H_local.inversed() * global_H_local1;

    /// 3. Computes the transformation from the DOF (in global frame) to the node's local frame DOF0global_H_Node0local and DOF1global_H_Node1local
    Transform DOF0_H_local0, DOF1_H_local1;
    l_interpolation->getDOFtoLocalTransform(beamIndex, DOF0_H_local0, DOF1_H_local1);

    /// 4. Computes the adequate transformation
    Transform global_R_DOF0(Vec3(0, 0, 0), x[node0Idx].getOrientation());
    Transform global_R_DOF1(Vec3(0, 0, 0), x[node1Idx].getOrientation());
    /// - rotation due to the optional transformation
    global_H_local0 = global_R_DOF0 * DOF0_H_local0;
    global_H_local1 = global_R_DOF1 * DOF1_H_local1;

    Transform DOF0global_H_Node0local, DOF1global_H_Node1local;

    DOF0global_H_Node0local.set(global_H_local0.getOrigin(), global_H_local.getOrientation());
    DOF1global_H_Node1local.set(global_H_local1.getOrigin(), global_H_local.getOrientation());

    //TODO(dmarchal 2017-05-17) Please specify who/when this will be done
    //TODO A verifier : global_H_local0.getOrigin() == x[node0Idx].getOrientation().rotate(DOF0_H_local0.getOrigin())

    /// compute the current local displacement of the beam (6dofs)
    /// 1. get the rest transformation from local to 0 and local to 1
    Transform local_H_local0_rest, local_H_local1_rest;
    l_interpolation->getSplineRestTransform(beamIndex, local_H_local0_rest, local_H_local1_rest);

    ///2. computes the local displacement of 0 and 1 in frame local:
    SpatialVector u0 = local_H_local0.CreateSpatialVector() - local_H_local0_rest.CreateSpatialVector();
    SpatialVector u1 = local_H_local1.CreateSpatialVector() - local_H_local1_rest.CreateSpatialVector();

    /// 3. put the result in a Vec6
    for (unsigned int i = 0; i < 3; i++)
    {
        U0Local[i] = u0.getLinearVelocity()[i];
        U0Local[i + 3] = u0.getAngularVelocity()[i];
        U1Local[i] = u1.getLinearVelocity()[i];
        U1Local[i + 3] = u1.getAngularVelocity()[i];
    }

    if (d_reinforceLength.getValue())
    {
        Vec3 P0, P1, P2, P3;
        Real length;
        Real rest_length = l_interpolation->getLength(beamIndex);
        l_interpolation->getSplinePoints(beamIndex, x, P0, P1, P2, P3);
        l_interpolation->computeActualLength(length, P0, P1, P2, P3);

        U0Local[0] = (-length + rest_length) / 2;
        U1Local[0] = (length - rest_length) / 2;
    }

    if (!d_useShearStressComputation.getValue())
    {
        /////////////////// TEST //////////////////////
        /// test: correction due to spline computation;
        Vec3 ResultNode0, ResultNode1;
        l_interpolation->computeStrechAndTwist(beamIndex, x, ResultNode0, ResultNode1);

        Real ux0 = -ResultNode0[0] + l_interpolation->getLength(beamIndex) / 2;
        Real ux1 = ResultNode1[0] - l_interpolation->getLength(beamIndex) / 2;

        U0Local[0] = ux0;
        U1Local[0] = ux1;

        U0Local[3] = -ResultNode0[2];
        U1Local[3] = ResultNode1[2];

        //////////////////////////////////////////////////
    }

    if (updateBeamMatrices)
    {
        ///find the beamMatrices:
        BeamLocalMatrices* beamMatrices = &m_localBeamMatrices[beamIndex];//new BeamLocalMatrices();

        /// compute Adjoint Matrices:
        beamMatrices->m_A0Ref = DOF0global_H_Node0local.inversed().getAdjointMatrix();
        beamMatrices->m_A1Ref = DOF1global_H_Node1local.inversed().getAdjointMatrix();

        /// compute the local mass matrices
        if (d_computeMass.getValue())
        {
            computeMass(beamIndex, (*beamMatrices));

        }

        /// compute the local stiffness matrices
        computeStiffness(beamIndex, (*beamMatrices));
    }

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



/////////////////////////////////////
/// GaussPoint3
/////////////////////////////////////

template <class DataTypes>
AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::GaussPoint3(Real x, Real y, Real z, Real w1, Real w2, Real w3)
{
    m_coordinates = { x, y, z };
    m_weights = { w1, w2, w3 };
    m_mechanicalState = MechanicalState::ELASTIC; //By default, before any deformation occurs
    m_prevStress = Vec9(); //By default, no deformation => 0 stress tensor
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::getGradN() const -> const Matrix9x12&
{
    return m_gradN;
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::setGradN(Matrix9x12 gradN)
{
    m_gradN = gradN;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::getMechanicalState() const -> const MechanicalState
{
    return m_mechanicalState;
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::setMechanicalState(MechanicalState newState)
{
    m_mechanicalState = newState;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::getPrevStress() const -> const Vec9&
{
    return m_prevStress;
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::setPrevStress(Vec9 newStress)
{
    m_prevStress = newStress;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::getWeights() const -> const Vec3&
{
    return m_weights;
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::setWeights(Vec3 weights)
{
    m_weights = weights;
}


/////////////////////////////////////
/// Interval3
/////////////////////////////////////

template <class DataTypes>
AdaptivePlasticBeamForceField<DataTypes>::Interval3::Interval3()
{
    //By default, integration is considered over [-1,1]*[-1,1]*[-1,1].
    m_a1 = m_a2 = m_a3 = -1;
    m_b1 = m_b2 = m_b3 = 1;
}

template <class DataTypes>
AdaptivePlasticBeamForceField<DataTypes>::Interval3::Interval3(Real a1, Real b1, Real a2, Real b2, Real a3, Real b3)
{
    m_a1 = a1;
    m_b1 = b1;
    m_a2 = a2;
    m_b2 = b2;
    m_a3 = a3;
    m_b3 = b3;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::Interval3::geta1() const -> Real
{
    return m_a1;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::Interval3::getb1() const -> Real
{
    return m_b1;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::Interval3::geta2() const -> Real
{
    return m_a2;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::Interval3::getb2() const -> Real
{
    return m_b2;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::Interval3::geta3() const -> Real
{
    return m_a3;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::Interval3::getb3() const -> Real
{
    return m_b3;
}

} // namespace _adaptiveplasticbeamforcefield_
} // namespace sofa::plugin::beamadapter::component::forcefield