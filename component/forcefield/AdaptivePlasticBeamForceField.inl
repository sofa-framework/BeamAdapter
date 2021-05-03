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
AdaptivePlasticBeamForceField<DataTypes>::AdaptivePlasticBeamForceField() :
    d_poissonRatio(initData(&d_poissonRatio, Real(0.3), "poissonRatio",
        "Value of the poisson ratio for the considered material")),
    d_youngModulus(initData(&d_youngModulus, "youngModulus",
        "Value of the young modulus for the considered material")),
    d_initialYieldStress(initData(&d_initialYieldStress, "initialYieldStress",
        "Yield stress of the considered material, prior to any elastic deformation")),
    d_plasticModulus(initData(&d_plasticModulus, "plasticModulus",
        "Approximation of the plastic modulus as a constant. Can be deduced from a generic law such as Ramberg-Osgood's"))
{
    if (d_initialYieldStress.getValue() < 0)
    {
        msg_error() << "yield Stress should be positive. Please provide a positive yield"
                    << "stress value for the considered material";
    }

}

template <class DataTypes>
AdaptivePlasticBeamForceField<DataTypes>::~AdaptivePlasticBeamForceField()
{
    //TO DO : should m_gaussPoints, m_integrationIntervals, m_beamMechanicalStates be deallocated explicitly ?
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

    //Initialisation of the lastPos field with the rest position
    m_lastPos = mstate->read(core::ConstVecCoordId::restPosition())->getValue();

    //Initialisation of generalised Hooke's law
    const Real E = d_youngModulus.getValue();
    const Real nu = d_poissonRatio.getValue();
    m_genHookesLaw[0][0] = m_genHookesLaw[4][4] = m_genHookesLaw[8][8] = 1 - nu;
    m_genHookesLaw[0][4] = m_genHookesLaw[0][8] = m_genHookesLaw[4][0] = m_genHookesLaw[8][0] = nu;
    m_genHookesLaw[4][8] = m_genHookesLaw[8][4] = nu;

    m_genHookesLaw[1][1] = m_genHookesLaw[2][2] = m_genHookesLaw[3][3] = (1 - 2*nu) / 2;
    m_genHookesLaw[5][5] = m_genHookesLaw[6][6] = m_genHookesLaw[7][7] = (1 - 2*nu) / 2;

    m_genHookesLaw[1][3] = m_genHookesLaw[2][6] = m_genHookesLaw[5][7] = (1 - 2*nu) / 2;
    m_genHookesLaw[3][1] = m_genHookesLaw[6][2] = m_genHookesLaw[7][5] = (1 - 2*nu) / 2;

    m_genHookesLaw *= E / ((1 + nu)*(1 - 2*nu));

    //Initialisation of Gauss points and integration intervals
    m_gaussPoints.clear();
    m_gaussPoints.resize(numBeams);
    m_integrationIntervals.clear();
    m_integrationIntervals.resize(numBeams);

    for (unsigned int b = 0; b < numBeams; b++)
    {
        BeamLocalMatrices* beamMatrices = &m_localBeamMatrices[b];

        initialiseInterval(b, m_integrationIntervals, beamGeometryParams[b]);
        initialiseGaussPoints(b, m_gaussPoints, m_integrationIntervals[b]);
        updateTangentStiffness(b, *beamMatrices); //Prior to plastic deformation, computes an elastic stiffness with Cep = C
    }

    //Initialisaiton of beam mechanical states to ELASTIC
    m_beamMechanicalStates.clear();
    m_beamMechanicalStates.resize(numBeams);
    std::fill(m_beamMechanicalStates.begin(), m_beamMechanicalStates.end(), MechanicalState::ELASTIC);

    //Initialisation of the comparison threshold for stress tensor norms to 0.
    // Plasticity computation requires to basically compare stress tensor norms to 0.
    // As stress norm values can vary of several orders of magnitude, depending on the
    // considered materials and/or applied forces, this comparison has to be carried out
    // carefully.
    // The idea here is to use the initialYieldStress of the material, and the
    // available precision limit (e.g. std::numeric_limits<double>::epsilon()).
    // We rely on the value of the initial Yield stress, as we can expect plastic
    // deformation to occur inside a relatively small intervl of stresses around this value.
    const int orderOfMagnitude = d_initialYieldStress.getValue(); //Should use std::abs, but d_initialYieldStress > 0
    m_stressComparisonThreshold = std::numeric_limits<double>::epsilon()*orderOfMagnitude;
}


template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::initialiseGaussPoints(int beam, vector<beamGaussPoints>& gaussPoints,
                                                                     const Interval3& integrationInterval)
{
    //Gaussian nodes coordinates and weights for a 1D integration on [-1,1]
    const double sqrt3_5 = helper::rsqrt(3.0 / 5);
    Vec3 canonical3NodesCoordinates = { -sqrt3_5, 0, sqrt3_5 };
    Vec3 canonical3NodesWeights = { 5.0 / 9, 8.0 / 9, 5.0 / 9 };

    const BeamPlasticInterpolation<DataTypes>::BeamSection beamSection = l_interpolation->getBeamSection(beam);
    double Iy = beamSection._Iy;
    double Iz = beamSection._Iz;
    double A = beamSection._A;
    const double L = l_interpolation->getLength(beam);

    Real E = d_youngModulus.getValue();
    Real nu = d_poissonRatio.getValue();

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
                newGaussPoint.setGradN( computeGradN(xChanged, yChanged, zChanged, L, A, Iy, Iz, E, nu) );
                newGaussPoint.setYieldStress( d_initialYieldStress.getValue() );
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
void AdaptivePlasticBeamForceField<DataTypes>::integrateBeam(beamGaussPoints& gaussPoints, LambdaType integrationFun)
{
    //Apply a generic (lambda) integration function to each Gauss point of a beam element
    for (unsigned int gp = 0; gp < gaussPoints.size(); gp++)
    {
        integrationFun(gaussPoints[gp]);
    }
}


template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::computeGradN(Real x, Real y, Real z, Real L,
                                                            Real A, Real Iy, Real Iz, Real E, Real nu,
                                                            Real kappaY, Real kappaZ) -> Matrix9x12
{
    Matrix9x12 gradN = Matrix9x12();
    Real xi = x / L;
    Real eta = y / L;
    Real zeta = z / L;

    Real L2 = L * L;
    Real G = E / (2.0 * (1.0 + nu));

    Real phiY, phiZ;
    if (A == 0)
    {
        phiY = 0.0;
        phiZ = 0.0;
    }
    else
    {
        phiY = (12.0 * E * Iy / (kappaZ * G * A * L2));
        phiZ = (12.0 * E * Iz / (kappaY * G * A * L2));
    }
    Real phiYInv = (1 / (1 + phiY));
    Real phiZInv = (1 / (1 + phiZ));

    //Row 0
    gradN[0][0] = -1 / L;
    gradN[0][1] = (phiZInv * 6 * eta * (1 - 2 * xi)) / L;
    gradN[0][2] = (phiYInv * 6 * zeta * (1 - 2 * xi)) / L;
    gradN[0][3] = 0;
    gradN[0][4] = phiYInv * zeta * (6 * xi - 4 - phiY);
    gradN[0][5] = phiZInv * eta * (4 - 6 * xi + phiZ);
    gradN[0][6] = 1 / L;
    gradN[0][7] = (phiZInv * 6 * eta * (2 * xi - 1)) / L;
    gradN[0][8] = (phiYInv * 6 * zeta * (2 * xi - 1)) / L;
    gradN[0][9] = 0;
    gradN[0][10] = phiYInv * zeta * (6 * xi - 2 + phiY);
    gradN[0][11] = phiZInv * eta * (2 - 6 * xi - phiZ);

    //Rows 1 and 3
    gradN[1][0] = 0.0;
    gradN[1][1] = -(phiZInv * phiZ) / (2 * L);
    gradN[1][2] = 0.0;
    gradN[1][3] = zeta / 2;
    gradN[1][4] = 0.0;
    gradN[1][5] = -(phiZInv * phiZ) / 4;
    gradN[1][6] = 0.0;
    gradN[1][7] = (phiZInv * phiZ) / (2 * L);
    gradN[1][8] = 0.0;
    gradN[1][9] = -zeta / 2;
    gradN[1][10] = 0.0;
    gradN[1][11] = -(phiZInv * phiZ) / 4;

    gradN[3][0] = 0.0;
    gradN[3][1] = -(phiZInv * phiZ) / (2 * L);
    gradN[3][2] = 0.0;
    gradN[3][3] = zeta / 2;
    gradN[3][4] = 0.0;
    gradN[3][5] = -(phiZInv * phiZ) / 4;
    gradN[3][6] = 0.0;
    gradN[3][7] = (phiZInv * phiZ) / (2 * L);
    gradN[3][8] = 0.0;
    gradN[3][9] = -zeta / 2;
    gradN[3][10] = 0.0;
    gradN[3][11] = -(phiZInv * phiZ) / 4;

    //Rows 2 and 6
    gradN[2][0] = gradN[2][1] = 0.0;
    gradN[2][2] = -(phiYInv * phiY) / (2 * L);
    gradN[2][3] = -eta / 2;
    gradN[2][4] = (phiYInv * phiY) / 4;
    gradN[2][5] = gradN[2][6] = gradN[2][7] = 0.0;
    gradN[2][8] = (phiYInv * phiY) / (2 * L);
    gradN[2][9] = eta / 2;
    gradN[2][10] = (phiYInv * phiY) / 4;
    gradN[2][11] = 0.0;

    gradN[6][0] = gradN[6][1] = 0.0;
    gradN[6][2] = -(phiYInv * phiY) / (2 * L);
    gradN[6][3] = -eta / 2;
    gradN[6][4] = (phiYInv * phiY) / 4;
    gradN[6][5] = gradN[6][6] = gradN[6][7] = 0.0;
    gradN[6][8] = (phiYInv * phiY) / (2 * L);
    gradN[6][9] = eta / 2;
    gradN[6][10] = (phiYInv * phiY) / 4;
    gradN[6][11] = 0.0;

    //Rows 4, 5, 7, 8
    gradN[4][0] = gradN[4][1] = gradN[4][2] = gradN[4][3] = 0.0;
    gradN[4][4] = gradN[4][5] = gradN[4][6] = gradN[4][7] = 0.0;
    gradN[4][8] = gradN[4][9] = gradN[4][10] = gradN[4][11] = 0.0;

    gradN[5][0] = gradN[5][1] = gradN[5][2] = gradN[5][3] = 0.0;
    gradN[5][4] = gradN[5][5] = gradN[5][6] = gradN[5][7] = 0.0;
    gradN[5][8] = gradN[5][9] = gradN[5][10] = gradN[5][11] = 0.0;

    gradN[7][0] = gradN[7][1] = gradN[7][2] = gradN[7][3] = 0.0;
    gradN[7][4] = gradN[7][5] = gradN[7][6] = gradN[7][7] = 0.0;
    gradN[7][8] = gradN[7][9] = gradN[7][10] = gradN[7][11] = 0.0;

    gradN[8][0] = gradN[8][1] = gradN[8][2] = gradN[8][3] = 0.0;
    gradN[8][4] = gradN[8][5] = gradN[8][6] = gradN[8][7] = 0.0;
    gradN[8][8] = gradN[8][9] = gradN[8][10] = gradN[8][11] = 0.0;

    return gradN;
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
    SOFA_UNUSED(v);

    VecDeriv& f = *dataf.beginEdit();
    const VecCoord& x = datax.getValue();

    f.resize(x.size()); // current content of the vector will remain the same (http://www.cplusplus.com/reference/vector/vector/resize/)

    unsigned int numBeams = l_interpolation->getNumBeams();
    m_localBeamMatrices.resize(numBeams);

    if (d_computeMass.getValue())
        computeGravityVector();

    const vector<beamGaussPoints>& gaussPoints = m_gaussPoints;

    /// Computes internal forces.
    for (unsigned int b = 0; b < numBeams; b++)
    {
        ///find the indices of the nodes
        unsigned int node0Idx, node1Idx;
        l_interpolation->getNodeIndices(b, node0Idx, node1Idx);

        /// find the beamMatrices:
        BeamLocalMatrices* beamMatrices = &m_localBeamMatrices[b];

        /// IF RIGIDIFICATION: no stiffness forces:
        if (node0Idx == node1Idx)
            continue;

        /// Computes local displacement increment
        Vec6 U0Local, U1Local, U0LocalLast, U1LocalLast, DeltaU0Local, DeltaU1Local;
        computeLocalDisplacement(x, U0Local, U1Local, b, node0Idx, node1Idx, true);
        computeLocalDisplacement(m_lastPos, U0LocalLast, U1LocalLast, b, node0Idx, node1Idx);
        DeltaU0Local = U0Local - U0LocalLast;
        DeltaU1Local = U1Local - U1LocalLast;

        /////////////////////////////COMPUTATION OF THE STIFFNESS FORCE

        // Computation of the new stress point, through material point iterations
        // Cf detail comment of the algorithm in computeStressIncrement

        bool isPlasticBeam = false;
        Vec12 localForces = Vec12();

        typedef std::function<void(GaussPoint3&)> IntegrationLambda;
        IntegrationLambda computeStress = [&](GaussPoint3& gp)
        {
            MechanicalState mechanicalState = gp.getMechanicalState();

            //Strain
            Vec12 displacementIncrement = Vec12();
            for (int i=0; i < 6; i++)
            {
                displacementIncrement[i] = DeltaU0Local[i];
                displacementIncrement[i+6] = DeltaU1Local[i];
            }
            Vec9 strainIncrement = gp.getGradN() * displacementIncrement;

            //Stress
            Vec9 newStressPoint = Vec9();
            computeStressIncrement(gp, newStressPoint,
                                   strainIncrement, mechanicalState);

            // Update the mechanical state information
            mechanicalState = gp.getMechanicalState();

            isPlasticBeam = isPlasticBeam || (mechanicalState == MechanicalState::PLASTIC);

            gp.setPrevStress(newStressPoint);

            Vec3 gpWeights = gp.getWeights();
            localForces += (gpWeights[1] * gpWeights[2] * gpWeights[3]) * gp.getGradN().transposed() * newStressPoint;
        };

        //Actually run gaussian integration on the beams Gauss points
        integrateBeam(m_gaussPoints[b], computeStress);

        // Updates the beam mechanical state information
        if (!isPlasticBeam)
            // TO DO: prior to any plastic deformation, this should be ELASTIC. But it makes no difference in the beam implementation
            m_beamMechanicalStates[b] = MechanicalState::POSTPLASTIC;

        //Update the tangent stiffness matrix with the new computed stresses
        //This matrix will then be used in addDForce and addKToMatrix methods
        if (isPlasticBeam)
            updateTangentStiffness(b, *beamMatrices);

        //Retrieve local forces from the gaussian integration
        //TO DO : better way to do this ? std::copy ?
        Vec6 f0 = { localForces[0], localForces[1], localForces[2], localForces[3], localForces[4], localForces[5] };
        Vec6 f1 = { localForces[6], localForces[7], localForces[8], localForces[9], localForces[10], localForces[11] };

        /// Compute the force in the global frame
        Vec6 F0_ref = beamMatrices->m_A0Ref.multTranspose(f0);
        Vec6 F1_ref = beamMatrices->m_A1Ref.multTranspose(f1);

        /// Add this force to vector f
        for (unsigned int i = 0; i < 6; i++)
        {
            f[node0Idx][i] -= F0_ref[i];
            f[node1Idx][i] -= F1_ref[i];
        }

        /// Compute local mass matrix (as in AdaptiveBeamForceFieldAndMass)
        if (d_computeMass.getValue())
            computeMass(b, (*beamMatrices));

    } // end for (unsigned int b = 0; b < numBeams; b++)

    /// Add gravity (as in AdaptiveBeamForceFieldAndMass)
    if (d_computeMass.getValue())
        addMDx(mparams, dataf, d_dataG, 1.0);

    // Save the current positions as a record for the next time step.
    //TO DO: check if this is copy operator
    m_lastPos = x;

    dataf.endEdit();
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


template< class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::computeStressIncrement(GaussPoint3& gp,
                                                                      Vec9& newStressPoint,
                                                                      const Vec9& strainIncrement,
                                                                      MechanicalState& pointMechanicalState)
{
    /// This method implements the radial return algorithm, as in "Numerical Implementation of
    /// Constitutive models: Rate-independent Deviatoric Plasticity", T.J.R. Hugues, 1984.
    /// The idea is to compute the stress increment in two steps : a purely elastic step, in
    /// which all deformation is considered elastic, and a plastic correction step, is
    /// deformation was actually important enough to generate plasticity.
    /// The plasticity model used in the computation is a Von Mises-Hill plasticity with
    /// linear mixed hardening.
    /// NB: we consider that the yield function and the plastic flow are equal (f=g). This
    /// corresponds to an associative flow rule (for plasticity).

    const Matrix9x9& C = m_genHookesLaw;

    // First step = computation of the elastic predictor, as if deformation was entirely elastic
    const MechanicalState mechanicalState = gp.getMechanicalState();

    Vec9 elasticIncrement = C*strainIncrement;
    Vec9 elasticPredictorStress = gp.getPrevStress() + elasticIncrement;

    const Vec9& backStress = gp.getBackStress();
    const Real yieldStress = gp.getYieldStress();

    if (vonMisesYield(elasticPredictorStress, backStress, yieldStress) < m_stressComparisonThreshold)
    {
        // The Gauss point is in elastic state: the back stress and yield stress
        // remain constant, and the new stress is equal to the trial stress.
        newStressPoint = elasticPredictorStress;

        // If the Gauss point was initially plastic, we update its mechanical state
        if (mechanicalState == MechanicalState::PLASTIC)
            gp.setMechanicalState(MechanicalState::POSTPLASTIC);
    }
    else
    {
        // If the Gauss point was initially elastic, we update its mechanical state
        if (mechanicalState == MechanicalState::POSTPLASTIC || mechanicalState == MechanicalState::ELASTIC)
            gp.setMechanicalState(MechanicalState::POSTPLASTIC);

        Vec9 shiftedDeviatoricElasticPredictor = deviatoricStress(elasticPredictorStress - backStress);

        // Gradient of the Von Mises yield function is colinear to the deviatoric stress tensor.
        // Thus we can compute the yield surface normal using the deviatoric stress.
        // For the Von Mises yield function, the normal direction to the yield surface doesn't
        // change between the elastic predictor and it's projection on the yield surface
        Real shiftDevElasticPredictorNorm = shiftedDeviatoricElasticPredictor.norm();
        Vec9 N = shiftedDeviatoricElasticPredictor / shiftDevElasticPredictorNorm;

        // Indicates the proportion of Kinematic vs isotropic hardening. beta=0 <=> kinematic, beta=1 <=> isotropic
        const Real beta = 0.5;

        const Real E = d_youngModulus.getValue();
        const Real nu = d_poissonRatio.getValue();
        const Real mu = E / ( 2*(1 + nu) ); // Lame coefficient

        // Plastic modulus
        const Real H = d_plasticModulus.getValue();

        // Computation of the plastic multiplier
        const double sqrt2 = helper::rsqrt(2.0);
        const double sqrt3 = helper::rsqrt(3.0);
        const double sqrt6 = sqrt2 * sqrt3;
        Real plasticMultiplier = (shiftDevElasticPredictorNorm - (sqrt2 / sqrt3) * yieldStress) / ( mu*sqrt6 *( 1 + H/(3*mu) ) );

        // Updating plastic variables
        newStressPoint = elasticPredictorStress - sqrt6*mu*plasticMultiplier * N;

        Real newYieldStress = yieldStress + beta * H * plasticMultiplier;
        gp.setYieldStress(newYieldStress);
        Vec9 newBackStress = backStress + (sqrt2 / sqrt3) * (1-beta) * H * plasticMultiplier * N;
        gp.setBackStress(newBackStress);

        Vec9 newPlasticStrain = gp.getPlasticStrain() + (sqrt3/sqrt2) * plasticMultiplier * N;
        gp.setPlasticStrain(newPlasticStrain);
        gp.setEffectivePlasticStrain(gp.getEffectivePlasticStrain() + plasticMultiplier);
    }
}


template< class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::updateTangentStiffness(unsigned int beamIndex,
                                                                      BeamLocalMatrices& beamLocalMatrices)
{

    // Switching to the 12x12 format for gaussian quadrature
    Matrix12x12 Kt = Matrix12x12();

    const Matrix9x9& C = m_genHookesLaw;
    const Real E = d_youngModulus.getValue();
    const Real nu = d_poissonRatio.getValue();

    typedef std::function<void(GaussPoint3&)> IntegrationLambda;
    IntegrationLambda computeTangentStiffness = [&](GaussPoint3& gp)
    {
        Vec9 currentStressPoint = gp.getPrevStress(); //Updated in computeStressIncrement

        Real H = d_plasticModulus.getValue();

        Matrix9x12 gradN = gp.getGradN();
        Matrix9x9 Cep = Matrix9x9();
        // Cep
        Vec9 gradient = vonMisesGradient(currentStressPoint);

        if ( gradient.norm() < m_stressComparisonThreshold
             || gp.getMechanicalState() != MechanicalState::PLASTIC )
            Cep = C; //TO DO: is that correct ?
        else
        {
            Vec9 N = helper::rsqrt(2.0 / 3.0) * gradient;
            Vec9 CN = C * N;
            //Conversion to matrix, TO DO: better way ? Eigen ? Use only matrices ?
            Mat<9, 1, Real> matCN = Mat<9, 1, Real>();
            for (int i = 0; i < CN.size(); i++)
                matCN[i][0] = CN[i];
            // NtC = (NC)t because of C symmetry
            Cep = C - ( matCN * matCN.transposed() ) / (N * CN + (2.0 / 3.0) * H); //TO DO: check that * operator is actually dot product
        }

        Vec3 gpWeights = gp.getWeights();
        Kt += gpWeights[1]*gpWeights[2]*gpWeights[3] * gradN.transposed()*Cep*gradN;
    }; //end computeTangentStiffness

    //Actually runs gaussian integration on the beams Gauss points
    integrateBeam(m_gaussPoints[beamIndex], computeTangentStiffness);

    //Update beam tangent stiffness
    Kt.getsub(0, 0, beamLocalMatrices.m_K00);
    Kt.getsub(0, 6, beamLocalMatrices.m_K01);
    Kt.getsub(6, 0, beamLocalMatrices.m_K10);
    Kt.getsub(6, 6, beamLocalMatrices.m_K11);
}

//----- Implementation of the Von Mises yield function -----//

template< class DataTypes>
typename AdaptivePlasticBeamForceField<DataTypes>::Vec9 AdaptivePlasticBeamForceField<DataTypes>::deviatoricStress(const Vec9& stressTensor)
{
    // Returns the deviatoric stress from a given stress tensor in Voigt notation

    Vec9 deviatoricStress = stressTensor;
    double mean = (stressTensor[0] + stressTensor[4] + stressTensor[8]) / 3.0;

    deviatoricStress[0] -= mean;
    deviatoricStress[4] -= mean;
    deviatoricStress[8] -= mean;

    return deviatoricStress;
}

template< class DataTypes>
typename AdaptivePlasticBeamForceField<DataTypes>::Real AdaptivePlasticBeamForceField<DataTypes>::equivalentStress(const Vec9& stressTensor)
{
    // Direct computation of the equivalent stress. We use the fact that the tensor is symmetric
    Real sigmaXX = stressTensor[0];
    Real sigmaXY = stressTensor[1];
    //Real sigmaXZ = stressTensor[2];
    //Real sigmaYX = stressTensor[3];
    Real sigmaYY = stressTensor[4];
    Real sigmaYZ = stressTensor[5];
    Real sigmaZX = stressTensor[6];
    //Real sigmaZY = stressTensor[7];
    Real sigmaZZ = stressTensor[8];

    double aux1 = 0.5 * ((sigmaXX - sigmaYY) * (sigmaXX - sigmaYY) + (sigmaYY - sigmaZZ) * (sigmaYY - sigmaZZ) + (sigmaZZ - sigmaXX) * (sigmaZZ - sigmaXX));
    double aux2 = 3.0 * (sigmaYZ * sigmaYZ + sigmaZX * sigmaZX + sigmaXY * sigmaXY);

    return helper::rsqrt(aux1 + aux2);
}

template< class DataTypes>
typename AdaptivePlasticBeamForceField<DataTypes>::Real AdaptivePlasticBeamForceField<DataTypes>::vonMisesYield(const Vec9& stressTensor,
                                                                                                                const Vec9& backStress,
                                                                                                                const Real yieldStress)
{
    return equivalentStress(stressTensor-backStress) - yieldStress;
}

template< class DataTypes>
typename AdaptivePlasticBeamForceField<DataTypes>::Vec9 AdaptivePlasticBeamForceField<DataTypes>::vonMisesGradient(const Vec9& stressTensor)
{
    // NB: this gradient represent the normal to the yield surface
    // in case the Von Mises yield criterion is used. The norm of the
    // gradient is sqrt(3/2): it has to be multiplied by sqrt(2/3)
    // to give the unit normal to the yield surface

    //Direct computation: TO DO: compare with deviatoric computation
    //Vec9 gradient = Vec9();

    //if (stressTensor.norm() < m_stressComparisonThreshold.getValue())
    //    return gradient; //TO DO: is that correct ?

    //Real sigmaXX = stressTensor[0];
    //Real sigmaXY = stressTensor[1];
    ////Real sigmaXZ = stressTensor[2];
    ////Real sigmaYX = stressTensor[3];
    //Real sigmaYY = stressTensor[4];
    //Real sigmaYZ = stressTensor[5];
    //Real sigmaZX = stressTensor[6];
    ////Real sigmaZY = stressTensor[7];
    //Real sigmaZZ = stressTensor[8];

    //gradient[0] = 2*sigmaXX - sigmaYY - sigmaZZ;
    //gradient[4] = 2*sigmaYY - sigmaZZ - sigmaXX;
    //gradient[8] = 2 * sigmaZZ - sigmaXX - sigmaYY;
    //gradient[1] = gradient[3] = 3 * sigmaXY;
    //gradient[2] = gradient[6] = 3 * sigmaZX;
    //gradient[5] = gradient[7] = 3 * sigmaYZ;

    //Real sigmaEq = equivalentStress(stressTensor);
    //gradient *= 1 / (2*sigmaEq);
    //return gradient;

    //Deviatoric-based computation: TO DO: compare with direct computation
    Real sigmaEq = equivalentStress(stressTensor);
    return (3 / (2 * sigmaEq)) * deviatoricStress(stressTensor);
}

//----------------------------------------------------------//

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::computeStiffness(int beam, BeamLocalMatrices& beamLocalMatrices)
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
    m_backStress = Vec9(); //By default, no plastic deformation => back stress is 0
    m_yieldStress = 0; //Changed by initialiseGaussPoints, depends on the material
    m_plasticStrain = Vec9(); //By default, no plastic deformation => no history
    m_effectivePlasticStrain = 0; //By default, no plastic deformation => no history
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

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::getBackStress() const -> const Vec9&
{
    return m_backStress;
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::setBackStress(Vec9 backStress)
{
    m_backStress = backStress;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::getYieldStress() const ->const Real
{
    return m_yieldStress;
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::setYieldStress(Real yieldStress)
{
    m_yieldStress = yieldStress;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::getPlasticStrain() const -> const Vec9&
{
    return m_plasticStrain;
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::setPlasticStrain(Vec9 plasticStrain)
{
    m_plasticStrain = plasticStrain;
}

template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::getEffectivePlasticStrain() const ->const Real
{
    return m_effectivePlasticStrain;
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::setEffectivePlasticStrain(Real effectivePlasticStrain)
{
    m_effectivePlasticStrain = effectivePlasticStrain;
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