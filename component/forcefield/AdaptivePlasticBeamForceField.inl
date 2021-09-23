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
using sofa::type::Quat;

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
        l_interpolation.set(dynamic_cast<BaseContext*>(this->getContext())->get<BInterpolation>(BaseContext::Local));

    if (!l_interpolation)
        msg_error() << "No Beam Interpolation found, the component can not work.";

    ForceField<DataTypes>::init();

    // The interpolation component has to be initiated before, in order to get
    // the information required for Gauss points initialisation. But with this method,
    // bwdInit() will be called twice on the interpolation component.
    // TO DO : is there a better way to guarantee the order of initialisation ?
    l_interpolation->bwdInit();

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

    //Initialisation of Gauss points, integration intervals, and local matrices
    m_localBeamMatrices.clear();
    m_localBeamMatrices.resize(numBeams);
    m_gaussPoints.clear();
    m_gaussPoints.resize(numBeams);
    m_integrationIntervals.clear(); // No need to resize as Intervals will be pushed_back

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

    const BeamInterpolation<DataTypes>::BeamSection beamSection = l_interpolation->getBeamSection(beam);
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
                newGaussPoint.setNx( computeNx(xChanged, yChanged, zChanged, L, A, Iy, Iz, E, nu) );
                newGaussPoint.setYieldStress( d_initialYieldStress.getValue() );
                gaussPoints[beam][gaussPointIt] = newGaussPoint;
                gaussPointIt++;
            }
        }
    }
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::initialiseInterval(int beam, vector<Interval3>& integrationIntervals,
                                                                  const BeamGeometry& beamGeometryParams)
{
    if (beamGeometryParams.sectionShape == "rectangular")
    {
        Real L = beamGeometryParams._L;
        Real Ly = beamGeometryParams._Ly;
        Real Lz = beamGeometryParams._Lz;

        // Integration interval definition for a local frame at the centre of the beam
        integrationIntervals.push_back(Interval3(-L/2, L/2, -Ly/2, Ly/2, -Lz/2, Lz/2));
    }
    else if (beamGeometryParams.sectionShape == "elliptic")
    {
        // TO DO: implement quadrature method for elliptic cross section
        msg_error() << "Quadrature method for " << beamGeometryParams.sectionShape
            << " shape cross section has not been implemented yet. Methods for rectangular cross sections are available";
    }
    else if (beamGeometryParams.sectionShape == "circular")
    {
        Real L = beamGeometryParams._L;
        Real r = beamGeometryParams._r;
        if (beamGeometryParams._rInner != 0) // NB: this test should work, as _rInner is set to 0 for every other case
        {
            //TO DO: implement quadrature method for a disc and a hollow-disc cross section
            msg_error() << "Quadrature method for " << beamGeometryParams.sectionShape
                << " shape cross section has not been implemented yet. Methods for rectangular and circular cross sections are available";
        }

        // Integration interval definition for a local frame at the centre of the beam
        integrationIntervals.push_back(Interval3(-L / 2, L / 2, -r / 2, r / 2, -r / 2, r / 2));
    }
    else
    {
        msg_error() << "Quadrature method for " << beamGeometryParams.sectionShape
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
auto AdaptivePlasticBeamForceField<DataTypes>::computeNx(Real x, Real y, Real z, Real L,
                                                         Real A, Real Iy, Real Iz, Real E, Real nu,
                                                         Real kappaY, Real kappaZ) -> Matrix3x12
{
    Matrix3x12 Nx = Matrix3x12(); // Sets each element to 0
    Real xi = x / L;
    Real eta = y / L;
    Real zeta = z / L;

    Real xi2 = xi * xi;
    Real xi3 = xi * xi * xi;

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
    Real phiYInv = ( 1 / (phiY - 1) );
    Real phiZInv = ( 1 / (1 + phiZ) );

    Nx[0][0] = 1 - xi;
    Nx[0][1] = ( 1 + 0.5*phiZInv - 6*phiZInv * xi2 ) * eta;
    Nx[0][2] = ( 1 - 0.5*phiYInv + 6*phiYInv * xi2 ) * zeta;
    //Nx[0][3] = 0;
    Nx[0][4] = ( (phiYInv / 4) - xi - 3*phiYInv * xi2 ) * z;
    Nx[0][5] = ( (phiZInv / 4) + xi - 3*phiZInv * xi2 ) * y;
    Nx[0][6] = xi;
    Nx[0][7] = -Nx[0][1];
    Nx[0][8] = -Nx[0][2];
    //Nx[0][9] = 0;
    Nx[0][10] = Nx[0][4] + 2*xi*z;
    Nx[0][11] = Nx[0][5] - 2*xi*y;

    //Nx[1][0] = 0;
    Nx[1][1] = 0.5 * ( 1 - (2 + phiZInv)*xi + 4*phiZInv*xi3 ) ;
    //Nx[1][2] = 0;
    Nx[1][3] = z * (xi - 1);
    //Nx[1][4] = 0;
    Nx[1][5] = (L / 8) * ( 1 - 2*phiZInv*xi - 4*xi2 + 8*phiZInv*xi3 );
    //Nx[1][6] = 0;
    Nx[1][7] = -Nx[1][1] + 1;
    //Nx[1][8] = 0;
    Nx[1][9] = - z * xi;
    //Nx[1][10] = 0;
    Nx[1][11] = (L / 8) * ( -1 - 2*phiZInv*xi + 4*xi2 + 8*phiZInv*xi3 );

    //Nx[2][0] = 0;
    //Nx[2][1] = 0;
    Nx[2][2] = 0.5 * ( 1 + (phiYInv - 2)*xi - 4*phiYInv*xi3 );
    Nx[2][3] = y * (1 - xi);
    Nx[2][4] = (L / 8) * ( -1 - 2*phiYInv*xi + 4*xi2 + 8*phiYInv*xi3 );
    //Nx[2][5] = 0;
    //Nx[2][6] = 0;
    //Nx[2][7] = 0;
    Nx[2][8] = -Nx[2][2] + 1;
    Nx[2][9] = y * xi;
    Nx[2][10] = (L / 8) * ( 1 - 2*phiYInv*xi - 4*xi2 + 8*phiYInv*xi3);
    //Nx[2][11] = 0;

    return Nx;
}


template <class DataTypes>
auto AdaptivePlasticBeamForceField<DataTypes>::computeGradN(Real x, Real y, Real z, Real L,
                                                            Real A, Real Iy, Real Iz, Real E, Real nu,
                                                            Real kappaY, Real kappaZ) -> Matrix9x12
{
    Matrix9x12 gradN = Matrix9x12(); // Sets each element to 0
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
    Real phiYInv = ( 1 / (phiY - 1) );
    Real phiZInv = ( 1 / (1 + phiZ) );

    //Row 0
    gradN[0][0] = -1 / L;
    gradN[0][1] = - ( 12 * phiZInv * xi * eta ) / L;
    gradN[0][2] = ( 12 * phiYInv * xi * zeta) / L;
    // gradN[0][3] = 0;
    gradN[0][4] = - ( 1 + 6 * phiYInv * xi ) * zeta;
    gradN[0][5] = ( 1 - 6 * phiZInv * xi) * eta;
    gradN[0][6] = 1 / L;
    gradN[0][7] = - gradN[0][1];
    gradN[0][8] = - gradN[0][2];
    // gradN[0][9] = 0;
    gradN[0][10] = gradN[0][4] + 2*zeta;
    gradN[0][11] = gradN[0][5] - 2*eta;

    //Rows 1 and 3
    gradN[1][3] = zeta / 2;
    gradN[1][9] = -zeta / 2;

    gradN[3][3] = zeta / 2;
    gradN[3][9] = -zeta / 2;

    //Rows 2 and 6
    gradN[2][3] = -eta / 2;
    gradN[2][9] = eta / 2;

    gradN[6][3] = -eta / 2;
    gradN[6][9] = eta / 2;

    //Rows 4, 5, 7, 8 are null

    return gradN;
}



template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::reinit()
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
            localForces += (gpWeights[0] * gpWeights[1] * gpWeights[2]) * gp.getGradN().transposed() * newStressPoint;
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
            gp.setMechanicalState(MechanicalState::PLASTIC);

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
        // NB: the change of variables to account for the right integration intervals
        // is done beforehand in the computation of gpWeights and gradN.
        // This include the multiplication of the result by the section dimensions.
        Kt += gpWeights[0]*gpWeights[1]*gpWeights[2] * gradN.transposed()*Cep*gradN;
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


template<class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields()) return;
    if (!this->mstate) return;

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    unsigned int numBeams = l_interpolation->getNumBeams();

    std::vector<defaulttype::Vector3> centrelinePointsSF;
    std::vector<defaulttype::Vector3> centrelinePointsSpline;
    std::vector<defaulttype::Vector3> visualGaussPoints;
    std::vector<RGBAColor> colours;

    for (unsigned int i = 0; i < numBeams; ++i)
    {
        drawElement(i, visualGaussPoints, centrelinePointsSF, centrelinePointsSpline, colours, x);

        ////****** Local frame ******//
        const BeamInterpolation<DataTypes>::BeamSection beamSection = l_interpolation->getBeamSection(i);
        const double L = l_interpolation->getLength(i);

        Vec3 localPos(0.0, 0.0, 0.0);
        Real baryX = 0.5;
        Transform globalHLocalInterpol;
        l_interpolation->InterpolateTransformUsingSpline(i, baryX, localPos, x, globalHLocalInterpol);

        Quat q = globalHLocalInterpol.getOrientation();
        q.normalize();

        Vec3 P1, x, y, z;
        P1 = globalHLocalInterpol.getOrigin();
        x = q.rotate(Vec3(L / 6.0, 0, 0));
        y = q.rotate(Vec3(0, L / 8.0, 0));
        z = q.rotate(Vec3(0, 0, L / 8.0));
        float radius_arrow = (float)L / 60.0f;

        vparams->drawTool()->drawArrow(P1, P1 + x, radius_arrow, sofa::type::RGBAColor(1, 0, 0, 1));
        vparams->drawTool()->drawArrow(P1, P1 + y, radius_arrow, sofa::type::RGBAColor(1, 0, 0, 1));
        vparams->drawTool()->drawArrow(P1, P1 + z, radius_arrow, sofa::type::RGBAColor(1, 0, 0, 1));
    }

    vparams->drawTool()->setPolygonMode(2, true);
    vparams->drawTool()->setLightingEnabled(true);
    vparams->drawTool()->drawPoints(visualGaussPoints, 3, colours);
    vparams->drawTool()->drawLines(centrelinePointsSF, 1.0, RGBAColor(0.24f, 0.72f, 0.96f, 1.0f));
    vparams->drawTool()->drawLines(centrelinePointsSpline, 1.0, RGBAColor(0.16f, 0.61f, 0.07f, 1.0f));
    vparams->drawTool()->setLightingEnabled(false);
    vparams->drawTool()->setPolygonMode(0, false);

    //if (node0Idx == node1Idx) /// rigidification case
    //{
    //    vparams->drawTool()->drawLines(points, 2, sofa::type::RGBAColor(0, 0, 1, 1));
    //    continue;
    //}
    //else
    //    vparams->drawTool()->drawLines(points, 2, sofa::type::RGBAColor(1, 0, 0, 1));
}

template<class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::drawElement(unsigned int beamIndex, std::vector< defaulttype::Vector3 >& visualGaussPoints,
                                                           std::vector< defaulttype::Vector3 >& centrelinePointsSF,
                                                           std::vector< defaulttype::Vector3 >& centrelinePointsSpline,
                                                           std::vector<RGBAColor>& colours,
                                                           const VecCoord& x)
{
    Transform globalH0Local, globalH1Local;
    l_interpolation->computeTransform2(beamIndex, globalH0Local, globalH1Local, x);
    unsigned int node0Idx, node1Idx;
    l_interpolation->getNodeIndices(beamIndex, node0Idx, node1Idx);

    Vec6 U0Local, U1Local;
    computeLocalDisplacement(x, U0Local, U1Local, beamIndex, node0Idx, node1Idx);

    //Displacement in Vec12 form
    Mat<12, 1, Real> discreteU = Mat<12, 1, Real>();
    for (int i = 0; i < 6; i++)
    {
        discreteU[i][0] = U0Local[i];
        discreteU[i + 6][0] = U1Local[i];
    }

    // Getting central point
    Vec3 localPos(0.0, 0.0, 0.0);
    Real baryX = 0.5;
    Transform globalHLocalInterpol;
    l_interpolation->InterpolateTransformUsingSpline(beamIndex, baryX, localPos, x, globalHLocalInterpol);

    Quat q = globalHLocalInterpol.getOrientation();
    q.normalize();
    Vec3 centreNode = globalHLocalInterpol.getOrigin();

    typedef std::function<void(GaussPoint3&)> IntegrationLambda;
    IntegrationLambda computeVisualGP = [&](GaussPoint3& gp)
    {
        Mat<3, 1, Real> continuousU = gp.getNx() * discreteU;
        Vec3 gpCoord = gp.getCoord();

        Vec3 beamVec = { continuousU[0][0] + gpCoord[0],
                         continuousU[1][0] + gpCoord[1],
                         continuousU[2][0] + gpCoord[2] };

        Vec3 visualGP = centreNode + q.rotate(beamVec);
        visualGaussPoints.push_back(visualGP);

        MechanicalState mechanicalState = gp.getMechanicalState();
        if (mechanicalState == MechanicalState::ELASTIC)
            colours.push_back({ 1.0f,0.015f,0.015f,1.0f }); //RED
        else if (mechanicalState == MechanicalState::PLASTIC)
            colours.push_back({ 0.051f,0.15f,0.64f,1.0f }); //BLUE
        else
            colours.push_back({ 0.078f,0.41f,0.078f,1.0f }); //GREEN
    };

    integrateBeam(m_gaussPoints[beamIndex], computeVisualGP);

    ////****** Centreline ******//

    const unsigned int NBSEG = 40; //number of segments descretising the centreline

    //---------- Version with "mechanical" shape functions ----------//

    const BeamInterpolation<DataTypes>::BeamSection beamSection = l_interpolation->getBeamSection(beamIndex);
    double Iy = beamSection._Iy;
    double Iz = beamSection._Iz;
    double A = beamSection._A;
    const double L = l_interpolation->getLength(beamIndex);

    Vec3 node0 = globalH0Local.getOrigin();
    Vec3 node1 = globalH1Local.getOrigin();

    Real E = d_youngModulus.getValue();
    Real nu = d_poissonRatio.getValue();

    Real l0 = m_integrationIntervals[beamIndex].geta1();
    Real l1 = m_integrationIntervals[beamIndex].getb1();
    Real segLength = (l1 - l0) / NBSEG;

    centrelinePointsSF.push_back(node0);

    for (unsigned int seg = 1; seg < NBSEG; seg++)
    {
        //Shape function of the centreline point
        Real x = l0 + seg*segLength;
        Matrix3x12 drawN = computeNx(x, 0, 0, L, A, Iy, Iz, E, nu);
        Mat<3, 1, Real> continuousU = drawN * discreteU;

        Vec3 beamVec = { continuousU[0][0] + x,
                         continuousU[1][0],
                         continuousU[2][0] };
        Vec3 clp = centreNode + q.rotate(beamVec);
        centrelinePointsSF.push_back(clp); //First time as the end of the former segment
        centrelinePointsSF.push_back(clp); //Second time as the beginning of the next segment
    }

    centrelinePointsSF.push_back(node1);

    //---------- Version with spline interpolation functions ----------//

    Vec3 pos = globalH0Local.getOrigin();

    for (double baryCoord = 0.0; baryCoord < 1.0 + 1.0 / NBSEG; baryCoord+= 1.0/NBSEG)
    {
        centrelinePointsSpline.push_back(pos);
        Vec3 localPos(0.0, 0.0, 0.0);
        this->l_interpolation->interpolatePointUsingSpline(beamIndex, baryCoord, localPos, x, pos);
        centrelinePointsSpline.push_back(pos);
    }
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
auto AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::getNx() const -> const Matrix3x12&
{
    return m_Nx;
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::setNx(Matrix3x12 Nx)
{
    m_Nx = Nx;
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
auto AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::getCoord() const -> const Vec3&
{
    return m_coordinates;
}

template <class DataTypes>
void AdaptivePlasticBeamForceField<DataTypes>::GaussPoint3::setCoord(Vec3 coord)
{
    m_coordinates = coord;
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