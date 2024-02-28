/******************************************************************************
*                              BeamAdapter plugin                             *
*                  (c) 2006 Inria, University of Lille, CNRS                  *
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
* Authors: see Authors.md                                                     *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
//
// C++ Implementation : WireBeamInterpolation / AdaptiveBeamForceFieldAndMass
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#pragma once

#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>
#include <sofa/core/behavior/BaseLocalMassMatrix.h>

#include <BeamAdapter/component/forcefield/AdaptiveBeamForceFieldAndMass.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/ScopedAdvancedTimer.h>


namespace sofa::component::forcefield
{

namespace _adaptivebeamforcefieldandmass_
{

/* ************* ADAPTIVE FORCEFIELD_AND_MASS ************** */
using sofa::core::behavior::ForceField ;
using sofa::core::objectmodel::BaseContext ;
using sofa::type::Vec3 ;
using sofa::type::Quat ;
using sofa::helper::ReadAccessor ;
using sofa::core::ConstVecCoordId ;
using std::set ;
using sofa::helper::ScopedAdvancedTimer;

template <class DataTypes>
AdaptiveBeamForceFieldAndMass<DataTypes>::AdaptiveBeamForceFieldAndMass()
    : d_computeMass(initData(&d_computeMass,true,"computeMass","if false, only compute the stiff elastic model"))
    , d_massDensity(initData(&d_massDensity,(Real)1.0,"massDensity", "Density of the mass (usually in kg/m^3)" ))
    , d_useShearStressComputation(initData(&d_useShearStressComputation, true, "shearStressComputation","if false, suppress the shear stress in the computation"))
    , d_reinforceLength(initData(&d_reinforceLength, false, "reinforceLength", "if true, a separate computation for the error in elongation is peformed"))
    , l_interpolation(initLink("interpolation","Path to the Interpolation component on scene"))
    , l_instrumentParameters(initLink("instrumentParameters", "link to an object specifying physical parameters based on abscissa"))
{
}

template <class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Loading);

    if(!l_interpolation)
        l_interpolation.set(dynamic_cast<BaseContext *>(this->getContext())->get<BInterpolation>(BaseContext::Local));

    if (!l_interpolation) {
        msg_error() << "No Beam Interpolation found, the component can not work.";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }

    ForceField<DataTypes>::init();
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}


template <class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::reinit()
{
    init();
}


template <class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::computeGravityVector()
{
    const Vec3& gravity = this->getContext()->getGravity();
    m_gravity = Vec6(gravity[0], gravity[1], gravity[2], 0, 0, 0);
}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::computeStiffness(int beamId, BeamLocalMatrices& beamLocalMatrices)
{
    Real x_curv = 0.0 ;
    Real _nu = 0.0 ;
    Real _E = 0.0 ;

    ///Get the curvilinear abscissa of the extremity of the beam
    l_interpolation->getYoungModulusAtX(beamId, x_curv, _E, _nu);

    /// material parameters
    Real _G = _E / (2.0 * (1.0 + _nu));

    Real phiy, phiz;
    Real L2 = (Real) (beamLocalMatrices._L * beamLocalMatrices._L);
    Real L3 = (Real) (L2 * beamLocalMatrices._L);
    Real EIy = (Real)(_E * beamLocalMatrices._Iy);
    Real EIz = (Real)(_E * beamLocalMatrices._Iz);

    /// Find shear-deformation parameters
    if (beamLocalMatrices._Asy == 0)
        phiy = 0.0;
    else
        phiy =(L2 ==0)? 0.0: (Real)(24.0*(1.0+_nu)* beamLocalMatrices._Iz/(beamLocalMatrices._Asy*L2));

    if (beamLocalMatrices._Asz == 0)
        phiz = 0.0;
    else
        phiz =(L2 ==0)? 0.0: (Real)(24.0*(1.0+_nu)* beamLocalMatrices._Iy/(beamLocalMatrices._Asz*L2));

    beamLocalMatrices.m_K00.clear(); beamLocalMatrices.m_K01.clear(); beamLocalMatrices.m_K10.clear(); beamLocalMatrices.m_K11.clear();

    /// diagonal values
    beamLocalMatrices.m_K00[0][0] = beamLocalMatrices.m_K11[0][0] = (beamLocalMatrices._L == 0.0)? 0.0 :_E* beamLocalMatrices._A/ beamLocalMatrices._L;
    beamLocalMatrices.m_K00[1][1] = beamLocalMatrices.m_K11[1][1] = (L3 == 0.0)? 0.0 :(Real)(12.0*EIz/(L3*(1.0+phiy)));
    beamLocalMatrices.m_K00[2][2] = beamLocalMatrices.m_K11[2][2] = (L3 == 0.0)? 0.0 : (Real)(12.0*EIy/(L3*(1.0+phiz)));
    beamLocalMatrices.m_K00[3][3] = beamLocalMatrices.m_K11[3][3] = (beamLocalMatrices._L == 0.0)? 0.0 : _G* beamLocalMatrices._J/ beamLocalMatrices._L;
    beamLocalMatrices.m_K00[4][4] = beamLocalMatrices.m_K11[4][4] = (beamLocalMatrices._L == 0.0)? 0.0 : (Real)((4.0+phiz)*EIy/(beamLocalMatrices._L*(1.0+phiz)));
    beamLocalMatrices.m_K00[5][5] = beamLocalMatrices.m_K11[5][5] = (beamLocalMatrices._L == 0.0)? 0.0 : (Real)((4.0+phiy)*EIz/(beamLocalMatrices._L*(1.0+phiy)));

    /// diagonal blocks
    beamLocalMatrices.m_K00[4][2]  =(L2 == 0.0)? 0.0 : (Real)(-6.0*EIy/(L2*(1.0+phiz)));
    beamLocalMatrices.m_K00[5][1]  =(L2 == 0.0)? 0.0 : (Real)( 6.0*EIz/(L2*(1.0+phiy)));
    beamLocalMatrices.m_K11[5][1]  = -beamLocalMatrices.m_K00[5][1];
    beamLocalMatrices.m_K11[4][2]  = -beamLocalMatrices.m_K00[4][2];

    /// lower non-diagonal blocks
    beamLocalMatrices.m_K10[0][0]   = -beamLocalMatrices.m_K00[0][0];
    beamLocalMatrices.m_K10[1][1]   = -beamLocalMatrices.m_K00[1][1];
    beamLocalMatrices.m_K10[1][5]   = -beamLocalMatrices.m_K00[5][1];
    beamLocalMatrices.m_K10[2][2]   = -beamLocalMatrices.m_K00[2][2];
    beamLocalMatrices.m_K10[2][4]   = -beamLocalMatrices.m_K00[4][2];
    beamLocalMatrices.m_K10[3][3]   = -beamLocalMatrices.m_K00[3][3];
    beamLocalMatrices.m_K10[4][2]  = beamLocalMatrices.m_K00[4][2];
    beamLocalMatrices.m_K10[4][4]  =(beamLocalMatrices._L == 0.0)? 0.0 : (Real)((2.0-phiz)*EIy/(beamLocalMatrices._L*(1.0+phiz)));
    beamLocalMatrices.m_K10[5][1]  = beamLocalMatrices.m_K00[5][1];
    beamLocalMatrices.m_K10[5][5]  =(beamLocalMatrices._L == 0.0)? 0.0 : (Real)((2.0-phiy)*EIz/(beamLocalMatrices._L*(1.0+phiy)));

    /// Make a symetric matrix with diagonal blocks
    for (int i=0; i<=5; i++)
    {
        for (int j=i+1; j<6; j++)
        {
            beamLocalMatrices.m_K00[i][j] =  beamLocalMatrices.m_K00[j][i];
            beamLocalMatrices.m_K11[i][j] =  beamLocalMatrices.m_K11[j][i];
        }
    }

    /// upper non-diagonal block : set k_loc10 as the transposed matrix of k_loc01
    beamLocalMatrices.m_K01 = beamLocalMatrices.m_K10.transposed();
}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::computeMass(int beamId, BeamLocalMatrices& beamLocalMatrix)
{
    SOFA_UNUSED(beamId);
    Real L2 = (Real) (beamLocalMatrix._L * beamLocalMatrix._L);
    beamLocalMatrix.m_M00.clear(); beamLocalMatrix.m_M01.clear(); beamLocalMatrix.m_M10.clear(); beamLocalMatrix.m_M11.clear();

    Real AL = beamLocalMatrix._A * beamLocalMatrix._L;
    Real Iz_A = (beamLocalMatrix._A == 0.0) ? 0.0 : beamLocalMatrix._Iz / beamLocalMatrix._A;
    Real Iy_A = (beamLocalMatrix._A == 0.0) ? 0.0 : beamLocalMatrix._Iy / beamLocalMatrix._A;

    /// diagonal values
    beamLocalMatrix.m_M00[0][0] = beamLocalMatrix.m_M11[0][0] = (Real)(1.0 / 3.0);
    beamLocalMatrix.m_M00[1][1] = beamLocalMatrix.m_M11[1][1] = (L2 == 0.0) || (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(13.0 / 35.0 + 6.0 * Iz_A / (5.0 * L2));
    beamLocalMatrix.m_M00[2][2] = beamLocalMatrix.m_M11[2][2] = (L2 == 0.0) || (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(13.0 / 35.0 + 6.0 * Iy_A / (5.0 * L2));
    beamLocalMatrix.m_M00[3][3] = beamLocalMatrix.m_M11[3][3] = (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(beamLocalMatrix._J / (3.0 * beamLocalMatrix._A));
    beamLocalMatrix.m_M00[4][4] = beamLocalMatrix.m_M11[4][4] = (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(L2 / 105.0 + 2 * Iy_A / 15.0);
    beamLocalMatrix.m_M00[5][5] = beamLocalMatrix.m_M11[5][5] = (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(L2 / 105.0 + 2 * Iz_A / 15.0);

    /// diagonal blocks
    beamLocalMatrix.m_M00[4][2] = (beamLocalMatrix._L == 0.0) || (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(-11.0 * beamLocalMatrix._L / 210.0 - beamLocalMatrix._Iy / (10 * AL));
    beamLocalMatrix.m_M00[5][1] = (beamLocalMatrix._L == 0.0) || (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(11.0 * beamLocalMatrix._L / 210.0 + beamLocalMatrix._Iz / (10 * AL));
    beamLocalMatrix.m_M11[5][1]  = -beamLocalMatrix.m_M00[5][1];
    beamLocalMatrix.m_M11[4][2]  = -beamLocalMatrix.m_M00[4][2];

    beamLocalMatrix.m_M00 *= beamLocalMatrix._rho * AL;
    beamLocalMatrix.m_M11 *= beamLocalMatrix._rho * AL;

    /// lower non-diagonal blocks
    beamLocalMatrix.m_M10[0][0] = (Real)(1.0 / 6.0);
    beamLocalMatrix.m_M10[1][1] = (L2 == 0.0) || (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(9.0 / 70.0 - 6.0 * Iz_A / (5.0 * L2));
    beamLocalMatrix.m_M10[2][2] = (L2 == 0.0) || (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(9.0 / 70.0 - 6.0 * Iy_A / (5.0 * L2));
    beamLocalMatrix.m_M10[3][3]  = (beamLocalMatrix._A == 0.0) ? 0.0: (Real)(beamLocalMatrix._J/(6.0*beamLocalMatrix._A));
    beamLocalMatrix.m_M10[4][4] = (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(-L2 / 140.0 - Iy_A / 30.0);
    beamLocalMatrix.m_M10[5][5] = (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(-L2 / 140.0 - Iz_A / 30.0);

    beamLocalMatrix.m_M10[1][5] = (beamLocalMatrix._L == 0.0) || (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(13 * beamLocalMatrix._L / 420.0 - beamLocalMatrix._Iz / (10.0 * AL));
    beamLocalMatrix.m_M10[2][4] = (beamLocalMatrix._L == 0.0) || (beamLocalMatrix._A == 0.0) ? 0.0 : (Real)(-13 * beamLocalMatrix._L / 420.0 + beamLocalMatrix._Iy / (10.0 * AL));
    beamLocalMatrix.m_M10[4][2]  = -beamLocalMatrix.m_M10[2][4];
    beamLocalMatrix.m_M10[5][1]  = -beamLocalMatrix.m_M10[1][5];

    beamLocalMatrix.m_M10 *= beamLocalMatrix._rho * AL;

    /// Make a symetric matrix with diagonal blocks
    for (int i=0; i<=5; i++)
    {
        for (int j=i+1; j<6; j++)
        {
            beamLocalMatrix.m_M00[i][j] =  beamLocalMatrix.m_M00[j][i];
            beamLocalMatrix.m_M11[i][j] =  beamLocalMatrix.m_M11[j][i];
        }
    }

    /// upper non-diagonal block : set k_loc10 as the transposed matrix of k_loc01
    beamLocalMatrix.m_M01 = beamLocalMatrix.m_M10.transposed();
}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::applyStiffnessLarge( VecDeriv& df, const VecDeriv& dx,
                                                                    int bIndex, Index nd0Id, Index nd1Id,
                                                                    SReal factor )
{
    if(nd0Id==nd1Id) /// Return in case of rigidification
        return;

    Vec6NoInit U0, U1, u0, u1, f0, f1, F0, F1;
    BeamLocalMatrices &beamLocalMatrix = m_localBeamMatrices[bIndex];

    for (unsigned int i=0; i<6; i++)
    {
        U0[i] = dx[nd0Id][i];
        U1[i] = dx[nd1Id][i];
    }

    /// displacement in local frame
    u0 = beamLocalMatrix.m_A0Ref*U0;
    u1 = beamLocalMatrix.m_A1Ref*U1;

    /// internal force in local frame
    f0 = beamLocalMatrix.m_K00*u0 +  beamLocalMatrix.m_K01*u1;
    f1 = beamLocalMatrix.m_K10*u0 +  beamLocalMatrix.m_K11*u1;

    /// force in global frame
    F0 = beamLocalMatrix.m_A0Ref.multTranspose(f0);
    F1 = beamLocalMatrix.m_A1Ref.multTranspose(f1);

    /// put the result in df
    for (unsigned int i=0; i<6; i++)
    {
        df[nd0Id][i] -= F0[i]*factor;
        df[nd1Id][i] -= F1[i]*factor;
    }
}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::applyMassLarge( VecDeriv& df, int bIndex, Index nd0Id, Index nd1Id, SReal factor)
{
    const BeamLocalMatrices &beamLocalMatrix = m_localBeamMatrices[bIndex];

    /// displacement in local frame (only gravity as external force)
    Vec6 a0 = beamLocalMatrix.m_A0Ref * m_gravity;
    Vec6 a1 = beamLocalMatrix.m_A1Ref * m_gravity;

    /// internal force in local frame
    Vec6 f0 = beamLocalMatrix.m_M00*a0 + beamLocalMatrix.m_M01*a1;
    Vec6 f1 = beamLocalMatrix.m_M10*a0 + beamLocalMatrix.m_M11*a1;

    /// force in global frame
    Vec6 F0 = beamLocalMatrix.m_A0Ref.multTranspose(f0);
    Vec6 F1 = beamLocalMatrix.m_A1Ref.multTranspose(f1);

    /// put the result in df
    for (unsigned int i=0; i<6; i++)
    {
        df[nd0Id][i] += F0[i]*factor;
        df[nd1Id][i] += F1[i]*factor;
    }
}



/////////////////////////////////////
/// Mass Interface
/////////////////////////////////////

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addMDx(const MechanicalParams* mparams , DataVecDeriv& dataf, const DataVecDeriv& datadx, SReal factor)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(datadx);

    auto f = sofa::helper::getWriteOnlyAccessor(dataf);

    auto size = l_interpolation->getStateSize();
    if (f.size() != size)
        f.resize(size);

    unsigned int numBeams = l_interpolation->getNumBeams();
    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        applyMassLarge( f.wref(), b, node0Idx, node1Idx, factor);
    }
}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addMToMatrix(const MechanicalParams *mparams,
                                                            const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(mstate);
    Real mFact = (Real)mparams->mFactor();

    unsigned int numBeams = l_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        BeamLocalMatrices &bLM = m_localBeamMatrices[b];
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        /// matrices in global frame
        Matrix6x6 M00 = bLM.m_A0Ref.multTranspose((bLM.m_M00 * bLM.m_A0Ref));
        Matrix6x6 M01 = bLM.m_A0Ref.multTranspose((bLM.m_M01 * bLM.m_A1Ref));
        Matrix6x6 M10 = bLM.m_A1Ref.multTranspose((bLM.m_M10 * bLM.m_A0Ref));
        Matrix6x6 M11 = bLM.m_A1Ref.multTranspose((bLM.m_M11 * bLM.m_A1Ref));

        int index0[6], index1[6];
        for (int i=0;i<6;i++)
            index0[i] = r.offset+node0Idx*6+i;
        for (int i=0;i<6;i++)
            index1[i] = r.offset+node1Idx*6+i;

        for (int i=0;i<6;i++)
        {
            for (int j=0;j<6;j++)
            {
                r.matrix->add(index0[i], index0[j],  M00(i,j)*mFact);
                r.matrix->add(index0[i], index1[j],  M01(i,j)*mFact);
                r.matrix->add(index1[i], index0[j],  M10(i,j)*mFact);
                r.matrix->add(index1[i], index1[j],  M11(i,j)*mFact);
            }
        }

    }
}

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::buildMassMatrix(sofa::core::behavior::MassMatrixAccumulator* matrices)
{
    const unsigned int numBeams = l_interpolation->getNumBeams();


    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        BeamLocalMatrices &bLM = m_localBeamMatrices[b];
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        /// matrices in global frame
        const Matrix6x6 M00 = bLM.m_A0Ref.multTranspose(bLM.m_M00 * bLM.m_A0Ref);
        const Matrix6x6 M01 = bLM.m_A0Ref.multTranspose(bLM.m_M01 * bLM.m_A1Ref);
        const Matrix6x6 M10 = bLM.m_A1Ref.multTranspose(bLM.m_M10 * bLM.m_A0Ref);
        const Matrix6x6 M11 = bLM.m_A1Ref.multTranspose(bLM.m_M11 * bLM.m_A1Ref);


        matrices->add(node0Idx * 6, node0Idx * 6,  M00);
        matrices->add(node0Idx * 6, node1Idx * 6,  M01);
        matrices->add(node1Idx * 6, node0Idx * 6,  M10);
        matrices->add(node1Idx * 6, node1Idx * 6,  M11);
    }
}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addMBKToMatrix(const MechanicalParams* mparams,
                                                              const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(mstate);
    Real kFact = (Real)mparams->kFactor();
    Real mFact = (Real)mparams->mFactor();

    Real totalMass = 0;
    unsigned int numBeams = l_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        const BeamLocalMatrices &bLM = m_localBeamMatrices[b];
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        int index0[6], index1[6];
        for (int i=0;i<6;i++)
            index0[i] = r.offset+node0Idx*6+i;
        for (int i=0;i<6;i++)
            index1[i] = r.offset+node1Idx*6+i;

        if(node0Idx!=node1Idx) // no rigidification
        {
            // matrices in global frame
            Matrix6x6 K00 = bLM.m_A0Ref.multTranspose((bLM.m_K00 * bLM.m_A0Ref));
            Matrix6x6 K01 = bLM.m_A0Ref.multTranspose((bLM.m_K01 * bLM.m_A1Ref));
            Matrix6x6 K10 = bLM.m_A1Ref.multTranspose((bLM.m_K10 * bLM.m_A0Ref));
            Matrix6x6 K11 = bLM.m_A1Ref.multTranspose((bLM.m_K11 * bLM.m_A1Ref));

            for (int i=0;i<6;i++)
            {
                for (int j=0;j<6;j++)
                {
                    r.matrix->add(index0[i], index0[j], - K00(i,j)*kFact);
                    r.matrix->add(index0[i], index1[j], - K01(i,j)*kFact);
                    r.matrix->add(index1[i], index0[j], - K10(i,j)*kFact);
                    r.matrix->add(index1[i], index1[j], - K11(i,j)*kFact);
                }
            }
        }

        // matrices in global frame
        Matrix6x6 M00 = bLM.m_A0Ref.multTranspose((bLM.m_M00 * bLM.m_A0Ref));
        Matrix6x6 M01 = bLM.m_A0Ref.multTranspose((bLM.m_M01 * bLM.m_A1Ref));
        Matrix6x6 M10 = bLM.m_A1Ref.multTranspose((bLM.m_M10 * bLM.m_A0Ref));
        Matrix6x6 M11 = bLM.m_A1Ref.multTranspose((bLM.m_M11 * bLM.m_A1Ref));

        for (int i=0;i<6;i++)
        {
            for (int j=0;j<6;j++)
            {
                totalMass += M00(i,j)*mFact + M11(i,j)*mFact;
                r.matrix->add(index0[i], index0[j],  M00(i,j)*mFact);
                r.matrix->add(index0[i], index1[j],  M01(i,j)*mFact);
                r.matrix->add(index1[i], index0[j],  M10(i,j)*mFact);
                r.matrix->add(index1[i], index1[j],  M11(i,j)*mFact);
            }
        }
    }
}

template <class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::buildDampingMatrix(core::behavior::DampingMatrix*)
{
    // No damping in this ForceField
}


/////////////////////////////////////
/// ForceField Interface
/////////////////////////////////////

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addForce (const MechanicalParams* mparams ,
                                                         DataVecDeriv& dataf,
                                                         const DataVecCoord& datax,
                                                         const DataVecDeriv& v)
{
    SCOPED_TIMER("AdaptiveBeamForceFieldAndMass_addForce");
    SOFA_UNUSED(v);

    auto f = sofa::helper::getWriteOnlyAccessor(dataf);
    const VecCoord& x = datax.getValue();

    f.resize(x.size()); // current content of the vector will remain the same (http://www.cplusplus.com/reference/vector/vector/resize/)

    unsigned int numBeams = l_interpolation->getNumBeams();
    m_localBeamMatrices.resize(numBeams);

    if(d_computeMass.getValue())
        computeGravityVector();

    /// TODO:
    ///* Redimentionner _localBeamMatrices
    ///* Calculer les rotation et les transformations
    ///* Calculer la matrice "locale"
    ///* Calculer la force exercée par chaque beam
    ///* Calculer la force exercée par la gravitée
    for (unsigned int beamId=0; beamId <numBeams; beamId++)
    {
        ///find the indices of the nodes
        sofa::Index node0Idx, node1Idx;
        l_interpolation->getNodeIndices(beamId, node0Idx, node1Idx);

        ///find the beamMatrices:
        BeamLocalMatrices& beamMatrices = m_localBeamMatrices[beamId];

        ///////////// new : Calcul du repère local de la beam & des transformations adequates///////////////
        Transform global_H_local0, global_H_local1;

        /// 1. get the current transform of the beam:
        l_interpolation->computeTransform(beamId, node0Idx, node1Idx, global_H_local0, global_H_local1, x);

        /// 2. Computes the frame of the beam based on the spline interpolation:
        Transform global_H_local;
        Real baryX = 0.5;
        Real L = l_interpolation->getLength(beamId);

        l_interpolation->InterpolateTransformUsingSpline(global_H_local, baryX, global_H_local0, global_H_local1, L);

        Transform local_H_local0 = global_H_local.inversed()*global_H_local0;
        Transform local_H_local1 = global_H_local.inversed()*global_H_local1;

        /// 3. Computes the transformation from the DOF (in global frame) to the node's local frame DOF0global_H_Node0local and DOF1global_H_Node1local
        Transform DOF0_H_local0, DOF1_H_local1;
        l_interpolation->getDOFtoLocalTransform(beamId, DOF0_H_local0, DOF1_H_local1);

        /// 4. Computes the adequate transformation
        Transform global_R_DOF0(Vec3(0,0,0), x[node0Idx].getOrientation());
        Transform global_R_DOF1(Vec3(0,0,0), x[node1Idx].getOrientation());
        /// - rotation due to the optional transformation
        global_H_local0 = global_R_DOF0*DOF0_H_local0;
        global_H_local1 = global_R_DOF1*DOF1_H_local1;

        Transform DOF0global_H_Node0local, DOF1global_H_Node1local;

        DOF0global_H_Node0local.set(global_H_local0.getOrigin(), global_H_local.getOrientation() );
        DOF1global_H_Node1local.set(global_H_local1.getOrigin(), global_H_local.getOrientation() );

        //TODO(dmarchal 2017-05-17) Please specify who/when this will be done
        //TODO A verifier : global_H_local0.getOrigin() == x[node0Idx].getOrientation().rotate(DOF0_H_local0.getOrigin())

        /// compute Adjoint Matrices:
        beamMatrices.m_A0Ref = DOF0global_H_Node0local.inversed().getAdjointMatrix();
        beamMatrices.m_A1Ref = DOF1global_H_Node1local.inversed().getAdjointMatrix();

        /////////////////////////////////////// COMPUTATION OF THE MASS AND STIFFNESS  MATRIX (LOCAL)

        /// Update Interpolation & geometrical parameters with current positions

        /// material parameters
        beamMatrices._rho = d_massDensity.getValue();

        /// Temp : we only overide values for which a Data has been set in the WireRestShape
        if (l_instrumentParameters.get())
        {
            Real x_curv = 0;
            l_interpolation->getAbsCurvXFromBeam(beamId, x_curv);

            /// The length of the beams is only known to the interpolation !
            l_instrumentParameters->getInterpolationParam(x_curv, beamMatrices._rho, beamMatrices._A, beamMatrices._Iy,
                beamMatrices._Iz, beamMatrices._Asy, beamMatrices._Asz, beamMatrices._J);
        }
        else
        {
            l_interpolation->getInterpolationParam(beamId, beamMatrices._L, beamMatrices._A, beamMatrices._Iy,
                beamMatrices._Iz, beamMatrices._Asy, beamMatrices._Asz, beamMatrices._J);
        }


        /// compute the local mass matrices
        if(d_computeMass.getValue())
        {
            computeMass(beamId, beamMatrices);
        }

        /// IF RIGIDIFICATION: no stiffness forces:
        if(node0Idx==node1Idx)
            continue;

        /// compute the local stiffness matrices
        computeStiffness(beamId, beamMatrices);

        /////////////////////////////COMPUTATION OF THE STIFFNESS FORCE
        /// compute the current local displacement of the beam (6dofs)
        /// 1. get the rest transformation from local to 0 and local to 1
        Transform local_H_local0_rest,local_H_local1_rest;
        l_interpolation->getSplineRestTransform(beamId, local_H_local0_rest, local_H_local1_rest);

        ///2. computes the local displacement of 0 and 1 in frame local:
        SpatialVector u0 = local_H_local0.CreateSpatialVector() - local_H_local0_rest.CreateSpatialVector();
        SpatialVector u1 = local_H_local1.CreateSpatialVector() - local_H_local1_rest.CreateSpatialVector();

        /// 3. put the result in a Vec6
        Vec6NoInit U0local, U1local;

        for (unsigned int i=0; i<3; i++)
        {
            U0local[i] = u0.getLinearVelocity()[i];
            U0local[i+3] = u0.getAngularVelocity()[i];
            U1local[i] = u1.getLinearVelocity()[i];
            U1local[i+3] = u1.getAngularVelocity()[i];
        }

        if(d_reinforceLength.getValue())
        {
            Vec3 P0,P1,P2,P3;
            Real length;
            Real rest_length = l_interpolation->getLength(beamId);
            l_interpolation->getSplinePoints(beamId, x, P0, P1, P2, P3);
            l_interpolation->computeActualLength(length, P0, P1, P2, P3);

            U0local[0]=(-length+rest_length)/2;
            U1local[0]=( length-rest_length)/2;
        }

        if (!d_useShearStressComputation.getValue())
        {
            /////////////////// TEST //////////////////////
            /// test: correction due to spline computation;
            Vec3 ResultNode0, ResultNode1;
            l_interpolation->computeStrechAndTwist(beamId, x, ResultNode0, ResultNode1);

            Real ux0 =-ResultNode0[0] + l_interpolation->getLength(beamId)/2;
            Real ux1 = ResultNode1[0] - l_interpolation->getLength(beamId)/2;

            U0local[0] = ux0;
            U1local[0] = ux1;

            U0local[3] =-ResultNode0[2];
            U1local[3] = ResultNode1[2];

            //////////////////////////////////////////////////
        }

        /// compute the force in the local frame:
        Vec6 f0 = beamMatrices.m_K00 * U0local + beamMatrices.m_K01 * U1local;
        Vec6 f1 = beamMatrices.m_K10 * U0local + beamMatrices.m_K11 * U1local;

        /// compute the force in the global frame
        Vec6 F0_ref = beamMatrices.m_A0Ref.multTranspose(f0);
        Vec6 F1_ref = beamMatrices.m_A1Ref.multTranspose(f1);

        /// Add this force to vector f
        for (unsigned int i=0; i<6; i++)
        {
            f[node0Idx][i]-=F0_ref[i];
            f[node1Idx][i]-=F1_ref[i];
        }
    }

    if(d_computeMass.getValue())
    {
        /// will add gravity directly using m_gravity:
        DataVecDeriv emptyVec;
        addMDx(mparams, dataf, emptyVec, 1.0);
    }
}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addDForce(const MechanicalParams* mparams,
                                                         DataVecDeriv& datadF, const DataVecDeriv& datadX )
{
    auto df = sofa::helper::getWriteOnlyAccessor(datadF);
    const VecDeriv& dx = datadX.getValue();
    const double kFactor = mparams->kFactor();

    df.resize(dx.size()); // current content of the vector will remain the same (http://www.cplusplus.com/reference/vector/vector/resize/)

    unsigned int numBeams = l_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        applyStiffnessLarge( df.wref(), dx, b, node0Idx, node1Idx, kFactor);
    }
}



template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addKToMatrix(const MechanicalParams* mparams,
                                                            const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef matrixRef = matrix->getMatrix(mstate);
    Real k = (Real)mparams->kFactor();

    unsigned int numBeams = l_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        const BeamLocalMatrices &beamLocalMatrix = m_localBeamMatrices[b];
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        if(node0Idx==node1Idx)
            continue;

        // matrices in global frame
        Matrix6x6 K00 = beamLocalMatrix.m_A0Ref.multTranspose((beamLocalMatrix.m_K00 * beamLocalMatrix.m_A0Ref));
        Matrix6x6 K01 = beamLocalMatrix.m_A0Ref.multTranspose((beamLocalMatrix.m_K01 * beamLocalMatrix.m_A1Ref));
        Matrix6x6 K10 = beamLocalMatrix.m_A1Ref.multTranspose((beamLocalMatrix.m_K10 * beamLocalMatrix.m_A0Ref));
        Matrix6x6 K11 = beamLocalMatrix.m_A1Ref.multTranspose((beamLocalMatrix.m_K11 * beamLocalMatrix.m_A1Ref));

        int index0[6], index1[6];
        for (int i=0;i<6;i++)
            index0[i] = matrixRef.offset+node0Idx*6+i;
        for (int i=0;i<6;i++)
            index1[i] = matrixRef.offset+node1Idx*6+i;

        for (int i=0;i<6;i++)
        {
            for (int j=0;j<6;j++)
            {
                matrixRef.matrix->add(index0[i], index0[j], - K00(i,j)*k);
                matrixRef.matrix->add(index0[i], index1[j], - K01(i,j)*k);
                matrixRef.matrix->add(index1[i], index0[j], - K10(i,j)*k);
                matrixRef.matrix->add(index1[i], index1[j], - K11(i,j)*k);
            }
        }
    }
}

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::buildStiffnessMatrix(core::behavior::StiffnessMatrix* matrix)
{
    const unsigned int numBeams = l_interpolation->getNumBeams();


    auto dfdx = matrix->getForceDerivativeIn(this->mstate)
                       .withRespectToPositionsIn(this->mstate);

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        const BeamLocalMatrices &beamLocalMatrix = m_localBeamMatrices[b];
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        if(node0Idx==node1Idx)
            continue;

        // matrices in global frame
        Matrix6x6 K00 = beamLocalMatrix.m_A0Ref.multTranspose(beamLocalMatrix.m_K00 * beamLocalMatrix.m_A0Ref);
        Matrix6x6 K01 = beamLocalMatrix.m_A0Ref.multTranspose(beamLocalMatrix.m_K01 * beamLocalMatrix.m_A1Ref);
        Matrix6x6 K10 = beamLocalMatrix.m_A1Ref.multTranspose(beamLocalMatrix.m_K10 * beamLocalMatrix.m_A0Ref);
        Matrix6x6 K11 = beamLocalMatrix.m_A1Ref.multTranspose(beamLocalMatrix.m_K11 * beamLocalMatrix.m_A1Ref);


        dfdx(node0Idx*6, node0Idx*6) += - K00;
        dfdx(node0Idx*6, node1Idx*6) += - K01;
        dfdx(node1Idx*6, node0Idx*6) += - K10;
        dfdx(node1Idx*6, node1Idx*6) += - K11;
    }
}


/////////////////////////////////////
/// Visualization
/////////////////////////////////////

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::draw(const VisualParams *vparams)
{
    if (!vparams->displayFlags().getShowForceFields() && !vparams->displayFlags().getShowBehaviorModels()) return;
    if (!mstate) return;

    vparams->drawTool()->saveLastState();

    ReadAccessor<Data<VecCoord> > x = mstate->read(ConstVecCoordId::position()) ;

    unsigned int numBeams = l_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        Transform globalH0Local,  globalH1Local;

        unsigned int node0Idx, node1Idx;
        l_interpolation->getNodeIndices(b, node0Idx, node1Idx);
        l_interpolation->computeTransform(b, node0Idx, node1Idx, globalH0Local, globalH1Local, x.ref());

        if (vparams->displayFlags().getShowBehaviorModels() && node0Idx!=node1Idx)
            drawElement(vparams, b, globalH0Local, globalH1Local);

        if(vparams->displayFlags().getShowForceFields())
        {
            // /  test ///
            type::vector<type::Vec3> points;
            Vec3 pos = globalH0Local.getOrigin();
            for (double i=0.0; i<1.00001; i+=0.02)
            {
                points.push_back(pos);
                Vec3 localPos(0.0,0.0,0.0);
                this->l_interpolation->interpolatePointUsingSpline(b, i, localPos, x.ref(), pos);
                points.push_back(pos);
            }

            if(node0Idx==node1Idx) /// rigidification case
            {
                vparams->drawTool()->drawLines(points,2, sofa::type::RGBAColor(0,0,1,1));
                continue;
            }
            else
                vparams->drawTool()->drawLines(points,2, sofa::type::RGBAColor(1,0,0,1));

            double length = (double) l_interpolation->getLength(b);


            Vec3 localPos(0.0,0.0,0.0);
            Real baryX = 0.5;
            Transform globalHLocalInterpol;
            l_interpolation->InterpolateTransformUsingSpline(b, baryX, localPos, x.ref(), globalHLocalInterpol);

            Quat q = globalHLocalInterpol.getOrientation();
            q.normalize();

            Vec3 P1, x,y,z;
            P1 = globalHLocalInterpol.getOrigin();
            x= q.rotate(Vec3(length/6.0,0,0));
            y= q.rotate(Vec3(0,length/8.0,0));
            z= q.rotate(Vec3(0,0,length/8.0));
            float radius_arrow = (float)length/60.0f;

            vparams->drawTool()->drawArrow(P1,P1 + x, radius_arrow, sofa::type::RGBAColor(1,0,0,1));
            vparams->drawTool()->drawArrow(P1,P1 + y, radius_arrow, sofa::type::RGBAColor(1,0,0,1));
            vparams->drawTool()->drawArrow(P1,P1 + z, radius_arrow, sofa::type::RGBAColor(1,0,0,1));
        }
    }

    vparams->drawTool()->restoreLastState();
}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::drawElement(const VisualParams *vparams, int beam,
                                                           Transform &global_H0_local, Transform &global_H1_local)
{
    double length = (double) l_interpolation->getLength(beam);

    /// ARROWS
    Vec3 sizeArrows (length/4., length/8., length/8.);

    vparams->drawTool()->drawFrame(global_H0_local.getOrigin(), global_H0_local.getOrientation(), sizeArrows);
    vparams->drawTool()->drawFrame(global_H1_local.getOrigin(), global_H1_local.getOrientation(), sizeArrows);
}


} /// namespace _adaptivebeamforcefieldandmass_

} /// namespace sofa::component::forcefield
