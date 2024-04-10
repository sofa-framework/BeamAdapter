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
// C++ Implementation : WireBeamInterpolation / AdaptiveInflatableBeamForceField
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

#include <sofa/core/behavior/ForceField.inl>
#include <BeamAdapter/component/forcefield/AdaptiveInflatableBeamForceField.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/helper/decompose.h>


#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/helper/OptionsGroup.h>

#include <sofa/simulation/Simulation.h>
#include <sofa/core/visual/VisualParams.h>


namespace sofa::component::forcefield::_AdaptiveInflatableBeamForceField_
{

/* ************* ADAPTIVE FORCEFIELD_AND_MASS ************** */
using sofa::core::behavior::ForceField ;
using sofa::core::objectmodel::BaseContext ;
using sofa::type::Vec3 ;
using sofa::type::Quat ;
using sofa::helper::ReadAccessor ;
using sofa::core::ConstVecCoordId ;
using std::set ;

template <class DataTypes>
AdaptiveInflatableBeamForceField<DataTypes>::AdaptiveInflatableBeamForceField()
    : d_computeMass(initData(&d_computeMass,true,"computeMass","if false, only compute the stiff elastic model"))
    , d_massDensity(initData(&d_massDensity,(Real)1.0,"massDensity", "Density of the mass (usually in kg/m^3)" ))
    , d_reinforceLength(initData(&d_reinforceLength, false, "reinforceLength", "if true, a separate computation for the error in elongation is peformed"))
    , d_pressure(initData(&d_pressure, (Real)0.0, "pressure", "pressure inside the inflatable Beam"))
    , d_dataG(initData(&d_dataG,"dataG","Gravity 3d vector"))

    , l_interpolation(initLink("interpolation","Path to the Interpolation component on scene"))
    , l_instrumentParameters(initLink("instrumentParameters", "link to an object specifying physical parameters based on abscissa"))
{
}

template <class DataTypes>
void AdaptiveInflatableBeamForceField<DataTypes>::init()
{
    if(!l_interpolation)
        l_interpolation.set(dynamic_cast<BaseContext *>(this->getContext())->get<BInterpolation>(BaseContext::Local));

    if(!l_interpolation)
        msg_error() << "No Beam Interpolation found, the component can not work.";

    ForceField<DataTypes>::init();
}


template <class DataTypes>
void AdaptiveInflatableBeamForceField<DataTypes>::reinit()
{
    init();
}


template <class DataTypes>
void AdaptiveInflatableBeamForceField<DataTypes>::computeGravityVector()
{
    Vec3 gravity = getContext()->getGravity();

    auto _G = sofa::helper::getWriteOnlyAccessor(d_dataG);
    _G.resize(l_interpolation->getStateSize());

    m_gravity = Vec3(gravity[0],gravity[1],gravity[2]);

    if(_G.size()==0)
    {
        dmsg_warning() << "_G.size = 0";
        return;
    }

    for (unsigned int i=0; i<_G.size(); i++)
    {
        _G[i][0]=gravity[0];
        _G[i][1]=gravity[1];
        _G[i][2]=gravity[2];

        _G[i][3]=(Real)0.0;
        _G[i][4]=(Real)0.0;
        _G[i][5]=(Real)0.0;
    }
}


template<class DataTypes>
void AdaptiveInflatableBeamForceField<DataTypes>::computeStiffness(int beam, BeamLocalMatrices& beamLocalMatrices)
{
    Real x_curv = 0.0 ;
    Real _nu = 0.0 ;
    Real _E = 0.0 ;

    ///Get the curvilinear abscissa of the extremity of the beam
    l_interpolation->getYoungModulusAtX(beam,x_curv, _E, _nu);

    /// material parameters
    Real _G;
    _G=_E/(2.0*(1.0+_nu));

    /// interpolation & geometrical parameters
    Real _A, _L, _Iy, _Iz, _Asy, _Asz, _J;
    l_interpolation->getInterpolationParam(beam, _L, _A, _Iy , _Iz, _Asy, _Asz, _J);



    /// Temp : we only overide values for which a Data has been set in the WireRestShape
    if(l_instrumentParameters.get())
    {
        Real x_curv = 0, _rho;
        l_interpolation->getAbsCurvXFromBeam(beam, x_curv);
        l_instrumentParameters->getInterpolationParam(x_curv, _rho, _A, _Iy , _Iz, _Asy, _Asz, _J);	// The length of the beams is only known to the interpolation !
    }

    /// correction for inflated beam (effective shear area for circular tubes with thin walls
    _Asy = 0.5;
    _Asz = 0.5;

    /// pressure
    Real P = this->d_pressure.getValue();



    Real phiy, phiz;
    Real L2 = (Real) (_L * _L);
    Real L3 = (Real) (L2 * _L);
    Real E_P_Iy = (Real)((_E+P) * _Iy);
    Real E_P_Iz = (Real)((_E+P) * _Iz);

    /// Find shear-deformation parameters
    if (_Asy == 0)
        phiy = 0.0;
    else
        phiy = (Real)12*((_E+P)*_Iz)/((P+_Asy*_G)*_A*L2); /// formula from "Mechanics of Inflatable Fabric Beams, C. Wielgosz et al.
                                                            /// except that they consider the Pressure as a force (P*_A)

    if (_Asz == 0)
        phiz = 0.0;
    else
        phiz = (Real)12*((_E+P)*_Iy)/((P+_Asy*_G)*_A*L2);

    beamLocalMatrices.m_K00.clear(); beamLocalMatrices.m_K01.clear(); beamLocalMatrices.m_K10.clear(); beamLocalMatrices.m_K11.clear();



    /// diagonal values
    beamLocalMatrices.m_K00[0][0] = beamLocalMatrices.m_K11[0][0] = _E*_A/_L;
    beamLocalMatrices.m_K00[1][1] = beamLocalMatrices.m_K11[1][1] = (Real)(12.0*E_P_Iz/(L3*(1.0+phiy)));
    beamLocalMatrices.m_K00[2][2] = beamLocalMatrices.m_K11[2][2]   = (Real)(12.0*E_P_Iy/(L3*(1.0+phiz)));
    beamLocalMatrices.m_K00[3][3] = beamLocalMatrices.m_K11[3][3]   = _G*_J/_L;
    beamLocalMatrices.m_K00[4][4] = beamLocalMatrices.m_K11[4][4]   = (Real)((4.0+phiz)*E_P_Iy/(_L*(1.0+phiz)));
    beamLocalMatrices.m_K00[5][5] = beamLocalMatrices.m_K11[5][5]   = (Real)((4.0+phiy)*E_P_Iz/(_L*(1.0+phiy)));

    /// diagonal blocks
    beamLocalMatrices.m_K00[4][2]   = (Real)(-6.0*E_P_Iy/(L2*(1.0+phiz)));
    beamLocalMatrices.m_K00[5][1]   = (Real)( 6.0*E_P_Iz/(L2*(1.0+phiy)));
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
    beamLocalMatrices.m_K10[4][4]  = (Real)((2.0-phiz)*E_P_Iy/(_L*(1.0+phiz)));
    beamLocalMatrices.m_K10[5][1]  = beamLocalMatrices.m_K00[5][1];
    beamLocalMatrices.m_K10[5][5]  = (Real)((2.0-phiy)*E_P_Iz/(_L*(1.0+phiy)));

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
void AdaptiveInflatableBeamForceField<DataTypes>::computeMass(int beam, BeamLocalMatrices& beamLocalMatrix)
{
    /// material parameters
    Real _rho;
    _rho = d_massDensity.getValue();

    /// interpolation & geometrical parameters
    Real _A, _L, _Iy, _Iz, _Asy, _Asz, _J;
    l_interpolation->getInterpolationParam(beam, _L, _A, _Iy , _Iz, _Asy, _Asz, _J);

    /// Temp : we only overide values for which a Data has been set in the WireRestShape
    if(l_instrumentParameters.get())
    {
        Real x_curv = 0;
        l_interpolation->getAbsCurvXFromBeam(beam, x_curv);

        /// The length of the beams is only known to the interpolation !
        l_instrumentParameters->getInterpolationParam(x_curv, _rho, _A, _Iy , _Iz, _Asy, _Asz, _J);
    }

    Real L2 = (Real) (_L * _L);
    beamLocalMatrix.m_M00.clear(); beamLocalMatrix.m_M01.clear(); beamLocalMatrix.m_M10.clear(); beamLocalMatrix.m_M11.clear();

    /// diagonal values
    beamLocalMatrix.m_M00[0][0] = beamLocalMatrix.m_M11[0][0] = (Real)(1.0/3.0);
    beamLocalMatrix.m_M00[1][1] = beamLocalMatrix.m_M11[1][1] = (Real)(13.0/35.0 + 6.0*_Iz/(5.0*_A*L2));
    beamLocalMatrix.m_M00[2][2] = beamLocalMatrix.m_M11[2][2] = (Real)(13.0/35.0 + 6.0*_Iy/(5.0*_A*L2));
    beamLocalMatrix.m_M00[3][3] = beamLocalMatrix.m_M11[3][3] = (Real)(_J/(3.0*_A));
    beamLocalMatrix.m_M00[4][4] = beamLocalMatrix.m_M11[4][4] = (Real)(L2/105.0 + 2*_Iy/(15.0*_A));
    beamLocalMatrix.m_M00[5][5] = beamLocalMatrix.m_M11[5][5] = (Real)(L2/105.0 + 2*_Iz/(15.0*_A));

    /// diagonal blocks
    beamLocalMatrix.m_M00[4][2]  = (Real)(-11.0*_L/210.0 - _Iy/(10*_A*_L)  );
    beamLocalMatrix.m_M00[5][1]  = (Real)( 11.0*_L/210.0 + _Iz/(10*_A*_L)  );
    beamLocalMatrix.m_M11[5][1]  = -beamLocalMatrix.m_M00[5][1];
    beamLocalMatrix.m_M11[4][2]  = -beamLocalMatrix.m_M00[4][2];

    beamLocalMatrix.m_M00 *= _rho*_A*_L;
    beamLocalMatrix.m_M11 *= _rho*_A*_L;

    /// lower non-diagonal blocks
    beamLocalMatrix.m_M10[0][0]  = (Real)(1.0/6.0);
    beamLocalMatrix.m_M10[1][1]  = (Real)(9.0/70.0 - 6.0*_Iz/(5.0*_A*L2));
    beamLocalMatrix.m_M10[2][2]  = (Real)(9.0/70.0 - 6.0*_Iy/(5.0*_A*L2));
    beamLocalMatrix.m_M10[3][3]  = (Real)(_J/(6.0*_A));
    beamLocalMatrix.m_M10[4][4]  = (Real)(-L2/140.0 - _Iy/(30.0*_A));
    beamLocalMatrix.m_M10[5][5]  = (Real)(-L2/140.0 - _Iz/(30.0*_A));

    beamLocalMatrix.m_M10[1][5]  = (Real)( 13*_L/420.0 - _Iz/(10.0*_A*_L));
    beamLocalMatrix.m_M10[2][4]  = (Real)(-13*_L/420.0 + _Iy/(10.0*_A*_L));
    beamLocalMatrix.m_M10[4][2]  = -beamLocalMatrix.m_M10[2][4];
    beamLocalMatrix.m_M10[5][1]  = -beamLocalMatrix.m_M10[1][5];

    beamLocalMatrix.m_M10 *= _rho*_A*_L;

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
void AdaptiveInflatableBeamForceField<DataTypes>::applyStiffnessLarge( VecDeriv& df, const VecDeriv& dx,
                                                                    int bIndex, Index nd0Id, Index nd1Id,
                                                                    SReal factor )
{
    if(nd0Id==nd1Id) /// Return in case of rigidification
        return;

    Vec6 U0, U1, u0, u1, f0, f1, F0, F1;
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
void AdaptiveInflatableBeamForceField<DataTypes>::applyMassLarge( VecDeriv& df, const VecDeriv& dx,
                                                               int bIndex, Index nd0Id, Index nd1Id,
                                                               SReal factor)
{
    Vec6 A0, A1, a0, a1, f0, f1, F0, F1;
    BeamLocalMatrices &beamLocalMatrix = m_localBeamMatrices[bIndex];

    for (unsigned int i=0; i<6; i++)
    {
        A0[i] = dx[nd0Id][i];
        A1[i] = dx[nd1Id][i];
    }

    /// displacement in local frame
    a0 = beamLocalMatrix.m_A0Ref*A0;
    a1 = beamLocalMatrix.m_A1Ref*A1;

    /// internal force in local frame
    f0 = beamLocalMatrix.m_M00*a0 + beamLocalMatrix.m_M01*a1;
    f1 = beamLocalMatrix.m_M10*a0 + beamLocalMatrix.m_M11*a1;

    /// force in global frame
    F0 = beamLocalMatrix.m_A0Ref.multTranspose(f0);
    F1 = beamLocalMatrix.m_A1Ref.multTranspose(f1);

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
void AdaptiveInflatableBeamForceField<DataTypes>::addMDx(const MechanicalParams* mparams , DataVecDeriv& dataf, const DataVecDeriv& datadx, SReal factor)
{
    SOFA_UNUSED(mparams);

    auto f = sofa::helper::getWriteOnlyAccessor(dataf);
    const VecDeriv& dx = datadx.getValue();

    unsigned int numBeams = l_interpolation->getNumBeams();

    if (f.size()!=dx.size())
        f.resize(dx.size()); // current content of the vector will remain the same (http://www.cplusplus.com/reference/vector/vector/resize/)

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        applyMassLarge( f.wref(), dx, b, node0Idx, node1Idx, factor);
    }
}


template<class DataTypes>
void AdaptiveInflatableBeamForceField<DataTypes>::addMToMatrix(const MechanicalParams *mparams,
                                                            const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    Real mFact = (Real)mparams->mFactor();

    unsigned int numBeams = l_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        BeamLocalMatrices &bLM = m_localBeamMatrices[b];
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        /// matrices in global frame
        Matrix6x6 M00, M01, M10, M11;

        M00=bLM.m_A0Ref.multTranspose( ( bLM.m_M00 * bLM.m_A0Ref  )  );
        M01=bLM.m_A0Ref.multTranspose( ( bLM.m_M01 * bLM.m_A1Ref  )  );
        M10=bLM.m_A1Ref.multTranspose( ( bLM.m_M10 * bLM.m_A0Ref  )  );
        M11=bLM.m_A1Ref.multTranspose( ( bLM.m_M11 * bLM.m_A1Ref  )  );

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
void AdaptiveInflatableBeamForceField<DataTypes>::addMBKToMatrix(const MechanicalParams* mparams,
                                                              const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    Real kFact = (Real)mparams->kFactor();
    Real mFact = (Real)mparams->mFactor();

    Real totalMass = 0;
    unsigned int numBeams = l_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        BeamLocalMatrices &bLM = m_localBeamMatrices[b];
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        int index0[6], index1[6];
        for (int i=0;i<6;i++)
            index0[i] = r.offset+node0Idx*6+i;
        for (int i=0;i<6;i++)
            index1[i] = r.offset+node1Idx*6+i;

        if(node0Idx!=node1Idx) // no rigidification
        {
            // matrices in global frame
            Matrix6x6 K00, K01, K10, K11;

            K00=bLM.m_A0Ref.multTranspose( ( bLM.m_K00 * bLM.m_A0Ref  )  );
            K01=bLM.m_A0Ref.multTranspose( ( bLM.m_K01 * bLM.m_A1Ref  )  );
            K10=bLM.m_A1Ref.multTranspose( ( bLM.m_K10 * bLM.m_A0Ref  )  );
            K11=bLM.m_A1Ref.multTranspose( ( bLM.m_K11 * bLM.m_A1Ref  )  );

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
        Matrix6x6 M00, M01, M10, M11;

        M00=bLM.m_A0Ref.multTranspose( ( bLM.m_M00 * bLM.m_A0Ref  )  );
        M01=bLM.m_A0Ref.multTranspose( ( bLM.m_M01 * bLM.m_A1Ref  )  );
        M10=bLM.m_A1Ref.multTranspose( ( bLM.m_M10 * bLM.m_A0Ref  )  );
        M11=bLM.m_A1Ref.multTranspose( ( bLM.m_M11 * bLM.m_A1Ref  )  );

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


/////////////////////////////////////
/// ForceField Interface
/////////////////////////////////////

template<class DataTypes>
void AdaptiveInflatableBeamForceField<DataTypes>::addForce (const MechanicalParams* mparams ,
                                                         DataVecDeriv& dataf,
                                                         const DataVecCoord& datax,
                                                         const DataVecDeriv& v)
{
    SOFA_UNUSED(v);

    auto f = sofa::helper::getWriteOnlyAccessor(dataf);
    const VecCoord& x = datax.getValue();

    f.resize(x.size()); // current content of the vector will remain the same (http://www.cplusplus.com/reference/vector/vector/resize/)

    unsigned int numBeams = l_interpolation->getNumBeams();
    m_localBeamMatrices.resize(numBeams);

    if(d_computeMass.getValue())
        computeGravityVector();

    /// TODO:
    ///* Resize _localBeamMatrices

    ///* Compute rotations and transformations

    ///* Compute local frame matrix

    ///* Compute force contribution of each beam

    ///* Calculer la force exercée par la gravitée
    for (unsigned int b=0; b<numBeams; b++)
    {
        ///find the indices of the nodes
        sofa::Index node0Idx, node1Idx;
        l_interpolation->getNodeIndices(b, node0Idx, node1Idx);

        ///find the beamMatrices:
        BeamLocalMatrices  *beamMatrices = &m_localBeamMatrices[b]  ;//new BeamLocalMatrices();

        ///////////// new : Calcul du repère local de la beam & des transformations adequates///////////////
        Transform global_H_local0, global_H_local1;

        /// 1. get the current transform of the beam:
        dmsg_info() << "in addForce";
        l_interpolation->computeTransform(b, node0Idx, node1Idx, global_H_local0, global_H_local1, x);

        /// 2. Computes the frame of the beam based on the spline interpolation:
        Transform global_H_local;
        Real baryX = 0.5;
        Real L = l_interpolation->getLength(b);

        l_interpolation->InterpolateTransformUsingSpline(global_H_local, baryX, global_H_local0, global_H_local1, L);

        Transform local_H_local0 = global_H_local.inversed()*global_H_local0;
        Transform local_H_local1 = global_H_local.inversed()*global_H_local1;

        /// 3. Computes the transformation from the DOF (in global frame) to the node's local frame DOF0global_H_Node0local and DOF1global_H_Node1local
        Transform DOF0_H_local0, DOF1_H_local1;
        l_interpolation->getDOFtoLocalTransform(b, DOF0_H_local0, DOF1_H_local1);

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
        beamMatrices->m_A0Ref = DOF0global_H_Node0local.inversed().getAdjointMatrix();
        beamMatrices->m_A1Ref = DOF1global_H_Node1local.inversed().getAdjointMatrix();

        /////////////////////////////////////// COMPUTATION OF THE MASS AND STIFFNESS  MATRIX (LOCAL)
        /// compute the local mass matrices
        if(d_computeMass.getValue())
        {
            computeMass(b, (*beamMatrices));

        }

        /// IF RIGIDIFICATION: no stiffness forces:
        if(node0Idx==node1Idx)
            continue;

        /// compute the local stiffness matrices
        computeStiffness(b, (*beamMatrices));

        /////////////////////////////COMPUTATION OF THE STIFFNESS FORCE
        /// compute the current local displacement of the beam (6dofs)
        /// 1. get the rest transformation from local to 0 and local to 1
        Transform local_H_local0_rest,local_H_local1_rest;
        l_interpolation->getSplineRestTransform(b,local_H_local0_rest, local_H_local1_rest);

        ///2. computes the local displacement of 0 and 1 in frame local:
        SpatialVector u0 = local_H_local0.CreateSpatialVector() - local_H_local0_rest.CreateSpatialVector();
        SpatialVector u1 = local_H_local1.CreateSpatialVector() - local_H_local1_rest.CreateSpatialVector();

        /// 3. put the result in a Vec6
        Vec6 U0local, U1local;

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
            Real rest_length = l_interpolation->getLength(b);
            l_interpolation->getSplinePoints(b,x,P0,P1,P2,P3);
            l_interpolation->computeActualLength(length, P0,P1,P2,P3);

            U0local[0]=(-length+rest_length)/2;
            U1local[0]=( length-rest_length)/2;
        }

        /// compute the force in the local frame:
        Vec6 f0 = beamMatrices->m_K00 * U0local + beamMatrices->m_K01 * U1local;
        Vec6 f1 = beamMatrices->m_K10 * U0local + beamMatrices->m_K11 * U1local;

        /// ADD the effect of the pressure in the axial direction
        // inner radius of the tube
        Real r = l_interpolation->d_innerRadius.getValue();

        if (r==(Real)0){
            msg_warning()<<" Inflatable Beam Force Field suppose that the interpolation is a tube ";
        }

        /// effect of pressure is added as an external force (- internal force)
        f1[0]-= r*r*M_PI*this->d_pressure.getValue();
        f0[0]+= r*r*M_PI*this->d_pressure.getValue();



        /// compute the force in the global frame
        Vec6 F0_ref = beamMatrices->m_A0Ref.multTranspose(f0);
        Vec6 F1_ref = beamMatrices->m_A1Ref.multTranspose(f1);

        /// Add this force to vector f (negative as it is internal force)
        for (unsigned int i=0; i<6; i++)
        {
            f[node0Idx][i]-=F0_ref[i];
            f[node1Idx][i]-=F1_ref[i];
        }
    }

    if(d_computeMass.getValue())
    {
        /// add gravity:
        addMDx(mparams, dataf, d_dataG, 1.0);
    }
}


template<class DataTypes>
void AdaptiveInflatableBeamForceField<DataTypes>::addDForce(const MechanicalParams* mparams,
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
void AdaptiveInflatableBeamForceField<DataTypes>::addKToMatrix(const MechanicalParams* mparams,
                                                            const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef matrixRef = matrix->getMatrix(mstate);
    Real k = (Real)mparams->kFactor();

    unsigned int numBeams = l_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        BeamLocalMatrices &beamLocalMatrix = m_localBeamMatrices[b];
        l_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        if(node0Idx==node1Idx)
            continue;

        // matrices in global frame
        Matrix6x6 K00, K01, K10, K11;
        K00=beamLocalMatrix.m_A0Ref.multTranspose( ( beamLocalMatrix.m_K00 * beamLocalMatrix.m_A0Ref  )  );
        K01=beamLocalMatrix.m_A0Ref.multTranspose( ( beamLocalMatrix.m_K01 * beamLocalMatrix.m_A1Ref  )  );
        K10=beamLocalMatrix.m_A1Ref.multTranspose( ( beamLocalMatrix.m_K10 * beamLocalMatrix.m_A0Ref  )  );
        K11=beamLocalMatrix.m_A1Ref.multTranspose( ( beamLocalMatrix.m_K11 * beamLocalMatrix.m_A1Ref  )  );

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



/////////////////////////////////////
/// Visualization
/////////////////////////////////////

template<class DataTypes>
void AdaptiveInflatableBeamForceField<DataTypes>::draw(const VisualParams *vparams)
{
    if (!vparams->displayFlags().getShowForceFields() && !vparams->displayFlags().getShowBehaviorModels()) return;
    if (!mstate) return;

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

}


template<class DataTypes>
void AdaptiveInflatableBeamForceField<DataTypes>::drawElement(const VisualParams *vparams, int beam,
                                                           Transform &global_H0_local, Transform &global_H1_local)
{
    double length = (double) l_interpolation->getLength(beam);

    /// ARROWS
    Vec3 sizeArrows (length/4., length/8., length/8.);

    vparams->drawTool()->drawFrame(global_H0_local.getOrigin(), global_H0_local.getOrientation(), sizeArrows);
    vparams->drawTool()->drawFrame(global_H1_local.getOrigin(), global_H1_local.getOrientation(), sizeArrows);
}


} // namespace sofa::component::forcefield::_AdaptiveInflatableBeamForceField_

