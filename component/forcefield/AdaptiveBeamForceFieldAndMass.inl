/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
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
#ifndef SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEAMFORCEFIELDANDMASS_INL
#define SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEAMFORCEFIELDANDMASS_INL

#include <sofa/core/behavior/ForceField.inl>
#include "AdaptiveBeamForceFieldAndMass.h"

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/helper/decompose.h>


#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/helper/OptionsGroup.h>

#include <sofa/helper/gl/Cylinder.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/helper/gl/Axis.h>
#include <sofa/core/visual/VisualParams.h>




namespace sofa
{

namespace component
{

namespace forcefield
{

namespace _adaptivebeamforcefieldandmass_
{

/* ************* ADAPTIVE FORCEFIELD_AND_MASS ************** */
using sofa::core::behavior::ForceField ;
using sofa::core::objectmodel::BaseContext ;
using sofa::defaulttype::Vector3 ;
using sofa::defaulttype::Quat ;
using sofa::helper::ReadAccessor ;
using sofa::core::ConstVecCoordId ;
using std::set ;

template <class DataTypes>
AdaptiveBeamForceFieldAndMass<DataTypes>::AdaptiveBeamForceFieldAndMass()
    : d_timoshenko(initData(&d_timoshenko,true,"timoshenko","use Timoshenko beam (non-null section shear area)"))
    , d_computeMass(initData(&d_computeMass,true,"computeMass","if false, only compute the stiff elastic model"))
    , d_massDensity(initData(&d_massDensity,(Real)1.0,"massDensity", "Density of the mass (usually in kg/m^3)" ))
    , d_shearStressComputation(initData(&d_shearStressComputation, true, "shearStressComputation","if false, suppress the shear stress in the computation"))
    , d_reinforceLength(initData(&d_reinforceLength, false, "reinforceLength", "if true, a separate computation for the error in elongation is peformed"))
    , m_interpolation(initLink("interpolation","Path to the Interpolation component on scene"))
    , m_instrumentParameters(initLink("instrumentParameters", "link to an object specifying physical parameters based on abscissa"))
{
    d_localBeamMatrices.resize(2);
}

template <class DataTypes>
AdaptiveBeamForceFieldAndMass<DataTypes>::~AdaptiveBeamForceFieldAndMass()
{
}

template <class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::init()
{
    if(!m_interpolation)
        m_interpolation.set(dynamic_cast<BaseContext *>(this->getContext())->get<BInterpolation>(BaseContext::Local));

    if(!m_interpolation)
        serr<<"No Beam Interpolation found, the component can not work!"<<sendl;

    this->ForceField<DataTypes>::init();

    //TODO(dmarchal 2017-05-17) Please specify who/when this will be done
    // TODO : Why do we resize this vector ?
    d_localBeamMatrices.resize(2);
}

template <class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::reinit()
{
    init();
}


template <class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::computeGravityVector()
{


    Vec3 g = this->getContext()->getGravity();

    VecDeriv& _G = *d_dataG.beginEdit();
    _G.resize(m_interpolation->getStateSize());


    this->m_gravity = Vec3(g[0],g[1],g[2]);


    if(_G.size()==0){
        serr<<"WARNING : _G.size = 0"<<sendl;
        return; }

    for (unsigned int i=0; i<_G.size(); i++)
    {
        _G[i][0]=g[0];      _G[i][1]=g[1];      _G[i][2]=g[2];
        _G[i][3]=(Real)0.0; _G[i][4]=(Real)0.0; _G[i][5]=(Real)0.0;
    }

    d_dataG.endEdit();
}

///////////////////////////// MASS INTERFACE

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addMDx(const MechanicalParams* /*mparams*/ , DataVecDeriv& dataf, const DataVecDeriv& datadx, double factor)
{
    VecDeriv& f = *dataf.beginEdit() ;
    const VecDeriv& dx = datadx.getValue();

    unsigned int numBeams = m_interpolation->getNumBeams();

    //TODO(dmarchal 2017-05-17) So what is the answer to the question in a comment ?
    if (f.size()!=dx.size())
        f.resize(dx.size()); // will reset the value ????

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        m_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        applyMassLarge( f, dx, b, node0Idx, node1Idx, factor );

    }

    dataf.endEdit() ;
}



///////////////////////////// FORCE - FIELD INTERFACE
template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addForce (const MechanicalParams* mparams ,
                                                         DataVecDeriv& dataf,
                                                         const DataVecCoord& datax, const DataVecDeriv& v)
{
    VecDeriv& f = *dataf.beginEdit() ;
    const VecCoord& x = datax.getValue();

    //TODO(dmarchal 2017-05-17) So what is the answer to the question in a comment ?
    f.resize(x.size()); // will reset the value ????

    unsigned int numBeams = m_interpolation->getNumBeams();
    d_localBeamMatrices.resize(numBeams);

    if(d_computeMass.getValue())
    {
        computeGravityVector();
    }

    /// TODO:
    ///* Redimentionner _localBeamMatrices
    ///* Calculer les rotation et les transformations
    ///* Calculer la matrice "locale"
    ///* Calculer la force exercée par chaque beam
    ///* Calculer la force exercée par la gravitée
    for (unsigned int b=0; b<numBeams; b++)
    {
        ///find the indices of the nodes
        unsigned int node0Idx, node1Idx;
        m_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        ///find the beamMatrices:
        BeamLocalMatrices  *beamMatrices = &d_localBeamMatrices[b]  ;//new BeamLocalMatrices();

        ///////////// new : Calcul du repère local de la beam & des transformations adequates///////////////
        Transform global_H_local0, global_H_local1;

        /// 1. get the current transform of the beam:
        sout << "in addForce" << sendl;
        m_interpolation->computeTransform2(b, global_H_local0, global_H_local1, x);

        /// 2. Computes the frame of the beam based on the spline interpolation:
        Transform global_H_local;
        Real baryX = 0.5;
        Real L = m_interpolation->getLength(b);

        m_interpolation->InterpolateTransformUsingSpline(global_H_local, baryX, global_H_local0, global_H_local1, L);

        Transform local_H_local0 = global_H_local.inversed()*global_H_local0;
        Transform local_H_local1 = global_H_local.inversed()*global_H_local1;

        /// 3. Computes the transformation from the DOF (in global frame) to the node's local frame DOF0global_H_Node0local and DOF1global_H_Node1local
        Transform DOF0_H_local0, DOF1_H_local1;
        m_interpolation->getDOFtoLocalTransform(b, DOF0_H_local0, DOF1_H_local1);

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
        beamMatrices->loc0_Ad_ref = DOF0global_H_Node0local.inversed().getAdjointMatrix();
        beamMatrices->loc1_Ad_ref = DOF1global_H_Node1local.inversed().getAdjointMatrix();



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
        m_interpolation->getSplineRestTransform(b,local_H_local0_rest, local_H_local1_rest);

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
            Real rest_length = m_interpolation->getLength(b);
            m_interpolation->getSplinePoints(b,x,P0,P1,P2,P3);
            m_interpolation->computeActualLength(length, P0,P1,P2,P3);

            U0local[0]=(-length+rest_length)/2;
            U1local[0]=( length-rest_length)/2;
        }

        if (!d_shearStressComputation.getValue())
        {
            /////////////////// TEST //////////////////////
            /// test: correction due to spline computation;
            Vec3 ResultNode0, ResultNode1;
            m_interpolation->computeStrechAndTwist(b, x, ResultNode0, ResultNode1);

            Real ux0 =-ResultNode0[0] + m_interpolation->getLength(b)/2;
            Real ux1 = ResultNode1[0] - m_interpolation->getLength(b)/2;

            U0local[0] = ux0;
            U1local[0] = ux1;

            U0local[3] =-ResultNode0[2];
            U1local[3] = ResultNode1[2];

            //////////////////////////////////////////////////
        }

        /// compute the force in the local frame:
        Vec6 f0 = beamMatrices->k_loc00 * U0local + beamMatrices->k_loc01 * U1local;
        Vec6 f1 = beamMatrices->k_loc10 * U0local + beamMatrices->k_loc11 * U1local;

        /// compute the force in the global frame
        Vec6 F0_ref = beamMatrices->loc0_Ad_ref.multTranspose(f0);
        Vec6 F1_ref = beamMatrices->loc1_Ad_ref.multTranspose(f1);

        /// Add this force to vector f
        for (unsigned int i=0; i<6; i++)
        {
            f[node0Idx][i]-=F0_ref[i];
            f[node1Idx][i]-=F1_ref[i];
        }
    }

    if(d_computeMass.getValue())
    {

        /// add gravity:
        this->addMDx(mparams , dataf,d_dataG,1.0);
    }
    dataf.endEdit() ;

}

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addDForce(const MechanicalParams* mparams,
                                                         DataVecDeriv& datadF, const DataVecDeriv& datadX )
{
    VecDeriv& df = *datadF.beginEdit();
    const VecDeriv& dx=datadX.getValue();
    const double kFactor=mparams->kFactor();

    df.resize(dx.size()); // will reset the value ????

    unsigned int numBeams = m_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        m_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        applyStiffnessLarge( df, dx, b, node0Idx, node1Idx, kFactor );
    }

    datadF.endEdit();
}

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::computeStiffness(int beam, BeamLocalMatrices& b)
{
    Real x_curv = 0.0 ;
    Real _nu = 0.0 ;
    Real _E = 0.0 ;

    ///Get the curvilinear abscissa of the extremity of the beam
    m_interpolation->getYoungModulusAtX(beam,x_curv, _E, _nu);

    /// material parameters
    Real _G;
    _G=_E/(2.0*(1.0+_nu));

    /// interpolation & geometrical parameters
    Real _A, _L, _Iy, _Iz, _Asy, _Asz, _J;
    m_interpolation->getInterpolationParam(beam, _L, _A, _Iy , _Iz, _Asy, _Asz, _J);

    /// Temp : we only overide values for which a Data has been set in the WireRestShape
    if(m_instrumentParameters.get())
    {
        Real x_curv = 0, _rho;
        m_interpolation->getAbsCurvXFromBeam(beam, x_curv);
        m_instrumentParameters->getInterpolationParam(x_curv, _rho, _A, _Iy , _Iz, _Asy, _Asz, _J);	// The length of the beams is only known to the interpolation !
    }
    Real   phiy, phiz;
    Real L2 = (Real) (_L * _L);
    Real L3 = (Real) (L2 * _L);
    Real EIy = (Real)(_E * _Iy);
    Real EIz = (Real)(_E * _Iz);

    /// Find shear-deformation parameters
    if (_Asy == 0)
        phiy = 0.0;
    else
        phiy = (Real)(24.0*(1.0+_nu)*_Iz/(_Asy*L2));

    if (_Asz == 0)
        phiz = 0.0;
    else
        phiz = (Real)(24.0*(1.0+_nu)*_Iy/(_Asz*L2));

    b.k_loc00.clear(); b.k_loc01.clear(); b.k_loc10.clear(); b.k_loc11.clear();

    /// diagonal values
    b.k_loc00[0][0] = b.k_loc11[0][0] = _E*_A/_L;
    b.k_loc00[1][1] = b.k_loc11[1][1] = (Real)(12.0*EIz/(L3*(1.0+phiy)));
    b.k_loc00[2][2] = b.k_loc11[2][2]   = (Real)(12.0*EIy/(L3*(1.0+phiz)));
    b.k_loc00[3][3] = b.k_loc11[3][3]   = _G*_J/_L;
    b.k_loc00[4][4] = b.k_loc11[4][4]   = (Real)((4.0+phiz)*EIy/(_L*(1.0+phiz)));
    b.k_loc00[5][5] = b.k_loc11[5][5]   = (Real)((4.0+phiy)*EIz/(_L*(1.0+phiy)));

    /// diagonal blocks
    b.k_loc00[4][2]   = (Real)(-6.0*EIy/(L2*(1.0+phiz)));
    b.k_loc00[5][1]   = (Real)( 6.0*EIz/(L2*(1.0+phiy)));
    b.k_loc11[5][1]  = -b.k_loc00[5][1];
    b.k_loc11[4][2]  = -b.k_loc00[4][2];

    /// lower non-diagonal blocks
    b.k_loc10[0][0]   = -b.k_loc00[0][0];
    b.k_loc10[1][1]   = -b.k_loc00[1][1];
    b.k_loc10[1][5]   = -b.k_loc00[5][1];
    b.k_loc10[2][2]   = -b.k_loc00[2][2];
    b.k_loc10[2][4]   = -b.k_loc00[4][2];
    b.k_loc10[3][3]   = -b.k_loc00[3][3];
    b.k_loc10[4][2]  = b.k_loc00[4][2];
    b.k_loc10[4][4]  = (Real)((2.0-phiz)*EIy/(_L*(1.0+phiz)));
    b.k_loc10[5][1]  = b.k_loc00[5][1];
    b.k_loc10[5][5]  = (Real)((2.0-phiy)*EIz/(_L*(1.0+phiy)));

    /// Make a symetric matrix with diagonal blocks
    for (int i=0; i<=5; i++)
    {
        for (int j=i+1; j<6; j++)
        {
            b.k_loc00[i][j] =  b.k_loc00[j][i];
            b.k_loc11[i][j] =  b.k_loc11[j][i];
        }
    }

    /// upper non-diagonal block : set k_loc10 as the transposed matrix of k_loc01
    b.k_loc01 = b.k_loc10.transposed();
}

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::computeMass(int beam,BeamLocalMatrices& b)
{
    /// material parameters
    Real _rho;
    _rho=this->d_massDensity.getValue();

    /// interpolation & geometrical parameters
    Real _A, _L, _Iy, _Iz, _Asy, _Asz, _J;
    m_interpolation->getInterpolationParam(beam, _L, _A, _Iy , _Iz, _Asy, _Asz, _J);

    /// Temp : we only overide values for which a Data has been set in the WireRestShape
    if(m_instrumentParameters.get())
    {
        Real x_curv = 0;
        m_interpolation->getAbsCurvXFromBeam(beam, x_curv);

        /// The length of the beams is only known to the interpolation !
        m_instrumentParameters->getInterpolationParam(x_curv, _rho, _A, _Iy , _Iz, _Asy, _Asz, _J);
    }

    Real L2 = (Real) (_L * _L);
    b.m_loc00.clear(); b.m_loc01.clear(); b.m_loc10.clear(); b.m_loc11.clear();

    /// diagonal values
    b.m_loc00[0][0] = b.m_loc11[0][0] = (Real)(1.0/3.0);
    b.m_loc00[1][1] = b.m_loc11[1][1] = (Real)(13.0/35.0 + 6.0*_Iz/(5.0*_A*L2));
    b.m_loc00[2][2] = b.m_loc11[2][2] = (Real)(13.0/35.0 + 6.0*_Iy/(5.0*_A*L2));
    b.m_loc00[3][3] = b.m_loc11[3][3] = (Real)(_J/(3.0*_A));
    b.m_loc00[4][4] = b.m_loc11[4][4] = (Real)(L2/105.0 + 2*_Iy/(15.0*_A));
    b.m_loc00[5][5] = b.m_loc11[5][5] = (Real)(L2/105.0 + 2*_Iz/(15.0*_A));

    /// diagonal blocks
    b.m_loc00[4][2]  = (Real)(-11.0*_L/210.0 - _Iy/(10*_A*_L)  );
    b.m_loc00[5][1]  = (Real)( 11.0*_L/210.0 + _Iz/(10*_A*_L)  );
    b.m_loc11[5][1]  = -b.m_loc00[5][1];
    b.m_loc11[4][2]  = -b.m_loc00[4][2];

    b.m_loc00 *= _rho*_A*_L;
    b.m_loc11 *= _rho*_A*_L;

    /// lower non-diagonal blocks
    b.m_loc10[0][0]  = (Real)(1.0/6.0);
    b.m_loc10[1][1]  = (Real)(9.0/70.0 - 6.0*_Iz/(5.0*_A*L2));
    b.m_loc10[2][2]  = (Real)(9.0/70.0 - 6.0*_Iy/(5.0*_A*L2));
    b.m_loc10[3][3]  = (Real)(_J/(6.0*_A));
    b.m_loc10[4][4]  = (Real)(-L2/140.0 - _Iy/(30.0*_A));
    b.m_loc10[5][5]  = (Real)(-L2/140.0 - _Iz/(30.0*_A));

    b.m_loc10[1][5]  = (Real)( 13*_L/420.0 - _Iz/(10.0*_A*_L));
    b.m_loc10[2][4]  = (Real)(-13*_L/420.0 + _Iy/(10.0*_A*_L));
    b.m_loc10[4][2]  = -b.m_loc10[2][4];
    b.m_loc10[5][1]  = -b.m_loc10[1][5];

    b.m_loc10 *= _rho*_A*_L;

    /// Make a symetric matrix with diagonal blocks
    for (int i=0; i<=5; i++)
    {
        for (int j=i+1; j<6; j++)
        {
            b.m_loc00[i][j] =  b.m_loc00[j][i];
            b.m_loc11[i][j] =  b.m_loc11[j][i];
        }
    }

    /// upper non-diagonal block : set k_loc10 as the transposed matrix of k_loc01
    b.m_loc01 = b.m_loc10.transposed();
}

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::applyStiffnessLarge( VecDeriv& df, const VecDeriv& dx,
                                                                    int bIndex, Index nd0Id, Index nd1Id,
                                                                    const double &factor )
{
    /// in case of rigidification:
    if(nd0Id==nd1Id)
        return;


    Vec6 U0, U1, u0, u1, f0, f1, F0, F1;
    BeamLocalMatrices &bLM = d_localBeamMatrices[bIndex];

    for (unsigned int i=0; i<6; i++)
    {
        U0[i] = dx[nd0Id][i];
        U1[i] = dx[nd1Id][i];
    }

    /// displacement in local frame
    u0 = bLM.loc0_Ad_ref*U0;
    u1 = bLM.loc1_Ad_ref*U1;

    /// internal force in local frame
    f0 = bLM.k_loc00*u0 +  bLM.k_loc01*u1;
    f1 = bLM.k_loc10*u0 +  bLM.k_loc11*u1;

    /// force in global frame
    F0 = bLM.loc0_Ad_ref.multTranspose(f0);
    F1 = bLM.loc1_Ad_ref.multTranspose(f1);

    /// put the result in df
    for (unsigned int i=0; i<6; i++)
    {
        df[nd0Id][i] -= F0[i]*factor;
        df[nd1Id][i] -= F1[i]*factor;
    }
}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::applyMassLarge( VecDeriv& df, const VecDeriv& dx,
                                                               int bIndex, Index nd0Id, Index nd1Id,
                                                               const double &factor)
{

    Vec6 A0, A1, a0, a1, f0, f1, F0, F1;
    BeamLocalMatrices &bLM = d_localBeamMatrices[bIndex];

    for (unsigned int i=0; i<6; i++)
    {
        A0[i] = dx[nd0Id][i];
        A1[i] = dx[nd1Id][i];
    }

    /// displacement in local frame
    a0 = bLM.loc0_Ad_ref*A0;
    a1 = bLM.loc1_Ad_ref*A1;

    /// internal force in local frame
    f0 = bLM.m_loc00*a0 +  bLM.m_loc01*a1;
    f1 = bLM.m_loc10*a0 +  bLM.m_loc11*a1;

    /// force in global frame
    F0 = bLM.loc0_Ad_ref.multTranspose(f0);
    F1 = bLM.loc1_Ad_ref.multTranspose(f1);

    /// put the result in df
    for (unsigned int i=0; i<6; i++)
    {
        df[nd0Id][i] += F0[i]*factor;
        df[nd1Id][i] += F1[i]*factor;
    }
}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addKToMatrix(const MechanicalParams* mparams,
                                                            const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    Real k = (Real)mparams->kFactor();

    unsigned int numBeams = m_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        BeamLocalMatrices &bLM = d_localBeamMatrices[b];
        m_interpolation->getNodeIndices( b,  node0Idx, node1Idx );


        if(node0Idx==node1Idx)
            continue;

        // matrices in global frame
        Matrix6x6 K00, K01, K10, K11;
        K00=bLM.loc0_Ad_ref.multTranspose( ( bLM.k_loc00 * bLM.loc0_Ad_ref  )  );
        K01=bLM.loc0_Ad_ref.multTranspose( ( bLM.k_loc01 * bLM.loc1_Ad_ref  )  );
        K10=bLM.loc1_Ad_ref.multTranspose( ( bLM.k_loc10 * bLM.loc0_Ad_ref  )  );
        K11=bLM.loc1_Ad_ref.multTranspose( ( bLM.k_loc11 * bLM.loc1_Ad_ref  )  );

        int index0[6], index1[6];
        for (int i=0;i<6;i++)
            index0[i] = r.offset+node0Idx*6+i;
        for (int i=0;i<6;i++)
            index1[i] = r.offset+node1Idx*6+i;

        for (int i=0;i<6;i++)
        {
            for (int j=0;j<6;j++)
            {
                r.matrix->add(index0[i], index0[j], - K00(i,j)*k);
                r.matrix->add(index0[i], index1[j], - K01(i,j)*k);
                r.matrix->add(index1[i], index0[j], - K10(i,j)*k);
                r.matrix->add(index1[i], index1[j], - K11(i,j)*k);
            }
        }
    }
}

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::addMToMatrix(const MechanicalParams *mparams,
                                                            const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    Real mFact = (Real)mparams->mFactor();

    unsigned int numBeams = m_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        BeamLocalMatrices &bLM = d_localBeamMatrices[b];
        m_interpolation->getNodeIndices( b,  node0Idx, node1Idx );

        /// matrices in global frame
        Matrix6x6 M00, M01, M10, M11;

        M00=bLM.loc0_Ad_ref.multTranspose( ( bLM.m_loc00 * bLM.loc0_Ad_ref  )  );
        M01=bLM.loc0_Ad_ref.multTranspose( ( bLM.m_loc01 * bLM.loc1_Ad_ref  )  );
        M10=bLM.loc1_Ad_ref.multTranspose( ( bLM.m_loc10 * bLM.loc0_Ad_ref  )  );
        M11=bLM.loc1_Ad_ref.multTranspose( ( bLM.m_loc11 * bLM.loc1_Ad_ref  )  );

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
void AdaptiveBeamForceFieldAndMass<DataTypes>::addMBKToMatrix(const MechanicalParams* mparams,
                                                              const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    Real kFact = (Real)mparams->kFactor();
    Real mFact = (Real)mparams->mFactor();

    Real totalMass = 0;
    unsigned int numBeams = m_interpolation->getNumBeams();

    for (unsigned int b=0; b<numBeams; b++)
    {
        unsigned int node0Idx, node1Idx;
        BeamLocalMatrices &bLM = d_localBeamMatrices[b];
        m_interpolation->getNodeIndices( b,  node0Idx, node1Idx );


        int index0[6], index1[6];
        for (int i=0;i<6;i++)
            index0[i] = r.offset+node0Idx*6+i;
        for (int i=0;i<6;i++)
            index1[i] = r.offset+node1Idx*6+i;


        if(node0Idx!=node1Idx) // no rigidification
        {
            // matrices in global frame
            Matrix6x6 K00, K01, K10, K11;

            K00=bLM.loc0_Ad_ref.multTranspose( ( bLM.k_loc00 * bLM.loc0_Ad_ref  )  );
            K01=bLM.loc0_Ad_ref.multTranspose( ( bLM.k_loc01 * bLM.loc1_Ad_ref  )  );
            K10=bLM.loc1_Ad_ref.multTranspose( ( bLM.k_loc10 * bLM.loc0_Ad_ref  )  );
            K11=bLM.loc1_Ad_ref.multTranspose( ( bLM.k_loc11 * bLM.loc1_Ad_ref  )  );

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

        M00=bLM.loc0_Ad_ref.multTranspose( ( bLM.m_loc00 * bLM.loc0_Ad_ref  )  );
        M01=bLM.loc0_Ad_ref.multTranspose( ( bLM.m_loc01 * bLM.loc1_Ad_ref  )  );
        M10=bLM.loc1_Ad_ref.multTranspose( ( bLM.m_loc10 * bLM.loc0_Ad_ref  )  );
        M11=bLM.loc1_Ad_ref.multTranspose( ( bLM.m_loc11 * bLM.loc1_Ad_ref  )  );

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

template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::draw(const VisualParams *vparams)
{
    if (!vparams->displayFlags().getShowForceFields() && !vparams->displayFlags().getShowBehaviorModels()) return;
    if (!this->mstate) return;

    ReadAccessor<Data<VecCoord> > x = this->mstate->read(ConstVecCoordId::position()) ;

    unsigned int numBeams = m_interpolation->getNumBeams();


    for (unsigned int b=0; b<numBeams; b++)
    {
        Transform global_H0_local,  global_H1_local;

        m_interpolation->computeTransform2(b, global_H0_local, global_H1_local, x.ref());
        unsigned int node0Idx, node1Idx;
        m_interpolation->getNodeIndices( b,  node0Idx, node1Idx );


        if (vparams->displayFlags().getShowBehaviorModels() && node0Idx!=node1Idx)
            drawElement(vparams, b, global_H0_local, global_H1_local);

        if(vparams->displayFlags().getShowForceFields())
        {

            // /  test ///
            std::vector<Vector3> points;
            Vec3 pos = global_H0_local.getOrigin();
            for (double i=0.0; i<1.00001; i+=0.02)
            {
                points.push_back(pos);
                Vec3 localPos(0.0,0.0,0.0);
                this->m_interpolation->interpolatePointUsingSpline(b, i, localPos, x.ref(), pos);
                points.push_back(pos);
            }


            if(node0Idx==node1Idx)
            {
                /// rigidification case !!
                vparams->drawTool()->drawLines(points,2, Vec<4,float>(0,0,1,1));
                continue;
            }
            else
            {
                /// other case
                vparams->drawTool()->drawLines(points,2, Vec<4,float>(1,0,0,1));
            }

            double Length = (double) m_interpolation->getLength(b);

            Vec3 localPos(0.0,0.0,0.0);
            Real baryX = 0.5;
            Transform global_H_localInterpol;
            this->m_interpolation->InterpolateTransformUsingSpline(b, baryX, localPos, x.ref(), global_H_localInterpol);


            Quat q =global_H_localInterpol.getOrientation();
            q.normalize();

            Vec3 P1, x,y,z;
            P1 = global_H_localInterpol.getOrigin();
            x= q.rotate(Vec3(Length/6.0,0,0));
            y= q.rotate(Vec3(0,Length/8.0,0));
            z= q.rotate(Vec3(0,0,Length/8.0));
            float radius_arrow = (float)Length/60.0f;

            vparams->drawTool()->drawArrow(P1,P1 + x, radius_arrow, Vec<4,float>(1,0,0,1));
            vparams->drawTool()->drawArrow(P1,P1 + y, radius_arrow, Vec<4,float>(1,0,0,1));
            vparams->drawTool()->drawArrow(P1,P1 + z, radius_arrow, Vec<4,float>(1,0,0,1));
        }

    }

}


template<class DataTypes>
void AdaptiveBeamForceFieldAndMass<DataTypes>::drawElement(const VisualParams *vparams, int beam,
                                                           Transform &global_H0_local, Transform &global_H1_local)
{


    double Length = (double) m_interpolation->getLength(beam);

    /// ARROWS
    Vec3 sizeArrows (Length/4, Length/8, Length/8);

    vparams->drawTool()->drawFrame(global_H0_local.getOrigin(), global_H0_local.getOrientation(), sizeArrows );
    vparams->drawTool()->drawFrame(global_H1_local.getOrigin(), global_H1_local.getOrientation(), sizeArrows );
}

} /// namespace _adaptivebeamforcefieldandmass_

} /// namespace forcefield

} /// namespace component

} /// namespace sofa

#endif /* SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEAMFORCEFIELDANDMASS_INL */
