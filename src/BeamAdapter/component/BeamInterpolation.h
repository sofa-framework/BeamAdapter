/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
 *                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
 *                                                                             *
 * This library is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *                                                                             *
 * This library is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                           *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this library; if not, write to the Free Software Foundation,     *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
 *******************************************************************************
 *                               SOFA :: Modules                               *
 *                                                                             *
 * Authors: The SOFA Team and external contributors (see Authors.txt)          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/
//
// C++ Implementation : BeamInterpolation / AdaptiveBeamForceFieldAndMass
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_FEM_BEAMINTERPOLATION_H
#define SOFA_COMPONENT_FEM_BEAMINTERPOLATION_H


#include <BeamAdapter/config.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/helper/logging/Messaging.h>
#include <sofa/helper/OptionsGroup.h>

namespace sofa
{

namespace component
{

namespace fem
{

namespace _beaminterpolation_
{


using sofa::type::vector;
using sofa::helper::OptionsGroup;
using sofa::core::topology::BaseMeshTopology;
using sofa::core::objectmodel::BaseObjectDescription ;
using sofa::core::objectmodel::BaseObject ;
using sofa::core::ConstVecCoordId ;
using sofa::defaulttype::SolidTypes ;
using sofa::type::Vec ;
using sofa::type::Quat ;
using sofa::defaulttype::Rigid3Types ;
using sofa::core::behavior::MechanicalState ;
using sofa::component::container::MechanicalObject ;

/*!
 * \class BeamInterpolation
 *
 * This class implements a Sofa Component that provide interpolation method to compute Finite Element elastic force and mass based on
 * Adaptive 6D beam elements.
 * - Adaptive beam interpolation
 * - Adaptive Force and Mass computation
 * - Adaptive Mapping
 *
 * AdaptiveBeam Interpolation provides the basis of the Beam computation
 * As the computation is adaptive, the interpolation can be modified at each time step.
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 *
 */
template<class DataTypes>
class BeamInterpolation : public virtual BaseObject
{
public:
    SOFA_CLASS( SOFA_TEMPLATE(BeamInterpolation, DataTypes) , BaseObject);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::VecReal VecReal;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef unsigned int Index;
    typedef BaseMeshTopology::EdgeID ElementID;
    typedef vector<BaseMeshTopology::EdgeID> VecElementID;
    typedef vector<BaseMeshTopology::Edge> VecEdges;
    typedef vector<unsigned int> VecIndex;

    typedef typename SolidTypes<Real>::Transform Transform;
    typedef typename SolidTypes<Real>::SpatialVector SpatialVector;

    typedef Vec<2, Real> Vec2;
    typedef Vec<3, Real> Vec3;
    typedef Vec<6, Real> Vec6;

    typedef vector<Vec<3, Real> > VectorVec3;

public:
    BeamInterpolation() ;
    virtual ~BeamInterpolation() override {}

    //////////////////////////////////// Exposing this object in the factory ///////////////////////
    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T* obj, sofa::core::objectmodel::BaseContext* context, BaseObjectDescription* arg)
    {
        if (dynamic_cast<MechanicalState<DataTypes>*>(context->getMechanicalState()) == nullptr)
        {
            return false;
        }
        return BaseObject::canCreate(obj, context, arg);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////// Inherited from Base ///////////////////////////////////////
    virtual void init() override ;
    virtual void bwdInit() override ;
    virtual void reinit() override ;
    virtual void reset() override ;

    //TODO(dmarchal@cduriez) Ca me semble détourner l'API pour faire des choses par surprise. A mon avis la bonne solution
    //est d'implémenter un vrai binding Python pour BeamInterpolation. Avec une fonction updateInterpolation
    /// In the context of beam interpolation, this function (easily access with Python) is used to update the interpolation (input / output)
    virtual void storeResetState() override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////


    void updateInterpolation();
    /**
     * @brief Returns true if the interpolation is specified in the scene file (case of saved executed scenes...)
     */
    bool interpolationIsAlreadyInitialized();
    bool verifyTopology();
    void computeCrossSectionInertiaMatrix();

    unsigned int getNumBeamsNotUnderControl(){return this->m_numBeamsNotUnderControl;}
    unsigned int getNumBeams(){return this->d_edgeList.getValue().size();}

    void getAbsCurvXFromBeam(int beam, Real& x_curv);
    void getAbsCurvXFromBeam(int beam, Real& x_curv_start, Real& x_curv_end);

    static void getControlPointsFromFrame(
                                const Transform& global_H_local0, const Transform& global_H_local1,
                                const Real& L,
                                Vec3& P0, Vec3& P1,
                                Vec3& P2, Vec3& P3);


    void getDOFtoLocalTransform(unsigned int edgeInList,Transform &DOF0_H_local0, Transform &DOF1_H_local1);
    void getDOFtoLocalTransformInGlobalFrame(unsigned int edgeInList, Transform &DOF0Global_H_local0, Transform &DOF1Global_H_local1, const VecCoord &x);


    int computeTransform(unsigned int edgeInList,  Transform &global_H0_local,  Transform &global_H1_local,
                         Transform &local0_H_local1,  Quat<Real>& local_R_local0, const VecCoord &x);

    int computeTransform2(unsigned int edgeInList,
                          Transform &global_H_local0,  Transform &global_H_local1, const VecCoord &x);

    void getTangent(Vec3& t, const Real& baryCoord,
                    const Transform &global_H_local0, const Transform &global_H_local1,const Real &L);

    int getNodeIndices(unsigned int edgeInList, unsigned int &node0Idx, unsigned int &node1Idx );

    void getInterpolationParam(unsigned int edgeInList, Real &_L, Real &_A, Real &_Iy , Real &_Iz,
                               Real &_Asy, Real &_Asz, Real &J);

    /// spline base interpolation of points and transformation
    void interpolatePointUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos, const VecCoord &x, Vec3& posResult) {
        interpolatePointUsingSpline(edgeInList,baryCoord,localPos,x,posResult,true, ConstVecCoordId::position());
    }

    void interpolatePointUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos,
                                     const VecCoord &x, Vec3& posResult, bool recompute, const ConstVecCoordId &vecXId);


    void getTangentUsingSplinePoints(unsigned int edgeInList, const Real& baryCoord, const ConstVecCoordId &vecXId, Vec3& t );

    void getSplinePoints(unsigned int edgeInList, const VecCoord &x, Vec3& P0, Vec3& P1, Vec3& P2, Vec3 &P3);

    ///vId_Out provides the id of the multiVecId which stores the position of the Bezier Points
    void updateBezierPoints( const VecCoord &x, sofa::core::VecCoordId &vId_Out);
    void updateBezierPoints( const VecCoord &x, unsigned int index, VectorVec3& v);

    /// getLength / setLength => provides the rest length of each spline
    Real getLength(unsigned int edgeInList) ;
    void setLength(unsigned int edgeInList, Real &length) ;

    /// computeActualLength => given the 4 control points of the spline, it provides an estimate of the length (using gauss points integration)
    void computeActualLength(Real &length, const Vec3& P0, const Vec3& P1, const Vec3& P2, const Vec3 &P3);

    void computeStrechAndTwist(unsigned int edgeInList, const VecCoord &x, Vec3 &ResultNodeO, Vec3 &ResultNode1);
    void InterpolateTransformUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos,
                                         const VecCoord &x, Transform &global_H_localInterpol);

    /// generic implementation of the interpolation =>TODO?  could:migrate to Solidtypes files ?
    void InterpolateTransformUsingSpline(Transform& global_H_localResult, const Real &baryCoord,
                                         const Transform &global_H_local0, const Transform &global_H_local1,const Real &L);

    void InterpolateTransformAndVelUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos,
                                               const VecCoord &x, const VecDeriv &v,
                                               Transform &global_H_localInterpol, Deriv &v_interpol);

    /// 3DOF mapping
    void MapForceOnNodeUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos,
                                   const VecCoord& x, const Vec3& finput,
                                   SpatialVector& FNode0output, SpatialVector& FNode1output );

    /// 6DoF mapping
    void MapForceOnNodeUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos,
                                   const VecCoord& x, const SpatialVector& f6DofInput,
                                   SpatialVector& FNode0output, SpatialVector& FNode1output );

    /// compute the total bending Rotation Angle while going through the Spline (to estimate the curvature)
    Real ComputeTotalBendingRotationAngle(const Real& dx_computation, const Transform &global_H_local0,
                                          const Transform &global_H_local1,const Real &L,
                                          const Real& baryCoordMin, const Real& baryCoordMax);

    void RotateFrameForAlignX(const Quat<Real> &input,  Vec3 &x, Quat<Real> &output);

    unsigned int getStateSize() const ;

    struct BeamSection{
        double _r; 			///<radius of the section
        double _rInner;		///<inner radius of the section if beam is hollow
        double _Iy;         /// < Iy and Iz are the cross-section moment of inertia (assuming mass ratio = 1) about the y and z axis;
        double _Iz; 		/// < see https://en.wikipedia.org/wiki/Second_moment_of_area
        double _J;  		///< Polar moment of inertia (J = Iy + Iz)
        double _A; 			///< A is the cross-sectional area;
        double _Asy; 		///< _Asy is the y-direction effective shear area =  10/9 (for solid circular section) or 0 for a non-Timoshenko beam
        double _Asz; 		///< _Asz is the z-direction effective shear area;
    };
    BeamSection &getBeamSection(int /*edgeIndex*/ ){return this->m_constantSection;}

    Data<helper::OptionsGroup>   crossSectionShape;

    /// Circular Cross Section
    Data<Real>          d_radius;
    Data<Real>          d_innerRadius;

    /// Square Cross Section
    Data<Real>          d_sideLength;

    /// Elliptic Cross Section
    Data<Real>          d_smallRadius;
    Data<Real>          d_largeRadius;

    /// Rectangular Cross Section
    Data<Real>          d_lengthY;
    Data<Real>          d_lengthZ;
    Data<bool>          d_dofsAndBeamsAligned;
    Data<Real>          d_defaultYoungModulus;
    Data<Real>          d_poissonRatio;
    Data<bool>          d_straight;

    virtual void clear();
    virtual void addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1, const Real &angle);
    virtual void getSamplingParameters(type::vector<Real>& xP_noticeable,
                                       type::vector<int>& nbP_density) ;
    virtual Real getRestTotalLength() ;
    virtual void getCollisionSampling(Real &dx, const Real& x_localcurv_abs) ;
    virtual void getNumberOfCollisionSegment(Real &dx, unsigned int &numLines) ;
    virtual void getYoungModulusAtX(int beamId,Real& x_curv, Real& youngModulus, Real& cPoisson) ;
    void setTransformBetweenDofAndNode(int beam, const Transform &DOF_H_Node, unsigned int zeroORone ) ;
    virtual void getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest);

    //TODO(dmarchal@cduriez) strange name... seems to be wire based...shouldn't it go to WireBeamInterpolation.
    virtual void getBeamAtCurvAbs(const Real& x_input, unsigned int &edgeInList_output, Real& baryCoord_output, unsigned int start=0);

    ///////// for AdaptiveControllers
    bool isControlled(){return m_isControlled;}
    void setControlled(bool value){m_isControlled=value;}

    /// Collision information
    void addCollisionOnBeam(unsigned int b) ;
    void clearCollisionOnBeam() ;

protected :
    /// DATA INPUT (that could change in real-time)
    MechanicalObject<sofa::defaulttype::Vec3Types>::SPtr m_StateNodes;

    ///1.m_edgeList : list of the edge in the topology that are concerned by the Interpolation
    Data< VecElementID >        d_edgeList;
    const VecEdges*             m_topologyEdges {nullptr};

    ///2.m_lengthList: list of the length of each beam
    Data< vector< double > >    d_lengthList;

    ///3. (optional) apply a rigid Transform between the degree of Freedom and the first node of the beam
    /// Indexation based on the num of Edge
    Data< vector< Transform > > d_DOF0TransformNode0;

    ///4. (optional) apply a rigid Transform between the degree of Freedom and the second node of the beam
    Data< vector< Transform > > d_DOF1TransformNode1;

    Data< vector< Vec2 > >      d_curvAbsList;

    ///5. (optional) list of the beams in m_edgeList that need to be considered for collision
    Data< sofa::type::vector<int> > d_beamCollision;

    /// INPUT / OUTPUT FOR DOING EXTERNAL COMPUTATION OF Beam Interpolation (use it as a kind of data engine)
    ///Input 1. VecID => (1) "current" Pos, Vel    (2) "free" PosFree, VelFree   (3) "rest" PosRest, V=0
    Data< OptionsGroup > d_vecID;
    ///Input 2. Vector of 2-tuples (indice of the beam   ,   barycentric between 0 and 1)
    Data< vector< Vec2 > > d_InterpolationInputs;

    ///Output
    Data< VecCoord > d_InterpolatedPos;
    Data< VecDeriv > d_InterpolatedVel;

    /// GEOMETRICAL COMPUTATION (for now we suppose that the radius of the beam do not vary in space / in time)
    BeamSection      m_constantSection;

    /// Topology

    /// pointer to the topology
    BaseMeshTopology* m_topology {nullptr};

    /// pointer on mechanical state
    MechanicalState<DataTypes>* m_mstate {nullptr} ;

    /// this->brokenInTwo = if true, the wire is in two separate parts
    bool  m_isControlled            {false} ;
    bool  m_brokenInTwo              {false} ;
    unsigned int m_numBeamsNotUnderControl {0} ;
};

#if !defined(SOFA_BEAMINTERPOLATION_CPP)
extern template class SOFA_BEAMADAPTER_API BeamInterpolation<Rigid3Types>;
#endif

} /// namespace _beaminterpolation_

using _beaminterpolation_::BeamInterpolation ;

} /// namespace fem

} /// namespace component

} /// namespace sofa

#endif  /*SOFA_COMPONENT_FEM_BEAMINTERPOLATION_H*/
