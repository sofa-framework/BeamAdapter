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
#pragma once

#include <BeamAdapter/config.h>
#include <BeamAdapter/utils/BeamSection.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/component/statecontainer/MechanicalObject.h>

namespace sofa::component::fem
{

namespace _beaminterpolation_
{

using sofa::helper::OptionsGroup;
using sofa::core::topology::BaseMeshTopology;
using sofa::core::ConstVecCoordId;
using sofa::core::behavior::MechanicalState;
using sofa::component::statecontainer::MechanicalObject;

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
class BeamInterpolation : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS( SOFA_TEMPLATE(BeamInterpolation, DataTypes) , sofa::core::objectmodel::BaseObject);

    using Coord = typename DataTypes::Coord;
    using VecCoord = typename DataTypes::VecCoord;
    using Real = typename Coord::value_type;

    using Deriv = typename DataTypes::Deriv;
    using VecDeriv = typename DataTypes::VecDeriv;

    using Vec2 = sofa::type::Vec<2, Real>;
    using Vec3 = sofa::type::Vec<3, Real>;
    using Vec3NoInit = sofa::type::VecNoInit<3, Real>;
    using Quat = sofa::type::Quat<Real>;
    using VectorVec3 = type::vector <Vec3>;

    using Transform = typename sofa::defaulttype::SolidTypes<Real>::Transform;
    using SpatialVector = typename sofa::defaulttype::SolidTypes<Real>::SpatialVector;

    using PointID = BaseMeshTopology::PointID;
    using ElementID = BaseMeshTopology::EdgeID;
    using VecElementID = type::vector<BaseMeshTopology::EdgeID>;
    using VecEdges = type::vector<BaseMeshTopology::Edge>;
    
    using BeamSection = sofa::beamadapter::BeamSection;

public:
    BeamInterpolation() ;
    virtual ~BeamInterpolation() override = default;

    //////////////////////////////////// Exposing this object in the factory ///////////////////////
    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T* obj, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
    {
        if (dynamic_cast<MechanicalState<DataTypes>*>(context->getMechanicalState()) == nullptr)
        {
            arg->logError(std::string("No mechanical state with the datatype '") + DataTypes::Name() +
                "' found in the context node.");
            return false;
        }
        return BaseObject::canCreate(obj, context, arg);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////// Inherited from Base ///////////////////////////////////////
    void init() override ;
    void bwdInit() override ;
    void reinit() override ;
    void reset() override ;

    //TODO(dmarchal@cduriez) Ca me semble détourner l'API pour faire des choses par surprise. A mon avis la bonne solution
    //est d'implémenter un vrai binding Python pour BeamInterpolation. Avec une fonction updateInterpolation
    /// In the context of beam interpolation, this function (easily access with Python) is used to update the interpolation (input / output)
    void storeResetState() override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////


    void updateInterpolation();
    /**
     * @brief Returns true if the interpolation is specified in the scene file (case of saved executed scenes...)
     */
    bool interpolationIsAlreadyInitialized();
    bool verifyTopology();
    void computeCrossSectionInertiaMatrix();

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

    
    int computeTransform(const ElementID edgeInList, Transform &global_H_local0,  Transform &global_H_local1, const VecCoord &x);
    int computeTransform(const ElementID edgeInList, const PointID node0Idx, const PointID node1Idx, Transform& global_H_local0, Transform& global_H_local1, const VecCoord& x);

    

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

    /// Method to rotate a Frame define by a Quat @param input around an axis @param x, x will be normalized in method. Output is return inside @param output
    void RotateFrameForAlignX(const Quat &input, Vec3 &x, Quat &output);
    
    /// Method to rotate a Frame define by a Quat @param input around an axis @param x , x has to be normalized. Output is return inside @param output
    void RotateFrameForAlignNormalizedX(const Quat& input, const Vec3& x, Quat& output);

    unsigned int getStateSize() const ;

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

    Real          m_defaultYoungModulus;
    Real          m_defaultPoissonRatio;
    Data<type::vector<Real>>          d_defaultYoungModulus;
    Data<type::vector<Real>>          d_poissonRatio;

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

    /////////////////////////// Deprecated Methods  ////////////////////////////////////////// 
    [[deprecated("Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06")]]
    unsigned int getNumBeamsNotUnderControl() {
        msg_warning() << "Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06";
        return 0;
    }

protected :
    /// DATA INPUT (that could change in real-time)
    using StateDataTypes = sofa::defaulttype::StdVectorTypes< sofa::type::Vec<DataTypes::spatial_dimensions, Real>, sofa::type::Vec<DataTypes::spatial_dimensions, Real>, Real >;
    typename MechanicalObject<StateDataTypes>::SPtr m_StateNodes;

    ///1.m_edgeList : list of the edge in the topology that are concerned by the Interpolation
    Data< VecElementID >        d_edgeList;
    const VecEdges*             m_topologyEdges {nullptr};

    ///2.m_lengthList: list of the length of each beam
    Data< type::vector< double > >    d_lengthList;

    ///3. (optional) apply a rigid Transform between the degree of Freedom and the first node of the beam
    /// Indexation based on the num of Edge
    Data< type::vector< Transform > > d_DOF0TransformNode0;

    ///4. (optional) apply a rigid Transform between the degree of Freedom and the second node of the beam
    Data< type::vector< Transform > > d_DOF1TransformNode1;

    Data< type::vector< Vec2 > >      d_curvAbsList;

    ///5. (optional) list of the beams in m_edgeList that need to be considered for collision
    Data< sofa::type::vector<int> > d_beamCollision;

    /// INPUT / OUTPUT FOR DOING EXTERNAL COMPUTATION OF Beam Interpolation (use it as a kind of data engine)
    ///Input 1. VecID => (1) "current" Pos, Vel    (2) "free" PosFree, VelFree   (3) "rest" PosRest, V=0
    Data< OptionsGroup > d_vecID;
    ///Input 2. Vector of 2-tuples (indice of the beam   ,   barycentric between 0 and 1)
    Data< type::vector< Vec2 > > d_InterpolationInputs;

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
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_BEAMINTERPOLATION_CPP)
extern template class SOFA_BEAMADAPTER_API BeamInterpolation<sofa::defaulttype::Rigid3Types>;
#endif

} /// namespace _beaminterpolation_

using _beaminterpolation_::BeamInterpolation ;

} /// namespace sofa::component::fem
