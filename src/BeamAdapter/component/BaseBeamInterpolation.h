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
#pragma once

#include <BeamAdapter/config.h>

#include <BeamAdapter/component/engine/WireRestShape.h>

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/type/vector.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/component/statecontainer/MechanicalObject.h>


namespace sofa::component::fem
{

namespace _basebeaminterpolation_
{
using sofa::core::topology::BaseMeshTopology ;
using sofa::type::Quat ;
using sofa::type::Vec ;
using sofa::type::Vec3d ;
using sofa::type::vector;
using sofa::core::behavior::MechanicalState;
using sofa::component::statecontainer::MechanicalObject;

/*!
 * \class BaseBeamInterpolation
 *
 */
template<class DataTypes>
class BaseBeamInterpolation : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BaseBeamInterpolation, DataTypes) ,
        sofa::core::objectmodel::BaseObject);

    using Coord = typename DataTypes::Coord;
    using VecCoord = typename DataTypes::VecCoord;
    using Real = typename Coord::value_type;

    using Deriv = typename DataTypes::Deriv;
    using VecDeriv = typename DataTypes::VecDeriv;

    using Transform = typename sofa::defaulttype::SolidTypes<Real>::Transform;
    using SpatialVector = typename sofa::defaulttype::SolidTypes<Real>::SpatialVector;

    using Vec2 = sofa::type::Vec<2, Real>;
    using Vec3 = sofa::type::Vec<3, Real>;
    using Vec3NoInit = sofa::type::VecNoInit<3, Real>;
    using Quat = sofa::type::Quat<Real>;
    using VectorVec3 = type::vector <Vec3>;

    using PointID = BaseMeshTopology::PointID;
    using EdgeID = BaseMeshTopology::EdgeID;
    using VecEdgeID = type::vector<BaseMeshTopology::EdgeID>;
    using VecEdges = type::vector<BaseMeshTopology::Edge>;

    using BeamSection = sofa::beamadapter::BeamSection;

    BaseBeamInterpolation(/*sofa::component::engine::WireRestShape<DataTypes> *_restShape = nullptr*/);

    virtual ~BaseBeamInterpolation() = default;

    void bwdInit() override;

    static void getControlPointsFromFrame(
        const Transform& global_H_local0, const Transform& global_H_local1,
        const Real& L,
        Vec3& P0, Vec3& P1,
        Vec3& P2, Vec3& P3);

    /// Method to rotate a Frame define by a Quat @param input around an axis @param x, x will be normalized in method. Output is return inside @param output
    static void RotateFrameForAlignX(const Quat& input, Vec3& x, Quat& output);

    /// Method to rotate a Frame define by a Quat @param input around an axis @param x , x has to be normalized. Output is return inside @param output
    static void RotateFrameForAlignNormalizedX(const Quat& input, const Vec3& x, Quat& output);

    virtual void clear();

public:
    virtual void addBeam(const BaseMeshTopology::EdgeID& eID, const Real& length, const Real& x0, const Real& x1, const Real& angle);
    unsigned int getNumBeams() { return this->d_edgeList.getValue().size(); }
    
    void getAbsCurvXFromBeam(int beam, Real& x_curv);
    void getAbsCurvXFromBeam(int beam, Real& x_curv_start, Real& x_curv_end);

    /// getLength / setLength => provides the rest length of each spline using @sa d_lengthList
    virtual Real getRestTotalLength() = 0;
    Real getLength(unsigned int edgeInList);
    void setLength(unsigned int edgeInList, Real& length);
    
    /// Collision information using @sa d_beamCollision
    virtual void getCollisionSampling(Real& dx, const Real& x_localcurv_abs) = 0;
    void addCollisionOnBeam(unsigned int b);
    void clearCollisionOnBeam();

    virtual void getSamplingParameters(type::vector<Real>& xP_noticeable,
        type::vector<int>& nbP_density) = 0;
    virtual void getNumberOfCollisionSegment(Real& dx, unsigned int& numLines) = 0;


    virtual void getCurvAbsAtBeam(const unsigned int& edgeInList_input, const Real& baryCoord_input, Real& x_output) = 0;
    virtual void getSplineRestTransform(unsigned int edgeInList, Transform& local_H_local0_rest, Transform& local_H_local1_rest) = 0;
    
    /// Returns the BeamSection @sa m_beamSection corresponding to the given beam
    virtual const BeamSection& getBeamSection(sofa::Index beamId) = 0;
    /// Returns the BeamSection data depending on the beam position at the given beam, similar to @getBeamSection
    virtual void getInterpolationParameters(sofa::Index beamId, Real& _L, Real& _A, Real& _Iy, Real& _Iz, Real& _Asy, Real& _Asz, Real& J) = 0;
    /// Returns the Young modulus, Poisson's ratio and massDensity coefficient of the section at the given curvilinear abscissa
    virtual void getMechanicalParameters(sofa::Index beamId, Real& youngModulus, Real& cPoisson, Real& massDensity) = 0;


    virtual void getBeamAtCurvAbs(const Real& x_input, unsigned int& edgeInList_output, Real& baryCoord_output, unsigned int start = 0);

    int computeTransform(const EdgeID edgeInList, Transform& global_H_local0, Transform& global_H_local1, const VecCoord& x);
    int computeTransform(const EdgeID edgeInList, const PointID node0Idx, const PointID node1Idx, Transform& global_H_local0, Transform& global_H_local1, const VecCoord& x);

    void getDOFtoLocalTransform(unsigned int edgeInList, Transform& DOF0_H_local0, Transform& DOF1_H_local1);
    void getDOFtoLocalTransformInGlobalFrame(unsigned int edgeInList, Transform& DOF0Global_H_local0, Transform& DOF1Global_H_local1, const VecCoord& x);
    void setTransformBetweenDofAndNode(int beam, const Transform& DOF_H_Node, unsigned int zeroORone);

    void getTangent(Vec3& t, const Real& baryCoord,
        const Transform& global_H_local0, const Transform& global_H_local1, const Real& L);

    int getNodeIndices(unsigned int edgeInList, unsigned int& node0Idx, unsigned int& node1Idx);

    void getSplinePoints(unsigned int edgeInList, const VecCoord& x, Vec3& P0, Vec3& P1, Vec3& P2, Vec3& P3);
    unsigned int getStateSize() const;

    void computeActualLength(Real& length, const Vec3& P0, const Vec3& P1, const Vec3& P2, const Vec3& P3);

    void computeStrechAndTwist(unsigned int edgeInList, const VecCoord& x, Vec3& ResultNodeO, Vec3& ResultNode1);


    
    ///vId_Out provides the id of the multiVecId which stores the position of the Bezier Points
    void updateBezierPoints(const VecCoord& x, sofa::core::VecCoordId& vId_Out);
    void updateBezierPoints(const VecCoord& x, unsigned int index, VectorVec3& v);


    /// spline base interpolation of points and transformation
    void interpolatePointUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos, const VecCoord& x, Vec3& posResult) {
        interpolatePointUsingSpline(edgeInList, baryCoord, localPos, x, posResult, true, sofa::core::vec_id::read_access::freePosition);
    }

    void interpolatePointUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos,
        const VecCoord& x, Vec3& posResult, bool recompute, const sofa::core::ConstVecCoordId& vecXId);


    void InterpolateTransformUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos,
        const VecCoord& x, Transform& global_H_localInterpol);

    void InterpolateTransformUsingSpline(Transform& global_H_localResult, const Real& baryCoord,
        const Transform& global_H_local0, const Transform& global_H_local1, const Real& L);

    void InterpolateTransformAndVelUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos,
        const VecCoord& x, const VecDeriv& v,
        Transform& global_H_localInterpol, Deriv& v_interpol);


    /// compute the total bending Rotation Angle while going through the Spline (to estimate the curvature)
    Real ComputeTotalBendingRotationAngle(const Real& dx_computation, const Transform& global_H_local0,
        const Transform& global_H_local1, const Real& L,
        const Real& baryCoordMin, const Real& baryCoordMax);


    /// 3DOF mapping
    void MapForceOnNodeUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos,
        const VecCoord& x, const Vec3& finput,
        SpatialVector& FNode0output, SpatialVector& FNode1output);

    /// 6DoF mapping
    void MapForceOnNodeUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos,
        const VecCoord& x, const SpatialVector& f6DofInput,
        SpatialVector& FNode0output, SpatialVector& FNode1output);

protected:
    


public:
    /// DATA INPUT (that could change in real-time)
    using StateDataTypes = sofa::defaulttype::StdVectorTypes< sofa::type::Vec<DataTypes::spatial_dimensions, Real>, sofa::type::Vec<DataTypes::spatial_dimensions, Real>, Real >;
    typename MechanicalObject<StateDataTypes>::SPtr m_StateNodes;

    Data< VecEdgeID > d_edgeList;
    const VecEdges* m_topologyEdges{ nullptr };

    ///2. Vector of length of each beam. Same size as @sa d_edgeList
    Data< type::vector< double > >    d_lengthList;

    ///3. (optional) apply a rigid Transform between the degree of Freedom and the first node of the beam. Indexation based on the num of Edge
    Data< type::vector< Transform > > d_DOF0TransformNode0;

    ///4. (optional) apply a rigid Transform between the degree of Freedom and the second node of the beam. Indexation based on the num of Edge
    Data< type::vector< Transform > > d_DOF1TransformNode1;

    /// Vector of 2 Real. absciss curv of first and second point defining the beam.
    Data< type::vector< Vec2 > >      d_curvAbsList;

    ///5. (optional) list of the beams in m_edgeList that need to be considered for collision
    Data< sofa::type::vector<int> > d_beamCollision;

    Data<bool>          d_dofsAndBeamsAligned;


    /// pointer to the topology
    BaseMeshTopology* m_topology{ nullptr };

    /// pointer on mechanical state
    MechanicalState<DataTypes>* m_mstate{ nullptr };
};


#if !defined(SOFA_PLUGIN_BEAMADAPTER_BaseBeamInterpolation_CPP)
extern template class SOFA_BEAMADAPTER_API BaseBeamInterpolation<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace _basebeaminterpolation_

/// Import the privately defined into the expected sofa namespace.
using _basebeaminterpolation_::BaseBeamInterpolation ;

} // namespace sofa::component::fem
