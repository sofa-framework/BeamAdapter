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


#include "initBeamAdapter.h"
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>

#include <sofa/core/objectmodel/BaseObject.h>


namespace sofa
{

namespace component
{

namespace fem
{

using sofa::helper::vector;
using namespace sofa::core::topology;
using namespace sofa::defaulttype;


/** Compute Finite Element elastic force and mass based on Adaptive 6D beam elements.
  - Adaptive beam interpolation
  - Adaptive Force and Mass computation
  - Adaptive Mapping

  TODO : put in a separate class what is specific to wire shape !
 */


/// AdaptiveBeam Interpolation provides the basis of the Beam computation
/// As the computation is adaptive, the interpolation can be modified at each time step.
template<class DataTypes>
class BeamInterpolation : public virtual sofa::core::objectmodel::BaseObject
{
public:

	SOFA_CLASS( SOFA_TEMPLATE(BeamInterpolation, DataTypes) , sofa::core::objectmodel::BaseObject);

	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::VecReal VecReal;
	typedef typename DataTypes::Coord Coord;
	typedef typename DataTypes::Deriv Deriv;
	typedef typename Coord::value_type Real;
	typedef unsigned int Index;
	typedef BaseMeshTopology::EdgeID ElementID;
	typedef sofa::helper::vector<BaseMeshTopology::EdgeID> VecElementID;
	typedef sofa::helper::vector<BaseMeshTopology::Edge> VecEdges;
	typedef helper::vector<unsigned int> VecIndex;


	typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
	typedef typename sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;

    typedef Vec<2, Real> Vec2;
	typedef Vec<3, Real> Vec3;
	typedef Vec<6, Real> Vec6;


	BeamInterpolation()
	: radius(initData(&radius, (Real)1.0f, "radius", "radius of the beam (for now only constant radius are used)"))
	, innerRadius(initData(&innerRadius, (Real)0.0f, "innerRadius", "inner radius of the beam if it applies"))
	, dofsAndBeamsAligned(initData(&dofsAndBeamsAligned, true, "dofsAndBeamsAligned", "if false, a transformation for each beam is computed between the DOF and the beam nodes"))
	, m_edgeList(initData(&m_edgeList, "edgeList", "list of the edge in the topology that are concerned by the Interpolation"))
	, m_lengthList(initData(&m_lengthList, "lengthList", "list of the length of each beam"))
	, m_DOF0TransformNode0(initData(&m_DOF0TransformNode0, "DOF0TransformNode0", "Optional rigid transformation between the degree of Freedom and the first node of the beam"))
	, m_DOF1TransformNode1(initData(&m_DOF1TransformNode1, "DOF1TransformNode1", "Optional rigid transformation between the degree of Freedom and the second node of the beam"))
	, m_curvAbsList(initData(&m_curvAbsList, "curvAbsList", ""))
	, _topology(NULL)
	, _mstate(NULL)
	{
		this->brokenInTwo=false;
		_isControlled=false;
		this->_numBeamsNotUnderControl = 0;
	}

	virtual ~BeamInterpolation(){}

	void init();
	void bwdInit();
	void reinit(){init(); bwdInit(); }
	void reset(){bwdInit(); this->_numBeamsNotUnderControl=0;}

	/**
	 * @brief Returns true if the interpolation is specified in the scene file (case of saved executed scenes...) 
	 */
	bool interpolationIsAlreadyInitialized();

	unsigned int getNumBeamsNotUnderControl(){return this->_numBeamsNotUnderControl;}
	unsigned int getNumBeams(){return this->m_edgeList.getValue().size();}

//	VecElementID &getEdgeList(){return this->Edge_List;}

	void getAbsCurvXFromBeam(int beam, Real& x_curv);
	void getAbsCurvXFromBeam(int beam, Real& x_curv_start, Real& x_curv_end);

	void getDOFtoLocalTransform(unsigned int edgeInList,Transform &DOF0_H_local0, Transform &DOF1_H_local1);

	void getDOFtoLocalTransformInGlobalFrame(unsigned int edgeInList, Transform &DOF0Global_H_local0, Transform &DOF1Global_H_local1, const VecCoord &x);


	void computeTransform(unsigned int edgeInList,  Transform &global_H0_local,  Transform &global_H1_local,
			Transform &local0_H_local1,  Quat& local_R_local0, const VecCoord &x);

	void computeTransform2(unsigned int edgeInList,  Transform &global_H_local0,  Transform &global_H_local1, const VecCoord &x);

	void getTangent(Vec3& t, const Real& baryCoord, const Transform &global_H_local0, const Transform &global_H_local1,const Real &L);

	void getNodeIndices(unsigned int edgeInList, unsigned int &node0Idx, unsigned int &node1Idx );

	void getInterpolationParam(unsigned int edgeInList, Real &_L, Real &_A, Real &_Iy , Real &_Iz,
			Real &_Asy, Real &_Asz, Real &J);


	// spline base interpolation of points and transformation
	void interpolatePointUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos, const VecCoord &x, Vec3& posResult);
	void getSplinePoints(unsigned int edgeInList, const VecCoord &x, Vec3& P0, Vec3& P1, Vec3& P2, Vec3 &P3);

	// getLength / setLength => provides the rest length of each spline
	Real getLength(unsigned int edgeInList)
	{
		Real _L = this->m_lengthList.getValue()[edgeInList];
		return _L;
	}

	void setLength(unsigned int edgeInList, Real &length)
	{
		vector<double> &lengthList = *m_lengthList.beginEdit();
		lengthList[edgeInList] = length;
		m_lengthList.endEdit();
	}

	// computeActualLength => given the 4 control points of the spline, it provides an estimate of the length (using gauss points integration)
	void computeActualLength(Real &length, const Vec3& P0, const Vec3& P1, const Vec3& P2, const Vec3 &P3);


	void computeStrechAndTwist(unsigned int edgeInList, const VecCoord &x, Vec3 &ResultNodeO, Vec3 &ResultNode1);
	void InterpolateTransformUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos, const VecCoord &x, Transform &global_H_localInterpol);
    // generic implementation of the interpolation =>TODO?  could:migrate to Solidtypes files ?
	void InterpolateTransformUsingSpline(Transform& global_H_localResult, const Real &baryCoord, const Transform &global_H_local0, const Transform &global_H_local1,const Real &L);

	void InterpolateTransformAndVelUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos, const VecCoord &x, const VecDeriv &v,
												Transform &global_H_localInterpol, Deriv &v_interpol);


	// 3DOF mapping
	void MapForceOnNodeUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos, const VecCoord& x, const Vec3& finput,
									SpatialVector& FNode0output, SpatialVector& FNode1output );

	// 6DoF mapping
	void MapForceOnNodeUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos, const VecCoord& x, const SpatialVector& f6DofInput,
										SpatialVector& FNode0output, SpatialVector& FNode1output );



	// compute the total bending Rotation Angle while going through the Spline (to estimate the curvature)
	void ComputeTotalBendingRotationAngle(Real& BendingAngle, const Real& dx_computation, const Transform &global_H_local0, const Transform &global_H_local1,const Real &L,
											const Real& baryCoordMin, const Real& baryCoordMax);


	void RotateFrameForAlignX(const Quat &input,  Vec3 &x, Quat &output);


	unsigned int getStateSize() const
	{
		if (this->_mstate==NULL)
		{
			serr << "WARNING no _mstate found" << sendl;
			return 0 ;
		}
		else
                {
//                    std::cout<<" get mstate named "<<this->_mstate->name<<"  size ="<<this->_mstate->getSize()<<std::endl;
			return this->_mstate->getSize();
		}
	}

	struct BeamSection{
		// double _L; //length
		double _r; //radius of the section
		double _rInner; //inner radius of the section if beam is hollow
		double _Iy;
		double _Iz; //Iz is the cross-section moment of inertia (assuming mass ratio = 1) about the z axis;
		double _J;  //Polar moment of inertia (J = Iy + Iz)
		double _A; // A is the cross-sectional area;
		double _Asy; //_Asy is the y-direction effective shear area =  10/9 (for solid circular section) or 0 for a non-Timoshenko beam
		double _Asz; //_Asz is the z-direction effective shear area;
	};
	BeamSection &getBeamSection(int /*edgeIndex*/ ){return this->_constantRadius;}

    //template< class TReal >
	struct CurvAbscissa
	{
		CurvAbscissa()
		{
		}

		CurvAbscissa(Real r1, Real r2)
			: m_abscissa(std::make_pair(r1, r2))
		{
		}

		/// write to an output stream
        inline friend std::ostream& operator << (std::ostream& out, const CurvAbscissa& a)
        {
            out << a.m_abscissa.first << " " << a.m_abscissa.second;
            return out;
        }

		/// read from an input stream
		inline friend std::istream& operator >> (std::istream& in, CurvAbscissa& a)
		{
			in >> a.m_abscissa.first >> a.m_abscissa.second;
			return in;
		}

		const Real &first() const {return m_abscissa.first;}
		const Real &second() const {return m_abscissa.second;}

        std::pair< Real, Real > m_abscissa;
	};

	Data<Real> radius;
	Data<Real> innerRadius;
	Data<bool> dofsAndBeamsAligned;

	///////// for AdaptiveControllers
	bool isControlled(){return _isControlled;}
	void setControlled(bool value){_isControlled=value;}



	virtual void clear();

	virtual void addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1, const Real &angle);

	virtual void getSamplingParameters(helper::vector<Real>& /*xP_noticeable*/, helper::vector< int>& /*nbP_density*/)
	{
		serr<<"getSamplingParameters is not implemented when _restShape== NULL : TODO !! "<<sendl;
	}

	virtual Real getRestTotalLength()
	{
		Real le(0.0);
		const vector< double > &lengthList = m_lengthList.getValue();

		for (unsigned int i = 0; i < lengthList.size(); i++)
			le += lengthList[i];

		return le;
	}

	virtual void getCollisionSampling(Real &dx, const Real& /*x_localcurv_abs*/)
	{
		unsigned int numLines = 30;
		dx = getRestTotalLength()/numLines;
	}

	virtual void getNumberOfCollisionSegment(Real &dx, unsigned int &numLines)
	{
		numLines = 30;
		dx = getRestTotalLength()/numLines;
	}


	virtual void getYoungModulusAtX(int /*beamId*/,Real& /*x_curv*/, Real& youngModulus, Real& cPoisson)
	{
		youngModulus = (Real) 1000000.0;
		cPoisson     = (Real) 0.4;
	}


    void setTransformBetweenDofAndNode(int beam, const Transform &DOF_H_Node, unsigned int zeroORone )
    {
        if (beam > (int) (m_DOF0TransformNode0.getValue().size()-1) || beam > (int) (m_DOF1TransformNode1.getValue().size()-1))
        {
            serr<<"WARNING setTransformBetweenDofAndNode on non existing beam"<<sendl;
            return;
        }

        if (!zeroORone)
        {
			vector<Transform> &DOF0_TransformNode0 = *m_DOF0TransformNode0.beginEdit();
            
			DOF0_TransformNode0[beam] = DOF_H_Node;

			m_DOF0TransformNode0.endEdit();
        }
        else
        {
			vector<Transform> &DOF1_TransformNode1 = *m_DOF1TransformNode1.beginEdit();

            DOF1_TransformNode1[beam] = DOF_H_Node;

			m_DOF1TransformNode1.endEdit();
        }
    }


	virtual void getRestTransform(unsigned int edgeInList, Transform &local0_H_local1_rest);
	virtual void getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest);
	virtual void getBeamAtCurvAbs(const Real& x_input, unsigned int &edgeInList_output, Real& baryCoord_output, unsigned int start=0);
	virtual bool breaksInTwo(const Real &x_min_out,  Real &x_break, int &numBeamsNotUnderControlled );


	/// Pre-construction check method called by ObjectFactory.
	/// Check that DataTypes matches the MechanicalState.
	template<class T>
	static bool canCreate(T* obj, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
	{
		if (dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes>*>(context->getMechanicalState()) == NULL)
		{
			return false;
		}
		return BaseObject::canCreate(obj, context, arg);
	}


	virtual std::string getTemplateName() const
	{
		return templateName(this);
	}

	static std::string templateName(const BeamInterpolation<DataTypes>* = NULL)
	{
		return DataTypes::Name();
	}


protected :
	// DATA INPUT (that could change in real-time)
	//1.m_edgeList : list of the edge in the topology that are concerned by the Interpolation
	Data< VecElementID > m_edgeList;
	const VecEdges *_topologyEdges;

	//2.m_lengthList: list of the length of each beam
	Data< vector< double > > m_lengthList;

	//3. (optional) apply a rigid Transform between the degree of Freedom and the first node of the beam
	// Indexation based on the num of Edge
	Data< vector< Transform > > m_DOF0TransformNode0;

	//4. (optional) apply a rigid Transform between the degree of Freedom and the second node of the beam
	Data< vector< Transform > > m_DOF1TransformNode1;

    Data< vector< CurvAbscissa > > m_curvAbsList;

	// GEOMETRICAL COMPUTATION (for now we suppose that the radius of the beam do not vary in space / in time)
	BeamSection _constantRadius;

	// Topology

	// pointer to the topology
	sofa::core::topology::BaseMeshTopology* _topology;

	// verify that the m_edgeList always contains existing edges
	bool verifyTopology();

	// pointer on mechanical state
	sofa::core::behavior::MechanicalState<DataTypes> *_mstate;

	// pointer on an external rest-shape
	//sofa::component::engine::WireRestShape<DataTypes> *_restShape;

	// bool => tells if the Beams are controlled :
	// an external controller is providing the list of edge and the length of the beams (only for a wire)


	// this->brokenInTwo = if true, the wire is in two separate parts
	bool _isControlled;
	bool brokenInTwo;
	unsigned int  _numBeamsNotUnderControl;
};

#if defined(WIN32) && !defined(SOFA_BUILD_BEAMADAPTER)
#pragma warning(disable : 4231)
#ifndef SOFA_FLOAT
extern template class SOFA_BEAMADAPTER_API BeamInterpolation<Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_BEAMADAPTER_API BeamInterpolation<Rigid3fTypes>;
#endif
#endif

} // namespace fem

} // namespace component

} // namespace sofa

#endif  /*SOFA_COMPONENT_FEM_BEAMINTERPOLATION_H*/
