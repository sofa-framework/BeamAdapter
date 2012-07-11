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
// C++ Implementation : AdaptiveBeamController
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
//
#ifndef SOFA_COMPONENT_ENGINE_BASERESTSHAPE_H
#define SOFA_COMPONENT_ENGINE_BASERESTSHAPE_H

#include "initBeamAdapter.h"
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/DataEngine.h>
#include <sofa/component/topology/EdgeSetTopologyModifier.h>
#include <sofa/component/topology/EdgeSetGeometryAlgorithms.h>
#include <sofa/component/topology/Edge2QuadTopologicalMapping.h>


namespace sofa
{

namespace component
{

namespace engine
{

/*! \class BaseRestShape
 * \brief Describe the shape functions in a single segment
 */
template <class DataTypes>
class BaseRestShape : public sofa::core::DataEngine
{
public:
	SOFA_CLASS(BaseRestShape,sofa::core::DataEngine);
	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::Coord    Coord   ;
	typedef typename DataTypes::Deriv    Deriv   ;
	typedef typename Coord::value_type   Real    ;
	typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
	typedef sofa::defaulttype::Vec<3, Real> Vec3;
	typedef sofa::defaulttype::Vec<2, Real> Vec2;
	typedef sofa::defaulttype::Quat Quat;
	typedef typename sofa::helper::vector<Vec2>::iterator vecIt;

	/**
	 * @brief Default Constructor.
	 */
	BaseRestShape()
	: length(initData(&length, (Real)1.0, "length", "total length of the wire instrument"))
	, density(initData(&density, "densityOfBeams", "density of beams between key points"))
	, keyPoints(initData(&keyPoints,"keyPoints","key points of the shape (curv absc)"))
	, _poissonRatio(initData(&_poissonRatio,(Real)0.49,"poissonRatio","Poisson Ratio")) 	// TODO => TABLE zith size = density.size() or keyPoints.size()-1
	, _youngModulus1(initData(&_youngModulus1,(Real)5000,"youngModulus","Young Modulus"))  // TODO => TABLE zith size = density.size() or keyPoints.size()-1
	{
	}

	virtual void init() ;
	virtual void update(){}

	/**
	 * @brief Default Destructor.
	 */
	~BaseRestShape(){}


	/*!
	 * for coils: a part of the coil instrument can be brokenIn2  (by default the point of release is the end of the straight length)
	 */
	virtual Real getReleaseCurvAbs(){return length.getValue();}


	// todo => topological change !
	virtual void releaseWirePart();


	/*!
	 * this function is called by the force field to evaluate the rest position of each beam
	 */
	virtual void getRestTransformOnX(Transform &global_H_local, const Real &x) ;

        virtual void getCollisionSampling(Real &dx, const Real & /*x_curv*/){dx=length.getValue();}

	/*!
	 * This function provides a vector with the curvilinear abscissa of the noticeable point(s)
	 * and the minimum density (number of points) between them.
	 */
	virtual void getSamplingParameters(helper::vector<Real>& xP_noticeable, helper::vector<int>& nbP_density);

	virtual void getNumberOfCollisionSegment(Real &dx, unsigned int &/*numLines*/){dx=length.getValue();}

	/*!
	 * this function gives the Young modulus and Poisson's coefficient of the beam depending on the beam position
	 */
	virtual void getYoungModulusAtX(Real& /*x_curv*/, Real& youngModulus, Real& cPoisson);


	Real getLength(){
		return length.getValue();
	}


protected:

	Data<Real> length;  							/*!< Total length of the wire instrument */
	Data< sofa::helper::vector<int> > density; 		/*!< Density of beams between key points */
	Data< sofa::helper::vector<Real> > keyPoints; 	/*!< Key points of the shape (curv absc) */
	Data<Real> _poissonRatio; 						/*!< Poisson Ratio */
	Data<Real> _youngModulus1; 						/*!< Young Modulus */


public :

	/// Construction method called by ObjectFactory.
	template<class T>
	static typename T::SPtr create(T* obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
	{
                return core::objectmodel::BaseObject::create(obj, context, arg);
	}

	virtual std::string getTemplateName() const
	{
		return templateName(this);
	}

	static std::string templateName(const BaseRestShape<DataTypes>* = NULL)
	{
		return DataTypes::Name();
	}



};


} // namespace engine

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_ENGINE_BASERESTSHAPE_H */
