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

#ifndef SOFA_COMPONENT_ENGINE_WIRERESTSHAPE_H
#define SOFA_COMPONENT_ENGINE_WIRERESTSHAPE_H

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


template <class DataTypes>
class SOFA_BEAMADAPTER_API WireRestShape : public sofa::core::DataEngine
{
public:
	SOFA_CLASS(WireRestShape,sofa::core::DataEngine);
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


	 WireRestShape():
	 procedural( initData(&procedural,(bool)true,"procedural","is the guidewire shape mathemetically defined ?") )
	, fileName( initData(&fileName, std::string(""), "filename","data of filename") )
	, NonProceduralScale( initData ( &NonProceduralScale, (Real)1.0, "nonProceduralScale", "scale of the model defined by file" ) )
	, length(initData(&length, (Real)1.0, "length", "total length of the wire instrument"))
	, straightLength(initData(&straightLength, (Real)1.0, "straightLength", "length of the initial straight shape"))
	, spireDiameter(initData(&spireDiameter, (Real)0.1, "spireDiameter", "diameter of the spire"))
	, spireHeight(initData(&spireHeight, (Real)0.01, "spireHeight", "height between each spire"))
	, density(initData(&density, "densityOfBeams", "density of beams between key points"))
	, keyPoints(initData(&keyPoints,"keyPoints","key points of the shape (curv absc)"))
	, numEdges(initData(&numEdges, 10, "numEdges","number of Edges for the visual model"))
	, numEdgesCollis(initData(&numEdgesCollis,"numEdgesCollis", "number of Edges for the collision model" ))
	, _poissonRatio(initData(&_poissonRatio,(Real)0.49,"poissonRatio","Poisson Ratio"))
	, _youngModulus1(initData(&_youngModulus1,(Real)5000,"youngModulus","Young Modulus"))
	, _youngModulus2(initData(&_youngModulus2,(Real)3000,"youngModulusExtremity","youngModulus for beams at the extremity\nonly if not straight"))
	, _radius1(initData(&_radius1,(Real)1.0f,"radius","radius"))
	, _radius2(initData(&_radius2,(Real)1.0f,"radiusExtremity","radius for beams at the extremity\nonly if not straight"))
	, _innerRadius1(initData(&_innerRadius1,(Real)0.0f,"innerRadius","inner radius if it applies"))
	, _innerRadius2(initData(&_innerRadius2,(Real)0.0f,"innerRadiusExtremity","inner radius for beams at the extremity\nonly if not straight"))
	, _massDensity1(initData(&_massDensity1,(Real)1.0,"massDensity", "Density of the mass (usually in kg/m^3)" ))
	, _massDensity2(initData(&_massDensity2,(Real)1.0,"massDensityExtremity", "Density of the mass at the extremity\nonly if not straight" ))
	,edge2QuadMap(NULL)
	{
		 brokenIn2=false;
	}

	 /**
	  * @brief Default Destructor.
	  */
	 ~WireRestShape(){}

	 void RotateFrameForAlignX(const sofa::defaulttype::Quat &input, Vec3 &x, sofa::defaulttype::Quat &output);

	 void init();
	 void reinit(){ }

	 void update(){ }

	 void bwdInit();


	 // for coils: a part of the coil instrument can be brokenIn2  (by default the point of release is the end of the straight length)
	 virtual Real getReleaseCurvAbs(){return straightLength.getValue();}
	 // todo => topological change !
	 virtual void releaseWirePart();



	 // this function is called by the force field to evaluate the rest position of each beam
	 virtual void getRestTransformOnX(Transform &global_H_local, const Real &x);


	 // this function provides a vector with the curviliar abscissa of the noticeable point(s)
	 // and the minimum density (number of points) between them
	 virtual void getSamplingParameters(helper::vector<Real>& xP_noticeable, helper::vector<int>& nbP_density);

	 //this function gives the Young modulus and Poisson's coefficient of the beam depending on the beam position
	 virtual void getYoungModulusAtX(const Real& x_curv, Real& youngModulus, Real& cPoisson);

	 //this function gives the mass density and the BeamSection data depending on the beam position
	 void getInterpolationParam(const Real& x_curv, Real &_rho, Real &_A, Real &_Iy , Real &_Iz, Real &_Asy, Real &_Asz, Real &_J);

	 //Functions enabling to load and use a geometry given from OBJ external file
	 void LoadFile();

	 void InitRestConfig();

	 void getRestPosNonProcedural(Real& abs, Coord &p);

	 void computeOrientation(const Vec3& AB, const Quat& Q, Quat &result);



	 virtual Real getLength(){
		 if(brokenIn2)
			 return straightLength.getValue();
		 else
			 return length.getValue();
	 }

	 virtual void getCollisionSampling(Real &dx, const Real &x_curv);

	 virtual void getNumberOfCollisionSegment(Real &dx, unsigned int &numLines)
	 {
		 numLines = 0;
		 for (unsigned i=0; i<numEdgesCollis.getValue().size(); i++)
		 {
			 numLines += (unsigned int)numEdgesCollis.getValue()[i];
		 }
		 dx=length.getValue()/numLines;
	 }




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

	 static std::string templateName(const WireRestShape<DataTypes>* = NULL)
	 {
		 return DataTypes::Name();
	 }

protected:

	 // Analitical creation of wire shape...
	 Data<bool> procedural;
	 Data<std::string> fileName;
	 Data<Real> NonProceduralScale;
	 Data<Real> length;
	 Data<Real> straightLength;
	 Data<Real> spireDiameter;
	 Data<Real> spireHeight;
	 Data< sofa::helper::vector<int> > density;
	 Data< sofa::helper::vector<Real> > keyPoints;
	 Data<int> numEdges;
	 Data< sofa::helper::vector<int> > numEdgesCollis;


	 //User Data about the Young modulus
	 Data<Real> _poissonRatio;
	 Data<Real> _youngModulus1;
	 Data<Real> _youngModulus2;

	 //Data required for the File loading
	 sofa::helper::vector<Vec3> localRestPositions;
	 sofa::helper::vector<Transform> localRestTransforms;
	 sofa::helper::vector<Real> curvAbs;
	 double absOfGeometry;

	 // Radius
	 Data<Real> _radius1, _radius2, _innerRadius1, _innerRadius2;
	 struct BeamSection{
		double _r; //radius of the section
		double _rInner; //inner radius of the section if beam is hollow
		double _Iy;
		double _Iz; //Iz is the cross-section moment of inertia (assuming mass ratio = 1) about the z axis;
		double _J;  //Polar moment of inertia (J = Iy + Iz)
		double _A; // A is the cross-sectional area;
		double _Asy; //_Asy is the y-direction effective shear area =  10/9 (for solid circular section) or 0 for a non-Timoshenko beam
		double _Asz; //_Asz is the z-direction effective shear area;
	 };
	 BeamSection beamSection1, beamSection2;
	 Data<Real> _massDensity1, _massDensity2;

	 // broken in 2 case
	 bool brokenIn2;
	 sofa::core::topology::TopologyContainer* _topology;
	 sofa::component::topology::EdgeSetGeometryAlgorithms<DataTypes>* edgeGeo;
	 sofa::component::topology::EdgeSetTopologyModifier* edgeMod;
	 sofa::component::topology::Edge2QuadTopologicalMapping* edge2QuadMap;
	 bool edgeSetInNode;

};


} // namespace engine

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_ENGINE_WIRERESTSHAPE_H */
