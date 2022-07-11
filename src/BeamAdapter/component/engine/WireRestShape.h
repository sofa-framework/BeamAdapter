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
// C++ Implementation : WireRestShape
//
// Description:
//
// Contributors:
//   - Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
#pragma once

#include <BeamAdapter/config.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyModifier.h>
#include <sofa/core/DataEngine.h>
#include <sofa/component/topology/mapping/Edge2QuadTopologicalMapping.h>
#include <sofa/component/topology/container/dynamic/EdgeSetGeometryAlgorithms.h>
#include <sofa/core/loader/MeshLoader.h>

namespace sofa::component::engine
{

namespace _wirerestshape_
{

using sofa::type::Quat;
using sofa::type::vector;
using sofa::core::topology::TopologyContainer;
using sofa::component::topology::container::dynamic::EdgeSetGeometryAlgorithms;
using sofa::component::topology::container::dynamic::EdgeSetTopologyModifier;
using sofa::component::topology::mapping::Edge2QuadTopologicalMapping;
using sofa::core::loader::MeshLoader;

/**
 * \class WireRestShape
 * \brief Describe the shape functions on multiple segments
 *
 *  Describe the shape functions on multiple segments using curvilinear abscissa
 */
template <class DataTypes>
class WireRestShape : public sofa::core::DataEngine
{
public:
    SOFA_CLASS(WireRestShape,sofa::core::DataEngine);
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;
    typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
    typedef sofa::type::Vec<3, Real> Vec3;
    typedef sofa::type::Vec<2, Real> Vec2;
    typedef typename vector<Vec2>::iterator vecIt;

    typedef typename core::visual::VisualParams VisualParams;
    typedef typename core::objectmodel::BaseContext BaseContext;
    typedef typename core::objectmodel::BaseObjectDescription BaseObjectDescription;

    /**
     * @brief Default Constructor.
     */
     WireRestShape() ;

     /*!
      * @brief Default Destructor.
      */
     ~WireRestShape(){}

     /////////////////////////// Inherited from Base //////////////////////////////////////////
     virtual std::string getTemplateName() const override
     {
         return templateName(this);
     }

     static std::string templateName(const WireRestShape<DataTypes>* = NULL)
     {
         return DataTypes::Name();
     }

     /// Construction method called by ObjectFactory.
     template<class T>
     static typename T::SPtr create(T* obj, BaseContext* context, BaseObjectDescription* arg)
     {
         return core::objectmodel::BaseObject::create(obj, context, arg);
     }


     /////////////////////////// Inherited from BaseObject //////////////////////////////////////////
     virtual void parse(BaseObjectDescription* arg) override;
     virtual void init() override ;
     virtual void reinit() override {}
     virtual void doUpdate() override {}
     virtual void bwdInit() override ;
     void draw(const VisualParams * vparams) override ;


     /////////////////////////// Methods of WireRestShape  //////////////////////////////////////////

     /// For coils: a part of the coil instrument can be brokenIn2  (by default the point of release is the end of the straight length)
     Real getReleaseCurvAbs() const {return d_straightLength.getValue();}

     /// This function is called by the force field to evaluate the rest position of each beam
     void getRestTransformOnX(Transform &global_H_local, const Real &x);

     /// This function gives the Young modulus and Poisson's coefficient of the beam depending on the beam position
     void getYoungModulusAtX(const Real& x_curv, Real& youngModulus, Real& cPoisson);

     /// This function gives the mass density and the BeamSection data depending on the beam position
     void getInterpolationParam(const Real& x_curv, Real &_rho, Real &_A, Real &_Iy , Real &_Iz, Real &_Asy, Real &_Asz, Real &_J);

     /**
      * This function provides a vector with the curviliar abscissa of the noticeable point(s)
      * and the minimum density (number of points) between them
      */
     void getSamplingParameters(vector<Real>& xP_noticeable, vector<int>& nbP_density) const ;


     /// Functions enabling to load and use a geometry given from OBJ external file
     void initRestConfig();
     void getRestPosNonProcedural(Real& abs, Coord &p);
     void computeOrientation(const Vec3& AB, const Quat<Real>& Q, Quat<Real> &result);
     void initFromLoader();
     bool checkTopology();

     Real getLength() ;
     void getCollisionSampling(Real &dx, const Real &x_curv) ;
     void getNumberOfCollisionSegment(Real &dx, unsigned int &numLines) ;

     //TODO(dmarchal 2017-05-17) Please specify who and when it will be done either a time after wich
     //we can remove the todo.
     // todo => topological change !
     void releaseWirePart();

     void rotateFrameForAlignX(const Quat<Real> &input, Vec3 &x, Quat<Real> &output);


protected:
     /// Analitical creation of wire shape...
     Data<bool> d_isAProceduralShape;
     Data<Real> d_nonProceduralScale;
     Data<Real> d_length;
     Data<Real> d_straightLength;
     Data<Real> d_spireDiameter;
     Data<Real> d_spireHeight;
     Data<vector<int> > d_density;
     Data<vector<Real> > d_keyPoints;
     Data< int > d_numEdges;
     Data<vector<int> > d_numEdgesCollis;

     /// User Data about the Young modulus
     Data<Real> d_poissonRatio;
     Data<Real> d_youngModulus1;
     Data<Real> d_youngModulus2;

     /// Radius
     Data<Real> d_radius1;
     Data<Real> d_radius2;
     Data<Real> d_innerRadius1;
     Data<Real> d_innerRadius2;

     Data<Real> d_massDensity1;
     Data<Real> d_massDensity2;

     /// broken in 2 case
     Data<bool> d_brokenIn2;
     Data<bool>	d_drawRestShape;

     /// Data required for the File loading
     vector<Vec3> 		m_localRestPositions;
     vector<Transform> 	m_localRestTransforms;
     vector<Real> 		m_curvAbs ;
     double 							m_absOfGeometry {0};

     struct BeamSection{
        double _r; 			///>radius of the section
        double _rInner; 	///>inner radius of the section if beam is hollow
        double _Iy;
        double _Iz; 		///>Iz is the cross-section moment of inertia (assuming mass ratio = 1) about the z axis;
        double _J;  		///>Polar moment of inertia (J = Iy + Iz)
        double _A; 			///> A is the cross-sectional area;
        double _Asy; 		///>_Asy is the y-direction effective shear area =  10/9 (for solid circular section) or 0 for a non-Timoshenko beam
        double _Asz; 		///>_Asz is the z-direction effective shear area;
     };
     BeamSection beamSection1;
     BeamSection beamSection2;

     TopologyContainer* _topology {nullptr} ;
     EdgeSetTopologyModifier* edgeMod {nullptr} ;
     Edge2QuadTopologicalMapping* edge2QuadMap ;
     MeshLoader* loader {nullptr};
     bool edgeSetInNode {false};
};



} // namespace _wirerestshape_

using _wirerestshape_::WireRestShape;

} // namespace sofa::component::engine
