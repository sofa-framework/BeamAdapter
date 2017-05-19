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
// C++ Implementation : WireRestShape
//
// Description:
//
// Contributors:
//   - Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef SOFA_COMPONENT_ENGINE_WIRERESTSHAPE_H
#define SOFA_COMPONENT_ENGINE_WIRERESTSHAPE_H

//TODO(dmarchal 2017-05-17) Do we really need so much include ?
#include "initBeamAdapter.h"
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/DataEngine.h>
#include <SofaBaseTopology/EdgeSetTopologyModifier.h>
#include <SofaBaseTopology/EdgeSetGeometryAlgorithms.h>
#include <SofaTopologyMapping/Edge2QuadTopologicalMapping.h>
#include <SofaLoader/MeshObjLoader.h>

namespace sofa
{

namespace component
{

namespace engine
{

/*!
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
    typedef sofa::defaulttype::Vec<3, Real> Vec3;
    typedef sofa::defaulttype::Vec<2, Real> Vec2;
    typedef sofa::defaulttype::Quat Quat;
    typedef typename sofa::helper::vector<Vec2>::iterator vecIt;

    /*!
     * @brief Default Constructor.
     */
     WireRestShape() ;

     /*!
      * @brief Default Destructor.
      */
     ~WireRestShape(){}

     void rotateFrameForAlignX(const sofa::defaulttype::Quat &input, Vec3 &x, sofa::defaulttype::Quat &output);

     /// These are virtual in-herited from BaseObeject,
     /// use the override keyword to make things clear http://en.cppreference.com/w/cpp/language/override
     void init() override ;
     void reinit() override{ }
     void update() override { }
     void bwdInit() override ;
     void draw(const core::visual::VisualParams* vparams) override ;


     /*!
      * For coils: a part of the coil instrument can be brokenIn2  (by default the point of release is the end of the straight length)
      */
     virtual Real getReleaseCurvAbs(){return d_straightLength.getValue();}

     //TODO(dmarchal 2017-05-17) Please specify who and when it will be done either a time after wich
     //we can remove the todo.
     // todo => topological change !
     virtual void releaseWirePart();

     /*!
      * This function is called by the force field to evaluate the rest position of each beam
      */
     virtual void getRestTransformOnX(Transform &global_H_local, const Real &x);

     /*!
      * This function provides a vector with the curviliar abscissa of the noticeable point(s)
      * and the minimum density (number of points) between them
      */
     virtual void getSamplingParameters(helper::vector<Real>& xP_noticeable, helper::vector<int>& nbP_density);

     /*!
      * This function gives the Young modulus and Poisson's coefficient of the beam depending on the beam position
      */
     virtual void getYoungModulusAtX(const Real& x_curv, Real& youngModulus, Real& cPoisson);

     /*!
      * this function gives the mass density and the BeamSection data depending on the beam position
      */
     void getInterpolationParam(const Real& x_curv, Real &_rho, Real &_A, Real &_Iy , Real &_Iz, Real &_Asy, Real &_Asz, Real &_J);

     /*!
      * Functions enabling to load and use a geometry given from OBJ external file
      */
     void initRestConfig();
     void getRestPosNonProcedural(Real& abs, Coord &p);
     void computeOrientation(const Vec3& AB, const Quat& Q, Quat &result);
     void initFromLoader();
     bool checkTopology();

     virtual Real getLength() ;
     virtual void getCollisionSampling(Real &dx, const Real &x_curv) ;
     virtual void getNumberOfCollisionSegment(Real &dx, unsigned int &numLines) ;

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
     /// Analitical creation of wire shape...
     Data<bool> d_procedural;
     Data<Real> d_nonProceduralScale;
     Data<Real> d_length;
     Data<Real> d_straightLength;
     Data<Real> d_spireDiameter;
     Data<Real> d_spireHeight;
     Data< sofa::helper::vector<int> > d_density;
     Data< sofa::helper::vector<Real> > d_keyPoints;
     Data< int > d_numEdges;
     Data< sofa::helper::vector<int> > d_numEdgesCollis;

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
     sofa::helper::vector<Vec3> 		m_localRestPositions;
     sofa::helper::vector<Transform> 	m_localRestTransforms;
     sofa::helper::vector<Real> 		m_curvAbs ;
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

     sofa::core::topology::TopologyContainer* _topology {nullptr} ;
     sofa::component::topology::EdgeSetGeometryAlgorithms<DataTypes>* edgeGeo {nullptr};
     sofa::component::topology::EdgeSetTopologyModifier* edgeMod {nullptr};
     sofa::component::topology::Edge2QuadTopologicalMapping* edge2QuadMap {nullptr};
     sofa::core::loader::MeshLoader* loader {nullptr};
     bool edgeSetInNode {false};
};


//TODO(dmarchal 2017-05-17 Use extern template here to reduce code bloat.

} // namespace engine

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_ENGINE_WIRERESTSHAPE_H */
