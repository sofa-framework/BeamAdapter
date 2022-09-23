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
#include <BeamAdapter/utils/BeamSection.h>
#include <BeamAdapter/component/model/WireSectionMaterial.h>

#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyContainer.h>
#include <sofa/core/loader/MeshLoader.h>

namespace sofa::component::engine
{

namespace _wirerestshape_
{

using sofa::core::topology::TopologyContainer;
using sofa::core::loader::MeshLoader;

using namespace sofa::beamadapter;

/**
 * \class WireRestShape
 * \brief Describe the shape functions on multiple segments
 *  
 *  Describe the full shape of a Wire with a given length and radius. The wire is discretized by a set of beams (given by the keyPoints and the relatives Beam density)
 *  This component compute the beam discretization and the shape functions on multiple segments using curvilinear abscissa.
 */
template <class DataTypes>
class WireRestShape : public core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(WireRestShape, core::objectmodel::BaseObject);

    using Coord = typename DataTypes::Coord;
    using Real = typename Coord::value_type;
    using Transform = typename sofa::defaulttype::SolidTypes<Real>::Transform;
    using Vec3 = sofa::type::Vec<3, Real>;
    using Quat = sofa::type::Quat<Real>;
   
    using BeamSection = sofa::beamadapter::BeamSection;

    /**
     * @brief Default Constructor.
     */
     WireRestShape() ;

     /*!
      * @brief Default Destructor.
      */
     virtual ~WireRestShape() = default;

     /////////////////////////// Inherited from BaseObject //////////////////////////////////////////
     void parse(core::objectmodel::BaseObjectDescription* arg) override;
     void init() override ;
       
     void draw(const core::visual::VisualParams * vparams) override ;


     /////////////////////////// Methods of WireRestShape  //////////////////////////////////////////

     /// This function is called by the force field to evaluate the rest position of each beam
     void getRestTransformOnX(Transform &global_H_local, const Real &x);

     /// This function gives the Young modulus and Poisson's coefficient of the beam depending on the beam position
     void getYoungModulusAtX(const Real& x_curv, Real& youngModulus, Real& cPoisson) const;

     /// This function gives the mass density and the BeamSection data depending on the beam position
     void getInterpolationParam(const Real& x_curv, Real &_rho, Real &_A, Real &_Iy , Real &_Iz, Real &_Asy, Real &_Asz, Real &_J) const;

     /**
      * This function provides a type::vector with the curviliar abscissa of the noticeable point(s) 
      * and the minimum density (number of points) between them. (Nb. nbP_density.size() == xP_noticeable.size() - 1)
      */
     void getSamplingParameters(type::vector<Real>& xP_noticeable, type::vector<int>& nbP_density) const ;


     /// Functions enabling to load and use a geometry given from OBJ external file
     void initRestConfig();
     void getRestPosNonProcedural(Real& abs, Coord &p);
     void computeOrientation(const Vec3& AB, const Quat& Q, Quat &result);     
     void initFromLoader();
     bool checkTopology();

     [[nodiscard]] bool fillTopology();
     Real getLength() ;
     void getCollisionSampling(Real &dx, const Real &x_curv) ;
     void getNumberOfCollisionSegment(Real &dx, unsigned int &numLines) ;

     void rotateFrameForAlignX(const Quat &input, Vec3 &x, Quat &output);


     /////////////////////////// Deprecated Methods  ////////////////////////////////////////// 

     /// For coils: a part of the coil instrument can be brokenIn2  (by default the point of release is the end of the straight length)
     Real getReleaseCurvAbs() const {
         msg_warning() << "Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06";
         return 0.0;
     }

     void releaseWirePart() {
         msg_warning() << "Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06";
     }

public:
     /// Analitical creation of wire shape...
     Data<bool> d_isAProceduralShape;
     Data<Real> d_nonProceduralScale;
     Data<Real> d_length;
     Data<Real> d_straightLength;
     Data<Real> d_spireDiameter;
     Data<Real> d_spireHeight;
     Data<type::vector<int> > d_density;
     Data<type::vector<Real> > d_keyPoints;
     Data< int > d_numEdges;
     Data<type::vector<int> > d_numEdgesCollis;

     /// broken in 2 case
     Data<bool>	d_drawRestShape;

     /// Link to be set to the topology container in the component graph.
     SingleLink<WireRestShape<DataTypes>, WireSectionMaterial, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_sectionMaterial1;

     /// Link to be set to the topology container in the component graph.
     SingleLink<WireRestShape<DataTypes>, WireSectionMaterial, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_sectionMaterial2;

private:
     /// Data required for the File loading
     type::vector<Vec3> 		m_localRestPositions;
     type::vector<Transform> 	m_localRestTransforms;
     type::vector<Real> 		m_curvAbs ;
     double 							m_absOfGeometry {0};
     
     /// Link to be set to the topology container in the component graph.
     SingleLink<WireRestShape<DataTypes>, TopologyContainer, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_topology;     
     /// Pointer to the topology container, should be set using @sa l_topology, otherwise will search for one in current Node.
     TopologyContainer* _topology{ nullptr }; 

     /// Link to be set to the topology container in the component graph.
     SingleLink<WireRestShape<DataTypes>, MeshLoader, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_loader;     
     /// Pointer to the MeshLoader, should be set using @sa l_loader, otherwise will search for one in current Node.
     MeshLoader* loader{ nullptr };
};


#if !defined(SOFA_PLUGIN_BEAMADAPTER_WIRERESTSHAPE_CPP)
extern template class SOFA_BEAMADAPTER_API WireRestShape<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace _wirerestshape_

using _wirerestshape_::WireRestShape;

} // namespace sofa::component::engine
