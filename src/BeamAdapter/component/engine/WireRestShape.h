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
#include <BeamAdapter/component/model/BaseRodSectionMaterial.h>

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
 * Describe the full shape of a Wire with a given set of @sa BaseRodSectionMaterial. The wire is discretized by a set of beams (given by the keyPoints and the relatives Beam density)
 * @sa d_keyPoints and @d_density are computed by method @sa initLengths using the set of rod sections description.
 * This component compute the beam discretization and the shape functions on multiple segments using curvilinear abscissa.
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
   
    /**
     * @brief Default Constructor.
     */
     WireRestShape() ;

     /*!
      * @brief Default Destructor.
      */
     virtual ~WireRestShape() = default;

     /////////////////////////// Inherited from BaseObject //////////////////////////////////////////
     void init() override ;


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
     void computeOrientation(const Vec3& AB, const Quat& Q, Quat &result);     
     
     
     Real getLength() ;
     void getCollisionSampling(Real &dx, const Real &x_curv);
     void getNumberOfCollisionSegment(Real &dx, unsigned int &numLines) ;



     /////////////////////////// Deprecated Methods  ////////////////////////////////////////// 

     /// For coils: a part of the coil instrument can be brokenIn2  (by default the point of release is the end of the straight length)
     [[deprecated("Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06")]]
     Real getReleaseCurvAbs() const {
         msg_warning() << "Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06";
         return 0.0;
     }

     [[deprecated("Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06")]]
     void releaseWirePart() {
         msg_warning() << "Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06";
     }

protected:
    /// Internal method to init Lengths vector @sa d_keyPoints using the length of each materials @sa l_sectionMaterials.
    void initLengths();
    /// Internal method to init Edge Topology @sa _topology using the list of materials @sa l_sectionMaterials. Returns false if init can't be performed.
    bool initTopology();


public:
     Data<type::vector<int> > d_density;
     Data<type::vector<Real> > d_keyPoints;
     
     /// Vector or links to the Wire section material. The order of the linked material will define the WireShape structure.
     MultiLink<WireRestShape<DataTypes>, BaseRodSectionMaterial<DataTypes>, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_sectionMaterials;

private:
     /// Link to be set to the topology container in the component graph.
     SingleLink<WireRestShape<DataTypes>, TopologyContainer, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_topology;     
     /// Pointer to the topology container, should be set using @sa l_topology, otherwise will search for one in current Node.
     TopologyContainer* _topology{ nullptr }; 
};


#if !defined(SOFA_PLUGIN_BEAMADAPTER_WIRERESTSHAPE_CPP)
extern template class SOFA_BEAMADAPTER_API WireRestShape<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace _wirerestshape_

using _wirerestshape_::WireRestShape;

} // namespace sofa::component::engine
