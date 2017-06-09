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
#ifndef SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEAMFORCEFIELDANDMASS_H
#define SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEAMFORCEFIELDANDMASS_H

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/visual/VisualParams.h>

#include "../initBeamAdapter.h"
#include "../BeamInterpolation.h"
#include "../engine/WireRestShape.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Forward declarations, see https://en.wikipedia.org/wiki/Forward_declaration
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa
{

namespace component
{

namespace forcefield
{

/////////////////////////////////// private namespace pattern //////////////////////////////////////
/// To avoid the lacking of names imported with with 'using' in the other's component namespace
/// you should use a private namespace and "export" only this one in the public namespace.
/// This is done at the end of this file, have a look if you are not used to this pattern.
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace _adaptivebeamforcefieldandmass_
{

using sofa::helper::vector;
using sofa::component::engine::WireRestShape ;
using sofa::component::fem::BeamInterpolation ;
using sofa::core::behavior::MultiMatrixAccessor ;
using sofa::core::visual::VisualParams ;
using sofa::core::behavior::Mass ;
using sofa::core::MechanicalParams ;
using sofa::defaulttype::Vec ;
using sofa::defaulttype::Mat ;

/*!
 * \class AdaptiveBeamForceFieldAndMass
 * @brief AdaptiveBeamForceFieldAndMass Class
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 */
template<class DataTypes>
class AdaptiveBeamForceFieldAndMass : public Mass<DataTypes>
{
public:

    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBeamForceFieldAndMass,DataTypes),
               SOFA_TEMPLATE(core::behavior::Mass,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::VecReal VecReal;
    typedef VecCoord Vector;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;
    typedef BeamInterpolation<DataTypes> BInterpolation;

    typedef unsigned int Index;
    typedef core::topology::BaseMeshTopology::Edge Element;
    typedef sofa::helper::vector<core::topology::BaseMeshTopology::Edge> VecElement;
    typedef helper::vector<unsigned int> VecIndex;

    typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
    typedef typename sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;

    /// remove ?
    typedef Vec<12, Real> Displacement;        ///< the displacement vector
    typedef Mat<3, 3, Real> Transformation; ///< matrix for rigid transformations like rotations
    typedef Mat<12, 12, Real> StiffnessMatrix;

    /// replace by:

    typedef Vec<3, Real> Vec3;
    typedef Vec<6, Real> Vec6;        ///< the displacement vector
    typedef Mat<6, 6, Real> Matrix6x6;

    /*!
     * \class BeamLocalMatrices
     * @brief BeamLocalMatrices Class
     */
protected:
    class BeamLocalMatrices{

    public:
        BeamLocalMatrices(){}
        BeamLocalMatrices(const BeamLocalMatrices &v)
        {
            this->k_loc00 = v.k_loc00;
            this->k_loc10 = v.k_loc10;
            this->k_loc11 = v.k_loc11;
            this->k_loc01 = v.k_loc01;

            this->m_loc00 = v.m_loc00;
            this->m_loc11 = v.m_loc11;
            this->m_loc01 = v.m_loc01;
            this->m_loc10 = v.m_loc10;

            this->loc0_Ad_ref = v.loc0_Ad_ref;
            this->loc1_Ad_ref = v.loc1_Ad_ref;
        }
        ~BeamLocalMatrices(){}
        // stiffness Matrices
        Matrix6x6 k_loc00, k_loc01, k_loc10, k_loc11;
        // mass Matrices
        Matrix6x6 m_loc00, m_loc01, m_loc10, m_loc11;
        // adjoint Matrices
        Matrix6x6 loc0_Ad_ref, loc1_Ad_ref;
    };

public:
    AdaptiveBeamForceFieldAndMass( ) ;
    virtual ~AdaptiveBeamForceFieldAndMass() ;

    /// This is inhereted from BaseObject
    virtual void init() override ;
    virtual void reinit() override ;
    virtual void draw(const VisualParams* vparams) override ;

    /// Mass Interface
    virtual  void addMDx(const MechanicalParams* mparams, DataVecDeriv& f, const DataVecDeriv& dx, double factor);
    virtual void addMToMatrix(const MechanicalParams *mparams, const MultiMatrixAccessor* matrix);
    virtual void addMBKToMatrix(const MechanicalParams* mparams, const MultiMatrixAccessor* matrix);

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    virtual  void accFromF(const MechanicalParams* mparams, DataVecDeriv& , const DataVecDeriv& )
    {serr<<"accFromF can not be implemented easily: It necessitates a solver because M^-1 is not available"<<sendl;}

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    virtual double getKineticEnergy(const MechanicalParams* mparams, const DataVecDeriv& )  const ///< vMv/2 using dof->getV()
    {serr<<"getKineticEnergy not yet implemented"<<sendl;return 0;}

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    virtual void addGravityToV(const MechanicalParams* mparams, DataVecDeriv& )
    {serr<<"addGravityToV not implemented yet"<<sendl;}

    /// ForceField Interface
    virtual void addForce (const MechanicalParams* mparams,
                           DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);

    virtual void  addDForce(const MechanicalParams* mparams,
                            DataVecDeriv&   datadF , const DataVecDeriv&   datadX );

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    virtual double getPotentialEnergy(const MechanicalParams* mparams, const DataVecCoord& )
    const {serr<<"getPotentialEnergy not yet implemented"<<sendl; return 0; }

    void addKToMatrix(const MechanicalParams* mparams,
                      const MultiMatrixAccessor* matrix);

    void computeStiffness(int beam, BeamLocalMatrices& beamMatrices);
    void computeMass(int beam, BeamLocalMatrices& beamMatrices);


    Data<bool> d_timoshenko; ///< use Timoshenko beam (non-null section shear area)
    Data<bool> d_computeMass; ///< if false, only compute the stiff elastic model
    Data<Real> d_massDensity; ///< Density of the mass (usually in kg/m^3)
    Data<bool> d_shearStressComputation; ///< if false, suppress the shear stress in the computation
    Data<bool> d_reinforceLength; ///<if true, perform a separate computation to evaluate the elongation
    DataVecDeriv d_dataG;

protected :
    /// pointer to the interpolation
    SingleLink<AdaptiveBeamForceFieldAndMass<DataTypes>, BInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_interpolation;

    /// pointer to the WireRestShape
    SingleLink<AdaptiveBeamForceFieldAndMass<DataTypes>, WireRestShape<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_instrumentParameters;

    void applyMassLarge( VecDeriv& df, const VecDeriv& dx, int bIndex, Index nd0Id, Index nd1Id, const double &factor);
    void applyStiffnessLarge( VecDeriv& df, const VecDeriv& dx, int beam, Index nd0Id, Index nd1Id, const double &factor );

    /// compute the gravity vector (if necessary)
    void computeGravityVector();
    Vec3 m_gravity;

    vector<BeamLocalMatrices> d_localBeamMatrices;

private:
    void drawElement(const VisualParams* vparams, int beam,
                     Transform &global_H0_local, Transform &global_H1_local) ;
};

/// Instantiate the templates so that they are not instiated in each translation unit (see )
#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_PLUGIN_BEAMADAPTER_ADAPTVEBEAMFORCEFIELD_CPP)
#ifdef SOFA_WITH_FLOAT
extern template class SOFA_BEAMADAPTER_API AdaptiveBeamForceFieldAndMass<sofa::defaulttype::Rigid3fTypes> ;
#endif
#ifdef SOFA_WITH_DOUBLE
extern template class SOFA_BEAMADAPTER_API AdaptiveBeamForceFieldAndMass<sofa::defaulttype::Rigid3dTypes> ;
#endif
#endif

} /// namespace _adaptivebeamforcefieldandmass_



////////////////////////////////// EXPORT NAMES IN SOFA NAMESPACE //////////////////////////////////
/// 'Export' the objects defined in the private namespace into the 'public' one.
////////////////////////////////////////////////////////////////////////////////////////////////////
using _adaptivebeamforcefieldandmass_::AdaptiveBeamForceFieldAndMass ;


} /// namespace forcefield

} /// namespace component

} /// namespace sofa



#endif  /* SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEAMFORCEFIELDANDMASS_H */
