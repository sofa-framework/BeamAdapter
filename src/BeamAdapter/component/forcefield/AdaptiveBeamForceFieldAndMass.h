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

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <BeamAdapter/config.h>
#include <BeamAdapter/component/BeamInterpolation.h>
#include <BeamAdapter/component/engine/WireRestShape.h>


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Forward declarations, see https://en.wikipedia.org/wiki/Forward_declaration
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa::component::forcefield
{

/////////////////////////////////// private namespace pattern //////////////////////////////////////
/// To avoid the lacking of names imported with with 'using' in the other's component namespace
/// you should use a private namespace and "export" only this one in the public namespace.
/// This is done at the end of this file, have a look if you are not used to this pattern.
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace _adaptivebeamforcefieldandmass_
{

using sofa::core::behavior::MultiMatrixAccessor;
using sofa::core::visual::VisualParams;
using sofa::core::MechanicalParams;

/*!
 * \class AdaptiveBeamForceFieldAndMass
 * @brief AdaptiveBeamForceFieldAndMass Class
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 */
template<class DataTypes>
class AdaptiveBeamForceFieldAndMass : virtual public sofa::core::behavior::Mass<DataTypes>, virtual public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS2(SOFA_TEMPLATE(AdaptiveBeamForceFieldAndMass,DataTypes),
               SOFA_TEMPLATE(core::behavior::Mass,DataTypes),
               SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

    using Coord = typename DataTypes::Coord;
    using VecCoord = typename DataTypes::VecCoord;
    using Real = typename Coord::value_type;

    using Deriv = typename DataTypes::Deriv;
    using VecDeriv = typename DataTypes::VecDeriv;

    using DataVecCoord = Data<VecCoord>;
    using DataVecDeriv = Data<VecDeriv>;

    using Vec3 = sofa::type::Vec<3, Real>;
    using Vec6 = sofa::type::Vec<6, Real>;
    using Vec6NoInit = sofa::type::VecNoInit<6, Real>;
    using Matrix6x6 = sofa::type::Mat<6, 6, Real>;
    using Matrix6x6NoInit = sofa::type::MatNoInit<6, 6, Real>;
    using Transform = typename sofa::defaulttype::SolidTypes<Real>::Transform;
    using SpatialVector = typename sofa::defaulttype::SolidTypes<Real>::SpatialVector;

    using BInterpolation = sofa::component::fem::BeamInterpolation<DataTypes>;
    using WireRestShape = sofa::component::engine::WireRestShape<DataTypes>;
    typedef sofa::core::behavior::Mass<DataTypes> MassT;

protected:

    /*!
     * \class BeamLocalMatrices
    * @brief BeamLocalMatrices Class
    */
    class BeamLocalMatrices{

    public:
        BeamLocalMatrices() = default;

        Matrix6x6NoInit m_K00, m_K01, m_K10, m_K11; /// stiffness Matrices
        Matrix6x6NoInit m_M00, m_M01, m_M10, m_M11; /// mass Matrices
        Matrix6x6NoInit m_A0Ref, m_A1Ref;   /// adjoint Matrices

        Real _A, _L, _Iy, _Iz, _Asy, _Asz, _J; ///< Interpolation & geometrical parameters
        Real _rho;
    };

public:
    AdaptiveBeamForceFieldAndMass( ) ;
    virtual ~AdaptiveBeamForceFieldAndMass() = default;


    /////////////////////////////////////
    /// This is inhereted from BaseObject
    /////////////////////////////////////
    void init() override ;
    void reinit() override ;
    void draw(const VisualParams* vparams) override ;


    /////////////////////////////////////
    /// Mass Interface
    /////////////////////////////////////
    void addMDx(const MechanicalParams* mparams, DataVecDeriv& f, const DataVecDeriv& dx, double factor) override;
    void addMToMatrix(const MechanicalParams *mparams, const MultiMatrixAccessor* matrix) override;
    void addMBKToMatrix(const MechanicalParams* mparams, const MultiMatrixAccessor* matrix) override;

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    void accFromF(const MechanicalParams* mparams, DataVecDeriv& , const DataVecDeriv& ) override
    {
        SOFA_UNUSED(mparams);
        msg_error()<<"accFromF can not be implemented easily: It necessitates a solver because M^-1 is not available";
    }

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    double getKineticEnergy(const MechanicalParams* mparams, const DataVecDeriv& )  const override ///< vMv/2 using dof->getV()
    {
        SOFA_UNUSED(mparams);
        msg_error() << "getKineticEnergy not yet implemented";
        return 0;
    }
    
    bool isDiagonal() const override { return false; }

    /////////////////////////////////////
    /// ForceField Interface
    /////////////////////////////////////
    void addForce(const MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;

    void addDForce(const MechanicalParams* mparams, DataVecDeriv&   datadF , const DataVecDeriv&   datadX ) override;

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    double getPotentialEnergy(const MechanicalParams* mparams, const DataVecCoord& )const override
    {
        SOFA_UNUSED(mparams);
        msg_error()<<"getPotentialEnergy not yet implemented"; 
        return 0; 
    }

    SReal getGravitationalPotentialEnergy(const core::MechanicalParams* mparams, const DataVecCoord& x, const Deriv& gravity) const override
    {
        SOFA_UNUSED(mparams);
        SOFA_UNUSED(x);
        SOFA_UNUSED(gravity);
        msg_error()<<"getGravitationalPotentialEnergy not yet implemented";
        return 0.0;
    }

    void addGravitationalForce(const core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v, const Deriv& gravity) override
    {
        // Force directly computed in addForce
        SOFA_UNUSED(mparams);
        SOFA_UNUSED(f);
        SOFA_UNUSED(x);
        SOFA_UNUSED(v);
        SOFA_UNUSED(gravity);
        return;
    }

    using sofa::core::behavior::ForceField<DataTypes>::addKToMatrix;
    void addKToMatrix(const MechanicalParams* mparams,
                      const MultiMatrixAccessor* matrix) override;

    void computeStiffness(int beam, BeamLocalMatrices& beamLocalMatrices);
    void computeMass(int beam, BeamLocalMatrices& beamMatrices);


    Data<bool> d_computeMass;               ///< if false, only compute the stiff elastic model
    Data<Real> d_massDensity;               ///< Density of the mass (usually in kg/m^3)
    Data<bool> d_useShearStressComputation; ///< if false, suppress the shear stress in the computation
    Data<bool> d_reinforceLength;           ///< if true, perform a separate computation to evaluate the elongation
    DataVecDeriv d_dataG;


    bool insertInNode( core::objectmodel::BaseNode* node ) override { return core::behavior::ForceField<DataTypes>::insertInNode(node); }
    bool removeInNode( core::objectmodel::BaseNode* node ) override { return core::behavior::ForceField<DataTypes>::removeInNode(node); }

protected :

    SingleLink<AdaptiveBeamForceFieldAndMass<DataTypes>, BInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_interpolation;
    SingleLink<AdaptiveBeamForceFieldAndMass<DataTypes>, WireRestShape, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_instrumentParameters;

    void applyMassLarge( VecDeriv& df, const VecDeriv& dx, int bIndex, Index nd0Id, Index nd1Id, const double &factor);
    void applyStiffnessLarge( VecDeriv& df, const VecDeriv& dx, int beam, Index nd0Id, Index nd1Id, const double &factor );
    void computeGravityVector();

    Vec3 m_gravity;
    type::vector<BeamLocalMatrices> m_localBeamMatrices;

private:

    void drawElement(const VisualParams* vparams, int beam,
                     Transform &global_H0_local, Transform &global_H1_local) ;
};

/// Instantiate the templates so that they are not instiated in each translation unit (see )
#if !defined(SOFA_PLUGIN_BEAMADAPTER_ADAPTIVEBEAMFORCEFIELD_CPP)
extern template class SOFA_BEAMADAPTER_API AdaptiveBeamForceFieldAndMass<sofa::defaulttype::Rigid3Types> ;
#endif

} /// namespace _adaptivebeamforcefieldandmass_


////////////////////////////////// EXPORT NAMES IN SOFA NAMESPACE //////////////////////////////////
/// 'Export' the objects defined in the private namespace into the 'public' one.
////////////////////////////////////////////////////////////////////////////////////////////////////
using _adaptivebeamforcefieldandmass_::AdaptiveBeamForceFieldAndMass ;


} /// namespace sofa::component::forcefield
