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

#include <BeamAdapter/config.h>

#include <sofa/type/Transform.h>
#include <sofa/type/SpatialVector.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <BeamAdapter/component/BaseBeamInterpolation.h>
#include <BeamAdapter/component/engine/WireRestShape.h>


namespace beamadapter
{

/*!
 * \class AdaptiveBeamForceFieldAndMass
 * @brief AdaptiveBeamForceFieldAndMass Class
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 */
template<class DataTypes>
class AdaptiveBeamForceFieldAndMass : public core::behavior::Mass<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBeamForceFieldAndMass,DataTypes),
               SOFA_TEMPLATE(core::behavior::Mass,DataTypes));

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
    using Transform = sofa::type::Transform<Real>;
    using SpatialVector = sofa::type::SpatialVector<Real>;

    using BInterpolation = BaseBeamInterpolation<DataTypes>;
    using core::behavior::Mass<DataTypes>::mstate;

protected:

    /*!
     * \class BeamLocalMatrices
    * @brief BeamLocalMatrices Class
    */
    class BeamLocalMatrices
    {
    public:
        BeamLocalMatrices() = default;

        Matrix6x6NoInit m_K00, m_K01, m_K10, m_K11; /// stiffness Matrices
        Matrix6x6NoInit m_M00, m_M01, m_M10, m_M11; /// mass Matrices
        Matrix6x6NoInit m_A0Ref, m_A1Ref;   /// adjoint Matrices

        Real _A, _L, _Iy, _Iz, _Asy, _Asz, _J; ///< Interpolation & geometrical parameters
        Real _rho;
        Real _youngM;
        Real _cPoisson;
    };

public:
    AdaptiveBeamForceFieldAndMass( ) ;
    virtual ~AdaptiveBeamForceFieldAndMass() = default;


    /////////////////////////////////////
    /// This is inhereted from BaseObject
    /////////////////////////////////////
    void init() override ;
    void reinit() override ;
    void draw(const sofa::core::visual::VisualParams* vparams) override ;


    /////////////////////////////////////
    /// Mass Interface
    /////////////////////////////////////
    void addMDx(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecDeriv& dx, SReal factor) override;
    void addMToMatrix(const sofa::core::MechanicalParams *mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix) override;
    void addMBKToMatrix(const sofa::core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix) override;

    void buildMassMatrix(sofa::core::behavior::MassMatrixAccumulator* matrices) override;
    void buildStiffnessMatrix(core::behavior::StiffnessMatrix* matrix) override;
    void buildDampingMatrix(core::behavior::DampingMatrix* matrices) override;

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    void accFromF(const sofa::core::MechanicalParams* mparams, DataVecDeriv& , const DataVecDeriv& ) override
    {
        SOFA_UNUSED(mparams);
        msg_error()<<"accFromF can not be implemented easily: It necessitates a solver because M^-1 is not available";
    }

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    SReal getKineticEnergy(const sofa::core::MechanicalParams* mparams, const DataVecDeriv& )  const override ///< vMv/2 using dof->getV()
    {
        SOFA_UNUSED(mparams);
        msg_error() << "getKineticEnergy not yet implemented";
        return 0;
    }

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    void addGravityToV(const sofa::core::MechanicalParams* mparams, DataVecDeriv& ) override
    {
        SOFA_UNUSED(mparams);
        msg_error() << "addGravityToV not implemented yet";
    }
    
    bool isDiagonal() const override { return false; }

    /////////////////////////////////////
    /// ForceField Interface
    /////////////////////////////////////
    void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;

    void addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv&   datadF , const DataVecDeriv&   datadX ) override;

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    SReal getPotentialEnergy(const sofa::core::MechanicalParams* mparams, const DataVecCoord& )const override
    {
        SOFA_UNUSED(mparams);
        msg_error()<<"getPotentialEnergy not yet implemented"; 
        return 0_sreal; 
    }

    using sofa::core::behavior::ForceField<DataTypes>::addKToMatrix;
    void addKToMatrix(const sofa::core::MechanicalParams* mparams,
                      const sofa::core::behavior::MultiMatrixAccessor* matrix) override;

    void computeStiffness(const sofa::Index beamID, BeamLocalMatrices& beamLocalMatrices);
    void computeMass(const sofa::Index beamID, BeamLocalMatrices& beamMatrices);

    Data<bool> d_computeMass;               ///< if false, only compute the stiff elastic model
    Real m_defaultMassDensity;
    Data<type::vector<Real>> d_massDensity; ///< Density of the mass
    Data<bool> d_useShearStressComputation; ///< if false, suppress the shear stress in the computation
    Data<bool> d_reinforceLength;           ///< if true, perform a separate computation to evaluate the elongation

protected :
    SingleLink<AdaptiveBeamForceFieldAndMass<DataTypes>, BInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_interpolation;

    void applyMassLarge( VecDeriv& df, const sofa::Index beamID, const sofa::Index nd0Id, const sofa::Index nd1Id, const SReal factor);
    void applyStiffnessLarge( VecDeriv& df, const VecDeriv& dx, const sofa::Index beamID, const sofa::Index nd0Id, const sofa::Index nd1Id, const SReal factor );
    void computeGravityVector();

private:
    type::vector<BeamLocalMatrices> m_localBeamMatrices;
    Vec6 m_gravity;

    void drawElement(const sofa::core::visual::VisualParams* vparams, const sofa::Index beamID,
                     const Transform &global_H0_local, const Transform &global_H1_local) ;
};

/// Instantiate the templates so that they are not instiated in each translation unit (see )
#if !defined(SOFA_PLUGIN_BEAMADAPTER_ADAPTIVEBEAMFORCEFIELD_CPP)
extern template class SOFA_BEAMADAPTER_API AdaptiveBeamForceFieldAndMass<sofa::defaulttype::Rigid3Types> ;
#endif

} // namespace beamadapter
