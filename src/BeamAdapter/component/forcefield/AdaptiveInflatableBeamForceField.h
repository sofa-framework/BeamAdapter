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
// C++ Implementation : BeamInterpolation / AdaptiveInflatableBeamForceField
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
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/type/vector.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/visual/VisualParams.h>

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
namespace _AdaptiveInflatableBeamForceField_
{

using sofa::type::vector;
using sofa::component::engine::WireRestShape ;
using sofa::component::fem::BeamInterpolation ;
using sofa::core::behavior::MultiMatrixAccessor ;
using sofa::core::visual::VisualParams ;
using sofa::core::behavior::Mass ;
using sofa::core::MechanicalParams ;
using sofa::type::Vec ;
using sofa::type::Mat ;
using sofa::defaulttype::Rigid3Types ;
using core::objectmodel::Data ;
using core::topology::BaseMeshTopology;
using sofa::defaulttype::SolidTypes;

/*!
 * \class AdaptiveInflatableBeamForceField
 * @brief AdaptiveInflatableBeamForceField Class
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 */
template<class DataTypes>
class AdaptiveInflatableBeamForceField : public Mass<DataTypes>
{
public:

    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveInflatableBeamForceField,DataTypes),
               SOFA_TEMPLATE(core::behavior::Mass,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::VecReal VecReal;
    typedef VecCoord Vector;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef Data<VecCoord> DataVecCoord;
    typedef Data<VecDeriv> DataVecDeriv;
    typedef BeamInterpolation<DataTypes> BInterpolation;

    typedef unsigned int Index;
    typedef BaseMeshTopology::Edge Element;
    typedef type::vector<BaseMeshTopology::Edge> VecElement;
    typedef type::vector<unsigned int> VecIndex;

    typedef typename SolidTypes<Real>::Transform Transform;
    typedef typename SolidTypes<Real>::SpatialVector SpatialVector;

    typedef Vec<3, Real> Vec3;
    typedef Vec<6, Real> Vec6;          ///< the displacement vector
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
            m_K00 = v.m_K00;
            m_K10 = v.m_K10;
            m_K11 = v.m_K11;
            m_K01 = v.m_K01;

            m_M00 = v.m_M00;
            m_M11 = v.m_M11;
            m_M01 = v.m_M01;
            m_M10 = v.m_M10;

            m_A0Ref = v.m_A0Ref;
            m_A1Ref = v.m_A1Ref;
        }
        ~BeamLocalMatrices(){}

        Matrix6x6 m_K00, m_K01, m_K10, m_K11; /// stiffness Matrices
        Matrix6x6 m_M00, m_M01, m_M10, m_M11; /// mass Matrices
        Matrix6x6 m_A0Ref, m_A1Ref;   /// adjoint Matrices
    };

public:

    AdaptiveInflatableBeamForceField( ) ;
    virtual ~AdaptiveInflatableBeamForceField() = default;


    /////////////////////////////////////
    /// This is inhereted from BaseObject
    /////////////////////////////////////
    virtual void init() override ;
    virtual void reinit() override ;
    virtual void draw(const VisualParams* vparams) override ;


    /////////////////////////////////////
    /// Mass Interface
    /////////////////////////////////////
    virtual void addMDx(const MechanicalParams* mparams, DataVecDeriv& f, const DataVecDeriv& dx, SReal factor) override;
    virtual void addMToMatrix(const MechanicalParams *mparams, const MultiMatrixAccessor* matrix) override;
    virtual void addMBKToMatrix(const MechanicalParams* mparams, const MultiMatrixAccessor* matrix) override;

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    virtual void accFromF(const MechanicalParams* mparams, DataVecDeriv& , const DataVecDeriv& ) override
    {
        SOFA_UNUSED(mparams);
        msg_error() << "accFromF can not be implemented easily: It necessitates a solver because M^-1 is not available";
    }

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    virtual SReal getKineticEnergy(const MechanicalParams* mparams, const DataVecDeriv& )  const override ///< vMv/2 using dof->getV()
    {
        SOFA_UNUSED(mparams);
        msg_error() << "getKineticEnergy not yet implemented";
        return 0;
    }

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    virtual void addGravityToV(const MechanicalParams* mparams, DataVecDeriv& ) override
    {
        SOFA_UNUSED(mparams);
        msg_error() << "addGravityToV not implemented yet";
    }


    bool isDiagonal() const override { return false; }

    /////////////////////////////////////
    /// ForceField Interface
    /////////////////////////////////////
    virtual void addForce(const MechanicalParams* mparams,
                          DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;

    virtual void addDForce(const MechanicalParams* mparams,
                           DataVecDeriv&   datadF , const DataVecDeriv&   datadX ) override;

    //TODO(dmarchal 2017-05-17) So what do we do ? For who is this message intended for ? How can we make this code "more" manageable.
    virtual SReal getPotentialEnergy(const MechanicalParams* mparams, const DataVecCoord& ) const override
    {
        SOFA_UNUSED(mparams);
        msg_error() << "getPotentialEnergy not yet implemented";
        return 0_sreal;
    }

    using sofa::core::behavior::ForceField<DataTypes>::addKToMatrix;
    void addKToMatrix(const MechanicalParams* mparams,
                      const MultiMatrixAccessor* matrix) override;

    void computeStiffness(int beam, BeamLocalMatrices& beamLocalMatrices);
    void computeMass(int beam, BeamLocalMatrices& beamMatrices);


    Data<bool> d_computeMass;               ///< if false, only compute the stiff elastic model
    Data<Real> d_massDensity;               ///< Density of the mass (usually in kg/m^3)
    Data<bool> d_reinforceLength;           ///< if true, perform a separate computation to evaluate the elongation
    Data<Real> d_pressure;                  ///< Pressure applied inside the inflatable beams
    DataVecDeriv d_dataG;

protected :

    SingleLink<AdaptiveInflatableBeamForceField<DataTypes>, BInterpolation          , BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_interpolation;
    SingleLink<AdaptiveInflatableBeamForceField<DataTypes>, WireRestShape<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_instrumentParameters;

    void applyMassLarge( VecDeriv& df, const VecDeriv& dx, int bIndex, Index nd0Id, Index nd1Id, SReal factor);
    void applyStiffnessLarge( VecDeriv& df, const VecDeriv& dx, int beam, Index nd0Id, Index nd1Id, SReal factor );
    void computeGravityVector();

    Vec3 m_gravity;
    type::vector<BeamLocalMatrices> m_localBeamMatrices;

    using Mass<DataTypes>::getContext;
    using Mass<DataTypes>::mstate;

private:

    void drawElement(const VisualParams* vparams, int beam,
                     Transform &global_H0_local, Transform &global_H1_local) ;
};

/// Instantiate the templates so that they are not instiated in each translation unit (see )
#if !defined(SOFA_PLUGIN_BEAMADAPTER_ADAPTIVEINFLATABLEBEAMFORCEFIELD_CPP)
extern template class SOFA_BEAMADAPTER_API AdaptiveInflatableBeamForceField<Rigid3Types> ;
#endif

} /// namespace _AdaptiveInflatableBeamForceField_


////////////////////////////////// EXPORT NAMES IN SOFA NAMESPACE //////////////////////////////////
/// 'Export' the objects defined in the private namespace into the 'public' one.
////////////////////////////////////////////////////////////////////////////////////////////////////
using _AdaptiveInflatableBeamForceField_::AdaptiveInflatableBeamForceField ;


} /// namespace sofa::component::forcefield
