/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
//
// C++ Implementation : AdaptivePlasticBeamForceField
//
// Description: Implementation of plasticity over the AdaptiveBeamForceFieldAndMass
// force field interface.
//
//
// Author: Camille Krewcuns, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#pragma once

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/visual/VisualParams.h>

#include <sofa/helper/fixed_array.h>

#include "../BeamPlasticInterpolation.h"
#include "AdaptiveBeamForceFieldAndMass.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa::plugin::beamadapter::component::forcefield
{

namespace _adaptiveplasticbeamforcefield_
{

using sofa::helper::vector;
using sofa::core::visual::VisualParams;
using sofa::core::MechanicalParams;
using sofa::defaulttype::Vec;
using sofa::defaulttype::Mat;
using sofa::core::behavior::MultiMatrixAccessor;
using sofa::defaulttype::Rigid3Types;
using sofa::defaulttype::SolidTypes;

using sofa::plugin::beamadapter::component::forcefield::_beamplasticinterpolation_::BeamPlasticInterpolation;
using sofa::plugin::beamadapter::component::forcefield::_beamplasticinterpolation_::BeamGeometry;
using sofa::component::forcefield::_adaptivebeamforcefieldandmass_::AdaptiveBeamForceFieldAndMass;

/*!
 * \class AdaptiveInflatableBeamForceField
 * @brief AdaptiveInflatableBeamForceField Class
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 */
template<class DataTypes>
class AdaptivePlasticBeamForceField : public AdaptiveBeamForceFieldAndMass<DataTypes>
{
public:

    SOFA_CLASS(SOFA_TEMPLATE(AdaptivePlasticBeamForceField, DataTypes),
        SOFA_TEMPLATE(AdaptiveBeamForceFieldAndMass, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename Coord::value_type Real;
    typedef Data<VecCoord> DataVecCoord;
    typedef Data<VecDeriv> DataVecDeriv;
    typedef AdaptiveBeamForceFieldAndMass<DataTypes>::BeamLocalMatrices BeamLocalMatrices;
    typedef BeamPlasticInterpolation<DataTypes> BPInterpolation;

    typedef typename SolidTypes<Real>::Transform Transform;
    typedef typename SolidTypes<Real>::SpatialVector SpatialVector;

    typedef Vec<3, Real> Vec3;
    typedef Vec<6, Real> Vec6;
    typedef Vec<9, Real> Vec9;          ///< Second-order tensor in vector notation
    typedef Mat<9, 9, Real> Matrix9x9;  ///< Fourth-order tensor in vector notation

    /** \enum MechanicalState
     *  \brief Types of mechanical state associated with the (Gauss) integration
     *  points. The POSTPLASTIC state corresponds to points which underwent plastic
     *  deformation, but on which constraints were released so that the plasticity
     *  process stopped.
     */
    enum MechanicalState {
        ELASTIC = 0,
        PLASTIC = 1,
        POSTPLASTIC = 2,
    };

    /*!
     * \class BeamPlasticVariables
     * @brief BeamPlasticVariables Class
     */
protected:

    class BeamPlasticVariables
    {
    public:
        BeamPlasticVariables() {}
        BeamPlasticVariables(const BeamPlasticVariables& v)
        {
            m_Kt00 = v.m_Kt00;
            m_Kt10 = v.m_Kt10;
            m_Kt11 = v.m_Kt11;
            m_Kt01 = v.m_Kt01;
        }
        ~BeamPlasticVariables() {}

        Matrix6x6 m_Kt00, m_Kt01, m_Kt10, m_Kt11; /// Local tangent stiffness matrix
    };

public:

    AdaptivePlasticBeamForceField();
    virtual ~AdaptivePlasticBeamForceField();

    /////////////////////////////////////
    /// This is inhereted from BaseObject
    /////////////////////////////////////
    virtual void init() override;
    virtual void reinit() override;
    virtual void draw(const VisualParams* vparams) override;


    /////////////////////////////////////
    /// ForceField Interface
    /////////////////////////////////////
    virtual void addForce(const MechanicalParams* mparams, DataVecDeriv& dataf,
                          const DataVecCoord& datax, const DataVecDeriv& v);

    virtual void addDForce(const MechanicalParams* mparams, DataVecDeriv& datadF,
                           const DataVecDeriv& datadX);

    void addKToMatrix(const MechanicalParams* mparams, const MultiMatrixAccessor* matrix);

    void computeStiffness(int beam, BeamLocalMatrices& beamLocalMatrices);

protected:

    SingleLink<AdaptiveBeamForceFieldAndMass<DataTypes>, BPInterpolation, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_interpolation;

    vector<BeamPlasticVariables> m_plasticVariables;

    void applyStiffnessLarge(VecDeriv& df, const VecDeriv& dx, int beam, Index nd0Id, Index nd1Id, const double& factor);

};


/// Instantiate the templates so that they are not instiated in each translation unit (see )
#if !defined(SOFA_PLUGIN_BEAMADAPTER_ADAPTVEPLASTICBEAMFORCEFIELD_CPP)
    extern template class SOFA_BEAMADAPTER_API AdaptivePlasticBeamForceField<Rigid3Types>;
#endif

} // namespace _adaptiveplasticbeamforcefield_

} // namespace sofa::plugin::beamadapter::component::forcefield