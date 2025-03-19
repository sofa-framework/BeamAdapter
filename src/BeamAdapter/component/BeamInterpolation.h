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
#include <BeamAdapter/utils/BeamSection.h>
#include <BeamAdapter/component/BaseBeamInterpolation.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/component/statecontainer/MechanicalObject.h>

namespace beamadapter
{

using sofa::helper::OptionsGroup;
using sofa::core::topology::BaseMeshTopology;
using sofa::core::behavior::MechanicalState;
using sofa::component::statecontainer::MechanicalObject;

/*!
 * \class BeamInterpolation
 *
 * This class implements a Sofa Component that provide interpolation method to compute Finite Element elastic force and mass based on
 * Adaptive 6D beam elements.
 * - Adaptive beam interpolation
 * - Adaptive Force and Mass computation
 * - Adaptive Mapping
 *
 * AdaptiveBeam Interpolation provides the basis of the Beam computation
 * As the computation is adaptive, the interpolation can be modified at each time step.
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 *
 */
template<class DataTypes>
class BeamInterpolation : public BaseBeamInterpolation<DataTypes>
{
public:
    SOFA_CLASS( SOFA_TEMPLATE(BeamInterpolation, DataTypes) , SOFA_TEMPLATE(BaseBeamInterpolation, DataTypes));

    using Coord = typename DataTypes::Coord;
    using VecCoord = typename DataTypes::VecCoord;
    using Real = typename Coord::value_type;

    using Deriv = typename DataTypes::Deriv;
    using VecDeriv = typename DataTypes::VecDeriv;

    using Vec2 = sofa::type::Vec<2, Real>;
    using Vec3 = sofa::type::Vec<3, Real>;
    using Vec3NoInit = sofa::type::VecNoInit<3, Real>;
    using Quat = sofa::type::Quat<Real>;
    using VectorVec3 = type::vector <Vec3>;

    using Transform = typename sofa::defaulttype::SolidTypes<Real>::Transform;
    using SpatialVector = typename sofa::defaulttype::SolidTypes<Real>::SpatialVector;

    using PointID = BaseMeshTopology::PointID;
    using ElementID = BaseMeshTopology::EdgeID;
    using VecElementID = type::vector<BaseMeshTopology::EdgeID>;
    using VecEdges = type::vector<BaseMeshTopology::Edge>;    

    using BaseBeamInterpolation<DataTypes>::d_componentState;
public:
    BeamInterpolation() ;
    virtual ~BeamInterpolation() override = default;

    //////////////////////////////////// Exposing this object in the factory ///////////////////////
    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T* obj, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
    {
        if (dynamic_cast<MechanicalState<DataTypes>*>(context->getMechanicalState()) == nullptr)
        {
            arg->logError(std::string("No mechanical state with the datatype '") + DataTypes::Name() +
                "' found in the context node.");
            return false;
        }
        return sofa::core::objectmodel::BaseObject::canCreate(obj, context, arg);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////// Inherited from Base ///////////////////////////////////////
    void init() override ;
    void bwdInit() override ;
    void reinit() override ;
    void reset() override ;

    //TODO(dmarchal@cduriez) Ca me semble détourner l'API pour faire des choses par surprise. A mon avis la bonne solution
    //est d'implémenter un vrai binding Python pour BeamInterpolation. Avec une fonction updateInterpolation
    /// In the context of beam interpolation, this function (easily access with Python) is used to update the interpolation (input / output)
    void storeResetState() override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////


    void updateInterpolation();
    /**
     * @brief Returns true if the interpolation is specified in the scene file (case of saved executed scenes...)
     */
    bool interpolationIsAlreadyInitialized();
    bool verifyTopology();
    void computeCrossSectionInertiaMatrix();

    const BeamSection& getBeamSection(sofa::Index beamId) override { 
        SOFA_UNUSED(beamId);
        return this->m_constantSection; 
    }
    void getInterpolationParameters(sofa::Index beamId, Real &_L, Real &_A, Real &_Iy , Real &_Iz,
                               Real &_Asy, Real &_Asz, Real &J) override;
    void getMechanicalParameters(sofa::Index beamId, Real& youngModulus, Real& cPoisson, Real& massDensity) override;

    void getTangentUsingSplinePoints(unsigned int edgeInList, const Real& baryCoord, const sofa::core::ConstVecCoordId &vecXId, Vec3& t );

  
    /// computeActualLength => given the 4 control points of the spline, it provides an estimate of the length (using gauss points integration)
   

    

    virtual void getCurvAbsAtBeam(const unsigned int& edgeInList_input, const Real& baryCoord_input, Real& x_output) {}
    virtual void getBeamAtCurvAbs(const Real& x_input, unsigned int& edgeInList_output, Real& baryCoord_output, unsigned int start = 0) {}


    Data<helper::OptionsGroup>   crossSectionShape;

    /// Circular Cross Section
    Data<Real>          d_radius;
    Data<Real>          d_innerRadius;

    /// Square Cross Section
    Data<Real>          d_sideLength;

    /// Elliptic Cross Section
    Data<Real>          d_smallRadius;
    Data<Real>          d_largeRadius;

    /// Rectangular Cross Section
    Data<Real>          d_lengthY;
    Data<Real>          d_lengthZ;
    Data<bool>          d_dofsAndBeamsAligned;

    Real          m_defaultYoungModulus;
    Real          m_defaultPoissonRatio;
    Data<type::vector<Real>>          d_defaultYoungModulus;
    Data<type::vector<Real>>          d_poissonRatio;

    Data<bool>          d_straight;

    virtual void clear();
    virtual void addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1, const Real &angle);
    virtual void getSamplingParameters(type::vector<Real>& xP_noticeable,
                                       type::vector<int>& nbP_density) ;
    Real getRestTotalLength() override;
    void getCollisionSampling(Real &dx, const Real& x_localcurv_abs) override;
    void getNumberOfCollisionSegment(Real &dx, unsigned int &numLines) override;

    void setTransformBetweenDofAndNode(int beam, const Transform &DOF_H_Node, unsigned int zeroORone );
    void getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest) override;

    /////////////////////////// Deprecated Methods  ////////////////////////////////////////// 
    [[deprecated("Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06")]]
    unsigned int getNumBeamsNotUnderControl() {
        msg_warning() << "Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06";
        return 0;
    }

protected :
    /// INPUT / OUTPUT FOR DOING EXTERNAL COMPUTATION OF Beam Interpolation (use it as a kind of data engine)
    ///Input 1. VecID => (1) "current" Pos, Vel    (2) "free" PosFree, VelFree   (3) "rest" PosRest, V=0
    Data< OptionsGroup > d_vecID;
    ///Input 2. Vector of 2-tuples (indice of the beam   ,   barycentric between 0 and 1)
    Data< type::vector< Vec2 > > d_InterpolationInputs;

    ///Output
    Data< VecCoord > d_InterpolatedPos;
    Data< VecDeriv > d_InterpolatedVel;

    /// GEOMETRICAL COMPUTATION (for now we suppose that the radius of the beam do not vary in space / in time)
    BeamSection      m_constantSection;
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_BEAMINTERPOLATION_CPP)
extern template class SOFA_BEAMADAPTER_API BeamInterpolation<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace beamadapter
