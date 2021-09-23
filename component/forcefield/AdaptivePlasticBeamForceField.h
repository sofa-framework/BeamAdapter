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

using sofa::component::fem::_beaminterpolation_::BeamInterpolation;
using sofa::component::fem::_beaminterpolation_::BeamGeometry;
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
    typedef sofa::type::RGBAColor RGBAColor;
    typedef AdaptiveBeamForceFieldAndMass<DataTypes>::BeamLocalMatrices BeamLocalMatrices;
    typedef BeamInterpolation<DataTypes> BInterpolation;

    typedef typename SolidTypes<Real>::Transform Transform;
    typedef typename SolidTypes<Real>::SpatialVector SpatialVector;

    typedef Vec<3, Real> Vec3;
    typedef Vec<6, Real> Vec6;
    typedef Vec<9, Real> Vec9;          ///< Second-order tensor in vector notation
    typedef Vec<12, Real> Vec12;
    typedef Mat<9, 9, Real> Matrix9x9;  ///< Fourth-order tensor in vector notation
    typedef Mat<3, 12, Real> Matrix3x12;
    typedef Mat<9, 12, Real> Matrix9x12;
    typedef Mat<12, 12, Real> Matrix12x12; ///< Tangent stiffness matrix

    /** \enum class MechanicalState
     *  \brief Types of mechanical state associated with the (Gauss) integration
     *  points. The POSTPLASTIC state corresponds to points which underwent plastic
     *  deformation, but on which constraints were released so that the plasticity
     *  process stopped.
     */
    enum class MechanicalState {
        ELASTIC = 0,
        PLASTIC = 1,
        POSTPLASTIC = 2,
    };

    ///<3-dimensional Gauss point for reduced integration
    class GaussPoint3
    {
    public:
        GaussPoint3() {}
        GaussPoint3(Real x, Real y, Real z, Real w1, Real w2, Real w3);
        ~GaussPoint3() {}

        auto getNx() const -> const Matrix3x12&;
        void setNx(Matrix3x12 Nx);

        auto getGradN() const -> const Matrix9x12&;
        void setGradN(Matrix9x12 gradN);

        auto getMechanicalState() const -> const MechanicalState;
        void setMechanicalState(MechanicalState newState);

        auto getPrevStress() const -> const Vec9&;
        void setPrevStress(Vec9 newStress);

        auto getWeights() const -> const Vec3&;
        void setWeights(Vec3 weights);

        auto getCoord() const -> const Vec3&;
        void setCoord(Vec3 coord);

        auto getBackStress() const -> const Vec9&;
        void setBackStress(Vec9 backStress);

        auto getYieldStress() const ->const Real;
        void setYieldStress(Real yieldStress);

        auto getPlasticStrain() const -> const Vec9&;
        void setPlasticStrain(Vec9 plasticStrain);

        auto getEffectivePlasticStrain() const ->const Real;
        void setEffectivePlasticStrain(Real effectivePlasticStrain);

    protected:
        Vec3 m_coordinates;
        Vec3 m_weights;

        Matrix3x12 m_Nx; /// Shape functions value for the Gauss point (matrix form)
        Matrix9x12 m_gradN; /// Small strain hypothesis deformation gradient, applied to the beam shape functions (matrix form)
        MechanicalState m_mechanicalState; /// State of the Gauss point deformation (elastic, plastic, or postplastic)
        Vec9 m_prevStress; /// Value of the stress tensor at previous time step
        Vec9 m_backStress; /// Centre of the yield surface, in stress space
        Real m_yieldStress; /// Elastic limit, varying if plastic deformation occurs

        //Plasticity history variables
        Vec9 m_plasticStrain;
        Real m_effectivePlasticStrain;
    };

    ///<3 Real intervals [a1,b1], [a2,b2] and [a3,b3], for 3D reduced integration
    class Interval3
    {
    public:
        //By default, integration is considered over [-1,1]*[-1,1]*[-1,1].
        Interval3();
        Interval3(Real a1, Real b1, Real a2, Real b2, Real a3, Real b3);
        ~Interval3() {}

        auto geta1() const -> Real;
        auto getb1() const -> Real;
        auto geta2() const -> Real;
        auto getb2() const -> Real;
        auto geta3() const -> Real;
        auto getb3() const -> Real;

    protected:
        Real m_a1, m_b1, m_a2, m_b2, m_a3, m_b3;
    };

public:

    typedef Vec<27, GaussPoint3> beamGaussPoints;

    AdaptivePlasticBeamForceField();
    virtual ~AdaptivePlasticBeamForceField();

    /////////////////////////////////////
    /// This is inhereted from BaseObject
    /////////////////////////////////////
    virtual void init() override;
    virtual void bwdInit() override;
    virtual void reinit() override;

    /////////////////////////////////////
    /// ForceField Interface
    /////////////////////////////////////
    virtual void addForce(const MechanicalParams* mparams, DataVecDeriv& dataf,
                          const DataVecCoord& datax, const DataVecDeriv& v);


protected:

    /// Position at the last time step, to handle increments for the plasticity resolution
    VecCoord m_lastPos;

    vector<MechanicalState> m_beamMechanicalStates;

    Data<Real> d_youngModulus;
    Data<Real> d_poissonRatio;
    Data<Real> d_initialYieldStress;
    Matrix9x9 m_genHookesLaw;

    // TO DO: Plastic modulus should not be a constant. Is this approximation relevant ?
    // Otherwise a generic law such as Ramberg-Osgood should be use. (Cf Plugin BeamPlastic)
    Data<Real> d_plasticModulus; ///Linearisation coefficient for the plastic behaviour law

    // Threshold used to compare stress tensor norms to 0. See detailed explanation
    // at the computation of the threshold in the init() method.
    Real m_stressComparisonThreshold;

    /////////////////////////////////////
    /// Gaussian integration
    /////////////////////////////////////

    vector<beamGaussPoints> m_gaussPoints;
    vector<Interval3> m_integrationIntervals;

    void initialiseGaussPoints(int beam, vector<beamGaussPoints> &gaussPoints, const Interval3& integrationInterval);
    void initialiseInterval(int beam, vector<Interval3> &integrationIntervals, const BeamGeometry& beamGeometryParams);

    auto computeNx(Real x, Real y, Real z, Real L, Real A, Real Iy, Real Iz,
                   Real E, Real nu, Real kappaY = 1.0, Real kappaZ = 1.0)-> Matrix3x12;

    auto computeGradN(Real x, Real y, Real z, Real L, Real A, Real Iy, Real Iz,
                      Real E, Real nu, Real kappaY = 1.0, Real kappaZ = 1.0) -> Matrix9x12;

    inline double changeCoordinate(double x, double a, double b)
    {
        return 0.5 * ((b - a) * x + a + b);
    }
    inline double changeWeight(double w, double a, double b)
    {
        return 0.5 * (b - a) * w;
    }

    template <typename LambdaType>
    void integrateBeam(beamGaussPoints& gaussPoints, LambdaType integrationFun);

    /////////////////////////////////////
    /// Methods for plasticity
    /////////////////////////////////////

    /// Computes local displacement of a beam element using the corotational model
    void computeLocalDisplacement(const VecCoord& x, Vec6& U0Local, Vec6& U1Local,
                                  unsigned int beamIndex, unsigned int node0Idx, unsigned int node1Idx,
                                  bool updateBeamMatrices=false);

    /// Computes stress increment on a single point of space, from previous stress and current strain
    void computeStressIncrement(GaussPoint3& gp, Vec9& newStressPoint,
                                const Vec9& strainIncrement, MechanicalState& pointMechanicalState);

    void updateTangentStiffness(unsigned int beamIndex, BeamLocalMatrices& beamLocalMatrices);

    //----- Implementation of the Von Mises yield function -----//

    /// Computes the deviatoric stress from a tensor
    Vec9 deviatoricStress(const Vec9& stressTensor);
    /// Computes the equivalent stress from a tensor
    Real equivalentStress(const Vec9& stressTensor);
    /// Evaluates the Von Mises yield function for given stress tensor and yield stress
    Real vonMisesYield(const Vec9& stressTensor, const Vec9& backStress, const Real yieldStress);
    /// Computes the Von Mises yield function gradient, for a given stress tensor
    Vec9 vonMisesGradient(const Vec9& stressTensor);

    //----- Visualisation -----//

    void draw(const core::visual::VisualParams* vparams) override;
    void drawElement(unsigned int beamIndex, std::vector<defaulttype::Vector3>& visualGaussPoints,
                     std::vector< defaulttype::Vector3 >& centrelinePointsSF,
                     std::vector< defaulttype::Vector3 >& centrelinePointsSpline,
                     std::vector<RGBAColor>& colours, const VecCoord& x);
};


/// Instantiate the templates so that they are not instiated in each translation unit (see )
#if !defined(SOFA_PLUGIN_BEAMADAPTER_ADAPTVEPLASTICBEAMFORCEFIELD_CPP)
    extern template class SOFA_BEAMADAPTER_API AdaptivePlasticBeamForceField<Rigid3Types>;
#endif

} // namespace _adaptiveplasticbeamforcefield_

} // namespace sofa::plugin::beamadapter::component::forcefield