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

#include <SofaSimpleFem/../../SofaSimpleFem_test/ForceFieldTestCreation.h>
#include "../../../component/forcefield/AdaptivePlasticBeamForceField.h"

namespace sofa
{

/**  Test suite for AdaptivePlasticBeamForceField: we check if the accurate forces are computed
*/
template <typename _AdaptivePlasticBeamForceField>
struct AdaptivePlasticBeamForceField_test : public ForceField_test<_AdaptivePlasticBeamForceField>
{

    typedef _AdaptivePlasticBeamForceField ForceType;
    typedef ForceField_test<_AdaptivePlasticBeamForceField> Inherited;
    typedef typename ForceType::DataTypes DataTypes;

    typedef typename ForceType::VecCoord VecCoord;
    typedef typename ForceType::VecDeriv VecDeriv;
    typedef typename ForceType::Coord Coord;
    typedef typename ForceType::Deriv Deriv;
    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef typename Coord::value_type Real;
    typedef sofa::defaulttype::Vec<3, Real> Vec3;
    typedef sofa::defaulttype::Vec<6, Real> Vec6;

    typedef component::container::MechanicalObject<DataTypes> DOF;

    VecCoord x;
    VecDeriv v, f;

    AdaptivePlasticBeamForceField_test() :Inherited::ForceField_test(std::string(BEAMADAPTER_TEST_SCENES_DIR) + "/" + "AdaptivePlasticBeamForceField.scn")
    {
        //Position
        x.resize(2);
        DataTypes::setCPos(x[0], { 0, 0, 1 } );
        DataTypes::setCRot(x[0], { 0, 0, 0, 1 } );
        DataTypes::setCPos(x[1], { 1, 0, 1 } );
        DataTypes::setCRot(x[1], { 0, 0, 0, 1 } );
        //Velocity
        v.resize(2);
        DataTypes::set(v[0], 0, 0, 0, 0, 0, 0);
        DataTypes::set(v[1], 0, 0, 0, 0, 0, 0);
        //Expected force
        f.resize(2);
        Vec6 fbegin(-0.25, 0, 0, 0, 0, 0);
        Vec6 fend(0.25, 0, 0, 0, 0, 0);
        DataTypes::set(f[0], fbegin[0], fbegin[1], fbegin[2], fbegin[3], fbegin[4], fbegin[5]);
        DataTypes::set(f[1], fend[0], fend[1], fend[2], fend[3], fend[4], fend[5]);

        // Init simulation
        sofa::simulation::getSimulation()->init(Inherited::node.get());
    }

    //Test the value of the force it should be equal for each vertex to Pressure*area/4
    void test_valueForce()
    {
        // run the forcefield_test
        Inherited::run_test(x, v, f);
    }
};

// ========= Define the list of types to instanciate.
//using ::testing::Types;
typedef ::testing::Types<
    sofa::plugin::beamadapter::component::forcefield::_adaptiveplasticbeamforcefield_::AdaptivePlasticBeamForceField<defaulttype::Rigid3Types>
> TestTypes; // the types to instanciate.



// ========= Tests to run for each instanciated type
TYPED_TEST_SUITE(AdaptivePlasticBeamForceField_test, TestTypes);

// test case
TYPED_TEST(AdaptivePlasticBeamForceField_test, extension)
{
    this->errorMax *= 100;
    this->deltaRange = std::make_pair(1, this->errorMax * 10);
    this->debug = false;

    // run test
    this->test_valueForce();
}

} // namespace sofa
