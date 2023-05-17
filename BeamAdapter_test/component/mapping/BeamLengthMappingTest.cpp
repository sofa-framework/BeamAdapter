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
#include <string>
using std::string ;
#include <sofa/component/mapping/testing/MappingTestCreation.h>
#include <sofa/simulation/graph/DAGSimulation.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/simulation/UpdateLinksVisitor.h>
#include <sofa/simulation/InitVisitor.h>

#include <sofa/component/statecontainer/MechanicalObject.h>

#include <BeamAdapter/component/mapping/BeamLengthMapping.h>
#include <BeamAdapter/component/BeamInterpolation.h>

#include <sofa/simulation/common/SceneLoaderXML.h>
using sofa::simulation::SceneLoaderXML ;
using sofa::simulation::Node ;
using sofa::component::statecontainer::MechanicalObject ;

namespace sofa {
  namespace { // anonymous namespace
using namespace core;
using namespace component;
using type::Vec;
using type::Mat;


/**  Test suite for RigidMapping.
The test cases are defined in the #Test_Cases member group.
  */
template <typename _BeamLengthMapping>
struct BeamLengthMappingTest : public sofa::mapping_test::Mapping_test<_BeamLengthMapping>
{

    typedef _BeamLengthMapping BeamLengthMapping;
    typedef sofa::mapping_test::Mapping_test<BeamLengthMapping> Inherit;

    typedef typename BeamLengthMapping::In InDataTypes;
    typedef typename InDataTypes::VecCoord InVecCoord;
    typedef typename InDataTypes::VecDeriv InVecDeriv;
    typedef typename InDataTypes::Coord InCoord;
    typedef typename InDataTypes::Deriv InDeriv;
    typedef sofa::component::statecontainer::MechanicalObject<InDataTypes> InMechanicalObject;
    typedef typename InMechanicalObject::ReadVecCoord  ReadInVecCoord;
    typedef typename InMechanicalObject::WriteVecCoord WriteInVecCoord;
    typedef typename InMechanicalObject::WriteVecDeriv WriteInVecDeriv;
    typedef typename InCoord::Pos Translation;
    typedef typename InCoord::Rot Rotation;
    typedef typename InDataTypes::Real InReal;
    typedef Mat<InDataTypes::spatial_dimensions,InDataTypes::spatial_dimensions,InReal> RotationMatrix;


    typedef typename BeamLengthMapping::Out OutDataTypes;
    typedef typename OutDataTypes::VecCoord OutVecCoord;
    typedef typename OutDataTypes::VecDeriv OutVecDeriv;
    typedef typename OutDataTypes::Coord OutCoord;
    typedef typename OutDataTypes::Deriv OutDeriv;
    typedef sofa::component::statecontainer::MechanicalObject<OutDataTypes> OutMechanicalObject;
    typedef typename OutMechanicalObject::WriteVecCoord WriteOutVecCoord;
    typedef typename OutMechanicalObject::WriteVecDeriv WriteOutVecDeriv;
    typedef typename OutMechanicalObject::ReadVecCoord ReadOutVecCoord;
    typedef typename OutMechanicalObject::ReadVecDeriv ReadOutVecDeriv;


    BeamLengthMapping* beamLengthMapping;

    BeamLengthMappingTest()
    {
        this->errorFactorDJ = 200;

        beamLengthMapping = static_cast<BeamLengthMapping*>( this->mapping );

        // beamLengthMapping::getJs is not yet implemented
        this->flags &= ~Inherit::TEST_getJs;

        // beamLengthMapping::getK is not yet implemented
        //this->flags &= ~Inherit::TEST_getK;

        // beamLengthMapping::applyDJT is not yet implemented
        //this->flags &= ~Inherit::TEST_applyDJT;

    }

    /** @name Test_Cases
      verify the computation of the beam lengths and the derivatives
      */
    ///@{
    /** Two frames are placed + line topology + and the mapping is constructed
    */

    bool testCase1()
    {
        const int Nout=2; // WARNING this number has to be changed to test with more than one beam !!
        const int Nin=3;

        string scene =
                "<?xml version='1.0'?>"
                ""
                "<Node 	name='Root' gravity='0 0 0' time='0' animate='0'>"
                "   			<EulerImplicit rayleighStiffness='0.08' rayleighMass='0.08' printLog='false' />"
                "               <CGLinearSolver iterations='100' threshold='1e-10' tolerance='1e-15' />"
                "               <Mesh name='meshSuture' edges='0 1 1 2' />"
                "               <MechanicalObject template='Rigid' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1   2 0 0 0 0 0 1' showObject='1'/>"
                "                <BeamInterpolation name='Interpol' radius='0.3' defaultYoungModulus='1e7'  DOF0TransformNode0='0.2 0.3 0.1 0 0 0 1   0.12 0.3 0.1 0 0 0.3826834 0.9238795 '  DOF1TransformNode1='-0.3 0.3 0.1 0 0 0 1  -0.2 0.25 0 0 0 0 1' dofsAndBeamsAligned='0' straight='0'/>"
                "               <Node name='Map' > "
                "                      <MechanicalObject template='Vec1d' name='mappedDOFs' position='0.5  0.8'  />"
                "                      <BeamLengthMapping name='beamLMap' geometricStiffness='1' interpolation='@Interpol' input='@DOFs' output='@mappedDOFs' />"
                "               </Node>"
                "</Node> " ;
        this->root = SceneLoaderXML::loadFromMemory ( "testCase1", scene.c_str());



        //std::cout<<"*******************  Get the mapping ";
        BeamLengthMapping* lengthMapping;
        this->root->getTreeObject(lengthMapping);
        this->mapping = lengthMapping;
        this->deltaRange.first= 1;
        this->deltaRange.second= 1000;

        MechanicalObject<defaulttype::Rigid3dTypes>* FromModel = nullptr;
        this->root->getTreeObject(FromModel);
        this->inDofs = FromModel;

        MechanicalObject<defaulttype::Vec1Types>* ToModel = nullptr;
        this->root->getTreeObject(ToModel);
        this->outDofs= ToModel;

        const Data<InVecCoord>& dataInX = *FromModel->read(VecCoordId::position());
        const InVecCoord& xin = dataInX.getValue();

        const Data<OutVecCoord>& dataOutX = *ToModel->read(VecCoordId::position());
        const OutVecCoord& xout = dataOutX.getValue();

        std::cout<<" x in  = "<<xin<<std::endl;

        std::cout<<" x out  = "<<xout<<std::endl;



        // new values of DOFs: exact same position
        InVecCoord xin_new(Nin);
        xin_new[1][0]=1.0;
        xin_new[2][0]=2.0;


        // expected mapped values
        OutVecCoord expectedChildCoords(Nout);
        expectedChildCoords[0][0] = 0.5;
        expectedChildCoords[1][0] = 0.7416993975182119;



        return this->runTest(xin,xout,xin_new,expectedChildCoords);


    }
    bool testCase2()
    {
        const int Nout=1; // WARNING this number has to be changed to test with more than one beam !!
        const int Nin=2;

        string scene =
                "<?xml version='1.0'?>"
                "<Node 	name='Root' gravity='0 0 0' time='0' animate='0'>"
                "    <RequiredPlugin name=\"Sofa.Component.ODESolver.Backward\"/>"
                "    <RequiredPlugin name=\"Sofa.Component.LinearSolver.Iterative\"/>"
                "    <EulerImplicitSolver rayleighStiffness='0.08' rayleighMass='0.08' printLog='false' />"
                "    <CGLinearSolver iterations='100' threshold='1e-10' tolerance='1e-15' />"
                "    <Mesh name='meshSuture' edges='0 1' />"
                "    <MechanicalObject template='Rigid3' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
                "    <BeamInterpolation name='Interpol' radius='0.1'/>"
                "    <Node name='Map' > "
                "           <MechanicalObject template='Vec1' name='mappedDOFs' position='1.0'  />"
                "           <BeamLengthMapping name='beamLMap' geometricStiffness='1' interpolation='@Interpol' input='@DOFs' output='@mappedDOFs' />"
                "    </Node>"
                "</Node> " ;
        this->root = SceneLoaderXML::loadFromMemory ( "testCase1", scene.c_str());



        //std::cout<<"*******************  Get the mapping ";
        BeamLengthMapping* lengthMapping;
        this->root->getTreeObject(lengthMapping);
        this->mapping = lengthMapping;
        this->deltaRange.first= 100;
        this->deltaRange.second= 100000;


        MechanicalObject<defaulttype::Rigid3dTypes>* FromModel = nullptr;
        this->root->getTreeObject(FromModel);
        this->inDofs = FromModel;

        MechanicalObject<defaulttype::Vec1Types>* ToModel = nullptr;
        this->root->getTreeObject(ToModel);
        this->outDofs= ToModel;

        const Data<InVecCoord>& dataInX = *FromModel->read(VecCoordId::position());
        const InVecCoord& xin = dataInX.getValue();

        const Data<OutVecCoord>& dataOutX = *ToModel->read(VecCoordId::position());
        const OutVecCoord& xout = dataOutX.getValue();

        std::cout<<" x in  = "<<xin<<std::endl;

        std::cout<<" x out  = "<<xout<<std::endl;



        // new values in (curve) pos=  0.912591 -0.419907 0 0 0 -0.283387 0.959006
        InVecCoord xin_new(Nin);
        xin_new[1][0]=0.912591;
        xin_new[1][1]=-0.419907;
        xin_new[1][5]= -0.283387;
        xin_new[1][6]= 0.959006;

        // expected mapped values
        OutVecCoord expectedChildCoords(Nout);
        expectedChildCoords[0][0] = 1.0197604809762477;



        return this->runTest(xin,xout,xin_new,expectedChildCoords);


    }


};

// Define the list of types to instanciate. We do not necessarily need to test all combinations.
using ::testing::Types;
typedef Types<
mapping::_beamlengthmapping_::BeamLengthMapping<defaulttype::Rigid3dTypes,defaulttype::Vec1dTypes>
//,mapping::_beamlengthmapping_::BeamLengthMapping<defaulttype::Rigid3fTypes,defaulttype::Vec1fTypes>
> DataTypes; // the types to instanciate.

// Test suite for all the instanciations
TYPED_TEST_CASE(BeamLengthMappingTest, DataTypes);
// first test case, failing
// Error is: Position of mapped particle 1 is wrong: 0.72042110251343161, expected: 0.7416993975182119. difference should be less than 2.2204460492503131e-15 (0.02127829500478029) [on MSVC]
// This failing test could be either 
// - reliable (something wrong has been introduced in SOFA or BeamAdapter) 
// - or the test itself was relying on something wrong (expected test results or the component)
TYPED_TEST( BeamLengthMappingTest , DISABLED_testCase1 )
{
    ASSERT_TRUE(this->testCase1());
}

// second test case
TYPED_TEST( BeamLengthMappingTest , DISABLED_testCase2 )
{
    ASSERT_TRUE(this->testCase2());
}

  }// anonymous namespace
  }//sofa
