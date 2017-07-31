#include <string>
using std::string ;

#include <SofaTest/Sofa_test.h>
#include <sofa/helper/BackTrace.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>

using sofa::helper::vector;
using namespace sofa::defaulttype;



#include "GeometricalUtils.cpp"

namespace sofa
{




template <class Real>
class GeometricalUtilsTest : public Sofa_test<>
{
public:

    GeometricalUtilsTest(){}
    
    void static testComputeSplineLength( ){
        Vec<3, Real> P0((Real)1.0, (Real)3.2, (Real)2.8 );
        Vec<3, Real> P1((Real)2.0, (Real)3.6, (Real)2.4 );
        Vec<3, Real> P2((Real)3.0, (Real)3.0, (Real)2.9 );
        Vec<3, Real> P3((Real)4.0, (Real)3.2, (Real)2.8 );


        // verification


        Real length_verif=0.0;
        Vec<3, Real> seg, pos;
        pos=P0;

        int i=0;
        for (Real bx=0.00001; bx<1.00000001; bx+=0.00001)
        {
            // compute length
            seg  = -pos;
            pos = P0*(1-bx)*(1-bx)*(1-bx) + P1*3*bx*(1-bx)*(1-bx) + P2*3*bx*bx*(1-bx) + P3*bx*bx*bx;
            seg += pos;
            length_verif += seg.norm();
            i++;
        }


        Real length=0.0;
        computeSplineLength(length,P0,P1,P2,P3);

        EXPECT_LE(fabs(length-length_verif), fabs(0.001*length_verif) );


        std::cout<<"length "<<length<<" length_verif "<<length_verif<< "  i ="<<i <<std::endl;

        return;
        
        
    }

};


typedef  GeometricalUtilsTest<double> GeometricalUtilsTestdouble;



TEST(GeometricalUtilsTest, TestSpline)
{
    GeometricalUtilsTestdouble::testComputeSplineLength();
}

/*
std::vector< std::vector<std::string> > unintvalues = {
    {"", "t"}, {"toto", "f"}, {"./mappedTopo/topoMap", "t" }
};


TEST_P(WireRestShapeTest, SimpleSceneValid) {
    ASSERT_NO_THROW(this->simpleSceneTest(GetParam() ) );
}

INSTANTIATE_TEST_CASE_P(checkRestShapeInits,
                        WireRestShapeTest,
                        ::testing::ValuesIn(unintvalues));

*/
} // namespace SOFA
