#include <gtest/gtest.h>
#include <sofa/testing/BaseSimulationTest.h>

#include "../../../component/forcefield/AdaptivePlasticBeamForceField.h"
#include "../../../component/forcefield/AdaptivePlasticBeamForceField.inl"

using sofa::plugin::beamadapter::component::forcefield::_adaptiveplasticbeamforcefield_::AdaptivePlasticBeamForceField;
using sofa::defaulttype::Rigid3dTypes;

typedef AdaptivePlasticBeamForceField<Rigid3dTypes> AdaptivePlaticBeamForceField_Rigid;
typedef AdaptivePlaticBeamForceField_Rigid::Interval3 Interval3;
typedef AdaptivePlaticBeamForceField_Rigid::GaussPoint3 GaussPoint3;
typedef AdaptivePlaticBeamForceField_Rigid::beamGaussPoints beamGaussPoints;
typedef AdaptivePlaticBeamForceField_Rigid::Vec3 Vec3;

inline double changeCoordinate(double x, double a, double b)
{
    return 0.5 * ((b - a) * x + a + b);
}
inline double changeWeight(double w, double a, double b)
{
    return 0.5 * (b - a) * w;
}

template <typename LambdaType>
void integrateBeam(beamGaussPoints& gaussPoints, LambdaType integrationFun)
{
    //Apply a generic (lambda) integration function to each Gauss point of a beam element
    for (unsigned int gp = 0; gp < gaussPoints.size(); gp++)
    {
        integrationFun(gaussPoints[gp]);
    }
}

TEST(Quadrature, BaseIntervalVolume)
{
	Interval3 interval = Interval3(-1, 1, -1, 1, -1, 1);

    const double sqrt3_5 = sofa::helper::rsqrt(3.0 / 5);
    Vec3 canonical3NodesCoordinates = { -sqrt3_5, 0, sqrt3_5 };
    Vec3 canonical3NodesWeights = { 5.0 / 9, 8.0 / 9, 5.0 / 9 };

	beamGaussPoints gaussPoints;

    unsigned int gaussPointIt = 0;
    for (unsigned int i = 0; i < 3; i++)
    {
        double x = canonical3NodesCoordinates[i];
        double w1 = canonical3NodesWeights[i];
        // Changing first coordinate and weight to adapt to the integration interval
        double a1 = interval.geta1();
        double b1 = interval.getb1();
        double xChanged = changeCoordinate(x, a1, b1);
        double w1Changed = changeWeight(w1, a1, b1);

        // with [-1,1] interval, changeCoordinate and changeWeight should not make any difference
        ASSERT_DOUBLE_EQ(x, xChanged);
        ASSERT_DOUBLE_EQ(w1, w1Changed);

        for (unsigned int j = 0; j < 3; j++)
        {
            double y = canonical3NodesCoordinates[j];
            double w2 = canonical3NodesWeights[j];
            // Changing second coordinate and weight to adapt to the integration interval
            double a2 = interval.geta2();
            double b2 = interval.getb2();
            double yChanged = changeCoordinate(y, a2, b2);
            double w2Changed = changeWeight(w2, a2, b2);

            // with [-1,1] interval, changeCoordinate and changeWeight should not make any difference
            ASSERT_DOUBLE_EQ(y, yChanged);
            ASSERT_DOUBLE_EQ(w2, w2Changed);

            for (unsigned int k = 0; k < 3; k++)
            {
                double z = canonical3NodesCoordinates[k];
                double w3 = canonical3NodesWeights[k];
                // Changing third coordinate and weight to adapt to the integration interval
                double a3 = interval.geta3();
                double b3 = interval.getb3();
                double zChanged = changeCoordinate(z, a3, b3);
                double w3Changed = changeWeight(w3, a3, b3);

                // with [-1,1] interval, changeCoordinate and changeWeight should not make any difference
                ASSERT_DOUBLE_EQ(z, zChanged);
                ASSERT_DOUBLE_EQ(w3, w3Changed);

                GaussPoint3 newGaussPoint = GaussPoint3(xChanged, yChanged, zChanged, w1Changed, w2Changed, w3Changed);
                gaussPoints[gaussPointIt] = newGaussPoint;
                gaussPointIt++;
            }
        }
    } // end for (unsigned int i = 0; i < 3; i++)

    double volume = 0;
    typedef std::function<void(GaussPoint3&)> IntegrationLambda;
    IntegrationLambda integrateVolume = [&](GaussPoint3& gp)
    {
        Vec3 gpWeights = gp.getWeights();
        double prod = gpWeights[0] * gpWeights[1] * gpWeights[2];
        volume += prod;
    };

    integrateBeam(gaussPoints, integrateVolume);
    ASSERT_DOUBLE_EQ(volume, 8.0); // Expected volume of a [-1,1]*[-1,1]*[-1,1] cube
}

TEST(Quadrature, ChangedIntervalVolume)
{
    Interval3 interval = Interval3(-2, 2, -1, 5, 0, 3);

    const double sqrt3_5 = sofa::helper::rsqrt(3.0 / 5);
    Vec3 canonical3NodesCoordinates = { -sqrt3_5, 0, sqrt3_5 };
    Vec3 canonical3NodesWeights = { 5.0 / 9, 8.0 / 9, 5.0 / 9 };

    beamGaussPoints gaussPoints;

    unsigned int gaussPointIt = 0;
    for (unsigned int i = 0; i < 3; i++)
    {
        double x = canonical3NodesCoordinates[i];
        double w1 = canonical3NodesWeights[i];
        // Changing first coordinate and weight to adapt to the integration interval
        double a1 = interval.geta1();
        double b1 = interval.getb1();
        double xChanged = changeCoordinate(x, a1, b1);
        double w1Changed = changeWeight(w1, a1, b1);

        for (unsigned int j = 0; j < 3; j++)
        {
            double y = canonical3NodesCoordinates[j];
            double w2 = canonical3NodesWeights[j];
            // Changing second coordinate and weight to adapt to the integration interval
            double a2 = interval.geta2();
            double b2 = interval.getb2();
            double yChanged = changeCoordinate(y, a2, b2);
            double w2Changed = changeWeight(w2, a2, b2);

            for (unsigned int k = 0; k < 3; k++)
            {
                double z = canonical3NodesCoordinates[k];
                double w3 = canonical3NodesWeights[k];
                // Changing third coordinate and weight to adapt to the integration interval
                double a3 = interval.geta3();
                double b3 = interval.getb3();
                double zChanged = changeCoordinate(z, a3, b3);
                double w3Changed = changeWeight(w3, a3, b3);

                GaussPoint3 newGaussPoint = GaussPoint3(xChanged, yChanged, zChanged, w1Changed, w2Changed, w3Changed);
                gaussPoints[gaussPointIt] = newGaussPoint;
                gaussPointIt++;
            }
        }
    } // end for (unsigned int i = 0; i < 3; i++)

    double volume = 0;
    typedef std::function<void(GaussPoint3&)> IntegrationLambda;
    IntegrationLambda integrateVolume = [&](GaussPoint3& gp)
    {
        Vec3 gpWeights = gp.getWeights();
        double prod = gpWeights[0] * gpWeights[1] * gpWeights[2];
        volume += prod;
    };

    integrateBeam(gaussPoints, integrateVolume);
    ASSERT_DOUBLE_EQ(volume, 72.0); // Expected volume of a [-2, 2]*[-1, 5]*[0, 3]; cube
}