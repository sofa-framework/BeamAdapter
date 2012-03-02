/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
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
#ifndef SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEAMFORCEFIELDANDMASS_H
#define SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEAMFORCEFIELDANDMASS_H


#include "initBeamAdapter.h"
#include "BeamInterpolation.h"
#include "WireRestShape.h"
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>

#include <sofa/core/objectmodel/BaseObject.h>


namespace sofa
{

namespace component
{

namespace forcefield
{
    using sofa::helper::vector;
    using namespace sofa::core::topology;
    using namespace sofa::component::fem;
	using engine::WireRestShape;

template<class DataTypes>
class SOFA_BEAMADAPTER_API AdaptiveBeamForceFieldAndMass : public core::behavior::Mass<DataTypes>
{
public:

    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBeamForceFieldAndMass,DataTypes),  SOFA_TEMPLATE(core::behavior::Mass,DataTypes));

   // SOFA_CLASS2(SOFA_TEMPLATE(AdaptiveBeamForceFieldAndMass,DataTypes),  SOFA_TEMPLATE(core::behavior::Mass,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::VecReal VecReal;
    typedef VecCoord Vector;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
	typedef core::objectmodel::Data<VecCoord> DataVecCoord;
	typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;
    typedef  BeamInterpolation<DataTypes> BInterpolation;

    typedef unsigned int Index;
    typedef core::topology::BaseMeshTopology::Edge Element;
    typedef sofa::helper::vector<core::topology::BaseMeshTopology::Edge> VecElement;
	typedef helper::vector<unsigned int> VecIndex;

    typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
    typedef typename sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;

     /// remove ?

    typedef Vec<12, Real> Displacement;        ///< the displacement vector
    typedef Mat<3, 3, Real> Transformation; ///< matrix for rigid transformations like rotations
    typedef Mat<12, 12, Real> StiffnessMatrix;

    /// replace by:

    typedef Vec<3, Real> Vec3;
    typedef Vec<6, Real> Vec6;        ///< the displacement vector
    typedef Mat<6, 6, Real> Matrix6x6;

protected :
    // pointer to the interpolation
    BInterpolation* m_interpolation;
    Data< helper::vector< std::string > > m_interpolationPath;
	SingleLink<AdaptiveBeamForceFieldAndMass<DataTypes>, WireRestShape<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_instrumentParameters;

public:
    AdaptiveBeamForceFieldAndMass( )
    : m_interpolation(NULL)
    , m_interpolationPath(initData(&m_interpolationPath,"interpolation","Path to the Interpolation component on scene"))
    , _timoshenko(initData(&_timoshenko,true,"timoshenko","use Timoshenko beam (non-null section shear area)"))
    , _computeMass(initData(&_computeMass,true,"computeMass","if false, only compute the stiff elastic model"))
    , _massDensity(initData(&_massDensity,(Real)1.0,"massDensity", "Density of the mass (usually in kg/m^3)" ))
    , _shearStressComputation(initData(&_shearStressComputation, true, "shearStressComputation","if false, suppress the shear stress in the computation"))
	, m_instrumentParameters(initLink("instrumentParameters", "link to an object specifying physical parameters based on abscissa"))
    {
        _localBeamMatrices.resize(2);
    }

    virtual void init();
    virtual void reinit();

    // Mass Interface
    virtual  void addMDx(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecDeriv& dx, double factor);

    virtual void addMToMatrix(const core::MechanicalParams *mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix);

    virtual void addMBKToMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix);

    virtual  void accFromF(const core::MechanicalParams* /* PARAMS FIRST */, DataVecDeriv& , const DataVecDeriv& )
    {serr<<"accFromF can not be implemented easily: It necessitates a solver because M^-1 is not available"<<sendl;}

    virtual double getKineticEnergy(const core::MechanicalParams* /* PARAMS FIRST */, const DataVecDeriv& )  const ///< vMv/2 using dof->getV()
    {serr<<"getKineticEnergy not yet implemented"<<sendl;return 0;}

    virtual void addGravityToV(const core::MechanicalParams* /* PARAMS FIRST */, DataVecDeriv& )
    {serr<<"addGravityToV not implemented yet"<<sendl;}



    // ForceField Interface
    virtual void addForce (const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);

    virtual void  addDForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX );

    virtual double getPotentialEnergy(const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord& ) const {serr<<"getPotentialEnergy not yet implemented"<<sendl; return 0; }

    void addKToMatrix(const core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix);

    void draw(const core::visual::VisualParams* vparams);

    void drawElement(const core::visual::VisualParams* vparams, int beam, Transform &global_H0_local, Transform &global_H1_local);


    Data<bool> _timoshenko;
    Data<bool> _computeMass;
    Data<Real> _massDensity;
    Data<bool> _shearStressComputation;


  //  void setPathToInterpolation(const std::string &o){m_interpolationPath.setValue(o);};
protected:


    void applyMassLarge( VecDeriv& df, const VecDeriv& dx, int bIndex, Index nd0Id, Index nd1Id, const double &factor);
    void applyStiffnessLarge( VecDeriv& df, const VecDeriv& dx, int beam, Index nd0Id, Index nd1Id, const double &factor );

//    void computeTransform(int i, Index a, Index b);
	//void computeRotationLarge( Transformation &r, const Vector &p, Index a, Index b);
//	void accumulateForceLarge( VecDeriv& f, const VecCoord& x, int i, Index a, Index b);
	//void accumulateDampingLarge( Vector& f, Index elementIndex );


    class BeamLocalMatrices{

    public:
        BeamLocalMatrices(){}
        BeamLocalMatrices(const BeamLocalMatrices &v)
        {
            this->k_loc00 = v.k_loc00;
            this->k_loc10 = v.k_loc10;
            this->k_loc11 = v.k_loc11;
            this->k_loc01 = v.k_loc01;

            this->m_loc00 = v.m_loc00;
            this->m_loc11 = v.m_loc11;
            this->m_loc01 = v.m_loc01;
            this->m_loc10 = v.m_loc10;

            this->loc0_Ad_ref = v.loc0_Ad_ref;
            this->loc1_Ad_ref = v.loc1_Ad_ref;
        }
        ~BeamLocalMatrices(){}
        // stiffness Matrices
        Matrix6x6 k_loc00, k_loc01, k_loc10, k_loc11;
        // mass Matrices
        Matrix6x6 m_loc00, m_loc01, m_loc10, m_loc11;
        // adjoint Matrices
        Matrix6x6 loc0_Ad_ref, loc1_Ad_ref;

    };

    void computeStiffness(int beam, BeamLocalMatrices& beamMatrices);
    void computeMass(int beam,BeamLocalMatrices& beamMatrices);

    helper::vector<BeamLocalMatrices> _localBeamMatrices;

    DataVecDeriv data_G;
    // compute the gravity vector (if necessary)
    void computeGravityVector();
    Vec3 gravity;





};


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEAMFORCEFIELDANDMASS_H */
