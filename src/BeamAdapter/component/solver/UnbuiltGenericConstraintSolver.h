/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_UNBUILTGENERICCONSTRAINTSOLVER_H
#define SOFA_COMPONENT_CONSTRAINTSET_UNBUILTGENERICCONSTRAINTSOLVER_H

#include <SofaConstraint/GenericConstraintSolver.h>
#include <SofaConstraint/ConstraintSolverImpl.h>
#include <sofa/core/behavior/BaseConstraint.h>
#include <sofa/core/behavior/ConstraintSolver.h>
#include <sofa/core/behavior/BaseConstraintCorrection.h>

#include <sofa/simulation/Node.h>
#include <sofa/simulation/MechanicalVisitor.h>

#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaBaseLinearSolver/SparseMatrix.h>

#include <sofa/helper/set.h>
#include <sofa/helper/map.h>

#include <BeamAdapter/initBeamAdapter.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

using namespace sofa::defaulttype;
using namespace sofa::component::linearsolver;
using namespace helper::system::thread;
using core::behavior::ConstraintResolution;

class UnbuiltGenericConstraintSolver;

class SOFA_BEAMADAPTER_API UnbuiltGenericConstraintProblem : public GenericConstraintProblem
{
public:
    // For unbuilt version :
    SparseMatrix<double> Wdiag;
    std::list<unsigned int> constraints_sequence;
    std::vector<core::behavior::BaseConstraintCorrection*> cclist_elem1, cclist_elem2;

    UnbuiltGenericConstraintProblem(const GenericConstraintProblem& cp)
        : GenericConstraintProblem(cp)
    {
        this->scaleTolerance=true;
        this->allVerified=true;
        this->sor=1.0;
    }
    ~UnbuiltGenericConstraintProblem() { this->freeConstraintResolutions(); }

    void unbuiltGaussSeidel(double timeout=0, UnbuiltGenericConstraintSolver* solver = NULL);
};



class SOFA_BEAMADAPTER_API UnbuiltGenericConstraintSolver : public GenericConstraintSolver
{
    typedef std::vector<core::behavior::BaseConstraintCorrection*> list_cc;
    typedef std::vector<list_cc> VecListcc;
    typedef sofa::core::MultiVecId MultiVecId;

public:
    SOFA_CLASS(UnbuiltGenericConstraintSolver, GenericConstraintSolver);

    bool buildSystem(const core::ConstraintParams * /*cParams*/, MultiVecId res1, MultiVecId res2=MultiVecId::null());
    bool solveSystem(const core::ConstraintParams * /*cParams*/, MultiVecId res1, MultiVecId res2=MultiVecId::null());

protected:
    UnbuiltGenericConstraintSolver();
    virtual ~UnbuiltGenericConstraintSolver();

    UnbuiltGenericConstraintProblem *unbuit_current_cp;
};

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
