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

#include "UnbuiltGenericConstraintSolver.h"

#include <sofa/core/visual/VisualParams.h>

#include <sofa/simulation/AnimateVisitor.h>
#include <sofa/simulation/BehaviorUpdatePositionVisitor.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/SolveVisitor.h>
#include <sofa/simulation/VectorOperations.h>

#include <sofa/simulation/Simulation.h>
#include <sofa/helper/gl/template.h>
#include <sofa/helper/gl/Axis.h>
#include <sofa/helper/gl/Cylinder.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/system/thread/CTime.h>
#include <math.h>
#include <iostream>

#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

UnbuiltGenericConstraintSolver::UnbuiltGenericConstraintSolver() : GenericConstraintSolver()
{
}


UnbuiltGenericConstraintSolver::~UnbuiltGenericConstraintSolver()
{
}

bool UnbuiltGenericConstraintSolver::buildSystem(const core::ConstraintParams *cParams, MultiVecId /*res1*/, MultiVecId /*res2*/)
{
    unbuit_current_cp = new UnbuiltGenericConstraintProblem( *(this->current_cp) );
    unsigned int numConstraints = 0;

    sofa::helper::AdvancedTimer::stepBegin("Accumulate Constraint");
    // mechanical action executed from root node to propagate the constraints
    simulation::MechanicalResetConstraintVisitor(cParams).execute(context);
    // calling buildConstraintMatrix
    //simulation::MechanicalAccumulateConstraint(&cparams /* PARAMS FIRST */, core::MatrixDerivId::constraintMatrix(), numConstraints).execute(context);

    MechanicalSetConstraint(cParams, core::MatrixDerivId::holonomicC(), numConstraints).execute(context);
    MechanicalAccumulateConstraint2(cParams, core::MatrixDerivId::holonomicC()).execute(context);

    // suppress the constraints that are on DOFS currently concerned by projective constraint
    core::MechanicalParams mparams = core::MechanicalParams(*cParams);
    simulation::MechanicalProjectJacobianMatrixVisitor(&mparams).execute(context);

    sofa::helper::AdvancedTimer::stepEnd  ("Accumulate Constraint");
    sofa::helper::AdvancedTimer::valSet("numConstraints", numConstraints);

    unbuit_current_cp->clear(numConstraints);

    sofa::helper::AdvancedTimer::stepBegin("Get Constraint Value");
    MechanicalGetConstraintViolationVisitor(cParams /* PARAMS FIRST */, &unbuit_current_cp->dFree).execute(context);
    sofa::helper::AdvancedTimer::stepEnd ("Get Constraint Value");

    sofa::helper::AdvancedTimer::stepBegin("Get Constraint Resolutions");
    MechanicalGetConstraintResolutionVisitor(cParams /* PARAMS FIRST */, unbuit_current_cp->constraintsResolutions).execute(context);
    sofa::helper::AdvancedTimer::stepEnd("Get Constraint Resolutions");

    if (this->f_printLog.getValue()) sout<<"UnbuiltGenericConstraintSolver: "<<numConstraints<<" constraints"<<sendl;

    for (unsigned int i=0;i<constraintCorrections.size();i++)
    {
        core::behavior::BaseConstraintCorrection* cc = constraintCorrections[i];
        cc->resetForUnbuiltResolution(unbuit_current_cp->getF(), unbuit_current_cp->constraints_sequence);
    }

    SparseMatrix<double>* Wdiag = &unbuit_current_cp->Wdiag;
    Wdiag->resize(numConstraints, numConstraints);

    // for each contact, the pair of constraint correction that is involved with the contact is memorized
    unbuit_current_cp->cclist_elem1.resize(numConstraints);
    unbuit_current_cp->cclist_elem2.resize(numConstraints);
            unsigned int nbObjects = 0;
    for(unsigned int i=0; i<numConstraints;)
    {
        nbObjects++;
        unsigned int l = unbuit_current_cp->constraintsResolutions[i]->nbLines;

        bool elem1 = false, elem2 = false;
        for (unsigned int j=0;j<constraintCorrections.size();j++)
        {
            core::behavior::BaseConstraintCorrection* cc = constraintCorrections[j];
            if(cc->hasConstraintNumber(i))
            {
                if(elem1)
                {
                    unbuit_current_cp->cclist_elem2[i] = cc;
                    elem2=true;
                }
                else
                {
                    unbuit_current_cp->cclist_elem1[i] = cc;
                    elem1=true;
                }

            }
        }

        if(elem1)	// for each contact, the pair of constraintcorrection is called to add the contribution
            unbuit_current_cp->cclist_elem1[i]->getBlockDiagonalCompliance(Wdiag, i, i+l-1);
        else
            serr<<"WARNING: no constraintCorrection found for constraint"<<i<<sendl;

        if(elem2)
            unbuit_current_cp->cclist_elem2[i]->getBlockDiagonalCompliance(Wdiag, i, i+l-1);
        else
            unbuit_current_cp->cclist_elem2[i] = (NULL);

        double** w =  unbuit_current_cp->getW();
        for(unsigned int m=i; m<i+l; m++)
            for(unsigned int n=i; n<i+l; n++)
                w[m][n] = Wdiag->element(m, n);

        i += l;
    }





    if ( displayTime.getValue() )
    {
        sout<<" build_LCP " << ( (double) timer.getTime() - time)*timeScale<<" ms" <<sendl;
        time = (double) timer.getTime();
    }
    return true;
}



bool UnbuiltGenericConstraintSolver::solveSystem(const core::ConstraintParams * /*cParams*/, MultiVecId /*res1*/, MultiVecId /*res2*/)
{
    unbuit_current_cp->tolerance = tolerance.getValue();
    unbuit_current_cp->maxIterations = maxIt.getValue();
    unbuit_current_cp->scaleTolerance = scaleTolerance.getValue();
    unbuit_current_cp->allVerified = allVerified.getValue();
    unbuit_current_cp->sor = sor.getValue();


    sofa::helper::AdvancedTimer::stepBegin("ConstraintsUnbuiltGaussSeidel");
    unbuit_current_cp->unbuiltGaussSeidel(0, this);


    sofa::helper::AdvancedTimer::stepEnd("ConstraintsUnbuiltGaussSeidel");


    if ( displayTime.getValue() )
    {
        sout<<" TOTAL solve_LCP " <<( (double) timer.getTime() - time)*timeScale<<" ms" <<sendl;
        time = (double) timer.getTime();
    }

    if(this->f_printLog.getValue())
    {
        std::cout<<"i \t \t dfree \t \t f"<<std::endl;
        for (unsigned int i=0; i<unbuit_current_cp->getDimension(); i++)
        {
            std::cout<<""<<i<<"\t \t "<<unbuit_current_cp->getDfree()[i]<<" \t \t"<<unbuit_current_cp->getF()[i]<<std::endl;
        }
        //afficheLCP(unbuit_current_cp->getDfree(), unbuit_current_cp->getW(), unbuit_current_cp->getF(), unbuit_current_cp->getDimension());
    }


    current_cp = dynamic_cast<GenericConstraintProblem *>( unbuit_current_cp );

    return true;
}



// unbuilt Gauss Seidel: the
void UnbuiltGenericConstraintProblem::unbuiltGaussSeidel(double timeout, UnbuiltGenericConstraintSolver* solver)
{

    if(!dimension)
        return;



    double t0 = (double)CTime::getTime() ;
    double timeScale = 1.0 / (double)CTime::getTicksPerSec();

    double *dfree = getDfree();
    double *force = getF();
    double **w = getW();
    double tol = tolerance;

    double *d = _d.ptr();

    int i, j, l, nb;

    double errF[6];
    double error=0.0;

    bool convergence = false;
    sofa::helper::vector<double> tempForces;
    if(sor != 1.0) tempForces.resize(dimension);

    if(scaleTolerance && !allVerified)
        tol *= dimension;

    for(i=0; i<dimension; )
    {
        constraintsResolutions[i]->init(i, w, force);
        int nb = constraintsResolutions[i]->nbLines;
        // TODO : previous forces don't work with multiple constraints
//		if(cclist_elem1[i]) cclist_elem1[i]->setConstraintDForce(force, i, i+nb-1, true);
//		if(cclist_elem2[i]) cclist_elem2[i]->setConstraintDForce(force, i, i+nb-1, true);



        i += nb;


    }
    memset(force, 0, dimension * sizeof(double));	// Erase previous forces for the time being

    sofa::helper::vector<double>* graph_residuals = NULL;
    sofa::helper::vector<double> tabErrors;

    if(solver)
    {
        graph_residuals = &(*solver->graphErrors.beginEdit())["Error"];
        graph_residuals->clear();

        tabErrors.resize(dimension);
    }

    bool warning=false;

    for(i=0; i<this->maxIterations; i++)
    {


        bool constraintsAreVerified = true;
        if(sor != 1.0)
        {
            for(j=0; j<dimension; j++)
                tempForces[j] = force[j];
        }

        error=0.0;
        for(j=0; j<dimension; ) // increment of j realized at the end of the loop
        {


            //1. nbLines provide the dimension of the constraint  (max=6)
            nb = constraintsResolutions[j]->nbLines;

            //1bis => if the compliance is null=> do not consider this contact...
            if(w[j][j]<1e-20)
            {
                warning=true;
                j += nb;
                continue;
            }

            //2. for each line we compute the actual value of d
            //   (a)d is set to dfree
            for(l=0; l<nb; l++)
            {
                errF[l] = force[j+l];
                d[j+l] = dfree[j+l];
            }
            //   (b) contribution of forces are added to d
            if(cclist_elem1[j]) cclist_elem1[j]->addConstraintDisplacement(d, j, j+nb-1);
            if(cclist_elem2[j]) cclist_elem2[j]->addConstraintDisplacement(d, j, j+nb-1);

            //3. the specific resolution of the constraint(s) is called




            constraintsResolutions[j]->resolution(j, w, d, force, dfree);

            //4. the error is measured (displacement due to the new resolution (i.e. due to the new force))
            double contraintError = 0.0;
            if(nb > 1)
            {
                for(l=0; l<nb; l++)
                {
                    double lineError = 0.0;
                    for (int m=0; m<nb; m++)
                    {
                        double dofError = w[j+l][j+m] * (force[j+m] - errF[m]);
                        lineError += dofError * dofError;
                    }
                    lineError = sqrt(lineError);
                    if(lineError > tol)
                        constraintsAreVerified = false;

                    contraintError += lineError;
                }
            }
            else
            {
                contraintError = fabs(w[j][j] * (force[j] - errF[0]));
                if(contraintError > tol)
                    constraintsAreVerified = false;
            }

            if(constraintsResolutions[j]->tolerance)
            {
                if(contraintError > constraintsResolutions[j]->tolerance)
                    constraintsAreVerified = false;
                contraintError *= tol / constraintsResolutions[j]->tolerance;
            }

            error += contraintError;
            if(solver)
                tabErrors[j] = contraintError;

            //5. the force is updated for the constraint corrections
            bool update = false;
            for(l=0; l<nb; l++)
                update |= (force[j+l] || errF[l]);

            if(update)
            {
                double tempF[6];
                for(l=0; l<nb; l++)
                {
                    tempF[l] = force[j+l];
                    force[j+l] -= errF[l]; // DForce
                }

                if(cclist_elem1[j]) cclist_elem1[j]->setConstraintDForce(force, j, j+nb-1, update);
                if(cclist_elem2[j]) cclist_elem2[j]->setConstraintDForce(force, j, j+nb-1, update);

                for(l=0; l<nb; l++)
                    force[j+l] = tempF[l];
            }

            j += nb;
        }

        if(solver)
            graph_residuals->push_back(error);

        if(sor != 1.0)
        {
            for(j=0; j<dimension; j++)
                force[j] = sor * force[j] + (1-sor) * tempForces[j];
        }

        double t1 = (double)CTime::getTime();
        double dt = (t1 - t0)*timeScale;

        if(timeout && dt > timeout)
        {
            if (warning)
                std::cerr<<"WARNING Some constraints were suppr because compliance was < 1e-20"<<std::endl;
            return;
        }
        else if(allVerified)
        {
            if(constraintsAreVerified)
            {
                convergence = true;
                break;
            }
        }
        else if(error < tol/* && i>0*/) // do not stop at the first iteration (that is used for initial guess computation)
        {
            convergence = true;
            break;
        }
    }
    if (warning)
        std::cerr<<" Some constraints were suppr because compliance was < 1e-20"<<std::endl;

    if(solver)
    {
        if(!convergence)
        {
            if(solver->f_printLog.getValue())
                solver->serr << "No convergence : error = " << error << solver->sendl;
            else
                solver->sout << "No convergence : error = " << error << solver->sendl;
        }
        else if(solver->displayTime.getValue())
            solver->sout<<" Convergence after " << i+1 << " iterations " << solver->sendl;
    }

    sofa::helper::AdvancedTimer::valSet("GS iterations", i+1);

    for(i=0; i<dimension; i += constraintsResolutions[i]->nbLines)
        constraintsResolutions[i]->store(i, force, convergence);

    if(solver)
    {
        solver->graphErrors.endEdit();

        sofa::helper::vector<double>& graph_constraints = (*solver->graphConstraints.beginEdit())["Constraints"];
        graph_constraints.clear();

        for(j=0; j<dimension; )
        {
            nb = constraintsResolutions[j]->nbLines;

            if(tabErrors[j])
                graph_constraints.push_back(tabErrors[j]);
            else if(constraintsResolutions[j]->tolerance)
                graph_constraints.push_back(constraintsResolutions[j]->tolerance);
            else
                graph_constraints.push_back(tol);

            j += nb;
        }
        solver->graphConstraints.endEdit();
    }



}

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(UnbuiltGenericConstraintSolver);

int UnbuiltGenericConstraintSolverClass = core::RegisterObject("A Generic Constraint Solver using the Linear Complementarity Problem formulation to solve Constraint based components")
.add< UnbuiltGenericConstraintSolver >();



} // namespace constraintset

} // namespace component

} // namespace sofa
