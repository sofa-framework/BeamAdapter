/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
#define SOFA_PLUGIN_BEAMADAPTER_FORCEFEEDBACKSPRINGFORCEFIELD_CPP

#include "ForceFeedbackSpringsForceField.inl"
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;


SOFA_DECL_CLASS(ForceFeedbackSpringsForceField)

///////////// SPECIALIZATION FOR RIGID TYPES //////////////


#ifndef SOFA_FLOAT

template<>
void ForceFeedbackSpringsForceField<Rigid3dTypes>::addForce(const core::MechanicalParams* /* mparams */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& /* v */)
{
    sofa::helper::WriteAccessor< DataVecDeriv > f1 = f;
    sofa::helper::ReadAccessor< DataVecCoord > p1 = x;

    sofa::helper::ReadAccessor< DataVecCoord > p0 = *this->getExtPosition();
    Deriv& forceFeedback = *d_forceFeedback.beginEdit();
    forceFeedback.clear();

    f1.resize(p1.size());

    if (this->recompute_indices.getValue())
    {
        this->recomputeIndices();
    }

    const VecReal& k = this->stiffness.getValue();
    const VecReal& k_a = this->angularStiffness.getValue();

    for (unsigned int i = 0; i < this->m_indices.size(); i++)
    {
        const unsigned int index = this->m_indices[i];
        unsigned int ext_index = this->m_indices[i];
        if(this->local_useRestMState)
            ext_index= this->m_ext_indices[i];

        // translation
        if (i >= this->m_pivots.size())
        {
            Vec3d dx = p1[index].getCenter() - p0[ext_index].getCenter();
            getVCenter(f1[index]) -=  dx * (i < k.size() ? k[i] : k[0]) ;
            getVCenter(forceFeedback) -=  dx * (i < k.size() ? k[i] : k[0]) ;
        }
        else
        {
            CPos localPivot = p0[ext_index].getOrientation().inverseRotate(this->m_pivots[i] - p0[ext_index].getCenter());
            CPos rotatedPivot = p1[index].getOrientation().rotate(localPivot);
            CPos pivot2 = p1[index].getCenter() + rotatedPivot;
            CPos dx = pivot2 - this->m_pivots[i];
            getVCenter(f1[index]) -= dx * (i < k.size() ? k[i] : k[0]) ;
            getVCenter(forceFeedback) -= dx * (i < k.size() ? k[i] : k[0]) ;

        }

        // rotation
        Quatd dq = p1[index].getOrientation() * p0[ext_index].getOrientation().inverse();
        Vec3d dir;
        double angle=0;
        dq.normalize();

        if (dq[3] < 0)
        {
            dq = dq * -1.0;
        }

        if (dq[3] < 0.999999999999999)
            dq.quatToAxis(dir, angle);

        getVOrientation(f1[index]) -= dir * angle * (i < k_a.size() ? k_a[i] : k_a[0]);
        getVOrientation(forceFeedback) -= dir * angle * (i < k_a.size() ? k_a[i] : k_a[0]);
    }
    d_forceFeedback.endEdit();
}


#endif // SOFA_FLOAT

#ifndef SOFA_DOUBLE

template<>
void ForceFeedbackSpringsForceField<Rigid3fTypes>::addForce(const core::MechanicalParams* /* mparams */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& /* v */)
{
    sofa::helper::WriteAccessor< DataVecDeriv > f1 = f;
    sofa::helper::ReadAccessor< DataVecCoord > p1 = x;

    sofa::helper::ReadAccessor< DataVecCoord > p0 = *getExtPosition();
    Deriv& forceFeedback = *d_forceFeedback.beginEdit();
    forceFeedback.clear();

    f1.resize(p1.size());

    if (recompute_indices.getValue())
    {
        recomputeIndices();
    }

    const VecReal& k = stiffness.getValue();
    const VecReal& k_a = angularStiffness.getValue();

    for (unsigned int i=0; i<m_indices.size(); i++)
    {
        const unsigned int index = m_indices[i];
        const unsigned int ext_index = m_ext_indices[i];

        // translation
        Vec3f dx = p1[index].getCenter() - p0[ext_index].getCenter();
        getVCenter(f1[index]) -=  dx * (i < k.size() ? k[i] : k[0]) ;
        getVCenter(forceFeedback) -=  dx * (i < k.size() ? k[i] : k[0]) ;

        // rotation
        Quatf dq = p1[index].getOrientation() * p0[ext_index].getOrientation().inverse();
        Vec3f dir;
        Real angle=0;
        dq.normalize();
        if (dq[3] < 0.999999999999999)
            dq.quatToAxis(dir, angle);
        dq.quatToAxis(dir, angle);

        //std::cout<<"dq : "<<dq <<"  dir :"<<dir<<"  angle :"<<angle<<std::endl;
        getVOrientation(f1[index]) -= dir * angle * (i < k_a.size() ? k_a[i] : k_a[0]) ;
        getVOrientation(forceFeedback) -= dir * angle * (i < k_a.size() ? k_a[i] : k_a[0]) ;
    }
    std::cout<<forceFeedback<<std::endl;
    d_forceFeedback.endEdit();
}
#endif // SOFA_FLOAT



int ForceFeedbackSpringsForceFieldClass = core::RegisterObject("Simple elastic springs applied to given degrees of freedom between their current and rest shape position")
#ifndef SOFA_FLOAT
        .add< ForceFeedbackSpringsForceField<Vec3dTypes> >()
        .add< ForceFeedbackSpringsForceField<Vec1dTypes> >()
        .add< ForceFeedbackSpringsForceField<Rigid3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< ForceFeedbackSpringsForceField<Vec3fTypes> >()
        .add< ForceFeedbackSpringsForceField<Vec1fTypes> >()
        .add< ForceFeedbackSpringsForceField<Rigid3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<Vec3dTypes>;
template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<Vec1dTypes>;
template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<Vec3fTypes>;
template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<Vec1fTypes>;
template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<Rigid3fTypes>;
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa
