/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
#ifndef SOFA_PLUGIN_BEAMADAPTER_FORCEFEEDBACKSPRINGFORCEFIELD_INL
#define SOFA_PLUGIN_BEAMADAPTER_FORCEFEEDBACKSPRINGFORCEFIELD_INL

#include "ForceFeedbackSpringsForceField.h"
#include <sofa/helper/system/config.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <SofaDeformable/RestShapeSpringsForceField.inl>

#include <assert.h>
#include <iostream>


namespace sofa
{

namespace component
{

namespace forcefield
{

template<class DataTypes>
ForceFeedbackSpringsForceField<DataTypes>::ForceFeedbackSpringsForceField()
    : d_forceFeedback(initData(&d_forceFeedback, "forceFeedback", "force that can be used for force feedback"))
{
}

template<class DataTypes>
void ForceFeedbackSpringsForceField<DataTypes>::init()
{
    sofa::component::forcefield::RestShapeSpringsForceField<DataTypes>::init();

    Deriv & forceFeedback = *d_forceFeedback.beginEdit();
    forceFeedback.clear();
    d_forceFeedback.endEdit();
}

template<class DataTypes>
void ForceFeedbackSpringsForceField<DataTypes>::bwdInit()
{
    sofa::component::forcefield::RestShapeSpringsForceField<DataTypes>::bwdInit();

    local_useRestMState = true;
}



template<class DataTypes>
void ForceFeedbackSpringsForceField<DataTypes>::addForce(const core::MechanicalParams*  mparams , DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv&  v )
{
    Deriv & forceFeedback = *d_forceFeedback.beginEdit();
    forceFeedback.clear();

    SOFA_UNUSED(mparams);
    SOFA_UNUSED(v);

    helper::WriteAccessor< DataVecDeriv > f1 = f;
    helper::ReadAccessor< DataVecCoord > p1 = x;
    helper::ReadAccessor< DataVecCoord > p0 = *this->getExtPosition();

    f1.resize(p1.size());

    if (this->recompute_indices.getValue())
    {
        this->recomputeIndices();
    }

    //Springs_dir.resize(m_indices.size() );
    if ( this->k.size()!= this->m_indices.size() )
    {
        const Real k0 = this->k[0];

        for (unsigned int i=0; i<this->m_indices.size(); i++)
        {
            const unsigned int index = this->m_indices[i];

            unsigned int ext_index = this->m_indices[i];
            if(local_useRestMState)
                ext_index= this->m_ext_indices[i];

            Deriv dx = p1[index] - p0[ext_index];
            f1[index] -=  dx * k0 ;
            forceFeedback -=  dx * k0 ;
        }
    }
    else
    {
        for (unsigned int i=0; i<this->m_indices.size(); i++)
        {
            const unsigned int index = this->m_indices[i];
            unsigned int ext_index = this->m_indices[i];
            if(local_useRestMState)
                ext_index= this->m_ext_indices[i];

            Deriv dx = p1[index] - p0[ext_index];
            f1[index] -=  dx * this->k[i];
            forceFeedback -=  dx * this->k[i];
        }
    }
    d_forceFeedback.endEdit();
}


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_PLUGIN_BEAMADAPTER_FORCEFEEDBACKSPRINGFORCEFIELD_INL



