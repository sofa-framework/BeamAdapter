/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
#ifndef SOFA_PLUGIN_BEAM_ADAPTER_H
#define SOFA_PLUGIN_BEAM_ADAPTER_H

#include <sofa/helper/system/config.h>


#ifdef SOFA_BUILD_BEAMADAPTER
#  define SOFA_BEAMADAPTER_API SOFA_EXPORT_DYNAMIC_LIBRARY
#else
#  define SOFA_BEAMADAPTER_API SOFA_IMPORT_DYNAMIC_LIBRARY
#endif


#endif /* SOFA_PLUGIN_BEAM_ADAPTER_H */
