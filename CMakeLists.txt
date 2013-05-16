cmake_minimum_required(VERSION 2.8)

project("BeamAdapter")

include(${SOFA_CMAKE_DIR}/preProject.cmake)

set(HEADER_FILES

	AdaptiveBeamConstraint.h           
	AdaptiveBeamConstraint.inl         
	AdaptiveBeamContactMapper.h         
	AdaptiveBeamContactMapper.inl         
	AdaptiveBeamController.h             
	AdaptiveBeamController.inl           
	AdaptiveBeamForceFieldAndMass.h       
	AdaptiveBeamForceFieldAndMass.inl    
	AdaptiveBeamLengthConstraint.h       
	AdaptiveBeamLengthConstraint.inl     
	AdaptiveBeamMapping.h               
	AdaptiveBeamMapping.inl             
	BaseRestShape.h                      
	BaseRestShape.inl                    
	BeamInterpolation.h                  
	BeamInterpolation.inl                
	ImplicitSurfaceAdaptiveConstraint.h   
	ImplicitSurfaceAdaptiveConstraint.inl 
	initBeamAdapter.h                     
	InterventionalRadiologyController.h   
	InterventionalRadiologyController.inl 
	MultiAdaptiveBeamContactMapper.h      
	MultiAdaptiveBeamContactMapper.inl    
	MultiAdaptiveBeamMapping.h            
	MultiAdaptiveBeamMapping.inl          
	SutureController.h                    
	SutureController.inl                  
	UnbuiltGenericConstraintSolver.h      
	WireBeamInterpolation.h               
	WireBeamInterpolation.inl             
	WireRestShape.h                       
	WireRestShape.inl  
	SteerableCatheter.h
	SteerableCatheter.inl

	)

set(SOURCE_FILES

	AdaptiveBeamConstraint.cpp            
	AdaptiveBeamContactMapper.cpp         
	AdaptiveBeamController.cpp            
	AdaptiveBeamForceFieldAndMass.cpp     
	AdaptiveBeamFrictionContact.cpp       
	AdaptiveBeamLengthConstraint.cpp      
	AdaptiveBeamMapping.cpp               
	BaseRestShape.cpp                     
	BeamInterpolation.cpp                 
	ImplicitSurfaceAdaptiveConstraint.cpp 
	initBeamAdapter.cpp                   
	InterventionalRadiologyController.cpp 
	MultiAdaptiveBeamContactMapper.cpp    
	MultiAdaptiveBeamMapping.cpp          
	SutureController.cpp                  
	UnbuiltGenericConstraintSolver.cpp    
	WireBeamInterpolation.cpp             
	WireRestShape.cpp                
	SteerableCatheter.cpp

	)


set(README_FILE "BeamAdapter.txt" )	

#message(status "BeamAdapter")
	
	
if(SOFA-PLUGIN_SOFAADVANCEDCONSTRAINT)
	#find_package(CUDA REQUIRED)
	#cuda_include_directories(${SOFA_CUDA_DIR})
	#cuda_include_directories(${SOFA_FRAMEWORK_DIR})
	#cuda_add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${CUDA_SOURCES})
	add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} )
else()
	add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} )
endif()

set(COMPILER_DEFINES "SOFA_BUILD_BEAMADAPTER")
set(LINKER_DEPENDENCIES SofaAdvancedConstraint SofaComponentAdvanced SofaVolumetricData )

if(SOFA-PLUGIN_SOFAADVANCEDCONSTRAINT)
    #list(APPEND LINKER_DEPENDENCIES SofaCUDA)
    #list(APPEND LINKER_DEPENDENCIES SofaCUDADev)
endif()

#
# Organize projects into folders
#

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
SET_PROPERTY(TARGET ${PROJECT_NAME}	PROPERTY FOLDER "plugins")


include(${SOFA_CMAKE_DIR}/postProject.cmake)