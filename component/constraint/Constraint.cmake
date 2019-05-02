list(APPEND HEADER_FILES
  component/constraint/AdaptiveBeamSlidingConstraint.h
  component/constraint/AdaptiveBeamSlidingConstraint.inl
  component/constraint/AdaptiveBeamLengthConstraint.h
  component/constraint/AdaptiveBeamLengthConstraint.inl
  )

list(APPEND SOURCE_FILES
  component/constraint/AdaptiveBeamSlidingConstraint.cpp
  component/constraint/AdaptiveBeamLengthConstraint.cpp
  )
  
if(SofaImplicitField_FOUND)
    list(APPEND HEADER_FILES
      component/constraint/ImplicitSurfaceAdaptiveConstraint.h
      component/constraint/ImplicitSurfaceAdaptiveConstraint.inl
    )

    list(APPEND SOURCE_FILES
      component/constraint/ImplicitSurfaceAdaptiveConstraint.cpp
    )
endif()

