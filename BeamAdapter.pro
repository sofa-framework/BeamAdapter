load(sofa/pre)
defineAsPlugin(BeamAdapter)

TEMPLATE = lib
TARGET = BeamAdapter


DEFINES += SOFA_BUILD_BEAMADAPTER


HEADERS += \
		AdaptiveBeamConstraint.h              \
		AdaptiveBeamConstraint.inl            \
		AdaptiveBeamController.h              \
		AdaptiveBeamController.inl            \
		AdaptiveBeamForceFieldAndMass.h       \
		AdaptiveBeamForceFieldAndMass.inl     \
		AdaptiveBeamLengthConstraint.h        \
		AdaptiveBeamLengthConstraint.inl      \
		AdaptiveBeamMapping.h                 \
		AdaptiveBeamMapping.inl               \
		BaseRestShape.h                       \
		BaseRestShape.inl                     \
		BeamInterpolation.h                   \
		BeamInterpolation.inl                 \
		initBeamAdapter.h                     \
		InterventionalRadiologyController.h   \
		InterventionalRadiologyController.inl \
		MultiAdaptiveBeamMapping.h            \
		MultiAdaptiveBeamMapping.inl          \
		SteerableCatheter.h                   \
		SteerableCatheter.inl                 \
		SutureController.h                    \
		SutureController.inl                  \
		UnbuiltGenericConstraintSolver.h      \
		WireBeamInterpolation.h               \
		WireBeamInterpolation.inl             \
		WireRestShape.h                       \
		WireRestShape.inl                     \


SOURCES += \
		AdaptiveBeamConstraint.cpp            \
		AdaptiveBeamController.cpp            \
		AdaptiveBeamForceFieldAndMass.cpp     \
		AdaptiveBeamLengthConstraint.cpp      \
		AdaptiveBeamMapping.cpp               \
		BaseRestShape.cpp                     \
		BeamInterpolation.cpp                 \
		initBeamAdapter.cpp                   \
		InterventionalRadiologyController.cpp \
		MultiAdaptiveBeamMapping.cpp          \
		SteerableCatheter.cpp                 \
		SutureController.cpp                  \
		UnbuiltGenericConstraintSolver.cpp    \
		WireBeamInterpolation.cpp             \
		WireRestShape.cpp                     \


contains(DEFINES,SOFA_DEV2){
HEADERS +=  \
		AdaptiveBeamContactMapper.h\
		AdaptiveBeamContactMapper.inl\
		AdaptiveBeamFrictionContact.h \
		AdaptiveBeamFrictionContact.inl \
		MultiAdaptiveBeamContactMapper.h\
		MultiAdaptiveBeamContactMapper.inl

SOURCES +=  \
		AdaptiveBeamContactMapper.cpp\
		AdaptiveBeamFrictionContact.cpp\
		MultiAdaptiveBeamContactMapper.cpp

}


contains(DEFINES,SOFA_HAVE_SOFAEVE){
HEADERS +=  \
		ImplicitSurfaceAdaptiveConstraint.h   \
		ImplicitSurfaceAdaptiveConstraint.inl

SOURCES +=  \
		ImplicitSurfaceAdaptiveConstraint.cpp
}


README_FILE = BeamAdapter.txt

unix : QMAKE_POST_LINK = cp $$SRC_DIR/$$README_FILE $$LIB_DESTDIR
win32 : QMAKE_POST_LINK = copy \"$$toWindowsPath($$SRC_DIR/$$README_FILE)\" \"$$LIB_DESTDIR\"

load(sofa/post)
