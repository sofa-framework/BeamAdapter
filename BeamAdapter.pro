load(sofa/pre)
defineAsPlugin(BeamAdapter)

TEMPLATE = lib
TARGET = BeamAdapter


DEFINES += SOFA_BUILD_BEAMADAPTER


HEADERS += \
		AdaptiveBeamConstraint.h              \
		AdaptiveBeamConstraint.inl            \
		AdaptiveBeamContactMapper.h           \
		AdaptiveBeamContactMapper.inl         \
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
		ImplicitSurfaceAdaptiveConstraint.h   \
		ImplicitSurfaceAdaptiveConstraint.inl \
		initBeamAdapter.h                     \
		InterventionalRadiologyController.h   \
		InterventionalRadiologyController.inl \
		MultiAdaptiveBeamContactMapper.h      \
		MultiAdaptiveBeamContactMapper.inl    \
		MultiAdaptiveBeamMapping.h            \
		MultiAdaptiveBeamMapping.inl          \
		SutureController.h                    \
		SutureController.inl                  \
		UnbuiltGenericConstraintSolver.h      \
		WireBeamInterpolation.h               \
		WireBeamInterpolation.inl             \
		WireRestShape.h                       \
		WireRestShape.inl                     \


SOURCES += \
		AdaptiveBeamConstraint.cpp            \
		AdaptiveBeamContactMapper.cpp         \
		AdaptiveBeamController.cpp            \
		AdaptiveBeamForceFieldAndMass.cpp     \
		AdaptiveBeamFrictionContact.cpp       \
		AdaptiveBeamLengthConstraint.cpp      \
		AdaptiveBeamMapping.cpp               \
		BaseRestShape.cpp                     \
		BeamInterpolation.cpp                 \
		ImplicitSurfaceAdaptiveConstraint.cpp \
		initBeamAdapter.cpp                   \
		InterventionalRadiologyController.cpp \
		MultiAdaptiveBeamContactMapper.cpp    \
		MultiAdaptiveBeamMapping.cpp          \
		SutureController.cpp                  \
		UnbuiltGenericConstraintSolver.cpp    \
		WireBeamInterpolation.cpp             \
		WireRestShape.cpp                     \


README_FILE = BeamAdapter.txt

unix : QMAKE_POST_LINK = cp $$SRC_DIR/$$README_FILE $$LIB_DESTDIR
win32 : QMAKE_POST_LINK = copy \"$$toWindowsPath($$SRC_DIR/$$README_FILE)\" \"$$LIB_DESTDIR\"

load(sofa/post)
