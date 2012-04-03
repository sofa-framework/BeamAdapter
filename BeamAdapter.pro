load(sofa/pre)
defineAsPlugin(BeamAdapter)

TEMPLATE = lib
TARGET = BeamAdapter


DEFINES += SOFA_BUILD_BEAMADAPTER


HEADERS += \
		initBeamAdapter.h                     \
		WireRestShape.h                       \
		WireRestShape.inl                     \
                BaseRestShape.h                       \
		BaseRestShape.inl                     \
		BeamInterpolation.h                   \
		BeamInterpolation.inl                 \
		WireBeamInterpolation.h               \
		WireBeamInterpolation.inl             \
		AdaptiveBeamForceFieldAndMass.h       \
		AdaptiveBeamForceFieldAndMass.inl     \
		AdaptiveBeamController.h              \
		AdaptiveBeamController.inl            \
		InterventionalRadiologyController.h   \
		InterventionalRadiologyController.inl \
		AdaptiveBeamMapping.h                 \
                AdaptiveBeamMapping.inl               \
		AdaptiveBeamContactMapper.h\
		AdaptiveBeamContactMapper.inl\
		MultiAdaptiveBeamContactMapper.h\
		MultiAdaptiveBeamContactMapper.inl\
		MultiAdaptiveBeamMapping.h            \
		MultiAdaptiveBeamMapping.inl           	\
                SutureController.h			\
                SutureController.inl		\
		AdaptiveBeamConstraint.h 		\
                AdaptiveBeamConstraint.inl             \
AdaptiveBeamLengthConstraint.h \
AdaptiveBeamLengthConstraint.inl \
ImplicitSurfaceAdaptiveConstraint.h \
ImplicitSurfaceAdaptiveConstraint.inl


SOURCES += \
		initBeamAdapter.cpp                   \
		WireRestShape.cpp                     \
				BaseRestShape.cpp                     \
		WireBeamInterpolation.cpp             \
		BeamInterpolation.cpp                 \
		AdaptiveBeamForceFieldAndMass.cpp     \
		AdaptiveBeamController.cpp            \
		InterventionalRadiologyController.cpp \
		AdaptiveBeamMapping.cpp               \   
AdaptiveBeamContactMapper.cpp\
AdaptiveBeamFrictionContact.cpp\
MultiAdaptiveBeamContactMapper.cpp\
		MultiAdaptiveBeamMapping.cpp 		\
                SutureController.cpp	\
                AdaptiveBeamConstraint.cpp          \
AdaptiveBeamLengthConstraint.cpp \
ImplicitSurfaceAdaptiveConstraint.cpp

	  
		  
README_FILE = BeamAdapter.txt

unix : QMAKE_POST_LINK = cp $$SRC_DIR/$$README_FILE $$LIB_DESTDIR
win32 : QMAKE_POST_LINK = copy \"$$toWindowsPath($$SRC_DIR/$$README_FILE)\" \"$$LIB_DESTDIR\"

load(sofa/post)
