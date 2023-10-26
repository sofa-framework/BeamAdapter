def createScene(rootnode):
    settings = rootnode.addChild('Settings')
    settings.addObject('RequiredPlugin', name='BeamAdapter')  # Needed to use components [BeamProjectionDifferenceMultiMapping]
    settings.addObject('RequiredPlugin', name='Sofa.Component.AnimationLoop')  # Needed to use components [FreeMotionAnimationLoop]
    settings.addObject('RequiredPlugin', name='Sofa.Component.Constraint.Lagrangian.Correction')  # Needed to use components [GenericConstraintCorrection]
    settings.addObject('RequiredPlugin', name='Sofa.Component.Constraint.Lagrangian.Solver')  # Needed to use components [GenericConstraintSolver]
    settings.addObject('RequiredPlugin', name='Sofa.Component.Constraint.Projective')  # Needed to use components [FixedConstraint]
    settings.addObject('RequiredPlugin', name='Sofa.Component.LinearSolver.Direct')  # Needed to use components [SparseLDLSolver]
    settings.addObject('RequiredPlugin', name='Sofa.Component.Mass')  # Needed to use components [UniformMass]  
    settings.addObject('RequiredPlugin', name='Sofa.Component.ODESolver.Backward')  # Needed to use components [EulerImplicitSolver]
    settings.addObject('RequiredPlugin', name='Sofa.Component.SolidMechanics.Spring')  # Needed to use components [RestShapeSpringsForceField]
    settings.addObject('RequiredPlugin', name='Sofa.Component.StateContainer')  # Needed to use components [MechanicalObject]
    settings.addObject('RequiredPlugin', name='Sofa.Component.Topology.Container.Dynamic')  # Needed to use components [EdgeSetTopologyContainer,PointSetTopologyContainer]
    settings.addObject('RequiredPlugin', name='Sofa.Component.Visual')  # Needed to use components [VisualStyle]  
    settings.addObject('RequiredPlugin', name='Sofa.GUI.Component')  # Needed to use components [AttachBodyButtonSetting]

    rootnode.addObject('VisualStyle', displayFlags='showBehavior showVisual')
    rootnode.addObject('AttachBodyButtonSetting', stiffness=0.1)
    rootnode.gravity.value = [0, -9810, 0]
    rootnode.dt.value = 0.01

    rootnode.addObject('FreeMotionAnimationLoop')
    rootnode.addObject('GenericConstraintSolver', maxIterations=1000, tolerance=1e-3)

    simulation = rootnode.addChild('Simulation')
    simulation.addObject('EulerImplicitSolver')
    simulation.addObject('SparseLDLSolver', template='CompressedRowSparseMatrixMat3x3d')
    simulation.addObject('GenericConstraintCorrection')

    # Beam
    nbEdges = 3
    beam = simulation.addChild('Beam')
    beam.addObject('EdgeSetTopologyContainer', edges=[[i, i + 1] for i in range(nbEdges)])
    beam.addObject('MechanicalObject', template='Rigid3',
                   position=[[100 * i, 0, 0, 0, 0, 0, 1] for i in range(nbEdges + 1)])
    beam.addObject('AdaptiveBeamForceFieldAndMass', massDensity=1e-6)
    beam.addObject('BeamInterpolation', straight=False, dofsAndBeamsAligned=False,
                   defaultYoungModulus=1e5)

    # Particles
    particles = simulation.addChild('Particles')
    particles.addObject('PointSetTopologyContainer')
    particles.addObject('MechanicalObject', template='Rigid3', showObject=True, showObjectScale=30, drawMode=1,
                        position=[[0, 0, 0, 0, 0, 0, 1], [80, 0, 0, 0, 0, 0, 1], [150, 0, 0, 0, 0, 0, 1]])
    particles.addObject('UniformMass', totalMass=0.005)
    particles.addObject('FixedConstraint', indices=[0, 2])  # Fix first and last particles

    # This will constrain the first particle and its projection on the beam to remain attached together
    fixing = particles.addChild('FixingConstraintParticle1')
    beam.addChild(fixing)
    fixing.addObject('MechanicalObject', template='Rigid3', position=[0, 0, 0, 0, 0, 0, 0])
    fixing.addObject('RestShapeSpringsForceField', stiffness=1e6, angularStiffness=1e6)
    fixing.addObject('BeamProjectionDifferenceMultiMapping',
                     directions=[1, 1, 1, 1, 1, 1, 1],  # The three positions and rotations
                     input1=particles.getMechanicalState().linkpath,
                     indicesInput1=[0],
                     input2=beam.getMechanicalState().linkpath,
                     interpolationInput2=beam.BeamInterpolation.linkpath,
                     output=fixing.getMechanicalState().linkpath)

    # This will constrain the second particle to slide along the beam, and it will also constrain the orientation
    # of the beam and the particle to remain the same
    sliding = particles.addChild('SlidingConstraintParticle2')
    beam.addChild(sliding)
    sliding.addObject('MechanicalObject', template='Rigid3', position=[0, 0, 0, 0, 0, 0, 0])
    sliding.addObject('RestShapeSpringsForceField', stiffness=1e6, angularStiffness=0)
    sliding.addObject('BeamProjectionDifferenceMultiMapping',
                      directions=[0, 1, 1, 1, 1, 1, 1],  # Only y, z positions to allow the particle to slide on the beam
                      # but this time we add the three rotations
                      input1=particles.getMechanicalState().linkpath,
                      indicesInput1=[1],
                      input2=beam.getMechanicalState().linkpath,
                      interpolationInput2=beam.BeamInterpolation.linkpath,
                      output=sliding.getMechanicalState().linkpath)

    # This will constrain the third particle to slide along the beam
    sliding = particles.addChild('SlidingConstraintParticle3')
    beam.addChild(sliding)
    sliding.addObject('MechanicalObject', template='Rigid3', position=[0, 0, 0, 0, 0, 0, 0])
    sliding.addObject('RestShapeSpringsForceField', stiffness=1e6, angularStiffness=1e6)
    sliding.addObject('BeamProjectionDifferenceMultiMapping',
                      directions=[0, 1, 1, 0, 0, 0, 0],  # Only y, z positions to allow the particle to slide on the beam
                      input1=particles.getMechanicalState().linkpath,
                      indicesInput1=[2],
                      input2=beam.getMechanicalState().linkpath,
                      interpolationInput2=beam.BeamInterpolation.linkpath,
                      output=sliding.getMechanicalState().linkpath)
