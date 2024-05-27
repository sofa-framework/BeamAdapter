/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "BeamAdapter", "index.html", [
    [ "Beam Adapter Documentation", "index.html", [
      [ "1. Introduction", "index.html#autotoc_md2", null ],
      [ "2. Documentation cover:", "index.html#autotoc_md3", null ],
      [ "3. Technical roadmap", "index.html#autotoc_md4", [
        [ "AdaptiveBeamMapping", "index.html#autotoc_md5", null ],
        [ "TO ADD List:", "index.html#autotoc_md6", null ]
      ] ],
      [ "4. Examples list:", "index.html#autotoc_md7", [
        [ "Simple scenario", "index.html#autotoc_md8", null ],
        [ "Component examples", "index.html#autotoc_md9", null ],
        [ "Complex scenario", "index.html#autotoc_md10", null ]
      ] ]
    ] ],
    [ "theory", "md_doc_modeling_theory.html", [
      [ "Modeling Theory", "md_doc_modeling_theory.html#autotoc_md11", [
        [ "1. Corotational beam model", "md_doc_modeling_theory.html#autotoc_md12", null ],
        [ "2. B-splines , Bézier Splines", "md_doc_modeling_theory.html#autotoc_md13", null ]
      ] ],
      [ "Implementation", "md_doc_modeling_theory.html#autotoc_md14", [
        [ "WireRestShape", "md_doc_modeling_theory.html#autotoc_md15", null ],
        [ "WireBeamInterpolation", "md_doc_modeling_theory.html#autotoc_md16", null ],
        [ "AdaptiveBeamForceFieldAndMass", "md_doc_modeling_theory.html#autotoc_md17", null ],
        [ "AdaptiveBeamMapping", "md_doc_modeling_theory.html#autotoc_md18", null ]
      ] ],
      [ "References", "md_doc_modeling_theory.html#autotoc_md19", null ]
    ] ],
    [ "Scene Implementation", "md_doc_modeling_implementation.html", [
      [ "All BeamAdapter components", "md_doc_modeling_implementation.html#autotoc_md21", [
        [ "1. The Components on Root Node", "md_doc_modeling_implementation.html#autotoc_md22", [
          [ "RequiredPlugin", "md_doc_modeling_implementation.html#autotoc_md23", null ],
          [ "VisualStyle", "md_doc_modeling_implementation.html#autotoc_md24", null ],
          [ "FreeMotionAnimationLoop", "md_doc_modeling_implementation.html#autotoc_md25", null ],
          [ "LCPConstraintSolver", "md_doc_modeling_implementation.html#autotoc_md26", null ],
          [ "CollisionPipeline", "md_doc_modeling_implementation.html#autotoc_md27", null ],
          [ "BruteForceDetection", "md_doc_modeling_implementation.html#autotoc_md28", null ],
          [ "LocalMinDistance", "md_doc_modeling_implementation.html#autotoc_md29", null ],
          [ "CollisionResponse", "md_doc_modeling_implementation.html#autotoc_md30", null ],
          [ "CollisionGroup", "md_doc_modeling_implementation.html#autotoc_md31", null ]
        ] ],
        [ "2. The Components for the shape of the catheter", "md_doc_modeling_implementation.html#autotoc_md32", [
          [ "SteerableCatheter", "md_doc_modeling_implementation.html#autotoc_md33", null ],
          [ "EdgeSetTopologyContainer", "md_doc_modeling_implementation.html#autotoc_md34", null ],
          [ "EdgeSetTopologyModifier, EdgeSetGeometryAlgorithms", "md_doc_modeling_implementation.html#autotoc_md35", null ],
          [ "MechanicalObject", "md_doc_modeling_implementation.html#autotoc_md36", null ]
        ] ],
        [ "3. The Components for simulation of the catheter", "md_doc_modeling_implementation.html#autotoc_md37", [
          [ "EulerImplicit", "md_doc_modeling_implementation.html#autotoc_md38", null ],
          [ "BTDLinearSolver", "md_doc_modeling_implementation.html#autotoc_md39", null ],
          [ "RegularGrid", "md_doc_modeling_implementation.html#autotoc_md40", null ],
          [ "MechanicalObject", "md_doc_modeling_implementation.html#autotoc_md41", null ],
          [ "WireBeamInterpolation", "md_doc_modeling_implementation.html#autotoc_md42", null ],
          [ "AdaptiveBeamForceFieldAndMass", "md_doc_modeling_implementation.html#autotoc_md43", null ],
          [ "InterventionalRadiologyController", "md_doc_modeling_implementation.html#autotoc_md44", null ],
          [ "LinearSolverConstraintCorrection", "md_doc_modeling_implementation.html#autotoc_md45", null ],
          [ "FixedConstraint", "md_doc_modeling_implementation.html#autotoc_md46", null ]
        ] ],
        [ "4. Collision", "md_doc_modeling_implementation.html#autotoc_md47", [
          [ "EdgeSetTopologyContainer, EdgeSetTopologyModifier", "md_doc_modeling_implementation.html#autotoc_md48", null ],
          [ "MechanicalObject", "md_doc_modeling_implementation.html#autotoc_md49", null ],
          [ "MultiAdaptiveBeamMapping", "md_doc_modeling_implementation.html#autotoc_md50", null ],
          [ "Point / Line", "md_doc_modeling_implementation.html#autotoc_md51", null ]
        ] ],
        [ "5. Visualization", "md_doc_modeling_implementation.html#autotoc_md52", [
          [ "Topology Components: QuadSetTopologyContainer, QuadSetTopologyModifier, etc...", "md_doc_modeling_implementation.html#autotoc_md53", null ],
          [ "MechanicalObject", "md_doc_modeling_implementation.html#autotoc_md54", null ],
          [ "AdaptiveBeamMapping", "md_doc_modeling_implementation.html#autotoc_md55", null ],
          [ "OglModel", "md_doc_modeling_implementation.html#autotoc_md56", null ],
          [ "IdentityMapping", "md_doc_modeling_implementation.html#autotoc_md57", null ]
        ] ]
      ] ],
      [ "References", "md_doc_modeling_implementation.html#autotoc_md58", null ]
    ] ],
    [ "Todo List", "todo.html", null ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", null ],
        [ "Variables", "namespacemembers_vars.html", null ]
      ] ]
    ] ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", null ],
        [ "Variables", "functions_vars.html", null ],
        [ "Typedefs", "functions_type.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"_adaptive_beam_contact_mapper_8cpp.html#a3d2088e3645f6e64b21c3b1cab13b761",
"classsofa_1_1component_1_1constraint_1_1__implicitsurfaceadaptiveconstraint___1_1_implicit_surface_adaptive_constraint.html#a1ee497138540ec7537efb26cc94a823c",
"classsofa_1_1component_1_1controller_1_1__interventionalradiologycontroller___1_1_interventional_radiology_controller.html#a8eabd560d3832185c58f5324d2c72c84",
"classsofa_1_1component_1_1fem_1_1__beaminterpolation___1_1_beam_interpolation.html#aebfcdbbd0fe4efbff10b180be50ca6bc",
"classsofa_1_1component_1_1mapping_1_1__adaptivebeammapping___1_1_adaptive_beam_mapping.html#aed2bbb4982de3add6ef6820d4ae5f2d6",
"structsofa_1_1beamadapter_1_1_beam_section.html#a3c7e2592540bfbfddb22de4d0bc0f4c3"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';