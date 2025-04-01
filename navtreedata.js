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
    [ "Beam Adapter Documentation", "index.html", "index" ],
    [ "theory", "md_doc_2modeling_2theory.html", [
      [ "Modeling Theory", "md_doc_2modeling_2theory.html#autotoc_md11", [
        [ "1. Corotational beam model", "md_doc_2modeling_2theory.html#autotoc_md12", null ],
        [ "2. B-splines , Bézier Splines", "md_doc_2modeling_2theory.html#autotoc_md13", null ]
      ] ],
      [ "Implementation", "md_doc_2modeling_2theory.html#autotoc_md14", [
        [ "WireRestShape", "md_doc_2modeling_2theory.html#autotoc_md15", null ],
        [ "WireBeamInterpolation", "md_doc_2modeling_2theory.html#autotoc_md16", null ],
        [ "AdaptiveBeamForceFieldAndMass", "md_doc_2modeling_2theory.html#autotoc_md17", null ],
        [ "AdaptiveBeamMapping", "md_doc_2modeling_2theory.html#autotoc_md18", null ]
      ] ],
      [ "References", "md_doc_2modeling_2theory.html#autotoc_md19", null ]
    ] ],
    [ "Scene Implementation", "md_doc_2modeling_2implementation.html", [
      [ "All BeamAdapter components", "md_doc_2modeling_2implementation.html#autotoc_md21", [
        [ "1. The Components on Root Node", "md_doc_2modeling_2implementation.html#autotoc_md22", [
          [ "RequiredPlugin", "md_doc_2modeling_2implementation.html#autotoc_md23", null ],
          [ "VisualStyle", "md_doc_2modeling_2implementation.html#autotoc_md24", null ],
          [ "FreeMotionAnimationLoop", "md_doc_2modeling_2implementation.html#autotoc_md25", null ],
          [ "LCPConstraintSolver", "md_doc_2modeling_2implementation.html#autotoc_md26", null ],
          [ "CollisionPipeline", "md_doc_2modeling_2implementation.html#autotoc_md27", null ],
          [ "BruteForceDetection", "md_doc_2modeling_2implementation.html#autotoc_md28", null ],
          [ "LocalMinDistance", "md_doc_2modeling_2implementation.html#autotoc_md29", null ],
          [ "CollisionResponse", "md_doc_2modeling_2implementation.html#autotoc_md30", null ],
          [ "CollisionGroup", "md_doc_2modeling_2implementation.html#autotoc_md31", null ]
        ] ],
        [ "2. The Components for the shape of the catheter", "md_doc_2modeling_2implementation.html#autotoc_md32", [
          [ "SteerableCatheter", "md_doc_2modeling_2implementation.html#autotoc_md33", null ],
          [ "EdgeSetTopologyContainer", "md_doc_2modeling_2implementation.html#autotoc_md34", null ],
          [ "EdgeSetTopologyModifier, EdgeSetGeometryAlgorithms", "md_doc_2modeling_2implementation.html#autotoc_md35", null ],
          [ "MechanicalObject", "md_doc_2modeling_2implementation.html#autotoc_md36", null ]
        ] ],
        [ "3. The Components for simulation of the catheter", "md_doc_2modeling_2implementation.html#autotoc_md37", [
          [ "EulerImplicit", "md_doc_2modeling_2implementation.html#autotoc_md38", null ],
          [ "BTDLinearSolver", "md_doc_2modeling_2implementation.html#autotoc_md39", null ],
          [ "RegularGrid", "md_doc_2modeling_2implementation.html#autotoc_md40", null ],
          [ "MechanicalObject", "md_doc_2modeling_2implementation.html#autotoc_md41", null ],
          [ "WireBeamInterpolation", "md_doc_2modeling_2implementation.html#autotoc_md42", null ],
          [ "AdaptiveBeamForceFieldAndMass", "md_doc_2modeling_2implementation.html#autotoc_md43", null ],
          [ "InterventionalRadiologyController", "md_doc_2modeling_2implementation.html#autotoc_md44", null ],
          [ "LinearSolverConstraintCorrection", "md_doc_2modeling_2implementation.html#autotoc_md45", null ],
          [ "FixedConstraint", "md_doc_2modeling_2implementation.html#autotoc_md46", null ]
        ] ],
        [ "4. Collision", "md_doc_2modeling_2implementation.html#autotoc_md47", [
          [ "EdgeSetTopologyContainer, EdgeSetTopologyModifier", "md_doc_2modeling_2implementation.html#autotoc_md48", null ],
          [ "MechanicalObject", "md_doc_2modeling_2implementation.html#autotoc_md49", null ],
          [ "MultiAdaptiveBeamMapping", "md_doc_2modeling_2implementation.html#autotoc_md50", null ],
          [ "Point / Line", "md_doc_2modeling_2implementation.html#autotoc_md51", null ]
        ] ],
        [ "5. Visualization", "md_doc_2modeling_2implementation.html#autotoc_md52", [
          [ "Topology Components: QuadSetTopologyContainer, QuadSetTopologyModifier, etc...", "md_doc_2modeling_2implementation.html#autotoc_md53", null ],
          [ "MechanicalObject", "md_doc_2modeling_2implementation.html#autotoc_md54", null ],
          [ "AdaptiveBeamMapping", "md_doc_2modeling_2implementation.html#autotoc_md55", null ],
          [ "OglModel", "md_doc_2modeling_2implementation.html#autotoc_md56", null ],
          [ "IdentityMapping", "md_doc_2modeling_2implementation.html#autotoc_md57", null ]
        ] ]
      ] ],
      [ "References", "md_doc_2modeling_2implementation.html#autotoc_md58", null ]
    ] ],
    [ "Todo List", "todo.html", null ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", null ],
        [ "Functions", "namespacemembers_func.html", null ],
        [ "Variables", "namespacemembers_vars.html", null ],
        [ "Enumerations", "namespacemembers_enum.html", null ]
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
"_adaptive_beam_controller_8h_source.html",
"classbeamadapter_1_1_base_rod_section_material.html#aa123ceac868b5e59ba5ce9c99ac559f4",
"classbeamadapter_1_1_wire_rest_shape.html#a7293515ade2b3aa4581c054e88d1953a"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';