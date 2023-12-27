Beam Adapter Documentation
=======================
Contributors: See https://github.com/sofa-framework/BeamAdapter/blob/master/Authors.md

### Introduction

The plugin Beam Adpater allows to have a dynamic implementation of the FEM beam elements.
The underlying mechanical fundations are based on beam theory (https://en.wikipedia.org/wiki/Euler–Bernoulli_beam_theory), that have been extended to large transformation using corotational formualation.
These elements are particularly useful for modeling 1D like deformable structures (cables, threads, flexible needles, catheters, etc...).


This tutorial presents the mechanical and implementation basis of the plugin Beam Adpater.
It also presents some examples of use and lists the current limitations.

#### This Documentation covers:
- [Mechanical basis](docs/modeling.md)

## Limitation and ongoing work on the plugin


## Implementation
#Shape function: BeamInterpolation#
In the following, we describe the main components of the implementation of BeamAdapter Plugin
One of the central component of the beams is the shape function. 
The description of the shape function relies on:
- One topology of edges
- A set of frames with Dofs (3 translations, 3 rotations) that corresponds to the points of the topology
- A spline (3d order) support for each edge

## Examples
