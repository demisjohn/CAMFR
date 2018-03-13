# CAMFR
Forked from Sourceforge project for maintenance
Originally written by Peter Bienstman

Better maintained fork by [demisjohn](https://github.com/demisjohn/CAMFR)

# Description from [SourceForge](http://camfr.sourceforge.net/)
## Introduction

CAMFR (CAvity Modelling FRamework) is a fast, flexible, friendly full-vectorial Maxwell solver. Its main focus is on applications in the field of nanophotonics, like
wavelength-scale microstructures (like photonic crystal devices)
lasers (like vertical-cavity surface-emitting lasers)
light-emitting diodes (like resonant-cavity LEDs)
It is based on a combination of eigenmode expansion and advanced boundary conditions like perfectly matched layers (PML).

## Features

CAMFR is an ongoing active research project, started at the photonics group of the Department of Information Technology (INTEC) at Ghent University in Belgium. This means that it contains many attractive features and algorithms currently not yet found in commercial modelling tools. CAMFR can be used to calculate
- the scattering matrix of a structure
- the field inside a structure, for any given excitation
- band diagrams of an infinite periodic structure
-threshold material gain and resonance wavelength of laser modes
- the response to a current source in an arbitrary cavity
- structures terminated by a semi-infinite repetition of another structure

This functionality is currently available for two types of geometries:
- 2D Cartesian structures
- 3D cylindrical symmetric structures

Additionally, there is code to model the extraction from light emitting diodes, either planar devices, or 3D devices which incorporate 2D periodic structures.

Defining structures is quite straightforward, either layer-by-layer, or using geometric primitive shapes. There are also integrated plotting routines for rapid simulation feedback.

## Framework character
CAMFR is conceived as a C++ framework, with all the algorithms implemented in terms of abstract waveguides and scatterers. This makes it extremely easy to extend CAMFR to new geometries.
The end user does not deal with this C++ code directly, but rather through bindings to the Python scripting language. This makes the code very clear and flexible, and allows e.g. to seamlessly integrate CAMFR with Python-aware visualistion tools.

## License and support
All the code is released under the GPL.
