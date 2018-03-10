# UTRI3
Standard 2D 3-node linear triangle User Element for Abaqus/Standard (UEL Subroutine).

Element formulation:
1. Can be set for either plane strain or plane stress.
1. Can be set for either isotropic or lamina material behavior.

This module:
1. Must be included in the user source code with an INCLUDE statement.
1. Must be accessed in the user source code with a USE statement within the UEL user subroutine.

The module-level element subroutine:
1. Must be called in the user source code with a CALL statement within the UEL user subroutine.

## USAGE
All instructions are given in the `utri3.f90` source code.
