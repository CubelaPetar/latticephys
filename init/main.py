#!/usr/bin/env python 3

###################################################################################
#
# Title:  Oscillation Modes of a one dimensional molecule
# Author: Petar Cubela
# Date of creation: 28.12.2022
# Description:
#
#
####################################################################################

"""
Illustration of the eigenmodes of an one dimensional molecule.
"""

#------# Modules #--------#

from math import sqrt
import numpy as np
import numpy.linalg as la
import matplotlib.pylab as plt

#-------------------------#

#-----# Constants #-------#

SPRING_CONSTANT     = 1.0              # newton per meters
MASS_SMALL          = 1.0
MASS_LARGE          = 2.0
OMEGA_NULL_SQUARED  = SPRING_CONSTANT / MASS_SMALL
ALPHA               = sqrt( MASS_SMALL / MASS_LARGE )

#-------------------------#



def dynamical_matrix():
    matrix = np.array( [ [ OMEGA_NULL_SQUARED, - ALPHA * OMEGA_NULL_SQUARED , 0.0   ],
                         [ - ALPHA * OMEGA_NULL_SQUARED,  2 * ALPHA ** 2 * OMEGA_NULL_SQUARED,   - ALPHA * OMEGA_NULL_SQUARED ],
                         [ 0.0 ,- ALPHA * OMEGA_NULL_SQUARED, OMEGA_NULL_SQUARED ]] )
    return matrix

def eigenh(A):
    """ diagonalize a hermitian matrix and return its eigenvalues and vectors ordered by size """
    eigenValues, eigenVectors = la.eigh(A)
    idx = np.argsort(eigenValues)
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]

    return (eigenValues, eigenVectors)

def draw():
    """ Draw the band structure along specified symmetry directions. """

def disperion_along_direction():
    """ compute the band structure along a specific firection 
    args:
    list of values

    start_point
    end_point
    qx
    qy

    Lx
    Ly
    """
def main():
    print("==============================")
    print("Hello, world!")
    print("==============================")
    A = np.zeros((6,6))
    A[0:3, 0:3] = np.zeros((3,3))
    A[0:3, 3:6] = np.eye(3)
    A[3:6, 3:6] = np.zeros((3,3))
    A[3:6, 0:3] = dynamical_matrix()

    eigenValues = eigenh(A)[0]
    print(eigenValues)

if __name__ == '__main__':
    main()


