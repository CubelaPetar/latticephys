
#!/usr/bin/env python 3

###################################################################################
#
# Title:  Atom in lattices class
# Author: Petar Cubela
# Date of creation: 28.12.2022
# Description:
#
#
####################################################################################

"""
Illustration of the eigenmodes of molecules or lattices.
"""

#------# Modules #--------#

#import math
import pygame

#-------------------------#

#-----# Constants #-------#

pygame.init()

WIDTH, HEIGHT = 800, 800
WIN           = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Planet Simulation")

WHITE     = (255, 255, 255)
YELLOW    = (255, 255, 0)
BLUE      = (100,149, 237)
RED       = (188, 39, 50)
DARK_GREY = (80, 78, 81)

FONT      = pygame.font.SysFont("comicsans", 20)



#-------------------------#

class Atoms():
    """ Create all atoms coupled by springs a molecule or lattice"""

    TIMESTEP = .01       
    SCALE    = 100

    def __init__(self, x: float, y: float, mass: float, spring_const: float, color: tuple, radius: float):
        self.x              = x
        self.y              = y

        self.mass           = mass
        self.spring_const   = spring_const

        self.color          = color
        self.radius         = radius

        self.trajectory     = []

        self.x_vel           = 0.0
        self.y_vel           = 0.0

    def draw(self, win):
        x                    = self.x * self.SCALE + WIDTH /2
        y                    = self.y * self.SCALE + HEIGHT /2

        if len(self.trajectory) > 2:
            updated_points      = []
            for point in self.trajectory:
                x, y             = point
                x                = x * self.SCALE + WIDTH /2
                y                = y * self.SCALE + HEIGHT /2
                updated_points.append((x, y))

            #pygame.draw.lines(win, self.color, False, updated_points, 2)

        pygame.draw.circle(win, self.color, (x, y), self.radius)

    def attraction(self, other):
        other_x, other_y     = other.x, other.y
        distance_x           = self.x - other_x 
        distance_y           = self.y - other_y

        force_x              = 0.
        force_y              = 0.
    
        force_x              = - self.spring_const * distance_x
        force_y              = - self.spring_const * distance_y
        return force_x, force_y

    def update_position(self, atoms) -> None:
        total_fx = total_fy  = 0.0
        for atom in atoms:
            if self == atom:
                continue

            fx, fy              = self.attraction(atom)
            total_fx           += fx
            total_fy           += fy

        self.x_vel          += total_fx / self.mass * self.TIMESTEP
        self.y_vel          += total_fy / self.mass * self.TIMESTEP

        self.x += self.x_vel * self.TIMESTEP
        self.y += self.y_vel * self.TIMESTEP
        self.trajectory.append((self.x, self.y))


def main():
    run           = True
    clock         = pygame.time.Clock()


    small_mass1   = Atoms(-1, 0, 1, 1,  BLUE, 20 )
    small_mass1.x_vel   = -.1

    large_mass    = Atoms( 0, 0, 2, 1,  RED, 30 )
    large_mass.x_vel   = 0

    small_mass2   = Atoms(+1, 0, 1, 1,  BLUE, 20 )
    small_mass2.x_vel   = .1



    atoms       = [ small_mass1, large_mass, small_mass2 ]

    while run:
        clock.tick(60)
        WIN.fill((0, 0, 0))
        #pygame.display.update()

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False

        for atom in atoms:
            atom.update_position(atoms)
            atom.draw(WIN)

        pygame.display.update()

    pygame.quit()



if __name__ == '__main__':
    main()
