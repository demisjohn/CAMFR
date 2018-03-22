# -*- coding: utf-8 -*-
"""

CAMFR example showing the solving and plotting of a simple multi-mode Silicon waveguide clad in SiO2.

Spyder `Run Configuration` should be either:
    Execute in Current Interpreter
        OR
        
    Execute in New Interpreter 
        AND
    Interact with shell after completion


@author: Demis D. John
@email: demis@ucsb.edu
"""

####################################################

''' Import the CAMFR module: '''
import camfr

from numpy import pi


####################################################

''' Some refractive indices ("n"), thicknesses and widths '''
wl = 1.550      # wavelength in microns

core = camfr.Material( 3.4 )   # Silicon, RIx
clad = camfr.Material( 1.446 ) # Thermal Oxide
air = camfr.Material( 1.0 )

tCore = 0.300   # core thickness in Microns
tClad = 2.0

wCore = 1.00    # core width
wSides = 2.0    # lateral width of simulation space



''' Setup some CAMFR simulation parameters '''

camfr.set_lambda(wl)    # all CAMFR functions/classes should be under the camfr object
lam0 = abs(camfr.get_lambda())
k0 = 2*pi/lam0  # propagation const.


''' Number of modes to find. Computation time scales logarithmically with the number of eigenmodes used. '''
camfr.set_N(4)



''' -- Construct Waveguide -- '''
''' Vertical, bottom-to-top '''
center = camfr.Slab( clad(tClad) + core(tCore) + clad(tClad) )
side = camfr.Slab( clad(2*tClad) + air(center.width() - 2*tClad)  )

''' Horizontal, left-to-right '''
wg = side(wSides) + center(wCore) + side(wSides)    


s = camfr.Section( wg, 16, 30 )    
'''
    Section(
        Waveguide, 
        <# plane waves in 1st stage estimation>, 
        <# 1D modes in each Slab during 2nd stage>
    )
    
    Num. of planewaves should improve mode ordering and finding
    Num. of Slab Modes should improve accuracy of final mode profiles
'''


## CAMFR mode-solver params
camfr.set_section_solver(camfr.L)
camfr.set_mode_correction(camfr.full) 
''' `full`: solve for all modes.  Other options include `guided_only`, `none` don't refine the modes after the initial estimate. See CAMFR.pdf page 52.'''

camfr.set_left_wall(camfr.E_wall)
camfr.set_right_wall(camfr.E_wall)
camfr.set_lower_wall(camfr.slab_E_wall)
camfr.set_upper_wall(camfr.slab_E_wall)

''' Perfectly-matched layers @ simulation edges, to simulate infinite layers. '''
camfr.set_right_PML(-0.04)
camfr.set_left_PML(-0.04)
camfr.set_upper_PML(-0.04)
camfr.set_lower_PML(-0.04)


s.calc()    # calculate the modes.

''' Plot some modes. '''
#fig = s.plot(dx=0.3, dy=0.3)
fig = s.plot(field=['P','Ex','Ey'], mode=[0,1,2], dx=0.05, dy=0.05)

''' Save the figure as an image. '''
#fig.savefig('Silicon_WG_-_Modesolver_example_v1.png')


print( 'done.' )



