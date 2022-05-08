#! /usr/bin/env python

###########################################################################
#
# File:     matplotlib_section.py
# Author:   demis@ucsb.edu
# Date:     20180311
# Version:  1.0
#
# Copyright (C) 2018 Demis D. John - Univ. of California Santa Barbara
#
############################################################################

from _camfr import *    # import the Section and Slab classes, in order to add functions to them.
import numpy as np
import matplotlib.pyplot as plt



############################################################################
#
# This file contains MatPlotLib-based plotting functions that will be added to
# various CAMFR classes, such as the Section and Slab classes.
# To add functions to classes defined externally, first create a function with 
# a placeholder name, typically beginning with two underscores __, as this 
# indicates that a user should generally not be calling this func directly.
# Then, add an attribute to the target class that points to this function.
# `self` will be the object, once this func is called as a method of that object
# 
# Demis D. John, 2018-03-11, Univ. of California Santa Barbara
#
############################################################################

############################################################################
#
# Choose some global plotting options:
#
############################################################################

## Colormap
colormap = plt.get_cmap('hot')

AxisBGColor = 'black'   # background color of every axis



############################################################################
#
# Functions
#
############################################################################

def __Section_plot(self, field="Ex", mode=0, dx=0.100, dy=0.100, annotations=True):
    '''
    Plot a 2D mode profile of the specified field, using MatPlotLib.
    
    Parameters
    ----------
    ModeObj: a CAMFR Mode object, often acquired via `SectionObj.mode(0)`.
    
    field : {'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'}
        Electric (E) or Magnetic (H) field, x/y/z component.  
        Case-insensitive.
        Defaults to 'Ex', x-component (horizontal) of the electric field, or `E1` in CAMFR parlance.
        Can pass an iterable of strings, such as 
            ['Ex', 'Ey']
        to plot multiple fields on a single figure.  Each field will be plotted along the rows of figure axes.
    
    mode : integer
        Which waveguide mode to plot, an integer from 0 to get_N(), however many modes are calculated for the Section
        Defaults to 0.
        Can pass an iterable of integers, such as
            [0, 1, 2]
        to plot multiple modes on a single figure.  Each mode will be plotted along the columns of figure axes.
    
    dx, dy: float
        x & y resolution for Mode Profiles & Gamma/field plotting. Default is 0.100um for both.
    
    annotations : boolean, optional
        If true, the effective index, mode number and field component will written on each mode plot.  True by default.
    
    Not Implemented:
        field = 'P' will plot normalized Power.

    
    
    Returns
    -------
    Returns the handle to the matplotlib figure object. This allows you to save or close the figure etc., via
        modeplot(mode=0).savefig('savedimage.png')
    or
        >>> fig = modeplot(mode=0, field='Ez')
        >>> fig.savefig('savedimage.png')
        >>> from matplotlib.pyplot import close as plot_close
        >>> plot_close(fig)     # close the figure during loops
    You can also retrieve the axis, pcolormesh, axes etc. objects from this handle as so:
        >>> [ax] = fig.get_axes()
        >>> QuadMeshObj = ax.get_children()[0]
        >>> XAxisObj = ax.get_children()[5]

    '''
    
    ## sanitize `field` argument, check if iterable
    if hasattr(field, '__iter__'):
        numfields = len(field) 
    else:
        numfields = 1;      # loop only once
        field = [field]     # make iterable with one item
        
    ## Convert `field` strings to corresponding CAMFR string
    #   See pg. 59 `Field` of CAMFR manual for available options.
    #   For a Section, axis 1 = X, horizontal, and axis 2 = Y, vertical.
    # user options :: CAMFR option as dictionary
    fieldopts = {}  # initialize dictionary
    fieldopts['ex'] = 'E1'
    fieldopts['ey'] = 'E2'
    fieldopts['ez'] = 'Ez'
    fieldopts['hx'] = 'H1'
    fieldopts['hy'] = 'H2'
    fieldopts['hz'] = 'Hz'
    fieldopts['p'] = 'abs_S'
    
    # create the complementary array of CAMFR fields
    try:
        cfield = [fieldopts[a] for a in    [b.lower() for b in field]  ]
    except KeyError as k:
        raise ValueError(  "Unrecognized value %s found for the `field` argument. " % k  + "The following options are valid: %s"  % fieldopts.keys()  )
    #end try(fieldopts)
    
        
    
    
    ## sanitize `mode` argument, check if iterable
    if hasattr(mode, '__iter__'):
        nummodes = len(mode)
    else:
        nummodes = 1;      # loop only once
        mode = [mode]      # make iterable with one item
    
    # get Section object
    obj = self
    
    
    # no harm in re-importing already imported modules, no extra processing time:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm      # colormaps
    import numpy as np
    plt.ion()   # interactive plotting mode, for live updating
    
    

    
    fig, ax = plt.subplots(nrows=nummodes, ncols=numfields, sharex=True, sharey=True, squeeze=False)

    
    #print(  "Mode Profile Resolution: (%f,%f)" %(dx,dy)  )
    Modes = np.copy(mode)   # convert to array
    #Fields = np.copy(field)      
    m = -1
    for modeN in Modes:
        f = -1
        m = m+1
        for camfrfield in cfield:
            f = f+1
            print( "Calculating fields for Mode %i: %s" %( modeN,field[f].title() )  )
            x = np.arange(0, obj.width(), dx)
            y = np.arange(0, obj.height(), dy)
            X,Y = np.meshgrid(x,y)
            
            
            intensity = np.zeros((np.size(y),np.size(x)), float)
            F = np.zeros((np.size(y),np.size(x)), float)

            
            for i in range(0, np.size(y)):
                for j in range(0, np.size(x)):
                    # Use getattr to allow user to specify the field to plot/use.
                    F[i,j] =  \
                        np.abs( 
                            getattr(    obj.mode(modeN).field(Coord(X[i,j], Y[i,j], 0)) , camfrfield)() 
                        ) ** 2.0
                #end for(x)
            #end for(y)
            
            #TEtot = np.sum(TEfield[:])
            #TMtot = np.sum(TMfield[:])
            #TEfrac = 100.*(TEtot/(TEtot+TMtot))
            
            #print( '* TE Fraction: %3.1f ' %TEfrac + 'n = %1.7f + i* ' % obj.mode(modeN).n_eff().real + '%1.7e' % obj.mode(modeN).n_eff().imag   )
            
            ## Plot the specified mode/field:
            axis = ax[ m , f ]  # which axis to plot on
            
            axis.pcolormesh( X, Y, F, cmap=colormap )
            if m==( len(Modes)-1 ):   axis.set_xlabel(r'x ($\mu{}m$)')  # LaTeX notation, overkill
            if f==0:    axis.set_ylabel(r'y ($\mu{}m$)')
            axis.set_xlim( axis.get_xlim()[0], obj.width() )
            axis.set_ylim( axis.get_ylim()[0], obj.height() )
            #axis.set_axis_bgcolor( AxisBGColor )   # this version works for matplotlib <v2.0
            axis.set_facecolor( AxisBGColor )       # only matplotlib >v2.0
            
            
            if annotations:
                titlestr = "Mode(" + str(modeN) + "): " + field[f].title()
                #axis.set_title(  titlestr  )
                axis.text( 0.05, 0.9, titlestr, transform=axis.transAxes, horizontalalignment='left', color='green', fontsize=9, fontweight='bold')
                
                n_str = "$\mathregular{n_{eff} =}$ %0.5f" % ( obj.mode(modeN).n_eff().real )
                if f==0: axis.text( 0.05, 0.05, n_str, transform=axis.transAxes, horizontalalignment='left', color='green', fontsize=9, fontweight='bold')
            #end if(annotations)
            
            ## update the plots:
            if (m==0) and (f==0): 
                fig.canvas.window().raise_()    # bring plot window to front (a hack - delete this if it causes trouble)
            fig.canvas.draw()   # update the figure
            plt.pause(0.05)     # allow GUI to update (may pop a warning)
            
        # end for Fields
    # end for Modes

    return fig
        
#end _Section_plot()


# add the above function to the Section class:
Section.plot   =   __Section_plot  
# This determines the real name of the function as a Section method, and points to this function.
