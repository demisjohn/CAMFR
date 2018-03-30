#! /usr/bin/env python

##############################################################################
#
# File:     camfr_PIL.py
# Authors:  Lieven.Vanholme@intec.ugent.be
#           Peter.Bienstman@ugent.be
#           Pieter.Dumon@intec.ugent.be
#
##############################################################################

from camfr import *
#from numpy import *
import numpy as np

# Colormap codes.

whiteblack  = 1
blackwhite  = 2
palet       = 3

# Colormap colors.

NEGCOLOR    = (0,0,255)         # Blue.
MIDCOLOR    = (255,255,255)     # White.
POSCOLOR    = (255,0,0)         # Red.

MAPCOLORS   = [(0  ,0  ,0  ),   # Black.
               (0  ,0  ,255),   # Blue.
               (255,0  ,255),   # Purple.
               (255,0  ,0  ),   # Red.
               (255,255,0  ),   # Yellow.
               (255,255,255)]   # White.

# PIL file formats.

formats = ['.gif', '.GIF', '.jpg', '.JPG', '.jpeg', '.JPEG', '.bmp', '.BMP',  
           '.dib', '.DIB', '.png', '.PNG', '.tiff', '.TIFF', '.tif', '.TIF',
           '.eps', '.EPS', '.ps', '.PS', '.pdf', '.PDF', '.xbm', '.XBM']

# Arrow vector variables.

ARROWSIZE = 10      # For each 10 data pixel, we get one arrow.

def set_arrowsize(arrsz=10):
   global ARROWSIZE
   ARROWSIZE = arrsz

#                   p3
#                    | \
#    p1-------------p2  \  
#     |     *           p4
#    p7-------------p6  /
#                    | /  
#                   p5
# poolcoordinates

p1 = [ 2.897, 0.344]
p2 = [ 0.245, 0.344]
p3 = [ 0.464, 0.377]
p4 = [ 0.000, 0.667]
p5 = [-0.464, 0.377]
p6 = [-0.245, 0.344]
p7 = [ 3.387, 0.344]    

ARROW = np.array([ p1, p2, p3, p4, p5, p6, p7 ])



##############################################################################
#
# Auxialiary draw function.
#
##############################################################################

def _create_window_and_draw(drawobject):
    
    from TkPlotCanvas import Frame, PlotCanvas, TOP, SUNKEN, BOTH, YES
    
    window = Frame()
    window.pack(fill=BOTH, expand=YES)

    def display(value):
        print(value)

    c = PlotCanvas( window, 500, 500, zoom=1, select=display,
                    relief=SUNKEN,border=2)
    c.pack(side=TOP, fill=BOTH, expand=YES)
    c.draw(drawobject, 'automatic', 'automatic')

    window.mainloop()


    
##############################################################################
#
# scatter_plot
#
##############################################################################

def scatter_plot(x, y):
    
    import TkPlotCanvas
    
    v = []
    for i in range(len(x)):
        v.append([x[i],y[i]])
        
    _create_window_and_draw(TkPlotCanvas.PolyMarker(v, fillcolor='white'))



##############################################################################
#
# plot_vector
#
##############################################################################

def plot_vector(v):

    import TkPlotCanvas
    
    pass
    try:
        is2d = v[0][0]
        _create_window_and_draw(TkPlotCanvas.PolyLine(v))
    except TypeError:
        w = [] 
        for i in range(len(v)):
            w.append([i,v[i]])
        _create_window_and_draw(TkPlotCanvas.PolyLine(w))



##############################################################################
#
# _create_scaled_matrix_plot
#
#  Low-level routine.
#  Given a colormap with N colors, returns a plot picture of matrix.
#  The matrix should contain elements in the range 0..N.
#  r_x is the range for the horizontal axis (matrix columns),
#  r_y for the vertical axis (rows).
#
##############################################################################

def _create_scaled_matrix_plot(colormap, z, r_x=0, r_y=0,
                               min_area = 100000, scale =1):
    import Image
    
    def round(x):
        return int(np.math.floor(x+.5))

    # Determine width and height of a pixel.

    [height, width, scale_x, scale_y] = \
    _scale_function( z, r_x, r_y, min_area, scale)

    # Prepare picture.
    
    pic = Image.new("RGB",(width*scale_x,height*scale_y))
    for x in range(width):
        for y in range(height):
            pic.paste(colormap[round(z[y,x])],
                      (scale_x* x,   scale_y* y,
                       scale_x*(x+1),scale_y*(y+1)))

    return pic


            
##############################################################################
#
#  _create_scaled_arrow_plot
#
#  Low-level routine.
#  Scales the plot and draws black arrows on a white canvas.
#
##############################################################################

def _create_scaled_arrow_plot(px, pz, r_x=0, r_y=0,
                              min_area = 100000, scale =1):
    import Image, ImageDraw
    global ARROWSIZE
    
    # Determine width and height of a vector.

    [height, width, scale_x, scale_y] = \
    _scale_function( pz, r_x, r_y, min_area, scale)

    # Prepare picture.
    
    len_x   = width*scale_x
    len_y   = height*scale_y

    pic     = Image.new("RGB",(len_x,len_y))  
    # This one is black, color is only defined in the manual.
    pic.paste( 0xFFFFFF,(0, 0,len_x,len_y))
    draw    = ImageDraw.Draw(pic)

    # Ofsety = len(r_y) - (height-1)*ARROWSIZE.
    pz*=scale_x
    px*=scale_y
    
    for x in np.array(range(width))[::ARROWSIZE]:
        for y in np.array(range(height))[::ARROWSIZE]:
            X = scale_x*x
            Y = scale_y*y
            _create_arrow(draw,(X,Y), pz[y,x], px[y,x])
            
    del draw

    return pic



##############################################################################
#
#  _scale_function
#
#  Low-level routine.
#  Returns the height, width and scale. Scale gives the size of one point in
#  canvas-pixels.
#
##############################################################################

def _scale_function(z, r_x, r_y, min_area, scale):

    def round(x):
        return int(np.math.floor(x+.5))
        
    height = z.shape[0]
    width  = z.shape[1]

    if type(r_x)!=np.ndarray or np.asarray(r_x).shape[0]==1: r_x = np.range(width)
    if type(r_y)!=np.ndarray or np.asarray(r_y).shape[0]==1: r_y = np.range(height)

    if len(r_x)>1:  d_x = r_x[1] - r_x[0]
    else :          d_x = 1
    if len(r_y)>1:  d_y = r_y[1] - r_y[0]
    else:           d_y = 1

    if (height * width * d_x * d_y  <  min_area):
        scale = np.math.sqrt(min_area/(height*width*d_x*d_y))

    if d_x < d_y:
        scale_x = round(scale*d_x)
        scale_y = scale_x*round(d_y/d_x)
    else:
        scale_y = round(scale*d_y)
        scale_x = scale_y*round(d_x/d_y)

    if scale_x < 1 or scale_y < 1:
        scale_x=1
        scale_y=1

    return [height, width, scale_x, scale_y]



##############################################################################
#
#  _create_arrow
#
#  Draws an arrow on the canvas. Size is dx and dy.
#  Center of the arrow is point p.
#
##############################################################################


def _create_arrow(draw, p, dx, dy):

    dq = np.sqrt(dx**2+dy**2)            # Distance.

    if (np.abs(dx) <= 1e-10): dx = 1e-10  # Don't set arc yet (+/-).
    arc = np.arctan(dy/dx)
    if dx < 0: arc += np.pi

    arrow = np.array(ARROW)

    # Turn the arrow.
    arrow[:,0] += arc
    # Size the arrow.
    arrow[:,1] *= dq

    points = []
    # Calc every point.
    for pt in arrow:
        dpx = pt[1] * np.cos(pt[0])
        dpy = pt[1] * np.sin(pt[0])
        points.append((p[0]+dpx,p[1]-dpy))

    draw.polygon(points, fill=0x000000)



##############################################################################
#
# _create_matrix_plot
#
#  Low-level routine. Scales matrix and sets colormap.
#
##############################################################################

def _create_matrix_plot(z_, r_x=0, r_y=0, colorcode=0,
                        min_area=100000, scale=1):

    import numpy.oldnumeric.mlab as ml
        
    # Scale z and find appropriate colormap.

    z = z_.copy()
    
    zmax = ml.max(ml.max(z))
    zmin = ml.min(ml.min(z))

    if (zmin < 0) and (0 < zmax) and (not colorcode):
        colormap = create_bipolar_colormap()
        zmax = ml.max(np.asarray([-zmin, zmax]))
        z += zmax
        z *= (len(colormap)-1)/(2*zmax)
    else:
        if colorcode == whiteblack:
            colormap = create_white_black_colormap()
        elif colorcode == blackwhite:
            colormap = create_black_white_colormap()    
        else:
            colormap = create_unipolar_colormap()

        z -= zmin
        if (zmax != zmin):
            z *= (len(colormap)-1)/(zmax-zmin)

    return _create_scaled_matrix_plot(colormap, z, r_x, r_y,
                                      min_area = min_area, scale = scale)



##############################################################################
#
#  _create_arrow_plot
#
#  Giving the arrows the right size.
#
##############################################################################

def _create_arrow_plot(px, pz, r_x=0, r_y=0,
                       min_area=100000, scale=1):
    import numpy.oldnumeric.mlab as ml
    
    # Scale pz & px
    
    pmax    = np.max( ml.max(ml.max(np.abs(pz))), ml.max(ml.max(np.abs(px))))
    if (pmax == 0): cst     = 0
    else:           cst     = ARROWSIZE/pmax
    pz      = pz*cst 
    px      = px*cst     

    return _create_scaled_arrow_plot(px, pz, r_x, r_y,
                                      min_area = min_area, scale = scale)



##############################################################################
#
# _output_pic
#
#  Low_level routine. Writes picture to a file or shows it on canvas.
#
##############################################################################

def _output_pic(pic, filename=0):

    import Tkinter, ImageTk, os, sys
    
    if filename:
        if '.' in filename:
            (name, suffix) = os.path.splitext(filename)
        else:
            suffix  = '.jpg'
            filename += suffix
        slash       = '/'
        script      = sys.argv[0]
        script_path = os.path.abspath(script)
        pic_path    = os.path.dirname(script_path)+slash+filename
        if not suffix in formats:
            print("File format not supported. Defaulting to jpg.")
            pic_path += ".jpg"

        pic.save(pic_path)    
        print("Created", pic_path)
    else:
        root     = Tkinter.Tk()    
        canvas   = Tkinter.Canvas(width=pic.size[0], height=pic.size[1])
        backdrop = ImageTk.PhotoImage(pic)
        canvas.create_image(0, 0, image=backdrop, anchor=Tkinter.NW)
        canvas.pack()

        root.mainloop()



##############################################################################
#
# _overlay_pictures
#
#  Low_level routine. If contour=1, overlays the contours of pic2 with
#  pic1. Otherwise the entire pic2 is overlayed with pic1.
#
##############################################################################

def _overlay_pictures(pic1, pic2, contour):

    if (contour):
        import ImageFilter, ImageChops
        return ImageChops.multiply(pic2.filter(ImageFilter.CONTOUR), pic1)
    else:
        import Image
        return Image.blend(pic2, pic1, 0.5)



##############################################################################
#
# Plot a matrix.
#
#  High-level routine.
#
##############################################################################

def plot_matrix(z, r_x=0, r_y=0, filename=0, colorcode=0):

    _output_pic(_create_matrix_plot(z, r_x, r_y, colorcode), filename)



##############################################################################
#
# Plot an arrow matrix.
#
#  High-level routine.
#
##############################################################################

def plot_arrow(px, pz, r_x=0, r_z=0, filename=0):

    _output_pic(_create_arrow_plot(px, pz, r_z, r_x), filename)



##############################################################################
#
# Create a phasor movie.
#
##############################################################################

def _create_phasor_movie(z_, r_x=0, r_y=0, min_area=100000, scale=1, ln=0):
    
    import numpy.oldnumeric.mlab as ml
    
    # Make movie memory.
    
    movie    = []
    frames   = 16

    z = z_.copy()

    if ln:
        colormap = create_unipolar_colormap()
        
        # Scale factors for z.
        # We're not sure about the lowest zmin value during the calculations.
        # So, we make it our own by adding a 'zcst' to each one.
        # Log(z**2) gives a (much) nicer plot than log(abs(z)).
        # We should get positive values anyhow.

        # Supplement to make sure the z isn't zero this really changes
        # the feeling of the image.. so not that correct.
        zcst    = 1e-8

        zmax    = 2 * np.log(ml.max(ml.max(np.abs(z)))+ zcst)
        zmin    = np.log(zcst)                 # Around -23.
        z_scale = (len(colormap)-1)/(zmax-zmin)
        
        # Calculate each frame.
        
        for Nr in np.range(0,frames):
            pic = _create_scaled_matrix_plot(colormap,
                        ((np.log(z.real**2 + zcst) - zmin)*z_scale),
                        r_x, r_y, min_area = min_area, scale = scale)
            movie.append(pic)
            z *= np.exp(2j*pi/frames)
            
    else:
        colormap = create_bipolar_colormap()
        
        # Scale factors for z.

        zmax    = ml.max(ml.max(np.abs(z)))
        if (zmax == 0):
            # in this case, z=0, the middle of the color palet
            #(z+zmax)*z_scale) should be = len(colormap)/2
            z_scale = (len(colormap)-1)/2
            zmax = 1 
        else:
            z_scale = (len(colormap)-1)/(2*zmax)

        # Calculate each frame.
        
        for Nr in np.range(0,frames):
            pic = \
              _create_scaled_matrix_plot(
                  colormap,((z+zmax)*z_scale).real,r_x,r_y,
                            min_area = min_area, scale = scale)
            movie.append(pic)
            z *= np.exp(2j*np.pi/frames)

    return movie



##############################################################################
#
# _output_movie
#
#  Low_level routine. Writes movie to a file or shows it on canvas.
#
##############################################################################

def _output_movie(movie, filename):

    import Tkinter, ImageTk, gifmaker, os, sys

    frames = len(movie)

    if filename:
        #for Nr in range(frames):
        #    fn = filename + str(Nr) + ".jpg"
        #    movie[Nr].save(fn)
        #return
        
        for Nr in range(0,frames):
            movie[Nr] = movie[Nr].convert("P")  
        if '.' in filename:
            (name, suffix) = os.path.splitext(filename)
            if suffix != '.gif':
                print("File format not supported. Defaulting to gif.")
                filename = filename[:filename.index('.')] + ".gif"
        else:
            filename += ".gif"
            
        slash     = '/'
        script    = sys.argv[0]
        totalpath = os.path.abspath(script)
        userpath  = os.path.dirname(totalpath)+slash+filename
        fp = open(userpath,"wb")
        
        gifmaker.makedelta(fp, movie)
        fp.close()
        print( "Created", userpath )

    else:
        root = Tkinter.Tk()
        
        # Close window procedure.

        stop = [0]    
        def callback(): stop[0] = 1 
        root.protocol("WM_DELETE_WINDOW", callback)
        
        # Animate picture.
        
        canvas = Tkinter.Canvas(width=movie[0].size[0],
                                height=movie[0].size[1])

        while not stop[0]:
            for x in np.range(frames):
                if stop[0]: break
                backdrop = ImageTk.PhotoImage(movie[x])
                canvas.create_image(0, 0, image=backdrop, anchor=Tkinter.NW)
                canvas.pack()
                root.update()
                root.after(int(100))

        root.destroy()       



##############################################################################
#
# Creates a movie starting from a complex matrix representing phasors.
#
##############################################################################

def phasormovie(z, r_x=0, r_y=0, filename=0, ln=0):

    _output_movie(_create_phasor_movie(z, r_x, r_y, ln=ln), filename)



##############################################################################
#
# Different colormaps.
#
# _create_color_range creates a range of colors, making a
# smooth trannsition between two colors.
#
##############################################################################

def create_bipolar_colormap():

    return _create_color_range(NEGCOLOR,MIDCOLOR ) +\
           _create_color_range(MIDCOLOR,POSCOLOR, 1 )


def create_unipolar_colormap():

    colormap, l, end = [], len(MAPCOLORS), 0

    if (l > 1) :
        for x in np.range(l-1):
            if (x == l-2): end = 1            
            colormap += _create_color_range(MAPCOLORS[x],MAPCOLORS[x+1], end)
            
    return colormap


def create_white_black_colormap():
    return  _create_color_range((255,255,255), (0,0,0), 1 )


def create_black_white_colormap():
    return  _create_color_range((0,0,0), (255,255,255), 1 )

def _create_color_range(c1=(0,0,0), c2=(255,255,255), ad_last = 0 ):

    # Calc max range.
    diff    = np.array(c2)-array(c1)  
    colors  = float(np.max(np.abs(diff)))

    # Dred,dgreen,dblue.
    dc = diff / colors

    # Start_to_end + dc.
    if ad_last: colors +=1

    return [sum( (c1 + (x*dc).astype(int))*[1, 0x100, 0x10000] )
            for x in np.range(int(colors))]



##############################################################################
#
# Backend-independent visualisation functions.
#
#  They are duplicated in each backend however, to allow flexible run time
#  switching of backends.
#
##############################################################################

##############################################################################
#
# Plot distribution of effective indices in complex plane.
#
##############################################################################
    
def plot_neff(waveguide):
    
    x,y = [],[]
    
    for i in np.range(waveguide.N()):
        n = waveguide.mode(i).n_eff()
        x.append(n.real)
        y.append(n.imag)

    scatter_plot(x,y)



##############################################################################
#
# Plot a complex function.
#
##############################################################################

def plot_f(f, r_x, r_y, filename=0, colormap=palet):
    
    fz = np.zeros([len(r_y),len(r_x)], float)

    for i_y in np.range(len(r_y)):
      for i_x in np.range(len(r_x)):
        fz[len(r_y)-1-i_y, i_x] = np.abs(f(r_x[i_x] + r_y[i_y]*1j))

    plot_matrix(fz, r_x, r_y, filename, colormap)



##############################################################################
#
# Plot the refractive index profile in a waveguide.
#
##############################################################################

def plot_n_waveguide(waveguide, r_x):
    
    v = []
    
    for i_x in np.range(len(r_x)):
        v.append((r_x[i_x], np.abs(waveguide.n(Coord(r_x[i_x],0,0)))))
        
    plot_vector(v)



##############################################################################
#
# Plot the refractive index profile in a stack.
#
##############################################################################

def plot_n_stack(stack, r_x, r_z, r_y=0, filename=0, colormap=whiteblack):
   
    rxrange = False
    ryrange = False
    rzrange = False
    try:
      rx = float(r_x)
      ax1 = r_y
      ax2 = r_z
      r_x = [r_x]
      axc = r_x
    except TypeError:
      rxrange = True
    try:
      ry = float(r_y)
      ax1 = r_x
      ax2 = r_z
      r_y = [r_y]
      axc = r_y
    except TypeError:
      ryrange = True
    try:
      rz = float(r_z)
      ax1 = r_y
      ax2 = r_x
      r_z = [r_z]
      axc = r_z
    except TypeError:
      rzrange = True
    if rxrange and ryrange and rzrange:
      print("Error: plot_n_stack can only make cross sections.")

    n = np.zeros([len(ax1),len(ax2)],float)
    if rzrange:
     _calc_n_stack(n, stack, r_x[::-1], r_y[::-1], r_z)
    else:
     _calc_n_stack(n, stack, r_x, r_y[::-1], r_z)
     
    plot_matrix(n,ax2,ax1,filename,colormap)


##############################################################################
#
# Plot the poynting arrows in a stack.
#
##############################################################################

def plot_arrow_stack(stack, r_x, r_z, r_y = 0, filename=0):
   
    rxrange = False
    ryrange = False
    rzrange = False
    try:
      rx = float(r_x)
      r_x = [r_x]
    except TypeError:
      rxrange = True
    try:
      ry = float(r_y)
      r_y = [r_y]
    except TypeError:
      ryrange = True
    try:
      rz = float(r_z)
      r_z = [r_z]
    except TypeError:
      rzrange = True
    if rxrange and ryrange and rzrange:
      print("Error: plot_n_stack can only make cross sections")
    px   = np.zeros([len(r_x),len(r_z)], float)
    py   = np.zeros([len(r_y),len(r_z)], float)
    pz   = np.zeros([len(r_x),len(r_z)], float)
    if rzrange:
      _calc_arrow_stack(px, py, pz, stack, r_x[::-1], r_y[::-1], r_z)
    else:
      _calc_arrow_stack(px, py, pz, stack, r_x, r_y[::-1], r_z)
    if not rxrange:
      plot_arrow(py, pz, r_y, r_z, filename)
    elif not ryrange:
      plot_arrow(px, pz, r_x, r_z, filename)
    else:
      plot_arrow(px, py, r_x, r_y, filename)



##############################################################################
#
# Plot the refractive index profile in a Section.
#
##############################################################################

def plot_n_section(stack, r_x, r_y, filename, colormap):
    
    n = np.zeros([len(r_y),len(r_x)], float)
    _calc_n_stack(n, stack, array(r_x)[::-1], r_y)

    plot_matrix(n, r_x, r_y, filename, colormap)    



##############################################################################
#
# Wrapper for plot_n.
#
##############################################################################

def plot_n(o, r1, r2=0, r3=0, filename=0, colormap=whiteblack):

    if type(r2)!=np.ndarray or np.asarray(r2).shape[0]==1:
        plot_n_waveguide(o, r1)
    elif type(o) == Stack or type(o) == BlochStack or type(o) == Cavity:
        if type(r3)!=np.ndarray or np.asarray(r3).shape[0]==1:
          plot_n_stack(o, r1, r2, 0.0, filename, colormap)
        else:
          plot_n_stack(o, r1, r3, r2, filename, colormap)
    elif type(o) == Section:
        plot_n_section(o, r1, r2, filename, colormap)
    else:
        print("Unsupported argument for plot_n.")



##############################################################################
#
# Plot the field profile of a waveguide mode.
#
##############################################################################

def plot_field_waveguide(mode, component, r_x):
    
    v = []
    
    for i_x in range(len(r_x)):
      v.append((r_x[i_x],component(mode.field(Coord(r_x[i_x],0,0)))))
        
    plot_vector(v)



##############################################################################
#
# Plot the field profile in a stack.
#
##############################################################################

def plot_field_stack(stack, component, r_x, r_z, r_y = 0, filename=0,
                     colormap=whiteblack, overlay_n=1, contour=1, arrow=0):

    rxrange = False
    ryrange = False
    rzrange = False
    try:
      rx = float(r_x)
      ax1 = r_y
      ax2 = r_z
      r_x = [r_x]
      axc = r_x
    except TypeError:
      rxrange = True
    try:
      ry = float(r_y)
      ax1 = r_x
      ax2 = r_z
      r_y = [r_y]
      axc = r_y
    except TypeError:
      ryrange = True
    try:
      rz = float(r_z)
      ax1 = r_y
      ax2 = r_x
      r_z = [r_z]
      axc = r_z
    except TypeError:
      rzrange = True
    if rxrange and ryrange and rzrange:
      print("Error: plot_n_stack can only make cross sections")

    f = np.zeros([len(ax1),len(ax2)], float)
    if rzrange:
      _calc_field_stack(f, stack, r_x[::-1], r_y[::-1], component, r_z)
    else:
      _calc_field_stack(f, stack, r_x, r_y[::-1], component, r_z)
    if not (overlay_n or arrow) :
        return plot_matrix(f, ax2, ax1, filename, colormap)

    # Overlay index/arrow profile.
    
    if overlay_n:
        n = np.zeros([len(ax1),len(ax2)], float)
        if rzrange:
          _calc_n_stack(n, stack, r_x[::-1], r_y[::-1], r_z)
        else:
          _calc_n_stack(n, stack, r_x, r_y[::-1], r_z)
        pic_n = _create_matrix_plot(n, ax2, ax1, whiteblack)

    if arrow:
        px   = np.zeros([len(r_x),len(r_z)], float)
        py   = np.zeros([len(r_y),len(r_z)], float)
        pz   = np.zeros([len(r_x),len(r_z)], float)
        if rzrange:
          _calc_arrow_stack(px, py, pz, stack, r_x[::-1], r_y[::-1], r_z)
        else:
          _calc_arrow_stack(px, py, pz, stack, r_x, r_y[::-1], r_z)
        if not rxrange:
          pic_p = _create_arrow_plot(py, pz, r_z, r_y)
        elif not ryrange:
          pic_p = _create_arrow_plot(px, pz, r_z, r_x)
        else:
          pic_p = _create_arrow_plot(px, py, r_x, r_y)

    if arrow:
        if overlay_n:
            pic_n = _overlay_pictures(pic_p, pic_n, 0)
        else:
            pic_n = pic_p

    pic_f = _create_matrix_plot(f, ax2, ax1, colormap)

    _output_pic(_overlay_pictures(pic_f, pic_n, contour), filename)



##############################################################################
#
# Calculates the field values for all the points of the matrix.
#
##############################################################################

def _calc_field_stack(f, stack, r_x, r_y, component, r_z=0):

   if type(r_z)!=np.ndarray or asarray(r_z).shape[0]==1:
        # 2D 
    for x in range(len(r_x)):
        for z in range(len(r_y)):
            f[x,z] = component(stack.field(Coord(r_x[x],0,r_y[z])))
            
   elif len(r_x)==1:
    for y in range(len(r_y)):
        for z in range(len(r_z)):
            f[y,z] = component(stack.field(Coord(r_x[0], r_y[y], r_z[z])))
   elif len(r_y)==1:
    for x in range(len(r_x)):
        for z in range(len(r_z)):
            f[x,z] = component(stack.field(Coord(r_x[x], r_y[0], r_z[z])))
   else:
    for x in range(len(r_x)):
        for y in range(len(r_y)):
            f[x,y] = component(stack.field(Coord(r_x[x], r_y[y], r_z[0])))



##############################################################################
#
# Asks the refractive index for all the points of the matrix.
#
##############################################################################

def _calc_n_stack(n, stack, r_x, r_y, r_z=0):

    if type(r_z)!=ndarray or asarray(r_z).shape[0]==1:
        # 2D
        for x in range(len(r_x)):
            for z in range(len(r_y)):
                n[x,z] = stack.n(Coord(r_x[x],0,r_y[z])).real
            
    elif len(r_x)==1:
      for y in range(len(r_y)):
         for z in range(len(r_z)):
            n[y,z] = stack.n(Coord(r_x[0], r_y[y], r_z[z])).real

    elif len(r_y)==1:
      for x in range(len(r_x)):
         for z in range(len(r_z)):
            n[x,z] = stack.n(Coord(r_x[x], r_y[0], r_z[z])).real

    else:
      for x in range(len(r_x)):
         for y in range(len(r_y)):
            n[y,x] = stack.n(Coord(r_x[x], r_y[y], r_z[0])).real



##############################################################################
#
# Calculates the poynting values for some points of the matrix.
# The number of points depends on the size of the arrow.
#
##############################################################################

def _calc_arrow_stack(px, pz, stack, r_x, r_z):
    
    global ARROWSIZE
    if get_polarisation() == TM:
        for x in np.array(range(len(r_x)))[::ARROWSIZE]:
            for z in np.array(range(len(r_z)))[::ARROWSIZE]:
                f       = stack.field(Coord(r_x[x], 0, r_z[z]))
                h2      = f.H2().conjugate()
                px[x,z] = (-f.Ez() * h2).real
                pz[x,z] = ( f.E1() * h2).real

    elif get_polarisation() == TE:
        for x in np.array(range(len(r_x)))[::ARROWSIZE]:
            for z in np.array(range(len(r_z)))[::ARROWSIZE]:
                f       = stack.field(Coord(r_x[x], 0, r_z[z]))
                e2      = f.E2()
                px[x,z] = ( e2 * f.Hz().conjugate()).real
                pz[x,z] = (-e2 * f.H1().conjugate()).real

    elif get_polarisation() == TE_TM: # 3D stack
        if len(r_x)==1:
          for y in np.array(range(len(r_y)))[::ARROWSIZE]:
             for z in np.array(range(len(r_z)))[::ARROWSIZE]:
                f = stack.field(Coord(x,r_y[y], r_z[z]))
                py[y,z] = (f.H1().conjugate()*f.Ez() \
                           -f.E1()*f.Hz().conjugate()).real
                pz[y,z] = (f.E1()*f.H2().conjugate() \
                           -f.H1().conjugate()*f.E2()).real
        elif len(r_y)==1:
          for x in np.array(range(len(r_x)))[::ARROWSIZE]:
             for z in np.array(range(len(r_z)))[::ARROWSIZE]:
                f = stack.field(Coord(r_x[x],y, r_z[z]))
                px[x,z] = (f.E2()*f.Hz().conjugate() \
                           -f.H2().conjugate()*f.Ez()).real
                pz[x,z] = (f.E1()*f.H2().conjugate() \
                           -f.H1().conjugate()*f.E2()).real
        else:
          for x in np.array(range(len(r_x)))[::ARROWSIZE]:
             for y in np.array(range(len(r_y)))[::ARROWSIZE]:
                f = stack.field(Coord(r_x[x],r_y[y], z))
                px[y,x] = (f.E2()*f.Hz().conjugate() \
                           -f.H2().conjugate()*f.Ez()).real
                py[y,x] = (f.E1()*f.H2().conjugate() \
                           -f.H1().conjugate()*f.E2()).real

    else:
        print("Error: no polarisation defined.")

    

##############################################################################
#
# Plot the field profile of a section mode.
#
##############################################################################

def plot_field_section_mode(mode, component, r_x, r_y, filename, colormap,
                            overlay_n=1, contour=1):
    
    f = np.zeros([len(r_y),len(r_x)], float)
    
    for i_x in np.range(len(r_x)):
      for i_y in np.range(len(r_y)):
        f[len(r_y)-1-i_y,i_x] = \
              component(mode.field(Coord(r_x[i_x], r_y[i_y], 0)))

    if not overlay_n:
        return plot_matrix(f, r_x, r_y, filename, colormap)

    # Overlay index profile.
        
    n = np.zeros([len(r_y),len(r_x)], float)
        
    for i_x in np.range(len(r_x)):
        for i_y in np.range(len(r_y)):
            n[len(r_y)-1-i_y,i_x] = mode.n(Coord(r_x[i_x], r_y[i_y], 0)).real

    pic_n = _create_matrix_plot(n, r_x, r_y, whiteblack)
    pic_f = _create_matrix_plot(f, r_x, r_y, colormap)

    _output_pic(_overlay_pictures(pic_f, pic_n, contour), filename)



##############################################################################
#
# Wrapper for plot_field.
#
##############################################################################

def plot_field(o, component, r1, r2=0, r3=0, filename=0,
               colormap=0, overlay_n=1, contour=1, arrow=0):

    if type(r2)!=np.ndarray or np.asarray(r2).shape[0]==1:
        plot_field_waveguide(o, component, r1)
    elif type(o) == Stack or type(o) == BlochMode or type(o) == Cavity:
        if type(r3)!=np.ndarray or np.asarray(r3).shape[0]==1: # 2D stack.
          plot_field_stack(o, component, r1, r2, 0, filename, colormap,
                         overlay_n,contour,arrow)
        else:
          plot_field_stack(o, component, r1, r3, r2, filename, colormap,
                         overlay_n,contour,arrow)

    elif type(o) == SectionMode:
        plot_field_section_mode(o,component,r1,r2,filename,colormap,overlay_n,
                                contour)
    else:
        print("Unsupported argument for plot_field.")
        


##############################################################################
#
# Animate the field profile in a stack.
#
##############################################################################

def animate_field_stack(stack, component, r_x, r_z, r_y = 0, filename=0,
                        overlay_n=1, contour=1, ln=0):

    rxrange = False
    ryrange = False
    rzrange = False
    try:
      rx = float(r_x)
      ax1 = r_y
      ax2 = r_z
      r_x = [r_x]
      axc = r_x
    except TypeError:
      rxrange = True
    try:
      ry = float(r_y)
      ax1 = r_x
      ax2 = r_z
      r_y = [r_y]
      axc = r_y
    except TypeError:
      ryrange = True
    try:
      rz = float(r_z)
      ax1 = r_y
      ax2 = r_x
      r_z = [r_z]
      axc = r_z
    except TypeError:
      rzrange = True
    if rxrange and ryrange and rzrange:
      print("Error: plot_n_stack can only make cross sections")

    f = np.zeros([len(ax1),len(ax2)], complex)
    if rzrange:
      _calc_field_stack(f, stack, r_x[::-1], r_y[::-1], component, r_z)
    else:
      _calc_field_stack(f, stack, r_x, r_y[::-1], component, r_z)
    if not overlay_n:
        return phasormovie(f, ax2, ax1, filename, ln)

    # Overlay index profile.

    n = np.zeros([len(ax1),len(ax2)], float)
    if rzrange:
      _calc_n_stack(n, stack, r_x[::-1], r_y[::-1], r_z)
    else:
      _calc_n_stack(n, stack, r_x, r_y[::-1], r_z)

    pic_n = _create_matrix_plot (n, ax2, ax1, whiteblack)
    mov_f = _create_phasor_movie(f, ax2, ax1, ln=ln)

    for Nr in np.range(len(mov_f)):
        mov_f[Nr] = _overlay_pictures(mov_f[Nr], pic_n, contour)

    _output_movie(mov_f, filename)



##############################################################################
#
# Animate the field profile of a section mode.
#
##############################################################################

def animate_field_section_mode(mode, component, r_x, r_y, filename=0,
                               overlay_n=1, contour=1):
    
    f = np.zeros([len(r_y),len(r_x)], complex)

    for i_x in np.range(len(r_x)):
      for i_y in np.range(len(r_y)):
        f[len(r_y)-1-i_y,i_x] = \
              component(mode.field(Coord(r_x[i_x], r_y[i_y], 0)))

    if not overlay_n:
        return phasormovie(f, r_x, r_y, filename)

    # Overlay index profile.
        
    n = np.zeros([len(r_y),len(r_x)], float)
        
    for i_x in np.range(len(r_x)):
        for i_y in np.range(len(r_y)):
            n[len(r_y)-1-i_y,i_x] = mode.n(Coord(r_x[i_x], r_y[i_y], 0)).real

    pic_n = _create_matrix_plot (n, r_x, r_y, whiteblack)
    mov_f = _create_phasor_movie(f, r_x, r_y)

    for Nr in np.range(len(mov_f)):
        mov_f[Nr] = _overlay_pictures(mov_f[Nr], pic_n, contour)

    _output_movie(mov_f, filename)   



##############################################################################
#
# Wrapper for animate_field.
#
##############################################################################

def animate_field(o, component, r1, r2, r3=0, filename=0, overlay_n=1,
                  contour=1, ln=0):

    if type(o) == Stack or type(o) == BlochMode or type(o) == Cavity:
      if type(r3)!=np.ndarray: # 2D
         animate_field_stack(o,component,r1,r2,0.0,filename,
                            overlay_n,contour,ln)
      else:
         animate_field_stack(o,component,r1,r3,r2,filename,
                            overlay_n,contour,ln)
    elif type(o) == SectionMode:
        animate_field_section_mode(o,component,r1,r2,filename,
                                   overlay_n,contour)
    else:
        print("Unsupported argument for animate_field.")



##############################################################################
#
# Inject plot functions into C++ classes.
#
##############################################################################

import slab_plot, stack_plot

Slab.plot       = lambda self : slab_plot.SlabPlot(self)
Stack.plot      = lambda self : stack_plot.StackPlot(self)
BlochStack.plot = lambda self : stack_plot.StackPlot(self)
Cavity.plot     = lambda self : stack_plot.StackPlot(self)

Slab.plot_n = plot_n
Circ.plot_n = plot_n

Mode.plot_field = plot_field

Stack.plot_n        = plot_n
Stack.plot_field    = plot_field
Stack.animate_field = animate_field

BlochStack.plot_n = plot_n

BlochMode.plot_field    = plot_field
BlochMode.animate_field = animate_field

Cavity.plot_n        = plot_n
Cavity.plot_field    = plot_field
Cavity.animate_field = animate_field

