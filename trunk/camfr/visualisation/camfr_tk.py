from camfr import *
from Numeric import *

# Colormap codes.

whiteblack  = 1
blackwhite  = 2
palet       = 3

##############################################################################
#
# Auxialiary draw function.
#
##############################################################################

def create_window_and_draw(drawobject):
    
    from TkPlotCanvas import *
    
    window = Frame()
    window.pack(fill=BOTH, expand=YES)

    def display(value):
	print value

    c = PlotCanvas(window,500,500,zoom=1,select=display,relief=SUNKEN,border=2)
    c.pack(side=TOP, fill=BOTH, expand=YES)
    c.draw(drawobject, 'automatic', 'automatic')

    window.mainloop()


    
##############################################################################
#
# scatter_plot.
#
##############################################################################

def scatter_plot(x, y):
    
    import TkPlotCanvas
    
    v = []
    for i in range(len(x)):
        v.append([x[i],y[i]])
        
    create_window_and_draw(TkPlotCanvas.PolyMarker(v, fillcolor='white'))



##############################################################################
#
# Plot a vector.
#
##############################################################################

def plot_vector(v):

    import TkPlotCanvas
    
    pass
    try:
        is2d = v[0][0]
        create_window_and_draw(TkPlotCanvas.PolyLine(v))
    except TypeError:
        w = [] 
        for i in range(len(v)):
            w.append([i,v[i]])
        create_window_and_draw(TkPlotCanvas.PolyLine(w))



##############################################################################
#
# plot_scaled_matrix
#
#  Given a colormap with N colors, plots a matrix containing
#  elements in the range 0..N.
#  r_x is the range for the horizontal axis (matrix columns),
#  r_y for the vertical axis (rows)
#
##############################################################################

def plot_scaled_matrix(root, colormap, z, r_x=0, r_y=0):

    import Tkinter
    
    def round(x):
        return int(math.floor(x+.5))

    # Determine width and height of a pixel.

    height = z.shape[0]
    width  = z.shape[1]

    if not r_x: r_x = range(width)
    if not r_y: r_y = range(height)

    d_x = r_x[1] - r_x[0]
    d_y = r_y[1] - r_y[0]
    
    min_area, scale = 100000, 1
    if (height*width*d_x*d_y < min_area):
        scale = math.sqrt(min_area/(height*width*d_x*d_y))

    if d_x < d_y:
        scale_x = round(scale*d_x)
        scale_y = scale_x*round(d_y/d_x)
    else:
        scale_y = round(scale*d_y)
        scale_x = scale_y*round(d_x/d_y)
    
    # Prepare picture.

    pic = Tkinter.Canvas(root, width=scale_x*width,
                         height=scale_y*height, bg="White")
    
    for x in range(width):
        for y in range(height):
            pic.create_rectangle(scale_x*x+2,     scale_y*y+2,
                                 scale_x*(x+1)+2, scale_y*(y+1)+2,
                                 fill=colormap[round(z[y,x])], width=0)
    
    return pic

            

##############################################################################
#
# Plot a matrix.
#
##############################################################################

def plot_matrix(z, r_x=0, r_y=0, filename=0, colorcode=0):

    if filename:
        print "Saving to file not supported."

    import MLab, Tkinter, ImageTk
        
    # Scale z and find appropriate colormap.
    
    zmax = MLab.max(MLab.max(z))
    zmin = MLab.min(MLab.min(z))

    if (zmin < 0) and (0 < zmax) :
        colormap = create_bipolar_colormap()
        zmax = MLab.max([-zmin, zmax])
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
            
    # Put picture on canvas.

    root = Tkinter.Tk()  
    pic = plot_scaled_matrix(root, colormap, z, r_x, r_y)
    pic.pack()
    root.mainloop()



##############################################################################
#
# Creates a movie starting from a complex matrix representing phasors.
#
##############################################################################

def phasormovie(z, r_x=0, r_y=0, filename=0):

    if filename:
        print "Saving to file not supported."
        
    import MLab, Tkinter
    
    # Make movie memory.
    
    movie    = []
    frames   = 16
    colormap = create_bipolar_colormap()

    # Scale factors for z.

    zmax = MLab.max(MLab.max(abs(z)))
    z_scale = (len(colormap)-1)/(2*zmax)

    # Calculate each frame.
    
    root = Tkinter.Tk()    
    for Nr in range(0,frames):
        pic = plot_scaled_matrix(root, colormap, ((z+zmax)*z_scale).real,
                                 r_x, r_y)
        movie.append(pic)
        z *= exp(2j*pi/frames)
        
    # Close window procedure.

    stop = [0]    
    def callback(): stop[0] = 1 
    root.protocol("WM_DELETE_WINDOW", callback)
            
    # Animate picture.

    while not stop[0]:
        for x in range(frames):
            if stop[0]: break
            movie[x].pack()
            root.update()
            root.after(int(100))
            movie[x].forget()

    root.destroy()  



##############################################################################
#
# Different colormaps.
#
##############################################################################

def create_bipolar_colormap():

    colormap = []
    
    # blue-white-1
    for j in range(0,255):
        colormap.append('#%02X%02XFF'%(j,j))
    # white-red 
    for j in arange(255,-1,-1):
        colormap.append('#FF%02X%02X'%(j,j))
  
    return colormap



def create_unipolar_colormap():

    colormap = []
    # black-blue-1
    for j in range(0,255):
        colormap.append('#0000%02X'%j)    
    # blue-purple-1
    for j in range(0,255):
        colormap.append('#%02X00FF'%j)
    # purple-red+1
    for j in arange(255,0,-1):
        colormap.append('#FF00%02X'%j)      
    # red-yellow-1
    for j in range(0,255):
        colormap.append('#FF%02X00'%j)
    # yellow-white       
    for j in range(0,256):
        colormap.append('#FFFF%02X'%j)

    return colormap



def create_white_black_colormap():

    colormap = []
    # white-black              -#FFFFFF-#000000
    for j in arange(255,-1,-1):
        colormap.append('#%02X%02X%02X'%(j,j,j))

    return colormap



def create_black_white_colormap():

    colormap = []
    # black-white              -#000000-#FFFFFF
    for j in range(0,256):
        colormap.append('#%02X%02X%02X'%(j,j,j))

    return colormap



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
    
    for i in range(N()):
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
    
    fz = zeros([len(r_y),len(r_x)], Float)

    for i_y in range(len(r_y)):
      for i_x in range(len(r_x)):
        fz[len(r_y)-1-i_y, i_x] = abs(f(r_x[i_x] + r_y[i_y]*1j))

    plot_matrix(fz, r_x, r_y, filename, colormap)



##############################################################################
#
# Plot the refractive index profile in a waveguide.
#
##############################################################################

def plot_n_waveguide(waveguide, r_x):
    
    v = []
    
    for i_x in range(len(r_x)):
      v.append((r_x[i_x], abs(waveguide.n(Coord(r_x[i_x],0,0)))))
        
    plot_vector(v)



##############################################################################
#
# Plot the refractive index profile in a stack.
#
##############################################################################

def plot_n_stack(stack, r_x, r_z, filename=0, colormap=whiteblack):
    
    n = zeros([len(r_x),len(r_z)], Float)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        n[len(r_x)-1-i_x,i_z] = stack.n(Coord(r_x[i_x], 0, r_z[i_z])).real

    plot_matrix(n, r_z, r_x, filename, colormap)



##############################################################################
#
# Plot the refractive index profile in a Section.
#
##############################################################################

def plot_n_section(stack, r_x, r_y, filename, colormap):
    
    n = zeros([len(r_y),len(r_x)], Float)

    for i_x in range(len(r_x)):
      for i_y in range(len(r_y)):
        n[len(r_y)-1-i_y,i_x] = stack.n(Coord(r_x[i_x], r_y[i_y], 0)).real

    plot_matrix(n, r_x, r_y, filename, colormap)    



##############################################################################
#
# Wrapper for plot_n.
#
##############################################################################

def plot_n(o, r1, r2=0, filename=0, colormap=whiteblack):

    if not r2:
        plot_n_waveguide(o, r1)
    elif type(o) == Stack or type(o) == BlochStack or type(o) == Cavity:
        plot_n_stack(o, r1, r2, filename, colormap)
    elif type(o) == Section:
        plot_n_section(o, r1, r2, filename, colormap)    
    else:
        print "Unsupported argument for plot_n."


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

def plot_field_stack(stack, component, r_x, r_z, filename, colormap):
    
    f = zeros([len(r_x),len(r_z)], Float)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        f[len(r_x)-1-i_x,i_z] = \
              component(stack.field(Coord(r_x[i_x], 0, r_z[i_z])))

    plot_matrix(f, r_z, r_x, filename, colormap)



##############################################################################
#
# Plot the field profile of a section mode.
#
##############################################################################

def plot_field_section_mode(mode, component, r_x, r_y, filename, colormap):
    
    f = zeros([len(r_y),len(r_x)], Float)

    for i_x in range(len(r_x)):
      for i_y in range(len(r_y)):
        f[len(r_y)-1-i_y,i_x] = \
              component(mode.field(Coord(r_x[i_x], r_y[i_y], 0)))

    plot_matrix(f, r_x, r_y, filename, colormap)    



##############################################################################
#
# Wrapper for plot_field.
#
##############################################################################

def plot_field(o, component, r1, r2=0, filename=0, colormap=0):

    if not r2:
        plot_field_waveguide(o, component, r1)
    elif type(o) == Stack or type(o) == BlochMode or type(o) == Cavity:
        plot_field_stack(o, component, r1, r2, filename, colormap)
    elif type(o) == SectionMode:
        plot_field_section_mode(o, component, r1, r2, filename, colormap)
    else:
        print "Unsupported argument for plot_field."
        


##############################################################################
#
# Animate the field profile in a stack.
#
##############################################################################

def animate_field_stack(stack, component, r_x, r_z, filename=0):
    
    f = zeros([len(r_x),len(r_z)], Complex)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        f[len(r_x)-1-i_x,i_z] = \
              component(stack.field(Coord(r_x[i_x], 0, r_z[i_z])))

    phasormovie(f, r_z, r_x, filename)



##############################################################################
#
# Animate the field profile of a section mode.
#
##############################################################################

def animate_field_section_mode(mode, component, r_x, r_y, filename=0):
    
    f = zeros([len(r_y),len(r_x)], Complex)

    for i_x in range(len(r_x)):
      for i_y in range(len(r_y)):
        f[len(r_y)-1-i_y,i_x] = \
              component(mode.field(Coord(r_x[i_x], r_y[i_y], 0)))

    phasormovie(f, r_x, r_y, filename)



##############################################################################
#
# Wrapper for animate_field.
#
##############################################################################

def animate_field(o, component, r1, r2, filename=0):

    if type(o) == Stack or type(o) == BlochMode or type(o) == Cavity:
        animate_field_stack(o, component, r1, r2, filename)
    elif type(o) == SectionMode:
        animate_field_section_mode(o, component, r1, r2, filename)
    else:
        print "Unsupported argument for animate_field."


