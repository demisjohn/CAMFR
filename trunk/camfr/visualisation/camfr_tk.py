from TkPlotCanvas import *
from camfr import *
from Numeric import *
from MLab import *
from cmath import *

##############################################################################
#
# Plot a vector.
#
##############################################################################

def create_window_and_draw(drawobject):
    window = Frame()
    window.pack(fill=BOTH, expand=YES)

    def display(value):
	print value

    c = PlotCanvas(window,500,500,zoom=1,select=display,relief=SUNKEN,border=2)
    c.pack(side=TOP, fill=BOTH, expand=YES)
    c.draw(drawobject, 'automatic', 'automatic')

    window.mainloop()

def plot_vector(v):
    create_window_and_draw(PolyLine(v))



##############################################################################
#
# Plot a matrix.
#
##############################################################################

def plot_matrix(z):

    # Scale z and find appropriate colormap.
    
    z = flipud(z)

    height = z.shape[0]
    width  = z.shape[1]
    
    zmax = max(max(z))
    zmin = min(min(z))
   
    if (zmin < 0) and (0 < zmax) :
        colormap = create_bipolar_color_map()
        zmax = max([-zmin, zmax])
        zmin = -zmax
    else:
        colormap = create_unipolar_colormap_2()

    interval = zmax - zmin

    z = z - zmin

    colors = len(colormap)
 
    if (interval != 0):
        z *= (colors-1)/interval

    # (re)size picture.

    area, scale = 100000, 1
    if (height*width < area):
        scale = int(math.sqrt(area/(height*width)))

    # Put picture on canvas.

    root = Tk()
    canvas = Canvas(root, width=scale*width, height=scale*height, bg="White")
    canvas.pack()

    for y in range(height):
        for x in range(width):
            canvas.create_rectangle(scale*x+2,       scale*y+2,          \
                                    scale*x+2+scale, scale*y+2+scale,    \
                                    fill=colormap[int(z[y,x])], width=0)
 
    root.mainloop()



##############################################################################
#
# Creates a movie starting from a complex matrix representing phasors.
#
##############################################################################

def phasormovie(z):

    # Make movie-memory.
    
    height = z.shape[0]
    width  = z.shape[1]    
    frames = 16
    movie  = zeros([frames, height, width], Float)

    colormap = create_bipolar_color_map()

    # Fill all frames, rescale the data.
    # The central color ends up where data = zero.

    z = flipud(z)

    for Nr in range(0,frames):
        movie[Nr] = z.real
        z = z * exp(2j*pi/frames)
 
    interval = 2*max(max(abs(z)))
    colors   = len(colormap)
    
    for Nr in range(0,frames):
        if (interval != 0): 
            movie[Nr] = movie[Nr] + interval/2.
            movie[Nr] = movie[Nr] * (colors-1)/interval
        else:
            movie[Nr] = movie[Nr] + colors/2.
    
    # (re)size picture.

    area, scale = 100000, 1
    if (height*width < area):
        scale = int(math.sqrt(area/(height*width)))

    # Put picture on canvas.

    root = Tk()

    canvas = Canvas(root, width=scale*width, height=scale*height, bg="White")
    canvas.pack()

    for y in range(height):
        for x in range(0,width):
            canvas.create_rectangle(scale*x+2,       scale*y+2,          \
                                    scale*x+2+scale, scale*y+2+scale,    \
                                    fill=colormap[int(movie[0,y,x])], width=0)

    # Close window procedure
    
    stop = [0]    
    def callback(): stop[0] = 1 
    root.protocol("WM_DELETE_WINDOW", callback)

    # Animate picture.
    
    while not stop[0]:
        for frame in movie:
            if stop[0]: break
            object = 1   
            for y in range(height):
                for x in range(0, width):
                    canvas.itemconfigure(object,fill=colormap[int(frame[y,x])])
                    object += 1  
            root.update()

    root.destroy()       


    
##############################################################################
#
# Different colormaps.
#
##############################################################################

def create_bipolar_color_map():

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
    # blue-sea-1
    for j in range(0,255):
        colormap.append('#00%02XFF'%j)
    # sea-green+1
    for j in arange(255,0,-1):
        colormap.append('#00FF%02X'%j)
    # green-yellow-1
    for j in range(0,255):
        colormap.append('#%02XFF00'%j)
    # yellow-red        
    for j in arange(255,-1,-1):
        colormap.append('#FF%02X00'%j)

    return colormap

def create_unipolar_colormap_2():

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



##############################################################################
#
# Plot distribution of effective indices in complex plane.
#
##############################################################################
    
def plot_neff(waveguide):
    
    v = []
    
    for i in range(N()):
	n = waveguide.mode(i).n_eff()
	v.append((n.real, n.imag))
        
    create_window_and_draw(PolyMarker(v, fillcolor='white'))


    
##############################################################################
#
# Plot a complex function.
#
##############################################################################

def plot_f(f, r_x, r_y):
    
    fz = zeros([len(r_x),len(r_y)], Float)

    for i_y in range(len(r_y)):
      for i_x in range(len(r_x)):
        fz[i_x,i_y] = abs(f(r_x[i_x] + r_y[i_y]*1j))

    plot_matrix(fz)



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

def plot_n_stack(stack, r_x, r_z):
    
    n = zeros([len(r_x),len(r_z)], Float)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        n[i_x,i_z] = stack.n(Coord(r_x[i_x], 0, r_z[i_z])).real

    plot_matrix(n)



##############################################################################
#
# Wrapper for plot_n.
#
##############################################################################

def plot_n(o, r_x, r_z=0):
    if not r_z:
        plot_n_waveguide(o, r_x)
    if type(o) == Stack:
        plot_n_stack(o, r_x, r_z)



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

def plot_field_stack(stack, component, r_x, r_z):
    
    f = zeros([len(r_x),len(r_z)], Float)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        f[i_x,i_z] = component(stack.field(Coord(r_x[i_x], 0, r_z[i_z])))

    plot_matrix(f)



##############################################################################
#
# Wrapper for plot_field.
#
##############################################################################

def plot_field(o, component, r_x, r_z=0):
    if not r_z:
        plot_field_waveguide(o, component, r_x)
    if type(o) == Stack:
        plot_field_stack(o, component, r_x, r_z)



##############################################################################
#
# Animate the field profile in a stack.
#
##############################################################################

def animate_field(stack, component, r_x, r_z):
    
    f = zeros([len(r_x),len(r_z)], Complex)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        f[i_x,i_z] = component(stack.field(Coord(r_x[i_x], 0, r_z[i_z])))

    phasormovie(f)
