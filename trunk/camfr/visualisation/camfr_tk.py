from TkPlotCanvas import *
from camfr import *
from Numeric import *
from MLab import *

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

def plot2D(z):

    # Create color map.
    
    colormap = [] 
    # black-blue
    for j in range(0,16):
        colormap.append('#0000%x%x'%(j,j))    
    # blue-sea
    for j in range(0,16):
        colormap.append('#00%x%xFF'%(j,j))
    # sea-green
    for j in range(0,16):
        colormap.append('#00FF%x%x'%(15-j,15-j))
    # green-yellow
    for j in range(0,16):
        colormap.append('#%x%xFF00'%(j,j))
    # yellow-red        
    for j in range(0,16):
        colormap.append('#FF%x%x00'%(15-j,15-j))

    # Scale z.

    height = z.shape[1]
    width  = z.shape[0]
    
    zmax = max(max(z))
    zmin = min(min(z))
    z -= zmin
    
    colors = len(colormap)
    
    interval =  zmax - zmin
    if (interval != 0):
        z *= (colors-1)/interval

    # Put picture on canvas.

    root = Tk()
    canvas = Canvas(root, width=width, height=height, bg="White")
    canvas.pack()
    
    for y in range(0,height):
        for x in range(0,width):
            canvas.create_line(x+2,y+2,x+3,y+2, fill=colormap[int(z[x,y])])

    root.mainloop()



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

    plot2D(fz)



##############################################################################
#
# Plot the refractive index profile in a waveguide.
#
##############################################################################

def plot_n(waveguide, steps=100):
    
    v = []
    
    for i in range(steps):
	x = i * waveguide.width() / (1.0*steps)
	v.append((x, abs(waveguide.n(Coord(x, 0, 0)))))
        
    plot_vector(v)



##############################################################################
#
# Plot the refractive index profile in a stack.
#
##############################################################################

def plot_n(stack, r_x, r_z):
    
    n = zeros([len(r_z),len(r_x)], Float)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        n[i_z,i_x] = stack.n(Coord(r_x[i_x], 0, r_z[i_z])).real

    plot2D(n)



##############################################################################
#
# Plot the field profile of a waveguide mode.
#
##############################################################################

def plot_field(waveguide, mode, function, steps=100):
    
    v = []
    
    for i in range(steps):
	x = i * waveguide.width() / (1.0 * steps)
	v.append((x,function(waveguide.mode(mode).field(Coord(x,0,0)))))
        
    plot_vector(v)



##############################################################################
#
# Plot the field profile in a stack.
#
##############################################################################

def plot_field(stack, r_x, r_z, component):
    
    f = zeros([len(r_z),len(r_x)], Float)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        f[i_z,i_x] = component(stack.field(Coord(r_x[i_x], 0, r_z[i_z])))

    plot2D(f)
