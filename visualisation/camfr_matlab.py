from Numeric import *
from MLab import *
from camfr import *
import pymat

##############################################################################
#
# Thin wrapper around pymat.
#
##############################################################################

class Matlab:
    def __init__(self):
	self.H = pymat.open()
	
    def __call__(self, s):
        pymat.eval(self.H, s)

    def eval(self, s):
        pymat.eval(self.H, s)

    def put(self, name, array):
	pymat.put(self.H, name, array)

    def get(self, name):
	return pymat.get(self.H, name)

    def close(self):
	pymat.close(self.H)

matlab = Matlab()

def figure():
    matlab("figure")

def hold_on():
    matlab("hold on")

	

##############################################################################
#
# Plot a vector.
#
##############################################################################

def plot_vector(v):
    matlab.put("v", array(v))
    matlab("plot(v(:,1),v(:,2))")



##############################################################################
#
# Plot a matrix.
#
##############################################################################

def plot_matrix(z, r_x=0, r_y=0):
    
    matlab.put('z', z)
    if r_x:
        matlab.put('x', r_x)
    if r_y:
        matlab.put('y', r_y)

    if not r_x and not r_y:
        matlab("pcolor(z)")
    else:
        matlab("pcolor(y,x,z)")
   
    if (min(min(z)) < 0) and (0 < max(max(z))):
        create_bipolar_color_map()
    else:
        matlab("colormap jet")
        
    matlab("axis equal")
    matlab("axis off")
    matlab("shading flat")



##############################################################################
#
# Colormap suited for field plot.
#
##############################################################################

def create_bipolar_color_map():

    matlab("""
    
    L=length(colormap)/2
    a = zeros(2*L,3)

    for I = 1:L,
       a(I,1)=1
       a(I,2)=I/L
       a(I,3)=I/L
       a(I+L,3)=1
       a(I+L,1)=1+(1/L)-I/L
       a(I+L,2)=1+(1/L)-I/L
    end

    colormap(a)
    
    """)



##############################################################################
#
# Creates a movie starting from a complex matrix representing phasors.
#
##############################################################################

def phasormovie(z, r_x=0, r_y=0):
    
    matlab.put('z', z)
    if r_x:
        matlab.put('x', r_x)
    if r_y:
        matlab.put('y', r_y)
        
    matlab("mx = max(max(abs(z)))")
    matlab("v = [-mx,mx]")

    matlab("k = 0")
    for k in range(1,17):
        matlab("k = k + 1")
        matlab("m = real(z*exp(i*(k/16)*2*pi))")
  
        if not r_x and not r_y:
            matlab("pcolor(m)")
        else:
            matlab("pcolor(y,x,m)")

        matlab("caxis(v)")
        
        create_bipolar_color_map()
                
        matlab("axis equal")
        matlab("axis off")
        matlab("shading interp")
        matlab("frames(k) = getframe")

    matlab("movie(frames,100)")



##############################################################################
#
# Plot distribution of effective indices in complex plane.
#
##############################################################################
    
def plot_neff(waveguide):
    v = []
    for i in range(N()):
	n_eff = waveguide.mode(i).n_eff() 
	v.append((n_eff.real, n_eff.imag))
    matlab.put("v", array(v))
    matlab("plot(v(:,1),v(:,2),'o')")



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

    plot_matrix(fz, r_x, r_y)



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

    plot_matrix(n, r_x, r_z)



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

    plot_matrix(f, r_x, r_z)



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

    phasormovie(f, r_x, r_z)
    

