from camfr_work import *
from numpy import *
import MLab, pymat

##############################################################################
#
# Thin wrapper around pymat.
#
# Some changes by Ross Stanley.
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
# Scatter plot.
#
##############################################################################

def scatter_plot(x, y):
    matlab.put("x", array(x))
    matlab.put("y", array(y)) 
    matlab("plot(x,y,'o')")



##############################################################################
#
# Plot a vector.
#
##############################################################################

def plot_vector(v):
        
    pass
    try:
        is2d = v[0][0]
        matlab.put("v", array(v))
        matlab("plot(v(:,1),v(:,2))")
    except TypeError:
        w = [] 
        for i in range(len(v)):
            w.append([i,v[i]])
        matlab.put("w", array(w))
        matlab("plot(w(:,1),w(:,2))")



##############################################################################
#
# Plot a matrix.
#
# Change 'pcolor' to 'imagesc' and reverse x and y
##############################################################################

def plot_matrix(z, r_x=0, r_y=0, filename=0, colorcode=0):

    if filename:
        print "Saving to file not supported."
    
    matlab.put('z', z)
    if r_x:
        matlab.put('x', r_x)
    if r_y:
        matlab.put('y', r_y)

    if not r_x and not r_y:
        matlab("imagesc(z)")
    else:
        matlab("imagesc(x,y,z)")
   
    if (MLab.min(MLab.min(z)) < 0) and (0 < MLab.max(MLab.max(z))):
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
# Change 'pcolor' to 'imagesc' and reverse x and y
##############################################################################

def phasormovie(z, r_x=0, r_y=0, filename=0):

    if filename:
        print "Saving to file not supported."
    
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
            matlab("imagesc(m)")
        else:
            matlab("imagesc(x,y,m)")

        matlab("caxis(v)")
        
        create_bipolar_color_map()
                
        matlab("axis equal")
        matlab("axis off")
        matlab("shading interp")
        matlab("frames(k) = getframe")

    matlab("movie(frames,100)")



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
    
    for i in range(waveguide.N()):
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
    
    fz = zeros([len(r_y),len(r_x)], float)

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
    
    n = zeros([len(r_x),len(r_z)], float)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        n[len(r_x)-1-i_x,i_z] = stack.n(Coord(r_x[i_x], 0, r_z[i_z])).real

    plot_matrix(n, r_z, r_x, filename, colormap)



##############################################################################
#
# Plot the refractive index profile in a Section.
#
# Added default values for filename and colormap
##############################################################################

def plot_n_section(stack, r_x, r_y, filename=0, colormap=whiteblack):
    
    n = zeros([len(r_y),len(r_x)], float)

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
# Added default values for filename and colormap
##############################################################################

def plot_field_stack(stack, component, r_x, r_z, filename=0, colormap=0,
                     overlay_n=1, contour=1):
    
    f = zeros([len(r_x),len(r_z)], float)

    for i_x in range(len(r_x)):
      for i_z in range(len(r_z)):
        f[len(r_x)-1-i_x,i_z] = \
              component(stack.field(Coord(r_x[i_x], 0, r_z[i_z])))

    plot_matrix(f, r_z, r_x, filename, colormap)



##############################################################################
#
# Plot the field profile of a section mode.
#
# Added default values for filename and colormap
##############################################################################

def plot_field_section_mode(mode, component, r_x, r_y, filename=0, colormap=0,
                            overlay_n=1, contour=1):
    
    f = zeros([len(r_y),len(r_x)], float)

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

def plot_field(o, component, r1, r2=0, filename=0, colormap=0,
               overlay_n=1, contour=1):

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

def animate_field_stack(stack, component, r_x, r_z, filename=0,
                        overlay_n=1, contour=1):
    
    f = zeros([len(r_x),len(r_z)], complex)

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

def animate_field_section_mode(mode, component, r_x, r_y, filename=0,
                               overlay_n=1, contour=1):
    
    f = zeros([len(r_y),len(r_x)], complex)

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

def animate_field(o, component, r1, r2, filename=0, overlay_n=1, contour=1):

    if type(o) == Stack or type(o) == BlochMode or type(o) == Cavity:
        animate_field_stack(o, component, r1, r2, filename)
    elif type(o) == SectionMode:
        animate_field_section_mode(o, component, r1, r2, filename)
    else:
        print "Unsupported argument for animate_field."

