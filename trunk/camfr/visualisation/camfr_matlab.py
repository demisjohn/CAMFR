from Numeric import *
from camfr import *
import pymat

# Thin wrapper around pymat

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

# Global Matlab Session object

matlab = Matlab()

# Plotting functions

def plot_vector(v):
    matlab.put("v", array(v))
    matlab("plot(v(:,1),v(:,2))")

def plot_neff(waveguide):
    v = []
    for i in range(N()):
	n_eff = waveguide.mode(i).n_eff() 
	v.append((n_eff.real, n_eff.imag))
    matlab.put("v", array(v))
    matlab("plot(v(:,1),v(:,2),'o')")

def plot_field(waveguide, mode, function, steps=100):
    v = []
    PML = waveguide.get_imag_start_thickness()
    for i in range(steps):
	x = i * waveguide.width().real / (1.0 * steps)
	v.append((x,function(waveguide.mode(mode).field(Coord(x+PML*1j,0,0)))))
    plot_vector(v)
    
def plot_n(waveguide, steps=100):
    v = []
    for i in range(steps):
	x = i * waveguide.width().real / (1.0*steps)
	v.append((x, abs(waveguide.n(Coord(x, 0, 0)))))
    plot_vector(v)

def figure():
    matlab("figure")

def hold_on():
    matlab("hold on")
	
