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

def plot_E1(waveguide, mode, PML=0):
    v = []
    for i in range(100):
	x = i * waveguide.width().real / 100.0
	v.append((x,abs(waveguide.mode(mode).field(Coord(x-PML*1j,0,0)).E1())))
    plot_vector(v)
    
def plot_E2(waveguide, mode, PML=0):
    v = []
    for i in range(100):
	x = i * waveguide.width().real / 100.0
	v.append((x,abs(waveguide.mode(mode).field(Coord(x-PML*1j,0,0)).E2())))
    plot_vector(v)

def plot_Ez(waveguide, mode, PML=0):
    v = []
    for i in range(100):
	x = i * waveguide.width().real / 100.0
	v.append((x,abs(waveguide.mode(mode).field(Coord(x-PML*1j,0,0)).Ez())))
    plot_vector(v)

def plot_H1(waveguide, mode, PML=0):
    v = []
    for i in range(100):
	x = i * waveguide.width().real / 100.0
	v.append((x,abs(waveguide.mode(mode).field(Coord(x-PML*1j,0,0)).H1())))
    plot_vector(v)
    
def plot_H2(waveguide, mode, PML=0):
    v = []
    for i in range(100):
	x = i * waveguide.width().real / 100.0
	v.append((x,abs(waveguide.mode(mode).field(Coord(x-PML*1j,0,0)).H2())))
    plot_vector(v)

def plot_Hz(waveguide, mode, PML=0):
    v = []
    for i in range(100):
	x = i * waveguide.width().real / 100.0
	v.append((x,abs(waveguide.mode(mode).field(Coord(x-PML*1j,0,0)).Hz())))
    plot_vector(v)
    
def plot_n(waveguide):
    v = []
    for i in range(100):
	x = i * waveguide.width().real / 100.0
	v.append((x, abs(waveguide.n(Coord(x, 0, 0)))))
    plot_vector(v)

def figure():
    matlab("figure")

def hold_on():
    matlab("hold on")
	
