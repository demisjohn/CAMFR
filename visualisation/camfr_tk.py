from TkPlotCanvas import *
from camfr import *

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

def plot_neff(waveguide):
    v = []
    for i in range(N()):
	n = waveguide.mode(i).n_eff()
	v.append((n.real, n.imag))
    create_window_and_draw(PolyMarker(v, fillcolor='white'))

def plot_E1(wg, mode, PML=0):
    v = []
    for i in range(100):
	x = i * wg.width().real / 100.0
	v.append((x, abs(wg.mode(mode).field(Coord(x-PML*1j, 0, 0)).E1())))
    plot_vector(v)
    
def plot_E2(wg, mode, PML=0):
    v = []
    for i in range(100):
	x = i * wg.width().real / 100.0
	v.append((x, abs(wg.mode(mode).field(Coord(x-PML*1j, 0, 0)).E2())))
    plot_vector(v)

def plot_Ez(wg, mode, PML=0):
    v = []
    for i in range(100):
	x = i * wg.width().real / 100.0
	v.append((x, abs(wg.mode(mode).field(Coord(x-PML*1j, 0, 0)).Ez())))
    plot_vector(v)

def plot_H1(wg, mode, PML=0):
    v = []
    for i in range(100):
	x = i * wg.width().real / 100.0
	v.append((x, abs(wg.mode(mode).field(Coord(x-PML*1j, 0, 0)).H1())))
    plot_vector(v)
    
def plot_H2(wg, mode, PML=0):
    v = []
    for i in range(100):
	x = i * wg.width().real / 100.0
	v.append((x, abs(wg.mode(mode).field(Coord(x-PML*1j, 0, 0)).H2())))
    plot_vector(v)

def plot_Hz(wg, mode, PML=0):
    v = []
    for i in range(100):
	x = i * wg.width().real / 100.0
	v.append((x, abs(wg.mode(mode).field(Coord(x-PML*1j, 0, 0)).Hz())))
    plot_vector(v)

def plot_n(waveguide):
    v = []
    for i in range(100):
	x = i * waveguide.width().real / 100.0
	v.append((x,abs(waveguide.n(Coord(x, 0, 0)))))
    plot_vector(v)
