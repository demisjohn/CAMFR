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

def plot_field(waveguide, mode, function, steps=100):
    v = []
    for i in range(steps):
	x = i * waveguide.width() / (1.0 * steps)
	v.append((x,function(waveguide.mode(mode).field(Coord(x,0,0)))))
    plot_vector(v)
    
def plot_n(waveguide, steps=100):
    v = []
    for i in range(steps):
	x = i * waveguide.width() / (1.0*steps)
	v.append((x, abs(waveguide.n(Coord(x, 0, 0)))))
    plot_vector(v)
