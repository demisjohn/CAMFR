from camfr_work import *
from Numeric import *

# Colormap codes.

whiteblack  = 1
blackwhite  = 2
palet       = 3

# PIL file formats.

formats = ['.gif', '.GIF', '.jpg', '.JPG', '.jpeg', '.JPEG', '.bmp', '.BMP',  
           '.dib', '.DIB', '.png', '.PNG', '.tiff', '.TIFF', '.tif', '.TIF',
           '.eps', '.EPS', '.ps', '.PS', '.pdf', '.PDF', '.xbm', '.XBM']

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
	print value

    c = PlotCanvas(window,500,500,zoom=1,select=display,relief=SUNKEN,border=2)
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

def _create_scaled_matrix_plot(colormap, z, r_x=0, r_y=0):
    
    import Image
    
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
    
    pic = Image.new("RGB",(width*scale_x,height*scale_y))
    for x in range(width):
        for y in range(height):
            pic.paste(colormap[round(z[y,x])],
                      (scale_x* x,   scale_y* y,
                       scale_x*(x+1),scale_y*(y+1)))

    return pic

            

##############################################################################
#
# _create_matrix_plot
#
#  Low-level routine. Scales matrix and sets colormap.
#
##############################################################################

def _create_matrix_plot(z, r_x=0, r_y=0, colorcode=0):

    import MLab
        
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

    return _create_scaled_matrix_plot(colormap, z, r_x, r_y)



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
            print "File format not supported. Defaulting to jpg."
            pic_path += ".jpg"

        pic.save(pic_path)    
        print "Created", pic_path
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
# Create a phasor movie.
#
##############################################################################

def _create_phasor_movie(z, r_x=0, r_y=0):
    
    import MLab
    
    # Make movie memory.
    
    movie    = []
    frames   = 16
    colormap = create_bipolar_colormap()

    # Scale factors for z.

    zmax = MLab.max(MLab.max(abs(z)))
    z_scale = (len(colormap)-1)/(2*zmax)

    # Calculate each frame.
    
    for Nr in range(0,frames):
        pic = \
          _create_scaled_matrix_plot(colormap,((z+zmax)*z_scale).real,r_x,r_y)
        movie.append(pic)
        z *= exp(2j*pi/frames)

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
                print "File format not supported. Defaulting to gif."
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
        print  "Created", userpath
        
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
            for x in range(frames):
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

def phasormovie(z, r_x=0, r_y=0, filename=0):

    _output_movie(_create_phasor_movie(z, r_x, r_y), filename)



##############################################################################
#
# Different colormaps.
#
##############################################################################

def create_bipolar_colormap():

    colormap = []
    
    # blue-white-1              -#FF0000-#FFFEFE
    for j in range(0,255): 
        colormap.append( 0xFF0000 + j*0x101 )
    # white-red                 -#FFFFFF-#0000FF
    for j in arange(255,-1,-1):
        colormap.append( 0xFF + j*0x10100 )
  
    return colormap



def create_unipolar_colormap():

    colormap = []
    # black-blue-1              -#000000-#FE0000
    for j in range(0,255):
        colormap.append( j*0x10000 )    
    # blue-purple-1             -#FF0000-#FF00FF
    for j in range(0,255):
        colormap.append( 0xFF0000 + j)
    # purple-red+1              -#FF00FF-#0100FF
    for j in arange(255,0,-1):
        colormap.append( 0x0000FF + j*0x10000)      
    # red-yellow-1              -#0000FF-#00FEFF
    for j in range(0,255):
        colormap.append( 0xFF + j*0x100)
    # yellow-white              -#00FFFF-#FFFFFF
    for j in range(0,256):
        colormap.append( 0xFFFF + j*0x10000)

    return colormap



def create_white_black_colormap():

    colormap = []
    # white-black              -#FFFFFF-#000000
    for j in arange(255,-1,-1):
        colormap.append( j*0x10101 )    

    return colormap



def create_black_white_colormap():

    colormap = []
    # black-white              -#000000-#FFFFFF
    for j in range(0,256):
        colormap.append( j*0x10101 )    

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

def plot_field_stack(stack, component, r_x, r_z, filename, colormap,
                     overlay_n=1, contour=1):

    f = zeros([len(r_x),len(r_z)], Float)

    for i_x in range(len(r_x)):
        for i_z in range(len(r_z)):
            f[len(r_x)-1-i_x,i_z] = \
                  component(stack.field(Coord(r_x[i_x], 0, r_z[i_z])))

    if not overlay_n:
        return plot_matrix(f, r_z, r_x, filename, colormap)

    # Overlay index profile.
        
    n = zeros([len(r_x),len(r_z)], Float)
        
    for i_x in range(len(r_x)):
        for i_z in range(len(r_z)):
            n[len(r_x)-1-i_x,i_z] = stack.n(Coord(r_x[i_x], 0, r_z[i_z])).real

    pic_n = _create_matrix_plot(n, r_z, r_x, whiteblack)
    pic_f = _create_matrix_plot(f, r_z, r_x, colormap)

    _output_pic(_overlay_pictures(pic_f, pic_n, contour), filename)



##############################################################################
#
# Plot the field profile of a section mode.
#
##############################################################################

def plot_field_section_mode(mode, component, r_x, r_y, filename, colormap,
                            overlay_n=1, contour=1):
    
    f = zeros([len(r_y),len(r_x)], Float)
    
    for i_x in range(len(r_x)):
      for i_y in range(len(r_y)):
        f[len(r_y)-1-i_y,i_x] = \
              component(mode.field(Coord(r_x[i_x], r_y[i_y], 0)))

    if not overlay_n:
        return plot_matrix(f, r_x, r_y, filename, colormap)

    # Overlay index profile.
        
    n = zeros([len(r_y),len(r_x)], Float)
        
    for i_x in range(len(r_x)):
        for i_y in range(len(r_y)):
            n[len(r_y)-1-i_y,i_x] = mode.n(Coord(r_x[i_x], r_y[i_y], 0)).real

    pic_n = _create_matrix_plot(n, r_x, r_y, whiteblack)
    pic_f = _create_matrix_plot(f, r_x, r_y, colormap)

    _output_pic(_overlay_pictures(pic_f, pic_n, contour), filename)



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
        plot_field_stack(o,component,r1,r2,filename,colormap,overlay_n,contour)
    elif type(o) == SectionMode:
        plot_field_section_mode(o,component,r1,r2,filename,colormap,overlay_n,
                                contour)
    else:
        print "Unsupported argument for plot_field."
        


##############################################################################
#
# Animate the field profile in a stack.
#
##############################################################################

def animate_field_stack(stack, component, r_x, r_z, filename=0,
                        overlay_n=1, contour=1):

    f = zeros([len(r_x),len(r_z)], Complex)

    for i_x in range(len(r_x)):
        for i_z in range(len(r_z)):
            f[len(r_x)-1-i_x,i_z] = \
                  component(stack.field(Coord(r_x[i_x], 0, r_z[i_z])))

    if not overlay_n:
        return phasormovie(f, r_z, r_x, filename)

    # Overlay index profile.
        
    n = zeros([len(r_x),len(r_z)], Float)
        
    for i_x in range(len(r_x)):
        for i_z in range(len(r_z)):
            n[len(r_x)-1-i_x,i_z] = stack.n(Coord(r_x[i_x], 0, r_z[i_z])).real

    pic_n = _create_matrix_plot (n, r_z, r_x, whiteblack)
    mov_f = _create_phasor_movie(f, r_z, r_x)

    for Nr in range(len(mov_f)):
        mov_f[Nr] = _overlay_pictures(mov_f[Nr], pic_n, contour)

    _output_movie(mov_f, filename)        



##############################################################################
#
# Animate the field profile of a section mode.
#
##############################################################################

def animate_field_section_mode(mode, component, r_x, r_y, filename=0,
                               overlay_n=1, contour=1):
    
    f = zeros([len(r_y),len(r_x)], Complex)

    for i_x in range(len(r_x)):
      for i_y in range(len(r_y)):
        f[len(r_y)-1-i_y,i_x] = \
              component(mode.field(Coord(r_x[i_x], r_y[i_y], 0)))

    if not overlay_n:
        return phasormovie(f, r_x, r_y, filename)

    # Overlay index profile.
        
    n = zeros([len(r_y),len(r_x)], Float)
        
    for i_x in range(len(r_x)):
        for i_y in range(len(r_y)):
            n[len(r_y)-1-i_y,i_x] = mode.n(Coord(r_x[i_x], r_y[i_y], 0)).real

    pic_n = _create_matrix_plot (n, r_x, r_y, whiteblack)
    mov_f = _create_phasor_movie(f, r_x, r_y)

    for Nr in range(len(mov_f)):
        mov_f[Nr] = _overlay_pictures(mov_f[Nr], pic_n, contour)

    _output_movie(mov_f, filename)   



##############################################################################
#
# Wrapper for animate_field.
#
##############################################################################

def animate_field(o, component, r1, r2, filename=0, overlay_n=1, contour=1):

    if type(o) == Stack or type(o) == BlochMode or type(o) == Cavity:
        animate_field_stack(o,component,r1,r2,filename,overlay_n,contour)
    elif type(o) == SectionMode:
        animate_field_section_mode(o,component,r1,r2,filename,overlay_n,
                                   contour)
    else:
        print "Unsupported argument for animate_field."
        
