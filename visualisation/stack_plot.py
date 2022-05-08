#! /usr/bin/env python

##############################################################################
#
# File:    stack_plot.py
# Author:  Lieven.Vanholme@intec.ugent.be
#
##############################################################################

from camfr import *
from numpy import *
from tkinter import *

# win bug?? of toch ergens code die ik mis
import camfr_PIL

HELPCONTENT = """\
CAMFR plot: some tips, oddities:

    Zoom in on the structure by drawing a rubber box in the
    index profile window.

    To save animations:
    
      - 'play' should be on.
      - Save with '.gif' extension to get an animated gif.
      - Save with '.jpg' extension to get a number of
        jpg files for each frame. 

    If image disapears of the canvas, the size of the area
    was to small. Click restore on the configuration frame.

    Choose mode of Blochstack.
    Only the release of mouse-button-2 will plot the
    described mode. It's hard to get the appropriate
    mode, this is a solution:
    Keep button-2 pressed.
    Click button-1.
    Move the pointer to the left (or to the right)
    Click button-1 again untill wanted mode is reached.
    Release button-2.
"""       

INFOCONTENT = """\
CAMFR plot 1.2
        
For information
contact: Lieven Vanholme
email:   Lieven.Vanholme@intec.ugent.be"""



##############################################################################
#
# class StackPlot
#
#   This class provides all the widgets around the plot canvas,
#   and increases the functionality of the canvas itself.
#
##############################################################################

class StackPlot:
    
    import camfr_PIL

    # Stack values.
    N               = 0         # Number of modes in blochstack.

    # Canvas values.
    MINVALUE        = 1         # The smallest size for the canvas.
                                # (checked by configure)
    MAXVALUE        = 800       # The maximum size of the canvas.
                                # (checked by configure)
    FIELDSIZE       = 220       # The start size for the field canvas.
                                # (used by restore)
    STACKSIZE       = 100       # The start size for the stack canvas.
                                # (used by restore)
    fieldRange      = [ FIELDSIZE, FIELDSIZE, 6, range(1), range(1)]
                                # Plotrange  x-axis, y-axis, pixel.
                                # (of field)
    stackRange      = [ STACKSIZE, STACKSIZE, 1, range(1), range(1)]
                                # Plotrange  x-axis, y-axis, pixel.
                                # (of refractiveindex)

    # Memory values.
    fm              = [[],[],[]]        # Array of field Matrices.
    fm_nr           = 0                 # Current field Matrix number.
    mat_f           = zeros(0,float)    # Field matrix.
    mat_n           = zeros(0,float)    # Index matrix.
    mat_px          = zeros(0,float)    # Poynting-x matrix.
    mat_py          = zeros(0,float)    # Poynting-y matrix.

    # Picture variables.
    PLAY            = 0                 # Stop-play-pause = 0-1-2.
    REDRAWMOVIE     = 0                 # Don't play while redrawing.
    REDRAWPICTURE   = 0
    REDRAWINDEX     = 0
    REDRAWPOYNTING  = 0

    def __init__(self, stack):

        # Tk var.
        self.window = Frame()   # This alows the make of Tk variables.
        
        # About stack.
        self.stack  = stack
        self.width  = stack.width()*(1 - 1e-16)
                                        # A numerical inacuracy between
                                        # blochstack.width and slab.width. 

        if type(stack) == Stack:
            self.top        = stack.length()
            self.bottom     = 0
            
        elif type(stack) == Cavity:
            self.top        =   stack.top_stack().length()
            self.bottom     = - stack.bot_stack().length()
            
        elif type(stack) == BlochStack:
            self.top        = stack.length()
            self.bottom     = 0
            self.blochStack = stack
            self.N          = stack.N()
            if self.N > 0 :
                self.stack  = stack.mode(0)
            self.mode       = IntVar()
            self.mode.set(0)

        self.length = self.top - self.bottom

        # Bipolarmap.
        camfr_PIL.NEGCOLOR  = (0,0,255)         # Blue.
        camfr_PIL.MIDCOLOR  = (255,255,255)     # White.
        camfr_PIL.POSCOLOR  = (255,0,0)         # Red.
        
        # Main window.
        self.window.master.wm_title("CAMFR stack plot")
        self.window.master.wm_iconname("CAMFR")

        self.indexFlag      = BooleanVar()
        self.poyntingFlag   = BooleanVar()
        self.horizontalFlag = BooleanVar()
        self.logPlot        = BooleanVar()
        
        self.window.grid()
        self._horizontalgrid() # Determinates whether horizontal
                               # or vertical layout.
        self._gridSupport()    # Fit the canvas around the image.

        self._makeRange(self.stackRange, self.bottom, self.top, self.width, 0)
        self._makeRange(self.fieldRange, self.bottom, self.top, self.width, 0)

        self.window.master.protocol("WM_DELETE_WINDOW", self._quit)
        self.window.master.bind("<Control_L>q",self._quit)

        # Make widgets.
        self._makeMenuBar()
        self._makeFieldPlotWindow()
        self._makeStackIndexWindow()
        self._makeFieldButtons()
        self._makeMovieButtons()
        
        # Griding.
        self.window.master.config(menu=self.menuBar)
        self._gridAll()

        # Start up.
        self._drawStack()
        self._drawNewRegion()
        self.window.mainloop()

   ###########################################################################
   #
   # All widgets get sticked on the grid, two options are supported,
   # field canvas at the right or at the top.
   #
   ###########################################################################

    def _gridAll(self):
        if self.horizontalFlag.get():
            self.fieldPlotFrame.grid(ipadx=10, ipady=10, row=0, \
                                     columnspan=3, rowspan=1, column=0)
                                # keep rowspan=1 otherwise rowspan stays 3.
            self.fbf.grid(ipadx=10, ipady=10, row=1, column=0)      
            self.wBf.grid(ipadx=10, ipady=10, row=1, column=1)
            self.stackIndexFrame.grid(ipadx=10, ipady=10, row=1, column=2)
        else:
            self.fieldPlotFrame.grid(ipadx=10, ipady=10, row=0, \
                                     columnspan=1, rowspan=3, column=1)
            self.fbf.grid(ipadx=10, ipady=10, row=0, column=0)
            self.wBf.grid(ipadx=10, ipady=10, row=1, column=0)
            self.stackIndexFrame.grid(ipadx=10, ipady=10, row=2, column=0)

   ###########################################################################
   #
   # Makes the top menubar.
   #
   ###########################################################################

    def _makeMenuBar(self):

        self.indexFlag.set(1)
        self.logPlot.set(0)
        
        self.menuBar = Menu()
        
        fileMenu = Menu(self.menuBar,tearoff=0)
        fileMenu.add_command(       label = "Save as raw data",
                             command = lambda:self._save("CAMFRRAWDATA"))
        fileMenu.add_command(       label="Save as picture",
                             command = lambda:self._save("CAMFRPICTURE"))
        fileMenu.add_command(       label="Save as movie",
                             command = lambda:self._save("CAMFRGIFMOVIE"))
        fileMenu.add_command(       label="Save as movieframes",
                             command = lambda:self._save("CAMFRFRAMESMOVIE"))
        fileMenu.add_command(       label="Quit",
                                    command = self._quit)
        
        configMenu = Menu(self.menuBar,tearoff=0)
        configMenu.add_checkbutton( label='show index window',
                                    command = self._drawField,
                                    variable=self.indexFlag)
        configMenu.add_checkbutton( label='show Poynting vector',
                                    command = self._drawField,
                                    variable=self.poyntingFlag)        
        configMenu.add_checkbutton( label='log plot',
                                    command = self._drawLogField,
                                    variable=self.logPlot)        
        configMenu.add_checkbutton( label='horizontal layout',
                                    command = self._gridAll,
                                    variable=self.horizontalFlag)
        configMenu.add_command(     label='configure...',
                                    command = self._configurationFrame)
        configMenu.add_command(     label='color...',
                                    command = self._colorFrame)
        
        helpMenu = Menu(self.menuBar, tearoff=0)
        helpMenu.add_command(       label="Help",
                                    command = self._help )
        helpMenu.add_command(       label="About",
                                    command = self._about )
        
        self.menuBar.add_cascade(   label="File",       menu=fileMenu)
        self.menuBar.add_cascade(   label="Config",     menu=configMenu)
        self.menuBar.add_cascade(   label="Help",       menu=helpMenu)


   ###########################################################################
   #
   # Makes a frame with the play-pause button.
   # More buttons (plain pause or full stop) can easily be added
   #
   ###########################################################################

    def _makeMovieButtons(self):      
        self.wBf    = Frame ()         # WorkButtonFrame.

        self.playpauseB  = Button(self.wBf,text='play ',
                             bd = 0, #border
                             font=( "Courier New", 10, ""),
                             activebackground="#FF5555" ,
                             command = self._moviePlay)

        # Griding.
        self.wBf.grid()
        self.playpauseB.grid(row=0, column=0)
     

   ###########################################################################
   #
   # Construction of a frame which holds the 'E1r, ..' buttons
   #
   ###########################################################################

    def _makeFieldButtons(self):
        
        self.fbf    = Frame()          # FieldButtonFrame.
        self.fbf.grid()
        self.field  = StringVar()

        if ( get_polarisation() == TE ):
            self.field.set("E2r")
            rEr = Radiobutton( self.fbf, text="E2.r", variable=self.field,
                               value="E2r", command = self._drawNewField,
                               anchor=W, height=-1, width=-1)
            rEi = Radiobutton( self.fbf, text="E2.i", variable=self.field,
                               value="E2i", command = self._drawNewField,
                               anchor=W)
            rHr = Radiobutton( self.fbf, text="H1.r", variable=self.field,
                               value="H1r", command = self._drawNewField,
                               anchor=W)
            rHi = Radiobutton( self.fbf, text="H1.i", variable=self.field,
                               value="H1i", command = self._drawNewField,
                               anchor=W)
            rzr = Radiobutton( self.fbf, text="Hz.r", variable=self.field,
                               value="Hzr", command = self._drawNewField,
                               anchor=W)
            rzi = Radiobutton( self.fbf, text="Hz.i", variable=self.field,
                               value="Hzi", command = self._drawNewField,
                               anchor=W)
        else:
            self.field.set("H2r")
            rEr = Radiobutton( self.fbf, text="E1.r", variable=self.field,
                               value="E1r", command = self._drawNewField,
                               anchor=W)
            rEi = Radiobutton( self.fbf, text="E1.i", variable=self.field,
                               value="E1i", command = self._drawNewField,
                               anchor=W)
            rHr = Radiobutton( self.fbf, text="H2.r", variable=self.field,
                               value="H2r", command = self._drawNewField,
                               anchor=W)
            rHi = Radiobutton( self.fbf, text="H2.i", variable=self.field,
                               value="H2i", command = self._drawNewField,
                               anchor=W)
            rzr = Radiobutton( self.fbf, text="Ez.r", variable=self.field,
                               value="Ezr", command = self._drawNewField,
                               anchor=W)
            rzi = Radiobutton( self.fbf, text="Ez.i", variable=self.field,
                               value="Ezi", command = self._drawNewField,
                               anchor=W)

        # Grid.
        rEr.grid(row=0 , column=0, padx=0,pady=0, sticky=W)
        rEi.grid(row=0 , column=1, padx=0,pady=0, sticky=W)
        rHr.grid(row=1 , column=0, padx=0,pady=0, sticky=W)
        rHi.grid(row=1 , column=1, padx=0,pady=0, sticky=W)
        rzr.grid(row=2 , column=0, padx=0,pady=0, sticky=W)
        rzi.grid(row=2 , column=1, padx=0,pady=0, sticky=W)

        # Mode button.
        if self.N:

            def modeVar(var):
                if (self.mode.get() < 0): self.mode.set(0)
                elif (self.mode.get() > (self.N - 1)): self.mode.set(self.N-1)
                self._drawNewMode()   

            def left():
                if (self.mode.get() > 0):
                    self.mode.set( self.mode.get() - 1)
                    self._drawNewMode()                    

            def right():
                if (self.mode.get() < self.N -1):
                    self.mode.set( self.mode.get() + 1)
                    self._drawNewMode()                    

            #mode chooser field button frame
            mcfbf = Frame(self.fbf)
            mcfbf.grid(row=3 , columnspan=2, sticky=W+E+N+S)
            
            leftButton      = Button(mcfbf, text="<",
                                     bd = 0, #border
                                     font=( "", 10, "bold"),
                                     command=left)
            rightButton     = Button(mcfbf, text=">",
                                     bd = 0,
                                     font=( "", 10, "bold"),
                                     command=right)
            
            modeEntry       = Entry(mcfbf, width=3)
            modeEntry.config(textvariable=self.mode)
            modeEntry.bind('<Key-Return>', modeVar)

            leftButton.grid(    row=0 , column=0)
            modeEntry.grid(     row=0 , column=1)
            rightButton.grid(   row=0 , column=2)
         

   ###########################################################################
   #
   # Make the Frame to plot the field in.
   #
   ###########################################################################

    def _makeFieldPlotWindow(self):

        from matrix_plot_canvas import MatrixPlotCanvas
        
        self.fieldPlotFrame = Frame()
        self.fieldPlotFrame.grid()
        
        def fieldZoomer(x0,y0,x1,y1):

            # Convert window coordinates to array numbers.
            xrange  = self.fieldRange[3]           
            yrange  = self.fieldRange[4]
            pixel   = self.fieldRange[2]

            x0, y0 = int(x0/pixel), int(y0/pixel)
            x1, y1 = int(x1/pixel), int(y1/pixel)

            #  Area selected.
            if (x0 != x1) and (y0 != y1): 
                self._makeRange( self.fieldRange,
                                 xrange[x0], xrange[x1],
                                 yrange[y0], yrange[y1])
                self._drawNewRegion()
     
        self.fieldCanvas = MatrixPlotCanvas(self.fieldPlotFrame,
                            self.fieldRange[0], self.fieldRange[1],
                            relief=SUNKEN, border=2, zoom=fieldZoomer)
        
        self.fieldCanvas.grid(row=0,column=0)

   ###########################################################################
   #
   # Make the frame to plot the refractive window in.
   #
   ###########################################################################

    def _makeStackIndexWindow(self):

        from matrix_plot_canvas import MatrixPlotCanvas
        
        self.stackIndexFrame = Frame()
        self.stackIndexFrame.grid()

        def stackZoomer(x0,y0,x1,y1):

            # Fit the window coordinates to array numbers.
            xrange  = self.stackRange[3]
            yrange  = self.stackRange[4]
            pixel   = self.stackRange[2]

            x0, y0 = int(x0/pixel), int(y0/pixel)
            x1, y1 = int(x1/pixel), int(y1/pixel)

            if (x0 != x1) and (y0 != y1): 
                self._makeRange(self.fieldRange,
                                xrange[x0], xrange[x1],
                                yrange[y0], yrange[y1])
                self._drawNewRegion()
     
        self.stackCanvas = MatrixPlotCanvas(self.stackIndexFrame,
                            self.stackRange[0], self.stackRange[1],
                            relief=SUNKEN, border=2, zoom = stackZoomer)
        
        self.stackCanvas.grid(row=0,column=0)

   ###########################################################################
   #
   # Stops the movie.
   #
   ###########################################################################

    def _movieStop(self):
        self.PLAY = 0
        self.playpauseB.config(text="play ")
        self._drawField()

   ###########################################################################
   #
   # Makes a movie and send it to _movie()
   #
   ###########################################################################

    def _moviePlay(self):

        if (self.PLAY == 0):    #stopped
            self._movieNewPlay()
            self.PLAY = 1
            self.playpauseB.config(text="pause")
            self._movie()
        elif(self.PLAY == 1):     #playing
            self.PLAY = 2
            self.playpauseB.config(text="play ")
        elif (self.PLAY == 2):    #paused
            self.PLAY = 1
            self.playpauseB.config(text="pause")


#    def _moviePlay(self):
#        #if new field..
#        self._movieNewPlay()
#
#        if (self.PLAY == 2):
#            self.pauseB.config(bg = "#CCCCCC")
#            
#        self.PLAY = 1
#        self._movie()


    def _movieNewPlay(self):
        #  New field?.
        if self.REDRAWMOVIE:
            self.REDRAWMOVIE = 0
            self.movie = camfr_PIL._create_phasor_movie(
                self.fm[self.fm_nr], min_area=0,
                scale = self.fieldRange[2],
                ln = self.logPlot.get())
        

   ###########################################################################
   #
   # Makes en plays the movie out of a matrix, index is added by _drawpic().
   #
   ###########################################################################

    def _movie(self):       
        import time

        frames , x = len(self.movie), 0
        while not (self.PLAY == 0):
            if (x == frames): x = 0
            if (self.PLAY == 0): break
            if (self.PLAY == 2): self.fieldCanvas.update()
                                # To fix pause bug so far.
            if (self.PLAY == 1):
                self.pic_f = self.movie[x]
                self._drawPic()
                x += 1
            
            time.sleep(0.1)

   ###########################################################################
   #
   # Draw the stack in the refractive index window
   #
   ###########################################################################
            
    def _drawStack(self):
        xr = self.stackRange[3]
        yr = self.stackRange[4]
        n = zeros([len(yr),len(xr)], float)

        camfr_PIL._calc_n_stack(n, self.stack, yr, xr)

        pic = camfr_PIL._create_matrix_plot(n, colorcode=camfr_PIL.whiteblack,
                                            min_area = 0,
                                            scale=self.stackRange[2] )
        self.stackCanvas.draw(pic)

   ###########################################################################
   #
   # Calc the field (depends on what user chose).
   #
   ###########################################################################
        
    def _calcField(self):
        self.REDRAWMOVIE = 1     # Movie pictures will have to be rewritten.
        self.REDRAWPICTURE = 1   # Picture will have to be rewritten.

        if self.PLAY==2 : self.PLAY = 0  # because we don't want a frame from
                                         # the old movie
        
        field   = self.field.get()
        z       = zeros(0,complex)

        def setfield(nr, comp):
            self.fm_nr = nr
            if not len(self.fm[nr]): self.fm[nr] = self._calc_field(comp)  

        if ( get_polarisation() == TE ):
            if (field == "E2r") or (field == "E2i") :
                setfield(0, lambda f : f.E2())
            if (field == "H1r") or (field == "H1i") :
                setfield(1, lambda f : f.H1())
            if (field == "Hzr") or (field == "Hzi") :
                setfield(2, lambda f : f.Hz())

            if (field == "E2r") or (field == "H1r") or (field == "Hzr"):
                z = array(self.fm[self.fm_nr].real)                 
            if (field == "E2i") or (field == "H1i") or (field == "Hzi"):
                z = array(self.fm[self.fm_nr].imag)
        else:    
            if (field == "E1r") or (field == "E1i") :
                setfield(0, lambda f : f.E1())
            if (field == "H2r") or (field == "H2i") :
                setfield(1, lambda f : f.H2())
            if (field == "Ezr") or (field == "Ezi") :
                setfield(2, lambda f : f.Ez())
            
            if (field == "E1r") or (field == "H2r") or (field == "Ezr"):
                z = array(self.fm[self.fm_nr].real)
            if (field == "E1i") or (field == "H2i") or (field == "Ezi"):
                z = array(self.fm[self.fm_nr].imag)  

        self.mat_f = array(z)

   ###########################################################################
   #
   # Calculatess the index for the fieldcanvas.
   #
   ###########################################################################

    def _calcIndex(self):      
        xr = self.fieldRange[3]
        yr = self.fieldRange[4]
        n = zeros([len(yr),len(xr)], float)

        camfr_PIL._calc_n_stack(n, self.stack, yr, xr)

        self.mat_n = array(n)

   ###########################################################################
   #
   # Calculatess the Poynting vector for the fieldcanvas.
   #
   ###########################################################################

    def _calcPoynting(self):    
        xr = self.fieldRange[3]
        yr = self.fieldRange[4]      
        px   = zeros([len(yr),len(xr)], float)
        py   = zeros([len(yr),len(xr)], float)

        if self._hasField():
            camfr_PIL._calc_arrow_stack(py, px, self.stack, yr, xr)

        self.mat_py = py      
        self.mat_px = px

   ###########################################################################
   #
   # If a new field (e.g. E1real) is chosen:
   # this wil calculate and draw the new field 
   #
   ###########################################################################

    def _drawNewField(self): 
        self._calcField()
        self._drawField()

       
   ###########################################################################
   #
   # If a new region is chosen in the stack canvas:
   # calculate new values for the fieldcanvas.
   #
   ###########################################################################

    def _drawNewRegion(self):
        self.REDRAWPOYNTING = 1
        self.REDRAWINDEX    = 1
        self.fm             = [[],[],[]]
        
        self._calcField()
        self._calcIndex()
        self._calcPoynting()
        self._drawField()


   ###########################################################################
   #
   # A new mode, a new poynting (BlochstackCode)
   #
   ###########################################################################
   
    def _drawNewMode(self):

        self.fm  = [[],[],[]] # remove old field
        self.stack = self.blochStack.mode(self.mode.get())

        self.REDRAWPOYNTING = 1
        self._calcPoynting()
        
        self._drawNewField()  # a new mode, a new field


   ###########################################################################
   #
   # Doesn't need to be recalculated, just take log(mat_f) and redraw.
   #
   ###########################################################################

    def _drawLogField(self):

        self.REDRAWPICTURE  = 1
        self.REDRAWMOVIE    = 1
        self._drawField()

   ###########################################################################
   #
   # Drawfield passes the drawing to _moviePlay() or _drawPic().
   #
   ###########################################################################

    def _drawField(self):

        # Make sure not to make again the same pic.
        # Control by REDRAWPICTURE and REDRAMMOVIE.
        
        # Make indexplot.
        if self.REDRAWINDEX and self.indexFlag.get():
                # Prevent self.mat_n to be be overwritten by camfr_PIL.
                mat_n = array(self.mat_n)
                self.pic_n = camfr_PIL._create_matrix_plot( \
                                mat_n, colorcode=camfr_PIL.whiteblack,\
                                min_area = 0, scale=self.fieldRange[2] )
                self.REDRAWINDEX = 0   

        # Make poyntingplot.
        if self.REDRAWPOYNTING and self.poyntingFlag.get():
                mat_px = array(self.mat_px)
                mat_py = array(self.mat_py)
                self.pic_p = camfr_PIL._create_arrow_plot( \
                                mat_py, mat_px, \
                                min_area = 0, scale=self.fieldRange[2] )
                self.REDRAWPOYNTING = 0

        # Movie or pic?.
        if (self.PLAY > 0) :
            if self.REDRAWMOVIE:
               self._movieNewPlay()
            self._movie()
        else:
            if self.REDRAWPICTURE:
                if self.logPlot.get():
                    mat_f   = log( abs(self.fm[self.fm_nr]) +1e-4)
                    # **2 instead of abs().
                    
                    #better should check if matrix is zero?
                    if self._hasField() :  cc = 3 # Palet.
                    else: cc = 1
                          
                else:
                    mat_f   = array(self.mat_f)
                    
                    #better should check if matrix is zero?
                    if self._hasField() :    cc = 0
                    else:                    cc = 1
                    
                self.pic_f  = camfr_PIL._create_matrix_plot( \
                                mat_f, min_area = 0, colorcode=cc, \
                                scale = self.fieldRange[2])
                self.REDRAWPICTURE = 0

            self._drawPic()

   ###########################################################################
   #
   # Passes the right picture with(out) index to the canvas.
   #
   ###########################################################################

    def _drawPic(self):
        if not (self.indexFlag.get() or self.poyntingFlag.get()):
            self.fieldCanvas.draw( self.pic_f)
        elif (self.indexFlag.get() and self.poyntingFlag.get()):
            self.fieldCanvas.draw(
                camfr_PIL._overlay_pictures(
                    self.pic_f, camfr_PIL._overlay_pictures(
                        self.pic_p, self.pic_n, 0), 1))
        elif self.indexFlag.get():
            self.fieldCanvas.draw(
                camfr_PIL._overlay_pictures(self.pic_f, self.pic_n, 1))
        else:
            self.fieldCanvas.draw(
                camfr_PIL._overlay_pictures(self.pic_f, self.pic_p, 0))

   ###########################################################################
   #
   # Calculates the field matrix.
   #
   ###########################################################################

    def _calc_field(self,component):

        xr = self.fieldRange[3]
        yr = self.fieldRange[4]
        f  = zeros([len(yr),len(xr)], complex)

        if self._hasField():
            #a 'popup' appears while calculating 
            calcInfoFrame = infoFrame("calculating...")
            self.stackCanvas.update()#needed to show calcInfoFrame(apparently)
            camfr_PIL._calc_field_stack(f, self.stack, yr, xr, component)
            calcInfoFrame.destroy()

        return f


   ###########################################################################
   #
   # stack:
   # Controls if (inc)field is set. This has consequences for plotting etc..
   # blochstack:
   # control if N is > 0
   #
   ###########################################################################

    def _hasField(self):
        hasField = 0
        tr = 0
        try:
            tr = len(self.stack.inc_field())
        except AttributeError:
            pass
        if ( tr or self.N >0 ): hasField = 1

        return hasField


            
   ###########################################################################
   #
   # Makes the fysical range:
   # Stack area which will be plotted
   #
   #           ^  x0y0--------|
   #           |  |           |
   #           y  |----------x1y1
   #              x -->   
   #
   # Not so consistent choice of axis.
   # Historical reasons.
   #
   ###########################################################################

    def _makeRange(self,plRange,x0,x1,y0,y1):   
        # Fit it the best plotsize.
        w_stack   =   y0-y1
        l_stack   =   x1-x0
        l_window  =   plRange[0]
        w_window  =   plRange[1]
        x_pixel   =   plRange[2]
        
        if (l_stack/l_window > w_stack/w_window):
            cst = l_stack/(l_window/x_pixel)
        else: cst = w_stack/(w_window/x_pixel)

        plRange[3] = arange(x0,x1,cst)
        plRange[4] = arange(y0,y1,-cst)



   ###########################################################################
   #
   # The canvas will be sticked around the stack,
   # keeping the same area as proposed.
   #
   ###########################################################################

    def _gridSupport(self):    
        import math
        # Fit plot window.
        cst = math.sqrt(
            (self.stackRange[0]*self.stackRange[1])/(self.width*self.length))
        self.stackRange[0] = cst * self.length 
        self.stackRange[1] = cst * self.width
        
        cst = math.sqrt(
            (self.fieldRange[0]*self.fieldRange[1])/(self.width*self.length))
        self.fieldRange[0] = cst * self.length 
        self.fieldRange[1] = cst * self.width

   ###########################################################################
   #
   # This function determinates is the stack is horizontal oriented.
   #
   ###########################################################################

    def _horizontalgrid(self):      
        if ((self.width/self.length) < 1):  self.horizontalFlag.set(1)
        else : self.horizontalFlag.set(0)

   ###########################################################################
   #
   # Makes an object of the configurationFrame class.
   #
   ###########################################################################

    def _configurationFrame(self):
        confr = configurationFrame(self)

   ###########################################################################
   #
   # Makes an object of the colorFrame class.
   #
   ###########################################################################

    def _colorFrame(self):
        colfr = colorFrame(self)

   ###########################################################################
   #
   # This object needs a special quit if the movie plays.
   #
   ###########################################################################

    def _quit(self, mesage=0):
        import sys    
        self.PLAY = 0
        self.window.master.destroy()

   ###########################################################################
   #
   # To save theFieldPlot as gif, jpg of xls (as values)
   # if a movie plays, the movie will be saved as gif.
   #
   ###########################################################################

    def _save(self, saveas="CAMFRPICTURE"):
        import tkFileDialog, tkMessageBox, os, gifmaker, Image, ImagePalette
        
        index   = self.indexFlag.get()

        if (saveas=="CAMFRPICTURE"):
            ifile = "myplot.jpg"
            ftypes = [("picture", ".jpg"),
                      ("picture", ".bmp"),
                      ("picture", ".png")]
            ttle = "save plot"
        elif (saveas=="CAMFRRAWDATA"):
            ifile = "myplot.xls"
            ftypes = [("numeric tabbed", ".xls")]
            ttle = "save data"
        elif (saveas=="CAMFRGIFMOVIE"):
            ifile = "mymovie.gif"
            ftypes = [("movie",".gif")]
            ttle = "save movie as gif"
        elif (saveas=="CAMFRFRAMESMOVIE"):
            ifile = "mymovie.jpg"
            ftypes = [("moviefiles",".jpg"),
                      ("moviefiles", ".bmp"),
                      ("moviefiles", ".png")]
            ttle = "save movie as jpg files"
                      
        #  Get file name.
        name    = tkFileDialog.asksaveasfilename(
                        filetypes   = ftypes ,
                        initialfile = ifile,
                        title       = ttle)
        
        (root, ext) = os.path.splitext(name)

        # Save options.
        if (saveas=="CAMFRPICTURE"):
            if not(ext in camfr_PIL.formats):
                ext = '.jpg'
            name = root + ext
            if not index: self.pic_f.save(name)
            else:
                    camfr_PIL._overlay_pictures(
                    self.pic_f, self.pic_n, 1).save(name)
            print("file saved as ", name)

        elif (saveas=="CAMFRRAWDATA"):
            if not (ext in ['.xls','.XLS']):
                ext = '.xls'
            name = root + ext
            self._saveMatrix( name, self.mat_f)
            print("file saved as", name)


        elif (saveas=="CAMFRGIFMOVIE"):
            if not (ext in ['.gif', '.GIF']):
                ext = '.gif'
            name = root + ext

            # copy from playbutton 
            self._movieNewPlay()
    
            movie2 = range(0,len(self.movie))
            for Nr in range(0,len(self.movie)):
                if not index: movie2[Nr] = self.movie[Nr].convert("P")
                else:
                    movie2[Nr] = camfr_PIL._overlay_pictures(
                        self.movie[Nr], self.pic_n, 1).convert("P")
            fp = open(name,"wb")
            gifmaker.makedelta(fp, movie2)
            fp.close()
            print("file saved as", name)

        elif (saveas=="CAMFRFRAMESMOVIE"):
            if not(ext in camfr_PIL.formats):
                ext = '.jpg'

            # this wil make a movie if nessecary
            self._movieNewPlay()

            for Nr in range(0,len(self.movie)):
                nameNR = root + "_%i" %Nr + ext
                if not index: self.movie[Nr].save(nameNR)
                else:
                    camfr_PIL._overlay_pictures(
                    self.movie[Nr], self.pic_n, 1).save(nameNR)

            print("file saved as %i "%len(self.movie),root,"_Nr",ext," frames")

        else:
            tkMessageBox.showwarning(
                title="fail",
                message="no file/data is saved")


   ###########################################################################
   #
   # Use this to save a matrix as a tabbed (xls) file.
   #
   ###########################################################################

    def _saveMatrix(self, name, matrix):     
        File = open(name,"w")
        for row in matrix:
            for column in row:
                File.writelines( [str(column) + '\t'] )
                #print >> File, column, '\t',
            File.writelines( [''] )
            #print >> File, ''
        File.close()

   ###########################################################################
   #
   # About box.
   #
   ###########################################################################

    def _about(self):
        import tkMessageBox
        tkMessageBox.showinfo('About CAMFR stack plot... ', INFOCONTENT )

   ###########################################################################
   #
   # Help file.
   #
   ###########################################################################

    def _help(self):
        contents = HELPCONTENT

        helpFrame = Toplevel()
        helpFrame.title('CAMFR stack plot help')
        text = Text(helpFrame, relief=SUNKEN, width=60, height=24)
        text.focus_set()
        text.insert(0.0, contents)
        scrollbar = Scrollbar(helpFrame)
        scrollbar.pack(fill=Y, side=RIGHT)
        text.pack(fill=BOTH, expand=YES)
        
        text.configure(yscrollcommand=(scrollbar, 'set'))
        scrollbar.configure(command=(text, 'yview'))



##############################################################################
#
# Simple frame, makes it easy to change colors from the bipolarcolormap.
#
##############################################################################
        
class colorFrame: 
    def __init__(self, stackPlotObject):
        import tkColorChooser

        spo = stackPlotObject
        ct = Toplevel()
        ct.title('configure color')
        
        def negcol():
            c = camfr_PIL.NEGCOLOR
            n = tkColorChooser.askcolor(
                color="#%02x%02x%02x"%(c[0],c[1],c[2]))[0]
            if n:
                camfr_PIL.NEGCOLOR = n
                negButton.config( bg="#%02x%02x%02x"%(n[0],n[1],n[2]) )
            ct.focus()
   
        def midcol():
            c = camfr_PIL.MIDCOLOR
            n = tkColorChooser.askcolor(
                color="#%02x%02x%02x"%(c[0],c[1],c[2]))[0]
            if n:
                camfr_PIL.MIDCOLOR = n
                midButton.config( bg="#%02x%02x%02x"%(n[0],n[1],n[2]) )
            ct.focus()

        def poscol():
            c = camfr_PIL.POSCOLOR
            n = tkColorChooser.askcolor(
                color="#%02x%02x%02x"%(c[0],c[1],c[2]))[0]
            if n:
                camfr_PIL.POSCOLOR = n
                posButton.config( bg="#%02x%02x%02x"%(n[0],n[1],n[2]) )
            ct.focus()
        
        def ok():
            spo.REDRAWMOVIE, spo.REDRAWPICTURE = 1,1
            spo._drawField()
            #  grr.. I don't know how to remove this if a movie plays..
            ct.destroy()


        n = camfr_PIL.NEGCOLOR
        m = camfr_PIL.MIDCOLOR
        p = camfr_PIL.POSCOLOR
        negButton       = Button(ct, text="neg color",
                                 bg="#%02x%02x%02x"%(n[0],n[1],n[2]),
                                 command=negcol)
        midButton       = Button(ct, text="mid color",
                                 bg="#%02x%02x%02x"%(m[0],m[1],m[2]),
                                 command=midcol)
        posButton       = Button(ct, text="pos color",
                                 bg="#%02x%02x%02x"%(p[0],p[1],p[2]),
                                 command=poscol)
        okButton        = Button(ct, text="OK", command=ok)
        negButton.pack()
        midButton.pack()        
        posButton.pack()
        okButton.pack()


##############################################################################
#
# class infoFrame
#
#   This class makes a small informationpanel.
#
##############################################################################

class infoFrame:
    def __init__(self, infostring="info"):
        infot       = Toplevel()
        self.infot  = infot
        infot.title('camfr')
        infoLabel   = Label(infot, text=infostring,
                            fg="black",bg="white",
                            font=("Helvetica", 12))
        infoLabel.pack()

    def destroy(self):
        self.infot.destroy()

##############################################################################
#
# class configurationFrame
#
#   This class makes a small panel in order to give the user the oppertunity
#   to change the pictures size.
#
##############################################################################

class configurationFrame:
    def __init__(self, stackPlotObject):
        spo = stackPlotObject
        self.spo = spo
        # Main window.
        ct = Toplevel()
        ct.title('configure plotting areas')
        cff = Frame(ct, borderwidth=2, relief=GROOVE)   # Config field.
        cfs = Frame(ct, borderwidth=2, relief=GROOVE)   # Config stack.
        cfp = Frame(ct, borderwidth=2, relief=GROOVE)   # Config poynting.
        cfb = Frame(ct, borderwidth=2)                  # Config buttons.

        # All entries.
        self.fieldRangeLength   = Entry(cff, width=4)
        self.fieldRangeLength.insert(0, int(spo.fieldRange[0]) )
        self.fieldRangeWidth    = Entry(cff, width=4)
        self.fieldRangeWidth.insert(0, int(spo.fieldRange[1]) )
        self.fieldRangePixel    = Entry(cff, width=4)
        self.fieldRangePixel.insert(0, int(spo.fieldRange[2]) )
        self.stackRangeLength   = Entry(cfs, width=4)
        self.stackRangeLength.insert(0, int(spo.stackRange[0]) )
        self.stackRangeWidth    = Entry(cfs, width=4)
        self.stackRangeWidth.insert(0, int(spo.stackRange[1]) )
        self.stackRangePixel    = Entry(cfs, width=4)
        self.stackRangePixel.insert(0, int(spo.stackRange[2]) )
        self.poyntingLength     = Entry(cfp, width=4)
        self.poyntingLength.insert(0, int(camfr_PIL.ARROWSIZE) )

        # Labels.
        fieldLabel     = Label(cff, text="field canvas")
        stackLabel     = Label(cfs, text="index canvas")
        flengthLabel   = Label(cff, text="Length")
        fwidthLabel    = Label(cff, text="Width")
        fpixelLabel    = Label(cff, text="Pixelsize")
        slengthLabel   = Label(cfs, text="Length")
        swidthLabel    = Label(cfs, text="Width")
        spixelLabel    = Label(cfs, text="Pixelsize")
        poyntingLabel  = Label(cfp, text="Poynting arrow length")
      
        ######################################################################
        #
        # Reconfigure the field and stack canvas, using the user-defined size
        #
        ######################################################################

        def ok():           
            # Read fields.
            mi, ma = StackPlot.MINVALUE, StackPlot.MAXVALUE
            if not self._validate(self.fieldRangeLength,
                                  self.spo.fieldRange, 0, mi, ma): return
            if not self._validate(self.fieldRangeWidth,
                                  self.spo.fieldRange, 1, mi, ma): return
            m = min( self.spo.fieldRange[0] ,self.spo.fieldRange[1])/3
            if not self._validate(self.fieldRangePixel,
                                  self.spo.fieldRange, 2, 1, m): return
            if not self._validate(self.stackRangeLength,
                                  self.spo.stackRange, 0, mi, ma): return
            if not self._validate(self.stackRangeWidth,
                                  self.spo.stackRange, 1, mi, ma): return
            m = min( self.spo.stackRange[0], self.spo.stackRange[1])/3   
            if not self._validate(self.stackRangePixel,
                                  self.spo.stackRange, 2, 1, m): return

            camfr_PIL.ARROWSIZE = int(self.poyntingLength.get())

            # Recalc canvas.
            self.spo._movieStop()
            self.spo.fieldPlotFrame.grid_forget()
            self.spo.stackIndexFrame.grid_forget()
            
            self.spo._makeRange(self.spo.fieldRange,
                                self.spo.fieldRange[3][0],
                                self.spo.fieldRange[3][-1],
                                self.spo.fieldRange[4][0],
                                self.spo.fieldRange[4][-1])
            self.spo._makeRange(self.spo.stackRange,
                                self.spo.stackRange[3][0],
                                self.spo.stackRange[3][-1],
                                self.spo.stackRange[4][0],
                                self.spo.stackRange[4][-1])
            self.spo._makeFieldPlotWindow()
            self.spo._makeStackIndexWindow()

            # New grid-plot.            
            self.spo._gridAll()
            self.spo._drawNewRegion()
            self.spo._drawStack()
            
            ct.destroy()

        ######################################################################
        #
        # Restore will resize the canvas to its original settings. The
        # original settings are StackPlot.FIELDSIZE and StackPlot.STACKSIZE.
        #
        ######################################################################

        def restore():
            self.spo.fieldPlotFrame.grid_forget()
            self.spo.stackIndexFrame.grid_forget()
            self.spo.fieldRange[0] = StackPlot.FIELDSIZE
            self.spo.fieldRange[1] = StackPlot.FIELDSIZE
            self.spo.stackRange[0] = StackPlot.STACKSIZE
            self.spo.stackRange[1] = StackPlot.STACKSIZE
            self.spo._gridSupport()
            self.spo._makeRange(self.spo.stackRange, 0,
                                self.spo.length, self.spo.width, 0)
            self.spo._makeRange(self.spo.fieldRange, 0,
                                self.spo.length, self.spo.width, 0)
            self.spo._gridAll()
            self.spo._drawNewRegion()
            self.spo._drawStack()
            ct.destroy()
            
        def cancel():
            ct.destroy()
        
        okButton        = Button(cfb, text="OK",    command=ok)
        cancelButton    = Button(cfb, text="cancel", command=cancel)
        restoreButton   = Button(cfb, text="restore", command=restore)

        # Grid.
        ct.grid()
        cff.grid(ipadx=10, ipady=10, row=0, column=0)
        cfs.grid(ipadx=10, ipady=10, row=0, column=1)
        cfp.grid(ipadx=10, ipady=10, row=1, column=0, columnspan=2)
        cfb.grid(ipadx=10, ipady=10, row=2, column=0, columnspan=2)
        
        cff.grid()
        fieldLabel.grid(        row=0, column=0, columnspan=2)
        flengthLabel.grid(      row=1, column=0)
        fwidthLabel.grid(       row=2, column=0)
        fpixelLabel.grid(       row=3, column=0)
        self.fieldRangeLength.grid( row=1, column=1)
        self.fieldRangeWidth.grid(  row=2, column=1)
        self.fieldRangePixel.grid(  row=3, column=1)

        cfs.grid()
        stackLabel.grid(        row=0, column=0, columnspan=2)        
        slengthLabel.grid(      row=1, column=0)
        swidthLabel.grid(       row=2, column=0)
        spixelLabel.grid(       row=3, column=0)
        self.stackRangeLength.grid( row=1, column=1)
        self.stackRangeWidth.grid(  row=2, column=1)
        self.stackRangePixel.grid(  row=3, column=1)

        cfp.grid()
        poyntingLabel.grid(        row=0, column=0)
        self.poyntingLength.grid(  row=0, column=1)

        cfb.grid()
        okButton.grid(      row=0, column=0)
        cancelButton.grid(  row=0, column=1)
        restoreButton.grid( row=0, column=2)

   ###########################################################################
   #
   # Validates the entry in the config frame.
   #
   ###########################################################################

    def _validate(self, entry, endResult, nr, mi=1, ma=800):   
        import tkMessageBox
        # Controls if it is integer.
        try:
            result = int(entry.get())
        except ValueError:
            tkMessageBox.showwarning(
                "Illegal value",
                "Not an integer." + "\nPlease try again",
                #parent = ct
            )
            return 0

        # Controls if it is too small.
        if self.spo.MINVALUE is not None and result < mi:
            tkMessageBox.showwarning(
                "Too small",
                "The allowed minimum value is %s. "
                "Please try again." % mi
            )
            entry.focus()
            return 0

        # Controls if it is too big.
        if self.spo.MAXVALUE is not None and result > ma:
            tkMessageBox.showwarning(
                "Too large",
                "The allowed maximum value is %s. "
                "Please try again." % ma,
            )
            entry.focus()
            return 0

        endResult[nr] = result
        
        return 1


