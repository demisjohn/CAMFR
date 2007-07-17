#! /usr/bin/env python

##############################################################################
#
# File:    slab_plot.py
# Author:  Lieven.Vanholme@intec.ugent.be
#
##############################################################################

from camfr import *
from numpy import *
from Tkinter import *
from TkPlotCanvas import *



##############################################################################
#
# class SlabPlot
#
#   Provides all the user interactivity of the mode-plot-screen.
#
##############################################################################

class SlabPlot:

   slab, width = 0, 0                      # TM
   colors  = ['black','red','#0000FF']     # Color for real, imag and neff.
   markers = [0,0]                         # Point to the plotlines.
   plRange = ()                            # Plotrange.
   n_effs  = []                            # Keep n_eff of al modes.

   def __init__(self,slab):

        # Main window.
        self.window = Frame()
        self.window.grid()
        self.window.master.title("CAMFR slab modeplot")

        # About slab.
        self.mode = IntVar()
        self.mode.set(0)
        self.slab = slab
        self.width = slab.width()
        self.plRange = (0, self.width)

        # Make widgets.
        self._makeMenuBar()
        self.makeSlider()
        self.makePlotWindows()
        cf = self.makeChooseField()
        wb = self.makeWorkButtons()
        self.makeValuesField()

        # Griding.
        self.window.master.config(menu=self.menuBar)
        
        cf.grid(row=3,rowspan=2, column=1, sticky=N)
        self.plotFrame.grid(row=0,rowspan=2,column=1, columnspan=2, sticky=NW)
        wb.grid(row=4,column=2)
        self.slider.grid(row=2,column=1,columnspan=2, sticky=NW)
        self.values.grid(row=3,column=2)

        self.window.mainloop()


   ###########################################################################
   #
   # Menubar
   #
   ###########################################################################

   def _makeMenuBar(self):
        self.menuBar = Menu()
        
        fileMenu = Menu(self.menuBar,tearoff=0)
        fileMenu.add_command( label="Quit",
                              command = self.window.master.destroy)
        fileMenu.add_command( label = "Save as raw data",
                              command = lambda:self._save())
                              
        configMenu = Menu(self.menuBar,tearoff=0)
        configMenu.add_checkbutton( label='show effective index distribution',
                                    command = self.changeViewMode)
   
        self.menuBar.add_cascade( label="File",   menu=fileMenu)
        self.menuBar.add_cascade( label="Config", menu=configMenu)


   ###########################################################################
   #
   # Builds the two canvases.
   # The canvas where one mode is plotted 'c'.
   # And a canvas that gives all the different modes 'd'.
   #
   ###########################################################################
   
   def makePlotWindows(self):
        
        self.plotFrame = Frame()
        
        def display(value):
            self.c.select(value)
            self.plRange = value
            self.draw()

        # Mode window.
        self.c = PlotCanvas(self.plotFrame, "300", "200",
                            relief=SUNKEN, border=2,
                            zoom = 1, select = display)
        self.c.grid(row=0,column=0)

        def takeXY(value):
            self.goToMode(array([value[0],value[1]]))
            self.draw()
            self.drawNeff()

        # Neff window.
        self.d = PlotCanvasXY(self.plotFrame, "300", "200", relief=SUNKEN,
                              border=2, zoom = 1, giveXY = takeXY)

   ###########################################################################
   #
   # Make the work buttons (more, reset, quit).
   #     
   ###########################################################################

   def makeWorkButtons(self):

        wB = Frame ()
        clearB = Button(wB, text='Reset Zoom', command = self.resetPlot)

        clearB.pack(side = LEFT)
        
        return wB

   ###########################################################################
   #
   # Reset the plot.
   #
   ###########################################################################

   def resetPlot(self):
        self.plRange = (0, self.width)
        self.draw()
        self.drawNeff()

   ###########################################################################
   #
   # Switch between one or two plot canvases.
   #
   ###########################################################################
        
   def changeViewMode(self):
        if (len(self.plotFrame.grid_slaves())==2):
            self.d.grid_forget()
        else:
            self.d.grid(row=1,column=0)
            self.drawNeff()

   ###########################################################################
   #
   # Make a mode-choose slider.
   #
   ###########################################################################

   def makeSlider(self):
        self.slider = Scale(from_=0, to=self.slab.N()-1 ,
                            orient  = HORIZONTAL,
                            length  = "300",
                            width   = "5", 
                            variable = self.mode,
                            command = self.setSlideVar)

   def setSlideVar(self, val):
        self.mode.set(int(val))
        self.drawNeff()
        self.draw()

   ###########################################################################
   #
   # Draw the Neff of all the modes on canvas 'd'.
   #
   ###########################################################################

   def drawNeff(self):
        if (self.markers[0]==0): self.getNeffMarkers()
        self.setMarker()

        self.d.clear()
        self.d.draw(PlotGraphics(self.markers), 'automatic', 'automatic')

   ###########################################################################
   #
   # Sets a Marker in canvas 'd'.
   # If user chooses a certain mode via the slider.
   #
   ###########################################################################

   def setMarker(self):
        n = self.slab.mode(self.mode.get()).n_eff()
        self.markers[1] = PolyMarker([[n.real,n.imag]], color='red',
                                     fillcolor='red')

        # N_eff-text.
        self.n_eff.set("%f  %fj" %(n.real, n.imag))

   ###########################################################################
   #
   # Get Neff markers of all modes.
   #
   ###########################################################################

   def getNeffMarkers(self):
        ne = self.n_effs

        # Find al n_eff.
        for i in range(self.slab.N()):
            n = self.slab.mode(i).n_eff()
            ne.append([n.real,n.imag])

        self.markers[0] = PolyMarker(ne, color=self.colors[2],
                                     fillcolor='white')

   ###########################################################################
   #
   # Make a frame which shows the values.
   #
   ###########################################################################

   def makeValuesField(self):
        vF = Frame ()
        
        self.n_eff  = StringVar()

        n_effLabel  = Label (vF,text="n_eff", justify=LEFT)
        self.n_eff.set(self.slab.mode(self.mode.get()).n_eff())
        self.n_effValue = Label(vF,textvariable=self.n_eff,
                                width =30, justify=LEFT)

        vF.grid()
        n_effLabel.grid(row=0,column=0, sticky=W)
        self.n_effValue.grid(row=0,column=1, sticky=W)

        self.values = vF

   ###########################################################################
   #
   # Make all the different checkbuttons where users can choose a field.
   #
   ###########################################################################
   
   def makeChooseField(self):
        fF = Frame ()
        fF.grid()
        
        self.varEr      = IntVar()
        self.varEi      = IntVar()
        self.varHr      = IntVar()
        self.varHi      = IntVar()
        self.varZr      = IntVar()
        self.varZi      = IntVar()
        self.varNr      = IntVar()
        self.varNi      = IntVar()

        if ( get_polarisation() == TE ):
            cEr     = Checkbutton(fF, text="E2.r",    variable=self.varEr,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[0])
            cEi     = Checkbutton(fF, text="E2.i",    variable=self.varEi,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[1])
            cHr     = Checkbutton(fF, text="H1.r",    variable=self.varHr,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[0])
            cHi     = Checkbutton(fF, text="H1.i",    variable=self.varHi,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[1])
            czr     = Checkbutton(fF, text="Hz.r",    variable=self.varZr,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[0])
            czi     = Checkbutton(fF, text="Hz.i",    variable=self.varZi,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[1])
            cEr.select()
        else:
            cEr     = Checkbutton(fF, text="E1.r",    variable=self.varEr,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[0])
            cEi     = Checkbutton(fF, text="E1.i",    variable=self.varEi,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[1])
            cHr     = Checkbutton(fF, text="H2.r",    variable=self.varHr,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[0])
            cHi     = Checkbutton(fF, text="H2.i",    variable=self.varHi,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[1])
            czr     = Checkbutton(fF, text="Ez.r",    variable=self.varZr,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[0])
            czi     = Checkbutton(fF, text="Ez.i",    variable=self.varZi,
                                  command=self.draw, padx=0, pady=0,
                                  foreground=self.colors[1])
            cHr.select()

        cNr      = Checkbutton(fF, text="n.r",     variable=self.varNr,
                               command=self.draw, padx=0, pady=0,
                               foreground=self.colors[0])
        cNi      = Checkbutton(fF, text="n.i",     variable=self.varNi,
                               command=self.draw, padx=0, pady=0,
                               foreground=self.colors[1])

        

        cEr.grid(row=0, column=0, padx=0, pady=0, sticky=W)
        cEi.grid(row=0, column=1, padx=0, pady=0, sticky=W)
        cHr.grid(row=1, column=0, padx=0, pady=0, sticky=W)
        cHi.grid(row=1, column=1, padx=0, pady=0, sticky=W)
        czr.grid(row=2, column=0, padx=0, pady=0, sticky=W)
        czi.grid(row=2, column=1, padx=0, pady=0, sticky=W)
        cNr.grid(row=3, column=0, padx=0, pady=0, sticky=W)
        cNi.grid(row=3, column=1, padx=0, pady=0, sticky=W)
        
        return fF


   ###########################################################################
   #
   # Draw fields (lines) on canvas 'c'.
   #
   ###########################################################################

   def draw(self):
        self.lines = []
        self.c.clear()
        Er, Ei = self.varEr.get(), self.varEi.get()
        Hr, Hi = self.varHr.get(), self.varHi.get()
        Zr, Zi = self.varZr.get(), self.varZi.get()
        Nr, Ni = self.varNr.get(), self.varNi.get()

        if Nr: self.getNlines(lambda f : f.real, 0)
        if Ni: self.getNlines(lambda f : f.imag, 0)

        if ( get_polarisation() == TE):
            if Er: self.getFieldLines( lambda f : f.E2().real,0)
            if Ei: self.getFieldLines( lambda f : f.E2().imag,1)
            if Hr: self.getFieldLines( lambda f : f.H1().real,0)
            if Hi: self.getFieldLines( lambda f : f.H1().imag,1)
            if Zr: self.getFieldLines( lambda f : f.Hz().real,0)
            if Zi: self.getFieldLines( lambda f : f.Hz().imag,1)
        else:
            if Er: self.getFieldLines( lambda f : f.E1().real,0)
            if Ei: self.getFieldLines( lambda f : f.E1().imag,1)
            if Hr: self.getFieldLines( lambda f : f.H2().real,0)
            if Hi: self.getFieldLines( lambda f : f.H2().imag,1)
            if Zr: self.getFieldLines( lambda f : f.Ez().real,0)
            if Zi: self.getFieldLines( lambda f : f.Ez().imag,1)
        
        if self.lines:         
            self.c.draw(PlotGraphics(self.lines), 'automatic', 'automatic')

   ###########################################################################
   #
   # Get the plot of the refractive index (as a number of points).
   #
   ###########################################################################
   
   def getNlines(self,cmp,col):
        r, pr = [], self.plRange
        
        for x in arange(pr[0],pr[1],(pr[1]-pr[0])/100.):
            r.append((x, cmp(self.slab.n(Coord(x,0,0)))  ))
            
        self.lines.append(PolyLine(r, color=self.colors[col]))

   ###########################################################################
   #
   # Get the plot of the field (as a number of points).
   #
   ###########################################################################

   def getFieldLines(self,cmp, col):        
        r, mo, sl, pr = [], self.mode.get(), self.slab, self.plRange

        for x in arange(pr[0],pr[1],(pr[1]-pr[0])/100.):
            r.append(( x, cmp(sl.mode(mo).field(Coord(x,0,0))) ))

        self.lines.append(PolyLine(r, color=self.colors[col]))     
  
   ###########################################################################
   #
   # Find mode with x y close to n_eff (in canvas 'd').
   #
   ###########################################################################
   
   def goToMode(self,xy):
        """Find mode with x y close to n_eff."""
        neff, mo = array(self.n_effs), self.mode.get()

        # Scan distance.
        d = sum((neff[0]-xy)**2)
        for n in range(0,len(neff)):
            d2 = sum((neff[n]-xy)**2)
            if (d >= d2): d, mo = d2, n

        self.mode.set(mo)

   ###########################################################################
   #
   # Save this mode
   #
   ###########################################################################

   def _save(self):
        File = open("mymode.xls","w")
        PlotGraphics(self.lines).writeToFile(File,'\t')
        File.close()

##############################################################################
#
# class PlotCanvasXY
#
#   Inherits from PlotCanvas, overriding the showvalue function.
#
##############################################################################

class PlotCanvasXY(PlotCanvas):
    
   def __init__(self, master, width, height, background='white', **attr):

         if attr.has_key('giveXY'):
            self.giveXY = attr['giveXY']
            del attr['giveXY']

         apply(PlotCanvas.__init__,
              (self, master, width, height, background), attr)
         
         #enables mouse selecting
         self.canvas.bind('<ButtonRelease-1>', self._mouseRelease2)


   def _showValue(self, event):
        scale, shift = self.transformation
        x       = self.canvas.canvasx(event.x)
        y       = self.canvas.canvasy(event.y)
        point   = numpy.array([x, y])
        point   = (point-shift)/scale
        self.giveXY(point)

   ###########################################################################
   #
   #  _mouseRelease2 covers _mouseRelease and enables object selecting
   #  with the left mouse button
   #   
   ###########################################################################

   def _mouseRelease2(self, event):
      if ( (self.startx == self.canvas.canvasx(event.x)) and
           (self.starty == self.canvas.canvasy(event.y))):
         self._showValue(event)
      else:
         self._mouseRelease(event)
