#! /usr/bin/env python

##############################################################################
#
# File:    matrix_plot_canvas.py
# Author:  Lieven.Vanholme@intec.ugent.be
#
##############################################################################

from Tkinter import *

##############################################################################
#
# class MatrixPlotCanvas
#
#  Provides the frame with a Canvas with extra functionality:
#  draw function and mouse bindings.
#
##############################################################################

class MatrixPlotCanvas(Frame):
    
    global myImage              # Global to keep picture on canvas.
    
    def __init__(self, master, width, height, **attr):

        self.zoom = 0
        if attr.has_key('zoom'):
            self.zoom = attr['zoom']
            del attr['zoom']
            
        self.selectfn = None
        if attr.has_key('select'):
            self.selectfn = attr['select']
            del attr['select']
            
        apply(Frame.__init__, (self, master), attr)

        self.canvas = Canvas(self, width=width, height=height)
        self.canvas.pack(fill=BOTH, expand=YES)
        border_w    = self.canvas.winfo_reqwidth() - \
                       int(self.canvas.cget('width'))
        border_h    = self.canvas.winfo_reqheight() - \
                       int(self.canvas.cget('height'))
        self.border = (border_w, border_h)

        if self.zoom or self.selectfn is not None:
            self.mouse_state = 0
            self.canvas.bind('<Button-1>', self._mousePressed)
            self.canvas.bind('<B1-Motion>', self._mouseMotion)
            self.canvas.bind('<ButtonRelease-1>', self._mouseRelease)

        self.canvas.bind('<Button-2>', self._showValue)
        self.canvas.bind('<ButtonRelease-2>', self._hideValue)
 
    	self.rubberband = None
    	self.rectangle  = None

    ##########################################################################
    #
    # draw
    #
    #  Puts the image on the canvas.
    #  If smaller, the image will be placed in the middle of the canvas.
    #
    ##########################################################################

    def draw(self,pic):

        import Tkinter, ImageTk

        # Old band has no connection with new draw.
        if self.rubberband:
            self.canvas.delete(self.rubberband) 

        self.myImage=ImageTk.PhotoImage(pic)

        # Center the image.
        im_w = self.myImage.width()        
        im_h = self.myImage.height()
        ca_w = int(self.canvas.cget('width'))
        ca_h = int(self.canvas.cget('height'))
        self.yborder = int(( ca_h - im_h )/ 2)
        self.xborder = int(( ca_w - im_w )/ 2)
        
        self.canvas.create_image(self.xborder, self.yborder,
                                 image=self.myImage, anchor=Tkinter.NW)
        self.canvas.update()

    ##########################################################################
    #
    # mousePressed
    #
    #  Used as callback function to position the user's mouse on the canvas.
    #
    ##########################################################################
        
    def _mousePressed(self,event):
        self.startx = self.canvas.canvasx(event.x)
        self.starty = self.canvas.canvasy(event.y)


    ##########################################################################
    #
    # mouseMotion
    #
    #  If the mouse moves over the canvas, a rectangle (rubberband)
    #  will appear.
    #
    ##########################################################################

    def _mouseMotion(self,event):
        x = self.canvas.canvasx(event.x)
        #same? x, event.x (float and int)
        y = self.canvas.canvasy(event.y)

        if (self.startx != event.x)  and (self.starty != event.y) : 
                self.canvas.delete(self.rubberband)
                self.rubberband = self.canvas.create_rectangle(
                    self.startx, self.starty, x, y)
        #self.update_idletasks() ?


    ##########################################################################
    #
    # mouseRelease
    #
    # As long as the mouse moves, the rubberband can be anywhere in the
    # canvas. At the end, the rubberband should only be part of the image.
    #
    ##########################################################################
   
    def _mouseRelease(self,event):
        """
        (imgx0,imgy0)-----------|
        |                       |
        |                       |
        ------------(imgx1, imgy1)
                                       (rectx0,recty0)------------|
                                        |                         |
                                        |                         |
                                        ------------(rectx1, recty1)

        It's somehow tricky to substract '1' from img heigth and width
        but an image with length '1' has just one pixel, say x0=0, x1=0
        """
        # Coordinates of image and rectangle.
        rectx0 = self.startx
        recty0 = self.starty
        rectx1 = self.canvas.canvasx(event.x)
        recty1 = self.canvas.canvasy(event.y)

        imgx0 = self.xborder
        imgy0 = self.yborder
        imgx1 = imgx0 + self.myImage.width()  - 1
        imgy1 = imgy0 + self.myImage.height() - 1
        
        #--order top corner, bottom corner--#
        if rectx1 < rectx0: rectx1, rectx0 = rectx0, rectx1 
        if recty1 < recty0: recty1, recty0 = recty0, recty1

        #--fit rectangel in image--#
        self.canvas.delete(self.rubberband)
        
        if (rectx0 < imgx1) and (recty0 < imgy1) and \
           (rectx1 > imgx0) and (recty1 > imgy0) :
            if (rectx0 < imgx0): rectx0 = imgx0
            if (recty0 < imgy0): recty0 = imgy0
            if (rectx1 > imgx1): rectx1 = imgx1
            if (recty1 > imgy1): recty1 = imgy1
            self.rubberband = self.canvas.create_rectangle(
                                rectx0, recty0, rectx1, recty1,
                                outline="red",width=1)
            self.zoom(rectx0-imgx0, recty0-imgy0, rectx1-imgx0, recty1-imgy0)
        else:        
            pass


    def _showValue(self,event):
        pass
    def _hideValue(self,event):
        pass
    def _setsize(self):
        pass
