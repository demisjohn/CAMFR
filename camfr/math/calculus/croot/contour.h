
/////////////////////////////////////////////////////////////////////////////
//
// File:     contour.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20010322
// Version:  1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef CONTOUR_H
#define CONTOUR_H

#include "../function.h"

/////////////////////////////////////////////////////////////////////////////
//
// CLASS: Contour
//
//   A rectangular contour in the complex plane bounded by
//   (bottom_left,top_right). Used for calculation of the values of the
//   integrals z^n / f(z) for n=0..M at its egdes and in the four
//   subcontours.
//
//   Points on the contour/subcontours are labeled as follows:
//
//                  tl  tc  tr
//                  cl  cc  cr
//                  bl  bc  br
//
//   t,b: top, bottom, l,r: left, right, c:center
//
//   eps, mu, max_k are numeric parameters for the Patterson quadrature.
//
/////////////////////////////////////////////////////////////////////////////

typedef enum {top_left, top_right, bottom_left, bottom_right} Subcontour;
typedef enum {br_cr, cr_tr, tr_tc, tc_tl, tl_cl, cl_bl, bl_bc, bc_br,
              cc_tc, cc_cr, cc_bc, cc_cl} Segment;

class Contour
{
  public:

    Contour(const Complex& bottom_left, const Complex& top_right,
            ComplexFunction& f, unsigned int M,
            Real eps=1e-4, Real mu=1e-4, unsigned int max_k=8);
    Contour(const Contour&);
    Contour& operator=(const Contour&);
    ~Contour() {}

    bool encloses(const Complex& z) const;

    // Returns the top left, top right, bottom left or bottom right
    // subcontour, without precalculating any line segment.

    const Contour subcontour(Subcontour s) const;

    // Returns the contour integrals z^n/f, either over the entire contour,
    // or over a subcontour. Integrals along common line segments are reused.

    vector<Complex> contour_integrals() const
      {return path_integrals(8, main_segments, main_signs);}
    
    vector<Complex> subcontour_integrals(Subcontour s) const
      {return path_integrals(4, sub_segments[s], sub_signs[s]);}

    // Returns equally-sized adjacent contour to the right. Integrals along
    // common line segments are precalculated and reused.
    
    Contour adjacent_r() const;

    // Returns three equally-sized adjacent contours to the right, upwards,
    // and diagonally upwards-right. Integrals along common line segments
    // are precalculated and reused.
    
    vector<Contour> adjacent_ur() const;

    // Return a contour twice the size of the original one, with the same
    // bottom right point. Integrals along common line segments
    // are precalculated and reused.

    Contour double_ur() const;

    void set_integrals(Segment segment, const vector<Complex>& ints)
      {integrals[segment] = ints; know_integrals[segment] = true;}

    vector<Complex> get_integrals(Segment segment) const;

    ComplexFunction* get_f()  const {return f;}
    Complex get_bottom_left() const {return bl;}
    Complex get_top_right()   const {return tr;}
    Complex get_center()      const {return cc;}
    
  protected:

    unsigned int M;
    ComplexFunction* f;
    Real eps, mu;
    unsigned int max_k;

    Complex tl, tc, tr;
    Complex cl, cc, cr;
    Complex bl, bc, br;
    
    mutable bool know_integrals[12];
    mutable vector<Complex> integrals[12];

  private:

    void create_internal_points();
    void copy_from(const Contour& c);
    
    vector<Complex> path_integrals
      (unsigned int no_segments, Segment segments[], int signs[]) const;

    // Lookup tables describing contours, line segments,
    // begin and end points.

    Complex* begin[12];
    Complex*   end[12];
    
    static Segment main_segments[8], sub_segments[4][4];
    static int        main_signs[8], sub_signs   [4][4];
};



#endif
