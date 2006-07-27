
/////////////////////////////////////////////////////////////////////////////
//
// File:     camfr_wrap_2.cpp
// Author:   Peter.Bienstman@UGent.be
//
// Copyright (C) 2002-2006 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <boost/python.hpp>

#include "camfr_wrap.h"
#include "cavity.h"
#include "bloch.h"
#include "infstack.h"
#include "primitives/planar/planar.h"
#include "primitives/circ/circ.h"
#include "primitives/slab/generalslab.h"
#include "primitives/slab/isoslab/slab.h"
#include "primitives/slab/isoslab/slabwall.h"
#include "primitives/slab/isoslab/slabdisp.h"
#include "primitives/section/section.h"
#include "primitives/section/sectiondisp.h"
#include "primitives/section/refsection.h"
#include "primitives/blochsection/blochsection.h"
#include "primitives/blochsection/blochsectionmode.h"

/////////////////////////////////////////////////////////////////////////////
//
// Python wrappers for the CAMFR classes.
// The wrappers are created with the Boost library.
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// Additional functions used in the Python interface.
//
/////////////////////////////////////////////////////////////////////////////

inline Complex planar_static_get_kt(Planar p) {return Planar::get_kt();}
inline void planar_static_set_kt(Planar p, Complex kt) {Planar::set_kt(kt);}

inline SectionMode* section_get_mode(const Section& s, int i)
{
  check_wg_index(s,i); 
  return dynamic_cast<SectionMode*>(s.get_mode(i+1));
}

inline BlochSectionMode* blochsection_get_mode(const BlochSection& s, int i)
{
  check_wg_index(s,i);
  return dynamic_cast<BlochSectionMode*>(s.get_mode(i+1));
}

inline BlochMode* blochstack_get_mode(BlochStack& b, int i)
{
  check_wg_index(b,i); 
  return dynamic_cast<BlochMode*>(b.get_mode(i+1));
}

inline void cavity_set_current_source(Cavity& c, Coord& pos, Coord& ori) 
  {c.set_source(pos,ori);}

inline void cavity_set_general_source(Cavity& c,
                                      const cVector& fw, const cVector& bw) 
  {c.set_source(fw,bw);}

inline boost::python::object blochmode_fw_bw(BlochMode& b, Real z)
{
  cVector fw(global.N,fortranArray);
  cVector bw(global.N,fortranArray);

  b.fw_bw_field(Coord(0,0,z), &fw, &bw);

  return boost::python::make_tuple(fw, bw);
}

inline boost::python::object blochmode_fw_bw_2(BlochMode& b, Real z, Limit l)
{
  cVector fw(global.N,fortranArray);
  cVector bw(global.N,fortranArray);

  b.fw_bw_field(Coord(0,0,z,Plus,Plus,l), &fw, &bw);

  return boost::python::make_tuple(fw, bw);
}

inline Real blochstack_length(BlochStack& bs) 
  {return real(bs.get_total_thickness());} 
inline Real blochstack_width(BlochStack& bs)
  {return real(bs.c1_size());}
inline Real cavity_length(Cavity& c) 
  {return real(c.get_bot()->get_total_thickness() 
             + c.get_top()->get_total_thickness());} 
inline Real cavity_width(Cavity& c)
  {return real(c.get_top()->get_inc()->c1_size());}
inline Real slab_width(Slab& s)
  {return real(s.get_width());}
inline Real section_width(Section& s)
  {return real(s.get_width());}
inline Real section_height(Section& s)
  {return real(s.get_height());}
inline Real blochsection_width(BlochSection& s)
  {return real(s.get_width());}
inline Real blochsection_height(BlochSection& s)
  {return real(s.get_height());}
inline Complex blochmode_n(BlochMode& m, Coord &c)
  {return m.get_geom()->n_at(c);}
inline Complex sectionmode_n(SectionMode& m, Coord &c)
  {return m.get_geom()->n_at(c);}


  
/////////////////////////////////////////////////////////////////////////////
//
// Functions dealing with handling of default arguments.
//
/////////////////////////////////////////////////////////////////////////////

Real cavity_calc_sigma(Cavity& c)
  {return c.calc_sigma();}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cav_find_mode, Cavity::find_mode,2,5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cav_find_modes, \
  Cavity::find_modes_in_region,3,7)



/////////////////////////////////////////////////////////////////////////////
//
// The following functions are used when expanding an abritrarily shaped 
// field in slabmodes.
//
/////////////////////////////////////////////////////////////////////////////

inline cVector slab_expand_field(Slab& s, PyObject* o, Real eps)
  {PythonFunction f(o); return s.expand_field(&f, eps);}

inline cVector slab_expand_gaussian
  (Slab& s, Complex height, Complex width, Complex pos, Real eps)
    {GaussianFunction f(height,width,pos); return s.expand_field(&f, eps);}

inline cVector slab_expand_plane_wave
  (Slab& s, const Complex& amplitude, const Complex& angle, Real eps)
{
  Complex index = s.get_core()->n();
  PlaneWaveFunction f(amplitude,angle, index);
  return s.expand_field(&f, eps);
}



/////////////////////////////////////////////////////////////////////////////
//
// More exported functions.
//
/////////////////////////////////////////////////////////////////////////////

void camfr_wrap_2()
{
  using namespace boost::python;

  // Wrap Cavity.

  class_<Cavity>("Cavity", init<Stack&, Stack&>())
    .def("find_mode",      &Cavity::find_mode, cav_find_mode())
    .def("find_all_modes", &Cavity::find_modes_in_region,cav_find_modes())
    .def("sigma",          cavity_calc_sigma)
    .def("set_source",     cavity_set_current_source)
    .def("set_source",     cavity_set_general_source)
    .def("length",         cavity_length)
    .def("width",          cavity_width)
    .def("field",          &Cavity::field)
    .def("n",              &Cavity::n_at)
    .def("bot_stack",      &Cavity::get_bot,
         return_value_policy<reference_existing_object>())
    .def("top_stack",      &Cavity::get_top,
         return_value_policy<reference_existing_object>())
    ;

  // Wrap BlochStack.

  class_<BlochStack, bases<MultiWaveguide> >
    ("BlochStack", init<const Expression&>())
    .def("mode",        blochstack_get_mode,
         return_value_policy<reference_existing_object>())
    .def("length",      blochstack_length)
    .def("width",       blochstack_width)
    .def("beta_vector", &BlochStack::get_beta_vector)
    .def("__repr__",    &BlochStack::repr)
    ;

  // Wrap BlochMode.

  class_<BlochMode, bases<Mode> >("BlochMode", no_init)
    .def("fw_field", &BlochMode::fw_field)
    .def("bw_field", &BlochMode::bw_field)
    .def("fw_bw",    blochmode_fw_bw)
    .def("fw_bw",    blochmode_fw_bw_2)
    .def("S_flux",   &BlochMode::S_flux)
    .def("n",        blochmode_n)
    ;

  // Wrap InfStack.

  class_<InfStack, bases<DenseScatterer> >
    ("InfStack", init<const Expression&>())
    .def("R12", &InfStack::get_R12,
         return_value_policy<reference_existing_object>())
    ;

  // Wrap RealFunction.

  class_<RealFunction, boost::noncopyable>("RealFunction", no_init)
    .def("times_called", &RealFunction::times_called)
    .def("__call__",     &RealFunction::operator())
    ;

  // Wrap ComplexFunction.

  class_<ComplexFunction, boost::noncopyable>("RealFunction", no_init)
    .def("times_called", &ComplexFunction::times_called)
    .def("__call__",     &ComplexFunction::operator())
    ;

  // Wrap Planar.

  class_<Planar, bases<MonoWaveguide> >("Planar", init<Material&>())
    .def("set_theta", &Planar::set_theta)    
    .def("set_kt",    planar_static_set_kt)
    .def("get_kt",    planar_static_get_kt)
    ;

  // Wrap Circ.

  class_<Circ, bases<MultiWaveguide> >("Circ", init<Term&>())
    .def(init<Expression&>())
    ;

  // Wrap SlabWall.

  class_<SlabWall, boost::noncopyable>("SlabWall", no_init)
    .def("R", &SlabWall::get_R12)
    ;

  // Wrap SlabWallMixed.

  class_<SlabWallMixed, bases<SlabWall> >
    ("SlabWallMixed", init<const Complex&, const Complex&>());

  scope().attr("slab_E_wall")  = SlabWallMixed(1.0,  1.0);
  scope().attr("slab_H_wall")  = SlabWallMixed(1.0, -1.0);
  scope().attr("slab_no_wall") = SlabWallMixed(0.0,  1.0);

  // Wrap SlabWall_TBC.

  class_<SlabWall_TBC, bases<SlabWall> >
    ("SlabWall_TBC", init<const Complex&, const Material&>());

  // Wrap SlabWall_PC.

  class_<SlabWall_PC, bases<SlabWall> >
    ("SlabWall_PC", init<const Expression&>());

  // Wrap SlabDisp.

  class_<SlabDisp, bases<ComplexFunction> >
    ("SlabDisp", init<Expression&, Real>())
    .def(init<Expression&, Real, SlabWall*, SlabWall*>())
    ;

  // Wrap Slab.

  class_<Slab, bases<MultiWaveguide> >("Slab", init<const Term&>())
    .def(init<const Expression&>())
    .def("set_lower_wall",    &Slab::set_lower_wall)
    .def("set_upper_wall",    &Slab::set_upper_wall)
    .def("width",             slab_width)
    //.def("disp",              &Slab::get_disp)
    .def("expand_field",      slab_expand_field)
    .def("expand_gaussian",   slab_expand_gaussian) 
    .def("expand_plane_wave", slab_expand_plane_wave)
    .def("set_dummy",         &Slab::set_dummy)    
    .def("add_kz2_estimate",  &Slab::add_kz2_estimate)
    ;

  // Wrap SectionDisp.

  class_<SectionDisp, bases<ComplexFunction> >
    ("SectionDisp", init<Stack&, Stack&, Real, int>());

  // Wrap Section.

  class_<Section, bases<MultiWaveguide> >
  ("Section", init<Expression&, optional<int, int> >())
    .def(init<Expression&, Expression&, optional<int, int> >())    
    .def(init<const Term&>())
    .def("mode",         section_get_mode,
         return_value_policy<reference_existing_object>())
    .def("disp",         &Section::get_disp)
    .def("width",        section_width)
    .def("height",       section_height)
    .def("eps",          &Section::eps_at)
    .def("mu",           &Section::mu_at)
    .def("n",            &Section::n_at)
    .def("set_sorting",  &Section::set_sorting)
    .def("set_estimate", &Section::set_estimate)
    ;

  // Wrap RefSection.

  class_<RefSection, bases<MultiWaveguide> >
  ("RefSection", init<Material&, const Complex&, const Complex&, int>())
    ;  

  // Wrap SectionMode.

  class_<SectionMode, boost::noncopyable, bases<Mode> >
    ("SectionMode", no_init)
    .def("n", sectionmode_n)
    ;

  // Wrap BlochSection.

  class_<BlochSection, bases<MultiWaveguide> >
  ("BlochSection", init<Expression& >())    
    .def(init<const Term&>())
    .def("mode",          blochsection_get_mode,
         return_value_policy<reference_existing_object>())
    .def("width",         blochsection_width)
    .def("height",        blochsection_height)
    .def("eps",           &BlochSection::eps_at)
    .def("mu",            &BlochSection::mu_at)
    .def("n",             &BlochSection::n_at)
    .def("order",         &BlochSection::order)
    .def("set_theta_phi", &BlochSection::set_theta_phi)    
    .def("set_kx0_ky0",   &BlochSection::set_kx0_ky0)    
    .def("get_kx0",       &BlochSection::get_kx0)    
    .def("get_ky0",       &BlochSection::get_ky0)
    ;

  // Wrap BlochSectionMode.

  class_<BlochSectionMode, boost::noncopyable, bases<Mode> >
    ("BlochSectionMode", no_init)
    .def("get_Mx",   &BlochSectionMode::get_Mx)    
    .def("get_My",   &BlochSectionMode::get_My)    
    .def("get_kx",   &BlochSectionMode::get_kx)    
    .def("get_ky",   &BlochSectionMode::get_ky)
    ;

}
