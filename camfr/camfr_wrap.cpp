
/////////////////////////////////////////////////////////////////////////////
//
// File:          camfr_wrap.cpp
// Author:        Peter.Bienstman@rug.ac.be
// Date:          20021119
// Version:       2.1
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <boost/python.hpp>
#include "Numeric/arrayobject.h"

#include "defs.h"
#include "coord.h"
#include "mode.h"
#include "field.h"
#include "material.h"
#include "waveguide.h"
#include "scatterer.h" 
#include "expression.h"
#include "stack.h"
#include "cavity.h"
#include "bloch.h"
#include "icache.h"
#include "infstack.h"
#include "math/calculus/function.h"
#include "primitives/planar/planar.h"
#include "primitives/circ/circ.h"
#include "primitives/slab/generalslab.h"
#include "primitives/slab/isoslab/slab.h"
#include "primitives/slab/isoslab/slabwall.h"
#include "primitives/slab/isoslab/slabdisp.h"
#include "primitives/section/section.h"
#include "primitives/section/sectiondisp.h"
#include "primitives/section/refsection.h"

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

inline void set_lambda(Complex l)
  {global.lambda = l;}

inline Complex get_lambda()
  {return global.lambda;}

inline void set_N(int n)
  {global.N = n;}

inline int get_N()
  {return global.N;}

inline void set_polarisation(Polarisation pol)
  {global.polarisation = pol;}

inline Polarisation get_polarisation()
  {return global.polarisation;}

inline void set_gain_material(Material* m)
  {global.gain_mat = m;} 

inline void set_solver(Solver s)
  {global.solver = s;} 

inline void set_stability(Stability s)
  {global.stability = s;}

inline void set_precision(unsigned int i)
  {global.precision = i;} 

inline void set_precision_enhancement(unsigned int i)
  {global.precision_enhancement = i;}

inline void set_dx_enhanced(Real d)
  {global.dx_enhanced = d;} 

inline void set_precision_rad(unsigned int i)
  {global.precision_rad = i;} 

inline void set_C_upperright(const Complex& c)
  {global.C_upperright = c;} 

inline void set_sweep_from_previous(bool b)
  {global.sweep_from_previous = b;} 

inline void set_sweep_steps(unsigned int i)
  {global.sweep_steps = i;} 

inline void set_eps_trace_coarse(Real d)
  {global.eps_trace_coarse = d;} 

inline void set_chunk_tracing(bool b)
  {global.chunk_tracing = b;} 

inline void set_unstable_exp_threshold(Real d)
  {global.unstable_exp_threshold = d;} 

inline void set_field_calc_heuristic(Field_calc_heuristic f)
  {global.field_calc_heuristic = f;} 

inline void set_bloch_calc(Bloch_calc s)
  {global.bloch_calc = s;}

inline void set_eigen_calc(Eigen_calc s)
  {global.eigen_calc = s;}

inline void set_orthogonal(bool b)
  {global.orthogonal = b;}

inline void set_degenerate(bool b)
  {global.degenerate = b;}

inline void set_circ_order(int n)
  {global_circ.order = n;}

inline void set_circ_fieldtype(Fieldtype f)
  {global_circ.fieldtype = f;}

inline void set_lower_PML(Real PML)
{
  if (PML > 0)
    py_print("Warning: gain in PML.");

  global_slab.lower_PML = PML;
}

inline void set_upper_PML(Real PML)
{
  if (PML > 0)
    py_print("Warning: gain in PML.");

  global_slab.upper_PML = PML;
}

inline void set_left_PML(Real PML)
{
  if (PML > 0)
    py_print("Warning: gain in PML.");

  global_section.left_PML = PML;
}

inline void set_right_PML(Real PML)
{
  if (PML > 0)
    py_print("Warning: gain in PML.");

  global_section.right_PML = PML;
}

inline void set_circ_PML(Real PML)
{
  if (PML > 0)
    py_print("Warning: gain in PML.");

  global_circ.PML = PML;
}

inline void set_left_wall(Section_wall_type w)
  {global_section.leftwall = w;}

inline void set_right_wall(Section_wall_type w)
  {global_section.rightwall = w;}

inline void set_lower_wall(SlabWall* w)
  {global_slab.lowerwall = w;}

inline void set_upper_wall(SlabWall* w)
  {global_slab.upperwall = w;}

inline void set_beta(const Complex& beta)
  {global.slab_ky = beta;}

inline void set_mode_surplus(Real l)
  {global.mode_surplus = l;}

inline void set_backward_modes(bool b)
  {global.backward_modes = b;}

inline void set_section_solver(Section_solver s)
  {global_section.section_solver = s;}

inline void set_mode_correction(Mode_correction c)
  {global_section.mode_correction = c;}

inline void set_keep_all_estimates(bool b)
  {global_section.keep_all_estimates = b;}

inline void set_eta_ASR(Real eta)
{
  if ((eta > 1.0) || (eta < 0.0))
    py_print("Warning: eta_ASR should be between 0 and 1.");

  global_slab.eta_ASR=eta;
}

inline void set_davy(bool b)
  {global.davy = b;}

inline int mode_pol(const Mode& m) {return m.pol;}

inline Complex field_E1(const Field& f) {return f.E1;}
inline Complex field_E2(const Field& f) {return f.E2;}
inline Complex field_Ez(const Field& f) {return f.Ez;}
inline Complex field_H1(const Field& f) {return f.H1;}
inline Complex field_H2(const Field& f) {return f.H2;}
inline Complex field_Hz(const Field& f) {return f.Hz;}

inline void check_index(int i)
{
  if ( (i<0) || (i>=int(global.N)) )
  {
    PyErr_SetString(PyExc_IndexError, "index out of bounds.");
    throw std::out_of_range("index out of bounds.");
  }
}

inline void check_wg_index(const Waveguide& w, int i)
{
  if ( (i<0) || (i>=int(w.N())) )
  {
    PyErr_SetString(PyExc_IndexError, "index out of bounds.");
    throw std::out_of_range("index out of bounds.");
  }
}

inline Mode* waveguide_get_mode(const Waveguide& w, int i)
  {check_wg_index(w,i); return w.get_mode(i+1);}
inline Mode* waveguide_get_fw_mode(const Waveguide& w, int i)
  {check_wg_index(w,i); return w.get_fw_mode(i+1);}
inline Mode* waveguide_get_bw_mode(const Waveguide& w, int i)
  {check_wg_index(w,i); return w.get_bw_mode(i+1);}

inline SectionMode* section_get_mode(const Section& s, int i)
  {check_wg_index(s,i); return dynamic_cast<SectionMode*>(s.get_mode(i+1));}

inline BlochMode* blochstack_get_mode(BlochStack& b, int i)
  {check_wg_index(b,i); return dynamic_cast<BlochMode*>(b.get_mode(i+1));}

inline Complex stack_R12(const Stack& s, int i, int j)
  {check_index(i); check_index(j); return s.R12(i+1,j+1);}
inline Complex stack_R21(const Stack& s, int i, int j)
  {check_index(i); check_index(j); return s.R21(i+1,j+1);}
inline Complex stack_T12(const Stack& s, int i, int j)
  {check_index(i); check_index(j); return s.T12(i+1,j+1);}
inline Complex stack_T21(const Stack& s, int i, int j)
  {check_index(i); check_index(j); return s.T21(i+1,j+1);}

inline Real stack_inc_S_flux(Stack& s, Real c1_start, Real c1_stop, Real eps)
  {return dynamic_cast<MultiWaveguide*>(s.get_inc())
     ->S_flux(s.inc_field_expansion(),c1_start,c1_stop,eps);}

inline Real stack_ext_S_flux(Stack& s, Real c1_start, Real c1_stop, Real eps)
  {return dynamic_cast<MultiWaveguide*>(s.get_ext())
     ->S_flux(s.ext_field_expansion(),c1_start,c1_stop,eps);}

inline void cavity_set_current_source(Cavity& c, Coord& pos, Coord& ori) 
  {c.set_source(pos,ori);}

inline void cavity_set_general_source(Cavity& c,
                                      const cVector& fw, const cVector& bw) 
  {c.set_source(fw,bw);}

inline boost::python::object stack_fw_bw(Stack& s, Real z)
{
  cVector fw(global.N,fortranArray);
  cVector bw(global.N,fortranArray);

  s.fw_bw_field(Coord(0,0,z), &fw, &bw);

  return boost::python::make_tuple(fw, bw);
}

inline boost::python::object stack_fw_bw_2(Stack& s, Real z, Limit l)
{
  cVector fw(global.N,fortranArray);
  cVector bw(global.N,fortranArray);

  s.fw_bw_field(Coord(0,0,z,Plus,Plus,l), &fw, &bw);

  return boost::python::make_tuple(fw, bw);
}

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

inline Real stack_length(Stack& s) 
  {return real(s.get_total_thickness());}
inline Real stack_width(Stack& s) 
  {return real(s.get_inc()->c1_size());}
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
inline Complex blochmode_n(BlochMode& m, Coord &c)
  {return m.get_geom()->n_at(c);}
inline Complex sectionmode_n(SectionMode& m, Coord &c)
  {return m.get_geom()->n_at(c);}

inline void free_tmp_interfaces(Waveguide& w)
  {interface_cache.deregister(&w);}



/////////////////////////////////////////////////////////////////////////////
//
// The following functions and classes are used when expanding an
// abritrarily shaped field in slabmodes.
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// PythonFunction
//
/////////////////////////////////////////////////////////////////////////////

class PythonFunction : public ComplexFunction
{
  public:

    PythonFunction(PyObject* f_): f(f_) {}

    Complex operator()(const Complex& z)
      {counter++; return boost::python::call<Complex>(f, z);}

  protected:

    PyObject* f;
};

inline cVector slab_expand_field(Slab& s, PyObject* o, Real eps)
  {PythonFunction f(o); return s.expand_field(&f, eps);}

inline void stack_set_inc_field_function(Stack& s, PyObject* o, Real eps)
{
  Slab* slab = dynamic_cast<Slab*>(s.get_inc());
  
  if (!slab)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "set_inc_field_function only implemented for slabs.");
    exit(-1); //throw boost::python::argument_error();
  }

  PythonFunction f(o);
  s.set_inc_field(slab->expand_field(&f, eps));
}



/////////////////////////////////////////////////////////////////////////////
//
// GaussianFunction
//
/////////////////////////////////////////////////////////////////////////////

class GaussianFunction : public ComplexFunction
{
  public:

    GaussianFunction (Complex height, Complex width, Complex position)
      : h(height), w(width), p(position) {}

    Complex operator()(const Complex& x)
	{counter++; return h*exp(-(x-p)*(x-p)/(w*w*2.0));}

  protected:

    Complex h, w, p;
};

inline cVector slab_expand_gaussian
  (Slab& s, Complex height, Complex width, Complex pos, Real eps)
    {GaussianFunction f(height,width,pos); return s.expand_field(&f, eps);}

inline void stack_set_inc_field_gaussian
  (Stack& s, Complex height, Complex width, Complex pos, Real eps)
{
  Slab* slab = dynamic_cast<Slab*>(s.get_inc());
  
  if (!slab)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "set_inc_field_gaussian only implemented for slabs.");
    exit(-1); //throw boost::python::argument_error();
  }

  GaussianFunction f(height,width,pos);
  s.set_inc_field(slab->expand_field(&f, eps));
}



/////////////////////////////////////////////////////////////////////////////
//
// PlaneWaveFunction
//
/////////////////////////////////////////////////////////////////////////////

class PlaneWaveFunction : public ComplexFunction
{
  public:

    PlaneWaveFunction (Complex amplitude, Complex angle, Complex index)
      : am(amplitude), an(angle), n(index) {}

    Complex operator()(const Complex& x)
      {counter++; return am*exp(-2.0*I*pi*sin(an)*n*x/global.lambda);}

  protected:

    Complex am, an, n;
};

inline cVector slab_expand_plane_wave
  (Slab& s, const Complex& amplitude, const Complex& angle, Real eps)
{
        Complex index = s.get_core()->n();
        PlaneWaveFunction f(amplitude,angle, index); 
        return s.expand_field(&f, eps);
}

inline void stack_set_inc_field_plane_wave
  (Stack& s, const Complex& amplitude, const Complex& angle, Real eps)
{
  Slab* slab = dynamic_cast<Slab*>(s.get_inc());
 
  if (!slab)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "set_inc_field_plane_Wave only implemented for slabs.");
    exit (-1); //throw boost::python::argument_error();
  }

  Complex index = slab->get_core()->n();
  PlaneWaveFunction f(amplitude,angle,index);
  s.set_inc_field(slab->expand_field(&f, eps));
}



/////////////////////////////////////////////////////////////////////////////
//
// Functions converting C++ objects to and from Python objects.
//
/////////////////////////////////////////////////////////////////////////////

struct cVector_to_python
{
  static PyObject* convert(const cVector& c)
  {    
    int dim[1]; dim[0] = c.rows();

    PyArrayObject* result 
      = (PyArrayObject*) PyArray_FromDims(1, dim, PyArray_CDOUBLE);

    for (int i=0; i<c.rows(); i++)
        *(Complex*)(result->data + i*result->strides[0]) = c(i+1);

    return PyArray_Return(result);
  }
};

struct register_cVector_from_python
{

  register_cVector_from_python()
  {
    boost::python::converter::registry::insert
      (&convertible, &construct, boost::python::type_id<cVector>());
  }

  static void* convertible(PyObject* o)
  {
    if (!PyArray_Check(o))
      return NULL;

    if (    ( ((PyArrayObject*)(o))->nd != 1) 
         || ( ((PyArrayObject*)(o))->dimensions[0] != int(global.N)) )
      return NULL;

    return o;
  }
    
    
  static void construct
    (PyObject* o, boost::python::converter::rvalue_from_python_stage1_data* 
     data)
  {
    void* storage = ((
      boost::python::converter::rvalue_from_python_storage<cVector>*)data)
        ->storage.bytes;

    new (storage) cVector(global.N, fortranArray);

    PyArrayObject* a = (PyArrayObject *)
      PyArray_ContiguousFromObject(o, PyArray_CDOUBLE, 1, 1);

    for (int i=0; i<global.N; i++)      
      (*(cVector*)(storage))(i+1) = *(Complex*)(a->data + i*a->strides[0]);
    
    Py_DECREF(a);

    data->convertible = storage;
  }

};

struct cMatrix_to_python
{
  static PyObject* convert(const cMatrix& c)
  {
    int dim[2]; dim[0] = c.rows(); dim[1] = c.columns();
  
    PyArrayObject* result 
      = (PyArrayObject*) PyArray_FromDims(2, dim, PyArray_CDOUBLE);
    
    for (int i=0; i<c.rows(); i++)
      for (int j=0; j<c.columns(); j++)
        *(Complex*)(result->data + i*result->strides[0] + j*result->strides[1])
          = c(i+1,j+1);

    return PyArray_Return(result);
  }
};



/////////////////////////////////////////////////////////////////////////////
//
// Wrapper functions warning about deprecated features.
//
/////////////////////////////////////////////////////////////////////////////

Term material_to_term(BaseMaterial& m, const Complex& d)
{
  if (real(d) < 0)
    py_print("Warning: negative real length of material.");

  if(abs(imag(d)) > 0)
    py_error("Error: complex thickness deprecated in CAMFR 1.0.");
  
  return Term(m(d));
} 

Term waveguide_to_term(Waveguide& w, const Complex& d)
{
  if (real(d) < 0)
    py_print("Warning: negative real length of waveguide.");

  if(abs(imag(d)) > 0)
    py_error("Error: complex thickness deprecated in CAMFR 1.0.");
  
  return Term(w(d));
}


  
/////////////////////////////////////////////////////////////////////////////
//
// Functions dealing with handling of default arguments.
//
/////////////////////////////////////////////////////////////////////////////

Complex stack_lateral_S_flux(Stack& s, Real c)
  {return s.lateral_S_flux(c);}

Complex stack_lateral_S_flux_2(Stack& s, Real c, int k)
  {std::vector<Complex> S_k; s.lateral_S_flux(c, &S_k); return S_k[k];}

void stack_set_inc_field(Stack& s, const cVector& f)
  {s.set_inc_field(f);}

void stack_set_inc_field_2(Stack& s, const cVector& f, const cVector& b) 
  {s.set_inc_field(f, &const_cast<cVector&>(b));}

Real cavity_calc_sigma(Cavity& c)
  {return c.calc_sigma();}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cav_find_mode, Cavity::find_mode,2,5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cav_find_modes, \
  Cavity::find_modes_in_region,3,7)



/////////////////////////////////////////////////////////////////////////////
//
// The CAMFR module itself
//
/////////////////////////////////////////////////////////////////////////////

#undef PyArray_Type
extern PyTypeObject PyArray_Type;

BOOST_PYTHON_MODULE(_camfr)
{
  using namespace boost::python;

  implicitly_convertible<Scatterer,Term>();
  implicitly_convertible<Stack,Term>();

  import_array();

  to_python_converter<cVector, cVector_to_python>();
  to_python_converter<cMatrix, cMatrix_to_python>();
  register_cVector_from_python();

  // Wrap Limit enum.

  enum_<Limit>("Limit")
    .value("Plus", Plus)
    .value("Min",  Min)
    ;

  scope().attr("Plus") = Plus;
  scope().attr("Min")  = Min;

  // Wrap Solver enum.

  enum_<Solver>("Solver")
    .value("ADR",    ADR)
    .value("track",  track)
    .value("series", series)    
    .value("ASR",    ASR)
    ;

  scope().attr("ADR")    = ADR;
  scope().attr("track")  = track;
  scope().attr("series") = series;
  scope().attr("ASR")    = ASR;

  // Wrap Stability enum.

  enum_<Stability>("Stability")
    .value("normal", normal)
    .value("extra",  extra)
    .value("SVD",    SVD)
    ;

  scope().attr("normal") = normal;
  scope().attr("extra")  = extra;
  scope().attr("SVD")    = SVD;

  // Wrap Field_calc_heuristic enum.

  enum_<Field_calc_heuristic>("Field_calc_heuristic")
    .value("identical", identical)
    .value("symmetric", symmetric)
    ;

  scope().attr("identical") = identical;
  scope().attr("symmetric") = symmetric;

  // Wrap Bloch_calc enum.

  enum_<Bloch_calc>("Bloch_calc")
    .value("GEV", GEV)
    .value("T",   T)
    ;

  scope().attr("GEV") = GEV;
  scope().attr("T")   = T;

  // Wrap Eigen_calc enum.

  enum_<Eigen_calc>("Eigen_calc")
    .value("lapack",  lapack)
    .value("arnoldi", arnoldi)
    ;

  scope().attr("lapack")  = lapack;
  scope().attr("arnoldi") = arnoldi;

  // Wrap Polarisation enum.

  enum_<Polarisation>("Polarisation")
    .value("unknown", unknown)
    .value("TEM",     TEM)
    .value("TE",      TE)
    .value("TM",      TM)
    .value("HE",      HE)
    .value("EH",      EH)
    .value("TE_TM",   TE_TM)
    ;

  scope().attr("unknown") = unknown;
  scope().attr("TEM")     = TEM;
  scope().attr("TE")      = TE;
  scope().attr("TM")      = TM;
  scope().attr("HE")      = HE;
  scope().attr("EH")      = EH;
  scope().attr("TE_TM")   = TE_TM;

  // Wrap Fieldtype enum.

  enum_<Fieldtype>("Fieldtype")
    .value("cos_type", cos_type)
    .value("sin_type", sin_type)
    ;

  scope().attr("cos_type") = cos_type;
  scope().attr("sin_type") = sin_type; 

 // Wrap Section_wall_type enum.

  enum_<Section_wall_type>("Section_wall_type")
    .value("E_wall", E_wall)
    .value("H_wall", H_wall)
    ;

  scope().attr("E_wall") = E_wall;
  scope().attr("H_wall") = H_wall;

  // Wrap Sort_type.

  enum_<Sort_type>("Sort_type")
    .value("highest_index", highest_index)
    .value("lowest_loss",   lowest_loss)
    ;
  
  scope().attr("highest_index") = highest_index;
  scope().attr("lowest_loss")   = lowest_loss;

  // Wrap Section_solver enum.

  enum_<Section_solver>("Section_solver")
    .value("OS", OS)
    .value("NT", NT)
    .value("L",  L)
    ;

  scope().attr("OS") = OS;
  scope().attr("NT") = NT; 
  scope().attr("L")  = L;   

  // Wrap Mode_correction enum.

  enum_<Mode_correction>("Mode_correction")
    .value("none",        none)
    .value("snap",        snap)
    .value("guided_only", guided_only)
    .value("full",        full)
    ;

  scope().attr("none")        = none;
  scope().attr("snap")        = snap;
  scope().attr("guided_only") = guided_only;
  scope().attr("full")        = full;

  // Wrap getters and setters for global parameters.

  def("set_lambda",                 set_lambda);
  def("get_lambda",                 get_lambda);
  def("set_N",                      set_N);
  def("N",                          get_N);
  def("set_polarisation",           set_polarisation);  
  def("get_polarisation",           get_polarisation);
  def("set_gain_material",          set_gain_material);
  def("set_solver",                 set_solver);
  def("set_stability",              set_stability);
  def("set_precision",              set_precision);
  def("set_precision_enhancement",  set_precision_enhancement);
  def("set_dx_enhanced",            set_dx_enhanced);
  def("set_precision_rad",          set_precision_rad);
  def("set_C_upperright",           set_C_upperright);
  def("set_sweep_from_previous",    set_sweep_from_previous);
  def("set_sweep_steps",            set_sweep_steps);
  def("set_eps_trace_coarse",       set_eps_trace_coarse);
  def("set_chunk_tracing",          set_chunk_tracing);
  def("set_unstable_exp_threshold", set_unstable_exp_threshold);
  def("set_field_calc_heuristic",   set_field_calc_heuristic);
  def("set_bloch_calc",             set_bloch_calc);
  def("set_eigen_calc",             set_eigen_calc);
  def("set_orthogonal",             set_orthogonal);
  def("set_degenerate",             set_degenerate);
  def("set_circ_order",             set_circ_order);
  def("set_circ_field_type",        set_circ_fieldtype);
  def("set_left_wall",              set_left_wall);
  def("set_right_wall",             set_right_wall);
  def("set_upper_wall",             set_upper_wall);
  def("set_lower_wall",             set_lower_wall);
  def("set_left_PML",               set_left_PML);
  def("set_right_PML",              set_right_PML);
  def("set_upper_PML",              set_upper_PML);
  def("set_lower_PML",              set_lower_PML);
  def("set_circ_PML",               set_circ_PML);  
  def("set_eta_ASR",                set_eta_ASR);
  def("set_beta",                   set_beta);
  def("set_section_solver",         set_section_solver);    
  def("set_keep_all_estimates",     set_keep_all_estimates);  
  def("set_mode_correction",        set_mode_correction);
  def("set_mode_surplus",           set_mode_surplus);
  def("set_backward_modes",         set_backward_modes); 
  def("set_davy",                   set_davy);
  def("free_tmps",                  free_tmps);
  def("free_tmp_interfaces",        free_tmp_interfaces);

  // Wrap Coord.

  class_<Coord>("Coord", init<const Real&, const Real&, const Real&,
         optional<Limit, Limit, Limit> >())
    .def("__repr__", &Coord::repr)
    ;

  // Wrap Field.

  class_<Field>("Field", no_init)
    .def("E1",       field_E1)
    .def("E2",       field_E2)
    .def("Ez",       field_Ez)
    .def("H1",       field_H1)
    .def("H2",       field_H2)
    .def("Hz",       field_Hz)
    .def("S1",       &Field::S1)
    .def("S2",       &Field::S2)
    .def("Sz",       &Field::Sz)
    .def("abs_E",    &Field::abs_E)
    .def("abs_H",    &Field::abs_H)
    .def("abs_S",    &Field::abs_S)
    .def("__repr__", &Field::repr)
    ;

  // Wrap FieldExpansion.

  class_<FieldExpansion>("FieldExpansion", no_init)
    .def("field",    &FieldExpansion::field)
    .def("__repr__", &FieldExpansion::repr)
    ;

  // Wrap BaseMaterial.

  class_<BaseMaterial, boost::noncopyable>("BaseMaterial", no_init);

  // Wrap Material.

  class_<Material, bases<BaseMaterial> >
    ("Material", init<const Complex&, optional<const Complex&> >())
    .def("__call__", material_to_term)
    .def("n",        &Material::n)
    .def("set_n",    &Material::set_n)
    .def("epsr",     &Material::epsr)
    .def("mur",      &Material::mur)
    .def("set_mur",  &Material::set_mur)
    .def("eps",      &Material::eps)
    .def("mu",       &Material::mu)
    .def("gain",     &Material::gain)
    .def("__repr__", &Material::repr)
    ;

  // Wrap Material_length.

  class_<Material_length>("Material_length", no_init);

  // Wrap Mode.

  class_<Mode>("Mode", no_init)
    .def("field",    &Mode::field)
    .def("n_eff",    &Mode::n_eff)
    .def("kz",       &Mode::get_kz)
    .def("pol",      mode_pol)
    .def("__repr__", &Mode::repr)
    ;

  // Wrap Waveguide.

  class_<Waveguide, boost::noncopyable>("Waveguide", no_init)
    .def("core",     &Waveguide::get_core,
         return_value_policy<reference_existing_object>())
    .def("eps",      &Waveguide::eps_at)
    .def("mu",       &Waveguide::mu_at)
    .def("n",        &Waveguide::n_at)
    .def("N",        &Waveguide::N)
    .def("mode",     waveguide_get_mode,
         return_value_policy<reference_existing_object>())
    .def("fw_mode",  waveguide_get_fw_mode,
         return_value_policy<reference_existing_object>())
    .def("bw_mode",  waveguide_get_bw_mode,
         return_value_policy<reference_existing_object>())
    .def("calc",     &Waveguide::find_modes)
    .def("__repr__", &Waveguide::repr)
    .def("__call__", waveguide_to_term)
    ;

  // Wrap Waveguide_length.

  class_<Waveguide_length>("Waveguide_length", no_init);

  // Wrap MultiWaveguide.

  class_<MultiWaveguide, bases<Waveguide>, boost::noncopyable>
    ("MultiWaveguide", no_init)
    .def("field_from_source", &MultiWaveguide::field_from_source)
    ;

  // Wrap MonoWaveguide.

  class_<MonoWaveguide, bases<Waveguide>, boost::noncopyable>
    ("MonoWaveguide", no_init);

  // Wrap Scatterer.

  class_<Scatterer, boost::noncopyable>("Scatterer", no_init)
    .def("calc", &Scatterer::calcRT)
    .def("free", &Scatterer::freeRT)
    .def("inc",  &Scatterer::get_inc,
         return_value_policy<reference_existing_object>())
    .def("ext",  &Scatterer::get_ext,
         return_value_policy<reference_existing_object>())
    .def(self + Expression())
    .def(self + Term())
    ;

  // Wrap MultiScatterer.

  class_<MultiScatterer, bases<Scatterer>, boost::noncopyable>
    ("MultiScatterer", no_init);

  // Wrap DenseScatterer.

  class_<DenseScatterer, bases<MultiScatterer>, boost::noncopyable>
    ("DenseScatterer", no_init);

  // Wrap DiagScatterer.

  class_<DiagScatterer, bases<MultiScatterer>, boost::noncopyable>
    ("DiagScatterer", no_init);

  // Wrap MonoScatterer.

  class_<MonoScatterer, bases<Scatterer>, boost::noncopyable>
    ("MonoScatterer", no_init);

  // Wrap SquashedScatterer.

  class_<SquashedScatterer, bases<DenseScatterer> >
    ("SquashedScatterer", init<DenseScatterer&>());

  // Wrap FlippedScatterer.

  class_<FlippedScatterer, bases<MultiScatterer> >
    ("FlippedScatterer", init<MultiScatterer&>());

  // Wrap E_Wall.

  class_<E_Wall, bases<DiagScatterer> >("E_Wall", init<Waveguide&>());

  // Wrap H_Wall.

  class_<H_Wall, bases<DiagScatterer> >("H_Wall", init<Waveguide&>());

  // Wrap Expression.

  class_<Expression>("Expression")
    .def(init<const Term&>())
    .def(init<const Expression&>())
    .def("flatten",  &Expression::flatten)
    .def("inc",  &Expression::get_inc,
         return_value_policy<reference_existing_object>())
    .def("ext",  &Expression::get_ext,
         return_value_policy<reference_existing_object>())
    .def("__repr__", &Expression::repr)
    .def("add",      &Expression::operator+=)
    .def(self += self)
    .def(self + self)
    .def(self + Term())
    .def(self * int())
    .def(int() * self)
    ;

  // Wrap Term.

  class_<Term>("Term", init<Scatterer&>())
    .def(init<Stack&>())
    .def(init<const Expression&>())
    .def("inc",  &Term::get_inc,
         return_value_policy<reference_existing_object>())
    .def("ext",  &Term::get_ext,
         return_value_policy<reference_existing_object>())
    .def("__repr__", &Term::repr)
    .def(self + self)
    .def(self + Expression())
    .def(self * int())
    .def(int() * self)
    ;

  // Wrap Stack.

  class_<Stack>("Stack", init<const Expression&, optional<int> >())
    .def(init<const Term&, optional<int> >())
    .def("calc",                     &Stack::calcRT)
    .def("free",                     &Stack::freeRT)
    .def("inc",                      &Stack::get_inc,
         return_value_policy<reference_existing_object>())
    .def("ext",                      &Stack::get_ext,
         return_value_policy<reference_existing_object>())
    .def("scatterer",                &Stack::as_multi,
         return_value_policy<reference_existing_object>())
    .def("length",                   stack_length)
    .def("width",                    stack_width)
    .def("set_inc_field",            stack_set_inc_field)
    .def("set_inc_field",            stack_set_inc_field_2)
    .def("set_inc_field_function",   stack_set_inc_field_function)
    .def("set_inc_field_gaussian",   stack_set_inc_field_gaussian)
    .def("set_inc_field_plane_wave", stack_set_inc_field_plane_wave)
    .def("inc_field",                &Stack::get_inc_field)
    .def("refl_field",               &Stack::get_refl_field)
    .def("trans_field",              &Stack::get_trans_field)
    .def("inc_S_flux",               stack_inc_S_flux)
    .def("ext_S_flux",               stack_ext_S_flux)
    .def("field",                    &Stack::field)
    .def("fw_bw",                    stack_fw_bw)
    .def("fw_bw",                    stack_fw_bw_2)
    .def("lateral_S_flux",           stack_lateral_S_flux)
    .def("lateral_S_flux",           stack_lateral_S_flux_2)
    .def("eps",                      &Stack::eps_at)
    .def("mu",                       &Stack::mu_at)
    .def("n",                        &Stack::n_at)
    .def("R12",                      &Stack::get_R12)
    .def("R21",                      &Stack::get_R21)
    .def("T12",                      &Stack::get_T12)
    .def("T21",                      &Stack::get_T21)
    .def("R12",                      stack_R12)
    .def("R21",                      stack_R21)
    .def("T12",                      stack_T12)
    .def("T21",                      stack_T21)
    .def(self + Expression())
    .def(self + Term())
    ;

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
    .def("set_kt",    &Planar::set_kt)
    .def("get_kt",    &Planar::get_kt)
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

  scope().attr("slab_E_wall")    = SlabWallMixed(1.0,  1.0);
  scope().attr("slab_H_wall")    = SlabWallMixed(1.0, -1.0);
  scope().attr("slab_open_wall") = SlabWallMixed(0.0,  1.0);

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
}

