
/////////////////////////////////////////////////////////////////////////////
//
// File:          camfr_wrap.cpp
// Author:        Peter.Bienstman@rug.ac.be
// Date:          20010705
// Version:       1.0
//
// Copyright (C) 2001 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string>
#include <iostream>
#include <stdexcept>
#include "boost/python/class_builder.hpp"
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

/////////////////////////////////////////////////////////////////////////////
//
// Python wrappers for the CAMFR classes.
// The wrappers are created with the Boost library.
//
/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace boost::python;



/////////////////////////////////////////////////////////////////////////////
//
// Additional functions used in the Python interface.
//
/////////////////////////////////////////////////////////////////////////////

inline void set_lambda(Real l)
  {global.lambda = l;}

inline Real get_lambda()
  {return global.lambda;}

inline void set_N(int n)
  {global.N = n;}

inline int get_N()
  {return global.N;}

inline void set_polarisation(long pol)
  {global.polarisation = Polarisation(pol);}

inline void set_gain_material(Material* m)
  {global.gain_mat = m;} 

inline void set_solver(long s)
  {global.solver = Solver(s);} 

inline void set_stability(long s)
  {global.stability = Stability(s);}

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

inline void set_field_calc(long f)
  {global.field_calc = Field_calc(f);} 

inline void set_bloch_calc(long s)
  {global.bloch_calc = Bloch_calc(s);}

inline void set_eigen_calc(long s)
  {global.eigen_calc = Eigen_calc(s);}

inline void set_orthogonal(bool b)
  {global.orthogonal = b;}

inline void set_circ_order(int n)
  {global_circ.order = n;}

inline void set_circ_fieldtype(long f)
  {global_circ.fieldtype = Fieldtype(f);}

inline void set_left_wall(SlabWall* w)
  {global_slab.leftwall = w;}

inline void set_right_wall(SlabWall* w)
  {global_slab.rightwall = w;}

inline void set_beta(const Complex& beta)
  {global_slab.beta = beta;}

template class enum_as_int_converters<Limit>;
template class enum_as_int_converters<Solver>;
template class enum_as_int_converters<Stability>;
template class enum_as_int_converters<Field_calc>;
template class enum_as_int_converters<Bloch_calc>;
template class enum_as_int_converters<Eigen_calc>;
template class enum_as_int_converters<Polarisation>;
template class enum_as_int_converters<Fieldtype>;

inline Complex mode_kz (const Mode& m) {return m.get_kz();}
inline Complex mode_pol(const Mode& m) {return m.pol;}

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

inline void check_bloch_index(int i)
{
  if ( (i<0) || (i>=int(2*global.N)) )
  {
    PyErr_SetString(PyExc_IndexError, "index out of bounds.");
    throw std::out_of_range("index out of bounds.");
  }
}

inline Mode* waveguide_get_mode(const Waveguide& w, int i)
  {check_index(i); return w.get_mode(i+1);}
inline Mode* waveguide_get_fw_mode(const Waveguide& w, int i)
  {check_index(i); return w.get_fw_mode(i+1);}
inline Mode* waveguide_get_bw_mode(const Waveguide& w, int i)
  {check_index(i); return w.get_bw_mode(i+1);}

inline BlochMode* blochstack_get_mode(BlochStack& b, int i)
  {check_bloch_index(i); return dynamic_cast<BlochMode*>(b.get_mode(i+1));}

inline cMatrix stack_get_R12(const Stack& s)
  {if (s.as_multi()) return s.as_multi()->get_R12();}
inline cMatrix stack_get_R21(const Stack& s)
  {if (s.as_multi()) return s.as_multi()->get_R21();}
inline cMatrix stack_get_T12(const Stack& s)
  {if (s.as_multi()) return s.as_multi()->get_T12();}
inline cMatrix stack_get_T21(const Stack& s)
  {if (s.as_multi()) return s.as_multi()->get_T21();}

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



/////////////////////////////////////////////////////////////////////////////
//
// Functions convertion C++ objects to and from Python objects.
//
// Note on memory management: as the Blitz ref-counting system might free
// data that is still in use in the Python environment, we have to make
// a copy of the array data. However, we do this only for vectors. For
// matrices, the problem does not occur in practice and the copying
// would incur more overhead.
//
/////////////////////////////////////////////////////////////////////////////

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

PyObject* to_python(cVector& c)
{
  char* copy = (char*) malloc(c.rows()*sizeof(Complex));
  memcpy(copy, c.data(), c.rows()*sizeof(Complex));                    
  
  int dim[1]; dim[0] = c.rows();
  
  PyArrayObject* result = (PyArrayObject*)
    PyArray_FromDimsAndData(1, dim, PyArray_CDOUBLE, copy);

  return PyArray_Return(result);
}

PyObject* to_python(const cVector& c)
  {return to_python(const_cast<cVector&>(c));}

cVector from_python(PyObject* o, type<const cVector&>)
{
  if (!PyArray_Check(o))
  {
    PyErr_SetString(PyExc_ValueError, "expected numerical array.");
    throw argument_error();
  }

  PyArrayObject* a = (PyArrayObject*) PyArray_Cast
    ((PyArrayObject*) o, PyArray_CDOUBLE);

  if ( (a->nd != 1) || (a->dimensions[0] != int(global.N)) )
  {
    PyErr_SetString(PyExc_ValueError, "array has wrong dimensions.");
    throw argument_error();
  }
  
  return cVector((Complex*) a->data, global.N, neverDeleteData, fortranArray);
}

PyObject* to_python(cMatrix& c)
{ 
  int dim[2]; dim[0] = global.N; dim[1] = global.N;

  cMatrix c_bak; c_bak = c.copy();
  
  PyArrayObject* result = (PyArrayObject*)
    PyArray_FromDimsAndData(2, dim, PyArray_CDOUBLE, (char*) c_bak.data());

  // Reflect Fortran storage orders.

  int tmp = result->strides[0];
  result->strides[0] = result->strides[1];
  result->strides[1] = tmp;

  return PyArray_Return(result);
}

PyObject* to_python(const cMatrix& c)
  {return to_python(const_cast<cMatrix&>(c));}

PyObject* to_python(Material* m)
  {return python_extension_class_converters<Material>::smart_ptr_to_python(m);}

PyObject* to_python(const Material* m)
  {return to_python(const_cast<Material*>(m));}

PyObject* to_python(Mode* m)
  {return python_extension_class_converters<Mode>::smart_ptr_to_python(m);}

PyObject* to_python(const Mode* m)
  {return to_python(const_cast<Mode*>(m));}

PyObject* to_python(BlochMode* m)
  {return python_extension_class_converters<BlochMode>
     ::smart_ptr_to_python(m);}

PyObject* to_python(const BlochMode* m)
  {return to_python(const_cast<BlochMode*>(m));}

PyObject* to_python(Waveguide* wg) {return
   python_extension_class_converters<Waveguide>::smart_ptr_to_python(wg);}

PyObject* to_python(const Waveguide* wg)
  {return to_python(const_cast<Waveguide*>(wg));}

BOOST_PYTHON_END_CONVERSION_NAMESPACE



/////////////////////////////////////////////////////////////////////////////
//
// These functions compensate for the lack of implicit type conversion
// in Python, mainly when converting to Term in expressions like
//
//   const Expression  operator+(const Term& L,       const Term& R);
//   const Expression& operator+(const Expression& L, const Term& R);
//
/////////////////////////////////////////////////////////////////////////////

Term material_to_term(BaseMaterial& m, const Complex& d)
{
  if (real(d) < 0)
    cout << "Warning: negative real length of material." << endl;
  
  return Term(m(d));
} 

Term waveguide_to_term(Waveguide& w, const Complex& d)
{
  if (real(d) < 0)
    cout << "Warning: negative real length of waveguide." << endl;
  
  return Term(w(d));
}

#define TERM_PLUS_TERM(t1, t2) inline const Expression \
 operator+(t1& L, t2& R) {return Term(L) + Term(R);}

#define TERM_PLUS_TERM_F(t1, t2) (const Expression (*) (t1&, t2&)) &operator+

TERM_PLUS_TERM(Stack,      Stack)
TERM_PLUS_TERM(Stack,      const Expression)
TERM_PLUS_TERM(Stack,      const Term)
TERM_PLUS_TERM(Stack,      Scatterer)
TERM_PLUS_TERM(const Term, Stack)
TERM_PLUS_TERM(const Term, const Expression)
TERM_PLUS_TERM(const Term, Scatterer)
TERM_PLUS_TERM(Scatterer,  Stack)
TERM_PLUS_TERM(Scatterer,  const Expression)
TERM_PLUS_TERM(Scatterer,  const Term)
TERM_PLUS_TERM(Scatterer,  Scatterer)
  
#define EX_PLUS_TERM(t) const Expression& \
 operator+(const Expression& L, t& R) {return L + Term(R);}

#define EX_PLUS_TERM_F(t) (const Expression& (*) (const Expression&, t&)) \
      &operator+
  
EX_PLUS_TERM(Stack)
EX_PLUS_TERM(const Expression)
EX_PLUS_TERM(Scatterer)


  
/////////////////////////////////////////////////////////////////////////////
//
// These functions compensate for the current lack of automatic handling
// of default arguments in Boost.
//
/////////////////////////////////////////////////////////////////////////////

Complex stack_lateral_S_flux(Stack& s, const Complex& c)
  {return s.lateral_S_flux(c);}

Complex stack_lateral_S_flux_2(Stack& s, const Complex& c, int k)
  {vector<Complex> S_k; s.lateral_S_flux(c, &S_k); return S_k[k];}

Real cavity_calc_sigma(Cavity& c)
  {return c.calc_sigma();}

void cavity_find_modes_in_region_3
  (Cavity& c, Real lambda_start, Real lambda_stop, Real delta_lambda)
    {return c.find_modes_in_region(lambda_start, lambda_stop, delta_lambda);} 

void cavity_find_modes_in_region_7
  (Cavity& c, Real lambda_start, Real lambda_stop,
   Real delta_lambda, unsigned int number,
   Real n_imag_start, Real n_imag_stop, unsigned int passes)
{
  return c.find_modes_in_region(lambda_start, lambda_stop, delta_lambda,
                                number, n_imag_start, n_imag_stop, passes);
}

void cavity_find_mode_2(Cavity& c, Real lambda_start, Real lambda_stop)
  {return c.find_mode(lambda_start, lambda_stop);}

void cavity_find_mode_5
  (Cavity& c, Real lambda_start, Real lambda_stop,
   Real n_imag_start, Real n_imag_stop, unsigned int passes)
{
  return
    c.find_mode(lambda_start, lambda_stop, n_imag_start, n_imag_stop, passes); 
}

void cavity_set_source(Cavity& c, Coord& pos, Coord& orientation)
  {c.set_source(pos,orientation);}

void stack_set_inc_field(Stack& s, const cVector& f) 
  {s.set_inc_field(f);}

void stack_set_inc_field_2(Stack& s, const cVector& f, const cVector& b) 
  {s.set_inc_field(f, &const_cast<cVector&>(b));}



/////////////////////////////////////////////////////////////////////////////
//
// The CAMFR module itself
//
/////////////////////////////////////////////////////////////////////////////

BOOST_PYTHON_MODULE_INIT(camfr_work)
{
  try
  {
     boost::python::module_builder camfr("camfr_work");

    import_array();
    
    // Splash screen.
    
    cout << endl
         << "CAMFR 1.0pre - "
         << "Copyright (C) 1998-2002 Peter Bienstman - Ghent University."
         << endl << endl;

    // Wrap Limit enum.
    
    camfr.add(boost::python::make_ref(Plus), "Plus");
    camfr.add(boost::python::make_ref(Min),  "Min");
    
    // Wrap Solver enum.
    
    camfr.add(boost::python::make_ref(ADR),   "ADR");
    camfr.add(boost::python::make_ref(track), "track");

    // Wrap Stability enum.
    
    camfr.add(boost::python::make_ref(normal), "normal");
    camfr.add(boost::python::make_ref(extra),  "extra");
    camfr.add(boost::python::make_ref(SVD),    "SVD");

    // Wrap Field_calc enum.
    
    camfr.add(boost::python::make_ref(T_T), "T_T");
    camfr.add(boost::python::make_ref(S_T), "S_T");
    camfr.add(boost::python::make_ref(S_S), "S_S");

    // Wrap Bloch_calc enum.

    camfr.add(boost::python::make_ref(GEV), "GEV");
    camfr.add(boost::python::make_ref(T),   "T");

    // Wrap Eigen_calc enum.

    camfr.add(boost::python::make_ref(lapack),  "lapack");
    camfr.add(boost::python::make_ref(arnoldi), "arnoldi");

    // Wrap Polarisation enum.
    
    camfr.add(boost::python::make_ref(unknown), "unknown");
    camfr.add(boost::python::make_ref(TEM),     "TEM");
    camfr.add(boost::python::make_ref(TE),      "TE");
    camfr.add(boost::python::make_ref(TM),      "TM");
    camfr.add(boost::python::make_ref(HE),      "HE");
    camfr.add(boost::python::make_ref(EH),      "EH");    
    camfr.add(boost::python::make_ref(TE_TM),   "TE_TM");

    // Wrap Fieldtype enum.

    camfr.add(boost::python::make_ref(cos_type), "cos_type");
    camfr.add(boost::python::make_ref(sin_type), "sin_type");    

    // Wrap getters and setters for global parameters.
      
    camfr.def(set_lambda,                 "set_lambda");
    camfr.def(get_lambda,                 "get_lambda");
    camfr.def(set_N,                      "set_N");
    camfr.def(get_N,                      "N");
    camfr.def(set_polarisation,           "set_polarisation");
    camfr.def(set_gain_material,          "set_gain_material");
    camfr.def(set_solver,                 "set_solver");
    camfr.def(set_stability,              "set_stability");
    camfr.def(set_precision,              "set_precision");
    camfr.def(set_precision_enhancement,  "set_precision_enhancement");
    camfr.def(set_dx_enhanced,            "set_dx_enhanced");
    camfr.def(set_precision_rad,          "set_precision_rad");
    camfr.def(set_C_upperright,           "set_C_upperright");
    camfr.def(set_sweep_from_previous,    "set_sweep_from_previous");
    camfr.def(set_sweep_steps,            "set_sweep_steps");
    camfr.def(set_eps_trace_coarse,       "set_eps_trace_coarse");
    camfr.def(set_chunk_tracing,          "set_chunk_tracing");
    camfr.def(set_unstable_exp_threshold, "set_unstable_exp_threshold");
    camfr.def(set_field_calc,             "set_field_calc");
    camfr.def(set_bloch_calc,             "set_bloch_calc");
    camfr.def(set_eigen_calc,             "set_eigen_calc");
    camfr.def(set_orthogonal,             "set_orthogonal");
    camfr.def(set_circ_order,             "set_circ_order");
    camfr.def(set_circ_fieldtype,         "set_circ_field_type");
    camfr.def(set_left_wall,              "set_left_wall");
    camfr.def(set_right_wall,             "set_right_wall");
    camfr.def(set_beta,                   "set_beta");
    camfr.def(free_tmps,                  "free_tmps");

    // Wrap Coord.

    boost::python::class_builder<Coord> Coord_(camfr, "Coord");

    Coord_.def(boost::python::constructor
               <const Complex&, const Complex&, const Complex&>());
    Coord_.def(boost::python::constructor
               <const Complex&, const Complex&, const Complex&,
                      Limit,          Limit,          Limit>());
    Coord_.def(&Coord::repr, "__repr__");

    // Wrap Field.

    boost::python::class_builder<Field> Field_(camfr, "Field");

    Field_.def(field_E1,      "E1");
    Field_.def(field_E2,      "E2");
    Field_.def(field_Ez,      "Ez");
    Field_.def(field_H1,      "H1");
    Field_.def(field_H2,      "H2");
    Field_.def(field_Hz,      "Hz");
    Field_.def(&Field::abs_E, "abs_E");
    Field_.def(&Field::abs_H, "abs_H");
    Field_.def(&Field::repr,  "__repr__");

    // Wrap FieldExpansion.

    boost::python::class_builder<FieldExpansion> 
      FieldExpansion_(camfr, "FieldExpansion");
    
    FieldExpansion_.def(&FieldExpansion::field, "field");
    FieldExpansion_.def(&FieldExpansion::repr,  "__repr__");
    
    // Wrap Material_length.

    boost::python::class_builder<Material_length>
      Material_length_(camfr, "Material_length");
    
    // Wrap BaseMaterial.

    boost::python::class_builder<BaseMaterial> 
      BaseMaterial_(camfr, "BaseMaterial");

    BaseMaterial_.def(&material_to_term, "__call__");
    
    // Wrap Material.

    boost::python::class_builder<Material> Material_(camfr, "Material");
    Material_.declare_base(BaseMaterial_);
    Material_.def(&material_to_term, "__call__");
    
    Material_.def(boost::python::constructor<const Complex&>());
    Material_.def(boost::python::constructor<const Complex&,const Complex&>());
    Material_.def(&Material::n,    "n");
    Material_.def(&Material::epsr, "epsr");
    Material_.def(&Material::mur,  "mur");
    Material_.def(&Material::eps,  "eps");
    Material_.def(&Material::mu,   "mu");
    Material_.def(&Material::gain, "gain");
    Material_.def(&Material::repr, "__repr__");

    // Wrap Waveguide_length.

    boost::python::class_builder<Waveguide_length>
      Waveguide_length_(camfr, "Waveguide_length");
    
    // Wrap Waveguide.

    boost::python::class_builder<Waveguide> Waveguide_(camfr, "Waveguide");

    Waveguide_.def(&Waveguide::get_core,    "core");
    Waveguide_.def(&Waveguide::eps_at,      "eps");
    Waveguide_.def(&Waveguide::mu_at,       "mu");
    Waveguide_.def(&Waveguide::n_at,        "n");   
    Waveguide_.def(&Waveguide::N,           "N");
    Waveguide_.def(waveguide_get_mode,      "mode");
    Waveguide_.def(waveguide_get_fw_mode,   "fw_mode");
    Waveguide_.def(waveguide_get_bw_mode,   "bw_mode");
    Waveguide_.def(&Waveguide::find_modes,  "calc");
    Waveguide_.def(&Waveguide::repr,        "__repr__");    
    Waveguide_.def(&waveguide_to_term,      "__call__");

    // Wrap MultiWaveguide.

    boost::python::class_builder<MultiWaveguide> 
      MultiWaveguide_(camfr, "MultiWaveguide");

    MultiWaveguide_.declare_base(Waveguide_);

    MultiWaveguide_.def(&MultiWaveguide::field_from_source,
                                                       "field_from_source");
    MultiWaveguide_.def(&waveguide_to_term,            "__call__"); // tmp
    
    // Wrap MonoWaveguide.

    boost::python::class_builder<MonoWaveguide> 
      MonoWaveguide_(camfr, "MonoWaveguide");

    MonoWaveguide_.declare_base(Waveguide_);
    MonoWaveguide_.def(&waveguide_to_term, "__call__"); // tmp  

    // Wrap Mode.

    boost::python::class_builder<Mode> Mode_(camfr, "Mode");

    Mode_.def(&Mode::field,    "field");
    Mode_.def(&Mode::n_eff,    "n_eff");
    Mode_.def(mode_kz,         "kz");
    Mode_.def(mode_pol,        "pol");
    Mode_.def(&Mode::repr,     "__repr__");

    // Wrap Scatterer.

    boost::python::class_builder<Scatterer> Scatterer_(camfr, "Scatterer");

    Scatterer_.def(&Scatterer::calcRT,  "calc");
    Scatterer_.def(&Scatterer::freeRT,  "free");
    Scatterer_.def(&Scatterer::get_inc, "inc");
    Scatterer_.def(&Scatterer::get_ext, "ext");
    Scatterer_.def(TERM_PLUS_TERM_F(Scatterer, const Term),       "__add__");
    Scatterer_.def(TERM_PLUS_TERM_F(Scatterer, Stack),            "__add__");
    Scatterer_.def(TERM_PLUS_TERM_F(Scatterer, const Expression), "__add__");
    Scatterer_.def(TERM_PLUS_TERM_F(Scatterer, Scatterer),        "__add__");
    
    // Wrap MultiScatterer.

    boost::python::class_builder<MultiScatterer> 
      MultiScatterer_(camfr, "MultiScatterer");
    MultiScatterer_.declare_base(Scatterer_);

    // Wrap DenseScatterer.

    boost::python::class_builder<DenseScatterer> 
      DenseScatterer_(camfr, "DenseScatterer");
    DenseScatterer_.declare_base(MultiScatterer_);

    // Wrap DiagScatterer.

    boost::python::class_builder<DiagScatterer> 
      DiagScatterer_(camfr, "DiagScatterer");
    DiagScatterer_.declare_base(MultiScatterer_);

    // Wrap MonoScatterer.

    boost::python::class_builder<MonoScatterer> 
      MonoScatterer_(camfr, "MonoScatterer");
    MonoScatterer_.declare_base(Scatterer_);

    // Wrap FlippedScatterer.

    boost::python::class_builder<FlippedScatterer>
      FlippedScatterer_(camfr, "FlippedScatterer");
    FlippedScatterer_.declare_base(MultiScatterer_);

    FlippedScatterer_.def(boost::python::constructor<MultiScatterer&>());

    // Wrap E_Wall.

    boost::python::class_builder<E_Wall> E_Wall_(camfr, "E_Wall");
    E_Wall_.declare_base(DiagScatterer_);   

    E_Wall_.def(boost::python::constructor<Waveguide&>());
      
    // Wrap H_Wall.

    boost::python::class_builder<H_Wall> H_Wall_(camfr, "H_Wall");
    H_Wall_.declare_base(DiagScatterer_);
    
    H_Wall_.def(boost::python::constructor<Waveguide&>());

    // Wrap Expression.

    boost::python::class_builder<Expression> Expression_(camfr, "Expression");

    Expression_.def(boost::python::constructor<>());
    Expression_.def(boost::python::constructor<const Term&>());
    Expression_.def(boost::python::constructor<const Expression&>());
    Expression_.def(&Expression::flatten,             "flatten");
    Expression_.def(&Expression::repr,                "__repr__");
    Expression_.def(&Expression::operator+=,          "add"); // todo: iadd
    Expression_.def(EX_PLUS_TERM_F(const Term),       "__add__");
    Expression_.def(EX_PLUS_TERM_F(Stack),            "__add__");
    Expression_.def(EX_PLUS_TERM_F(const Expression), "__add__");
    Expression_.def(EX_PLUS_TERM_F(Scatterer),        "__add__");
    Expression_.def(boost::python::operators<boost::python::op_mul>(), 
                    boost::python::right_operand<unsigned int>());
    Expression_.def(boost::python::operators<boost::python::op_mul>(),  
                    boost::python::left_operand<unsigned int>());
    
    // Wrap Term.

    boost::python::class_builder<Term> Term_(camfr, "Term");
    
    Term_.def(boost::python::constructor<Scatterer&>());
    Term_.def(boost::python::constructor<Stack&>());
    Term_.def(boost::python::constructor<const Expression&>());
    Term_.def(&Term::get_inc, "inc");
    Term_.def(&Term::get_ext, "ext");
    Term_.def(&Term::repr,    "__repr__");
    Term_.def(TERM_PLUS_TERM_F(const Term, const Term),       "__add__");
    Term_.def(TERM_PLUS_TERM_F(const Term, Stack),            "__add__");
    Term_.def(TERM_PLUS_TERM_F(const Term, const Expression), "__add__");
    Term_.def(TERM_PLUS_TERM_F(const Term, Scatterer),        "__add__");
    Term_.def(boost::python::operators<boost::python::op_mul>(), 
              boost::python::right_operand<unsigned int>());
    Term_.def(boost::python::operators<boost::python::op_mul>(),  
              boost::python::left_operand<unsigned int>());    
    
    // Wrap Stack.
    
    boost::python::class_builder<Stack> Stack_(camfr, "Stack");
    
    Stack_.def(boost::python::constructor<const Expression&>());
    Stack_.def(&Stack::calcRT,              "calc");
    Stack_.def(&Stack::freeRT,              "free");
    Stack_.def(&Stack::get_inc,             "inc");
    Stack_.def(&Stack::get_ext,             "ext");
    Stack_.def(&Stack::get_total_thickness, "length");
    Stack_.def(stack_set_inc_field,         "set_inc_field");
    Stack_.def(stack_set_inc_field_2,       "set_inc_field");
    Stack_.def(&Stack::get_inc_field,       "inc_field");
    Stack_.def(&Stack::get_refl_field,      "refl_field");
    Stack_.def(&Stack::get_trans_field,     "trans_field");
    Stack_.def(stack_inc_S_flux,            "inc_S_flux");
    Stack_.def(stack_ext_S_flux,            "ext_S_flux");    
    Stack_.def(&Stack::field,               "field");
    Stack_.def(stack_lateral_S_flux,        "lateral_S_flux");
    Stack_.def(stack_lateral_S_flux_2,      "lateral_S_flux");
    Stack_.def(&Stack::eps_at,              "eps");
    Stack_.def(&Stack::mu_at,               "mu");
    Stack_.def(&Stack::n_at,                "n");
    Stack_.def(stack_get_R12,               "R12");
    Stack_.def(stack_get_R21,               "R21");
    Stack_.def(stack_get_T12,               "T12");
    Stack_.def(stack_get_T21,               "T21");
    Stack_.def(stack_R12,                   "R12");
    Stack_.def(stack_R21,                   "R21");
    Stack_.def(stack_T12,                   "T12");
    Stack_.def(stack_T21,                   "T21");
    
    Stack_.def(TERM_PLUS_TERM_F(Stack, const Term),       "__add__");
    Stack_.def(TERM_PLUS_TERM_F(Stack, Stack),            "__add__");
    Stack_.def(TERM_PLUS_TERM_F(Stack, const Expression), "__add__");
    Stack_.def(TERM_PLUS_TERM_F(Stack, Scatterer),        "__add__");
    
    // Wrap Cavity.

    boost::python::class_builder<Cavity> Cavity_(camfr, "Cavity");

    Cavity_.def(boost::python::constructor<Stack&, Stack&>());
    Cavity_.def(cavity_find_modes_in_region_3, "find_modes_in_region");
    Cavity_.def(cavity_find_modes_in_region_7, "find_modes_in_region");
    Cavity_.def(cavity_find_mode_2,            "find_mode");
    Cavity_.def(cavity_find_mode_5,            "find_mode");
    Cavity_.def(cavity_calc_sigma,             "sigma");
    Cavity_.def(cavity_set_source,             "set_source");
    Cavity_.def(&Cavity::field,                "field");

    // Wrap BlochStack.

    boost::python::class_builder<BlochStack> BlochStack_(camfr, "BlochStack");
    BlochStack_.declare_base(MultiWaveguide_);

    BlochStack_.def(boost::python::constructor<const Expression&>());
    BlochStack_.def(blochstack_get_mode,              "mode");
    BlochStack_.def(&BlochStack::get_total_thickness, "length");
    BlochStack_.def(&BlochStack::get_beta_vector,     "beta_vector");
    BlochStack_.def(&BlochStack::repr,                "__repr__");
    
    // Wrap BlochMode.

    boost::python::class_builder<BlochMode> BlochMode_(camfr, "BlochMode");
    BlochMode_.declare_base(Mode_);

    BlochMode_.def(&BlochMode::fw_field, "fw_field");
    BlochMode_.def(&BlochMode::bw_field, "bw_field");
    BlochMode_.def(&BlochMode::S_flux,   "S_flux");

    // Wrap InfStack.

    boost::python::class_builder<InfStack> InfStack_(camfr, "InfStack");
    InfStack_.declare_base(DenseScatterer_);

    InfStack_.def(boost::python::constructor<const Expression&>());
    InfStack_.def(boost::python::constructor
                  <const Expression&,const Complex&>()); // tmp
    InfStack_.def(&InfStack::get_R12,  "R12");

    // Wrap RealFunction.

    boost::python::class_builder<RealFunction> 
      RealFunction_(camfr, "RealFunction");

    RealFunction_.def(&RealFunction::times_called, "times_called");
    RealFunction_.def(&RealFunction::operator(),   "__call__");

    // Wrap ComplexFunction.

    boost::python::class_builder<ComplexFunction> 
      ComplexFunction_(camfr, "ComplexFunction");

    ComplexFunction_.def(&ComplexFunction::times_called, "times_called");
    ComplexFunction_.def(&ComplexFunction::operator(),   "__call__");

    // Wrap Planar.

    boost::python::class_builder<Planar> Planar_(camfr, "Planar");
    Planar_.declare_base(MonoWaveguide_);
    
    Planar_.def(boost::python::constructor<Material&>());
    Planar_.def(&Planar::set_theta, "set_theta");
    Planar_.def(&Planar::repr,     "__repr__"); // tmp
    Planar_.def(waveguide_to_term, "__call__"); // tmp 
    
    // Wrap Circ.

    boost::python::class_builder<Circ> Circ_(camfr, "Circ");
    Circ_.declare_base(MultiWaveguide_);
    
    Circ_.def(constructor<const Term&>());
    Circ_.def(constructor<const Expression&>());    
    Circ_.def(&Circ::repr,       "__repr__"); // tmp
    Circ_.def(waveguide_to_term, "__call__"); // tmp

    // Wrap SlabWall.

    boost::python::class_builder<SlabWall> SlabWall_(camfr, "SlabWall");

    SlabWall_.def(&SlabWall::get_R12, "R");

    // Wrap SlabWallMixed.

    boost::python::class_builder<SlabWallMixed> 
      SlabWallMixed_(camfr, "SlabWallMixed");
    SlabWallMixed_.declare_base(SlabWall_);

    SlabWallMixed_.def(constructor<const Complex&, const Complex&>());

    // Wrap SlabWall_TBC.

    boost::python::class_builder<SlabWall_TBC> 
      SlabWall_TBC_(camfr, "SlabWall_TBC");
    SlabWall_TBC_.declare_base(SlabWall_);

    SlabWall_TBC_.def(boost::python::constructor
                      <const Complex&, const Material&>());

    camfr.add(make_ref(slab_E_wall),    "slab_E_wall");
    camfr.add(make_ref(slab_H_wall),    "slab_H_wall");
    camfr.add(make_ref(slab_open_wall), "slab_open_wall");

    // Wrap SlabWall_PC.

    boost::python::class_builder<SlabWall_PC> 
      SlabWall_PC_(camfr, "SlabWall_PC");
    SlabWall_PC_.declare_base(SlabWall_);

    SlabWall_PC_.def(constructor<const Expression&>());

    // Wrap SlabDisp.

    boost::python::class_builder<SlabDisp> SlabDisp_(camfr, "SlabDisp");
    SlabDisp_.declare_base(ComplexFunction_);

    SlabDisp_.def(boost::python::constructor<Expression&, Real>());
    SlabDisp_.def(boost::python::constructor
                  <Expression&, Real, SlabWall*, SlabWall*>());
    SlabDisp_.def(&SlabDisp::operator(), "__call__"); // tmp   

    // Wrap Slab.

    boost::python::class_builder<Slab> Slab_(camfr, "Slab");
    Slab_.declare_base(MultiWaveguide_);

    Slab_.def(boost::python::constructor<const Term&>());
    Slab_.def(boost::python::constructor<const Expression&>());

    Slab_.def(&Slab::set_left_wall,  "set_left_wall");
    Slab_.def(&Slab::set_right_wall, "set_right_wall");
    Slab_.def(&Slab::get_width,      "width");
    Slab_.def(&Slab::repr,           "__repr__"); // tmp
    Slab_.def(waveguide_to_term,     "__call__"); // tmp

    // Wrap SectionDisp.

    boost::python::class_builder<SectionDisp> 
      SectionDisp_(camfr, "SectionDisp");
    SectionDisp_.declare_base(ComplexFunction_);

    SectionDisp_.def(boost::python::constructor<Stack&, Stack&, Real, int>());
    SectionDisp_.def(&SectionDisp::operator(), "__call__"); // tmp

    // Wrap Section.

    boost::python::class_builder<Section> Section_(camfr, "Section");
    Section_.declare_base(MultiWaveguide_);

    Section_.def(boost::python::constructor<const Expression&>());
    Section_.def(boost::python::constructor<const Expression&,int>());
    Section_.def(boost::python::constructor
                 <const Expression&,const Expression&>());
    Section_.def(boost::python::constructor
                 <const Expression&,const Expression&,int>());

    Section_.def(&Section::get_width, "width");
    Section_.def(&Section::repr,      "__repr__"); // tmp
  }
  catch(...)
  {
    //handle_exception();
  }
}
