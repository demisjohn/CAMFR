 
/////////////////////////////////////////////////////////////////////////////
//
// File:     sectionoverlap.h
// Author:   Peter.Bienstman@rug.ac.be
// Date:     20020225
// Version:  1.0
//
// Copyright (C) 2002 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SECTIONOVERLAP_H
#define SECTIONOVERLAP_H

#include "../../field.h"

/////////////////////////////////////////////////////////////////////////////
//
// Calculates overlapintegral Int(E_I x H_II) between mode_I and mode_II,
// for the slice which ranges from z_start to z_stop.
//
// The field expansion of these modes can be passed in advance in field_I
// and field_II, and the overlapintegrals can be passed in a cache, indexed
// by I_index and II_index.
//
/////////////////////////////////////////////////////////////////////////////

class SectionMode;     // forward declaration - see sectionmode.h
class Section2D_Mode;  // forward declaration - see sectionmode.h
class OverlapMatrices; // forward declaration - see slabmatrixcache.h

Complex overlap_slice(const SectionMode* mode_I,
                      const SectionMode* mode_II,
                      const Complex& z_start,
                      const Complex& z_stop,
                      FieldExpansion* field_I=NULL,
                      FieldExpansion* field_II=NULL,
                      OverlapMatrices* m=NULL,
                      int I_index=1, int II_index=1);

Complex overlap_numeric(const SectionMode* mode_I,
                        const SectionMode* mode_II);

Complex overlap_pw(const Section2D_Mode* sec_I_mode, 
                   const Section2D_Mode* sec_II_mode);

Complex overlap(const Section2D_Mode* sec_I_mode,
                const Section2D_Mode* sec_II_mode);

#endif



