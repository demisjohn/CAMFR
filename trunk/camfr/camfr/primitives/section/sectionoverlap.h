
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
// Calculates overlapintegral Int(E_I x H_II) between mode_I and mode_II.
//
/////////////////////////////////////////////////////////////////////////////

class SectionMode; // forward declaration - see sectionmode.h

Complex overlap_slice(SectionMode* mode_I,
                      SectionMode* mode_II,
                      const Complex& z_start,
                      const Complex& z_stop,
                      FieldExpansion* field_I=NULL,
                      FieldExpansion* field_II=NULL);

Complex overlap_numeric(const SectionMode* mode_I,
                        const SectionMode* mode_II);



#endif



