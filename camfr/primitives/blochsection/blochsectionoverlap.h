 
/////////////////////////////////////////////////////////////////////////////
//
// File:     blochsectionoverlap.h
// Author:   Peter.Bienstman@UGent.be
// Date:     20050518
// Version:  1.0
//
// Copyright (C) 2005 Peter Bienstman - Ghent University
//
/////////////////////////////////////////////////////////////////////////////

#ifndef BLOCHSECTIONOVERLAP_H
#define BLOCHSECTIONOVERLAP_H

#include "../../field.h"

/////////////////////////////////////////////////////////////////////////////
//
// Calculates overlapintegral Int(E_I x H_II) between mode_I and mode_II,
//
/////////////////////////////////////////////////////////////////////////////

class BlochSectionMode;     // forward declaration - see blochsectionmode.h
class BlochSection2D_Mode;  // forward declaration - see blochsectionmode.h

Complex overlap(const BlochSection2D_Mode* sec_I_mode,
                const BlochSection2D_Mode* sec_II_mode);

#endif



