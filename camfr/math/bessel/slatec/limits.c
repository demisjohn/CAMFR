/**************************************************************************/
/*                                                                        */
/* File:     limits.c                                                     */
/* Author:   Peter.Bienstman@rug.ac.be                                    */
/* Date:     20000103                                                     */
/* Version:  1.0                                                          */
/*                                                                        */
/* Copyright (C) 2000 Peter Bienstman <Peter.Bienstman@rug.ac.be>         */
/*                                                                        */
/**************************************************************************/

/**************************************************************************/
/*                                                                        */
/* Wrapper to generate both single and double precision                   */
/* versions of machar.                                                    */
/*                                                                        */
/**************************************************************************/

#define SP
#include "machar.c"
#undef SP

#define DP
#include "machar.c"
