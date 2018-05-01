///////////////////////////////////////////////////////////////////////////
// LibSBTL_vu_HE - SBTL library for gaseous helium based on:
//
// Ortiz-Vega, D.O., Hall, K.R., Holste, J.C., Arp, V.D., Harvey, A.H., and Lemmon, E.W.:
//
//   final equation of state, to be published in J. Phys. Chem. Ref. Data, 2018.
//
// Arp, V.D., McCarty, R.D., and Friend, D.G.:
//
//  "Thermophysical Properties of Helium-4 from 0.8 to 1500 K with
//   Pressures to 2000 MPa,"
//   NIST Technical Note 1334 (revised), 1998.
//
// Arp, V.D., McCarty, R.D., and Friend, D.G.:
//
//  "Thermophysical Properties of Helium-4 from 0.8 to 1500 K with
//   Pressures to 2000 MPa,"
//   NIST Technical Note 1334 (revised), 1998.
//
// Hands, B.A. and Arp, V.D.:
//
//  "A Correlation of Thermal Conductivity Data for Helium,"
//   Cryogenics, 21(12):697-703, 1981."
//
// Copyright (C) Idaho National Laboratory.
// All rights reserved.
//
// Disclaimer:
// The Idaho National Laboratory (INL) uses its best efforts to deliver a high-quality software and to verify that the computed information is correct.
// However, INL makes no warranties to that effect, and INL shall not be liable for any damage that may result from errors or omissions in the software.
//
// Version: 0.9.0
//
// SBTL_def.h
//
#pragma once

#ifdef __SSE2__
//
#include "emmintrin.h"
#define IROUND(d) (_mm_cvttsd_si32(_mm_load_sd(&d)))
//
#else
//
#pragma message ("IROUND: casting operation. Enable SSE2, if possible.")
#define IROUND(d) ((unsigned int)(d))
//
#endif
