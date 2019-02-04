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
// The Idaho National Laboratory (INL) uses its best efforts to deliver a high-quality software and
// to verify that the computed information is correct. However, INL makes no warranties to that
// effect, and INL shall not be liable for any damage that may result from errors or omissions in
// the software.
//
// Version: 0.9.0
//
// PT_FLASH
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_HE.h"
#include "SBTL_call_conv.h"
//
#define ITMAX 10
//
// initial guess from auxiliary splines
extern "C" void __stdcall VU_TP_HE_INI(double t, double p, double & vt, double & u);
//
// forward functions with derivatives
extern "C" void __stdcall DIFF_P_VU_HE_T(
    double vt, double v, double u, double & p, double & dpdv, double & dpdu, double & dudv);
extern "C" void __stdcall DIFF_T_VU_HE_T(
    double vt, double v, double u, double & t, double & dtdv, double & dtdu, double & dudv);
extern "C" void __stdcall DIFF_P_VU_HE_TT(
    double vt, double u, double & p, double & dpdv, double & dpdu, double & dudv);
extern "C" void __stdcall DIFF_T_VU_HE_TT(
    double vt, double u, double & t, double & dtdv, double & dtdu, double & dudv);
//
SBTLAPI int __stdcall PT_FLASH_HE(double p, double t, double & v, double & vt, double & u) throw()
{
  static const double df_p = 1.e-10; // rel. deviation in p
  static const double df_t = 1.e-10; // abs. deviation in t

  double tx, px, den;
  double dtdv_u, dtdu_v, dudv_t;
  double dpdv_u, dpdu_v, dudv_p;

  // calculate initial guesses
  VU_TP_HE_INI(t, p, vt, u);

  // newtons method
  double f_p = -1., f_t = -1., p_inv = 1. / p;
  int icount = 0;
  while (fabs(f_p * p_inv) > df_p || fabs(f_t) > df_t)
  {
    DIFF_P_VU_HE_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p); // px, transformed derivatives
    DIFF_T_VU_HE_TT(vt, u, tx, dtdv_u, dtdu_v, dudv_t); // tx, transformed derivatives
    f_p = px - p;
    f_t = tx - t;
    den = dtdu_v * dpdv_u - dtdv_u * dpdu_v;
    vt = vt + (-dtdu_v * f_p + f_t * dpdu_v) / den;
    u = u + (-f_t * dpdv_u + dtdv_u * f_p) / den;
    if (icount++ > ITMAX)
    {
      return I_ERR;
    }
  }
  v = exp(vt);
  return I_OK;
}
//
SBTLAPI int __stdcall PT_FLASH_DERIV_HE(double p,
                                        double t,
                                        double & v,
                                        double & vt,
                                        double & dvdp_t,
                                        double & dvdt_p,
                                        double & dpdt_v,
                                        double & u,
                                        double & dudp_t,
                                        double & dudt_p,
                                        double & dpdt_u) throw()
{
  static const double df_p = 1.e-10; // rel. deviation in p
  static const double df_t = 1.e-10; // abs. deviation in t

  double tx, px, den;
  double dtdv_u, dtdu_v, dudv_t;
  double dpdv_u, dpdu_v, dudv_p;
  double p_, t_;

  // calculate initial guesses
  VU_TP_HE_INI(t, p, vt, u);

  // newtons method
  double f_p = -1., f_t = -1., p_inv = 1. / p;
  int icount = 0;
  while (fabs(f_p * p_inv) > df_p || fabs(f_t) > df_t)
  {
    DIFF_P_VU_HE_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p); // px, transformed derivatives
    DIFF_T_VU_HE_TT(vt, u, tx, dtdv_u, dtdu_v, dudv_t); // tx, transformed derivatives
    f_p = px - p;
    f_t = tx - t;
    den = dtdu_v * dpdv_u - dtdv_u * dpdu_v;
    vt = vt + (-dtdu_v * f_p + f_t * dpdu_v) / den;
    u = u + (-f_t * dpdv_u + dtdv_u * f_p) / den;
    if (icount++ > ITMAX)
    {
      return I_ERR;
    }
  }
  v = exp(vt);
  // derivatives
  DIFF_P_VU_HE_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
  DIFF_T_VU_HE_T(vt, v, u, t_, dtdv_u, dtdu_v, dudv_t);
  //
  dvdp_t = 1. / (dpdv_u + dpdu_v * dudv_t);
  dvdt_p = 1. / (dtdv_u + dtdu_v * dudv_p);
  dpdt_v = -dvdt_p / dvdp_t;
  //
  dudp_t = 1. / (dpdu_v + dpdv_u / dudv_t);
  dudt_p = 1. / (dtdu_v + dtdv_u / dudv_p);
  dpdt_u = -dudt_p / dudp_t;
  return I_OK;
}
//
SBTLAPI int
PT_FLASH_HE_T(double p, double t, double & vt, double & u) throw()
{
  static const double df_p = 1.e-10; //-8   //rel. deviation in p
  static const double df_t = 1.e-10; //-8   //abs. deviation in t

  double tx, px, den;
  double dtdv_u, dtdu_v, dudv_t;
  double dpdv_u, dpdu_v, dudv_p;

  // calculate initial guesses
  VU_TP_HE_INI(t, p, vt, u);

  // newtons method
  double f_p = -1., f_t = -1., p_inv = 1. / p;
  int icount = 0;
  while (fabs(f_p * p_inv) > df_p || fabs(f_t) > df_t)
  {
    DIFF_P_VU_HE_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p); // px, transformed derivatives
    DIFF_T_VU_HE_TT(vt, u, tx, dtdv_u, dtdu_v, dudv_t); // tx, transformed derivatives
    f_p = px - p;
    f_t = tx - t;
    den = dtdu_v * dpdv_u - dtdv_u * dpdu_v;
    vt = vt + (-dtdu_v * f_p + f_t * dpdu_v) / den;
    u = u + (-f_t * dpdv_u + dtdv_u * f_p) / den;
    if (icount++ > ITMAX)
    {
      return I_ERR;
    }
  }
  return I_OK;
}
