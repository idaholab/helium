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
// HS_FLASH
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_HE.h"
#include "SBTL_call_conv.h"
//
#define ITMAX 10

// initial guess from auxiliary splines
extern "C" void __stdcall VU_SH_HE_INI(double s, double h, double & vt, double & u);
//
// forward functions with derivatives
extern "C" void __stdcall DIFF_P_VU_HE_T(
    double vt, double v, double u, double & p, double & dpdv, double & dpdu, double & dudv);
extern "C" void __stdcall DIFF_T_VU_HE_T(
    double vt, double v, double u, double & t, double & dtdv, double & dtdu, double & dudv);
extern "C" void __stdcall DIFF_S_VU_HE_T(
    double vt, double v, double u, double & s, double & dsdv, double & dsdu, double & dudv);
extern "C" void __stdcall DIFF_S_VU_HE_TT(
    double vt, double u, double & s, double & dsdv, double & dsdu, double & dudv);
extern "C" void __stdcall DIFF_P_VU_HE_TT(
    double vt, double u, double & p, double & dpdv, double & dpdu, double & dudv);
//
SBTLAPI int __stdcall HS_FLASH_HE(double h, double s, double & v, double & vt, double & u) throw()
{
  static const double df_h = 1.e-8;  //-8   //abs. deviation in h
  static const double df_s = 1.e-10; //-8   //abs. deviation in s

  double hx, sx, px, den;
  double dhdv_u, dhdu_v;
  double dpdv_u, dpdu_v, dudv_p;
  double dsdv_u, dsdu_v, dudv_s;

  // calculate initial guesses
  VU_SH_HE_INI(s, h, vt, u);
  v = exp(vt);

  // newtons method
  double f_h = -1., f_s = -1.;
  int icount = 0;
  while (fabs(f_h) > df_h || fabs(f_s) > df_s)
  {
    DIFF_P_VU_HE_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p); // px, transformed derivatives
    DIFF_S_VU_HE_TT(vt, u, sx, dsdv_u, dsdu_v, dudv_s); // sx, transformed derivatives
    hx = u + px * v * 1.e3;
    dhdv_u = (dpdv_u * v + px * v) * 1.e3;
    dhdu_v = 1. + dpdu_v * v * 1.e3;
    f_h = hx - h;
    f_s = sx - s;
    den = dhdv_u * dsdu_v - dhdu_v * dsdv_u;
    vt = vt + (-dsdu_v * f_h + f_s * dhdu_v) / den;
    u = u + (-f_s * dhdv_u + dsdv_u * f_h) / den;
    v = exp(vt);
    if (icount++ > ITMAX)
    {
      return I_ERR;
    }
  }
  return I_OK;
}
//
SBTLAPI void __stdcall HS_FLASH_DERIV_HE(double v,
                                         double vt,
                                         double u,
                                         double & dvdh_s,
                                         double & dvds_h,
                                         double & dhds_v,
                                         double & dudh_s,
                                         double & duds_h,
                                         double & dhds_u) throw()
{
  double dhdv_u, dhdu_v, dudv_h;
  double dpdv_u, dpdu_v, dudv_p;
  double dsdv_u, dsdu_v, dudv_s;
  double p_, s_;

  // derivatives
  DIFF_P_VU_HE_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
  DIFF_S_VU_HE_T(vt, v, u, s_, dsdv_u, dsdu_v, dudv_s);
  dhdv_u = (dpdv_u * v + p_) * 1.e3;
  dhdu_v = 1. + dpdu_v * v * 1.e3;
  dudv_h = -dhdv_u / dhdu_v;
  //
  dvdh_s = 1. / (dhdv_u + dhdu_v * dudv_s);
  dvds_h = 1. / (dsdv_u + dsdu_v * dudv_h);
  dhds_v = -dvds_h / dvdh_s;
  //
  dudh_s = 1. / (dhdu_v + dhdv_u / dudv_s);
  duds_h = 1. / (dsdu_v + dsdv_u / dudv_h);
  dhds_u = -duds_h / dudh_s;
}
//
SBTLAPI void __stdcall HS_PT_FLASH_DERIV_HE(double v,
                                            double vt,
                                            double u,
                                            double & p,
                                            double & dpdh_s,
                                            double & dpds_h,
                                            double & dhds_p,
                                            double & t,
                                            double & dtdh_s,
                                            double & dtds_h,
                                            double & dhds_t) throw()
{
  double dhdv_u, dhdu_v;
  double dpdv_u, dpdu_v, dudv_p;
  double dtdv_u, dtdu_v, dudv_t;
  double dsdv_u, dsdu_v, dudv_s;
  double s_;

  // derivatives
  DIFF_P_VU_HE_T(vt, v, u, p, dpdv_u, dpdu_v, dudv_p);
  DIFF_T_VU_HE_T(vt, v, u, t, dtdv_u, dtdu_v, dudv_t);
  DIFF_S_VU_HE_T(vt, v, u, s_, dsdv_u, dsdu_v, dudv_s);
  dhdv_u = (dpdv_u * v + p) * 1.e3;
  dhdu_v = 1. + dpdu_v * v * 1.e3;
  //
  dpdh_s = (dpdv_u * dsdu_v - dpdu_v * dsdv_u) / (dhdv_u * dsdu_v - dhdu_v * dsdv_u);
  dpds_h = (dpdv_u * dhdu_v - dpdu_v * dhdv_u) / (dsdv_u * dhdu_v - dsdu_v * dhdv_u);
  dhds_p = -dpds_h / dpdh_s;
  //
  dtdh_s = (dtdv_u * dsdu_v - dtdu_v * dsdv_u) / (dhdv_u * dsdu_v - dhdu_v * dsdv_u);
  dtds_h = (dtdv_u * dhdu_v - dtdu_v * dhdv_u) / (dsdv_u * dhdu_v - dsdu_v * dhdv_u);
  dhds_t = -dtds_h / dtdh_s;
}
