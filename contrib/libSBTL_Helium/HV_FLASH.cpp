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
// HV_FLASH
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_HE.h"
#include "SBTL_call_conv.h"
//
#define ITMAX 10

//initial guess from auxiliary splines
extern "C" double __stdcall U_VH_HE_INI_T(double vt, double h);
//
//forward functions with derivatives
extern "C" void __stdcall DIFF_P_VU_HE_T(double vt, double v, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_HE_TT(double vt, double u, double& p, double& dpdv, double& dpdu, double& dudv);
//
SBTLAPI int __stdcall FLASH_VH_HE(double v, double h, double& u) throw()
{
    double vt;
    static const double df_h=1.e-8;  //-8   //abs. deviation in h

    double hx,px;
    double dhdu_v;
    double dpdv_u, dpdu_v, dudv_p;

    vt=log(v);

    //calculate initial guess
    u=U_VH_HE_INI_T(vt,h);

    //newtons method
    double f_h=-1.;
    int icount=0;
    while(fabs(f_h)>df_h) {
        DIFF_P_VU_HE_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        hx=u+px*v*1.e3;
        dhdu_v=1.+dpdu_v*v*1.e3;
        f_h=hx-h;
        u=u-f_h/dhdu_v;
        if(icount++>ITMAX) {
            u=ERR_VAL;
            return I_ERR;
        }
    }
    return I_OK;
}
//
SBTLAPI int __stdcall FLASH_VH_HE_T(double v, double vt, double h, double& u) throw()
{

    static const double df_h=1.e-8;  //-8   //abs. deviation in h

    double hx,px;
    double dhdu_v;
    double dpdv_u, dpdu_v, dudv_p;

    //calculate initial guess
    u=U_VH_HE_INI_T(vt,h);

    //newtons method
    double f_h=-1.;
    int icount=0;
    while(fabs(f_h)>df_h) {
        DIFF_P_VU_HE_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        hx=u+px*v*1.e3;
        dhdu_v=1.+dpdu_v*v*1.e3;
        f_h=hx-h;
        u=u-f_h/dhdu_v;
        if(icount++>ITMAX) {
            u=ERR_VAL;
            return I_ERR;
        }
    }
    return I_OK;
}
//
SBTLAPI void __stdcall VH_FLASH_DERIV_HE(double v, double vt, double u, double& dudv_h, double& dudh_v, double& dvdh_u) throw()
{
    double dpdv_u, dpdu_v, dudv_p;
    double dhdv_u, dhdu_v;
    double p_;

    //derivatives
    DIFF_P_VU_HE_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
    dhdv_u=(dpdv_u*v+p_)*1.e3;
    dhdu_v=1.+dpdu_v*v*1.e3;
    dudv_h=-dhdv_u/dhdu_v;
    //
    dudh_v=1./dhdu_v;
    dvdh_u=1./dhdv_u;
}
//
