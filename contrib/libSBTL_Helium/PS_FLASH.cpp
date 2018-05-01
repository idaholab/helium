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
// PS_FLASH
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_HE.h"
#include "SBTL_call_conv.h"
//
#define ITMAX 10
//
//initial guess from auxiliary splines
extern "C" void  __stdcall VU_SP_HE_INI(double s, double p, double& vt, double& u);
//
//forward functions with derivatives
extern "C" void __stdcall DIFF_P_VU_HE_T(double vt, double v, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_HE_T(double vt, double v, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_HE_TT(double vt, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_HE_TT(double vt, double u, double& s, double& dsdv, double& dsdu, double& dudv);
//
SBTLAPI int __stdcall PS_FLASH_HE(double p, double s, double& v, double& vt, double& u) throw()
{
    static const double df_p=1.e-10; //-8   //rel. deviation in p
    static const double df_s=1.e-10; //-8   //abs. deviation in t

    double sx,px,den;
    double dsdv_u, dsdu_v, dudv_s;
    double dpdv_u, dpdu_v, dudv_p;

    //calculate initial guesses
    VU_SP_HE_INI(s, p, vt, u);

    //newtons method
    double f_p=-1.,f_s=-1.,p_inv=1./p;
    int icount=0;
    while(fabs(f_p*p_inv)>df_p || fabs(f_s)>df_s) {
        DIFF_P_VU_HE_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        DIFF_S_VU_HE_TT(vt, u, sx, dsdv_u, dsdu_v, dudv_s);      // sx, transformed derivatives
        f_p=px-p;
        f_s=sx-s;
        den=dsdu_v*dpdv_u-dsdv_u*dpdu_v;
        vt=vt+(-dsdu_v*f_p+f_s*dpdu_v)/den;
        u =u +(-f_s*dpdv_u+dsdv_u*f_p)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    v=exp(vt);
    return I_OK;
}
//
SBTLAPI void __stdcall PS_FLASH_DERIV_HE(double v, double vt, double u, double& dvdp_s, double& dvds_p, double& dpds_v, double& dudp_s, double& duds_p, double& dpds_u) throw()
{
    double dsdv_u, dsdu_v, dudv_s;
    double dpdv_u, dpdu_v, dudv_p;
    double p_,s_;

    //derivatives
    DIFF_P_VU_HE_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
    DIFF_S_VU_HE_T(vt, v, u, s_, dsdv_u, dsdu_v, dudv_s);
    //
    dvdp_s=1./(dpdv_u+dpdu_v*dudv_s);
    dvds_p=1./(dsdv_u+dsdu_v*dudv_p);
    dpds_v=-dvds_p/dvdp_s;
    //
    dudp_s=1./(dpdu_v+dpdv_u/dudv_s);
    duds_p=1./(dsdu_v+dsdv_u/dudv_p);
    dpds_u=-duds_p/dudp_s;
}
//
SBTLAPI int PS_FLASH_HE_T(double p, double s, double& vt, double& u) throw()
{
    static const double df_p=1.e-10; //-8   //rel. deviation in p
    static const double df_s=1.e-10; //-8   //abs. deviation in s

    double sx,px,den;
    double dsdv_u, dsdu_v, dudv_s;
    double dpdv_u, dpdu_v, dudv_p;

    //calculate initial guesses
    VU_SP_HE_INI(s, p, vt, u);

    //newtons method
    double f_p=-1.,f_s=-1.,p_inv=1./p;
    int icount=0;
    while(fabs(f_p*p_inv)>df_p || fabs(f_s)>df_s) {
        DIFF_P_VU_HE_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        DIFF_S_VU_HE_TT(vt, u, sx, dsdv_u, dsdu_v, dudv_s);      // sx, transformed derivatives
        f_p=px-p;
        f_s=sx-s;
        den=dsdu_v*dpdv_u-dsdv_u*dpdu_v;
        vt=vt+(-dsdu_v*f_p+f_s*dpdu_v)/den;
        u =u +(-f_s*dpdv_u+dsdv_u*f_p)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    return I_OK;
}
