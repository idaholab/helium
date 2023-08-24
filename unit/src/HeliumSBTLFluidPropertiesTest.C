//* This file is part of helium
//* https://github.com/idaholab/helium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/helium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeliumSBTLFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(HeliumSBTLFluidPropertiesTest, test)
{
  const Real T = 120.0 + 273.15;
  const Real p = 101325;

  const Real rho_from_p_T = _fp->rho_from_p_T(p, T);
  const Real rho = rho_from_p_T;

  const Real h_from_p_T = _fp->h_from_p_T(p, T);
  const Real h = h_from_p_T;

  const Real e_from_p_rho = _fp->e_from_p_rho(p, rho);
  const Real e = e_from_p_rho;

  const Real v = 1 / rho;

  const Real s_from_v_e = _fp->s_from_v_e(v, e);
  const Real s = s_from_v_e;

  // p
  REL_TEST(_fp->p_from_v_e(v, e), p, REL_TOL_CONSISTENCY);
  REL_TEST(_fp->p_from_h_s(h, s), p, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->p_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->p_from_h_s, h, s, REL_TOL_DERIVATIVE);

  // T
  REL_TEST(_fp->T_from_v_e(v, e), T, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->T_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // rho
  // TODO: REL_TEST(rho, rho_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(rho_from_p_T, 0.12402618805565359, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->rho_from_p_s(p, s), rho_from_p_T, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->rho_from_p_T, p, T, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->rho_from_p_s, p, s, REL_TOL_DERIVATIVE);

  // e
  // TODO: REL_TEST(e, e_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(e_from_p_rho, 1230100.5223275987, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->e_from_v_h(v, h), e, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->e_from_p_rho, p, rho, 5e-4);
  DERIV_TEST(_fp->e_from_v_h, v, h, 2e-6);

  // c
  const Real c = _fp->c_from_v_e(v, e);
  // TODO: REL_TEST(c, c_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(c, 1167.0498683516942, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->c_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // cp
  const Real cp = _fp->cp_from_v_e(v, e);
  // TODO: REL_TEST(cp, cp_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(cp, 5193.0901602885488, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->cp_from_p_T(p, T), 5193.0901602885488, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->cp_from_v_e, v, e, 0.001); // allow 0.1% here (numerical derivative)
  DERIV_TEST(_fp->cp_from_p_T, p, T, 0.001); // allow 0.1% here (numerical derivative)

  // cv
  const Real cv = _fp->cv_from_v_e(v, e);
  // TODO: REL_TEST(cv, cv_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(cv, 3116.0786782056684, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->cv_from_p_T(p, T), 3116.0786782056684, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->cv_from_v_e, v, e, 1e-4); // allow 0.01% here (numerical derivative)

  // mu
  Real mu = _fp->mu_from_v_e(v, e);
  // TODO: REL_TEST(mu, mu_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(mu, 0.000024003810057301473, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->mu_from_p_T(p, T), 0.000024003810057301473, REL_TOL_SAVED_VALUE);
  Real dmu_dv, dmu_de;
  _fp->mu_from_v_e(v, e, mu, dmu_de, dmu_dv);
  REL_TEST(mu, 0.000024003810057301473, REL_TOL_SAVED_VALUE);
  REL_TEST(dmu_dv, 0., REL_TOL_SAVED_VALUE);
  REL_TEST(dmu_de, 0., REL_TOL_SAVED_VALUE);

  // k
  const Real k = _fp->k_from_v_e(v, e);
  // TODO: REL_TEST(k, k_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(k, 0.18809848357335232, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->k_from_p_T(p, T), 0.18809848357335232, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->k_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->k_from_p_T, p, T, REL_TOL_DERIVATIVE);

  // s
  REL_TEST(s, 29384.437673339806, REL_TOL_EXTERNAL_VALUE);
  // TODO: REL_TEST(s, s_saved, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->s_from_h_p(h, p), s, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->s_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->s_from_h_p, h, p, REL_TOL_DERIVATIVE);

  // g
  const Real g = _fp->g_from_v_e(v, e);
  REL_TEST(g, -9505426.5901812277, REL_TOL_EXTERNAL_VALUE);
  // TODO: REL_TEST(g, g_saved, REL_TOL_SAVED_VALUE);

  // h
  // TODO: REL_TEST(h, h_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(h_from_p_T, 2047065.0810911956, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->h_from_p_T, p, T, REL_TOL_DERIVATIVE);

  // beta
  const Real beta = _fp->beta_from_p_T(p, T);
  REL_TEST(beta, 0.00254252762578, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(beta, 0.0025424810406517082, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->beta_from_p_T, p, T, REL_TOL_DERIVATIVE);

  // molar mass
  REL_TEST(_fp->molarMass(), 4.002602e-3, REL_TOL_SAVED_VALUE);
}
