#ifndef BR_H
#define BR_H

#include "cfg_decoding.h"

extern int br_poly_init();
extern int br_poly_exit();
extern int br_poly_clear();
extern int br_g_poly_gen();
extern int br_lag_poly_construct();
extern int br_k_poly_construct(long long tv_idx);
extern int br_k_poly_ded(long long base_tv_idx, long long ded_tv_idx);
extern int br_m_poly_construct(long long tv_idx);
extern int br_v_matric_gen(long long tv_idx);
extern int is_popov_form(long long tv_idx);
extern int ms_reduction(long long tv_idx);
extern int br_g_poly_ret(long long tv_idx);
extern int br_q_poly_ret(long long tv_idx);
extern int br_y_high_deg_trans(unsigned char *tras_poly);
extern int br_v_matric_gen(long long tv_idx);
extern int br_test();

#endif
