#ifndef FACTORIZTION_H
#define FACTORIZTION_H

#include "cfg_decoding.h"

extern long long dev_to_cnt;
extern long long sing_era_cnt;
extern unsigned char erasure_flag[CODEWORD_LEN];
extern unsigned char dev_to_flag[CODEWORD_LEN];
extern unsigned char store_q0_dev[MAX_POLY_TERM_SIZE];
extern unsigned char store_q1_dev[MAX_POLY_TERM_SIZE];
extern unsigned char store_v_dev[MAX_POLY_TERM_SIZE];

extern int g_poly_trans(unsigned char *poly, unsigned char *g_poly, long long pole_idx, long long z_degree);
extern long long lt_get(unsigned char *poly);
extern long long chien_search_x_y(unsigned char *poly);
extern unsigned char recur_root(unsigned char *poly, unsigned char *g_poly);
extern int recur_poly_update(unsigned char *poly, long long pole_idx, unsigned char root_prev_coef);
extern int factorization_recur(unsigned char *poly, unsigned char *est_msg);
extern int poly_dev_cal(unsigned char *poly, unsigned char *x_dev_poly, unsigned char *y_dev_poly);
extern int factorization_free();
extern unsigned char poly_term_dev_cal(long long term_idx, unsigned char term_coef, unsigned char x_val, unsigned char y_val);
extern int poly_dev_build(unsigned char *poly, unsigned char *poly_dev);
extern int fac_her_lagrange_poly_construct();
extern int fac_ret_poly_construct();
extern int fac_ret_encoding();
extern int fac_dev_init();

#endif
