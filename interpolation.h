#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "cfg_decoding.h"

extern long long x_term_degree_table[MAX_POLY_TERM_SIZE];
extern long long y_term_degree_table[MAX_POLY_TERM_SIZE];
extern long long z_term_degree_table[MAX_POLY_TERM_SIZE];

extern unsigned char intp_poly_coef[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
extern unsigned char intp_poly_tmp[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
extern unsigned char q0_poly_coef[MAX_POLY_TERM_SIZE], q1_poly_coef[MAX_POLY_TERM_SIZE];
extern unsigned char intp_poly_degree[KOT_INTP_POLY_NUM];
extern unsigned char min_intp_poly[MAX_POLY_TERM_SIZE];

extern unsigned char poly_eva_x_y(unsigned char *poly_x_y, unsigned char x_val, unsigned char y_val);
extern int term_degree_table_init();
extern long long poly_degree_cal(unsigned char *poly);
extern long long poly_z_degree_get(unsigned char *poly);
extern int poly_init();
extern int poly_normal_update(unsigned char *poly_tmp, unsigned char *poly_update, unsigned char *poly_min, unsigned char hs_dev_min, unsigned char hs_dev_self);
extern int poly_min_update(unsigned char *poly_tmp, unsigned char *poly_min, unsigned char x_point_val, unsigned char hs_dev_min);
extern unsigned char hasse_dev_cal(long long point_idx, long long poly_idx, unsigned char test_sym);
extern int poly_dev_test(unsigned char *test_poly_seq);

#endif