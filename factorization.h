#ifndef FACTORIZTION_H
#define FACTORIZTION_H

#include "cfg_decoding.h"

extern int g_poly_trans(unsigned char *poly, unsigned char *g_poly, long long pole_idx, long long z_degree);
extern long long lt_get(unsigned char *poly);
extern long long chien_search_x_y(unsigned char *poly);
extern unsigned char recur_root(unsigned char *poly, unsigned char *g_poly);
extern int recur_poly_update(unsigned char *poly, long long pole_idx, unsigned char root_prev_coef);
extern int factorization_recur(unsigned char *poly, unsigned char *est_msg);
extern int poly_dev_cal(unsigned char *poly, unsigned char *x_dev_poly, unsigned char *y_dev_poly);
extern int factorization_free();

#endif