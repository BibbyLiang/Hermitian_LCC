#ifndef LCC_DECODING_H
#define LCC_DECODING_H

#include "cfg_decoding.h"

extern float chnl_rel_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
extern float chnl_rel_order[CODEWORD_LEN];
extern long long chnl_rel_order_idx[CODEWORD_LEN];
extern long long chnl_rel_max_id[CODEWORD_LEN];
extern long long chnl_rel_scd_id[CODEWORD_LEN];
extern long long tst_vct_num;

extern unsigned char poly_eva_x_y(unsigned char *poly_x_y, unsigned char x_val, unsigned char y_val);
extern int chnl_rel_cal(float **input_seq, long long input_len);
extern int tst_vct_num_cal();
extern int tst_vct_init();
extern int tst_vct_exit();
extern int tst_vct_form();
extern int term_degree_table_init();
extern long long poly_degree_cal(unsigned char *poly);
extern long long poly_z_degree_get(unsigned char *poly);
extern int koetter_interpolation_hermitian();
extern int poly_init();
extern int poly_normal_update(unsigned char *poly_tmp, unsigned char *poly_update, unsigned char *poly_min, unsigned char hs_dev_min, unsigned char hs_dev_self);
extern int poly_min_update(unsigned char *poly_tmp, unsigned char *poly_min, unsigned char x_point_val, unsigned char hs_dev_min);
extern unsigned char hasse_dev_cal(long long point_idx, long long poly_idx, unsigned char test_sym);
extern int poly_dev_test(unsigned char *test_poly_seq);
extern int g_poly_trans(unsigned char *poly, unsigned char *g_poly, long long pole_idx, long long z_degree);
extern long long lt_get(unsigned char *poly);
extern long long chien_search_x_y(unsigned char *poly);
extern unsigned char recur_root(unsigned char *poly, unsigned char *g_poly);
extern int recur_poly_update(unsigned char *poly, long long pole_idx, unsigned char root_prev_coef);
extern int factorization_recur(unsigned char *poly, unsigned char *est_msg);
extern int factorization_free();

#endif