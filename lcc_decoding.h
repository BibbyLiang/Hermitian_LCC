#ifndef LCC_DECODING_H
#define LCC_DECODING_H

#include "cfg_decoding.h"

extern float chnl_rel_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
extern float chnl_rel_order[CODEWORD_LEN];
extern long long chnl_rel_order_idx[CODEWORD_LEN];
extern long long chnl_rel_max_id[CODEWORD_LEN];
extern long long chnl_rel_scd_id[CODEWORD_LEN];
extern long long tst_vct_num;

extern int chnl_rel_cal(float **input_seq, long long input_len);
extern int tst_vct_num_cal();
extern int tst_vct_init();
extern int tst_vct_exit();
extern int tst_vct_form();
extern int term_degree_table_init();
extern int koetter_interpolation_hermitian();
extern int poly_init();

#endif