#ifndef LCC_DECODING_H
#define LCC_DECODING_H

#include "cfg_decoding.h"

extern clock_t start, stop;
extern float runtime;
extern long long err_cnt;

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
extern int koetter_interpolation_hermitian();
extern int recover_lag_poly_init();
extern int cwd2msg(unsigned char *cwd, unsigned char *msg);
extern int check_result_cwd(unsigned char *cwd, unsigned char *est_cwd);
extern int check_result_msg(unsigned char *msg, unsigned char *est_msg);
extern int koetter_intp_her_lcc(unsigned char *test_poly_seq,
                                	long long intp_start_idx,
                                	unsigned char input_poly[][MAX_POLY_TERM_SIZE],
                                	unsigned char *output_poly);
extern int her_fac(unsigned char *poly,
	                 unsigned char *est_msg,
	                 unsigned char *est_cwd);
extern int her_fac_free(unsigned char *poly, unsigned char *est_cwd);

#endif