#ifndef LCC_DECODING_H
#define LCC_DECODING_H

#include "cfg_decoding.h"

extern FILE *frc;
extern char log_name[255];

extern clock_t start, stop;
extern float runtime;
extern long long err_cnt;

extern float chnl_rel_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
extern float chnl_rel_order[CODEWORD_LEN];
extern long long chnl_rel_order_idx[CODEWORD_LEN];
extern long long chnl_rel_max_id[CODEWORD_LEN];
extern long long chnl_rel_scd_id[CODEWORD_LEN];
extern long long tst_vct_num;
extern unsigned char **tst_vct;
extern unsigned char **ret_tst_vct;
extern long long *tv_err;
extern unsigned char **tv_est_msg;
extern unsigned char **tv_est_cwd;
extern unsigned char *tv_dec_output_flag;

extern long long dev_1_avg_cnt;
extern long long dev_2_avg_cnt;
extern long long dev_out_cnt;

extern int chnl_rel_cal(float **input_seq, long long input_len);
extern int chnl_rel_order_pair_adj();
extern int tst_vct_num_cal();
extern int tst_vct_init();
extern int tst_vct_exit();
extern int tst_vct_form();
extern int koetter_interpolation_hermitian();
extern int recover_lag_poly_init();
extern int cwd2msg(unsigned char *cwd, unsigned char *msg);
extern long long check_cwd_dimension(unsigned char *cwd);
extern int check_result_cwd(unsigned char *cwd, unsigned char *est_cwd, long long tv_idx);
extern int check_result_msg(unsigned char *msg, unsigned char *est_msg, long long tv_idx);
extern int koetter_intp_her_lcc(unsigned char *test_poly_seq,
									long long tv_idx,
                                	long long intp_start_idx,
                                	long long intp_end_idx,
                                	unsigned char input_poly[][MAX_POLY_TERM_SIZE],
                                	unsigned char *output_poly);
extern int her_fac(unsigned char *poly,
	                 unsigned char *est_msg,
	                 unsigned char *est_cwd);
extern int her_fac_free(unsigned char *poly, unsigned char *est_cwd);
extern int kot_node_init();
extern int kot_node_clear();
extern int kot_node_exit();
extern int kot_node_save(unsigned char poly[][MAX_POLY_TERM_SIZE], long long layer_idx, long long place_idx);
extern long long kot_node_load(unsigned char poly[][MAX_POLY_TERM_SIZE], long long tv_idx);
extern int her_cwd_check(unsigned char *cwd);
extern int MLcriterion(unsigned char est_cwd[], unsigned char recv_cwd[]);
extern int ReTMLcriterion(unsigned char est_cwd[], unsigned char recv_cwd[]);
extern int her_cmm_intp();
extern int her_ucm_proc(long long tv_idx);
extern int check_result_tv(unsigned char *tv, unsigned char *est_cwd, unsigned char *est_msg);
extern int her_lcc_check_result();
extern int her_lcc_dec();

#endif