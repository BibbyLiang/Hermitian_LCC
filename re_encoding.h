#ifndef RE_ENCODING_H
#define RE_ENCODING_H

extern unsigned char lag_poly[CODEWORD_LEN][MAX_POLY_TERM_SIZE];
extern unsigned char ret_poly[MAX_POLY_TERM_SIZE];
extern unsigned char ret_cwd_poly[CODEWORD_LEN];
extern unsigned char keep_sym[CODEWORD_LEN];
extern long long keep_flag[CODEWORD_LEN];
extern long long keep_cnt;
extern unsigned char ret_est_msg[MESSAGE_LEN];
extern unsigned char ret_trans_cwd[CODEWORD_LEN];
extern unsigned char ret_est_msg[MESSAGE_LEN];
extern unsigned char ret_est_cwd[CODEWORD_LEN];
extern unsigned char v_poly[MAX_POLY_TERM_SIZE];
extern unsigned char v_val[CODEWORD_LEN];
extern unsigned char t_val[GF_FIELD][GF_FIELD][GF_FIELD];
#if (1 == CFG_Y_RET_STORE)
extern unsigned char y_t_val[CODEWORD_LEN][CODEWORD_LEN];
#endif

extern int her_lagrange_poly_construct();
extern int ret_poly_construct();
extern int ret_encoding();
extern int ret_trans();
extern int ret_cwd_recover();
extern int v_poly_construct();
extern int keep_position_set(long long *keep_poition);
extern unsigned char ret_fac_free_dev(long long sym_idx, unsigned char dev_flag, unsigned char dev_poly_v);
extern int re_encoding_transform();
extern int ret_fac_free(unsigned char *q0_poly, unsigned char *q1_poly, unsigned char *v_poly);
extern int ret_t_val_cal();
extern int era_cwd_gen();
extern int test_poly_dev();

#endif