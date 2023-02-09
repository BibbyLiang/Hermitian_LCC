#ifndef RE_ENCODING_H
#define RE_ENCODING_H

extern unsigned char lag_poly[GF_FIELD][MAX_POLY_TERM_SIZE];
extern unsigned char ret_poly[MAX_POLY_TERM_SIZE];
extern unsigned char ret_cwd_poly[CODEWORD_LEN];
extern unsigned char keep_sym[CODEWORD_LEN];
extern long long keep_flag[CODEWORD_LEN];
extern unsigned char ret_est_msg[MESSAGE_LEN];
extern unsigned char ret_trans_cwd[CODEWORD_LEN];
extern unsigned char ret_est_msg[MESSAGE_LEN];
extern unsigned char ret_est_cwd[CODEWORD_LEN];
extern unsigned char v_poly[MAX_POLY_TERM_SIZE];

extern int her_lagrange_poly_construct();
extern int ret_poly_construct();
extern int ret_encoding();
extern int ret_trans();
extern int ret_cwd_recover();
extern int v_poly_construct();

#endif