#ifndef ENCODING_H
#define ENCODING_H

#include "cfg_decoding.h"

extern unsigned char af_pnt[CODEWORD_LEN][2];
extern long long pole_basis_pow[MESSAGE_LEN][2];
extern unsigned char msg_poly[MESSAGE_LEN];
extern unsigned char cwd_poly[CODEWORD_LEN];
extern unsigned char recv_poly[CODEWORD_LEN];
extern unsigned char est_msg_poly[MESSAGE_LEN];
extern unsigned char est_cwd_poly[CODEWORD_LEN];
extern long long pow_trans_order[CODEWORD_LEN];

extern int affine_points_cal();
extern int pole_basis_cal();
extern int rnd_msg_gen();
extern int her_encoding(unsigned char *msg, unsigned char *cwd);
extern int her_convert(unsigned char *poly);
extern unsigned char her_degree_check(unsigned char *poly);

#endif