#ifndef GF_CAL_H
#define GF_CAL_H

#include "cfg_decoding.h"

extern unsigned char power_polynomial_table[GF_FIELD][2];

#if (1 == GF_CAL_COUNT)
extern unsigned char cnt_switch;
extern long long add_cnt;
extern long long mul_cnt;
extern long long div_cnt;
extern long long real_cbm_cnt;
extern long long real_mul_ff_cnt;
extern long long pow_cnt;
extern long long add_cnt_prev;
extern long long mul_cnt_prev;
extern long long div_cnt_prev;
extern long long real_cbm_cnt_prev;
extern long long real_mul_ff_cnt_prev;
extern long long pow_cnt_prev;
#endif

extern unsigned char gf_pow2poly(unsigned char val_in_pow);
extern unsigned char gf_poly2pow(unsigned char val_in_poly);
extern unsigned char gf_location(unsigned char val);
extern unsigned char gf_add(unsigned char a, unsigned char b);
extern unsigned char gf_multp(unsigned char a, unsigned char b);
extern unsigned char gf_div(unsigned char a, unsigned char b);
extern unsigned char gf_mod_single_term(unsigned char a, unsigned char b);
extern long long gf_degree(unsigned char* a, long long len_a);
extern unsigned char gf_div_q_r(unsigned char* dividend, long long len_dividend,
								   unsigned char* divisor, long long len_divisor,
								   unsigned char* quotien, long long len_quotien,
								   unsigned char* remainder, long long len_remainder);
extern unsigned char gf_multp_poly(unsigned char* a, long long len_a,
									   unsigned char* b, long long len_b,
									   unsigned char* product, long long len_product);

extern int gf_multp_poly_hw(unsigned char* a, unsigned char len_a,
				 				  unsigned char* b, unsigned char len_b,
				 				  unsigned char* product, unsigned char len_product);
extern long long real_combine(long long n, long long k);
extern unsigned char gf_real_mutp_ff(long long n, unsigned char ff);
unsigned char gf_pow_cal(unsigned char ff, long long n);
extern unsigned char phase_trans(unsigned char phase);
#if (1 == GF_CAL_COUNT)
extern int gf_count_hist(long long err_cnt);
#endif
extern void gf_count_switch(unsigned char count_switch);
extern void BubbleSort4(float *A, int len, long long *A_idx);
extern void BubbleSort5(float *A, int len, long long *A_idx);
extern void gf_count_reset();
extern long long term_search(long long x_degree, long long y_degree, long long z_degree);

#endif
