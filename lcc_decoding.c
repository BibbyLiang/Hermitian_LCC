#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "time.h"
#include "debug_info.h"
#include "math.h"
#include "gf_cal.h"
#include "encoding.h"
#include "channel.h"
#include "mod.h"
#include "interpolation.h"
#include "factorization.h"
#include "re_encoding.h"
#include "br.h"
#include "sys_gen.h"
#include "lcc_decoding.h"

#define PI	3.1415926

FILE *frc;//use to output log
char log_name[255];

clock_t start, stop;
float runtime;
long long err_cnt = 0;

/*col(locator, 0xff~GF_FIELD-2)-row(mesg, CODEWORD_LEN), same as matlab, contrary to most papers*/
float chnl_rel_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
float chnl_rel_order[CODEWORD_LEN];
long long chnl_rel_order_idx[CODEWORD_LEN];
long long chnl_rel_order_pair_idx[CODEWORD_LEN];
long long chnl_rel_max_id[CODEWORD_LEN];
long long chnl_rel_scd_id[CODEWORD_LEN];

long long tst_vct_num = 0;
unsigned char **tst_vct;
unsigned char **ret_tst_vct;
unsigned char **tv_est_msg;
unsigned char **tv_est_cwd;
long long *tv_err;
unsigned char *tv_dec_output_flag;
float *tst_vct_rel;
long long *tst_vct_seq;

unsigned char recover_lag_poly[CODEWORD_LEN][MAX_POLY_TERM_SIZE];
unsigned char recover_lag_flag[CODEWORD_LEN];

unsigned char ***node_poly;
unsigned char *node_tore_flag;
long long node_num = 0;
unsigned char cmm_intp_poly[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
long long cmm_intp_poly_degree[KOT_INTP_POLY_NUM];
long long intp_cnt = 0;

long long dev_1_avg_cnt = 0;
long long dev_2_avg_cnt = 0;
long long dev_out_cnt = 0;

long long new_erasure_val = (CODEWORD_LEN - (((MESSAGE_LEN + (GF_Q) * (GF_Q - 1) / 2) / GF_Q + 1) * GF_Q));

long long sys_ret_ok_cnt = 0;
long long lag_ret_ok_cnt = 0;

int chnl_rel_seq_order()
{
	long long i = 0, j = 0;
	float tmp = 0, max_val = 0, scd_val = 0;
	long long max_idx = 0, scd_idx = 0;
	float chnl_rel[CODEWORD_LEN];

	memset(chnl_rel_order, 0, sizeof(float) * CODEWORD_LEN);
	memset(chnl_rel, 0, sizeof(float) * CODEWORD_LEN);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		chnl_rel_order_idx[i] = i;
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		max_val = 0;
		scd_val = 0;
	
		for(j = 0; j < (CODEWORD_LEN + 1); j++)
		{
			if(max_val < chnl_rel_matrix[j][i])
			{
				max_val = chnl_rel_matrix[j][i];
				max_idx = j;
			}

			chnl_rel_max_id[i] = max_idx;
		}
		for(j = 0; j < (CODEWORD_LEN + 1); j++)
		{
			if((scd_val < chnl_rel_matrix[j][i])
				&& (max_val > chnl_rel_matrix[j][i]))
			{
				scd_val = chnl_rel_matrix[j][i];
				scd_idx = j;
			}

			chnl_rel_scd_id[i] = scd_idx;
		}

		if(0 == scd_val)
		{
			scd_val = 0.00001;
		}
		chnl_rel[i] = max_val / scd_val;
		DEBUG_NOTICE("chnl_rel_order: %f %f %f %ld %ld\n",
		             max_val,
		             scd_val,
		             chnl_rel[i],
		             chnl_rel_max_id[i],
		             chnl_rel_scd_id[i]);
	}

	memcpy(chnl_rel_order, chnl_rel, sizeof(float) * CODEWORD_LEN);

	BubbleSort4(chnl_rel_order, CODEWORD_LEN, chnl_rel_order_idx);//�ɿ�������
	//BubbleSort5(chnl_rel_order, CODEWORD_LEN, chnl_rel_order_idx);//�ɿ��Ƚ���

#if (1 == TEST_MODE)
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("chnl_rel_order_idx: %ld %f %ld\n",
					 i,
					 chnl_rel_order[i],
		             chnl_rel_order_idx[i]);
	}
#endif

	//memcpy(chnl_rel_order_pair_idx, chnl_rel_order_idx, sizeof(long long) * CODEWORD_LEN);

	return 0;
}

#if (1 == CFG_RET_ETA_OPT)
int chnl_rel_order_pair_adj()
{
	long long i = 0, j = 0;
	long long pair_cnt = 0, sym_cnt = 0;
	long long sym_idx = 0;
	float sym_rel_val = 0;
	long long pair_hist[GF_FIELD];
	memset(pair_hist, 0, sizeof(long long) * GF_FIELD);
	unsigned char tmp_flag = 0;
	unsigned char sym_choose_flag[CODEWORD_LEN];
	memset(sym_choose_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		sym_idx = chnl_rel_order_idx[i];
		for(j = 0; j < GF_FIELD; j++)
		{
			if(af_pnt[sym_idx][0] == power_polynomial_table[j][0])
			{
				if((1 == tmp_flag)
					&& (0 == pair_hist[j]))
				{
					sym_choose_flag[i] = 2;
					break;
				}

				pair_hist[j]++;
				sym_cnt++;
				if(0 == tmp_flag)
				{
					sym_choose_flag[i] = 1;
				}
				else
				{
					sym_choose_flag[i] = 3;
				}
				if(1 == pair_hist[j])
				{
					pair_cnt++;
				}
				break;
			}
		}
		DEBUG_NOTICE("sym_choose: %ld | %ld %ld | %ld | %x\n", i, sym_cnt, pair_cnt, sym_idx, sym_choose_flag[i]);
		
		if((GF_FIELD - (MESSAGE_LEN - GF_Q * (GF_Q - 1) / 2) / GF_Q) <= pair_cnt)
		{
			tmp_flag = 1;
		}
		if(ETA <= sym_cnt)
		{
			break;
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("chnl_rel_order_idx_prv: %ld %f %ld\n",
					 i,
					 chnl_rel_order[i],
					 chnl_rel_order_idx[i]);
	}
#endif

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(3 == sym_choose_flag[i])
		{
			for(j = 0; j < CODEWORD_LEN; j++)
			{
				if(2 == sym_choose_flag[j])
				{
					DEBUG_NOTICE("sym_rel_exc: %ld %ld | %ld %ld\n", i, j, chnl_rel_order_idx[i], chnl_rel_order_idx[j]);
					sym_idx = chnl_rel_order_idx[i];
					sym_rel_val = chnl_rel_order[i];
					chnl_rel_order_idx[i] = chnl_rel_order_idx[j];
					chnl_rel_order[i] = chnl_rel_order[j];
					chnl_rel_order_idx[j] = sym_idx;
					chnl_rel_order[j] = sym_rel_val;
					sym_choose_flag[i] = 2;
					sym_choose_flag[j] = 3;
					break;
				}
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("chnl_rel_order_idx_adj: %ld %f %ld\n",
					 i,
					 chnl_rel_order[i],
		             chnl_rel_order_idx[i]);
	}
#endif

	return 0;
}
#endif

int chnl_rel_cal(float **input_seq,
				    long long input_len)
{
	long long i = 0, j = 0, k = 0;
	float n0, d0 = 0, d1 = 0, temp = 0;
	unsigned char tmp_bit = 0;
	/*for BPSK*/
	float map[input_len][2];

	for(i = 0; i < (CODEWORD_LEN + 1); i++)
	{
		memset(chnl_rel_matrix[i], 0, sizeof(float) * CODEWORD_LEN);
	}

	n0 = 1 / ((float)MESSAGE_LEN / (float)CODEWORD_LEN) / (pow(10, eb2n0 / 10) / 2);

	long long seq_idx = 0;
	float base_symbol = 0.0;
	float d_sum = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		d_sum = 0;

		for(j = 0; j < GF_FIELD; j++)
		{
			d0 = 0;
			d1 = 0;

			for(k = 0; k < GF_Q; k++)
			{
				seq_idx = i * GF_Q + k;
				base_symbol = (float)((power_polynomial_table[j][1] >> k) & 0x1);
				base_symbol = 1 - 2 * base_symbol;
				
				d0 = (input_seq[seq_idx][0] - base_symbol) * (input_seq[seq_idx][0] - base_symbol)
					  + (input_seq[seq_idx][1] - (0.0)) * (input_seq[seq_idx][1] - (0.0));
				if(d0 >= d1)
				{
					d1 = d0;
				}

			}
			
			if(d_sum < d1)
			{
				d_sum = d1;
			}
		}
		DEBUG_INFO("d0: %ld | %f %f| %f\n", i, d0, d1, d_sum);
	}

	for(i = 0; i < input_len; i++)
	{
		d0 = (input_seq[i][0] - (1.0)) * (input_seq[i][0] - (1.0))
				+ (input_seq[i][1] - (0.0)) * (input_seq[i][1] - (0.0));
		d1 = (input_seq[i][0] - (-1.0)) * (input_seq[i][0] - (-1.0))
				+ (input_seq[i][1] - (0.0)) * (input_seq[i][1] - (0.0));

#if 1/*for float precision problem*/
		d0 = d0 / 1e4;
		d1 = d1 / 1e4;
#endif

		/*for BPSK*/
		map[i][0] = 1 / (PI * n0) * exp((-d0) / n0);
        map[i][1] = 1 / (PI * n0) * exp((-d1) / n0);
        DEBUG_INFO("map: %d %f | %f %f | %f %f\n",
        		   i,
        		   n0,
        		   d0,
        		   d1,
        		   map[i][0],
        		   map[i][1]);
        if((0 == map[i][0])
        	|| (0 == map[i][1]))
        {
        	DEBUG_INFO("map: %d %f | %f %f | %f %f\n",
        			  i,
        			  n0,
        			  map[i][0],
        			  map[i][0],
        			  d0,
        			  d1);
        }
	}

	for(i = 0; i < GF_FIELD; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			chnl_rel_matrix[i][j] = 1;

			for(k = 0; k < GF_Q; k++)
			{
				/*for BPSK*/
				tmp_bit = ((power_polynomial_table[i][1] >> k) & 0x1);
				DEBUG_INFO("chnl_rel_matrix: %d %d | %d | %d | %f %f\n",
		        		   i,
		        		   j,
		        		   k,
		        		   tmp_bit,
		        		   chnl_rel_matrix[i][j],
		        		   map[j * GF_Q + k][tmp_bit]);

				chnl_rel_matrix[i][j] = chnl_rel_matrix[i][j]
									   * map[j * GF_Q + k][tmp_bit];
			}
		}
	}

	for (i = 0; i < CODEWORD_LEN; i++)
    {
		DEBUG_NOTICE("recv: %x\n", recv_poly[i]);

        temp = 0;
        for (j = 0; j < GF_FIELD; j++)
        {
            temp = temp + chnl_rel_matrix[j][i];
        }    

        for (j = 0; j < GF_FIELD; j++)
        {
            chnl_rel_matrix[j][i] = chnl_rel_matrix[j][i] / temp;
            if(chnl_rel_matrix[j][i] != chnl_rel_matrix[j][i])
            {
            	DEBUG_SYS("NAN: %d %d | %f %f | %f %f %f %f %f %f %f %f\n",
            			   i,
            	           j,
            	           temp,
            	           chnl_rel_matrix[j][i],
            	           input_seq[i * GF_Q][0],
            	           input_seq[i * GF_Q + 1][0],
            	           input_seq[i * GF_Q + 2][0],
            	           input_seq[i * GF_Q + 3][0],
            	           input_seq[i * GF_Q + 4][0],
            	           input_seq[i * GF_Q + 5][0],
            	           input_seq[i * GF_Q + 6][0],
            	           input_seq[i * GF_Q + 7][0]);
            }
			DEBUG_NOTICE("chnl_rel: %ld %ld | %f\n", i, power_polynomial_table[j][0], chnl_rel_matrix[j][i]);
        }
    }

#if 0//(1 == TEST_MODE)
	for (i = 0; i < CODEWORD_LEN; i++)
    {
        for (j = 0; j < (CODEWORD_LEN + 1); j++)
        {
            chnl_rel_matrix[j][i] = 0;
        }
    }
    chnl_rel_matrix[1][1] = 1;
    chnl_rel_matrix[1][7] = 1;
    chnl_rel_matrix[3][6] = 1;
    chnl_rel_matrix[2][0] = 1;

	chnl_rel_matrix[0][5] = 0.9;
	chnl_rel_matrix[1][5] = 0.1;

	chnl_rel_matrix[2][2] = 0.8;
	chnl_rel_matrix[0][2] = 0.2;

	chnl_rel_matrix[1][3] = 0.7;
	chnl_rel_matrix[0][3] = 0.3;

	chnl_rel_matrix[0][4] = 0.6;
	chnl_rel_matrix[2][4] = 0.4;
#endif

	chnl_rel_seq_order();

#if (1 == CFG_RET_ETA_OPT)
	chnl_rel_order_pair_adj();
#endif

	return 0;
}

int tst_vct_num_cal()
{
	tst_vct_num = (long long)pow(2, ETA);
	DEBUG_NOTICE("tst_vct_num: %ld\n", tst_vct_num);

	return 0;
}

int tst_vct_init()
{
	long long i = 0;

	tst_vct = (unsigned char**)malloc(sizeof(unsigned char*) * tst_vct_num);
	for (i = 0; i < tst_vct_num; i++)
	{
		tst_vct[i] = (unsigned char*)malloc(sizeof(unsigned char) * CODEWORD_LEN);
	}
	for (i = 0; i < tst_vct_num; i++)
	{
		memset(tst_vct[i], 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	}

	ret_tst_vct = (unsigned char**)malloc(sizeof(unsigned char*) * tst_vct_num);
	for (i = 0; i < tst_vct_num; i++)
	{
		ret_tst_vct[i] = (unsigned char*)malloc(sizeof(unsigned char) * CODEWORD_LEN);
	}
	for (i = 0; i < tst_vct_num; i++)
	{
		memset(ret_tst_vct[i], 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	}

	tv_est_msg = (unsigned char**)malloc(sizeof(unsigned char*) * tst_vct_num);
	for (i = 0; i < tst_vct_num; i++)
	{
		tv_est_msg[i] = (unsigned char*)malloc(sizeof(unsigned char) * MESSAGE_LEN);
	}
	for (i = 0; i < tst_vct_num; i++)
	{
		memset(tv_est_msg[i], 0xFF, sizeof(unsigned char) * MESSAGE_LEN);
	}

	tv_est_cwd = (unsigned char**)malloc(sizeof(unsigned char*) * tst_vct_num);
	for (i = 0; i < tst_vct_num; i++)
	{
		tv_est_cwd[i] = (unsigned char*)malloc(sizeof(unsigned char) * CODEWORD_LEN);
	}
	for (i = 0; i < tst_vct_num; i++)
	{
		memset(tv_est_cwd[i], 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	}

	tv_err = (long long*)malloc(sizeof(long long) * tst_vct_num);
	for(i = 0; i < tst_vct_num; i++)
	{
		tv_err[i] = 0;
	}
	
	tv_dec_output_flag = (unsigned char*)malloc(sizeof(unsigned char) * tst_vct_num);
	for(i = 0; i < tst_vct_num; i++)
	{
		tv_dec_output_flag[i] = 0;
	}
	
	tst_vct_rel = (float*)malloc(sizeof(float) * tst_vct_num);
	for(i = 0; i < tst_vct_num; i++)
	{
		tst_vct_rel[i] = 0;
	}
	
	tst_vct_seq = (long long*)malloc(sizeof(long long) * tst_vct_num);
	for(i = 0; i < tst_vct_num; i++)
	{
		tst_vct_seq[i] = i;
	}

	return 0;
}

int tst_vct_exit()
{
	long long i = 0;

	for (i = 0; i < tst_vct_num; i++)
	{
  		free(tst_vct[i]);
		tst_vct[i] = NULL;
  	}
	free(tst_vct);
	tst_vct = NULL;

	for (i = 0; i < tst_vct_num; i++)
	{
  		free(ret_tst_vct[i]);
		ret_tst_vct[i] = NULL;
  	}
	free(ret_tst_vct);
	ret_tst_vct = NULL;
	
	for (i = 0; i < tst_vct_num; i++)
	{
  		free(tv_est_msg[i]);
		tv_est_msg[i] = NULL;
  	}
	free(tv_est_msg);
	tv_est_msg = NULL;
	
	for (i = 0; i < tst_vct_num; i++)
	{
  		free(tv_est_cwd[i]);
		tv_est_cwd[i] = NULL;
  	}
	free(tv_est_cwd);
	tv_est_cwd = NULL;

	free(tv_err);
	tv_err = NULL;
	
	free(tv_dec_output_flag);
	tv_dec_output_flag = NULL;
	
	free(tst_vct_rel);
	tst_vct_rel = NULL;
	
	free(tst_vct_seq);
	tst_vct_seq = NULL;

	return 0;
}

int tst_vct_form()
{
	long long i = 0, j = 0, pickout_place = 0;
	long long cwd_idx = 0;
	long long pickout_sym = 0xFF;

	for(i = 0; i < tst_vct_num; i++)
	{
		tst_vct_rel[i] = 0;
	}	
	for(i = 0; i < tst_vct_num; i++)
	{
		tst_vct_seq[i] = i;
	}

	for(i = 0; i < tst_vct_num; i++)
	{
		tst_vct_rel[i] = 1;
	
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			cwd_idx = chnl_rel_order_idx[j];

			if(ETA < j)
			{
				pickout_place = 0;
#if 0//there may be some problems because of the precision of float			
				pickout_sym = chnl_rel_max_id[cwd_idx];
#else
				if(0xFF == recv_poly[cwd_idx])
				{
					pickout_sym = 0;
				}
				else
				{
					pickout_sym = recv_poly[cwd_idx] + 1;
				}
#endif				
			}
			else
			{
				pickout_place = (i >> j) & 0x1;
				if(0 == pickout_place)
				{
					pickout_sym = chnl_rel_max_id[cwd_idx];
				}
				else
				{
					pickout_sym = chnl_rel_scd_id[cwd_idx];
				}
				
				tst_vct_rel[i] = tst_vct_rel[i] * chnl_rel_matrix[pickout_sym][cwd_idx];
			}
			DEBUG_NOTICE("tst_vct_cal: %ld %ld | %ld %ld | %ld\n",
			             i,
			             j,
			             cwd_idx,
			             pickout_place,
			             pickout_sym);
			tst_vct[i][cwd_idx] = power_polynomial_table[pickout_sym][0];
		}
		
		DEBUG_NOTICE("tst_vct_rel: %d | %f\n",
		              i,
		              tst_vct_rel[i]);
	}

	/*for test*/
	//memcpy(tst_vct[1], tst_vct[0], sizeof(unsigned char) * CODEWORD_LEN);
	//memcpy(tst_vct[2], tst_vct[0], sizeof(unsigned char) * CODEWORD_LEN);
	//memcpy(tst_vct[3], tst_vct[0], sizeof(unsigned char) * CODEWORD_LEN);

	for(i = 0; i < tst_vct_num; i++)
	{
		tv_err[i] = 0;
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if(tst_vct[i][j] != cwd_poly[j])
			{
				DEBUG_NOTICE("tv_err_idx: %ld %ld | %x %x\n", i, j, tst_vct[i][j], cwd_poly[j]);
				tv_err[i]++;
			}
		}
		DEBUG_NOTICE("tv_err: %ld | %ld\n", i, tv_err[i]);
		
		tv_dec_output_flag[i] = 0;
	}

#if (1 == CFG_DEBUG_NOTICE)
	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			DEBUG_NOTICE("tst_vct: %ld %ld | %x\n",
			             i,
			             j,
			             tst_vct[i][j]);
			if(tst_vct[i][j] != recv_poly[j])
			{
				DEBUG_NOTICE("tst_vct_diff: %ld %ld | %x %x\n",
				             i,
				             j,
				             tst_vct[i][j],
				             recv_poly[j]);
			}
		}
	}
#endif

#if 0
	tst_vct[0][0] = 0xff;
	tst_vct[0][1] = 0xff;
	tst_vct[0][2] = 0xff;
	tst_vct[0][3] = 0xff;
	tst_vct[0][4] = 0xff;
	tst_vct[0][5] = 0xff;
	tst_vct[0][6] = 0xff;
	tst_vct[0][7] = 0xff;
	tst_vct[0][8] = 0xff;
	tst_vct[0][9] = 0xff;
	tst_vct[0][10] = 0xff;
	tst_vct[0][11] = 0xff;
	tst_vct[0][12] = 0xff;
	tst_vct[0][13] = 0xff;
	tst_vct[0][14] = 0xff;
	tst_vct[0][15] = 0xff;
	tst_vct[0][16] = 0xff;
	tst_vct[0][17] = 0xff;
	tst_vct[0][18] = 0xff;
	tst_vct[0][19] = 0xff;
	tst_vct[0][20] = 0xff;
	tst_vct[0][21] = 0xff;
	tst_vct[0][22] = 0xff;
	tst_vct[0][23] = 0xff;
	tst_vct[0][24] = 0xff;
	tst_vct[0][25] = 0xff;
	tst_vct[0][26] = 0xff;
	tst_vct[0][27] = 0xff;
	tst_vct[0][28] = 0xff;
	tst_vct[0][29] = 0xff;
	tst_vct[0][30] = 0xff;
	tst_vct[0][31] = 0xff;
	tst_vct[0][32] = 0x8;
	tst_vct[0][33] = 0x2;
	tst_vct[0][34] = 0x0;
	tst_vct[0][35] = 0x4;
	tst_vct[0][36] = 0x0;
	tst_vct[0][37] = 0x5;
	tst_vct[0][38] = 0x1;
	tst_vct[0][39] = 0x5;
	tst_vct[0][40] = 0x5;
	tst_vct[0][41] = 0xe;
	tst_vct[0][42] = 0x0;
	tst_vct[0][43] = 0x8;
	tst_vct[0][44] = 0xb;
	tst_vct[0][45] = 0xff;
	tst_vct[0][46] = 0x8;
	tst_vct[0][47] = 0x8;
	tst_vct[0][48] = 0x2;
	tst_vct[0][49] = 0x6;
	tst_vct[0][50] = 0xff;
	tst_vct[0][51] = 0xff;
	tst_vct[0][52] = 0x4;
	tst_vct[0][53] = 0x3;
	tst_vct[0][54] = 0x9;
	tst_vct[0][55] = 0x0;
	tst_vct[0][56] = 0x0;
	tst_vct[0][57] = 0xd;
	tst_vct[0][58] = 0x6;
	tst_vct[0][59] = 0x1;
	tst_vct[0][60] = 0x2;
	tst_vct[0][61] = 0x1;
	tst_vct[0][62] = 0xff;
	tst_vct[0][63] = 0x1;
#endif

	BubbleSort5(tst_vct_rel, tst_vct_num, tst_vct_seq);
	for(i = 0; i < tst_vct_num; i++)
	{
		DEBUG_INFO("all_rel_sort: %ld | %f | %ld\n",
		           i,
		           tst_vct_rel[i],
		           tst_vct_seq[i]);
	}

	return 0;
}

int koetter_interpolation_hermitian()
{
	poly_init();
#if (0 == CFG_RET)
#if 1
	poly_dev_test(recv_poly);
#else
#if 0
	koetter_intp_her_lcc(recv_poly,
						 0,
						 0,
						 CODEWORD_LEN - 1,
						 intp_poly_coef,
						 min_intp_poly);
	her_fac(min_intp_poly, est_msg_poly, est_cwd_poly);
#endif
	poly_init();
	koetter_intp_her_lcc(recv_poly,
						 0,
						 0,
						 CODEWORD_LEN - 1 - ETA,
						 intp_poly_coef,
						 min_intp_poly);
#if 0
	koetter_intp_her_lcc(recv_poly,
						 0,
						 CODEWORD_LEN - ETA,
						 CODEWORD_LEN - 1,
						 intp_poly_coef,
						 min_intp_poly);
#endif
#endif
#else
	poly_dev_test(ret_trans_cwd);
#endif

#if (1 == CFG_FAC_FREE)

#if (1 == CFG_RET)
	//ret_fac_free(q0_poly_coef, q1_poly_coef, v_poly);
	//ret_cwd_recover();
	//cwd2msg(est_cwd_poly, est_msg_poly);

	her_fac(min_intp_poly, est_msg_poly, est_cwd_poly);
#else
	//factorization_free();

	her_fac(min_intp_poly, est_msg_poly, est_cwd_poly);
#endif

#else

#if (1 == CFG_RET)
	//factorization_recur(min_intp_poly, ret_est_msg);
	//ret_cwd_recover();

	her_fac(min_intp_poly, est_msg_poly, est_cwd_poly);
#else
	//factorization_recur(min_intp_poly, est_msg_poly);
	
	her_fac(min_intp_poly, est_msg_poly, est_cwd_poly);
#endif

#endif

	return 0;
}

int recover_lag_poly_init()
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char x_poly[MAX_POLY_TERM_SIZE], y_poly[MAX_POLY_TERM_SIZE];
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE], tmp_poly_y[MAX_POLY_TERM_SIZE];
	unsigned char y_her_curve_find_flag = 0;
	unsigned char x_div = 0x0, y_div = 0x0;
	unsigned char check_val = 0xFF;

	unsigned ret_x_flag = 0;
	long long lag_cnt = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		memset(recover_lag_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	memset(recover_lag_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);

	memset(erasure_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
#if 1
	memset(erasure_flag, 1, sizeof(unsigned char) * new_erasure_val);
	DEBUG_SYS("new_erasure: %ld\n", new_erasure_val);
#endif

	/*cal lag. poly., this can be opt.*/
	/*part of lag. poly. can be cal. and stored in ReT*/
	for(i = 0; i < CODEWORD_LEN; i++)//af. point
	{
		if(1 == erasure_flag[i])
		{
			continue;
		}
		lag_cnt++;
		recover_lag_flag[i] = 1;

		memset(x_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		x_poly[0] = 0x0;
		memset(y_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		y_poly[0] = 0x0;
		x_div = 0x0;
		y_div = 0x0;

		for(j = 0; j < GF_FIELD; j++)//locator
		{
			ret_x_flag = 0;
			for(k = 0; k < CODEWORD_LEN; k++)
			{
#if 1
				if((0 == erasure_flag[k])
					&& (af_pnt[k][0] == power_polynomial_table[j][0]))
#else
				if(af_pnt[k][0] == power_polynomial_table[j][0])
#endif
				{
					ret_x_flag = 1;
					break;
				}
			}

			if(af_pnt[i][0] == power_polynomial_table[j][0])
			{				             
				//continue;
			}
			else if(0 == ret_x_flag)
			{
				DEBUG_NOTICE("ret_x_flag_0: %ld | %x | %x %x\n",
							 i,
							 power_polynomial_table[j][0],
							 af_pnt[i][0],
							 af_pnt[i][1]);
			}
			else
			{
				memset(tmp_poly_x, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

				/*this can be cal. offline*/
				x_div = gf_multp(x_div, gf_add(af_pnt[i][0], power_polynomial_table[j][0]));

				for(k = (MAX_POLY_TERM_SIZE - 1); k >= 0; k--)
				{
					if(0xFF != x_poly[k])
					{
						tmp_poly_x[k] = gf_multp(x_poly[k], power_polynomial_table[j][0]);
						for(l = 0; l < MAX_POLY_TERM_SIZE; l++)
						{
							if((x_term_degree_table[l] == (x_term_degree_table[k] + 1))
								&& (y_term_degree_table[l] == y_term_degree_table[k])
								&& (z_term_degree_table[l] == z_term_degree_table[k]))
							{
								x_poly[l] = x_poly[k];
								x_poly[k] = 0xFF;
								break;
							}
						}
					}
				}
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if(0xFF == tmp_poly_x[k])
					{
						continue;
					}
					else if((0xFF == x_poly[k])
							&& (0xFF != tmp_poly_x[k]))
					{
						x_poly[k] = tmp_poly_x[k];
					}
					else
					{
						x_poly[k] = gf_add(x_poly[k], tmp_poly_x[k]);
					}
				}
			}

			if(af_pnt[i][1] == power_polynomial_table[j][0])
			{
				//continue;
			}
			else
			{
				memset(tmp_poly_y, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				y_her_curve_find_flag = 0;
				for(k = 0; k < CODEWORD_LEN; k++)
				{
					if((af_pnt[k][1] == power_polynomial_table[j][0])
						&& (af_pnt[k][0] == af_pnt[i][0]))
					{
						y_her_curve_find_flag = 1;
						break;
					}
				}

				/*this can be cal. offline*/
				y_div = gf_multp(y_div, gf_add(af_pnt[i][1], power_polynomial_table[j][0]));

				for(k = (MAX_POLY_TERM_SIZE - 1); k >= 0; k--)
				{
					if((0xFF != y_poly[k])
						&& (1 == y_her_curve_find_flag))//notice this
					{
						tmp_poly_y[k] = gf_multp(y_poly[k], power_polynomial_table[j][0]);
						for(l = 0; l < MAX_POLY_TERM_SIZE; l++)
						{
							if((x_term_degree_table[l] == x_term_degree_table[k])
								&& (y_term_degree_table[l] == (y_term_degree_table[k] + 1))
								&& (z_term_degree_table[l] == z_term_degree_table[k]))
							{
								y_poly[l] = y_poly[k];
								y_poly[k] = 0xFF;
								break;
							}
						}
					}
				}
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if(0xFF == tmp_poly_y[k])
					{
						continue;
					}
					else if((0xFF == y_poly[k])
							&& (0xFF != tmp_poly_y[k]))
					{
						y_poly[k] = tmp_poly_y[k];
					}
					else
					{
						y_poly[k] = gf_add(y_poly[k], tmp_poly_y[k]);
					}
				}
			}
		}

		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)//for y
		{
			if(0xFF != y_poly[j])
			{
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)//for x
				{
					if(0xFF != x_poly[k])
					{
						for(l = 0; l < MAX_POLY_TERM_SIZE; l++)//for new term
						{
							if((x_term_degree_table[l] == x_term_degree_table[k])
								&& (y_term_degree_table[l] == (y_term_degree_table[k] + y_term_degree_table[j]))
								&& (z_term_degree_table[l] == z_term_degree_table[k]))
							{
								recover_lag_poly[i][l] = gf_add(recover_lag_poly[i][l],
														gf_multp(x_poly[k], y_poly[j]));
								recover_lag_poly[i][l] = gf_div(recover_lag_poly[i][l], x_div);
								recover_lag_poly[i][l] = gf_div(recover_lag_poly[i][l], y_div);
							}
						}
					}
				}
			}
		}

		her_convert(recover_lag_poly[i]);

#if (1 == TEST_MODE)
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != recover_lag_poly[i][j])
			{
				DEBUG_NOTICE("lag_poly: %ld | %ld %ld %ld | %x\n",
							  i,
							  x_term_degree_table[j],
							  y_term_degree_table[j],
							  z_term_degree_table[j],
							  recover_lag_poly[i][j]);
			}
		}
#endif		
	}

	return 0;
}

int cwd2msg(unsigned char *cwd, unsigned char *msg)
{
	long long i = 0, j = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(ret_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
#if 1
		if(0 == recover_lag_flag[i])
		{
			DEBUG_NOTICE("recover_erasure_skip: %ld\n", i);
			continue;
		}
		else
#endif		
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(0xFF != recover_lag_poly[i][j])
				{
					if((0xFF == cwd[i])
						|| (0xFF == recover_lag_poly[i][j]))
					{
						tmp_poly[j] = 0xFF;
					}
					else if(0x0 == cwd[i])
					{
						tmp_poly[j] = recover_lag_poly[i][j];
					}
					else if(0x0 == recover_lag_poly[i][j])
					{
						tmp_poly[j] = cwd[i];
					}
					else
					{
						tmp_poly[j] = gf_multp(cwd[i], recover_lag_poly[i][j]);
					}
					if((0xFF == ret_poly[j])
						&& (0xFF == tmp_poly[j]))
					{
						ret_poly[j] = 0xFF;
					}
					else if(0xFF == ret_poly[j])
					{
						ret_poly[j] = tmp_poly[j];
					}
					else if(0xFF == tmp_poly[j])
					{
						ret_poly[j] = ret_poly[j];
					}
					else
					{
						ret_poly[j] = gf_add(ret_poly[j], tmp_poly[j]);
					}

					//her_convert(ret_poly);

					DEBUG_NOTICE("ret_poly_cal: %ld | %x %x | %ld %ld %ld | %x\n",
								  i,
								  cwd[i],
								  recover_lag_poly[i][j],
					              x_term_degree_table[j],
					              y_term_degree_table[j],
					              z_term_degree_table[j],
					              ret_poly[j]);
				}
			}
		}
	}
	her_convert(ret_poly);

#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != ret_poly[i])
		{
			DEBUG_NOTICE("ret_poly: %ld %ld %ld | %x\n",
			              x_term_degree_table[i],
			              y_term_degree_table[i],
			              z_term_degree_table[i],
			              ret_poly[i]);
		}
	}
#endif

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if((pole_basis_pow[i][0] == x_term_degree_table[j])
				&& (pole_basis_pow[i][1] == y_term_degree_table[j])
				&& (0 == z_term_degree_table[j]))
			{
				msg[i] = ret_poly[j];
				break;
			}
		}
		DEBUG_NOTICE("cwd2msg: %ld | %x %x\n", i, msg[i], msg_poly[i]);
	}

	return 0;
}

long long check_cwd_dimension(unsigned char *cwd)
{
	long long i = 0, j = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(ret_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	long long cwd_dim = 0, term_degree = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
#if 1
		if(0 == recover_lag_flag[i])
		{
			DEBUG_NOTICE("recover_erasure_skip: %ld\n", i);
			continue;
		}
		else
#endif		
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				/*just check those high part*/
				if(MESSAGE_LEN >=
				   (x_term_degree_table[j] * GF_Q
				    + y_term_degree_table[j] * (GF_Q + 1)
				    - (GF_Q) * (GF_Q - 1) / 2
				    + 1))
				{
					continue;
				}

				if(0xFF != recover_lag_poly[i][j])
				{
					if((0xFF == cwd[i])
						|| (0xFF == recover_lag_poly[i][j]))
					{
						tmp_poly[j] = 0xFF;
					}
					else if(0x0 == cwd[i])
					{
						tmp_poly[j] = recover_lag_poly[i][j];
					}
					else if(0x0 == recover_lag_poly[i][j])
					{
						tmp_poly[j] = cwd[i];
					}
					else
					{
						tmp_poly[j] = gf_multp(cwd[i], recover_lag_poly[i][j]);
					}
					if((0xFF == ret_poly[j])
						&& (0xFF == tmp_poly[j]))
					{
						ret_poly[j] = 0xFF;
					}
					else if(0xFF == ret_poly[j])
					{
						ret_poly[j] = tmp_poly[j];
					}
					else if(0xFF == tmp_poly[j])
					{
						ret_poly[j] = ret_poly[j];
					}
					else
					{
						ret_poly[j] = gf_add(ret_poly[j], tmp_poly[j]);
					}

					her_convert(ret_poly);
#if 1
					DEBUG_NOTICE("ret_poly_cal: %ld | %x %x | %ld %ld %ld | %x\n",
								  i,
								  cwd[i],
								  recover_lag_poly[i][j],
					              x_term_degree_table[j],
					              y_term_degree_table[j],
					              z_term_degree_table[j],
					              ret_poly[j]);
#endif
				}
			}
		}
	}
	her_convert(ret_poly);

#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != ret_poly[i])
		{
			DEBUG_NOTICE("ret_poly: %ld %ld %ld | %x\n",
			              x_term_degree_table[i],
			              y_term_degree_table[i],
			              z_term_degree_table[i],
			              ret_poly[i]);
		}
	}
#endif
	term_degree = poly_degree_cal(ret_poly);
	cwd_dim = term_degree + 1 - (GF_Q) * (GF_Q - 1) / 2;

	return cwd_dim;
}

int check_result_cwd(unsigned char *cwd, unsigned char *est_cwd, long long tv_idx)
{
	long long i = 0;
	long long genus = 0, radius = 0;
	long long cwd_err = 0;
	int val = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(cwd[i] != est_cwd[i])
		{
			DEBUG_IMPORTANT("check_result_cwd: %ld | %x %x\n",
			               i,
			               cwd[i],
			               est_cwd[i]);
			val++;
		}
	}
	if(0 == val)
	{
		DEBUG_IMPORTANT("cwd_ok\n");
		if(GF_Q > min_intp_idx)
		{
			DEBUG_IMPORTANT("min_intp_idx: %ld\n", min_intp_idx);
		}
	}
	else
	{
		//err_cnt++;
		genus = GF_Q * (GF_Q - 1) / 2;
		//radius = (CODEWORD_LEN - MESSAGE_LEN - genus + 1 - 1) / 2;
		/*for GS (m = 1), that is the GS decoding bound*/
		radius = (CODEWORD_LEN - (MESSAGE_LEN + genus - 1) - 2 * genus) / 2;
		if(0 == ((MESSAGE_LEN - genus + 1) % 2))
		{
			radius = radius - 1;
		}
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			//if(cwd_poly[i] != recv_poly[i])
			if(cwd_poly[i] != tst_vct[tv_idx][i])
			{
				cwd_err++;
			}
		}
		if(cwd_err > radius)
		{
			DEBUG_IMPORTANT("err_out_of_radius_cwd: %ld %ld %ld\n",
			          	   cwd_err,
			          	   radius,
			          	   val);
		}
		else
		{
			DEBUG_SYS("cwd decoding err: %ld %ld\n",
			          cwd_err,
			          radius);
#if 1
			for(i = 0; i < MESSAGE_LEN; i++)
			{
				DEBUG_SYS("check_result_msg_err: %ld | %x %x\n",
				          i,
				          msg_poly[i],
				          est_msg_poly[i]);
			}
#endif			
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				DEBUG_SYS("check_result_cwd_err: %ld | %ld | %x %x %x\n",
						  tv_idx,
			              i,
			              cwd[i],
			              est_cwd[i],
			              tst_vct[tv_idx][i]);
			}
			for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
			{
				if(0xFF != min_intp_poly[i])
				{
					DEBUG_SYS("min_intp_poly_err: %ld %ld %ld | %x\n",
					          x_term_degree_table[i],
					          y_term_degree_table[i],
					          z_term_degree_table[i],
					          min_intp_poly[i]);
				}
			}
			for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
			{
				if(0xFF != q0_poly_coef[i])
				{
					DEBUG_SYS("q0_poly_err: %ld %ld %ld | %x\n",
					          x_term_degree_table[i],
					          y_term_degree_table[i],
					          z_term_degree_table[i],
					          q0_poly_coef[i]);
				}
			}
			for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
			{
				if(0xFF != q1_poly_coef[i])
				{
					DEBUG_SYS("q1_poly_err: %ld %ld %ld | %x\n",
					          x_term_degree_table[i],
					          y_term_degree_table[i],
					          z_term_degree_table[i],
					          q1_poly_coef[i]);
				}
			}

			for(i = 0; i < CODEWORD_LEN; i++)
			{
				DEBUG_SYS("keep_poition: %ld | %ld %ld\n", i, keep_cnt, keep_flag[i]);
			}
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				DEBUG_SYS("chnl_rel_order_idx: %ld %f %ld\n",
							 i,
							 chnl_rel_order[i],
				             chnl_rel_order_idx[i]);
			}
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				DEBUG_SYS("chnl_rel_order_idx: %ld %f %ld\n",
							 i,
							 chnl_rel_order[i],
				             chnl_rel_order_idx[i]);
			}
			DEBUG_SYS("Transmission over Channel:\n");
			for(i = 0; i < symbol_num; i++)
			{
				DEBUG_SYS("rx_symbol[%ld][0] = %f;\n", i, recv_seq[i][0]);
				DEBUG_SYS("rx_symbol[%ld][1] = %f;\n", i, recv_seq[i][1]);
			}
			for(i = 0; i < MESSAGE_LEN; i++)
			{
				DEBUG_SYS("msg_poly[%ld] = 0x%x;\n",
				          i,
				          msg_poly[i]);
			}
			DEBUG_SYS("\n");

#if (1 == OUTPUT_LOG)
			frc = fopen(log_name, "a+");
			fprintf(frc, "cwd decoding err: %ld %ld\n",
						  cwd_err,
						  radius);
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				fprintf(frc, "check_result_cwd_err: %ld | %x %x %x\n",
							 i,
							 cwd[i],
							 est_cwd[i],
							 tst_vct[tv_idx][i]);
			}
			for(i = 0; i < MESSAGE_LEN; i++)
			{
				fprintf(frc, "msg_poly[%ld] = 0x%x;\n",
							 i,
							 msg_poly[i]);
			}
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				fprintf(frc, "cwd_poly[%ld] = 0x%x;\n",
							 i,
							 cwd_poly[i]);
			}
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				fprintf(frc, "tst_vct[tv_idx][%ld] = 0x%x;\n",
							 i,
							 tst_vct[tv_idx][i]);
			}
			for(i = 0; i < symbol_num; i++)
			{
				fprintf(frc, "recv_seq[%ld][0] = %f;\n",
							 i,
							 recv_seq[i][0]);
				fprintf(frc, "recv_seq[%ld][1] = %f;\n",
							 i,
							 recv_seq[i][1]);
			}
			DEBUG_SYS("Transmission over Channel:\n");
			for(i = 0; i < symbol_num; i++)
			{
				fprintf(frc, "recv_seq[%ld][0] = %f;\n", i, recv_seq[i][0]);
				fprintf(frc, "recv_seq[%ld][1] = %f;\n", i, recv_seq[i][1]);
			}
			for(i = 0; i < MESSAGE_LEN; i++)
			{
				fprintf(frc, "msg_poly[%ld] = 0x%x;\n",
			            i,
			            msg_poly[i]);
			}
			DEBUG_SYS("\n");
			fclose(frc);
			frc = NULL;
#endif
		}		
	}

	return val;
}

int check_result_msg(unsigned char *msg, unsigned char *est_msg, long long tv_idx)
{
	long long i = 0;
	long long genus = 0, radius = 0, cwd_err = 0;
	int val = 0;

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		if(msg[i] != est_msg[i])
		{
			DEBUG_IMPORTANT("check_result_msg: %ld | %x %x\n",
			               i,
			               msg[i],
			               est_msg[i]);
			val++;
		}
	}
	if(0 == val)
	{
		DEBUG_IMPORTANT("msg_ok\n");
		if(GF_Q > min_intp_idx)
		{
			DEBUG_IMPORTANT("min_intp_idx: %ld\n", min_intp_idx);
		}
	}
	else
	{
		//err_cnt++;
		genus = GF_Q * (GF_Q - 1) / 2;
		//radius = (CODEWORD_LEN - MESSAGE_LEN - genus + 1 - 1) / 2;
		/*for GS (m = 1), that is the GS decoding bound*/
		radius = (CODEWORD_LEN - (MESSAGE_LEN + genus - 1) - 2 * genus) / 2;
		if(0 == ((MESSAGE_LEN - genus + 1) % 2))
		{
			radius = radius - 1;
		}
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			//if(cwd_poly[i] != recv_poly[i])
			if(cwd_poly[i] != tst_vct[tv_idx][i])
			{
				cwd_err++;
			}
		}
		if(cwd_err > radius)
		{
			DEBUG_IMPORTANT("err_out_of_radius_msg: %ld %ld %ld\n",
			               cwd_err,
			          	   radius,
			          	   val);
		}
		else
		{
			DEBUG_SYS("msg decoding err: %ld %ld\n",
			          cwd_err,
			          radius);
			for(i = 0; i < MESSAGE_LEN; i++)
			{
				DEBUG_SYS("check_result_msg_err: %ld | %x %x\n",
				          i,
				          msg[i],
				          est_msg[i]);
			}
#if 0			
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				DEBUG_SYS("check_result_cwd_err: %ld | %x %x\n",
			              i,
			              cwd_poly[i],
			              est_cwd_poly[i]);
			}
#endif			
			for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
			{
				if(0xFF != min_intp_poly[i])
				{
					DEBUG_SYS("min_intp_poly_err: %ld %ld %ld | %x\n",
					          x_term_degree_table[i],
					          y_term_degree_table[i],
					          z_term_degree_table[i],
					          min_intp_poly[i]);
				}
			}
		}
	}

	return val;
}

int koetter_intp_her_lcc(unsigned char *test_poly_seq,
							 long long tv_idx,
                             long long intp_start_idx,
                             long long intp_end_idx,
                             unsigned char input_poly[][MAX_POLY_TERM_SIZE],
                             unsigned char *output_poly)
{
	long long i = 0, j = 0, k = 0;
	unsigned char hasse_dev[KOT_INTP_POLY_NUM], q0_val = 0xFF, q1_val = 0xFF;
	memset(hasse_dev, 0xFF, sizeof(unsigned char) * KOT_INTP_POLY_NUM);
	long long min_idx = -1, min_degree = 65536, degree_val = 0, tmp_val = 0;
	unsigned char z_flag = 0;
	long long intp_idx = 0;
	long long l_idx = 0, p_idx = 0;
#if 1
	/*for eta == 0 case*/
	if(intp_start_idx > intp_end_idx)
	{
		min_idx = -1;
		min_degree = 65536;
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			intp_poly_degree[j] = poly_degree_cal(input_poly[j]);
			
#if (1 == CFG_FIND_LATER_POLY)			
			if(min_degree >= intp_poly_degree[j])
#else
			if(min_degree > intp_poly_degree[j])
#endif		
			{
				min_degree = intp_poly_degree[j];
				min_idx = j;
			}
		}
		//DEBUG_SYS("min_idx: %ld\n", min_idx);
	}
#endif
	/*intp_poly_tmp is used to store intp. med. result here*/
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		intp_poly_degree[i] = poly_degree_cal(input_poly[i]);
		memcpy(intp_poly_tmp[i], input_poly[i], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

#if (1 == CFG_INTP_ORDER_ERR)//for debug
	long long tmp_idx = 0;
	for(i = 0; i < 26; i++)
	{
		tmp_idx = chnl_rel_order_idx[CODEWORD_LEN - 1 - i];
		chnl_rel_order_idx[CODEWORD_LEN - 1 - i] = i;
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if(i == chnl_rel_order_idx[j])
			{
				break;
			}
		}
		chnl_rel_order_idx[j] = tmp_idx;
	}
	tmp_idx = chnl_rel_order_idx[2];
	for(i = 2; i < 37; i++)
	{
		chnl_rel_order_idx[i] = chnl_rel_order_idx[i + 1];
	}
	chnl_rel_order_idx[37] = tmp_idx;
	long long tmp_iter = 28;
	for(i = 0; i < tmp_iter; i++)
	{
		tmp_idx = chnl_rel_order_idx[37 - i];
		chnl_rel_order_idx[37 - i] = chnl_rel_order_idx[36 - i];
		chnl_rel_order_idx[36 - i] = tmp_idx;
	}
#endif	

	//for(i = intp_start_idx; i < CODEWORD_LEN; i++)
	for(i = intp_start_idx; i <= intp_end_idx; i++)
	{
		intp_idx = chnl_rel_order_idx[CODEWORD_LEN - 1 - i];
		//intp_idx = i;//there is a bug

		DEBUG_NOTICE("intp point: %ld %ld | %x %x | %x %x\n",
		             i,
					 intp_idx,
					 af_pnt[intp_idx][0],
					 af_pnt[intp_idx][1],
					 test_poly_seq[intp_idx],
					 cwd_poly[intp_idx]);
		/**notice that these two val should be clear*/
		min_idx = -1;
		min_degree = 65536;

#if (1 == CFG_RET)
		if(1 == keep_flag[intp_idx])
		{
			DEBUG_NOTICE("ReT Keep Skip: %ld %ld\n", i, intp_idx);
			continue;
		}
#endif

		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			/*compute hasse dev.*/
#if 0
			poly_q0q1_get(j);
			q0_val = poly_eva_x_y(q0_poly_coef, af_pnt[i][0], af_pnt[i][1]);
			q1_val = poly_eva_x_y(q1_poly_coef, af_pnt[i][0], af_pnt[i][1]);
			q1_val = gf_multp(q1_val, test_poly_seq[i]);
			hasse_dev[j] = gf_add(q0_val, q1_val);
			DEBUG_NOTICE("hasse_dev: %ld | %ld | %x %x | %x | %x\n",
						 i,
						 j,
						 q0_val,
						 q1_val,
						 test_poly_seq[i],
						 hasse_dev[j]);
#else
			hasse_dev[j] = hasse_dev_cal(intp_idx, j, test_poly_seq[intp_idx]);
#endif

			/*find the min poly whose hasse dev is not zero*/
			if(0xFF != hasse_dev[j])
			{
				degree_val = poly_degree_cal(input_poly[j]);
				if(degree_val < min_degree)
				{
					min_degree = degree_val;
					min_idx = j;
				}
				if((degree_val == min_degree)
					&& (-1 != min_idx))
				{
					tmp_val = poly_degree_compare(input_poly[j], input_poly[min_idx]);
					if(2 == tmp_val)
					{
						min_degree = degree_val;
						min_idx = j;
					}
					DEBUG_NOTICE("equal_poly_degree: %d\n", tmp_val);
				}
			}
		}
		DEBUG_NOTICE("min_poly: %ld | %ld\n",
					 min_idx,
					 min_degree);

		/*update poly*/
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			if(0xFF != hasse_dev[j])
			{
				if(min_idx == j)
				{
					poly_min_update(intp_poly_tmp[j], input_poly[j], af_pnt[intp_idx][0], hasse_dev[j]);
				}
				else
				{
					poly_normal_update(intp_poly_tmp[j], input_poly[j], input_poly[min_idx], hasse_dev[min_idx], hasse_dev[j]);
				}
			}

			/*convert poly to avoid high x-degree*/
			her_convert(intp_poly_tmp[j]);
		}

		min_idx = -1;
		min_degree = 65536;
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			z_flag = 0;
			/*notice that the prev poly should be stored for updating properly*/
			memcpy(input_poly[j], intp_poly_tmp[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			intp_poly_degree[j] = poly_degree_cal(input_poly[j]);
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != input_poly[j][k])
				{
					DEBUG_NOTICE("intp_poly_coef_update: %ld | %ld %ld %ld | %ld | %x\n",
								 j,
								 x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 intp_poly_degree[j],
								 input_poly[j][k]);
					if(0 != z_term_degree_table[k])
					{
						z_flag = 1;
					}
				}
			}
#if 0//to avoid non-z-term poly.
			if((min_degree > intp_poly_degree[j])
				&& (0 != z_flag))
#else
#if 0
			if(((min_degree > intp_poly_degree[j])
				&& ((CODEWORD_LEN - 1) != i))
					||
				((min_degree > intp_poly_degree[j])
				&& ((CODEWORD_LEN - 1) == i)
				&& (0 != z_flag)))
#else

#if (1 == CFG_FIND_LATER_POLY)
			if(min_degree >= intp_poly_degree[j])
#else
			if(min_degree > intp_poly_degree[j])
#endif
			
#endif
#endif
			{
#if 0
				if((GF_Q > j)
					&& ((CODEWORD_LEN - 1) == i))
				{
					continue;
				}
#endif				
				min_degree = intp_poly_degree[j];
				min_idx = j;
			}
		}
		
		if(i > (CODEWORD_LEN - 1 - ETA))
		{
			l_idx = i - (CODEWORD_LEN - 1 - ETA);
			p_idx = tv_idx / ((long long)(pow(2, ETA - l_idx)));
			kot_node_save(input_poly, l_idx, p_idx);
		}
		intp_cnt++;
		DEBUG_NOTICE("intp: %ld | %ld\n", tv_idx, intp_cnt);
	}

#if 1
	her_convert(input_poly[min_idx]);
#endif

	//min_idx = 5;//for debug

	memcpy(output_poly, input_poly[min_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	/*notice that this function should be kept*/
	//poly_q0q1_get_new(output_poly);
	//min_intp_idx = min_idx;

#if (0 == CFG_FAC_FREE)
#if (1 == CFG_RET)
	if(0 != intp_start_idx)/*this step should not be excuted in cmm intp*/
	{
		unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
		poly_q0q1_get_new(output_poly);
		memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		poly_mul(v_poly, q0_poly_coef, tmp_poly);
		memcpy(q0_poly_coef, tmp_poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		memset(output_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		memcpy(output_poly, q0_poly_coef, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0xFF != q1_poly_coef[i])
			{
#if (0 == CFG_QUICK_POLY_SEARCH)
				for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
				{
					if((x_term_degree_table[j] == x_term_degree_table[i])
						&& (y_term_degree_table[j] == y_term_degree_table[i])
						&& (1 == z_term_degree_table[j]))
					{
						output_poly[j] = gf_add(output_poly[j], q1_poly_coef[i]);
						break;
					}
				}
#else
				j = term_search(x_term_degree_table[i], y_term_degree_table[i], 1);
				output_poly[j] = gf_add(output_poly[j], q1_poly_coef[i]);
#endif
			}
		}
	}
#endif
#endif

#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != output_poly[i])
		{
			DEBUG_NOTICE("min_intp_poly: %ld | %ld %ld %ld | %x\n",
						 min_idx,
						 x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 output_poly[i]);
		}
	}
#endif

	return 0;
}

int fac_msg2cwd(unsigned char *msg,
              		 unsigned char *cwd)
{
	her_encoding(msg, cwd);
#if (1 == CFG_RET)
	long long i = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		cwd[i] = gf_add(cwd[i], ret_cwd_poly[i]);
	}
	//cwd2msg(cwd, msg);
#endif

	return 0;
}

int her_era_cwd_gen(unsigned char *output_cwd)
{
	long long i = 0, j = 0, k = 0, l = 0, m = 0;
	long long k_idx = 0;
	unsigned char l_val = 0x0, r_val = 0xFF;
	unsigned char y_her_curve_find_flag = 0;
	unsigned char x_i_sel_flag = 0, y_i_sel_flag = 0;
	unsigned char x_j_sel_flag = 0, y_j_sel_flag = 0;
	unsigned char tmp_val = 0xFF;

	memcpy(keep_sym, output_cwd, sizeof(unsigned char) * CODEWORD_LEN);

	for(k = 0; k < CODEWORD_LEN; k++)
	{
		if((0 == erasure_flag[k])
			|| (0 == dev_to_flag[k]))
		{
			output_cwd[k] = output_cwd[k];//useless
			continue;
		}
		for(k_idx = 0; k_idx < GF_FIELD; k_idx++)
		{
			if(power_polynomial_table[k_idx][0] == af_pnt[k][0])
			{
				break;
			}
		}
		r_val = 0xFF;
		for(m = 0; m < CODEWORD_LEN; m++)
		{
			if(1 == erasure_flag[m])
			{
				continue;
			}
			l_val = 0x0;
			x_j_sel_flag = 0;
			y_j_sel_flag = 0;
			
			j = m / GF_Q;

			if(af_pnt[m][0] == power_polynomial_table[j][0])
			{
				x_j_sel_flag = 1;
				for(i = 0; i < GF_FIELD; i++)
				{
					x_i_sel_flag = 0;						
					l = i * GF_Q;

					if((0 == erasure_flag[l])
						&& (af_pnt[m][0] != power_polynomial_table[i][0]))
					{
						x_i_sel_flag = 1;								
					}

					if(1 == x_i_sel_flag)
					{
						l_val = gf_multp(l_val, t_val[i][j][k_idx]);
						DEBUG_NOTICE("x_r_val_cal: %ld %ld | %x %x | %x %x %x %x\n",
						             k,
						             m,
						             power_polynomial_table[j][0],
						             power_polynomial_table[i][0],
						             t_val[i][j][k_idx],
						             l_val,
						             keep_sym[j],
						             r_val);
					}
				} 
			}

			if((0xFF == y_t_val[m][k])
				|| (0xFF == l_val))
			{
				l_val = 0xFF;
			}
			else if(0x0 == y_t_val[m][k])
			{
				l_val = l_val;
			}
			else if(0x0 == l_val)
			{
				l_val = y_t_val[m][k];
			}
			else
			{
				l_val = gf_multp(l_val, y_t_val[m][k]);
			}

			if((0xFF == keep_sym[m])
				|| (0xFF == l_val))
			{
				tmp_val = 0xFF;
			}
			else if(0x0 == keep_sym[m])
			{
				tmp_val = l_val;
			}
			else if(0x0 == l_val)
			{
				tmp_val = keep_sym[m];
			}
			else
			{
				tmp_val = gf_multp(keep_sym[m], l_val);
			}
			if((0xFF == r_val)
				&& (0xFF == tmp_val))
			{
				r_val = 0xFF;
			}
			else if((0xFF != r_val)
					&& (0xFF == tmp_val))
			{
				r_val = r_val;
			}
			else if((0xFF == r_val)
					&& (0xFF != tmp_val))
			{
				r_val = tmp_val;
			}
			else
			{
				r_val = gf_add(r_val, tmp_val);
			}
		}
		DEBUG_NOTICE("ret_cwd_gen: %ld | %x %x\n",
		             k,
		             r_val,
		             output_cwd[k]);

		output_cwd[k] = r_val;
	}

	return 0;
}

int her_ret_fac_free(unsigned char *q0_poly,
						unsigned char *q1_poly,
						unsigned char *v_poly,
						unsigned char *output_cwd)
{
	long long i = 0, j = 0;
	unsigned char q0_val = 0xFF, q1_val = 0xFF;
	unsigned char dev_init_flag = 0;
	long long erasure_keep_cnt = 0;
	memset(erasure_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
	memset(dev_to_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
	dev_to_cnt = 0;
	sing_era_cnt = 0;
	long long fac_free_sing_limit = 0;
	keep_cnt = 0;

	unsigned char q0_all_zero = 1;
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q0_poly_coef[i])
		{
			q0_all_zero = 0;
			break;
		}
	}
	if(1 == q0_all_zero)
	{
		memset(output_cwd, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
		return 2;
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		q1_val = poly_eva_x_y(q1_poly, af_pnt[i][0], af_pnt[i][1]);
		if(0xFF != q1_val)
		{
			if(0xFF == v_val[i])
			{
				output_cwd[i] = 0xFF;
			}
			else
			{
				q0_val = poly_eva_x_y(q0_poly, af_pnt[i][0], af_pnt[i][1]);
				if(0xFF == q0_val)
				{
					output_cwd[i] = 0xFF;
				}
				else
				{
					output_cwd[i] = gf_multp(v_val[i], q0_val);
					output_cwd[i] = gf_div(output_cwd[i], q1_val);
				}
			}
		}
		else
		{
			if(0xFF == v_val[i])
			{
				output_cwd[i] = ret_fac_free_dev(i, dev_init_flag, 1);
				q0_val = poly_eva_x_y(q0_poly, af_pnt[i][0], af_pnt[i][1]);
				output_cwd[i] = gf_multp(output_cwd[i], q0_val);
				dev_init_flag = 1;
			}
			else
			{
				q0_val = poly_eva_x_y(q0_poly, af_pnt[i][0], af_pnt[i][1]);
				if(0xFF != q0_val)
				{
					output_cwd[i] = 0xFF;

					return 1;
				}
				else
				{
					output_cwd[i] = ret_fac_free_dev(i, dev_init_flag, 0);
					output_cwd[i] = gf_multp(output_cwd[i], v_val[i]);
					dev_init_flag = 1;
				}
			}

			if(1 == erasure_flag[i])
			{
				dev_to_cnt++;
			}
		}
#if (1 == TEST_MODE)
		q0_val = poly_eva_x_y(q0_poly, af_pnt[i][0], af_pnt[i][1]);
#endif
		DEBUG_NOTICE("ret_fac_free_cal: %ld | %d | %x %x %x | %x\n",
				     i,
				     erasure_flag[i],
				     q0_val,
				     q1_val,
				     v_val[i],
				     output_cwd[i]);
#if (1 == TEST_MODE)
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if(1 == erasure_flag[j])
			{
				DEBUG_IMPORTANT("erasure_flag_check: %d %ld\n", j, dev_to_cnt);
			}
		}
#endif
	}

	if(0 != dev_to_cnt)
	{
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			if(1 == erasure_flag[i])
			{
				DEBUG_IMPORTANT("erasure_flag_find: %ld\n", i);
				for(j = 0; j < CODEWORD_LEN; j++)
				{
					if(af_pnt[i][0] == af_pnt[j][0])
					{
						erasure_flag[j] = 1;
						DEBUG_IMPORTANT("erasure_flag_set: %ld\n", j);
					}
				}
			}
		}
		
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			if(1 == erasure_flag[i])
			{
				sing_era_cnt++;
			}
		}

		DEBUG_IMPORTANT("dev_to_cnt: %ld\n", dev_to_cnt, sing_era_cnt);

		fac_free_sing_limit = CODEWORD_LEN
		                    - (GF_Q
		                       * (long long)(ceil((MESSAGE_LEN + (GF_Q) * (GF_Q - 1) / 2.0) / GF_Q)));
		DEBUG_IMPORTANT("fac_free_sing_limit: %ld\n", fac_free_sing_limit);
		if(fac_free_sing_limit >= sing_era_cnt)
		{
			her_era_cwd_gen(output_cwd);
		}
		else
		{
			if(CODEWORD_LEN != sing_era_cnt)
			{
				DEBUG_SYS("Call Recur. Fac.: %ld\n", sing_era_cnt, fac_free_sing_limit);
#if (1 == DEV_RECORD)
				dev_out_cnt++;
#endif				
			}
#if (1 == DEV_RECORD)			
			else
			{
				dev_2_avg_cnt = dev_2_avg_cnt - CODEWORD_LEN;
			}
#endif			
		}
	}
	
#if (1 == CFG_PRG_DECODING)	
#if (1 == CFG_CWD_DIM_CHECK)
	long long cwd_dim = check_cwd_dimension(output_cwd);
	DEBUG_NOTICE("output_cwd_dim: %ld\n", cwd_dim);
	if(MESSAGE_LEN < cwd_dim)
#endif
	{
#if 0	
		DEBUG_NOTICE("force erasure decoding\n");
		memset(erasure_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
		memset(dev_to_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
#if 0		
		memset(erasure_flag, 1, sizeof(unsigned char) * GF_Q);
		memset(dev_to_flag, 1, sizeof(unsigned char) * GF_Q);
#else
		memset(erasure_flag, 1, sizeof(unsigned char) * keep_cnt);
		memset(dev_to_flag, 1, sizeof(unsigned char) * keep_cnt);
#endif
		her_era_cwd_gen(output_cwd);
		cwd_dim = check_cwd_dimension(output_cwd);
		DEBUG_NOTICE("cwd_dim_after: %ld\n", cwd_dim);
		if(MESSAGE_LEN < cwd_dim)
		{
			memset(output_cwd, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
		}
#endif
		/*clear this illegal est_cwd*/
		memset(output_cwd, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	}
#endif	

	return 0;
}

int her_fac_free(unsigned char *poly, unsigned char *est_cwd)
{
#if (1 == CFG_RET)
	poly_q0q1_get_new(poly);
	her_ret_fac_free(q0_poly_coef,
					 q1_poly_coef,
					 v_poly,
					 est_cwd);
#else

#endif

	return 0;
}

int her_fac(unsigned char *poly,
              unsigned char *est_msg,
              unsigned char *est_cwd)
{
	//start = clock();

#if (1 == CFG_FAC_FREE)
	poly_q0q1_get_new(poly);
	her_ret_fac_free(q0_poly_coef, q1_poly_coef, v_poly, est_cwd);

	if(((CODEWORD_LEN
		 - (GF_Q * (long long)(ceil((MESSAGE_LEN + (GF_Q) * (GF_Q - 1) / 2.0) / GF_Q))))
		 < sing_era_cnt)
		 && (CODEWORD_LEN != sing_era_cnt))
	{
		DEBUG_SYS("Recall Recur. Fac. to Gen. Msg: %ld\n", sing_era_cnt);
#if (1 == CFG_RET)/*notice that this q0 should be mul by v in this case*/
		long long i = 0, j = 0;
		unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
		memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		poly_mul(v_poly, q0_poly_coef, tmp_poly);
		memcpy(q0_poly_coef, tmp_poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		memset(poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		memcpy(poly, q0_poly_coef, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0xFF != q1_poly_coef[i])
			{
				j = term_search(x_term_degree_table[i], y_term_degree_table[i], 1);
				poly[j] = gf_add(poly[j], q1_poly_coef[i]);
			}
		}
#endif
		factorization_recur(poly, est_msg);
		fac_msg2cwd(est_msg, est_cwd);

		return 0;
	}

#if (1 == CFG_RET)	
	long long i = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		est_cwd[i] = gf_add(est_cwd[i], ret_cwd_poly[i]);
		DEBUG_NOTICE("est_cwd: %ld | %x\n",
		             i,
		             est_cwd[i]);
	}
#endif

	//cwd2msg(est_cwd, est_msg);
#else

	factorization_recur(poly, est_msg);
	fac_msg2cwd(est_msg, est_cwd);

#endif

	//stop = clock();
	//runtime = runtime + (stop - start) / 1000.0000;

	return 0;
}

int kot_node_init()
{
	long long i = 0, j = 0;
	node_num = 0;
	for(i = 0; i <= ETA; i++)
	{
		node_num = node_num + (long long)(pow(2, i));
	}
	DEBUG_NOTICE("node_num: %ld\n", node_num);

	node_poly = (unsigned char***)malloc(sizeof(unsigned char**) * node_num);
	for(i = 0; i < node_num; i++)
	{
		node_poly[i] = (unsigned char**)malloc(sizeof(unsigned char*) * KOT_INTP_POLY_NUM);
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			node_poly[i][j] = (unsigned char*)malloc(sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
	}

	for(i = 0; i < node_num; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			memset(node_poly[i][j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
	}

	node_tore_flag = (unsigned char*)malloc(sizeof(unsigned char) * node_num);
	memset(node_tore_flag, 0, sizeof(unsigned char) * node_num);

	return 0;
}

int kot_node_exit()
{
	long long i = 0, j = 0;

	for(i = 0; i < node_num; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			free(node_poly[i][j]);
			node_poly[i][j] = NULL;
		}
		free(node_poly[i]);
		node_poly[i] = NULL;
	}
	free(node_poly);
	node_poly = NULL;

	free(node_tore_flag);
	node_tore_flag = NULL;

	return 0;
}

int kot_node_clear()
{
	long long i = 0, j = 0;

	for(i = 0; i < node_num; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			memset(node_poly[i][j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
	}

	memset(node_tore_flag, 0, sizeof(unsigned char) * node_num);

	return 0;
}

long long kot_node_layer_get(long long node_idx)
{
	long long i = 0;
	long long layer_idx = 0;
	long long tmp_val_left = 0, tmp_val_right = 0;

	for(i = 0; i < ETA; i++)
	{
		tmp_val_left = (long long)(pow(2, i)) - 1;
		tmp_val_right = tmp_val_left + (long long)(pow(2, i)) - 1;

		if((tmp_val_left <= node_idx)
			&& (tmp_val_right >= node_idx))
		{
			layer_idx = i;
			break;
		}
	}

	return layer_idx;
}

long long kot_node_place_get(long long node_idx)
{
	long long i = 0;
	long long place_idx = 0;
	long long tmp_val_left = 0, tmp_val_right = 0;

	for(i = 0; i < ETA; i++)
	{
		tmp_val_left = (long long)(pow(2, i)) - 1;
		tmp_val_right = tmp_val_left + (long long)(pow(2, i)) - 1;

		if((tmp_val_left <= node_idx)
			&& (tmp_val_right >= node_idx))
		{
			place_idx = node_idx - tmp_val_left;
			break;
		}
	}

	return place_idx;
}

int kot_node_save(unsigned char poly[][MAX_POLY_TERM_SIZE], long long layer_idx, long long place_idx)
{
	long long i = 0;
	long long node_idx = 0;

	node_idx = (long long)(pow(2, layer_idx)) - 1 + place_idx;
	DEBUG_NOTICE("kot_node_save: %ld %ld | %ld | %d\n",
	             layer_idx,
	             place_idx,
	             node_idx,
	             node_tore_flag[node_idx]);
	
	if(0 == node_tore_flag[node_idx])
	{
		for(i = 0; i < KOT_INTP_POLY_NUM; i++)
		{
			memcpy(node_poly[node_idx][i], poly[i], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
		node_tore_flag[node_idx] = 1;
		DEBUG_NOTICE("kot_node_save_OK: %ld %ld | %ld | %d\n",
	             	layer_idx,
	             	place_idx,
	             	node_idx,
	             	node_tore_flag[node_idx]);
	}

	return 0;
}

long long kot_node_load(unsigned char poly[][MAX_POLY_TERM_SIZE], long long tv_idx)
{
	long long i = 0;
	long long l_idx = 0, p_idx = 0, node_idx = 0;
	long long div_num = 0;

#if 1
	if(0 == ETA)// notice this, eta can be zero
	{
		return l_idx;
	}
#endif

	for(l_idx = ETA - 1; l_idx > 0; l_idx--)
	{
		div_num = (long long)(pow(2, (ETA - l_idx)));

		p_idx = tv_idx / div_num;
		node_idx = (long long)(pow(2, l_idx)) - 1 + p_idx;
		if(0 != node_tore_flag[node_idx])
		{
			for(i = 0; i < KOT_INTP_POLY_NUM; i++)
			{
				memcpy(poly[i], node_poly[node_idx][i], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
			DEBUG_NOTICE("kot_node_load_OK: %ld %ld | %ld | %ld\n",
			             l_idx,
			             p_idx,
			             node_idx,
			             tv_idx);
			break;
		}
	}

	return l_idx;
}

int her_lcc_tv_init(long long tv_idx)
{
	poly_init();

	return 0;
}

int her_lcc_tv_process(long long tv_idx)
{
	long long layer_idx = 0;

	her_lcc_tv_init(tv_idx);

	//layer_idx = kot_node_load(intp_poly_coef, tv_idx);

#if (1 == CFG_RET)
	koetter_intp_her_lcc(ret_tst_vct[tv_idx],
						 tv_idx,
						 layer_idx,
						 CODEWORD_LEN - 1,
						 intp_poly_coef,
						 min_intp_poly);
#else
	koetter_intp_her_lcc(tst_vct[tv_idx],
						 tv_idx,
						 layer_idx,
						 CODEWORD_LEN - 1,
						 intp_poly_coef,
						 min_intp_poly);
#endif

	her_fac(min_intp_poly, tv_est_msg[tv_idx], tv_est_cwd[tv_idx]);

#if (1 == TEST_MODE)
	long long i = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("est_cwd: %ld | %ld | %x %x\n",
					 tv_idx,
					 i,
					 tv_est_cwd[tv_idx][i],
					 cwd_poly[i]);
	}
#endif

	return 0;
}

int check_result_tv(unsigned char *tv, unsigned char *est_cwd, unsigned char *est_msg)
{
	long long i = 0;
	long long genus = 0, radius = 0, cwd_err = 0, msg_err = 0, tv_diff = 0;
	int val = 0;

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		if(msg_poly[i] != est_msg[i])
		{
			DEBUG_IMPORTANT("check_result_msg: %ld | %x %x\n",
			               i,
			               msg_poly[i],
			               est_msg[i]);
			msg_err++;
		}
	}
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(cwd_poly[i] != est_cwd[i])
		{
			DEBUG_IMPORTANT("check_result_cwd: %ld | %x %x\n",
			               i,
			               cwd_poly[i],
			               est_cwd[i]);
			cwd_err++;
		}
	}

	if((0 == msg_err)
		&& (0 == cwd_err))
	{
		DEBUG_IMPORTANT("dec_ok\n");
		val = 0;
	}
	else
	{
		genus = GF_Q * (GF_Q - 1) / 2;
		//radius = (CODEWORD_LEN - MESSAGE_LEN - genus + 1 - 1) / 2;
		radius = (CODEWORD_LEN - (MESSAGE_LEN + genus - 1) - 2 * genus) / 2;
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			if(cwd_poly[i] != tv[i])
			{
				tv_diff++;
			}
		}
		if(tv_diff > radius)
		{
			val = 1;
			DEBUG_IMPORTANT("err_out_of_radius: %ld %ld\n",
			                tv_diff,
			          	    radius);
		}
		else
		{
			val = 2;
			DEBUG_SYS("msg decoding err: %ld %ld\n",
			          tv_diff,
			          radius);
			for(i = 0; i < MESSAGE_LEN; i++)
			{
				DEBUG_SYS("check_result_msg_err: %ld | %x %x\n",
				          i,
				          msg_poly[i],
				          est_msg[i]);
			}
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				DEBUG_SYS("check_result_cwd_err: %ld | %x %x\n",
			              i,
			              cwd_poly[i],
			              est_cwd[i]);
			}
		}
	}

	return val;
}

int her_lcc_check_result()
{
	long long i = 0, j = 0;
	unsigned char err_this_frame = 1;
	int cwd_err_val = CODEWORD_LEN, msg_err_val = MESSAGE_LEN;
	int check_val = 0;
	long long genus = GF_Q * (GF_Q - 1) / 2;
	//long long radius = (CODEWORD_LEN - MESSAGE_LEN - genus + 1 - 1) / 2;
	/*for GS (m = 1), that is the GS decoding bound*/
	long long radius = (CODEWORD_LEN - (MESSAGE_LEN + genus - 1) - 2 * genus) / 2;
	if(0 == ((MESSAGE_LEN - genus + 1) % 2))
	{
		radius = radius - 1;
	}
	unsigned char err_dump = 0;

#if (1 == CFG_SYS_GEN)
	check_val = 0;
	check_val = sys_check_ret_cwd();
	if(0 == check_val)
	{
		sys_ret_ok_cnt++;
		err_this_frame = 0;
		DEBUG_NOTICE("sys_ret_ok_cnt: %ld\n", sys_ret_ok_cnt);

		//return 0;
	}
	
	check_val = 0;
	check_val = sys_ret_try_cwd_check();
	if(0 == check_val)
	{
		err_this_frame = 0;
	}
#endif
#if (1 == CFG_RET)/*add this to check the ret_cwd_poly*/
	check_val = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(ret_cwd_poly[i] != cwd_poly[i])
		{
			//DEBUG_SYS("lag_ret_diff: %ld | %x %x\n", i, ret_cwd_poly[i], cwd_poly[i]);
			check_val = 1;
			break;
		}
	}
	if(0 == check_val)
	{
		lag_ret_ok_cnt++;
		memcpy(est_cwd_poly, ret_cwd_poly, sizeof(unsigned char) * CODEWORD_LEN);
		cwd2msg(ret_cwd_poly, est_msg_poly);
		err_this_frame = 0;
		DEBUG_NOTICE("lag_ret_ok_cnt: %ld\n", lag_ret_ok_cnt);

		//return 0;
	}
#endif
	
	for(i = 0; i < tst_vct_num; i++)
	{
		if(0 == tv_dec_output_flag[i])
		{
			/*this tv. has no output, no need to check*/
			DEBUG_NOTICE("tv_skip_check: %ld\n", i);
			continue;
		}

		cwd_err_val = CODEWORD_LEN;
		msg_err_val = MESSAGE_LEN;

		cwd_err_val = check_result_cwd(cwd_poly, tv_est_cwd[i], i);
#if (1 == CFG_FAC_FREE)
		if(0 == cwd_err_val)
		{
			cwd2msg(tv_est_cwd[i], tv_est_msg[i]);
		}
		else
		{
			continue;
		}
#else

#if (1 == CFG_RET)
		if(0 == cwd_err_val)
		{
			cwd2msg(tv_est_cwd[i], tv_est_msg[i]);
		}
		else
		{
			continue;
		}
#endif

#endif

		msg_err_val = check_result_msg(msg_poly, tv_est_msg[i], i);
		if((0 == cwd_err_val)
			&& (0 == msg_err_val))
		{
			err_this_frame = 0;
			DEBUG_NOTICE("dec_suc_tv: %ld\n", i);
			break;
		}
	}

	if(0 != err_this_frame)
	{
		for(i = 0; i < tst_vct_num; i++)
		{
			if(tv_err[i] <= radius)
			{
				DEBUG_SYS("dec_err: %ld | %ld\n", i, tv_err[i]);
#if (1 == OUTPUT_LOG)
				frc = fopen(log_name, "a+");
				fprintf(frc, "dec_err: %ld | %ld\n", i, tv_err[i]);

				if(0 == err_dump)
				{
					for(j = 0; j < symbol_num; j++)
					{
						fprintf(frc, "recv_seq[%ld][0] = %f;\n", j, recv_seq[j][0]);
						fprintf(frc, "recv_seq[%ld][1] = %f;\n", j, recv_seq[j][1]);
					}
					for(j = 0; j < MESSAGE_LEN; j++)
					{
						fprintf(frc, "msg_poly[%ld] = 0x%x;\n",
					            j,
					            msg_poly[j]);
					}
					err_dump = 1;
				}

				fclose(frc);
				frc = NULL;
#endif
			}
		}
		err_cnt++;
	}

#if (1 == TEST_MODE)
	for(i = 0; i < tst_vct_num; i++)
	{
		if(0 == tv_dec_output_flag[i])
		{
			/*this tv. has no output, no need to check*/
			continue;
		}
		check_val = check_result_tv(tst_vct[i], tv_est_cwd[i], tv_est_msg[i]);
		DEBUG_NOTICE("dec_tv_check: %ld | %ld | %ld\n", i, tv_err[i], check_val);
	}
#endif

	return 0;
}

int her_cwd_check(unsigned char *cwd)
{
	long long i = 0, j = 0, k = 0;
	unsigned char x_val = 0xFF, y_val = 0xFF, val = 0xFF;

	for(i = 0; i < GF_FIELD; i++)
	{
		for(j = 0; j < GF_FIELD; j++)
		{
			val = 0xFF;
			for(k = 0; k < CODEWORD_LEN; k++)
			{
				x_val = gf_pow_cal(af_pnt[i][0], i);
				y_val = gf_pow_cal(af_pnt[i][1], j);
				val = gf_add(val,
							 gf_multp(gf_add(cwd[k], cwd[k]),
									  gf_multp(x_val, y_val)));
			}
			DEBUG_NOTICE("ret_cwd_check: %ld %ld | %x\n",
						 i,
						 j,
						 val);
		}
	}

	return 0;
}

/*this should be modified for AG codes*/
int MLcriterion(unsigned char est_cwd[], unsigned char recv_cwd[])
{
    long long i, j;
    float gama1, gama2;
    unsigned char same_pos_vector[CODEWORD_LEN], differ_pos_vector[CODEWORD_LEN];
    long long count1 = 0, count2 = 0;
    long long Hammin_dis, differ_pos, same_pos;
    float symbol_reliabilty_metric[CODEWORD_LEN];
    float temp, temp_metric[CODEWORD_LEN];
    float max_symbol_reliability, sec_max_symbol_reliability;
    long long max_row, sec_max_row;
    unsigned char cwd[CODEWORD_LEN];
    long long code_min_dis = CODEWORD_LEN - MESSAGE_LEN + 1 - (GF_Q) * (GF_Q - 1) / 2;

    for (i = 0; i < CODEWORD_LEN; i++)
    {
        same_pos_vector[i] = 0xFF;
        differ_pos_vector[i] = 0xFF;
        temp_metric[i] = 0;
        cwd[i] = est_cwd[i];
    }

    for (i = 0; i < CODEWORD_LEN; i++)
    {
        if (cwd[i] == recv_cwd[i])
        {
            same_pos_vector[count1] = i;
            count1++;
        }
        else
        {
            differ_pos_vector[count2] = i;
            count2++;
        }
    }
    Hammin_dis = count2;
    gama1 = 1;
    for (j = 0; j < count2; j++)
    {
        differ_pos = differ_pos_vector[j];
        max_symbol_reliability = 0;
        for (i = 0; i < GF_FIELD; i++)
        {
            if (chnl_rel_matrix[i][differ_pos] > max_symbol_reliability)
                max_symbol_reliability = chnl_rel_matrix[i][differ_pos];
        }
        for (i = 0; i < GF_FIELD; i++)
        {
            if (power_polynomial_table[i][0] == cwd[differ_pos])
            {
                //temp=i;
                break;
            }
        }
        DEBUG_INFO("gama1: %d %d | %f %f %f\n",
        		   differ_pos,
        		   i,
	    		   gama1,
	    		   max_symbol_reliability,
	    		   chnl_rel_matrix[i][differ_pos]);
        gama1 = gama1 * max_symbol_reliability / chnl_rel_matrix[i][differ_pos];
    }
    for (i = 0; i < count1; i++)
    {
        same_pos = same_pos_vector[i];
        max_symbol_reliability = 0;
        sec_max_symbol_reliability = 0;
        max_row = -1;
        sec_max_row = -1;
        for (j = 0; j < GF_FIELD; j++)
        {
            if (chnl_rel_matrix[j][same_pos] > max_symbol_reliability)
            {
                sec_max_symbol_reliability = max_symbol_reliability;
                sec_max_row = max_row;
                max_symbol_reliability = chnl_rel_matrix[j][same_pos];
                max_row = j;
            }
            else if (chnl_rel_matrix[j][same_pos] > sec_max_symbol_reliability)
            {
                sec_max_symbol_reliability = chnl_rel_matrix[j][same_pos];
                sec_max_row = j;
            }
        }
        symbol_reliabilty_metric[i] = max_symbol_reliability / sec_max_symbol_reliability;
        temp_metric[i] = symbol_reliabilty_metric[i];
    }
    //insertion sort the metric from small to big
    for (j = 1; j < count1; j++)
    {
        temp = temp_metric[j];
        i = j - 1;
        while (i >= 0 && temp_metric[i] > temp)
        {
            temp_metric[i + 1] = temp_metric[i];
            i = i - 1;
        }
        temp_metric[i + 1] = temp;
    }
    gama2 = 1;
    //for (i = 0; i < (CODEWORD_LEN - MESSAGE_LEN + 1 - Hammin_dis); i++)
    for (i = 0; i < (code_min_dis - Hammin_dis); i++)
    {
        gama2 = gama2 * temp_metric[i];
        DEBUG_INFO("gama2: %f %f\n",
    		        gama2,
    		        temp_metric[i]);
    }
    DEBUG_INFO("gama: %f %f\n",
    		   gama1,
    		   gama2);
    if (gama1 <= gama2)
        return 1;
    else
        return 0;
}

int ReTMLcriterion(unsigned char est_cwd[], unsigned char recv_cwd[])
{
    long long i, j;
    float gama1, gama2;
    unsigned char same_pos_vector[CODEWORD_LEN], differ_pos_vector[CODEWORD_LEN];
    long long count1 = 0, count2 = 0;
    long long Hammin_dis, differ_pos, same_pos;
    float symbol_reliabilty_metric[CODEWORD_LEN];
    float temp, temp_metric[CODEWORD_LEN];
    float max_symbol_reliability, sec_max_symbol_reliability;
    long long max_row, sec_max_row;
    unsigned char cwd[CODEWORD_LEN];
    //long long code_min_dis = CODEWORD_LEN - MESSAGE_LEN + 1 - (GF_Q) * (GF_Q - 1) / 2;
    //long long ret_msg_k = (CODEWORD_LEN - (GF_Q) * (GF_Q - 1) / 2) / GF_Q * GF_Q;
    long long ret_msg_k = (long long)(ceil(((float)MESSAGE_LEN + (GF_Q) * (GF_Q - 1) / 2.0) / GF_Q) * GF_Q);
    long long code_min_dis = CODEWORD_LEN - ret_msg_k + 1 - (GF_Q) * (GF_Q - 1) / 2;
    DEBUG_NOTICE("ret_len: %ld %ld\n", ret_msg_k, code_min_dis);

    for (i = 0; i < CODEWORD_LEN; i++)
    {
        same_pos_vector[i] = 0xFF;
        differ_pos_vector[i] = 0xFF;
        temp_metric[i] = 0;
        cwd[i] = est_cwd[i];
    }

    for (i = 0; i < CODEWORD_LEN; i++)
    {
        if (cwd[i] == recv_cwd[i])
        {
            same_pos_vector[count1] = i;
            count1++;
        }
        else
        {
            differ_pos_vector[count2] = i;
            count2++;
        }
    }
    Hammin_dis = count2;
    gama1 = 1;
    for (j = 0; j < count2; j++)
    {
        differ_pos = differ_pos_vector[j];
        max_symbol_reliability = 0;
        for (i = 0; i < GF_FIELD; i++)
        {
            if (chnl_rel_matrix[i][differ_pos] > max_symbol_reliability)
                max_symbol_reliability = chnl_rel_matrix[i][differ_pos];
        }
        for (i = 0; i < GF_FIELD; i++)
        {
            if (power_polynomial_table[i][0] == cwd[differ_pos])
            {
                //temp=i;
                break;
            }
        }
        DEBUG_INFO("gama1: %d %d | %f %f %f\n",
        		   differ_pos,
        		   i,
	    		   gama1,
	    		   max_symbol_reliability,
	    		   chnl_rel_matrix[i][differ_pos]);
        gama1 = gama1 * max_symbol_reliability / chnl_rel_matrix[i][differ_pos];
    }
    for (i = 0; i < count1; i++)
    {
        same_pos = same_pos_vector[i];
        max_symbol_reliability = 0;
        sec_max_symbol_reliability = 0;
        max_row = -1;
        sec_max_row = -1;
        for (j = 0; j < GF_FIELD; j++)
        {
            if (chnl_rel_matrix[j][same_pos] > max_symbol_reliability)
            {
                sec_max_symbol_reliability = max_symbol_reliability;
                sec_max_row = max_row;
                max_symbol_reliability = chnl_rel_matrix[j][same_pos];
                max_row = j;
            }
            else if (chnl_rel_matrix[j][same_pos] > sec_max_symbol_reliability)
            {
                sec_max_symbol_reliability = chnl_rel_matrix[j][same_pos];
                sec_max_row = j;
            }
        }
        symbol_reliabilty_metric[i] = max_symbol_reliability / sec_max_symbol_reliability;
        temp_metric[i] = symbol_reliabilty_metric[i];
    }
    //insertion sort the metric from small to big
    for (j = 1; j < count1; j++)
    {
        temp = temp_metric[j];
        i = j - 1;
        while (i >= 0 && temp_metric[i] > temp)
        {
            temp_metric[i + 1] = temp_metric[i];
            i = i - 1;
        }
        temp_metric[i + 1] = temp;
    }
    gama2 = 1;
    //for (i = 0; i < (CODEWORD_LEN - MESSAGE_LEN + 1 - Hammin_dis); i++)
    for (i = 0; i < (code_min_dis - Hammin_dis); i++)
    {
        gama2 = gama2 * temp_metric[i];
        DEBUG_INFO("gama2: %f %f\n",
    		        gama2,
    		        temp_metric[i]);
    }
    DEBUG_INFO("gama: %f %f\n",
    		   gama1,
    		   gama2);
    if (gama1 <= gama2)
        return 1;
    else
        return 0;
}

int her_cmm_intp()
{
	long long i = 0, j = 0;

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		memset(cmm_intp_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	poly_init();

#if (1 == CFG_RET)
	koetter_intp_her_lcc(ret_tst_vct[0],
						 0,
						 0,//keep_cnt,//notice this, there are some problems
						 CODEWORD_LEN - 1 - ETA,
						 intp_poly_coef,
						 min_intp_poly);
#else
	koetter_intp_her_lcc(tst_vct[0],
						 0,
						 0,
						 CODEWORD_LEN - 1 - ETA,
						 intp_poly_coef,
						 min_intp_poly);
#endif

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		memcpy(cmm_intp_poly[i], intp_poly_coef[i], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	memcpy(cmm_intp_poly_degree, intp_poly_degree, sizeof(long long) * KOT_INTP_POLY_NUM);

#if (1 == TEST_MODE)
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != cmm_intp_poly[i][j])
			{
				DEBUG_NOTICE("cmm_intp_poly: %ld | %ld %ld %ld | %ld | %x\n",
							 i,
							 x_term_degree_table[j],
							 y_term_degree_table[j],
							 z_term_degree_table[j],
							 cmm_intp_poly_degree[i],
							 cmm_intp_poly[i][j]);
			}
		}
	}
#endif

	return 0;
}

int her_ucm_proc(long long tv_idx)
{
	long long i = 0;
	long long layer_idx = CODEWORD_LEN - 1 - ETA + 1;

	//start = clock();

	poly_init();
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		memcpy(intp_poly_coef[i], cmm_intp_poly[i], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	memcpy(intp_poly_degree, cmm_intp_poly_degree, sizeof(long long) * KOT_INTP_POLY_NUM);

	/*load*/
	layer_idx = kot_node_load(intp_poly_coef, tv_idx);
	layer_idx = layer_idx + (CODEWORD_LEN - 1 - ETA) + 1;
	DEBUG_NOTICE("layer_idx: %ld\n", layer_idx);
	//layer_idx = 0;

#if (1 == CFG_RET)
	koetter_intp_her_lcc(ret_tst_vct[tv_idx],
						 tv_idx,
						 layer_idx,
						 CODEWORD_LEN - 1,
						 intp_poly_coef,
						 min_intp_poly);
#else
	koetter_intp_her_lcc(tst_vct[tv_idx],
						 tv_idx,
						 layer_idx,
						 CODEWORD_LEN - 1,
						 intp_poly_coef,
						 min_intp_poly);
#endif

	//stop = clock();
	//runtime = runtime + (stop - start) / 1000.0000;

	her_fac(min_intp_poly, tv_est_msg[tv_idx], tv_est_cwd[tv_idx]);

#if (1 == TEST_MODE)
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("est_cwd: %ld | %ld | %x %x\n",
					 tv_idx,
					 i,
					 tv_est_cwd[tv_idx][i],
					 cwd_poly[i]);
	}
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_NOTICE("est_msg: %ld | %ld | %x %x\n",
					 tv_idx,
					 i,
					 tv_est_msg[tv_idx][i],
					 msg_poly[i]);
	}
#endif

	return 0;
}

int dump_err_info()
{
	return 0;
}

int her_lcc_dec()
{
	long long i = 0, j = 0;
	int ml_check_val = 0;
	long long tv_recv_diff_num = 0, tv_est_err_num = 0;
	long long tv_idx = 0;

#if (0 == CFG_BR)

	kot_node_clear();

	her_cmm_intp();

	intp_cnt = 0;
	for(i = 0; i < tst_vct_num; i++)
	{
		tv_idx = tst_vct_seq[i];
	
		her_ucm_proc(tv_idx);

		tv_dec_output_flag[tv_idx] = 1;

#if (1 == CFG_PRG_DECODING)
		if(0 == ml_check_val)
		{
			//her_cwd_check(tv_est_cwd[i]);
			ml_check_val = MLcriterion(tv_est_cwd[tv_idx], recv_poly);
			if(0 != ml_check_val)
			{
				tv_recv_diff_num = 0;
				tv_est_err_num = 0;
				for(j = 0; j < CODEWORD_LEN; j++)
				{
					if(recv_poly[j] != tv_est_cwd[tv_idx][j])
					{
						tv_recv_diff_num++;
					}
					if(cwd_poly[j] != tv_est_cwd[tv_idx][j])
					{
						tv_est_err_num++;
					}
				}
				DEBUG_NOTICE("ML pass: %ld %ld %ld %ld\n", tv_idx, tv_err[tv_idx], tv_recv_diff_num, tv_est_err_num);
				break;
			}
		}
#endif

	}
	DEBUG_NOTICE("intp_cnt: %ld\n", intp_cnt);

#else

	int popov_flag = 1;
	//long long genus = GF_Q * (GF_Q - 1) / 2;
	//long long radius = (CODEWORD_LEN - MESSAGE_LEN - genus + 1 - 1) / 2;
	long long ms_to_cnt = 0;
	//unsigned bug_flag = 0;
	br_poly_clear();

	for(i = 0; i < tst_vct_num; i++)
	{
		tv_idx = tst_vct_seq[i];
	
		//br_poly_clear();
		//bug_flag = 0;
		ms_to_cnt = 0;
		
		if(0 == tv_idx)
		{
			br_k_poly_construct(0);
		}
		else
		{
#if (0 == CFG_ACD_BR)		
			br_k_poly_ded(0, tv_idx);
#else
			br_k_poly_construct(tv_idx);
#endif			
		}
		
		br_m_poly_construct(tv_idx);
		br_v_matric_gen(tv_idx);

		//start = clock();

		popov_flag = is_popov_form(tv_idx);
		while(0 == popov_flag)
		{
			ms_reduction(tv_idx);
			popov_flag = is_popov_form(tv_idx);
			ms_to_cnt++;
			if(CODEWORD_LEN < ms_to_cnt)
			{
				DEBUG_SYS("br_ms_to: %ld\n", ms_to_cnt);
				break;
			}
		}

		//stop = clock();
		//runtime = runtime + (stop - start) / 1000.0000;

		br_g_poly_ret(tv_idx);
		br_q_poly_ret(tv_idx);

		her_fac(min_intp_poly, tv_est_msg[tv_idx], tv_est_cwd[tv_idx]);

		tv_dec_output_flag[tv_idx] = 1;

#if 0		
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			DEBUG_NOTICE("br_cwd_check: %ld %ld | %x %x\n", i, j, tv_est_cwd[i][j], cwd_poly[j]);
			if((tv_est_cwd[i][j] != cwd_poly[j])
				&& (radius >= tv_err[i]))
			{
				DEBUG_SYS("br_test_err: %ld | %ld\n", i, tv_err[i]);
				bug_flag = 1;
			}
		}
		if(1 == bug_flag)
		{
			DEBUG_SYS("Transmission over Channel:\n");
			for(i = 0; i < symbol_num; i++)
			{
				DEBUG_SYS("rx_symbol[%ld][0] = %f;\n", i, recv_seq[i][0]);
				DEBUG_SYS("rx_symbol[%ld][1] = %f;\n", i, recv_seq[i][1]);
			}
			for(i = 0; i < MESSAGE_LEN; i++)
			{
				DEBUG_SYS("msg_poly[%ld] = 0x%x;\n",
				          i,
				          msg_poly[i]);
			}
		}
#endif		
	}

#endif

	her_lcc_check_result();

	return 0;
}
