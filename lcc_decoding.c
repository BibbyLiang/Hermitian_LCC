#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "time.h"
#include "debug_info.h"
#include "math.h"
#include "gf_cal.h"
#include "encoding.h"
#include "channel.h"
#include "interpolation.h"
#include "factorization.h"
#include "re_encoding.h"
#include "lcc_decoding.h"

#define PI	3.1415926

clock_t start, stop;
float runtime;
long long err_cnt = 0;

/*col(locator, 0xff~GF_FIELD-2)-row(mesg, CODEWORD_LEN), same as matlab, contrary to most papers*/
float chnl_rel_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
float chnl_rel_order[CODEWORD_LEN];
long long chnl_rel_order_idx[CODEWORD_LEN];
long long chnl_rel_max_id[CODEWORD_LEN];
long long chnl_rel_scd_id[CODEWORD_LEN];
long long tst_vct_num = 0;
unsigned char **tst_vct;

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

	BubbleSort4(chnl_rel_order, CODEWORD_LEN, chnl_rel_order_idx);//¿É¿¿¶ÈÉýÐò
#if (1 == TEST_MODE)
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("chnl_rel_order_idx: %ld %f %ld\n",
					 i,
					 chnl_rel_order[i],
		             chnl_rel_order_idx[i]);
	}
#endif
	return 0;
}

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
		d0 = d0 / 1e2;
		d1 = d1 / 1e2;
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

	DEBUG_NOTICE("tst_vct_init OK\n");

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

	return 0;
}

int tst_vct_form()
{
	long long i = 0, j = 0, pickout_place = 0;
	long long cwd_idx = 0;
	unsigned char pickout_sym = 0xFF;

	for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			cwd_idx = chnl_rel_order_idx[j];

			if(ETA < j)
			{
				pickout_place = 0;
				pickout_sym = chnl_rel_max_id[cwd_idx];
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
			}
			DEBUG_NOTICE("tst_vct_cal: %ld %ld | %ld %ld | %ld\n",
			             i,
			             j,
			             cwd_idx,
			             pickout_place,
			             pickout_sym);
			tst_vct[i][cwd_idx] = power_polynomial_table[pickout_sym][0];
		}
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
		}
	}
#endif	

	return 0;
}

int koetter_interpolation_hermitian()
{
	poly_init();
#if (0 == CFG_RET)	
	poly_dev_test(recv_poly);
#else
	poly_dev_test(ret_trans_cwd);
#endif
#if (1 == CFG_FAC_FREE)
#if (1 == CFG_RET)
	ret_fac_free(q0_poly_coef, q1_poly_coef, v_poly);
	ret_cwd_recover();
#else
	factorization_free();
#endif
#else
#if (1 == CFG_RET)
	factorization_recur(min_intp_poly, ret_est_msg);
	ret_cwd_recover();
#else
	factorization_recur(min_intp_poly, est_msg_poly);
#endif	
#endif

	return 0;
}

int check_result_cwd(unsigned char *cwd, unsigned char *est_cwd)
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
			DEBUG_SYS("min_intp_idx: %ld\n", min_intp_idx);
		}
	}
	else
	{
		err_cnt++;
		genus = GF_Q * (GF_Q - 1) / 2;
		radius = (CODEWORD_LEN - MESSAGE_LEN - genus + 1 - 1) / 2;
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			if(cwd_poly[i] != recv_poly[i])
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
				DEBUG_SYS("check_result_cwd_err: %ld | %x %x %x\n",
			              i,
			              cwd[i],
			              est_cwd[i],
			              recv_poly[i]);
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
		}
	}

	return val;
}

int check_result_msg(unsigned char *msg, unsigned char *est_msg)
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
			DEBUG_SYS("min_intp_idx: %ld\n", min_intp_idx);
		}
	}
	else
	{
		err_cnt++;
		genus = GF_Q * (GF_Q - 1) / 2;
		radius = (CODEWORD_LEN - MESSAGE_LEN - genus + 1 - 1) / 2;
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			if(cwd_poly[i] != recv_poly[i])
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
