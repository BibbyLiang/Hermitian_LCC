#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "debug_info.h"
#include "math.h"
#include "gf_cal.h"
#include "encoding.h"
#include "channel.h"
#include "lcc_decoding.h"

#define PI	3.1415926

/*col(locator, 0xff~GF_FIELD-2)-row(mesg, CODEWORD_LEN), same as matlab, contrary to most papers*/
float chnl_rel_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
float chnl_rel_order[CODEWORD_LEN];
long long chnl_rel_order_idx[CODEWORD_LEN];
long long chnl_rel_max_id[CODEWORD_LEN];
long long chnl_rel_scd_id[CODEWORD_LEN];
long long tst_vct_num = 0;
unsigned char **tst_vct;

long long x_term_degree_table[MAX_POLY_TERM_SIZE];
long long y_term_degree_table[MAX_POLY_TERM_SIZE];
long long z_term_degree_table[MAX_POLY_TERM_SIZE];

unsigned char intp_poly_coef[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
unsigned char intp_poly_tmp[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
unsigned char q0_poly_coef[MAX_POLY_TERM_SIZE], q1_poly_coef[MAX_POLY_TERM_SIZE];
unsigned char intp_poly_degree[KOT_INTP_POLY_NUM];
unsigned char min_intp_poly[MAX_POLY_TERM_SIZE];

unsigned char poly_eva_x_y(unsigned char *poly_x_y, unsigned char x_val, unsigned char y_val)
{
	long long i = 0;
	unsigned char val = 0xFF, tmp_x = 0xFF, tmp_y = 0xFF, tmp_prod = 0x0;
	
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0 != z_term_degree_table[i])
		{
			continue;
		}
		if(0xFF == poly_x_y[i])
		{
			continue;
		}

		tmp_x = gf_pow_cal(x_val, x_term_degree_table[i]);
		tmp_y = gf_pow_cal(y_val, y_term_degree_table[i]);
		tmp_prod = gf_multp(tmp_x, tmp_y);
		tmp_prod = gf_multp(tmp_prod, poly_x_y[i]);
		val = gf_add(val, tmp_prod);
#if 0
		DEBUG_NOTICE("poly_eva_x_y_cal: %x %x | %ld %ld %ld | %x | %x\n",
				     x_val,
				     y_val,
				     x_term_degree_table[i],
				     y_term_degree_table[i],
				     z_term_degree_table[i],
				     poly_x_y[i],
	                 val);
#endif	                 
	}
#if 0
	DEBUG_NOTICE("poly_eva_x_y: %x %x | %x\n",
				 x_val,
				 y_val,
	             val);
#endif	             

	return val;
}

int poly_eva_x_y_test()
{
	/*for poly_eva_x_y test*/
	long long i = 0, j = 0;
	unsigned char tmp_msg[MAX_POLY_TERM_SIZE];
	memset(tmp_msg, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if((x_term_degree_table[j] == pole_basis_pow[i][0])
				&& (y_term_degree_table[j] == pole_basis_pow[i][1])
				&& (0 == z_term_degree_table[j]))
			{
				tmp_msg[j] = msg_poly[i];
				DEBUG_NOTICE("tmp_msg: %d | %d %d | %x\n",
				             j,
			                 x_term_degree_table[j],
			                 y_term_degree_table[j],
			                 tmp_msg[j]);
			}
		}
	}
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		poly_eva_x_y(tmp_msg, af_pnt[i][0], af_pnt[i][1]);
	}
	
	return 0;
}

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

	BubbleSort4(chnl_rel_order, CODEWORD_LEN, chnl_rel_order_idx);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("chnl_rel_order_idx: %ld %f %ld\n",
					 i,
					 chnl_rel_order[i],
		             chnl_rel_order_idx[i]);
	}

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
    chnl_rel_matrix[3][5] = 1;
    chnl_rel_matrix[5][1] = 1;
    chnl_rel_matrix[0][3] = 1;

	chnl_rel_matrix[1][6] = 0.9;
	chnl_rel_matrix[7][6] = 0.1;

	chnl_rel_matrix[3][4] = 0.8;
	chnl_rel_matrix[5][4] = 0.2;

	chnl_rel_matrix[7][2] = 0.7;
	chnl_rel_matrix[6][2] = 0.3;

	chnl_rel_matrix[3][0] = 0.6;
	chnl_rel_matrix[7][0] = 0.4;
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

int term_degree_table_init()
{
	long long i = 0;
	
	for(i = 0; i < ((MAX_DEGREE + 1) * (MAX_DEGREE + 1) * (MAX_DEGREE + 1)); i++)
	{
		x_term_degree_table[i] = i / ((MAX_DEGREE + 1) * (MAX_DEGREE + 1));
		y_term_degree_table[i] = (i / (MAX_DEGREE + 1)) % (MAX_DEGREE + 1);
		z_term_degree_table[i] = i % (MAX_DEGREE + 1);
#if 0		
		DEBUG_NOTICE("term_degree_table: %ld | %ld %ld %ld\n",
		             i,
		             x_term_degree_table[i],
		             y_term_degree_table[i],
		             z_term_degree_table[i]);
#endif		             
	}

	return 0;
}

long long poly_degree_cal(unsigned char *poly)
{
	long long i = 0;
	long long max_degree = 0, degree_val = 0;

	long long w_z = GF_Q * pole_basis_pow[MESSAGE_LEN - 1][0]
				  + (GF_Q + 1) * pole_basis_pow[MESSAGE_LEN - 1][1];

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly[i])
		{
			degree_val = GF_Q * x_term_degree_table[i]
			           + (GF_Q + 1) * y_term_degree_table[i]
			           + w_z * z_term_degree_table[i];
			if(degree_val > max_degree)
			{
				max_degree = degree_val;
			}
		}
	}
	//DEBUG_NOTICE("max_degree: %ld\n", max_degree);

	return max_degree;
}

long long poly_z_degree_get(unsigned char *poly)
{
	long long i = 0;
	long long max_degree = 0, degree_val = 0, z_degree = 0;

	long long w_z = GF_Q * pole_basis_pow[MESSAGE_LEN - 1][0]
				  + (GF_Q + 1) * pole_basis_pow[MESSAGE_LEN - 1][1];

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly[i])
		{
			degree_val = GF_Q * x_term_degree_table[i]
			           + (GF_Q + 1) * y_term_degree_table[i]
			           + w_z * z_term_degree_table[i];
			if(degree_val > max_degree)
			{
				max_degree = degree_val;
				z_degree = z_term_degree_table[i];
			}
		}
	}
	DEBUG_NOTICE("z_degree: %ld\n", z_degree);

	return z_degree;
}

/*compare the degrees of poly_1 and poly_2*/
unsigned char poly_degree_compare(unsigned char *poly_1, unsigned char *poly_2)
{
	long long i = 0;
	long long poly_degree_1 = 0, poly_degree_2 = 0;
	long long z_degree_1 = 0, z_degree_2 = 0;

	poly_degree_1 = poly_degree_cal(poly_1);
	poly_degree_2 = poly_degree_cal(poly_2);
	z_degree_1 = poly_z_degree_get(poly_1);
	z_degree_2 = poly_z_degree_get(poly_2);

	if(poly_degree_1 == poly_degree_2)
	{
		if(z_degree_1 > z_degree_2)
		{
			return 1;
		}
		else if(z_degree_1 < z_degree_2)
		{
			return 2;
		}
		else
		{
			return 0;
		}
	}
	else if(poly_degree_1 > poly_degree_2)
	{
		return 1;
	}
	else
	{
		return 2;
	}
}

int koetter_interpolation_hermitian()
{
	term_degree_table_init();
	poly_init();
	poly_dev_test(recv_poly);

	factorization_recur(min_intp_poly, est_msg_poly);
	factorization_free();
	
	return 0;
}

int poly_init()
{
	long long i = 0, j = 0;
	long long y_degree = 0, z_degree = 0;
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		memset(intp_poly_coef[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		y_degree = i % GF_Q;
		z_degree = i / GF_Q;
		
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if((0 == x_term_degree_table[j])
				&& (y_degree == y_term_degree_table[j])
				&& (z_degree == z_term_degree_table[j]))
			{
				intp_poly_coef[i][j] = 0x0;
				DEBUG_NOTICE("intp_poly_coef: %ld | %ld %ld %ld | %x\n",
				             i,
				             x_term_degree_table[j],
				             y_term_degree_table[j],
				             z_term_degree_table[j],
				             intp_poly_coef[i][j]);
				break;
			}
		}
		
		intp_poly_degree[i] = poly_degree_cal(intp_poly_coef[i]);
		memcpy(intp_poly_tmp[i], intp_poly_coef[i], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	memset(min_intp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	
	return 0;
}

int poly_q0q1_get(long long poly_idx)
{
	long long i = 0, j = 0;
	
	memset(q0_poly_coef, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(q1_poly_coef, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		//if(0 != z_term_degree_table[i])
		if(1 == z_term_degree_table[i])
		{
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if((x_term_degree_table[j] == x_term_degree_table[i])
					&& (y_term_degree_table[j] == y_term_degree_table[i])
					&& (0 == z_term_degree_table[j]))
				{
					break;
				}
			}
			q1_poly_coef[j] = intp_poly_coef[poly_idx][i];
		}
		else
		{
			q0_poly_coef[i] = intp_poly_coef[poly_idx][i];
		}
	}
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q0_poly_coef[i])
		{
			DEBUG_NOTICE("q0_poly: %ld %ld | %x\n",
			             x_term_degree_table[i],
				         y_term_degree_table[i],
				         q0_poly_coef[i]);
		}
	}
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q1_poly_coef[i])
		{
			DEBUG_NOTICE("q1_poly: %ld %ld | %x\n",
			             x_term_degree_table[i],
				         y_term_degree_table[i],
				         q1_poly_coef[i]);
		}
	}

	return 0;
}

int poly_normal_update(unsigned char *poly_tmp, unsigned char *poly_update, unsigned char *poly_min, unsigned char hs_dev_min, unsigned char hs_dev_self)
{
	long long i = 0;
	unsigned char poly_update_val = 0xFF, poly_min_val = 0xFF;

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		/*notice that these two vals should be clear*/
		poly_update_val = 0xFF, poly_min_val = 0xFF;

		/*is it different form the TCOM papaer?*/
		if(0xFF != poly_update[i])
		{
			//poly_update_val = gf_multp(hs_dev_min, poly_update[i]);
			poly_update_val = poly_update[i];
		}
		if(0xFF != poly_min[i])
		{
			poly_min_val = gf_multp(hs_dev_self, poly_min[i]);
			poly_min_val = gf_div(poly_min_val, hs_dev_min);
		}
		if(0xFF == poly_update_val)
		{
			poly_tmp[i] = poly_min_val;
		}
		else if(0xFF == poly_min_val)
		{
			poly_tmp[i] = poly_update_val;
		}
		else
		{
			poly_tmp[i] = gf_add(poly_update_val, poly_min_val);
		}
	}

	return 0;
}

int poly_min_update(unsigned char *poly_tmp, unsigned char *poly_min, unsigned char x_point_val, unsigned char hs_dev_min)
{
	long long i = 0, j = 0;
	/*warn that there may be stack overflow*/
	unsigned char x_up_poly_min[MAX_POLY_TERM_SIZE];
	memset(x_up_poly_min, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char poly_min_val = 0xFF;

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly_min[i])
		{
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if((x_term_degree_table[j] == (x_term_degree_table[i] + 1))
					&& (y_term_degree_table[j] == y_term_degree_table[i])
					&& (z_term_degree_table[j] == z_term_degree_table[i]))
				{
					x_up_poly_min[j] = poly_min[i];
					break;
				}
			}
			DEBUG_NOTICE("x_up_poly_min: %ld %ld %ld | %x\n",
			             x_term_degree_table[j],
			             y_term_degree_table[j],
			             z_term_degree_table[j],
			             x_up_poly_min[j]);
		}
	}
	
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly_min[i])
		{
			poly_min_val = gf_multp(poly_min[i], x_point_val);
		}
		else
		{
			/*notice that this val should be clear*/
			poly_min_val = 0xFF;
		}
		
		if((0xFF == x_up_poly_min[i])
			&& (0xFF == poly_min_val))
		{
			poly_tmp[i] = 0xFF;
		}
		else if((0xFF == x_up_poly_min[i])
			   && (0xFF != poly_min_val))
		{
			poly_tmp[i] = poly_min_val;
		}
		else if((0xFF != x_up_poly_min[i])
			   && (0xFF == poly_min_val))
		{
			poly_tmp[i] = x_up_poly_min[i];
		}
		else
		{
			poly_tmp[i] = gf_add(x_up_poly_min[i], poly_min_val);
		}
		
		if(0xFF != poly_tmp[i])
		{
			/*is it different form the TCOM papaer?*/
			//poly_min[i] = gf_multp(poly_min[i], hs_dev_min);
			DEBUG_NOTICE("poly_min: %ld %ld %ld | %x %x | %x\n",
			             x_term_degree_table[i],
			             y_term_degree_table[i],
			             z_term_degree_table[i],
			             poly_min_val,
			             x_up_poly_min[i],
			             poly_tmp[i]);
		}
	}

	return 0;
}

unsigned char hasse_dev_cal(long long point_idx, long long poly_idx, unsigned char test_sym)
{
	unsigned char val = 0xFF, q0_val = 0xFF, q1_val = 0xFF;

	poly_q0q1_get(poly_idx);
	q0_val = poly_eva_x_y(q0_poly_coef, af_pnt[point_idx][0], af_pnt[point_idx][1]);
	q1_val = poly_eva_x_y(q1_poly_coef, af_pnt[point_idx][0], af_pnt[point_idx][1]);
	q1_val = gf_multp(q1_val, test_sym);
	val = gf_add(q0_val, q1_val);
	DEBUG_NOTICE("hasse_dev: %ld | %ld | %x %x | %x | %x\n",
				 point_idx,
	             poly_idx,
	             q0_val,
	             q1_val,
	             test_sym,
	             val);

	return val;
}

int poly_dev_test(unsigned char *test_poly_seq)
{
	long long i = 0, j = 0, k = 0;
	unsigned char hasse_dev[KOT_INTP_POLY_NUM], q0_val = 0xFF, q1_val = 0xFF;
	memset(hasse_dev, 0xFF, sizeof(unsigned char) * KOT_INTP_POLY_NUM);
	long long min_idx = -1, min_degree = 65536, degree_val = 0, tmp_val = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("intp point: %ld | %x %x | %x\n",
		             i,
		             af_pnt[i][0],
		             af_pnt[i][1],
		             test_poly_seq[i]);
		/**notice that these two val should be clear*/
		min_idx = -1;
		min_degree = 65536;

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
			hasse_dev[j] = hasse_dev_cal(i, j, test_poly_seq[i]);
#endif
			
			/*find the min poly whose hasse dev is not zero*/
			if(0xFF != hasse_dev[j])
			{
				degree_val = poly_degree_cal(intp_poly_coef[j]);
				if(degree_val < min_degree)
				{
					min_degree = degree_val;
					min_idx = j;
				}
				if((degree_val == min_degree)
					&& (-1 != min_idx))
				{
					tmp_val = poly_degree_compare(intp_poly_coef[j], intp_poly_coef[min_idx]);
					if(1 == tmp_val)
					{
						min_degree = degree_val;
						min_idx = j;
					}
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
					poly_min_update(intp_poly_tmp[j], intp_poly_coef[j], af_pnt[i][0], hasse_dev[j]);
				}
				else
				{
					poly_normal_update(intp_poly_tmp[j], intp_poly_coef[j], intp_poly_coef[min_idx], hasse_dev[min_idx], hasse_dev[j]);
				}
			}
		}

		min_idx = -1;
		min_degree = 65536;
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			/*notice that the prev poly should be stored for updating properly*/
			memcpy(intp_poly_coef[j], intp_poly_tmp[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			intp_poly_degree[j] = poly_degree_cal(intp_poly_coef[j]);
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != intp_poly_coef[j][k])
				{
					DEBUG_NOTICE("intp_poly_coef_update: %ld | %ld %ld %ld | %ld | %x\n",
					             j,
					             x_term_degree_table[k],
					             y_term_degree_table[k],
					             z_term_degree_table[k],
					             intp_poly_degree[j],
					             intp_poly_coef[j][k]);
				}             
			}
			
			if(min_degree > intp_poly_degree[j])
			{
				min_degree = intp_poly_degree[j];
				min_idx = j;
			}
		}
	}

	memcpy(min_intp_poly, intp_poly_coef[min_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	poly_q0q1_get(min_idx);
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != min_intp_poly[i])
		{
			DEBUG_NOTICE("min_intp_poly: %ld | %ld %ld %ld | %x\n",
						 min_idx,
			             x_term_degree_table[i],
			             y_term_degree_table[i],
			             z_term_degree_table[i],
			             min_intp_poly[i]);
		}
	}

	return 0;
}

int g_poly_trans(unsigned char *poly, unsigned char *g_poly, long long pole_idx, long long z_degree)
{
	long long i = 0, j = 0;
	/*warn that there may be stack overflow*/
	long long pole_x_degree = 0, pole_y_degree = 0;

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if((0 != z_term_degree_table[i])
			&& (0xFF != poly[i]))
		{
			pole_x_degree = pole_basis_pow[pole_idx][0] * z_degree + x_term_degree_table[i];
			pole_y_degree = pole_basis_pow[pole_idx][1] * z_degree + y_term_degree_table[i];
		
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if((pole_x_degree == x_term_degree_table[j])
					&& (pole_y_degree == y_term_degree_table[j])
					&& (0 == z_term_degree_table[j]))//notice this
				{
					g_poly[j] = poly[i];
			        break;     		 
				}
			}

			poly[i] = 0xFF;
		}
	}

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != g_poly[i])
		{
			DEBUG_NOTICE("g_poly: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
			             y_term_degree_table[i],
			             z_term_degree_table[i],
			             g_poly[i]);
		}
	}
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly[i])
		{
			DEBUG_NOTICE("poly_rest: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
			             y_term_degree_table[i],
			             z_term_degree_table[i],
			             poly[i]);
		}
	}

	return 0;
}

long long lt_get(unsigned char *poly)
{
	long long i = 0;
	long long degree_val = 0, max_degree = -1, lt_idx = -1;
	
	long long w_z = GF_Q * pole_basis_pow[MESSAGE_LEN - 1][0]
				  + (GF_Q + 1) * pole_basis_pow[MESSAGE_LEN - 1][1];

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly[i])
		{
			degree_val = GF_Q * x_term_degree_table[i]
			           + (GF_Q + 1) * y_term_degree_table[i]
			           + w_z * z_term_degree_table[i];
			if(degree_val > max_degree)
			{
				max_degree = degree_val;
				lt_idx = i;
			}
		}
	}
	
	DEBUG_NOTICE("lt: %ld %ld %ld | %ld | %x\n",
	             x_term_degree_table[lt_idx],
	             y_term_degree_table[lt_idx],
	             z_term_degree_table[lt_idx],
	             max_degree,
	             poly[lt_idx]);
	
	return lt_idx;
}

long long chien_search_x_y(unsigned char *poly)
{
	long long i = 0;
	long long root_idx = -1;
	unsigned char val = 0xFF, x_val = 0xFF, y_val = 0xFF;

	for(i = 0; i < GF_FIELD; i++)
	{
		x_val = af_pnt[i][0];
		y_val = af_pnt[i][1];
		val = poly_eva_x_y(poly, x_val, y_val);
		/*notice that we assume that there is only one root*/
		if(0xFF == val)
		{
			root_idx = i;
			break;
		}
	}

	return root_idx;
}

unsigned char recur_root(unsigned char *poly, unsigned char *g_poly)
{
	long long i = 0;
	unsigned char root_val = 0xFF;
	long long poly_degree = -1, g_poly_degree = -1, poly_lt_idx = -1;

	poly_degree = poly_degree_cal(poly);
	g_poly_degree = poly_degree_cal(g_poly);
	if(poly_degree > g_poly_degree)
	{
		/*is there any root?*/
		root_val = 0xFF;
	}
	if(poly_degree < g_poly_degree)
	{
		root_val = 0xFF;
	}
	if(poly_degree == g_poly_degree)
	{
		poly_lt_idx = lt_get(poly);
		root_val = poly[poly_lt_idx];
		DEBUG_NOTICE("root_val_cal: %x %x\n",
		             root_val,
		             g_poly[poly_lt_idx]);
		root_val = gf_div(root_val, g_poly[poly_lt_idx]);
		//root_val = gf_add(poly[poly_lt_idx], g_poly[poly_lt_idx]);
	}

	DEBUG_NOTICE("root_val: %ld %ld | %x\n",
	             poly_degree,
	             g_poly_degree,
	             root_val);

	return root_val;
}

int recur_poly_update(unsigned char *poly, long long pole_idx, unsigned char root_prev_coef)
{
	long long i = 0, j = 0;
	long long x_degree = 0, y_degree = 0;
	unsigned char tmp_val = 0xFF;

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if((0xFF != poly[i])
			&& (0 != z_term_degree_table[i]))
		{
			x_degree = x_term_degree_table[i] + pole_basis_pow[pole_idx][0];
			y_degree = y_term_degree_table[i] + pole_basis_pow[pole_idx][1];

			DEBUG_NOTICE("poly_update_degree: %ld %ld | %ld %ld | %ld %ld\n",
						 x_term_degree_table[i],
						 y_term_degree_table[i],
						 pole_basis_pow[pole_idx][0],
						 pole_basis_pow[pole_idx][1],
			             x_degree,
			             y_degree);

			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if((x_degree == x_term_degree_table[j])
					&& (y_degree == y_term_degree_table[j])
					&& (0 == z_term_degree_table[j]))
				{
					tmp_val = gf_multp(root_prev_coef, poly[i]);
					/*notice that it should be an addition operation*/
					poly[j] = gf_add(poly[j], tmp_val);
				}
			}
		}
	}

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly[i])
		{
			DEBUG_NOTICE("recur_poly_update: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
			             y_term_degree_table[i],
			             z_term_degree_table[i],
			             poly[i]);
		}
	}

	return 0;
}

int factorization_recur(unsigned char *poly, unsigned char *est_msg)
{
	long long i = 0;
	/*warn that there may be stack overflow*/
	unsigned char poly_rest[MAX_POLY_TERM_SIZE];
	memcpy(poly_rest, poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char g_poly[MAX_POLY_TERM_SIZE];
	memset(g_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	long long z_degree = 0;

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		/*warn that we assume that there is only one z term*/
		if((0xFF != poly[i])
			&& (0 != z_term_degree_table[i]))
		{
			z_degree = z_term_degree_table[i];
			DEBUG_NOTICE("z_degree: %ld\n", z_degree);
			break;
		}
	}

	/*recur factorization refer to wanyq's thesis algorithm 2.3*/
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		memcpy(poly_rest, poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		memset(g_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		g_poly_trans(poly_rest, g_poly, MESSAGE_LEN - 1 - i, z_degree);
		est_msg[MESSAGE_LEN - 1 - i] = recur_root(poly_rest, g_poly);
		recur_poly_update(poly, MESSAGE_LEN - 1 - i, est_msg[MESSAGE_LEN - 1 - i]);
	}

	return 0;	
}

int factorization_free()
{
	long long i = 0;
	unsigned char q0_val = 0xFF, q1_val = 0xFF;

	/*q0 and q1 have been got in the last part of fun:poly_dev_test*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		q0_val = poly_eva_x_y(q0_poly_coef, af_pnt[i][0], af_pnt[i][1]);
		q1_val = poly_eva_x_y(q1_poly_coef, af_pnt[i][0], af_pnt[i][1]);
		est_cwd_poly[i] = gf_div(q0_val, q1_val);
		DEBUG_NOTICE("est_cwd_poly: %ld | %x %x | %x\n",
		             i,
		             q0_val,
		             q1_val,
		             est_cwd_poly[i]);
	}

	return 0;
}