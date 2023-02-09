#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "debug_info.h"
#include "math.h"
#include "gf_cal.h"
#include "encoding.h"
#include "channel.h"
#include "interpolation.h"
#include "factorization.h"
#include "lcc_decoding.h"

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

int term_degree_table_init()
{
	long long i = 0;

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
#if (0 == CFG_DYM_SIZE)
		x_term_degree_table[i] = i / ((MAX_DEGREE + 1) * (MAX_DEGREE + 1));
		y_term_degree_table[i] = (i / (MAX_DEGREE + 1)) % (MAX_DEGREE + 1);
		z_term_degree_table[i] = i % (MAX_DEGREE + 1);
#else
		x_term_degree_table[i] = i / (Y_MAX_SIZE * Z_MAX_SIZE);
		y_term_degree_table[i] = (i / Z_MAX_SIZE) % Y_MAX_SIZE;
		z_term_degree_table[i] = i % Z_MAX_SIZE;	
#endif
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
	/*notice these degrees should be started form -1*/
	long long max_degree = -1, degree_val = -1;

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
	unsigned char z_flag = 0;

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
			
			/*convert poly to avoid high x-degree*/
			her_convert(intp_poly_tmp[j]);
		}

		min_idx = -1;
		min_degree = 65536;
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			z_flag = 0;
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
					if(0 != z_term_degree_table[k])
					{
						z_flag = 1;
					}
				}
			}
			
			if((min_degree > intp_poly_degree[j])
				&& (0 != z_flag))
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
	
#if 0
	long long zero_root_cnt = 0;
	unsigned char q1_eva_val = 0xFF;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		q1_eva_val = poly_eva_x_y(q1_poly_coef, af_pnt[i][0], af_pnt[i][1]);
		if(0xFF == q1_eva_val)
		{
			zero_root_cnt++;
		}
	}
	if(0 != zero_root_cnt)
	{
		DEBUG_SYS("zero_root_cnt: %ld\n", zero_root_cnt);
	}
#endif

	return 0;
}