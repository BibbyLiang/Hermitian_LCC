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
#include "re_encoding.h"
#include "lcc_decoding.h"

long long dev_to_cnt = 0;
unsigned char erasure_flag[CODEWORD_LEN];
unsigned char dev_to_flag[CODEWORD_LEN];

unsigned char store_q0_dev[MAX_POLY_TERM_SIZE];
unsigned char store_q1_dev[MAX_POLY_TERM_SIZE];
unsigned char store_v_dev[MAX_POLY_TERM_SIZE];

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
#if (0 == CFG_QUICK_POLY_SEARCH)		
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
#else
			j = term_search(pole_x_degree, pole_y_degree, 0);
			g_poly[j] = poly[i];
#endif

			poly[i] = 0xFF;
		}
	}

	//cnt_switch = 0;
	her_convert(poly);
	her_convert(g_poly);
	//cnt_switch = 1;
#if (1 == TEST_MODE)
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
#endif
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
#if (0 == CFG_QUICK_POLY_SEARCH)
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
#else
			j = term_search(x_degree, y_degree, 0);
			tmp_val = gf_multp(root_prev_coef, poly[i]);
			/*notice that it should be an addition operation*/
			poly[j] = gf_add(poly[j], tmp_val);
#endif
		}
	}

	//cnt_switch = 0;
	her_convert(poly);
	//cnt_switch = 1;
#if (1 == TEST_MODE)
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
#endif
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

#if 1
	her_encoding(est_msg, est_cwd_poly);
#endif

	return 0;	
}

int poly_dev_cal(unsigned char *poly, unsigned char *x_dev_poly, unsigned char *y_dev_poly)
{
	long long i = 0, j = 0;
	memset(x_dev_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(y_dev_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 1; i < MAX_POLY_TERM_SIZE; i++)
	{
		if((0 != (x_term_degree_table[i] % 2))
			&& (0xFF != poly[i]))
		{
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(((x_term_degree_table[j] + 1) == x_term_degree_table[i])
					&& (y_term_degree_table[j] == y_term_degree_table[i])
					&& (z_term_degree_table[j] == z_term_degree_table[i]))
				{
					break;
				}
			}
			x_dev_poly[j] = poly[i];
		}
		
		if((0 != (y_term_degree_table[i] % 2))
			&& (0xFF != poly[i]))
		{
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(((y_term_degree_table[j] + 1) == y_term_degree_table[i])
					&& (x_term_degree_table[j] == x_term_degree_table[i])
					&& (z_term_degree_table[j] == z_term_degree_table[i]))
				{
					break;
				}
			}
			y_dev_poly[j] = poly[i];
		}
	}
#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != x_dev_poly[i])
		{
			DEBUG_NOTICE("x_dev: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
			             y_term_degree_table[i],
			             z_term_degree_table[i],
			             x_dev_poly[i]);
		}
	}
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != y_dev_poly[i])
		{
			DEBUG_NOTICE("y_dev: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
			             y_term_degree_table[i],
			             z_term_degree_table[i],
			             y_dev_poly[i]);
		}
	}
#endif
	return 0;
}

int factorization_free()
{
	long long i = 0, j = 0;
	unsigned char q0_val = 0xFF, q1_val = 0xFF;
	/*warn that there may be stack overflow*/
	unsigned char q0_dev[MAX_POLY_TERM_SIZE], q1_dev[MAX_POLY_TERM_SIZE], tmp_dev[MAX_POLY_TERM_SIZE];
	memset(q0_dev, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(q1_dev, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(tmp_dev, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char x_dev_val = 0xFF, y_dev_val = 0xFF;
	long long dev_cnt = 0, dev_flag = 0;;

	dev_to_cnt = 0;
	memset(erasure_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
	memset(dev_to_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);

	/*q0 and q1 have been got in the last part of fun:poly_dev_test*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		q0_val = poly_eva_x_y(q0_poly_coef, af_pnt[i][0], af_pnt[i][1]);
		q1_val = poly_eva_x_y(q1_poly_coef, af_pnt[i][0], af_pnt[i][1]);

		if((0xFF == q0_val)
			&& (0xFF == q1_val))
		{
#if 0		
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(0xFF != q0_poly_coef[j])
				{
					DEBUG_NOTICE("q0_poly: %ld %ld | %x\n",
					             x_term_degree_table[j],
						         y_term_degree_table[j],
						         q0_poly_coef[j]);
				}
			}
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(0xFF != q1_poly_coef[j])
				{
					DEBUG_NOTICE("q1_poly: %ld %ld | %x\n",
					             x_term_degree_table[j],
						         y_term_degree_table[j],
						         q1_poly_coef[j]);
				}
			}
			poly_dev_cal(q0_poly_coef, x_dev, y_dev);
			x_dev_val = poly_eva_x_y(x_dev, af_pnt[i][0], af_pnt[i][1]);
			y_dev_val = poly_eva_x_y(y_dev, af_pnt[i][0], af_pnt[i][1]);
			q0_val = gf_add(x_dev_val, y_dev_val);
			DEBUG_NOTICE("q0_val_cal: %x %x %x\n",
			             x_dev_val,
			             y_dev_val,
			             q0_val);
			poly_dev_cal(q1_poly_coef, x_dev, y_dev);
			x_dev_val = poly_eva_x_y(x_dev, af_pnt[i][0], af_pnt[i][1]);
			y_dev_val = poly_eva_x_y(y_dev, af_pnt[i][0], af_pnt[i][1]);
			q1_val = gf_add(x_dev_val, y_dev_val);
			DEBUG_NOTICE("q1_val_cal: %x %x %x\n",
			             x_dev_val,
			             y_dev_val,
			             q1_val);
			est_cwd_poly[i] = gf_div(q0_val, q1_val);
#else

#if (1 == CFG_STORE_DEV)//store first-dev
			if(0 == dev_flag)
			{
				fac_dev_init();
				dev_flag = 1;
			}
			q0_val = poly_eva_x_y(store_q0_dev, af_pnt[i][0], af_pnt[i][1]);
			q1_val = poly_eva_x_y(store_q1_dev, af_pnt[i][0], af_pnt[i][1]);
			if((0xFF == q0_val)
				&& (0xFF == q1_val))
			{
				dev_to_cnt++;
				erasure_flag[i] = 1;
				dev_to_flag[i] = 1;
			}
#else
			dev_cnt = 0;
			memcpy(q0_dev, q0_poly_coef, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			memcpy(q1_dev, q1_poly_coef, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			while((0xFF == q0_val)
				  && (0xFF == q1_val))
			{
				poly_dev_build(q0_dev, tmp_dev);
				q0_val = poly_eva_x_y(tmp_dev, af_pnt[i][0], af_pnt[i][1]);
				memcpy(q0_dev, tmp_dev, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

				poly_dev_build(q1_dev, tmp_dev);
				q1_val = poly_eva_x_y(tmp_dev, af_pnt[i][0], af_pnt[i][1]);
				memcpy(q1_dev, tmp_dev, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				dev_cnt++;
				DEBUG_NOTICE("dev_cnt: %ld\n", dev_cnt);
				if(2 < dev_cnt)
				{
					DEBUG_NOTICE("dev_timeout\n");
					dev_to_cnt++;
					erasure_flag[i] = 1;
					dev_to_flag[i] = 1;
					break;
				}
			}
#endif

			est_cwd_poly[i] = gf_div(q0_val, q1_val);
#endif
		}
		else
		{
			est_cwd_poly[i] = gf_div(q0_val, q1_val);
		}
		DEBUG_NOTICE("est_cwd_poly: %ld | %x %x | %x\n",
		             i,
		             q0_val,
		             q1_val,
		             est_cwd_poly[i]);
	}
	
#if 0//test
	unsigned char dev_test_poly[MAX_POLY_TERM_SIZE];
	poly_dev_build(q0_poly_coef, dev_test_poly);
#endif

	if(0 != dev_to_cnt)
	{
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			if(1 == erasure_flag[i])
			{
				for(j = 0; j < CODEWORD_LEN; j++)
				{
					if(af_pnt[i][0] == af_pnt[j][0])
					{
						erasure_flag[j] = 1;
					}
				}
			}
		}

		DEBUG_IMPORTANT("dev_to_cnt: %ld\n", dev_to_cnt);
		fac_her_lagrange_poly_construct();
		fac_ret_poly_construct();
		fac_ret_encoding();
	}

	return 0;
}

unsigned char poly_term_dev_cal(long long term_idx, unsigned char term_coef, unsigned char x_val, unsigned char y_val)
{
	long long i = 0;
	long long x_degree = 0, y_degree = 0;
	unsigned char tmp_val_x_dev = 0xFF, tmp_val_y_dev = 0xFF, tmp_val = 0xFF;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char val = 0xFF;

	if(0 != (x_term_degree_table[term_idx] % 2))
	{
		x_degree = x_term_degree_table[term_idx] - 1;
		y_degree = y_term_degree_table[term_idx];
		if(0 <= x_degree)
		{
			tmp_val = gf_pow_cal(x_val, x_degree);
		}
		else
		{
			tmp_val = 0xFF;
		}
		tmp_val_x_dev = gf_multp(tmp_val, 
								 gf_pow_cal(y_val, y_degree));
	}
	if(0 != (y_term_degree_table[term_idx] % 2))
	{
		x_degree = x_term_degree_table[term_idx] + GF_Q;
		y_degree = y_term_degree_table[term_idx] - 1;
		tmp_val = gf_pow_cal(x_val, x_degree);
		if(0 <= y_degree)
		{
			tmp_val_y_dev = gf_multp(tmp_val, 
									 gf_pow_cal(y_val, y_degree));
		}
		else
		{
			tmp_val_y_dev = 0xFF;
		}
	}
	val = gf_add(tmp_val_x_dev, tmp_val_y_dev);
	val = gf_multp(val, term_coef);
	
	DEBUG_NOTICE("poly_term_dev_cal: %x %x | %x\n",
	             tmp_val_x_dev,
	             tmp_val_y_dev,
	             val);

	return val;
}

int poly_dev_build(unsigned char *poly, unsigned char *poly_dev)
{
	long long i = 0, j = 0;
	long long x_degree = 0, y_degree = 0;

	memset(poly_dev, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly[i])
		{
			if(0 != (x_term_degree_table[i] % 2))
			{
				x_degree = x_term_degree_table[i] - 1;
				y_degree = y_term_degree_table[i];
				for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
				{
					if((x_degree == x_term_degree_table[j])
						&& (y_degree == y_term_degree_table[j])
						&& (0 == z_term_degree_table[j]))
					{
						poly_dev[j] = gf_add(poly_dev[j], poly[i]);
#if 1
						DEBUG_NOTICE("poly_dev_cal_x: %ld %ld | %ld %ld | %x\n",
						             x_term_degree_table[i],
						             y_term_degree_table[i],
						             x_degree,
						             y_degree,
						             poly_dev[j]);
#endif						
						break;
					}
				}
			}
			if(0 != (y_term_degree_table[i] % 2))
			{
				x_degree = x_term_degree_table[i] + GF_Q;
				y_degree = y_term_degree_table[i] - 1;
				for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
				{
					if((x_degree == x_term_degree_table[j])
						&& (y_degree == y_term_degree_table[j])
						&& (0 == z_term_degree_table[j]))
					{
						poly_dev[j] = gf_add(poly_dev[j], poly[i]);
#if 1
						DEBUG_NOTICE("poly_dev_cal_y: %ld %ld | %ld %ld | %x\n",
						             x_term_degree_table[i],
						             y_term_degree_table[i],
						             x_degree,
						             y_degree,
						             poly_dev[j]);
#endif						
						break;
					}
				}
			}
		}
	}

	her_convert(poly_dev);
#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly_dev[i])
		{
			DEBUG_NOTICE("poly_dev: %ld %ld %ld | %x\n",
			              x_term_degree_table[i],
			              y_term_degree_table[i],
			              z_term_degree_table[i],
			              poly_dev[i]);
		}
	}
#endif
	return 0;
}

int fac_her_lagrange_poly_construct()
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char x_poly[MAX_POLY_TERM_SIZE], y_poly[MAX_POLY_TERM_SIZE];
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE], tmp_poly_y[MAX_POLY_TERM_SIZE];
	unsigned char y_her_curve_find_flag = 0;
	unsigned char x_div = 0x0, y_div = 0x0;
	unsigned char check_val = 0xFF;

	unsigned ret_x_flag = 0;

	for(i = 0; i < CODEWORD_LEN; i++)//af. point
	{
		if(1 == erasure_flag[i])
		{
			continue;
		}

		memset(lag_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
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
				if((0 == erasure_flag[k])
					&& (af_pnt[k][0] == power_polynomial_table[j][0]))
				{
					ret_x_flag = 1;
					break;
				}
			}

			if(af_pnt[i][0] == power_polynomial_table[j][0])
			{
#if 0			
				DEBUG_NOTICE("x_skip: %ld | %x | %x %x\n",
				             i,
				             power_polynomial_table[j][0],
				             af_pnt[i][0],
				             af_pnt[i][1]);
#endif				             
				//continue;
			}
#if 1
			else if(0 == ret_x_flag)
			{
				DEBUG_NOTICE("ret_x_flag_0: %ld | %x | %x %x\n",
				             i,
				             power_polynomial_table[j][0],
				             af_pnt[i][0],
				             af_pnt[i][1]);
			}
#endif
			else
			{
				memset(tmp_poly_x, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

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
#if 0				
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if((0xFF != x_poly[k])
						|| (0xFF != tmp_poly_x[k]))
					{
						DEBUG_NOTICE("x_poly_cal: %ld %ld | %ld %ld %ld | %x %x\n",
									  i,
									  j,
						              x_term_degree_table[k],
						              y_term_degree_table[k],
						              z_term_degree_table[k],
						              x_poly[k],
						              tmp_poly_x[k]);
					}
				}
#endif
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
#if 0
		check_val = gf_div(poly_eva_x_y(x_poly, af_pnt[i][0], af_pnt[i][1]), x_div);
		DEBUG_NOTICE("x_check_val: %ld | %x %x | %x\n", i, af_pnt[i][0], af_pnt[i][1], check_val);
		check_val = gf_div(poly_eva_x_y(y_poly, af_pnt[i][0], af_pnt[i][1]), y_div);
		DEBUG_NOTICE("y_check_val: %ld | %x %x | %x\n", i, af_pnt[i][0], af_pnt[i][1], check_val);
#endif		
		//memcpy(lag_poly[i], x_poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
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
								lag_poly[i][l] = gf_add(lag_poly[i][l],
														gf_multp(x_poly[k], y_poly[j]));
								lag_poly[i][l] = gf_div(lag_poly[i][l], x_div);
								lag_poly[i][l] = gf_div(lag_poly[i][l], y_div);
							}
						}
					}
				}
			}
		}
		//cnt_switch = 0;
		her_convert(lag_poly[i]);
		//cnt_switch = 1;
#if (1 == TEST_MODE)
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != lag_poly[i][j])
			{
				DEBUG_NOTICE("lag_poly: %ld | %ld %ld %ld | %x\n",
							  i,
				              x_term_degree_table[j],
				              y_term_degree_table[j],
				              z_term_degree_table[j],
				              lag_poly[i][j]);
			}
		}
#endif		
#if 0
		check_val = poly_eva_x_y(lag_poly[i], af_pnt[i][0], af_pnt[i][1]);
		DEBUG_NOTICE("check_val: %ld | %x %x | %x\n", i, af_pnt[i][0], af_pnt[i][1], check_val);
#endif
	}

	return 0;
}

int fac_ret_poly_construct()
{
	long long i = 0, j = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(ret_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
#if (0 == CFG_RET)
	memcpy(keep_sym, est_cwd_poly, sizeof(unsigned char) * CODEWORD_LEN);
#else
	memcpy(keep_sym, ret_est_cwd, sizeof(unsigned char) * CODEWORD_LEN);
#endif
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(1 == erasure_flag[i])
		{
			DEBUG_NOTICE("erasure_skip: %ld\n", i);
			continue;
		}
		else
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(0xFF != lag_poly[i][j])
				{
					tmp_poly[j] = gf_multp(keep_sym[i], lag_poly[i][j]);
					ret_poly[j] = gf_add(ret_poly[j], tmp_poly[j]);
					//cnt_switch = 0;
					her_convert(ret_poly);
					//cnt_switch = 1;
					DEBUG_NOTICE("ret_poly_cal: %ld | %x %x | %ld %ld %ld | %x\n",
								  i,
								  keep_sym[i],
								  lag_poly[i][j],
					              x_term_degree_table[j],
					              y_term_degree_table[j],
					              z_term_degree_table[j],
					              ret_poly[j]);
				}
			}
		}
	}
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
	return 0;
}

int fac_ret_encoding()
{
	long long i = 0, j = 0, k = 0;
	unsigned char x_val = 0xFF, y_val = 0xFF;;
	unsigned char val = 0xFF;
#if (0 == CFG_RET)
	memset(ret_cwd_poly, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
#endif
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == dev_to_flag[i])
		{
#if (0 == CFG_RET)
			ret_cwd_poly[i] = est_cwd_poly[i];
#else
			ret_est_cwd[i] = ret_est_cwd[i];
#endif
		}
		else
		{
#if (0 == CFG_RET)		
			ret_cwd_poly[i] = poly_eva_x_y(ret_poly, af_pnt[i][0], af_pnt[i][1]);
#else
			ret_est_cwd[i] = poly_eva_x_y(ret_poly, af_pnt[i][0], af_pnt[i][1]);
#endif
			DEBUG_NOTICE("ret_cwd_poly: %ld | %x %x %x\n",
			             i,
			             cwd_poly[i],
			             ret_est_cwd[i],
			             ret_cwd_poly[i]);
		}             
	}
#if 0
	/*check ret-cwd*/
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
							 gf_multp(ret_cwd_poly[k],
							 		  gf_multp(x_val, y_val)));
			}
			DEBUG_NOTICE("ret_cwd_check: %ld %ld | %x\n",
			             i,
			             j,
			             val);
		}
	}
#endif
#if (0 == CFG_RET)
	//memcpy(ret_est_cwd, ret_cwd_poly, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(est_cwd_poly, ret_cwd_poly, sizeof(unsigned char) * CODEWORD_LEN);
#endif
	return 0;
}

int fac_dev_init()
{
	memset(store_q0_dev, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(store_q1_dev, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(store_v_dev, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	poly_dev_build(q0_poly_coef, store_q0_dev);
	poly_dev_build(q1_poly_coef, store_q1_dev);
#if (1 == CFG_RET)	
	poly_dev_build(v_poly, store_v_dev);
#endif

#if 0//test
	DEBUG_NOTICE("test_poly_dev\n");
	test_poly_dev();
#endif

	return 0;
}