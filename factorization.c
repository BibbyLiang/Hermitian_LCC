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

	her_convert(poly);
	her_convert(g_poly);

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

	her_convert(poly);

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
	
	return 0;
}

int factorization_free()
{
	long long i = 0, j = 0;
	unsigned char q0_val = 0xFF, q1_val = 0xFF;
	/*warn that there may be stack overflow*/
	unsigned char x_dev[MAX_POLY_TERM_SIZE], y_dev[MAX_POLY_TERM_SIZE];
	memset(x_dev, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(y_dev, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char x_dev_val = 0xFF, y_dev_val = 0xFF;

	/*q0 and q1 have been got in the last part of fun:poly_dev_test*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		q0_val = poly_eva_x_y(q0_poly_coef, af_pnt[i][0], af_pnt[i][1]);
		q1_val = poly_eva_x_y(q1_poly_coef, af_pnt[i][0], af_pnt[i][1]);

		if((0xFF == q0_val)
			&& (0xFF == q1_val))
		{
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

	return 0;
}