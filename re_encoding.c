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

unsigned char lag_poly[CODEWORD_LEN][MAX_POLY_TERM_SIZE];
unsigned char ret_poly[MAX_POLY_TERM_SIZE];
unsigned char ret_cwd_poly[CODEWORD_LEN];
unsigned char keep_sym[CODEWORD_LEN];
long long keep_flag[CODEWORD_LEN];
unsigned char ret_trans_cwd[CODEWORD_LEN];
unsigned char ret_est_msg[MESSAGE_LEN];
unsigned char ret_est_cwd[CODEWORD_LEN];
unsigned char v_poly[MAX_POLY_TERM_SIZE];
unsigned char v_val[CODEWORD_LEN];

int her_lagrange_poly_construct()
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char x_poly[MAX_POLY_TERM_SIZE], y_poly[MAX_POLY_TERM_SIZE];
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE], tmp_poly_y[MAX_POLY_TERM_SIZE];
	unsigned char y_her_curve_find_flag = 0;
	unsigned char x_div = 0x0, y_div = 0x0;
	
	unsigned ret_x_flag = 0;

	memset(keep_flag, 0, sizeof(long long) * CODEWORD_LEN);
#if 1//(1 == TEST_MODE)	
	keep_flag[2] = 1;
	keep_flag[3] = 1;
	keep_flag[4] = 1;
	keep_flag[5] = 1;
#endif	

	for(i = 0; i < CODEWORD_LEN; i++)//af. point
	{
#if 1			
		if(0 == keep_flag[i])
		{
			continue;
		}
#endif	

		memset(lag_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		memset(x_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		x_poly[0] = 0x0;
		memset(y_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		y_poly[0] = 0x0;
		x_div = 0x0;
		y_div = 0x0;
		
		for(j = 0; j < GF_FIELD; j++)//locator
		{
#if 1
			ret_x_flag = 0;
			for(k = 0; k < CODEWORD_LEN; k++)
			{
				if((1 == keep_flag[k])
					&& (af_pnt[k][0] == power_polynomial_table[j][0]))
				{
					ret_x_flag = 1;
					break;
				}
			}
#endif
			if(af_pnt[i][0] == power_polynomial_table[j][0])
			{
#if 1
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
#if 0
						DEBUG_NOTICE("x_poly_cal_process: %ld %ld | %ld %ld %ld | %x %x\n",
									  i,
									  j,
						              x_term_degree_table[k],
						              y_term_degree_table[k],
						              z_term_degree_table[k],
						              x_poly[k],
						              tmp_poly_x[k]);
#endif
						for(l = 0; l < MAX_POLY_TERM_SIZE; l++)
						{
#if 0						
							DEBUG_NOTICE("x_poly_cal_find: %ld %ld | %ld %ld %ld | %x %x\n",
										  i,
										  j,
							              x_term_degree_table[l],
							              y_term_degree_table[l],
							              z_term_degree_table[l],
							              x_poly[l],
							              x_poly[k]);
#endif							              
							if((x_term_degree_table[l] == (x_term_degree_table[k] + 1))
								&& (y_term_degree_table[l] == y_term_degree_table[k])
								&& (z_term_degree_table[l] == z_term_degree_table[k]))
							{
								x_poly[l] = x_poly[k];
								x_poly[k] = 0xFF;
#if 0								
								DEBUG_NOTICE("x_poly_cal_set: %ld %ld | %ld %ld %ld | %x %x\n",
											  i,
											  j,
								              x_term_degree_table[l],
								              y_term_degree_table[l],
								              z_term_degree_table[l],
								              x_poly[l],
								              x_poly[k]);
#endif
								break;
							}
						}
					}
				}
#if 1
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if(0xFF != x_poly[k])
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
#if 1
				DEBUG_NOTICE("y_skip: %ld | %x | %x %x\n",
				             i,
				             power_polynomial_table[j][0],
				             af_pnt[i][0],
				             af_pnt[i][1]);
#endif				             
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
#if 1
						DEBUG_NOTICE("y_her_curve_find_flag: %ld | %x | %x %x\n",
						             i,
						             power_polynomial_table[j][0],
						             af_pnt[k][0],
						             af_pnt[k][1]);
#endif						
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
#if 1
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if(0xFF != y_poly[k])
					{
						DEBUG_NOTICE("y_poly_cal: %ld %ld | %ld %ld %ld | %x %x\n",
									  i,
									  j,
						              x_term_degree_table[k],
						              y_term_degree_table[k],
						              z_term_degree_table[k],
						              y_poly[k],
						              tmp_poly_y[k]);
					}
				}
#endif				
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
#if 1
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != x_poly[j])
			{
				DEBUG_NOTICE("x_poly: %ld | %ld %ld %ld | %x\n",
							  i,
				              x_term_degree_table[j],
				              y_term_degree_table[j],
				              z_term_degree_table[j],
				              x_poly[j]);
			}
		}
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != y_poly[j])
			{
				DEBUG_NOTICE("y_poly: %ld | %ld %ld %ld | %x\n",
							  i,
				              x_term_degree_table[j],
				              y_term_degree_table[j],
				              z_term_degree_table[j],
				              y_poly[j]);
			}
		}
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
		
		her_convert(lag_poly[i]);
		
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
	}

	return 0;
}

int ret_poly_construct()
{
	long long i = 0, j = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(ret_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	memcpy(keep_sym, recv_poly, sizeof(unsigned char) * CODEWORD_LEN);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == keep_flag[i])
		{
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

	return 0;
}

int ret_encoding()
{
	long long i = 0, j = 0, k = 0;
	unsigned char x_val = 0xFF, y_val = 0xFF;;
	unsigned char val = 0xFF;

	memset(ret_cwd_poly, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		ret_cwd_poly[i] = poly_eva_x_y(ret_poly, af_pnt[i][0], af_pnt[i][1]);
		DEBUG_NOTICE("ret_cwd_poly: %ld | %x\n",
		             i,
		             ret_cwd_poly[i]);
	}
	
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

	return 0;
}

int ret_trans()
{
	long long i = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		ret_trans_cwd[i] = gf_add(recv_poly[i], ret_cwd_poly[i]);
		DEBUG_NOTICE("ret_trans_cwd: %ld | %x %x | %x\n",
		             i,
		             recv_poly[i],
		             ret_cwd_poly[i],
		             ret_trans_cwd[i]);
	}

	return 0;
}

int ret_cwd_recover()
{
	long long i = 0;

	memset(ret_est_cwd, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	her_encoding(ret_est_msg, ret_est_cwd);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		est_cwd_poly[i] = gf_add(ret_est_cwd[i], ret_cwd_poly[i]);
		DEBUG_NOTICE("ret_est_cwd: %ld | %x %x | %x\n",
		             i,
		             ret_est_cwd[i],
		             ret_cwd_poly[i],
		             est_cwd_poly[i]);
	}

	return 0;
}

int v_poly_construct()
{
	long long i = 0, j = 0, k = 0;
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE];
	unsigned char keep_gf_flag[GF_FIELD];
	memset(keep_gf_flag, 1, sizeof(unsigned char) * GF_FIELD);
	
	for(i = 0; i < GF_FIELD; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if((af_pnt[j][0] == power_polynomial_table[i][0])
				&& (0 == keep_flag[j]))
			{
				keep_gf_flag[i] = 0;
				break;
			}
		}
	}

	memset(v_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	v_poly[0] = 0x0;
	memset(v_val, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < GF_FIELD; i++)
	{
		if(0 == keep_gf_flag[i])
		{
			continue;
		}
		else
		{
			memset(tmp_poly_x, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			
			for(j = (MAX_POLY_TERM_SIZE - 1); j >= 0; j--)
			{
				if(0xFF != v_poly[j])
				{
					tmp_poly_x[j] = gf_multp(v_poly[j], power_polynomial_table[i][0]);

					for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
					{
						if((x_term_degree_table[k] == (x_term_degree_table[j] + 1))
							&& (y_term_degree_table[k] == y_term_degree_table[j])
							&& (z_term_degree_table[k] == z_term_degree_table[j]))
						{
							v_poly[k] = v_poly[j];
							v_poly[j] = 0xFF;
							break;
						}
					}
				}
			}
			
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(0xFF == tmp_poly_x[j])
				{
					continue;
				}
				else if((0xFF == v_poly[j])
						&& (0xFF != tmp_poly_x[j]))
				{
					v_poly[j] = tmp_poly_x[j];
				}
				else
				{
					v_poly[j] = gf_add(v_poly[j], tmp_poly_x[j]);
				}
			}
		}
	}

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != v_poly[i])
		{
			DEBUG_NOTICE("v_poly: %ld %ld %ld | %x\n",
			              x_term_degree_table[i],
			              y_term_degree_table[i],
			              z_term_degree_table[i],
			              v_poly[i]);
		}
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		v_val[i] = poly_eva_x_y(v_poly, af_pnt[i][0], af_pnt[i][1]);
		DEBUG_NOTICE("v_val: %ld | %x %x | %x\n",
					 i,
					 af_pnt[i][0],
					 af_pnt[i][1],
		             v_val[i]);
	}

	return 0;
}

int keep_position_set(long long *keep_poition)
{
	long long i = 0;



	return 0;
}