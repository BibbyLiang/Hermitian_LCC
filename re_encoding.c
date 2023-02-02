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

int her_lagrange_poly_construct()
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char x_poly[MAX_POLY_TERM_SIZE], y_poly[MAX_POLY_TERM_SIZE];
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE], tmp_poly_y[MAX_POLY_TERM_SIZE];
	unsigned char y_her_curve_find_flag = 0;

	for(i = 0; i < CODEWORD_LEN; i++)//af. point
	{
		memset(lag_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		memset(x_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		x_poly[0] = 0x0;
		memset(y_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		y_poly[0] = 0x0;
		for(j = 0; j < GF_FIELD; j++)//locator
		{
			if(af_pnt[i][0] == power_polynomial_table[j][0])
			{
				DEBUG_NOTICE("x_skip: %ld | %x | %x %x\n",
				             i,
				             power_polynomial_table[j][0],
				             af_pnt[i][0],
				             af_pnt[i][1]);
				//continue;
			}
			else
			{
				memset(tmp_poly_x, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				
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
								DEBUG_NOTICE("x_poly_cal_set: %ld %ld | %ld %ld %ld | %x %x\n",
											  i,
											  j,
								              x_term_degree_table[l],
								              y_term_degree_table[l],
								              z_term_degree_table[l],
								              x_poly[l],
								              x_poly[k]);
								break;
							}
						}
					}
				}
				
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
				DEBUG_NOTICE("y_skip: %ld | %x | %x %x\n",
				             i,
				             power_polynomial_table[j][0],
				             af_pnt[i][0],
				             af_pnt[i][1]);
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
						DEBUG_NOTICE("y_her_curve_find_flag: %ld | %x | %x %x\n",
						             i,
						             power_polynomial_table[j][0],
						             af_pnt[k][0],
						             af_pnt[k][1]);
						break;
					}
				}
				
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
							}
						}
					}
				}
			}
		}
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