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

#define K_PLUS_G	CODEWORD_LEN//(MESSAGE_LEN + (GF_Q) * (GF_Q - 1) / 2)

long long sys_ret_z_deg = 0;

long long sys_intp_seq[K_PLUS_G];
unsigned char sys_intp_point[K_PLUS_G][2];
unsigned char sys_intp_poly[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];
unsigned char sys_poly_tmp[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];
unsigned char sys_q0[MAX_POLY_TERM_SIZE], sys_q1[MAX_POLY_TERM_SIZE];
unsigned char sys_prev_poly[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];

unsigned char sys_proc_col_flag[CODEWORD_LEN];
unsigned char **sys_proc_mat;
unsigned char **sys_sqr_mat;
unsigned char sys_gen_mat[MESSAGE_LEN][CODEWORD_LEN];
unsigned char sys_tmp_gau_mat[MESSAGE_LEN][CODEWORD_LEN];
long long this_new_col_idx = 0;

unsigned char aj_factor_flag[GF_FIELD][CFG_SYS_GEN_POLY];
unsigned char aj_factor_flag_store[GF_FIELD][CFG_SYS_GEN_POLY];
unsigned char aj_factor_a_set_flag[CFG_SYS_GEN_POLY];

unsigned char sys_elm_poly[CODEWORD_LEN][MAX_POLY_TERM_SIZE];
unsigned char sys_act_poly[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];
unsigned char sys_act_poly_store[CFG_SYS_GEN_POLY][MAX_POLY_TERM_SIZE];

unsigned char sys_ret_poly[MAX_POLY_TERM_SIZE];

unsigned char sys_ret_cwd[CODEWORD_LEN];

unsigned char sys_intp_min_poly[MAX_POLY_TERM_SIZE];
unsigned char sys_mul_poly[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
unsigned char sys_break_mat_poly[KOT_INTP_POLY_NUM][KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];

unsigned char sys_g_store_val[KOT_INTP_POLY_NUM][KOT_INTP_POLY_NUM][CODEWORD_LEN];
long long sys_g_store_deg[KOT_INTP_POLY_NUM][KOT_INTP_POLY_NUM];

unsigned char sys_intp_coef_mat[KOT_INTP_POLY_NUM][KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
//unsigned char sys_intp_coef_mat_tmp[KOT_INTP_POLY_NUM][KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
/*note that the arrays have overflowed, then change them frome 3D to 2D*/
unsigned char sys_intp_coef_mat_tmp[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];
unsigned char sys_intp_coef_mat_base_tmp[KOT_INTP_POLY_NUM][MAX_POLY_TERM_SIZE];

unsigned char sys_tv_trans[CODEWORD_LEN];

int sys_intp_point_gen(long long *rel_order_idx)
{
	long long i = 0;
	long long idx = 0;

	for(i = 0; i < K_PLUS_G; i++)
	{
		memset(sys_intp_point[i], 0xFF, sizeof(unsigned char) * 2);
	}
	memset(sys_intp_seq, 0, sizeof(long long) * K_PLUS_G);

#if 1
	for(i = 0; i < K_PLUS_G; i++)
	{
		idx = rel_order_idx[CODEWORD_LEN - 1 - i];
		sys_intp_point[i][0] = af_pnt[idx][0];
		sys_intp_point[i][1] = af_pnt[idx][1];
		
		sys_intp_seq[i] = idx;
		
		DEBUG_NOTICE("sys_intp_point: %ld | %ld | %x %x\n",
		             i,
		             idx,
		             sys_intp_point[i][0],
		             sys_intp_point[i][1]);

		if(MESSAGE_LEN > i)
		{
			//sys_proc_col_flag[idx] = 1;
		}
	}
#else//for information set selection test
	for(i = 0; i < K_PLUS_G; i++)
	{
		idx = i;
	
		sys_intp_point[i][0] = af_pnt[idx][0];
		sys_intp_point[i][1] = af_pnt[idx][1];

		sys_intp_seq[i] = idx;

		DEBUG_NOTICE("sys_intp_point: %ld | %ld | %x %x\n",
		             i,
		             idx,
		             sys_intp_point[i][0],
		             sys_intp_point[i][1]);

		if(MESSAGE_LEN > i)
		{
			//sys_proc_col_flag[K_PLUS_G - 1 - i] = 1;
		}
	}
#endif

	return 0;
}

int sys_intp_poly_init()
{
	long long i = 0, j = 0;

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		memset(sys_intp_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		memset(sys_prev_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		j = term_search(0, i, 0);
		sys_intp_poly[i][j] = 0x0;

#if (1 == TEST_MODE)
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != sys_intp_poly[i][j])
			{
				DEBUG_NOTICE("sys_intp_poly: %ld %ld %ld | %x\n",
				             x_term_degree_table[j],
					         y_term_degree_table[j],
					         z_term_degree_table[j],
					         sys_intp_poly[i][j]);
			}
		}
#endif		
	}

	return 0;
}

int sys_q0q1_get(unsigned char *poly)
{
	long long i = 0, j = 0;
	
	memset(sys_q0, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(sys_q1, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if((0 != y_term_degree_table[i])
			&& (0 == z_term_degree_table[i]))
		{
			j = term_search(x_term_degree_table[i], y_term_degree_table[i], 0);
			sys_q1[j] = poly[i];
		}
		else
		{
			sys_q0[i] = poly[i];
		}
	}

#if (1 == TEST_MODE)	
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != sys_q0[i])
		{
			DEBUG_NOTICE("sys_q0_poly: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
				         y_term_degree_table[i],
				         z_term_degree_table[i],
				         sys_q0[i]);
		}
	}
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != sys_q1[i])
		{
			DEBUG_NOTICE("sys_q1_poly: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
				         y_term_degree_table[i],
				         z_term_degree_table[i],
				         sys_q1[i]);
		}
	}
#endif

	return 0;
}

unsigned char sys_hasse_dev_cal(unsigned char *poly, unsigned char *input_point)
{
	unsigned char val = 0xFF;

	val = poly_eva_x_y(poly, input_point[0], input_point[1]);

	return val;
}

int sys_poly_update(unsigned char *poly_tmp,//output in this iter
                         unsigned char *poly_update,//to be update
                         unsigned char *poly_min,//min poly
                         unsigned char hs_dev_min,
                         unsigned char hs_dev_self,
                         unsigned char x_point_val,
                         unsigned char hs_flag)
{
	long long i = 0, j = 0;
	/*warn that there may be stack overflow*/
	unsigned char x_up_poly_min[MAX_POLY_TERM_SIZE];
	memset(x_up_poly_min, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char poly_update_val = 0xFF, poly_min_val = 0xFF;

	if(1 == hs_flag)
	{
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0 != z_term_degree_table[i])
			{
				continue;
			}
		
			if(0xFF != poly_min[i])
			{
				j = term_search((x_term_degree_table[i] + 1), y_term_degree_table[i], z_term_degree_table[i]);
				x_up_poly_min[j] = poly_min[i];
#if 0				
				DEBUG_NOTICE("x_up_poly_min: %ld %ld %ld | %x\n",
				             x_term_degree_table[j],
				             y_term_degree_table[j],
				             z_term_degree_table[j],
				             x_up_poly_min[j]);
#endif
			}
		}
		
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0 != z_term_degree_table[i])
			{
				continue;
			}
		
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
#if 0//(1 == TEST_MODE)		
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
#endif		
		}
	}
	else
	{
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0 != z_term_degree_table[i])
			{
				continue;
			}
		
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
	}

#if 0//(1 == TEST_MODE)	
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly_tmp[i])
		{
			DEBUG_NOTICE("poly_tmp: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
			             y_term_degree_table[i],
			             z_term_degree_table[i],
			             poly_tmp[i]);
		}
	}
#endif

	return 0;
}

int sys_has_dev_test()
{
	long long i = 0, j = 0;
	unsigned dev_val = 0xFF;

	for(i = 0; i < K_PLUS_G; i++)
	{
		if(0 == sys_proc_col_flag[i])
		{
			continue;
		}

		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			dev_val = sys_hasse_dev_cal(sys_intp_poly[j], af_pnt[i]);
			DEBUG_NOTICE("sys_has_dev_test: %ld | %ld | %x %x | %x\n",
			             i,
			             j,
			             af_pnt[i][0],
			             af_pnt[i][1],
			             dev_val);
		}
	}

	return 0;
}

int sys_mat_init()
{
	long long i = 0, j = 0, k = 0;
	unsigned char val = 0xFF;
	unsigned char term_tmp[MAX_POLY_TERM_SIZE];

	memset(sys_proc_col_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		memset(sys_gen_mat[i], 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
		memset(sys_tmp_gau_mat[i], 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	}

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			memset(term_tmp, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			k = term_search(pole_basis_pow[i][0], pole_basis_pow[i][1], 0);
			term_tmp[k] = 0x0;
			val = poly_eva_x_y(term_tmp, af_pnt[j][0], af_pnt[j][1]);
			sys_gen_mat[i][j] = val;
			DEBUG_NOTICE("sys_gen_mat: %ld %ld | %x\n", i, j, sys_gen_mat[i][j]);
		}
	}

	return 0;
}

int sys_par_col_init(long long row, long long col)
{
	long long i = 0;

	sys_proc_mat = (unsigned char**)malloc(sizeof(unsigned char*) * row);
	for (i = 0; i < row; i++)
	{
		sys_proc_mat[i] = (unsigned char*)malloc(sizeof(unsigned char) * col);
	}

	sys_sqr_mat = (unsigned char**)malloc(sizeof(unsigned char*) * row);
	for (i = 0; i < row; i++)
	{
		sys_sqr_mat[i] = (unsigned char*)malloc(sizeof(unsigned char) * row);
	}

	return 0;
}

int sys_par_col_exit(long long row, long long col)
{
	long long i = 0;

	for (i = 0; i < row; i++)
	{
  		free(sys_proc_mat[i]);
		sys_proc_mat[i] = NULL;
  	}
	free(sys_proc_mat);
	sys_proc_mat = NULL;

	for (i = 0; i < row; i++)
	{
  		free(sys_sqr_mat[i]);
		sys_sqr_mat[i] = NULL;
  	}
	free(sys_sqr_mat);
	sys_sqr_mat = NULL;

	return 0;
}

int sys_par_col_gen(long long row, long long col)
{
	long long i = 0, j = 0;
	long long row_cnt = 0, col_cnt = 0;

	for (i = 0; i < row; i++)
	{
		memset(sys_proc_mat[i], 0xFF, sizeof(unsigned char) * col);
	}

#if (0 == CFG_FAST_FULL_GAU)
	col_cnt = 0;
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 != sys_proc_col_flag[i])
		{
			row_cnt = 0;
			for(j = (MESSAGE_LEN - row); j < MESSAGE_LEN; j++)
			{
				sys_proc_mat[row_cnt][col_cnt] = sys_gen_mat[j][i];
				//DEBUG_NOTICE("sys_proc_mat_cal: %ld %ld | %ld %ld | %x\n", i, j, col_cnt, row_cnt, sys_proc_mat[row_cnt][col_cnt]);
				row_cnt++;
			}
			col_cnt++;
			if(col <= col_cnt)
			{
				break;
			}
		}
	}
#else
	if((col <= ((MESSAGE_LEN - (GF_Q) * (GF_Q - 1) / 2) + 1))
		&& (row == MESSAGE_LEN))
	{
		for(i = 0; i < MESSAGE_LEN; i++)
		{
			memcpy(sys_proc_mat[i], sys_tmp_gau_mat[i], sizeof(unsigned char) * (col - 1));
		}
		for(i = 0; i < MESSAGE_LEN; i++)
		{
			sys_proc_mat[i][col - 1] = sys_gen_mat[i][this_new_col_idx];
		}
	}
	else
	{
		col_cnt = 0;
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			if(0 != sys_proc_col_flag[i])
			{
				row_cnt = 0;
				for(j = (MESSAGE_LEN - row); j < MESSAGE_LEN; j++)
				{
					sys_proc_mat[row_cnt][col_cnt] = sys_gen_mat[j][i];
					row_cnt++;
				}
				col_cnt++;
				if(col <= col_cnt)
				{
					break;
				}
			}
		}
	}
#endif

#if 0//(1 == TEST_MODE)
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			DEBUG_NOTICE("sys_proc_mat: %ld %ld | %x\n", i, j, sys_proc_mat[i][j]);
		}
	}
#endif

	return 0;
}

int sys_par_mat_trans(long long row, long long col)
{
	long long i = 0, j = 0, k = 0;
	unsigned char val = 0xFF;

	for (i = 0; i < row; i++)
	{
		memset(sys_sqr_mat[i], 0xFF, sizeof(unsigned char) * row);
	}

	for(i = 0; i < row; i++)
	{
		for(j = 0; j < row; j++)
		{
			val = 0xFF;
			for(k = 0; k < col; k++)
			{
				val = gf_add(val,
				             gf_multp(sys_proc_mat[i][k], sys_proc_mat[j][k]));
				//DEBUG_NOTICE("sys_sqr_mat_cal: %ld %ld %ld | %x %x\n", i, j, k, sys_proc_mat[i][k], sys_proc_mat[j][k]);
			}
			sys_sqr_mat[i][j] = val;
			DEBUG_NOTICE("sys_sqr_mat: %ld %ld | %x\n", i, j, sys_sqr_mat[i][j]);
		}
	}

	return 0;
}

int sys_par_mat_gau_elm(long long row)
{
	long long i = 0, j = 0, k = 0;
	unsigned char div_factor = 0xFF;
	int val = 0, rank = 0;

	for(i = 0; i < row; i++)
	{
		for(j = 0; j < row; j++)
		{
			if(j == i)
			{
				continue;
			}

			div_factor = gf_div(sys_sqr_mat[i][i], sys_sqr_mat[j][i]);
			for(k = 0; k < row; k++)
			{
				sys_sqr_mat[j][k] = gf_add(sys_sqr_mat[j][k],
				                           gf_div(sys_sqr_mat[i][k], div_factor));
				DEBUG_NOTICE("sys_sqr_mat_gau_elm_cal: %ld %ld %ld | %x\n", i, j, k, sys_sqr_mat[i][j]);
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < row; j++)
		{
			DEBUG_NOTICE("sys_sqr_mat_gau_elm: %ld %ld | %x\n", i, j, sys_sqr_mat[i][j]);
		}
	}
#endif

	for(j = 0; j < row; j++)
	{
		for(i = 0; i < row; i++)
		{
			if(0xFF != sys_sqr_mat[i][j])
			{
				rank++;
				break;
			}
		}
	}

	if(row <= rank)
	{
		val = 0;
	}
	else
	{
		val = 1;
	}

	return val;
}

int sys_proc_mat_gau_elm(long long row, long long col)
{
	long long i = 0, j = 0, k = 0;
	long long sel_row_idx = 0;
	unsigned char div_factor = 0xFF;
	int val = 0, rank = 0;

	for(i = 0; i < col; i++)
	{
		for(j = 0; j < row; j++)
		{
			if(0xFF != sys_proc_mat[j][i])
			{
				sel_row_idx = j;
				break;
			}
		}

		if(j >= row)
		{
			continue;
		}

		for(j = 0; j < col; j++)
		{
			if((i == j)
				|| (0xFF == sys_proc_mat[sel_row_idx][j]))
			{
				continue;
			}

			div_factor = gf_div(sys_proc_mat[sel_row_idx][i], sys_proc_mat[sel_row_idx][j]);
			for(k = 0; k < row; k++)
			{
				sys_proc_mat[k][j] = gf_add(sys_proc_mat[k][j],
				                            gf_div(sys_proc_mat[k][i], div_factor));
				DEBUG_NOTICE("sys_sqr_mat_gau_elm_cal: %ld %ld %ld | %ld | %x | %x\n", i, j, k, sel_row_idx, div_factor, sys_proc_mat[k][j]);
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			DEBUG_NOTICE("sys_proc_mat_gau_elm: %ld %ld | %x\n", i, j, sys_proc_mat[i][j]);
		}
	}
#endif

	for(j = 0; j < col; j++)
	{
		for(i = 0; i < row; i++)
		{
			if(0xFF != sys_proc_mat[i][j])
			{
				rank++;
				break;
			}
		}
	}

	if(col <= rank)
	{
		val = 0;
	}
	else
	{
		val = 1;
	}

	return val;
}

/*this function apply to sys_proc_mat*/
int sys_mat_gau_elm(long long row, long long col)
{
	long long i = 0, j = 0, k = 0;
	long long sel_row_idx = 0;
	unsigned char div_factor = 0xFF;
	int val = 0, rank = 0, zero_col = 0;;

	for(i = 0; i < col; i++)
	{
		for(j = 0; j < row; j++)
		{
			if(0xFF != sys_proc_mat[j][i])
			{
				sel_row_idx = j;
				break;
			}
		}

		if(j >= row)
		{
			continue;
		}

		for(j = 0; j < col; j++)
		{
			if((i == j)
				|| (0xFF == sys_proc_mat[sel_row_idx][j]))
			{
				continue;
			}

			div_factor = gf_div(sys_proc_mat[sel_row_idx][i], sys_proc_mat[sel_row_idx][j]);
			for(k = 0; k < row; k++)
			{
				sys_proc_mat[k][j] = gf_add(sys_proc_mat[k][j],
				                     gf_div(sys_proc_mat[k][i], div_factor));
				//DEBUG_NOTICE("sys_mat_gau_elm_cal: %ld %ld %ld | %ld | %x | %x\n", i, j, k, sel_row_idx, div_factor, sys_proc_mat[k][j]);
			}
		}
#if (1 == CFG_ZERO_COL_CHECK_GAU_ELM)/*fast full gau based on zero-col check*/
		//if(MESSAGE_LEN == row)
		{
			for(j = 0; j < row; j++)
			{
				if(0xFF != sys_proc_mat[j][i])
				{
					break;
				}
			}
			if(row <= j)
			{
				zero_col++;
			}
			if(0 < zero_col)
			{
				val = 1;
				return val;
			}
		}
#endif		
	}

#if (1 == TEST_MODE)
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			//DEBUG_NOTICE("sys_mat_gau_elm: %ld %ld | %x\n", i, j, sys_proc_mat[i][j]);
		}
	}
#endif

	for(j = 0; j < col; j++)
	{
		for(i = 0; i < row; i++)
		{
			if(0xFF != sys_proc_mat[i][j])
			{
				//DEBUG_NOTICE("sel_col_non_zero_idx: %ld %ld\n", i, j);
				rank++;
				break;
			}
		}
	}

	if(col <= rank)
	{
		val = 0;
	}
	else
	{
		val = 1;
	}
	DEBUG_NOTICE("rank: %ld %ld\n", col, rank);

	return val;
}

int sys_check_poly_w(long long *poly_deg)
{
	long long i = 0;
	long long max_deg = 0, min_deg = 65536;
	int val = 0;

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		if(max_deg < poly_deg[i])
		{
			max_deg = poly_deg[i];
		}
		if(min_deg > poly_deg[i])
		{
			min_deg = poly_deg[i];
		}
	}

	if(GF_Q >= (max_deg - min_deg))
	{
		val = 0;
	}
	else
	{
		val = 1;
	}

	return val;
}

int sys_check_part_dep(long long col)
{
	long long i = 0;
	int val = 0;

	sys_par_col_init(MESSAGE_LEN - col, col);
	sys_par_col_gen(MESSAGE_LEN - col, col);
	//sys_par_mat_trans(MESSAGE_LEN - col, col);
	//val = sys_par_mat_gau_elm(MESSAGE_LEN - col);
	//val = sys_proc_mat_gau_elm(MESSAGE_LEN - col, col);
	val = sys_mat_gau_elm(MESSAGE_LEN - col, col);
	sys_par_col_exit(MESSAGE_LEN - col, col);

	return val;
}

int sys_check_full_dep(long long col)
{
	long long i = 0;
	int val = 0;

	sys_par_col_init(MESSAGE_LEN, col);
	sys_par_col_gen(MESSAGE_LEN, col);
	val = sys_mat_gau_elm(MESSAGE_LEN, col);
#if (1 == CFG_FAST_FULL_GAU)
	if(0 == val)
	{
		for(i = 0; i < MESSAGE_LEN; i++)
		{
			memcpy(sys_tmp_gau_mat[i], sys_proc_mat[i], sizeof(unsigned char) * col);
		}
	}
#endif	
	sys_par_col_exit(MESSAGE_LEN, col);

	return val;
}

int check_aj_factor(long long locator_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char val = 0xFF;

	memset(aj_factor_flag[locator_idx], 0, sizeof(unsigned char) * CFG_SYS_GEN_POLY);

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{	
		aj_factor_flag[locator_idx][i] = 1;
		for(j = 1; j < GF_Q; j++)
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0 != z_term_degree_table[k])
				{
					continue;
				}
			
				if((0xFF != sys_intp_poly[i][k])
					&& (j == y_term_degree_table[k])
					&& (0 == z_term_degree_table[k]))
				{
					l = term_search(x_term_degree_table[k], 0, 0);
					//tmp_poly[l] = gf_add(tmp_poly[l], sys_intp_poly[i][k]);
					tmp_poly[l] = sys_intp_poly[i][k];
				}
			}

#if (1 == TEST_MODE)
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != tmp_poly[k])
				{
					DEBUG_NOTICE("tmp_poly: %ld %ld %ld | %x\n",
								 x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 tmp_poly[k]);
				}
			}
#endif

			val = poly_eva_x_y(tmp_poly, power_polynomial_table[locator_idx][0], 0xFF);
			if(0xFF != val)
			{
				aj_factor_flag[locator_idx][i] = 0;
				break;
			}
		}

		DEBUG_NOTICE("aj_factor_flag: %x %d | %d\n", power_polynomial_table[locator_idx][0], i, aj_factor_flag[locator_idx][i]);
	}

	return 0;
}

int poly_aj_factor_div(unsigned char poly_seq[][MAX_POLY_TERM_SIZE], long long poly_idx, long long locator_idx)
{
	long long i = 0, j = 0;
	unsigned char q0_tmp[MAX_POLY_TERM_SIZE];
	memset(q0_tmp, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char res_tmp[MAX_POLY_TERM_SIZE];
	memset(res_tmp, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char all_clean_flag = 0;

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0 != z_term_degree_table[i])
		{
			continue;
		}

		if((0xFF != poly_seq[poly_idx][i])
			&& (0 == z_term_degree_table[i]))
		{
			j = term_search(x_term_degree_table[i], y_term_degree_table[i], 0);
			q0_tmp[j] = poly_seq[poly_idx][i];
		}
	}
	
#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q0_tmp[i])
		{
			DEBUG_NOTICE("q0_tmp: %ld %ld %ld | %x\n",
						 x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 q0_tmp[i]);
		}
	}
#endif

	all_clean_flag = 0;
	memset(res_tmp, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	while(0 == all_clean_flag)
	{
		for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
		{
			if(0 != z_term_degree_table[i])
			{
				continue;
			}

			if(0xFF != q0_tmp[i])
			{
#if 1
				DEBUG_NOTICE("q0_tmp_res: %ld %ld %ld | %x\n",
							 x_term_degree_table[i],
							 y_term_degree_table[i],
							 z_term_degree_table[i],
							 q0_tmp[i]);
#endif				
				all_clean_flag = 0;
				break;
			}
		}
		if(MAX_POLY_TERM_SIZE <= i)
		{
			all_clean_flag = 1;
			break;
		}

		for(i = (MAX_POLY_TERM_SIZE - 1); i >= 0; i--)
		{
			if(0 != z_term_degree_table[i])
			{
				continue;
			}

			if(0xFF != q0_tmp[i])
			{
				j = term_search(x_term_degree_table[i] - 1, y_term_degree_table[i], 0);
				res_tmp[j] = gf_add(q0_tmp[i], res_tmp[j]);

				q0_tmp[i] = 0xFF;
				q0_tmp[j] = gf_add(q0_tmp[j],
								   gf_multp(res_tmp[j], power_polynomial_table[locator_idx][0]));
#if (1 == TEST_MODE)
				DEBUG_NOTICE("q0_tmp_res_cal: %ld %ld %ld | %x || %ld %ld %ld | %x || %ld %ld %ld | %x\n",
							 x_term_degree_table[i],
							 y_term_degree_table[i],
							 z_term_degree_table[i],
							 q0_tmp[i],
							 x_term_degree_table[j],
							 y_term_degree_table[j],
							 z_term_degree_table[j],
							 q0_tmp[j],
							 x_term_degree_table[j],
							 y_term_degree_table[j],
							 z_term_degree_table[j],
							 res_tmp[j]);
#endif			
			}
		}
	}
	//memset(q0_tmp, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memcpy(q0_tmp, res_tmp, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q0_tmp[i])
		{
			DEBUG_NOTICE("q0_tmp_div: %x | %ld %ld %ld | %x\n",
						 power_polynomial_table[locator_idx][0],
						 x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 q0_tmp[i]);
		}
	}
#endif

	//memset(sys_act_poly[poly_idx], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memcpy(sys_act_poly[poly_idx], q0_tmp, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != sys_act_poly[poly_idx][i])
		{
			DEBUG_NOTICE("sys_act_poly_div: %x | %d | %ld %ld %ld | %x\n",
						 power_polynomial_table[locator_idx][0],
						 poly_idx,
						 x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 sys_act_poly[poly_idx][i]);
		}
	}
#endif

	return 0;
}

int poly_aj_factor_trans(long long locator_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char aj_tmp_flag[CFG_SYS_GEN_POLY];
	memset(aj_tmp_flag, 0, sizeof(unsigned char) * CFG_SYS_GEN_POLY);
	unsigned char aj_tmp_b_flag = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	unsigned char tmp_val[CFG_SYS_GEN_POLY];
	memset(tmp_val, 0xFF, sizeof(unsigned char) * CFG_SYS_GEN_POLY);
	long long l_star_idx = 0;
	long long deg_star = 65536;
	long long deg_tmp[CFG_SYS_GEN_POLY];
	memset(deg_tmp, 0, sizeof(long long) * CFG_SYS_GEN_POLY);

	memset(aj_factor_a_set_flag, 0, sizeof(unsigned char) * CFG_SYS_GEN_POLY);

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		memset(aj_tmp_flag, 0, sizeof(unsigned char) * CFG_SYS_GEN_POLY);
		aj_tmp_b_flag = 0;
		l_star_idx = 0;
	    deg_star = 65536;
	    memset(deg_tmp, 0, sizeof(long long) * CFG_SYS_GEN_POLY);
	    memset(tmp_val, 0xFF, sizeof(unsigned char) * CFG_SYS_GEN_POLY);

		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(1 == aj_factor_a_set_flag[j])
			{
				/*this aj has been trans.*/
				continue;
			}

			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0 != z_term_degree_table[k])
				{
					continue;
				}
			
				if((0xFF != sys_act_poly[j][k])
					&& ((CFG_SYS_GEN_POLY - 1 - i) == y_term_degree_table[k])
					&& (0 == z_term_degree_table[k]))
				{
					l = term_search(x_term_degree_table[k], 0, 0);
					tmp_poly[l] = sys_act_poly[j][k];
				}
			}

			tmp_val[j] = poly_eva_x_y(tmp_poly, power_polynomial_table[locator_idx][0], 0xFF);
			if(0xFF != tmp_val[j])
			{
				aj_tmp_flag[j] = 1;
				aj_tmp_b_flag = 1;
			}
			//deg_tmp[j] = poly_degree_cal(tmp_poly);
			deg_tmp[j] = poly_degree_cal(sys_act_poly[j]);
			if((deg_star > deg_tmp[j])
				&& (0 < deg_tmp[j])
				&& (0xFF != tmp_val[j]))
			{
				DEBUG_NOTICE("l_star: %ld %ld | %ld | %d\n", deg_star, deg_tmp[j], j, aj_tmp_flag[j]);
				deg_star = deg_tmp[j];
				l_star_idx = j;
			}
		}

		//memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);	

#if (1 == TEST_MODE)
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			DEBUG_NOTICE("b_trans: %ld | %ld | %ld | %d | %x\n",
						  j,
						  l_star_idx,
						  aj_tmp_b_flag,
						  aj_tmp_flag[j],
						  tmp_val[j]);
		}
#endif

		if(0 != aj_tmp_b_flag)
		{
			aj_factor_a_set_flag[l_star_idx] = 1;

			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				if(0 == aj_tmp_flag[j])
				{
					continue;
				}
				if(j == l_star_idx)
				{
					continue;
				}
				
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if(0 != z_term_degree_table[k])
					{
						continue;
					}
				
#if 0				
					sys_intp_poly[j][k] = gf_add(gf_multp(tmp_val[l_star_idx], sys_intp_poly[j][k]),
					                             gf_multp(tmp_val[j], sys_intp_poly[l_star_idx][k]));
#endif
					sys_act_poly[j][k] = gf_add(gf_multp(tmp_val[l_star_idx], sys_act_poly[j][k]),
					                     		gf_multp(tmp_val[j], sys_act_poly[l_star_idx][k]));
				}
				//aj_factor_flag[locator_idx][j] = 1;
			}
		}

#if (1 == TEST_MODE)
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			deg_tmp[j] = poly_degree_cal(sys_act_poly[j]);
		
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_act_poly[j][k])
				{
					DEBUG_NOTICE("sys_intp_poly_trans: %d | %ld | %ld | %ld %ld %ld | %x\n",
								 locator_idx,
								 j,
								 deg_tmp[j],
					             x_term_degree_table[k],
						         y_term_degree_table[k],
						         z_term_degree_table[k],
						         sys_act_poly[j][k]);
				}
			}
		}
#endif
	}

	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		aj_factor_flag[locator_idx][i] = 1;
		aj_tmp_b_flag = 0;

		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0 != z_term_degree_table[k])
				{
					continue;
				}

				if((0xFF != sys_act_poly[i][k])
					&& (j == y_term_degree_table[k])
					&& (0 == z_term_degree_table[k]))
				{
					l = term_search(x_term_degree_table[k], 0, 0);
					tmp_poly[l] = sys_act_poly[i][k];
				}
			}
#if (1 == TEST_MODE)
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != tmp_poly[k])
				{
					DEBUG_NOTICE("tmp_poly: %ld %ld %ld | %x\n",
					             x_term_degree_table[k],
						         y_term_degree_table[k],
						         z_term_degree_table[k],
						         tmp_poly[k]);
				}
			}
#endif
			tmp_val[j] = poly_eva_x_y(tmp_poly, power_polynomial_table[locator_idx][0], 0xFF);
			DEBUG_NOTICE("tmp_val: %x %x\n", power_polynomial_table[locator_idx][0], tmp_val[j]);
			
			if(0xFF != tmp_val[j])
			{
				aj_tmp_b_flag = 1;
			}
			
			if(1 == aj_tmp_b_flag)
			{
				break;
			}
		}
		
		if(1 == aj_tmp_b_flag)
		{
			aj_factor_flag[locator_idx][i] = 0;
		}
#if (1 == TEST_MODE)		
		deg_tmp[i] = poly_degree_cal(sys_act_poly[i]);
#endif
		DEBUG_NOTICE("sys_intp_poly_trans_result: %ld %ld | %ld | %ld\n",
					 locator_idx,
					 i,
					 deg_tmp[i],
					 aj_factor_flag[locator_idx][i]);
	}

	return 0;
}

int sys_elm_intp_poly_proc(unsigned char poly_seq[][MAX_POLY_TERM_SIZE], long long pnt_idx)
{
	long long i = 0, j = 0;
	long long sel_poly_idx = 0;
	unsigned char eva_val = 0xFF, min_eva_val = 0xFF;;
	long long poly_deg[CFG_SYS_GEN_POLY], min_deg = 65536;
	memset(poly_deg, 0, sizeof(long long) * CFG_SYS_GEN_POLY);

	min_deg = 65536;
	memset(poly_deg, 0, sizeof(long long) * CFG_SYS_GEN_POLY);
	for(j = 0; j < CFG_SYS_GEN_POLY; j++)
	{
		poly_deg[j] = poly_degree_cal(poly_seq[j]);
		if(poly_deg[j] <= min_deg)
		{
#if 0		
			eva_val = sys_hasse_dev_cal(poly_seq[j], af_pnt[pnt_idx]);
			if(0xFF != eva_val)
			{
				min_deg = poly_deg[j];
				sel_poly_idx = j;
				min_eva_val = eva_val;
			}
#else
			min_deg = poly_deg[j];
			sel_poly_idx = j;
#endif			
		}
	}
	min_eva_val = sys_hasse_dev_cal(poly_seq[sel_poly_idx], af_pnt[pnt_idx]);

	memcpy(sys_elm_poly[pnt_idx], poly_seq[sel_poly_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
	{
		if(0 != z_term_degree_table[j])
		{
			continue;
		}

		if(0xFF != sys_elm_poly[pnt_idx][j])
		{
			sys_elm_poly[pnt_idx][j] = gf_div(sys_elm_poly[pnt_idx][j], min_eva_val);
			DEBUG_NOTICE("sys_elm_poly: %d | %ld | %x | %ld %ld %ld | %x\n",
						 pnt_idx,
						 sel_poly_idx,
						 min_eva_val,
			             x_term_degree_table[j],
				         y_term_degree_table[j],
				         z_term_degree_table[j],
				         sys_elm_poly[pnt_idx][j]);
		}
	}
#if (1 == TEST_MODE)
	eva_val = sys_hasse_dev_cal(sys_elm_poly[pnt_idx], af_pnt[pnt_idx]);
	DEBUG_NOTICE("sys_elm_poly_sel_pos: %ld | %x | %ld\n", pnt_idx, eva_val, min_deg);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if((0 == sys_proc_col_flag[i])
			|| (pnt_idx == i))
		{
			continue;
		}

		eva_val = sys_hasse_dev_cal(sys_elm_poly[pnt_idx], af_pnt[i]);
		DEBUG_NOTICE("sys_elm_poly_other_pos: %ld | %x\n", i, eva_val);
	}
#endif
	return 0;
}

int sys_elm_poly_gen()
{
	long long i = 0, j = 0, k = 0;
	long long pnt_idx = 0, locator_idx = 0;
	long long sel_poly_idx = 0;
	long long val = 0;
	long long deg_tmp[CFG_SYS_GEN_POLY];
	memset(deg_tmp, 0, sizeof(long long) * CFG_SYS_GEN_POLY);
	unsigned char this_gf_x_trans_flag = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		memset(sys_elm_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	for(i = 0; i < GF_FIELD; i++)
	{
		val = 0;
		for(j = 0; j < GF_Q; j++)
		{
			if(1 == sys_proc_col_flag[i * GF_Q + j])
			{
				val = 1;
				break;
			}
		}
		if(1 == val)
		{
			check_aj_factor(i);
		}
	}

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == (i % GF_Q))
		{
			this_gf_x_trans_flag = 0;
		}

		if(0 == sys_proc_col_flag[i])
		{
			continue;
		}

		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			memcpy(sys_act_poly[j], sys_intp_poly[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}

		for(j = 0; j < GF_FIELD; j++)
		{
			if(power_polynomial_table[j][0] == af_pnt[i][0])
			{
				locator_idx = j;
				break;
			}
		}

		val = 0;
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(1 == aj_factor_flag[locator_idx][j])
			{
				val = 1;
				break;
			}
		}

#if 0
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(0 == aj_factor_flag[locator_idx][j])
			{
				poly_aj_factor_trans(locator_idx);
				break;
			}
		}
#else
		if(0 == this_gf_x_trans_flag)
		{
			DEBUG_NOTICE("trans: %ld\n", i);
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				if(0 == aj_factor_flag[locator_idx][j])
				{
					poly_aj_factor_trans(locator_idx);
					break;
				}
			}

			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				memcpy(sys_act_poly_store[j], sys_act_poly[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
			for(j = 0; j < GF_FIELD; j++)
			{
				memcpy(aj_factor_flag_store[j], aj_factor_flag[j], sizeof(unsigned char) * CFG_SYS_GEN_POLY);
			}
			
			//this_gf_x_trans_flag = 1;
		}
		else
		{
			DEBUG_NOTICE("no_trans: %ld\n", i);
#if 0
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				memcpy(sys_act_poly[j], sys_act_poly_store[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
#endif			
			for(j = 0; j < GF_FIELD; j++)
			{
				memcpy(aj_factor_flag[j], aj_factor_flag_store[j], sizeof(unsigned char) * CFG_SYS_GEN_POLY);
			}
		}
#endif

#if 0
		val = 65536;
		sel_poly_idx = 0;
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(1 == aj_factor_flag[locator_idx][j])
			{
				deg_tmp[j] = poly_degree_cal(sys_act_poly[j]);
				if((val > deg_tmp[j])
					&& (0 < deg_tmp[j]))
				{
					val = deg_tmp[j];
					sel_poly_idx = j;
				}
			}
		}
#endif		

#if 0
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(1 == aj_factor_flag[locator_idx][j])
			{
				poly_aj_factor_div(sys_act_poly, j, locator_idx);
			}
		}
#else
		if(0 == this_gf_x_trans_flag)
		{
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				if(1 == aj_factor_flag[locator_idx][j])
				{
					poly_aj_factor_div(sys_act_poly, j, locator_idx);
				}
			}
			
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				memcpy(sys_act_poly_store[j], sys_act_poly[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
			
			this_gf_x_trans_flag = 1;
		}
		else
		{
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				memcpy(sys_act_poly[j], sys_act_poly_store[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
		}
#endif

		for(k = 0; k < CODEWORD_LEN; k++)
		{
			if((k == i)
				|| (af_pnt[k][0] != af_pnt[i][0])
				|| (0 == sys_proc_col_flag[k]))
			{
				continue;
			}
			sys_kot_intp(sys_act_poly, af_pnt[k]);
		}

#if (1 == TEST_MODE)
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_act_poly[j][k])
				{
					DEBUG_NOTICE("sys_intp_poly_final_get: %d %d | %ld %ld %ld | %x\n",
								 i,
								 j,
					             x_term_degree_table[k],
						         y_term_degree_table[k],
						         z_term_degree_table[k],
						         sys_act_poly[j][k]);
				}
			}
		}
#endif
		sys_elm_intp_poly_proc(sys_act_poly, i);
	}

	return 0;
}

int sys_ret_poly_gen()
{
	long long i = 0, j = 0;
	unsigned char eva_val = 0xFF;
	int val = 0;
	long long cwd_dim = CODEWORD_LEN;

	memset(sys_ret_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == sys_proc_col_flag[i])
		{
			continue;
		}

		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0 != z_term_degree_table[j])
			{
				continue;
			}

			if(0xFF != sys_elm_poly[i][j])
			{
				sys_ret_poly[j] = gf_add(sys_ret_poly[j], 
										 gf_multp(sys_elm_poly[i][j], recv_poly[i]));
			}
		}
	}

#if 1//(1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != sys_ret_poly[i])
		{
			DEBUG_NOTICE("sys_ret_poly: %ld %ld %ld | %x\n",
						 x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 sys_ret_poly[i]);
		}
	}

	memset(sys_ret_cwd, 0xFF ,sizeof(unsigned char) * CODEWORD_LEN);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
#if 0	
		if(0 == sys_proc_col_flag[i])
		{
			continue;
		}
#endif
		if(0 == sys_proc_col_flag[i])
		{
			eva_val = sys_hasse_dev_cal(sys_ret_poly, af_pnt[i]);
			DEBUG_NOTICE("sys_ret_pos_check: %ld | %x %x %x\n", i, eva_val, cwd_poly[i], recv_poly[i]);
			sys_ret_cwd[i] = eva_val;
		}
		else
		{
			sys_ret_cwd[i] = recv_poly[i];
		}
#if 0
		if((eva_val != recv_poly[i])
			&& (1 == sys_proc_col_flag[i]))
		{
			DEBUG_SYS("sys_ret_cwd_err: %ld | %x %x\n", i, eva_val, recv_poly[i]);
			val = 1;
		}
#endif
	}
	
	cwd_dim = poly_degree_cal(sys_ret_poly) + 1 - (GF_Q) * (GF_Q - 1) / 2;
	//cwd_dim = check_cwd_dimension(sys_ret_cwd);

	if((0 == val)
		&& (MESSAGE_LEN >= cwd_dim))
	{
		DEBUG_NOTICE("sys_ret_ok: %d\n", cwd_dim);
	}
	else
	{
		DEBUG_SYS("sys_ret_err: %d\n", cwd_dim);
	}
#endif

	return 0;
}

int sys_kot_intp(unsigned char poly_seq[][MAX_POLY_TERM_SIZE], unsigned char *input_point)
{
	long long j = 0, k = 0;
	unsigned char dev_val = 0xFF;
	long long min_idx = -1, min_deg = 65536, deg_tmp = 0;
	unsigned char dev_store[CFG_SYS_GEN_POLY];

	for(j = 0; j < CFG_SYS_GEN_POLY; j++)
	{
		dev_val = sys_hasse_dev_cal(poly_seq[j], input_point);
		dev_store[j] = dev_val;
		DEBUG_NOTICE("dev_val: %ld | %x %x | %x\n",
		             j,
		             input_point[0],
		             input_point[1],
		             dev_val);

		if(0xFF != dev_val)
		{
			deg_tmp = poly_degree_cal(poly_seq[j]);
			if(min_deg > deg_tmp)
			{
				min_deg = deg_tmp;
				min_idx = j;
			}
		}
	}

	DEBUG_NOTICE("min_idx: %x %x | %ld\n", input_point[0], input_point[1], min_idx);
	
	/*update poly.*/
	if(-1 != min_idx)
	{
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			memset(sys_poly_tmp[j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

			if(0xFF == dev_store[j])
			{
				continue;
			}
			if(j == min_idx)
			{
				sys_poly_update(sys_poly_tmp[j],
				                poly_seq[j],
				                poly_seq[min_idx],
				                dev_store[min_idx],
				                dev_store[j],
				                input_point[0],
				                1);
			}
			else
			{
				sys_poly_update(sys_poly_tmp[j],
				                poly_seq[j],
				                poly_seq[min_idx],
				                dev_store[min_idx],
				                dev_store[j],
				                input_point[0],
				                0);
			}
		}

		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			if(0xFF == dev_store[j])
			{
				continue;
			}

			memcpy(poly_seq[j], sys_poly_tmp[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);			
		
#if (1 == TEST_MODE)
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != poly_seq[j][k])
				{
					DEBUG_NOTICE("sys_act_poly_add: %d | %ld %ld %ld | %x\n",
								 j,
					             x_term_degree_table[k],
						         y_term_degree_table[k],
						         z_term_degree_table[k],
						         poly_seq[j][k]);
				}
			}
#endif
		}
	}

	return 0;
}

int sys_gen_test()
{
	long long i = 0, j = 0, k = 0;
	long long point_idx = 0, locator_idx = 0;
	unsigned char dev_val = 0xFF;
	long long min_idx = -1, min_deg = 65536, deg_tmp = 0;
	unsigned char dev_store[CFG_SYS_GEN_POLY];
	memset(dev_store, 0xFF, sizeof(unsigned char) * CFG_SYS_GEN_POLY);
	long long poly_deg[CFG_SYS_GEN_POLY];
	memset(poly_deg, 0, sizeof(long long) * CFG_SYS_GEN_POLY);
	int val = 0;

	memset(sys_proc_col_flag, 0, sizeof(unsigned char) * CODEWORD_LEN);
	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		memset(sys_poly_tmp[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	sys_intp_point_gen(chnl_rel_order_idx);
	sys_intp_poly_init();

	point_idx = 0;

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		if(1 == sys_proc_col_flag[sys_intp_seq[point_idx]])
		{
			continue;
		}
#if 0	
		min_idx = -1, min_deg = 65536, deg_tmp = 0;
		
		/*cal. hasse dev.*/
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			dev_val = sys_hasse_dev_cal(sys_intp_poly[j], sys_intp_point[point_idx]);
			dev_store[j] = dev_val;
			DEBUG_NOTICE("dev_val: %ld %ld | %x %x | %x\n",
			             i,
			             j,
			             sys_intp_point[point_idx][0],
			             sys_intp_point[point_idx][1],
			             dev_val);

			if(0xFF != dev_val)
			{
				deg_tmp = poly_degree_cal(sys_intp_poly[j]);
				if(min_deg > deg_tmp)
				{
					min_deg = deg_tmp;
					min_idx = j;
				}
			}
		}
		
		DEBUG_NOTICE("min_idx: %ld %ld | %ld\n", i, point_idx, min_idx);
		
		/*update poly.*/
		if(-1 != min_idx)
		{
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				memset(sys_poly_tmp[j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

				if(0xFF == dev_store[j])
				{
					continue;
				}

				if(j == min_idx)
				{
					sys_poly_update(sys_poly_tmp[j],
					                sys_intp_poly[j],
					                sys_intp_poly[min_idx],
					                dev_store[min_idx],
					                dev_store[j],
					                sys_intp_point[point_idx][0],
					                1);
				}
				else
				{
					sys_poly_update(sys_poly_tmp[j],
					                sys_intp_poly[j],
					                sys_intp_poly[min_idx],
					                dev_store[min_idx],
					                dev_store[j],
					                sys_intp_point[point_idx][0],
					                0);
				}
			}

			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				if(0xFF == dev_store[j])
				{
					continue;
				}

				memcpy(sys_intp_poly[j], sys_poly_tmp[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);			
			}
		}
#else
		sys_kot_intp(sys_intp_poly, sys_intp_point[point_idx]);
#endif

		/*set intp. point flag*/
		sys_proc_col_flag[sys_intp_seq[point_idx]] = 1;
		this_new_col_idx = sys_intp_seq[point_idx];

		/*cal. deg. for Inf. set check*/
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			poly_deg[j] = poly_degree_cal(sys_intp_poly[j]);
		}

#if 1/*check and determine the next intp. point*/
		/*equal or smaller? only smaller?*/
		if((MESSAGE_LEN - (GF_Q) * (GF_Q - 1) / 2) < i)
		{
			val = sys_check_poly_w(poly_deg);

			if(0 != val)
			{
				val = sys_check_part_dep(i + 1);
			}

			if(0 != val)
			{
				val = sys_check_full_dep(i + 1);
			}

			if(0 != val)
			{
				for(j = 0; j < CFG_SYS_GEN_POLY; j++)
				{
					memcpy(sys_intp_poly[j], sys_prev_poly[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				}
				sys_proc_col_flag[sys_intp_seq[point_idx]] = 0;
				i = i - 1;
				
				/**notice this to avoid lock loop*/
				point_idx++;
				if(K_PLUS_G <= point_idx)
				{
					DEBUG_NOTICE("K_PLUS_G out: %ld\n", point_idx);
					point_idx = 0;
					continue;
					//break;
				}

				continue;
			}
			else
			{
				for(j = 0; j < CFG_SYS_GEN_POLY; j++)
				{
					memcpy(sys_prev_poly[j], sys_intp_poly[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				}
			}
		}		
		else
		{
			for(j = 0; j < CFG_SYS_GEN_POLY; j++)
			{
				memcpy(sys_prev_poly[j], sys_intp_poly[j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
		}		
#endif

		/*avoid overflow*/
		point_idx++;
		if(K_PLUS_G <= point_idx)
		{
			DEBUG_NOTICE("K_PLUS_G out: %ld\n", point_idx);
			//break;
			point_idx = 0;
			continue;
		}

#if (1 == TEST_MODE)
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			//poly_deg[j] = poly_degree_cal(sys_intp_poly[j]);
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_intp_poly[j][k])
				{
					DEBUG_NOTICE("sys_intp_poly: %ld | %ld | %ld %ld %ld | %x\n",
								 j,
								 poly_deg[j],
					             x_term_degree_table[k],
					             y_term_degree_table[k],
					             z_term_degree_table[k],
					             sys_intp_poly[j][k]);
				}
			}
		}
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			DEBUG_NOTICE("sys_poly_deg: %ld | %ld\n", j, poly_deg[j]);
		}

		//sys_has_dev_test();
#endif
	}

	min_idx = -1;
	min_deg = 65536;
	for(i = 0; i < CFG_SYS_GEN_POLY; i++)
	{
		poly_deg[i] = poly_degree_cal(sys_intp_poly[i]);
		if(poly_deg[i] < min_deg)
		{
			min_deg = poly_deg[i];
			min_idx = i;
		}
	}
	memcpy(sys_intp_min_poly, sys_intp_poly[min_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	/*check if dev. are zero*/
#if (1 == TEST_MODE)	
	sys_has_dev_test();
#endif

	/*check if there is aj_factor*/
#if 0	
	for(i = 0; i < GF_FIELD; i++)
	{
		check_aj_factor(i);
	}
#endif	
	
	/*process the poly. according to aj_factor*/
#if 0	
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < GF_FIELD; j++)
		{
			if(power_polynomial_table[j][0] == af_pnt[i][0])
			{
				locator_idx = j;
			}
		}

		if((1 == aj_factor_flag[locator_idx][0])
			&& (1 == sys_proc_col_flag[i]))
		{
			poly_aj_factor_div(0, locator_idx);
		}
		
		if((0 == aj_factor_flag[locator_idx][0])
			&& (1 == sys_proc_col_flag[i]))
		{
			poly_aj_factor_trans(locator_idx);
		}
	}
#else
	//cnt_switch = 1;
	sys_elm_poly_gen();
	//cnt_switch = 0;
	sys_ret_poly_gen();
	//cnt_switch = 0;
#endif

	return 0;
}


/*the poly_array for br dec. should be borrowed.*/
int sys_br_g_poly_gen()
{
	long long i = 0, j = 0, k = 0;
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE];
	int val = 0;

	memset(br_g_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	br_g_poly[0] = 0x0;

	for(i = 0; i < GF_FIELD; i++)
	{
		val = 0;
		for(j = 0; j < GF_Q; j++)
		{
			if(0 == sys_proc_col_flag[i * GF_Q + j])
			{
				val = 1;
				break;
			}
		}
		if(0 == val)
		{
			continue;
		}
	
		memset(tmp_poly_x, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		
		for(j = (MAX_POLY_TERM_SIZE - 1); j >= 0; j--)
		{
			if(0xFF != br_g_poly[j])
			{
				tmp_poly_x[j] = gf_multp(br_g_poly[j], power_polynomial_table[i][0]);

				k = term_search((x_term_degree_table[j] + 1), y_term_degree_table[j], z_term_degree_table[j]);
				br_g_poly[k] = br_g_poly[j];
				br_g_poly[j] = 0xFF;
			}
		}
		
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF == tmp_poly_x[j])
			{
				continue;
			}
			else if((0xFF == br_g_poly[j])
					&& (0xFF != tmp_poly_x[j]))
			{
				br_g_poly[j] = tmp_poly_x[j];
			}
			else
			{
				br_g_poly[j] = gf_add(br_g_poly[j], tmp_poly_x[j]);
			}
		}
	}
	
	memset(br_g_gamma_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	br_g_gamma_poly[0] = 0x0;

	for(i = 0; i < GF_FIELD; i++)
	{
		val = 0;
		for(j = 0; j < GF_Q; j++)
		{
			if(0 == sys_proc_col_flag[i * GF_Q + j])
			{
				val = 1;
				break;
			}
		}
		if(1 == val)
		{
			continue;
		}
	
		memset(tmp_poly_x, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		
		for(j = (MAX_POLY_TERM_SIZE - 1); j >= 0; j--)
		{
			if(0xFF != br_g_gamma_poly[j])
			{
				tmp_poly_x[j] = gf_multp(br_g_gamma_poly[j], power_polynomial_table[i][0]);

				k = term_search((x_term_degree_table[j] + 1), y_term_degree_table[j], z_term_degree_table[j]);
				br_g_gamma_poly[k] = br_g_gamma_poly[j];
				br_g_gamma_poly[j] = 0xFF;
			}
		}
		
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF == tmp_poly_x[j])
			{
				continue;
			}
			else if((0xFF == br_g_gamma_poly[j])
					&& (0xFF != tmp_poly_x[j]))
			{
				br_g_gamma_poly[j] = tmp_poly_x[j];
			}
			else
			{
				br_g_gamma_poly[j] = gf_add(br_g_gamma_poly[j], tmp_poly_x[j]);
			}
		}
	}

	//her_convert(br_g_poly);

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != br_g_poly[i])
		{
			DEBUG_NOTICE("sys_br_g_poly: %ld %ld %ld | %x\n",
						  x_term_degree_table[i],
						  y_term_degree_table[i],
						  z_term_degree_table[i],
						  br_g_poly[i]);
		}
	}
	
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != br_g_gamma_poly[i])
		{
			DEBUG_NOTICE("br_g_gamma_poly: %ld %ld %ld | %x\n",
						  x_term_degree_table[i],
						  y_term_degree_table[i],
						  z_term_degree_table[i],
						  br_g_gamma_poly[i]);
		}
	}

	return 0;
}

int sys_br_lag_part_gen()
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char x_poly[MAX_POLY_TERM_SIZE], y_poly[MAX_POLY_TERM_SIZE];
	unsigned char tmp_poly_x[MAX_POLY_TERM_SIZE], tmp_poly_y[MAX_POLY_TERM_SIZE];
	unsigned char y_her_curve_find_flag = 0;
	unsigned char x_div = 0x0, y_div = 0x0;
	unsigned char check_val = 0xFF;

	unsigned char ret_x_flag = 0;
	long long aj_proc_cnt = 0;

	for(i = 0; i < CODEWORD_LEN; i++)//af. point
	{	
		if(1 == sys_proc_col_flag[i])
		{
			DEBUG_NOTICE("sys_gamma skip: %d\n", i);
			continue;
		}
		
		memset(br_lag_part_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		
		memset(x_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		x_poly[0] = 0x0;
		memset(y_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		y_poly[0] = 0x0;
		x_div = 0x0;
		y_div = 0x0;

		for(j = 0; j < GF_FIELD; j++)//locator
		{
			ret_x_flag = 0;
			aj_proc_cnt = 0;

			for(k = 0; k < CODEWORD_LEN; k++)
			{
				if((af_pnt[k][0] == power_polynomial_table[j][0])
					&& (0 == sys_proc_col_flag[k]))
				{
					ret_x_flag = 1;
					break;
				}
				
				if((1 == sys_proc_col_flag[k])
					&& (af_pnt[k][0] == power_polynomial_table[j][0]))
				{
					DEBUG_NOTICE("sys_aj_gamma_skip: %x\n", power_polynomial_table[j][0]);
					aj_proc_cnt++;
				}
			}
			/*notice this!*/
			if(GF_Q == aj_proc_cnt)
			{
				ret_x_flag = 0;
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
								br_lag_part_poly[i][l] = gf_add(br_lag_part_poly[i][l],
														gf_multp(x_poly[k], y_poly[j]));
								br_lag_part_poly[i][l] = gf_div(br_lag_part_poly[i][l], x_div);
								br_lag_part_poly[i][l] = gf_div(br_lag_part_poly[i][l], y_div);
							}
						}
					}
				}
			}
		}

		//her_convert(br_lag_poly[i]);

		/*if there is no common factor for gamma, the part lag. poly. shpuld be 1, not 0.*/
#if (1 == TEST_MODE)
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != br_lag_part_poly[i][j])
			{
				DEBUG_NOTICE("br_lag_part_poly: %ld | %ld %ld %ld | %x\n",
							  i,
							  x_term_degree_table[j],
							  y_term_degree_table[j],
							  z_term_degree_table[j],
							  br_lag_part_poly[i][j]);
			}
		}
#endif		
	}

	return 0;
}

int sts_br_vct_trans(long long tv_idx)
{
	long long i = 0, j = 0;
	unsigned char g_val = 0xFF;

	memset(sys_tv_trans, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(1 == sys_proc_col_flag[i])
		{
			continue;
		}
		
		g_val = poly_eva_x_y(br_g_gamma_poly, af_pnt[i][0], af_pnt[i][1]);
		sys_tv_trans[i] = gf_add(tst_vct[tv_idx][i], sys_ret_cwd[i]);
		sys_tv_trans[i] = gf_div(sys_tv_trans[i], g_val);
		DEBUG_NOTICE("sys_tv_trans: %ld | %x\n", i, sys_tv_trans[i]);
	}

	return 0;
}

int sys_br_k_poly_construct(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	long long intp_idx = 0;

	//for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			intp_idx = chnl_rel_order_idx[CODEWORD_LEN - 1 - j];

			if(0xFF == sys_tv_trans[intp_idx])
			{
				continue;
			}
#if 1
			if(1 == sys_proc_col_flag[intp_idx])
			{
				continue;
			}
#endif

			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_lag_part_poly[intp_idx][k])
				{
					if((0xFF == sys_tv_trans[intp_idx])
						|| (0xFF == br_lag_part_poly[intp_idx][k]))
					{
						tmp_poly[k] = 0xFF;
					}
					else if(0x0 == sys_tv_trans[intp_idx])
					{
						tmp_poly[k] = br_lag_part_poly[intp_idx][k];
					}
					else if(0x0 == br_lag_part_poly[intp_idx][k])
					{
						tmp_poly[k] = sys_tv_trans[intp_idx];
					}
					else
					{
						tmp_poly[k] = gf_multp(sys_tv_trans[intp_idx], br_lag_part_poly[intp_idx][k]);
					}

					if((0xFF == br_k_poly[tv_idx][k])
						&& (0xFF == tmp_poly[k]))
					{
						br_k_poly[tv_idx][k] = 0xFF;
					}
					else if(0xFF == br_k_poly[tv_idx][k])
					{
						br_k_poly[tv_idx][k] = tmp_poly[k];
					}
					else if(0xFF == tmp_poly[k])
					{
						br_k_poly[tv_idx][k] = br_k_poly[tv_idx][k];
					}
					else
					{
						br_k_poly[tv_idx][k] = gf_add(br_k_poly[tv_idx][k], tmp_poly[k]);
					}

#if 1
					DEBUG_NOTICE("sys_br_k_poly_cal: %ld %ld | %x %x | %ld %ld %ld | %x\n",
								  tv_idx,
								  intp_idx,
								  sys_tv_trans[intp_idx],
								  br_lag_part_poly[intp_idx][k],
					              x_term_degree_table[k],
					              y_term_degree_table[k],
					              z_term_degree_table[k],
					              br_k_poly[tv_idx][k]);
#endif
				}
			}

			if((CODEWORD_LEN - 1 - ETA) == j)
			{
				memcpy(br_k_com_poly, br_k_poly[tv_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}
		}

#if (1 == TEST_MODE)
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != br_k_poly[tv_idx][k])
			{
				DEBUG_NOTICE("sys_br_k_poly: %d | %ld %ld %ld | %x\n",
				              tv_idx,
				              x_term_degree_table[k],
				              y_term_degree_table[k],
				              z_term_degree_table[k],
				              br_k_poly[tv_idx][k]);
			}
		}
#endif

		memset(br_zk_poly[tv_idx], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		memcpy(br_zk_poly[tv_idx], br_k_poly[tv_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		l = term_search(0, 0, 1);
		br_zk_poly[tv_idx][l] = 0x0;

#if (1 == TEST_MODE)
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != br_zk_poly[tv_idx][k])
			{
				DEBUG_NOTICE("sys_br_zk_poly: %d | %ld %ld %ld | %x\n",
				              tv_idx,
				              x_term_degree_table[k],
				              y_term_degree_table[k],
				              z_term_degree_table[k],
				              br_zk_poly[tv_idx][k]);
			}
		}
#endif		
	}

	return 0;
}

int sys_br_ret_deg_down()
{
	long long w_z = GF_Q * pole_basis_pow[MESSAGE_LEN - 1][0]
				  + (GF_Q + 1) * pole_basis_pow[MESSAGE_LEN - 1][1];

	long long g_gamma_deg = poly_degree_cal(br_g_gamma_poly);
	
	sys_ret_z_deg = w_z - g_gamma_deg;
	DEBUG_NOTICE("sys_ret_z_deg: %ld - %ld = %ld\n", w_z, g_gamma_deg, sys_ret_z_deg);

	return 0;
}

int sys_br_m_poly_construct(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;

	//for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < (BR_BASIS_NUM / 2); j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_g_poly[k])
				{
					l = term_search(x_term_degree_table[k], (y_term_degree_table[k] + j), z_term_degree_table[k]);
					br_m_poly[tv_idx][j][l] = br_g_poly[k];
				}
			}
		}
		
		for(j = (BR_BASIS_NUM / 2); j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_zk_poly[tv_idx][k])
				{
					l = term_search(x_term_degree_table[k], (y_term_degree_table[k] + j - (BR_BASIS_NUM / 2)), z_term_degree_table[k]);
					br_m_poly[tv_idx][j][l] = br_zk_poly[tv_idx][k];
				}
			}
		}
	}

	//for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			br_y_high_deg_trans(br_m_poly[tv_idx][j]);
		}		
	}

#if (1 == TEST_MODE)
	//for(i = 0; i < tst_vct_num; i++)
	{
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_m_poly[tv_idx][j][k])
				{
					DEBUG_NOTICE("sys_br_m_poly: %d %d | %ld %ld %ld | %x\n",
					              tv_idx,
					              j,
								  x_term_degree_table[k],
								  y_term_degree_table[k],
								  z_term_degree_table[k],
								  br_m_poly[tv_idx][j][k]);
				}
			}
		}		
	}
#endif	

	return 0;
}

int sys_br_v_matric_gen(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	long long w_z = sys_ret_z_deg;
	long long up_x_deg = 0;

	for(i = 0; i < BR_BASIS_NUM; i++)//for row
	{
		for(j = 0; j < BR_BASIS_NUM; j++)//for col
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if((0xFF != br_m_poly[tv_idx][i][k])
					&& ((j % GF_Q) == y_term_degree_table[k])
					&& ((j / GF_Q) == z_term_degree_table[k]))
				{
					l = term_search(x_term_degree_table[k], 0, 0);
					br_v_matrix_poly[tv_idx][i][j][l] = br_m_poly[tv_idx][i][k];
				}
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < BR_BASIS_NUM; i++)//for row
	{
		for(j = 0; j < BR_BASIS_NUM; j++)//for col
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					DEBUG_NOTICE("sys_br_v_matrix_poly: %ld | %ld %ld | %ld %ld %ld | %x\n",
					             tv_idx,
					             i,
					             j,
					             x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 br_v_matrix_poly[tv_idx][i][j][k]);
				}
			}
		}
	}
#endif

	for(i = 0; i < BR_BASIS_NUM; i++)//for row
	{
		for(j = 0; j < BR_BASIS_NUM; j++)//for col
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

			up_x_deg = (((GF_Q + 1) * (j % GF_Q)) + (j / GF_Q * w_z)) / GF_Q;

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					l = term_search((x_term_degree_table[k] + up_x_deg), 0, 0);
					tmp_poly[l] = br_v_matrix_poly[tv_idx][i][j][k];
					br_v_matrix_poly[tv_idx][i][j][k] = 0xFF;
				}
			}
			memcpy(br_v_matrix_poly[tv_idx][i][j], tmp_poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
#if (1 == TEST_MODE)
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					DEBUG_NOTICE("sys_br_v_matrix_poly_up: %ld | %ld %ld | %ld | %ld %ld %ld | %x\n",
					             tv_idx,
					             i,
					             j,
					             up_x_deg,
					             x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 br_v_matrix_poly[tv_idx][i][j][k]);
				}
			}
#endif
		}
	}

	return 0;
}

int sys_br_g_poly_ret(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	long long w_z = sys_ret_z_deg;
	long long down_x_deg = 0;

	for(i = 0; i < BR_BASIS_NUM; i++)//for row
	{
		for(j = 0; j < BR_BASIS_NUM; j++)//for col
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

			down_x_deg = (((GF_Q + 1) * (j % GF_Q)) + (j / GF_Q * w_z)) / GF_Q;

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					if(0 > (x_term_degree_table[k] - down_x_deg))
					{
						DEBUG_SYS("down_x_deg_err: %ld %ld | %ld %ld\n", i, j, x_term_degree_table[k], down_x_deg);
					}
				
					l = term_search((x_term_degree_table[k] - down_x_deg), 0, 0);
					tmp_poly[l] = br_v_matrix_poly[tv_idx][i][j][k];
					br_v_matrix_poly[tv_idx][i][j][k] = 0xFF;
				}
			}
			memcpy(br_v_matrix_poly[tv_idx][i][j], tmp_poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
#if (1 == TEST_MODE)
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					DEBUG_NOTICE("br_v_matrix_poly_down: %ld | %ld %ld | %ld | %ld %ld %ld | %x\n",
					             tv_idx,
					             i,
					             j,
					             down_x_deg,
					             x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 br_v_matrix_poly[tv_idx][i][j][k]);
				}
			}
#endif
		}
	}

	return 0;
}

long long sys_poly_degree_cal(unsigned char *poly)
{
	long long i = 0;
	/*notice these degrees should be started form -1*/
	long long max_degree = -100, degree_val = -100;

	long long w_z = sys_ret_z_deg;

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

	DEBUG_NOTICE("max_degree: %ld %ld\n", w_z, max_degree);

	return max_degree;
}

int sys_br_q_poly_ret(long long tv_idx)
{
	long long i = 0, j = 0, k = 0, l = 0;
	long long min_q_deg = 65536, tmp_deg = 0;
	long long min_idx = -1;

	for(i = 0; i < BR_BASIS_NUM; i++)
	{
		memset(br_q_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		
		for(j = 0; j < BR_BASIS_NUM; j++)
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != br_v_matrix_poly[tv_idx][i][j][k])
				{
					l = term_search(x_term_degree_table[k],
					                y_term_degree_table[k] + (j % GF_Q),
					                z_term_degree_table[k] + (j / GF_Q));
					br_q_poly[i][l] = br_v_matrix_poly[tv_idx][i][j][k];
				}
			}
		}
		
		tmp_deg = sys_poly_degree_cal(br_q_poly[i]);
		if(tmp_deg <= min_q_deg)
		{
			min_q_deg = tmp_deg;
			min_idx = i;
		}
		
#if (1 == TEST_MODE)
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != br_q_poly[i][k])
			{
				DEBUG_NOTICE("br_q_poly: %ld | %ld %ld %ld | %x\n",
				             i,
				             x_term_degree_table[k],
							 y_term_degree_table[k],
							 z_term_degree_table[k],
							 br_q_poly[i][k]);
			}
		}
#endif		
	}

	memcpy(min_intp_poly, br_q_poly[min_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	min_intp_idx = min_idx;
	DEBUG_NOTICE("br_min_idx: %ld\n", min_intp_idx);

	return 0;
}

int sys_her_fac(unsigned char *poly,
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

#if 0//(1 == CFG_RET)	
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

int sys_check_ret_cwd()
{
	long long i = 0;
	int val = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(sys_ret_cwd[i] != cwd_poly[i])
		{
			val = 1;
			break;
		}
	}
	if(0 == val)
	{
		memcpy(est_cwd_poly, sys_ret_cwd, sizeof(unsigned char) * CODEWORD_LEN);
		cwd2msg(sys_ret_cwd, est_msg_poly);
	}

	return val;
}

int sys_br_test()
{
	long long i = 0, j = 0, k = 0;
	int popov_flag = 1;
	long long ms_to_cnt = 0;
	unsigned char test_poly[MAX_POLY_TERM_SIZE];
	memset(test_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	int val = 0;
	long long genus = GF_Q * (GF_Q - 1) / 2;
	/*for GS (m = 1), that is the GS decoding bound*/
	long long radius = (CODEWORD_LEN - (MESSAGE_LEN - genus + 1) - 2 * genus) / 2;
	if(0 == ((MESSAGE_LEN - genus + 1) % 2))
	{
		radius = radius - 1;
	}
	DEBUG_NOTICE("GS radius: %ld\n", radius);

	sys_br_g_poly_gen();
	sys_br_lag_part_gen();
	sys_br_ret_deg_down();
	memcpy(v_poly, br_g_gamma_poly, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		v_val[i] = poly_eva_x_y(v_poly, af_pnt[i][0], af_pnt[i][1]);
		DEBUG_NOTICE("sys_v_val: %ld | %x %x | %x\n",
					 i,
					 af_pnt[i][0],
					 af_pnt[i][1],
		             v_val[i]);
	}

	for(k = 0; k < tst_vct_num; k++)
	{
		sts_br_vct_trans(k);
		sys_br_k_poly_construct(k);
		sys_br_m_poly_construct(k);
		sys_br_v_matric_gen(k);

		ms_to_cnt = 0;/*notice this, this should be clear*/
		popov_flag = is_popov_form(k);
		while(0 == popov_flag)
		{
			ms_reduction(k);
			popov_flag = is_popov_form(k);
			ms_to_cnt++;
			if(CODEWORD_LEN < ms_to_cnt)
			{
				DEBUG_SYS("br_ms_to: %ld\n", ms_to_cnt);
				break;
			}
		}

		sys_br_g_poly_ret(k);
		sys_br_q_poly_ret(k);

		sys_her_fac(min_intp_poly, tv_est_msg[k], tv_est_cwd[k]);
		tv_dec_output_flag[k] = 1;

		for(i = 0; i < CODEWORD_LEN; i++)
		{
			tv_est_cwd[k][i] = gf_add(tv_est_cwd[k][i], sys_ret_cwd[i]);
			DEBUG_NOTICE("br_sys_tv_est_cwd: %d | %ld | %x %x %x\n", tv_err[k], i, tv_est_cwd[k][i], cwd_poly[i], recv_poly[i]);
			if((tv_est_cwd[k][i] != cwd_poly[i])
				&& (tv_err[k] <= radius))
			{
				DEBUG_SYS("sys_dec_err: %ld | %ld | %ld | %x %x | %x %x\n",
						  k,
				          tv_err[k],
				          i,
				          gf_add(sys_ret_cwd[i], tst_vct[k][i]),
				          gf_add(sys_ret_cwd[i], cwd_poly[i]),
				          tv_est_cwd[k][i],
				          cwd_poly[i]);
			}
		}
	}
	br_poly_clear();

	her_lcc_check_result();

	return 0;
}

int sys_poly_q0q1_get(unsigned char *poly_seq)
{
	long long i = 0, j = 0;
	
	memset(q0_poly_coef, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	memset(q1_poly_coef, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != poly_seq[i])
		{
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
				q1_poly_coef[j] = poly_seq[i];
			}
			else
			{
				q0_poly_coef[i] = poly_seq[i];
			}
		}
	}
#if (1 == TEST_MODE)	
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != q0_poly_coef[i])
		{
			DEBUG_NOTICE("q0_poly: %ld %ld %ld | %x\n",
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
			DEBUG_NOTICE("q1_poly: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
				         y_term_degree_table[i],
				         z_term_degree_table[i],
				         q1_poly_coef[i]);
		}
	}
#endif
	return 0;
}

int sys_g_poly_all_val_cal()
{
	long long i = 0, j = 0, k = 0, l = 0;
	long long y_deg = 0;
	unsigned char x_val = 0xFF, y_val = 0xFF, tmp_val = 0xFF;

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			memset(sys_g_store_val[i][j], 0, sizeof(unsigned char) * CODEWORD_LEN);
		}
	}

	for(l = 0; l < CODEWORD_LEN; l++)
	{
		for(i = 0; i < KOT_INTP_POLY_NUM; i++)
		{
			for(j = 0; j < KOT_INTP_POLY_NUM; j++)
			{
				y_deg = j % GF_Q;

				x_val = 0xFF;
				y_val = 0x0;//notice this
				for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
				{
					if((0 != y_term_degree_table[k])
						|| (0 != z_term_degree_table[k]))
					{
						continue;
					}
					
					if(0xFF != sys_break_mat_poly[i][j][k])
					{
						tmp_val = gf_pow_cal(af_pnt[l][0], x_term_degree_table[k]);
						tmp_val = gf_multp(tmp_val, sys_break_mat_poly[i][j][k]);
						x_val = gf_add(x_val, tmp_val);

						y_val = gf_pow_cal(af_pnt[l][1], y_deg);
					}
				}
				
				sys_g_store_val[i][j][l] = gf_multp(x_val, y_val);
				if(0xFF != sys_g_store_val[i][j][l])
				{
					DEBUG_NOTICE("sys_g_store_val: %ld %ld | %x %x | %x %x | %x\n",
								 i,
								 j,
								 af_pnt[l][0],
								 af_pnt[l][1],
								 x_val,
								 y_val,
								 sys_g_store_val[i][j][l]);
				}
			}
		}
	}

	return 0;
}

unsigned char sys_hasse_full_dev_cal(long long poly_coef_idx, long long pnt_idx, unsigned char test_sym)
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char tmp_val = 0xFF, ele_val = 0xFF, val = 0xFF;
	long long z_deg = 0;
	unsigned char q_val[KOT_INTP_POLY_NUM];
	memset(q_val, 0xFF, sizeof(unsigned char) * KOT_INTP_POLY_NUM);
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != sys_intp_coef_mat[poly_coef_idx][i][k])
			{
				tmp_val = gf_pow_cal(af_pnt[pnt_idx][0], x_term_degree_table[k]);
				tmp_val = gf_multp(tmp_val, sys_intp_coef_mat[poly_coef_idx][i][k]);
				q_val[i] = gf_add(q_val[i], tmp_val);
			}
		}
	}

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		z_deg = i / GF_Q;
		ele_val = 0xFF;

		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			tmp_val = gf_multp(q_val[j], sys_g_store_val[j][i][pnt_idx]);
			ele_val = gf_add(tmp_val, ele_val);
		}
		
		tmp_val = gf_pow_cal(test_sym, z_deg);
		ele_val = gf_multp(tmp_val, ele_val);
		
		val = gf_add(val, ele_val);
	}

	return val;
}

int sys_kot_poly_init()
{
	long long i = 0, j = 0, k = 0;

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			memset(sys_intp_coef_mat[i][j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
		k = term_search(0, 0, 0);
		sys_intp_coef_mat[i][i][k] = 0x0;
		//memcpy(sys_intp_coef_mat_tmp[i][j], sys_intp_coef_mat[i][j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}

	return 0;
}

int sys_br_v_matric_break()
{
	long long i = 0, j = 0, k = 0, l = 0;
	long long deg_tmp = 0, max_deg = 0;
	long long y_deg = 0, z_deg = 0;
	long long w_z = GF_Q * pole_basis_pow[MESSAGE_LEN - 1][0]
				  + (GF_Q + 1) * pole_basis_pow[MESSAGE_LEN - 1][1];

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			memset(sys_break_mat_poly[i][j], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
		}
		memset(sys_g_store_deg[i], -65536, sizeof(long long) * KOT_INTP_POLY_NUM);
	}

	for(i = 0; i < (KOT_INTP_POLY_NUM / 2); i++)//for row
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)//for col
		{
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if((0xFF != sys_intp_poly[i][k])
					&& ((j % GF_Q) == y_term_degree_table[k])
					&& ((j / GF_Q) == z_term_degree_table[k]))
				{
					l = term_search(x_term_degree_table[k], 0, 0);
					sys_break_mat_poly[i][j][l] = gf_add(sys_intp_poly[i][k], sys_break_mat_poly[i][j][l]);
				}
			}
		}
	}
	
	for(i = (KOT_INTP_POLY_NUM / 2); i < KOT_INTP_POLY_NUM; i++)
	{
		l = term_search(0, 0, 0);
		sys_break_mat_poly[i][i][l] = 0x0;
	}

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)//for row
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)//for col
		{
			y_deg = j % GF_Q;
			z_deg = j / GF_Q;
			max_deg = -65536;

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_break_mat_poly[i][j][k])
				{
					deg_tmp = x_term_degree_table[k] * GF_Q
							+ y_deg * (GF_Q + 1)
							+ z_deg * w_z;
					
					if(deg_tmp > max_deg)
					{
						max_deg = deg_tmp;
					}

					DEBUG_NOTICE("sys_break_mat_poly: %ld %ld | %ld %ld %ld | %x\n",
					             i,
					             j,
					             x_term_degree_table[k],
								 y_term_degree_table[k],
								 z_term_degree_table[k],
								 sys_break_mat_poly[i][j][k]);
				}
			}
			
			sys_g_store_deg[i][j] = max_deg;
			DEBUG_NOTICE("sys_g_store_deg: %ld %ld | %ld\n", i, j, sys_g_store_deg[i][j]);
		}
	}

	return 0;
}

int sys_poly_mat_mul(long long sel_idx)
{
	long long i = 0, j = 0, k = 0, l = 0, s = 0;
	unsigned char tmp_poly[MAX_POLY_TERM_SIZE];
	memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		memset(sys_mul_poly[i], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
	}
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			memset(tmp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_intp_coef_mat[sel_idx][i][k])
				{
					for(l = 0; l < MAX_POLY_TERM_SIZE; l++)
					{
						if(0xFF != sys_break_mat_poly[i][j][l])
						{
							s = term_search(x_term_degree_table[k] + x_term_degree_table[l], 0, 0);
							sys_mul_poly[j][s] = gf_add(sys_mul_poly[j][s],
														gf_multp(sys_intp_coef_mat[sel_idx][i][k], sys_break_mat_poly[i][j][l]));
						}
					}
				}
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
		{
			if(0xFF != sys_mul_poly[i][k])
			{
				DEBUG_NOTICE("sys_mul_poly: %ld | %ld %ld %ld | %x\n",
				             i,
				             x_term_degree_table[k],
							 y_term_degree_table[k],
							 z_term_degree_table[k],
							 sys_mul_poly[i][k]);
			}
		}
	}
#endif

	return 0;
}

int sys_mul_poly_recover()
{
	long long i = 0, j = 0, k = 0;
	long long y_deg = 0, z_deg = 0;

	memset(min_intp_poly, 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		y_deg = i % GF_Q;
		z_deg = i / GF_Q;
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != sys_mul_poly[i][j])
			{
				k = term_search(x_term_degree_table[j], y_deg, z_deg);
				min_intp_poly[k] = gf_add(min_intp_poly[k], sys_mul_poly[i][j]);
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < MAX_POLY_TERM_SIZE; i++)
	{
		if(0xFF != min_intp_poly[i])
		{
			DEBUG_NOTICE("sys_min_intp_poly: %ld %ld %ld | %x\n",
			             x_term_degree_table[i],
						 y_term_degree_table[i],
						 z_term_degree_table[i],
						 min_intp_poly[i]);
		}
	}
#endif

	return 0;
}

long long sys_min_mul_poly_idx_determine()
{
	long long i = 0, j = 0, k = 0;
	long long min_idx = -1, min_deg = 65536, deg_tmp = 0, deg_this_ele = 0;
	long long deg_mul[KOT_INTP_POLY_NUM][KOT_INTP_POLY_NUM];
	long long deg_mat[KOT_INTP_POLY_NUM][KOT_INTP_POLY_NUM];
	long long deg_this_row[KOT_INTP_POLY_NUM];
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			deg_mul[i][j] = -65536;
			deg_mat[i][j] = -65536;
		}
	}
	memset(deg_this_row, 0, sizeof(long long) * KOT_INTP_POLY_NUM);
	long long y_deg = 0, z_deg = 0;
	long long w_z = GF_Q * pole_basis_pow[MESSAGE_LEN - 1][0]
				  + (GF_Q + 1) * pole_basis_pow[MESSAGE_LEN - 1][1];
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			y_deg = j % GF_Q;
			z_deg = j / GF_Q;
			
			deg_this_ele = -65536;
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if((0xFF != intp_poly_coef[i][k])
					&& (y_deg == y_term_degree_table[k])
					&& (z_deg == z_term_degree_table[k]))
				{
					deg_tmp = x_term_degree_table[k] * GF_Q
							+ y_deg * (GF_Q + 1)
							+ z_deg * w_z;
					
					if(deg_tmp > deg_this_ele)
					{
						deg_this_ele = deg_tmp;
					}
				}
			}
			
			deg_mul[i][j] = deg_this_ele;
			DEBUG_NOTICE("deg_mul: %ld %ld | %ld\n", i, j, deg_mul[i][j]);
		}
	}
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			deg_this_ele = -65536;
			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_break_mat_poly[i][j][k])
				{
					deg_tmp = x_term_degree_table[k] * GF_Q;
					
					if(deg_tmp > deg_this_ele)
					{
						deg_this_ele = deg_tmp;
					}
				}
			}
			
			deg_mat[i][j] = deg_this_ele;
			DEBUG_NOTICE("deg_mat: %ld %ld | %ld\n", i, j, deg_mat[i][j]);
		}
	}
	
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)//row in mul
	{
		deg_this_ele = 0;
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)//col in mat
		{
			for(k = 0; k < KOT_INTP_POLY_NUM; k++)//row in mat and col in mul
			{
				deg_tmp = deg_mul[i][k] + deg_mat[k][j];
				
				if(deg_tmp > deg_this_ele)
				{
					deg_this_ele = deg_tmp;
				}
			}
		}
		deg_this_row[i] = deg_this_ele;
		
		if(deg_this_row[i] <= min_deg)
		{
			min_deg = deg_this_row[i];
			min_idx = i;
		}
	}

	DEBUG_NOTICE("min_idx: %ld | %ld\n", min_deg, min_idx);

	return min_idx;
}

int sys_final_has_dev_test()
{
	long long i = 0, j = 0;
	unsigned dev_val = 0xFF;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < CFG_SYS_GEN_POLY; j++)
		{
			dev_val = sys_hasse_full_dev_cal(j, i, gf_add(tst_vct[0][i], sys_ret_cwd[i]));
			DEBUG_NOTICE("sys_final_has_dev_test: %ld | %x %x | %x %x | %x | %x | %x\n",
						 j,
						 af_pnt[i][0],
						 af_pnt[i][1],
						 tst_vct[0][i],
						 sys_ret_cwd[i],
						 recv_poly[i],
						 gf_add(tst_vct[0][i], sys_ret_cwd[i]),
						 dev_val);
		}
	}

	return 0;
}

long long sys_intp_coef_mat_deg_cal(long long poly_coef_idx)
{
	long long i = 0, j = 0, k = 0;
	long long y_deg = 0, z_deg = 0;
	long long w_z = GF_Q * pole_basis_pow[MESSAGE_LEN - 1][0]
				  + (GF_Q + 1) * pole_basis_pow[MESSAGE_LEN - 1][1];
	long long max_deg = -65536, deg_this_col = -65536, deg_this_ele = -65536, deg_tmp = 0;
	long long val = 0;

	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		y_deg = i % GF_Q;
		z_deg = i / GF_Q;
		deg_this_col = -65536;
		
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			deg_this_ele = -65536;

			for(k = 0; k < MAX_POLY_TERM_SIZE; k++)
			{
				if(0xFF != sys_intp_coef_mat[poly_coef_idx][j][k])
				{
					deg_tmp = x_term_degree_table[k] * GF_Q;
					if(deg_tmp > deg_this_ele)
					{
						deg_this_ele = deg_tmp;
					}
				}
			}
			DEBUG_NOTICE("deg_this_ele: %ld %ld | %ld\n", i, j, deg_this_ele);
			deg_tmp = deg_this_ele + sys_g_store_deg[i][j];
			if(deg_tmp > deg_this_col)
			{
				deg_this_col = deg_tmp;
			}
		}
		DEBUG_NOTICE("deg_this_col: %ld | %ld\n", i, deg_this_col);
		//deg_tmp = deg_this_col + y_deg * (GF_Q + 1) + z_deg * w_z;
		deg_tmp = deg_this_col;
		if(deg_tmp > max_deg)
		{
			max_deg = deg_tmp;
		}
	}
	
	DEBUG_NOTICE("max_deg: %ld\n", max_deg);
	val = max_deg;

	return val;
}

int sys_kot_update(unsigned char poly_mat_update[][MAX_POLY_TERM_SIZE],
					   unsigned char poly_mat_base[][MAX_POLY_TERM_SIZE],
					   unsigned char poly_mat_min[][MAX_POLY_TERM_SIZE],
					   unsigned char dev_self,
					   unsigned char dev_min,
					   unsigned char alpha,
					   unsigned char min_flag)
{
	long long i = 0, j = 0, k = 0;

	if(1 == min_flag)
	{
		for(i = 0; i < KOT_INTP_POLY_NUM; i++)
		{
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(0xFF != poly_mat_base[i][j])
				{
					k = term_search(x_term_degree_table[j] + 1, 0, 0);
					poly_mat_update[i][k] = gf_add(poly_mat_update[i][k], poly_mat_base[i][j]);
					
					poly_mat_update[i][j] = gf_add(poly_mat_update[i][j],
										 	gf_multp(poly_mat_base[i][j], alpha));
				}
			}
		}
	}
	else
	{
		for(i = 0; i < KOT_INTP_POLY_NUM; i++)
		{
			for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
			{
				if(0xFF != poly_mat_base[i][j])
				{
					poly_mat_update[i][j] = gf_add(poly_mat_update[i][j],
										 	gf_multp(poly_mat_base[i][j], dev_min));
				}
				if(0xFF != poly_mat_min[i][j])
				{
					poly_mat_update[i][j] = gf_add(poly_mat_update[i][j],
										 	gf_multp(poly_mat_min[i][j], dev_self));
				}
			}
		}
	}

#if (1 == TEST_MODE)
	for(i = 0; i < KOT_INTP_POLY_NUM; i++)
	{
		for(j = 0; j < MAX_POLY_TERM_SIZE; j++)
		{
			if(0xFF != poly_mat_update[i][j])
			{
				DEBUG_NOTICE("poly_mat_update: %ld | %ld %ld %ld | %x\n",
							 i,
				             x_term_degree_table[j],
							 y_term_degree_table[j],
							 z_term_degree_table[j],
							 poly_mat_update[i][j]);
			}
		}
	}
#endif

	return 0;
}

int sys_kot_test()
{
	long long i = 0, j = 0, k = 0, l = 0;
	unsigned char dev_val = 0xFF;
	long long min_idx = -1, min_deg = 65536, deg_tmp = 0;
	unsigned char dev_store[KOT_INTP_POLY_NUM];
	memset(dev_store, 0xFF, sizeof(unsigned char) * KOT_INTP_POLY_NUM);

	sys_kot_poly_init();
	sys_br_v_matric_break();
	sys_g_poly_all_val_cal();

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(1 == sys_proc_col_flag[i])
		{
			continue;
		}
		
		min_deg = 65536;
		
		for(j = 0; j < KOT_INTP_POLY_NUM; j++)
		{
			dev_val = sys_hasse_full_dev_cal(j, i, gf_add(tst_vct[0][i], sys_ret_cwd[i]));
			dev_store[j] = dev_val;
			DEBUG_NOTICE("dev_val: %ld | %x %x | %x %x | %x | %x %x | %x\n",
						 j,
						 af_pnt[i][0],
						 af_pnt[i][1],
						 tst_vct[0][i],
						 sys_ret_cwd[i],
						 recv_poly[i],
						 gf_add(tst_vct[0][i], sys_ret_cwd[i]),
						 gf_add(recv_poly[i], sys_ret_cwd[i]),
						 dev_val);

			if(0xFF != dev_val)
			{
				deg_tmp = sys_intp_coef_mat_deg_cal(j);
				if(min_deg > deg_tmp)
				{
					min_deg = deg_tmp;
					min_idx = j;
				}
			}
		}

		DEBUG_NOTICE("min_idx: %ld | %x %x | %ld | %ld\n", i, af_pnt[i][0], af_pnt[i][1], min_idx, min_deg);

		if(-1 != min_idx)
		{
			for(j = 0; j < KOT_INTP_POLY_NUM; j++)
			{
				memcpy(sys_intp_coef_mat_base_tmp[j], sys_intp_coef_mat[min_idx][j], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
			}

			for(j = 0; j < KOT_INTP_POLY_NUM; j++)
			{
#if 1
				for(k = 0; k < KOT_INTP_POLY_NUM; k++)
				{
					memset(sys_intp_coef_mat_tmp[k], 0xFF, sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				}
#endif

				if(0xFF == dev_store[j])
				{
					continue;
				}
#if 1
				DEBUG_NOTICE("sys_kot_update: %ld %ld\n", j, min_idx);

				if(min_idx == j)
				{
					sys_kot_update(sys_intp_coef_mat_tmp,
								   sys_intp_coef_mat[j],
								   sys_intp_coef_mat_base_tmp,
								   dev_store[j],
								   dev_store[min_idx],
								   af_pnt[i][0],
								   1);
				}
				else
				{
					sys_kot_update(sys_intp_coef_mat_tmp,
								   sys_intp_coef_mat[j],
								   sys_intp_coef_mat_base_tmp,
								   dev_store[j],
								   dev_store[min_idx],
								   af_pnt[i][0],
								   0);
				}
#endif
				for(k = 0; k < KOT_INTP_POLY_NUM; k++)
				{
					memcpy(sys_intp_coef_mat[j][k], sys_intp_coef_mat_tmp[k], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				}
			}

			for(j = 0; j < KOT_INTP_POLY_NUM; j++)
			{			
				if(0xFF == dev_store[j])
				{
					continue;
				}			
#if 0
				for(k = 0; k < KOT_INTP_POLY_NUM; k++)
				{
					memcpy(sys_intp_coef_mat[j][k], sys_intp_coef_mat_tmp[j][k], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);
				}
#endif

#if (1 == TEST_MODE)
				for(k = 0; k < KOT_INTP_POLY_NUM; k++)
				{
					for(l = 0; l < MAX_POLY_TERM_SIZE; l++)
					{
						if(0xFF != sys_intp_coef_mat[j][k][l])
						{
							DEBUG_NOTICE("sys_intp_coef_mat: %d %d | %ld %ld %ld | %x\n",
										 j,
										 k,
							             x_term_degree_table[l],
								         y_term_degree_table[l],
								         z_term_degree_table[l],
								         sys_intp_coef_mat[j][k][l]);
						}
					}
				}
#endif
			}
		}		
	}

#if 1
	min_deg = 65536;
	for(j = 0; j < KOT_INTP_POLY_NUM; j++)
	{
		deg_tmp = sys_intp_coef_mat_deg_cal(j);
		dev_store[j] = deg_tmp;
		if(min_deg > deg_tmp)
		{
			min_deg = deg_tmp;
			min_idx = j;
		}

#if (1 == TEST_MODE)		
		for(k = 0; k < KOT_INTP_POLY_NUM; k++)
		{
			for(l = 0; l < MAX_POLY_TERM_SIZE; l++)
			{
				if(0xFF != sys_intp_coef_mat[j][k][l])
				{
					DEBUG_NOTICE("sel_sys_intp_coef_mat: %d %d | %ld %ld %ld | %x\n",
								 j,
								 k,
					             x_term_degree_table[l],
						         y_term_degree_table[l],
						         z_term_degree_table[l],
						         sys_intp_coef_mat[j][k][l]);
				}
			}
		}
#endif		
	}
	//memcpy(min_intp_poly, intp_poly_coef[min_idx], sizeof(unsigned char) * MAX_POLY_TERM_SIZE);

	//min_idx = sys_min_mul_poly_idx_determine();

	for(j = 0; j < KOT_INTP_POLY_NUM; j++)
	{
		if(min_idx != j)
		{
			continue;
		}

		sys_poly_mat_mul(min_idx);
		sys_mul_poly_recover();

		for(i = 0; i < CODEWORD_LEN; i++)
		{
			dev_val = sys_hasse_full_dev_cal(min_idx, i, gf_add(tst_vct[0][i], sys_ret_cwd[i]));
			DEBUG_NOTICE("sys_final_has_dev_test: %ld | %x %x | %x %x | %x | %x | %x\n",
						 j,
						 af_pnt[i][0],
						 af_pnt[i][1],
						 tst_vct[0][i],
						 sys_ret_cwd[i],
						 recv_poly[i],
						 gf_add(tst_vct[0][i], sys_ret_cwd[i]),
						 dev_val);
		}

		her_fac(min_intp_poly, tv_est_msg[0], tv_est_cwd[0]);

		for(i = 0; i < CODEWORD_LEN; i++)
		{
			tv_est_cwd[0][i] = gf_add(tv_est_cwd[0][i], sys_ret_cwd[i]);
			DEBUG_NOTICE("kot_sys_tv_est_cwd: %d | %d | %ld | %x %x %x\n", j, tv_err[0], i, tv_est_cwd[0][i], cwd_poly[i], recv_poly[i]);
			if((tv_est_cwd[0][i] != cwd_poly[i])
				&& (tv_err[0] <= ((CODEWORD_LEN - MESSAGE_LEN - (GF_Q * (GF_Q - 1) / 2) + 1 - 1) / 2)))
			{
				DEBUG_SYS("sys_dec_err: %ld | %ld | %x %x | %x %x\n",
				          tv_err[0],
				          i,
				          gf_add(sys_ret_cwd[i], tst_vct[0][i]),
				          gf_add(sys_ret_cwd[i], cwd_poly[i]),
				          tv_est_cwd[0][i],
				          cwd_poly[i]);
			}
		}
	}
#endif
	return 0;
}
