#include <stdlib.h>
#include <string.h>
#include "cfg_decoding.h"
#include "debug_info.h"
#include "gf_cal.h"
#include "encoding.h"
#include "rnd.h"
#include "mod.h"
#include "time.h"
#include "interpolation.h"
#include "factorization.h"
#include "re_encoding.h"
#include "lcc_decoding.h"

void init_simulation()
{
	srand(time(NULL));
	init_genrand((long)(time(NULL)));

	DEBUG_SYS("init_simulation OK\n");

	return;
}

void main()
{
	long long i = 0;
	long long sim_cnt = 1, monitor_cnt = 1;

	init_simulation();

	mod_init();
	tst_vct_num_cal();
	tst_vct_init();

	for(i = 0; i < sim_cnt; i++)
	{
		if(0 == (i % monitor_cnt))
		{
			DEBUG_SYS("sim: %ld / %ld, err: %ld\n", i, sim_cnt, err_cnt);
		}
	
		affine_points_cal();
		term_degree_table_init();
		pole_basis_cal();
		rnd_msg_gen();
		her_encoding(msg_poly);
		trans_over_chnl();
		chnl_rel_cal(recv_seq, symbol_num);

		/*test*/
		her_lagrange_poly_construct();

		tst_vct_form();

		cnt_switch = 1;

		koetter_interpolation_hermitian();

		cnt_switch = 0;

#if (1 == CFG_FAC_FREE)
		check_result_cwd(cwd_poly, est_cwd_poly);
		//check_result_msg(msg_poly, est_msg_poly);
#else
		check_result_msg(msg_poly, est_msg_poly);
#endif
	}
	DEBUG_SYS("sim: %ld, err: %ld\n", sim_cnt, err_cnt);
	DEBUG_SYS("Decoding Complexity: %f %f\n", (float)(add_cnt / sim_cnt), (float)(mul_cnt / sim_cnt));

	tst_vct_exit();
	mod_exit();
	
	return;
}
