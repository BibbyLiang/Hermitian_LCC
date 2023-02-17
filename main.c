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

	runtime = 0;
	err_cnt = 0;

	DEBUG_SYS("init_simulation OK\n");

	return;
}

void main()
{
	long long i = 0;
	long long sim_cnt = 100, monitor_cnt = 1;

	init_simulation();

	mod_init();
	tst_vct_num_cal();
	tst_vct_init();
	
	affine_points_cal();
	term_degree_table_init();
	pole_basis_cal();

	for(i = 0; i < sim_cnt; i++)
	{
		if(0 == (i % monitor_cnt))
		{
			DEBUG_SYS("sim: %ld / %ld, err: %ld\n", i, sim_cnt, err_cnt);
		}

		rnd_msg_gen();
		her_encoding(msg_poly, cwd_poly);
		trans_over_chnl();
		chnl_rel_cal(recv_seq, symbol_num);

		cnt_switch = 1;
		start = clock();

#if (1 == CFG_RET)
		re_encoding_transform();
#endif		

		tst_vct_form();

		koetter_interpolation_hermitian();

		stop = clock();
		runtime = runtime + (stop - start) / 1000.0000;
		cnt_switch = 0;

#if (1 == CFG_FAC_FREE)
		check_result_cwd(cwd_poly, est_cwd_poly);
#else
#if (0 == CFG_RET)
		check_result_msg(msg_poly, est_msg_poly);
#else
		check_result_cwd(cwd_poly, est_cwd_poly);
#endif
#endif
	}
	DEBUG_SYS("sim: %ld, err: %ld\n", sim_cnt, err_cnt);
	DEBUG_SYS("Decoding Complexity: %f %f\n", (float)(add_cnt / sim_cnt), (float)(mul_cnt / sim_cnt));
	DEBUG_SYS("Decoding Measured Time: %f\n", runtime / ((float)sim_cnt));

	tst_vct_exit();
	mod_exit();
	
	return;
}
