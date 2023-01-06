#include <stdlib.h>
#include <string.h>
#include "cfg_decoding.h"
#include "debug_info.h"
#include "gf_cal.h"
#include "encoding.h"
#include "rnd.h"
#include "mod.h"
#include "time.h"
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
	init_simulation();

	mod_init();

	affine_points_cal();
	pole_basis_cal();
	rnd_msg_gen();
	her_encoding(msg_poly);
	trans_over_chnl();
	chnl_rel_cal(recv_seq, symbol_num);
	
	tst_vct_num_cal();
	tst_vct_init();
	tst_vct_form();
	
	koetter_interpolation_hermitian();
	
	tst_vct_exit();
	mod_exit();
	
	return;
}
