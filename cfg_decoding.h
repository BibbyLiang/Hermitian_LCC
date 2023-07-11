#ifndef CFG_DECODING_H
#define CFG_DECODING_H

#define EARLY_TERMINATION	0
#define ET_NUM				150

#define GF_Q			4
#if (2 == GF_Q)
#define GF_FIELD		4
#define MESSAGE_LEN     5
#endif
#if (3 == GF_Q)
#define GF_FIELD        8
#define MESSAGE_LEN     4
#endif
#if (4 == GF_Q)
#define GF_FIELD        16
#define MESSAGE_LEN     39
#endif
#if (6 == GF_Q)
#define GF_FIELD        64
#define MESSAGE_LEN     409
#endif
#if (8 == GF_Q)
#define GF_FIELD        256
#define MESSAGE_LEN     223
#endif
#define CODEWORD_LEN    (GF_Q * GF_Q * GF_Q)   

#define GF_CAL_COUNT	 1

#define SYS_ENC			 0

#define CFG_PRG_DECODING 0

#define CFG_PRG_RET_ET	 0

#define TEST_MODE		 0

#define FIX_INPUT_DBG	 0

#define TV_TEST			 0

#define OUTPUT_LOG		 0

#define ETA				 0

#define MAX_DEGREE		 GF_FIELD
#define CFG_DYM_SIZE	 1
#if (0 == CFG_DYM_SIZE)
#define MAX_POLY_TERM_SIZE	((MAX_DEGREE + 1) * (MAX_DEGREE + 1) * (MAX_DEGREE + 1))
#else
#define X_MAX_SIZE	     (POLY_TERM_SIZE_FACTOR * (MAX_DEGREE + 1))
#define Y_MAX_SIZE		 (POLY_TERM_SIZE_FACTOR * (MAX_DEGREE + 1))
#define Z_MAX_SIZE		  2
#define MAX_POLY_TERM_SIZE	(X_MAX_SIZE * Y_MAX_SIZE * Z_MAX_SIZE)
#endif
#define CFG_QUICK_POLY_SEARCH	1
#define KOT_INTP_POLY_NUM	   (2 * GF_Q)
#define BR_BASIS_NUM		   (2 * GF_Q)

#define CFG_FIND_LATER_POLY 1

#define CFG_FAC_FREE	 	1
#define CFG_STORE_DEV	 	CFG_FAC_FREE//it is efficient
#define FAC_FREE_ERR	 	0
#define DEV_RECORD		 	0
#define CFG_CWD_DIM_CHECK	1

#define CFG_RET			 1
#define CFG_FAST_RET	 1
#define CFG_Y_RET_STORE  CFG_FAST_RET//there are some bugs, it cannot be set as 0 independently
#define CFG_RET_ETA_OPT  0
#define CFG_RET_L		 ((MESSAGE_LEN - (GF_Q) * (GF_Q - 1) / 2 + 0) / GF_Q * GF_Q + 0)
#define CFG_NREL_NO_RET	 1 //it must set to 1, or some bug may be occured.
#define CFG_P_NREL_SORT  0//there are some bugs.

#define CFG_INTP_ORDER_ERR 0

#define CFG_BR			  1
#define CFG_ACD_BR		  0
#define CFG_TV_K_POLY_IND 1//this must set to 1, or some bug may be occured.

#define CFG_SYS_GEN		  1
#define CFG_SYS_GEN_POLY  GF_Q
#define CFG_FAST_FULL_GAU 1
#define CFG_ZERO_COL_CHECK_GAU_ELM 1

#if (0 == CFG_BR)
#define POLY_TERM_SIZE_FACTOR	1
#else
#define POLY_TERM_SIZE_FACTOR	2
#endif

#endif
