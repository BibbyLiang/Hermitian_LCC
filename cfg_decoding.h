#ifndef CFG_DECODING_H
#define CFG_DECODING_H

#define GF_Q			4
#if (2 == GF_Q)
#define GF_FIELD		4
#define MESSAGE_LEN     5
#endif
#if (3 == GF_Q)
#define GF_FIELD        8
#define MESSAGE_LEN     5
#endif
#if (4 == GF_Q)
#define GF_FIELD        16
#define MESSAGE_LEN     42
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

#define TEST_MODE		 0

#define OUTPUT_LOG		 0

#define ETA				 1

#define MAX_DEGREE		 GF_FIELD
#define CFG_DYM_SIZE	 1
#if (0 == CFG_DYM_SIZE)
#define MAX_POLY_TERM_SIZE	((MAX_DEGREE + 1) * (MAX_DEGREE + 1) * (MAX_DEGREE + 1))
#else
#define X_MAX_SIZE	     (MAX_DEGREE + 1)
#define Y_MAX_SIZE		 (MAX_DEGREE + 1)
#define Z_MAX_SIZE		 2
#define MAX_POLY_TERM_SIZE	(X_MAX_SIZE * Y_MAX_SIZE * Z_MAX_SIZE)
#endif
#define CFG_QUICK_POLY_SEARCH	1
#define KOT_INTP_POLY_NUM	(2 * GF_Q)

#define CFG_FAC_FREE	 0
#define CFG_STORE_DEV	 1

#define CFG_RET			 0

#endif
