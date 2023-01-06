#ifndef CFG_DECODING_H
#define CFG_DECODING_H

#define GF_Q			2
#if (2 == GF_Q)
#define GF_FIELD		4
#define MESSAGE_LEN     5
#endif
#if (3 == GF_Q)
#define GF_FIELD        8
#define MESSAGE_LEN     3
#endif
#if (4 == GF_Q)
#define GF_FIELD        16
#define MESSAGE_LEN     7
#endif
#if (6 == GF_Q)
#define GF_FIELD        64
#define MESSAGE_LEN     21
#endif
#if (8 == GF_Q)
#define GF_FIELD        256
#define MESSAGE_LEN     223
#endif
#define CODEWORD_LEN    (GF_Q * GF_Q * GF_Q)   

#define GF_CAL_COUNT	 1

#define SYS_ENC			 0

#define TEST_MODE		 1

#define OUTPUT_LOG		 0

#define ETA				 2

#define MAX_DEGREE		 8
#define MAX_POLY_TERM_SIZE	((MAX_DEGREE + 1) * (MAX_DEGREE + 1) * (MAX_DEGREE + 1))
#define KOT_INTP_POLY_NUM	(2 * GF_Q)

#endif
