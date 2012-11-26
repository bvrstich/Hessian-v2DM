//nog enkele definities:
#ifdef PQ

#define __Q_CON

#endif

#ifdef PQG

#define __Q_CON
#define __G_CON

#include "PHM.h"

#endif

#ifdef PQGT1

#define __Q_CON
#define __G_CON
#define __T1_CON

#include "PHM.h"
#include "DPM.h"

#endif

#ifdef PQGT2

#define __Q_CON
#define __G_CON
#define __T2_CON

#include "PHM.h"
#include "PPHM.h"

#endif

#ifdef PQGT

#define __Q_CON
#define __G_CON
#define __T1_CON
#define __T2_CON

#include "PHM.h"
#include "DPM.h"
#include "PPHM.h"

#endif

#include "lapack.h"
#include "Matrix.h"
#include "RecMat.h"
#include "Vector.h"
#include "TPM.h"
#include "SPM.h"
#include "Hessian.h"
#include "HessBar.h"

#include "SUP.h"
#include "EIG.h"

#include "Tools.h"

#include "Newton.h"
