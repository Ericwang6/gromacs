#include "gromacs/mdlib/qmmm.h"
#include "deepmd/DeepPot.h"

#ifndef GMX_MDLIB_QM_DP_H
#define GMX_MDLIB_QM_DP_H

extern deepmd::DeepPot dp;
void init_dp(t_QMrec* qm);
real call_dp(const t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[]);

#endif