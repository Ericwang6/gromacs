#include "gmxpre.h"

#include "config.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/ns.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "deepmd/DeepPot.h"

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

const float ANG2NM = 0.1;
const float EV2KJ = 96.48533132;
const float EVANG2KJNM = 964.8533132;

#ifdef HIGH_PREC
typedef double VALUETYPE;
#else
typedef float VALUETYPE;
#endif

deepmd::DeepPot dp;

void init_dp(t_QMrec* qm)
{
    std::cout << "Init DeepPotential..." << std::endl;
    snew(qm->dp_model_dir, 200);
    qm->dp_model_dir = getenv("GMX_DP_MODEL_DIR");
    std::cout << qm->dp_model_dir << std::endl;
    dp.init(qm->dp_model_dir);
}

real call_dp(const t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
    t_QMMMrec* QMMMrec;
    QMMMrec = fr->qr;

    /* Index init */
    int ii, jj;

    /* DeePMD */
    std::vector<VALUETYPE > dbox = {}; // no pbc
    std::vector<VALUETYPE > dcoord;
    std::vector<VALUETYPE > dforce;
    std::vector<VALUETYPE > dvirial;
    std::vector<int > dtype;
    dcoord.resize(qm->nrQMatoms * 3);
    dforce.resize(qm->nrQMatoms * 3);
    dvirial.resize(9);

    /* read type.raw file */
    std::ifstream typeFile("type.raw");
    int val;
    if (typeFile.is_open())
    {
        for (ii = 0; ii < qm->nrQMatoms; ii++)
        {
            typeFile >> val;
            dtype.push_back(val);
        }
    }
    else
    {
        gmx_fatal(FARGS, "Fail to read type file");
    }
    /* read type.raw file */

    /* Init coords */
    for (ii = 0; ii < qm->nrQMatoms; ii++)
    {
        dcoord[ii * DIM]     = qm->xQM[ii][XX] / ANG2NM;
        dcoord[ii * DIM + 1] = qm->xQM[ii][YY] / ANG2NM;
        dcoord[ii * DIM + 2] = qm->xQM[ii][ZZ] / ANG2NM;
    }

    /* call dp */
    double ener;
    dp.compute(ener, dforce, dvirial, dcoord, dtype, dbox);
    real QMener = ener * EV2KJ;

    /* for debugger */
    std::cout << "QM energy:" << QMener << std::endl;
    std::cout << "QM coords (" << qm->nrQMatoms << " atoms): " << std::endl;
    for (ii = 0; ii < qm->nrQMatoms; ii++)
    {
        for (jj = 0; jj < DIM; jj++)
        {
            std::cout << dcoord[ii * DIM + jj] << " ";
        }
        std::cout << std::endl;
    }

    /* fake forces */
    // for (ii = 0; ii < qm->nrQMatoms; ii++)
    // {
    //     for (jj = 0; jj < DIM; jj++)
    //     {
    //         dforce[ii * DIM + jj] = 0.0;
    //     }
    // }

    /* DP forces */
    for (ii = 0; ii < qm->nrQMatoms; ii++)
    {
        for (jj = 0; jj < DIM; jj++)
        {
            f[ii][jj]      = dforce[ii * DIM + jj] * EVANG2KJNM;
            fshift[ii][jj] = dforce[ii * DIM + jj] * EVANG2KJNM;
        }
    }

    for (ii = 0; ii < mm->nrMMatoms; ii++)
    {
        for (jj = 0; jj < DIM; jj++)
        {
            f[ii + qm->nrQMatoms][jj]      = 0.0;
            fshift[ii + qm->nrQMatoms][jj] = 0.0;
        }
    }
    return (QMener);
}