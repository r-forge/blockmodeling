#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void critFun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void optPar(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void optParMulti(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void parArr2Vec(void *, void *, void *, void *, void *);
extern void parVec2Arr(void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(regd)(void *, void *, void *, void *, void *);
extern void F77_NAME(regdne)(void *, void *, void *, void *, void *);
extern void F77_NAME(regdow)(void *, void *, void *, void *, void *);
extern void F77_NAME(regdowne)(void *, void *, void *, void *, void *);
extern void F77_NAME(rege)(void *, void *, void *, void *, void *);
extern void F77_NAME(regene)(void *, void *, void *, void *, void *);
extern void F77_NAME(regenm)(void *, void *, void *, void *, void *);
extern void F77_NAME(regenmdiag)(void *, void *, void *, void *, void *);
extern void F77_NAME(regenmne)(void *, void *, void *, void *, void *);
extern void F77_NAME(regeow)(void *, void *, void *, void *, void *);
extern void F77_NAME(regeowne)(void *, void *, void *, void *, void *);
extern void F77_NAME(regeownm)(void *, void *, void *, void *, void *);
extern void F77_NAME(regeownmdiag)(void *, void *, void *, void *, void *);
extern void F77_NAME(regeownmne)(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"critFun",     (DL_FUNC) &critFun,     30},
    {"optPar",      (DL_FUNC) &optPar,      37},
    {"optParMulti", (DL_FUNC) &optParMulti, 42},
    {"parArr2Vec",  (DL_FUNC) &parArr2Vec,   5},
    {"parVec2Arr",  (DL_FUNC) &parVec2Arr,   5},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"regd",         (DL_FUNC) &F77_NAME(regd),         5},
    {"regdne",       (DL_FUNC) &F77_NAME(regdne),       5},
    {"regdow",       (DL_FUNC) &F77_NAME(regdow),       5},
    {"regdowne",     (DL_FUNC) &F77_NAME(regdowne),     5},
    {"rege",         (DL_FUNC) &F77_NAME(rege),         5},
    {"regene",       (DL_FUNC) &F77_NAME(regene),       5},
    {"regenm",       (DL_FUNC) &F77_NAME(regenm),       5},
    {"regenmdiag",   (DL_FUNC) &F77_NAME(regenmdiag),   5},
    {"regenmne",     (DL_FUNC) &F77_NAME(regenmne),     5},
    {"regeow",       (DL_FUNC) &F77_NAME(regeow),       5},
    {"regeowne",     (DL_FUNC) &F77_NAME(regeowne),     5},
    {"regeownm",     (DL_FUNC) &F77_NAME(regeownm),     5},
    {"regeownmdiag", (DL_FUNC) &F77_NAME(regeownmdiag), 5},
    {"regeownmne",   (DL_FUNC) &F77_NAME(regeownmne),   5},
    {NULL, NULL, 0}
};

void R_init_blockmodelingRoxygen2(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
