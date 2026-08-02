/* Calls into a Fortran archive from C with no Fortran driver on the
 * link line, to prove the package supplied the Fortran runtime itself.
 *
 * ddot_ is called for a checkable result. dgemm_ is referenced only for
 * its address, and that reference is the part that matters: dgemm calls
 * XERBLA, and xerbla's WRITE and STOP pull in _gfortran_st_write and
 * _gfortran_stop_string. Without libgfortran on the line the link fails
 * there. (ddot alone would not — it touches no runtime entry point.)
 *
 * As in ../c/main.c, the address must live in a global array with
 * external linkage; a local pointer compared against null is folded away
 * and takes the relocation — and therefore the archive member, and
 * therefore the whole test — with it.
 */
#include <stdio.h>

extern double ddot_(const int *n, const double *x, const int *incx,
                    const double *y, const int *incy);
extern void dgemm_();

void *eplinalg_blas_syms[] = {
    (void *)dgemm_,
};

int main(void)
{
    const int n = 3, inc = 1;
    const double x[3] = {1.0, 2.0, 3.0};
    const double y[3] = {4.0, 5.0, 6.0};
    double d;

    if (eplinalg_blas_syms[0] == 0) {
        printf("FAIL: dgemm_ resolved to NULL\n");
        return 1;
    }

    d = ddot_(&n, x, &inc, y, &inc);
    if (d != 32.0) {
        printf("FAIL: ddot_ returned %g, expected 32\n", d);
        return 1;
    }

    printf("OK: C-only consumer used a Fortran archive via "
           "EPLINALG_FORTRAN_TAG (ddot=%g)\n", d);
    return 0;
}
