/* Referencing METIS_MUMPS_PartGraphKway is what makes this a real test:
 * it lives in kmetis.c.o, and `nm -u libmetis_mumps.a` shows that object
 * undefined-referencing log and pow (and the objects it drags in adding
 * sqrt and powf). Pulling it out of the archive is therefore enough to
 * fail the link if libm is missing from eplinalg::metis's exported
 * interface. Calling it is not necessary and would need a real graph.
 *
 * The address must be taken from a *global* array with external linkage.
 * The obvious spelling -- a local `void *p = (void *)f; if (!p) ...` --
 * is folded away entirely: the compiler knows a function address is never
 * null, drops the comparison, and with it the only relocation naming the
 * symbol, so the archive member is never pulled in and the test silently
 * passes with libm missing. A global array cannot be elided that way.
 */
#include <stdio.h>

extern int METIS_MUMPS_PartGraphKway();
extern int METIS_MUMPS_SetDefaultOptions();

void *eplinalg_metis_syms[] = {
    (void *)METIS_MUMPS_PartGraphKway,
    (void *)METIS_MUMPS_SetDefaultOptions,
};

int main(void)
{
    size_t n = sizeof(eplinalg_metis_syms) / sizeof(eplinalg_metis_syms[0]);
    size_t i;

    for (i = 0; i < n; ++i) {
        if (eplinalg_metis_syms[i] == 0) {
            printf("FAIL: METIS symbol %zu resolved to NULL\n", i);
            return 1;
        }
    }

    printf("OK: C-only consumer linked and ran (%zu METIS symbols)\n", n);
    return 0;
}
