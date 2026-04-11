/* blacs_pinfo_.c -- multifloats override of the BLACS pinfo entry point.
 *
 * This file replaces the upstream blacs_pinfo_.c verbatim-plus-one-line:
 * after MPI_Init has been called (either by the user or by the existing
 * lazy-init here), we also call multifloats_mpi_init(), which registers
 * the MPI_FLOAT64X2 / MPI_COMPLEX128X2 datatypes and the MPI_DD_SUM /
 * MPI_ZZ_SUM / MPI_DD_AMX / MPI_DD_AMN / MPI_ZZ_AMX / MPI_ZZ_AMN
 * user-defined ops that the migrated BLACS reductions rely on.
 *
 * Every ScaLAPACK user program calls blacs_pinfo (directly or via the
 * lazy call inside blacs_gridmap_) before any BLACS communication, so
 * this single hook is sufficient.  multifloats_mpi_init is idempotent
 * -- repeat calls after the first are a one-branch no-op -- so mixing
 * multifloats BLACS with single-precision or other precision variants
 * in the same program is safe.
 */
#include "Bdef.h"
#include "multifloats_c.h"

#if (INTFACE == C_CALL)
void Cblacs_pinfo(Int *mypnum, Int *nprocs)
#else
F_VOID_FUNC blacs_pinfo_(Int *mypnum, Int *nprocs)
#endif
{
   Int ierr;
   extern Int BI_Iam, BI_Np;
   MpiInt flag, Iam = BI_Iam, Np = BI_Np;
   MpiInt argc=0;
   char **argv=NULL;
   if (BI_COMM_WORLD == NULL)
   {
      MPI_Initialized(&flag);

      if (!flag)
         ierr = MPI_Init(&argc,&argv);  // call Init and ignore argc and argv

      BI_COMM_WORLD = (Int *) malloc(sizeof(Int));
      *BI_COMM_WORLD = MPI_Comm_c2f(MPI_COMM_WORLD);
   }
   MPI_Comm_size(MPI_COMM_WORLD, &Np);
   MPI_Comm_rank(MPI_COMM_WORLD, &Iam);
   *mypnum = BI_Iam = Iam;
   *nprocs = BI_Np  = Np;

   /* Multifloats one-time MPI handle registration. Idempotent. */
   multifloats_mpi_init();
}
