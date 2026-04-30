/*
 * C-side MPI stubs that supplement upstream's libmpiseq for fully serial
 * (no mpiexec, no real MPI runtime) test executables.
 *
 * Upstream's `external/MUMPS_5.8.2/libseq/mpic.c` only stubs four
 * routines: `MPI_Init`, `MPI_Comm_rank`, `MPI_Comm_size`,
 * `MPI_Comm_f2c`. The migrated standard archives in this build
 * (`libblacs.a`, `libqblacs.a`, `libpblas.a`, `libqpblas.a`) are C code
 * that calls many more — `MPI_Send`, `MPI_Recv`, `MPI_Isend`,
 * `MPI_Type_match_size`, `MPI_Comm_split`, etc. Linking those archives
 * against libmpiseq alone leaves the C-MPI calls unresolved. This file
 * provides single-process implementations of every additional symbol
 * the BLACS / PBLAS C side references.
 *
 * All stubs assume rank=0, size=1, and either no-op (for collectives,
 * since there's nothing to coordinate with) or print a "should not be
 * called" error and abort (for point-to-point, which a single-rank
 * solver should never invoke). Matches the contract upstream's
 * Fortran-side mpi.f stubs use.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>

#include "mpi.h"  /* libmpiseq's mpi.h — not real MPI */

#ifndef LIBSEQ_CALL
#define LIBSEQ_CALL
#endif

/* libmpiseq's mpi.h only defines MPI_Comm. Provide the additional
 * opaque types BLACS / PBLAS reference. All are integer-shaped so
 * single-rank stubs can pass them through trivially. */
typedef int        MPI_Datatype;
typedef int        MPI_Op;
typedef int        MPI_Request;
typedef intptr_t   MPI_Aint;
typedef struct { int MPI_SOURCE; int MPI_TAG; int MPI_ERROR; int count; } MPI_Status;
typedef void (MPI_User_function)(void *invec, void *inoutvec, int *len, MPI_Datatype *dt);

static void mpiseq_should_not_be_called(const char *fn)
{
    fprintf(stderr,
            "Error. %s should not be called from a single-rank MUMPS "
            "test built against libmpiseq.\n", fn);
    abort();
}

/* ── Lifecycle ──────────────────────────────────────────────────────── */

/* MPI_Init / MPI_Finalize / MPI_Comm_rank / MPI_Wtime / MPI_Comm_f2c are
 * already supplied by upstream's libmpiseq mpic.c — don't redefine. */

int LIBSEQ_CALL MPI_Initialized(int *flag)               { *flag = 1; return 0; }
int LIBSEQ_CALL MPI_Comm_size(MPI_Comm comm, int *size)  { (void) comm; *size = 1; return 0; }
/* MPI_Comm_c2f: convert C handle to Fortran handle. Single-rank: any
 * non-zero handle maps to the upstream USE_COMM_WORLD constant
 * (-987654) which the F77 layer interprets as MPI_COMM_WORLD; the
 * pass-through (LIBSEQ_INT)comm is also valid because libseq's
 * MPI_Comm IS LIBSEQ_INT. */
int LIBSEQ_CALL MPI_Comm_c2f(MPI_Comm comm)              { return (int) comm; }
int LIBSEQ_CALL MPI_Abort(MPI_Comm comm, int errorcode)
{
    (void) comm;
    fprintf(stderr, "MPI_Abort called with errorcode=%d\n", errorcode);
    exit(errorcode);
    return 0;
}

/* ── Communicator / group queries ────────────────────────────────────── */

int LIBSEQ_CALL MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)         { *newcomm = comm; return 0; }
int LIBSEQ_CALL MPI_Comm_free(MPI_Comm *comm)                          { (void) comm; return 0; }
int LIBSEQ_CALL MPI_Comm_create(MPI_Comm comm, void *group, MPI_Comm *newcomm)
{
    (void) group; *newcomm = comm; return 0;
}
int LIBSEQ_CALL MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
    (void) color; (void) key; *newcomm = comm; return 0;
}
int LIBSEQ_CALL MPI_Comm_group(MPI_Comm comm, void *group)             { (void) comm; (void) group; return 0; }
int LIBSEQ_CALL MPI_Comm_get_attr(MPI_Comm c, int kv, void *v, int *f) { (void) c; (void) kv; (void) v; *f = 0; return 0; }
int LIBSEQ_CALL MPI_Group_incl(void *g, int n, int *r, void *ng)       { (void) g; (void) n; (void) r; (void) ng; return 0; }
int LIBSEQ_CALL MPI_Group_free(void *g)                                { (void) g; return 0; }

/* ── Collectives — single rank, so the "operation" is a no-op or copy ─ */

int LIBSEQ_CALL MPI_Barrier(MPI_Comm comm)                             { (void) comm; return 0; }

int LIBSEQ_CALL MPI_Bcast(void *buf, int count, MPI_Datatype dt, int root, MPI_Comm comm)
{
    (void) buf; (void) count; (void) dt; (void) root; (void) comm; return 0;
}

/* Single-rank reduce SEMANTICALLY equals `memcpy(recvbuf, sendbuf,
 * count * sizeof(element))`, but libmpiseq's mpi.h doesn't expose
 * datatype-size constants — MPI_Datatype is just an opaque int and
 * MPI_INT / MPI_DOUBLE / MPI_DOUBLE_COMPLEX / MPI_DOUBLE_PRECISION /
 * etc. integer values aren't defined on the C side. Using a
 * hard-coded element size truncates quad-precision reductions
 * (REAL(KIND=16) = 16 B, COMPLEX(KIND=16) = 32 B) and over-copies
 * 4-byte INT reductions past sendbuf.
 *
 * The honest stub aborts on any call. tests/mumps's standard build
 * (mpiexec + real MPI) never reaches here. The fully-sequential
 * libmpiseq link path (B3 in tests/mumps/TODO.md) would invoke this
 * via BLACS/PBLAS C code, which is currently blocked anyway by the
 * Q-prefixed scalapack symbol gap. When B3 lands, this stub needs a
 * datatype→size table built by mirroring mpich's MPI_INT etc. values
 * or by extending libmpiseq's mpi.h with size-encoding macros. */
int LIBSEQ_CALL MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
                           MPI_Datatype dt, MPI_Op op, int root, MPI_Comm comm)
{
    (void) sendbuf; (void) recvbuf; (void) count;
    (void) dt; (void) op; (void) root; (void) comm;
    mpiseq_should_not_be_called(
        "MPI_Reduce (libmpiseq stub cannot infer datatype size — see "
        "tests/mumps/c/mpiseq_c_stubs.c)");
    return 0;
}

int LIBSEQ_CALL MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                              MPI_Datatype dt, MPI_Op op, MPI_Comm comm)
{
    (void) sendbuf; (void) recvbuf; (void) count;
    (void) dt; (void) op; (void) comm;
    mpiseq_should_not_be_called(
        "MPI_Allreduce (libmpiseq stub cannot infer datatype size — see "
        "tests/mumps/c/mpiseq_c_stubs.c)");
    return 0;
}

/* ── Point-to-point — illegal for a single rank ───────────────────────── */

int LIBSEQ_CALL MPI_Send(const void *b, int c, MPI_Datatype d, int dst, int tag, MPI_Comm cm)
{
    (void) b; (void) c; (void) d; (void) dst; (void) tag; (void) cm;
    mpiseq_should_not_be_called("MPI_Send");
    return 0;
}
int LIBSEQ_CALL MPI_Rsend(const void *b, int c, MPI_Datatype d, int dst, int tag, MPI_Comm cm)
{
    (void) b; (void) c; (void) d; (void) dst; (void) tag; (void) cm;
    mpiseq_should_not_be_called("MPI_Rsend");
    return 0;
}
int LIBSEQ_CALL MPI_Recv(void *b, int c, MPI_Datatype d, int src, int tag,
                         MPI_Comm cm, MPI_Status *st)
{
    (void) b; (void) c; (void) d; (void) src; (void) tag; (void) cm; (void) st;
    mpiseq_should_not_be_called("MPI_Recv");
    return 0;
}
int LIBSEQ_CALL MPI_Isend(const void *b, int c, MPI_Datatype d, int dst, int tag,
                          MPI_Comm cm, MPI_Request *rq)
{
    (void) b; (void) c; (void) d; (void) dst; (void) tag; (void) cm;
    if (rq) *rq = 0;
    mpiseq_should_not_be_called("MPI_Isend");
    return 0;
}
int LIBSEQ_CALL MPI_Irecv(void *b, int c, MPI_Datatype d, int src, int tag,
                          MPI_Comm cm, MPI_Request *rq)
{
    (void) b; (void) c; (void) d; (void) src; (void) tag; (void) cm;
    if (rq) *rq = 0;
    mpiseq_should_not_be_called("MPI_Irecv");
    return 0;
}
int LIBSEQ_CALL MPI_Sendrecv(const void *sb, int sc, MPI_Datatype sd, int dst, int stag,
                             void *rb, int rc, MPI_Datatype rd, int src, int rtag,
                             MPI_Comm cm, MPI_Status *st)
{
    (void) sb; (void) sc; (void) sd; (void) dst; (void) stag;
    (void) rb; (void) rc; (void) rd; (void) src; (void) rtag;
    (void) cm; (void) st;
    mpiseq_should_not_be_called("MPI_Sendrecv");
    return 0;
}
int LIBSEQ_CALL MPI_Waitall(int count, MPI_Request *reqs, MPI_Status *stats)
{
    (void) count; (void) reqs; (void) stats; return 0;
}
int LIBSEQ_CALL MPI_Testall(int count, MPI_Request *reqs, int *flag, MPI_Status *stats)
{
    (void) count; (void) reqs; (void) stats; *flag = 1; return 0;
}

/* ── Datatype handles — opaque integers for libmpiseq's purposes ─────── */

int LIBSEQ_CALL MPI_Type_vector(int c, int bl, int st, MPI_Datatype old, MPI_Datatype *new_)
{
    (void) c; (void) bl; (void) st; (void) old; *new_ = 0; return 0;
}
int LIBSEQ_CALL MPI_Type_indexed(int c, int *bls, int *disps, MPI_Datatype old, MPI_Datatype *new_)
{
    (void) c; (void) bls; (void) disps; (void) old; *new_ = 0; return 0;
}
int LIBSEQ_CALL MPI_Type_create_struct(int c, int *bls, MPI_Aint *disps,
                                       MPI_Datatype *types, MPI_Datatype *new_)
{
    (void) c; (void) bls; (void) disps; (void) types; *new_ = 0; return 0;
}
int LIBSEQ_CALL MPI_Type_commit(MPI_Datatype *t)                       { (void) t; return 0; }
int LIBSEQ_CALL MPI_Type_free(MPI_Datatype *t)                         { (void) t; return 0; }
int LIBSEQ_CALL MPI_Type_match_size(int typeclass, int size, MPI_Datatype *t)
{
    (void) typeclass; (void) size; *t = 0; return 0;
}

int LIBSEQ_CALL MPI_Pack_size(int incount, MPI_Datatype dt, MPI_Comm cm, int *size)
{
    (void) dt; (void) cm; *size = incount; return 0;
}
int LIBSEQ_CALL MPI_Pack(const void *inbuf, int incount, MPI_Datatype dt,
                         void *outbuf, int outsize, int *position, MPI_Comm cm)
{
    (void) dt; (void) outsize; (void) cm;
    if (inbuf && outbuf && incount > 0)
        memcpy((char *) outbuf + *position, inbuf, (size_t) incount);
    if (position) *position += incount;
    return 0;
}
int LIBSEQ_CALL MPI_Unpack(const void *inbuf, int insize, int *position,
                           void *outbuf, int outcount, MPI_Datatype dt, MPI_Comm cm)
{
    (void) insize; (void) dt; (void) cm;
    if (inbuf && outbuf && outcount > 0)
        memcpy(outbuf, (const char *) inbuf + *position, (size_t) outcount);
    if (position) *position += outcount;
    return 0;
}

int LIBSEQ_CALL MPI_Op_create(MPI_User_function *fn, int commute, MPI_Op *op)
{
    (void) fn; (void) commute; *op = 0; return 0;
}
int LIBSEQ_CALL MPI_Op_free(MPI_Op *op)                                { *op = 0; return 0; }

int LIBSEQ_CALL MPI_Error_class(int errcode, int *errclass)            { *errclass = errcode; return 0; }

/* MPI_Wtime is already supplied by libmpiseq mpic.c. */
