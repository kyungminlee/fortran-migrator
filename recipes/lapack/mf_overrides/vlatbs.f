*> \brief \b VLATBS solves a triangular banded system of equations.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download VLATBS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/vlatbs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/vlatbs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/vlatbs.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE VLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB, X,
*                          SCALE, CNORM, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, NORMIN, TRANS, UPLO
*       INTEGER            INFO, KD, LDAB, N
*       TYPE(real64x2)   SCALE
*       ..
*       .. Array Arguments ..
*       TYPE(real64x2)   CNORM( * )
*       TYPE(cmplx64x2)         AB( LDAB, * ), X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> VLATBS solves one of the triangular systems
*>
*>    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
*>
*> with scaling to prevent overflow, where A is an upper or lower
*> triangular band matrix.  Here A**T denotes the transpose of A, x and b
*> are n-element vectors, and s is a scaling factor, usually less than
*> or equal to 1, chosen so that the components of x will be less than
*> the overflow threshold.  If the unscaled problem will not cause
*> overflow, the Level 2 BLAS routine VTBSV is called.  If the matrix A
*> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
*> non-trivial solution to A*x = 0 is returned.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the matrix A is upper or lower triangular.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the operation applied to A.
*>          = 'N':  Solve A * x = s*b     (No transpose)
*>          = 'T':  Solve A**T * x = s*b  (Transpose)
*>          = 'C':  Solve A**H * x = s*b  (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          Specifies whether or not the matrix A is unit triangular.
*>          = 'N':  Non-unit triangular
*>          = 'U':  Unit triangular
*> \endverbatim
*>
*> \param[in] NORMIN
*> \verbatim
*>          NORMIN is CHARACTER*1
*>          Specifies whether CNORM has been set or not.
*>          = 'Y':  CNORM contains the column norms on entry
*>          = 'N':  CNORM is not set on entry.  On exit, the norms will
*>                  be computed and stored in CNORM.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] KD
*> \verbatim
*>          KD is INTEGER
*>          The number of subdiagonals or superdiagonals in the
*>          triangular matrix A.  KD >= 0.
*> \endverbatim
*>
*> \param[in] AB
*> \verbatim
*>          AB is TYPE(cmplx64x2) array, dimension (LDAB,N)
*>          The upper or lower triangular band matrix A, stored in the
*>          first KD+1 rows of the array. The j-th column of A is stored
*>          in the j-th column of the array AB as follows:
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= KD+1.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is TYPE(cmplx64x2) array, dimension (N)
*>          On entry, the right hand side b of the triangular system.
*>          On exit, X is overwritten by the solution vector x.
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is TYPE(real64x2)
*>          The scaling factor s for the triangular system
*>             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.
*>          If SCALE = 0, the matrix A is singular or badly scaled, and
*>          the vector x is an exact or approximate solution to A*x = 0.
*> \endverbatim
*>
*> \param[in,out] CNORM
*> \verbatim
*>          CNORM is TYPE(real64x2) array, dimension (N)
*>
*>          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
*>          contains the norm of the off-diagonal part of the j-th column
*>          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
*>          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
*>          must be greater than or equal to the 1-norm.
*>
*>          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
*>          returns the 1-norm of the offdiagonal part of the j-th column
*>          of A.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -k, the k-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup latbs
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  A rough bound on x is computed; if that is less than overflow, VTBSV
*>  is called, otherwise, specific code is used which checks for possible
*>  overflow or divide-by-zero at every operation.
*>
*>  A columnwise scheme is used for solving A*x = b.  The basic algorithm
*>  if A is lower triangular is
*>
*>       x[1:n] := b[1:n]
*>       for j = 1, ..., n
*>            x(j) := x(j) / A(j,j)
*>            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
*>       end
*>
*>  Define bounds on the components of x after j iterations of the loop:
*>     M(j) = bound on x[1:j]
*>     G(j) = bound on x[j+1:n]
*>  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
*>
*>  Then for iteration j+1 we have
*>     M(j+1) <= G(j) / | A(j+1,j+1) |
*>     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
*>            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
*>
*>  where CNORM(j+1) is greater than or equal to the infinity-norm of
*>  column j+1 of A, not counting the diagonal.  Hence
*>
*>     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
*>                  1<=i<=j
*>  and
*>
*>     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
*>                                   1<=i< j
*>
*>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine VTBSV if the
*>  reciprocal of the largest M(j), j=1,..,n, is larger than
*>  max(underflow, 1/overflow).
*>
*>  The bound on x(j) is also used to determine when a step in the
*>  columnwise method can be performed without fear of overflow.  If
*>  the computed bound is greater than a large constant, x is scaled to
*>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
*>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
*>
*>  Similarly, a row-wise scheme is used to solve A**T *x = b  or
*>  A**H *x = b.  The basic algorithm for A upper triangular is
*>
*>       for j = 1, ..., n
*>            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
*>       end
*>
*>  We simultaneously compute two bounds
*>       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
*>       M(j) = bound on x(i), 1<=i<=j
*>
*>  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
*>  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
*>  Then the bound on x(j) is
*>
*>       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
*>
*>            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
*>                      1<=i<=j
*>
*>  and we can safely call VTBSV if 1/M(n) and 1/G(n) are both greater
*>  than max(underflow, 1/overflow).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE VLATBS( UPLO, TRANS, DIAG, NORMIN, N, KD, AB, LDAB,
     $                   X,
     $                   SCALE, CNORM, INFO )
      USE multifloats, only: abs, aimag, cmplx64x2, conjg, max, min,
     +real64x2, dd_half, dd_one, dd_two, dd_zero, operator(+),
     +operator(-), operator(*), operator(/), operator(**), operator(==),
     +operator(/=), operator(<), operator(>), operator(<=),
     +operator(>=), assignment(=)
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, KD, LDAB, N
      TYPE(real64x2)   SCALE
*     ..
*     .. Array Arguments ..
      TYPE(real64x2)   CNORM( * )
      TYPE(cmplx64x2)         AB( LDAB, * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST, JLEN, MAIND
      TYPE(real64x2)   BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCFAC,
     $                   XBND, XJ, XMAX
      TYPE(cmplx64x2)         CSUMJ, TJJS, USCAL, ZDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ITAMAX, IVAMAX
      TYPE(real64x2)   TLAMCH, TVASUM
      TYPE(cmplx64x2)         VDOTC, VDOTU, VLADIV
      EXTERNAL           LSAME, ITAMAX, IVAMAX, TLAMCH, TVASUM,
     $                   VDOTC,
     $                   VDOTU, VLADIV
*     ..
*     .. External Subroutines ..
      EXTERNAL           TSCAL, XERBLA, VAXPY, VTSCAL, VTBSV
*     ..
*     .. Intrinsic Functions ..

*     ..
*     .. Statement Functions ..
      TYPE(real64x2)   CABS1, CABS2
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( real64x2( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      CABS2( ZDUM ) = ABS( real64x2( ZDUM ) / real64x2(limbs=[2.D0, 
     +0.0_8]) ) +
     $                ABS( AIMAG( ZDUM ) / real64x2(limbs=[2.D0, 
     +0.0_8]) )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Test the input parameters.
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT.
     $         LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( KD.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'VLATBS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      SCALE = DD_ONE
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine machine dependent parameters to control overflow.
*
      SMLNUM = TLAMCH( 'Safe minimum' ) / TLAMCH( 'Precision' )
      BIGNUM = DD_ONE / SMLNUM
*
      IF( LSAME( NORMIN, 'N' ) ) THEN
*
*        Compute the 1-norm of each column, not including the diagonal.
*
         IF( UPPER ) THEN
*
*           A is upper triangular.
*
            DO 10 J = 1, N
               JLEN = MIN( KD, J-1 )
               CNORM( J ) = TVASUM( JLEN, AB( KD+1-JLEN, J ), 1 )
   10       CONTINUE
         ELSE
*
*           A is lower triangular.
*
            DO 20 J = 1, N
               JLEN = MIN( KD, N-J )
               IF( JLEN.GT.0 ) THEN
                  CNORM( J ) = TVASUM( JLEN, AB( 2, J ), 1 )
               ELSE
                  CNORM( J ) = DD_ZERO
               END IF
   20       CONTINUE
         END IF
      END IF
*
*     Scale the column norms by TSCFAC if the maximum element in CNORM is
*     greater than BIGNUM/2.
*
      IMAX = ITAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM*DD_HALF ) THEN
         TSCFAC = DD_ONE
      ELSE
         TSCFAC = DD_HALF / ( SMLNUM*TMAX )
         CALL TSCAL( N, TSCFAC, CNORM, 1 )
      END IF
*
*     Compute a bound on the computed solution vector to see if the
*     Level 2 BLAS routine VTBSV can be used.
*
      XMAX = DD_ZERO
      DO 30 J = 1, N
         XMAX = MAX( XMAX, CABS2( X( J ) ) )
   30 CONTINUE
      XBND = XMAX
      IF( NOTRAN ) THEN
*
*        Compute the growth in A * x = b.
*
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = KD + 1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = 1
         END IF
*
         IF( TSCFAC.NE.DD_ONE ) THEN
            GROW = DD_ZERO
            GO TO 60
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, G(0) = max{x(i), i=1,...,n}.
*
            GROW = DD_HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 40 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 60
*
               TJJS = AB( MAIND, J )
               TJJ = CABS1( TJJS )
*
               IF( TJJ.GE.SMLNUM ) THEN
*
*                 M(j) = G(j-1) / abs(A(j,j))
*
                  XBND = MIN( XBND, MIN( DD_ONE, TJJ )*GROW )
               ELSE
*
*                 M(j) could overflow, set XBND to 0.
*
                  XBND = DD_ZERO
               END IF
*
               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN
*
*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
*
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
*
*                 G(j) could overflow, set GROW to 0.
*
                  GROW = DD_ZERO
               END IF
   40       CONTINUE
            GROW = XBND
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( DD_ONE, DD_HALF / MAX( XBND, SMLNUM ) )
            DO 50 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 60
*
*              G(j) = G(j-1)*( 1 + CNORM(j) )
*
               GROW = GROW*( DD_ONE / ( DD_ONE+CNORM( J ) ) )
   50       CONTINUE
         END IF
   60    CONTINUE
*
      ELSE
*
*        Compute the growth in A**T * x = b  or  A**H * x = b.
*
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
            MAIND = KD + 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
            MAIND = 1
         END IF
*
         IF( TSCFAC.NE.DD_ONE ) THEN
            GROW = DD_ZERO
            GO TO 90
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, M(0) = max{x(i), i=1,...,n}.
*
            GROW = DD_HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 70 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 90
*
*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
*
               XJ = DD_ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
*
               TJJS = AB( MAIND, J )
               TJJ = CABS1( TJJS )
*
               IF( TJJ.GE.SMLNUM ) THEN
*
*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
*
                  IF( XJ.GT.TJJ )
     $               XBND = XBND*( TJJ / XJ )
               ELSE
*
*                 M(j) could overflow, set XBND to 0.
*
                  XBND = DD_ZERO
               END IF
   70       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( DD_ONE, DD_HALF / MAX( XBND, SMLNUM ) )
            DO 80 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 90
*
*              G(j) = ( 1 + CNORM(j) )*G(j-1)
*
               XJ = DD_ONE + CNORM( J )
               GROW = GROW / XJ
   80       CONTINUE
         END IF
   90    CONTINUE
      END IF
*
      IF( ( GROW*TSCFAC ).GT.SMLNUM ) THEN
*
*        Use the Level 2 BLAS solve if the reciprocal of the bound on
*        elements of X is not too small.
*
         CALL VTBSV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, X, 1 )
      ELSE
*
*        Use a Level 1 BLAS solve, scaling intermediate results.
*
         IF( XMAX.GT.BIGNUM*DD_HALF ) THEN
*
*           Scale X so that its components are less than or equal to
*           BIGNUM in absolute value.
*
            SCALE = ( BIGNUM*DD_HALF ) / XMAX
            CALL VTSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         ELSE
            XMAX = XMAX*DD_TWO
         END IF
*
         IF( NOTRAN ) THEN
*
*           Solve A * x = b
*
            DO 120 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
*
               XJ = CABS1( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = AB( MAIND, J )*TSCFAC
               ELSE
                  TJJS = TSCFAC
                  IF( TSCFAC.EQ.DD_ONE )
     $               GO TO 110
               END IF
               TJJ = CABS1( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN
*
*                    abs(A(j,j)) > SMLNUM:
*
                  IF( TJJ.LT.DD_ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by 1/b(j).
*
                        REC = DD_ONE / XJ
                        CALL VTSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = VLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               ELSE IF( TJJ.GT.DD_ZERO ) THEN
*
*                    0 < abs(A(j,j)) <= SMLNUM:
*
                  IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
*                       to avoid overflow when dividing by A(j,j).
*
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.DD_ONE ) THEN
*
*                          Scale by 1/CNORM(j) to avoid overflow when
*                          multiplying x(j) times column j.
*
                        REC = REC / CNORM( J )
                     END IF
                     CALL VTSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = VLADIV( X( J ), TJJS )
                  XJ = CABS1( X( J ) )
               ELSE
*
*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to A*x = 0.
*
                  DO 100 I = 1, N
                     X( I ) = DD_ZERO
  100             CONTINUE
                  X( J ) = DD_ONE
                  XJ = DD_ONE
                  SCALE = DD_ZERO
                  XMAX = DD_ZERO
               END IF
  110          CONTINUE
*
*              Scale x if necessary to avoid overflow when adding a
*              multiple of column j of A.
*
               IF( XJ.GT.DD_ONE ) THEN
                  REC = DD_ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
*
*                    Scale x by 1/(2*abs(x(j))).
*
                     REC = REC*DD_HALF
                     CALL VTSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
*
*                 Scale x by 1/2.
*
                  CALL VTSCAL( N, DD_HALF, X, 1 )
                  SCALE = SCALE*DD_HALF
               END IF
*
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
*
*                    Compute the update
*                       x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
*                                             x(j)* A(max(1,j-kd):j-1,j)
*
                     JLEN = MIN( KD, J-1 )
                     CALL VAXPY( JLEN, -X( J )*TSCFAC,
     $                           AB( KD+1-JLEN, J ), 1, X( J-JLEN ), 1 )
                     I = IVAMAX( J-1, X, 1 )
                     XMAX = CABS1( X( I ) )
                  END IF
               ELSE IF( J.LT.N ) THEN
*
*                 Compute the update
*                    x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
*                                          x(j) * A(j+1:min(j+kd,n),j)
*
                  JLEN = MIN( KD, N-J )
                  IF( JLEN.GT.0 )
     $               CALL VAXPY( JLEN, -X( J )*TSCFAC, AB( 2, J ), 1,
     $                           X( J+1 ), 1 )
                  I = J + IVAMAX( N-J, X( J+1 ), 1 )
                  XMAX = CABS1( X( I ) )
               END IF
  120       CONTINUE
*
         ELSE IF( LSAME( TRANS, 'T' ) ) THEN
*
*           Solve A**T * x = b
*
            DO 170 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               XJ = CABS1( X( J ) )
               USCAL = TSCFAC
               REC = DD_ONE / MAX( XMAX, DD_ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  REC = REC*DD_HALF
                  IF( NOUNIT ) THEN
                     TJJS = AB( MAIND, J )*TSCFAC
                  ELSE
                     TJJS = TSCFAC
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.DD_ONE ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     REC = MIN( DD_ONE, REC*TJJ )
                     USCAL = VLADIV( USCAL, TJJS )
                  END IF
                  IF( REC.LT.DD_ONE ) THEN
                     CALL VTSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
*
               CSUMJ = DD_ZERO
               IF( USCAL.EQ.cmplx64x2( DD_ONE ) ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call VDOTU to perform the dot product.
*
                  IF( UPPER ) THEN
                     JLEN = MIN( KD, J-1 )
                     CSUMJ = VDOTU( JLEN, AB( KD+1-JLEN, J ), 1,
     $                       X( J-JLEN ), 1 )
                  ELSE
                     JLEN = MIN( KD, N-J )
                     IF( JLEN.GT.1 )
     $                  CSUMJ = VDOTU( JLEN, AB( 2, J ), 1, X( J+1 ),
     $                          1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( UPPER ) THEN
                     JLEN = MIN( KD, J-1 )
                     DO 130 I = 1, JLEN
                        CSUMJ = CSUMJ + ( AB( KD+I-JLEN, J )*USCAL )*
     $                          X( J-JLEN-1+I )
  130                CONTINUE
                  ELSE
                     JLEN = MIN( KD, N-J )
                     DO 140 I = 1, JLEN
                        CSUMJ = CSUMJ + ( AB( I+1, J )*USCAL )*X( J+I )
  140                CONTINUE
                  END IF
               END IF
*
               IF( USCAL.EQ.cmplx64x2( TSCFAC ) ) THEN
*
*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  IF( NOUNIT ) THEN
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                     TJJS = AB( MAIND, J )*TSCFAC
                  ELSE
                     TJJS = TSCFAC
                     IF( TSCFAC.EQ.DD_ONE )
     $                  GO TO 160
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( TJJ.LT.DD_ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           REC = DD_ONE / XJ
                           CALL VTSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = VLADIV( X( J ), TJJS )
                  ELSE IF( TJJ.GT.DD_ZERO ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL VTSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = VLADIV( X( J ), TJJS )
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0 and compute a solution to A**T *x = 0.
*
                     DO 150 I = 1, N
                        X( I ) = DD_ZERO
  150                CONTINUE
                     X( J ) = DD_ONE
                     SCALE = DD_ZERO
                     XMAX = DD_ZERO
                  END IF
  160             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
*                 product has already been divided by 1/A(j,j).
*
                  X( J ) = VLADIV( X( J ), TJJS ) - CSUMJ
               END IF
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
  170       CONTINUE
*
         ELSE
*
*           Solve A**H * x = b
*
            DO 220 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               XJ = CABS1( X( J ) )
               USCAL = TSCFAC
               REC = DD_ONE / MAX( XMAX, DD_ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  REC = REC*DD_HALF
                  IF( NOUNIT ) THEN
                     TJJS = CONJG( AB( MAIND, J ) )*TSCFAC
                  ELSE
                     TJJS = TSCFAC
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.DD_ONE ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     REC = MIN( DD_ONE, REC*TJJ )
                     USCAL = VLADIV( USCAL, TJJS )
                  END IF
                  IF( REC.LT.DD_ONE ) THEN
                     CALL VTSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
*
               CSUMJ = DD_ZERO
               IF( USCAL.EQ.cmplx64x2( DD_ONE ) ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call VDOTC to perform the dot product.
*
                  IF( UPPER ) THEN
                     JLEN = MIN( KD, J-1 )
                     CSUMJ = VDOTC( JLEN, AB( KD+1-JLEN, J ), 1,
     $                       X( J-JLEN ), 1 )
                  ELSE
                     JLEN = MIN( KD, N-J )
                     IF( JLEN.GT.1 )
     $                  CSUMJ = VDOTC( JLEN, AB( 2, J ), 1, X( J+1 ),
     $                          1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( UPPER ) THEN
                     JLEN = MIN( KD, J-1 )
                     DO 180 I = 1, JLEN
                        CSUMJ = CSUMJ + ( CONJG( AB( KD+I-JLEN, J ) )*
     $                          USCAL )*X( J-JLEN-1+I )
  180                CONTINUE
                  ELSE
                     JLEN = MIN( KD, N-J )
                     DO 190 I = 1, JLEN
                        CSUMJ = CSUMJ + ( CONJG( AB( I+1, J ) )*USCAL )
     $                          *X( J+I )
  190                CONTINUE
                  END IF
               END IF
*
               IF( USCAL.EQ.cmplx64x2( TSCFAC ) ) THEN
*
*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  X( J ) = X( J ) - CSUMJ
                  XJ = CABS1( X( J ) )
                  IF( NOUNIT ) THEN
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                     TJJS = CONJG( AB( MAIND, J ) )*TSCFAC
                  ELSE
                     TJJS = TSCFAC
                     IF( TSCFAC.EQ.DD_ONE )
     $                  GO TO 210
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( TJJ.LT.DD_ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           REC = DD_ONE / XJ
                           CALL VTSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = VLADIV( X( J ), TJJS )
                  ELSE IF( TJJ.GT.DD_ZERO ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL VTSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = VLADIV( X( J ), TJJS )
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0 and compute a solution to A**H *x = 0.
*
                     DO 200 I = 1, N
                        X( I ) = DD_ZERO
  200                CONTINUE
                     X( J ) = DD_ONE
                     SCALE = DD_ZERO
                     XMAX = DD_ZERO
                  END IF
  210             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
*                 product has already been divided by 1/A(j,j).
*
                  X( J ) = VLADIV( X( J ), TJJS ) - CSUMJ
               END IF
               XMAX = MAX( XMAX, CABS1( X( J ) ) )
  220       CONTINUE
         END IF
         SCALE = SCALE / TSCFAC
      END IF
*
*     Scale the column norms by 1/TSCFAC for return.
*
      IF( TSCFAC.NE.DD_ONE ) THEN
         CALL TSCAL( N, DD_ONE / TSCFAC, CNORM, 1 )
      END IF
*
      RETURN
*
*     End of VLATBS
*
      END
