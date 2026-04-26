*> \brief \b TLATRS solves a triangular system of equations with the scale factor set to prevent overflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download TLATRS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/tlatrs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/tlatrs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/tlatrs.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE TLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
*                          CNORM, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, NORMIN, TRANS, UPLO
*       INTEGER            INFO, LDA, N
*       TYPE(real64x2)   SCALE
*       ..
*       .. Array Arguments ..
*       TYPE(real64x2)   A( LDA, * ), CNORM( * ), X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> TLATRS solves one of the triangular systems
*>
*>    A *x = s*b  or  A**T *x = s*b
*>
*> with scaling to prevent overflow.  Here A is an upper or lower
*> triangular matrix, A**T denotes the transpose of A, x and b are
*> n-element vectors, and s is a scaling factor, usually less than
*> or equal to 1, chosen so that the components of x will be less than
*> the overflow threshold.  If the unscaled problem will not cause
*> overflow, the Level 2 BLAS routine TTRSV is called.  If the matrix A
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
*>          = 'N':  Solve A * x = s*b  (No transpose)
*>          = 'T':  Solve A**T* x = s*b  (Transpose)
*>          = 'C':  Solve A**T* x = s*b  (Conjugate transpose = Transpose)
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
*> \param[in] A
*> \verbatim
*>          A is TYPE(real64x2) array, dimension (LDA,N)
*>          The triangular matrix A.  If UPLO = 'U', the leading n by n
*>          upper triangular part of the array A contains the upper
*>          triangular matrix, and the strictly lower triangular part of
*>          A is not referenced.  If UPLO = 'L', the leading n by n lower
*>          triangular part of the array A contains the lower triangular
*>          matrix, and the strictly upper triangular part of A is not
*>          referenced.  If DIAG = 'U', the diagonal elements of A are
*>          also not referenced and are assumed to be 1.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max (1,N).
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is TYPE(real64x2) array, dimension (N)
*>          On entry, the right hand side b of the triangular system.
*>          On exit, X is overwritten by the solution vector x.
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is TYPE(real64x2)
*>          The scaling factor s for the triangular system
*>             A * x = s*b  or  A**T* x = s*b.
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
*> \ingroup latrs
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  A rough bound on x is computed; if that is less than overflow, TTRSV
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
*>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine TTRSV if the
*>  reciprocal of the largest M(j), j=1,..,n, is larger than
*>  max(underflow, 1/overflow).
*>
*>  The bound on x(j) is also used to determine when a step in the
*>  columnwise method can be performed without fear of overflow.  If
*>  the computed bound is greater than a large constant, x is scaled to
*>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
*>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
*>
*>  Similarly, a row-wise scheme is used to solve A**T*x = b.  The basic
*>  algorithm for A upper triangular is
*>
*>       for j = 1, ..., n
*>            x(j) := ( b(j) - A[1:j-1,j]**T * x[1:j-1] ) / A(j,j)
*>       end
*>
*>  We simultaneously compute two bounds
*>       G(j) = bound on ( b(i) - A[1:i-1,i]**T * x[1:i-1] ), 1<=i<=j
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
*>  and we can safely call TTRSV if 1/M(n) and 1/G(n) are both greater
*>  than max(underflow, 1/overflow).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE TLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X,
     $                   SCALE,
     $                   CNORM, INFO )
      USE multifloats, only: abs, max, min, real64x2, dd_half, dd_one,
     +dd_zero, operator(+), operator(-), operator(*), operator(/),
     +operator(**), operator(==), operator(/=), operator(<),
     +operator(>), operator(<=), operator(>=), assignment(=)
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, LDA, N
      TYPE(real64x2)   SCALE
*     ..
*     .. Array Arguments ..
      TYPE(real64x2)   A( LDA, * ), CNORM( * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST
      TYPE(real64x2)   BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS,
     $                   TMAX, TSCFAC, USCAL, XBND, XJ, XMAX
*     ..
*     .. Local Arrays ..
      TYPE(real64x2)   WORK(1)
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ITAMAX
      TYPE(real64x2)   TASUM, TDOT, TLAMCH, TLANGE
      EXTERNAL           LSAME, ITAMAX, TASUM, TDOT, TLAMCH,
     $                   TLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           TAXPY, TSCAL, TTRSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
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
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TLATRS', -INFO )
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
               CNORM( J ) = TASUM( J-1, A( 1, J ), 1 )
   10       CONTINUE
         ELSE
*
*           A is lower triangular.
*
            DO 20 J = 1, N - 1
               CNORM( J ) = TASUM( N-J, A( J+1, J ), 1 )
   20       CONTINUE
            CNORM( N ) = DD_ZERO
         END IF
      END IF
*
*     Scale the column norms by TSCFAC if the maximum element in CNORM is
*     greater than BIGNUM.
*
      IMAX = ITAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM ) THEN
         TSCFAC = DD_ONE
      ELSE
*
*        Avoid NaN generation if entries in CNORM exceed the
*        overflow threshold
*
         IF( TMAX.LE.TLAMCH('Overflow') ) THEN
*           Case 1: All entries in CNORM are valid floating-point numbers
            TSCFAC = DD_ONE / ( SMLNUM*TMAX )
            CALL TSCAL( N, TSCFAC, CNORM, 1 )
         ELSE
*           Case 2: At least one column norm of A cannot be represented
*           as floating-point number. Find the offdiagonal entry A( I, J )
*           with the largest absolute value. If this entry is not +/- Infinity,
*           use this value as TSCFAC.
            TMAX = DD_ZERO
            IF( UPPER ) THEN
*
*              A is upper triangular.
*
               DO J = 2, N
                  TMAX = MAX( TLANGE( 'M', J-1, 1, A( 1, J ), 1,
     $                        WORK ), TMAX )
               END DO
            ELSE
*
*              A is lower triangular.
*
               DO J = 1, N - 1
                  TMAX = MAX( TLANGE( 'M', N-J, 1, A( J+1, J ), 1,
     $                        WORK ), TMAX )
               END DO
            END IF
*
            IF( TMAX.LE.TLAMCH('Overflow') ) THEN
               TSCFAC = DD_ONE / ( SMLNUM*TMAX )
               DO J = 1, N
                  IF( CNORM( J ).LE.TLAMCH('Overflow') ) THEN
                     CNORM( J ) = CNORM( J )*TSCFAC
                  ELSE
*                    Recompute the 1-norm without introducing Infinity
*                    in the summation
                     CNORM( J ) = DD_ZERO
                     IF( UPPER ) THEN
                        DO I = 1, J - 1
                           CNORM( J ) = CNORM( J ) +
     $                                  TSCFAC * ABS( A( I, J ) )
                        END DO
                     ELSE
                        DO I = J + 1, N
                           CNORM( J ) = CNORM( J ) +
     $                                  TSCFAC * ABS( A( I, J ) )
                        END DO
                     END IF
                  END IF
               END DO
            ELSE
*              At least one entry of A is not a valid floating-point entry.
*              Rely on TRSV to propagate Inf and NaN.
               CALL TTRSV( UPLO, TRANS, DIAG, N, A, LDA, X, 1 )
               RETURN
            END IF
         END IF
      END IF
*
*     Compute a bound on the computed solution vector to see if the
*     Level 2 BLAS routine TTRSV can be used.
*
      J = ITAMAX( N, X, 1 )
      XMAX = ABS( X( J ) )
      XBND = XMAX
      IF( NOTRAN ) THEN
*
*        Compute the growth in A * x = b.
*
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
         END IF
*
         IF( TSCFAC.NE.DD_ONE ) THEN
            GROW = DD_ZERO
            GO TO 50
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, G(0) = max{x(i), i=1,...,n}.
*
            GROW = DD_ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 30 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              M(j) = G(j-1) / abs(A(j,j))
*
               TJJ = ABS( A( J, J ) )
               XBND = MIN( XBND, MIN( DD_ONE, TJJ )*GROW )
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
   30       CONTINUE
            GROW = XBND
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( DD_ONE, DD_ONE / MAX( XBND, SMLNUM ) )
            DO 40 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              G(j) = G(j-1)*( 1 + CNORM(j) )
*
               GROW = GROW*( DD_ONE / ( DD_ONE+CNORM( J ) ) )
   40       CONTINUE
         END IF
   50    CONTINUE
*
      ELSE
*
*        Compute the growth in A**T * x = b.
*
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
         END IF
*
         IF( TSCFAC.NE.DD_ONE ) THEN
            GROW = DD_ZERO
            GO TO 80
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, M(0) = max{x(i), i=1,...,n}.
*
            GROW = DD_ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 60 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
*
               XJ = DD_ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
*
*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
*
               TJJ = ABS( A( J, J ) )
               IF( XJ.GT.TJJ )
     $            XBND = XBND*( TJJ / XJ )
   60       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( DD_ONE, DD_ONE / MAX( XBND, SMLNUM ) )
            DO 70 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = ( 1 + CNORM(j) )*G(j-1)
*
               XJ = DD_ONE + CNORM( J )
               GROW = GROW / XJ
   70       CONTINUE
         END IF
   80    CONTINUE
      END IF
*
      IF( ( GROW*TSCFAC ).GT.SMLNUM ) THEN
*
*        Use the Level 2 BLAS solve if the reciprocal of the bound on
*        elements of X is not too small.
*
         CALL TTRSV( UPLO, TRANS, DIAG, N, A, LDA, X, 1 )
      ELSE
*
*        Use a Level 1 BLAS solve, scaling intermediate results.
*
         IF( XMAX.GT.BIGNUM ) THEN
*
*           Scale X so that its components are less than or equal to
*           BIGNUM in absolute value.
*
            SCALE = BIGNUM / XMAX
            CALL TSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         END IF
*
         IF( NOTRAN ) THEN
*
*           Solve A * x = b
*
            DO 110 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
*
               XJ = ABS( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = A( J, J )*TSCFAC
               ELSE
                  TJJS = TSCFAC
                  IF( TSCFAC.EQ.DD_ONE )
     $               GO TO 100
               END IF
               TJJ = ABS( TJJS )
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
                        CALL TSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
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
                     CALL TSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE
*
*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to A*x = 0.
*
                  DO 90 I = 1, N
                     X( I ) = DD_ZERO
   90             CONTINUE
                  X( J ) = DD_ONE
                  XJ = DD_ONE
                  SCALE = DD_ZERO
                  XMAX = DD_ZERO
               END IF
  100          CONTINUE
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
                     CALL TSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
*
*                 Scale x by 1/2.
*
                  CALL TSCAL( N, DD_HALF, X, 1 )
                  SCALE = SCALE*DD_HALF
               END IF
*
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
*
*                    Compute the update
*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
*
                     CALL TAXPY( J-1, -X( J )*TSCFAC, A( 1, J ), 1, X,
     $                           1 )
                     I = ITAMAX( J-1, X, 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               ELSE
                  IF( J.LT.N ) THEN
*
*                    Compute the update
*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
*
                     CALL TAXPY( N-J, -X( J )*TSCFAC, A( J+1, J ), 1,
     $                           X( J+1 ), 1 )
                     I = J + ITAMAX( N-J, X( J+1 ), 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               END IF
  110       CONTINUE
*
         ELSE
*
*           Solve A**T * x = b
*
            DO 160 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               XJ = ABS( X( J ) )
               USCAL = TSCFAC
               REC = DD_ONE / MAX( XMAX, DD_ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  REC = REC*DD_HALF
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCFAC
                  ELSE
                     TJJS = TSCFAC
                  END IF
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.DD_ONE ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     REC = MIN( DD_ONE, REC*TJJ )
                     USCAL = USCAL / TJJS
                  END IF
                  IF( REC.LT.DD_ONE ) THEN
                     CALL TSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
*
               SUMJ = DD_ZERO
               IF( USCAL.EQ.DD_ONE ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call TDOT to perform the dot product.
*
                  IF( UPPER ) THEN
                     SUMJ = TDOT( J-1, A( 1, J ), 1, X, 1 )
                  ELSE IF( J.LT.N ) THEN
                     SUMJ = TDOT( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( UPPER ) THEN
                     DO 120 I = 1, J - 1
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
  120                CONTINUE
                  ELSE IF( J.LT.N ) THEN
                     DO 130 I = J + 1, N
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
  130                CONTINUE
                  END IF
               END IF
*
               IF( USCAL.EQ.TSCFAC ) THEN
*
*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  X( J ) = X( J ) - SUMJ
                  XJ = ABS( X( J ) )
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCFAC
                  ELSE
                     TJJS = TSCFAC
                     IF( TSCFAC.EQ.DD_ONE )
     $                  GO TO 150
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  TJJ = ABS( TJJS )
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
                           CALL TSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE IF( TJJ.GT.DD_ZERO ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL TSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0, and compute a solution to A**T*x = 0.
*
                     DO 140 I = 1, N
                        X( I ) = DD_ZERO
  140                CONTINUE
                     X( J ) = DD_ONE
                     SCALE = DD_ZERO
                     XMAX = DD_ZERO
                  END IF
  150             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot
*                 product has already been divided by 1/A(j,j).
*
                  X( J ) = X( J ) / TJJS - SUMJ
               END IF
               XMAX = MAX( XMAX, ABS( X( J ) ) )
  160       CONTINUE
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
*     End of TLATRS
*
      END
