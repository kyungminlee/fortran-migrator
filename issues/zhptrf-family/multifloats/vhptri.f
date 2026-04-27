*> \brief \b VHPTRI
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download VHPTRI + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/vhptri.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/vhptri.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/vhptri.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE VHPTRI( UPLO, N, AP, IPIV, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       TYPE(cmplx64x2)         AP( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> VHPTRI computes the inverse of a complex Hermitian indefinite matrix
*> A in packed storage using the factorization A = U*D*U**H or
*> A = L*D*L**H computed by VHPTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U*D*U**H;
*>          = 'L':  Lower triangular, form is A = L*D*L**H.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] AP
*> \verbatim
*>          AP is TYPE(cmplx64x2) array, dimension (N*(N+1)/2)
*>          On entry, the block diagonal matrix D and the multipliers
*>          used to obtain the factor U or L as computed by VHPTRF,
*>          stored as a packed triangular matrix.
*>
*>          On exit, if INFO = 0, the (Hermitian) inverse of the original
*>          matrix, stored as a packed triangular matrix. The j-th column
*>          of inv(A) is stored in the array AP as follows:
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;
*>          if UPLO = 'L',
*>             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by VHPTRF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is TYPE(cmplx64x2) array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
*>               inverse could not be computed.
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
*> \ingroup hptri
*
*  =====================================================================
      SUBROUTINE VHPTRI( UPLO, N, AP, IPIV, WORK, INFO )
      USE multifloats, only: abs, cmplx64x2, conjg, real, real64x2,
     +dd_one, operator(+), operator(-), operator(*), operator(/),
     +operator(**), operator(==), operator(/=), operator(<),
     +operator(>), operator(<=), operator(>=), assignment(=)
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      TYPE(cmplx64x2)         AP( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      TYPE(cmplx64x2)         CONE, ZERO
      ! Converted to assignments below: PARAMETER          ( ONE = real64x2(limbs=[1.0D0, 0.0_8]), CONE = cmplx64x2(real64x2(limbs=[1.0D0, 0.0_8]),real64x2(limbs=[0.0D0, 0.0_8])),                    ZERO = cmplx64x2(real64x2(limbs=[0.0D0, 0.0_8]),real64x2(limbs=[0.0D0, 0.0_8])) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, K, KC, KCNEXT, KP, KPC, KSTEP, KX, NPP
      TYPE(real64x2)   AK, AKP1, D, T
      TYPE(cmplx64x2)         AKKP1, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      TYPE(cmplx64x2)         VDOTC
      EXTERNAL           LSAME, VDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, VCOPY, VHPMV, VSWAP
*     ..
*     .. Intrinsic Functions ..

*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      CONE = cmplx64x2(real64x2(limbs=[1.0D0, 0.0_8]),
     +real64x2(limbs=[0.0D0, 0.0_8]))
      ZERO = cmplx64x2(real64x2(limbs=[0.0D0, 0.0_8]),
     +real64x2(limbs=[0.0D0, 0.0_8]))

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'VHPTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check that the diagonal matrix D is nonsingular.
*
      IF( UPPER ) THEN
*
*        Upper triangular storage: examine D from bottom to top
*
         KP = N*( N+1 ) / 2
         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. AP( KP ).EQ.ZERO )
     $         RETURN
            KP = KP - INFO
   10    CONTINUE
      ELSE
*
*        Lower triangular storage: examine D from top to bottom.
*
         KP = 1
         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. AP( KP ).EQ.ZERO )
     $         RETURN
            KP = KP + N - INFO + 1
   20    CONTINUE
      END IF
      INFO = 0
*
      IF( UPPER ) THEN
*
*        Compute inv(A) from the factorization A = U*D*U**H.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
         KC = 1
   30    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N )
     $      GO TO 50
*
         KCNEXT = KC + K
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            AP( KC+K-1 ) = DD_ONE / real64x2( AP( KC+K-1 ) )
*
*           Compute column K of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL VCOPY( K-1, AP( KC ), 1, WORK, 1 )
               CALL VHPMV( UPLO, K-1, -CONE, AP, WORK, 1, ZERO,
     $                     AP( KC ), 1 )
               AP( KC+K-1 ) = AP( KC+K-1 ) -
     $                        REAL( VDOTC( K-1, WORK, 1, AP( KC ),
     $                              1 ) )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( AP( KCNEXT+K-1 ) )
            AK = real64x2( AP( KC+K-1 ) ) / T
            AKP1 = real64x2( AP( KCNEXT+K ) ) / T
            AKKP1 = AP( KCNEXT+K-1 ) / T
            D = T*( AK*AKP1-DD_ONE )
            AP( KC+K-1 ) = AKP1 / D
            AP( KCNEXT+K ) = AK / D
            AP( KCNEXT+K-1 ) = -AKKP1 / D
*
*           Compute columns K and K+1 of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL VCOPY( K-1, AP( KC ), 1, WORK, 1 )
               CALL VHPMV( UPLO, K-1, -CONE, AP, WORK, 1, ZERO,
     $                     AP( KC ), 1 )
               AP( KC+K-1 ) = AP( KC+K-1 ) -
     $                        REAL( VDOTC( K-1, WORK, 1, AP( KC ),
     $                              1 ) )
               AP( KCNEXT+K-1 ) = AP( KCNEXT+K-1 ) -
     $                            VDOTC( K-1, AP( KC ), 1,
     $                                   AP( KCNEXT ),
     $                            1 )
               CALL VCOPY( K-1, AP( KCNEXT ), 1, WORK, 1 )
               CALL VHPMV( UPLO, K-1, -CONE, AP, WORK, 1, ZERO,
     $                     AP( KCNEXT ), 1 )
               AP( KCNEXT+K ) = AP( KCNEXT+K ) -
     $                          REAL( VDOTC( K-1, WORK, 1,
     $                                AP( KCNEXT ),
     $                          1 ) )
            END IF
            KSTEP = 2
            KCNEXT = KCNEXT + K + 1
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the leading
*           submatrix A(1:k+1,1:k+1)
*
            KPC = ( KP-1 )*KP / 2 + 1
            CALL VSWAP( KP-1, AP( KC ), 1, AP( KPC ), 1 )
            KX = KPC + KP - 1
            DO 40 J = KP + 1, K - 1
               KX = KX + J - 1
               TEMP = CONJG( AP( KC+J-1 ) )
               AP( KC+J-1 ) = CONJG( AP( KX ) )
               AP( KX ) = TEMP
   40       CONTINUE
            AP( KC+KP-1 ) = CONJG( AP( KC+KP-1 ) )
            TEMP = AP( KC+K-1 )
            AP( KC+K-1 ) = AP( KPC+KP-1 )
            AP( KPC+KP-1 ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = AP( KC+K+K-1 )
               AP( KC+K+K-1 ) = AP( KC+K+KP-1 )
               AP( KC+K+KP-1 ) = TEMP
            END IF
         END IF
*
         K = K + KSTEP
         KC = KCNEXT
         GO TO 30
   50    CONTINUE
*
      ELSE
*
*        Compute inv(A) from the factorization A = L*D*L**H.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         NPP = N*( N+1 ) / 2
         K = N
         KC = NPP
   60    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 )
     $      GO TO 80
*
         KCNEXT = KC - ( N-K+2 )
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            AP( KC ) = DD_ONE / real64x2( AP( KC ) )
*
*           Compute column K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL VCOPY( N-K, AP( KC+1 ), 1, WORK, 1 )
               CALL VHPMV( UPLO, N-K, -CONE, AP( KC+N-K+1 ), WORK, 1,
     $                     ZERO, AP( KC+1 ), 1 )
               AP( KC ) = AP( KC ) - REAL( VDOTC( N-K, WORK, 1,
     $                    AP( KC+1 ), 1 ) )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( AP( KCNEXT+1 ) )
            AK = real64x2( AP( KCNEXT ) ) / T
            AKP1 = real64x2( AP( KC ) ) / T
            AKKP1 = AP( KCNEXT+1 ) / T
            D = T*( AK*AKP1-DD_ONE )
            AP( KCNEXT ) = AKP1 / D
            AP( KC ) = AK / D
            AP( KCNEXT+1 ) = -AKKP1 / D
*
*           Compute columns K-1 and K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL VCOPY( N-K, AP( KC+1 ), 1, WORK, 1 )
               CALL VHPMV( UPLO, N-K, -CONE, AP( KC+( N-K+1 ) ),
     $                     WORK,
     $                     1, ZERO, AP( KC+1 ), 1 )
               AP( KC ) = AP( KC ) - REAL( VDOTC( N-K, WORK, 1,
     $                    AP( KC+1 ), 1 ) )
               AP( KCNEXT+1 ) = AP( KCNEXT+1 ) -
     $                          VDOTC( N-K, AP( KC+1 ), 1,
     $                          AP( KCNEXT+2 ), 1 )
               CALL VCOPY( N-K, AP( KCNEXT+2 ), 1, WORK, 1 )
               CALL VHPMV( UPLO, N-K, -CONE, AP( KC+( N-K+1 ) ),
     $                     WORK,
     $                     1, ZERO, AP( KCNEXT+2 ), 1 )
               AP( KCNEXT ) = AP( KCNEXT ) -
     $                        REAL( VDOTC( N-K, WORK, 1,
     $                              AP( KCNEXT+2 ),
     $                        1 ) )
            END IF
            KSTEP = 2
            KCNEXT = KCNEXT - ( N-K+3 )
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the trailing
*           submatrix A(k-1:n,k-1:n)
*
            KPC = NPP - ( N-KP+1 )*( N-KP+2 ) / 2 + 1
            IF( KP.LT.N )
     $         CALL VSWAP( N-KP, AP( KC+KP-K+1 ), 1, AP( KPC+1 ), 1 )
            KX = KC + KP - K
            DO 70 J = K + 1, KP - 1
               KX = KX + N - J + 1
               TEMP = CONJG( AP( KC+J-K ) )
               AP( KC+J-K ) = CONJG( AP( KX ) )
               AP( KX ) = TEMP
   70       CONTINUE
            AP( KC+KP-K ) = CONJG( AP( KC+KP-K ) )
            TEMP = AP( KC )
            AP( KC ) = AP( KPC )
            AP( KPC ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = AP( KC-N+K-1 )
               AP( KC-N+K-1 ) = AP( KC-N+KP-1 )
               AP( KC-N+KP-1 ) = TEMP
            END IF
         END IF
*
         K = K - KSTEP
         KC = KCNEXT
         GO TO 60
   80    CONTINUE
      END IF
*
      RETURN
*
*     End of VHPTRI
*
      END
