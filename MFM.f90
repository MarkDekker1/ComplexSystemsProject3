!----------------------------------------------------------------------
!   MFM model subroutines
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----


      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION X,Y,lda, c1, c2, beta,alpha

	lda = PAR(1)
	c1 = PAR(2)
	c2 = PAR(3)
	beta = PAR(4)
	alpha = PAR(5)

	X = U(1)
	Y = U(2)

	F(1) = alpha - lda * X + X**3/3
	F(2) = beta - (c1+c2*X)*Y+Y**3/3

      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,Z)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: Z

!      Set the parameters
       PAR(1)=1.	!lda
       PAR(2)=1.	!c1
       PAR(3)=1.	!c2
       PAR(4)=0.5	!beta
       PAR(5)=0.	!alpha


!      Set the variables equilibria
       U(1)=0.
       U(2)=0.558
       !U(1)=SQRT(3.)
       !U(2)=0.324

      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
