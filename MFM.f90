!----------------------------------------------------------------------
!   MFM model subroutines
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----


      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION fi,z           ,sigma,N,r,  w

	sigma = PAR(1)
	N = PAR(2)
	r = PAR(3)

	fi = U(1)
	z = U(2)

	w = 1+ (r-1)*z**(N-1)-r/N*(1-z**N)/(1-z)

	F(1) = -w*fi*(1-fi)
	F(2) = (sigma-fi*(r-1))*z*(1-z)*(1-z**(N-1))

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
       PAR(1)=1.	!sigma
       PAR(2)=25.	!N
       PAR(3)=25.	!r


!      Set the variables equilibria
       U(1)=0.23331466895097552
       U(2)=0.0

      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
