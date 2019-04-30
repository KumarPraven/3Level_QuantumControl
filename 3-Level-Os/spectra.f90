!.......Analyse optimal field into frequency components.
	PROGRAM power_spectrum
        IMPLICIT DOUBLE PRECISION (a-h,o-z)
        INTEGER, PARAMETER :: nt=(64*1024)-1
        DOUBLE PRECISION, DIMENSION(:) :: E(0:nt), Ek(0:nt), EE(0:nt)
        DOUBLE PRECISION :: t, Wp
!
	Totaltime=4*1024.D0
        dt=Totaltime/DFLOAT(nt)
	pi=4.D0*DATAN(1.D0)
!
        OPEN(8,FILE='optfield-o.dat')
        DO i=0,nt
           t=DFLOAT(i)*dt
           READ(8,*) dummy, E(i), EE(i)
        END DO
!        OPEN(7,FILE='power-pump.dat')
!        CALL spectra(E,nt,Totaltime)
        OPEN(7,FILE='power-stoke.dat')
        CALL spectra(EE,nt,Totaltime)
!
	END PROGRAM power_spectrum
!
!========================================================================
!.......This subroutine analyses a time dependent electric field into its
!.......frequency components
!
        SUBROUTINE spectra(E,nt,Totaltime)
        IMPLICIT DOUBLE PRECISION (a-h, o-z)
        INTEGER, INTENT(IN) :: nt
        DOUBLE PRECISION, INTENT(IN) :: Totaltime
        DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: E(0:nt)
        DOUBLE PRECISION, DIMENSION(:) :: FPRINT(0:nt)
        COMPLEX*16, DIMENSION(:) :: FF2(0:nt)
!
!.......atomic unit of time :    t = 2.4189D-17 s
!............................hertz = 4.13412763D+16
!.........................au2hertz = 6.5796839D+15
!
	pi=4.0D0*DATAN(1.0D0)
        au2cm1=219475.797D0
        dt=Totaltime/DFLOAT(nt)
!
!.......Analyse field into frequency components
        DO j=0,nt
           FF2(j)=DCMPLX(E(j))
        END DO
!
        ISGN=+1
!       WRITE(6,*)' Fourier analysis of E(t)'
        CALL FOUR1(FF2,nt+1,ISGN)
!
!.......Scale to get correct f_l after Fourier Transform,
!.......should multiply by a factor of sqrt(nt)*dt.
        factor=DSQRT(DFLOAT(nt))*dt
        FF2(:)=factor*FF2(:)
!
!.......Write out f_l-- weights of field as a frequency spectrum
!.......According to dfreq=1/(N*dt), the unit of frequency is a.u.
!.......dfreq=1.00d00/Totaltime, the unit of frequency is cm inverse.
!
!       dfreq=2.0D0*pi*au2cm1/Totaltime
        dfreq=2.D0*pi/Totaltime
        FPRINT(:)=CDABS(FF2(:))
!
        DO j=0,nt/2
           freq=DFLOAT(j)*dfreq
           WRITE(7,*) freq, FPRINT(j), DREAL(FF2(j)), AIMAG(FF2(j))
        END DO
!        freq=-DFLOAT(nt/2)*dfreq
!        DO j=(nt/2+1),nt
!           WRITE(7,*) freq, FPRINT(j)  !freq printed in cm inverse
!           freq=freq+dfreq
!        END DO
!
!203     FORMAT(D15.8,2X,D15.8)
!
        END SUBROUTINE spectra
!
!========================================================================
!.......Fourier transform pf the data's.
	SUBROUTINE FOUR1(DATA,NN,ISGN)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION DATA(2*NN)
        DATA PI/3.1415926535897932384626433D0/
!
        TPI=2.0D0*PI
        N=2*NN
        C=1.0D0/DSQRT(DFLOAT(NN))
        J=1
        DO 011 I=1,N,2
         IF(J.GT.I)THEN
           TEMPR=DATA(J)
           TEMPI=DATA(J+1)
           DATA(J)=DATA(I)
           DATA(J+1)=DATA(I+1)
           DATA(I)=TEMPR
           DATA(I+1)=TEMPI
         ENDIF
        M=N/2
001     IF((M.GE.2).AND.(J.GT.M))THEN
          J=J-M
          M=M/2
          GO TO 001
        ENDIF
        J=J+M
011     CONTINUE
        MMAX=2
002     IF(N.GT.MMAX)THEN
          ISTEP=2*MMAX
          THETA=TPI/(ISGN*MMAX)
          WPR=-2.0D0*SIN(0.5D0*THETA)**2
          WPI=SIN(THETA)
          WR=1.
          WI=0.
          DO 013 M=1,MMAX,2
           DO 012 I=M,N,ISTEP
              J=I+MMAX
              TEMPR= WR*DATA(J)- WI*DATA(J+1)
              TEMPI= WR*DATA(J+1)+ WI*DATA(J)
              DATA(J)=DATA(I)-TEMPR
              DATA(J+1)=DATA(I+1)-TEMPI
              DATA(I)=DATA(I)+TEMPR
              DATA(I+1)=DATA(I+1)+TEMPI
012        CONTINUE
         WTEMP=WR
         WR=WR*WPR-WI*WPI+WR
         WI=WI*WPR+WTEMP*WPI+WI
013     CONTINUE
        MMAX=ISTEP
        GO TO 002
        ENDIF
        DO 014 I=1,N
           IR=2*I-1
           II=2*1
           DATA(I)=DATA(I)*C
014     CONTINUE
        RETURN
        END SUBROUTINE FOUR1

