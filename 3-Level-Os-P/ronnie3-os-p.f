C.....Krotov-Ronnie method for Three-level system.............*
C.....Without using RWA.......................................*
C.....Off-Resonance case and penalty on 2nd state population..*
C
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(n=3,nt=16*1024)
C
      COMMON/block1/ Tc, Tau, Ttime
      COMMON/block2/ alpha, beta
      COMMON/block3/ hbar, mu1, mu2
      COMMON/block4/ W1, W2, W3, Wp, Ws
      COMMON/block5/ Ep, Es
      COMMON/block6/ Epi, Esi
      COMMON/block7/ a2tmp
C
      DOUBLE PRECISION Tc, Tau, Ttime
      DOUBLE PRECISION alpha, beta
      DOUBLE PRECISION hbar, mu1, mu2
      DOUBLE PRECISION W1, W2, W3, Wp, Ws
      DOUBLE PRECISION, DIMENSION(0:nt):: Ep, Es
      DOUBLE PRECISION, DIMENSION(0:nt):: Epi, Esi
      COMPLEX*16 a2tmp(0:nt)
C
      INTEGER i, idt, itr
      DOUBLE PRECISION t, dt
      DOUBLE PRECISION sum1, sum2, TotI, TotF
      DOUBLE PRECISION F11, F22, F33
C
      COMPLEX*16 ai(n), af(n), X1
      COMPLEX*16 a(n), b(n), a1(n), b1(n)
      COMPLEX*16 ait(n,0:nt), aft(n,0:nt)
C
      DOUBLE PRECISION Jo, Jp, Jc
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: Jfi
C
      EXTERNAL derivs
C
      Ttime=1000.D0
      Tc=Ttime/2.D0
      Tau=125.D0
C
      dt=Ttime/DFLOAT(nt)
      PRINT*,"DT=", dt, "NT=", nt
C
      alpha=1.0D03
      beta=0.01D0
C
      hbar=1.D0
      mu1=1.D0
      mu2=1.D0
C
      W1=0.D0
      W2=0.517D0 
      W3=0.517D0-0.375D0
C
      Wp=0.517D0-0.1D0
      Ws=0.375D0-0.1D0
C
C.....Preparation of wave packets
      sum1=0.D0
      sum2=0.D0
      DO i=1,n
         READ(1,*) ai(i)
         READ(2,*) af(i)
         sum1=sum1+(CDABS(ai(i)))**2
         sum2=sum2+(CDABS(af(i)))**2
      END DO
      WRITE(*,*) sum1, sum2
C
C.....Makes the normalization
      DO i=1,n
         ai(i)=ai(i)/DSQRT(sum1)
         af(i)=af(i)/DSQRT(sum2)
      END DO
C
C.....Normalization constant varification
      sum1=0.D0
      sum2=0.D0
      DO i=1,n
         sum1=sum1+(CDABS(ai(i)))**2
         sum2=sum2+(CDABS(af(i)))**2
      END DO
      WRITE(*,*) sum1, sum2
C
      OPEN(unit=21, file='optfield.dat')
      OPEN(unit=22, file='progress.dat')
      OPEN(unit=23, file='population.dat')
C
C.....Optimal Control Iterations
C
      DO i=1,n
         a(i)=ai(i)
      END DO 
C
C.....Forward propagation of Initial wave packet (a)
      DO i=0,nt
         t=DFLOAT(i)*dt
C
         DO i1=1,n
            ait(i1,i)=a(i1)
         END DO
         a2tmp(i)=ait(2,i)
C
         idt=1
C
         CALL RTSV(i,t)
         CALL DERIVS(idt,i,t,a,b)
         CALL RK4(idt,i,a,b,n,t,dt,a,derivs)
C
         F11=0.D0
         F22=0.D0
         F33=0.D0
         F11=CDABS(a(1))**2
         F22=CDABS(a(2))**2
         F33=CDABS(a(3))**2
C
      END DO
      TotI=0.D0
      TotI=TotI+F11+F22+F33
      PRINT*, "TotI=", TotI
C
      ALLOCATE( Jfi(15005) )
      DO itr=1,15000
C
C.....Calculation of the cost functional
C.....Jo is the objective term
      Jo=0.D0
      DO i=1,n
         Jo=Jo+CDABS(DCONJG(ait(i,nt))*af(i))
      END DO
      Jo=Jo**2
C
      Jp=0.D0
      Jc=0.D0
      DO i=0,nt
         t=DFLOAT(i)*dt
C
         Ep0=0.D0
         Es0=0.D0
         Ep0=DEXP(-(t-Tc)**2/(2.D0*Tau**2))
         Es0=DEXP(-(t-Tc)**2/(2.D0*Tau**2))
C
         Jp=Jp + ( (dt/Ep0)*(Ep(i)-Epi(i))**2 + 
     &             (dt/Es0)*(Es(i)-Esi(i))**2 )
         Jc=Jc + dt*(DCONJG(a2tmp(i))*a2tmp(i))
         IF (itr==15000) WRITE(21,*) t, Ep(i), Es(i)
C
      END DO
      Jp=alpha*Jp
      Jc=beta*Jc
C
C.....Jfi is the final optimized cost functional.
      Jfi(itr)=Jo-Jp-Jc
      WRITE(22,*) itr, Jo, Jp, Jc, Jfi(itr)
C
C.....The generation of Lagrange multiplier (a1)
C
      X1=DCMPLX(0.D0,0.D0)
      DO i=1,n
         X1=X1+DCONJG(af(i))*ait(i,nt)
      END DO
C
      DO i=1,n
         a1(i)=DCMPLX(0.D0,0.D0)
      END DO 
C
      DO i=1,n
         a1(i)=af(i)*X1
      END DO
C
C.....Backward propagation of Lagrange multiplier (a1)
      DO i=nt,0,-1
         t=DFLOAT(i)*dt
C
         DO i1=1,n
            aft(i1,i)=a1(i1)
         END DO
C
         idt=-1
C
         CALL DERIVS(idt,i,t,a1,b1)
         CALL RK4(idt,i,a1,b1,n,t,-dt,a1,derivs)
C
         F11=0.D0
         F22=0.D0
         F33=0.D0
         F11=CDABS(a1(1))**2
         F22=CDABS(a1(2))**2
         F33=CDABS(a1(3))**2
C
      END DO
      TotF=0.D0
      TotF=TotF+F11+F22+F33
      PRINT*, "itr=", itr, "TotF=", TotF
C
C.....Forward propagation of Lagrange multiplier (a1)
C      DO i=0,nt
C         t=DFLOAT(i)*dt
CC
C         DO i1=1,n
C            aft(i1,i)=a1(i1)
C         END DO
CC
          idt=-1
CC
C         CALL DERIVS(idt,i,t,a1,b1)
C         CALL RK4(idt,i,a1,b1,n,t,dt,a1,derivs)
CC
C      END DO
C
C.....Forward propagation of Initial wave packet (a)
      DO i=1,n
         a(i)=ai(i)
      END DO
C
      DO i=0,nt
         t=DFLOAT(i)*dt
C
         DO i1=1,n
            ait(i1,i)=a(i1)
         END DO
         a2tmp(i)=ait(2,i)
C
         idt=1
C
         CALL FIELD(i,t,ait,aft)
         CALL DERIVS(idt,i,t,a,b)
         CALL RK4(idt,i,a,b,n,t,dt,a,derivs)
C
         F11=0.D0
         F22=0.D0
         F33=0.D0
         F11=CDABS(a(1))**2
         F22=CDABS(a(2))**2
         F33=CDABS(a(3))**2
         IF (itr==15000) WRITE(23,7779) t, F11, F22, F33
 7779	format(2x,58F15.7)
	 
C
      END DO
      TotI=0.D0
      TotI=TotI+F11+F22+F33
      PRINT*, "itr=", itr, "TotI=", TotI
C
      END DO
      END
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C.....See Numerical Recipes, p553 (1st edition)
      SUBROUTINE RK4(idt,i2,y,dydx,n,x,h,yout,derivs)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
C
      INTEGER idt, i2, j, n, nmax
      PARAMETER (nmax=10)
C
      DOUBLE PRECISION h, x, h6, hh, xh
      COMPLEX*16 y(n), dydx(n), yout(n)
      COMPLEX*16 yt(nmax), dyt(nmax), dym(nmax)
C
      hh=h*0.5D0
      h6=h/6.D0
      xh=x+hh
C
      DO j=1,n
         yt(j)=y(j)+hh*dydx(j)
      END DO
      CALL DERIVS(idt,i2,xh,yt,dyt)
C
      DO j=1,n
         yt(j)=y(j)+hh*dyt(j)
      END DO
      CALL DERIVS(idt,i2,xh,yt,dym)
C
      DO j=1,n
         yt(j)=y(j)+h*dym(j)
         dym(j)=dyt(j)+dym(j)
      END DO
      CALL DERIVS(idt,i2,x+h,yt,dyt)
C
      DO j=1,n
         yout(j)=y(j)+h6*(dydx(j)+dyt(j)+2.D0*dym(j))
      END DO
C
      END SUBROUTINE RK4
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      SUBROUTINE RTSV(k,t)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(n=3,nt=16*1024)
C
      COMMON/block1/ Tc, Tau, Ttime
      COMMON/block2/ alpha, beta
      COMMON/block3/ hbar, mu1, mu2
      COMMON/block4/ W1, W2, W3, Wp, Ws
      COMMON/block5/ Ep, Es
      COMMON/block6/ Epi, Esi
      COMMON/block7/ a2tmp
C
      DOUBLE PRECISION Tc, Tau, Ttime
      DOUBLE PRECISION alpha, beta
      DOUBLE PRECISION hbar, mu1, mu2
      DOUBLE PRECISION W1, W2, W3, Wp, Ws
      DOUBLE PRECISION, DIMENSION(0:nt):: Ep, Es
      DOUBLE PRECISION, DIMENSION(0:nt):: Epi, Esi
      COMPLEX*16 a2tmp(0:nt)
C
      INTEGER k
      DOUBLE PRECISION t, Epr, Esr, Ep0, Es0
C
      Epr=5.0D-02 !....per fs
      Esr=5.0D-02 !....per fs
C
      Ep0=Epr*DEXP(-(t-Tc)**2/(2.D0*Tau**2))
      Es0=Esr*DEXP(-(t-Tc)**2/(2.D0*Tau**2))
C
      Ep(k)=Ep0
      Es(k)=Es0
C
      OPEN(unit=30, file='initfield.dat')
      WRITE(30,*) t, Ep(k), Es(k)
      Epi(k)=Ep(k)
      Esi(k)=Es(k)
C
      END SUBROUTINE RTSV
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      SUBROUTINE FIELD(m,t,ait,aft)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(n=3,nt=16*1024)
C
      COMMON/block1/ Tc, Tau, Ttime
      COMMON/block2/ alpha, beta
      COMMON/block3/ hbar, mu1, mu2
      COMMON/block4/ W1, W2, W3, Wp, Ws
      COMMON/block5/ Ep, Es
      COMMON/block6/ Epi, Esi
      COMMON/block7/ a2tmp
C
      DOUBLE PRECISION Tc, Tau, Ttime
      DOUBLE PRECISION alpha, beta
      DOUBLE PRECISION hbar, mu1, mu2
      DOUBLE PRECISION W1, W2, W3, Wp, Ws
      DOUBLE PRECISION, DIMENSION(0:nt):: Ep, Es
      DOUBLE PRECISION, DIMENSION(0:nt):: Epi, Esi
      COMPLEX*16 a2tmp(0:nt)
C
      INTEGER m
      DOUBLE PRECISION t, Ep0, Es0
      COMPLEX*16 ait(n,0:nt), aft(n,0:nt)
      COMPLEX*16 za1, za2, czi
      DATA czi/(0.D0,1.D0)/
C
      Ep0=0.D0
      Es0=0.D0
      Ep0=DEXP(-(t-Tc)**2/(2.D0*Tau**2))
      Es0=DEXP(-(t-Tc)**2/(2.D0*Tau**2))
C
      za1=DCMPLX(0.D0,0.D0)
      za1=( DCONJG(aft(1,m))*ait(2,m) + DCONJG(aft(2,m))*ait(1,m) )*
     &(1.D0+DCOS(2.D0*Wp*t)) - DCMPLX(0.D0,1.D0)*DSIN(2.D0*Wp*t)*(
     &DCONJG(aft(1,m))*ait(2,m) - DCONJG(aft(2,m))*ait(1,m) )
C
      za2=DCMPLX(0.D0,0.D0)
      za2=( DCONJG(aft(2,m))*ait(3,m) + DCONJG(aft(3,m))*ait(2,m) )*
     &(1.D0+DCOS(2.D0*Ws*t)) + DCMPLX(0.D0,1.D0)*DSIN(2.D0*Ws*t)*(
     &DCONJG(aft(2,m))*ait(3,m) - DCONJG(aft(3,m))*ait(2,m) )
C
      Epi(m)=Ep(m)
      Esi(m)=Es(m)
C
      Ep(m)=Epi(m)-(Ep0/(2.D0*alpha))*AIMAG(za1)
      Es(m)=Esi(m)-(Es0/(2.D0*alpha))*AIMAG(za2)
C
      END SUBROUTINE FIELD
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      SUBROUTINE DERIVS(idt,l,x,y,dydx)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(n=3,nt=16*1024)
C
      COMMON/block1/ Tc, Tau, Ttime
      COMMON/block2/ alpha, beta
      COMMON/block3/ hbar, mu1, mu2
      COMMON/block4/ W1, W2, W3, Wp, Ws
      COMMON/block5/ Ep, Es
      COMMON/block6/ Epi, Esi
      COMMON/block7/ a2tmp
C
      DOUBLE PRECISION Tc, Tau, Ttime
      DOUBLE PRECISION alpha, beta
      DOUBLE PRECISION hbar, mu1, mu2
      DOUBLE PRECISION W1, W2, W3, Wp, Ws
      DOUBLE PRECISION, DIMENSION(0:nt):: Ep, Es
      DOUBLE PRECISION, DIMENSION(0:nt):: Epi, Esi
      COMPLEX*16 a2tmp(0:nt)
C
      INTEGER l, idt
      DOUBLE PRECISION x, Delta
      COMPLEX*16 czi, y(n), dydx(n)
      DATA czi/(0.D0,1.D0)/
C
      Delta=0.1D0

      dydx(1)= (czi/2.D0)*Ep(l)*y(2)*(1.D0+DCOS(2.D0*Wp*x) -
     &czi*DSIN(2.D0*Wp*x) )
     &+ czi*Delta*y(1)

      IF (idt==1) THEN
      dydx(2)= (czi/2.D0)* ( Ep(l)*y(1)*(1.D0+DCOS(2.D0*Wp*x) +
     &czi*DSIN(2.D0*Wp*x) ) + Es(l)*y(3)*(1.D0+DCOS(2.D0*Ws*x) +
     &czi*DSIN(2.D0*Ws*x) ) ) 
C    &- czi*Delta*y(2)
      ELSE 
      dydx(2)= (czi/2.D0)* ( Ep(l)*y(1)*(1.D0+DCOS(2.D0*Wp*x) +
     &czi*DSIN(2.D0*Wp*x) ) + Es(l)*y(3)*(1.D0+DCOS(2.D0*Ws*x) +
     &czi*DSIN(2.D0*Ws*x) ) ) + beta*a2tmp(l)
C    &- czi*Delta*y(2)
      END IF

      dydx(3)= (czi/2.D0)*Es(l)*y(2)*(1.D0+DCOS(2.D0*Ws*x) -
     &czi*DSIN(2.D0*Ws*x) ) 
     &+ czi*Delta*y(3)
C
      END SUBROUTINE DERIVS
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
