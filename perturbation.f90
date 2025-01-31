module globe
implicit none
double complex one
double precision pi, G, kB, m_H
double precision M_star, M_planet, M_moon, R_star, R_planet, R_moon
integer nr, nt, nl, m, l0, m0 ! nt physical space, nl spectral space, m azimuthal mode, l0 and m0 tidal potential wavenumbers
parameter (nr=100, nt=80, nl=60, m=2, l0=2, m0=2)
double precision r(0:nr), theta(nt)
double precision rho0(0:nr), p0(0:nr), T0(0:nr), cs(0:nr), Psi0(0:nr), g0(0:nr)
double precision B01(0:nr,nt), B02(0:nr,nt), B03(0:nr,nt)
double precision rho1_s(0:nr,nl), h1_s(0:nr,nl), Psi1_s(0:nr,nl) ! spectral space
double precision u1_s(0:nr,nl), u2_s(0:nr,nl), u3_s(0:nr,nl) ! spectral space
double precision B1_s(0:nr,nl), B2_s(0:nr,nl), B3_s(0:nr,nl) ! spectral space
double precision rho1_p(0:nr,nt), h1_p(0:nr,nt), Psi1_p(0:nr,nt) ! physical space
double precision u1_p(0:nr,nt), u2_p(0:nr,nt), u3_p(0:nr,nt) ! physical space
double precision B1_p(0:nr,nt), B2_p(0:nr,nt), B3_p(0:nr,nt) ! physical space
double precision gamma_adia(0:nr), N2(0:nr)
double precision vis(0:nr), eta(0:nr)
end module globe

program main
use globe
implicit none
call constant ! give constant
call grid ! give grid
call readin ! read equilibrium rho0, p0, T0, Psi0, g0
call field ! give equilibrium field B01, B02, B03
call diffusivity ! give vis, eta
call linear ! calculate linear perturbation equations
call energy ! calculate kinetic energy, magnetic energy, energy flux
call dissipation ! calculate viscous and ohmic dissipations
call output ! output variables
end program main

subroutine constant
use globe
implicit none
one=(0.d0,1.d0)
pi=acos(-1.d0)
G=6.67d-8
kB=1.38d-16
m_H=1.67d-24
M_star=2.d33
M_planet=2.d30
R_star=7.d10
R_planet=7.d9
end subroutine constant

subroutine grid
use globe
implicit none
integer i, j
end subroutine grid

subroutine readin
use globe
implicit none
integer i
end subroutine readin

subroutine field
use globe
implicit none
integer i, j
end subroutine field

subroutine diffusivity
use globe
implicit none
integer i
end subroutine diffusivity

subroutine linear
use globe
implicit none
integer i, j, l
double precision a(0:nr,nl), b(0:nr,nl), c(0:nr,nl)
end subroutine linear

subroutine energy
use globe
implicit none
integer i, j, l
end subroutine energy

subroutine dissipation
use globe
implicit none
integer i, j, l
end subroutine dissipation

subroutine output
use globe
implicit none
integer i, j, l
end subroutine output 

!     PPP(N,M,C) is the N,M-th associated Legendre function evaluated at
!     C, using the recursion relation below from Abramowitz & Stegun.
      DOUBLE PRECISION FUNCTION PPP(N,M,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION P(0:10000)
      IF ( N .LT. M ) THEN
       PPP=0.D0
       RETURN
      END IF
      S=DSQRT(1.D0 - C**2)
      P(M)=1.D0
      DO 10 I=1,M
       P(M)=P(M) * (-0.5D0) * S * DFLOAT(I+M)
10    CONTINUE
      P(M+1)=DFLOAT(2*M+1) * C * P(M)
      DO 20 I=M+2,N
       P(I)=(DFLOAT(2*I-1)*C*P(I-1) - DFLOAT(I+M-1)*P(I-2))/DFLOAT(I-M)
20    CONTINUE
      PPP=P(N)
      RETURN
      END FUNCTION PPP
