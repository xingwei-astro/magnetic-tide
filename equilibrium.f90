program main
implicit none
integer i, n, n0    ! n0 for radius
parameter (n=1000000)
double precision pi, G, Rg, Ms, Rs, X_H, X_He, mu_m, g_s, rho_m, rho_c, p_c, tem_c
double precision r(0:n), rho(0:n), p(0:n), q(0:n) ! q polytropic index
double precision x(0:n), a, a0
double precision dr, radius, r0, sigma

pi=acos(-1.0d0)
G=6.67d-8
Rg=8.314d7
Ms=1.99d33
Rs=6.96d10
g_s=G*Ms/Rs**2
rho_m=3.0d0*Ms/(4.0d0*pi*Rs**3)
X_H=0.7d0
X_He=0.3d0
mu_m=1.0d0/(0.5d0+1.5d0*X_H+0.25d0*X_He)

dr=1.d-4    	! integration step
r0=13.45d0     	! interface position by manual test
sigma=1.d-3    	! overshoot thickness = pressure scale height = 0.05 Rs
do i=0,n
 r(i)=dr*dble(i)
! q(i)=4.0d0
 q(i)=2.75d0-1.25d0*tanh((r(i)-r0)/sigma)
enddo

r(0)=1.d-6
rho(0)=1.d0
p(0)=1.d0
do i=0,n-1
 x(i)=r(i)**2*rho(i)
 call integral(0,i,r(0:i),x(0:i),a)
 p(i+1)=p(i)-rho(i)*a*dr/r(i)**2
 rho(i+1)=rho(i)*(p(i+1)/p(i))**(1.d0/(1.d0+1.d0/q(i)))
 if(rho(i+1).lt.1.d-7) then
  n0=i
  radius=r(n0)
  exit
 endif
enddo
write(6,*) 'n0=', n0, 'r_max=', radius, 'interface', r0/radius    ! see if r0/radius=0.7 of Sun

call integral(0,n0,r(0:n0),x(0:n0),a0)
write(6,*) '3rho_c/rho_m=', radius**3/a0, 'p_c/rho_c*R/GM=', radius/a0
rho_c=radius**3/a0*rho_m/3.d0
p_c=Rs**2/radius**2*4.d0*pi*G*rho_c**2
tem_c=mu_m/Rg*p_c/rho_c
write(6,*) 'rho_c=', rho_c, 'p_c=', p_c, 'tem_c=', tem_c

do i=0,n0
 x(i)=r(i)**2*rho(i)
 call integral(0,i,r(0:i),x(0:i),a)
 write(2,'(4E15.6)') r(i)/radius, rho(i), 5.d0/3.d0*p(i)/rho(i), 4.d0*pi*G*Rs*rho_c/radius**3*a/g_s
enddo
end program main

subroutine integral(n1,n2,x,y,a)
implicit none        
integer i, n1, n2
double precision x(n1:n2), y(n1:n2), a
a=0.0
do i=n1, n2-1
 a=a+0.5*(x(i+1)-x(i))*(y(i)+y(i+1))
enddo
end subroutine integral
