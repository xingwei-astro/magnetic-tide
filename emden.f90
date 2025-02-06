      program main
      implicit none
      double precision pi
      parameter (pi=acos(-1.d0))
      integer n, i, m
      parameter (n=10000)
      double precision x(n), y(n)
      double precision dx, p

      p=1.5d0
      dx=1.d-3
      do i=1, n
       x(i)=dx*(i-1)
      enddo
      y(1)=1.d0
      y(2)=y(1)
      do i=2, n-1
       y(i+1)=((dx-x(i))*y(i-1)+2.d0*x(i)*y(i)-dx**2*x(i)*y(i)**p)/(x(i)+dx)
       if(y(i+1).lt.0.d0) then
        m=i
        exit
       endif
      enddo
      do i=2, m
       write(1,'(3E15.6)') x(i)/x(m), y(i), sin(pi*x(i)/x(m))/(pi*x(i)/x(m))
      enddo
      end program main

