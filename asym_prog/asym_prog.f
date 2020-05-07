

       program asym_prog
       implicit real*8(a-h,o-z)
       dimension :: den(100), e0(100), e1(100)
       dimension :: esym(100), xkf(100), alp(100)
       dimension :: ea(100)

       
       open(400,file='par.don')
       open(500,file='eos.don')

       open(777,file='parab.don')

       read(400,*) n
       pi=3.14159d0
       pi2 = pi*pi

       do i=1,n
          read(500,*) xkf(i), e0(i), e1(i) 
          esym(i) = e1(i) - e0(i) 
          den(i) = 2.d0*(xkf(i)**3)/(3.d0*pi2)
       end do
       
       do i=1,6 
          alp(i) = 0.d0 + 0.2d0 * (i-1)
       end do

       do i=1,n
          ea(1) = e0(i) + esym(i)*alp(1)**2 
          ea(2) = e0(i) + esym(i)*alp(2)**2
          ea(3) = e0(i) + esym(i)*alp(3)**2
          ea(4) = e0(i) + esym(i)*alp(4)**2
          ea(5) = e0(i) + esym(i)*alp(5)**2
          ea(6) = e0(i) + esym(i)*alp(6)**2
          write(777,1414) den(i), ea(1), ea(2), 
     1                    ea(3), ea(4), ea(5), ea(6)
       end do

1414   format(2x,F8.4,2x,F8.4,2x,F8.4,2x,F8.4,2x,F8.4,
     1        2x,F8.4,2x,F8.4)

       end program  


       function asym_spread(n, xkf, e0, e1)
       implicit real*8(a-h,o-z)
       dimension :: den(100), e0(100), e1(100)
       dimension :: esym(100), xkf(100), alp(100)
       dimension :: ea(100)
        
       pi=3.14159d0
       pi2 = pi*pi

       do i=1,n
          esym(i) = e1(i) - e0(i) 
          den(i) = 2.d0*(xkf(i)**3)/(3.d0*pi2)
       end do
       
       do i=1,6 
          alp(i) = 0.d0 + 0.2d0 * (i-1)
       end do

       do i=1,n
          ea(1) = e0(i) + esym(i)*alp(1)*alp(1) 
          ea(2) = e0(i) + esym(i)*alp(2)*alp(2)
          ea(3) = e0(i) + esym(i)*alp(3)*alp(3)
          ea(4) = e0(i) + esym(i)*alp(4)*alp(4)
          ea(5) = e0(i) + esym(i)*alp(5)*alp(5)
          ea(6) = e0(i) + esym(i)*alp(6)*alp(6)
!          write(777,1414) den(i), ea(1), ea(2), 
!     1                    ea(3), ea(4), ea(5), ea(6)  
       end do

!1414   format(2x,F8.4,2x,F8.4,2x,F8.4,2x,F8.4,2x,F8.4,
!     1        2x,F8.4,2x,F8.4)

       return ea  

