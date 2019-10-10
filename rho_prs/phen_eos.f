       
       program phen_eos
       implicit real*8(a-h,o-z)
        
       dimension :: den(100),xkf(100)
       
       open(404,file='dens.don') 
       open(505,file='phen_pars.don')
       open(999,file='e0_phens.don')
       
       read(505,*) n, nswitch
       
       pi = 3.14159265
       pi2 = pi*pi 

       if(nswitch .EQ. 0) then              
          do i=1,n
             read(404,*) den(i)
          end do
       else if(nswitch .EQ. 1) then
          do i=1,n
             read(404,*) xkf(i)
             den(i) = 2.d0*xkf(i)*xkf(i)*xkf(i)/(3.d0*pi2)
          end do
       end if
       
c      parameatrization of empirical EoS for symmmetric nuclear matter
        
       xk0 = 220
       xk1 = 260

       rho0 = 0.16d0
       
       fact=(3.d0*pi2/2.d0)**(2.d0/3.d0) 
       hbc=197.327d0
       hbc2=hbc**2
       xm=938.926d0 
       tfact=(3.d0*hbc2/10.d0/xm)
       totfact=fact*tfact
c      
       alpha=-29.47-46.74*(xk0+44.21)/(xk0-166.11)
       beta=23.37*(xk0+254.53)/(xk0-166.11)
       sigma=(xk0+44.21)/210.32
       
       alpha1=-29.47-46.74*(xk1+44.21)/(xk1-166.11)
       beta1=23.37*(xk1+254.53)/(xk1-166.11)
       sigma1=(xk1+44.21)/210.32
       
        
       gam=0.72d0
       alph=0.2d0 
       a1=119.14d0
       b1=-816.95d0
       c1=724.51d0
       d1=-32.99d0
       d2=891.15d0
       ff1=a1*2.d0*(0.5d0)**(5.d0/3.d0)
       ff2=d1*2.d0*(0.5d0)**(5.d0/3.d0) + 
     1     d2*2.d0*(0.5d0)**(8.d0/3.d0) 
       
c         phenom_eos_section
         
       do i=1,n 
          datapt = den(i)
          rat=datapt/rho0
          ee=ff1*(datapt)**(2.d0/3.d0) + b1*datapt +
     1    c1*(datapt)**(alph+1.d0) + ff2*(datapt)**(5.d0/3.d0)  
          ee2=totfact*(datapt)**(2.d0/3.d0)+(alpha/2.d0)*(rat)+
     1    (beta/(sigma + 1.d0))*(rat)**(sigma)                
          ee3=totfact*(datapt)**(2.d0/3.d0) + (alpha1/2.d0)*(rat)+
     1    (beta1/(sigma1 + 1.d0))*(rat)**(sigma1)
          pt=22.d0*rat**gam + 12.d0*rat**(2.d0/3.d0)                
          write(999,1010) den(i), ee, ee2, ee3, pt
       end do      
1010   format(2x,F7.3,2x,F8.4,2x,F8.4,2x,F8.4,2x,F8.4)
       end
       
