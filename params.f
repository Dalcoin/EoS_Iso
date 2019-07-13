       program denx

       implicit real*8 (a-h,o-z)
       dimension der(100), esym(100)
       dimension outpt1(100), outpt2(100), coeff1(4,100), coeff2(4,100)
       dimension outpt3(100), outpt4(100), coeff3(4,100), coeff4(4,100)
       dimension desym(100), de0(100), de1(100), prse0(100), prse1(100)
       dimension P(100), den(100), rho(100), rhog(100)
       dimension xkf(100), e0(100), e1(100), prsesym(100)
       dimension e0d(100), e1d(100)
       dimension xkf0(100), xkf1(100), den0(100), den1(100)
       dimension de0dr(100), ee(100), ee2(100)
       dimension outpt5(100), coeff5(100)

       pi = 3.14159265
       pi2 = pi*pi 

       open(unit=120, file='values.don')
       open(unit=140, file='data.srt')

       read(120,*) n, n_opt,mic,isnm,isym_emp,k0,rho0,fff
       n2=n-1

       den_s = 0.0d0
       den_e = 1.8d0



       if(n_opt .eq. 0) then
          open(unit=510, file='ex.don')
       else if(n_opt .eq. 1) then
          open(unit=500, file='e0.don')
          open(unit=501, file='e1.don')
          open(unit=550, file='et.don') 
       end if

       if(n_opt .eq. 0) then
          do i = 1,n
c             read(510,*) den(i), e0(i), e1(i)
             read(510,*) xkf(i), e0(i), e1(i)
             den(i) = (2.d0*xkf(i)**3)/(3.d0*pi2)
             esym(i) = e1(i) - e0(i)
          end do
       else if(n_opt .eq. 1) then
          do i = 1,n
             read(500,*) xkf0(i), den0(i), e0(i)
             read(501,*) xkf1(i), den1(i), e1(i) 
          end do
       end if     

c  parametrization of empirical EoS for symmmetric nuclear matter

       fact=(3.d0*pi2/2.d0)**(2.d0/3.d0) 
       hbc=197.327d0
       hbc2=hbc**2
       xm=938.926d0 
       tfact=(3.d0*hbc2/10.d0/xm)
       totfact=fact*tfact
c
       alpha=-29.47-46.74*(k0+44.21)/(k0-166.11)
       beta=23.37*(k0+254.53)/(k0-166.11)
       sigma=(k0+44.21)/210.32
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
          ee(i)=ff1*(datapt)**(2.d0/3.d0) + b1*datapt +
     1    c1*(datapt)**(alph+1.d0) + ff2*(datapt)**(5.d0/3.d0) 
          
          ee2(i)=totfact*(datapt)**(2.d0/3.d0) + (alpha/2.d0)*(rat)+
     1    (beta/(sigma + 1.d0))*(rat)**(sigma) 
           
          
          if(mic.eq.1) go to 3355 
          
             if(isnm.eq.1) then
                e0(i)=ee(i) 
             else 
                e0(i)=ee2(i) 
             end if 
          
c             if(datapt.le.0.0019.and.e0.gt.0.d0) then
c                e0(i)=0.d0 
c             end if
               
             pt=22.d0*rat**gam + 12.d0*rat**(2.d0/3.d0)      
             if(isym_emp.eq.0) then
                esym(i)= e1(i)-e0(i) 
             else 
                esym(i)=pt      
             end if      
             go to 3355 
          
3355      continue
       end do

       if(n_opt .eq. 0) then
          call dcsakm(n,den,esym,outpt1,coeff1)
          call dcsakm(n,den,e0,outpt2,coeff2)
          call dcsakm(n,den,e1,outpt3,coeff3)
       else if(n_opt .eq. 1) then
          call dcsakm(n,den0,e0,outpt2,coeff2)
          call dcsakm(n,den1,e1,outpt3,coeff3)
          do i=1,n
             den_step = (den_e - den_s)/n*(i-1.d0) + den_s   
             den(i) = den_step
             e0d(i) = dcsval(den(i),n2,outpt2,coeff2)     
             e1d(i) = dcsval(den(i),n2,outpt3,coeff3)  
             esym(i) = e1d(i) - e0d(i) 
          end do
          call dcsakm(n,den,esym,outpt1,coeff1)
       end if

       do i=1,n
         desym(i) = dcsder(1,den(i),n2,outpt1,coeff1)
         de0(i) = dcsder(1,den(i),n2,outpt2,coeff2)
         de1(i) = dcsder(1,den(i),n2,outpt3,coeff3)
         prse0(i) = den(i)*den(i)*de0(i)   
         prse1(i) = den(i)*den(i)*de1(i)
         prsesym(i) = den(i)*den(i)*desym(i) 
       end do

       call dcsakm(n,de0,den,outpt4,coeff4)
       call dcsakm(n,den,prse0,outpt5,coeff5)
       rho0 = dcsval(0.d0,n2,outpt4,coeff4)
       e0o = dcsval(rho0,n2,outpt2,coeff2)
       rho1 = 0.1d0
       prs0 = dcsval(rho0,n2,outpt5,coeff5)

       esym0 = dcsval(rho0,n2,outpt1,coeff1)
       esym1 = dcsval(rho1,n2,outpt1,coeff1)

       bigL = 3.d0*rho0*dcsder(1,rho0,n2,outpt1,coeff1) 
       bigK = 9.d0*rho0*rho0*dcsder(2,rho0,n2,outpt1,coeff1)
       bigKD = 9.d0*rho1*rho1*dcsder(2,rho1,n2,outpt1,coeff1)
       bigK0 = 9.d0*rho0*rho0*dcsder(2,rho0,n2,outpt2,coeff2)

         
       if(n_opt .eq. 1) then 
          write(550,6734)
          write(550,1000)

          write(550,8629)

          do i=1,n
             write(550,1400) den(i), e0d(i), e1d(i), esym(i)
          end do
          write(550,1000)
          write(550,4545)
          write(550,1000)
          write(550,*) bigL
          write(550,1000)
          write(500,4646)
          write(550,1000)
          write(550,*) bigK


       end if
       
1400   format(3x,F7.3,2x,F7.3,2x,F7.3,2x,F7.3,2x,F7.3,2x,F7.3,2x,F7.3)
1500   format(2x,F8.4,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3
     1        ,2x,F8.3,2x,F8.3)

       
1001   format('-----------------------------------------------|')
1000   format('                                         ')         
6734   format('---------------------------Pressures---------------------
     1-------|')
4545   format('-------------------------L-Value-------------------------
     1-------|')
4646   format('----------------K-Value------------------|')
3636   format('---------------K-0-Value------------------|')
3737   format('---------------K-D-Value------------------|')
4747   format('--------------E-Sat-Value----------------|')
4848   format('--------------R-Sat-Value----------------|')
4949   format('--------------P-Sat-Value----------------|')

5050   format('------------Esym-Sat-Value----------------|')
5151   format('------------Esym-Den-Value----------------|')


8629   format('     Den      P0       P1',       
     1       '       Psym      E0       E1       Esym')   

0001   format('------------Isoscalar-Values-------------|')
1111   format('------------Isovector-Values-------------|')
2222   format('------------Tabulated-Values-------------|')
        
       if(n_opt .eq. 0) then
          write(140,6734)
          write(140,1000)
          write(140,8629)
          do i=1,n
             write(140,1400) den(i), prse0(i), prse1(i), prsesym(i),
     1                               e0(i), e1(i), esym(i) 
          end do


          write(140,1000)
          write(140,0001)
          write(140,1000)

          write(140,1000)
          write(140,4848)
          write(140,1000)
          write(140,*) rho0
          write(140,1000)
          write(140,4747)
          write(140,1000)
          write(140,*) e0o
          write(140,1000)
          write(140,4949) 
          write(140,1000) 
          write(140,*) prs0
          write(140,1000)
          write(140,3636)
          write(140,1000)
          write(140,*) bigK0    
          write(140,1000)      

          write(140,1000)
          write(140,1111)
          write(140,1000)
       
          write(140,1000)        
          write(140,4545)
          write(140,1000)
          write(140,*) bigL       
          write(140,1000)
          write(140,4646)
          write(140,1000)
          write(140,*) bigK
          write(140,1000)
          write(140,3737)
          write(140,1000)
          write(140,*) bigKD
          write(140,1000)
          write(140,5050)
          write(140,1000)
          write(140,*) esym0
          write(140,1000)
          write(140,5151)
          write(140,1000)
          write(140,*) esym1
          write(140,1000) 

          write(140,1000)
          write(140,2222)
          write(140,1000)


          write(140,*) "   rho_o     E_o      K_o        L_o    ",  
     1                 "   K         K_d       E_sym_o   E_sym_den " 
          write(140,1500) rho0, e0o, bigK0, bigL, bigK, bigKD,
     1                    esym0, esym1

       end if
       stop
       
       end program

