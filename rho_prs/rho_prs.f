
       program rho_prs
       implicit real*8(a-h,o-z)
       
       dimension outpt1(100), coeff1(4,100)        
       dimension :: e(100),rho(100),den(100)
       dimension :: prs(100)

       open(444,file='par.don')
       open(555,file='ex.don')
       open(999,file='out.don')

 
       read(444,*) n, rho_0, nbeta

       n2=n-1
       rho2_0 = 2.d0*rho_0

       do i=1,n
          read(555,*) rho(i), e(i)
       end do

       call dcsakm(n,rho,e,outpt1,coeff1)       
       
       prho0  = rho_0*rho_0*dcsder(1,rho_0,n2,outpt1,coeff1)
       p2rho0 = rho2_0*rho2_0*dcsder(1,rho2_0,n2,outpt1,coeff1)

       if(nbeta .EQ. 1) then
          prho0  = dcsval(rho_0,n2,outpt1,coeff1)
          p2rho0 = dcsval(rho2_0,n2,outpt1,coeff1)
       end if

       write(999,1010) prho0, p2rho0
1010   format(2x,F8.4,6x,F8.4)
       end

