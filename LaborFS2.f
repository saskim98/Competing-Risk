	program main
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3
      integer, parameter ::  initial = 1000, initer = 100000
      integer, parameter :: niter =  100000, niter1 = 20000 

	real*8 yold, told, dtime, age, labor, labor1, labor2

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc), sigmay

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 sigmas

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb), omega(nb,nb) 
	
	real*8 sigmab, sigmat, sigmad, sigmap, sigmam
	real*8 a0, b0, a1, b1, d0, s0
	
	real*8 mean_betaa(nc)
	real*8 mean_sigmay
	real*8 mean_theta(nt)
	real*8 mean_sigmas
	real*8 mean_phi(np,nl)  
	real*8 mean_bmu
	real*8 mean_omega(nb,nb)
      
	real*8 bardic, dicbar, dic
      
      common /vij/ij
	
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      common /vsigmay/sigmay
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      common /vsigmas/sigmas

      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega

      common /vsigmab/sigmab
      common /vsigmat/sigmat
      common /vsigmad/sigmad
      common /vsigmap/sigmap
      common /vsigmam/sigmam
      common /va0/a0
      common /vb0/b0
      common /va1/a1
      common /vb1/b1
      common /vd0/d0
      common /vs0/s0
      
      external gen_dic
      
    	open(unit = 5, file = 'Labor1.txt')    	
    	open(unit = 6, file = 'LaborFS2_Initial.txt')
    	    	
    	open(unit = 11, file = 'LaborFS2_Output1.txt')
    	open(unit = 12, file = 'LaborFS2_Output2.txt')
    	open(unit = 13, file = 'LaborFS2_Output3.txt')
    	open(unit = 14, file = 'LaborFS2_Output4.txt')
    	open(unit = 15, file = 'LaborFS2_Output5.txt')
    	open(unit = 16, file = 'LaborFS2_Output6.txt')
    	            
    	open(unit = 21, file = 'LaborFS2_Output7.txt')
      
      iseed = 9999999

      do i = 1, ni
        do j = 1, nj
            yij(i,j) = -999.d0
        enddo
    	enddo
    	  
	do id = 1, 1835
	
	  read(5,*) i, j, yold, told, dtime, 
     +                  age, labor, labor1, labor2

        age = (age - 28.6603139d0)/5.04093971d0
        labor1 = (labor1 - 0.1591928251d0)/0.3662666476d0
        labor2 = (labor2 - 0.1457399103d0)/0.353241509d0
        
        yij(i,j) = yold
        
        tij(i,j) = told
        
        ti(i) = dtime
        
        zi(i) = labor
        	  
	  ui(i,1) = 1.d0
	  ui(i,2) = age
	  ui(i,3) = labor1
	  ui(i,4) = labor2

	  xi(i,1) = 1.d0
	  xi(i,2) = age
        
	enddo    

	do i = 1, ni
	  ij(i) = 0
	  do j = 1, nj	  
            if (yij(i,j) .ne. -999.d0) then
                ij(i) = ij(i) + 1
            endif
        enddo
      enddo  
            	      	      	
c     Hyperparameter for prior distribution

	sigmab = 100.d0; sigmat = 100.d0
      sigmap = 100.d0; sigmam = 100.d0

      d0 = dfloat(nb) + 0.1d0; s0 = 0.1d0     
      a0 = 2.0d0 ; b0 = 1.0d0  
      a1 = 2.0d0 ; b1 = 1.0d0

c     Initial values 

      read(6,*) betaa
      read(6,*) sigmay
      read(6,*) theta
      read(6,*) sigmas
      read(6,*) phi 
      read(6,*) bmu
      read(6,*) omega
                  
      do k = 1, nb
        do i = 1, ni
            bi(i,k) = 0.d0
        enddo
      enddo
                  
      call rnset(iseed)          
            
c     Initial sampling for latent variables
            
      call gibbs_bi1_MLE(iseed)
      call gibbs_bi2_MLE(iseed)
      call gibbs_bi3_MLE(iseed)
      call gibbs_bi4_MLE(iseed)
      do ir = 1, initial
          write(*,*) ir
          call gibbs_bi1(iseed)
          call gibbs_bi2(iseed)
          call gibbs_bi3(iseed)
          call gibbs_bi4(iseed)
	enddo

      
c     Gibbs sampler
            
      icount = 0
	do ir = 1, initer
	
          call gibbs(iseed)	

          write(*,1) ir, betaa, sigmay
          write(*,1) ir, theta
          write(*,1) ir, sigmas
          do l = 1, nl-1
              write(*,1) ir, (phi(j,l),j=1,np)
          enddo
          write(*,1) ir, bmu
          write(*,*)
        
          if (ir - ir/5*5 .eq. 1) then
        
              icount = icount + 1
        
              write(11,2) icount, betaa, sigmay
              write(12,2) icount, theta
              write(13,2) icount, sigmas
              write(14,2) icount, phi 
              write(15,2) icount, bmu
              write(16,2) icount, omega
            
          endif
        
      enddo           

      do jj = 1, nc
          mean_betaa(jj) = 0.d0
      enddo
      
	mean_sigmay = 0.d0

      do jj = 1, nt
          mean_theta(jj) = 0.d0
      enddo
            
	mean_sigmas = 0.d0
      
      do l = 1, nl
          do jj = 1, np
              mean_phi(jj,l) = 0.d0
          enddo
      enddo
      
	mean_bmu = 0.d0

      do j1 = 1, nb
          do j2 = 1, nb
              mean_omega(j1,j2) = 0.d0
          enddo
      enddo
            
      bardic = 0.0d0      
      
      icount = 0
	do ir = 1, niter
	
          call gibbs(iseed)	

          write(*,1) ir, betaa, sigmay
          write(*,1) ir, theta
          write(*,1) ir, sigmas
          do l = 1, nl-1
              write(*,1) ir, (phi(j,l),j=1,np)
          enddo
          write(*,1) ir, bmu
          write(*,*)
        
          if (ir - ir/5*5 .eq. 1) then
        
              icount = icount + 1
        
              write(11,2) icount, betaa, sigmay
              write(12,2) icount, theta
              write(13,2) icount, sigmas
              write(14,2) icount, phi 
              write(15,2) icount, bmu
              write(16,2) icount, omega
            
              do jj = 1, nc
                  mean_betaa(jj) = mean_betaa(jj) 
     +                           + betaa(jj)/dfloat(niter1)
              enddo
      
	        mean_sigmay = mean_sigmay + sigmay/dfloat(niter1)

              do jj = 1, nt
                  mean_theta(jj) = mean_theta(jj) 
     +                           + theta(jj)/dfloat(niter1)
              enddo
            
	        mean_sigmas = mean_sigmas + sigmas/dfloat(niter1)
      
              do l = 1, nl
                  do jj = 1, np
                      mean_phi(jj,l) = mean_phi(jj,l) 
     +                               + phi(jj,l)/dfloat(niter1)
                  enddo
              enddo
      
	        mean_bmu = mean_bmu + bmu/dfloat(niter1)

              do j1 = 1, nb
                  do j2 = 1, nb
                      mean_omega(j1,j2) = mean_omega(j1,j2) 
     +                                  + omega(j1,j2)
     +                                    /dfloat(niter1)
                  enddo
              enddo
              
              
              bardic = bardic 
     +               + gen_dic(betaa, sigmay, theta,  
     +                         sigmas, phi, bmu, omega)
     +                 /dfloat(niter1)
              
          endif
        
      enddo           
              
c     Posterior estimates      
      
      do jj = 1, nc
          betaa(jj) = mean_betaa(jj)
      enddo
      
	sigmay = mean_sigmay

      do jj = 1, nt
          theta(jj) = mean_theta(jj)
      enddo
            
	mean_sigmas = sigmas
      
      do l = 1, nl
          do jj = 1, np
              mean_phi(jj,l) = phi(jj,l)
          enddo
      enddo
      
	mean_bmu = bmu

      do j1 = 1, nb
          do j2 = 1, nb
              mean_omega(j1,j2) = omega(j1,j2)
          enddo
      enddo
                  
      do k = 1, nb
          do i = 1, ni
              bi(i,k) = 0.d0
          enddo
      enddo
                                          
      call gibbs_bi1_MLE(iseed)
      call gibbs_bi2_MLE(iseed)
      call gibbs_bi3_MLE(iseed)
      call gibbs_bi4_MLE(iseed)
      do ir = 1, initial
          write(*,*) ir
          call gibbs_bi1(iseed)
          call gibbs_bi2(iseed)
          call gibbs_bi3(iseed)
          call gibbs_bi4(iseed)
	enddo
      
	do ir = 1, initer
	
          write(*,1) ir
          call gibbs_bi1(iseed)      
          call gibbs_bi2(iseed)  
          call gibbs_bi3(iseed)    
          call gibbs_bi4(iseed)      
        
      enddo           
      
      dicbar = 0.d0
      icount = 0
	do ir = 1, niter
	
          write(*,1) ir
          call gibbs_bi1(iseed) 
          call gibbs_bi2(iseed)  
          call gibbs_bi3(iseed) 
          call gibbs_bi4(iseed)      

          if (ir - ir/5*5 .eq. 1) then
        
              icount = icount + 1
                      
              dicbar = dicbar
     +               + gen_dic(betaa, sigmay, theta, 
     +                         sigmas, phi, bmu, omega)
     +                 /dfloat(niter1)
              
          endif
        
      enddo           
      
      call rnget(iseed)

      dic = -4.d0*bardic + 2.d0*dicbar
      
      write(21,3) dicbar, bardic, dic   
      
    1 format(i6, 100f10.5)   
    2 format(i6, 100f20.10)   
    3 format(100f20.10)   

      end program		

    
      real*8 function gen_dic(betaa, sigmay, theta, 
     +                        sigmas, phi, bmu, omega)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3
      
      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)

	real*8 ti(ni), ui(ni,nu)

	real*8 zi(ni), xi(ni,nx)
      
	real*8 betaa(nc), sigmay
	real*8 theta(nt) 
	real*8 sigmas
	real*8 phi(np,nl) 
	real*8 bmu, omega(nb,nb) 
      
	real*8 bi(ni,nb)
      
	real*8 api, omegai(nb,nb)
      real*8 sij, ystar		
	real*8 si, xxt(nt), xt, xp
      real*8 xxp(np), bstar(nb)	
      real*8 suma, sumb, bAb

      integer ipvt(nb)
      real*8 fac(nb,nb), det1, det2, det

	real*8 pdf1, pdf2, pdf3, pdf4, dic	
      
      common /vij/ij
	
      common /vyij/yij
      common /vtij/tij
      
      common /vti/ti
      common /vui/ui

      common /vzi/zi
      common /vxi/xi
      
      common /vbi/bi
                  
      external dconst, dlinrg, dblinf, dlftrg, dlfdrg
      
	api = dconst('PI')
      
      call dlinrg(nb, omega, nb, omegai, nb)
      
      call dlftrg(nb,omega,nb,fac,nb,ipvt)
	call dlfdrg(nb,fac,nb,ipvt,det1,det2)
	det = det1*10.d0**det2
      
      dic = 0.d0
      do i = 1, ni
      
          pdf1 = 0.d0
          do j = 1, ij(i)
        
              sij = tij(i,j) - bmu - bi(i,nb)
            
              if (sij .ge. 0.d0) then
                
                  ystar = yij(i,j) 
     +                  - (betaa(1) + bi(i,1))*sij   
     +                  - (betaa(2) + bi(i,2))*sij**2
            
              else
            
                  ystar = yij(i,j) 
     +                  - (betaa(3) + bi(i,3))*sij   
            
              endif
            
              pdf1 = pdf1
     +             - dlog(2.d0*api*sigmay)/2.d0 
     +             - ystar**2/(2.d0*sigmay)
            
          enddo
                  
              
          si = ti(i) - bmu - bi(i,nb)

          do j = 1, nu
              xxt(j) = ui(i,j)
          enddo
          do j = 1, nb-1
              xxt(nu+j) = bi(i,j)
          enddo
                
          xt = 0.d0
          do j = 1, nt
              xt = xt + xxt(j)*theta(j)
          enddo

          ystar = si - xt
          
          pdf2 = -dlog(2.d0*api*sigmas)/2.d0 
     +         - ystar**2/(2.d0*sigmas)
          

          do j = 1, nx
              xxp(j) = xi(i,j)
          enddo
          do j = 1, nb-1
              xxp(nx+j) = bi(i,j)
          enddo
          
          suma = 0.d0
          do l = 1, nl-1              
        
              xp = 0.d0
              do j = 1, np
                  xp = xp + xxp(j)*phi(j,l)
              enddo
            
              suma = suma + dexp(xp)
            
              sumb = 0.d0
              if (zi(i) .eq. dfloat(l)) then
        
                  sumb = 0.d0
                  do j = 1, np
                      sumb = sumb + xxp(j)*phi(j,l)
                  enddo
                        
              endif  
              
          enddo
          
          pdf3 = sumb - dlog(1.d0 + suma)
        
          
          do jj = 1, nb
              bstar(jj) = bi(i,jj)
          enddo
          
	    bAb = dblinf(nb, nb, omegai, nb, bstar, bstar)
          
          pdf4 = -dfloat(nb)*dlog(2.d0*api)/2.d0 
     +         - dlog(det)/2.d0 
     +         - bAb/2.d0
          
          dic = dic + pdf1 + pdf2 + pdf3 + pdf4 

      enddo

	gen_dic = dic

      end function
       
	include 'LaborFS2_Gibbs.f'
	include 'optim1.f'
	include 'gilks2.f'
	include 'tnorm.f'
