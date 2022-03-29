	program main
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 191, nj = 12, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3
      integer, parameter :: initial = 200, initer = 20000
      integer, parameter :: niter = 20000, nniter = 1

	real*8 yold, told, dtime, age
      real*8 labor, labor1, labor2

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc), sigmay

	real*8 ti(ni), tti(ni)
      real*8 ui(ni,nu), theta(nt)      
	real*8 delta, ezeta
      real*8 sigmas, v1

	real*8 zi(ni), zzi(ni)
      real*8 xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

      real*8 pp(ni,nl)
      real*8 error, temp
      
	real*8 tmean(ni)
      real*8 pmean(ni,nl)
      real*8 rmse
      real*8 rmse1
      real*8 rmse2
      real*8 rmse3
      
      real*8 mean_tmean(ni)
      real*8 mean_pmean(ni,nl)
      real*8 mean_rmse
      real*8 mean_rmse1
      real*8 mean_rmse2
      real*8 mean_rmse3

	real*8 std_tmean(ni)
      real*8 std_pmean(ni,nl)
      real*8 std_rmse
      real*8 std_rmse1
      real*8 std_rmse2
      real*8 std_rmse3

	real*8 lower_tmean(ni)
      real*8 lower_pmean(ni,nl)
      real*8 lower_rmse
      real*8 lower_rmse1
      real*8 lower_rmse2
      real*8 lower_rmse3
      
	real*8 upper_tmean(ni)
      real*8 upper_pmean(ni,nl)
      real*8 upper_rmse
      real*8 upper_rmse1
      real*8 upper_rmse2
      real*8 upper_rmse3

	real*8 seq_tmean(niter,ni)
      real*8 seq_pmean(niter,ni,nl)
      real*8 seq_rmse(niter)
      real*8 seq_rmse1(niter)
      real*8 seq_rmse2(niter)
      real*8 seq_rmse3(niter)
      
      real*8 ahpd(niter), alow(2), aupp(2)
      
      common /vij/ij
	
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      common /vsigmay/sigmay
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      common /vdelta/delta
      common /vezeta/ezeta
      common /vsigmas/sigmas
      common /vv1/v1

      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega

      common /vpp/pp
      
    	open(unit = 5, file = 'Labor2.txt')    	
    	open(unit = 6, file = 'LaborFE_Initial.txt')
    	    	
    	open(unit = 11, file = 'LaborFE_Output1.txt')
    	open(unit = 12, file = 'LaborFE_Output2.txt')
    	open(unit = 13, file = 'LaborFE_Output3.txt')
    	open(unit = 14, file = 'LaborFE_Output4.txt')
    	open(unit = 15, file = 'LaborFE_Output5.txt')
    	open(unit = 16, file = 'LaborFE_Output6.txt')

    	open(unit = 21, file = 'PLaborFE_Output1.txt')
    	open(unit = 22, file = 'PLaborFE_Output2.txt')
    	open(unit = 23, file = 'PLaborFE_Output3.txt')
      
      iseed = 9999999

      do i = 1, ni
          do j = 1, nj
              yij(i,j) = -999.d0
          enddo
    	enddo
    	  
	do id = 1, 714
	
	    read(5,*) i, j, yold, told, dtime, 
     +                    age, labor, labor1, labor2

          age = (age - 28.6603139d0)/5.04093971d0
          labor1 = (labor1 - 0.1591928251d0)/0.3662666476d0
          labor2 = (labor2 - 0.1457399103d0)/0.353241509d0
        
          yij(i,j) = yold
        
          tij(i,j) = told
        
          tti(i) = dtime
        
c         labor = 1:vacuum, labor = 2:C-section, labor = 3:spontaneous       
        
          zzi(i) = labor
        	  
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
            	      	      	
      ezeta = 1.d0

c     Initial values 

      read(6,*) betaa
      read(6,*) sigmay
      read(6,*) theta
      read(6,*) delta
      read(6,*) sigmas
      read(6,*) v1
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
      call gibbs_bi12(i,iseed)
      call gibbs_bi3(i,iseed)
      call gibbs_bi4_MLE(i,iseed)
      do ir = 1, initial
          write(*,*) ir
          call gibbs_bi12(i,iseed)
          call gibbs_bi3(i,iseed)
          call gibbs_bi4(i,iseed)
      enddo
      
      
c     Gibbs sampler
            
	do ir = 1, initer
	          
          write(*,*) ir
          read(11,*) irr, betaa, sigmay
          read(12,*) irr, theta
          read(13,*) irr, delta, sigmas, v1
          read(14,*) irr, phi 
          read(15,*) irr, bmu
          read(16,*) irr, omega
                                  
          call gibbs(iseed)	
                                  
      enddo           

      
      do i = 1, ni
          
          mean_tmean(i) = 0.d0
          std_tmean(i) = 0.d0
          
          do l = 1, nl
              mean_pmean(i,l) = 0.d0
              std_pmean(i,l) = 0.d0
          enddo
          
      enddo 
            
      mean_rmse = 0.d0
      std_rmse = 0.d0

      mean_rmse1 = 0.d0
      std_rmse1 = 0.d0

      mean_rmse2 = 0.d0
      std_rmse2 = 0.d0

      mean_rmse3 = 0.d0
      std_rmse3 = 0.d0
      
	do ir = 1, niter
	          
          write(*,*) ir
          read(11,*) irr, betaa, sigmay
          read(12,*) irr, theta
          read(13,*) irr, delta, sigmas, v1
          read(14,*) irr, phi 
          read(15,*) irr, bmu
          read(16,*) irr, omega

          do irr = 1, initial              
              call gibbs(iseed)	
          enddo
          
          do i = 1, ni
              
              tmean(i) = 0.d0
              
              do l = 1, nl
                  pmean(i,l) = 0.d0
              enddo
              
          enddo
          do irr = 1, nniter
              
              call gibbs(iseed)	
              
              do i = 1, ni
                  
                  tmean(i) = tmean(i) + ti(i)/dfloat(nniter)
                  
                  do l = 1, nl
                      pmean(i,l) = pmean(i,l) 
     +                           + pp(i,l)/dfloat(nniter)
                  enddo
                  
              enddo
              
          enddo
          
          do i = 1, ni
              
              mean_tmean(i) = mean_tmean(i) + tmean(i)/dfloat(niter)
              std_tmean(i) = std_tmean(i) + tmean(i)**2
              seq_tmean(ir,i) = tmean(i)

              do l = 1, nl
                  mean_pmean(i,l) = mean_pmean(i,l) 
     +                            + pmean(i,l)/dfloat(niter)
                  std_pmean(i,l) = std_pmean(i,l) + pmean(i,l)**2
                  seq_pmean(ir,i,l) = pmean(i,l)
              enddo
              
          enddo
          
          rmse = 0.d0
          rmse1 = 0.d0; icount1 = 0
          rmse2 = 0.d0; icount2 = 0
          rmse3 = 0.d0; icount3 = 0
          do i = 1, ni
              
              error = tti(i) - tmean(i)
              
              rmse = rmse + error**2/dfloat(ni)
              
              if (zzi(i) .eq. 1.d0) then
                  rmse1 = rmse1 + error**2
                  icount1 = icount1 + 1                  
              endif
                  
              if (zzi(i) .eq. 2.d0) then
                  rmse2 = rmse2 + error**2
                  icount2 = icount2 + 1                  
              endif

              if (zzi(i) .eq. 3.d0) then
                  rmse3 = rmse3 + error**2
                  icount3 = icount3 + 1                  
              endif
              
          enddo
          
          rmse1 = rmse1/dfloat(icount1)
          rmse2 = rmse2/dfloat(icount2)
          rmse3 = rmse3/dfloat(icount3)
          
          mean_rmse = mean_rmse + rmse/dfloat(niter)
          std_rmse = std_rmse + rmse**2
          seq_rmse(ir) = rmse
                    
          mean_rmse1 = mean_rmse1 + rmse1/dfloat(niter)
          std_rmse1 = std_rmse1 + rmse1**2
          seq_rmse1(ir) = rmse1

          mean_rmse2 = mean_rmse2 + rmse2/dfloat(niter)
          std_rmse2 = std_rmse2 + rmse2**2
          seq_rmse2(ir) = rmse2
          
          mean_rmse3 = mean_rmse3 + rmse3/dfloat(niter)
          std_rmse3 = std_rmse3 + rmse3**2
          seq_rmse3(ir) = rmse3
          
      enddo           

      do i = 1, ni
          
		temp = (std_tmean(i) - dfloat(niter)*mean_tmean(i)**2)
     +         / dfloat(niter-1)
		std_tmean(i) = dsqrt(temp)     

          do l = 1, nl

		    temp = (std_pmean(i,l) - dfloat(niter)
     +                                 *mean_pmean(i,l)**2)
     +             / dfloat(niter-1)
		    std_pmean(i,l) = dsqrt(temp)     
              
          enddo
          
      enddo     
      
	temp = (std_rmse - dfloat(niter)*mean_rmse**2)
     +     / dfloat(niter-1)
	std_rmse = dsqrt(temp)
      
	temp = (std_rmse1 - dfloat(niter)*mean_rmse1**2)
     +     / dfloat(niter-1)
	std_rmse1 = dsqrt(temp)

	temp = (std_rmse2 - dfloat(niter)*mean_rmse2**2)
     +     / dfloat(niter-1)
	std_rmse2 = dsqrt(temp)

	temp = (std_rmse3 - dfloat(niter)*mean_rmse3**2)
     +     / dfloat(niter-1)
	std_rmse3 = dsqrt(temp)
      
      do i = 1, ni
          
          do ir = 1, niter
              ahpd(ir) = seq_tmean(ir,i)
          enddo
          call hpd(niter, 0.05d0, ahpd, alow, aupp)
          lower_tmean(i) = alow(1); upper_tmean(i) = aupp(1)
          
          do l = 1, nl

              do ir = 1, niter
                  ahpd(ir) = seq_pmean(ir,i,l)
              enddo
              call hpd(niter, 0.05d0, ahpd, alow, aupp)
              lower_pmean(i,l) = alow(1)
              upper_pmean(i,l) = aupp(1)
                            
          enddo
          
      enddo	
      
      do ir = 1, niter
          ahpd(ir) = seq_rmse(ir)
      enddo
      call hpd(niter, 0.05d0, ahpd, alow, aupp)
      lower_rmse = alow(1); upper_rmse = aupp(1)
      
      do ir = 1, niter
          ahpd(ir) = seq_rmse1(ir)
      enddo
      call hpd(niter, 0.05d0, ahpd, alow, aupp)
      lower_rmse1 = alow(1); upper_rmse1 = aupp(1)

      do ir = 1, niter
          ahpd(ir) = seq_rmse2(ir)
      enddo
      call hpd(niter, 0.05d0, ahpd, alow, aupp)
      lower_rmse2 = alow(1); upper_rmse2 = aupp(1)
      
      do ir = 1, niter
          ahpd(ir) = seq_rmse3(ir)
      enddo
      call hpd(niter, 0.05d0, ahpd, alow, aupp)
      lower_rmse3 = alow(1); upper_rmse3 = aupp(1)
      
      do i = 1, ni
          write(21,1) i, tti(i), zzi(i), 
     +                   mean_tmean(i), std_tmean(i), 
     +                   lower_tmean(i), upper_tmean(i)
      enddo

      write(22,2) mean_rmse, std_rmse, 
     +            lower_rmse, upper_rmse
      write(22,*)

      write(22,2) mean_rmse1, std_rmse1, 
     +            lower_rmse1, upper_rmse1
      write(22,*)

      write(22,2) mean_rmse2, std_rmse2, 
     +            lower_rmse2, upper_rmse2
      write(22,*)

      write(22,2) mean_rmse3, std_rmse3, 
     +            lower_rmse3, upper_rmse3
      write(22,*)
      
      do i = 1, ni
          write(23,1) i, tti(i), zzi(i), 
     +                   (mean_pmean(i,l), std_pmean(i,l), 
     +                    lower_pmean(i,l), upper_pmean(i,l),
     +                    l = 1, nl)
      enddo
            
    1 format(i6, 100f20.10)   
    2 format(100f20.10)   

	end program		

	include 'PLaborFE_Gibbs.f'
	include 'optim1.f'
	include 'gilks2.f'
	include 'tnorm.f'
      include 'hpd.f'                        
