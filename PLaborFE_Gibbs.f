	subroutine gibbs(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 191, nj = 12, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc), sigmay

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 delta, ezeta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb), omega(nb,nb) 
		
      real*8 pp(ni,nl)
      
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
      
      call gibbs_bi12(iseed)
      
      call gibbs_bi3(iseed)
      
      call gibbs_bi4(iseed)  
      
	call gibbs_zi(iseed)  
      
	call gibbs_ti(iseed)      
                      
      end subroutine
     

	subroutine gibbs_ti(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 191, nj = 12, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 delta, ezeta, sigmas, v1
	real*8 zi(ni), bmu, bi(ni,nb)
      
      real*8 labor1, labor2
	real*8 si, xxt(nt), xt	
      real*8 zeta, epsilon
      real*8 shape, scale, rv
	      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
            
      common /vdelta/delta
      common /vezeta/ezeta
      
      common /vsigmas/sigmas
      common /vv1/v1

      common /vzi/zi
      
      common /vbmu/bmu
      common /vbi/bi

      external drnexp, drnstt
      
      do i = 1, ni
      
          if (zi(i) .eq. 1.d0) then

              labor1 = 1.d0; labor2 = 0.d0
          
          else if (zi(i) .eq. 2.d0) then

              labor1 = 0.d0; labor2 = 1.d0
          
          else 

              labor1 = 0.d0; labor2 = 0.d0

          endif
          
          labor1 = (labor1 - 0.1591928251d0)/0.3662666476d0
          labor2 = (labor2 - 0.1457399103d0)/0.353241509d0
            
          ui(i,3) = labor1
          ui(i,4) = labor2
            
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
      
          call rnset(iseed)
          call drnexp(1,zeta)
          call rnget(iseed)

          call rnset(iseed)
          call drnstt(1,v1,rv)
          call rnget(iseed)
                    
          epsilon = dsqrt(sigmas/v1)*rv

          si = xt + delta*(zeta - ezeta) + epsilon
      
          ti(i) = si + bmu + bi(i,nb)    
          
      enddo
     
	end subroutine
      
            
	subroutine gibbs_zi(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 191, nj = 12, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 
	real*8 bi(ni,nb)

      real*8 pp(ni,nl)
      
      real*8 xxp(np), xb, prob(nl)
      real*8 sump, cprob(nl), u
      		
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi

      common /vbi/bi

      common /vpp/pp
      
      external drnunf
      
      do i = 1, ni
          
          do l = 1, nl
              prob(l) = 1.d0
          enddo      
          sump = 0.d0
          do l = 1, nl-1
        
              do j = 1, nx
                  xxp(j) = xi(i,j)
              enddo
              do j = 1, nb-1
                  xxp(nx+j) = bi(i,j)
              enddo
                   
              xb = 0.d0
              do j = 1, np
                  xb = xb + xxp(j)*phi(j,l)
              enddo
            
              prob(l) = dexp(xb)
              sump = sump + prob(l)
          
          enddo
                          
          do l = 1, nl          
              prob(l) = prob(l)/(1.d0 + sump)
              pp(i,l) = prob(l)
          enddo
      
          cprob(1) = prob(1)
          do l = 2, nl          
              cprob(l) = cprob(l-1) + prob(l)
          enddo
                       
          call rnset(iseed)
          u = drnunf()
          call rnget(iseed)
                		                        
          if (u .le. cprob(1)) then
              
              zi(i) = 1.d0
              
          else if (u .le. cprob(2)) then
              
              zi(i) = 2.d0
              
          else
              
              zi(i) = 3.d0
              
          endif
          
      enddo
            
      end subroutine

      
	subroutine gibbs_bi4(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 191, nj = 12, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3
	integer, parameter :: npoint = 5, cpoint = 3, ndim = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc), sigmay

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      real*8 epsilon, astar, error
      real*8 yystar(npoint), xxstar(npoint,ndim)
      real*8 xx(ndim,ndim), xy(ndim)
      real*8 xxi(ndim,ndim), xhat(ndim)

	integer idum
	real*8 omegai(nb,nb)

      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      common /vsigmay/sigmay
            
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega
      
      common /vidum/idum
      common /vomegai/omegai

      external fbi4, dlinrg, drnnof, drnunf

      epsilon = 0.1d0

      call dlinrg(nb, omega, nb, omegai, nb)
      
      do i = 1, ni
      
          idum = i
        
          bold = bi(i,4)

          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(fbi4, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)

  11      do kk = 1, npoint
                  
              error = dfloat(kk - cpoint)
              astar = amean + error*epsilon
                  
              yystar(kk) = fbi4(astar)
              xxstar(kk,1) = astar**2
              xxstar(kk,2) = astar
              xxstar(kk,3) = 1.d0
                  
          enddo
              
          do kk = 1, npoint
              if (kk .ne. cpoint) then                  
                  if (yystar(kk) .le. ynewlo) then
                      epsilon = epsilon*1.2d0
                      goto 11
                  endif
              endif
          enddo
        
          do kk1 = 1, ndim
              do kk2 = 1, ndim
                  xx(kk1,kk2) = 0.d0
                  do kk = 1, npoint
                      xx(kk1,kk2) = xx(kk1,kk2) 
     +                            + xxstar(kk,kk1)
     +                              *xxstar(kk,kk2)
                  enddo
              enddo
              xy(kk1) = 0.d0
              do kk2 = 1, npoint
                  xy(kk1) = xy(kk1) + xxstar(kk2,kk1)
     +                                *yystar(kk2)
              enddo
          enddo

          call dlinrg(ndim, xx, ndim, xxi, ndim)	
                        
          do kk1 = 1, ndim
              xhat(kk1) = 0.d0
              do kk2 = 1, ndim
                  xhat(kk1) = xhat(kk1)
     +                      + xxi(kk1,kk2)*xy(kk2)
              enddo
          enddo
                                                       
          asigma = 1.0d0/(xhat(1)*2.d0)            
          if (asigma .lt. 0.0d0) asigma = -asigma

          bpdf = -fbi4(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)

          do ii = 1, 50
	  
              call rnset(iseed)
	        rv = drnnof()
	        call rnget(iseed)
		    anew = amean + rv*dsqrt(asigma)

		    apdf = -fbi4(anew) 
     +             + (anew - amean)**2/(2.d0*asigma)
		    ratio = apdf - bpdf 
     
              if (ratio .ge. 0.0d0) then
		        bpdf = apdf
		        bold = anew
	        else
		        call rnset(iseed)
		        u = drnunf()
		        call rnget(iseed)
		        if (dlog(u) .le. ratio) then
		            bpdf = apdf
		            bold = anew
		        endif
		    endif
		
          enddo     
       
          bi(i,4) = bold
                                   
      enddo
              
      end subroutine

    
	subroutine gibbs_bi4_mle(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 191, nj = 12, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc), sigmay

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	integer idum
	real*8 omegai(nb,nb)

      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      common /vsigmay/sigmay
            
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega
      
      common /vidum/idum
      common /vomegai/omegai

      external fbi4, dlinrg

      call dlinrg(nb, omega, nb, omegai, nb)
      
      do i = 1, ni
      
          idum = i
        
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bi(i,4)
          call nelmin(fbi4, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    bi(i,4) = xmin(1)
                                   
      enddo
              
      end subroutine
      
           
      real*8 function fbi4(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 191, nj = 12, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc), sigmay

	real*8 bmu, bi(ni,nb)

	integer idum
	real*8 omegai(nb,nb)

	real*8 star, sij, ystar
	real*8 bstar(nb), temp
      real*8 sum1, sum2, pdf
            
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      common /vsigmay/sigmay
            
      common /vbmu/bmu
      common /vbi/bi
      
      common /vidum/idum
      common /vomegai/omegai

      external dblinf, dlftrg, dlfdrg

      i = idum
      
      bi(i,4) = star
      
      sum1 = 0.d0             
      if (ij(i) .lt. 3) then
          iji = 2
      else
          iji = 3
      endif
      do j = 1, iji
        
          sij = tij(i,j) - bmu - bi(i,nb)
            
          if (sij .ge. 0.d0) then
                
              ystar = yij(i,j) 
     +              - (betaa(1) + bi(i,1))*sij   
     +              - (betaa(2) + bi(i,2))*sij**2   
            
          else

              ystar = yij(i,j) 
     +              - (betaa(3) + bi(i,3))*sij   
            
          endif
                      
          sum1 = sum1 - ystar**2/(2.d0*sigmay)
                                              
      enddo

      do j = 1, nb
          bstar(j) = bi(i,j)
      enddo
      temp = dblinf(nb,nb,omegai,nb,bstar,bstar)
      sum2 = -temp/2.d0
      
      pdf = sum1 + sum2
                        
      fbi4 = -pdf
     
      end function  
      
            
	subroutine gibbs_bi3(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 191, nj = 12, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc), sigmay

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 omegai(nb,nb), sij, ystar, summ
      real*8 hmean, hsigma, amean, asigma, rv

      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      common /vsigmay/sigmay
            
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega
      
      external dlinrg, drnnof

      call dlinrg(nb, omega, nb, omegai, nb)
   
      do i = 1, ni
      
          summ = 0.d0 
          do jj = 1, nb
              if (jj .ne. 3) then
                  summ = summ + omegai(3,jj)*bi(i,jj)
              endif
          enddo
      
          hmean = -summ
          hsigma = omegai(3,3)
          
          if (ij(i) .lt. 3) then
              iji = 2
          else
              iji = 3
          endif
          do j = 1, iji
        
              sij = tij(i,j) - bmu - bi(i,nb)
            
              if (sij .lt. 0.d0) then
                
                  ystar = yij(i,j) - betaa(3)*sij   
            
                  hmean = hmean + sij*ystar/sigmay
                  hsigma = hsigma + sij**2/sigmay
                                        
              endif
                    
          enddo
      
          amean = hmean/hsigma
          asigma = 1.d0/hsigma

          call rnset(iseed)
	    rv = drnnof()
	    call rnget(iseed)

          bi(i,3) = amean + rv*sqrt(asigma)
                                                
      enddo
      
      end subroutine
          
		      
	subroutine gibbs_bi12(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 191, nj = 12, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc), sigmay
	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 omegai(nb,nb)
	real*8 omegai1(2,2), omegai2(2,2)
	real*8 wij(2), bstar(2)
      real*8 sij, ystar

      real*8 hmean(2), hsigma(2,2)
      real*8 amean(2), asigma(2,2)
	real*8 tol, rsig(2,2), rv(2)
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      common /vsigmay/sigmay
            
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega
      
      external dlinrg, dmach, dchfac, drnmvn

      call dlinrg(nb, omega, nb, omegai, nb)
      
      do j1 = 1, 2
          do j2 = 1, 2
              omegai1(j1,j2) = omegai(j1,j2)
              omegai2(j1,j2) = omegai(j1,2+j2)
          enddo
      enddo
      
      do i = 1, ni
      
          bstar(1) = bi(i,3)
          bstar(2) = bi(i,4)
          do j1 = 1, 2
          
              hmean(j1) = 0.d0
              do j2 = 1, 2
              
                  hmean(j1) = hmean(j1) 
     +                      - omegai2(j1,j2)*bstar(j2)
              
                  hsigma(j1,j2) = omegai1(j1,j2)
              
              enddo
          
          enddo
      
          if (ij(i) .lt. 3) then
              iji = 2
          else
              iji = 3
          endif
          do j = 1, iji
        
              sij = tij(i,j) - bmu - bi(i,nb)
            
              if (sij .ge. 0.d0) then
                
                  ystar = yij(i,j) 
     +                  - betaa(1)*sij   
     +                  - betaa(2)*sij**2   
              
                  wij(1) = sij
                  wij(2) = sij**2                            
                                                
                  do j1 = 1, 2
        
                      hmean(j1) = hmean(j1) 
     +                          + wij(j1)*ystar/sigmay
     
                      do j2 = 1, 2
                          hsigma(j1,j2) = hsigma(j1,j2) 
     +                                  + wij(j1)*wij(j2)/sigmay
                      enddo
            
                  enddo
                    
              endif
          
          enddo

	    call dlinrg(2, hsigma, 2, asigma, 2)
	
	    do j1 = 1, 2
	        amean(j1) = 0.d0
		    do j2 = 1, 2
		        amean(j1) = amean(j1) + asigma(j1,j2)*hmean(j2)
		    enddo
	    enddo

	    tol = 100.d0*dmach(4)
	    call dchfac(2, asigma, 2, tol, irank, rsig, 2)
	    call rnset(iseed)
	    call drnmvn(1, 2, rsig, 2, rv, 1)
	    call rnget(iseed)
        
          bi(i,1) = amean(1) + rv(1)
          bi(i,2) = amean(2) + rv(2)
          
      enddo
                                                        
      end subroutine