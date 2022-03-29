	subroutine gibbs(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc), sigmay

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 lambda(ni), zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb), omega(nb,nb) 
	
	real*8 sigmab, sigmat, sigmad, sigmap, sigmam
	real*8 a0, b0, a1, b1, a2, b2, d0, s0
	
      common /vij/ij
	
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      common /vsigmay/sigmay
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      common /vlambda/lambda
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1

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
      common /va2/a2
      common /vb2/b2
      common /vd0/d0
      common /vs0/s0
            
      call gibbs_zeta(iseed)

      call gibbs_bi1(iseed) 
      
      call gibbs_bi2(iseed)  

      call gibbs_bi3(iseed) 
      
      call gibbs_bi4(iseed)      
      
	call gibbs_bmu(iseed)


	call gibbs_beta12(iseed)

	call gibbs_beta3(iseed)

	call gibbs_sigmay(iseed)

      
      call gibbs_theta(iseed)
      
      call gibbs_delta(iseed)
            
	call gibbs_sigmas(iseed)
      
	call gibbs_v1(iseed)

	call gibbs_lambda(iseed)
      
                    
	call gibbs_phi(iseed)      

      
	call gibbs_omega(iseed)
            
	end subroutine


	subroutine gibbs_bi3(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3
	integer, parameter :: npoint = 5, cpoint = 3, ndim = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 a0, b0

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
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
 
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external fbi3, dlinrg, drnnof, drnunf

      epsilon = 0.1d0

      call dlinrg(nb, omega, nb, omegai, nb)

      do i = 1, ni
      
          idum = i
        
          bold = bi(i,3)

          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(fbi3, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)

  11      do kk = 1, npoint
                  
              error = dfloat(kk - cpoint)
              astar = amean + error*epsilon
                  
              yystar(kk) = fbi3(astar)
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

          bpdf = -fbi3(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)

          do ii = 1, 20
	  
              call rnset(iseed)
	        rv = drnnof()
	        call rnget(iseed)
		    anew = amean + rv*dsqrt(asigma)

		    apdf = -fbi3(anew) 
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
       
          bi(i,3) = bold
                              
      enddo
                             
      end subroutine


	subroutine gibbs_bi3_MLE(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 a0, b0

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	integer idum
	real*8 omegai(nb,nb)

      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
 
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external fbi3, dlinrg, drnnof, drnunf

      call dlinrg(nb, omega, nb, omegai, nb)

      do i = 1, ni
      
          idum = i
        
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bi(i,3)
          call nelmin(fbi3, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    bi(i,3) = xmin(1)
                 
      enddo
                             
      end subroutine

      
      real*8 function fbi3(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb)

	real*8 a0, b0

	integer idum
	real*8 omegai(nb,nb)

	real*8 star, sij, ystar
	real*8 bstar(nb), temp
	real*8 si, xxt(nt), xt
      real*8 xxp(np), xp
      real*8 suma, sumb, sumc
      real*8 sum1, sum2      
      real*8 pdf1, pdf2, pdf
            
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
 
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external dlinrg, dblinf, dlftrg, dlfdrg

      i = idum
      
      bi(i,3) = star
      
      suma = a0; sumb = b0
      do ii = 1, ni
          do j = 1, ij(ii)
        
              sij = tij(ii,j) - bmu - bi(ii,nb)
            
              if (sij .ge. 0.d0) then
                
                  ystar = yij(ii,j) 
     +                  - (betaa(1) + bi(ii,1))*sij   
     +                  - (betaa(2) + bi(ii,2))*sij**2   
            
              else

                  ystar = yij(ii,j) 
     +                  - (betaa(3) + bi(ii,3))*sij   
            
              endif
            
              suma = suma + 1.d0/2.d0
              sumb = sumb + ystar**2/2.d0 
            
          enddo
      enddo

      do j = 1, nb
          bstar(j) = bi(i,j)
      enddo
      temp = dblinf(nb,nb,omegai,nb,bstar,bstar)
      sumc = -temp/2.d0
      
      do j = 1, nx
          xxp(j) = xi(i,j)
      enddo
      do j = 1, nb-1
          xxp(nx+j) = bi(i,j)
      enddo

      sum1 = 0.d0
      do l = 1, nl-1              
        
          xp = 0.d0
          do j = 1, np
              xp = xp + xxp(j)*phi(j,l)
          enddo
           
          sum1 = sum1 + dexp(xp)

          sum2 = 0.d0
          if (zi(i) .eq. dfloat(l)) then
        
              sum2 = 0.d0
              do j = 1, np
                  sum2 = sum2 + xxp(j)*phi(j,l)
              enddo
                        
          endif            
            
      enddo

      pdf1 = -suma*dlog(sumb) + sumc 
     +     + sum2 - dlog(1.d0 + sum1)
                        

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

      ystar = si - xt - delta*(zeta(i) - ezeta)
          
      pdf2 = -(v1 + 1.d0)/2.d0*dlog(1.d0 + ystar**2/sigmas)

      
      pdf = pdf1 + pdf2
                             
      fbi3 = -pdf
     
      end function
      
      

	subroutine gibbs_bi2(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3
	integer, parameter :: npoint = 5, cpoint = 3, ndim = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 a0, b0

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
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
 
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external fbi2, dlinrg, drnnof, drnunf

      epsilon = 0.1d0

      call dlinrg(nb, omega, nb, omegai, nb)

      do i = 1, ni
      
          idum = i
        
          bold = bi(i,2)

          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(fbi2, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)

  11      do kk = 1, npoint
                  
              error = dfloat(kk - cpoint)
              astar = amean + error*epsilon
                  
              yystar(kk) = fbi2(astar)
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

          bpdf = -fbi2(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)

          do ii = 1, 20
	  
              call rnset(iseed)
	        rv = drnnof()
	        call rnget(iseed)
		    anew = amean + rv*dsqrt(asigma)

		    apdf = -fbi2(anew) 
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
       
          bi(i,2) = bold
                                   
      enddo
                       
      end subroutine
       
   
	subroutine gibbs_bi2_MLE(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 a0, b0

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	integer idum
	real*8 omegai(nb,nb)

      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
 
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external fbi2, dlinrg, drnnof, drnunf

      call dlinrg(nb, omega, nb, omegai, nb)

      do i = 1, ni
      
          idum = i
        
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bi(i,2)
          call nelmin(fbi2, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    bi(i,2) = xmin(1)
                          
      enddo
                       
      end subroutine
       

      real*8 function fbi2(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb)

	real*8 a0, b0

	integer idum
	real*8 omegai(nb,nb)

	real*8 star, sij, ystar
	real*8 bstar(nb), temp
	real*8 si, xxt(nt), xt
      real*8 xxp(np), xp
      real*8 suma, sumb, sumc
      real*8 sum1, sum2      
      real*8 pdf1, pdf2, pdf
            
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
 
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external dlinrg, dblinf, dlftrg, dlfdrg

      i = idum
      
      bi(i,2) = star
      
      suma = a0; sumb = b0
      do ii = 1, ni
      
          do j = 1, ij(ii)
        
              sij = tij(ii,j) - bmu - bi(ii,nb)
            
              if (sij .ge. 0.d0) then
                
                  ystar = yij(ii,j) 
     +                  - (betaa(1) + bi(ii,1))*sij   
     +                  - (betaa(2) + bi(ii,2))*sij**2   
            
              else

                  ystar = yij(ii,j) 
     +                  - (betaa(3) + bi(ii,3))*sij   
            
              endif
            
              suma = suma + 1.d0/2.d0
              sumb = sumb + ystar**2/2.d0 
            
          enddo
        
      enddo
      
      do j = 1, nb
          bstar(j) = bi(i,j)
      enddo
      temp = dblinf(nb,nb,omegai,nb,bstar,bstar)
      sumc = -temp/2.d0
      

      do j = 1, nx
          xxp(j) = xi(i,j)
      enddo
      do j = 1, nb-1
          xxp(nx+j) = bi(i,j)
      enddo

      sum1 = 0.d0
      do l = 1, nl-1              
        
          xp = 0.d0
          do j = 1, np
              xp = xp + xxp(j)*phi(j,l)
          enddo
           
          sum1 = sum1 + dexp(xp)

          sum2 = 0.d0
          if (zi(i) .eq. dfloat(l)) then
        
              sum2 = 0.d0
              do j = 1, np
                  sum2 = sum2 + xxp(j)*phi(j,l)
              enddo
                        
          endif            
            
      enddo

      pdf1 = -suma*dlog(sumb) + sumc 
     +     + sum2 - dlog(1.d0 + sum1)
                        

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

      ystar = si - xt - delta*(zeta(i) - ezeta)
          
      pdf2 = -(v1 + 1.d0)/2.d0*dlog(1.d0 + ystar**2/sigmas)

      
      pdf = pdf1 + pdf2
                             
      fbi2 = -pdf
     
      end function
      
      

	subroutine gibbs_bi1(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3
	integer, parameter :: npoint = 5, cpoint = 3, ndim = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 a0, b0

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
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
 
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external fbi1, dlinrg, drnnof, drnunf

      epsilon = 0.1d0

      call dlinrg(nb, omega, nb, omegai, nb)

      do i = 1, ni
      
          idum = i
        
          bold = bi(i,1)

          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(fbi1, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)

  11      do kk = 1, npoint
                  
              error = dfloat(kk - cpoint)
              astar = amean + error*epsilon
                  
              yystar(kk) = fbi1(astar)
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

          bpdf = -fbi1(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)

          do ii = 1, 20
	  
              call rnset(iseed)
	        rv = drnnof()
	        call rnget(iseed)
		    anew = amean + rv*dsqrt(asigma)

		    apdf = -fbi1(anew) 
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
       
          bi(i,1) = bold
                                            
      enddo
                             
      end subroutine
      

	subroutine gibbs_bi1_MLE(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 a0, b0

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	integer idum
	real*8 omegai(nb,nb)

      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
 
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external fbi1, dlinrg, drnnof, drnunf

      call dlinrg(nb, omega, nb, omegai, nb)

      do i = 1, ni
      
          idum = i
        
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bi(i,1)
          call nelmin(fbi1, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    bi(i,1) = xmin(1)
                           
      enddo
                             
      end subroutine
  
      
      real*8 function fbi1(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 

	real*8 bmu, bi(ni,nb)

	real*8 a0, b0

	integer idum
	real*8 omegai(nb,nb)

	real*8 star, sij, ystar
	real*8 bstar(nb), temp
	real*8 si, xxt(nt), xt
      real*8 xxp(np), xp
      real*8 suma, sumb, sumc
      real*8 sum1, sum2      
      real*8 pdf1, pdf2, pdf
            
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
 
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi
      
      common /vbmu/bmu
      common /vbi/bi

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external dlinrg, dblinf, dlftrg, dlfdrg

      i = idum
      
      bi(i,1) = star
      
      suma = a0; sumb = b0
      do ii = 1, ni
      
          do j = 1, ij(ii)
        
              sij = tij(ii,j) - bmu - bi(ii,nb)
            
              if (sij .ge. 0.d0) then
                
                  ystar = yij(ii,j) 
     +                  - (betaa(1) + bi(ii,1))*sij   
     +                  - (betaa(2) + bi(ii,2))*sij**2   
            
              else

                  ystar = yij(ii,j) 
     +                  - (betaa(3) + bi(ii,3))*sij   
            
              endif
            
              suma = suma + 1.d0/2.d0
              sumb = sumb + ystar**2/2.d0 
            
          enddo
        
      enddo

      do j = 1, nb
          bstar(j) = bi(i,j)
      enddo
      temp = dblinf(nb,nb,omegai,nb,bstar,bstar)
      sumc = -temp/2.d0
      
      
      do j = 1, nx
          xxp(j) = xi(i,j)
      enddo
      do j = 1, nb-1
          xxp(nx+j) = bi(i,j)
      enddo

      sum1 = 0.d0
      do l = 1, nl-1              
        
          xp = 0.d0
          do j = 1, np
              xp = xp + xxp(j)*phi(j,l)
          enddo
           
          sum1 = sum1 + dexp(xp)

          sum2 = 0.d0
          if (zi(i) .eq. dfloat(l)) then
        
              sum2 = 0.d0
              do j = 1, np
                  sum2 = sum2 + xxp(j)*phi(j,l)
              enddo
                        
          endif            
            
      enddo

      pdf1 = -suma*dlog(sumb) + sumc 
     +     + sum2 - dlog(1.d0 + sum1)
                        
              
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

      ystar = si - xt - delta*(zeta(i) - ezeta)
          
      pdf2 = -(v1 + 1.d0)/2.d0*dlog(1.d0 + ystar**2/sigmas)
                                        
      
      pdf = pdf1 + pdf2
                             
      fbi1 = -pdf
     
      end function
      
      

	subroutine gibbs_bi4(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3
	integer, parameter :: npoint = 5, cpoint = 3, ndim = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 a0, b0

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
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
       
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega

      common /va0/a0
      common /vb0/b0
      
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

          do ii = 1, 20
	  
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

  
	subroutine gibbs_bi4_MLE(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 bmu, bi(ni,nb), omega(nb,nb) 

	real*8 a0, b0

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	integer idum
	real*8 omegai(nb,nb)

      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
       
      common /vbmu/bmu
      common /vbi/bi
      common /vomega/omega

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external fbi4, dlinrg, drnnof, drnunf

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
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1

	real*8 bmu, bi(ni,nb)

	real*8 a0, b0

	integer idum
	real*8 omegai(nb,nb)

	real*8 star, sij, ystar
	real*8 bstar(nb), temp
	real*8 si, xxt(nt), xt
      real*8 suma, sumb, sumc
      real*8 pdf1, pdf2, pdf
            
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1
       
      common /vbmu/bmu
      common /vbi/bi

      common /va0/a0
      common /vb0/b0
      
      common /vidum/idum
      common /vomegai/omegai

      external dlinrg, dblinf, dlftrg, dlfdrg

      i = idum
      
      bi(i,4) = star
      
      suma = a0; sumb = b0
      do ii = 1, ni
      
          do j = 1, ij(ii)
        
              sij = tij(ii,j) - bmu - bi(ii,nb)
            
              if (sij .ge. 0.d0) then
                
                  ystar = yij(ii,j) 
     +                  - (betaa(1) + bi(ii,1))*sij   
     +                  - (betaa(2) + bi(ii,2))*sij**2   
            
              else

                  ystar = yij(ii,j) 
     +                  - (betaa(3) + bi(ii,3))*sij   
            
              endif
            
              suma = suma + 1.d0/2.d0
              sumb = sumb + ystar**2/2.d0 
            
          enddo
        
      enddo
      
      do j = 1, nb
          bstar(j) = bi(i,j)
      enddo
      temp = dblinf(nb,nb,omegai,nb,bstar,bstar)
      sumc = -temp/2.d0
      
      pdf1 = -suma*dlog(sumb) + sumc 
                        
      
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

      ystar = si - xt - delta*(zeta(i) - ezeta)
          
      pdf2 = -(v1 + 1.d0)/2.d0*dlog(1.d0 + ystar**2/sigmas)

      
      pdf = pdf1 + pdf2
                             
      fbi4 = -pdf
     
      end function
      
      
	subroutine gibbs_bmu(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3
	integer, parameter :: npoint = 5, cpoint = 3, ndim = 3
      
      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1
	
	real*8 bmu, bi(ni,nb)

	real*8 sigmam, a0, b0
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      real*8 epsilon, astar, error
      real*8 yystar(npoint), xxstar(npoint,ndim)
      real*8 xx(ndim,ndim), xy(ndim)
      real*8 xxi(ndim,ndim), xhat(ndim)
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1

      common /vbmu/bmu
      common /vbi/bi

      common /vsigmam/sigmam

      common /va0/a0
      common /vb0/b0

      external fbmu, dlinrg, drnnof, drnunf

      epsilon = 0.1d0
      
      bold = bmu

      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fbmu, nopt, start, xmin, ynewlo,  
     +            reqmin, step, konvge, kcount, 
     +            icount, numres, ifault)
	amean = xmin(1)
      
  11  do kk = 1, npoint
                  
          error = dfloat(kk - cpoint)
          astar = amean + error*epsilon
                  
          yystar(kk) = fbmu(astar)
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
     +                        + xxstar(kk,kk1)
     +                          *xxstar(kk,kk2)
              enddo
          enddo
          xy(kk1) = 0.d0
          do kk2 = 1, npoint
              xy(kk1) = xy(kk1) + xxstar(kk2,kk1)
     +                            *yystar(kk2)
          enddo
      enddo

      call dlinrg(ndim, xx, ndim, xxi, ndim)	
                        
      do kk1 = 1, ndim
          xhat(kk1) = 0.d0
          do kk2 = 1, ndim
              xhat(kk1) = xhat(kk1)
     +                  + xxi(kk1,kk2)*xy(kk2)
          enddo
      enddo
                                                       
      asigma = 1.0d0/(xhat(1)*2.d0)            
      if (asigma .lt. 0.0d0) asigma = -asigma

      bpdf = -fbmu(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)

      do ii = 1, 20
	  
          call rnset(iseed)
	    rv = drnnof()
	    call rnget(iseed)
	    anew = amean + rv*dsqrt(asigma)

	    apdf = -fbmu(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
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
        
      bmu = bold
            
	end subroutine


      real*8 function fbmu(bmu)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1
	
	real*8 bmu, bi(ni,nb)

	real*8 sigmam, a0, b0

	real*8 sij, ystar
	real*8 si, xxt(nt), xt
      real*8 suma, sumb
      real*8 pdf1, pdf2, pdf
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1

      common /vbi/bi

      common /vsigmam/sigmam

      common /va0/a0
      common /vb0/b0

      external dlinrg, dblinf, dlftrg, dlfdrg
      
      suma = a0; sumb = b0
      do i = 1, ni
      
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
            
              suma = suma + 1.d0/2.d0
              sumb = sumb + ystar**2/2.d0 
            
          enddo
        
      enddo      
      pdf1 = -suma*dlog(sumb) 
       
      pdf2 = 0.d0                
      do i = 1, ni
              
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

          ystar = si - xt - delta*(zeta(i) - ezeta)
          
          pdf2 = pdf2 - (v1 + 1.d0)/2.d0
     +                  *dlog(1.d0 + ystar**2/sigmas)
          
      enddo
                  
      pdf = pdf1 + pdf2 - bmu**2/(2.d0*sigmam)
                                                    
      fbmu = -pdf
    
      end function

     
	subroutine gibbs_beta3(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7 
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)
	real*8 bmu, bi(ni,nb)
	real*8 sigmab, a0, b0

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt
      
	real*8 bold, anew, amean, asigma
      real*8 rv, u, bpdf, apdf, ratio

	real*8 sij, ystar, suma, sumb, sumc, sumd
      
      common /vij/ij
	
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vbmu/bmu
      common /vbi/bi

      common /vsigmab/sigmab
      common /va0/a0
      common /vb0/b0

      external fbeta3, dlinrg, dblinf, drnnof, drnunf

      bold = betaa(3)

      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 2000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fbeta3, nopt, start, xmin, ynewlo,  
     +            reqmin, step, konvge, kcount, 
     +            icount, numres, ifault)
	amean = xmin(1)

      betaa(3) = amean

      suma = a0; sumb = b0
      sumc = 0.d0; sumd = 0.d0      
      do i = 1, ni
      
          do j = 1, ij(i)
        
              sij = tij(i,j) - bmu - bi(i,nb)
            
              if (sij .ge. 0.d0) then

                  ystar = yij(i,j) 
     +                  - (betaa(1) + bi(i,1))*sij   
     +                  - (betaa(2) + bi(i,2))*sij**2   
                                
              else
            
                  ystar = yij(i,j) 
     +                  - (betaa(3) + bi(i,3))*sij   

                  sumc = sumc + sij**2
                  
                  sumd = sumd + ystar*sij
            
              endif
            
              suma = suma + 1.d0/2.d0
              sumb = sumb + ystar**2/2.d0 
            
          enddo
                                  
      enddo

      asigma = -suma*(sumc*sumb - sumd**2)/sumb**2
     +       - 1.d0/sigmab

      asigma = -1.d0/asigma
      if (asigma .lt. 0.0d0) asigma = -asigma      
            
      bpdf = -fbeta3(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)
                
      do ii = 1, 20
	  
          call rnset(iseed)
	    rv = drnnof()
	    call rnget(iseed)
		anew = amean + rv*dsqrt(asigma)

		apdf = -fbeta3(anew) 
     +       + (anew - amean)**2/(2.d0*asigma)
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

      betaa(3) = bold
      
	end subroutine


      real*8 function fbeta3(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7 
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)
	real*8 bmu, bi(ni,nb)
	real*8 sigmab, a0, b0

	real*8 star, sij, ystar		
      real*8 suma, sumb, sumc, pdf
      
      common /vij/ij
	
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vbmu/bmu
      common /vbi/bi

      common /vsigmab/sigmab
      common /va0/a0
      common /vb0/b0

      betaa(3) = star
      
      suma = a0; sumb = b0
      do i = 1, ni
      
        do j = 1, ij(i)
        
            sij = tij(i,j) - bmu - bi(i,nb)
            
            if (sij .ge. 0.d0) then
                
                ystar = yij(i,j) 
     +                - (betaa(1) + bi(i,1))*sij   
     +                - (betaa(2) + bi(i,2))*sij**2   
            
            else
            
                ystar = yij(i,j) 
     +                - (betaa(3) + bi(i,3))*sij   
            
            endif
            
            suma = suma + 1.d0/2.d0
            sumb = sumb + ystar**2/2.d0 
            
        enddo
        
      enddo
      
      sumc = -betaa(3)**2/(2.d0*sigmab)

      pdf = -suma*dlog(sumb) + sumc
                       
      fbeta3 = -pdf
    
	end function


	subroutine gibbs_beta12(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7 
	integer, parameter :: nx = 2, np = 5, nl = 3
	integer, parameter :: nc1 = nc-1

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)
	real*8 bmu, bi(ni,nb)
	real*8 sigmab, a0, b0

	real*8 start(nc1), xmin(nc1), ynewlo, reqmin, step(nc1)
	integer konvge, kcount, icount, numres, ifault, nopt

      real*8 bold(nc1), anew(nc1)
      real*8 amean(nc1), asigma(nc1,nc1), der2(nc1,nc1)
	real*8 tol, rsig(nc1,nc1), hmean(nc1)
	real*8 hpdf, bpdf, apdf, ratio, u

	real*8 sij, ystar, wws(nc1)      
      real*8 suma, sumb, sumc(nc1,nc1), sumd(nc1)
      
      common /vij/ij
	
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vbmu/bmu
      common /vbi/bi

      common /vsigmab/sigmab
      common /va0/a0
      common /vb0/b0

      external fbeta1, dlinrg, dchfac, dblinf, drnmvn, drnunf

      bold(1) = betaa(1)
      bold(2) = betaa(2)

	nopt = nc1
	reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	do j = 1, nc1
	  step(j) = 0.2d0 ; start(j) = bold(j)
	enddo	
	call nelmin(fbeta1, nopt, start, xmin, 
     +            ynewlo, reqmin, step, konvge, 
     +            kcount, icount, numres, ifault)       
      do j = 1, nc1
        amean(j) = xmin(j)
      enddo

      betaa(1) = amean(1)
      betaa(2) = amean(2)

      do j1 = 1, nc1
        do j2 = 1, nc1
            der2(j1,j2) = 0.d0
        enddo
      enddo
      do j1 = 1, nc1
        do j2 = j1, nc1

            suma = a0; sumb = b0
            do jj1 = 1, nc1
                do jj2 = 1, nc1
                    sumc(jj1,jj2) = 0.d0
                enddo
                sumd(jj1) = 0.d0
            enddo
            do i = 1, ni
      
                do j = 1, ij(i)
        
                    sij = tij(i,j) - bmu - bi(i,nb)
            
                    if (sij .ge. 0.d0) then

                        wws(1) = sij 
                        wws(2) = sij**2

                        ystar = yij(i,j) 
     +                        - (betaa(1) + bi(i,1))*sij   
     +                        - (betaa(2) + bi(i,2))*sij**2   

                        do jj1 = 1, nc1
                            do jj2 = 1, nc1
                                sumc(jj1,jj2) = sumc(jj1,jj2) 
     +                                        + wws(jj1)*wws(jj2)
                            enddo
                            sumd(jj1) = sumd(jj1) + ystar*wws(jj1)
                        enddo
                                
                    else
            
                        ystar = yij(i,j) 
     +                        - (betaa(3) + bi(i,3))*sij   
            
                    endif
            
                    suma = suma + 1.d0/2.d0
                    sumb = sumb + ystar**2/2.d0 
            
                enddo
                                  
            enddo

            der2(j1,j2) = -suma*(sumc(j1,j2)*sumb 
     +                           - sumd(j1)*sumd(j2))
     +                     /sumb**2

        enddo
      enddo

      der2(1,1) = der2(1,1) - 1.d0/sigmab    
      der2(2,2) = der2(2,2) - 1.d0/sigmab

      do j1 = 1, nc1
        do j2 = j1, nc1
            der2(j2,j1) = der2(j1,j2)
        enddo
      enddo                

      do j1 = 1, nc1
        do j2 = 1, nc1
            der2(j1,j2) = -der2(j1,j2)
        enddo
      enddo                

	call dlinrg(nc1, der2, nc1, asigma, nc1)	
	  	  
	tol = 100.0d0*dmach(4)
	call dchfac(nc1, asigma, nc1, tol, irank, rsig, nc1)

	do j = 1, nc1
	  hmean(j) = bold(j) - amean(j)
	enddo
	hpdf = dblinf(nc1, nc1, der2, nc1, hmean, hmean)
	bpdf = -fbeta1(bold) + hpdf/2.d0

	do ii = 1, 20
	
	    call rnset(iseed)
		call drnmvn(1, nc1, rsig, nc1, anew, 1)
		call rnget(iseed)		
		do j = 1, nc1
	      anew(j) = amean(j) + anew(j)
		enddo
		
		do j = 1, nc1
		    hmean(j) = anew(j) - amean(j)
		enddo
		hpdf = dblinf(nc1, nc1, der2, nc1, hmean, hmean)
		apdf = -fbeta1(anew) + hpdf/2.d0

		ratio = apdf - bpdf 

		if (ratio .ge. 0.0d0) then
		    bpdf = apdf
			do j = 1, nc1
				bold(j) = anew(j)
			enddo
		else
		    call rnset(iseed)
			u = drnunf()
			call rnget(iseed)
			if (dlog(u) .le. ratio) then
			    bpdf = apdf
				do j = 1, nc1
					bold(j) = anew(j)
				enddo
			endif
		endif
		
	enddo      

      betaa(1) = bold(1)
      betaa(2) = bold(2)

	end subroutine


      real*8 function fbeta1(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7 
	integer, parameter :: nx = 2, np = 5, nl = 3
	integer, parameter :: nc1 = nc-1

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc)
	real*8 bmu, bi(ni,nb)
	real*8 sigmab, a0, b0

	real*8 star(nc1), sij, ystar		
      real*8 suma, sumb, sumc, pdf
      
      common /vij/ij
	
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      
      common /vbmu/bmu
      common /vbi/bi

      common /vsigmab/sigmab
      common /va0/a0
      common /vb0/b0

      betaa(1) = star(1)
      betaa(2) = star(2)
      
      suma = a0; sumb = b0
      do i = 1, ni
      
        do j = 1, ij(i)
        
            sij = tij(i,j) - bmu - bi(i,nb)
            
            if (sij .ge. 0.d0) then
                
                ystar = yij(i,j) 
     +                - (betaa(1) + bi(i,1))*sij   
     +                - (betaa(2) + bi(i,2))*sij**2   
            
            else
            
                ystar = yij(i,j) 
     +                - (betaa(3) + bi(i,3))*sij   
            
            endif
            
            suma = suma + 1.d0/2.d0
            sumb = sumb + ystar**2/2.d0 
            
        enddo
        
      enddo
      
      sumc = -betaa(1)**2/(2.d0*sigmab)
     +     - betaa(2)**2/(2.d0*sigmab)

      pdf = -suma*dlog(sumb) + sumc
                       
      fbeta1 = -pdf
    
	end function


	subroutine gibbs_sigmay(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
	real*8 betaa(nc), sigmay
	real*8 bmu, bi(ni,nb)
	real*8 a0, b0

	real*8 sij, ystar		
      real*8 shape, scale, rv
      
      common /vij/ij
	
      common /vyij/yij
      common /vtij/tij
      common /vbetaa/betaa
      common /vsigmay/sigmay
      
      common /vbmu/bmu
      common /vbi/bi

      common /va0/a0
      common /vb0/b0

      external drngam

      shape = a0; scale = b0
      do i = 1, ni
      
        do j = 1, ij(i)
        
            sij = tij(i,j) - bmu - bi(i,nb)
            
            if (sij .ge. 0.d0) then
                
                ystar = yij(i,j) 
     +                - (betaa(1) + bi(i,1))*sij   
     +                - (betaa(2) + bi(i,2))*sij**2   
            
            else
            
                ystar = yij(i,j) 
     +                - (betaa(3) + bi(i,3))*sij   
            
            endif
            
            shape = shape + 1.d0/2.d0
            scale = scale + ystar**2/2.d0 
            
        enddo
        
      enddo

      call rnset(iseed)
	call drngam(1, shape, rv)
	call rnget(iseed)
	sigmay = scale/rv
     
      end subroutine


	subroutine gibbs_zeta(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 lambda(ni), zeta(ni), ezeta
	real*8 delta, sigmas
	real*8 bmu, bi(ni,nb)

	real*8 si, xxt(nt), xt, ystar		

      real*8 hmean, hsigma, amean, asigma, trim, rv

      logical la, lb
	
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vlambda/lambda
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      
      common /vsigmas/sigmas

      common /vbmu/bmu
      common /vbi/bi

      external ytuvn

      do i = 1, ni
              
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

          ystar = si - xt + delta*ezeta

          hmean = lambda(i)*ystar*delta/sigmas - 1.d0
          hsigma = lambda(i)*delta**2/sigmas

          amean = hmean/hsigma
          asigma = 1.d0/hsigma
                            
          trim = -amean/dsqrt(asigma)
            
          la = .false. ; lb = .true.

		call rnset(iseed)
		rv = ytuvn(trim, trim, la, lb, iseed)
		call rnget(iseed)
		    
		zeta(i) = amean + rv*dsqrt(asigma)         
                    
      enddo
    
      end subroutine
        
      
	subroutine gibbs_theta(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1
	
	real*8 bmu, bi(ni,nb)
      
	real*8 sigmat

	real*8 start(nt), xmin(nt), ynewlo, reqmin, step(nt)
	integer konvge, kcount, icount, numres, ifault, nopt

      real*8 bold(nt), anew(nt), amean(nt), asigma(nt,nt)
	real*8 der2(nt,nt), tol, rsig(nt,nt), hmean(nt)
	real*8 hpdf, bpdf, apdf, ratio, u
      
	real*8 si, xxt(nt), xt, ystar, zstar     
      real*8 temp1, temp2
      	      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1

      common /vbmu/bmu
      common /vbi/bi

      common /vsigmat/sigmat

      external ftheta, dlinrg, dchfac, dblinf, drnmvn, drnunf
      
      do j = 1, nt
          bold(j) = theta(j) 
      enddo
      
      nopt = nt
      reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
      do j = 1, nt
          step(j) = 0.2d0 ; start(j) = bold(j)
      enddo	
      call nelmin(ftheta, nopt, start, xmin, 
     +            ynewlo, reqmin, step, konvge, 
     +            kcount, icount, numres, ifault)       
      do j = 1, nt
          amean(j) = xmin(j)
          theta(j) = amean(j)
      enddo

      do j1 = 1, nt
          do j2 = 1, nt
              der2(j1,j2) = 0.d0
          enddo
      enddo
      do j1 = 1, nt
          do j2 = j1, nt
              
              do i = 1, ni
              
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
          
                  zstar = zeta(i) - ezeta  

                  ystar = si - xt - delta*zstar
          
                  temp1 = 1.d0 + ystar**2/sigmas

                  temp2 = 1.d0 - ystar**2/sigmas
                  
                  der2(j1,j2) = der2(j1,j2) 
     +                        - (v1 + 1.d0)*(1.d0/sigmas)
     +                          *temp2/temp1**2
     +                          *xxt(j1)*xxt(j2)
                    
              enddo      
                  
          enddo
      enddo

      do j = 1, nt
          der2(j,j) = der2(j,j) - 1.d0/sigmat     
      enddo
      
      do j1 = 1, nt
          do j2 = j1, nt
              der2(j2,j1) = der2(j1,j2)
          enddo
      enddo                
      do j1 = 1, nt
          do j2 = 1, nt
              der2(j1,j2) = -der2(j1,j2)
          enddo
      enddo                

      call dlinrg(nt, der2, nt, asigma, nt)	
	  	  
      tol = 100.0d0*dmach(4)
      call dchfac(nt, asigma, nt, tol, irank, rsig, nt)

      do j = 1, nt
          hmean(j) = bold(j) - amean(j)
      enddo
      hpdf = dblinf(nt, nt, der2, nt, hmean, hmean)
      bpdf = -ftheta(bold) + hpdf/2.d0

      do ii = 1, 20
	
          call rnset(iseed)
          call drnmvn(1, nt, rsig, nt, anew, 1)
          call rnget(iseed)		
          do j = 1, nt
              anew(j) = amean(j) + anew(j)
          enddo
		
          do j = 1, nt
              hmean(j) = anew(j) - amean(j)
          enddo
          hpdf = dblinf(nt, nt, der2, nt, hmean, hmean)
          apdf = -ftheta(anew) + hpdf/2.d0

          ratio = apdf - bpdf 

          if (ratio .ge. 0.0d0) then
              bpdf = apdf
              do j = 1, nt
                  bold(j) = anew(j)
              enddo
          else
              call rnset(iseed)
              u = drnunf()
              call rnget(iseed)
              if (dlog(u) .le. ratio) then
                  bpdf = apdf
                  do j = 1, nt
                      bold(j) = anew(j)
                  enddo
              endif
          endif
		
      enddo      
	
	do j = 1, nt
	  theta(j) = bold(j)
	enddo            
      
      end subroutine
      
     
      real*8 function ftheta(theta)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1
	
	real*8 bmu, bi(ni,nb)
      
	real*8 sigmat
      
	real*8 si, xxt(nt), xt, ystar, zstar     
      real*8 sum1, sum2, pdf
      	      
      common /vti/ti
      common /vui/ui
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1

      common /vbmu/bmu
      common /vbi/bi

      common /vsigmat/sigmat

      sum1 = 0.d0
      do i = 1, ni
              
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
          
          zstar = zeta(i) - ezeta  

          ystar = si - xt - delta*zstar
          
          sum1 = sum1 - (v1 + 1.d0)/2.d0
     +                  *dlog(1.d0 + ystar**2/sigmas)

      enddo

      sum2 = 0.d0
      do j = 1, nt
          sum2 = sum2 - theta(j)**2/(2.d0*sigmat)
      enddo
      
      pdf = sum1 + sum2
      
      ftheta = -pdf
     
	end function
       
      
	subroutine gibbs_delta(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1
	
	real*8 bmu, bi(ni,nb)
	real*8 sigmad

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
	real*8 si, xxt(nt), xt, ystar, zstar     
      real*8 temp1, temp2
	      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1

      common /vbmu/bmu
      common /vbi/bi

      common /vsigmad/sigmad
      
      external fdelta, dlinrg, drnnof, drnunf

      bold = delta

      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fdelta, nopt, start, xmin, ynewlo,  
     +            reqmin, step, konvge, kcount, 
     +            icount, numres, ifault)
	amean = xmin(1)
	delta = amean
      
      asigma = 0.d0
      do i = 1, ni
              
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
          
          zstar = zeta(i) - ezeta  

          ystar = si - xt - delta*zstar
          
          temp1 = 1.d0 + ystar**2/sigmas

          temp2 = 1.d0 - ystar**2/sigmas

          asigma = asigma 
     +           - (v1 + 1.d0)*(zstar**2/sigmas)
     +             *temp2/temp1**2
                    
      enddo      
      asigma = asigma - 1.d0/sigmad
      
      asigma = -1.d0/asigma
      if (asigma .lt. 0.0d0) asigma = -asigma      

      bpdf = -fdelta(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)
                
      do ii = 1, 20
	  
          call rnset(iseed)
	    rv = drnnof()
	    call rnget(iseed)
		anew = amean + rv*dsqrt(asigma)

		apdf = -fdelta(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
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

      delta = bold
            
      end subroutine
                       

      real*8 function fdelta(delta)      
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1
	
	real*8 bmu, bi(ni,nb)
	real*8 sigmad
      
	real*8 si, xxt(nt), xt, ystar, zstar     
      real*8 summ, pdf
	      
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vsigmas/sigmas
      common /vv1/v1

      common /vbmu/bmu
      common /vbi/bi

      common /vsigmad/sigmad
            
      summ = 0.d0
      do i = 1, ni
              
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
          
          zstar = zeta(i) - ezeta  

          ystar = si - xt - delta*zstar
          
          summ = summ - (v1 + 1.d0)/2.d0
     +                  *dlog(1.d0 + ystar**2/sigmas)
                    
      enddo    
      
      pdf = summ - delta**2/(2.d0*sigmad)
      
      fdelta = -pdf
    
      end function

      
	subroutine gibbs_sigmas(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1
	
	real*8 bmu, bi(ni,nb)

	real*8 a1, b1

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

	real*8 si, xxt(nt), xt, ystar		
      real*8 sumi(ni), summ
	
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1

      common /vbmu/bmu
      common /vbi/bi

      common /va1/a1
      common /vb1/b1

      common /vsumi/sumi

      external fsigmas, dlinrg, drnnof, drnunf

      do i = 1, ni
              
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

          ystar = si - xt - delta*(zeta(i) - ezeta)
                                
          sumi(i) = ystar**2
                    
      enddo

      bold = dlog(sigmas)

      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fsigmas, nopt, start, xmin, ynewlo,  
     +            reqmin, step, konvge, kcount, 
     +            icount, numres, ifault)
	amean = xmin(1)

      summ = 0.d0
      do i = 1, ni
          summ = summ 
     +         - (v1 + 1.d0)/2.d0*sumi(i)*dexp(-amean)
     +           /(1.d0 + sumi(i)*dexp(-amean))**2
      enddo
      
      asigma = summ - b1*dexp(-amean)
	
      asigma = -1.d0/asigma
      if (asigma .lt. 0.0d0) asigma = -asigma      

      bpdf = -fsigmas(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)
                
      do ii = 1, 20
	  
          call rnset(iseed)
	    rv = drnnof()
	    call rnget(iseed)
		anew = amean + rv*dsqrt(asigma)

		apdf = -fsigmas(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
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

      sigmas = dexp(bold)
           
	end subroutine


      real*8 function fsigmas(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 v1
	real*8 a1, b1
	
      real*8 sumi(ni)
      
      real*8 star, summ, pdf
            
      common /vv1/v1

      common /va1/a1
      common /vb1/b1

      common /vsumi/sumi
                  
      summ = 0.d0
      do i = 1, ni
          summ = summ 
     +         - star/2.d0
     +         - (v1 + 1.d0)/2.d0
     +           *dlog(1.d0 + sumi(i)*dexp(-star))
      enddo
      pdf = summ - a1*star - b1*dexp(-star)
                 
      fsigmas = -pdf
    
	end function

      
	subroutine gibbs_v1(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3
	integer, parameter :: npoint = 5, cpoint = 3, ndim = 3

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 zeta(ni), ezeta
	real*8 delta, sigmas, v1
	
	real*8 bmu, bi(ni,nb)

	real*8 a2, b2

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      real*8 epsilon, astar, error
      real*8 yystar(npoint), xxstar(npoint,ndim)
      real*8 xx(ndim,ndim), xy(ndim)
      real*8 xxi(ndim,ndim), xhat(ndim)

	real*8 si, xxt(nt), xt, ystar		
      real*8 suma, sumb
	
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      common /vsigmas/sigmas
      common /vv1/v1

      common /vbmu/bmu
      common /vbi/bi

      common /va2/a2
      common /vb2/b2

      common /vsuma/suma
      common /vsumb/sumb

      external fv1, dlinrg, drnnof, drnunf

      epsilon = 0.1d0

      suma = 0.d0; sumb = 0.d0
      do i = 1, ni
              
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

        ystar = si - xt - delta*(zeta(i) - ezeta)

        suma = suma + 1.d0                                
        sumb = sumb + dlog(1.d0 + ystar**2/sigmas)
                    
      enddo 
      
      bold = dlog(v1)

      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fv1, nopt, start, xmin, ynewlo,  
     +            reqmin, step, konvge, kcount, 
     +            icount, numres, ifault)
	amean = xmin(1)
	
  11  do kk = 1, npoint
                  
          error = dfloat(kk - cpoint)
          astar = amean + error*epsilon
                  
          yystar(kk) = fv1(astar)
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
     +                        + xxstar(kk,kk1)
     +                          *xxstar(kk,kk2)
              enddo
          enddo
          xy(kk1) = 0.d0
          do kk2 = 1, npoint
              xy(kk1) = xy(kk1) + xxstar(kk2,kk1)
     +                            *yystar(kk2)
          enddo
      enddo

      call dlinrg(ndim, xx, ndim, xxi, ndim)	
                        
      do kk1 = 1, ndim
          xhat(kk1) = 0.d0
          do kk2 = 1, ndim
              xhat(kk1) = xhat(kk1)
     +                  + xxi(kk1,kk2)*xy(kk2)
          enddo
      enddo
                                                       
      asigma = 1.0d0/(xhat(1)*2.d0)            
      if (asigma .lt. 0.0d0) asigma = -asigma

      bpdf = -fv1(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)
                
      do ii = 1, 20
	  
          call rnset(iseed)
	    rv = drnnof()
	    call rnget(iseed)
		anew = amean + rv*dsqrt(asigma)

		apdf = -fv1(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
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

      v1 = dexp(bold)
                 
	end subroutine


      real*8 function fv1(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 a2, b2
      real*8 suma, sumb
      
      real*8 star, v1, summ, pdf
            
      common /va2/a2
      common /vb2/b2

      common /vsuma/suma
      common /vsumb/sumb

      external dlngam
            
      v1 = dexp(star)
      
      summ = suma*(dlngam((v1 + 1.d0)/2.d0) - dlngam(v1/2.d0))
     +     - v1/2.d0*sumb

      pdf = summ + a2*star - b2*dexp(star)
                 
      fv1 = -pdf
    
	end function


	subroutine gibbs_lambda(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 ti(ni), ui(ni,nu), theta(nt) 
	real*8 lambda(ni), zeta(ni), ezeta
	real*8 delta, sigmas, v1
	real*8 bmu, bi(ni,nb)

	real*8 si, xxt(nt), xt, ystar		
      real*8 shape, scale, rv
	
      common /vti/ti
      common /vui/ui
      common /vtheta/theta
      
      common /vlambda/lambda
      
      common /vzeta/zeta
      common /vezeta/ezeta
      common /vdelta/delta
      
      common /vsigmas/sigmas
      common /vv1/v1

      common /vbmu/bmu
      common /vbi/bi

      external drngam

      do i = 1, ni
              
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

          ystar = si - xt - delta*(zeta(i) - ezeta)
            
          shape = (v1 + 1.d0)/2.d0
          scale = (1.d0 + ystar**2/sigmas)/2.d0 

          call rnset(iseed)
	    call drngam(1, shape, rv)
	    call rnget(iseed)
	    lambda(i) = rv/scale
                    
      enddo
     
	end subroutine


	subroutine gibbs_phi(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 zi(ni), xi(ni,nx), phi(np,nl) 
	real*8 bi(ni,nb)
	real*8 sigmap

	real*8 start(np), xmin(np), ynewlo, reqmin, step(np)
	integer konvge, kcount, icount, numres, ifault, nopt

      real*8 bold(np), anew(np), amean(np), asigma(np,np)
	real*8 der2(np,np), tol, rsig(np,np), hmean(np)
	real*8 hpdf, bpdf, apdf, ratio, u

      integer ldum
      real*8 xxp(ni,np)
      real*8 suma, sumb, xb
      		
      common /vzi/zi
      common /vxi/xi
      common /vphi/phi

      common /vbi/bi

      common /vsigmap/sigmap

      common /vldum/ldum
      common /vxxp/xxp

      external fphi, dlinrg, dchfac, dblinf, drnmvn, drnunf

      do l = 1, nl-1

        ldum = l
        
        do i = 1, ni
            do j = 1, nx
                xxp(i,j) = xi(i,j)
            enddo
            do j = 1, nb-1
                xxp(i,nx+j) = bi(i,j)
            enddo
        enddo      

        do j = 1, np
            bold(j) = phi(j,l) 
        enddo
      
	  nopt = np
	  reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	  do j = 1, np
	      step(j) = 0.2d0 ; start(j) = bold(j)
	  enddo	
	  call nelmin(fphi, nopt, start, xmin, 
     +              ynewlo, reqmin, step, konvge, 
     +              kcount, icount, numres, ifault)       
        do j = 1, np
            amean(j) = xmin(j)
            phi(j,l) = amean(j)
        enddo

        do j1 = 1, np
            do j2 = 1, np
                der2(j1,j2) = 0.d0
            enddo
        enddo
        do j1 = 1, np
            do j2 = j1, np

                do i = 1, ni
              
                    suma = 0.d0
                    do k = 1, nl-1                    
              
                        xb = 0.d0
                        do j = 1, np
                            xb = xb + xxp(i,j)*phi(j,k)
                        enddo
                        
                        suma = suma + dexp(xb)
                        
                    enddo

                    xb = 0.d0
                    do j = 1, np
                        xb = xb + xxp(i,j)*phi(j,l)
                    enddo
                    
                    sumb = 1.d0 + suma - dexp(xb)
                
                    der2(j1,j2) = der2(j1,j2) 
     +                          - dexp(xb)*sumb
     +                            /(1.d0 + suma)**2
     +                            *xxp(i,j1)*xxp(i,j2)
        
                enddo

            enddo
        enddo

        do j = 1, np
            der2(j,j) = der2(j,j) - 1.d0/sigmap       
        enddo

        do j1 = 1, np
            do j2 = j1, np
                der2(j2,j1) = der2(j1,j2)
            enddo
        enddo                
        do j1 = 1, np
            do j2 = 1, np
                der2(j1,j2) = -der2(j1,j2)
            enddo
        enddo                

	  call dlinrg(np, der2, np, asigma, np)	
	  	  
	  tol = 100.0d0*dmach(4)
	  call dchfac(np, asigma, np, tol, irank, rsig, np)

	  do j = 1, np
	      hmean(j) = bold(j) - amean(j)
	  enddo
	  hpdf = dblinf(np, np, der2, np, hmean, hmean)
	  bpdf = -fphi(bold) + hpdf/2.d0

	  do ii = 1, 20
	
	      call rnset(iseed)
		    call drnmvn(1, np, rsig, np, anew, 1)
		    call rnget(iseed)		
		    do j = 1, np
	          anew(j) = amean(j) + anew(j)
		    enddo
		
		    do j = 1, np
		        hmean(j) = anew(j) - amean(j)
		    enddo
		    hpdf = dblinf(np, np, der2, np, hmean, hmean)
		    apdf = -fphi(anew) + hpdf/2.d0

		    ratio = apdf - bpdf 

		    if (ratio .ge. 0.0d0) then
		        bpdf = apdf
			    do j = 1, np
				    bold(j) = anew(j)
			    enddo
		    else
		        call rnset(iseed)
			    u = drnunf()
			    call rnget(iseed)
			    if (dlog(u) .le. ratio) then
			        bpdf = apdf
				    do j = 1, np
					    bold(j) = anew(j)
				    enddo
			    endif
		    endif
		
	  enddo      

        do j = 1, np
            phi(j,l) = bold(j)
        enddo

      enddo
                                         
	end subroutine


      real*8 function fphi(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 zi(ni), phi(np,nl) 
	real*8 sigmap

      integer ldum
      real*8 xxp(ni,np)
      
      real*8 star(np)
      real*8 suma, sumb, xb
      real*8 sum1, sum2, pdf
      		
      common /vzi/zi
      common /vphi/phi

      common /vsigmap/sigmap

      common /vldum/ldum
      common /vxxp/xxp

      l = ldum
      
      do j = 1, np
        phi(j,l) = star(j)
      enddo 
            
      sum1 = 0.d0
      do i = 1, ni

        suma = 0.d0
        do k = 1, nl-1              
        
            xb = 0.d0
            do j = 1, np
                xb = xb + xxp(i,j)*phi(j,k)
            enddo
            
            suma = suma + dexp(xb)
            
        enddo

        sumb = 0.d0
        if (zi(i) .eq. dfloat(l)) then
        
            sumb = 0.d0
            do j = 1, np
                sumb = sumb + xxp(i,j)*phi(j,l)
            enddo
                        
        endif            
                        
        sum1 = sum1 + sumb - dlog(1.d0 + suma)
        
      enddo
      
      sum2 = 0.d0
      do j = 1, np
        sum2 = sum2 - phi(j,l)**2/(2.d0*sigmap)
      enddo
      
      pdf = sum1 + sum2            

      fphi = -pdf
     
	end function

		
	
	subroutine gibbs_omega(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

	real*8 bi(ni,nb), omega(nb,nb) 
	real*8 d0, s0
	
	real*8 shape, qq(nb,nb), bstar(nb)
		
      common /vbi/bi
      common /vomega/omega

      common /vd0/d0
      common /vs0/s0

	shape = d0	      
      do k1 = 1, nb
        do k2 = 1, nb
		    qq(k1,k2) = 0.d0
        enddo
        qq(k1,k1) = d0*s0
      enddo

      do i = 1, ni            

	  shape = shape + 1.d0

        do k = 1, nb
            bstar(k) = bi(i,k)
        enddo
                                       
        do k1 = 1, nb
            do k2 = 1, nb
                qq(k1,k2) = qq(k1,k2) + bstar(k1)*bstar(k2)
            enddo
        enddo
                            
      enddo
      
	call gen_iwish(omega,shape,qq,nb,iseed)
      
	end subroutine


	subroutine gen_iwish(sigma,shape,qq,nn,iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      real*8 sigma(nn,nn), qq(nn,nn), shape
      
	real*8 qqinv(nn,nn)
	real*8 df, chi(nn), anv(nn,nn), bb(nn,nn)
	real*8 tol, rsig(nn,nn), suma
	real*8 btemp(nn,nn), aainv(nn,nn)
            
	external drnchi, drnnof, dlinrg, 
     +         dchfac, dmach, dmxytf, dmrrrr 
                  
      do j = 1, nn
      
		df = shape - dfloat(j) + 1.d0
		
		call rnset(iseed)
		call drnchi(1, df, chi(j))
		call rnget(iseed)
		
		do i = 2, nn
			call rnset(iseed)		
			anv(j,i) = drnnof()
			call rnget(iseed)
		enddo
		
	enddo

      bb(1,1) = chi(1)
      do j = 2, nn

		suma = 0.0d0
		do i = 1, j-1         
			suma = suma + anv(i,j)**2
		enddo
		bb(j,j) = chi(j) + suma

		bb(1,j) = anv(1,j)*dsqrt(chi(1))
          bb(j,1) = bb(1,j)
        
        if (j .gt. 2) then
		    do i = 2, j-1         
			    suma = 0.0d0
			    do k = 1, i-1
				    suma = suma + anv(k,j)*anv(k,i)
			    enddo
			    bb(i,j) = anv(i,j)*dsqrt(chi(i)) + suma
		        bb(j,i) = bb(i,j)
		    enddo
	  endif
	  
	enddo
      
	call dlinrg(nn, qq, nn, qqinv, nn)	

	tol = 100.0d0*dmach(4)	
	call dchfac(nn, qqinv, nn, tol, irank, rsig, nn)
	
      call dmxtyf(nn,nn,rsig,nn,nn,nn,bb,
     +            nn,nn,nn,btemp,nn)
      call dmrrrr(nn,nn,btemp,nn,nn,nn,rsig,
     +            nn,nn,nn,aainv,nn)
		
	call dlinrg(nn, aainv, nn, sigma, nn)	      
                                                        
      end subroutine      


	subroutine gen_wish(sigma,shape,qq,nn,iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 446, nj = 15, nb = 4
	integer, parameter :: nc = 3, nu = 4, nt = 7
	integer, parameter :: nx = 2, np = 5, nl = 3

      real*8 sigma(nn,nn), qq(nn,nn), shape
      
	real*8 qqinv(nn,nn)
	real*8 df, chi(nn), anv(nn,nn), bb(nn,nn)
	real*8 tol, rsig(nn,nn), suma
	real*8 btemp(nn,nn), aainv(nn,nn)
            
	external drnchi, drnnof, dlinrg, 
     +         dchfac, dmach, dmxytf, dmrrrr 
                  
      do j = 1, nn
      
		df = shape - dfloat(j) + 1.d0
		
		call rnset(iseed)
		call drnchi(1, df, chi(j))
		call rnget(iseed)
		
		do i = 2, nn
			call rnset(iseed)		
			anv(j,i) = drnnof()
			call rnget(iseed)
		enddo
		
	enddo

      bb(1,1) = chi(1)
      do j = 2, nn

		suma = 0.0d0
		do i = 1, j-1         
			suma = suma + anv(i,j)**2
		enddo
		bb(j,j) = chi(j) + suma

		bb(1,j) = anv(1,j)*dsqrt(chi(1))
          bb(j,1) = bb(1,j)
        
          if (j .gt. 2) then
		    do i = 2, j-1         
			    suma = 0.0d0
			    do k = 1, i-1
				    suma = suma + anv(k,j)*anv(k,i)
			    enddo
			    bb(i,j) = anv(i,j)*dsqrt(chi(i)) + suma
		        bb(j,i) = bb(i,j)
		    enddo
	    endif
	  
	enddo
      
	call dlinrg(nn, qq, nn, qqinv, nn)	

	tol = 100.0d0*dmach(4)	
	call dchfac(nn, qqinv, nn, tol, irank, rsig, nn)
	
      call dmxtyf(nn,nn,rsig,nn,nn,nn,bb,
     +            nn,nn,nn,btemp,nn)
      call dmrrrr(nn,nn,btemp,nn,nn,nn,rsig,
     +            nn,nn,nn,sigma,nn)
		                                                        
      end subroutine      
	