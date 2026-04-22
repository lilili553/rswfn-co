      subroutine rks_x_pbe
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew, K. Burke, and M. Ernzerhof
c     Generalized gradient approximation made simple
c     Phys. Rev. Lett. 77 (1996) 3865-3868
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t4 = rho**2
      t5 = t2**2
      zk(i) = -0.7385587663820224D0*t2*rho*(0.1804D1-0.804D0/(1.D0
     #+0.7131826587600489D-2*sigma/t5/t4))
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = rho**2
      t5 = t2**2
      t10 = 1.D0+0.7131826587600489D-2*sigma/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.7385587663820224D0*t3*t13
      t20 = t10**2
      t21 = 1/t20
      vrhoa(i) = -0.9847450218426965D0*t2*t13+0.1129303341188623D-1/t2
     #/t4*t21*sigma
      vsigmaaa(i) = -0.1693955011782934D-1/t3*t21
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = rho**2
      t5 = t2**2
      t10 = 1.D0+0.7131826587600489D-2*sigma/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.7385587663820224D0*t3*t13
      t20 = t10**2
      t21 = 1/t20
      t22 = 1/t2/t4*t21
      vrhoa(i) = -0.9847450218426965D0*t2*t13+0.1129303341188623D-1
     #*t22*sigma
      vsigmaaa(i) = -0.1693955011782934D-1/t3*t21
      t37 = t4**2
      t41 = 1/t20/t10
      t43 = sigma**2
      v2rhoa2(i) = -0.6564966812284644D0/t5*t13-0.2258606682377246D-1
     #/t2/t4/rho*t21*sigma+0.8590928633765426D-3/t37/t4*t41*t43
      v2rhoasigmaaa(i) = 0.2258606682377246D-1*t22
     #-0.644319647532407D-3/t37/rho*t41*sigma
      v2sigmaaa2(i) = 0.9664794712986105D-3/t37*t41
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      
      endif ! ideriv
      return
      end