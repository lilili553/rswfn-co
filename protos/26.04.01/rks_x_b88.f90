      subroutine rks_x_b88
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     A.D. Becke
c     Density-functional exchange-energy approximation with correct 
c     asymptotic behaviour
c     Phys. Rev. A38 (1988) 3098-3100
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
      t3 = t2*rho
      t5 = 1/t3
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      zk(i) = -0.7385587663820224D0*t3-0.5291668409558467D-2*t5*sigma/
     #(1.D0+0.317500104573508D-1*t8*t10)
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
      t5 = 1/t3
      t6 = t5*sigma
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      t13 = 1.D0+0.317500104573508D-1*t8*t10
      t14 = 1/t13
      zk(i) = -0.7385587663820224D0*t3-0.5291668409558467D-2*t6*t14
      t18 = rho**2
      t20 = 1/t2/t18
      t24 = t13**2
      t25 = 1/t24
      t30 = t2**2
      t35 = 1/t30/t18
      t39 = dsqrt(1.D0+0.1587401051968199D1*sigma*t35)
      t40 = 1/t39
      vrhoa(i) = -0.9847450218426965D0*t2+0.705555787941129D-2*t20
     #*sigma*t14+0.2645834204779234D-2*t6*t25*(-0.8466669455293548D-1
     #*t7*t20*t10-0.106673350692263D0*sigma/t30/t18/rho*t40)
      vsigmaaa(i) = -0.2116667363823387D-1*t5*t14
     #+0.5291668409558467D-2*t6*t25*(0.6350002091470161D-1/t7*t5*t10
     #+0.8000501301919725D-1*t35*t40)
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
      t5 = 1/t3
      t6 = t5*sigma
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      t13 = 1.D0+0.317500104573508D-1*t8*t10
      t14 = 1/t13
      zk(i) = -0.7385587663820224D0*t3-0.5291668409558467D-2*t6*t14
      t18 = rho**2
      t20 = 1/t2/t18
      t21 = t20*sigma
      t24 = t13**2
      t25 = 1/t24
      t29 = t18*rho
      t30 = t2**2
      t32 = 1/t30/t29
      t35 = 1/t30/t18
      t38 = 1.D0+0.1587401051968199D1*sigma*t35
      t39 = dsqrt(t38)
      t40 = 1/t39
      t43 = -0.8466669455293548D-1*t7*t20*t10-0.106673350692263D0
     #*sigma*t32*t40
      t44 = t25*t43
      vrhoa(i) = -0.9847450218426965D0*t2+0.705555787941129D-2*t21*t14
     #+0.2645834204779234D-2*t6*t44
      t49 = 1/t7
      t55 = 0.6350002091470161D-1*t49*t5*t10+0.8000501301919725D-1*t35
     #*t40
      t56 = t25*t55
      vsigmaaa(i) = -0.2116667363823387D-1*t5*t14
     #+0.5291668409558467D-2*t6*t56
      t62 = 1/t2/t29
      t69 = 1/t24/t13
      t70 = t43**2
      t77 = t18**2
      t83 = sigma**2
      t89 = 1/t39/t38
      v2rhoa2(i) = -0.6564966812284644D0/t30-0.3292593677058602D-1*t62
     #*sigma*t14-0.1411111575882258D-1*t21*t44-0.5291668409558467D-2
     #*t6*t69*t70+0.2645834204779234D-2*t6*t25*(0.3951112412470322D0
     #*t7*t62*t10+0.106673350692263D1*sigma/t30/t77*t40
     #-0.4515557042823225D0*t83/t2/t77/t29*t89)
      t100 = t5*t25
      v2rhoasigmaaa(i) = 0.2822223151764516D-1*t20*t14
     #-0.705555787941129D-2*t21*t56+0.1058333681911693D-1*t100*t43
     #-0.5291668409558467D-2*t6*t69*t43*t55+0.2645834204779234D-2*t6
     #*t25*(-0.169333389105871D0*t49*t20*t10-0.640040104153578D0*t32
     #*t40+0.3386667782117419D0*sigma/t2/t77/t18*t89)
      t124 = t55**2
      v2sigmaaa2(i) = 0.4233334727646774D-1*t100*t55
     #-0.1058333681911693D-1*t6*t69*t124+0.5291668409558467D-2*t6*t25*
     #(-0.1270000418294032D0/t7/sigma*t5*t10+0.1600100260383945D0
     #/sigma*t35*t40-0.2540000836588064D0/t2/t77/rho*t89)
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

