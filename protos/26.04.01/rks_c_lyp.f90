      subroutine rks_c_lyp
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     C. Lee, W. Yang, and R.G. Parr
c     Development of the Colle-Salvetti correlation-energy formula into
c     a functional of the electron density
c     Phys. Rev. B37 (1988) 785-789
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
      t3 = 1/t2
      t6 = 1/(1.D0+0.349D0*t3)
      t9 = 0.2533D0*t3
      t10 = dexp(-t9)
      t12 = rho**2
      t14 = t2**2
      t20 = t3*t6
      zk(i) = -0.4918D-1*t6*rho-0.649176D-2*t10*t6/t14/t12/rho*(0.25D0
     #*t12*(0.1148493600075277D2*t14*t12+(0.2611111111111111D1
     #-0.9850555555555556D-1*t3-0.1357222222222222D0*t20)*sigma-0.5D0*
     #(0.25D1-0.1407222222222222D-1*t3-0.1938888888888889D-1*t20)
     #*sigma-0.2777777777777778D-1*(t9+0.349D0*t20-11.D0)*sigma)
     #-0.4583333333333333D0*t12*sigma)
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
      t3 = 1/t2
      t5 = 1.D0+0.349D0*t3
      t6 = 1/t5
      t9 = 0.2533D0*t3
      t10 = dexp(-t9)
      t11 = t10*t6
      t12 = rho**2
      t14 = t2**2
      t16 = 1/t14/t12/rho
      t20 = t3*t6
      t22 = 0.2611111111111111D1-0.9850555555555556D-1*t3
     #-0.1357222222222222D0*t20
      t30 = t9+0.349D0*t20-11.D0
      t33 = 0.1148493600075277D2*t14*t12+t22*sigma-0.5D0*(0.25D1
     #-0.1407222222222222D-1*t3-0.1938888888888889D-1*t20)*sigma
     #-0.2777777777777778D-1*t30*sigma
      t38 = 0.25D0*t12*t33-0.4583333333333333D0*t12*sigma
      zk(i) = -0.4918D-1*t6*rho-0.649176D-2*t11*t16*t38
      t45 = t14*rho
      t49 = t30/rho*sigma
      t54 = rho*sigma
      t60 = t5**2
      t61 = 1/t60
      t64 = t12**2
      t66 = 1/t64/rho
      t81 = 1/t2/rho
      t83 = t81*t6
      t85 = 1/t45
      t86 = t85*t61
      vrhoa(i) = -0.4918D-1*t6-0.649176D-2*t11*t16*(0.5D0*rho*t33
     #+0.25D0*t12*(0.3062649600200738D2*t45-0.2777777777777778D-1*t49)
     #-0.25D0*t54)-0.5721273333333333D-2*t61*t3-0.548120936D-3*t66*t10
     #*t6*t38-0.75520808D-3*t10*t61*t66*t38+0.2380312D-1*t11/t14/t64
     #*t38-0.649176D-2*t11*t16*(0.25D0*t12*((0.3283518518518519D-1*t81
     #+0.4524074074074074D-1*t83-0.1578901851851852D-1*t86)*sigma
     #-0.5D0*(0.4690740740740741D-2*t81+0.6462962962962963D-2*t83
     #-0.2255574074074074D-2*t86)*sigma-0.2777777777777778D-1*(
     #-0.8443333333333333D-1*t81-0.1163333333333333D0*t83
     #+0.4060033333333333D-1*t86)*sigma+0.2777777777777778D-1*t49)
     #-0.6666666666666667D0*t54)
      vsigmaaa(i) = 0.7213066666666667D-3*t11*t85-0.2596704D-1*t11*t16
     #*(0.25D0*t12*t22-0.6666666666666667D0*t12)
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
      t3 = 1/t2
      t5 = 1.D0+0.349D0*t3
      t6 = 1/t5
      t9 = 0.2533D0*t3
      t10 = dexp(-t9)
      t11 = t10*t6
      t12 = rho**2
      t13 = t12*rho
      t14 = t2**2
      t16 = 1/t14/t13
      t17 = t14*t12
      t18 = 0.1148493600075277D2*t17
      t20 = t3*t6
      t22 = 0.2611111111111111D1-0.9850555555555556D-1*t3
     #-0.1357222222222222D0*t20
      t23 = t22*sigma
      t28 = 0.5D0*(0.25D1-0.1407222222222222D-1*t3
     #-0.1938888888888889D-1*t20)*sigma
      t30 = t9+0.349D0*t20-11.D0
      t32 = 0.2777777777777778D-1*t30*sigma
      t33 = t18+t23-t28-t32
      t38 = 0.25D0*t12*t33-0.4583333333333333D0*t12*sigma
      zk(i) = -0.4918D-1*t6*rho-0.649176D-2*t11*t16*t38
      t45 = t14*rho
      t47 = 1/rho
      t48 = t30*t47
      t49 = t48*sigma
      t51 = 0.3062649600200738D2*t45-0.2777777777777778D-1*t49
      t54 = rho*sigma
      t56 = 0.5D0*rho*t33+0.25D0*t12*t51-0.25D0*t54
      t60 = t5**2
      t61 = 1/t60
      t64 = t12**2
      t65 = t64*rho
      t66 = 1/t65
      t67 = t66*t10
      t68 = t6*t38
      t71 = t10*t61
      t76 = 1/t14/t64
      t81 = 1/t2/rho
      t83 = t81*t6
      t85 = 1/t45
      t86 = t85*t61
      t88 = 0.3283518518518519D-1*t81+0.4524074074074074D-1*t83
     #-0.1578901851851852D-1*t86
      t99 = -0.8443333333333333D-1*t81-0.1163333333333333D0*t83
     #+0.4060033333333333D-1*t86
      t103 = t88*sigma-0.5D0*(0.4690740740740741D-2*t81
     #+0.6462962962962963D-2*t83-0.2255574074074074D-2*t86)*sigma
     #-0.2777777777777778D-1*t99*sigma+0.2777777777777778D-1*t49
      t107 = 0.25D0*t12*t103-0.6666666666666667D0*t54
      vrhoa(i) = -0.4918D-1*t6-0.649176D-2*t11*t16*t56
     #-0.5721273333333333D-2*t61*t3-0.548120936D-3*t67*t68
     #-0.75520808D-3*t71*t66*t38+0.2380312D-1*t11*t76*t38-0.649176D-2
     #*t11*t16*t107
      t116 = 0.25D0*t12*t22-0.6666666666666667D0*t12
      vsigmaaa(i) = 0.7213066666666667D-3*t11*t85-0.2596704D-1*t11*t16
     #*t116
      t121 = 1/t60/t5
      t126 = t64*t12
      t127 = 1/t126
      t137 = t99*t47*sigma
      t141 = t30/t12*sigma
      t150 = rho*t51
      t165 = 1/t2/t126
      t166 = t165*t10
      t188 = 1/t2/t12
      t190 = t188*t6
      t192 = 1/t17
      t193 = t192*t61
      t195 = 1/t13
      t196 = t195*t121
      s1 = -0.2662299191111111D-2*t121*t85-0.7628364444444444D-2*t61
     #*t81+0.9500762890666667D-2*t127*t10*t68+0.9521248D-1*t11*t76*t56
     #-0.2596704D-1*t11*t16*(0.5D0*rho*t103+0.25D0*t12*(
     #-0.2777777777777778D-1*t137+0.2777777777777778D-1*t141))
     #-0.649176D-2*t11*t16*(t18+t23-t28-t32+t150)+0.9521248D-1*t11*t76
     #*t107-0.2192483744D-2*t67*t6*t56-0.302083232D-2*t71*t66*t56
      s2 = s1-0.2550589422186667D-3*t166*t61*t38-0.302083232D-2*t71
     #*t66*t107-0.3514234932266667D-3*t10*t121*t165*t38
     #-0.2192483744D-2*t67*t6*t107
      v2rhoa2(i) = s2-0.9255935539253333D-4*t166*t68
     #-0.2221624533333333D0*t11/t14/t65*t38-0.1298352D-1*t11*t16*
     #(0.25D0*t12*((-0.4378024691358025D-1*t188-0.6032098765432099D-1
     #*t190+0.3157803703703704D-1*t193-0.3673578308641975D-2*t196)
     #*sigma-0.5D0*(-0.6254320987654321D-2*t188-0.8617283950617284D-2
     #*t190+0.4511148148148148D-2*t193-0.5247969012345679D-3*t196)
     #*sigma-0.2777777777777778D-1*(0.1125777777777778D0*t188
     #+0.1551111111111111D0*t190-0.8120066666666667D-1*t193
     #+0.9446344222222222D-2*t196)*sigma+0.5555555555555556D-1*t137
     #-0.5555555555555556D-1*t141)-0.6666666666666667D0*sigma)
     #-0.649176D-2*t11*t16*(0.1D1*t150+0.2552208000167282D2*t17-0.5D0
     #*sigma)+0.1309027338666667D-1*t71*t127*t38
      v2rhoasigmaaa(i) = -0.649176D-2*t11*t16*(-0.9444444444444444D0
     #*rho-0.2777777777777778D-1*rho*t30)+0.6090232622222222D-4*t195
     #*t10*t6+0.8391200888888889D-4*t71*t195+0.9978075555555556D-2*t11
     #*t192-0.1298352D-1*t11*t16*(0.25D0*t12*(-0.3D-22*t83+0.1D-22*t86
     #+0.5555555555555556D-1*t48)+0.1333333333333333D1*rho)
     #-0.1298352D-1*t11*t192*t22-0.2192483744D-2*t67*t6*t116
     #-0.302083232D-2*t71*t66*t116+0.9521248D-1*t11*t76*t116
     #-0.2596704D-1*t11*t16*(0.25D0*t12*t88-0.1333333333333333D1*rho)
      v2sigmaaa2(i) = 0.D0
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
      
