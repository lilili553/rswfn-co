      subroutine rks_x_pw91
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
c     M.R. Pederson, D.J. Singh, C. Fiolhais
c     Atoms, molecules, solids and surfaces:
c     Applications of the generalized gradient approximation
c     for exchange and correlation
c     Phys. Rev. B 46 (1992) 6671--6687
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
      t4 = dsqrt(sigma)
      t6 = t4/t3
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t17 = dexp(-0.261211729852336D1*sigma*t14)
      t25 = sigma**2
      t26 = t11**2
      zk(i) = -0.7385587663820224D0*t3*(1.D0+t10+0.261211729852336D-1*
     #(0.2743D0-0.1508D0*t17)*sigma*t14)/(1.D0+t10
     #+0.272926271249799D-5*t25/t2/t26/rho)
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
      t4 = dsqrt(sigma)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigma*t14
      t17 = dexp(-0.261211729852336D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigma
      t23 = 1.D0+t10+0.261211729852336D-1*t20*t14
      t24 = t3*t23
      t25 = sigma**2
      t26 = t11**2
      t29 = 1/t2/t26/rho
      t32 = 1.D0+t10+0.272926271249799D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.7385587663820224D0*t24*t33
      t43 = 0.8466757147455043D-1*t4/t2/t11*t8
      t46 = 1/t12/t11/rho
      t50 = dsqrt(1.D0+0.1587419712813815D1*t15)
      t51 = 1/t50
      t53 = 0.1066750825533289D0*sigma*t46*t51
      t57 = t25/t2/t26/t11
      t66 = t32**2
      t67 = 1/t66
      vrhoa(i) = -0.9847450218426965D0*t2*t23*t33-0.3692793831910112D0
     #*t3*(-t43-t53-0.5487637560595959D-1*t57*t17-0.1393129225879125D0
     #*t20*t46)*t33+0.3692793831910112D0*t24*t67*(-t43-t53
     #-0.2911213559997856D-4*t57)
      t76 = 0.6350067860591282D-1/t4*t5*t8
      t78 = 0.8000631191499664D-1*t14*t51
      vsigmaaa(i) = -0.7385587663820224D0*t3*(t76+t78
     #+0.411572817044697D-1*t29*t17*sigma+0.1044846919409344D0*t19*t14
     #)*t33+0.7385587663820224D0*t24*t67*(t76+t78
     #+0.2183410169998392D-4*sigma*t29)
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
      t4 = dsqrt(sigma)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigma*t14
      t17 = dexp(-0.261211729852336D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigma
      t23 = 1.D0+t10+0.261211729852336D-1*t20*t14
      t24 = t3*t23
      t25 = sigma**2
      t26 = t11**2
      t29 = 1/t2/t26/rho
      t32 = 1.D0+t10+0.272926271249799D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.7385587663820224D0*t24*t33
      t36 = t2*t23
      t40 = 1/t2/t11
      t43 = 0.8466757147455043D-1*t4*t40*t8
      t44 = t11*rho
      t46 = 1/t12/t44
      t49 = 1.D0+0.1587419712813815D1*t15
      t50 = dsqrt(t49)
      t51 = 1/t50
      t53 = 0.1066750825533289D0*sigma*t46*t51
      t56 = 1/t2/t26/t11
      t57 = t25*t56
      t62 = -t43-t53-0.5487637560595959D-1*t57*t17
     #-0.1393129225879125D0*t20*t46
      t63 = t3*t62
      t66 = t32**2
      t67 = 1/t66
      t69 = -t43-t53-0.2911213559997856D-4*t57
      t70 = t67*t69
      vrhoa(i) = -0.9847450218426965D0*t36*t33-0.3692793831910112D0
     #*t63*t33+0.3692793831910112D0*t24*t70
      t73 = 1/t4
      t76 = 0.6350067860591282D-1*t73*t5*t8
      t78 = 0.8000631191499664D-1*t14*t51
      t79 = t29*t17
      t84 = t76+t78+0.411572817044697D-1*t79*sigma
     #+0.1044846919409344D0*t19*t14
      t85 = t3*t84
      t90 = t76+t78+0.2183410169998392D-4*sigma*t29
      t91 = t67*t90
      vsigmaaa(i) = -0.7385587663820224D0*t85*t33+0.7385587663820224D0
     #*t24*t91
      t107 = 0.395115333547902D0*t4/t2/t44*t8
      t109 = 1/t12/t26
      t112 = 0.1066750825533289D1*sigma*t109*t51
      t116 = t25/t2/t26/t44
      t118 = 1/t50/t49
      t120 = 0.4515683437631874D0*t116*t118
      t124 = t26**2
      t139 = 1/t66/t32
      t140 = t69**2
      v2rhoa2(i) = -0.6564966812284644D0/t12*t23*t33
     #-0.1969490043685393D1*t2*t62*t33+0.1969490043685393D1*t36*t70
     #-0.3692793831910112D0*t3*(t107+t112-t120+0.9877747609072727D0
     #*t116*t17-0.764498826669826D0*t25*sigma/t124/t11*t17
     #+0.1021628098978025D1*t20*t109)*t33+0.7385587663820224D0*t63*t70
     #-0.7385587663820224D0*t24*t139*t140+0.3692793831910112D0*t24*t67
     #*(t107+t112-t120+0.3687537175997285D-3*t116)
      t156 = 0.1693351429491009D0*t73*t40*t8
      t158 = 0.6400504953199731D0*t46*t51
      t159 = sigma*t56
      t161 = 0.3386762578223905D0*t159*t118
      v2rhoasigmaaa(i) = -0.9847450218426965D0*t2*t84*t33
     #+0.9847450218426965D0*t36*t91-0.3692793831910112D0*t3*(-t156
     #-t158+t161-0.6585165072715151D0*t56*t17*sigma
     #+0.5733741200023695D0*t25/t124/rho*t17-0.5572516903516501D0*t19
     #*t46)*t33+0.3692793831910112D0*t63*t91+0.3692793831910112D0*t85
     #*t70-0.7385587663820224D0*t24*t139*t69*t90+0.3692793831910112D0
     #*t24*t67*(-t156-t158+t161-0.2328970847998285D-3*t159)
      t193 = 0.1270013572118256D0/t4/sigma*t5*t8
      t197 = 0.1600126238299933D0/sigma*t14*t51
      t199 = 0.2540071933667929D0*t29*t118
      t211 = t90**2
      v2sigmaaa2(i) = -0.7385587663820224D0*t3*(-t193+t197-t199
     #-0.4300305900017772D0/t124*t17*sigma+0.3292582536357576D0*t79)
     #*t33+0.1477117532764045D1*t85*t91-0.1477117532764045D1*t24*t139
     #*t211+0.7385587663820224D0*t24*t67*(-t193+t197-t199
     #+0.8733640679993569D-4*t29)
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

      