      subroutine rks_c_p86
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew
c     Density-functional approximation for the correlation energy of
c     the inhomogeneous electron gas
c     Phys. Rev. B33 (1986) 8822-8824
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
      logical t5
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t6 = t2**(1.D0/6.D0)
      t12 = dlog(t4)
      t18 = piecewise(1.D0 .le. t4,-0.1423D0/(1.D0
     #+0.8292885914166397D0*t6+0.20682485366586D0*t3),0.311D-1*t12
     #-0.48D-1+0.12407009817988D-2*t3*t12-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t31 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t3
     #+0.2843543831490386D-5*t22)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t22+0.1763993811759022D-1*t2)
      t34 = rho**(1.D0/6.D0)
      t39 = dexp(-0.81290825D-3*t20/t31/t34/rho)
      t41 = rho**(1.D0/3.D0)
      zk(i) = rho*t18+0.1D1*t39*t31*sigma/t41/rho
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t5 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.8292885914166397D0*t6+0.20682485366586D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t5,-0.1423D0/t9,0.311D-1*t12-0.48D-1
     #+0.12407009817988D-2*t3*t12-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t33 = t20/t31
      t34 = rho**(1.D0/6.D0)
      t36 = 1/t34/rho
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rho**(1.D0/3.D0)
      t43 = 1/t41/rho
      t44 = sigma*t43
      zk(i) = rho*t18+0.1D1*t40*t44
      t47 = piecewise(t5,0.D0,0.D0)
      t49 = t9**2
      t51 = t6**2
      t52 = t51**2
      t55 = rho**2
      t56 = 1/t55
      t59 = 1/t22
      t60 = t59*t56
      t73 = piecewise(t5,0.1423D0/t49*(-0.1382147652361066D0/t52/t6
     #*t56-0.6894161788861999D-1*t60),-0.1036666666666667D-1*t2
     #-0.4135669939329333D-3*t59*t12*t56-0.4135669939329333D-3*t3*t2
     #+0.2398688564811013D-2*t60)
      t75 = t31**2
      t80 = 1/t3*t56
      t84 = t28**2
      t93 = (-0.4811024840421814D-2*t60-0.1895695887660258D-5*t80)*t29
     #-1.D0*t24/t84*(-0.1803772444038489D1*t60-0.1210946621972718D0
     #*t80-0.1763993811759022D-1*t56)
      vrhoa(i) = rho*t47+t18+rho*t73+0.1D1*(0.81290825D-3*t20/t75*t36
     #*t93+0.9483929583333333D-3*t33/t34/t55)*t39*t31*sigma*t43+0.1D1
     #*t39*t93*t44-0.1333333333333333D1*t40*sigma/t41/t55
      t115 = dsqrt(rho)
      vsigmaaa(i) = -0.16258165D-2*t20/t115/t55*t39+0.4D1*t40*t43
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
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t5 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.8292885914166397D0*t6+0.20682485366586D0*t3
      t11 = 0.1423D0/t9
      t12 = dlog(t4)
      t14 = t3*t12
      t18 = piecewise(t5,-t11,0.311D-1*t12-0.48D-1+0.12407009817988D-2
     #*t14-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t32 = 1/t31
      t33 = t20*t32
      t34 = rho**(1.D0/6.D0)
      t36 = 1/t34/rho
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rho**(1.D0/3.D0)
      t43 = 1/t41/rho
      t44 = sigma*t43
      zk(i) = rho*t18+0.1D1*t40*t44
      t47 = piecewise(t5,0.D0,0.D0)
      t48 = rho*t47
      t49 = t9**2
      t50 = 1/t49
      t51 = t6**2
      t52 = t51**2
      t53 = t52*t6
      t54 = 1/t53
      t55 = rho**2
      t56 = 1/t55
      t59 = 1/t22
      t60 = t59*t56
      t62 = -0.1382147652361066D0*t54*t56-0.6894161788861999D-1*t60
      t66 = t59*t12
      t69 = t3*t2
      t73 = piecewise(t5,0.1423D0*t50*t62,-0.1036666666666667D-1*t2
     #-0.4135669939329333D-3*t66*t56-0.4135669939329333D-3*t69
     #+0.2398688564811013D-2*t60)
      t75 = t31**2
      t76 = 1/t75
      t77 = t20*t76
      t79 = 1/t3
      t80 = t79*t56
      t82 = -0.4811024840421814D-2*t60-0.1895695887660258D-5*t80
      t84 = t28**2
      t85 = 1/t84
      t86 = t24*t85
      t90 = -0.1803772444038489D1*t60-0.1210946621972718D0*t80
     #-0.1763993811759022D-1*t56
      t93 = t82*t29-1.D0*t86*t90
      t94 = t36*t93
      t98 = 1/t34/t55
      t101 = 0.81290825D-3*t77*t94+0.9483929583333333D-3*t33*t98
      t102 = t101*t39
      t103 = t31*sigma
      t104 = t103*t43
      t107 = t39*t93
      t111 = 1/t41/t55
      t112 = sigma*t111
      vrhoa(i) = t48+t18+rho*t73+0.1D1*t102*t104+0.1D1*t107*t44
     #-0.1333333333333333D1*t40*t112
      t115 = dsqrt(rho)
      t117 = 1/t115/t55
      vsigmaaa(i) = -0.16258165D-2*t20*t117*t39+0.4D1*t40*t43
      t129 = (-0.843D-1/(1.D0+0.1101176160755631D1*t6
     #+0.1619735131738333D0*t3)+t11)*t56
      t135 = (-0.1555D-1*t12+0.211D-1-0.80645563816922D-3*t14
     #+0.421838333811592D-2*t3)*t56
      t137 = piecewise(t5,0.1709920934161366D1
     #*t129,0.1709920934161366D1*t135)
      t139 = t55*rho
      t147 = t62**2
      t152 = t55**2
      t153 = 1/t152
      t156 = 1/t139
      t160 = 1/t22/t2
      t161 = t160*t153
      t163 = t59*t156
      t180 = piecewise(t5,-0.2846D0/t49/t9*t147+0.1423D0*t50*(
     #-0.1151789710300888D0/t53/t2*t153+0.2764295304722132D0*t54*t156
     #-0.4596107859241333D-1*t161+0.13788323577724D0*t163
     #),0.1036666666666667D-1*t56-0.2757113292886222D-3*t160*t12*t153
     #-0.4521665800333405D-2*t163+0.8271339878658667D-3*t66*t156
     #+0.4135669939329333D-3*t3*t56+0.1599125709874009D-2*t161)
      t185 = piecewise(t5,-0.1709920934161366D1*t129,
     #-0.1709920934161366D1*t135)
      t193 = t93**2
      t203 = 1/t69*t153
      t205 = t79*t156
      t215 = t90**2
      t226 = (-0.3207349893614542D-2*t161+0.9622049680843627D-2*t163
     #-0.6318986292200858D-6*t203+0.3791391775320515D-5*t205)*t29-2.D0
     #*t82*t85*t90+2.D0*t24/t84/t28*t215-1.D0*t86*(
     #-0.1202514962692326D1*t161+0.3607544888076978D1*t163
     #-0.4036488739909061D-1*t203+0.2421893243945437D0*t205
     #+0.3527987623518044D-1*t156)
      t247 = t101**2
      v2rhoa2(i) = rho*t137+0.6222222222222222D1*t40*sigma/t41/t139
     #+2.D0*rho*t180+rho*t185-0.5333333333333333D1*t102*t103*t111
     #+0.2D1*(-0.16258165D-2*t20/t75/t31*t36*t193
     #-0.1896785916666667D-2*t77*t98*t93+0.81290825D-3*t77*t36*t226
     #-0.2054851409722222D-2*t33/t34/t139)*t39*t104+0.4D1*t102*t93
     #*sigma*t43+0.2D1*t39*t226*t44-0.5333333333333333D1*t107*t112
     #+0.2D1*t247*t39*t104+4.D0*t48+4.D0*t47+4.D0*t73
      t254 = 1/t20
      t266 = t117*t39
      v2rhoasigmaaa(i) = 0.4D1*(0.406454125D-3*t254*t76*t94
     #+0.4741964791666667D-3*t254*t32*t98)*t39*t104-0.16258165D-2*t101
     #*t20*t266+0.4D1*t102*t31*t43-0.16258165D-2*t33*t266*t93+0.4D1
     #*t107*t43+0.2167755333333333D-2*t20/t115/t139*t39
     #-0.5333333333333333D1*t40*t111
      t287 = t41**2
      v2sigmaaa2(i) = -0.9754899D-2*t254*t117*t39+0.264327929167225D-5
     #/t287/t139*t32*t39
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

