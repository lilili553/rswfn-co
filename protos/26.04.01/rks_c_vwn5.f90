      subroutine rks_c_vwn5
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     S.H. Vosko, L. Wilk, and M. Nusair
c     Accurate spin-dependent electron liquid correlation energies for
c     local spin density calculations: a critical analysis
c     Can. J. Phys. 58 (1980) 1200-1211
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
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t4 = t1**(1.D0/6.D0)
      t7 = 1/(0.6203504908994D0*t2+0.2935818660072219D1*t4+0.129352D2)
      t10 = dlog(0.6203504908994D0*t2*t7)
      t16 = datan(0.615199081975908D1/(0.1575246635799487D1*t4
     #+0.372744D1))
      t20 = (0.7876233178997433D0*t4+0.10498D0)**2
      t22 = dlog(t20*t7)
      zk(i) = rho*(0.310907D-1*t10+0.3878329487811301D-1*t16
     #+0.9690227711544374D-3*t22)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t4 = t1**(1.D0/6.D0)
      t6 = 0.6203504908994D0*t2+0.2935818660072219D1*t4+0.129352D2
      t7 = 1/t6
      t10 = dlog(0.6203504908994D0*t2*t7)
      t11 = 0.310907D-1*t10
      t13 = 0.1575246635799487D1*t4+0.372744D1
      t16 = datan(0.615199081975908D1/t13)
      t17 = 0.3878329487811301D-1*t16
      t19 = 0.7876233178997433D0*t4+0.10498D0
      t20 = t19**2
      t22 = dlog(t20*t7)
      t23 = 0.9690227711544374D-3*t22
      zk(i) = rho*(t11+t17+t23)
      t25 = t2**2
      t26 = 1/t25
      t28 = rho**2
      t29 = 1/t28
      t32 = t6**2
      t33 = 1/t32
      t37 = t4**2
      t38 = t37**2
      t40 = 1/t38/t4
      t41 = t40*t29
      t43 = -0.2067834969664667D0*t26*t29-0.4893031100120365D0*t41
      t51 = t13**2
      t52 = 1/t51
      vrhoa(i) = t11+t17+t23+rho*(0.5011795824473985D-1*(
     #-0.2067834969664667D0*t26*t7*t29-0.6203504908994D0*t2*t33*t43)
     #/t2*t6+0.626408570946439D-1*t52*t40*t29/(1.D0+0.378469910464D2
     #*t52)+0.9690227711544374D-3*(-0.2625411059665811D0*t19*t7*t41
     #-1.D0*t20*t33*t43)/t20*t6)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t4 = t1**(1.D0/6.D0)
      t6 = 0.6203504908994D0*t2+0.2935818660072219D1*t4+0.129352D2
      t7 = 1/t6
      t10 = dlog(0.6203504908994D0*t2*t7)
      t11 = 0.310907D-1*t10
      t13 = 0.1575246635799487D1*t4+0.372744D1
      t16 = datan(0.615199081975908D1/t13)
      t17 = 0.3878329487811301D-1*t16
      t19 = 0.7876233178997433D0*t4+0.10498D0
      t20 = t19**2
      t22 = dlog(t20*t7)
      t23 = 0.9690227711544374D-3*t22
      zk(i) = rho*(t11+t17+t23)
      t25 = t2**2
      t26 = 1/t25
      t27 = t26*t7
      t28 = rho**2
      t29 = 1/t28
      t32 = t6**2
      t33 = 1/t32
      t34 = t2*t33
      t37 = t4**2
      t38 = t37**2
      t39 = t38*t4
      t40 = 1/t39
      t41 = t40*t29
      t43 = -0.2067834969664667D0*t26*t29-0.4893031100120365D0*t41
      t46 = -0.2067834969664667D0*t27*t29-0.6203504908994D0*t34*t43
      t47 = 1/t2
      t48 = t46*t47
      t49 = t48*t6
      t51 = t13**2
      t52 = 1/t51
      t53 = t52*t40
      t55 = 1.D0+0.378469910464D2*t52
      t56 = 1/t55
      t58 = t53*t29*t56
      t60 = t19*t7
      t63 = t20*t33
      t66 = -0.2625411059665811D0*t60*t41-1.D0*t63*t43
      t67 = 1/t20
      t68 = t66*t67
      t69 = t68*t6
      vrhoa(i) = t11+t17+t23+rho*(0.5011795824473985D-1*t49
     #+0.626408570946439D-1*t58+0.9690227711544374D-3*t69)
      t77 = 1/t25/t1
      t79 = t28**2
      t80 = 1/t79
      t81 = t77*t7*t80
      t88 = 1/t28/rho
      t92 = 1/t32/t6
      t94 = t43**2
      t102 = 1/t39/t1
      t103 = t102*t80
      t105 = t40*t88
      t107 = -0.1378556646443111D0*t77*t80+0.4135669939329333D0*t26
     #*t88-0.4077525916766971D0*t103+0.978606220024073D0*t105
      t125 = t80*t56
      t134 = t51**2
      t138 = t55**2
      s1 = 0.2004718329789594D0*t49+0.2505634283785756D0*t58
      v2rhoa2(i) = s1+0.387609108461775D-2*t69+2.D0*rho*
     #(0.5011795824473985D-1*(-0.1378556646443111D0*t81
     #+0.4135669939329333D0*t26*t33*t29*t43+0.4135669939329333D0*t27
     #*t88+0.12407009817988D1*t2*t92*t94-0.6203504908994D0*t34*t107)
     #*t47*t6+0.1670598608157995D-1*t46/t2/t1*t6*t29
     #+0.5011795824473985D-1*t48*t43+0.3289159980064473D-1/t51/t13*t77
     #*t125+0.5220071424553658D-1*t52*t102*t125-0.1252817141892878D0
     #*t53*t88*t56-0.1244848083156773D1/t134/t13*t77*t80/t138
     #+0.9690227711544374D-3*(0.3446391616107778D-1*t81
     #+0.5250822119331622D0*t19*t33*t41*t43-0.2187842549721509D0*t60
     #*t103+0.5250822119331622D0*t60*t105+2.D0*t20*t92*t94-1.D0*t63
     #*t107)*t67*t6+0.2544083100456872D-3*t66/t20/t19*t6*t40*t29
     #+0.9690227711544374D-3*t68*t43)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end

