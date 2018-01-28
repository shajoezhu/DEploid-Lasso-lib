      subroutine elnet2(beta,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,ulam,thr,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(nlam),xv(ni)
      real cl(2,ni)
      integer ju(ni),ia(nx),kin(nlam)
      real, dimension (:), allocatable :: a,g
      integer, dimension (:), allocatable :: mm,ix
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)
      allocate(a(1:ni),stat=jerr)
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(g(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(ix(1:ni),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      bta=beta
      omb=1.0-bta
      ix=0
      if(flmin .ge. 1.0)goto 10771
      eqs=max(eps,flmin)
      alf=eqs**(1.0/(nlam-1))
10771 continue
      rsq=0.0
      a=0.0
      mm=0
      nlp=0
      nin=nlp                                                              1035
      iz=0                                                                 1035
      mnl=min(mnlam,nlam)                                                  1035
      alm=0.0                                                              1036
10780 do 10781 j=1,ni                                                      1036
      if(ju(j).eq.0)goto 10781                                             1036
      g(j)=abs(dot_product(y,x(:,j)))                                      1036
10781 continue                                                             1037
10782 continue                                                             1037
10790 do 10791 m=1,nlam                                                    1037
      alm0=alm                                                             1038
      if(flmin .lt. 1.0)goto 10811                                         1038
      alm=ulam(m)                                                          1038
      goto 10801                                                           1039
10811 if(m .le. 2)goto 10821                                               1039
      alm=alm*alf                                                          1039
      goto 10801                                                           1040
10821 if(m .ne. 1)goto 10831                                               1040
      alm=big                                                              1040
      goto 10841                                                           1041
10831 continue                                                             1041
      alm0=0.0                                                             1042
10850 do 10851 j=1,ni                                                      1042
      if(ju(j).eq.0)goto 10851                                             1042
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1042
10851 continue                                                             1043
10852 continue                                                             1043
      alm0=alm0/max(bta,1.0e-3)                                            1043
      alm=alf*alm0                                                         1044
10841 continue                                                             1045
10801 continue                                                             1045
      dem=alm*omb                                                          1045
      ab=alm*bta                                                           1045
      rsq0=rsq                                                             1045
      jz=1                                                                 1046
      tlam=bta*(2.0*alm-alm0)                                              1047
10860 do 10861 k=1,ni                                                      1047
      if(ix(k).eq.1)goto 10861                                             1047
      if(ju(k).eq.0)goto 10861                                             1048
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       1049
10861 continue                                                             1050
10862 continue                                                             1050
10870 continue                                                             1050
10871 continue                                                             1050
      if(iz*jz.ne.0) go to 10360                                           1051
10880 continue                                                             1051
      nlp=nlp+1                                                            1051
      dlx=0.0                                                              1052
10890 do 10891 k=1,ni                                                      1052
      if(ix(k).eq.0)goto 10891                                             1052
      gk=dot_product(y,x(:,k))                                             1053
      ak=a(k)                                                              1053
      u=gk+ak*xv(k)                                                        1053
      v=abs(u)-vp(k)*ab                                                    1053
      a(k)=0.0                                                             1055
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)))
      if(a(k).eq.ak)goto 10891                                             1057
      if(mm(k) .ne. 0)goto 10911                                           1057
      nin=nin+1                                                            1057
      if(nin.gt.nx)goto 10892                                              1058
      mm(k)=nin                                                            1058
      ia(nin)=k                                                            1059
10911 continue                                                             1060
      del=a(k)-ak                                                          1060
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1061
      y=y-del*x(:,k)                                                       1061
      dlx=max(xv(k)*del**2,dlx)                                            1062
10891 continue                                                             1063
10892 continue                                                             1063
      if(nin.gt.nx)goto 10872                                              1064
      if(dlx .ge. thr)goto 10931                                           1064
      ixx=0                                                                1065
10940 do 10941 k=1,ni                                                      1065
      if(ix(k).eq.1)goto 10941                                             1065
      if(ju(k).eq.0)goto 10941                                             1066
      g(k)=abs(dot_product(y,x(:,k)))                                      1067
      if(g(k) .le. ab*vp(k))goto 10961                                     1067
      ix(k)=1                                                              1067
      ixx=1                                                                1067
10961 continue                                                             1068
10941 continue                                                             1069
10942 continue                                                             1069
      if(ixx.eq.1) go to 10880                                             1070
      goto 10872                                                           1071
10931 continue                                                             1072
      if(nlp .le. maxit)goto 10981                                         1072
      jerr=-m                                                              1072
      return                                                               1072
10981 continue                                                             1073
10360 continue                                                             1073
      iz=1                                                                 1074
10990 continue                                                             1074
10991 continue                                                             1074
      nlp=nlp+1                                                            1074
      dlx=0.0                                                              1075
11000 do 11001 l=1,nin                                                     1075
      k=ia(l)                                                              1075
      gk=dot_product(y,x(:,k))                                             1076
      ak=a(k)                                                              1076
      u=gk+ak*xv(k)                                                        1076
      v=abs(u)-vp(k)*ab                                                    1076
      a(k)=0.0                                                             1078
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)))
      if(a(k).eq.ak)goto 11001                                             1080
      del=a(k)-ak                                                          1080
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1081
      y=y-del*x(:,k)                                                       1081
      dlx=max(xv(k)*del**2,dlx)                                            1082
11001 continue                                                             1083
11002 continue                                                             1083
      if(dlx.lt.thr)goto 10992                                             1083
      if(nlp .le. maxit)goto 11021                                         1083
      jerr=-m                                                              1083
      return                                                               1083
11021 continue                                                             1084
      goto 10991                                                           1085
10992 continue                                                             1085
      jz=0                                                                 1086
      goto 10871                                                           1087
10872 continue                                                             1087
      if(nin .le. nx)goto 11041                                            1087
      jerr=-10000-m                                                        1087
      goto 10792                                                           1087
11041 continue                                                             1088
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1088
      kin(m)=nin                                                           1089
      rsqo(m)=rsq                                                          1089
      almo(m)=alm                                                          1089
      lmu=m                                                                1090
      if(m.lt.mnl)goto 10791                                               1090
      if(flmin.ge.1.0)goto 10791                                           1091
      me=0                                                                 1091
11050 do 11051 j=1,nin                                                     1091
      if(ao(j,m).ne.0.0) me=me+1                                           1091
11051 continue                                                             1091
11052 continue                                                             1091
      if(me.gt.ne)goto 10792                                               1092
      if(rsq-rsq0.lt.sml*rsq)goto 10792                                    1092
      if(rsq.gt.rsqmax)goto 10792                                          1093
10791 continue                                                             1094
10792 continue                                                             1094
      deallocate(a,mm,g,ix)                                                1095
      return                                                               1096
      end                                                                  1097
