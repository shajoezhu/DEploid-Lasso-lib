      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni)                  959
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         960
      integer jd(*),ia(nx),nin(nlam)                                        961
      real, dimension (:), allocatable :: xm,xs,xv,vlam
      integer, dimension (:), allocatable :: ju
      allocate(xm(1:ni),stat=jerr)                                          966
      allocate(xs(1:ni),stat=ierr)                                          966
      jerr=jerr+ierr                                                        967
      allocate(ju(1:ni),stat=ierr)                                          967
      jerr=jerr+ierr                                                        968
      allocate(xv(1:ni),stat=ierr)                                          968
      jerr=jerr+ierr                                                        969
      allocate(vlam(1:nlam),stat=ierr)                                      969
      jerr=jerr+ierr                                                        970
      if(jerr.ne.0) return                                                  971
      call chkvars(no,ni,x,ju)                                              972
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  973
      if(maxval(ju) .gt. 0)goto 10581                                       973
      jerr=7777                                                             973
      return                                                                973
10581 continue                                                              974
      call standard1(no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)           975
      if(jerr.ne.0) return                                                  976
      cl=cl/ys                                                              976
      if(isd .le. 0)goto 10601                                              976
10610 do 10611 j=1,ni                                                       976
      cl(:,j)=cl(:,j)*xs(j)                                                 976
10611 continue                                                              976
10612 continue                                                              976
10601 continue                                                              977
      if(flmin.ge.1.0) vlam=ulam/ys                                         978
      call elnet2(parm,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  981
10620 do 10621 k=1,lmu                                                      981
      alm(k)=ys*alm(k)                                                      981
      nk=nin(k)                                                             982
10630 do 10631 l=1,nk                                                       982
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          982
10631 continue                                                              982
10632 continue                                                              982
      a0(k)=0.0                                                             983
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))           984
10621 continue                                                              985
10622 continue                                                              985
      deallocate(xm,xs,ju,xv,vlam)                                          986
      return                                                                987
      end
