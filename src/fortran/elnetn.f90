subroutine elnetn (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni)
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)
      integer jd(*),ia(nx),nin(nlam)
      real, dimension (:), allocatable :: xm,xs,xv,vlam
      integer, dimension (:), allocatable :: ju
      allocate(xm(1:ni),stat=jerr)
      allocate(xs(1:ni),stat=ierr)

      jerr=jerr+ierr
      allocate(ju(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xv(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(vlam(1:nlam),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return

      call chkvars(no,ni,x,ju)

      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0

      if(maxval(ju) .gt. 0)goto 10581
            jerr=7777
            return
10581 continue

      call standard1(no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)

      if(jerr.ne.0) return
      cl=cl/ys

      if(isd .le. 0)goto 10601
      10610 do 10611 j=1,ni
                  cl(:,j)=cl(:,j)*xs(j)
            10611 continue
      10612 continue
10601 continue

      if(flmin.ge.1.0) vlam=ulam/ys

      call elnet2(parm,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)

      if(jerr.gt.0) return

10620 do 10621 k=1,lmu
            alm(k)=ys*alm(k)
            nk=nin(k)
      10630 do 10631 l=1,nk
                  ca(l,k)=ys*ca(l,k)/xs(ia(l))
            10631 continue
      10632 continue
            a0(k)=0.0
            if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))
10621 continue

10622 continue

      deallocate(xm,xs,ju,xv,vlam)
      return
      end
