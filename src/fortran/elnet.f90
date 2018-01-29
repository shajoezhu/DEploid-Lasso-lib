subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni)
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)
      integer jd(*),ia(nx),nin(nlam)
      real, dimension (:), allocatable :: vq;
      if(maxval(vp) .gt. 0.0)goto 10021
            jerr=10000
            return
10021 continue

      allocate(vq(1:ni),stat=jerr)
      if(jerr.ne.0) return

      vq=max(0.0,vp)
      vq=vq*ni/sum(vq)

      if(ka .ne. 1)goto 10041
            call elnetu  (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
            goto 10051
10041 continue
            call elnetn (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue

10031 continue
      deallocate(vq)
      return
      end
