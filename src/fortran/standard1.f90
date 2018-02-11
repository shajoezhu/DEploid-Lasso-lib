subroutine standard1 (no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)
      integer ju(ni)
      real, dimension (:), allocatable :: v
      allocate(v(1:no),stat=jerr)
      if(jerr.ne.0) return
      w=w/sum(w)
      v=sqrt(w)

      if(intr .ne. 0)goto 10651
            ym=0.0
            y=v*y
            ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)
            y=y/ys

      10660 do 10661 j=1,ni
                  if(ju(j).eq.0)goto 10661
                        xm(j)=0.0
                        x(:,j)=v*x(:,j)
                        xv(j)=dot_product(x(:,j),x(:,j))

                        if(isd .eq. 0)goto 10681
                              xbq=dot_product(v,x(:,j))**2
                              vc=xv(j)-xbq
                              xs(j)=sqrt(vc)
                              x(:,j)=x(:,j)/xs(j)
                              xv(j)=1.0+xbq/vc
                              goto 10691
                  10681 continue
                              xs(j)=1.0
                  10691 continue
                  10671 continue
            10661 continue
      10662 continue
            go to 10700
10651 continue

10710 do 10711 j=1,ni
            if(ju(j).eq.0)goto 10711
            xm(j)=dot_product(w,x(:,j))
            x(:,j)=v*(x(:,j)-xm(j))
            xv(j)=dot_product(x(:,j),x(:,j))
            if(isd.gt.0) xs(j)=sqrt(xv(j))
      10711 continue
      10712 continue

      if(isd .ne. 0)goto 10731
            xs=1.0
            goto 10741
      10731 continue

      10750 do 10751 j=1,ni
            if(ju(j).eq.0)goto 10751
            x(:,j)=x(:,j)/xs(j)
      10751 continue
      10752 continue
            xv=1.0
      10741 continue
10721 continue

      ym=dot_product(w,y)
      y=v*(y-ym)
      ys=sqrt(dot_product(y,y))
      y=y/ys
10700 continue

      deallocate(v)
      return
      end
