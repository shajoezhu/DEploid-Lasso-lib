      subroutine standard1 (no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)    989
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        989
      integer ju(ni)                                                        990
      real, dimension (:), allocatable :: v
      allocate(v(1:no),stat=jerr)                                           993
      if(jerr.ne.0) return                                                  994
      w=w/sum(w)                                                            994
      v=sqrt(w)                                                             995
      if(intr .ne. 0)goto 10651                                             995
      ym=0.0                                                                995
      y=v*y                                                                 996
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                         996
      y=y/ys                                                                997
10660 do 10661 j=1,ni                                                       997
      if(ju(j).eq.0)goto 10661                                              997
      xm(j)=0.0                                                             997
      x(:,j)=v*x(:,j)                                                       998
      xv(j)=dot_product(x(:,j),x(:,j))                                      999
      if(isd .eq. 0)goto 10681                                              999
      xbq=dot_product(v,x(:,j))**2                                          999
      vc=xv(j)-xbq                                                         1000
      xs(j)=sqrt(vc)                                                       1000
      x(:,j)=x(:,j)/xs(j)                                                  1000
      xv(j)=1.0+xbq/vc                                                     1001
      goto 10691                                                           1002
10681 continue                                                             1002
      xs(j)=1.0                                                            1002
10691 continue                                                             1003
10671 continue                                                             1003
10661 continue                                                             1004
10662 continue                                                             1004
      go to 10700                                                          1005
10651 continue                                                             1006
10710 do 10711 j=1,ni                                                      1006
      if(ju(j).eq.0)goto 10711                                             1007
      xm(j)=dot_product(w,x(:,j))                                          1007
      x(:,j)=v*(x(:,j)-xm(j))                                              1008
      xv(j)=dot_product(x(:,j),x(:,j))                                     1008
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1009
10711 continue                                                             1010
10712 continue                                                             1010
      if(isd .ne. 0)goto 10731                                             1010
      xs=1.0                                                               1010
      goto 10741                                                           1011
10731 continue                                                             1011
10750 do 10751 j=1,ni                                                      1011
      if(ju(j).eq.0)goto 10751                                             1011
      x(:,j)=x(:,j)/xs(j)                                                  1011
10751 continue                                                             1012
10752 continue                                                             1012
      xv=1.0                                                               1013
10741 continue                                                             1014
10721 continue                                                             1014
      ym=dot_product(w,y)                                                  1014
      y=v*(y-ym)                                                           1014
      ys=sqrt(dot_product(y,y))                                            1014
      y=y/ys                                                               1015
10700 continue                                                             1015
      deallocate(v)                                                        1016
      return                                                               1017
      end                                                                  1018
