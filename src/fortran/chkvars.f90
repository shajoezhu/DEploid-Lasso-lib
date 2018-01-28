subroutine chkvars(no,ni,x,ju)
      real x(no,ni)
      integer ju(ni)
11060 do 11061 j=1,ni
            ju(j)=0
            t=x(1,j)
      11070 do 11071 i=2,no
                  if(x(i,j).eq.t)goto 11071
                  ju(j)=1
                  goto 11072
            11071 continue
      11072 continue
      11061 continue
11062 continue
      return
end
