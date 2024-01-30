      integer function scan(char,charset)
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      character * (*) char,charset
      ilen1 =len(char)
      ilen2 = len(charset)
      do i = 1,ilen1
        if (index(charset,char(i:i)) .ne. 0) then
          scan = i
          return
        endif
      enddo
      scan = 0
      return
      END
