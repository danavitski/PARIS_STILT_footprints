      subroutine prmprt(stcprm)
! CHG(03/12/03)
      IMPLICIT REAL*8 (A-H,O-Z)

      real*8 stcprm(15)
      ifind(l,k) = -2 + 3*l +k
      write(*,*)' Gamma = ',stcprm(1)
      do k=1,3
        write(*,'('' line '',i1,'':  '',$)')k-1
        write(*,'(3f10.6)') (stcprm(ifind(l,k)),l=1,3)
      enddo
      write(*,*)'x0 = ',stcprm(11),', y0 = ',stcprm(12)
      write(*,*)'cos(th) = ',stcprm(13),', sin(th) = ',stcprm(14)
      write(*,*)'gridsize = ',stcprm(15)
      write(*,*)
      return
      END

      integer function ans_pt(prompt,values)
      integer scan
      character * (*) prompt,values
      character*80 instr
        write(*,'(1x,a,$)') prompt
        read(*,'(a)') instr
        ans_pt = (scan(values,instr)+1)/2
      return
      END

      PROGRAM dmapftst
      PARAMETER (REARTH=6367.47)
      real*8 stcprm(15)
      integer ano_pt,ans_pt,choice,ano_pj
      alat = 0.
      ano_pj = 1
      do  while (ano_pj .gt. 0)
        write(*,*) 'O - Oblique Stereographic'
        write(*,*) 'L - Lambert Polar'
        write(*,*) 'T - Transverse Mercator'
        write(*,*) 'U - Oblique Mercator'
        write(*,*) 'C - Oblique Lambert'
        choice = ans_pt('Enter Choice: ','OoLlTtUuCc')
        if (choice .eq. 1) then
          write(*,'('' Enter Latitude and Longitude: '',$)')
          read(*,*)alat,along
          call sobstr(stcprm,alat,along)
        else if (choice .eq. 2) then
          write(*,'('' Enter Two Ref. Latitudes and One Longitude: '',  &
     &            $)')
          read(*,*)alat1,alat2,along
          write(*,*)'Eqvlat = ',eqvlat(alat1,alat2)
          call stlmbr(stcprm,eqvlat(alat1,alat2),along)
        else if (choice .eq. 3) then
          write(*,'('' Enter Latitude and Longitude: '',$)')
          read(*,*)alat,along
          call stvmrc(stcprm,alat,along)
        else if (choice .eq. 4) then
          write(*,'(a,$)') 'Enter Central latitude and longitude: '
          read(*,*) alat1,along1
          write(*,'(a,$)') 'Enter Secondary latitude and longitude: '
          read(*,*) alat2,along2
          call sobmrc(stcprm, alat1, along1, alat2, along2)
        else if (choice .eq. 5) then
          write(*,'(a,$)') 'Enter Central latitude and longitude: '
          read(*,*) alat,along
          write(*,'(a,$)') 'Enter lat & lon of second point on circle: '
          read(*,*) alat1,along1
          write(*,'(a,$)') 'Enter lat & lon of third point on circle: '
          read(*,*) alat2,along2
          call soblmb(stcprm, alat1,along1, alat,along, alat2,along2)
        endif
        choice = ans_pt('1-point or 2-point scaling? ','1o2t')
        if (choice .eq. 1) then
          write(*,'(a,$)') 'Enter x,y of anchor point: '
          read(*,*) x,y
          write(*,'(a,$)') 'Enter lat,long of anchor point: '
          read(*,*) alat,along
          write(*,'(a,$)') 'Enter lat,long of reference point: '
          read (*,*) rlat,rlong
          write(*,'(a,$)') 'Enter grid size at reference point: '
          read(*,*) grdsiz
          write(*,'(a,$)')                                              &
     &       'Enter y-axis orientation at reference point: '
          read(*,*) orient
          call stcm1p(stcprm, x,y, alat,along,                          &
     &        rlat,rlong, grdsiz,orient)
        else
          write(*,'(a,$)') 'Enter x,y of first anchor point: '
          read(*,*) x1,y1
          write(*,'(a,$)') 'Enter lat,long of first anchor point: '
          read(*,*) alat1,along1
          write(*,'(a,$)') 'Enter x,y of second anchor point: '
          read(*,*) x2,y2
          write(*,'(a,$)') 'Enter lat,long of second anchor point: '
          read(*,*) alat2,along2
          call stcm2p(stcprm,x1,y1,alat1,along1, x2,y2,alat2,along2)
        endif
        call prmprt(stcprm)
        ano_pt = 1
        do while (ano_pt .ne. 0)
          ano_pt = ans_pt('Translate x,y point? (y/n) ','yY')
          if (ano_pt .ne. 0) then
            write(*,'(1x,a,$)') 'Enter x,y: '
            read(*,*) x,y
            call cxy2ll(stcprm, x,y, alat,along)
            call cll2xy(stcprm, alat,along, oldx,oldy)
            write(*,*) 'x,y = (',oldx,',',oldy,'), lat,long = (',alat,  &
     &         ',',along,').'
            write(*,*) 'gridsize(x,y) = ', cgszxy(stcprm,x,y),          &
     &         'gridsize(l,l) = ', cgszll(stcprm,alat,along)
            call cpolxy(stcprm, x,y, enx,eny,enz)
            write(*,*) 'Polar axis from x,y = (',                       &
     &         enx,',',eny,',',enz,')'
            call cpolll(stcprm, alat,along, enx,eny,enz)
            write(*,*) 'Polar axis from lat,long = (',                  &
     &         enx,',',eny,',',enz,')'
            call cgrnxy(stcprm, x,y, enx,eny,enz)
            write(*,*) 'Greenwich axis from x,y = (',                   &
     &         enx,',',eny,',',enz,')'
            call cgrnll(stcprm, alat,along, enx,eny,enz)
            write(*,*) 'Greenwich axis from lat,long = (',              &
     &            enx,',',eny,',',enz,')'
            call cc2gxy(stcprm, x,y, 0.,10., ug,vg)
            call cg2cxy(stcprm, x,y, ug,vg, ue,vn)
            write(*,*) 'x,y winds from (E,N) to (Ug,Vg):(',ue,',',vn,   &
     &         ') to (',ug,',',vg,')'
            call cc2gll(stcprm, alat,along, 0.,10., ug,vg)
            call cg2cll(stcprm, alat,along, ug,vg, ue,vn)
            write(*,*) 'l,l winds from (E,N) to (Ug,Vg):(',ue,',',vn,   &
     &            ') to (',ug,',',vg,')'
            call ccrvxy(stcprm, x,y, gx,gy)
            write(*,*) 'x,y curvature vector (gx,gy):(',gx,',',gy,')'
            call ccrvll(stcprm, alat,along, gx,gy)
            write(*,*) 'lat,long curvature vector (gx,gy):('            &
     &         ,gx,',',gy,')'
           endif
         enddo
      ano_pj = ans_pt('Another Projection? (y/n) ','yY')
      enddo
      stop
      END PROGRAM dmapftst
