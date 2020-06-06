! Программа моделирования негидростатики 3D
! Разделение переменных на 2 части: переменные мелкой воды и дельта-переменные

program main
  !use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
  use variables
  implicit none

  !print '(3a)', 'This file was compiled by ', compiler_version(), ' using the options:'
  !print '(2a)', '  ', compiler_options()

  open(17, file="debout.csv")

  call Init
  write(*,*) "Init done."

  !call DebTest

  !call PrintInitData
  !       call TPoutPltX
  !       call TPoutPltY

  call Gout

  !write(*,*) "PrintInitData"

  call Steps

end program main

!-----------------------------------------------------------------------------------------------

subroutine info(startTime)
  use variables
  implicit none
  real(4) :: startTime
  integer(4) :: curTime, difTime, estimateTime, rest
  integer(4) :: td, th, tm, ts, tde, the, tme, tse
  integer(4) :: ix, iy, iz, isl, isr, isb, ist
  real(R8) :: vel, velxy, maxVel3, maxVel2, Et2, Et3, Ek2, Ek3, Ep, Ek2i, Ek3i, mci, wci, mtot, z1, z2
  real(R8) :: uci, vci, ucni, vcni, voli, zci, zcni, zl, zr, zb, zt, hmin, hmax, hmin0, hmax0, tetaMin, tetaMax
  real(R8) :: fT, fB, hsl, hsr, hsb, hst, hi, tci, tmci, tTot, tmTot
  real(R8), save :: Ep0, vol0
  real(R8) :: umx(nx,ny), umy(nx,ny)                 ! средние по столбцу потоковые скорости
  real(R8) :: f(nx, ny, nz)
  real(8):: soundk, dtLim_x, dtLim_y
  real(8):: z                                         ! число Зайцев
  logical, save :: isFirstTimeHere = .true.
  character(len=100), save:: fname

  curTime = second()
  difTime = curTime - int(startTime)
  ts = mod(difTime,60)
  rest = (difTime - ts) / 60
  tm = mod(rest, 60)
  th = (rest - tm) / 60
  if(istep==0) then
    estimateTime = 0
  else
    if(nt>0) then
      estimateTime = (real(difTime) / real(istep)) * nt
    else
      estimateTime = (real(difTime) / real(ttime)) * maxt
    end if
  end if

  if(istep==0.) then
    z = 0
  else
    z = difTime / ((nx-1.d0)*(ny-1.d0)*(nz-1.d0)*istep) * 1.d6
  end if

  tse = mod(estimateTime,60)                    ! прогноз секунды
  rest = (estimateTime - tse) / 60              ! остаток в минутах
  tme = mod(rest, 60)                           ! прогноз минуты
  the = (rest - tme) / 60                       ! прогноз часы
  if(th<100 .and. the<100) then
    !                                       дни         dt step,iout    зайцы    cf
    !                                        V          V      V          V      V
    write(*,"(6(i2.2,a),i0, a,1p,g12.4, a,0p,f0.2, a,1p,g12.4, 2(a,i0), a,0p,f0.2,a,  a,1p,g12.4)") &
      th,':',tm,':',ts,'/',the,':',tme,':',tse,'  ', &
      printproc, &
      '%. t=', GetCTime(ttime), ' (', GetCTime(ttime)/(3600.*24.), ' дн) dt=', dt, &
      ' step: ', istep, ' iout=', iout, ' (', z, ' Зайцев)'
  else
    rest = th                                   ! прошло часы
    th = mod(rest, 24)                          ! часы (0..24)
    td = (rest - th) / 24                       ! дни
    rest = the                                  ! остаток часы
    the = mod(rest, 24)                         ! прогноз часы (0..24)
    tde = (rest - tme) / 24                     ! прогноз дни
    write(*,"(2(i0,a,3(i2.2,a)),i0, a,1p,g12.4, a,0p,f0.2, a,1p,g12.4, 2(a,i0), a,0p,f0.2,a, a,1p,g12.4)") &
      td,':',th,':',tm,':',ts,'/',tde,':',the,':',tme,':',tse,'  ', &
      printproc, &
      '%. t=', GetCTime(ttime), ' (', GetCTime(ttime)/(3600.*24.), ' дн) dt=', dt, &
      ' step: ', istep, ' iout=', iout, ' (', z, ' зайцев)'
  end if

  ! статистика:
  vol = 0.
  mtot = 0.
  tTot = 0
  tmTot = 0
  maxVel3 = 0.
  maxVel2 = 0.
  Ek2 = 0.
  Ek3 = 0.
  Ep = 0.

  hmax = -1.d50; hmax0 = hmax; tetaMax = hmax
  hmin =  1.d50; hmin0 = hmin; tetaMin = hmin

  do ix=1, ncx
    isl = ix                                                            ! индекс левой грани
    isr = ix + 1                                                        ! индекс правой грани
    do iy=1, ncy
      if(c_type(ix,iy)<=CELL_DELETED) cycle
      isb = iy
      ist = iy + 1
      do iz=1, ncz

        wci = c_w(ix, iy, iz)                               ! вертикальная скорость
        hi = fz_z(ix, iy, iz) - fz_z(ix, iy, iz+1)          ! толщина слоя в ячейке

        voli = c_dx(ix) * c_dy(iy) * hi                     ! объем слой-ячейки
        mci = c_rho(ix, iy, iz) * voli                      ! масса слой-ячейки
        tci = c_teta(ix, iy, iz) * voli                     ! псевдо-масса слой-ячейки
        tmci = c_teta(ix, iy, iz) * c_rho(ix, iy, iz) * voli
        zcni = (fz_z(ix, iy, iz) + fz_z(ix, iy, iz+1)) / 2.

        ! кинетическая энергия, с учетом вертикальной скорости
        Ek3i = sqrt(c_u(ix, iy, iz)**2 + c_v(ix, iy, iz)**2 + wci**2) * mci
        Ek3 = Ek3 + Ek3i

        ! кинетическая энергия, без учета вертикальной скорости
        Ek2i = sqrt(c_u(ix, iy, iz)**2 + c_v(ix, iy, iz)**2) * mci
        Ek2 = Ek2 + Ek2i

        ! потенциальная энергия:
        Ep = Ep + g * mci * zcni

        ! объем и масса
        vol = vol + voli
        mTot = mTot + mci
        tTot = tTot + tci
        tmTot = tmTot + tmci

        velxy = sqrt(c_u(ix, iy, iz)**2 + c_v(ix, iy, iz)**2)
        if(velxy>maxVel2) maxVel2 = velxy
        vel = sqrt(c_u(ix, iy, iz)**2 + c_v(ix, iy, iz)**2 + c_w(ix,iy,iz)**2)
        if(vel>maxVel3) maxVel3 = velxy

        if(c_teta(ix,iy,iz)>tetaMax) tetaMax = c_teta(ix,iy,iz)
        if(c_teta(ix,iy,iz)<tetaMin) tetaMin = c_teta(ix,iy,iz)

      end do

      hi = fz_z(ix, iy, 1) - fz_z(ix, iy, nz)           ! реальная глубина
      if(hi>hmax) hmax = hi
      if(hi<hmin) hmin = hi

      hi = c_h0(ix, iy)                                 ! глубина мелкой воды
      if(hi>hmax0) hmax0 = hi
      if(hi<hmin0) hmin0 = hi

    end do
  end do

  if(isFirstTimeHere) then
    Ep0 = Ep
    vol0 = vol
  end if

  Et2 = (Ek2 + Ep - Ep0) / vol0
  Et3 = (Ek3 + Ep - Ep0) / vol0

  if(isFirstTimeHere) then
    write(fname,"(a,i0,a)") "./stat-", taskNum, ".csv"
    open(117, file=fname)
    write(117,*) "t(sec);t(day);wind;maxVel2;maxVel3;vol;M;t-M;tm-M;Et2;Et3;Ek2;Ek3;Ep;hmin;hmax;hmin_sw;hmax_sw;minTeta;maxTeta"
    isFirstTimeHere=.false.
  else
    open(117, file=fname, access="append")
  end if

  write(117,"(50(1p,g20.12,';'))") ttime, ttime/(3600.*24.), cfSurfFrictionU, maxvel2, maxvel3, &
    vol, mTot, tTot, tmTot, Et2, Et3, Ek2, Ek3, Ep - Ep0, hmin, hmax, hmin0, hmax0, tetaMin, tetaMax
  close(117)

end subroutine info

!-----------------------------------------------------------------------------------------------

subroutine Steps
  use variables
  real(4), save :: curTime, lastTime, startTime, esimateTime
  integer(4) :: i, j, k
  real(R8) :: dx, dy, dxMin, dyMin
  integer:: proc1 = 0, proc2 = 0, proc
  integer(4), parameter :: nts = 100
  logical:: needPltOut

  call SaveTaskParams

  !call CheckConstX

  call BuildFluxesSw                                            ! начальные потоки для МВ
#if 0
  call BuildFluxes
#endif

  do while(.true.)
    call CalcDt                                                 ! вычисление шага по времени

    istep=istep+1                                               ! Счетчик временных шагов
    ttime=ttime+dt

    if(IsDeb(-1)) then
      write(17,"(/a,i0,20(a,1p,g20.12))") &
        "------------;---;---;---;step:;",istep,";time:;",ttime,";dt:;",dt
      flush(17)
    end if

    if (istep .eq. 1) then
      startTime = second()
      call info(startTime)
      lastTime = startTime
    endif

    !if((NSTEP.ge.NT).or.(NSTEP.eq.1)) then
    !  minHeight = 10e10
    !  do J=1,NY-1
    !    do I=1,NX-1
    !      if (c_z(I,J,1)-c_z(I,J,NZ) .lt. minHeight) then
    !        minHeight =c_z(I,J,1)-c_z(I,J,NZ)
    !      endif
    !    enddo
    !  enddo
    !  write(*,*) "NSTEP=",NSTEP," minHeight=",minHeight
    !endif

    !if(ttime>4.e5) then
    !  cfSurfFrictionU=0.
    !  cfViscVert = 0
    !end if

    call Step

    !if(MOD(NSTEP,NPRINT*20) .eq. 0) then
    ! call SaveAllData
    !endif

    ! разбираемся, нужен ли вывод на этом шаге:
    needPltOut = .false.
    if(nprint>0) then
      if((istep-1)/nprint /= istep/nprint) needPltOut = .true.
    end if
    if(tprint>0) then
      if(floor((ttime-dt)/tprint) /= floor(ttime/tprint)) needPltOut = .true.
    end if

    if(needPltOut) then
      !call PrintData
      !call TPoutPltX
      !call TPoutPltY

      call Gout

    endif

    if(nt>0) proc1 = istep * 100 / NT
    if(maxt>0.) proc2 = ttime * 100. / maxt
    proc = max(proc1, proc2)

    if (proc>printproc) then
      call info(starttime)
      printproc = proc
      lasttime = second()
    else
      curTime = second()
      if(curTime-lastTime>20) then
        call info(startTime)
        lastTime = curTime
      !else if(mod(nstep,25)==0) then
      !  call info(startTime)
      !  lastTime = curTime
      end if
    endif

    if((nt>0 .and. istep.ge.nt) .or. (maxt>0 .and. ttime>maxt)) then
#if DABC>0
      close(17,status='keep')
#endif
      close(27,status='keep')
      close(21,status='keep')
      close(22,status='keep')
      close(23,status='keep')
      close(24,status='keep')
      close(25,status='keep')
      close(26,status='keep')
      exit
    endif

  !      if(nts>0) then
  !        if(mod(nstep,nts)==0) call TestSym
  !      end if

  end do

  write(*, *) "-- Done --------------------------------------------"

end subroutine Steps

!-----------------------------------------------------------------------------------------------

subroutine Step
  use variables
  implicit none
  integer :: i,j,k
  logical isFirstDeb

  call SetupVariableBCs                                         ! обработка данных для ГУ, переменных во времени

  call ShWater                                                  ! расчет 0-переменных по мелкой воде

  call Phase1
  call Phase2X
  call Phase2Y
  call Phase2Z_in
  call Phase2Z_top
  call Phase3

  !call DebTest

  !call CheckConstX
  !call TestSym                      ! тест на симметрию относительно Y=0: h(y)=h(-y), u(y)=u(-y), v(y)=-v(-y)

# if 0
  do i=1,ncx; do j=1,ncy; do k=1,ncz
    if(c_dteta(i,j,k)/=cs_dteta(k,j,i)) then
      debdum=0
    endif
    if(c_drho(i,j,k)/=cs_drho(i,j,k)) then
      debdum=1
    endif
    if(c_du(i,j,k)/=cs_du(i,j,k)) then
      debdum=2
    endif
    if(c_dv(i,j,k)/=cs_dv(i,j,k)) then
      debdum=3
    endif
    if(c_dw(i,j,k)/=cs_dw(i,j,k)) then
      debdum=4
    endif
  end do; end do; end do
#endif

#if 0
  do i=1,nxx; do j=1,nxy; do k=1,ncz
    if(fx_dteta(i,j,k)/=fxn_dteta(i,j,k)) then
      debdum=0
    endif
    if(fx_drho(i,j,k)/=fxn_drho(i,j,k)) then
      debdum=1
    endif
    if(fx_du(i,j,k)/=fxn_du(i,j,k)) then
      debdum=2
    endif
    if(fx_dv(i,j,k)/=fxn_dv(i,j,k)) then
      debdum=3
    endif
    if(fx_dw(i,j,k)/=fxn_dw(i,j,k)) then
      debdum=4
    endif
  end do; end do; end do
#endif

  call NextStep

end subroutine Step

!-----------------------------------------------------------------------------------------------------

! Вычисление величины шага по времени
subroutine CalcDt
  use variables
  integer(4) :: i, j, k
  real(R8):: dx, dy, dtMin, sound, uc, vc, wc, dtk

  dtMin = 1.d100
  maxSound = 0.

  do i=1,ncx
    dx = c_dx(i)
    if(i==31) &
      debdum = 1

    do j=1,ncy
      if(c_type(i,j)<=CELL_DELETED) cycle                ! пропускаем удалённые ячейки

      dy = c_dy(j)

      ! ограничение на шаг мелкой воды:
      sound = sqrt(g * c_h0(i,j))
      dtk = dx / (abs(c_u0(i,j)) + sound)
      if(dtk<dtmin) then
        dtmin = dtk
      end if
      dtk = dy / (abs(c_v0(i,j)) + sound)
      if(dtk<dtmin) then
        dtmin = dtk
      end if

      ! ограничения от негидростатики:
      do k=1,ncz                                          ! цикл по слоям

        ! по направлению X:
        uc = c_u0(i,j) + c_du(i,j,k)                      ! реальная X-скорость в слой-ячейке
        dtk = dx / (abs(uc) + cs_sound)
        if(dtk<dtmin) then
          dtmin = dtk
        endif

        ! по направлению Y:
        vc = c_v0(i,j) + c_dv(i,j,k)                      ! реальная Y-скорость в слой-ячейке
        dtk = dy / (abs(vc) + cs_sound)
        if(dtk<dtmin) then
          dtmin = dtk
        endif

        ! по направлению Z:
        wc = c_GetW0(i,j,k) + c_dw(i,j,k)                  ! реальная W-скорость в слой-ячейке
        dtk = (fz_z(i,j,k) - fz_z(i,j,k+1)) / (abs(wc) + cs_sound)
        if(dtk<dtmin) then
          dtmin = dtk
        endif

      enddo

    enddo
  enddo

  !if(dtmin>dtVisc) dtmin = dtVisc                        ! ограничение от вязкости

  dt = CFL * dtmin

end subroutine CalcDt

!-----------------------------------------------------------------------------------------------

! для каждой грани переписываем все параметры ufn -> uf
subroutine NextStep
  use variables
  implicit none

  ! Переход к новому слою на гранях. Вместо копирования массивов меняем местами указатели:
  call ExchangePointers

  c_zt(:,:)=fz_z(:,:,1)

  ! мелкая вода:
  c_z0 = cn_z0
  c_h0 = cn_h0
  c_u0 = cn_u0
  c_v0 = cn_v0
  c_rho0 = cn_rho0

end subroutine NextStep

!-----------------------------------------------------------------------------------------------------------------

! сохраняет параметры задачи в файл в директорию out для идентификации расчета
subroutine SaveTaskParams
  use variables
  implicit none

  return

  open(27,file="./out/info.txt")

  write(27,*) "taskNum=",taskNum
  write(27,*) "nx=", nx
  write(27,*) "ny=", ny
  write(27,*) "nz=", nz
  write(27,*) "dl=", dl
  write(27,*) "dw=", dw
  write(27,*) "bcTypeX=", bcTypeX
  write(27,*) "bcTypeY=", bcTypeY
  write(27,*) "bcTypeTopT=", bcTypeTopT
  write(27,*) "bcTypeTopU=", bcTypeTopU
  write(27,*) "CFL=", CFL
  write(27,*) "PAN=", PAN
  write(27,*) "Coriolis=", Coriolis
  write(27,*) "cfViscVert=", cfViscVert
  write(27,*) "cfViscVert=", cfViscVert
  write(27,*) "cfViscVertT=", cfViscVertT
  write(27,*) "cf3ViscVertT=", cf3ViscVertT
  write(27,*) "cfBottomFriction=", cfBottomFriction
  write(27,*) "cfSurfFrictionU=", cfSurfFrictionU
  write(27,*) "cfSurfFrictionV=", cfSurfFrictionV
  write(27,*) "TSurf=", TSurf(1,1) ! возможно удалить!!!
  write(27,*) "uAtm=", uAtm
  write(27,*) "SIG_0=", SIG_0
  write(27,*) "sig_ex=", sig_ex
  write(27,*) "CFD_L=", CFD_L
  write(27,*) "CFD_H=", CFD_H
  write(27,*) "CFD_T=", CFD_T
  write(27,*) "CFD_TEMP=", CFD_TEMP
  write(27,*) "#-----------------------------------------------------------------------"
  write(27,*) "isFixCover=", isFixCover
  write(27,*) "isExchangePh1=", isExchangePh1
  write(27,*) "isImplExchange=", isImplExchange
  write(27,*) "isImplExchangeS2=", isImplExchangeS2
  write(27,*) "isImplExchangeS3=", isImplExchangeS3
  write(27,*) "isImplCoriolis=", isImplCoriolis
  write(27,*) "isConsFilterZ0=", isConsFilterZ0
  write(27,*) "#-----------------------------------------------------------------------"
  write(27,*) "# Options:"
  write(27,*) "optLimitMinRhoC=", optLimitMinRhoC
  write(27,*) "optLimitMinHC=", optLimitMinHC
  write(27,*) "optLimitMinRhoF=", optLimitMinRhoF
  write(27,*) "optLimitMinHF=", optLimitMinHF
  write(27,*) "optNoAssert=", optNoAssert

  close(27)
end subroutine
