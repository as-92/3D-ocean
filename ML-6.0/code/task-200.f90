! начальные данные для задачи динамики Белого моря
subroutine SetupTask_200
  use variables
  implicit none

  interface
    subroutine GeoCoords2IJ(lat, lon, i, j)
      real(R8):: lat, lon
      integer:: i, j
    end subroutine
  end interface

  integer :: i, j, k, ic, jc
  integer:: bcTypeDef
  logical:: delL, delR
  real(R8):: gymid
  integer:: nxWS, nyWS, tpL, tpR
  real(R8):: r0, ri, z0, tMaxD, h

  integer:: file_WSdepth, file_WSmap                            ! файл с глубинами и карта БМ
  integer(1), pointer:: map(:,:)                                ! карта
  real(4), pointer:: depth(:,:)                                 ! глубины
  real(R8):: modu(10)
  real(R8):: nbx(10), nby(10)                                   ! координаты нормали к границе втока
  real(R8):: cspx(10), cspy(10)                                 ! средняя точка участка втока
  real(R8):: lb(10)                                             ! длина участка втока
  real(R8):: sb(10)                                             ! поперечное сечение втока
  integer:: n(10)                                               ! счётчики граней на границе с данным номером (для осреднения)
  integer:: ibc, cnt
  real(R8):: sideLen, l
  type(TBcData) datai                                           ! данные на границе
  integer:: file_topex_u, file_topex_v, file_topex_z            ! файлы с данными прилива
  integer:: nTopex, i0, j0, it
  real(R8):: d, dist, lat, lon, ux, uy, xc, yc, xt, yt

  ncx=133; ncy=110; ncz=2;                                      ! количество ячеек сетки

  ! время счёта:
#if 0
  nt=20000                                                      ! полное число шагов по времени
  nprint=max(1,nt/200)                                          ! интервал печати
#else
  tMaxD = 10.                                                   ! длительность расчёта (дни)
  tMaxD = 1. / 24. * 2.                                         ! длительность расчёта (дни)
  maxt = tMaxD * 60. * 60. * 24.                                ! длительность расчёта (сек)
  tprint=maxt/400.
#endif

  CFL = 0.3 ! Число Куранта
  pan=0.0 !- параметр Паниковского

  ! константы:
  g = 9.80616d0
  rho0 = 1.
  teta0 = 1.
  sound0 = 100. !0.15 * sqrt(g * dh)             ! псевдо-скорость звука

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  open(newunit=file_WSdepth, file='data/WS/bottom_depth.dat',  action='READ', err=100)
  open(newunit=file_WSmap,   file='data/WS/white_sea_map-abc.dat', action='READ', err=101)

  ! размеры рабочей области прочитаем из файла:
  read(file_WSdepth, *) gxmin, gxmax, gymax, gymin              ! в градусах (широта - от полюса)
  gymin = 90. - gymin; gymax = 90. - gymax                      ! в градусах (широта - от экватора)
  ! еще в файле есть глубины, но они будут прочитаны позже

  ! середина БМ:
  r0 = 6.37d6                                                   ! радиус Земли
  gymid = (gymin + gymax) / 2.                                  ! широта середины рабочей области (град)
  ri = r0 * cos(gymid / 180.d0 * PI)                            ! радиус параллели, проходящей по середине БМ
  lgradx = 2. * PI * ri / 360.                                  ! длина одного градуса параллели БМ
  lgrady = 2. * PI * r0 / 360.                                  ! длина одного градуса мередиана
  dl = lgradx * (gxmax - gxmin)
  dw = lgrady * (gymax - gymin)

  nx = ncx+1; ny = ncy+1; nz = ncz+1
  call allocate

  ! читаем и запоминаем карту и глубины:
  nxWS = 133; nyWS = 110                                        ! фиксированный размер карты в файле
  allocate (map(nxWS,nyWS), depth(nxWS,nyWS))
  do j=1, nyWS
    read(file_WSmap, '(10000(i1))') (map(i,j),i=1,nxWS)
  end do
  read(file_WSdepth,*) depth

  close(file_WSdepth)
  close(file_WSmap)

  ! координаты узлов:
  xmin=-dl/2; xmax=-xmin; do i=1,nxx; x(i) = xmin + (i-1.) * dl / (nxx-1); end do
  ymin=-dw/2; ymax=-ymin; do j=1,nyy; y(j) = ymin + (j-1.) * dw / (nyy-1); end do

  ! строим сетку (ячейки):
  z0 = 0.                                                       ! поверхность
  do i=1,ncx
    do j=1,ncy
      select case(map(i,ny-j))
        case(0)                                                 ! вода
          c_type(i,j) = CELL_INNER
        case(1)                                                 ! берег
          c_type(i,j) = CELL_DELETED
        case(3)                                                 ! водяные граничные ячейки
          c_type(i,j) = CELL_DELETED_SPEC
        case default
          write(*,*) "ERROR in SetupTask_200. i,j,map(i,j)", i, j, map(i,j)
      end select
      h = depth(i,ny-j)                                         ! глубина в ячейке
      c_b(i,j) = z0 - h
      ! z-сетка в ячейке:
      do k=1,nz
        fz_z(i,j,k) = z0 - (k-1.) * h / (nz-1.)
      end do
    end do
  end do
  deallocate (map, depth)

  !-- грани - внутренние, удалённые, граничные -------------------------------------------------------------

  bcTypeDef = BC_STICK                                          ! тип границы вода-берег
  fx_bc = 0; fy_bc = 0

  nBcData = 2                                                   ! резервируем 2 места для ГУ рек

  ! X-грани:
  do i=1,nxx
    do j=1,nxy
      if(i==1) then                                             ! левый край области
        if(c_type(1,j)<=CELL_DELETED) then                      ! если смежная ячейка была удалена
          fx_type(i,j) = BC_DELETED                             ! .. грань тоже удаляем
        else
          fx_type(i,j) = bcTypeDef                              ! .. иначе это граница вода-берег
        end if
      else if(i==nxx) then                                      ! правый край области
        if(c_type(ncx,j)<=CELL_DELETED) then                    ! если смежная ячейка была удалена
          fx_type(i,j) = BC_DELETED                             ! .. грань тоже удаляем
        else
          fx_type(i,j) = bcTypeDef                              ! .. иначе это граница вода-берег
        end if
      else                                                      ! внутренняя грань, есть две смежные ячейки
        tpL = c_type(i-1,j)
        tpR = c_type(i  ,j)
        delL = (tpL <= CELL_DELETED)                            ! была ли удалена ячейка слева
        delR = (tpR <= CELL_DELETED)                            ! была ли удалена ячейка справа
        if(delL .and. delR) then                                ! если обе смежные ячейки были удалены
          fx_type(i,j) = BC_DELETED                             ! .. грань тоже удаляем
        else if(delL .or. delR) then                            ! если только одна смежная ячейка была удалена
          if((.not.delL .and. tpR==CELL_DELETED_SPEC) .or. &
             (.not.delR .and. tpL==CELL_DELETED_SPEC)) then     ! .. если удалённая ячейка - это спец-ячейка
            call SetupWSxBC_200(i,j)                            ! .. .. устанавливаем специальные ГУ
          else                                                  ! .. иначе это грань вода-берег
            fx_type(i,j) = bcTypeDef
          end if
        else                                                    ! грань вода-вода
          call assertn(c_type(i-1,j)>=0, "SetupTask_200", i, j, int(c_type(i-1,j)))
          call assertn(c_type(i  ,j)>=0, "SetupTask_200", i, j, int(c_type(i  ,j)))
          fx_type(i,j) = BC_INNER
        end if
      end if
    end do
  end do

  ! Y-грани:
  do i=1,nyx
    do j=1,nyy
      if(j==1) then                                             ! левый край области
        if(c_type(j,1)==CELL_DELETED) then                      ! если смежная ячейка была удалена
          fy_type(i,j) = BC_DELETED                             ! .. грань тоже удаляем
        else
          fy_type(i,j) = bcTypeDef                              ! .. иначе это граница вода-берег
        end if
      else if(j==nyy) then                                      ! правый край области
        if(c_type(i,ncy)==CELL_DELETED) then                    ! если смежная ячейка была удалена
          fy_type(i,j) = BC_DELETED                             ! .. грань тоже удаляем
        else
          fy_type(i,j) = bcTypeDef                              ! .. иначе это граница вода-берег
        end if
      else                                                      ! внутренняя грань, есть две смежные ячейки
        tpL = c_type(i,j-1)
        tpR = c_type(i,j  )
        delL = (tpL <= CELL_DELETED)                            ! была ли удалена ячейка слева
        delR = (tpR <= CELL_DELETED)                            ! была ли удалена ячейка справа
        if(delL .and. delR) then                                ! если обе смежные ячейки были удалены
          fy_type(i,j) = BC_DELETED                             ! .. грань тоже удаляем
        else if(delL .or. delR) then                            ! если только одна смежная ячейка была удалена
          if((.not.delL .and. tpR==CELL_DELETED_SPEC) .or. &
             (.not.delR .and. tpL==CELL_DELETED_SPEC)) then     ! .. если удалённая ячейка - это спец-ячейка
            call SetupWSyBC_200(i,j)                            ! .. .. устанавливаем специальные ГУ
          else                                                  ! .. иначе это грань вода-берег
            fy_type(i,j) = bcTypeDef
          end if
        else                                                    ! грань вода-вода
          call assertn(c_type(i,j-1)>=0, "SetupTask_200", i, j, int(c_type(i,j-1)))
          call assertn(c_type(i,j-1)>=0, "SetupTask_200", i, j, int(c_type(i,j  )))
          fy_type(i,j) = BC_INNER
        end if
      end if
    end do
  end do

  ! начальные параметры:
  c_rho = rho0
  c_teta = teta0
  c_u = 0.
  c_v = 0.
  c_w = 0.

  ! мелкая вода:
  c_h0(:,:) = fz_z(:,:,1) - fz_z(:,:,nz)
  c_rho0 = rho0
  c_u0 = 0.
  c_v0 = 0.

  !-- ГУ стока рек ----------------------------------------------------------------------

  allocate(bcData(nBcData))

  ! для участков с ГУ вход-выход надо найти направление втока и поперечное сечение
  nbx=0; nby=0; cspx=0; cspy=0; lb=0; sb=0; n=0

  ! находим усреднённые нормали к границам втока:
  !
  ! пока в fx_bc и fy_bc записаны индексы водных границ (1, 2 и 3). Позже вместо индекса 3 (Горло)
  ! будут записаны индивидуальные индексы в bcData с индивидуальными, переменными во времени граничными данными

  ! цикл по X-граням:
  do i=1,nxx; do j=1,nxy
    ibc = fx_bc(i,j)
    if(ibc==0) cycle                                            ! пропускаем грань, если на ней нет спец. условия

    sideLen = (y(j+1) - y(j))                                   ! длина грани
    lb(ibc) = lb(ibc) + sideLen
    ! нормаль должна смотреть в сторону рабочей ячейки:
    if(fx_IsOwnRight(i,j)) then
      nbx(ibc) = nbx(ibc) + sideLen
    else
      nbx(ibc) = nbx(ibc) - sideLen
    end if
    cspx(ibc) = cspx(ibc) + x(i); cspy(ibc) = cspy(ibc) + (y(j)+y(j+1))/2.
    n(ibc) = n(ibc) + 1
  end do; end do

  ! цикл по Y-граням:
  do i=1,nyx; do j=1,nyy
    ibc = fy_bc(i,j)
    if(ibc==0) cycle                                            ! пропускаем грань, если на ней нет спец. условия

    sideLen = (x(i+1) - x(i))                                   ! длина грани
    lb(ibc) = lb(ibc) + sideLen
    ! нормаль должна смотреть в сторону рабочей ячейки:
    if(fy_IsOwnRight(i,j)) then
      nby(ibc) = nby(ibc) + sideLen
    else
      nby(ibc) = nby(ibc) - sideLen
    end if
    cspx(ibc) = cspx(ibc) + (x(i)+x(i+1))/2.; cspy(ibc) = cspy(ibc) + y(j)
    n(ibc) = n(ibc) + 1
  end do; end do

  ! осреднение нормали и координат центра:
  do ibc=1,10                                                   ! мы закладывали ограничение на 10 спец-границ
    if(n(ibc)==0) cycle                                         ! для этого типа границы ни одной грани не найдено
    ! вектор (nbx,nby) должен иметь единичную длину:
    l = sqrt(nbx(ibc)**2 + nby(ibc)**2)
    nbx(ibc) = nbx(ibc) / l; nby(ibc) = nby(ibc) / l
    cspx(ibc) = cspx(ibc) / real(n(ibc)); cspy(ibc) = cspy(ibc) / real(n(ibc))
  end do

  ! находим поперечное сечение S=h*l (по направлению усреднённой нормали):

  cnt = 2

  ! цикл по X-граням:
  do i=1,nxx; do j=1,nxy
    ibc = fx_bc(i,j)
    if(ibc==0) cycle                                            ! пропускаем грань, если на ней нет спец. условия

    ! с какой стороны от грани лежит "живая" ячейка?
    if(fx_IsOwnRight(i,j)) then
      ic = i
    else
      ic = i-1
    end if

    h = fz_z(ic,j,1) - fz_z(ic,j,nz)                            ! глубина в приграничной ячейке
    sideLen = (y(j+1) - y(j)) * abs(nbx(ibc))                   ! длина грани поперёк нормали
    sb(ibc) = sb(ibc) + h * sideLen

    if(ibc==3) then                                             ! для граней в Горле
      cnt = cnt + 1
      fx_bc(i,j) = cnt                                          ! .. записывем уникальный индекс в массиве bcData
    end if
  end do; end do

  ! цикл по Y-граням:
  do i=1,nyx; do j=1,nyy
    ibc = fy_bc(i,j)
    if(ibc==0) cycle                                            ! пропускаем грань, если на ней нет спец. условия

    ! с какой стороны от грани лежит "живая" ячейка?
    if(fy_IsOwnRight(i,j)) then
      jc = j
    else
      jc = j-1
    end if

    h = fz_z(i,jc,1) - fz_z(i,jc,nz)                            ! глубина в приграничной ячейке
    sideLen = (x(i+1) - x(i)) * abs(nby(ibc))                   ! длина грани поперёк нормали
    sb(ibc) = sb(ibc) + h * sideLen

    if(ibc==3) then                                             ! для граней в Горле
      cnt = cnt + 1
      fx_bc(i,j) = cnt                                          ! .. записывем уникальный индекс в массиве bcData
    end if
  end do; end do

  call assert(cnt==nBcData)

  ! сток рек (среднегодовой и месячный) [м^3/сек]:
  !                годовой янв   фев   мар   апр   май   июн   июл   авг   сен   окт   ноя   дек
  ! Северная Двина:  3331  1167  928   847  3646  14722 5623  3078  2238  2058  2432  2012  1320
  ! Онега:           505   178   145   133   480  1937   819   430   318   376   484   399   253

  ! скорости течения:
  modu = 0.
  modu(1) = 1973. / sb(1)     ! май
  modu(2) = 14722. / sb(2)

  ! для первых двух спец-ГУ (реки) строим вектор скорости течения (ux,uy) по нормали к границе:
  do ibc=1,2
    ux = modu(ibc) * nbx(ibc)                                   ! компоненты скорости
    uy = modu(ibc) * nby(ibc)

    datai.u = ux
    datai.v = uy
    datai.w = 0.
    datai.rho = rho0
    datai.teta = teta0

    bcData(ibc) = datai
  end do

  !-- ГУ прилива в горле ------------------------------------------------------------------------------------------

  open(newunit=file_topex_u, file="data/WS/gorlo_u.dat", action='read', err=102)
  open(newunit=file_topex_v, file="data/WS/gorlo_v.dat", action='read', err=102)
  open(newunit=file_topex_z, file="data/WS/gorlo_z.dat", action='read', err=102)

  ! пропускаем 3 строки заголовка:
  do i=1,3
    read(file_topex_z,"(A)")
    read(file_topex_u,"(A)")
    read(file_topex_v,"(A)")
  end do

  ! чтение topex-данных:
  nTopex = 0
  do while (.true.)
    nTopex = nTopex + 1
    read(file_topex_z,*,end=1) topex(nTopex).lat, topex(nTopex).lon, topex(nTopex).mAmp(1), topex(nTopex).mPh(1), topex(nTopex).sAmp(1), topex(nTopex).sPh(1)
    read(file_topex_u,*      ) topex(nTopex).lat, topex(nTopex).lon, topex(nTopex).mAmp(2), topex(nTopex).mPh(2), topex(nTopex).sAmp(2), topex(nTopex).sPh(2)
    read(file_topex_v,*      ) topex(nTopex).lat, topex(nTopex).lon, topex(nTopex).mAmp(3), topex(nTopex).mPh(3), topex(nTopex).sAmp(3), topex(nTopex).sPh(3)

    ! в какую ячейку попала точка topex (ближайшая живая ячейка):
    call GeoCoords2IJ(topex(nTopex).lat, topex(nTopex).lon, i0, j0)   ! ячейка (i0,j0) не обязательно существует

    ! ищем живую ячейку (ic,jc), ближайшую к (i0,j0)
    ic = -1; jc = -1
    do i=1, ncx; do j=1,ncy
      if(c_type(i,j)<=CELL_DELETED) cycle
      d = (i-i0)**2 + (j-j0)**2
      if(ic<0 .or. d<dist) then
        dist = d
        ic = i
        jc = j
      end if
    end do; end do

    topex(nTopex).i = ic
    topex(nTopex).j = jc

  end do
1 nTopex = nTopex - 1                                           ! последнее чтение было неудачным - откатываем назад

  close(file_topex_z)
  close(file_topex_u)
  close(file_topex_v)

  ! ищем топексы, относящиеся к гранями type=BC_IN_T

  do i=1,nxx; do j=1,nxy
    if(fx_type(i,j)/=BC_IN_T) cycle                             ! ниже обрабатываем только грани ГУ Горла
    xc = x(i); yc = (y(j) + y(j+1)) / 2.                        ! координата центра грани

    ! ищем ближайший топекс:
    it = -1
    do k=1, nTopex                                              ! цикл по топексам
      lat = topex(k).lat; lon = topex(k).lon                    ! точка топекса
      xt = xmin + (lon - gxmin) * lgradx
      yt = ymin + (lat - gymin) * lgrady

      d = sqrt((xt-xc)**2 + (yt-yc)**2)
      if(it<0 .or. d<dist) then
        dist = d
        it = k
      end if
    end do

    ! нашли topex(it) - ближайший к ячейке (ic,j)
    ibc = fx_bc(i,j)                                            ! данные ГУ для грани
    bcData(ibc).ind = it                                        ! для грани записали индекс соответствующего топекса
  end do; end do

  do i=1,nyx; do j=1,nyy
    if(fy_type(i,j)/=BC_IN_T) cycle                             ! ниже обрабатываем только грани ГУ Горла
    xc = (x(i) + x(i+1)) / 2.; yc = y(j)                        ! координата центра грани

    ! ищем ближайший топекс:
    it = -1
    do k=1, nTopex                                              ! цикл по топексам
      lat = topex(k).lat; lon = topex(k).lon                    ! точка топекса
      xt = xmin + (lon - gxmin) * lgradx
      yt = ymin + (lat - gymin) * lgrady

      d = sqrt((xt-xc)**2 + (yt-yc)**2)
      if(it<0 .or. d<dist) then
        dist = d
        it = k
      end if
    end do

    ! нашли topex(it) - ближайший к ячейке (ic,j)
    ibc = fy_bc(i,j)                                            ! данные ГУ для грани
    bcData(ibc).ind = it                                        ! для грани записали индекс соответствующего топекса
  end do; end do

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  !                                                                 !
  !                        О Т Л А Д К А                            !
  !                                                                 !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  debstartstep = 0; debendstep = 2
  debstartz = 28; debendz = 30

  !debcells(7, 102, :) = 1
  !debcells(24, 1, 29) = 1
  !debcells(28, 1, 24) = 1
  !debcells(29, 1, 24) = 1
  !debSidesX(24, 1) = 1
  debSidesY(19, 96, :) = 1
  !debSidesZ(24, 1, 29) = 1

  return

100 write(*,*) "ERROR. file 'data/WS/bottom_depth.dat' not found"; call avost
101 write(*,*) "ERROR. file 'data/WS/white_sea_map-abc.dat not found"; call avost
102 write(*,*) "ERROR. topex file not found"; call avost

end subroutine SetupTask_200

!-----------------------------------------------------------------------------------------------------------------

! установка спец-ГУ га X-грань
subroutine SetupWSxBC_200(i,j)
  use variables
  implicit none
  integer:: i, j

  integer:: tpl, tpr

  tpl = c_type(i-1,j)                                           ! тип ячейки слева
  tpr = c_type(i  ,j)                                           ! тип ячейки справа

  if(tpl==CELL_DELETED_SPEC .or. tpr==CELL_DELETED_SPEC) then   ! удалённая спец-ячейка
    ! грань принадлежит одному из трёх втоков; разбираемся - какому
    if(j<10) then                                               ! р. Онега
      fx_type(i,j) = BC_IN                                      ! .. постоянные по времени ГУ
      fx_type(i,j) = BC_STICK                                   ! .. заглушка - стена
      fx_bc(i,j) = 1                                            ! .. индекс данных ГУ
    else if(j<50) then                                          ! р. Северная Двина
      fx_type(i,j) = BC_IN                                      ! .. постоянные по времени ГУ
      fx_type(i,j) = BC_STICK                                   ! .. заглушка - стена
      fx_bc(i,j) = 2                                            ! .. индекс данных ГУ
    else                                                        ! Горло БМ
      fx_type(i,j) = BC_IN_T                                    ! .. переменные по времени ГУ
      fx_bc(i,j) = 3                                            ! .. индекс данных ГУ
      nBcData = nBcData + 1                                     ! счётчик граней в Горле
    end if
  else
    write(*,*) "ERROR in SetupWSxBC_200. i,j:", i, j
    call avost
  end if
end subroutine

!-----------------------------------------------------------------------------------------------------------------

! установка спец-ГУ га Y-грань
subroutine SetupWSyBC_200(i,j)
  use variables
  implicit none
  integer:: i, j

  integer:: tpl, tpr

  tpl = c_type(i,j-1)                                           ! тип ячейки слева
  tpr = c_type(i,j  )                                           ! тип ячейки справа

  if(tpl==CELL_DELETED_SPEC .or. tpr==CELL_DELETED_SPEC) then   ! удалённая спец-ячейка
    ! грань принадлежит одному из трёх втоков; разбираемся - какому
    if(j<10) then                                               ! р. Онега
      fy_type(i,j) = BC_IN                                      ! .. постоянные по времени ГУ
      fy_type(i,j) = BC_STICK                                   ! .. заглушка - стена
      fy_bc(i,j) = 1                                            ! .. номер ГУ
    else if(j<50) then                                          ! р. Северная Двина
      fy_type(i,j) = BC_IN                                      ! .. постоянные по времени ГУ
      fy_type(i,j) = BC_STICK                                   ! .. заглушка - стена
      fy_bc(i,j) = 2                                            ! .. номер ГУ
    else                                                        ! Горло БМ
      fy_type(i,j) = BC_IN_T                                    ! .. переменные по времени ГУ
      fy_bc(i,j) = 3                                            ! .. индекс данных ГУ
      nBcData = nBcData + 1                                     ! счётчик граней в Горле
    end if
  else
    write(*,*) "ERROR in SetupWSyBC_200. i,j:", i, j
    call avost
  end if
end subroutine

!--------------------------------------------------------------------------------------------

subroutine GeoCoords2IJ(lat, lon, i, j)
  use variables
  implicit none

  real(R8):: lat, lon
  integer:: i, j

  i = floor((lat - gxmin) / lgradx) + 1
  j = floor((lon - gymin) / lgrady) + 1
end subroutine
