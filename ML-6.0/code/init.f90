!Тип ГУ на боковых границах. 1: периодичные, 2: жесткая стенка, 3: in/out для теста Steady solution with shock
subroutine Init
  use variables
  implicit none

  logical:: isabcini1 = .true.   ! экспонентациальное увеличение толщин слоёв в тесте Залесного
  logical:: isusesavingdata = .false.!.true.   ! использовать ли сохраненные данные
  logical:: isseriadtest  !В Лабораторных тестах серии D используется 3 слоя
  logical:: isfirsttimehere = .true.
  integer:: ic, jc, is, iz, i, j, k, i1, i2, j1, j2, il, ir, jl, jr, ccr, ccl, kb, kt
  integer:: loadstep, loadiout
  real(R8):: loadtime
  real(R8):: c, zi, hi, ui, vi
  real(R8), allocatable :: htmp(:), lam(:)
  real(R8):: xci, hci, bci, tci, zl, zr, xl, xr, xi, dhi, dx, dy, ht
  real(R8):: alfah, betah, arg, cosa, deltop, delbot, dispan, dmu_0, cf
  real(R8):: dzbot, dztop, dzbotcoef, dztopcoef, rbot, rtop, xbot, xtop, yci, ybot, ytop
  real(R8):: p_0, t_0, z0, r0, r_0, ri, rho1, rho2, u_0, h_0, stepen, cx, cy, y1, y2, initrho, t_surf
  real(R8):: xlockgate, zclayer, zlevel, zlevel2, rhomid, zpt, rho_0

  cfSurfFrictionU = 0.

  ! *************************************************************
  ! Выбор теста
  ! -1) 2D тест Залесного про экватор c кориолисом
  ! -2) Steady solution with shock
  ! -4) Тест 5.1 из статьи J.Sainte-Marie Analytical solutions... 2013
  ! -5) Проточный тест на неровном дне с периодическими ГУ по оси X
  ! -6) Проточный тест на неровном дне с периодическими ГУ по оси Y
  ! -7) Проточный тест на неровном дне с периодическими ГУ по двум осям
  ! -8) Вихрь
  ! -9) Курганов
  ! 10 Бугор на поверхности/дне по X с нулевыми скоростями, период. ГУ
  !
  ! -12 Струя в периодической области с Кориолисом
  ! -13 Лабораторные тесты из статьи C. GLADSTONE "An experimental investigation of density-stratified inertial gravity currents"

  ! -200 Белое море

  ! -998 тест мелой воды
  ! -999 тест мелой воды (вихрь)

  test_numb = 8

  ! что печатать: 1 - печать, 0 - не печатать
  print_height = 1
  print_velocity = 1
  print_velocity_by_color = 1
  print_density = 1
  print_debug = 0

  isexchangeph1 = .true.                  ! делать ли обмен после фазы 1
  isimplexchange = .false.                ! неявный обмен консервативных переменных
  isimplexchanges2 = .false.               ! делать ли обмен потоковых переменных
  sig_ex = 0.
  isfixcover = .false.                    ! расчет с фиксированной крышкой
  isimplcoriolis = .true.
  isconsfilterz0 = .false.                ! фильтр консервативных уровней поверхности

  coriolis = 0.
  sound0 = 100.

  cft_rho0 = 1.d0; cft_alpha = 0.002d0    ! параметры в уравнении состояния rho(T)

  sig_0 = 0; sig = 0.

  maxt = 0.
  tprint = 0
  nt = 0
  nprint = 0

  select case (test_numb)
    case ( 8);  call SetupTask_8          ! вихрь
    case (10);  call SetupTask_10         ! горбы

    case default
      write(*,*) "Wrong test_numb"
      stop
  end select

  write(*,"(4(a,i0),3(a,1p,g12.4))") &
    "Test ",test_numb," cells: ",ncx,"x",ncy,"x",ncz,". Size: ",dl,"x",dw,". CFL=",cfl

  call dims   ! обезразмеривание

  cs_sound = sound0 / sqrt(rho0)
  sound0_2 = sound0**2

  !-- готовим списки объектов сетки: -------------------------------------------------

  ! подсчет количества рабочих ячеек:
  ncells = 0
  do i=1, ncx
    do j=1, ncy
      if(c_type(i,j)>=0) ncells = ncells + 1
    end do
  end do
  ! список ячеек:
  allocate(cells.ptr(ncells))
  ic = 1
  do i=1, ncx
    do j=1, ncy
      if(c_type(i,j)>=0) then
        cells.ptr(ic).i = i
        cells.ptr(ic).j = j
        ic = ic + 1
      end if
    end do
  end do

  ! подсчет количества рабочих X-граней:
  nfx_sides = 0
  do i=1, nxx
    do j=1, nxy
      if(fx_type(i,j)>=0) then
        k = fx_type(i,j)                                ! тип ГУ на грани
        nfx_sides(k) = nfx_sides(k) + 1                 ! счётчик граней с данным ГУ
      end if
    end do
  end do
  ! список X-граней:
  do k = 0, nbc                                         ! цикл по типам ГУ
    if(nfx_sides(k)==0) cycle                           ! если такого типа ГУ в задаче нет
    allocate(fx_sides(k).ptr(nfx_sides(k)))             ! выделяем место под список

    is = 1
    do i=1, nxx
      do j=1, nxy
        if(fx_type(i,j)==k) then                       ! если грань нужного типа
          fx_sides(k).ptr(is).i = i                     ! заносим индексы грани в соответствующий список
          fx_sides(k).ptr(is).j = j
          is = is + 1
        end if
      end do
    end do
  end do

  ! подсчет количества рабочих Y-граней:
  nfy_sides = 0
  do i=1, nyx
    do j=1, nyy
      if(fy_type(i,j)>=0) then
        k = fy_type(i,j)                               ! тип ГУ на грани
        nfy_sides(k) = nfy_sides(k) + 1                 ! счётчик граней с данным ГУ
      end if
    end do
  end do
  ! список Y-граней:
  do k = 0, nbc                                         ! цикл по типам ГУ
    if(nfy_sides(k)==0) cycle                           ! если такого типа ГУ в задаче нет
    allocate(fy_sides(k).ptr(nfy_sides(k)))             ! выделяем место под список

    is = 1
    do i=1, nyx
      do j=1, nyy
        if(i==44 .and. j==1) &
          debdum=1
        if(fy_type(i,j)==k) then                       ! если грань нужного типа
          fy_sides(k).ptr(is).i = i                     ! заносим индексы грани в соответствующий список
          fy_sides(k).ptr(is).j = j
          is = is + 1
        end if
      end do
    end do
  end do

  !.......................................................................

  ! шаги сетки:
  do i=1,ncx; c_dx(i) = x(i+1) - x(i); end do
  do i=1,ncy; c_dy(i) = y(i+1) - y(i); end do

  ! довычисляем параметры в ячейках:
  c_zt(:,:) = fz_z(:,:,1)
  c_b(:,:) = fz_z(:,:,nz)
  c_z0(:,:) = c_b(:,:) + c_h0(:,:)                              ! уровень поверхности в МВ
  fzn_z(:,:,nz) = fz_z(:,:,nz)

  ! на гранях:
  fxn_z(:,:,nz) = fx_z(:,:,nz)
  fyn_z(:,:,nz) = fy_z(:,:,nz)

  ! параметры на гранях, как полусумма значений в ячейках:
  do i=1,nxx; do j=1,nxy                                        ! для X-граней
    ccl = i - 1                                                  ! ячейка слева
    ccr = i                                                      ! ячейка справа

    ! на граничных гранях корректируем индексы ячеек:
    if(fx_type(i,j)==BC_PERIODIC) then                          ! для периодических граней
      if(i==1) ccl = ncx
      if(i==nxx) ccr = 1
    else if(fx_type(i,j)>0) then                                ! любая другая граничная грань
      if(i==1) then; ccl = ccr
      else if(c_type(ccl,j)<0) then; ccl = ccr; end if
      if(i==nxx) then; ccr = ccl
      else if(c_type(ccr,j)<0) then; ccr = ccl; end if
    end if

    call assertn(c_type(ccl,j)>=0 .and. c_type(ccr,j)>=0, "Init-sides-X. Ошибка сетки", i, j, -1)

    fx_z(i,j,:) = (fz_z(ccl,j,:) + fz_z(ccr,j,:)) / 2.
    fx_b(i,j) = fx_z(i,j,nz)
    !fx_zt(i,j) = fx_z(i,j,1)

    fx_u(i,j,:) = (c_u(ccl,j,:) + c_u(ccr,j,:)) / 2.
    fx_v(i,j,:) = (c_v(ccl,j,:) + c_v(ccr,j,:)) / 2.
    fx_w(i,j,:) = (c_w(ccl,j,:) + c_w(ccr,j,:)) / 2.
    fx_rho(i,j,:) = (c_rho(ccl,j,:) + c_rho(ccr,j,:)) / 2.
    fx_teta(i,j,:) = (c_teta(ccl,j,:) + c_teta(ccr,j,:)) / 2.

    ! для мелкой воды:
    fx_h0(i,j) = (c_h0(ccl,j) + c_h0(ccr,j)) / 2.
    fx_u0(i,j) = (c_u0(ccl,j) + c_u0(ccr,j)) / 2.
    fx_v0(i,j) = (c_v0(ccl,j) + c_v0(ccr,j)) / 2.
    fx_rho0(i,j) = (c_rho0(ccl,j) + c_rho0(ccr,j)) / 2.

    fx_z0(i,j) = fx_b(i,j) + fx_h0(i,j)

  end do; end do

  do i=1,nyx; do j=1,nyy                                        ! для Y-граней
    ccl = j - 1                                                  ! ячейка слева
    ccr = j                                                      ! ячейка справа

    ! на граничных гранях корректируем индексы ячеек:
    if(fy_type(i,j)==BC_PERIODIC) then                          ! для периодических граней
      if(j==1) ccl = ncy
      if(j==nyy) ccr = 1
    else if(fy_type(i,j)>0) then                                ! любая другая граничная грань
      if(j==1) then; ccl = ccr
      else if(c_type(i,ccl)<0) then; ccl = ccr; end if
      if(j==nyy) then; ccr = ccl
      else if(c_type(i,ccr)<0) then; ccr = ccl; end if
    end if

    call assertn(c_type(i,ccl)>=0 .and. c_type(i,ccr)>=0, "Init-sides-Y. Ошибка сетки", i, j, -1)

    fy_z(i,j,:) = (fz_z(i,ccl,:) + fz_z(i,ccr,:)) / 2.
    fy_b(i,j) = fy_z(i,j,nz)

    fy_u(i,j,:) = (c_u(i,ccl,:) + c_u(i,ccr,:)) / 2.
    fy_v(i,j,:) = (c_v(i,ccl,:) + c_v(i,ccr,:)) / 2.
    fy_w(i,j,:) = (c_w(i,ccl,:) + c_w(i,ccr,:)) / 2.
    fy_rho(i,j,:) = (c_rho(i,ccl,:) + c_rho(i,ccr,:)) / 2.
    fy_teta(i,j,:) = (c_teta(i,ccl,:) + c_teta(i,ccr,:)) / 2.

    ! для мелкой воды:
    fy_h0(i,j) = (c_h0(i,ccl) + c_h0(i,ccr)) / 2.
    fy_u0(i,j) = (c_u0(i,ccl) + c_u0(i,ccr)) / 2.
    fy_v0(i,j) = (c_v0(i,ccl) + c_v0(i,ccr)) / 2.
    fy_rho0(i,j) = (c_rho0(i,ccl) + c_rho0(i,ccr)) / 2.

    fy_z0(i,j) = fy_b(i,j) + fy_h0(i,j)

  end do; end do

  ! на горизонтальных гранях:
  do ic=1,ncells
    i = cells.ptr(ic).i
    j = cells.ptr(ic).j

    do k=1,nz
      kt = k - 1                          ! слой выше грани
      kb = k                              ! слой ниже грани
      if(k==1) kt = kb
      if(k==nz) kb = kt

      fz_u(i,j,k) = (c_u(i,j,kt) + c_u(i,j,kb)) / 2.
      fz_v(i,j,k) = (c_v(i,j,kt) + c_v(i,j,kb)) / 2.
      fz_w(i,j,k) = (c_w(i,j,kt) + c_w(i,j,kb)) / 2.
      fz_rho(i,j,k) = (c_rho(i,j,kt) + c_rho(i,j,kb)) / 2.
      fz_dteta(i,j,k) = (c_teta(i,j,kt) + c_teta(i,j,kb)) / 2.
    end do
  end do

  !.......................................................................

  ! для вычисления потоков, определяем переменные на грани для t=t[n+1]
  fxn_h0 = fx_h0; fxn_u0 = fx_u0; fxn_v0 = fx_v0; fxn_rho0 = fx_rho0
  fyn_h0 = fy_h0; fyn_u0 = fy_u0; fyn_v0 = fy_v0; fyn_rho0 = fy_rho0

  ! вертикальная скорость в мелкой воде:
  call CalcSwW

  ! вычисляем z с точкой; вертикальные компоненты заданы в Н/У:
  do ic=1,ncells
    i = cells.ptr(ic).i
    j = cells.ptr(ic).j

    ! z с точкой на поверхности:
    zpt = c_dw(i,j,1) - &
               (fx_du(i+1,j,1) + fx_du(i,j,1)) / 2. * (fx_z(i+1,j,1) - fx_z(i,j,1)) / c_dx(i) - &
               (fy_dv(i,j+1,1) + fy_dv(i,j,1)) / 2. * (fy_z(i,j+1,1) - fy_z(i,j,1)) / c_dy(j)

    ! z с точкой на разделах слоёв в ячейках:
    do k = 1, nz                                              ! цикл по разделам слоёв
      zi = fz_z(i,j,k)
      fz_zp(i,j,k) = zpt * (zi - c_b(i,j)) / (c_zt(i,j) - c_b(i,j))
    end do
  end do

  ! дельта-переменные:

  do i=1,nxx; do j=1,nxy                          ! цикл по X-граням
    if(fx_type(i,j)<0) cycle                      ! пропускаем не рабочие грани
    do k=1,ncz
      fx_w0(i,j,k) = fx_GetW0(i,j,k)              ! Вертикальная скорость МВ на X-гранях
      fx_dteta(i,j,k) = fx_teta(i,j,k) - teta0
      fx_drho(i,j,k) = fx_rho(i,j,k) - rho0
      fx_du(i,j,k) = fx_u(i,j,k) - fx_u0(i,j)
      fx_dv(i,j,k) = fx_v(i,j,k) - fx_v0(i,j)
      fx_dw(i,j,k) = fx_w(i,j,k) - fx_w0(i,j,k)
    end do
  end do; end do

  do i=1,nyx; do j=1,nyy                          ! цикл по Y-граням
    if(fy_type(i,j)<0) cycle                      ! пропускаем не рабочие грани
    do k=1,ncz
      fy_w0(i,j,k) = fy_GetW0(i,j,k)              ! Вертикальная скорость МВ на X-гранях
      fy_dteta(i,j,k) = fy_teta(i,j,k) - teta0
      fy_drho(i,j,k) = fy_rho(i,j,k) - rho0
      fy_du(i,j,k) = fy_u(i,j,k) - fy_u0(i,j)
      fy_dv(i,j,k) = fy_v(i,j,k) - fy_v0(i,j)
      fy_dw(i,j,k) = fy_w(i,j,k) - fy_w0(i,j,k)
    end do
  end do; end do

  do i=1,ncx; do j=1,ncy                          ! цикл по ячейкам
    if(c_type(i,j)<0) cycle                       ! пропускаем не рабочие ячейки
    do k=1,nz
      fz_w0(i,j,k) = fz_GetW0(i,j,k)              ! Вертикальная скорость МВ на Z-гранях
      fz_dteta(i,j,k) = fz_teta(i,j,k) - teta0
      fz_drho(i,j,k) = fz_rho(i,j,k) - rho0
      fz_du(i,j,k) = fz_u(i,j,k) - c_u0(i,j)
      fz_dv(i,j,k) = fz_v(i,j,k) - c_v0(i,j)
      fz_dw(i,j,k) = fz_w(i,j,k) - fz_w0(i,j,k)
    end do
    do k=1,ncz
      c_w0(i,j,k) = c_GetW0(i,j,k)               ! Вертикальная скорость МВ в слой-ячейках
      c_dteta(i,j,k) = c_teta(i,j,k) - teta0
      c_drho(i,j,k) = c_rho(i,j,k) - rho0
      c_du(i,j,k) = c_u(i,j,k) - c_u0(i,j)
      c_dv(i,j,k) = c_v(i,j,k) - c_v0(i,j)
      c_dw(i,j,k) = c_w(i,j,k) - c_w0(i,j,k)
    end do
  end do; end do

  printproc = 0
  ttime=0
  istep=0

  !     *****************************************************************
  !     Вычисление интегральных параметров
  vol=0  ! Начальный полный объем
  totmass_0=0  ! Начальная полная масса

  ekin_0=0  ! Начальная суммарная кинетическая энергия
  epot_0=0  ! Начальная суммарная потенциальная энергия
  etot_0=0   ! Начальная полная энергия

  do k=1,nz-1
    do j=1,ny-1
      do i=1,nx-1
        dx=x(i+1)-x(i)
        dy=y(j+1)-y(j)

        vol_0=vol_0+dx*(fz_z(i,j,k)-fz_z(i,j,k+1))
        totmass_0= totmass_0+(rho0+c_drho(i,j,k))*dx*dy*(fz_z(i,j,k)-fz_z(i,j,k+1))

        ekin_0=ekin_0+0.5*(rho0+c_drho(i,j,k))*dx*dy*(fz_z(i,j,k)-fz_z(i,j,k+1))*(c_du(i,j,k)**2+c_dv(i,j,k)**2)

        epot_0=epot_0+0.5*dx*dy*g*(rho0+c_drho(i,j,k))*(fz_z(i,j,k)**2-fz_z(i,j,k+1)**2)
      enddo
    enddo
  enddo
  etot_0=ekin_0+epot_0

end subroutine Init
