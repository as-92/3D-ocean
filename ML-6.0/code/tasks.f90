! начальные данные, спецефичные для задач
!---------------------------------------------------------------------------------------------------------
! Ниже для каждой из задач необходимо задать следующие значения в ячейках:
!   c_b(:,:)    - z-координата дна
! для мелкой воды:
!   c_h0(:,:)   - глубина (расстояние от поверхности до дна)
!   c_u0(:,:)   - X-скорость
!   c_v0(:,:)   - Y-скорость
!   c_rho0(:,:)
!   c_wb0(:,:)     - вычисляемый параметр
!   c_wt0(:,:)     - вычисляемый параметр
! для слоёв, сверху вниз:
!   fz_z(:,:,:) - уровни раздела слоёв
!   c_dteta(:,:,:) - отклонение фиктивной плотности от teta0
!   c_drho(:,:,:)  - отклонение реальной плотности от rho0
!   c_du(:,:,:)    - отклонение X-скорости от c_u0(:,:)
!   c_dv(:,:,:)    - отклонение Y-скорости от c_v0(:,:)
!   c_dw(:,:,:)    - отклонение Z-скорости от c_GetW(:,:,:)
! Также, для каждой ячейки и грани необходимо задать c_type, fx_type, fy_type

! Вихрь
subroutine SetupTask_8
  use variables
  implicit none

  integer:: i, ic, il, ir, iz, j, jc, jl, jr, k
  real(R8):: alfah, bci, betah, cosa, delbot, deltop, dmu_0, dx, dy
  real(R8):: dzbot, dzbotcoef, dztop, dztopcoef, h_0
  real(R8):: p_0, r0, r_0, rbot, ri, rtop, stepen, tci, u_0
  real(R8):: ui, vi, zi, hi, xbot, xci, xtop, yci, z0

  dl=2.; dw=2.                  ! размер области
  nx=51; ny=51; nz=11;           ! количество узлов сетки
  bctypex = 1; bctypey = 1      ! тип ГУ на боковых границах

  call allocate

  ! отладка:
  debstartstep = 0; debendstep = -5
  debstartz = 1; debendz = 400

  !debcells(1, 1) = 1
  !debcells(ncx, 1) = 1
  !debSidesX(1, 1) = 1
  !debSidesX(nxx, 1) = 1

  !debSidesY(44, 1) = 1
  !debsidesY(44, 2) = 1

  !debcells(31, 1) = 1
  !debSidesX(31, 1) = 1
  !debSidesZ(31, 1) = 1

#if 0
  nt=10 ! полное число шагов по времени
  nprint=1 ! интервал печати
#else
  maxt=10.
  tprint=maxt/200.
#endif

  CFL = 0.3 ! Число Куранта

  pan=0.0 !- параметр Паниковского

  g = 1
  rho0 = 1.

  ! Расчет начальной сетки
  x(1) = -dl / 2.
  dx = dl / (nx-1)
  do i=1,ncx; x(i+1) = x(i) + dx; end do                    ! координаты узлов сетки

  y(1) = -dw / 2.
  dy = dw / (ny-1)
  do j=1,ny-1; y(j+1) = y(j) + dy; end do                   ! координаты узлов сетки

  z0=1.                                   ! фоновый уровень поверхности
  b0=0.                                   ! фоновый уровень дна
  h_0 = z0 - b0

  alfaH = 0.404
  betaH = 0.4
  r_0= 0.1

  ! начальные данные на сетке; мелкая вода:
  do jc=1, ncy
    do ic=1, ncx
      xci = (x(ic) + x(ic+1)) / 2.               ! x-координата ячейки
      yci = (y(jc) + y(jc+1)) / 2.               ! y-координата ячейки
      ri = sqrt(xci**2+yci**2)
      stepen=betaH*(1.-(ri**2/r_0**2))
      ui =  alfaH/r_0*exp(stepen)*yci
      vi = -alfaH/r_0*exp(stepen)*xci
      hi = h_0 - alfaH**2*exp(2*stepen)/(4.*g*betaH)    ! глубина

      c_u0(ic,jc) = ui                            ! скорость
      c_v0(ic,jc) = vi
      c_b(ic,jc) = b0                             ! дно
      zi = b0 + hi                                ! поверхность
      ! равномерная сетка по вертикали от поверхности до дна:
      do iz = 1, nz
        fz_z(ic, jc, iz) = zi - (iz - 1.) * hi / (nz - 1.)
      end do
    end do
  end do

  c_h0(:,:) = fz_z(:,:,1) - fz_z(:,:,nz)
  c_rho0 = rho0

  ! константы:
  rho0 = 1.
  teta0 = 1.
  sound0 = 10.             ! псевдо-скорость звука:

  ! начальные данные на сетке; слои:
  c_rho = rho0
  c_teta = teta0
  do i=1,ncx; do j=1,ncy; do k=1,ncz
    c_u(i,j,k) = c_u0(i,j)
    c_v(i,j,k) = c_v0(i,j)
    c_w(i,j,k) = 0.
  end do; end do; end do

  ! For WB test
#if 0
  do k=1,nz-1
    do j=1,ny-1
      do i=1,nx-1
        c_du(i,j,k)=0.25*( fx_du(i,j,k)+fx_du(i+1,j,k)+fy_du(i,j,k)+fy_du(i,j+1,k))
        c_dv(i,j,k)=0.25*(fx_dv(i,j,k)+fx_dv(i+1,j,k)+fy_dv(i,j,k)+fy_dv(i,j+1,k))
      enddo
    enddo
  enddo
#endif

end

!================================================================================================================

! Бугор на дне и на поверхности:
subroutine SetupTask_10
  use variables
  implicit none

  integer:: i, ic, il, ir, iz, j, jc, jl, jr, k
  real(R8):: alfah, bci, betah, cosa, delbot, deltop, dmu_0, dx, dy
  real(R8):: dzbot, dzbotcoef, dztop, dztopcoef, h_0, hi
  real(R8):: p_0, r0, r_0, rbot, ri, rtop, stepen, tci, u_0
  real(R8):: ui, vi, xbot, ybot, xci, xtop, ytop, yci, z0

#if 0
  dl=20.; dw=1.                 ! размер области
  nx=65; ny=2; nz=2;           ! количество узлов сетки
#elif 1
  dl=1.; dw=20.                 ! размер области
  nx=2; ny=65; nz=2;           ! количество узлов сетки
#else
  dl=20.; dw=20.                 ! размер области
  nx=45; ny=45; nz=2;           ! количество узлов сетки
#endif
  bctypex = 1; bctypey = 1      ! тип ГУ на боковых границах

  call allocate

  ! отладка:
  debstartstep = 0; debendstep = 5
  debstartz = 1; debendz = 400

  !debcells(1, 1) = 1
  !debcells(ncx, 1) = 1
  !debSidesX(1, 1) = 1
  !debSidesX(nxx, 1) = 1

  !debCells(20, 12) = 1
  !debCells(12, 20) = 1
  !debSidesX(20, 12) = 1
  !debsidesY(12, 20) = 1

  !debcells(31, 1) = 1
  !debSidesX(31, 1) = 1
  !debSidesZ(31, 1) = 1

#if 0
  nt=200 ! полное число шагов по времени
  nprint=nt/200 ! интервал печати
#else
  maxt=50.
  tprint=maxt/200.
#endif

  CFL = 0.3 ! Число Куранта

  pan=0.0 !- параметр Паниковского

  g = 1
  rho0 = 1.

  ! Расчет начальной сетки
  x(1) = -dl / 2.
  dx = dl / (nx-1)
  do i=1,ncx; x(i+1) = x(i) + dx; end do                    ! координаты узлов сетки

  y(1) = -dw / 2.
  dy = dw / (ny-1)
  do j=1,ny-1; y(j+1) = y(j) + dy; end do                   ! координаты узлов сетки

  z0=1.                                   ! фоновый уровень поверхности
  b0=0.                                   ! фоновый уровень дна

  ! для поверхности:
  dztopcoef = 1. !0.01                        ! коэффициент подъема поверхности (0<cf<1)
  rtop = max(dl,dw)/4.                            ! ширина возмущения поверхности
  xtop = x(1) + 0.5 * dl                  ! координата центра  на поверхности
  ytop = y(1) + 0.5 * dw
  dztop = dztopcoef * (z0 - b0)           ! вертикальный размер возмущения

  ! возмущение дна:
  dzbotcoef = 0.                          ! коэффициент подъема дна (0<cf<1)
  rbot = max(dl,dw)/4.                            ! ширина возмущения дна
  xbot = x(1) + 0.5 * dl                  ! координата возмущения
  ybot = y(1) + 0.5 * dw
  dzbot = dzbotcoef * (z0 - b0)           ! вертикальная величина колебания дна

  do jc=1, ncy
    yci = (y(jc) + y(jc+1)) / 2.                 ! y-координата ячейки
    do ic=1, ncx
      xci = (x(ic) + x(ic+1)) / 2.               ! x-координата ячейки

      ! бугор на поверхности:
      deltop = 0.
      r0 = max(dl, dw)
      do i=-1,1; do j=-1,1
        ri = sqrt((xci+dl*i-xtop)**2 + (yci+dw*i-ytop)**2)
        if(ri<r0) r0 = ri
      end do; end do

      if(r0<rtop) then                                                 ! если в пределах возмущения
        ri = r0 / rtop                                                ! нормируем на 0..1
        cosa = cos(pi * ri)
        deltop = dztop * ((1.d0 + cosa) / 2.d0)
      end if

      ! бугор на дне:
      delbot = 0.
      r0 = max(dl, dw)
      do i=-1,1; do j=-1,1
        ri = sqrt((xci+dl*i-xbot)**2 + (yci+dw*i-ybot)**2)
        if(ri<r0) r0 = ri
      end do; end do

      if(r0<rbot) then                                                 ! если в пределах возмущения
        ri = r0 / rbot                                                 ! нормируем на 0..1
        delbot = dzbot * ((1. + cos(pi * ri)) / 2.)
      end if

      bci = b0 + delbot                             ! уровень дна в ячейке
      tci = z0 + deltop                             ! уровень поверхности в ячейке

      ! равномерная сетка по вертикали от поверхности до дна:
      do iz = 1, nz
        fz_z(ic, jc, iz) = tci - (iz - 1.) * (tci - bci)/(nz - 1.)
      end do
    end do
  end do

  ! начальные данные на сетке; мелкая вода:
  u_0 = 0 !0.1d0                                                            ! скорость потока
  c_h0(:,:) = fz_z(:,:,1) - fz_z(:,:,nz)
  c_rho0 = rho0
  dmu_0 = u_0 * c_h0(1,1)                                             ! импульс в одной из ячеек
  c_u0 = dmu_0 / c_h0                                                 ! скорость с учётом равного импульса
  c_v0 = 0.

  ! константы:
  rho0 = 1.
  teta0 = 1.
  sound0 = 2.             ! псевдо-скорость звука:

  ! начальные данные на сетке; слои:
  c_rho = rho0
  c_teta = teta0
  do i=1,ncx; do j=1,ncy; do k=1,ncz
    c_u(i,j,k) = c_u0(i,j)
    c_v(i,j,k) = c_v0(i,j)
    c_w(i,j,k) = 0.
  end do; end do; end do

  ! For WB test
#if 0
  do k=1,nz-1
    do j=1,ny-1
      do i=1,nx-1
        c_du(i,j,k)=0.25*( fx_du(i,j,k)+fx_du(i+1,j,k)+fy_du(i,j,k)+fy_du(i,j+1,k))
        c_dv(i,j,k)=0.25*(fx_dv(i,j,k)+fx_dv(i+1,j,k)+fy_dv(i,j,k)+fy_dv(i,j+1,k))
      enddo
    enddo
  enddo
#endif

end

!------------------------------------------------------------------------------------------------

subroutine SetupTask_998
  use variables
  implicit none

#if 0
  ! сетка:
  nx=51  ! количество узлов сетки по горизонтали
  ny=51
  nz=2     ! число расчетных узлов по в вертикали
!      nxcells = nx - 1;   nycells = ny - 1
!      nxsides = nx    ;   nysides = ny

  dl=2.;              dw=2.                            ! размер расчётной области
  xmin = -dl/2.;      xmax = dl/2. ! [cм]
  ymin = -dw/2.;      ymax = dw/2.

  bctypex = 2 ! Тип ГУ на боковых границах. 1: периодичные, 2: жесткая стенка, 3:
  bctypey = 2 ! Тип ГУ на боковых границах. 1: периодичные

  call allocate

  ! отладка:
  debstartstep = 0; debendstep = -10
  debstartz = 1; debendz = 400
  debcells(25, 24) = 1
  debcells(25, 25) = 1
  debcells(25, 26) = 1
  debcells(25, 27) = 1
  !debSidesX(25, 48) = 1
  !debSidesX(26, 48) = 1
  !debSidesY(26, 27) = 1
  debsidesy(25, 26) = 1

  pan=0.0 !- параметр Паниковского

  alfah = 0.404
  betah = 0.4
  r_0= 0.1
  p_0 = 1.
  h_0 = 1.                  ! фоновая глубина
  g= 1.!9.81 !! гравитационная постоянная

  ! ***************************************************************
  !     ЧИСЛО ШАГОВ ПО ВРЕМЕНИ И ИНТЕРВАЛ ПЕЧАТИ

  nt = 400                   ! полное число шагов по времени
  nprint = 1
  !maxt = 50.; tprint = maxt / 250.

  cfl=0.3d0 ! Число Куранта

  sigma2ph=1.

  !     *******************************************************
  !     Расчет начальной сетки
  !     *******************************************************

  do i=1,nx
    x(i) = xmin + (i-1.) * (xmax-xmin) / (nx-1.)
  enddo

  do i=1,ny
    y(i) = ymin + (i-1.) * (ymax-ymin) / (ny-1.)
  enddo

  kface=nz

  z0=0.
  b0 = -h_0

  do jc=1, ncy
    do ic=1, ncx
      xci = (x(ic) + x(ic+1)) / 2.               ! x-координата ячейки
      yci = (y(jc) + y(jc+1)) / 2.               ! y-координата ячейки
      ui =  0.
      vi = 0
      hi = h_0                                    ! глубина
      if(xci<(xmin+xmax)/2.) ui = 0.1d0
      !if(yci<(ymin+ymax)/2.) vi = 0.1d0

      c_du(ic,jc,:) = ui                           ! скорость
      c_dv(ic,jc,:) = vi
      c_b(ic,jc) = b0                             ! дно
      c_zt(ic, jc) = b0 + hi                      ! поверхность
      ! равномерная сетка по вертикали от поверхности до дна:
      do iz = 1, nz
        fz_z(ic, jc, iz) = c_zt(ic, jc) - (iz - 1.) * hi / (nz - 1.)
      end do

      ! начальные данные для мелкой воды:
      c_h0 = hi
      c_u0(ic,jc) = ui
      c_v0(ic,jc) = vi
      c_rho0(ic,jc)=1.

    end do
  end do

  ! сетка на гранях по оси x (поверхность и дно):
  ! внутренние грани:
  do j=1,ny-1
    do i=2,nx-1
      fx_z(i,j,1)=0.5*(fz_z(i,j,1)+fz_z(i-1,j,1))
      fx_z(i,j,nz)=0.5*(fz_z(i,j,nz)+fz_z(i-1,j,nz))
    enddo

    ! граничные грани:
    fx_z(1,j, nz) = 0.5*(fz_z(1,j,nz)+fz_z(nx-1,j,nz))
    fx_z(1,j, 1) = 0.5*(fz_z(1,j,1)+fz_z(nx-1,j,1))
    fx_z(nx,j, nz) = fx_z(1,j, nz)
    fx_z(nx,j, 1) = fx_z(1,j, 1)

    ! промежуточные разделы слоёв:
    do k = 2, nz-1
      fx_z(:,j,k) = fx_z(:,j,1)-(fx_z(:,j,1)-fx_z(:,j,nz))/(nz-1.)*(k-1.)
    end do

  enddo
  !hx = fx_z(:,:, 1)
  fx_b = fx_z(:,:, nz)
  ! сетка на гранях по оси y (поверхность и дно):
  ! внутренние грани:
  do i=1,nx-1
    do j=2,ny-1
      fy_z(i,j,1)=0.5*(fz_z(i,j,1)+fz_z(i,j-1,1))
      fy_z(i,j,nz)=0.5*(fz_z(i,j,nz)+fz_z(i,j-1,nz))
    enddo

    ! граничные грани:
    fy_z(i,1, nz) = 0.5*(fz_z(i,1,nz)+fz_z(i,ny-1,nz))
    fy_z(i,1, 1) = 0.5*(fz_z(i,1,1)+fz_z(i,ny-1,1))
    fy_z(i,ny, nz) = fy_z(i,1, nz)
    fy_z(i,ny, 1) = fy_z(i,1, 1)

    ! промежуточные разделы слоёв:
    do k = 2, nz-1
      fy_z(i,:,k) = fy_z(i,:,1)-(fy_z(i,:,1)-fy_z(i,:,nz))/(nz-1.)*(k-1.)
    end do

  enddo
  !HY = fy_z(:,:, 1)
  fy_b = fy_z(:,:, nz)


  ! For WB test
  do k=1,nz
    do j=1,ny-1
      do i=1,nx-1
        fz_z(i,j,k)=0.25*( fx_z(i,j,k)+fx_z(i+1,j,k)+fy_z(i,j,k)+fy_z(i,j+1,k))
      enddo
    enddo
  enddo

  do j=1,ny-1
    do i=1,nx-1
      c_zt(i,j)=fz_z(i,j,1)
      c_b(i,j)=fz_z(i,j,nz)
    enddo
  enddo

  ! потоковые скорости по оси X
  do j=1,ny-1
    do k=1,nz-1
      do i=2,nx-1
        fx_du(i,j,k)=0.5*(c_du(i,j,k)+c_du(i-1,j,k))
        fx_dv(i,j,k)=0.5*(c_dv(i,j,k)+c_dv(i-1,j,k))
      enddo
      fx_du(1,j,k)=0.5*(c_du(1,j,k)+c_du(nx-1,j,k))
      fx_du(nx,j,k)=fx_du(1,j,k)
      fx_dv(1,j,k)=0.5*(c_dv(1,j,k)+c_dv(nx-1,j,k))
      fx_dv(nx,j,k)=fx_dv(1,j,k)
    enddo
  enddo

  ! потоковые скорости по оси Y
  do i=1,nx-1
    do k=1,nz-1
      do j=2,ny-1
        fy_du(i,j,k)=0.5*(c_du(i,j,k)+c_du(i,j-1,k))
        fy_dv(i,j,k)=0.5*(c_dv(i,j,k)+c_dv(i,j-1,k))
      enddo
      fy_du(i,1,k)=0.5*(c_du(i,1,k)+c_du(i,ny-1,k))
      fy_du(i,ny,k)=fy_du(i,1,k)
      fy_dv(i,1,k)=0.5*(c_dv(i,1,k)+c_dv(i,ny-1,k))
      fy_dv(i,ny,k)=fy_dv(i,1,k)
    enddo
  enddo

  ! For WB test
  do k=1,nz-1
    do j=1,ny-1
      do i=1,nx-1
        c_du(i,j,k)=0.25*( fx_du(i,j,k)+fx_du(i+1,j,k)+fy_du(i,j,k)+fy_du(i,j+1,k))
        c_dv(i,j,k)=0.25*(fx_dv(i,j,k)+fx_dv(i+1,j,k)+fy_dv(i,j,k)+fy_dv(i,j+1,k))
      enddo
    enddo
  enddo

  ! на гранях: потоковые переменные мелкой воды:
  do i=1,nxx
    do j=1,nxy
      il = max(i-1,1)
      ir = min(i,ncx)
      fx_h0(i,j) = (c_h0(il,j) + c_h0(ir,j)) / 2.
      fx_u0(i,j) = (c_u0(il,j) + c_u0(ir,j)) / 2.
      fx_v0(i,j) = (c_v0(il,j) + c_v0(ir,j)) / 2.
      fx_rho0(i,j) = (c_rho0(il,j) + c_rho0(ir,j)) / 2.
    end do
  end do
  do i=1,nyx
    do j=1,nyy
      jl = max(j-1,1)
      jr = min(j,ncy)
      fy_h0(i,j) = (c_h0(i,jl) + c_h0(i,jr)) / 2.
      fy_u0(i,j) = (c_u0(i,jl) + c_u0(i,jr)) / 2.
      fy_v0(i,j) = (c_v0(i,jl) + c_v0(i,jr)) / 2.
      fy_rho0(i,j) = (c_rho0(i,jl) + c_rho0(i,jr)) / 2.
    end do
  end do
#endif

end

!------------------------------------------------------------------------------------------------

subroutine SetupTask_999
  use variables
  implicit none

#if 0
  ! сетка:
  nx=51  ! количество узлов сетки по горизонтали
  ny=51
  nz=2     ! число расчетных узлов по в вертикали
!      nxcells = nx - 1;   nycells = ny - 1
!      nxsides = nx    ;   nysides = ny

  dl=2.;              dw=2.                            ! размер расчётной области
  xmin = -dl/2.;      xmax = dl/2. ! [cм]
  ymin = -dw/2.;      ymax = dw/2.

  bctypex = 2 ! Тип ГУ на боковых границах. 1: периодичные, 2: жесткая стенка, 3:
  bctypey = 2 ! Тип ГУ на боковых границах. 1: периодичные

  call allocate

  pan=0.0 !- параметр Паниковского

  alfah = 0.404
  betah = 0.4
  r_0= 0.1
  p_0 = 1.
  h_0 = 1.                  ! фоновая глубина
  g= 1.!9.81 !! гравитационная постоянная

  ! ***************************************************************
  !     ЧИСЛО ШАГОВ ПО ВРЕМЕНИ И ИНТЕРВАЛ ПЕЧАТИ

  !nprint = 1; nt = 400 * nprint                   ! полное число шагов по времени
  maxt = 50.; tprint = maxt / 250.

  cfl=0.4 ! Число Куранта

  sigma2ph=1.

  !     *******************************************************
  !     Расчет начальной сетки
  !     *******************************************************

  do i=1,nx
    x(i) = xmin + (i-1.) * (xmax-xmin) / (nx-1.)
  enddo

  do i=1,ny
    y(i) = ymin + (i-1.) * (ymax-ymin) / (ny-1.)
  enddo

  kface=nz

  z0=0.
  b0 = -h_0

  do jc=1, ncy
    do ic=1, ncx
      xci = (x(ic) + x(ic+1)) / 2.               ! x-координата ячейки
      yci = (y(jc) + y(jc+1)) / 2.               ! y-координата ячейки
      ri = sqrt(xci**2+yci**2)
      stepen=betah*(1.-(ri**2/r_0**2))
      ui =  alfah/r_0*exp(stepen)*yci
      vi = -alfah/r_0*exp(stepen)*xci
      hi = h_0 - alfah**2*exp(2*stepen)/(4.*g*betah)    ! глубина

      c_du(ic,jc,:) = ui                           ! скорость
      c_dv(ic,jc,:) = vi
      c_b(ic,jc) = b0                             ! дно
      c_zt(ic, jc) = b0 + hi                      ! поверхность
      ! равномерная сетка по вертикали от поверхности до дна:
      do iz = 1, nz
        fz_z(ic, jc, iz) = c_zt(ic, jc) - (iz - 1.) * hi / (nz - 1.)
      end do

      ! начальные данные для мелкой воды:
      c_h0 = hi
      c_u0(ic,jc) = ui
      c_v0(ic,jc) = vi
      c_rho0(ic,jc)=1.

    end do
  end do

  ! сетка на гранях по оси x (поверхность и дно):
  ! внутренние грани:
  do j=1,ny-1
    do i=2,nx-1
      fx_z(i,j,1)=0.5*(fz_z(i,j,1)+fz_z(i-1,j,1))
      fx_z(i,j,nz)=0.5*(fz_z(i,j,nz)+fz_z(i-1,j,nz))
    enddo

    ! граничные грани:
    fx_z(1,j, nz) = 0.5*(fz_z(1,j,nz)+fz_z(nx-1,j,nz))
    fx_z(1,j, 1) = 0.5*(fz_z(1,j,1)+fz_z(nx-1,j,1))
    fx_z(nx,j, nz) = fx_z(1,j, nz)
    fx_z(nx,j, 1) = fx_z(1,j, 1)

    ! промежуточные разделы слоёв:
    do k = 2, nz-1
      fx_z(:,j,k) = fx_z(:,j,1)-(fx_z(:,j,1)-fx_z(:,j,nz))/(nz-1.)*(k-1.)
    end do

  enddo
  !hx = fx_z(:,:, 1)
  fx_b = fx_z(:,:, nz)
  ! сетка на гранях по оси y (поверхность и дно):
  ! внутренние грани:
  do i=1,nx-1
    do j=2,ny-1
      fy_z(i,j,1)=0.5*(fz_z(i,j,1)+fz_z(i,j-1,1))
      fy_z(i,j,nz)=0.5*(fz_z(i,j,nz)+fz_z(i,j-1,nz))
    enddo

    ! граничные грани:
    fy_z(i,1, nz) = 0.5*(fz_z(i,1,nz)+fz_z(i,ny-1,nz))
    fy_z(i,1, 1) = 0.5*(fz_z(i,1,1)+fz_z(i,ny-1,1))
    fy_z(i,ny, nz) = fy_z(i,1, nz)
    fy_z(i,ny, 1) = fy_z(i,1, 1)

    ! промежуточные разделы слоёв:
    do k = 2, nz-1
      fy_z(i,:,k) = fy_z(i,:,1)-(fy_z(i,:,1)-fy_z(i,:,nz))/(nz-1.)*(k-1.)
    end do

  enddo
  !HY = fy_z(:,:, 1)
  fy_b = fy_z(:,:, nz)


  ! For WB test
  do k=1,nz
    do j=1,ny-1
      do i=1,nx-1
        fz_z(i,j,k)=0.25*( fx_z(i,j,k)+fx_z(i+1,j,k)+fy_z(i,j,k)+fy_z(i,j+1,k))
      enddo
    enddo
  enddo

  do j=1,ny-1
    do i=1,nx-1
      c_zt(i,j)=fz_z(i,j,1)
      c_b(i,j)=fz_z(i,j,nz)
    enddo
  enddo

  ! потоковые скорости по оси X
  do j=1,ny-1
    do k=1,nz-1
      do i=2,nx-1
        fx_du(i,j,k)=0.5*(c_du(i,j,k)+c_du(i-1,j,k))
        fx_dv(i,j,k)=0.5*(c_dv(i,j,k)+c_dv(i-1,j,k))
      enddo
      fx_du(1,j,k)=0.5*(c_du(1,j,k)+c_du(nx-1,j,k))
      fx_du(nx,j,k)=fx_du(1,j,k)
      fx_dv(1,j,k)=0.5*(c_dv(1,j,k)+c_dv(nx-1,j,k))
      fx_dv(nx,j,k)=fx_dv(1,j,k)
    enddo
  enddo

  ! потоковые скорости по оси Y
  do i=1,nx-1
    do k=1,nz-1
      do j=2,ny-1
        fy_du(i,j,k)=0.5*(c_du(i,j,k)+c_du(i,j-1,k))
        fy_dv(i,j,k)=0.5*(c_dv(i,j,k)+c_dv(i,j-1,k))
      enddo
      fy_du(i,1,k)=0.5*(c_du(i,1,k)+c_du(i,ny-1,k))
      fy_du(i,ny,k)=fy_du(i,1,k)
      fy_dv(i,1,k)=0.5*(c_dv(i,1,k)+c_dv(i,ny-1,k))
      fy_dv(i,ny,k)=fy_dv(i,1,k)
    enddo
  enddo

  ! For WB test
  do k=1,nz-1
    do j=1,ny-1
      do i=1,nx-1
        c_du(i,j,k)=0.25*( fx_du(i,j,k)+fx_du(i+1,j,k)+fy_du(i,j,k)+fy_du(i,j+1,k))
        c_dv(i,j,k)=0.25*(fx_dv(i,j,k)+fx_dv(i+1,j,k)+fy_dv(i,j,k)+fy_dv(i,j+1,k))
      enddo
    enddo
  enddo

  ! на гранях: потоковые переменные мелкой воды:
  do i=1,nxx
    do j=1,nxy
      il = max(i-1,1)
      ir = min(i,ncx)
      fx_h0(i,j) = (c_h0(il,j) + c_h0(ir,j)) / 2.
      fx_u0(i,j) = (c_u0(il,j) + c_u0(ir,j)) / 2.
      fx_v0(i,j) = (c_v0(il,j) + c_v0(ir,j)) / 2.
      fx_rho0(i,j) = (c_rho0(il,j) + c_rho0(ir,j)) / 2.
    end do
  end do
  do i=1,nyx
    do j=1,nyy
      jl = max(j-1,1)
      jr = min(j,ncy)
      fy_h0(i,j) = (c_h0(i,jl) + c_h0(i,jr)) / 2.
      fy_u0(i,j) = (c_u0(i,jl) + c_u0(i,jr)) / 2.
      fy_v0(i,j) = (c_v0(i,jl) + c_v0(i,jr)) / 2.
      fy_rho0(i,j) = (c_rho0(i,jl) + c_rho0(i,jr)) / 2.
    end do
  end do

#endif

end
