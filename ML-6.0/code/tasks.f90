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

!============================================================================================================

! Проточный тест с аналитическим решением
subroutine SetupTask_4
  use variables
  implicit none

  integer:: i, ic, iz, j, jc, k
  real(R8):: dx, dy
  real(R8):: tci, bci
  real(R8):: ui, vi, xi, yi, zi, z0

  dl=20.; dw=1.                                                 ! размер области
  nx=65; ny=2; nz=21                                             ! количество узлов сетки

  bcTypeX = BC_IN_Z; bcTypeY = BC_SLIDE

  call allocate

  ! отладка:
  debstartstep = 0; debendstep = 5
  debstartz = 1; debendz = 400

  !debcells(51, 1, 20) = 20
  !debcells(ncx, 1) = 1
  debSidesX(1,1,:) = 1
  !debSidesX(nxx, 1) = 1

  !debCells(20, 12) = 1
  !debCells(12, 20) = 1
  !debSidesX(20, 12) = 1
  !debsidesY(12, 20) = 1

  !debcells(31, 1) = 1
  !debSidesX(31, 1) = 1
  debSidesZ(1,1,:) = 1

#if 1
  nt=20000                                                          ! полное число шагов по времени
  nprint=max(1,nt/200)                                           ! интервал печати
#else
  maxt=100.
  tprint=maxt/200.
#endif

  CFL = 0.3                                                     ! Число Куранта

  pan = 0.0                                                     !- параметр Паниковского

  ! константы:
  g = 9.8d0
  rho0 = 1.
  teta0 = 1.
  sound0 = 1.             ! псевдо-скорость звука:
  talpha = 0.1
  tbeta = 1.

  ! Расчет начальной сетки
  x(1) = 0.
  dx = dl / (nx-1)
  do i=1,ncx; x(i+1) = x(i) + dx; end do                        ! координаты узлов сетки

  y(1) = -dw / 2.
  dy = dw / (ny-1)
  do j=1,ny-1; y(j+1) = y(j) + dy; end do                       ! координаты узлов сетки

  z0=task4_h(x(1))                                              ! фоновый уровень поверхности
  b0=0.                                                         ! фоновый уровень дна

  do jc=1, ncy
    do ic=1, ncx
      xi = (x(ic) + x(ic+1)) / 2.                               ! x-координата ячейки

      bci = task4_b(xi)                                         ! уровень дна в ячейке
      tci = task4_h(xi) + task4_b(xi)                           ! уровень поверхности в ячейке

      ! равномерная сетка по вертикали от поверхности до дна:
      fz_z(ic, jc, 1) = tci
      fz_z(ic, jc, nz) = bci
      do iz = 2, nz-1
        fz_z(ic, jc, iz) = bci + (nz - iz) * (tci - bci)/(nz - 1.)
      end do
    end do
  end do

  ! начальные данные на сетке; слои:
  c_rho = rho0
  c_teta = teta0
  do i=1,ncx; do j=1,ncy; do k=1,ncz
    xi = (x(i) + x(i+1)) / 2.                                   ! x-координата ячейки
    zi = (fz_z(i,j,k) + fz_z(i,j,k+1)) / 2.
    c_u(i,j,k) = task4_u(xi, zi)
    c_v(i,j,k) = 0.
    c_w(i,j,k) = 0.
  end do; end do; end do

  ! начальные данные на сетке; мелкая вода:
  do i=1,ncx; do j=1,ncy;
    ui = 0.
    do k=1,ncz
      ui = ui + c_u(i,j,k) * (fz_z(i,j,k) - fz_z(i,j,k+1))
    end do
    c_h0(i,j) = fz_z(i,j,1) - fz_z(i,j,nz)
    c_u0(i,j) = ui / c_h0(i,j)
    c_v0(i,j) = 0.
    c_rho0 = rho0
  end do; end do

  if(bctypex==BC_IN) then                                             ! ГУ "вход" слева и справа

    allocate(bcData(2))                                              ! данные на левом и правом торце

    ! левая грань:
    xi = x(1)
    bcData(1).h = task4_h(xi)
    bcData(1).u = 0.
    bcData(1).v = 0.
    bcData(1).w = 0.
    bcData(1).rho = rho0
    bcData(1).teta = teta0

    i = 1
    do j=1,ncy
      fx_bc(i,j) = 1
    end do

    ! правая грань:
    bcData(2).h = z0 - b0
    bcData(2).u = 0.
    bcData(2).v = 0.
    bcData(2).w = 0.
    bcData(2).rho = rho0
    bcData(2).teta = teta0

    i = nxx
    do j=1,ncy
      fx_bc(i,j) = 2
    end do
  end if
end

!---------------------------------------------------------------------------------------------------------------------

! Вихрь X-Y
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

!============================================================================================================

! Вихревая пара X-Z (всплытие)
subroutine SetupTask_9
  use variables
  implicit none

  interface
    subroutine CalcVertexXZ(xi, zi, x0, z0, alpha, beta, r0, ui, wi, tetai)
      real(R8):: xi, zi, x0, z0, alpha, beta, r0, ui, wi, tetai
    end subroutine
  end interface

  integer:: i, ic, il, ir, iz, j, jc, jl, jr, k
  real(R8):: alphah, bci, betah, cosa, delbot, deltop, dmu_0, dx, dy
  real(R8):: dzbot, dzbotcoef, dztop, dztopcoef, h_0
  real(R8):: p_0, r0, r_0, rbot, ri, rtop, a, tci, u_0
  real(R8):: ui, vi, zi, hi, xi, yi, z0, xc0, zc0, dcx

  dl=5.; dw=0.1; dh=5.                                          ! размер области
  ncx=50; ncy=1; ncz=50;                                        ! количество ячеек сетки
  bctypex = BC_PERIODIC; bctypey = BC_PERIODIC                  ! тип ГУ на боковых границах

  nx = ncx+1; ny = ncy+1; nz = ncz+1
  call allocate

  ! отладка:
  debstartstep = 0; debendstep = -2
  debstartz = 28; debendz = 30

  debcells(24, 1, 28) = 1
  debcells(24, 1, 29) = 1
  debcells(28, 1, 24) = 1
  debcells(29, 1, 24) = 1
  !debSidesX(24, 1) = 1
  debSidesZ(24, 1, 29) = 1

  !debSidesY(44, 1) = 1
  !debsidesY(44, 2) = 1

  !debcells(23, 1) = 1
  !debSidesX(31, 1) = 1
  !debSidesZ(31, 1) = 1

#if 0
  nt=2000 ! полное число шагов по времени
  nprint=10 ! интервал печати
#else
  maxt=100.
  tprint=maxt/400.
#endif

  CFL = 0.3 ! Число Куранта
  pan=0.0 !- параметр Паниковского

  ! константы:
  g = 1
  rho0 = 1.
  teta0 = 1.
  sound0 = 1. !0.15 * sqrt(g * dh)             ! псевдо-скорость звука

  ! Расчет начальной сетки
  x(1) = -dl / 2.
  dx = dl / (nx-1)
  do i=1,ncx; x(i+1) = x(i) + dx; end do                    ! координаты узлов сетки

  y(1) = -dw / 2.
  dy = dw / (ny-1)
  do j=1,ny-1; y(j+1) = y(j) + dy; end do                   ! координаты узлов сетки

  z0 = 0.                                                   ! фоновый уровень поверхности
  b0 = z0 - dh                                              ! фоновый уровень дна

  ! параметры вихря:
  alphaH = 0.404
  betaH = 0.4
  r_0= 0.25

  xc0 = 0.                                                  ! центр вихревой пары
  dcx = r_0 * 3.                                            ! расстояние между центрами вихрей
  zc0 = b0 + r_0 * 8.                                       ! уровень вихрей по z

  !dcx = 0.
  !zc0 = (b0 + z0) / 2.

  ! равномерная сетка по вертикали от поверхности до дна:
  do iz = 1, nz
    fz_z(:, :, iz) = z0 - (iz - 1.) * dh / (nz - 1.)
  end do

  ! фоновые параметры:
  c_rho = rho0
  c_teta = teta0
  c_u = 0.
  c_v = 0.
  c_w = 0.

  ! начальные данные на сетке:
  do i=1, ncx
    do j=1, ncy
      do k = 1, ncz

        xi = (x(i) + x(i+1)) / 2.                             ! x-координата ячейки
        yi = (y(j) + y(j+1)) / 2.                             ! y-координата ячейки
        zi = (fz_z(i, j, k) + fz_z(i, j, k+1)) / 2.             ! z-кордината слой-ячейки

        ! два вихря:
        call CalcVertexXZ(xi, zi, xc0 - dcx/2., zc0, -alphaH, betaH, r_0, c_u(i,j,k), c_w(i,j,k), c_teta(i,j,k))
        call CalcVertexXZ(xi, zi, xc0 + dcx/2., zc0,  alphaH, betaH, r_0, c_u(i,j,k), c_w(i,j,k), c_teta(i,j,k))

      end do
    end do
  end do

  ! мелкая вода:
  c_h0 = dh
  c_rho0 = rho0
  c_u0 = 0.
  c_v0 = 0.

end

!================================================================================================================

! Бугор на дне и на поверхности:
subroutine SetupTask_10
  use variables
  implicit none

  integer:: i, ic, il, ir, iz, j, jc, jl, jr, k
  real(R8):: alfah, bci, betah, cosa, delbot, deltop, dmu_0, dx, dy
  real(R8):: dzbot, dzbotcoef, dztop, dztopcoef, h_0, hi
  real(R8):: p_0, r0, r_0, rbot, ri, rtop, stepen, tci, u_0, ul, ur
  real(R8):: ui, vi, xbot, ybot, xci, xtop, ytop, yci, z0

#if 1         /* 1D-X */

  dl=20.; dw=1.                                                 ! размер области
  nx=65; ny=2; nz=21                                             ! количество узлов сетки

#elif 0       /* 1D-Y */

  dl=1.; dw=20.                                                 ! 1D-Y: размер области
  nx=2; ny=65; nz=2;                                            ! количество узлов сетки

#else         /* 2D */

  dl=20.; dw=20.                                                ! 2D: размер области
  nx=45; ny=45; nz=2;                                           ! количество узлов сетки

#endif

  bctypex = BC_PERIODIC; bctypey = BC_PERIODIC                  ! тип ГУ на боковых границах - периодика
  bctypex = BC_SLIDE; bctypey = BC_SLIDE                        ! тип ГУ на боковых границах - скольжение
  bctypex = BC_IN; bctypey = BC_SLIDE

  call allocate

  ! отладка:
  debstartstep = 0; debendstep = 5
  debstartz = 1; debendz = 400

  debcells(51, 1, 20) = 20
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

#if 1
  nt=20000                                                          ! полное число шагов по времени
  nprint=max(1,nt/200)                                           ! интервал печати
#else
  maxt=100.
  tprint=maxt/200.
#endif

  CFL = 0.3                                                     ! Число Куранта

  pan = 0.1                                                     !- параметр Паниковского

  g = 1

  ! константы:
  rho0 = 1.
  teta0 = 1.
  sound0 = 1.             ! псевдо-скорость звука:

  ! Расчет начальной сетки
  x(1) = -dl / 2.
  dx = dl / (nx-1)
  do i=1,ncx; x(i+1) = x(i) + dx; end do                        ! координаты узлов сетки

  y(1) = -dw / 2.
  dy = dw / (ny-1)
  do j=1,ny-1; y(j+1) = y(j) + dy; end do                       ! координаты узлов сетки

  z0=1.                                                         ! фоновый уровень поверхности
  b0=0.                                                         ! фоновый уровень дна

  ! для поверхности:
  dztopcoef = 0!.5 !0.05                                              ! коэффициент подъема поверхности (0<cf<1)
  rtop = 0.2 * max(dl,dw)                                      ! ширина возмущения поверхности
  xtop = x(1) + 0.3 * dl                                        ! координата центра  на поверхности
  ytop = y(1) + 0.5 * dw
  dztop = dztopcoef * (z0 - b0)                                 ! вертикальный размер возмущения

  ! возмущение дна:
  dzbotcoef = 0.5                                               ! коэффициент подъема дна (0<cf<1)
  rbot = 0.2 * max(dl,dw)                                      ! ширина возмущения дна
  xbot = x(1) + 0.8 * dl                                        ! координата возмущения
  ybot = y(1) + 0.5 * dw
  dzbot = dzbotcoef * (z0 - b0)                                 ! вертикальная величина колебания дна

  do jc=1, ncy
    yci = (y(jc) + y(jc+1)) / 2.                                ! y-координата ячейки
    do ic=1, ncx
      xci = (x(ic) + x(ic+1)) / 2.                              ! x-координата ячейки

      ! бугор на поверхности:
      deltop = 0.
      r0 = max(dl, dw)
      do i=-1,1; do j=-1,1
        ri = sqrt((xci+dl*i-xtop)**2 + (yci+dw*i-ytop)**2)
        if(ri<r0) r0 = ri
      end do; end do

      if(r0<rtop) then                                          ! если в пределах возмущения
        ri = r0 / rtop                                          ! нормируем на 0..1
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

#if 0
      ! равномерная сетка по вертикали от поверхности до дна:
#else
      ! НЕравномерная сетка по вертикали:
      h_0 = 0.                                                   ! сумма для нормировки
      do k=1,ncz
        if(k==1) then
          hi = 1.                                               ! шаг в слое
        else
          hi = hi * 1.1
        end if
        cfh(k) = hi
        h_0 = h_0 + hi
      end do

      ! нормировка на общую толщину tci-bci:
      do k=1, ncz
        cfh(k) = cfh(k) / h_0
      end do

      fz_z(ic, jc, 1) = tci
      fz_z(ic, jc, nz) = bci

#endif

      ! координаты разделов слоёв:
      fz_z(ic, jc, 1) = tci
      fz_z(ic, jc, nz) = bci
      do k = 2, nz-1
        !fz_z(ic, jc, iz) = tci - (iz - 1.) * (tci - bci)/(nz - 1.)
        fz_z(ic, jc, k) = fz_z(ic, jc, k-1) - (tci - bci) * cfh(k-1)
      end do
    end do
  end do

  ! начальные данные на сетке; мелкая вода:
  u_0 = 0.1                                                            ! скорость потока
  ul = 0.1
  ur = 0.1
  c_h0(:,:) = fz_z(:,:,1) - fz_z(:,:,nz)
  c_rho0 = rho0
  dmu_0 = u_0 * c_h0(1,1)                                             ! импульс в одной из ячеек
  c_u0 = dmu_0 / c_h0                                                 ! скорость с учётом равного импульса
  c_v0 = 0.

  ! начальные данные на сетке; слои:
  c_rho = rho0
  c_teta = teta0
  do i=1,ncx; do j=1,ncy; do k=1,ncz
    c_u(i,j,k) = c_u0(i,j)
    c_v(i,j,k) = c_v0(i,j)
    c_w(i,j,k) = 0.
  end do; end do; end do

  if(bctypex==BC_IN) then                                             ! ГУ "вход" слева и справа

    allocate(bcData(2))                                              ! данные на левом и правом торце

    ! левая грань:
    bcData(1).h = z0 - b0
    bcData(1).u = ul
    bcData(1).v = 0.
    bcData(1).w = 0.
    bcData(1).rho = rho0
    bcData(1).teta = teta0

    i = 1
    do j=1,ncy
      fx_bc(i,j) = 1
    end do

    ! правая грань:
    bcData(2).h = z0 - b0
    bcData(2).u = ur
    bcData(2).v = 0.
    bcData(2).w = 0.
    bcData(2).rho = rho0
    bcData(2).teta = teta0

    i = nxx
    do j=1,ncy
      fx_bc(i,j) = 2
    end do
  end if

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

!===============================================================================================

subroutine CalcVertexXZ(xi, zi, x0, z0, alpha, beta, r0, ui, wi, tetai)
  use variables
  implicit none
  real(R8):: xi, zi, x0, z0, alpha, beta, r0, ui, wi, tetai
  real(R8):: ri, a, delu, delw, delteta

  ri = sqrt((xi-x0)**2 + (zi-z0)**2)                    ! расстояние до центра (xc0, zc0)
  a = beta * (1. - (ri**2 / r0**2))
  delu = alpha / r0 * exp(a) * (zi - z0)
  delw = -alpha / r0 * exp(a) * (xi - x0)
  delteta = -alpha**2 * exp(2 * a) / (4. * sound0**2 * beta)

  ui = ui + delu
  wi = wi + delw
  tetai = tetai + delteta
end subroutine
