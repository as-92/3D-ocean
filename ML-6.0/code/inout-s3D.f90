subroutine Gout
  call TPoutPlt
end subroutine Gout

!--------------------------------------------------------------------------------------------

subroutine TPoutPlt
  use variables
  implicit none

  include 'tecio.h'

  character(1) :: c0
  character(len=200) :: fname, buf
  real(8) :: zt(nx), zi, bi, hi,z1, z2, ti
  real(4) :: data(nx * ny * nz)
  integer(4) :: ic, jc, js, is, iz, isr, isl, jsb, jsf, rc, cnt, ndata
  integer(4) :: i1, i2, i3, i4, i5, i6, i, j, k
  ! переменные:                           x  y  z  u  v  w |u| h rho T fz_z z0 rot id
  integer(4) :: ValueLocation(14) =    (/ 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0 /)         ! 0: point, 1: cell
  integer(4) :: PassiveVarList(14) =   (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0 /)
  integer(4) :: ShareVarFromZone(14) = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0 /)
  ! integer(4) :: ValueLocation(4) = (/ 1, 1, 1, 0 /)         ! 0: point, 1: cell
  ! integer(4) :: PassiveVarList(4) = (/ 0, 0, 0, 0 /)
  ! integer(4) :: ShareVarFromZone(4) = (/ 0, 0, 0, 0 /)
  integer(4) :: ValueLocationA(3) = (/ 1, 1, 1 /)         ! 0: point, 1: cell
  integer(4) :: PassiveVarListA(3) = (/ 0, 0, 0 /)
  integer(4) :: ShareVarFromZoneA(3) = (/ 0, 0, 0 /)
  ! для вычисления вертикальных скоростей:
  real(8) :: f(nx, ny, nz), w(nx-1, ny-1, nz-1), h_x, h_y, fB, fT, zci, zcni, zl, zr, zb, zf, hl, hr, hf, uci, ucni, val
  real(8) :: umx(nxx,nxy),umy(nyx,nyy)                                   ! средние по столбцу потоковые скорости
  real(8) :: Znode(nx, ny, nz)
  real(8) :: hb, dvdx, dudy, rot(ncx, ncy, ncz)

  c0 = char(0)                                                          ! терминальный 0 для c-strings
  ! высоты в узлах сетки
  do iz = nz, 1, -1
    do js = 2, ny-1
      do is = 2, nx-1
        Znode(is,js,iz) = 0.25 * (fx_z(is-1,js,iz)+fx_z(is,js,iz)+fy_z(is,js-1,iz)+fy_z(is,js,iz))
      end do
    end do
    js = 1
    do is = 2, nx-1
      Znode(is,js,iz) = (fx_z(is,js,iz)+fy_z(is-1,js,iz)+fy_z(is,js,iz)) / 3.0
    end do
    js = ny
    do is = 2, nx-1
      Znode(is,js,iz) = (fx_z(is,js-1,iz)+fy_z(is-1,js,iz)+fy_z(is,js,iz)) / 3.0
    end do
    is = 1
    do js = 2, ny-1
      Znode(is,js,iz) = (fx_z(is,js-1,iz)+fx_z(is,js,iz)+fy_z(is,js,iz)) / 3.0
    end do
    is = nx
    do js = 2, ny-1
      Znode(is,js,iz) = (fx_z(is,js-1,iz)+fx_z(is,js,iz)+fy_z(is-1,js,iz)) / 3.0
    end do
    Znode(1,1,iz) = 0.5 * (fx_z(1,1,iz)+fy_z(1,1,iz))
    Znode(1,ny,iz) = 0.5 * (fx_z(1,ny-1,iz)+fy_z(1,ny,iz))
    Znode(nx,1,iz) = 0.5 * (fx_z(nx,1,iz)+fy_z(nx-1,1,iz))
    Znode(nx,ny,iz) = 0.5 * (fx_z(nx,ny-1,iz)+fy_z(nx-1,ny,iz))
  end do

  ! средние скорости на гранях оси X:
  umx = 0.
  do iz = 1, nz-1
    umx = umx + fx_du(:, :, iz) * (fx_z(:, :, iz) - fx_z(:, :, iz+1))
  end do
  do i=1, nx
    do j=1, ny-1
      if(fx_type(i,j)<0) then             ! удаленая грань
        umx(i,j) = 0
      else
        z1 = fx_z(i, j, 1)
        z2 = fx_z(i, j, nz)
        hi = z1 - z2
        umx(i,j) = umx(i,j) / hi
      end if
    end do
  end do

  ! средние скорости на гранях оси Y:
  umy = 0.
  do iz = 1, nz-1
    umy = umy + fy_du(:, :, iz) * (fy_z(:, :, iz) - fy_z(:, :, iz+1))
  end do
  do i=1, nx-1
    do j=1,ny
      if(fy_type(i,j)<0) then             ! удаленая грань
        umy(i,j) = 0
      else
        umy(i,j) = umy(i,j) / (fy_z(i, j, 1) - fy_z(i, j, nz))
      end if
    end do
  end do

  ! расчёт вертикальных компонент скоростей:
#if 0
  if(iout>0) then
    do jc=1, ncy
      jsb = jc
      jsf = jc + 1
      h_y = y(jsf) - y(jsb)
      do ic=1, ncx
        if(c_type(ic,jc)<=CELL_DELETED) cycle
        isl = ic                                                            ! индекс левой грани
        isr = ic + 1                                                        ! индекс правой грани
        h_x = x(isr) - x(isl)
        f(ic, jc, 1) = 0.
        do iz=1, nz-1
          hl = fx_z(isl, jc, iz) - fx_z(isl, jc, iz+1)
          hr = fx_z(isr, jc, iz) - fx_z(isr, jc, iz+1)
          hb = fy_z(ic, jsb, iz) - fy_z(ic, jsb, iz+1)
          hf = fy_z(ic, jsf, iz) - fy_z(ic, jsf, iz+1)
          f(ic, jc, iz+1) = f(ic, jc, iz) + &
            (hr * (fx_du(isr, jc, iz) - umx(isr, jc)) - hl * (fx_du(isl, jc, iz) - umx(isl, jc))) / h_x + &
            (hf * (fy_du(ic, jsf, iz) - umx(ic, jsf)) - hb * (fy_du(ic, jsb, iz) - umx(ic, jsb))) / h_y
        end do
        do iz=1, nz-1
          ! вертикальная скорость:
          zl = (fx_z(isl, jc, iz) + fx_z(isl, jc, iz+1)) / 2.
          zr = (fx_z(isr, jc, iz) + fx_z(isr, jc, iz+1)) / 2.
          zb = (fy_z(ic, jsb, iz) + fy_z(ic, jsb, iz+1)) / 2.
          zf = (fy_z(ic, jsf, iz) + fy_z(ic, jsf, iz+1)) / 2.
          zci = (cs_z(ic, jc, iz) + cs_z(ic, jc, iz+1)) / 2.
          zcni = (fz_z(ic, jc, iz) + fz_z(ic, jc, iz+1)) / 2.
          uci = cs_du(ic, jc, iz)
          ucni = c_du(ic, jc, iz)
          fT = f(ic, jc, iz)
          fB = f(ic, jc, iz+1)

          w(ic, jc, iz) = (zcni - zci) / (dt/2.) + &
                (uci + ucni) / 2. * (zr - zl) / h_x + &
                (uci + ucni) / 2. * (zf - zb) / h_y + &!!! перепроверить!!!!
                (fT + fB) / 2.
        end do
      end do
    end do
  else
    w = 0.
  end if
#endif

  ! ротор скорости:
  do jc=1, ncy
    do ic=1, ncx
      if(c_type(ic,jc)<=CELL_DELETED) cycle
      ! по граням:
      !rot(ic, jc, :) = (vx(ic+1, jc, :) - vx(ic, jc, :)) / dx - (fy_u(ic, jc+1, :) - fy_u(ic, jc, :)) / dy
      ! по ячейкам:
      do iz=1, nz-1
#if 0
        if(ic==1) then
          dvdx = (c_dv(2, jc, iz) - c_dv(1, jc, iz)) / c_dx(1)
        else if(ic==nx-1) then
          dvdx = (c_dv(ic, jc, iz) - c_dv(ic-1, jc, iz)) / c_dx(nx-1)
        else
          dvdx = (c_dv(ic+1, jc, iz) - c_dv(ic-1, jc, iz)) / (c_dx(ic-1)+c_dx(ic))
        end if
        if(jc==1) then
          dudy = (c_du(ic, 2, iz) - c_du(ic, 1, iz)) / c_dy(1)
        else if(jc==ny-1) then
          dudy = (c_du(ic, jc, iz) - c_du(ic, jc-1, iz)) / c_dy(ny-1)
        else
          dudy = (c_du(ic, jc+1, iz) - c_du(ic, jc-1, iz)) / (c_dy(jc-1)+c_dy(jc))
        end if
        rot(ic, jc, iz) = dvdx - dudy
#else
        rot(ic, jc, iz) = 0.
#endif
      end do
    end do
  end do

  ! конструируем имя файла:
  write(fname, "(a,i6.6,a)") "./out/c", iout, ".plt"//c0

  ! удаляем существующие файлы:
  if(iout==0) then
    call execute_command_line("mkdir -p ./out/")
    call execute_command_line("rm -f ./out/*.*")
    call execute_command_line("rm -f ./out/*")
  end if

  write(buf,"(a,i0,a)") "TITLE = ""TEST_NO: ", test_numb, """"
  rc = tecini142( &
          buf//c0, &                                                    ! заголовок
          'x y z u v w absVel h rho T fz_z z0 rot id'//c0, &                         ! переменные
          fname//c0, &                                                  ! имя файла
          './out'//c0, &                                                ! scratch dir
          0, &                                                          ! формат файла: 0=.plt 1=.szplt
          0, &                                                          ! тип файла: 0=full, 1=grid, 2=data
          0, &                                                          ! debug
          1 &                                                           ! 0=real(4), 1=real(8)
  )

  ! общие параметры:
  if(iout==0) then
    write(buf, "(i0,a)") test_numb, c0;   rc = tecauxstr142("task"//c0, buf)
    write(buf, "(f4.2,a)") cfl, c0;       rc = tecauxstr142("cfl"//c0, buf)
    write(buf, "(i0,a)") nx, c0;          rc = tecauxstr142("nx"//c0, buf)
    write(buf, "(i0,a)") ny, c0;          rc = tecauxstr142("ny"//c0, buf)
    write(buf, "(i0,a)") nz, c0;          rc = tecauxstr142("nz"//c0, buf)
    write(buf, "(f4.2,a)") pan, c0;       rc = tecauxstr142("pan"//c0, buf)
    write(buf, "(i0,a)") istep, c0;       rc = tecauxstr142("nstep"//c0, buf)
    write(buf, "(1p,g0,a)") Coriolis, c0;  rc = tecauxstr142("coriolis"//c0, buf)
    write(buf, "(1p,e8.1,a)") sig_0 * CFD_T / CFD_L**2, c0;  rc = tecauxstr142("sig"//c0, buf)
  end if

  !write(buf, "(i0,a)") istep, c0;         rc = tecauxstr142("stepmax"//c0, buf)    ! наибольший достигнутый шаг
  !write(buf, "(1p,g0.4,a)") ttime, c0;    rc = tecauxstr142("tmax"//c0, buf)       ! наибольшее достигнутое время

  ! время:
  ti = GetCTime(ttime)
  !if(test_numb==1) ti = ti / 3600. / 24.              ! сек -> дни

  write(buf,"(i0, a,i0)") iout, ". step ", istep
  rc = teczne142( &
          trim(buf)//c0, &      ! ZoneTitle, &
          0, &                  ! ZoneType, &
          nx, &                 ! IMxOrNumPts, &
          ny, &                 ! JMxOrNumElements, &
          nz, &                 ! KMxOrNumFaces, &
          0, &                  ! ICellMax, &
          0, &                  ! JCellMax, &
          0, &                  ! KCellMax, &
          ti, &                 ! SolutionTime, &
          1, &                  ! StrandID, &
          0, &                  ! ParentZone, &
          1, &                  ! IsBlock, &
          0, &                  ! NumFaceConnections, &
          0, &                  ! FaceNeighborMode, &
          0, &                  ! TotalNumFaceNodes, &
          0, &                  ! NumConnectedBoundaryFaces, &
          0, &                  ! TotalNumBoundaryConnections, &
          PassiveVarList, &               ! PassiveVarList, &
          ValueLocation, &      ! ValueLocation, &
          ShareVarFromZone, &               ! ShareVarFromZone, &
          0 &                   ! ShareConnectivityFromZone &
  )

  ! общие параметры:
  write(buf, "(i0,a)") istep, c0;       rc = teczauxstr142("step"//c0, buf)
  write(buf, "(1p,g0.4,a)") ti, c0;     rc = teczauxstr142("t"//c0, buf)
  ti = ti / 3600. / 24.              ! сек -> дни
  write(buf, "(1p,g0.2,a)") ti, c0;     rc = teczauxstr142("tday"//c0, buf)

  ndata = nx * ny * nz

  ! координата x:
  cnt = 1
  do iz = nz, 1, -1
    do js = 1, ny
      do is = 1, nx
        data(cnt) = real(GetCLength(x(is)), 4)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! координата y:
  cnt = 1
  do iz = nz, 1, -1
    do js = 1, ny
      do is = 1, nx
        data(cnt) = real(GetCLength(y(js)), 4)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! координата Z:
  cnt = 1
  do iz = nz, 1, -1
    do js = 1, ny
      do is = 1, nx
        data(cnt) = real(GetCDepth(Znode(is, js, iz)), 4)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! данные в слой-ячейках:
  ndata = ncx * ncy * ncz

  ! скорость в слое u:
  cnt = 1
  data = 0.
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          !data(cnt) = GetCVelocity(c_u(ic, jc, iz))
          data(cnt) = c_u0(ic,jc)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! скорость в слое v:
  cnt = 1
  data = 0.
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          !data(cnt) = GetCVelocity(c_v(ic, jc, iz))
          data(cnt) = c_v0(ic,jc)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! скорость в слое w:
  cnt = 1
  data = 0.
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          !data(cnt) = w(ic, jc, iz) / CH * CT
          data(cnt)=0.
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! скорость в слое absVel:
  cnt = 1
  data = 0.
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          !data(cnt) = GetCVelocity(SQRT(c_u(ic, jc, iz)**2+c_v(ic, jc, iz)**2))
          data(cnt) = GetCVelocity(SQRT(c_u0(ic, jc)**2+c_v0(ic, jc)**2))
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! высота слоя h:
  cnt = 1
  data = 0.
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          data(cnt) = GetCDepth(fz_z(ic, jc, iz) - fz_z(ic, jc, iz+1))
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! плотность:
  cnt = 1
  data = 0.
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          data(cnt) = GetCDensity(rho0 + c_drho(ic, jc, iz))
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! температура:
  cnt = 1
  data = 0.
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          data(cnt) = c_dteta(ic,jc,iz)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! z
  cnt = 1
  data = 0.
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          data(cnt) = GetCDepth(fz_z(ic, jc, iz) + fz_z(ic, jc, iz+1)) / 2.
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! z0
  cnt = 1
  data = 0.
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          data(cnt) = GetCDepth(fz_z(ic, jc, 1))
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! rot
  cnt = 1
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          data(cnt) = rot(ic, jc, iz)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! id:
  cnt = 1
  do iz = ncz, 1, -1
    do jc = 1, ncy
      do ic = 1, ncx
        if(c_type(ic, jc)>CELL_DELETED) &
          data(cnt) = ic*1000 + jc + iz * 0.01
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  rc = tecend142()

  !...........................................................................................................

  if(test_numb==4) then     ! точное решение

    ! ! удаляем существующие файлы:
    ! if(iout==0) then
    !   call execute_command_line("rm -f ./out2/*.*")
    !   call execute_command_line("rm -f ./out2/*")
    ! end if

    ! ! конструируем имя файла:
    ! write(fname, "(a,i6.6,a)") "./out2/a", iout, ".plt"//c0

    ! write(buf,"(a,i0,a)") "TITLE = ""TEST_NO: ", test_numb, """"
    ! rc = tecini142( &
    !         buf//c0, &                                                    ! заголовок
    !         'x h ha'//c0, &                                               ! переменные
    !         fname//c0, &                                  ! имя файла
    !         './out'//c0, &                                                ! scratch dir
    !         0, &                                                          ! формат файла: 0=.plt 1=.szplt
    !         0, &                                                          ! тип файла: 0=full, 1=grid, 2=data
    !         0, &                                                          ! debug
    !         1 &                                                           ! 0=real(4), 1=real(8)
    ! )
    ! rc = teczne142( &
    !         "analytics"//c0, &    ! ZoneTitle, &
    !         0, &                  ! ZoneType, &
    !         nxx, &             ! IMxOrNumPts, &
    !         1, &                  ! JMxOrNumElements, &
    !         1, &                  ! KMxOrNumFaces, &
    !         0, &                  ! ICellMax, &
    !         0, &                  ! JCellMax, &
    !         0, &                  ! KCellMax, &
    !         ttime, &              ! SolutionTime, &
    !         2, &                  ! StrandID, &
    !         0, &                  ! ParentZone, &
    !         1, &                  ! IsBlock, &
    !         0, &                  ! NumFaceConnections, &
    !         0, &                  ! FaceNeighborMode, &
    !         0, &                  ! TotalNumFaceNodes, &
    !         0, &                  ! NumConnectedBoundaryFaces, &
    !         0, &                  ! TotalNumBoundaryConnections, &
    !         PassiveVarListA, &    ! PassiveVarList, &
    !         ValueLocationA, &     ! ValueLocation, &
    !         ShareVarFromZoneA, &  ! ShareVarFromZone, &
    !         0 &                   ! ShareConnectivityFromZone &
    ! )

    ! ndata = nxx

    ! ! координата x:
    ! cnt = 1
    ! do is = 1, nxx
    !   data(cnt) = real(x(is), 4)
    !   cnt = cnt + 1
    ! end do
    ! rc = tecdat142(ndata, data, 0)

    ! ! координата ZX(1):
    ! cnt = 1
    ! do is = 1, nxx
    !   data(cnt) = real(ZX(is,1), 4)                    ! z-координата
    !   cnt = cnt + 1
    ! end do
    ! rc = tecdat142(ndata, data, 0)

    ! ! точное решение по высоте:
    ! cnt = 1
    ! do is = 1, nxx
    !   bi = task4_b(x(is))                            ! дно
    !   hi = task4_h(x(is))                            ! глубина
    !   data(cnt) = real(bi + hi, 4)                   ! z-координата
    !   cnt = cnt + 1
    ! end do
    ! rc = tecdat142(ndata, data, 0)

    ! rc = tecend142()

  end if

  iout = iout + 1

  return
end subroutine TPoutPlt

!--------------------------------------------------------------------------------------------------------------------

! В отличие от TPoutPlt делит каждую ячейку по средним линиям на 4 прямоугольника. Таким образом точно выводится и
! положение сетки в центре ячейки (консервативные значения h), и в центрах граней (потоковые значения).
! Отличия:
! - Вместо координат одной ячейки выводятся координаты верщин 4 ячеек
! - Вместо ячеечно-центрированных значений, выводятся узловые значения (консервативные, два потоковых и усредненные по 4 граням потоковые)
subroutine TPoutPlt2
  use variables
  implicit none

  interface
    function GetTPoutPlt2Data(i,j,k,vc,vfx,vfy)
      use variables
      implicit none
      integer:: i,j,k
      real(R8):: vc(:,:,:),vfx(:,:,:),vfy(:,:,:)
      real(8):: GetTPoutPlt2Data
    end function GetTPoutPlt2Data
  end interface

  include 'tecio.h'

  character(1) :: c0
  character(len=200) :: fname, buf
  real(8) :: zt(nx), xi, yi, zi, bi, hi, z1, z2, ti
  real(4),pointer:: data(:)
  integer(4) :: ic, jc, js, is, iz, isr, isl, jsb, jsf, rc, cnt, ndata
  integer(4) :: i1, i2, i3, i4, i5, i6, i, j, k, kk
  ! переменные:                           x  y  z  u  v  w  h rho T z0 rot id
  integer(4) :: ValueLocation(12) =    (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0 /)         ! 0: cell, 1: node
  integer(4) :: PassiveVarList(12) =   (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  integer(4) :: ShareVarFromZone(12) = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  !integer(4) :: ValueLocationA(3) = (/ 1, 1, 1 /)         ! 0: point, 1: cell
  !integer(4) :: PassiveVarListA(3) = (/ 0, 0, 0 /)
  !integer(4) :: ShareVarFromZoneA(3) = (/ 0, 0, 0 /)
  ! для вычисления вертикальных скоростей:
  real(8) :: f(nx, ny, nz), w(nx-1, ny-1, nz-1), h_x, h_y, fB, fT, zci, zcni, zl, zr, zb, zf, hl, hr, hf, uci, ucni
  real(8) :: umx(nxx,nxy),umy(nyx,nyy)                                   ! средние по столбцу потоковые скорости
  real(8) :: Znode(nx, ny, nz)
  real(8) :: hb, dvdx, dudy, val
  real(8):: c_rot(nx-1, ny-1, nz), fx_rot(nx, ny-1, nz), fy_rot(nx-1, ny, nz)

  c0 = char(0)                                                          ! терминальный 0 для c-strings
  ! высоты в узлах сетки
  do iz = nz, 1, -1
    do js = 2, nyy-1
      do is = 2, nxx-1
        Znode(is,js,iz) = 0.25 * (fx_z(is-1,js,iz)+fx_z(is,js,iz)+fy_z(is,js-1,iz)+fy_z(is,js,iz))
      end do
    end do
    js = 1
    do is = 2, nxx-1
      Znode(is,js,iz) = (fx_z(is,js,iz)+fy_z(is-1,js,iz)+fy_z(is,js,iz)) / 3.0
    end do
    js = nyy
    do is = 2, nxx-1
      Znode(is,js,iz) = (fx_z(is,js-1,iz)+fy_z(is-1,js,iz)+fy_z(is,js,iz)) / 3.0
    end do
    is = 1
    do js = 2, nyy-1
      Znode(is,js,iz) = (fx_z(is,js-1,iz)+fx_z(is,js,iz)+fy_z(is,js,iz)) / 3.0
    end do
    is = nxx
    do js = 2, nyy-1
      Znode(is,js,iz) = (fx_z(is,js-1,iz)+fx_z(is,js,iz)+fy_z(is-1,js,iz)) / 3.0
    end do
    Znode(1,1,iz) = 0.5 * (fx_z(1,1,iz)+fy_z(1,1,iz))
    Znode(1,nyy,iz) = 0.5 * (fx_z(1,nyy-1,iz)+fy_z(1,nyy,iz))
    Znode(nxx,1,iz) = 0.5 * (fx_z(nxx,1,iz)+fy_z(nxx-1,1,iz))
    Znode(nxx,nyy,iz) = 0.5 * (fx_z(nxx,nyy-1,iz)+fy_z(nxx-1,nyy,iz))
  end do

  ! средние скорости на гранях оси X:
  umx = 0.
  do iz = 1, nz-1
    umx = umx + fx_du(:, :, iz) * (fx_z(:, :, iz) - fx_z(:, :, iz+1))
  end do
  do i=1, nx
    do j=1, ny-1
      z1 = fx_z(i, j, 1)
      z2 = fx_z(i, j, nz)
      hi = z1 - z2
      umx(i,j) = umx(i,j) / hi
    end do
  end do
  ! средние скорости на гранях оси Y:
  umy = 0.
  do iz = 1, nz-1
    umy = umy + fy_du(:, :, iz) * (fy_z(:, :, iz) - fy_z(:, :, iz+1))
  end do
  do i=1, nx-1
    do j=1,ny
      umy(i,j) = umy(i,j) / (fy_z(i, j, 1) - fy_z(i, j, nz))
    end do
  end do

  ! расчёт вертикальных компонент скоростей:
  if(iout>0) then
    do jc=1, ncy
      jsb = jc
      jsf = jc + 1
      h_y = y(jsf) - y(jsb)
      do ic=1, ncx
        isl = ic                                                            ! индекс левой грани
        isr = ic + 1                                                        ! индекс правой грани
        h_x = x(isr) - x(isl)
        f(ic, jc, 1) = 0.
        do iz=1, nz-1
          hl = fx_z(isl, jc, iz) - fx_z(isl, jc, iz+1)
          hr = fx_z(isr, jc, iz) - fx_z(isr, jc, iz+1)
          hb = fy_z(ic, jsb, iz) - fy_z(ic, jsb, iz+1)
          hf = fy_z(ic, jsf, iz) - fy_z(ic, jsf, iz+1)
          f(ic, jc, iz+1) = f(ic, jc, iz) + &
            (hr * (fx_du(isr, jc, iz) - umx(isr, jc)) - hl * (fx_du(isl, jc, iz) - umx(isl, jc))) / h_x + &
            (hf * (fy_du(ic, jsf, iz) - umx(ic, jsf)) - hb * (fy_du(ic, jsb, iz) - umx(ic, jsb))) / h_y
        end do
        do iz=1, nz-1
          ! вертикальная скорость:
          zl = (fx_z(isl, jc, iz) + fx_z(isl, jc, iz+1)) / 2.
          zr = (fx_z(isr, jc, iz) + fx_z(isr, jc, iz+1)) / 2.
          zb = (fy_z(ic, jsb, iz) + fy_z(ic, jsb, iz+1)) / 2.
          zf = (fy_z(ic, jsf, iz) + fy_z(ic, jsf, iz+1)) / 2.
          zci = (fzs_z(ic, jc, iz) + fzs_z(ic, jc, iz+1)) / 2.
          zcni = (fz_z(ic, jc, iz) + fz_z(ic, jc, iz+1)) / 2.
          uci = cs_du(ic, jc, iz)
          ucni = c_du(ic, jc, iz)
          fT = f(ic, jc, iz)
          fB = f(ic, jc, iz+1)

          w(ic, jc, iz) = (zcni - zci) / (dt/2.) + &
                (uci + ucni) / 2. * (zr - zl) / h_x + &
                (uci + ucni) / 2. * (zf - zb) / h_y + &!!! перепроверить!!!!
                (fT + fB) / 2.
        end do
      end do
    end do
  else
    w = 0.
  end if

  ! ротор скорости в ячейках:
  do jc=1, ny-1
    do ic=1, nx-1
      ! по граням:
      !rot(ic, jc, :) = (vx(ic+1, jc, :) - vx(ic, jc, :)) / dx - (fy_u(ic, jc+1, :) - fy_u(ic, jc, :)) / dy
      ! по ячейкам:
      do iz=1, nz-1
        if(ic==1) then
          dvdx = (c_dv(2, jc, iz) - c_dv(1, jc, iz)) / c_dx(1)
        else if(ic==nx-1) then
          dvdx = (c_dv(ic, jc, iz) - c_dv(ic-1, jc, iz)) / c_dx(nx-1)
        else
          dvdx = (c_dv(ic+1, jc, iz) - c_dv(ic-1, jc, iz)) / (c_dx(ic-1)+c_dx(ic))
        end if
        if(jc==1) then
          dudy = (c_du(ic, 2, iz) - c_du(ic, 1, iz)) / c_dy(1)
        else if(jc==ny-1) then
          dudy = (c_du(ic, jc, iz) - c_du(ic, jc-1, iz)) / c_dy(ny-1)
        else
          dudy = (c_du(ic, jc+1, iz) - c_du(ic, jc-1, iz)) / (c_dy(jc-1)+c_dy(jc))
        end if
        c_rot(ic, jc, iz) = dvdx - dudy
      end do
    end do
  end do

  ! ротор на X-гранях:
  fx_rot(2:nx-1,:,:) = (c_rot(1:nx-2,:,:) + c_rot(2:nx-1,:,:)) / 2.     ! внутренние грани
  if(bcTypeX==1) then   ! периодичесткие ГУ
    fx_rot(1,:,:) = (c_rot(1,:,:) + c_rot(nx-1,:,:)) / 2.
    fx_rot(nx,:,:) = fx_rot(1,:,:)
  else
    fx_rot(1,:,:) = c_rot(1,:,:)
    fx_rot(nx,:,:) = c_rot(nx-1,:,:)
  end if

  ! ротор на Y-гранях:
  fy_rot(:,2:ny-1,:) = (c_rot(:,1:ny-2,:) + c_rot(:,2:ny-1,:)) / 2.     ! внутренние грани
  if(bcTypeY==1) then   ! периодичесткие ГУ
    fy_rot(:,1,:) = (c_rot(:,1,:) + c_rot(:,ny-1,:)) / 2.
    fy_rot(:,ny,:) = fy_rot(:,1,:)
  else
    fy_rot(:,1,:) = c_rot(:,1,:)
    fy_rot(:,ny,:) = c_rot(:,ny-1,:)
  end if

  ! конструируем имя файла:
  write(fname, "(a,i6.6,a)") "./out/c", iout, ".plt"//c0

  ! удаляем существующие файлы:
  if(iout==0) then
    call execute_command_line("mkdir -p ./out/")
    call execute_command_line("rm -f ./out/*.*")
    call execute_command_line("rm -f ./out/*")
  end if

  write(buf,"(a,i0,a)") "TITLE = ""TEST_NO: ", test_numb, """"
  rc = tecini142( &
          buf//c0, &                                                    ! заголовок
          'x y z u v w h rho T z0 rot id'//c0, &                        ! переменные
          fname//c0, &                                                  ! имя файла
          './out'//c0, &                                                ! scratch dir
          0, &                                                          ! формат файла: 0=.plt 1=.szplt
          0, &                                                          ! тип файла: 0=full, 1=grid, 2=data
          0, &                                                          ! debug
          1 &                                                           ! 0=real(4), 1=real(8)
  )

  ! общие параметры:
  if(iout==0) then
    write(buf, "(i0,a)") test_numb, c0;   rc = tecauxstr142("task"//c0, buf)
    write(buf, "(f4.2,a)") cfl, c0;       rc = tecauxstr142("cfl"//c0, buf)
    write(buf, "(i0,a)") nx, c0;          rc = tecauxstr142("nx"//c0, buf)
    write(buf, "(i0,a)") ny, c0;          rc = tecauxstr142("ny"//c0, buf)
    write(buf, "(i0,a)") nz, c0;          rc = tecauxstr142("nz"//c0, buf)
    write(buf, "(f4.2,a)") pan, c0;       rc = tecauxstr142("pan"//c0, buf)
    write(buf, "(i0,a)") istep, c0;       rc = tecauxstr142("nstep"//c0, buf)
    write(buf, "(1p,g0,a)") Coriolis, c0;  rc = tecauxstr142("coriolis"//c0, buf)
    write(buf, "(1p,e8.1,a)") sig_0 * CFD_T / CFD_L**2, c0;  rc = tecauxstr142("sig"//c0, buf)
  end if

  !write(buf, "(i0,a)") istep, c0;         rc = tecauxstr142("stepmax"//c0, buf)    ! наибольший достигнутый шаг
  !write(buf, "(1p,g0.4,a)") ttime, c0;    rc = tecauxstr142("tmax"//c0, buf)       ! наибольшее достигнутое время

  ! время:
  ti = GetCTime(ttime)
  !if(test_numb==1) ti = ti / 3600. / 24.              ! сек -> дни

  write(buf,"(i0, a,i0)") iout, ". step ", istep
  rc = teczne142( &
          trim(buf)//c0, &      ! ZoneTitle, &
          0, &                  ! ZoneType, &
          2*nxx-1, &             ! IMxOrNumPts, &
          2*nyy-1, &                 ! JMxOrNumElements, &
          nz, &                  ! KMxOrNumFaces, &
          0, &                  ! ICellMax, &
          0, &                  ! JCellMax, &
          0, &                  ! KCellMax, &
          ti, &                 ! SolutionTime, &
          1, &                  ! StrandID, &
          0, &                  ! ParentZone, &
          1, &                  ! IsBlock, &
          0, &                  ! NumFaceConnections, &
          0, &                  ! FaceNeighborMode, &
          0, &                  ! TotalNumFaceNodes, &
          0, &                  ! NumConnectedBoundaryFaces, &
          0, &                  ! TotalNumBoundaryConnections, &
          PassiveVarList, &               ! PassiveVarList, &
          ValueLocation, &      ! ValueLocation, &
          ShareVarFromZone, &               ! ShareVarFromZone, &
          0 &                   ! ShareConnectivityFromZone &
  )

  ! общие параметры:
  write(buf, "(i0,a)") istep, c0;       rc = teczauxstr142("step"//c0, buf)
  write(buf, "(1p,g0.4,a)") ti, c0;     rc = teczauxstr142("t"//c0, buf)
  ti = ti / 3600. / 24.              ! сек -> дни
  write(buf, "(1p,g0.2,a)") ti, c0;     rc = teczauxstr142("tday"//c0, buf)

  ! число узлов сетки:
  ndata = (nxx*2-1) * (nyy*2-1) * nz
  allocate (data(ndata))

  ! координата x:
  cnt = 1
  do iz = nz, 1, -1
    do js = 1, 2*nyy-1
      do is = 1, nxx
        data(cnt) = real(GetCLength(x(is)), 4)
        cnt = cnt + 1
        if(is/=nxx) then
          xi = (x(is)+x(is+1))/2.
          data(cnt) = real(GetCLength(xi), 4)
          cnt = cnt + 1
        end if
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! координата y:
  cnt = 1
  do iz = nz, 1, -1
    do js = 1, 2*nyy-1
      if(mod(js,2)==1) then
        yi = y((js+1)/2)
      else
        yi = (y(js/2) + y(js/2+1)) / 2.
      end if
      do is = 1, 2*nxx-1
        data(cnt) = real(GetCLength(yi), 4)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! координата Z:
  !ndata = nxx * nyy * nz
  cnt = 1
  do iz = nz, 1, -1
    do js = 1, 2*nyy-1          ! счетчик маленьких узлов
      do is = 1, 2*nxx-1
        zi = GetTPoutPlt2Data(is,js,iz,fz_z,fx_z,fy_z)
        data(cnt) = real(GetCDepth(zi), 4)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! данные в слой-ячейках:
  !ndata = 2*ncx * 2*ncy * (nz-1)

  ! скорость в слое u:
  cnt = 1
  do k = nz, 1, -1
    do j = 1, 2*ncy-1
      do i = 1, 2*ncx-1
        kk = k
        if(kk==nz) kk = nz - 1
        val = GetTPoutPlt2Data(i,j,kk,c_du,fx_du,fy_du)
        data(cnt) = GetCVelocity(val)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! скорость в слое v:
  cnt = 1
  do k = nz, 1, -1
    do j = 1, 2*ncy-1
      do i = 1, 2*ncx-1
        kk = k
        if(kk==nz) kk = nz - 1
        val = GetTPoutPlt2Data(i,j,kk,c_dv,fx_dv,fy_dv)
        data(cnt) = GetCVelocity(val)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! скорость в слое w:
  cnt = 1
  do k = nz, 1, -1
    do j = 1, 2*ncy-1
      do i = 1, 2*ncx-1
        data(cnt) = 0.
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! высота слоя h:
  cnt = 1
  do k = nz, 1, -1
    do j = 1, 2*ncy-1
      do i = 1, 2*ncx-1
        kk = k
        if(kk==nz) kk = nz - 1
        val = 0.
        data(cnt) = GetCDepth(val)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! плотность:
  cnt = 1
  do k = nz, 1, -1
    do j = 1, 2*ncy-1
      do i = 1, 2*ncx-1
        kk = k
        if(kk==nz) kk = nz - 1
        val = rho0+GetTPoutPlt2Data(i,j,kk,c_drho,fx_drho,fy_drho)
        data(cnt) = GetCDensity(val)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! температура:
  cnt = 1
  do k = nz, 1, -1
    do j = 1, 2*ncy-1
      do i = 1, 2*ncx-1
        kk = k
        if(kk==nz) kk = nz - 1
        val = GetTPoutPlt2Data(i,j,kk,c_dteta,fx_dteta,fy_dteta)
        data(cnt) = GetCDensity(val)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! z0
  cnt = 1
  do k = nz, 1, -1
    do j = 1, 2*ncy-1
      do i = 1, 2*ncx-1
        val = GetTPoutPlt2Data(i,j,1,c_dv,fx_dv,fy_dv)
        data(cnt) = GetCDepth(val)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! rot
  cnt = 1
  do k = nz, 1, -1
    do j = 1, 2*ncy-1
      do i = 1, 2*ncx-1
        kk = k
        if(kk==nz) kk = nz - 1
        val = GetTPoutPlt2Data(i,j,kk,c_rot,fx_rot,fy_rot)
        data(cnt) = GetCVelocity(val)
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  ! переменные в ячейках:
  ndata = 2*(nxx-1) * 2*(nyy-1) * (nz-1)

  ! id
  cnt = 1
  do k = nz-1, 1, -1
    do j = 1, 2*(ncy-1)
      do i = 1, 2*(ncx-1)
        val = real(((i+1)/2)*1000 + (j+1)/2)
        data(cnt) = val
        cnt = cnt + 1
      end do
    end do
  end do
  rc = tecdat142(ndata, data, 0)

  rc = tecend142()
  deallocate(data)

  iout = iout + 1

  return

end subroutine TPoutPlt2

!-----------------------------------------------------------------------------------------------

function GetTPoutPlt2Data(i,j,k,vc,vfx,vfy)
  use variables
  implicit none
  integer:: i,j,k
  real(R8):: vc(:,:,:),vfx(:,:,:),vfy(:,:,:)
  real(8):: GetTPoutPlt2Data
  integer:: cs,il,ir,jb,jt,i0,j0
  real(8):: v, vl, vr, vt, vb

  if(k==0) then
    debdum = 0
  end if

  cs = 2*mod(j,2) + mod(i,2)
  select case(cs)
    case(0)                                   ! центр большой ячейки
      v = vc(i/2,j/2,k)
    case(1)                                   ! центр X-грани
      v = vfx((i+1)/2,j/2,k)
    case(2)                                   ! центр Y-грани
      v = vfy(i/2,(j+1)/2,k)
    case(3)                                   ! узел большой ячейки

      ! нижняя X-грань:
      jb = (j-1)/2
      if(jb<1) jb = 1                         ! находимся на ymin, нижней грани нет

      ! верхняя X-грань:
      jt = (j-1)/2 + 1
      if(jt>=ny) jt = ny - 1                  ! находимся на ymax, верхней грани нет

      ! левая Y-грань
      il = (i-1)/2
      if(il<1) il = 1

      ! правая Y-грань
      ir = (i-1)/2 + 1
      if(ir>=nx) ir = nx-1

      i0 = (i+1)/2
      j0 = (j+1)/2
      v = ( vfx(i0,jb,k) + vfx(i0,jt,k) + &
            vfy(il,j0,k) + vfy(ir,j0,k) ) / 4.
  end select

  GetTPoutPlt2Data = v
end function GetTPoutPlt2Data

