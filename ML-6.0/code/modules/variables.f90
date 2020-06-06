module variables
  use, intrinsic:: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
  ! опции:
  logical, parameter:: optLimitMinRhoC = .false.        ! для консервативных переменных: if(rho<eps) rho = eps
  logical, parameter:: optLimitMinHC = .false.          ! для консервативных переменных: if(h<eps) h = eps
  logical, parameter:: optLimitMinRhoF = .false.        ! для потоковых переменных: if(rho<eps) rho = eps
  logical, parameter:: optLimitMinHF = .false.          ! для потоковых переменных: if(h<eps) h = eps
  logical, parameter:: optNoAssert = .false.            ! отмена диагностики и авоста по ассертам

  integer(1), parameter:: BC_DELETED = -1                           ! уничтоженная грань
  integer(1), parameter:: BC_INNER = 0                              ! внутренняя грань
  integer(1), parameter:: BC_PERIODIC = 1                           ! периодический ГУ
  integer(1), parameter:: BC_SLIDE = 2                              ! скольжение
  integer(1), parameter:: BC_STICK = 5                              ! прилипание
  integer(1), parameter:: BC_IN = 6                                 ! вток
  integer(1), parameter:: BC_IN_T = 7                               ! вток, переменный по времени, индивидуальный для каждой грани
  integer(1), parameter:: BC_FIX_U = 8                              ! на входе задана скорость
  integer(1), parameter:: BC_H_VAR_T = 9
  integer(1), parameter:: NBC = 9                                   ! максимальный номер ГУ для грани

  integer(1), parameter:: CELL_DELETED_SPEC = -2                    ! уничтоженная специальная ячейка
  integer(1), parameter:: CELL_DELETED = -1                         ! уничтоженная ячейка
  integer(1), parameter:: CELL_INNER = 0                            ! внутренняя ячейка
  integer(1), parameter:: CELL_INNER_SPEC = 1                       ! внутренняя спец-ячейка

  integer(4), parameter:: T_INVR=1, T_INVQ=2, T_INVU=3, T_INVW=5, T_INVD=4    ! номера инвариантов
  integer(4), parameter:: DIRM=0, DIRP=1                            ! направления

  real(R8):: NAN

  type TBcData                                                      ! данные для расчета ГУ "вток"
    real(R8):: h                                                    ! толщина слоя
    real(R8):: u                                                    ! компоненты скорости
    real(R8):: v
    real(R8):: w
    real(R8):: rho                                                  ! реальная плотность
    real(R8):: teta
    real(R8):: T                                                    ! температура
    integer:: ind                                                   ! ссылка на дополнительные данные
  end type TBcData

  type Vector2
    real(R8):: x, y
  end type

  type Vector
    real(R8):: x, y, z
  end type

  type TopexData
    integer:: i, j                                                  ! индекс ячейки, ближайшей к точке данных
    real(R8):: lat, lon                                             ! географические координаты точки данных
    real(R8):: mAmp(3), mPh(3), sAmp(3), sPh(3)                     ! характеристики прилива
  end type

  integer:: bcTypeTopT,bcTypeTopU, bcTypeX, bcTypeY
  integer:: iout, ioutX, ioutY
  !, nxcells, nycells, nxsides, nysides
  integer:: ncx, ncy                  ! число ячеек по x и y (массивы c_..)
  integer:: nxx, nxy                  ! число X-граней по x и y (массивы fx_..)
  integer:: nyx, nyy                  ! число Y-граней по x и y (массивы fy_..)
  integer:: ncz                       ! число слоёв

  real(R8),pointer::Tsurf(:,:)   ! Температура на поверхности

  real(R8):: CFD_T, CFD_L, CFD_H, CFD_TEMP, CFD_RHO      ! размерные коэффициенты

  real(R8):: cfViscVertT,uMod,cfViscVert,cf3ViscVertT
  real(R8):: cfBottomFriction,uAtm,vAtm
  real(R8):: cfSurfFrictionU,cfSurfFrictionV,tmp,minHeight
  real(R8):: dtvisc                                  ! максимально допустимый шаг - ограничение от горизонтальной вязкости
  real(R8):: rho0                                     ! плотность - константаё
  real(R8):: teta0
  real(R8):: sound0, sound0_2                         ! модельная "скорость звука" и её квадрат
  real(R8):: cs_sound                                 ! "скорость звука" в характеристических уравнениях: = sound0 / sqrt(rho0)
  real(R8):: cfT_rho0, cfT_alpha                      ! параметры в уравнении состояния rho(T)

  real(R8):: Coriolis                                 ! сила Кориолиса

  logical:: isImplExchange                            ! выбор способа обмена между слоями: true - явно-неявный метод АВС (см. sig_ex)
  logical:: isImplExchangeS2                          ! надо ли выполнить обмен на гранях после фазы 2?
  logical:: isImplExchangeS3                          ! надо ли выполнить обмен на гранях после фазы 3?
  logical:: isExchangePh1                             ! делать ли обмен после фазы 1
  logical:: isFixCover                                ! расчет с фиксированной крышкой
  logical:: isConsFilterZ0                            ! фильтрация консервативного уровня поверхности
  logical:: isImplCoriolis                            ! расчет по полунеявной формуле Кориолиса
  real(R8):: sig_ex

  integer(4) nx,ny,nz,nprint,istep,nt,printproc,kface
  integer(4) print_debug,print_density,print_height
  integer(4) print_velocity,print_velocity_by_color
  integer:: taskNum
  !INTEGER(4) isNanHeight,isNanVelocity,isNanDensity,test_numb,mI,mJ,mK,loadStep,loadIout

  real(R8):: dl, dw, dh, xmin, xmax, ymin, ymax
  real(R8):: gxmax, gxmin, gymax, gymin                             ! мин-мах в градусах
  real(R8):: lgradx, lgrady                                         ! длина одного градуса в метрах
  real(R8):: cfl, dt, ttime, g, pan, sig_0, sig, sigma, sigma2Ph

  real(R8):: ekin_0         ! начальная кинетическая энергия
  real(R8):: epot_0         ! начальная потенциальная энергия
  real(R8):: etot_0         ! начальная полная энергия
  real(R8):: vol_0
  real(R8):: totmass_0
  real(R8):: vol
  real(R8):: maxt           ! ограничение на расчетное время (или 0)
  real(R8):: tprint         ! частота выдачи графики по времени (или 0)

  real(R8) :: talpha, tbeta, b0, q0, maxSound
  real(R8), pointer :: ubcl0(:), ubcr0(:)

  real(R8),pointer:: x(:), y(:)                       ! координаты прямоугольной сетки
  real(R8),pointer:: c_b(:,:),fx_b(:,:),fy_b(:,:)     ! уровень дна
  integer(1),pointer:: c_map(:,:)                     ! карта
  !real(R8),pointer:: c_z(:,:,:),cs_z(:,:,:)   ! z-координата раздела между слоями
  real(R8),pointer:: fx_z(:,:,:),fy_z(:,:,:),fxn_z(:,:,:),fyn_z(:,:,:)
  real(R8),pointer:: fy_tmp(:,:,:)

  ! переменные мелкой воды (в названии объекта присутствует 0):
  real(R8), pointer:: c_z0(:,:),   cs_z0(:,:),   cn_z0(:,:)                     ! уровень поверхности МВ
  real(R8), pointer:: c_h0(:,:),   cs_h0(:,:),   cn_h0(:,:)                     ! глубина жидкости в ячейках
  real(R8), pointer:: c_rho0(:,:), cs_rho0(:,:), cn_rho0(:,:)                 ! плотность в ячейках
  real(R8), pointer:: c_u0(:,:),   cs_u0(:,:),   cn_u0(:,:)                     ! скорость в направлении X в ячейках
  real(R8), pointer:: c_v0(:,:),   cs_v0(:,:),   cn_v0(:,:)                     ! скорость в направлении Y в ячейках
  real(R8), pointer:: c_w0(:,:,:), cs_w0(:,:,:), cn_w0(:,:,:)
  real(R8), pointer:: c_wt0(:,:),  cs_wt0(:,:),  cn_wt0(:,:)                    ! вертикальная скорость на поверхности в ячейках
  real(R8), pointer:: c_wb0(:,:),  cs_wb0(:,:),  cn_wb0(:,:)                    ! вертикальная скорость на дне в ячейках
  real(R8), pointer:: cs_G0(:,:)                                  ! G = sqrt(g/h)

  real(R8), pointer:: c_sumFluxH0(:,:)                          ! суммарный поток в ячейке для h
  real(R8), pointer:: c_sumFluxM0(:,:)                          ! суммарный поток в ячейке для h*rho
  real(R8), pointer:: c_sumFluxIx0(:,:)                         ! суммарный поток в ячейке для h*rho*u
  real(R8), pointer:: c_sumFluxIy0(:,:)                         ! суммарный поток в ячейке для h*rho*v

  real(R8), pointer:: fx_z0(:,:), fxn_z0(:,:)
  real(R8), pointer:: fx_h0(:,:), fxn_h0(:,:)
  real(R8), pointer:: fx_rho0(:,:), fxn_rho0(:,:)
  real(R8), pointer:: fx_u0(:,:), fxn_u0(:,:)
  real(R8), pointer:: fx_v0(:,:), fxn_v0(:,:)
  real(R8), pointer:: fx_w0(:,:,:), fxn_w0(:,:,:)
  real(R8), pointer:: fx_wb0(:,:), fxn_wb0(:,:)
  real(R8), pointer:: fx_wt0(:,:), fxn_wt0(:,:)

  real(R8), pointer:: fy_z0(:,:), fyn_z0(:,:)
  real(R8), pointer:: fy_h0(:,:), fyn_h0(:,:)
  real(R8), pointer:: fy_rho0(:,:), fyn_rho0(:,:)
  real(R8), pointer:: fy_u0(:,:), fyn_u0(:,:)
  real(R8), pointer:: fy_v0(:,:), fyn_v0(:,:)
  real(R8), pointer:: fy_w0(:,:,:), fyn_w0(:,:,:)
  real(R8), pointer:: fy_wb0(:,:), fyn_wb0(:,:)
  real(R8), pointer:: fy_wt0(:,:), fyn_wt0(:,:)

  real(R8), pointer:: fx_fluxH0(:,:),  fy_fluxH0(:,:)
  real(R8), pointer:: fx_fluxM0(:,:),  fy_fluxM0(:,:)
  real(R8), pointer:: fx_fluxIx0(:,:), fy_fluxIx0(:,:)
  real(R8), pointer:: fx_fluxIy0(:,:), fy_fluxIy0(:,:)

  real(R8),pointer:: c_zt(:,:),cs_zt(:,:)                         ! z-координата реальной поверхности
  real(R8), pointer:: fz_w0(:,:,:), fzn_w0(:,:,:)

  real(R8),pointer::f_tmp3(:,:,:)                                 ! "пустой" указатель - только для обмена 3D указателей

  ! консервативные переменные в слой-ячейках (по z - в середине слоя):
  real(R8),pointer:: c_du(:,:,:), cs_du(:,:,:), c_u(:,:,:), cs_u(:,:,:)
  real(R8),pointer:: c_dv(:,:,:), cs_dv(:,:,:), c_v(:,:,:), cs_v(:,:,:)
  real(R8),pointer:: c_dw(:,:,:), cs_dw(:,:,:), c_w(:,:,:), cs_w(:,:,:)
  real(R8),pointer:: c_drho(:,:,:), cs_drho(:,:,:), c_rho(:,:,:), cs_rho(:,:,:)
  real(R8),pointer:: c_dteta(:,:,:), cs_dteta(:,:,:), c_teta(:,:,:), cs_teta(:,:,:)

  real(R8),pointer:: cs_G(:,:,:)                                       ! константа в инвариантах R и Q
  real(R8),pointer:: cs_zp(:,:,:)                         ! z с точкой

  ! данные на X-слой-гранях:
  real(R8),pointer:: fx_du(:,:,:), fxn_du(:,:,:), fx_u(:,:,:), fxn_u(:,:,:)
  real(R8),pointer:: fx_dv(:,:,:), fxn_dv(:,:,:), fx_v(:,:,:), fxn_v(:,:,:)
  real(R8),pointer:: fx_dw(:,:,:), fxn_dw(:,:,:), fx_w(:,:,:), fxn_w(:,:,:)
  real(R8),pointer:: fx_dteta(:,:,:), fxn_dteta(:,:,:), fx_teta(:,:,:), fxn_teta(:,:,:)
  real(R8),pointer:: fx_drho(:,:,:), fxn_drho(:,:,:), fx_rho(:,:,:), fxn_rho(:,:,:)
  real(R8),pointer:: fx_dH(:,:), fxn_dH(:,:)                                      ! для поверхности fx_dH = fx_z(1) - fx_z0

  ! данные на Y-слой-гранях:
  real(R8),pointer:: fy_du(:,:,:), fyn_du(:,:,:), fy_u(:,:,:), fyn_u(:,:,:)
  real(R8),pointer:: fy_dv(:,:,:), fyn_dv(:,:,:), fy_v(:,:,:), fyn_v(:,:,:)
  real(R8),pointer:: fy_dw(:,:,:), fyn_dw(:,:,:), fy_w(:,:,:), fyn_w(:,:,:)
  real(R8),pointer:: fy_dteta(:,:,:), fyn_dteta(:,:,:), fy_teta(:,:,:), fyn_teta(:,:,:)
  real(R8),pointer:: fy_drho(:,:,:), fyn_drho(:,:,:), fy_rho(:,:,:), fyn_rho(:,:,:)
  real(R8),pointer:: fy_dH(:,:), fyn_dH(:,:)                                      ! для поверхности fy_dH = fy_z(1) - fy_z0

  ! данные на Z-гранях:
  real(R8),pointer:: fz_zp(:,:,:), fzn_zp(:,:,:)                ! z с точкой
  real(R8),pointer:: fz_z(:,:,:), fzn_z(:,:,:), fzs_z(:,:,:)    ! z-координата раздела слоёв
  real(R8),pointer:: fz_du(:,:,:), fzn_du(:,:,:), fz_u(:,:,:), fzn_u(:,:,:)
  real(R8),pointer:: fz_dv(:,:,:), fzn_dv(:,:,:), fz_v(:,:,:), fzn_v(:,:,:)
  real(R8),pointer:: fz_dw(:,:,:), fzn_dw(:,:,:), fz_w(:,:,:), fzn_w(:,:,:)
  real(R8),pointer:: fz_dteta(:,:,:), fzn_dteta(:,:,:), fz_teta(:,:,:), fzn_teta(:,:,:)
  real(R8),pointer:: fz_drho(:,:,:), fzn_drho(:,:,:), fz_rho(:,:,:), fzn_rho(:,:,:)

  real(R8),pointer:: c_tmp(:,:,:), DIMUC(:,:,:), DIMVC(:,:,:)

  type(TBcData), pointer:: bcData(:)                            ! данные для ГУ "вход". Индекс = fx_bc или fy_bc
  integer:: nBcData                                             ! количество элементов в bcData
  integer(4), pointer:: fx_bc(:,:), fy_bc(:,:)                  ! индекс данных с граничными условиями (индекс в bcData) или 0

  ! для задачи БМ:
  real(R8):: dzg0, dzg1                                         ! колебания уровня в двух точках горла
  real(R8), pointer:: ksiWS(:,:)                                ! относительное положение граней на линии поперёк Горла
  logical:: needCalc_BC_H_VAR_T

  type(TopexData) topex(100)                                    ! данные приливов
  type(Vector2):: gspos(2)                                      ! кординаты точек на линии поперёк Горла

  real(R8), pointer :: c_dx(:), c_dy(:)                         ! шаги сетки (размеры ячеек в прямоугольной сетке)
  integer(1), pointer:: c_type(:,:)                             ! типы ячеек
  integer(1), pointer:: fx_type(:,:), fy_type(:,:)              ! типы граней

  ! реальное размещение потоковых переменных:
  real(R8), allocatable, target:: fx_dteta1(:,:,:), fx_du1(:,:,:), fx_dv1(:,:,:), fx_dw1(:,:,:), fx_drho1(:,:,:), &
                                  fx_dteta2(:,:,:), fx_du2(:,:,:), fx_dv2(:,:,:), fx_dw2(:,:,:), fx_drho2(:,:,:)
  real(R8), allocatable, target:: fy_dteta1(:,:,:), fy_du1(:,:,:), fy_dv1(:,:,:), fy_dw1(:,:,:), fy_drho1(:,:,:), &
                                  fy_dteta2(:,:,:), fy_du2(:,:,:), fy_dv2(:,:,:), fy_dw2(:,:,:), fy_drho2(:,:,:)
  real(R8), allocatable, target:: fz_dteta1(:,:,:), fz_du1(:,:,:), fz_dv1(:,:,:), fz_dw1(:,:,:), fz_drho1(:,:,:), &
                                  fz_dteta2(:,:,:), fz_du2(:,:,:), fz_dv2(:,:,:), fz_dw2(:,:,:), fz_drho2(:,:,:)

  real(R8), allocatable, target:: fx_teta1(:,:,:), fx_u1(:,:,:), fx_v1(:,:,:), fx_w1(:,:,:), fx_rho1(:,:,:), &
                                  fx_teta2(:,:,:), fx_u2(:,:,:), fx_v2(:,:,:), fx_w2(:,:,:), fx_rho2(:,:,:)
  real(R8), allocatable, target:: fy_teta1(:,:,:), fy_u1(:,:,:), fy_v1(:,:,:), fy_w1(:,:,:), fy_rho1(:,:,:), &
                                  fy_teta2(:,:,:), fy_u2(:,:,:), fy_v2(:,:,:), fy_w2(:,:,:), fy_rho2(:,:,:)
  real(R8), allocatable, target:: fz_teta1(:,:,:), fz_u1(:,:,:), fz_v1(:,:,:), fz_w1(:,:,:), fz_rho1(:,:,:), &
                                  fz_teta2(:,:,:), fz_u2(:,:,:), fz_v2(:,:,:), fz_w2(:,:,:), fz_rho2(:,:,:)

  real(R8), allocatable, target:: fx_z1(:,:,:), fx_zp1(:,:,:), &
                                  fx_z2(:,:,:), fx_zp2(:,:,:)
  real(R8), allocatable, target:: fy_z1(:,:,:), fy_zp1(:,:,:), &
                                  fy_z2(:,:,:), fy_zp2(:,:,:)
  real(R8), allocatable, target:: fz_z1(:,:,:), fz_zp1(:,:,:), &
                                  fz_z2(:,:,:), fz_zp2(:,:,:)

  ! реальное размещение потоковых переменных МЕЛКОЙ ВОДЫ:
  real(R8), allocatable, target:: fx_z01(:,:), fx_h01(:,:), fx_rho01(:,:), fx_u01(:,:), fx_v01(:,:), fx_wb01(:,:), fx_wt01(:,:), &
                                  fx_z02(:,:), fx_h02(:,:), fx_rho02(:,:), fx_u02(:,:), fx_v02(:,:), fx_wb02(:,:), fx_wt02(:,:)
  real(R8), allocatable, target:: fy_z01(:,:), fy_h01(:,:), fy_rho01(:,:), fy_u01(:,:), fy_v01(:,:), fy_wb01(:,:), fy_wt01(:,:), &
                                  fy_z02(:,:), fy_h02(:,:), fy_rho02(:,:), fy_u02(:,:), fy_v02(:,:), fy_wb02(:,:), fy_wt02(:,:)
  real(R8), allocatable, target:: fx_w01(:,:,:), fy_w01(:,:,:), fz_w01(:,:,:), &
                                  fx_w02(:,:,:), fy_w02(:,:,:), fz_w02(:,:,:)

  ! описание сетки:
  type PAIR; integer:: i, j; end type PAIR                          ! тип: пара целых индексов
  type PAIRS; type(PAIR), pointer:: ptr(:); end type PAIRS          ! тип: список пар индексов

  type(PAIRS):: cells                                               ! список индексов внутренних ячеек
  type(PAIRS):: fx_sides(0:NBC)                                     ! fx_sides(ibc) - список X-граней с ГУ ibc
  type(PAIRS):: fy_sides(0:NBC)                                     ! fy_sides(ibc) - список Y-граней с ГУ ibc
  integer:: ncells                                                  ! число внутренних ячеек
  integer:: nfx_sides(0:NBC)                                        ! nfx_sides(ibc) - число X-граней с ГУ ibc
  integer:: nfy_sides(0:NBC)                                        ! nfy_sides(ibc) - число Y-граней с ГУ ibc

  ! отладочная информация:
  real(R8) :: debdum
  real(R8), pointer :: deby(:), debz1(:), debz2(:)

  integer(1), pointer :: debCells(:,:,:), debSidesX(:,:,:), debSidesY(:,:,:), debSidesZ(:,:,:)  ! ячейки и грани под отладкой
  integer(4) :: debStartStep, debEndStep, debStartZ, debEndZ  ! ограничения на шаги и слои под отладкой

  character(len=:),allocatable::frm619,frm9997,frm9998,frm9999
  character(len=:),allocatable::frm623,frm622,frm621,frm620,frm624
  character(len=:),allocatable::frm617,frm618,frm9622

  !character (len=120) :: name
  character (len=10), allocatable :: colorArray(:)

  real(R8), parameter:: eps = 1.d-5
  real(R8), parameter:: eps1 = 1.d-7
  real(R8), parameter:: PI = 3.14159265358979d0

contains !-----------------------------------------------------------------------------------------------------

  subroutine allocate

    if(NX<=1 .or. NY<=1 .or. NZ<=0 ) then
      write(*,*) "ERROR Allocate(). NX, NY, NZ:", NX, NY, NZ
      stop
    end if

    ncx = nx - 1; ncy = ny - 1
    nxx = nx    ; nxy = ny - 1
    nyx = nx - 1; nyy = ny
    ncz = nz -1

    colorArray = [character (len=10) :: "BLACK", "RED", "GREEN","BLUE", "CYAN", "YELLOW", "PURPLE"]

    !-- переменные мелкой воды (в названии объекта присутствует0): ---------------------
    allocate(c_sumFluxH0(ncx,ncy))                          ! суммарный поток в ячейке для h
    allocate(c_sumFluxM0(ncx,ncy))                          ! суммарный поток в ячейке для h*rho
    allocate(c_sumFluxIx0(ncx,ncy))                         ! суммарный поток в ячейке для h*rho*u
    allocate(c_sumFluxIy0(ncx,ncy))                         ! суммарный поток в ячейке для h*rho*v
    allocate(cs_G0(ncx,ncy))                                  ! G = sqrt(g/h)

    allocate(fx_fluxH0(nxx,nxy),  fy_fluxH0(nyx,nyy))
    allocate(fx_fluxM0(nxx,nxy),  fy_fluxM0(nyx,nyy))
    allocate(fx_fluxIx0(nxx,nxy), fy_fluxIx0(nyx,nyy))
    allocate(fx_fluxIy0(nxx,nxy), fy_fluxIy0(nyx,nyy))

    allocate(c_z0(ncx,ncy),   cs_z0(ncx,ncy),   cn_z0(ncx,ncy))                     ! уровень поверхности МВ
    allocate(c_h0(ncx,ncy),   cs_h0(ncx,ncy),   cn_h0(ncx,ncy))                     ! глубина жидкости в ячейках
    allocate(c_rho0(ncx,ncy), cs_rho0(ncx,ncy), cn_rho0(ncx,ncy))                   ! плотность в ячейках
    allocate(c_u0(ncx,ncy),   cs_u0(ncx,ncy),   cn_u0(ncx,ncy))                     ! скорость в направлении X в ячейках
    allocate(c_v0(ncx,ncy),   cs_v0(ncx,ncy),   cn_v0(ncx,ncy))                     ! скорость в направлении Y в ячейках
    allocate(c_wt0(ncx,ncy),  cs_wt0(ncx,ncy),  cn_wt0(ncx,ncy))                    ! вертикальная скорость на поверхности в ячейках
    allocate(c_wb0(ncx,ncy),  cs_wb0(ncx,ncy),  cn_wb0(ncx,ncy))                    ! вертикальная скорость на дне в ячейках

    allocate(c_w0(ncx,ncy,ncz), cs_w0(ncx,ncy,ncz), cn_w0(ncx,ncy,ncz))

    !----------------------------------------------------------------------------------------

    allocate(c_b(ncx,ncy), c_type(ncx, ncy), c_tmp(ncx,ncy,ncz), c_dx(ncx), c_dy(ncy),          &
        c_dteta(ncx,ncy,ncz), cs_dteta(ncx,ncy,ncz), c_teta(ncx,ncy,ncz), cs_teta(ncx,ncy,ncz), &
        c_drho(ncx,ncy,ncz),  cs_drho(ncx,ncy,ncz),  c_rho(ncx,ncy,ncz),  cs_rho(ncx,ncy,ncz),  &
        c_du(ncx,ncy,ncz),    cs_du(ncx,ncy,ncz),    c_u(ncx,ncy,ncz),    cs_u(ncx,ncy,ncz),    &
        c_dv(ncx,ncy,ncz),    cs_dv(ncx,ncy,ncz),    c_v(ncx,ncy,ncz),    cs_v(ncx,ncy,ncz),    &
        c_dw(ncx,ncy,ncz),    cs_dw(ncx,ncy,ncz),    c_w(ncx,ncy,ncz),    cs_w(ncx,ncy,ncz),    &
        c_zt(ncx,ncy),        cs_zt(ncx,ncy),   &
                              cs_G(ncx,ncy,ncz) &
    )

    allocate(fxn_dH(nxx,nxy), fyn_dH(nyx,nyy), fzs_z(ncx,ncy,nz))

    allocate(x(nxx), fx_b(nxx,nxy), fx_type(nxx,nxy), fx_bc(nx,nxy), &
        fx_z1(nxx,nxy,nz), fx_zp1(nxx,nxy,ncz), fx_dteta1(nxx,nxy,ncz), fx_du1(nxx,nxy,ncz), fx_dv1(nxx,nxy,ncz), fx_dw1(nxx,nxy,ncz), fx_drho1(nxx,nxy,ncz), &
        fx_z2(nxx,nxy,nz), fx_zp2(nxx,nxy,ncz), fx_dteta2(nxx,nxy,ncz), fx_du2(nxx,nxy,ncz), fx_dv2(nxx,nxy,ncz), fx_dw2(nxx,nxy,ncz), fx_drho2(nxx,nxy,ncz), &
        fx_teta1(nxx,nxy,ncz), fx_u1(nxx,nxy,ncz), fx_v1(nxx,nxy,ncz), fx_w1(nxx,nxy,ncz), fx_rho1(nxx,nxy,ncz), &
        fx_teta2(nxx,nxy,ncz), fx_u2(nxx,nxy,ncz), fx_v2(nxx,nxy,ncz), fx_w2(nxx,nxy,ncz), fx_rho2(nxx,nxy,ncz) &
    )

    allocate(y(nyy), fy_b(nyx,nyy), fy_type(nyx,nyy), fy_bc(nyx,nyy), &
        fy_z1(nyx,nyy,nz), fy_dteta1(nyx,nyy,ncz), fy_du1(nyx,nyy,ncz), fy_dv1(nyx,nyy,ncz), fy_dw1(nyx,nyy,ncz), fy_drho1(nyx,nyy,ncz), &
        fy_z2(nyx,nyy,nz), fy_dteta2(nyx,nyy,ncz), fy_du2(nyx,nyy,ncz), fy_dv2(nyx,nyy,ncz), fy_dw2(nyx,nyy,ncz), fy_drho2(nyx,nyy,ncz), &
        fy_teta1(nyx,nyy,ncz), fy_u1(nyx,nyy,ncz), fy_v1(nyx,nyy,ncz), fy_w1(nyx,nyy,ncz), fy_rho1(nyx,nyy,ncz), &
        fy_teta2(nyx,nyy,ncz), fy_u2(nyx,nyy,ncz), fy_v2(nyx,nyy,ncz), fy_w2(nyx,nyy,ncz), fy_rho2(nyx,nyy,ncz) &
    )

    allocate( &
        fz_z1(ncx,ncy,nz), fz_dteta1(ncx,ncy,nz), fz_du1(ncx,ncy,nz), fz_dv1(ncx,ncy,nz), fz_dw1(ncx,ncy,nz), fz_drho1(ncx,ncy,nz), &
        fz_z2(ncx,ncy,nz), fz_dteta2(ncx,ncy,nz), fz_du2(ncx,ncy,nz), fz_dv2(ncx,ncy,nz), fz_dw2(ncx,ncy,nz), fz_drho2(ncx,ncy,nz), &
        fz_zp1(ncx,ncy,nz), fz_teta1(ncx,ncy,nz), fz_u1(ncx,ncy,nz), fz_v1(ncx,ncy,nz), fz_w1(ncx,ncy,nz), fz_rho1(ncx,ncy,nz), &
        fz_zp2(ncx,ncy,nz), fz_teta2(ncx,ncy,nz), fz_u2(ncx,ncy,nz), fz_v2(ncx,ncy,nz), fz_w2(ncx,ncy,nz), fz_rho2(ncx,ncy,nz) &
    )

    ! для мелкой воды:
    allocate(fx_z01(nxx,nxy), fx_h01(nxx,nxy), fx_rho01(nxx,nxy), fx_u01(nxx,nxy), fx_v01(nxx,nxy), &
             fx_z02(nxx,nxy), fx_h02(nxx,nxy), fx_rho02(nxx,nxy), fx_u02(nxx,nxy), fx_v02(nxx,nxy))
    allocate(fy_z01(nyx,nyy), fy_h01(nyx,nyy), fy_rho01(nyx,nyy), fy_u01(nyx,nyy), fy_v01(nyx,nyy), &
             fy_z02(nyx,nyy), fy_h02(nyx,nyy), fy_rho02(nyx,nyy), fy_u02(nyx,nyy), fy_v02(nyx,nyy))

    allocate(fx_wb01(nxx,nxy), fx_wt01(nxx,nxy), fx_w01(nxx,nxy,ncz), &
             fx_wb02(nxx,nxy), fx_wt02(nxx,nxy), fx_w02(nxx,nxy,ncz))
    allocate(fy_wb01(nyx,nyy), fy_wt01(nyx,nyy), fy_w01(nyx,nyy,ncz), &
             fy_wb02(nyx,nyy), fy_wt02(nyx,nyy), fy_w02(nyx,nyy,ncz))
    allocate(fz_w01(ncx,ncy,nz), &
             fz_w02(ncx,ncy,nz))

    !allocate(c_types(ncx,ncy), fx_types(nxx,nxy), fy_types(nyx,nyy))  ! типы объектов

    allocate (fy_tmp(nyx,nyy,nz))

    allocate(Tsurf(ncx,ncy))

    allocate(ubcl0(nz-1),ubcr0(nz-1))

    allocate(deby(1:ny),debz1(1:nz),debz2(1:nz))
    allocate(debCells(ncx,ncy,ncz), debSidesX(nxx,nxy,ncz), debSidesY(nyx,nyy,ncz), debSidesZ(ncx,ncy,nz))

    call ExchangePointers

    NAN = IEEE_VALUE(1., IEEE_QUIET_NAN)

    ! размерные коэффициенты:
    CFD_T = 1.; CFD_L = 1.; CFD_H = 1.; CFD_TEMP = 1.

    sigma2Ph=1.
    sigma = 0.5

    debCells=0; debSidesX=0; debSidesY=0
    debStartStep = 0; debEndStep = -1
    debStartZ = 0; debEndZ = 99999

    ! стандартные типы ячеек и граней для прямоугольника:
    c_type = CELL_INNER

    fx_type(1,:) = bcTypeX
    fx_type(2:nx-1,:) = BC_INNER
    fx_type(nx,:) = bcTypeX
    fx_bc = 0                             ! нет ГУ с данными на границе

    fy_type(:,1) = bcTypeY
    fy_type(:,2:ny-1) = BC_INNER
    fy_type(:,ny) = bcTypeY
    fy_bc = 0                             ! нет ГУ с данными на границе

    needCalc_BC_H_VAR_T = .false.

  end subroutine allocate

  !----------------------------------------------------------------------------------------------------

  subroutine ExchangePointers
    implicit none
    logical, save:: isOddEx2 = .true.

    if(isOddEx2) then

      fx_dteta => fx_dteta1; fxn_dteta => fx_dteta2; fx_teta => fx_teta1; fxn_teta => fx_teta2
      fx_drho  => fx_drho1;  fxn_drho  => fx_drho2;  fx_rho  => fx_rho1;  fxn_rho  => fx_rho2
      fx_du    => fx_du1;    fxn_du    => fx_du2;    fx_u    => fx_u1;    fxn_u    => fx_u2
      fx_dv    => fx_dv1;    fxn_dv    => fx_dv2;    fx_v    => fx_v1;    fxn_v    => fx_v2
      fx_dw    => fx_dw1;    fxn_dw    => fx_dw2;    fx_w    => fx_w1;    fxn_w    => fx_w2
      fx_z     => fx_z1;     fxn_z     => fx_z2

      fy_dteta => fy_dteta1; fyn_dteta => fy_dteta2; fy_teta => fy_teta1; fyn_teta => fy_teta2
      fy_drho  => fy_drho1;  fyn_drho  => fy_drho2;  fy_rho  => fy_rho1;  fyn_rho  => fy_rho2
      fy_du    => fy_du1;    fyn_du    => fy_du2;    fy_u    => fy_u1;    fyn_u    => fy_u2
      fy_dv    => fy_dv1;    fyn_dv    => fy_dv2;    fy_v    => fy_v1;    fyn_v    => fy_v2
      fy_dw    => fy_dw1;    fyn_dw    => fy_dw2;    fy_w    => fy_w1;    fyn_w    => fy_w2
      fy_z     => fy_z1;     fyn_z     => fy_z2

      fz_dteta => fz_dteta1; fzn_dteta => fz_dteta2; fz_teta => fz_teta1; fzn_teta => fz_teta2
      fz_drho  => fz_drho1;  fzn_drho  => fz_drho2;  fz_rho  => fz_rho1;  fzn_rho  => fz_rho2
      fz_du    => fz_du1;    fzn_du    => fz_du2;    fz_u    => fz_u1;    fzn_u    => fz_u2
      fz_dv    => fz_dv1;    fzn_dv    => fz_dv2;    fz_v    => fz_v1;    fzn_v    => fz_v2
      fz_dw    => fz_dw1;    fzn_dw    => fz_dw2;    fz_w    => fz_w1;    fzn_w    => fz_w2
      fz_z     => fz_z1;     fzn_z     => fz_z2
      fz_zp    => fz_zp1;    fzn_zp    => fz_zp2

      ! переменные мелкой воды:
      fx_z0   => fx_z01;   fxn_z0   => fx_z02
      fx_h0   => fx_h01;   fxn_h0   => fx_h02
      fx_rho0 => fx_rho01; fxn_rho0 => fx_rho02
      fx_u0   => fx_u01;   fxn_u0   => fx_u02
      fx_v0   => fx_v01;   fxn_v0   => fx_v02
      fx_wb0  => fx_wb01;  fxn_wb0  => fx_wb02
      fx_wt0  => fx_wt01;  fxn_wt0  => fx_wt02
      fx_w0   => fx_w01;   fxn_w0   => fx_w02

      fy_z0   => fy_z01;   fyn_z0   => fy_z02
      fy_h0   => fy_h01;   fyn_h0   => fy_h02
      fy_rho0 => fy_rho01; fyn_rho0 => fy_rho02
      fy_u0   => fy_u01;   fyn_u0   => fy_u02
      fy_v0   => fy_v01;   fyn_v0   => fy_v02
      fy_wb0  => fy_wb01;  fyn_wb0  => fy_wb02
      fy_wt0  => fy_wt01;  fyn_wt0  => fy_wt02
      fy_w0   => fy_w01;   fyn_w0   => fy_w02

      fz_w0   => fz_w01;   fzn_w0   => fz_w02

    else

      fx_dteta => fx_dteta2; fxn_dteta => fx_dteta1; fx_teta => fx_teta2; fxn_teta => fx_teta1
      fx_drho  => fx_drho2;  fxn_drho  => fx_drho1;  fx_rho  => fx_rho2;  fxn_rho  => fx_rho1
      fx_du    => fx_du2;    fxn_du    => fx_du1;    fx_u    => fx_u2;    fxn_u    => fx_u1
      fx_dv    => fx_dv2;    fxn_dv    => fx_dv1;    fx_v    => fx_v2;    fxn_v    => fx_v1
      fx_dw    => fx_dw2;    fxn_dw    => fx_dw1;    fx_w    => fx_w2;    fxn_w    => fx_w1
      fx_z     => fx_z2;     fxn_z     => fx_z1

      fy_dteta => fy_dteta2; fyn_dteta => fy_dteta1; fy_teta => fy_teta2; fyn_teta => fy_teta1
      fy_drho  => fy_drho2;  fyn_drho  => fy_drho1;  fy_rho  => fy_rho2;  fyn_rho  => fy_rho1
      fy_du    => fy_du2;    fyn_du    => fy_du1;    fy_u    => fy_u2;    fyn_u    => fy_u1
      fy_dv    => fy_dv2;    fyn_dv    => fy_dv1;    fy_v    => fy_v2;    fyn_v    => fy_v1
      fy_dw    => fy_dw2;    fyn_dw    => fy_dw1;    fy_w    => fy_w2;    fyn_w    => fy_w1
      fy_z     => fy_z2;     fyn_z     => fy_z1

      fz_dteta => fz_dteta2; fzn_dteta => fz_dteta1; fz_teta => fz_teta2; fzn_teta => fz_teta1
      fz_drho  => fz_drho2;  fzn_drho  => fz_drho1;  fz_rho  => fz_rho2;  fzn_rho  => fz_rho1
      fz_du    => fz_du2;    fzn_du    => fz_du1;    fz_u    => fz_u2;    fzn_u    => fz_u1
      fz_dv    => fz_dv2;    fzn_dv    => fz_dv1;    fz_v    => fz_v2;    fzn_v    => fz_v1
      fz_dw    => fz_dw2;    fzn_dw    => fz_dw1;    fz_w    => fz_w2;    fzn_w    => fz_w1
      fz_z     => fz_z2;     fzn_z     => fz_z1
      fz_zp    => fz_zp2;    fzn_zp    => fz_zp1

      ! переменные мелкой воды:
      fx_z0   => fx_z02;   fxn_z0   => fx_z01
      fx_h0   => fx_h02;   fxn_h0   => fx_h01
      fx_rho0 => fx_rho02; fxn_rho0 => fx_rho01
      fx_u0   => fx_u02;   fxn_u0   => fx_u01
      fx_v0   => fx_v02;   fxn_v0   => fx_v01
      fx_wb0  => fx_wb02;  fxn_wb0  => fx_wb01
      fx_wt0  => fx_wt02;  fxn_wt0  => fx_wt01
      fx_w0   => fx_w02;   fxn_w0   => fx_w01

      fy_z0   => fy_z02;   fyn_z0   => fy_z01
      fy_h0   => fy_h02;   fyn_h0   => fy_h01
      fy_rho0 => fy_rho02; fyn_rho0 => fy_rho01
      fy_u0   => fy_u02;   fyn_u0   => fy_u01
      fy_v0   => fy_v02;   fyn_v0   => fy_v01
      fy_wb0  => fy_wb02;  fyn_wb0  => fy_wb01
      fy_wt0  => fy_wt02;  fyn_wt0  => fy_wt01
      fy_w0   => fy_w02;   fyn_w0   => fy_w01

      fz_w0   => fz_w02;   fzn_w0   => fz_w01

    end if

    isOddEx2 = .not. isOddEx2

  end subroutine ExchangePointers

  !----------------------------------------------------------------------------------------------------

  ! Распределение температуры в тесте Залесного(test_numb == 1)
  function getT(z) result(T)
    real(R8), intent(in) :: z ! input
    real(R8)             :: T ! output
    T = (7.5 *(1. - tanh((z-80.)/30.)) + (5000.-z) / 5000. * 10.)
  end function getT

  !----------------------------------------------------------------------------------------------------

  ! Зависимость плотности от температуры:
  function getRHO(T) result(RHO)
    real(R8), intent(in) :: T ! input
    real(R8)             :: RHO ! output
    RHO = cfT_rho0 * ( 1. - cfT_alpha * T )
  end function getRHO

  !----------------------------------------------------------------------------------------------------

  function getTfromRHO(RHO) result(T)
    real(R8), intent(in) :: RHO ! input
    real(R8)             :: T ! output
    T = (1. - RHO/cfT_rho0) / cfT_alpha
  end function getTfromRHO

  !----------------------------------------------------------------------------------------------------

  ! для теста 4

  ! глубина в ячейке:
  function task4_h(xi)
    implicit none
    real(R8) :: xi, task4_h

    task4_h = 1. + (pi/2. + atan(50./20.*(xi-20./3.)) ) / pi -exp(-50./20.**2*(xi-3./4.*20.)**2)

  end function

  !----------------------------------------------------------------------------------------------------

  ! уровень дна в ячейке:
  function task4_b(xi)
    implicit none
    real(R8) :: xi, task4_b
    real(R8) :: hi

    hi = task4_h(xi)
    task4_b = b0 - hi  - (talpha*tbeta)**2 / (2. * g * (sin(tbeta*hi))**2)
  end function

  !----------------------------------------------------------------------------------------------------

  ! скорость (решение в тесте 4)
  function task4_u(xi, zi)
    implicit none
    real(R8) :: xi, zi, task4_u
    real(R8) :: hi, bi

    hi = task4_h(xi)
    bi = task4_b(xi)
    task4_u = talpha * tbeta  / sin(tbeta * hi) * cos(tbeta * (zi - bi))
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! есть ли справа от грани рабочая ячейка
  function fx_IsOwnRight(i, j)
    integer(4):: i, j
    logical:: fx_IsOwnRight
    if(i==ncx) then                                             ! если грань на левом краю области
      fx_IsOwnRight = .false.                                   ! .. справа ячейки быть не может
    else
      fx_IsOwnRight = c_type(i, j) > CELL_DELETED
    end if
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! есть ли справа от грани рабочая ячейка
  function fy_IsOwnRight(i, j)
    integer(4):: i, j
    logical:: fy_IsOwnRight
    if(j==ncy) then                                             ! если грань на левом краю области
      fy_IsOwnRight = .false.                                   ! .. справа ячейки быть не может
    else
      fy_IsOwnRight = c_type(i, j) > CELL_DELETED
    end if
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  function cs_GetG(i,j,k)
    integer(4):: i, j, k
    real(R8):: cs_GetG

    cs_GetG = cs_sound / (teta0 + cs_dteta(i,j,k))
  end

  !----------------------------------------------------------------------------------------------------------------------------

  ! вертикальная скорость в слой-ячейке (i,j,k) для t=t[n]
  function c_GetW0(i,j,k)
    integer(4):: i, j, k
    real(R8):: z, c_GetW0, a
    z = (fz_z(i,j,k) + fz_z(i,j,k+1)) / 2.
    a = (z - c_b(i,j)) / c_h0(i,j)
    c_GetW0 = c_wt0(i,j) * a + c_wb0(i,j) * (1. - a)
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! вертикальная скорость в слой-ячейке (i,j,k) для t=t[n+1/2]
  function cs_GetW0(i,j,k)
    integer(4):: i, j, k
    real(R8):: z, cs_GetW0, a, wt, wb
    z = (fzs_z(i,j,k) + fzs_z(i,j,k+1)) / 2.      ! середина слоя t[n+1/2]
    a = (z - c_b(i,j)) / cs_h0(i,j)
    wt = cs_wt0(i,j)
    wb = cs_wb0(i,j)
    cs_GetW0 = wt * a + wb * (1. - a)
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! вертикальная скорость в слой-ячейке (i,j,k) для t=t[n+1]
  function cn_GetW0(i,j,k)
    integer(4):: i, j, k
    real(R8):: z, cn_GetW0, a
    z = (fzn_z(i,j,k) + fzn_z(i,j,k+1)) / 2.      ! середина слоя t[n+1]
    a = (z - c_b(i,j)) / cn_h0(i,j)
    cn_GetW0 = cn_wt0(i,j) * a + cn_wb0(i,j) * (1. - a)
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! вертикальная скорость на X-слой-грани (i,j,k) для t=t[n]
  function fx_GetW0(i,j,k)
    integer(4):: i, j, k
    real(R8):: z, fx_GetW0, a
    z = (fx_z(i,j,k) + fx_z(i,j,k+1)) / 2.
    a = (z - fx_b(i,j)) / fx_h0(i,j)
    fx_GetW0 = fx_wt0(i,j) * a + fx_wb0(i,j) * (1. - a)
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! вертикальная скорость на X-слой-грани (i,j,k) для t=t[n+1]
  function fxn_GetW0(i,j,k)
    integer(4):: i, j, k
    real(R8):: z, fxn_GetW0, a
    z = (fxn_z(i,j,k) + fxn_z(i,j,k+1)) / 2.
    a = (z - fx_b(i,j)) / fxn_h0(i,j)
    fxn_GetW0 = fxn_wt0(i,j) * a + fxn_wb0(i,j) * (1. - a)
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! вертикальная скорость на Y-слой-грани для t=t[n]
  function fy_GetW0(i,j,k)
    integer(4):: i, j, k
    real(R8):: z, fy_GetW0, a
    z = (fy_z(i,j,k) + fy_z(i,j,k+1)) / 2.
    a = (z - fy_b(i,j)) / fy_h0(i,j)
    fy_GetW0 = fy_wt0(i,j) * a + fy_wb0(i,j) * (1. - a)
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! вертикальная скорость на Y-грани для t=t[n+1]
  function fyn_GetW0(i,j,k)
    integer(4):: i, j, k
    real(R8):: z, fyn_GetW0, a
    z = (fyn_z(i,j,k) + fyn_z(i,j,k+1)) / 2.
    a = (z - fy_b(i,j)) / fyn_h0(i,j)
    fyn_GetW0 = fyn_wt0(i,j) * a + fyn_wb0(i,j) * (1. - a)
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! вертикальная скорость на Z-слой-грани для t=t[n]
  function fz_GetW0(i,j,k)
    integer(4):: i, j, k
    real(R8):: z, fz_GetW0, a
    z = fz_z(i,j,k)
    a = (z - c_b(i,j)) / c_h0(i,j)
    fz_GetW0 = c_wt0(i,j) * a + c_wb0(i,j) * (1. - a)
    debdum=0
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! вертикальная скорость на Z-грани для t=t[n+1]
  function fzn_GetW0(i,j,k)
    integer(4):: i, j, k
    real(R8):: z, fzn_GetW0, a
    z = fzn_z(i,j,k)
    a = (z - c_b(i,j)) / cn_h0(i,j)
    fzn_GetW0 = cn_wt0(i,j) * a + cn_wb0(i,j) * (1. - a)
    debdum=0
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  function c_GetRho(i,j,k,in)
    integer(4):: i, j, k, in
    real(R8):: z, c_GetRho
    select case (in)
      case(0); c_GetRho = c_drho(i,j,k) + rho0
      case(1); c_GetRho = cs_drho(i,j,k) + rho0
      case(2); c_GetRho = c_drho(i,j,k) + rho0
      case default; call avost
    end select
  end

  !----------------------------------------------------------------------------------------------------------------------------

  function c_GetTeta(i,j,k,in)
    integer(4):: i, j, k, in
    real(R8):: z, c_GetTeta
    select case (in)
      case(0); c_GetTeta = c_dteta(i,j,k) + teta0
      case(1); c_GetTeta = cs_dteta(i,j,k) + teta0
      case(2); c_GetTeta = c_dteta(i,j,k) + teta0
      case default; call avost
    end select
  end

  !----------------------------------------------------------------------------------------------------------------------------

  function c_GetU(i,j,k,in)
    integer(4):: i, j, k, in
    real(R8):: z, c_GetU
    select case (in)
      case(0); c_GetU = c_du(i,j,k) + c_u0(i,j)
      case(1); c_GetU = cs_du(i,j,k) + cs_u0(i,j)
      case(2); c_GetU = c_du(i,j,k) + cn_u0(i,j)
      case default; call avost
    end select
  end

  !----------------------------------------------------------------------------------------------------------------------------

  function c_GetV(i,j,k,in)
    integer(4):: i, j, k, in
    real(R8):: z, c_GetV
    select case (in)
      case(0); c_GetV = c_dv(i,j,k) + c_v0(i,j)
      case(1); c_GetV = cs_dv(i,j,k) + cs_v0(i,j)
      case(2); c_GetV = c_dv(i,j,k) + cn_v0(i,j)
      case default; call avost
    end select
  end

  !----------------------------------------------------------------------------------------------------------------------------

  function c_GetW(i,j,k,in)
    integer(4):: i, j, k, in
    real(R8):: z, c_GetW
    select case (in)
      case(0); c_GetW = c_dw(i,j,k) + c_GetW0(i,j,k)
      case(1); c_GetW = cs_dw(i,j,k) + cs_GetW0(i,j,k)
      case(2); c_GetW = c_dw(i,j,k) + cn_GetW0(i,j,k)
      case default; call avost
    end select
  end

  !----------------------------------------------------------------------------------------------------------------------------

  ! Вычисление безразмерной гравитационной постоянной по размерной [Си]: g = gC * CFD_L^2 / (CFD_T^2 * CFD_H)
  function Get0Gravity(gC)
    real(R8) gC, Get0Gravity
    Get0Gravity = gC * CFD_L * CFD_L / (CFD_T * CFD_T * CFD_H)
  end function

  ! вычисление безразмерной длины по длине в Си и коэффициентам обезразмеривания: l = lC * CFD_L
  function Get0Length(lC)
    real(R8) lC, Get0Length
    Get0Length = lC * CFD_L
  end function

  ! вычисление безразмерного времени по времени в Си и коэффициентам обезразмеривания: t = tC * CFD_T
  function Get0Time(tC)
    real(R8) tC, Get0Time
    Get0Time = tC * CFD_T
  end function

  ! вычисление безразмерной скорости по скорости в Си и коэффициентам обезразмеривания: u = uC * CFD_L / CFD_T
  function Get0Velocity(uC)
    real(R8) uC, Get0Velocity
    Get0Velocity = uC * CFD_L / CFD_T
  end function

  ! вычисление безразмерной глубины по глубине в Си и коэффициентам обезразмеривания: h = hC * CFD_H
  function Get0Depth(hC)
    real(R8) hC, Get0Depth
    Get0Depth = hC * CFD_H
  end function

  ! вычисление безразмерной площади по площади в Си и коэффициентам обезразмеривания: s = sC * CFD_L^2
  function Get0Square(sC)
    real(R8) sC, Get0Square
    Get0Square = sC * CFD_L * CFD_L
  end function

  ! вычисление безразмерной плотности по плотности в Си и коэффициентам обезразмеривания: rho = rhoC * CFD_RHO
  function Get0Density(tC)
    real(R8) tC, Get0Density
    Get0Density = tC * CFD_RHO
  end function

  ! вычисление безразмерной температуры: t = tC * CFD_TEMP
  function Get0Temperature(tC)
    real(R8) tC, Get0Temperature
    Get0Temperature = tC * CFD_TEMP
  end function

  ! вычисление безразмерной частоты вращения из частоты [rad/sec]
  function Get0Omega(oC)
    real(R8) oC, Get0Omega
    Get0Omega = oC / CFD_T
  end function

  !----------------------------------------------------------------------------------------------------------------------------

  ! Вычисление размерной гравитационной постоянной по размерной [Си]: gC = gC * (CFD_T^2 * CFD_H) / CFD_L^2
  function GetCGravity(g0)
    real(R8) g0, GetCGravity
    GetCGravity = g0 / (CFD_L * CFD_L) * (CFD_T * CFD_T * CFD_H)
  end function

  ! вычисление размерной [Си] длины по безразмерной длине
  function GetCLength(l0)
    real(R8) l0, GetCLength
    GetCLength = l0 / CFD_L
  end function

  ! вычисление размерного времени по безразмерному времени
  function GetCTime(t0)
    real(R8) t0, GetCTime
    GetCTime = t0 / CFD_T
  end function

  ! вычисление размерной скорости по безразмерной скорости
  function GetCVelocity(u0)
    real(R8) u0, GetCVelocity
    GetCVelocity = u0 / CFD_L * CFD_T
  end function

  ! вычисление размерной скорости по безразмерной скорости
  function GetCVelocityZ(u0)
    real(R8) u0, GetCVelocityZ
    GetCVelocityZ = u0 / CFD_H * CFD_T
  end function

  ! вычисление размерной глубины по безразмерной глубине
  function GetCDepth(h0)
    real(R8) h0, GetCDepth
    GetCDepth = h0 / CFD_H
  end function

  ! вычисление размерной площади по безразмерной площади
  function GetCSquare(s0)
    real(R8) s0, GetCSquare
    GetCSquare = s0 / CFD_L**2
  end function

  function GetCDensity(t0)
    real(R8) t0, GetCDensity
    GetCDensity = t0 / CFD_RHO
  end function

  ! вычисление размерной температуры по безразмерной
  function GetCTemperature(t0)
    real(R8) t0, GetCTemperature
    GetCTemperature = t0 / CFD_TEMP
  end function

  ! вычисление размерной [rad/sec] частоты вращения из безразмерной частоты
  function GetCOmega(oC)
    real(R8) oC, GetCOmega
    GetCOmega = oC * CFD_T
  end function

  !---------------------------------------------------------------------------------------------------------------------

  ! обезразмеривание коэфициентов и сеточных переменных
  subroutine dims
    implicit none

    Coriolis = Coriolis / CFD_T / CFD_L
    cfViscVert = cfViscVert * CFD_H**2 / CFD_T
    cfViscVertT = cfViscVertT * CFD_H**2 / CFD_T
    cf3ViscVertT = cf3ViscVertT * CFD_H / CFD_T  ! почему так ???!!!
    cfBottomFriction = cfBottomFriction * CFD_H / CFD_L
    cfSurfFrictionU = cfSurfFrictionU * CFD_H / CFD_L
    cfSurfFrictionV = cfSurfFrictionV * CFD_H / CFD_L
    sig_0 = sig_0 * CFD_L**2 / CFD_T

    uAtm = Get0Velocity(uAtm)
    g = Get0Gravity(g)
    dl = Get0Length(dl)
    dw = Get0Length(dw)

    x = x * CFD_L
    y = y * CFD_L

    fx_z = fx_z * CFD_H
    fy_z = fy_z * CFD_H
!    hx = hx * CFD_H
!    hy = hy * CFD_H
    c_zt = c_zt * CFD_H
    fx_b = fx_b * CFD_H
    fy_b = fy_b * CFD_H
    c_b = c_b * CFD_H

    !dx = Get0Length(dx)
    !dy = Get0Length(dy)

    c_du = c_du * (CFD_L / CFD_T)
    fx_du = fx_du * (CFD_L / CFD_T)
    fy_du = fy_du * (CFD_L / CFD_T)
  end subroutine

  !---------------------------------------------------------------------------------------------------------------------

  function GradToRad(grad)
    real(8):: grad, GradToRad
    GradToRad = grad / 180.d0 * pi
  end function GradToRad

  !---------------------------------------------------------------------------------------------------------------------

  subroutine DebTest
    integer:: i,j,k
      do i=1,nyx; do j=1,nyy; do k=1,ncz
        if(abs(fy_teta(i,j,k)-fx_teta(j,i,k))/=0.) then
          debdum=0
        end if
        if(abs(fy_u(i,j,k)-fx_v(j,i,k))/=0.) then
          debdum=1
        end if
        if(abs(fy_v(i,j,k)-fx_u(j,i,k))/=0.) then
          debdum=2
        end if
        if(abs(fy_z(i,j,k)-fx_z(j,i,k))/=0.) then
          debdum=3
        end if
      end do; end do; end do
  end subroutine

  !---------------------------------------------------------------------------------------------------------------------

  function IsDeb(iz)
    integer(4) :: iz
    logical :: IsDeb

    IsDeb = (istep>=debStartStep .and. istep<=debEndStep)
  end function IsDeb

  !---------------------------------------------------------------------------------------------------------------------

  function IsDebCell(ix, iy, iz)
    integer(4) :: ix, iy, iz, k
    logical :: IsDebCell, rc

    if(.not. IsDeb(-1)) then
      IsDebCell = .false.
      return
    end if

    if(iz<=0) then
      do k=1, ncz
        if(debCells(ix, iy, k)/=0) then
          IsDebCell = .true.
          return
        end if
      end do
      IsDebCell = .false.
    else
      IsDebCell = debCells(ix, iy, iz) /= 0
    end if
  end function IsDebCell

  !---------------------------------------------------------------------------------------------------------------------

  function IsDebSideX(ix, iy, iz)
    integer(4) :: ix, iy, iz, k
    logical :: IsDebSideX

    if(.not. IsDeb(-1)) then
      IsDebSideX = .false.
      return
    end if

    if(iz<=0) then
      do k=1, ncz
        if(debSidesX(ix, iy, k)/=0) then
          IsDebSideX = .true.
          return
        end if
      end do
      IsDebSideX = .false.
    else
      if(iy>ubound(debSidesX,2)) then
        debdum=0
      endif
      IsDebSideX = debSidesX(ix, iy, iz) /= 0
    end if
  end function IsDebSideX

  !---------------------------------------------------------------------------------------------------------------------

  function IsDebSideY(ix, iy, iz)
    integer(4) :: ix, iy, iz, k
    logical :: IsDebSideY

    if(.not. IsDeb(-1)) then
      IsDebSideY = .false.
      return
    end if

    if(iz<=0) then
      do k=1, ncz
        if(debSidesY(ix, iy, k)/=0) then
          IsDebSideY = .true.
          return
        end if
      end do
      IsDebSideY = .false.
    else
      IsDebSideY = debSidesY(ix, iy, iz) /= 0
    end if
  end function IsDebSideY

  !---------------------------------------------------------------------------------------------------------------------

  function IsDebSideZ(ix, iy, iz)
    integer(4) :: ix, iy, iz, k
    logical :: IsDebSideZ

    if(.not. IsDeb(-1)) then
      IsDebSideZ = .false.
      return
    end if

    if(iz<=0) then
      do k=1, nz
        if(debSidesZ(ix, iy, k)/=0) then
          IsDebSideZ = .true.
          return
        end if
      end do
      IsDebSideZ = .false.
    else
      IsDebSideZ = debSidesZ(ix, iy, iz) /= 0
    end if
  end function IsDebSideZ

  !---------------------------------------------------------------------------------------------------------------------

  subroutine avost(text, i)
    character(*), optional:: text
    integer, optional:: i

    if (present(text) .and. present(i)) then
      write(*,"(a,i0,3a,i0)") "step=", istep, ". AVOST: ", text, ":", i
    else if (present(text)) then
      write(*,"(a,i0,2a)") "step=", istep, ". AVOST: ", text
    else
      write(*,"(a,i0,2a)") "step=", istep, ". AVOST"
    end if
    stop
  end subroutine

  !---------------------------------------------------------------------------------------------------------------------

  subroutine assert(cond, text)
    logical:: cond
    character(*), optional:: text
    logical, save:: isFirst = .true.

    if(.not. cond) then
      if(isFirst) then
        if (present(text)) write(*,"(a,a)") "ASSERT! ", text
        isFirst = .false.
      end if
      if(optNoAssert) return
      if(nprint>0 .or. tprint>0.) then
        iout = 999999
#if 0
        call Gout
#endif
      end if
      call avost
    end if
  end subroutine

  !---------------------------------------------------------------------------------------------------------------------

  subroutine assertn(cond, text, i, j, k)
    logical cond
    character(*) text
    integer i, j, k
    logical, save:: isFirst = .true.

    if(.not. cond) then
      if(isFirst) then
        write(*,"(a,a,3I5)") "ASSERT! ", text, i, j, k
        isFirst = .false.
      end if
      if(optNoAssert) return
      if(nprint>0 .or. tprint>0.) then
        iout = 999999
#if 0
        call Gout
#endif
      end if
      call avost
    end if
  end subroutine

end module
