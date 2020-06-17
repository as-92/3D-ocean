! модуль расчёта 0-переменных по системе уравнений мелкой воды

subroutine ShWater
  use variables
  implicit none

  call PhaseSw1
  call PhaseSw2X
  call PhaseSw2Y
  call SumFluxesSw
  call PhaseSw3
  call CalcSwW

end subroutine ShWater

!------------------------------------------------------------------------------------------------------------------------------------

! Фаза 1 мелкой воды
subroutine PhaseSw1
  use variables
  implicit none

  integer(4):: i, j, ic
  real(8):: h_old,    h_new
  real(8):: m_old,    m_new
  real(8):: impx_old, impx_new
  real(8):: impy_old, impy_new

  do ic=1,ncells
    i=cells.ptr(ic).i
    j=cells.ptr(ic).j

    ! переменные на t=t[n]:
    h_old = c_h0(i,j)                                         ! глубина h для t=t[n]
    m_old = c_h0(i,j) * c_rho0(i,j)                           ! масса h*rho для t=t[n]
    impx_old = m_old * c_u0(i,j)                              ! импульс h*rho*u для t=t[n]
    impy_old = m_old * c_v0(i,j)                              ! импульс h*rho*v для t=t[n]

    ! переменные на t=t[n+1/2]:
    h_new = h_old - 0.5 * dt * c_sumFluxH0(i,j)               ! глубина h для t=t[n+1/2]
    m_new = m_old - 0.5 * dt * c_sumFluxM0(i,j)               ! масса h*rho для t=t[n+1/2]
    impx_new = impx_old - 0.5 * dt * c_sumFluxIx0(i,j)        ! импульс h*rho*u для t=t[n+1/2]
    impy_new = impy_old - 0.5 * dt * c_sumFluxIy0(i,j)        ! импульс h*rho*v для t=t[n+1/2]

    cs_h0(i,j) = h_new                                        ! глубина t=t[n+1/2]
    cs_rho0(i,j) = m_new / h_new                              ! плотность
    cs_u0(i,j) = impx_new / m_new                             ! компонента скорости по оси X
    cs_v0(i,j) = impy_new / m_new                             ! компонента скорости по оси Y

    cs_G0(i,j) = sqrt(g / h_new)

    if(IsDebCell(i, j, -1)) then
      write(17,"(/'Phase Sw1',2(';',i0),';',100(';',f0.5))") i, j, (x(i)+x(i+1))/2., (y(j)+y(j+1))/2.
      write(17,"(a)") ";i;j;;h;rho;u;v;fluxH;fluxM;fluxU;fluxV"
      write(17,"('side-L:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fx_h0(i  ,j  ), fx_rho0(i  ,j  ), fx_u0(i  ,j  ), fx_v0(i  ,j  ), &
            fx_fluxH0(i,j), fx_fluxM0(i,j), fx_fluxIx0(i,j), fx_fluxIy0(i,j)
      write(17,"('side-R:',2(';',i0),';',100(';',1p,g0.16))") &
          i+1, j, &
            fx_h0(i+1,j  ), fx_rho0(i+1,j  ), fx_u0(i+1,j  ), fx_v0(i+1,j  ), &
            fx_fluxH0(i+1,j), fx_fluxM0(i+1,j), fx_fluxIx0(i+1,j), fx_fluxIy0(i+1,j)
      write(17,"('side-B:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fy_h0(i  ,j  ), fy_rho0(i  ,j  ), fy_u0(i  ,j  ), fy_v0(i  ,j  ), &
            fy_fluxH0(i,j), fy_fluxM0(i,j), fy_fluxIx0(i,j), fy_fluxIy0(i,j)
      write(17,"('side-T:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j+1, &
            fy_h0(i  ,j+1), fy_rho0(i  ,j+1), fy_u0(i  ,j+1), fy_v0(i  ,j+1), &
            fy_fluxH0(i,j+1), fy_fluxM0(i,j+1), fy_fluxIx0(i,j+1), fy_fluxIy0(i,j+1)
      write(17,"('SumFlux:' ,2(';',i0),';;;;;',100(';',1p,g0.16))") i, j, &
            c_sumFluxH0(i,j), c_sumFluxM0(i,j), c_sumFluxIx0(i,j), c_sumFluxIy0(i,j)
      write(17,"('UC:' ,2(';',i0),';',100(';',1p,g0.16))") i, j, c_h0(i,j),  c_rho0(i,j),  c_u0(i,j),  c_v0(i,j)
      write(17,"('UCS:',2(';',i0),';',100(';',1p,g0.16))") i, j, cs_h0(i,j), cs_rho0(i,j), cs_u0(i,j), cs_v0(i,j), cs_G0(i,j)
    end if

  end do

  cs_z0(:,:) = c_b(:,:) + cs_h0(:,:)

end subroutine PhaseSw1

!------------------------------------------------------------------------------------------------------------------------------------

! Фаза 3 мелкой воды
subroutine PhaseSw3
  use variables
  implicit none

  integer(4):: i, j, ic
  real(8):: h_old,    h_new
  real(8):: m_old,    m_new
  real(8):: impx_old, impx_new
  real(8):: impy_old, impy_new

  do ic=1,ncells
    i=cells.ptr(ic).i
    j=cells.ptr(ic).j

    ! переменные на t=t[n]:
    h_old = cs_h0(i,j)                                        ! глубина h для t=t[n]
    m_old = cs_h0(i,j) * cs_rho0(i,j)                         ! масса h*rho для t=t[n]
    impx_old = m_old * cs_u0(i,j)                             ! импульс h*rho*u для t=t[n]
    impy_old = m_old * cs_v0(i,j)                             ! импульс h*rho*v для t=t[n]

    ! переменные на t=t[n+1/2]:
    h_new = h_old - 0.5 * dt * c_sumFluxH0(i,j)               ! глубина h для t=t[n+1/2]
    m_new = m_old - 0.5 * dt * c_sumFluxM0(i,j)               ! масса h*rho для t=t[n+1/2]
    impx_new = impx_old - 0.5 * dt * c_sumFluxIx0(i,j)        ! импульс h*rho*u для t=t[n+1/2]
    impy_new = impy_old - 0.5 * dt * c_sumFluxIy0(i,j)        ! импульс h*rho*v для t=t[n+1/2]

    cn_h0(i,j) = h_new                                        ! глубина t=t[n+1/2]
    cn_rho0(i,j) = m_new / h_new                              ! плотность
    cn_u0(i,j) = impx_new / m_new                             ! компонента скорости по оси X
    cn_v0(i,j) = impy_new / m_new                             ! компонента скорости по оси Y

    if(IsDebCell(i, j, -1)) then
      write(17,"(/'Phase Sw3',2(';',i0),';',100(';',f0.5))") i, j, (x(i)+x(i+1))/2., (y(j)+y(j+1))/2.
      write(17,"(a)") ";i;j;;h;rho;u;v;fluxH;fluxM;fluxU;fluxV"
      write(17,"('side-XL:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fxn_h0(i  ,j  ), fxn_rho0(i  ,j  ), fxn_u0(i  ,j  ), fxn_v0(i  ,j  ), &
            fx_fluxH0(i,j), fx_fluxM0(i,j), fx_fluxIx0(i,j), fx_fluxIy0(i,j)
      write(17,"('side-XR:',2(';',i0),';',100(';',1p,g0.16))") &
          i+1, j, &
            fxn_h0(i+1,j  ), fxn_rho0(i+1,j  ), fxn_u0(i+1,j  ), fxn_v0(i+1,j  ), &
            fx_fluxH0(i+1,j), fx_fluxM0(i+1,j), fx_fluxIx0(i+1,j), fx_fluxIy0(i+1,j)
      write(17,"('side-YL:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fyn_h0(i  ,j  ), fyn_rho0(i  ,j  ), fyn_u0(i  ,j  ), fyn_v0(i  ,j  ), &
            fy_fluxH0(i,j), fy_fluxM0(i,j), fy_fluxIx0(i,j), fy_fluxIy0(i,j)
      write(17,"('side-YR:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j+1, &
            fyn_h0(i  ,j+1), fyn_rho0(i  ,j+1), fyn_u0(i  ,j+1), fyn_v0(i  ,j+1), &
            fy_fluxH0(i,j+1), fy_fluxM0(i,j+1), fy_fluxIx0(i,j+1), fy_fluxIy0(i,j+1)
      write(17,"('SumFlux:' ,2(';',i0),';;;;;',100(';',1p,g0.16))") i, j, &
            c_sumFluxH0(i,j), c_sumFluxM0(i,j), c_sumFluxIx0(i,j), c_sumFluxIy0(i,j)
      write(17,"('UCS:',2(';',i0),';',100(';',1p,g0.16))") i, j, cs_h0(i,j), cs_rho0(i,j), cs_u0(i,j), cs_v0(i,j)
      write(17,"('UCN:',2(';',i0),';',100(';',1p,g0.16))") i, j, cn_h0(i,j), cn_rho0(i,j), cn_u0(i,j), cn_v0(i,j)
    end if

  end do

  cn_z0(:,:) = c_b(:,:) + cn_h0(:,:)

end subroutine PhaseSw3

!------------------------------------------------------------------------------------------------------------------------------------

! Фаза 2 мелкой воды на X-гранях
subroutine PhaseSw2X
  use variables
  implicit none

  interface
    subroutine TransportSwInvX(itype, i, j, off, inv); integer(4):: itype, i, j, off; real(R8):: inv; end subroutine TransportSwInvX
    function GetCellInvSwX(type, i, j, isN); integer(4):: type, i, j; logical:: isN; real(8):: GetCellInvSwX; end
    subroutine BuildFluxesSwX(i, j); integer(4):: i, j; end subroutine
  end interface

  integer(4):: i, j, il, ir, ic, is, ns, ibc
  real(R8):: invR, invQ, invU, invD                             ! новые инварианты
  real(R8):: soundr, soundl, sound, uf, M, Gr, Gq
  real(R8):: hi, ui, vi, rhoi
  real(R8):: bcH, bcU, bcV, bcRho
  logical:: isOwnRight                                          ! признак "ячейка справа"

  ! внутренние грани:
  ns = nfx_sides(BC_INNER)                                      ! количество внутренних граней
  do is=1, ns
    i = fx_sides(BC_INNER).ptr(is).i                            ! индексы грани на сетке
    j = fx_sides(BC_INNER).ptr(is).j
    call assertn(i>1 .and. i<nxx, "PhaseSW2_X. Ошибка сетки (BC_INNER)", i, j, -1)

    if(IsDebSideX(i, j, -1)) then
      write(17,"(/'Phase Sw2-X',2(';',i0),';',';inner')") i, j
    endif

    il = i-1
    ir = i

    ! вычисляем Мах на грани:
    soundr = sqrt(g * cs_h0(ir,j))
    soundl = sqrt(g * cs_h0(il,j))
    sound = (soundl + soundr) / 2.
    uf=0.5*(cs_u0(ir,j) + cs_u0(il,j))                          ! скорость на грани
    M = uf / sound                                              ! Мах на грани

    if(abs(M)>0.9) then
      write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
      call avost
    endif

    ! инвариант R:
    call TransportSwInvX(T_INVR, i, j, DIRP, invR)              ! перенос через левую ячейку вправо
    Gr = cs_G0(il,j)

    ! инвариант Q:
    call TransportSwInvX(T_INVQ, i, j, DIRM, invQ)              ! перенос через правую ячейку влево
    Gq = -cs_G0(ir,j)

    ! инварианты U и Rho:
    if(M>eps) then
      call TransportSwInvX(T_INVU, i, j, DIRP, invU)      ! перенос через левую ячейку вправо
      call TransportSwInvX(T_INVD, i, j, DIRP, invD)      ! перенос через левую ячейку вправо
    else if(M<-eps) then
      call TransportSwInvX(T_INVU, i, j, DIRM, invU)      ! перенос через правую ячейку влево
      call TransportSwInvX(T_INVD, i, j, DIRM, invD)      ! перенос через правую ячейку влево
    else
      invU = (cs_v0(il,j) + cs_v0(ir,j)) / 2.
      invD = (cs_rho0(il,j) + cs_rho0(ir,j)) / 2.
    end if

    ! восстанавливаем значения из инвариантов:
    rhoi = invD                                        ! плотность
    vi = invU                                          ! тангенцальная скорость

    ! invR = u + Gr * (b + h)
    ! invQ = u + Gq * (b + h)
    !----------------------------------
    ! invR - invQ = (Gr - Gq)*(b + h)
    ! invR*Gq - invQ*Gr = (Gq - Gr)*u

    hi = (invR - invQ) / (Gr - Gq) - fx_b(i,j)         ! глубина
    ui = (invR*Gq - invQ*Gr) / (Gq - Gr)               ! нормальная скорость

    fxn_h0(i,j) = hi
    fxn_u0(i,j) = ui
    fxn_v0(i,j) = vi                                         ! тангенцальная скорость
    fxn_rho0(i,j) = rhoi                                        ! плотность

    if(IsDebSideX(i, j, -1)) then
      write(17,"(a)") ";i;j;;h;rho;u;v;R;Q;U;D;Gr;Gq"
      write(17,"('side:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fx_h0(i,j), fx_rho0(i,j), fx_u0(i,j), fx_v0(i,j), &
            invR, invQ, invU, invD, Gr, Gq
      write(17,"('UFN:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fxn_h0(i,j), fxn_rho0(i,j), fxn_u0(i,j), fxn_v0(i,j)
    end if

    call BuildFluxesSwX(i,j)                                      ! построение потоков на грани

  end do

  !.. граничные условия ...................................................................

  ! периодические ГУ:
  ns = nfx_sides(BC_PERIODIC)                                   ! количество граней
  do is=1, ns
    i = fx_sides(BC_PERIODIC).ptr(is).i                         ! индексы грани на сетке
    j = fx_sides(BC_PERIODIC).ptr(is).j
    call assertn(i==1 .or. i==nxx, "Phase2_X, периодические ГУ", i,j,-1)

    if(i/=1) cycle                                              ! только для левой границы

    il = ncx                                                    ! ячейка слева (крайняя правая)
    ir = 1                                                      ! ячейка справа

    ! вычисляем Мах на грани:
    soundr = sqrt(g * cs_h0(ir,j))
    soundl = sqrt(g * cs_h0(il,j))
    sound = (soundl + soundr) / 2.
    uf=0.5*(cs_u0(ir,j) + cs_u0(il,j))                          ! скорость на грани
    M = uf / sound                                              ! Мах на грани

    if(abs(M)>0.9) then
      write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
      call avost
    endif

    ! инвариант R:
    call TransportSwInvX(T_INVR, nxx, j, DIRP, invR)        ! перенос через левую ячейку вправо
    Gr = cs_G0(il,j)

    ! инвариант Q:
    call TransportSwInvX(T_INVQ, 1, j, DIRM, invQ)          ! перенос через правую ячейку влево
    Gq = -cs_G0(ir,j)

    ! инварианты U и D:
    if(M>eps) then
      call TransportSwInvX(T_INVU, nxx, j, DIRP, invU)        ! перенос через левую ячейку вправо
      call TransportSwInvX(T_INVD, nxx, j, DIRP, invD)        ! перенос через левую ячейку вправо
    else if(M<-eps) then
      call TransportSwInvX(T_INVU, 1, j, DIRM, invU)          ! перенос через правую ячейку влево
      call TransportSwInvX(T_INVD, 1, j, DIRM, invD)          ! перенос через правую ячейку влево
    else
      invU = (cs_v0(il,j) + cs_v0(ir,j)) / 2.
      invD = (cs_rho0(il,j) + cs_rho0(ir,j)) / 2.
    end if

    ! восстанавливаем значения из инвариантов:
    fxn_rho0(1,j) = invD
    fxn_rho0(nxx,j) = fxn_rho0(1,j)
    fxn_v0(1,j) = invU
    fxn_v0(nxx,j) = fxn_v0(1,j)

    ! invR = u + Gr * (b + h)
    ! invQ = u + Gq * (h + h)
    !----------------------------------
    ! invR - Q = (Gr - Gq)*(b + h)
    ! invR*Gq - invQ*Gr = (Gq - Gr)*u

    fxn_h0(1,j) = (invR - invQ) / (Gr - Gq) - fx_b(i,j)
    fxn_h0(nxx,j) = fxn_h0(1,j)
    fxn_u0(1,j) = (invR*Gq - invQ*Gr) / (Gq - Gr)
    fxn_u0(nxx,j) = fxn_u0(1,j)

    call BuildFluxesSwX(1,j)                                    ! построение потоков на грани
    call BuildFluxesSwX(nxx,j)

    if(IsDebSideX(i, j, -1)) then
      write(17,"(a)") ""
      write(17,"('Phase Sw2-X',2(';',i0))") i, j
      write(17,"(a)") "BC-period;i;j;;h;rho;u;v;R;Q;U;D;Gr;Gq"
      write(17,"('side:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fx_h0(i,j), fx_rho0(i,j), fx_u0(i,j), fx_v0(i,j), &
            invR, invQ, invU, invD, Gr, Gq
    end if

  end do  ! do is

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! ГУ "стенка с проскальзыванием":
  ns = nfx_sides(BC_SLIDE)                                      ! количество граней
  do is=1, ns
    i = fx_sides(BC_SLIDE).ptr(is).i                            ! индексы грани на сетке
    j = fx_sides(BC_SLIDE).ptr(is).j

    if(IsDebSideX(i, j, -1)) then
      debdum = 0
    endif

    ! с какой стороны от стенки расчетная область?
    if(i==1) then
      isOwnRight = .true.                                       ! справа, если грань на левой стенке
    else if(i==nxx) then
      isOwnRight = .false.
    else
      isOwnRight = (c_type(i,j)>=0)
    end if

    if(isOwnRight) then                                         ! живая ячейка справа

      ! вычисляем Мах на грани:
      soundr = sqrt(g * cs_h0(i,j))
      sound = soundr
      uf = cs_u0(i,j)
      M = uf / sound                                            ! Мах на грани

      if(abs(M)>0.9) then
        write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
        call avost
      endif

      ! инвариант Q:
      call TransportSwInvX(T_INVQ, i, j, DIRM, invQ)            ! перенос через правую ячейку влево
      Gq = -cs_G0(ir,j)

      ! инварианты U и T:
      if(M<-eps) then
        call TransportSwInvX(T_INVU, i, j, DIRM, invU)        ! перенос через правую ячейку влево
        call TransportSwInvX(T_INVD, i, j, DIRM, invD)
      else
        invU = GetCellInvSwX(T_INVU, i, j, .true.)                ! .. вычисляем инвариант в центре ячейки
        invD = GetCellInvSwX(T_INVD, i, j, .true.)
      end if

      ! восстанавливаем значения из инвариантов:
      fxn_rho0(i,j) = invD
      fxn_v0(i,j) = invU

      ! u = 0
      ! invQ = 0 + Gq * (b + h)
      !----------------------------------
      ! u = 0
      ! h = invQ / Gq - b

      fxn_h0(i,j) = invQ / Gq - fx_b(i,j)
      fxn_u0(i,j) = 0.

    else    ! if(c_type(i,j)>=0)                                ! живая ячейка слева, стенка справа

      ! вычисляем число Маха на грани:
      il = i - 1
      sound = sqrt(g * cs_h0(il,j))
      uf = cs_u0(il,j)                                          ! скорость на грани
      M = uf / sound                                            ! Мах на грани

      if(abs(M)>0.9) then
        write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
        call avost
      endif

      ! инвариант R:
      call TransportSwInvX(T_INVR, i, j, DIRP, invR)          ! перенос через левую ячейку вправо
      Gr = cs_G0(il,j)

      ! инварианты U и T:
      if(M>eps) then                                            ! течение на стенку   >O|
        call TransportSwInvX(T_INVU, i, j, DIRP, invU)        ! .. перенос вправо
        call TransportSwInvX(T_INVD, i, j, DIRP, invD)
      else                                                      ! течение от стенки   <O|
        invU = GetCellInvSwX(T_INVU, il, j, .true.)               ! .. вычисляем инвариант в центре ячейки
        invD = GetCellInvSwX(T_INVD, il, j, .true.)
      end if

      ! восстанавливаем значения из инвариантов:
      fxn_rho0(i,j) = invD
      fxn_v0(i,j) = invU

      ! u = 0
      ! invR = 0 + Gr * (b + h)
      !----------------------------------
      ! u = 0
      ! h = invR / Gr - b

      fxn_h0(i,j) = invR / Gr - fx_b(i,j)
      fxn_u0(i,j) = 0.

    end if  ! if(c_type(i,j)>=0)

    call BuildFluxesSwX(i,j)                                      ! построение потоков на грани

  end do  ! do is

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! ГУ "стенка с прилипанием":
  ns = nfx_sides(BC_STICK)                                      ! количество граней
  do is=1, ns
    i = fx_sides(BC_STICK).ptr(is).i                            ! индексы грани на сетке
    j = fx_sides(BC_STICK).ptr(is).j

    if(IsDebSideX(i, j, -1)) then
      debdum = 0
    endif

    ! с какой стороны от стенки расчетная область?
    if(i==1) then
      isOwnRight = .true.                                       ! справа, если грань на левой стенке
    else if(i==nxx) then
      isOwnRight = .false.
    else
      isOwnRight = (c_type(i,j)>=0)
    end if

    if(isOwnRight) then                                         ! живая ячейка справа

      ! вычисляем Мах на грани:
      soundr = sqrt(g * cs_h0(i,j))
      sound = soundr
      uf = cs_u0(i,j)
      M = uf / sound                                            ! Мах на грани

      if(abs(M)>0.9) then
        write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
        call avost
      endif

      ! инвариант Q:
      call TransportSwInvX(T_INVQ, i, j, DIRM, invQ)            ! перенос через правую ячейку влево
      Gq = -cs_G0(i, j)

      ! инвариант D:
      if(M<-eps) then                                           ! течение на стенку
        call TransportSwInvX(T_INVD, i, j, DIRM, invD)          ! .. перенос инварианта чараз ячейку
      else                                                      ! течение от стенки
        invD = GetCellInvSwX(T_INVD, i, j, .true.)              ! .. вычисляем инвариант в центре ячейки
      end if

      ! восстанавливаем значения из инвариантов:
      fxn_rho0(i,j) = invD
      fxn_v0(i,j) = 0.

      ! u = 0
      ! invQ = 0 + Gq * (b + h)
      !----------------------------------
      ! u = 0
      ! h = invQ / Gq - b

      fxn_h0(i,j) = invQ / Gq - fx_b(i,j)
      fxn_u0(i,j) = 0.

    else    ! if(c_type(i,j)>=0)                                ! живая ячейка слева, стенка справа

      ! вычисляем число Маха на грани:
      il = i - 1
      sound = sqrt(g * cs_h0(il,j))
      uf = cs_u0(il,j)                                         ! скорость на грани
      M = uf / sound                                            ! Мах на грани

      if(abs(M)>0.9) then
        write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
        call avost
      endif

      ! инвариант R:
      call TransportSwInvX(T_INVR, i, j, DIRP, invR)          ! перенос через левую ячейку вправо
      Gr = cs_G0(il, j)

      ! инвариант D:
      if(M>eps) then                                            ! течение на стенку
        call TransportSwInvX(T_INVD, i, j, DIRP, invD)        ! .. перенос инварианта чараз ячейку
      else                                                      ! течение от стенки
        invD = GetCellInvSwX(T_INVD, il, j, .true.)                ! .. вычисляем инвариант в центре ячейки
      end if

      ! восстанавливаем значения из инвариантов:
      fxn_rho0(i,j) = invD
      fxn_v0(i,j) = 0.

      ! u = 0
      ! invR = 0 + Gr * (b + h)
      !----------------------------------
      ! u = 0
      ! h = invR / Gr - b

      fxn_h0(i,j) = invR / Gr - fx_b(i,j)
      fxn_u0(i,j) = 0.

    end if  ! if(c_type(i,j)>=0)

    call BuildFluxesSwX(i,j)                                    ! построение потоков на грани

  end do  ! do is

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! ГУ "постоянный вток":
  ns = nfx_sides(BC_IN)                                         ! количество граней
  do is=1, ns
    i = fx_sides(BC_IN).ptr(is).i                             ! индексы грани на сетке
    j = fx_sides(BC_IN).ptr(is).j
    il = i - 1
    ir = i

    ibc = fx_bc(i,j)                                            ! индекс индивидуальных параметров ГУ

    ! заграничные данные:
    bcH = bcData(ibc).h
    bcU = bcData(ibc).u
    bcV = bcData(ibc).v
    bcRho = bcData(ibc).rho

    if(IsDebSideX(i, j, -1)) then
      write(17,"(/'Phase Sw2-X',2(';',i0),';',';BC_IN')") i, j
    endif

    isOwnRight = fx_IsOwnRight(i,j)

    ic = i - 1; if(isOwnRight) ic = i

    ! вычисляем Мах на грани:
    soundr = sqrt(g * bcH)
    soundl = sqrt(g * cs_h0(ic,j))
    sound = (soundl + soundr) / 2.
    uf = 0.5*(cs_u0(ic,j) + bcU)                                ! скорость на грани
    M = uf / sound                                              ! Мах на грани

    if(abs(M)>0.9) then
      write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
      call avost
    endif

    ! инвариант R (слева направо):
    if(isOwnRight) then                                         ! ячейка справа, граница слева
      Gr = sqrt(g / bcH)
      invR = bcU + Gr * bcH
    else                                                        ! слева есть ячейка
      call TransportSwInvX(T_INVR, i, j, DIRP, invR)            ! .. перенос через левую ячейку вправо
      Gr = cs_G0(ic,j)
    end if

    ! инвариант Q (справа налево):
    if(isOwnRight) then                                         ! справа есть ячейка
      call TransportSwInvX(T_INVQ, i, j, DIRM, invQ)            ! перенос через правую ячейку влево
      Gq = -cs_G0(ir,j)
    else
      Gq = -sqrt(g / bcH)
      invQ = bcU + Gq * bcH
    end if

    ! инварианты U и Rho:
      ! транспортные инварианты:
    if(isOwnRight .and. M<-eps) then                            ! течение из ячейки на границу справа налево  |<O
      call TransportSwInvX(T_INVU, i, j, DIRM, invU)            ! .. перенос через правую ячейку влево
      call TransportSwInvX(T_INVD, i, j, DIRM, invD)            ! .. перенос через правую ячейку влево
    else if(.not.isOwnRight .and. M>eps) then                   ! течение из ячейки на границу слева направо  O>|
      call TransportSwInvX(T_INVU, i, j, DIRP, invU)            ! .. перенос через левую ячейку вправо
      call TransportSwInvX(T_INVD, i, j, DIRP, invD)            ! .. перенос через левую ячейку вправо
    else                                                        ! течение из-за границы
      invU = bcV
      invD = bcRho
    end if

    ! восстанавливаем значения из инвариантов:
    rhoi = invD                                                 ! плотность
    vi = invU                                                   ! тангенцальная скорость

    ! invR = u + Gr * (b + h)
    ! invQ = u + Gq * (b + h)
    !----------------------------------
    ! invR - invQ = (Gr - Gq)*(b + h)
    ! invR*Gq - invQ*Gr = (Gq - Gr)*u

    hi = (invR - invQ) / (Gr - Gq) - fx_b(i,j)                  ! глубина
    ui = (invR*Gq - invQ*Gr) / (Gq - Gr)                        ! нормальная скорость

    fxn_h0(i,j) = hi
    fxn_u0(i,j) = ui
    fxn_v0(i,j) = vi                                            ! тангенцальная скорость
    fxn_rho0(i,j) = rhoi                                        ! плотность

    if(IsDebSideX(i, j, -1)) then
      write(17,"(a)") ";i;j;;h;rho;u;v;R;Q;U;D;Gr;Gq"
      write(17,"('side:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            bcH, bcRho, bcU, bcV
      write(17,"('UF:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fx_h0(i,j), fx_rho0(i,j), fx_u0(i,j), fx_v0(i,j), &
            invR, invQ, invU, invD, Gr, Gq
      write(17,"('UFN:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fxn_h0(i,j), fxn_rho0(i,j), fxn_u0(i,j), fxn_v0(i,j)
    end if

    call BuildFluxesSwX(i,j)                                    ! построение потоков на грани

  end do  ! do is

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! ГУ "переменный вток":
  ns = nfx_sides(BC_IN_T)                                       ! количество граней
  do is=1, ns
    i = fx_sides(BC_IN_T).ptr(is).i                             ! индексы грани на сетке
    j = fx_sides(BC_IN_T).ptr(is).j
    il = i - 1
    ir = i

    ibc = fx_bc(i,j)                                            ! индекс индивидуальных параметров ГУ

    ! заграничные данные:
    bcH = bcData(ibc).h
    bcU = bcData(ibc).u
    bcV = bcData(ibc).v
    bcRho = bcData(ibc).rho

    if(IsDebSideX(i, j, -1)) then
      write(17,"(/'Phase Sw2-X',2(';',i0),';',';BC_IN_T')") i, j
    endif

    isOwnRight = fx_IsOwnRight(i,j)

    ic = i - 1; if(isOwnRight) ic = i

    ! вычисляем Мах на грани:
    soundr = sqrt(g * bcH)
    soundl = sqrt(g * cs_h0(ic,j))
    sound = (soundl + soundr) / 2.
    uf = 0.5*(cs_u0(ic,j) + bcU)                                ! скорость на грани
    M = uf / sound                                              ! Мах на грани

    if(abs(M)>0.9) then
      write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
      call avost
    endif

    ! инвариант R (слева направо):
    if(isOwnRight) then                                         ! ячейка справа, граница слева
      Gr = sqrt(g / bcH)
      invR = bcU + Gr * bcH
    else                                                        ! слева есть ячейка
      call TransportSwInvX(T_INVR, i, j, DIRP, invR)            ! .. перенос через левую ячейку вправо
      Gr = cs_G0(ic,j)
    end if

    ! инвариант Q (справа налево):
    if(isOwnRight) then                                         ! справа есть ячейка
      call TransportSwInvX(T_INVQ, i, j, DIRM, invQ)            ! перенос через правую ячейку влево
      Gq = -cs_G0(ir,j)
    else
      Gq = -sqrt(g / bcH)
      invQ = bcU + Gq * bcH
    end if

    ! инварианты U и Rho:
      ! транспортные инварианты:
    if(isOwnRight .and. M<-eps) then                            ! течение из ячейки на границу справа налево  |<O
      call TransportSwInvX(T_INVU, i, j, DIRM, invU)            ! .. перенос через правую ячейку влево
      call TransportSwInvX(T_INVD, i, j, DIRM, invD)            ! .. перенос через правую ячейку влево
    else if(.not.isOwnRight .and. M>eps) then                   ! течение из ячейки на границу слева направо  O>|
      call TransportSwInvX(T_INVU, i, j, DIRP, invU)            ! .. перенос через левую ячейку вправо
      call TransportSwInvX(T_INVD, i, j, DIRP, invD)            ! .. перенос через левую ячейку вправо
    else                                                        ! течение из-за границы
      invU = bcV
      invD = bcRho
    end if

    ! восстанавливаем значения из инвариантов:
    rhoi = invD                                                 ! плотность
    vi = invU                                                   ! тангенцальная скорость

    ! invR = u + Gr * (b + h)
    ! invQ = u + Gq * (b + h)
    !----------------------------------
    ! invR - invQ = (Gr - Gq)*(b + h)
    ! invR*Gq - invQ*Gr = (Gq - Gr)*u

    hi = (invR - invQ) / (Gr - Gq) - fx_b(i,j)                  ! глубина
    ui = (invR*Gq - invQ*Gr) / (Gq - Gr)                        ! нормальная скорость

    fxn_h0(i,j) = hi
    fxn_u0(i,j) = ui
    fxn_v0(i,j) = vi                                            ! тангенцальная скорость
    fxn_rho0(i,j) = rhoi                                        ! плотность

    if(IsDebSideX(i, j, -1)) then
      write(17,"(a)") ";i;j;;h;rho;u;v;R;Q;U;D;Gr;Gq"
      write(17,"('side:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            bcH, bcRho, bcU, bcV
      write(17,"('UF:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fx_h0(i,j), fx_rho0(i,j), fx_u0(i,j), fx_v0(i,j), &
            invR, invQ, invU, invD, Gr, Gq
      write(17,"('UFN:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fxn_h0(i,j), fxn_rho0(i,j), fxn_u0(i,j), fxn_v0(i,j)
    end if

    call BuildFluxesSwX(i,j)                                    ! построение потоков на грани

  end do  ! do is

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! ГУ "фиксированная скорость":
  ns = nfx_sides(BC_FIX_U)                                       ! количество граней
  do is=1, ns
    i = fx_sides(BC_FIX_U).ptr(is).i                             ! индексы грани на сетке
    j = fx_sides(BC_FIX_U).ptr(is).j
    il = i - 1
    ir = i
    call assertn(.false., "Phase2_X. Не реализовано BC_FIX_U", -1, -1, -1)

    call BuildFluxesSwX(i,j)                                    ! построение потоков на грани

  end do  ! do is

end subroutine PhaseSw2X

!------------------------------------------------------------------------------------------------------------------------------------

! Фаза 2 мелкой воды на Y-гранях
subroutine PhaseSw2Y
  use variables
  implicit none

  interface
    subroutine TransportSwInvY(itype, i, j, off, inv); integer(4):: itype, i, j, off; real(R8):: inv; end subroutine
    function GetCellInvSwY(type, i, j, isN); integer(4):: type, i, j; logical:: isN; real(8):: GetCellInvSwY; end
    subroutine BuildFluxesSwY(i, j); integer(4):: i, j; end subroutine
  end interface

  integer(4):: i, j, jl, jr, jc, is, ns, ibc
  real(R8):: invR, invQ, invU, invD                             ! новые инварианты
  real(R8):: soundr, soundl, sound, uf, M, Gr, Gq
  real(R8):: hi, ui, vi, rhoi
  real(R8):: bcH, bcU, bcV, bcRho
  logical:: isOwnRight                                          ! признак "ячейка справа"

  ! внутренние грани:
  ns = nfy_sides(BC_INNER)                                      ! количество внутренних граней
  do is=1, ns
    i = fy_sides(BC_INNER).ptr(is).i                            ! индексы грани на сетке
    j = fy_sides(BC_INNER).ptr(is).j
    call assertn(j>1 .and. j<nyy, "PhaseSW2_Y. Ошибка сетки (BC_INNER)", i, j, -1)

    if(IsDebSideY(i, j, -1)) then
      write(17,"(/'Phase Sw2-Y',2(';',i0),';;inner')") i, j
    endif

    jl = j-1
    jr = j

    ! вычисляем Мах на грани:
    soundr = sqrt(g * cs_h0(i,jr))
    soundl = sqrt(g * cs_h0(i,jl))
    sound = (soundl + soundr) / 2.
    uf=0.5*(cs_v0(i,jr) + cs_v0(i,jl))                          ! скорость на грани
    M = uf / sound                                              ! Мах на грани

    if(abs(M)>0.9) then
      write(*,*) "ERROR in PhaseSw2Y: |M|>1. side (i,j):", i, j
      call avost
    endif

    ! инвариант R:
    call TransportSwInvY(T_INVR, i, j, DIRP, invR)          ! перенос через левую ячейку вправо
    Gr = cs_G0(i,jl)

    ! инвариант Q:
    call TransportSwInvY(T_INVQ, i, j, DIRM, invQ)          ! перенос через правую ячейку влево
    Gq = -cs_G0(i,jr)

    ! инварианты U и Rho:
    if(M>eps) then
      call TransportSwInvY(T_INVU, i, j, DIRP, invU)      ! перенос через левую ячейку вправо
      call TransportSwInvY(T_INVD, i, j, DIRP, invD)      ! перенос через левую ячейку вправо
    else if(M<-eps) then
      call TransportSwInvY(T_INVU, i, j, DIRM, invU)      ! перенос через правую ячейку влево
      call TransportSwInvY(T_INVD, i, j, DIRM, invD)      ! перенос через правую ячейку влево
    else
      invU = (cs_u0(i,jl) + cs_u0(i,jr)) / 2.
      invD = (cs_rho0(i,jl) + cs_rho0(i,jr)) / 2.
    end if

    ! восстанавливаем значения из инвариантов:
    rhoi = invD                                        ! плотность
    ui = invU                                          ! тангенцальная скорость

    ! invR = u + Gr * (b + h)
    ! invQ = u + Gq * (b + h)
    !----------------------------------
    ! invR - invQ = (Gr - Gq)*(b + h)
    ! invR*Gq - invQ*Gr = (Gq - Gr)*u

    hi = (invR - invQ) / (Gr - Gq) - fy_b(i,j)         ! глубина
    vi = (invR*Gq - invQ*Gr) / (Gq - Gr)               ! нормальная скорость

    fyn_h0(i,j) = hi
    fyn_u0(i,j) = ui                                          ! тангенцальная скорость
    fyn_v0(i,j) = vi
    fyn_rho0(i,j) = rhoi                                        ! плотность

    if(IsDebSideY(i, j, -1)) then
      write(17,"(a)") ";i;j;;h;rho;u;v;R;Q;U;D;Gr;Gq"
      write(17,"('side:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fy_h0(i,j), fy_rho0(i,j), fy_u0(i,j), fy_v0(i,j), &
            invR, invQ, invU, invD, Gr, Gq
      write(17,"('UFN:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fyn_h0(i,j), fyn_rho0(i,j), fyn_u0(i,j), fyn_v0(i,j)
    end if

    call BuildFluxesSwY(i,j)                                      ! построение потоков на грани

  end do

  !.. граничные условия ...................................................................

  ! периодические ГУ:
  ns = nfy_sides(BC_PERIODIC)                                   ! количество граней
  do is=1, ns
    i = fy_sides(BC_PERIODIC).ptr(is).i                         ! индексы грани на сетке
    j = fy_sides(BC_PERIODIC).ptr(is).j
    call assertn(j==1 .or. j==nyy, "Phase2_Y, периодические ГУ", i,j,-1)

    if(j/=1) cycle                                              ! только для левой границы

    jl = ncy                                                    ! ячейка слева (крайняя правая)
    jr = 1                                                      ! ячейка справа

    ! вычисляем Мах на грани:
    soundr = sqrt(g * cs_h0(i,jr))
    soundl = sqrt(g * cs_h0(i,jl))
    sound = (soundl + soundr) / 2.
    uf=0.5*(cs_v0(i,jr) + cs_v0(i,jl))                          ! скорость на грани
    M = uf / sound                                              ! Мах на грани

    if(abs(M)>0.9) then
      write(*,*) "ERROR in PhaseSw2Y: |M|>1. side (i,j):", i, j
      call avost
    endif

    ! инвариант R:
    call TransportSwInvY(T_INVR, i, nyy, DIRP, invR)        ! перенос через левую ячейку вправо
    Gr = cs_G0(i,jl)

    ! инвариант Q:
    call TransportSwInvY(T_INVQ, i, 1, DIRM, invQ)          ! перенос через правую ячейку влево
    Gq = -cs_G0(i,jr)

    ! инварианты U и D:
    if(M>eps) then
      call TransportSwInvY(T_INVU, i, nyy, DIRP, invU)        ! перенос через левую ячейку вправо
      call TransportSwInvY(T_INVD, i, nyy, DIRP, invD)        ! перенос через левую ячейку вправо
    else if(M<-eps) then
      call TransportSwInvY(T_INVU, i, 1, DIRM, invU)          ! перенос через правую ячейку влево
      call TransportSwInvY(T_INVD, i, 1, DIRM, invD)          ! перенос через правую ячейку влево
    else
      invU = (cs_u0(i,jl) + cs_u0(i,jr)) / 2.
      invD = (cs_rho0(i,jl) + cs_rho0(i,jr)) / 2.
    end if

    ! восстанавливаем значения из инвариантов:
    fyn_rho0(i,1) = invD
    fyn_rho0(i,nyy) = fyn_rho0(i,1)
    fyn_u0(i,1) = invU
    fyn_u0(i,nyy) = fyn_u0(i,1)

    ! invR = u + Gr * (b + h)
    ! invQ = u + Gq * (b + h)
    !----------------------------------
    ! invR - Q = (Gr - Gq) * (b + h)
    ! invR*Gq - invQ*Gr = (Gq - Gr)*u

    fyn_h0(i,1) = (invR - invQ) / (Gr - Gq) - fy_b(i,j)
    fyn_h0(i,nyy) = fyn_h0(i,1)
    fyn_v0(i,1) = (invR*Gq - invQ*Gr) / (Gq - Gr)
    fyn_v0(i,nyy) = fyn_v0(i,1)

    call BuildFluxesSwY(i,1)                                      ! построение потоков на грани
    call BuildFluxesSwY(i,nyy)

  end do  ! do is

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! ГУ "стенка с проскальзыванием":
  ns = nfy_sides(BC_SLIDE)                                      ! количество граней
  do is=1, ns
    i = fy_sides(BC_SLIDE).ptr(is).i                            ! индексы грани на сетке
    j = fy_sides(BC_SLIDE).ptr(is).j

    if(IsDebSideY(i, j, -1)) then
      write(17,"(/'Phase Sw2-Y',2(';',i0),';;BC-slide')") i, j
    endif

    ! с какой стороны от стенки расчетная область?
    if(j==1) then
      isOwnRight = .true.                                       ! справа, если грань на левой стенке
    else if(j==nyy) then
      isOwnRight = .false.
    else
      isOwnRight = (c_type(i,j)>=0)
    end if

    if(isOwnRight) then                                         ! живая ячейка справа

      ! вычисляем Мах на грани:
      soundr = sqrt(g * cs_h0(i,j))
      sound = soundr
      uf = cs_v0(i,j)
      M = uf / sound                                            ! Мах на грани

      if(abs(M)>0.9) then
        write(*,*) "ERROR in PhaseSw2-Y: |M|>1. side (i,j):", i, j
        call avost
      endif

      invR = NaN; Gr = NaN

      ! инвариант Q:
      call TransportSwInvY(T_INVQ, i, j, DIRM, invQ)        ! перенос через правую ячейку влево
      Gq = -cs_G0(i,j)

      ! инварианты U и T:
      if(M<-eps) then
        call TransportSwInvY(T_INVU, i, j, DIRM, invU)        ! перенос через правую ячейку влево
        call TransportSwInvY(T_INVD, i, j, DIRM, invD)
      else
        invU = GetCellInvSwY(T_INVU, i, j, .true.)                ! .. вычисляем инвариант в центре ячейки
        invD = GetCellInvSwY(T_INVD, i, j, .true.)
      end if

      ! восстанавливаем значения из инвариантов:
      fyn_rho0(i,j) = invD
      fyn_u0(i,j) = invU

      ! u = 0
      ! invQ = 0 + Gq * (b + h)
      !----------------------------------
      ! u = 0
      ! h = invQ / Gq - b

      fyn_h0(i,j) = invQ / Gq - fy_b(i,j)
      fyn_v0(i,j) = 0.

      if(IsDebSideY(i, j, -1)) then
        write(17,"(a)") "|O<-;i;j;;h;rho;u;v;R;Q;U;D;Gr;Gq"
        write(17,"('side:',2(';',i0),';',100(';',1p,g0.16))") &
            i, j, &
              fy_h0(i,j), fy_rho0(i,j), fy_u0(i,j), fy_v0(i,j), &
              invR, invQ, invU, invD, Gr, Gq
        write(17,"('UFN:',2(';',i0),';',100(';',1p,g0.16))") &
            i, j, &
              fyn_h0(i,j), fyn_rho0(i,j), fyn_u0(i,j), fyn_v0(i,j)
      end if

    else    ! if(c_type(i,j)>=0)                                ! живая ячейка слева, стенка справа

      ! вычисляем число Маха на грани:
      jl = j - 1
      sound = sqrt(g * cs_h0(i,jl))
      uf = cs_v0(i,jl)                                          ! скорость на грани
      M = uf / sound                                            ! Мах на грани

      if(abs(M)>0.9) then
        write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
        call avost
      endif

      ! инвариант R:
      call TransportSwInvY(T_INVR, i, j, DIRP, invR)        ! перенос через левую ячейку вправо
      Gr = cs_G0(i,jl)

      invQ = NaN; Gq = NaN

      ! инварианты U и T:
      if(M>eps) then
        call TransportSwInvY(T_INVU, i, j, DIRP, invU)        ! перенос вправо
        call TransportSwInvY(T_INVD, i, j, DIRP, invD)
      else
        invU = GetCellInvSwY(T_INVU, i, jl, .true.)               ! .. вычисляем инвариант в центре ячейки
        invD = GetCellInvSwY(T_INVD, i, jl, .true.)
      end if

      ! восстанавливаем значения из инвариантов:
      fyn_rho0(i,j) = invD
      fyn_u0(i,j) = invU

      ! u = 0
      ! invR = 0 + Gr * (b + h)
      !----------------------------------
      ! u = 0
      ! h = invR / Gr - b

      fyn_h0(i,j) = invR / Gr - fy_b(i,j)
      fyn_v0(i,j) = 0.

      if(IsDebSideY(i, j, -1)) then
        write(17,"(a)") "->O|;i;j;;h;rho;u;v;R;Q;U;D;Gr;Gq"
        write(17,"('side:',2(';',i0),';',100(';',1p,g0.16))") &
            i, j, &
              fy_h0(i,j), fy_rho0(i,j), fy_u0(i,j), fy_v0(i,j), &
              invR, invQ, invU, invD, Gr, Gq
        write(17,"('UFN:',2(';',i0),';',100(';',1p,g0.16))") &
            i, j, &
              fyn_h0(i,j), fyn_rho0(i,j), fyn_u0(i,j), fyn_v0(i,j)
      end if

    end if  ! if(c_type(i,j)>=0)

    call BuildFluxesSwY(i,j)                                      ! построение потоков на грани

  end do  ! do is

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! ГУ "стенка с прилипанием":
  ns = nfy_sides(BC_STICK)                                      ! количество граней
  do is=1, ns
    i = fy_sides(BC_STICK).ptr(is).i                            ! индексы грани на сетке
    j = fy_sides(BC_STICK).ptr(is).j

    ! с какой стороны от стенки расчетная область?
    if(j==1) then
      isOwnRight = .true.                                       ! справа, если грань на левой стенке
    else if(j==nyy) then
      isOwnRight = .false.
    else
      isOwnRight = (c_type(i,j)>=0)
    end if

    if(isOwnRight) then                                         ! живая ячейка справа

      ! вычисляем Мах на грани:
      soundr = sqrt(g * cs_h0(i,j))
      sound = soundr
      uf = cs_v0(i,j)
      M = uf / sound                                            ! Мах на грани

      if(abs(M)>0.9) then
        write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
        call avost
      endif

      ! инвариант Q:
      call TransportSwInvY(T_INVQ, i, j, DIRM, invQ)            ! перенос через правую ячейку влево
      Gq = -cs_G0(i,j)

      ! инвариант D:
      if(M<-eps) then
        call TransportSwInvY(T_INVD, i, j, DIRM, invD)
      else
        invD = GetCellInvSwY(T_INVD, i, j, .true.)                ! .. вычисляем инвариант в центре ячейки
      end if

      ! восстанавливаем значения из инвариантов:
      fyn_rho0(i,j) = invD
      fyn_u0(i,j) = 0.

      ! u = 0
      ! invQ = 0 + Gq * (b + h)
      !----------------------------------
      ! u = 0
      ! h = invQ / Gq - b

      fyn_h0(i,j) = invQ / Gq - fy_b(i,j)
      fyn_v0(i,j) = 0.

    else    ! if(c_type(i,j)>=0)                                ! живая ячейка слева, стенка справа

      ! вычисляем число Маха на грани:
      jl = j - 1
      sound = sqrt(g * cs_h0(i,jl))
      uf = cs_v0(i,jl)                                         ! скорость на грани
      M = uf / sound                                            ! Мах на грани

      if(abs(M)>0.9) then
        write(*,*) "ERROR in PhaseSw2X: |M|>1. side (i,j):", i, j
        call avost
      endif

      ! инвариант R:
      call TransportSwInvY(T_INVR, i, j, DIRP, invR)        ! перенос через левую ячейку вправо
      Gr = cs_G0(i, jl)

      ! инвариант D:
      if(M>eps) then
        call TransportSwInvY(T_INVD, i, j, DIRP, invD)
      else
        invD = GetCellInvSwY(T_INVD, i, jl, .true.)                ! .. вычисляем инвариант в центре ячейки
      end if

      ! восстанавливаем значения из инвариантов:
      fyn_rho0(i,j) = invD
      fyn_u0(i,j) = 0.

      ! u = 0
      ! invR = 0 + Gr * (b + h)
      !----------------------------------
      ! u = 0
      ! h = invR / Gr - b

      fyn_h0(i,j) = invR / Gr - fy_b(i,j)
      fyn_v0(i,j) = 0.

    end if  ! if(c_type(i,j)>=0)

    call BuildFluxesSwY(i,j)                                      ! построение потоков на грани

  end do  ! do is

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! ГУ "постоянный вток":
  ns = nfy_sides(BC_IN)                                         ! количество граней
  do is=1, ns
    i = fy_sides(BC_IN).ptr(is).i                             ! индексы грани на сетке
    j = fy_sides(BC_IN).ptr(is).j
    jl = j - 1
    jr = j

    ibc = fy_bc(i,j)                                            ! индекс индивидуальных параметров ГУ

    ! заграничные данные:
    bcH = bcData(ibc).h
    bcU = bcData(ibc).u
    bcV = bcData(ibc).v
    bcRho = bcData(ibc).rho

    if(IsDebSideY(i, j, -1)) then
      write(17,"(/'Phase Sw2-Y',2(';',i0),';',';BC_IN')") i, j
    endif

    isOwnRight = fy_IsOwnRight(i,j)

    jc = j - 1; if(isOwnRight) jc = j

    ! вычисляем Мах на грани:
    soundr = sqrt(g * bcH)
    soundl = sqrt(g * cs_h0(i,jc))
    sound = (soundl + soundr) / 2.
    uf = 0.5*(cs_u0(i,jc) + bcU)                                ! скорость на грани
    M = uf / sound                                              ! Мах на грани

    if(abs(M)>0.9) then
      write(*,*) "ERROR in PhaseSw2Y: |M|>1. side (i,j):", i, j
      call avost
    endif

    ! инвариант R (слева направо):
    if(isOwnRight) then                                         ! ячейка справа, граница слева
      Gr = sqrt(g / bcH)
      invR = bcU + Gr * bcH
    else                                                        ! слева есть ячейка
      call TransportSwInvY(T_INVR, i, j, DIRP, invR)            ! .. перенос через левую ячейку вправо
      Gr = cs_G0(i,jc)
    end if

    ! инвариант Q (справа налево):
    if(isOwnRight) then                                         ! справа есть ячейка
      call TransportSwInvY(T_INVQ, i, j, DIRM, invQ)            ! перенос через правую ячейку влево
      Gq = -cs_G0(i,jc)
    else
      Gq = -sqrt(g / bcH)
      invQ = bcU + Gq * bcH
    end if

    ! инварианты U и Rho:
      ! транспортные инварианты:
    if(isOwnRight .and. M<-eps) then                            ! течение из ячейки на границу справа налево  |<O
      call TransportSwInvY(T_INVU, i, j, DIRM, invU)            ! .. перенос через правую ячейку влево
      call TransportSwInvY(T_INVD, i, j, DIRM, invD)            ! .. перенос через правую ячейку влево
    else if(.not.isOwnRight .and. M>eps) then                   ! течение из ячейки на границу слева направо  O>|
      call TransportSwInvY(T_INVU, i, j, DIRP, invU)            ! .. перенос через левую ячейку вправо
      call TransportSwInvY(T_INVD, i, j, DIRP, invD)            ! .. перенос через левую ячейку вправо
    else                                                        ! течение из-за границы
      invU = bcV
      invD = bcRho
    end if

    ! восстанавливаем значения из инвариантов:
    rhoi = invD                                                 ! плотность
    ui = invU                                                   ! тангенцальная скорость

    ! invR = u + Gr * (b + h)
    ! invQ = u + Gq * (b + h)
    !----------------------------------
    ! invR - invQ = (Gr - Gq) * (b + h)
    ! invR*Gq - invQ*Gr = (Gq - Gr)*u

    hi = (invR - invQ) / (Gr - Gq) - fy_b(i,j)                  ! глубина
    vi = (invR*Gq - invQ*Gr) / (Gq - Gr)                        ! нормальная скорость

    fyn_h0(i,j) = hi
    fyn_u0(i,j) = ui
    fyn_v0(i,j) = vi                                            ! тангенцальная скорость
    fyn_rho0(i,j) = rhoi                                        ! плотность

    if(IsDebSideY(i, j, -1)) then
      write(17,"(a)") ";i;j;;h;rho;u;v;R;Q;U;D;Gr;Gq"
      write(17,"('side:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            bcH, bcRho, bcU, bcV
      write(17,"('UF:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fy_h0(i,j), fy_rho0(i,j), fy_u0(i,j), fy_v0(i,j), &
            invR, invQ, invU, invD, Gr, Gq
      write(17,"('UFN:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fyn_h0(i,j), fyn_rho0(i,j), fyn_u0(i,j), fyn_v0(i,j)
    end if

    call BuildFluxesSwY(i,j)                                      ! построение потоков на грани

  end do  ! do is

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! ГУ "переменный вток":
  ns = nfy_sides(BC_IN_T)                                       ! количество граней
  do is=1, ns
    i = fy_sides(BC_IN_T).ptr(is).i                             ! индексы грани на сетке
    j = fy_sides(BC_IN_T).ptr(is).j
    jl = j - 1
    jr = j

    ibc = fy_bc(i,j)                                            ! индекс индивидуальных параметров ГУ

    ! заграничные данные:
    bcH = bcData(ibc).h
    bcU = bcData(ibc).u
    bcV = bcData(ibc).v
    bcRho = bcData(ibc).rho

    if(IsDebSideY(i, j, -1)) then
      write(17,"(/'Phase Sw2-Y',2(';',i0),';',';BC_IN_T')") i, j
    endif

    isOwnRight = fy_IsOwnRight(i,j)

    jc = j - 1; if(isOwnRight) jc = j

    ! вычисляем Мах на грани:
    soundr = sqrt(g * bcH)
    soundl = sqrt(g * cs_h0(i,jc))
    sound = (soundl + soundr) / 2.
    uf = 0.5*(cs_u0(i,jc) + bcU)                                ! скорость на грани
    M = uf / sound                                              ! Мах на грани

    if(abs(M)>0.9) then
      write(*,*) "ERROR in PhaseSw2Y: |M|>1. side (i,j):", i, j
      call avost
    endif

    ! инвариант R (слева направо):
    if(isOwnRight) then                                         ! ячейка справа, граница слева
      Gr = sqrt(g / bcH)
      invR = bcU + Gr * bcH
    else                                                        ! слева есть ячейка
      call TransportSwInvY(T_INVR, i, j, DIRP, invR)            ! .. перенос через левую ячейку вправо
      Gr = cs_G0(i,jc)
    end if

    ! инвариант Q (справа налево):
    if(isOwnRight) then                                         ! справа есть ячейка
      call TransportSwInvY(T_INVQ, i, j, DIRM, invQ)            ! перенос через правую ячейку влево
      Gq = -cs_G0(i,jc)
    else
      Gq = -sqrt(g / bcH)
      invQ = bcU + Gq * bcH
    end if

    ! инварианты U и Rho:
      ! транспортные инварианты:
    if(isOwnRight .and. M<-eps) then                            ! течение из ячейки на границу справа налево  |<O
      call TransportSwInvY(T_INVU, i, j, DIRM, invU)            ! .. перенос через правую ячейку влево
      call TransportSwInvY(T_INVD, i, j, DIRM, invD)            ! .. перенос через правую ячейку влево
    else if(.not.isOwnRight .and. M>eps) then                   ! течение из ячейки на границу слева направо  O>|
      call TransportSwInvY(T_INVU, i, j, DIRP, invU)            ! .. перенос через левую ячейку вправо
      call TransportSwInvY(T_INVD, i, j, DIRP, invD)            ! .. перенос через левую ячейку вправо
    else                                                        ! течение из-за границы
      invU = bcV
      invD = bcRho
    end if

    ! восстанавливаем значения из инвариантов:
    rhoi = invD                                                 ! плотность
    ui = invU                                                   ! тангенцальная скорость

    ! invR = u + Gr * (b + h)
    ! invQ = u + Gq * (b + h)
    !----------------------------------
    ! invR - invQ = (Gr - Gq) * (b + h)
    ! invR*Gq - invQ*Gr = (Gq - Gr)*u

    hi = (invR - invQ) / (Gr - Gq) - fy_b(i,j)                  ! глубина
    vi = (invR*Gq - invQ*Gr) / (Gq - Gr)                        ! нормальная скорость

    fyn_h0(i,j) = hi
    fyn_u0(i,j) = ui
    fyn_v0(i,j) = vi                                            ! тангенцальная скорость
    fyn_rho0(i,j) = rhoi                                        ! плотность

    if(IsDebSideY(i, j, -1)) then
      write(17,"(a)") ";i;j;;h;rho;u;v;R;Q;U;D;Gr;Gq"
      write(17,"('side:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            bcH, bcRho, bcU, bcV
      write(17,"('UF:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fy_h0(i,j), fy_rho0(i,j), fy_u0(i,j), fy_v0(i,j), &
            invR, invQ, invU, invD, Gr, Gq
      write(17,"('UFN:',2(';',i0),';',100(';',1p,g0.16))") &
          i, j, &
            fyn_h0(i,j), fyn_rho0(i,j), fyn_u0(i,j), fyn_v0(i,j)
    end if

    call BuildFluxesSwY(i,j)                                      ! построение потоков на грани

  end do  ! do is

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! ГУ "фиксированная скорость":
  ns = nfy_sides(BC_FIX_U)                                       ! количество граней
  do is=1, ns
    i = fy_sides(BC_FIX_U).ptr(is).i                             ! индексы грани на сетке
    j = fy_sides(BC_FIX_U).ptr(is).j
    jl = j - 1
    jr = j
    call avost("Не реализовано. " // __FILE__, __LINE__)

    call BuildFluxesSwY(i,j)                                    ! построение потоков на грани

  end do  ! do is

end subroutine PhaseSw2Y

!------------------------------------------------------------------------------------------------------------------------------------

! перенос инварианта type через ячейку на грань (i,j) в направлении off
subroutine TransportSwInvX(itype, i, j, off, inv)
  use variables
  implicit none

  interface
    function GetCellInvSwX(itype, i, j, isN); integer(4)::itype,i,j; logical::isN; real(8)::GetCellInvSwX; end
    function GetSideInvSwX(itype, i, j, GG); integer(4)::itype,i,j; real(8)::GetSideInvSwX, GG; end
  end interface

  integer(4):: itype, i, j, off
  real(R8):: inv
  integer(4):: io, ic
  real(R8):: invo, invt, invc, invs, invn, sound, lam, dx, q, Imin, Imax
  character(1):: iname(4)

  call assertn(itype>0 .and. itype<5, "TransportSwInvX. i,j,itype:", i, j, itype)

  io = i - (2*off-1)                                            ! индекс заднего ребра
  ic = i - off                                                  ! индекс ячейки

  invo = GetSideInvSwX(itype, io, j, cs_G0(ic, j))                             ! инвариант на задней грани
  invt = GetSideInvSwX(itype, i, j, cs_G0(ic, j))                              ! инвариант на передней грани
  invc = GetCellInvSwX(itype, ic, j, .false.)                    ! инвариант в ячейке t[n]
  invs = GetCellInvSwX(itype, ic, j, .true.)                     ! инвариант в ячейке t[n+1/2]

  invn = 2.*invs - invo                                         ! новый инвариант предвариательно

  sound = sqrt(g * cs_h0(ic,j))                                 ! скорость звука в ячейке

  ! вычисляем характеристическую скорость lambda:
  select case(itype)
    case (T_INVR)
      lam = cs_u0(ic,j) + sound                                 ! для инварианта R: lam = u + c
    case (T_INVQ)
      lam = cs_u0(ic,j) - sound                                 ! для инварианта Q: lam = u - c
    case (T_INVU,T_INVD)
      lam = cs_u0(ic,j)                                         ! для транспортных инвариантов: lam = u
  end select
  if(off==0) lam = -lam                                         ! если перенос влево
  dx = c_dx(ic)

  q = 2. * (invs - invc) + lam * (invt - invo) / dx * dt              ! аппроксимация правой части
  Imin = min(invo, invc, invt) + q                              ! нижня мажоранта
  Imax = max(invo, invc, invt) + q                              ! верхняя мажоранта

  inv = min(Imax, max(Imin, invn))                              ! лимитирование

  if(IsDebSideX(i, j, -1)) then
    iname = ['R', 'Q', 'U', 'D']
    write(17,"(/'TranspI-SW-X:',2(';',i0),';;off=',i0)") i, j, off
    write(17,"(a)") ";i;j;;I;rho;u;v;G"
    write(17,"(a,a,2(';',i0),';',100(';',1p,g0.16))") &
      iname(itype), 'o', io, j, invo, fx_rho0(io,j), fx_u0(io,j), fx_v0(io,j)
    write(17,"(a,a,2(';',i0),';',100(';',1p,g0.16))") &
      iname(itype), 't', i , j, invt, fx_rho0(i ,j), fx_u0(i ,j), fx_v0(i ,j)
    write(17,"(a,a,2(';',i0),';',100(';',1p,g0.16))") &
      iname(itype), 'c', ic, j, invc, c_rho0(ic,j), c_u0(ic,j), c_v0(ic,j), cs_G0(ic, j)
    write(17,"(a,a,2(';',i0),';',100(';',1p,g0.16))") &
      iname(itype), 's', ic, j, invs, cs_rho0(ic,j), cs_u0(ic,j), cs_v0(ic,j), cs_G0(ic, j)
    write(17,"(a,a,';;;;',1p,g0.16,a)") iname(itype), 'n', invn, &
          ";u;sound;lam;dx;dt;q;Imin;Imax"
    write(17,"(a,a,2(';',i0),';',100(';',1p,g0.16))") iname(itype), ' ', i , j, inv, &
          cs_u0(ic,j), sound, lam, dx, dt, q, Imin, Imax
  end if

end subroutine TransportSwInvX

!------------------------------------------------------------------------------------------------------------------------------------

! вычисление инварианта типа type в ячейке (i,j) на слое isN
function GetCellInvSwX(type, i, j, isN)
  use variables
  implicit none
  integer(4):: type, i, j
  logical:: isN
  real(8):: GetCellInvSwX, inv, GG

  select case(type)
    case (T_INVR)
      GG = cs_G0(i,j)
      if(isN) then; inv = cs_u0(i,j) + GG * cs_z0(i,j); else; inv = c_u0(i,j) + GG * c_z0(i,j); endif
    case (T_INVQ)
      GG = cs_G0(i,j)
      if(isN) then; inv = cs_u0(i,j) - GG * cs_z0(i,j); else; inv = c_u0(i,j) - GG * c_z0(i,j); endif
    case (T_INVU)
      if(isN) then; inv = cs_v0(i,j); else; inv = c_v0(i,j); endif
    case (T_INVD)
      if(isN) then; inv = cs_rho0(i,j); else; inv = c_rho0(i,j); endif
    case default
      write(*,*) "ERROR GetCellInvSwX. type:", type
      call Exit
  end select

  GetCellInvSwX = inv

end function

!------------------------------------------------------------------------------------------------------------------------------------

! вычисление инварианта типа type на грани-X (i,j) на слое t[n]
function GetSideInvSwX(type, i, j, GG)
  use variables
  implicit none
  integer(4):: type, i, j
  real(8):: GetSideInvSwX, inv, GG

  select case(type)
    case (T_INVR)
      inv = fx_u0(i,j) + GG * fx_z0(i,j)
    case (T_INVQ)
      inv = fx_u0(i,j) - GG * fx_z0(i,j)
    case (T_INVU)
      inv = fx_v0(i,j)
    case (T_INVD)
      inv = fx_rho0(i,j)
    case default
      write(*,*) "ERROR GetCellInvSwX. type:", type
      call Exit
  end select

  GetSideInvSwX = inv

end function

!------------------------------------------------------------------------------------------------------------------------------------

subroutine BuildFluxesSwX(i, j)
  use variables
  implicit none
  integer(4):: i, j
  real(R8):: fH, fM, fIx, fIy

  fxn_z0(i,j) = fx_b(i,j) + fxn_h0(i,j)                         ! уровень поверхности

  fH = fxn_h0(i,j) * fxn_u0(i,j)                                ! поток объема
  fM = fH * fxn_rho0(i,j)                                       ! поток массы

  ! ниже используется приближение Буссинеска (постоянная плотность):
  fIx = fH * fxn_u0(i,j) + g/2. * fxn_z0(i,j)**2
  !fIx = fH * fxn_u0(i,j) + g/2. * fxn_h0(i,j)**2
  fIy = fH * fxn_v0(i,j)

  fx_fluxH0(i,j) = fH
  fx_fluxM0(i,j) = fM
  fx_fluxIx0(i,j) = fIx
  fx_fluxIy0(i,j) = fIy

  if(IsDebSideX(i, j, -1)) then
    debdum = 0
!    write(17,"(a)") "BuildFluxesSwX;i;j;h;z;rho;u;v;fH;fM;fIx;fIy"
!    write(17,"(2(';',i0),100(';',1p,g0.16))") i, j, fxn_h0(i,j), fxn_z0(i,j), fxn_rho0(i,j), fxn_u0(i,j), fxn_v0(i,j), &
!                   fH, fM, fIx, fIy
  endif
end subroutine

!------------------------------------------------------------------------------------------------------------------------------------

! перенос инварианта type через ячейку на грань (i,j) в направлении off
subroutine TransportSwInvY(itype, i, j, off, inv)
  use variables
  implicit none
  integer(4):: itype, i, j, off
  real(R8):: inv

  interface
    function GetCellInvSwY(itype, i, j, isN); integer(4)::itype,i,j; logical::isN; real(8)::GetCellInvSwY; end
    function GetSideInvSwY(itype, i, j, GG); integer(4)::itype,i,j; real(8)::GetSideInvSwY, GG; end
  end interface

  integer(4):: jo, jc
  real(R8):: invo, invt, invc, invs, invn, sound, lam, dy, q, Imin, Imax
  character(1):: iname(4)

  jo = j - (2*off-1)                                            ! индекс заднего ребра
  jc = j - off                                                  ! индекс ячейки

  invo = GetSideInvSwY(itype, i, jo, cs_G0(i, jc))              ! инвариант на задней грани
  invt = GetSideInvSwY(itype, i, j, cs_G0(i, jc))               ! инвариант на передней грани
  invc = GetCellInvSwY(itype, i, jc, .false.)                   ! инвариант в ячейке t[n]
  invs = GetCellInvSwY(itype, i, jc, .true.)                    ! инвариант в ячейке t[n+1/2]

  invn = 2.*invs - invo                                         ! новый инвариант предвариательно

  sound = sqrt(g * cs_h0(i,jc))                                 ! скорость звука в ячейке

  ! вычисляем характеристическую скорость lambda:
  select case(itype)
    case (T_INVR)
      lam = cs_v0(i,jc) + sound                                 ! для инварианта R: lam = u + c
    case (T_INVQ)
      lam = cs_v0(i,jc) - sound                                 ! для инварианта Q: lam = u - c
    case (T_INVU,T_INVD)
      lam = cs_v0(i,jc)                                         ! для транспортных инвариантов: lam = u
  end select
  if(off==0) lam = -lam                                         ! если перенос влево
  dy = c_dy(jc)

  q = 2. * (invs - invc) + lam * (invt - invo) / dy * dt        ! аппроксимация правой части
  Imin = min(invo, invc, invt) + q                              ! нижня мажоранта
  Imax = max(invo, invc, invt) + q                              ! верхняя мажоранта

  inv = min(Imax, max(Imin, invn))                              ! лимитирование

  if(IsDebSideY(i, j, -1)) then
    iname = ['R', 'Q', 'U', 'D']
    write(17,"(/'TranspI-SW-Y:',2(';',i0),';;off=',i0)") i, j, off
    write(17,"(a)") ";i;j;;I;rho;u;v;G"
    write(17,"(a,a,2(';',i0),';',100(';',1p,g0.16))") &
      iname(itype), 'o', i, jo, invo, fy_rho0(i,jo), fy_u0(i,jo), fy_v0(i,jo)
    write(17,"(a,a,2(';',i0),';',100(';',1p,g0.16))") &
      iname(itype), 't', i , j, invt, fy_rho0(i ,j), fy_u0(i ,j), fy_v0(i ,j)
    write(17,"(a,a,2(';',i0),';',100(';',1p,g0.16))") &
      iname(itype), 'c', i, jc, invc, c_rho0(i,jc), c_u0(i,jc), c_v0(i,jc), cs_G0(i, jc)
    write(17,"(a,a,2(';',i0),';',100(';',1p,g0.16))") &
      iname(itype), 's', i, jc, invs, cs_rho0(i,jc), cs_u0(i,jc), cs_v0(i,jc), cs_G0(i, jc)
    write(17,"(a,a,';;;;',1p,g0.16,a)") iname(itype), 'n', invn, &
          ";u;sound;lam;dy;dt;q;Imin;Imax"
    write(17,"(a,a,2(';',i0),';',100(';',1p,g0.16))") iname(itype), ' ', i , j, inv, &
          cs_v0(i,jc), sound, lam, dy, dt, q, Imin, Imax
  end if

end subroutine TransportSwInvY

!------------------------------------------------------------------------------------------------------------------------------------

! вычисление инварианта типа type в ячейке (i,j) на слое isN
function GetCellInvSwY(type, i, j, isN)
  use variables
  implicit none
  integer(4):: type, i, j
  logical:: isN
  real(8):: GetCellInvSwY, inv, GG

  select case(type)
    case (T_INVR)
      GG = cs_G0(i,j)
      if(isN) then; inv = cs_v0(i,j) + GG * cs_z0(i,j); else; inv = c_v0(i,j) + GG * c_z0(i,j); endif
    case (T_INVQ)
      GG = cs_G0(i,j)
      if(isN) then; inv = cs_v0(i,j) - GG * cs_z0(i,j); else; inv = c_v0(i,j) - GG * c_z0(i,j); endif
    case (T_INVU)
      if(isN) then; inv = cs_u0(i,j); else; inv = c_u0(i,j); endif
    case (T_INVD)
      if(isN) then; inv = cs_rho0(i,j); else; inv = c_rho0(i,j); endif
    case default
      write(*,*) "ERROR GetCellInvSwY. type:", type
      call Exit
  end select

  GetCellInvSwY = inv

end function

!------------------------------------------------------------------------------------------------------------------------------------

! вычисление инварианта типа type на грани-Y (i,j) на слое t[n]
function GetSideInvSwY(type, i, j, GG)
  use variables
  implicit none
  integer(4):: type, i, j
  real(8):: GetSideInvSwY, inv, GG

  select case(type)
    case (T_INVR)
      inv = fy_v0(i,j) + GG * fy_z0(i,j)
    case (T_INVQ)
      inv = fy_v0(i,j) - GG * fy_z0(i,j)
    case (T_INVU)
      inv = fy_u0(i,j)
    case (T_INVD)
      inv = fy_rho0(i,j)
    case default
      write(*,*) "ERROR GetCellInvSwY. type:", type
      call Exit
  end select

  GetSideInvSwY = inv

end function

!------------------------------------------------------------------------------------------------------------------------------------

subroutine BuildFluxesSwY(i, j)
  use variables
  implicit none
  integer(4):: i, j
  real(8):: fH, fM, fIx, fIy

  fyn_z0(i,j) = fy_b(i,j) + fyn_h0(i,j)                         ! уровень поверхности

  fH = fyn_h0(i,j) * fyn_v0(i,j)                                ! поток объема
  fM = fH * fyn_rho0(i,j)                                       ! поток массы

  ! ниже используется приближение Буссинеска (постоянная плотность):
  fIx = fH * fyn_u0(i,j)
  fIy = fH * fyn_v0(i,j) + g/2. * fyn_z0(i,j)**2
  !fIy = fH * fyn_v0(i,j) + g/2. * fyn_h0(i,j)**2

  fy_fluxH0(i,j) = fH
  fy_fluxM0(i,j) = fM
  fy_fluxIx0(i,j) = fIx
  fy_fluxIy0(i,j) = fIy

  if(IsDebSideY(i, j, -1)) then
    write(17,"(a)") "BuildFluxesSwY;i;j;;h;z;rho;u;v;fH;fM;fIx;fIy"
    write(17,"(2(';',i0),';',100(';',1p,g0.16))") i, j, fyn_h0(i,j), fyn_z0(i,j), fyn_rho0(i,j), fyn_u0(i,j), fyn_v0(i,j), &
                   fH, fM, fIx, fIy
  endif
end subroutine

!------------------------------------------------------------------------------------------------------------------------------------

! суммирование потоков с граней в ячейку
subroutine SumFluxesSw
  use variables
  implicit none
  integer(4):: i, j, ic
  real(8):: dx, dy

  do ic=1,ncells
    i=cells.ptr(ic).i
    j=cells.ptr(ic).j

    dx = c_dx(i)
    dy = c_dy(j)
    c_sumFluxH0(i,j)  = (fx_fluxH0(i+1,j) - fx_fluxH0(i,j)) / dx + &
                        (fy_fluxH0(i,j+1) - fy_fluxH0(i,j)) / dy
    c_sumFluxM0(i,j)  = (fx_fluxM0(i+1,j) - fx_fluxM0(i,j)) / dx + &
                        (fy_fluxM0(i,j+1) - fy_fluxM0(i,j)) / dy
    c_sumFluxIx0(i,j) = (fx_fluxIx0(i+1,j) - fx_fluxIx0(i,j)) / dx + &
                        (fy_fluxIx0(i,j+1) - fy_fluxIx0(i,j)) / dy - &
!                          0
                        g * c_b(i,j) * (fxn_z0(i+1,j) - fxn_z0(i,j)) / dx
    c_sumFluxIy0(i,j) = (fx_fluxIy0(i+1,j) - fx_fluxIy0(i,j)) / dx + &
                        (fy_fluxIy0(i,j+1) - fy_fluxIy0(i,j)) / dy - &
!                          0
                        g * c_b(i,j) * (fyn_z0(i,j+1) - fyn_z0(i,j)) / dy

    if(IsDebCell(i, j, -1)) then
      write(17,"(/'SumFluxesSw',2(';',i0),';',100(';',1p,g0.16))") i, j
      write(17,"(a)") "side;i;j;;fH;fM;fIx;fIy"
      write(17,"('L',2(';',i0),';',100(';',1p,g0.16))") i  , j  , fx_fluxH0(i  ,j  ), fx_fluxM0(i  ,j  ), fx_fluxIx0(i  ,j  ), fx_fluxIy0(i  ,j  )
      write(17,"('R',2(';',i0),';',100(';',1p,g0.16))") i+1, j  , fx_fluxH0(i+1,j  ), fx_fluxM0(i+1,j  ), fx_fluxIx0(i+1,j  ), fx_fluxIy0(i+1,j  )
      write(17,"('B',2(';',i0),';',100(';',1p,g0.16))") i  , j  , fy_fluxH0(i  ,j  ), fy_fluxM0(i  ,j  ), fy_fluxIx0(i  ,j  ), fy_fluxIy0(i  ,j  )
      write(17,"('T',2(';',i0),';',100(';',1p,g0.16))") i  , j+1, fy_fluxH0(i  ,j+1), fy_fluxM0(i  ,j+1), fy_fluxIx0(i  ,j+1), fy_fluxIy0(i  ,j+1)
      write(17,"('Sum:;;;',100(';',1p,g0.16))") c_sumFluxH0(i,j), c_sumFluxM0(i,j), c_sumFluxIx0(i,j), c_sumFluxIy0(i,j)
    endif

  end do
end subroutine SumFluxesSw

!-------------------------------------------------------------------------------------------------------------

! начальное построение потоков на гранях по начальным данным
subroutine BuildFluxesSw
  use variables
  implicit none

  integer(4):: i, j

  do i=1,nxx
    do j=1, nxy
      if(fx_type(i,j)>=0) then
        call BuildFluxesSwX(i,j)
      end if
    end do
  end do

  do i=1,nyx
    do j=1, nyy
      if(fy_type(i,j)>=0) then
        call BuildFluxesSwY(i,j)
      end if
    end do
  end do

  call SumFluxesSw

end subroutine

!----------------------------------------------

subroutine CalcSwW
  use variables
!  , only: assertn, avost, c_b, c_dx, c_dy, c_h0, c_rho0, c_type, c_u0, c_v0, c_wb0, c_wt0, &
!                       cells, cn_h0, cn_rho0, cn_u0, cn_v0, cn_wb0, cn_wt0, cs_h0, cs_rho0, cs_u0, cs_v0, &
!                       cs_wb0, cs_wt0, debdum, dt, eps, fx_b, fx_h0, fx_rho0, fx_sides, fx_type, fx_u0, fx_v0, &
!                       fx_wb0, fx_wt0, fxn_h0, fxn_rho0, fxn_u0, fxn_v0, fxn_wb0, fxn_wt0, fxn_z0, fy_b, &
!                       fy_h0, fy_rho0, fy_sides, fy_type, fy_u0, fy_v0, fy_wb0, fy_wt0, fyn_h0, fyn_rho0, &
!                       fyn_u0, fyn_v0, fyn_wb0, fyn_wt0, fyn_z0, g, ncells, ncx, ncy, nfx_sides, nfy_sides, &
!                       nxx, nxy, nyx, nyy, x, y, BC_PERIODIC
  implicit none

  integer(4):: ic, i, j, cl, cr
  real(8):: dx, dy, h, h4, dudx, dvdy

  ! вертикальная скорость в ячейках:
  do ic=1,ncells
    i=cells.ptr(ic).i
    j=cells.ptr(ic).j

    ! вертикальная скорость t[n]:
#if 0
    c_wb0(i,j) = (fx_u0(i+1,j) + fx_u0(i,j)) / 2 * (fx_b(i+1,j) - fx_b(i,j)) / c_dx(i) + &
                 (fy_v0(i,j+1) + fy_v0(i,j)) / 2 * (fy_b(i,j+1) - fy_b(i,j)) / c_dy(j)
#else
    c_wb0(i,j) = c_u0(i,j) * (fx_b(i+1,j) - fx_b(i,j)) / c_dx(i) + &
                 c_v0(i,j) * (fy_b(i,j+1) - fy_b(i,j)) / c_dy(j)
    cs_wb0(i,j) = cs_u0(i,j) * (fx_b(i+1,j) - fx_b(i,j)) / c_dx(i) + &
                  cs_v0(i,j) * (fy_b(i,j+1) - fy_b(i,j)) / c_dy(j)
#endif
    c_wt0(i,j) = c_wb0(i,j) - &
                 (fx_h0(i+1,j) + fx_h0(i,j)) / 2 * (fx_u0(i+1,j) - fx_u0(i,j)) / c_dx(i) - &
                 (fy_h0(i,j+1) + fy_h0(i,j)) / 2 * (fy_v0(i,j+1) - fy_v0(i,j)) / c_dy(j)
    ! вертикальная скорость t[n+1]:
#if 0
    cn_wb0(i,j) = (fxn_u0(i+1,j) + fxn_u0(i,j)) / 2 * (fx_b(i+1,j) - fx_b(i,j)) / c_dx(i) + &
                  (fyn_v0(i,j+1) + fyn_v0(i,j)) / 2 * (fy_b(i,j+1) - fy_b(i,j)) / c_dy(j)
#else
    cn_wb0(i,j) = cn_u0(i,j) * (fx_b(i+1,j) - fx_b(i,j)) / c_dx(i) + &
                  cn_v0(i,j) * (fy_b(i,j+1) - fy_b(i,j)) / c_dy(j)
#endif
    cn_wt0(i,j) = cn_wb0(i,j) - &
                  (fxn_h0(i+1,j) + fxn_h0(i,j)) / 2 * (fxn_u0(i+1,j) - fxn_u0(i,j)) / c_dx(i) - &
                  (fyn_h0(i,j+1) + fyn_h0(i,j)) / 2 * (fyn_v0(i,j+1) - fyn_v0(i,j)) / c_dy(j)
    ! вертикальная скорость t[n+1/2]:
#if 0
    cs_wb0(i,j) = (c_wb0(i,j) + cn_wb0(i,j)) / 2.
    cs_wt0(i,j) = (c_wt0(i,j) + cn_wt0(i,j)) / 2.
#else
    dx = c_dx(i)
    dy = c_dy(j)
#  if 1
    h = ( fx_h0(i+1,j)+fxn_h0(i+1,j) + fx_h0(i,j)+fxn_h0(i,j) + &
          fy_h0(i,j+1)+fyn_h0(i,j+1) + fy_h0(i,j)+fyn_h0(i,j) ) / 8.
#  else
!#warning ERROR!!!
    h = ( fx_h0(i+1,j)+fxn_h0(i+1,j) + fx_h0(i,j)+fxn_h0(i,j) ) / 4.
    !h4 = ( fx_h0(i+1,j)+fxn_h0(i+1,j) + fx_h0(i,j)+fxn_h0(i,j) )
#  endif
    dudx = ( (fx_u0(i+1,j)+fxn_u0(i+1,j))/2. - (fx_u0(i,j)+fxn_u0(i,j))/2. ) / dx
    dvdy = ( (fy_v0(i,j+1)+fyn_v0(i,j+1))/2. - (fy_v0(i,j)+fyn_v0(i,j))/2. ) / dy
    !cs_wt0(i,j) = cs_wb0(i,j) - h * (dudx + dvdy)
    !cs_wt0(i,j) = cs_wb0(i,j) - 0.125 * h4 * (fx_u0(i+1,j)+fxn_u0(i+1,j)-fx_u0(i,j)-fxn_u0(i,j)) / dx
    !cs_wt0(i,j) = cs_wt0(i,j) - 0.125 * h4 * (fy_v0(i,j+1)+fyn_v0(i,j+1)-fy_v0(i,j)-fyn_v0(i,j)) / dy
#endif

    if(IsDebCell(i, j, -1)) then
      write(17,"(/'CalcSwW',2(';',i0),';',100(';',f0.5))") i, j, (x(i)+x(i+1))/2., (y(j)+y(j+1))/2.
      write(17,"(a)") ";i;j;;c_wb0;c_wt0;cs_wb0;cs_wt0;cn_wb0;cn_wt0;h;dudx;dvdy"
      write(17,"(2(';',i0),';',100(';',1p,g0.16))") i, j, c_wb0(i,j), c_wt0(i,j), cs_wb0(i,j), cs_wt0(i,j), cn_wb0(i,j), cn_wt0(i,j), h, dudx, dvdy
    end if

    debdum=1
  end do

  ! вертикальная скорость на X-гранях:
  do i=1,nxx
    do j=1, nxy

      if(fx_type(i,j)==BC_DELETED) cycle

      cl = i - 1                                               ! ячейка слева
      cr = i                                                   ! ячейка справа

      ! на граничных гранях корректируем индексы ячеек:
      if(fx_type(i,j)==BC_PERIODIC) then                        ! для периодических граней
        if(i==1) cl = ncx
        if(i==nxx) cr = 1
      else if(fx_type(i,j)>0) then                              ! любая другая граничная грань
        if(i==1) then; cl = cr
        else if(c_type(cl,j)<0) then; cl = cr; end if
        if(i==nxx) then; cr = cl
        else if(c_type(cr,j)<0) then; cr = cl; end if
      end if

      call assertn(c_type(cl,j)>=0 .and. c_type(cr,j)>=0, "CalcW_X. Ошибка сетки", i, j, -1)
      fx_wb0(i,j) = (c_wb0(cl,j) + c_wb0(cr,j)) / 2.          ! вертикальная скорость на дне t[n]
      fx_wt0(i,j) = (c_wt0(cl,j) + c_wt0(cr,j)) / 2.          ! вертикальная скорость на поверхности t[n]
      fxn_wb0(i,j) = (cn_wb0(cl,j) + cn_wb0(cr,j)) / 2.       ! вертикальная скорость на дне t[n+1]
      fxn_wt0(i,j) = (cn_wt0(cl,j) + cn_wt0(cr,j)) / 2.       ! вертикальная скорость на поверхности t[n+1]
    end do
  end do

  ! вертикальная скорость на Y-гранях:
  do i=1,nyx
    do j=1, nyy

      if(fy_type(i,j)==BC_DELETED) cycle

      cl = j - 1                                               ! ячейка слева
      cr = j                                                   ! ячейка справа

      ! на граничных гранях корректируем индексы ячеек:
      if(fy_type(i,j)==BC_PERIODIC) then                        ! для периодических граней
        if(j==1) cl = ncy
        if(j==nyy) cr = 1
      else if(fy_type(i,j)>0) then                              ! любая другая граничная грань
        if(j==1) then; cl = cr
        else if(c_type(i,cl)<0) then; cl = cr; end if
        if(j==nyy) then; cr = cl
        else if(c_type(i,cr)<0) then; cr = cl; end if
      end if

      call assertn(c_type(i,cl)>=0 .and. c_type(i,cr)>=0, "CalcW_Y. Ошибка сетки", i, j, -1)
      fy_wb0(i,j) = (c_wb0(i,cl) + c_wb0(i,cr)) / 2.          ! вертикальная скорость на дне t[n]
      fy_wt0(i,j) = (c_wt0(i,cl) + c_wt0(i,cr)) / 2.          ! вертикальная скорость на поверхности t[n]
      fyn_wb0(i,j) = (cn_wb0(i,cl) + cn_wb0(i,cr)) / 2.       ! вертикальная скорость на дне t[n+1]
      fyn_wt0(i,j) = (cn_wt0(i,cl) + cn_wt0(i,cr)) / 2.       ! вертикальная скорость на поверхности t[n+1]
    end do
  end do

end subroutine
