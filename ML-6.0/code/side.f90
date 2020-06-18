! Фаза 2 на X-гранях
subroutine Phase2X
  use variables
  implicit none

  interface
    subroutine TransportInvX(itype, i, j, k, off, inv); integer(4):: itype, i, j, k, off; real(R8):: inv; end subroutine TransportInvX
    subroutine TransportHX(i, j, off, inv); integer(4):: i, j, off; real(R8):: inv; end
  end interface

  integer(4):: i, j, k, il, ir, is, ic, ns
  real(R8):: invR, invQ, invU, invW, invD                             ! новые инварианты
  real(R8):: soundr, soundl, sound, ul, ur, uf, M, Gr, Gq
  real(R8):: drhoi, dtetai, dui, dvi, dwi
  real(R8):: w0i, zi, z0
  real(R8):: H0l, H0r, Hl, Hr, invH
  logical:: isFirstOut

  ! внутренние грани:
  ns = nfx_sides(BC_INNER)                                      ! количество внутренних граней
  do is=1, ns
    i = fx_sides(BC_INNER).ptr(is).i                            ! индексы грани на сетке
    j = fx_sides(BC_INNER).ptr(is).j
    call assertn(i>1 .and. i<nxx, "Phase2_X. Ошибка сетки (BC_INNER)", i, j, -1)

#if DLEV>0
    if(IsDebSideX(i, j, -1)) then
      write(17,"(/'Phase 2-X',2(';',i0),';;inner')") i, j
    endif
#endif

    il = i-1
    ir = i

    isFirstOut = .true.
    do k=1, ncz

      ! скорость на грани:
      !uf = fx_u0(i,j) + fx_du(i,j,k)                    ! скорость на грани
      ul = cs_u(il,j,k) !c_GetU(il,j,k,1)
      ur = cs_u(ir,j,k) !c_GetU(ir,j,k,1)
      uf=(ul + ur) / 2.
      M = uf / cs_sound                                 ! число Маха для переноса инвариантов

      if(abs(M)>0.9) then
        write(*,*) "ERROR in Phase2X: |M|>1. side (i,j,k):", i, j, k
        call avost
      endif

      ! инвариант R:
      call TransportInvX(T_INVR, i, j, k, DIRP, invR)              ! перенос через левую ячейку вправо
      Gr = cs_G(il,j,k)

      ! инвариант Q:
      call TransportInvX(T_INVQ, i, j, k, DIRM, invQ)              ! перенос через правую ячейку влево
      Gq = -cs_G(ir,j,k)

      ! инварианты U, W и Rho:
!      if(M>eps) then
      if(M>0) then
        call TransportInvX(T_INVU, i, j, k, DIRP, invU)      ! перенос через левую ячейку вправо
        call TransportInvX(T_INVW, i, j, k, DIRP, invW)      ! перенос через левую ячейку вправо
        call TransportInvX(T_INVD, i, j, k, DIRP, invD)      ! перенос через левую ячейку вправо
!      else if(M<-eps) then
      else if(M<=0) then
        call TransportInvX(T_INVU, i, j, k, DIRM, invU)      ! перенос через правую ячейку влево
        call TransportInvX(T_INVW, i, j, k, DIRM, invW)      ! перенос через правую ячейку влево
        call TransportInvX(T_INVD, i, j, k, DIRM, invD)      ! перенос через правую ячейку влево
      else
        invU = (cs_dv(il,j,k) + cs_dv(ir,j,k)) / 2.
        invW = (cs_dw(il,j,k) + cs_dw(ir,j,k)) / 2.
        invD = (cs_drho(il,j,k) + cs_drho(ir,j,k)) / 2.
      end if

      ! восстанавливаем значения из инвариантов:
      drhoi = invD                                        ! дельта плотности
      dvi = invU                                          ! дельта тангенцальной скорости
      dwi = invW                                          ! дельта вертикальной скорости

      ! invR = du + Gr * dteta
      ! invQ = du + Gq * dteta
      !----------------------------------
      ! invR - invQ = (Gr - Gq) * dteta
      ! invR*Gq - invQ*Gr = (Gq - Gr) * du

      dtetai = (invR - invQ) / (Gr - Gq)                  ! дельта псевдоплотности
      dui = (invR*Gq - invQ*Gr) / (Gq - Gr)               ! дельта нормальной скорости

!      w0i = fxn_GetW0(i,j,k)
      w0i = fxn_w0(i,j,k)

      fxn_dteta(i,j,k) = dtetai                           ! псевдоплотность
      fxn_du(i,j,k) = dui                                 ! нормальная скорость
      fxn_dv(i,j,k) = dvi                                 ! тангенцальная скорость
      fxn_dw(i,j,k) = dwi                                 ! вертикальная скорость
      fxn_drho(i,j,k) = drhoi                             ! плотность

#if DLEV>0
      if(IsDebSideX(i, j, k)) then
        !if(isFirstOut) &
          write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
        write(17,"('side:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
              teta0+fx_dteta(i,j,k), rho0+fx_drho(i,j,k), fx_u0(i,j)+fx_du(i,j,k), fx_v0(i,j)+fx_dv(i,j,k), fx_w0(i,j,k)+fx_dw(i,j,k), &
              invR, invQ, invU, invW, invD, Gr, Gq
        write(17,"('UFN:;;;',100(';',1p,g0.16))") &
              teta0+fxn_dteta(i,j,k), rho0+fxn_drho(i,j,k), fxn_u0(i,j)+fxn_du(i,j,k), fxn_v0(i,j)+fxn_dv(i,j,k), fxn_w0(i,j,k)+fxn_dw(i,j,k)
        isFirstOut = .false.
      end if
#endif

    end do      ! k

  end do   ! i,j

  !.. граничные условия ...................................................................

  call fx_BndPeriodic                                           ! периодические ГУ
  call fx_BndSlide                                              ! ГУ "стенка с проскальзыванием"
  call fx_BndStick                                              ! ГУ "стенка с прилипанием"
  call fx_BndInoutFix                                           ! ГУ "постоянный вток"
  call fx_BndInoutVarT                                          ! ГУ "переменный по времени вток"
  call fx_BndInoutVarZ                                          ! ГУ "переменный по Z вток"
  call fx_BndFixU                                               ! ГУ "фиксированная скорость"

end subroutine Phase2X

!------------------------------------------------------------------------------------------------------------------------------------

! Фаза 2 на Y-гранях
subroutine Phase2Y
  use variables
  implicit none

  interface
    subroutine TransportInvY(itype, i, j, k, off, inv); integer(4):: itype, i, j, k, off; real(R8):: inv; end
    subroutine TransportHY(i, j, off, inv); integer(4):: i, j, off; real(R8):: inv; end
  end interface

  integer(4):: i, j, k, jl, jr, is, ns
  real(R8):: invR, invQ, invU, invW, invD, invH                             ! новые инварианты
  real(R8):: soundr, soundl, sound, ul, ur, uf, M, Gr, Gq
  real(R8):: drhoi, dtetai, dui, dvi, dwi
  real(R8) :: w0i, zi, Hl, Hr, H0l, H0r, z0
  logical:: isOwnRight                                          ! признак "ячейка справа"
  logical:: isFirstOut

  ! внутренние грани:
  ns = nfy_sides(BC_INNER)                                      ! количество внутренних граней
  do is=1, ns
    i = fy_sides(BC_INNER).ptr(is).i                            ! индексы грани на сетке
    j = fy_sides(BC_INNER).ptr(is).j

    if(IsDebSideY(i, j, -1)) then
      write(17,"(/'Phase 2-Y',2(';',i0),';inner')") i, j
    endif

    jl = j-1
    jr = j

    isFirstOut = .true.
    do k=1, ncz

      ! скорость на грани:
      !uf = fy_v0(i,j) + fy_dv(i,j,k)                    ! скорость на грани
      ul = cs_v(i,jl,k)
      ur = cs_v(i,jr,k)
      uf=(ul + ur) / 2.
      M = uf / cs_sound                                 ! число Маха для переноса инвариантов

      if(abs(M)>0.9) then
        write(*,*) "ERROR in Phase2Y: |M|>1. side (i,j,k):", i, j, k
        call avost
      endif

      ! инвариант R:
      call TransportInvY(T_INVR, i, j, k, DIRP, invR)              ! перенос через левую ячейку вправо
      Gr = cs_G(i,jl,k)

      ! инвариант Q:
      call TransportInvY(T_INVQ, i, j, k, DIRM, invQ)              ! перенос через правую ячейку влево
      Gq = -cs_G(i,jr,k)

      ! инварианты U, W и Rho:
!      if(M>eps) then
      if(M>0) then
        call TransportInvY(T_INVU, i, j, k, DIRP, invU)      ! перенос через левую ячейку вправо
        call TransportInvY(T_INVW, i, j, k, DIRP, invW)      ! перенос через левую ячейку вправо
        call TransportInvY(T_INVD, i, j, k, DIRP, invD)      ! перенос через левую ячейку вправо
!      else if(M<-eps) then
      else if(M<=0) then
        call TransportInvY(T_INVU, i, j, k, DIRM, invU)      ! перенос через правую ячейку влево
        call TransportInvY(T_INVW, i, j, k, DIRM, invW)      ! перенос через правую ячейку влево
        call TransportInvY(T_INVD, i, j, k, DIRM, invD)      ! перенос через правую ячейку влево
      else
        invU = (cs_du(i,jl,k) + cs_du(i,jr,k)) / 2.
        invW = (cs_dw(i,jl,k) + cs_dw(i,jr,k)) / 2.
        invD = (cs_drho(i,jl,k) + cs_drho(i,jr,k)) / 2.
      end if

      ! восстанавливаем значения из инвариантов:
      drhoi = invD                                        ! дельта плотности
      dui = invU                                          ! дельта тангенцальной скорости
      dwi = invW                                          ! дельта вертикальной скорости

      ! invR = du + Gr * dteta
      ! invQ = du + Gq * dteta
      !----------------------------------
      ! invR - invQ = (Gr - Gq) * dteta
      ! invR*Gq - invQ*Gr = (Gq - Gr) * du

      dtetai = (invR - invQ) / (Gr - Gq)                  ! дельта псевдоплотности
      dvi = (invR*Gq - invQ*Gr) / (Gq - Gr)               ! дельта нормальной скорости

      w0i = fyn_w0(i,j,k)

      fyn_dteta(i,j,k) = dtetai
      fyn_du(i,j,k) = dui
      fyn_dv(i,j,k) = dvi                            ! тангенцальная скорость
      fyn_dw(i,j,k) = dwi                            ! вертикальная скорость
      fyn_drho(i,j,k) = drhoi                ! плотность

      if(IsDebSideY(i, j, k)) then
        if(isFirstOut) &
          write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
        write(17,"('side:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
              teta0+fy_dteta(i,j,k), rho0+fy_drho(i,j,k), fy_u0(i,j)+fy_du(i,j,k), fy_v0(i,j)+fy_dv(i,j,k), fy_w0(i,j,k)+fy_dw(i,j,k), &
              invR, invQ, invU, invW, invD, Gr, Gq
        write(17,"('UFN:;;;',100(';',1p,g0.16))") &
              teta0+fyn_dteta(i,j,k), rho0+fyn_drho(i,j,k), fyn_u0(i,j)+fyn_du(i,j,k), fyn_v0(i,j)+fyn_dv(i,j,k), fyn_w0(i,j,k)+fyn_dw(i,j,k)
        isFirstOut = .false.
      end if

    end do      ! k

  end do   ! i,j

  !.. граничные условия ...................................................................

  call fy_BndPeriodic                                           ! периодические ГУ
  call fy_BndSlide                                              ! ГУ "стенка с проскальзыванием"
  call fy_BndStick                                              ! ГУ "стенка с прилипанием"
  call fy_BndInoutFix                                           ! ГУ "постоянный вток"
  call fy_BndInoutVarT                                          ! ГУ "переменный по времени вток"
  !call fy_BndInoutVarZ                                          ! ГУ "переменный по Z вток"
  call fy_BndFixU                                               ! ГУ "фиксированная скорость"

end

!------------------------------------------------------------------------------------------------------------------------------------

! Первая часть фазы 2-Z. На горизонтальных гранях (кроме верхней) определяются величины fzn_du, fzn_dv, fzn_dw, fzn_dteta, fzn_drho.
! На верхней грани определяются только транспортные инварианты и величины: fzn_du, fzn_dv, fzn_drho.
subroutine Phase2Z_in
  use variables
  implicit none

  interface
    subroutine TransportInvZ(itype, i, j, k, off, inv); integer(4):: itype, i, j, k, off; real(R8):: inv; end
  end interface

  integer(4):: i, j, k, jl, jr, ic, nc
  real(R8):: invR, invQ, invU, invW, invD                             ! новые инварианты
  real(R8):: soundr, soundl, sound, uf, ul, ur, M, Gr, Gq
  real(R8):: drhoi, dtetai, dui, dvi, dwi
  real(R8):: w0i, zi
  real(R8):: bxp, bxm, byp, bym, uc, vc, wc
  logical:: isOwnRight                                          ! признак "ячейка справа"
  logical:: isFirstOut

  isFirstOut = .true.

  nc = ncells                                                  ! количество рабочих ячеек
  do ic=1, nc
    i = cells.ptr(ic).i                                         ! индексы ячеек на сетке
    j = cells.ptr(ic).j

    if(IsDebSideZ(i, j, -1)) then
      write(17,"(/'Phase 2-Z',2(';',i0))") i, j
    endif

    isFirstOut = .true.
    do k=2, nz-1                                                ! внутренние горизонтальные грани

      !uf = fz_GetW0(i,j,k) + fz_dw(i,j,k)                        ! скорость на грани
      uf = fz_w0(i,j,k) + fz_dw(i,j,k)                        ! скорость на грани
      ul = cs_w0(i,j,k-1) + cs_dw(i,j,k-1)
      ur = cs_w0(i,j,k) + cs_dw(i,j,k)
      uf = (ul + ur) / 2. - fz_zp(i,j,k)
      M = uf / cs_sound                                         ! число Маха

      if(abs(M)>0.9) then
        write(*,*) "ERROR in Phase2Z: |M|>1. side (i,j,k):", i, j, k
        call avost
      endif

      ! инвариант R (снизу вверх):
      call TransportInvZ(T_INVR, i, j, k, DIRM, invR)              ! перенос через нижнюю ячейку вверх
      Gr = cs_G(i,j,k)

      ! инвариант Q(сверху вниз):
      call TransportInvZ(T_INVQ, i, j, k, DIRP, invQ)              ! перенос через верхнюю ячейку вниз
      Gq = -cs_G(i,j,k-1)

      ! инварианты U, W и Rho:
      !if(M>eps) then                                           ! скорость вверх
      if(M>0.) then                                           ! скорость вверх
        call TransportInvZ(T_INVU, i, j, k, DIRM, invU)         ! перенос через нижнюю ячейку вверх
        call TransportInvZ(T_INVW, i, j, k, DIRM, invW)
        call TransportInvZ(T_INVD, i, j, k, DIRM, invD)
      !else if(M<-eps) then
      else if(M<=0) then
        call TransportInvZ(T_INVU, i, j, k, DIRP, invU)      ! перенос через верхнюю ячейку вниз
        call TransportInvZ(T_INVW, i, j, k, DIRP, invW)
        call TransportInvZ(T_INVD, i, j, k, DIRP, invD)
      else
        invU = (cs_du(i,j,k-1) + cs_du(i,j,k)) / 2.
        invW = (cs_dv(i,j,k-1) + cs_dv(i,j,k)) / 2.
        invD = (cs_drho(i,j,k-1) + cs_drho(i,j,k)) / 2.
      end if

      ! восстанавливаем значения из инвариантов:
      drhoi = invD                                        ! дельта плотности
      dui = invU                                          ! дельта тангенцальной скорости U
      dvi = invW                                          ! дельта тангенцальной скорости V

      ! invR = du + Gr * dteta
      ! invQ = du + Gq * dteta
      !----------------------------------
      ! invR - invQ = (Gr - Gq) * dteta
      ! invR*Gq - invQ*Gr = (Gq - Gr) * du

      dtetai = (invR - invQ) / (Gr - Gq)                  ! дельта псевдоплотности
      dwi = (invR*Gq - invQ*Gr) / (Gq - Gr)               ! дельта нормальной скорости

      fzn_dteta(i,j,k) = dtetai
      fzn_du(i,j,k) = dui
      fzn_dv(i,j,k) = dvi                            ! тангенцальная скорость
      fzn_dw(i,j,k) = dwi                            ! вертикальная скорость
      fzn_drho(i,j,k) = drhoi                ! плотность

      ! вычисление абсолютных переменных:
      fzn_teta(i,j,k) = dtetai + teta0
      fzn_rho(i,j,k) = drhoi + rho0
      fzn_u(i,j,k) = dui + cn_u0(i,j)
      fzn_v(i,j,k) = dvi + cn_v0(i,j)

      if(IsDebSideZ(i, j, k)) then
        !if(isFirstOut) &
          write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
        write(17,"('side:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
              fzn_dteta(i,j,k), fzn_drho(i,j,k), fzn_du(i,j,k), fzn_dv(i,j,k), fzn_dw(i,j,k), &
              invR, invQ, invU, invW, invD, Gr, Gq
        write(17,"('UF:;;;',100(';',1p,g0.16))") &
              fz_teta(i,j,k), fz_rho(i,j,k), fz_u(i,j,k), fz_v(i,j,k)
        write(17,"('UFN:;;;',100(';',1p,g0.16))") &
              fzn_teta(i,j,k), fzn_rho(i,j,k), fzn_u(i,j,k), fzn_v(i,j,k)
        isFirstOut = .false.
      end if

    end do      ! k

    !-- отдельно для дна (горизонтальная грань k=nz): --------------------------------------

    k = nz

    uf = cs_w0(i,j,k-1) + cs_dw(i,j,k-1)                        ! скорость в придонной ячейке
    M = uf / cs_sound                                         ! число Маха

    ! инвариант Q(сверху вниз):
    call TransportInvZ(T_INVQ, i, j, k, DIRP, invQ)              ! перенос через верхнюю ячейку вниз
    Gq = -cs_G(i,j,k-1)

    ! инварианты U, W и Rho:
    if(M<-eps) then
      call TransportInvZ(T_INVU, i, j, k, DIRP, invU)             ! перенос через верхнюю ячейку вниз
      call TransportInvZ(T_INVW, i, j, k, DIRP, invW)
      call TransportInvZ(T_INVD, i, j, k, DIRP, invD)
    else
      invU = cs_du(i,j,ncz)
      invW = cs_dv(i,j,ncz)
      invD = cs_drho(i,j,ncz)
    end if

    ! восстанавливаем значения из инвариантов:
    drhoi = invD                                        ! дельта плотности
    dui = invU                                          ! дельта тангенцальной скорости U
    dvi = invW                                          ! дельта тангенцальной скорости V

    ! dB       dB
    ! -- * u + -- * v = w           -- ГУ для нижней границы
    ! dx       dy
    !
    ! Q = dw + G * dteta

    bxp = fxn_z(i+1,j,nz)
    bxm = fxn_z(i,j,nz)
    byp = fyn_z(i,j+1,nz)
    bym = fyn_z(i,j,nz)
    uc = cn_u0(i,j) + dui
    vc = cn_v0(i,j) + dvi

    wc = uc * (bxp - bxm) / c_dx(i) + vc * (byp - bym) / c_dy(j)
    dwi = wc - cn_wb0(i,j)
    dtetai = (invQ - dwi) / Gq

    fzn_dteta(i,j,k) = dtetai
    fzn_du(i,j,k) = dui
    fzn_dv(i,j,k) = dvi                            ! тангенцальная скорость
    fzn_dw(i,j,k) = dwi                            ! вертикальная скорость
    fzn_drho(i,j,k) = drhoi                ! плотность

    ! вычисление абсолютных переменных:
    fzn_teta(i,j,k) = dtetai + teta0
    fzn_rho(i,j,k) = drhoi + rho0
    fzn_u(i,j,k) = dui + cn_u0(i,j)
    fzn_v(i,j,k) = dvi + cn_v0(i,j)

    if(IsDebSideZ(i, j, k)) then
!      if(isFirstOut) &
        write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
      write(17,"('side:',3(';',i0),100(';',1p,g0.16))") &
          i, j, k, &
            fzn_dteta(i,j,k), fzn_drho(i,j,k), fzn_du(i,j,k), fzn_dv(i,j,k), fzn_dw(i,j,k), &
            invR, invQ, invU, invW, invD, Gr, Gq
      write(17,"('UF:;;;',100(';',1p,g0.16))") &
            fz_teta(i,j,k), fz_rho(i,j,k), fz_u(i,j,k), fz_v(i,j,k)
      write(17,"('UFN:;;;',100(';',1p,g0.16))") &
            fzn_teta(i,j,k), fzn_rho(i,j,k), fzn_u(i,j,k), fzn_v(i,j,k)
      isFirstOut = .false.
    end if

    !-- для верхней границы только трансплртные инварианты: --------------------------------

    k = 1

    ! инварианты U, W и Rho:
    if(M>eps) then                                           ! скорость вверх
      call TransportInvZ(T_INVU, i, j, 1, DIRM, invU)         ! перенос через нижнюю ячейку вверх
      call TransportInvZ(T_INVW, i, j, 1, DIRM, invW)
      call TransportInvZ(T_INVD, i, j, 1, DIRM, invD)
    else
      invU = cs_du(i,j,1)
      invW = cs_dv(i,j,1)
      invD = cs_drho(i,j,1)
    end if

    drhoi = invD                                        ! дельта плотности
    dui = invU                                          ! дельта тангенцальной скорости U
    dvi = invW                                          ! дельта тангенцальной скорости V

    fzn_du(i,j,k) = dui
    fzn_dv(i,j,k) = dvi                            ! тангенцальная скорость
    fzn_drho(i,j,k) = drhoi                ! плотность

    ! вычисление абсолютных переменных:
    fzn_rho(i,j,k) = drhoi + rho0
    fzn_u(i,j,k) = dui + cn_u0(i,j)
    fzn_v(i,j,k) = dvi + cn_v0(i,j)

    if(IsDebSideZ(i, j, k)) then
      !if(isFirstOut) &
        write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
        write(17,"('side:',3(';',i0),';',3(';',1p,g0.16),';',100(';',1p,g0.16))") &
            i, j, k, &
              fzn_drho(i,j,k), fzn_du(i,j,k), fzn_dv(i,j,k), &
              invR, invQ, invU, invW, invD, Gr, Gq
        write(17,"('UF:;;;;',100(';',1p,g0.16))") &
              fz_rho(i,j,k), fz_u(i,j,k), fz_v(i,j,k)
        write(17,"('UFN:;;;;',100(';',1p,g0.16))") &
              fzn_rho(i,j,k), fzn_u(i,j,k), fzn_v(i,j,k)
      isFirstOut = .false.
    end if

  end do   ! i,j

end

!-----------------------------------------------------------------------------------------------------------------------------------------------

! Вторая часть фазы 2-Z, определение верхней границы. Определяются величины fzn_dteta(1), fzn_dw(1)
subroutine Phase2Z_top
  use variables
  implicit none

  integer:: i, j, k, il, ir, jl, jr, nc, ic, ibc
  real(R8):: z0, zp, Gr, QQ, uxl, uxr, vyl, vyr, uc, vc, uf, M, H0l, H0r, Hl, Hr, w0, dzdx, dzdy, dwdz, bcU, h
  real(R8):: invR, invH
  logical:: isOwnRight

  !-- построение z-сетки на X-гранях: ---------------------------------------

  ! перенос dH:
  do i=1,nxx
    do j=1,nxy
      if(fx_type(i,j)<0) cycle                                  ! пропускаем не рабочие грани
      select case(fx_type(i,j))                                 ! в зависимости от типа ГУ

        case(BC_INNER) !.. внутренние грани ....................................................

          il = i-1
          ir = i
          ! скорость на грани:
          !uf = fx_u0(i,j) + fx_du(i,j,1)                      ! X-скорость на поверхности
          uf = (fzn_du(il,j,1)+c_u0(il,j) + fzn_du(ir,j,1)+c_u0(ir,j)) / 2.
          M = uf / cs_sound                                     ! число Маха для переноса инвариантов

          if(abs(M)>0.9) then
            write(*,*) "ERROR in CalcH_X: |M|>1. side (i,j):", i, j
            call avost
          endif

          ! перенос H:
          if(M>eps) then
            call TransportHX(i, j, DIRP, invH)                  ! перенос через левую ячейку вправо
          else if(M<-eps) then
            call TransportHX(i, j, DIRM, invH)                  ! перенос через правую ячейку влево
          else
            H0l = c_b(il,j) + cs_h0(il,j)                       ! поверхность МВ в ячейке слева
            H0r = c_b(ir,j) + cs_h0(ir,j)                       ! поверхность МВ в ячейке справа
            Hl = fzs_z(il,j,1)                                  ! реальная поверхность в ячейке слева
            Hr = fzs_z(ir,j,1)                                  ! реальная поверхность в ячейке справа
            invH = ((Hl-H0l) + (Hr-H0r)) / 2.
          end if
          fxn_dH(i,j) = invH

        case(BC_PERIODIC) !.. периодические ГУ ...........................................................

          if(i/=1) cycle                                        ! только для левой границы
          il = ncx                                              ! ячейка слева (крайняя правая)
          ir = 1                                                ! ячейка справа
          ! скорость на грани:
          uf = fx_u0(i,j) + fx_du(i,j,1)                        ! X-скорость на поверхности
          M = uf / cs_sound                                     ! число Маха для переноса инвариантов

          if(abs(M)>0.9) then
            write(*,*) "ERROR in CalcH_X: |M|>1. side (i,j):", i, j
            call avost
          endif

          ! перенос H:
          if(M>eps) then
            call TransportHX(nxx, j, DIRP, invH)                ! перенос через левую ячейку вправо
          else if(M<-eps) then
            call TransportHX(1, j, DIRM, invH)                  ! перенос через правую ячейку влево
          else
            H0l = c_b(il,j) + cs_h0(il,j)                       ! поверхность МВ в ячейке слева
            H0r = c_b(ir,j) + cs_h0(ir,j)                       ! поверхность МВ в ячейке справа
            Hl = fzs_z(il,j,1)                                  ! реальная поверхность в ячейке слева
            Hr = fzs_z(ir,j,1)                                  ! реальная поверхность в ячейке справа
            invH = ((Hl-H0l) + (Hr-H0r)) / 2.
          end if
          fxn_dH(1,j) = invH
          fxn_dH(nxx,j) = invH

        case(BC_STICK, BC_SLIDE) !.. стена ...................................................................

          il = i-1
          ir = i
          isOwnRight = fx_IsOwnRight(i,j)

          ! скорость на грани:
          if(isOwnRight) then                                   ! ячейка справа
            uf = fzn_du(ir,j,1)+c_u0(ir,j)
          else
            uf = fzn_du(il,j,1)+c_u0(il,j)
          end if
          M = uf / cs_sound                                     ! число Маха для переноса инвариантов

          if(isOwnRight) then                                   ! ячейка справа                         | O

            if(M<-eps) then                                     ! течение справа налево (через ячейку)  |<O-
              call TransportHX(i, j, DIRM, invH)                ! перенос через правую ячейку влево
            else
              H0r = c_b(ir,j) + cs_h0(ir,j)                     ! поверхность МВ в ячейке справа
              Hr = fzs_z(ir,j,1)                                ! реальная поверхность в ячейке справа
              invH = Hr - H0r
            end if

          else                                                  ! ячейка слева                           O |

            if(M>eps) then                                      ! течение слева направо (через ячейку)  -O>|
              call TransportHX(i, j, DIRP, invH)                ! перенос через левую ячейку вправо
            else
              H0l = c_b(il,j) + cs_h0(il,j)                     ! поверхность МВ в ячейке слева
              Hl = fzs_z(il,j,1)                                ! реальная поверхность в ячейке слева
              invH = Hl - H0l
            end if

          end if

          fxn_dH(i,j) = invH

        case(BC_IN, BC_IN_T, BC_IN_Z) !.. вход ...........................................................

          il = i-1
          ir = i

          ibc = fx_bc(i,j)
          isOwnRight = fx_IsOwnRight(i,j)

          if(taskNum==4) then
            bcU = task4_u0(x(i))
          else
            bcU = bcData(ibc).u
          endif

          ! скорость на грани:
          !uf = bcData(ibc).u
          if(isOwnRight) then                                   ! ячейка справа
            uf = (bcU + fzn_du(ir,j,1)+c_u0(ir,j)) / 2.
          else
            uf = (fzn_du(il,j,1)+c_u0(il,j) + bcU) / 2.
          end if

          M = uf / cs_sound                                     ! число Маха для переноса инвариантов

          if(abs(M)>0.9) then
            write(*,*) "ERROR in CalcH_X: |M|>1. side (i,j):", i, j
            call avost
          endif

          ! перенос H:
          if(isOwnRight) then                                   ! ячейка справа

            if(M>eps) then                                      ! течение слева направо (от границы)
              invH = 0.
            else if(M<-eps) then                                ! течение справа налево (через ячейку)
              call TransportHX(i, j, DIRM, invH)                ! перенос через правую ячейку влево
            else if(M>eps) then
              invH = 0.                                         ! за границей h и h0 совпадают
            else
              H0r = c_b(ir,j) + cs_h0(ir,j)                     ! поверхность МВ в ячейке справа
              Hr = fzs_z(ir,j,1)                                ! реальная поверхность в ячейке справа
              invH = (0. + (Hr-H0r)) / 2.
            end if

          else                                                  ! ячейка слева

            if(M>eps) then                                      ! течение слева направо (через ячейку)
              call TransportHX(i, j, DIRP, invH)                ! перенос через левую ячейку вправо
            else if(M<-eps) then                                ! течение справа налево (от границы)
              invH = 0.
            else
              H0l = c_b(il,j) + cs_h0(il,j)                     ! поверхность МВ в ячейке слева
              Hl = fzs_z(il,j,1)                                ! реальная поверхность в ячейке слева
              invH = ((Hl-H0l) + 0.) / 2.
            end if

          end if

          fxn_dH(i,j) = invH

        case default !.. неизвестный тип ГУ ..........................................................

          call avost("Не рализовано. "//__FILE__, __LINE__)

      end select
    end do
  end do

  do i=1,nxx
    do j=1,nxy
      if(fx_type(i,j)<0) cycle                                  ! пропускаем не рабочие грани
      ! Выше на каждой рабочей грани вычислен уровень поверхности.
      ! Определяем z-координаты разделов между слоями
      z0 = fxn_z0(i,j) + fxn_dH(i,j)                            ! реальная поверхность
      h = z0 - fx_b(i,j)                                        ! текущая глубина
      fxn_z(i,j,1) = z0
      do k=2,nz-1
        fxn_z(i,j,k) = fxn_z(i,j,k-1) - h * cfh(k-1)
      end do
      fxn_z(i,j,nz) = fx_b(i,j)

      ! после построения сетки можно вычислить абсолютные значения потоковых переменных:
      do k=1,ncz
        fxn_w0(i,j,k) = fxn_GetW0(i,j,k)                        ! Вертикальная скорость МВ на X-гранях
        fxn_teta(i,j,k) = teta0       + fxn_dteta(i,j,k)
        fxn_rho(i,j,k)  = rho0        + fxn_drho(i,j,k)
        fxn_u(i,j,k)    = fxn_u0(i,j) + fxn_du(i,j,k)
        fxn_v(i,j,k)    = fxn_v0(i,j) + fxn_dv(i,j,k)
        fxn_w(i,j,k)    = fxn_w0(i,j,k) + fxn_dw(i,j,k)
      end do

      if(IsDebSideX(i, j, -1)) then
        write(17,"(/a,1000(';',i0))") "Phase 2-Z-top-X;i;j;;fxn_dH;fxn_z0;"
        write(17,"(2(';',i0),';',2(';',1p,g24.16),a,1000(';',1p,g24.16))") i, j, fxn_dH(i,j), fxn_z0(i,j)
        write(17,"(a)") ";;;k;fxn_z;fxn_w0;fxn_teta;fxn_rho;fxn_u;fxn_v;fxn_w"
        do k=1,ncz
          if(IsDebSideX(i,j,k)) then
            write(17,"('UFN:;;;',i0,100(';',1p,g0.16))") k, &
                  fxn_z(i,j,k), fxn_w0(i,j,k), fxn_teta(i,j,k), fxn_rho(i,j,k), fxn_u(i,j,k), fxn_v(i,j,k), fxn_w(i,j,k)
          end if
        end do
      end if
    end do
  end do

  !-- построение z-сетки на Y-гранях: ---------------------------------------

  ! перенос dH:
  do i=1,nyx
    do j=1,nyy
      if(fy_type(i,j)<0) cycle                                  ! пропускаем не рабочие грани
      select case(fy_type(i,j))                                 ! в зависимости от типа ГУ

        case(BC_INNER) !.. внутренние грани ....................................................
          jl = j-1
          jr = j
          ! скорость на грани:
          uf = fy_v0(i,j) + fy_dv(i,j,1)                        ! Y-скорость на поверхности
          M = uf / cs_sound                                     ! число Маха для переноса инвариантов

          if(abs(M)>0.9) then
            write(*,*) "ERROR in CalcH_Y (inner): |M|>1. side (i,j):", i, j
            call avost
          endif

          ! перенос H:
          if(M>eps) then
            call TransportHY(i, j, DIRP, invH)                  ! перенос через левую ячейку вправо
          else if(M<-eps) then
            call TransportHY(i, j, DIRM, invH)                  ! перенос через правую ячейку влево
          else
            H0l = c_b(i,jl) + cs_h0(i,jl)                       ! поверхность МВ в ячейке слева
            H0r = c_b(i,jr) + cs_h0(i,jr)                       ! поверхность МВ в ячейке справа
            Hl = fzs_z(i,jr,1)                                  ! реальная поверхность в ячейке слева
            Hr = fzs_z(i,jr,1)                                  ! реальная поверхность в ячейке справа
            invH = ((Hl-H0l) + (Hr-H0r)) / 2.
          end if

          fyn_dH(i,j) = invH

        case(BC_PERIODIC) !.. периодические ГУ ...........................................................
          if(j/=1) cycle                                        ! только для левой границы
          jl = ncy                                              ! ячейка слева (крайняя правая)
          jr = 1                                                ! ячейка справа
          ! скорость на грани:
          uf = fy_v0(i,j) + fy_dv(i,j,1)                        ! Y-скорость на поверхности
          M = uf / cs_sound                                     ! число Маха для переноса инвариантов

          if(abs(M)>0.9) then
            write(*,*) "ERROR in CalcH_Y (periodic): |M|>1. side (i,j):", i, j
            call avost
          endif

          ! перенос H:
          if(M>eps) then
            call TransportHY(i, nyy, DIRP, invH)                ! перенос через левую ячейку вправо
          else if(M<-eps) then
            call TransportHY(i, 1, DIRM, invH)                  ! перенос через правую ячейку влево
          else
            H0l = c_b(i,jl) + cs_h0(i,jl)                       ! поверхность МВ в ячейке слева
            H0r = c_b(i,jr) + cs_h0(i,jr)                       ! поверхность МВ в ячейке справа
            Hl = fzs_z(i,jr,1)                                  ! реальная поверхность в ячейке слева
            Hr = fzs_z(i,jr,1)                                  ! реальная поверхность в ячейке справа
            invH = ((Hl-H0l) + (Hr-H0r)) / 2.
          end if

          fyn_dH(i,j) = invH
          fyn_dH(i,nyy) = invH

        case(BC_STICK, BC_SLIDE) !.. стена ...................................................................

          jl = j-1
          jr = j
          isOwnRight = fy_IsOwnRight(i,j)

          ! скорость на грани:
          if(isOwnRight) then                                   ! ячейка справа
            uf = fzn_dv(i,jr,1)+c_v0(i,jr)
          else
            uf = fzn_dv(i,jl,1)+c_v0(i,jl)
          end if
          M = uf / cs_sound                                     ! число Маха для переноса инвариантов

          if(isOwnRight) then                                   ! ячейка справа                         | O

            if(M<-eps) then                                     ! течение справа налево (через ячейку)  |<O-
              call TransportHY(i, j, DIRM, invH)                ! перенос через правую ячейку влево
            else
              H0r = c_b(i,jr) + cs_h0(i,jr)                     ! поверхность МВ в ячейке справа
              Hr = fzs_z(i,jr,1)                                ! реальная поверхность в ячейке справа
              invH = Hr - H0r
            end if

          else                                                  ! ячейка слева                           O |

            if(M>eps) then                                      ! течение слева направо (через ячейку)  -O>|
              call TransportHY(i, j, DIRP, invH)                ! перенос через левую ячейку вправо
            else
              H0l = c_b(i,jl) + cs_h0(i,jl)                     ! поверхность МВ в ячейке слева
              Hl = fzs_z(i,jl,1)                                ! реальная поверхность в ячейке слева
              invH = Hl - H0l
            end if

          end if

          fyn_dH(i,j) = invH

        case(BC_IN, BC_IN_T) !.. вход ...........................................................
          jl = j-1
          jr = j

          ibc = fy_bc(i,j)
          isOwnRight = fy_IsOwnRight(i,j)

          ! скорость на грани:
          if(isOwnRight) then                                   ! ячейка справа
            uf = (bcData(ibc).v + fzn_dv(i,jr,1)+c_v0(i,jr)) / 2.
          else
            uf = (fzn_dv(i,jl,1)+c_v0(i,jl) + bcData(ibc).v) / 2.
          end if
          M = uf / cs_sound                                     ! число Маха для переноса инвариантов

          if(abs(M)>0.9) then
            write(*,*) "ERROR in CalcH_Y: |M|>1. side (i,j):", i, j
            call avost
          endif

          ! перенос H:
          if(isOwnRight) then                                   ! ячейка справа

            if(M>eps) then                                      ! течение слева направо (от границы)
              invH = 0.
            else if(M<-eps) then                                ! течение справа налево (через ячейку)
              call TransportHY(i, j, DIRM, invH)                ! перенос через правую ячейку влево
            else
              H0r = c_b(i,jr) + cs_h0(i,jr)                     ! поверхность МВ в ячейке справа
              Hr = fzs_z(i,jr,1)                                ! реальная поверхность в ячейке справа
              invH = (0. + (Hr-H0r)) / 2.
            end if

          else                                                  ! ячейка слева

            if(M>eps) then                                      ! течение слева направо (через ячейку)
              call TransportHY(i, j, DIRP, invH)                ! перенос через левую ячейку вправо
            else if(M<-eps) then                                ! течение справа налево (от границы)
              invH = 0.
            else
              H0l = c_b(i,jl) + cs_h0(i,jl)                     ! поверхность МВ в ячейке слева
              Hl = fzs_z(i,jl,1)                                ! реальная поверхность в ячейке слева
              invH = ((Hl-H0l) + 0.) / 2.
            end if

          end if

          fyn_dH(i,j) = invH

        case default !.. неизвестный тип ГУ ..........................................................

          call avost("Не рализовано. "//__FILE__, __LINE__)

      end select
    end do
  end do

  do i=1,nyx
    do j=1,nyy
      if(fy_type(i,j)<0) cycle                                  ! пропускаем не рабочие грани
      ! Выше на каждой рабочей грани вычислен уровень поверхности.
      ! Определяем z-координаты разделов между слоями
      z0 = fyn_z0(i,j) + fyn_dH(i,j)                            ! реальная поверхность
      h = z0 - fy_b(i,j)                                        ! текущая глубина
      fyn_z(i,j,1) = z0
      do k=2,nz-1
        fyn_z(i,j,k) = fyn_z(i,j,k-1) - h * cfh(k-1)
      end do
      fyn_z(i,j,nz) = fy_b(i,j)

      ! после построения сетки можно вычислить абсолютные значения потоковых переменных:
      do k=1,ncz
        fyn_w0(i,j,k) = fyn_GetW0(i,j,k)                        ! Вертикальная скорость МВ на Y-гранях
        fyn_teta(i,j,k) = teta0       + fyn_dteta(i,j,k)
        fyn_rho(i,j,k)  = rho0        + fyn_drho(i,j,k)
        fyn_u(i,j,k)    = fyn_u0(i,j) + fyn_du(i,j,k)
        fyn_v(i,j,k)    = fyn_v0(i,j) + fyn_dv(i,j,k)
        fyn_w(i,j,k)    = fyn_w0(i,j,k) + fyn_dw(i,j,k)
      end do

      if(IsDebSideY(i, j, -1)) then
        write(17,"(a,1000(';',i0))") "/Phase 2-Z-top-Y;i;j;;fyn_dH;fyn_z0"
        write(17,"(2(';',i0),';',2(';',1p,g24.16),a,1000(';',1p,g24.16))") i, j, fyn_dH(i,j), fyn_z0(i,j)
        write(17,"(a)") ";;;k;fyn_z;fyn_w0;fyn_teta;fyn_rho;fyn_u;fyn_v;fyn_w"
        do k=1,ncz
          if(IsDebSideY(i,j,k)) then
            write(17,"('UFN:;;;',i0,100(';',1p,g0.16))") k, &
                  fyn_z(i,j,k), fyn_w0(i,j,k), fyn_teta(i,j,k), fyn_rho(i,j,k), fyn_u(i,j,k), fyn_v(i,j,k), fyn_w(i,j,k)
          end if
        end do
      end if

    end do
  end do

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  nc = ncells                                                   ! количество рабочих ячеек
  do ic=1, nc
    i = cells.ptr(ic).i                                         ! индексы ячеек на сетке
    j = cells.ptr(ic).j
    k = 1

    ir = i + 1                                                  ! индексы граней вокруг ячейки
    il = i
    jr = j+1
    jl = j

    if(IsDebSideZ(i, j, -1)) then
      write(17,"(/'Phase 2-Z-top-Z',2(';',i0))") i, j
    endif

    ! инвариант R (снизу вверх):
    call TransportInvZ(T_INVR, i, j, k, DIRM, invR)             ! перенос через нижнюю ячейку вверх
    Gr = cs_G(i,jl,k)

    QQ = rho0 * g / sound0_2

    z0 = c_b(i,j) + cn_h0(i,j)                                  ! уровень поверхности МВ t[n+1]
    uxr = fxn_u0(ir,j) + fxn_du(ir,j,1)
    uxl = fxn_u0(il,j) + fxn_du(il,j,1)
    vyr = fyn_v0(i,jr) + fyn_dv(i,jr,1)
    vyl = fyn_v0(i,jl) + fyn_dv(i,jl,1)
    uc = (uxl + uxr) / 2.
    vc = (vyl + vyr) / 2.

    w0 = (cn_wb0(i,j) * cn_z0(i,j) - cn_wt0(i,j) * c_b(i,j)) / cn_h0(i,j)
    dwdz = (cn_wt0(i,j) - cn_wb0(i,j)) / cn_h0(i,j)

    fzn_z(i,j,1) = ( fzs_z(i,j,1) + dt/2. * &
      ( &
        invR + Gr * QQ * z0 + &
        w0 - &
        uc * (fxn_z(ir,j,1) - fxn_z(il,j,1)) / c_dx(i) - &
        vc * (fyn_z(i,jr,1) - fyn_z(i,jl,1)) / c_dy(j) &
      ) &
                   ) / &
        (1. + dt/2. * (Gr * QQ - dwdz))

    fzn_dteta(i,j,1) = -QQ * &
        (z0 - fzn_z(i,j,1))

    fzn_dw(i,j,1) = invR - Gr * fzn_dteta(i,j,1)

    ! новая сетка по вертикали в ячейке (i,j):
    h = fzn_z(i,j,1) - c_b(i,j)
    do k=2,nz-1
      fzn_z(i,j,k) = fzn_z(i,j,k-1) - h * cfh(k-1)
    end do
    fzn_z(i,j,nz) = c_b(i,j)

    ! после построения сетки можно вычислить абсолютные значения потоковых переменных:
    do k=1,nz
      fzn_w0(i,j,k) = fzn_GetW0(i,j,k)                          ! Вертикальная скорость МВ на Z-гранях
      fzn_teta(i,j,k) = teta0       + fzn_dteta(i,j,k)
      fzn_w(i,j,k)    = fzn_w0(i,j,k) + fzn_dw(i,j,k)
      if(k/=nz) then
        cn_w0(i,j,k) = cn_GetW0(i,j,k)                          ! Вертикальная скорость МВ в слой-ячейках
      end if
    end do

    !fzn_zp(i,j,1) = (fzn_z(i,j,1) - fzs_z(i,j,1)) / (dt/2.)
    dzdx = (fxn_z(i+1,j,1) - fxn_z(i,j,1)) / c_dx(i)
    dzdy = (fyn_z(i,j+1,1) - fyn_z(i,j,1)) / c_dy(j)
    fzn_zp(i,j,1) = fzn_w(i,j,1) - (fzn_u(i,j,1) * dzdx + fzn_v(i,j,1) * dzdy)
    do k=2,nz-1
      fzn_zp(i,j,k) = fzn_zp(i,j,1) * (fzn_z(i,j,k) - fzn_z(i,j,nz)) / (fzn_z(i,j,1) - fzn_z(i,j,nz))
    end do
    fzn_zp(i,j,nz) = 0

    if(IsDebSideZ(i, j, -1)) then
      write(17,"(a)") ";i;j;k;fzn_z;fzn_w0;fzn_teta;fzn_rho;fzn_u;fzn_v;fzn_w"
      do k=1,nz
        if(IsDebSideZ(i, j, k)) then
          write(17,"('UFN:',3(';',i0),100(';',1p,g0.16))") &
              i, j, k, &
              fzn_z(i,j,k),fzn_w0(i,j,k),fzn_teta(i,j,k), fzn_rho(i,j,k), fzn_u(i,j,k), fzn_v(i,j,k), fzn_w(i,j,k)
        end if
      end do
    end if

  end do   ! i,j


end


