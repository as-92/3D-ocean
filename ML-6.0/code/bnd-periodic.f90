! периодические ГУ на X-гранях:
subroutine fx_BndPeriodic
  use variables
  implicit none

  real(R8) :: drhoi, dtetai, dui, dvi, dwi
  real(R8) :: gq, gr, M
  integer :: i, j, k, il, ir, is, ns
  real(R8) :: invD, invQ, invR, invU, invW
  real(R8) :: uf, ul, ur

  ns = nfx_sides(BC_PERIODIC)                                   ! количество граней
  do is=1, ns
    i = fx_sides(BC_PERIODIC).ptr(is).i                         ! индексы грани на сетке
    j = fx_sides(BC_PERIODIC).ptr(is).j
    call assertn(i==1 .or. i==nxx, "Phase2_X, периодические ГУ", i,j,-1)

    if(i/=1) cycle                                              ! только для левой границы

#if DLEV>0
    if(IsDebSideX(i, j, -1)) then
      write(17,"(/'Phase 2-X',2(';',i0),';;periodic')") i, j
    endif
#endif

    il = ncx                                                    ! ячейка слева (крайняя правая)
    ir = 1                                                      ! ячейка справа

    do k=1, ncz

      ! скорость на грани:
      !uf = fx_u0(i,j) + fx_du(i,j,k)                    ! скорость на грани
      ul = cs_u(il,j,k) !c_GetU(il,j,k,1)
      ur = cs_u(ir,j,k) !c_GetU(ir,j,k,1)
      uf=(ul + ur) / 2.
      M = uf / cs_sound

      if(abs(M)>0.9) then
        write(*,*) "ERROR in Phase2X: |M|>1. side (i,j,k):", i, j, k
        call avost
      endif

      ! инвариант R:
      call TransportInvX(T_INVR, nxx, j, k, DIRP, invR)        ! перенос через левую ячейку вправо
      Gr = cs_G(il,j,k)

      ! инвариант Q:
      call TransportInvX(T_INVQ, 1, j, k, DIRM, invQ)          ! перенос через правую ячейку влево
      Gq = -cs_G(ir,j,k)

      ! инварианты U и D:
      if(M>eps) then
        call TransportInvX(T_INVU, nxx, j, k, DIRP, invU)        ! перенос через левую ячейку вправо
        call TransportInvX(T_INVW, nxx, j, k, DIRP, invW)      ! перенос через левую ячейку вправо
        call TransportInvX(T_INVD, nxx, j, k, DIRP, invD)        ! перенос через левую ячейку вправо
      else if(M<-eps) then
        call TransportInvX(T_INVU, 1, j, k, DIRM, invU)          ! перенос через правую ячейку влево
        call TransportInvX(T_INVW, 1, j, k, DIRM, invW)      ! перенос через правую ячейку влево
        call TransportInvX(T_INVD, 1, j, k, DIRM, invD)          ! перенос через правую ячейку влево
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

      fxn_dteta(1,j,k) = dtetai
      fxn_du(1,j,k) = dui
      fxn_dv(1,j,k) = dvi                                 ! тангенцальная скорость
      fxn_dw(1,j,k) = dwi                                 ! вертикальная скорость
      fxn_drho(1,j,k) = drhoi                             ! плотность

      fxn_dteta(nxx,j,k) = fxn_dteta(1,j,k)
      fxn_du(nxx,j,k) = fxn_du(1,j,k)
      fxn_dv(nxx,j,k) = fxn_dv(1,j,k)
      fxn_dw(nxx,j,k) = fxn_dw(1,j,k)
      fxn_drho(nxx,j,k) = fxn_drho(1,j,k)

#if DLEV>0
      if(IsDebSideX(i, j, k)) then
        write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
        write(17,"('UF:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
              teta0+fx_dteta(i,j,k), rho0+fx_drho(i,j,k), fx_u0(i,j)+fx_du(i,j,k), fx_v0(i,j)+fx_dv(i,j,k), fx_w0(i,j,k)+fx_dw(i,j,k), &
              invR, invQ, invU, invW, invD, Gr, Gq
        write(17,"('UFN:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
              teta0+fxn_dteta(i,j,k), rho0+fxn_drho(i,j,k), fxn_u0(i,j)+fxn_du(i,j,k), fxn_v0(i,j)+fxn_dv(i,j,k), fxn_w0(i,j,k)+fxn_dw(i,j,k)
        i = nxx
        write(17,"('UFN:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
              teta0+fxn_dteta(i,j,k), rho0+fxn_drho(i,j,k), fxn_u0(i,j)+fxn_du(i,j,k), fxn_v0(i,j)+fxn_dv(i,j,k), fxn_w0(i,j,k)+fxn_dw(i,j,k)
      end if
#endif

    end do   ! do k

  end do  ! do is

end subroutine fx_BndPeriodic

!============================================================================================================

! периодические ГУ на Y-гранях:
subroutine fy_BndPeriodic
  use variables
  implicit none

  real(R8) :: drhoi, dtetai, dui, dvi, dwi
  real(R8) :: gq, gr, M
  integer :: i, j, k, jl, jr, is, ns
  real(R8) :: invD, invQ, invR, invU, invW
  real(R8) :: uf, ul, ur

  ns = nfy_sides(BC_PERIODIC)                                   ! количество граней
  do is=1, ns
    i = fy_sides(BC_PERIODIC).ptr(is).i                         ! индексы грани на сетке
    j = fy_sides(BC_PERIODIC).ptr(is).j
    call assertn(j==1 .or. j==nyy, "Phase2-Y, периодические ГУ", i,j,-1)

    if(j/=1) cycle                                              ! только для левой границы

    if(IsDebSideY(i, j, -1)) then
      write(17,"(/'Phase 2-Y',2(';',i0),';;periodic')") i, j
    endif

    jl = ncy                                                    ! ячейка слева (крайняя правая)
    jr = 1                                                      ! ячейка справа

    do k=1, ncz

      ! скорость на грани:
      !uf = fy_v0(i,j) + fy_dv(i,j,k)                    ! скорость на грани
      ul = cs_v(i,jl,k)
      ur = cs_v(i,jr,k)
      uf=(ul + ur) / 2.
      M = uf / cs_sound

      if(abs(M)>0.9) then
        write(*,*) "ERROR in Phase2Y-periodic: |M|>1. side (i,j,k):", i, j, k
        call avost
      endif

      ! инвариант R:
      call TransportInvY(T_INVR, i, nyy, k, DIRP, invR)        ! перенос через левую ячейку вправо
      Gr = cs_G(i,jl,k)

      ! инвариант Q:
      call TransportInvY(T_INVQ, i, 1, k, DIRM, invQ)          ! перенос через правую ячейку влево
      Gq = -cs_G(i,jr,k)

      ! инварианты U и D:
      if(M>eps) then
        call TransportInvY(T_INVU, i, nyy, k, DIRP, invU)        ! перенос через левую ячейку вправо
        call TransportInvY(T_INVW, i, nyy, k, DIRP, invW)      ! перенос через левую ячейку вправо
        call TransportInvY(T_INVD, i, nyy, k, DIRP, invD)        ! перенос через левую ячейку вправо
      else if(M<-eps) then
        call TransportInvY(T_INVU, i, 1, k, DIRM, invU)          ! перенос через правую ячейку влево
        call TransportInvY(T_INVW, i, 1, k, DIRM, invW)      ! перенос через правую ячейку влево
        call TransportInvY(T_INVD, i, 1, k, DIRM, invD)          ! перенос через правую ячейку влево
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

      !w0i = fyn_GetW0(i,j,k)

      fyn_dteta(i,1,k) = dtetai
      fyn_du(i,1,k) = dui
      fyn_dv(i,1,k) = dvi                            ! тангенцальная скорость
      fyn_dw(i,1,k) = dwi                            ! вертикальная скорость
      fyn_drho(i,1,k) = drhoi                ! плотность

      fyn_dteta(i,nyy,k) = fyn_dteta(i,1,k)
      fyn_du(i,nyy,k) = fyn_du(i,1,k)
      fyn_dv(i,nyy,k) = fyn_dv(i,1,k)
      fyn_dw(i,nyy,k) = fyn_dw(i,1,k)
      fyn_drho(i,nyy,k) = fyn_drho(i,1,k)

      if(IsDebSideY(i, j, k)) then
        write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
        write(17,"('side:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
              teta0+fy_dteta(i,j,k), rho0+fy_drho(i,j,k), fy_u0(i,j)+fy_du(i,j,k), fy_v0(i,j)+fy_dv(i,j,k), fy_w0(i,j,k)+fy_dw(i,j,k), &
              invR, invQ, invU, invW, invD, Gr, Gq
        write(17,"('UFN:;;;',100(';',1p,g0.16))") &
              teta0+fyn_dteta(i,j,k), rho0+fyn_drho(i,j,k), fyn_u0(i,j)+fyn_du(i,j,k), fyn_v0(i,j)+fyn_dv(i,j,k), fyn_w0(i,j,k)+fyn_dw(i,j,k)
      end if

    end do   ! do k

  end do  ! do is

end subroutine fy_BndPeriodic
