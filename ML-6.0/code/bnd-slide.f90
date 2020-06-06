! ГУ 'стенка с проскальзыванием' на X-гранях:
subroutine fx_BndSlide
  use variables
  implicit none

  real(R8) :: drhoi, dtetai, dui, dvi, dwi
  real(R8) :: gq, gr, M
  integer :: i, j, k, is, ic, ns
  real(R8) :: invD, invQ, invR, invU, invW
  real(R8) :: uf, ul, ur

  ns = nfx_sides(BC_SLIDE)                                      ! количество граней
  do is=1, ns
    i = fx_sides(BC_SLIDE).ptr(is).i                            ! индексы грани на сетке
    j = fx_sides(BC_SLIDE).ptr(is).j

    ! с какой стороны от стенки расчетная область?
    if(fx_IsOwnRight(i, j)) then                                ! живая ячейка справа  |<-Q

      ic = i                                                    ! индекс i ячейки
      do k=1, ncz

        ! вычисляем Мах на грани:
        uf = cs_u(ic, j, k)
        M = uf / cs_sound                                       ! Мах на грани

        if(abs(M)>0.9) then
          write(*,*) "ERROR in Phase2X-slide: |M|>1. side (i,j,k):", i, j, k
          call avost
        endif

        ! инвариант Q:
        call TransportInvX(T_INVQ, i, j, k, DIRM, invQ)         ! перенос через правую ячейку влево
        Gq = -cs_G(ic,j,k)

        ! транспортные инварианты:
        if(M<-eps) then
          call TransportInvX(T_INVU, i, j, k, DIRM, invU)       ! перенос через правую ячейку влево
          call TransportInvX(T_INVW, i, j, k, DIRM, invW)
          call TransportInvX(T_INVD, i, j, k, DIRM, invD)
        else
          invU = cs_dv(ic,j,k)                                  ! инвариант в центре ячейки
          invW = cs_dw(ic,j,k)
          invD = cs_drho(ic,j,k)
        end if

        ! восстанавливаем значения из инвариантов:
        drhoi = invD                                            ! дельта плотности
        dvi = invU                                              ! дельта тангенцальной скорости
        dwi = invW                                              ! дельта вертикальной скорости

        ! du = 0
        ! invQ = 0 + Gq * dteta
        !----------------------------------
        ! du = 0
        ! dteta = invQ / Gq

        dtetai = invQ / Gq                                      ! дельта псевдоплотности
        dui = 0.                                                ! дельта нормальной скорости

        fxn_dteta(i,j,k) = dtetai
        fxn_du(i,j,k) = dui
        fxn_dv(i,j,k) = dvi                                     ! тангенцальная скорость
        fxn_dw(i,j,k) = dwi                                     ! вертикальная скорость
        fxn_drho(i,j,k) = drhoi                                 ! плотность

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
        end if
#endif

      end do   ! k

    else                                                        ! живая ячейка слева, стенка справа: R->|

      ic = i - 1                                                ! индекс i ячейки
      do k = 1, ncz

        ! вычисляем число Маха на грани:
        uf = cs_u0(ic,j)                                          ! скорость на грани
        M = uf / cs_sound                                            ! Мах на грани

        if(abs(M)>0.9) then
          write(*,*) "ERROR in Phase2X-slide: |M|>1. side (i,j,k):", i, j, k
          call avost
        endif

        ! инвариант R:
        call TransportInvX(T_INVR, i, j, k, DIRP, invR)          ! перенос через левую ячейку вправо
        Gr = cs_G(ic,j,k)

        ! транспортные инварианты:
        if(M>eps) then                                            ! течение на стенку   >O|
          call TransportInvX(T_INVU, i, j, k, DIRP, invU)        ! .. перенос вправо
          call TransportInvX(T_INVW, i, j, k, DIRP, invW)
          call TransportInvX(T_INVD, i, j, k, DIRP, invD)
        else                                                      ! течение от стенки   <O|
          invU = cs_dv(ic,j,k)
          invW = cs_dw(ic,j,k)
          invD = cs_drho(ic,j,k)
        end if

        ! восстанавливаем значения из инвариантов:
        drhoi = invD                                        ! дельта плотности
        dvi = invU                                          ! дельта тангенцальной скорости
        dwi = invW                                          ! дельта вертикальной скорости

        ! du = 0
        ! invR = 0 + Gr * dteta
        !----------------------------------
        ! du = 0
        ! dteta = invR / Gr

        dtetai = invR / Gr
        dui = 0.

        fxn_dteta(i,j,k) = dtetai
        fxn_du(i,j,k) = dui
        fxn_dv(i,j,k) = dvi                                     ! тангенцальная скорость
        fxn_dw(i,j,k) = dwi                                     ! вертикальная скорость
        fxn_drho(i,j,k) = drhoi                                 ! плотность

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
        end if
#endif

      end do      ! k

    end if  ! if(c_type(i,j)>=0)

  end do  ! do is
end subroutine fx_BndSlide

!===============================================================================================================================

! ГУ 'стенка с проскальзыванием' на Y-гранях:
subroutine fy_BndSlide
  use variables
  implicit none

  real(R8) :: drhoi, dtetai, dui, dvi, dwi
  real(R8) :: gq, gr, M
  integer :: i, j, k, is, jc, ns
  real(R8) :: invD, invQ, invR, invU, invW
  real(R8) :: uf, ul, ur

  ns = nfx_sides(BC_SLIDE)                                      ! количество граней
  do is=1, ns
    i = fx_sides(BC_SLIDE).ptr(is).i                            ! индексы грани на сетке
    j = fx_sides(BC_SLIDE).ptr(is).j

    ! с какой стороны от стенки расчетная область?
    if(fy_IsOwnRight(i, j)) then                                ! живая ячейка справа  |<-Q

      jc = j                                                    ! индекс j ячейки
      do k=1, ncz

        ! вычисляем Мах на грани:
        uf = cs_u(i, jc, k)
        M = uf / cs_sound                                       ! Мах на грани

        if(abs(M)>0.9) then
          write(*,*) "ERROR in Phase2Y-slide: |M|>1. side (i,j,k):", i, j, k
          call avost
        endif

        ! инвариант Q:
        call TransportInvY(T_INVQ, i, j, k, DIRM, invQ)         ! перенос через правую ячейку влево
        Gq = -cs_G(i,jc,k)

        ! транспортные инварианты:
        if(M<-eps) then
          call TransportInvY(T_INVU, i, j, k, DIRM, invU)       ! перенос через правую ячейку влево
          call TransportInvY(T_INVW, i, j, k, DIRM, invW)
          call TransportInvY(T_INVD, i, j, k, DIRM, invD)
        else
          invU = cs_dv(i,jc,k)                                  ! инвариант в центре ячейки
          invW = cs_dw(i,jc,k)
          invD = cs_drho(i,jc,k)
        end if

        ! восстанавливаем значения из инвариантов:
        drhoi = invD                                            ! дельта плотности
        dvi = invU                                              ! дельта тангенцальной скорости
        dwi = invW                                              ! дельта вертикальной скорости

        ! du = 0
        ! invQ = 0 + Gq * dteta
        !----------------------------------
        ! du = 0
        ! dteta = invQ / Gq

        dtetai = invQ / Gq                                      ! дельта псевдоплотности
        dui = 0.                                                ! дельта нормальной скорости

        fyn_dteta(i,j,k) = dtetai
        fyn_du(i,j,k) = dui
        fyn_dv(i,j,k) = dvi                                     ! тангенцальная скорость
        fyn_dw(i,j,k) = dwi                                     ! вертикальная скорость
        fyn_drho(i,j,k) = drhoi                                 ! плотность

#if DLEV>0
        if(IsDebSideY(i, j, k)) then
          write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
          write(17,"('UF:',3(';',i0),100(';',1p,g0.16))") &
              i, j, k, &
                teta0+fy_dteta(i,j,k), rho0+fy_drho(i,j,k), fy_u0(i,j)+fy_du(i,j,k), fy_v0(i,j)+fy_dv(i,j,k), fy_w0(i,j,k)+fy_dw(i,j,k), &
                invR, invQ, invU, invW, invD, Gr, Gq
          write(17,"('UFN:',3(';',i0),100(';',1p,g0.16))") &
              i, j, k, &
                teta0+fyn_dteta(i,j,k), rho0+fyn_drho(i,j,k), fyn_u0(i,j)+fyn_du(i,j,k), fyn_v0(i,j)+fyn_dv(i,j,k), fyn_w0(i,j,k)+fyn_dw(i,j,k)
        end if
#endif

      end do   ! k

    else                                                        ! живая ячейка слева, стенка справа: R->|

      jc = j - 1                                                ! индекс j ячейки
      do k = 1, ncz

        ! вычисляем число Маха на грани:
        uf = cs_u0(i,jc)                                          ! скорость на грани
        M = uf / cs_sound                                         ! Мах на грани

        if(abs(M)>0.9) then
          write(*,*) "ERROR in Phase2Y-slide: |M|>1. side (i,j,k):", i, j, k
          call avost
        endif

        ! инвариант R:
        call TransportInvY(T_INVR, i, j, k, DIRP, invR)          ! перенос через левую ячейку вправо
        Gr = cs_G(i,jc,k)

        ! транспортные инварианты:
        if(M>eps) then                                            ! течение на стенку   >O|
          call TransportInvY(T_INVU, i, j, k, DIRP, invU)        ! .. перенос вправо
          call TransportInvY(T_INVW, i, j, k, DIRP, invW)
          call TransportInvY(T_INVD, i, j, k, DIRP, invD)
        else                                                      ! течение от стенки   <O|
          invU = cs_dv(i,jc,k)
          invW = cs_dw(i,jc,k)
          invD = cs_drho(i,jc,k)
        end if

        ! восстанавливаем значения из инвариантов:
        drhoi = invD                                        ! дельта плотности
        dvi = invU                                          ! дельта тангенцальной скорости
        dwi = invW                                          ! дельта вертикальной скорости

        ! du = 0
        ! invR = 0 + Gr * dteta
        !----------------------------------
        ! du = 0
        ! dteta = invR / Gr

        dtetai = invR / Gr
        dui = 0.

        fyn_dteta(i,j,k) = dtetai
        fyn_du(i,j,k) = dui
        fyn_dv(i,j,k) = dvi                                     ! тангенцальная скорость
        fyn_dw(i,j,k) = dwi                                     ! вертикальная скорость
        fyn_drho(i,j,k) = drhoi                                 ! плотность

#if DLEV>0
        if(IsDebSideY(i, j, k)) then
          write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
          write(17,"('UF:',3(';',i0),100(';',1p,g0.16))") &
              i, j, k, &
                teta0+fy_dteta(i,j,k), rho0+fy_drho(i,j,k), fy_u0(i,j)+fy_du(i,j,k), fy_v0(i,j)+fy_dv(i,j,k), fy_w0(i,j,k)+fy_dw(i,j,k), &
                invR, invQ, invU, invW, invD, Gr, Gq
          write(17,"('UFN:',3(';',i0),100(';',1p,g0.16))") &
              i, j, k, &
                teta0+fyn_dteta(i,j,k), rho0+fyn_drho(i,j,k), fyn_u0(i,j)+fyn_du(i,j,k), fyn_v0(i,j)+fyn_dv(i,j,k), fyn_w0(i,j,k)+fyn_dw(i,j,k)
        end if
#endif

      end do      ! k

    end if  ! if(c_type(i,j)>=0)

  end do  ! do is
end subroutine fy_BndSlide

