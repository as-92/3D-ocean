! ГУ "переменный вток" на X-гранях:
subroutine fx_BndInoutVar
  use variables
  implicit none

  integer:: i, j, k, is, ic, ns, ibc
  real(R8):: bcH, bcU, bcV, bcW, bcRho, bcTeta
  real(R8):: u1, u2, uf, M
  real(R8):: invR, invQ, invU, invW, invD, Gr, Gq
  real(R8):: dui, dvi, dwi, drhoi, dtetai
  logical:: ownRight

  ns = nfx_sides(BC_IN_T)                                       ! количество граней
  do is=1, ns
    i = fx_sides(BC_IN_T).ptr(is).i                             ! индексы грани на сетке
    j = fx_sides(BC_IN_T).ptr(is).j

    ibc = fx_bc(i,j)                                            ! индекс индивидуальных параметров ГУ

#if DLEV>0
    if(IsDebSideX(i, j, -1)) then
      write(17,"(/'Phase 2-X',2(';',i0),';;BC_IN_T',i0)") i, j, ibc
    endif
#endif

    ! заграничные данные:
    bcH = bcData(ibc).h
    bcU = bcData(ibc).u
    bcV = bcData(ibc).v
    bcW = bcData(ibc).w
    bcRho = bcData(ibc).rho
    bcTeta = bcData(ibc).teta

    ownRight = fx_IsOwnRight(i,j)                                  ! лежит ли ячейка справа

    ic = i; if(ownRight) ic = i + 1

    ! вычисляем инварианты из-за границы "как есть", без изменения глубины

    do k=1, ncz

      ! скорость на грани:
      u1 = cs_u(ic,j,k)                                         ! скорость в ячейке
      u2 = bcU                                                  ! скорость извне
      uf = (u1 + u2) / 2.
      M = uf / cs_sound                                         ! число Маха для переноса инвариантов

      if(abs(M)>0.9) then
        write(*,*) "ERROR in Phase2Y: |M|>1. side (i,j,k):", i, j, k
        call avost
      endif

      ! инвариант R (слева направо):
      if(ownRight) then                                         ! ячейка справа, граница слева
        Gr = cs_sound / bcTeta
        invR = 0.                                               ! .. параметры мелкой воды и многослойные совпадают
      else                                                      ! слева есть ячейка
        call TransportInvX(T_INVR, i, j, k, DIRP, invR)         ! .. перенос через левую ячейку вправо
        Gr = cs_G(i-1,j,k)
      end if

      ! инвариант Q (справа налево):
      if(ownRight) then                                         ! справа есть ячейка
        call TransportInvX(T_INVQ, i, j, k, DIRM, invQ)         ! .. перенос через правую ячейку влево
        Gq = -cs_G(i,j,k)
      else
        Gq = -cs_sound / bcTeta
        invQ = 0.
      end if

      ! транспортные инварианты:
      if((ownRight .and. M<-eps) .or. &                         ! течение из ячейки на границу   |<O
         (.not.ownRight .and. M>eps)) then                      ! течение из ячейки на границу   O>|
        invU = cs_du(ic,j,k)
        invW = cs_dw(ic,j,k)
        invD = cs_drho(ic,j,k)
      else                                                      ! течение от границы в ячейку
        invU = 0.
        invW = 0.
        invD = 0.
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

      fxn_dteta(i,j,k) = dtetai                           ! псевдоплотность
      fxn_du(i,j,k) = dui                                 ! нормальная скорость
      fxn_dv(i,j,k) = dvi                                 ! тангенцальная скорость
      fxn_dw(i,j,k) = dwi                                 ! вертикальная скорость
      fxn_drho(i,j,k) = drhoi                             ! плотность

#if DLEV>0
      if(IsDebSideX(i, j, k)) then
        write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
        write(17,"('BC:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
            bcTeta, bcRho, bcU, bcV, bcW
        write(17,"('side:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
              teta0+fx_dteta(i,j,k), rho0+fx_drho(i,j,k), fx_u0(i,j)+fx_du(i,j,k), fx_v0(i,j)+fx_dv(i,j,k), fx_w0(i,j,k)+fx_dw(i,j,k), &
              invR, invQ, invU, invW, invD, Gr, Gq
        write(17,"('UFN:;;;',100(';',1p,g0.16))") &
              teta0+fxn_dteta(i,j,k), rho0+fxn_drho(i,j,k), fxn_u0(i,j)+fxn_du(i,j,k), fxn_v0(i,j)+fxn_dv(i,j,k), fxn_w0(i,j,k)+fxn_dw(i,j,k)
      end if
#endif

    end do      ! k

  end do  ! do is

end subroutine fx_BndInoutVar

!===================================================================================================================

! ГУ "переменный вток" на Y-гранях:
subroutine fy_BndInoutVar
  use variables
  implicit none

  integer:: i, j, k, is, jc, ns, ibc
  real(R8):: bcH, bcU, bcV, bcW, bcRho, bcTeta
  real(R8):: u1, u2, uf, M
  real(R8):: invR, invQ, invU, invW, invD, Gr, Gq
  real(R8):: dui, dvi, dwi, drhoi, dtetai
  logical:: ownRight

  ns = nfy_sides(BC_IN_T)                                         ! количество граней
  do is=1, ns
    i = fy_sides(BC_IN_T).ptr(is).i                               ! индексы грани на сетке
    j = fy_sides(BC_IN_T).ptr(is).j

    ibc = fy_bc(i,j)                                            ! индекс индивидуальных параметров ГУ

#if DLEV>0
    if(IsDebSideY(i, j, -1)) then
      write(17,"(/'Phase 2-Y',2(';',i0),';;BC_IN_T',i0)") i, j, ibc
    endif
#endif

    ! заграничные данные:
    bcH = bcData(ibc).h
    bcU = bcData(ibc).u
    bcV = bcData(ibc).v
    bcW = bcData(ibc).w
    bcRho = bcData(ibc).rho
    bcTeta = bcData(ibc).teta

    ownRight = fy_IsOwnRight(i,j)                                  ! лежит ли ячейка справа

    jc = j; if(ownRight) jc = j + 1

    ! вычисляем инварианты из-за границы "как есть", без изменения глубины

    do k=1, ncz

      ! скорость на грани:
      u1 = cs_u(i,jc,k)                                         ! скорость в ячейке
      u2 = bcU                                                  ! скорость извне
      uf = (u1 + u2) / 2.
      M = uf / cs_sound                                         ! число Маха для переноса инвариантов

      if(abs(M)>0.9) then
        write(*,*) "ERROR in Phase2Y: |M|>1. side (i,j,k):", i, j, k
        call avost
      endif

      ! инвариант R (слева направо):
      if(ownRight) then                                         ! ячейка справа, граница слева
        Gr = cs_sound / bcTeta
        invR = 0.                                               ! .. параметры мелкой воды и многослойные совпадают
      else                                                      ! слева есть ячейка
        call TransportInvY(T_INVR, i, j, k, DIRP, invR)         ! .. перенос через левую ячейку вправо
        Gr = cs_G(i,jc,k)
      end if

      ! инвариант Q (справа налево):
      if(ownRight) then                                         ! справа есть ячейка
        call TransportInvY(T_INVQ, i, j, k, DIRM, invQ)         ! .. перенос через правую ячейку влево
        Gq = -cs_G(i,jc,k)
      else
        Gq = -cs_sound / bcTeta
        invQ = 0.
      end if

      ! транспортные инварианты:
      if((ownRight .and. M<-eps) .or. &                         ! течение из ячейки на границу   |<O
         (.not.ownRight .and. M>eps)) then                      ! течение из ячейки на границу   O>|
        invU = cs_du(i,jc,k)
        invW = cs_dw(i,jc,k)
        invD = cs_drho(i,jc,k)
      else                                                      ! течение от границы в ячейку
        invU = 0.
        invW = 0.
        invD = 0.
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

      fyn_dteta(i,j,k) = dtetai                           ! псевдоплотность
      fyn_du(i,j,k) = dui                                 ! нормальная скорость
      fyn_dv(i,j,k) = dvi                                 ! тангенцальная скорость
      fyn_dw(i,j,k) = dwi                                 ! вертикальная скорость
      fyn_drho(i,j,k) = drhoi                             ! плотность

#if DLEV>0
      if(IsDebSideY(i, j, k)) then
        write(17,"(a)") ";i;j;k;teta;rho;u;v;w;InvR;InvQ;InvU;InvW;InvD;Gr;Gq"
        write(17,"('BC:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
            bcTeta, bcRho, bcU, bcV, bcW
        write(17,"('side:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, &
              teta0+fy_dteta(i,j,k), rho0+fy_drho(i,j,k), fy_u0(i,j)+fy_du(i,j,k), fy_v0(i,j)+fy_dv(i,j,k), fy_w0(i,j,k)+fy_dw(i,j,k), &
              invR, invQ, invU, invW, invD, Gr, Gq
        write(17,"('UFN:;;;',100(';',1p,g0.16))") &
              teta0+fyn_dteta(i,j,k), rho0+fyn_drho(i,j,k), fyn_u0(i,j)+fyn_du(i,j,k), fyn_v0(i,j)+fyn_dv(i,j,k), fyn_w0(i,j,k)+fyn_dw(i,j,k)
      end if
#endif

    end do      ! k

  end do  ! do is

end subroutine fy_BndInoutVar

