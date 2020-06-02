
!-------------------------------------------------------------------------------------------------------------------------------------------------

subroutine TransportInvX(itype, i, j, k, off, inv)
  use variables
  implicit none

  interface
    function GetCellInvX(itype, i, j, k, isN); integer(4)::itype,i,j,k; logical::isN; real(8)::GetCellInvX; end
    function GetSideInvX(itype, i, j, k, GG); integer(4)::itype,i,j,k; real(8)::GetSideInvX, GG; end
  end interface

  integer(4):: itype, i, j, k, off
  real(R8):: inv
  integer(4):: io, ic
  real(R8):: invo, invt, invc, invs, invn, sound, uc, lam, dx, q, Imin, Imax
  character(1):: iname(5)

  io = i - (2*off-1)                                            ! индекс заднего ребра
  ic = i - off                                                  ! индекс ячейки

  invo = GetSideInvX(itype, io, j, k, cs_G(ic, j, k))                             ! инвариант на задней грани
  invt = GetSideInvX(itype, i, j, k, cs_G(ic, j, k))                              ! инвариант на передней грани
  invc = GetCellInvX(itype, ic, j, k, .false.)                    ! инвариант в ячейке t[n]
  invs = GetCellInvX(itype, ic, j, k, .true.)                     ! инвариант в ячейке t[n+1/2]

  invn = 2.*invs - invo                                         ! новый инвариант предвариательно

  uc = cs_u0(ic,j) + cs_du(ic, j, k)                            ! скорость в ячейке
  sound = cs_sound                                              ! скорость звука в ячейке

  ! вычисляем характеристическую скорость lambda:
  select case(itype)
    case (T_INVR)
      lam = uc + sound                                 ! для инварианта R: lam = u + c
    case (T_INVQ)
      lam = uc - sound                                 ! для инварианта Q: lam = u - c
    case default
      lam = uc                                         ! для транспортных инвариантов: lam = u
  end select
  if(off==0) lam = -lam                                         ! если перенос влево
  dx = c_dx(ic)

  q = 2. * (invs - invc) + lam * (invt - invo) / dx * dt              ! аппроксимация правой части
  Imin = min(invo, invc, invt) + q                              ! нижня мажоранта
  Imax = max(invo, invc, invt) + q                              ! верхняя мажоранта

  inv = min(Imax, max(Imin, invn))                              ! лимитирование

  if(IsDebSideX(i, j, -1)) then
    iname = ['R', 'Q', 'U', 'D', 'W']
    write(17,"(/'TranspI-X:',3(';',i0))") i, j, k
    write(17,"(a)") ";i;j;k;dteta;drho;du;dv;dw;I;G"
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") &
      iname(itype), 'o', io, j, k, fx_dteta(io,j,k), fx_drho(io,j,k), fx_du(io,j,k), fx_dv(io,j,k), fx_dw(io,j,k), invo
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") &
      iname(itype), 't', i , j, k, fx_dteta(i ,j,k), fx_drho(i ,j,k), fx_du(i ,j,k), fx_dv(i ,j,k), fx_dw(i ,j,k), invt
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") &
      iname(itype), 'c', ic, j, k, c_dteta(ic,j,k), c_drho(ic,j,k), c_du(ic,j,k), c_dv(ic,j,k), c_dw(ic,j,k), invc, cs_G(ic, j, k)
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") &
      iname(itype), 's', ic, j, k, cs_dteta(ic,j,k), cs_drho(ic,j,k), cs_du(ic,j,k), cs_dv(ic,j,k), cs_dw(ic,j,k), invs, cs_G(ic, j, k)
    write(17,"(a,a,';;;;',1p,g0.16,a)") iname(itype), 'n', invn, &
          ";u;sound;lam;dx;q;Imin;Imax"
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") iname(itype), ' ', i , j, k, inv, &
          uc, sound, lam, dx, q, Imin, Imax
  end if

end subroutine TransportInvX

!-------------------------------------------------------------------------------------------------------------------------------------------------

subroutine TransportInvY(itype, i, j, k, off, inv)
  use variables
  implicit none

  interface
    function GetCellInvY(itype, i, j, k, isN); integer(4)::itype,i,j,k; logical::isN; real(8)::GetCellInvY; end
    function GetSideInvY(itype, i, j, k, GG); integer(4)::itype,i,j,k; real(8)::GetSideInvY, GG; end
  end interface

  integer(4):: itype, i, j, k, off
  real(R8):: inv
  integer(4):: jo, jc
  real(R8):: invo, invt, invc, invs, invn, sound, uc, lam, dy, q, Imin, Imax
  character(1):: iname(5)

  jo = j - (2*off-1)                                            ! индекс заднего ребра
  jc = j - off                                                  ! индекс ячейки

  invo = GetSideInvY(itype, i, jo, k, cs_G(i, jc, k))                             ! инвариант на задней грани
  invt = GetSideInvY(itype, i, j , k, cs_G(i, jc, k))                              ! инвариант на передней грани
  invc = GetCellInvY(itype, i, jc, k, .false.)                    ! инвариант в ячейке t[n]
  invs = GetCellInvY(itype, i, jc, k, .true.)                     ! инвариант в ячейке t[n+1/2]

  invn = 2.*invs - invo                                         ! новый инвариант предвариательно

  uc = cs_v0(i,jc) + cs_dv(i, jc, k)                            ! скорость в ячейке
  sound = cs_sound                                 ! скорость звука в ячейке

  ! вычисляем характеристическую скорость lambda:
  select case(itype)
    case (T_INVR)
      lam = uc + sound                                 ! для инварианта R: lam = u + c
    case (T_INVQ)
      lam = uc - sound                                 ! для инварианта Q: lam = u - c
    case default
      lam = uc                                         ! для транспортных инвариантов: lam = u
  end select
  if(off==0) lam = -lam                                         ! если перенос влево
  dy = c_dy(jc)

  q = 2. * (invs - invc) + lam * (invt - invo) / dy * dt              ! аппроксимация правой части
  Imin = min(invo, invc, invt) + q                              ! нижня мажоранта
  Imax = max(invo, invc, invt) + q                              ! верхняя мажоранта

  inv = min(Imax, max(Imin, invn))                              ! лимитирование

  if(IsDebSideY(i, j, -1)) then
    iname = ['R', 'Q', 'U', 'D', 'W']
    write(17,"(/'TranspI-Y:',3(';',i0))") i, j, k
    write(17,"(a)") ";i;j;k;dteta;drho;du;dv;dw;I;G"
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") &
      iname(itype), 'o', i, jo, k, fy_dteta(i,jo,k), fy_drho(i,jo,k), fy_du(i,jo,k), fy_dv(i,jo,k), fy_dw(i,jo,k), invo
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") &
      iname(itype), 't', i , j, k, fy_dteta(i ,j,k), fy_drho(i ,j,k), fy_du(i ,j,k), fy_dv(i ,j,k), fy_dw(i ,j,k), invt
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") &
      iname(itype), 'c', i, jc, k, c_dteta(i,jc,k), c_drho(i,jc,k), c_du(i,jc,k), c_dv(i,jc,k), c_dw(i,jc,k), invc, cs_G(i, jc, k)
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") &
      iname(itype), 's', i, jc, k, cs_dteta(i,jc,k), cs_drho(i,jc,k), cs_du(i,jc,k), cs_dv(i,jc,k), cs_dw(i,jc,k), invs, cs_G(i, jc, k)
    write(17,"(a,a,';;;;',1p,g0.16,a)") iname(itype), 'n', invn, &
          ";u;sound;lam;dx;q;Imin;Imax"
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") iname(itype), ' ', i , j, k, inv, &
          uc, sound, lam, dy, q, Imin, Imax
  end if

end subroutine TransportInvY

!-------------------------------------------------------------------------------------------------------------------------------------------------

subroutine TransportInvZ(itype, i, j, k, off, inv)
  use variables
  implicit none

  interface
    function GetCellInvZ(itype, i, j, k, isN); integer(4)::itype,i,j,k; logical::isN; real(8)::GetCellInvZ; end
    function GetSideInvZ(itype, i, j, k, GG); integer(4)::itype,i,j,k; real(8)::GetSideInvZ, GG; end
  end interface

  integer(4):: itype, i, j, k, off
  real(R8):: inv
  integer(4):: ko, kc
  real(R8):: invo, invt, invc, invs, invn, sound, wi, wc, czp, lam, dz, q, Imin, Imax
  character(1):: iname(5)

  ko = k - (2*off-1)                                            ! индекс заднего ребра
  kc = k - off                                                  ! индекс ячейки

  invo = GetSideInvZ(itype, i, j, ko, cs_G(i, j, kc))           ! инвариант на задней грани
  invt = GetSideInvZ(itype, i, j, k, cs_G(i, j, kc))            ! инвариант на передней грани
  invc = GetCellInvZ(itype, i, j, kc, .false.)                  ! инвариант в ячейке t[n]
  invs = GetCellInvZ(itype, i, j, kc, .true.)                   ! инвариант в ячейке t[n+1/2]

  invn = 2.*invs - invo                                         ! новый инвариант предвариательно

  wi = cs_GetW0(i,j,kc) + cs_dw(i,j,kc)                         ! эйлерова вертикальная скорость в слой-ячейке
  wc = wi - (fz_zp(i, j, k) + fz_zp(i, j, k)) / 2.              ! характеристическая скорость в ячейке
  sound = cs_sound                                              ! скорость звука в ячейке

  ! вычисляем характеристическую скорость lambda:
  select case(itype)
    case (T_INVR)
      lam = wc + sound                                 ! для инварианта R: lam = u + c
    case (T_INVQ)
      lam = wc - sound                                 ! для инварианта Q: lam = u - c
    case default
      lam = wc                                         ! для транспортных инвариантов: lam = u
  end select
  dz = fz_z(i,j,k) - fz_z(i,j,ko)
  !dz = fzs_z(i,j,k) - fzs_z(i,j,ko)

  q = 2. * (invs - invc) + lam * (invt - invo) / dz * dt              ! аппроксимация правой части
  Imin = min(invo, invc, invt) + q                              ! нижня мажоранта
  Imax = max(invo, invc, invt) + q                              ! верхняя мажоранта

  inv = min(Imax, max(Imin, invn))                              ! лимитирование

  if(IsDebSideZ(i, j, -1)) then
    iname = ['R', 'Q', 'U', 'D', 'W']
    write(17,"(/'TranspI-Z:',3(';',i0))") i, j, k
    write(17,"(a)") ";i;j;I;G"
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") iname(itype), 'o', i, j, ko, invo
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") iname(itype), 't', i, j, k , invt
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") iname(itype), 'c', i, j, kc, invc, cs_G(i, j, kc)
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") iname(itype), 's', i, j, kc, invs, cs_G(i, j, kc)
    write(17,"(a,a,';;;;',1p,g0.16,a)") iname(itype), 'n', invn, &
          ";wc;sound;lam;dz;q;Imin;Imax"
    write(17,"(a,a,3(';',i0),100(';',1p,g0.16))") iname(itype), ' ', i , j, k, inv, &
          wc, sound, lam, dz, q, Imin, Imax
  end if

end subroutine TransportInvZ

!------------------------------------------------------------------------------------------------------------

function GetCellInvX(itype, i, j, k, isN)
  use variables
  implicit none

  integer(4)::itype,i,j,k
  logical::isN
  real(8)::GetCellInvX, inv

  select case(itype)

    case (T_INVR)

      if(isN) then
        inv = cs_du(i,j,k) + cs_G(i,j,k) * cs_dteta(i,j,k)
      else
        inv = c_du(i,j,k) + cs_G(i,j,k) * c_dteta(i,j,k)
      end if

    case (T_INVQ)

      if(isN) then
        inv = cs_du(i,j,k) - cs_G(i,j,k) * cs_dteta(i,j,k)
      else
        inv = c_du(i,j,k) - cs_G(i,j,k) * c_dteta(i,j,k)
      end if

    case (T_INVU)
      if(isN) then; inv = cs_dv(i,j,k); else; inv = c_dv(i,j,k); endif

    case (T_INVW)
      if(isN) then; inv = cs_dw(i,j,k); else; inv = c_dw(i,j,k); endif

    case (T_INVD)
      if(isN) then; inv = cs_drho(i,j,k); else; inv = c_drho(i,j,k); endif

  end select

  GetCellInvX = inv

end

!------------------------------------------------------------------------------------------------------------

function GetSideInvX(itype, i, j, k, GG)
  use variables
  implicit none

  integer(4)::itype,i,j,k
  logical::isN
  real(8)::GetSideInvX, inv,GG

  select case(itype)

    case (T_INVR)
      inv = fx_du(i,j,k) + GG * fx_dteta(i,j,k)

    case (T_INVQ)
      inv = fx_du(i,j,k) - GG * fx_dteta(i,j,k)

    case (T_INVU)
      inv = fx_dv(i,j,k)

    case (T_INVW)
      inv = fx_dw(i,j,k)

    case (T_INVD)
      inv = fx_drho(i,j,k)

  end select

  GetSideInvX = inv
end

!------------------------------------------------------------------------------------------------------------

function GetCellInvY(itype, i, j, k, isN)
  use variables
  implicit none

  integer(4)::itype,i,j,k
  logical::isN
  real(8)::GetCellInvY, inv

  select case(itype)

    case (T_INVR)

      if(isN) then
        inv = cs_dv(i,j,k) + cs_G(i,j,k) * cs_dteta(i,j,k)
      else
        inv = c_dv(i,j,k) + cs_G(i,j,k) * c_dteta(i,j,k)
      end if

    case (T_INVQ)

      if(isN) then
        inv = cs_dv(i,j,k) - cs_G(i,j,k) * cs_dteta(i,j,k)
      else
        inv = c_dv(i,j,k) - cs_G(i,j,k) * c_dteta(i,j,k)
      end if

    case (T_INVU)
      if(isN) then; inv = cs_du(i,j,k); else; inv = c_dv(i,j,k); endif

    case (T_INVW)
      if(isN) then; inv = cs_dw(i,j,k); else; inv = c_dw(i,j,k); endif

    case (T_INVD)
      if(isN) then; inv = cs_drho(i,j,k); else; inv = c_drho(i,j,k); endif

  end select

  GetCellInvY = inv

end

!------------------------------------------------------------------------------------------------------------

function GetSideInvY(itype, i, j, k, GG)
  use variables
  implicit none

  integer(4)::itype,i,j,k
  logical::isN
  real(8)::GetSideInvY, inv,GG

  select case(itype)

    case (T_INVR)
      inv = fy_dv(i,j,k) + GG * fy_dteta(i,j,k)

    case (T_INVQ)
      inv = fy_dv(i,j,k) - GG * fy_dteta(i,j,k)

    case (T_INVU)
      inv = fy_du(i,j,k)

    case (T_INVW)
      inv = fy_dw(i,j,k)

    case (T_INVD)
      inv = fy_drho(i,j,k)

  end select

  GetSideInvY = inv
end

!------------------------------------------------------------------------------------------------------------

function GetCellInvZ(itype, i, j, k, isN)
  use variables
  implicit none

  integer(4)::itype,i,j,k
  logical::isN
  real(8)::GetCellInvZ, inv

  select case(itype)

    case (T_INVR)

      if(isN) then
        inv = cs_dw(i,j,k) + cs_G(i,j,k) * cs_dteta(i,j,k)
      else
        inv = c_dw(i,j,k) + cs_G(i,j,k) * c_dteta(i,j,k)
      end if

    case (T_INVQ)

      if(isN) then
        inv = cs_dw(i,j,k) - cs_G(i,j,k) * cs_dteta(i,j,k)
      else
        inv = c_dw(i,j,k) - cs_G(i,j,k) * c_dteta(i,j,k)
      end if

    case (T_INVU)
      if(isN) then; inv = cs_du(i,j,k); else; inv = c_dv(i,j,k); endif

    case (T_INVW)
      if(isN) then; inv = cs_dv(i,j,k); else; inv = c_dv(i,j,k); endif

    case (T_INVD)
      if(isN) then; inv = cs_drho(i,j,k); else; inv = c_drho(i,j,k); endif

  end select

  GetCellInvZ = inv

end

!------------------------------------------------------------------------------------------------------------

function GetSideInvZ(itype, i, j, k, GG)
  use variables
  implicit none

  integer(4)::itype,i,j,k
  logical::isN
  real(8)::GetSideInvZ, inv,GG

  select case(itype)

    case (T_INVR)
      inv = fz_dw(i,j,k) + GG * fz_dteta(i,j,k)

    case (T_INVQ)
      inv = fz_dw(i,j,k) - GG * fz_dteta(i,j,k)

    case (T_INVU)
      inv = fz_du(i,j,k)

    case (T_INVW)
      inv = fz_dv(i,j,k)

    case (T_INVD)
      inv = fz_drho(i,j,k)

  end select

  GetSideInvZ = inv
end

!-------------------------------------------------------------------------------------------------------------------------------------------------

! перенос dH
subroutine TransportHX(i, j, off, inv)
  use variables
  implicit none

  integer(4):: i, j, off
  real(R8):: inv
  integer(4):: io, ic
  real(R8):: invo, invt, invc, invs, invn, sound, uc, lam, dx, q, Imin, Imax

  io = i - (2*off-1)                                            ! индекс заднего ребра
  ic = i - off                                                  ! индекс ячейки

  ! инварианты inv = z[1] - H0:
  invo = fx_z(io,j,1) - fx_z0(io,j)                             ! на задней грани t[n]
  invt = fx_z(i,j,1) - fx_z0(i,j)                               ! на целевой грани t[n]
#if 0
!#warning ERROR !!!
  invc = fz_z(ic,j,1) - cn_z0(ic,j)                              ! в ячейке t[n+1]
#else
  invc = fz_z(ic,j,1) - c_z0(ic,j)                              ! в ячейке t[n]
#endif
  invs = fzs_z(ic,j,1) - cs_z0(ic,j)                            ! в ячейке t[n+1/2]

  invn = 2.*invs - invo                                         ! новый инвариант предвариательно

  uc = cs_u0(ic,j) + cs_du(ic, j, 1)                            ! скорость в ячейке
  lam = uc                                                      ! для транспортных инвариантов: lam = u
  if(off==0) lam = -lam                                         ! если перенос влево
  dx = c_dx(ic)

  q = 2. * (invs - invc) + lam * (invt - invo) / dx * dt        ! аппроксимация правой части
  Imin = min(invo, invc, invt) + q                              ! нижня мажоранта
  Imax = max(invo, invc, invt) + q                              ! верхняя мажоранта

  inv = min(Imax, max(Imin, invn))                              ! лимитирование

  if(IsDebSideX(i, j, -1)) then
    write(17,"(/'TranspH-X:',3(';',i0))") i, j
    write(17,"(a)") ";i;j;I"
    write(17,"(a,2(';',i0),';',100(';',1p,g0.16))") 'Ho', io, j, invo
    write(17,"(a,2(';',i0),';',100(';',1p,g0.16))") 'Ht', i , j, invt
    write(17,"(a,2(';',i0),';',100(';',1p,g0.16))") 'Hc', ic, j, invc
    write(17,"(a,2(';',i0),';',100(';',1p,g0.16))") 'Hs', ic, j, invs
    write(17,"(a,';;;;',1p,g0.16,a)") 'Hn', invn, &
          ";u;lam;dx;q;Imin;Imax"
    write(17,"(a,2(';',i0),';',100(';',1p,g0.16))") 'H', i , j, inv, &
          uc, lam, dx, q, Imin, Imax
  end if

end subroutine TransportHX

!-------------------------------------------------------------------------------------------------------------------------------------------------

! перенос dH
subroutine TransportHY(i, j, off, inv)
  use variables
  implicit none

  integer(4):: i, j, off
  real(R8):: inv
  integer(4):: jo, jc
  real(R8):: invo, invt, invc, invs, invn, sound, vc, lam, dy, q, Imin, Imax

  jo = j - (2*off-1)                                            ! индекс заднего ребра
  jc = j - off                                                  ! индекс ячейки

  ! инварианты inv = z[1] - H0:
  invo = fy_z(i,jo,1) - fy_z0(i,jo)                             ! на задней грани t[n]
  invt = fy_z(i,j,1) - fy_z0(i,j)                               ! на целевой грани t[n]
#if 0
!#  warning ERROR !!!
  invc = fz_z(i,jc,1) - cn_z0(i,jc)                              ! в ячейке t[n]
#else
  invc = fz_z(i,jc,1) - c_z0(i,jc)                              ! в ячейке t[n]
#endif
  invs = fzs_z(i,jc,1) - cs_z0(i,jc)                            ! в ячейке t[n+1/2]

  invn = 2.*invs - invo                                         ! новый инвариант предвариательно

  vc = cs_v0(i,jc) + cs_dv(i, jc,1)                            ! скорость в ячейке
  lam = vc                                                      ! для транспортных инвариантов: lam = u
  if(off==0) lam = -lam                                         ! если перенос влево
  dy = c_dy(jc)

  q = 2. * (invs - invc) + lam * (invt - invo) / dy * dt        ! аппроксимация правой части
  Imin = min(invo, invc, invt) + q                              ! нижня мажоранта
  Imax = max(invo, invc, invt) + q                              ! верхняя мажоранта

  inv = min(Imax, max(Imin, invn))                              ! лимитирование

  if(IsDebSideY(i, j, -1)) then
    write(17,"(/'TranspH-Y:',3(';',i0))") i, j
    write(17,"(a)") ";i;j;I"
    write(17,"(a,2(';',i0),';',100(';',1p,g0.16))") 'Ho', i, jo, invo
    write(17,"(a,2(';',i0),';',100(';',1p,g0.16))") 'Ht', i, j , invt
    write(17,"(a,2(';',i0),';',100(';',1p,g0.16))") 'Hc', i, jc, invc
    write(17,"(a,2(';',i0),';',100(';',1p,g0.16))") 'Hs', i, jc, invs
    write(17,"(a,';;;;',1p,g0.16,a)") 'Hn', invn, &
          ";v;lam;dy;q;Imin;Imax"
    write(17,"(a,2(';',i0),';',100(';',1p,g0.16))") 'H', i , j, inv, &
          vc, lam, dy, q, Imin, Imax
  end if

end subroutine TransportHY

