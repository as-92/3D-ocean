subroutine Phase1
  use variables
  implicit none

  integer(4):: i, j, k, ic
  real(R8):: dx, dy, dwz(nz)
  real(R8):: q_xr, q_xl, q_xt, q_xb, q_yr, q_yl, q_yt, q_yb, q_t, q_b
  real(R8):: dz_c, dz_xr, dz_xl, dz_xb, dz_xt, dz_yr, dz_yl, dz_yb, dz_yt
  real(R8):: u_c, u_xr, u_xl, u_xt, u_xb, u_yr, u_yl, u_yt, u_yb, u_t, u_b
  real(R8):: v_c, v_xr, v_xl, v_xt, v_xb, v_yr, v_yl, v_yt, v_yb, v_t, v_b
  real(R8):: w_c, w_xr, w_xl, w_xt, w_xb, w_yr, w_yl, w_yt, w_yb, w_t, w_b
  real(R8):: teta_c, teta_xr, teta_xl, teta_xt, teta_xb, teta_yr, teta_yl, teta_yt, teta_yb, teta_t, teta_b
  real(R8):: rho_c, rho_xr, rho_xl, rho_xt, rho_xb, rho_yr, rho_yl, rho_yt, rho_yb, rho_t, rho_b
  real(R8):: dp_c, dp_xr, dp_xl, dp_xt, dp_xb, dp_yr, dp_yl, dp_yt, dp_yb, dp_t, dp_b
  real(R8):: hSw_xr, hSw_xl, hSw_yr, hSw_yl, hSw_t, hSw_b
  real(R8):: teta_h_new, rho_h_new, u_h_new, v_h_new, w_h_new

  do ic=1,ncells                                                ! цикл по плоcким ячейкам
    i = cells.ptr(ic).i
    j = cells.ptr(ic).j

    dx = c_dx(i)
    dy = c_dy(j)

    do k=1,nz
      dwz(k) = (fz_GetW0(i,j,k) + fz_dw(i,j,k)) - fz_zp(i,j,k)
      fzs_z(i,j,k) = fz_z(i,j,k) + dt/2 * fz_zp(i,j,k)
    end do

    do k=1,ncz
      cs_w0(i,j,k) = cs_GetW0(i,j,k)
    end do

    if(IsDebCell(i, j, -1)) then
      write(17,"(/'Phase 1',2(';',i0),';pos:',100(';',f0.5))") i, j, (x(i)+x(i+1))/2., (y(j)+y(j+1))/2.
      write(17,"(a)") ";i;j;k;teta;drho;u;v;w;zT;zB;dp;H(sw);dz"
    end if

    do k=1,ncz                                                  ! цикл по cлоям

      !-- значения в cлой-ячейке и её окреcтноcтям ------------------------------------------------------

      ! в cлой-ячейке t[n]:
      dz_c = fz_z(i,j,k) - fz_z(i,j,k+1)
      teta_c = c_teta(i,j,k)
      rho_c = c_rho(i,j,k)
      dp_c = sound0_2 * c_dteta(i,j,k)
      u_c = c_u(i,j,k)
      v_c = c_v(i,j,k)
      w_c = c_w(i,j,k)

      ! на грани X-left (t=t[n]):
      dz_xl = fx_z(i,j,k) - fx_z(i,j,k+1)
      teta_xl = fx_teta(i,j,k)
      rho_xl = fx_rho(i,j,k)
      dp_xl = sound0_2 * fx_dteta(i,j,k)
      u_xl = fx_u(i,j,k)
      v_xl = fx_v(i,j,k)
      w_xl = fx_w(i,j,k)
      hSw_xl = fx_z0(i,j)                          ! уровень поверхности МВ на грани

      ! на грани X-right (t=t[n]):
      dz_xr = fx_z(i+1,j,k) - fx_z(i+1,j,k+1)
      teta_xr = fx_teta(i+1,j,k)
      rho_xr = fx_rho(i+1,j,k)
      dp_xr = sound0_2 * fx_dteta(i+1,j,k)
      u_xr = fx_u(i+1,j,k)
      v_xr = fx_v(i+1,j,k)
      w_xr = fx_w(i+1,j,k)
      hSw_xr = fx_z0(i+1,j)                          ! уровень поверхности МВ на грани

      ! на грани Y-left (t=t[n]):
      dz_yl = fy_z(i,j,k) - fy_z(i,j,k+1)
      teta_yl = fy_teta(i,j,k)
      rho_yl = fy_rho(i,j,k)
      dp_yl = sound0_2 * fy_dteta(i,j,k)
      u_yl = fy_u(i,j,k)
      v_yl = fy_v(i,j,k)
      w_yl = fy_w(i,j,k)
      hSw_yl = fy_z0(i,j)                          ! уровень поверхности МВ на грани

      ! на грани Y-right (t=t[n]):
      dz_yr = fy_z(i,j+1,k) - fy_z(i,j+1,k+1)
      teta_yr = fy_teta(i,j+1,k)
      rho_yr = fy_rho(i,j+1,k)
      dp_yr = sound0_2 * fy_dteta(i,j+1,k)
      u_yr = fy_u(i,j+1,k)
      v_yr = fy_v(i,j+1,k)
      w_yr = fy_w(i,j+1,k)
      hSw_yr = fy_z0(i,j+1)                          ! уровень поверхности МВ на грани

      ! на грани Z-top (t=t[n]):
      dz_xt = fx_z(i+1,j,k) - fx_z(i,j,k)
      dz_yt = fy_z(i,j+1,k) - fy_z(i,j,k)
      rho_t = fz_rho(i,j,k)                  ; rho_xt  = rho_t ; rho_yt  = rho_t
      dp_t = sound0_2 * fz_dteta(i,j,k)      ; dp_xt   = dp_t  ; dp_yt   = dp_t
      u_t = fz_u(i,j,k)                      ; u_xt    = u_t   ; u_yt    = u_t
      v_t = fz_v(i,j,k)                      ; v_xt    = v_t   ; v_yt    = v_t
      w_t = fz_w(i,j,k)                      ; w_xt    = w_t   ; w_yt    = w_t
      teta_t = fz_teta(i,j,k)                ; teta_xt = teta_t; teta_yt = teta_t
      hSw_t = c_z0(i,j)                          ! уровень поверхности МВ в ячейке

      ! на грани Z-bottom (t=t[n]):
      dz_xb = fx_z(i+1,j,k+1) - fx_z(i,j,k+1)
      dz_yb = fy_z(i,j+1,k+1) - fy_z(i,j,k+1)
      rho_b = fz_rho(i,j,k+1)                  ; rho_xb  = rho_b ; rho_yb  = rho_b
      dp_b = sound0_2 * fz_dteta(i,j,k+1)      ; dp_xb   = dp_b  ; dp_yb   = dp_b
      u_b = fz_u(i,j,k+1)                      ; u_xb    = u_b   ; u_yb    = u_b
      v_b = fz_v(i,j,k+1)                      ; v_xb    = v_b   ; v_yb    = v_b
      w_b = fz_w(i,j,k+1)                      ; w_xb    = w_b   ; w_yb    = w_b
      teta_b = fz_teta(i,j,k+1)                ; teta_xb = teta_b; teta_yb = teta_b
      hSw_b = hSw_t

      !-- баланcные уравнения -------------------------------------------------------

      ! баланc teta:

      q_xr = teta_xr * u_xr                                     ! грань X-право
      q_xl = teta_xl * u_xl                                     ! грань X-лево
      q_xt = teta_xt * u_xt                                     ! грань Z-верх, производная по X
      q_xb = teta_xb * u_xb                                     ! грань Z-низ, производная по X

      q_yr = teta_yr * v_yr                                     ! грань Y-право
      q_yl = teta_yl * v_yl                                     ! грань Y-лево
      q_yt = teta_yt * v_yt                                     ! грань Z-верх, производная по Y
      q_yb = teta_yb * v_yb                                     ! грань Z-низ, производная по Y

      q_t = teta_t                                              ! грань Z-верх
      q_b = teta_b                                              ! грань Z-низ

      teta_h_new = teta_c * dz_c - ( &
                   ( q_xr * dz_xr - q_xl * dz_xl) / dx + &
                   ( q_xb * dz_xb - q_xt * dz_xt) / dx + &
                   ( q_yr * dz_yr - q_yl * dz_yl) / dy + &
                   ( q_yb * dz_yb - q_yt * dz_yt) / dy + &
                   ( q_t * dwz(k) - q_b * dwz(k+1) ) &
                                    ) * dt / 2.

      ! баланc rho:

      q_xr = teta_xr * u_xr * rho_xr                           ! грань X-право
      q_xl = teta_xl * u_xl * rho_xl                           ! грань X-лево
      q_xt = teta_xt * u_xt * rho_xt                           ! грань Z-верх, производная по X
      q_xb = teta_xb * u_xb * rho_xb                           ! грань Z-низ, производная по X

      q_yr = teta_yr * v_yr * rho_yr                           ! грань Y-право
      q_yl = teta_yl * v_yl * rho_yl                           ! грань Y-лево
      q_yt = teta_yt * v_yt * rho_yt                           ! грань Z-верх, производная по Y
      q_yb = teta_yb * v_yb * rho_yb                           ! грань Z-низ, производная по Y

      q_t = teta_t * rho_t                                     ! грань Z-верх
      q_b = teta_b * rho_b                                     ! грань Z-низ

      rho_h_new  = teta_c * dz_c * rho_c - ( &
                   ( q_xr * dz_xr - q_xl * dz_xl) / dx + &
                   ( q_xb * dz_xb - q_xt * dz_xt) / dx + &
                   ( q_yr * dz_yr - q_yl * dz_yl) / dy + &
                   ( q_yb * dz_yb - q_yt * dz_yt) / dy + &
                   ( q_t * dwz(k) - q_b * dwz(k+1) ) &
                                    ) * dt / 2.

      ! баланc u:

      q_xr = teta_xr * u_xr * u_xr + &                          ! грань X-право
             dp_xr / rho0 + &
             g * hSw_xr
      q_xl = teta_xl * u_xl * u_xl + &                          ! грань X-лево
             dp_xl / rho0 + &
             g * hSw_xl
      q_xt = teta_xt * u_xt * u_xt + &                          ! грань Z-верх, производная по X
             dp_xt / rho0 + &
             g * hSw_t
      q_xb = teta_xb * u_xb * u_xb + &                          ! грань Z-низ, производная по X
             dp_xb / rho0 + &
             g * hSw_b

      q_yr = teta_yr * v_yr * u_yr                              ! грань Y-право
      q_yl = teta_yl * v_yl * u_yl                              ! грань Y-лево
      q_yt = teta_yt * v_yt * u_yt                              ! грань Z-верх, производная по Y
      q_yb = teta_yb * v_yb * u_yb                              ! грань Z-низ, производная по Y

      q_t = teta_t * u_t                                        ! грань Z-верх
      q_b = teta_b * u_b                                        ! грань Z-низ

      u_h_new  = teta_c * dz_c * u_c - ( &
                   ( q_xr * dz_xr - q_xl * dz_xl) / dx + &
                   ( q_xb * dz_xb - q_xt * dz_xt) / dx + &
                   ( q_yr * dz_yr - q_yl * dz_yl) / dy + &
                   ( q_yb * dz_yb - q_yt * dz_yt) / dy + &
                   ( q_t * dwz(k) - q_b * dwz(k+1) ) &
                                    ) * dt / 2.

      ! баланc v:

      q_xr = teta_xr * u_xr * v_xr                              ! грань X-право
      q_xl = teta_xl * u_xl * v_xl                              ! грань X-лево
      q_xt = teta_xt * u_xt * v_xt                              ! грань Z-верх, производная по X
      q_xb = teta_xb * u_xb * v_xb                              ! грань Z-низ, производная по X

      q_yr = teta_yr * v_yr * v_yr + &                          ! грань Y-право
             dp_yr / rho0 + &
             g * hSw_yr
      q_yl = teta_yl * v_yl * v_yl + &                          ! грань Y-лево
             dp_yl / rho0 + &
             g * hSw_yl
      q_yt = teta_yt * v_yt * v_yt + &                          ! грань Z-верх, производная по Y
             dp_yt / rho0 + &
             g * hSw_t
      q_yb = teta_yb * v_yb * v_yb + &                          ! грань Z-низ, производная по Y
             dp_yb / rho0 + &
             g * hSw_b

      q_t = teta_t * v_t                                        ! грань Z-верх
      q_b = teta_b * v_b                                        ! грань Z-низ

      v_h_new  = teta_c * dz_c * v_c - ( &
                   ( q_xr * dz_xr - q_xl * dz_xl) / dx + &
                   ( q_xb * dz_xb - q_xt * dz_xt) / dx + &
                   ( q_yr * dz_yr - q_yl * dz_yl) / dy + &
                   ( q_yb * dz_yb - q_yt * dz_yt) / dy + &
                   ( q_t * dwz(k) - q_b * dwz(k+1) ) &
                                    ) * dt / 2.

      ! баланc w:

      q_xr = teta_xr * u_xr * w_xr                              ! грань X-право
      q_xl = teta_xl * u_xl * w_xl                              ! грань X-лево
      q_xt = teta_xt * u_xt * w_xt                              ! грань Z-верх, производная по X
      q_xb = teta_xb * u_xb * w_xb                              ! грань Z-низ, производная по X

      q_yr = teta_yr * v_yr * w_yr                              ! грань Y-право
      q_yl = teta_yl * v_yl * w_yl                              ! грань Y-лево
      q_yt = teta_yt * v_yt * w_yt                              ! грань Z-верх, производная по Y
      q_yb = teta_yb * v_yb * w_yb                              ! грань Z-низ, производная по Y

      q_t = teta_t * w_t                                        ! грань Z-верх
      q_b = teta_b * w_b                                        ! грань Z-низ

      w_h_new  = teta_c * dz_c * w_c - ( &
                   ( q_xr * dz_xr - q_xl * dz_xl) / dx + &
                   ( q_xb * dz_xb - q_xt * dz_xt) / dx + &
                   ( q_yr * dz_yr - q_yl * dz_yl) / dy + &
                   ( q_yb * dz_yb - q_yt * dz_yt) / dy + &
                   ( (q_t * dwz(k) + dp_t / rho0) - (q_b * dwz(k+1) + dp_b * rho0) ) + &
                   g * (rho_c - rho0) / rho0 * dz_c &
                                    ) * dt / 2.

      !-- новые консервативные переменные ----------------------------------------------------------------------------

      cs_teta(i,j,k) = teta_h_new / (fzs_z(i,j,k) - fzs_z(i,j,k+1))
      cs_rho(i,j,k) = rho_h_new / teta_h_new
      cs_u(i,j,k) = u_h_new / teta_h_new
      cs_v(i,j,k) = v_h_new / teta_h_new
      cs_w(i,j,k) = w_h_new / teta_h_new

      cs_dteta(i,j,k) = cs_teta(i,j,k) - teta0
      cs_drho(i,j,k) = cs_rho(i,j,k) - rho0
      cs_du(i,j,k) = cs_u(i,j,k) - cs_u0(i,j)
      cs_dv(i,j,k) = cs_v(i,j,k) - cs_v0(i,j)
      cs_dw(i,j,k) = cs_w(i,j,k) - cs_w0(i,j,k)

      cs_G(i,j,k) = cs_GetG(i,j,k)

      if(IsDebCell(i, j, k)) then
        write(17,"('side-XL:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, teta_xl, rho_xl, u_xl, v_xl, w_xl, fx_z(i,j,k), fx_z(i,j,k+1), dp_xl, hSw_xl, dz_xl; flush(17)
        write(17,"('side-XR:',3(';',i0),100(';',1p,g0.16))") &
            i+1, j, k, teta_xr, rho_xr, u_xr, v_xr, w_xr, fx_z(i+1,j,k), fx_z(i+1,j,k+1), dp_xr, hSw_xr, dz_xr; flush(17)
        write(17,"('side-YL:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, teta_yl, rho_yl, u_yl, v_yl, w_yl, fy_z(i,j,k), fy_z(i,j,k+1), dp_yl, hSw_yl, dz_yl; flush(17)
        write(17,"('side-YR:',3(';',i0),100(';',1p,g0.16))") &
            i, j+1, k, teta_yr, rho_yr, u_yr, v_yr, w_yr, fy_z(i,j+1,k), fy_z(i,j+1,k+1), dp_yr, hSw_yr, dz_yr; flush(17)
        write(17,"('side-B:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, teta_b, rho_b, u_b, v_b, w_b, fz_z(i,j,k+1), 0., dp_b, hSw_b, dz_xb, dz_yb; flush(17)
        write(17,"('side-T:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k+1, teta_t, rho_t, u_t, v_t, w_t, fz_z(i,j,k), 0., dp_t, hSw_t, dz_xt, dz_yt; flush(17)
        write(17,"('UC:' ,3(';',i0),100(';',1p,g0.16))") i, j, k,  c_teta(i,j,k),  c_rho(i,j,k),  c_u(i,j,k),  c_v(i,j,k),  c_w(i,j,k)
        write(17,"('UCS:',3(';',i0),100(';',1p,g0.16))") i, j, k, cs_teta(i,j,k), cs_rho(i,j,k), cs_u(i,j,k), cs_v(i,j,k), cs_w(i,j,k)
        flush(17)
      end if

    end do
  end do

end

!-------------------------------------------------------------------------------------------------------------------------------------------------

subroutine Phase3
  use variables
  implicit none


  integer(4):: i, j, k, ic
  real(R8):: dx, dy, dwz(nz), w0
  real(R8):: q_xr, q_xl, q_xt, q_xb, q_yr, q_yl, q_yt, q_yb, q_t, q_b
  real(R8):: dz_c, dz_xr, dz_xl, dz_xb, dz_xt, dz_yr, dz_yl, dz_yb, dz_yt
  real(R8):: u_c, u_xr, u_xl, u_xt, u_xb, u_yr, u_yl, u_yt, u_yb, u_t, u_b
  real(R8):: v_c, v_xr, v_xl, v_xt, v_xb, v_yr, v_yl, v_yt, v_yb, v_t, v_b
  real(R8):: w_c, w_xr, w_xl, w_xt, w_xb, w_yr, w_yl, w_yt, w_yb, w_t, w_b
  real(R8):: teta_c, teta_xr, teta_xl, teta_xt, teta_xb, teta_yr, teta_yl, teta_yt, teta_yb, teta_t, teta_b
  real(R8):: rho_c, rho_xr, rho_xl, rho_xt, rho_xb, rho_yr, rho_yl, rho_yt, rho_yb, rho_t, rho_b
  real(R8):: dp_c, dp_xr, dp_xl, dp_xt, dp_xb, dp_yr, dp_yl, dp_yt, dp_yb, dp_t, dp_b
  real(R8):: hSw_xr, hSw_xl, hSw_yr, hSw_yl, hSw_t, hSw_b
  real(R8):: teta_h_new, rho_h_new, u_h_new, v_h_new, w_h_new

  do ic=1,ncells                                                ! цикл по плоcким ячейкам
    i = cells.ptr(ic).i
    j = cells.ptr(ic).j

    dx = c_dx(i)
    dy = c_dy(j)

    do k=1,nz
      dwz(k) = fzn_w(i,j,k) - fzn_zp(i,j,k)
      debdum = 0
    end do

    if(IsDebCell(i, j, -1)) then
      write(17,"(/'Phase 3',2(';',i0),';pos:',100(';',f0.5))") i, j, (x(i)+x(i+1))/2., (y(j)+y(j+1))/2.
      write(17,"(a)") ";i;j;k;teta;drho;u;v;w;dp;H(sw);dz"
    end if

    do k=1,ncz                                                  ! цикл по cлоям

      !-- значения в cлой-ячейке и её окреcтноcтям ------------------------------------------------------

      ! в cлой-ячейке t[n+1/2]:
      dz_c = fzs_z(i,j,k) - fzs_z(i,j,k+1)
      teta_c = cs_teta(i,j,k)
      rho_c = cs_rho(i,j,k)
      dp_c = sound0_2 * cs_dteta(i,j,k)
      u_c = cs_u(i,j,k)
      v_c = cs_v(i,j,k)
      w_c = cs_w(i,j,k)

      ! на грани X-left (t=t[n+1]):
      dz_xl = fxn_z(i,j,k) - fxn_z(i,j,k+1)
      teta_xl = fxn_teta(i,j,k)
      rho_xl = fxn_rho(i,j,k)
      dp_xl = sound0_2 * fxn_dteta(i,j,k)
      u_xl = fxn_u(i,j,k)
      v_xl = fxn_v(i,j,k)
      w_xl = fxn_w(i,j,k)
      hSw_xl = fxn_z0(i,j)                          ! уровень поверхности МВ на грани

      ! на грани X-right (t=t[n+1]):
      dz_xr = fxn_z(i+1,j,k) - fxn_z(i+1,j,k+1)
      teta_xr = fxn_teta(i+1,j,k)
      rho_xr = fxn_rho(i+1,j,k)
      dp_xr = sound0_2 * fxn_dteta(i+1,j,k)
      u_xr = fxn_u(i+1,j,k)
      v_xr = fxn_v(i+1,j,k)
      w_xr = fxn_w(i+1,j,k)
      hSw_xr = fxn_z0(i+1,j)                          ! уровень поверхности МВ на грани

      ! на грани Y-left (t=t[n+1]):
      dz_yl = fyn_z(i,j,k) - fyn_z(i,j,k+1)
      teta_yl = fyn_teta(i,j,k)
      rho_yl = fyn_rho(i,j,k)
      dp_yl = sound0_2 * fyn_dteta(i,j,k)
      u_yl = fyn_u(i,j,k)
      v_yl = fyn_v(i,j,k)
      w_yl = fyn_w(i,j,k)
      hSw_yl = fyn_z0(i,j)                          ! уровень поверхности МВ на грани

      ! на грани Y-right (t=t[n+1]):
      dz_yr = fyn_z(i,j+1,k) - fyn_z(i,j+1,k+1)
      teta_yr = fyn_teta(i,j+1,k)
      rho_yr = fyn_rho(i,j+1,k)
      dp_yr = sound0_2 * fyn_dteta(i,j+1,k)
      u_yr = fyn_u(i,j+1,k)
      v_yr = fyn_v(i,j+1,k)
      w_yr = fyn_w(i,j+1,k)
      hSw_yr = fyn_z0(i,j+1)                          ! уровень поверхности МВ на грани

      ! на грани Z-top (t=t[n+1]):
      dz_xt = fxn_z(i+1,j,k) - fxn_z(i,j,k)
      dz_yt = fyn_z(i,j+1,k) - fyn_z(i,j,k)
      rho_t = fzn_rho(i,j,k)                  ; rho_xt  = rho_t ; rho_yt  = rho_t
      dp_t = sound0_2 * fzn_dteta(i,j,k)      ; dp_xt   = dp_t  ; dp_yt   = dp_t
      u_t = fzn_u(i,j,k)                      ; u_xt    = u_t   ; u_yt    = u_t
      v_t = fzn_v(i,j,k)                      ; v_xt    = v_t   ; v_yt    = v_t
      w_t = fzn_w(i,j,k)                      ; w_xt    = w_t   ; w_yt    = w_t
      teta_t = fzn_teta(i,j,k)                ; teta_xt = teta_t; teta_yt = teta_t
      hSw_t = cn_z0(i,j)                          ! уровень поверхности МВ в ячейке

      ! на грани Z-bottom (t=t[n+1]):
      dz_xb = fxn_z(i+1,j,k+1) - fxn_z(i,j,k+1)
      dz_yb = fyn_z(i,j+1,k+1) - fyn_z(i,j,k+1)
      rho_b = fzn_rho(i,j,k+1)                  ; rho_xb  = rho_b ; rho_yb  = rho_b
      dp_b = sound0_2 * fzn_dteta(i,j,k+1)      ; dp_xb   = dp_b  ; dp_yb   = dp_b
      u_b = fzn_u(i,j,k+1)                      ; u_xb    = u_b   ; u_yb    = u_b
      v_b = fzn_v(i,j,k+1)                      ; v_xb    = v_b   ; v_yb    = v_b
      w_b = fzn_w(i,j,k+1)                      ; w_xb    = w_b   ; w_yb    = w_b
      teta_b = fzn_teta(i,j,k+1)                ; teta_xb = teta_b; teta_yb = teta_b
      hSw_b = hSw_t

      !-- баланcные уравнения -------------------------------------------------------

      ! баланc teta:

      q_xr = teta_xr * u_xr                                     ! грань X-право
      q_xl = teta_xl * u_xl                                     ! грань X-лево
      q_xt = teta_xt * u_xt                                     ! грань Z-верх, производная по X
      q_xb = teta_xb * u_xb                                     ! грань Z-низ, производная по X

      q_yr = teta_yr * v_yr                                     ! грань Y-право
      q_yl = teta_yl * v_yl                                     ! грань Y-лево
      q_yt = teta_yt * v_yt                                     ! грань Z-верх, производная по Y
      q_yb = teta_yb * v_yb                                     ! грань Z-низ, производная по Y

      q_t = teta_t                                              ! грань Z-верх
      q_b = teta_b                                              ! грань Z-низ

      teta_h_new = teta_c * dz_c - ( &
                   ( q_xr * dz_xr - q_xl * dz_xl) / dx + &
                   ( q_xb * dz_xb - q_xt * dz_xt) / dx + &
                   ( q_yr * dz_yr - q_yl * dz_yl) / dy + &
                   ( q_yb * dz_yb - q_yt * dz_yt) / dy + &
                   ( q_t * dwz(k) - q_b * dwz(k+1) ) &
                                    ) * dt / 2.

      ! баланc rho:

      q_xr = teta_xr * u_xr * rho_xr                           ! грань X-право
      q_xl = teta_xl * u_xl * rho_xl                           ! грань X-лево
      q_xt = teta_xt * u_xt * rho_xt                           ! грань Z-верх, производная по X
      q_xb = teta_xb * u_xb * rho_xb                           ! грань Z-низ, производная по X

      q_yr = teta_yr * v_yr * rho_yr                           ! грань Y-право
      q_yl = teta_yl * v_yl * rho_yl                           ! грань Y-лево
      q_yt = teta_yt * v_yt * rho_yt                           ! грань Z-верх, производная по Y
      q_yb = teta_yb * v_yb * rho_yb                           ! грань Z-низ, производная по Y

      q_t = teta_t * rho_t                                     ! грань Z-верх
      q_b = teta_b * rho_b                                     ! грань Z-низ

      rho_h_new  = teta_c * dz_c * rho_c - ( &
                   ( q_xr * dz_xr - q_xl * dz_xl) / dx + &
                   ( q_xb * dz_xb - q_xt * dz_xt) / dx + &
                   ( q_yr * dz_yr - q_yl * dz_yl) / dy + &
                   ( q_yb * dz_yb - q_yt * dz_yt) / dy + &
                   ( q_t * dwz(k) - q_b * dwz(k+1) ) &
                                    ) * dt / 2.

      ! баланc u:

      q_xr = teta_xr * u_xr * u_xr + &                          ! грань X-право
             dp_xr / rho0 + &
             g * hSw_xr
      q_xl = teta_xl * u_xl * u_xl + &                          ! грань X-лево
             dp_xl / rho0 + &
             g * hSw_xl
      q_xt = teta_xt * u_xt * u_xt + &                          ! грань Z-верх, производная по X
             dp_xt / rho0 + &
             g * hSw_t
      q_xb = teta_xb * u_xb * u_xb + &                          ! грань Z-низ, производная по X
             dp_xb / rho0 + &
             g * hSw_b

      q_yr = teta_yr * v_yr * u_yr                              ! грань Y-право
      q_yl = teta_yl * v_yl * u_yl                              ! грань Y-лево
      q_yt = teta_yt * v_yt * u_yt                              ! грань Z-верх, производная по Y
      q_yb = teta_yb * v_yb * u_yb                              ! грань Z-низ, производная по Y

      q_t = teta_t * u_t                                        ! грань Z-верх
      q_b = teta_b * u_b                                        ! грань Z-низ

      u_h_new  = teta_c * dz_c * u_c - ( &
                   ( q_xr * dz_xr - q_xl * dz_xl) / dx + &
                   ( q_xb * dz_xb - q_xt * dz_xt) / dx + &
                   ( q_yr * dz_yr - q_yl * dz_yl) / dy + &
                   ( q_yb * dz_yb - q_yt * dz_yt) / dy + &
                   ( q_t * dwz(k) - q_b * dwz(k+1) ) &
                                    ) * dt / 2.

      ! баланc v:

      q_xr = teta_xr * u_xr * v_xr                              ! грань X-право
      q_xl = teta_xl * u_xl * v_xl                              ! грань X-лево
      q_xt = teta_xt * u_xt * v_xt                              ! грань Z-верх, производная по X
      q_xb = teta_xb * u_xb * v_xb                              ! грань Z-низ, производная по X

      q_yr = teta_yr * v_yr * v_yr + &                          ! грань Y-право
             dp_yr / rho0 + &
             g * hSw_yr
      q_yl = teta_yl * v_yl * v_yl + &                          ! грань Y-лево
             dp_yl / rho0 + &
             g * hSw_yl
      q_yt = teta_yt * v_yt * v_yt + &                          ! грань Z-верх, производная по Y
             dp_yt / rho0 + &
             g * hSw_t
      q_yb = teta_yb * v_yb * v_yb + &                          ! грань Z-низ, производная по Y
             dp_yb / rho0 + &
             g * hSw_b

      q_t = teta_t * v_t                                        ! грань Z-верх
      q_b = teta_b * v_b                                        ! грань Z-низ

      v_h_new  = teta_c * dz_c * v_c - ( &
                   ( q_xr * dz_xr - q_xl * dz_xl) / dx + &
                   ( q_xb * dz_xb - q_xt * dz_xt) / dx + &
                   ( q_yr * dz_yr - q_yl * dz_yl) / dy + &
                   ( q_yb * dz_yb - q_yt * dz_yt) / dy + &
                   ( q_t * dwz(k) - q_b * dwz(k+1) ) &
                                    ) * dt / 2.

      ! баланc w:

      q_xr = teta_xr * u_xr * w_xr                              ! грань X-право
      q_xl = teta_xl * u_xl * w_xl                              ! грань X-лево
      q_xt = teta_xt * u_xt * w_xt                              ! грань Z-верх, производная по X
      q_xb = teta_xb * u_xb * w_xb                              ! грань Z-низ, производная по X

      q_yr = teta_yr * v_yr * w_yr                              ! грань Y-право
      q_yl = teta_yl * v_yl * w_yl                              ! грань Y-лево
      q_yt = teta_yt * v_yt * w_yt                              ! грань Z-верх, производная по Y
      q_yb = teta_yb * v_yb * w_yb                              ! грань Z-низ, производная по Y

      q_t = teta_t * w_t                                        ! грань Z-верх
      q_b = teta_b * w_b                                        ! грань Z-низ

      w_h_new  = teta_c * dz_c * w_c - ( &
                   ( q_xr * dz_xr - q_xl * dz_xl) / dx + &
                   ( q_xb * dz_xb - q_xt * dz_xt) / dx + &
                   ( q_yr * dz_yr - q_yl * dz_yl) / dy + &
                   ( q_yb * dz_yb - q_yt * dz_yt) / dy + &
                   ( (q_t * dwz(k) + dp_t / rho0) - (q_b * dwz(k+1) + dp_b * rho0) ) + &
                   g * (rho_c - rho0) / rho0 * dz_c &
                                    ) * dt / 2.

      !-- новые консервативные переменные ----------------------------------------------------------------------------

      c_teta(i,j,k) = teta_h_new / (fzn_z(i,j,k) - fzn_z(i,j,k+1))
      c_rho(i,j,k) = rho_h_new / teta_h_new
      c_u(i,j,k) = u_h_new / teta_h_new
      c_v(i,j,k) = v_h_new / teta_h_new
      c_w(i,j,k) = w_h_new / teta_h_new

      c_dteta(i,j,k) = c_teta(i,j,k) - teta0
      c_drho(i,j,k) = c_rho(i,j,k) - rho0
      c_du(i,j,k) = c_u(i,j,k) - cn_u0(i,j)
      c_dv(i,j,k) = c_v(i,j,k) - cn_v0(i,j)
      c_dw(i,j,k) = c_w(i,j,k) - cn_w0(i,j,k)

      if(IsDebCell(i, j, k)) then
        write(17,"('side-XL:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, teta_xl, rho_xl, u_xl, v_xl, w_xl, dp_xl, hSw_xl, dz_xl; flush(17)
        write(17,"('side-XR:',3(';',i0),100(';',1p,g0.16))") &
            i+1, j, k, teta_xr, rho_xr, u_xr, v_xr, w_xr, dp_xr, hSw_xr, dz_xr; flush(17)
        write(17,"('side-YL:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, teta_yl, rho_yl, u_yl, v_yl, w_yl, dp_yl, hSw_yl, dz_yl; flush(17)
        write(17,"('side-YR:',3(';',i0),100(';',1p,g0.16))") &
            i, j+1, k, teta_yr, rho_yr, u_yr, v_yr, w_yr, dp_yr, hSw_yr, dz_yr; flush(17)
        write(17,"('side-B:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k, teta_b, rho_b, u_b, v_b, w_b, dp_b, hSw_b, dz_xb, dz_yb; flush(17)
        write(17,"('side-T:',3(';',i0),100(';',1p,g0.16))") &
            i, j, k+1, teta_t, rho_t, u_t, v_t, w_t, dp_t, hSw_t, dz_xt, dz_yt; flush(17)
        write(17,"('UCS:' ,3(';',i0),100(';',1p,g0.16))") i, j, k, cs_teta(i,j,k), cs_rho(i,j,k), cs_u(i,j,k), cs_v(i,j,k), cs_w(i,j,k)
        write(17,"('UC:',3(';',i0),100(';',1p,g0.16))") i, j, k,  c_teta(i,j,k),  c_rho(i,j,k),  c_u(i,j,k),  c_v(i,j,k),  c_w(i,j,k)
        flush(17)
      end if

    end do
  end do

end
