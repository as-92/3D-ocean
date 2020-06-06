! ГУ "фиксированная скорость"
subroutine fx_BndFixU
  use variables
  implicit none

  integer i, j, k, is, ic, ns

  ns = nfx_sides(BC_FIX_U)                                       ! количество граней
  do is=1, ns
    i = fx_sides(BC_FIX_U).ptr(is).i                               ! индексы грани на сетке
    j = fx_sides(BC_FIX_U).ptr(is).j

    call assertn(.false., "Phase2_X. Не реализовано BC_FIX_U", -1, -1, -1)

  end do  ! do is

end subroutine fx_BndFixU

!===========================================================================================================

subroutine fy_BndFixU
  use variables
  implicit none

  integer i, j, k, is, jc, ns

  ns = nfx_sides(BC_FIX_U)                                         ! количество граней
  do is=1, ns
    i = fx_sides(BC_FIX_U).ptr(is).i                               ! индексы грани на сетке
    j = fx_sides(BC_FIX_U).ptr(is).j

    call assertn(.false., "Phase2_Y. Не реализовано BC_FIX_U", -1, -1, -1)

  end do  ! do is

end subroutine fy_BndFixU
