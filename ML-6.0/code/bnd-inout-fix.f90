! ГУ "постоянный вток" на X-гранях:
subroutine fx_BndInoutFix
  use variables
  implicit none

  integer i, j, k, is, ic, ns

  ns = nfx_sides(BC_IN)                                         ! количество граней
  do is=1, ns
    i = fx_sides(BC_IN).ptr(is).i                               ! индексы грани на сетке
    j = fx_sides(BC_IN).ptr(is).j

    call assertn(.false., "Phase2_X. Не реализовано BC_IN", -1, -1, -1)

  end do  ! do is

end subroutine fx_BndInoutFix

!===========================================================================================================

subroutine fy_BndInoutFix
  use variables
  implicit none

  integer i, j, k, is, jc, ns

  ns = nfx_sides(BC_IN)                                         ! количество граней
  do is=1, ns
    i = fx_sides(BC_IN).ptr(is).i                               ! индексы грани на сетке
    j = fx_sides(BC_IN).ptr(is).j

    call assertn(.false., "Phase2_Y. Не реализовано BC_IN", -1, -1, -1)

  end do  ! do is

end subroutine fy_BndInoutFix
