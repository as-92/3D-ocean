! настройка ГУ, переменных во времени
subroutine SetupVariableBCs
  use variables
  use datetime_module

  type(datetime) :: timeBase, startTime
  type(timedelta) :: tdelta

  real(R8), parameter:: om = 2. * 360. / (60. * 60. * 24.)
  integer:: i, j, k, is, ns, it, ibc
  real(R8), save:: delta0 = -1                                    ! с такой разницы стартуем
  real(R8):: delta, deltag, timh, psx, psy, pax, pay, pbx, pby
  real(R8):: dsax, dsay, dbax, dbay
  real(R8):: ma, mp, sa, sp, f(3)

  ! для задачи "Белое море" вычисляем высоту прилива на крайних точках Горла:
  if(taskNum==200) then

    ! timh=3216. + timestep(ns)/3600.
    ! уравнение сделано по Гринвичу
    ! 3216 (часы) - прямой сдвиг от начала года на 00 часов 15 мая т.е. на 134 дня

    ! этот кусок кода исполняется один раз для задачи 200:
    if(delta0<0) then                                             ! здесь первый раз, вычисляем стартовую разницу во времени

      ! 1 января 2005, полночь по Гринвичу:
      timeBase = datetime(2005,1, 1,0,0,0) ! 2005-01-01 00:00:00
      ! стартовая дата расчёта:
      timeBase = datetime(2005,5,15,0,0,0) ! 2005-05-15 00:00:00

      ! дельта времени от даты расчёта до стартовой даты (сек)
      tdelta = startTime - timeBase
      delta0 = tdelta.total_seconds()

      ! Известны высоты уровня на концах отрезка. Концы отрезка заданы координатами gspos(1) и gspos(2). Для каждой из граней
      ! с типом ГУ BC_H_VAR_T необходимо найти проекцию на этот отрезок и запомнить коэффициент 0<ksi<1. В дальнейшем для таких граней будет
      ! выполняться линейная интерполяция: h[side] = dzg0 + ksi * (dzg1-dzg0)

      ns = max(nfx_sides(BC_H_VAR_T), nfy_sides(BC_H_VAR_T))    ! количество граней с ГУ 'BC_H_VAR_T'
      if(ns>0) then
        allocate(ksiWS(ns,2))                                   ! по количеству граней
        needCalc_BC_H_VAR_T = .true.

        ns = nfx_sides(BC_H_VAR_T)                              ! количество граней
        do is=1, ns
          i = fx_sides(BC_H_VAR_T).ptr(is).i                    ! индексы грани на сетке
          j = fx_sides(BC_H_VAR_T).ptr(is).j

          psx = x(i); psy = (y(j)+y(j+1)) / 2.                  ! координата центра грани (точка S)
          pax = gspos(1).x; pay = gspos(1).y                      ! координата опорной точки А
          pbx = gspos(2).x; pby = gspos(2).y                      ! координата опорной точки В
          ! ищем проекцию pk точки ps на отрезок (А,В):
          dsax = psx - pax; dsay = psy - pay                    ! вектор AS
          dbax = pbx - pax; dbay = pby - pay                    ! вектор AB
          ksiWS(is,1) = (dsax*dbax + dsay*dbay) / (dbax*dbax + dbay*dbay)
          call assertn(ksiWS(is,1)>=0 .and. ksiWS(is,1)<=1., "SetupVariableBCs-ksi-X", i, j, 0)
        end do

        ns = nfy_sides(BC_H_VAR_T)                              ! количество граней
        do is=1, ns
          i = fy_sides(BC_H_VAR_T).ptr(is).i                    ! индексы грани на сетке
          j = fy_sides(BC_H_VAR_T).ptr(is).j
          psx = (x(i)+x(i+1)) / 2.; psy = y(j)                  ! координата центра грани (точка S)
          pax = gspos(1).x; pay = gspos(1).y                      ! координата опорной точки А
          pbx = gspos(2).x; pby = gspos(2).y                      ! координата опорной точки В
          ! ищем проекцию pk точки ps на отрезок (А,В):
          dsax = psx - pax; dsay = psy - pay                    ! вектор AS
          dbax = pbx - pax; dbay = pby - pay                    ! вектор AB
          ksiWS(is,2) = (dsax*dbax + dsay*dbay) / (dbax*dbax + dbay*dbay)
          call assertn(ksiWS(is,2)>=0 .and. ksiWS(is,2)<=1., "SetupVariableBCs-ksi-Y", i, j, 0)
        end do

      end if

    end if      ! if(delta0<0)

    ! для задачи 200 (БМ) каждый шаг:

    delta = delta0 + GetCTime(ttime)                            ! текущая разница во времени (сек)
    deltag = delta + 4. * 3600.                                 ! + 4 часа разницы с Гринвичем

    if(needCalc_BC_H_VAR_T) then                                ! если применяется метод Семёнова
      timh = delta / 3600.                                      ! часы от начала 2005 года (=3216 для 15.05.2005 00:00:00)

      ! отклонение от равновесия в двух береговых точках на время расчётное ttime:
      dzg0 = .8088154793d0 * 5.1d0 * cos(GradToRad((13.9430356d0 * timh) + 274.2344961d0 - 216.8d0)) + &
          .8834757209d0 * 14.8d0 * cos(GradToRad((15.0410686d0 * timh) + 187.4798548d0 - 330.1d0)) + &
          1.037439823d0 * 209.7d0 * cos(GradToRad((28.9841042d0 * timh) + 100.8260815d0 - 334.0d0)) + &
          1.d0 * 58.9d0 * cos(GradToRad((30.0410686d0 * timh) - 28.6d0))
      dzg0 = -dzg0 / 100.                                       ! [см] -> [м]

      dzg1 = .8088154793d0 * 1.4d0 * cos(GradToRad((13.9430356d0 * timh) + 274.2344961d0 - 272.7d0)) + &
          .8834757209d0 * 9.7d0 * cos(GradToRad((15.0410686d0 * timh) + 187.4798548d0 - 344.1d0)) + &
          1.037439823d0 * 237.9d0 * cos(GradToRad((28.9841042d0 * timh) + 100.8260815d0 - 358.9d0)) + &
          1. * 45.3 * cos(GradToRad((30.0410686 * timh) - 62.0))
      dzg1 = -dzg1 / 100.                                       ! [см] -> [м]
    end if
  end if

  ! для каждой грани с переменными по времени ГУ вычисляем параметры:

  !-- вход ---------------------------------------------------------------

  ns = nfx_sides(BC_IN_T)
  do is=1,ns
    i = fx_sides(BC_IN_T).ptr(is).i                             ! индексы грани на сетке
    j = fx_sides(BC_IN_T).ptr(is).j
    ibc = fx_bc(i,j)
    it = bcData(ibc).ind                                        ! номер топекса, соответствующего грани

    ! читаем данные из топекса (3 параметра):
    do k=1,3
      ma = topex(it).mAmp(k) * 1.d-2                            ! амплитуды в TOPEX указаны в см/сек
      mp = topex(it).mPh(k)
      sa = topex(it).sAmp(k) * 1.e-2
      sp = topex(it).sPh(k)
      f(k) = ma * cos(GradToRad(om*delta - mp)) + sa * cos(GradToRad(om*delta - sp))
    end do

    ! расчитанные данные кладём в ГУ грани:
    bcData(ibc).h = -fx_b(i,j) + Get0Depth(f(1))                ! глубина с возмущением -> безразмерная
    bcData(ibc).u = Get0Velocity(f(2))
    bcData(ibc).v = Get0Velocity(f(3))
    bcData(ibc).w = 0.
    bcData(ibc).rho = rho0
    bcData(ibc).teta = teta0
  end do      ! is

  ns = nfy_sides(BC_IN_T)
  do is=1,ns
    i = fy_sides(BC_IN_T).ptr(is).i                             ! индексы грани на сетке
    j = fy_sides(BC_IN_T).ptr(is).j
    ibc = fy_bc(i,j)
    it = bcData(ibc).ind                                        ! номер топекса, соответствующего грани

    ! читаем данные из топекса (3 параметра):
    do k=1,3
      ma = topex(it).mAmp(k) * 1.d-2                            ! амплитуды в TOPEX указаны в см/сек
      mp = topex(it).mPh(k)
      sa = topex(it).sAmp(k) * 1.e-2
      sp = topex(it).sPh(k)
      f(k) = ma * cos(GradToRad(om*delta - mp)) + sa * cos(GradToRad(om*delta - sp))
    end do

    ! расчитанные данные кладём в ГУ грани:
    bcData(ibc).h = -fy_b(i,j) + Get0Depth(f(1))                ! глубина с возмущением -> безразмерная
    bcData(ibc).u = Get0Velocity(f(2))
    bcData(ibc).v = Get0Velocity(f(3))
    bcData(ibc).w = 0.
    bcData(ibc).rho = rho0
    bcData(ibc).teta = teta0
  end do      ! is

  !-- заданная глубина -----------------------------------------------------------------------

  !
  !    Н Е    Р Е А Л И З О В А Н О
  !

  call assert(nfx_sides(BC_H_VAR_T)==0 .and. nfy_sides(BC_H_VAR_T)==0, &
        "SetupVariableBCs - НЕ РЕАЛИЗОВАНО")

end

