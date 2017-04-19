module gsw_mod
  !! a stripped down copy of the Gibbs-SeaWater (GSW) Oceanographic Toolbox, using just the
  !! density calculation function [[gsw_rho]]

  implicit none
  private

  ! public methods
  !------------------------------------------------------------
  public :: gsw_rho


  ! private module variables
  !------------------------------------------------------------
  integer,   parameter :: r8 = selected_real_kind(14,30)
  real (r8), parameter :: gsw_cp0 = 3991.86795711963_r8
  real (r8), parameter :: gsw_sfac = 0.0248826675584615_r8



contains



  !================================================================================
  !================================================================================



  elemental function gsw_ct_from_pt (sa, pt)
    !! Calculates Conservative Temperature from potential temperature of seawater
    implicit none

    real (r8), intent(in) :: sa
    !! Absolute Salinity [g/kg]

    real (r8), intent(in) :: pt
    !! potential temperature with reference pressure of 0 dbar [deg C]

    real (r8) :: gsw_ct_from_pt
    !! Conservative Temperature  [deg C]

    real (r8) :: pot_enthalpy, x2, x, y

    x2 = gsw_sfac*sa
    x = sqrt(x2)
    y = pt*0.025_r8        ! normalize for F03 and F08

    pot_enthalpy =  61.01362420681071_r8 + y*(168776.46138048015_r8 + &
         y*(-2735.2785605119625_r8 + y*(2574.2164453821433_r8 + &
         y*(-1536.6644434977543_r8 + y*(545.7340497931629_r8 + &
         (-50.91091728474331_r8 - 18.30489878927802_r8*y)*y))))) + &
         x2*(268.5520265845071_r8 + y*(-12019.028203559312_r8 + &
         y*(3734.858026725145_r8 + y*(-2046.7671145057618_r8 + &
         y*(465.28655623826234_r8 + (-0.6370820302376359_r8 - &
         10.650848542359153_r8*y)*y)))) + &
         x*(937.2099110620707_r8 + y*(588.1802812170108_r8 + &
         y*(248.39476522971285_r8 + (-3.871557904936333_r8 - &
         2.6268019854268356_r8*y)*y)) + &
         x*(-1687.914374187449_r8 + x*(246.9598888781377_r8 + &
         x*(123.59576582457964_r8 - 48.5891069025409_r8*x)) + &
         y*(936.3206544460336_r8 + &
         y*(-942.7827304544439_r8 + y*(369.4389437509002_r8 + &
         (-33.83664947895248_r8 - 9.987880382780322_r8*y)*y))))))

    gsw_ct_from_pt = pot_enthalpy/gsw_cp0

    return
  end function gsw_ct_from_pt



  !================================================================================
  !================================================================================



  elemental function gsw_rho (s, pt, p)
    !!  Calculates specific volume from Absolute Salinity, Conservative
    !!  Temperature and pressure, using the computationally-efficient
    !!  polynomial expression for specific volume (Roquet et al., 2014).

    implicit none
    real (r8), intent(in) :: s, pt, p
    real (r8) :: sa
    real (r8) :: gsw_rho
    real (r8) :: ct
    real (r8) :: gsw_specvol
    real (r8) :: xs, ys, z

    !! Absolute Salinity                                        [ g/kg ]
    !  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
    !  p   =  sea pressure                                             [ dbar ]
    !         ( i.e. absolute pressure - 10.1325 dbar )
    !
    !  specvol  =  specific volume                                   [ m^3/kg ]
    !--------------------------------------------------------------------------

    real (r8), parameter :: offset = 5.971840214030754e-1_r8

    real (r8), parameter :: v000 =  1.0769995862e-3_r8
    real (r8), parameter :: v001 = -6.0799143809e-5_r8
    real (r8), parameter :: v002 =  9.9856169219e-6_r8
    real (r8), parameter :: v003 = -1.1309361437e-6_r8
    real (r8), parameter :: v004 =  1.0531153080e-7_r8
    real (r8), parameter :: v005 = -1.2647261286e-8_r8
    real (r8), parameter :: v006 =  1.9613503930e-9_r8
    real (r8), parameter :: v010 = -3.1038981976e-4_r8
    real (r8), parameter :: v011 =  2.4262468747e-5_r8
    real (r8), parameter :: v012 = -5.8484432984e-7_r8
    real (r8), parameter :: v013 =  3.6310188515e-7_r8
    real (r8), parameter :: v014 = -1.1147125423e-7_r8
    real (r8), parameter :: v020 =  6.6928067038e-4_r8
    real (r8), parameter :: v021 = -3.4792460974e-5_r8
    real (r8), parameter :: v022 = -4.8122251597e-6_r8
    real (r8), parameter :: v023 =  1.6746303780e-8_r8
    real (r8), parameter :: v030 = -8.5047933937e-4_r8
    real (r8), parameter :: v031 =  3.7470777305e-5_r8
    real (r8), parameter :: v032 =  4.9263106998e-6_r8
    real (r8), parameter :: v040 =  5.8086069943e-4_r8
    real (r8), parameter :: v041 = -1.7322218612e-5_r8
    real (r8), parameter :: v042 = -1.7811974727e-6_r8
    real (r8), parameter :: v050 = -2.1092370507e-4_r8
    real (r8), parameter :: v051 =  3.0927427253e-6_r8
    real (r8), parameter :: v060 =  3.1932457305e-5_r8
    real (r8), parameter :: v100 = -1.5649734675e-5_r8
    real (r8), parameter :: v101 =  1.8505765429e-5_r8
    real (r8), parameter :: v102 = -1.1736386731e-6_r8
    real (r8), parameter :: v103 = -3.6527006553e-7_r8
    real (r8), parameter :: v104 =  3.1454099902e-7_r8
    real (r8), parameter :: v110 =  3.5009599764e-5_r8
    real (r8), parameter :: v111 = -9.5677088156e-6_r8
    real (r8), parameter :: v112 = -5.5699154557e-6_r8
    real (r8), parameter :: v113 = -2.7295696237e-7_r8
    real (r8), parameter :: v120 = -4.3592678561e-5_r8
    real (r8), parameter :: v121 =  1.1100834765e-5_r8
    real (r8), parameter :: v122 =  5.4620748834e-6_r8
    real (r8), parameter :: v130 =  3.4532461828e-5_r8
    real (r8), parameter :: v131 = -9.8447117844e-6_r8
    real (r8), parameter :: v132 = -1.3544185627e-6_r8
    real (r8), parameter :: v140 = -1.1959409788e-5_r8
    real (r8), parameter :: v141 =  2.5909225260e-6_r8
    real (r8), parameter :: v150 =  1.3864594581e-6_r8
    real (r8), parameter :: v200 =  2.7762106484e-5_r8
    real (r8), parameter :: v201 = -1.1716606853e-5_r8
    real (r8), parameter :: v202 =  2.1305028740e-6_r8
    real (r8), parameter :: v203 =  2.8695905159e-7_r8
    real (r8), parameter :: v210 = -3.7435842344e-5_r8
    real (r8), parameter :: v211 = -2.3678308361e-7_r8
    real (r8), parameter :: v212 =  3.9137387080e-7_r8
    real (r8), parameter :: v220 =  3.5907822760e-5_r8
    real (r8), parameter :: v221 =  2.9283346295e-6_r8
    real (r8), parameter :: v222 = -6.5731104067e-7_r8
    real (r8), parameter :: v230 = -1.8698584187e-5_r8
    real (r8), parameter :: v231 = -4.8826139200e-7_r8
    real (r8), parameter :: v240 =  3.8595339244e-6_r8
    real (r8), parameter :: v300 = -1.6521159259e-5_r8
    real (r8), parameter :: v301 =  7.9279656173e-6_r8
    real (r8), parameter :: v302 = -4.6132540037e-7_r8
    real (r8), parameter :: v310 =  2.4141479483e-5_r8
    real (r8), parameter :: v311 = -3.4558773655e-6_r8
    real (r8), parameter :: v312 =  7.7618888092e-9_r8
    real (r8), parameter :: v320 = -1.4353633048e-5_r8
    real (r8), parameter :: v321 =  3.1655306078e-7_r8
    real (r8), parameter :: v330 =  2.2863324556e-6_r8
    real (r8), parameter :: v400 =  6.9111322702e-6_r8
    real (r8), parameter :: v401 = -3.4102187482e-6_r8
    real (r8), parameter :: v402 = -6.3352916514e-8_r8
    real (r8), parameter :: v410 = -8.7595873154e-6_r8
    real (r8), parameter :: v411 =  1.2956717783e-6_r8
    real (r8), parameter :: v420 =  4.3703680598e-6_r8
    real (r8), parameter :: v500 = -8.0539615540e-7_r8
    real (r8), parameter :: v501 =  5.0736766814e-7_r8
    real (r8), parameter :: v510 = -3.3052758900e-7_r8
    real (r8), parameter :: v600 =  2.0543094268e-7_r8

    sa = (35.16504 / 35.0) * s

    ct = pt
    ct = gsw_ct_from_pt (sa, pt)

    xs = sqrt(gsw_sfac*sa + offset)
    ys = ct*0.025_r8
    z = p*1e-4_r8

    gsw_specvol = v000 + xs*(v010 + xs*(v020 + xs*(v030 + xs*(v040 + xs*(v050 &
         + v060*xs))))) + ys*(v100 + xs*(v110 + xs*(v120 + xs*(v130 + xs*(v140 &
         + v150*xs)))) + ys*(v200 + xs*(v210 + xs*(v220 + xs*(v230 + v240*xs))) &
         + ys*(v300 + xs*(v310 + xs*(v320 + v330*xs)) + ys*(v400 + xs*(v410 &
         + v420*xs) + ys*(v500 + v510*xs + v600*ys))))) + z*(v001 + xs*(v011 &
         + xs*(v021 + xs*(v031 + xs*(v041 + v051*xs)))) + ys*(v101 + xs*(v111 &
         + xs*(v121 + xs*(v131 + v141*xs))) + ys*(v201 + xs*(v211 + xs*(v221 &
         + v231*xs)) + ys*(v301 + xs*(v311 + v321*xs) + ys*(v401 + v411*xs &
         + v501*ys)))) + z*(v002 + xs*(v012 + xs*(v022 + xs*(v032 + v042*xs))) &
         + ys*(v102 + xs*(v112 + xs*(v122 + v132*xs)) + ys*(v202 + xs*(v212 &
         + v222*xs) + ys*(v302 + v312*xs + v402*ys))) + z*(v003 + xs*(v013 &
         + v023*xs) + ys*(v103 + v113*xs + v203*ys) + z*(v004 + v014*xs + v104*ys &
         + z*(v005 + v006*z)))))
    gsw_rho = 1.0_r8/gsw_specvol
    return
  end function



  !================================================================================

end module gsw_mod
