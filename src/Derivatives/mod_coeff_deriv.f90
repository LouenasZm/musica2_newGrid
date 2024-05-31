!==============================================================================
module mod_coeff_deriv
!==============================================================================
  !> Module for coefficients of derivative schemes
!==============================================================================
  use precision
  implicit none
  ! ---------------------------------------------------------------------------
  ! Coefficients for first derivatives
  real(wp) :: a3(-1:1),a5(-2:2),a7(-3:3),a9(-4:4),a11(-5:5)
  real(wp), dimension(2)  :: a01
  real(wp), dimension(3)  :: a02,a20
  real(wp), dimension(5)  :: a13,a04,a31,a40
  real(wp), dimension(7)  :: a24,a15,a06,a42,a51,a60
  real(wp), dimension(9)  :: a35,a26,a17,a08
  real(wp), dimension(11) :: a46,a37,a28,a19,a010
  ! ---------------------------------------------------------------------------
  ! Coefficients for second derivatives
  real(wp) :: b3(-1:1),b5(-2:2)
  real(wp), dimension(4)  :: b03,b30
  real(wp), dimension(5)  :: b13,b04,b31,b40
  ! ---------------------------------------------------------------------------
  ! Coefficients SBP (Summation By Parts)
  real(wp) :: as4(-2:2)
  real(wp) :: as4p0(4),as4p1(4),as4p2(5),as4p3(6)
  real(wp) :: as4m0(4),as4m1(4),as4m2(5),as4m3(6)
  real(wp) :: as6(-3:3)
  real(wp) :: as6p0(6),as6p1(6),as6p2(6),as6p3(7),as6p4(8),as6p5(9)
  real(wp) :: as6m0(6),as6m1(6),as6m2(6),as6m3(7),as6m4(8),as6m5(9)
  !real(wp) :: h11,h22,h33,h44,h55,h66
  ! ---------------------------------------------------------------------------
  ! Extrapolation coefficients for bc_wall (in cartesian)
  real(wp), dimension(3,2)  :: cextrp2,cextrp3
  ! ---------------------------------------------------------------------------
  integer, parameter :: iFD=0
  
contains

  !===============================================================
  subroutine init_coeff_deriv(stencil,is_DRP)
  !===============================================================
    !> Initialize derivative scheme coefficients
  !===============================================================
    use warnstop
    implicit none
    ! ------------------------------------------------------------
    ! Scheme parameters
    integer, intent(in)  :: stencil
    logical, intent(in)  :: is_DRP
    ! ------------------------------------------------------------

    ! Second-order derivatives
    ! ========================
    
    ! Standard schemes on 3-point stencil
    ! -----------------------------------
    ! order 2 (centered)
    b3( 1)= 1.0_wp
    b3( 0)=-2.0_wp
    b3(-1)= b3(1)
    
    ! order 2 (decentered 0-3)
    b03(1)= 2.0_wp
    b03(2)=-5.0_wp
    b03(3)= 4.0_wp
    b03(4)=-1.0_wp    
    ! order 2 (decentered 3-0)
    b30 = -b03

    ! Standard schemes on 5-point stencil
    ! -----------------------------------
    ! order 4 (centered)
    b5( 2)=-1.0_wp/12.0_wp
    b5( 1)= 4.0_wp/3.0_wp
    b5( 0)=-5.0_wp/2.0_wp
    b5(-1)= b5(1)
    b5(-2)= b5(2)

    ! order 4 (decentered 1-3)
    b13(1)=11.0_wp/12.0_wp
    b13(2)=-5.0_wp/3.0_wp
    b13(3)= 1.0_wp/2.0_wp
    b13(4)= 1.0_wp/3.0_wp
    b13(5)=-1.0_wp/12.0_wp
    ! order 4 (decentered 3-1)
    b31 = -b13

    ! order 4 (decentered 0-4)
    b04(1)= 35.0_wp/12.0_wp
    b04(2)=-26.0_wp/3.0_wp
    b04(3)= 19.0_wp/2.0_wp
    b04(4)=-14.0_wp/3.0_wp
    b04(5)= 11.0_wp/12.0_wp
    ! order 4 (decentered 4-0)    
    b40 = -b04

    ! First-order derivatives
    ! =======================

    ! Boundary scheme o1 / 2pts
    ! -------------------------
    a01(1)=-1.0_wp
    a01(2)= 1.0_wp

    ! Standard schemes on 3-point stencil
    ! -----------------------------------
    ! order 2 (centered)
    a3( 1)=  0.5_wp
    a3( 0)=  0.0_wp
    a3(-1)= -a3(1)

    ! order 2 (decentered 0-2)
    a02(1)= -3.0_wp/2.0_wp
    a02(2)=  2.0_wp
    a02(3)= -1.0_wp/2.0_wp
    ! order 2 (decentered 2-0)
    a20 = -a02

    ! Standard schemes on 5-point stencil
    ! -----------------------------------
    ! order 4 (centered)
    a5( 2)= -1.0_wp/12.0_wp
    a5( 1)=  2.0_wp/3.0_wp
    a5( 0)=  0.0_wp
    a5(-1)= -a5(1)
    a5(-2)= -a5(2)

    ! order 4 (decentered 1-3)
    a13(1)= -1.0_wp/4.0_wp
    a13(2)= -5.0_wp/6.0_wp
    a13(3)=  3.0_wp/2.0_wp
    a13(4)= -1.0_wp/2.0_wp
    a13(5)=  1.0_wp/12.0_wp
    ! order 4 (decentered 3-1)
    a31 = -a13

    ! order 4 (decentered 0-4)
    a04(1)= -25.0_wp/12.0_wp
    a04(2)=  4.0_wp
    a04(3)= -3.0_wp
    a04(4)=  4.0_wp/3.0_wp
    a04(5)= -1.0_wp/4.0_wp
    ! order 4 (decentered 4-0)    
    a40 = -a04

    if (stencil<7) return
    
    ! Standard schemes on 7-point stencil
    ! -----------------------------------
    ! order 6 (centered)
    a7( 1) =  3.0_wp/4.0_wp
    a7( 2) = -3.0_wp/20.0_wp
    a7( 3) =  1.0_wp/60.0_wp
    a7( 0) =  0.0_wp
    a7(-1) = -a7(1)
    a7(-2) = -a7(2)
    a7(-3) = -a7(3)

    ! order 6 (decentered 2-4)
    a24(1) =  1.0_wp/30.0_wp
    a24(2) = -2.0_wp/5.0_wp
    a24(3) = -7.0_wp/12.0_wp
    a24(4) =  4.0_wp/3.0_wp
    a24(5) = -1.0_wp/2.0_wp
    a24(6) =  2.0_wp/15.0_wp
    a24(7) = -1.0_wp/60.0_wp
    ! order 6 (decentered 4-2)
    a42 = -a24

    ! order 6 (decentered 1-5)
    a15(1) = -1.0_wp/6.0_wp
    a15(2) =-77.0_wp/60.0_wp
    a15(3) =  5.0_wp/2.0_wp
    a15(4) = -5.0_wp/3.0_wp
    a15(5) =  5.0_wp/6.0_wp
    a15(6) = -1.0_wp/4.0_wp
    a15(7) =  1.0_wp/30.0_wp
    ! order 6 (decentered 5-1)
    a51 = -a15
 
    ! order 6 (decentered 0-6)
    a06(1) =-49.0_wp/20.0_wp
    a06(2) =  6.0_wp
    a06(3) =-15.0_wp/2.0_wp
    a06(4) = 20.0_wp/3.0_wp
    a06(5) =-15.0_wp/4.0_wp
    a06(6) =  6.0_wp/5.0_wp
    a06(7) = -1.0_wp/6.0_wp
    ! order 6 (decentered 6-0)
    a60 = -a06

    ! if DRP (Dispersion Relation Preserving) and 7-pt stencil,
    ! replace by optimized scheme
    ! ~> Tam's scheme [Tam & Webb JCP 1993]
    if ((stencil==7).and.(is_DRP)) then
       ! Tam's scheme (centered)
       a7( 1) =  0.770882380518_wp
       a7( 2) = -0.166705904415_wp
       a7( 3) =  0.020843142770_wp
       a7( 0) =  0.0_wp
       a7(-1) = -a7(1)
       a7(-2) = -a7(2)
       a7(-3) = -a7(3)
       ! Tam's scheme (decentered 2-4)
       a24(1) =  0.049041958_wp
       a24(2) = -0.468840357_wp
       a24(3) = -0.474760914_wp
       a24(4) =  1.273274737_wp
       a24(5) = -0.518484526_wp
       a24(6) =  0.166138533_wp
       a24(7) = -0.026369431_wp
       ! Tam's scheme (decentered 1-5)
       a15(1) = -0.209337622_wp
       a15(2) = -1.084875676_wp
       a15(3) =  2.147776050_wp
       a15(4) = -1.388928322_wp
       a15(5) =  0.768949766_wp
       a15(6) = -0.281814650_wp
       a15(7) =  0.048230454_wp
       ! Tam's scheme (decentered 0-6)
       a06(1) = -2.192280339_wp
       a06(2) =  4.748611401_wp
       a06(3) = -5.108851915_wp
       a06(4) =  4.461567104_wp
       a06(5) = -2.833498741_wp
       a06(6) =  1.128328861_wp
       a06(7) = -0.203876371_wp
       ! Tam's scheme (decentered 4-2)
       a42 = -a24
       ! Tam's scheme (decentered 5-1)
       a51 = -a15
       ! Tam's scheme (decentered 6-0)
       a60 = -a06
    end if
    
    if (stencil<9) return
    
    ! Standard schemes on 9-point stencil
    ! -----------------------------------
    ! order 8 (centered)
    a9( 1) =  4.0_wp/5.0_wp
    a9( 2) = -1.0_wp/5.0_wp
    a9( 3) =  4.0_wp/105.0_wp
    a9( 4) = -1.0_wp/280.0_wp
    a9( 0) =  0.0_wp
    a9(-1) = -a9(1)
    a9(-2) = -a9(2)
    a9(-3) = -a9(3)
    a9(-4) = -a9(4)

    ! order 8 (decentered 3-5)
    a35(1) = -1.0_wp/168.0_wp
    a35(2) =  1.0_wp/14.0_wp
    a35(3) = -1.0_wp/2.0_wp
    a35(4) = -9.0_wp/20.0_wp
    a35(5) =  5.0_wp/4.0_wp
    a35(6) = -1.0_wp/2.0_wp
    a35(7) =  1.0_wp/6.0_wp
    a35(8) = -1.0_wp/28.0_wp
    a35(9) =  1.0_wp/280.0_wp

    ! order 8 (decentered 2-6)
    a26(1) =  1.0_wp/56.0_wp
    a26(2) = -2.0_wp/7.0_wp
    a26(3) =-19.0_wp/20.0_wp
    a26(4) =  2.0_wp/1.0_wp
    a26(5) = -5.0_wp/4.0_wp
    a26(6) =  2.0_wp/3.0_wp
    a26(7) = -1.0_wp/4.0_wp
    a26(8) =  2.0_wp/35.0_wp
    a26(9) = -1.0_wp/168.0_wp

    ! order 8 (decentered 1-7)
    a17(1) = -1.0_wp/8.0_wp
    a17(2) =-223.0_wp/140.0_wp
    a17(3) =  7.0_wp/2.0_wp
    a17(4) = -7.0_wp/2.0_wp
    a17(5) = 35.0_wp/12.0_wp
    a17(6) = -7.0_wp/4.0_wp
    a17(7) =  7.0_wp/10.0_wp
    a17(8) = -1.0_wp/6.0_wp
    a17(9) =  1.0_wp/56.0_wp

    ! order 8 (decentered 0-8)
    a08(1) =-761.0_wp/280.0_wp
    a08(2) =  8.0_wp
    a08(3) =-14.0_wp
    a08(4) = 56.0_wp/3.0_wp
    a08(5) =-35.0_wp/2.0_wp
    a08(6) = 56.0_wp/5.0_wp
    a08(7) =-14.0_wp/3.0_wp
    a08(8) =  8.0_wp/7.0_wp
    a08(9) = -1.0_wp/8.0_wp

    ! if DRP (Dispersion Relation Preserving) and 9-pt stencil,
    ! replace by optimized scheme
    ! ~> [Bogey & Bailly JCP 2004] / boundary schemes not defined ~> standard
    if ((stencil==9).and.(is_DRP)) then
       ! 9-points optimized scheme o(4)
       a9(1) =  0.841570125482_wp
       a9(2) = -0.244678631765_wp
       a9(3) =  0.059463584768_wp
       a9(4) = -0.007650904064_wp
       a9(0) =  0.0_wp
       a9(-1) = -a9(1)
       a9(-2) = -a9(2)
       a9(-3) = -a9(3)
       a9(-4) = -a9(4)
    end if
    
    if (stencil<9) return
    
    ! Standard schemes on 11-point stencil
    ! -----------------------------------
    ! order 10 (centered)
    a11( 1) =  5.0_wp/6.0_wp
    a11( 2) = -5.0_wp/21.0_wp
    a11( 3) =  5.0_wp/84.0_wp
    a11( 4) = -5.0_wp/504.0_wp
    a11( 5) =  1.0_wp/1260.0_wp
    a11(-1) = -a11(1)
    a11(-2) = -a11(2)
    a11(-3) = -a11(3)
    a11(-4) = -a11(4)
    a11(-5) = -a11(5)
! Att! decentered order to be recalculated (WRONG COEFF)
!!$    ! order 10 (decentered 4-6)
!!$    a46( 1) =  1.0_wp/840.0_wp
!!$    a46( 2) = -1.0_wp/63.0_wp
!!$    a46( 3) =  3.0_wp/28.0_wp
!!$    a46( 4) = -4.0_wp/7.0_wp
!!$    a46( 5) =-11.0_wp/30.0_wp
!!$    a46( 6) =  6.0_wp/5.0_wp
!!$    a46( 7) = -1.0_wp/2.0_wp
!!$    a46( 8) =  4.0_wp/21.0_wp
!!$    a46( 9) = -3.0_wp/56.0_wp
!!$    a46(10) =  1.0_wp/105.0_wp
!!$    a46(11) = -1.0_wp/1260.0_wp
!!$    ! order 10 (decentered 3-7)
!!$    a37( 1) = -1.0_wp/360.0_wp
!!$    a37( 2) =  1.0_wp/24.0_wp
!!$    a37( 3) = -3.0_wp/8.0_wp
!!$    a37( 4) =-319.0_wp/420.0_wp
!!$    a37( 5) =  7.0_wp/4.0_wp
!!$    a37( 6) =-21.0_wp/20.0_wp
!!$    a37( 7) =  7.0_wp/12.0_wp
!!$    a37( 8) = -1.0_wp/4.0_wp
!!$    a37( 9) =  3.0_wp/40.0_wp
!!$    a37(10) = -1.0_wp/72.0_wp
!!$    a37(11) =  1.0_wp/840.0_wp
!!$    ! order 10 (decentered 2-8)
!!$    a28( 1) =  1.0_wp/90.0_wp
!!$    a28( 2) = -2.0_wp/9.0_wp
!!$    a28( 3) = -341.0_wp/280.0_wp
!!$    a28( 4) =  8.0_wp/3.0_wp
!!$    a28( 5) = -7.0_wp/3.0_wp
!!$    a28( 6) = 28.0_wp/15.0_wp
!!$    a28( 7) = -7.0_wp/6.0_wp
!!$    a28( 8) =  8.0_wp/15.0_wp
!!$    a28( 9) = -1.0_wp/6.0_wp
!!$    a28(10) =  2.0_wp/63.0_wp
!!$    a28(11) = -1.0_wp/360.0_wp
!!$    ! order 10 (decentered 1-9)
!!$    a19( 1) = -1.0_wp/10.0_wp
!!$    a19( 2) =-139.0_wp/76.0_wp
!!$    a19( 3) =  9.0_wp/2.0_wp
!!$    a19( 4) = -6.0_wp
!!$    a19( 5) =  7.0_wp
!!$    a19( 6) =-63.0_wp/10.0_wp
!!$    a19( 7) = 21.0_wp/5.0_wp
!!$    a19( 8) = -2.0_wp
!!$    a19( 9) =  9.0_wp/14.0_wp
!!$    a19(10) = -1.0_wp/8.0_wp
!!$    a19(11) =  1.0_wp/90.0_wp
!!$    ! order 10 (decentered 0-10)
!!$    a010( 1) =-536.0_wp/183.0_wp
!!$    a010( 2) = 10.0_wp
!!$    a010( 3) =-45.0_wp/2.0_wp
!!$    a010( 4) = 40.0_wp
!!$    a010( 5) =-105.0_wp/2.0_wp
!!$    a010( 6) =252.0_wp/5.0_wp
!!$    a010( 7) =-35.0_wp
!!$    a010( 8) =120.0_wp/7.0_wp
!!$    a010( 9) =-45.0_wp/8.0_wp
!!$    a010(10) = 10.0_wp/9.0_wp
!!$    a010(11) = -1.0_wp/10.0_wp
! Att! replace by optimized ones
    ! 11-points optimized scheme decentered 4-6
    ! [Berland, Marsden, Bogey & Bailly JCP 2007]
    a46( 1)= 0.016756572303_wp
    a46( 2)=-0.117478455239_wp
    a46( 3)= 0.411034935097_wp
    a46( 4)=-1.130286765151_wp
    a46( 5)= 0.341435872100_wp
    a46( 6)= 0.556396830543_wp
    a46( 7)=-0.082525734207_wp
    a46( 8)= 0.003565834658_wp
    a46( 9)= 0.001173034777_wp
    a46(10)=-0.000071772607_wp
    a46(11)=-0.000000352272_wp
    ! 11-points optimized scheme decentered 3-7
    ! [Berland, Marsden, Bogey & Bailly JCP 2007]
    a37( 1)=-0.013277273810_wp
    a37( 2)= 0.115976072920_wp
    a37( 3)=-0.617479187931_wp
    a37( 4)=-0.274113948206_wp
    a37( 5)= 1.086208764655_wp
    a37( 6)=-0.402951626982_wp
    a37( 7)= 0.131066986242_wp
    a37( 8)=-0.028154858354_wp
    a37( 9)= 0.002596328316_wp
    a37(10)= 0.000128743150_wp
    a37(11)= 0.0_wp
    ! 11-points optimized scheme decentered 2-8
    ! [Berland, Marsden, Bogey & Bailly JCP 2007]
    a28( 1)= 0.046246319744_wp
    a28( 2)=-0.462989982072_wp
    a28( 3)=-0.459203180244_wp
    a28( 4)= 1.205900619436_wp
    a28( 5)=-0.423956587692_wp
    a28( 6)= 0.102329382027_wp
    a28( 7)=-0.006253229685_wp
    a28( 8)=-0.002025942780_wp
    a28( 9)=-0.000016793609_wp
    a28(10)=-0.000015302561_wp
    a28(11)=-0.000015302561_wp
    ! 11-points optimized scheme decentered 1-9
    ! [Berland, Marsden, Bogey & Bailly JCP 2007]
    a19( 1)=-0.180022054228_wp
    a19( 2)=-1.237550583044_wp
    a19( 3)= 2.484731692990_wp
    a19( 4)=-1.810320814061_wp
    a19( 5)= 1.112990048440_wp
    a19( 6)=-0.481086916514_wp
    a19( 7)= 0.126598690230_wp
    a19( 8)=-0.015510730165_wp
    a19( 9)= 0.000021609059_wp
    a19(10)= 0.000156447571_wp
    a19(11)=-0.000007390277_wp
    ! 11-points optimized scheme decentered 0-10
    ! [Berland, Marsden, Bogey & Bailly JCP 2007]
    a010( 1)=-2.391602219538_wp
    a010( 2)= 5.832490322294_wp
    a010( 3)=-7.650218001182_wp
    a010( 4)= 7.907810563576_wp
    a010( 5)=-5.922599052629_wp
    a010( 6)= 3.071037015445_wp
    a010( 7)=-1.014956769726_wp
    a010( 8)= 0.170022256519_wp
    a010( 9)= 0.002819958377_wp
    a010(10)=-0.004791009708_wp
    a010(11)=-0.000013063429_wp

    ! if DRP (Dispersion Relation Preserving) and 11-pt stencil,
    ! replace by optimized scheme
    if ((stencil==11).and.(is_DRP)) then
       ! 11-points optimized scheme o(4)
       ! [Bogey & Bailly JCP 2004]
       a11( 1) =  0.872756993962_wp
       a11( 2) = -0.286511173973_wp
       a11( 3) =  0.090320001280_wp
       a11( 4) = -0.020779405824_wp
       a11( 5) =  0.002484594688_wp
       a11( 0) =  0.0_wp
       a11(-1) = -a11(1)
       a11(-2) = -a11(2)
       a11(-3) = -a11(3)
       a11(-4) = -a11(4)
       a11(-5) = -a11(5)
       ! 11-points optimized scheme decentered 4-6
       ! [Berland, Marsden, Bogey & Bailly JCP 2007]
       a46( 1)= 0.016756572303_wp
       a46( 2)=-0.117478455239_wp
       a46( 3)= 0.411034935097_wp
       a46( 4)=-1.130286765151_wp
       a46( 5)= 0.341435872100_wp
       a46( 6)= 0.556396830543_wp
       a46( 7)=-0.082525734207_wp
       a46( 8)= 0.003565834658_wp
       a46( 9)= 0.001173034777_wp
       a46(10)=-0.000071772607_wp
       a46(11)=-0.000000352272_wp
       ! 11-points optimized scheme decentered 3-7
       ! [Berland, Marsden, Bogey & Bailly JCP 2007]
       a37( 1)=-0.013277273810_wp
       a37( 2)= 0.115976072920_wp
       a37( 3)=-0.617479187931_wp
       a37( 4)=-0.274113948206_wp
       a37( 5)= 1.086208764655_wp
       a37( 6)=-0.402951626982_wp
       a37( 7)= 0.131066986242_wp
       a37( 8)=-0.028154858354_wp
       a37( 9)= 0.002596328316_wp
       a37(10)= 0.000128743150_wp
       a37(11)= 0.0_wp
       ! 11-points optimized scheme decentered 2-8
       ! [Berland, Marsden, Bogey & Bailly JCP 2007]
       a28( 1)= 0.046246319744_wp
       a28( 2)=-0.462989982072_wp
       a28( 3)=-0.459203180244_wp
       a28( 4)= 1.205900619436_wp
       a28( 5)=-0.423956587692_wp
       a28( 6)= 0.102329382027_wp
       a28( 7)=-0.006253229685_wp
       a28( 8)=-0.002025942780_wp
       a28( 9)=-0.000016793609_wp
       a28(10)=-0.000015302561_wp
       a28(11)=-0.000015302561_wp
       ! 11-points optimized scheme decentered 1-9
       ! [Berland, Marsden, Bogey & Bailly JCP 2007]
       a19( 1)=-0.180022054228_wp
       a19( 2)=-1.237550583044_wp
       a19( 3)= 2.484731692990_wp
       a19( 4)=-1.810320814061_wp
       a19( 5)= 1.112990048440_wp
       a19( 6)=-0.481086916514_wp
       a19( 7)= 0.126598690230_wp
       a19( 8)=-0.015510730165_wp
       a19( 9)= 0.000021609059_wp
       a19(10)= 0.000156447571_wp
       a19(11)=-0.000007390277_wp
       ! 11-points optimized scheme decentered 0-10
       ! [Berland, Marsden, Bogey & Bailly JCP 2007]
       a010( 1)=-2.391602219538_wp
       a010( 2)= 5.832490322294_wp
       a010( 3)=-7.650218001182_wp
       a010( 4)= 7.907810563576_wp
       a010( 5)=-5.922599052629_wp
       a010( 6)= 3.071037015445_wp
       a010( 7)=-1.014956769726_wp
       a010( 8)= 0.170022256519_wp
       a010( 9)= 0.002819958377_wp
       a010(10)=-0.004791009708_wp
       a010(11)=-0.000013063429_wp
    end if
    
    if (iFD.eq.1) then
       ! DF 3 points
       a11(1) =  0.5_wp
       a11(2) =  0.0_wp
       a11(3) =  0.0_wp
       a11(4) =  0.0_wp
       a11(5) =  0.0_wp
       a11(0) =  0.0_wp
       a11(-1) = -a11(1)
       a11(-2) = -a11(2)
       a11(-3) = -a11(3)
       a11(-4) = -a11(4)
       a11(-5) = -a11(5)
       a9(1) =  0.5_wp
       a9(2) = -0.0_wp
       a9(3) =  0.0_wp
       a9(4) = -0.0_wp
       a9(0) =  0.0_wp
       a9(-1) = -a9(1)
       a9(-2) = -a9(2)
       a9(-3) = -a9(3)
       a9(-4) = -a9(4)
       a7( 1) = 0.5_wp
       a7( 2) = 0.0_wp
       a7( 3) = 0.0_wp
       a7( 0) = 0.0_wp
       a7(-1) = -a7(1)
       a7(-2) = -a7(2)
       a7(-3) = -a7(3)
       a5(1) =  0.5
       a5(2) =  0.
       a5(0) =  0.
       a5(-1) = -a5(1)
       a5(-2) = -a5(2)
    else if (iFD.eq.2) then
       ! DF 5 points
       a11(1) =  2.0_wp/3.0_wp
       a11(2) = -1.0_wp/12.0_wp
       a11(3) =  0.0_wp
       a11(4) =  0.0_wp
       a11(5) =  0.0_wp
       a11(0) =  0.0_wp
       a11(-1) = -a11(1)
       a11(-2) = -a11(2)
       a11(-3) = -a11(3)
       a11(-4) = -a11(4)
       a11(-5) = -a11(5)
       a9(1) =  2.0_wp/3.0_wp
       a9(2) = -1.0_wp/12.0_wp
       a9(3) =  0.0_wp
       a9(4) =  0.0_wp
       a9(0) =  0.0_wp
       a9(-1) = -a9(1)
       a9(-2) = -a9(2)
       a9(-3) = -a9(3)
       a9(-4) = -a9(4)
       a7( 1) =  2.0_wp/3.0_wp
       a7( 2) = -1.0_wp/12.0_wp
       a7( 3) = 0.0_wp
       a7( 0) = 0.0_wp
       a7(-1) = -a7(1)
       a7(-2) = -a7(2)
       a7(-3) = -a7(3)
    else if (iFD.eq.3) then
       ! DF 7 points
       a11(1) =  3.0_wp/4.0_wp
       a11(2) = -3.0_wp/20.0_wp
       a11(3) =  1.0_wp/60.0_wp
       a11(4) =  0.
       a11(5) =  0.
       a11(0) =  0.
       a11(-1) = -a11(1)
       a11(-2) = -a11(2)
       a11(-3) = -a11(3)
       a11(-4) = -a11(4)
       a11(-5) = -a11(5)
       a9(1) =  3.0_wp/4.0_wp
       a9(2) = -3.0_wp/20.0_wp
       a9(3) =  1.0_wp/60.0_wp
       a9(4) =  0.0_wp
       a9(0) =  0.0_wp
       a9(-1) = -a9(1)
       a9(-2) = -a9(2)
       a9(-3) = -a9(3)
       a9(-4) = -a9(4)
    else if (iFD.eq.4) then
       ! DF 9 points
       a11(1) =  4.0_wp/5.0_wp
       a11(2) = -1.0_wp/5.0_wp
       a11(3) =  4.0_wp/105.0_wp
       a11(4) = -1.0_wp/280.0_wp
       a11(5) =  0.0_wp
       a11(0) =  0.0_wp
       a11(-1) = -a11(1)
       a11(-2) = -a11(2)
       a11(-3) = -a11(3)
       a11(-4) = -a11(4)
       a11(-5) = -a11(5)
    endif
  end subroutine init_coeff_deriv

  !===============================================================
  subroutine init_coeff_SBP(order_sbp)
  !===============================================================
    !> Initialize derivative scheme coefficients
    !> SBP (Summation By Parts) schemes
  !===============================================================
    use warnstop
    implicit none
    ! ------------------------------------------------------------
    ! Scheme parameters
    integer, intent(in) :: order_sbp
    ! ------------------------------------------------------------
    real(wp) :: h11,h22,h33,h44,h55,h66
    real(wp) :: q11,q12,q13,q14,q15,q16,q23,q24,q25,q26
    real(wp) :: q34,q35,q36,q45,q46,q47,q56,q57,q58
    ! ------------------------------------------------------------

    ! order 4 (from Gustafsson book)
    ! =======

    as4p0(1)=-24.0_wp/17.0_wp
    as4p0(2)= 59.0_wp/34.0_wp
    as4p0(3)=-4.0_wp/17.0_wp
    as4p0(4)=-3.0_wp/34.0_wp

    as4p1(1)=-1.0_wp/2.0_wp
    as4p1(2)= 0.0_wp
    as4p1(3)= 1.0_wp/2.0_wp
    as4p1(4)= 0.0_wp

    as4p2(1)= 4.0_wp/43.0_wp
    as4p2(2)=-59.0_wp/86.0_wp
    as4p2(3)= 0.0_wp
    as4p2(4)= 59.0_wp/86.0_wp
    as4p2(5)=-4.0_wp/43.0_wp

    as4p3(1)= 3.0_wp/98.0_wp
    as4p3(2)= 0.0_wp
    as4p3(3)=-59.0_wp/98.0_wp
    as4p3(4)= 0.0_wp
    as4p3(5)= 32.0_wp/49.0_wp
    as4p3(6)=-4.0_wp/49.0_wp

    as4m0=-as4p0
    as4m1=-as4p1
    as4m2=-as4p2
    as4m3=-as4p3

    ! Standard schemes order 4 (centered)
    as4( 2) = -1.0_wp/12.0_wp
    as4( 1) =  2.0_wp/3.0_wp
    as4( 0) =  0.0_wp
    as4(-1) = -as4(1)
    as4(-2) = -as4(2)
    
!!$    ! order 4 (from Fernandez et al. C&F 2014)
!!$    ! =======
!!$    
!!$    h11 = 17.0_wp/48.0_wp
!!$    h22 = 59.0_wp/48.0_wp
!!$    h33 = 43.0_wp/48.0_wp
!!$    h44 = 49.0_wp/48.0_wp
!!$    
!!$    q11=-1.0_wp/2.0_wp
!!$    q12= 59.0_wp/96.0_wp
!!$    q13=-1.0_wp/12.0_wp
!!$    q14=-1.0_wp/32.0_wp
!!$    q23= 59.0_wp/96.0_wp
!!$    q24= 0.0_wp
!!$    q34= 59.0_wp/96.0_wp
!!$    
!!$!!    ! setting up H
!!$!!    H = eye(n,n);
!!$!!    
!!$!!    H(1,1)= h11
!!$!!    H(2,2)= h22
!!$!!    H(3,3)= h33
!!$!!    H(4,4)= h44
!!$!!    
!!$!!    H(n  ,n  )= h11
!!$!!    H(n-1,n-1)= h22
!!$!!    H(n-2,n-2)= h33
!!$!!    H(n-3,n-3)= h44
!!$
!!$    !setting up coefficients
!!$    as4p0(1)= q11
!!$    as4p0(2)= q12
!!$    as4p0(3)= q13
!!$    as4p0(4)= q14
!!$    as4p0=as4p0/h11
!!$    
!!$    as4p1(1)=-q12
!!$    as4p1(2)= 0.0_wp
!!$    as4p1(3)= q23
!!$    as4p1(4)= q24
!!$    as4p1=as4p1/h22
!!$    
!!$    as4p2(1)=-q13
!!$    as4p2(2)=-q23
!!$    as4p2(3)= 0.0_wp
!!$    as4p2(4)= q34
!!$    as4p2(5)=-1.0_wp/12.0_wp
!!$    as4p2=as4p2/h33
!!$    
!!$    as4p3(1)=-q14
!!$    as4p3(2)=-q24
!!$    as4p3(3)=-q34
!!$    as4p3(4)= 0.0_wp
!!$    as4p3(5)= 2.0_wp/3.0_wp
!!$    as4p3(6)=-1.0_wp/12.0_wp
!!$    as4p3=as4p3/h44
!!$    
!!$    as4m0=-as4p0
!!$    as4m1=-as4p1
!!$    as4m2=-as4p2
!!$    as4m3=-as4p3
    
    ! order 6 (from Gustafsson book)
    ! =======
    
    as6p0(1)=-21600.0_wp/13649.0_wp
    as6p0(2)= 104009.0_wp/54596.0_wp
    as6p0(3)= 30443.0_wp/81894.0_wp
    as6p0(4)=-33311.0_wp/27298.0_wp
    as6p0(5)= 16863.0_wp/27298.0_wp
    as6p0(6)=-15025.0_wp/163788.0_wp

    as6p1(1)=-104009.0_wp/240260.0_wp
    as6p1(2)= 0.0_wp
    as6p1(3)=-311.0_wp/72078.0_wp
    as6p1(4)= 20229.0_wp/24026.0_wp
    as6p1(5)=-24337.0_wp/48052.0_wp
    as6p1(6)= 36661.0_wp/360390.0_wp

    as6p2(1)=-30443.0_wp/162660.0_wp
    as6p2(2)= 311.0_wp/32532.0_wp
    as6p2(3)= 0.0_wp
    as6p2(4)=-11155.0_wp/16266.0_wp
    as6p2(5)= 41287.0_wp/32532.0_wp
    as6p2(6)=-21999.0_wp/54220.0_wp

    as6p3(1)= 33311.0_wp/107180.0_wp
    as6p3(2)=-20229.0_wp/21436.0_wp
    as6p3(3)= 485.0_wp/1398.0_wp
    as6p3(4)= 0.0_wp
    as6p3(5)= 4147.0_wp/21436.0_wp
    as6p3(6)= 25427.0_wp/321540.0_wp
    as6p3(7)= 72.0_wp/5359.0_wp

    as6p4(1)=-16863.0_wp/78770.0_wp
    as6p4(2)= 24337.0_wp/31508.0_wp
    as6p4(3)=-41287.0_wp/47262.0_wp
    as6p4(4)=-4147.0_wp/15754.0_wp
    as6p4(5)= 0.0_wp
    as6p4(6)= 342523.0_wp/472620.0_wp
    as6p4(7)=-1296.0_wp/7877.0_wp
    as6p4(8)= 144.0_wp/7877.0_wp

    as6p5(1)= 15025.0_wp/525612.0_wp
    as6p5(2)=-36661.0_wp/262806.0_wp
    as6p5(3)= 21999.0_wp/87602.0_wp
    as6p5(4)=-25427.0_wp/262806.0_wp
    as6p5(5)=-342523.0_wp/525612.0_wp
    as6p5(6)= 0.0_wp
    as6p5(7)= 32400.0_wp/43801.0_wp
    as6p5(8)=-6480.0_wp/43801.0_wp
    as6p5(9)= 720.0_wp/43801.0_wp

    as6m0=-as6p0
    as6m1=-as6p1
    as6m2=-as6p2
    as6m3=-as6p3
    as6m4=-as6p4
    as6m5=-as6p5

    ! Standard schemes order 6 (centered)
    as6( 3) =  1.0_wp/60.0_wp
    as6( 2) = -3.0_wp/20.0_wp
    as6( 1) =  3.0_wp/4.0_wp
    as6( 0) =  0.0_wp
    as6(-1) = -as6(1)
    as6(-2) = -as6(2)
    as6(-3) = -as6(3)

    ! order 6 (from Fernandez et al. C&F 2014)
    ! =======
    
    ! optimized value comment out to use other values
    q56 =5591070156686698065364559.0_wp/7931626489314500743872000.0_wp
    
    h11=0.13649e5_wp/0.43200e5_wp
    h22=0.12013e5_wp/0.8640e4_wp
    h33= 0.2711e4_wp/0.4320e4_wp
    h44= 0.5359e4_wp/0.4320e4_wp
    h55= 0.7877e4_wp/0.8640e4_wp
    h66=0.43801e5_wp/0.43200e5_wp

    q11=-1.0_wp/2.0_wp
    q12=-0.953e3_wp/0.16200e5_wp + q56;
    q13= 0.715489e6_wp/0.259200e6_wp - 4.0_wp*q56
    q14=-0.62639e5_wp/0.14400e5_wp + 6.0_wp*q56
    q15= 0.147127e6_wp/0.51840e5_wp - 4.0_wp*q56
    q16=-0.89387e5_wp/0.129600e6_wp + q56
    
    q23=-0.57139e5_wp/0.8640e4_wp + 10.0_wp*q56
    q24= 0.745733e6_wp/0.51840e5_wp - 20.0_wp*q56
    q25=-0.18343e5_wp/0.1728e4_wp + 15.0_wp*q56
    q26= 0.240569e6_wp/0.86400e5_wp - 4.0_wp*q56
    
    q34=-0.176839e6_wp/0.12960e5_wp + 20.0_wp*q56
    q35= 0.242111e6_wp/0.17280e5_wp - 20.0_wp*q56
    q36=-0.182261e6_wp/0.43200e5_wp + 6.0_wp*q56
    
    q45=-0.165041e6_wp/0.25920e5_wp + 10.0_wp*q56
    q46= 0.710473e6_wp/0.259200e6_wp - 4.0_wp*q56
    q47= 1.0_wp/6.0e1_wp
    
    q57=-3.0_wp/2.0e1_wp
    q58= 1.0_wp/6.0e1_wp
    
!!    ! setting up H
!!    H = eye(n,n);
!!    
!!    H(1,1)= h11
!!    H(2,2)= h22
!!    H(3,3)= h33
!!    H(4,4)= h44
!!    H(5,5)= h55
!!    H(6,6)= h66
!!    
!!    H(n  ,n  )= h11
!!    H(n-1,n-1)= h22
!!    H(n-2,n-2)= h33
!!    H(n-3,n-3)= h44
!!    H(n-4,n-4)= h55
!!    H(n-5,n-5)= h66
    
    ! setting up coefficients
    as6p0(1)= q11
    as6p0(2)= q12
    as6p0(3)= q13
    as6p0(4)= q14
    as6p0(5)= q15
    as6p0(6)= q16
    as6p0=as6p0/h11
    
    as6p1(1)=-q12
    as6p1(2)= 0.0_wp
    as6p1(3)= q23
    as6p1(4)= q24
    as6p1(5)= q25
    as6p1(6)= q26
    as6p1=as6p1/h22
    
    as6p2(1)=-q13
    as6p2(2)=-q23
    as6p2(3)= 0.0_wp
    as6p2(4)= q34
    as6p2(5)= q35
    as6p2(6)= q36
    as6p2=as6p2/h33
    
    as6p3(1)=-q14
    as6p3(2)=-q24
    as6p3(3)=-q34
    as6p3(4)= 0.0_wp
    as6p3(5)= q45
    as6p3(6)= q46
    as6p3(7)= q47
    as6p3=as6p3/h44
   
    as6p4(1)=-q15
    as6p4(2)=-q25
    as6p4(3)=-q35
    as6p4(4)=-q45
    as6p4(5)= 0.0_wp
    as6p4(6)= q56
    as6p4(7)= q57
    as6p4(8)= q58
    as6p4=as6p4/h55

    as6p5(1)=-q16
    as6p5(2)=-q26
    as6p5(3)=-q36
    as6p5(4)=-q46
    as6p5(5)=-q56    
    as6p5(6)= 0.0_wp 
    as6p5(7)= 3.0_wp/4.0_wp
    as6p5(8)=-3.0_wp/20.0_wp;
    as6p5(9)= 1.0_wp/60.0_wp;
    as6p5=as6p5/h66

    as6m0=-as6p0
    as6m1=-as6p1
    as6m2=-as6p2
    as6m3=-as6p3
    as6m4=-as6p4
    as6m5=-as6p5
        
  end subroutine init_coeff_SBP

  !===============================================================
  subroutine test_coeff_SBP(order_sbp)
  !===============================================================
    !> Test derivative scheme coefficients
    !> SBP (Summation By Parts) schemes
  !===============================================================
    use warnstop
    implicit none
    ! ------------------------------------------------------------
    ! Scheme parameters
    integer, intent(in) :: order_sbp
    ! ------------------------------------------------------------
    integer, parameter :: n=20
    integer :: i
    real(wp) :: tpi,dxl
    real(wp), dimension(n) :: xl,vl,dvl
    ! ------------------------------------------------------------

    tpi=2.0_wp*acos(-1.0_wp)

    ! grid
    ! ====
    dxl=1.0_wp/dble(n-1)
    xl(1)=0.0_wp
    do i=2,n
       xl(i)=xl(i-1)+dxl
    enddo
    
    ! test function
    ! =============
    vl=sin(tpi*xl)

    ! write results (single proc)
    ! =============    
    open(194,file='dvl.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) n
    write(194) (xl(i),i=1,n)
    write(194) (vl(i),i=1,n)
    
    ! derivative
    ! ==========
    
    ! order 4
    ! -------
    i=1
    dvl(i)= as4p0(1)*vl(1)+as4p0(2)*vl(2) &
          + as4p0(3)*vl(3)+as4p0(4)*vl(4)
    i=2
    dvl(i)= as4p1(1)*vl(1)+as4p1(2)*vl(2) &
          + as4p1(3)*vl(3)+as4p1(4)*vl(4)
    i=3
    dvl(i)= as4p2(1)*vl(1)+as4p2(2)*vl(2) &
          + as4p2(3)*vl(3)+as4p2(4)*vl(4) &
          + as4p2(5)*vl(5)
    i=4
    dvl(i)= as4p3(1)*vl(1)+as4p3(2)*vl(2) &
          + as4p3(3)*vl(3)+as4p3(4)*vl(4) &
          + as4p3(5)*vl(5)+as4p3(6)*vl(6)
    
    do i=5,n-4
       dvl(i)= as4(1)* (vl(i+1)-vl(i-1)) &
             + as4(2)* (vl(i+2)-vl(i-2))
    enddo
    
    i=n-3
    dvl(i)= as4m3(1)*vl(n  )+as4m3(2)*vl(n-1) &
          + as4m3(3)*vl(n-2)+as4m3(4)*vl(n-3) &
          + as4m3(5)*vl(n-4)+as4m3(6)*vl(n-5)
    i=n-2
    dvl(i)= as4m2(1)*vl(n  )+as4m2(2)*vl(n-1) &
          + as4m2(3)*vl(n-2)+as4m2(4)*vl(n-3) &
          + as4m2(5)*vl(n-4)
    i=n-1
    dvl(i)= as4m1(1)*vl(n  )+as4m1(2)*vl(n-1) &
          + as4m1(3)*vl(n-2)+as4m1(4)*vl(n-3)
    i=n
    dvl(i)= as4m0(1)*vl(n  )+as4m0(2)*vl(n-1) &
          + as4m0(3)*vl(n-2)+as4m0(4)*vl(n-3)
       
    dvl=dvl/dxl/tpi
    
    ! write results (single proc)
    write(194) (dvl(i),i=1,n)
    
    ! order 6
    ! -------
    i=1
    dvl(i)= as6p0(1)*vl(1)+as6p0(2)*vl(2) &
          + as6p0(3)*vl(3)+as6p0(4)*vl(4) &
          + as6p0(5)*vl(5)+as6p0(6)*vl(6)
    i=2
    dvl(i)= as6p1(1)*vl(1)+as6p1(2)*vl(2) &
          + as6p1(3)*vl(3)+as6p1(4)*vl(4) &
          + as6p1(5)*vl(5)+as6p1(6)*vl(6)
    i=3
    dvl(i)= as6p2(1)*vl(1)+as6p2(2)*vl(2) &
          + as6p2(3)*vl(3)+as6p2(4)*vl(4) &
          + as6p2(5)*vl(5)+as6p2(6)*vl(6)
    i=4
    dvl(i)= as6p3(1)*vl(1)+as6p3(2)*vl(2) &
          + as6p3(3)*vl(3)+as6p3(4)*vl(4) &
          + as6p3(5)*vl(5)+as6p3(6)*vl(6) &
          + as6p3(7)*vl(7)
    i=5
    dvl(i)= as6p4(1)*vl(1)+as6p4(2)*vl(2) &
          + as6p4(3)*vl(3)+as6p4(4)*vl(4) &
          + as6p4(5)*vl(5)+as6p4(6)*vl(6) &
          + as6p4(7)*vl(7)+as6p4(8)*vl(8)
    i=6
    dvl(i)= as6p5(1)*vl(1)+as6p5(2)*vl(2) &
          + as6p5(3)*vl(3)+as6p5(4)*vl(4) &
          + as6p5(5)*vl(5)+as6p5(6)*vl(6) &
          + as6p5(7)*vl(7)+as6p5(8)*vl(8) &
          + as6p5(9)*vl(9)
    
    do i=7,n-6
       dvl(i)= as6(1)* (vl(i+1)-vl(i-1)) &
             + as6(2)* (vl(i+2)-vl(i-2)) &
             + as6(3)* (vl(i+3)-vl(i-3))
    enddo
    
    i=n-5
    dvl(i)= as6m5(1)*vl(n  )+as6m5(2)*vl(n-1) &
          + as6m5(3)*vl(n-2)+as6m5(4)*vl(n-3) &
          + as6m5(5)*vl(n-4)+as6m5(6)*vl(n-5) &
          + as6m5(7)*vl(n-6)+as6m5(8)*vl(n-7) &
          + as6m5(9)*vl(n-8)
    i=n-4
    dvl(i)= as6m4(1)*vl(n  )+as6m4(2)*vl(n-1) &
          + as6m4(3)*vl(n-2)+as6m4(4)*vl(n-3) &
          + as6m4(5)*vl(n-4)+as6m4(6)*vl(n-5) &
          + as6m4(7)*vl(n-6)+as6m4(8)*vl(n-7)
    i=n-3
    dvl(i)= as6m3(1)*vl(n  )+as6m3(2)*vl(n-1) &
          + as6m3(3)*vl(n-2)+as6m3(4)*vl(n-3) &
          + as6m3(5)*vl(n-4)+as6m3(6)*vl(n-5) &
          + as6m3(7)*vl(n-6)
    i=n-2
    dvl(i)= as6m2(1)*vl(n  )+as6m2(2)*vl(n-1) &
          + as6m2(3)*vl(n-2)+as6m2(4)*vl(n-3) &
          + as6m2(5)*vl(n-4)+as6m2(6)*vl(n-5)
    i=n-1
    dvl(i)= as6m1(1)*vl(n  )+as6m1(2)*vl(n-1) &
          + as6m1(3)*vl(n-2)+as6m1(4)*vl(n-3) &
          + as6m1(5)*vl(n-4)+as6m1(6)*vl(n-5)
    i=n
    dvl(i)= as6m0(1)*vl(n  )+as6m0(2)*vl(n-1) &
          + as6m0(3)*vl(n-2)+as6m0(4)*vl(n-3) &
          + as6m0(5)*vl(n-4)+as6m0(6)*vl(n-5)
       
    dvl=dvl/dxl/tpi
    
    ! write results (single proc)
    write(194) (dvl(i),i=1,n)
    close(194)
    
  end subroutine test_coeff_SBP
  

  !===============================================================
  subroutine init_coeff_wall_cart
  !===============================================================
    !> Initialize extrapolation coefficients for bc_wall
    !> for cartesian solver
  !===============================================================
    use mod_bc
    use mod_grid
    use mod_constant
    use mod_mpi
    implicit none
    ! ------------------------------------------------------------
    real(wp) :: ry
    ! ------------------------------------------------------------

    ! imin
    if (is_bc_wall(1,1)) then
      ry = (x(3)-x(2))/(x(2)-x(1))
      cextrp2(1,1) = (1.0_wp+ry)**2/(ry*(2.0_wp+ry))
      cextrp3(1,1) =  1.0_wp       /(ry*(2.0_wp+ry))
    endif
    ! imax
    if (is_bc_wall(1,2)) then
      ry = (x(nx-2)-x(nx-1))/(x(nx-1)-x(nx))
      cextrp2(1,2) = (1.0_wp+ry)**2/(ry*(2.0_wp+ry))
      cextrp3(1,2) =  1.0_wp       /(ry*(2.0_wp+ry))
    endif
    ! jmin
    if (is_bc_wall(2,1)) then
      ry = (y(3)-y(2))/(y(2)-y(1))
      cextrp2(2,1) = (1.0_wp+ry)**2/(ry*(2.0_wp+ry))
      cextrp3(2,1) =  1.0_wp       /(ry*(2.0_wp+ry))
    endif
    ! jmax
    if (is_bc_wall(2,2)) then
      ry = (y(ny-2)-y(ny-1))/(y(ny-1)-y(ny))
      cextrp2(2,2) = (1.0_wp+ry)**2/(ry*(2.0_wp+ry))
      cextrp3(2,2) =  1.0_wp       /(ry*(2.0_wp+ry))
    endif
    ! if 3D
    if (nz.gt.1) then
      ! kmin
      if (is_bc_wall(3,1)) then
        ry = (z(3)-z(2))/(z(2)-z(1))
        cextrp2(3,1) = (1.0_wp+ry)**2/(ry*(2.0_wp+ry))
        cextrp3(3,1) =  1.0_wp       /(ry*(2.0_wp+ry))
      endif
      ! kmax
      if (is_bc_wall(3,2)) then
        ry = (z(nz-2)-z(nz-1))/(z(nz-1)-z(nz))
        cextrp2(3,2) = (1.0_wp+ry)**2/(ry*(2.0_wp+ry))
        cextrp3(3,2) =  1.0_wp       /(ry*(2.0_wp+ry))
      endif
    endif

    if ((iproc.eq.0).and.(verbose)) then
     write(*,*) 'c2 for border_wall (imin to kmax):', cextrp2
     write(*,*) 'c3 for border_wall (imin to kmax):', cextrp3
    endif

  end subroutine init_coeff_wall_cart

end module mod_coeff_deriv
