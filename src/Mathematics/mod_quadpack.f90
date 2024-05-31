!=================================================================================
module mod_quadpack
!=================================================================================
  !> Module for numerical quadratures using Gauss-Kronrod rule
  !> adapted from QUADPACK routines (extended to double precision)
!=================================================================================
  use precision
  implicit none
  !-------------------------------------------------------------------------------
  ! Gauss-Kronrod abscissas and weight pairs
  real(wp), dimension(:), allocatable :: xgk,wgk,wg
  !-------------------------------------------------------------------------------

contains
  
  !===============================================================================
  subroutine init_GK(n_gk)
  !===============================================================================
    !> Initialization of Gauss-Kronrod abscissae and weight pairs
    !
    !  The abscissae and weights are given for the interval (-1,1).
    !  Because of symmetry only the positive abscissae and their corresponding
    !  weights are given.
    !
    !  xgk: abscissae of the Kronrod rule
    !       -> xgk(2),xgk(4)... abscissae of the Gauss rule
    !       -> xgk(1),xgk(3)... optimally added abscissae to the Gauss rule
    !
    !  wgk: weights of the Kronrod rule
    !
    !  wg : weigths of the Gauss rule
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: n_gk ! # of Gauss-Kronrod points
    ! ----------------------------------------------------------------------------        

    if (n_gk==15) then

       allocate(xgk(8),wgk(8),wg(4))
       ! xgk: abscissae of the 15-point Kronrod rule       
       xgk=[0.9914553711208126_wp,0.9491079123427585_wp,0.8648644233597691_wp,0.7415311855993944_wp,&
            0.5860872354676911_wp,0.4058451513773972_wp,0.2077849550078985_wp,0.0_wp ]
       ! wgk: weights of the 15-point Kronrod rule
       wgk=[2.293532201052922e-2_wp,6.309209262997855e-2_wp,0.1047900103222502_wp,0.1406532597155259_wp,&
              0.1690047266392679_wp, 0.1903505780647854_wp, 0.2044329400752989_wp,0.2094821410847278_wp]
       ! wg: weigths of the 8-point Gauss rule
       wg =[0.1294849661688697_wp,0.2797053914892767_wp,0.3818300505051189_wp,0.4179591836734694_wp]

    elseif (n_gk==21) then

       allocate(xgk(11),wgk(11),wg(5))
       ! xgk: abscissae of the 21-point Kronrod rule       
       xgk=[0.9956571630258081_wp,0.9739065285171717_wp,0.9301574913557082_wp,0.8650633666889845_wp,&
            0.7808177265864169_wp,0.6794095682990244_wp,0.5627571346686047_wp,0.4333953941292472_wp,&
            0.2943928627014602_wp,0.1488743389816312_wp,0.0_wp]
       ! wgk: weights of the 21-point Kronrod rule
       wgk=[1.169463886737187e-2_wp,3.255816230796473e-2_wp,5.475589657435200e-2_wp,7.503967481091995e-2_wp,&
            9.312545458369761e-2_wp,1.093871588022976e-1_wp,1.234919762620659e-1_wp,1.347092173114733e-1_wp,&
            1.427759385770601e-1_wp,1.477391049013385e-1_wp,1.494455540029169e-1_wp]
       ! wg: weigths of the 10-point Gauss rule
       wg =[0.06667134430868814_wp,0.1494513491505806_wp,0.2190863625159820_wp,0.2692667193099964_wp,&
            0.2955242247147529_wp]
       
    elseif (n_gk==31) then

       allocate(xgk(16),wgk(16),wg(8))
       ! xgk: abscissae of the 31-point Kronrod rule       
       xgk=[0.9980022986933971_wp,0.9879925180204854_wp,0.9677390756791391_wp,0.9372733924007059_wp,&
            0.8972645323440819_wp,0.8482065834104272_wp,0.7904185014424659_wp,0.7244177313601700_wp,&
            0.6509967412974170_wp,0.5709721726085388_wp,0.4850818636402397_wp,0.3941513470775634_wp,&
            0.2991800071531688_wp,0.2011940939974345_wp,0.1011420669187175_wp,0.0_wp]
       ! wgk: weights of the 31-point Kronrod rule
       wgk=[5.377479872923349e-3_wp,1.500794732931612e-2_wp,2.546084732671532e-2_wp,3.534636079137585e-2_wp,&
            4.458975132476488e-2_wp,5.348152469092809e-2_wp,6.200956780067064e-2_wp,6.985412131872826e-2_wp,&
            7.684968075772038e-2_wp,8.308050282313302e-2_wp,8.856444305621177e-2_wp,9.312659817082532e-2_wp,&
            9.664272698362368e-2_wp,9.917359872179196e-2_wp,1.007698455238756e-1_wp,1.013300070147915e-1_wp]
       ! wg: weigths of the 15-point Gauss rule
       wg =[3.075324199611727e-2_wp,7.036604748810812e-2_wp,1.071592204671719e-1_wp,1.395706779261543e-1_wp,&
            1.662692058169939e-1_wp,1.861610000155622e-1_wp,1.984314853271116e-1_wp,2.025782419255613e-1_wp]
       
    elseif (n_gk==41) then
       
       allocate(xgk(21),wgk(21),wg(10))
       ! xgk: abscissae of the 41-point Kronrod rule       
       xgk=[0.9988590315882777_wp,0.9931285991850949_wp,0.9815078774502503_wp,0.9639719272779138_wp,&
            0.9408226338317548_wp,0.9122344282513259_wp,0.8782768112522820_wp,0.8391169718222188_wp,&
            0.7950414288375512_wp,0.7463319064601508_wp,0.6932376563347514_wp,0.6360536807265150_wp,&
            0.5751404468197103_wp,0.5108670019508271_wp,0.4435931752387251_wp,0.3737060887154196_wp,&
            0.3016278681149130_wp,0.2277858511416451_wp,0.1526054652409227_wp,0.07652652113349733_wp,&
            0.0_wp]
       ! wgk: weights of the 41-point Kronrod rule
       wgk=[3.073583718520532e-3_wp,8.600269855642942e-3_wp,1.462616925697125e-2_wp,2.038837346126652e-2_wp,&
            2.588213360495116e-2_wp,3.128730677703280e-2_wp,3.660016975820080e-2_wp,4.166887332797369e-2_wp,&
            4.643482186749767e-2_wp,5.094457392372869e-2_wp,5.519510534828599e-2_wp,5.911140088063957e-2_wp,&
            6.265323755478117e-2_wp,6.583459713361842e-2_wp,6.864867292852162e-2_wp,7.105442355344407e-2_wp,&
            7.303069033278667e-2_wp,7.458287540049919e-2_wp,7.570449768455667e-2_wp,7.637786767208074e-2_wp,&
            7.660071191799966e-2_wp]
       ! wg: weigths of the 20-point Gauss rule
       wg =[1.761400713915212e-2_wp,4.060142980038694e-2_wp,6.267204833410906e-2_wp,8.327674157670475e-2_wp,&
            1.019301198172404e-1_wp,1.181945319615184e-1_wp,1.316886384491766e-1_wp,1.420961093183821e-1_wp,&
            1.491729864726037e-1_wp,1.527533871307259e-1_wp]

    elseif (n_gk==51) then
       
       allocate(xgk(26),wgk(26),wg(13))
       ! xgk: abscissae of the 51-point Kronrod rule
       xgk=[0.9992621049926098_wp,0.9955569697904981_wp,0.9880357945340772_wp,0.9766639214595175_wp,&
            0.9616149864258425_wp,0.9429745712289743_wp,0.9207471152817016_wp,0.8949919978782754_wp,&
            0.8658470652932756_wp,0.8334426287608340_wp,0.7978737979985001_wp,0.7592592630373576_wp,&
            0.7177664068130844_wp,0.6735663684734684_wp,0.6268100990103174_wp,0.5776629302412230_wp,&
            0.5263252843347192_wp,0.4730027314457150_wp,0.4178853821930377_wp,0.3611723058093878_wp,&
            0.3030895389311078_wp,0.2438668837209884_wp,0.1837189394210489_wp,0.1228646926107104_wp,&
            0.06154448300568508_wp,0.0_wp]
       ! wgk: weights of the 51-point Kronrod rule
       wgk=[1.987383892330316e-3_wp,5.561932135356714e-3_wp,9.473973386174152e-3_wp,1.323622919557167e-2_wp,&
            1.684781770912830e-2_wp,2.043537114588284e-2_wp,2.400994560695322e-2_wp,2.747531758785174e-2_wp,&
            3.079230016738749e-2_wp,3.400213027432934e-2_wp,3.711627148341554e-2_wp,4.008382550403238e-2_wp,&
            4.287284502017005e-2_wp,4.550291304992179e-2_wp,4.798253713883671e-2_wp,5.027767908071567e-2_wp,&
            5.236288580640748e-2_wp,5.425112988854549e-2_wp,5.595081122041232e-2_wp,5.743711636156783e-2_wp,&
            5.868968002239421e-2_wp,5.972034032417406e-2_wp,6.053945537604586e-2_wp,6.112850971705305e-2_wp,&
            6.147118987142532e-2_wp,6.158081806783294e-2_wp]
       ! wg: weigths of the 25-point Gauss rule
       wg =[1.139379850102629e-2_wp,2.635498661503214e-2_wp,4.093915670130631e-2_wp,5.490469597583519e-2_wp,&
            6.803833381235692e-2_wp,8.014070033500102e-2_wp,9.102826198296365e-2_wp,1.005359490670506e-1_wp,&
            1.085196244742637e-1_wp,1.148582591457116e-1_wp,1.194557635357848e-1_wp,1.222424429903100e-1_wp,&
            1.231760537267155e-1_wp]

    elseif (n_gk==61) then
       
       allocate(xgk(31),wgk(31),wg(15))
       ! xgk: abscissae of the 61-point Kronrod rule
       xgk=[0.9994844100504906_wp,0.9968934840746495_wp,0.9916309968704046_wp,0.9836681232797472_wp,&
            0.9731163225011263_wp,0.9600218649683075_wp,0.9443744447485600_wp,0.9262000474292743_wp,&
            0.9055733076999078_wp,0.8825605357920527_wp,0.8572052335460611_wp,0.8295657623827684_wp,&
            0.7997278358218391_wp,0.7677774321048262_wp,0.7337900624532268_wp,0.6978504947933158_wp,&
            0.6600610641266270_wp,0.6205261829892429_wp,0.5793452358263617_wp,0.5366241481420199_wp,&
            0.4924804678617786_wp,0.4470337695380892_wp,0.4004012548303944_wp,0.3527047255308781_wp,&
            0.3040732022736251_wp,0.2546369261678898_wp,0.2045251166823099_wp,0.1538699136085835_wp,&
            0.1028069379667370_wp,0.05147184255531770_wp,0.0_wp]
       ! wgk: weights of the 61-point Kronrod rule
       wgk=[1.389013698677008e-3_wp,3.890461127099884e-3_wp,6.630703915931292e-3_wp,9.273279659517763e-3_wp,&
            1.182301525349634e-2_wp,1.436972950704580e-2_wp,1.692088918905327e-2_wp,1.941414119394238e-2_wp,&
            2.182803582160919e-2_wp,2.419116207808060e-2_wp,2.650995488233310e-2_wp,2.875404876504129e-2_wp,&
            3.090725756238776e-2_wp,3.298144705748373e-2_wp,3.497933802806002e-2_wp,3.688236465182123e-2_wp,&
            3.867894562472759e-2_wp,4.037453895153596e-2_wp,4.196981021516425e-2_wp,4.345253970135607e-2_wp,&
            4.481480013316266e-2_wp,4.605923827100699e-2_wp,4.718554656929915e-2_wp,4.818586175708713e-2_wp,&
            4.905543455502978e-2_wp,4.979568342707421e-2_wp,5.040592140278235e-2_wp,5.088179589874961e-2_wp,&
            5.122154784925877e-2_wp,5.142612853745903e-2_wp,5.149472942945157e-2_wp]
       ! wg: weigths of the 30-point Gauss rule
       wg =[7.968192496166606e-3_wp,1.846646831109096e-2_wp,2.878470788332337e-2_wp,3.879919256962705e-2_wp,&
            4.840267283059405e-2_wp,5.749315621761907e-2_wp,6.597422988218050e-2_wp,7.375597473770521e-2_wp,&
            8.075589522942022e-2_wp,8.689978720108298e-2_wp,9.212252223778613e-2_wp,9.636873717464426e-2_wp,&
            9.959342058679527e-2_wp,1.017623897484055e-1_wp,1.028526528935588e-1_wp]
    endif

  end subroutine init_GK

  !===============================================================================
  subroutine qag(f,a,b,epsabs,epsrel,key,resu,abserr,neval,ier)
  !===============================================================================
    !> QAG approximates an integral over a finite interval [QUADPACK]
    !
    ! The routine calculates an approximation "resu" to a definite integral   
    !      I=integral of f over (a,b),
    !  hopefully satisfying
    !    ||I-resu|| <= max(epsabsS,epsrel*||I||)
    !
    ! QAG is a simple globally adaptive integrator using the strategy of 
    ! Aind (Piessens,1973). It is possible to choose between 6 pairs of
    ! Gauss-Kronrod quadrature formulae for the rule evaluation component. 
    ! The pairs of high degree of precision are suitable for handling
    ! integration difficulties due to a strongly oscillating integrand.
    !
    !  KEY chooses the order of the local integration rule:
    !  1:  7 Gauss points, 15 Gauss-Kronrod points,
    !  2: 10 Gauss points, 21 Gauss-Kronrod points,
    !  3: 15 Gauss points, 31 Gauss-Kronrod points,
    !  4: 20 Gauss points, 41 Gauss-Kronrod points,
    !  5: 25 Gauss points, 51 Gauss-Kronrod points,
    !  6: 30 Gauss points, 61 Gauss-Kronrod points.
    !
    !  IER return code.
    !  0: normal and reliable termination of the routine.  It is assumed that the 
    !     requested accuracy has been achieved.
    !  1: maximum number of subdivisions allowed has been achieved.  One can 
    !     allow more subdivisions by increasing the value of LIMIT in QAG. 
    !     However,if this yields no improvement it is advised to analyze the
    !     integrand to determine the integration difficulties.  If the position
    !     of a local difficulty can be determined,such as a singularity or
    !     discontinuity within the interval) one will probably gain from 
    !     splitting up the interval at this point and calling the integrator 
    !     on the subranges.  If possible,an appropriate special-purpose 
    !     integrator should be used which is designed for handling the type 
    !     of difficulty involved.
    !  2: the occurrence of roundoff error is detected,which prevents the
    !     requested tolerance from being achieved.
    !  3: extremely bad integrand behavior occurs at some points of the
    !     integration interval.
    !  6: the input is invalid,because epsabs<0 and epsrel<0.
    !
    !  Authors:
    !  Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
    !  QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
    !
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: key      ! order of the local integration rule
    integer, intent(out) :: neval   ! number of times the integral was evaluated
    integer, intent(out) :: ier     ! return code
    real(wp), external :: f         ! name of the integrand function routine
    real(wp), intent(in) :: a,b     ! integration bounds
    real(wp), intent(in) :: epsabs  ! absolute accuracy requested
    real(wp), intent(in) :: epsrel  ! relative accuracy requested
    real(wp), intent(out) :: resu   ! estimated value of the integral
    real(wp), intent(out) :: abserr ! estimate of ||I-resu||
    ! ----------------------------------------------------------------------------           
    integer, parameter :: limit=500 ! maximum number of subintervals allowed
                                    ! in the subdivision process of QAGE
    integer :: last
    integer, dimension(limit) :: iord
    real(wp), dimension(limit) :: alist,blist,elist,rlist
    ! ----------------------------------------------------------------------------           

    ! Gauss-Kronrod abscissas and weight pairs
    ! ========================================
    select case(key)
    case(1)
       ! xgk: abscissae of the 15-point Kronrod rule
       ! wgk: weights of the 15-point Kronrod rule
       !  wg: weights of the 7-point Gauss rule
       if (.not.allocated(xgk)) call init_GK(15)      
    case(2)
       ! xgk: abscissae of the 21-point Kronrod rule
       ! wgk: weights of the 21-point Kronrod rule
       !  wg: weights of the 10-point Gauss rule
       if (.not.allocated(xgk)) call init_GK(21)        
    case(3)
       ! xgk: abscissae of the 31-point Kronrod rule
       ! wgk: weights of the 31-point Kronrod rule
       !  wg: weights of the 15-point Gauss rule
       if (.not.allocated(xgk)) call init_GK(31)      
    case(4)       
       ! xgk: abscissae of the 41-point Kronrod rule
       ! wgk: weights of the 41-point Kronrod rule
       !  wg: weights of the 20-point Gauss rule
       if (.not.allocated(xgk)) call init_GK(41)
    case(5)      
       ! xgk: abscissae of the 51-point Kronrod rule
       ! wgk: weights of the 51-point Kronrod rule
       !  wg: weights of the 25-point Gauss rule
       if (.not.allocated(xgk)) call init_GK(51)
    case(6)      
       ! xgk: abscissae of the 61-point Kronrod rule
       ! wgk: weights of the 61-point Kronrod rule
       !  wg: weights of the 30-point Gauss rule
       if (.not.allocated(xgk)) call init_GK(61)
    end select

    call qage(f,a,b,epsabs,epsrel,key,limit,resu,abserr,neval, &
              ier,alist,blist,rlist,elist,iord,last)

  end subroutine qag

  !===============================================================================
  subroutine qage(f,a,b,epsabs,epsrel,key,limit,resu,abserr,neval, &
                  ier,alist,blist,rlist,elist,iord,last)
  !===============================================================================
    !> QAGE estimates a definite integral [QUADPACK]
    !
    ! The routine calculates an approximation "resu" to a definite integral   
    !      I=integral of f over (a,b),
    !  hopefully satisfying
    !    ||I-resu|| <= max(epsabsS,epsrel*||I||)
    !
    !  Authors:
    !  Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
    !  QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
    !
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: limit                ! maximum number of subintervals allowed
    integer, intent(in) :: key      ! order of the local integration rule
    integer, intent(out) :: neval   ! number of times the integral was evaluated
    integer, intent(out) :: ier     ! return code
    integer :: last ! number of subintervals actually produced in the subdivision
    real(wp), external :: f         ! name of the integrand function routine
    real(wp), intent(in) :: a,b     ! integration bounds
    real(wp), intent(in) :: epsrel  ! relative accuracy requested
    real(wp), intent(in) :: epsabs  ! absolute accuracy requested
    real(wp), intent(out) :: resu   ! estimated value of the integral
    real(wp), intent(out) :: abserr ! estimate of ||I-resu||
    integer, intent(out), dimension(limit) :: iord ! pointers to the error estimates
    real(wp), dimension(limit) :: alist,blist,elist,rlist ! workspace
    ! alist,blist contains in entries 1 through LAST the left and right ends of the partition subintervals
    ! rlist contains in entries 1 through LAST the integral approximations on the subintervals
    ! elist contains in entries 1 through LAST the absolute error estimates on the subintervals
    ! ----------------------------------------------------------------------------           
    !  Local parameters:
    !  -----------------
    !    alist     - list of left end points of all subintervals considered up to now
    !    blist     - list of right end points of all subintervals considered up to now
    !    elist(i)  - error estimate applying to rlist(i)
    !    maxerr    - pointer to the interval with largest error estimate
    !    errmax    - elist(maxerr)
    !    area      - sum of the integrals over the subintervals
    !    errsum    - sum of the errors over the subintervals
    !    errbnd    - requested accuracy max(epsabs,epsrel*abs(resu))
    !    *****1    - variable for the left subinterval
    !    *****2    - variable for the right subinterval
    !    last      - index for subdivision
    integer :: iroff1,iroff2
    integer :: keyf,maxerr,nrmax
    real(wp) :: area,area1,area12,area2
    real(wp) :: a1,a2,b1,b2,c,resabs
    real(wp) :: defabs,defab1,defab2
    real(wp) :: errbnd,errmax,errsum
    real(wp) :: error1,error2,erro12
    ! ----------------------------------------------------------------------------           

    ! Test on validity of parameters
    ! ==============================
    ier=0
    neval=0
    last=0
    resu=0.0_wp
    abserr=0.0_wp
    alist(1)=a
    blist(1)=b
    rlist(1)=0.0_wp
    elist(1)=0.0_wp
    iord(1)=0

    if (epsabs<0.0_wp .and. epsrel<0.0_wp ) then
       ier=6
       return
    endif

    ! First approximation to the integral
    ! ===================================
    keyf=key
    keyf=max(keyf,1)
    keyf=min(keyf,6)

    c=keyf
    neval=0

    if (keyf==1) then
       call qk15(f,a,b,resu,abserr,defabs,resabs)
    else if (keyf==2) then
       call qk21(f,a,b,resu,abserr,defabs,resabs)
    else if (keyf==3) then
       call qk31(f,a,b,resu,abserr,defabs,resabs)
    else if (keyf==4) then
       call qk41(f,a,b,resu,abserr,defabs,resabs)
    else if (keyf==5) then
       call qk51(f,a,b,resu,abserr,defabs,resabs)
    else if (keyf==6) then
       call qk61(f,a,b,resu,abserr,defabs,resabs)
    endif

    last=1
    rlist(1)=resu
    elist(1)=abserr
    iord(1)=1

    ! Test on accuracy
    ! ================

    errbnd=max(epsabs,epsrel*abs(resu))

    if (abserr<=50.0_wp*epsilon(defabs)*defabs .and. &
         errbnd<abserr) then
       ier=2
    endif

    if (limit==1) then
       ier=1
    endif

    if (ier/=0 .or. (abserr<=errbnd.and.abserr/=resabs) &
        .or. abserr==0.0_wp ) then

       if (keyf/=1) then
          neval=(10*keyf+1)*(2*neval+1)
       else
          neval=30*neval+15
       endif

       return

    endif

    ! Initialization
    ! ==============
    errmax=abserr
    maxerr=1
    area=resu
    errsum=abserr
    nrmax=1
    iroff1=0
    iroff2=0

    do last=2,limit
   
       ! Bisect the subinterval with the largest error estimate
       ! ======================================================
   
       a1=alist(maxerr)
       b1=0.5_wp*(alist(maxerr)+blist(maxerr))
       a2=b1
       b2=blist(maxerr)

       if (keyf==1) then
          call qk15(f,a1,b1,area1,error1,resabs,defab1)
       else if (keyf==2) then
          call qk21(f,a1,b1,area1,error1,resabs,defab1)
       else if (keyf==3) then
          call qk31(f,a1,b1,area1,error1,resabs,defab1)
       else if (keyf==4) then
          call qk41(f,a1,b1,area1,error1,resabs,defab1)
       else if (keyf==5) then
          call qk51(f,a1,b1,area1,error1,resabs,defab1)
       else if (keyf==6) then
          call qk61(f,a1,b1,area1,error1,resabs,defab1)
       endif

       if (keyf==1) then
          call qk15(f,a2,b2,area2,error2,resabs,defab2)
       else if (keyf==2) then
          call qk21(f,a2,b2,area2,error2,resabs,defab2)
       else if (keyf==3) then
          call qk31(f,a2,b2,area2,error2,resabs,defab2)
       else if (keyf==4) then
          call qk41(f,a2,b2,area2,error2,resabs,defab2)
       else if (keyf==5) then
          call qk51(f,a2,b2,area2,error2,resabs,defab2)
       else if (keyf==6) then
          call qk61(f,a2,b2,area2,error2,resabs,defab2)
       endif
   
       ! Improve previous approximations to integral and error and test for accuracy
       ! ===========================================================================   
       neval=neval+1
       area12=area1+area2
       erro12=error1+error2
       errsum=errsum+erro12 - errmax
       area=area+area12 - rlist(maxerr)

       if (defab1/=error1 .and. defab2/=error2) then

          if (abs(rlist(maxerr)-area12) <= 1.0e-5_wp*abs(area12) &
               .and. 0.99_wp*errmax<=erro12) then
             iroff1=iroff1+1
          endif

          if (10<last .and. errmax<erro12) then
             iroff2=iroff2+1
          endif

       endif

       rlist(maxerr)=area1
       rlist(last)=area2
       errbnd=max(epsabs,epsrel*abs(area))
   
       ! Test for roundoff error and eventually set error flag
       ! =====================================================
   
       if (errbnd<errsum) then

          if (6<=iroff1 .or. 20<=iroff2) then
             ier=2
          endif
      
          ! Set error flag in the case that the number of subintervals equals limit
          ! =======================================================================
          if (last==limit) then
             ier=1
          endif
      
          ! Set error flag in the case of bad integrand behavior
          ! at a point of the integration range
          ! ====================================================
          if (max(abs(a1),abs(b2))<=(1.0_wp+c*1.0e3_wp*epsilon(a1)) &
               *(abs(a2)+1.0e4_wp*tiny(a2))) then
             ier=3
          endif

       endif
   
       ! Append the newly-created intervals to the list
       ! ==============================================  
       if (error2 <= error1 ) then
          alist(last)=a2
          blist(maxerr)=b1
          blist(last)=b2
          elist(maxerr)=error1
          elist(last)=error2
       else
          alist(maxerr)=a2
          alist(last)=a1
          blist(last)=b1
          rlist(maxerr)=area2
          rlist(last)=area1
          elist(maxerr)=error2
          elist(last)=error1
       endif
   
       ! Call QSORT to maintain the descending ordering
       ! in the list of error estimates and select the subinterval
       ! with the largest error estimate (to be bisected next)
       ! =========================================================
       call qsort(limit,last,maxerr,errmax,elist,iord,nrmax)

       if (ier/=0 .or. errsum<=errbnd) then
          exit
       endif

    enddo

    ! Compute final resu
    ! ==================
    resu=sum(rlist(1:last))

    abserr=errsum

    if (keyf/= 1) then
       neval=(10*keyf+1)*(2*neval+1)
    else
       neval=30*neval+15
    endif

  end subroutine qage

  !===============================================================================
  subroutine  qk15(f,a,b,resu,abserr,resabs,resasc)
  !===============================================================================
    !> QK15 carries out a 15 point Gauss-Kronrod quadrature rule [QUADPACK]
    !
    ! This routine approximates I=integral(a<=x<=b)f(x)dx
    ! with an error estimate, and J=integral(a<=x<=b)|f(x)|dx
    !
    !  Authors:
    !  Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
    !  QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
    !
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), external :: f         ! name of the integrand function routine
    real(wp), intent(in) :: a,b     ! integration bounds
    real(wp), intent(out) :: resu   ! estimated value of the integral
    real(wp), intent(out) :: abserr ! estimate of |I-resu|
    real(wp), intent(out) :: resabs ! approximation to the integral of the absolute value of f
    real(wp), intent(out) :: resasc ! approximation to the integral |f-I/(b-a)| over [a,b]
    ! ----------------------------------------------------------------------------        
    integer :: j,jtw,jtwm1
    real(wp) :: centr ! mid point of the interval
    real(wp) :: hlgth,dhlgth
    real(wp) :: absc  ! abscissa
    real(wp) :: fc,fsum
    real(wp) :: fval1,fval2 ! function values
    real(wp) :: fv1(7),fv2(7)
    real(wp) :: resg  ! resu of the 7-point Gauss formula
    real(wp) :: resk  ! resu of the 15-point Kronrod formula
    real(wp) :: reskh ! approximation to the mean value of f over (a,b), i.e. to i/(b-a)
    ! ----------------------------------------------------------------------------        

    ! Bound transformation
    ! ====================
    centr=0.5_wp*(a+b)
    hlgth=0.5_wp*(b-a)
    dhlgth=abs(hlgth)

    ! Compute the 15-point Kronrod approximation to the integral
    ! ==========================================================
    ! and estimate the absolute error

    fc=f(centr)
    resg=fc*wg(4)
    resk=fc*wgk(8)
    resabs=abs(resk)

    do j=1,3
       jtw=j*2
       absc=hlgth*xgk(jtw)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtw)=fval1
       fv2(jtw)=fval2
       fsum=fval1+fval2
       resg=resg+wg(j)*fsum
       resk=resk+wgk(jtw)*fsum
       resabs=resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    enddo

    do j=1,4
       jtwm1=j*2-1
       absc=hlgth*xgk(jtwm1)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtwm1)=fval1
       fv2(jtwm1)=fval2
       fsum=fval1+fval2
       resk=resk+wgk(jtwm1)*fsum
       resabs=resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    enddo

    reskh=resk * 0.5_wp
    resasc=wgk(8)*abs(fc-reskh)

    do j=1,7
       resasc=resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    enddo

    resu=resk*hlgth
    resabs=resabs*dhlgth
    resasc=resasc*dhlgth
    abserr=abs((resk-resg)*hlgth)

    if (resasc/=0.0_wp.and.abserr/=0.0_wp) then
       abserr=resasc*min(1.0_wp,(200.0_wp*abserr/resasc)**1.5_wp)
    endif

    if (resabs>tiny(resabs)/(50.0_wp*epsilon(resabs))) then
       abserr=max((epsilon(resabs)*50.0_wp)*resabs,abserr)
    endif

  end subroutine qk15

  !===============================================================================
  subroutine  qk21(f,a,b,resu,abserr,resabs,resasc)
  !===============================================================================
    !> QK21 carries out a 21 point Gauss-Kronrod quadrature rule [QUADPACK]
    !
    ! This routine approximates I=integral(a<=x<=b)f(x)dx
    ! with an error estimate, and J=integral(a<=x<=b)|f(x)|dx
    !
    !  Authors:
    !  Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
    !  QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
    !
    !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), external :: f         ! name of the integrand function routine
    real(wp), intent(in) :: a,b     ! integration bounds
    real(wp), intent(out) :: resu   ! estimated value of the integral
    real(wp), intent(out) :: abserr ! estimate of |I-resu|
    real(wp), intent(out) :: resabs ! approximation to the integral of the absolute value of f
    real(wp), intent(out) :: resasc ! approximation to the integral |f-I/(b-a)| over [a,b]
    ! ----------------------------------------------------------------------------        
    integer :: j,jtw,jtwm1
    real(wp) :: centr ! mid point of the interval
    real(wp) :: hlgth,dhlgth
    real(wp) :: absc  ! abscissa
    real(wp) :: fc,fsum
    real(wp) :: fval1,fval2 ! function values
    real(wp) :: fv1(10),fv2(10)
    real(wp) :: resg  ! resu of the 7-point Gauss formula
    real(wp) :: resk  ! resu of the 15-point Kronrod formula
    real(wp) :: reskh ! approximation to the mean value of f over (a,b), i.e. to i/(b-a)
    ! ----------------------------------------------------------------------------        

    ! Bound transformation
    ! ====================
    centr=0.5_wp*(a+b)
    hlgth=0.5_wp*(b-a)
    dhlgth=abs(hlgth)

    ! Compute the 21-point Kronrod approximation to the integral
    ! ==========================================================
    ! and estimate the absolute error

    resg=0.0_wp
    fc=f(centr)
    resk=wgk(11)*fc
    resabs=abs(resk)

    do j=1,5
       jtw=2*j
       absc=hlgth*xgk(jtw)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtw)=fval1
       fv2(jtw)=fval2
       fsum=fval1+fval2
       resg=resg+wg(j)*fsum
       resk=resk+wgk(jtw)*fsum
       resabs=resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    enddo

    do j=1,5
       jtwm1=2*j-1
       absc=hlgth*xgk(jtwm1)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtwm1)=fval1
       fv2(jtwm1)=fval2
       fsum=fval1+fval2
       resk=resk+wgk(jtwm1)*fsum
       resabs=resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    enddo

    reskh=resk*0.5_wp
    resasc=wgk(11)*abs(fc-reskh)

    do j=1,10
       resasc=resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    enddo

    resu=resk*hlgth
    resabs=resabs*dhlgth
    resasc=resasc*dhlgth
    abserr=abs((resk-resg)*hlgth)

    if (resasc/=0.0_wp.and.abserr/=0.0_wp) then
       abserr=resasc*min(1.0_wp,(200.0_wp*abserr/resasc)**1.5_wp)
    endif

    if (resabs>tiny(resabs)/(50.0_wp*epsilon(resabs))) then
       abserr=max((epsilon(resabs)*50.0_wp)*resabs,abserr)
    endif

  end subroutine qk21

  !===============================================================================
  subroutine  qk31(f,a,b,resu,abserr,resabs,resasc)
  !===============================================================================
    !> QK31 carries out a 31 point Gauss-Kronrod quadrature rule [QUADPACK]
    !
    ! This routine approximates I=integral(a<=x<=b)f(x)dx
    ! with an error estimate, and J=integral(a<=x<=b)|f(x)|dx
    !
    !  Authors:
    !  Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
    !  QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
    !
    !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), external :: f         ! name of the integrand function routine
    real(wp), intent(in) :: a,b     ! integration bounds
    real(wp), intent(out) :: resu   ! estimated value of the integral
    real(wp), intent(out) :: abserr ! estimate of |I-resu|
    real(wp), intent(out) :: resabs ! approximation to the integral of the absolute value of f
    real(wp), intent(out) :: resasc ! approximation to the integral |f-I/(b-a)| over [a,b]
    ! ----------------------------------------------------------------------------        
    integer :: j,jtw,jtwm1
    real(wp) :: centr ! mid point of the interval
    real(wp) :: hlgth,dhlgth
    real(wp) :: absc  ! abscissa
    real(wp) :: fc,fsum
    real(wp) :: fval1,fval2 ! function values
    real(wp) :: fv1(15),fv2(15)
    real(wp) :: resg  ! resu of the 7-point Gauss formula
    real(wp) :: resk  ! resu of the 15-point Kronrod formula
    real(wp) :: reskh ! approximation to the mean value of f over (a,b), i.e. to i/(b-a)
    ! ----------------------------------------------------------------------------        

    ! Bound transformation
    ! ====================
    centr=0.5_wp*(a+b)
    hlgth=0.5_wp*(b-a)
    dhlgth=abs(hlgth)

    ! Compute the 31-point Kronrod approximation to the integral
    ! ==========================================================
    ! and estimate the absolute error

    fc=f(centr)
    resg=wg(8)*fc
    resk=wgk(16)*fc
    resabs=abs(resk)

    do j=1,7
       jtw=j*2
       absc=hlgth*xgk(jtw)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtw)=fval1
       fv2(jtw)=fval2
       fsum=fval1+fval2
       resg=resg+wg(j)*fsum
       resk=resk+wgk(jtw)*fsum
       resabs=resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    enddo

    do j=1,8
       jtwm1=j*2-1
       absc=hlgth*xgk(jtwm1)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtwm1)=fval1
       fv2(jtwm1)=fval2
       fsum=fval1+fval2
       resk=resk+wgk(jtwm1)*fsum
       resabs=resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    enddo

    reskh=resk*0.5_wp
    resasc=wgk(16)*abs(fc-reskh)

    do j=1,15
       resasc=resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    enddo

    resu=resk*hlgth
    resabs=resabs*dhlgth
    resasc=resasc*dhlgth
    abserr=abs((resk-resg)*hlgth)

    if (resasc/=0.0_wp.and.abserr/=0.0_wp) &
         abserr=resasc*min(1.0_wp,(200.0_wp*abserr/resasc)**1.5_wp)

    if (resabs>tiny(resabs)/(50.0_wp*epsilon(resabs))) then
       abserr=max((epsilon(resabs)*50.0_wp)*resabs,abserr)
    endif

  end subroutine qk31

  !===============================================================================
  subroutine  qk41(f,a,b,resu,abserr,resabs,resasc)
  !===============================================================================
    !> QK41 carries out a 41 point Gauss-Kronrod quadrature rule [QUADPACK]
    !
    ! This routine approximates I=integral(a<=x<=b)f(x)dx
    ! with an error estimate, and J=integral(a<=x<=b)|f(x)|dx
    !
    !  Authors:
    !  Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
    !  QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
    !
    !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), external :: f         ! name of the integrand function routine
    real(wp), intent(in) :: a,b     ! integration bounds
    real(wp), intent(out) :: resu   ! estimated value of the integral
    real(wp), intent(out) :: abserr ! estimate of |I-resu|
    real(wp), intent(out) :: resabs ! approximation to the integral of the absolute value of f
    real(wp), intent(out) :: resasc ! approximation to the integral |f-I/(b-a)| over [a,b]
    ! ----------------------------------------------------------------------------        
    integer :: j,jtw,jtwm1
    real(wp) :: centr ! mid point of the interval
    real(wp) :: hlgth,dhlgth
    real(wp) :: absc  ! abscissa
    real(wp) :: fc,fsum
    real(wp) :: fval1,fval2 ! function values
    real(wp) :: fv1(20),fv2(20)
    real(wp) :: resg  ! resu of the 7-point Gauss formula
    real(wp) :: resk  ! resu of the 15-point Kronrod formula
    real(wp) :: reskh ! approximation to the mean value of f over (a,b), i.e. to i/(b-a)
    ! ----------------------------------------------------------------------------        

    ! Bound transformation
    ! ====================
    centr=0.5_wp*(a+b)
    hlgth=0.5_wp*(b-a)
    dhlgth=abs(hlgth)

    ! Compute the 41-point Kronrod approximation to the integral
    ! ==========================================================
    ! and estimate the absolute error

    resg=0.0_wp
    fc=f(centr)
    resk=wgk(21)*fc
    resabs=abs(resk)

    do j=1,10
       jtw=j*2
       absc=hlgth*xgk(jtw)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtw)=fval1
       fv2(jtw)=fval2
       fsum=fval1+fval2
       resg=resg+wg(j)*fsum
       resk=resk+wgk(jtw)*fsum
       resabs=resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    enddo

    do j=1,10
       jtwm1=j*2-1
       absc=hlgth*xgk(jtwm1)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtwm1)=fval1
       fv2(jtwm1)=fval2
       fsum=fval1+fval2
       resk=resk+wgk(jtwm1)*fsum
       resabs=resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    enddo

    reskh=resk*0.5_wp
    resasc=wgk(21)*abs(fc-reskh)

    do j=1,20
       resasc=resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    enddo

    resu=resk*hlgth
    resabs=resabs*dhlgth
    resasc=resasc*dhlgth
    abserr=abs((resk-resg)*hlgth)

    if (resasc/=0.0_wp.and.abserr/=0.0_wp) &
         abserr=resasc*min(1.0_wp,(200.0_wp*abserr/resasc)**1.5_wp)

    if (resabs>tiny(resabs)/(50.0_wp*epsilon(resabs))) then
       abserr=max((epsilon(resabs)*50.0_wp)*resabs,abserr)
    endif

  end subroutine qk41

  !===============================================================================
  subroutine  qk51(f,a,b,resu,abserr,resabs,resasc)
  !===============================================================================
    !> QK51 carries out a 51 point Gauss-Kronrod quadrature rule [QUADPACK]
    !
    ! This routine approximates I=integral(a<=x<=b)f(x)dx
    ! with an error estimate, and J=integral(a<=x<=b)|f(x)|dx
    !
    !  Authors:
    !  Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
    !  QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
    !
    !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), external :: f         ! name of the integrand function routine
    real(wp), intent(in) :: a,b     ! integration bounds
    real(wp), intent(out) :: resu   ! estimated value of the integral
    real(wp), intent(out) :: abserr ! estimate of |I-resu|
    real(wp), intent(out) :: resabs ! approximation to the integral of the absolute value of f
    real(wp), intent(out) :: resasc ! approximation to the integral |f-I/(b-a)| over [a,b]
    ! ----------------------------------------------------------------------------        
    integer :: j,jtw,jtwm1
    real(wp) :: centr ! mid point of the interval
    real(wp) :: hlgth,dhlgth
    real(wp) :: absc  ! abscissa
    real(wp) :: fc,fsum
    real(wp) :: fval1,fval2 ! function values
    real(wp) :: fv1(25),fv2(25)
    real(wp) :: resg  ! resu of the 7-point Gauss formula
    real(wp) :: resk  ! resu of the 15-point Kronrod formula
    real(wp) :: reskh ! approximation to the mean value of f over (a,b), i.e. to i/(b-a)
    ! ----------------------------------------------------------------------------        

    ! Bound transformation
    ! ====================
    centr=0.5_wp*(a+b)
    hlgth=0.5_wp*(b-a)
    dhlgth=abs(hlgth)

    ! Compute the 51-point Kronrod approximation to the integral
    ! ==========================================================
    ! and estimate the absolute error

    fc=f(centr)
    resg=wg(13)*fc
    resk=wgk(26)*fc
    resabs=abs(resk)

    do j=1,12
       jtw=j*2
       absc=hlgth*xgk(jtw)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtw)=fval1
       fv2(jtw)=fval2
       fsum=fval1+fval2
       resg=resg+wg(j)*fsum
       resk=resk+wgk(jtw)*fsum
       resabs=resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    enddo

    do j=1,13
       jtwm1=j*2-1
       absc=hlgth*xgk(jtwm1)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtwm1)=fval1
       fv2(jtwm1)=fval2
       fsum=fval1+fval2
       resk=resk+wgk(jtwm1)*fsum
       resabs=resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    enddo

    reskh=resk*0.5_wp
    resasc=wgk(26)*abs(fc-reskh)

    do j=1,25
       resasc=resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    enddo

    resu=resk*hlgth
    resabs=resabs*dhlgth
    resasc=resasc*dhlgth
    abserr=abs((resk-resg)*hlgth)

    if (resasc/=0.0_wp.and.abserr/=0.0_wp) then
       abserr=resasc*min(1.0_wp,(200.0_wp*abserr/resasc)**1.5_wp)
    endif

    if (resabs>tiny(resabs)/(50.0_wp*epsilon(resabs))) then
       abserr=max((epsilon(resabs)*50.0_wp)*resabs,abserr)
    endif

  end subroutine qk51

  !===============================================================================
  subroutine  qk61(f,a,b,resu,abserr,resabs,resasc)
  !===============================================================================
    !> QK61 carries out a 61 point Gauss-Kronrod quadrature rule [QUADPACK]
    !
    ! This routine approximates I=integral(a<=x<=b)f(x)dx
    ! with an error estimate, and J=integral(a<=x<=b)|f(x)|dx
    !
    !  Authors:
    !  Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
    !  QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
    !
    !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), external :: f         ! name of the integrand function routine
    real(wp), intent(in) :: a,b     ! integration bounds
    real(wp), intent(out) :: resu   ! estimated value of the integral
    real(wp), intent(out) :: abserr ! estimate of |I-resu|
    real(wp), intent(out) :: resabs ! approximation to the integral of the absolute value of f
    real(wp), intent(out) :: resasc ! approximation to the integral |f-I/(b-a)| over [a,b]
    ! ----------------------------------------------------------------------------        
    integer :: j,jtw,jtwm1
    real(wp) :: centr ! mid point of the interval
    real(wp) :: hlgth,dhlgth
    real(wp) :: absc  ! abscissa
    real(wp) :: fc,fsum
    real(wp) :: fval1,fval2 ! function values
    real(wp) :: fv1(30),fv2(30)
    real(wp) :: resg  ! resu of the 7-point Gauss formula
    real(wp) :: resk  ! resu of the 15-point Kronrod formula
    real(wp) :: reskh ! approximation to the mean value of f over (a,b), i.e. to i/(b-a)
    ! ----------------------------------------------------------------------------        

    ! Bound transformation
    ! ====================
    centr=0.5_wp*(b+a)
    hlgth=0.5_wp*(b-a)
    dhlgth=abs(hlgth)

    ! Compute the 61-point Kronrod approximation to the integral
    ! ==========================================================
    ! and estimate the absolute error

    resg=0.0_wp
    fc=f(centr)
    resk=wgk(31)*fc
    resabs=abs(resk)

    do j=1,15
       jtw=j*2
       absc=hlgth*xgk(jtw)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtw)=fval1
       fv2(jtw)=fval2
       fsum=fval1+fval2
       resg=resg+wg(j)*fsum
       resk=resk+wgk(jtw)*fsum
       resabs=resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    enddo

    do j=1,15
       jtwm1=j*2-1
       absc=hlgth*xgk(jtwm1)
       fval1=f(centr-absc)
       fval2=f(centr+absc)
       fv1(jtwm1)=fval1
       fv2(jtwm1)=fval2
       fsum=fval1+fval2
       resk=resk+wgk(jtwm1)*fsum
       resabs=resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    enddo

    reskh=resk * 0.5_wp
    resasc=wgk(31)*abs(fc-reskh)

    do j=1,30
       resasc=resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    enddo

    resu=resk*hlgth
    resabs=resabs*dhlgth
    resasc=resasc*dhlgth
    abserr=abs((resk-resg)*hlgth)

    if (resasc/=0.0_wp.and.abserr/=0.0_wp) then
       abserr=resasc*min(1.0_wp,(200.0_wp*abserr/resasc)**1.5_wp)
    endif

    if (resabs>tiny(resabs)/(50.0_wp* epsilon(resabs))) then
       abserr=max((epsilon(resabs)*50.0_wp)*resabs,abserr)
    endif

  end subroutine qk61

  !===============================================================================
  subroutine qsort(limit,last,maxerr,ermax,elist,iord,nrmax)
  !===============================================================================
    !> QSORT maintains the order of a list of local error estimates [QUADPACK]
    !
    ! This routine maintains the descending ordering in the list of the 
    ! local error estimates resuing from the interval subdivision process. 
    ! At each call two error estimates are inserted using the sequential 
    ! search top-down for the largest error estimate and bottom-up for the
    ! smallest error estimate.
    !
    ! iord(last): The first K elements contain pointers to the error estimates
    ! such that elist(iord(1)) through elist(iord(k)) form a decreasing sequence,
    ! with k=last if  last<=(limit/2+2), and otherwise k=limit+1-last.
    !
    !  Authors:
    !  Robert Piessens, Elise de Doncker-Kapenger, Christian Ueberhuber, David Kahaner,
    !  QUADPACK, a Subroutine Package for Automatic Integration, Springer Verlag, 1983
    !
  !===============================================================================
    implicit none
    integer, intent(in) :: limit ! maximum number of error estimates the list can contain
    integer, intent(in) :: last          ! current number of error estimates
    integer, intent(inout) :: maxerr     ! index in the list of the nrmax-th largest error
    real(wp), intent(out) :: ermax       ! nrmax-th largest error=elist(maxerr)
    real(wp), intent(in) :: elist(last)  ! contains the error estimates
    integer, intent(inout) :: iord(last) ! pointers to the error estimates
    integer, intent(inout) :: nrmax      ! number of err
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    integer :: ibeg,isucc,jbnd,jupbn
    real(wp) :: errmin,errmax
    ! ----------------------------------------------------------------------------

    ! Check whether the list contains more than two error estimates
    ! =============================================================
    if (last<=2) then
       iord(1)=1
       iord(2)=2
       go to 90
    endif

    ! This part of the routine is only executed if, due to a
    ! difficult integrand, subdivision increased the error
    ! estimate. In the normal case the insert procedure should
    ! start after the nrmax-th largest error estimate.

    errmax=elist(maxerr)

    do i=1,nrmax-1

       isucc=iord(nrmax-1)

       if (errmax<=elist(isucc)) then
          exit
       endif

       iord(nrmax)=isucc
       nrmax=nrmax-1

    enddo

    ! Compute the number of elements in the list to be maintained
    ! in descending order. This number depends on the number of
    ! subdivisions still allowed.
 
    jupbn=last

    if ((limit/2+2) < last ) then
       jupbn=limit+3-last
    endif

    errmin=elist(last)

    ! Insert errmax by traversing the list top-down,starting
    ! comparison from the element elist(iord(nrmax+1))
    ! ======================================================

    jbnd=jupbn-1
    ibeg=nrmax+1

    do i=ibeg,jbnd
       isucc=iord(i)
       if (elist(isucc) <= errmax ) then
          go to 60
       endif
       iord(i-1)=isucc
    enddo

    iord(jbnd)=maxerr
    iord(jupbn)=last
    go to 90

    ! Insert errmin by traversing the list bottom-up
    ! ==============================================
60  continue

    iord(i-1)=maxerr
    k=jbnd

    do j=i,jbnd
       isucc=iord(k)
       if (errmin < elist(isucc)) then
          go to 80
       endif
       iord(k+1)=isucc
       k=k-1
    enddo

    iord(i)=last
    go to 90

80  continue

    iord(k+1)=last

    ! Set maxerr and ermax
    ! ====================
90  continue

    maxerr=iord(nrmax)
    ermax=elist(maxerr)

  end subroutine qsort

end module mod_quadpack
