!=================================================================================
module mod_fluid
!=================================================================================
  !> Module to define fluid properties & constants [use feos_xxx.ini]
!=================================================================================
  use precision
  implicit none
  ! ------------------------------------------------------------------------------
  character(len=10) :: fluidname
  ! ------------------------------------------------------------------------------
  character(len=3) :: eos_type   ! read in param.ini
  character(len=1) :: visc_type  ! read in param.ini
  ! ------------------------------------------------------------------------------
  real(wp) :: cpfg,cvfg
  real(wp) :: gam1,igm1
  ! ------------------------------------------------------------------------------
  ! critical quantities (temperature, pressure, density & compressibility factor)
  real(wp) :: Tc,pc,roc,zc,roc0
  ! fluid constants
  real(wp) :: pmol,mdm,teb,nexp,cvinf,om,rg,Pr,gam
  ! parameters of Span-Wagner model
  real(wp) :: n01,n02,n03,n04,n05,n06,n07,n08,n09,n10,n11,n12
  real(wp) :: eta1,eta2,eta3,eta4
  ! ------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine read_eos
  !===============================================================================
    !> Read file feos_xxx.ini
  !===============================================================================
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer, parameter :: uni=60
    logical :: iexist
    ! ----------------------------------------------------------------------------

    ! Check file
    ! ==========
    inquire(file='feos_'//trim(fluidname)//'.ini', exist=iexist)

    if (.not.iexist) call mpistop('Feosdon file does not exist!',0)

    open(uni,file='feos_'//trim(fluidname)//'.ini', form='formatted')

    ! Read file
    ! =========
    read(uni,*) ! Header
    read(uni,'(45X,E13.6)') Tc
    read(uni,'(45X,E13.6)') pc
    read(uni,'(45X,E13.6)') roc
    read(uni,'(45X,E13.6)') zc
    read(uni,'(45X,E13.6)') pmol
    read(uni,'(45X,E13.6)') mdm
    read(uni,'(45X,E13.6)') teb
    read(uni,'(45X,E13.6)') om
    read(uni,'(45X,E13.6)') cvinf
    read(uni,*)
    read(uni,'(45X,E13.6)') nexp
    read(uni,'(45X,E13.6)') rg
    read(uni,'(45X,E13.6)') Pr
    read(uni,'(45X,E13.6)') gam

    if (eos_type(1:2).eq.'sw') then
       read(uni,*) ! Header
       read(uni,'(45X,E15.6)') n01
       read(uni,'(45X,E15.6)') n02
       read(uni,'(45X,E15.6)') n03
       read(uni,'(45X,E15.6)') n04
       read(uni,'(45X,E15.6)') n05
       read(uni,'(45X,E15.6)') n06
       read(uni,'(45X,E15.6)') n07
       read(uni,'(45X,E15.6)') n08
       read(uni,'(45X,E15.6)') n09
       read(uni,'(45X,E15.6)') n10
       read(uni,'(45X,E15.6)') n11
       read(uni,'(45X,E15.6)') n12
       read(uni,'(45X,E15.6)') tc
       read(uni,'(45X,E15.6)') pc
       read(uni,'(45X,E15.6)') roc
       read(uni,*) ! Header
       read(uni,'(45X,E15.6)') eta1
       read(uni,'(45X,E15.6)') eta2
       read(uni,'(45X,E15.6)') eta3
       read(uni,'(45X,E15.6)') eta4
    endif
    close(uni)
    
  end subroutine read_eos
  
end module mod_fluid
