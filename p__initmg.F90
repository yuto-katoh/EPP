!-----------------------------------------------------------------------
module p__initmg
  !-----------------------------------------------------------------------
  use c__eprm  ,only: lxa1, ialt, resf, naopt
  use v__in    ,only: pi,tpi
  implicit none
  !
  real(kind=8)   ,save :: delx
  !
  private
  public :: initmg,delx
  !
contains
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine initmg( bf,alt,fce,lat,hh, pn_all )
    !-----------------------------------------------------------------------
    real(kind=8),intent(out) ::  bf(lxa1)
    real(kind=8),intent(out) ::  alt(lxa1),fce(lxa1),lat(lxa1),hh(lxa1)
    real(kind=8),intent(inout) ::  pn_all(3,lxa1)
    !
    real(kind=8),parameter ::  me = 9.10939d-31
    real(kind=8),parameter ::  RE = 6.378d6
    real(kind=8),parameter ::  R0 = RE + dble(ialt)* 1.d3

    real(kind=8),parameter ::  bb_eq = 1.1590d-7   ! [nT]@L=6.45
    real(kind=8),parameter ::  L = 6.45d0    

    !    real(kind=8),parameter ::  bb_eq = 4.8594d-7   ! [nT]@L=4
    !    real(kind=8),parameter ::  L = 4.d0
    !    real(kind=8),parameter ::  bb_eq = 1.4398d-7   ! [nT]@L=6
    !    real(kind=8),parameter ::  L = 6.d0
    real(kind=8),parameter ::  qq = 1.602d-19
    real(kind=8),parameter ::  vc = 2.99792458d8
    real(kind=8)    ::  Brad,Bram,fact1,rr0,cos2_th0,th0,yy0,sm
    real(kind=8)    ::  bb0,wce,slen,ds,th,rr,factp,dth,ss
    integer(kind=4) ::  ii,jj,ind,lxb1,lxb2

    ! file import
    character (LEN=75) :: temp
    real(kind=8), dimension(1000) ::  in_hei, in_O, in_N2, in_O2
    real(kind=8), dimension(3,1000) ::  in_all

    ! profile calc
    real(kind=8) :: hf

    !interpolation
    real(kind=8) :: ds_,m
    integer :: a,b,zpos
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    factp = 180.d0/pi
    rr0 = R0/RE
    !-----------------------------------------------------------------------
    cos2_th0 = rr0/L
    th0 = acos( sqrt(cos2_th0) )
    !### th0 : latitude (0 deg at the magnetic equator)
    !-----------------------------------------------------------------------
    yy0 = sqrt(3.d0)*sin(th0)
    fact1 = log(yy0+sqrt(1.d0+yy0**2)) + yy0*sqrt(1.d0+yy0**2)
    sm = rr0/2.d0/sqrt(3.d0)/cos2_th0*fact1
    !### rr0 : altitude of the base point of the field line (500 km)
    !### sm  : length of field line from rr0 to magnetic equator
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    fact1 = cos(th0)**6
    fact1 = 1.d0/fact1
    Brad = -2.d0*sin(th0)*fact1
    Bram = cos(th0)*fact1
    bb0 = sqrt( Brad**2 + Bram**2 )*bb_eq
    wce = qq*bb0/me
    slen = vc/wce
    !-----------------------------------------------------------------------
    rr0 = rr0*RE/slen
    sm  = sm*RE/slen
    !-----------------------------------------------------------------------
    ds = 1.d2/slen    *1e1/resf !Not sure about this change
    delx = ds
    th = th0
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    do ii=1,lxa1
       !-----------------------------------------------------------------------
       ss = sm + ds*dble(ii)
       dth = ds*cos2_th0/rr0/sqrt(1.d0+3.d0*sin(th)**2)/cos(th)
       dth = dth
       if(ii.eq.1) dth = 0.d0
       th = th + dth
       !-----------------------------------------------------------------------
       rr = rr0*cos(th)**2/cos2_th0
       !-----------------------------------------------------------------------
       fact1 = cos(th)**6
       fact1 = 1.d0/fact1
       !-----------------------------------------------------------------------
       Brad = -2.d0*sin(th)*fact1
       Bram = cos(th)*fact1
       bf(ii) = sqrt( Brad**2 + Bram**2 )*bb_eq
       fce(ii) = qq*bf(ii)/me/tpi
       lat(ii) = th*factp
       !      hh(ii)  = ds*dble(ii)*slen
       hh(ii)  = ds*dble(ii)
       alt(ii) = (rr*slen - RE)*1.d-3
       !-----------------------------------------------------------------------
    end do
    !-----------------------------------------------------------------------

    !BUILD NEUTRAL DENSITY PROFILES using msise model file input (must have 1000 values)

    !we will now read in the file with the number densities!
    open(18, file="data/msis.lst")
    do ii=1,27
       read(18,*) temp !discard the first 27 lines, then save rest
    end do
    do ii=1,1000
       read(18,*) in_hei(ii), &
            &in_all(1,ii), in_all(2,ii), in_all(3,ii)
    end do

    ds_ = in_hei(2)-in_hei(1)
    zpos = 1d0 - in_hei(1)/ds_ 

    !now perform linear interpolation to get the new (fine) profiles
    do ii=1,lxa1
       hf = alt(ii) !get the height
       m = zpos + hf/ds_ !calculate index
       a = int(m)
       a = min(max(a,1),999)
       b = a+1

       do jj=1,3
          pn_all(jj,ii) = (in_all(jj,a) + &
               (in_all(jj,b)-in_all(jj,a))*(m-dble(a))) * 1d6 !cm-3 to m-3
       end do
    end do

    pn_all(:,:) = naopt * pn_all(:,:)

    !-----------------------------------------------------------------------
    bb0 = bf(1)
    do ii=1,lxa1
       bf(ii) = bf(ii)/bb0
       write(91,'(I10,7(1PE12.4))') & 
            &       ii,alt(ii),bf(ii),fce(ii),lat(ii),hh(ii),pn_all(1,ii) 
    end do
    !-----------------------------------------------------------------------
    return
    !-----------------------------------------------------------------------
  end subroutine initmg
  !-----------------------------------------------------------------------
end module p__initmg
