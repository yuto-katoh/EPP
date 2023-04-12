!-----------------------------------------------------------------------
module p__init
  !-----------------------------------------------------------------------
  use c__eprm,  only: lx,lxa1,npa1, ialt
  use v__in,    only: q,tpi,pi, t_live
  use p__initmg,only: initmg,delx
  use p__initcs,only: initcs
  implicit none
  !
  real(kind=8), save ::T,alp,sigmas_ion(3,10000), sigmas_dif(1000,10000)
  real(kind=8), save ::sigmas_elas(3,10000),sigmas_tot(3,10000)
  !
  private
  public :: T, init,alp, sigmas_ion,sigmas_elas,sigmas_tot,sigmas_dif
  !
contains
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine init(xx,vp1,pm1,lor,bf,pind,alt,fce,lat,hh,pn_all)
    !-----------------------------------------------------------------------
    real(kind=8)   ,intent(out) ::  xx(npa1),vp1(3,npa1)
    real(kind=8)   ,intent(out) ::  pm1(3,npa1),lor(npa1)
    real(kind=8)   ,intent(out) ::  bf(lxa1),alt(lxa1)
    real(kind=8)   ,intent(out) ::  fce(lxa1),lat(lxa1),hh(lxa1)
    real(kind=8)   ,intent(out) ::  pn_all(3,lxa1)
    integer(kind=4),intent(out) ::  pind(npa1)
    !
    real(kind=8),parameter ::  kb = 1.381d-23 , vc = 2.9979d8
    real(kind=8),parameter ::  ev = 1.16d4 , me = 9.11d-31
    
    !
    integer(kind=4),external :: mt_init
    real(kind=8)   ,external :: mtrand_real8
    !
    real(kind=8) ::  phs,vpp,aa,factp,vges2,aaa,Tc,tot,phi,scale
    integer      ::  ii,jj,kk,ll,iseed,istat
    real(kind=8) ::  a,b,c

    !-----------------------------------------------------------------------
    factp = 180.d0/pi
    iseed = 634
    istat = mt_init(iseed)
    !-----------------------------------------------------------------------
    call initmg( bf,alt,fce,lat,hh,pn_all)
    call initcs(sigmas_ion, sigmas_elas, sigmas_dif, sigmas_tot)
    !-----------------------------------------------------------------------

    T = 10.d0   !Energy in keV
!    print*, "Enter Energy [keV]:"
!    read (*,*) T
    print*, "Energy [keV]:",T
    T = T*ev*1d3

    
    !alp = 72.33
    !alp = 30.0
!    print*, "Enter initial pitch angle [°]:"
!    read (*,*) alp
    alp = 70.0
    print*, "Initial pitch angle [°]:",alp
    alp = alp/factp
    
       do ii=1,npa1

       aa = me*vc**2/( kb*T + me*vc**2 )
       vpp = sqrt(1.d0 - aa**2)

       xx(ii) = 2.d0*delx
       vp1(1,ii) = vpp*cos(alp)
       vp1(2,ii) = vpp*sin(alp)
       vp1(3,ii) = vp1(2,ii)**2/2.d0

       lor(ii) = 1.d0/sqrt(1.d0-vp1(1,ii)**2-vp1(2,ii)**2)
       pm1(1,ii) = vp1(1,ii)*lor(ii)
       pm1(2,ii) = vp1(2,ii)*lor(ii)
       pm1(3,ii) = 0.d0

       pind(ii) = 0
    end do





    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    return
    !-----------------------------------------------------------------------
  end subroutine init
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
end module p__init
