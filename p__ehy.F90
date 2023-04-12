!-----------------------------------------------------------------------
module ehy_main
  !-----------------------------------------------------------------------
  use c__eprm  , only: istep,lxa1,npa1,ialt,lx,mfopt,resfout,numtout,live_res,live_max
  use v__in    , only: v__in__ini,pi,delt,t_live
  use p__init  , only: init, T, alp
  use p__initmg, only: delx
  use p__passv , only: passv
  use p__passx , only: passx
  implicit none
  !-----------------------------------------------------------------------
  real(kind=8)    ::  bf(lxa1)
  real(kind=8)    ::  xx(npa1),vp1(3,npa1),pm1(3,npa1),lor(npa1),mask(npa1)
  integer(kind=4) ::  pind(npa1)
  real(kind=8)    ::  alt(lxa1),fce(lxa1),lat(lxa1),hh(lxa1)
  real(kind=8)    ::  pn_all(3,lxa1)
  integer(kind=4) ::  np_alt(ialt*resfout,numtout)

  real(kind=8),parameter ::  ev = 1.16d4
  !
  real(kind=8)    ::  vpara,vperp,vpp,pxx,palt,palp
  real(kind=8)    ::  factp,aa,inen
  integer(kind=4) ::  ii,jj,ip0,ifile,ic,id
  !
contains
  !-----------------------------------------------------------------------
  subroutine ehy
    !-----------------------------------------------------------------------
    np_alt(:,:) = 0
    !-----------------------------------------------------------------------
    call v__in__ini
    !-----------------------------------------------------------------------

    call init( xx,vp1,pm1,lor,bf,pind,alt,fce,lat,hh, pn_all)
    !-----------------------------------------------------------------------

    factp = 180.d0/pi
    inen = sum(vp1(1,:)**2 + vp1(2,:)**2)

    do
       do ii = 1,istep
          call passv( xx,vp1,pm1,lor,bf,pind )
          call passx( xx,vp1,alt,np_alt,pind, pn_all, fce )
         
       end do

        

       !if the simulation is not finished now, run another cycle
       if(ANY(pind==0)) then
          mask = 0
          where (pind ==0) mask = 1
          print *, int(sum(mask)), npa1/sum(mask)*sum(mask*vp1(1,:)**2 &
               + mask*vp1(2,:)**2)/inen * 1e2
       else
          print *, char(7)
          exit
       end if
    end do

    !save the result
    id = int(2e8+1e8*int(mfopt*1.1)+1e6*int(alp*factp)+int(T/ev))
    !id = 60
    write(id,*)  int(T/ev), npa1, mfopt, int(alp*factp), numtout, live_max, resfout !header
    do ii = 1,ialt*resfout
       write(id,*) dble(ii)/resfout,np_alt(ii,:)
    end do
    !-----------------------------------------------------------------------
    return
    !-----------------------------------------------------------------------
  end subroutine ehy
  !-----------------------------------------------------------------------
end module ehy_main
