!-----------------------------------------------------------------------
module p__passx
  !-----------------------------------------------------------------------
  use c__eprm,  only: lx1,lx2,npa1,lxa1,ialt,resfout,live_max,live_res,numtout
  use v__in,    only: delt, tpi, pi, t_live
  use p__initmg,only: delx
  use p__init,  only: sigmas_ion, sigmas_elas, sigmas_tot, sigmas_dif
  implicit none
  !
  private
  public :: passx
  !
contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine passx( xp1,vp1,alt,np_alt,pind, pn_all, fce )
    !-----------------------------------------------------------------------

    real(kind=8)   ,intent(inout) ::  vp1(3,npa1),xp1(npa1)
    real(kind=8)   ,intent(inout) ::  alt(lxa1), pn_all(3,lxa1)
    real(kind=8)   ,intent(in)    ::  fce(lxa1)
    integer(kind=4),intent(inout) ::  pind(npa1)
    integer(kind=4),intent(inout) ::  np_alt(ialt*resfout,numtout)
    !
    real(kind=8),external :: mtrand_real8

    real(kind=8),parameter ::  kb = 1.381d-23 , vc = 2.9979d8
    real(kind=8),parameter ::  ev = 1.16d4 , me = 9.11d-31

    real(kind=8),parameter :: E_pair(3) = (/ 35.0, 37.0, 33.0 /)
    character(len=2),parameter :: Species_N(3) = (/ "O ","N2","O2" /)

    real(kind=8) ::   aa,bb,raa1,rba1,half,hdelt
    real(kind=8) ::   p_dice,p_i,p_tot
    real(kind=8) ::   sigma_ion(3), sigma_elas(3), ratio, ratios(3)
    real(kind=8) ::   sigma, sigma_tot(3), mac_sigma, wce, vges2, T
    real(kind=8) ::   Eind, aaa,  vfact, Tloss
    logical      ::   mask(1000)
    integer      ::   ii,ic,id,id2,jj,kk,ll,ind,dectf, pos
    real(kind=8) ::   phi, phi_ch, b_dice

    !-----------------------------------------------------------------------
    half = 0.5d0
    wce = fce(1) *tpi !scaling factor for delta_t
    !-----------------------------------------------------------------------
    do ii = 1,npa1
       
       !-----------------------------------------------------------------------
       if( pind(ii).eq.0 ) then
          !-----------------------------------------------------------------------
          
          xp1(ii) = xp1(ii) + vp1(1,ii)*delt(ii)
          t_live(ii) = t_live(ii) + delt(ii)/wce*1.d3 !increase livetime (ms)

          !-----------------------------------------------------------------------
          aa   = xp1(ii)/delx
          raa1 = dint( aa+half )
          ic   = raa1

          !output place and time of the first 10 particles
!          if((ii<=10).and.(alt(ic)<ialt)) then
!             write(100+ii,*) alt(ic), t_live(ii)
!          endif

          
          !-----------------------------------------------------------------------
          if( ic.le.lx1-1 ) then
             pind(ii) = 1
          else if( ic.gt.lx2+1 ) then
             pind(ii) = 2
          endif

          !destroy particle if lifetime gets too big
          if(t_live(ii).ge.live_max) then
             pind(ii) = 3
          endif
          !-----------------------------------------------------------------------
          vges2 = vp1(1,ii)**2 + vp1(2,ii)**2
          aaa = sqrt(1-vges2)
          T = me*vc**2/kb*(1/aaa-1)  !get kin. energy
          Eind = 500*log(T/ev) !index in array of sigma (logarithmic scale)
          ind = int(Eind)
          phi = acos(vp1(1,ii)/sqrt(vges2)) !pitch angle
          
          !check if some collision is happening
          sigma_tot = sigmas_tot(:,ind)  + &
               (Eind-dble(ind))*(sigmas_tot(:,ind+1)-sigmas_tot(:,ind))
          mac_sigma = dot_product(sigma_tot, pn_all(:,ic))
          p_tot = 1 - exp(-mac_sigma * sqrt(vges2)*vc * delt(ii)/wce)
          p_dice = mtrand_real8()

          
          if((p_dice.le.p_tot)) then

             sigma_elas = sigmas_elas(:,ind)  + &
                  (Eind-dble(ind))*(sigmas_elas(:,ind+1)-sigmas_elas(:,ind))
             sigma_ion  = sigmas_ion(:,ind)  + &
                  (Eind-dble(ind))*(sigmas_ion(:,ind+1)-sigmas_ion(:,ind))

             !now see if its some ionization
             ratio = dot_product(sigma_ion, pn_all(:,ic))/mac_sigma
             p_dice = mtrand_real8()

             
             if((p_dice.le.ratio)) then
                ratios = sigma_ion/sum(sigma_ion)
                p_dice = mtrand_real8()

                
                if(p_dice.le.sum(ratios(1:3))) jj = 3
                if(p_dice.le.sum(ratios(1:2))) jj = 2
                if(p_dice.le.sum(ratios(1:1))) jj = 1

                !record collision
                id = int(alt(ic)*resfout)
                id2 = max(min(int(t_live(ii)/live_res)+1,numtout),1)
                np_alt(id, id2) = np_alt(id, id2) + 1
                

                if(ii.eq.1) then  !print information about collision of particle 1
                   print*, alt(ic),int(T/ev),Species_N(jj),sigmas_ion(jj,ind)* &
                        1d4,pn_all(jj,ic)*1d-6,int(-E_pair(jj)), phi*360/tpi, &
                        t_live(ii)
                end if

                T = T - E_pair(jj)*ev  !subtract ionization energy
                if(T/ev.le.20) then  !destroy particle if energy <20eV
                   pind(ii) = 4
                   exit
                endif
                
                !change velocity
                aaa = me*vc**2/(kb*T+me*vc**2)
                vfact = sqrt(1-aaa**2)/sqrt(vges2)
                vp1(1,ii) = vp1(1,ii)*vfact
                vp1(2,ii) = vp1(2,ii)*vfact
                vges2 = vp1(1,ii)**2 + vp1(2,ii)**2

             end if
             !end of ionoization -> now elastic
             
             !change pitch angle!
             b_dice = mtrand_real8()*tpi
             p_dice = mtrand_real8()

             mask = sigmas_dif(:,ind).ge.p_dice
             pos =  minloc(sigmas_dif(:,ind), DIM=1, MASK=mask)
             phi_ch = (dble(pos)-1)/1000*pi
             phi = acos(-sin(phi)*cos(b_dice)*sin(phi_ch)+cos(phi)*cos(phi_ch))
             vp1(1,ii) = cos(phi)*sqrt(vges2)
             vp1(2,ii) = sin(phi)*sqrt(vges2)


             !if(ii.eq.1) then  !print information about collision of particle 1
             !   print*,alt(ic),int(T/ev),"--",sigmas_elas(jj,ind)*1d4,pn_all(jj,ic)*  &
             !        1d-6,int(phi_ch*360/tpi),phi*360/tpi, t_live(ii)
             !end if

             
          end if
          !at the end, recalculate timestep
          delt(ii) =  min(-wce/(mac_sigma*sqrt(vges2)*vc) * log(1-0.005), 1000*wce/(abs(vp1(1,ii))*vc))

          
    end if

       !-----------------------------------------------------------------------
    end do
    !-----------------------------------------------------------------------
    return
    !-----------------------------------------------------------------------
  end subroutine passx
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
end module p__passx


