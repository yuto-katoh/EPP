!-----------------------------------------------------------------------
module p__initcs

  use v__in, only: pi,tpi
  implicit none

  private
  public :: initcs

contains
  real(kind=8) function Get_sigma_ion(Molecule,Orbital,T)
    implicit none  
    real(kind=8), intent(in):: T
    integer, intent(in)  :: Molecule,Orbital
    real(KIND=8) :: U(3,5),N(3,5),Qq(3,5), Bbb(3,5)
    real(KIND=8) :: sigtot, tt,uu,S, a_0,R,nn
    integer :: ii
    !declare constants
    a_0 = 0.529d-10 !bohr-radius
    R = 13.6057d0 !Rydberg-constant
    !factors for the cross section model:
    ! O-Atom
    U(1,:) = (/796.1,84.762,68.5050,69.6520,   0./)
    N(1,:) = (/2.,2.,1.6748,2.3252,   0./)
    Qq(1,:) = (/1.,1.,1.,1.,0./)
    ! N2-Molekül
    U(2,:)  = (/71.13,63.18,44.3,54.91,   0./)
    N(2,:)  = (/2.,2.,4.,2.,   0./)
    Qq(2,:)  = (/1.,1.,1.,1.,   0./)
    ! O2-Molekül
    U(3,:) = (/79.73,90.92,59.89,71.84,84.88/)
    N(3,:) = (/2.,2.,4.,2.,2./)
    Qq(3,:) = (/1.,1.,1.,1.,1./)
    !ionization Energys
    Bbb(1,:) = (/562.878,33.913,16.603,13.618,   1./)
    Bbb(2,:) =  (/41.72,21.  ,17.07,15.58,   1./)
    Bbb(3,:)  = (/46.19,29.82,19.64,19.79,12.07/)

    !get the correct sigma
    ii=Orbital
    tt = T/Bbb(Molecule,ii)
    uu = U(Molecule,ii)/Bbb(Molecule,ii)
    S = 2d0*tpi*(a_0)**2*N(Molecule,ii)*(R/Bbb(Molecule,ii))**2
    Get_sigma_ion = S/(tt+(uu+1))*(Qq(Molecule,ii)*log(tt)/2 &
         *(1-1/tt**2)+(2-Qq(Molecule,ii))*(1-1/tt-log(tt)/(tt+1)))
  end function Get_sigma_ion

  real(KIND=8) function Get_sigma_dif(Energy, phi)
    implicit none
    real(KIND=8), intent(in) :: Energy, phi
    real(KIND=8) :: gamma_e, e_rest_mass, epsilon, gamma, NN
    gamma_e = 0.6
    e_rest_mass = 510998.9  !eV
    epsilon = Energy/e_rest_mass
    gamma = gamma_e * (6.22d-5)/(epsilon*(epsilon+2))
    NN = 0.25* 1/(gamma*(1+gamma))
    Get_sigma_dif = 1/NN*(1/(1+2*gamma-cos(phi))**2)  *sin(phi)

  end function Get_sigma_dif
  
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  subroutine initcs(sigmas_ion, sigmas_elas, sigmas_dif, sigmas_tot)


    real(kind=8), intent(out) ::sigmas_ion(3,10000), sigmas_dif(1000,10000)
    real(kind=8), intent(out) ::sigmas_elas(3,10000),sigmas_tot(3,10000)


    real(kind=8) ::  phs,vpp,aa,factp,vges2,aaa,Tc,tot,phi,scale
    integer      ::  ii,jj,kk,ll,iseed,istat
    real(kind=8) ::  a,b,c

    !calculate sigmas:
    !Energy index  x means  E = exp(x/500) Ev  !
    do ll=1,10000 !!for all energys    do ll=1,10000
       Tc = exp(dble(ll)/500) !!get an exponential scaling ([Tc] = ev)

       !ionization
       do jj=1,3 !!for all gases
          tot = 0d0
          do kk=1,5 !!for all orbitals
             tot = tot + Get_sigma_ion(jj,kk,Tc)
          end do
          sigmas_ion(jj,ll)=tot
       end do

       !elastic scattering (way too simple model)
       Tc = exp(dble(ll)/500) !!get an exponential scaling ([Tc] = ev) 
       sigmas_elas(1,ll) = dble(10)**(-18.3641) * Tc**(-0.63593)
       sigmas_elas(2,ll) = dble(10)**(-18.05689) * Tc**(-0.63211)
       sigmas_elas(3,ll) = dble(10)**(-18.81805) * Tc**(-0.7278020)

       !differential cross section
       do jj = 1,1000
          phi = (dble(jj)-1)/1000*pi
          sigmas_dif(jj,ll) = Get_sigma_dif(Tc,phi)
       end do
       scale = sum(sigmas_dif(:,ll))
       sigmas_dif(:,ll) = sigmas_dif(:,ll)/scale
       do kk = 2,1000
          sigmas_dif(kk,ll) =  sigmas_dif(kk-1,ll) +  sigmas_dif(kk,ll)
       end do


    end do

    sigmas_tot =  sigmas_ion +  sigmas_elas

    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
 


  end subroutine initcs

end module p__initcs
