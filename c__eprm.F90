!-----------------------------------------------------------------------
module c__eprm
  implicit none
!-----------------------------------------------------------------------
  integer(kind=4) :: istep = 200
!  integer(kind=4),parameter :: istep = 10000000
!-----------------------------------------------------------------------

    
  integer(kind=4),parameter :: npa1 = 10000
  
  integer(kind=4),parameter :: ialt = 400
  integer(kind=4),parameter :: resf = 1000 !(resolution of atmosphere 1/km)
  integer(kind=4),parameter :: resfout = 5 !(resolution of ionization profile 1/km)

  real(kind=8),parameter :: live_max = 200.0  !max lifetime of the particles
  
  !real(kind=8),parameter :: live_res = 0.05   !time resolution of output file
  !set these equal for integrated profile!
  real(kind=8),parameter :: live_res = live_max
  
  integer(kind=4),parameter :: numtout = int(live_max/live_res+0.4999999)
  
  integer(kind=4),parameter :: mfopt = 1 !factor for mirror force; 0:OFF 1:ON
  integer(kind=4),parameter :: naopt = 1 !factor for neutral atmosphere; 0:OFF 1:ON

  integer(kind=4),parameter :: lx = resf*ialt , lxa1 = lx
  integer(kind=4),parameter :: lpn1 = 256

  !  integer(kind=4),parameter :: npa1 = lpn1*lx
  integer(kind=4),parameter :: nsp = 2 , nqa = 2
  integer(kind=4),parameter :: lx1 = 2 , lx2 = lx-2

!-----------------------------------------------------------------------
end module c__eprm
