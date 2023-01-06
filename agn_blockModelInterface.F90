!*******************************************************************************

!  Module:     agn_blockModelInterface

!  Description: All routines here map subgrid models onto a given block.
!               All routines take a parameter list for the models.
!               
!               These routines allow us to abstract away the physical
!                 models and geometry of feedback regions. We use these
!                 routines to choose appropriate physics models and 
!                 geometries based on the parameter list.


module agn_blockModelInterface

  use agn_data
  use agn_geometryTools
  use agn_miscTools
  use agn_physicsModels
  use ieee_arithmetic

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

contains

!------------------------------------------------------------------------------
! samples the gas around a point to estimate accretion rate
subroutine getAccretionRate(lb, parmList, dt, mgdot, Rho_surr, alpha, pmass)

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, & 
                             Grid_getBlkIndexLimits, Grid_getBlkPhysicalSize, &
                             Grid_releaseBlkPtr, Grid_getBlkCenterCoords, &
                             Grid_getCellCoords

  use Multispecies_interface, ONLY : Multispecies_getSum
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Grid_data, ONLY : gr_smalle
 
  implicit none

#include "Multispecies.h"

  real, parameter :: Tq = 5.e5

  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(IN) :: parmList
  real, intent(IN) :: dt, pmass
  real, intent(OUT) :: Rho_surr, mgdot, alpha
 
  real :: pos(3), mass
  real :: dV, dx(3), size(3), coord(3), blkCenter(3), zoneCoord(3), dist2, dist(3)
  real :: faceCoord(3), faceDist(3)
  real, dimension(GRID_KHI_GC):: zoneZ
  real, dimension(GRID_JHI_GC):: zoneY
  real, dimension(GRID_IHI_GC):: zoneX
  integer :: i, j, k
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  logical :: isInside, isInInnerShell, isInOuterShell
  real :: phi, theta
  real :: tag, Cs_surr, Vol_surr, Mdot_edd, cs
  real :: r_in, r_out, s_in, s_out, Vol_in, Vol_out
  real, pointer :: solnData(:,:,:,:)
  real :: mh_surr, mh, nh, m_cold
  real, save :: mu
  logical :: firstCall = .true.

  if (firstCall) then
    call PhysicalConstants_get("proton mass", mu)
    firstCall =.false.
  end if

  alpha  = parmList(PARM_LIST_ALPHA)
  pos(:) = parmList(PARM_LIST_CENTERX:PARM_LIST_CENTERZ)
  mass   = parmList(PARM_LIST_MASS)
  call getAngleFromTag(parmList(PARM_LIST_TAG), phi, theta)

  mgdot = 0.
  
  call Grid_getBlkPhysicalSize(lb ,size)    ! dimensions of block
  call Grid_getBlkCenterCoords(lb, blkCenter)
  coord = blkCenter - 0.5*size                 ! block edge coordinates
  call Grid_getDeltas(lb, dx)
  dV    = product(dx(:))                  ! cell volume

  r_in  = 2.*dx(1)
  r_out = 4.*dx(1)

  call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)
  call Grid_getBlkPtr(lb, solnData, CENTER)

  Rho_surr = 0.
  Cs_surr  = 0.
  Vol_surr = 0.
  mh_surr  = 0.
  s_in     = 0.
  s_out    = 0.
  Vol_in   = 0.
  Vol_out  = 0.
  m_cold   = 0.

  call Grid_getCellCoords(IAXIS, lb, CENTER, .true., zoneX, GRID_IHI_GC)
  call Grid_getCellCoords(JAXIS, lb, CENTER, .true., zoneY, GRID_JHI_GC)
  call Grid_getCellCoords(KAXIS, lb, CENTER, .true., zoneZ, GRID_KHI_GC)

  ! average accretion rate over surrounding zones
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
  do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
  do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
    zoneCoord(1) = zoneX(i)
    zoneCoord(2) = zoneY(j)
    zoneCoord(3) = zoneZ(k)
  
    ! build a cylinder
    isInside = isInRegion(REGION_CUBE, pos, agn_rAccrete*dx(1), &
               agn_rAccrete*dx(1), &
               phi, theta, zoneCoord, agn_Lbox)
     
    if (isInside) then
      if(trim(agn_AccretionMode) == "cold") then ! sum all gas with T < Tq = 5.e5K
        if(solnData(TEMP_VAR,i,j,k) <= Tq) then
          m_cold = m_cold + solnData(DENS_VAR,i,j,k)*dV
        endif
      else
        if(solnData(EINT_VAR,i,j,k) < gr_smalle) &   
          solnData(EINT_VAR,i,j,k) = gr_smalle 
        Rho_surr = Rho_surr + solnData(DENS_VAR,i,j,k)
        cs       = sqrt( solnData(GAMC_VAR,i,j,k) &
                         * (solnData(GAMC_VAR,i,j,k)-1) &
                         * solnData(EINT_VAR,i,j,k) )
        if(IEEE_IS_NAN(cs)) then 
          print*, 'gamc, eint, cs = ', solnData(GAMC_VAR,i,j,k), solnData(EINT_VAR,i,j,k), cs
        endif
        Cs_surr  = Cs_surr + cs
        !! TODO : test getsum
        !call Multispecies_getSum(A, mh)
        !mh_surr  = mh_surr + mh
        Vol_surr = Vol_surr + 1 ! cell size is uniform in one block 
      endif
    end if

    ! Computing entropies in shells surrouding the BH
    isInInnerShell = isInRegion(REGION_SHELL, pos, r_in, 0., &
                     phi, theta, zoneCoord, agn_Lbox)
    isInOuterShell = isInRegion(REGION_SHELL, pos, r_out, r_in, &
                     phi, theta, zoneCoord, agn_Lbox)
    if (isInInnerShell) then
      s_in = s_in + solnData(TEMP_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)**(2./3.)
      Vol_in = Vol_in + 1
    end if
    if (isInOuterShell) then
      s_out = s_out + solnData(TEMP_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)**(2./3.)
      Vol_out = Vol_out + 1
    end if

  end do
  end do
  end do ! search zones
  call Grid_releaseBlkPtr(lb, solnData)

  Rho_surr = Rho_surr / Vol_surr
  Cs_surr  = Cs_surr  / Vol_surr
  !mh_surr  = mh_surr / Vol_surr * mu
  nh       = Rho_surr / mu
  s_in     = s_in / Vol_in
  s_out    = s_out / Vol_out

  if (trim(agn_AccretionMode) == "pope") then
    call popeAccretion(mass, Rho_surr, Cs_surr, mgdot)
  else if (trim(agn_AccretionMode) == "alpha") then
    call alphaAccretion(mass, Rho_surr, Cs_surr, mgdot)
  else if (trim(agn_AccretionMode) == "beta") then
    call betaAccretion(mass, Rho_surr, Cs_surr, nh, mgdot)
  else if (trim(agn_AccretionMode) == "entropy") then
    call entropyAccretion(mass, Rho_surr, Cs_surr, r_in, r_out, s_in, s_out, alpha, mgdot)
  else if (trim(agn_AccretionMode) == "cold") then
    mgdot = (m_cold + pmass) / ACCRETION_TIMESCALE
    !!! DEBUG
    !if(abs(mgdot) > HUGE(1.) .or. mgdot .lt. 0.) &   
    !  print*, 'mgdot, m_cold, pmass = ', mgdot, m_cold, pmass
  else
    call bondiAccretion(mass, Rho_surr, Cs_surr, mgdot)
  end if

  ! Apply Eddington upper limit
  Mdot_edd = 4.*agn_pi*agn_G*mass*agn_mp/agn_sigmaT/agn_c/agn_Efff
  if (trim(agn_AccretionMode) /= "cold" .and. mgdot > Mdot_edd) mgdot = Mdot_edd

  if(trim(agn_feedbackMode) == "episodic jet" .or. &
     trim(agn_feedbackMode) == "cosmic ray") mgdot = 0.

  return
end subroutine getAccretionRate


!------------------------------------------------------------------------------
! get the resolution near the AGN
subroutine getLocalResolution(lb, parmList, localDx)

  use Grid_interface, ONLY : Grid_getDeltas

  implicit none

  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(IN) :: parmList
  real, intent(OUT) :: localDx

  real :: dx(3)

  call Grid_getDeltas(lb, dx)

  localDx = dx(1)

  return
end subroutine getLocalREsolution


!------------------------------------------------------------------------------
! computes a radius so that we don't remove more than 10% of the gas
! we will also enforce a minimum of 2 zones
subroutine getDepletionRadius(lb, parmList, dt, localDensity, rDeplete)

  implicit none

  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(IN) :: parmList
  real, intent(IN) :: dt, localDensity
  real, intent(OUT) :: rDeplete

  rDeplete = parmList(PARM_LIST_MGDOT) * dt
  rDeplete = rDeplete / localDensity / agn_43pi / agn_depleteFrac
  rDeplete = rDeplete**(1/3.)
 
  if(agn_rDeplete == 0) then
    rDeplete = agn_rDeplete
  else 
    rDeplete = max(rDeplete, agn_rDeplete*parmList(PARM_LIST_LOCALDX))
  endif

  !rDeplete = 2.0*parmList(PARM_LIST_LOCALDX)
  !rDeplete = parmList(PARM_LIST_LOCALDX)
  !rDeplete = 2.0*agn_kpc2cm
  !rDeplete = max(rDeplete, 2.0*agn_kpc2cm)
  !rDeplete = max(rDeplete, 2.0*parmList(PARM_LIST_LOCALDX))
  !rDeplete = max(rDeplete, 2.0*agn_kpc2cm, parmList(PARM_LIST_LOCALDX))

  return
end subroutine getDepletionRadius


!------------------------------------------------------------------------------
! removes gas around the BH (due to accretion)
! Note: This function is not parallel and only removes from neighboring cells
! Currently using removeAllGas instead
subroutine removeGas(lb, parmList, dt)

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, & 
                             Grid_getBlkIndexLimits, Grid_getBlkPhysicalSize, &
                             Grid_releaseBlkPtr, Grid_getBlkCenterCoords
  use Grid_data, ONLY : gr_smallrho

  use Eos_interface, ONLY : Eos_wrapped

  implicit none

  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(IN) :: parmList
  real, intent(IN) :: dt

  real :: pos(3), dvol, mgdot, massInCell,massToRemove, dx(3), size(3), coord(3)
  real :: blkCenter(3)
  integer :: iPos, jPos, kPos, i, j, k
  integer :: blkLimits(2,MDIM), blkLimitsGC(2,MDIM)
  real, pointer :: solnData(:,:,:,:)
  logical, save :: firstCall = .true.

  if (firstCall) then
    firstCall = .false.
  end if

  pos(:) = parmList(PARM_LIST_CENTERX:PARM_LIST_CENTERZ)
  mgdot  = parmList(PARM_LIST_MGDOT)

  massToRemove = mgdot*dt/27. ! shared amongst nearest cells

  call Grid_getBlkPhysicalSize(lb ,size)    ! dimensions of block
  call Grid_getBlkCenterCoords(lb, blkCenter)
  coord = blkCenter - 0.5*size                 ! block edge coordinates
  call Grid_getDeltas(lb, dx)
  dvol = product(dx(1:3))
  
  call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)

  iPos = blkLimits(1,1) + int((pos(1)-coord(1))/dx(1))
  jPos = blkLimits(1,2) + int((pos(2)-coord(2))/dx(2))
  kPos = blkLimits(1,3) + int((pos(3)-coord(3))/dx(3))

  call Grid_getBlkPtr(lb, solnData, CENTER)

  do k = kPos-1, kPos+1
  do j = jPos-1, jPos+1
  do i = iPos-1, iPos+1
    massInCell = solnData(DENS_VAR,i,j,k)*dvol
    solnData(DENS_VAR,i,j,k) = max((massInCell-massToRemove)/dvol, gr_smallrho)
  end do
  end do
  end do ! scan zones

  call Eos_wrapped(MODE_DENS_EI, blkLimits, lb)

  call Grid_releaseBlkPtr(lb, solnData)

  return
end subroutine removeGas

!------------------------------------------------------------------------------
! Compare accretion rate to Eddington limit. Uses this comparison and 
!   user parameters to determine appropriate feedback mode

function getFeedbackMode(lb, parmList, dt, simTime)

  implicit none
  
  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(IN) :: parmList
  real, intent(IN) :: dt, simTime
  real :: oldBHMass, feedbackMode, getFeedbackMode, jetClock

  real :: mgdot, mbh, Mdot_edd, delM

  mgdot = parmList(PARM_LIST_MGDOT)
  mbh   = parmList(PARM_LIST_MASS)
  oldBHMass = parmList(PARM_LIST_OLDMASS)

  Mdot_edd = 4.*agn_pi*agn_G*mbh*agn_mp/agn_sigmaT/agn_c/agn_Efff

  if (trim(agn_AccretionMode) == "cold" .and. &
      mgdot == 0.) then ! this is possible in the cold gas feedback mode
    getFeedbackMode = FEEDBACK_MODE_NONE
    return
  endif

  if ( mgdot/Mdot_edd > agn_threshold .or. & 
       trim(agn_feedbackMode) == "pure quasar")  then
    feedbackMode = FEEDBACK_MODE_QUASAR
  else
    if (trim(agn_feedbackMode) == "jet") then
      feedbackMode = FEEDBACK_MODE_JET

    else if (trim(agn_feedbackMode) == "episodic jet") then
      feedbackMode = FEEDBACK_MODE_NONE
      if(simTime >= agn_jetStartTime) then
        jetClock = mod(simTime - agn_jetStartTime, agn_jetCycle)
        if (jetClock <= agn_jetDuration) then
          feedbackMode = FEEDBACK_MODE_JET
        endif
      endif

    else if (trim(agn_feedbackMode) == "cosmic ray") then
      feedbackMode = FEEDBACK_MODE_NONE
      if(simTime >= agn_jetStartTime) then
        jetClock = mod(simTime - agn_jetStartTime, agn_jetCycle)
        if (jetClock <= agn_jetDuration) then
          feedbackMode = FEEDBACK_MODE_COSMICRAY
        endif
      endif

    else
      delM = (mbh - oldBHMass)/oldBHMass
      if (delM >= agn_minDelM) then
        feedbackMode = FEEDBACK_MODE_BUBBLE
      else
        feedbackMode = FEEDBACK_MODE_NONE
      end if
    end if
  end if

  getFeedbackMode = feedbackMode

  return
end function getFeedbackMode


!-----------------------------------------------------------------------------
! for a given feedback model, chooses an appropriate feedback region 
!   (i.e. Bubbles go in spheres, jets go in cylinders)

function getFeedbackRegion(lb, parmList)

  implicit none

  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(IN) :: parmList
  real :: getFeedbackRegion, regionShape

  if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_JET) then
    regionShape = REGION_CYLINDER
  else if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_COSMICRAY) then
    regionShape = REGION_CYLINDER
  else if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_QUASAR) then
    regionShape = REGION_SPHERE
  else if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_BUBBLE) then
    regionShape = REGION_SPHERE
  else
    regionShape = REGION_POINT
  end if

  getFeedbackRegion = regionShape

  return
end function getFeedbackRegion


!-----------------------------------------------------------------------------
! Returns the center of the feedback region.
!   Jet centers are on the central BH.
!   Bubbles get randomly placed within a sphere twice the bubble diameter.
subroutine getFeedbackCenter(lb, parmList, centerBH, center)

  use healpix_types
  use rngmod
  use Driver_data, ONLY : dr_globalMe
 
  implicit none
 
  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(IN) :: parmList
  real, intent(IN) :: centerBH(3)
  real, intent(OUT) :: center(3)
  real :: dist, flip
  integer :: iDir
  type(planck_rng), save :: pos_seed
  logical, save :: firstCall = .true.

  if (firstCall) then
    call rand_init(pos_seed, dr_globalMe+1, dr_globalMe+2, dr_globalMe+3, dr_globalMe+4)

    firstCall = .false.
  endif

  if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_BUBBLE .and. &
      trim(agn_feedbackMode) /= "fixbub" ) then
    do iDir = 1, 3
      dist = parmList(PARM_LIST_REJ) * rand_uni(pos_seed)
      flip = rand_uni(pos_seed)
      if (flip <= 0.5) then
        flip = -1.
      else
        flip = 1.
      end if

      center(iDir) = centerBH(iDir) + dist * flip

      ! assume periodic boundary conditions
      if (center(iDir) >= agn_xmax) then
        center(iDir) = center(iDir) - agn_Lbox
      else if (center(iDir) < agn_xmin) then
        center(iDir) = center(iDir) + agn_Lbox
      end if
    end do ! iDir

  else
    center(:) = centerBH(:)
  end if

  return
end subroutine getFeedbackCenter


!------------------------------------------------------------------------------
! Returns the angle of the feedback region.
!   For jets, this is the angle assigned when the BH was created.
!     If agn_jetPrecessAngle/=0., assign a random theta from the BH jet axis.
!   For bubbles, which have randomly-assigned positions, 
!     the bubble axis will point back to the central BH.
subroutine getFeedbackAngle(lb, parmList, centerBH, center, phi, theta)

  use healpix_types
  use rngmod
  use Driver_data, ONLY: dr_globalMe, dr_simTime

  implicit none

  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(IN) :: parmList
  real, intent(IN) :: centerBH(3), center(3)
  real, intent(OUT) :: phi, theta
  integer :: tag, iDir
  real :: phiBH, thetaBH, dist(3), r, theta_p, phi_p
  type(planck_rng), save :: angle_seed
  logical, save :: firstCall = .true.

  if(firstCall) then
    call rand_init(angle_seed, dr_globalMe+1, dr_globalMe+2, dr_globalMe+3, dr_globalMe+4)
    firstCall = .false.
  endif

  call getAngleFromTag(parmList(PARM_LIST_TAG), phiBH, thetaBH)

  if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_BUBBLE .and. &
      trim(agn_feedbackMode) /= "fixbub") then

    do iDir = 1, 3
      dist(iDir) = center(iDir) - centerBH(iDir)
      if (dist(iDir) >= agn_xmax) then
        dist(iDir) = dist(iDir) - agn_Lbox
      else if (dist(iDir) < agn_xmin) then
        dist(iDir) = dist(iDir) + agn_Lbox
      end if
    end do

    r = sqrt( sum(dist(:)**2.) )
    theta = acos(dist(3)/r)
  
    r = sqrt( sum(dist(1:2)**2.) )
    if (dist(2) >= 0.) then 
      phi = acos(dist(1)/r)
    else
      phi = agn_pi + (dist(1)/r)
    end if

    if (phi /= phi) phi = 0.

  else
!    theta_p = rand_uni(angle_seed) * agn_jetPrecessAngle / 180. * agn_pi
    theta_p = agn_jetPrecessAngle / 180. * agn_pi
    if(agn_jetPrecessAngle == 0.) then
      phi_p = 0.
    else
!      phi_p   = rand_uni(angle_seed) * 2. * agn_pi
      phi_p   = mod(dr_simTime/PRECESSION_PERIOD, 2.*agn_pi) + 0.5*agn_pi
    endif
!    theta   = thetaBH + theta_p
!    phi     = phiBH + phi_p
    theta   = theta_p
    phi     = phi_p
  end if
  
  return
end subroutine getFeedbackAngle


!------------------------------------------------------------------------------
! returns the radius and height of injection region
subroutine getFeedbackScale(lb, parmList, dt, localDensity, rej, hej)

  use Grid_interface, ONLY : Grid_getDeltas

  implicit none

  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(IN) :: parmList
  real, intent(IN) :: dt, localDensity
  real, intent(OUT) :: rej, hej
  real :: Ebub, mass, oldMass, dx(3)

  if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_JET) then
    rej = agn_Rej
    hej = agn_Hej

  else if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_COSMICRAY) then
    rej = agn_Rej
    hej = agn_Hej

  else if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_QUASAR) then
    call Grid_getDeltas(lb, dx)
    rej = 8.0*dx(1)
    hej = 8.0*dx(1)

  else if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_BUBBLE) then

    if ( trim(agn_feedbackMode) == "bubble" ) then
      mass    = parmList(PARM_LIST_MASS)
      oldMass = parmList(PARM_LIST_OLDMASS)
      call getFeedbackEnergy(parmList(PARM_LIST_MGDOT), mass, oldMass, dt, &
                            FEEDBACK_MODE_BUBBLE, Ebub)
      Ebub = Ebub * dt
      rej = agn_Rej * (Ebub/agn_Ebub0/localDensity*agn_rhoBub0)**(1./5.)

! Fixed radius that scaled with Mvir (Battaglia 2010)
!      rej = 70.8 * agn_kpc2cm 
! Fixed small radius to compare with jets
!      rej = 3.0 * agn_kpc2cm
      
      hej = rej
    else
      rej = agn_Rej
      hej = agn_Hej
    end if
  
  else
    rej = 0.
    hej = 0.
  end if

  return
end subroutine getFeedbackScale


!------------------------------------------------------------------------------
! this routine returns the volume of a cube which completely
!   encloses the requested region. This is used to aid parallel processing.

function getFeedbackVol(lb, parmList)

  implicit none
  
  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(IN) :: parmList
  real :: getFeedbackVol
  real :: phi, theta, regionShape, rej, hej, tag

  tag = parmList(PARM_LIST_TAG)
  call getAngleFromTag(tag, phi, theta)
  regionShape = parmList(PARM_LIST_REGION)
  rej = parmList(PARM_LIST_REJ)
  hej = parmList(PARM_LIST_HEJ)

  getfeedbackVol = getVolume(regionShape, hej, rej, phi, theta)
  
  return
end function getFeedbackVol


!-----------------------------------------------------------------------------
! tests if a given feedback region will fit entirely on the local processor
subroutine testFit(parmList, willFit)

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getListOfBlocks, &
                             Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords

  implicit none

  real, dimension(MAX_PARM), intent(IN) :: parmList
  logical, intent(OUT) :: willFit
  integer :: b, lb, blkCount, blkList(MAXBLOCKS)
  real :: regionShape, tag, pos(3), regionRadius, regionVolume
  real :: volRemaining, volIntersect, phi, theta, rej, hej, rdepl
  real :: size(3), blkRadius, blkCenter(3)

  call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  pos(:)      = parmList(PARM_LIST_CENTERX:PARM_LIST_CENTERZ)
  tag = parmList(PARM_LIST_TAG)
  call getAngleFromTag(tag, phi, theta)
  regionShape = parmList(PARM_LIST_REGION)
  rej = parmList(PARM_LIST_REJ)
  hej = parmList(PARM_LIST_HEJ)
  rdepl = parmList(PARM_LIST_RDEPLETE)
  if(rdepl .gt. rej .or. rdepl .gt. hej) then
    rej = rdepl
    hej = rdepl
    regionShape = REGION_CUBE
  endif
  regionRadius = getRadius(regionShape, hej, rej, phi, theta)
  regionVolume = getVolume(regionShape, hej, rej, phi, theta)
  volRemaining = parmList(PARM_LIST_VOL)

  if (regionShape == REGION_POINT) then
    willFit = .true.
  else 
    do b = 1, blkCount
      lb = blkList(b)
      call Grid_getBlkPhysicalSize(lb, size)      ! dimensions of block
      call Grid_getBlkCenterCoords(lb, blkCenter) ! block center
      blkRadius = 0.5*size(1)                     ! half-width of block
    
      volIntersect = getVolIntersect(pos, regionRadius, blkCenter, & 
                                     blkRadius, agn_Lbox)

      volRemaining = volRemaining - volIntersect
    end do

    willFit = (abs(volRemaining/regionVolume) <= VOL_TOL)
  end if

  return
end subroutine testFit

!-----------------------------------------------------------------------------
! removes gas around all AGN
subroutine removeAllGas(lb, parmList, doneApplying)

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getCellCoords, &
                             Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr, &
                             Grid_getBlkIndexLimits
  use Grid_data, ONLY: gr_smallrho
 
  use Driver_interface, ONLY : Driver_getDt 
  use Grid_data, ONLY : gr_smallrho
  use Eos_interface, ONLY : Eos_wrapped
  use agn_miscTools
  use agn_injectTracersModule
  use Cool_data
  use Driver_data, ONLY : dr_simTime, dr_globalNumProcs, dr_globalMe
  use Particles_data, ONLY : pt_numLocal

  implicit none

  real, parameter :: Tq = 5.e5

  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(INOUT) :: parmList
  logical, intent(OUT) :: doneApplying

  real, pointer :: solnData(:,:,:,:)
  real :: zoneCoord(3), coord(3), size(3), dx(3), dV, blkCenter(3), blkRadius
  real :: displ(3), del(3)
  integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
  real, dimension(GRID_KHI_GC) :: zoneZ
  real, dimension(GRID_JHI_GC) :: zoneY
  real, dimension(GRID_IHI_GC) :: zoneX      
  logical :: isInside, isInsideMag, magResolved
  integer :: i, j, k, ii, jj, kk, neq, ipar
  integer, save :: nstep, oldNstep

  real :: Eej, Mej, Pej, vxej, vyej, vzej, dt, magEnergy, magStrength
  real :: mass, oldMass, time

  real :: pos(3), mgdot, mode, regionShape,theta,phi, regionRadius, regionVolume
  real :: volIntersect, volRemaining, volFracRemaining, rej, hej, massToRemove
  real :: massInCell, eid, dedt, tcool, q, rpar, tag, pv
  logical, save :: firstCall = .true.
  integer, save :: ptagmult

  pos(:)      = parmList(PARM_LIST_CENTERX:PARM_LIST_CENTERZ)
  mgdot       = parmList(PARM_LIST_MGDOT)
  regionShape = REGION_CUBE
  rej         = parmList(PARM_LIST_RDEPLETE)
  hej         = parmList(PARM_LIST_RDEPLETE)

  if (firstCall) then
    ptagmult = 10**(int(alog10(real(dr_globalNumProcs)))+1)
    firstCall = .false.
  endif

  doneApplying = .false.

  time = dr_simTime
  call Grid_getBlkPhysicalSize(lb,size)    ! dimensions of block
  call Grid_getBlkCenterCoords(lb, blkCenter)
  coord = blkCenter - 0.5*size                 ! block edge coordinates
  call Grid_getDeltas(lb,del)
  dx    = del(1)                           ! cell dimensions
  dV    = product(del(:))                  ! cell volume
  blkRadius = 0.5*size(1)                  ! half-width of block

  call Driver_getDt(dt)
  massToRemove = mgdot*dt/(agn_43pi*rej**3)
  neq = 0
  ipar = 0
  rpar = 0.

  ! check if injection region intersects current block
  regionRadius = getRadius(regionShape, hej, rej, phi, theta)
  regionVolume = getVolume(regionShape, hej, rej, phi, theta)
  volIntersect =getVolIntersect(pos,regionRadius,blkCenter,blkRadius,agn_Lbox)
  volRemaining = parmList(PARM_LIST_VOL) - volIntersect
  volFracRemaining = volRemaining/regionVolume

  doneApplying = (parmList(PARM_LIST_VOL)/regionVolume <= VOL_TOL)
  if (.not. doneApplying .and. volIntersect > 0.) then

    call Grid_getBlkPtr(lb, solnData)

    ! unpack zones
    call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)

    call Grid_getCellCoords(KAXIS, lb, CENTER, .true., zoneZ, GRID_KHI_GC)
    call Grid_getCellCoords(JAXIS, lb, CENTER, .true., zoneY, GRID_JHI_GC)
    call Grid_getCellCoords(IAXIS, lb, CENTER, .true., zoneX, GRID_IHI_GC)

    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
    do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
      zoneCoord(1) = zoneX(i)
      zoneCoord(2) = zoneY(j)
      zoneCoord(3) = zoneZ(k)

      isInside = isInRegion(regionShape, pos, hej, rej, &
                            phi, theta, zoneCoord, agn_Lbox)

      if (parmList(PARM_LIST_MODE) /= FEEDBACK_MODE_NONE .and. isInside) then
        ! Remove accreted hot gas
        massInCell = solnData(DENS_VAR,i,j,k)*dV
        solnData(DENS_VAR,i,j,k) = max((massInCell-massToRemove)/dV, &
                                  (1.0-agn_depleteFrac)*massInCell/dV)
      end if ! we are in the injection region
    end do
    end do
    end do ! loop over zones

    parmList(PARM_LIST_VOL) = volRemaining
    doneApplying = (volFracRemaining <= VOL_TOL)

    call Grid_releaseBlkPtr(lb, solnData)
  
    call Eos_wrapped(MODE_DENS_EI, blkLimits, lb)

  end if ! sum in this block

  return
end subroutine removeAllGas


!-----------------------------------------------------------------------------
! sums the mass in the given injection region 
subroutine getTotalMass(lb, parmList, doneApplying)

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getCellCoords, &
                             Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr, &
                             Grid_getBlkIndexLimits
  
  use agn_miscTools

  implicit none

  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(INOUT) :: parmList
  logical, intent(OUT) :: doneApplying

  real, pointer :: solnData(:,:,:,:)
  real :: zoneCoord(3), coord(3), size(3), dx(3), dV, blkCenter(3), blkRadius
  real :: displ(3), del(3)
  integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
  real, dimension(GRID_KHI_GC) :: zoneZ
  real, dimension(GRID_JHI_GC) :: zoneY
  real, dimension(GRID_IHI_GC) :: zoneX      
  logical :: isInside, isInsideMag, magResolved
  integer :: i, j, k, ii, jj, kk
  integer, save :: nstep, oldNstep

  real :: Eej, Mej, Pej, vxej, vyej, vzej, dt, magEnergy, magStrength
  real :: mass, oldMass, time

  real :: pos(3), mgdot, mode, regionShape,theta,phi, regionRadius, regionVolume
  real :: volIntersect, volRemaining, volFracRemaining, rej, hej, totalMass
  real :: tag
  logical, save :: firstCall = .true.

  pos(:)      = parmList(PARM_LIST_CENTERX:PARM_LIST_CENTERZ)
  mgdot       = parmList(PARM_LIST_MGDOT)
  mass        = parmList(PARM_LIST_MASS)
  oldMass     = parmList(PARM_LIST_OLDMASS)
  regionShape = parmList(PARM_LIST_REGION)
  mode        = parmList(PARM_LIST_MODE)
  rej         = parmList(PARM_LIST_REJ)
  hej         = parmList(PARM_LIST_HEJ)
  tag         = int(parmList(PARM_LIST_TAG))
  call getAngleFromTag(tag, phi, theta)

  doneApplying = .false.

  call Grid_getBlkPhysicalSize(lb,size)    ! dimensions of block
  call Grid_getBlkCenterCoords(lb, blkCenter)
  coord = blkCenter - 0.5*size                 ! block edge coordinates
  call Grid_getDeltas(lb,del)
  dx    = del(1)                           ! cell dimensions
  dV    = product(del(:))                  ! cell volume
  blkRadius = 0.5*size(1)                  ! half-width of block

  totalMass = 0.

  if (regionShape == REGION_POINT) then
    doneApplying = .true.
  else
    ! check if injection region intersects current block
    regionRadius = getRadius(regionShape, hej, rej, phi, theta)
    regionVolume = getVolume(regionShape, hej, rej, phi, theta)
    volIntersect =getVolIntersect(pos,regionRadius,blkCenter,blkRadius,agn_Lbox)
    volRemaining = parmList(PARM_LIST_VOL) - volIntersect
    volFracRemaining = volRemaining/regionVolume

    doneApplying = (parmList(PARM_LIST_VOL)/regionVolume <= VOL_TOL)
    if (.not. doneApplying .and. volIntersect > 0.) then

      call Grid_getBlkPtr(lb, solnData)

      ! unpack zones
      call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)

      call Grid_getCellCoords(KAXIS, lb, CENTER, .true., zoneZ, GRID_KHI_GC)
      call Grid_getCellCoords(JAXIS, lb, CENTER, .true., zoneY, GRID_JHI_GC)
      call Grid_getCellCoords(IAXIS, lb, CENTER, .true., zoneX, GRID_IHI_GC)

      do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
      do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        zoneCoord(1) = zoneX(i)
        zoneCoord(2) = zoneY(j)
        zoneCoord(3) = zoneZ(k)

        isInside = isInRegion(regionShape, pos, hej, rej, &
                              phi, theta, zoneCoord, agn_Lbox)

        if (isInside) then
          totalMass = totalMass + solnData(DENS_VAR,i,j,k)*dV
        end if ! we are in the injection region
      end do
      end do
      end do ! loop over zones

      parmList(PARM_LIST_REGMASS) = parmList(PARM_LIST_REGMASS) + totalMass
      parmList(PARM_LIST_VOL) = volRemaining
      doneApplying = (volFracRemaining <= VOL_TOL)

      call Grid_releaseBlkPtr(lb, solnData)
    end if ! sum in this block
  end if ! point versus extended sources

  return
end subroutine getTotalMass


!------------------------------------------------------------------------------
! unpacks blocks into zones and applies feedback to gas
subroutine applyFeedback(lb, parmList, doneApplying)

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getCellCoords, &
                             Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr, &
                             Grid_getBlkIndexLimits

  use Logfile_interface, ONLY : Logfile_stamp
  use agn_injectTracersModule
  use agn_physicsModels, ONLY : getFeedbackEnergy
  use Driver_interface, ONLY : Driver_getDt, Driver_getNStep, Driver_getSimTime
  use Driver_data, ONLY : dr_globalMe, dr_simTime
#ifdef CRAY_MSCALAR 
  use eos_gammaCRData
#endif

  implicit none
  
  integer, intent(IN) :: lb
  real, dimension(MAX_PARM), intent(INOUT) :: parmList
  logical, intent(OUT) :: doneApplying

  real, pointer :: solnData(:,:,:,:)
  real :: zoneCoord(3), coord(3), size(3), dx(3), dV, blkCenter(3), blkRadius
  real :: displ(3), del(3)
  real, dimension(GRID_KHI_GC) :: zoneZ
  real, dimension(GRID_JHI_GC) :: zoneY
  real, dimension(GRID_IHI_GC) :: zoneX
  logical :: isInside, isInsideMag, magResolved
  integer :: partBlk, blkList(1)
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: i, j, k, ipos(3)
  integer, save :: nstep, oldNstep = huge(1)

  real :: Eej, Edot, Mej, Pej, vxej, vyej, vzej, Dej, Ecrej
  real :: dt, magEnergy, magStrength, egas, pgas, ecr, pcr
  real :: mass, oldMass, totalMass, time
  real :: magRegionShape, magRegionRadius, magRegionScale

  real :: pos(3), mgdot, mode, regionShape,theta,phi, regionRadius, regionVolume
  real :: tag, volIntersect, volRemaining, volFracRemaining, rej, hej
  real :: Bx, By, Bz, Bxej, Byej, Bzej, mag
 
  real :: totalE = 0., ener   , totalVol = 0., totalB = 0., totalEej

  integer :: nsubzones = 4, ii, jj, kk
  real :: subCoord(3), subDispl(3), subdx(3), subdV, cellLeftEdge(3)

  pos(:)      = parmList(PARM_LIST_CENTERX:PARM_LIST_CENTERZ)
  mgdot       = parmList(PARM_LIST_MGDOT) 
  mass        = parmList(PARM_LIST_MASS) 
  oldMass     = parmList(PARM_LIST_OLDMASS) 
  totalMass   = parmList(PARM_LIST_REGMASS) 
  regionShape = parmList(PARM_LIST_REGION)
  mode        = parmList(PARM_LIST_MODE)
  rej         = parmList(PARM_LIST_REJ)
  hej         = parmList(PARM_LIST_HEJ)
  tag         = parmList(PARM_LIST_TAG)
  call getAngleFromTag(tag, phi, theta)

  call Driver_getDt(dt)

  doneApplying = .false.

  call Grid_getBlkPhysicalSize(lb,size)    ! dimensions of block
  call Grid_getBlkCenterCoords(lb, blkCenter)
  coord = blkCenter - 0.5*size                 ! block edge coordinates
  call Grid_getDeltas(lb,del)
  dx(:)    = del(:)                           ! cell dimensions
  dV    = product(del(:))                  ! cell volume
  blkRadius = 0.5*size(1)                  ! half-width of block

  subdx(:) = del(:)/nsubzones
  subdV    = product(subdx(:))

  ! magnetic injection parameters
  magRegionShape = REGION_SPHERE
  magRegionRadius = min(hej, rej)
  magRegionScale = magRegionRadius/2.


!------------------------------------------------------------------------
  ! apply point sources here
  if (regionShape == REGION_POINT) then

    if (parmList(PARM_LIST_MODE) == FEEDBACK_MODE_QUASAR) then

      partBlk = UNKNOWN
      blkList(1) = lb
      call gr_findBlock(blkList, 1, pos, partBlk)
  
      if (partBlk == lb) then
    
        call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)
        i = blkLimits(LOW, IAXIS) + int((pos(1)-coord(1))/dx(1))
        j = blkLimits(LOW, JAXIS) + int((pos(2)-coord(2))/dx(2))
        k = blkLimits(LOW, KAXIS) + int((pos(3)-coord(3))/dx(3))
   
        call quasarFeedback(mgdot, dt, (/ 0.0, 0.0, 0.0 /), 0.5, Eej)

        call Grid_getBlkPtr(lb, solnData)
         
        solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) * & 
                                   solnData(DENS_VAR,i,j,k) + Eej/dV
        solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) / & 
                                   solnData(DENS_VAR,i,j,k)
        solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) * & 
                                   solnData(DENS_VAR,i,j,k) + Eej/dV
        solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) / & 
                                   solnData(DENS_VAR,i,j,k)
        call Grid_releaseBlkPtr(lb, solnData) 

        doneApplying = .true.
      end if ! point is actually on block
    end if ! quasars

  ! apply extended sources here
  else
    
    ! check if injection region intersects current block
    regionRadius = getRadius(regionShape, hej, rej, phi, theta)
    regionVolume = getVolume(regionShape, hej, rej, phi, theta)
    volIntersect =getVolIntersect(pos,regionRadius,blkCenter,blkRadius,agn_Lbox)
    volRemaining = parmList(PARM_LIST_VOL) - volIntersect
    volFracRemaining = volRemaining/regionVolume

    doneApplying = (parmList(PARM_LIST_VOL)/regionVolume <= VOL_TOL)

    if (.not. doneApplying .and. volIntersect > 0.) then

      call Grid_getBlkPtr(lb, solnData)
    
      call getMagFieldStrength(dr_simTime, dt, mgdot, mass, oldMass, &
                               mode, magRegionRadius, magStrength)

      ! avoid magnetic monopoles - only apply magnetic field if we 
      !   are resolved enough (only warn once per timestep)
      magResolved = (magRegionRadius > dx(1))
      call Driver_getNStep(nstep)
      if (dr_globalMe == MASTER_PE .and. nstep /= oldNstep .and. & 
          .not. magResolved) then
        call Logfile_stamp(dr_globalMe, "WARNING: Not enough resolution to apply magnetic field injection!")
        oldNstep = nstep
      end if

      ! unpack zones
      call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)

      call Grid_getCellCoords(KAXIS, lb, CENTER, .true., zoneZ, GRID_KHI_GC)
      call Grid_getCellCoords(JAXIS, lb, CENTER, .true., zoneY, GRID_JHI_GC)
      call Grid_getCellCoords(IAXIS, lb, CENTER, .true., zoneX, GRID_IHI_GC)

      do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
      do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        zoneCoord(1) = zoneX(i)
        zoneCoord(2) = zoneY(j)
        zoneCoord(3) = zoneZ(k)
        displ(:) = getDistance(pos, zoneCoord, phi, theta, agn_Lbox)

        ! apply magnetic field injection
#ifdef MAGP_VAR
        if (magResolved .and. mode /= FEEDBACK_MODE_QUASAR) then
          Bx = solnData(MAGX_VAR,i,j,k)
          By = solnData(MAGY_VAR,i,j,k)
          Bz = solnData(MAGZ_VAR,i,j,k)

          call magFeedback(dt, magRegionScale, magStrength, & 
                           displ, Bxej, Byej, Bzej, magEnergy)

          call rotateBack(phi, theta, Bxej, Byej, Bzej, Bx, By, Bz) 

          solnData(MAGX_VAR,i,j,k) = solnData(MAGX_VAR,i,j,k) + Bx
          solnData(MAGY_VAR,i,j,k) = solnData(MAGY_VAR,i,j,k) + By
          solnData(MAGZ_VAR,i,j,k) = solnData(MAGZ_VAR,i,j,k) + Bz

#ifdef MAG_FACE_VAR
#endif
          totalB = totalB + magEnergy*dV
        end if ! apply magnetic field
#endif

        ! apply thermal and kinetic injection
        Eej = 0.
        Mej = 0.
        Pej = 0.
        Dej = 0.
        vzej = 0.
        Ecrej = 0.

        if (mode == FEEDBACK_MODE_COSMICRAY) then
          isInside = isInRegion(regionShape, pos, hej, rej, &
                                phi, theta, zoneCoord, agn_Lbox)

#ifdef CRAY_MSCALAR 

!! If inject CR jets using inflow BC from hy_uhd_unsplit, mute here
!!          isInside = .false.

          if(isInside) then
            call cosmicrayFeedback(displ, rej, hej, dt, Mej, Pej, Eej, Ecrej)

            solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) + 0.
            solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) + 0.
            solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) + Pej
            solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)*solnData(DENS_VAR,i,j,k) + Eej + Ecrej
            solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k) = solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k)*solnData(DENS_VAR,i,j,k) + Ecrej

            solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k) + Mej

            ecr = sum(solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k))
            pcr = ecr*(eos_gamma_cr-1.)
            egas = max(solnData(EINT_VAR,i,j,k) - ecr, 0.)
            pgas = egas*(eos_gamma_gas-1.)

            solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) / &
                                       solnData(DENS_VAR,i,j,k)
            solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) / & 
                                       solnData(DENS_VAR,i,j,k)
            solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) / &
                                       solnData(DENS_VAR,i,j,k)
            solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) / &
                                       solnData(DENS_VAR,i,j,k)
            solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k) = solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k) / solnData(DENS_VAR,i,j,k)

            solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                                       0.5*(solnData(VELX_VAR,i,j,k)**2. + &
                                            solnData(VELY_VAR,i,j,k)**2. + &
                                            solnData(VELZ_VAR,i,j,k)**2.)

            solnData(GAMC_VAR,i,j,k) = (eos_gamma_gas*pgas + eos_gamma_cr*pcr) / &
                                       (pgas + pcr)
            solnData(GAME_VAR,i,j,k) = (eos_gamma_gas*(eos_gamma_cr-1.)*pgas + &
                                        eos_gamma_cr*(eos_gamma_gas-1.)*pcr) / &
                                      ((eos_gamma_cr-1.)*pgas + (eos_gamma_gas-1.)*pcr)

          endif
#endif

        else

         if (mode == FEEDBACK_MODE_QUASAR) then
          totalEej = 0
          do ii = 1, nsubzones
          do jj = 1, nsubzones
          do kk = 1, nsubzones
            cellLeftEdge(:) = zoneCoord(:) - dx(:)/2.
            subCoord(1) = cellLeftEdge(1) + (ii-0.5) * subdx(1)
            subCoord(2) = cellLeftEdge(2) + (jj-0.5) * subdx(2)
            subCoord(3) = cellLeftEdge(3) + (kk-0.5) * subdx(3)
            subDispl(:) = getDistance(pos, subCoord, phi, theta, agn_Lbox)
            call quasarFeedback(mgdot, dt, subDispl, rej, Eej)

            totalEej = totalEej + Eej * subdV / dV
          end do
          end do
          end do
            
          solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) * & 
                                     solnData(DENS_VAR,i,j,k) + totalEej
         else if (mode == FEEDBACK_MODE_JET) then
          call jetFeedback(mgdot, dt, displ, rej, hej, Mej, Pej, Eej, Ecrej)
          solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) * &
                                     solnData(DENS_VAR,i,j,k) + Eej + Ecrej
#ifdef CRAY_MSCALAR
          if(NMASS_SCALARS > 0. .and. Ecrej > 0.) then
            solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k) = &
              solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k) * &
              solnData(DENS_VAR,i,j,k) + Ecrej
          endif
#endif
          totalE = totalE + (Eej+Ecrej)*dV
         else if (mode == FEEDBACK_MODE_BUBBLE) then
          call bubbleFeedback(mgdot, mass, oldMass, dt, displ, &
                              totalMass, rej, Eej, Mej, Pej)
          solnData(EINT_VAR,i,j,k) = (solnData(EINT_VAR,i,j,k) + Eej) * &
                                      solnData(DENS_VAR,i,j,k)
          totalE = totalE + Eej*solnData(DENS_VAR,i,j,k)*dV
         end if

          vxej = Pej * sin(theta) * cos(phi) 
          vyej = Pej * sin(theta) * sin(phi)
          vzej = Pej * cos(theta)

          solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) * & 
                                     solnData(DENS_VAR,i,j,k) + vxej
          solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) * & 
                                     solnData(DENS_VAR,i,j,k) + vyej
          solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) * & 
                                     solnData(DENS_VAR,i,j,k) + vzej

          Mej = Mej / solnData(DENS_VAR,i,j,k)
          solnData(DENS_VAR,i,j,k) = solnData(DENS_VAR,i,j,k) * (1. + Mej)

#ifdef CRAY_MSCALAR
          if(NMASS_SCALARS > 0. .and. Ecrej > 0.) then
            ecr = sum(solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k))
            pcr = ecr * (eos_gamma_cr-1.)
            egas = max(solnData(EINT_VAR,i,j,k) - ecr, 0.)
            pgas = egas * (eos_gamma_gas-1.)
            solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k) = &
              solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k) / &
              solnData(DENS_VAR,i,j,k)
            solnData(GAMC_VAR,i,j,k) = (eos_gamma_gas*pgas + eos_gamma_cr*pcr)/&
                                       (pgas + pcr)
            solnData(GAME_VAR,i,j,k) = (eos_gamma_gas*(eos_gamma_cr-1.)*pgas + &
                                        eos_gamma_cr*(eos_gamma_gas-1.)*pcr) / &
                                      ((eos_gamma_cr-1.)*pgas + (eos_gamma_gas-1.)*pcr)
          endif
#endif

          solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) / & 
                                     solnData(DENS_VAR,i,j,k)
          solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) / & 
                                     solnData(DENS_VAR,i,j,k)
          solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) / & 
                                     solnData(DENS_VAR,i,j,k)
          solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) / & 
                                     solnData(DENS_VAR,i,j,k)

          solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + & 
                                     0.5*(solnData(VELX_VAR,i,j,k)**2. + &
                                          solnData(VELY_VAR,i,j,k)**2. + &
                                          solnData(VELZ_VAR,i,j,k)**2.)

! inject tracer fluid for jet material
#ifdef FJET_MSCALAR
          if (agn_useTracers) then
            isInside = isInRegion(regionShape, pos, hej, rej, &
                                  phi, theta, zoneCoord, agn_Lbox)
            if(isInside) &   
              solnData(FJET_MSCALAR,i,j,k) = 1.d0
          endif
#endif

        endif ! mode=cosmicray
      end do
      end do
      end do ! loop over zones

      parmList(PARM_LIST_VOL) = volRemaining
      doneApplying = (volFracRemaining <= VOL_TOL)
      
      call Grid_releaseBlkPtr(lb, solnData)
      
      ! inject particle tracers into this block  
      call Driver_getSimTime(time)
      if (agn_useTracers) then
        if (mode == FEEDBACK_MODE_JET .and. time > agn_nextInjectionTime) then
          call InjectTracers(lb, parmList, volIntersect/regionVolume)
        else if (mode == FEEDBACK_MODE_BUBBLE) then
          call InjectTracers(lb, parmList, volIntersect/regionVolume)
        end if
      end if

    end if ! inject in this block

  end if ! point versus extended sources
 
  return
end subroutine applyFeedback

end module
