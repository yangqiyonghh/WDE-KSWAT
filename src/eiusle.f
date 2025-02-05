      subroutine eiusle
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine computes the USLE erosion index (EI)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    amp_r(:,:)  |none        |alpha factor for rain (mo max 0.5h rain)
!!    idg(:)      |none        |array location of random number seed
!!                             |used for a given process
!!    ihru        |none        |HRU number
!!    i_mo        |none        |month being simulated
!!    ovrlnd(:)   |mm H2O      |overland flow onto HRU from upstream
!!                             |routing unit
!!    precipday   |mm H2O      |amount of water reaching soil surface
!!    rndseed(:,:)|none        |random number generator seed
!!    snomlt      |mm H2O      |amount of snow melt in HRU on current day
!!    tconc(:)    |hr          |time of concentration
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units                  |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    usle_ei     |100(ft-tn in)/(acre-hr)|USLE rainfall erosion index
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units       |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ab          |mm H2O      |lowest value xa can have
!!    ajp         |mm H2O      |highest value xa can have
!!    j           |none        |HRU number
!!    pkrf        |none        |intermediate calculation
!!    pkrf30      |mm/hr       |maximum 30-min. storm intensity
!!    preceff     |mm H2O      |amount of rainfall impacting ground surface
!!                             |on day in HRUHRU当天影响地面的降雨量
!!    xa          |none        |fraction of daily rainfall occuring during
!!                             |half-hour of max intensity rain最大强度降雨半小时内的日降雨量
!!    xb          |none        |intermediate calculation中间计算
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Log, Log10
!!    SWAT: Expo, Atri

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: j
      real :: ab, xa, preceff, ajp, xb, pkrf, pkrf30

      j = 0
      j = ihru

      ab = 0.02083
      pkrf = 0.
      pkrf30 = 0.
      preceff = 0.
      xa = 0.
      xb = 0.

      preceff = precipday - snomlt - ovrlnd(j)
      if (preceff > 1.e-4) then
        ajp = 1. - Expo(-125. / preceff)
        xa = Atri(ab, amp_r(i_mo,hru_sub(j)), ajp, rndseed(idg(4),j))
        xb = -2. * Log(1. - xa)!式子4-25
        pkrf30 = 2. * preceff * xa!式子4-26
        pkrf = xb * preceff!式子4-25
        usle_ei = preceff * (12.1 + 8.9 * (Log10(pkrf) - .4343)) *     !式子4-21 
     &    pkrf30 / 1000.
!    *    peakr / 10.
        if (usle_ei < 1.e-4) usle_ei = 0.
        usle_eifac(j) = usle_ei
      endif

      return
      end