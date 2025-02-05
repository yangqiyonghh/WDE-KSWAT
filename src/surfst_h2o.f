      subroutine surfst_h2o

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine determines the net surface runoff reaching the 
!!    main channel on a given day. The net amount of water reaching
!!    the main channel can include water in surface runoff from the 
!!    previous day and will exclude surface runoff generated on the
!!    current day which takes longer than one day to reach the main
!!    channel

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    brt(:)      |none          |fraction of surface runoff that takes
!!                               |one day or less to reach the subbasin
!!                               |outlet需要一天或更短时间才能到达子流域出口的地表径流部分
!!    ihru        |none          |HRU number
!!    surf_bs(1,:)|mm H2O        |amount of surface runoff lagged over one
!!                               |day
!!    surfq(:)    |mm H2O        |surface runoff generated in HRU on the
!!                               |current day
!!    hhqday(:)   |mm H2O        |surface runoff generated in HRU on the
!!                               |current hour at current day  
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bsprev      |mm H2O        |surface runoff lagged from prior day
!!    qday        |mm H2O        |surface runoff loading to main channel
!!                               |from HRU on current day
!!    surf_bs(1,:)|mm H2O        |amount of surface runoff lagged over one
!!                               |day
!!    hhqday(:)   |mm H2O        |surface runoff generated in HRU on the
!!                               |current hour at current day  当前小时当天HRU中产生的地表径流
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


      use parm

      integer :: j
      
      j = 0
      j = ihru

      if (ievent == 0) then
          
         !
         !inflpcp= max((1-sum(sol_st(:,j))/sum(sol_fc(:,j))),0.0)*surf_bs(1,j) 
         ! surf_bs(1,j) =surf_bs(1,j)-inflpcp

         bsprev = surf_bs(1,j)!滞后一天的地表径流量
	   surf_bs(1,j) = Max(1.e-6, surf_bs(1,j) + surfq(j))
         qday = surf_bs(1,j) * brt(j)
         surf_bs(1,j) = surf_bs(1,j) - qday
	    !if (iyr==2008 .and.(i==164 .or. i==165 .or. i==163).and. (J==1))print*, qday ,surf_bs(1,j) , brt(j)
      else
		!亚日径流滞后（仅适用于透水表面的径流）
         bsprev = hhsurf_bs(1,j,nstep)		! 与前一天相比的滞后 J.Jeong 4/06/2009

	   do k=1,nstep

	   !! Left-over (previous timestep) + inflow (current  timestep)
           hhsurf_bs(1,j,k) = Max(0., bsprev + hhqday(k))
   	
	   !! 到达主河道的径流和泥沙的新估计
	     hhqday(k) = hhsurf_bs(1,j,k) * brt(j)
	     hhsurf_bs(1,j,k) = hhsurf_bs(1,j,k) - hhqday(k)
   	  
	   !! lagged at the end of time step  
	     bsprev = hhsurf_bs(1,j,k) 
	   end do

	   !! daily total yield from the HRU
	   qday = sum(hhqday) + sum(ubnrunoff)  !!!!这里把一天的进行相加了，小时尺度开发
	end if
     
      return
      end