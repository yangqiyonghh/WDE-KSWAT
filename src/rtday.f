      subroutine rtday

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine routes the daily flow through the reach using a 
!!    variable storage coefficient

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ch_d(:)     |m             |average depth of main channel
!!    ch_k(2,:)   |mm/hr         |effective hydraulic conductivity of
!!                               |main channel alluvium
!!    ch_l2(:)    |km            |length of main channel
!!    ch_n(2,:)   |none          |Manning's "n" value for the main channel
!!    ch_s(2,:)   |m/m           |average slope of main channel
!!    ch_w(2,:)   |m             |average width of main channel
!!    chside(:)   |none          |change in horizontal distance per unit河道边坡上每单位垂直距离变化的水平距离变化；始终设置为2（斜率=1/2）
!!                               |change in vertical distance on channel side
!!                               |slopes; always set to 2 (slope=1/2)河道边坡上每单位垂直距离变化的水平距离变化；始终设置为2（斜率=1/2）
!!    evrch       |none          |Reach evaporation adjustment factor.
!!                               |Evaporation from the reach is multiplied by
!!                               |EVRCH. This variable was created to limit the
!!                               |evaporation predicted in arid regions.
!!    inum1       |none          |reach number
!!    inum2       |none          |inflow hydrograph storage location number流入过程线存储位置编号
!!    pet_day     |mm H2O        |potential evapotranspiration
!!    phi(1,:)    |m^2           |cross-sectional area of flow in channel at河道满岸深度处的水流截面积
!!                               |bankfull depth
!!    phi(6,:)    |m             |bottom width of main channel主渠道底部宽度
!!    rnum1       |none          |fraction of overland flow
!!    rchstor(:)   |m^3 H2O       |water stored in reach
!!    varoute(2,:)|m^3 H2O       |water flowing into reach on day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    rcharea     |m^2           |cross-sectional area of flow
!!    rchdep      |m             |depth of flow on day
!!    rtevp       |m^3 H2O       |evaporation from reach on day
!!    rttime      |hr            |reach travel time
!!    rttlc       |m^3 H2O       |transmission losses from reach on day
!!    rtwtr       |m^3 H2O       |water leaving reach on day
!!    sdti        |m^3/s         |average flow on day in reach
!!    rchstor(:)   |m^3 H2O       |water stored in reach
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    det         |hr            |time step (24 hours)
!!    c           |none          |inverse of channel side slope
!!    jrch        |none          |reach number
!!    p           |m             |wetted perimeter润湿周边
!!    rh          |m             |hydraulic radius
!!    scoef       |none          |Storage coefficient (fraction of water in 
!!                               |reach flowing out on day)
!!    topw        |m             |top width of main channel
!!    vol         |m^3 H2O       |volume of water in reach at beginning of
!!                               |day一天开始时可到达的水量
!!    wtrin       |m^3 H2O       |amount of water flowing into reach on day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sqrt, Min
!!    SWAT: Qman

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
!!    Modified by Balaji Narasimhan
!!    Spatial Sciences Laboratory, Texas A&M University
      use parm

      integer :: jrch
      real :: wtrin, scoef, p, tbase, topw, vol, c, rh
	  real :: volrt, maxrt, adddep, addp, addarea, vc, aaa
	  real :: rttlc1, rttlc2, rtevp1, rtevp2, det

      jrch = 0
      jrch = inum1

      wtrin = 0.
      wtrin = varoute(2,inum2) * (1. - rnum1)

!! calculate volume of water in reach计算河段的水量
      vol = 0.
      vol = wtrin + rchstor(jrch)
!! Find average flowrate in a day求出一天内的平均流量
      volrt = vol / 86400.

!! Find maximum flow capacity of the channel at bank full查找满岸时河道的最大流量
      c = 0.
      c = chside(jrch)
	p = phi(6,jrch) + 2. * ch_d(jrch) * Sqrt(1. + c * c)!式子7-5
	rh = phi(1,jrch) / p!式子7-6水力半径
	maxrt = Qman(phi(1,jrch), rh, ch_n(2,jrch), ch_s(2,jrch))

      sdti = 0.
	rchdep = 0.
	p = 0.
	rh = 0.
	vc = 0.

!! If average flowrate is greater than than the channel capacity at bank full
!! then simulate flood plain flow else simulate the regular channel flow
!! 如果平均流量大于满岸时的河道容量，则模拟洪泛平原流量，否则模拟常规河道流量
      
      if (volrt > maxrt) then
	  rcharea = phi(1,jrch)
	  rchdep = ch_d(jrch)
	  p = phi(6,jrch) + 2. * ch_d(jrch) * Sqrt(1. + c * c)!式子7-5
	  rh = phi(1,jrch) / p!式子7-6水力半径
	  sdti = maxrt
	  adddep = 0
        !PRINT*,volrt,jrch,I
!! 用迭代法在1cm的间隔深度找到volrt的横截面积和深度找到深度，直到放电速率等于volrt
        
	  Do While (sdti < volrt)
          adddep = adddep + 0.01
          addarea = rcharea + ((ch_w(2,jrch) * 5) + 4 * adddep) * adddep!式子7-9
          addp = p + (ch_w(2,jrch) * 4) + 2. * adddep * Sqrt(1. + 4 * 4)!式子7-10
	    rh = addarea / addp
          sdti = Qman(addarea, rh, ch_n(2,jrch), ch_s(2,jrch))
          
        end do

	  rcharea = addarea
	  rchdep = ch_d(jrch) + adddep
	  p = addp
	  sdti = volrt
	else
!! 用迭代法在1cm的间隔深度处找到volrt的横截面积和深度找到深度，直到放电速率等于volrt
	  Do While (sdti < volrt)
	    rchdep = rchdep + 0.01
	    rcharea = (phi(6,jrch) + c * rchdep) * rchdep
	    p = phi(6,jrch) + 2. * rchdep * Sqrt(1. + c * c)
	    rh = rcharea / p
          sdti = Qman(rcharea, rh, ch_n(2,jrch), ch_s(2,jrch))
	  end do
	  sdti = volrt
	end if

!!计算水位处渠道的顶部宽度
      topw = 0.
      if (rchdep <= ch_d(jrch)) then
        topw = phi(6,jrch) + 2. * rchdep * c
      else
        topw = 5 * ch_w(2,jrch) + 2. * (rchdep - ch_d(jrch)) * 4.
      end if

!!	模拟的时间步长（小时）
        det = 24.

      if (sdti > 0.) then
        !! 计算速度和行程时间
  	    vc = sdti / rcharea  
        vel_chan(jrch) = vc
	  rttime = ch_l2(jrch) * 1000. / (3600. * vc)


        !! 计算当天离开河段的水量
        scoef = 0.
 	  rtwtr = 0.
        !scoef = 2. * det / (2. * rttime + det)
        scoef =  det / (rttime + det)
        if (scoef > 1.) scoef = 1.
        rtwtr = scoef * (wtrin + rchstor(jrch))
        !新存储系数替换
        rtwtr = vc * rcharea * 86400.
        rtwtr = amin1 (rtwtr, wtrin)

!! 计算一天结束时河道中的水量
      rchstor(jrch) = rchstor(jrch) + wtrin - rtwtr

!! 添加if语句以防止rchstor变为负数
      if (rchstor(jrch) < 0.0) rchstor(jrch) = 0.0

!! 传输和蒸发损失按比例从通道储存和流出的体积中提取

       !! 计算传输损耗
	  rttlc = 0.

	  if (rtwtr > 0.) then

	!!  Total time in hours to clear the water净水总时间（小时）

          rttlc = det * ch_k(2,jrch) * ch_l2(jrch) * p
          rttlc2 = rttlc * rchstor(jrch) / (rtwtr + rchstor(jrch))

	    if (rchstor(jrch) <= rttlc2) then
	      rttlc2 = min(rttlc2, rchstor(jrch))
	      rchstor(jrch) = rchstor(jrch) - rttlc2
	      rttlc1 = rttlc - rttlc2
	      if (rtwtr <= rttlc1) then
	        rttlc1 = min(rttlc1, rtwtr)
	        rtwtr = rtwtr - rttlc1
	      else
	        rtwtr = rtwtr - rttlc1
	      end if
	    else
	      rchstor(jrch) = rchstor(jrch) - rttlc2
	      rttlc1 = rttlc - rttlc2
	      if (rtwtr <= rttlc1) then
	        rttlc1 = min(rttlc1, rtwtr)
	        rtwtr = rtwtr - rttlc1
	      else
	        rtwtr = rtwtr - rttlc1
	      end if
	    end if
	  rttlc = rttlc1 + rttlc2
        end if


        !! calculate evaporation计算蒸发量
	  rtevp = 0.
       if (rtwtr > 0.) then

          aaa = evrch * pet_day / 1000.

	    if (rchdep <= ch_d(jrch)) then
            rtevp = aaa * ch_l2(jrch) * 1000. * topw
	    else
		  if (aaa <=  (rchdep - ch_d(jrch))) then
              rtevp = aaa * ch_l2(jrch) * 1000. * topw
	      else
	        rtevp = (rchdep - ch_d(jrch)) 
	        rtevp = rtevp + (aaa - (rchdep - ch_d(jrch))) 
              topw = phi(6,jrch) + 2. * ch_d(jrch) * c           
	        rtevp = rtevp * ch_l2(jrch) * 1000. * topw
	      end if
	    end if

	    rtevp2 = rtevp * rchstor(jrch) / (rtwtr + rchstor(jrch))

	    if (rchstor(jrch) <= rtevp2) then
	      rtevp2 = min(rtevp2, rchstor(jrch))
	      rchstor(jrch) = rchstor(jrch) - rtevp2
	      rtevp1 = rtevp - rtevp2
	      if (rtwtr <= rtevp1) then
	        rtevp1 = min(rtevp1, rtwtr)
	        rtwtr = rtwtr - rtevp1
	      else
	        rtwtr = rtwtr - rtevp1
	      end if
	    else
	      rchstor(jrch) = rchstor(jrch) - rtevp2
	      rtevp1 = rtevp - rtevp2
	      if (rtwtr <= rtevp1) then
	        rtevp1 = min(rtevp1, rtwtr)
	        rtwtr = rtwtr - rtevp1
	      else
	        rtwtr = rtwtr - rtevp1
	      end if
	    end if
	  rtevp = rtevp1 + rtevp2
        end if

      else
        rtwtr = 0.
        sdti = 0.
	  rchstor(jrch) = 0.
	  vel_chan(jrch) = 0.
        flwin(jrch) = 0.
        flwout(jrch) = 0.
      end if

!! 未计算河段的降水量，因为子流域中的HRU面积加起来等于整个子流域面积（包括河道面积），因此降水量计入子流域环流

!!      volinprev(jrch) = wtrin
!!	qoutprev(jrch) = rtwtr

      if (rtwtr < 0.) rtwtr = 0.
      if (rchstor(jrch) < 0.) rchstor(jrch) = 0.

      if (rchstor(jrch) < 10.) then
        rtwtr = rtwtr + rchstor(jrch)
        rchstor(jrch) = 0.
      end if

      return
      end