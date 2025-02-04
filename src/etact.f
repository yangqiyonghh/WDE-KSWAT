      subroutine etact
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates potential plant transpiration for Priestley-
!!    Taylor and Hargreaves ET methods, and potential and actual soil
!!    evaporation. NO3 movement into surface soil layer due to evaporation
!!    is also calculated.


!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    canstor(:)   |mm H2O        |amount of water held in canopy storage
!!    elevb(:,:)   |m             |elevation at center of band in subbasin
!!    elevb_fr(:,:)|none          |fraction of subbasin area within elevation 
!!                                |band
!!    ep_max       |mm H2O        |maximum amount of transpiration (plant et)
!!                                |that can occur on current day in HRU 
!!    esco(:)      |none          |soil evaporation compensation factor
!!    ihru         |none          |HRU number
!!    ipet         |none          |code for potential ET method
!!                                |0 Priestley-Taylor method
!!                                |1 Penman/Monteith method
!!                                |2 Hargreaves method
!!    laiday(:)    |m**2/m**2     |leaf area index
!!    pet_day      |mm H2O        |potential evapotranspiration on current day
!!                                |in HRU
!!    pot_vol(:)   |m**3 H2O      |current volume of water stored in the
!!                                |depression/impounded area
!!    sno_hru(:)   |mm H2O        |amount of water in snow in HRU on current day当天HRU的雪中水量
!!    snoeb(:,:)   |mm H2O        |snow water content in elevation band on 
!!                                |current day
!!    sol_cov(:)   |kg/ha         |amount of residue on soil surface
!!    sol_fc(:,:)  |mm H2O        |amount of water available to plants in soil
!!                                |layer at field capacity (fc - wp water)
!!    sol_nly(:)   |none          |number of soil layers in profile
!!    sol_no3(:,:) |kg N/ha       |amount of nitrogen stored in the nitrate
!!                                |pool 
!!    sol_st(:,:)  |mm H2O        |amount of water stored in the soil layer on
!!                                |current day
!!    sol_z(:,:)   |mm            |depth to bottom of soil layer
!!    tavband(:,:) |deg C         |average temperature for the day in band in HRUHRU频带内当天的平均温度
!!    tmpav(:)     |deg C         |average air temperature on current day for
!!                                |HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    canev        |mm H2O        |amount of water evaporated from canopy
!!                                |storage
!!    ep_max       |mm H2O        |maximum amount of transpiration (plant et)
!!                                |that can occur on current day in HRU
!!    es_day       |mm H2O        |actual amount of evaporation (soil et) that
!!                                |occurs on day in HRU
!!    sno_hru(:)   |mm H2O        |amount of water in snow in HRU on current day
!!    sno3up       |kg N/ha       |amount of nitrate moving upward in the soil
!!                                |profile in watershed
!!    snoeb(:,:)   |mm H2O        |snow water content in elevation band on 
!!                                |current day
!!    snoev        |mm H2O        |amount of water in snow lost through
!!                                |sublimation on current day当天通过升华损失的雪中水量
!!    sol_st(:,:)  |mm H2O        |amount of water stored in the soil layer on
!!                                |current day
!!    sol_sw(:)    |mm H2O        |amount of water stored in the soil profile
!!                                |on current day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cej          |
!!    dep          |mm            |soil depth from which evaporation will occur
!!                                |in current soil layer当前土层将发生蒸发的土壤深度
!!    eaj          |none          |weighting factor to adjust PET for impact of
!!                                |plant cover
!!    effnup       |
!!    eos1         |none          |variable to hold intermediate calculation
!!                                |result
!!    eosl         |mm H2O        |maximum amount of evaporation that can occur
!!                                |from soil profile
!!    es_max       |mm H2O        |maximum amount of evaporation (soil et)
!!                                |that can occur on current day in HRU
!!    esd          |mm            |maximum soil depth from which evaporation允许发生蒸发的最大土壤深度
!!                                |is allowed to occur
!!    esleft       |mm H2O        |potenial soil evap that is still available仍然可用的潜在土壤蒸发量
!!    etco         |
!!    evz          |
!!    evzp         |
!!    ib           |none          |counter
!!    j            |none          |HRU number
!!    ly           |none          |counter
!!    no3up        |kg N/ha       |amount of nitrate moving upward in profile
!!    pet          |mm H2O        |amount of PET remaining after water stored
!!                                |in canopy is evaporated
!!    sev          |mm H2O        |amount of evaporation from soil layer
!!    sumsnoeb     |mm H2O        |amount of snow in elevation bands whose air
!!                                |temperature is greater than 0 degrees C
!!    xx           |none          |variable to hold intermediate calculation 
!!                                |result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Min, Max
!!    SWAT: Expo

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: j, ib, ly
!!    real, parameter :: esd = 500., etco = 0.80, effnup = 0.1
      real :: esd, etco, effnup
      real :: no3up, es_max, eos1, xx, cej, eaj, pet, esleft
      real :: sumsnoeb, evzp, eosl, dep, evz, sev

      j = 0
      j = ihru
      pet = 0.
      pet = pet_day
!!   添加语句以测试上面的真实语句
	esd = 500.
	etco = 0.80
	effnup = 0.1

!! evaporate canopy storage first
!! canopy storage is calculated by the model only if the Green & Ampt
!! method is used to calculate surface runoff. The curve number methods
!! take canopy effects into account in the equations. For either of the
!! CN methods, canstor will always equal zero.
!! 只有在使用Green&Ampt方法计算地表径流的情况下，该模型才能计算蒸发冠层蓄水量第一冠层蓄水量。
!曲线数方法在方程中考虑了冠层效应。对于任意一种CN方法，canstor将始终等于零。
      pet = pet - canstor(j)
      if (pet < 0.) then
        canstor(j) = -pet
        canev = pet_day
        pet = 0.
        ep_max = 0.
        es_max = 0.
      else
        canev = canstor(j)
        canstor(j) = 0.
      endif

      if (pet > 1.0e-6) then

        !! compute potential plant evap for methods other that Penman-Monteith用Penman-Monteith以外的方法计算潜在的电厂蒸发量
        if (ipet /= 1) then
          if (laiday(j) <= 3.0) then
            ep_max = laiday(j) * pet / 3.!式子2-81
          else
            ep_max = pet!式子2-82
          end if
          if (ep_max < 0.) ep_max = 0.
        end if

        !! compute potential soil evaporation计算潜在土壤蒸发量
        cej = -5.e-5
        eaj = 0.
        es_max = 0.
        eos1 = 0.
        if (sno_hru(j) >= 0.5) then
          eaj = 0.5
        else
          eaj = Exp(cej * (sol_cov(j)+ 0.1))!式子2-84
        end if
        es_max = pet * eaj!式子2-83
        eos1 = pet / (es_max + ep_max + 1.e-10)!式子2-85
        eos1 = es_max * eos1!式子2-85
        es_max = Min(es_max, eos1)
        es_max = Max(es_max, 0.)
!        if (pot_vol(j) > 1.e-4) es_max = 0.

        !! make sure maximum plant and soil ET doesn't exceed potential ET确保最大植物和土壤ET不超过潜在ET
        !!if (pet_day < es_max + ep_max) then
          !!es_max = pet_day - ep_max
          if (pet < es_max + ep_max) then
            es_max = pet * es_max / (es_max + ep_max)
            ep_max = pet * ep_max / (es_max + ep_max)
          end if
          if (pet < es_max + ep_max) then
            es_max = pet - ep_max - 1.0e-6
          end if
        !!end if

        !! initialize soil evaporation variables初始化土壤蒸发变量
        esleft = 0.
        esleft = es_max

        !! compute sublimation计算升华
        if (elevb_fr(1,hru_sub(j)) <= 0.) then
          !! compute sublimation without elevation bands计算无高程带的升华
          if (tmpav(j) > 0.) then
            if (sno_hru(j) >= esleft) then
              !! take all soil evap from snow cover清除积雪中的所有土壤
              sno_hru(j) = sno_hru(j) - esleft
              snoev = snoev + esleft
              esleft = 0.
            else
              !! take all soil evap from snow cover then start taking from soil从积雪中提取所有土壤蒸发物，然后开始从土壤中提取
              esleft = esleft - sno_hru(j)
              snoev = snoev + sno_hru(j)
              sno_hru(j) = 0.
            endif
          endif
        else
          !! elevation bands
          sumsnoeb = 0.
          !! calculate air temp in elevation bands and sum snow
          !! for elevation bands with temp > 0 deg C计算海拔带中的空气温度，并对温度>0摄氏度的海拔带求和降雪量
          do ib = 1, 10
            if (elevb_fr(ib,hru_sub(j)) <= 0.) exit
            if (tavband(ib,j) > 0.) then 
              sumsnoeb = sumsnoeb +                                     
     &                             snoeb(ib,j) * elevb_fr(ib,hru_sub(j))
            end if
          end do
          
          !! compute sublimation from elevation bands从高程带计算升华
          if (sumsnoeb >= esleft .and. sumsnoeb > 0.01) then
            do ib = 1, 10
              if (elevb_fr(ib,hru_sub(j)) <= 0.) exit
              if (tavband(ib,j) > 0.) then
                snoev = snoev + snoeb(ib,j) * (esleft / sumsnoeb) *     
     &                                           elevb_fr(ib,hru_sub(j))
                snoeb(ib,j) = snoeb(ib,j) - snoeb(ib,j) * (esleft /     
     &                                                         sumsnoeb)
              end if
            end do
          else
            do ib = 1, 10
              if (elevb_fr(ib,hru_sub(j)) <= 0.) exit
              if (tavband(ib,j) > 0.) then
                snoev = snoev + snoeb(ib,j) * elevb_fr(ib,hru_sub(j))
                snoeb(ib,j) = 0.
              end if
            end do
          end if
          esleft = esleft - snoev
          sno_hru(j) = sno_hru(j) - snoev
        endif

!! take soil evap from each soil layer从每个土层中提取土壤蒸发量
      evzp = 0.
      eosl = 0.
      eosl = esleft
      do ly = 1, sol_nly(j)

        !! depth exceeds max depth for soil evap (esd)深度超过土壤蒸发的最大深度
        dep = 0.
        if (ly == 1) then
          dep = sol_z(1,j)
        else
          dep = sol_z(ly-1,j)
        endif
        
        if (dep < esd) then
          !! calculate evaporation from soil layer计算土层蒸发量
          evz = 0.
          sev = 0. 
          xx = 0.
          evz = eosl * sol_z(ly,j) / (sol_z(ly,j) + Exp(2.374 -   !式子2-92
     &       .00713 * sol_z(ly,j)))
          sev = evz - evzp * esco(j)!式子2-94
          evzp = evz
          if (sol_st(ly,j) < sol_fc(ly,j)) then
            xx =  2.5 * (sol_st(ly,j) - sol_fc(ly,j)) / sol_fc(ly,j)!式子2-95
            sev = sev * Expo(xx)!式子2-95
          end if
          sev = Min(sev, sol_st(ly,j) * etco)!式子2-97

          if (sev < 0.) sev = 0.
          if (sev > esleft) sev = esleft
			
          !! adjust soil storage, potential evap调整土壤储量、潜在蒸发量
          if (surf_bs(1,j)>sev)then
              surf_bs(1,j)=surf_bs(1,j)-sev
              sev=0.

          else
              sev=sev-surf_bs(1,j)
              surf_bs(1,j)=0.
           end if
          if (sol_st(ly,j) > sev) then
            esleft = esleft - sev
            sol_st(ly,j) = Max(1.e-6, sol_st(ly,j) - sev)
          else
            esleft = esleft - sol_st(ly,j)
            sol_st(ly,j) = 0.
          endif
        endif

        !! compute no3 flux from layer 2 to 1 by soil evaporation通过土壤蒸发调节土壤储存，从第2层到第1层的潜在蒸发计算no3通量
        if (ly == 2) then
          no3up = 0.
          no3up = effnup * sev * sol_no3(2,j) / (sol_st(2,j) + 1.e-6)
          no3up = Min(no3up, sol_no3(2,j))
          sno3up = sno3up + no3up * hru_dafr(j)
          sol_no3(2,j) = sol_no3(2,j) - no3up
          sol_no3(1,j) = sol_no3(1,j) + no3up
        endif

      end do

      !! update total soil water content更新土壤总含水量
      sol_sw(j) = 0.
      do ly = 1, sol_nly(j)
        sol_sw(j) = sol_sw(j) + sol_st(ly,j)
        
      end do

      !! calculate actual amount of evaporation from soil计算土壤的实际蒸发量
      es_day = es_max - esleft
      if (es_day < 0.) es_day = 0.

      end if
  
      return
      end
