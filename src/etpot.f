      subroutine etpot
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates potential evapotranspiration using one
!!    of three methods. If Penman-Monteith is being used, potential plant
!!    transpiration is also calculated.

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name       |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    laiday(:)  |m**2/m**2      |leaf area index
!!    albday     |none           |albedo for the day in HRUHRU当天的反照率
!!    cht(:)     |m              |canopy height冠层高度
!!    co2(:)     |ppmv           |CO2 concentration
!!    gsi(:)     |m/s            |maximum stomatal conductance
!!    hru_ra(:)  |MJ/m^2         |solar radiation for the day in HRUHRU当天的太阳辐射
!!    hru_rmx(:) |MJ/m^2         |maximum possible radiation for the day in HRU
!!    hru_sub(:) |none           |subbasin in which HRU is located
!!    icr(:)     |none           |sequence number of crop grown within the
!!                               |current year
!!    idplt(:)   |none           |land cover code from crop.dat
!!    igro(:)    |none           |land cover status code
!!                               |0 no land cover currently growing
!!                               |1 land cover growing
!!    ihru       |none           |HRU number
!!    ipet       |none           |code for potential ET method
!!                               |0 Priestley-Taylor method
!!                               |1 Penman/Monteith method
!!                               |2 Hargreaves method
!!                               |3 read in PET values
!!    nro(:)     |none           |sequence number of year in rotation
!!    petmeas    |mm H2O         |potential ET value read in for day
!!    rhd(:)     |none           |relative humidity for the day in HRU
!!    sno_hru(:) |mm H2O         |amount of water in snow in HRU on current day当天HRU的雪中水量
!!    sub_elev(:)|m              |elevation of HRU
!!    tmn(:)     |deg C          |minimum air temperature on current day in HRU
!!    tmpav(:)   |deg C          |average air temperature on current day in HRU
!!    tmx(:)     |deg C          |maximum air temperature on current day for HRU
!!    u10(:)     |m/s            |wind speed (measured at 10 meters above 
!!                               |surface)
!!    vpd2(:)    |(m/s)*(1/kPa)  |rate of decline in stomatal conductance per
!!               |               |unit increase in vapor pressure deficit
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ep_max      |mm H2O        |maximum amount of transpiration (plant et) HRU中当天可能发生的最大蒸腾量（植物et）
!!                               |that can occur on current day in HRU
!!    pet_day     |mm H2O        |potential evapotranspiration on current day inHRU当天的潜在蒸散量
!!                               |HRU
!!    vpd         |kPa           |vapor pressure deficitpor压力不足
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    chz         |cm            |vegetation height
!!    d           |cm            |displacement height for plant type
!!    dlt         |kPa/deg C     |slope of the saturation vapor pressure-
!!                               |temperature curve
!!    ea          |kPa           |saturated vapor pressure
!!    ed          |kPa           |actual vapor pressure
!!    fvpd        |kPa           |amount of vapro pressure deficit over vapro压力不足量超过阈值
!!                               |threshold value
!!    gma         |kPa/deg C     |psychrometric constant
!!    j           |none          |HRU number
!!    pb          |kPa           |mean atmospheric pressure
!!    pet_alpha  |none           |alpha factor in Priestley-Taylor PET equation
!!    ralb        |MJ/m2         |net incoming radiation for PETPET净入射辐射
!!    ralb1       |MJ/m2         |net incoming radiation
!!    ramm        |MJ/m2         |extraterrestrial radiation
!!    rbo         |none          |net emissivity
!!    rc          |s/m           |canopy resistance
!!    rho         |MJ/(m3*kPa)   |K1*0.622*xl*rho/pb
!!    rn          |MJ/m2         |net radiation净辐射
!!    rn_pet      |MJ/m2         |net radiation for continuous crop cover连续作物覆盖的净辐射
!!    rout        |MJ/m2         |outgoing radiation
!!    rto         |none          |cloud cover factor
!!    rv          |s/m           |aerodynamic resistance to sensible heat and
!!                               |vapor transfer
!!    tk          |deg K         |average air temperature on current day for HRU
!!    uzz         |m/s           |wind speed at height zz
!!    xl          |MJ/kg         |latent heat of vaporization
!!    xx          |kPa           |difference between vpd and vpthreshold
!!    zom         |cm            |roughness length for momentum transfer
!!    zov         |cm            |roughness length for vapor transfer
!!    zz          |cm            |height at which wind speed is determined确定风速的高度
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Log, Sqrt, Max, Min
!!    SWAT: Ee

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: j
      real :: tk, pb, gma, xl, ea, ed, dlt, ramm, ralb1, ralb, xx
      real :: rbo, rto, rn, uzz, zz, zom, zov, rv, rn_pet, fvpd
      real :: rc, rho, rout, d, chz, gsi_adj, pet_alpha

      !! initialize local variables
      j = 0
      j = ihru

      tk = 0.
      tk = tmpav(j) + 273.15
	
      !! 计算平均大气压
      pb = 0.
      pb = 101.3 - sub_elev(hru_sub(j)) *                               
     &                       (0.01152 - 0.544e-6 * sub_elev(hru_sub(j)))!式子1-51

      !!计算蒸发潜热
      xl = 0.
      xl = 2.501 - 2.361e-3 * tmpav(j)!式子1-49

      !! 计算湿度常数
      gma = 0.
      gma = 1.013e-3 * pb / (0.622 * xl)!式子1-50

!计算饱和蒸汽压、实际蒸汽压和蒸汽压不足
      ea = 0.
      ed = 0.
      vpd = 0.
      ea = Ee(tmpav(j))
      ed = ea * rhd(j)!式子1-46
      vpd = ea - ed!式子1-48

      !!计算饱和蒸汽压力曲线的斜率
      dlt = 0.
      dlt = 4098. * ea / (tmpav(j) + 237.3)**2!式子1-47
	


!! 确定电位ET

      select case (ipet)

       case (0)   !! 普里斯特利-泰勒潜在蒸散发法
     
!! 净辐射计算PET的净短波辐射
          ralb = 0.
          if (sno_hru(j) <= .5) then
            ralb = hru_ra(j) * (1.0 - 0.23)
          else
            ralb = hru_ra(j) * (1.0 - 0.8)!式子1-19，19页
          end if

        !! 计算长波净辐射

          !! SWAT手册中的净发射率方程2.2.20
          rbo = 0.
          rbo = -(0.34 - 0.139 * Sqrt(ed))!式子1-28

          !! 云量因子方程2.2.19
          rto = 0.
            if (hru_rmx(j) < 1.e-4) then
		    rto = 0.
            else
              rto = 0.9 * (hru_ra(j) / hru_rmx(j)) + 0.1!式子1-28
            end if

          !! 长波净辐射方程2.2.21
          rout = 0.
          rout = rbo * rto * 4.9e-9 * (tk**4)!式子1-28

          !! 计算净辐射
          rn_pet = 0.
          rn_pet = ralb + rout!式子1-18
       !! 净辐射

          pet_alpha = 1.28
          pet_day = pet_alpha * (dlt / (dlt + gma)) * rn_pet / xl
          pet_day = Max(0., pet_day)


       case (1)   !! 彭曼蒙特潜在蒸散发法

       !! 净辐射
         !! 计算PET的净短波辐射
          ralb = 0.
          if (sno_hru(j) <= .5) then!式子1-19，19页
            ralb = hru_ra(j) * (1.0 - 0.23) !式子1-18
          else
            ralb = hru_ra(j) * (1.0 - 0.8) 
          end if
         !! 计算最大工厂ET的净短波辐射
          ralb1 = 0.
          ralb1 = hru_ra(j) * (1.0 - albday) !1-18

!! 计算SWAT手册中的净长波辐射净发射率方程2.2.20
          rbo = 0.
          rbo = -(0.34 - 0.139 * Sqrt(ed))!式子1-28

          !! 云量因子方程2.2.19
          rto = 0.
            if (hru_rmx(j) < 1.e-4) then
		    rto = 0.
            else
              rto = 0.9 * (hru_ra(j) / hru_rmx(j)) + 0.1!式子1-28
            end if

          !!长波净辐射方程2.2.21
          rout = 0.
          rout = rbo * rto * 4.9e-9 * (tk**4)!式子1-28!式子1-24

          !!计算净辐射
          rn = 0.
          rn_pet = 0.
          rn = ralb1 + rout
          rn_pet = ralb + rout
       !! 净辐射

          rho = 0.
          rho = 1710. - 6.85 * tmpav(j)!式子2-71

          if (u10(j) < 0.01) u10(j) = 0.01

        !! 潜在ET：40厘米高的参考作物苜蓿
           rv = 0.
           rc = 0.
           rv = 114. / (u10(j) * (170./1000.)**0.2)!式子2-72
           rc = 49. / (1.4 - 0.4 * co2(hru_sub(j)) / 330.)!式子2-74
           pet_day = (dlt * rn_pet + gma * rho * vpd / rv) /            
     &                               (xl * (dlt + gma * (1. + rc / rv)))!式子2-54

           pet_day = Max(0., pet_day)
 
        !! 最大装置ET
          if (igro(j) <= 0) then!判断是否有植被
            ep_max = 0.0
          else
!确定风速和风速测量高度，如有必要，调整至雨棚上方100厘米（1米）
            uzz = 0.
            zz = 0.
            if (cht(j) <= 1.0) then!冠层高度<1M
              zz = 170.0!确定风速的高度24页
            else
              zz = cht(j) * 100. + 100.!式子1-42
            end if
            uzz = u10(j) * (zz/1000.)**0.2!式子1-43

            !! 计算冠层高度（cm）
            chz = 0.
            if (cht(j) < 0.01) then
              chz = 1.!植被高度cm
            else
              chz = cht(j) * 100.
            end if

            !!计算动量传递的粗糙度长度
            zom = 0.
            if (chz <= 200.) then
              zom = 0.123 * chz!式子2-56
            else
              zom = 0.058 * chz**1.19!式子2-57
            end if
 
            !! calculate roughness length for vapor transfer计算蒸汽传输的粗糙度长度
            zov = 0.
            zov = 0.1 * zom

            !! calculate zero-plane displacement of wind profile计算风廓线的零平面位移
            d = 0.
            d = 0.667 * chz

            !! calculate aerodynamic resistance计算空气动力学阻力
            rv = Log((zz - d) / zom) * Log((zz - d) / zov)!式子2-55
            rv = rv / ((0.41)**2 * uzz)!式子2-55

            !! adjust stomatal conductivity for low vapor pressure
            !! this adjustment will lower maximum plant ET for plants
            !! sensitive to very low vapor pressure!! 调节气孔导度以适应低蒸气压这种调节将降低对极低蒸气压敏感的植物的最大ET
            xx = 0.
            fvpd = 0.
            xx = vpd - 1.
            if (xx > 0.0) then
              fvpd = Max(0.1,1.0 - vpd2(idplt(j)) * xx)
            else
              fvpd = 1.0
            end if
            gsi_adj = gsi(idplt(j)) * fvpd
            
            if (gsi_adj > 1.e-6) then
            !! calculate canopy resistance计算冠层阻力
            rc = 1. / gsi_adj                    !single leaf resistance单叶片阻力
            rc = rc / (0.5 * (laiday(j) + 0.01)                         
     &                           * (1.4 - 0.4 * co2(hru_sub(j)) / 330.))!式子2-66，2-67

            !! calculate maximum plant ET计算最大PET
            ep_max = (dlt * rn + gma * rho * vpd / rv) /                
     &                               (xl * (dlt + gma * (1. + rc / rv)))!式子2-54
            if (ep_max < 0.) ep_max = 0.
            ep_max = Min(ep_max, pet_day)
            else
              ep_max = 0.
            end if
          end if
       
       case (2)   !! HARGREAVES POTENTIAL EVAPOTRANSPIRATION METHOD哈格里夫斯潜在蒸散发法

        !! extraterrestrial radiation
        !! 37.59 is coefficient in equation 2.2.6 !!extraterrestrial
        !! 30.00 is coefficient in equation 2.2.7 !!max at surface
		!! 地外辐射
!! 37.59是方程2.2.6中的系数！！外星的
!! 30.00是方程2.2.7中的系数！！表面最大值
        ramm = 0.
        ramm = hru_rmx(j) * 37.59 / 30. 

        if (tmx(j) > tmn(j)) then
         pet_day = harg_petco(hru_sub(j))*(ramm / xl)*(tmpav(j) + 17.8)*
     &                                            (tmx(j) - tmn(j))**0.5
         pet_day = Max(0., pet_day)
        else
          pet_day = 0.
        endif
          
       case (3)  !! READ IN PET VALUES
        pet_day = petmeas
  
      end select

      return
      end