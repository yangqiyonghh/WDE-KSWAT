      subroutine swu
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine distributes potential plant evaporation through
!!    the root zone and calculates actual plant water use based on soil
!!    water availability. Also estimates water stress factor.     
	 !这个子程序通过根区分配潜在的植物蒸发量，并根据土壤水分的可用性计算实际的植物用水量。还估算了水分胁迫因子
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ep_max      |mm H2O        |maximum amount of transpiration (plant et)
!!                               |that can occur on current day in HRU
!!    epco(:)     |none          |plant water uptake compensation factor (0-1)
!!    icr(:)      |none          |sequence number of crop grown within the
!!                               |current year
!!    idc(:)      |none          |crop/landcover category:
!!                               |1 warm season annual legume
!!                               |2 cold season annual legume
!!                               |3 perennial legume
!!                               |4 warm season annual
!!                               |5 cold season annual
!!                               |6 perennial
!!                               |7 trees
!!    idplt(:)    |none          |land cover code from crop.dat
!!    ihru        |none          |HRU number
!!    iwatable    |none          |high water table code:
!!                               |0 no high water table
!!                               |1 high water table
!!    nro(:)      |none          |sequence number of year in rotation
!!    phuacc(:)   |none          |fraction of plant heat units accumulated工厂累计热量单位比例
!!    sol_fc(:,:) |mm H2O        |amount of water available to plants in soil
!!                               |layer at field capacity (fc - wp water)
!!    sol_nly(:)  |none          |number of soil layers in profile
!!    sol_st(:,:) |mm H2O        |amount of water stored in the soil layer on
!!                               |current day
!!    sol_ul(:,:) |mm H2O        |amount of water held in the soil layer at
!!                               |saturation
!!    sol_z(:,:)  |mm            |depth to bottom of soil layer
!!    sol_zmx(:)  |mm            |maximum rooting depth
!!    stsol_rd(:) |mm            |storing last soil root depth for use in harvestkillop/killop储存最后的土壤根深以用于收获
!!    ubw         |none          |water uptake distribution parameter
!!                               |This parameter controls the amount of
!!                               |water removed from the different soil layers
!!                               |by the plant. In particular, this parameter
!!                               |allows the amount of water removed from
!!                               |the surface layer via plant uptake to be
!!                               |controlled. While the relationship between
!!                               |UBW and H2O removed from the surface layer is
!!                               |affected by the depth of the soil profile, in
!!                               |general, as UBW increases the amount of water
!!                               |removed from the surface layer relative to the
!!                               |amount removed from the entire profile
!!                               |increases
!!    uobw        |none          |water uptake normalization parameter
!!                               |This variable normalizes the water uptake so
!!                               |that the model can easily verify that uptake
!!                               |from the different soil layers sums to 1.0
									!吸水归一化参数这个变量使水的吸收标准化，因此模型可以很容易地验证从
									!不同土壤层的吸收总和为1.0
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ep_day      |mm H2O        |actual amount of transpiration that occurs
!!                               |on day in HRU
!!    sol_rd      |mm            |current rooting depth
!!    sol_st(:,:) |mm H2O        |amount of water stored in the soil layer on
!!                               |current day
!!    sol_sw(:)   |mm H2O        |amount of water stored in soil profile on
!!                               |current day
!!    strsw(:)    |none          |fraction of potential plant growth achieved
!!                               |on the day where the reduction is caused by
!!                               |water stress水分胁迫导致植物生长减少的当天的潜在植物生长比例
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    gx          |
!!    ir          |
!!    j           |none          |HRU number
!!    k           |none          |counter (soil layer)
!!    reduc       |none          |fraction of water uptake by plants achieved
!!                               |where the reduction is caused by low water
!!                               |content
!!    sum         |
!!    sump        |
!!    wuse(:)     |mm H2O        |water uptake by plants in each soil layer土壤各层植物的水分吸收
!!    xx          |mm H2O        |water uptake by plants from all layers各层植物对水分的吸收
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Max

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: j, k, ir
      real, dimension(mlyr) :: wuse
      real :: sum, xx, gx, reduc, sump

      j = 0
      j = ihru

      select case (idc(idplt(j)))
        case (1, 2, 4, 5)
          sol_rd = 2.5 * phuacc(j) * sol_zmx(j)!式子5-31
          if (sol_rd > sol_zmx(j)) sol_rd = sol_zmx(j)
          if (sol_rd < 10.) sol_rd = 10.
        case default
          sol_rd = sol_zmx(j)
      end select

	  stsol_rd(j) = sol_rd ! cole armen 26 Feb

      if (ep_max <= 0.01) then
        strsw(j) = 1.
      else
        !! initialize variables
        gx = 0.
        ir = 0
        sump = 0.
        wuse = 0.
        xx = 0.
 
!!  compute aeration stress计算通气应力
        if (sol_sw(j) > sol_sumfc(j)) then
          satco = (sol_sw(j) - sol_sumfc(j)) / (sol_sumul(j) - 
     &                                                 sol_sumfc(j))
          pl_aerfac = .85
          scparm = 100. * (satco - pl_aerfac) / (1.0001 - pl_aerfac)
          if (scparm > 0.) then
            strsa(j) = 1. - (scparm / (scparm + Exp(2.9014 - .03867 *
     &                                                      scparm)))
          else
            strsa(j) = 1.
          end if
        end if

        do k = 1, sol_nly(j)
          if (ir > 0) exit

          if (sol_rd <= sol_z(k,j)) then
            gx = sol_rd
            ir = k
          else
            gx = sol_z(k,j)
          end if

          sum = 0.
          if (sol_rd <= 0.01) then
            sum = ep_max / uobw
          else
            sum = ep_max * (1. - Exp(-ubw * gx / sol_rd)) / uobw!5-33
          end if

          !! don't allow compensation for aeration stress不允许对通气压力进行补偿
!          if (strsa(j) > .99) then
!           yy = 0.
!          else
!            yy= sump - xx
!          end if
          wuse(k) = sum - sump + 1. * epco(j)!5-35
          wuse(k) = sum - sump + (sump - xx) * epco(j)!5-35
          sump = sum


          !! adjust uptake if sw is less than 25% of plant available water如果SW小于植物可用水量的25%，调整吸收
          reduc = 0.
          if (sol_st(k,j) < sol_fc(k,j)/4.) then
            reduc = Exp(5. * (4. * sol_st(k,j) / sol_fc(k,j) - 1.))!式子5-36
          else
            reduc = 1.
          endif
          wuse(k) = wuse(k) * reduc

          if (sol_st(k,j) < wuse(k)) then
            wuse(k) = sol_st(k,j)
          end if

          sol_st(k,j) = Max(1.e-6, sol_st(k,j) - wuse(k))
          xx = xx + wuse(k)
          end do

        !! update total soil water in profile更新剖面土壤总水量
        sol_sw(j) = 0.
        do k = 1, sol_nly(j)
          sol_sw(j) = sol_sw(j) + sol_st(k,j)
        end do

        strsw(j) = xx / ep_max
        ep_day = xx
      end if

      return
      end