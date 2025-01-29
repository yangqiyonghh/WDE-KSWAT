      subroutine gwmod
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine estimates groundwater contribution to
!!    streamflow

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    alpha_bf(:) |1/days        |alpha factor for groundwater recession curve
!!    alpha_bfe(:)|none          |Exp(-alpha_bf(:))
!!    deepst(:)   |mm H2O        |depth of water in deep aquifer
!!    ihru        |none          |HRU number
!!    gw_delaye(:)|none          |Exp(-1./(delay(:)) where delay(:) is the 
!!                               |groundwater delay (time required for water
!!                               |leaving the bottom of the root zone to reach
!!                               |the shallow aquifer; units-days)
!!    gw_revap(:) |none          |revap coeff: this variable controls the amount
!!                               |of water moving from the shallow aquifer to
!!                               |the root zone as a result of soil moisture
!!                               |depletion
!!    gw_spyld(:) |m**3/m**3     |specific yield for shallow aquifer
!!    gwht(:)     |m             |groundwater height
!!    gwqmn(:)    |mm H2O        |threshold depth of water in shallow aquifer
!!                               |required before groundwater flow will occur
!!    pet_day     |mm H2O        |potential evapotranspiration on current day
!!                               |in HRU
!!    rchrg(:)    |mm H2O        |amount of water entering shallow aquifer on
!!                               |previous day in HRU
!!    rchrg_dp(:) |none          |recharge to deep aquifer: the fraction of
!!                               |root zone percolation that reaches the deep向深层含水层的补给：到达深层含水层根区渗流的分数
!!                               |aquifer
!!    revapmn(:)  |mm H2O        |threshold depth of water in shallow aquifer
!!                               |required to allow revap to occur
!!    sepbtm(:)   |mm H2O        |percolation from bottom of soil profile for
!!                               |the day in HRU
!!    shallst(:)  |mm H2O        |depth of water in shallow aquifer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    deepst(:)   |mm H2O        |depth of water in deep aquifer
!!    gw_q(:)     |mm H2O        |groundwater contribution to streamflow from
!!                               |HRU on current day
!!    gwht(:)     |m             |groundwater height
!!    gwseep      |mm H2O        |amount of water recharging deep aquifer on
!!                               |current day in HRUHRU中当天深层含水层的补给水量
!!    rchrg(:)    |mm H2O        |amount of water recharging both aquifers on
!!                               |current day in HRUHRU当日补给两含水层的水量
!!    revapday    |mm H2O        |amount of water moving from the shallow 
!!                               |aquifer into the soil profile or being taken
!!                               |up by plant roots in the shallow aquifer
!!    shallst(:)  |mm H2O        |depth of water in shallow aquifer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    rchrg1      |mm H2O        |amount of water entering shallow aquifer on
!!                               |previous day前一天进入浅层含水层的水量
!!    rchrg_karst |mm H2O        |amount of water from secondary channels,
!!                               |ponds, and wetlands
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
!!    revap is subtracted and rchrg is delayed (johnson, 1977)减去revap，延迟rchrg

      use parm

      integer :: j
      real :: rchrg1, rchrg_karst
      real :: xx, r2
      j = 0
      j = ihru

      rchrg1 = 0.
      rchrg_karst = 0.
      rchrg1 = rchrg(j) + rchrg_src(j)
      if (curyr==1 .and. i==1 .and. karst_hru(j)>=1)shallst(j)=karst_hs!初始水深
      
!! add seepage from secondary channels, ponds, and wetlands;增加二级渠道、池塘和湿地的渗流；
      rchrg_karst = tloss + twlpnd + twlwet
!! compute shallow aquifer level for current day, assumes karst losses 
!! infiltrate at the same speed as what goes through the soil profile.计算当天的浅层含水层水位，假设岩溶损失以与土壤剖面相同的速度渗透。
      if (karst_hru(j)>=1)sepbtm(j) =Q_karst_S(j)
      if (karst_hru(j)>=1)gw_delaye(j)=Exp(-gw_karst_delaye)
      if (karst_hru(j)>=1)alpha_bfe(j)=Exp(-gw_alpha_bfe)
      rchrg(j) = 0.
      rchrg(j) = (1.-gw_delaye(j)) * (sepbtm(j) + gwq_ru(j) +           
     &                             rchrg_karst) + gw_delaye(j) * rchrg1!式子2-132
      if (rchrg(j) < 1.e-6) rchrg(j) = 0.
      gwq_ru(j) = 0.
      
!! compute deep aquifer level for day计算一天的深层含水层水位
      gwseep = rchrg(j) * rchrg_dp(j)!式子2-134
      deepst(j) = deepst(j) + gwseep

      shallst(j) = shallst(j) + (rchrg(j) - gwseep)!式子2-135
      gwht(j) = gwht(j) * alpha_bfe(j) + rchrg(j) * (1. - alpha_bfe(j)) !式子2-150
     &    / (800. * gw_spyld(j) * alpha_bf(j) + 1.e-6)
      gwht(j) = Max(1.e-6, gwht(j))

      
      if (gw_karst_hru(inum1)==1 .and. karst_hru(j)>=1)then
          
          if (curyr==1 .and. i==1)karst_FM_height(j)=karst_hf!初始水深
          karst_FM_height(j)=karst_FM_height(j)+Q_karst_F(j)          
          if (karst_FM_height(j)>shallst(j) )then!地下管道向基质流动
          
              karst_sc(j)=max(k_sc*(karst_FM_height(j)-shallst(j) )**alpha_sc,0.010)
              !karst_sc(j)=min(karst_sc(j),0.5*(karst_FM_height(j)-shallst(j)))
              if (karst_sc(j)>karst_FM_height(j))then
                  karst_sc(j)=karst_FM_height(j)
                  shallst(j)=shallst(j)+karst_sc(j)
              else
                  shallst(j)=shallst(j)+karst_sc(j)
                  
              endif
              Q_karst_F1(j)=max(Q_karst_F(j)-karst_sc(j),0.0)
              Q_karst_S1(j)=max(Q_karst_S(j)+karst_sc(j),0.0)              
              karst_sc1(j)=-1.0*karst_sc(j)
              !if (karst_sc1(j)>100 .or. karst_sc1(j)<-100)print*,karst_sc(j),karst_sc1(j),1
              karst_FM_height(j)=max(karst_FM_height(j)-karst_sc(j),0.0)
          else!基质向地下管道流动
              !if ( j==172)print*,k_sc,alpha_sc,shallst(j),i,1
              karst_sc(j)=max(k_sc*(shallst(j) -karst_FM_height(j))**alpha_sc,0.01)
              !karst_sc(j)=min(karst_sc(j),0.5*(shallst(j)-karst_FM_height(j)))
              if (karst_sc(j)>shallst(j))then
                  shallst(j)=0.
                  karst_sc(j)=shallst(j)
                  
              else
                  shallst(j)=shallst(j)-karst_sc(j)
                  
              endif
              Q_karst_F1(j)=max(Q_karst_F(j)+karst_sc(j),0.0)
              Q_karst_S1(j)=max(Q_karst_S(j)-karst_sc(j),0.0)
              karst_sc1(j)=1.0*karst_sc(j)
              !if (karst_sc1(j)>100 .or. karst_sc1(j)<-100)print*,karst_sc(j),karst_sc1(j)
              karst_FM_height(j)=karst_FM_height(j)+karst_sc(j)
              
          endif 

              Q_karst_FM(j)=k_fm*(karst_FM_height(j))**alpha_FM!管道出流量
              Q_karst_FM(j)=min(Q_karst_FM(j),karst_FM_height(j))
              !if ( j==172)print*,Q_karst_FM(j),karst_FM_height(j),karst_sc(j),i
              if (karst_FM_height(j) >= gw_karst_h2) then
                karst_FM_height(j) = karst_FM_height(j) - Q_karst_FM(j)
                if (karst_FM_height(j)< gw_karst_h2) then
                 Q_karst_FM(j)= karst_FM_height(j) +Q_karst_FM(j) -gw_karst_h2
                   karst_FM_height(j) = gw_karst_h2
                end if
               else
                Q_karst_FM(j) = 0.
               end if   
      endif
      
      
      
!! compute groundwater contribution to streamflow for day计算地下水对一天流量的贡献
      if (karst_hru(j)==0)then
      if (shallst(j) > gwqmn(j)) then
        gw_q(j) = gw_q(j) * alpha_bfe(j) + (rchrg(j) - gwseep ) *       
     &                                               (1. - alpha_bfe(j))
      else
        gw_q(j) = 0.
      end if
      end if
      if (karst_hru(j)>=1)gw_q(j) = gw_q(j) * alpha_bfe(j) + (rchrg(j) - gwseep ) *       
     &                                               (1. - alpha_bfe(j))
      !! bmp adjustmentbmp调整
      gw_q(j) = gw_q(j) * bmp_flos(j)

!! compute revap to soil profile/plant roots计算土壤剖面/植物根系的revap
      revapday = gw_revap(j) * pet_day!式子2-145
      if (shallst(j) < revapmn(j)) then
        revapday = 0.!式子2-146
      else
        shallst(j) = shallst(j) - revapday
        if (shallst(j) < revapmn(j)) then
          revapday = shallst(j) + revapday - revapmn(j)!式子2-147
          shallst(j) = revapmn(j)
        end if
      end if

!! remove ground water flow from shallow aquifer storage清除浅层蓄水层中的地下水
      if (karst_hru(j)==0)then!非岩溶
      if (shallst(j) >= gwqmn(j)) then
        shallst(j) = shallst(j) - gw_q(j)
        if (shallst(j) < gwqmn(j)) then
          gw_q(j) = shallst(j) + gw_q(j) - gwqmn(j)
          shallst(j) = gwqmn(j)
        end if
       else
        gw_q(j) = 0.
       end if
       
      else!岩溶基质

      if (shallst(j) >= gw_karst_h1) then
        shallst(j) = shallst(j) - gw_q(j)
        if (shallst(j) < gw_karst_h1) then
          gw_q(j) = shallst(j) + gw_q(j) -gw_karst_h1
          shallst(j) = gw_karst_h1
        end if
       else
        gw_q(j) = 0.
       end if          
          
       end if
       
      return
      end