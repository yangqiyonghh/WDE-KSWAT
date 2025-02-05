      subroutine percmain
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine is the master soil percolation component.

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    drainmod tile equations   08/2006
!!	dep_imp(:)	|mm            |depth to impervious layer
!!    drainmod tile equations   08/2006

!!    icrk        |none          |crack flow code
!!                               |1 simulate crack flow in watershed
!!    inflpcp     |mm H2O        |amount of precipitation that infiltrates
!!                               |into soil (enters soil)
!!    ihru        |none          |HRU number
!!    drainmod tile equations   01/2006
!!    itdrn       |none          |tile drainage equations flag/code
!!                               |1 simulate tile flow using subroutine drains(wt_shall)
!!                               |0 simulate tile flow using subroutine origtile(wt_shall,d) 
!!    iwtdn       |none          |water table depth algorithms flag/code
!!                               |1 simulate wt_shall using subroutine new water table depth routine
!!                               |0 simulate wt_shall using subroutine original water table depth routine  
!!    ismax       |none          |maximum depressional storage selection flag/code
!!                               |1 dynamic stmaxd computed as a function of random roughness and rain intensity 
!!                               |by depstor.f
!!                               |0 static stmaxd read from .bsn for the global value or .sdr for specific hrus                 
!!    drainmod tile equations   01/2006
!!    sol_fc(:,:) |mm H2O        |amount of water available to plants in soil田间持水量下土壤/土层中植物可利用的水量 (fc-wp)
!!                               |layer at field capacity (fc - wp)
!!    sol_nly(:)  |none          |number of layers in soil profile!土壤层数
!!    sol_st(:,:) |mm H2O        |amount of water stored in the soil layer on
!!                               |the current day (less wp water)当天土层中储存的水量
!!    sol_ul(:,:) |mm H2O        |amount of water held in the soil layer at
!!                               |saturation
!!    voltot      |mm            |total volume of cracks expressed as depth
!!                               |per unit area
!!    wat_tbl(:)  |mm            |water table based on depth from soil surface
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    new water table depth  equations   01/2009
!!	c			|none		   |a factor used to convert airvol to wtd
!!	deep_p      |mm			   |total thickness of soil profile in HRU
!!    dg          |mm			   |soil layer thickness in HRU
!!    new water table depth  equations   01/2009
!!    flat(:,:)   |mm H2O        |lateral flow storage array横向流存储阵列
!!    latlyr      |mm H2O        |lateral flow in soil layer for the day当天土层中的横向流量
!!    latq(:)     |mm H2O        |total lateral flow in soil profile for the 
!!                               |day in HRUHRU中当天土壤剖面中的总横向流量
!!    lyrtile     |mm H2O        |drainage tile flow in soil layer for day排水砖在土层中的流动新的地下水位深度方程
!!    new water table depth  equations   01/2009
!!	ne_p		|mm/hr		   |effective porosity in HRU for all soil profile layers 
!!	ne_w		|mm/hr		   |effective porosity in HRU for soil layers above wtd 
!!    new water table depth  equations   01/2009
!!    qtile       |mm H2O        |drainage tile flow in soil profile for the day当天土壤剖面中的排水砖流量
!!    sepday      |mm H2O        |micropore percolation from soil layer土层微孔渗流
!!    sepbtm(:)   |mm H2O        |percolation from bottom of soil profile for
!!                               |the day in HRUHRU中当天土壤剖面底部的渗流
!!    sol_prk(:,:)|mm H2O        |percolation storage array渗流存储阵列
!!    sol_st(:,:) |mm H2O        |amount of water stored in the soil layer on
!!                               |the current day (less wp water)
!!    sol_sw(:)   |mm H2O        |amount of water stored in the soil profile
!!                               |on current day当天土壤剖面中储存的水量
!!    sw_excess   |mm H2O        |amount of water in excess of field capacity
!!                               |stored in soil layer on the current day
!!    new water table depth  equations   01/2009
!!    wat		    |mm H2O        |土壤表面以下至防渗层的浅地下水位深度
!!    new water table depth  equations   01/2009
!!    wt_shall    |mm H2O        |防渗层上方的浅水位深度
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    j1          |none          |counter
!!	w2			|mm			   |
!!	y1			|mm 		   |dummy variable for wat
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max
!!    SWAT: percmacro, percmicro, drains(wt_shall), origtile(wt_shall,d)

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: j, j1, nn, k, sb,isp, ii
      real :: lid_cuminf_total
      
      lid_cuminf_total = 0.

      j = 0
      j = ihru
      sb = inum1
      isp = isep_typ(j) 	   !! J.Jeong 6/25/14

      !! 初始化进入第一土层的水

      if (icrk == 1) then
        sepday = Max(0., inflpcp - voltot)
      else
        sepday = inflpcp
      end if

!!  增加灌溉用水
	if (aird(j)>0) then
	  j=j
	end if
      sepday = inflpcp + aird(j) + pot_seep(j)
      pot_seep(j) = 0.
      
!!如果未安装，或在此处重置地下水来源（否则为虚拟）
!!  change per JGA 10/12/11 irrigation problem with reach
!!	if (irrsc(j) > 2)  aird(j) = 0.
!!      aird(j) = 0.

      !! calculate crack flow 
      if (icrk == 1) then 
	    call percmacro
	    sepday = sepday - sepcrktot
	  endif

      do j1 = 1, sol_nly(j)
        !! 添加从上覆层进入土层的水
        sol_st(j1,j) = sol_st(j1,j) + sepday

        !!将LID改良土层（雨水花园和多孔路面）的渗滤土壤水添加到相应HRU的第一个土层
        if (j1 == 1.and.urblu(j) > 0) then
!          lid_cuminf_total = 0.
          do ii = 1, 4
            if (ii == 2) then  ! 2: rain garden, 3: cistern and 4: porous pavement
              lid_cuminf_total = lid_cuminf_total + 
     &        lid_sw_add(j,ii) * lid_farea(j,ii) * 
     &        fcimp(urblu(j)) * rg_sarea(sb,urblu(j))
            else if (ii == 3) then
              lid_cuminf_total = lid_cuminf_total + 
     &        lid_sw_add(j,ii) * lid_farea(j,ii) * 
     &        fcimp(urblu(j))
            else if (ii == 4) then
              lid_cuminf_total = lid_cuminf_total + 
     &        lid_sw_add(j,ii) * lid_farea(j,ii) * 
     &        fcimp(urblu(j))
            end if
          end do
          sol_st(j1,j) = sol_st(j1,j) + lid_cuminf_total
        end if
        
 	  !! 化粪池流入生物带层  J.Jeong
	  ! 如果土壤温度高于零度，则在生物带层中添加STE。
	  if(j1==i_sep(j).and.sol_tmp(j1,j) > 0. .and. isep_opt(j) /= 0) then
		  sol_st(j1,j) = sol_st(j1,j) + qstemm(j)  ! in mm
	    qvol = qstemm(j) * hru_ha(j) * 10.
		  xx = qvol / hru_ha(j) / 1000.
          sol_no3(j1,j) = sol_no3(j1,j) + xx *(sptno3concs(isp) 
     &                    + sptno2concs(isp))  
          sol_nh3(j1,j) = sol_nh3(j1,j) + xx * sptnh4concs(isp) 
          sol_orgn(j1,j) = sol_orgn(j1,j) + xx * sptorgnconcs(isp)*0.5
          sol_fon(j1,j) = sol_fon(j1,j) + xx * sptorgnconcs(isp) * 0.5
          sol_orgp(j1,j) = sol_orgp(j1,j) + xx * sptorgps(isp) * 0.5
          sol_fop(j1,j) = sol_fop(j1,j) + xx * sptorgps(isp) * 0.5
          sol_solp(j1,j) = sol_solp(j1,j) + xx * sptminps(isp)  
          bio_bod(j)=bio_bod(j)+xx*sptbodconcs(isp)   ! J.Jeong 4/03/09
        end if

       !! 确定分层重力排水
        sw_excess = 0.
        sw_excess = sol_st(j1,j) - sol_fc(j1,j)

        !! 初始化当前图层的变量
        sepday = 0.
        latlyr = 0.
        lyrtile = 0.
        lyrtilex = 0.

        if (sw_excess > 1.e-5) then
		!计算瓦片流量（lyrtile）、侧向流量（latlyr）和渗流（sepday）
          call percmicro(j1)

          sol_st(j1,j) = sol_st(j1,j) - sepday - latlyr - lyrtile!更新土壤含水量
          sol_st(j1,j) = Max(1.e-6,sol_st(j1,j))

          !! 如果超过田间容量（高地下水位），则重新分配土壤水分
          call sat_excess(j1)
!         sol_st(j1,j) = sol_st(j1,j) - lyrtilex
!         sol_st(j1,j) = Max(1.e-6,sol_st(j1,j))
        end if

        !! 汇总计算
        if (j1 == sol_nly(j)) then
          sepbtm(j) = sepbtm(j) + sepday
          !if (j==378)print*,sepbtm(j) ,inflpcp,sepday,i,0
        endif
        latq(j) = latq(j) + latlyr
        qtile = qtile + lyrtile
        flat(j1,j) = latlyr + lyrtile
        sol_prk(j1,j) = sol_prk(j1,j) + sepday
	  if (latq(j) < 1.e-6) latq(j) = 0.
        if (qtile < 1.e-6) qtile = 0.
        if (flat(j1,j) < 1.e-6) flat(j1,j) = 0.
      end do
          
        !! bmp adjustment
        latq(j) = latq(j) * bmp_flos(j)
        qtile = qtile * bmp_flot(j)
      
        !! 城市分布式bmps对渗流的贡献
        if (ievent > 0) then
          sepbtm(j) = sepbtm(j) + bmp_recharge(sb) 
        endif

      !! 更新土壤剖面水分
      sol_sw(j) = 0.
      do j1 = 1, sol_nly(j)
        sol_sw(j) = sol_sw(j) + sol_st(j1,j)

      end do

      !! 计算浅水表深度和瓦片流量
      qtile = 0.
      wt_shall = 0.    !CB 8/24/09
      wt_shall = dep_imp(j)
      !! drainmod tile equations排水模瓷砖方程   08/11/2006
      if (sol_tmp(2,j) > 0.) then   !Daniel 1/29/09
        por_air = 0.5
        d = dep_imp(j) - ddrain(j)
        !! drainmod wt_shall equations  排水模型wt_shall方程 10/23/2006
        if (iwtdn == 0) then !compute wt_shall using original eq-Daniel使用原始方程Daniel计算wt_shall 10/23/06
          if (sol_sw(j) > sol_sumfc(j)) then
            yy = sol_sumul(j) * por_air
            if (yy < 1.1 * sol_sumfc(j)) then
              yy = 1.1 * sol_sumfc(j)
            end if
            xx = (sol_sw(j) - sol_sumfc(j)) / (yy - sol_sumfc(j))
            if (xx > 1.) xx = 1.
            wt_shall = xx * dep_imp(j)
		    wat = dep_imp(j) - wt_shall
			if(wat > dep_imp(j)) wat = dep_imp(j)
          end if
        else
          !compute water table depth using Daniel's modifications
   !       Updated water table depth D.Moriasi 4/8/2014使用Daniel的修改计算地下水位深度更新的地下水位深度D.Moriasi 2014年8月4日
           swst_del = 0.
           sw_del = 0.
           wt_del = 0.
           wtst_del = 0.
          do j1 = 1, sol_nly(j)
!            if (wat_tbl(j) < sol_z(j1,j)) then
             swst_del = sol_stpwt(j1,j) - sol_st(j1,j)
             sw_del = sol_swpwt(j) - sol_sw(j)
             wt_del = sw_del * vwt(j1,j)
             wtst_del = swst_del * vwt(j1,j)
 !            wat_tbl(j) = wat_tbl(j) + wt_del 
             wat_tbl(j) = wat_tbl(j) + wtst_del  
             if(wat_tbl(j) < 0.0) wat_tbl(j) = 0.0
	       if(wat_tbl(j) > dep_imp(j)) wat_tbl(j) = dep_imp(j)
	       wt_shall = dep_imp(j) - wat_tbl(j)
	       sol_swpwt(j) = sol_sw(j)
	       sol_stpwt(j1,j) = sol_st(j1,j)
!	        exit
!	      end if
!       更新的地下水位深度D.Moriasi 4/8/2014
	    end do
        end if
        !! drainmod wt_shall equations   排水模型wt_shall方程10/23/2006
        
        if (ddrain(j) > 0.) then
          if (wt_shall <= d) then
            qtile = 0.
          else
            !! Start Daniel's tile equations modifications 开始修改Daniel的tile方程 01/2006
            if (itdrn == 1) then
              call drains     ! compute tile flow using drainmod tile equations 使用draimod瓷砖方程计算瓷砖流量
              !! drainmod tile equations  排水模瓷砖方程 01/2006
            else !! compute tile flow using existing tile equations使用现有瓦片方程计算瓦片流量
              call origtile(d)! existing tile equations 现有平铺方程式
	        if(qtile < 0.) qtile=0.
            end if 
          end if
        end if
      end if
      !! End Daniel's tile equations modifications 结束Daniel的tile方程修改 01/2006

      if (qtile > 0.) then
        !! update soil profile water after tile drainage瓷砖排水后更新土壤剖面水
        sumqtile = qtile
        do j1 = 1, sol_nly(j)
          xx = sol_st(j1,j) - sol_fc(j1,j)
          if (xx > 0.) then
            if (xx > sumqtile) then
              sol_st(j1,j) = sol_st(j1,j) - sumqtile
              sumqtile = 0.
            else
              sumqtile = sumqtile - xx
              sol_st(j1,j) = sol_fc(j1,j)
            end if
          end if
        end do
        if (sumqtile > 0.) then
          qtile = qtile - sumqtile
          qtile = amax1(0., qtile)
        end if
      end if

      !! update soil profile water更新土壤剖面水分
      sol_sw(j) = 0.
      do j1 = 1, sol_nly(j)
        sol_sw(j) = sol_sw(j) + sol_st(j1,j)
      end do

      return
      end