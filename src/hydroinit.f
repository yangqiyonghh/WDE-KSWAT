      subroutine hydroinit 

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    该子程序计算与流域水文相关的变量：子流域的集中时间、滞后地表径流、峰值径流率方程的系数和横向流行进时间。

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~1
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ch_l1(:)    |km            |longest tributary channel length in subbasin
!!    ch_l2(:)    |km            |main channel length in subbasin
!!    ch_n(1,:)   |none          |Manning's "n" value for the tributary channels
!!    ch_s(1,:)   |m/m           |average slope of tributary channels
!!    da_km       |km2           |area of the watershed in square kilometers
!!    gdrain(:)   |hrs           |drain tile lag time: the amount of time
!!                               |between the transfer of water from the soil
!!                               |to the drain tile and the release of the
!!                               |water from the drain tile to the reach.
!!    hru_dafr(:) |km2/km2       |fraction of total watershed area contained
!!                               |in HRU
!!    hru_km(:)   |km2           |area of HRU in square kilometers
!!    hru_slp(:)  |m/m           |average slope steepness
!!    hru_sub(:)  |none          |subbasin in which HRU is located
!!    lat_ttime(:)|days          |lateral flow travel time
!!   tile_ttime(:)|none          |Exponential of the tile flow travel time
!!    ldrain(:)   |none          |soil layer where drainage tile is located排水砖所在的土层
!!    nhru        |none          |number of HRUs in watershed
!!    ov_n(:)     |none          |Manning's "n" value for overland flow
!!    slsoil(:)   |m             |slope length for lateral subsurface flow
!!    slsubbsn(:) |m             |average slope length for subbasin
!!    sol_k(:,:)  |mm/hr         |saturated hydraulic conductivity of soil
!!                               |layer
!!    sol_nly(:)  |none          |number of layers in soil profile
!!    sub_fr(:)   |none          |fraction of total watershed area contained in
!!                               |subbasin
!!    surlag      |days          |Surface runoff lag time.
!!                               |This parameter is needed in subbasins where
!!                               |the time of concentration is greater than 1 
!!                               |day. SURLAG is used to create a "storage" for
!!                               |surface runoff to allow the runoff to take 
!!                               |longer than 1 day to reach the subbasin outlet
!!    tconc(:)     |hr           |time of concentration
!!    uhalpha      |             |alpha coefficient for estimating unit hydrograph
!!                               |using a gamma function (*.bsn)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 


!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    brt(:)      |none          |fraction of surface runoff within the subbasin
!!                               |which takes 1 day or less to reach the 
!!                               |subbasin outlet
!!    lat_ttime(:)|none          |Exponential of the lateral flow travel time
!!   tile_ttime(:)|none          |Exponential of the tile flow travel time
!!    sub_tc(:)   |hr            |time of concentration for subbasin次盆地集中时间
!!    t_ov(:)     |hr            |time for flow from farthest point in subbasin
!!                               |to enter a channel从子盆地最远点流入河道的时间
!!    tconc(:)    |hr            |time of concentration for HRUHRU集中时间
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    j           |none          |counter
!!    l           |none          |counter
!!    scmx        |mm/hr         |maximum soil hydraulic conductivity
!!    t_ch        |hr            |time for flow entering the farthest upstream 水流进入最远上游河道到达底池出水口的时间
!!                               |channel to reach the subbasin outlet
!!    xx          |none          |variable to hold calculation result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~  
!!    SWAT: Ttcoef

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


      use parm

      integer :: j, l
      real :: t_ch, scmx, xx

      do j = 1, nhru

!! subbasin !!
!!   计算集中时间（陆上和航道时间之和）
        t_ch = 0
        t_ov(j) = .0556 * (slsubbsn(j)*ov_n(j)) ** .6 / hru_slp(j) ** .3!式子2-25
        t_ch = .62 * ch_l1(j) * ch_n(1,hru_sub(j)) ** .75 /             
     &      ((da_km * sub_fr(hru_sub(j)))**.125 *                       
     &                                         ch_s(1,hru_sub(j))**.375)!式子2-33
        sub_tc(hru_sub(j)) = t_ov(j) + t_ch
!! end subbasin !!


!! HRU !!
!!    计算集中时间（陆上和航道时间之和）
        t_ch = 0
        ch_l1(j) = ch_l1(j) * hru_dafr(j) / sub_fr(hru_sub(j))
        t_ov(j) = .0556 * (slsubbsn(j)*ov_n(j)) ** .6 / hru_slp(j) ** .3
        t_ch = .62 * ch_l1(j) * ch_n(1,hru_sub(j)) ** .75 /             
     &              ((da_km*hru_dafr(j))**.125*ch_s(1,hru_sub(j))**.375)
        tconc(j) = t_ov(j) + t_ch

          surlag_option=1
          if (surlag_option==1)then
!修改---------------------hru-----------------------------------------------------------------------------------
        t_ov(j) = .2525 * (slsubbsn(j)*(ov_n(j)) ** 2 )/ hru_slp(j) !式子2-25
        t_ch = 0.33* ch_l1(j) *(slsubbsn(j)**0.25)*(ov_n(j)**0.874)* !式子2-33            
     &     (ch_n(1,hru_sub(j)) ** .75) /(((da_km*hru_dafr(j))**.125)*
     &             (ch_s(1,hru_sub(j))**.375)*(hru_slp(j) **0.437))
         tconc(j) = t_ov(j) + t_ch
         !print*,tconc(j),j,22
!修改-------------------------hru-----------------------------------------------------------------------------        
!修改surlag============================================================
        v_ov =0
        v_ch =0
        v_minutes =0
        v_ov =0.0011*hru_slp(j)/(ov_n(j)**2)
        v_ch =0.594*(hru_slp(j) **0.437)*( (da_km*hru_dafr(j))**.125)*
     &       (ch_s(1,hru_sub(j))**.375)/((slsubbsn(j)**0.25)*
     &                  (ov_n(j)**0.874)*(ch_n(1,hru_sub(j)) ** .75))
        v_minutes =(v_ov+v_ch)/2
        surlag(j) =0.87* ch_l1(j)/v_minutes

!修改surlag============================================================     
        end if        
        
        
        
!!    compute delivery ratio计算交付比率
      rto = tconc(j) / sub_tc(hru_sub(j))
      dr_sub(j) = amin1(.95,rto ** .5)


!!    compute fraction of surface runoff that is reaching the main channel计算到达主河道的地表径流的分数
        if (ievent > 0) then
           brt(j) = 1. - Exp(-surlag(j) / (tconc(j) / (idt / 60.)))	!! urban modeling by J.Jeong
        else
           brt(j) = 1. - Exp(-surlag(j) / tconc(j))
        endif
        !PRINT*, brt(j) 
        if (brt(j)<0.3) brt(j) =0.3
        if (isproj == 2) brt(j) = 1.
        
        
!!    compute lateral flow travel time计算横向流行程时间
        if (lat_ttime(j) <= 0.) then
            scmx = 0.
            do l = 1, sol_nly(j)
              if (sol_k(l,j) > scmx) then
                scmx = sol_k(l,j)
              endif
            end do
            !! unit conversion:
            !! xx = m/(mm/h) * 1000.(mm/m)/24.(h/d) / 4.
            xx = 0.
            xx = 10.4 * slsoil(j) / scmx
            lat_ttime(j) = 1. - Exp(-1./xx)
        else
          lat_ttime(j) = 1. - Exp(-1./lat_ttime(j))
        end if

        if (ldrain(j) > 0 .and. gdrain(j) > 0.01) then
            tile_ttime(j) = 1. - Exp(-24. / gdrain(j))
        else
            tile_ttime(j) = 0.
        end if

!!    compute routing coefficients for main channel计算主信道的路由系数
        if (ch_l2(hru_sub(j)) > 0.) call ttcoef(hru_sub(j))
        if (j == hru1(hru_sub(j))) then
          if (alpha_bnk(hru_sub(j)) <= 0.) then
            alpha_bnk(hru_sub(j)) = alpha_bf(j)
          end if
          alpha_bnke(hru_sub(j)) = Exp(-alpha_bnk(hru_sub(j)))
        end if

      end do

      if (ievent > 0) then
!!    compute unit hydrograph for computing subbasin hydrograph from direct runoff
!计算单位过程线，用于根据直接径流计算次流域过程线
      do isb = 1, msub 
        ql = 0.
        sumq = 0.

        tb = .5 + .6 * sub_tc(isb) + 0.0  !baseflow time, hr基本流量时间，先生

        if (tb > 48.) tb = 48.			   !maximum 48hrs
        tp = .375 * tb						! time to peak flow达到峰值流量的时间
	  
	  !! convert to time step (from hr)转换为时间步长（从小时开始）, J.Jeong March 2009
	  tb = ceiling(tb * 60./ real(idt))
	  tp = int(tp * 60./ real(idt))         
	  
	  if(tp==0) tp = 1
	  if(tb==tp) tb = tb + 1
	  itb(isb) = int(tb) 
        
	  ! Triangular Unit Hydrograph三角形单位过程线
	  if (iuh==1) then
	    do i = 1, itb(isb)
            xi = float(i)
 	      if (xi < tp) then           !! rising limb of hydrograph过程线上升段
              q = xi / tp
            else                        !! falling limb of hydrograph过程线下降段
              q = (tb - xi) / (tb - tp)
            end if
            q = Max(0.,q)
            uh(isb,i) = (q + ql) / 2.
            ql = q
            sumq = sumq + uh(isb,i)
          end do
          
		do i = 1, itb(isb)
            uh(isb,i) = uh(isb,i) / sumq
          end do
	  
	  ! Gamma Function Unit Hydrograph伽马函数单位过程线
	  elseif (iuh==2) then
          i = 1; q=1.
		do while (q>0.0001)
            xi = float(i)
		   q = (xi / tp) ** uhalpha * exp((1.- xi / tp) * uhalpha)
            q = Max(0.,q)
            uh(isb,i) = (q + ql) / 2.
            ql = q
            sumq = sumq + uh(isb,i)
	      i = i + 1
	      if (i>3.*nstep) exit
	    end do
	    itb(isb) = i - 1
		  do i = 1, itb(isb)
            uh(isb,i) = uh(isb,i) / sumq
          end do
	  endif 

      end do
      end if

      return
      end