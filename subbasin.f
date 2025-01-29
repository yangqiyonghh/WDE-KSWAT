      subroutine subbasin
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine controls the simulation of the land phase of the 
!!    hydrologic cycle

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name           |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    auto_wstr(:)   |none          |water stress factor which triggers auto
!!                                  |irrigation
!!    bio_e(:)       |(kg/ha)/      |biomass-energy ratio
!!                   |     (MJ/m**2)|The potential (unstressed) growth rate per
!!                                  |unit of intercepted photosynthetically
!!                                  |active radiation.
!!    canev          |mm H2O        |amount of water evaporated from canopy
!!                                  |storage
!!    ep_day         |mm H2O        |actual amount of transpiration that occurs
!!                                  |on day in HRUHRU中当天发生的实际蒸腾量
!!    es_day         |mm H2O        |actual amount of evaporation (soil et) that
!!                                  |occurs on day in HRUHRU中当天发生的实际蒸发量（土壤et）
!!    gw_q(:)        |mm H2O        |groundwater contribution to streamflow from
!!                                  |HRU on current day
!!    hru_ra(:)      |MJ/m^2        |solar radiation for the day in HRU
!!    iida           |julian date   |day being simulated (current julian date)
!!    idplt(:)       |none          |land cover code from crop.dat
!!    igro(:)        |none          |land cover status code
!!                                  |0 no land cover currently growing
!!                                  |1 land cover growing
!!    inum1          |none          |subbasin number
!!    imp_trig(:)    |none          |release/impound action code:
!!                                  |0 begin impounding water
!!                                  |1 release impounded water
!!    irrsc(:)       |none          |irrigation source code:
!!                                  |1 divert water from reach
!!                                  |2 divert water from reservoir
!!                                  |3 divert water from shallow aquifer
!!                                  |4 divert water from deep aquifer
!!                                  |5 divert water from source outside
!!                                  |  watershed
!!    iurban(:)      |none          |urban simulation code:
!!                                  |0  no urban sections in HRU
!!                                  |1  urban sections in HRU, simulate using
!!                                  |   USGS regression equations
!!                                  |2  urban sections in HRU, simulate using
!!                                  |   build up/wash off algorithm
!!    latq(:)        |mm H2O        |total lateral flow in soil profile for the
!!                                  |day in HRU
!!    nafert(:)      |none          |sequence number of auto-fert application
!!                                  |within the year
!!    nair(:)        |none          |sequence number of auto-irrigation 
!!                                  |application within the year
!!    nfert(:)       |none          |sequence number of fertilizer application
!!                                  |within the year
!!    nirr(:)        |none          |sequence number of irrigation application
!!                                  |within the year
!!    nrelease(:)    |none          |sequence number of impound/release
!!                                  |operation within the year
!!    nro(:)         |none          |sequence number of year in rotation
!!    peakr          |m^3/s         |peak runoff rate
!!    pet_day        |mm H2O        |potential evapotranspiration on current
!!                                  |day in HRU
!!    phuacc(:)      |none          |fraction of plant heat units accumulated
!!    phubase(:)     |heat units    |base zero total heat units (used when no
!!                                  |land cover is growing)
!!                                  |pesticide application occurs
!!    pot_fr(:)      |km2/km2       |fraction of HRU area that drains into
!!                                  |pothole
!!    pot_vol(:)     |m**3 H2O      |current volume of water stored in the 
!!                                  |depression/impounded area
!!    precipday      |mm H2O        |precipitation for the day in HRU
!!    qday           |mm H2O        |surface runoff loading to main channel from
!!                                  |HRU for day
!!    qtile          |mm H2O        |drainage tile flow in soil layer for the 
!!                                  |day
!!    sci(:)         |none          |retention coefficient for CN method based
!!                                  |on plant ET
!!    sedyld(:)      |metric tons   |soil loss for day in HRU
!!    smx(:)         |none          |retention coefficient for CN method based
!!                                  |on soil moisture
!!    surfq(:)       |mm H2O        |surface runoff generated on day in HRU
!!    tmn(:)         |deg C         |minimum temperature for the day in HRU
!!    tmpav(:)       |deg C         |average temperature for the day in HRU
!!    tmx(:)         |deg C         |maximum temperature for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    albday      |none          |albedo, the fraction of the solar radiation
!!                               |reflected at the soil surface back into
!!                               |space
!!    etday       |mm H2O        |actual evapotranspiration occuring on day
!!                               |in HRU
!!    ihru        |none          |HRU number
!!    inflpcp     |mm H2O        |amount of precipitation that infiltrates
!!                               |into soil (enters soil)
!!    nafert(:)   |none          |sequence number of auto-fert application
!!                               |within the year
!!    nair(:)     |none          |sequence number of auto-irrigation 
!!                               |application within the year
!!    qdfr        |none          |fraction of water yield that is surface
!!                               |runoff
!!    qdr(:)      |mm H2O        |total amount of water entering main channel
!!                               |for day from HRU从HRU进入主通道的总水量
!!    sci(:)      |none          |retention coefficient for CN method based
!!                               |on plant ET
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    d           |
!!    gma         |kPa/deg C     |psychrometric constant
!!    ho          |              |net radiation
!!    j           |none          |HRU number
!!    pet_alpha   |none          |alpha factor in Priestley-Taylor ET 
!!                               |equation
!!    tmpk        |deg K         |average temperature for the day in the HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Max
!!    SWAT: varinit, albedo, solt, surface, percmain, etpot, etact, fert
!!    SWAT: confert, graze, plantmod, nminrl, nitvol, pminrl, gwmod, apply, gwmod_deep
!!    SWAT: washp, decay, pestlch, enrsb, pesty, orgn, psed, nrain, nlch
!!    SWAT: solp, subwq, bacteria, urban, pothole, latsed, surfstor
!!    SWAT: substor, wetland, hrupond, irrsub, autoirr, watuse, watbal
!!    SWAT: sumv, virtual

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: j,sb,kk
      real :: tmpk, d, gma, ho, pet_alpha, aphu, phuop,obs,sim,wapi_all1

      ihru = 0
      ihru = hru1(inum1) 
      
      call sub_subbasin
      
      do iihru = 1, hrutot(inum1)
      
      j = 0
      j = ihru
      
      !if (station_flag/=1 )then
      !if (alpha_delay1/=sftmp)print*,alpha_delay1,sftmp,1,i,j
      !if (alpha_delay2/=smtmp)print*,alpha_delay2,smtmp,2,i,j
      !if (sur_karst_h1/=smfmn)print*,sur_karst_h1,smfmn,3,i,j
      !if (sur_karst_h2/=smfmn+smfmx)print*,sur_karst_h2,smfmn+smfmx,3,i,j
      !if (k_sec/=timp)print*,k_sec,timp,4,i,j
      !if (k_eh/=snocovmx)print*,k_eh,snocovmx,5,i,j
      !if (alpha_eh/=sno50cov)print*,alpha_eh,sno50cov,6,i,j
      !if (k_f/=spcon_bsn)print*,k_f,spcon_bsn,7,i,j
      !if (gw_karst_delaye/=spexp_bsn)print*,gw_karst_delaye,spexp_bsn,8,i,j
      !if (gw_alpha_bfe/=n_updis)print*,gw_alpha_bfe,n_updis,9,i,j
      !if (k_fm/=p_updis)print*,k_fm,p_updis,10,i,j
      !if (alpha_FM/=pperco_bsn)print*,alpha_FM,pperco_bsn,11,i,j
      !if (k_sc/=phoskd_bsn)print*,k_sc,phoskd_bsn,12,i,j
      !if (karst_he/=RCN11)print*,karst_he,RCN11,14,i,j
      !if (alpha_sc/=psp_bsn)print*,alpha_sc,psp_bsn,13,i,j
      !if (karst_he/=RCN11)print*,karst_he,RCN11,14,i,j
      !if (karst_hs/=nperco_bsn)print*,karst_hs,nperco_bsn,15,i,j
      !if (karst_hf/=eros_spl)print*,karst_hf,eros_spl,16,i,j
      !if (gw_karst_h1/=wdpq)print*,gw_karst_h1,wdpq,17,i,j
      !if (gw_karst_h2/=wgpq)print*,gw_karst_h2,wgpq,18,i,j
      !end if

      if (cswat == 2) then
          if (tillage_switch(ihru) .eq. 1) then
              if (tillage_days(ihru) .ge. 30) then
                    tillage_switch(ihru) = 0
                    tillage_days(ihru) = 0
              else
                    tillage_days(ihru) = tillage_days(ihru) + 1
              end if                
              !tillage_depth(ihru) = dtil
              !tillage_switch(ihru) = .TRUE. 
          end if
      end if
      !!by zhang DSSAT tillage  
      !!====================== 

      call varinit
      if (icr(j) <= 0) icr(j) = 1
      
      i_wtrhru = 0
      idplrot(icr(j),ihru) = idplt(j)
      if (idplt(j) /= 0) then
          if (cpnm(idplt(j)) == "WATR") then
              i_wtrhru = 1
          end if
      endif
      

      !计算干湿事件的指标
        !wdps=50.3
        !wgps=0.9
      if (api_idum==1)then
      if (station_flag==1 .and. (inum1==1 .or. inum1==2 .or.inum1==3 .or.inum1==4 .or.inum1==24 .or.inum1==25 .or.inum1==26 .or.
     &      inum1==27 .or.inum1==28 .or.inum1==29 .or.inum1==30 .or.inum1==31 .or.inum1==32.or.inum1==33 .or.inum1==34 .or.inum1==37
     & .or.inum1==43.or.inum1==44 .or.inum1==46 .or.inum1==23))then
          
          IF (curyr<=1)then
              IWAPI(j)=0.
          wapi1=0.
          wapi2=0.


          pre_pcp(j,2:int(karst_parm(20))) =pre_pcp(j,1:int(karst_parm(20))-1)
          pre_pcp(j,1)=precipday
          !print*,int(wdlpq)
          do kkk=2,int(karst_parm(20))
              wapi1=wapi1+(karst_parm(21)**kkk)*pre_pcp(j,kkk)
              wapi2=wapi2+karst_parm(21)**kkk
          enddo
      
          wapi_all(j)=wapi_all(j)+wapi1/wapi2
          if (i==365)then
              WAPI(j)=wapi_all(j)/365
              api_id=365
          endif

          !另外一种改进
          !wapi_all_1=0.
          wapi1=0.
          wapi2=0.

          pre_pcp(j,2:int(karst_parm(20))) =pre_pcp(j,1:int(karst_parm(20))-1)
          pre_pcp(j,1)=precipday
          !print*,int(wdlpq)
          do kkk=1,int(karst_parm(20))
              wapi1=wapi1+(karst_parm(21)**kkk)*pre_pcp(j,kkk)
              wapi2=wapi2+karst_parm(21)**kkk
          enddo
      
          wapi_all_1(j)=wapi_all_1(j)+wapi1/wapi2
          if (i==365)then
              WAPI_1(j)=wapi_all_1(j)/365
              api_id=365
          endif          
          !另外一种改进          
          
          
          
          else
          if (j==1)api_id=1+api_id
          wapi_all1=0.
          wapi1=0.
          wapi2=0.

          pre_pcp(j,2:int(karst_parm(20))) =pre_pcp(j,1:int(karst_parm(20))-1)
          pre_pcp(j,1)=precipday

          do kkk=2,int(karst_parm(20))
              wapi1=wapi1+(karst_parm(21)**kkk)*pre_pcp(j,kkk)
              wapi2=wapi2+karst_parm(21)**kkk
          enddo
      
          wapi_all1=wapi1/wapi2
          if (wapi_all1<=0.)wapi_all1=0.
          IWAPI(j)=(wapi_all1-WAPI(j))/(WAPI(j)+wapi_all1)     
          !IWAPI(j)=(wapi_all-WAPI(j))/(wapi_all)      
          WAPI(j)=(wapi_all1  + WAPI(j)*api_id)/(1+api_id)
          
          
       !改进形式开始
          wapi_all1=0.
          wapi1=0.
          wapi2=0.

          pre_pcp(j,2:int(karst_parm(20))) =pre_pcp(j,1:int(karst_parm(20))-1)
          pre_pcp(j,1)=precipday

          do kkk=1,int(karst_parm(20))
              wapi1=wapi1+(karst_parm(21)**kkk)*pre_pcp(j,kkk)
              wapi2=wapi2+karst_parm(21)**kkk
          enddo
      
          wapi_all1=wapi1/wapi2
          if (wapi_all1<=0)wapi_all1=0.
          IWAPI_1(j)=(wapi_all1-WAPI_1(j))/(WAPI_1(j)+wapi_all1)    
          WAPI_1(j)=(wapi_all1  + WAPI_1(j)*api_id)/(1+api_id)

          !该进形式结束          
          
          
          
          endif          
          
      else
          
          IF (curyr<=1)then
              IWAPI(j)=0.
          wapi1=0.
          wapi2=0.

          pre_pcp(j,2:int(wdps)) =pre_pcp(j,1:int(wdps)-1)
          pre_pcp(j,1)=precipday
          !print*,int(wdlpq)
          do kkk=2,int(wdps)
              wapi1=wapi1+(wgps**kkk)*pre_pcp(j,kkk)
              wapi2=wapi2+wgps**kkk
          enddo
      
          wapi_all(j)=wapi_all(j)+wapi1/wapi2
          if (i==365)then
              WAPI(j)=wapi_all(j)/365
              api_id=365
          endif
          
          
          !另外一种改进开始
          !wapi_all_1=0.
          wapi1=0.
          wapi2=0.
          pre_pcp(j,2:int(wdps)) =pre_pcp(j,1:int(wdps)-1)
          pre_pcp(j,1)=precipday
          !print*,int(wdlpq)
          do kkk=1,int(wdps)
              wapi1=wapi1+(wgps**kkk)*pre_pcp(j,kkk)
              wapi2=wapi2+wgps**kkk
          enddo
      
          wapi_all_1(j)=wapi_all_1(j)+wapi1/wapi2
          if (i==365)then
              WAPI_1(j)=wapi_all_1(j)/365
              api_id=365
          endif          
          !另外一种改进结束
          
          
          else
          
          wapi_all1=0.
          wapi1=0.
          wapi2=0.
          if (j==1)api_id=1+api_id
          pre_pcp(j,2:int(wdps)) =pre_pcp(j,1:int(wdps)-1)
          pre_pcp(j,1)=precipday

          do kkk=2,int(wdps)
              wapi1=wapi1+(wgps**kkk)*pre_pcp(j,kkk)
              wapi2=wapi2+wgps**kkk
          enddo
      
          wapi_all1=wapi1/wapi2
          if (wapi_all1<=0)wapi_all1=0.
          IWAPI(j)=(wapi_all1-WAPI(j))/(WAPI(j)+wapi_all1)    
          WAPI(j)=(wapi_all1  + WAPI(j)*api_id)/(1+api_id)
          
          
          
          
       !改进形式开始
          wapi_all1=0.
          wapi1=0.
          wapi2=0.

          pre_pcp(j,2:int(wdps)) =pre_pcp(j,1:int(wdps)-1)
          pre_pcp(j,1)=precipday

          do kkk=1,int(wdps)
              wapi1=wapi1+(wgps**kkk)*pre_pcp(j,kkk)
              wapi2=wapi2+wgps**kkk
          enddo
          !if (J==1)print*,iyr,i
          wapi_all1=wapi1/wapi2
          if (wapi_all1<=0)wapi_all=0.
          IWAPI_1(j)=(wapi_all1-WAPI_1(j))/(WAPI_1(j)+wapi_all1)    
          WAPI_1(j)=(wapi_all1  + WAPI_1(j)*api_id)/(1+api_id)
          !该进形式结束
          if (IWAPI_1(j)<0)print*,IWAPI_1(j)
          endif
      endif
      
      else
          IWAPI(j)=0.
          IWAPI_1(j)=0.
      endif
      !计算干湿事件的指标      
      !IWAPI(j)=0.
      
      if (curyr==1 .and. i==1)origin_sol_k_parm(:,j)=sol_k(:,j)
      if (curyr==1 .and. i==1)origin_gw1_parm(j)=alpha_bfe(j)
      if (curyr==1 .and. i==1)origin_gw2_parm(j)=gw_delaye(j)
      if (curyr==1 .and. i==1)origin_esco_parm(j)=esco(j)
      if (curyr==1 .and. i==1)origin_brt_parm(j)=origin_brt_parm(j)
      if (curyr==1 .and. i==1)origin_Ia=0.2
      if (curyr==1 .and. i==1)origin_ch_n(inum1)=ch_n(2,inum1)
      if (curyr==1 .and. i==1)origin_ch_k(inum1)=ch_k(2,inum1)
      !原始模型的干湿事件计算
      sol_k(:,j)=(2.71**IWAPI_1(j))*origin_sol_k_parm(:,j)
      alpha_bfe(j)=(2.71**IWAPI(j))*origin_gw1_parm(j)
      esco(j)=(2.71**IWAPI_1(j))*origin_esco_parm(j)
      !brt(j)=(2.71**IWAPI_1(j))*origin_brt_parm(j)
      esco(j)=min(esco(j),1.0)
      alpha_bfe(j)=min(alpha_bfe(j),1.0)
      !brt(j)=min(brt(j),1.0)
      if (iihru==1)then
          ch_n(2,inum1)=(2.71**(-IWAPI_1(j)))*origin_ch_n(inum1)
          ch_k(2,inum1)=(2.71**IWAPI_1(j))*origin_ch_k(inum1)
          ch_n(2,inum1)=max(0.01,ch_n(2,inum1))
          ch_n(2,inum1)=min(0.7,ch_n(2,inum1))
      endif
      
      surf_Ia=(2.71**(-IWAPI_1(j)))*origin_Ia
      surf_Ia=min(surf_Ia,0.5)
      surf_Ia=max(surf_Ia,0.00001)
      
      
      !原始模型的干湿事件计算
      if (station_flag==1 .and. (inum1==1 .or. inum1==2 .or.inum1==3 .or.inum1==4 .or.inum1==24 .or.inum1==25 .or.inum1==26 .or.
     &      inum1==27 .or.inum1==28 .or.inum1==29 .or.inum1==30 .or.inum1==31 .or.inum1==32.or.inum1==33 .or.inum1==34 .or.inum1==37
     & .or.inum1==43.or.inum1==44 .or.inum1==46 .or.inum1==23))then
      if (api_idum==1)then
      
          
         
      alpha_delay1=(2.71**(-IWAPI(j)))*karst_parm(1)!(0.1-5)下渗给表层岩溶带的延迟系数，就是比如土壤下渗量并不是一下子就补给到岩溶带中，存在滞时
      alpha_delay2=(2.71**(-IWAPI(j)))*karst_parm(2)!(0.1-5)表层岩溶带传递给基质的过程，,存在滞时
      sur_karst_h1=karst_parm(3)!(1-50)表层岩溶带发生延迟流的阈值水深
      sur_karst_h2=(karst_parm(3)+karst_parm(4))!(1-20)表层岩溶带发生二次流的阈值水深
      k_sec=(2.71**IWAPI(j))*karst_parm(5)!(0.0001-0.5)表层岩溶带二次流的出流比流量系数(线性出流)
      k_eh=(2.71**IWAPI(j))*karst_parm(6)!(0.0001-0.5)表层岩溶带延迟出流比流量系数
      alpha_eh=(2.71**IWAPI(j))*karst_parm(7)!(0.5-2)表层岩溶带延迟出流指数系数
      k_f=(2.71**IWAPI(j))*karst_parm(8)!(0.0001-0.5)表层岩溶带传递给管道流，这里与基质流传递方式不一样，线性传递
      gw_karst_delaye=karst_parm(9)!(2.71**(-IWAPI(j)))*karst_parm(9)!(0.1-5)补给地下岩溶基质的延迟系数
      gw_alpha_bfe=karst_parm(10)!(2.71**(-IWAPI(j)))*karst_parm(10)!(0.1-5)地下岩溶基质出流的延迟系数
      k_fm=(2.71**IWAPI(j))*karst_parm(11)!(0.0001-0.5)管道流的出流系数
      alpha_FM=(2.71**IWAPI(j))*karst_parm(12)!(0.5-3)管道流的比流量指数系数?
      k_sc=(2.71**IWAPI(j))*karst_parm(13)!(0.0001-0.5)基质和管道的交换系数
      alpha_sc=(2.71**IWAPI(j))*karst_parm(14)!(0.1-1)基质和管道的交换指数
      karst_he=karst_parm(15) !(1-25)表层岩溶带初始水深
      karst_hs=karst_parm(16)!(1-25)基质初始水深
      karst_hf=karst_parm(17)!(1-25)管道初始水深
      gw_karst_h1=karst_parm(18)!岩溶基质阈值水深(1-50)
      gw_karst_h2=karst_parm(19)!岩溶管道阈值水深(1-50) 
      else
      alpha_delay1=karst_parm(1)!(0.1-5)下渗给表层岩溶带的延迟系数，就是比如土壤下渗量并不是一下子就补给到岩溶带中，存在滞时
      alpha_delay2=karst_parm(2)!(0.1-5)表层岩溶带传递给基质的过程，,存在滞时
      sur_karst_h1=karst_parm(3)!(1-50)表层岩溶带发生延迟流的阈值水深
      sur_karst_h2=(karst_parm(3)+karst_parm(4))!(1-20)表层岩溶带发生二次流的阈值水深
      k_sec=karst_parm(5)!(0.0001-0.5)表层岩溶带二次流的出流比流量系数(线性出流)
      k_eh=karst_parm(6)!(0.0001-0.5)表层岩溶带延迟出流比流量系数
      alpha_eh=karst_parm(7)!(0.5-2)表层岩溶带延迟出流指数系数
      k_f=karst_parm(8)!(0.0001-0.5)表层岩溶带传递给管道流，这里与基质流传递方式不一样，线性传递
      gw_karst_delaye=karst_parm(9)!(2.71**(-IWAPI(j)))*karst_parm(9)!(0.1-5)补给地下岩溶基质的延迟系数
      gw_alpha_bfe=karst_parm(10)!(2.71**(-IWAPI(j)))*karst_parm(10)!(0.1-5)地下岩溶基质出流的延迟系数
      k_fm=karst_parm(11)!(0.0001-0.5)管道流的出流系数
      alpha_FM=karst_parm(12)!(0.5-3)管道流的比流量指数系数?
      k_sc=karst_parm(13)!(0.0001-0.5)基质和管道的交换系数
      alpha_sc=karst_parm(14)!(0.1-1)基质和管道的交换指数
      karst_he=karst_parm(15) !(1-25)表层岩溶带初始水深
      karst_hs=karst_parm(16)!(1-25)基质初始水深
      karst_hf=karst_parm(17)!(1-25)管道初始水深
      gw_karst_h1=karst_parm(18)!岩溶基质阈值水深(1-50)
      gw_karst_h2=karst_parm(19)!岩溶管道阈值水深(1-50)           
          
          
      endif
      
      if (k_sec>=0.3)k_sec=0.3
      if (k_eh>=0.3)k_eh=0.3
      if (k_f>=0.3)k_f=0.3
      if (k_fm>=0.3)k_fm=0.3
      if (k_sc>=0.3)k_sc=0.3
      if (alpha_eh>=5)alpha_eh=5.0
      if (alpha_FM>=5)alpha_FM=5.0
      if (alpha_sc>=5)alpha_sc=5.0            

          
      else
         
      if (api_idum==1)then  

      alpha_delay1=(2.71**(-IWAPI(j)))*sftmp!(0.1-5)下渗给表层岩溶带的延迟系数，就是比如土壤下渗量并不是一下子就补给到岩溶带中，存在滞时
      alpha_delay2=(2.71**(-IWAPI(j)))*smtmp!(0.1-5)表层岩溶带传递给基质的过程，,存在滞时
      sur_karst_h1=smfmn!(2.71**(-IWAPI(j)))*smfmn!(1-50)表层岩溶带发生延迟流的阈值水深
      sur_karst_h2=(smfmn+smfmx)!(2.71**(-IWAPI(j)))*(smfmn+smfmx)!(1-20)表层岩溶带发生二次流的阈值水深
      k_sec=(2.71**IWAPI(j))*timp!(0.0001-0.5)表层岩溶带二次流的出流比流量系数(线性出流)
      k_eh=(2.71**IWAPI(j))*snocovmx!(0.0001-0.5)表层岩溶带延迟出流比流量系数
      alpha_eh=(2.71**IWAPI(j))*sno50cov!(0.5-2)表层岩溶带延迟出流指数系数
      k_f=(2.71**IWAPI(j))*spcon_bsn!(0.0001-0.5)表层岩溶带传递给管道流，这里与基质流传递方式不一样，线性传递
      gw_karst_delaye=spexp_bsn!(2.71**(-IWAPI(j)))*spexp_bsn!(2.71**(-IWAPI(j)))*spexp_bsn!(0.1-5)补给地下岩溶基质的延迟系数
      gw_alpha_bfe=n_updis!(2.71**(-IWAPI(j)))*n_updis!(2.71**(-IWAPI(j)))*n_updis!(0.1-5)地下岩溶基质出流的延迟系数
      k_fm=(2.71**IWAPI(j))*p_updis!(0.0001-0.5)管道流的出流系数
      alpha_FM=(2.71**IWAPI(j))*pperco_bsn!(0.5-3)管道流的比流量指数系数?
      k_sc=(2.71**IWAPI(j))*phoskd_bsn!(0.0001-0.5)基质和管道的交换系数
      alpha_sc=(2.71**IWAPI(j))*psp_bsn!(0.1-1)基质和管道的交换指数
      karst_he=RCN11 !(1-25)表层岩溶带初始水深
      karst_hs=nperco_bsn!(1-25)基质初始水深
      karst_hf=eros_spl!(1-25)管道初始水深
      gw_karst_h1=wdpq!(2.71**(-IWAPI(j)))*wdpq!岩溶基质阈值水深(1-50)
      gw_karst_h2=wgpq!(2.71**(-IWAPI(j)))*wgpq!岩溶管道阈值水深(1-50)           
          

      
      else
      alpha_delay1=sftmp!(0.1-5)下渗给表层岩溶带的延迟系数，就是比如土壤下渗量并不是一下子就补给到岩溶带中，存在滞时
      alpha_delay2=smtmp!(0.1-5)表层岩溶带传递给基质的过程，,存在滞时
      sur_karst_h1=smfmn!(1-50)表层岩溶带发生延迟流的阈值水深
      sur_karst_h2=(smfmn+smfmx)!(1-20)表层岩溶带发生二次流的阈值水深
      k_sec=timp!(0.0001-0.5)表层岩溶带二次流的出流比流量系数(线性出流)
      k_eh=snocovmx!(0.0001-0.5)表层岩溶带延迟出流比流量系数
      alpha_eh=sno50cov!(0.5-2)表层岩溶带延迟出流指数系数
      k_f=spcon_bsn!(0.0001-0.5)表层岩溶带传递给管道流，这里与基质流传递方式不一样，线性传递
      gw_karst_delaye=spexp_bsn!(2.71**(-IWAPI(j)))*spexp_bsn!(0.1-5)补给地下岩溶基质的延迟系数
      gw_alpha_bfe=n_updis!(2.71**(-IWAPI(j)))*n_updis!(0.1-5)地下岩溶基质出流的延迟系数
      k_fm=p_updis!(0.0001-0.5)管道流的出流系数
      alpha_FM=pperco_bsn!(0.5-3)管道流的比流量指数系数?
      k_sc=phoskd_bsn!(0.0001-0.5)基质和管道的交换系数
      alpha_sc=psp_bsn!(0.1-1)基质和管道的交换指数
      karst_he=RCN11 !(1-25)表层岩溶带初始水深
      karst_hs=nperco_bsn!(1-25)基质初始水深
      karst_hf=eros_spl!(1-25)管道初始水深
      gw_karst_h1=wdpq!岩溶基质阈值水深(1-50)
      gw_karst_h2=wgpq!岩溶管道阈值水深(1-50)             
          
          
      endif
      
      
      if (k_sec>=0.3)k_sec=0.3
      if (k_eh>=0.3)k_eh=0.3
      if (k_f>=0.3)k_f=0.3
      if (k_fm>=0.3)k_fm=0.3
      if (k_sc>=0.3)k_sc=0.3
      if (alpha_eh>=5)alpha_eh=5.0
      if (alpha_FM>=5)alpha_FM=5.0
      if (alpha_sc>=5)alpha_sc=5.0      
      end if
      

      
	if (i_wtrhru == 1) then
         call water_hru
      else 

        !! Simulate land covers other than water

        !! update base zero total heat units
        if (tmpav(j) > 0. .and. phutot(hru_sub(j)) > 0.01) then
           phubase(j) = phubase(j) + tmpav(j) / phutot(hru_sub(j))
        end if
        
        call schedule_ops

        !! calculate albedo for day
        call albedo

        !! calculate soil temperature for soil layers
        call solt
        !if (j==1 )print*,i,kkk,sol_sw(1)
!       if (ipot(j) /= j .and. imp_trig(nro(j),nrelease(j),j)==1)       &  Srini pothole
!
!     &        then             
          !! calculate surface runoff if HRU is not impounded or an 
          !! undrained depression--
          call surface

          !! add surface flow that was routed across the landscape on the previous day
       !!   qday = qday + surfq_ru(j)
       !!   surfq_ru(j) = 0.
          
          !! 计算有效降雨量(percs进入土壤的量)

          inflpcp = Max(0.,precipday - surfq(j))!+inflpcp)
          !if (j==378)print*,precipday,inflpcp,surfq(j),-2

!        end if
         
        !! 执行管理操作
        if (yr_skip(j) == 0) call operatn
          
        if (auto_wstr(j) > 1.e-6 .and. irrsc(j) > 2) call autoirr    
        if (karst_hru(j)==2)inflpcp=precipday+inflpcp+surf_bs(1,j) !岩溶洼地
        if (hru_slp(j)*100>50 .and. karst_hru(j)==1)go to 888!高坡度，石漠化区域，忽略土壤层的计算，直接到表层岩溶带
        !! 执行土壤水分路线
         
        call percmain
888     continue
        !! 计算蒸发蒸腾量
        call etpot
        if (hru_slp(j)*100>50 .and.karst_hru(j)==1 )go to 889!高坡度，石漠化区域，忽略土壤层的计算，直接到表层岩溶带        
!        if (pot_vol(j) < 1.e-6) call etact
        call etact!土壤蒸发
889     continue
        
        !岩溶洼地的计算
        !if (karst_hru(j)==2)call karst_pond!岩溶洼地hru
        
        !表层岩溶带的计算
        if (karst_hru(j)>=1)call surface_karst!岩溶hru
        !if (karst_hru(j)>=1 .and. j==546)print*,i,sur_karst_height(j),Q_sur_sec(j),Q_sur_eh(j),Q_karst_S(j),Q_karst_F(j)
        !! 使用气候驱动因素计算地下水位深度
        call wattable

        !! new CN method
        if (icn == 1) then 
        sci(j) = sci(j) + pet_day*exp(-cncoef_sub(hru_sub(j))*sci(j)/   
     &    smx(j)) - precipday + qday + qtile + latq(j) + sepbtm(j)
        else if (icn == 2) then 
        sci(j) = sci(j) + pet_day*exp(-cncoef_sub(hru_sub(j))*sci(j)/   
     &    smx(j)) - precipday + qday + latq(j) + sepbtm(j) + qtile
        sci(j) = amin1(sci(j),smxco * smx(j))
        end if 
        
        !! 在连续施肥作业中施用肥料
        if (icfrt(j) == 1) then
          ndcfrt(j) = ndcfrt(j) + 1
          call confert
        end if
        
        !! 农药在害虫连续作业中的应用
        if (icpst(j) == 1) then 
          ndcpst(j) = ndcpst(j) + 1
          call conapply
        end if 
        
        !! 清除放牧生物量并施用肥料
        if (igrz(j) == 1) then
          ndeat(j) = ndeat(j) + 1
          call graze
        end if
       
        !! 计算作物生长
        call plantmod
        
        !! 检查休眠
        if (igro(j) == 1) call dormant
        !! compute actual ET for day in HRU计算HRU中当天的实际ET
        etday = ep_day + es_day + canev!在etact计算

!! 编写每日气温和土壤温度文件
!! 如果用户需要，也可以在readfile.f中取消注释

!      write (120,12112) i,j,tmx(j),tmn(j),(sol_tmp(k,j),k=1,sol_nly(j))
!12112  format (2i4,12f8.2)

        !! compute nitrogen and phosphorus mineralization 

      if (cswat == 0) then
        call nminrl
	end if
	if (cswat == 1) then
		call carbon
	end if
	
	!! Add by zhang
	!!=================
	if (cswat == 2) then
	  call carbon_zhang2
	end if
	!! Add by zhang
	!!=================	

        call nitvol
        if (sol_P_model == 1) then
            call pminrl
        else
            call pminrl2
        end if
!!    compute biozone processes in septic HRUs
!!    if 1)current is septic hru and 2)  soil temperature is above zero
!!    计算脓毒症重症监护室中的生物区过程
!!    如果1）电流为化粪池，2）土壤温度高于零度
	  if (isep_opt(j)/=0.and.iyr>=isep_iyr(j)) then
	   if (sol_tmp(i_sep(j),j) > 0.) call biozone     
	  endif

        !! 计算地下水贡献
        call gwmod
        call gwmod_deep
        !if (gw_karst_hru(inum1)==1 .and. karst_hru(j)>=1)call gw_karst!地下河  ,这里加karst_hru(j)>=1的目的是为了保证该子流域为岩溶区域

        !! 计算农药冲刷  
        if (precipday >= 2.54) call washp

        !! 计算农药降解
        call decay

        !!计算农药在土壤中的移动
        call pestlch

        if (surfq(j) > 0. .and. peakr > 1.e-6) then
          if (precipday > 0.) then
            call enrsb(0)
            if (sedyld(j) > 0.) call pesty(0)

		  if (cswat == 0) then
			call orgn(0)
	    end if
	    if (cswat == 1) then
	    
		    call orgncswat(0)
		  end if
		  
		  !! Add by zhang
		  !! ====================
		  if (cswat == 2) then
		    call orgncswat2(0)
		  end if
		  !! Add by zhang
		  !! ====================

            call psed(0)
          end if
        end if

        !! 在降雨中向土壤剖面添加硝酸盐
        call nrain

        !! 计算硝酸盐移动浸出
        call nlch

        !! 计算磷的运动
        call solp

        !! 计算细菌运输
        call bacteria

        !! 计算城市地区的负荷
        if (urblu(j) > 0) then
	     if(ievent == 0) then
	        call urban ! daily simulation
	     else
		     call urbanhr ! subdaily simulation J.Jeong 4/20/2009
	     endif
	  endif
	  
!! Srini Pothole
        !! 计算不排水洼地/蓄水面积（如水稻）过程
        !if (pot_fr(j) > 0.) then
        !   if (ievent == 0) then   
        !  call pothole
        !   else
        !      call potholehr
        !   endif
        !endif
        
        !! 计算横向流中的泥沙荷载并添加到sedyld中
        call latsed

        !! 计算地下水流中的养分负荷
        call gwnutr
        call gw_no3

        !! 地表径流中的滞后养分和沉积物
        call surfstor

        !! 滞后地下水流和地下水流中的硝酸盐

        call substor

        !! 添加前一天穿过景观的横向流量
      !!  latq(j) = latq(j) + latq_ru(j)
      !!  latq_ru(j) = 0.
        
        !! 计算由于现场过滤带边缘造成的污染物减少
        if (vfsi(j) >0.)then
          call filter
          if (filterw(j) > 0.) call buffer
        end if
              if (vfsi(j) == 0. .and. filterw(j) > 0.) then 
                call filtw
                call buffer
              end if

	 !! 计算田间草地水道造成的污染物减少量
         if (grwat_i(j) == 1) then
          call grass_wway
        end if

	 !! 计算固定BMP eff导致的污染物减少
	 !  if (bmp_flag(j) == 1) then
       !   call bmpfixed
       ! end if
          if (sub_karst==1)then!计算刁江
              if (karst_hru(j)==1)then!为岩溶hru
              !这些子流域在原始流域中是不存在的，因为根据dem水文分析没有包含这些，如果是开始进行岩溶流域计算，只需把以下子流域的地下水加上就行
                  if (inum1==13.or.inum1==14 .or.inum1==15 .or.inum1==16.or.inum1==17.or.inum1==18.or.inum1==20 .or. inum1==21)then 
                     qdr(j) =  latq(j) + gw_q(j)+ gw_qdeep(j)+Q_karst_FM(j)+ Q_sur_sec(j)+ Q_sur_eh(j)
      !               if (iyr==2008 .and.(i==164 .or. i==165))print*,1
      !               if (iyr==2008 .and.(i==164 .or. i==165))print*,precipday,latq(j) , gw_q(j)+ gw_qdeep(j)+Q_karst_FM(j)+
      !&           Q_sur_sec(j)+Q_sur_eh(j)
                  else
                     qdr(j) = qday + latq(j) + gw_q(j) + qtile + gw_qdeep(j) +Q_karst_FM(j)+ Q_sur_sec(j)+ Q_sur_eh(j)
      !               if (iyr==2008 .and.(i==164 .or. i==165 .or. i==163).and. (J==1))print*,i,J
      !               if (iyr==2008 .and.(i==164 .or. i==165 .or. i==163).and. (J==1))print*,precipday,qday,
      !&               latq(j) ,  gw_q(j) + qtile + gw_qdeep(j) +Q_karst_FM(j)+ Q_sur_sec(j)+ Q_sur_eh(j)
      !               if (iyr==2008 .and.(i==164 .or. i==165 .or. i==163).and. (J==1))print*,inflpcp,sepbtm(j)     
      !               if (iyr==2008 .and.(i==164 .or. i==165 .or. i==163).and. (J==1))print*,  ""                    
                  endif
                  
              elseif (karst_hru(j)==2)then!岩溶洼地
                     qdr(j) =  latq(j) + gw_q(j)+ gw_qdeep(j)+Q_karst_FM(j)+ Q_sur_sec(j)+ Q_sur_eh(j)
      !               if (iyr==2008 .and.(i==164 .or. i==165))print*,3
      !                if (iyr==2008 .and.(i==164 .or. i==165))print*,precipday,qday , latq(j) ,gw_q(j)+ gw_qdeep(j)+Q_karst_FM(j)+
      !&                 Q_sur_sec(j)+ Q_sur_eh(j)
              else!不是岩溶hru
                 qdr(j) = qday + latq(j) + gw_q(j) + qtile + gw_qdeep(j) 
                 !if (iyr==2008 .and.(i==164 .or. i==165))print*,4
                 !if (iyr==2008 .and.(i==164 .or. i==165))print*,precipday,qday , latq(j) , gw_q(j), gw_qdeep(j) 
                 !
              endif
              

          else!计算滂江
              if (karst_hru(j)==1)then!为岩溶hru
                  if (inum1==48.or.inum1==49 )then!子流域在隔壁流域的
                     qdr(j) = latq(j) + gw_q(j)+ gw_qdeep(j)+Q_karst_FM(j)+ Q_sur_sec(j)+ Q_sur_eh(j)
                  else 
                     qdr(j) = qday + latq(j) + gw_q(j) + qtile + gw_qdeep(j) +Q_karst_FM(j)+ Q_sur_sec(j)+ Q_sur_eh(j)
                  endif
              elseif (karst_hru(j)==2)then!为岩溶洼地hru
                   qdr(j) =  latq(j) + gw_q(j)+ gw_qdeep(j)+Q_karst_FM(j)+ Q_sur_sec(j)+ Q_sur_eh(j)
              else!不是岩溶hru
                  qdr(j) = qday + latq(j) + gw_q(j) + qtile + gw_qdeep(j) 
              endif
              
          end if
              !if (i==155 .and. inum1==27)print*,qdr(j),karst_hru(j),Q_karst_FM(j),Q_sur_sec(j), Q_sur_eh(j)
        !! compute water yield for HRU
        
        if (qdr(j) < 0.) qdr(j) = 0.
        if (qdr(j) > 0.) then
          qdfr = qday / qdr(j)
        else
          qdfr = 0.
        end if

        !! 计算chl-a、CBOD和溶解氧负荷
        call subwq

        !! 计算湿地过程
        call wetlan

        !! 计算池过程
        if (ievent == 0) then
           call hrupond
        else
           call hrupondhr
        endif
        
!       斯里尼坑洞       
        if (pot_fr(j) > 0.) call pothole
                
        xx = sed_con(j)+soln_con(j)+solp_con(j)+orgn_con(j)+orgp_con(j)
        if (xx > 1.e-6) then
          call urb_bmp
        end if

        !! 消耗性用水（池塘、浅层含水层、深层含水层）
        call watuse

        !! 进行水平衡
        call watbal
        
        !! qdayout是在湿地、池塘和坑洞之后离开hru的地表径流
        qdayout(j) = qday

      endif

      !! 执行输出摘要
      call sumv

!! 汇总每个子流域多个人力资源单元的输出
!! 一种新的fig方法的存储到达载荷
777            continue         
      call virtual
      aird(j) = 0.
      
      ihru = ihru + 1
      
      end do
		
      !!路线2景观单元
      if (ils2flag(inum1) > 0) then
      isub = inum1                        ! save the subbasin number
      
      !! 从山坡计算输出
      ihout1 = mhyd_bsn + (inum1 - 1) * 4 ! first outflow hyd number
      ihout = ihout1                      ! outflow hyd number
      inum1 = 1                           ! landscape unit number
      inum2 = isub                        ! subbasin number
      call routeunit                      ! hillslope unit
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout
      
      !! 从谷底计算输出
      inum1 = 2                           ! landscape unit number
      ihout = ihout + 1                   ! outflow hyd number
      sumdaru = 0.
      do j = 1, hrutot(isub)
        sumdaru = sumdaru + hru_km(j)
      end do 
      daru_km(inum2,inum1) = sumdaru
      call routeunit                      ! valley bottom unit
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout
      
      !! 从山坡到谷底的路线输出
      ihout = ihout + 1                   ! outflow hyd number
      inum1 = 2                           ! valley bottom landscape unit
      inum2 = ihout1                      ! inflow hyd=outlfow from hillslope
      inum3 = isub                        ! subbasin number
      rnum1 = 1.                          ! fraction overland flow
      iru_sub = 1                         ! route across landscape unit
      !! 计算输沙能力的加权K系数
      sumk = 0.
      ovsl = 0.
      ovs = 0.
      do j = 1, hrutot(isub)
        sumk = sumk + usle_k(j) * hru_rufr(inum1,j)
        ovsl = ovsl + slsubbsn(j)
        ovs = ovs + hru_slp(j)
      end do 
      ovsl = ovsl / hrutot(isub)
      ovs = ovs / hrutot(isub)
      ru_k(isub,inum1) = sumk
      ru_ovsl(isub,inum1) = ovsl
      ru_ovs(isub,inum1) = ovs
      ru_ktc(isub,inum1) = 50.
      ru_a(isub,inum1) = daru_km(isub,1) / ru_ovsl(isub,inum1)
      call routels(iru_sub)               ! route across valley bottom
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      inum3s(ihout) = inum3
      ihouts(ihout) = ihout
      
      !! 添加带有谷底负载的路线
      inum1 = ihout                       ! hyd from routed 
      inum2 = ihout - 1                   ! hyd from loading
      ihout = ihout + 1                   ! outflow hyd number
      call addh                           ! add hyd's
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout
      
      !! save landscape routed output in place of subbasin output for routing
      varoute(isub,:) = varoute(ihout,:)
      end if
      
 1000 format(4i10,a10)
      return
      end