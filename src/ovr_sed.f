      subroutine ovr_sed()
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine computes splash erosion by raindrop impact and flow erosion by overland flow

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cht(:)      |m             |canopy height
!!    fimp(:)     |fraction      |fraction of HRU area that is
!!                               |impervious (both directly and
!!                               |indirectly connected)
!!    hhqday(:)   |mm H2O        |surface runoff generated each timestep 
!!                               |of day in HRU
!!    hru_km(:)   |km2           |area of HRU in square kilometers
!!    idt         |minutes       |length of time step used to report
!!    inum1       |none          |subbasin number
!!    laiday(:)   |m2/m2         |leaf area index叶面积指数
!!    rainsub(:,:)|mm H2O        |precipitation for the time step during the
!!                               |day in HRU
!!    eros_spl	  |none          |coefficient of splash erosion varing 0.9-3.1
!!    urblu(:)    |none          |urban land type identification number from
!!                               |urban.dat
!!    usle_k(:)   |              |USLE equation soil erodibility (K) factor
!!
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhsedy(:,:)|tons           |sediment yield from HRU drung a time step
!!                               |applied to HRU

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!	  bed_shear		|N/m2		   |shear stress b/w stream bed and flow	
!!	  erod_k		|g/J		   |soil detachability value	土壤可剥离性值
!!    jj			|none          |HRU number
!!    kk			|none          |time step of the day
!!	  ke_direct		|J/m2/mm	   |rainfall kinetic energy of direct throughfall
!!	  ke_leaf		|J/m2/mm	   |rainfall kinetic energy of leaf drainage
!!	  ke_total		|J/m2   	   |total kinetic energy of rainfall
!!	  percent_clay	|percent	   |percent clay
!!	  percent_sand	|percent	   |percent sand
!!	  percent_silt	|percent	   |percent silt
!!	  pheff     	|m			   |effective plant height
!!	  rdepth_direct	|mm			   |rainfall depth of direct throughfall|直接穿透降雨深度
!!	  rdepth_leaf	|mm			   |rainfall depth of leaf drainage叶片排水的降雨深度
!!	  rdepth_tot	|mm			   |total rainfall depth 总降雨深度
!!    rintnsty	    |mm/hr         |rainfall intensity降雨强度
!!	  sedspl		|tons		   |时间步长内降雨影响的产沙量
!!	  sedov 		|tons		   |时间步长内陆流产沙量
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: log10, Exp, Real
!!
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!!  Splash erosion model is adopted from EUROSEM model developed by Morgan (2001).
!!	Rill/interill erosion model is adoped from Modified ANSWERS model by Park et al.(1982)
!!  Code developed by J. Jeong and N. Kannan, BRC.

      use parm

	integer :: k, j
	real :: percent_clay, percent_silt, percent_sand, erod_k
	real :: ke_direct, ke_leaf, ke_total,pheff, c
	real :: rdepth_direct, rdepth_leaf, rdepth_tot, canopy_cover
	real :: bed_shear, sedov, sedspl, rain_d50, rintnsty

	j = ihru
	
!! Fraction of sand
      percent_clay = sol_clay(1,j)
	percent_silt = sol_silt(1,j) 
	percent_sand = 100. - percent_clay - percent_silt

!! 土壤可剥离性值采用EUROSEM用户指南（表1）
	if ((percent_clay>=40.) .and. (percent_sand>=20.) .and. 
     &(percent_sand<=45.)) then
	  erod_k = 2.0 !clay
      elseif ((percent_clay>=27.) .and. (percent_sand>=20.) .and.
     &(percent_sand<=45.)) then
	  erod_k = 1.7 !Clay loam
      elseif ((percent_silt<=40.).and.(percent_sand<=20.)) then
	  erod_k = 2.0 !Clay
      elseif ((percent_silt>40.).and.(percent_clay>=40.)) then
	  erod_k = 1.6 !Silty clay
      elseif ((percent_clay>=35.).and.(percent_sand>=45.)) then
	  erod_k = 1.9 !Sandy clay
      elseif ((percent_clay>=27.).and.(percent_sand<20.)) then
	  erod_k = 1.6 !Silty clay loam
      elseif ((percent_clay<=10.).and.(percent_silt>=80.)) then
	  erod_k = 1.2 !Silt
      elseif (percent_silt>=50.) then
	  erod_k = 1.5 !Silt loam
      elseif ((percent_clay>=7.) .and. (percent_sand<=52.) .and. 
     &(percent_silt>=28.)) then
	  erod_k = 2.0 !Loam
      elseif (percent_clay>=20.) then
	  erod_k = 2.1 !Sandy clay loam
      elseif (percent_clay>=percent_sand-70.) then
	  erod_k = 2.6 !Sandy loam
      elseif (percent_clay>=(2. * percent_sand) - 170.) then
	  erod_k = 3.0 !Loamy sand
      else
	  erod_k = 1.9 !Sand 
      end if
	
 !!	基于叶面积指数的冠层覆盖
!!	如果LAI>=1，则假设冠层覆盖率为100%
	  if(laiday(j)>=1.) then
	    canopy_cover = 1.
	  else
	    canopy_cover = laiday(j)
	  end if

      do k=1,nstep
	  rintnsty = 60. * rainsub(j,k) / Real(idt) 
	  rain_d50 = 0.188 * rintnsty ** 0.182

	  if(rintnsty>0) then
	
      !!直接穿透产生的降雨动能（J/m^2/mm）
	    ke_direct = 8.95 + 8.44 * log10(rintnsty)
        if(ke_direct<0.) ke_direct = 0.
	  !! 叶片排水产生的降雨动能（J/m^2）
	    pheff = 0.5 * cht(j)
	    ke_leaf = 15.8 * pheff ** 0.5 - 5.87
	    if (ke_leaf<0) ke_leaf = 0.

	  !! 雨量
	    rdepth_tot = rainsub(j,k) / (idt * 60.)
	    rdepth_leaf = rdepth_tot * canopy_cover 
	    rdepth_direct = rdepth_tot - rdepth_leaf 
	  else
	    ke_direct = 0.
	    ke_leaf = 0.
		rdepth_tot = 0.
	    rdepth_leaf = 0.
	    rdepth_direct = 0.
	  endif

	!! 降雨总动能（J/m^2）
	  ke_total = 0.001 * (rdepth_direct * ke_direct + rdepth_leaf * 
     &	  ke_leaf)

	!!雨滴撞击引起的土壤总剥离
	  sedspl = erod_k * ke_total * exp(-1 * hhqday(k) / 1000.) * !时间步长内降雨影响的产沙量
     &	  hru_km(j) ! tons

	!! Impervious area of HRU
	  if(urblu(j)>0) sedspl = sedspl * (1.- fimp(urblu(j)))

	!! 允许飞溅侵蚀的最大水深
	  if(hhqday(k)>=3.* rain_d50.or.hhqday(k)<=1.e-3) sedspl = 0.


	!! 地表径流侵蚀
    !! usle方程（ysed.f）中使用的覆盖和管理因子
        if (idplt(j) > 0) then     
          c = Exp((-.2231 - cvm(idplt(j))) *                            
     &	      Exp(-.00115 * sol_cov(j)) + cvm(idplt(j)))         !式子4-10     
        else
          if (sol_cov(j) > 1.e-4) then
            c = Exp(-.2231 * Exp(-.00115 * sol_cov(j)))               
          else
            c = .8
          end if
	  end if
	!! 5厘格栅处的水的比重=9807N/m3
	  bed_shear = 9807 * (hhqday(k) / 1000.) * hru_slp(j) ! N/m2
	  sedov = 11.02 * rill_mult * usle_k(j) * c_factor * c * 
     &	  bed_shear ** eros_expo ! kg/hour/m2
	  if(ievent > 0) then
	    sedov = 16.667 * sedov * hru_km(j) * idt ! tons per time step
	  else
	    sedov = 24000. * sedov * hru_km(j)	! tons per day
	  endif

	!! Impervious area of HRU
	  if(urblu(j)>0) sedov = sedov * (1.- fimp(urblu(j)))

	  hhsedy(j,k) = dratio(inum1) * (sedspl + sedov)
	  if (hhsedy(j,k) < 1.e-10) hhsedy(j,k) = 0.

	end do
	return
      end