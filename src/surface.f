      subroutine surface

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine models surface hydrology at any desired time step

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ihru        |none          |HRU number
!!    ovrlnd(:)   |mm H2O        |overland flow onto HRU from upstream
!!                               |routing unit
!!    peakr       |mm/hr         |peak runoff rate
!!    precipday   |mm H2O        |effective precipitation for the day in HRU
!!    qday        |mm H2O        |surface runoff loading to main channel 
!!                               |for day
!!    surfq(:)    |mm H2O        |surface runoff generated in HRU during
!!                               |the day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    precipday   |mm H2O        |effective precipitation for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Max
!!    SWAT: canopyint, snom, crackvol, dailycn, volq, crackflow, surfst_h2o,
!!    SWAT: alph, pkq, tran, eiusle, ysed

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


      use parm

      integer :: j,sb,kk
      real :: precip_fr
      real :: irfr,hruvirr

      j = 0
      j = ihru
      sb = hru_sub(j)
      hruirrday = 0.
      irmmdt = 0.

      !! compute canopy interception
      if (idplt(j) > 0) then
        call canopyint
      end if

      !! compute snow melt
      call snom
      !if (j==1)print*,i,kkk,precipdt(:)
      !! output by elevation band to output.snw
      if (isnow == 1) then    
        write(115,1010) i, iyr, subnum(j), hruno(j),                   
     &                (snoeb(ib,j), ib = 1,10)
      end if
      
      if (isnow ==1) then
        write (116,1010) i, iyr, subnum(j), hruno(j),
     &    (tavband(ib,j), ib = 1, 10)
      end if
      
      !! compute crack volume
      if (icrk == 1) call crackvol

      !! add overland flow from upstream routing unit
      precipday = precipday + ovrlnd(j)
      if (nstep > 0) then
        do ii = 1, nstep
          precipdt(ii+1) = precipdt(ii+1) + ovrlnd_dt(j,ii)
        end do
      end if
      !if (j==1)print*,i,kkk,precipdt(:),tmx(j)
      !! add irrigation from retention-irrigation ponds to soil water
      if (ri_luflg(j)==1) then
        irfr = hru_km(j)* (1.-fimp(urblu(j))) / ri_subkm(sb) 
        do ii=1,nstep
          !amount irrigated in hru
          hruvirr = ri_totpvol(ii) * irfr !m3
          irmmdt(ii) = hruvirr / (hru_km(j) 
     &       * (1.- fimp(urblu(j))) * 1000.) !mm/dt
          
          !add irrigated water to soil water content
          do kk=1,sol_nly(j)
            if(irmmdt(ii)<sol_ul(kk,j)-sol_st(kk,j)) then
               sol_st(kk,j) = sol_st(kk,j) + irmmdt(ii)
               exit
            else
               sol_st(kk,j) = sol_ul(kk,j)
               irmmdt(ii) = irmmdt(ii) - (sol_ul(kk,j)-sol_st(kk,j))
            end if
          end do
         
        end do
      end if

      !!calculate subdaily curve number value
      call dailycn

        !! compute runoff - surfq in mm H2O
      if (precipday > 0.1) then
        call volq 
        !! bmp adjustment
        surfq(j) = surfq(j) * bmp_flo(j)
        !! adjust runoff for loss into crack volume
        if (surfq(j) > 0. .and. icrk == 1) call crackflow
      end if

      surfq(j) = surfq(j) + qird(j)
      qird(j) = 0.

      !! !! 计算白天（qday）到达主河道的地表径流量，并存储剩余量
      call surfst_h2o

      !! calculate half-hour rainfall计算半小时降雨量
      if (precipday > 0.01) call alph(0)

      if (qday > 0.0001) then
        !! compute peak rate - peakr in m3/s  计算洪峰速度
        call pkq(0)  
      end if  

      if (qday > 0.0001 .and. peakr > 0.) then
        !! 计算非HUMUS数据集的传输损耗
        call tran
        call eiusle

	!! 计算降雨和地表径流对泥沙的侵蚀
		call ovr_sed
      end if

      call cfactor
      if (surfq(j) > 1.e-6 .and. peakr > 1.e-6) call ysed(0)

      if (qday < 0.) qday = 0.
      
1010  format (2(i4,1x),a5,a4,1x,10f8.1)
      return
      end