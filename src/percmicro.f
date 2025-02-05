      subroutine percmicro(ly1)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine computes percolation and lateral subsurface flow
!!    from a soil layer when field capacity is exceeded

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hru_slp(:)   |m/m           |average slope steepness
!!    ihru         |none          |HRU number
!!    iwatable     |none          |high water table code:
!!                                |0 no high water table
!!                                |1 high water table
!!    ldrain(:)    |none          |soil layer where drainage tile is located
!!    slsoil(:)    |m             |slope length for lateral subsurface flow侧向地下水流的斜坡长度
!!    sol_fc(:,:)  |mm H2O        |amount of water available to plants in soil
!!                                |layer at field capacity (fc - wp water)
!!    sol_hk(:,:)  |none          |beta coefficient to calculate hydraulic 用于计算水力传导率的β系数
!!                                |conductivity
!!    sol_k(:,:)   |mm/hr         |saturated hydraulic conductivity of soil
!!                                |layer
!!    sol_nly(:)   |none          |number of soil layers in HRU
!!    sol_st(:,:)  |mm H2O        |amount of water stored in the soil layer
!!                                |on any given day
!!    sol_sumfc(:) |mm H2O        |amount of water held in the soil profile
!!                                |at field capacity
!!    sol_sw(:)    |mm H2O        |amount of water stored in the soil profile
!!                                |on any given day
!!    sol_tmp(:,:) |deg C         |daily average temperature of soil layer
!!    sol_ul(:,:)  |mm H2O        |amount of water held in the soil layer at
!!                                |saturation (sat - wp water)饱和时土层中的含水量（sat-wp水）
!!    sol_z(:,:)   |mm            |depth to bottom of soil layer至土层底部的深度
!!    sw_excess    |mm H2O        |amount of water in soil that exceeds field 
!!                                |capacity (gravity drained water)
!!    tdrain(:)    |hrs           |time to drain soil to field capacity
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    latlyr       |mm H2O        |lateral subsurface flow in layer
!!    lyrtile      |mm H2O        |drainage tile flow in layer for day in HRUHRU中一天的排水砖分层流动
!!    sepday       |mm H2O        |percolation from soil layer土层渗流
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    adjf         |none          |adjustment factor for lateral flow横向流量调整系数
!!    dg           |mm            |depth of soil layer
!!    ho           |none          |variable to hold intermediate calculation
!!                                |result
!!    j            |none          |HRU number
!!    ly1          |none          |soil layer number
!!    ratio        |none          |ratio of seepage to (latq + sepday)
!!    yy           |mm            |depth to top of soil layer至土层顶部的深度
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer, intent (in) :: ly1
      integer :: j
      real :: adjf, yy, dg, ho, ratio, sol_k_sep

      j = 0
      j = ihru

      adjf = 1.

!! 如果层的温度为0摄氏度或更低，则没有水流
      if (sol_tmp(ly1,j) <= 0.) then
        sepday = 0.
        return
      end if

        !! 利用山坡蓄水法计算横向流量
        if (ly1 == 1) then
          yy = 0.
        else
          yy = 0.
          yy = sol_z(ly1-1,j)
        end if

        dg = 0.
        ho = 0.
        latlyr = 0.
        dg = sol_z(ly1,j) - yy
        if (sol_ul(ly1,j) - sol_fc(ly1,j)==0.) then
          ho=0.
        else
          ho = 2. * sw_excess / ((sol_ul(ly1,j) - sol_fc(ly1,j)) /  dg)
        end if
        latlyr = adjf * ho * sol_k(ly1,j) * hru_slp(j) / slsoil(j)      
     &                                                            * .024!式子2-127

      if (latlyr < 0.) latlyr = 0. 
      if (latlyr > sw_excess) latlyr = sw_excess

      sol_hk(ly1,j) = (sol_ul(ly1,j) - sol_fc(ly1,j)) / sol_k(ly1,j)!式子2-107

!!  septic changes 1/28/09 
      if (ly1 == i_sep(j)) then
         if (isep_opt(j) == 1) then !active system
           sol_k_sep = sol_k(ly1,j)* 
     &              (sol_st(ly1,j) - sol_fc(ly1,j))/
     &              (sol_ul(ly1,j) - sol_fc(ly1,j))
           sol_k_sep = Max(1.e-6, sol_k_sep)
           sol_k_sep = Min(sol_k(ly1,j), sol_k_sep)
           
           sol_hk(ly1,j) = (sol_ul(ly1,j) - sol_fc(ly1,j)) 
     &      / sol_k_sep
         
         elseif (isep_opt(j) == 2) then !failing system
           sol_hk(ly1,j) = 1.e10
         endif
      endif 
!!  septic changes 1/28/09
      
      sol_hk(ly1,j) = Max(2., sol_hk(ly1,j))

      !! compute seepage to the next layer
      sepday = 0.
      sepday = sw_excess * (1. - Exp(-24. / sol_hk(ly1,j)))!式子2-106
      !if (karst_hru(j)>=1)sepday = sw_excess * (1. - Exp(-24. / (wdlpq*sol_hk(ly1,j))))
      !! 将生物带层的最大渗流限制在潜在perc量以下
	if(ly1 == i_sep(j).and.isep_opt(j)==1) then
	   sepday = min(sepday,sol_k_sep *24.)
	   bz_perc(j) = sepday
	end if
      
      !! 如果下一层饱和，则限制渗流
      if (ly1 == sol_nly(j)) then
        xx = (dep_imp(j) - sol_z(ly1,j)) / 1000.
        
        if (xx < 1.e-4) then
          sepday = 0.
        else
          sepday = sepday * xx / (xx + Exp(8.833 - 2.598 * xx))!式子2-117
        end if
      end if

      !! 检查质量平衡
      if (sepday + latlyr > sw_excess) then
        ratio = 0.
        ratio = sepday / (latlyr + sepday)
        sepday = 0.
        latlyr = 0.
        sepday = sw_excess * ratio
        latlyr = sw_excess * (1. - ratio)
      endif
      if (sepday + lyrtile > sw_excess) then
        sepday = 0.
        sepday = sw_excess - lyrtile
      endif
      !if (j==378)print*,sepday,ly1 ,i,-1
      return
      end