       subroutine zero_karst


      use parm
      
      alpha_delay1=0
      sftmp=1!下渗给表层岩溶带的延迟系数，就是比如土壤下渗量并不是一下子就补给到岩溶带中，存在滞时
      alpha_delay2=0
      smtmp=1!表层岩溶带传递给基质的过程，存在滞时
      sur_karst_h1=0
      smfmn=1!表层岩溶带发生延迟流的阈值水深
      sur_karst_h2=0
      smfmx=2!表层岩溶带发生二次流的阈值水深
      k_sec=0
      timp=0.5!表层岩溶带二次流的出流比流量系数(线性出流)
      k_eh=0.
      SNOCOVMX=0.5!表层岩溶带延迟出流比流量系数
      alpha_eh=0
      SNO50COV=0.5!表层岩溶带延迟出流指数系数
      k_f=0.
      spcon_bsn=0.5!表层岩溶带传递给管道流，这里与基质流传递方式不一样，线性传递
      gw_karst_delaye=0
      spexp_bsn=0.5!补给地下岩溶基质的延迟系数
      gw_alpha_bfe=0
      n_updis=0.5!!地下岩溶基质出流的延迟系数
      k_fm=0
      p_updis=0.5!管道流的出流系数
      alpha_FM=0
      pperco_bsn=0.5!管道流的比流量指数系数
      k_sc=0
      phoskd_bsn=0.5!基质和管道的交换系数
      alpha_sc=0
      psp_bsn=0.5!基质和管道的交换指数      
      karst_he=1
      rcn_sub_bsn=1 !表层岩溶带初始水深
      karst_hs=1
      nperco_bsn=1!基质初始水深
      karst_hf=1
      eros_spl=1!管道初始水深
      IC_karst=0.1
      karst_type=0
      sur_rchrg_karst=0.
      sur_karst_height=0.
      sur_connect=0.
      Q_sur_sec=0.
      Q_sur_eh=0.
      Q_karst_S=0.
      gw_karst_hru=0
      Q_karst_F=0.
      karst_hru=1
      gw_karst_h1=0.
      gw_karst_h2=0.
      karst_FM_height=0.
      pre_pcp=0.
      IWAPI=0.
      WAPI=0.
      IWAPI_1=0.
      WAPI_1=0.      
      wapi_all=0.  
      wapi_all_1=0.
      !wapi_all=0.
      api_id=0
      origin_parm=0.
      origin_gw1_parm=0.
      origin_gw2_parm   =0. 
      origin_esco_parm=0.
      origin_brt_parm=0.
      origin_Ia=0.
      surf_Ia=0.
      origin_ch_n=0.
      origin_ch_k=0.
      karst_sc1=0.
      Matrix_conduct_flow_exchange=0.
      Q_karst_S1=0.
      Q_karst_F1=0.
      return
      end