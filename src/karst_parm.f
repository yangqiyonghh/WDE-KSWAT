      subroutine karst_parm
      !分配给不同岩溶的参数
      use parm
      integer :: j


      j = 0


      do j = 1, nhru
      !!可用的率定参数
      !!初始参数的暂时设置位置
      
      !1代表碳酸盐岩裂隙溶洞水：石灰岩裂隙溶洞水
      !2代表碳酸盐岩夹碎屑岩裂隙溶洞水
      !3代表碳酸盐岩裂隙溶洞水：白云岩裂隙孔洞裂隙水
      karst_type(j)=1
      if (karst_type(j)==1)alpha_delay1=sftmp!(0.1-5)下渗给表层岩溶带的延迟系数，就是比如土壤下渗量并不是一下子就补给到岩溶带中，存在滞时,karst_type=1
      if (karst_type(j)==2)alpha_delay1=smtmp!(0.1-5)下渗给表层岩溶带的延迟系数，就是比如土壤下渗量并不是一下子就补给到岩溶带中，存在滞时,karst_type=2
      if (karst_type(j)==3)alpha_delay1=smfmx!(0.1-5)下渗给表层岩溶带的延迟系数，就是比如土壤下渗量并不是一下子就补给到岩溶带中，存在滞时,karst_type=3
      
      if (karst_type(j)==1)alpha_delay2=smfmn!(0.1-5)表层岩溶带传递给基质的过程，存在滞时,karst_type=1
      if (karst_type(j)==2)alpha_delay2=timp!(0.1-5)表层岩溶带传递给基质的过程，存在滞时,karst_type=2
      if (karst_type(j)==3)alpha_delay2=snocovmx!(0.1-5)表层岩溶带传递给基质的过程，存在滞时,karst_type=3
      
      if (karst_type(j)==1)sur_karst_h1=sno50cov!(25-100)表层岩溶带发生延迟流的阈值水深,karst_type=1
      if (karst_type(j)==2)sur_karst_h1=spcon_bsn!(25-100)表层岩溶带发生延迟流的阈值水深,karst_type=2
      if (karst_type(j)==3)sur_karst_h1=spexp_bsn!(25-100)表层岩溶带发生延迟流的阈值水深,karst_type=3
      
      if (karst_type(j)==1)sur_karst_h2=sno50cov+RCN11!(1-50)表层岩溶带发生二次流的阈值水深,karst_type=1
      if (karst_type(j)==2)sur_karst_h2=spcon_bsn+cmn_bsn!(1-50)表层岩溶带发生二次流的阈值水深,karst_type=2
      if (karst_type(j)==3)sur_karst_h2=spexp_bsn+n_updis!(1-50)表层岩溶带发生二次流的阈值水深,karst_type=3
      
      if (karst_type(j)==1)k_sec=p_updis!(0.0001-0.1)表层岩溶带二次流的出流比流量系数(线性出流),karst_type=1
      if (karst_type(j)==2)k_sec=nperco_bsn!(0.0001-0.1)表层岩溶带二次流的出流比流量系数(线性出流),karst_type=2
      if (karst_type(j)==3)k_sec=pperco_bsn!(0.0001-0.1)表层岩溶带二次流的出流比流量系数(线性出流),karst_type=3
      
      
      if (karst_type(j)==1)k_eh=phoskd_bsn!(0.0001-0.1)表层岩溶带延迟出流比流量系数,karst_type=1
      if (karst_type(j)==2)k_eh=psp_bsn!(0.0001-0.1)表层岩溶带延迟出流比流量系数,karst_type=2
      if (karst_type(j)==3)k_eh=rsdco!(0.0001-0.1)表层岩溶带延迟出流比流量系数,karst_type=3
      
      if (karst_type(j)==1)alpha_eh=percop!(0.5-3)表层岩溶带延迟出流指数系数,karst_type=1
      if (karst_type(j)==2)alpha_eh=wdpq!(0.5-3)表层岩溶带延迟出流指数系数,karst_type=2
      if (karst_type(j)==3)alpha_eh=wgpq!(0.5-3)表层岩溶带延迟出流指数系数,karst_type=3
      
      if (karst_type(j)==1)k_f=wdlpq!(0.0001-0.1)表层岩溶带传递给管道流，这里与基质流传递方式不一样，线性传递,karst_type=1
      if (karst_type(j)==2)k_f=wglpq!(0.0001-0.1)表层岩溶带传递给管道流，这里与基质流传递方式不一样，线性传递,karst_type=2      
      if (karst_type(j)==3)k_f=wdps!(0.0001-0.1)表层岩溶带传递给管道流，这里与基质流传递方式不一样，线性传递,karst_type=3
      
      if (karst_type(j)==1)gw_karst_delaye=wgps!(0.1-5)补给地下岩溶基质的延迟系数,karst_type=1
      if (karst_type(j)==2)gw_karst_delaye=wdlps!(0.1-5)补给地下岩溶基质的延迟系数,karst_type=2      
      if (karst_type(j)==3)gw_karst_delaye=wglps!(0.1-5)补给地下岩溶基质的延迟系数,karst_type=3
      
      if (karst_type(j)==1)gw_alpha_bfe=bactkdq!(0.1-5)地下岩溶基质出流的延迟系数,karst_type=1
      if (karst_type(j)==2)gw_alpha_bfe=thbact!(0.1-5)地下岩溶基质出流的延迟系数,karst_type=2      
      if (karst_type(j)==3)gw_alpha_bfe=wof_p!(0.1-5)地下岩溶基质出流的延迟系数,karst_type=3
       
      if (karst_type(j)==1)k_fm=wof_lp!(0.0001-0.1)管道流的出流系数,karst_type=1
      if (karst_type(j)==2)k_fm=wdpf!(0.0001-0.1)管道流的出流系数,karst_type=2      
      if (karst_type(j)==3)k_fm=wgpf!(0.0001-0.1)管道流的出流系数,karst_type=3
      
      if (karst_type(j)==1)alpha_FM=wdlpf!(0.5-3)管道流的比流量指数系数?,karst_type=1
      if (karst_type(j)==2)alpha_FM=wglpf!(0.5-3)管道流的比流量指数系数?,karst_type=2      
      if (karst_type(j)==2)alpha_FM=cncoef!(0.5-3)管道流的比流量指数系数?,karst_type=3
      
      if (karst_type(j)==1)k_sc=cdn_bsn!(0.0001-0.05)基质和管道的交换系数,karst_type=1
      if (karst_type(j)==2)k_sc=sdnco_bsn!(0.0001-0.05)基质和管道的交换系数,karst_type=2      
      if (karst_type(j)==3)k_sc=bact_swf!(0.0001-0.05)基质和管道的交换系数,karst_type=3
      
      if (karst_type(j)==1)alpha_sc=bactmx!(0.1-1)基质和管道的交换指数,karst_type=1
      if (karst_type(j)==2)alpha_sc=bactminlp!(0.1-1)基质和管道的交换指数,karst_type=2      
      alpha_sc=bactminp!(0.1-1)基质和管道的交换指数,karst_type=3
      
      if (karst_type(j)==1)karst_he=wdlprch !(1-25)表层岩溶带初始水深,karst_type=1
      if (karst_type(j)==2)karst_he=wdprch !(1-25)表层岩溶带初始水深,karst_type=2      
      if (karst_type(j)==3)karst_he=wdlpres !(1-25)表层岩溶带初始水深,karst_type=3
      
      if (karst_type(j)==1)karst_hs=wdpres!(1-25)基质初始水深,karst_type=1
      if (karst_type(j)==2)karst_hs=tb_adj!(1-25)基质初始水深,karst_type=2
      if (karst_type(j)==3)karst_hs=depimp_bsn!(1-25)基质初始水深,karst_type=3      
      
      if (karst_type(j)==1)karst_hf=ddrain_bsn!(1-25)管道初始水深,karst_type=1
      if (karst_type(j)==2)karst_hf=tdrain_bsn!(1-25)管道初始水深,karst_type=2
      if (karst_type(j)==3)karst_hf=gdrain_bsn!(1-25)管道初始水深,karst_type=3      
            
      end do
      

      return
      end