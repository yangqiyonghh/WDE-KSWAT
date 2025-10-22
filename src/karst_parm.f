      subroutine karst_parm
      !分配给不同岩溶的参数
      use parm
      integer :: j


      j = 0


      do j = 1, nhru

      karst_type(j)=1
      if (karst_type(j)==1)alpha_delay1=sftmp!(0.1-5)
      if (karst_type(j)==2)alpha_delay1=smtmp!(0.1-5)
      if (karst_type(j)==3)alpha_delay1=smfmx!(0.1-5)
      
      if (karst_type(j)==1)alpha_delay2=smfmn
      if (karst_type(j)==2)alpha_delay2=timp
      if (karst_type(j)==3)alpha_delay2=snocovmx
      
      if (karst_type(j)==1)sur_karst_h1=sno50cov
      if (karst_type(j)==2)sur_karst_h1=spcon_bsn
      if (karst_type(j)==3)sur_karst_h1=spexp_bsn
      
      if (karst_type(j)==1)sur_karst_h2=sno50cov+RCN11
      if (karst_type(j)==2)sur_karst_h2=spcon_bsn+cmn_bsn
      if (karst_type(j)==3)sur_karst_h2=spexp_bsn+n_updis
      
      if (karst_type(j)==1)k_sec=p_updis
      if (karst_type(j)==2)k_sec=nperco_bsn
      if (karst_type(j)==3)k_sec=pperco_bsn
      
      
      if (karst_type(j)==1)k_eh=phoskd_bsn
      if (karst_type(j)==2)k_eh=psp_bsn
      if (karst_type(j)==3)k_eh=rsdco
      
      if (karst_type(j)==1)alpha_eh=percop
      if (karst_type(j)==2)alpha_eh=wdpq
      if (karst_type(j)==3)alpha_eh=wgpq
      
      if (karst_type(j)==1)k_f=wdlpq
      if (karst_type(j)==2)k_f=wglpq     
      if (karst_type(j)==3)k_f=wdps
      
      if (karst_type(j)==1)gw_karst_delaye=wgps
      if (karst_type(j)==2)gw_karst_delaye=wdlps   
      if (karst_type(j)==3)gw_karst_delaye=wglps
      
      if (karst_type(j)==1)gw_alpha_bfe=bactkdq
      if (karst_type(j)==2)gw_alpha_bfe=thbact     
      if (karst_type(j)==3)gw_alpha_bfe=wof_p
       
      if (karst_type(j)==1)k_fm=wof_lp
      if (karst_type(j)==2)k_fm=wdpf    
      if (karst_type(j)==3)k_fm=wgpf
      
      if (karst_type(j)==1)alpha_FM=wdlpf
      if (karst_type(j)==2)alpha_FM=wglpf     
      if (karst_type(j)==2)alpha_FM=cncoef
      
      if (karst_type(j)==1)k_sc=cdn_bsn
      if (karst_type(j)==2)k_sc=sdnco_bsn     
      if (karst_type(j)==3)k_sc=bact_swf
      
      if (karst_type(j)==1)alpha_sc=bactmx
      if (karst_type(j)==2)alpha_sc=bactminlp  
      alpha_sc=bactminp
      
      if (karst_type(j)==1)karst_he=wdlprch 
      if (karst_type(j)==2)karst_he=wdprch     
      if (karst_type(j)==3)karst_he=wdlpres 
      
      if (karst_type(j)==1)karst_hs=wdpres
      if (karst_type(j)==2)karst_hs=tb_adj
      if (karst_type(j)==3)karst_hs=depimp_bsn    
      
      if (karst_type(j)==1)karst_hf=ddrain_bsn
      if (karst_type(j)==2)karst_hf=tdrain_bsn
      if (karst_type(j)==3)karst_hf=gdrain_bsn     
            
      end do
      

      return
      end