       subroutine zero_karst


      use parm
      
      alpha_delay1=0
      sftmp=1!������������ܴ����ӳ�ϵ�������Ǳ�������������������һ���ӾͲ��������ܴ��У�������ʱ
      alpha_delay2=0
      smtmp=1!������ܴ����ݸ����ʵĹ��̣�������ʱ
      sur_karst_h1=0
      smfmn=1!������ܴ������ӳ�������ֵˮ��
      sur_karst_h2=0
      smfmx=2!������ܴ���������������ֵˮ��
      k_sec=0
      timp=0.5!������ܴ��������ĳ���������ϵ��(���Գ���)
      k_eh=0.
      SNOCOVMX=0.5!������ܴ��ӳٳ���������ϵ��
      alpha_eh=0
      SNO50COV=0.5!������ܴ��ӳٳ���ָ��ϵ��
      k_f=0.
      spcon_bsn=0.5!������ܴ����ݸ��ܵ�������������������ݷ�ʽ��һ�������Դ���
      gw_karst_delaye=0
      spexp_bsn=0.5!�����������ܻ��ʵ��ӳ�ϵ��
      gw_alpha_bfe=0
      n_updis=0.5!!�������ܻ��ʳ������ӳ�ϵ��
      k_fm=0
      p_updis=0.5!�ܵ����ĳ���ϵ��
      alpha_FM=0
      pperco_bsn=0.5!�ܵ����ı�����ָ��ϵ��
      k_sc=0
      phoskd_bsn=0.5!���ʺ͹ܵ��Ľ���ϵ��
      alpha_sc=0
      psp_bsn=0.5!���ʺ͹ܵ��Ľ���ָ��      
      karst_he=1
      rcn_sub_bsn=1 !������ܴ���ʼˮ��
      karst_hs=1
      nperco_bsn=1!���ʳ�ʼˮ��
      karst_hf=1
      eros_spl=1!�ܵ���ʼˮ��
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