      subroutine karst_parm
      !�������ͬ���ܵĲ���
      use parm
      integer :: j


      j = 0


      do j = 1, nhru
      !!���õ��ʶ�����
      !!��ʼ��������ʱ����λ��
      
      !1����̼��������϶�ܶ�ˮ��ʯ������϶�ܶ�ˮ
      !2����̼�����Ҽ���м����϶�ܶ�ˮ
      !3����̼��������϶�ܶ�ˮ����������϶�׶���϶ˮ
      karst_type(j)=1
      if (karst_type(j)==1)alpha_delay1=sftmp!(0.1-5)������������ܴ����ӳ�ϵ�������Ǳ�������������������һ���ӾͲ��������ܴ��У�������ʱ,karst_type=1
      if (karst_type(j)==2)alpha_delay1=smtmp!(0.1-5)������������ܴ����ӳ�ϵ�������Ǳ�������������������һ���ӾͲ��������ܴ��У�������ʱ,karst_type=2
      if (karst_type(j)==3)alpha_delay1=smfmx!(0.1-5)������������ܴ����ӳ�ϵ�������Ǳ�������������������һ���ӾͲ��������ܴ��У�������ʱ,karst_type=3
      
      if (karst_type(j)==1)alpha_delay2=smfmn!(0.1-5)������ܴ����ݸ����ʵĹ��̣�������ʱ,karst_type=1
      if (karst_type(j)==2)alpha_delay2=timp!(0.1-5)������ܴ����ݸ����ʵĹ��̣�������ʱ,karst_type=2
      if (karst_type(j)==3)alpha_delay2=snocovmx!(0.1-5)������ܴ����ݸ����ʵĹ��̣�������ʱ,karst_type=3
      
      if (karst_type(j)==1)sur_karst_h1=sno50cov!(25-100)������ܴ������ӳ�������ֵˮ��,karst_type=1
      if (karst_type(j)==2)sur_karst_h1=spcon_bsn!(25-100)������ܴ������ӳ�������ֵˮ��,karst_type=2
      if (karst_type(j)==3)sur_karst_h1=spexp_bsn!(25-100)������ܴ������ӳ�������ֵˮ��,karst_type=3
      
      if (karst_type(j)==1)sur_karst_h2=sno50cov+RCN11!(1-50)������ܴ���������������ֵˮ��,karst_type=1
      if (karst_type(j)==2)sur_karst_h2=spcon_bsn+cmn_bsn!(1-50)������ܴ���������������ֵˮ��,karst_type=2
      if (karst_type(j)==3)sur_karst_h2=spexp_bsn+n_updis!(1-50)������ܴ���������������ֵˮ��,karst_type=3
      
      if (karst_type(j)==1)k_sec=p_updis!(0.0001-0.1)������ܴ��������ĳ���������ϵ��(���Գ���),karst_type=1
      if (karst_type(j)==2)k_sec=nperco_bsn!(0.0001-0.1)������ܴ��������ĳ���������ϵ��(���Գ���),karst_type=2
      if (karst_type(j)==3)k_sec=pperco_bsn!(0.0001-0.1)������ܴ��������ĳ���������ϵ��(���Գ���),karst_type=3
      
      
      if (karst_type(j)==1)k_eh=phoskd_bsn!(0.0001-0.1)������ܴ��ӳٳ���������ϵ��,karst_type=1
      if (karst_type(j)==2)k_eh=psp_bsn!(0.0001-0.1)������ܴ��ӳٳ���������ϵ��,karst_type=2
      if (karst_type(j)==3)k_eh=rsdco!(0.0001-0.1)������ܴ��ӳٳ���������ϵ��,karst_type=3
      
      if (karst_type(j)==1)alpha_eh=percop!(0.5-3)������ܴ��ӳٳ���ָ��ϵ��,karst_type=1
      if (karst_type(j)==2)alpha_eh=wdpq!(0.5-3)������ܴ��ӳٳ���ָ��ϵ��,karst_type=2
      if (karst_type(j)==3)alpha_eh=wgpq!(0.5-3)������ܴ��ӳٳ���ָ��ϵ��,karst_type=3
      
      if (karst_type(j)==1)k_f=wdlpq!(0.0001-0.1)������ܴ����ݸ��ܵ�������������������ݷ�ʽ��һ�������Դ���,karst_type=1
      if (karst_type(j)==2)k_f=wglpq!(0.0001-0.1)������ܴ����ݸ��ܵ�������������������ݷ�ʽ��һ�������Դ���,karst_type=2      
      if (karst_type(j)==3)k_f=wdps!(0.0001-0.1)������ܴ����ݸ��ܵ�������������������ݷ�ʽ��һ�������Դ���,karst_type=3
      
      if (karst_type(j)==1)gw_karst_delaye=wgps!(0.1-5)�����������ܻ��ʵ��ӳ�ϵ��,karst_type=1
      if (karst_type(j)==2)gw_karst_delaye=wdlps!(0.1-5)�����������ܻ��ʵ��ӳ�ϵ��,karst_type=2      
      if (karst_type(j)==3)gw_karst_delaye=wglps!(0.1-5)�����������ܻ��ʵ��ӳ�ϵ��,karst_type=3
      
      if (karst_type(j)==1)gw_alpha_bfe=bactkdq!(0.1-5)�������ܻ��ʳ������ӳ�ϵ��,karst_type=1
      if (karst_type(j)==2)gw_alpha_bfe=thbact!(0.1-5)�������ܻ��ʳ������ӳ�ϵ��,karst_type=2      
      if (karst_type(j)==3)gw_alpha_bfe=wof_p!(0.1-5)�������ܻ��ʳ������ӳ�ϵ��,karst_type=3
       
      if (karst_type(j)==1)k_fm=wof_lp!(0.0001-0.1)�ܵ����ĳ���ϵ��,karst_type=1
      if (karst_type(j)==2)k_fm=wdpf!(0.0001-0.1)�ܵ����ĳ���ϵ��,karst_type=2      
      if (karst_type(j)==3)k_fm=wgpf!(0.0001-0.1)�ܵ����ĳ���ϵ��,karst_type=3
      
      if (karst_type(j)==1)alpha_FM=wdlpf!(0.5-3)�ܵ����ı�����ָ��ϵ��?,karst_type=1
      if (karst_type(j)==2)alpha_FM=wglpf!(0.5-3)�ܵ����ı�����ָ��ϵ��?,karst_type=2      
      if (karst_type(j)==2)alpha_FM=cncoef!(0.5-3)�ܵ����ı�����ָ��ϵ��?,karst_type=3
      
      if (karst_type(j)==1)k_sc=cdn_bsn!(0.0001-0.05)���ʺ͹ܵ��Ľ���ϵ��,karst_type=1
      if (karst_type(j)==2)k_sc=sdnco_bsn!(0.0001-0.05)���ʺ͹ܵ��Ľ���ϵ��,karst_type=2      
      if (karst_type(j)==3)k_sc=bact_swf!(0.0001-0.05)���ʺ͹ܵ��Ľ���ϵ��,karst_type=3
      
      if (karst_type(j)==1)alpha_sc=bactmx!(0.1-1)���ʺ͹ܵ��Ľ���ָ��,karst_type=1
      if (karst_type(j)==2)alpha_sc=bactminlp!(0.1-1)���ʺ͹ܵ��Ľ���ָ��,karst_type=2      
      alpha_sc=bactminp!(0.1-1)���ʺ͹ܵ��Ľ���ָ��,karst_type=3
      
      if (karst_type(j)==1)karst_he=wdlprch !(1-25)������ܴ���ʼˮ��,karst_type=1
      if (karst_type(j)==2)karst_he=wdprch !(1-25)������ܴ���ʼˮ��,karst_type=2      
      if (karst_type(j)==3)karst_he=wdlpres !(1-25)������ܴ���ʼˮ��,karst_type=3
      
      if (karst_type(j)==1)karst_hs=wdpres!(1-25)���ʳ�ʼˮ��,karst_type=1
      if (karst_type(j)==2)karst_hs=tb_adj!(1-25)���ʳ�ʼˮ��,karst_type=2
      if (karst_type(j)==3)karst_hs=depimp_bsn!(1-25)���ʳ�ʼˮ��,karst_type=3      
      
      if (karst_type(j)==1)karst_hf=ddrain_bsn!(1-25)�ܵ���ʼˮ��,karst_type=1
      if (karst_type(j)==2)karst_hf=tdrain_bsn!(1-25)�ܵ���ʼˮ��,karst_type=2
      if (karst_type(j)==3)karst_hf=gdrain_bsn!(1-25)�ܵ���ʼˮ��,karst_type=3      
            
      end do
      

      return
      end