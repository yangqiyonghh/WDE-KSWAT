      subroutine gw_karst_num!���º�������������

      use parm

      
      if (sub_karst==1 )then !1��������,
          gw_karst_hru(2)=1!�𶴵��º�
          gw_karst_hru(25)=1!���̵��º�
          gw_karst_hru(27)=1!���ϵ��º�
          gw_karst_hru(15)=1!�˿����ºӶ�
          gw_karst_hru(16)=1!�˿����ºӶ�
          gw_karst_hru(50)=1!�˿����ºӶ�
          gw_karst_hru(51)=1!�˿����ºӶ�
          gw_karst_hru(55)=1!�˿����º�+��ͷ���º�
          gw_karst_hru(53)=1!�˿����ºӶ�
          gw_karst_hru(52)=1!�˿����ºӶ�
          gw_karst_hru(66)=1!��ͷ���º�
          gw_karst_hru(20)=1!������º�
          gw_karst_hru(73)=1!��˳���º�
          gw_karst_hru(21)=1!��˳���º�
          gw_karst_hru(6)=1!��ɽ���º�+�Ŷɵ��º�
          gw_karst_hru(59)=1!��ĥ���º�
          gw_karst_hru(17)=1!��ĥ���º�
          gw_karst_hru(74)=1!��̶���º�
          gw_karst_hru(13)=1!��̶���º�
          gw_karst_hru(18)=1!��̶���º�
          gw_karst_hru(77)=1!�������º�
          gw_karst_hru(80)=1!�������º�
          gw_karst_hru(8)=1!���ҵ��º�
          gw_karst_hru(84)=1!ͥ�ɵ��º�
          gw_karst_hru(14)=1 !�������º�
          gw_karst_hru(85)=1 !�������º�
      else!0:�轭����
          gw_karst_hru(85)=1 
      endif
      
    
      return
      end