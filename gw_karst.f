      subroutine gw_karst!���º�

      use parm

      
      !ֻ�и��������ڴ��ڵ��ºӲŴ��ڽ���
      real :: xx, r2
      integer::j
      j = ihru
      if (curyr==1 .and. i==1)karst_FM_height(j)=karst_hf!��ʼˮ��
      karst_FM_height(j)=karst_FM_height(j)+Q_karst_F(j)
      
      Q_karst_FM(j)=k_fm*(karst_FM_height(j))**alpha_FM!�ܵ�������
      Q_karst_FM(j)=min(Q_karst_FM(j),karst_FM_height(j))
      if (karst_FM_height(j) >= gw_karst_h2) then
        karst_FM_height(j) = karst_FM_height(j) - Q_karst_FM(j)
        if (karst_FM_height(j)< gw_karst_h2) then
         Q_karst_FM(j)= karst_FM_height(j) +Q_karst_FM(j) -gw_karst_h2
           karst_FM_height(j) = gw_karst_h2
        end if
       else
        Q_karst_FM(j) = 0.
       end if         
      
      !karst_FM_height(j)=max(karst_FM_height(j)-Q_karst_FM(j),0.0)
      !���ǻ��ʺ͹ܵ��Ľ���
      !if (j==13)print*,i,karst_FM_height(j),shallst(j),Q_karst_FM(j),11
      if (karst_FM_height(j)>shallst(j) )then!���¹ܵ����������
          
          karst_sc(j)=max(k_sc*(karst_FM_height(j)-shallst(j) )**alpha_sc,0.010)
          karst_sc(j)=min(karst_sc(j),0.5*(karst_FM_height(j)-shallst(j)))
          shallst(j)=shallst(j)+karst_sc(j)
          karst_FM_height(j)=max(karst_FM_height(j)-karst_sc(j),0.0)
      else!��������¹ܵ�����
          karst_sc(j)=max(k_sc*(shallst(j) -karst_FM_height(j))**alpha_sc,0.01)
          karst_sc(j)=min(karst_sc(j),0.5*(shallst(j)-karst_FM_height(j)))
          shallst(j)=max(shallst(j)-karst_sc(j),0.0)
          karst_FM_height(j)=karst_FM_height(j)+karst_sc(j)
          
      endif   
      !if (j==13)print*,i,karst_FM_height(j),karst_sc(j),22
      return
      end