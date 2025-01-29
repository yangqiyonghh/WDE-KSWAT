      subroutine surface_karst!������ܴ�

      use parm

      
      
      real :: xx, r2,sur_karst_height1
      integer::j
      j = ihru
      if (curyr==1 .and. i==1)sur_karst_height(j)=karst_he!��ʼˮ��
          sur_connect=min((sur_karst_height(j)-sur_karst_h1)/(sur_karst_h2-sur_karst_h1),1.0)
          sur_connect=max(sur_connect,0.01) 
      
      if (hru_slp(j)>0.5)then!ʯĮ������
          !���������ܴ��Ĳ�����
          sur_rchrg_karst(j)=sur_rchrg_karst(j)*2.71**(-1/alpha_delay1)+(1-2.71**(-1/alpha_delay1))*inflpcp!alpha_delay1Ϊ�ӳ�ϵ��
          sur_karst_height(j)=sur_karst_height(j)+sur_rchrg_karst(j)!���±�����ܴ��Ĵ�ˮ��
        
          if (sur_karst_height(j)>sur_karst_h2)then!�ж��������ӳٳ���
              Q_sur_sec(j)=k_sec*(max(sur_karst_height(j)-sur_karst_h2,0.0))!������
              !sur_connect������ܴ�ˮ����ͨ��
              !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!�ӳٳ���
              Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!�ӳٳ���
              Q_sur_eh(j)=min(Q_sur_eh(j),sur_karst_height(j)-sur_karst_h1)
                    if (Q_sur_eh(j)>sur_karst_height(j)-sur_karst_h1 .and. sur_karst_height(j)-sur_karst_h1>0)Q_sur_eh(j)=
     &         sur_karst_height(j)-sur_karst_h1              
              !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!���³����������ܴ��Ĵ�ˮ��
          else!�޶�����
               if (sur_karst_height(j)>sur_karst_h1)then!ֻ���ӳٳ���
                   Q_sur_sec(j)=0
                    !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!�ӳٳ���
                   Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!�ӳٳ��� 
              
                   !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!���³����������ܴ��Ĵ�ˮ��
               else!��ˮ��
                   Q_sur_sec(j)=0
                   Q_sur_eh(j)=0
                   !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!���³����������ܴ��Ĵ�ˮ��
               endif
          endif
      !   IF (Q_sur_eh(j)>10000) print*,i,j,Q_sur_eh(j),sur_connect,k_eh,sur_karst_height(j),sur_karst_h1,sur_karst_h2,alpha_eh,
      !&     karst_he,sur_rchrg_karst(j),1
          !���ʵĲ�����
          !Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)-sur_karst_h1),0.0)
          Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)),0.0)          
          !���¹ܵ��Ĳ��������е��ºӵĲŽ��м���
          !if (gw_karst_hru(inum1)==1)Q_karst_F(j)=k_f*max(sur_karst_height(j)-sur_karst_h1,0.0)
          if (gw_karst_hru(inum1)==1)Q_karst_F(j)=k_f*max(sur_karst_height(j),0.0)
          sur_karst_height1=sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j)-Q_karst_S(j)-Q_karst_F(j)!���³����������ܴ��Ĵ�ˮ��
          if (sur_karst_height1-sur_karst_h1<=0 )then
              Q_sur_sec(j)=max((sur_karst_height(j)-sur_karst_h1)*Q_sur_sec(j)/(Q_sur_sec(j)+Q_sur_eh(j)+Q_karst_S(j)+Q_karst_F(j)
     &         +1e-6),0.0)
              Q_sur_eh(j)=max((sur_karst_height(j)-sur_karst_h1)*Q_sur_eh(j)/(Q_sur_sec(j)+Q_sur_eh(j)+Q_karst_S(j)+Q_karst_F(j)
     &              +1e-6),0.)
              Q_karst_S(j)=max((sur_karst_height(j)-sur_karst_h1)*Q_karst_S(j)/(Q_sur_sec(j)+Q_sur_eh(j)+Q_karst_S(j)+Q_karst_F(j)
     &              +1e-6),0.)
              Q_karst_F(j)=max((sur_karst_height(j)-sur_karst_h1)*Q_karst_F(j)/(Q_sur_sec(j)+Q_sur_eh(j)+Q_karst_S(j)+Q_karst_F(j)
     &              +1e-6),0.)  
              sur_karst_height(j)=sur_karst_h1
              
          else
              sur_karst_height(j)=sur_karst_height1
          end if
          
      else!��ʯĮ������
          !if (j==378)print*,sur_rchrg_karst(j),sepbtm(j),i,1
          !���������ܴ��Ĳ�����
          sur_rchrg_karst(j)=sur_rchrg_karst(j)*2.71**(-1/alpha_delay1)+(1-2.71**(-1/alpha_delay1))*sepbtm(j)!alpha_delayΪ�ӳ�ϵ��
          sur_karst_height(j)=sur_karst_height(j)+sur_rchrg_karst(j)!���±�����ܴ��Ĵ�ˮ��
         !if (j==378)print*,sur_rchrg_karst(j),sur_karst_height(j),i,2
          if (sur_karst_height(j)>sur_karst_h2)then!�ж��������ӳٳ���
              Q_sur_sec(j)=k_sec*(max(sur_karst_height(j)-sur_karst_h2,0.0))!������
              !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!�ӳٳ���
              Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!�ӳٳ���  
                    if (Q_sur_eh(j)>sur_karst_height(j)-sur_karst_h1 .and. sur_karst_height(j)-sur_karst_h1>=0)Q_sur_eh(j)=
     &         sur_karst_height(j)-sur_karst_h1              
              !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!���³����������ܴ��Ĵ�ˮ��
          else!�޶�����
               if (sur_karst_height(j)>sur_karst_h1)then!ֻ���ӳٳ���
                   Q_sur_sec(j)=0
                    !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!�ӳٳ���
                    Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!�ӳٳ���
                    Q_sur_eh(j)=min(Q_sur_eh(j),sur_karst_height(j)-sur_karst_h1)

                    !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!���³����������ܴ��Ĵ�ˮ��
               else!��ˮ��
                   Q_sur_sec(j)=0
                   Q_sur_eh(j)=0
                   !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!���³����������ܴ��Ĵ�ˮ��
               endif

          endif 
          !���ʵĲ�����

          !Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)-sur_karst_h1),0.0)
          Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)),0.0)  
          
          !���¹ܵ��Ĳ��������е��ºӵĲŽ��м���
          if (gw_karst_hru(inum1)==1)Q_karst_F(j)=k_f*max(sur_karst_height(j),0.0)
          sur_karst_height1=sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j)-Q_karst_S(j)-Q_karst_F(j)!���³����������ܴ��Ĵ�ˮ��
          if (sur_karst_height1-sur_karst_h1<=0)then
              Q_sur_sec(j)=max((sur_karst_height(j)-sur_karst_h1)*Q_sur_sec(j)/(Q_sur_sec(j)+Q_sur_eh(j)+Q_karst_S(j)+Q_karst_F(j)
     &         +1e-6),0.0)
              Q_sur_eh(j)=max((sur_karst_height(j)-sur_karst_h1)*Q_sur_eh(j)/(Q_sur_sec(j)+Q_sur_eh(j)+Q_karst_S(j)+Q_karst_F(j)
     &              +1e-6),0.)
              Q_karst_S(j)=max((sur_karst_height(j)-sur_karst_h1)*Q_karst_S(j)/(Q_sur_sec(j)+Q_sur_eh(j)+Q_karst_S(j)+Q_karst_F(j)
     &              +1e-6),0.)
              Q_karst_F(j)=max((sur_karst_height(j)-sur_karst_h1)*Q_karst_F(j)/(Q_sur_sec(j)+Q_sur_eh(j)+Q_karst_S(j)+Q_karst_F(j)
     &              +1e-6),0.)  
              sur_karst_height(j)=sur_karst_h1
          else
              sur_karst_height(j)=sur_karst_height1
          end if
      !    IF (Q_sur_eh(j)>10000) print*,i,j,Q_sur_eh(j),sur_connect,k_eh,sur_karst_height(j),sur_karst_h1,sur_karst_h2,alpha_eh,
      !&  karst_he,sur_rchrg_karst(j),2         
      endif

      !if (i==22 .and. j==542)
      
      return
      end