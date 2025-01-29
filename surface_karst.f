      subroutine surface_karst!表层岩溶带

      use parm

      
      
      real :: xx, r2,sur_karst_height1
      integer::j
      j = ihru
      if (curyr==1 .and. i==1)sur_karst_height(j)=karst_he!初始水深
          sur_connect=min((sur_karst_height(j)-sur_karst_h1)/(sur_karst_h2-sur_karst_h1),1.0)
          sur_connect=max(sur_connect,0.01) 
      
      if (hru_slp(j)>0.5)then!石漠化区域
          !计算表层岩溶带的补给量
          sur_rchrg_karst(j)=sur_rchrg_karst(j)*2.71**(-1/alpha_delay1)+(1-2.71**(-1/alpha_delay1))*inflpcp!alpha_delay1为延迟系数
          sur_karst_height(j)=sur_karst_height(j)+sur_rchrg_karst(j)!更新表层岩溶带的储水量
        
          if (sur_karst_height(j)>sur_karst_h2)then!有二次流和延迟出流
              Q_sur_sec(j)=k_sec*(max(sur_karst_height(j)-sur_karst_h2,0.0))!二次流
              !sur_connect表层岩溶带水文连通性
              !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!延迟出流
              Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!延迟出流
              Q_sur_eh(j)=min(Q_sur_eh(j),sur_karst_height(j)-sur_karst_h1)
                    if (Q_sur_eh(j)>sur_karst_height(j)-sur_karst_h1 .and. sur_karst_height(j)-sur_karst_h1>0)Q_sur_eh(j)=
     &         sur_karst_height(j)-sur_karst_h1              
              !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!更新出流后表层岩溶带的储水量
          else!无二次流
               if (sur_karst_height(j)>sur_karst_h1)then!只有延迟出流
                   Q_sur_sec(j)=0
                    !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!延迟出流
                   Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!延迟出流 
              
                   !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!更新出流后表层岩溶带的储水量
               else!无水流
                   Q_sur_sec(j)=0
                   Q_sur_eh(j)=0
                   !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!更新出流后表层岩溶带的储水量
               endif
          endif
      !   IF (Q_sur_eh(j)>10000) print*,i,j,Q_sur_eh(j),sur_connect,k_eh,sur_karst_height(j),sur_karst_h1,sur_karst_h2,alpha_eh,
      !&     karst_he,sur_rchrg_karst(j),1
          !基质的补给量
          !Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)-sur_karst_h1),0.0)
          Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)),0.0)          
          !地下管道的补给量，有地下河的才进行计算
          !if (gw_karst_hru(inum1)==1)Q_karst_F(j)=k_f*max(sur_karst_height(j)-sur_karst_h1,0.0)
          if (gw_karst_hru(inum1)==1)Q_karst_F(j)=k_f*max(sur_karst_height(j),0.0)
          sur_karst_height1=sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j)-Q_karst_S(j)-Q_karst_F(j)!更新出流后表层岩溶带的储水量
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
          
      else!非石漠化区域
          !if (j==378)print*,sur_rchrg_karst(j),sepbtm(j),i,1
          !计算表层岩溶带的补给量
          sur_rchrg_karst(j)=sur_rchrg_karst(j)*2.71**(-1/alpha_delay1)+(1-2.71**(-1/alpha_delay1))*sepbtm(j)!alpha_delay为延迟系数
          sur_karst_height(j)=sur_karst_height(j)+sur_rchrg_karst(j)!更新表层岩溶带的储水量
         !if (j==378)print*,sur_rchrg_karst(j),sur_karst_height(j),i,2
          if (sur_karst_height(j)>sur_karst_h2)then!有二次流和延迟出流
              Q_sur_sec(j)=k_sec*(max(sur_karst_height(j)-sur_karst_h2,0.0))!二次流
              !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!延迟出流
              Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!延迟出流  
                    if (Q_sur_eh(j)>sur_karst_height(j)-sur_karst_h1 .and. sur_karst_height(j)-sur_karst_h1>=0)Q_sur_eh(j)=
     &         sur_karst_height(j)-sur_karst_h1              
              !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!更新出流后表层岩溶带的储水量
          else!无二次流
               if (sur_karst_height(j)>sur_karst_h1)then!只有延迟出流
                   Q_sur_sec(j)=0
                    !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!延迟出流
                    Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!延迟出流
                    Q_sur_eh(j)=min(Q_sur_eh(j),sur_karst_height(j)-sur_karst_h1)

                    !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!更新出流后表层岩溶带的储水量
               else!无水流
                   Q_sur_sec(j)=0
                   Q_sur_eh(j)=0
                   !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!更新出流后表层岩溶带的储水量
               endif

          endif 
          !基质的补给量

          !Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)-sur_karst_h1),0.0)
          Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)),0.0)  
          
          !地下管道的补给量，有地下河的才进行计算
          if (gw_karst_hru(inum1)==1)Q_karst_F(j)=k_f*max(sur_karst_height(j),0.0)
          sur_karst_height1=sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j)-Q_karst_S(j)-Q_karst_F(j)!更新出流后表层岩溶带的储水量
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