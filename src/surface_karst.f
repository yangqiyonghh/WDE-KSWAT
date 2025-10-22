      subroutine surface_karst!表层岩溶带

      use parm

      
      
      real :: xx, r2,sur_karst_height1
      integer::j
      j = ihru
      if (curyr==1 .and. i==1)sur_karst_height(j)=karst_he!初始水深
          sur_connect=min((sur_karst_height(j)-sur_karst_h1)/(sur_karst_h2-sur_karst_h1),1.0)
          sur_connect=max(sur_connect,0.01) 
      
      if (hru_slp(j)>0.5)then

          sur_rchrg_karst(j)=sur_rchrg_karst(j)*2.71**(-1/alpha_delay1)+(1-2.71**(-1/alpha_delay1))*inflpcp
          sur_karst_height(j)=sur_karst_height(j)+sur_rchrg_karst(j)
        
          if (sur_karst_height(j)>sur_karst_h2)then
              Q_sur_sec(j)=k_sec*(max(sur_karst_height(j)-sur_karst_h2,0.0))

              !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh
              Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!延迟出流
              Q_sur_eh(j)=min(Q_sur_eh(j),sur_karst_height(j)-sur_karst_h1)
                    if (Q_sur_eh(j)>sur_karst_height(j)-sur_karst_h1 .and. sur_karst_height(j)-sur_karst_h1>0)Q_sur_eh(j)=
     &         sur_karst_height(j)-sur_karst_h1              
              !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!更新出流后表层岩溶带的储水量
          else
               if (sur_karst_height(j)>sur_karst_h1)then!只有延迟出流
                   Q_sur_sec(j)=0
                    !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!延迟出流
                   Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!延迟出流 
              
                   !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!更新出流后表层岩溶带的储水量
               else
                   Q_sur_sec(j)=0
                   Q_sur_eh(j)=0
                   !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!更新出流后表层岩溶带的储水量
               endif
          endif

          !Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)-sur_karst_h1),0.0)
          Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)),0.0)          

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
          
      else


          sur_rchrg_karst(j)=sur_rchrg_karst(j)*2.71**(-1/alpha_delay1)+(1-2.71**(-1/alpha_delay1))*sepbtm(j)!
          sur_karst_height(j)=sur_karst_height(j)+sur_rchrg_karst(j)!

          if (sur_karst_height(j)>sur_karst_h2)then!
              Q_sur_sec(j)=k_sec*(max(sur_karst_height(j)-sur_karst_h2,0.0))!二次流
              !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!
              Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh! 
                    if (Q_sur_eh(j)>sur_karst_height(j)-sur_karst_h1 .and. sur_karst_height(j)-sur_karst_h1>=0)Q_sur_eh(j)=
     &         sur_karst_height(j)-sur_karst_h1              
              !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!
          else!
               if (sur_karst_height(j)>sur_karst_h1)then!
                   Q_sur_sec(j)=0
                    !Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0)/(sur_karst_h2-sur_karst_h1))**alpha_eh!
                    Q_sur_eh(j)=k_eh*sur_connect*(max(sur_karst_height(j)-sur_karst_h1,0.0))**alpha_eh!
                    Q_sur_eh(j)=min(Q_sur_eh(j),sur_karst_height(j)-sur_karst_h1)

                    !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!
               else!
                   Q_sur_sec(j)=0
                   Q_sur_eh(j)=0
                   !sur_karst_height(j)=max(sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j),0.0)!
               endif

          endif 


          !Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)-sur_karst_h1),0.0)
          Q_karst_S(j)=max(Q_karst_S(j)*2.71**(-1/alpha_delay2)+(1-2.71**(-1/alpha_delay2))*(sur_karst_height(j)),0.0)  
          
          !
          if (gw_karst_hru(inum1)==1)Q_karst_F(j)=k_f*max(sur_karst_height(j),0.0)
          sur_karst_height1=sur_karst_height(j)-Q_sur_sec(j)-Q_sur_eh(j)-Q_karst_S(j)-Q_karst_F(j)!
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
   
      endif

   
      
      return
      end