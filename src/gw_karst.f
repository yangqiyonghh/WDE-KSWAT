      subroutine gw_karst!地下河

      use parm

      

      real :: xx, r2
      integer::j
      j = ihru
      if (curyr==1 .and. i==1)karst_FM_height(j)=karst_hf!初始水深
      karst_FM_height(j)=karst_FM_height(j)+Q_karst_F(j)
      
      Q_karst_FM(j)=k_fm*(karst_FM_height(j))**alpha_FM!管道出流量
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


      if (karst_FM_height(j)>shallst(j) )then!
          
          karst_sc(j)=max(k_sc*(karst_FM_height(j)-shallst(j) )**alpha_sc,0.010)
          karst_sc(j)=min(karst_sc(j),0.5*(karst_FM_height(j)-shallst(j)))
          shallst(j)=shallst(j)+karst_sc(j)
          karst_FM_height(j)=max(karst_FM_height(j)-karst_sc(j),0.0)
      else!
          karst_sc(j)=max(k_sc*(shallst(j) -karst_FM_height(j))**alpha_sc,0.01)
          karst_sc(j)=min(karst_sc(j),0.5*(shallst(j)-karst_FM_height(j)))
          shallst(j)=max(shallst(j)-karst_sc(j),0.0)
          karst_FM_height(j)=karst_FM_height(j)+karst_sc(j)
          
      endif   

      return
      end