      subroutine gw_karst_num!地下河所在子流域编号

      use parm

      
      if (sub_karst==1 )then !1：刁江流域,
          gw_karst_hru(2)=1!金洞地下河
          gw_karst_hru(25)=1!拉盘地下河
          gw_karst_hru(27)=1!拉南地下河
          gw_karst_hru(15)=1!八况地下河段
          gw_karst_hru(16)=1!八况地下河段
          gw_karst_hru(50)=1!八况地下河段
          gw_karst_hru(51)=1!八况地下河段
          gw_karst_hru(55)=1!八况地下河+龙头地下河
          gw_karst_hru(53)=1!八况地下河段
          gw_karst_hru(52)=1!八况地下河段
          gw_karst_hru(66)=1!龙头地下河
          gw_karst_hru(20)=1!拉竹地下河
          gw_karst_hru(73)=1!永顺地下河
          gw_karst_hru(21)=1!永顺地下河
          gw_karst_hru(6)=1!尖山地下河+九渡地下河
          gw_karst_hru(59)=1!九磨地下河
          gw_karst_hru(17)=1!九磨地下河
          gw_karst_hru(74)=1!索潭地下河
          gw_karst_hru(13)=1!索潭地下河
          gw_karst_hru(18)=1!索潭地下河
          gw_karst_hru(77)=1!加若地下河
          gw_karst_hru(80)=1!加若地下河
          gw_karst_hru(8)=1!拉烈地下河
          gw_karst_hru(84)=1!庭律地下河
          gw_karst_hru(14)=1 !百赖地下河
          gw_karst_hru(85)=1 !百赖地下河
      else!0:滂江流域
          gw_karst_hru(85)=1 
      endif
      
    
      return
      end