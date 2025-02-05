      subroutine karst_pond!岩溶洼地

      use parm

      
      
      real :: xx, r2
      integer::j
      j = ihru
      
      if (hru_slp(j)>0.5)then!石漠化区域
          (inflpcp+qday)
          
          
      else!非石漠化区域
          (inflpcp)
          
      endif
      
          

      return
      end