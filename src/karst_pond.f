      subroutine karst_pond!�����ݵ�

      use parm

      
      
      real :: xx, r2
      integer::j
      j = ihru
      
      if (hru_slp(j)>0.5)then!ʯĮ������
          (inflpcp+qday)
          
          
      else!��ʯĮ������
          (inflpcp)
          
      endif
      
          

      return
      end