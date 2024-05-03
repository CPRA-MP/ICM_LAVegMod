subroutine reset_coverages
    ! global arrays updated by subroutine:
    !      bare_old
    !      bare_new
    !      dead_Flt
    !
    ! This subroutine converts the previous model year's new bareground (bare_new) to old bareground (bare_old) and removes the whole Morph pixel portion of the dead flotant (dead_Flt). 
    
   
    use params
    implicit none
!
!    ! local variables --> none
!  
!    ! Add all new bareground to old bareground 
!    bare_old(:,2) = bare_old(:,1) + bare_new(:,1)
!
!    ! Dump previous values of old bareground and move new values to column 1 for the next processing step
!    bare_old(:,1) = bare_old(:,2)
!    bare_old(:,2) = 0
!
!    ! Remove the whole Morph pixel portion of the dead flotant 
!    dead_Flt(:,2) = dead_Flt(:,1) - (dead_Flt(:,1)/((1/256)*(1/256)))   ! 256 is the number of Morph 30 m pixels in each 480 m grid cell 
!
!    ! Dump previous values of dead flotant and move new values to column 1 for the next processing step
!    dead_Flt(:,1) = dead_Flt(:,2)
!    dead_Flt(:,2) = 0


end