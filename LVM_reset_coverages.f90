subroutine reset_coverages
    ! global arrays updated by subroutine:
    !      coverages(ngrid,ncov,2)
    !
    ! This subroutine converts the previous model year's new bareground to old bareground and removes the whole Morph pixel portion of the dead flotant. 
    ! It should be the first thing run (before any land/water updates from Morph)

    use params
    implicit none
!
!    ! local variables

!    ! Add all new bareground to old bareground 
    coverages(:,boi,2) = coverages(:,boi,1) + coverages(:,bni,1)
!
!    ! Dump previous values of old bareground and shift new values to state 1 for the next processing step
    coverages(:,boi,1) = coverages(:,boi,2)
    coverages(:,boi,2) = 0
    
!    ! Remove the whole Morph pixel portion of the dead flotant
    coverages(:,dfi,2) = coverages(:,dfi,1) - floor(coverages(:,dfi,1)/dem_pixel_proportion(:))
!    dead_Flt(:,2) = dead_Flt(:,1) - (dead_Flt(:,1)/((1/256)*(1/256)))   ! 256 is the number of Morph 30 m pixels in each 480 m grid cell 
!
!    ! Dump previous values of dead flotant and shift new values to state 1 for the next processing step
    coverages(:,dfi,1) = coverages(:,dfi,2)
    coverages(:,dfi,2) = 0 

end