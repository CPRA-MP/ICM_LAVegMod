subroutine reset_coverages
    ! global arrays updated by subroutine:
    !      coverages(ngrid,ncov,2)
    !
    ! This subroutine converts the previous model year's new bareground to old bareground and removes the whole Morph pixel portion of the dead flotant. 
    ! It should be the first thing run (before any land/water updates from Morph)

    use params
    implicit none

    ! local variables

    ! Add all new bareground to old bareground and reset the new bareground
    coverages(:,boi) = coverages(:,boi) + coverages(:,bni)
    coverages(:,boi) = 0
    
    ! Remove the whole Morph pixel portion of the dead flotant
    coverages(:,dfi) = coverages(:,dfi) - dem_pixel_proportion(:)*floor(coverages(:,dfi)/dem_pixel_proportion(:))

end