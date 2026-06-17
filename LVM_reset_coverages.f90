subroutine reset_coverages
    ! global arrays updated by subroutine:
    !      coverages
    !       boi
    !       bni
    !       dfi
    !       dem_pixel_proportion
    
    ! This subroutine converts the previous model year's new bareground to old bareground and removes the whole Morph pixel portion of the dead flotant. 
    ! It should be the first thing run (before any land/water updates from Morph)

    use params
    implicit none
    real(sp) :: residual_dead_flt_post_morph    ! portion of previous year's dead flotant coverage that was unable to be converted to water in ICM-Morph due to grid resolution 

    ! local variables
    ! none


    ! Add all new bareground to old bareground and reset the new bareground
    coverages(:,boi) = coverages(:,boi) + coverages(:,bni)
    coverages(:,bni) = 0
    
    ! Remove the whole Morph pixel portion of the dead flotant
    ! The current year updated water area, `water_from_morph`, is read in from the Morph output file.
    ! This water area will have been updated to convert the previous year's dead flotant to water, but only whole (30-m) pixels of dead flotant were able to be converted to water.
    ! This step reduces previous year's dead flotant by the portion that Morph was able to convert to water.
    ! The residual dead_flotant will remain as the starting amount of dead_flotant for the current model year.
    residual_dead_flt_post_morph = dem_pixel_proportion(:)*floor(coverages(:,dfi)/dem_pixel_proportion(:))
    coverages(:,dfi) = max(0.0,coverages(:,dfi) - residual_dead_flt_post_morph )

end