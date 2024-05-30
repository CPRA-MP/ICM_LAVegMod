subroutine land_change
    ! global arrays updated by subroutine:
    !       water_from_morph
    !       coverages
    !
    ! This subroutine modifies the land/waver coverage in Veg cells based on processes in Morph.  
    ! The change is based on the % of the cell that is water, as determined by Morph. 
    ! Possible outcomes are: 1) no change; 2) land gain; 3) land loss 
    ! If there is land gain, then the % water in the cell is decreased and the % new bareground is increased 
    ! If there is land loss, then the % water in the cell is increased and the % land is decreased in the following order: 
    ! old bareground is reduced first.

    use params
    implicit none
    
    ! local variables                                                                        
    integer :: ig                                                                                   ! iterator over Veg grid cells
    integer :: ic                                                                                   ! iterator over coverage types
    integer :: fi                                                                                   ! iterator of flotant coverage indices
    integer :: fli                                                                                  ! local varaible to store index of flotants in coverages(ngrid,ncov)
    real :: delta_water                                                                             ! difference between Morph water and Veg water
    real :: total_flotant                                                                           ! sum of all flotant coverage types
    real :: veg_land                                                                                ! total land (decimal %) in the Veg cell
    real :: morph_land                                                                              ! total land (decimal %) in the Morph cell
    real :: scale_land                                                                              ! ratio of Morph land to Veg land (0.0 <= scale_land <= 1.0)
    real :: no_change_threshold                                                                     ! ratio of DEM pixel to Veg grid cell indicating no change at DEM-pixel level was detected

    
    
    water_from_morph = water_from_morph * (1 - coverages(:,nmi))                                  ! Adjust the Morph water to be comparable to the Veg water (Morph treats NOTMOD as NoData, so the %water does not account for NOTMOD. Morph Water + Morph Land + Veg NOTMOD = 100%)
   
    do ig = 1, ngrid                                                                                ! Loop through every grid cell comparing the amount of water and making the needed changes
        no_change_threshold = dem_pixel_proportion(ig)                                              ! portion of LAVegMod grid cell that is one DEM pixel - can't have land change less than one pixel
        delta_water = coverages(ig,wti) - water_from_morph(ig)                                    ! Compare Morph water to Veg water and proceed accordingly
                                                                                                    ! Check if change in water area is greater *in magnitude* than the no_change_threshold; the first IF is met when land gain occurs, the ELSE IF is land loss
                                                                                                    ! 1) NO LAND CHANGE
        if (delta_water > no_change_threshold) then                                                 ! 2) LAND GAIN; delta_water is positive and greater than change threshold
            coverages(ig,bni) = coverages(ig,bni) + delta_water                                 !   - add delta_water to new bareground and subtract water 
            coverages(ig,wti) = coverages(ig,wti) - delta_water                                 !   - subtract delta_water from water
        
        else if (delta_water < -1.0*no_change_threshold) then                                       ! 3) LAND LOSS;  delta_water is negative and less than change_threshold
            delta_water = -1.0*delta_water                                                          ! redefining delta water so we're adding a positive rather than having to subtract a negative    
            ! total_flotant = ELBA2_Flt(ig) + PAHE2_Flt(ig) + bare_Flt(ig) + dead_Flt(ig)
            total_flotant = coverages(ig,dfi)                                                     ! add DEAD_Flt to count of total flotant coverage
            total_flotant = total_flotant + coverages(ig,bfi)                                     ! add BARE_Flt to count of total flotant coverage
            do fi = 1,flt_thn_cnt                                                                   ! loop over living thin mat flotant coverages and add to total flotant coverage
                fli = flt_thn_indices(fi)
                total_flotant = total_flotant + coverages(ig,fli)
            end do
            do fi = 1,flt_thk_cnt                                                                   ! loop over living thick mat flotant coverages and add to total flotant coverage
                fli = flt_thk_indices(fi)
                total_flotant = total_flotant + coverages(ig,fli)
            end do
           
            morph_land = 1.0 - (water_from_morph(ig) + coverages(ig,nmi) + total_flotant)         ! calculate land area from last ICM-Morph outputs (ignoring NotMod)
            veg_land = 1.0 - (coverages(ig,wti)+ coverages(ig,bni) + total_flotant)             ! calculate land area from last ICM-LAVegMod outputs
            
            if (morph_land < 0.0) then
                write(*,*) '*****************WARNING************************'                       ! print a warning message -- why is morph land less than 0?
                write(*,*) 'Grid Cell ID ',ig,' has negative land area.'
                write(*,*) '      morph_land: ', morph_land
                write(*,*) 'water_from_morph: ', water_from_morph(ig)
                write(*,*) '   notmod/upland: ', coverages(ig,nmi)
                write(*,*) '   total_flotant: ', total_flotant
                do fi = 1,flt_thn_cnt
                    fli = flt_thn_indices(fi)
                    write(*,*) '       ',cov_symbol(fli),': ', coverages(ig,fli)
                end do
                do fi = 1,flt_thk_cnt
                    fli = flt_thk_indices(fi)
                    write(*,*) '       ',cov_symbol(fli),': ', coverages(ig,fli)
                end do
                write(*,*) '        dead_Flt: ', coverages(ig,dfi)
                write(*,*) '************************************************'
            
            else if (veg_land > 0) then                                                             ! there is land available to be lost
                if (coverages(ig,boi) >= delta_water) then                                        !   old bareground is large enough to cover all of the land loss
                    coverages(ig,boi) = coverages(ig,boi) - delta_water                         !   reduce old bareground area by amount of new water area needed
                    coverages(ig,wti) = coverages(ig,wti) + delta_water                         !   update water area
                    scale_land = 1.0                                                                !   vegetated land will not change due to old bareground available
                
                else if (coverages(ig,boi) > 0.0) then                                            ! there is old bareground present, but it's not enough to cover the full land loss 
                    delta_water = delta_water - coverages(ig,boi)                                 !   determine residual amount of new water area needed after old bareground is exhausted
                    coverages(ig,wti) = coverages(ig,wti) + coverages(ig,boi)                 !   add old bareground area to water
                    veg_land = veg_land - coverages(ig,boi)                                       !   calculate the residual vegetated land to be lost after all of old bareground was exhausted and converted to water
                    coverages(ig,boi) = 0.0                                                       !   exhaust all of the old bareground area
                    if (veg_land > 0) then                                                          !   there is still vegetated land to be lost (reduced)
                        coverages(ig,wti) = coverages(ig,wti) + delta_water                     !       add residual amount of new water area to water
                        scale_land = morph_land/veg_land                                            !       calculate the factor by which to reduce all veg coverages (decimal less than 1.0)
                    else                                                                            !   there is not enough vegetated land available to be lost (meaning veg_land is 0)
                        write(*,*) '*****************WARNING************************'
                        write(*,*) 'Grid Cell ID ',ig,' does not have enough land area available to meet calculated land loss area.'
                        write(*,*) '************************************************'
                    end if 
                else                                                                                ! there's no old bareground to reduce first 
                    coverages(ig,wti) = coverages(ig,wti) + delta_water                         !   add to new water area to water
                    scale_land = morph_land/veg_land                                                !   calculate the factor by which to reduce all veg coverages (decimal less than 1.0)
                end if  
                
                do ic=1,ncov
                    if (cov_grp(ic) > 7) then                                                       ! skip coverage groups for water, bareground, not modeled, and flotant
                        coverages(ig,ic) = coverages(ig,ic) * scale_land
                    end if
                end do
            else 
                
                write(*,*) '*****************WARNING************************'                       ! print error message --  veg land is 0, meaning there's nothing to be reduced
                write(*,*) 'Grid Cell ID ',ig,' has no land area availalb and therefore does not have enough land area available to meet calculated land loss area.'
                write(*,*) '************************************************'
            end if  
        end if
    end do
    return
end