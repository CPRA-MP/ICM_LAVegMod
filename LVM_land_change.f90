subroutine land_change
    ! global arrays updated by subroutine:
    !       water_from_morph
    !       all coverages
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
    integer :: ig                                                       ! iterator over Veg grid cells
    real :: delta_water                                                 ! difference between Morph water and Veg water
    real :: total_flotant                                               ! sum of all flotant coverage types
    real :: veg_land                                                    ! total land (decimal %) in the Veg cell
    real :: morph_land                                                  ! total land (decimal %) in the Morph cell
    real :: scale_land                                                  ! ratio of Morph land to Veg land
    real :: no_change_threshold                                         ! ratio of DEM pixel to Veg grid cell indicating no change at DEM-pixel level was detected
    
    ! Adjust the Morph water to be comparable to the Veg water (Morph treats NOTMOD as NoData, so the %water does not account for NOTMOD. Morph Water + Morph Land + Veg NOTMOD = 100%)
    water_from_morph = water_from_morph * (1 - upland(:,1))             ! NOTMOD is now called upland
    no_change_threshold = dem_res**2 / grid_res**2                      ! portion of LAVegMod grid cell that is one DEM pixel !! MOVE TO PARAMS AND MOVE CALCULATION TO SET_IO IF NEEDED IN OTHER SUBROUTINES !!
    
    ! Loop through every grid cell comparing the amount of water and making the needed changes
    do ig = 1, ngrid
        delta_water = water(ig,1) - water_from_morph(ig)                ! Compare Morph water to Veg water and proceed accordingly
                                                                        ! 1) No change; 0.00390625 = (30*30)/(480*480) is the fraction of one morph cell in a veg grid cell
        if (delta_water > no_change_threshold) then                     ! 2) Land gain; add new bareground and subtract water 
            bare_new(ig,2) = bare_new(ig,1) + delta_water
            water(ig,2) = water(ig,1) - delta_water
        
        else if (delta_water < -1.0*no_change_threshold) then           ! 3) land loss; add to water, convert bareground to water, then convert vegetation to water
            delta_water = -1.0*delta_water                              ! redefining delta water so we're adding a positive rather than having to subtract a negative    
            ! calculate land in Veg cell and Morph cell
            total_flotant = ELBA2_Flt(ig,1) + PAHE2_Flt(ig,1) + bare_Flt(ig,1) + dead_Flt(ig,1)                    !
            morph_land = 1.0 - (water_from_morph(ig) + upland(ig,1) + total_flotant)
            veg_land = 1.0 - (water(ig,1) + upland(ig,1) + total_flotant)
            
            if (morph_land < 0.0) then
                ! print a warning message -- why is morph land less than 0?
                write(*,*) 'Grid Cell ID ',ig,' has negative land area.'
                write(*,*) '      morph_land: ', morph_land
                write(*,*) 'water_from_morph: ', water_from_morph(ig)
                write(*,*) '          upland: ', upland
                write(*,*) '   total_flotant: ', total_flotant
                write(*,*) '       ELBA2_Flt: ', ELBA2_Flt(ig,1)
                write(*,*) '       PAHE2_Flt: ', PAHE2_Flt(ig,1)
                write(*,*) '        bare_Flt: ', bare_Flt(ig,1)
                write(*,*) '        dead_Flt: ', dead_Flt(ig,1)
                
            else if (veg_land > 0) then                                           ! there is land to be lost
                
                if (bare_old(ig,1) >= delta_water) then                       ! old bareground is large enough to cover the change
                    bare_old(ig,2) = bare_old(ig,1) - delta_water
                    water(ig,2) = water(ig,1) + delta_water
                
                else if (bare_old(ig,1) > 0.0) then                                   ! there is old bareground present, but it's not enough to cover the full land loss 
                    delta_water = delta_water - bare_old(ig,1)
                    water(ig,2) = water(ig,1) + bare_old(ig,1)
                    veg_land = veg_land - bare_old(ig,1)
                    bare_old(ig,2) = 0.0 
                    
                    if (veg_land > 0) then                                        ! there is still vegetated land to be lost (reduced)
                        water(ig,2) = water(ig,2) + delta_water
                        ! Calculate the factor by which to reduce all veg coverages
                        scale_land = morph_land/veg_land                        ! decimal less than 1
                    else    !meaning veg_land is 0
                        ! print an error message that the full land loss could not be accounted for
                    end if 
                else    ! there's no old bareground to reduce first 
                    water(ig,2) = water(ig,1) + delta_water
                    ! Calculate the factor by which to reduce all veg coverages
                    scale_land = morph_land/veg_land                            ! decimal less than 1
                end if  
                ! Apply the scale to all veg coverages
                QULA3(ig,2) = QULA3(ig,1) * scale_land           ! QULA3
                QULE(ig,2)  = QULE(ig,1)  * scale_land           ! QULE
                QUNI(ig,2)  = QUNI(ig,1)  * scale_land           ! QUNI
                QUTE(ig,2)  = QUTE(ig,1)  * scale_land           ! QUTE
                QUVI(ig,2)  = QUVI(ig,1)  * scale_land           ! QUVI
                ULAM(ig,2)  = ULAM(ig,1)  * scale_land           ! ULAM
                NYAQ2(ig,2) = NYAQ2(ig,1) * scale_land           ! NYAQ2
                SANI(ig,2)  = SANI(ig,1)  * scale_land           ! SANI
                TADI2(ig,2) = TADI2(ig,1) * scale_land           ! TADI2
                COES(ig,2)  = COES(ig,1)  * scale_land           ! COES
                MOCE2(ig,2) = MOCE2(ig,1) * scale_land           ! MOCE2
                PAHE2(ig,2) = PAHE2(ig,1) * scale_land           ! PAHE2
                SALA2(ig,2) = SALA2(ig,1) * scale_land           ! SALA2
                ZIMI(ig,2)  = ZIMI(ig,1)  * scale_land           ! ZIMI
                CLMA10(ig,2)= CLMA10(ig,1)* scale_land           ! CLMA10
                ELCE(ig,2)  = ELCE(ig,1)  * scale_land           ! ELCE
                IVFR(ig,2)  = IVFR(ig,1)  * scale_land           ! IVFR
                PAVA(ig,2)  = PAVA(ig,1)  * scale_land           ! PAVA
                PHAU7(ig,2) = PHAU7(ig,1) * scale_land           ! PHAU7
                POPU5(ig,2) = POPU5(ig,1) * scale_land           ! POPU5
                SALA(ig,2)  = SALA(ig,1)  * scale_land           ! SALA
                SCCA11(ig,2)= SCCA11(ig,1)* scale_land           ! SCCA11
                TYDO(ig,2)  = TYDO(ig,1)  * scale_land           ! TYDO
                SCAM6(ig,2) = SCAM6(ig,1) * scale_land           ! SCAM6
                SCRO5(ig,2) = SCRO5(ig,1) * scale_land           ! SCRO5
                SPCY(ig,2)  = SPCY(ig,1)  * scale_land           ! SPCY
                SPPA(ig,2)  = SPPA(ig,1)  * scale_land           ! SPPA
                AVGE(ig,2)  = AVGE(ig,1)  * scale_land           ! AVGE
                DISP(ig,2)  = DISP(ig,1)  * scale_land           ! DISP
                JURO(ig,2)  = JURO(ig,1)  * scale_land           ! JURO
                SPAL(ig,2)  = SPAL(ig,1)  * scale_land           ! SPAL
                BAHABI(ig,2)= BAHABI(ig,1)* scale_land           ! BAHABI
                DISPBI(ig,2)= DISPBI(ig,1)* scale_land           ! DISPBI
                PAAM2(ig,2) = PAAM2(ig,1) * scale_land           ! PAAM2
                SOSE(ig,2)  = SOSE(ig,1)  * scale_land           ! SOSE
                SPPABI(ig,2)= SPPABI(ig,1)* scale_land           ! SPPABI
                SPVI3(ig,2) = SPVI3(ig,1) * scale_land           ! SPVI3
                STHE9(ig,2) = STHE9(ig,1) * scale_land           ! STHE9
                UNPA(ig,2)  = UNPA(ig,1)  * scale_land           ! UNPA                
            else 
                !print error message --  veg land is 0, meaning there's nothing to be reduced   and/or  ! print error message : print('WARNING: Cell ID = ' + str(loc) + '. Ignored a change in land loss after reducing bareground because current land is ' +  str(curLand))

            end if  
        end if
    end do
    return
end