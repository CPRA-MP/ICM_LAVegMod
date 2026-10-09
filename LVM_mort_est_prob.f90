subroutine mort_est_prob
    ! global arrays updated by subroutine:
    !   establish_P
    !   mortality_P
    
    ! global arrays used by subroutine:
    !   ngrid
    !   ncov
    !   cov_grp
    !   grid_elev
    !   establish_tables
    !   est_Y_bins
    !   n_Y_bins
    !   mortality_tables
    !   mort_Y_bins
    !   n_Y_bins
    !   sal_av_yr
    !   grid_comp
    !   wlv_yr
    !   est_X_bins
    !   mort_X_bins
    !   n_X_bins
    !   barrier_island
    !   tree_establishment

    ! This subroutine calculates the establishment and mortality probability for each species/coverage in each grid cell
   
    use params
    implicit none

    ! local variables
    integer :: ig                                                                                                                       ! iterator over veg grid cells
    integer :: ic                                                                                                                       ! iterator over veg grid coverages (columns)
    integer :: cover_group                                                                                                              ! cover group value;  e.g., cover_group = 13 is saline emergent wetland vegetation
    real(sp) :: tol                                                                                                                     ! level of tolerance allowed in the sum; values +/- this value are considered in range
    real(sp) :: minY                                                                                                                    ! minimum value included in the establishment/mortality input tables
    real(sp) :: maxY                                                                                                                    ! maximum value included in the establishment/mortality input tables
    real(sp) :: var1                                                                                                                    ! filtered value of varialbe1 passed into oneway_interp to boound the variable by the extreme min/max values in the establishment/mortality input tables
    real(sp) :: wlv_ptile                                                                                                               ! interpolated value representing the calculated WLV percentile for the given species - used in the universal mortality table
    real(sp) :: sal_ptile                                                                                                               ! interpolated value representing the calculated SAL percentile for the given species  - used in the universal mortality table
    real(sp) :: sal_ptile_est														! shifted salinity percentile to pass in to the mortality table to return the corresponding establishment conditions
    real(sp) :: mortP_for_establish											                ! temp value to store the probability of mortality which will be transformed into the probabilty of establishment
    establish_P = 0                                                                                                                     ! initialize with all 0s; for coverages that do not calculate an establishment probability, the value will remain 0 (nothing will establish)
    mortality_P = 0                                                                                                                     ! initialize with all 0s; for coverages that do not calculate a mortality probability, the value will remain 0 (nothing will die)
	
    do ig = 1,ngrid                                                                                                                     ! Loop through every grid cell
        if (grid_comp(ig) > 0) then                                                                                                     ! check that grid cell has an allowable ICM-Hydro compartment ID
            if (grid_comp(ig)<= ncomp) then                                                                                             ! check that grid cell has an allowable ICM-Hydro compartment ID
                do ic = 1, ncov                                                                                                         ! Loop through every coverage (column)
                    cover_group = cov_grp(ic)                                                                                           ! Identify which coverage group this coverage (column) belongs to

                    if (cover_group == 8 .or. cover_group == 14) then                                                                   ! For bottomland hardwood forest and barrier island species (coverage group 8 and 14), 
                                                                                                                                        !    calculate establishment probability from elevation above mean water level
                        minY = minval(est_Y_bins(:,ic))
                        maxY = maxval(est_Y_bins(:,ic))

                        if (grid_comp(ig) > 0) then
                            var1 = max(min(grid_elev(ig)-stg_av_yr(grid_comp(ig)),maxY),minY)                                           ! apply low/high pass filter to limit variable1 to be set to extreme values located in the input table
                        else
                            var1 = minY
                        endif

                        call oneway_interp(var1, establish_tables(1,:,ic), est_Y_bins(:,ic), n_Y_bins, establish_P(ig,ic))

                        minY = minval(mort_Y_bins(:,ic))
                        maxY = maxval(mort_Y_bins(:,ic))

                        if (grid_comp(ig) > 0) then
                            var1 = max(min(grid_elev(ig)-stg_av_yr(grid_comp(ig)),maxY),minY)                                           ! apply low/high pass filter to limit variable1 to be set to extreme values located in the input table
                        else
                            var1 = minY
                        endif      

                        call oneway_interp(var1, mortality_tables(1,:,ic), mort_Y_bins(:,ic), n_Y_bins, mortality_P(ig,ic))
                    
                    elseif (cover_group == 4 .or. cover_group == 5 .or. cover_group >= 9) then                                          ! For swamp forest, thick and thin floating marsh, emergent wetland (fresh, intermediate, brackish, and saline) (coverage groups 4-5, 9-13), calculate establishment probability from wlv and annual salinity  
                        
                        if (prob_table_type == 1) then                                                                                  ! Use universal probability tables based on sal and wlv percentiles:
                                                                                                                                        !     This setting utilizes one unified probability of mortality table:
                                                                                                                                        !        -  0% mortality within the inner-quartile range
                                                                                                                                        !        - linearly increasing mortality from 25th/75th percentiles, until reaching 
                                                                                                                                        !        - 80% mortality at the 5th/95th percentiles
																	!        - linearly increasing mortality from the 5th/95th percentiles, until reaching		
                                                                                                                                        !        - 97.5% mortality at the 1st/99th percentiles, and
                                                                                                                                        !        - 100% mortality at the 0th/100th percentils
                                                                                                                                        !
                                                                                                                                        !     Fresh marsh, flotant, and swamp forest are fresh tolerant and therefore have adjusted salinity niche percentile values
                                                                                                                                        !     passed into the model via 'coverage_attributes.csv'. These species salinity percentiles were adjusted so that
                                                                                                                                        !     CRMS-derived calculated 25th percentile salinity is overwritten by the CRMS-derived 5th percentile salinity.
                                                                                                                                        !     And the original 5th-percentile, 99th, and 0th percentile salinity values are set to 0.01, 0.001, and 0.0001, respectively.
                                                                                                                                        !     This adjusts the above logic such that mortatlity probability for fresh species:
                                                                                                                                        !        -  0% mortality between the 5th/75th percentiles
                                                                                                                                        !        - linearly increasing mortality from 5th/75th percentiles, until reaching 
                                                                                                                                        !        - 80% mortality at the 0.01 ppt and the 95th percentiles
                                                                                                                                        !        - linearly increasing mortality from 0.01 ppt and 95th percentile, until reaching
                                                                                                                                        !        - 97.% mortality at 0.001 ppt and 99th percentile, and
                                                                                                                                        !        - 100% mortality at 0.0001 ppt and the 100th percentile
                                                                                                                                        !        
                                                                                                                                        !     The above adjustment is also made on the WLV for species that are tolerant of still/stagnant water.
                                                                                                                                        !     
                                                                                                                                        !     The universal probability of mortality table is also the basis for establishment probabilities.
                                                                                                                                        !     In general, the probability of establishment one minus the probability of mortality (75% mortality = 25% establishment).
                                                                                                                                        !     However, since young vegetation is considered more vulnerable to non-ideal conditions than mature veg,
                                                                                                                                        !     The salinity and WLV niches (defined by the CRMS-derived percentiles for mortality conditions) are shifted down,
                                                                                                                                        !     to less optimum conditions.  
                                                                                                                                        !
                                                                                                                                        !     For example, if a species began to senenesce with: 
                                                                                                                                        !         25th percentile salinity = 0% mortality  100% establishment
                                                                                                                                        !         20th percentile salinity = 20% mortality 80% establishment
                                                                                                                                        !         15th percentile salinity = 40% mortality 60% establishment
                                                                                                                                        !         10th percentile salinity = 60% mortality 40% establishment
                                                                                                                                        !          5th percentile salinity = 80% mortality 20% establishment
                                                                                                                                        
                                                                                                                                        !     after adjusting for growing conditions with a 5% establishment shift, that species while growing, would now establish with:
                                                                                                                                        !         25th percentile salinity = 0% mortality  80% establishment
                                                                                                                                        !         20th percentile salinity = 20% mortality 60% establishment
                                                                                                                                        !         15th percentile salinity = 40% mortality 40% establishment
                                                                                                                                        !         10th percentile salinity = 60% mortality 20% establishment
                                                                                                                                        !          5th percentile salinity = 80% mortality  0% establishment

                                                                                                                                        !     The magnitude in which to shift between mortality and establishment as shown above is defined with the input parameter: 'mort_est_shift' 
                                                                                                                                        !

                            if (cover_group == 4 .or. cover_group == 5) then								! If flotant, set WLV to median so that WLV is not a factor in the probability of mortality/establishment of flotant marsh
                                wlv_ptile = 0.5
                            else                                                                                                        ! If not flotant, convert annual wlv into a percentile via interpolation
				call oneway_interp(wlv_yr(grid_comp(ig)),wlv_niche(ic),niche_ptiles,n_niche_ptiles,wlv_ptile)	! Calculate percentile for the year's WLV based on the input WLV niche percentile values for the species
                            endif
                            
                            call oneway_interp(sal_av_yr(grid_comp(ig)),sal_niche(ic),niche_ptiles,n_niche_ptiles,sal_ptile)		! Calculate percentile for the year's salinity based on the input salinity niche percentile values for the species

                            sal_ptile_est = max(0.0, min(1.0,sal_ptile + mort_est_shift) )						! Shift the salinity percentile used in mortality lookup for use in determining establishment probability
                            
                            call twoway_interp(sal_ptile, wlv_ptile, universal_mortality_table, univ_ptile_Y_bins, n_ptile, univ_ptile_X_bins, n_ptile, mortality_P(ig,ic)) 
                            call twoway_interp(sal_ptile_est, wlv_ptile, universal_mortality_table, univ_ptile_Y_bins, n_ptile, univ_ptile_X_bins, n_ptile, mortP_for_establish)
                            
                            establish_P(ig,ic) = max(0.0001,1.0 - mortP_for_establish)                                                  ! establishment = (1.0-mortality), after applying the establishment percentile shift
                                                                                                                                        ! do not let establishment equal true zero since this is about relative competition for establishment cover if no species have good establishment criteria, all species will have ability to compete equally

                        elseif (prob_table_type == 2) then                                                                              ! Use species-specific probability tables based on sal and wlv values instead of universal table logic above

                            call twoway_interp(sal_av_yr(grid_comp(ig)), wlv_yr(grid_comp(ig)), mortality_tables(:,:,ic), mort_Y_bins(:,ic), n_Y_bins, mort_X_bins(:,ic), n_X_bins, mortality_P(ig,ic))
                            call twoway_interp(sal_av_yr(grid_comp(ig)), wlv_yr(grid_comp(ig)), establish_tables(:,:,ic), est_Y_bins(:,ic), n_Y_bins, est_X_bins(:,ic), n_X_bins, establish_P(ig,ic))

			endif
                                                                                                                                        ! For water, not mod, new bareground, old bareground, bareground flotant, dead flotant (coverage groups 0-3, 6-7), do nothing
                    endif
                end do 
            end if
        end if
    end do

    ! Zero-out establish_P for barrier island species not in barrier island cells and keep it in barrier island cells
    do ic=1,ncov
        if (cov_grp(ic) == 14) then                                                                                                     ! if it's a barrier island species (cover group 14), then the expansion liklihood stays in the barrier island cells (multiplied by 1) and is removed from non-barrier island cells (multiplied by 0)
            establish_P(:,ic) = establish_P(:,ic) * barrier_island                                                                      ! barrier island is a 1D array of size ngrid (1 if island; 0 if not)
        end if
    end do


    ! Zero-out establish_P for non-barrier island species in barrier island cells 
    do ig=1,ngrid
        if (barrier_island(ig) > 0) then
            do ic=1,ncov
                cover_group = cov_grp(ic)                                                                                               ! Identify which coverage group this coverage (column) belongs to
                if (cover_group == 14 .or. cover_group < 4) then
                    ! do nothing, leave as is
                else
                    establish_P(ig,ic) = 0.0                                                                                            ! The initial map has been processed so all non-barrier island species are removed from the BI areas. However if the initial map has non-barrier island species on barrier island cells, then the current code only stops them from expanding and does not remove them
                end if
            end do
        end if
    end do


    ! Zero-out establish_P for swampforest and bottomland hardwood without the tree establishment condition 
    do ic=1,ncov
        cover_group = cov_grp(ic)                                                                                                       ! Identify which coverage group this coverage (column) belongs to
        if (cover_group == 9 .or. cover_group == 8) then                                                                                                      ! Cover group 9 is swamp forest
            establish_P(:,ic) = establish_P(:,ic) * tree_establishment                                                                  ! tree establishment is a 1D array of size ngrid (1 if conditions met; 0 if not)
        end if
    end do

    ! Zero-out establish_P for bottomland hardwood with salinity greater than 1 ppt; increase mortality_P to 100% (1.0) if salinity is greater than 1 ppt
    do ig=1,ngrid
        do ic=1,ncov
            cover_group = cov_grp(ic)                                                                                                       ! Identify which coverage group this coverage (column) belongs to
            if (cover_group == 8) then  
                if (grid_comp(ig) > 0) then                                                                                                 ! check that grid cell has an allowable ICM-Hydro compartment ID
                    if (sal_av_yr(grid_comp(ig)) > 1.0) then                                                                                ! Check if the salinity threshold was crossed. Bottomland hardwood does not have salinty in the est/mort tables, so this criteria is necessary
                        mortality_P(ig,ic) = 1.0                                                                                            ! Removes all bottomland hardwood species coverage from that cell
                        establish_P(ig,ic) = 0.0                                                                                            ! Stops any establishement of bottomland hardwood species in that cell 
                    end if
                end if                                                          
            end if
        end do
    end do



end 