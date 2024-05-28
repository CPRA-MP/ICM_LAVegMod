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
    !   wlv_smr
    !   est_X_bins
    !   mort_X_bins
    !   n_X_bins
    !   barrier_island
    !   tree_establishment

    ! This subroutine calculates the establishment and mortality probability for each species/coverage in each grid cell
   
    use params
    implicit none

    ! local variables
    integer :: ig                       ! iterator over veg grid cells
    integer :: ic                       ! iterator over veg grid coverages (columns)
    integer :: cover_group              ! cover group value;  e.g., cover_group = 13 is saline emergent wetland vegetation
    real(sp) :: tol                     ! level of tolerance allowed in the sum; values +/- this value are considered in range

    establish_P = 0                                     ! initialize with all 0s; for coverages that do not calculate an establishment probability, the value will remain 0 (nothing will establish)
    mortality_P = 0                                     ! initialize with all 0s; for coverages that do not calculate an mortality probability, the value will remain 0 (nothing will die)

    do ig = 1,ngrid                                     ! Loop through every grid cell 
        do ic = 1, ncov                                 ! Loop through every coverage (column)
            cover_group = cov_grp(ic)                   ! Identify which coverage group this coverage (column) belongs to
            if (cover_group == 8 .or. cover_group == 14) then                             ! For bottomland hardwood forest and barrier island species (coverage group 8 and 14), calculate establishment probability from elevation 
                call oneway_interp(grid_elev(ig), establish_tables(:,:,ic), est_Y_bins(:,ic), n_Y_bins, establish_P(ig,ic))    ! oneway_interp(variable1,table,variable1bins, var1bin_n, yint)
                call oneway_interp(grid_elev(ig), mortality_tables(:,:,ic), mort_Y_bins(:,ic), n_Y_bins, mortality_P(ig,ic))   ! oneway_interp(variable1,table,variable1bins, var1bin_n, yint)
            elseif (cover_group == 4 .or. cover_group == 5 .or. cover_group >= 9) then    ! For swamp forest, thick and thin floating marsh, emergent wetland (fresh, intermediate, brackish, and saline) (coverage groups 4-5, 9-13), calculate establishment probability from wlv and annual salinity
                call twoway_interp(sal_av_yr(grid_comp(ig)), wlv_smr(grid_comp(ig)), establish_tables(:,:,ic), est_Y_bins(:,ic), n_Y_bins, est_X_bins(:,ic), n_X_bins, establish_P(ig,ic))             ! twoway_interp(variable1, variable2, table, variable1bins, var1bin_n, variable2bins, var2bin_n, yint)
                call twoway_interp(sal_av_yr(grid_comp(ig)), wlv_smr(grid_comp(ig)), mortality_tables(:,:,ic), mort_Y_bins(:,ic), n_Y_bins, mort_X_bins(:,ic), n_X_bins, mortality_P(ig,ic))            ! twoway_interp(variable1, variable2, table, variable1bins, var1bin_n, variable2bins, var2bin_n, yint)
            else                                                                           ! For water, not mod, new bareground, old bareground, bareground flotant, dead flotant (coverage groups 0-3, 6-7), do nothing
                return ! is that right? 
            endif
        end do 
    end do

    ! Zero-out establish_P for barrier island species not in barrier island cells and keep it in barrier island cells
    do ic=1,ncov
        if (cov_grp == 14) then                          ! if it's a barrier island species (cover group 14), then the expansion liklihood stays in the barrier island cells (multiplied by 1) and is removed from non-barrier island cells (multiplied by 0)
            establish_P(:,ic) = establish_P(:,ic) * barrier_land ! barrier island is a 1D array of size ngrid (1 if island; 0 if not)
        end if
    end do


    ! Zero-out establish_P for non-barrier island species in barrier island cells 
    do ig=1,ngrid
        if (barrier_island(ig) > 0) then
            do ic=1,ncov
                cover_group = cov_grp(ic)                   ! Identify which coverage group this coverage (column) belongs to
                if (cover_group == 14 .or. cover_group < 4) then
                    ! do nothing, leave as is
                else
                    establish_P(ig,ic) = 0.0
                end if
            end do
        end if
    end do


    ! Zero-out establish_P for swampforest model without the tree establishment condition 
    do ic=1,ncov
        cover_group = cov_grp(ic)                                     ! Identify which coverage group this coverage (column) belongs to
        if (cover_group == 9) then                                    ! Cover group 9 is swamp forest
            establish_P(:,ic) = establish_P(:,ic) * tree_establishment        ! tree establishment is a 1D array of size ngrid (1 if conditions met; 0 if not)
        endif
    end do




end 