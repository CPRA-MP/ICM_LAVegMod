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
    integer :: ig                                                                                                                       ! iterator over veg grid cells
    integer :: ic                                                                                                                       ! iterator over veg grid coverages (columns)
    integer :: cover_group                                                                                                              ! cover group value;  e.g., cover_group = 13 is saline emergent wetland vegetation
    real(sp) :: tol                                                                                                                     ! level of tolerance allowed in the sum; values +/- this value are considered in range
    real(sp) :: minY                                                                                                                    ! minimum value included in the establishment/mortality input tables
    real(sp) :: maxY                                                                                                                    ! maximum value included in the establishment/mortality input tables
    real(sp) :: var1                                                                                                                    ! filtered value of varialbe1 passed into oneway_interp to boound the variable by the extreme min/max values in the establishment/mortality input tables

    
    establish_P = 0                                                                                                                     ! initialize with all 0s; for coverages that do not calculate an establishment probability, the value will remain 0 (nothing will establish)
    mortality_P = 0                                                                                                                     ! initialize with all 0s; for coverages that do not calculate a mortality probability, the value will remain 0 (nothing will die)

    do ig = 1,ngrid                                                                                                                     ! Loop through every grid cell
        if (grid_comp(ig) > 0) then                                                                                                     ! check that grid cell has an allowable ICM-Hydro compartment ID
            if (grid_comp(ig)<= ncomp) then                                                                                             ! check that grid cell has an allowable ICM-Hydro compartment ID
                do ic = 1, ncov                                                                                                         ! Loop through every coverage (column)
                    cover_group = cov_grp(ic)                                                                                           ! Identify which coverage group this coverage (column) belongs to
                    if (cover_group == 8 .or. cover_group == 14) then                                                                   ! For bottomland hardwood forest and barrier island species (coverage group 8 and 14), calculate establishment probability from elevation 
                                                                                                                                        !   - oneway_interp(variable1,table,variable1bins, var1bin_n, yint)
                        minY = min(est_Y_bins(:,ic))
                        maxY = max(est_Y_bins(:,ic))
                        var1 = max(min(grid_elev(ig),maxY),minY)                                                                        ! apply low/high pass filter to limit variable1 to be set to extreme values located in the input table
                        call oneway_interp(var1, establish_tables(:,:,ic), est_Y_bins(:,ic), n_Y_bins, establish_P(ig,ic))

                        minY = min(mort_Y_bins(:,ic))
                        maxY = max(mort_Y_bins(:,ic))
                        var1 = max(min(grid_elev(ig),maxY),minY)                                                                        ! apply low/high pass filter to limit variable1 to be set to extreme values located in the input table
                        call oneway_interp(var1, mortality_tables(:,:,ic), mort_Y_bins(:,ic), n_Y_bins, mortality_P(ig,ic))
                        
                    elseif (cover_group == 4 .or. cover_group == 5 .or. cover_group >= 9) then                                          ! For swamp forest, thick and thin floating marsh, emergent wetland (fresh, intermediate, brackish, and saline) (coverage groups 4-5, 9-13), calculate establishment probability from wlv and annual salinity
                                                                                                                                        !   - twoway_interp(variable1, variable2, table, variable1bins, var1bin_n, variable2bins, var2bin_n, yint)
                        call twoway_interp(sal_av_yr(grid_comp(ig)), wlv_smr(grid_comp(ig)), mortality_tables(:,:,ic), mort_Y_bins(:,ic), n_Y_bins, mort_X_bins(:,ic), n_X_bins, mortality_P(ig,ic))
                        call twoway_interp(sal_av_yr(grid_comp(ig)), wlv_smr(grid_comp(ig)), establish_tables(:,:,ic), est_Y_bins(:,ic), n_Y_bins, est_X_bins(:,ic), n_X_bins, establish_P(ig,ic))
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


    ! Zero-out establish_P for swampforest model without the tree establishment condition 
    do ic=1,ncov
        cover_group = cov_grp(ic)                                                                                                       ! Identify which coverage group this coverage (column) belongs to
        if (cover_group == 9) then                                                                                                      ! Cover group 9 is swamp forest
            establish_P(:,ic) = establish_P(:,ic) * tree_establishment                                                                  ! tree establishment is a 1D array of size ngrid (1 if conditions met; 0 if not)
        endif
    end do




end 