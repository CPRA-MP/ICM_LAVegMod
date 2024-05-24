subroutine update_coverages
    ! global arrays updated by subroutine:
    !    to be filled in
    !
    ! This subroutine 
   
    use params
    implicit none

    ! local variables 
    real :: unoccupied_land          !come back to make this an array    ! total land on which vegetation can establish
    integer :: ig                                                        ! iterator over Veg grid cells

! To-Do:
! add intermediate print outs
! figure out the 1s and 2s
! add all local variables 
! acute salinity 
! tree and barrier island things
! flotant
! check sum


    ! Calculate the establishment and mortality probabilities for every species and grid cell 
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


    ! Allow class three dispersal species ("weedy") to establish on any new bareground
    do ig=1,ngrid
        if (coverages(ig,bni,2) > 0) then
            total_est_P = 0.0
            do ic=1,ncov
                dispersal_class = cov_disp_class(ic)
                if (dispersal_class == 3) then
                    total_est_P = total_est_P + establish_P(ig,ic)          ! sum the estalishment Ps of the weedy species 
                endif
            end do

            do ic=1,ncov
                dispersal_class = cov_disp_class(ic)
                if (dispersal_class == 3) then
                    coverages(ig,ic,2) = coverages(ig,ic,1) +  ((establish_P(ig,ic)/total_est_P)*coverages(ig,bni,2))         ! sum the estalishment Ps of the weedy species 
                endif
            end do          
            coverages(ig,bni,2) = 0.0 ! reset the new bareground 
            end do 
        else 
        endif
    end do

    ! Apply the mortality probabilty 
    coverages(:,:,2) = coverages(:,:,1) * (1 - mortality_p) ! Need to come back to this -- this is where the 1 and 2s get tricky -- 2 has been updated for bareground
  
    ! Sum the total unoccupied land 
    do ic=1,ncov
        if (FFIBS(ic)>0)
            newly_unoccupied_lnd = newly_unoccupied_lnd + coverages(:,ic,2)             ! If statement excludes flotant and non-veg coverages
        end if
    end do
    total_unoccupied_lnd = newly_unoccupied_lnd + coverages(:,bni,2) + coverages(:,boi,2)   ! add the new barground and old bareground to it 
   
    ! Sum the unoccupied flotant 
    do il=i,size(flt_thn_indices)
        newly_unoccupied_flt = newly_unoccupied_flt + coverages(:,flt_thn_indices(il),2)
    end do
    do il=i,size(flt_thk_indices)
        newly_unoccupied_flt = newly_unoccupied_flt + coverages(:,flt_thk_indices(il),2)
    end do
    total_unoccupied_flt = newly_unoccupied_flt + coverages(:,bfi,2)


    ! Calculate the dispersal coverage for each species ! D_i = total coverage of that vegetation in those cells / the area of those cells (remember those cells may not be the same size)
    do ig=i,ngid
        do ic=1,ncov
            disp_cov(ig,ic)= coverages(ig,ic,2)                            ! dispersal coverage for the coverages within the central cell 

            numerator = 0
            denominator = 0            
            do inb=1,max_neighbors
                neighbor = nearest_neighbors(ig,inb)                    ! neighbor is a grid cell ID 
                if (neighbor < 0) then                                    ! if neighbor index is -9999, then it has reached the end of nearest neighbors
                    return
                elseif
                    numerator = numerator + (coverages(neighbor,ic,2) * grid_a(neighbor))
                    denominator = denominator + grid_a(neighbor)
                endif
                disp_cov(ig,ic) = disp_cov(ig,ic) + (numerator/denominator)       ! add to it the dispersal coverage for the coverages in the surrounding cells (nearest neighbors)
            end do

            numerator = 0
            denominator = 0
            dispersal_class = cov_disp_class(ic)          
            if (dispersal_class == 2 .or. dispersal_class == 3) then
                do inb=1,max_neighbors
                    neighbor = near_neighbors(ig,inb)                       ! neighbor is a grid cell ID 
                    if (neighbor < 0) then                                    ! if neighbor index is -9999, then it has reached the end of near neighbors
                        return
                    elseif
                        numerator = numerator + (coverages(neighbor,ic,2) * grid_a(neighbor))
                        denominator = denominator + grid_a(neighbor)
                    endif
                    disp_cov(ig,ic) = disp_cov(ig,ic) + (numerator/denominator)   ! add to it the dispersal coverage for the coverages in the surrounding cells (nearest neighbors)
                end do
            endif
        end do   
    end do 



    ! Calculate expansion likelihood 
    exp_lkd = establish_P * disp_cov                                               ! includes all species with a establishment probabilty (e.g., flotant)

    ! Sum expansion likelihood across all non-flotant species
    do ic=1,ncov
        if (FFIBS(ic)>0) then
            exp_lkd_total = exp_lkd_total + exp_lkd(:,ic)
        end if
    end do
    
    ! If exp_lkd_total is 0, then allow dispersal class 3 (weedy) species to establish; otherwise, allow all non-flotant veg to establish
    do ig=1,ngrid
        if (exp_lkd_total(ig) == 0.0) then
            total_est_P = 0.0
            do ic=1,ncov
                dispersal_class = cov_disp_class(ic)
                if (dispersal_class == 3) then
                    total_est_P = total_est_P + establish_P(ig,ic)          ! sum the estalishment Ps of the weedy species 
                end if
            end do
            if (total_est_P>0) then
                do ic=1,ncov
                    dispersal_class = cov_disp_class(ic)
                    if (dispersal_class == 3) then
                        coverages(ig,ic,2) = coverages(ig,ic,1) +  ((establish_P(ig,ic)/total_est_P)*total_unoccupied_lnd)         ! sum the estalishment Ps of the weedy species 
                    end if
                end do          
                coverages(ig,bni,2) = 0.0                                   ! reset the new bareground 
                coverages(ig,boi,2) = 0.0                                   ! reset the old bareground                 
            else
                ! leave bareground intact 
            end if
        else
            coverages(ig,bni,2) = 0.0                                       ! reset the new bareground 
            coverages(ig,boi,2) = 0.0                                       ! reset the old bareground 
            do ic=1,ncov
                if (FFIBS(ic)>0) then
                    coverages(ig,ic,2) = coverages(ig,ic,2) + ((exp_lkd(ig,ic)/exp_lkd_total(ig))*total_unoccupied_lnd)
                end if
            end do
        end if
    end do

    ! Sum expansion liklihood across flotant species
    do il=i,size(flt_thn_indices)
        exp_lkd_total_flt = exp_lkd_total_flt + exp_lkd(:,flt_thn_indices(il))
    end do
    do il=i,size(flt_thk_indices)
        exp_lkd_total_flt = exp_lkd_total_flt + exp_lkd(:,flt_thk_indices(il))
    end do
    
    ! Do the flotant stuff
    ! call acute salinity ? or leave it here
    ! call check_sum and round to nearest 1 m2
end