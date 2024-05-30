subroutine update_coverages
    ! global arrays updated by subroutine:
    !   coverages
    !   exp_lkd
    !   exp_lkd_total
    
    ! global arrays used by subroutine:
    !   ngrid
    !   mortality_p
    !   establish_P
    !   cov_grp
    !   cov_disp_class
    !   bni
    !   boi

    ! This subroutine updates the coverage values of all the coverages in every grid cell 
   
    use params
    implicit none

    ! local variables 
    integer :: ig                                                                   ! iterator over Veg grid cells
    integer :: ic                                                                   ! iterator over coverage types (columns)
    integer :: cover_group
    integer :: dispersal_class
    real(sp) :: total_est_P

    
    
    ! Calculate the establishment and mortality probabilities for every species and grid cell 
    call mort_est_prob

    ! Allow class three dispersal species ("weedy") to establish on any new bareground
    call weedy_establishment

    ! Sum the total unoccupied land
    call sum_unoccupied_lnd
      
    ! Sum the unoccupied flotant, keeping track of the different types (dead thin, dead thick, and bareground flotant)
    
    call sum_unoccupied_flt

    ! Apply the mortality probabilty to the coverages
    coverages(:,:) = coverages(:,:) * (1 - mortality_p) 

    ! Calculate the dispersal coverage for each species ! D_i = total coverage of that vegetation in those cells / the area of those cells (remember those cells may not be the same size)
    call calc_dispersal_coverage

    ! Calculate expansion likelihood 
    exp_lkd = 0.0                                                                   ! initialize array to zero before first used
    exp_lkd = establish_P * disp_cov                                                ! includes all species with a establishment probabilty (e.g., flotant)

    ! Sum expansion likelihood across all non-flotant species
    exp_lkd_total = 0.0                                                             ! initialize array to zero before first used
    
    do ic=1,ncov
        cover_group = cov_grp(ic)
        if (cover_group > 7) then                                                   ! excludes all flotant and non-veg coverages (cover groups 1-7)
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
                        coverages(ig,ic) = coverages(ig,ic) +  ((establish_P(ig,ic)/total_est_P)*total_unoccupied_lnd(ig))         ! sum the estalishment Ps of the weedy species 
                    end if
                end do          
                coverages(ig,bni) = 0.0                                   ! reset the new bareground 
                coverages(ig,boi) = 0.0                                   ! reset the old bareground                 
            else
                ! leave bareground intact 
            end if
        else
            coverages(ig,bni) = 0.0                                       ! reset the new bareground 
            coverages(ig,boi) = 0.0                                       ! reset the old bareground 
            do ic=1,ncov
            cover_group = cov_grp(ic)                                     ! Identify which coverage group this coverage (column) belongs to
                if (cover_group > 7) then                                 ! Excludes all flotant types and non-veg coverages (cover groups 1-7)
                    coverages(ig,ic) = coverages(ig,ic) + ((exp_lkd(ig,ic)/exp_lkd_total(ig))*total_unoccupied_lnd(ig))
                end if
            end do
        end if
    end do

    ! Apply flotant changes 
    call update_flotant

    ! Apply acute salinity coverage changes
    call acute_salinity

    ! Apply aculte salinity covage changes for flotant 
    call acute_salinity_flt

    ! Decide what this function will do and describe it here
    call round_coverages

    ! Check that the sum of all coverages in each cell is 1.0 +/- given tolerance 
    call check_sum

    ! Check sum needs to be called before adding the whole pixel portion of dead flotant to water so it is not double counted in dead flotant and water

    ! Add the whole Morph pixel portion of the dead flotant to water
    coverages(:,wti) = coverages(:,wti) + (dem_pixel_proportion*floor(coverages(:,dfi)/dem_pixel_proportion))

    
end