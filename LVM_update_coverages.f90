subroutine update_coverages
    ! global arrays updated by subroutine:
    !   coverages
    !   exp_lkd_total
    
    ! global arrays used by subroutine:
    !   exp_lkd
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
    integer :: dispersal_class
    real(sp) :: total_est_P

 
    ! Sum expansion likelihood across all non-flotant species
    exp_lkd_total = 0.0                                                             ! initialize array to zero before first used
    
    do ic=1,ncov
        if (cov_grp(ic) > 7) then                                                   ! excludes all flotant and non-veg coverages (cover groups 1-7)
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
                if (cov_grp(ic) > 7) then                                 ! Excludes all flotant types and non-veg coverages (cover groups 1-7)
                    coverages(ig,ic) = coverages(ig,ic) + ((exp_lkd(ig,ic)/exp_lkd_total(ig))*total_unoccupied_lnd(ig))
                end if
            end do
        end if
    end do


    
end