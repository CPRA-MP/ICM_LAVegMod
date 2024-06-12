subroutine acute_salinity_lnd
    ! global arrays updated by subroutine:
    !   coverages
    
    ! global arrays used by subroutine:
    !   ngrid
    !   sal_mx_14d_yr
    !   grid_comp
    !   ncov
    !   cov_grp
    !   bni
    
    ! This subroutine changes coverages based on acute salinity stress (the threshold level is hard coded as 5.5 ppt). It only impacts fresh emergent wetland vegetation (cov_grp = 10). 
   
    use params
    implicit none

    ! local variables
    integer :: ig                       ! iterator over veg grid cells
    integer :: ic                       ! iterator over veg grid coverages (columns)
    real(sp) :: acute_sal_threshold     ! threshold for max 14 day salinity at which acute salinity stress is or is not applied; units are ppt

    acute_sal_threshold = 5.5           ! ppt
    do ig=1,ngrid
        if (sal_mx_14d_yr(grid_comp(ig)) >= acute_sal_threshold) then
            do ic=1,ncov
                if (cov_grp(ic) == 10) then                                      ! cov_grp = 10; fresh emergent wetland vegetation
                    coverages(ig,bni) = coverages(ig,bni) + coverages(ig,ic)
                    coverages(ig,ic) = 0.0
                end if
            end do
        end if 
    end do



end 