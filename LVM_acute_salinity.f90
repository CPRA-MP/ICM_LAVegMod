subroutine acute_salinity
    ! global arrays updated by subroutine:
    !   coverages
    
    ! global arrays used by subroutine:
    !   ngrid
    !   sal_mx_14d_yr
    !   grid_comp
    !   ncov
    !   cov_grp
    !   bni
    
    ! This subroutine changes coverages based on acute salinity stress (the threshold level is hard coded as 5.5 ppt)
   
    use params
    implicit none

    ! local variables
    integer :: ig                       ! iterator over veg grid cells
    integer :: ic                       ! iterator over veg grid coverages (columns)
    integer :: coverage_group           ! cover group value;  e.g., cover_group = 13 is saline emergent wetland vegetation
    real(sp) :: acute_sal_threshold     ! threshold for max 14 day salinity at which acute salinity stress is or is not applied; units are ppt

    acute_sal_threshold = 5.5           ! ppt
    do ig=1,ngrid
        if (sal_mx_14d_yr(grid_comp(ig)) >= acute_sal_threshold) then
            do ic=1,ncov
                coverage_group=cov_grp(ic)
                if (coverage_group == 10) then
                    coverages(ig,bni) = coverages(ig,bni) + coverages(ig,ic)
                    coverages(ig,ic) = 0.0
                end if
            end do
        end if 
    end do



end 