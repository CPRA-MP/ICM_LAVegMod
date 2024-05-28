subroutine acute_salinity_flt
    ! global arrays updated by subroutine:
    !   coverages
    
    ! global arrays used by subroutine:
    !   ngrid
    !   sal_mx_14d_yr
    !   grid_comp
    !   flt_thn_cnt
    !   flt_thk_cnt 
    !   bfi
    !   dfi
    !   flt_thn_indices
    !   flt_thk_indices

    ! This subroutine changes coverages of flotant based on acute salinity stress (the threshold level is hard coded as 5.5 ppt)
   
    use params
    implicit none

    ! local variables
    integer :: ig                       ! iterator over veg grid cells
    integer :: il                       ! iterator thin or thick mat flotant species
    real(sp) :: acute_sal_threshold     ! threshold for max 14 day salinity at which acute salinity stress is or is not applied; units are ppt

    acute_sal_threshold = 5.5           ! ppt
    do ig=1,ngrid
        if (sal_mx_14d_yr(grid_comp(ig)) >= acute_sal_threshold) then
            do il=i,flt_thn_cnt
                coverages(ig,dfi,2) = coverages(ig,dfi,2) + coverages(ig,flt_thn_indices(il),2)        ! Thin mat is added to dead flotant
                coverages(ig,flt_thn_indices(il),2) = 0.0                                              ! Zero-out thin mat
            end do
            do il=i,flt_thk_cnt
                coverages(ig,bfi,2) = coverages(ig,bfi,2) + coverages(ig,flt_thk_indices(il),2)        ! Thick mat is added to bareground flotant
                coverages(ig,flt_thk_indices(il),2) = 0.0                                              ! Zero-out thick mat            
            end do            
        end if 
    end do
end 