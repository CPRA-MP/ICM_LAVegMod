subroutine sum_unoccupied_flt
    ! global arrays updated by subroutine:
    !   total_unoccupied_flt
    !   newly_unoccupied_thn_flt
    !   newly_unoccupied_thk_flt
    
    ! global arrays used by subroutine:
    !   ngrid
    !   ncov
    !   cov_grp
    !   coverages
    !   bni
    !   boi
    !   mortality_p
    !   flt_thn_cnt
    !   flt_thk_cnt


    ! This subroutine sums the amount of unoccupied flotant area in each veg grid cell by type (dead thin, dead thick, bareground)
   
    use params
    implicit none

    
    integer :: il                               ! iterator over flotant species within the thin and thick mat categories
    integer :: cover_group                      ! cover group value;  e.g., cover_group = 13 is saline emergent wetland vegetation

    total_unoccupied_flt = 0.0                                                      ! initialize array to zero before first used
    newly_unoccupied_thn_flt = 0.0                                                  ! initialize array to zero before first used
    newly_unoccupied_thk_flt = 0.0                                                  ! initialize array to zero before first used
    
    do il=1,flt_thn_cnt
        newly_unoccupied_thn_flt = newly_unoccupied_thn_flt + (coverages(:,flt_thn_indices(il))*mortality_p(:,flt_thn_indices(il)))
        coverages(:,flt_thn_indices(il)) = coverages(:,flt_thn_indices(il)) * (1 - mortality_p(:,flt_thn_indices(il)))
    end do
    do il=1,flt_thk_cnt
        newly_unoccupied_thk_flt = newly_unoccupied_thk_flt + (coverages(:,flt_thk_indices(il))*mortality_p(:,flt_thk_indices(il)))
        coverages(:,flt_thk_indices(il)) = coverages(:,flt_thk_indices(il)) * (1 - mortality_p(:,flt_thk_indices(il)))
    end do
    
    coverages(:,bfi) = coverages(:,bfi) + newly_unoccupied_thn_flt + newly_unoccupied_thk_flt   ! Add the newly unoccupied thin and thick mat flotant to the bareground flotant coverage (then it can be seen in the intermediate output)
    total_unoccupied_flt = coverages(:,bfi)                                                     ! total unoccupied flotant is the same now as barground flotant


end 