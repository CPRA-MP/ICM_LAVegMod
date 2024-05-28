subroutine sum_unoccupied_flt(total_unoccupied_flt, newly_unoccupied_thn_flt, newly_unoccupied_thk_flt)
    ! global arrays updated by subroutine:
    !   
    
    ! global arrays used by subroutine:
    !   ngrid
    !   ncov
    !   cov_grp
    !   coverages
    !   bni
    !   boi
    !   mortality_p

    ! This subroutine sums the amount of unoccupied flotant area in each veg grid cell by type (dead thin, dead thick, bareground)
   
    use params
    implicit none

     ! dummy local variables populated with arrays passed into subroutine
    real(sp),dimension(ngrid), intent(inout) :: total_unoccupied_flt               ! dummy variable to hold the total amount of unoccupied land in each veg grid cell; it is returned to parent subroutine
    real(sp),dimension(ngrid), intent(inout) :: total_unoccupied_thn_flt           ! dummy variable to hold the total amount of unoccupied land in each veg grid cell; it is returned to parent subroutine
    real(sp),dimension(ngrid), intent(inout) :: total_unoccupied_thk_flt           ! dummy variable to hold the total amount of unoccupied land in each veg grid cell; it is returned to parent subroutine


    ! local variables
    integer :: il                               ! iterator over flotant species within the thin and thick mat categories
    integer :: cover_group                      ! cover group value;  e.g., cover_group = 13 is saline emergent wetland vegetation
    real(sp) :: newly_unoccupied_lnd(ngrid)     ! portion of each grid cell that is unoccupied due to veg-mortality


    do il=i,size(flt_thn_indices)
        newly_unoccupied_thn_flt = newly_unoccupied_thn_flt + (coverages(:,flt_thn_indices(il),2)*mortality_p(:,ic))
    end do
    do il=i,size(flt_thk_indices)
        newly_unoccupied_flt = newly_unoccupied_flt + (coverages(:,flt_thk_indices(il),2)*mortality_p(:,ic))
    end do
    total_unoccupied_flt = newly_unoccupied_thn_flt + newly_unoccupied_thk_flt + coverages(:,bfi,2)

end 