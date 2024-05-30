subroutine sum_unoccupied_lnd(total_unoccupied_lnd)
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

    ! This subroutine sums the amount of unoccupied land in each veg grid cell
   
    use params
    implicit none

     ! dummy local variables populated with arrays passed into subroutine
     real(sp),dimension(ngrid), intent(inout) :: total_unoccupied_lnd                       ! dummy variable to hold the total amount of unoccupied land in each veg grid cell; it is returned to parent subroutine


    ! local variables
    integer :: ic                                                                           ! iterator over veg grid coverages (columns)
    integer :: cover_group                                                                  ! cover group value;  e.g., cover_group = 13 is saline emergent wetland vegetation
    real(sp) :: newly_unoccupied_lnd(ngrid)                                                 ! portion of each grid cell that is unoccupied due to veg-mortality


    do ic=1,ncov
        cover_group = cov_grp(ic)                                                           ! Excludes flotant and non-veg coverages (cover groups 1-7)
        if (cover_group > 7) then
            newly_unoccupied_lnd = newly_unoccupied_lnd + (coverages(:,ic,2)*mortality_p(:,ic))
        end if
    end do
    total_unoccupied_lnd = newly_unoccupied_lnd + coverages(:,bni,2) + coverages(:,boi,2)   ! add the new barground and old bareground to the newly unoccupied land to give the total 

end 