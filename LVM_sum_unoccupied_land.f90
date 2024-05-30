subroutine sum_unoccupied_lnd
    ! global arrays updated by subroutine:
    !   newly_unoccupied_lnd
    !   total_unoccupied_lnd   
    
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


    ! local variables
    integer :: ic                                                                           ! iterator over veg grid coverages (columns)
    integer :: cover_group                                                                  ! cover group value;  e.g., cover_group = 13 is saline emergent wetland vegetation

    total_unoccupied_lnd = 0.0                                                              ! initialize array to zero before first used
    newly_unoccupied_lnd = 0.0                                                              ! initialize array to zero before first used
    
    do ic=1,ncov
        cover_group = cov_grp(ic)                                                           ! Excludes flotant and non-veg coverages (cover groups 1-7)
        if (cover_group > 7) then
            newly_unoccupied_lnd = newly_unoccupied_lnd + (coverages(:,ic)*mortality_p(:,ic))
        end if
    end do
    total_unoccupied_lnd = newly_unoccupied_lnd + coverages(:,bni) + coverages(:,boi)   ! add the new barground and old bareground to the newly unoccupied land to give the total 

end 