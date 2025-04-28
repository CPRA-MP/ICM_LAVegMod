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
 
    total_unoccupied_lnd = 0.0                                                              ! initialize 1D array to zero before first used
    newly_unoccupied_lnd = 0.0                                                              ! initialize 1D array to zero before first used
    
    ! Apply the mortality probability 
    do ic=1,ncov                                  
        if (cov_grp(ic)  > 7) then                                                          ! Excludes flotant and non-veg coverages (cover groups 1-7)
            newly_unoccupied_lnd = newly_unoccupied_lnd + (coverages(:,ic)*mortality_p(:,ic))
            coverages(:,ic) = coverages(:,ic)*(1 - mortality_p(:,ic))
        end if
    end do
    coverages(:,bni) = coverages(:,bni) + newly_unoccupied_lnd                          ! add the newly unoccupied land to the new bareground category (then it can be seen in the intermediate output)
    total_unoccupied_lnd = coverages(:,bni) + coverages(:,boi)                          ! add the new barground and old bareground for the total unuccupied land in each cell 

end 