subroutine calc_dispersal_coverage(disp_cov)
    ! global arrays updated by subroutine:
    !   
    
    ! global arrays used by subroutine:
    !   ngrid
    !   ncov
    !   coverages
    !   max_neighbors
    !   nearest_neighbors
    !   near_neighbors
    !   grid_a
    !   cov_disp_class
    

    ! This subroutine calculates the dispersal coverage for every coverage in every cell 
   
    use params
    implicit none

     ! dummy local variables populated with arrays passed into subroutine
    real(sp),dimension(ngrid,ncov), intent(inout) :: disp_cov               ! dummy variable to hold the dispersal coverage of each coverage in each veg cell; it is returned to parent subroutine

    ! local variables
    integer :: ig                               ! iterator over vegetation grid 
    integer :: ic                               ! iterator over coverage type (columns of coverages)
    integer :: il                               ! iterator over flotant species within the thin and thick mat categories
    integer :: cover_group                      ! cover group value;  e.g., cover_group = 13 is saline emergent wetland vegetation
    real(sp) :: numerator                       ! numerator in the dispersal calc (sum of the area of that coverage in the surrounding cells)
    real(sp) :: denominator                     ! denominator in the dispersal calc (sum of the total area in the surrounding cells)
    integer :: inb                              ! iterator over the list of neighbor cells (near or nearest neighbor)
    integer :: neighbor                         ! the grid cell ID of the current neighbor 
    integer :: dispersal_class                  ! the dispersal class to which the current coverage belongs

    do ig=1,ngrid
        do ic=1,ncov
            disp_cov(ig,ic)= coverages(ig,ic)                            ! dispersal coverage for the coverages within the central cell 

            numerator = 0
            denominator = 0            
            do inb=1,max_neighbors
                neighbor = nearest_neighbors(ig,inb)                    ! neighbor is a grid cell ID 
                if (neighbor < 0) then                                    ! if neighbor index is -9999, then it has reached the end of nearest neighbors
                    return
                else
                    numerator = numerator + (coverages(neighbor,ic) * grid_a(neighbor))
                    denominator = denominator + grid_a(neighbor)
                endif
                disp_cov(ig,ic) = disp_cov(ig,ic) + (numerator/denominator)       ! add to it the dispersal coverage for the coverages in the surrounding cells (nearest neighbors)
            end do

            numerator = 0
            denominator = 0
            dispersal_class = cov_disp_class(ic)          
            if (dispersal_class == 2 .or. dispersal_class == 3) then
                do inb=1,max_neighbors
                    neighbor = near_neighbors(ig,inb)                       ! neighbor is a grid cell ID 
                    if (neighbor < 0) then                                    ! if neighbor index is -9999, then it has reached the end of near neighbors
                        return
                    else
                        numerator = numerator + (coverages(neighbor,ic) * grid_a(neighbor))
                        denominator = denominator + grid_a(neighbor)
                    endif
                    disp_cov(ig,ic) = disp_cov(ig,ic) + (numerator/denominator)   ! add to it the dispersal coverage for the coverages in the surrounding cells (nearest neighbors)
                end do
            endif
        end do   
    end do 


end 