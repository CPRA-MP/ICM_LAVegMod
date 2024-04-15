subroutine neighbors

    use params
    implicit none

    ! PARAMS variables used:
    !   grid_x
    !   grid_y
    !   build_neighbors
    !   nearest_neighbors
    !   nearest_neighbors_file
    !   nearest_neighbors_dist
    !   near neighbors
    !   near_neighbors_file
    !   near_neighbors_dist

    ! local variables
    integer :: g0                           ! iterator to store the grid cell of interest
    integer :: gi                           ! iterator to store all other grid cells (not g0)
    integer :: maxcount                     ! count of the largest number of neighbors for any grid cell (used in output file formatting)
    integer :: count                        ! count of the number of neighbors for current grid cell
    integer :: closest_index                ! local pointer to index of minimum distance value
    integer :: current_dist                 ! distance from g0 to currrent gi
    integer :: arbitrary_max                ! initial nonzero distance that is larger than largest neighbor search radius
    real(sp),dimension(:),allocatable :: d  ! array of distances from each grid cell of interest (g0) to all other grid cells (gi)
    
    allocate(d(ngrid))
    
    arbitrary_max = near_neighbors_dist**2
    
    if (build_neighbors == 1) then
        write(  *,*) ' - calculating centroid-to-centroid distances for all grid cells and finding neighbors'
        write(000,*) ' - calculating centroid-to-centroid distances for all grid cells and finding neighbors'

        maxcount = 0
        
        ! determine which grid cells are nearest and near neighbors to each grid cell of interest (g0)
        do g0 = 1,ngrid
            d = arbitrary_max
            
            ! check to see if maximum number of neighboring cells has already been found - if so, report out error message but continue on
            if (maxcount == max_neighbors) then
                write(  *,*) ' - found the maximum number of neighbors (set in input_params). Stopping neighbor analysis.'
                write(000,*) ' - found the maximum number of neighbors (set in input_params). Stopping neighbor analysis.'
                
                open(unit=999, file=trim(adjustL('veg/__NEIGHBORING_GRID_CELL_ERRORS__.txt')))
                write(999,*) 'Found the maximum number of neighbors (set in input_params). Stopping neighbor analysis. Run continued.'
                write(999,*) 'Maximum number of neighbors was set to: ', max_neighbors
                close(999)
                
                exit
            end if
            
            
            ! calculate distance from current grid cell of interest (g0) to all other grid cells (gi)
            do gi = 1,ngrid
                d(gi) = ( ( grid_x(g0) - grid_x(gi) )**2 + ( ( grid_y(g0) - grid_y(gi) )**2 )**0.5
            end do
            
            ! determine list of all NEAREST NEIGHBOR grid cells to current grid cell of interest (g0)
            count = 1
            do gi = 1,ngrid
                closest_index = MINLOC( d,DIM=1,MASK=(d<=nearest_neighbors_dist) )
                current_dist = d(closest_index)
                if (current_dist <= nearest_neighbors_dist) then
                    if (gi /= g0) then                  ! do not include g0 in the list of neighbors
                        nearest_neighbors(gi,count) = closest_index
                        count = count + 1
                    end if
                    d(closest_index) = arbitrary_max    ! replace current minimum distance with max so it is no longer a candidate for MINLOC
                else    
                    exit                                ! this ELSE will be triggered when the minimum distance returned from MINLOC is greater than nearest_neighbors_dist
                end if
            end do
            maxcount = MAX(count,maxcount)
            
            ! determine list of all NEAR NEIGHBOR grid cells to current grid cell of interest (g0)
            count = 1
            do gi = 1,ngrid
                closest_index = MINLOC( d,DIM=1,MASK=(d<=near_neighbors_dist) )
                current_dist = d(closest_index)
                if (current_dist <= near_neighbors_dist) then
                    if (gi /= g0) then                  ! do not include g0 in the list of neighbors
                        near_neighbors(gi,count) = closest_index
                    end if
                    count = count + 1
                    d(closest_index) = arbitrary_max    ! replace current minimum distance with max so it is no longer a candidate for MINLOC
                else    
                    exit                                ! this ELSE will be triggered when the minimum distance returned from MINLOC is greater than near_neighbors_dist
                end if
            end do
            maxcount = MAX(count,maxcount)
        end do
                    

        
    else
        write(  *,*) ' - reading list of neighboring cells from file'
        write(000,*) ' - reading list of neighboring cells from file'
        
        
        
    
    end if


 

2345    format(I0,54(',',F0.4))

    return
end
