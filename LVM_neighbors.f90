subroutine neighbors

    use params
    implicit none

    ! PARAMS variables used:
    !   build_neighbors
    !   grid_x
    !   grid_y
    !   max_neighbors
    !   nearest_neighbors
    !   nearest_neighbors_file
    !   nearest_neighbors_dist
    !   near neighbors
    !   near_neighbors_file
    !   near_neighbors_dist

    ! local variables
    integer :: g                            ! iterator
    integer :: g0                           ! iterator to store the grid cell of interest
    integer :: gi                           ! iterator to store all other grid cells (not g0)
    integer :: count                        ! count of the number of neighbors for current grid cell
    integer :: closest_index                ! local pointer to index of minimum distance value
    integer :: current_dist                 ! distance from g0 to currrent gi
    integer :: arbitrary_max                ! initial nonzero distance that is larger than largest neighbor search radius
    real(sp),dimension(:),allocatable :: d  ! array of distances from each grid cell of interest (g0) to all other grid cells (gi)
    
    allocate(d(ngrid))
    
    arbitrary_max = near_neighbors_dist**2
    
    ! initialize arrays of nearest and near neighbors to negative value so each grid cell can have a unique number of neighbors - will require >0 filter when looking up neighboring grid cell IDs
    nearest_neighbors = -9999
    near_neighbors = -9999
    
    if (build_neighbors == 1) then
        write(  *,'(A)') ' - calculating centroid-to-centroid distances for all grid cells and finding neighbors'
        write(000,'(A)') ' - calculating centroid-to-centroid distances for all grid cells and finding neighbors'

        ! determine which grid cells are nearest and near neighbors to each grid cell of interest (g0)
        do g0 = 1,ngrid
            d = arbitrary_max
            
            ! calculate distance from current grid cell of interest (g0) to all other grid cells (gi)
            do gi = 1,ngrid
                d(gi) = ( ( grid_x(g0) - grid_x(gi) )**2 + ( ( grid_y(g0) - grid_y(gi) )**2 ) )**0.5
            end do
            
            ! determine list of all NEAREST NEIGHBOR grid cells to current grid cell of interest (g0)
            count = 0
            do gi = 1,max_neighbors
                closest_index = MINLOC( d,DIM=1,MASK=(d<=nearest_neighbors_dist) )
                current_dist = d(closest_index)
                if (current_dist <= nearest_neighbors_dist) then
                    if (closest_index /= g0) then                  ! do not include g0 in the list of neighbors (count=0 when current_dist = 0 = d(closest_index)
                        nearest_neighbors(g0,count) = closest_index
                    end if
                    count = count + 1
                    d(closest_index) = arbitrary_max    ! replace current minimum distance with max so it is no longer a candidate for MINLOC
                else    
                    exit                                ! this ELSE will be triggered when the minimum distance returned from MINLOC is greater than nearest_neighbors_dist
                end if
            end do
            
            ! determine list of all NEAR NEIGHBOR grid cells to current grid cell of interest (g0)
            count = 0
            do gi = 1,max_neighbors
                closest_index = MINLOC( d,DIM=1,MASK=(d<=near_neighbors_dist) )
                current_dist = d(closest_index)
                if (current_dist <= near_neighbors_dist) then
                    if (closest_index /= g0) then                  ! do not include g0 in the list of neighbors (count=0 when current_dist = 0 = d(closest_index)
                        near_neighbors(g0,count) = closest_index
                    end if
                    count = count + 1
                    d(closest_index) = arbitrary_max    ! replace current minimum distance with max so it is no longer a candidate for MINLOC
                else    
                    exit                                ! this ELSE will be triggered when the minimum distance returned from MINLOC is greater than near_neighbors_dist
                end if
            end do
        end do
        
        ! write lists of NEAREST and NEAR NEIGHBOR grid cells to file
        write(  *,'(A,A)') ' - writing nearest neighbor file: ', trim(adjustL(nearest_neighbors_file))
        write(000,'(A,A)') ' - writing nearest neighbor file: ', trim(adjustL(nearest_neighbors_file))
        
        open(unit=202, file=trim(adjustL(nearest_neighbors_file)))
        write(202,'(A,I0,A)') 'gridID,nearest neighbors within ',nearest_neighbors_dist,' m:'
        do g0 = 1,ngrid
            write(202,2345) g0, nearest_neighbors(g0,:)
        end do
        close(202)
        
        write(  *,'(A,A)') ' - writing near neighbor file: ', trim(adjustL(near_neighbors_file))
        write(000,'(A,A)') ' - writing near neighbor file: ', trim(adjustL(near_neighbors_file))
    
        open(unit=203, file=trim(adjustL(near_neighbors_file)))
        write(202,'(A,I0,A)') 'gridID,near neighbors within ',near_neighbors_dist,' m:'
        do g0 = 1,ngrid
            write(203,2345) g0, near_neighbors(g0,:)
        end do
        close(203)
    
    else
        ! read lists of NEAREST and NEAR NEIGHBOR grid cells to file
        write(  *,'(A,A)') ' - reading list of nearest neighboring cells from file: ', trim(adjustL(nearest_neighbors_file))
        write(000,'(A,A)') ' - reading list of nearest neighboring cells from file: ', trim(adjustL(nearest_neighbors_file))
        
        open(unit=202, file=trim(adjustL(nearest_neighbors_file)))
        read(202,*) dump_txt
        do g = 1,ngrid
            read(202,*)
            read(202,2345) g0, nearest_neighbors(g0,:)   ! read in grid cell of interest from file and assign list of nearest neighbors to array
        end do
        close(202)
        
        write(  *,*) ' - reading list of near neighboring cells from file: ', trim(adjustL(near_neighbors_file))
        write(000,*) ' - reading list of near neighboring cells from file: ', trim(adjustL(near_neighbors_file))
    
        open(unit=203, file=trim(adjustL(near_neighbors_file)))
        write(202,'(A,I0,A)') 'gridID,near neighbors within ',near_neighbors_dist,' m:'
        do g = 1,ngrid
            read(203,2345) g0, near_neighbors(g0,:)   ! read in grid cell of interest from file and assign list of nearest neighbors to array
        end do
        close(203)
        
    
    end if


 
2345 format(I0,<max_neighbors>(',',I0))


    return
end
