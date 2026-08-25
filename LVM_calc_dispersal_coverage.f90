subroutine calc_dispersal_coverage
    ! global arrays updated by subroutine:
    !   disp_cov
    
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
    ! D_i = total coverage of that vegetation in those cells / the area of those cells (remember those cells may not be the same size)
    
    use params
    implicit none


    ! local variables
    integer :: ig                                                               ! iterator over vegetation grid 
    integer :: ic                                                               ! iterator over coverage type (columns of coverages)
    integer :: il                                                               ! iterator over flotant species within the thin and thick mat categories
    integer :: cover_group                                                      ! cover group value;  e.g., cover_group = 13 is saline emergent wetland vegetation
    real(sp) :: numerator                                                       ! numerator in the dispersal calc (sum of the area of that coverage in the surrounding cells)
    real(sp) :: denominator                                                     ! denominator in the dispersal calc (sum of the total area in the surrounding cells)
    integer :: inb                                                              ! iterator over the list of neighbor cells (near or nearest neighbor)
    integer :: neighbor                                                         ! the grid cell ID of the current neighbor 
    integer :: dispersal_class                                                  ! the dispersal class to which the current coverage belongs

    disp_cov = 0.0                                                              ! initialize array before first use
    
    

    open(unit=666, file='veg/'//trim(adjustL(fnc_tag))//'_'//year//'_V_veg_high_disp_est_debug.csv')
    write(666,'A') 'grid,numerator,denominator,sal_av_yr,wlv_smr'    


    do ig=1,ngrid
        do ic=1,ncov
            disp_cov(ig,ic)= coverages(ig,ic)                                   ! dispersal coverage for the coverages within the central cell 

            numerator = 0
            denominator = 0            
            do inb=1,max_neighbors
                neighbor = nearest_neighbors(ig,inb)                            ! neighbor is a grid cell ID 
                if (neighbor > 0) then                                          ! if neighbor index is -9999, then it has reached the end of nearest neighbors
                    numerator = numerator + (coverages(neighbor,ic) * grid_a(neighbor))     ! units of real area
                    denominator = denominator + grid_a(neighbor)                            ! units of real area                                                                       
                endif    
            end do
            if (denominator > 0) then
                disp_cov(ig,ic) = disp_cov(ig,ic) + (numerator/denominator)     ! add to it the dispersal coverage for the coverages in the surrounding cells (nearest neighbors)
            endif

            numerator = 0
            denominator = 0
            dispersal_class = cov_disp_class(ic)          
            if (dispersal_class == 2 .or. dispersal_class == 3) then
                do inb=1,max_neighbors
                    neighbor = near_neighbors(ig,inb)                           ! neighbor is a grid cell ID 
                    if (neighbor > 0) then                                      ! if neighbor index is -9999, then it has reached the end of near neighbors
                        numerator = numerator + (coverages(neighbor,ic) * grid_a(neighbor))
                        denominator = denominator + grid_a(neighbor)
                    endif
                end do
            endif
            if (denominator > 0) then
                disp_cov(ig,ic) = disp_cov(ig,ic) + (numerator/denominator) ! add to it the dispersal coverage for the coverages in the surrounding cells (nearest neighbors)
            endif 
        end do   

        write(666,3333) ig,numerator,denominator,sal_av_yr(grid_comp(ig)),wlv_smr(grid_comp(ig))
    end do 

    close(666)


    open(unit=666, file='veg/'//trim(adjustL(fnc_tag))//'_'//year//'_V_veg_high_disp_est_exp_lkd.csv')
    write(666,2222) 'exp_lkd_total',trim(adjustL(veg_coverage_file_header))
    do ig = 1,ngrid
        write(666,4444) ig,exp_lkd(ig,:),exp_lkd_total(ig)
    end do
    close(666)

    open(unit=666, file='veg/'//trim(adjustL(fnc_tag))//'_'//year//'_V_veg_high_disp_est_est_p.csv')
    write(666,'(A)') trim(adjustL(veg_coverage_file_header))
    do ig = 1,ngrid
        write(666,5555) ig,establish_P(ig,:)
    end do
    close(666)

    open(unit=666, file='veg/'//trim(adjustL(fnc_tag))//'_'//year//'_V_veg_high_disp_est_mort_p.csv')
    write(666,'(A)') trim(adjustL(veg_coverage_file_header))
    do ig = 1,ngrid
        write(666,5555) ig,mortality_P(ig,:)
    end do
    close(666)

    open(unit=666, file='veg/'//trim(adjustL(fnc_tag))//'_'//year//'_V_veg_high_disp_est_disp_cov.csv')
    write(666,'(A)') trim(adjustL(veg_coverage_file_header))
    do ig = 1,ngrid
        write(666,5555) ig,disp_cov(ig,:)
    end do
    close(666)



2222    format(A,','A)
3333    format(I0,',',F0.4,',',F0.4,',',F0.4,',',F0.4)
4444    format(I0,',',F0.4,<ncov>(',',F0.4))
5555    format(I0,',',<ncov>(',',F0.4))

    return
end 