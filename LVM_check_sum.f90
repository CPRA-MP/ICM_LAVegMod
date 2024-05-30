subroutine check_sum
    ! global arrays updated by subroutine:
    !   none

    ! This subroutine checks that the coverages in each cell sum to 1.0 +/- tolerance 
    ! The tolerance level is hard coded below 
   
    use params
    implicit none

    ! local variables
    integer :: ig                       ! iterator over veg grid cells
    real(sp) :: tol                     ! level of tolerance allowed in the sum; values +/- this value are considered in range


    tol = 0.005

    do ig=1,ngrid
        if (((1.0 - tol) > sum(coverages(ig,:))) .or. ((1.0 + tol) < sum(coverages(ig,:)))) then
           ! print an error message that it's out of bounds 
            write(*,*) '*****************WARNING************************'                       
            write(*,*) 'Check sum error: Grid Cell ID ',ig,' is ',sum(coverages(ig,:)), ', which is beyond the tolerated total sum coverage.'
            write(*,*) '************************************************'
        end if 
    end do




end 