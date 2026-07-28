subroutine twoway_interp(y, x, table, Yrows, nYrows, Xcols, nXcols, VALxy)
    
    ! inputs to the subroutine:
    !   x      : the (x) value of the input variable specific to the grid cell. For MP29, X is water level variability. 
    !   y      : the (y) value of the input variable specific to the grid cell. For MP29, Y is mean annual salinity. 
    !   table  : establishment or mortality table for one species (2D) being interpolated over in the X- and Y-dimensions
    !   Xcols  : the 1d array of column header values, in ascending order, for the X-dimension of (table). For MP29, Xcols is (est_X_bins), or (mort_X_bins), for one species (1D)
    !   nXcols : the number of 'bins' used to discretizethe the interpolation table in the X-dimension
    !   Yrows  : the 1d array of row header values, in ascending order, for the Y-dimension of (table). For MP29, Yrows is (est_Y_bins), or (mort_Y_bins), for one species (1D)
    !   nYrows : the number of 'bins' used to discretizethe the interpolation table in the Y-dimension
    !
    ! output of the subroutine:
    !   VALxy  : the output interpolated value for the input x and y values
    !
    ! global arrays updated by subroutine for MP29:
    !   establish_P : this subroutine returns the output (VALxy) and can return that into a variable, or in the case of LVM_mort_est_prob.f90 into a speficic location of a pre-existing array, specifically the (mortality_P) and (establishment_P) probablity arrays for a given grid cell(ig) and coverage type (ic)
    !   mortality_P : this subroutine returns the output (VALxy) and can return that into a variable, or in the case of LVM_mort_est_prob.f90 into a speficic location of a pre-existing array, specifically the (mortality_P) and (establishment_P) probablity arrays for a given grid cell(ig) and coverage type (ic)
    !
    ! global arrays used by subroutine:
    !   refer to the parent subroutines which call this subroutine, global arrays are passed into this subroutine and locally stored. Refer to the intent(in) and intent(out) variables defined below.
    !   
    !
    ! This subroutine interpolates the establishment or mortality probability based on two inputs. For each grid cell, it is iteratively called across all species and interpolates one species' establishment/mortality table at a time. 
    ! This subroutine is called for any LAVegMod cover group that uses two variables to define the probability of mortality and establishment.
    ! For MP23 and MP29, the cover groups that use two variables are: swamp forest, thick and thin floating marsh, and emergent wetland (fresh, intermediate, brackish, and saline). 
    ! For MP23 and MP29, the est and mort probabilities of those species are a function of mean annual salinity (sal_av_yr) and water level variabilty (wlv_smr). 
    ! 
    ! The logic/pseudocode for this 2-dimensional interpolation is based on the structure of the MP29 establishment and mortality tables that are ordered:
    !    - in ascending order in the X-dimension from left to right, and
    !    - in ascending order in the Y-dimension from top to bottom.
    !
    ! In the schematic below, to find the interpolated value at location (x,y) the boundaries surrounding the desired location at (x,y) are defined as:
    !    - the upper left bound is to the (l)eft of (x), and (a)bove (y)
    !    - the upper right bound is to the (r)ight of (x), and (a)bove (y)
    !    - the lower left bound is to the (l)eft of (x), and (b)elow (y)
    !    - the lower right bound is to the (r)ight of (x), and (b)elow (y)
    !               
    !               
    !                              left       x         right
    !                         (-) ----------------------------------------(+)
    !                above     | VAL(l,a)   VAL(x,a)   VAL(r,a)     _ _           
    !                          |                                     |   (y-above)
    !                  y       | VAL(l,y)  *VAL(x,y)*  VAL(r,y)     _|_              
    !                          |           
    !                below     | VAL(l,b)   VAL(x,b)   VAL(r,b)           
    !                          |
    !                          |  |--(x-left)--|
    !                          |
    !                         (+)
    !               
    !               
    !               

              
               
    use params
    implicit none

    real(sp),intent(in) :: x                                     !   value of input variable for the X-dimension of the 2d interpolation table (X-dimension is water level variability for MP29 establishment and mortality tables) 
    real(sp),intent(in) :: y                                     !   value of input variable for the Y-dimension of the 2d interpolation table (Y-dimension is salinity for MP29 establishment and mortality tables)  
    real(sp),intent(out) :: VALxy                                !   final interpolated value for input variables X=(x) and Y=(y)

    integer,intent(in) :: nXcols                                 !   number of columns for the X-dimension of the 2d interpolation table
    integer,intent(in) :: nYrows                                 !   number of rows for the Y-dimension of the 2d interpolation table
    real(sp),dimension(nXcols),intent(in) :: Xcols               !   column header values for the X-dimension of the 2d interpolation table - ASCENDING ORDER FROM LEFT TO RIGHT
    real(sp),dimension(nYrows),intent(in) :: Yrows               !   row header values for the Y-dimension of the 2d interpolation table - ASCENDING ORDER FROM TOP TO BOTTOM
    real(sp),dimension(nXcols,nYrows),intent(in) :: table        !   2d table of values being interpolated across in the X and Y dimensions (establishment and mortality tables for each species in MP29)
    
    integer :: ib                                                !   iterator over: (nXcols) or (nYrows)
    integer :: closest_index                                     !   index corresponding to: the closest row to (x), or the closet column to (y)
    real(sp) :: min_dif                                          !   the smallest difference between: (x) and the closest column header value, or (y) and the closest row header value
    real(sp) :: dif                                              !   the difference between: (x) and any iterable column header value, or (y) and any iterable row header value
    
    integer :: left                                              !   index of column to the (left) of (x) in the 2d interpolation table     
    integer :: right                                             !   index of column to the (right) of (x) in the 2d interpolation table 
    integer :: above                                             !   index of row (above) (y) in the 2d interpolation table 
    integer :: below                                             !   index of row (below) (y) in the 2d interpolation table 
    real(sp) :: r                                                !   column header value to the (r)ight of (x) in the X-dimension
    real(sp) :: l                                                !   column header value to the (l)eft of (x) in the X-dimension
    real(sp) :: a                                                !   row header value (a)bove (y) in the Y-dimension
    real(sp) :: b                                                !   row header value (b)elow (y) in the Y-dimension
    real(sp) :: VALla                                            !   table lookup value for interpolation boundary to the (l)eft of (x)  & (a)bove (y)
    real(sp) :: VALra                                            !   table lookup value for interpolation boundary to the (r)ight of (x) & (a)bove (y)
    real(sp) :: VALlb                                            !   table lookup value for interpolation boundary to the (l)eft of (x)  & (b)elow (y)   
    real(sp) :: VALrb                                            !   table lookup value for interpolation boundary to the (r)ight of (x) & (b)elow (y)   
    real(sp) :: VALxa                          	                 !   VAL(x,a) = interpolated value at (x), at the upper boundary (a)bove (y) 
    real(sp) :: VALxb                                            !   VAL(x,b) = interpolated value at (x), at the lower boundary (b)elow (y)
    real(sp) :: x_int_wgt                                        !   interpolation weighting factor in the X-dimension
    real(sp) :: y_int_wgt                                        !   interpolation weighting factor in the Y-dimension



                                                                 ! Read Y-dimension of input 2d interpolation table and determine which rows are above and below the value of y and will be used as the bounds for interpolation in the Y-dimension
    closest_index = -9999                                        ! arbitary large negative initial value
    dif = 0                                                      ! set initial difference to zero
    min_dif = 9999                                               ! arbitary large initial value

    do ib = 1,nYrows                                             ! iterate over the rows in the Y-dimension of the input 2d interpolation table					    
        dif = abs(Yrows(ib) - y)                                 ! calculate the distance between y and the header value of the current row
        if (dif < min_dif) then                                  ! check if the value of the current row header is closer to y than previously closest row header
            closest_index = ib                                   ! if closer, store the current index as the closest to y
            min_dif = dif                                        ! calculate the magnitude of the distance between y and the row header value in the closest_index
        end if 
    end do

    if ( Yrows(closest_index) > y )  then                        ! if the min_dif row header value is greater than y:
        below = closest_index                                    !    - then the closest_index is below y when ascending the row header values from top to bottom
        above = below - 1                                        !    - set index for row above by moving back up the table (which is -1 in index values) (above index is always -1 value less than below index)
    elseif( Yrows(closest_index) < y ) then                      ! if the min_dif row header value is less than y
        above = closest_index                                    !    - then the closest_index is above y when ascending the row header values from top to bottom
        below = above + 1                                        !    - set index for row below by  moving down the table (which is +1 in index values) (below index is always +1 value greater than above index)
    else                                                         ! if the min_dif row header value is equal to y
        below = closest_index                                    !    - then the closest_index is at y and set to below
        above = below                                            !    - above and below are equal
    end if 


                                                                 ! Read X-dimension of input 2d interpolation table and determine which columns are to the right and to the left of the value of x and will be used as the bounds for interpolation in the X-dimension
    closest_index = -9999                                        ! arbitary large negative initial value
    dif = 0                                                      ! set initial difference to zero
    min_dif = 9999                                               ! arbitary large initial value

    do ib = 1,nXcols                                             ! iterate over the columns in the X-dimension of the input 2d interpolation table
        dif = abs(Xcols(ib) - x)                                 ! calculate the distance between x and the header value of the current column
        if (dif < min_dif) then                                  ! check if the value of the current column header is closer to x than previously closest column header
            closest_index = ib                                   ! if closer, store the current index as the closest to x
            min_dif = dif                                        ! calculate the magnitude of the distance between x and the column header value in the closest_index
        end if 
    end do

    if ( Xcols(closest_index) > x )  then                        ! if the min_dif column header value is greater than x:
        right = closest_index                                    !    - then the closest_index is to the right of x when ascending the column header values from left to right
        left = right - 1                                         !    - set index for column to the left by moving back one step (left index is always -1 value less than right index)
    elseif( Xcols(closest_index) < x ) then                      ! if the min_dif column header value is less than x
        left = closest_index                                     !    - then the closest_index is to the left of x when ascending the column header values from left to right
        right = left + 1                                         !    - set index for column to the right by moving forward one index step (right index is always +1 value greater than left index)
    else                                                         ! if the min_dif column header value is equal to x
        left = closest_index                                     !    - then the closest_index is at x and set to left
        right = left                                             !    - right and left are equal
    end if 


    r = Xcols(right)                                             ! set column header value to the right of the interpolation bounds in X-dimension
    l = Xcols(left)                                              ! set column header value to the left of the interpolation bounds in X-dimension
    b = Yrows(below)                                             ! set row header value below the interpolation bounds in Y-dimension
    a = Yrows(above)                                             ! set row header value above the interpolation bounds in Y-dimension
    
    VALla = table(left,above)                                    ! lookup value from table at the upper left interpolation bound to the left of and above x
    VALra = table(right,above)                                   ! lookup value from table at the upper right interpolation bound to the right of and above x
    VALlb = table(left,below)                                    ! lookup value from table at the lower left interpolation bound to the left of and below x
    VALrb = table(right,below)                                   ! lookup value from table at the lower right interpolation bound to the right of and below x
    
    if (r==l) then                                               ! if the right and left bounding values for interpolation in the X-dimension are the same
        x_int_wgt = 0.0                                          !    - set the X-dimensional interpolation weighting factor to zero
    else                                                         ! if the right and left bounding values for interpolation in the X-dimension are not the same
        x_int_wgt = (x-l)/(r-l)                                  !     - scale the magnitude from left to right by the distance between the left bound and x
    endif

    if (a==b) then                                               ! if the below and above bounding values for interpolation in the Y-dimension are the same
        y_int_wgt = 0.0	                                         !    - set the Y-dimensional interpolation weighting factor to zero
    else                                                         ! if the below and above bounding values for interpolation in the Y-dimension are not the same
        y_int_wgt = (y-a)/(b-a)                                  !     - scale the magnitude from above to below by the distance between the above bound and y
    endif

    VALxa = VALla + (VALra - VALla) * x_int_wgt                 ! interpolate in the X-dimension for the row above y: VAL(x,a)
    VALxb = VALlb + (VALrb - VALlb) * x_int_wgt                 ! interpolate in the X-dimension for the row below y: VAL(x,b)
    VALxy = VALxa + (VALxb - VALxa) * y_int_wgt                 ! interpolate in the Y-dimension between the above value at x [VAL(x,a)] and the below value [VAL(x,b)] at x

end