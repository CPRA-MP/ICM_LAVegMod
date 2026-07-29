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
    !    - the upper left bound is (l)ess than (x), and (l)ess than (y)
    !    - the upper right bound is (g)reater than (x), and (l)ess than (y)
    !    - the lower left bound is (l)ess than (x), and (g)reater than (y)
    !    - the lower right bound is (g)reater than (x), and (g)reater than (y)
    !               
    !               
    !                              lesser      x       greater
    !                         (-) ----------------------------------------(+)
    !               lesser     | VAL(l,l)   VAL(x,l)   VAL(g,l)     _ _           
    !                          |                                     |   (y-lesserY)
    !                  y       | VAL(l,y)  *VAL(x,y)*  VAL(g,y)     _|_              
    !                          |           
    !               greater    | VAL(l,g)   VAL(x,g)   VAL(g,g)           
    !                          |
    !                          | |--(x-lesserX)--|
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
    
    integer :: ilesserX                                          !   index of column with header value less than (x) in the 2d interpolation table     
    integer :: igreaterX                                         !   index of column with header value greater than (x) in the 2d interpolation table 
    integer :: ilesserY                                          !   index of row with header value less than (y) in the 2d interpolation table     
    integer :: igreaterY                                         !   index of row with header value greater than (x) in the 2d interpolation table
    real(sp) :: lesserX                                          !   column header value less than (x) in the X-dimension
    real(sp) :: greaterX                                         !   column  header value greater than (x) in the X-dimension
    real(sp) :: lesserY                                          !   row header value less than (y) in the Y-dimension
    real(sp) :: greaterY                                         !   row header value greater than (y) in the Y-dimension
    real(sp) :: VALll                                            !   table lookup value for interpolation boundary lesser than (x)  & lesser than (y)
    real(sp) :: VALgl                                            !   table lookup value for interpolation boundary greater than (x) & lesser than (y)
    real(sp) :: VALlg                                            !   table lookup value for interpolation boundary lesser than (x)  & greater than (y)   
    real(sp) :: VALgg                                            !   table lookup value for interpolation boundary greater than (x) & greater than (y)   
    real(sp) :: VALxl                          	                 !   VAL(x,a) = interpolated value at (x), at the lesser Y-dimension boundary less than (y) 
    real(sp) :: VALxg                                            !   VAL(x,b) = interpolated value at (x), at the greater Y-dimension boundary greater than (y)
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
        igreaterY = closest_index                                !    - then the closest_index is greater than y when ascending the row header values from top to bottom
        ilesserY = igreaterY - 1                                 !    - lesser index is always -1 value less than greater index
    elseif( Yrows(closest_index) < y ) then                      ! if the min_dif row header value is less than y
        ilesserY = closest_index                                 !    - then the closest_index is less than y when ascending the row header values from top to bottom
        igreaterY = ilesserY + 1                                 !    -  greater index is always +1 value more than lesser index
    else                                                         ! if the min_dif row header value is equal to y
        igreaterY = closest_index                                !    - then the closest_index is at y and set to lesser index
        ilesserY = igreaterY                                     !    - lesser and greater indices are equal
    end if 


                                                                 ! Read X-dimension of input 2d interpolation table and determine which columns are greater and less than the value of x and will be used as the bounds for interpolation in the X-dimension
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
        igreaterX = closest_index                                !    - then the closest_index is greater than x when ascending the column header values from left to right
        ilesserX = igreaterX - 1                                 !    - lesser index is always -1 value less than greater index
    elseif( Xcols(closest_index) < x ) then                      ! if the min_dif column header value is less than x
        ilesserX = closest_index                                 !    - then the closest_index is to the left of x when ascending the column header values from left to right
        igreaterX = ilesserX + 1                                 !    - greater index is always +1 value more than lesser index
    else                                                         ! if the min_dif column header value is equal to x
        ilesserX = closest_index                                 !    - then the closest_index is at x and set to lesser index
        igreaterX = ilesserX                                     !    - lesser and greater indices are equal
    end if 


    greaterX = Xcols(igreaterX)                                  ! set column header value that is the greater/upper of the interpolation bounds in X-dimension
    lesserX = Xcols(ilesserX)                                    ! set column header value that is the lesser/lower of the interpolation bounds in X-dimension
    greaterY = Yrows(igreaterY)                                  ! set row header value that is the greater/upper of the interpolation bounds in Y-dimension
    lesserY = Yrows(ilesserY)                                    ! set row header value  that is the lesser/lower of the interpolation bounds in Y-dimension
    
    VALll = table(ilesserX,ilesserY)                               ! lookup value from table at the upper left interpolation bound lesser than (x) and lesser than (y)
    VALgl = table(igreaterX,ilesserY)                              ! lookup value from table at the upper right interpolation bound greater than (x) and lesser than (y)
    VALlg = table(ilesserX,igreaterY)                              ! lookup value from table at the lower left interpolation bound lesser than (x) and greater than (y)
    VALgg = table(igreaterX,igreaterY)                             ! lookup value from table at the lower right interpolation bound greater than (x) and greater than (y)
    
    if (greaterX==lesserX) then                                  ! if the right and left bounding values for interpolation in the X-dimension are the same
        x_int_wgt = 0.0                                          !    - set the X-dimensional interpolation weighting factor to zero
    else                                                         ! if the right and left bounding values for interpolation in the X-dimension are not the same
        x_int_wgt = (x-lesserX)/(greaterX-lesserX)               !     - scale the magnitude from left to right by the distance between the left bound and x
    endif

    if (greaterY==lesserY) then                                  ! if the below and above bounding values for interpolation in the Y-dimension are the same
        y_int_wgt = 0.0	                                         !    - set the Y-dimensional interpolation weighting factor to zero
    else                                                         ! if the below and above bounding values for interpolation in the Y-dimension are not the same
        y_int_wgt = (y-lesserY)/(greaterY-lesserY)               !     - scale the magnitude from above to below by the distance between the above bound and y
    endif

    VALxl = VALll + (VALgl - VALll) * x_int_wgt                 ! interpolate in the X-dimension for the row index lesser than y: VAL(x,l)
    VALxg = VALlg + (VALgg - VALlg) * x_int_wgt                 ! interpolate in the X-dimension for the row index greater than y: VAL(x,g)
    VALxy = VALxl + (VALxg - VALxl) * y_int_wgt                 ! interpolate in the Y-dimension between the bounding values calculated at x, [VAL(x,l)] and [VAL(x,g)]

end