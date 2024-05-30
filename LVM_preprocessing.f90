subroutine preprocessing

    use params
    implicit none

    ! local variables
    integer :: i                                                                                    ! iterator
    integer :: g                                                                                    ! local ICM-LAVegMod grid ID variable
    integer :: ic                                                                                   ! iterator over coverage types
    integer :: y                                                                                    ! iterator over establishment/mortality table columns

    
    ! initialize grid data arrays to zero before reading in
    FFIBS_score = 0.0
    pct_vglnd_BLHF = 0.0
    pct_vglnd_SWF = 0.0
    pct_vglnd_FM = 0.0
    pct_vglnd_IM = 0.0
    pct_vglnd_BM = 0.0
    pct_vglnd_SM = 0.0
    
    ! read ICM-LAVegMod model attributes for each coverage type
    write(  *,*) ' - reading in model attributes for coverage types'
    write(000,*) ' - reading in model attributes for coverage types'
    cov_grp = 0
    cov_disp_class = 0
    FFIBS = 0
    open(unit=100, file=trim(adjustL(coverage_attribute_file)))
    read(100,'(A)') dump_txt                                                                        ! dump column header row 
    do i = 1,ncov
        read(100,*) cov_symbol(i), dump_txt, dump_txt, dump_txt, dump_txt, cov_grp(i), cov_disp_class(i), FFIBS(i)
    end do
    close(100)
    
    ! process coverage groups and save indices for specifc coverage groups
    flt_thn_cnt = 0                                                                                 ! initialize counter for finding thin mat flotant coverages
    flt_thn_indices = 0                                                                             ! initialize array to store coverage group index of FLOTANT coverage types in coverages(ngrid,ncov)
    flt_thk_cnt = 0                                                                                 ! initialize counter for finding thin mat flotant coverages
    flt_thk_indices = 0                                                                             ! initialize array to store coverage group index of FLOTANT coverage types in coverages(ngrid,ncov)
    
    do ic = 1,ncov
        if (cov_grp(ic) == 0) then
            wti = ic                                                                                ! coverage group index of WAT in coverages(ngrid,ncov)
        else if (cov_grp(ic) == 1) then
            nmi = ic                                                                                ! coverage group index of NOTMOD in coverages(ngrid,ncov)
        else if (cov_grp(ic) == 2) then
            boi = ic                                                                                ! coverage group index of BAREGRND_OLD in coverages(ngrid,ncov)
        else if (cov_grp(ic) == 3) then
            bni = ic                                                                                ! coverage group index of BAREGRND_NEW in coverages(ngrid,ncov)
        else if (cov_grp(ic) == 4) then                                                             ! coverage group dimension indices for flotants in coverages(ngrid,ncov)
            flt_thn_cnt = flt_thn_cnt + 1
            flt_thn_indices(flt_thn_cnt) = ic                                                       ! store coverage group indices for thin flotants in new array
        else if (cov_grp(ic) == 5) then                                                             ! coverage group dimension indices for flotants in coverages(ngrid,ncov)
            flt_thk_cnt = flt_thk_cnt + 1
            flt_thk_indices(flt_thk_cnt) = ic                                                       ! store coverage group indices for thick flotants in new array
        else if (cov_grp(ic) == 6) then
            bfi = ic                                                                                ! coverage group index of BARE_Flt in coverages(ngrid,ncov)
        else if (cov_grp(ic) == 7) then
            dfi = ic                                                                                ! coverage group index of DEAD_Flt in coverages(ngrid,ncov)
        end if
    end do
    
    ! read ICM-LAVegMod grid cell attributes
    grid_x = 0.0                                                                                    ! initialize data array to zero before reading in
    grid_y = 0.0                                                                                    ! initialize data array to zero before reading in
    grid_a = 0.0                                                                                    ! initialize data array to zero before reading in
    grid_comp = 0                                                                                   ! initialize data array to zero before reading in
    
    write(  *,*) ' - reading in grid cell attributes'
    write(000,*) ' - reading in grid cell attributes'
    open(unit=101, file=trim(adjustL(grid_file)))
    read(101,1234) dump_txt                                                                         ! dump column header row
    do i = 1,ngrid
        read(101,*) g, grid_x(g), grid_y(g), grid_a(g), grid_comp(g)
        dem_pixel_proportion(g) = dem_res**2 / grid_a(g)
    end do
    close(101)
    

    
    ! read ICM-LAVegMod grid output file into arrays - read in last year output veg coverages into initial values
    coverages = 0.0                                                                                 ! initialize grid data arrays to zero before reading in

    write(  *,*) ' - reading in LAVegMod grid-level output from previous year'
    write(000,*) ' - reading in LAVegMod grid-level output from previous year'
    open(unit=102, file=trim(adjustL(veg_in_file)))
    read(102,1234) veg_coverage_file_header                                                         ! dump column header row ! format 1234 must match structure of veg_out_file column headers
    do i = 1,ngrid
        read(102,*) g, coverages(g,:)
    end do
    close(102)
    
    read(veg_coverage_file_header,*) dump_txt,cov_symbol_check
    
    ! check that the input files are in the same species coverage column/row order - abort model run if orders do not match
    do i=1,ncov
        if (cov_symbol(i) /= cov_symbol_check(i)) then
            write(*,'(4A)') '    - column-order of ',trim(adjustL(veg_in_file)),' DOES NOT match row-order of ', trim(adjustL(coverage_attribute_file))
            write(*,'(A)') '    - EXITING SIMULATION NOW! Correct input files and re-submit run.'
            stop
        end if
    end do
    write(*,'(4A)') '    - column-order of ',trim(adjustL(veg_in_file)),' matches row-order of ', trim(adjustL(coverage_attribute_file))
    write(*,'(A)')  '    - continuing on with run'
    
   ! read ICM-Hydro compartment hydro output data in from file
    stg_mx_yr = 0.0                                                                                 ! initialize data array to zero before reading in
    stg_av_yr = 0.0                                                                                 ! initialize data array to zero before reading in
    stg_av_smr = 0.0                                                                                ! initialize data array to zero before reading in
    wlv_smr = 0.0                                                                                   ! initialize data array to zero before reading in
    sal_av_yr = 0.0                                                                                 ! initialize data array to zero before reading in
    sal_av_yr = 0.0                                                                                 ! initialize data array to zero before reading in
    sal_mx_14d_yr = 0.0                                                                             ! initialize data array to zero before reading in
    tmp_av_yr = 0.0                                                                                 ! initialize data array to zero before reading in   
    tmp_av_smr = 0.0                                                                                ! initialize data array to zero before reading in  

    write(  *,*) ' - reading in annual ICM-Hydro compartment-level output'
    write(000,*) ' - reading in annual ICM-Hydro compartment-level output'
    open(unit=103, file=trim(adjustL(hydro_comp_out_file)))
    read(103,*) dump_txt                                                                            ! dump header
    do i = 1,ncomp
        read(103,*) dump_txt,               &
   &         stg_mx_yr(i),                  &
   &         stg_av_yr(i),                  &
   &         stg_av_smr(i),                 &
   &         wlv_smr(i),                    &
   &         sal_av_yr(i),                  &
   &         sal_av_smr(i),                 &
   &         sal_mx_14d_yr(i),              &
   &         tmp_av_yr(i),                  &
   &         tmp_av_smr(i),                 &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt,                      &
   &         dump_flt
    end do
    close(103)
    
    ! read ICM-Morph landscape data for grid in from file
    grid_elev = 0.0                                                                                 ! initialize data array to zero before reading in
    water_from_morph = 0.0                                                                          ! initialize data array to zero before reading in
    
    write(  *,*) ' - reading in annual ICM-Morph landscape data output'
    write(000,*) ' - reading in annual ICM-Morph landscape data output'
    open(unit=104, file=trim(adjustL(morph_grid_out_file)))
    read(104,*) dump_txt        ! dump header
    do i = 1,ngrid
        read(104,*) g,                          &                                                   ! grid cell ID
   &            dump_flt,                       &                                                   ! elevation of water bottom portion of grid cell, as calculated in ICM-Morph    
   &            grid_elev(g),                   &                                                   ! elevation of land portion of grid cell, as calculated in ICM-Morph
   &            dump_flt,                       &                                                   ! percent land (wetland+upland) of grid cell, as calculated in ICM-Morph, in the output grid file this is from 0-100%
   &            dump_flt,                       &                                                   ! percent wetland of grid cell, as calculated in ICM-Morph, in the output grid file this is from 0-100%
   &            water_from_morph(g)                                                                 ! percent water of grid cell, as calculated in ICM-Morph, in the output grid file this is from 0-100%
    end do
    close(104)
    
    water_from_morph = water_from_morph/100.0                                                       ! convert from 0-100% to 0.0-1.0 proportion

    ! read barrier island domain
    barrier_island = 0                                                                              ! initialize data array to zero before reading in
    
    write(  *,*) ' - reading in barrier island domain map'
    write(000,*) ' - reading in barrier island domain map'
    open(unit=105, file=trim(adjustL(morph_grid_out_file)))
    read(105,*) dump_txt                                                                            ! dump header
    do i = 1,ngrid
        read(105,*) g,                          &                                                   ! grid cell ID
   &            barrier_island(g)                                                                   ! barrier island flag (1 if grid cell is in barrier island domain; 0, if not)
    end do
    close(105)
    
    ! read tree establishment criteria
    tree_establishment = 0                                                                           ! initialize data array to zero before reading in
    
    write(  *,*) ' - reading in tree establishment criteria'
    write(000,*) ' - reading in tree establishment criteria'
    open(unit=106, file=trim(adjustL(morph_grid_out_file)))
    read(106,*) dump_txt                                                                            ! dump header
    do i = 1,ngrid
        read(106,*) g,                          &                                                   ! grid cell ID
   &            tree_establishment(g)                                                               ! tree establishment criteria (1 if conditions are met; 0, if not)
    end do
    close(106)
    
    ! read species establishment probability tables
    est_X_bins = 0.0                                                                                ! initialize data array to zero before reading in
    est_Y_bins = 0.0                                                                                ! initialize data array to zero before reading in
    establish_tables = 0.0                                                                          ! initialize data array to zero before reading in
    
    write(  *,*) ' - reading in probability of establishment tables'
    write(000,*) ' - reading in probability of establishment tables'
    do ic = 1,ncov
        open(unit=107, file=trim(adjustL('veg/tables/establishment_'//trim(adjustL(cov_symbol(ic)))))//'.csv' )
        read(107,*) dump_txt,est_X_bins(:,ic)                                                       ! first row of establisment table is the X-variable values defining each establishment 'bin'
        do y = 1,n_Y_bins
            read(107,*) est_Y_bins(y,ic),       &                                                   ! first column of establishment table is Y-variable value defining each establishment 'bin'
   &            establish_tables(:,y,ic)                                                            ! read the rest of the columns into the establishment probability array for given y value and coverage type, ic
        end do
        close(107)
    end do
    
    ! read species mortality probability tables
    mort_X_bins = 0.0                                                                                ! initialize data array to zero before reading in
    mort_Y_bins = 0.0                                                                                ! initialize data array to zero before reading in
    mortality_tables = 0.0                                                                          ! initialize data array to zero before reading in
    
    write(  *,*) ' - reading in probability of mortality tables'
    write(000,*) ' - reading in probability of mortality tables'
    do ic = 1,ncov
        open(unit=108, file=trim(adjustL('veg/tables/mortality_'//trim(adjustL(cov_symbol(ic)))))//'.csv' )
        read(108,*) dump_txt,mort_X_bins(:,ic)                                                      ! first row of mortality table is the X-variable values defining each mortality 'bin'
        do y = 1,n_Y_bins
            read(108,*) mort_Y_bins(y,ic),       &                                                  ! first column of mortality table is Y-variable value defining each mortality 'bin'
   &            mortality_tables(:,y,ic)                                                            ! read the rest of the columns into the mortality probability array for given y value and coverage type, ic
        end do
        close(107)
    end do
    
    
1234    format(A,<ncov>(',',A))

    return
end
