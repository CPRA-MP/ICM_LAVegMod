subroutine set_io
    
    !  subroutine that reads in configuration file that contains input/output files and settings
    !  refer to PARAMS for description of all variables.
    
    use params
    implicit none
    
    call params_alloc_io
    
    open(unit=001, file=trim(adjustL('veg/LAVegMod_input_params.csv')))
 
    ! settings and input files that are static for all years
    read(001,*) start_year,dump_txt
    read(001,*) elapsed_year,dump_txt
    read(001,*) ncov,dump_txt    
    read(001,*) ngrid,dump_txt
    read(001,*) ncomp,dump_txt
    read(001,*) grid_res,dump_txt
    read(001,*) dem_res,dump_txt
    read(001,*) grid_file,dump_txt
    read(001,*) build_neighbors
    read(001,*) nearest_neighbors_file
    read(001,*) nearest_neighbors_dist
    read(001,*) near_neighbors_file
    read(001,*) near_neighbors_dist
    read(001,*) max_neighbors
    read(001,*) coverage_attribute_file
    
    ! input files that change for each model year
    read(001,*) veg_in_file,dump_txt        ! we can have this set automatically via the elapsed_year variable internal to the code instead of having it written to an input file
    read(001,*) hydro_comp_out_file 
    read(001,*) morph_grid_out_file
    ! output files     
    read(001,*) veg_out_file,dump_txt       ! we can have this set automatically via the elapsed_year variable internal to the code instead of having it written to an input file
    read(001,*) veg_summary_file,dump_txt   ! we can have this set automatically via the elapsed_year variable internal to the code instead of having it written to an input file
    ! filenaming convention
    read(001,*) fnc_tag
    
    fnc_tag =trim(adjustL(fnc_tag))
    mterm = fnc_tag(1:6)
    sterm = fnc_tag(8:10)
    gterm = fnc_tag(12:15)
    cterm = fnc_tag(17:20)
    uterm = fnc_tag(22:24)
    vterm = fnc_tag(26:28)
    

    
    close(001)

    return
end
