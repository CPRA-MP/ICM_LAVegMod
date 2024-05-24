module params

!! module to define parameter types for all global variables

    implicit none

    ! generic variables used across all subroutines
    integer,parameter :: sp=selected_real_kind(p=6)                 ! determine compiler KIND value for 4-byte (single precision) floating point real numbers
    integer,parameter :: dp=selected_real_kind(p=15)                ! determine compiler KIND value for 8-byte (double precision) floating point real numbers
    integer,parameter :: fn_len=300                                 ! maximum character length allowed for filename character strings read in from input_params.csv
    character*3000 :: dump_txt                                      ! dummy variable to use for skipping lines in input files
    integer :: dump_int                                             ! dummy variable to use for data in input files
    real(sp) :: dump_flt                                            ! dummy variable to use for data in input files
    integer :: tp                                                   ! flag to indicate which timeperiod to use for inundation calculations (1-12=month; 13 = annual)

    ! I/O settings in subroutine: SET_IO
    integer :: start_year                                           ! first year of model run
    integer :: elapsed_year                                         ! elapsed year of model simulation
    integer :: ncov                                                 ! number of grid cell coverages included in the model (e.g., water, new_bareground, old_bareground, UNPA, SPPA, etc.)
    integer :: ngrid                                                ! number of ICM-LAVegMod grid cells
    integer :: ncomp                                                ! number of ICM-Hydro compartments
    integer :: dem_res                                              ! XY resolution of DEM (meters)
    integer :: build_neighbors                                      ! flag - set to 1 if near and nearest neighbor lists need to be built from Grid XY data, set to 0 if neighbor files already exist
    integer :: nearest_neighbors_dist                               ! distance in which a neighboring grid cell is considered a nearest neighbor (meters) *must be smaller magnitude than "near_neighbor_dist"*
    integer :: near_neighbors_dist                                  ! distance in which a neighboring grid cell is considered a near neighbor (meters) *must be larger magnitude than "nearest_neighbor_dist"*
    integer :: max_neighbors                                        ! maximum number of grid cells that will be allowed in the near and nearest neighbor lists
    integer :: n_X_bins                                             ! number of bins definiing the X-axis of the establishment and mortability input tables
    integer :: n_Y_bins                                             ! number of bins definiing the Y-axis of the establishment and mortability input tables
    
    ! input files in subroutine: SET_IO
    character*fn_len :: coverage_attribute_file                     ! file name, with relative path, to csv with model attributes for each coverage type - this file row-order must match the column-order of veg_in_file, below
    character*fn_len :: grid_file                                   ! file name, with relative path, to csv with X and Y coordinates (UTM meters) of grid cell centroids - also includes the grid cell area (sq meters) and the overlaid ICM-Hydro compartment ID
    character*fn_len :: nearest_neighbors_file                      ! file name, with relative path, to csv with list of grid cells that are the defined nearest neighbors
    character*fn_len :: near_neighbors_file                         ! file name, with relative path, to csv with list of grid cells that are the defined near neighbors
    character*fn_len :: veg_in_file                                 ! file name, with relative path, to *vegty.csv file from previous model year read in to set initial conditions for the current model year
    character*fn_len :: hydro_comp_out_file                         ! file name, with relative path, to *compartment_out.csv from current model year's ICM-Hydro simulation
    character*fn_len :: morph_grid_out_file                         ! file name, with relative path, to *grid_data.csv file from previous model year's ICM-Morph simulation
    
    ! output files in subroutine: SET_IO
    character*fn_len :: veg_out_file                                ! file name, with relative path, to *vegty.csv file for current year written to disk for final landscape of the current model year
    character*fn_len :: veg_summary_file                            ! file name, with relative path, to *vegsm.csv file for current year written to disk for final landscape of the current model year
    
    ! QAQC save point information in subroutine: SET_IO
    character*fn_len :: fnc_tag                                     ! file naming convention tag
    character*6 :: mterm                                            ! file naming convention model name term
    character*3 :: sterm                                            ! file naming convention scenario term
    character*4 :: gterm                                            ! file naming convention group term
    character*4 :: cterm                                            ! file naming convention CLARA scenario term
    character*3 :: uterm                                            ! file naming convention uncertainty term
    character*3 :: vterm                                            ! file naming convention variance term

    ! define variables read in or calculated from files in subroutine: PREPROCESSING
    integer,dimension(:),allocatable ::  grid_comp                  ! ICM-Hydro compartment ID overlaying ICM-LAVegMod grid (-)
    real(sp),dimension(:),allocatable ::  grid_x                    ! X coordinate of ICM-LAVegMod grid cell centroid (UTM Zone 15N meters)
    real(sp),dimension(:),allocatable ::  grid_y                    ! Y coordinate of ICM-LAVegMod grid cell centroid (UTM Zone 15N meters)
    real(sp),dimension(:),allocatable :: grid_a                     ! area of ICM_LAVegMod grid cell (sq meters)
    real(sp),dimension(:),allocatable :: dem_pixel_proportion       ! proportion of each ICM-LAVegGrid cell that is occupied by ONE ICM-Morph DEM pixel (-); for the 2023 grid this was equal to 1/256 = (30*30)/(480*480)
    
    ! define coverage attribute variables read in from input attribute table in subroutine: PREPROCESSING
    character*20,dimension(:),allocatable ::  cov_symbol            ! USDA code/symbol for each vegetation coverage type (e.g., SPPA, SPAL, etc.) - from *coverage_attribute_file*
    character*20,dimension(:),allocatable ::  cov_symbol_check      ! USDA code/symbol for each vegetation coverage type used to check consistency in input files - from *veg_in_file*
    integer,dimension(:),allocatable ::  cov_grp                    ! model group ID for LAVegMod process model assigned to each respective coverage type
                                                                    !       cov_grp =  0; water
                                                                    !       cov_grp =  1; not modeled/developed
                                                                    !       cov_grp =  2; old bareground
                                                                    !       cov_grp =  3; new bareground
                                                                    !       cov_grp =  4; flotant marsh - thin mat
                                                                    !       cov_grp =  5; flotant marsh - thick mat
                                                                    !       cov_grp =  6; flotant marsh - bare mat
                                                                    !       cov_grp =  7; flotant marsh - dead
                                                                    !       cov_grp =  8; bottomland hardwood forest
                                                                    !       cov_grp =  9; swamp forest
                                                                    !       cov_grp = 10; fresh emergent wetland vegetation
                                                                    !       cov_grp = 11; intermediate emergent wetland vegetation
                                                                    !       cov_grp = 12; brackish emergent wetland vegetation
                                                                    !       cov_grp = 13; saline emergent wetland vegetation
                                                                    !       cov_grp = 14; barrier island vegetation
    integer,dimension(:),allocatable ::  cov_disp_class             ! disperal class ID for each respective coverage type
                                                                    !       cov_disp_class =  0; no dispersal - used for non-vegetative coverage types (e.g., water)
                                                                    !       cov_disp_class =  1; nearest disperal - species can disperse only from "nearest neighboring" areas (distance is assigned in SET_IO from "veg/LAVegMod_input_params.csv")
                                                                    !       cov_disp_class =  2; near disperal - species can disperse only from "near and nearest neighboring" areas (distances are assigned in SET_IO from "veg/LAVegMod_input_params.csv")
                                                                    !       cov_disp_class =  3; always available - "weedy" species that are assumed always available for establishment/infinite dispersal distance
    real(sp),dimension(:),allocatable ::  FFIBS                     ! FFIBS score assigned to each respective coverage type
    integer :: wti                                                  ! index in coverages(ngrid,ncov,2) for water coverage group
    integer :: nmi                                                  ! index in coverages(ngrid,ncov,2) for NotMod coverage group
    integer :: boi                                                  ! index in coverages(ngrid,ncov,2) for old bareground coverage group
    integer :: bni                                                  ! index in coverages(ngrid,ncov,2) for new bareground coverage group
    integer :: dfi                                                  ! index in coverages(ngrid,ncov,2) for dead flotant marsh coverage group
    integer :: bfi                                                  ! index in coverages(ngrid,ncov,2) for bare mat flotant marsh coverage group
    integer :: flt_thn_cnt                                          ! count of species included in the thin mat flotant marsh coverage group, excluding dead flotant
    integer,dimension(:),allocatable ::  flt_thn_indices            ! 1D array that stores the coverage group indices of thin mat flotant marsh coverage types in coverages(ngrid,ncov,2), excluding dead flotant
    integer :: flt_thk_cnt                                          ! count of species included in the thick mat flotant marsh coverage group, excluding dead flotant
    integer,dimension(:),allocatable ::  flt_thk_indices            ! 1D array that stores the coverage group indices of thin mat flotant marsh coverage types in coverages(ngrid,ncov,2), excluding dead flotant


    ! species coverage grid in: PREPROCESSING
    ! define variables used to define the vegetation species coverage at each grid cell that are read in from file
    ! these variables are 3D arrays [i,j,k] where the ith dimension represents the grid cell ID and the jth dimension represents the coverage type column, and the kth dimension represents the coverage value of type j for the previous coverage state [k=1] and for the current coverage state [j=2]
    character*3000 :: veg_coverage_file_header                      ! text string that saves the first row of the veg input file to use as a header in the output file
    real(sp),dimension(:,:,:),allocatable :: coverages              ! percent of ICM_LAVegMod grid cell that is each coverage type (k=1 will be the previous state ICM-LAVegMod % value, where as k=2 will store state % value as used by ICM-LAVegMod; not to be confused with 'water_from_morph' variable)

    
    ! define ICM-Hydro variables read in from compartment_out summary file in subroutine: PREPROCESSING
    real(sp),dimension(:),allocatable :: stg_mx_yr                  ! Maximum water surface elevation (stage) during the year (m NAVD88)
    real(sp),dimension(:),allocatable :: stg_av_yr                  ! Mean water surface elevation (stage) during the year (m NAVD88)
    real(sp),dimension(:),allocatable :: stg_av_smr                 ! Mean water surface elevation (stage) during growing season/summer  (m NAVD88) - growing season/summer defined as May 1 through Aug 31 (inclusive) defined in ICM-Hydro/2D_ICM_summaries.f 
    real(sp),dimension(:),allocatable :: wlv_smr                    ! Water level variability during growing season/summer (m) - water level variability calculated from tidal range (see ICM-Hydro/2D_ICM_summaries.f) - growing season/summer defined as May 1 through Aug 31 (inclusive) defined in ICM-Hydro (2D_ICM_summaries.f)   
    real(sp),dimension(:),allocatable :: sal_av_yr                  ! Mean water salinity during the year (ppt)    
    real(sp),dimension(:),allocatable :: sal_av_smr                 ! Mean water salinity during growing season/summer (ppt) - growing season/summer defined as May 1 through Aug 31 (inclusive) defined in ICM-Hydro/2D_ICM_summaries.f
    real(sp),dimension(:),allocatable :: sal_mx_14d_yr              ! Maximum 2-week mean salinity during the year (ppt)
    real(sp),dimension(:),allocatable :: tmp_av_yr                  ! Mean water temperture during the year (deg C)  
    real(sp),dimension(:),allocatable :: tmp_av_smr                 ! Mean water temperture during growing season/summer (deg C) - growing season/summer defined as May 1 through Aug 31 (inclusive) defined in ICM-Hydro/2D_ICM_summaries.f

    ! define variables read in from ICM-Morph output files in subroutine: PREPROCESSING
    real(sp),dimension(:),allocatable :: grid_elev                  ! elevation of land portion of grid cell, as calculated at end of previous year's ICM-Morph run (m, NAVD88)
    real(sp),dimension(:),allocatable :: water_from_morph           ! proportion of ICM-LAVegMod grid cell that is water, as calculated at end of previous year's ICM-Morph run (0.0 - 1.0)
    
    ! define variables read in in subroutine: PREPROCESSING
    integer,dimension(:),allocatable ::  barrier_island             ! flag indicating whether the grid cell is located in a barrier island domain or not (1 if island; 0 if not)
    integer,dimension(:),allocatable ::  tree_establishment         ! flag indicating whether the grid cell has met hydrologic tree establishment criteria for the year (1 if conditions met; 0 if not)
    
    ! define variables read in for establishment and mortability tables in subroutine: PREPROCESSING
    real(sp),dimension(:,:),allocatable :: est_X_bins               ! array holding the values used to define the X-axis of each species establishent tables - the first dimension is the location in the X-axis, the second dimension is the coverage index, ic
    real(sp),dimension(:,:),allocatable :: est_Y_bins               ! array holding the values used to define the Y-axis of each species establishent tables - the first dimension is the location in the X-axis, the second dimension is the coverage index, ic
    real(sp),dimension(:,:,:),allocatable :: establish_tables       ! 2-dimensional establishment probablity table for each species - first dimension is X value of table, second dimension is Y value, third dimension is the coverage index, ic
   
    
    
    ! these variables are 2D arrays [i,j] where the ith dimension represents the grid cell ID and the jth dimension represents the species coverage for the previous year [j=1] and for the current model year [j=2]
    real(sp),dimension(:,:),allocatable :: FFIBS_score              ! weighted FFIBS score of ICM-LAVegMod grid cell - used for accretion
    real(sp),dimension(:,:),allocatable :: pct_vglnd_BLHF           ! percent of vegetated land that is bottomland hardwood forest
    real(sp),dimension(:,:),allocatable :: pct_vglnd_SWF            ! percent of vegetated land that is swamp forest
    real(sp),dimension(:,:),allocatable :: pct_vglnd_FM             ! percent of vegetated land that is fresh (attached) marsh
    real(sp),dimension(:,:),allocatable :: pct_vglnd_IM             ! percent of vegetated land that is intermediate marsh
    real(sp),dimension(:,:),allocatable :: pct_vglnd_BM             ! percent of vegetated land that is brackish marsh
    real(sp),dimension(:,:),allocatable :: pct_vglnd_SM             ! percent of vegetated land that is saline marsh

    ! define variables read in or calculated from files in subroutine: NEIGHBORS
    integer,dimension(:,:),allocatable ::  nearest_neighbors        ! list of grid cell IDs that are the nearest neighbors to each grid cell
    integer,dimension(:,:),allocatable ::  near_neighbors           ! list of grid cell IDs that are the near neighbors to each grid cell
    
    ! define variables read in or calculated from files in subroutine: establishment_Pcalc
    real(sp),dimension(:,:),allocatable :: establish_p              ! the establishment probability for every cell and species in this model year based on the species-specific variables

end module params
