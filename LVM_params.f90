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
    integer :: ngrid                                                ! number of ICM-LAVegMod grid cells
    integer :: ncomp                                                ! number of ICM-Hydro compartments
    integer :: grid_res                                             ! XY resolution of ICM-LAVegMod grid (meters) - only applicable for regular Cartesian grid
    integer :: dem_res                                              ! XY resolution of DEM (meters)

    ! input files in subroutine: SET_IO
    character*fn_len :: veg_in_file                                 ! file name, with relative path, to *vegty.csv file from previous model year read in to set initial conditions for the current model year
    character*fn_len :: hydro_comp_out_file                         ! file name, with relative path, to *compartment_out.csv from current model year's ICM-Hydro simulation
    
    ! output files in subroutine: SET_IO
    character*fn_len :: veg_out_file                                ! file name, with relative path, to *vegty.csv file for current year written to disc for final landscape of the current model year
    
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
    integer,dimension(:),allocatable ::  grid_x                     ! X coordinate of ICM-LAVegMod grid cell centroid (UTM Zone 15N meters)
    integer,dimension(:),allocatable ::  grid_y                     ! Y coordinate of ICM-LAVegMod grid cell centroid (UTM Zone 15N meters) 

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
    real(sp),dimension(:),allocatable :: water_from_morph           ! percent of ICM-LAVegMod grid cell that is water, as calculated at end of previous year's ICM-Morph run
    
    ! species coverage grid in: PREPROCESSING
    ! define variables used to define the vegetation species coverage at each grid cell that are read in from file
    ! these variables are 2D arrays [i,j] where the ith dimension represents the grid cell ID and the jth dimension represents the species coverage for the previous year [j=1] and for the current model year [j=2]
    character*3000 :: veg_coverage_file_header                          ! text string that saves the first row of the veg input file to use as a header in the output file
    real(sp),dimension(:,:),allocatable :: water                        ! percent of ICM_LAVegMod grid cell that is water (j=1 will be last year's ICM-LAVegMod % value, where as j=2 will store current year's % value as used by ICM-LAVegMod; not to be confused with 'water_from_morph' variable)
    real(sp),dimension(:,:),allocatable :: upland                       ! percent of ICM-LAVegMod grid cell that is upland/developed (e.g., NotMod) and is too high and dry for wetland vegetation
    real(sp),dimension(:,:),allocatable :: bare_old                     ! percent of ICM-LAVegMod grid cell that is non-vegetated wetland and was bare in previous year (old bare ground)
    real(sp),dimension(:,:),allocatable :: bare_new                     ! percent of ICM-LAVegMod grid cell that is non-vegetated wetland and is newly bare during current year (new bare ground)
    real(sp),dimension(:,:),allocatable :: QULA3                        !
    real(sp),dimension(:,:),allocatable :: QULE                         !
    real(sp),dimension(:,:),allocatable :: QUNI                         !
    real(sp),dimension(:,:),allocatable :: QUTE                         !
    real(sp),dimension(:,:),allocatable :: QUVI                         !
    real(sp),dimension(:,:),allocatable :: ULAM                         !
    real(sp),dimension(:,:),allocatable :: NYAQ2                        !
    real(sp),dimension(:,:),allocatable :: SANI                         !
    real(sp),dimension(:,:),allocatable :: TADI2                        !
    real(sp),dimension(:,:),allocatable :: ELBA2_Flt                    !
    real(sp),dimension(:,:),allocatable :: PAHE2_Flt                    !
    real(sp),dimension(:,:),allocatable :: bare_Flt                     !
    real(sp),dimension(:,:),allocatable :: dead_Flt                     !
    real(sp),dimension(:,:),allocatable :: COES                         !
    real(sp),dimension(:,:),allocatable :: MOCE2                        !
    real(sp),dimension(:,:),allocatable :: PAHE2                        !
    real(sp),dimension(:,:),allocatable :: SALA2                        !
    real(sp),dimension(:,:),allocatable :: ZIMI                         !
    real(sp),dimension(:,:),allocatable :: CLMA10                       !
    real(sp),dimension(:,:),allocatable :: ELCE                         !
    real(sp),dimension(:,:),allocatable :: IVFR                         !
    real(sp),dimension(:,:),allocatable :: PAVA                         !
    real(sp),dimension(:,:),allocatable :: PHAU7                        !
    real(sp),dimension(:,:),allocatable :: POPU5                        !
    real(sp),dimension(:,:),allocatable :: SALA                         !
    real(sp),dimension(:,:),allocatable :: SCCA11                       !
    real(sp),dimension(:,:),allocatable :: TYDO                         !
    real(sp),dimension(:,:),allocatable :: SCAM6                        !
    real(sp),dimension(:,:),allocatable :: SCRO5                        !
    real(sp),dimension(:,:),allocatable :: SPCY                         !
    real(sp),dimension(:,:),allocatable :: SPPA                         !
    real(sp),dimension(:,:),allocatable :: AVGE                         !
    real(sp),dimension(:,:),allocatable :: DISP                         !
    real(sp),dimension(:,:),allocatable :: JURO                         !
    real(sp),dimension(:,:),allocatable :: SPAL                         !
    real(sp),dimension(:,:),allocatable :: BAHABI                       !
    real(sp),dimension(:,:),allocatable :: DISPBI                       !
    real(sp),dimension(:,:),allocatable :: PAAM2                        !
    real(sp),dimension(:,:),allocatable :: SOSE                         !
    real(sp),dimension(:,:),allocatable :: SPPABI                       !
    real(sp),dimension(:,:),allocatable :: SPVI3                        !
    real(sp),dimension(:,:),allocatable :: STHE9                        !
    real(sp),dimension(:,:),allocatable :: UNPA                         !
    real(sp),dimension(:,:),allocatable :: FFIBS_score                  ! weighted FFIBS score of ICM-LAVegMod grid cell - used for accretion
    real(sp),dimension(:,:),allocatable :: pct_vglnd_BLHF               ! percent of vegetated land that is bottomland hardwood forest
    real(sp),dimension(:,:),allocatable :: pct_vglnd_SWF                ! percent of vegetated land that is swamp forest
    real(sp),dimension(:,:),allocatable :: pct_vglnd_FM                 ! percent of vegetated land that is fresh (attached) marsh
    real(sp),dimension(:,:),allocatable :: pct_vglnd_IM                 ! percent of vegetated land that is intermediate marsh
    real(sp),dimension(:,:),allocatable :: pct_vglnd_BM                 ! percent of vegetated land that is brackish marsh
    real(sp),dimension(:,:),allocatable :: pct_vglnd_SM                 ! percent of vegetated land that is saline marsh

end module params
