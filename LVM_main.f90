!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                  
!       ICM-LAVegMod - Vegetative Coverage Model   
!                                                  
!                                                  
!   Fortran version of ICM-LAVegMod developed         
!   for 2029 Coastal Master Plan - LA CPRA         
!                                                  
!   original model (2012 MP): Visser and Duke-Sylvester (2012)
!   revised model (2017 MP): Visser and Duke-Sylvester (2017)
!   revised model (2023 MP): Foster-Martinez et al. (2023)
!   current model (2029 MP): Foster-Martinez and White (in prep)
!                                              
!   Questions: eric.white@la.gov, mrfoster@uno.edu
!   last update: 12/15/2023
!                                                     
!   project site: https://github.com/CPRA-MP      
!   documentation: http://coastal.la.gov/our-plan  
!                                                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
program main
    
    use params
    use sort
    
    implicit none
    
    ! local variables
    integer,dimension(8) :: dtvalues                ! variable to store date time values
    
    character*17 :: dtstrf                          ! string to hold formatted datetime
    character*19 :: dtstr                           ! string to hold formatted datetime

    
    call date_and_time(VALUES=dtvalues)             ! grab simulation start time
    write(dtstrf,8888) dtvalues(1),dtvalues(2),dtvalues(3),'_',dtvalues(5),'.',dtvalues(6),'.',dtvalues(7)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    open(unit=000, file=trim(adjustL('log.log')))!'_ICM-LAVegMod_runlog_')//dtstr//trim('.log')))
    
    write(  *,*)
    write(  *,*) '*************************************************************'
    write(  *,*) '****                                                     ****'
    write(  *,*) '****    ****   STARTING ICM-LAVegMod SIMULATION   ***    ****'
    write(  *,*) '****                                                     ****'
    write(  *,*) '*************************************************************'
    write(  *,*)
    write(  *,*) 'Started ICM-LAVegMod simulation at: ',dtstr
    write(  *,*)

    write(000,*)
    write(000,*) '*************************************************************'
    write(000,*) '****                                                     ****'
    write(000,*) '****    ****   STARTING ICM-LAVegMod SIMULATION   ***    ****'
    write(000,*) '****                                                     ****'
    write(000,*) '*************************************************************'    
    write(000,*)
    write(000,*) 'Started ICM-LAVegMod simulation at: ',dtstr
    write(000,*)

    call set_io                                     ! input/output settings - must be run BEFORE parameter allocation   
    call params_alloc
   
    call preprocessing
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),'_',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Preprocessing subroutine ended at: ',dtstr
    write(000,*) 'Preprocessing subroutine ended at: ',dtstr

    call neighbors
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),'_',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Nearest neighbors subroutine ended at: ',dtstr
    write(000,*) 'Nearest neighbors subroutine ended at: ',dtstr

    ! Reset new bareground and adjust dead flotant 
    !call reset_coverages
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Reset Coverages subroutine ended at: ',dtstr
    write(000,*) 'Reset Coverages subroutine ended at: ',dtstr

    ! Adjust the total land and water in each cell based on changes from ICM-Morph
    !call land_change
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Land Change subroutine ended at: ',dtstr
    write(000,*) 'Land Change subroutine ended at: ',dtstr

    ! write intermediate coverage file for post updates from ICM-Morph - also write summary output file
    call write_output('imo_morphu',1)   !currently this is a 10*character string being passed in to write_output - currently needs to be padded with spaces
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Write Intermediate Output subroutine ended at: ',dtstr
    write(000,*) 'Write Intermediate Output subroutine ended at: ',dtstr 

    ! Calculate the establishment and mortality probabilities for every species and grid cell for this model year's environmental conditions
    !call mort_est_prob
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Mortality Establishment Probability subroutine ended at: ',dtstr
    write(000,*) 'Mortality Establishment Probability subroutine ended at: ',dtstr

    ! Allow high dispersal species (class three "weedy" species) to establish on any new bareground
    !call hhigh_disp_est_nbg
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'High Dispersal Establishment on New Bareground subroutine ended at: ',dtstr
    write(000,*) 'High Dispersal Establishment on New Bareground subroutine ended at: ',dtstr

    ! write intermediate coverage file for post high dispersal establishment on new bareground (1 = first time high dispersal is called)- also write summary output file
    call write_output('imo_hdest1',1)   !currently this is a 10*character string being passed in to write_output - currently needs to be padded with spaces
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Write Intermediate Output subroutine ended at: ',dtstr
    write(000,*) 'Write Intermediate Output subroutine ended at: ',dtstr 

    ! Apply mortality to non-flotant species and sum the total unoccupied land in each grid cell 
    !call sum_unoccupied_lnd
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Sum Unoccupied Land subroutine ended at: ',dtstr
    write(000,*) 'Sum Unoccupied Land ended at: ',dtstr

    ! Apply mortality to flotant species and sum the unoccupied flotant, keeping track of the different types (dead thin, dead thick, and bareground flotant)
    !call sum_unoccupied_flt
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Sum Unoccupied Flotant subroutine ended at: ',dtstr
    write(000,*) 'Sum Unoccupied Flotant ended at: ',dtstr

    ! write intermediate coverage file for post mortality - also write summary output file
    call write_output('imo_mort  ',1)   !currently this is a 10*character string being passed in to write_output - currently needs to be padded with spaces
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Write Intermediate Output subroutine ended at: ',dtstr
    write(000,*) 'Write Intermediate Output subroutine ended at: ',dtstr 

    ! Calculate the dispersal coverage for each species 
    !call calc_dispersal_coverage
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Calculate Dispersal Coverage subroutine ended at: ',dtstr
    write(000,*) 'Calculate Dispersal Coverage subroutine ended at: ',dtstr

    ! Calculate expansion likelihood, which accounts for the species establishment probability and abdundance in the area 
    exp_lkd = 0.0                                                                   ! initialize array to zero before first used
    exp_lkd = establish_P * disp_cov                                                ! includes all species with a establishment probabilty (e.g., flotant)

    ! Update the coverages using the expansion liklihood to determine what establishes on the unoccupied land
    !call update_coverages
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Update Coverages subroutine ended at: ',dtstr
    write(000,*) 'Update Coverages subroutine ended at: ',dtstr


    ! Update flotant coverages using the expansion liklihood to determine what establishes on unoccupied flotant
    !call update_flotant
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Update Flotant subroutine ended at: ',dtstr
    write(000,*) 'Update Flotant subroutine ended at: ',dtstr

    ! write intermediate coverage file for post standard establishment - also write summary output file
    call write_output('imo_stest ',1)   !currently this is a 10*character string being passed in to write_output - currently needs to be padded with spaces
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Write Intermediate Output subroutine ended at: ',dtstr
    write(000,*) 'Write Intermediate Output subroutine ended at: ',dtstr 

    ! Allow high dispersal species (class three "weedy" species) to establish on any remaining bareground
    !call high_disp_stest
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'High Dispersal Establishment on Remaining Bareground subroutine ended at: ',dtstr
    write(000,*) 'High Dispersal Establishment on Remaining Bareground subroutine ended at: ',dtstr

    ! write intermediate coverage file for post standard establishment - also write summary output file
    call write_output('imo_hdest2',1)   !currently this is a 10*character string being passed in to write_output - currently needs to be padded with spaces
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Write Intermediate Output subroutine ended at: ',dtstr
    write(000,*) 'Write Intermediate Output subroutine ended at: ',dtstr 


    ! Apply coverage changes caused by acute salinity
    !call acute_salinity_lnd
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Acute Salinity Land subroutine ended at: ',dtstr
    write(000,*) 'Acute Salinity Land subroutine ended at: ',dtstr

    ! Apply coverage changes in flotant species caused by acute salinity 
    !call acute_salinity_flt
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Acute Salinity Flotant subroutine ended at: ',dtstr
    write(000,*) 'Acute Salinity Flotant subroutine ended at: ',dtstr

    ! Round all coverages to the nearest 1 m^2 and check if anything is negative
    !call round_coverages
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Round Coverages subroutine ended at: ',dtstr
    write(000,*) 'Round Coverages subroutine ended at: ',dtstr

    ! Check that the sum of all coverages in each cell is 1.0 +/- given tolerance 
    !call check_sum
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Check Sum subroutine ended at: ',dtstr
    write(000,*) 'Check Sum subroutine ended at: ',dtstr

    ! Note, check sum needs to be called before adding the whole pixel portion of dead flotant to water so it is not double counted in dead flotant and water

    ! Add the whole Morph pixel portion of the dead flotant to water
    coverages(:,wti) = coverages(:,wti) + (dem_pixel_proportion*floor(coverages(:,dfi)/dem_pixel_proportion))

    
    ! write final coverage file for End of Year landscape - also write summary output file
    call write_output('eoy       ',1)   !currently this is a 10*character string being passed in to write_output - currently needs to be padded with spaces
    call date_and_time(VALUES=dtvalues)
    write(dtstr,8889) dtvalues(1),'-',dtvalues(2),'-',dtvalues(3),' ',dtvalues(5),':',dtvalues(6),':',dtvalues(7)
    write(  *,*) 'Write Output subroutine ended at: ',dtstr
    write(000,*) 'Write Output subroutine ended at: ',dtstr    


    write(  *,*)
    write(  *,*) 'Ended ICM-LAVegMod simulation at: ',dtstr
    write(  *,*)
    
    write(000,*)
    write(000,*) 'Ended ICM-LAVegMod simulation at: ',dtstr
    write(000,*)
    close(000)


8888    format(I4.4,I2.2,I2.2,a,I2.2,a,I2.2,a,I2.2)
8889    format(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2,a,I2.2)
    
end program
