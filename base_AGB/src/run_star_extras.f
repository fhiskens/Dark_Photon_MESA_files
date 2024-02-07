!***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib

      implicit none

    real(dp) :: original_diffusion_dt_limit
    real(dp) :: postAGB_check = 0.0
    real(dp) :: rot_set_check = 0.0
    logical :: wd_diffusion = .false.
    real(dp) :: X_C_init, X_N_init

    integer :: time0, time1, clock_rate


       ! Parameters for the inclusion of energy-loss to the production of dark photons (MJD, FJH, RRV)
    real(dp) :: mixing_chi, m_dp, log10mDP
    ! Arrays for importing and evaluating data:
    real(dp) :: T_array(82) 
    real(dp) :: wpl_arrays(11, 22) !# of different wpl arrays, dim of each wpl array
    real(dp) :: gam_arrays(11, 7) !# of gam arrays, dim of each temp_array
    real(dp) :: table_data(82, 22, 7, 11, 11) !#Temperatures, wpl dim of tables, gam dim of tables, 10x10 tables
    real(dp) :: internal_workspace(82, 22, 7, 3, 11, 11) ! Additional 3 for internal processing
    real(dp) :: wpl_boundary_array(23)
    real(dp) :: gam_boundary_array(8)
    
    ! Physical constants
    real(dP) :: hbar_FH, c_FH, kb_FH, eVtoerg_FH, const1, const2, alpha_FH, me_FH, u_FH
                
             
    integer :: bsm_load, table_load, Ngam, Nwpl, NT, NTarray(1), HB_check, wpl_check, RGB_check, &
                early_check
    character(10) :: mass_dir
    character(256) :: path_to_file    
    
    
      ! these routines are called by the standard run_star check_model
      contains



            subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).
         s% other_neu => other_neu
         s% other_mesh_delta_coeff_factor => other_mesh_delta_coeff_factor


         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      
      
      ! Function to find the correct temperature indices to use.
          subroutine find_nearest_Temps(log10_T, T_array, T_indices)
         real(dp), intent(in) :: T_array(82), log10_T
         integer, intent(inout) :: T_indices(2)
         integer :: Tind(1)
         
         Tind = minloc(abs(log10_T-T_array))
         T_indices(1) = Tind(1) ! Finding the closest table to the temperature of the cell

                     ! Finding the second-closest table to the temperature of the cell
         if (log10_T-T_array(T_indices(1))>0) then
             T_indices(2) = T_indices(1) + 1
         else if (log10_T - T_array(T_indices(1)) < 0) then
             T_indices(2) = T_indices(1) - 1
             if (T_indices(2) > 82) then
                 write(*,*) "Tried to access table 83 or higher", log10_T
             end if
         end if
    end subroutine
    
    
    ! Function to find the correct plasma frequency indices to use.
    subroutine find_which_wpl_index(wpl, wpl_boundary_array, wpl_index)
        real(dp), intent(in) :: wpl, wpl_boundary_array(23)
        integer, intent(inout) :: wpl_index
        integer :: wpl_ind(1)        
        
        wpl_ind = minloc(abs(wpl-wpl_boundary_array))
        if (wpl-wpl_boundary_array(wpl_ind(1))>0) then
            wpl_index = wpl_ind(1)
        else if (wpl-wpl_boundary_array(wpl_ind(1))<=0) then
            wpl_index = wpl_ind(1)-1
        end if
    end subroutine
    
    ! Function to find the correct gam indices to use.
    subroutine find_which_gam_index(gam, gam_boundary_array, gam_index)
        real(dp), intent(in) :: gam, gam_boundary_array(8)
        integer, intent(inout) :: gam_index
        integer :: gam_ind(1)
        
        gam_ind = minloc(abs(gam-gam_boundary_array))
        if (gam-gam_boundary_array(gam_ind(1))>0) then
            gam_index = gam_ind(1)
        else if (gam-gam_boundary_array(gam_ind(1))<0) then
            gam_index = gam_ind(1)-1
        end if
        
        if (gam_index < 1) then
            gam_index = 1
            write(*,*) 'Gamma too low, defaulting to index 1', gam
        end if
    
    end subroutine
    
    ! Function for evaluating G in 'evaluate_energy_loss' subroutine.
    subroutine evaluate_G_function(Tind, wpl_ind, gam_ind, table_data, internal_workspace, gam_arrays, &
                                    wpl_arrays, gam_value, wpl_value, G_value, ierr)
        use interp_2d_lib_db, only: interp_rgbi3p_db
        real(dp), intent(in) :: table_data(82, 22, 7, 11, 11)
        real(dp), intent(in) :: internal_workspace(82, 22, 7, 3, 11, 11)
        real(dp), intent(in) :: gam_arrays(11, 7)
        real(dp), intent(in) :: wpl_arrays(11, 22)
        real(dp), intent(in) :: gam_value(1), wpl_value(1)
        integer, intent(in) :: Tind, wpl_ind, gam_ind
        real(dp), intent(inout) :: G_value(1)
        integer, intent(out) :: ierr
        
        !Internal parameters
        real(dp) :: table(11,11), wk_array(3, 11, 11), gam_array(11), wpl_array(11)
        integer :: i
        
        table = table_data(Tind, wpl_ind, gam_ind, :, :)
        !write(*,*) "Check table"
        !write(*,*) table
        wk_array = internal_workspace(Tind, wpl_ind, gam_ind, :, :, :)
        gam_array = gam_arrays(:, gam_ind)
        !write(*,*) "Check gam array"
        !write(*,*) gam_array
        wpl_array = wpl_arrays(:, wpl_ind)
        !write(*,*) "Check wpl array"
        !write(*,*) wpl_array
        
        !write(*,*) "Check work array"
        !do i=1,3
        !    write(*,*) wk_array(i, :, :)
        !end do
        
        call interp_rgbi3p_db(2, 11, 11, gam_array, wpl_array, table, &
        1, gam_value, wpl_value, G_value, ierr, wk_array)
        
    end subroutine
    
    
    ! Function for interpolating between temperature indices to find the final value of G.
    subroutine find_G_final(Tind1, Tind2, log10_T, T_array, g1, g2, g_fin)
        real(dp), intent(in) :: T_array(82), log10_T
        real(dp), intent(in) :: g1(1), g2(1)
        integer, intent(in) :: Tind1, Tind2
        real(dp), intent(inout) :: g_fin
        
        g_fin = (log10_T - T_array(Tind1))/(T_array(Tind2)-T_array(Tind1))*(g2(1)-g1(1)) + g1(1)
        ! Take 10^(g_fin) and use this as our final value (our interpolation was in log10(g))
        g_fin = 10**(g_fin)
    end subroutine
    
    
    
    ! Function for evaluating the energy-loss to transverse and longitudinal dark photons.
    subroutine evaluate_energy_loss(log10_T, wpl1, gam1, T_array, gam_arrays, wpl_arrays, &
                                    gam_boundary_array, wpl_boundary_array, table_data, &
                                    internal_workspace, g_fin, ierr)
        real(dp), intent(in) :: log10_T, wpl1(1), gam1(1)
        real(dp), intent(in) :: T_array(82), gam_arrays(11,7), wpl_arrays(11,22)
        real(dp), intent(in) :: wpl_boundary_array(23), gam_boundary_array(8)
        real(dp), intent(in) :: table_data(82, 22, 7, 11, 11)
        real(dp), intent(in) :: internal_workspace(82, 22, 7, 3, 11, 11)
        real(dp), intent(inout) :: g_fin
        integer, intent(out) :: ierr
        
        
        integer :: T_indices(2), wpl_index_int, gam_index_int
        real(dp) :: g1(1), g2(1)
    
        !if (wpl_check /= 1) then
        !    find_which_wpl_index(wpl1(1), wpl_boundary_array, wpl_index_int)
        !    write(*,*) 
        !end if
        !wpl_check = 1
    
        call find_nearest_Temps(log10_T, T_array, T_indices)
        call find_which_wpl_index(wpl1(1), wpl_boundary_array, wpl_index_int)
        call find_which_gam_index(gam1(1), gam_boundary_array, gam_index_int)
        
        
        call evaluate_G_function(T_indices(1), wpl_index_int, gam_index_int, table_data, internal_workspace, gam_arrays, &
                                    wpl_arrays, gam1, wpl1, g1, ierr)
        call evaluate_G_function(T_indices(2), wpl_index_int, gam_index_int, table_data, internal_workspace, gam_arrays, &
                                    wpl_arrays, gam1, wpl1, g2, ierr)
        
        call find_G_final(T_indices(1), T_indices(2), log10_T, T_array, g1, g2, g_fin)
        
        !write(*,*) g1, g2, log10(g_fin), log10_T, T_array(T_indices(1)), T_array(T_indices(2))
        
     end subroutine

	! Modified neutrino energy-loss function. Includes additional energy-loss to dark photons.
      subroutine other_neu(  &
                  id, k, T, log10_T, Rho, log10_Rho, abar, zbar, z2bar, log10_Tlim, flags, loss, sources, ierr)
               use neu_lib, only: neu_get
               use neu_def
               use interp_2d_lib_db, only: interp_rgbi3p_db

               integer, intent(in) :: id ! id for star
               integer, intent(in) :: k ! cell number or 0 if not for a particular cell
               real(dp), intent(in) :: T ! temperature
               real(dp), intent(in) :: log10_T ! log10 of temperature
               real(dp), intent(in) :: Rho ! density
               real(dp), intent(in) :: log10_Rho ! log10 of density
               real(dp), intent(in) :: abar ! mean atomic weight
               real(dp), intent(in) :: zbar ! mean charge
               real(dp), intent(in) :: z2bar ! mean charge squared
               real(dp), intent(in) :: log10_Tlim
               logical, intent(inout) :: flags(num_neu_types) ! true if should include the type of loss
               real(dp), intent(inout) :: loss(num_neu_rvs) ! total from all sources
               real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
               integer, intent(out) :: ierr
               character(10) :: T_file, wpl_index, gam_index
             character(256) :: full_path, path_to_T_dir


            ! Additional parameters to help calculate additional energy-loss
             real(dp) :: ans(1), gam1(1), wpl1(1), ne, g1(1), g2(1), g_fin, gam(1), wpl(1), &
                           eps_T, eps_L, eps_total, log10_T_test, log10_T_internal, ye
      	     integer :: i, j, l, T_indices(2), wpl_index_int, gam_index_int
             type (star_info), pointer :: s

             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return

             ! Calls the standard neu_get routine
             call neu_get(  &
                T, log10_T, Rho, log10_Rho, abar, zbar, z2bar, log10_Tlim, flags, &
                loss, sources, ierr)
             if (ierr /= 0) return

             ! Step One: Reading in dark photon parameters from the inlist, then loading and setting up the relevant interpolation tables
             
             if (bsm_load /= 1) then
                 bsm_load = 1 ! To guarantee that this step only happens once.
                 
                 log10mDP = s% job% extras_rpar(4) ! log10(m_dp) - useful to have loaded
                 mixing_chi = 10**(s% job% extras_rpar(5)) ! The dark photon mixing parameter
                 m_dp = 10**(log10mDP) ! Dark photon mass
                 
                 write(mass_dir, '(f10.2)') log10mDP ! Convert log10mDP to a string
                 path_to_file = '/path_to_file/Dark_Photon_MESA_files/Data/m_' // trim(adjustl(mass_dir)) ! Defining the path to the relevant directory
                 !write(*,*) path_to_file
                 
                 open(2, file = trim(path_to_file) // "/file_list.txt", status = 'old', access = 'sequential', &
                                    form = 'formatted', action = 'read')
                 read(2,*) T_array
                 close(2)
                 
                 ! Load in information about the spacing of our grids in log10(wpl) and log10(gam1) as defined in the paper
                 open(3, file = trim(path_to_file) // "/wpl_arrays.txt", status = 'old', access = 'sequential', &
                                    form = 'formatted', action = 'read')
                 read(3,*) wpl_arrays !wpl_array.txt needs to be formatted as 20 rows with 10 columns
                 close(3)

                 open(4, file = trim(path_to_file) // "/gam_arrays.txt", status = 'old', access = 'sequential', &
                                     form = 'formatted', action = 'read')
                 read(4,*) gam_arrays !gam_array.txt needs to be formatted as 7 rows with 10 columns
                 close(4)
                 
                 !  These values are set - we have tested accuracy of interpolation across grids of this size and spaced as per the files
                 Ngam = 11 
                 Nwpl = 11
                 
                 
                 ! Now we read in the data to the first NT entries of the array import_data
                 do i=1, 82
                            write(T_file, '(f10.2)') T_array(i) ! Convert to string with 2dp
                            path_to_T_dir = trim(path_to_file) // '/T_' // trim(adjustl(T_file)) // '/'
                            !write(*,*) path_to_T_dir
                            
                            do j=1, 22
                                do l=1,7
                                    write(wpl_index, '(I2)') j !Convert wpl index to string
                                    write(gam_index, '(I2)') l !Same
                                    full_path = trim(path_to_T_dir) // trim(adjustl(wpl_index)) // '_' // &
                                     trim(adjustl(gam_index)) // '.txt'
                                     
                                    !write(*,*) full_path
                                    open(1, file = trim(full_path), status = 'old', access = 'sequential', &
                                           form = 'formatted', action = 'read')
                                    read(1, *) table_data(i, j, l, :, :)
                                    close(1)    
                                end do
                            end do
                 end do
                 
                 
                 ! Now we set up the arrays wk_array - a set of NT 3 x Ngam x Nwpl arrays which act as the internal workspace for
                 ! the interpolation routines. These cannot be initialised without also computing and returning an interpolated value
                 ! so we also specify a value of gam1 (x coord) and wpl1 (y coord) below. The computed value is stored in the global
                 ! paramter ans, and is never used.
                 gam1 = 3. 
                 wpl1 = 2.
                 
                 do i=1,82
                     do j=1,22
                         do l=1,7
                             !write(*,*) i, j, l
                             !write(*,*) gam_arrays(:,l)
                             !write(*,*) wpl_arrays(:,j)
                             !table_data(i,j,l,:,:)
                             call interp_rgbi3p_db(1, Ngam, Nwpl, gam_arrays(:, l), wpl_arrays(:, j), table_data(i, j, l, :, :), & 
                             1, gam1, wpl1, ans, ierr, internal_workspace(i, j, l, :, :, :))
                         end do
                     end do
                 end do
                 
                 ! Finally, we define a set of physical constants. This is performed once.
                 hbar_FH = 6.582d-16 ! hbar in eVs
                 c_FH = 3.0d10 ! c in cgs
                 kb_FH = 8.617d-5 ! Boltzmann's constant in eV/K
                 eVtoerg_FH = 1.60218d-12 ! conversion btw eV and erg
                 alpha_FH = 1.0/137.0 ! em fine structure constant
                 me_FH = 0.511d6 ! Electron mass in eV/c^2
                 u_FH = 1.66054d-24 ! atomic mass unit in g

                 const1 = hbar_FH*hbar_FH*hbar_FH*c_FH*c_FH*c_FH
                 const2 = eVtoerg_FH/(hbar_FH*hbar_FH*hbar_FH*hbar_FH*c_FH*c_FH*c_FH)
                 
                 ! To do: Write arrays of wpl and gam values at boundaries of tables
                 wpl_boundary_array(1) = wpl_arrays(1, 1)
                 do i=1, 22
                     wpl_boundary_array(i+1) =wpl_arrays(11, i)
                 end do
                 
                 gam_boundary_array(1) = gam_arrays(1, 1)
                 do i=1, 7
                     gam_boundary_array(i+1) = gam_arrays(11, i)
                 end do
                
                 write(*,*) "Dark Photon data has been initialised"
                 table_load = 1 ! Tells other threads that the tables and constants defined above have been loaded.
                 
             end if
             
             ! Step Two: Determine the magnitude of energy-loss to transverse dark photons
             ne = Rho*(s% ye(k))/u_FH ! Electron number density in cell
             ye = zbar/abar
             !wpl(1) = (sqrt(const1*4*pi*alpha_FH*ne/me_FH)) ! plasma frequency of cell
             wpl(1) = (3.33d5*8.617d-5)*sqrt((Rho)*(ye))/ &
                  (1+(1.019d-6*(Rho)*(ye))**(2./3.))**(1./4.)
             gam(1) = (const1*(4*pi*alpha_FH**2/(3*me_FH))*sqrt((2*pi*me_FH)/(3*kb_FH*T))*(wpl(1))**2*(z2bar/abar)*(Rho)/u_FH) ! value of gam1 in cell
             
             if (wpl(1) < 1.0d6) then
                 if (table_load == 1) then
             
                 ! Specify the appropriate values of log10_wpl and log10_gam - our interpolation variables
                 wpl1(1) = log10(wpl(1))
                 gam1(1) = log10(gam(1))
                 
                 log10_T_internal = log10_T
                 
                 if (log10_T >= 9.0) then
                     write(*,*) log10_T
                     log10_T_internal = 8.999
                     write(*,*) "log10_T is greater than 9.0"
                 end if
                 
                 if (gam(1) < 1.0d-12) then
                     gam(1) = -11.999
                 end if
              
                 
                 call evaluate_energy_loss(log10_T_internal, wpl1, gam1, T_array, gam_arrays, wpl_arrays, &
                                    gam_boundary_array, wpl_boundary_array, table_data, &
                                    internal_workspace, g_fin, ierr)
                 else
                     g_fin = 0. ! If tables haven't yet been loaded, we set g_fin as 0
                 end if
             else
                 !write(*,*) "plasma frequency greater than upper value"
                 g_fin = 0.
             end if
             
             !write(*,*) g1, g2, log10(g_fin)
             
             
             
             ! Set the value as 0, then if tables have been loaded we replace it with the actual value
             eps_T = 0
             eps_L = 0
             
             if (table_load == 1) then 
                 eps_T = (const2/(pi**2))*mixing_chi**2*m_dp**4*g_fin/Rho
                 
                 if (wpl(1)>m_dp) then
                     ! Longitudinal energy-loss (resonant) - from Ayala et al. 2020
                     eps_L = (1/kb_FH)**3*(1.62d4*mixing_chi**2*m_dp**2/Rho)*wpl(1)**2*sqrt(wpl(1)**2-m_dp**2) & 
                     /(exp(wpl(1)/(kb_FH*T))-1)
                 end if
             end if


             eps_total = eps_T + eps_L

             loss(ineu) = loss(ineu)+eps_total ! Benign if commented out
            
      end subroutine other_neu


      ! Function for reducing mesh_delta_coeff near dark photon RPR
      subroutine other_mesh_delta_coeff_factor(id, eps_h, eps_he, eps_z, ierr)
        !use const_def
        !use chem_def
        integer, intent(in) :: id
        real(dp), intent(in), dimension(:) :: eps_h, eps_he, eps_z
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        real(dp) :: wpl(50000), ne, c, hbar1, mH, alpha, me, coef
        integer :: k, nz, opt

        include 'formats'

        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        coef = 0.01


        opt = 1 !1 for basic

        c = 3.d10
        me = 0.511d6
        hbar1 = 6.582d-16
        mH = 1.66054d-24
        alpha = 1.0/137.0

        nz = s% nz

        do k=1, nz
                 ne = (s% rho(k))*(s% ye(k))/mH
                 !wpl(k) = sqrt(hbar1**3*c**3*4*3.14159265359*alpha*ne/me)
                 wpl(k) = (3.33d5*8.617d-5)*sqrt((s% rho(k))*(s% ye(k)))/ &
                  (1+(1.019d-6*(s% rho(k))*(s% ye(k)))**(2./3.))**(1./4.)
                 
                 if (abs(log10(wpl(k)) - s% job% extras_rpar(4)) < 0.01) then
                      s% mesh_delta_coeff_factor(k) = coef
                 end if
        end do



      end subroutine other_mesh_delta_coeff_factor
     


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
	 character(len=256) :: photosphere_summary, tau100_summary


         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

! set the correct summary file for the BC tables depending on [a/Fe]
        !photosphere_summary = 'table_' // trim(s% job% extras_cpar(2)) // 'summary.txt'
        !tau100_summary = 'table100_' // trim(s% job% extras_cpar(2)) // 'summary.txt'
        !call table_atm_init(.true., tau100_summary, photosphere_summary, ierr)





      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 7
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n), c, hbar1, me, mH, wpl(50000), alpha, ne
         integer, intent(out) :: ierr
         integer :: i, nz, k(1), j
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         c = 3.d10
         me = 0.511d6
         hbar1 = 6.582d-16
         mH = 1.66054d-24
         alpha = 1.0/137.0

         nz = s% nz

         do i=1, nz
                  wpl(i) = (3.33d5*8.617d-5)*sqrt((s% rho(i))*(s% ye(i)))/ &
                  (1+(1.019d-6*(s% rho(i))*(s% ye(i)))**(2./3.))**(1./4.)
         end do

         names(1) = 'max_he_ye'

         if (s%   max_eps_he_k > 0) then
                vals(1) = s% ye(s% max_eps_he_k)
         else
                vals(1) = 0
         end if


         names(2) = 'central_wpl'
         vals(2) = wpl(nz)

         k = minloc(abs(wpl-10**(s% job% extras_rpar(4))))
         j = k(1)

         names(3) = 'nearest_wpl'
         vals(3) = wpl(j)

         names(4) = 'nearest_wpl_m'
         vals(4) = s% m(j)

         names(5) = 'nearest_wpl_T'
         vals(5) = s% T(j)

         names(6) = 'nearest_wpl_Rho'
         vals(6) = s% Rho(j)
         
         names(7) = 'nearest_wpl_deg'
         vals(7) = s% eta(j)

         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.


      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
         
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n), c, me, hbar1, mH, alpha, wpl(50000), kb, const1, gam(50000), ne
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, i
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         c = 3.d10
         me = 0.511d6
         hbar1 = 6.582d-16
         mH = 1.66054d-24
         alpha = 1.0/137.0
         kb = 8.617d-5
         const1 = hbar1*hbar1*hbar1*c*c*c
         
         
         
         
         do i=1, nz
             wpl(i) = (3.33d5*8.617d-5)*sqrt((s% rho(i))*(s% ye(i)))/ &
             (1+(1.019d-6*(s% rho(i))*(s% ye(i)))**(2./3.))**(1./4.)
             !gam(i) = const1*(4*3.14159265359*alpha**2/(3*me))*sqrt((2*3.14159265359*me)/(3*kb*(s% T(i))))* &
             !(wpl(i))**2*(s% z2bar(i))*(s% Rho(i))/mH
         end do
         
         
         

         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
         names(1) = 'plasma_frequency'
         vals(:, 1) = wpl(1:nz)
         
         

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
	 real(dp) :: envelope_mass_fraction, L_He, L_tot, min_center_h1_for_diff, &
            critmass, feh, rot_full_off, rot_full_on, frac2, mass_difference
         real(dp), parameter :: huge_dt_limit = 3.15d16 ! ~1 Gyr
         real(dp), parameter :: new_varcontrol_target = 1d-3
         real(dp), parameter :: Zsol = 0.0142
         type (star_info), pointer :: s
	    logical :: diff_test1, diff_test2, diff_test3


         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going


! set ROTATION: extra param are set in inlist: star_job - this has been carried over from the MIST example run_star_extras files
        rot_full_off = s% job% extras_rpar(1) !1.2
        rot_full_on = s% job% extras_rpar(2) !1.8

        if (rot_set_check == 0) then
            if ((s% job% extras_rpar(3) > 0.0) .and. (s% initial_mass > rot_full_off)) then
                !check if ZAMS is achieved, then set rotation
                if ((abs(log10(s% power_h_burn * Lsun / s% L(1))) < 0.01) .and. (s% star_age > 1d2)) then
                    if (s% initial_mass <= rot_full_on) then
                        frac2 = (s% initial_mass - rot_full_off) / (rot_full_on - rot_full_off)
                        frac2 = 0.5d0*(1 - cos(pi*frac2))
                    else
                        frac2 = 1.0
                    end if
                    write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                    write(*,*) 'new omega_div_omega_crit, fraction', s% job% extras_rpar(3) * frac2, frac2
                    write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                    s% job% new_omega_div_omega_crit = s% job% extras_rpar(3) * frac2
                    s% job% set_near_zams_omega_div_omega_crit_steps = 10
                    rot_set_check = 1
                end if
            else
                write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                write(*,*) 'no rotation'
                write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                rot_set_check = 1
            end if

        end if

! End run at TP AGB
        mass_difference = s% he_core_mass - s% c_core_mass
        if (s% center_he4 < 1d-4) then
            if (mass_difference < 0.5*(s% he_core_mass)) then
                if  (s% luminosity_by_category(2, 1) > s% luminosity_by_category(3, 1)) then
                    write(*,*) 'Star at TP-AGB phase - stopping run'
                    call star_write_model(id, s% job% save_model_filename, ierr)
                    extras_finish_step = terminate
                end if
            end if
        end if

! check DIFFUSION: to determine whether or not diffusion should happen
! no diffusion for fully convective, post-MS, and mega-old models
! do diffusion during the WD phase
	    min_center_h1_for_diff = 1d-10
	    diff_test1 = abs(s% mass_conv_core - s% star_mass) < 1d-2 !fully convective
	    diff_test2 = s% star_age > 5d10 !really old
	    diff_test3 = s% center_h1 < min_center_h1_for_diff !past the main sequence
	    if( diff_test1 .or. diff_test2 .or. diff_test3 )then
            s% diffusion_dt_limit = huge_dt_limit
        else
            s% diffusion_dt_limit = original_diffusion_dt_limit
	    end if

        if (wd_diffusion) then
            s% diffusion_dt_limit = original_diffusion_dt_limit
        end if



         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve



      end module run_star_extras
