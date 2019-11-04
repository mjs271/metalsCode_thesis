program FD_hMetal
use FD_hMetal_mod ! the module with all the subroutines
use PhreeqcRM ! calculates the geochemistry
implicit none

! ==============================================================================
!                           SIMULATION PARAMETERS
! ==============================================================================

integer, parameter          :: nthreads = 0
    ! number of OpenMP threads for reaction module. <= 0 implies the number of
    ! threads equals the number of processors of the computer
double precision, parameter :: maxtime = 1.577d8 ! 5 years = 1.577d8 seconds
double precision, parameter :: Omega = 0.4d0 ! length of domain [m]
double precision, parameter :: dx = 1d-2 ! spatial discretization parameter
integer, parameter          :: ncell = nint(Omega/dx)
    ! number of chemistry cells
    ! subtract 1 because won't be calculating chemistry for boundary cell
integer, parameter          :: ntrans = ncell + 1
    ! number of cells for transport
double precision, parameter :: dt = 25920.0_dp ! time step [s] (25920 seconds = 0.3 days)
integer, parameter          :: nsteps = nint(maxtime/dt) ! number of time steps

! ==============================================================================
!                           PHYSICAL PARAMETERS
! ==============================================================================

double precision, parameter        :: porosity = 0.47d0 ! porosity of medium
double precision, parameter        :: init_D = 4.27e-10 ! effective diffusion coefficient
double precision                   :: D(ntrans - 1) = init_D ! vector for diffusion coefficient

! ==============================================================================
!                       GENERAL/PHREEQCRM VARIABLES
! ==============================================================================

double precision, dimension(ncell) :: den_0, prs_0, tmp_0, sat_0, por_0, vol_0 ! vectors for passing initial parameters to PHREEQCRM
double precision                   :: cur_time ! current simulation time
integer                            :: i, j, m ! loop iterators
integer                            :: id, status ! id for the reaction module calls and status variable
integer                            :: ngrd, so_col ! number of grid points in the reaction module and number of columns in the selected output writes
integer                            :: ncomp, nchem ! number of aqueous components to transport and number of chemistry cells (get this from the RM for a quick error catching exercise)
integer                            :: ic1(ncell, 7) ! initial condition array
integer, dimension(ncell)          :: mask ! detailed chemistry printing mask for passing to RM
double precision, allocatable      :: bc_conc(:, :) ! boundary condition array
double precision, allocatable      :: comp_conc(:, :) ! concentrations array for RM
double precision, allocatable      :: sout(:, :) ! selected output array
double precision, allocatable      :: concs(:, :) ! concentrations array for transport calculations
integer                            :: bc1(1) ! this is for one boundary condition, to be passed to RM

! ==============================================================================
!                             READ/WRITE VARIABLES
! ==============================================================================

double precision, parameter   :: save_dt = dt * 1e3_dp
integer, parameter            :: save_steps = nint(maxtime/(save_dt)) + 1
    ! time step for saving concentrations for plotting and number of writes to be done
double precision, allocatable :: fig3_plot_concs(: , :), fig4_plot_concs(: , :)
    ! arrays for holding the data for generating the figure 3 and 4 plots
double precision, allocatable :: plot_times(:) ! array of times for writing out the plot data

! some more plotting-related variables
integer, parameter          :: fig3_numel = 6
integer                     :: fig3_elements(fig3_numel)
integer, parameter          :: fig3Plot_unit = 12
character(*), parameter     :: fig3Plot_name = 'fig3_concs.txt'

integer, parameter          :: fig4_numel = 3
integer                     :: fig4_elements(fig4_numel)
integer, parameter          :: fig4Plot_unit = 14
character(*), parameter     :: fig4Plot_name = 'fig4_concs.txt'

! ==============================================================================
! ==============================================================================

! initialize the simulation time to 0
cur_time = 0.0d0

! print some summary info, to begin with
print *, '======================================='
! print *, 'grid_Pe =', (dx * init_v) / init_D
print *, 'CFL = ', (init_D / porosity * dt) / dx**2
print *, '1.0d0 / dx = ', 1.0d0 / dx
print *, 'dt = ', dt
print *, 'ncell (chemistry) = ', ncell
print *, '======================================='

! ==============================================================================
!                   DEFINE PHYSICAL CONDITIONS (for PHREEQCRM)
! ==============================================================================
den_0 = 1.0d0 ! Density
prs_0 = 0.986923d0 ! Pressure 0.986923 atm = 1 bar
tmp_0 = 25.0d0 ! Temperature
sat_0 = 1.0d0 ! Saturation
por_0 = porosity ! porosity
vol_0 = 1.0d0 ! reaction cell volume

! print chemistry mask to print detailed output (headings) for only first element (this should speed things up marginally)
mask = 0
mask(1) = 1

! create the reaction module (RM)
id = RM_Create(ncell, nthreads)
! location of database file when running from main directory
status = RM_LoadDatabase(id, './phreeqc_files/PHREEQC-database.txt')
if (status < 0) then
    print *, 'Database Load Error'
    call exit(status)
endif

! initialize all the settings for the reaction module
! see the PhreeqcRM docs for more explanation: https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqcrm/classphreeqcrm.html
status = RM_OpenFiles(id)
status = RM_SetRepresentativeVolume(id, vol_0)
status = RM_SetSaturation(id, sat_0)
status = RM_SetPorosity(id, por_0)
status = RM_SetTime(id, cur_time)
status = RM_SetTimeStep(id, dt)
status = RM_SetDensity(id, den_0)
status = RM_SetTemperature(id, tmp_0)
status = RM_SetPressure(id, prs_0)
status = RM_SetSelectedOutputOn(id, 1) ! turn on selected output
status = RM_SetUnitsSolution(id, 2) ! 2 = mol/L
status = RM_SetUnitsPPassemblage(id, 0) ! 0 = mol/L of representative volume
status = RM_SetUnitsKinetics(id, 0) ! 0 = mol/L of representative volume
status = RM_SetUnitsExchange(id, 0) ! 0 = mol/L of representative volume
status = RM_SetUnitsSurface(id, 0) ! 0 = mol/L of representative volume
status = RM_SetUnitsSSassemblage(id, 1) ! 1 = mol/L of representative volume
status = RM_SetPrintChemistryMask(id, mask)
status = RM_SetPrintChemistryOn(id, 0, 0, 0) ! this toggles detailed output about reaction calculations
status = RM_SetScreenOn(id, 0) ! turn off messages about rebalancing
status = RM_SetComponentH2O(id, 1)
    ! 0 implies that H2O is not given as a separate component--just H and O
    ! 1 is default with H2O as separate component
! location of phreeqc input file when running from main directory
status = RM_RunFile(id, 1, 1, 1, './phreeqc_files/PHREEQC-Case1-Input.txt')

! the following isn't strictly necessary if you already know what your problem
! looks like, but it's quick and only done once.

ncomp = RM_FindComponents(id)
! get number of grid cells in model (this can change, so it needs to be checked)
ngrd = RM_GetGridCellCount(id)
if (ngrd /= ncell) then
    print *, '****** Cell count differs from particle number ******'
    call exit(status)
endif

nchem = RM_GetChemistryCellCount(id)
if (nchem /= ncell) then
    print *, '****** Chem count differs from particle number ******'
    call exit(status)
endif

! ==============================================================================
!                       INITIAL/BOUNDARY CONDITIONS
! ==============================================================================
! ===  (1) SOLUTIONS, (2) EQUILIBRIUM_PHASES, (3) EXCHANGE,
! ===  (4) SURFACE, (5) GAS_PHASE, (6) SOLID_SOLUTIONS, and (7) KINETICS

ic1 = -1
ic1(:, 1) = 1
ic1(:, 2) = 1
ic1(:, 4) = 1
ic1(:, 7) = 1

allocate(bc_conc(1, ncomp)) ! must be 2D array for module--requires singleton dimension in this case
bc1 = 0 ! corresponds to solution zero in input file

status = RM_InitialPhreeqc2Module(id, ic1)
status = RM_RunCells(id)

! get aqueous ion concentrations from the reaction module
allocate(comp_conc(ncell, ncomp))
status = RM_GetConcentrations(id, comp_conc)
status = RM_InitialPhreeqc2Concentrations(id, bc_conc, 1, bc1)

! ==============================================================================
!               SELECTED OUTPUT (based on the PHREEQC input file)
! ==============================================================================

! the names of the fields in sOut can be found in the file: phreeqc_files/sOut_fields.txt

! assign the indices of sOut that will be used for re-creating Fig. 3 of Arora, et al.
! These are: pH, Alkalinity, N(5), Sulfate, Feii, Pb
fig3_elements = (/ 7, 50, 10, 14, 12, 28 /)

! assign the indices of sOut that will be used for re-creating Fig. 4 of Arora, et al.
! These are: Ferrihydrite, FeS, and Siderite
fig4_elements = (/ 21, 22, 23  /)

! ==============================================================================
!        CONCS (what is being transported and given to the reaction module)
! ==============================================================================
! dimension of concs is ntrans (# chemistry cells + 1 for boundary) x ncomp

! For this specific problem:

! If (RM_SetComponentH2O(id, 0)):
! (1) H, (2) O, (3) Charge, (4) Acetate, (5) C, (6) Ca, (7) Cl, (8) Fe, (9) Feii,
! (10) K, (11) Mg, (12) N, (13) Na, (14) Naq, (15) Oxy, (16) Pb, (17) S,
! (18) Sulfate, (19) Zn

! If (RM_SetComponentH2O(id, 1)):
! (1) H2O, (2) H, (3) O, (4) Charge, (5) Acetate, (6) C, (7) Ca, (8) Cl, (9) Fe,
! (10) Feii, (11) K, (12) Mg, (13) N, (14) Na, (15) Naq, (16) Oxy, (17) Pb,
! (18) S, (19) Sulfate, (20) Zn

! ==============================================================================

! get the number of columns in the selected output, allocate the array, and get the output
so_col = RM_GetSelectedOutputcolumncount(id)
allocate(sout(ncell, so_col))
status = RM_GetSelectedOutput(id, sout)

! allocate the concs and plotting arrays
allocate(concs(ntrans, ncomp), fig3_plot_concs(ncell, fig3_numel),&
               fig4_plot_concs(ncell, fig4_numel), plot_times(save_steps))

! initialize the concentrations array used for transport calculations
concs(1, :) = bc_conc(1, :)
concs(2 : ntrans, :) = comp_conc

! get the specific quantities needed for Fig. 3 and 4 from the full sOut
fig3_plot_concs(:, 1 : fig3_numel) = sOut(:, fig3_elements)
fig4_plot_concs(:, :) = sOut(:, fig4_elements)

! array used later for deciding whether to get sOut and write the plot information
plot_times = (/((i - 1) * save_dt, i = 1, save_steps)/)

! this is the file that stores what comes from fig3_plot_concs
! the header gives the shape of the full array and the number of time steps
open (unit = fig3Plot_unit, file = fig3Plot_name, action = 'write')
write (fig3Plot_unit, *) shape(fig3_plot_concs), save_steps
write (fig3Plot_unit, *) fig3_plot_concs
close (unit = fig3Plot_unit, status = 'keep')
open (unit = fig4Plot_unit, file = fig4Plot_name, action = 'write')
write (fig4Plot_unit, *) shape(fig4_plot_concs), save_steps
write (fig4Plot_unit, *) fig4_plot_concs
close (unit = fig4Plot_unit, status = 'keep')

! counter for plot_times
j = 2

! time stepping loop
do m = 1, nsteps

    print *, 'step = ', m, ' of ', nsteps

    ! conduct the diffusion calculations
    call diffuse(concs(:, :), D, dx, dt, porosity, ntrans)

    ! update simulation time and send it to the RM
    cur_time = cur_time + dt
    status = RM_SetTime(id, cur_time)

    ! send the concs to the RM, run the geochem, and get it back
    comp_conc = concs(2 : ntrans, :)
    status = RM_SetConcentrations(id, comp_conc)
    status = RM_RunCells(id)
    status = RM_GetConcentrations(id, comp_conc)
    concs(2 : ntrans, :) = comp_conc

    ! if it's time, write out the plotting data
    if (cur_time >= plot_times(j)) then
        status = RM_GetSelectedOutput(id, sout)

        fig3_plot_concs(:, 1 : fig3_numel) = sOut(:, fig3_elements)
        fig4_plot_concs(:, :) = sOut(:, fig4_elements)

        open (unit = fig3Plot_unit, file = fig3Plot_name, action = 'write',&
              status = 'old', access = 'append')
        write (fig3Plot_unit, *) fig3_plot_concs
        close (unit = fig3Plot_unit, status = 'keep')
        open (unit = fig4Plot_unit, file = fig4Plot_name, action = 'write',&
              status = 'old', access = 'append')
        write (fig4Plot_unit, *) fig4_plot_concs
        close (unit = fig4Plot_unit, status = 'keep')
        j = j + 1
    endif

enddo

! destroy the reaction module
status = RM_Destroy(id)
deallocate(concs, fig3_plot_concs, fig4_plot_concs, plot_times, sout, comp_conc, bc_conc)

end program FD_hMetal
