module mod_PT_hMetal
use kdtree2_module ! this is for the kD-tree fixed radius search
implicit none

! global constants
integer, parameter          :: sp = kind(1.0), dp = kind(1.0d0)
double precision, parameter :: pi = 4.0d0 * atan(1.0d0), kappaD = 0.5d0,&
                               kappa = (1.0d0 - kappaD), kappaM = 0.5d0,&
                               kappaI = 1.0d0 - kappaM
    ! these kappas are the constants that scale how much of the total diffusion
    ! is simulated by random-walk vs. mass-transfer, and then by each "direction"
    ! of mass transfer
integer, parameter          :: nspec = 20, numsd = 3, num_alloc = 4e3
    ! NOTE: nspec (the number of chemical species in the model) is hard-coded for
        ! this specific problem. a more elegant solution would be better
    ! numsd is distances over which particle reactions are considered
    ! num_alloc is the maximum number of nearby particles expected to be found in the fixed radius search
        ! a larger number here avoids crashes and is better if you aren't expecting memory problems

! mobile particle type
! NOTE: these are only for a 1D problem
type mparticle
    double precision :: loc ! real-valued spatial location
    double precision :: concs(nspec) ! vector of chemical concentrations
    logical          :: active ! indicates whether particle is active and within the domain
end type

! immobile particle type
type iparticle
    double precision :: loc
    double precision :: concs(nspec)
end type

! a couple of derived types for the kD tree search
! holds indices of nearby particles
type index_array
    integer, allocatable :: indices(:)
end type
! holds the distances to the corresponding particle held by index_array
type dist_array
    double precision, allocatable :: dists(:)
end type

contains

! subroutine to initialize the random number generator seed from clock time
! from: https://gnu.huihoo.org/gcc/gcc-4.2.4/gfortran/RANDOM_005fSEED.html
subroutine init_random_seed()
    integer              :: i, n, clock
    integer, allocatable :: seed(:)

    call random_seed(size = n)
    allocate (seed(n))
    call system_clock(count = clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
end subroutine init_random_seed

! moves active particles via random-walk diffusion
subroutine diffuse(p, np, D, dt, alive)
    type(mparticle),  intent(inout) :: p(:) ! mobile particle array
    double precision, intent(in   ) :: D, dt ! diffusion coefficient and time step
    integer,          intent(in   ) :: np, alive(:) ! number and array of indices of active particles
    double precision                :: normvec(np) ! vector which will hold Normal(0, 1) values

    ! call N(0, 1) generator
    call box_mullerp(np, normvec)

    p(alive)%loc = p(alive)%loc + sqrt(2.0d0 * kappaD * D * dt) * normvec
end subroutine diffuse

! reflective lower boundary
subroutine reflectlow(p, low, alive)
    type(mparticle),  intent(inout) :: p(:) ! mobile particle array
    double precision, intent(in   ) :: low ! lower spatial boundary
    integer,          intent(in   ) :: alive(:) ! array of indices of active particles

    ! if particle has exited lower boundary, flip the negative location
    ! to the same positive value
    where (p(alive)%loc < low) p(alive)%loc = -p(alive)%loc
end subroutine reflectlow

! reflective upper boundary
subroutine reflecthigh(p, high, alive)
    type(mparticle),  intent(inout) :: p(:) ! mobile particle array
    double precision, intent(in   ) :: high ! upper spatial boundary
    integer,          intent(in   ) :: alive(:) ! array of indices of active particles

    ! if particle has exited upper boundary, reflect it back the same distance
    ! by which it overstepped
    where (p(alive)%loc > high) p(alive)%loc = 2.0d0 * high - p(alive)%loc
end subroutine reflecthigh

! since PHREEQCRM can't accept the particle array as input, this subroutine
! assigns the values in the 2D concs array to its corresponding immobile
! particle
subroutine concs_to_iparts(c, p, np)
    double precision, intent(in   ) :: c(:, :) ! 2D concentration array used by PHREEQCRM
    type(iparticle),  intent(inout) :: p(:) ! immobile particle array
    integer,          intent(in   ) :: np ! number of immobile particles
    integer                         :: i ! iteration variable

    do i = 1, np
        p(i)%concs = c(i, :)
    enddo
end subroutine concs_to_iparts

! this does the opposite of above
subroutine iparts_to_concs(c, p, np)
    double precision, intent(inout) :: c(:, :) ! 2D concentration array used by PHREEQCRM
    type(iparticle),  intent(in   ) :: p(:) ! immobile particle array
    integer,          intent(in   ) :: np ! number of immobile particles
    integer                         :: i ! iteration variable

    do i = 1, np
        c(i, :) = p(i)%concs
    enddo
end subroutine iparts_to_concs

! this subroutine builds the distance matrix for all mobile-immobile pairwise
! distances less than the cutoff radius
subroutine build_Distmat(ip, mp, na, alive, ni, D, dt, Distmat)
    type(iparticle),  intent(in   ) :: ip(:) ! immobile particle array
    type(mparticle),  intent(in   ) :: mp(:) ! mobile particle array
    integer,          intent(in   ) :: na, alive(:), ni
        ! number and array of indices of active mobile particles and number
        ! of immobile particles
    double precision, intent(in   ) :: D, dt ! diffusion coefficient and time step
    double precision, intent(  out) :: Distmat(ni, na) ! pairwise distance matrix
    type(kdtree2), pointer          :: tree ! this is the KD tree
    integer                         :: ntot, dim = 1
        ! total number of particles (active mobile + immobile)
        !****Note: hard coded one spatial dimension
    integer                         :: bindex, i, j ! indexing and looping variables
    real(kdkind)                    :: locs(na + ni), r2
        ! array holding locations of inactive and active immobile particles
        ! and value of squared search radius for KD search
    type(index_array), allocatable  :: closeguys(:) ! this holds the indices of nearby particles
    type(dist_array), allocatable   :: close_dists(:) ! this holds the distances to the corresponding nearby particle

    Distmat = 0.0d0

    ! calculate total number of particles to be considered for mass balance
    ntot = na + ni
    ! build locs array--immobile particles will be at the beginning of the
    ! array, and mobile will be at the end
    locs(1 : ni) = real(ip%loc, kdkind)
    locs(ni + 1 : ntot) = real(mp(alive)%loc, kdkind)
    ! calculate interaction distance to be numsd standard deviations of the
    ! Brownian Motion process--r2 is this distance squared
    ! ****NOTE: numsd is a global variable that is hard-coded above
    r2 = (real(numsd, kdkind) * sqrt(4.0_kdkind * real(D, kdkind) *&
                                     real(dt, kdkind)))**2

    ! build the KD tree and search it
    ! ****NOTE: num_alloc is a global variable that is hard-coded above
    call maketree(tree, dim, ntot, locs)

    allocate (closeguys(ni), close_dists(ni))
    ! this finds the closest mobile particles to each immobile particle
    call search(1, ni, tree, r2, num_alloc, closeguys, close_dists)
    ! NOTE: this search returns the SQUARED distance between two points
        ! also, the point itself is included in the closeguys list
    call kdtree2_destroy(tree)

    ! loop over immobile particles to build distance matrix
    do i = 1, ni ! immobile particle loop
        do j = 1, size(closeguys(i)%indices) ! mobile particle loop
            ! this is the mobile particle loop
            ! note that the closeguys array is indexed to the loc array,
            ! and thus its true index in the mp array is calculated below

            ! current mobile particle's index in locs array
            if (closeguys(i)%indices(j) <= ni) cycle
                ! if B is an immobile particle, skip this loop
                ! also prevents distance with self
            bindex = closeguys(i)%indices(j) - ni
            ! NOTE: Distmat is indexed to the mobile alive array
            Distmat(i, bindex) = close_dists(i)%dists(j)
        enddo
    enddo
    deallocate (closeguys, close_dists)
end subroutine build_Distmat

! algorithm set forth in mobile-immobile JCP paper using explicit matrix forward
! solve
subroutine immobile2mobile(ip, mp, na, alive, ni, D, dt, nc, Distmat)
    type(iparticle),  intent(inout) :: ip(:) ! immobile particle array
    type(mparticle),  intent(inout) :: mp(:) ! mobile particle array
    integer,          intent(in   ) :: na, alive(:), ni, nc
        ! number and array of indices of active mobile particles, number
        ! of immobile particles, and number of species in concs array
    double precision, intent(in   ) :: D, dt ! diffusion coefficient and time step
    double precision, intent(in   ) :: Distmat(ni, na) ! pairwise distance matrix
    integer                         :: i, j ! looping variables
    double precision                :: WImat(na, ni) ! mass-transfer matrix
    double precision                :: denom(na, ni) ! denom is used in colocation probability calculation but is pre-calculated for efficiency
    double precision                :: DTmat(na, ni) ! holds the transpose of the distance matrix
    double precision                :: colsum(ni) ! holds the column sums of the MT matrix

    DTmat = transpose(Distmat)

    ! calculate the denominator of the exponential bit of the colocation probability, scaled by the kappas
    do i = 1, ni
        denom(:, i) = -kappa * kappaI * 4.0d0 *  D * dt
    enddo

    ! zero out the very small values and exponentiate DTmat
    where (DTmat < 1.0d-16) WImat = 0.0d0
    where (DTmat /= 0.0d0) WImat = exp(DTmat / denom)

    ! iterate the normalization--10 times is an ad hoc guess, but seems to keep
        ! the error on the order of 10^-15
    do j = 1, 10
        colsum = sum(WImat, 1)

        do i = 1, ni
            if (colsum(i) /= 0.0d0) then
                WImat(:, i) = WImat(:, i) / colsum(i)
            endif
        enddo
    enddo

    ! conduct the mass transfers for all species
    do i = 1, nc
        mp(alive)%concs(i) = matmul(WImat, ip%concs(i))
    enddo
end subroutine immobile2mobile

! this does the opposite of the above
subroutine mobile2immobile(ip, mp, na, alive, ni, D, dt, nc, Distmat)
    type(iparticle),  intent(inout) :: ip(:) ! immobile particle array
    type(mparticle),  intent(inout) :: mp(:) ! mobile particle array
    integer,          intent(in   ) :: na, alive(:), ni, nc
        ! number and array of indices of active mobile particles, number
        ! of immobile particles, and number of species in concs array
    double precision, intent(in   ) :: D, dt ! diffusion coefficient and time step
    double precision, intent(in   ) :: Distmat(ni, na) ! pairwise distance matrix
    integer                         :: i, j ! looping variables
    double precision                :: WMmat(ni, na) ! mass-transfer matrix
    double precision                :: denom(ni, na) ! denom is used in colocation probability calculation but is pre-calculated for efficiency
    double precision                :: colsum(na) ! holds the column sums of the MT matrix

    ! calculate the denominator of the exponential bit of the colocation probability, scaled by the kappas
    do i = 1, ni
        denom(i, :) = -kappa * kappaM * 4.0d0 * D * dt
    enddo

    ! zero out the very small values and exponentiate DTmat
    where (Distmat < 1.0d-16) WMmat = 0.0d0
    where (Distmat /= 0.0d0) WMmat = exp(Distmat / denom)

    ! iterate the normalization--10 times is an ad hoc guess, but seems to keep
        ! the error on the order of 10^-15
    do j = 1, 10
        colsum = sum(WMmat, 1)
        do i = 1, na
            if (colsum(i) /= 0.0d0) then
                WMmat(:, i) = WMmat(:, i) / colsum(i)
            endif
        enddo
    enddo

    ! conduct the mass transfers for all species
    do i = 1, nc
        ip%concs(i) = matmul(WMmat, mp(alive)%concs(i))
    enddo
end subroutine mobile2immobile

! this builds a KD tree
subroutine maketree(tree2, d, n, locs)
    type(kdtree2), pointer, intent(  out) :: tree2 ! this is the KD tree
    integer,                intent(in   ) :: d, n ! number of spatial dimensions, number of particles
    real(kdkind),           intent(in   ) :: locs(d, n)
        ! location array for particles, with dimension d x n (number of
        ! spatial dimensions x number of particles)

    ! build the tree
    tree2 => kdtree2_create(locs, dim=d, sort=.false., rearrange=.true.)
        ! currently don't see a need to sort, as false is quicker, while
        ! rearrange = true is quicker
end subroutine maketree

! this searches an already built KD tree
subroutine search(start, end, tree, r2, num_alloc, closeguys, close_dists)
    integer,                intent(in   ) :: start, end ! start and end indices to search
    integer,                intent(in   ) :: num_alloc ! how large to to preallocate results array within KD tree module (hard-coded above)
    type(kdtree2), pointer, intent(in   ) :: tree ! the KD tree
    real(kdkind),           intent(in   ) :: r2 ! squared search radius
    type(index_array),      intent(  out) :: closeguys(:) ! this holds the indices of nearby particles
    type(dist_array),       intent(  out) :: close_dists(:) ! this holds the distances to the corresponding nearby particle
    integer                               :: i, n, nf ! loop iterator, number of particles in search, and number of particles found by search
    type(kdtree2_result), allocatable     :: results(:) ! results array from KD tree module

    allocate (results(num_alloc))
    n = end - start + 1

    ! loop over all particles
    do i = 1, n
        ! the type of search used here finds all the particles within
        ! squared distance r2 from the i^th particle in the list
        ! the hard-coded 0 is the 'correlation time' of the search
        call kdtree2_r_nearest_around_point(tree, i + start - 1, 0, r2, nf, num_alloc, results)

        ! allocate these based on how many nearby particles were found
        allocate (closeguys(i)%indices(nf), close_dists(i)%dists(nf))

        closeguys(i)%indices = results(1 : nf)%idx
        close_dists(i)%dists = results(1 : nf)%dis
    enddo

    deallocate (results)
end subroutine search

! these next two subroutines use the Box-Muller transform to generate N(0,1)
! random numbers from U(0,1)
! https://goo.gl/DQgmMu
! Note: this polar formulation seems to be consistently ~20% faster than the
! version below that uses trig functions
! reference for polar version (and standard version):
! https://www.taygeta.com/random/gaussian.html
subroutine box_mullerp(n, z)
    integer,          intent(in   ) :: n ! size of random array to be generated
    double precision, intent(  out) :: z(n) ! array of random variables

    ! ========================== LOCAL VARIABLES ===============================
    integer          :: j
    double precision :: w, x1, x2
    double precision :: rand(2)

    ! initialize the random seed, just in case
    call init_random_seed()

    do j = 1, n/2
        w = 1.0d0
        do while (w >= 1.0d0)
            call random_number(rand)
            x1 = 2.0d0 * rand(1) - 1.0d0
            x2 = 2.0d0 * rand(2) - 1.0d0
            w = x1**2 + x2**2
        enddo
        w = sqrt((-2.0d0 * log(w)) / w)
        z(2 * j - 1 : 2 * j) = (/x1 * w, x2 * w/)
    enddo

    if (mod(n, 2) /= 0) then
        w = 1.0d0
        do while (w >= 1.0d0)
            call random_number(rand)
            x1 = 2.0d0 * rand(1) - 1.0d0
            x2 = 2.0d0 * rand(2) - 1.0d0
            w = x1**2 + x2**2
        enddo
        w = sqrt((-2.0d0 * log(w)) / w)
        z(n) = x1 * w
    endif
end subroutine box_mullerp

! subroutine box_muller()
!     integer, parameter :: n = 1e8
!     integer            :: i, j
!     double precision   :: x1, x2, y1, y2, z(n)

!     call init_random_seed()
!     i = 1

!     do j = 1, n/2
!         x1 = rand()
!         x2 = rand()
!         y1 = sqrt(-2.0d0 * log(x1)) * cos(2.0d0 * pi * x2)
!         y2 = sqrt(-2.0d0 * log(x1)) * sin(2.0d0 * pi * x2)
!         z(i : i + 1) = (/y1, y2/)
!         i = i + 2
!     enddo
! end subroutine box_muller

! quicksort subroutine to sort a 1-dimensional array
! from: https://gist.github.com/t-nissie/479f0f16966925fa29ea#file-sortreal8-f-L4
recursive subroutine quicksort(a, first, last)
    implicit none
    double precision :: a(:), x, t
    integer          :: first, last
    integer          :: i, j

    x = a( (first + last) / 2 )
    i = first
    j = last
    do
        do while (a(i) < x)
            i=i+1
        enddo
        do while (x < a(j))
            j=j-1
        enddo
        if (i >= j) exit
        t = a(i);  a(i) = a(j);  a(j) = t
        i=i+1
        j=j-1
    enddo
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine

end module mod_PT_hMetal
