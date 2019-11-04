module FD_hMetal_mod
implicit none

! global constants
double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
integer, parameter          :: sp = kind(1.0), dp = kind(1.d0)
integer, parameter          :: nspec = 20 ! this is hard-coded for this problem. a more elegant solution would be better

contains

! diffusion subroutine using explicit Euler finite-differences
subroutine diffuse(mat, D, dx, dt, phi, ncell)
    double precision, intent(inout) :: mat(:, :)
    double precision, intent(in   ) :: D(:), dx, dt, phi
    integer,          intent(in   ) :: ncell
    double precision                :: Dmat(ncell - 2, nspec), temp(nspec)
    integer                         :: i

    ! stop concentrations from diffusing in across lower boundary (i.e., the final cell)
        ! by eliminating concentration gradient between cells end and end - 1,
        ! then replace after
    temp = mat(ncell, :)
    mat(ncell, :) = mat(ncell - 1, :)
    do i = 1, nspec
        Dmat(:, i) = D(2 : ncell - 1)
    enddo

    ! finite-difference evaluation
    mat(2 : ncell - 1, :) = mat(2 : ncell - 1, :) + ((Dmat * dt)/(phi * dx**2)) *&
                            (mat(3 : ncell, :) - 2 * mat(2 : ncell - 1, :)&
                                + mat(1 : ncell - 2, :))
end subroutine diffuse

end module FD_hMetal_mod
