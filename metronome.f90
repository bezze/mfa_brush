! Routine that introduces an external driving periodical force on the grafted chains
subroutine metronome(mode_metro)
#include "control_simulation.h"
use commons
implicit none
integer, intent(in) :: mode_metro
real(kind=8) :: r_last(3), r_last_mod, beta=25., cos_lim, cos_th
integer :: l

select case (mode_metro)

case(0)  ! Init  variables 
print*,""
print*," * Simulation with active polymer"
print*,""

case(1)

cos_lim = abs( cos(3.1416/2.- beta*3.1416/180.) ) ! tidy up please

Do l=1,n_chain !n_chain , loop over chains
        r_last(:) = r0(:, l*n_mon)-r0(:, 1+(l-1)*n_mon)
        r_last_mod = sqrt(dot_product(r_last,r_last))
        cos_th = r_last(1)/r_last_mod
        if( cos_th.le.cos_lim) then 
                force(1, 2+(l-1)*n_mon) = force(2, 2+(l-1)*n_mon) + sign(1.,cos_th)*50.
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  TEST TO CHECK IF THE ROUTINE IS CALCULATING FORCES AND ENERGY
        !  CORRECTLY. 
        !print*, r0, "positions"
        !print*, k_bend, "bending constant"
        !print*, F_bend, "bending FORCE"
        !print*, v_bend, "bending ENERGY"
        ! ERASE AFTER CHECKING
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
End do
end select
end subroutine metronome
