! Routine that introduces an external driving periodical force on the grafted chains
subroutine metronome(mode_metro)
#include "control_simulation.h"
use commons
implicit none
integer, intent(in) :: mode_metro
real(kind=8) :: r_last(3), r_last_mod, r_2rel(3) , beta=65., cos_lim, cos_th
real(kind=8), SAVE, ALLOCATABLE ::  cos_mem(:)
integer :: l

select case (mode_metro)

case(0)  ! Init  variables 
print*,""
print*," * Simulation with active polymer"
print*,""
ALLOCATE(cos_mem(n_chain) )
cos_mem = 0.
!case(1)  ! Memory storing. Saves 2nd bead position at time t, for comparing with t+dt 

case(1)  ! Apply forces where angle condition requires
! We'll compare 'cos(theta)' rather than 'theta'
! It's biyective in [0,pi]->[-1,1], and cheaper to calclulate
!cos_lim = abs( cos(pi/2.- beta*pi/180.) ) 
cos_lim = cos(beta*pi/180.)

Do l=1,n_chain !n_chain , loop over chains
        !r_last(:) = r0(:, l*n_mon)-r0(:, 1+(l-1)*n_mon)
        !r_last_mod = sqrt(dot_product(r_last,r_last))
        !cos_th = r_last(1)/r_last_mod        
        r_2rel = r0(:, 2+(l-1)*n_mon)- r0(:, 1+(l-1)*n_mon)         
        cos_th = r_2rel(1)/SQRT(DOT_PRODUCT(r_2rel,r_2rel))

        if( (cos_mem(l).le.cos_lim ).and.(cos_th.ge.cos_lim) ) then 
                force(1, 2+(l-1)*n_mon) = force(1, 2+(l-1)*n_mon) + sign(1.,cos_th)*50.
        end if

        cos_mem(l) = cos_th  ! Store new memory for next step

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
