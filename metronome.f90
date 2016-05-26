! Routine that introduces an external driving periodical force on the grafted chains
subroutine metronome(mode_metro)
#include "control_simulation.h"
use commons
implicit none
integer, intent(in) :: mode_metro
real(kind=8) :: r_last(3), r_last_mod, r_2rel(3) , beta=85., cos_lim, sin_lim, cos_th, f_ext
real(kind=8), SAVE, ALLOCATABLE ::  cos_mem(:)
integer :: l, kick(n_chain) 
integer, SAVE :: count_local 
select case (mode_metro)

case(0)  ! Init  variables 
print*,""
print*," * Simulation with active polymer"
print*,""
ALLOCATE(cos_mem(n_chain) )
cos_mem = 0.
count_local = 0
!case(1)  ! Memory storing. Saves 2nd bead position at time t, for comparing with t+dt 

case(1)  ! Apply forces where angle condition requires
! We'll compare 'cos(theta)' rather than 'theta'
! It's biyective in [0,pi]->[-1,1], and cheaper to calclulate
!cos_lim = abs( cos(pi/2.- beta*pi/180.) ) 
cos_lim = cos(beta*pi/180.)
sin_lim = sin(beta*pi/180.)
f_ext = 10000. 
kick = 0

do l = 1, n_chain
    !r_last(:) = r0(:, l*n_mon)-r0(:, 1+(l-1)*n_mon)
    !r_last_mod = sqrt(dot_product(r_last,r_last))
    !cos_th = r_last(1)/r_last_mod        
    r_2rel = r0(:, 3+(l-1)*n_mon)- r0(:, 1+(l-1)*n_mon)         !REDO!!!
    cos_th = r_2rel(1)/SQRT(DOT_PRODUCT(r_2rel,r_2rel))
    if( (cos_mem(l).le.cos_lim ).and.(cos_th.ge.cos_lim) ) then 
        kick(l) = 15
        if(l.le.n_chain/2) then  ! top wall 
            force(1, 2+(l-1)*n_mon) = force(1, 2+(l-1)*n_mon) + f_ext*(sin_lim)!*sign(1.,cos_th)
            force(1, 3+(l-1)*n_mon) = force(1, 3+(l-1)*n_mon) + f_ext*(sin_lim)!*sign(1.,cos_th)
            force(3, 2+(l-1)*n_mon) = force(3, 2+(l-1)*n_mon) + f_ext*(cos_lim)!*sign(1.,cos_th)
            force(3, 3+(l-1)*n_mon) = force(3, 3+(l-1)*n_mon) + f_ext*(cos_lim)!*sign(1.,cos_th)
        else if (l.gt.n_chain/2) then !bottom wall
            force(1, 2+(l-1)*n_mon) = force(1, 2+(l-1)*n_mon) + f_ext*(sin_lim)!*sign(1.,cos_th)
            force(1, 3+(l-1)*n_mon) = force(1, 3+(l-1)*n_mon) + f_ext*(sin_lim)!*sign(1.,cos_th)
            force(3, 2+(l-1)*n_mon) = force(3, 2+(l-1)*n_mon) - f_ext*(cos_lim)!*sign(1.,cos_th)
            force(3, 3+(l-1)*n_mon) = force(3, 3+(l-1)*n_mon) - f_ext*(cos_lim)!*sign(1.,cos_th)
        end if

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
End do ! end do l

count_local = count_local + 1

if(count_local .eq.10) then
    count_local = 0
    write(79,*) i_time, kick
end if

end select
end subroutine metronome
