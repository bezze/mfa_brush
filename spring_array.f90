subroutine spring_array()
! Routine that introduces an elastic (Hook) force between neighboring chains
! The force is applied to the second bead of the chain

#include "control_simulation.h"
use commons

cols=NINT(SQRT(0.5*n_chain*boundary(1)/boundary(2)))
rows=NINT(SQRT(0.5*n_chain*boundary(2)/boundary(1)))

do j_chain = 0, 1 !0 is top wall. 1 is for bottom wall
    l_half =  j_chain * n_mon * n_chain / 2 ! starting particle to set it's location
    do l = 0, n_chain-1
        l = l + l_half  ! top/bot correction

        ! esto no me gusta nada, encontrar una mejor forma de hacer el loop

        r_graft = r0(:, 1+(l-1)*n_mon)- r0(:, 1+(l-1)*n_mon)  
        ! r_graft es el vector que une los bead #1
        d = SQRT(DOT_PRODUCT(r_graft,r_graft))
        ! d es la magnitud de r_graft. Será la distancia de equil.



    end do

end do
end subroutine spring_array
