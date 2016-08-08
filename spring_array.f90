subroutine spring_array()
! Routine that introduces an elastic (Hook) force between neighboring chains
! The force is applied to the second bead of the chain

#include "control_simulation.h"
use commons

cols= NINT(SQRT(0.5*n_chain*boundary(1)/boundary(2)))
rows= NINT(SQRT(0.5*n_chain*boundary(2)/boundary(1)))
neig= (1, -1, rows, -rows)
f_array= 0.d0


do j_chain = 0, 1 !0 is top wall. 1 is for bottom wall
    l_half =  j_chain * n_mon * n_chain / 2 ! starting particle to set it's location
    do l = 0, n_chain-1
        l = l + l_half  ! top/bot correction
        
        ! NEIGHBOURS:
        ! Each bead has 4 neigh. . According to gen_brush case(5) they're 
        ! generated first in rows and then in cols (row -> fast, col -> slow),
        ! starting from the lower left origin towards up and right.
        ! This means that each (closest) neighboring bead is of the form:
        ! -- Horizontal) ( i_chain +/- n_rows )*n_mon
        ! -- Vertical)   ( i_chain +/- 1 )*n_mon
        
        do i = 1, 4
            
            r_neigh = r0(:, 2+l*n_mon) - r0(:, 2+(l+neigh(i))*n_mon)
            ! r_neigh is the vector that joins the 2nd beads of chains
            r_graft = r0(:, 1+l*n_mon) - r0(:, 1+(l+neigh(i))*n_mon)  
            ! r_graft is the vector that joins the 1st beads of chains
            d = SQRT(DOT_PRODUCT(r_graft,r_graft))
            ! d es la magnitud de r_graft. Será la distancia de equil.
            dist = SQRT(DOT_PRODUCT(r_neigh,r_neigh))
            ! dist is the distance between the 2nd beads of neighbors
            f_array(:,2+l*n_mon) = f_array(:,2+l*n_mon) - k_array*(dist-d)

        end do ! i

    end do ! l

end do ! j_chain
end subroutine spring_array
