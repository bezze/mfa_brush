subroutine spring_array(mode)
#include "control_simulation.h"
#ifdef SPRING_ARRAY
    ! Routine that introduces an elastic (Hook) force between neighboring chains
    ! The force is applied to the second bead of the chain
    use commons

    INTEGER,  intent(in) :: mode

    INTEGER, SAVE :: cols,rows, bead_site
    INTEGER :: i,j,ni,nj,ir0,igraft,ineigh
    REAL (KIND=8) ::  r_neigh(3), r_graft(3), d, dist
    REAL (KIND=8), ALLOCATABLE :: f_array(:,:)

    select case (mode)

    case(0)

        cols= NINT(SQRT(0.5*n_chain*boundary(1)/boundary(2)))
        rows= NINT(SQRT(0.5*n_chain*boundary(2)/boundary(1)))
        k_spr_x= 80
        k_spr_y= 0
        f_array= 0.d0
        bead_site= 4 ! This is the number of the monomer in each chain, with grafted =0. It should be ranged
        ! between 1 and n_mon-1. If not, weird array springs will be formed.

        ALLOCATE(Mindex(0:rows+1,0:cols+1))

        Mindex=-1 ! Initializing to -1.

        ! This matrix should contain all indexes of the 1st bead of every chain
        do j=1,cols
            do i=1,rows
                Mindex(i,j) =  (j*rows-i)*n_mon + 1 ! +1 is the grafted bead
            end do
        end do

        ! This is the trick for periodic boundary conditions
        Mindex(:,0) = Mindex(:,cols)
        Mindex(:,cols+1) = Mindex(:,1)
        Mindex(0,:) = Mindex(rows,:)
        Mindex(rows+1,:) = Mindex(1,:)

        print *, ""
        print *, " * (Spring Array ON) Chains are linked by springs"
        print *, " -- k_spr_x=", k_spr_x 
        print *, " -- k_spr_y=", k_spr_y 

    case(1)
        ALLOCATE(f_array(3,n_part))
        f_array=0.d0
        v_array=0.d0

        do l=0,1 ! Bot/Top loop
            do j=1,cols
                do i=1,rows
                    igraft= Mindex(i,j)+l*n_chain/2 
                    ir0= igraft + bead_site

                    do nj=-1,1,2
                        ineigh= Mindex(i,j+nj)+ bead_site +l*n_chain/2
                        r_neigh = r0(:, ir0) - r0(:, ineigh)
                        ! r_neigh is the vector that joins the 2nd beads of chains

                        r_graft = r0(:, igraft) - r0(:, ineigh-bead_site) 
                        ! r_graft is the vector that joins the 1st beads of chains

                        d = SQRT(DOT_PRODUCT(r_graft,r_graft))
                        ! d es la magnitud de r_graft. Será la distancia de equil.

                        dist = SQRT(DOT_PRODUCT(r_neigh,r_neigh))
                        ! dist is the distance between the 2nd beads of neighbors

                        f_array(:, ir0) = f_array(:, ir0) - k_spr_x*(dist-d)*r_neigh/dist
                        v_array = v_array + .5*k_spr_x*(dist-d)**2
                        !print *, 'dist=',dist,'||  d=', d
                        !print *, k_spr_x*(dist-d)*r_neigh/dist
                    end do !nj

                    do ni=-1,1,2
                        ineigh= Mindex(i,j+ni)+ bead_site +l*n_chain/2
                        r_neigh = r0(:, ir0) - r0(:, ineigh)
                        ! r_neigh is the vector that joins the 2nd beads of chains

                        r_graft = r0(:, igraft) - r0(:, ineigh-bead_site) 
                        ! r_graft is the vector that joins the 1st beads of chains

                        d = SQRT(DOT_PRODUCT(r_graft,r_graft))
                        ! d es la magnitud de r_graft. Será la distancia de equil.

                        dist = SQRT(DOT_PRODUCT(r_neigh,r_neigh))
                        ! dist is the distance between the 2nd beads of neighbors

                        f_array(:, ir0) = f_array(:, ir0) - k_spr_y*(dist-d)*r_neigh/dist
                        v_array = v_array + .5*k_spr_y*(dist-d)**2
                        !print *, 'dist=',dist,'||  d=', d
                        !print *, k_spr_y*(dist-d)*r_neigh/dist
                    end do !ni
                end do !i
            end do !j
!print *, f_array
        end do !l
!print *, f_array(:,Mindex(2,1)),'|' ,f_array(:,Mindex(2,2)),'|',f_array(:,Mindex(2,3))
!do i=0,rows+1
!print *, i,')',Mindex(i,:)
!enddo
        force = force + f_array ! f_array is a vector FULL of zeros, may be expensive

        ! NEIGHBOURS:
        ! Each bead has 4 neigh. . According to gen_brush case(5) they're 
        ! generated first in rows and then in cols (row -> fast, col -> slow),
        ! starting from the lower left origin towards up and right.
        ! This means that each (closest) neighboring bead is of the form:
        ! -- Horizontal) ( i_chain +/- n_rows )*n_mon
        ! -- Vertical)   ( i_chain +/- 1 )*n_mon
    end select

#endif
end subroutine spring_array
