module allocation

IMPLICIT NONE
public :: allocation_3pg

contains

    subroutine allocation_3pg(m0, fertility, pRx, pRn, f_phys, pFS, m, npp_fract_root, npp_fract_stem, npp_fract_foliage, n_sp) 
        ! calculate allocation parameters

        ! input
        integer:: n_sp
        real(kind=8), dimension(n_sp), intent(in) :: m0, fertility, pRx, pRn, f_phys, pFS

        ! output
        real(kind=8),dimension(n_sp), intent(out) ::  npp_fract_root, npp_fract_stem, npp_fract_foliage, m

        ! Get growing conditions modifier (assumes more C to roots in harsher growing conditions)
        m(:) = m0(:) + (1.d0 - m0(:)) * fertility(:)

        ! Get allocation fraction to roots, foliage and stems
        npp_fract_root(:) = pRx(:) * pRn(:) / (pRn(:) + (pRx(:) - pRn(:)) * f_phys(:) * m(:))
        npp_fract_stem(:) = (1.d0 - npp_fract_root(:)) / (1.d0 + pFS(:))
        npp_fract_foliage(:) = 1.d0 - npp_fract_root(:) - npp_fract_stem(:)

    end subroutine allocation_3pg

end module allocation