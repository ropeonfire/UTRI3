!! UTRI3.f90
!!
!! Copyright (c) 2016-2018 William Matthew Peterson
!! This code is distributed under the LGPLv3 license.
!! 
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!! Element Type:
!! U21: plane strain
!! U25: plane stress
!! 
!! Element Properties:
!! prop(1) = E1
!! prop(2) = E2
!! prop(3) = nu12
!! prop(4) = G12
!! prop(5) = density
!! prop(6) = thickness
!! iprop(1) = matbehav  (1:ISOTROPIC, 2:LAMINA)
!! 
!! Input File Usage:
!! *user element, type=U2X, nodes=3, coordinates=2, properties=6, iproperties=1
!!  1,2
!! *user property, elset=utri3_elset
!!  <E1>, <E2>, <nu12>, <G12>, <rho>, <thk>, <matbehav>
!! *element, type=U2X, elset=utri3_elset
!!  elem#, n1,n2,n3
!!
!! Subroutine Usage:
!! At the beginning of the user source code (e.g. "usub.for"), include the UTRI3.F90 module:
!!  include 'utri3.f90'
!! 
!! Within the UEL subroutine, USE the UTRI3 module and call the appropriate module subroutine:
!!  UEL_code:   use utri3
!!  UEL_code:   
!!  UEL_code:   if ((jtype==21).or.(jtype==25)) then
!!  UEL_code:       call tri3(rhs, amatrx, svars, ndofel, nsvars, props, nprops, coords, mcrd, &
!!  UEL_code:           nnodes, u, jtype, time, dtime, kstep, kinc, jelem, predef, npredf,     &
!!  UEL_code:           pnewdt, jprops, njprops, nthreads, ithread)
!!  UEL_code:   endif
!! 
!! Note 1: nthreads and ithread are not currently used and may be left undefined.
!! Note 2: If fixed-form is used in the main subroutine source-code, you must use the correct
!!         indentation levels, line lengths, and line continuation markers for the 'UEL_code:'
!!         in the example section above.
!! 
!! Analysis Usage (command line):
!!    abaqus job=Job-1 input=input.inp user=usub.for
!! 
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!DIR$ FREEFORM

module utri3
    !! Standard 2D 3-node linear triangle User Element for Abaqus/Standard.
    
    use uparams
    implicit none
    
    private
    public :: tri3
    
    type :: tri3_type
        integer        :: nnodes = NELEMNODES_TRI3                              !! No. of nodes on element
        integer        :: nnface = NFACENODES_TRI3                              !! No. of nodes on face of the element adjacent to the interface
        integer        :: ngp    = NINTEGRPTS_TRI3                              !! No. of gauss integation points
        integer        :: ndof   = NDOF_TRI3                                    !! No. of degrees of freedom
        integer        :: elemtype                                              !! Integer type encoding (see also clipped6.externaldb.getIntegerTypes)
        integer        :: label                                                 !! Integer Part-level element label
        integer        :: matbehav                                              !! [ISOTROPIC | LAMINA | ENGINEERING_CONSTANTS]
        integer        :: pstype                                                !! 2D plane strain=1, 2D plane stress=2
        real(8)        :: thick                                                 !! 2D=thickness
        real(8)        :: density                                               !! mass density (M/L^3)
        integer, dimension(NELEMNODES_TRI3)          :: nodes                   !! nodal connectivity list
        integer, dimension(NFACENODES_TRI3)          :: facenodes               !! list of nodes on interface
        real(8), dimension(NPROPS_2D)                :: props                   !! mechanical properties
        real(8), dimension(DIMS_2D, NELEMNODES_TRI3) :: curcrds                 !! current coordinates for each node on element
        real(8), dimension(NDOF_TRI3)                :: u                       !! nodal displacement vector
        real(8), dimension(NDOF_TRI3)                :: rhs                     !! nodal force vector
        real(8), dimension(STRS_2D, NFACENODES_TRI3) :: exstrs                  !! stresses computed at or extrapolated to interface gauss points
        real(8), dimension(STRS_2D, NINTEGRPTS_TRI3) :: gpstrs                  !! stresses calculated at internal gauss points
        real(8), dimension(NDOF_TRI3, NDOF_TRI3)     :: stiff                   !! stiffness matrix
        real(8), dimension(NINTEGRPTS_TRI3)          :: gpwgt                   !! numerical weights for internal gauss points
        real(8), dimension(DIMS_2D, NINTEGRPTS_TRI3) :: gpcrds                  !! natural (r,s) coords for internal gauss points
        contains
        procedure, pass :: initElement
        procedure, pass :: computeCurrentCoords
        procedure, pass :: computeCmat
        procedure, pass :: integrateElement
        procedure, pass :: assembleElement
    end type tri3_type
    
    
contains
    subroutine tri3(rhs, amatrx, svars, ndofel, nsvars, props, nprops, coords,   &
        mcrd, nnodes, u, jtype, time, dtime, kstep, kinc, jelem, predef, npredf, &
        pnewdt, jprops, njprops, nthreads, ithread)
        
        implicit none
        
        !! input
        real(8), dimension(ndofel,1)            :: rhs
        real(8), dimension(ndofel,ndofel)       :: amatrx
        real(8), dimension(nsvars)              :: svars
        real(8), dimension(nprops)              :: props
        real(8), dimension(mcrd,nnodes)         :: coords
        real(8), dimension(ndofel)              :: u
        real(8), dimension(2)                   :: time
        real(8), dimension(2,npredf,nnodes)     :: predef
        integer, dimension(njprops)             :: jprops
        real(8) :: dtime, pnewdt
        integer :: ndofel, nsvars, nprops, mcrd, nnodes, npredf, njprops
        integer :: jtype, kstep, kinc, jelem
        integer :: nthreads, ithread
        
        !! local
        type(tri3_type) :: e
        real(8), dimension(STRS_2D,STRS_2D) :: cmat
        
        !! begin
        call e%initElement(jelem,jtype,props,nprops,jprops,njprops,svars,nsvars,u,ndofel)
        call e%computeCurrentCoords(coords)
        call e%integrateElement(jelem, svars, nsvars)
        call e%assembleElement(amatrx, rhs, ndofel)
        
    end subroutine tri3
    
    
    ! ##############################################################################################
    
    
    subroutine initElement(self,jelem,jtype,props,nprops,jprops,njprops,svars,nsvars,u,ndofel)
        !! Initialize element properties and settings.
        implicit none
        
        !! input
        class(tri3_type), intent(inout)         :: self
        integer, intent(in)                     :: jelem
        integer, intent(in)                     :: jtype
        real(8), dimension(nprops), intent(in)  :: props
        integer, dimension(njprops), intent(in) :: jprops
        real(8), dimension(nsvars), intent(in)  :: svars
        real(8), dimension(ndofel), intent(in)  :: u
        integer, intent(in)                     :: nprops
        integer, intent(in)                     :: njprops
        integer, intent(in)                     :: nsvars
        integer, intent(in)                     :: ndofel
        
        !! begin
        if (jtype == 21) then
            self%elemtype = jtype
            self%pstype   = PSTYPE_STRN
        elseif (jtype == 25) then
            self%elemtype = jtype
            self%pstype   = PSTYPE_STRS
        else
            print "(1x,a,i0)", "Error, unrecognized user element type: ", self%elemtype
            print "(1x,a,i0)", "jelem: ", jelem
            call exit()
        endif
        
        if (ndofel /= NDOF_TRI3) then
            print "(1x,a,i0)", "Error, user element is incompatible with NDOF: ", ndofel
            print "(1x,a,i0)", "jelem: ", jelem
            call exit()
        endif
        
        if ( any([1,2] == jprops(1)) ) then
            !! 1:Isotropic | 2:Lamina
            self%matbehav = jprops(1)
        else
            print "(1x,a,i0)", "Error, unexpected material behavior: ", jprops(1)
            print "(1x,a,i0)", "jelem: ", jelem
            call exit()
        endif
        
        self%props(1:4) = props(1:4)    !! (E1, E2, nu12, G12)
        self%density    = props(5)
        self%thick      = props(6)
        self%nodes      = [1:self%nnodes]
        self%u          = u(1:ndofel)
        
        self%gpwgt      = [one]
        self%gpcrds     = reshape([one/three, one/three], shape=[DIMS_2D, NINTEGRPTS_TRI3])
        
    end subroutine initElement
    
    
    ! ##############################################################################################
    
    
    subroutine computeCurrentCoords(self, coords)
        !! Compute the current nodal coordinates.
        implicit none
        
        !! input
        class(tri3_type), intent(inout)     :: self
        real(8), dimension(:,:), intent(in) :: coords
        
        !! local
        integer :: inode, icrd
        
        !! begin
        do inode = 1,self%nnodes                                                                    !! for each node,
            do icrd = 1,DIMS_2D                                                                     !! for each spatial coordinate direction,
                self%curcrds(icrd,inode) = coords(icrd,inode) + self%u(icrd + DIMS_2D*(inode-1))    !! compute the current coordinates & assign to the corresponding node on the subelement
            enddo                                                                                   !! example: inode=2; icrd=2; self%curcrds(2,2) = coords(2,2) + self%u(4) = YCURRENT @ node 2
        enddo
        
    end subroutine computeCurrentCoords
    
    
    ! ##############################################################################################
    
    
    subroutine computeCmat(self, jelem, cmat, svars, nsvars)
        !! Build constitutive/elasticity matrix in local coordinates. If isotropic or aligned with 
        !! global csys, rotation not needed.
        
        implicit none
        
        !! input
        class(tri3_type), intent(inout)        :: self
        integer, intent(in)                    :: jelem
        real(8), dimension(:,:), intent(out)   :: cmat
        real(8), dimension(nsvars), intent(in) :: svars
        integer, intent(in)                    :: nsvars
        
        !! local
        real(8) :: e1, e2, nu12, g12
        real(8) :: nu21
        
        !! begin
        cmat = zero
        
        if (self%matbehav == 1) then
            !! ISOTROPIC
            e1   = self%props(1)
            nu12 = self%props(3)
            
            if (self%pstype == 1) then
                !! Plane strain
                cmat(1,1) = (one-nu12)
                cmat(1,2) = nu12
                cmat(2,1) = cmat(1,2)
                cmat(2,2) = cmat(1,1)
                cmat(3,3) = (one-two*nu12) / two
                cmat = (e1 / ((one+nu12)*(one-two*nu12))) * cmat
                
            elseif (self%pstype == 2) then
                !! Plane stress
                cmat(1,1) = e1 / (one-nu12*nu12)
                cmat(1,2) = e1*nu12 / (one-nu12*nu12)
                cmat(2,1) = cmat(1,2)
                cmat(2,2) = cmat(1,1)
                cmat(3,3) = e1*(one-nu12) / (two*(one-nu12*nu12))
                
            else
                print "(1x,a,i0)", "Error, unsupported plane stress/strain type: ", self%pstype
                print "(1x,a,i0)", "jelem: ", jelem
                call exit()
            endif
            
        elseif (self%matbehav == 2) then
            !! LAMINA
            e1   = self%props(1)
            e2   = self%props(2)
            nu12 = self%props(3)
            g12  = self%props(4)
            nu21 = nu12 * e2/e1
            
            if (self%pstype == 1) then
                !! Plane strain
                print *, "Error, computeCmat not implemented for jprops(1)=2 (matbehav=LAMINA) in plane strain."
                print "(1x,a,i0)", "jelem: ", jelem
                call exit()
            
            elseif (self%pstype == 2) then
                !! Plane stress
                cmat(1,1) = e1/(one-nu12*nu21)
                cmat(1,2) = (e1*nu21)/(one-nu12*nu21)
                cmat(2,1) = cmat(1,2)
                cmat(2,2) = e2/(one-nu12*nu21)
                cmat(3,3) = g12
                
            else
                print "(1x,a,i0)", "Error, unsupported plane stress/strain type: ", self%pstype
                print "(1x,a,i0)", "jelem: ", jelem
                call exit()
            endif
            
        elseif (self%matbehav == 3) then
            !! ENGINEERING_CONSTANTS
            print *, "Error, computeCmat not implemented for jprops(1)=3 (matbehav=ENGINEERING_CONSTANTS)."
            print "(1x,a,i0)", "jelem: ", jelem
            call exit()
            
        else
            print "(1x,a,i0)", "Error, unsupported material behavior: ", self%matbehav
            print "(1x,a,i0)", "jelem: ", jelem
            call exit()
        endif
        
    end subroutine computeCmat
    
    
    ! ##############################################################################################
    
    
    subroutine integrateElement(self, jelem, svars, nsvars)
        !! Compute element strains, stresses, RHS vector, and stiffness matrix.
        
        implicit none
        
        !! input
        class(tri3_type)                         :: self
        integer, intent(in)                      :: jelem
        real(8), dimension(nsvars), intent(in)   :: svars
        integer, intent(in)                      :: nsvars
        
        !! local
        integer                                  :: igp, j, k                                       !! counters (specializations for TRI3)
        real(8)                                  :: gpwgt, y23, y31, y12, x32, x13, x21, x31, y21   !! needed by TRI3 only
        real(8)                                  :: detjac, halfdetjac                              !! "1/2" factor for triangular elems.
        
        real(8), dimension(STRS_2D)              :: strain
        real(8), dimension(STRS_2D, STRS_2D)     :: cmat
        real(8), dimension(STRS_2D, self%ndof)   :: bmat
        real(8), dimension(self%ndof, STRS_2D)   :: bt
        real(8), dimension(self%ndof, STRS_2D)   :: btc
        real(8), dimension(self%ndof, self%ndof) :: btcb
        real(8), dimension(self%ndof)            :: bts
        
        
        !! begin 
        self%stiff  = zero
        self%rhs    = zero
        self%gpstrs = zero
        
        !! Build the stress/strain constitutive matrix (local coords).
        call self%computeCmat(jelem, cmat, svars, nsvars)
        
        !! Compute strain-displacement matrix, strain, stress, and contrib to stiffness matrix and rhs,
        !! at each GP (tri3 has only 1).
        gploop: do igp = 1,self%ngp
            !! Reset values for each GP.
            bmat   = zero       !! strain-disp, b = invjac * dfncs
            bt     = zero       !! transpose(b)
            btc    = zero       !! transpose(b) * cmat
            btcb   = zero       !! transpose(b) * cmat * b
            bts    = zero       !! transpose(b) * stress vector 
            detjac = zero
            
            gpwgt = self%gpwgt(igp)
            
            !! TRI3 allows simplified calcs using only curcrds(icrd,inode).
            !! Notation example: y23 = y2 - y3
            y23 = self%curcrds(2,2) - self%curcrds(2,3)
            y31 = self%curcrds(2,3) - self%curcrds(2,1)
            y12 = self%curcrds(2,1) - self%curcrds(2,2)
            y21 = self%curcrds(2,2) - self%curcrds(2,1)
            
            x32 = self%curcrds(1,3) - self%curcrds(1,2)
            x13 = self%curcrds(1,1) - self%curcrds(1,3)
            x21 = self%curcrds(1,2) - self%curcrds(1,1)
            x31 = self%curcrds(1,3) - self%curcrds(1,1)
            
            detjac = x21*y31 - x31*y21
            
            !! Compute the strain-displacement matrix.
            bmat(1,1) = y23/detjac
            bmat(3,1) = x32/detjac
            
            bmat(2,2) = bmat(3,1) ! x32/detjac
            bmat(3,2) = bmat(1,1) ! y23/detjac
            
            bmat(1,3) = y31/detjac
            bmat(3,3) = x13/detjac
            
            bmat(2,4) = bmat(3,3) ! x13/detjac
            bmat(3,4) = bmat(1,3) ! y31/detjac
            
            bmat(1,5) = y12/detjac
            bmat(3,5) = x21/detjac
            
            bmat(2,6) = bmat(3,5) ! x21/detjac
            bmat(3,6) = bmat(1,5) ! y12/detjac
            
            bt = transpose(bmat)
            btc = matmul(bt, cmat)
            btcb = matmul(btc,bmat)
            
            !! Add gp contribution to stiffness.
            halfdetjac = (one/two) * detjac
            self%stiff = self%stiff + (btcb * halfdetjac * gpwgt * self%thick)
            
            do j = 1,STRS_2D
                strain(j) = zero
                do k = 1,self%ndof
                    strain(j) = strain(j) + (bmat(j,k) * self%u(k))
                enddo
            enddo
            
            do j = 1,STRS_2D
                self%gpstrs(j,igp) = zero
                do k = 1,STRS_2D
                    self%gpstrs(j,igp) = self%gpstrs(j,igp) + (cmat(j,k) * strain(k))
                enddo
            enddo
            
            do j = 1,self%ndof
                do k = 1,STRS_2D
                    bts(j) = bts(j) + (bmat(k,j) * self%gpstrs(k,igp))
                enddo
            enddo
            
            self%rhs = self%rhs + (bts * halfdetjac * gpwgt * self%thick)
            
        enddo gploop
        
    end subroutine integrateElement
    
    
    ! ##############################################################################################
    
    
    subroutine assembleElement(self, amatrx, rrhs, ndofel)
        !! Assemble to a) tangent stiffness matrix (amatrx), and b) residual force vector (rhs).
        
        implicit none
        
        !! input
        class(tri3_type), intent(in)                :: self
        real(8), dimension(:,:), intent(inout)      :: amatrx
        real(8), dimension(ndofel,1), intent(inout) :: rrhs
        integer, intent(in)                         :: ndofel
        
        !! local
        integer                                     :: inode, icrd, n, i, j
        integer                                     :: GlobalI, GlobalJ
        integer, dimension(self%ndof)               :: GlobalDof
        
        !! begin
        do inode = 1, self%nnodes
            n = self%nodes(inode)
            do icrd = 1, DIMS_2D
                GlobalDof(icrd + DIMS_2D*(inode-1)) = icrd + DIMS_2D*(n-1)
            end do
        end do
        
        do i = 1, self%ndof
            GlobalI = GlobalDof(i)
            rrhs(GlobalI,1) = rrhs(GlobalI,1) - self%rhs(i)
            do j = 1, self%ndof
                GlobalJ = GlobalDof(j)
                amatrx(GlobalI,GlobalJ) = amatrx(GlobalI,GlobalJ) + self%stiff(i,j)
            end do
        end do
        
    end subroutine assembleElement
    
    
    ! ##############################################################################################
    
    
end module utri3
