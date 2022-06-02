module DarkEnergyInterface
    use precision
    use interpolation
    use classes
    implicit none

private

    type, extends(TCambComponent) :: TDarkEnergyModel
        logical :: is_cosmological_constant = .true.
        integer :: num_perturb_equations = 0
        real(dl) :: w_lam = -1_dl !p/rho for the dark energy (an effective value, used e.g. for halofit)
        real(dl) :: wa = 0._dl !may not be used, just for compatibility with e.g. halofit
        real(dl) :: cs2_lam = 1_dl !rest-frame sound speed, though may not be used
        logical :: no_perturbations = .false. !Don't change this, no perturbations is unphysical
    contains
        procedure :: Init
        procedure :: ReadParams 
        procedure :: w_de
        procedure :: grho_de
        procedure :: grhot_de
        procedure :: Effective_w_wa  ! Used as approximate values for non-linear corrections
        procedure :: PrintFeedback
        procedure :: PerturbationInitial

        procedure :: PerturbedStressEnergy !Get density perturbation and heat flux for sources
        procedure :: diff_rhopi_Add_Term
        procedure :: PerturbationEvolve
        
        procedure, nopass :: SelfPointer
        procedure, nopass :: PythonClass
        procedure, private :: setcgammappf
    end type TDarkEnergyModel

    public TDarkEnergyModel
    
contains

    subroutine Init(this, State)
        
        use classes
        class(TDarkEnergyModel), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State

        this%is_cosmological_constant = .false.
        this%num_perturb_equations = 1

    end subroutine Init

    subroutine ReadParams(this, Ini)
        
        use IniObjects
        use FileUtils
        class(TDarkEnergyModel) :: this
        class(TIniFile), intent(in) :: Ini

        this%w_lam = Ini%Read_Double('w', -1.d0)
        this%wa = Ini%Read_Double('wa', 0.d0)
        
        if (this%w_lam + this%wa > 0) then
            error stop 'w + wa > 0, giving w>0 at high redshift'
        endif

    end subroutine ReadParams

    function w_de(this, a)

        class(TDarkEnergyModel) :: this
        real(dl) :: w_de, al
        real(dl), intent(IN) :: a

        w_de = this%w_lam+ this%wa*(1._dl-a)

    end function w_de  

    function grho_de(this, a) result(res) ! relative density (8 pi G a^4 rho_de / grhov)
        
        class(TDarkEnergyModel) :: this
        real(dl), intent(IN) :: a
        real(dl) :: res

        res = a ** (1._dl - 3. * this%w_lam - 3. * this%wa)
        if (this%wa /= 0) then
            res = res*exp(-3. * this%wa * (1._dl - a))
        endif
        
    end function grho_de

    function grhot_de(this, a) result(res) ! Get grhov_t = 8*pi*rho_de*a**2
        
        class(TDarkEnergyModel), intent(inout) :: this
        real(dl), intent(IN) :: a
        real(dl) :: res

        if (a > 1e-10) then ! Ensure a valid result
            res = this%grho_de(a) / (a * a)
        else
            res = 0._dl
        end if

    end function grhot_de

    subroutine Effective_w_wa(this, w, wa)
        
        class(TDarkEnergyModel), intent(inout) :: this
        real(dl), intent(out) :: w, wa
        w = this%w_lam
        wa = this%wa
    
    end subroutine Effective_w_wa

    subroutine PrintFeedback(this, FeedbackLevel)
        
        class(TDarkEnergyModel) :: this
        integer, intent(in) :: FeedbackLevel

        if (FeedbackLevel >0) then
            write(*,'("(w0, wa) = (", f8.5,", ", f8.5, ")")') this%w_lam, this%wa
        endif 

    end subroutine PrintFeedback

    subroutine PerturbationInitial(this, y, a, tau, k) !Get intinitial values for perturbations at a (or tau)
        
        class(TDarkEnergyModel), intent(in) :: this
        real(dl), intent(out) :: y(:)
        real(dl), intent(in) :: a, tau, k
        y = 0 !For standard adiabatic perturbations can usually just set to zero to good accuracy
    
    end subroutine PerturbationInitial

    subroutine PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe=0
    dgqe=0

    end subroutine PerturbedStressEnergy

    function diff_rhopi_Add_Term(this, dgrhoe, dgqe,grho, gpres, w, grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, grhok, w, adotoa, &
        k, grhov_t, z, k2, yprime(:), y(:), Kf1
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot

    ! Ensure, that the result is set, when the function is not implemented by
    ! subclasses
    ppiedot = 0._dl

    end function diff_rhopi_Add_Term

    subroutine PerturbationEvolve(this, ayprime, w, w_ix, a, adotoa, k, z, y)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a,adotoa, k, z, y(:), w
    integer, intent(in) :: w_ix
    end subroutine PerturbationEvolve

    function TDarkEnergyModel_PythonClass()
        
        character(LEN=:), allocatable :: TDarkEnergyModel_PythonClass
        TDarkEnergyModel_PythonClass = 'DarkEnergyModel'
    
    end function TDarkEnergyPPF_PythonClass

    subroutine TDarkEnergyPPF_SelfPointer(cptr,P)
    
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TDarkEnergyModel), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    
    end subroutine TDarkEnergyPPF_SelfPointer

   subroutine setcgammappf(this)
    
        class(TDarkEnergyModel) :: this
        this%c_Gamma_ppf = 0.4_dl 
    
    end subroutine setcgammappf

end module DarkEnergyInterface
