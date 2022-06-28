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
        real(dl) :: wa = 0._dl 
        real(dl) :: w0 = 0._dl
        real(dl) :: w1 = 0._dl
        real(dl) :: w2 = 0._dl
        real(dl) :: w3 = 0._dl
        real(dl) :: abound1 = 0._dl
        real(dl) :: abound2 = 0._dl
        real(dl) :: abound3 = 0._dl
        real(dl) :: amid = 0._dl
        real(dl), dimension (5) :: alpha
        integer :: sim = 1
        real(dl) :: c_Gamma_ppf = 0.4_dl
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
        procedure, nopass :: SelfPointer => TDarkEnergyModel_SelfPointer
        procedure, nopass :: PythonClass => TDarkEnergyModel_PythonClass
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
        this%w0 = Ini%Read_Double('w0', 0.d0)
        this%w1 = Ini%Read_Double('w1', 0.d0)
        this%w2 = Ini%Read_Double('w2', 0.d0)
        this%w3 = Ini%Read_Double('w3', 0.d0)
        this%abound1 = Ini%Read_Double('abound1', 0.d0)
        this%abound2 = Ini%Read_Double('abound2', 0.d0)
        this%abound3 = Ini%Read_Double('abound3', 0.d0)
        this%amid = Ini%Read_Double('amid', 0.d0)
        this%alpha = Ini%Read_Double_Array('alpha', 1)
        this%sim = Ini%Read_Int('sim', 1)
        
       ! if (this%w_lam + this%wa > 0) then
        !    error stop 'w + wa > 0, giving w>0 at high redshift'
        !endif

    end subroutine ReadParams


    function w_de(this, a)

        
        class(TDarkEnergyModel) :: this
        real(dl) :: w_de, al, i, j, temp, a_min
        real(dl), intent(IN) :: a
        
        a_min=0
        
        if (this%sim==0) then !wa, w0 model (original)
            if (a<0.5) then
            w_de = this%w_lam
            else 
            w_de = this%wa
            end if
            
            if (a<a_min) then
                w_de=-1
            end if
            
        end if
        
        if (this%sim==1)then  !two-step model 
            
            
            if (a>this%amid) then
                w_de = this%w_lam
            else 
                w_de = this%wa
            end if
            
            if (a<a_min) then
                w_de=-1
            end if
            
        end if
        
        if (this%sim==2)then   !four-step model 
                     
            
            if (this%abound3<a) then
                w_de = this%w3
            end if
            
            if (a<this%abound3) then
                w_de = this%w2
            end if
            
            if (a<this%abound2) then
                w_de = this%w1
             
            end if
         
            if (a<this%abound1) then
                w_de = this%w0
                
            end if
            
            if (a<a_min) then
                w_de=-1
            end if
            
         end if
         
         if (this%sim==3)then  !natural log model 
         
             if (a<a_min) then
                  w_de=-1
             end if
             
              i=1
              j=5
              temp=0
              do while (i<=j)
                  temp=temp + (this%alpha(i)*(log(a**-1))**(i)) 
                  i=i+1
              end do  
               
              w_de = -1 + temp
              
              
            
          end if

    end function w_de  
    
           

    function grho_de(this, a) result(res) ! relative density (8 pi G a^4 rho_de / grhov)
        
        class(TDarkEnergyModel) :: this
        real(dl), intent(IN) :: a
        real(dl) :: res, i, j, temp, a_min
        
        a_min=0
        
        if (this%sim==0)then !wa, w0 model (original)
            if (this%wa /= 0) then
                res = res*exp(-3. * this%wa * (1._dl - a))
            
            endif
            
        end if
        
        if (this%sim==1)then !two-step model
           
          
        
            if (a > this%amid) then
                res = a**(-3*this%w_lam + 1)
            else
                res = a**(-3*this%wa + 1) * this%amid**(-3*(this%w_lam - this%wa))
            end if
            
            if (a<a_min) then
                res=a**4
            end if
        
        end if
        
        if (this%sim==2)then !four-step model
            
                      
            
            
            if (this%abound3<a) then
                 res = a**(-3*this%w3 + 1)
                 
            end if
            
           
            if (a<this%abound3) then
                res = a**(-3*this%w2 + 1) * this%abound3**(-3*(this%w3 - this%w2))
            end if
            
            if (a<this%abound2) then
                res = a**(-3*this%w1 + 1) * this%abound3**(-3*(this%w3 - this%w2))* this%abound2**(-3*(this%w2-this%w1))
            end if
         
            if (a<this%abound1) then
                res = a**(-3*this%w0 + 1) *this%abound3**(-3*(this%w3 - this%w2)) * &
                     this%abound2**(-3*(this%w2 - this%w1)) * this%abound1**(-3*(this%w1 - this%w0))
            end if
          
            
            if (a<a_min) then
                res=a**4
            end if
        
        end if
        
        if (this%sim==3)then !natural log model
            
            if (a<a_min) then
                  res = a**(4)
            end if
              
              i=1
              j=5
              temp=0
              do while (i<=j)
                  temp=temp + (this%alpha(i)*(log(a**-1))**(i))/(i+1) 
                  i=i+1
              end do  
               
              res = a**(-3*temp+4)  
         end if 
      
         
     
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

    subroutine PerturbedStressEnergy(this, dgrhoe, dgqe, a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        
        class(TDarkEnergyModel), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix
        real(dl) :: Gamma, S_Gamma, ckH, Gammadot, Fa, sigma
        real(dl) :: vT, grhoT, k2

        k2=k**2
        !ppf
        grhoT = grho - grhov_t
        vT = dgq / (grhoT + gpres_noDE)
        Gamma = ay(w_ix)

        !sigma for ppf
        sigma = (etak + (dgrho + 3 * adotoa / k * dgq) / 2._dl / k) / kf1 - &
            k * Gamma
        sigma = sigma / adotoa

        S_Gamma = grhov_t * (1 + w) * (vT + sigma) * k / adotoa / 2._dl / k2
        ckH = this%c_Gamma_ppf * k / adotoa

        if (ckH * ckH .gt. 3.d1) then ! ckH^2 > 30 
            Gamma = 0
            Gammadot = 0.d0
        else
            Gammadot = S_Gamma / (1 + ckH * ckH) - Gamma - ckH * ckH * Gamma
            Gammadot = Gammadot * adotoa
        endif
        ayprime(w_ix) = Gammadot !Set this here, and don't use PerturbationEvolve

        Fa = 1 + 3 * (grhoT + gpres_noDE) / 2._dl / k2 / kf1
        dgqe = S_Gamma - Gammadot / adotoa - Gamma
        dgqe = -dgqe / Fa * 2._dl * k * adotoa + vT * grhov_t * (1 + w)
        dgrhoe = -2 * k2 * kf1 * Gamma - 3 / k * adotoa * dgqe

    end subroutine PerturbedStressEnergy

    function diff_rhopi_Add_Term(this, dgrhoe, dgqe,grho, gpres, w, grhok, adotoa, Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
        
        class(TDarkEnergyModel), intent(in) :: this
        real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, w, grhok, adotoa 
        real(dl), intent(in) :: k, grhov_t, z, k2, yprime(:), y(:), Kf1
        integer, intent(in) :: w_ix
        real(dl) :: ppiedot, hdotoh

        hdotoh = (-3._dl * grho - 3._dl * gpres - 2._dl * grhok) / 6._dl / adotoa
        
        ppiedot = 3._dl * dgrhoe + dgqe * &
            (12._dl / k * adotoa + k / adotoa - 3._dl / k * (adotoa + hdotoh)) + &
            grhov_t * (1 + w) * k * z / adotoa - 2._dl * k2 * Kf1 * &
            (yprime(w_ix) / adotoa - 2._dl * y(w_ix))
        
        ppiedot = ppiedot * adotoa / Kf1

    end function diff_rhopi_Add_Term

    function TDarkEnergyModel_PythonClass()
        
        character(LEN=:), allocatable :: TDarkEnergyModel_PythonClass
        TDarkEnergyModel_PythonClass = 'DarkEnergyModel'
    
    end function TDarkEnergyModel_PythonClass

    subroutine TDarkEnergyModel_SelfPointer(cptr, P)
    
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TDarkEnergyModel), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    
    end subroutine TDarkEnergyModel_SelfPointer

    subroutine setcgammappf(this)
    
        class(TDarkEnergyModel) :: this
        this%c_Gamma_ppf = 0.4_dl 
    
    end subroutine setcgammappf

end module DarkEnergyInterface
