from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np, \
    AllocatableArrayDouble, f_pointer
from ctypes import c_int, c_double, byref, POINTER, c_bool

@fortran_class
# class DarkEnergyModel(F2003Class):
    
    _fortran_class_module_ = 'DarkEnergyInterface'
    _fortran_class_name_ = 'TDarkEnergyModel'
    
    _fields_ = [
        ("__is_cosmological_constant", c_bool),
        ("__num_perturb_equations", c_int),
        ("w", c_double, "w(0)"),
        ("wa", c_double, "-dw/da(0)"),
        ("c_Gamma_ppf", c_double, "-dw/da(0)"),
        ("__no_perturbations", c_bool, "turn off perturbations (unphysical, so hidden in Python)")
    ]
                
    def validate_params(self):
        if self.wa + self.w > 0:
            raise CAMBError('dark energy model has w + wa > 0, giving w>0 at high redshift')
        return True

    def set_params(self, w=-1.0, wa=0):
        self.w = w
        self.wa = wa
        self.validate_params()

# short names for models that support w/wa
F2003Class._class_names.update({'ppf': DarkEnergyModel})
