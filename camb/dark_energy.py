from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np, \
    AllocatableArrayDouble, f_pointer
from ctypes import c_int, c_double, byref, POINTER, c_bool

@fortran_class
class DarkEnergyModel(F2003Class):
    
    _fortran_class_module_ = 'DarkEnergyInterface'
    _fortran_class_name_ = 'TDarkEnergyModel'
    
    _fields_ = [
        ("__is_cosmological_constant", c_bool),
        ("__num_perturb_equations", c_int),
        ("w", c_double, "w(0)"),
        ("wa", c_double, "-dw/da(0)"),
        ("w0", c_double),
        ("w1", c_double),
        ("w2", c_double),
        ("w3", c_double),
        ("abound1", c_double),
        ("abound2", c_double),
        ("abound3", c_double),
        ("amid", c_double),
        ("alpha", c_double*5),
        ("sim", c_int),
        ("c_Gamma_ppf", c_double, "-dw/da(0)"),
        ("__no_perturbations", c_bool, "turn off perturbations (unphysical, so hidden in Python)")
    ]
                
    def validate_params(self):
        if self.wa + self.w > 0:
            raise CAMBError('dark energy model has w + wa > 0, giving w>0 at high redshift')
        return True

    def set_params(self, w=-1.0, wa=0, w0=0, w1=0, w2=0, w3=0, sim=1, alpha = [.0, 0.0, 0.0, 0.0, 0.0]):
        self.w = w
        self.wa = wa
        self.w0 = w0
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
        self.abound1 = abound1
        self.abound2 = abound2
        self.abound3 = abound3
        self.amid = amid
        self.alpha = alpha
        self.sim = sim
        #self.validate_params()

# short names for models that support w/wa
F2003Class._class_names.update({'ppf': DarkEnergyModel})
