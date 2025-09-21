import numpy as np
from properties import PropertyPackageInterface
from dataclasses import dataclass, field    

@dataclass
class Stream: 

    """ State container to store properties of streams in between unit ops. """

    pressure_Pa: float                  
    temperature_K: float                     
    molar_composition: np.ndarray       # order according to the backend property package
    molar_flow_mol_s: float             
    backend: PropertyPackageInterface = field(default=None)        

    def Z(self) -> float:
        return self.backend.get_Z(self.temperature_K, self.pressure_Pa, self.molar_composition)
    
    def density_SI(self) -> float:
        return self.backend.get_density_SI(self.temperature_K, self.pressure_Pa)

    