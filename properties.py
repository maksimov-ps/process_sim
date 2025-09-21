import numpy as np
from typing import Protocol


class PropertyPackageInterface(Protocol): 
    
    """Common interface for equation-of-state backends."""
    name: str
    components: tuple[str, ...]       
    
    def get_Z(self, 
              temperature_K: float, 
              pressure_Pa: float, 
              molar_composition: np.ndarray) -> float: ...


class IdealGasBackend:
    name = 'IdealGas'
    
    def __init__(self, 
                 components: tuple[str, ...]):
        self.components = components

    def get_Z(self, 
              temperature_K: float, 
              pressure_Pa: float, 
              molar_composition: np.ndarray) -> float:
        return 1.0
    
    def get_density_SI(self,
                       temperature_K: float,
                       pressure_Pa: float) -> float:
        R = 8.314462618 # J/(mol·K)
        return pressure_Pa / (R * temperature_K)
    

class SRK:
    name = 'SRK'
    
    def __init__(self, components: tuple[str, ...]):
        self.components = components

    def get_Z(self, 
              temperature_K: float, 
              pressure_Pa: float, 
              molar_composition: np.ndarray) -> float:
        # Placeholder implementation
        return 0.9
