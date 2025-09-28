import numpy as np
from typing import Protocol
from chemicals.elements import molecular_weight as MW_data_import 
from chemicals.elements import nested_formula_parser as MW_data_parser

class PropertyPackageInterface(Protocol): 
    
    """Common interface for equation-of-state backends."""
    name: str
    components: tuple[str, ...]       
    
    def get_compressibility_factor(self, 
              temperature_K: float, 
              pressure_Pa: float, 
              molar_composition: np.ndarray) -> float: ...
    
    def get_density_SI(self, 
              temperature_K: float, 
              pressure_Pa: float, 
              molar_composition: np.ndarray) -> float: ...


class IdealGasBackend:

    """ Ideal gas property package backend. """

    name = 'IdealGas'
    
    def __init__(self, 
                 components: tuple[str, ...]):
        self.components = components

    def get_compressibility_factor(self, 
              temperature_K: float, 
              pressure_Pa: float, 
              molar_composition: np.ndarray,) -> float:
        return 1.0
    
    def get_density_SI(self,
                       temperature_K: float,
                       pressure_Pa: float,
                       molar_composition: np.ndarray) -> float:
        
        molecular_weight_vec_g_mol = np.array([MW_data_import(MW_data_parser(f)) for f in self.components], dtype=float) 
        
        R = 8.314462618 # J/(mol·K)
        P = pressure_Pa
        T = temperature_K
        MW = np.dot(molecular_weight_vec_g_mol, molar_composition) / 1000 # kg/mol
        rho = P * MW / (R * T) 

        density_SI = rho 

        return density_SI
    

class SRK:
    name = 'SRK'
    
    def __init__(self, components: tuple[str, ...]):
        self.components = components

    def get_compressibility_factor(self, 
              temperature_K: float, 
              pressure_Pa: float, 
              molar_composition: np.ndarray) -> float:
        # Placeholder implementation
        return 0.9
