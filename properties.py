import numpy as np
from typing import Protocol
from chemicals.elements import molecular_weight as MW_data_import 
from chemicals.elements import nested_formula_parser as MW_data_parser

import thermo

class EquationOfStateInterface(Protocol): 
    
    """ Common interface for equation-of-state backends. """
    name: str
    components: tuple[str, ...]       
    
    def get_fugacity_coefs(self,
                           temperature_K: float, 
                           pressure_Pa: float, 
                           molar_composition: np.ndarray) -> np.ndarray: ...

    def get_compressibility_factor(self, 
                                   temperature_K: float,
                                   pressure_Pa: float,
                                   molar_composition: np.ndarray) -> float: ...
    
    def get_density_SI(self,
                       temperature_K: float,
                       pressure_Pa: float,
                       molar_composition: np.ndarray) -> float: ...
    

class ActivityModelInterface(Protocol):

    """ Common interface for activity model backends. """
    name: str
    components: tuple[str, ...]       
    
    def get_activity(self, 
                     temperature_K: float,
                     pressure_Pa: float,
                     molar_composition: np.ndarray) -> np.ndarray: ...





class GammaPhiPackage():

    """ Property package to calculate VLE using activity and EOS models. """
    def __init__(self, 
                 eos_backend: EquationOfStateInterface,
                 activity_model_backend: ActivityModelInterface):
        self.eos_backend = eos_backend
        self.activity_model_backend = activity_model_backend
        self.components = eos_backend.components

    def TP_flash(self): 

        " Isothermal-isobaric flash calculation. "




        return None



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
        
        molecular_weight_vec_g_mol = np.array([MW_data_import(MW_data_parser(f)) for f in self.components]) 
        
        R: float = 8.314462618 # J/(mol·K)
        P: float = pressure_Pa
        T: float = temperature_K
        MW: float = np.dot(molecular_weight_vec_g_mol, molar_composition) / 1000 # kg/mol
        rho: float = P * MW / (R * T) 

        density_SI: float = rho 

        return density_SI
    

class SoaveRedlichKwongBackend:
    name = 'SRK'
    
    def __init__(self, components: tuple[str, ...]):
        self.components = components

    def get_compressibility_factor(self, 
              temperature_K: float, 
              pressure_Pa: float, 
              molar_composition: np.ndarray) -> float:
        # Placeholder implementation
        return 0.9
