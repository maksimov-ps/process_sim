import numpy as np
from typing import Protocol

from chemicals import CAS_from_any
from chemicals import MW as MW_data_import
from chemicals.acentric import omega as acentric_factor_data_import
from chemicals.critical import Tc as critical_temperature_data_import, Pc as critical_pressure_data_import



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
    



class PropertyPackageInterface(Protocol):

    """ Common interface for property packages. """
    eos_backend: EquationOfStateInterface
    activity_model_backend: ActivityModelInterface

    def _TPD_stability_test(self): ...







class GammaPhiPackage():

    """ Property package to calculate VLE using activity and EOS models. """
    def __init__(self, 
                 eos_backend: EquationOfStateInterface,
                 activity_model_backend: ActivityModelInterface):
        self.eos_backend = eos_backend
        self.activity_model_backend = activity_model_backend


    def _TPD_stability_test(self, 
                            temperature_K: float,
                            pressure_Pa: float,
                            components: tuple[str, ...],
                            molar_composition: np.ndarray) -> bool:
        " Performs a Tangent Plane Stability test for a given conditions (T, P, z). "
        " The method is based on Michelsen (1982)"

        fugacity_coefs: np.ndarray = self.eos_backend.get_fugacity_coefs(temperature_K = temperature_K,
                                                                         pressure_Pa = pressure_Pa,
                                                                         components = components,
                                                                         molar_composition = molar_composition)




        return True


    def TP_flash(self, 
                 temperature_K: float,
                 pressure_Pa: float,
                 components: tuple[str, ...],
                 molar_composition: np.ndarray ) -> None: 

        " Isothermal-isobaric flash calculation. "
        TPD_stability: bool = self._TPD_stability_test(temperature_K = temperature_K,
                                                       pressure_Pa = pressure_Pa,
                                                       components = components,
                                                       molar_composition = molar_composition)



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
        
        # placeholder implementation
        return 0.9
    

class SoaveRedlichKwongBackend():
    name = 'SRK'
    

    def get_compressibility_factor(self, 
              temperature_K: float, 
              pressure_Pa: float, 
              molar_composition: np.ndarray) -> float:
        # Placeholder implementation
        return 0.9
    

    def get_fugacity_coefs(temperature_K: float, 
                           pressure_Pa: float, 
                           components: tuple[str, ...],
                           molar_composition: np.ndarray) -> np.ndarray:
        
        
        acentric_factor_data = np.zeros_like(molar_composition)
        critical_temperature_K_data = np.zeros_like(molar_composition)
        critical_pressure_Pa_data = np.zeros_like(molar_composition)

        for k in range(len(components)):
            compoent_CASRN = CAS_from_any(components[k])
            acentric_factor_data[k] = acentric_factor_data_import(CASRN = compoent_CASRN)  # Placeholder value
            critical_temperature_K_data[k] = critical_temperature_data_import(CASRN = compoent_CASRN)  # Placeholder value
            critical_pressure_Pa_data[k] = critical_pressure_data_import(CASRN = compoent_CASRN)  # Placeholder value

        
        
        
        # Placeholder implementation
        return np.ones_like(molar_composition)


class ActivityModelBackend():
    name = 'ActivityModel'
    
    def __init__(self):
        pass

    def get_activity(self, 
                     temperature_K: float,
                     pressure_Pa: float,
                     molar_composition: np.ndarray) -> np.ndarray:
        # Placeholder implementation
        return np.ones_like(molar_composition)