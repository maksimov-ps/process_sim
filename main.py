import numpy as np 
from core import Stream
from properties import GammaPhiPackage
from properties import SoaveRedlichKwongBackend, ActivityModelBackend


if __name__ == "__main__":
    
    feed = Stream(pressure_Pa = 1e5,
                  temperature_K = 300.0,
                  components = ('H2', 'N2', 'NH3', 'CH3OH','C2H5OH'),
                  molar_composition = np.array([0.7, 0.3, 0.0, 0.0, 0.0]),
                  molar_flow_mol_s = 1.0,
                  property_package_backend = GammaPhiPackage(eos_backend = SoaveRedlichKwongBackend, 
                                                             activity_model_backend = ActivityModelBackend))
    
    print(feed.vapour_fraction())

    



