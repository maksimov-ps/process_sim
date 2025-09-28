import numpy as np 
from core import Stream
from properties import IdealGasBackend, SRK

if __name__ == "__main__":

    feed = Stream(backend = IdealGasBackend(components=('H2', 'N2', 'NH3')),
                  pressure_Pa = 1e5,
                  temperature_K = 300.0,
                  molar_composition = np.array([0.7, 0.3, 0.0]),
                  molar_flow_mol_s = 1.0)
    
    print(feed.density_SI())

    



