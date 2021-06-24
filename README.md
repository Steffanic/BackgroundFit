# BackgroundFit
This code implements the reaction plane method proposed by C. Nattrass et al. in Phys. Rev. C 93, 044915.

This codes uses as input a ROOT file that has to contain histograms of Dphi distribution in Signal and Background regions. 
Check the example fInput.root for reference.  

You need to install iminuit:
http://iminuit.readthedocs.io/en/latest/

## Steffanic Fork
This fork adds a new non-symmetrized version of the fit function in a generally extensible form as more data become available. It is implemented in PatModel.py *name should be changed*. 

Here's a code snippet on how to use the background function. 

    from PatModel import background
    import numpy as np

    params = {"v_1t":0, "v_2t":0.1, "v_3t":0.1, "v_4t":0.02, "v_5t":0, "v_6t":0, "v_1a":0, "v_2a":0.1, "v_3a":0.1, "v_4a":0.02, "v_5a":0, "v_6a":0, "N^t":np.pi/(2*np.pi/6), "N^a":2*np.pi}
    reaction_plane_order=2
    azimuth_width = np.pi / 6
    azimuth_center = 0
    delta_phi = np.linspace(-np.pi/2, 3*np.pi/2, 100)

    BG = background(delta_phi, params, reaction_plane_order, azimuth_width, azimuth_center)
  
From there, BG is a numpy array of 100 bin values, background counts.
  
There are two important global variables: TRUNCATION_ORDER and MAX_FLOW_ORDER

TRUNCATION_ORDER: This is the order that the summations inside \tilde{v}_n and \tilde{w}_n will terminate at.

MAX_FLOW_ORDER: This is the highest order flow coefficient that we want to fit to.
  
  
  
