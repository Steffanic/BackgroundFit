import numpy as np

def reaction_plane_resolution(n:int, j:int) -> float:
    if n==2 and j==2:
        return 0.8
    if n==3 and j==3:
        return 0.6
    if n==4 and j==2:
        return 0.4

def reaction_plane_correlations(n:int, j:int) -> float:
    if n in [2, 4, 6] and j==2:
        return 1
    if n in [3,5] and j==2:
        return 0

    

def background_level( delta_phi:float
                    , params:dict
                    , reaction_plane_order:int
                    , azimuth_region_width:float
                    , azimuth_region_center:float):
    r''' Calculates background level returns \tilde{B}

        From the paper: https://arxiv.org/pdf/1802.01668.pdf. This function computes the value of equation A.30.

        Keyword Arguments:
        delta_phi -- phi_t - phi_a
        params -- Dictionary containing values for {"v_1t":, "v_2t":, "v_3t":, "v_4t":, "v_5t":, "v_6t", "v_1a":, "v_2a":, "v_3a":, "v_4a":, "v_5a":, "v_6a":, "N^t":, "N^a"}
            v_n[t|a] -- The v_n for trigger or associated particles.
            trigger_yield -- number of triggers in this azimuthal region
            associated_yield -- number of associated particles in this azimuthal region
        reaction_plane_order -- the reaction plane order that correlations will be calculated relative to. j in the literature. [Charles/Christine help me out here]
        azimuth_region_width -- the width of this azimuthal region, usually denoted as c in the literature
        azimuth_region_center -- the center of this azimuthal region, usually denoted as phi_s in the literature
    '''

    TRUNCATION_ORDER = 3

    pre_factor = (params["N^t"] * params["N^a"] * reaction_plane_order * azimuth_region_width) / ( 2 * (np.pi)**2 )

    sum_factor = lambda k: (params[f"v_{reaction_plane_order*k}t"] \
                            * 1 / (reaction_plane_order * k * azimuth_region_width)) \
                            * np.sin( reaction_plane_order * k * azimuth_region_width) \
                            * reaction_plane_resolution(reaction_plane_order * k, reaction_plane_order) \
                            * reaction_plane_correlations(reaction_plane_order * k, reaction_plane_order) \
                            * np.cos(reaction_plane_order * k * azimuth_region_center)

    summation = np.sum([sum_factor(k) for k in range(1, TRUNCATION_ORDER+1)])

    return pre_factor*(1 + 2*summation)


                    
