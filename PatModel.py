import numpy as np

TRUNCATION_ORDER = 2


def reaction_plane_resolution(n:int, j:int) -> float:
    if n==2 and j==2:
        return 0.8
    if n==3 and j==3:
        return 0.6
    if n==4 and j==2:
        return 0.4
    else:
        return 0

def reaction_plane_correlations(n:int, m:int, j:int) -> float:
    if n in [2, 4, 6] and m in [2, 4, 6] and j==2:
        return 1
    if (n in [3,5] or m in [3,5]) and j==2:
        return 0
    else:
        return 0

    

def background_level( params:dict
                    , reaction_plane_order:int
                    , azimuth_region_width:float
                    , azimuth_region_center:float):
    r''' Calculates background level returns \tilde{B}

        From the paper: https://arxiv.org/pdf/1802.01668.pdf. This function computes the value of equation A.31.1.

        Keyword Arguments:
        params -- Dictionary containing values for {"v_1t":, "v_2t":, "v_3t":, "v_4t":, "v_5t":, "v_6t", "v_1a":, "v_2a":, "v_3a":, "v_4a":, "v_5a":, "v_6a":, "N^t":, "N^a"}
            v_n[t|a] -- The v_n for trigger or associated particles.
            trigger_yield -- number of triggers in this azimuthal region
            associated_yield -- number of associated particles in this azimuthal region
        reaction_plane_order -- the reaction plane order that correlations will be calculated relative to. j in the literature. [Charles/Christine help me out here]
        azimuth_region_width -- the width of this azimuthal region, usually denoted as c in the literature
        azimuth_region_center -- the center of this azimuthal region, usually denoted as phi_s in the literature
    '''

    print(f"\n\n\n\n\n\n\n{params['N^t'] * params['N^a']}")
    pre_factor = (params["N^t"] * params["N^a"] * reaction_plane_order * azimuth_region_width) / ( 2 * (np.pi)**2 )

    sum_factor = lambda k: params[f"v_{int(reaction_plane_order*k)}t"] \
                            * (1 / (reaction_plane_order * k * azimuth_region_width)) \
                            * np.sin( reaction_plane_order * k * azimuth_region_width) \
                            * reaction_plane_resolution(reaction_plane_order * k, reaction_plane_order) \
                            * reaction_plane_correlations(reaction_plane_order * k, 0, reaction_plane_order) \
                            * np.cos(reaction_plane_order * k * azimuth_region_center)

    summation = np.sum([sum_factor(k) for k in range(1, TRUNCATION_ORDER+1)])

    return pre_factor*(1 + 2*summation)


def effective_v_nt( params:dict
                  , flow_order:int
                  , reaction_plane_order:int
                  , azimuth_region_width:float
                  , azimuth_region_center:float ):
    r''' Calculates effective v_n for the trigger particles, returns \tilde{v_n^t}

        From the paper: https://arxiv.org/pdf/1802.01668.pdf. This function computes the value of equation A.31.2

        Keyword Arguments:
        params -- Dictionary containing values for {"v_1t":, "v_2t":, "v_3t":, "v_4t":, "v_5t":, "v_6t", "v_1a":, "v_2a":, "v_3a":, "v_4a":, "v_5a":, "v_6a":, "N^t":, "N^a":}
            v_n[t|a] -- The v_n for trigger or associated particles.
            trigger_yield -- number of triggers in this azimuthal region
            associated_yield -- number of associated particles in this azimuthal region
        flow_order -- the order of the effective flow coefficient to be computed, the n in v_n
        reaction_plane_order -- the reaction plane order that correlations will be calculated relative to. j in the literature. [Charles/Christine help me out here]
        azimuth_region_width -- the width of this azimuthal region, usually denoted as c in the literature
        azimuth_region_center -- the center of this azimuthal region, usually denoted as phi_s in the literature
    '''
    non_sum_term = 0
    if flow_order / reaction_plane_order == flow_order // reaction_plane_order:
        non_sum_term = 1 / (flow_order * azimuth_region_width) \
                    * np.sin(flow_order * azimuth_region_width) \
                    * reaction_plane_resolution(flow_order, reaction_plane_order) \
                    * reaction_plane_correlations(flow_order, 0, reaction_plane_order) \
                    * np.cos(flow_order * azimuth_region_center)
    
    summand_prefactor = lambda k: params[f"v_{reaction_plane_order * k + flow_order}t"] \
            * reaction_plane_correlations(np.abs(reaction_plane_order * k + flow_order), flow_order, reaction_plane_order) \
                +  (0 if np.abs(reaction_plane_order * k - flow_order)==0 else params[f"v_{np.abs(reaction_plane_order * k - flow_order)}t"]) \
            * reaction_plane_correlations(np.abs(reaction_plane_order * k - flow_order), flow_order, reaction_plane_order) \
                    
    summand_postfactor = lambda k: (np.sin(reaction_plane_order * k * azimuth_region_width) \
                                * np.cos(reaction_plane_order * k * azimuth_region_center) \
                                * reaction_plane_resolution(reaction_plane_order * k, reaction_plane_order) )\
                                / (reaction_plane_order * k * azimuth_region_width)

    summand = lambda k: summand_prefactor(k) * summand_postfactor(k)

    summation = np.sum([summand(k) for k in range(1, TRUNCATION_ORDER + 1)])

    numerator = params[f"v_{flow_order}t"] + non_sum_term + summation

    denominator_sum_factor = lambda k: (params[f"v_{reaction_plane_order*k}t"] \
                            * 1 / (reaction_plane_order * k * azimuth_region_width)) \
                            * np.sin( flow_order * azimuth_region_width) \
                            * reaction_plane_resolution(reaction_plane_order * k, reaction_plane_order) \
                            * reaction_plane_correlations(reaction_plane_order * k, 0, reaction_plane_order) \
                            * np.cos(reaction_plane_order * k * azimuth_region_center)

    denominator_summation = np.sum([denominator_sum_factor(k) for k in range(1, TRUNCATION_ORDER+1)])

    denominator = 1 + 2 * denominator_summation

    return numerator / denominator

def effective_w_nt( params:dict
                  , flow_order:int
                  , reaction_plane_order:int
                  , azimuth_region_width:float
                  , azimuth_region_center:float ):
    r''' Calculates effective w_n for the trigger particles, returns \tilde{w_n^t}

        From the paper: https://arxiv.org/pdf/1802.01668.pdf. This function computes the value of equation A.31.3

        Keyword Arguments:
        params -- Dictionary containing values for {"v_1t":, "v_2t":, "v_3t":, "v_4t":, "v_5t":, "v_6t", "v_1a":, "v_2a":, "v_3a":, "v_4a":, "v_5a":, "v_6a":, "N^t":, "N^a"}
            v_n[t|a] -- The v_n for trigger or associated particles.
            trigger_yield -- number of triggers in this azimuthal region
            associated_yield -- number of associated particles in this azimuthal region
        flow_order -- the order of the effective flow coefficient to be computed, the n in v_n
        reaction_plane_order -- the reaction plane order that correlations will be calculated relative to. j in the literature. [Charles/Christine help me out here]
        azimuth_region_width -- the width of this azimuthal region, usually denoted as c in the literature
        azimuth_region_center -- the center of this azimuthal region, usually denoted as phi_s in the literature
    '''
    non_sum_term = 0
    if flow_order / reaction_plane_order == flow_order // reaction_plane_order:
        non_sum_term = 1 / (flow_order * azimuth_region_width) \
                    * np.sin(flow_order * azimuth_region_width) \
                    * reaction_plane_resolution(flow_order, reaction_plane_order) \
                    * reaction_plane_correlations(flow_order, 0, reaction_plane_order) \
                    * np.sin(flow_order * azimuth_region_center)
    
    summand_prefactor = lambda k: params[f"v_{reaction_plane_order * k + flow_order}t"] \
            * reaction_plane_correlations(np.abs(reaction_plane_order * k + flow_order), flow_order, reaction_plane_order) \
                +  (0 if np.abs(reaction_plane_order * k - flow_order)==0 else params[f"v_{np.abs(reaction_plane_order * k - flow_order)}t"]) \
            * reaction_plane_correlations(np.abs(reaction_plane_order * k - flow_order), flow_order, reaction_plane_order) \
                    
    summand_postfactor = lambda k: (np.sin(reaction_plane_order * k * azimuth_region_width) \
                                * np.cos(reaction_plane_order * k * azimuth_region_center) \
                                * reaction_plane_resolution(reaction_plane_order * k, reaction_plane_order) )\
                                / (reaction_plane_order * k * azimuth_region_width)

    summand = lambda k: summand_prefactor(k) * summand_postfactor(k)

    summation = np.sum([summand(k) for k in range(1, TRUNCATION_ORDER + 1)])

    numerator = non_sum_term + summation

    denominator_sum_factor = lambda k: (params[f"v_{reaction_plane_order*k}t"] \
                            * 1 / (reaction_plane_order * k * azimuth_region_width)) \
                            * np.sin( flow_order * azimuth_region_width) \
                            * reaction_plane_resolution(reaction_plane_order * k, reaction_plane_order) \
                            * reaction_plane_correlations(reaction_plane_order * k, 0, reaction_plane_order) \
                            * np.cos(reaction_plane_order * k * azimuth_region_center)

    denominator_summation = np.sum([denominator_sum_factor(k) for k in range(1, TRUNCATION_ORDER+1)])

    denominator = 1 + 2 * denominator_summation

    return numerator / denominator


def background(   delta_phi:float
                , params:dict
                , reaction_plane_order:int
                , azimuth_region_width:float
                , azimuth_region_center:float ):
    r''' Calculates background distribution over delta-phi, returns B(\Delta \phi)

        From the paper: https://arxiv.org/pdf/1802.01668.pdf. This function computes the value of equation A.30

        Keyword Arguments:
        delta_phi -- Azimuthal correlation observable
        params -- Dictionary containing values for {"v_1t":, "v_2t":, "v_3t":, "v_4t":, "v_5t":, "v_6t", "v_1a":, "v_2a":, "v_3a":, "v_4a":, "v_5a":, "v_6a":, "N^t":, "N^a"}
            v_n[t|a] -- The v_n for trigger or associated particles.
            trigger_yield -- number of triggers in this azimuthal region
            associated_yield -- number of associated particles in this azimuthal region
        reaction_plane_order -- the reaction plane order that correlations will be calculated relative to. j in the literature. [Charles/Christine help me out here]
        azimuth_region_width -- the width of this azimuthal region, usually denoted as c in the literature
        azimuth_region_center -- the center of this azimuthal region, usually denoted as phi_s in the literature
    '''
    B_tilde = background_level(params, reaction_plane_order, azimuth_region_width, azimuth_region_center)

    summand = lambda n: params[f"v_{n}a"] \
                        * (effective_v_nt(params, n, reaction_plane_order, azimuth_region_width, azimuth_region_center) \
                            *np.cos(n * delta_phi) + effective_w_nt(params, n, reaction_plane_order, azimuth_region_width, azimuth_region_center) \
                            *np.sin(n * delta_phi))

    summation = np.sum([summand(n) for n in range(1, TRUNCATION_ORDER+1)], axis=0)

    return B_tilde * (1 + 2 * summation)