import numpy as np


def background_level( delta_phi:float
                    , params:dict
                    , trigger_yield:int
                    , associated_yield:int
                    , reaction_plane_order:int
                    , azimuth_region_width:float
                    , azimuth_region_center:float):
                    r''' Calculates background level \tilde{B} and effective \tilde{v}^t_n, \tilde{w}^t_n; returns B(\Delta \phi)

                        From the paper: https://arxiv.org/pdf/1802.01668.pdf. This function computes the value of equation A.30.

                        Keyword Arguments:
                        delta_phi -- phi_t - phi_a
                        params -- Dictionary containing values for {"v_1t":, "v_2t":, "v_3t":, "v_4t":, "v_6t", "v_1a":, "v_2a":, "v_3a":, "v_4a":, "v_6a"}
                        trigger_yield -- number of triggers in this azimuthal region
                        associated_yield -- number of associated particles in this azimuthal region
                        reaction_plane_order -- the reaction plane order that correlations will be calculated relative to. j in the literature. [Charles/Christine help me out here]
                        azimuth_region_width -- the width of this azimuthal region, usually denoted as c in the literature
                        azimuth_region_center -- the center of this azimuthal region, usually denoted as phi_s in the literature
                    '''

                    
