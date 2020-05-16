## GCEN
GCEN++: A C++ Tool for Gene Co-Expression Networks Analysis and Function annotation

#### Implementation  
**RWR**  
$$
p_k = αp_0 + (1-α)Wp_{k - 1}
$$

Where  $$p_0$$ represents our initial, or prior, information on genes. $$W$$ is a stochastic matrix, that is, its columns sum to 1. The parameter $$α$$ describes the trade-off between prior information and network smoothing. Repeated iteration of this equation converges to a steady-state.  
