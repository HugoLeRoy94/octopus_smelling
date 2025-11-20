def Popen(c,epsilon,Ko,Kc):
    """
    epsilon: energy difference between the open and closed state, without ligand. typically <0
    Kc: ligand dissociation constant when channel in closed state
    Ko: open-state dissociation constant : Ko_ = Ko exp(beta -epsilon/5)
    c: ligand concentration
    """
    c = c + 1e-18
    return (1+c/Ko)**5/((1+c/Ko)**5+np.exp(-epsilon)*(1+c/Kc)**5)
def Pmin(epsilon):
    return 1/(1+np.exp(-epsilon))
def Pmax(epsilon,Ko,Kc):
    return 1/(1+(Ko/Kc)**5*np.exp(-epsilon))
def PopenNorm(c,epsilon,Ko,Kc):
    return (Popen(c,epsilon,Ko,Kc) - Pmin(epsilon))/(Pmax(epsilon,Ko,Kc) - Pmin(epsilon))

def PmaxHetero(epsilon,Kos,Kcs):
    return 1/(1+np.prod([Ko/Kc for Ko,Kc in zip(Kos,Kcs)]) * np.exp(-epsilon))
def PoHetero(c,epsilon,Kos,Kcs):
    """
    Kos and Kcs are both list of Kc, and Ko. reapeating values will increase the number of monomer of a type.
    """
    return np.prod([1 + c/Ko for Ko in Kos],axis=0)/(np.prod([1+c/Ko for Ko in Kos],axis=0) + np.exp(-epsilon)*np.prod([1+c/Kc for Kc in Kcs],axis=0))
def PoHeteroNorm_vectorized(c, epsilon, Kos, Kcs):
    """
    Vectorized calculation of the Normalized Open Probability for Hetero-pentamers.
    Uses the MWC state ratio form: P = 1 / (1 + L_iso * R)
    """
    # 1. Isomerization constant L (Closed/Open equilibrium at c=0)
    # epsilon shape: (L, Nr, 1, 1) -> L_iso: (L, Nr, 1)
    L_iso = np.exp(-epsilon).squeeze(-1)
    
    # 2. State Ratios R(c) = Prod_subunits [ (1+c/Kc) / (1+c/Ko) ]
    # Calculates the relative stability of Closed vs Open state for all c
    # Shape: (L, Nr, Nc)
    R_c = np.prod((1 + c / Kcs) / (1 + c / Kos), axis=2)
    
    # 3. Limiting Ratio R(inf) = Prod(Ko/Kc)
    # This defines the saturation limit as c -> infinity
    # Kos/Kcs shape: (L, Nr, 5, 1) -> prod axis 2 -> (L, Nr, 1)
    # We DO NOT squeeze the last dim so it broadcasts with L_iso (L, Nr, 1)
    R_inf = np.prod(Kos / Kcs, axis=2)
    
    # 4. Calculate Probabilities: P = 1 / (1 + L * R)
    P_c = 1.0 / (1.0 + L_iso * R_c)
    P_min = 1.0 / (1.0 + L_iso)       # At c=0, R=1, shape (L, Nr, 1)
    P_max = 1.0 / (1.0 + L_iso * R_inf) # Shape (L, Nr, 1)
    
    # 5. Normalize
    # Use np.divide with 'where' to handle cases where P_max == P_min (flat response)
    delta_P = P_max - P_min
    return np.divide(P_c - P_min, delta_P, where=(delta_P != 0), out=np.zeros_like(P_c))

def PoHeteroNorm(c,epsilon,Kos,Kcs):
    return (PoHetero(c,epsilon,Kos,Kcs) - Pmin(epsilon))/(PmaxHetero(epsilon,Kos,Kcs) - Pmin(epsilon))

def generate_discrete_curve(epsilon, Kos, Kcs, Nc, c_min, c_max, Na):
    """
    Generates a discrete activation curve from the continuous Popen model.

    Args:
        epsilon, Ko, Kc: Parameters for the Popen function.
        Nc (int): Number of concentration bins.
        c_min (float): Minimum concentration for binning.
        c_max (float): Maximum concentration for binning.
        Na (int): Number of activity bins.

    Returns:
        np.array: A 1D array of length Nc containing the discrete
                  activity curve (e.g., [0, 0, 1] for Na=2, Nc=3).
    """
    
    # --- Concentration Bins (Logarithmic) ---
    # Create Nc+1 bin edges from c_min to c_max
    c_edges = np.logspace(np.log10(c_min), np.log10(c_max), Nc + 1)
    
    # Get Nc bin representative centers (geometric mean)
    c_centers = np.sqrt(c_edges[:-1] * c_edges[1:])
    
    # --- Activity Bins (Linear) ---
    # Create Na+1 bin edges from 0 to 1
    a_edges = np.linspace(0, 1, Na + 1)
    a_centers = 0.5*(a_edges[1:] + a_edges[:-1])
    a_centers[0] = 0.
    a_centers[-1] = 1.
    
    # --- Discretization ---
    # 1. Get the continuous Popen value at each concentration center
    #continuous_values = PopenNorm(c_centers, epsilon, Ko, Kc)
    continuous_values = PoHeteroNorm(c_centers, epsilon, Kos, Kcs)
    
    # 2. Digitize these values into bins
    discrete_curve = np.digitize(continuous_values, a_edges)
    return c_edges,a_centers[discrete_curve-1]

def generate_activity_matrix(epsilon, Kos, Kcs, c_centers, a_edges):
    """
    Generates the discrete activity matrix for a set of ligands.
    
    Args:
        epsilon: (L, Nr, 5)
        Kos:     (L, Nr, 5)
        Kcs:     (L, Nr, 5)
        c_centers: (Nc,)
        a_edges:   (Na+1,)
        
    Returns:
        discrete_A: Integer array of shape (L, Nr, Nc)
    """
    L, Nr, penta = epsilon.shape
    Nc = len(c_centers)
    
    # --- 1. Prepare Input Shapes for Broadcasting ---
    
    # c: Add axes for L, Nr, and Subunit -> (1, 1, 1, Nc)
    c_grid = c_centers[np.newaxis, np.newaxis, np.newaxis, :]
    
    # Kos, Kcs: Add axis for Concentration -> (L, Nr, 5, 1)
    Kos_grid = Kos[:, :, :, np.newaxis]
    Kcs_grid = Kcs[:, :, :, np.newaxis]
    
    # epsilon: All 5 subunits share the same epsilon for the complex global transition.
    # We take the first value and add axes -> (L, Nr, 1, 1)
    eps_grid = epsilon[:, :, 0][:, :, np.newaxis, np.newaxis]
    
    # --- 2. Calculate Continuous Response ---
    continuous_A = PoHeteroNorm_vectorized(c_grid, eps_grid, Kos_grid, Kcs_grid)
    
    # --- 3. Discretize ---
    # Clip to ensure numerical stability doesn't push values < 0 or > 1
    continuous_A = np.clip(continuous_A, 0, 1)
    
    # Digitize using inner edges as thresholds
    discrete_A = np.digitize(continuous_A, a_edges[1:-1])
    
    return discrete_A