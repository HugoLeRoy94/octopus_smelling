import numpy as np

def Pmin(epsilon):
    return 1/(1+np.exp(-epsilon))
def Pmax(epsilon,Ko,Kc):
    return 1/(1+(Ko/Kc)**5*np.exp(-epsilon))
def PopenNorm(c,epsilon,Ko,Kc):
    return (Popen(c,epsilon,Ko,Kc) - Pmin(epsilon))/(Pmax(epsilon,Ko,Kc) - Pmin(epsilon))

def PmaxHetero(epsilon, Kos, Kcs):
    # Using list comprehension for element-wise division if inputs are lists/arrays
    ratio_prod = np.prod([Ko/Kc for Ko, Kc in zip(Kos, Kcs)])
    return 1/(1 + ratio_prod * np.exp(-epsilon))
def PoHetero(c, epsilon, Kos, Kcs):
    """
    Calculates absolute open probability for a heteropentamer.
    c: float or np.array of ligand concentrations
    epsilon: float, total gating energy
    Kos, Kcs: lists/arrays of length 5 containing dissociation constants
    """
    # Ensure c is handled correctly if it's a scalar or array
    term_o = np.prod([1 + c/Ko for Ko in Kos], axis=0)
    term_c = np.prod([1 + c/Kc for Kc in Kcs], axis=0)
    term_eps = np.sum(epsilon,axis=0)
    return term_o / (term_o + np.exp(-term_eps) * term_c)
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

def PoHeteroNorm(c, epsilon, Kos, Kcs):
    p_abs = PoHetero(c, epsilon, Kos, Kcs)
    p_min = Pmin(epsilon)
    p_max = PmaxHetero(epsilon, Kos, Kcs)
    # Normalize
    return (p_abs - p_min) / (p_max - p_min)

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

def generate_continuous_A(epsilon,Kos,Kcs,c_centers,a_edges,indices_to_avg=None):
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

    if indices_to_avg is not None:
        all_indices = np.arange(continuous_A.shape[1])
        remaining_indices = np.delete(all_indices, indices_to_avg) # [0, 1]
        # 2. Extract the "untouched" rows
        untouched_part = continuous_A[:, remaining_indices, :]
        # 3. Calculate the mean of the "target" rows
        # We use axis=1 to average across the row dimension
        avg_part = np.mean(continuous_A[:, indices_to_avg, :], axis=1, keepdims=True)
        # 4. Combine them
        continuous_A = np.concatenate([untouched_part, avg_part], axis=1)

    return continuous_A
    

def generate_activity_matrix(epsilon, Kos, Kcs, c_centers, a_edges,indices_to_avg=None):
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
    
    # --- 2. Calculate Continuous Response ---
    continuous_A = generate_continuous_A(epsilon,Kos,Kcs,c_centers,a_edges,indices_to_avg)
    
    # --- 3. Discretize ---
    # Clip to ensure numerical stability doesn't push values < 0 or > 1
    continuous_A = np.clip(continuous_A, 0, 1)

    # Digitize using inner edges as thresholds
    discrete_A = np.digitize(continuous_A, a_edges[1:-1])
    
    return discrete_A

import itertools
import random

def generate_Ks(L=100, n_types=2, n_subunits=5, n_hetero_sample='all',Kc_Amp=10.,eps_Amp=-10.):
    """
    Generates Ks, Kos, and Epsilon for N types with flexible heteromer selection.
    
    Args:
        L (int): Number of simulations.
        n_types (int): Number of base subunit types.
        n_subunits (int): Number of subunits per receptor.
        n_hetero_sample (int or str): 
            - 'all': Generate all possible hetero-combinations.
            - 0: Generate only Homomers.
            - int > 0: Generate all Homomers + n randomly selected Heteromers.
        Kc_Amp (float): amplitude of the uniform distribution from which the highest Kc is drawn
        eps_Amp (float): amplitude of the uniform distribution from which epsilon is drawn
    """
    
    # 1. Generate all base combinations
    all_combos_iter = itertools.combinations_with_replacement(range(n_types), n_subunits)
    all_combos = list(all_combos_iter)
    
    # 2. Separate Homomers and Heteromers
    homomers = [c for c in all_combos if len(set(c)) == 1]
    heteromers = [c for c in all_combos if len(set(c)) > 1]
    
    # 3. Select Heteromers based on input
    if n_hetero_sample == 'all':
        selected_hetero = heteromers
    elif isinstance(n_hetero_sample, int):
        # Clamp the number to the maximum available heteromers to avoid errors
        n_to_pick = min(n_hetero_sample, len(heteromers))
        if n_to_pick > 0:
            selected_hetero = random.sample(heteromers, n_to_pick)
            # Optional: Sort them to keep output somewhat deterministic/tidy
            selected_hetero.sort() 
        else:
            selected_hetero = []
    else:
        # Fallback if bad input provided
        print(f"Warning: Unknown n_hetero_sample '{n_hetero_sample}'. returning only Homomers.")
        selected_hetero = []

    # 4. Final Receptor List: Homomers first, then selected Heteromers
    # Sorting homomers ensures Index 0 is Type 0, Index 1 is Type 1, etc.
    homomers.sort() 
    final_combos = homomers + selected_hetero
    
    Nr = len(final_combos)
    
    # 5. Pre-allocate arrays
    Kcs = np.zeros((L, Nr, n_subunits), dtype=float)
    Kos = np.zeros((L, Nr, n_subunits), dtype=float)
    epsilon = np.zeros((L, Nr, n_subunits), dtype=float)
    
    # 6. Generate Parameters
    for l in range(L):
        # --- A. Generate Base Parameters for N types ---
        base_Kc = []
        base_Ko = []
        base_eps = []
        
        prev_Kc = 0
        prev_Ko = None 
        
        for i in range(n_types):
            # Kc increases
            current_Kc = (prev_Kc if i > 0 else 0) + np.random.random() * Kc_Amp
            base_Kc.append(current_Kc)
            
            # Ko decreases (and Ko <= Kc)
            if i == 0:
                current_Ko = np.random.random() * current_Kc
            else:
                current_Ko = np.random.random() * prev_Ko
            base_Ko.append(current_Ko)
            
            base_eps.append(np.random.random() * eps_Amp)
            
            prev_Kc = current_Kc
            prev_Ko = current_Ko

        # --- B. Map to Receptor Structure ---
        for r_idx, combo_indices in enumerate(final_combos):
            Kcs[l, r_idx, :] = [base_Kc[i] for i in combo_indices]
            Kos[l, r_idx, :] = [base_Ko[i] for i in combo_indices]
            
            combo_eps_values = [base_eps[i] for i in combo_indices]
            epsilon[l, r_idx, :] = [np.mean(combo_eps_values)] * n_subunits
            
    return Kcs, Kos, epsilon

#def generate_Ks(L=100, n_types=2, Nhetero=3):
#    """
#    Generates Ks, Kos, and Epsilon for N types.
#    
#    Args:
#        L (int): Number of simulations.
#        n_types (int): Number of base subunit types.
#        n_subunits (int): Number of subunits per receptor.
#        include_hetero (bool): If True, generates all mixed combinations. 
#                               If False, generates only the N homomers.
#    """
#    n_subunits=5
#    # 1. Define Receptors to build
#    if Nhetero==0:
#        all_combos = [(i,) * n_subunits for i in range(n_types)]
#    else:
#        all_combos = list(itertools.combinations_with_replacement(range(n_types), n_subunits))
#        # Sort: Homomers first, then Heteromers
#        all_combos.sort(key=lambda x: (len(set(x)) > 1, x))
#        
#            
#    if include_hetero:
#        # Generate all combinations (homo + hetero)
#        
#        
#    elif Nhetero>=all_combos.__len__():
#        # Generate ONLY Homomers [(0,0,0,0,0), (1,1,1,1,1)...]
#        all_combos = [(i,) * n_subunits for i in range(n_types)]
#    
#    Nr = len(all_combos)
#    
#    # 2. Pre-allocate arrays
#    Kcs = np.zeros((L, Nr, n_subunits), dtype=float)
#    Kos = np.zeros((L, Nr, n_subunits), dtype=float)
#    epsilon = np.zeros((L, Nr, n_subunits), dtype=float)
#    
#    for l in range(L):
#        # --- A. Generate Base Parameters for N types ---
#        # (This remains the same regardless of include_hetero)
#        base_Kc = []
#        base_Ko = []
#        base_eps = []
#        
#        prev_Kc = 0
#        prev_Ko = None 
#        
#        for i in range(n_types):
#            # Kc increases with type index
#            current_Kc = (prev_Kc if i > 0 else 0) + np.random.random() * 10
#            base_Kc.append(current_Kc)
#            
#            # Ko decreases with type index (and is < Kc)
#            if i == 0:
#                current_Ko = np.random.random() * current_Kc
#            else:
#                current_Ko = np.random.random() * prev_Ko
#            base_Ko.append(current_Ko)
#            
#            base_eps.append(np.random.random() * -10)
#            
#            prev_Kc = current_Kc
#            prev_Ko = current_Ko
#
#        # --- B. Map Base Params to Receptor Structure ---
#        for r_idx, combo_indices in enumerate(all_combos):
#            # Map parameters based on the indices in the combination
#            Kcs[l, r_idx, :] = [base_Kc[i] for i in combo_indices]
#            Kos[l, r_idx, :] = [base_Ko[i] for i in combo_indices]
#            
#            # Epsilon is the average of the subunits
#            combo_eps_values = [base_eps[i] for i in combo_indices]
#            epsilon[l, r_idx, :] = [np.mean(combo_eps_values)] * n_subunits
#            
#    return Kcs, Kos, epsilon

#def generate_Ks_3hetero():
#    # draw random parameters # 0 is homo1, 2 is homo2, 3 is hetero
#    # Ko << Kc
#    L = 100
#    Nr = 3
#    Kcs_3hetero = np.zeros((L,Nr,5),dtype=float)
#    Kos_3hetero = np.zeros((L,Nr,5),dtype=float)
#    epsilon_3hetero = np.zeros((L,Nr,5),dtype=float)
#    for l in range(L):
#        # select Kc1 \in [0, 10], Ko1 \in [0,Kc1] -> Ko1<=Kc1
#        Kc1 = np.random.random() * 10
#        Ko1 = np.random.random() * Kc1
#        # select Kc2 \in [Kc1, Kc1+10] Kc2>=Kc1 & Ko2 \in [0,Ko1] Ko2<=Ko1 | Ko2<=Kc2
#        Kc2 = Kc1 + np.random.random()*10
#        Ko2 = np.random.random()*Ko1
#        # generate three receptors:
#        Kcs_3hetero[l,0,:] = [Kc1] * 5 # Homo-1
#        Kcs_3hetero[l,1,:] = [Kc1] * 3 + [Kc2] * 2 # Hetero
#        Kcs_3hetero[l,2,:] = [Kc2] * 5 # Homo-2
#
#        Kos_3hetero[l,0,:] = [Ko1] * 5 # Homo-1
#        Kos_3hetero[l,1,:] = [Ko1] * 3 + [Ko2] * 2 # Hetero
#        Kos_3hetero[l,2,:] = [Ko2] * 5 # Homo-2
#
#        # epsilon \in [-10,0] -> epsilon <<1    
#        eps1 = np.random.random() * -10
#        eps2 = np.random.random() * -10
#        epsilon_3hetero[l,0,:] = [eps1] * 5
#        epsilon_3hetero[l,2,:] = [eps2] * 5
#        epsilon_3hetero[l,1,:] = [eps1 * 3./5. + eps2*2./5.]*5
#    return Kcs_3hetero,Kos_3hetero,epsilon_3hetero