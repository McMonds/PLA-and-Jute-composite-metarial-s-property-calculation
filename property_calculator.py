
import math
import numpy as np # For potential array operations later

# --- Constants ---
PI = math.pi
STEFAN_BOLTZMANN = 5.67e-8 # W/(m^2 K^4)
GRAVITY = 9.81 # m/s^2

# --- Helper Functions (Placeholders - Need Actual Models) ---

def get_pla_properties(pla_grade_name):
    # TODO: Implement lookup or require detailed input
    # Example placeholder values - REPLACE WITH ACTUAL DATA
    print(f"--- WARNING: Using placeholder properties for PLA grade: {pla_grade_name} ---")
    props = {
        'rho_pla': 1240, # kg/m^3
        'Cp_pla': 1800, # J/(kg K) - Assume constant for simplicity now
        'k_pla': 0.13, # W/(m K) - Assume constant
        'Tg_pla': 60, # deg C
        'Tm_pla': 155, # deg C
        'K_pla': 15000, # Pa.s^n (Power Law Consistency Index - NEEDS REAL DATA)
        'n_pla': 0.45, # (Power Law Index - NEEDS REAL DATA)
        'Ea_pla': 50000, # J/mol (Arrhenius Activation Energy - NEEDS REAL DATA)
        'yield_strength_pla': lambda T_C: max(1e6, 50e6 * math.exp(-0.05 * max(0, T_C - 25))) # VERY ROUGH EXAMPLE MODEL (Pa vs deg C)
    }
    return props

def get_jute_properties(jute_type):
    # TODO: Implement lookup or require detailed input
    print(f"--- WARNING: Using placeholder properties for Jute type: {jute_type} ---")
    props = {
        'rho_jute': 1450, # kg/m^3
        'Cp_jute': 1400, # J/(kg K)
        'k_jute': 0.05, # W/(m K)
        'T_degradation_jute': 220 # deg C
    }
    return props

def calculate_composite_properties(pla_props, jute_props, jute_loading_wt):
    print("\n--- Calculating Composite Properties ---")
    # Convert wt% to volume fraction
    wt_pla = 1.0 - jute_loading_wt
    vol_pla_rel = wt_pla / pla_props['rho_pla']
    vol_jute_rel = jute_loading_wt / jute_props['rho_jute']
    phi_jute = vol_jute_rel / (vol_pla_rel + vol_jute_rel)
    print(f"Jute Volume Fraction (phi_jute): {phi_jute:.3f}")

    # Rule of Mixtures / Simple Estimates (PLACEHOLDERS)
    rho_comp = pla_props['rho_pla'] * (1 - phi_jute) + jute_props['rho_jute'] * phi_jute
    Cp_comp = pla_props['Cp_pla'] * wt_pla + jute_props['Cp_jute'] * jute_loading_wt
    # Simple conductivity estimate (often poor for composites)
    k_comp = pla_props['k_pla'] * (1 - phi_jute) + jute_props['k_jute'] * phi_jute
    Tg_comp = pla_props['Tg_pla'] # Assume unchanged for now
    Tm_comp = pla_props['Tm_pla'] # Assume unchanged for now

    # --- CRITICAL: Composite Viscosity Estimation ---
    # This is highly complex. Using PLA values as a starting point (BAD ASSUMPTION!)
    # Real model needed (e.g., Mooney, Krieger-Dougherty, or measured data)
    K_comp = pla_props['K_pla'] * (1 + 2.5 * phi_jute + 10 * phi_jute**2) # Example filler effect (NEEDS VALIDATION)
    n_comp = pla_props['n_pla'] # Assume n unchanged (Often not true)
    Ea_comp = pla_props['Ea_pla']
    print(f"--- WARNING: Composite viscosity params are ROUGH ESTIMATES (K={K_comp:.1f}, n={n_comp:.2f}) ---")

    # Composite Yield Strength (Assume matrix dominates - NEEDS BETTER MODEL/DATA)
    yield_strength_comp = pla_props['yield_strength_pla']
    print(f"--- WARNING: Using PLA yield strength model for composite ---")

    comp_props = {
        'phi_jute': phi_jute, 'rho_comp': rho_comp, 'Cp_comp': Cp_comp, 'k_comp': k_comp,
        'Tg_comp': Tg_comp, 'Tm_comp': Tm_comp, 'K_comp': K_comp, 'n_comp': n_comp,
        'Ea_comp': Ea_comp, 'yield_strength_comp': yield_strength_comp,
        'emissivity_comp': 0.9 # Assume
    }
    print(f"Calculated rho_comp: {rho_comp:.1f} kg/m^3")
    print(f"Calculated Cp_comp: {Cp_comp:.1f} J/(kg K)")
    print(f"Calculated k_comp: {k_comp:.3f} W/(m K)")
    return comp_props

def get_viscosity(T_K, shear_rate, K, n, Ea, R=8.314):
    """Calculates viscosity using Power Law and Arrhenius models."""
    if shear_rate < 1e-6: shear_rate = 1e-6 # Avoid division by zero / issues at low shear
    # Arrhenius temperature dependence for K
    # We need a reference temperature T_ref_K where K was defined. Assume 200C = 473.15K for now
    T_ref_K = 200 + 273.15
    K_eff = K * math.exp((Ea / R) * (1 / T_K - 1 / T_ref_K))
    eta = K_eff * (shear_rate**(n - 1))
    return eta

def estimate_extrusion(rpm, screw_D, channel_H, die_D, die_L, T_zones_C, comp_props):
    print("\n--- Estimating Extrusion Process ---")
    # --- Simplified Flow Rate Estimation ---
    # Very rough model: Q proportional to RPM and screw geometry
    # Needs calibration factor based on actual measurements
    CALIBRATION_FACTOR = 0.5 # COMPLETE GUESS - Needs real data
    Q_drag = CALIBRATION_FACTOR * PI**2 * screw_D**2 * channel_H * (rpm / 60) * math.sin(math.radians(17.7)) * math.cos(math.radians(17.7)) # Simplified drag flow
    Q = Q_drag # Assume pressure flow is negligible for now (ANOTHER GUESS)
    print(f"Estimated Volumetric Flow Rate (Q): {Q * 1e6:.2f} cm^3/s (HIGHLY APPROXIMATE!)")
    if Q <= 0:
        print("ERROR: Estimated flow rate is zero or negative. Check RPM/Geometry.")
        return None

    # --- Die Calculations ---
    die_R = die_D / 2.0
    avg_melt_T_K = (T_zones_C[-1] + T_zones_C[-2]) / 2.0 + 273.15 # Guess initial melt temp = avg of last two zones
    shear_rate_die = (4 * Q) / (PI * die_R**3) # Apparent wall shear rate
    print(f"Estimated Shear Rate in Die: {shear_rate_die:.1f} 1/s")

    # Iterate to find consistent melt temp including shear heating
    T_melt_exit_K = avg_melt_T_K
    print("Iterating for melt temperature (incl. shear heating):")
    for i in range(5): # Simple iteration
        eta_die = get_viscosity(T_melt_exit_K, shear_rate_die,
                                comp_props['K_comp'], comp_props['n_comp'], comp_props['Ea_comp'])
        tau_die = eta_die * shear_rate_die # Wall shear stress

        # Estimate shear heating in the die land
        power_shear_die = tau_die * shear_rate_die * (PI * die_R**2 * die_L) # Power ~ stress*rate*volume
        delta_T_shear_die = power_shear_die / (Q * comp_props['rho_comp'] * comp_props['Cp_comp'])

        # Assume exit temp is roughly die set temp + shear heating in die
        # A better model would integrate heat transfer and shear heating along the whole barrel
        T_melt_exit_K_new = T_zones_C[-1] + 273.15 + delta_T_shear_die
        print(f"  Iter {i+1}: Assumed T={T_melt_exit_K-273.15:.1f}C -> eta={eta_die:.1f} Pa.s -> DeltaT_shear={delta_T_shear_die:.1f}C -> New T={T_melt_exit_K_new-273.15:.1f}C")

        if abs(T_melt_exit_K_new - T_melt_exit_K) < 0.5:
            T_melt_exit_K = T_melt_exit_K_new
            break
        T_melt_exit_K = T_melt_exit_K_new
    else:
         print("--- WARNING: Melt temperature iteration did not fully converge. ---")

    print(f"Estimated Melt Exit Temperature: {T_melt_exit_K - 273.15:.1f} °C")
    print(f"Estimated Viscosity at Die Exit: {eta_die:.1f} Pa·s")

    # --- Pressure Drop Estimation (Simplified Power Law) ---
    n = comp_props['n_comp']
    K_eff_die = comp_props['K_comp'] * math.exp((comp_props['Ea_comp'] / 8.314) * (1 / T_melt_exit_K - 1 / (200+273.15)))
    try:
        delta_P_die = (2 * K_eff_die * die_L / die_R) * ((3*n + 1) / (n * PI * die_R**3) * Q)**n
        print(f"Estimated Pressure Drop across Die: {delta_P_die / 1e6:.2f} MPa")
    except OverflowError:
        print("--- WARNING: Pressure drop calculation overflowed (likely very high viscosity/flow) ---")
        delta_P_die = float('inf')

    return {'Q': Q, 'T_melt_exit_K': T_melt_exit_K, 'eta_die': eta_die, 'delta_P_die': delta_P_die}

def model_cooling_and_stress(Q, T_melt_exit_K, target_D, T_air_K, h_conv, L_max_hang, comp_props):
    print("\n--- Modeling Filament Cooling and Gravity Stress ---")
    filament_R = target_D / 2.0
    A_cross_section = PI * filament_R**2
    v_filament = Q / A_cross_section
    print(f"Filament Exit Velocity: {v_filament * 1000:.2f} mm/s")

    # --- Simplified Cooling Model (Lumped Capacitance along length) ---
    # Ignores radial temperature gradients - assumes filament cools uniformly
    # More accurate model needs Finite Difference Method for radial conduction
    perimeter = 2 * PI * filament_R
    rho = comp_props['rho_comp']
    Cp = comp_props['Cp_comp']
    emissivity = comp_props['emissivity_comp']
    Tg_K = comp_props['Tg_comp'] + 273.15

    z = 0.0 # Distance from die
    t = 0.0 # Time
    T_curr_K = T_melt_exit_K
    dt = 0.01 # Time step for simulation (s)
    z_solidify = -1.0 # Distance where T drops below Tg

    print("Simulating cooling along filament length:")
    print("  Dist (mm) | Time (s) | Temp (°C)")
    print(f"  {z*1000:9.1f} | {t:8.3f} | {T_curr_K-273.15:8.1f} (Start)")

    max_steps = int((L_max_hang / v_filament) / dt) + 100
    if max_steps > 5000:
        print("--- WARNING: Simulation might take time due to low velocity / long length ---")
        max_steps = 5000 # Limit steps

    for step in range(max_steps):
        if T_curr_K < T_air_K + 0.1: # Stop if close to ambient
             print("  Filament cooled to near ambient.")
             break

        T_surf_K = T_curr_K # Lumped capacitance assumption
        q_conv = h_conv * perimeter * (T_surf_K - T_air_K)
        q_rad = emissivity * STEFAN_BOLTZMANN * perimeter * (T_surf_K**4 - T_air_K**4)
        q_total = q_conv + q_rad # Heat loss rate per unit length (W/m)

        # Temperature change rate dT/dt = (Heat Loss Rate) / (mass_per_length * Cp)
        mass_per_length = rho * A_cross_section
        dT_dt = -q_total / (mass_per_length * Cp)

        T_curr_K += dT_dt * dt
        t += dt
        z = v_filament * t

        if step % (max(1, int(0.05 / dt))) == 0: # Print every ~0.05s
            print(f"  {z*1000:9.1f} | {t:8.3f} | {T_curr_K-273.15:8.1f}")

        if T_curr_K < Tg_K and z_solidify < 0:
            z_solidify = z
            print(f"  SOLIDIFIED BELOW Tg ({Tg_K-273.15:.1f}°C) at z = {z_solidify*1000:.1f} mm")

    if z_solidify < 0 :
         print(f"--- WARNING: Filament did not cool below Tg within {L_max_hang:.2f} m ---")
         z_solidify = L_max_hang # Consider it not solidified for stress check

    # --- Gravity Stress Check ---
    max_stress = rho * GRAVITY * L_max_hang
    print(f"\nMax Gravitational Stress (at L={L_max_hang:.1f}m): {max_stress/1e3:.2f} kPa")

    # Evaluate yield strength at the point it should be solid enough
    T_at_solidify_point = T_curr_K if z_solidify >= L_max_hang else Tg_K # Use Tg if solidified earlier
    # Need temp in C for the example yield strength function
    yield_strength_at_T = comp_props['yield_strength_comp'](T_at_solidify_point - 273.15)
    print(f"Estimated Yield Strength at Solidification Temp ({T_at_solidify_point-273.15:.1f}°C): {yield_strength_at_T/1e3:.2f} kPa")

    PASS = max_stress < yield_strength_at_T and z_solidify < 0.05 # Example criteria: Solidify fast and be strong enough
    print("\n--- Feasibility Check (Gravity Pull Only) ---")
    print(f"Solidification Distance (z < Tg): {z_solidify*1000:.1f} mm")
    print(f"Gravity Stress < Yield Strength at Solidification Point?: {max_stress < yield_strength_at_T} ({max_stress/1e3:.1f} < {yield_strength_at_T/1e3:.1f} kPa)")
    print(f"Likely Feasible (Needs Fast Solidification & Strength > Stress): {PASS}")
    print("--- NOTE: This ignores die swell, which gravity alone cannot correct! ---")

    return {'z_solidify_m': z_solidify, 'max_stress_Pa': max_stress, 'yield_strength_Pa': yield_strength_at_T, 'PASS': PASS}


# --- Main Script Logic ---
if __name__ == "__main__":
    print("----------------------------------------------------")
    print(" Jute/PLA Filament Extrusion Modeler (Gravity Pull) ")
    print("      --- THEORETICAL & HIGHLY APPROXIMATE ---      ")
    print("----------------------------------------------------")

    # --- Get User Inputs ---
    print("\n--- Please Provide Material Information ---")
    pla_grade = input("Enter PLA Grade Name (e.g., NatureWorks 4043D): ")
    jute_type_in = input("Enter Jute Type/Treatment (e.g., Raw Corchorus O.): ")
    jute_load_in = float(input("Enter Jute Loading by Weight (%): ")) / 100.0

    print("\n--- Please Provide Extruder Geometry (meters) ---")
    screw_D_in = float(input("Enter Screw Diameter (e.g., 0.016): "))
    channel_H_in = float(input("Enter Approx. Screw Channel Depth (e.g., 0.002): "))
    die_D_in = float(input("Enter Target Filament Diameter (e.g., 0.00175): "))
    die_L_in = float(input("Enter Die Land Length (e.g., 0.005): "))

    print("\n--- Please Provide Operating Conditions ---")
    rpm_in = float(input("Enter Screw Speed (RPM): "))
    T_air_in_C = float(input("Enter Ambient Air Temperature (°C): "))
    print("Enter Heater Zone Set Temperatures (°C):")
    T_feed_in = float(input("  Feed Zone Temp (°C): "))
    T_back_in = float(input("  Back Zone Temp (°C): "))
    T_middle_in = float(input(" Middle Zone Temp (°C): "))
    T_front_in = float(input("  Front/Die Zone Temp (°C): "))
    T_zones_in_C = [T_feed_in, T_back_in, T_middle_in, T_front_in]

    print("\n--- Please Provide Cooling & Setup Estimates ---")
    h_conv_in = float(input("Estimate Convective Cooling Coefficient 'h' (W/m^2 K) (e.g., 15=calm air, 50-100=good fan): "))
    L_max_hang_in = float(input("Enter Max Filament Hanging Length (m) (e.g., 1.0): "))

    # --- Run Calculations ---
    try:
        # 1. Get Material Properties
        pla_data = get_pla_properties(pla_grade)
        jute_data = get_jute_properties(jute_type_in)
        composite_data = calculate_composite_properties(pla_data, jute_data, jute_load_in)

        # Check Degradation Temp
        if T_front_in > jute_data['T_degradation_jute']:
            print(f"\n*** CRITICAL WARNING: Front zone temp ({T_front_in}°C) exceeds estimated Jute degradation temp ({jute_data['T_degradation_jute']}°C)! ***")

        # 2. Estimate Extrusion Process
        extrusion_results = estimate_extrusion(rpm_in, screw_D_in, channel_H_in, die_D_in, die_L_in, T_zones_in_C, composite_data)

        if extrusion_results:
            # 3. Model Cooling & Stress
            cooling_results = model_cooling_and_stress(
                extrusion_results['Q'],
                extrusion_results['T_melt_exit_K'],
                die_D_in,
                T_air_in_C + 273.15,
                h_conv_in,
                L_max_hang_in,
                composite_data
            )

            # 4. Final Summary (using results)
            print("\n\n--- FINAL PREDICTION SUMMARY ---")
            print(f"Conditions: {jute_load_in*100:.1f}% Jute, Zones=[{T_feed_in}, {T_back_in}, {T_middle_in}, {T_front_in}]°C, RPM={rpm_in}")
            print(f"Predicted Melt Exit Temp: {extrusion_results['T_melt_exit_K'] - 273.15:.1f} °C")
            print(f"Predicted Solidification Distance (below Tg): {cooling_results['z_solidify_m']*1000:.1f} mm")
            print(f"Max Gravity Stress / Est. Yield Strength: {cooling_results['max_stress_Pa']/1e3:.1f} kPa / {cooling_results['yield_strength_Pa']/1e3:.1f} kPa")
            print(f"Feasibility Prediction (Gravity Pull Only): {'PASS (Potentially)' if cooling_results['PASS'] else 'FAIL (Likely Stretching/Inconsistency)'}")
            print("Reminder: Model uses many assumptions (viscosity, flow rate, cooling) and ignores die swell.")
            print("Experimental validation is ESSENTIAL.")

    except Exception as e:
        print("\n--- AN ERROR OCCURRED DURING CALCULATION ---")
        print(e)
        import traceback
        traceback.print_exc()


    print("\n----------------------------------------------------")
    print(" End of Simulation ")
    print("----------------------------------------------------")
