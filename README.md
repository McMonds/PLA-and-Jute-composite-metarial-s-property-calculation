# PLA-and-Jute-composite-metarial-s-property-calculation
Python script for calculating physical and mathematical properties of PLA and Jute composite materials and predict
a fine tune temperature for Fliabot EX6 extruder.


It asks for information related to material properties and extruder tools details in precision 
then it calculates and gives a temperature combination.


its a prototype.

Here is the Key Physics & Math Concepts Employed

    Material Science: Rule of Mixtures (for density, Cp), Nielsen/Maxwell models (for k_comp), 
	Polymer Rheology (Viscosity models: Newtonian, Power-Law, Arrhenius, effects of fillers), Solid Mechanics (Yield Strength vs Temperature).

    Thermodynamics: Heat Capacity, Latent Heat (if modeling melting precisely), Glass Transition, Melting Point, Thermal Degradation.

    Heat Transfer:

        Conduction (Fourier's Law: q = -k∇T) - Heating from barrel, heat flow within melt/filament.

        Convection (Newton's Law of Cooling: q = hA(T_s - T_f)) - Cooling by air/fan.

        Radiation (Stefan-Boltzmann Law: q = εσA(T_s⁴ - T_surr⁴)) - Cooling by radiation.

        Shear Heating Generation (Q_shear ≈ η * γ̇² * Volume) - Frictional heat within the melt.

        Transient Heat Equation (ρCp ∂T/∂t = ∇·(k∇T) + Q_gen) - Governing equation for temperature change over time/distance.

    Fluid Mechanics (Rheology):

        Viscosity (η): Resistance to flow.

        Shear Rate (γ̇): Rate of deformation (related to RPM, geometry, Q).

        Shear Stress (τ = η * γ̇).

        Flow Rate (Q): Volume of material moved per unit time (related to screw characteristics, RPM, pressure).

        Pressure Drop (ΔP): Pressure difference needed to push melt through the die (Hagen-Poiseuille for Newtonian, more complex for Power-Law).

        Die Swell: Expansion of polymer upon exiting the die.

    Statics/Solid Mechanics:

        Stress (σ = Force/Area).

        Gravitational Force (F = mg = ρ * Volume * g).

        Condition for No Yielding: σ_gravity < yield_strength(T).

    Mathematics: Algebra, Calculus (Derivatives for rates, Gradients for heat flow, Integrals if averaging), 
	Differential Equations (Transient Heat Equation), Numerical Methods (Finite Difference/Element for solving PDE), Unit Conversions.
