# Arbitrary-phase neutral-solute partitioning foundation

The phase partitioning kernel solves an **independent, analytically tractable**
mass-conserving case: each neutral dilute solute has a user-specified
constant ratio K(ip)=concentration in phase p / concentration in reference.
All phase volumes are fixed, all phases are present, and partitioning
does not change other solute activity coefficients.

For total moles n_i, reference concentration is

  c_i,ref = n_i / sum_p (K_i,p V_p)

and c_i,p = K_i,p c_i,ref. This supports an arbitrary finite number of
liquid and gas phases and an arbitrary finite number of neutral solutes.
All concentrations are mol/m3 of the specific phase; no molarity-versus-
molality conversion is implied. Gas phase concentrations are not pressures.

**Limits:** This is not rigorous nonideal liquid-liquid equilibrium,
does not solve chemical speciation, does not handle electrolytes,
phase disappearance, reactive species, gas EOS, saturation or temperature-
dependent partition constants. Coefficients need independent experimental
validation. This kernel has not yet been integrated with the radioactive
closed-vessel model; invoking it requires an explicit input mapping.
