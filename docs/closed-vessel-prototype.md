# Restricted zero-dimensional vessel demonstration

The vessel simulation advances branched radioactive decay, prescribed phase
radiation deposition, gross radiolysis yields, a 0D heat balance, and ideal
gas headspace pressure in a single step.

Gross radiolysis products are **not** automatically added to gas inventory.
Only an explicitly supplied *net* gas transfer may alter headspace moles.
Recombination, phase partitioning, kinetics, chemistry and gas release are
not resolved, and the pressure calculation uses fixed volume and the ideal
gas law. The nuclear data in tests are strictly synthetic.

This is not a validated multiphase equilibrium or process safety solver.
