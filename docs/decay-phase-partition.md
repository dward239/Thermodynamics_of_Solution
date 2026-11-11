# Decay-to-neutral-partition end-to-end demonstration

The standalone bridge first advances a generic first-order decay network,
converts remaining atom counts into element-total moles, and then applies
the analytically mass-conserving linear partitioning kernel to distribute
**one neutral pseudo-solute per chemical element** across fixed volumes.

This is an architectural coupling test, NOT a mechanistic radiochemical
simulation. Actual radionuclides can have multiple oxidation states,
complexes, and solids with very different phase behavior; those cannot
be represented by this one-solute-per-element assumption. No aqueous
electroneutrality, kinetic source, or real equilibrium constant is solved.

The test uses deliberately synthetic 10-second half-life and user-defined
aqueous/organic partition ratios. It validates bookkeeping and rejection
when a daughter element has no partition model.
