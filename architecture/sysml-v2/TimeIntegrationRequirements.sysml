// Candidate requirements and test allocation in SysML v2 notation.
// These are architecture contracts, not validated numerical proof in SysML.
package RadioThermochemicalTimeRequirements {
    private import ScalarValues::*;

    part def AdaptiveVesselSubject {
        attribute acceptedStepCount : Integer;
        attribute rejectedStepCount : Integer;
        attribute endTimeSeconds : Real;
        attribute nuclearClosureAtoms : Real;
        attribute hydrogenClosureMol : Real;
    }
    requirement def TimeProgressRequirement {
        subject model : AdaptiveVesselSubject;
        require constraint { model.acceptedStepCount >= 1 }
    }
    requirement def NonnegativeStepCountRequirement {
        subject model : AdaptiveVesselSubject;
        require constraint { model.rejectedStepCount >= 0 }
    }
    requirement def FiniteEndTimeRequirement {
        subject model : AdaptiveVesselSubject;
        require constraint { model.endTimeSeconds >= 0.0 }
    }
    requirement def NuclearConservationRequirement {
        subject model : AdaptiveVesselSubject;
        require constraint { model.nuclearClosureAtoms == 0.0 }
    }
    requirement def HydrogenConservationRequirement {
        subject model : AdaptiveVesselSubject;
        require constraint { model.hydrogenClosureMol == 0.0 }
    }

    requirement timeProgress : TimeProgressRequirement;
    requirement nonnegativeStepCounts : NonnegativeStepCountRequirement;
    requirement finiteEndTime : FiniteEndTimeRequirement;
    requirement nuclearConservation : NuclearConservationRequirement;
    requirement hydrogenConservation : HydrogenConservationRequirement;
    // Finite-precision acceptance criteria and actual test evidence are
    // specified in architecture/sysml-v2/requirements-traceability.md.
}
