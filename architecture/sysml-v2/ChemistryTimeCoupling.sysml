// Accepted-step chemical projection: *software* architecture, not a DAE.
// This is an experimental SysML v2 textual view; not parsed in CI.
package RadioThermochemicalChemistryCoupling {
    private import ScalarValues::*;
    private import RadioThermochemicalConnections::*;

    item def ElementInventoryData {
        attribute componentCount : Integer;
        attribute totalElementMol : Real;
    }
    item def EquilibriumProjectionData {
        attribute speciesCount : Integer;
        attribute solidPhaseCount : Integer;
        attribute maximumScaledMassResidual : Real;
    }
    port def ElementInventoryPort {
        out item inventory : ElementInventoryData;
    }
    port def EquilibriumProjectionPort {
        out item state : EquilibriumProjectionData;
    }

    part def AcceptedStepController {
        port nuclideOut : NuclidePort;
        port chemistryIn : ~EquilibriumProjectionPort;
    }
    part def NuclideElementMapper {
        port nuclideIn : ~NuclidePort;
        port elementsOut : ElementInventoryPort;
    }
    part def ReducedChemicalEquilibrium {
        port elementsIn : ~ElementInventoryPort;
        port equilibriumOut : EquilibriumProjectionPort;
    }
    part def AcceptedStepChemicalIBD {
        part controller : AcceptedStepController;
        part elementMapping : NuclideElementMapper;
        part equilibrium : ReducedChemicalEquilibrium;
        connect controller.nuclideOut to elementMapping.nuclideIn;
        connect elementMapping.elementsOut to equilibrium.elementsIn;
        connect equilibrium.equilibriumOut to controller.chemistryIn;
    }

    action def MapNuclideAtomsToElementMol;
    action def AddBackgroundComponentInventories;
    action def SuppressZeroInventorySpeciesAndSolids;
    action def SolveNeutralMultiphaseEquilibrium;
    action def AuditFluidAndSolidElementConservation;
    action def CommitCandidateTimeState;

    action def AcceptedStepChemistryProjection {
        action map : MapNuclideAtomsToElementMol;
        action background : AddBackgroundComponentInventories;
        action removeAbsent : SuppressZeroInventorySpeciesAndSolids;
        action solve : SolveNeutralMultiphaseEquilibrium;
        action audit : AuditFluidAndSolidElementConservation;
        action commit : CommitCandidateTimeState;
        first map then background;
        first background then removeAbsent;
        first removeAbsent then solve;
        first solve then audit;
        first audit then commit;
    }

    item def AcceptedCandidateEvent;
    item def MappingCompletedEvent;
    item def EquilibriumConvergedEvent;
    item def ConservationPassedEvent;
    item def ProjectionFailedEvent;

    state def AcceptedChemistryProjectionStates {
        state awaitingCandidate;
        state mappingInventory;
        state solvingEquilibrium;
        state auditingConservation;
        state readyToCommit;
        state failed;

        transition start first awaitingCandidate
            accept AcceptedCandidateEvent then mappingInventory;
        transition mapped first mappingInventory
            accept MappingCompletedEvent then solvingEquilibrium;
        transition solved first solvingEquilibrium
            accept EquilibriumConvergedEvent then auditingConservation;
        transition conserved first auditingConservation
            accept ConservationPassedEvent then readyToCommit;
        transition failedMap first mappingInventory
            accept ProjectionFailedEvent then failed;
        transition failedSolve first solvingEquilibrium
            accept ProjectionFailedEvent then failed;
        transition failedAudit first auditingConservation
            accept ProjectionFailedEvent then failed;
    }

    part chemicalConnection : AcceptedStepChemicalIBD;
    action chemistryOnAccept : AcceptedStepChemistryProjection;
    state chemistryLifecycle : AcceptedChemistryProjectionStates;
}
