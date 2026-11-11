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
        port candidateTemperatureOut : ThermalPort;
        port chemistryIn : ~EquilibriumProjectionPort;
    }
    part def NuclideElementMapper {
        port nuclideIn : ~NuclidePort;
        port elementsOut : ElementInventoryPort;
    }
    part def ReducedChemicalEquilibrium {
        port elementsIn : ~ElementInventoryPort;
        port candidateTemperatureIn : ~ThermalPort;
        port equilibriumOut : EquilibriumProjectionPort;
    }
    part def AcceptedStepChemicalIBD {
        part controller : AcceptedStepController;
        part elementMapping : NuclideElementMapper;
        part equilibrium : ReducedChemicalEquilibrium;
        connect controller.nuclideOut to elementMapping.nuclideIn;
        connect elementMapping.elementsOut to equilibrium.elementsIn;
        connect controller.candidateTemperatureOut to equilibrium.candidateTemperatureIn;
        connect equilibrium.equilibriumOut to controller.chemistryIn;
    }

    action def MapNuclideAtomsToElementMol;
    action def AddBackgroundComponentInventories;
    action def SuppressZeroInventorySpeciesAndSolids;
    action def EvaluateConfiguredStandardStates;
    action def SolveNeutralMultiphaseEquilibrium;
    action def AuditFluidAndSolidElementConservation;
    action def CommitCandidateTimeState;

    action def AcceptedStepChemistryProjection {
        action map : MapNuclideAtomsToElementMol;
        action background : AddBackgroundComponentInventories;
        action removeAbsent : SuppressZeroInventorySpeciesAndSolids;
        action standardStates : EvaluateConfiguredStandardStates;
        action solve : SolveNeutralMultiphaseEquilibrium;
        action audit : AuditFluidAndSolidElementConservation;
        action commit : CommitCandidateTimeState;
        first map then background;
        first background then removeAbsent;
        first removeAbsent then standardStates;
        first standardStates then solve;
        first solve then audit;
        first audit then commit;
    }

    item def AcceptedCandidateEvent;
    item def MappingCompletedEvent;
    item def StandardStatesReadyEvent;
    item def EquilibriumConvergedEvent;
    item def ConservationPassedEvent;
    item def ProjectionFailedEvent;

    state def AcceptedChemistryProjectionStates {
        state awaitingCandidate;
        state mappingInventory;
        state evaluatingStandardStates;
        state solvingEquilibrium;
        state auditingConservation;
        state readyToCommit;
        state failed;

        transition start first awaitingCandidate
            accept AcceptedCandidateEvent then mappingInventory;
        transition mapped first mappingInventory
            accept MappingCompletedEvent then evaluatingStandardStates;
        transition propertiesReady first evaluatingStandardStates
            accept StandardStatesReadyEvent then solvingEquilibrium;
        transition solved first solvingEquilibrium
            accept EquilibriumConvergedEvent then auditingConservation;
        transition conserved first auditingConservation
            accept ConservationPassedEvent then readyToCommit;
        transition failedMap first mappingInventory
            accept ProjectionFailedEvent then failed;
        transition failedProperties first evaluatingStandardStates
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
