// Software architecture baseline in SysML v2 textual notation.
// Part definitions/usages represent computational components and their states.
// An external SysML v2 parser is not integrated into CI yet.
package RadioThermochemicalArchitecture {
    private import ScalarValues::*;
    part def ChemicalSystemModel {
        attribute phaseCount : Integer;
        attribute speciesCount : Integer;
        attribute reactionCount : Integer;
    }
    part def NuclideInventoryModel {
        attribute nuclideCount : Integer;
        attribute atomInventory : Real;
    }
    part def NuclearDecayModel {
        attribute elapsedSeconds : Real;
        attribute emittedEnergyJ : Real;
    }
    part def RadiationDepositionModel {
        attribute absorbedEnergyJ : Real;
    }
    part def RadiolysisSourceModel {
        attribute grossGeneratedMoles : Real;
    }
    part def ThermalBalanceModel {
        attribute temperatureK : Real;
        attribute sensibleEnergyJ : Real;
    }
    part def GeneralIdealFluidEquilibriumModel {
        attribute componentCount : Integer;
        attribute fluidSpeciesCount : Integer;
        attribute normalizedMassResidual : Real;
    }
    part def ChemicalEquilibriumModel {
        attribute normalizedResidual : Real;
    }
    part def NuclearChemicalInventoryBridge {
        attribute elementInventoryMol : Real;
    }
    part def DilutePhasePartitionModel {
        attribute partitionedSoluteMol : Real;
    }
    part def HeadspacePressureModel {
        attribute pressurePa : Real;
    }
    part def NeutralDecayPartitionBridge {
        attribute daughterElementMoles : Real;
    }
    part def PureSolidPhaseSelector {
        attribute candidateSolidCount : Integer;
        attribute activeSolidCount : Integer;
    }
    part def SinglePhaseChargeBalance {
        attribute chargedPhaseCount : Integer;
        attribute freeCounterionTotalMol : Real;
        attribute chargeResidualMol : Real;
    }
    part def DaviesActivityCoupling {
        attribute ionicStrengthMolL : Real;
        attribute outerConsistencyResidual : Real;
    }
    part def RadiolyticHydrogenKinetics {
        attribute dissolvedHydrogenMol : Real;
        attribute headspaceHydrogenMol : Real;
        attribute consumedHydrogenMol : Real;
    }
    part def AdaptiveTimeIntegratorModel {
        attribute acceptedStepCount : Integer;
        attribute rejectedStepCount : Integer;
        attribute currentTimeSeconds : Real;
        attribute localErrorEstimate : Real;
    }
    part def AcceptedChemicalProjectionModel {
        attribute activeComponentCount : Integer;
        attribute activeSolidCount : Integer;
        attribute maximumScaledComponentClosure : Real;
    }
    part def SimulationController {
        part chemicalSystem : ChemicalSystemModel;
        part inventory : NuclideInventoryModel;
        part decay : NuclearDecayModel;
        part deposition : RadiationDepositionModel;
        part radiolysis : RadiolysisSourceModel;
        part thermal : ThermalBalanceModel;
        part equilibrium : ChemicalEquilibriumModel;
        part idealFluid : GeneralIdealFluidEquilibriumModel;
        part solidSelection : PureSolidPhaseSelector;
        part chargeBalance : SinglePhaseChargeBalance;
        part daviesCoupling : DaviesActivityCoupling;
        part hydrogenKinetics : RadiolyticHydrogenKinetics;
        part adaptiveIntegrator : AdaptiveTimeIntegratorModel;
        part acceptedChemistry : AcceptedChemicalProjectionModel;
        part inventoryBridge : NuclearChemicalInventoryBridge;
        part partition : DilutePhasePartitionModel;
        part headspace : HeadspacePressureModel;
        part decayPartitionBridge : NeutralDecayPartitionBridge;
    }
    part simulator : SimulationController;
}
