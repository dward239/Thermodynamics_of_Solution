// Logical software IBD in SysML v2 textual notation.
// This shows FUNCTION-CALL / DATA contracts, not physical fluid piping.
// No SysML v2 parser is currently part of CI; graphical versions live in views/.
package RadioThermochemicalConnections {
    private import ScalarValues::*;

    item def StepRequestData {
        attribute startTimeSeconds : Real;
        attribute requestedDurationSeconds : Real;
    }
    item def StepResultData {
        attribute endTimeSeconds : Real;
        attribute temperatureK : Real;
        attribute pressurePa : Real;
    }
    item def NuclideStateData {
        attribute atomCount : Real;
    }
    item def RadiationEnergyData {
        attribute depositedEnergyJ : Real;
    }
    item def HydrogenInventoryData {
        attribute dissolvedMol : Real;
        attribute headspaceMol : Real;
    }
    item def ThermalStateData {
        attribute temperatureK : Real;
    }

    port def StepRequestPort { out item request : StepRequestData; }
    port def StepResultPort { out item result : StepResultData; }
    port def NuclidePort { out item inventory : NuclideStateData; }
    port def RadiationPort { out item deposition : RadiationEnergyData; }
    port def HydrogenPort { out item inventory : HydrogenInventoryData; }
    port def ThermalPort { out item temperature : ThermalStateData; }

    part def AdaptiveController {
        port requestOut : StepRequestPort;
        port resultIn : ~StepResultPort;
    }
    part def VesselStepKernel {
        port requestIn : ~StepRequestPort;
        port resultOut : StepResultPort;
        port nuclearStateOut : NuclidePort;
        port nuclearStateIn : ~NuclidePort;
        port radiationOut : RadiationPort;
        port radiationIn : ~RadiationPort;
        port hydrogenOut : HydrogenPort;
        port hydrogenIn : ~HydrogenPort;
        port thermalOut : ThermalPort;
        port thermalIn : ~ThermalPort;
    }
    part def DecayKernel {
        port inventoryIn : ~NuclidePort;
        port inventoryOut : NuclidePort;
        port radiationOut : RadiationPort;
    }
    part def RadiolysisHydrogenKernel {
        port depositionIn : ~RadiationPort;
        port hydrogenOut : HydrogenPort;
    }
    part def ThermalBalanceKernel {
        port energyIn : ~RadiationPort;
        port temperatureOut : ThermalPort;
    }
    // The connections express software data dependencies of implemented
    // modules; an actual C++ call may exchange a richer struct than
    // the simplified scalar attributes in these conceptual data items.
    part def AdaptiveVesselSoftwareIBD {
        part controller : AdaptiveController;
        part vessel : VesselStepKernel;
        part decay : DecayKernel;
        part hydrogen : RadiolysisHydrogenKernel;
        part thermal : ThermalBalanceKernel;

        connect controller.requestOut to vessel.requestIn;
        connect vessel.resultOut to controller.resultIn;

        connect vessel.nuclearStateOut to decay.inventoryIn;
        connect decay.inventoryOut to vessel.nuclearStateIn;
        connect decay.radiationOut to vessel.radiationIn;

        connect vessel.radiationOut to hydrogen.depositionIn;
        connect hydrogen.hydrogenOut to vessel.hydrogenIn;
        connect vessel.radiationOut to thermal.energyIn;
        connect thermal.temperatureOut to vessel.thermalIn;
    }

    part connectedAdaptiveVessel : AdaptiveVesselSoftwareIBD;
}
