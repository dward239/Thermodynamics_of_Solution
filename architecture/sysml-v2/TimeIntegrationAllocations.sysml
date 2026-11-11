// Functional allocation view: conceptual mapping of behavior to software parts.
// This does not auto-generate C++ or assert runtime SysML execution.
package RadioThermochemicalTimeAllocations {
    private import RadioThermochemicalConnections::*;
    private import RadioThermochemicalTimeBehavior::*;

    allocation def SoftwareFunctionAllocation;

    part system : AdaptiveVesselSoftwareIBD;
    action adaptiveTrial : EvaluateAdaptiveTrial;
    action acceptedTrial : AcceptedBranch;
    action rejectedTrial : RejectedBranch;

    allocation trialExecution : SoftwareFunctionAllocation
        allocate adaptiveTrial to system.controller;
    allocation acceptedHandling : SoftwareFunctionAllocation
        allocate acceptedTrial to system.controller;
    allocation rejectedHandling : SoftwareFunctionAllocation
        allocate rejectedTrial to system.controller;
}
