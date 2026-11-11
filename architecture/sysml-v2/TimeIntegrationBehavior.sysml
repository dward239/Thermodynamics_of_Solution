// SysML v2 action and state views corresponding to rendered ACT and STM.
// Temporal control here is conceptual: the implementation is authoritative
// for guards, retry policies, step ordering and fail-fast behavior.
package RadioThermochemicalTimeBehavior {
    private import ScalarValues::*;

    item def TrialAcceptedEvent;
    item def TrialRejectedEvent;
    item def IntegrationCompletedEvent;
    item def IntegrationFailedEvent;
    item def StartIntegrationEvent;

    action def ProposeTrialStep;
    action def EvaluateFullStep;
    action def EvaluateFirstHalf;
    action def EvaluateSecondHalf;
    action def CompareLocalError;
    action def CheckConservation;
    action def CommitAcceptedState;
    action def ReduceStepSize;

    action def EvaluateAdaptiveTrial {
        action propose : ProposeTrialStep;
        action full : EvaluateFullStep;
        action half1 : EvaluateFirstHalf;
        action half2 : EvaluateSecondHalf;
        action estimate : CompareLocalError;
        first propose then full;
        first full then half1;
        first half1 then half2;
        first half2 then estimate;
        // Acceptance/rejection is a decision of the state machine below.
    }
    action def AcceptedBranch {
        action audit : CheckConservation;
        action commit : CommitAcceptedState;
        first audit then commit;
    }
    action def RejectedBranch {
        action retry : ReduceStepSize;
    }

    state def AdaptiveTimeControllerStates {
        state ready;
        state evaluating;
        state auditing;
        state adjusting;
        state committed;
        state completed;
        state failed;

        transition start first ready
            accept StartIntegrationEvent then evaluating;
        transition acceptable first evaluating
            accept TrialAcceptedEvent then auditing;
        transition rejected first evaluating
            accept TrialRejectedEvent then adjusting;
        transition accepted first auditing
            accept TrialAcceptedEvent then committed;
        transition auditFailure first auditing
            accept IntegrationFailedEvent then failed;
        transition retry first adjusting
            accept StartIntegrationEvent then evaluating;
        transition retryFailure first adjusting
            accept IntegrationFailedEvent then failed;
        transition next first committed
            accept StartIntegrationEvent then evaluating;
        transition finished first committed
            accept IntegrationCompletedEvent then completed;
    }
    state simulationLifeCycle : AdaptiveTimeControllerStates;
}
