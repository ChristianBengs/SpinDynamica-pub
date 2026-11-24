(* ::Package:: *)

(* SpinDynamica Optimal Control — Time-Optimal Pulse Search, TOP Curves, Joint Constraints & PMP Diagnostics *)

BeginPackage["OptimalControl`TimeOptimal`", {
    "Global`",
    "SpinDynamica`",
    "QM`Hilbert`",
    "QM`SpinOperators`",
    "QM`Liouville`",
    "OptimalControl`Kernels`",
    "OptimalControl`Objectives`",
    "OptimalControl`GRAPE`",
    "OptimalControl`Ensemble`"
}];

TimeOptimalControlOptimize::usage = "TimeOptimalControlOptimize[problemOrEnsemble, opts] searches for the shortest numerically identified pulse duration satisfying TargetFidelity under experimental constraints, supporting \"Progressive\" duration reduction, \"Joint\" fidelity-constrained Augmented Lagrangian optimization, \"JointPenalized\" DurationPenalty sweeps, and \"Hybrid\" refinement.";
TimeOptimalPulseCurve::usage = "TimeOptimalPulseCurve[problemOrEnsemble, durations, opts] optimizes controls across a grid of pulse durations using warm-start continuation and returns a structured numerical time-optimal pulse (TOP) curve.";
VerifyOptimalControlStationarity::usage = "VerifyOptimalControlStationarity[result, opts] evaluates continuous-time Pontryagin Maximum Principle switching functions s_k(t) on a fine sub-slice grid and reports integrated and discrete stationarity residuals, quadrature consistency, and an explicit scope disclaimer.";

InitialDuration::usage = "InitialDuration is an option for TimeOptimalControlOptimize specifying the initial pulse duration in seconds (default Automatic, taken from the input problem timing).";
MinimumDuration::usage = "MinimumDuration is an option for TimeOptimalControlOptimize specifying the lower duration bound in seconds for the search (default Automatic, set to 0.05 * InitialDuration).";
DurationTolerance::usage = "DurationTolerance is an option for TimeOptimalControlOptimize specifying the bracket stopping tolerance in seconds (default Automatic).";
DurationReductionFactor::usage = "DurationReductionFactor is an option for TimeOptimalControlOptimize specifying the multiplicative reduction factor r in (0, 1) during progressive duration reduction (default 0.8).";
FidelityCriterion::usage = "FidelityCriterion is an option for TimeOptimalControlOptimize and TimeOptimalPulseCurve selecting the feasibility metric: Automatic, \"Fidelity\", \"Mean\", \"SoftMin\", or \"Minimum\".";
DiscretizationStrategy::usage = "DiscretizationStrategy is an option for TimeOptimalControlOptimize and TimeOptimalPulseCurve selecting \"FixedSlices\", \"FixedSliceDuration\", or {\"FixedSliceDuration\", dtRef} when rescaling pulse duration.";
WaveformRescaling::usage = "WaveformRescaling is an option for TimeOptimalControlOptimize and TimeOptimalPulseCurve controlling initial-guess continuation between durations: \"NormalizedTime\", \"AreaPreserving\", or \"None\".";
AttemptsPerDuration::usage = "AttemptsPerDuration is an option for TimeOptimalControlOptimize and TimeOptimalPulseCurve specifying {coarseAttempts, boundaryAttempts} or an integer number of optimization restarts per duration.";
ValidationEnsemble::usage = "ValidationEnsemble is an option for TimeOptimalControlOptimize specifying a denser out-of-sample OptimalControlEnsemble on which the final pulse is independently evaluated.";
StationaritySubdivisions::usage = "StationaritySubdivisions is an option for VerifyOptimalControlStationarity specifying the number of sub-slice time points per control slice (default 8).";

TimeOptimalControlOptimize::badprob = "Expected an OptimalControlProblem, OptimalControlEnsemble, CooperativeObjective, or problem list; got `1`.";
TimeOptimalControlOptimize::baddur = "Invalid duration settings `1`.";
TimeOptimalControlOptimize::badopt = "Invalid option specification `1`: `2`.";
TimeOptimalControlOptimize::infeasibleinit = "Initial duration `1` s did not numerically reach TargetFidelity `2` (best achieved criterion fidelity `3`). Failure to converge does not prove physical infeasibility.";
TimeOptimalPulseCurve::baddurs = "Durations must be a non-empty list of positive real numbers; got `1`.";
VerifyOptimalControlStationarity::badinput = "VerifyOptimalControlStationarity expects a valid OptimalControlResult; got `1`.";

Begin["`Private`"];

$OptimalityDisclaimer = "Numerically identified shortest feasible duration within the chosen control discretization, optimizer tolerances, and restart budget; neither failure of GRAPE to converge at shorter durations nor small local residuals constitute a rigorous proof of physical or global time optimality.";
$PMPDisclaimer = "Pontryagin Maximum Principle switching-function residuals measure local first-order constrained stationarity at fixed total duration; small residuals are necessary rather than sufficient for global or time optimality, and free-final-time optimality further requires Hamiltonian transversality.";

Options[TimeOptimalControlOptimize] = {
    Method -> "Progressive",
    TargetFidelity -> 0.999,
    InitialDuration -> Automatic,
    MinimumDuration -> Automatic,
    DurationTolerance -> Automatic,
    DurationReductionFactor -> 0.8,
    FidelityCriterion -> Automatic,
    DiscretizationStrategy -> "FixedSlices",
    WaveformRescaling -> "NormalizedTime",
    AttemptsPerDuration -> {1, 3},
    ValidationEnsemble -> None,
    MaxIterations -> 300,
    RandomSeed -> Automatic,
    Verbose -> False
};

Options[TimeOptimalPulseCurve] = {
    TargetFidelity -> 0.999,
    FidelityCriterion -> Automatic,
    DiscretizationStrategy -> "FixedSlices",
    WaveformRescaling -> "NormalizedTime",
    AttemptsPerDuration -> {1, 2},
    ValidationEnsemble -> None,
    MaxIterations -> 300,
    RandomSeed -> Automatic,
    Verbose -> False
};

Options[VerifyOptimalControlStationarity] = {
    StationaritySubdivisions -> 8,
    Tolerance -> 10^-3
};

$TOOnlyOptionNames = {
    InitialDuration, MinimumDuration, DurationTolerance, DurationReductionFactor,
    FidelityCriterion, DiscretizationStrategy, WaveformRescaling,
    AttemptsPerDuration, ValidationEnsemble
};

(* Normalize accepted problem/ensemble wrappers into canonical OptimalControlProblem or OptimalControlEnsemble *)
normalizeTOInput[input_] := Which[
    MatchQ[input, _OptimalControlProblem] && Length[input] >= 4,
        {"Single", input, List @@ input[[5 ;;]]},
    MatchQ[input, _OptimalControlEnsemble],
        With[{own = Which[
            Length[input] >= 2 && ListQ[input[[1]]] && !MatrixQ[input[[1]]], List @@ input[[2 ;;]],
            Length[input] >= 5, List @@ input[[5 ;;]],
            True, {}
        ]},
            {"Ensemble", input, own}
        ],
    MatchQ[input, {(_OptimalControlProblem | {_?NumericQ, _OptimalControlProblem}) ..}],
        {"Ensemble", OptimalControlEnsemble[input], {}},
    MatchQ[input, CooperativeObjective[{(_OptimalControlProblem | {_?NumericQ, _OptimalControlProblem}) ..}, ___Rule]],
        With[{cOpts = List @@ Rest[input]},
            {"Ensemble", OptimalControlEnsemble[First[input], Sequence @@ cOpts], cOpts}
        ],
    True,
        $Failed
];

hasSliceListDriftQ["Single", prob_OptimalControlProblem] :=
    ListQ[prob[[1]]] && !MatrixQ[prob[[1]], NumericQ];

hasSliceListDriftQ["Ensemble", OptimalControlEnsemble[items_List, ___]] :=
    AnyTrue[items, Function[item,
        With[{p = If[ListQ[item], item[[2]], item]},
            ListQ[p[[1]]] && !MatrixQ[p[[1]], NumericQ]
        ]
    ]];

hasSliceListDriftQ["Ensemble", OptimalControlEnsemble[drift_, ___]] :=
    ListQ[drift] && !MatrixQ[drift, NumericQ];

extractBaselineTiming["Single", prob_OptimalControlProblem, defaultNt_Integer] :=
    ($OCSolverRegistry["ResolveTimingVector"])[prob[[4]], defaultNt];

extractBaselineTiming["Ensemble", ens_OptimalControlEnsemble, defaultNt_Integer] :=
    ($OCSolverRegistry["ExtractEnsembleTiming"])[ens, defaultNt];

retargetObject["Single", prob_OptimalControlProblem, dtVec_List] :=
    OptimalControlProblem[prob[[1]], prob[[2]], prob[[3]], dtVec, Sequence @@ (List @@ prob[[5 ;;]])];

retargetObject["Ensemble", ens_OptimalControlEnsemble, dtVec_List] :=
    ($OCSolverRegistry["RetargetEnsembleTiming"])[ens, dtVec];

resolveEffectiveCriterion["Single", critOpt_, ownAndUserOpts_List] :=
    If[MemberQ[{"Fidelity", "Mean", "SoftMin", "Minimum"}, critOpt],
        If[critOpt === "Minimum" || critOpt === "Mean" || critOpt === "SoftMin", "Fidelity", critOpt],
        If[critOpt === Automatic, "Fidelity", $Failed]
    ];

resolveEffectiveCriterion["Ensemble", critOpt_, ownAndUserOpts_List] := Module[{agg},
    agg = EnsembleAggregation /. ownAndUserOpts /. Options[OptimalControlEnsemble];
    Which[
        critOpt === Automatic,
            If[MatchQ[agg, {"SoftMin", _?NumericQ}], "SoftMin", "Mean"],
        MemberQ[{"Mean", "SoftMin", "Minimum"}, critOpt],
            critOpt,
        critOpt === "Fidelity",
            If[MatchQ[agg, {"SoftMin", _?NumericQ}], "SoftMin", "Mean"],
        True,
            $Failed
    ]
];

validateDiscretization[strat_, dt0_List, isSliceListQ_] := Module[{dtRef},
    Which[
        strat === "FixedSlices",
            <|"Mode" -> "FixedSlices", "RefDt" -> Mean[dt0], "SliceCount0" -> Length[dt0], "SliceListQ" -> isSliceListQ|>,
        strat === "FixedSliceDuration",
            <|"Mode" -> "FixedSliceDuration", "RefDt" -> Mean[dt0], "SliceCount0" -> Length[dt0], "SliceListQ" -> isSliceListQ|>,
        MatchQ[strat, {"FixedSliceDuration", _?NumericQ}] && TrueQ[Im[N[strat[[2]]]] == 0 && Re[strat[[2]]] > 0.],
            dtRef = N[Re[strat[[2]]]];
            <|"Mode" -> "FixedSliceDuration", "RefDt" -> dtRef, "SliceCount0" -> Length[dt0], "SliceListQ" -> isSliceListQ|>,
        True,
            Message[TimeOptimalControlOptimize::badopt, DiscretizationStrategy -> strat,
                "expected \"FixedSlices\", \"FixedSliceDuration\", or {\"FixedSliceDuration\", dtRef}"];
            $Failed
    ]
];

buildDurationGrid[dur_?NumericQ, dt0_List, discSpec_Association] := Module[{ntNew, tot0},
    If[discSpec["Mode"] === "FixedSlices",
        tot0 = Max[Total[dt0], 10^-15];
        Return[Developer`ToPackedArray[N[(dur / tot0) * dt0]]]
    ];
    ntNew = Max[1, Round[N[dur] / discSpec["RefDt"]]];
    If[TrueQ[discSpec["SliceListQ"]] && ntNew =!= discSpec["SliceCount0"],
        Message[TimeOptimalControlOptimize::badopt, DiscretizationStrategy -> "FixedSliceDuration",
            "cannot alter slice count when Drift is an explicit SliceList"];
        Return[$Failed]
    ];
    Developer`ToPackedArray[ConstantArray[N[dur] / ntNew, ntNew]]
];

parseAttemptBudget[attOpt_] := Which[
    IntegerQ[attOpt] && attOpt >= 1, {attOpt, attOpt},
    MatchQ[attOpt, {_Integer?Positive, _Integer?Positive}], attOpt,
    True,
        Message[TimeOptimalControlOptimize::badopt, AttemptsPerDuration -> attOpt,
            "expected positive integer or {coarseAttempts, boundaryAttempts}"];
        $Failed
];

buildAltRescale[mode_String] := Switch[mode,
    "NormalizedTime", "AreaPreserving",
    "AreaPreserving", "NormalizedTime",
    _, "NormalizedTime"
];

(* Rescale previous waveform to new duration/discretization and enforce static control constraints *)
buildWarmStartGuess[prevRes_OptimalControlResult, newDt_List, grapeOpts_List, mode_String, mCtrl_Integer] := Module[
    {oldDt, oldDur, newDur, ntNew, uOld, uRes, uScaled,
     ctrlBnds, ampBnd, normBnd, cSpec, bIdent, pCtx, uProj, paramOpt, winOpt, ctrlBasis},
    If[mode === "None", Return[Automatic]];
    oldDt = prevRes["TimeSteps"];
    oldDur = Max[Total[oldDt], 10^-15];
    newDur = Max[Total[newDt], 10^-15];
    ntNew = Length[newDt];

    If[mode === "NormalizedTime",
        (* GRAPE buildInitialVector reuses prevRes["Parameters"] if shapes match, else interpolates Controls *)
        Return[prevRes]
    ];

    (* AreaPreserving: resample u(t), scale by T_old / T_new, project onto bounds/disks, map to basis if needed *)
    uOld = prevRes["Controls"];
    uRes = ($OCSolverRegistry["ResampleControls"])[uOld, oldDt, newDt];
    uScaled = N[uRes * (oldDur / newDur)];

    {ctrlBnds, ampBnd, normBnd, paramOpt, winOpt} =
        {ControlBounds, AmplitudeBound, AmplitudeNormBound, ControlParameterization, ControlWindow} /. grapeOpts /. Options[GRAPEOptimize];
    bIdent = BuildControlBasis["PiecewiseConstant", newDt, None];
    cSpec = ($OCSolverRegistry["ResolveConstraintSpec"])[ctrlBnds, ampBnd, normBnd, None, newDt, ntNew, mCtrl];
    If[AssociationQ[cSpec] && AssociationQ[bIdent],
        pCtx = ($OCSolverRegistry["CompileProxContext"])[bIdent, newDt, cSpec, 0., 0.];
        If[AssociationQ[pCtx], uProj = ($OCSolverRegistry["ApplyProximal"])[pCtx, uScaled, 0.], uProj = uScaled],
        uProj = uScaled
    ];
    ctrlBasis = BuildControlBasis[paramOpt, newDt, winOpt];
    If[AssociationQ[ctrlBasis] && !TrueQ[ctrlBasis["IsIdentityQ"]],
        Quiet[Check[LeastSquares[ctrlBasis["BasisMatrix"], ToPackedRealMatrix[uProj]], ToPackedRealMatrix[uProj]]],
        ToPackedRealMatrix[uProj]
    ]
];

estimateMatrixExponentials[res_OptimalControlResult] := Module[{kind, evals, comp, compEns, nObj},
    kind = Lookup[res[[1]], "Kind", "Single"];
    evals = Max[Lookup[res[[1]], "Evaluations", 0],
        Lookup[res[[1]], "ObjectiveEvaluations", 0] + Lookup[res[[1]], "GradientEvaluations", 0], 0];
    If[kind === "Ensemble",
        compEns = Lookup[res[[1]], "CompiledEnsemble", <||>];
        nObj = Max[evals, 1];
        nObj * Total[Map[Length[Lookup[#, "TimeSteps", {}]] * Lookup[#, "NumStages", 1] &, Lookup[compEns, "Members", {}]]],
        comp = Lookup[res[[1]], "CompiledProblem", <||>];
        Max[evals, 1] * Length[Lookup[comp, "TimeSteps", res["TimeSteps"]]] * Lookup[comp, "NumStages", 1]
    ]
];

summarizeResult[res_OptimalControlResult, kind_String, effCrit_String, targetFid_?NumericQ] := Module[
    {assoc = res[[1]], fid, meanF, softMinF, minF, maxViol, critF, feasQ, durVal, objEvals, gradEvals, nExp},
    fid = N[assoc["FinalFidelity"]];
    If[kind === "Ensemble",
        meanF = N[Lookup[assoc, "MeanMemberFidelity", fid]];
        softMinF = N[Lookup[assoc, "SoftMinMemberFidelity", fid]];
        minF = N[Lookup[assoc, "MinMemberFidelity", fid]],
        meanF = fid;
        softMinF = fid;
        minF = fid
    ];
    maxViol = N[Lookup[assoc, "MaxConstraintViolation", 0.]];
    critF = Switch[effCrit,
        "Minimum", minF,
        "Mean", meanF,
        "SoftMin", softMinF,
        _, fid
    ];
    feasQ = TrueQ[NumericQ[critF] && critF >= targetFid - 10^-12 && maxViol <= 10^-5];
    durVal = N[assoc["TotalDuration"]];
    objEvals = Lookup[assoc, "ObjectiveEvaluations", Lookup[assoc, "Evaluations", 0]];
    gradEvals = Lookup[assoc, "GradientEvaluations", Lookup[assoc, "Evaluations", 0]];
    nExp = estimateMatrixExponentials[res];
    <|
        "Duration" -> durVal,
        "NumSlices" -> Length[assoc["TimeSteps"]],
        "Fidelity" -> fid,
        "MeanFidelity" -> meanF,
        "SoftMinFidelity" -> softMinF,
        "MinFidelity" -> minF,
        "CriterionFidelity" -> critF,
        "MaxConstraintViolation" -> maxViol,
        "FeasibleQ" -> feasQ,
        "Iterations" -> Lookup[assoc, "Iterations", 0],
        "Evaluations" -> Lookup[assoc, "Evaluations", 0],
        "ObjectiveEvaluations" -> objEvals,
        "GradientEvaluations" -> gradEvals,
        "MatrixExponentials" -> nExp,
        "ElapsedTime" -> N[Lookup[assoc, "ElapsedTime", 0.]],
        "Status" -> Lookup[assoc, "Status", "Unknown"],
        "Result" -> res
    |>
];

evaluateValidationEnsemble[None | False, ___] := None;
evaluateValidationEnsemble[valEns_, res_OptimalControlResult, effCrit_String, targetFid_?NumericQ] := Module[
    {retargeted, comp, assoc, inArg, ev, meanF, minF, softMinF, maxF, critF, feasQ},
    If[res === None || !MatchQ[res, _OptimalControlResult], Return[None]];
    If[!MatchQ[valEns, _OptimalControlEnsemble],
        Message[TimeOptimalControlOptimize::badopt, ValidationEnsemble -> valEns, "expected OptimalControlEnsemble or None"];
        Return[$Failed]
    ];
    retargeted = ($OCSolverRegistry["RetargetEnsembleTiming"])[valEns, res["TimeSteps"]];
    comp = ($OCSolverRegistry["CompileIfEnsemble"])[retargeted];
    If[comp === $Failed || !AssociationQ[comp], Return[$Failed]];
    assoc = res[[1]];
    inArg = If[TrueQ[Lookup[comp, "HasPhaseFreedom", False]] && AssociationQ[Lookup[assoc, "OptimalPhases", None]],
        <|"Controls" -> res["Controls"], "OptimalPhases" -> assoc["OptimalPhases"]|>,
        res["Controls"]
    ];
    ev = ($OCSolverRegistry["EvaluateEnsemble"])[comp, inArg];
    If[ev === $Failed || !AssociationQ[ev], Return[$Failed]];
    meanF = N[ev["MeanFidelity"]];
    minF = N[ev["MinFidelity"]];
    softMinF = N[Lookup[ev, "SoftMinFidelity", meanF]];
    maxF = N[ev["MaxFidelity"]];
    critF = Switch[effCrit,
        "Minimum", minF,
        "SoftMin", softMinF,
        "Mean", meanF,
        _, meanF
    ];
    feasQ = TrueQ[NumericQ[critF] && critF >= targetFid - 10^-12];
    <|
        "NumMembers" -> comp["NumMembers"],
        "MeanFidelity" -> meanF,
        "SoftMinFidelity" -> softMinF,
        "MinFidelity" -> minF,
        "MaxFidelity" -> maxF,
        "CriterionFidelity" -> critF,
        "FeasibleQ" -> feasQ,
        "MemberFidelities" -> ev["MemberFidelities"],
        "Weights" -> comp["Weights"],
        "ParameterPoints" -> comp["ParameterPoints"],
        "GridShape" -> comp["GridShape"]
    |>
];

runDurationAttempt[kind_String, obj_, dtVec_List, prevRes_, grapeOpts_List,
                   rescaleMode_String, nAttempts_Integer, effCrit_String, targetFid_?NumericQ,
                   seedBase_, isAnchorPhase_?BooleanQ] := Module[
    {targetObj, newDur = Total[dtVec], critKey = $OCSolverRegistry["FidelityCriterionKey"],
     mCtrl, bestSumm = None, attempts = 0, totObj = 0, totGrad = 0, totExp = 0,
     aIdx, initChoice, seedChoice, attemptOpts, res, summ, altMode},
    targetObj = retargetObject[kind, obj, dtVec];
    mCtrl = If[prevRes =!= None, Dimensions[prevRes["Controls"]][[2]], 1];

    Do[
        attempts++;
        {initChoice, seedChoice} = Which[
            aIdx === 1 && prevRes =!= None && rescaleMode =!= "None",
                {buildWarmStartGuess[prevRes, dtVec, grapeOpts, rescaleMode, mCtrl],
                 If[seedBase === Automatic, 1000 + aIdx, seedBase + 100 * aIdx]},
            aIdx === 1 && prevRes === None,
                {InitialControls /. grapeOpts /. Options[GRAPEOptimize],
                 If[seedBase === Automatic, 1000 + aIdx, seedBase]},
            aIdx === 2 && prevRes =!= None && rescaleMode =!= "None",
                altMode = buildAltRescale[rescaleMode];
                {buildWarmStartGuess[prevRes, dtVec, grapeOpts, altMode, mCtrl],
                 If[seedBase === Automatic, 2000 + aIdx, seedBase + 200 * aIdx]},
            True,
                {Automatic,
                 If[seedBase === Automatic, 3000 + Round[10^6 * N[newDur]] + 37 * aIdx, seedBase + 300 * aIdx + Round[10^5 * N[newDur]]]}
        ];
        attemptOpts = Join[
            {InitialControls -> initChoice,
             RandomSeed -> seedChoice,
             TargetFidelity -> targetFid,
             critKey -> effCrit},
            If[!isAnchorPhase && aIdx === 1 && prevRes =!= None && kind === "Ensemble",
                {SoftMinSchedule -> None},
                {}
            ],
            FilterRules[grapeOpts, Except[InitialControls | RandomSeed | TargetFidelity | DurationPenalty | OptimizeDurations | FixedTotalDuration]]
        ];
        res = Quiet[Check[GRAPEOptimize[targetObj, Sequence @@ attemptOpts], $Failed]];
        If[MatchQ[res, _OptimalControlResult],
            summ = summarizeResult[res, kind, effCrit, targetFid];
            totObj += summ["ObjectiveEvaluations"];
            totGrad += summ["GradientEvaluations"];
            totExp += summ["MatrixExponentials"];
            If[bestSumm === None || (summ["FeasibleQ"] && !bestSumm["FeasibleQ"]) ||
               (summ["FeasibleQ"] === bestSumm["FeasibleQ"] && summ["CriterionFidelity"] > bestSumm["CriterionFidelity"]),
                bestSumm = summ
            ];
            If[summ["FeasibleQ"], Break[]]
        ],
        {aIdx, nAttempts}
    ];

    <|
        "Duration" -> newDur,
        "NumSlices" -> Length[dtVec],
        "BestSummary" -> bestSumm,
        "AttemptsUsed" -> attempts,
        "ObjectiveEvaluations" -> totObj,
        "GradientEvaluations" -> totGrad,
        "MatrixExponentials" -> totExp
    |>
];

(* Build initial parameter vector for Joint / JointPenalized / Hybrid solvers *)
extractJointInitVector[comp_Association, prevRes_] := If[MatchQ[prevRes, _OptimalControlResult],
    ($OCSolverRegistry["BuildInitialVector"])[comp, prevRes, 0.5, Automatic],
    None
];

compileForJoint[kind_String, obj_, dtVec_List, durMode_, durBnds_, durPen_, mergedOpts_List, effCrit_String: "Automatic"] := Module[
    {cleanOpts, targetObj, aggCur},
    aggCur = EnsembleAggregation /. mergedOpts /. Options[OptimalControlEnsemble];
    cleanOpts = Join[
        {OptimizeDurations -> durMode,
         DurationBounds -> durBnds,
         DurationPenalty -> durPen,
         FixedTotalDuration -> False},
        If[kind === "Ensemble" && effCrit === "Minimum" && !MatchQ[aggCur, {"SoftMin", _?NumericQ}],
            {EnsembleAggregation -> {"SoftMin", 20.}},
            {}
        ],
        FilterRules[mergedOpts, Except[OptimizeDurations | DurationBounds | DurationPenalty | FixedTotalDuration |
            MultiStart | TimeStepSchedule | SoftMinSchedule | StepMonitor | MonitorInterval |
            If[kind === "Ensemble" && effCrit === "Minimum" && !MatchQ[aggCur, {"SoftMin", _?NumericQ}], EnsembleAggregation, Unevaluated[Sequence[]]]]]
    ];
    targetObj = retargetObject[kind, obj, dtVec];
    If[kind === "Single",
        ($OCSolverRegistry["CompileSingle"])[
            OptimalControlProblem[targetObj[[1]], targetObj[[2]], targetObj[[3]], targetObj[[4]],
                Sequence @@ Join[FilterRules[cleanOpts, Options[OptimalControlProblem]], FilterRules[List @@ targetObj[[5 ;;]], Except[OptimizeDurations | DurationBounds | DurationPenalty | FixedTotalDuration]]]]
        ],
        ($OCSolverRegistry["CompileIfEnsemble"])[targetObj, Sequence @@ FilterRules[cleanOpts, Options[OptimalControlEnsemble]]]
    ]
];

checkJointCompatibility[comp_Association, kind_String] := Module[{members, hasRC},
    members = If[kind === "Ensemble", comp["Members"], {comp}];
    If[AnyTrue[members, Lookup[#, "DriftKind", "Static"] === "SliceList" &],
        Message[TimeOptimalControlOptimize::badopt, Method -> "Joint",
            "joint duration optimization does not support SliceList drift; use Method -> \"Progressive\""];
        Return[$Failed]
    ];
    If[TrueQ[comp["HasSlewBounds"]],
        Message[TimeOptimalControlOptimize::badopt, Method -> "Joint",
            "SlewRateBound depends coupledly on slice durations; use Method -> \"Progressive\""];
        Return[$Failed]
    ];
    hasRC = AnyTrue[members, Lookup[#, "RunningCostSpec", None] =!= None &];
    If[hasRC,
        Message[TimeOptimalControlOptimize::badopt, Method -> "Joint",
            "RunningCost cannot be combined with Method -> \"Joint\" AL; use Method -> \"Progressive\" or \"JointPenalized\""];
        Return[$Failed]
    ];
    True
];

runJointAugmentedLagrangian[kind_String, obj_, dt0_List, discSpec_Association, t0_?NumericQ, tMin_?NumericQ,
                            tTol_?NumericQ, targetFid_?NumericQ, effCrit_String, mergedOpts_List, seedOpt_] := Module[
    {tStart = AbsoluteTime[], optDurOpt, durBndsOpt, durMode, tMax, tRef, durBndResolved,
     comp, xVec, warmupRes = None, warmupAttempt = 0, totObj = 0, totGrad = 0, totExp = 0,
     hist = {}, bestFeasRes = None, bestFeasDur = Infinity, bestOverall = None, lastALSumm = None,
     evalBase, assembleFn, dispatchOpt, targetFidEff, mu = 2., rho = 50., prevG = -Infinity, prevDur,
     maxOuter = 14, stallCount = 0, evalALFn, maxIterUser, optsIter, optRes, resIter, summIter,
     durCur, fidCur, gCur, tFinal, polishRes, polishSumm},

    optDurOpt = OptimizeDurations /. mergedOpts /. Options[GRAPEOptimize];
    durBndsOpt = DurationBounds /. mergedOpts /. Options[GRAPEOptimize];
    durMode = If[optDurOpt === False || optDurOpt === None || optDurOpt === Automatic, "Total", optDurOpt];
    tMax = If[MatchQ[durBndsOpt, {_?NumericQ, _?NumericQ}] && TrueQ[durBndsOpt[[2]] >= t0],
        N[Re[durBndsOpt[[2]]]],
        Max[1.25 * t0, 10^-6]
    ];
    durBndResolved = Which[
        durMode === "Total", {tMin, tMax},
        durBndsOpt === Automatic, {Max[10^-12, 0.05 * (tMin / Length[dt0])], Max[2. * Mean[dt0], tMax]},
        True, durBndsOpt
    ];

    comp = compileForJoint[kind, obj, dt0, durMode, durBndResolved, 0., mergedOpts, effCrit];
    If[comp === $Failed || !AssociationQ[comp], Return[$Failed]];
    If[checkJointCompatibility[comp, kind] === $Failed, Return[$Failed]];

    evalBase = If[kind === "Ensemble",
        Function[v, ($OCSolverRegistry["EvaluateEnsemble"])[comp, v]],
        Function[v, EvaluateOptimalControl[comp, v]]
    ];
    assembleFn = If[kind === "Ensemble", $OCSolverRegistry["AssembleEnsembleResult"], $OCSolverRegistry["AssembleSingleResult"]];
    dispatchOpt = $OCSolverRegistry["DispatchOptimizer"];

    (* Warm start at initial duration so AL starts inside or near the feasible basin *)
    Module[{warmAttempt, warmDt},
        warmDt = buildDurationGrid[t0, dt0, discSpec];
        If[warmDt === $Failed, Return[$Failed]];
        warmAttempt = runDurationAttempt[kind, obj, warmDt, None, mergedOpts, "None", 2, effCrit, targetFid, seedOpt, True];
        totObj += warmAttempt["ObjectiveEvaluations"];
        totGrad += warmAttempt["GradientEvaluations"];
        totExp += warmAttempt["MatrixExponentials"];
        warmupAttempt = warmAttempt["AttemptsUsed"];
        If[AssociationQ[warmAttempt["BestSummary"]],
            warmupRes = warmAttempt["BestSummary", "Result"];
            bestOverall = warmAttempt["BestSummary"];
            AppendTo[hist, Join[KeyDrop[bestOverall, "Result"], <|"Phase" -> "Warmup", "OuterIteration" -> 0, "AttemptsUsed" -> warmupAttempt|>]];
            If[bestOverall["FeasibleQ"],
                bestFeasRes = warmupRes;
                bestFeasDur = bestOverall["Duration"]
            ]
        ]
    ];

    xVec = If[MatchQ[warmupRes, _OptimalControlResult],
        extractJointInitVector[comp, warmupRes],
        ($OCSolverRegistry["BuildInitialVector"])[comp, InitialControls /. mergedOpts /. Options[GRAPEOptimize],
            InitialScale /. mergedOpts /. Options[GRAPEOptimize], seedOpt]
    ];
    If[xVec === $Failed || !VectorQ[xVec, NumericQ], Return[$Failed]];

    tRef = Max[t0, 10^-12];
    targetFidEff = targetFid;
    prevDur = t0;
    maxIterUser = Max[80, MaxIterations /. mergedOpts /. Options[TimeOptimalControlOptimize]];

    Do[
        optsIter = Join[
            {TargetFidelity -> None,
             MaxIterations -> Which[outer <= 3, Min[maxIterUser, 60], outer <= 6, Min[maxIterUser, 100], True, maxIterUser],
             GradientTolerance -> Which[outer <= 3, 10^-4, outer <= 6, 10^-5, True, 10^-6],
             FunctionTolerance -> 10^-12,
             StepTolerance -> 10^-11,
             RecordHistory -> True,
             MultiStart -> 1,
             TimeStepSchedule -> None,
             SoftMinSchedule -> None},
            FilterRules[mergedOpts, Except[TargetFidelity | MaxIterations | GradientTolerance | FunctionTolerance | StepTolerance | MultiStart | TimeStepSchedule | SoftMinSchedule | DurationPenalty]]
        ];
        evalALFn = With[{muK = mu, rhoK = rho, tfK = targetFidEff, trK = tRef, eB = evalBase, hasPhi = TrueQ[comp["HasPhaseFreedom"]]},
            Function[v,
                Module[{base = eB[v], fidB, gB, tVal, pS, pNS, qVal, sObjAL, objAL, gCtrlAL, gSliceAL, gDurAL, gPhiAL, nDurA},
                    If[base === $Failed || !AssociationQ[base], Return[$Failed]];
                    fidB = N[base["Fidelity"]];
                    gB = fidB - tfK;
                    tVal = N[base["TotalDuration"]];
                    pS = N[Lookup[base, "SmoothPenalty", 0.]];
                    pNS = N[Lookup[base, "NonSmoothPenalty", 0.]];
                    qVal = Max[0., muK - rhoK * gB];
                    sObjAL = -(tVal / trK) - pS - (qVal^2 - muK^2) / (2. * rhoK);
                    objAL = sObjAL - pNS;
                    gCtrlAL = Developer`ToPackedArray[qVal * base["FidelityControlGradient"] + base["SmoothPenaltyControlGradient"]];
                    gSliceAL = Developer`ToPackedArray[qVal * base["FidelitySliceGradient"] + (base["SliceGradient"] - base["FidelitySliceGradient"])];
                    nDurA = Length[base["FidelityActiveDurationGradient"]];
                    gDurAL = Developer`ToPackedArray[ConstantArray[-1. / trK, nDurA] + qVal * base["FidelityActiveDurationGradient"] + base["SmoothPenaltyActiveDurationGradient"]];
                    gPhiAL = If[hasPhi, Developer`ToPackedArray[qVal * base["PhaseGradient"]], {}];
                    Join[base, <|
                        "Objective" -> objAL,
                        "SmoothObjective" -> sObjAL,
                        "Gradient" -> Developer`ToPackedArray[Join[Flatten[gCtrlAL], gDurAL, gPhiAL]],
                        "ControlGradient" -> gCtrlAL,
                        "SliceGradient" -> gSliceAL,
                        "ActiveDurationGradient" -> gDurAL,
                        "PhaseGradient" -> gPhiAL
                    |>]
                ]
            ]
        ];
        optRes = dispatchOpt[evalALFn, xVec, comp, optsIter];
        If[optRes === $Failed || !AssociationQ[optRes], Break[]];
        xVec = optRes["ParametersVec"];
        resIter = assembleFn[comp, optRes];
        If[resIter === $Failed, Break[]];
        summIter = summarizeResult[resIter, kind, effCrit, targetFid];
        lastALSumm = summIter;
        totObj += Lookup[optRes, "ObjectiveEvaluations", summIter["Evaluations"]];
        totGrad += Lookup[optRes, "GradientEvaluations", summIter["Evaluations"]];
        totExp += estimateMatrixExponentials[resIter];
        If[bestOverall === None || (summIter["FeasibleQ"] && !bestOverall["FeasibleQ"]) ||
           (summIter["FeasibleQ"] && bestOverall["FeasibleQ"] && summIter["Duration"] < bestOverall["Duration"]) ||
           (!bestOverall["FeasibleQ"] && summIter["CriterionFidelity"] > bestOverall["CriterionFidelity"]),
            bestOverall = summIter
        ];
        If[summIter["FeasibleQ"] && summIter["Duration"] < bestFeasDur,
            bestFeasDur = summIter["Duration"];
            bestFeasRes = resIter
        ];
        durCur = summIter["Duration"];
        fidCur = summIter["Fidelity"];
        gCur = fidCur - targetFidEff;
        AppendTo[hist, Join[KeyDrop[summIter, "Result"], <|
            "Phase" -> "JointAL",
            "OuterIteration" -> outer,
            "Mu" -> mu,
            "Rho" -> rho,
            "EffectiveTargetFidelity" -> targetFidEff
        |>]];

        If[summIter["FeasibleQ"] && outer >= 2 && Abs[durCur - prevDur] <= tTol && gCur >= -10^-6,
            Break[]
        ];
        If[!summIter["FeasibleQ"] && outer >= 3 && Abs[durCur - prevDur] <= 0.1 * tTol && Abs[gCur - prevG] <= 10^-6,
            stallCount++,
            stallCount = 0
        ];
        If[stallCount >= 3, Break[]];

        mu = Max[0., mu - rho * gCur];
        If[gCur < -10^-7 || (Abs[prevG] < Infinity && Abs[gCur] > 0.3 * Abs[prevG]),
            rho = Min[10^9, 4. * rho]
        ];
        If[kind === "Ensemble" && effCrit === "Minimum",
            With[{upperCap = 1. - 0.15 * (1. - targetFid)},
                Which[
                    summIter["SoftMinFidelity"] >= targetFid - 10^-4 && summIter["MinFidelity"] < targetFid,
                        With[{gap = Max[0., targetFid - summIter["MinFidelity"]]},
                            targetFidEff = Min[upperCap, targetFidEff + 0.6 * gap];
                            rho = Min[Max[rho, 200.], 2500.];
                            mu = Min[mu, 20.]
                        ],
                    summIter["MinFidelity"] > targetFid + 0.002 && targetFidEff > targetFid,
                        targetFidEff = Max[targetFid, 0.5 * (targetFidEff + targetFid + (summIter["SoftMinFidelity"] - summIter["MinFidelity"]))]
                ]
            ]
        ];
        prevDur = durCur;
        prevG = gCur,
        {outer, maxOuter}
    ];

    (* Seal near-boundary final AL iterate with a fixed-duration polish at FinalALDuration (or +0.1*tTol bump) *)
    If[AssociationQ[lastALSumm] && !TrueQ[lastALSumm["FeasibleQ"]] && lastALSumm["CriterionFidelity"] >= targetFid - 2.*^-3,
        tFinal = N[lastALSumm["Duration"]];
        Module[{candDurs = DeleteDuplicates[{tFinal, Min[t0, tFinal + 0.1 * tTol]}], pAtt},
            Do[
                pAtt = runDurationAttempt[kind, obj, buildDurationGrid[dTry, dt0, discSpec],
                    lastALSumm["Result"], mergedOpts, "NormalizedTime", 1, effCrit, targetFid, seedOpt, False];
                totObj += pAtt["ObjectiveEvaluations"];
                totGrad += pAtt["GradientEvaluations"];
                totExp += pAtt["MatrixExponentials"];
                If[AssociationQ[pAtt["BestSummary"]],
                    polishSumm = pAtt["BestSummary"];
                    AppendTo[hist, Join[KeyDrop[polishSumm, "Result"], <|"Phase" -> "Polish", "AttemptsUsed" -> pAtt["AttemptsUsed"]|>]];
                    If[polishSumm["FeasibleQ"] && polishSumm["Duration"] < bestFeasDur,
                        bestFeasDur = polishSumm["Duration"];
                        bestFeasRes = polishSumm["Result"];
                        bestOverall = polishSumm;
                        lastALSumm = polishSumm;
                        Break[]
                    ]
                ],
                {dTry, candDurs}
            ]
        ]
    ];

    <|
        "BestFeasibleResult" -> bestFeasRes,
        "BestFeasibleDuration" -> If[bestFeasRes =!= None, bestFeasDur, None],
        "BestOverallSummary" -> bestOverall,
        "FinalALDuration" -> If[AssociationQ[lastALSumm], lastALSumm["Duration"], None],
        "FinalALSummary" -> lastALSumm,
        "History" -> hist,
        "TotalAttempts" -> warmupAttempt + Length[Select[hist, #["Phase"] === "JointAL" &]],
        "TotalObjectiveEvaluations" -> totObj,
        "TotalGradientEvaluations" -> totGrad,
        "TotalMatrixExponentials" -> totExp,
        "ElapsedTime" -> AbsoluteTime[] - tStart
    |>
];

formatTOOutput[methodName_String, targetFid_?NumericQ, effCrit_String, discOpt_, rescaleOpt_,
               bestFeasRes_, bestFeasDur_, bestOverall_, bracket_, hist_List,
               totAtt_Integer, totObj_Integer, totGrad_Integer, totExp_Integer, elapsed_?NumericQ,
               valEns_, explicitStatus_: None] := Module[
    {feasQ = bestFeasRes =!= None, bestRes, summ, valRes, statusStr},
    bestRes = If[feasQ, bestFeasRes, If[AssociationQ[bestOverall], bestOverall["Result"], None]];
    summ = Which[
        feasQ && MatchQ[bestFeasRes, _OptimalControlResult],
            summarizeResult[bestFeasRes, Lookup[bestFeasRes[[1]], "Kind", "Single"], effCrit, targetFid],
        AssociationQ[bestOverall],
            bestOverall,
        True,
            <|"Fidelity" -> Indeterminate, "MeanFidelity" -> Indeterminate, "SoftMinFidelity" -> Indeterminate,
              "MinFidelity" -> Indeterminate, "CriterionFidelity" -> Indeterminate|>
    ];
    valRes = If[valEns =!= None && valEns =!= False && MatchQ[bestRes, _OptimalControlResult],
        evaluateValidationEnsemble[valEns, bestRes, effCrit, targetFid],
        None
    ];
    statusStr = Which[
        explicitStatus =!= None, explicitStatus,
        feasQ, "FeasiblePulseFound",
        True, "TargetNotAttained"
    ];
    <|
        "Method" -> methodName,
        "ShortestFeasibleDuration" -> If[feasQ, N[bestFeasDur], None],
        "ShortestFeasibleResult" -> bestFeasRes,
        "BestResult" -> bestRes,
        "FeasibleQ" -> feasQ,
        "Status" -> statusStr,
        "TargetFidelity" -> N[targetFid],
        "FidelityCriterion" -> effCrit,
        "DiscretizationStrategy" -> discOpt,
        "WaveformRescaling" -> rescaleOpt,
        "AchievedFidelity" -> summ["Fidelity"],
        "MeanFidelity" -> summ["MeanFidelity"],
        "SoftMinFidelity" -> summ["SoftMinFidelity"],
        "MinFidelity" -> summ["MinFidelity"],
        "CriterionFidelity" -> summ["CriterionFidelity"],
        "Bracket" -> bracket,
        "History" -> hist,
        "TotalAttempts" -> totAtt,
        "TotalObjectiveEvaluations" -> totObj,
        "TotalGradientEvaluations" -> totGrad,
        "TotalMatrixExponentials" -> totExp,
        "ElapsedTime" -> N[elapsed],
        "ValidationResult" -> valRes,
        "OptimalityDisclaimer" -> $OptimalityDisclaimer
    |>
];

TimeOptimalControlOptimize[input_, opts___Rule] := Module[
    {norm, kind, obj, ownOpts, mergedOpts, grapeOpts,
     methodOpt, targetFidOpt, initDurOpt, minDurOpt, tolOpt, redOpt,
     critOpt, discOpt, rescaleOpt, attOpt, valEnsOpt, seedOpt, verboseOpt,
     defaultNt, dt0, t0, tMin, tTol, rFac, effCrit, isSliceListQ, discSpec, attSpec,
     kBase, kBound, tStart = AbsoluteTime[],
     totAtt = 0, totObj = 0, totGrad = 0, totExp = 0, hist = {},
     bestFeasRes = None, bestFeasDur = Infinity, bestOverall = None,
     dtCand, attRes, summ, prevRes, tLow, tHigh, tMid, nBisect,
     jointOut, hybridRef, penOpt, penList, lamVal, curRes, optDurOpt, durMode, compCheck},

    norm = normalizeTOInput[input];
    If[norm === $Failed,
        Message[TimeOptimalControlOptimize::badprob, input]; Return[$Failed]
    ];
    {kind, obj, ownOpts} = norm;
    mergedOpts = Join[{opts}, ownOpts, Options[TimeOptimalControlOptimize], Options[GRAPEOptimize]];
    grapeOpts = Prepend[
        FilterRules[mergedOpts, Except[Alternatives @@ Append[$TOOnlyOptionNames, Method]]],
        Method -> (Method /. ownOpts /. Options[GRAPEOptimize])
    ];

    {methodOpt, targetFidOpt, initDurOpt, minDurOpt, tolOpt, redOpt,
     critOpt, discOpt, rescaleOpt, attOpt, valEnsOpt, seedOpt, verboseOpt} =
        {Method, TargetFidelity, InitialDuration, MinimumDuration, DurationTolerance, DurationReductionFactor,
         FidelityCriterion, DiscretizationStrategy, WaveformRescaling, AttemptsPerDuration,
         ValidationEnsemble, RandomSeed, Verbose} /. mergedOpts;
    If[methodOpt === Automatic || !MemberQ[{"Progressive", "Joint", "JointPenalized", "Hybrid"}, methodOpt],
        If[MemberQ[{"Spectral", "BlockMatrixExp", "FirstOrder", "StateVector", "Liouville"}, methodOpt],
            grapeOpts = Prepend[FilterRules[grapeOpts, Except[Method]], Method -> methodOpt];
            methodOpt = "Progressive"
        ]
    ];

    If[!MemberQ[{"Progressive", "Joint", "JointPenalized", "Hybrid"}, methodOpt],
        Message[TimeOptimalControlOptimize::badopt, Method -> methodOpt,
            "expected \"Progressive\", \"Joint\", \"JointPenalized\", or \"Hybrid\""];
        Return[$Failed]
    ];
    If[!NumericQ[targetFidOpt] || !TrueQ[0. < Re[targetFidOpt] <= 1.],
        Message[TimeOptimalControlOptimize::badopt, TargetFidelity -> targetFidOpt, "expected number in (0, 1]"];
        Return[$Failed]
    ];
    If[!MemberQ[{"NormalizedTime", "AreaPreserving", "None"}, rescaleOpt],
        Message[TimeOptimalControlOptimize::badopt, WaveformRescaling -> rescaleOpt,
            "expected \"NormalizedTime\", \"AreaPreserving\", or \"None\""];
        Return[$Failed]
    ];
    effCrit = resolveEffectiveCriterion[kind, critOpt, mergedOpts];
    If[effCrit === $Failed,
        Message[TimeOptimalControlOptimize::badopt, FidelityCriterion -> critOpt,
            "expected Automatic, \"Fidelity\", \"Mean\", \"SoftMin\", or \"Minimum\""];
        Return[$Failed]
    ];

    defaultNt = TimeSteps /. mergedOpts /. Options[GRAPEOptimize];
    dt0 = extractBaselineTiming[kind, obj, defaultNt];
    If[dt0 === $Failed || !VectorQ[dt0, NumericQ],
        Message[TimeOptimalControlOptimize::baddur, dt0]; Return[$Failed]
    ];
    t0 = If[initDurOpt === Automatic, N[Total[dt0]], N[Re[initDurOpt]]];
    If[!NumericQ[t0] || !TrueQ[t0 > 0.],
        Message[TimeOptimalControlOptimize::baddur, InitialDuration -> initDurOpt]; Return[$Failed]
    ];
    tMin = If[minDurOpt === Automatic, 0.05 * t0, N[Re[minDurOpt]]];
    If[!NumericQ[tMin] || !TrueQ[0. < tMin <= t0],
        Message[TimeOptimalControlOptimize::baddur, MinimumDuration -> minDurOpt]; Return[$Failed]
    ];
    tTol = If[tolOpt === Automatic, Max[10^-6, 0.005 * t0], N[Re[tolOpt]]];
    If[!NumericQ[tTol] || !TrueQ[tTol > 0.],
        Message[TimeOptimalControlOptimize::baddur, DurationTolerance -> tolOpt]; Return[$Failed]
    ];
    rFac = N[Re[redOpt]];
    If[!NumericQ[rFac] || !TrueQ[0. < rFac < 1.],
        Message[TimeOptimalControlOptimize::badopt, DurationReductionFactor -> redOpt, "expected number in (0, 1)"];
        Return[$Failed]
    ];

    isSliceListQ = hasSliceListDriftQ[kind, obj];
    discSpec = validateDiscretization[discOpt, dt0, isSliceListQ];
    If[discSpec === $Failed, Return[$Failed]];
    attSpec = parseAttemptBudget[attOpt];
    If[attSpec === $Failed, Return[$Failed]];
    {kBase, kBound} = attSpec;

    Switch[methodOpt,
        "Progressive",
            (* Phase 1: Establish feasible pulse at InitialDuration *)
            dtCand = buildDurationGrid[t0, dt0, discSpec];
            If[dtCand === $Failed, Return[$Failed]];
            attRes = runDurationAttempt[kind, obj, dtCand, None, grapeOpts, "None",
                Max[kBase, kBound], effCrit, targetFidOpt, seedOpt, True];
            totAtt += attRes["AttemptsUsed"];
            totObj += attRes["ObjectiveEvaluations"];
            totGrad += attRes["GradientEvaluations"];
            totExp += attRes["MatrixExponentials"];
            summ = attRes["BestSummary"];
            If[AssociationQ[summ],
                bestOverall = summ;
                AppendTo[hist, Join[KeyDrop[summ, "Result"], <|
                    "Phase" -> "InitialAnchor", "AttemptsUsed" -> attRes["AttemptsUsed"],
                    "TotalObjectiveEvaluations" -> attRes["ObjectiveEvaluations"],
                    "TotalGradientEvaluations" -> attRes["GradientEvaluations"],
                    "TotalMatrixExponentials" -> attRes["MatrixExponentials"]
                |>]]
            ];
            If[!AssociationQ[summ] || !TrueQ[summ["FeasibleQ"]],
                Message[TimeOptimalControlOptimize::infeasibleinit, t0, N[targetFidOpt],
                    If[AssociationQ[summ], summ["CriterionFidelity"], Indeterminate]];
                Return[formatTOOutput["Progressive", targetFidOpt, effCrit, discOpt, rescaleOpt,
                    None, None, bestOverall, {tMin, t0}, hist,
                    totAtt, totObj, totGrad, totExp, AbsoluteTime[] - tStart, valEnsOpt, "InitialDurationInfeasible"]]
            ];
            bestFeasRes = summ["Result"];
            bestFeasDur = summ["Duration"];
            prevRes = bestFeasRes;
            tHigh = bestFeasDur;
            tLow = tMin;

            (* Phase 2: Adaptive duration reduction using warm starts + boundary retry escalation *)
            While[bestFeasDur > tMin + tTol,
                tMid = Max[tMin, rFac * bestFeasDur];
                If[tMid >= bestFeasDur - 10^-15 * Max[1., bestFeasDur], Break[]];
                dtCand = buildDurationGrid[tMid, dt0, discSpec];
                If[dtCand === $Failed, Return[$Failed]];
                attRes = runDurationAttempt[kind, obj, dtCand, prevRes, grapeOpts, rescaleOpt,
                    Max[kBase, kBound], effCrit, targetFidOpt, seedOpt, False];
                totAtt += attRes["AttemptsUsed"];
                totObj += attRes["ObjectiveEvaluations"];
                totGrad += attRes["GradientEvaluations"];
                totExp += attRes["MatrixExponentials"];
                summ = attRes["BestSummary"];
                If[AssociationQ[summ],
                    AppendTo[hist, Join[KeyDrop[summ, "Result"], <|
                        "Phase" -> "Reduction", "AttemptsUsed" -> attRes["AttemptsUsed"],
                        "TotalObjectiveEvaluations" -> attRes["ObjectiveEvaluations"],
                        "TotalGradientEvaluations" -> attRes["GradientEvaluations"],
                        "TotalMatrixExponentials" -> attRes["MatrixExponentials"]
                    |>]];
                    If[!bestOverall["FeasibleQ"] && summ["CriterionFidelity"] > bestOverall["CriterionFidelity"],
                        bestOverall = summ
                    ];
                    If[TrueQ[summ["FeasibleQ"]],
                        If[summ["Duration"] < bestFeasDur,
                            bestFeasDur = summ["Duration"];
                            bestFeasRes = summ["Result"];
                            bestOverall = summ
                        ];
                        prevRes = bestFeasRes;
                        tHigh = bestFeasDur;
                        If[tMid <= tMin + 10^-13 * Max[1., tMin], tLow = tMin; Break[]],
                        tLow = tMid;
                        Break[]
                    ],
                    tLow = tMid;
                    Break[]
                ]
            ];

            (* Phase 3: Safeguarded bracket refinement (never assume strict monotonicity) *)
            nBisect = 0;
            While[tHigh - tLow > tTol && nBisect < 25,
                nBisect++;
                tMid = 0.5 * (tLow + tHigh);
                If[tMid <= tLow + 10^-14 * Max[1., tMid] || tMid >= tHigh - 10^-14 * Max[1., tMid], Break[]];
                dtCand = buildDurationGrid[tMid, dt0, discSpec];
                If[dtCand === $Failed, Return[$Failed]];
                attRes = runDurationAttempt[kind, obj, dtCand, bestFeasRes, grapeOpts, rescaleOpt,
                    kBound, effCrit, targetFidOpt, seedOpt, False];
                totAtt += attRes["AttemptsUsed"];
                totObj += attRes["ObjectiveEvaluations"];
                totGrad += attRes["GradientEvaluations"];
                totExp += attRes["MatrixExponentials"];
                summ = attRes["BestSummary"];
                If[AssociationQ[summ],
                    AppendTo[hist, Join[KeyDrop[summ, "Result"], <|
                        "Phase" -> "BracketRefinement", "AttemptsUsed" -> attRes["AttemptsUsed"],
                        "TotalObjectiveEvaluations" -> attRes["ObjectiveEvaluations"],
                        "TotalGradientEvaluations" -> attRes["GradientEvaluations"],
                        "TotalMatrixExponentials" -> attRes["MatrixExponentials"]
                    |>]];
                    If[TrueQ[summ["FeasibleQ"]],
                        tHigh = Min[tHigh, summ["Duration"]];
                        If[summ["Duration"] < bestFeasDur,
                            bestFeasDur = summ["Duration"];
                            bestFeasRes = summ["Result"];
                            bestOverall = summ
                        ],
                        tLow = Max[tLow, summ["Duration"]]
                    ],
                    tLow = tMid
                ]
            ];
            formatTOOutput["Progressive", targetFidOpt, effCrit, discOpt, rescaleOpt,
                bestFeasRes, bestFeasDur, bestOverall, {tLow, bestFeasDur}, hist,
                totAtt, totObj, totGrad, totExp, AbsoluteTime[] - tStart, valEnsOpt],

        "Joint",
            jointOut = runJointAugmentedLagrangian[kind, obj, dt0, discSpec, t0, tMin, tTol,
                targetFidOpt, effCrit, grapeOpts, seedOpt];
            If[jointOut === $Failed, Return[$Failed]];
            formatTOOutput["Joint", targetFidOpt, effCrit, discOpt, rescaleOpt,
                jointOut["BestFeasibleResult"], jointOut["BestFeasibleDuration"],
                jointOut["BestOverallSummary"], {tMin, If[jointOut["BestFeasibleDuration"] =!= None, jointOut["BestFeasibleDuration"], t0]},
                jointOut["History"], jointOut["TotalAttempts"],
                jointOut["TotalObjectiveEvaluations"], jointOut["TotalGradientEvaluations"],
                jointOut["TotalMatrixExponentials"], jointOut["ElapsedTime"], valEnsOpt],

        "Hybrid",
            jointOut = runJointAugmentedLagrangian[kind, obj, dt0, discSpec, t0, tMin, tTol,
                targetFidOpt, effCrit, grapeOpts, seedOpt];
            If[jointOut === $Failed, Return[$Failed]];
            bestFeasRes = jointOut["BestFeasibleResult"];
            bestFeasDur = If[jointOut["BestFeasibleDuration"] =!= None, jointOut["BestFeasibleDuration"], Infinity];
            bestOverall = jointOut["BestOverallSummary"];
            hist = jointOut["History"];
            totAtt = jointOut["TotalAttempts"];
            totObj = jointOut["TotalObjectiveEvaluations"];
            totGrad = jointOut["TotalGradientEvaluations"];
            totExp = jointOut["TotalMatrixExponentials"];
            Module[{tCenter, finalALSumm = Lookup[jointOut, "FinalALSummary", None]},
                tCenter = If[NumericQ[Lookup[jointOut, "FinalALDuration", None]],
                    jointOut["FinalALDuration"],
                    If[bestFeasRes =!= None, bestFeasDur, t0]
                ];
                prevRes = If[AssociationQ[finalALSumm] && MatchQ[finalALSumm["Result"], _OptimalControlResult],
                    finalALSumm["Result"],
                    If[bestFeasRes =!= None, bestFeasRes, If[AssociationQ[bestOverall], bestOverall["Result"], None]]
                ];
                tHigh = If[bestFeasRes =!= None, bestFeasDur, t0];
                tLow = Max[tMin, 0.85 * tCenter];
                If[(bestFeasRes === None || tCenter < bestFeasDur - 0.5 * tTol) && NumericQ[tCenter],
                    Module[{probeDur = Min[If[bestFeasRes =!= None, bestFeasDur, t0], Max[tMin, 1.02 * tCenter]]},
                        dtCand = buildDurationGrid[probeDur, dt0, discSpec];
                        attRes = runDurationAttempt[kind, obj, dtCand, prevRes, grapeOpts, rescaleOpt,
                            kBound, effCrit, targetFidOpt, seedOpt, False];
                        totAtt += attRes["AttemptsUsed"]; totObj += attRes["ObjectiveEvaluations"];
                        totGrad += attRes["GradientEvaluations"]; totExp += attRes["MatrixExponentials"];
                        If[AssociationQ[attRes["BestSummary"]],
                            summ = attRes["BestSummary"];
                            AppendTo[hist, Join[KeyDrop[summ, "Result"], <|"Phase" -> "HybridProbe", "AttemptsUsed" -> attRes["AttemptsUsed"]|>]];
                            If[summ["FeasibleQ"],
                                If[summ["Duration"] < bestFeasDur,
                                    bestFeasDur = summ["Duration"]; bestFeasRes = summ["Result"]; bestOverall = summ
                                ];
                                tHigh = bestFeasDur; tLow = Max[tMin, 0.85 * bestFeasDur],
                                tLow = Max[tLow, summ["Duration"]]
                            ]
                        ]
                    ]
                ];
                If[bestFeasRes =!= None && tLow < tHigh - tTol,
                    (* Verify tLow is below threshold; if surprisingly feasible, step down once more *)
                    Do[
                        dtCand = buildDurationGrid[tLow, dt0, discSpec];
                        attRes = runDurationAttempt[kind, obj, dtCand, bestFeasRes, grapeOpts, rescaleOpt,
                            Max[1, kBase], effCrit, targetFidOpt, seedOpt, False];
                        totAtt += attRes["AttemptsUsed"]; totObj += attRes["ObjectiveEvaluations"];
                        totGrad += attRes["GradientEvaluations"]; totExp += attRes["MatrixExponentials"];
                        If[AssociationQ[attRes["BestSummary"]],
                            summ = attRes["BestSummary"];
                            AppendTo[hist, Join[KeyDrop[summ, "Result"], <|"Phase" -> "HybridProbe", "AttemptsUsed" -> attRes["AttemptsUsed"]|>]];
                            If[summ["FeasibleQ"],
                                If[summ["Duration"] < bestFeasDur,
                                    bestFeasDur = summ["Duration"]; bestFeasRes = summ["Result"]; bestOverall = summ
                                ];
                                tHigh = bestFeasDur;
                                If[tLow <= tMin + tTol, Break[]];
                                tLow = Max[tMin, 0.8 * tLow],
                                Break[]
                            ],
                            Break[]
                        ],
                        {4}
                    ];
                    nBisect = 0;
                    While[tHigh - tLow > tTol && nBisect < 15,
                        nBisect++;
                        tMid = 0.5 * (tLow + tHigh);
                        If[tMid <= tLow + 10^-14 * Max[1., tMid] || tMid >= tHigh - 10^-14 * Max[1., tMid], Break[]];
                        dtCand = buildDurationGrid[tMid, dt0, discSpec];
                        attRes = runDurationAttempt[kind, obj, dtCand, bestFeasRes, grapeOpts, rescaleOpt,
                            kBound, effCrit, targetFidOpt, seedOpt, False];
                        totAtt += attRes["AttemptsUsed"]; totObj += attRes["ObjectiveEvaluations"];
                        totGrad += attRes["GradientEvaluations"]; totExp += attRes["MatrixExponentials"];
                        If[AssociationQ[attRes["BestSummary"]],
                            summ = attRes["BestSummary"];
                            AppendTo[hist, Join[KeyDrop[summ, "Result"], <|"Phase" -> "HybridRefinement", "AttemptsUsed" -> attRes["AttemptsUsed"]|>]];
                            If[summ["FeasibleQ"],
                                tHigh = Min[tHigh, summ["Duration"]];
                                If[summ["Duration"] < bestFeasDur,
                                    bestFeasDur = summ["Duration"]; bestFeasRes = summ["Result"]; bestOverall = summ
                                ],
                                tLow = Max[tLow, summ["Duration"]]
                            ],
                            tLow = tMid
                        ]
                    ]
                ]
            ];
            formatTOOutput["Hybrid", targetFidOpt, effCrit, discOpt, rescaleOpt,
                bestFeasRes, If[bestFeasRes =!= None, bestFeasDur, None], bestOverall, {tLow, tHigh}, hist,
                totAtt, totObj, totGrad, totExp, AbsoluteTime[] - tStart, valEnsOpt],

        "JointPenalized",
            optDurOpt = OptimizeDurations /. mergedOpts /. Options[GRAPEOptimize];
            durMode = If[optDurOpt === False || optDurOpt === None || optDurOpt === Automatic, "Total", optDurOpt];
            penOpt = DurationPenalty /. mergedOpts /. Options[GRAPEOptimize];
            penList = Which[
                ListQ[penOpt] && VectorQ[penOpt, NumericQ] && AllTrue[penOpt, # >= 0. &], N[penOpt],
                NumericQ[penOpt] && TrueQ[penOpt > 0.], {N[Re[penOpt]]},
                True, {0.001, 0.005, 0.02, 0.08, 0.25}
            ];
            dtCand = buildDurationGrid[t0, dt0, discSpec];
            If[dtCand === $Failed, Return[$Failed]];
            compCheck = compileForJoint[kind, obj, dtCand, durMode, {tMin, Max[1.5 * t0, 10^-6]}, First[penList], grapeOpts];
            If[compCheck === $Failed, Return[$Failed]];
            If[TrueQ[compCheck["HasSlewBounds"]] || AnyTrue[If[kind === "Ensemble", compCheck["Members"], {compCheck}], Lookup[#, "DriftKind", "Static"] === "SliceList" &],
                Message[TimeOptimalControlOptimize::badopt, Method -> "JointPenalized",
                    "SlewRateBound and SliceList drift cannot be combined with OptimizeDurations"];
                Return[$Failed]
            ];
            prevRes = None;
            Do[
                totAtt++;
                curRes = Quiet[Check[GRAPEOptimize[retargetObject[kind, obj, dtCand],
                    OptimizeDurations -> durMode,
                    DurationBounds -> If[durMode === "Total", {tMin, Max[1.5 * t0, 10^-6]}, DurationBounds /. mergedOpts /. Options[GRAPEOptimize]],
                    DurationPenalty -> {lamVal, t0},
                    FixedTotalDuration -> False,
                    TargetFidelity -> None,
                    InitialControls -> If[prevRes =!= None, prevRes, InitialControls /. mergedOpts /. Options[GRAPEOptimize]],
                    Sequence @@ FilterRules[grapeOpts, Except[OptimizeDurations | DurationBounds | DurationPenalty | FixedTotalDuration | TargetFidelity | InitialControls]]
                ], $Failed]];
                If[MatchQ[curRes, _OptimalControlResult],
                    summ = summarizeResult[curRes, kind, effCrit, targetFidOpt];
                    totObj += summ["ObjectiveEvaluations"];
                    totGrad += summ["GradientEvaluations"];
                    totExp += summ["MatrixExponentials"];
                    AppendTo[hist, Join[KeyDrop[summ, "Result"], <|"Phase" -> "JointPenalized", "DurationPenalty" -> lamVal, "AttemptsUsed" -> 1|>]];
                    prevRes = curRes;
                    If[bestOverall === None || (summ["FeasibleQ"] && !bestOverall["FeasibleQ"]) ||
                       (summ["FeasibleQ"] && bestOverall["FeasibleQ"] && summ["Duration"] < bestOverall["Duration"]) ||
                       (!bestOverall["FeasibleQ"] && summ["CriterionFidelity"] > bestOverall["CriterionFidelity"]),
                        bestOverall = summ
                    ];
                    If[summ["FeasibleQ"] && summ["Duration"] < bestFeasDur,
                        bestFeasDur = summ["Duration"]; bestFeasRes = curRes
                    ]
                ],
                {lamVal, penList}
            ];
            formatTOOutput["JointPenalized", targetFidOpt, effCrit, discOpt, rescaleOpt,
                bestFeasRes, If[bestFeasRes =!= None, bestFeasDur, None], bestOverall,
                {tMin, If[bestFeasRes =!= None, bestFeasDur, t0]}, hist,
                totAtt, totObj, totGrad, totExp, AbsoluteTime[] - tStart, valEnsOpt]
    ]
];

TimeOptimalControlOptimize[
    drift:Except[_Rule | _OptimalControlProblem | _OptimalControlEnsemble | _CooperativeObjective | _List],
    controls:Except[_Rule], target_, timing:Except[_Rule], opts___Rule] :=
    TimeOptimalControlOptimize[
        OptimalControlProblem[drift, controls, target, timing,
            Sequence @@ FilterRules[FilterRules[{opts}, Except[Method]], Options[OptimalControlProblem]]],
        opts
    ];

TimeOptimalPulseCurve[input_, durations_List, opts___Rule] := Module[{tStart = AbsoluteTime[], norm, kind, obj, ownOpts, mergedOpts, grapeOpts,
     targetFidOpt, critOpt, discOpt, rescaleOpt, attOpt, valEnsOpt, seedOpt,
     defaultNt, dt0, effCrit, isSliceListQ, discSpec, attSpec, kBase, kBound,
     dursDesc, dursAsc, prevRes = None, recRev = {}, recAsc,
     totAtt = 0, totObj = 0, totGrad = 0, totExp = 0,
     bestFeasRes = None, bestFeasDur = None, dVal, dtCand, attRes, summ, valRes},

    If[durations === {} || !VectorQ[durations, NumericQ] || !AllTrue[durations, TrueQ[Im[N[#]] == 0 && Re[#] > 0.] &],
        Message[TimeOptimalPulseCurve::baddurs, durations]; Return[$Failed]
    ];
    norm = normalizeTOInput[input];
    If[norm === $Failed,
        Message[TimeOptimalControlOptimize::badprob, input]; Return[$Failed]
    ];
    {kind, obj, ownOpts} = norm;
    mergedOpts = Join[{opts}, ownOpts, Options[TimeOptimalPulseCurve], Options[GRAPEOptimize]];
    grapeOpts = FilterRules[mergedOpts, Except[Alternatives @@ $TOOnlyOptionNames]];

    {targetFidOpt, critOpt, discOpt, rescaleOpt, attOpt, valEnsOpt, seedOpt} =
        {TargetFidelity, FidelityCriterion, DiscretizationStrategy, WaveformRescaling,
         AttemptsPerDuration, ValidationEnsemble, RandomSeed} /. mergedOpts;

    If[!MemberQ[{"NormalizedTime", "AreaPreserving", "None"}, rescaleOpt],
        Message[TimeOptimalControlOptimize::badopt, WaveformRescaling -> rescaleOpt,
            "expected \"NormalizedTime\", \"AreaPreserving\", or \"None\""];
        Return[$Failed]
    ];
    effCrit = resolveEffectiveCriterion[kind, critOpt, mergedOpts];
    If[effCrit === $Failed,
        Message[TimeOptimalControlOptimize::badopt, FidelityCriterion -> critOpt,
            "expected Automatic, \"Fidelity\", \"Mean\", \"SoftMin\", or \"Minimum\""];
        Return[$Failed]
    ];

    defaultNt = TimeSteps /. mergedOpts /. Options[GRAPEOptimize];
    dt0 = extractBaselineTiming[kind, obj, defaultNt];
    If[dt0 === $Failed, Return[$Failed]];
    isSliceListQ = hasSliceListDriftQ[kind, obj];
    discSpec = validateDiscretization[discOpt, dt0, isSliceListQ];
    If[discSpec === $Failed, Return[$Failed]];
    attSpec = parseAttemptBudget[attOpt];
    If[attSpec === $Failed, Return[$Failed]];
    {kBase, kBound} = attSpec;

    (* Traverse largest-to-smallest so warm-start continuation starts from an easy regime *)
    dursAsc = Sort[DeleteDuplicates[N[Re[durations]]]];
    dursDesc = Reverse[dursAsc];

    Do[
        dtCand = buildDurationGrid[dVal, dt0, discSpec];
        If[dtCand === $Failed, Return[$Failed, Module]];
        attRes = runDurationAttempt[kind, obj, dtCand,
            If[rescaleOpt === "None", None, prevRes],
            grapeOpts, rescaleOpt,
            If[prevRes === None, Max[kBase, kBound], kBase],
            effCrit, targetFidOpt, seedOpt, prevRes === None];
        totAtt += attRes["AttemptsUsed"];
        totObj += attRes["ObjectiveEvaluations"];
        totGrad += attRes["GradientEvaluations"];
        totExp += attRes["MatrixExponentials"];
        summ = attRes["BestSummary"];
        If[AssociationQ[summ],
            If[rescaleOpt =!= "None" && (TrueQ[summ["FeasibleQ"]] || prevRes === None),
                prevRes = summ["Result"]
            ];
            valRes = If[valEnsOpt =!= None && valEnsOpt =!= False,
                evaluateValidationEnsemble[valEnsOpt, summ["Result"], effCrit, targetFidOpt],
                None
            ];
            AppendTo[recRev, Join[summ, <|
                "RequestedDuration" -> dVal,
                "AttemptsUsed" -> attRes["AttemptsUsed"],
                "ObjectiveEvaluations" -> attRes["ObjectiveEvaluations"],
                "GradientEvaluations" -> attRes["GradientEvaluations"],
                "MatrixExponentials" -> attRes["MatrixExponentials"],
                "ValidationResult" -> valRes
            |>]];
            If[TrueQ[summ["FeasibleQ"]],
                If[bestFeasDur === None || summ["Duration"] < bestFeasDur,
                    bestFeasDur = summ["Duration"];
                    bestFeasRes = summ["Result"]
                ]
            ],
            AppendTo[recRev, <|
                "Duration" -> dVal, "RequestedDuration" -> dVal, "NumSlices" -> Length[dtCand],
                "Fidelity" -> 0., "MeanFidelity" -> 0., "SoftMinFidelity" -> 0., "MinFidelity" -> 0.,
                "CriterionFidelity" -> 0., "MaxConstraintViolation" -> Infinity,
                "FeasibleQ" -> False, "AttemptsUsed" -> attRes["AttemptsUsed"],
                "ObjectiveEvaluations" -> attRes["ObjectiveEvaluations"],
                "GradientEvaluations" -> attRes["GradientEvaluations"],
                "MatrixExponentials" -> attRes["MatrixExponentials"],
                "ElapsedTime" -> 0., "Status" -> "Failed", "Result" -> None, "ValidationResult" -> None
            |>]
        ],
        {dVal, dursDesc}
    ];

    recAsc = Reverse[recRev];
    <|
        "Durations" -> (#["Duration"] & /@ recAsc),
        "Fidelities" -> (#["Fidelity"] & /@ recAsc),
        "MeanFidelities" -> (#["MeanFidelity"] & /@ recAsc),
        "SoftMinFidelities" -> (#["SoftMinFidelity"] & /@ recAsc),
        "MinFidelities" -> (#["MinFidelity"] & /@ recAsc),
        "CriterionFidelities" -> (#["CriterionFidelity"] & /@ recAsc),
        "FeasibleFlags" -> (#["FeasibleQ"] & /@ recAsc),
        "AttemptCounts" -> (#["AttemptsUsed"] & /@ recAsc),
        "TargetFidelity" -> N[targetFidOpt],
        "FidelityCriterion" -> effCrit,
        "DiscretizationStrategy" -> discOpt,
        "WaveformRescaling" -> rescaleOpt,
        "ShortestFeasibleDuration" -> bestFeasDur,
        "ShortestFeasibleResult" -> bestFeasRes,
        "Records" -> recAsc,
        "CurveData" -> Transpose[{#["Duration"] & /@ recAsc, #["CriterionFidelity"] & /@ recAsc}],
        "TotalAttempts" -> totAtt,
        "TotalObjectiveEvaluations" -> totObj,
        "TotalGradientEvaluations" -> totGrad,
        "TotalMatrixExponentials" -> totExp,
        "ElapsedTime" -> AbsoluteTime[] - tStart,
        "OptimalityDisclaimer" -> $OptimalityDisclaimer
    |>
];

TimeOptimalPulseCurve[
    drift:Except[_Rule | _OptimalControlProblem | _OptimalControlEnsemble | _CooperativeObjective | _List],
    controls:Except[_Rule], target_, durations_List, opts___Rule] :=
    If[durations === {} || !NumericQ[First[durations]],
        Message[TimeOptimalPulseCurve::baddurs, durations]; $Failed,
        TimeOptimalPulseCurve[
            OptimalControlProblem[drift, controls, target, First[durations],
                Sequence @@ FilterRules[{opts}, Options[OptimalControlProblem]]],
            durations, opts
        ]
    ];

(* TO5: Sub-slice Pontryagin Maximum Principle Switching Functions & Stationarity Diagnostic *)
classifySimpleConstraintBlocks[cSpec_Association, mCtrl_Integer] := Module[
    {hasBox = TrueQ[cSpec["HasBox"]], hasDisk = TrueQ[cSpec["HasDisk"]],
     disj = TrueQ[cSpec["DisjointDisks"]], groups, lims, lbRow, ubRow,
     diskChans, blocks = {}, overlapQ = False, cType},
    groups = Lookup[cSpec, "NormGroups", {}];
    lims = Lookup[cSpec, "NormLimits", {}];
    lbRow = If[hasBox, (cSpec["LowerBoundsMat"])[[1]], ConstantArray[-Infinity, mCtrl]];
    ubRow = If[hasBox, (cSpec["UpperBoundsMat"])[[1]], ConstantArray[Infinity, mCtrl]];
    If[!hasBox && !hasDisk,
        Return[<|"ValidQ" -> True, "Type" -> "Unconstrained", "Blocks" -> {<|"Kind" -> "Unconstrained", "Channels" -> Range[mCtrl]|>}|>]
    ];
    If[hasDisk && !disj,
        Return[<|"ValidQ" -> False, "Reason" -> "OverlappingDiskConstraints"|>]
    ];
    diskChans = Flatten[groups];
    If[hasBox && AnyTrue[diskChans, (lbRow[[#]] > -Infinity || ubRow[[#]] < Infinity) &],
        overlapQ = True
    ];
    If[overlapQ,
        Return[<|"ValidQ" -> False, "Reason" -> "OverlappingBoxAndDiskConstraints"|>]
    ];
    Do[
        AppendTo[blocks, <|"Kind" -> "Disk", "Channels" -> groups[[g]], "Limit" -> lims[[g]]|>],
        {g, Length[groups]}
    ];
    Do[
        If[!MemberQ[diskChans, k],
            If[NumericQ[lbRow[[k]]] && lbRow[[k]] > -Infinity && NumericQ[ubRow[[k]]] && ubRow[[k]] < Infinity,
                AppendTo[blocks, <|"Kind" -> "Box", "Channel" -> k, "Lower" -> lbRow[[k]], "Upper" -> ubRow[[k]]|>],
                If[lbRow[[k]] === -Infinity && ubRow[[k]] === Infinity,
                    AppendTo[blocks, <|"Kind" -> "Unconstrained", "Channels" -> {k}|>],
                    Return[<|"ValidQ" -> False, "Reason" -> "SemiInfiniteBoxBounds"|>, Module]
                ]
            ]
        ],
        {k, mCtrl}
    ];
    cType = Which[
        hasDisk && !hasBox, "Disk",
        hasBox && !hasDisk, "Box",
        True, "DisjointBoxAndDisk"
    ];
    <|"ValidQ" -> True, "Type" -> cType, "Blocks" -> blocks|>
];

subSliceProps[hRaw_, dtS_?NumericQ] := Module[{hJ, vals, vecs, vecsH, dH, dS},
    hJ = 0.5 * (hRaw + ConjugateTranspose[hRaw]);
    {vals, vecs} = Eigensystem[hJ];
    vals = Re[vals];
    vecs = Transpose[vecs];
    vecsH = ConjugateTranspose[vecs];
    dH = Exp[(-0.5 I dtS) * vals];
    dS = dH * dH;
    {ToPackedComplexMatrix[vecs . (dH * vecsH)], ToPackedComplexMatrix[vecs . (dS * vecsH)]}
];

evalMemberSubSliceSwitching[mem_Association, uMat_, dtVec_List, nSub_Integer, optPhases_] := Module[
    {nt = Length[dtVec], mCtrl = Dimensions[uMat][[2]], dim = mem["Dimension"],
     driftKind = Lookup[mem, "DriftKind", "Static"], h0 = mem["Drift"], hCtrls = mem["Controls"],
     dEval = Lookup[mem, "DriftEvaluator", None], obj = mem["Objective"],
     t0 = Lookup[mem, "InitialTimePoint", 0.], nTotal = Length[dtVec] * nSub,
     idMat, dtFine, tMidFine, tEndSlice, hJ, pair, uHalf, uSub, uFwd, uFinal,
     cotRes, fidVal, cDag, pfSpec, xMat, uMid, pMid, j, q, idx, sFine},

    idMat = IdentityMatrix[dim, WorkingPrecision -> MachinePrecision];
    dtFine = ConstantArray[0., nTotal];
    tMidFine = ConstantArray[0., nTotal];
    uHalf = Table[Null, {nTotal}];
    uSub = Table[Null, {nTotal}];

    tEndSlice = t0;
    idx = 0;
    Do[
        With[{dtS = dtVec[[j]] / nSub},
            If[driftKind === "Static",
                hJ = h0 + Sum[uMat[[j, k]] * hCtrls[[k]], {k, mCtrl}];
                pair = subSliceProps[hJ, dtS];
                Do[
                    idx++;
                    dtFine[[idx]] = dtS;
                    tMidFine[[idx]] = tEndSlice + (q - 0.5) * dtS;
                    uHalf[[idx]] = pair[[1]];
                    uSub[[idx]] = pair[[2]],
                    {q, nSub}
                ],
                Do[
                    idx++;
                    dtFine[[idx]] = dtS;
                    tMidFine[[idx]] = tEndSlice + (q - 0.5) * dtS;
                    hJ = If[driftKind === "SliceList", h0[[j]], dEval[tMidFine[[idx]]]] +
                        Sum[uMat[[j, k]] * hCtrls[[k]], {k, mCtrl}];
                    pair = subSliceProps[hJ, dtS];
                    uHalf[[idx]] = pair[[1]];
                    uSub[[idx]] = pair[[2]],
                    {q, nSub}
                ]
            ];
            tEndSlice += dtVec[[j]]
        ],
        {j, nt}
    ];

    uFwd = Table[Null, {nTotal + 1}];
    uFwd[[1]] = ToPackedComplexMatrix[idMat];
    Do[uFwd[[l + 1]] = uSub[[l]] . uFwd[[l]], {l, nTotal}];
    uFinal = uFwd[[nTotal + 1]];

    pfSpec = Lookup[obj, "PhaseFreedomSpec", None];
    cotRes = If[AssociationQ[pfSpec] && TrueQ[pfSpec["ParametricQ"]] && AssociationQ[optPhases],
        (pfSpec["EvaluateWithPhases"])[uFinal, optPhases["Vector"]],
        (obj["EvaluateCotangent"])[uFinal]
    ];
    fidVal = N[cotRes[[1]]];
    cDag = ToPackedComplexMatrix[cotRes[[2]]];

    xMat = cDag;
    sFine = ConstantArray[0., {nTotal, mCtrl}];
    Do[
        pMid = xMat . uHalf[[l]];
        uMid = uHalf[[l]] . uFwd[[l]];
        With[{bMid = uMid . pMid},
            Do[
                sFine[[l, k]] = Im[Tr[bMid . hCtrls[[k]]]],
                {k, mCtrl}
            ]
        ];
        xMat = xMat . uSub[[l]],
        {l, nTotal, 1, -1}
    ];

    <|
        "Fidelity" -> fidVal,
        "DtFine" -> Developer`ToPackedArray[dtFine],
        "TimeMidpoints" -> Developer`ToPackedArray[tMidFine],
        "Switching" -> Developer`ToPackedArray[sFine]
    |>
];

evalPointwiseResidual[uRow_List, sRow_List, blocks_List, uRefVec_List] := Module[{res = 0., denom = 0.},
    Do[
        Switch[blk["Kind"],
            "Box",
                With[{k = blk["Channel"], lo = blk["Lower"], hi = blk["Upper"]},
                    With[{sk = sRow[[k]], uk = uRow[[k]]},
                        res += Max[0., Max[hi * sk, lo * sk] - uk * sk];
                        denom += Max[Abs[lo], Abs[hi]] * Abs[sk]
                    ]
                ],
            "Disk",
                With[{ch = blk["Channels"], w = blk["Limit"]},
                    With[{sg = sRow[[ch]], ug = uRow[[ch]]},
                        res += Max[0., w * Norm[sg] - ug . sg];
                        denom += w * Norm[sg]
                    ]
                ],
            "Unconstrained",
                Do[
                    With[{ur = uRefVec[[k]], sk = sRow[[k]]},
                        res += ur * Abs[sk];
                        denom += ur * Abs[sk]
                    ],
                    {k, blk["Channels"]}
                ]
        ],
        {blk, blocks}
    ];
    {res, denom}
];

VerifyOptimalControlStationarity[input_, opts___Rule] := Module[
    {nSub, tol, assoc, kind, comp, members, weights, aggKind, beta, reasons = {},
     uMat, dtVec, nt, mCtrl, cSpec, blockSpec, optPhases,
     memOut, fids, fMin, expW, effW, sFineAgg, tFine, dtFine,
     sliceGradQuad, evalCheck, sliceGradExact, quadErr,
     blocks, uRefVec, uFine, rFine, dFine, rInt, dInt, relInt,
     rDisc = 0., dDisc = 0., sBarSlice, relDisc, maxPw, meanPw,
     bndCount = 0, bndDenom = 0, zeroCross, satQ},

    If[!MatchQ[input, _OptimalControlResult],
        Message[VerifyOptimalControlStationarity::badinput, input]; Return[$Failed]
    ];
    {nSub, tol} = {StationaritySubdivisions, Tolerance} /. {opts} /. Options[VerifyOptimalControlStationarity];
    If[!IntegerQ[nSub] || nSub < 1 || !NumericQ[tol] || !TrueQ[tol > 0.],
        Message[VerifyOptimalControlStationarity::badinput, {opts}]; Return[$Failed]
    ];

    assoc = input[[1]];
    kind = Lookup[assoc, "Kind", "Single"];
    comp = If[kind === "Ensemble", Lookup[assoc, "CompiledEnsemble", None], Lookup[assoc, "CompiledProblem", None]];
    If[!AssociationQ[comp], Message[VerifyOptimalControlStationarity::badinput, input]; Return[$Failed]];

    If[assoc["Space"] =!= "Hilbert",
        AppendTo[reasons, "UnsupportedSpace: continuous PMP switching residual only applies to Hilbert space (not Liouville)."]
    ];
    If[!TrueQ[comp["ControlBasis", "IsIdentityQ"]],
        AppendTo[reasons, "NonIdentityControlParameterization: ControlParameterization/ControlWindow must be identity PiecewiseConstant without windowing."]
    ];
    If[TrueQ[comp["HasSlewBounds"]],
        AppendTo[reasons, "CoupledSlewRateConstraint: SlewRateBound couples adjacent slices outside simple pointwise PMP."]
    ];
    members = If[kind === "Ensemble", comp["Members"], {comp}];
    If[AnyTrue[members, Lookup[#, "RunningCostSpec", None] =!= None &],
        AppendTo[reasons, "NonZeroRunningCost: intermediate-time RunningCost alters continuous costate dynamics."]
    ];
    If[Lookup[comp, "PowerPenalty", 0.] > 0. || Lookup[comp, "SmoothnessPenalty", 0.] > 0. ||
       Lookup[comp, "L1Penalty", 0.] > 0. || Lookup[comp, "TVPenalty", 0.] > 0.,
        AppendTo[reasons, "NonZeroControlRegularizationPenalty: PowerPenalty, SmoothnessPenalty, L1Penalty, or TVPenalty is active."]
    ];

    uMat = assoc["Controls"];
    dtVec = assoc["TimeSteps"];
    nt = Length[dtVec];
    mCtrl = Dimensions[uMat][[2]];
    cSpec = Lookup[comp, "ConstraintSpec", <|"HasBox" -> False, "HasDisk" -> False, "DisjointDisks" -> True|>];
    blockSpec = classifySimpleConstraintBlocks[cSpec, mCtrl];
    If[!TrueQ[blockSpec["ValidQ"]],
        AppendTo[reasons, "UnsupportedConstraintCombination: " <> blockSpec["Reason"]]
    ];

    If[reasons =!= {},
        Return[<|
            "ApplicableQ" -> False,
            "InapplicabilityReasons" -> reasons,
            "ConstraintType" -> Lookup[blockSpec, "Type", Lookup[blockSpec, "Reason", "Inapplicable"]],
            "QuadratureGradientError" -> None,
            "DiscreteResidual" -> None,
            "IntegratedResidual" -> None,
            "RelativeDiscreteResidual" -> None,
            "RelativeIntegratedResidual" -> None,
            "SwitchingScale" -> None,
            "MaxPointwiseResidual" -> None,
            "MeanPointwiseResidual" -> None,
            "IntraSliceResidualGap" -> None,
            "BoundaryFraction" -> None,
            "SwitchingZeroCrossings" -> None,
            "StationaritySatisfiedQ" -> False,
            "TimeGrid" -> {},
            "ControlTimeSeries" -> {},
            "SwitchingFunctions" -> {},
            "PointwiseResiduals" -> {},
            "Disclaimer" -> $PMPDisclaimer
        |>]
    ];

    weights = If[kind === "Ensemble", comp["Weights"], {1.}];
    aggKind = If[kind === "Ensemble", comp["Aggregation"], "Mean"];
    beta = If[kind === "Ensemble", comp["SoftMinBeta"], 0.];
    optPhases = Lookup[assoc, "OptimalPhases", None];

    memOut = Map[evalMemberSubSliceSwitching[#, uMat, dtVec, nSub, optPhases] &, members];
    fids = #["Fidelity"] & /@ memOut;
    effW = If[kind === "Ensemble" && aggKind === "SoftMin",
        fMin = Min[fids];
        expW = weights * Exp[-beta * (fids - fMin)];
        expW / Total[expW],
        weights
    ];
    tFine = memOut[[1, "TimeMidpoints"]];
    dtFine = memOut[[1, "DtFine"]];
    sFineAgg = Developer`ToPackedArray[Sum[effW[[m]] * memOut[[m, "Switching"]], {m, Length[members]}]];

    (* Quadrature consistency check against GRAPE FidelitySliceGradient *)
    sliceGradQuad = Developer`ToPackedArray[MapThread[#1 * Total[#2] &, {dtVec / nSub, Partition[sFineAgg, nSub]}]];
    evalCheck = If[kind === "Ensemble",
        ($OCSolverRegistry["EvaluateEnsemble"])[comp,
            If[TrueQ[comp["HasDurationOptimization"]] || TrueQ[comp["HasPhaseFreedom"]], assoc["FullParameters"], assoc["Parameters"]]],
        EvaluateOptimalControl[comp,
            If[TrueQ[comp["HasDurationOptimization"]] || TrueQ[comp["HasPhaseFreedom"]], assoc["FullParameters"], assoc["Parameters"]]]
    ];
    sliceGradExact = evalCheck["FidelitySliceGradient"];
    quadErr = Max[Abs[sliceGradQuad - sliceGradExact]] / Max[1., Max[Abs[sliceGradExact]]];

    blocks = blockSpec["Blocks"];
    uRefVec = Map[Max[Max[Abs[uMat[[All, #]]]], 1.] &, Range[mCtrl]];
    uFine = Flatten[ConstantArray[#, nSub] & /@ uMat, 1];
    {rFine, dFine} = Transpose[MapThread[evalPointwiseResidual[#1, #2, blocks, uRefVec] &, {uFine, sFineAgg}]];
    rInt = dtFine . rFine;
    dInt = dtFine . dFine;
    relInt = If[dInt > 10^-7, rInt / dInt, rInt];
    maxPw = Max[rFine];
    meanPw = rInt / Max[Total[dtFine], 10^-16];

    sBarSlice = sliceGradExact / dtVec;
    Do[
        With[{pair = evalPointwiseResidual[uMat[[j]], sBarSlice[[j]], blocks, uRefVec]},
            rDisc += dtVec[[j]] * pair[[1]];
            dDisc += dtVec[[j]] * pair[[2]]
        ],
        {j, nt}
    ];
    relDisc = If[dDisc > 10^-7, rDisc / dDisc, rDisc];

    Do[
        Switch[blk["Kind"],
            "Box",
                bndDenom += nt;
                With[{lo = blk["Lower"], hi = blk["Upper"], col = uMat[[All, blk["Channel"]]]},
                    With[{tolB = 10^-3 * Max[hi - lo, 1.]},
                        bndCount += Count[col, v_ /; (v >= hi - tolB || v <= lo + tolB)]
                    ]
                ],
            "Disk",
                bndDenom += nt;
                With[{w = blk["Limit"], sub = uMat[[All, blk["Channels"]]]},
                    With[{nrms = Norm /@ sub},
                        bndCount += Count[nrms, v_ /; v >= w * (1. - 10^-3)]
                    ]
                ]
        ],
        {blk, blocks}
    ];

    zeroCross = Table[
        With[{nz = Select[sFineAgg[[All, k]], Abs[#] > 10^-12 * Max[1., Max[Abs[sFineAgg[[All, k]]]]] &]},
            If[Length[nz] <= 1, 0, Count[Most[nz] * Rest[nz], v_ /; v < 0.]]
        ],
        {k, mCtrl}
    ];

    satQ = With[{qTol = If[AnyTrue[members, Lookup[#, "DriftKind", "Static"] === "Continuous" &], 10^-3, 10^-5]},
        If[blockSpec["Type"] === "Unconstrained",
            TrueQ[quadErr <= qTol && relDisc <= tol],
            TrueQ[quadErr <= qTol && relDisc <= tol && relInt <= Max[tol, 5. * tol]]
        ]
    ];

    <|
        "ApplicableQ" -> True,
        "InapplicabilityReasons" -> {},
        "ConstraintType" -> blockSpec["Type"],
        "QuadratureGradientError" -> N[quadErr],
        "DiscreteResidual" -> N[rDisc],
        "IntegratedResidual" -> N[rInt],
        "RelativeDiscreteResidual" -> N[relDisc],
        "RelativeIntegratedResidual" -> N[relInt],
        "SwitchingScale" -> N[dInt],
        "MaxPointwiseResidual" -> N[maxPw],
        "MeanPointwiseResidual" -> N[meanPw],
        "IntraSliceResidualGap" -> N[rInt - rDisc],
        "BoundaryFraction" -> If[bndDenom > 0, N[bndCount / bndDenom], 0.],
        "SwitchingZeroCrossings" -> zeroCross,
        "StationaritySatisfiedQ" -> satQ,
        "TimeGrid" -> tFine,
        "ControlTimeSeries" -> Developer`ToPackedArray[uFine],
        "SwitchingFunctions" -> sFineAgg,
        "PointwiseResiduals" -> Developer`ToPackedArray[rFine],
        "Disclaimer" -> $PMPDisclaimer
    |>
];

End[];
EndPackage[];
