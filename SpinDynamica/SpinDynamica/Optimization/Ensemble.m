(* ::Package:: *)

(* SpinDynamica Optimal Control — Robust Ensemble Optimal Control, Affine Hoisting & Parallel Evaluation *)

BeginPackage["OptimalControl`Ensemble`", {
    "Global`",
    "SpinDynamica`",
    "General`Matrices`",
    "QM`Hilbert`",
    "QM`SpinOperators`",
    "QM`Liouville`",
    "Events`Events`",
    "Ensemble`EnsembleAverage`",
    "OptimalControl`Kernels`",
    "OptimalControl`Objectives`",
    "OptimalControl`GRAPE`"
}];

OptimalControlEnsemble::usage = "OptimalControlEnsemble[drift, controls, target, timing, opts] or OptimalControlEnsemble[problems, opts] defines a robust or cooperative ensemble optimal-control problem across varying Hamiltonians, RF inhomogeneity, orientational grids, or heterogeneous spin systems.";
EnsembleParameters::usage = "EnsembleParameters is an option for OptimalControlEnsemble specifying symbolic parameter sweeps such as {delta -> {-100, 0, 100}} or weighted grids.";
RFInhomogeneity::usage = "RFInhomogeneity is an option for OptimalControlEnsemble specifying relative RF scaling factors (e.g. {0.8, 1.0, 1.2} or weighted rules) applied to all control channels.";
EnsembleAggregation::usage = "EnsembleAggregation is an option for OptimalControlEnsemble selecting \"Mean\" (weighted average) or {\"SoftMin\", beta} (smooth worst-case maximization).";
EvaluateEnsembleFidelity::usage = "EvaluateEnsembleFidelity[ensembleOrResult, opts] or EvaluateEnsembleFidelity[ensemble, controls] evaluates member-by-member fidelities, weights, mean, minimum, and maximum across the ensemble.";

OptimalControlEnsemble::inval = "Invalid ensemble specification `1`.";
OptimalControlEnsemble::incompat = "All ensemble members must share identical NumControls, TimeSteps, and ControlBasis.";
OptimalControlEnsemble::badagg = "Invalid EnsembleAggregation `1`; expected \"Mean\" or {\"SoftMin\", beta} with beta > 0.";
OptimalControlEnsemble::badschedule = "Invalid SoftMinSchedule `1`; expected None, Automatic, or a non-empty list of non-negative numbers.";

Begin["`Private`"];

ensembleOnlyOptions = {
    EnsembleParameters -> {},
    RFInhomogeneity -> None,
    NormalizeWeights -> True,
    EnsembleAggregation -> "Mean",
    Parallel -> False
};

Options[OptimalControlEnsemble] = Join[Options[OptimalControlProblem], ensembleOnlyOptions];

(* Parse a single axis specification into {vars, valueTuples, weights} *)
parseAxisSpec[Rule[varSpec_, listSpec_List]] := Module[{vals, ws, cleanItem},
    If[listSpec === {}, Return[$Failed]];
    cleanItem = Map[
        Function[entry,
            Which[
                MatchQ[entry, Rule[_, _?NumericQ]], {entry[[1]], N[entry[[2]]]},
                ListQ[varSpec] && MatchQ[entry, {_List, _?NumericQ}] && Length[entry[[1]]] === Length[varSpec],
                    {entry[[1]], N[entry[[2]]]},
                !ListQ[varSpec] && MatchQ[entry, {_?NumericQ, _?NumericQ}],
                    {entry[[1]], N[entry[[2]]]},
                !ListQ[varSpec] && MatchQ[entry, {_List, _?NumericQ}],
                    {entry[[1]], N[entry[[2]]]},
                True,
                    {entry, 1.}
            ]
        ],
        listSpec
    ];
    vals = cleanItem[[All, 1]];
    ws = cleanItem[[All, 2]];
    If[!AllTrue[ws, NumericQ] || Total[ws] <= 0., Return[$Failed]];
    {varSpec, vals, ws}
];
parseAxisSpec[_] := $Failed;

parseRFSpec[None | False | Automatic] := {{1.}, {1.}, False};
parseRFSpec[spec_List] := Module[{parsed, rfScaleSym},
    If[spec === {}, Return[$Failed]];
    parsed = parseAxisSpec[Rule[rfScaleSym, spec]];
    If[parsed === $Failed || !AllTrue[parsed[[2]], NumericQ], Return[$Failed]];
    {N[parsed[[2]]], N[parsed[[3]]], True}
];
parseRFSpec[_] := $Failed;

buildMemberGrid[epSpec_, rfSpec_, normW_] := Module[
    {epList, axes, rfVals, rfWs, hasRF, allAxesVars, allAxesVals, allAxesWs,
     nAxes, tuples, rawWeights, rulesList, rfList, paramPoints, gridShape},
    epList = Which[
        epSpec === {} || epSpec === None || epSpec === False, {},
        MatchQ[epSpec, _Rule], {epSpec},
        ListQ[epSpec] && AllTrue[epSpec, MatchQ[#, _Rule] &], epSpec,
        True, $Failed
    ];
    If[epList === $Failed, Return[$Failed]];
    axes = parseAxisSpec /@ epList;
    If[!FreeQ[axes, $Failed], Return[$Failed]];
    {rfVals, rfWs, hasRF} = parseRFSpec[rfSpec];
    If[rfVals === $Failed, Return[$Failed]];

    If[axes === {},
        rawWeights = If[TrueQ[normW], rfWs / Total[rfWs], rfWs];
        rulesList = ConstantArray[{}, Length[rfVals]];
        rfList = rfVals;
        paramPoints = If[hasRF, List /@ rfVals, None];
        gridShape = If[hasRF, {Length[rfVals]}, {1}];
        Return[<|
            "Rules" -> rulesList,
            "RFScales" -> rfList,
            "Weights" -> Developer`ToPackedArray[N[rawWeights]],
            "ParameterPoints" -> paramPoints,
            "GridShape" -> gridShape,
            "AxisVars" -> {}
        |>]
    ];

    allAxesVars = axes[[All, 1]];
    allAxesVals = axes[[All, 2]];
    allAxesWs = axes[[All, 3]];
    If[hasRF,
        AppendTo[allAxesVals, rfVals];
        AppendTo[allAxesWs, rfWs]
    ];
    nAxes = Length[allAxesVals];
    gridShape = Length /@ allAxesVals;
    tuples = Tuples[allAxesVals];
    rawWeights = Apply[Times, Tuples[allAxesWs], {1}];
    If[TrueQ[normW], rawWeights = rawWeights / Total[rawWeights]];

    rulesList = Map[
        Function[tup,
            Flatten[Table[
                If[ListQ[allAxesVars[[i]]],
                    Thread[allAxesVars[[i]] -> tup[[i]]],
                    {allAxesVars[[i]] -> tup[[i]]}
                ],
                {i, Length[allAxesVars]}
            ]]
        ],
        tuples
    ];
    rfList = If[hasRF, tuples[[All, -1]], ConstantArray[1., Length[tuples]]];
    paramPoints = Map[
        Function[tup,
            Flatten[Table[tup[[i]], {i, nAxes}]]
        ],
        tuples
    ];
    If[Length[gridShape] === 1 && Length[First[paramPoints]] === 1,
        paramPoints = Flatten[paramPoints]
    ];
    <|
        "Rules" -> rulesList,
        "RFScales" -> Developer`ToPackedArray[N[rfList]],
        "Weights" -> Developer`ToPackedArray[N[rawWeights]],
        "ParameterPoints" -> paramPoints,
        "GridShape" -> gridShape,
        "AxisVars" -> DeleteDuplicates[Flatten[allAxesVars]]
    |>
];

(* Try affine symbolic hoisting on drift/control matrix expressions *)
tryAffineDecomp[expr_, scalarSyms_List, sampleRules_List, space_String, basis_, opbasis_, dim_Integer, unitFactor_] := Module[
    {symMat, rawMat, decomp, testVals, mTest, exactMat, diff, evalList},
    If[scalarSyms === {} || Length[sampleRules] < 2 || FreeQ[expr, Alternatives @@ scalarSyms],
        Return[$Failed]
    ];
    If[MatchQ[expr, _NullOperator | _NullSuperoperator] || (NumericQ[expr] && TrueQ[expr == 0]),
        Return[$Failed]
    ];
    symMat = Quiet[Check[
        If[space === "Hilbert",
            If[OperatorQ[expr], MatrixRepresentation[expr, basis], expr],
            Which[
                SuperoperatorQ[expr], SuperoperatorMatrixRepresentation[expr, opbasis, basis],
                OperatorQ[expr], SuperoperatorMatrixRepresentation[CommutationSuperoperator[-I expr], opbasis, basis],
                True, expr
            ]
        ],
        $Failed
    ]];
    If[symMat === $Failed || !MatrixQ[symMat] || Dimensions[symMat] =!= {dim, dim},
        Return[$Failed]
    ];
    rawMat = Normal[symMat];
    decomp = Quiet[Check[General`Matrices`Private`decomposeSymbolicMatrix[rawMat, scalarSyms], $Failed]];
    If[!AssociationQ[decomp], Return[$Failed]];

    (* Verify exactness on multiple sample points (guards against multi-axis nonlinear expressions) *)
    Do[
        testVals = N[scalarSyms /. sampleRules[[idx]]];
        If[!VectorQ[testVals, NumericQ], Return[$Failed, Module]];
        mTest = Normal[General`Matrices`Private`evaluateAffineMatrix[decomp, testVals]];
        exactMat = N[rawMat /. sampleRules[[idx]]];
        If[!MatrixQ[mTest, NumericQ] || !MatrixQ[exactMat, NumericQ], Return[$Failed, Module]];
        diff = Max[Abs[mTest - exactMat]];
        If[!TrueQ[diff <= 10^-11 * Max[1., Max[Abs[exactMat]]]], Return[$Failed, Module]],
        {idx, DeleteDuplicates[{1, Ceiling[Length[sampleRules] / 2], Length[sampleRules]}]}
    ];

    evalList = Map[
        Function[rule,
            ToPackedComplexMatrix[unitFactor * Normal[General`Matrices`Private`evaluateAffineMatrix[decomp, N[scalarSyms /. rule]]]]
        ],
        sampleRules
    ];
    If[space === "Hilbert",
        evalList = HermitianSymmetrize /@ evalList
    ];
    evalList
];

(* Extract lightweight worker payload so subkernels only need OptimalControl`Kernels` *)
packMemberForKernel[comp_Association] := <|
    "Space" -> comp["Space"],
    "Method" -> comp["Method"],
    "Dimension" -> comp["Dimension"],
    "Drift" -> comp["Drift"],
    "FineTimeSteps" -> Lookup[comp, "FineTimeSteps", comp["TimeSteps"]],
    "NumStages" -> Lookup[comp, "NumStages", 1],
    "Controls" -> comp["Controls"],
    "ControlRows" -> comp["ControlRows"],
    "ControlFlatT" -> comp["ControlFlatT"],
    "TimeSteps" -> comp["TimeSteps"],
    "RunningCostSpec" -> Lookup[comp, "RunningCostSpec", None],
    "Objective" -> <|
        "EvaluateCotangent" -> comp["Objective", "EvaluateCotangent"],
        "EvaluateCotangentDirectional" -> Lookup[comp["Objective"], "EvaluateCotangentDirectional", None],
        "ExactHessianSupportedQ" -> TrueQ[Lookup[comp["Objective"], "ExactHessianSupportedQ", False]],
        "UnsupportedSecondOrderReason" -> Lookup[comp["Objective"], "UnsupportedSecondOrderReason", None],
        "EvaluateWithPhases" -> If[AssociationQ[Lookup[comp["Objective"], "PhaseFreedomSpec", None]],
            Lookup[comp["Objective", "PhaseFreedomSpec"], "EvaluateWithPhases", None], None],
        "PureKetFastPath" -> comp["Objective", "PureKetFastPath"],
        "InitialMatrix" -> Lookup[comp["Objective"], "InitialMatrix", None],
        "TargetMatrix" -> Lookup[comp["Objective"], "TargetMatrix", None],
        "Weights" -> Lookup[comp["Objective"], "Weights", None],
        "Modes" -> Lookup[comp["Objective"], "Modes", None],
        "Scales" -> Lookup[comp["Objective"], "Scales", None]
    |>
|>;

evalMemberKernel[mem_Association, uMat_, stepOverride_: None, phiVec_: None] := With[
    {h0 = If[AssociationQ[stepOverride] && Lookup[stepOverride, "Drift", None] =!= None,
         stepOverride["Drift"], mem["Drift"]],
     dtFine = If[AssociationQ[stepOverride], stepOverride["FineTimeSteps"], Lookup[mem, "FineTimeSteps", mem["TimeSteps"]]],
     nStages = Lookup[mem, "NumStages", 1],
     rcSpec = Lookup[mem, "RunningCostSpec", None],
     timingSpec = If[AssociationQ[stepOverride], stepOverride["TimingSpec"], None],
     cotFn = With[{ewp = Lookup[mem["Objective"], "EvaluateWithPhases",
             With[{pfs = Lookup[mem["Objective"], "PhaseFreedomSpec", None]},
                 If[AssociationQ[pfs], Lookup[pfs, "EvaluateWithPhases", None], None]
             ]]},
         If[phiVec =!= None && ewp =!= None,
             With[{fn = ewp, pv = phiVec},
                 Function[uFinal, With[{r = fn[uFinal, pv]}, {r[[1]], r[[2]], r[[3]]}]]
             ],
             mem["Objective", "EvaluateCotangent"]
         ]
     ]},
    Switch[mem["Method"],
        "Spectral",
            EvaluateHilbertSpectral[
                h0, mem["ControlRows"], mem["ControlFlatT"],
                dtFine, uMat, cotFn, False, nStages, rcSpec, timingSpec
            ],
        "BlockMatrixExp",
            EvaluateHilbertBlockExp[
                h0, mem["Controls"], dtFine, uMat,
                cotFn, False, nStages, rcSpec, timingSpec
            ],
        "FirstOrder",
            EvaluateHilbertFirstOrder[
                h0, mem["Controls"], dtFine, uMat,
                cotFn, False, nStages, rcSpec, timingSpec
            ],
        "StateVector",
            With[{pk = mem["Objective", "PureKetFastPath"]},
                EvaluateHilbertStateVector[
                    h0, mem["ControlRows"], mem["ControlFlatT"],
                    dtFine, uMat, pk["PsiMat0"], pk["TarMat"], pk["Weights"], pk["PhaseInsensitive"], nStages, rcSpec, timingSpec
                ]
            ],
        "Liouville",
            With[{ls = mem["Objective"]},
                EvaluateLiouvilleExact[
                    h0, mem["Controls"], dtFine, uMat,
                    ls["InitialMatrix"], ls["TargetMatrix"], ls["Weights"], ls["Modes"], ls["Scales"], nStages, rcSpec, timingSpec
                ]
            ]
    ]
];

evalMemberObjectiveOnlyKernel[mem_Association, uMat_, stepOverride_: None, phiVec_: None] := With[
    {h0 = If[AssociationQ[stepOverride] && Lookup[stepOverride, "Drift", None] =!= None,
         stepOverride["Drift"], mem["Drift"]],
     dtFine = If[AssociationQ[stepOverride], stepOverride["FineTimeSteps"], Lookup[mem, "FineTimeSteps", mem["TimeSteps"]]],
     nStages = Lookup[mem, "NumStages", 1],
     rcSpec = Lookup[mem, "RunningCostSpec", None],
     cotFn = With[{ewp = Lookup[mem["Objective"], "EvaluateWithPhases",
             With[{pfs = Lookup[mem["Objective"], "PhaseFreedomSpec", None]},
                 If[AssociationQ[pfs], Lookup[pfs, "EvaluateWithPhases", None], None]
             ]]},
         If[phiVec =!= None && ewp =!= None,
             With[{fn = ewp, pv = phiVec},
                 Function[uFinal, With[{r = fn[uFinal, pv]}, {r[[1]], r[[2]], r[[3]]}]]
             ],
             mem["Objective", "EvaluateCotangent"]
         ]
     ]},
    Switch[mem["Method"],
        "Spectral" | "BlockMatrixExp" | "FirstOrder",
            OptimalControl`Kernels`Private`EvaluateHilbertSpectralObjectiveOnly[
                h0, mem["ControlRows"], dtFine, uMat, cotFn, nStages, rcSpec
            ],
        "StateVector",
            With[{pk = mem["Objective", "PureKetFastPath"]},
                OptimalControl`Kernels`Private`EvaluateHilbertStateVectorObjectiveOnly[
                    h0, mem["ControlRows"], dtFine, uMat, pk["PsiMat0"], pk["TarMat"], pk["Weights"], pk["PhaseInsensitive"], nStages, rcSpec
                ]
            ],
        "Liouville",
            With[{ls = mem["Objective"]},
                OptimalControl`Kernels`Private`EvaluateLiouvilleObjectiveOnly[
                    h0, mem["Controls"], dtFine, uMat,
                    ls["InitialMatrix"], ls["TargetMatrix"], ls["Weights"], ls["Modes"], ls["Scales"], nStages, rcSpec
                ]
            ]
    ]
];

evalMemberKernelWithCache[mem_Association, uMat_] := With[
    {h0 = mem["Drift"],
     dtFine = Lookup[mem, "FineTimeSteps", mem["TimeSteps"]],
     nStages = Lookup[mem, "NumStages", 1],
     rcSpec = Lookup[mem, "RunningCostSpec", None],
     cotFn = mem["Objective", "EvaluateCotangent"]},
    Switch[mem["Method"],
        "Spectral",
            EvaluateHilbertSpectral[
                h0, mem["ControlRows"], mem["ControlFlatT"],
                dtFine, uMat, cotFn, False, nStages, rcSpec, None, True
            ],
        "StateVector",
            With[{pk = mem["Objective", "PureKetFastPath"]},
                EvaluateHilbertStateVector[
                    h0, mem["ControlRows"], mem["ControlFlatT"],
                    dtFine, uMat, pk["PsiMat0"], pk["TarMat"], pk["Weights"], pk["PhaseInsensitive"], nStages, rcSpec, None, True
                ]
            ],
        "Liouville",
            With[{ls = mem["Objective"]},
                EvaluateLiouvilleExact[
                    h0, mem["Controls"], dtFine, uMat,
                    ls["InitialMatrix"], ls["TargetMatrix"], ls["Weights"], ls["Modes"], ls["Scales"], nStages, rcSpec, None, True
                ]
            ],
        _,
            $Failed
    ]
];

evalMemberHVPKernel[mem_Association, mState_Association, vSliceMat_] := Switch[mem["Method"],
    "Spectral",
        OptimalControl`Kernels`Private`EvaluateHilbertSpectralHVP[
            mState, vSliceMat, mem["Objective", "EvaluateCotangentDirectional"]
        ],
    "StateVector",
        OptimalControl`Kernels`Private`EvaluateHilbertStateVectorHVP[mState, vSliceMat],
    "Liouville",
        OptimalControl`Kernels`Private`EvaluateLiouvilleHVP[mState, vSliceMat],
    _,
        $Failed
];

evalMemberHessianKernel[mem_Association, mState_Association] := Switch[mem["Method"],
    "Spectral",
        OptimalControl`Kernels`Private`EvaluateHilbertSpectralHessian[
            mState, mem["Controls"], mem["Objective", "EvaluateCotangentDirectional"]
        ],
    _,
        $Failed
];

(* Setup worker-resident partitions so only uMat moves during optimization *)
setupParallelWorkers[compEns_Association, parSpec_] := Module[
    {reqK, curK, kernels, kernelIds, nK, nM, members, chunks, workerChunks, evalFn, kNames, cacheId},
    reqK = If[IntegerQ[parSpec] && parSpec >= 2, Min[parSpec, Length[compEns["Members"]]],
        Max[2, Min[$ProcessorCount, Length[compEns["Members"]]]]];
    curK = Kernels[];
    If[Length[curK] < reqK,
        Quiet[LaunchKernels[reqK - Length[curK]]];
        curK = Kernels[]
    ];
    If[Length[curK] === 0, Return[compEns]];
    nK = Min[Length[curK], Length[compEns["Members"]]];
    kernels = Take[curK, nK];
    kernelIds = ParallelEvaluate[$KernelID, kernels];
    kNames = {
        OptimalControl`Kernels`EvaluateHilbertSpectral,
        OptimalControl`Kernels`Private`EvaluateHilbertSpectralObjectiveOnly,
        OptimalControl`Kernels`Private`EvaluateHilbertSpectralHVP,
        OptimalControl`Kernels`Private`EvaluateHilbertSpectralHessian,
        OptimalControl`Kernels`EvaluateHilbertBlockExp,
        OptimalControl`Kernels`EvaluateHilbertFirstOrder,
        OptimalControl`Kernels`EvaluateHilbertStateVector,
        OptimalControl`Kernels`Private`EvaluateHilbertStateVectorObjectiveOnly,
        OptimalControl`Kernels`Private`EvaluateHilbertStateVectorHVP,
        OptimalControl`Kernels`EvaluateLiouvilleExact,
        OptimalControl`Kernels`Private`EvaluateLiouvilleObjectiveOnly,
        OptimalControl`Kernels`Private`EvaluateLiouvilleHVP,
        OptimalControl`Kernels`Private`EnsureSpectralG2Cache,
        OptimalControl`Kernels`Private`buildG2Tensor,
        OptimalControl`Kernels`Private`contractG2,
        OptimalControl`Kernels`Private`$compiledG2TensorFn,
        OptimalControl`Kernels`Private`$compiledG2ContractFn,
        OptimalControl`Kernels`Private`ToPackedComplexTensor,
        OptimalControl`Kernels`HermitianSymmetrize,
        OptimalControl`Kernels`ToPackedComplexMatrix,
        OptimalControl`Kernels`ToPackedComplexVector,
        OptimalControl`Kernels`ToPackedRealMatrix,
        OptimalControl`Kernels`Private`prepareRunningRuntime,
        OptimalControl`Kernels`Private`evalRunningBoundary,
        OptimalControl`Kernels`Private`evalRunningTimeGradient,
        OptimalControl`Kernels`Private`foldTimingGradient,
        evalMemberKernel,
        evalMemberObjectiveOnlyKernel,
        evalMemberKernelWithCache,
        evalMemberHVPKernel,
        evalMemberHessianKernel
    };
    DistributeDefinitions @@ kNames;

    members = packMemberForKernel /@ compEns["Members"];
    nM = Length[members];
    chunks = Table[Range[i, nM, nK], {i, nK}];
    workerChunks = Map[members[[#]] &, chunks];
    cacheId = CreateUUID[];
    ParallelEvaluate[
        If[!AssociationQ[OptimalControl`Ensemble`Private`$WorkerMemberStore],
            OptimalControl`Ensemble`Private`$WorkerMemberStore = <||>
        ];
        If[!AssociationQ[OptimalControl`Ensemble`Private`$WorkerStateCache],
            OptimalControl`Ensemble`Private`$WorkerStateCache = <||>
        ],
        kernels
    ];
    Do[
        With[{cid = cacheId, chunk = workerChunks[[w]]},
            ParallelEvaluate[
                OptimalControl`Ensemble`Private`$WorkerMemberStore[cid] = chunk;
                Null,
                kernels[[w]]
            ]
        ],
        {w, nK}
    ];
    Join[compEns, <|
        "UseParallel" -> True,
        "WorkerKernels" -> kernels,
        "WorkerKernelIDs" -> kernelIds,
        "WorkerChunks" -> chunks,
        "WorkerCacheId" -> cacheId
    |>]
];

computeEnsembleParameterScales[memberComps_List, dtVec_List, basisMat_, mCtrl_Integer] := Module[
    {chanNorms, colScales},
    chanNorms = Table[
        Max[Mean[Map[Norm[Flatten[(#["Controls"])[[k]]]] / Sqrt[N[#["Dimension"]]] &, memberComps]], 10^-6],
        {k, mCtrl}
    ];
    colScales = Map[Max[Norm[dtVec * #], 10^-12] &, Transpose[basisMat]];
    <|
        "ChannelNorms" -> Developer`ToPackedArray[chanNorms],
        "ParameterScales" -> Developer`ToPackedArray[Flatten[Outer[Times, colScales, chanNorms]]]
    |>
];

(* Compile a parametric OptimalControlEnsemble *)
compileParametricEnsemble[ens:OptimalControlEnsemble[driftExpr_, controls_, targetExpr_, timing_, opts___Rule], extraOpts___Rule] := Module[
    {mergedOpts, epOpt, rfOpt, normWOpt, aggOpt, parOpt, basisOpt, opbasisOpt, spaceOpt, methodOpt,
     paramOpt, winOpt, freqUnit, ctrlUnit, defaultNt, t0Opt, driftIntOpt, runCostOpt, optDurOpt,
     durBndOpt, durPenOpt, fixedTotOpt, initOpt, powerPen, smoothPen, l1Pen, tvPen,
     ctrlBounds, ampBound, ampNormBound, slewBound,
     grid, nM, rulesList, rfList, weights, aggKind, aggBeta, basis, opbasis, ctrlList, hasSop, space, dim,
     dtVec, totalDur, nt, hamFactor, ctrlFactor, ctrlBasis, nParams, mCtrl,
     resolveUnit, toGen, resolveTim, resolveCSpec, resolveDur, parseDurPen, compileProx, compileDriftFn, hasSopDriftFn,
     cSpecRaw, durSpec, dpSpec, hasDur, cSpec, proxCtx, scalarSyms, staticScalarDriftQ, hoistedDrift, hoistedCtrls,
     targetStatic, compiledStaticObj, memberComps, baseComp, res, mIdx, hermErr,
     baseScales, ctrlScales, fullScales, lbOut, ubOut, hasBoundsOut, t0First},

    mergedOpts = Join[{extraOpts}, {opts}, Options[OptimalControlEnsemble]];
    {epOpt, rfOpt, normWOpt, aggOpt, parOpt, basisOpt, opbasisOpt, spaceOpt, methodOpt,
     paramOpt, winOpt, freqUnit, ctrlUnit, defaultNt, t0Opt, driftIntOpt, runCostOpt,
     optDurOpt, durBndOpt, durPenOpt, fixedTotOpt, initOpt, powerPen, smoothPen, l1Pen, tvPen,
     ctrlBounds, ampBound, ampNormBound, slewBound} =
        {EnsembleParameters, RFInhomogeneity, NormalizeWeights, EnsembleAggregation, Parallel,
         Basis, OperatorBasis, PropagateSpace, Method, ControlParameterization, ControlWindow,
         FrequencyUnits, ControlUnits, TimeSteps, InitialTimePoint, DriftIntegration, RunningCost,
         OptimizeDurations, DurationBounds, DurationPenalty, FixedTotalDuration, InitialControls,
         PowerPenalty, SmoothnessPenalty, L1Penalty, TVPenalty,
         ControlBounds, AmplitudeBound, AmplitudeNormBound, SlewRateBound} /. mergedOpts;

    Which[
        aggOpt === "Mean" || aggOpt === Automatic, aggKind = "Mean"; aggBeta = 0.,
        MatchQ[aggOpt, {"SoftMin", _?NumericQ}] && TrueQ[aggOpt[[2]] > 0],
            aggKind = "SoftMin"; aggBeta = N[aggOpt[[2]]],
        True, Message[OptimalControlEnsemble::badagg, aggOpt]; Return[$Failed]
    ];

    grid = buildMemberGrid[epOpt, rfOpt, normWOpt];
    If[grid === $Failed, Message[OptimalControlEnsemble::inval, ens]; Return[$Failed]];
    rulesList = grid["Rules"];
    rfList = grid["RFScales"];
    weights = grid["Weights"];
    nM = Length[weights];

    resolveTim = $OCSolverRegistry["ResolveTimingVector"];
    resolveUnit = $OCSolverRegistry["ResolveUnitFactor"];
    resolveCSpec = $OCSolverRegistry["ResolveConstraintSpec"];
    resolveDur = $OCSolverRegistry["ResolveDurationSpec"];
    parseDurPen = $OCSolverRegistry["ParseDurationPenalty"];
    compileProx = $OCSolverRegistry["CompileProxContext"];
    compileDriftFn = $OCSolverRegistry["CompileDrift"];
    hasSopDriftFn = $OCSolverRegistry["DriftContainsSuperoperatorQ"];
    toGen = $OCSolverRegistry["ToGenMatrix"];

    basis = If[basisOpt === Automatic, ZeemanBasis[], basisOpt];
    opbasis = If[opbasisOpt === Automatic, ZeemanKetBraOperatorBasis[Sorted -> "Canonical"], opbasisOpt];
    ctrlList = If[ListQ[controls] && !MatrixQ[controls], controls, {controls}];
    If[ctrlList === {}, Message[OptimalControlProblem::badctrl]; Return[$Failed]];

    dtVec = resolveTim[timing, defaultNt];
    If[dtVec === $Failed, Return[$Failed]];
    nt = Length[dtVec];
    totalDur = Total[dtVec];

    hamFactor = resolveUnit[freqUnit];
    ctrlFactor = If[ctrlUnit === Automatic, hamFactor, resolveUnit[ctrlUnit]];
    If[hamFactor === $Failed || ctrlFactor === $Failed, Return[$Failed]];

    t0First = With[{tCheck = t0Opt /. First[rulesList]},
        If[NumericQ[tCheck] && FreeQ[tCheck, Complex | DirectedInfinity | Indeterminate], N[Re[tCheck]], 0.]
    ];
    hasSop = hasSopDriftFn[driftExpr /. First[rulesList], t0First, totalDur] || AnyTrue[ctrlList /. First[rulesList], SuperoperatorQ];
    space = If[spaceOpt === Automatic,
        If[hasSop, "Liouville", "Hilbert"],
        If[MemberQ[{"Hilbert", "Liouville"}, spaceOpt], spaceOpt, $Failed]
    ];
    If[space === $Failed, Return[$Failed]];
    dim = If[space === "Hilbert", BasisDimension[basis], OperatorBasisDimension[opbasis]];
    If[!IntegerQ[dim] || dim < 1, Return[$Failed]];

    ctrlBasis = BuildControlBasis[paramOpt, dtVec, winOpt];
    If[ctrlBasis === $Failed, Return[$Failed]];
    nParams = ctrlBasis["NumParams"];
    mCtrl = Length[ctrlList];

    cSpecRaw = resolveCSpec[ctrlBounds, ampBound, ampNormBound, slewBound, dtVec, nParams, mCtrl];
    If[cSpecRaw === $Failed, Return[$Failed]];

    scalarSyms = Select[grid["AxisVars"], MatchQ[#, _Symbol] &];
    staticScalarDriftQ = FreeQ[driftExpr, _PeriodicFunction | _TimeDependentGeneratorFunction | _Function | _InterpolatingFunction] &&
        (!ListQ[driftExpr] || MatrixQ[driftExpr, NumericQ]);
    hoistedDrift = If[staticScalarDriftQ,
        tryAffineDecomp[driftExpr, scalarSyms, rulesList, space, basis, opbasis, dim, hamFactor],
        $Failed
    ];
    hoistedCtrls = Table[
        tryAffineDecomp[ctrlList[[k]], scalarSyms, rulesList, space, basis, opbasis, dim, ctrlFactor],
        {k, mCtrl}
    ];

    targetStatic = FreeQ[targetExpr, Alternatives @@ scalarSyms] || scalarSyms === {};
    compiledStaticObj = If[targetStatic,
        CompileObjective[targetExpr, space, basis, opbasis],
        None
    ];
    If[targetStatic && compiledStaticObj === $Failed, Return[$Failed]];

    memberComps = Table[
        Module[{rSub = rulesList[[mIdx]], rf = rfList[[mIdx]], dSym, cSyms, t0Raw, t0M,
                dComp, h0, dtFineM, nStagesM, isTDM, driftKindM, driftIntM, driftIntSpecM,
                dEvalM, dDerivM, stageMetaM, hCs, obj, rcSpecM, svAllowed, meth},
            dSym = driftExpr /. rSub;
            cSyms = ctrlList /. rSub;
            t0Raw = t0Opt /. rSub;
            If[!NumericQ[t0Raw] || !FreeQ[t0Raw, Complex | DirectedInfinity | Indeterminate] || Im[N[t0Raw]] != 0,
                Message[OptimalControlProblem::badt0, t0Raw]; Return[$Failed, Module]
            ];
            t0M = N[Re[t0Raw]];
            If[ListQ[hoistedDrift],
                h0 = hoistedDrift[[mIdx]];
                dtFineM = Developer`ToPackedArray[N[dtVec]];
                nStagesM = 1;
                isTDM = False;
                driftKindM = "Static";
                driftIntM = "Static";
                driftIntSpecM = None;
                dEvalM = None;
                dDerivM = None;
                stageMetaM = None,
                dComp = compileDriftFn[dSym, space, basis, opbasis, dim, hamFactor, dtVec, t0M, driftIntOpt];
                If[dComp === $Failed, Return[$Failed, Module]];
                h0 = dComp["Drift"];
                dtFineM = dComp["FineTimeSteps"];
                nStagesM = dComp["NumStages"];
                isTDM = dComp["IsTimeDependentDrift"];
                driftKindM = dComp["DriftKind"];
                driftIntM = dComp["DriftIntegration"];
                driftIntSpecM = Lookup[dComp, "DriftIntegrationSpec", None];
                dEvalM = dComp["DriftEvaluator"];
                dDerivM = dComp["DriftDerivEvaluator"];
                stageMetaM = dComp["StageMeta"]
            ];
            hCs = Table[
                Module[{baseC},
                    baseC = If[ListQ[hoistedCtrls[[k]]],
                        hoistedCtrls[[k, mIdx]],
                        With[{rawC = toGen[cSyms[[k]], space, basis, opbasis, dim]},
                            If[rawC === $Failed, Return[$Failed, Table]];
                            If[space === "Hilbert", HermitianSymmetrize[ToPackedComplexMatrix[ctrlFactor * rawC]],
                                ToPackedComplexMatrix[ctrlFactor * rawC]]
                        ]
                    ];
                    If[rf == 1., baseC, ToPackedComplexMatrix[rf * baseC]]
                ],
                {k, mCtrl}
            ];
            If[!FreeQ[hCs, $Failed], Return[$Failed, Module]];
            If[space === "Hilbert",
                hermErr = Max[Map[Max[Abs[# - ConjugateTranspose[#]]] / Max[1., Max[Abs[#]]] &, hCs]];
                If[hermErr > 10^-8, Message[OptimalControlProblem::nonherm, hermErr]; Return[$Failed, Module]]
            ];
            obj = If[targetStatic, compiledStaticObj, CompileObjective[targetExpr /. rSub, space, basis, opbasis]];
            If[obj === $Failed, Return[$Failed, Module]];
            rcSpecM = CompileRunningCost[runCostOpt /. rSub, obj, space, basis, opbasis, dtVec, t0M];
            If[rcSpecM === $Failed, Return[$Failed, Module]];
            If[AssociationQ[rcSpecM], rcSpecM = Join[rcSpecM, <|"InitialTimePoint" -> t0M|>]];
            svAllowed = obj["PureKetFastPath"] =!= None &&
                (rcSpecM === None || TrueQ[rcSpecM["CompatibleWithStateVectorQ"]]);
            meth = Which[
                methodOpt =!= Automatic, methodOpt,
                space === "Liouville", "Liouville",
                svAllowed && dim >= 16, "StateVector",
                True, "Spectral"
            ];
            If[(space === "Hilbert" && !MemberQ[{"Spectral", "BlockMatrixExp", "FirstOrder", "StateVector"}, meth]) ||
               (space === "Liouville" && meth =!= "Liouville") ||
               (meth === "StateVector" && !svAllowed),
                Message[OptimalControlProblem::badmethod, meth, space]; Return[$Failed, Module]
            ];
            <|
                "Space" -> space,
                "Method" -> meth,
                "Dimension" -> dim,
                "Basis" -> basis,
                "OperatorBasis" -> opbasis,
                "Drift" -> h0,
                "FineTimeSteps" -> dtFineM,
                "NumStages" -> nStagesM,
                "IsTimeDependentDrift" -> isTDM,
                "DriftKind" -> driftKindM,
                "DriftIntegration" -> driftIntM,
                "DriftIntegrationSpec" -> driftIntSpecM,
                "InitialTimePoint" -> t0M,
                "DriftEvaluator" -> dEvalM,
                "DriftDerivEvaluator" -> dDerivM,
                "DriftStageMeta" -> stageMetaM,
                "Controls" -> hCs,
                "ControlRows" -> ToPackedComplexMatrix[Flatten[#] & /@ hCs],
                "ControlFlatT" -> ToPackedComplexMatrix[Flatten[Transpose[#]] & /@ hCs],
                "SymbolicDrift" -> dSym,
                "SymbolicControls" -> If[rf == 1., cSyms, rf * cSyms],
                "DriftUnitFactor" -> hamFactor,
                "ControlUnitFactor" -> ctrlFactor,
                "FrequencyUnits" -> freqUnit,
                "ControlUnits" -> If[ctrlUnit === Automatic, freqUnit, ctrlUnit],
                "TimeSteps" -> dtVec,
                "TotalDuration" -> totalDur,
                "NumSlices" -> nt,
                "NumControls" -> mCtrl,
                "ControlBasis" -> ctrlBasis,
                "NumParams" -> nParams,
                "Objective" -> obj,
                "RunningCostSpec" -> rcSpecM,
                "RawTarget" -> (targetExpr /. rSub),
                "PowerPenalty" -> 0.,
                "SmoothnessPenalty" -> 0.
            |>
        ],
        {mIdx, nM}
    ];
    If[MemberQ[memberComps, $Failed], Return[$Failed]];

    If[optDurOpt =!= False && optDurOpt =!= None &&
       AnyTrue[memberComps, Lookup[#, "DriftKind", "Static"] === "SliceList" &],
        Message[OptimalControlProblem::nodrifttime]; Return[$Failed]
    ];
    durSpec = resolveDur[optDurOpt, durBndOpt, fixedTotOpt, dtVec, initOpt, ctrlBasis, cSpecRaw, "Static", mCtrl];
    If[durSpec === $Failed, Return[$Failed]];
    dpSpec = parseDurPen[durPenOpt, totalDur];
    If[dpSpec === $Failed, Return[$Failed]];
    cSpec = durSpec["ModifiedConstraintSpec"];
    hasDur = TrueQ[durSpec["HasDurationOptimization"]];
    proxCtx = compileProx[ctrlBasis, dtVec, cSpec, N[Max[0., l1Pen]], N[Max[0., tvPen]]];

    baseComp = First[memberComps];
    Module[{pfSpec, hasPhi, nPhi},
        pfSpec = Lookup[baseComp["Objective"], "PhaseFreedomSpec", None];
        hasPhi = AssociationQ[pfSpec] && TrueQ[pfSpec["ParametricQ"]] && TrueQ[pfSpec["NumPhases"] > 0];
        nPhi = If[hasPhi, pfSpec["NumPhases"], 0];
        If[!AllTrue[memberComps,
            With[{mPf = Lookup[#["Objective"], "PhaseFreedomSpec", None]},
                If[AssociationQ[mPf] && TrueQ[mPf["ParametricQ"]], Lookup[mPf, "NumPhases", 0], 0]
            ] === nPhi &],
            Message[OptimalControlEnsemble::incompat]; Return[$Failed]
        ];
        baseScales = computeEnsembleParameterScales[memberComps, dtVec, ctrlBasis["BasisMatrix"], mCtrl];
        ctrlScales = baseScales["ParameterScales"];
        fullScales = Developer`ToPackedArray[Join[
            ctrlScales,
            If[hasDur, ConstantArray[durSpec["DurationScale"], durSpec["NumActive"]], {}],
            If[hasPhi, ConstantArray[1., nPhi], {}]
        ]];
        lbOut = Developer`ToPackedArray[Join[
            cSpec["LowerBounds"],
            If[hasDur, durSpec["DurationLowerBounds"], {}],
            If[hasPhi, ConstantArray[-Infinity, nPhi], {}]
        ]];
        ubOut = Developer`ToPackedArray[Join[
            cSpec["UpperBounds"],
            If[hasDur, durSpec["DurationUpperBounds"], {}],
            If[hasPhi, ConstantArray[Infinity, nPhi], {}]
        ]];
        hasBoundsOut = TrueQ[cSpec["HasBox"] || hasDur];

        res = <|
            "Kind" -> "Ensemble",
            "Members" -> memberComps,
            "Weights" -> weights,
            "NumMembers" -> nM,
            "Aggregation" -> aggKind,
            "SoftMinBeta" -> aggBeta,
            "ParameterPoints" -> grid["ParameterPoints"],
            "GridShape" -> grid["GridShape"],
            "Space" -> space,
            "Method" -> baseComp["Method"],
            "Basis" -> basis,
            "OperatorBasis" -> opbasis,
            "TimeSteps" -> dtVec,
            "TotalDuration" -> totalDur,
            "NumSlices" -> nt,
            "NumControls" -> mCtrl,
            "ControlBasis" -> ctrlBasis,
            "NumParams" -> nParams,
            "PhaseFreedomSpec" -> pfSpec,
            "HasPhaseFreedom" -> hasPhi,
            "NumPhases" -> nPhi,
        "PowerPenalty" -> N[Max[0., powerPen]],
        "SmoothnessPenalty" -> N[Max[0., smoothPen]],
        "L1Penalty" -> N[Max[0., l1Pen]],
        "TVPenalty" -> N[Max[0., tvPen]],
        "ControlParameterScales" -> ctrlScales,
        "ParameterScales" -> fullScales,
        "ChannelNorms" -> baseScales["ChannelNorms"],
        "ControlLowerBounds" -> cSpec["LowerBounds"],
        "ControlUpperBounds" -> cSpec["UpperBounds"],
        "LowerBounds" -> lbOut,
        "UpperBounds" -> ubOut,
        "HasControlBounds" -> cSpec["HasBox"],
        "HasBounds" -> hasBoundsOut,
        "HasDiskBounds" -> cSpec["HasDisk"],
        "HasSlewBounds" -> cSpec["HasSlew"],
        "HasConstraints" -> TrueQ[cSpec["HasAnyConstraint"] || hasDur],
        "ConstraintSpec" -> cSpec,
        "ProxContext" -> proxCtx,
        "HasDurationOptimization" -> hasDur,
        "DurationMode" -> durSpec["DurationMode"],
        "DurationWeights" -> durSpec["DurationWeights"],
        "OptimizeDurations" -> optDurOpt,
        "ActiveDurationIndices" -> durSpec["ActiveIndices"],
        "FreezeControlSlices" -> durSpec["FreezeControlSlices"],
        "NumActiveDurations" -> durSpec["NumActive"],
        "DurationLowerBounds" -> durSpec["DurationLowerBounds"],
        "DurationUpperBounds" -> durSpec["DurationUpperBounds"],
        "FixedTotalDuration" -> durSpec["FixedTotalDuration"],
        "FixedTotalTarget" -> durSpec["FixedTotalTarget"],
        "ActiveDurationTargetSum" -> durSpec["ActiveTargetSum"],
        "DurationScale" -> durSpec["DurationScale"],
        "DurationPenaltyLambda" -> dpSpec["Lambda"],
        "DurationPenaltyReference" -> dpSpec["Reference"],
        "SymbolicDrift" -> baseComp["SymbolicDrift"],
        "SymbolicControls" -> baseComp["SymbolicControls"],
        "DriftUnitFactor" -> hamFactor,
        "ControlUnitFactor" -> ctrlFactor,
        "FrequencyUnits" -> freqUnit,
        "ControlUnits" -> baseComp["ControlUnits"],
        "InitialTimePoint" -> Lookup[baseComp, "InitialTimePoint", 0.],
        "DriftIntegration" -> Lookup[baseComp, "DriftIntegration", "Static"],
        "IsTimeDependentDrift" -> AnyTrue[memberComps, TrueQ[Lookup[#, "IsTimeDependentDrift", False]] &],
        "UseParallel" -> False,
        "OriginEnsemble" -> ens
    |>] ;
    If[TrueQ[parOpt] || (IntegerQ[parOpt] && parOpt >= 2),
        setupParallelWorkers[res, parOpt],
        res
    ]
];

mergeMemberRunningCost[rc1_, rc2_] := Which[
    rc1 === None, rc2,
    rc2 === None, rc1,
    True,
        Module[{svOK = TrueQ[rc1["CompatibleWithStateVectorQ"]] && TrueQ[rc2["CompatibleWithStateVectorQ"]]},
            <|
                "Space" -> rc1["Space"],
                "Dimension" -> rc1["Dimension"],
                "NumTerms" -> rc1["NumTerms"] + rc2["NumTerms"],
                "CompatibleWithStateVectorQ" -> svOK,
                "HasTimeDependence" -> TrueQ[rc1["HasTimeDependence"]] || TrueQ[rc2["HasTimeDependence"]],
                "TWSpecs" -> Join[rc1["TWSpecs"], rc2["TWSpecs"]],
                "HilbertTermEvaluators" -> Join[rc1["HilbertTermEvaluators"], rc2["HilbertTermEvaluators"]],
                "StateVectorTermEvaluators" -> If[svOK,
                    Join[rc1["StateVectorTermEvaluators"], rc2["StateVectorTermEvaluators"]], None],
                "LiouvilleTermEvaluators" -> Join[rc1["LiouvilleTermEvaluators"], rc2["LiouvilleTermEvaluators"]],
                "CustomEvaluators" -> Join[rc1["CustomEvaluators"], rc2["CustomEvaluators"]],
                "InitialTimePoint" -> Lookup[rc1, "InitialTimePoint", 0.]
            |>
        ]
];

(* Compile a heterogeneous list of OptimalControlProblem instances *)
compileHeterogeneousEnsemble[probsOrWeighted_List, opts___Rule] := Module[
    {mergedOpts, normWOpt, aggOpt, parOpt, runCostOpt, optDurOpt, durBndOpt, durPenOpt, fixedTotOpt, initOpt,
     powerPen, smoothPen, l1Pen, tvPen, ctrlBounds, ampBound, ampNormBound, slewBound,
     rawW, rawProbs, weights, aggKind, aggBeta, memberComps, baseComp, nM, dt0, basis0,
     mCtrl0, nParams0, nt0, resolveCSpec, resolveDur, parseDurPen, compileProx, cSpecRaw, durSpec, dpSpec, hasDur,
     cSpec, proxCtx, baseScales, ctrlScales, fullScales, lbOut, ubOut, hasBoundsOut, res},

    If[probsOrWeighted === {}, Message[OptimalControlEnsemble::inval, probsOrWeighted]; Return[$Failed]];
    mergedOpts = Join[{opts}, Options[OptimalControlEnsemble]];
    {normWOpt, aggOpt, parOpt, runCostOpt, optDurOpt, durBndOpt, durPenOpt, fixedTotOpt, initOpt,
     powerPen, smoothPen, l1Pen, tvPen, ctrlBounds, ampBound, ampNormBound, slewBound} =
        {NormalizeWeights, EnsembleAggregation, Parallel, RunningCost,
         OptimizeDurations, DurationBounds, DurationPenalty, FixedTotalDuration, InitialControls,
         PowerPenalty, SmoothnessPenalty, L1Penalty, TVPenalty,
         ControlBounds, AmplitudeBound, AmplitudeNormBound, SlewRateBound} /. mergedOpts;

    Which[
        aggOpt === "Mean" || aggOpt === Automatic, aggKind = "Mean"; aggBeta = 0.,
        MatchQ[aggOpt, {"SoftMin", _?NumericQ}] && TrueQ[aggOpt[[2]] > 0],
            aggKind = "SoftMin"; aggBeta = N[aggOpt[[2]]],
        True, Message[OptimalControlEnsemble::badagg, aggOpt]; Return[$Failed]
    ];

    {rawW, rawProbs} = Transpose[Map[
        Function[item,
            Which[
                MatchQ[item, _OptimalControlProblem], {1., item},
                MatchQ[item, {_?NumericQ, _OptimalControlProblem}], {N[First[item]], item[[2]]},
                True, {$Failed, $Failed}
            ]
        ],
        probsOrWeighted
    ]];
    If[!FreeQ[rawProbs, $Failed] || Total[rawW] <= 0.,
        Message[OptimalControlEnsemble::inval, probsOrWeighted]; Return[$Failed]
    ];
    weights = Developer`ToPackedArray[If[TrueQ[normWOpt], N[rawW / Total[rawW]], N[rawW]]];
    memberComps = CompileOptimalControlProblem /@ rawProbs;
    If[MemberQ[memberComps, $Failed], Return[$Failed]];

    If[runCostOpt =!= None && runCostOpt =!= False && runCostOpt =!= {},
        memberComps = Map[
            Function[mem,
                Module[{t0M = Lookup[mem, "InitialTimePoint", 0.], rcEns, mergedRC, newMeth = mem["Method"]},
                    rcEns = CompileRunningCost[runCostOpt, mem["Objective"], mem["Space"],
                        mem["Basis"], mem["OperatorBasis"], mem["TimeSteps"], t0M];
                    If[rcEns === $Failed, Return[$Failed]];
                    If[AssociationQ[rcEns], rcEns = Join[rcEns, <|"InitialTimePoint" -> t0M|>]];
                    mergedRC = mergeMemberRunningCost[Lookup[mem, "RunningCostSpec", None], rcEns];
                    If[newMeth === "StateVector" && mergedRC =!= None && !TrueQ[mergedRC["CompatibleWithStateVectorQ"]],
                        If[MemberQ[Lookup[mem, "Options", {}], Method -> "StateVector"],
                            Message[OptimalControlProblem::badmethod, newMeth, mem["Space"]]; Return[$Failed]
                        ];
                        newMeth = "Spectral"
                    ];
                    Join[mem, <|"RunningCostSpec" -> mergedRC, "Method" -> newMeth|>]
                ]
            ],
            memberComps
        ];
        If[MemberQ[memberComps, $Failed], Return[$Failed]]
    ];

    nM = Length[memberComps];
    baseComp = First[memberComps];
    mCtrl0 = baseComp["NumControls"];
    nt0 = baseComp["NumSlices"];
    nParams0 = baseComp["NumParams"];
    dt0 = baseComp["TimeSteps"];
    basis0 = baseComp["ControlBasis", "BasisMatrix"];

    If[!AllTrue[memberComps,
        #["NumControls"] === mCtrl0 && #["NumSlices"] === nt0 && #["NumParams"] === nParams0 &&
        Max[Abs[#["TimeSteps"] - dt0]] <= 10^-12 * Max[dt0] &&
        Max[Abs[#["ControlBasis", "BasisMatrix"] - basis0]] <= 10^-12 &],
        Message[OptimalControlEnsemble::incompat]; Return[$Failed]
    ];

    If[(optDurOpt === False || optDurOpt === None) && TrueQ[Lookup[baseComp, "HasDurationOptimization", False]],
        {optDurOpt, durBndOpt, fixedTotOpt} =
            {OptimizeDurations, DurationBounds, FixedTotalDuration} /. Lookup[baseComp, "Options", {}] /. Options[OptimalControlEnsemble]
    ];
    If[(durPenOpt === 0 || durPenOpt === 0. || durPenOpt === None || durPenOpt === False) &&
       TrueQ[Lookup[baseComp, "DurationPenaltyLambda", 0.] > 0.],
        durPenOpt = DurationPenalty /. Lookup[baseComp, "Options", {}] /. Options[OptimalControlEnsemble]
    ];
    If[optDurOpt =!= False && optDurOpt =!= None &&
       AnyTrue[memberComps, Lookup[#, "DriftKind", "Static"] === "SliceList" &],
        Message[OptimalControlProblem::nodrifttime]; Return[$Failed]
    ];

    resolveCSpec = $OCSolverRegistry["ResolveConstraintSpec"];
    resolveDur = $OCSolverRegistry["ResolveDurationSpec"];
    parseDurPen = $OCSolverRegistry["ParseDurationPenalty"];
    compileProx = $OCSolverRegistry["CompileProxContext"];
    cSpecRaw = resolveCSpec[ctrlBounds, ampBound, ampNormBound, slewBound, dt0, nParams0, mCtrl0];
    If[cSpecRaw === $Failed, Return[$Failed]];
    durSpec = resolveDur[optDurOpt, durBndOpt, fixedTotOpt, dt0, initOpt, baseComp["ControlBasis"], cSpecRaw, "Static", mCtrl0];
    If[durSpec === $Failed, Return[$Failed]];
    dpSpec = parseDurPen[durPenOpt, baseComp["TotalDuration"]];
    If[dpSpec === $Failed, Return[$Failed]];
    cSpec = durSpec["ModifiedConstraintSpec"];
    hasDur = TrueQ[durSpec["HasDurationOptimization"]];
    proxCtx = compileProx[baseComp["ControlBasis"], dt0, cSpec, N[Max[0., l1Pen]], N[Max[0., tvPen]]];

    Module[{pfSpec, hasPhi, nPhi},
        pfSpec = Lookup[baseComp["Objective"], "PhaseFreedomSpec", None];
        hasPhi = AssociationQ[pfSpec] && TrueQ[pfSpec["ParametricQ"]] && TrueQ[pfSpec["NumPhases"] > 0];
        nPhi = If[hasPhi, pfSpec["NumPhases"], 0];
        If[!AllTrue[memberComps,
            With[{mPf = Lookup[#["Objective"], "PhaseFreedomSpec", None]},
                If[AssociationQ[mPf] && TrueQ[mPf["ParametricQ"]], Lookup[mPf, "NumPhases", 0], 0]
            ] === nPhi &],
            Message[OptimalControlEnsemble::incompat]; Return[$Failed]
        ];
        baseScales = computeEnsembleParameterScales[memberComps, dt0, basis0, mCtrl0];
        ctrlScales = baseScales["ParameterScales"];
        fullScales = Developer`ToPackedArray[Join[
            ctrlScales,
            If[hasDur, ConstantArray[durSpec["DurationScale"], durSpec["NumActive"]], {}],
            If[hasPhi, ConstantArray[1., nPhi], {}]
        ]];
        lbOut = Developer`ToPackedArray[Join[
            cSpec["LowerBounds"],
            If[hasDur, durSpec["DurationLowerBounds"], {}],
            If[hasPhi, ConstantArray[-Infinity, nPhi], {}]
        ]];
        ubOut = Developer`ToPackedArray[Join[
            cSpec["UpperBounds"],
            If[hasDur, durSpec["DurationUpperBounds"], {}],
            If[hasPhi, ConstantArray[Infinity, nPhi], {}]
        ]];
        hasBoundsOut = TrueQ[cSpec["HasBox"] || hasDur];

        res = <|
            "Kind" -> "Ensemble",
            "Members" -> memberComps,
            "Weights" -> weights,
            "NumMembers" -> nM,
            "Aggregation" -> aggKind,
            "SoftMinBeta" -> aggBeta,
            "ParameterPoints" -> None,
            "GridShape" -> {nM},
            "Space" -> baseComp["Space"],
            "Method" -> baseComp["Method"],
            "Basis" -> baseComp["Basis"],
            "OperatorBasis" -> baseComp["OperatorBasis"],
            "TimeSteps" -> dt0,
            "TotalDuration" -> baseComp["TotalDuration"],
            "NumSlices" -> nt0,
            "NumControls" -> mCtrl0,
            "ControlBasis" -> baseComp["ControlBasis"],
            "NumParams" -> nParams0,
            "PhaseFreedomSpec" -> pfSpec,
            "HasPhaseFreedom" -> hasPhi,
            "NumPhases" -> nPhi,
        "PowerPenalty" -> N[Max[0., powerPen]],
        "SmoothnessPenalty" -> N[Max[0., smoothPen]],
        "L1Penalty" -> N[Max[0., l1Pen]],
        "TVPenalty" -> N[Max[0., tvPen]],
        "ControlParameterScales" -> ctrlScales,
        "ParameterScales" -> fullScales,
        "ChannelNorms" -> baseScales["ChannelNorms"],
        "ControlLowerBounds" -> cSpec["LowerBounds"],
        "ControlUpperBounds" -> cSpec["UpperBounds"],
        "LowerBounds" -> lbOut,
        "UpperBounds" -> ubOut,
        "HasControlBounds" -> cSpec["HasBox"],
        "HasBounds" -> hasBoundsOut,
        "HasDiskBounds" -> cSpec["HasDisk"],
        "HasSlewBounds" -> cSpec["HasSlew"],
        "HasConstraints" -> TrueQ[cSpec["HasAnyConstraint"] || hasDur],
        "ConstraintSpec" -> cSpec,
        "ProxContext" -> proxCtx,
        "HasDurationOptimization" -> hasDur,
        "DurationMode" -> durSpec["DurationMode"],
        "DurationWeights" -> durSpec["DurationWeights"],
        "OptimizeDurations" -> optDurOpt,
        "ActiveDurationIndices" -> durSpec["ActiveIndices"],
        "FreezeControlSlices" -> durSpec["FreezeControlSlices"],
        "NumActiveDurations" -> durSpec["NumActive"],
        "DurationLowerBounds" -> durSpec["DurationLowerBounds"],
        "DurationUpperBounds" -> durSpec["DurationUpperBounds"],
        "FixedTotalDuration" -> durSpec["FixedTotalDuration"],
        "FixedTotalTarget" -> durSpec["FixedTotalTarget"],
        "ActiveDurationTargetSum" -> durSpec["ActiveTargetSum"],
        "DurationScale" -> durSpec["DurationScale"],
        "DurationPenaltyLambda" -> dpSpec["Lambda"],
        "DurationPenaltyReference" -> dpSpec["Reference"],
        "SymbolicDrift" -> baseComp["SymbolicDrift"],
        "SymbolicControls" -> baseComp["SymbolicControls"],
        "DriftUnitFactor" -> baseComp["DriftUnitFactor"],
        "ControlUnitFactor" -> baseComp["ControlUnitFactor"],
        "FrequencyUnits" -> baseComp["FrequencyUnits"],
        "ControlUnits" -> baseComp["ControlUnits"],
        "InitialTimePoint" -> Lookup[baseComp, "InitialTimePoint", 0.],
        "DriftIntegration" -> Lookup[baseComp, "DriftIntegration", "Static"],
        "IsTimeDependentDrift" -> AnyTrue[memberComps, TrueQ[Lookup[#, "IsTimeDependentDrift", False]] &],
        "UseParallel" -> False,
        "OriginEnsemble" -> OptimalControlEnsemble[probsOrWeighted, opts]
    |>] ;
    If[TrueQ[parOpt] || (IntegerQ[parOpt] && parOpt >= 2),
        setupParallelWorkers[res, parOpt],
        res
    ]
];

compileEnsembleIfNeeded[ens:OptimalControlEnsemble[items_List, opts___Rule], extraOpts___Rule] :=
    compileHeterogeneousEnsemble[items, Sequence @@ {extraOpts, opts}];

compileEnsembleIfNeeded[ens:OptimalControlEnsemble[drift_, controls_, target_, timing_, opts___Rule], extraOpts___Rule] :=
    compileParametricEnsemble[ens, extraOpts];

compileEnsembleIfNeeded[other_, ___] := other;

(* Evaluate ensemble fidelity + aggregated exact analytical gradient *)
evaluateCompiledEnsemble[compEns_Association, input_] := Module[
    {nt, mCtrl, nParams, nCtrl, basisAssoc, isIdent, hasDur, durMode, activeIdx, freezeIdx, nAct,
     hasPhi, nPhi, phiVec, dtVec, jointFlag = False, pMat, uMat, members, weights, nM, tTot,
     aggKind, beta, buildSched, buildTimingTensors, stepOverrides = None,
     staticTimingCache = <||>, memberResults, fids, runVals, grads, durGrads, phiGrads, leakVals,
     memberScores, sMin, fMin, expW, effWeights, aggScore, aggTermFid, aggRunVal, softMinFid,
     gradU, gradUFid, gDtKern = None, gPhi = None, lamP, lamS, lamT, tRef, pen, gRegU, gDtPen, diffDt, du, duSq,
     gU, gParam, gParamFid, gRegParam, gDt, gDurActive, gDurFidActive, gRegDurActive, gTotDur, gradOut,
     kernels, chunks, cid, workerOut, ovMap},

    nt = compEns["NumSlices"];
    mCtrl = compEns["NumControls"];
    nParams = compEns["NumParams"];
    nCtrl = nParams * mCtrl;
    basisAssoc = compEns["ControlBasis"];
    isIdent = basisAssoc["IsIdentityQ"];
    hasDur = TrueQ[Lookup[compEns, "HasDurationOptimization", False]];
    durMode = Lookup[compEns, "DurationMode", If[hasDur, "Slices", "None"]];
    activeIdx = Lookup[compEns, "ActiveDurationIndices", {}];
    freezeIdx = Lookup[compEns, "FreezeControlSlices", {}];
    nAct = If[hasDur, Lookup[compEns, "NumActiveDurations", 0], 0];
    hasPhi = TrueQ[Lookup[compEns, "HasPhaseFreedom", False]];
    nPhi = If[hasPhi, Lookup[compEns, "NumPhases", 0], 0];
    dtVec = Developer`ToPackedArray[N[compEns["TimeSteps"]]];
    phiVec = Developer`ToPackedArray[ConstantArray[0., nPhi]];

    pMat = Which[
        (hasDur || hasPhi) && VectorQ[input, NumericQ] && Length[input] === nCtrl + nAct + nPhi,
            jointFlag = True;
            If[hasDur,
                If[durMode === "Total",
                    tTot = N[Re[input[[nCtrl + 1]]]];
                    If[!TrueQ[tTot > 0.], Return[$Failed]];
                    dtVec = Developer`ToPackedArray[tTot * compEns["DurationWeights"]],
                    dtVec[[activeIdx]] = Developer`ToPackedArray[N[Re[input[[nCtrl + 1 ;; nCtrl + nAct]]]]]
                ]
            ];
            If[hasPhi, phiVec = Developer`ToPackedArray[N[Re[input[[nCtrl + nAct + 1 ;; -1]]]]]];
            Partition[Developer`ToPackedArray[N[Re[input[[1 ;; nCtrl]]]]], mCtrl],
        (hasDur || hasPhi) && AssociationQ[input] &&
        (KeyExistsQ[input, "TimeSteps"] || KeyExistsQ[input, "TotalDuration"] || KeyExistsQ[input, "Phases"] || KeyExistsQ[input, "OptimalPhases"]),
            Module[{cIn = Lookup[input, "Parameters", Lookup[input, "Controls", $Failed]],
                    dtIn = Lookup[input, "TimeSteps", None],
                    tIn = Lookup[input, "TotalDuration", None],
                    phIn = Lookup[input, "Phases", If[AssociationQ[Lookup[input, "OptimalPhases", None]], input[["OptimalPhases", "Vector"]], phiVec]]},
                If[cIn === $Failed, Return[$Failed]];
                If[hasDur,
                    Which[
                        VectorQ[dtIn, NumericQ] && Length[dtIn] === nt,
                            dtVec = Developer`ToPackedArray[N[Re[dtIn]]],
                        durMode === "Total" && NumericQ[tIn] && TrueQ[Re[tIn] > 0.],
                            dtVec = Developer`ToPackedArray[N[Re[tIn]] * compEns["DurationWeights"]],
                        dtIn =!= None || tIn =!= None,
                            Return[$Failed]
                    ]
                ];
                If[hasPhi,
                    If[!VectorQ[phIn, NumericQ] || Length[phIn] =!= nPhi, Return[$Failed]];
                    phiVec = Developer`ToPackedArray[N[Re[phIn]]]
                ];
                jointFlag = True;
                Which[
                    MatrixQ[cIn, NumericQ] && Dimensions[cIn] === {nParams, mCtrl}, ToPackedRealMatrix[cIn],
                    VectorQ[cIn, NumericQ] && Length[cIn] === nCtrl, Partition[Developer`ToPackedArray[N[Re[cIn]]], mCtrl],
                    !isIdent && MatrixQ[cIn, NumericQ] && Dimensions[cIn] === {nt, mCtrl},
                        LeastSquares[basisAssoc["BasisMatrix"], ToPackedRealMatrix[cIn]],
                    !isIdent && VectorQ[cIn, NumericQ] && Length[cIn] === nt * mCtrl,
                        LeastSquares[basisAssoc["BasisMatrix"], Partition[Developer`ToPackedArray[N[Re[cIn]]], mCtrl]],
                    True, $Failed
                ]
            ],
        VectorQ[input, NumericQ] && Length[input] === nCtrl,
            Partition[Developer`ToPackedArray[N[Re[input]]], mCtrl],
        MatrixQ[input, NumericQ] && Dimensions[input] === {nParams, mCtrl},
            ToPackedRealMatrix[input],
        !isIdent && MatrixQ[input, NumericQ] && Dimensions[input] === {nt, mCtrl},
            LeastSquares[basisAssoc["BasisMatrix"], ToPackedRealMatrix[input]],
        !isIdent && VectorQ[input, NumericQ] && Length[input] === nt * mCtrl,
            LeastSquares[basisAssoc["BasisMatrix"], Partition[Developer`ToPackedArray[N[Re[input]]], mCtrl]],
        True,
            $Failed
    ];
    If[pMat === $Failed, Return[$Failed]];
    uMat = If[isIdent, pMat, Developer`ToPackedArray[basisAssoc["BasisMatrix"] . pMat]];
    If[freezeIdx =!= {} && !isIdent, Return[$Failed]];
    If[freezeIdx =!= {},
        uMat = ToPackedRealMatrix[uMat];
        uMat[[freezeIdx, All]] = 0.
    ];
    If[hasDur && !AllTrue[dtVec, TrueQ[# > 0.] &], Return[$Failed]];

    members = compEns["Members"];
    weights = compEns["Weights"];
    nM = compEns["NumMembers"];
    aggKind = compEns["Aggregation"];
    beta = compEns["SoftMinBeta"];

    If[hasDur,
        buildSched = $OCSolverRegistry["BuildFineDriftSchedule"];
        buildTimingTensors = $OCSolverRegistry["BuildTimingDerivTensors"];
        stepOverrides = Table[
            Module[{mem = members[[mIdx]], t0M, sched, tSpec},
                t0M = Lookup[mem, "InitialTimePoint", 0.];
                If[Lookup[mem, "DriftKind", "Static"] === "Continuous",
                    sched = buildSched[mem["DriftEvaluator"], dtVec, t0M, mem["DriftIntegrationSpec"], mem["Space"]];
                    If[sched === $Failed, Return[$Failed, Module]];
                    tSpec = buildTimingTensors[mem, dtVec, sched];
                    <|"Drift" -> sched["Drift"], "FineTimeSteps" -> sched["FineTimeSteps"], "TimingSpec" -> tSpec|>,
                    tSpec = If[KeyExistsQ[staticTimingCache, t0M],
                        staticTimingCache[t0M],
                        staticTimingCache[t0M] = buildTimingTensors[mem, dtVec, None]
                    ];
                    <|"Drift" -> None, "FineTimeSteps" -> dtVec, "TimingSpec" -> tSpec|>
                ]
            ],
            {mIdx, nM}
        ];
        If[MemberQ[stepOverrides, $Failed], Return[$Failed]]
    ];

    If[TrueQ[compEns["UseParallel"]],
        kernels = compEns["WorkerKernels"];
        chunks = compEns["WorkerChunks"];
        cid = compEns["WorkerCacheId"];
        ovMap = If[hasDur,
            AssociationThread[compEns["WorkerKernelIDs"], Map[stepOverrides[[#]] &, chunks]],
            <||>
        ];
        workerOut = With[{wId = cid, uIn = uMat, ovById = ovMap, hd = hasDur, hp = hasPhi, pv = If[hasPhi, phiVec, None]},
            ParallelEvaluate[
                With[{mList = OptimalControl`Ensemble`Private`$WorkerMemberStore[wId]},
                    With[{oList = If[hd, ovById[$KernelID], ConstantArray[None, Length[mList]]]},
                        MapThread[
                            With[{r = evalMemberKernel[#1, uIn, #2, pv]},
                                With[{ext = Lookup[r, "Extra", None]},
                                    {r["Fidelity"], Lookup[r, "RunningCost", 0.], r["Gradient"],
                                     Lookup[r, "DurationGradient", None],
                                     If[hp && AssociationQ[ext], Lookup[ext, "PhaseGradient", None], None],
                                     If[AssociationQ[ext], Lookup[ext, "Leakage", None], None]}
                                ]
                            ] &,
                            {mList, oList}
                        ]
                    ]
                ],
                kernels
            ]
        ];
        fids = ConstantArray[0., nM];
        runVals = ConstantArray[0., nM];
        grads = ConstantArray[0., {nM, nt, mCtrl}];
        If[hasDur, durGrads = ConstantArray[0., {nM, nt}]];
        If[hasPhi, phiGrads = ConstantArray[0., {nM, nPhi}]];
        leakVals = ConstantArray[None, nM];
        Do[
            Do[
                fids[[chunks[[w, pos]]]] = workerOut[[w, pos, 1]];
                runVals[[chunks[[w, pos]]]] = workerOut[[w, pos, 2]];
                grads[[chunks[[w, pos]]]] = workerOut[[w, pos, 3]];
                If[hasDur, durGrads[[chunks[[w, pos]]]] = workerOut[[w, pos, 4]]];
                If[hasPhi, phiGrads[[chunks[[w, pos]]]] = workerOut[[w, pos, 5]]];
                leakVals[[chunks[[w, pos]]]] = workerOut[[w, pos, 6]],
                {pos, Length[chunks[[w]]]}
            ],
            {w, Length[chunks]}
        ],
        memberResults = With[{pv = If[hasPhi, phiVec, None]},
            If[hasDur,
                MapThread[evalMemberKernel[#1, uMat, #2, pv] &, {members, stepOverrides}],
                Map[evalMemberKernel[#, uMat, None, pv] &, members]
            ]
        ];
        If[!AllTrue[memberResults, AssociationQ], Return[$Failed]];
        fids = Developer`ToPackedArray[N[#["Fidelity"] & /@ memberResults]];
        runVals = Developer`ToPackedArray[N[Lookup[#, "RunningCost", 0.] & /@ memberResults]];
        grads = Developer`ToPackedArray[#["Gradient"] & /@ memberResults];
        If[hasDur, durGrads = Developer`ToPackedArray[#["DurationGradient"] & /@ memberResults]];
        If[hasPhi,
            phiGrads = Developer`ToPackedArray[Map[
                Function[mr, With[{ext = Lookup[mr, "Extra", None]},
                    If[AssociationQ[ext], Lookup[ext, "PhaseGradient", ConstantArray[0., nPhi]], ConstantArray[0., nPhi]]
                ]],
                memberResults
            ]]
        ];
        leakVals = Map[
            Function[mr, With[{ext = Lookup[mr, "Extra", None]},
                If[AssociationQ[ext], Lookup[ext, "Leakage", None], None]
            ]],
            memberResults
        ]
    ];

    memberScores = fids + runVals;
    If[aggKind === "SoftMin",
        sMin = Min[memberScores];
        expW = weights * Exp[-beta * (memberScores - sMin)];
        effWeights = expW / Total[expW];
        aggScore = sMin - (1. / beta) * Log[Total[expW]];
        fMin = Min[fids];
        aggTermFid = fMin - (1. / beta) * Log[Total[weights * Exp[-beta * (fids - fMin)]]];
        aggRunVal = aggScore - aggTermFid,
        effWeights = weights;
        aggTermFid = N[weights . fids];
        aggRunVal = N[weights . runVals];
        aggScore = aggTermFid + aggRunVal
    ];
    gradU = Developer`ToPackedArray[effWeights . grads];
    If[hasDur, gDtKern = Developer`ToPackedArray[effWeights . durGrads]];
    If[hasPhi, gPhi = Developer`ToPackedArray[effWeights . phiGrads]];

    lamP = compEns["PowerPenalty"];
    lamS = compEns["SmoothnessPenalty"];
    lamT = Lookup[compEns, "DurationPenaltyLambda", 0.];
    tRef = Max[Lookup[compEns, "DurationPenaltyReference", 1.], 10^-12];
    pen = 0.;
    gRegU = ConstantArray[0., {nt, mCtrl}];
    gDtPen = ConstantArray[0., nt];
    If[lamP > 0.,
        With[{uSq = Total[uMat * uMat, {2}]},
            pen += lamP * Total[dtVec * uSq];
            gRegU += (-2. * lamP) * (dtVec * uMat);
            If[hasDur, gDtPen -= lamP * uSq]
        ]
    ];
    If[lamS > 0. && nt > 1,
        diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]);
        du = Rest[uMat] - Most[uMat];
        duSq = Total[du * du, {2}];
        pen += lamS * Total[duSq / diffDt];
        gRegU[[1 ;; nt - 1]] += (2. * lamS) * (du / diffDt);
        gRegU[[2 ;; nt]] -= (2. * lamS) * (du / diffDt);
        If[hasDur,
            gDtPen[[1 ;; nt - 1]] += (0.5 * lamS) * (duSq / (diffDt^2));
            gDtPen[[2 ;; nt]] += (0.5 * lamS) * (duSq / (diffDt^2))
        ]
    ];
    If[lamT > 0.,
        pen += lamT * (Total[dtVec] / tRef);
        If[hasDur, gDtPen -= ConstantArray[lamT / tRef, nt]]
    ];
    softMinFid = If[aggKind === "SoftMin", aggTermFid, N[weights . fids]];
    With[{penNS = If[KeyExistsQ[compEns, "ProxContext"], ($OCSolverRegistry["EvalNonSmoothPenalty"])[compEns["ProxContext"], uMat], 0.]},
        gradUFid = gradU;
        If[freezeIdx =!= {}, gradUFid[[freezeIdx, All]] = 0.];
        gParamFid = If[isIdent, gradUFid, Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gradUFid]];
        gU = Developer`ToPackedArray[gradU + gRegU];
        If[freezeIdx =!= {}, gU[[freezeIdx, All]] = 0.];
        gParam = If[isIdent, gU, Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gU]];
        gRegParam = Developer`ToPackedArray[gParam - gParamFid];
        gDt = If[hasDur, Developer`ToPackedArray[gDtKern + gDtPen], None];
        gDurActive = Which[
            !hasDur, {},
            durMode === "Total", Developer`ToPackedArray[{(compEns["DurationWeights"]) . gDt}],
            True, Developer`ToPackedArray[gDt[[activeIdx]]]
        ];
        gDurFidActive = Which[
            !hasDur, {},
            durMode === "Total", Developer`ToPackedArray[{(compEns["DurationWeights"]) . gDtKern}],
            True, Developer`ToPackedArray[gDtKern[[activeIdx]]]
        ];
        gRegDurActive = Developer`ToPackedArray[gDurActive - gDurFidActive];
        gTotDur = If[hasDur && durMode === "Total", First[gDurActive], None];
        gradOut = If[(hasDur || hasPhi) && jointFlag,
            Developer`ToPackedArray[Join[
                Flatten[gParam],
                gDurActive,
                If[hasPhi, gPhi, {}]
            ]],
            gParam
        ];

        <|
            "Fidelity" -> aggTermFid,
            "TerminalFidelity" -> aggTermFid,
            "RunningCost" -> aggRunVal,
            "Leakage" -> If[AllTrue[leakVals, NumericQ], N[weights . leakVals], None],
            "SmoothObjective" -> aggScore - pen,
            "SmoothPenalty" -> pen,
            "NonSmoothPenalty" -> penNS,
            "Objective" -> aggScore - pen - penNS,
            "Penalty" -> pen + penNS,
            "Gradient" -> gradOut,
            "ControlGradient" -> gParam,
            "SliceGradient" -> gU,
            "FidelitySliceGradient" -> gradUFid,
            "DurationGradient" -> gDt,
            "ActiveDurationGradient" -> gDurActive,
            "TotalDurationGradient" -> gTotDur,
            "FidelityControlGradient" -> gParamFid,
            "FidelityActiveDurationGradient" -> gDurFidActive,
            "SmoothPenaltyControlGradient" -> gRegParam,
            "SmoothPenaltyActiveDurationGradient" -> gRegDurActive,
            "PhaseGradient" -> gPhi,
            "Phases" -> If[hasPhi, phiVec, None],
            "Parameters" -> pMat,
            "Controls" -> uMat,
            "TimeSteps" -> dtVec,
            "TotalDuration" -> Total[dtVec],
            "MemberFidelities" -> fids,
            "MemberRunningCosts" -> runVals,
            "MeanFidelity" -> N[weights . fids],
            "SoftMinFidelity" -> softMinFid,
            "MinFidelity" -> Min[fids],
            "MaxFidelity" -> Max[fids]
        |>
    ]
];

(* Phase C Ensemble Second-Order Capability, Cached Evaluator, HVP, and Explicit Hessian *)
checkEnsembleSecondOrderCapabilities[compEns_Association] := Module[
    {hasDur, hasPhi, hasNS, members, memberCaps, badCap, badReason, allSpec},
    hasDur = TrueQ[Lookup[compEns, "HasDurationOptimization", False]];
    hasPhi = TrueQ[Lookup[compEns, "HasPhaseFreedom", False]];
    hasNS = TrueQ[Lookup[compEns, "L1Penalty", 0.] > 0. || Lookup[compEns, "TVPenalty", 0.] > 0.];
    members = Lookup[compEns, "Members", {}];
    memberCaps = Table[
        With[{m = members[[i]],
              hasRC = AssociationQ[Lookup[members[[i]], "RunningCostSpec", None]],
              method = Lookup[members[[i]], "Method", "Unknown"],
              objOK = TrueQ[Lookup[members[[i]], "Objective", <||>] // (Lookup[#, "ExactHessianSupportedQ", False] &)]},
            Which[
                hasRC, <|"OK" -> False, "Spec" -> False, "Reason" -> "Second-order derivatives are not yet enabled for RunningCost / StatePenalty."|>,
                !MemberQ[{"Spectral", "StateVector", "Liouville"}, method],
                    <|"OK" -> False, "Spec" -> False, "Reason" -> "Second-order derivatives require Method -> \"Spectral\", \"StateVector\", or \"Liouville\" (got " <> ToString[method] <> ")."|>,
                !objOK,
                    <|"OK" -> False, "Spec" -> False, "Reason" -> Lookup[members[[i]]["Objective"], "UnsupportedSecondOrderReason", "Ensemble member objective does not support exact second derivatives."]|>,
                True,
                    <|"OK" -> True, "Spec" -> (method === "Spectral"), "Reason" -> None|>
            ]
        ],
        {i, Length[members]}
    ];
    badCap = SelectFirst[memberCaps, !TrueQ[#["OK"]] &, None];
    allSpec = AllTrue[memberCaps, TrueQ[#["Spec"]] &];
    badReason = Which[
        hasDur, "Second-order derivatives are not yet enabled for OptimizeDurations -> True.",
        hasPhi, "Second-order derivatives are not yet enabled for parametric PhaseFreedom.",
        hasNS, "Nonsmooth L1Penalty / TVPenalty problems do not have a smooth second-order Hessian.",
        AssociationQ[badCap], badCap["Reason"],
        True, None
    ];
    <|
        "SupportsSecondOrderQ" -> (badReason === None),
        "SupportsFastExplicitHessianQ" -> (badReason === None && allSpec),
        "Reason" -> badReason,
        "NumMembers" -> Length[members]
    |>
];

evaluateEnsembleWithCache[compEns_Association, input_] := Module[
    {nt, mCtrl, basisAssoc, isIdent, freezeIdx, dtVec, parsed, pMat, uMat, phiVec, jointFlag,
     caps, needCache, members, weights, nM, aggKind, beta, kernels, chunks, cid, workerOut,
     memberResults, fids, runVals, grads, leakVals, mStates,
     memberScores, sMin, fMin, expW, effWeights, aggScore, aggTermFid, aggRunVal, softMinFid,
     gradU, gradUFid, lamP, lamS, pen, gRegU, penNS, gU, gParam, gParamFid, gRegParam, cacheState},
    parsed = ($OCSolverRegistry["UnpackControlMatrix"])[compEns, input];
    If[parsed === $Failed, Return[$Failed]];
    {pMat, uMat, dtVec, phiVec, jointFlag} = parsed;
    nt = compEns["NumSlices"];
    mCtrl = compEns["NumControls"];
    basisAssoc = compEns["ControlBasis"];
    isIdent = basisAssoc["IsIdentityQ"];
    freezeIdx = Lookup[compEns, "FreezeControlSlices", {}];
    If[freezeIdx =!= {} && !isIdent, Return[$Failed]];
    If[freezeIdx =!= {},
        uMat = ToPackedRealMatrix[uMat];
        uMat[[freezeIdx, All]] = 0.
    ];
    caps = checkEnsembleSecondOrderCapabilities[compEns];
    needCache = TrueQ[caps["SupportsSecondOrderQ"]];
    If[!needCache,
        With[{baseRes = evaluateCompiledEnsemble[compEns, input]},
            If[!AssociationQ[baseRes], Return[$Failed]];
            Return[Join[baseRes, <|
                "SupportsSecondOrderQ" -> False,
                "UnsupportedSecondOrderReason" -> caps["Reason"],
                "EvalState" -> <|"SupportsSecondOrderQ" -> False, "Reason" -> caps["Reason"]|>
            |>]]
        ]
    ];

    members = compEns["Members"];
    weights = compEns["Weights"];
    nM = compEns["NumMembers"];
    aggKind = compEns["Aggregation"];
    beta = compEns["SoftMinBeta"];

    If[TrueQ[compEns["UseParallel"]],
        kernels = compEns["WorkerKernels"];
        chunks = compEns["WorkerChunks"];
        cid = compEns["WorkerCacheId"];
        workerOut = With[{wId = cid, uIn = uMat},
            ParallelEvaluate[
                With[{mList = OptimalControl`Ensemble`Private`$WorkerMemberStore[wId]},
                    Module[{resList, stateList, outList},
                        resList = Map[evalMemberKernelWithCache[#, uIn] &, mList];
                        stateList = Lookup[#, "EvalState", None] & /@ resList;
                        OptimalControl`Ensemble`Private`$WorkerStateCache[wId] = stateList;
                        outList = Map[
                            With[{r = #, ext = Lookup[#, "Extra", None]},
                                {r["Fidelity"], Lookup[r, "RunningCost", 0.], r["Gradient"],
                                 If[AssociationQ[ext], Lookup[ext, "Leakage", None], None]}
                            ] &,
                            resList
                        ];
                        outList
                    ]
                ],
                kernels
            ]
        ];
        If[!ListQ[workerOut] || Length[workerOut] =!= Length[chunks], Return[$Failed]];
        fids = ConstantArray[0., nM];
        runVals = ConstantArray[0., nM];
        grads = ConstantArray[0., {nM, nt, mCtrl}];
        leakVals = ConstantArray[None, nM];
        mStates = None;
        Do[
            Do[
                fids[[chunks[[w, pos]]]] = workerOut[[w, pos, 1]];
                runVals[[chunks[[w, pos]]]] = workerOut[[w, pos, 2]];
                grads[[chunks[[w, pos]]]] = workerOut[[w, pos, 3]];
                leakVals[[chunks[[w, pos]]]] = workerOut[[w, pos, 4]],
                {pos, Length[chunks[[w]]]}
            ],
            {w, Length[chunks]}
        ],
        memberResults = Map[evalMemberKernelWithCache[#, uMat] &, members];
        If[!AllTrue[memberResults, AssociationQ], Return[$Failed]];
        fids = Developer`ToPackedArray[N[#["Fidelity"] & /@ memberResults]];
        runVals = Developer`ToPackedArray[N[Lookup[#, "RunningCost", 0.] & /@ memberResults]];
        grads = Developer`ToPackedArray[#["Gradient"] & /@ memberResults];
        mStates = #["EvalState"] & /@ memberResults;
        leakVals = Map[
            Function[mr, With[{ext = Lookup[mr, "Extra", None]},
                If[AssociationQ[ext], Lookup[ext, "Leakage", None], None]
            ]],
            memberResults
        ]
    ];

    memberScores = fids + runVals;
    If[aggKind === "SoftMin",
        sMin = Min[memberScores];
        expW = weights * Exp[-beta * (memberScores - sMin)];
        effWeights = Developer`ToPackedArray[N[expW / Total[expW]]];
        aggScore = sMin - (1. / beta) * Log[Total[expW]];
        fMin = Min[fids];
        aggTermFid = fMin - (1. / beta) * Log[Total[weights * Exp[-beta * (fids - fMin)]]];
        aggRunVal = aggScore - aggTermFid,
        effWeights = Developer`ToPackedArray[N[weights]];
        aggTermFid = N[weights . fids];
        aggRunVal = N[weights . runVals];
        aggScore = aggTermFid + aggRunVal
    ];
    gradU = Developer`ToPackedArray[effWeights . grads];

    lamP = compEns["PowerPenalty"];
    lamS = compEns["SmoothnessPenalty"];
    pen = ($OCSolverRegistry["ComputeSmoothPenaltyValue"])[compEns, uMat, dtVec];
    gRegU = ConstantArray[0., {nt, mCtrl}];
    If[lamP > 0., gRegU += (-2. * lamP) * (dtVec * uMat)];
    If[lamS > 0. && nt > 1,
        With[{diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]), du = Rest[uMat] - Most[uMat]},
            gRegU[[1 ;; nt - 1]] += (2. * lamS) * (du / diffDt);
            gRegU[[2 ;; nt]] -= (2. * lamS) * (du / diffDt)
        ]
    ];
    softMinFid = If[aggKind === "SoftMin", aggTermFid, N[weights . fids]];
    penNS = If[KeyExistsQ[compEns, "ProxContext"], ($OCSolverRegistry["EvalNonSmoothPenalty"])[compEns["ProxContext"], uMat], 0.];
    gradUFid = gradU;
    If[freezeIdx =!= {}, gradUFid[[freezeIdx, All]] = 0.];
    gParamFid = If[isIdent, gradUFid, Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gradUFid]];
    gU = Developer`ToPackedArray[gradU + gRegU];
    If[freezeIdx =!= {}, gU[[freezeIdx, All]] = 0.];
    gParam = If[isIdent, gU, Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gU]];
    gRegParam = Developer`ToPackedArray[gParam - gParamFid];

    cacheState = <|
        "Kernel" -> "Ensemble",
        "SupportsSecondOrderQ" -> True,
        "SupportsFastExplicitHessianQ" -> TrueQ[caps["SupportsFastExplicitHessianQ"]],
        "UseParallel" -> TrueQ[compEns["UseParallel"]],
        "WorkerCacheId" -> Lookup[compEns, "WorkerCacheId", None],
        "MemberStates" -> mStates,
        "MemberSliceGrads" -> grads,
        "AggSliceGrad" -> gradU,
        "EffWeights" -> effWeights,
        "Aggregation" -> aggKind,
        "SoftMinBeta" -> beta,
        "BasisAssoc" -> basisAssoc,
        "ControlsMat" -> uMat,
        "ParametersMat" -> pMat,
        "TimeSteps" -> dtVec
    |>;

    <|
        "Fidelity" -> aggTermFid,
        "TerminalFidelity" -> aggTermFid,
        "RunningCost" -> aggRunVal,
        "Leakage" -> If[AllTrue[leakVals, NumericQ], N[weights . leakVals], None],
        "SmoothObjective" -> aggScore - pen,
        "SmoothPenalty" -> pen,
        "NonSmoothPenalty" -> penNS,
        "Objective" -> aggScore - pen - penNS,
        "Penalty" -> pen + penNS,
        "Gradient" -> gParam,
        "ControlGradient" -> gParam,
        "SliceGradient" -> gU,
        "FidelitySliceGradient" -> gradUFid,
        "FidelityControlGradient" -> gParamFid,
        "SmoothPenaltyControlGradient" -> gRegParam,
        "Phases" -> None,
        "Parameters" -> pMat,
        "Controls" -> uMat,
        "TimeSteps" -> dtVec,
        "TotalDuration" -> Total[dtVec],
        "MemberFidelities" -> fids,
        "MemberRunningCosts" -> runVals,
        "MeanFidelity" -> N[weights . fids],
        "SoftMinFidelity" -> softMinFid,
        "MinFidelity" -> Min[fids],
        "MaxFidelity" -> Max[fids],
        "SupportsSecondOrderQ" -> True,
        "UnsupportedSecondOrderReason" -> None,
        "EvalState" -> cacheState
    |>
];

evaluateEnsembleHVP[compEns_Association, evalState_Association, dirInput_] := Module[
    {caps, basisAssoc, nParams, mCtrl, nt, freezeIdx, isVec, vParamMat, vSliceMat,
     members, nM, effWeights, aggKind, beta, mGrads, aggGrad, kernels, chunks, cid, workerOut,
     hMemSlices, hFidSlice, dots, meanDot, covTerm, hPenSlice, hTotSlice, hParamMat},
    caps = checkEnsembleSecondOrderCapabilities[compEns];
    If[!TrueQ[caps["SupportsSecondOrderQ"]] || !TrueQ[Lookup[evalState, "SupportsSecondOrderQ", False]],
        Return[$Failed]
    ];
    basisAssoc = compEns["ControlBasis"];
    nParams = compEns["NumParams"];
    mCtrl = compEns["NumControls"];
    nt = compEns["NumSlices"];
    freezeIdx = Lookup[compEns, "FreezeControlSlices", {}];
    isVec = VectorQ[dirInput, NumericQ];
    vParamMat = Which[
        isVec && Length[dirInput] === nParams * mCtrl,
            Partition[Developer`ToPackedArray[N[Re[dirInput]]], mCtrl],
        MatrixQ[dirInput, NumericQ] && Dimensions[dirInput] === {nParams, mCtrl},
            ToPackedRealMatrix[dirInput],
        True,
            Return[$Failed]
    ];
    vSliceMat = If[TrueQ[basisAssoc["IsIdentityQ"]],
        vParamMat,
        Developer`ToPackedArray[basisAssoc["BasisMatrix"] . vParamMat]
    ];
    If[freezeIdx =!= {}, vSliceMat[[freezeIdx, All]] = 0.];

    members = compEns["Members"];
    nM = compEns["NumMembers"];
    effWeights = evalState["EffWeights"];
    aggKind = evalState["Aggregation"];
    beta = evalState["SoftMinBeta"];
    mGrads = evalState["MemberSliceGrads"];
    aggGrad = evalState["AggSliceGrad"];

    If[TrueQ[evalState["UseParallel"]],
        kernels = compEns["WorkerKernels"];
        chunks = compEns["WorkerChunks"];
        cid = evalState["WorkerCacheId"];
        workerOut = With[{wId = cid, vIn = vSliceMat},
            ParallelEvaluate[
                MapThread[
                    evalMemberHVPKernel[#1, #2, vIn] &,
                    {OptimalControl`Ensemble`Private`$WorkerMemberStore[wId],
                     OptimalControl`Ensemble`Private`$WorkerStateCache[wId]}
                ],
                kernels
            ]
        ];
        If[!ListQ[workerOut] || Length[workerOut] =!= Length[chunks], Return[$Failed]];
        hMemSlices = ConstantArray[0., {nM, nt, mCtrl}];
        Do[
            Do[
                hMemSlices[[chunks[[w, pos]]]] = workerOut[[w, pos]],
                {pos, Length[chunks[[w]]]}
            ],
            {w, Length[chunks]}
        ],
        hMemSlices = MapThread[
            evalMemberHVPKernel[#1, #2, vSliceMat] &,
            {members, evalState["MemberStates"]}
        ];
        If[!AllTrue[hMemSlices, MatrixQ[#, NumericQ] &], Return[$Failed]];
        hMemSlices = Developer`ToPackedArray[hMemSlices]
    ];

    hFidSlice = Developer`ToPackedArray[effWeights . hMemSlices];
    If[aggKind === "SoftMin",
        dots = Map[Total[Flatten[# * vSliceMat]] &, mGrads];
        meanDot = Total[Flatten[aggGrad * vSliceMat]];
        covTerm = Sum[effWeights[[a]] * dots[[a]] * mGrads[[a]], {a, nM}] - meanDot * aggGrad;
        hFidSlice = Developer`ToPackedArray[hFidSlice - beta * covTerm]
    ];
    hPenSlice = ($OCSolverRegistry["ApplySmoothPenaltyHVP"])[compEns, vSliceMat, evalState["TimeSteps"]];
    hTotSlice = Developer`ToPackedArray[hFidSlice + hPenSlice];
    If[freezeIdx =!= {}, hTotSlice[[freezeIdx, All]] = 0.];
    hParamMat = If[TrueQ[basisAssoc["IsIdentityQ"]],
        hTotSlice,
        Developer`ToPackedArray[basisAssoc["BasisTranspose"] . hTotSlice]
    ];
    If[isVec, Developer`ToPackedArray[Flatten[hParamMat]], hParamMat]
];

evaluateEnsembleHessian[compEns_Association, evalState_Association, mode_: Automatic] := Module[
    {caps, basisAssoc, nParams, mCtrl, nt, nSliceTot, nDim, freezeIdx, dtVec,
     members, nM, effWeights, aggKind, beta, mGrads, aggGrad, kernels, chunks, cid, workerOut,
     hMemList, hSliceFid, gFlatMats, gSMFlat, covMat, hSlicePen, hSliceTot, i, j, k, idx, eVec, cols},
    caps = checkEnsembleSecondOrderCapabilities[compEns];
    If[!TrueQ[caps["SupportsSecondOrderQ"]] || !TrueQ[Lookup[evalState, "SupportsSecondOrderQ", False]],
        Return[$Failed]
    ];
    basisAssoc = compEns["ControlBasis"];
    nParams = compEns["NumParams"];
    mCtrl = compEns["NumControls"];
    nt = compEns["NumSlices"];
    nSliceTot = nt * mCtrl;
    nDim = nParams * mCtrl;
    freezeIdx = Lookup[compEns, "FreezeControlSlices", {}];
    dtVec = evalState["TimeSteps"];

    If[TrueQ[caps["SupportsFastExplicitHessianQ"]] && mode =!= "ColumnHVP",
        members = compEns["Members"];
        nM = compEns["NumMembers"];
        effWeights = evalState["EffWeights"];
        aggKind = evalState["Aggregation"];
        beta = evalState["SoftMinBeta"];
        mGrads = evalState["MemberSliceGrads"];
        aggGrad = evalState["AggSliceGrad"];
        If[TrueQ[evalState["UseParallel"]],
            kernels = compEns["WorkerKernels"];
            chunks = compEns["WorkerChunks"];
            cid = evalState["WorkerCacheId"];
            workerOut = With[{wId = cid},
                ParallelEvaluate[
                    MapThread[
                        evalMemberHessianKernel[#1, #2] &,
                        {OptimalControl`Ensemble`Private`$WorkerMemberStore[wId],
                         OptimalControl`Ensemble`Private`$WorkerStateCache[wId]}
                    ],
                    kernels
                ]
            ];
            If[!ListQ[workerOut] || Length[workerOut] =!= Length[chunks], Return[$Failed]];
            hMemList = ConstantArray[0., {nM, nSliceTot, nSliceTot}];
            Do[
                Do[
                    hMemList[[chunks[[w, pos]]]] = workerOut[[w, pos]],
                    {pos, Length[chunks[[w]]]}
                ],
                {w, Length[chunks]}
            ],
            hMemList = MapThread[
                evalMemberHessianKernel[#1, #2] &,
                {members, evalState["MemberStates"]}
            ];
            If[!AllTrue[hMemList, MatrixQ[#, NumericQ] &], Return[$Failed]]
        ];
        hSliceFid = Sum[effWeights[[i]] * hMemList[[i]], {i, nM}];
        If[aggKind === "SoftMin",
            gFlatMats = Developer`ToPackedArray[Flatten /@ mGrads];
            gSMFlat = Developer`ToPackedArray[Flatten[aggGrad]];
            covMat = Transpose[effWeights * gFlatMats] . gFlatMats - Outer[Times, gSMFlat, gSMFlat];
            hSliceFid = hSliceFid - beta * covMat
        ];
        hSlicePen = ($OCSolverRegistry["BuildSmoothPenaltyHessian"])[compEns, dtVec];
        hSliceTot = Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ (hSliceFid + hSlicePen)];
        If[freezeIdx =!= {},
            Do[
                Do[
                    idx = (j - 1) * mCtrl + k;
                    hSliceTot[[idx, All]] = 0.;
                    hSliceTot[[All, idx]] = 0.,
                    {k, mCtrl}
                ],
                {j, freezeIdx}
            ]
        ];
        Return[($OCSolverRegistry["PullbackSliceHessianToBasis"])[basisAssoc, hSliceTot, mCtrl]]
    ];

    eVec = ConstantArray[0., nDim];
    cols = Table[
        eVec[[i]] = 1.;
        With[{c = evaluateEnsembleHVP[compEns, evalState, eVec]},
            eVec[[i]] = 0.;
            If[!VectorQ[c, NumericQ], Return[$Failed, Module]];
            c
        ],
        {i, nDim}
    ];
    Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ Transpose[cols]]
];

(* Public EvaluateOptimalControl on OptimalControlEnsemble *)
EvaluateOptimalControl[ens_OptimalControlEnsemble, input_, returnProps_: False] :=
    With[{comp = compileEnsembleIfNeeded[ens]},
        If[comp === $Failed, $Failed, evaluateCompiledEnsemble[comp, input]]
    ];

(* Forward-only ensemble objective evaluation (no backward passes) *)
evaluateCompiledEnsembleObjective[compEns_Association, input_] := Module[
    {nt, mCtrl, basisAssoc, isIdent, hasDur, freezeIdx, dtVec, phiVec, parsed, pMat, uMat, jointFlag,
     members, weights, nM, aggKind, beta, buildSched, stepOverrides = None,
     kernels, chunks, cid, ovMap, workerOut, memberResults, fids, runVals, leakVals,
     memberScores, sMin, fMin, expW, aggScore, aggTermFid, aggRunVal, softMinFid,
     penSmooth, penNS, hasPhi},
    parsed = ($OCSolverRegistry["UnpackControlMatrix"])[compEns, input];
    If[parsed === $Failed, Return[$Failed]];
    {pMat, uMat, dtVec, phiVec, jointFlag} = parsed;
    nt = compEns["NumSlices"];
    mCtrl = compEns["NumControls"];
    basisAssoc = compEns["ControlBasis"];
    isIdent = basisAssoc["IsIdentityQ"];
    hasDur = TrueQ[Lookup[compEns, "HasDurationOptimization", False]];
    freezeIdx = Lookup[compEns, "FreezeControlSlices", {}];
    hasPhi = TrueQ[Lookup[compEns, "HasPhaseFreedom", False]];
    If[freezeIdx =!= {} && !isIdent, Return[$Failed]];
    If[freezeIdx =!= {},
        uMat = ToPackedRealMatrix[uMat];
        uMat[[freezeIdx, All]] = 0.
    ];
    If[hasDur && !AllTrue[dtVec, TrueQ[# > 0.] &], Return[$Failed]];

    members = compEns["Members"];
    weights = compEns["Weights"];
    nM = compEns["NumMembers"];
    aggKind = compEns["Aggregation"];
    beta = compEns["SoftMinBeta"];

    If[hasDur,
        buildSched = $OCSolverRegistry["BuildFineDriftSchedule"];
        stepOverrides = Table[
            Module[{mem = members[[mIdx]], t0M, sched},
                t0M = Lookup[mem, "InitialTimePoint", 0.];
                If[Lookup[mem, "DriftKind", "Static"] === "Continuous",
                    sched = buildSched[mem["DriftEvaluator"], dtVec, t0M, mem["DriftIntegrationSpec"], mem["Space"]];
                    If[sched === $Failed, Return[$Failed, Module]];
                    <|"Drift" -> sched["Drift"], "FineTimeSteps" -> sched["FineTimeSteps"]|>,
                    <|"Drift" -> None, "FineTimeSteps" -> dtVec|>
                ]
            ],
            {mIdx, nM}
        ];
        If[MemberQ[stepOverrides, $Failed], Return[$Failed]]
    ];

    If[TrueQ[compEns["UseParallel"]],
        kernels = compEns["WorkerKernels"];
        chunks = compEns["WorkerChunks"];
        cid = compEns["WorkerCacheId"];
        ovMap = If[hasDur,
            AssociationThread[compEns["WorkerKernelIDs"], Map[stepOverrides[[#]] &, chunks]],
            <||>
        ];
        workerOut = With[{wId = cid, uIn = uMat, ovById = ovMap, hd = hasDur, pv = If[hasPhi, phiVec, None]},
            ParallelEvaluate[
                With[{mList = OptimalControl`Ensemble`Private`$WorkerMemberStore[wId]},
                    With[{oList = If[hd, ovById[$KernelID], ConstantArray[None, Length[mList]]]},
                        MapThread[
                            With[{r = evalMemberObjectiveOnlyKernel[#1, uIn, #2, pv]},
                                With[{ext = Lookup[r, "Extra", None]},
                                    {r["Fidelity"], Lookup[r, "RunningCost", 0.],
                                     If[AssociationQ[ext], Lookup[ext, "Leakage", None], None]}
                                ]
                            ] &,
                            {mList, oList}
                        ]
                    ]
                ],
                kernels
            ]
        ];
        If[!ListQ[workerOut] || Length[workerOut] =!= Length[chunks], Return[$Failed]];
        fids = ConstantArray[0., nM];
        runVals = ConstantArray[0., nM];
        leakVals = ConstantArray[None, nM];
        Do[
            Do[
                fids[[chunks[[w, pos]]]] = workerOut[[w, pos, 1]];
                runVals[[chunks[[w, pos]]]] = workerOut[[w, pos, 2]];
                leakVals[[chunks[[w, pos]]]] = workerOut[[w, pos, 3]],
                {pos, Length[chunks[[w]]]}
            ],
            {w, Length[chunks]}
        ],
        memberResults = With[{pv = If[hasPhi, phiVec, None]},
            If[hasDur,
                MapThread[evalMemberObjectiveOnlyKernel[#1, uMat, #2, pv] &, {members, stepOverrides}],
                Map[evalMemberObjectiveOnlyKernel[#, uMat, None, pv] &, members]
            ]
        ];
        If[!AllTrue[memberResults, AssociationQ], Return[$Failed]];
        fids = Developer`ToPackedArray[N[#["Fidelity"] & /@ memberResults]];
        runVals = Developer`ToPackedArray[N[Lookup[#, "RunningCost", 0.] & /@ memberResults]];
        leakVals = Map[
            Function[mr, With[{ext = Lookup[mr, "Extra", None]},
                If[AssociationQ[ext], Lookup[ext, "Leakage", None], None]
            ]],
            memberResults
        ]
    ];

    memberScores = fids + runVals;
    If[aggKind === "SoftMin",
        sMin = Min[memberScores];
        expW = weights * Exp[-beta * (memberScores - sMin)];
        aggScore = sMin - (1. / beta) * Log[Total[expW]];
        fMin = Min[fids];
        aggTermFid = fMin - (1. / beta) * Log[Total[weights * Exp[-beta * (fids - fMin)]]];
        aggRunVal = aggScore - aggTermFid,
        aggTermFid = N[weights . fids];
        aggRunVal = N[weights . runVals];
        aggScore = aggTermFid + aggRunVal
    ];

    penSmooth = ($OCSolverRegistry["ComputeSmoothPenaltyValue"])[compEns, uMat, dtVec];
    softMinFid = If[aggKind === "SoftMin", aggTermFid, N[weights . fids]];
    With[{penNS = If[KeyExistsQ[compEns, "ProxContext"], ($OCSolverRegistry["EvalNonSmoothPenalty"])[compEns["ProxContext"], uMat], 0.]},
        <|
            "Fidelity" -> aggTermFid,
            "TerminalFidelity" -> aggTermFid,
            "RunningCost" -> aggRunVal,
            "Leakage" -> If[AllTrue[leakVals, NumericQ], N[weights . leakVals], None],
            "SmoothObjective" -> aggScore - penSmooth,
            "SmoothPenalty" -> penSmooth,
            "NonSmoothPenalty" -> penNS,
            "Objective" -> aggScore - penSmooth - penNS,
            "Penalty" -> penSmooth + penNS,
            "Phases" -> If[hasPhi, phiVec, None],
            "Parameters" -> pMat,
            "Controls" -> uMat,
            "TimeSteps" -> dtVec,
            "TotalDuration" -> Total[dtVec],
            "MemberFidelities" -> fids,
            "MemberRunningCosts" -> runVals,
            "MeanFidelity" -> N[weights . fids],
            "SoftMinFidelity" -> softMinFid,
            "MinFidelity" -> Min[fids],
            "MaxFidelity" -> Max[fids]
        |>
    ]
];

(* Phase C Ensemble Second-Order Capability, Cached Evaluator, HVP, and Explicit Hessian *)
checkEnsembleSecondOrderCapabilities[compEns_Association] := Module[
    {hasDur, hasPhi, hasNS, members, memberCaps, badCap, badReason, allSpec},
    hasDur = TrueQ[Lookup[compEns, "HasDurationOptimization", False]];
    hasPhi = TrueQ[Lookup[compEns, "HasPhaseFreedom", False]];
    hasNS = TrueQ[Lookup[compEns, "L1Penalty", 0.] > 0. || Lookup[compEns, "TVPenalty", 0.] > 0.];
    members = Lookup[compEns, "Members", {}];
    memberCaps = Table[
        With[{m = members[[i]],
              hasRC = AssociationQ[Lookup[members[[i]], "RunningCostSpec", None]],
              method = Lookup[members[[i]], "Method", "Unknown"],
              objOK = TrueQ[Lookup[members[[i]], "Objective", <||>] // (Lookup[#, "ExactHessianSupportedQ", False] &)]},
            Which[
                hasRC, <|"OK" -> False, "Spec" -> False, "Reason" -> "Second-order derivatives are not yet enabled for RunningCost / StatePenalty."|>,
                !MemberQ[{"Spectral", "StateVector", "Liouville"}, method],
                    <|"OK" -> False, "Spec" -> False, "Reason" -> "Second-order derivatives require Method -> \"Spectral\", \"StateVector\", or \"Liouville\" (got " <> ToString[method] <> ")."|>,
                !objOK,
                    <|"OK" -> False, "Spec" -> False, "Reason" -> Lookup[members[[i]]["Objective"], "UnsupportedSecondOrderReason", "Ensemble member objective does not support exact second derivatives."]|>,
                True,
                    <|"OK" -> True, "Spec" -> (method === "Spectral"), "Reason" -> None|>
            ]
        ],
        {i, Length[members]}
    ];
    badCap = SelectFirst[memberCaps, !TrueQ[#["OK"]] &, None];
    allSpec = AllTrue[memberCaps, TrueQ[#["Spec"]] &];
    badReason = Which[
        hasDur, "Second-order derivatives are not yet enabled for OptimizeDurations -> True.",
        hasPhi, "Second-order derivatives are not yet enabled for parametric PhaseFreedom.",
        hasNS, "Nonsmooth L1Penalty / TVPenalty problems do not have a smooth second-order Hessian.",
        AssociationQ[badCap], badCap["Reason"],
        True, None
    ];
    <|
        "SupportsSecondOrderQ" -> (badReason === None),
        "SupportsFastExplicitHessianQ" -> (badReason === None && allSpec),
        "Reason" -> badReason,
        "NumMembers" -> Length[members]
    |>
];

evaluateEnsembleWithCache[compEns_Association, input_] := Module[
    {nt, mCtrl, basisAssoc, isIdent, freezeIdx, dtVec, parsed, pMat, uMat, phiVec, jointFlag,
     caps, needCache, members, weights, nM, aggKind, beta, kernels, chunks, cid, workerOut,
     memberResults, fids, runVals, grads, leakVals, mStates,
     memberScores, sMin, fMin, expW, effWeights, aggScore, aggTermFid, aggRunVal, softMinFid,
     gradU, gradUFid, lamP, lamS, pen, gRegU, penNS, gU, gParam, gParamFid, gRegParam, cacheState},
    parsed = ($OCSolverRegistry["UnpackControlMatrix"])[compEns, input];
    If[parsed === $Failed, Return[$Failed]];
    {pMat, uMat, dtVec, phiVec, jointFlag} = parsed;
    nt = compEns["NumSlices"];
    mCtrl = compEns["NumControls"];
    basisAssoc = compEns["ControlBasis"];
    isIdent = basisAssoc["IsIdentityQ"];
    freezeIdx = Lookup[compEns, "FreezeControlSlices", {}];
    If[freezeIdx =!= {} && !isIdent, Return[$Failed]];
    If[freezeIdx =!= {},
        uMat = ToPackedRealMatrix[uMat];
        uMat[[freezeIdx, All]] = 0.
    ];
    caps = checkEnsembleSecondOrderCapabilities[compEns];
    needCache = TrueQ[caps["SupportsSecondOrderQ"]];
    If[!needCache,
        With[{baseRes = evaluateCompiledEnsemble[compEns, input]},
            If[!AssociationQ[baseRes], Return[$Failed]];
            Return[Join[baseRes, <|
                "SupportsSecondOrderQ" -> False,
                "UnsupportedSecondOrderReason" -> caps["Reason"],
                "EvalState" -> <|"SupportsSecondOrderQ" -> False, "Reason" -> caps["Reason"]|>
            |>]]
        ]
    ];

    members = compEns["Members"];
    weights = compEns["Weights"];
    nM = compEns["NumMembers"];
    aggKind = compEns["Aggregation"];
    beta = compEns["SoftMinBeta"];

    If[TrueQ[compEns["UseParallel"]],
        kernels = compEns["WorkerKernels"];
        chunks = compEns["WorkerChunks"];
        cid = compEns["WorkerCacheId"];
        workerOut = With[{wId = cid, uIn = uMat},
            ParallelEvaluate[
                With[{mList = OptimalControl`Ensemble`Private`$WorkerMemberStore[wId]},
                    Module[{resList, stateList, outList},
                        resList = Map[evalMemberKernelWithCache[#, uIn] &, mList];
                        stateList = Lookup[#, "EvalState", None] & /@ resList;
                        OptimalControl`Ensemble`Private`$WorkerStateCache[wId] = stateList;
                        outList = Map[
                            With[{r = #, ext = Lookup[#, "Extra", None]},
                                {r["Fidelity"], Lookup[r, "RunningCost", 0.], r["Gradient"],
                                 If[AssociationQ[ext], Lookup[ext, "Leakage", None], None]}
                            ] &,
                            resList
                        ];
                        outList
                    ]
                ],
                kernels
            ]
        ];
        If[!ListQ[workerOut] || Length[workerOut] =!= Length[chunks], Return[$Failed]];
        fids = ConstantArray[0., nM];
        runVals = ConstantArray[0., nM];
        grads = ConstantArray[0., {nM, nt, mCtrl}];
        leakVals = ConstantArray[None, nM];
        mStates = None;
        Do[
            Do[
                fids[[chunks[[w, pos]]]] = workerOut[[w, pos, 1]];
                runVals[[chunks[[w, pos]]]] = workerOut[[w, pos, 2]];
                grads[[chunks[[w, pos]]]] = workerOut[[w, pos, 3]];
                leakVals[[chunks[[w, pos]]]] = workerOut[[w, pos, 4]],
                {pos, Length[chunks[[w]]]}
            ],
            {w, Length[chunks]}
        ],
        memberResults = Map[evalMemberKernelWithCache[#, uMat] &, members];
        If[!AllTrue[memberResults, AssociationQ], Return[$Failed]];
        fids = Developer`ToPackedArray[N[#["Fidelity"] & /@ memberResults]];
        runVals = Developer`ToPackedArray[N[Lookup[#, "RunningCost", 0.] & /@ memberResults]];
        grads = Developer`ToPackedArray[#["Gradient"] & /@ memberResults];
        mStates = #["EvalState"] & /@ memberResults;
        leakVals = Map[
            Function[mr, With[{ext = Lookup[mr, "Extra", None]},
                If[AssociationQ[ext], Lookup[ext, "Leakage", None], None]
            ]],
            memberResults
        ]
    ];

    memberScores = fids + runVals;
    If[aggKind === "SoftMin",
        sMin = Min[memberScores];
        expW = weights * Exp[-beta * (memberScores - sMin)];
        effWeights = Developer`ToPackedArray[N[expW / Total[expW]]];
        aggScore = sMin - (1. / beta) * Log[Total[expW]];
        fMin = Min[fids];
        aggTermFid = fMin - (1. / beta) * Log[Total[weights * Exp[-beta * (fids - fMin)]]];
        aggRunVal = aggScore - aggTermFid,
        effWeights = Developer`ToPackedArray[N[weights]];
        aggTermFid = N[weights . fids];
        aggRunVal = N[weights . runVals];
        aggScore = aggTermFid + aggRunVal
    ];
    gradU = Developer`ToPackedArray[effWeights . grads];

    lamP = compEns["PowerPenalty"];
    lamS = compEns["SmoothnessPenalty"];
    pen = ($OCSolverRegistry["ComputeSmoothPenaltyValue"])[compEns, uMat, dtVec];
    gRegU = ConstantArray[0., {nt, mCtrl}];
    If[lamP > 0., gRegU += (-2. * lamP) * (dtVec * uMat)];
    If[lamS > 0. && nt > 1,
        With[{diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]), du = Rest[uMat] - Most[uMat]},
            gRegU[[1 ;; nt - 1]] += (2. * lamS) * (du / diffDt);
            gRegU[[2 ;; nt]] -= (2. * lamS) * (du / diffDt)
        ]
    ];
    softMinFid = If[aggKind === "SoftMin", aggTermFid, N[weights . fids]];
    penNS = If[KeyExistsQ[compEns, "ProxContext"], ($OCSolverRegistry["EvalNonSmoothPenalty"])[compEns["ProxContext"], uMat], 0.];
    gradUFid = gradU;
    If[freezeIdx =!= {}, gradUFid[[freezeIdx, All]] = 0.];
    gParamFid = If[isIdent, gradUFid, Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gradUFid]];
    gU = Developer`ToPackedArray[gradU + gRegU];
    If[freezeIdx =!= {}, gU[[freezeIdx, All]] = 0.];
    gParam = If[isIdent, gU, Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gU]];
    gRegParam = Developer`ToPackedArray[gParam - gParamFid];

    cacheState = <|
        "Kernel" -> "Ensemble",
        "SupportsSecondOrderQ" -> True,
        "SupportsFastExplicitHessianQ" -> TrueQ[caps["SupportsFastExplicitHessianQ"]],
        "UseParallel" -> TrueQ[compEns["UseParallel"]],
        "WorkerCacheId" -> Lookup[compEns, "WorkerCacheId", None],
        "MemberStates" -> mStates,
        "MemberSliceGrads" -> grads,
        "AggSliceGrad" -> gradU,
        "EffWeights" -> effWeights,
        "Aggregation" -> aggKind,
        "SoftMinBeta" -> beta,
        "BasisAssoc" -> basisAssoc,
        "ControlsMat" -> uMat,
        "ParametersMat" -> pMat,
        "TimeSteps" -> dtVec
    |>;

    <|
        "Fidelity" -> aggTermFid,
        "TerminalFidelity" -> aggTermFid,
        "RunningCost" -> aggRunVal,
        "Leakage" -> If[AllTrue[leakVals, NumericQ], N[weights . leakVals], None],
        "SmoothObjective" -> aggScore - pen,
        "SmoothPenalty" -> pen,
        "NonSmoothPenalty" -> penNS,
        "Objective" -> aggScore - pen - penNS,
        "Penalty" -> pen + penNS,
        "Gradient" -> gParam,
        "ControlGradient" -> gParam,
        "SliceGradient" -> gU,
        "FidelitySliceGradient" -> gradUFid,
        "FidelityControlGradient" -> gParamFid,
        "SmoothPenaltyControlGradient" -> gRegParam,
        "Phases" -> None,
        "Parameters" -> pMat,
        "Controls" -> uMat,
        "TimeSteps" -> dtVec,
        "TotalDuration" -> Total[dtVec],
        "MemberFidelities" -> fids,
        "MemberRunningCosts" -> runVals,
        "MeanFidelity" -> N[weights . fids],
        "SoftMinFidelity" -> softMinFid,
        "MinFidelity" -> Min[fids],
        "MaxFidelity" -> Max[fids],
        "SupportsSecondOrderQ" -> True,
        "UnsupportedSecondOrderReason" -> None,
        "EvalState" -> cacheState
    |>
];

evaluateEnsembleHVP[compEns_Association, evalState_Association, dirInput_] := Module[
    {caps, basisAssoc, nParams, mCtrl, nt, freezeIdx, isVec, vParamMat, vSliceMat,
     members, nM, effWeights, aggKind, beta, mGrads, aggGrad, kernels, chunks, cid, workerOut,
     hMemSlices, hFidSlice, dots, meanDot, covTerm, hPenSlice, hTotSlice, hParamMat},
    caps = checkEnsembleSecondOrderCapabilities[compEns];
    If[!TrueQ[caps["SupportsSecondOrderQ"]] || !TrueQ[Lookup[evalState, "SupportsSecondOrderQ", False]],
        Return[$Failed]
    ];
    basisAssoc = compEns["ControlBasis"];
    nParams = compEns["NumParams"];
    mCtrl = compEns["NumControls"];
    nt = compEns["NumSlices"];
    freezeIdx = Lookup[compEns, "FreezeControlSlices", {}];
    isVec = VectorQ[dirInput, NumericQ];
    vParamMat = Which[
        isVec && Length[dirInput] === nParams * mCtrl,
            Partition[Developer`ToPackedArray[N[Re[dirInput]]], mCtrl],
        MatrixQ[dirInput, NumericQ] && Dimensions[dirInput] === {nParams, mCtrl},
            ToPackedRealMatrix[dirInput],
        True,
            Return[$Failed]
    ];
    vSliceMat = If[TrueQ[basisAssoc["IsIdentityQ"]],
        vParamMat,
        Developer`ToPackedArray[basisAssoc["BasisMatrix"] . vParamMat]
    ];
    If[freezeIdx =!= {}, vSliceMat[[freezeIdx, All]] = 0.];

    members = compEns["Members"];
    nM = compEns["NumMembers"];
    effWeights = evalState["EffWeights"];
    aggKind = evalState["Aggregation"];
    beta = evalState["SoftMinBeta"];
    mGrads = evalState["MemberSliceGrads"];
    aggGrad = evalState["AggSliceGrad"];

    If[TrueQ[evalState["UseParallel"]],
        kernels = compEns["WorkerKernels"];
        chunks = compEns["WorkerChunks"];
        cid = evalState["WorkerCacheId"];
        workerOut = With[{wId = cid, vIn = vSliceMat},
            ParallelEvaluate[
                MapThread[
                    evalMemberHVPKernel[#1, #2, vIn] &,
                    {OptimalControl`Ensemble`Private`$WorkerMemberStore[wId],
                     OptimalControl`Ensemble`Private`$WorkerStateCache[wId]}
                ],
                kernels
            ]
        ];
        If[!ListQ[workerOut] || Length[workerOut] =!= Length[chunks], Return[$Failed]];
        hMemSlices = ConstantArray[0., {nM, nt, mCtrl}];
        Do[
            Do[
                hMemSlices[[chunks[[w, pos]]]] = workerOut[[w, pos]],
                {pos, Length[chunks[[w]]]}
            ],
            {w, Length[chunks]}
        ],
        hMemSlices = MapThread[
            evalMemberHVPKernel[#1, #2, vSliceMat] &,
            {members, evalState["MemberStates"]}
        ];
        If[!AllTrue[hMemSlices, MatrixQ[#, NumericQ] &], Return[$Failed]];
        hMemSlices = Developer`ToPackedArray[hMemSlices]
    ];

    hFidSlice = Developer`ToPackedArray[effWeights . hMemSlices];
    If[aggKind === "SoftMin",
        dots = Map[Total[Flatten[# * vSliceMat]] &, mGrads];
        meanDot = Total[Flatten[aggGrad * vSliceMat]];
        covTerm = Sum[effWeights[[a]] * dots[[a]] * mGrads[[a]], {a, nM}] - meanDot * aggGrad;
        hFidSlice = Developer`ToPackedArray[hFidSlice - beta * covTerm]
    ];
    hPenSlice = ($OCSolverRegistry["ApplySmoothPenaltyHVP"])[compEns, vSliceMat, evalState["TimeSteps"]];
    hTotSlice = Developer`ToPackedArray[hFidSlice + hPenSlice];
    If[freezeIdx =!= {}, hTotSlice[[freezeIdx, All]] = 0.];
    hParamMat = If[TrueQ[basisAssoc["IsIdentityQ"]],
        hTotSlice,
        Developer`ToPackedArray[basisAssoc["BasisTranspose"] . hTotSlice]
    ];
    If[isVec, Developer`ToPackedArray[Flatten[hParamMat]], hParamMat]
];

evaluateEnsembleHessian[compEns_Association, evalState_Association, mode_: Automatic] := Module[
    {caps, basisAssoc, nParams, mCtrl, nt, nSliceTot, nDim, freezeIdx, dtVec,
     members, nM, effWeights, aggKind, beta, mGrads, aggGrad, kernels, chunks, cid, workerOut,
     hMemList, hSliceFid, gFlatMats, gSMFlat, covMat, hSlicePen, hSliceTot, i, j, k, idx, eVec, cols},
    caps = checkEnsembleSecondOrderCapabilities[compEns];
    If[!TrueQ[caps["SupportsSecondOrderQ"]] || !TrueQ[Lookup[evalState, "SupportsSecondOrderQ", False]],
        Return[$Failed]
    ];
    basisAssoc = compEns["ControlBasis"];
    nParams = compEns["NumParams"];
    mCtrl = compEns["NumControls"];
    nt = compEns["NumSlices"];
    nSliceTot = nt * mCtrl;
    nDim = nParams * mCtrl;
    freezeIdx = Lookup[compEns, "FreezeControlSlices", {}];
    dtVec = evalState["TimeSteps"];

    If[TrueQ[caps["SupportsFastExplicitHessianQ"]] && mode =!= "ColumnHVP",
        members = compEns["Members"];
        nM = compEns["NumMembers"];
        effWeights = evalState["EffWeights"];
        aggKind = evalState["Aggregation"];
        beta = evalState["SoftMinBeta"];
        mGrads = evalState["MemberSliceGrads"];
        aggGrad = evalState["AggSliceGrad"];
        If[TrueQ[evalState["UseParallel"]],
            kernels = compEns["WorkerKernels"];
            chunks = compEns["WorkerChunks"];
            cid = evalState["WorkerCacheId"];
            workerOut = With[{wId = cid},
                ParallelEvaluate[
                    MapThread[
                        evalMemberHessianKernel[#1, #2] &,
                        {OptimalControl`Ensemble`Private`$WorkerMemberStore[wId],
                         OptimalControl`Ensemble`Private`$WorkerStateCache[wId]}
                    ],
                    kernels
                ]
            ];
            If[!ListQ[workerOut] || Length[workerOut] =!= Length[chunks], Return[$Failed]];
            hMemList = ConstantArray[0., {nM, nSliceTot, nSliceTot}];
            Do[
                Do[
                    hMemList[[chunks[[w, pos]]]] = workerOut[[w, pos]],
                    {pos, Length[chunks[[w]]]}
                ],
                {w, Length[chunks]}
            ],
            hMemList = MapThread[
                evalMemberHessianKernel[#1, #2] &,
                {members, evalState["MemberStates"]}
            ];
            If[!AllTrue[hMemList, MatrixQ[#, NumericQ] &], Return[$Failed]]
        ];
        hSliceFid = Sum[effWeights[[i]] * hMemList[[i]], {i, nM}];
        If[aggKind === "SoftMin",
            gFlatMats = Developer`ToPackedArray[Flatten /@ mGrads];
            gSMFlat = Developer`ToPackedArray[Flatten[aggGrad]];
            covMat = Transpose[effWeights * gFlatMats] . gFlatMats - Outer[Times, gSMFlat, gSMFlat];
            hSliceFid = hSliceFid - beta * covMat
        ];
        hSlicePen = ($OCSolverRegistry["BuildSmoothPenaltyHessian"])[compEns, dtVec];
        hSliceTot = Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ (hSliceFid + hSlicePen)];
        If[freezeIdx =!= {},
            Do[
                Do[
                    idx = (j - 1) * mCtrl + k;
                    hSliceTot[[idx, All]] = 0.;
                    hSliceTot[[All, idx]] = 0.,
                    {k, mCtrl}
                ],
                {j, freezeIdx}
            ]
        ];
        Return[($OCSolverRegistry["PullbackSliceHessianToBasis"])[basisAssoc, hSliceTot, mCtrl]]
    ];

    eVec = ConstantArray[0., nDim];
    cols = Table[
        eVec[[i]] = 1.;
        With[{c = evaluateEnsembleHVP[compEns, evalState, eVec]},
            eVec[[i]] = 0.;
            If[!VectorQ[c, NumericQ], Return[$Failed, Module]];
            c
        ],
        {i, nDim}
    ];
    Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ Transpose[cols]]
];

parseSoftMinSchedule[schedOpt_, aggOpt_] := Module[{betaFinal, raw},
    Which[
        schedOpt === None || schedOpt === False, {},
        schedOpt === Automatic,
            If[MatchQ[aggOpt, {"SoftMin", _?NumericQ}] && TrueQ[aggOpt[[2]] > 0],
                betaFinal = N[aggOpt[[2]]];
                If[betaFinal > 4.,
                    Developer`ToPackedArray[N[{2., Sqrt[2. * betaFinal], betaFinal}]],
                    {betaFinal}
                ],
                {}
            ],
        ListQ[schedOpt] && Length[schedOpt] >= 1 && AllTrue[schedOpt, NumericQ[#] && TrueQ[# >= 0] &],
            raw = N[schedOpt];
            If[MatchQ[aggOpt, {"SoftMin", _?NumericQ}] && TrueQ[aggOpt[[2]] > 0],
                betaFinal = N[aggOpt[[2]]];
                If[Abs[Last[raw] - betaFinal] > 10^-9, AppendTo[raw, betaFinal]]
            ];
            N[raw],
        True,
            Message[OptimalControlEnsemble::badschedule, schedOpt]; $Failed
    ]
];

extractEnsembleTiming[ens:OptimalControlEnsemble[items_List, ___], defaultNt_Integer] := Module[{firstProb, resolveTim},
    firstProb = First[items];
    If[ListQ[firstProb], firstProb = firstProb[[2]]];
    resolveTim = $OCSolverRegistry["ResolveTimingVector"];
    resolveTim[firstProb[[4]], defaultNt]
];
extractEnsembleTiming[ens:OptimalControlEnsemble[_, _, _, timing_, ___], defaultNt_Integer] :=
    ($OCSolverRegistry["ResolveTimingVector"])[timing, defaultNt];

retargetEnsembleTiming[OptimalControlEnsemble[items_List, opts___Rule], timingSpec_] := Module[{newItems},
    newItems = Map[
        Function[item,
            If[MatchQ[item, {_?NumericQ, _OptimalControlProblem}],
                {item[[1]], OptimalControlProblem[item[[2, 1]], item[[2, 2]], item[[2, 3]], timingSpec,
                    Sequence @@ (List @@ item[[2, 5 ;;]])]},
                OptimalControlProblem[item[[1]], item[[2]], item[[3]], timingSpec,
                    Sequence @@ (List @@ item[[5 ;;]])]
            ]
        ],
        items
    ];
    OptimalControlEnsemble[newItems, opts]
];
retargetEnsembleTiming[OptimalControlEnsemble[drift_, controls_, target_, _, opts___Rule], timingSpec_] :=
    OptimalControlEnsemble[drift, controls, target, timingSpec, opts];
retargetEnsembleTiming[ens_OptimalControlEnsemble, totalDur_?NumericQ, ntNew_Integer] :=
    retargetEnsembleTiming[ens, {totalDur, ntNew}];

assembleEnsembleResult[comp_Association, optRes_Association, extraFields_Association: <||>] := Module[
    {finalEval, maxViol},
    finalEval = evaluateCompiledEnsemble[comp, optRes["ParametersVec"]];
    If[finalEval === $Failed, Return[$Failed]];
    maxViol = ($OCSolverRegistry["ComputeConstraintViolation"])[comp["ProxContext"], finalEval["Controls"]];
    OptimalControlResult[Join[<|
        "Kind" -> "Ensemble",
        "Controls" -> finalEval["Controls"],
        "Parameters" -> finalEval["Parameters"],
        "FinalFidelity" -> finalEval["Fidelity"],
        "TerminalFidelity" -> Lookup[finalEval, "TerminalFidelity", finalEval["Fidelity"]],
        "RunningCost" -> Lookup[finalEval, "RunningCost", 0.],
        "Leakage" -> Lookup[finalEval, "Leakage", None],
        "OptimalPhases" -> ($OCSolverRegistry["ExtractOptimalPhases"])[Lookup[comp, "PhaseFreedomSpec", None], Lookup[finalEval, "Phases", None]],
        "FinalObjective" -> finalEval["Objective"],
        "FinalPenalty" -> finalEval["Penalty"],
        "SmoothPenalty" -> finalEval["SmoothPenalty"],
        "NonSmoothPenalty" -> finalEval["NonSmoothPenalty"],
        "MaxConstraintViolation" -> maxViol,
        "Optimizer" -> Lookup[optRes, "Optimizer", "LBFGS"],
        "MemberFidelities" -> finalEval["MemberFidelities"],
        "MemberRunningCosts" -> Lookup[finalEval, "MemberRunningCosts", ConstantArray[0., comp["NumMembers"]]],
        "MemberWeights" -> comp["Weights"],
        "MeanMemberFidelity" -> finalEval["MeanFidelity"],
        "SoftMinMemberFidelity" -> Lookup[finalEval, "SoftMinFidelity", finalEval["MeanFidelity"]],
        "MinMemberFidelity" -> finalEval["MinFidelity"],
        "MaxMemberFidelity" -> finalEval["MaxFidelity"],
        "ParameterPoints" -> comp["ParameterPoints"],
        "GridShape" -> comp["GridShape"],
        "Gradient" -> finalEval["ControlGradient"],
        "ControlGradient" -> finalEval["ControlGradient"],
        "DurationGradient" -> Lookup[finalEval, "DurationGradient", None],
        "ActiveDurationGradient" -> Lookup[finalEval, "ActiveDurationGradient", None],
        "TotalDurationGradient" -> Lookup[finalEval, "TotalDurationGradient", None],
        "DurationMode" -> Lookup[comp, "DurationMode", "None"],
        "PhaseGradient" -> Lookup[finalEval, "PhaseGradient", None],
        "GradientNorm" -> Norm[Flatten[finalEval["Gradient"]]],
        "ScaledGradientNorm" -> optRes["GradientNorm"],
        "Iterations" -> optRes["Iterations"],
        "Evaluations" -> optRes["Evaluations"],
        "ObjectiveEvaluations" -> optRes["ObjectiveEvaluations"],
        "GradientEvaluations" -> optRes["GradientEvaluations"],
        "LineSearchMethod" -> optRes["LineSearchMethod"],
        "StrongWolfeSatisfiedCount" -> optRes["StrongWolfeSatisfiedCount"],
        "LastStepLength" -> optRes["LastStepLength"],
        "ElapsedTime" -> optRes["ElapsedTime"],
        "Status" -> optRes["Status"],
        "ConvergenceHistory" -> <|
            "Objective" -> optRes["ObjectiveHistory"],
            "Fidelity" -> optRes["FidelityHistory"],
            "GradientNorm" -> optRes["GradientNormHistory"]
        |>,
        "FinalPropagator" -> None,
        "FinalState" -> None,
        "CompiledEnsemble" -> comp,
        "CompiledProblem" -> First[comp["Members"]],
        "TimeSteps" -> Lookup[finalEval, "TimeSteps", comp["TimeSteps"]],
        "TotalDuration" -> Lookup[finalEval, "TotalDuration", comp["TotalDuration"]],
        "InitialTimeSteps" -> comp["TimeSteps"],
        "OptimizeDurations" -> Lookup[comp, "OptimizeDurations", False],
        "FullParameters" -> optRes["ParametersVec"],
        "Space" -> comp["Space"],
        "Method" -> comp["Method"],
        "Basis" -> comp["Basis"],
        "OperatorBasis" -> comp["OperatorBasis"],
        "SymbolicDrift" -> comp["SymbolicDrift"],
        "SymbolicControls" -> comp["SymbolicControls"],
        "DriftUnitFactor" -> comp["DriftUnitFactor"],
        "ControlUnitFactor" -> comp["ControlUnitFactor"],
        "FrequencyUnits" -> comp["FrequencyUnits"],
        "ControlUnits" -> comp["ControlUnits"],
        "InitialTimePoint" -> Lookup[comp, "InitialTimePoint", 0.],
        "DriftIntegration" -> Lookup[comp, "DriftIntegration", "Static"],
        "IsTimeDependentDrift" -> Lookup[comp, "IsTimeDependentDrift", False]
    |>, KeyTake[optRes, {
        "FISTARestarts", "LipschitzEstimate", "DynamicsTime", "OptimizerTime",
        "AcceptedSteps", "RejectedSteps", "HVPEvaluations", "HessianEvaluations",
        "SkippedCurvaturePairs", "ApproxWolfeCount", "ExperimentalConstrainedQ",
        "CompactLBFGSFinalSigma",
        "CGIterations", "NegativeCurvatureCount", "HybridSwitchIteration",
        "HybridSwitchReason", "TrustRegionFinalRadius", "ARCFinalSigma"
    }], extraFields]]
];

optimizeEnsembleCompiledStart[comp_Association, optsList_List, initOverride_: None, seedOverride_: None] := Module[
    {initOpt, initScale, seedOpt, buildInit, dispatchOpt, x0Vec, optRes},
    {initOpt, initScale, seedOpt} = {InitialControls, InitialScale, RandomSeed} /. optsList;
    If[initOverride =!= None, initOpt = initOverride];
    If[seedOverride =!= None, seedOpt = seedOverride];
    buildInit = $OCSolverRegistry["BuildInitialVector"];
    dispatchOpt = $OCSolverRegistry["DispatchOptimizer"];
    x0Vec = buildInit[comp, initOpt, initScale, seedOpt];
    If[x0Vec === $Failed, Return[$Failed]];
    optRes = dispatchOpt[
        Function[v, evaluateCompiledEnsemble[comp, v]],
        x0Vec, comp, optsList,
        Function[v, evaluateCompiledEnsembleObjective[comp, v]]
    ];
    If[optRes === $Failed, Return[$Failed]];
    assembleEnsembleResult[comp, optRes]
];

runEnsembleWithSchedules[ens_OptimalControlEnsemble, preComp_, tsList_List, smList_List, mergedOpts_List] := Module[
    {tAllStart = AbsoluteTime[], monKey = $OCSolverRegistry["MonitorStageKey"],
     maxIterOpt, tLimit, defaultNt, dt0, totalDur, prevRes = None, compCur = preComp,
     contHist = {}, totIters = 0, totEvals = 0, totObjEvals = 0, totGradEvals = 0,
     hasTS = Length[tsList] >= 1, hasSM = Length[smList] >= 1, stageIdx = 0,
     betaInit, s, ntK, remTime, ensK, compK, stageOpts, stageRes, betaVal, smRem,
     totalStages, appendStage},
    maxIterOpt = MaxIterations /. mergedOpts;
    tLimit = TimeConstraint /. mergedOpts;
    totalStages = Length[tsList] + If[hasTS && hasSM, Max[0, Length[smList] - 1], Length[smList]];

    appendStage[kind_String, metaRules_List, res_OptimalControlResult] := Module[{},
        totIters += res["Iterations"];
        totEvals += res["Evaluations"];
        totObjEvals += res["ObjectiveEvaluations"];
        totGradEvals += res["GradientEvaluations"];
        AppendTo[contHist, Join[<|"Stage" -> stageIdx, "Kind" -> kind|>, Association[metaRules], <|
            "Fidelity" -> res["FinalFidelity"],
            "Objective" -> res["FinalObjective"],
            "MeanFidelity" -> res["MeanMemberFidelity"],
            "MinFidelity" -> res["MinMemberFidelity"],
            "Iterations" -> res["Iterations"],
            "Evaluations" -> res["Evaluations"],
            "ElapsedTime" -> res["ElapsedTime"],
            "Status" -> res["Status"]
        |>]]
    ];

    If[!hasTS && !hasSM,
        Return[optimizeEnsembleCompiledStart[compCur, mergedOpts]]
    ];

    If[hasTS,
        defaultNt = TimeSteps /. mergedOpts;
        dt0 = extractEnsembleTiming[ens, defaultNt];
        If[dt0 === $Failed, Return[$Failed]];
        totalDur = Total[dt0];
        betaInit = If[hasSM, First[smList], None];
        Do[
            stageIdx++;
            ntK = tsList[[s]];
            remTime = If[NumericQ[tLimit], Max[0.01, tLimit - (AbsoluteTime[] - tAllStart)], Infinity];
            ensK = retargetEnsembleTiming[ens, totalDur, ntK];
            compK = compileEnsembleIfNeeded[ensK,
                Sequence @@ FilterRules[mergedOpts, Except[SoftMinSchedule | TimeStepSchedule | MultiStart | InitialControls | RandomSeed | StepMonitor | MonitorInterval]]];
            If[compK === $Failed, Return[$Failed, Module]];
            If[betaInit =!= None,
                compK = Join[compK, <|"Aggregation" -> If[betaInit == 0., "Mean", "SoftMin"], "SoftMinBeta" -> N[betaInit]|>]
            ];
            stageOpts = Join[
                If[s > 1, {InitialControls -> prevRes}, {}],
                If[stageIdx < totalStages, {MaxIterations -> Max[15, Ceiling[0.5 * maxIterOpt]]}, {}],
                {TimeConstraint -> remTime,
                 monKey -> <|"Kind" -> "TimeStep", "StageIndex" -> stageIdx, "TotalStages" -> totalStages, "NumSlices" -> ntK|>},
                mergedOpts
            ];
            stageRes = optimizeEnsembleCompiledStart[compK, stageOpts];
            If[stageRes === $Failed, Return[$Failed, Module]];
            prevRes = stageRes;
            compCur = compK;
            appendStage["TimeStep", {"NumSlices" -> ntK}, stageRes],
            {s, Length[tsList]}
        ]
    ];

    smRem = Which[!hasSM, {}, hasTS, Rest[smList], True, smList];
    Do[
        stageIdx++;
        betaVal = smRem[[s]];
        remTime = If[NumericQ[tLimit], Max[0.01, tLimit - (AbsoluteTime[] - tAllStart)], Infinity];
        compK = Join[compCur, <|
            "Aggregation" -> If[betaVal == 0., "Mean", "SoftMin"],
            "SoftMinBeta" -> N[betaVal]
        |>];
        stageOpts = Join[
            If[prevRes =!= None, {InitialControls -> prevRes}, {}],
            If[stageIdx < totalStages, {MaxIterations -> Max[15, Ceiling[0.5 * maxIterOpt]]}, {}],
            {TimeConstraint -> remTime,
             monKey -> <|"Kind" -> "SoftMin", "StageIndex" -> stageIdx, "TotalStages" -> totalStages,
                        "Beta" -> If[betaVal == 0., "Mean", betaVal]|>},
            mergedOpts
        ];
        stageRes = optimizeEnsembleCompiledStart[compK, stageOpts];
        If[stageRes === $Failed, Return[$Failed, Module]];
        prevRes = stageRes;
        compCur = compK;
        appendStage["SoftMin", {"Beta" -> If[betaVal == 0., "Mean", betaVal]}, stageRes],
        {s, Length[smRem]}
    ];

    OptimalControlResult[Join[prevRes[[1]], <|
        "Iterations" -> totIters,
        "Evaluations" -> totEvals,
        "ObjectiveEvaluations" -> totObjEvals,
        "GradientEvaluations" -> totGradEvals,
        "ElapsedTime" -> AbsoluteTime[] - tAllStart,
        "ContinuationHistory" -> contHist
    |>]]
];

(* Optimization driver for OptimalControlEnsemble & lists/CooperativeObjective of OptimalControlProblem *)
optimizeEnsembleInternal[ens_OptimalControlEnsemble, opts___Rule] := Module[
    {ensOwnOpts, mergedOpts, aggOpt, tsOpt, smOpt, msOpt, monOpt, monIntOpt, durPenOpt,
     explicitTFQ, tsList, smList, msSpec, preComp = None, singleRunner},
    ensOwnOpts = Which[
        Length[ens] >= 2 && ListQ[ens[[1]]] && !MatrixQ[ens[[1]]], List @@ ens[[2 ;;]],
        Length[ens] >= 5, List @@ ens[[5 ;;]],
        True, {}
    ];
    mergedOpts = Join[{opts}, ensOwnOpts, Options[OptimalControlEnsemble], Options[GRAPEOptimize]];
    durPenOpt = DurationPenalty /. mergedOpts;
    explicitTFQ = MemberQ[Keys[Join[{opts}, ensOwnOpts]], TargetFidelity];
    If[!explicitTFQ && (NumericQ[durPenOpt] && TrueQ[durPenOpt > 0.] || MatchQ[durPenOpt, {_?NumericQ, _?NumericQ}] && TrueQ[durPenOpt[[1]] > 0.]),
        mergedOpts = Prepend[mergedOpts, TargetFidelity -> None]
    ];
    {aggOpt, tsOpt, smOpt, msOpt, monOpt, monIntOpt} =
        {EnsembleAggregation, TimeStepSchedule, SoftMinSchedule, MultiStart, StepMonitor, MonitorInterval} /. mergedOpts;

    If[($OCSolverRegistry["ParseMonitorConfig"])[monOpt, monIntOpt] === $Failed, Return[$Failed]];
    tsList = ($OCSolverRegistry["ParseTimeStepSchedule"])[tsOpt];
    msSpec = ($OCSolverRegistry["ParseMultiStartSpec"])[msOpt];
    If[tsList === $Failed || msSpec === $Failed, Return[$Failed]];

    smList = parseSoftMinSchedule[smOpt, aggOpt];
    If[smList === $Failed, Return[$Failed]];

    (* Validate EnsembleAggregation eagerly even when SoftMinSchedule overrides early stages *)
    If[!MatchQ[aggOpt, "Mean" | Automatic] && !(MatchQ[aggOpt, {"SoftMin", _?NumericQ}] && TrueQ[aggOpt[[2]] > 0]),
        Message[OptimalControlEnsemble::badagg, aggOpt]; Return[$Failed]
    ];

    (* Precompile ensemble once whenever TimeStepSchedule is inactive *)
    If[tsList === {},
        preComp = compileEnsembleIfNeeded[ens, opts];
        If[preComp === $Failed, Return[$Failed]]
    ];

    If[msSpec["Mode"] === "Single",
        Return[runEnsembleWithSchedules[ens, preComp, tsList, smList, mergedOpts]]
    ];

    singleRunner = Function[{sInit, sSeed, sIdx, nTot},
        runEnsembleWithSchedules[ens, preComp, tsList, smList,
            Join[
                {InitialControls -> sInit, RandomSeed -> sSeed,
                 ($OCSolverRegistry["MonitorStageKey"]) -> <|"Kind" -> "MultiStart", "StartIndex" -> sIdx, "TotalStarts" -> nTot|>},
                FilterRules[mergedOpts, Except[MultiStart | InitialControls | RandomSeed]]
            ]
        ]
    ];
    ($OCSolverRegistry["RunMultiStart"])[singleRunner, msSpec, mergedOpts]
];

GRAPEOptimize[ens_OptimalControlEnsemble, opts___Rule] := optimizeEnsembleInternal[ens, opts];

GRAPEOptimize[probs:{(_OptimalControlProblem | {_?NumericQ, _OptimalControlProblem}) ..}, opts___Rule] :=
    optimizeEnsembleInternal[OptimalControlEnsemble[probs], opts];

GRAPEOptimize[CooperativeObjective[probs:{(_OptimalControlProblem | {_?NumericQ, _OptimalControlProblem}) ..}, cOpts___Rule], opts___Rule] :=
    optimizeEnsembleInternal[OptimalControlEnsemble[probs, cOpts], opts];

EvaluateEnsembleFidelity[ens_OptimalControlEnsemble, input_] := Module[{comp, inArg, ev},
    comp = compileEnsembleIfNeeded[ens];
    If[comp === $Failed, Return[$Failed]];
    inArg = If[MatchQ[input, _OptimalControlResult],
        If[TrueQ[Lookup[comp, "HasDurationOptimization", False]] || TrueQ[Lookup[comp, "HasPhaseFreedom", False]] ||
           (Lookup[input[[1]], "OptimizeDurations", False] =!= False && Lookup[input[[1]], "OptimizeDurations", False] =!= None) ||
           AssociationQ[Lookup[input[[1]], "OptimalPhases", None]],
            <|"Parameters" -> input["Parameters"], "Controls" -> input["Controls"], "TimeSteps" -> input["TimeSteps"],
              "OptimalPhases" -> Lookup[input[[1]], "OptimalPhases", None]|>,
            input["Controls"]
        ],
        input
    ];
    ev = evaluateCompiledEnsemble[comp, inArg];
    If[ev === $Failed, Return[$Failed]];
    <|
        "MemberFidelities" -> ev["MemberFidelities"],
        "Weights" -> comp["Weights"],
        "MeanFidelity" -> ev["MeanFidelity"],
        "MinFidelity" -> ev["MinFidelity"],
        "MaxFidelity" -> ev["MaxFidelity"],
        "AggregatedFidelity" -> ev["Fidelity"],
        "ParameterPoints" -> comp["ParameterPoints"],
        "GridShape" -> comp["GridShape"]
    |>
];

EvaluateEnsembleFidelity[res:OptimalControlResult[assoc_Association], ___Rule] := Module[{comp},
    If[assoc["Kind"] =!= "Ensemble",
        <|"MemberFidelities" -> {assoc["FinalFidelity"]}, "Weights" -> {1.},
          "MeanFidelity" -> assoc["FinalFidelity"], "MinFidelity" -> assoc["FinalFidelity"],
          "MaxFidelity" -> assoc["FinalFidelity"], "ParameterPoints" -> None, "GridShape" -> {1}|>,
        <|
            "MemberFidelities" -> assoc["MemberFidelities"],
            "Weights" -> assoc["MemberWeights"],
            "MeanFidelity" -> assoc["MeanMemberFidelity"],
            "MinFidelity" -> assoc["MinMemberFidelity"],
            "MaxFidelity" -> assoc["MaxMemberFidelity"],
            "AggregatedFidelity" -> assoc["FinalFidelity"],
            "ParameterPoints" -> assoc["ParameterPoints"],
            "GridShape" -> assoc["GridShape"]
        |>
    ]
];

(* Independent ensemble verification against SpinDynamica native propagators *)
verifyEnsembleSimulation[res:OptimalControlResult[assoc_Association], opts___Rule] := Module[
    {tol, chronOpt, compEns, members, weights, sdFids, ocFids, memberRes, dummyRes,
     aggKind, beta, fMin, expW, sdAgg, ocAgg},
    {tol, chronOpt} = {Tolerance, Chronology} /. {opts} /. Options[VerifyOptimalControlSimulation];
    compEns = assoc["CompiledEnsemble"];
    members = compEns["Members"];
    weights = compEns["Weights"];
    ocFids = assoc["MemberFidelities"];
    ocAgg = assoc["FinalFidelity"];

    sdFids = Table[
        dummyRes = OptimalControlResult[Join[assoc, <|
            "Kind" -> "Single",
            "FinalFidelity" -> ocFids[[m]],
            "FinalPropagator" -> None,
            "FinalState" -> None,
            "CompiledProblem" -> members[[m]],
            "Space" -> members[[m, "Space"]],
            "Basis" -> members[[m, "Basis"]],
            "OperatorBasis" -> members[[m, "OperatorBasis"]],
            "SymbolicDrift" -> members[[m, "SymbolicDrift"]],
            "SymbolicControls" -> members[[m, "SymbolicControls"]],
            "DriftUnitFactor" -> members[[m, "DriftUnitFactor"]],
            "ControlUnitFactor" -> members[[m, "ControlUnitFactor"]],
            "TimeSteps" -> Lookup[assoc, "TimeSteps", members[[m, "TimeSteps"]]],
            "TotalDuration" -> Lookup[assoc, "TotalDuration", members[[m, "TotalDuration"]]],
            "InitialTimePoint" -> Lookup[members[[m]], "InitialTimePoint", 0.]
        |>]];
        memberRes = VerifyOptimalControlSimulation[dummyRes, Tolerance -> tol, Chronology -> chronOpt];
        memberRes["SpinDynamicaFidelity"],
        {m, Length[members]}
    ];

    aggKind = compEns["Aggregation"];
    beta = compEns["SoftMinBeta"];
    sdAgg = If[aggKind === "SoftMin",
        fMin = Min[sdFids];
        expW = weights * Exp[-beta * (sdFids - fMin)];
        fMin - (1. / beta) * Log[Total[expW]],
        N[weights . sdFids]
    ];

    <|
        "FidelityMatchQ" -> TrueQ[Max[Abs[sdFids - ocFids]] <= tol && Abs[sdAgg - ocAgg] <= tol],
        "OptimalControlFidelity" -> ocAgg,
        "SpinDynamicaFidelity" -> sdAgg,
        "FidelityError" -> Abs[sdAgg - ocAgg],
        "MaxMemberError" -> Max[Abs[sdFids - ocFids]],
        "MemberFidelityErrors" -> Abs[sdFids - ocFids]
    |>
];

(* Register Ensemble hooks into GRAPE.m $OCSolverRegistry *)
$OCSolverRegistry["IsEnsemble"] = Function[x,
    MatchQ[x, _OptimalControlEnsemble] || (AssociationQ[x] && Lookup[x, "Kind", None] === "Ensemble")
];
$OCSolverRegistry["CompileIfEnsemble"] = compileEnsembleIfNeeded;
$OCSolverRegistry["EvaluateEnsemble"] = Function[{comp, input}, evaluateCompiledEnsemble[comp, input]];
$OCSolverRegistry["EvaluateEnsembleObjective"] = Function[{comp, input}, evaluateCompiledEnsembleObjective[comp, input]];
$OCSolverRegistry["CheckEnsembleSecondOrderCapabilities"] = checkEnsembleSecondOrderCapabilities;
$OCSolverRegistry["EvaluateEnsembleWithCache"] = Function[{comp, input}, evaluateEnsembleWithCache[comp, input]];
$OCSolverRegistry["EvaluateEnsembleHVP"] = Function[{comp, state, dir}, evaluateEnsembleHVP[comp, state, dir]];
$OCSolverRegistry["EvaluateEnsembleHessian"] = evaluateEnsembleHessian;
$OCSolverRegistry["VerifyEnsemble"] = Function[{res, optsList}, verifyEnsembleSimulation[res, Sequence @@ optsList]];
$OCSolverRegistry["EvaluateEnsembleFidelity"] = Function[res, EvaluateEnsembleFidelity[res]];
$OCSolverRegistry["RetargetEnsembleTiming"] = retargetEnsembleTiming;
$OCSolverRegistry["ExtractEnsembleTiming"] = extractEnsembleTiming;
$OCSolverRegistry["AssembleEnsembleResult"] = assembleEnsembleResult;
$OCSolverRegistry["ParseSoftMinSchedule"] = parseSoftMinSchedule;
$OCSolverRegistry["OptimizeCooperativeOrEnsemble"] = Function[{drift, ctrls, target, timing, optsList},
    Module[{items = If[MatchQ[target, _CooperativeObjective], First[target], target]},
        If[MatchQ[items, {(_OptimalControlProblem | {_?NumericQ, _OptimalControlProblem}) ..}],
            optimizeEnsembleInternal[OptimalControlEnsemble[items], Sequence @@ optsList],
            $Failed
        ]
    ]
];

End[];
EndPackage[];
