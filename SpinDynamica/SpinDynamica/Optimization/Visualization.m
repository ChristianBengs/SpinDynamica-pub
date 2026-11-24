(* ::Package:: *)

(* SpinDynamica Optimal Control — Waveform, Convergence, Trajectory & Ensemble Visualization *)

BeginPackage["OptimalControl`Visualization`", {
    "Global`",
    "SpinDynamica`",
    "QM`Hilbert`",
    "QM`SpinOperators`",
    "QM`Liouville`",
    "OptimalControl`Kernels`",
    "OptimalControl`Objectives`",
    "OptimalControl`GRAPE`",
    "OptimalControl`Ensemble`",
    "OptimalControl`TimeOptimal`"
}];

PlotOptimalControlWaveform::usage = "PlotOptimalControlWaveform[result, opts] plots the optimized piecewise-constant control waveform channels as a function of time.";
PlotOptimalControlConvergence::usage = "PlotOptimalControlConvergence[result, opts] plots the fidelity (or infidelity) and objective history across optimization iterations.";
PlotOptimalControlTrajectory::usage = "PlotOptimalControlTrajectory[result, initialState, observables, opts] propagates an initial state or operator under the optimized pulse and plots observable expectation trajectories.";
PlotEnsembleFidelity::usage = "PlotEnsembleFidelity[ensembleResultOrEval, opts] plots the member fidelity distribution across a 1D slice, 2D parameter grid, or ensemble index.";
PlotTimeOptimalPulseCurve::usage = "PlotTimeOptimalPulseCurve[curveOrOptResult, opts] plots fidelity, log-infidelity, ensemble statistics, shortest feasible waveform, or waveform evolution versus pulse duration (Method -> \"Fidelity\" | \"Infidelity\" | \"Ensemble\" | \"ShortestFeasible\" | \"WaveformEvolution\").";
PlotSwitchingFunctions::usage = "PlotSwitchingFunctions[stationarityAssocOrResult, opts] plots the control waveforms, continuous Pontryagin switching functions s_k(t), and pointwise stationarity residual r(t) over physical time.";

Begin["`Private`"];

Options[PlotOptimalControlWaveform] = {
    Method -> "Cartesian",
    ControlUnits -> Automatic,
    PlotLegends -> Automatic,
    FrameLabel -> Automatic,
    PlotLabel -> Automatic,
    ImageSize -> Medium
};

Options[PlotOptimalControlConvergence] = {
    Method -> "Infidelity",
    PlotLegends -> Automatic,
    FrameLabel -> Automatic,
    PlotLabel -> Automatic,
    ImageSize -> Medium
};

Options[PlotOptimalControlTrajectory] = {
    TargetNormalization -> "Frobenius",
    PlotLegends -> Automatic,
    FrameLabel -> Automatic,
    PlotLabel -> Automatic,
    ImageSize -> Medium
};

Options[PlotEnsembleFidelity] = {
    PlotLegends -> Automatic,
    FrameLabel -> Automatic,
    PlotLabel -> Automatic,
    ImageSize -> Medium
};

Options[PlotTimeOptimalPulseCurve] = {
    Method -> "Fidelity",
    PlotLegends -> Automatic,
    FrameLabel -> Automatic,
    PlotLabel -> Automatic,
    ImageSize -> Medium
};

Options[PlotSwitchingFunctions] = {
    StationaritySubdivisions -> 8,
    Tolerance -> 10^-3,
    PlotLegends -> Automatic,
    FrameLabel -> Automatic,
    PlotLabel -> Automatic,
    ImageSize -> Medium
};

buildStepSeries[tEdges_List, vals_List] := Module[{nt = Length[vals], pts},
    pts = Table[Null, {2 nt}];
    Do[
        pts[[2 j - 1]] = {tEdges[[j]], vals[[j]]};
        pts[[2 j]] = {tEdges[[j + 1]], vals[[j]]},
        {j, nt}
    ];
    pts
];

PlotOptimalControlWaveform[OptimalControlResult[assoc_Association], opts___Rule] := Module[
    {modeOpt, unitOpt, legOpt, frameOpt, labelOpt, sizeOpt,
     uMat, dtVec, tEdges, mCtrl, unitFactor, baseFactor, plotUnits, scaleRatio,
     uPlot, seriesList, legends, defFrame, amp, phase, pOpts},
    {modeOpt, unitOpt, legOpt, frameOpt, labelOpt, sizeOpt} =
        {Method, ControlUnits, PlotLegends, FrameLabel, PlotLabel, ImageSize} /. {opts} /. Options[PlotOptimalControlWaveform];

    uMat = assoc["Controls"];
    dtVec = assoc["TimeSteps"];
    tEdges = N[Lookup[assoc, "InitialTimePoint", 0.]] + Prepend[Accumulate[N[dtVec]], 0.];
    mCtrl = Dimensions[uMat][[2]];

    baseFactor = assoc["ControlUnitFactor"];
    If[unitOpt === Automatic,
        plotUnits = assoc["ControlUnits"];
        scaleRatio = 1.,
        plotUnits = unitOpt;
        unitFactor = ($OCSolverRegistry["ResolveUnitFactor"])[unitOpt];
        If[!NumericQ[unitFactor] || unitFactor == 0.,
            scaleRatio = 1.; plotUnits = assoc["ControlUnits"],
            scaleRatio = N[baseFactor / unitFactor]
        ]
    ];
    uPlot = scaleRatio * uMat;
    pOpts = FilterRules[{opts}, Complement[Options[ListLinePlot], {FrameLabel, PlotLabel, PlotLegends, ImageSize}]];

    If[modeOpt === "AmplitudePhase" && mCtrl === 2,
        amp = Sqrt[uPlot[[All, 1]]^2 + uPlot[[All, 2]]^2];
        phase = MapThread[If[#1 == 0. && #2 == 0., 0., ArcTan[#1, #2] / Degree] &, {uPlot[[All, 1]], uPlot[[All, 2]]}];
        Return[GraphicsColumn[{
            ListLinePlot[buildStepSeries[tEdges, amp],
                Frame -> True,
                FrameLabel -> {"Time (s)", "Amplitude (" <> ToString[plotUnits] <> ")"},
                PlotLabel -> If[labelOpt === Automatic, "Optimal Control Amplitude", labelOpt],
                ImageSize -> sizeOpt, PlotRange -> All, Sequence @@ pOpts],
            ListLinePlot[buildStepSeries[tEdges, phase],
                Frame -> True,
                FrameLabel -> {"Time (s)", "Phase (deg)"},
                PlotLabel -> "Optimal Control Phase",
                ImageSize -> sizeOpt, PlotRange -> {-185, 185}, Sequence @@ pOpts]
        }]],
        seriesList = Table[buildStepSeries[tEdges, uPlot[[All, k]]], {k, mCtrl}];
        legends = Which[
            legOpt =!= Automatic, legOpt,
            mCtrl === 1, {"Control 1"},
            mCtrl === 2, {"Control 1 (x)", "Control 2 (y)"},
            True, Table["Control " <> ToString[k], {k, mCtrl}]
        ];
        defFrame = {"Time (s)", "Amplitude (" <> ToString[plotUnits] <> ")"};
        ListLinePlot[seriesList,
            Frame -> True,
            FrameLabel -> If[frameOpt === Automatic, defFrame, frameOpt],
            PlotLabel -> If[labelOpt === Automatic, "Optimal Control Waveform", labelOpt],
            PlotLegends -> legends,
            ImageSize -> sizeOpt,
            PlotRange -> All,
            Sequence @@ pOpts
        ]
    ]
];

PlotOptimalControlConvergence[OptimalControlResult[assoc_Association], opts___Rule] := Module[
    {modeOpt, legOpt, frameOpt, labelOpt, sizeOpt, hist, fidHist, objHist, infidHist, pOpts},
    {modeOpt, legOpt, frameOpt, labelOpt, sizeOpt} =
        {Method, PlotLegends, FrameLabel, PlotLabel, ImageSize} /. {opts} /. Options[PlotOptimalControlConvergence];
    hist = assoc["ConvergenceHistory"];
    fidHist = N[hist["Fidelity"]];
    objHist = N[hist["Objective"]];
    If[fidHist === {} || !VectorQ[fidHist, NumericQ],
        Return[Graphics[Text["No convergence history recorded."], ImageSize -> sizeOpt]]
    ];
    pOpts = FilterRules[{opts}, Complement[Options[ListLinePlot], {FrameLabel, PlotLabel, PlotLegends, ImageSize}]];
    Switch[modeOpt,
        "Fidelity",
            ListLinePlot[{Transpose[{Range[0, Length[fidHist] - 1], fidHist}],
                          Transpose[{Range[0, Length[objHist] - 1], objHist}]},
                Frame -> True,
                FrameLabel -> If[frameOpt === Automatic, {"Iteration", "Value"}, frameOpt],
                PlotLabel -> If[labelOpt === Automatic, "Optimization Convergence", labelOpt],
                PlotLegends -> If[legOpt === Automatic, {"Fidelity", "Objective"}, legOpt],
                PlotRange -> All,
                ImageSize -> sizeOpt,
                Sequence @@ pOpts
            ],
        _,
            infidHist = Max[1. - #, 10^-16] & /@ fidHist;
            ListLogPlot[Transpose[{Range[0, Length[infidHist] - 1], infidHist}],
                Joined -> True,
                Frame -> True,
                FrameLabel -> If[frameOpt === Automatic, {"Iteration", "Infidelity (1 - F)"}, frameOpt],
                PlotLabel -> If[labelOpt === Automatic, "Infidelity Convergence", labelOpt],
                PlotRange -> All,
                ImageSize -> sizeOpt,
                Sequence @@ pOpts
            ]
    ]
];

toStateOpMat[expr_, basis_, dim_Integer] := Module[{v, m},
    Which[
        KetQ[expr],
            v = ToPackedComplexVector[Flatten[Normal[VectorRepresentation[expr, basis]]]];
            v = v / Max[Norm[v], 10^-16];
            ToPackedComplexMatrix[Outer[Times, v, Conjugate[v]]],
        OperatorQ[expr],
            ToPackedComplexMatrix[MatrixRepresentation[expr, basis]],
        MatrixQ[expr, NumericQ] && Dimensions[expr] === {dim, dim},
            ToPackedComplexMatrix[expr],
        VectorQ[expr, NumericQ] && Length[expr] === dim,
            v = ToPackedComplexVector[expr / Max[Norm[expr], 10^-16]];
            ToPackedComplexMatrix[Outer[Times, v, Conjugate[v]]],
        True,
            $Failed
    ]
];

toLiouvVec[expr_, opbasis_, basis_, dl_Integer] := Module[{v},
    Which[
        OperatorQ[expr],
            ToPackedComplexVector[Flatten[Normal[OperatorVectorRepresentation[expr, opbasis, basis]]]],
        KetQ[expr],
            ToPackedComplexVector[Flatten[Normal[OperatorVectorRepresentation[expr . Adjoint[expr], opbasis, basis]]]],
        VectorQ[expr, NumericQ] && Length[expr] === dl,
            ToPackedComplexVector[expr],
        MatrixQ[expr, NumericQ] && Length[expr]^2 === dl,
            v = Operator[expr, basis, CheckOperatorType -> False];
            ToPackedComplexVector[Flatten[Normal[OperatorVectorRepresentation[v, opbasis, basis]]]],
        True,
            $Failed
    ]
];

PlotOptimalControlTrajectory[OptimalControlResult[assoc_Association], initState_, obsSpec_, opts___Rule] := Module[
    {normOpt, legOpt, frameOpt, labelOpt, sizeOpt, comp, compEff, propRes, dtVec, dtComp, t0Val, schedNew,
     space, basis, opbasis, dim, obsList, nObs, rho0, obsMats, uCum, states, rMat,
     v0, obsVecs, pSlice, vStates, scaleVals, traces, legends, pOpts, j, r},
    {normOpt, legOpt, frameOpt, labelOpt, sizeOpt} =
        {TargetNormalization, PlotLegends, FrameLabel, PlotLabel, ImageSize} /. {opts} /. Options[PlotOptimalControlTrajectory];

    comp = If[assoc["Kind"] === "Ensemble",
        First[assoc["CompiledEnsemble", "Members"]],
        assoc["CompiledProblem"]
    ];
    dtVec = Developer`ToPackedArray[N[Lookup[assoc, "TimeSteps", comp["TimeSteps"]]]];
    dtComp = Developer`ToPackedArray[N[comp["TimeSteps"]]];
    t0Val = N[Lookup[assoc, "InitialTimePoint", Lookup[comp, "InitialTimePoint", 0.]]];
    compEff = If[Length[dtVec] =!= Length[dtComp] || Max[Abs[dtVec - dtComp]] > 10^-13 * Max[dtComp, 10^-12],
        If[Lookup[comp, "DriftKind", "Static"] === "Continuous" && Lookup[comp, "DriftEvaluator", None] =!= None,
            schedNew = ($OCSolverRegistry["BuildFineDriftSchedule"])[
                comp["DriftEvaluator"], dtVec, t0Val, comp["DriftIntegrationSpec"], comp["Space"]
            ];
            If[AssociationQ[schedNew],
                Join[comp, <|"TimeSteps" -> dtVec, "FineTimeSteps" -> schedNew["FineTimeSteps"], "Drift" -> schedNew["Drift"]|>],
                Join[comp, <|"TimeSteps" -> dtVec, "FineTimeSteps" -> dtVec|>]
            ],
            Join[comp, <|"TimeSteps" -> dtVec, "FineTimeSteps" -> dtVec|>]
        ],
        comp
    ];
    propRes = PropagateControls[compEff, assoc["Controls"]];
    space = compEff["Space"];
    basis = compEff["Basis"];
    opbasis = compEff["OperatorBasis"];
    dim = compEff["Dimension"];
    tGrid = t0Val + Prepend[Accumulate[dtVec], 0.];

    obsList = If[ListQ[obsSpec] && !MatrixQ[obsSpec] && !VectorQ[obsSpec, NumericQ], obsSpec, {obsSpec}];
    nObs = Length[obsList];

    If[space === "Hilbert",
        rho0 = toStateOpMat[initState, basis, dim];
        obsMats = Map[toStateOpMat[#, basis, dim] &, obsList];
        If[rho0 === $Failed || !FreeQ[obsMats, $Failed], Return[$Failed]];
        uCum = propRes["CumulativePropagators"];
        states = Table[uCum[[j]] . rho0 . ConjugateTranspose[uCum[[j]]], {j, Length[tGrid]}];
        scaleVals = Table[
            Which[
                normOpt === "TargetNorm", Max[Total[Abs[Flatten[obsMats[[r]]]]^2], 10^-16],
                normOpt === "None" || normOpt === False, 1.,
                True, Max[Norm[Flatten[rho0]] * Norm[Flatten[obsMats[[r]]]], 10^-16]
            ],
            {r, nObs}
        ];
        traces = Table[
            Table[{tGrid[[j]], Re[Tr[ConjugateTranspose[obsMats[[r]]] . states[[j]]]] / scaleVals[[r]]}, {j, Length[tGrid]}],
            {r, nObs}
        ],
        (* Liouville propagation *)
        v0 = toLiouvVec[initState, opbasis, basis, dim];
        obsVecs = Map[toLiouvVec[#, opbasis, basis, dim] &, obsList];
        If[v0 === $Failed || !FreeQ[obsVecs, $Failed], Return[$Failed]];
        pSlice = propRes["Propagators"];
        vStates = Table[Null, {Length[tGrid]}];
        vStates[[1]] = v0;
        Do[vStates[[j + 1]] = pSlice[[j]] . vStates[[j]], {j, Length[pSlice]}];
        scaleVals = Table[
            Which[
                normOpt === "TargetNorm", Max[Total[Abs[obsVecs[[r]]]^2], 10^-16],
                normOpt === "None" || normOpt === False, 1.,
                True, Max[Norm[v0] * Norm[obsVecs[[r]]], 10^-16]
            ],
            {r, nObs}
        ];
        traces = Table[
            Table[{tGrid[[j]], Re[Conjugate[obsVecs[[r]]] . vStates[[j]]] / scaleVals[[r]]}, {j, Length[tGrid]}],
            {r, nObs}
        ]
    ];

    legends = If[legOpt === Automatic,
        Table[ToString[obsList[[r]], InputForm], {r, nObs}],
        legOpt
    ];
    pOpts = FilterRules[{opts}, Complement[Options[ListLinePlot], {FrameLabel, PlotLabel, PlotLegends, ImageSize}]];
    ListLinePlot[traces,
        Frame -> True,
        FrameLabel -> If[frameOpt === Automatic, {"Time (s)", "Expectation"}, frameOpt],
        PlotLabel -> If[labelOpt === Automatic, "Expectation Trajectory", labelOpt],
        PlotLegends -> legends,
        PlotRange -> All,
        ImageSize -> sizeOpt,
        Sequence @@ pOpts
    ]
];

PlotEnsembleFidelity[input_, opts___Rule] := Module[
    {legOpt, frameOpt, labelOpt, sizeOpt, eval, fids, pts, gShape, mat, xVals, yVals, xTicks, yTicks},
    {legOpt, frameOpt, labelOpt, sizeOpt} =
        {PlotLegends, FrameLabel, PlotLabel, ImageSize} /. {opts} /. Options[PlotEnsembleFidelity];
    eval = Which[
        MatchQ[input, _OptimalControlResult] && KeyExistsQ[$OCSolverRegistry, "EvaluateEnsembleFidelity"],
            ($OCSolverRegistry["EvaluateEnsembleFidelity"])[input],
        AssociationQ[input] && KeyExistsQ[input, "MemberFidelities"],
            input,
        True,
            $Failed
    ];
    If[!AssociationQ[eval], Return[$Failed]];
    fids = N[eval["MemberFidelities"]];
    pts = Lookup[eval, "ParameterPoints", None];
    gShape = Lookup[eval, "GridShape", None];

    Which[
        ListQ[gShape] && Length[gShape] === 2 && Times @@ gShape === Length[fids] && gShape[[1]] >= 2 && gShape[[2]] >= 2,
            mat = ArrayReshape[fids, gShape];
            MatrixPlot[mat,
                ColorFunction -> "TemperatureMap",
                ColorFunctionScaling -> False,
                PlotRange -> {0., 1.},
                Frame -> True,
                FrameLabel -> If[frameOpt === Automatic, {"Parameter 2 Index", "Parameter 1 Index"}, frameOpt],
                PlotLabel -> If[labelOpt === Automatic, "Ensemble Fidelity Surface", labelOpt],
                PlotLegends -> If[legOpt === Automatic, Automatic, legOpt],
                ImageSize -> sizeOpt
            ],
        ListQ[pts] && Length[pts] === Length[fids] && VectorQ[pts, NumericQ],
            ListLinePlot[Transpose[{pts, fids}],
                Mesh -> All,
                Frame -> True,
                FrameLabel -> If[frameOpt === Automatic, {"Ensemble Parameter", "Fidelity"}, frameOpt],
                PlotLabel -> If[labelOpt === Automatic, "Ensemble Fidelity Profile", labelOpt],
                PlotRange -> {0., 1.02},
                ImageSize -> sizeOpt
            ],
        ListQ[pts] && MatrixQ[pts, NumericQ] && Length[pts] === Length[fids] && Dimensions[pts][[2]] === 1,
            ListLinePlot[Transpose[{pts[[All, 1]], fids}],
                Mesh -> All,
                Frame -> True,
                FrameLabel -> If[frameOpt === Automatic, {"Ensemble Parameter", "Fidelity"}, frameOpt],
                PlotLabel -> If[labelOpt === Automatic, "Ensemble Fidelity Profile", labelOpt],
                PlotRange -> {0., 1.02},
                ImageSize -> sizeOpt
            ],
        True,
            ListPlot[Transpose[{Range[Length[fids]], fids}],
                Joined -> True,
                Filling -> Axis,
                Frame -> True,
                FrameLabel -> If[frameOpt === Automatic, {"Member Index", "Fidelity"}, frameOpt],
                PlotLabel -> If[labelOpt === Automatic, "Ensemble Member Fidelities", labelOpt],
                PlotRange -> {0., 1.02},
                ImageSize -> sizeOpt
            ]
    ]
];

extractCurveRecords[input_Association] := Module[{recs, durs, fids, means, softs, mins, tarFid, bestRes},
    tarFid = Lookup[input, "TargetFidelity", None];
    bestRes = Lookup[input, "ShortestFeasibleResult", Lookup[input, "BestResult", None]];
    Which[
        KeyExistsQ[input, "Records"] && ListQ[input["Records"]],
            recs = SortBy[Select[input["Records"], NumericQ[Lookup[#, "Duration", Indeterminate]] &], #["Duration"] &];
            durs = #["Duration"] & /@ recs;
            fids = #["CriterionFidelity"] & /@ recs;
            means = Lookup[#, "MeanFidelity", #["CriterionFidelity"]] & /@ recs;
            softs = Lookup[#, "SoftMinFidelity", #["CriterionFidelity"]] & /@ recs;
            mins = Lookup[#, "MinFidelity", #["CriterionFidelity"]] & /@ recs,
        KeyExistsQ[input, "History"] && ListQ[input["History"]],
            recs = SortBy[Select[input["History"], NumericQ[Lookup[#, "Duration", Indeterminate]] &], #["Duration"] &];
            durs = #["Duration"] & /@ recs;
            fids = Lookup[#, "CriterionFidelity", Lookup[#, "Fidelity", 0.]] & /@ recs;
            means = Lookup[#, "MeanFidelity", #["Fidelity"]] & /@ recs;
            softs = Lookup[#, "SoftMinFidelity", #["Fidelity"]] & /@ recs;
            mins = Lookup[#, "MinFidelity", #["Fidelity"]] & /@ recs,
        True,
            Return[$Failed]
    ];
    <|"Records" -> recs, "Durations" -> durs, "Fidelities" -> fids, "Means" -> means,
      "SoftMins" -> softs, "Mins" -> mins, "TargetFidelity" -> tarFid, "BestResult" -> bestRes,
      "ShortestFeasibleDuration" -> Lookup[input, "ShortestFeasibleDuration", None]|>
];

PlotTimeOptimalPulseCurve[input_Association, opts___Rule] := Module[
    {modeOpt, legOpt, frameOpt, labelOpt, sizeOpt, parsed, durs, fids, means, softs, mins,
     tarFid, bestRes, tShort, recs, validRecs, pickIdx, waveSeries, waveLegends, tEdgesNorm,
     pOpts, infidPts, gridLines},
    {modeOpt, legOpt, frameOpt, labelOpt, sizeOpt} =
        {Method, PlotLegends, FrameLabel, PlotLabel, ImageSize} /. {opts} /. Options[PlotTimeOptimalPulseCurve];
    parsed = extractCurveRecords[input];
    If[parsed === $Failed, Return[$Failed]];
    {durs, fids, means, softs, mins, tarFid, bestRes, tShort, recs} =
        parsed /@ {"Durations", "Fidelities", "Means", "SoftMins", "Mins", "TargetFidelity",
                   "BestResult", "ShortestFeasibleDuration", "Records"};
    pOpts = FilterRules[{opts}, Complement[Options[ListLinePlot], {Method, FrameLabel, PlotLabel, PlotLegends, ImageSize}]];
    gridLines = {If[NumericQ[tShort], {{tShort, Directive[Dashed, Gray]}}, None],
                 If[NumericQ[tarFid], {{tarFid, Directive[Dashed, Red]}}, None]};

    Switch[modeOpt,
        "Fidelity",
            If[durs === {}, Return[Graphics[Text["No duration records available."], ImageSize -> sizeOpt]]];
            ListLinePlot[Transpose[{durs, fids}],
                Mesh -> All,
                Frame -> True,
                FrameLabel -> If[frameOpt === Automatic, {"Pulse Duration (s)", "Fidelity"}, frameOpt],
                PlotLabel -> If[labelOpt === Automatic, "Time-Optimal Pulse Curve", labelOpt],
                PlotLegends -> If[legOpt === Automatic, None, legOpt],
                GridLines -> gridLines,
                PlotRange -> All,
                ImageSize -> sizeOpt,
                Sequence @@ pOpts
            ],
        "Infidelity",
            If[durs === {}, Return[Graphics[Text["No duration records available."], ImageSize -> sizeOpt]]];
            infidPts = Transpose[{durs, Max[1. - #, 10^-16] & /@ fids}];
            ListLogPlot[infidPts,
                Joined -> True,
                Mesh -> All,
                Frame -> True,
                FrameLabel -> If[frameOpt === Automatic, {"Pulse Duration (s)", "Infidelity (1 - F)"}, frameOpt],
                PlotLabel -> If[labelOpt === Automatic, "Time-Optimal Infidelity Curve", labelOpt],
                GridLines -> {If[NumericQ[tShort], {{tShort, Directive[Dashed, Gray]}}, None],
                              If[NumericQ[tarFid], {{Max[1. - tarFid, 10^-16], Directive[Dashed, Red]}}, None]},
                PlotRange -> All,
                ImageSize -> sizeOpt,
                Sequence @@ pOpts
            ],
        "Ensemble",
            If[durs === {}, Return[Graphics[Text["No ensemble duration records available."], ImageSize -> sizeOpt]]];
            ListLinePlot[{Transpose[{durs, means}], Transpose[{durs, softs}], Transpose[{durs, mins}]},
                Mesh -> All,
                Frame -> True,
                FrameLabel -> If[frameOpt === Automatic, {"Pulse Duration (s)", "Ensemble Fidelity"}, frameOpt],
                PlotLabel -> If[labelOpt === Automatic, "Ensemble Time-Optimal Trade-off", labelOpt],
                PlotLegends -> If[legOpt === Automatic, {"Mean", "SoftMin", "Minimum"}, legOpt],
                GridLines -> gridLines,
                PlotRange -> All,
                ImageSize -> sizeOpt,
                Sequence @@ pOpts
            ],
        "ShortestFeasible",
            If[MatchQ[bestRes, _OptimalControlResult],
                PlotOptimalControlWaveform[bestRes, PlotLabel -> If[labelOpt === Automatic,
                    "Shortest Feasible Waveform (T = " <> ToString[NumberForm[bestRes["TotalDuration"], 5]] <> " s)", labelOpt],
                    ImageSize -> sizeOpt, Sequence @@ FilterRules[{opts}, Except[Method | PlotLabel | ImageSize]]],
                Graphics[Text["No feasible OptimalControlResult recorded."], ImageSize -> sizeOpt]
            ],
        "WaveformEvolution",
            validRecs = Select[recs, MatchQ[Lookup[#, "Result", None], _OptimalControlResult] &];
            If[validRecs === {},
                Return[If[MatchQ[bestRes, _OptimalControlResult],
                    PlotOptimalControlWaveform[bestRes, ImageSize -> sizeOpt],
                    Graphics[Text["No intermediate waveforms recorded."], ImageSize -> sizeOpt]
                ]]
            ];
            pickIdx = If[Length[validRecs] <= 1,
                {1},
                DeleteDuplicates[Round[Subdivide[1, Length[validRecs], Min[5, Length[validRecs] - 1]]]]
            ];
            validRecs = validRecs[[pickIdx]];
            waveSeries = Table[
                With[{resK = r["Result"], dtK = r["Result"]["TimeSteps"], uMatK = r["Result"]["Controls"]},
                    tEdgesNorm = Prepend[Accumulate[dtK], 0.] / Max[Total[dtK], 10^-15];
                    buildStepSeries[tEdgesNorm, If[Dimensions[uMatK][[2]] >= 2,
                        Sqrt[Total[uMatK[[All, 1 ;; 2]]^2, {2}]],
                        uMatK[[All, 1]]]]
                ],
                {r, validRecs}
            ];
            waveLegends = Table[
                "T = " <> ToString[NumberForm[r["Duration"], 4]] <> " s (F=" <> ToString[NumberForm[r["CriterionFidelity"], 4]] <> ")",
                {r, validRecs}
            ];
            ListLinePlot[waveSeries,
                Frame -> True,
                FrameLabel -> If[frameOpt === Automatic, {"Normalized Time (t / T)", "Control Amplitude"}, frameOpt],
                PlotLabel -> If[labelOpt === Automatic, "Waveform Evolution During Duration Reduction", labelOpt],
                PlotLegends -> If[legOpt === Automatic, waveLegends, legOpt],
                PlotRange -> All,
                ImageSize -> sizeOpt,
                Sequence @@ pOpts
            ],
        _,
            $Failed
    ]
];

PlotSwitchingFunctions[input_, opts___Rule] := Module[
    {legOpt, frameOpt, labelOpt, sizeOpt, diag, tGrid, uFine, sFine, rFine, mCtrl,
     ctrlSeries, swSeries, ctrlLegends, swLegends},
    {legOpt, frameOpt, labelOpt, sizeOpt} =
        {PlotLegends, FrameLabel, PlotLabel, ImageSize} /. {opts} /. Options[PlotSwitchingFunctions];
    diag = Which[
        AssociationQ[input] && KeyExistsQ[input, "ApplicableQ"], input,
        MatchQ[input, _OptimalControlResult],
            VerifyOptimalControlStationarity[input, Sequence @@ FilterRules[{opts}, Options[VerifyOptimalControlStationarity]]],
        True, $Failed
    ];
    If[diag === $Failed || !AssociationQ[diag], Return[$Failed]];
    If[!TrueQ[diag["ApplicableQ"]],
        Return[Graphics[
            Text[Style[StringRiffle[Prepend[diag["InapplicabilityReasons"],
                "Continuous PMP Switching Diagnostic Inapplicable:"], "\n"], 11]],
            Frame -> True, ImageSize -> sizeOpt
        ]]
    ];
    tGrid = diag["TimeGrid"];
    uFine = diag["ControlTimeSeries"];
    sFine = diag["SwitchingFunctions"];
    rFine = diag["PointwiseResiduals"];
    mCtrl = Dimensions[sFine][[2]];
    ctrlSeries = Table[Transpose[{tGrid, uFine[[All, k]]}], {k, mCtrl}];
    swSeries = Table[Transpose[{tGrid, sFine[[All, k]]}], {k, mCtrl}];
    ctrlLegends = If[legOpt === Automatic, Table["u_" <> ToString[k], {k, mCtrl}], legOpt];
    swLegends = Table["s_" <> ToString[k], {k, mCtrl}];
    GraphicsColumn[{
        ListLinePlot[ctrlSeries,
            Frame -> True,
            FrameLabel -> {"Time (s)", "Control"},
            PlotLabel -> If[labelOpt === Automatic, "Controls u(t)", labelOpt],
            PlotLegends -> ctrlLegends,
            PlotRange -> All,
            ImageSize -> sizeOpt],
        ListLinePlot[swSeries,
            Frame -> True,
            FrameLabel -> {"Time (s)", "s_k(t)"},
            PlotLabel -> "Pontryagin Switching Functions",
            PlotLegends -> swLegends,
            GridLines -> {None, {0.}},
            PlotRange -> All,
            ImageSize -> sizeOpt],
        ListLinePlot[Transpose[{tGrid, rFine}],
            Frame -> True,
            Filling -> Axis,
            FrameLabel -> If[frameOpt === Automatic, {"Time (s)", "Residual r(t)"}, frameOpt],
            PlotLabel -> "Pointwise PMP Residual r(t) (relInt = " <>
                ToString[ScientificForm[diag["RelativeIntegratedResidual"], 3]] <> ")",
            PlotRange -> All,
            ImageSize -> sizeOpt]
    }]
];

End[];
EndPackage[];
