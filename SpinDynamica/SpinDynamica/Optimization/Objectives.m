(* ::Package:: *)

(* SpinDynamica Optimal Control — Objective Constructors, Normalization & Analytical Cotangents *)

BeginPackage["OptimalControl`Objectives`", {
    "Global`",
    "SpinDynamica`",
    "General`Matrices`",
    "QM`Hilbert`",
    "QM`SpinOperators`",
    "QM`Liouville`",
    "OptimalControl`Kernels`"
}];

UnitaryTarget::usage = "UnitaryTarget[U, opts] defines a unitary-gate optimal-control target. Options include ObjectiveType -> \"UnitaryFidelity\" | \"UnitaryOverlap\" and PhaseFreedom.";
SubspaceUnitaryTarget::usage = "SubspaceUnitaryTarget[Usub, subspace, opts] defines a unitary-gate optimal-control target restricted to an embedded k-dimensional subspace with optional PhaseFreedom and LeakageWeight.";
LeakageWeight::usage = "LeakageWeight is an option for SubspaceUnitaryTarget specifying the penalty weight (>= 0) on subspace leakage L(U) = 1 - Tr[W^H W]/k.";
PhaseFreedom::usage = "PhaseFreedom is an option for UnitaryTarget and SubspaceUnitaryTarget specifying invariant phase freedom: Automatic, None, \"Global\", \"Diagonal\"/\"SubspaceDiagonal\", or parametric virtual rotations (\"LeftZ\", \"RightZ\", \"BothZ\", or \"Generators\").";
StateTarget::usage = "StateTarget[init -> target, opts] defines a state-to-state transfer target for pure kets or density operators. Options include ObjectiveType -> \"StateFidelity\" | \"StateOverlap\" | \"StateDistance\".";
ObservableTarget::usage = "ObservableTarget[init -> obs, opts] defines an observable transfer target maximizing Re Tr[obs^H rho(T)] (or |Tr[obs^H rho(T)]|^2). Options include ObjectiveType -> \"ObservableTransfer\" | \"ObservableMagnitude\", TargetNormalization -> \"Frobenius\" | \"TargetNorm\" | \"None\" | numeric.";
CoherenceTarget::usage = "CoherenceTarget[init -> coh, opts] defines a coherence-transfer optimal-control target supporting non-Hermitian spin coherence operators and phase-sensitive or phase-insensitive objectives.";
CooperativeObjective::usage = "CooperativeObjective[{spec1, spec2, ...}, opts] combines multiple transfer or unitary objectives with optional weights {weight, spec}. Options include NormalizeWeights -> True | False.";
ObjectiveType::usage = "ObjectiveType is an option for optimal-control targets selecting the mathematical fidelity or overlap formula.";
TargetNormalization::usage = "TargetNormalization is an option for ObservableTarget and CoherenceTarget specifying \"Frobenius\", \"TargetNorm\", \"None\", or a positive numeric scale factor.";
NormalizeWeights::usage = "NormalizeWeights is an option for CooperativeObjective and OptimalControlEnsemble controlling whether weights are normalized to sum to 1.";
PhaseInsensitive::usage = "PhaseInsensitive is an option for CoherenceTarget and ObservableTarget; True maximizes |Tr[O^H rho(T)]|^2 instead of Re Tr[O^H rho(T)].";
RunningCost::usage = "RunningCost is an option for OptimalControlProblem, OptimalControlEnsemble, and GRAPEOptimize specifying intermediate-time state penalties or custom trajectory costs evaluated at coarse slice boundaries.";
StatePenalty::usage = "StatePenalty[opOrState, weight, opts] specifies an intermediate-time penalty or reward on trajectory states with weight >= 0, TargetMode (Automatic, \"Expectation\", \"Magnitude\", \"StateFidelity\", \"Projector\"), TimeWeighting, Sign (-1 or +1), and TargetNormalization.";
TargetMode::usage = "TargetMode is an option for StatePenalty selecting \"Automatic\", \"Expectation\", \"Magnitude\", \"StateFidelity\", or \"Projector\".";
TimeWeighting::usage = "TimeWeighting is an option for StatePenalty specifying normalized temporal weights (\"Uniform\", \"FinalHalf\", \"LinearRamp\", explicit vector, or callable function).";

NormalizeTargetSpec::usage = "NormalizeTargetSpec[target] converts shorthand targets (unitary matrices/operators, rules, or lists of rules/targets) into canonical target objects.";
CompileObjective::usage = "CompileObjective[targetSpec, space, basis, opbasis, opts] compiles symbolic SpinDynamica targets into packed numerical evaluators and exact cotangent functions.";
CompileRunningCost::usage = "CompileRunningCost[runningCostSpec, compiledObj, space, basis, opbasis, dtVec, t0] compiles StatePenalty specifications and custom running-cost functions into fast boundary evaluators and adjoint cotangent injectors.";

ObjectiveTarget::usage = "ObjectiveTarget is a symbolic message head used for optimal-control target validation messages.";
ObjectiveTarget::invalspec = "Cannot interpret optimal-control target specification `1`.";
ObjectiveTarget::badnorm = "Invalid TargetNormalization setting `1`; expected \"Frobenius\", \"TargetNorm\", \"None\", or a positive number.";
StateTarget::notpsd = "StateFidelity requires positive semidefinite operators with positive trace; received `1`. Use ObservableTarget or ObjectiveType -> \"ObservableTransfer\" for traceless NMR spin operators.";
StateTarget::baddim = "Target dimension `1` does not match basis dimension `2`.";
StatePenalty::badspec = "Invalid RunningCost or StatePenalty specification `1`.";
ObjectiveTarget::badsub = "Invalid subspace target specification `1` or PhaseFreedom `2`.";

Begin["`Private`"];

Options[UnitaryTarget] = {ObjectiveType -> "UnitaryFidelity", PhaseFreedom -> Automatic};
Options[SubspaceUnitaryTarget] = {ObjectiveType -> Automatic, PhaseFreedom -> Automatic, LeakageWeight -> 0., PhaseInsensitive -> Automatic};
Options[StateTarget] = {ObjectiveType -> Automatic};
Options[ObservableTarget] = {ObjectiveType -> "ObservableTransfer", TargetNormalization -> "Frobenius", PhaseInsensitive -> False};
Options[CoherenceTarget] = {ObjectiveType -> Automatic, TargetNormalization -> "Frobenius", PhaseInsensitive -> False};
Options[CooperativeObjective] = {NormalizeWeights -> True};
Options[StatePenalty] = {TargetMode -> Automatic, TimeWeighting -> "Uniform", Sign -> -1, TargetNormalization -> Automatic, InitialState -> Automatic};

isKetOrVecQ[x_] := KetQ[x] || VectorQ[x, NumericQ];

toHilbertOperatorMatrix[expr_, basis_, dim_Integer] := Module[{mat},
    If[MatchQ[expr, _NullOperator | _NullSuperoperator] || (NumericQ[expr] && TrueQ[expr == 0]),
        Return[ToPackedComplexMatrix[ConstantArray[0., {dim, dim}]]]
    ];
    mat = Which[
        OperatorQ[expr], MatrixRepresentation[expr, basis],
        MatrixQ[expr, NumericQ], expr,
        True, $Failed
    ];
    If[mat === $Failed || !MatrixQ[mat, NumericQ] || Dimensions[mat] =!= {dim, dim},
        $Failed,
        ToPackedComplexMatrix[mat]
    ]
];

toKetVector[expr_, basis_, dim_Integer] := Module[{vec},
    vec = Which[
        KetQ[expr], Flatten[Normal[VectorRepresentation[expr, basis]]],
        VectorQ[expr, NumericQ] && Length[expr] === dim, expr,
        True, $Failed
    ];
    If[vec === $Failed || Length[vec] =!= dim || !VectorQ[vec, NumericQ],
        $Failed,
        ToPackedComplexVector[vec]
    ]
];

toLiouvilleVector[expr_, opbasis_, basis_, dl_Integer] := Module[{vec},
    If[MatchQ[expr, _NullOperator | _NullSuperoperator] || (NumericQ[expr] && TrueQ[expr == 0]),
        Return[ToPackedComplexVector[ConstantArray[0., dl]]]
    ];
    vec = Which[
        OperatorQ[expr], Flatten[Normal[OperatorVectorRepresentation[expr, opbasis, basis]]],
        KetQ[expr], Flatten[Normal[OperatorVectorRepresentation[expr . Adjoint[expr], opbasis, basis]]],
        VectorQ[expr, NumericQ] && Length[expr] === dl, expr,
        MatrixQ[expr, NumericQ] && (Length[expr]^2 === dl),
            Flatten[Normal[OperatorVectorRepresentation[Operator[expr, basis, CheckOperatorType -> False], opbasis, basis]]],
        True, $Failed
    ];
    If[vec === $Failed || Length[vec] =!= dl || !VectorQ[vec, NumericQ],
        $Failed,
        ToPackedComplexVector[vec]
    ]
];

psdAnalyzeMatrix[mat_, tol_: 10^-8] := Module[
    {sym, hermErr, ev, vecs, vec, tr, evClip, rhoNorm, purity, maxIdx},
    sym = HermitianSymmetrize[mat];
    hermErr = Max[Abs[mat - sym]] / Max[1., Max[Abs[mat]]];
    If[hermErr > tol, Return[<|"ValidQ" -> False, "Reason" -> "NonHermitian"|>]];
    {ev, vecs} = Eigensystem[sym];
    ev = Re[ev];
    tr = Total[ev];
    If[!TrueQ[tr > tol] || TrueQ[Min[ev] < -tol * Max[1., tr]],
        Return[<|"ValidQ" -> False, "Reason" -> "NotPSD"|>]
    ];
    evClip = Max[#, 0.] & /@ ev;
    tr = Total[evClip];
    vec = Transpose[vecs];
    rhoNorm = ToPackedComplexMatrix[vec . ((evClip / tr) * ConjugateTranspose[vec])];
    purity = Re[Total[(evClip / tr)^2]];
    maxIdx = First[Ordering[evClip, -1]];
    <|
        "ValidQ" -> True,
        "DensityMatrix" -> rhoNorm,
        "Purity" -> purity,
        "PureQ" -> (purity > 1. - 10^-9),
        "DominantKet" -> ToPackedComplexVector[vec[[All, maxIdx]]],
        "RootMatrix" -> ToPackedComplexMatrix[vec . (Sqrt[evClip / tr] * ConjugateTranspose[vec])]
    |>
];

resolveScaleFactor[normSpec_, matInit_, matTar_] := Module[{ni, ntNorm},
    Which[
        normSpec === "Frobenius" || normSpec === Automatic,
            ni = Norm[Flatten[matInit]];
            ntNorm = Norm[Flatten[matTar]];
            If[ni * ntNorm <= 0., $Failed, ni * ntNorm],
        normSpec === "TargetNorm",
            ntNorm = Total[Abs[Flatten[matTar]]^2];
            If[ntNorm <= 0., $Failed, ntNorm],
        normSpec === "None" || normSpec === False,
            1.,
        NumericQ[normSpec] && TrueQ[normSpec > 0],
            N[normSpec],
        True,
            Message[ObjectiveTarget::badnorm, normSpec]; $Failed
    ]
];

(* Normalize shorthand target syntax into canonical target head *)
NormalizeTargetSpec[t:(_UnitaryTarget | _SubspaceUnitaryTarget | _StateTarget | _ObservableTarget | _CoherenceTarget | _CooperativeObjective)] := t;

NormalizeTargetSpec[Rule[init_, tar_]] := Module[{},
    If[isKetOrVecQ[init] && isKetOrVecQ[tar],
        StateTarget[init -> tar, ObjectiveType -> "StateFidelity"],
        StateTarget[init -> tar, ObjectiveType -> Automatic]
    ]
];

NormalizeTargetSpec[items_List] := Module[{normItems},
    normItems = Map[
        Function[entry,
            Switch[entry,
                {_?NumericQ, _}, {entry[[1]], NormalizeTargetSpec[entry[[2]]]},
                _, NormalizeTargetSpec[entry]
            ]
        ],
        items
    ];
    If[!FreeQ[normItems, $Failed], $Failed, CooperativeObjective[normItems]]
];

NormalizeTargetSpec[uTar_ /; OperatorQ[uTar] || MatrixQ[uTar, NumericQ]] :=
    UnitaryTarget[uTar, ObjectiveType -> "UnitaryFidelity"];

NormalizeTargetSpec[other_] := (Message[ObjectiveTarget::invalspec, other]; $Failed);

(* Helpers for SubspaceUnitaryTarget isometry & parametric/non-parametric PhaseFreedom *)
resolveSubspaceIsometry[subSpec_, basis_, dim_Integer] := Module[{vecs, vRaw, vMat, k},
    vRaw = Which[
        ListQ[subSpec] && subSpec =!= {} && AllTrue[subSpec, IntegerQ[#] && 1 <= # <= dim &] &&
                Length[DeleteDuplicates[subSpec]] === Length[subSpec],
            (IdentityMatrix[dim, WorkingPrecision -> MachinePrecision])[[All, subSpec]],
        ListQ[subSpec] && subSpec =!= {} && AllTrue[subSpec, KetQ],
            vecs = toKetVector[#, basis, dim] & /@ subSpec;
            If[MemberQ[vecs, $Failed], Return[$Failed]];
            Transpose[vecs],
        MatrixQ[subSpec, NumericQ] && Dimensions[subSpec][[1]] === dim && 1 <= Dimensions[subSpec][[2]] <= dim,
            subSpec,
        ListQ[subSpec] && subSpec =!= {} && AllTrue[subSpec, isKetOrVecQ],
            vecs = toKetVector[#, basis, dim] & /@ subSpec;
            If[MemberQ[vecs, $Failed], Return[$Failed]];
            Transpose[vecs],
        True,
            $Failed
    ];
    If[vRaw === $Failed || !MatrixQ[vRaw, NumericQ], Return[$Failed]];
    vMat = ToPackedComplexMatrix[vRaw];
    If[Dimensions[vMat][[1]] =!= dim || Dimensions[vMat][[2]] < 1 || Dimensions[vMat][[2]] > dim,
        Return[$Failed]
    ];
    k = Dimensions[vMat][[2]];
    If[Max[Abs[ConjugateTranspose[vMat] . vMat - IdentityMatrix[k, WorkingPrecision -> MachinePrecision]]] > 10^-8,
        Return[$Failed]
    ];
    vMat
];

activeSpinLabelsFromBasis[basis_] := Module[{ss},
    ss = Quiet[Check[SpinSystem[basis], $Failed]];
    If[!ListQ[ss] || ss === {}, ss = Quiet[Check[SpinSystem[], {}]]];
    If[ListQ[ss], Map[If[ListQ[#], First[#], #] &, ss], {}]
];

diagonalizePhaseGenerators[genList_List, sideTag_String, basis_, dim_Integer, allSpins_List] := Module[
    {out, g, gMat, hermErr, ev, vec, vecH, lbl},
    out = Table[
        g = genList[[idx]];
        lbl = g;
        If[MatchQ[g, {"SpinZ", _}],
            With[{s = g[[2]]},
                If[!MemberQ[allSpins, s], Return[$Failed, Table]];
                g = opI[s, "z"];
                lbl = sideTag <> ":Iz[" <> ToString[s] <> "]"
            ],
            lbl = sideTag <> ":Gen" <> ToString[idx]
        ];
        gMat = toHilbertOperatorMatrix[g, basis, dim];
        If[gMat === $Failed, Return[$Failed, Table]];
        hermErr = Max[Abs[gMat - ConjugateTranspose[gMat]]] / Max[1., Max[Abs[gMat]]];
        If[hermErr > 10^-8, Return[$Failed, Table]];
        gMat = HermitianSymmetrize[gMat];
        {ev, vec} = Eigensystem[gMat];
        ev = Re[ev];
        vec = Transpose[vec];
        vecH = ConjugateTranspose[vec];
        <|"Mat" -> gMat, "Evals" -> Developer`ToPackedArray[ev],
          "Vec" -> ToPackedComplexMatrix[vec], "VecH" -> ToPackedComplexMatrix[vecH],
          "Label" -> lbl|>
        ,
        {idx, Length[genList]}
    ];
    If[!ListQ[out] || MemberQ[out, $Failed], $Failed, out]
];

parsePhaseFreedomSpec[pfOpt_, defaultMode_String, basis_, dim_Integer] := Module[
    {allSpins = activeSpinLabelsFromBasis[basis], effOpt, leftSpecs = {}, rightSpecs = {},
     leftGen, rightGen, modeStr},
    effOpt = If[pfOpt === Automatic, defaultMode, pfOpt];
    If[effOpt === None || effOpt === "None" || effOpt === False,
        Return[<|"ParametricQ" -> False, "BaseMode" -> "None", "Mode" -> "None"|>]
    ];
    If[effOpt === "Global",
        Return[<|"ParametricQ" -> False, "BaseMode" -> "Global", "Mode" -> "Global"|>]
    ];
    If[MemberQ[{"Diagonal", "SubspaceDiagonal"}, effOpt],
        Return[<|"ParametricQ" -> False, "BaseMode" -> "Diagonal", "Mode" -> "Diagonal"|>]
    ];

    Which[
        effOpt === "LeftZ",
            If[allSpins === {}, Return[$Failed]];
            leftSpecs = {"SpinZ", #} & /@ allSpins; modeStr = "LeftZ",
        MatchQ[effOpt, {"LeftZ", _}],
            leftSpecs = {"SpinZ", #} & /@ Flatten[{effOpt[[2]]}]; modeStr = "LeftZ",
        effOpt === "RightZ",
            If[allSpins === {}, Return[$Failed]];
            rightSpecs = {"SpinZ", #} & /@ allSpins; modeStr = "RightZ",
        MatchQ[effOpt, {"RightZ", _}],
            rightSpecs = {"SpinZ", #} & /@ Flatten[{effOpt[[2]]}]; modeStr = "RightZ",
        effOpt === "BothZ",
            If[allSpins === {}, Return[$Failed]];
            leftSpecs = {"SpinZ", #} & /@ allSpins;
            rightSpecs = leftSpecs; modeStr = "BothZ",
        MatchQ[effOpt, {"BothZ", _}],
            leftSpecs = {"SpinZ", #} & /@ Flatten[{effOpt[[2]]}];
            rightSpecs = leftSpecs; modeStr = "BothZ",
        MatchQ[effOpt, {"BothZ", _, _}],
            leftSpecs = {"SpinZ", #} & /@ Flatten[{effOpt[[2]]}];
            rightSpecs = {"SpinZ", #} & /@ Flatten[{effOpt[[3]]}]; modeStr = "BothZ",
        MatchQ[effOpt, {"Generators", {"Left", _List}}],
            leftSpecs = effOpt[[2, 2]]; modeStr = "GeneratorsLeft",
        MatchQ[effOpt, {"Generators", {"Right", _List}}],
            rightSpecs = effOpt[[2, 2]]; modeStr = "GeneratorsRight",
        MatchQ[effOpt, {"Generators", {"Both", _List, _List}}],
            leftSpecs = effOpt[[2, 2]]; rightSpecs = effOpt[[2, 3]]; modeStr = "GeneratorsBoth",
        MatchQ[effOpt, {"Generators", _List}],
            leftSpecs = effOpt[[2]]; modeStr = "Generators",
        True,
            Return[$Failed]
    ];
    If[leftSpecs === {} && rightSpecs === {}, Return[$Failed]];
    leftGen = diagonalizePhaseGenerators[leftSpecs, "Left", basis, dim, allSpins];
    rightGen = diagonalizePhaseGenerators[rightSpecs, "Right", basis, dim, allSpins];
    If[leftGen === $Failed || rightGen === $Failed, Return[$Failed]];
    <|
        "ParametricQ" -> True,
        "BaseMode" -> If[defaultMode === "None", "None", "Global"],
        "Mode" -> modeStr,
        "Left" -> leftGen,
        "Right" -> rightGen,
        "NumLeft" -> Length[leftGen],
        "NumRight" -> Length[rightGen],
        "NumPhases" -> Length[leftGen] + Length[rightGen],
        "LeftLabels" -> (#["Label"] & /@ leftGen),
        "RightLabels" -> (#["Label"] & /@ rightGen)
    |>
];

buildUnitaryLeafAssoc[uSubMat_, vMatOrNone_, dim_Integer, pfSpec_Association, lamL_?NumericQ] := Module[
    {kDim = Length[uSubMat], uSubD = ToPackedComplexMatrix[ConjugateTranspose[uSubMat]],
     hasV = (vMatOrNone =!= None), vMat, vDag, baseMode = pfSpec["BaseMode"],
     fullSubConst, zeroDimMat, baseEvalFn, baseDirFn, evalWithPhiFn, pfOut},
    vMat = If[hasV, vMatOrNone, None];
    vDag = If[hasV, ToPackedComplexMatrix[ConjugateTranspose[vMat]], None];
    fullSubConst = If[hasV, ToPackedComplexMatrix[vMat . uSubD . vDag], uSubD];
    zeroDimMat = ToPackedComplexMatrix[ConstantArray[0., {dim, dim}]];

    baseEvalFn = With[{uD = uSubD, v = vMat, vD = vDag, hasIsom = hasV, fConst = fullSubConst,
                       kN = N[kDim], bm = baseMode, lW = N[lamL]},
        Function[uFinal,
            Module[{wMat, zMat, fidVal, cFidDag, leakVal = 0., cLeakDag = 0., zSc, zDiag, zProd},
                wMat = If[hasIsom, vD . uFinal . v, uFinal];
                zMat = uD . wMat;
                Switch[bm,
                    "None",
                        fidVal = Re[Tr[zMat]] / kN;
                        cFidDag = (1. / kN) * fConst,
                    "Global",
                        zSc = Tr[zMat] / kN;
                        fidVal = Re[zSc * Conjugate[zSc]];
                        cFidDag = ((2. * Conjugate[zSc]) / kN) * fConst,
                    "Diagonal",
                        zDiag = Diagonal[zMat];
                        fidVal = Total[Re[zDiag * Conjugate[zDiag]]] / kN;
                        zProd = Conjugate[zDiag] * uD;
                        cFidDag = (2. / kN) * If[hasIsom, v . zProd . vD, zProd]
                ];
                If[hasIsom,
                    leakVal = 1. - Re[Total[Abs[Flatten[wMat]]^2]] / kN;
                    If[lW > 0.,
                        cLeakDag = (2. * lW / kN) * (v . ConjugateTranspose[wMat] . vD)
                    ]
                ];
                {fidVal - lW * leakVal,
                 ToPackedComplexMatrix[cFidDag + cLeakDag],
                 If[hasIsom, <|"Leakage" -> leakVal, "SubspaceFidelity" -> fidVal|>, <||>]}
            ]
        ]
    ];

    baseDirFn = With[{uD = uSubD, v = vMat, vD = vDag, hasIsom = hasV, fConst = fullSubConst,
                      kN = N[kDim], bm = baseMode, lW = N[lamL], zZero = zeroDimMat},
        Function[{uFinal, dU},
            Module[{dW, dzMat, dCFid = zZero, dCLeak = zZero, dzSc, dzDiag, dzProd},
                dW = If[hasIsom, vD . dU . v, dU];
                dzMat = uD . dW;
                Switch[bm,
                    "None",
                        dCFid = zZero,
                    "Global",
                        dzSc = Tr[dzMat] / kN;
                        dCFid = ((2. * Conjugate[dzSc]) / kN) * fConst,
                    "Diagonal",
                        dzDiag = Diagonal[dzMat];
                        dzProd = Conjugate[dzDiag] * uD;
                        dCFid = (2. / kN) * If[hasIsom, v . dzProd . vD, dzProd]
                ];
                If[hasIsom && lW > 0.,
                    dCLeak = (2. * lW / kN) * (v . ConjugateTranspose[dW] . vD)
                ];
                ToPackedComplexMatrix[dCFid + dCLeak]
            ]
        ]
    ];

    If[!TrueQ[pfSpec["ParametricQ"]],
        Return[<|
            "LeafKind" -> "Unitary",
            "Evaluate" -> baseEvalFn,
            "EvaluateDirectional" -> baseDirFn,
            "ExactHessianSupportedQ" -> True,
            "UnsupportedSecondOrderReason" -> None,
            "PhaseFreedomSpec" -> None
        |>]
    ];

    evalWithPhiFn = With[
        {bEval = baseEvalFn,
         gLeft = (Lookup[#, {"Mat", "Evals", "Vec", "VecH"}] & /@ pfSpec["Left"]),
         gRight = (Lookup[#, {"Mat", "Evals", "Vec", "VecH"}] & /@ pfSpec["Right"]),
         nL = pfSpec["NumLeft"], nR = pfSpec["NumRight"],
         idMat = Developer`ToPackedArray[IdentityMatrix[dim, WorkingPrecision -> MachinePrecision] + 0. I]},
        Function[{uFinal, phiVec},
            Module[{evalSide, pLDag, dPL, pRDag, dPR, uEff, bRes, phiVal, cEffDag, bExtra,
                    cUDag, kLeft, kRight, gL, gR, gPhi},
                evalSide = Function[{genSeq, angles}, Module[{n = Length[genSeq], mats, pref, suff, derivs, a, g, v, vH, e},
                    If[n === 0,
                        {idMat, {}},
                        mats = Table[
                            g = genSeq[[a]]; v = g[[3]]; vH = g[[4]]; e = g[[2]];
                            v . (Exp[I * angles[[a]] * e] * vH),
                            {a, n}
                        ];
                        pref = Table[idMat, {n + 1}];
                        Do[pref[[a + 1]] = pref[[a]] . mats[[a]], {a, n}];
                        suff = Table[idMat, {n + 1}];
                        Do[suff[[a]] = mats[[a]] . suff[[a + 1]], {a, n, 1, -1}];
                        derivs = Table[
                            pref[[a]] . ((I * genSeq[[a, 1]]) . mats[[a]]) . suff[[a + 1]],
                            {a, n}
                        ];
                        {pref[[n + 1]], derivs}
                    ]
                ]];
                {pLDag, dPL} = evalSide[gLeft, If[nL > 0, phiVec[[1 ;; nL]], {}]];
                {pRDag, dPR} = evalSide[gRight, If[nR > 0, phiVec[[nL + 1 ;; nL + nR]], {}]];
                uEff = ToPackedComplexMatrix[pLDag . uFinal . pRDag];
                {phiVal, cEffDag, bExtra} = bEval[uEff];
                cUDag = ToPackedComplexMatrix[pRDag . cEffDag . pLDag];
                gL = If[nL > 0,
                    kLeft = Transpose[uFinal . pRDag . cEffDag];
                    Table[Re[Total[Flatten[kLeft * dPL[[a]]]]], {a, nL}],
                    {}
                ];
                gR = If[nR > 0,
                    kRight = Transpose[cEffDag . pLDag . uFinal];
                    Table[Re[Total[Flatten[kRight * dPR[[b]]]]], {b, nR}],
                    {}
                ];
                gPhi = Developer`ToPackedArray[N[Join[gL, gR]]];
                {phiVal, cUDag, Join[bExtra, <|"PhaseVector" -> phiVec, "PhaseGradient" -> gPhi|>], gPhi}
            ]
        ]
    ];

    pfOut = Join[pfSpec, <|"EvaluateWithPhases" -> evalWithPhiFn|>];
    <|
        "LeafKind" -> "Unitary",
        "PhaseFreedomSpec" -> pfOut,
        "EvaluateDirectional" -> None,
        "ExactHessianSupportedQ" -> False,
        "UnsupportedSecondOrderReason" -> "Parametric PhaseFreedom second derivatives are not yet enabled.",
        "Evaluate" -> With[{ewp = evalWithPhiFn, zeros = ConstantArray[0., pfSpec["NumPhases"]]},
            Function[uFinal,
                Module[{r = ewp[uFinal, zeros]}, {r[[1]], r[[2]], r[[3]]}]
            ]
        ]
    |>
];

(* Compile single leaf target for Hilbert space *)
compileHilbertLeaf[t:UnitaryTarget[uTar_, opts___Rule], basis_, dim_Integer] := Module[
    {oType, pfOpt, defMode, pfSpec, uMat},
    {oType, pfOpt} = {ObjectiveType, PhaseFreedom} /. {opts} /. Options[UnitaryTarget];
    defMode = Which[
        oType === "UnitaryFidelity" || oType === Automatic, "Global",
        oType === "UnitaryOverlap", "None",
        True, Message[ObjectiveTarget::invalspec, t]; Return[$Failed]
    ];
    uMat = toHilbertOperatorMatrix[uTar, basis, dim];
    If[uMat === $Failed, Message[StateTarget::baddim, uTar, dim]; Return[$Failed]];
    pfSpec = parsePhaseFreedomSpec[pfOpt, defMode, basis, dim];
    If[pfSpec === $Failed, Message[ObjectiveTarget::badsub, uTar, pfOpt]; Return[$Failed]];
    buildUnitaryLeafAssoc[uMat, None, dim, pfSpec, 0.]
];

compileHilbertLeaf[t:SubspaceUnitaryTarget[uSub_, subSpec_, opts___Rule], basis_, dim_Integer] := Module[
    {oType, pfOpt, lamL, phaseIns, vMat, kDim, uSubMat, defMode, pfSpec},
    {oType, pfOpt, lamL, phaseIns} =
        {ObjectiveType, PhaseFreedom, LeakageWeight, PhaseInsensitive} /. {opts} /. Options[SubspaceUnitaryTarget];
    If[!NumericQ[lamL] || !TrueQ[Re[lamL] >= 0.] || Im[N[lamL]] != 0,
        Message[ObjectiveTarget::badsub, subSpec, LeakageWeight -> lamL]; Return[$Failed]
    ];
    vMat = resolveSubspaceIsometry[subSpec, basis, dim];
    If[vMat === $Failed, Message[ObjectiveTarget::badsub, subSpec, pfOpt]; Return[$Failed]];
    kDim = Dimensions[vMat][[2]];
    uSubMat = Which[
        MatrixQ[uSub, NumericQ] && Dimensions[uSub] === {kDim, kDim}, ToPackedComplexMatrix[uSub],
        kDim === dim && OperatorQ[uSub], toHilbertOperatorMatrix[uSub, basis, dim],
        True, $Failed
    ];
    If[uSubMat === $Failed, Message[StateTarget::baddim, uSub, kDim]; Return[$Failed]];
    defMode = Which[
        oType === "UnitaryOverlap" || phaseIns === False, "None",
        True, "Global"
    ];
    pfSpec = parsePhaseFreedomSpec[pfOpt, defMode, basis, dim];
    If[pfSpec === $Failed, Message[ObjectiveTarget::badsub, subSpec, pfOpt]; Return[$Failed]];
    buildUnitaryLeafAssoc[uSubMat, vMat, dim, pfSpec, N[Re[lamL]]]
];

compileHilbertLeaf[t:StateTarget[Rule[init_, tar_], opts___Rule], basis_, dim_Integer] := Module[
    {oType, k1, k2, m1, m2, m1H, m2H, a1, a2, sRoot, sNorm,
     zeroDimMat = ToPackedComplexMatrix[ConstantArray[0., {dim, dim}]]},
    oType = ObjectiveType /. {opts} /. Options[StateTarget];
    If[isKetOrVecQ[init] && isKetOrVecQ[tar],
        k1 = toKetVector[init, basis, dim];
        k2 = toKetVector[tar, basis, dim];
        If[k1 === $Failed || k2 === $Failed || Norm[k1] == 0 || Norm[k2] == 0,
            Message[StateTarget::baddim, init -> tar, dim]; Return[$Failed]];
        k1 = ToPackedComplexVector[k1 / Norm[k1]];
        k2 = ToPackedComplexVector[k2 / Norm[k2]];
        If[oType === Automatic, oType = "StateFidelity"];
        Return[Switch[oType,
            "StateFidelity",
                m1 = ToPackedComplexMatrix[Outer[Times, k1, Conjugate[k2]]];
                With[{p0 = k1, pT = k2, mK = m1},
                    <|
                        "LeafKind" -> "PureKet",
                        "Psi0" -> p0,
                        "PsiTar" -> pT,
                        "PhaseInsensitive" -> True,
                        "Evaluate" -> Function[uFinal, Module[{psiT, z},
                            psiT = uFinal . p0;
                            z = Conjugate[pT] . psiT;
                            {Re[z * Conjugate[z]], ToPackedComplexMatrix[(2. Conjugate[z]) * mK]}
                        ]],
                        "EvaluateDirectional" -> Function[{uFinal, dU},
                            ToPackedComplexMatrix[(2. * Conjugate[Conjugate[pT] . (dU . p0)]) * mK]
                        ],
                        "ExactHessianSupportedQ" -> True,
                        "UnsupportedSecondOrderReason" -> None
                    |>
                ],
            "StateOverlap",
                m1 = ToPackedComplexMatrix[Outer[Times, k1, Conjugate[k2]]];
                With[{p0 = k1, pT = k2, mK = m1, zZero = zeroDimMat},
                    <|
                        "LeafKind" -> "PureKet",
                        "Psi0" -> p0,
                        "PsiTar" -> pT,
                        "PhaseInsensitive" -> False,
                        "Evaluate" -> Function[uFinal, Module[{psiT, z},
                            psiT = uFinal . p0;
                            z = Conjugate[pT] . psiT;
                            {Re[z], mK}
                        ]],
                        "EvaluateDirectional" -> Function[{uFinal, dU}, zZero],
                        "ExactHessianSupportedQ" -> True,
                        "UnsupportedSecondOrderReason" -> None
                    |>
                ],
            "StateDistance",
                m1 = ToPackedComplexMatrix[Outer[Times, k1, Conjugate[k1]]];
                m2 = ToPackedComplexMatrix[Outer[Times, k2, Conjugate[k2]]];
                With[{mat1 = m1, mat2 = m2},
                    <|
                        "LeafKind" -> "Density",
                        "Evaluate" -> Function[uFinal, Module[{uH, rhoT, val},
                            uH = ConjugateTranspose[uFinal];
                            rhoT = uFinal . mat1 . uH;
                            val = Re[Tr[mat2 . rhoT]];
                            {val, ToPackedComplexMatrix[2. * (mat1 . uH . mat2)]}
                        ]],
                        "EvaluateDirectional" -> Function[{uFinal, dU},
                            ToPackedComplexMatrix[2. * (mat1 . ConjugateTranspose[dU] . mat2)]
                        ],
                        "ExactHessianSupportedQ" -> True,
                        "UnsupportedSecondOrderReason" -> None
                    |>
                ],
            _,
                Message[ObjectiveTarget::invalspec, t]; $Failed
        ]]
    ];
    m1 = toHilbertOperatorMatrix[init, basis, dim];
    m2 = toHilbertOperatorMatrix[tar, basis, dim];
    If[m1 === $Failed || m2 === $Failed, Message[StateTarget::baddim, init -> tar, dim]; Return[$Failed]];
    If[oType === Automatic,
        a1 = psdAnalyzeMatrix[m1];
        a2 = psdAnalyzeMatrix[m2];
        If[a1["ValidQ"] && a2["ValidQ"],
            oType = "StateFidelity",
            oType = "ObservableTransfer"
        ]
    ];
    Switch[oType,
        "StateFidelity",
            If[!AssociationQ[a1], a1 = psdAnalyzeMatrix[m1]];
            If[!AssociationQ[a2], a2 = psdAnalyzeMatrix[m2]];
            If[!a1["ValidQ"] || !a2["ValidQ"],
                Message[StateTarget::notpsd, init -> tar]; Return[$Failed]
            ];
            If[a1["PureQ"] && a2["PureQ"],
                k1 = a1["DominantKet"]; k2 = a2["DominantKet"];
                m1 = ToPackedComplexMatrix[Outer[Times, k1, Conjugate[k2]]];
                Return[With[{p0 = k1, pT = k2, mK = m1},
                    <|
                        "LeafKind" -> "PureKet",
                        "Psi0" -> p0,
                        "PsiTar" -> pT,
                        "PhaseInsensitive" -> True,
                        "Evaluate" -> Function[uFinal, Module[{psiT, z},
                            psiT = uFinal . p0;
                            z = Conjugate[pT] . psiT;
                            {Re[z * Conjugate[z]], ToPackedComplexMatrix[(2. Conjugate[z]) * mK]}
                        ]],
                        "EvaluateDirectional" -> Function[{uFinal, dU},
                            ToPackedComplexMatrix[(2. * Conjugate[Conjugate[pT] . (dU . p0)]) * mK]
                        ],
                        "ExactHessianSupportedQ" -> True,
                        "UnsupportedSecondOrderReason" -> None
                    |>
                ]]
            ];
            m1 = a1["DensityMatrix"]; m2 = a2["DensityMatrix"];
            If[a1["PureQ"] || a2["PureQ"],
                Return[With[{mat1 = m1, mat2 = m2},
                    <|
                        "LeafKind" -> "Density",
                        "Evaluate" -> Function[uFinal, Module[{uH, rhoT, val},
                            uH = ConjugateTranspose[uFinal];
                            rhoT = uFinal . mat1 . uH;
                            val = Re[Tr[mat2 . rhoT]];
                            {val, ToPackedComplexMatrix[2. * (mat1 . uH . mat2)]}
                        ]],
                        "EvaluateDirectional" -> Function[{uFinal, dU},
                            ToPackedComplexMatrix[2. * (mat1 . ConjugateTranspose[dU] . mat2)]
                        ],
                        "ExactHessianSupportedQ" -> True,
                        "UnsupportedSecondOrderReason" -> None
                    |>
                ]]
            ];
            sRoot = a2["RootMatrix"];
            With[{mat1 = m1, sr = sRoot},
                <|
                    "LeafKind" -> "MixedDensity",
                    "EvaluateDirectional" -> None,
                    "ExactHessianSupportedQ" -> False,
                    "UnsupportedSecondOrderReason" -> "Mixed-state Uhlmann-Jozsa fidelity second derivatives are not yet supported.",
                    "Evaluate" -> Function[uFinal, Module[{uH, rhoT, mMat, evM, vecM, vM, mu, fVal, mInvH, gRho},
                        uH = ConjugateTranspose[uFinal];
                        rhoT = uFinal . mat1 . uH;
                        mMat = HermitianSymmetrize[sr . rhoT . sr];
                        {evM, vecM} = Eigensystem[mMat];
                        vM = Transpose[vecM];
                        mu = Max[#, 10^-14] & /@ Re[evM];
                        fVal = Total[Sqrt[mu]];
                        mInvH = vM . ((1. / Sqrt[mu]) * ConjugateTranspose[vM]);
                        gRho = fVal * (sr . mInvH . sr);
                        {fVal^2, ToPackedComplexMatrix[2. * (mat1 . uH . gRho)]}
                    ]]
                |>
            ],
        "StateOverlap" | "HilbertSchmidt",
            sNorm = resolveScaleFactor["Frobenius", m1, m2];
            If[sNorm === $Failed, Return[$Failed]];
            m1H = ToPackedComplexMatrix[ConjugateTranspose[m1]];
            m2H = ToPackedComplexMatrix[ConjugateTranspose[m2]];
            With[{mat1 = m1, mat2 = m2, mat1H = m1H, mat2H = m2H, sn = sNorm},
                <|
                    "LeafKind" -> "Observable",
                    "Evaluate" -> Function[uFinal, Module[{uH, rhoT, val},
                        uH = ConjugateTranspose[uFinal];
                        rhoT = uFinal . mat1 . uH;
                        val = Re[Tr[mat2H . rhoT]] / sn;
                        {val, ToPackedComplexMatrix[(mat1 . uH . mat2H + mat1H . uH . mat2) / sn]}
                    ]],
                    "EvaluateDirectional" -> Function[{uFinal, dU},
                        With[{dUH = ConjugateTranspose[dU]},
                            ToPackedComplexMatrix[(mat1 . dUH . mat2H + mat1H . dUH . mat2) / sn]
                        ]
                    ],
                    "ExactHessianSupportedQ" -> True,
                    "UnsupportedSecondOrderReason" -> None
                |>
            ],
        "StateDistance",
            sNorm = 0.5 * (Total[Abs[Flatten[m1]]^2] + Total[Abs[Flatten[m2]]^2]);
            If[sNorm <= 0., Return[$Failed]];
            m1H = ToPackedComplexMatrix[ConjugateTranspose[m1]];
            m2H = ToPackedComplexMatrix[ConjugateTranspose[m2]];
            With[{mat1 = m1, mat2 = m2, mat1H = m1H, mat2H = m2H, sn = sNorm},
                <|
                    "LeafKind" -> "Observable",
                    "Evaluate" -> Function[uFinal, Module[{uH, rhoT, val},
                        uH = ConjugateTranspose[uFinal];
                        rhoT = uFinal . mat1 . uH;
                        val = 1. - Total[Abs[Flatten[rhoT - mat2]]^2] / (2. * sn);
                        {val, ToPackedComplexMatrix[(mat1 . uH . mat2H + mat1H . uH . mat2) / sn]}
                    ]],
                    "EvaluateDirectional" -> Function[{uFinal, dU},
                        With[{dUH = ConjugateTranspose[dU]},
                            ToPackedComplexMatrix[(mat1 . dUH . mat2H + mat1H . dUH . mat2) / sn]
                        ]
                    ],
                    "ExactHessianSupportedQ" -> True,
                    "UnsupportedSecondOrderReason" -> None
                |>
            ],
        "ObservableTransfer",
            compileHilbertLeaf[ObservableTarget[init -> tar, ObjectiveType -> "ObservableTransfer", TargetNormalization -> "Frobenius"], basis, dim],
        _,
            Message[ObjectiveTarget::invalspec, t]; $Failed
    ]
];

compileHilbertLeaf[t:ObservableTarget[Rule[init_, obs_], opts___Rule], basis_, dim_Integer] := Module[
    {oType, normSpec, phaseIns, m1, m2, m1H, m2H, sNorm},
    {oType, normSpec, phaseIns} = {ObjectiveType, TargetNormalization, PhaseInsensitive} /. {opts} /. Options[ObservableTarget];
    If[TrueQ[phaseIns] && oType === "ObservableTransfer", oType = "ObservableMagnitude"];
    m1 = toHilbertOperatorMatrix[init, basis, dim];
    m2 = toHilbertOperatorMatrix[obs, basis, dim];
    If[m1 === $Failed || m2 === $Failed, Message[StateTarget::baddim, init -> obs, dim]; Return[$Failed]];
    sNorm = resolveScaleFactor[normSpec, m1, m2];
    If[sNorm === $Failed, Return[$Failed]];
    m1H = ToPackedComplexMatrix[ConjugateTranspose[m1]];
    m2H = ToPackedComplexMatrix[ConjugateTranspose[m2]];
    Switch[oType,
        "ObservableTransfer" | "CoherenceTransfer" | Automatic,
            With[{mat1 = m1, mat2 = m2, mat1H = m1H, mat2H = m2H, sn = sNorm},
                <|
                    "LeafKind" -> "Observable",
                    "Evaluate" -> Function[uFinal, Module[{uH, rhoT, val},
                        uH = ConjugateTranspose[uFinal];
                        rhoT = uFinal . mat1 . uH;
                        val = Re[Tr[mat2H . rhoT]] / sn;
                        {val, ToPackedComplexMatrix[(mat1 . uH . mat2H + mat1H . uH . mat2) / sn]}
                    ]],
                    "EvaluateDirectional" -> Function[{uFinal, dU},
                        With[{dUH = ConjugateTranspose[dU]},
                            ToPackedComplexMatrix[(mat1 . dUH . mat2H + mat1H . dUH . mat2) / sn]
                        ]
                    ],
                    "ExactHessianSupportedQ" -> True,
                    "UnsupportedSecondOrderReason" -> None
                |>
            ],
        "ObservableMagnitude",
            With[{mat1 = m1, mat2 = m2, mat1H = m1H, mat2H = m2H, sn = sNorm},
                <|
                    "LeafKind" -> "Observable",
                    "Evaluate" -> Function[uFinal, Module[{uH, rhoT, z},
                        uH = ConjugateTranspose[uFinal];
                        rhoT = uFinal . mat1 . uH;
                        z = Tr[mat2H . rhoT] / sn;
                        {Re[z * Conjugate[z]],
                         ToPackedComplexMatrix[(2. / sn) * (Conjugate[z] * (mat1 . uH . mat2H) + z * (mat1H . uH . mat2))]}
                    ]],
                    "EvaluateDirectional" -> Function[{uFinal, dU},
                        Module[{uH = ConjugateTranspose[uFinal], dUH = ConjugateTranspose[dU], rhoT, dRhoT, z, dz},
                            rhoT = uFinal . mat1 . uH;
                            dRhoT = dU . mat1 . uH + uFinal . mat1 . dUH;
                            z = Tr[mat2H . rhoT] / sn;
                            dz = Tr[mat2H . dRhoT] / sn;
                            ToPackedComplexMatrix[(2. / sn) * (
                                Conjugate[dz] * (mat1 . uH . mat2H) + Conjugate[z] * (mat1 . dUH . mat2H) +
                                dz * (mat1H . uH . mat2) + z * (mat1H . dUH . mat2)
                            )]
                        ]
                    ],
                    "ExactHessianSupportedQ" -> True,
                    "UnsupportedSecondOrderReason" -> None
                |>
            ],
        _,
            Message[ObjectiveTarget::invalspec, t]; $Failed
    ]
];

compileHilbertLeaf[CoherenceTarget[Rule[init_, coh_], opts___Rule], basis_, dim_Integer] := Module[
    {oType, normSpec, phaseIns},
    {oType, normSpec, phaseIns} = {ObjectiveType, TargetNormalization, PhaseInsensitive} /. {opts} /. Options[CoherenceTarget];
    If[oType === Automatic,
        oType = If[TrueQ[phaseIns], "ObservableMagnitude", "ObservableTransfer"]
    ];
    compileHilbertLeaf[
        ObservableTarget[init -> coh, ObjectiveType -> oType, TargetNormalization -> normSpec, PhaseInsensitive -> phaseIns],
        basis, dim
    ]
];

compileHilbertLeaf[other_, _, _] := (Message[ObjectiveTarget::invalspec, other]; $Failed);

(* Compile Liouville leaf target specification *)
compileLiouvilleLeaf[t:StateTarget[Rule[init_, tar_], opts___Rule], opbasis_, basis_, dl_Integer] := Module[
    {oType, r1, r2, sNorm, mode},
    oType = ObjectiveType /. {opts} /. Options[StateTarget];
    r1 = toLiouvilleVector[init, opbasis, basis, dl];
    r2 = toLiouvilleVector[tar, opbasis, basis, dl];
    If[r1 === $Failed || r2 === $Failed, Message[StateTarget::baddim, init -> tar, dl]; Return[$Failed]];
    If[oType === Automatic, oType = "ObservableTransfer"];
    {mode, sNorm} = Switch[oType,
        "StateFidelity" | "ObservableTransfer" | "StateOverlap" | "HilbertSchmidt",
            {"Linear", resolveScaleFactor["Frobenius", r1, r2]},
        "StateDistance",
            {"Distance", Total[Abs[r1]^2] + Total[Abs[r2]^2]},
        _,
            {"Linear", resolveScaleFactor["Frobenius", r1, r2]}
    ];
    If[sNorm === $Failed || sNorm <= 0., Return[$Failed]];
    <|"Rho0" -> r1, "Obs" -> r2, "Mode" -> mode, "Scale" -> sNorm|>
];

compileLiouvilleLeaf[t:ObservableTarget[Rule[init_, obs_], opts___Rule], opbasis_, basis_, dl_Integer] := Module[
    {oType, normSpec, phaseIns, r1, r2, sNorm, mode},
    {oType, normSpec, phaseIns} = {ObjectiveType, TargetNormalization, PhaseInsensitive} /. {opts} /. Options[ObservableTarget];
    If[TrueQ[phaseIns] && oType === "ObservableTransfer", oType = "ObservableMagnitude"];
    r1 = toLiouvilleVector[init, opbasis, basis, dl];
    r2 = toLiouvilleVector[obs, opbasis, basis, dl];
    If[r1 === $Failed || r2 === $Failed, Message[StateTarget::baddim, init -> obs, dl]; Return[$Failed]];
    sNorm = resolveScaleFactor[normSpec, r1, r2];
    If[sNorm === $Failed, Return[$Failed]];
    mode = If[oType === "ObservableMagnitude", "Magnitude", "Linear"];
    <|"Rho0" -> r1, "Obs" -> r2, "Mode" -> mode, "Scale" -> sNorm|>
];

compileLiouvilleLeaf[CoherenceTarget[Rule[init_, coh_], opts___Rule], opbasis_, basis_, dl_Integer] := Module[
    {oType, normSpec, phaseIns},
    {oType, normSpec, phaseIns} = {ObjectiveType, TargetNormalization, PhaseInsensitive} /. {opts} /. Options[CoherenceTarget];
    If[oType === Automatic,
        oType = If[TrueQ[phaseIns], "ObservableMagnitude", "ObservableTransfer"]
    ];
    compileLiouvilleLeaf[
        ObservableTarget[init -> coh, ObjectiveType -> oType, TargetNormalization -> normSpec, PhaseInsensitive -> phaseIns],
        opbasis, basis, dl
    ]
];

compileLiouvilleLeaf[other_, _, _, _] := (Message[ObjectiveTarget::invalspec, other]; $Failed);

(* Main CompileObjective dispatcher *)
extractPrimaryInitSpec[targetSpec_] := Which[
    MatchQ[targetSpec, _StateTarget | _ObservableTarget | _CoherenceTarget] && MatchQ[targetSpec[[1]], _Rule],
        targetSpec[[1, 1]],
    MatchQ[targetSpec, _CooperativeObjective],
        Module[{leaves = targetSpec[[1]], firstLeaf, inits},
            firstLeaf = First[leaves];
            If[MatchQ[firstLeaf, {_?NumericQ, _}], firstLeaf = firstLeaf[[2]]];
            If[MatchQ[firstLeaf, _StateTarget | _ObservableTarget | _CoherenceTarget] && MatchQ[firstLeaf[[1]], _Rule],
                inits = DeleteDuplicates[Map[
                    With[{el = If[MatchQ[#, {_?NumericQ, _}], #[[2]], #]},
                        If[MatchQ[el, _StateTarget | _ObservableTarget | _CoherenceTarget] && MatchQ[el[[1]], _Rule], el[[1, 1]], None]
                    ] &,
                    leaves
                ]];
                If[Length[inits] === 1 && First[inits] =!= None, First[inits], None],
                None
            ]
        ],
    True, None
];

CompileObjective[targetSpec_, space_String, basis_, opbasis_, opts___Rule] := Module[
    {normSpec, normW, subEntries, rawW, rawTargets, weights, dim, dl, compiledLeaves,
     allPureKets, psiMat0, tarMat, phaseInsList, evalFns, primaryInit},
    normSpec = NormalizeTargetSpec[targetSpec];
    If[normSpec === $Failed, Return[$Failed]];
    primaryInit = extractPrimaryInitSpec[normSpec];
    If[MatchQ[normSpec, _CooperativeObjective],
        normW = NormalizeWeights /. (List @@ Rest[normSpec]) /. Options[CooperativeObjective];
        subEntries = First[normSpec];
        If[!ListQ[subEntries] || subEntries === {},
            Message[ObjectiveTarget::invalspec, targetSpec]; Return[$Failed]
        ];
        {rawW, rawTargets} = Transpose[Map[
            Function[item,
                If[MatchQ[item, {_?NumericQ, _}],
                    {N[First[item]], NormalizeTargetSpec[item[[2]]]},
                    {1., NormalizeTargetSpec[item]}
                ]
            ],
            subEntries
        ]];
        If[!FreeQ[rawTargets, $Failed], Return[$Failed]];
        weights = If[TrueQ[normW],
            If[Total[rawW] == 0, Return[$Failed], N[rawW / Total[rawW]]],
            N[rawW]
        ],
        rawTargets = {normSpec};
        weights = {1.}
    ];
    If[space === "Hilbert",
        dim = BasisDimension[basis];
        compiledLeaves = compileHilbertLeaf[#, basis, dim] & /@ rawTargets;
        If[MemberQ[compiledLeaves, $Failed], Return[$Failed]];
        allPureKets = AllTrue[compiledLeaves, #["LeafKind"] === "PureKet" &];
        evalFns = #["Evaluate"] & /@ compiledLeaves;
        If[Length[compiledLeaves] === 1,
            With[{evFn = First[evalFns], pf = Lookup[First[compiledLeaves], "PhaseFreedomSpec", None],
                  dirFn = Lookup[First[compiledLeaves], "EvaluateDirectional", None],
                  exactQ = TrueQ[Lookup[First[compiledLeaves], "ExactHessianSupportedQ", False]],
                  badReason = Lookup[First[compiledLeaves], "UnsupportedSecondOrderReason", None]},
                <|
                    "Kind" -> compiledLeaves[[1, "LeafKind"]],
                    "Space" -> "Hilbert",
                    "Dimension" -> dim,
                    "PureKetFastPath" -> If[allPureKets,
                        <|
                            "PsiMat0" -> ToPackedComplexMatrix[Transpose[{compiledLeaves[[1, "Psi0"]]}]],
                            "TarMat" -> ToPackedComplexMatrix[Transpose[{compiledLeaves[[1, "PsiTar"]]}]],
                            "Weights" -> Developer`ToPackedArray[{1.}],
                            "PhaseInsensitive" -> {compiledLeaves[[1, "PhaseInsensitive"]]}
                        |>,
                        None
                    ],
                    "PrimaryInitialState" -> primaryInit,
                    "PhaseFreedomSpec" -> pf,
                    "ExactHessianSupportedQ" -> exactQ,
                    "UnsupportedSecondOrderReason" -> badReason,
                    "EvaluateCotangentDirectional" -> dirFn,
                    "EvaluateCotangent" -> If[AssociationQ[pf],
                        Function[uFinal,
                            With[{res = evFn[uFinal]}, {res[[1]], res[[2]], Join[<|"SubFidelities" -> {res[[1]]}|>, res[[3]]]}]
                        ],
                        Function[uFinal,
                            With[{res = evFn[uFinal]},
                                {res[[1]], res[[2]],
                                 If[Length[res] >= 3 && AssociationQ[res[[3]]],
                                     Join[<|"SubFidelities" -> {res[[1]]}|>, res[[3]]],
                                     <|"SubFidelities" -> {res[[1]]}|>
                                 ]}
                            ]
                        ]
                    ],
                    "RawTarget" -> normSpec
                |>
            ],
            If[AnyTrue[compiledLeaves, AssociationQ[Lookup[#, "PhaseFreedomSpec", None]] &],
                Message[ObjectiveTarget::badsub, normSpec, "Parametric PhaseFreedom cannot be nested in CooperativeObjective"];
                Return[$Failed]
            ];
            With[{fns = evalFns, wVec = Developer`ToPackedArray[weights], nSub = Length[evalFns], dVal = dim,
                  dirFns = Lookup[#, "EvaluateDirectional", None] & /@ compiledLeaves,
                  exactQ = AllTrue[compiledLeaves, TrueQ[Lookup[#, "ExactHessianSupportedQ", False]] &],
                  badReason = First[Select[Lookup[#, "UnsupportedSecondOrderReason", None] & /@ compiledLeaves, # =!= None &],
                      "One or more leaves in CooperativeObjective do not support exact second derivatives."]},
                <|
                    "Kind" -> If[allPureKets, "PureKets", "Cooperative"],
                    "Space" -> "Hilbert",
                    "Dimension" -> dim,
                    "PureKetFastPath" -> If[allPureKets,
                        <|
                            "PsiMat0" -> ToPackedComplexMatrix[Transpose[#["Psi0"] & /@ compiledLeaves]],
                            "TarMat" -> ToPackedComplexMatrix[Transpose[#["PsiTar"] & /@ compiledLeaves]],
                            "Weights" -> wVec,
                            "PhaseInsensitive" -> (#["PhaseInsensitive"] & /@ compiledLeaves)
                        |>,
                        None
                    ],
                    "PrimaryInitialState" -> primaryInit,
                    "PhaseFreedomSpec" -> None,
                    "ExactHessianSupportedQ" -> exactQ,
                    "UnsupportedSecondOrderReason" -> If[exactQ, None, badReason],
                    "EvaluateCotangentDirectional" -> If[exactQ,
                        With[{dfns = dirFns, zZero = ToPackedComplexMatrix[ConstantArray[0., {dVal, dVal}]]},
                            Function[{uFinal, dU}, Module[{acc = zZero, idx},
                                Do[acc = acc + wVec[[idx]] * (dfns[[idx]])[uFinal, dU], {idx, nSub}];
                                ToPackedComplexMatrix[acc]
                            ]]
                        ],
                        None
                    ],
                    "EvaluateCotangent" -> Function[uFinal, Module[{vals, cSum, idx, pair, extraOut = <||>},
                        vals = ConstantArray[0., nSub];
                        cSum = ToPackedComplexMatrix[ConstantArray[0., {dVal, dVal}]];
                        Do[
                            pair = (fns[[idx]])[uFinal];
                            vals[[idx]] = pair[[1]];
                            cSum = cSum + wVec[[idx]] * pair[[2]];
                            If[Length[pair] >= 3 && AssociationQ[pair[[3]]],
                                extraOut = Join[extraOut, pair[[3]]]
                            ],
                            {idx, nSub}
                        ];
                        {wVec . vals, ToPackedComplexMatrix[cSum], Join[extraOut, <|"SubFidelities" -> vals|>]}
                    ]],
                    "RawTarget" -> normSpec
                |>
            ]
        ],
        (* Liouville space *)
        dl = OperatorBasisDimension[opbasis];
        compiledLeaves = compileLiouvilleLeaf[#, opbasis, basis, dl] & /@ rawTargets;
        If[!FreeQ[compiledLeaves, $Failed], Return[$Failed]];
        <|
            "Kind" -> "Liouville",
            "Space" -> "Liouville",
            "Dimension" -> dl,
            "InitialMatrix" -> ToPackedComplexMatrix[Transpose[#["Rho0"] & /@ compiledLeaves]],
            "TargetMatrix" -> ToPackedComplexMatrix[Transpose[#["Obs"] & /@ compiledLeaves]],
            "Weights" -> Developer`ToPackedArray[N[weights]],
            "Modes" -> (#["Mode"] & /@ compiledLeaves),
            "Scales" -> Developer`ToPackedArray[N[#["Scale"] & /@ compiledLeaves]],
            "ExactHessianSupportedQ" -> True,
            "UnsupportedSecondOrderReason" -> None,
            "PrimaryInitialState" -> primaryInit,
            "RawTarget" -> normSpec
        |>
    ]
];

(* Intermediate-time RunningCost & StatePenalty compilation *)
buildStaticTimeWeights[twSpec_, nt_Integer] := Module[{raw},
    Which[
        twSpec === "Uniform" || twSpec === Automatic,
            Developer`ToPackedArray[ConstantArray[1. / nt, nt]],
        twSpec === "FinalHalf",
            raw = Table[If[nt === 1 || j > Floor[nt / 2], 1., 0.], {j, nt}];
            Developer`ToPackedArray[raw / Max[Total[raw], 1.]],
        twSpec === "LinearRamp",
            raw = N[Range[nt]];
            Developer`ToPackedArray[raw / Total[raw]],
        ListQ[twSpec] && Length[twSpec] === nt && VectorQ[twSpec, NumericQ] && AllTrue[twSpec, TrueQ[# >= 0.] &],
            raw = N[Re[twSpec]];
            Developer`ToPackedArray[If[Total[raw] > 0., raw / Total[raw], ConstantArray[0., nt]]],
        True,
            None
    ]
];

makeTimeWeightEvaluator[twSpec_, nt_Integer] := Module[{staticW},
    staticW = buildStaticTimeWeights[twSpec, nt];
    If[VectorQ[staticW, NumericQ],
        Return[<|"StaticQ" -> True, "Weights" -> staticW, "EvalRaw" -> None|>]
    ];
    If[MatchQ[twSpec, _Symbol | _Function] && !MemberQ[{None, False, True}, twSpec],
        Return[<|"StaticQ" -> False, "Weights" -> None, "EvalRaw" -> With[{fn = twSpec},
            Function[{dtVec, t0},
                Module[{tEnds = N[t0] + Accumulate[N[dtVec]], raw},
                    raw = Table[
                        With[{v = Quiet[Check[N[Re[fn[tEnds[[j]], j]]], Quiet[Check[N[Re[fn[tEnds[[j]]]]], $Failed]]]]},
                            If[!NumericQ[v], $Failed, v]
                        ],
                        {j, Length[dtVec]}
                    ];
                    If[!FreeQ[raw, $Failed] || !AllTrue[raw, TrueQ[# >= 0.] &], $Failed,
                        Developer`ToPackedArray[raw]
                    ]
                ]
            ]
        ]|>]
    ];
    $Failed
];

CompileRunningCost[runSpec_, compiledObj_Association, space_String, basis_, opbasis_, dtVec_List, t0_: 0.] := Module[
    {nt = Length[dtVec], dim = compiledObj["Dimension"], items, compiledItems,
     svOK = True, hasTimeDep = False, termEvalFns, termSvFns, termLiouvFns, twSpecs},
    If[runSpec === None || runSpec === False || runSpec === {}, Return[None]];
    items = If[ListQ[runSpec], runSpec, {runSpec}];
    If[items === {}, Return[None]];

    compiledItems = Table[
        With[{item = items[[idx]]},
            If[MatchQ[item, _StatePenalty],
                Module[{spec, weight, optsList, modeOpt, twOpt, signOpt, normOpt, initOpt,
                        initExpr, tarExpr, coeff, twComp, m1, k1, m2, k2, m1H, m2H, hermErr,
                        effMode, a1, a2, sn, mOuter, vInit, vTar, rho0Cols, colMatches, rCol,
                        svCompatible = False, pkFast, svColIdx = None,
                        evalHilb, evalSV, evalLiouv},
                    If[Length[item] < 1, Message[StatePenalty::badspec, item]; Return[$Failed, Table]];
                    spec = item[[1]];
                    {weight, optsList} = If[Length[item] >= 2 && !MatchQ[item[[2]], _Rule],
                        {item[[2]], List @@ (item[[3 ;;]])},
                        {1., List @@ (item[[2 ;;]])}
                    ];
                    If[!NumericQ[weight] || !TrueQ[Re[weight] >= 0.] || Im[N[weight]] != 0,
                        Message[StatePenalty::badspec, item]; Return[$Failed, Table]
                    ];
                    {modeOpt, twOpt, signOpt, normOpt, initOpt} =
                        {TargetMode, TimeWeighting, Sign, TargetNormalization, InitialState} /. optsList /. Options[StatePenalty];
                    If[!NumericQ[signOpt] || Im[N[signOpt]] != 0,
                        Message[StatePenalty::badspec, item]; Return[$Failed, Table]
                    ];
                    coeff = N[Re[signOpt]] * N[Re[weight]];
                    twComp = makeTimeWeightEvaluator[twOpt, nt];
                    If[twComp === $Failed, Message[StatePenalty::badspec, item]; Return[$Failed, Table]];
                    If[!twComp["StaticQ"], hasTimeDep = True];

                    If[MatchQ[spec, _Rule],
                        {initExpr, tarExpr} = {spec[[1]], spec[[2]]},
                        tarExpr = spec;
                        initExpr = If[initOpt =!= Automatic, initOpt, Lookup[compiledObj, "PrimaryInitialState", None]]
                    ];
                    If[initExpr === None, Message[StatePenalty::badspec, item]; Return[$Failed, Table]];

                    If[space === "Hilbert",
                        (* Resolve initial state m1 (and optional pure ket k1) *)
                        k1 = If[isKetOrVecQ[initExpr],
                            With[{kv = toKetVector[initExpr, basis, dim]},
                                If[kv === $Failed || Norm[kv] == 0, Message[StateTarget::baddim, initExpr, dim]; Return[$Failed, Table]];
                                ToPackedComplexVector[kv / Norm[kv]]
                            ],
                            None
                        ];
                        m1 = If[k1 =!= None,
                            ToPackedComplexMatrix[Outer[Times, k1, Conjugate[k1]]],
                            With[{mat1 = toHilbertOperatorMatrix[initExpr, basis, dim]},
                                If[mat1 === $Failed, Message[StateTarget::baddim, initExpr, dim]; Return[$Failed, Table]];
                                mat1
                            ]
                        ];
                        (* Resolve target operator or ket m2 (and optional pure ket k2) *)
                        k2 = If[isKetOrVecQ[tarExpr],
                            With[{kv = toKetVector[tarExpr, basis, dim]},
                                If[kv === $Failed || Norm[kv] == 0, Message[StateTarget::baddim, tarExpr, dim]; Return[$Failed, Table]];
                                ToPackedComplexVector[kv / Norm[kv]]
                            ],
                            None
                        ];
                        m2 = If[k2 =!= None,
                            ToPackedComplexMatrix[Outer[Times, k2, Conjugate[k2]]],
                            With[{mat2 = toHilbertOperatorMatrix[tarExpr, basis, dim]},
                                If[mat2 === $Failed, Message[StateTarget::baddim, tarExpr, dim]; Return[$Failed, Table]];
                                mat2
                            ]
                        ];

                        If[k1 === None,
                            a1 = psdAnalyzeMatrix[m1];
                            If[a1["ValidQ"] && a1["PureQ"], k1 = a1["DominantKet"]]
                        ];
                        If[k2 === None,
                            a2 = psdAnalyzeMatrix[m2];
                            If[a2["ValidQ"] && a2["PureQ"], k2 = a2["DominantKet"]]
                        ];

                        hermErr = Max[Abs[m2 - ConjugateTranspose[m2]]] / Max[1., Max[Abs[m2]]];
                        effMode = Which[
                            modeOpt === Automatic,
                                If[k2 =!= None && k1 =!= None,
                                    "StateFidelity",
                                    If[hermErr > 10^-8, "Magnitude", "Expectation"]
                                ],
                            MemberQ[{"Expectation", "Magnitude", "StateFidelity", "Projector"}, modeOpt],
                                modeOpt,
                            True,
                                Message[StatePenalty::badspec, item]; Return[$Failed, Table]
                        ];

                        If[effMode === "StateFidelity" && k1 =!= None && k2 =!= None,
                            m1 = ToPackedComplexMatrix[Outer[Times, k1, Conjugate[k1]]];
                            m2 = ToPackedComplexMatrix[Outer[Times, k2, Conjugate[k2]]];
                            mOuter = ToPackedComplexMatrix[Outer[Times, k1, Conjugate[k2]]];
                            evalHilb = With[{cF = coeff, p0 = k1, pT = k2, mK = mOuter, dV = dim},
                                Function[{uJ, tJ, j, wJ},
                                    Module[{z = Conjugate[pT] . (uJ . p0), s = cF * wJ},
                                        {s * Re[z * Conjugate[z]], ToPackedComplexMatrix[(2. * s * Conjugate[z]) * mK]}
                                    ]
                                ]
                            ],
                            (* Operator Expectation / Magnitude / Projector in Hilbert Prop mode *)
                            If[effMode === "StateFidelity",
                                If[!AssociationQ[a1], a1 = psdAnalyzeMatrix[m1]];
                                If[!AssociationQ[a2], a2 = psdAnalyzeMatrix[m2]];
                                If[!a1["ValidQ"] || !a2["ValidQ"], Message[StateTarget::notpsd, item]; Return[$Failed, Table]];
                                m1 = a1["DensityMatrix"];
                                m2 = a2["DensityMatrix"]
                            ];
                            sn = Which[
                                effMode === "StateFidelity", 1.,
                                effMode === "Projector",
                                    If[normOpt === Automatic || normOpt === "None" || normOpt === False, 1.,
                                        resolveScaleFactor[normOpt, m1, m2]],
                                True,
                                    resolveScaleFactor[If[normOpt === Automatic, "Frobenius", normOpt], m1, m2]
                            ];
                            If[sn === $Failed || sn <= 0., Return[$Failed, Table]];
                            m1H = ToPackedComplexMatrix[ConjugateTranspose[m1]];
                            m2H = ToPackedComplexMatrix[ConjugateTranspose[m2]];
                            evalHilb = If[effMode === "Magnitude",
                                With[{cF = coeff, mat1 = m1, mat2 = m2, mat1H = m1H, mat2H = m2H, scaleN = sn},
                                    Function[{uJ, tJ, j, wJ},
                                        Module[{uH = ConjugateTranspose[uJ], rhoT, z, s = cF * wJ},
                                            rhoT = uJ . mat1 . uH;
                                            z = Tr[mat2H . rhoT] / scaleN;
                                            {s * Re[z * Conjugate[z]],
                                             ToPackedComplexMatrix[(2. * s / scaleN) * (Conjugate[z] * (mat1 . uH . mat2H) + z * (mat1H . uH . mat2))]}
                                        ]
                                    ]
                                ],
                                With[{cF = coeff, mat1 = m1, mat2 = m2, mat1H = m1H, mat2H = m2H, scaleN = sn},
                                    Function[{uJ, tJ, j, wJ},
                                        Module[{uH = ConjugateTranspose[uJ], rhoT, val, s = cF * wJ},
                                            rhoT = uJ . mat1 . uH;
                                            val = Re[Tr[mat2H . rhoT]] / scaleN;
                                            {s * val, ToPackedComplexMatrix[(s / scaleN) * (mat1 . uH . mat2H + mat1H . uH . mat2)]}
                                        ]
                                    ]
                                ]
                            ]
                        ];

                        (* Check StateVector fast-path compatibility *)
                        pkFast = compiledObj["PureKetFastPath"];
                        If[AssociationQ[pkFast] && k1 =!= None,
                            svColIdx = FirstPosition[
                                Transpose[pkFast["PsiMat0"]],
                                col_ /; Norm[col - k1] <= 10^-10 * Max[1., Norm[k1]],
                                None, {1}
                            ];
                            If[ListQ[svColIdx],
                                svColIdx = First[svColIdx];
                                svCompatible = True;
                                Which[
                                    effMode === "StateFidelity" && k2 =!= None,
                                        evalSV = With[{cF = coeff, pT = k2, colI = svColIdx, dV = dim, rCols = Length[pkFast["Weights"]]},
                                            Function[{psiMatJ, tJ, j, wJ},
                                                Module[{psi = psiMatJ[[All, colI]], z, s = cF * wJ, chiOut},
                                                    z = Conjugate[pT] . psi;
                                                    chiOut = ConstantArray[0. + 0. I, {dV, rCols}];
                                                    chiOut[[All, colI]] = (2. * s * z) * pT;
                                                    {s * Re[z * Conjugate[z]], Developer`ToPackedArray[chiOut]}
                                                ]
                                            ]
                                        ],
                                    effMode === "Magnitude",
                                        evalSV = With[{cF = coeff, mat2 = m2, mat2H = m2H, scaleN = sn, colI = svColIdx, dV = dim, rCols = Length[pkFast["Weights"]]},
                                            Function[{psiMatJ, tJ, j, wJ},
                                                Module[{psi = psiMatJ[[All, colI]], z, s = cF * wJ, chiOut},
                                                    z = (Conjugate[psi] . (mat2H . psi)) / scaleN;
                                                    chiOut = ConstantArray[0. + 0. I, {dV, rCols}];
                                                    chiOut[[All, colI]] = (2. * s / scaleN) * (Conjugate[z] * (mat2H . psi) + z * (mat2 . psi));
                                                    {s * Re[z * Conjugate[z]], Developer`ToPackedArray[chiOut]}
                                                ]
                                            ]
                                        ],
                                    True,
                                        evalSV = With[{cF = coeff, mat2Sum = ToPackedComplexMatrix[m2 + m2H], scaleN = sn, colI = svColIdx, dV = dim, rCols = Length[pkFast["Weights"]]},
                                            Function[{psiMatJ, tJ, j, wJ},
                                                Module[{psi = psiMatJ[[All, colI]], val, s = cF * wJ, chiOut},
                                                    val = 0.5 * Re[Conjugate[psi] . (mat2Sum . psi)] / scaleN;
                                                    chiOut = ConstantArray[0. + 0. I, {dV, rCols}];
                                                    chiOut[[All, colI]] = (s / scaleN) * (mat2Sum . psi);
                                                    {s * val, Developer`ToPackedArray[chiOut]}
                                                ]
                                            ]
                                        ]
                                ],
                                svCompatible = False
                            ],
                            svCompatible = False
                        ];
                        If[!svCompatible, evalSV = None];
                        <|"Kind" -> "StatePenalty", "Space" -> "Hilbert", "TW" -> twComp,
                          "SVCompatible" -> svCompatible, "EvalHilbert" -> evalHilb, "EvalSV" -> evalSV, "EvalLiouv" -> None|>,

                        (* Liouville StatePenalty *)
                        vInit = toLiouvilleVector[initExpr, opbasis, basis, dim];
                        vTar = toLiouvilleVector[tarExpr, opbasis, basis, dim];
                        If[vInit === $Failed || vTar === $Failed,
                            Message[StateTarget::baddim, initExpr -> tarExpr, dim]; Return[$Failed, Table]
                        ];
                        rho0Cols = Transpose[compiledObj["InitialMatrix"]];
                        colMatches = FirstPosition[rho0Cols, col_ /; Max[Abs[col - vInit]] <= 10^-9 * Max[1., Norm[vInit]], None, {1}];
                        If[!ListQ[colMatches],
                            Message[StatePenalty::badspec, item]; Return[$Failed, Table]
                        ];
                        rCol = First[colMatches];
                        m2 = Quiet[Check[toHilbertOperatorMatrix[tarExpr, basis, Sqrt[dim]], $Failed]];
                        hermErr = If[MatrixQ[m2, NumericQ], Max[Abs[m2 - ConjugateTranspose[m2]]] / Max[1., Max[Abs[m2]]], 0.];
                        effMode = Which[
                            modeOpt === Automatic, If[hermErr > 10^-8, "Magnitude", "Expectation"],
                            MemberQ[{"Expectation", "Magnitude", "StateFidelity", "Projector"}, modeOpt], modeOpt,
                            True, Message[StatePenalty::badspec, item]; Return[$Failed, Table]
                        ];
                        sn = Which[
                            effMode === "Projector",
                                If[normOpt === Automatic || normOpt === "None" || normOpt === False, 1.,
                                    resolveScaleFactor[normOpt, vInit, vTar]],
                            True,
                                resolveScaleFactor[If[normOpt === Automatic, "Frobenius", normOpt], vInit, vTar]
                        ];
                        If[sn === $Failed || sn <= 0., Return[$Failed, Table]];
                        evalLiouv = If[effMode === "Magnitude",
                            With[{cF = coeff, oVec = vTar, scaleN = sn, rC = rCol, dlV = dim, rCount = Length[rho0Cols]},
                                Function[{rhoMatJ, tJ, j, wJ},
                                    Module[{z = (Conjugate[oVec] . rhoMatJ[[All, rC]]) / scaleN, s = cF * wJ, csOut},
                                        csOut = ConstantArray[0. + 0. I, {dlV, rCount}];
                                        csOut[[All, rC]] = ((2. * s * z) / scaleN) * oVec;
                                        {s * Re[z * Conjugate[z]], Developer`ToPackedArray[csOut]}
                                    ]
                                ]
                            ],
                            With[{cF = coeff, oVec = vTar, scaleN = sn, rC = rCol, dlV = dim, rCount = Length[rho0Cols]},
                                Function[{rhoMatJ, tJ, j, wJ},
                                    Module[{z = (Conjugate[oVec] . rhoMatJ[[All, rC]]) / scaleN, s = cF * wJ, csOut},
                                        csOut = ConstantArray[0. + 0. I, {dlV, rCount}];
                                        csOut[[All, rC]] = (s / scaleN) * oVec;
                                        {s * Re[z], Developer`ToPackedArray[csOut]}
                                    ]
                                ]
                            ]
                        ];
                        <|"Kind" -> "StatePenalty", "Space" -> "Liouville", "TW" -> twComp,
                          "SVCompatible" -> False, "EvalHilbert" -> None, "EvalSV" -> None, "EvalLiouv" -> evalLiouv|>
                    ]
                ],
                (* Custom running-cost callable: fn[stateOrU, t, j] -> {val, cot} *)
                If[!MatchQ[item, _Function | _Symbol | _CompiledFunction],
                    Message[StatePenalty::badspec, item]; Return[$Failed, Table]
                ];
                hasTimeDep = True;
                With[{cf = item, tw = <|"StaticQ" -> True, "Weights" -> ConstantArray[1., nt], "EvalRaw" -> None|>,
                      dV = dim, rCount = If[space === "Liouville", Length[compiledObj["Weights"]], 0]},
                    If[space === "Hilbert",
                        <|"Kind" -> "Custom", "Space" -> "Hilbert", "TW" -> tw, "SVCompatible" -> False,
                          "EvalHilbert" -> Function[{uJ, tJ, j, wJ},
                              Module[{res = cf[uJ, tJ, j]},
                                  If[!ListQ[res] || Length[res] < 2 || !NumericQ[res[[1]]],
                                      $Failed,
                                      {N[Re[res[[1]]]], Developer`ToPackedArray[N[Normal[res[[2]]]] + 0. I]}
                                  ]
                              ]
                          ],
                          "EvalSV" -> None, "EvalLiouv" -> None|>,
                        <|"Kind" -> "Custom", "Space" -> "Liouville", "TW" -> tw, "SVCompatible" -> False,
                          "EvalHilbert" -> None, "EvalSV" -> None,
                          "EvalLiouv" -> Function[{rhoMatJ, tJ, j, wJ},
                              Module[{res = cf[rhoMatJ, tJ, j], cotM},
                                  If[!ListQ[res] || Length[res] < 2 || !NumericQ[res[[1]]],
                                      $Failed,
                                      cotM = If[VectorQ[res[[2]], NumericQ] && rCount === 1,
                                          Transpose[{res[[2]]}],
                                          Normal[res[[2]]]
                                      ];
                                      {N[Re[res[[1]]]], Developer`ToPackedArray[N[cotM] + 0. I]}
                                  ]
                              ]
                          ]|>
                    ]
                ]
            ]
        ],
        {idx, Length[items]}
    ];

    If[!ListQ[compiledItems] || MemberQ[compiledItems, $Failed], Return[$Failed]];
    svOK = AllTrue[compiledItems, TrueQ[#["SVCompatible"]] &];
    twSpecs = #["TW"] & /@ compiledItems;
    termEvalFns = #["EvalHilbert"] & /@ compiledItems;
    termSvFns = #["EvalSV"] & /@ compiledItems;
    termLiouvFns = #["EvalLiouv"] & /@ compiledItems;

    <|
        "Space" -> space,
        "Dimension" -> dim,
        "NumTerms" -> Length[compiledItems],
        "CompatibleWithStateVectorQ" -> svOK,
        "HasTimeDependence" -> hasTimeDep,
        "TWSpecs" -> twSpecs,
        "HilbertTermEvaluators" -> termEvalFns,
        "StateVectorTermEvaluators" -> If[svOK, termSvFns, None],
        "LiouvilleTermEvaluators" -> termLiouvFns,
        "CustomEvaluators" -> Cases[items, Except[_StatePenalty]]
    |>
];

End[];
EndPackage[];
