(* ::Package:: *)

BeginPackage["Optimization`Solvers`",
{"Global`",
 "SpinDynamica`",
 "OptimalControl`Kernels`",
 "OptimalControl`Objectives`",
 "OptimalControl`GRAPE`",
 "OptimalControl`Optimizers`",
 "OptimalControl`SecondOrderOptimizers`"
}];

ClearAll[
    RunGeneralOptimizer,
    MaxEvaluations,
    TargetObjective,
    ParameterScaling,
    BoundTransform,
    InitialSimplexScale,
    OptimizerOptions
];

RunGeneralOptimizer::usage = "RunGeneralOptimizer[oracles, x0, bounds, opts] optimizes a real-valued objective using a physics-independent solver backend.\n" <>
"oracles is an Association containing any available subset of {\"EvalObj\", \"EvalGrad\", \"EvalWithCache\", \"EvalHVP\", \"EvalHessian\"}.\n" <>
"x0 is the initial real parameter vector, and bounds is {lbVec, ubVec} (with -Infinity / Infinity for unbounded coordinates).\n" <>
"Supported Optimizer values include \"LBFGS\" (default), \"BFGS\", \"CG\", \"SteepestAscent\", \"RegularizedLBFGS\", " <>
"\"HybridLBFGSNewtonCG\", \"HybridLBFGSNewtonRFO\", \"TrustRegionNewtonCG\", \"TrustRegionGLTR\", \"AdaptiveCubicRegularization\", " <>
"\"NewtonRFO\", \"NewtonShift\", \"NewtonModifiedCholesky\", \"NewtonSpectral\", \"NelderMead\", \"Powell\", and \"FindMinimum\".\n" <>
"ParameterScaling -> Automatic | None | {m1, m2, ...} specifies characteristic parameter magnitudes m_i so internal coordinates x_i / m_i are O(1).";

MaxEvaluations::usage = "MaxEvaluations is an option for general optimization routines specifying the maximum number of objective evaluations.";
TargetObjective::usage = "TargetObjective is an option for general optimization routines specifying a target objective threshold for early termination (or None to disable target-based stopping).";
ParameterScaling::usage = "ParameterScaling is an option specifying characteristic parameter magnitudes m_i so internal coordinates x_i / m_i are O(1), Automatic to infer scales, or None for unit scaling.";
BoundTransform::usage = "BoundTransform is an option specifying smooth bound reparameterization (\"None\", \"Sigmoid\", \"Log\", or Automatic) for unconstrained second-order solvers.";
InitialSimplexScale::usage = "InitialSimplexScale is an option for derivative-free Nelder-Mead optimization setting the relative initial simplex step size.";
OptimizerOptions::usage = "OptimizerOptions is an option specifying solver-specific sub-options as an Association or list of rules.";

RunGeneralOptimizer::badoracle = "The required evaluation oracle (`1`) is missing for Optimizer -> \"`2`\".";
RunGeneralOptimizer::badopt = "Invalid option or specification for RunGeneralOptimizer: `1`.";
RunGeneralOptimizer::badbounds = "Initial point or bounds dimension mismatch (`1`).";
RunGeneralOptimizer::badx0 = "Initial point cannot satisfy incompatible bounds.";

Options[RunGeneralOptimizer] = {
    Optimizer -> Automatic,
    LineSearch -> Automatic,
    WolfeParameters -> {1.*^-4, 0.9},
    LBFGSMemory -> 15,
    MaxIterations -> 200,
    MaxEvaluations -> Automatic,
    MaxLineSearchEvaluations -> 25,
    GradientTolerance -> 1.*^-8,
    FunctionTolerance -> 1.*^-12,
    StepTolerance -> 1.*^-12,
    TimeConstraint -> Infinity,
    TargetObjective -> None,
    TargetFidelity -> None,
    ParameterScaling -> Automatic,
    BoundTransform -> "None",
    OptimizerOptions -> <||>,
    InitialSimplexScale -> 0.05,
    RecordHistory -> True,
    StepMonitor -> None,
    MonitorInterval -> 1,
    Verbose -> False
};

Begin["`Private`"];

$FirstOrderGRAPEOptimizers = {"LBFGS", "BFGS", "CG", "SteepestAscent", "GradientAscent"};
$SecondOrderGRAPEOptimizers = {
    "HybridLBFGSNewtonCG", "HybridLBFGSNewtonRFO",
    "TrustRegionNewtonCG", "TrustRegionGLTR", "AdaptiveCubicRegularization",
    "NewtonRFO", "NewtonShift", "NewtonModifiedCholesky", "NewtonSpectral"
};
$DenseHessianOptimizers = {"HybridLBFGSNewtonRFO", "NewtonRFO", "NewtonShift", "NewtonModifiedCholesky", "NewtonSpectral"};
$DerivativeFreeOptimizers = {"NelderMead", "Powell"};

finiteNumQ[v_] := NumericQ[v] && TrueQ[Abs[N[v]] < Infinity];

(* Resolve characteristic parameter magnitudes -> internal multiplier scaleVec = 1 / magVec *)
resolveScaleVector[x0Vec_List, lbVec_List, ubVec_List, scaleSpec_] := Module[
    {n = Length[x0Vec], mags},
    mags = Which[
        scaleSpec === None || scaleSpec === False,
            ConstantArray[1., n],
        NumericQ[scaleSpec] && TrueQ[N[scaleSpec] > 0.],
            ConstantArray[N[scaleSpec], n],
        ListQ[scaleSpec] && Length[scaleSpec] === n && AllTrue[scaleSpec, NumericQ[#] && TrueQ[N[#] > 0.] &],
            N[scaleSpec],
        scaleSpec === Automatic,
            Table[
                With[{x0 = Abs[N[x0Vec[[i]]]], l = lbVec[[i]], u = ubVec[[i]]},
                    Which[
                        x0 > 1.*^-12, x0,
                        finiteNumQ[l] && finiteNumQ[u] && Abs[N[u - l]] > 1.*^-12, 0.5 * Abs[N[u - l]],
                        finiteNumQ[l] && Abs[N[l]] > 1.*^-12, Abs[N[l]],
                        finiteNumQ[u] && Abs[N[u]] > 1.*^-12, Abs[N[u]],
                        True, 1.0
                    ]
                ],
                {i, n}
            ],
        True,
            Message[RunGeneralOptimizer::badopt, "Invalid ParameterScaling specification."];
            Return[$Failed]
    ];
    Developer`ToPackedArray[1. / Max[#, 1.*^-12] & /@ mags]
];

(* Normalize bounds to {lbVec, ubVec, boundTypes, hasBoundsQ} *)
parseBounds[n_Integer, bounds_] := Module[{lb, ub, ok = True},
    If[bounds === None || bounds === Automatic,
        Return[{ConstantArray[-Infinity, n], ConstantArray[Infinity, n], False}]
    ];
    If[!ListQ[bounds] || Length[bounds] =!= 2,
        Message[RunGeneralOptimizer::badbounds, "Bounds must be {lbVec, ubVec} or None."];
        Return[$Failed]
    ];
    {lb, ub} = bounds;
    If[!ListQ[lb] || !ListQ[ub] || Length[lb] =!= n || Length[ub] =!= n,
        Message[RunGeneralOptimizer::badbounds, "Bounds vectors length must match x0."];
        Return[$Failed]
    ];
    lb = Replace[lb, {None -> -Infinity, x_?NumericQ :> N[Re[x]]}, {1}];
    ub = Replace[ub, {None -> Infinity, x_?NumericQ :> N[Re[x]]}, {1}];
    Do[
        If[!(lb[[i]] === -Infinity || finiteNumQ[lb[[i]]]) || !(ub[[i]] === Infinity || finiteNumQ[ub[[i]]]),
            ok = False; Break[]
        ];
        If[TrueQ[lb[[i]] >= ub[[i]]], ok = False; Break[]],
        {i, n}
    ];
    If[!ok,
        Message[RunGeneralOptimizer::badbounds, "Lower bound must be strictly less than upper bound and numeric or +/-Infinity."];
        Return[$Failed]
    ];
    {lb, ub, AnyTrue[Join[lb, ub], finiteNumQ]}
];

(* Smooth 1D diffeomorphisms phi_i(u_i) mapping R -> open (lb_i, ub_i) for unconstrained 2nd-order solvers *)
sigmoidSafe[u_?NumericQ] := With[{uc = Clip[N[u], {-40., 40.}]},
    If[uc >= 0., 1. / (1. + Exp[-uc]), With[{e = Exp[uc]}, e / (1. + e)]]
];

buildBoundTransformSpec[x0Vec_List, lbVec_List, ubVec_List, scaleVec_List, mode_String] := Module[
    {n = Length[x0Vec], kinds, cVec, u0Vec, eps},
    If[mode === "None", Return[None]];
    kinds = Table[
        Which[
            finiteNumQ[lbVec[[i]]] && finiteNumQ[ubVec[[i]]], "TwoSided",
            finiteNumQ[lbVec[[i]]], "LowerOnly",
            finiteNumQ[ubVec[[i]]], "UpperOnly",
            True, "Unbounded"
        ],
        {i, n}
    ];
    cVec = ConstantArray[1., n];
    u0Vec = ConstantArray[0., n];
    Do[
        Switch[kinds[[i]],
            "Unbounded",
                cVec[[i]] = 1. / scaleVec[[i]];
                u0Vec[[i]] = N[x0Vec[[i]]] * scaleVec[[i]],
            "TwoSided",
                With[{l = lbVec[[i]], u = ubVec[[i]], span = ubVec[[i]] - lbVec[[i]]},
                    eps = 1.*^-7 * span;
                    With[{xc = Clip[N[x0Vec[[i]]], {l + eps, u - eps}]},
                        cVec[[i]] = span;
                        u0Vec[[i]] = Log[(xc - l) / (u - xc)]
                    ]
                ],
            "LowerOnly",
                With[{l = lbVec[[i]], c0 = Max[N[x0Vec[[i]] - lbVec[[i]]], 1. / scaleVec[[i]], 1.*^-4]},
                    cVec[[i]] = c0;
                    eps = 1.*^-7 * c0;
                    With[{xc = Max[N[x0Vec[[i]]], l + eps]},
                        u0Vec[[i]] = If[mode === "Log", Log[(xc - l) / c0], Log[(xc - l) / c0]]
                    ]
                ],
            "UpperOnly",
                With[{u = ubVec[[i]], c0 = Max[N[ubVec[[i]] - x0Vec[[i]]], 1. / scaleVec[[i]], 1.*^-4]},
                    cVec[[i]] = c0;
                    eps = 1.*^-7 * c0;
                    With[{xc = Min[N[x0Vec[[i]]], u - eps]},
                        u0Vec[[i]] = -Log[(u - xc) / c0]
                    ]
                ]
        ],
        {i, n}
    ];
    <|"Kinds" -> kinds, "C" -> Developer`ToPackedArray[cVec], "Lb" -> lbVec, "Ub" -> ubVec,
      "U0" -> Developer`ToPackedArray[u0Vec], "Dim" -> n|>
];

evalBoundTransform[spec_Association, uVec_List] := Module[
    {n = spec["Dim"], kinds = spec["Kinds"], c = spec["C"], lb = spec["Lb"], ub = spec["Ub"],
     x = ConstantArray[0., spec["Dim"]], dp = ConstantArray[1., spec["Dim"]], dpp = ConstantArray[0., spec["Dim"]],
     s, uc, ex},
    Do[
        Switch[kinds[[i]],
            "Unbounded",
                x[[i]] = c[[i]] * uVec[[i]];
                dp[[i]] = c[[i]];
                dpp[[i]] = 0.,
            "TwoSided",
                s = sigmoidSafe[uVec[[i]]];
                x[[i]] = lb[[i]] + c[[i]] * s;
                dp[[i]] = c[[i]] * s * (1. - s);
                dpp[[i]] = dp[[i]] * (1. - 2. * s),
            "LowerOnly",
                uc = Clip[N[uVec[[i]]], {-40., 40.}];
                ex = c[[i]] * Exp[uc];
                x[[i]] = lb[[i]] + ex;
                dp[[i]] = ex;
                dpp[[i]] = ex,
            "UpperOnly",
                uc = Clip[N[uVec[[i]]], {-40., 40.}];
                ex = c[[i]] * Exp[-uc];
                x[[i]] = ub[[i]] - ex;
                dp[[i]] = ex;
                dpp[[i]] = -ex
        ],
        {i, n}
    ];
    {Developer`ToPackedArray[x], Developer`ToPackedArray[dp], Developer`ToPackedArray[dpp]}
];

(* Normalize status keys to general-optimization vocabulary while retaining GRAPE status metadata *)
normalizeSolverStatus[st_String, hasTargetQ_] := Which[
    st === "TargetFidelityReached" && !TrueQ[hasTargetQ], "Converged",
    st === "TargetFidelityReached", "TargetReached",
    True, st
];

(* Built-in Box-Aware Derivative-Free Solver 1: Adaptive Nelder-Mead (Gao-Han 2012) *)
runNelderMeadSolver[evalObjFn_, x0Vec_List, scaleVec_List, lbVec_List, ubVec_List, hasBoundsQ_, optsAssoc_Association] := Module[
    {n = Length[x0Vec], invS, lbT, ubT, clipT, maxIter, maxEvals, funTol, stepTol, tLimit, targetObj,
     recHist, monOpt, monInt, simScale, alphaNM, betaNM, gammaNM, deltaNM,
     tStart, dynTime = 0., objEvals = 0, nIter = 0, nRestarts = 0, maxRestarts = 4, status = "MaxIterationsReached",
     verts, fVals, order, thBest, fBest, fBestNew, xBest, th0, evalT, objHist = {}, gHist = {},
     centroid, thR, fR, thE, fE, thC, fC, shrinkQ, fRange, stepSpread, step0, stepR, simpRank, i},
    invS = Developer`ToPackedArray[1. / scaleVec];
    lbT = If[hasBoundsQ, lbVec * scaleVec, ConstantArray[-Infinity, n]];
    ubT = If[hasBoundsQ, ubVec * scaleVec, ConstantArray[Infinity, n]];
    clipT[th_] := If[hasBoundsQ, Developer`ToPackedArray[MapThread[Clip[#1, {#2, #3}] &, {th, lbT, ubT}]], Developer`ToPackedArray[N[th]]];

    maxIter = Max[1, Round[Lookup[optsAssoc, MaxIterations, 200]]];
    maxEvals = With[{me = Lookup[optsAssoc, MaxEvaluations, Automatic]}, If[IntegerQ[me], Max[n + 2, me], Max[250, 40 * n * maxIter]]];
    funTol = Max[0., N[Lookup[optsAssoc, FunctionTolerance, 1.*^-12]]];
    stepTol = Max[0., N[Lookup[optsAssoc, StepTolerance, 1.*^-12]]];
    tLimit = N[Lookup[optsAssoc, TimeConstraint, Infinity]];
    targetObj = Lookup[optsAssoc, TargetObjective, None];
    If[!NumericQ[targetObj], targetObj = Lookup[optsAssoc, TargetFidelity, None]];
    recHist = TrueQ[Lookup[optsAssoc, RecordHistory, True]];
    simScale = Max[1.*^-4, N[Lookup[optsAssoc, InitialSimplexScale, 0.05]]];
    monOpt = Lookup[optsAssoc, StepMonitor, None];
    monInt = Max[1, Round[Lookup[optsAssoc, MonitorInterval, 1]]];

    alphaNM = 1.0;
    betaNM = 1.0 + 2.0 / Max[n, 1];
    gammaNM = 0.75 - 1.0 / (2.0 * Max[n, 1]);
    deltaNM = 1.0 - 1.0 / Max[n, 1];

    evalT[th_] := Module[{t0 = AbsoluteTime[], xc, r, val},
        xc = Developer`ToPackedArray[clipT[th] * invS];
        r = evalObjFn[xc];
        dynTime += AbsoluteTime[] - t0;
        objEvals++;
        val = Which[
            NumericQ[r], N[Re[r]],
            AssociationQ[r] && NumericQ[r["Objective"]], N[Re[r["Objective"]]],
            True, $Failed
        ];
        val
    ];

    tStart = AbsoluteTime[];
    th0 = clipT[N[x0Vec] * scaleVec];
    verts = Table[th0, {n + 1}];
    Do[
        step0 = simScale * Max[Abs[th0[[i]]], 1.0];
        If[finiteNumQ[ubT[[i]]] && th0[[i]] + step0 > ubT[[i]],
            If[finiteNumQ[lbT[[i]]] && th0[[i]] - step0 < lbT[[i]],
                step0 = If[(ubT[[i]] - th0[[i]]) >= (th0[[i]] - lbT[[i]]), 0.5 * (ubT[[i]] - th0[[i]]), -0.5 * (th0[[i]] - lbT[[i]])],
                step0 = -step0
            ]
        ];
        verts[[i + 1, i]] = th0[[i]] + step0;
        verts[[i + 1]] = clipT[verts[[i + 1]]],
        {i, n}
    ];
    fVals = evalT /@ verts;
    If[MemberQ[fVals, $Failed], Return[$Failed]];

    order = Ordering[fVals, All, Greater];
    verts = verts[[order]];
    fVals = fVals[[order]];
    thBest = verts[[1]];
    fBest = fVals[[1]];
    If[recHist, AppendTo[objHist, fBest]];

    While[True,
        If[NumericQ[targetObj] && TrueQ[fBest >= N[Re[targetObj]]], status = "TargetReached"; Break[]];
        If[nIter >= maxIter, status = "MaxIterationsReached"; Break[]];
        If[objEvals >= maxEvals, status = "MaxEvaluationsReached"; Break[]];
        If[AbsoluteTime[] - tStart >= tLimit, status = "TimeConstraintReached"; Break[]];
        fRange = Abs[fVals[[1]] - fVals[[-1]]];
        stepSpread = Max[Norm[# - verts[[1]], Infinity] & /@ Rest[verts]];
        simpRank = If[hasBoundsQ && n > 1 && nIter > 0 && nRestarts < maxRestarts,
            MatrixRank[(# - verts[[1]]) & /@ Rest[verts], Tolerance -> 1.*^-10],
            n
        ];
        If[(nIter > 0 && fRange <= funTol * Max[1., Abs[fBest]] && stepSpread <= stepTol) ||
           (hasBoundsQ && nIter > 0 && simpRank < n && nRestarts < maxRestarts),
            If[hasBoundsQ && nRestarts < maxRestarts,
                nRestarts++;
                verts = Table[thBest, {n + 1}];
                fVals[[1]] = fBest;
                Do[
                    stepR = 0.05 * Max[Abs[thBest[[i]]], 1.0] * (0.5^(nRestarts - 1));
                    If[finiteNumQ[ubT[[i]]] && thBest[[i]] + stepR > ubT[[i]],
                        If[finiteNumQ[lbT[[i]]] && thBest[[i]] - stepR < lbT[[i]],
                            stepR = If[(ubT[[i]] - thBest[[i]]) >= (thBest[[i]] - lbT[[i]]), 0.5 * (ubT[[i]] - thBest[[i]]), -0.5 * (thBest[[i]] - lbT[[i]])],
                            stepR = -stepR
                        ]
                    ];
                    verts[[i + 1, i]] = thBest[[i]] + stepR;
                    verts[[i + 1]] = clipT[verts[[i + 1]]];
                    fVals[[i + 1]] = evalT[verts[[i + 1]]];
                    If[fVals[[i + 1]] === $Failed, status = "EvaluationFailed"; Break[]],
                    {i, n}
                ];
                If[status === "EvaluationFailed", Break[]];
                order = Ordering[fVals, All, Greater];
                verts = verts[[order]];
                fVals = fVals[[order]];
                fBestNew = fVals[[1]];
                If[fBestNew <= fBest + funTol * Max[1., Abs[fBest]] &&
                   fRange <= funTol * Max[1., Abs[fBest]] && stepSpread <= stepTol,
                    thBest = verts[[1]]; fBest = fBestNew;
                    status = "Converged"; Break[]
                ];
                thBest = verts[[1]];
                fBest = fBestNew;
                Continue[],
                status = "Converged"; Break[]
            ]
        ];

        nIter++;
        centroid = Developer`ToPackedArray[Mean[Most[verts]]];
        thR = clipT[centroid + alphaNM * (centroid - verts[[-1]])];
        fR = evalT[thR];
        If[fR === $Failed, status = "EvaluationFailed"; Break[]];

        Which[
            fR > fVals[[1]],
                thE = clipT[centroid + betaNM * (thR - centroid)];
                fE = evalT[thE];
                If[fE === $Failed, status = "EvaluationFailed"; Break[]];
                If[fE > fR,
                    verts[[-1]] = thE; fVals[[-1]] = fE,
                    verts[[-1]] = thR; fVals[[-1]] = fR
                ],
            fR >= fVals[[-2]],
                verts[[-1]] = thR; fVals[[-1]] = fR,
            True,
                shrinkQ = False;
                If[fR > fVals[[-1]],
                    thC = clipT[centroid + gammaNM * (thR - centroid)];
                    fC = evalT[thC];
                    If[fC === $Failed, status = "EvaluationFailed"; Break[]];
                    If[fC >= fR, verts[[-1]] = thC; fVals[[-1]] = fC, shrinkQ = True],
                    thC = clipT[centroid - gammaNM * (centroid - verts[[-1]])];
                    fC = evalT[thC];
                    If[fC === $Failed, status = "EvaluationFailed"; Break[]];
                    If[fC > fVals[[-1]], verts[[-1]] = thC; fVals[[-1]] = fC, shrinkQ = True]
                ];
                If[shrinkQ,
                    Do[
                        verts[[j]] = clipT[verts[[1]] + deltaNM * (verts[[j]] - verts[[1]])];
                        fVals[[j]] = evalT[verts[[j]]];
                        If[fVals[[j]] === $Failed, status = "EvaluationFailed"; Break[]],
                        {j, 2, n + 1}
                    ];
                    If[status === "EvaluationFailed", Break[]]
                ]
        ];

        order = Ordering[fVals, All, Greater];
        verts = verts[[order]];
        fVals = fVals[[order]];
        thBest = verts[[1]];
        fBest = fVals[[1]];
        If[recHist, AppendTo[objHist, fBest]];
        If[monOpt =!= None && Mod[nIter, monInt] === 0,
            Quiet[Check[monOpt[<|"Iteration" -> nIter, "Objective" -> fBest, "ParametersVec" -> thBest * invS|>], Null]]
        ]
    ];

    xBest = Developer`ToPackedArray[thBest * invS];
    With[{elapsed = AbsoluteTime[] - tStart},
        <|
            "ParametersVec" -> xBest,
            "Objective" -> fBest,
            "GradientNorm" -> None,
            "Iterations" -> nIter,
            "Evaluations" -> objEvals,
            "ObjectiveEvaluations" -> objEvals,
            "GradientEvaluations" -> 0,
            "HVPEvaluations" -> 0,
            "HessianEvaluations" -> 0,
            "Status" -> status,
            "ElapsedTime" -> elapsed,
            "DynamicsTime" -> dynTime,
            "OptimizerTime" -> Max[0., elapsed - dynTime],
            "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
            "GradientNormHistory" -> {},
            "Optimizer" -> "NelderMead",
            "LineSearchMethod" -> "None",
            "ScaleVector" -> scaleVec,
            "Extra" -> <|"SimplexSpread" -> Max[Norm[# - verts[[1]], Infinity] & /@ Rest[verts]]|>
        |>
    ]
];

(* Built-in Box-Aware Derivative-Free Solver 2: Powell Direction-Set + Brent 1D Search *)
runPowellSolver[evalObjFn_, x0Vec_List, scaleVec_List, lbVec_List, ubVec_List, hasBoundsQ_, optsAssoc_Association] := Module[
    {n = Length[x0Vec], invS, lbT, ubT, clipT, maxIter, maxEvals, funTol, stepTol, tLimit, targetObj,
     recHist, monOpt, monInt, tStart, dynTime = 0., objEvals = 0, nIter = 0, status = "MaxIterationsReached",
     th, fVal, dirs, thStart, fStart, thExtrap, fExtrap, fMaxDrop, idxDrop, evalT, lineSearch1D,
     df, objHist = {}, thNew, fNew, stopQ = False},
    invS = Developer`ToPackedArray[1. / scaleVec];
    lbT = If[hasBoundsQ, lbVec * scaleVec, ConstantArray[-Infinity, n]];
    ubT = If[hasBoundsQ, ubVec * scaleVec, ConstantArray[Infinity, n]];
    clipT[v_] := If[hasBoundsQ, Developer`ToPackedArray[MapThread[Clip[#1, {#2, #3}] &, {v, lbT, ubT}]], Developer`ToPackedArray[N[v]]];

    maxIter = Max[1, Round[Lookup[optsAssoc, MaxIterations, 200]]];
    maxEvals = With[{me = Lookup[optsAssoc, MaxEvaluations, Automatic]}, If[IntegerQ[me], Max[n + 2, me], Max[250, 40 * n * maxIter]]];
    funTol = Max[0., N[Lookup[optsAssoc, FunctionTolerance, 1.*^-12]]];
    stepTol = Max[0., N[Lookup[optsAssoc, StepTolerance, 1.*^-12]]];
    tLimit = N[Lookup[optsAssoc, TimeConstraint, Infinity]];
    targetObj = Lookup[optsAssoc, TargetObjective, None];
    If[!NumericQ[targetObj], targetObj = Lookup[optsAssoc, TargetFidelity, None]];
    recHist = TrueQ[Lookup[optsAssoc, RecordHistory, True]];
    monOpt = Lookup[optsAssoc, StepMonitor, None];
    monInt = Max[1, Round[Lookup[optsAssoc, MonitorInterval, 1]]];

    evalT[v_] := Module[{t0 = AbsoluteTime[], xc, r},
        xc = Developer`ToPackedArray[clipT[v] * invS];
        r = evalObjFn[xc];
        dynTime += AbsoluteTime[] - t0;
        objEvals++;
        Which[
            NumericQ[r], N[Re[r]],
            AssociationQ[r] && NumericQ[r["Objective"]], N[Re[r["Objective"]]],
            True, $Failed
        ]
    ];

    (* Exact box-constrained bracket + golden/parabolic 1D maximization along dir from thIn *)
    lineSearch1D[thIn_List, fIn_?NumericQ, dir_List] := Module[
        {aLo = -Infinity, aHi = Infinity, dNorm = Norm[dir], step, a, b, c, fa = -Infinity, fb = -Infinity, fc = -Infinity,
         phiG = 0.6180339887498949, x1, x2, fx1, fx2, evalLine, k, tolA, lim = 20.0},
        If[dNorm <= 1.*^-15, Return[{thIn, fIn}]];
        If[hasBoundsQ,
            Do[
                With[{di = dir[[i]], xi = thIn[[i]], li = lbT[[i]], ui = ubT[[i]]},
                    Which[
                        di > 1.*^-15,
                            If[finiteNumQ[li], aLo = Max[aLo, (li - xi) / di]];
                            If[finiteNumQ[ui], aHi = Min[aHi, (ui - xi) / di]],
                        di < -1.*^-15,
                            If[finiteNumQ[li], aHi = Min[aHi, (li - xi) / di]];
                            If[finiteNumQ[ui], aLo = Max[aLo, (ui - xi) / di]]
                    ]
                ],
                {i, n}
            ]
        ];
        aLo = Max[aLo, -lim / dNorm];
        aHi = Min[aHi, lim / dNorm];
        If[aLo >= aHi - 1.*^-15, Return[{thIn, fIn}]];

        evalLine[alpha_] := evalT[thIn + alpha * dir];
        step = Min[0.1 * Max[Norm[thIn], 1.0] / dNorm, 0.25 * (aHi - aLo)];
        If[aHi <= 1.*^-12, step = -step];
        b = 0.; fb = fIn;

        (* Bracket maximum within [aLo, aHi] *)
        If[step > 0 && b + step <= aHi,
            c = b + step; fc = evalLine[c];
            If[fc === $Failed, Return[$Failed]];
            If[fc < fb && b - step >= aLo,
                a = c; fa = fc;
                step = -step;
                c = b + step; fc = evalLine[c];
                If[fc === $Failed, Return[$Failed]],
                a = Max[aLo, b - step]
            ],
            step = -Min[Abs[step], 0.25 * (aHi - aLo)];
            c = Max[aLo, b + step]; fc = evalLine[c];
            If[fc === $Failed, Return[$Failed]];
            a = b; fa = fb;
        ];

        (* Expand bracket toward uphill direction while inside [aLo, aHi] *)
        Do[
            If[!(fc >= fb), Break[]];
            a = b; fa = fb;
            b = c; fb = fc;
            step = 1.618 * step;
            c = Clip[b + step, {aLo, aHi}];
            If[Abs[c - b] < 1.*^-14, Break[]];
            fc = evalLine[c];
            If[fc === $Failed, Return[$Failed]];
            If[(step > 0 && c >= aHi) || (step < 0 && c <= aLo), Break[]],
            {14}
        ];

        If[a > c, {a, c} = {c, a}; {fa, fc} = {fc, fa}];
        a = Max[a, aLo]; c = Min[c, aHi];
        If[c - a <= 1.*^-14,
            With[{bestA = First[Ordering[{fa, fb, fc}, -1]] /. {1 -> a, 2 -> b, 3 -> c}},
                Return[{clipT[thIn + bestA * dir], Max[fa, fb, fc]}]
            ]
        ];

        (* Golden-section refinement on [a, c] (robust at boundaries and kinks) *)
        x1 = c - phiG * (c - a);
        x2 = a + phiG * (c - a);
        fx1 = evalLine[x1]; fx2 = evalLine[x2];
        If[fx1 === $Failed || fx2 === $Failed, Return[$Failed]];
        tolA = Max[1.*^-12, 0.25 * stepTol / dNorm];
        For[k = 0, k < 35 && (c - a) > tolA * Max[1., Abs[x1] + Abs[x2]] && objEvals < maxEvals, k++,
            If[fx1 < fx2,
                a = x1; fa = fx1; x1 = x2; fx1 = fx2;
                x2 = a + phiG * (c - a); fx2 = evalLine[x2];
                If[fx2 === $Failed, Return[$Failed]],
                c = x2; fc = fx2; x2 = x1; fx2 = fx1;
                x1 = c - phiG * (c - a); fx1 = evalLine[x1];
                If[fx1 === $Failed, Return[$Failed]]
            ]
        ];
        With[{cands = {{0., fIn}, {a, If[NumericQ[fa], fa, -Infinity]}, {x1, fx1}, {x2, fx2}, {c, If[NumericQ[fc], fc, -Infinity]}}},
            With[{best = First[MaximalBy[cands, Last]]},
                {clipT[thIn + best[[1]] * dir], best[[2]]}
            ]
        ]
    ];

    tStart = AbsoluteTime[];
    th = clipT[N[x0Vec] * scaleVec];
    fVal = evalT[th];
    If[fVal === $Failed, Return[$Failed]];
    If[recHist, AppendTo[objHist, fVal]];
    dirs = N[IdentityMatrix[n]];

    While[True,
        If[NumericQ[targetObj] && TrueQ[fVal >= N[Re[targetObj]]], status = "TargetReached"; Break[]];
        If[nIter >= maxIter, status = "MaxIterationsReached"; Break[]];
        If[objEvals >= maxEvals, status = "MaxEvaluationsReached"; Break[]];
        If[AbsoluteTime[] - tStart >= tLimit, status = "TimeConstraintReached"; Break[]];

        nIter++;
        thStart = th;
        fStart = fVal;
        fMaxDrop = 0.;
        idxDrop = 1;

        Do[
            With[{fPrev = fVal, lsRes = lineSearch1D[th, fVal, dirs[[i]]]},
                If[lsRes === $Failed, status = "EvaluationFailed"; stopQ = True; Break[]];
                {th, fVal} = lsRes;
                df = fVal - fPrev;
                If[df > fMaxDrop, fMaxDrop = df; idxDrop = i];
                If[objEvals >= maxEvals || AbsoluteTime[] - tStart >= tLimit, Break[]]
            ],
            {i, n}
        ];
        If[stopQ, Break[]];

        (* Powell criterion for updating direction set *)
        thExtrap = clipT[2. * th - thStart];
        fExtrap = evalT[thExtrap];
        If[fExtrap === $Failed, status = "EvaluationFailed"; Break[]];

        If[fExtrap > fStart &&
            2. * (fStart - 2. * fVal + fExtrap) * (fVal - fStart - fMaxDrop)^2 < fMaxDrop * (fExtrap - fStart)^2,
            With[{newDir = th - thStart},
                If[Norm[newDir] > 1.*^-14,
                    With[{ls2 = lineSearch1D[th, fVal, newDir]},
                        If[ls2 === $Failed, status = "EvaluationFailed"; Break[]];
                        {thNew, fNew} = ls2;
                        If[fNew > fVal, {th, fVal} = {thNew, fNew}];
                        dirs[[idxDrop]] = dirs[[-1]];
                        dirs[[-1]] = newDir / Norm[newDir]
                    ]
                ]
            ]
        ];
        If[Mod[nIter, Max[n, 2]] === 0, dirs = N[IdentityMatrix[n]]];

        If[recHist, AppendTo[objHist, fVal]];
        If[monOpt =!= None && Mod[nIter, monInt] === 0,
            Quiet[Check[monOpt[<|"Iteration" -> nIter, "Objective" -> fVal, "ParametersVec" -> th * invS|>], Null]]
        ];
        If[Abs[fVal - fStart] <= funTol * Max[1., Abs[fVal]] && Norm[th - thStart] <= stepTol * Max[1., Norm[th]],
            status = "Converged"; Break[]
        ]
    ];

    With[{elapsed = AbsoluteTime[] - tStart},
        <|
            "ParametersVec" -> Developer`ToPackedArray[th * invS],
            "Objective" -> fVal,
            "GradientNorm" -> None,
            "Iterations" -> nIter,
            "Evaluations" -> objEvals,
            "ObjectiveEvaluations" -> objEvals,
            "GradientEvaluations" -> 0,
            "HVPEvaluations" -> 0,
            "HessianEvaluations" -> 0,
            "Status" -> status,
            "ElapsedTime" -> elapsed,
            "DynamicsTime" -> dynTime,
            "OptimizerTime" -> Max[0., elapsed - dynTime],
            "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
            "GradientNormHistory" -> {},
            "Optimizer" -> "Powell",
            "LineSearchMethod" -> "BrentGolden",
            "ScaleVector" -> scaleVec,
            "Extra" -> <||>
        |>
    ]
];

(* Compiled FindMinimum Bridge (supports compiled Obj + optional analytical/FD Grad & box bounds) *)
runFindMinimumBridge[evalObjFn_, evalGradFn_, x0Vec_List, scaleVec_List, lbVec_List, ubVec_List, hasBoundsQ_, optsAssoc_Association] := Module[
    {n = Length[x0Vec], invS, th0, lbT, ubT, subOpts, fmMeth, maxIter, accGoal, precGoal, tLimit,
     recHist, vars, constraints, fmOpts, tStart, dynTime = 0., objEvals = 0, gradEvals = 0,
     objHist = {}, gHist = {}, nIter = 0, minRes, minVal, solRules, xSol, fSol, fLast = None, gLast = None,
     status = "Converged", monOpt, monInt, objWrap, gradWrap},
    invS = Developer`ToPackedArray[1. / scaleVec];
    th0 = Developer`ToPackedArray[N[x0Vec] * scaleVec];
    lbT = If[hasBoundsQ, lbVec * scaleVec, ConstantArray[-Infinity, n]];
    ubT = If[hasBoundsQ, ubVec * scaleVec, ConstantArray[Infinity, n]];
    subOpts = Lookup[optsAssoc, OptimizerOptions, <||>];
    If[ListQ[subOpts], subOpts = Association[subOpts]];
    If[!AssociationQ[subOpts], subOpts = <||>];
    fmMeth = Lookup[subOpts, "Method", Lookup[subOpts, Method, Automatic]];
    maxIter = Max[1, Round[Lookup[optsAssoc, MaxIterations, 200]]];
    accGoal = -Log10[Max[1.*^-15, Min[1.*^-2, N[Lookup[optsAssoc, FunctionTolerance, 1.*^-12]]]]];
    precGoal = -Log10[Max[1.*^-15, Min[1.*^-2, N[Lookup[optsAssoc, StepTolerance, 1.*^-12]]]]];
    tLimit = N[Lookup[optsAssoc, TimeConstraint, Infinity]];
    recHist = TrueQ[Lookup[optsAssoc, RecordHistory, True]];
    monOpt = Lookup[optsAssoc, StepMonitor, None];
    monInt = Max[1, Round[Lookup[optsAssoc, MonitorInterval, 1]]];

    vars = Table[Unique["goptFMVar"], {n}];

    objWrap[args___?NumericQ] := Module[{t0 = AbsoluteTime[], xc, r, val},
        xc = Developer`ToPackedArray[{args} * invS];
        r = evalObjFn[xc];
        dynTime += AbsoluteTime[] - t0;
        objEvals++;
        val = Which[
            NumericQ[r], -N[Re[r]],
            AssociationQ[r] && NumericQ[r["Objective"]], -N[Re[r["Objective"]]],
            True, $Failed
        ];
        If[NumericQ[val], fLast = -val];
        val
    ];

    gradWrap[args___?NumericQ] := Module[{t0 = AbsoluteTime[], xc, r, gX},
        xc = Developer`ToPackedArray[{args} * invS];
        r = evalGradFn[xc];
        dynTime += AbsoluteTime[] - t0;
        gradEvals++;
        gX = Which[
            AssociationQ[r] && KeyExistsQ[r, "Gradient"], Flatten[r["Gradient"]],
            VectorQ[r, NumericQ], r,
            True, Return[$Failed]
        ];
        gLast = Developer`ToPackedArray[-(N[gX] * invS)];
        gLast
    ];

    constraints = If[hasBoundsQ,
        DeleteCases[
            Table[
                Which[
                    finiteNumQ[lbT[[i]]] && finiteNumQ[ubT[[i]]], lbT[[i]] <= vars[[i]] <= ubT[[i]],
                    finiteNumQ[lbT[[i]]], lbT[[i]] <= vars[[i]],
                    finiteNumQ[ubT[[i]]], vars[[i]] <= ubT[[i]],
                    True, True
                ],
                {i, n}
            ],
            True
        ],
        {}
    ];

    fmOpts = Join[
        {
            MaxIterations -> maxIter,
            AccuracyGoal -> Min[12, Max[4, accGoal]],
            PrecisionGoal -> Min[12, Max[4, precGoal]],
            Method -> fmMeth
        },
        If[evalGradFn =!= None && fmMeth =!= "PrincipalAxis", {Gradient -> (gradWrap @@ vars), Evaluated -> False}, {}],
        {
            StepMonitor :> (
                nIter++;
                If[recHist,
                    If[NumericQ[fLast], AppendTo[objHist, fLast]];
                    If[VectorQ[gLast, NumericQ], AppendTo[gHist, Norm[gLast]]]
                ];
                If[monOpt =!= None && Mod[nIter, monInt] === 0,
                    Quiet[Check[monOpt[<|
                        "Iteration" -> nIter,
                        "Objective" -> fLast,
                        "GradientNorm" -> If[VectorQ[gLast, NumericQ], Norm[gLast], None],
                        "ParametersVec" -> Developer`ToPackedArray[vars * invS]
                    |>], Null]]
                ]
            )
        }
    ];

    tStart = AbsoluteTime[];
    minRes = TimeConstrained[
        Quiet[
            If[constraints === {},
                FindMinimum[objWrap @@ vars, Transpose[{vars, th0}], Evaluate[Sequence @@ fmOpts]],
                FindMinimum[Prepend[constraints, objWrap @@ vars] // Flatten, Transpose[{vars, th0}], Evaluate[Sequence @@ fmOpts]]
            ]
        ],
        tLimit,
        $Aborted
    ];
    If[minRes === $Aborted || minRes === $Failed || !MatchQ[minRes, {_?NumericQ, {_Rule ..}}],
        Return[$Failed]
    ];
    {minVal, solRules} = minRes;
    xSol = Developer`ToPackedArray[(vars /. solRules) * invS];
    fSol = -N[Re[minVal]];
    If[recHist && objHist === {}, objHist = {fSol}];

    With[{elapsed = AbsoluteTime[] - tStart},
        <|
            "ParametersVec" -> xSol,
            "Objective" -> fSol,
            "GradientNorm" -> If[VectorQ[gLast, NumericQ], Norm[gLast], None],
            "Iterations" -> nIter,
            "Evaluations" -> objEvals + gradEvals,
            "ObjectiveEvaluations" -> objEvals,
            "GradientEvaluations" -> gradEvals,
            "HVPEvaluations" -> 0,
            "HessianEvaluations" -> 0,
            "Status" -> status,
            "ElapsedTime" -> elapsed,
            "DynamicsTime" -> dynTime,
            "OptimizerTime" -> Max[0., elapsed - dynTime],
            "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
            "GradientNormHistory" -> Developer`ToPackedArray[gHist],
            "Optimizer" -> "FindMinimum",
            "LineSearchMethod" -> ToString[fmMeth],
            "ScaleVector" -> scaleVec,
            "Extra" -> <||>
        |>
    ]
];

(* Unified Physics-Independent Dispatcher *)
RunGeneralOptimizer[oracles_Association, x0In_List, boundsIn_: None, opts:(_Rule|_RuleDelayed)...] := Module[
    {optsAssoc, optsClean, monOpt, nDim, x0Vec, bParsed, lbVec, ubVec, hasBoundsQ, scaleVec, optName, lsOpt,
     targetVal, btMode, btSpec, subOpts, inlineSubOpts = <||>, evalObjFn, evalGradFn, evalWithCacheFn, evalHVPFn, evalHessFn,
     evalFlat, evalObjWrapped, hasTargetQ, optsForGRAPE, res, gHist, gNormFinal},
    monOpt = First[
        Join[
            Cases[{opts}, HoldPattern[StepMonitor -> v_] :> v, {1}],
            Cases[{opts}, HoldPattern[StepMonitor :> rhs_] :>
                If[MatchQ[HoldComplete[rhs], HoldComplete[_Function]], rhs, Function[{goptMonPayload}, rhs]], {1}]
        ],
        None
    ];
    optsClean = FilterRules[{opts}, Except[StepMonitor]];
    optsAssoc = Join[Association[Options[RunGeneralOptimizer]], Association[optsClean], <|StepMonitor -> monOpt|>];
    If[!VectorQ[x0In, NumericQ] || Length[x0In] === 0,
        Message[RunGeneralOptimizer::badopt, "Initial point x0 must be a non-empty numeric vector."];
        Return[$Failed]
    ];
    x0Vec = Developer`ToPackedArray[N[Re[x0In]]];
    nDim = Length[x0Vec];
    bParsed = parseBounds[nDim, boundsIn];
    If[bParsed === $Failed, Return[$Failed]];
    {lbVec, ubVec, hasBoundsQ} = bParsed;

    (* Check initial feasibility if bounds provided *)
    If[hasBoundsQ,
        Do[
            If[(finiteNumQ[lbVec[[i]]] && TrueQ[x0Vec[[i]] < lbVec[[i]] - 1.*^-12]) ||
               (finiteNumQ[ubVec[[i]]] && TrueQ[x0Vec[[i]] > ubVec[[i]] + 1.*^-12]),
                Message[RunGeneralOptimizer::badx0];
                Return[$Failed]
            ],
            {i, nDim}
        ];
        x0Vec = Developer`ToPackedArray[MapThread[Clip[#1, {#2, #3}] &, {x0Vec, lbVec, ubVec}]]
    ];

    scaleVec = resolveScaleVector[x0Vec, lbVec, ubVec, optsAssoc[ParameterScaling]];
    If[scaleVec === $Failed, Return[$Failed]];

    evalObjFn = Lookup[oracles, "EvalObj", None];
    evalGradFn = Lookup[oracles, "EvalGrad", None];
    evalWithCacheFn = Lookup[oracles, "EvalWithCache", None];
    evalHVPFn = Lookup[oracles, "EvalHVP", None];
    evalHessFn = Lookup[oracles, "EvalHessian", None];

    optName = optsAssoc[Optimizer];
    If[MatchQ[optName, {_String, (_Rule|_RuleDelayed)...}],
        inlineSubOpts = Association[Rest[optName]];
        optName = First[optName]
    ];
    If[optName === Automatic,
        optName = Which[
            evalGradFn =!= None || evalWithCacheFn =!= None, "LBFGS",
            True, "NelderMead"
        ]
    ];
    If[optName === "GradientAscent", optName = "SteepestAscent"];

    targetVal = optsAssoc[TargetObjective];
    hasTargetQ = NumericQ[targetVal];
    If[!hasTargetQ && NumericQ[optsAssoc[TargetFidelity]],
        targetVal = optsAssoc[TargetFidelity];
        hasTargetQ = True
    ];

    subOpts = optsAssoc[OptimizerOptions];
    If[ListQ[subOpts], subOpts = Association[subOpts]];
    If[!AssociationQ[subOpts], subOpts = <||>];
    subOpts = Join[inlineSubOpts, subOpts];

    (* Dispatch Derivative-Free Solvers *)
    If[MemberQ[$DerivativeFreeOptimizers, optName],
        If[evalObjFn === None,
            Message[RunGeneralOptimizer::badoracle, "EvalObj", optName];
            Return[$Failed]
        ];
        Return[
            If[optName === "NelderMead",
                runNelderMeadSolver[evalObjFn, x0Vec, scaleVec, lbVec, ubVec, hasBoundsQ, optsAssoc],
                runPowellSolver[evalObjFn, x0Vec, scaleVec, lbVec, ubVec, hasBoundsQ, optsAssoc]
            ]
        ]
    ];

    (* Dispatch FindMinimum Bridge *)
    If[optName === "FindMinimum",
        If[evalObjFn === None,
            Message[RunGeneralOptimizer::badoracle, "EvalObj", "FindMinimum"];
            Return[$Failed]
        ];
        Return[runFindMinimumBridge[evalObjFn, evalGradFn, x0Vec, scaleVec, lbVec, ubVec, hasBoundsQ, optsAssoc]]
    ];

    (* Build unified GRAPE objective+gradient oracle *)
    If[evalGradFn === None && evalWithCacheFn === None,
        Message[RunGeneralOptimizer::badoracle, "EvalGrad/EvalWithCache", optName];
        Return[$Failed]
    ];

    evalFlat = If[evalGradFn =!= None,
        Function[{xIn},
            Module[{r = evalGradFn[xIn]},
                If[r === $Failed || !AssociationQ[r] || !NumericQ[r["Objective"]] || !KeyExistsQ[r, "Gradient"],
                    $Failed,
                    <|"Objective" -> N[Re[r["Objective"]]],
                      "Fidelity" -> N[Re[Lookup[r, "Fidelity", r["Objective"]]]],
                      "Penalty" -> N[Re[Lookup[r, "Penalty", 0.]]],
                      "Gradient" -> Developer`ToPackedArray[N[Flatten[r["Gradient"]]]],
                      "Extra" -> Lookup[r, "Extra", <||>]|>
                ]
            ]
        ],
        Function[{xIn},
            Module[{r = evalWithCacheFn[xIn]},
                If[r === $Failed || !AssociationQ[r] || !NumericQ[r["Objective"]] || !KeyExistsQ[r, "Gradient"],
                    $Failed,
                    <|"Objective" -> N[Re[r["Objective"]]],
                      "Fidelity" -> N[Re[Lookup[r, "Fidelity", r["Objective"]]]],
                      "Penalty" -> N[Re[Lookup[r, "Penalty", 0.]]],
                      "Gradient" -> Developer`ToPackedArray[N[Flatten[r["Gradient"]]]],
                      "Extra" -> Lookup[r, "Extra", <||>]|>
                ]
            ]
        ]
    ];

    evalObjWrapped = If[evalObjFn === None, None,
        Function[{xIn},
            Module[{r = evalObjFn[xIn]},
                Which[
                    r === $Failed, $Failed,
                    NumericQ[r], <|"Objective" -> N[Re[r]], "Fidelity" -> N[Re[r]], "Penalty" -> 0., "Extra" -> <||>|>,
                    AssociationQ[r] && NumericQ[r["Objective"]],
                        <|"Objective" -> N[Re[r["Objective"]]],
                          "Fidelity" -> N[Re[Lookup[r, "Fidelity", r["Objective"]]]],
                          "Penalty" -> N[Re[Lookup[r, "Penalty", 0.]]],
                          "Extra" -> Lookup[r, "Extra", <||>]|>,
                    True, $Failed
                ]
            ]
        ]
    ];

    lsOpt = optsAssoc[LineSearch];
    If[lsOpt === "ZhangHager", lsOpt = "ZhangHagerNonmonotone"];
    optsForGRAPE = {
        Method -> optName,
        Optimizer -> optName,
        LineSearch -> lsOpt,
        WolfeParameters -> optsAssoc[WolfeParameters],
        LBFGSMemory -> Max[1, Round[optsAssoc[LBFGSMemory]]],
        MaxIterations -> Max[1, Round[optsAssoc[MaxIterations]]],
        MaxLineSearchEvaluations -> Max[1, Round[optsAssoc[MaxLineSearchEvaluations]]],
        GradientTolerance -> N[optsAssoc[GradientTolerance]],
        FunctionTolerance -> N[optsAssoc[FunctionTolerance]],
        StepTolerance -> N[optsAssoc[StepTolerance]],
        TimeConstraint -> optsAssoc[TimeConstraint],
        TargetFidelity -> If[hasTargetQ, N[Re[targetVal]], None],
        RecordHistory -> TrueQ[optsAssoc[RecordHistory]],
        StepMonitor -> optsAssoc[StepMonitor],
        MonitorInterval -> Max[1, Round[optsAssoc[MonitorInterval]]],
        Verbose -> TrueQ[optsAssoc[Verbose]],
        OptimizerOptions -> Normal[subOpts]
    };

    (* Dispatch First-Order GRAPE Solvers via $OCSolverRegistry["RunScaledOptimizer"] *)
    If[MemberQ[$FirstOrderGRAPEOptimizers, optName],
        res = ($OCSolverRegistry["RunScaledOptimizer"])[
            evalFlat, x0Vec, scaleVec, lbVec, ubVec, hasBoundsQ, optsForGRAPE, None, False, evalObjWrapped
        ];
        If[res === $Failed || !AssociationQ[res], Return[$Failed]];
        Return[
            Join[res, <|
                "Status" -> normalizeSolverStatus[Lookup[res, "Status", "Completed"], hasTargetQ],
                "HVPEvaluations" -> 0,
                "HessianEvaluations" -> 0,
                "ScaleVector" -> scaleVec
            |>]
        ]
    ];

    (* Dispatch Registered First-Order Extension (RegularizedLBFGS) *)
    If[optName === "RegularizedLBFGS",
        res = ($OCSolverRegistry["Optimizers"]["RegularizedLBFGS"])[<|
            "EvalFn" -> evalFlat,
            "EvalObjFn" -> evalObjWrapped,
            "X0Vec" -> x0Vec,
            "ScaleVec" -> scaleVec,
            "LBVec" -> lbVec,
            "UBVec" -> ubVec,
            "HasBounds" -> hasBoundsQ,
            "CustomProjFn" -> None,
            "HasCustomConstr" -> False,
            "AnyConstraints" -> hasBoundsQ,
            "OptsList" -> optsForGRAPE,
            "SubOpts" -> subOpts,
            "Comp" -> <||>
        |>];
        If[res === $Failed || !AssociationQ[res], Return[$Failed]];
        Return[
            Join[res, <|
                "Status" -> normalizeSolverStatus[Lookup[res, "Status", "Completed"], hasTargetQ],
                "HVPEvaluations" -> 0,
                "HessianEvaluations" -> 0,
                "ScaleVector" -> scaleVec
            |>]
        ]
    ];

    (* Dispatch Second-Order & Hybrid Solvers via $OCSolverRegistry["Optimizers"] *)
    If[MemberQ[$SecondOrderGRAPEOptimizers, optName],
        If[evalWithCacheFn === None || evalHVPFn === None,
            Message[RunGeneralOptimizer::badoracle, "EvalWithCache/EvalHVP", optName];
            Return[$Failed]
        ];
        If[MemberQ[$DenseHessianOptimizers, optName] && evalHessFn === None,
            Message[RunGeneralOptimizer::badoracle, "EvalHessian", optName];
            Return[$Failed]
        ];

        btMode = Lookup[optsAssoc, BoundTransform, "None"];
        If[btMode === Automatic, btMode = If[hasBoundsQ, "Sigmoid", "None"]];
        btSpec = If[hasBoundsQ && MemberQ[{"Sigmoid", "Log"}, btMode],
            buildBoundTransformSpec[x0Vec, lbVec, ubVec, scaleVec, btMode],
            None
        ];

        If[btSpec =!= None,
            (* Smooth reparameterization u -> x = phi(u) for unconstrained 2nd-order solvers *)
            Module[{uEvalFlat, uEvalObj, uEvalWithCache, uEvalHVP, uEvalHess, uScale = ConstantArray[1., nDim],
                    xFin, gFinX, rFin},
                uEvalObj = If[evalObjWrapped === None, None,
                    Function[{uIn},
                        With[{xMapped = First[evalBoundTransform[btSpec, uIn]]},
                            evalObjWrapped[xMapped]
                        ]
                    ]
                ];
                uEvalFlat = Function[{uIn},
                    Module[{xMapped, dp, dpp, r},
                        {xMapped, dp, dpp} = evalBoundTransform[btSpec, uIn];
                        r = evalFlat[xMapped];
                        If[r === $Failed, $Failed,
                            Append[r, "Gradient" -> Developer`ToPackedArray[r["Gradient"] * dp]]
                        ]
                    ]
                ];
                uEvalWithCache = Function[{uIn},
                    Module[{xMapped, dp, dpp, r, gX},
                        {xMapped, dp, dpp} = evalBoundTransform[btSpec, uIn];
                        r = evalWithCacheFn[xMapped];
                        If[r === $Failed || !AssociationQ[r] || !NumericQ[r["Objective"]] || !KeyExistsQ[r, "Gradient"],
                            $Failed,
                            gX = Developer`ToPackedArray[N[Flatten[r["Gradient"]]]];
                            <|"Objective" -> N[Re[r["Objective"]]],
                              "Fidelity" -> N[Re[Lookup[r, "Fidelity", r["Objective"]]]],
                              "Penalty" -> N[Re[Lookup[r, "Penalty", 0.]]],
                              "Gradient" -> Developer`ToPackedArray[gX * dp],
                              "Extra" -> Lookup[r, "Extra", <||>],
                              "EvalState" -> <|"InnerState" -> r["EvalState"], "GradX" -> gX, "Dp" -> dp, "Dpp" -> dpp|>|>
                        ]
                    ]
                ];
                uEvalHVP = Function[{stU, vU},
                    Module[{dp = stU["Dp"], dpp = stU["Dpp"], gX = stU["GradX"], hX},
                        hX = evalHVPFn[stU["InnerState"], Developer`ToPackedArray[dp * vU]];
                        If[!ArrayQ[hX, _, NumericQ], $Failed,
                            Developer`ToPackedArray[dp * Flatten[hX] + gX * dpp * vU]
                        ]
                    ]
                ];
                uEvalHess = If[evalHessFn === None, None,
                    Function[{stU},
                        Module[{dp = stU["Dp"], dpp = stU["Dpp"], gX = stU["GradX"], hMat},
                            hMat = evalHessFn[stU["InnerState"]];
                            If[!MatrixQ[hMat, NumericQ], $Failed,
                                Developer`ToPackedArray[
                                    Outer[Times, dp, dp] * hMat + DiagonalMatrix[gX * dpp]
                                ]
                            ]
                        ]
                    ]
                ];
                res = ($OCSolverRegistry["Optimizers"][optName])[<|
                    "EvalFn" -> uEvalFlat,
                    "EvalObjFn" -> uEvalObj,
                    "EvalWithCacheFn" -> uEvalWithCache,
                    "EvalHVPFn" -> uEvalHVP,
                    "EvalHessianFn" -> uEvalHess,
                    "X0Vec" -> btSpec["U0"],
                    "ScaleVec" -> uScale,
                    "LBVec" -> ConstantArray[-Infinity, nDim],
                    "UBVec" -> ConstantArray[Infinity, nDim],
                    "HasBounds" -> False,
                    "CustomProjFn" -> None,
                    "HasCustomConstr" -> False,
                    "AnyConstraints" -> False,
                    "OptsList" -> optsForGRAPE,
                    "SubOpts" -> subOpts,
                    "Comp" -> <||>
                |>];
                If[res === $Failed || !AssociationQ[res], Return[$Failed]];
                xFin = First[evalBoundTransform[btSpec, res["ParametersVec"]]];
                rFin = evalFlat[xFin];
                gNormFinal = If[AssociationQ[rFin], Norm[rFin["Gradient"] * (1. / scaleVec)], res["GradientNorm"]];
                Return[
                    Join[res, <|
                        "ParametersVec" -> xFin,
                        "GradientNorm" -> gNormFinal,
                        "Status" -> normalizeSolverStatus[Lookup[res, "Status", "Completed"], hasTargetQ],
                        "ScaleVector" -> scaleVec,
                        "Extra" -> Join[Lookup[res, "Extra", <||>], <|"BoundTransform" -> btMode|>]
                    |>]
                ]
            ],
            (* Unbounded OR active bounds cleanly rejected by SecondOrderOptimizers *)
            res = ($OCSolverRegistry["Optimizers"][optName])[<|
                "EvalFn" -> evalFlat,
                "EvalObjFn" -> evalObjWrapped,
                "EvalWithCacheFn" -> Function[{xIn},
                    Module[{r = evalWithCacheFn[xIn]},
                        If[r === $Failed || !AssociationQ[r] || !NumericQ[r["Objective"]] || !KeyExistsQ[r, "Gradient"],
                            $Failed,
                            <|"Objective" -> N[Re[r["Objective"]]],
                              "Fidelity" -> N[Re[Lookup[r, "Fidelity", r["Objective"]]]],
                              "Penalty" -> N[Re[Lookup[r, "Penalty", 0.]]],
                              "Gradient" -> Developer`ToPackedArray[N[Flatten[r["Gradient"]]]],
                              "Extra" -> Lookup[r, "Extra", <||>],
                              "EvalState" -> r["EvalState"]|>
                        ]
                    ]
                ],
                "EvalHVPFn" -> evalHVPFn,
                "EvalHessianFn" -> evalHessFn,
                "X0Vec" -> x0Vec,
                "ScaleVec" -> scaleVec,
                "LBVec" -> lbVec,
                "UBVec" -> ubVec,
                "HasBounds" -> hasBoundsQ,
                "CustomProjFn" -> None,
                "HasCustomConstr" -> False,
                "AnyConstraints" -> hasBoundsQ,
                "OptsList" -> optsForGRAPE,
                "SubOpts" -> subOpts,
                "Comp" -> <||>
            |>];
            If[res === $Failed || !AssociationQ[res], Return[$Failed]];
            Return[
                Join[res, <|
                    "Status" -> normalizeSolverStatus[Lookup[res, "Status", "Completed"], hasTargetQ],
                    "ScaleVector" -> scaleVec
                |>]
            ]
        ]
    ];

    Message[RunGeneralOptimizer::badopt, "Unknown Optimizer -> " <> ToString[optName]];
    $Failed
];

End[];
EndPackage[];
