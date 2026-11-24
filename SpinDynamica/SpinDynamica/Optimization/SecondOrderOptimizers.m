(* ::Package:: *)

(* SpinDynamica Optimal Control — Explicit Second-Order, Matrix-Free Krylov & Hybrid Optimizers (Phase D) *)

BeginPackage["OptimalControl`SecondOrderOptimizers`", {
    "Global`",
    "SpinDynamica`",
    "OptimalControl`Kernels`",
    "OptimalControl`GRAPE`"
}];

Begin["`Private`"];

(* Compiled Gill-Murray-Wright (1981) pivoted modified Cholesky solver for (P^T L D L^T P) p = -gMin.
   Returns packed vector Join[p, {modFlag}] where modFlag > 0.5 indicates negative curvature modification. *)
modifiedCholeskySolveCompiled = Compile[
    {{hIn, _Real, 2}, {gIn, _Real, 1}, {deltaIn, _Real}},
    Block[{n, cMat, lMat, dVec, perm, yVec, zVec, wVec, rhsPerm, pOut,
           gamma = 0., xi = 0., betaSq, delta, modFlag = 0.,
           j, i, k, q, maxDiag, thetaJ, tmpI, tmpF, sumVal, aij, dJ},
        n = Length[gIn];
        cMat = hIn;
        lMat = Table[0., {n}, {n}];
        dVec = Table[0., {n}];
        perm = Table[i, {i, n}];

        Do[
            tmpF = Abs[hIn[[i, i]]];
            If[tmpF > gamma, gamma = tmpF];
            Do[
                If[k != i,
                    aij = Abs[hIn[[i, k]]];
                    If[aij > xi, xi = aij]
                ],
                {k, n}
            ],
            {i, n}
        ];
        betaSq = Max[gamma, xi / Sqrt[Max[1., N[n * n - 1]]], 1.*^-15];
        delta = Max[deltaIn, 1.*^-14 * Max[1., gamma + xi]];

        Do[
            q = j;
            maxDiag = Abs[cMat[[j, j]]];
            Do[
                If[Abs[cMat[[k, k]]] > maxDiag,
                    maxDiag = Abs[cMat[[k, k]]];
                    q = k
                ],
                {k, j + 1, n}
            ];
            If[q != j,
                tmpI = perm[[j]]; perm[[j]] = perm[[q]]; perm[[q]] = tmpI;
                Do[
                    tmpF = cMat[[j, k]]; cMat[[j, k]] = cMat[[q, k]]; cMat[[q, k]] = tmpF,
                    {k, n}
                ];
                Do[
                    tmpF = cMat[[k, j]]; cMat[[k, j]] = cMat[[k, q]]; cMat[[k, q]] = tmpF,
                    {k, n}
                ];
                Do[
                    tmpF = lMat[[j, k]]; lMat[[j, k]] = lMat[[q, k]]; lMat[[q, k]] = tmpF,
                    {k, j - 1}
                ]
            ];

            lMat[[j, j]] = 1.;
            thetaJ = 0.;
            Do[
                aij = Abs[cMat[[i, j]]];
                If[aij > thetaJ, thetaJ = aij],
                {i, j + 1, n}
            ];
            dJ = Max[delta, Abs[cMat[[j, j]]], (thetaJ * thetaJ) / betaSq];
            If[cMat[[j, j]] < delta || dJ > cMat[[j, j]] + 1.*^-14 * Max[1., dJ],
                modFlag = 1.
            ];
            dVec[[j]] = dJ;

            Do[
                lMat[[i, j]] = cMat[[i, j]] / dJ,
                {i, j + 1, n}
            ];
            Do[
                Do[
                    cMat[[i, k]] -= lMat[[i, j]] * cMat[[k, j]];
                    cMat[[k, i]] = cMat[[i, k]],
                    {k, i, n}
                ],
                {i, j + 1, n}
            ],
            {j, n}
        ];

        rhsPerm = Table[-gIn[[perm[[i]]]], {i, n}];
        yVec = Table[0., {n}];
        Do[
            sumVal = rhsPerm[[i]];
            Do[sumVal -= lMat[[i, k]] * yVec[[k]], {k, i - 1}];
            yVec[[i]] = sumVal,
            {i, n}
        ];
        zVec = Table[yVec[[i]] / dVec[[i]], {i, n}];
        wVec = Table[0., {n}];
        Do[
            sumVal = zVec[[i]];
            Do[sumVal -= lMat[[k, i]] * wVec[[k]], {k, i + 1, n}];
            wVec[[i]] = sumVal,
            {i, n, 1, -1}
        ];
        pOut = Table[0., {n + 1}];
        Do[pOut[[perm[[i]]]] = wVec[[i]], {i, n}];
        pOut[[n + 1]] = modFlag;
        pOut
    ],
    CompilationTarget -> "WVM",
    RuntimeOptions -> "Speed"
];

ModifiedCholeskySolve[hMat_?MatrixQ, gMin_List, deltaIn_: 1.*^-8] := Module[
    {hSym, out, delta = N[deltaIn]},
    If[!NumericQ[delta], Return[$Failed]];
    hSym = Developer`ToPackedArray[N[Re[0.5 * (hMat + Transpose[hMat])]]];
    out = modifiedCholeskySolveCompiled[hSym, Developer`ToPackedArray[N[Re[gMin]]], N[delta]];
    <|
        "Direction" -> Developer`ToPackedArray[Most[out]],
        "Modified" -> TrueQ[Last[out] > 0.5]
    |>
];

(* Exact Rational Function Optimization (RFO) secular solver:
   Finds mu >= max(0, -lamMin) solving mu - s^2 Sum_i gTilde_i^2 / (evals_i + mu) = 0,
   equivalent to the lowest eigenpair of the (n+1)x(n+1) Goodwin-Kuprov augmented Hessian.
   Uses shifted variable u = mu - muLB > 0 so (evals + mu) = dLam + u never cancels near lamMin < 0. *)
SolveRFOSecular[evalsIn_List, evecsIn_?MatrixQ, gMinIn_List, sScaleIn_?NumericQ] := Module[
    {evals, evecs, gMin, s, sSq, gTilde, gNorm, lamMin, lamMaxAbs, muLB, dLam, minTol,
     restIdx, minIdx, restDLam, restG, sumRest, uLo, uHi, u, mu, fVal, dfVal, uNew,
     pRest, qMin, tau, p},
    evals = Developer`ToPackedArray[N[Re[evalsIn]]];
    evecs = Developer`ToPackedArray[N[Re[evecsIn]]];
    gMin = Developer`ToPackedArray[N[Re[gMinIn]]];
    s = Max[N[sScaleIn], 1.*^-12];
    sSq = s * s;
    gTilde = Developer`ToPackedArray[evecs . gMin];
    gNorm = Norm[gMin];
    If[TrueQ[gNorm == 0.] && TrueQ[Min[evals] >= 0.],
        Return[<|"Direction" -> ConstantArray[0., Length[gMin]], "Shift" -> 0., "HardCase" -> False, "NegativeCurvature" -> False|>]
    ];

    lamMin = Min[evals];
    lamMaxAbs = Max[Abs[evals]];
    muLB = Max[0., -lamMin];
    dLam = Developer`ToPackedArray[If[TrueQ[lamMin < 0.], evals - lamMin, evals]];
    minTol = 1.*^-10 * Max[1., lamMaxAbs];
    restIdx = Flatten[Position[dLam, x_ /; TrueQ[x > minTol + If[TrueQ[lamMin >= 0.], lamMin, 0.]], {1}]];
    If[TrueQ[lamMin >= 0.], restIdx = Flatten[Position[evals, x_ /; TrueQ[x > lamMin + minTol], {1}]]];
    minIdx = First[Ordering[evals, 1]];

    (* Saddle / hard case when lamMin < 0 and gMin is orthogonal to min eigenspace *)
    If[TrueQ[lamMin < -1.*^-12] && TrueQ[Sqrt[Total[gTilde[[Complement[Range[Length[evals]], restIdx]]]^2]] <= 1.*^-12 * Max[1., gNorm]],
        restDLam = Clip[dLam[[restIdx]], {1.*^-300, Infinity}];
        restG = gTilde[[restIdx]];
        sumRest = If[restIdx === {}, 0., sSq * Total[restG^2 / restDLam]];
        If[TrueQ[sumRest < muLB],
            pRest = If[restIdx === {}, ConstantArray[0., Length[gMin]], -((restG / restDLam) . evecs[[restIdx]])];
            qMin = evecs[[minIdx]];
            If[TrueQ[gMin . qMin > 0.], qMin = -qMin];
            tau = Clip[Sqrt[Max[0., (muLB - sumRest) / Max[1.*^-12, sSq * Abs[lamMin]]]], {0.05, 1.0}];
            p = Developer`ToPackedArray[pRest + tau * qMin];
            Return[<|"Direction" -> p, "Shift" -> muLB, "HardCase" -> True, "NegativeCurvature" -> True|>]
        ]
    ];

    If[TrueQ[gNorm == 0.] && TrueQ[lamMin < -1.*^-12],
        qMin = evecs[[minIdx]];
        Return[<|"Direction" -> Developer`ToPackedArray[0.05 * qMin], "Shift" -> muLB, "HardCase" -> True, "NegativeCurvature" -> True|>]
    ];

    (* Cancellation-free monotone secular root in u = mu - muLB > 0 *)
    uLo = 0.;
    uHi = s * gNorm + 1.*^-11 * Max[1., lamMaxAbs];
    u = 0.5 * (uLo + uHi);
    Do[
        With[{denom = Clip[dLam + u, {1.*^-300, Infinity}]},
            fVal = (muLB + u) - sSq * Total[(gTilde * gTilde) / denom];
            dfVal = 1. + sSq * Total[(gTilde * gTilde) / (denom * denom)]
        ];
        If[TrueQ[Abs[fVal] <= 1.*^-14 * Max[1., muLB + u]] || TrueQ[uHi - uLo <= 1.*^-14 * Max[u, 1.*^-300]],
            Break[]
        ];
        If[TrueQ[fVal < 0.], uLo = u, uHi = u];
        uNew = u - fVal / dfVal;
        If[TrueQ[uLo < uNew < uHi],
            u = uNew,
            u = If[TrueQ[uLo > 0.], Sqrt[uLo * uHi], 0.1 * uHi]
        ],
        {100}
    ];
    mu = muLB + u;
    p = Developer`ToPackedArray[-((gTilde / Clip[dLam + u, {1.*^-300, Infinity}]) . evecs)];
    <|
        "Direction" -> p,
        "Shift" -> mu,
        "HardCase" -> False,
        "NegativeCurvature" -> TrueQ[lamMin < -1.*^-12]
    |>
];

(* Exact More-Sorensen (1983) trust-region subproblem solver on a small j x j symmetric matrix:
   min_{||h|| <= delta} gSmall^T h + 1/2 h^T tMat h. Uses u = lam - lamLB > 0 to avoid cancellation. *)
SolveSmallTrustRegionSubproblem[tMatIn_?MatrixQ, gSmallIn_List, deltaIn_?NumericQ] := Module[
    {tSym, gSmall, delta, evals, evecs, gTilde, gNorm, lamMin, lamMaxAbs, hInt,
     lamLB, dLam, minTol, restIdx, minIdx, restDLam, restG, hRest, tau, qMin, h,
     uLo, uHi, u, lam, s2, s3, hNorm, phiMS, uNew, modelVal},
    tSym = Developer`ToPackedArray[N[Re[0.5 * (tMatIn + Transpose[tMatIn])]]];
    gSmall = Developer`ToPackedArray[N[Re[gSmallIn]]];
    delta = Max[N[deltaIn], 1.*^-16];
    gNorm = Norm[gSmall];
    {evals, evecs} = Eigensystem[tSym];
    evals = Developer`ToPackedArray[N[Re[evals]]];
    evecs = Developer`ToPackedArray[N[Re[evecs]]];
    gTilde = Developer`ToPackedArray[evecs . gSmall];
    lamMin = Min[evals];
    lamMaxAbs = Max[Abs[evals]];

    If[TrueQ[lamMin > 1.*^-13 * Max[1., lamMaxAbs]],
        hInt = Developer`ToPackedArray[-((gTilde / evals) . evecs)];
        If[TrueQ[Norm[hInt] <= delta],
            modelVal = gSmall . hInt + 0.5 * (hInt . (tSym . hInt));
            Return[<|
                "Step" -> hInt,
                "Lambda" -> 0.,
                "OnBoundary" -> False,
                "ModelValue" -> modelVal,
                "NegativeCurvature" -> False
            |>]
        ]
    ];

    lamLB = Max[0., -lamMin];
    dLam = Developer`ToPackedArray[If[TrueQ[lamMin < 0.], evals - lamMin, evals]];
    minTol = 1.*^-10 * Max[1., lamMaxAbs];
    restIdx = Flatten[Position[evals, x_ /; TrueQ[x > lamMin + minTol], {1}]];
    minIdx = First[Ordering[evals, 1]];

    (* Check hard case when lamMin <= 0 and min eigenspace is orthogonal to gSmall *)
    If[TrueQ[lamMin <= -1.*^-12] && TrueQ[Norm[gTilde[[Complement[Range[Length[evals]], restIdx]]]] <= 1.*^-12 * Max[1., gNorm]],
        restDLam = Clip[dLam[[restIdx]], {1.*^-300, Infinity}];
        restG = gTilde[[restIdx]];
        hRest = If[restIdx === {}, ConstantArray[0., Length[gSmall]], -((restG / restDLam) . evecs[[restIdx]])];
        If[TrueQ[Norm[hRest] <= delta],
            tau = Sqrt[Max[0., delta * delta - (hRest . hRest)]];
            qMin = evecs[[minIdx]];
            If[TrueQ[gSmall . qMin > 0.], qMin = -qMin];
            h = Developer`ToPackedArray[hRest + tau * qMin];
            modelVal = gSmall . h + 0.5 * (h . (tSym . h));
            Return[<|
                "Step" -> h,
                "Lambda" -> lamLB,
                "OnBoundary" -> True,
                "ModelValue" -> modelVal,
                "NegativeCurvature" -> True
            |>]
        ]
    ];

    If[TrueQ[gNorm == 0.],
        qMin = evecs[[minIdx]];
        h = Developer`ToPackedArray[delta * qMin];
        modelVal = 0.5 * (h . (tSym . h));
        Return[<|
            "Step" -> h,
            "Lambda" -> lamLB,
            "OnBoundary" -> True,
            "ModelValue" -> modelVal,
            "NegativeCurvature" -> TrueQ[lamMin < -1.*^-12]
        |>]
    ];

    uLo = 0.;
    uHi = gNorm / delta + 1.*^-10 * Max[1., lamMaxAbs];
    u = 0.5 * (uLo + uHi);
    Do[
        With[{denom = Clip[dLam + u, {1.*^-300, Infinity}]},
            s2 = Total[(gTilde * gTilde) / (denom * denom)];
            s3 = Total[(gTilde * gTilde) / (denom * denom * denom)]
        ];
        hNorm = Sqrt[s2];
        If[TrueQ[Abs[hNorm - delta] <= 1.*^-13 * delta] || TrueQ[uHi - uLo <= 1.*^-14 * Max[u, 1.*^-300]],
            Break[]
        ];
        phiMS = (1. / Max[hNorm, 1.*^-30]) - (1. / delta);
        If[TrueQ[phiMS < 0.], uLo = u, uHi = u];
        uNew = u + (s2 / Max[s3, 1.*^-300]) * ((hNorm / delta) - 1.);
        If[TrueQ[uLo < uNew < uHi],
            u = uNew,
            u = If[TrueQ[uLo > 0.], Sqrt[uLo * uHi], 0.1 * uHi]
        ],
        {100}
    ];
    lam = lamLB + u;
    h = Developer`ToPackedArray[-((gTilde / Clip[dLam + u, {1.*^-300, Infinity}]) . evecs)];
    modelVal = gSmall . h + 0.5 * (h . (tSym . h));
    <|
        "Step" -> h,
        "Lambda" -> lam,
        "OnBoundary" -> True,
        "ModelValue" -> modelVal,
        "NegativeCurvature" -> TrueQ[lamMin < -1.*^-12]
    |>
];

(* Exact Cartis-Gould-Toint (2011) Adaptive Cubic Regularization subproblem solver:
   min_{h} gSmall^T h + 1/2 h^T tMat h + (sigma / 3) ||h||^3. Uses u = lam - lamLB > 0. *)
SolveSmallARCSubproblem[tMatIn_?MatrixQ, gSmallIn_List, sigmaIn_?NumericQ] := Module[
    {tSym, gSmall, sigma, evals, evecs, gTilde, gNorm, lamMin, lamMaxAbs,
     lamLB, dLam, minTol, restIdx, minIdx, restDLam, restG, hRest, sRest, tau, qMin, h,
     uLo, uHi, u, lam, s2, s3, hNorm, fARC, dfARC, uNew, qMod, cMod},
    tSym = Developer`ToPackedArray[N[Re[0.5 * (tMatIn + Transpose[tMatIn])]]];
    gSmall = Developer`ToPackedArray[N[Re[gSmallIn]]];
    sigma = Max[N[sigmaIn], 1.*^-12];
    gNorm = Norm[gSmall];
    {evals, evecs} = Eigensystem[tSym];
    evals = Developer`ToPackedArray[N[Re[evals]]];
    evecs = Developer`ToPackedArray[N[Re[evecs]]];
    gTilde = Developer`ToPackedArray[evecs . gSmall];
    lamMin = Min[evals];
    lamMaxAbs = Max[Abs[evals]];
    lamLB = Max[0., -lamMin];
    dLam = Developer`ToPackedArray[If[TrueQ[lamMin < 0.], evals - lamMin, evals]];
    minTol = 1.*^-10 * Max[1., lamMaxAbs];
    restIdx = Flatten[Position[evals, x_ /; TrueQ[x > lamMin + minTol], {1}]];
    minIdx = First[Ordering[evals, 1]];

    If[TrueQ[lamMin <= -1.*^-12] && TrueQ[Norm[gTilde[[Complement[Range[Length[evals]], restIdx]]]] <= 1.*^-12 * Max[1., gNorm]],
        restDLam = Clip[dLam[[restIdx]], {1.*^-300, Infinity}];
        restG = gTilde[[restIdx]];
        hRest = If[restIdx === {}, ConstantArray[0., Length[gSmall]], -((restG / restDLam) . evecs[[restIdx]])];
        sRest = hRest . hRest;
        If[TrueQ[sigma * Sqrt[sRest] <= lamLB],
            tau = Sqrt[Max[0., (lamLB / sigma)^2 - sRest]];
            qMin = evecs[[minIdx]];
            If[TrueQ[gSmall . qMin > 0.], qMin = -qMin];
            h = Developer`ToPackedArray[hRest + tau * qMin];
            qMod = gSmall . h + 0.5 * (h . (tSym . h));
            cMod = qMod + (sigma / 3.) * (Norm[h]^3);
            Return[<|
                "Step" -> h,
                "Lambda" -> lamLB,
                "QuadModelValue" -> qMod,
                "CubicModelValue" -> cMod,
                "NegativeCurvature" -> True
            |>]
        ]
    ];

    If[TrueQ[gNorm == 0.],
        If[TrueQ[lamMin < -1.*^-12],
            qMin = evecs[[minIdx]];
            h = Developer`ToPackedArray[(lamLB / sigma) * qMin];
            qMod = 0.5 * (h . (tSym . h));
            cMod = qMod + (sigma / 3.) * (Norm[h]^3);
            Return[<|"Step" -> h, "Lambda" -> lamLB, "QuadModelValue" -> qMod, "CubicModelValue" -> cMod, "NegativeCurvature" -> True|>],
            Return[<|"Step" -> ConstantArray[0., Length[gSmall]], "Lambda" -> 0., "QuadModelValue" -> 0., "CubicModelValue" -> 0., "NegativeCurvature" -> False|>]
        ]
    ];

    uLo = 0.;
    uHi = Sqrt[sigma * gNorm] + 1.*^-10 * Max[1., lamMaxAbs];
    u = 0.5 * (uLo + uHi);
    Do[
        With[{denom = Clip[dLam + u, {1.*^-300, Infinity}]},
            s2 = Total[(gTilde * gTilde) / (denom * denom)];
            s3 = Total[(gTilde * gTilde) / (denom * denom * denom)]
        ];
        hNorm = Sqrt[s2];
        fARC = (lamLB + u) - sigma * hNorm;
        If[TrueQ[Abs[fARC] <= 1.*^-13 * Max[1., lamLB + u]] || TrueQ[uHi - uLo <= 1.*^-14 * Max[u, 1.*^-300]],
            Break[]
        ];
        If[TrueQ[fARC < 0.], uLo = u, uHi = u];
        dfARC = 1. + sigma * s3 / Max[hNorm, 1.*^-300];
        uNew = u - fARC / dfARC;
        If[TrueQ[uLo < uNew < uHi],
            u = uNew,
            u = If[TrueQ[uLo > 0.], Sqrt[uLo * uHi], 0.1 * uHi]
        ],
        {100}
    ];
    lam = lamLB + u;
    h = Developer`ToPackedArray[-((gTilde / Clip[dLam + u, {1.*^-300, Infinity}]) . evecs)];
    qMod = gSmall . h + 0.5 * (h . (tSym . h));
    cMod = qMod + (sigma / 3.) * (Norm[h]^3);
    <|
        "Step" -> h,
        "Lambda" -> lam,
        "QuadModelValue" -> qMod,
        "CubicModelValue" -> cMod,
        "NegativeCurvature" -> TrueQ[lamMin < -1.*^-12]
    |>
];

(* Standard L-BFGS two-loop inverse-Hessian-vector product for Morales-Nocedal TR-CG preconditioning *)
lbfgsTwoLoop[sList_List, yList_List, rhoList_List, vec_List] := Module[
    {mCur = Length[sList], alpha, qVec, syLast, yyLast, gamma, rVec, idx, beta},
    If[mCur === 0, Return[Developer`ToPackedArray[N[vec]]]];
    qVec = Developer`ToPackedArray[N[vec]];
    alpha = ConstantArray[0., mCur];
    Do[
        alpha[[idx]] = rhoList[[idx]] * (sList[[idx]] . qVec);
        qVec = qVec - alpha[[idx]] * yList[[idx]],
        {idx, mCur, 1, -1}
    ];
    syLast = Last[sList] . Last[yList];
    yyLast = Last[yList] . Last[yList];
    gamma = If[TrueQ[yyLast > 1.*^-24 && syLast > 0.], syLast / yyLast, 1.];
    rVec = gamma * qVec;
    Do[
        beta = rhoList[[idx]] * (yList[[idx]] . rVec);
        rVec = rVec + sList[[idx]] * (alpha[[idx]] - beta),
        {idx, 1, mCur}
    ];
    Developer`ToPackedArray[rVec]
];

(* Shared second-order context builder with capability guards, analytical oracle hooks, and scaled wrappers *)
buildSecondOrderContext[ctx_Association, optName_String, isDenseHessian_?BooleanQ] := Module[
    {optsList = ctx["OptsList"], subOpts = Lookup[ctx, "SubOpts", <||>],
     x0Vec = ctx["X0Vec"], scaleVec, invS, nDim, comp = Lookup[ctx, "Comp", <||>],
     hasCustomOracles, caps, maxHessDim, allowFastFallback},
    nDim = Length[x0Vec];
    scaleVec = Developer`ToPackedArray[N[ctx["ScaleVec"]]];
    invS = Developer`ToPackedArray[1. / scaleVec];
    hasCustomOracles = KeyExistsQ[ctx, "EvalWithCacheFn"] && KeyExistsQ[ctx, "EvalHVPFn"];

    If[!hasCustomOracles,
        If[!AssociationQ[comp] || !KeyExistsQ[$OCSolverRegistry, "SecondOrderCapabilities"],
            Message[OptimalControlProblem::unsupported2nd, optName, "Compiled second-order capability inspector unavailable."];
            Return[$Failed]
        ];
        caps = ($OCSolverRegistry["SecondOrderCapabilities"])[comp];
        If[!TrueQ[Lookup[caps, "SupportsSecondOrderQ", False]],
            Message[OptimalControlProblem::unsupported2nd, optName, Lookup[caps, "Reason", "Unsupported problem specification."]];
            Return[$Failed]
        ];
        If[TrueQ[Lookup[ctx, "AnyConstraints", False]],
            Message[OptimalControlProblem::unconstronly, optName];
            Return[$Failed]
        ];
        If[isDenseHessian,
            maxHessDim = Max[1, Round[Lookup[subOpts, "MaxHessianDim", 500]]];
            If[nDim > maxHessDim,
                Message[OptimalControlProblem::unsupported2nd, optName,
                    "Parameter dimension " <> ToString[nDim] <> " exceeds MaxHessianDim -> " <> ToString[maxHessDim] <>
                    ". Use \"TrustRegionNewtonCG\", \"TrustRegionGLTR\", or \"AdaptiveCubicRegularization\" (or raise \"MaxHessianDim\")."];
                Return[$Failed]
            ];
            allowFastFallback = Lookup[subOpts, "AllowColumnHVPHessian", True];
            If[!TrueQ[allowFastFallback] && !TrueQ[Lookup[caps, "SupportsFastExplicitHessianQ", False]],
                Message[OptimalControlProblem::unsupported2nd, optName,
                    "Fast explicit Goodwin-Kuprov Hessian is available only for Method -> \"Spectral\"; set \"AllowColumnHVPHessian\" -> True or use matrix-free Krylov solvers."];
                Return[$Failed]
            ]
        ],
        (
            If[TrueQ[Lookup[ctx, "AnyConstraints", False]],
                Message[OptimalControlProblem::unconstronly, optName];
                Return[$Failed]
            ];
            If[isDenseHessian,
                maxHessDim = Max[1, Round[Lookup[subOpts, "MaxHessianDim", 500]]];
                If[nDim > maxHessDim,
                    Message[OptimalControlProblem::unsupported2nd, optName,
                        "Parameter dimension " <> ToString[nDim] <> " exceeds MaxHessianDim -> " <> ToString[maxHessDim] <>
                        ". Use \"TrustRegionNewtonCG\", \"TrustRegionGLTR\", or \"AdaptiveCubicRegularization\" (or raise \"MaxHessianDim\")."];
                    Return[$Failed]
                ];
                If[!KeyExistsQ[ctx, "EvalHessianFn"] || ctx["EvalHessianFn"] === None,
                    Message[OptimalControlProblem::unsupported2nd, optName, "Custom analytical Hessian oracle missing."];
                    Return[$Failed]
                ]
            ]
        )
    ];

    <|
        "HasCustomOracles" -> hasCustomOracles,
        "ScaleVec" -> scaleVec,
        "InvScale" -> invS,
        "NumDim" -> nDim,
        "Comp" -> comp,
        "OptsList" -> optsList,
        "SubOpts" -> subOpts
    |>
];

(* Shared line search helper for explicit Hessian Newton methods *)
runNewtonLineSearch[evalCacheScaled_, th_, fVal_?NumericQ, gMax_, dir_, lsCfg_Association, tolScale_?NumericQ, nmStateIn_] := Module[
    {lsMode = lsCfg["Mode"], c1 = lsCfg["C1"], c2 = lsCfg["C2"], maxLs = lsCfg["MaxEvals"],
     lsExtra = lsCfg["ExtraOpts"], regLS, step = 1.0, phi0 = -fVal, dphi0,
     thTrial, resTrial, fTrial, phiTrial, tolEff, d0, denom, stepQ, trials = 0, swOut, nmOut},
    dphi0 = -(gMax . dir);
    If[!TrueQ[dphi0 < 0.], Return[$Failed]];

    If[lsMode === "Armijo",
        Do[
            trials++;
            thTrial = Developer`ToPackedArray[th + step * dir];
            resTrial = evalCacheScaled[thTrial];
            If[resTrial === $Failed || !NumericQ[resTrial["Objective"]] || !VectorQ[resTrial["Gradient"], NumericQ],
                step *= 0.5;
                If[step < 1.*^-13, Break[]];
                Continue[]
            ];
            fTrial = resTrial["Objective"];
            phiTrial = -fTrial;
            tolEff = Min[tolScale, 0.25 * Abs[c1 * step * dphi0]];
            If[TrueQ[phiTrial <= phi0 + c1 * step * dphi0 + tolEff],
                Return[<|
                    "Accepted" -> True,
                    "Theta" -> thTrial,
                    "Res" -> resTrial,
                    "StepLength" -> step,
                    "StrongWolfeSatisfied" -> False,
                    "ApproxWolfeSatisfied" -> False,
                    "Trials" -> trials,
                    "NMState" -> nmStateIn
                |>, Module]
            ];
            denom = 2. * (phiTrial - phi0 - dphi0 * step);
            stepQ = If[TrueQ[denom > 1.*^-24], -dphi0 * step^2 / denom, 0.5 * step];
            step = Clip[stepQ, {0.1 * step, 0.5 * step}];
            If[step < 1.*^-13, Break[]],
            {Max[maxLs, 25]}
        ];
        Return[$Failed]
    ];

    If[lsMode === "StrongWolfe",
        swOut = ($OCSolverRegistry["StrongWolfeLineSearch"])[
            evalCacheScaled, th, fVal, gMax, dir, 1.0, c1, c2, maxLs, tolScale
        ];
        If[AssociationQ[swOut],
            Return[Append[swOut, "NMState" -> nmStateIn]],
            Return[$Failed]
        ]
    ];

    regLS = Lookup[$OCSolverRegistry, "LineSearches", <||>];
    If[KeyExistsQ[regLS, lsMode],
        nmOut = (regLS[lsMode])[<|
            "EvalScaled" -> evalCacheScaled,
            "EvalObjScaled" -> None,
            "ProjBox" -> Identity,
            "HasConstraints" -> False,
            "Theta" -> th,
            "FVal" -> fVal,
            "GVal" -> gMax,
            "Dir" -> dir,
            "InitStep" -> 1.0,
            "C1" -> c1,
            "C2" -> c2,
            "MaxEvals" -> maxLs,
            "TolScale" -> tolScale,
            "PrevDeltaF" -> 0.,
            "NMState" -> nmStateIn,
            "ExtraOpts" -> lsExtra
        |>];
        Return[nmOut]
    ];

    $Failed
];

(* Core driver for Explicit Hessian Newton family: NewtonRFO, NewtonShift, NewtonModifiedCholesky, NewtonSpectral *)
runExplicitNewtonFamily[ctx_Association, variant_String] := Module[
    {soCtx, evalFlatFn = ctx["EvalFn"], evalWithCacheOracle, evalHessianOracle,
     x0Vec = ctx["X0Vec"], optsList = ctx["OptsList"], subOpts = Lookup[ctx, "SubOpts", <||>],
     scaleVec, invS, comp, hasOracles, nDim,
     maxIter, targetFid, fidCrit, gradTol, funTol, stepTol, tLimit, recHist, verbose,
     lsOpt, wolfeOpt, maxLsOpt, monOpt, monIntOpt, monCfg, stageTag, lastMonTime = 0., lastMonIter = -1,
     parseMonFn, emitMonFn, evalCritFn, resolveLSFn, lsCfg, explicitLSQ, nmState = <||>,
     tStart, dynTime = 0., evals = 0, objEvals = 0, gradEvals = 0, hvpEvals = 0, hessEvals = 0,
     acceptedSteps = 0, rejectedSteps = 0, negCurvCount = 0, swCount = 0, approxWolfeCnt = 0,
     evalCacheScaled, computeHessMinScaled, th, res0, fVal, phiVal, fidVal, penVal, gMax, gMin,
     gNorm, extra, bestFid, evalState, objHist, fidHist, gHist, nIter = 0, status = "MaxIterationsReached",
     lastStep = 0., fOld, stopNow, tolScale, tolEff,
     sRFO, sRFOMin, sRFOMax, muShift, muMin, muMax, deltaMC, eigFloorScale, maxStepNorm,
     hMin, eigSys, evalsH, evecsH, gTilde, minIdx, lamMin, maxAbsH,
     accepted, thNew, resNew, fNew, phiNew, dir, p, predDec, ared, rho,
     sVec, chol, hScale, lamAbs, mcOut, rfoOut, shiftTry, tauEsc, qMin, lsOut, maxInner},

    soCtx = buildSecondOrderContext[ctx, variant, True];
    If[soCtx === $Failed, Return[$Failed]];
    scaleVec = soCtx["ScaleVec"];
    invS = soCtx["InvScale"];
    nDim = soCtx["NumDim"];
    comp = soCtx["Comp"];
    hasOracles = soCtx["HasCustomOracles"];

    evalWithCacheOracle = If[hasOracles, ctx["EvalWithCacheFn"], $OCSolverRegistry["EvaluateWithCache"]];
    evalHessianOracle = If[hasOracles, ctx["EvalHessianFn"], $OCSolverRegistry["EvaluateHessian"]];

    {maxIter, targetFid, gradTol, funTol, stepTol, tLimit, recHist, verbose, lsOpt, wolfeOpt, maxLsOpt, monOpt, monIntOpt} =
        {MaxIterations, TargetFidelity, GradientTolerance, FunctionTolerance, StepTolerance,
         TimeConstraint, RecordHistory, Verbose, LineSearch, WolfeParameters, MaxLineSearchEvaluations,
         StepMonitor, MonitorInterval} /. optsList /. Options[GRAPEOptimize];
    fidCrit = Lookup[FilterRules[optsList, {$OCSolverRegistry["FidelityCriterionKey"]}], $OCSolverRegistry["FidelityCriterionKey"], Automatic];
    stageTag = Lookup[FilterRules[optsList, {$OCSolverRegistry["MonitorStageKey"]}], $OCSolverRegistry["MonitorStageKey"], "Optimization"];
    parseMonFn = $OCSolverRegistry["ParseMonitorConfig"];
    emitMonFn = $OCSolverRegistry["EmitMonitorIfDue"];
    evalCritFn = $OCSolverRegistry["EvalCriterionFidelity"];
    resolveLSFn = $OCSolverRegistry["ResolveLineSearchConfig"];

    monCfg = parseMonFn[monOpt, monIntOpt];
    If[monCfg === $Failed, Return[$Failed]];
    explicitLSQ = lsOpt =!= Automatic;
    lsCfg = resolveLSFn[If[lsOpt === Automatic, "Armijo", lsOpt], wolfeOpt, maxLsOpt, "BFGS", False];
    If[lsCfg === $Failed, Return[$Failed]];

    sRFOMin = Max[N[Lookup[subOpts, "MinRFOScale", 1.0]], 1.*^-4];
    sRFOMax = Max[sRFOMin, N[Lookup[subOpts, "MaxRFOScale", 1.*^8]]];
    sRFO = Clip[N[Lookup[subOpts, "InitialRFOScale", 1.0]], {sRFOMin, sRFOMax}];
    muMin = Max[N[Lookup[subOpts, "MinShift", 1.*^-12]], 1.*^-16];
    muMax = N[Lookup[subOpts, "MaxShift", 1.*^12]];
    muShift = Max[0., N[Lookup[subOpts, "InitialShift", 0.]]];
    deltaMC = Max[N[Lookup[subOpts, "ModifiedCholeskyDelta", 1.*^-8]], 1.*^-14];
    eigFloorScale = Max[N[Lookup[subOpts, "EigenvalueFloor", 1.*^-8]], 1.*^-14];
    maxInner = Max[5, Round[Lookup[subOpts, "MaxInnerTrials", 12]]];
    maxStepNorm = N[Lookup[subOpts, "MaxStepNorm", Infinity]];

    evalCacheScaled[thIn_] := Module[{t0Dyn = AbsoluteTime[], xIn, r, gFlat},
        xIn = Developer`ToPackedArray[thIn * invS];
        r = If[hasOracles, evalWithCacheOracle[xIn], evalWithCacheOracle[comp, xIn]];
        dynTime += AbsoluteTime[] - t0Dyn;
        evals++; objEvals++;
        If[r === $Failed || !AssociationQ[r] || !NumericQ[r["Objective"]] || !KeyExistsQ[r, "Gradient"],
            Return[$Failed]
        ];
        gradEvals++;
        gFlat = Developer`ToPackedArray[Flatten[r["Gradient"]] * invS];
        <|
            "Objective" -> N[Re[r["Objective"]]],
            "Fidelity" -> N[Re[r["Fidelity"]]],
            "Penalty" -> N[Re[Lookup[r, "Penalty", 0.]]],
            "Gradient" -> gFlat,
            "Extra" -> Lookup[r, "Extra", <||>],
            "EvalState" -> Lookup[r, "EvalState", None]
        |>
    ];

    computeHessMinScaled[evalSt_] := Module[{t0Dyn = AbsoluteTime[], hX, scaledH},
        hX = If[hasOracles,
            evalHessianOracle[evalSt],
            evalHessianOracle[comp, evalSt]
        ];
        dynTime += AbsoluteTime[] - t0Dyn;
        hessEvals++;
        If[!MatrixQ[hX, NumericQ], Return[$Failed]];
        scaledH = -(invS * ToPackedRealMatrix[hX]);
        Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ (invS * Transpose[scaledH])]
    ];

    tStart = AbsoluteTime[];
    th = Developer`ToPackedArray[N[x0Vec * scaleVec]];
    res0 = evalCacheScaled[th];
    If[res0 === $Failed, Return[$Failed]];
    {fVal, fidVal, penVal, gMax, extra, evalState} =
        {res0["Objective"], res0["Fidelity"], res0["Penalty"], res0["Gradient"], res0["Extra"], res0["EvalState"]};
    phiVal = -fVal;
    gMin = Developer`ToPackedArray[-gMax];
    gNorm = Norm[gMin];
    bestFid = fidVal;
    If[TrueQ[maxStepNorm === Infinity],
        maxStepNorm = Max[50., 25. * Max[1., Norm[th]]]
    ];

    objHist = If[recHist, {fVal}, {}];
    fidHist = If[recHist, {fidVal}, {}];
    gHist = If[recHist, {gNorm}, {}];

    With[{critFid0 = evalCritFn[fidVal, extra, fidCrit]},
        status = Which[
            NumericQ[targetFid] && TrueQ[critFid0 >= targetFid], "TargetFidelityReached",
            TrueQ[gNorm <= gradTol], "GradientToleranceReached",
            True, "Running"
        ]
    ];
    {lastMonTime, lastMonIter} = emitMonFn[monCfg, 0, evals, fVal, fidVal, bestFid, gNorm, 0., 0., status, extra, stageTag, lastMonTime, lastMonIter];
    If[status =!= "Running",
        With[{elapsed = AbsoluteTime[] - tStart},
            Return[<|
                "ParametersVec" -> Developer`ToPackedArray[th * invS],
                "Objective" -> fVal, "Fidelity" -> fidVal, "Penalty" -> penVal,
                "GradientNorm" -> gNorm, "ScaledGradient" -> gMax, "Optimizer" -> variant,
                "Iterations" -> 0, "Evaluations" -> evals, "ObjectiveEvaluations" -> objEvals,
                "GradientEvaluations" -> gradEvals, "HVPEvaluations" -> hvpEvals, "HessianEvaluations" -> hessEvals,
                "CGIterations" -> 0, "NegativeCurvatureCount" -> negCurvCount,
                "AcceptedSteps" -> 0, "RejectedSteps" -> 0, "SkippedCurvaturePairs" -> 0,
                "LineSearchMethod" -> lsCfg["Mode"], "StrongWolfeSatisfiedCount" -> 0, "ApproxWolfeCount" -> 0,
                "LastStepLength" -> 0., "Status" -> status,
                "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
                "FidelityHistory" -> Developer`ToPackedArray[fidHist],
                "GradientNormHistory" -> Developer`ToPackedArray[gHist],
                "ElapsedTime" -> elapsed, "DynamicsTime" -> dynTime, "OptimizerTime" -> Max[0., elapsed - dynTime],
                "Extra" -> extra
            |>]
        ]
    ];
    status = "MaxIterationsReached";

    Do[
        hMin = computeHessMinScaled[evalState];
        If[hMin === $Failed, status = "LineSearchTerminated"; Break[]];
        tolScale = 1.*^-14 * Max[1., Abs[fVal]];
        accepted = False;

        Which[
            (* 1. Rational Function Optimisation (Goodwin-Kuprov 2016) *)
            variant === "NewtonRFO",
                eigSys = Eigensystem[hMin];
                evalsH = Developer`ToPackedArray[N[Re[eigSys[[1]]]]];
                evecsH = Developer`ToPackedArray[N[Re[eigSys[[2]]]]];
                minIdx = First[Ordering[evalsH, 1]];
                lamMin = evalsH[[minIdx]];
                If[TrueQ[lamMin < -1.*^-12], negCurvCount++];

                If[explicitLSQ,
                    rfoOut = SolveRFOSecular[evalsH, evecsH, gMin, sRFO];
                    p = rfoOut["Direction"];
                    If[TrueQ[Norm[p] > maxStepNorm], p = Developer`ToPackedArray[(maxStepNorm / Norm[p]) * p]];
                    If[!TrueQ[gMax . p > 1.*^-14 * gNorm * Max[Norm[p], 1.*^-30]], p = Developer`ToPackedArray[-gMin]];
                    lsOut = runNewtonLineSearch[evalCacheScaled, th, fVal, gMax, p, lsCfg, tolScale, nmState];
                    If[AssociationQ[lsOut] && TrueQ[lsOut["Accepted"]],
                        accepted = True;
                        thNew = lsOut["Theta"]; resNew = lsOut["Res"]; lastStep = Norm[thNew - th];
                        If[TrueQ[lsOut["StrongWolfeSatisfied"]], swCount++];
                        If[TrueQ[lsOut["ApproxWolfeSatisfied"]], approxWolfeCnt++];
                        If[KeyExistsQ[lsOut, "NMState"], nmState = lsOut["NMState"]];
                        rejectedSteps += Max[0, Lookup[lsOut, "Trials", 1] - 1]
                    ],
                    Do[
                        rfoOut = SolveRFOSecular[evalsH, evecsH, gMin, sRFO];
                        p = rfoOut["Direction"];
                        If[TrueQ[Norm[p] > maxStepNorm] && TrueQ[sRFO < 0.25 * sRFOMax],
                            sRFO = Min[sRFOMax, Max[4. * sRFO, sRFO * Sqrt[Norm[p] / maxStepNorm]]];
                            Continue[]
                        ];
                        If[TrueQ[Norm[p] <= 1.*^-15 * Max[1., Norm[th]]],
                            Break[]
                        ];
                        resNew = evalCacheScaled[Developer`ToPackedArray[th + p]];
                        If[resNew === $Failed,
                            rejectedSteps++;
                            sRFO = Min[sRFOMax, 4. * sRFO];
                            Continue[]
                        ];
                        phiNew = -resNew["Objective"];
                        ared = phiVal - phiNew;
                        predDec = -(gMin . p + 0.5 * (p . (hMin . p)));
                        If[!TrueQ[predDec > 0.], predDec = -(gMin . p)];
                        tolEff = Min[tolScale, 0.25 * 1.*^-4 * Max[0., predDec]];
                        rho = If[TrueQ[predDec > 1.*^-24], ared / predDec, If[TrueQ[ared > -tolEff], 1., -1.]];
                        If[TrueQ[ared >= 1.*^-4 * Max[0., predDec] - tolEff] && TrueQ[phiNew <= phiVal + tolScale],
                            accepted = True;
                            thNew = Developer`ToPackedArray[th + p];
                            lastStep = Norm[p];
                            If[TrueQ[rho >= 0.75], sRFO = Max[sRFOMin, 0.5 * sRFO]];
                            Break[],
                            rejectedSteps++;
                            sRFO = Min[sRFOMax, Max[2.5 * sRFO, gNorm / Max[1.*^-12, 0.5 * Norm[p]]]]
                        ],
                        {inIter, maxInner}
                    ];
                    If[!accepted && VectorQ[p, NumericQ] && TrueQ[gMax . p > 0.],
                        lsOut = runNewtonLineSearch[evalCacheScaled, th, fVal, gMax, p,
                            <|"Mode" -> "Armijo", "C1" -> 1.*^-4, "C2" -> 0.9, "MaxEvals" -> maxLsOpt, "ExtraOpts" -> <||>|>,
                            tolScale, nmState];
                        If[AssociationQ[lsOut] && TrueQ[lsOut["Accepted"]],
                            accepted = True;
                            thNew = lsOut["Theta"]; resNew = lsOut["Res"]; lastStep = Norm[thNew - th];
                            rejectedSteps += Max[0, Lookup[lsOut, "Trials", 1] - 1]
                        ]
                    ]
                ],

            (* 2. Adaptive Diagonal-Shift Levenberg-Marquardt Newton (Cholesky) *)
            variant === "NewtonShift",
                hScale = Max[1., Max[Abs[Diagonal[hMin]]]];
                Do[
                    chol = If[TrueQ[muShift == 0.],
                        Quiet[Check[CholeskyDecomposition[hMin], $Failed]],
                        Quiet[Check[CholeskyDecomposition[hMin + muShift * IdentityMatrix[nDim, WorkingPrecision -> MachinePrecision]], $Failed]]
                    ];
                    If[chol === $Failed,
                        negCurvCount++;
                        If[shiftTry <= 2,
                            muShift = Min[muMax, Max[1.*^-5 * hScale, 8. * muShift]];
                            Continue[],
                            lamMin = Min[Re[Eigenvalues[hMin]]];
                            muShift = Min[muMax, Max[muShift, -lamMin + 1.*^-4 * Max[1., Abs[lamMin]]]];
                            chol = Quiet[Check[CholeskyDecomposition[hMin + muShift * IdentityMatrix[nDim, WorkingPrecision -> MachinePrecision]], $Failed]];
                            If[chol === $Failed, muShift = Min[muMax, Max[1., 4. * muShift]]; Continue[]]
                        ]
                    ];
                    p = Developer`ToPackedArray[LinearSolve[chol, LinearSolve[Transpose[chol], -gMin]]];
                    If[!VectorQ[p, NumericQ] || !TrueQ[gMax . p > 0.],
                        rejectedSteps++;
                        muShift = Min[muMax, Max[1.*^-4 * hScale, 4. * muShift]];
                        Continue[]
                    ];
                    If[TrueQ[Norm[p] <= 1.*^-15 * Max[1., Norm[th]]], Break[]];
                    resNew = evalCacheScaled[Developer`ToPackedArray[th + p]];
                    If[resNew === $Failed,
                        rejectedSteps++;
                        muShift = Min[muMax, Max[1.*^-4 * hScale, 4. * muShift]];
                        Continue[]
                    ];
                    phiNew = -resNew["Objective"];
                    ared = phiVal - phiNew;
                    predDec = -0.5 * (gMin . p) + 0.5 * muShift * (p . p);
                    tolEff = Min[tolScale, 0.25 * 1.*^-4 * Max[0., predDec]];
                    rho = If[TrueQ[predDec > 1.*^-24], ared / predDec, -1.];
                    If[TrueQ[rho >= 1.*^-4] && TrueQ[ared > -tolEff],
                        accepted = True;
                        thNew = Developer`ToPackedArray[th + p];
                        lastStep = Norm[p];
                        Which[
                            TrueQ[rho >= 0.75], muShift = If[TrueQ[muShift <= 1.*^-8 * hScale], 0., Max[0., 0.25 * muShift]],
                            TrueQ[rho < 0.25], muShift = Min[muMax, Max[1.*^-6 * hScale, 2. * muShift]]
                        ];
                        Break[],
                        rejectedSteps++;
                        muShift = Min[muMax, Max[1.*^-4 * hScale, 4. * muShift]]
                    ],
                    {shiftTry, maxInner}
                ],

            (* 3. Pivoted Gill-Murray-Wright Modified Cholesky + Line Search *)
            variant === "NewtonModifiedCholesky",
                mcOut = ModifiedCholeskySolve[hMin, gMin, deltaMC];
                p = mcOut["Direction"];
                If[TrueQ[mcOut["Modified"]], negCurvCount++];
                If[!VectorQ[p, NumericQ] || !TrueQ[gMax . p > 1.*^-14 * gNorm * Max[Norm[p], 1.*^-30]],
                    p = Developer`ToPackedArray[-gMin]
                ];
                If[TrueQ[Norm[p] > maxStepNorm], p = Developer`ToPackedArray[(maxStepNorm / Norm[p]) * p]];
                lsOut = runNewtonLineSearch[evalCacheScaled, th, fVal, gMax, p, lsCfg, tolScale, nmState];
                If[AssociationQ[lsOut] && TrueQ[lsOut["Accepted"]],
                    accepted = True;
                    thNew = lsOut["Theta"]; resNew = lsOut["Res"]; lastStep = Norm[thNew - th];
                    If[TrueQ[lsOut["StrongWolfeSatisfied"]], swCount++];
                    If[TrueQ[lsOut["ApproxWolfeSatisfied"]], approxWolfeCnt++];
                    If[KeyExistsQ[lsOut, "NMState"], nmState = lsOut["NMState"]];
                    rejectedSteps += Max[0, Lookup[lsOut, "Trials", 1] - 1]
                ],

            (* 4. Spectral Absolute-Eigenvalue Saddle-Free Newton (Dauphin / Nocedal) + Line Search *)
            variant === "NewtonSpectral",
                eigSys = Eigensystem[hMin];
                evalsH = Developer`ToPackedArray[N[Re[eigSys[[1]]]]];
                evecsH = Developer`ToPackedArray[N[Re[eigSys[[2]]]]];
                gTilde = Developer`ToPackedArray[evecsH . gMin];
                minIdx = First[Ordering[evalsH, 1]];
                lamMin = evalsH[[minIdx]];
                maxAbsH = Max[Abs[evalsH]];
                If[TrueQ[lamMin < -1.*^-12], negCurvCount++];
                With[{eFloor = Max[eigFloorScale, 1.*^-10 * Max[1., maxAbsH]]},
                    lamAbs = Developer`ToPackedArray[Map[Max[Abs[#], eFloor] &, evalsH]]
                ];
                p = Developer`ToPackedArray[-((gTilde / lamAbs) . evecsH)];
                (* Hard-case saddle escape when gMin is nearly orthogonal to negative-curvature eigendirection *)
                If[TrueQ[lamMin < -1.*^-8] && TrueQ[gNorm < 1.*^-2] && TrueQ[Abs[gTilde[[minIdx]]] < 1.*^-2 * Max[gNorm, 1.*^-15]],
                    qMin = evecsH[[minIdx]];
                    If[TrueQ[gMin . qMin > 0.], qMin = -qMin];
                    tauEsc = Min[0.5, Max[0.05, Sqrt[Max[gNorm, 1.*^-6]] / Max[1., Abs[lamMin]]]];
                    p = Developer`ToPackedArray[p + tauEsc * qMin]
                ];
                If[TrueQ[Norm[p] > maxStepNorm], p = Developer`ToPackedArray[(maxStepNorm / Norm[p]) * p]];
                If[!TrueQ[gMax . p > 0.], p = Developer`ToPackedArray[-gMin]];
                lsOut = runNewtonLineSearch[evalCacheScaled, th, fVal, gMax, p, lsCfg, tolScale, nmState];
                If[AssociationQ[lsOut] && TrueQ[lsOut["Accepted"]],
                    accepted = True;
                    thNew = lsOut["Theta"]; resNew = lsOut["Res"]; lastStep = Norm[thNew - th];
                    If[TrueQ[lsOut["StrongWolfeSatisfied"]], swCount++];
                    If[TrueQ[lsOut["ApproxWolfeSatisfied"]], approxWolfeCnt++];
                    If[KeyExistsQ[lsOut, "NMState"], nmState = lsOut["NMState"]];
                    rejectedSteps += Max[0, Lookup[lsOut, "Trials", 1] - 1]
                ]
        ];

        If[!accepted,
            status = If[gNorm <= Sqrt[gradTol], "GradientToleranceReached", "LineSearchTerminated"];
            Break[]
        ];

        acceptedSteps++;
        sVec = Developer`ToPackedArray[thNew - th];
        fOld = fVal;
        th = thNew;
        {fVal, fidVal, penVal, gMax, extra, evalState} =
            {resNew["Objective"], resNew["Fidelity"], resNew["Penalty"], resNew["Gradient"], resNew["Extra"], resNew["EvalState"]};
        phiVal = -fVal;
        gMin = Developer`ToPackedArray[-gMax];
        gNorm = Norm[gMin];
        bestFid = Max[bestFid, fidVal];
        nIter++;

        If[TrueQ[verbose] && (Mod[nIter, 5] === 0 || nIter === 1),
            Print["Iter ", nIter, " (", variant, "): Fidelity=", NumberForm[fidVal, 6],
                " Objective=", NumberForm[fVal, 6], " |step|=", NumberForm[lastStep, 4],
                " |grad|=", ScientificForm[gNorm, 3]]
        ];
        If[TrueQ[recHist],
            AppendTo[objHist, fVal];
            AppendTo[fidHist, fidVal];
            AppendTo[gHist, gNorm]
        ];

        With[{critFidNew = evalCritFn[fidVal, extra, fidCrit]},
            stopNow = Which[
                TrueQ[Abs[fVal - fOld] <= funTol * Max[1., Abs[fOld]]] && TrueQ[gNorm <= Sqrt[gradTol]],
                    If[NumericQ[targetFid] && TrueQ[critFidNew >= targetFid], "TargetFidelityReached", "FunctionToleranceReached"],
                TrueQ[Norm[sVec] <= stepTol * Max[1., Norm[th]]],
                    If[NumericQ[targetFid] && TrueQ[critFidNew >= targetFid], "TargetFidelityReached", "StepToleranceReached"],
                NumericQ[targetFid] && TrueQ[critFidNew >= targetFid],
                    "TargetFidelityReached",
                TrueQ[gNorm <= gradTol],
                    "GradientToleranceReached",
                AbsoluteTime[] - tStart >= tLimit,
                    "TimeConstraintReached",
                nIter >= maxIter,
                    "MaxIterationsReached",
                True,
                    None
            ]
        ];
        {lastMonTime, lastMonIter} = emitMonFn[
            monCfg, nIter, evals, fVal, fidVal, bestFid, gNorm, lastStep,
            AbsoluteTime[] - tStart, If[stopNow === None, "Running", stopNow],
            extra, stageTag, lastMonTime, lastMonIter
        ];
        If[stopNow =!= None, status = stopNow; Break[]],
        {maxIter}
    ];
    If[lastMonIter < nIter || (lastMonIter < 0 && nIter === 0),
        {lastMonTime, lastMonIter} = emitMonFn[monCfg, nIter, evals, fVal, fidVal, bestFid, gNorm, lastStep, AbsoluteTime[] - tStart, status, extra, stageTag, lastMonTime, lastMonIter]
    ];

    With[{elapsed = AbsoluteTime[] - tStart},
        <|
            "ParametersVec" -> Developer`ToPackedArray[th * invS],
            "Objective" -> fVal,
            "Fidelity" -> fidVal,
            "Penalty" -> penVal,
            "GradientNorm" -> gNorm,
            "ScaledGradient" -> gMax,
            "Optimizer" -> variant,
            "Iterations" -> nIter,
            "Evaluations" -> evals,
            "ObjectiveEvaluations" -> objEvals,
            "GradientEvaluations" -> gradEvals,
            "HVPEvaluations" -> hvpEvals,
            "HessianEvaluations" -> hessEvals,
            "CGIterations" -> 0,
            "NegativeCurvatureCount" -> negCurvCount,
            "AcceptedSteps" -> acceptedSteps,
            "RejectedSteps" -> rejectedSteps,
            "SkippedCurvaturePairs" -> 0,
            "LineSearchMethod" -> lsCfg["Mode"],
            "StrongWolfeSatisfiedCount" -> swCount,
            "ApproxWolfeCount" -> approxWolfeCnt,
            "LastStepLength" -> lastStep,
            "Status" -> status,
            "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
            "FidelityHistory" -> Developer`ToPackedArray[fidHist],
            "GradientNormHistory" -> Developer`ToPackedArray[gHist],
            "ElapsedTime" -> elapsed,
            "DynamicsTime" -> dynTime,
            "OptimizerTime" -> Max[0., elapsed - dynTime],
            "Extra" -> extra
        |>
    ]
];

(* Matrix-Free Steihaug-Toint Truncated Newton-CG inside a Trust Region with optional Morales-Nocedal L-BFGS Preconditioner *)
RunTrustRegionNewtonCG[ctx_Association] := Module[
    {soCtx, evalWithCacheOracle, evalObjOracle, evalHVPOracle,
     x0Vec = ctx["X0Vec"], optsList = ctx["OptsList"], subOpts = Lookup[ctx, "SubOpts", <||>],
     scaleVec, invS, comp, hasOracles, nDim,
     maxIter, targetFid, fidCrit, gradTol, funTol, stepTol, tLimit, recHist, verbose,
     monOpt, monIntOpt, monCfg, stageTag, lastMonTime = 0., lastMonIter = -1,
     parseMonFn, emitMonFn, evalCritFn,
     precMode, precMem, forcingExp, maxCG, etaAccept, deltaMax, deltaMin, useObjTrial,
     tStart, dynTime = 0., evals = 0, objEvals = 0, gradEvals = 0, hvpEvals = 0, cgIters = 0,
     acceptedSteps = 0, rejectedSteps = 0, negCurvCount = 0, skippedCurv = 0,
     sList = {}, yList = {}, rhoList = {}, initPrecHist,
     evalCacheScaled, evalObjScaled, evalHVPMinScaled, applyPrec, solveSubproblem,
     th, res0, fVal, phiVal, fidVal, penVal, gMax, gMin, gNorm, extra, bestFid, evalState,
     delta, objHist, fidHist, gHist, nIter = 0, status = "MaxIterationsReached",
     lastStep = 0., fOld, stopNow, tolScale, tolEff, trRetry,
     subOut, p, predDec, pMNorm, hitBoundary, thTrial, resTrial, objResTrial, phiTrial,
     ared, rho, accepted, sVec, yVec, sy, gradRes},

    soCtx = buildSecondOrderContext[ctx, "TrustRegionNewtonCG", False];
    If[soCtx === $Failed, Return[$Failed]];
    scaleVec = soCtx["ScaleVec"];
    invS = soCtx["InvScale"];
    nDim = soCtx["NumDim"];
    comp = soCtx["Comp"];
    hasOracles = soCtx["HasCustomOracles"];

    evalWithCacheOracle = If[hasOracles, ctx["EvalWithCacheFn"], $OCSolverRegistry["EvaluateWithCache"]];
    evalObjOracle = Lookup[ctx, "EvalObjFn", None];
    evalHVPOracle = If[hasOracles, ctx["EvalHVPFn"], $OCSolverRegistry["EvaluateHVP"]];

    {maxIter, targetFid, gradTol, funTol, stepTol, tLimit, recHist, verbose, monOpt, monIntOpt} =
        {MaxIterations, TargetFidelity, GradientTolerance, FunctionTolerance, StepTolerance,
         TimeConstraint, RecordHistory, Verbose, StepMonitor, MonitorInterval} /. optsList /. Options[GRAPEOptimize];
    fidCrit = Lookup[FilterRules[optsList, {$OCSolverRegistry["FidelityCriterionKey"]}], $OCSolverRegistry["FidelityCriterionKey"], Automatic];
    stageTag = Lookup[FilterRules[optsList, {$OCSolverRegistry["MonitorStageKey"]}], $OCSolverRegistry["MonitorStageKey"], "Optimization"];
    parseMonFn = $OCSolverRegistry["ParseMonitorConfig"];
    emitMonFn = $OCSolverRegistry["EmitMonitorIfDue"];
    evalCritFn = $OCSolverRegistry["EvalCriterionFidelity"];
    monCfg = parseMonFn[monOpt, monIntOpt];
    If[monCfg === $Failed, Return[$Failed]];

    precMode = Lookup[subOpts, "Preconditioner", "None"];
    If[!MemberQ[{"None", "LBFGS"}, precMode],
        Message[OptimalControlProblem::badopt, precMode]; Return[$Failed]
    ];
    precMem = Max[1, Round[Lookup[subOpts, "PreconditionerMemory", 10]]];
    forcingExp = Clip[N[Lookup[subOpts, "ForcingExponent", 0.5]], {0.1, 1.0}];
    maxCG = With[{rawMax = Lookup[subOpts, "MaxCGIterations", Automatic]},
        If[rawMax === Automatic, Min[nDim, 50], Max[1, Round[rawMax]]]
    ];
    etaAccept = Clip[N[Lookup[subOpts, "EtaAccept", 0.01]], {1.*^-6, 0.24}];
    deltaMax = Max[1., N[Lookup[subOpts, "MaxRadius", 100.]]];
    deltaMin = Max[1.*^-15, N[Lookup[subOpts, "MinRadius", 1.*^-14]]];
    useObjTrial = TrueQ[Lookup[subOpts, "UseObjectiveOnlyTrial", False]] && (evalObjOracle =!= None);

    initPrecHist = Lookup[subOpts, "InitialPreconditionerHistory", None];
    If[precMode === "LBFGS" && ListQ[initPrecHist] && Length[initPrecHist] === 3,
        {sList, yList, rhoList} = initPrecHist;
        If[Length[sList] > precMem,
            sList = Take[sList, -precMem];
            yList = Take[yList, -precMem];
            rhoList = Take[rhoList, -precMem]
        ]
    ];

    evalCacheScaled[thIn_] := Module[{t0Dyn = AbsoluteTime[], xIn, r, gFlat},
        xIn = Developer`ToPackedArray[thIn * invS];
        r = If[hasOracles, evalWithCacheOracle[xIn], evalWithCacheOracle[comp, xIn]];
        dynTime += AbsoluteTime[] - t0Dyn;
        evals++; objEvals++;
        If[r === $Failed || !AssociationQ[r] || !NumericQ[r["Objective"]] || !KeyExistsQ[r, "Gradient"],
            Return[$Failed]
        ];
        gradEvals++;
        gFlat = Developer`ToPackedArray[Flatten[r["Gradient"]] * invS];
        <|
            "Objective" -> N[Re[r["Objective"]]],
            "Fidelity" -> N[Re[r["Fidelity"]]],
            "Penalty" -> N[Re[Lookup[r, "Penalty", 0.]]],
            "Gradient" -> gFlat,
            "Extra" -> Lookup[r, "Extra", <||>],
            "EvalState" -> Lookup[r, "EvalState", None]
        |>
    ];

    evalObjScaled[thIn_] := Module[{t0Dyn = AbsoluteTime[], xIn, r},
        xIn = Developer`ToPackedArray[thIn * invS];
        r = evalObjOracle[xIn];
        dynTime += AbsoluteTime[] - t0Dyn;
        evals++; objEvals++;
        r
    ];

    evalHVPMinScaled[evalSt_, vTh_] := Module[{t0Dyn = AbsoluteTime[], vX, hX},
        vX = Developer`ToPackedArray[vTh * invS];
        hX = If[hasOracles,
            evalHVPOracle[evalSt, vX],
            evalHVPOracle[comp, evalSt, vX]
        ];
        dynTime += AbsoluteTime[] - t0Dyn;
        hvpEvals++;
        If[!ArrayQ[hX, _, NumericQ], Return[$Failed]];
        Developer`ToPackedArray[-(Flatten[hX] * invS)]
    ];

    applyPrec[vec_] := If[precMode === "LBFGS" && Length[sList] > 0,
        lbfgsTwoLoop[sList, yList, rhoList, vec],
        Developer`ToPackedArray[N[vec]]
    ];

    solveSubproblem[evalSt_, deltaCur_] := Module[
        {z, u, r, y, d, v, ry0, r0NormM, cgTol, ry, qVal = 0., j, hd, dHd, aQuad, bQuad, cQuad,
         disc, tau, pOut, predOut, mNormOut, alpha, zNext, uNext, zNextMSq, rNext, yNext, ryNext, beta, rd},
        z = ConstantArray[0., nDim];
        u = ConstantArray[0., nDim];
        r = gMin;
        y = applyPrec[r];
        ry = r . y;
        If[!TrueQ[ry > 0.],
            sList = {}; yList = {}; rhoList = {};
            y = r; ry = r . r
        ];
        d = Developer`ToPackedArray[-y];
        v = Developer`ToPackedArray[-r];
        ry0 = ry;
        r0NormM = Sqrt[Max[ry0, 0.]];
        cgTol = Min[0.5, gNorm^forcingExp] * r0NormM;

        Do[
            hd = evalHVPMinScaled[evalSt, d];
            If[hd === $Failed, Return[$Failed, Module]];
            cgIters++;
            dHd = d . hd;
            aQuad = d . v;

            (* Negative or zero curvature along conjugate direction: jump to trust boundary *)
            If[TrueQ[dHd <= 1.*^-15 * Max[aQuad, d . d * 1.*^-15]],
                negCurvCount++;
                bQuad = z . v;
                cQuad = (z . u) - deltaCur * deltaCur;
                disc = Max[0., bQuad * bQuad - aQuad * cQuad];
                tau = If[TrueQ[bQuad >= 0.],
                    (-cQuad) / Max[bQuad + Sqrt[disc], 1.*^-30],
                    (-bQuad + Sqrt[disc]) / Max[aQuad, 1.*^-30]
                ];
                pOut = Developer`ToPackedArray[z + tau * d];
                rd = r . d;
                predOut = -(qVal + tau * rd + 0.5 * tau * tau * dHd);
                If[!TrueQ[predOut > 0.], predOut = -(gMin . pOut)];
                Return[<|"Step" -> pOut, "Pred" -> predOut, "MNorm" -> deltaCur, "Boundary" -> True|>, Module]
            ];

            alpha = ry / dHd;
            zNext = Developer`ToPackedArray[z + alpha * d];
            uNext = Developer`ToPackedArray[u + alpha * v];
            zNextMSq = zNext . uNext;

            (* Exceeded trust-region ellipsoid: intersect boundary along d *)
            If[TrueQ[zNextMSq >= deltaCur * deltaCur],
                bQuad = z . v;
                cQuad = (z . u) - deltaCur * deltaCur;
                disc = Max[0., bQuad * bQuad - aQuad * cQuad];
                tau = If[TrueQ[bQuad >= 0.],
                    (-cQuad) / Max[bQuad + Sqrt[disc], 1.*^-30],
                    (-bQuad + Sqrt[disc]) / Max[aQuad, 1.*^-30]
                ];
                pOut = Developer`ToPackedArray[z + tau * d];
                rd = r . d;
                predOut = -(qVal + tau * rd + 0.5 * tau * tau * dHd);
                If[!TrueQ[predOut > 0.], predOut = -(gMin . pOut)];
                Return[<|"Step" -> pOut, "Pred" -> predOut, "MNorm" -> deltaCur, "Boundary" -> True|>, Module]
            ];

            rd = r . d;
            qVal = qVal + alpha * rd + 0.5 * alpha * alpha * dHd;
            z = zNext;
            u = uNext;
            rNext = Developer`ToPackedArray[r + alpha * hd];
            yNext = applyPrec[rNext];
            ryNext = rNext . yNext;
            If[!TrueQ[ryNext > 0.],
                Return[<|"Step" -> z, "Pred" -> -qVal, "MNorm" -> Sqrt[Max[0., z . u]], "Boundary" -> False|>, Module]
            ];
            If[TrueQ[Sqrt[ryNext] <= cgTol],
                Return[<|"Step" -> z, "Pred" -> -qVal, "MNorm" -> Sqrt[Max[0., z . u]], "Boundary" -> False|>, Module]
            ];
            beta = ryNext / ry;
            d = Developer`ToPackedArray[-yNext + beta * d];
            v = Developer`ToPackedArray[-rNext + beta * v];
            r = rNext;
            y = yNext;
            ry = ryNext,
            {j, maxCG}
        ];

        mNormOut = Sqrt[Max[0., z . u]];
        If[TrueQ[mNormOut == 0.],
            (* Fallback to Cauchy point along -y if CG exited before moving *)
            z = Developer`ToPackedArray[-Clip[deltaCur / Max[r0NormM, 1.*^-30], {0., 1.}] * y];
            predOut = -(gMin . z);
            Return[<|"Step" -> z, "Pred" -> predOut, "MNorm" -> Norm[z], "Boundary" -> True|>]
        ];
        <|"Step" -> z, "Pred" -> -qVal, "MNorm" -> mNormOut, "Boundary" -> False|>
    ];

    tStart = AbsoluteTime[];
    th = Developer`ToPackedArray[N[x0Vec * scaleVec]];
    res0 = evalCacheScaled[th];
    If[res0 === $Failed, Return[$Failed]];
    {fVal, fidVal, penVal, gMax, extra, evalState} =
        {res0["Objective"], res0["Fidelity"], res0["Penalty"], res0["Gradient"], res0["Extra"], res0["EvalState"]};
    phiVal = -fVal;
    gMin = Developer`ToPackedArray[-gMax];
    gNorm = Norm[gMin];
    bestFid = fidVal;

    delta = With[{rawD = Lookup[subOpts, "InitialRadius", Automatic]},
        If[NumericQ[rawD], Clip[N[rawD], {1.*^-6, deltaMax}], Clip[gNorm, {0.5, 10.}]]
    ];

    objHist = If[recHist, {fVal}, {}];
    fidHist = If[recHist, {fidVal}, {}];
    gHist = If[recHist, {gNorm}, {}];

    With[{critFid0 = evalCritFn[fidVal, extra, fidCrit]},
        status = Which[
            NumericQ[targetFid] && TrueQ[critFid0 >= targetFid], "TargetFidelityReached",
            TrueQ[gNorm <= gradTol], "GradientToleranceReached",
            True, "Running"
        ]
    ];
    {lastMonTime, lastMonIter} = emitMonFn[monCfg, 0, evals, fVal, fidVal, bestFid, gNorm, 0., 0., status, extra, stageTag, lastMonTime, lastMonIter];
    If[status =!= "Running",
        With[{elapsed = AbsoluteTime[] - tStart},
            Return[<|
                "ParametersVec" -> Developer`ToPackedArray[th * invS],
                "Objective" -> fVal, "Fidelity" -> fidVal, "Penalty" -> penVal,
                "GradientNorm" -> gNorm, "ScaledGradient" -> gMax, "Optimizer" -> "TrustRegionNewtonCG",
                "Iterations" -> 0, "Evaluations" -> evals, "ObjectiveEvaluations" -> objEvals,
                "GradientEvaluations" -> gradEvals, "HVPEvaluations" -> hvpEvals, "HessianEvaluations" -> 0,
                "CGIterations" -> cgIters, "NegativeCurvatureCount" -> negCurvCount,
                "TrustRegionFinalRadius" -> delta,
                "AcceptedSteps" -> 0, "RejectedSteps" -> 0, "SkippedCurvaturePairs" -> skippedCurv,
                "LineSearchMethod" -> "TrustRegion", "StrongWolfeSatisfiedCount" -> 0, "ApproxWolfeCount" -> 0,
                "LastStepLength" -> 0., "Status" -> status,
                "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
                "FidelityHistory" -> Developer`ToPackedArray[fidHist],
                "GradientNormHistory" -> Developer`ToPackedArray[gHist],
                "ElapsedTime" -> elapsed, "DynamicsTime" -> dynTime, "OptimizerTime" -> Max[0., elapsed - dynTime],
                "Extra" -> extra
            |>]
        ]
    ];
    status = "MaxIterationsReached";

    Do[
        tolScale = 1.*^-14 * Max[1., Abs[fVal]];
        accepted = False;
        Do[
            subOut = solveSubproblem[evalState, delta];
            If[subOut === $Failed, Break[]];
            p = subOut["Step"];
            predDec = subOut["Pred"];
            pMNorm = subOut["MNorm"];
            hitBoundary = subOut["Boundary"];
            If[TrueQ[Norm[p] <= 1.*^-15 * Max[1., Norm[th]]],
                Break[]
            ];
            thTrial = Developer`ToPackedArray[th + p];
            If[useObjTrial,
                objResTrial = evalObjScaled[thTrial];
                If[objResTrial === $Failed || !NumericQ[objResTrial["Objective"]],
                    rejectedSteps++;
                    delta = 0.25 * Max[Min[pMNorm, delta], 1.*^-30];
                    If[TrueQ[delta < deltaMin * Max[1., Norm[th]]], Break[]];
                    Continue[]
                ];
                phiTrial = -objResTrial["Objective"],
                resTrial = evalCacheScaled[thTrial];
                If[resTrial === $Failed,
                    rejectedSteps++;
                    delta = 0.25 * Max[Min[pMNorm, delta], 1.*^-30];
                    If[TrueQ[delta < deltaMin * Max[1., Norm[th]]], Break[]];
                    Continue[]
                ];
                phiTrial = -resTrial["Objective"]
            ];

            ared = phiVal - phiTrial;
            tolEff = Min[tolScale, 0.25 * etaAccept * Max[0., predDec]];
            rho = If[TrueQ[predDec > 1.*^-24], ared / predDec, -1.];

            If[TrueQ[rho < 0.25],
                delta = 0.25 * Max[Min[pMNorm, delta], 1.*^-30],
                If[TrueQ[rho > 0.75 && hitBoundary],
                    delta = Min[deltaMax, 2.0 * delta]
                ]
            ];

            If[TrueQ[rho > etaAccept] && TrueQ[ared > -tolEff],
                gradRes = If[useObjTrial, evalCacheScaled[thTrial], resTrial];
                If[gradRes === $Failed, Return[$Failed, Module]];
                accepted = True;
                Break[],
                rejectedSteps++;
                If[TrueQ[delta < deltaMin * Max[1., Norm[th]]], Break[]]
            ],
            {trRetry, 18}
        ];

        If[!accepted,
            status = If[gNorm <= Sqrt[gradTol], "GradientToleranceReached", "StepToleranceReached"];
            Break[]
        ];

        acceptedSteps++;
        sVec = Developer`ToPackedArray[thTrial - th];
        lastStep = Norm[sVec];
        fOld = fVal;
        th = thTrial;
        {fVal, fidVal, penVal, gMax, extra, evalState} =
            {gradRes["Objective"], gradRes["Fidelity"], gradRes["Penalty"], gradRes["Gradient"], gradRes["Extra"], gradRes["EvalState"]};
        phiVal = -fVal;
        If[precMode === "LBFGS",
            yVec = Developer`ToPackedArray[(-gMax) - gMin];
            sy = sVec . yVec;
            If[TrueQ[sy > 1.*^-12 * Norm[sVec] * Norm[yVec]],
                If[Length[sList] >= precMem,
                    sList = Rest[sList]; yList = Rest[yList]; rhoList = Rest[rhoList]
                ];
                AppendTo[sList, sVec]; AppendTo[yList, yVec]; AppendTo[rhoList, 1. / sy],
                skippedCurv++
            ]
        ];
        gMin = Developer`ToPackedArray[-gMax];
        gNorm = Norm[gMin];
        bestFid = Max[bestFid, fidVal];
        nIter++;

        If[TrueQ[verbose] && (Mod[nIter, 5] === 0 || nIter === 1),
            Print["Iter ", nIter, " (TR-NewtonCG): Fidelity=", NumberForm[fidVal, 6],
                " Objective=", NumberForm[fVal, 6], " delta=", ScientificForm[delta, 2],
                " cg=", cgIters, " |grad|=", ScientificForm[gNorm, 3]]
        ];
        If[TrueQ[recHist],
            AppendTo[objHist, fVal];
            AppendTo[fidHist, fidVal];
            AppendTo[gHist, gNorm]
        ];

        With[{critFidNew = evalCritFn[fidVal, extra, fidCrit]},
            stopNow = Which[
                TrueQ[Abs[fVal - fOld] <= funTol * Max[1., Abs[fOld]]] && TrueQ[gNorm <= Sqrt[gradTol]],
                    If[NumericQ[targetFid] && TrueQ[critFidNew >= targetFid], "TargetFidelityReached", "FunctionToleranceReached"],
                TrueQ[Norm[sVec] <= stepTol * Max[1., Norm[th]]],
                    If[NumericQ[targetFid] && TrueQ[critFidNew >= targetFid], "TargetFidelityReached", "StepToleranceReached"],
                NumericQ[targetFid] && TrueQ[critFidNew >= targetFid],
                    "TargetFidelityReached",
                TrueQ[gNorm <= gradTol],
                    "GradientToleranceReached",
                AbsoluteTime[] - tStart >= tLimit,
                    "TimeConstraintReached",
                nIter >= maxIter,
                    "MaxIterationsReached",
                True,
                    None
            ]
        ];
        {lastMonTime, lastMonIter} = emitMonFn[
            monCfg, nIter, evals, fVal, fidVal, bestFid, gNorm, lastStep,
            AbsoluteTime[] - tStart, If[stopNow === None, "Running", stopNow],
            extra, stageTag, lastMonTime, lastMonIter
        ];
        If[stopNow =!= None, status = stopNow; Break[]],
        {maxIter}
    ];
    If[lastMonIter < nIter || (lastMonIter < 0 && nIter === 0),
        {lastMonTime, lastMonIter} = emitMonFn[monCfg, nIter, evals, fVal, fidVal, bestFid, gNorm, lastStep, AbsoluteTime[] - tStart, status, extra, stageTag, lastMonTime, lastMonIter]
    ];

    With[{elapsed = AbsoluteTime[] - tStart},
        <|
            "ParametersVec" -> Developer`ToPackedArray[th * invS],
            "Objective" -> fVal,
            "Fidelity" -> fidVal,
            "Penalty" -> penVal,
            "GradientNorm" -> gNorm,
            "ScaledGradient" -> gMax,
            "Optimizer" -> "TrustRegionNewtonCG",
            "Iterations" -> nIter,
            "Evaluations" -> evals,
            "ObjectiveEvaluations" -> objEvals,
            "GradientEvaluations" -> gradEvals,
            "HVPEvaluations" -> hvpEvals,
            "HessianEvaluations" -> 0,
            "CGIterations" -> cgIters,
            "NegativeCurvatureCount" -> negCurvCount,
            "TrustRegionFinalRadius" -> delta,
            "AcceptedSteps" -> acceptedSteps,
            "RejectedSteps" -> rejectedSteps,
            "SkippedCurvaturePairs" -> skippedCurv,
            "LineSearchMethod" -> If[precMode === "LBFGS", "TrustRegion_LBFGSPrec", "TrustRegion"],
            "StrongWolfeSatisfiedCount" -> 0,
            "ApproxWolfeCount" -> 0,
            "LastStepLength" -> lastStep,
            "Status" -> status,
            "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
            "FidelityHistory" -> Developer`ToPackedArray[fidHist],
            "GradientNormHistory" -> Developer`ToPackedArray[gHist],
            "ElapsedTime" -> elapsed,
            "DynamicsTime" -> dynTime,
            "OptimizerTime" -> Max[0., elapsed - dynTime],
            "Extra" -> extra
        |>
    ]
];

(* Shared Krylov-Lanczos solver for TrustRegionGLTR and AdaptiveCubicRegularization:
   Builds reorthogonalised Lanczos basis (Q_j, T_j) with zero-HVP reuse on subproblem rejections. *)
runKrylovSubspaceSolver[ctx_Association, mode_String] := Module[
    {soCtx, evalWithCacheOracle, evalHVPOracle,
     x0Vec = ctx["X0Vec"], optsList = ctx["OptsList"], subOpts = Lookup[ctx, "SubOpts", <||>],
     scaleVec, invS, comp, hasOracles, nDim,
     maxIter, targetFid, fidCrit, gradTol, funTol, stepTol, tLimit, recHist, verbose,
     monOpt, monIntOpt, monCfg, stageTag, lastMonTime = 0., lastMonIter = -1,
     parseMonFn, emitMonFn, evalCritFn,
     maxKrylov, forcingExp, eta1, eta2, deltaMax, deltaMin, sigmaMin, sigmaMax,
     tStart, dynTime = 0., evals = 0, objEvals = 0, gradEvals = 0, hvpEvals = 0, cgIters = 0,
     acceptedSteps = 0, rejectedSteps = 0, negCurvCount = 0,
     evalCacheScaled, evalHVPMinScaled,
     th, res0, fVal, phiVal, fidVal, penVal, gMax, gMin, gNorm, extra, bestFid, evalState,
     delta, sigma, objHist, fidHist, gHist, nIter = 0, status = "MaxIterationsReached",
     lastStep = 0., fOld, stopNow, tolScale, tolEff,
     qList, diagT, subDiagT, betaPrev, qCur, qPrev, hv, alphaJ, wVec, qMat, coeffs, betaNext,
     jDim, tMat, gSmall, subSol, hVec, resNorm, krylovTol, p, predDec, hitBoundary, pNorm,
     thTrial, resTrial, phiTrial, ared, rho, accepted, sVec, retry},

    soCtx = buildSecondOrderContext[ctx, mode, False];
    If[soCtx === $Failed, Return[$Failed]];
    scaleVec = soCtx["ScaleVec"];
    invS = soCtx["InvScale"];
    nDim = soCtx["NumDim"];
    comp = soCtx["Comp"];
    hasOracles = soCtx["HasCustomOracles"];

    evalWithCacheOracle = If[hasOracles, ctx["EvalWithCacheFn"], $OCSolverRegistry["EvaluateWithCache"]];
    evalHVPOracle = If[hasOracles, ctx["EvalHVPFn"], $OCSolverRegistry["EvaluateHVP"]];

    {maxIter, targetFid, gradTol, funTol, stepTol, tLimit, recHist, verbose, monOpt, monIntOpt} =
        {MaxIterations, TargetFidelity, GradientTolerance, FunctionTolerance, StepTolerance,
         TimeConstraint, RecordHistory, Verbose, StepMonitor, MonitorInterval} /. optsList /. Options[GRAPEOptimize];
    fidCrit = Lookup[FilterRules[optsList, {$OCSolverRegistry["FidelityCriterionKey"]}], $OCSolverRegistry["FidelityCriterionKey"], Automatic];
    stageTag = Lookup[FilterRules[optsList, {$OCSolverRegistry["MonitorStageKey"]}], $OCSolverRegistry["MonitorStageKey"], "Optimization"];
    parseMonFn = $OCSolverRegistry["ParseMonitorConfig"];
    emitMonFn = $OCSolverRegistry["EmitMonitorIfDue"];
    evalCritFn = $OCSolverRegistry["EvalCriterionFidelity"];
    monCfg = parseMonFn[monOpt, monIntOpt];
    If[monCfg === $Failed, Return[$Failed]];

    maxKrylov = With[{rawK = Lookup[subOpts, "MaxKrylovDim", Automatic]},
        If[rawK === Automatic, Min[nDim, 40], Max[1, Min[nDim, Round[rawK]]]]
    ];
    forcingExp = Clip[N[Lookup[subOpts, "ForcingExponent", 0.5]], {0.1, 1.0}];
    eta1 = Clip[N[Lookup[subOpts, "Eta1", 0.01]], {1.*^-6, 0.24}];
    eta2 = Clip[N[Lookup[subOpts, "Eta2", 0.75]], {eta1 + 0.05, 0.99}];
    deltaMax = Max[1., N[Lookup[subOpts, "MaxRadius", 100.]]];
    deltaMin = Max[1.*^-15, N[Lookup[subOpts, "MinRadius", 1.*^-14]]];
    sigmaMin = Max[1.*^-12, N[Lookup[subOpts, "MinSigma", 1.*^-8]]];
    sigmaMax = Min[1.*^14, N[Lookup[subOpts, "MaxSigma", 1.*^12]]];

    evalCacheScaled[thIn_] := Module[{t0Dyn = AbsoluteTime[], xIn, r, gFlat},
        xIn = Developer`ToPackedArray[thIn * invS];
        r = If[hasOracles, evalWithCacheOracle[xIn], evalWithCacheOracle[comp, xIn]];
        dynTime += AbsoluteTime[] - t0Dyn;
        evals++; objEvals++;
        If[r === $Failed || !AssociationQ[r] || !NumericQ[r["Objective"]] || !KeyExistsQ[r, "Gradient"],
            Return[$Failed]
        ];
        gradEvals++;
        gFlat = Developer`ToPackedArray[Flatten[r["Gradient"]] * invS];
        <|
            "Objective" -> N[Re[r["Objective"]]],
            "Fidelity" -> N[Re[r["Fidelity"]]],
            "Penalty" -> N[Re[Lookup[r, "Penalty", 0.]]],
            "Gradient" -> gFlat,
            "Extra" -> Lookup[r, "Extra", <||>],
            "EvalState" -> Lookup[r, "EvalState", None]
        |>
    ];

    evalHVPMinScaled[evalSt_, vTh_] := Module[{t0Dyn = AbsoluteTime[], vX, hX},
        vX = Developer`ToPackedArray[vTh * invS];
        hX = If[hasOracles,
            evalHVPOracle[evalSt, vX],
            evalHVPOracle[comp, evalSt, vX]
        ];
        dynTime += AbsoluteTime[] - t0Dyn;
        hvpEvals++;
        If[!ArrayQ[hX, _, NumericQ], Return[$Failed]];
        Developer`ToPackedArray[-(Flatten[hX] * invS)]
    ];

    tStart = AbsoluteTime[];
    th = Developer`ToPackedArray[N[x0Vec * scaleVec]];
    res0 = evalCacheScaled[th];
    If[res0 === $Failed, Return[$Failed]];
    {fVal, fidVal, penVal, gMax, extra, evalState} =
        {res0["Objective"], res0["Fidelity"], res0["Penalty"], res0["Gradient"], res0["Extra"], res0["EvalState"]};
    phiVal = -fVal;
    gMin = Developer`ToPackedArray[-gMax];
    gNorm = Norm[gMin];
    bestFid = fidVal;

    delta = With[{rawD = Lookup[subOpts, "InitialRadius", Automatic]},
        If[NumericQ[rawD], Clip[N[rawD], {1.*^-6, deltaMax}], Clip[gNorm, {0.5, 10.}]]
    ];
    sigma = Clip[N[Lookup[subOpts, "InitialSigma", 1.0]], {sigmaMin, sigmaMax}];

    objHist = If[recHist, {fVal}, {}];
    fidHist = If[recHist, {fidVal}, {}];
    gHist = If[recHist, {gNorm}, {}];

    With[{critFid0 = evalCritFn[fidVal, extra, fidCrit]},
        status = Which[
            NumericQ[targetFid] && TrueQ[critFid0 >= targetFid], "TargetFidelityReached",
            TrueQ[gNorm <= gradTol], "GradientToleranceReached",
            True, "Running"
        ]
    ];
    {lastMonTime, lastMonIter} = emitMonFn[monCfg, 0, evals, fVal, fidVal, bestFid, gNorm, 0., 0., status, extra, stageTag, lastMonTime, lastMonIter];
    If[status =!= "Running",
        With[{elapsed = AbsoluteTime[] - tStart},
            Return[<|
                "ParametersVec" -> Developer`ToPackedArray[th * invS],
                "Objective" -> fVal, "Fidelity" -> fidVal, "Penalty" -> penVal,
                "GradientNorm" -> gNorm, "ScaledGradient" -> gMax, "Optimizer" -> mode,
                "Iterations" -> 0, "Evaluations" -> evals, "ObjectiveEvaluations" -> objEvals,
                "GradientEvaluations" -> gradEvals, "HVPEvaluations" -> hvpEvals, "HessianEvaluations" -> 0,
                "CGIterations" -> cgIters, "NegativeCurvatureCount" -> negCurvCount,
                "TrustRegionFinalRadius" -> If[mode === "TrustRegionGLTR", delta, None],
                "ARCFinalSigma" -> If[mode === "AdaptiveCubicRegularization", sigma, None],
                "AcceptedSteps" -> 0, "RejectedSteps" -> 0, "SkippedCurvaturePairs" -> 0,
                "LineSearchMethod" -> If[mode === "AdaptiveCubicRegularization", "CubicRegularization", "LanczosGLTR"],
                "StrongWolfeSatisfiedCount" -> 0, "ApproxWolfeCount" -> 0,
                "LastStepLength" -> 0., "Status" -> status,
                "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
                "FidelityHistory" -> Developer`ToPackedArray[fidHist],
                "GradientNormHistory" -> Developer`ToPackedArray[gHist],
                "ElapsedTime" -> elapsed, "DynamicsTime" -> dynTime, "OptimizerTime" -> Max[0., elapsed - dynTime],
                "Extra" -> extra
            |>]
        ]
    ];
    status = "MaxIterationsReached";

    Do[
        tolScale = 1.*^-14 * Max[1., Abs[fVal]];
        krylovTol = Max[1.*^-12 * gNorm, Min[0.5, gNorm^forcingExp] * gNorm];
        qList = {Developer`ToPackedArray[gMin / Max[gNorm, 1.*^-30]]};
        diagT = {};
        subDiagT = {};
        betaPrev = 0.;
        qPrev = ConstantArray[0., nDim];
        subSol = None;

        (* Build reorthogonalised Lanczos Krylov subspace (Q_j, T_j) incrementally *)
        Do[
            qCur = qList[[j]];
            hv = evalHVPMinScaled[evalState, qCur];
            If[hv === $Failed, subSol = $Failed; Break[]];
            cgIters++;
            alphaJ = qCur . hv;
            AppendTo[diagT, alphaJ];
            wVec = hv - alphaJ * qCur - If[j > 1, betaPrev * qPrev, 0.];
            qMat = Developer`ToPackedArray[qList];
            (* Two-pass full reorthogonalisation against Q_j *)
            coeffs = qMat . wVec;
            wVec = wVec - (coeffs . qMat);
            coeffs = qMat . wVec;
            wVec = Developer`ToPackedArray[wVec - (coeffs . qMat)];
            betaNext = Norm[wVec];

            jDim = Length[diagT];
            tMat = Developer`ToPackedArray[If[jDim === 1,
                {{diagT[[1]]}},
                DiagonalMatrix[diagT] + DiagonalMatrix[subDiagT, 1] + DiagonalMatrix[subDiagT, -1]
            ]];
            gSmall = Developer`ToPackedArray[Prepend[ConstantArray[0., jDim - 1], gNorm]];

            subSol = If[mode === "TrustRegionGLTR",
                SolveSmallTrustRegionSubproblem[tMat, gSmall, delta],
                SolveSmallARCSubproblem[tMat, gSmall, sigma]
            ];
            hVec = subSol["Step"];
            resNorm = betaNext * Abs[Last[hVec]];
            With[{tolEffK = If[mode === "AdaptiveCubicRegularization",
                    Max[1.*^-12 * gNorm, Min[krylovTol, 0.1 * Min[1., Norm[hVec]] * gNorm]],
                    krylovTol
                 ]},
                If[TrueQ[resNorm <= tolEffK] || TrueQ[betaNext <= 1.*^-13 * Max[1., gNorm]] || j === maxKrylov,
                    Break[]
                ]
            ];
            AppendTo[subDiagT, betaNext];
            betaPrev = betaNext;
            qPrev = qCur;
            AppendTo[qList, Developer`ToPackedArray[wVec / betaNext]],
            {j, maxKrylov}
        ];

        If[subSol === $Failed || !AssociationQ[subSol],
            status = "LineSearchTerminated"; Break[]
        ];
        qMat = Developer`ToPackedArray[qList];
        If[TrueQ[subSol["NegativeCurvature"]], negCurvCount++];

        accepted = False;
        Do[
            If[retry > 1,
                subSol = If[mode === "TrustRegionGLTR",
                    SolveSmallTrustRegionSubproblem[tMat, gSmall, delta],
                    SolveSmallARCSubproblem[tMat, gSmall, sigma]
                ]
            ];
            hVec = subSol["Step"];
            p = Developer`ToPackedArray[hVec . qMat];
            pNorm = Norm[p];
            predDec = If[mode === "TrustRegionGLTR", -subSol["ModelValue"], -subSol["QuadModelValue"]];
            hitBoundary = Lookup[subSol, "OnBoundary", True];
            If[!TrueQ[predDec > 0.], predDec = -(gMin . p)];
            If[TrueQ[pNorm <= 1.*^-15 * Max[1., Norm[th]]],
                Break[]
            ];
            thTrial = Developer`ToPackedArray[th + p];
            resTrial = evalCacheScaled[thTrial];
            If[resTrial === $Failed,
                rejectedSteps++;
                If[mode === "TrustRegionGLTR",
                    delta = 0.25 * Max[Min[pNorm, delta], 1.*^-30];
                    If[TrueQ[delta < deltaMin * Max[1., Norm[th]]], Break[]],
                    sigma = Min[sigmaMax, 2.0 * sigma];
                    If[TrueQ[sigma >= sigmaMax], Break[]]
                ];
                Continue[]
            ];
            phiTrial = -resTrial["Objective"];
            ared = phiVal - phiTrial;
            tolEff = Min[tolScale, 0.25 * eta1 * Max[0., predDec]];
            rho = If[TrueQ[predDec > 1.*^-24], ared / predDec, -1.];

            If[TrueQ[rho >= eta1] && TrueQ[ared > -tolEff],
                accepted = True;
                If[mode === "TrustRegionGLTR",
                    If[TrueQ[rho < 0.25],
                        delta = 0.25 * Max[Min[pNorm, delta], 1.*^-30],
                        If[TrueQ[rho > eta2 && hitBoundary], delta = Min[deltaMax, 2.0 * delta]]
                    ],
                    If[TrueQ[rho >= eta2],
                        sigma = Max[sigmaMin, 0.5 * sigma]
                    ]
                ];
                Break[],
                rejectedSteps++;
                If[mode === "TrustRegionGLTR",
                    delta = 0.25 * Max[Min[pNorm, delta], 1.*^-30];
                    If[TrueQ[delta < deltaMin * Max[1., Norm[th]]], Break[]],
                    sigma = Min[sigmaMax, 2.0 * sigma];
                    If[TrueQ[sigma >= sigmaMax], Break[]]
                ]
            ],
            {retry, 12}
        ];

        If[!accepted,
            status = If[gNorm <= Sqrt[gradTol], "GradientToleranceReached", "StepToleranceReached"];
            Break[]
        ];

        acceptedSteps++;
        sVec = Developer`ToPackedArray[thTrial - th];
        lastStep = Norm[sVec];
        fOld = fVal;
        th = thTrial;
        {fVal, fidVal, penVal, gMax, extra, evalState} =
            {resTrial["Objective"], resTrial["Fidelity"], resTrial["Penalty"], resTrial["Gradient"], resTrial["Extra"], resTrial["EvalState"]};
        phiVal = -fVal;
        gMin = Developer`ToPackedArray[-gMax];
        gNorm = Norm[gMin];
        bestFid = Max[bestFid, fidVal];
        nIter++;

        If[TrueQ[verbose] && (Mod[nIter, 5] === 0 || nIter === 1),
            Print["Iter ", nIter, " (", mode, "): Fidelity=", NumberForm[fidVal, 6],
                " Objective=", NumberForm[fVal, 6], " hvp=", hvpEvals, " |grad|=", ScientificForm[gNorm, 3]]
        ];
        If[TrueQ[recHist],
            AppendTo[objHist, fVal];
            AppendTo[fidHist, fidVal];
            AppendTo[gHist, gNorm]
        ];

        With[{critFidNew = evalCritFn[fidVal, extra, fidCrit]},
            stopNow = Which[
                TrueQ[Abs[fVal - fOld] <= funTol * Max[1., Abs[fOld]]] && TrueQ[gNorm <= Sqrt[gradTol]],
                    If[NumericQ[targetFid] && TrueQ[critFidNew >= targetFid], "TargetFidelityReached", "FunctionToleranceReached"],
                TrueQ[Norm[sVec] <= stepTol * Max[1., Norm[th]]],
                    If[NumericQ[targetFid] && TrueQ[critFidNew >= targetFid], "TargetFidelityReached", "StepToleranceReached"],
                NumericQ[targetFid] && TrueQ[critFidNew >= targetFid],
                    "TargetFidelityReached",
                TrueQ[gNorm <= gradTol],
                    "GradientToleranceReached",
                AbsoluteTime[] - tStart >= tLimit,
                    "TimeConstraintReached",
                nIter >= maxIter,
                    "MaxIterationsReached",
                True,
                    None
            ]
        ];
        {lastMonTime, lastMonIter} = emitMonFn[
            monCfg, nIter, evals, fVal, fidVal, bestFid, gNorm, lastStep,
            AbsoluteTime[] - tStart, If[stopNow === None, "Running", stopNow],
            extra, stageTag, lastMonTime, lastMonIter
        ];
        If[stopNow =!= None, status = stopNow; Break[]],
        {maxIter}
    ];
    If[lastMonIter < nIter || (lastMonIter < 0 && nIter === 0),
        {lastMonTime, lastMonIter} = emitMonFn[monCfg, nIter, evals, fVal, fidVal, bestFid, gNorm, lastStep, AbsoluteTime[] - tStart, status, extra, stageTag, lastMonTime, lastMonIter]
    ];

    With[{elapsed = AbsoluteTime[] - tStart},
        <|
            "ParametersVec" -> Developer`ToPackedArray[th * invS],
            "Objective" -> fVal,
            "Fidelity" -> fidVal,
            "Penalty" -> penVal,
            "GradientNorm" -> gNorm,
            "ScaledGradient" -> gMax,
            "Optimizer" -> mode,
            "Iterations" -> nIter,
            "Evaluations" -> evals,
            "ObjectiveEvaluations" -> objEvals,
            "GradientEvaluations" -> gradEvals,
            "HVPEvaluations" -> hvpEvals,
            "HessianEvaluations" -> 0,
            "CGIterations" -> cgIters,
            "NegativeCurvatureCount" -> negCurvCount,
            "TrustRegionFinalRadius" -> If[mode === "TrustRegionGLTR", delta, None],
            "ARCFinalSigma" -> If[mode === "AdaptiveCubicRegularization", sigma, None],
            "AcceptedSteps" -> acceptedSteps,
            "RejectedSteps" -> rejectedSteps,
            "SkippedCurvaturePairs" -> 0,
            "LineSearchMethod" -> If[mode === "AdaptiveCubicRegularization", "CubicRegularization", "LanczosGLTR"],
            "StrongWolfeSatisfiedCount" -> 0,
            "ApproxWolfeCount" -> 0,
            "LastStepLength" -> lastStep,
            "Status" -> status,
            "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
            "FidelityHistory" -> Developer`ToPackedArray[fidHist],
            "GradientNormHistory" -> Developer`ToPackedArray[gHist],
            "ElapsedTime" -> elapsed,
            "DynamicsTime" -> dynTime,
            "OptimizerTime" -> Max[0., elapsed - dynTime],
            "Extra" -> extra
        |>
    ]
];

RunTrustRegionGLTR[ctx_Association] := runKrylovSubspaceSolver[ctx, "TrustRegionGLTR"];
RunAdaptiveCubicRegularization[ctx_Association] := runKrylovSubspaceSolver[ctx, "AdaptiveCubicRegularization"];

(* Hybrid Phase 1 L-BFGS + Principled Phase 2 Second-Order Handoff (Newton-CG or Newton-RFO) *)
runHybridSecondOrder[ctx_Association, phase2Opt_String] := Module[
    {optName, isDenseP2, soCtx, evalFlatFn = ctx["EvalFn"], evalObjOnlyFn = Lookup[ctx, "EvalObjFn", None],
     x0Vec = ctx["X0Vec"], optsList = ctx["OptsList"], subOpts = Lookup[ctx, "SubOpts", <||>],
     scaleVec, invS, nDim,
     mMem, maxIter, targetFid, fidCrit, gradTol, funTol, stepTol, tLimit, recHist, verbose,
     wolfeOpt, maxLsOpt, monOpt, monIntOpt, monCfg, stageTag, lastMonTime = 0., lastMonIter = -1,
     parseMonFn, emitMonFn, evalCritFn, resolveLSFn, swLSFn, lsModeP1, lsCfg,
     minP1, maxP1, switchGradTol, switchRelGrad, carryPrec,
     tStart, dynTime = 0., evals = 0, objEvals = 0, gradEvals = 0, acceptedSteps = 0, rejectedSteps = 0,
     skippedCurv = 0, swCount = 0, approxWolfeCnt = 0, nmState = <||>, struggleStreak = 0,
     evalScaledP1, evalObjScaledP1, th, res0, fVal, fidVal, penVal, gVal, gNorm, g0Norm, extra, bestFid,
     objHist, fidHist, gHist, recentObjHist, nIter = 0, status = "Running", stopNow, fOld, lastStep = 0., prevDeltaF = 0.,
     sList = {}, yList = {}, rhoList = {}, dir, initStep, tolScale, lsOut, accepted,
     thNew, fNew, fidNew, penNew, gNew, extraNew, sVec, yVec, sy,
     switchReason = None, switchIter = None, isFidelityScaleQ, basinReadyQ, slowProgQ,
     p2RemIter, p2RemTime, p2SubOpts, p2Ctx, p2Opts, p2Res, mergedObj, mergedFid, mergedG, totalElapsed},

    optName = If[phase2Opt === "TrustRegionNewtonCG", "HybridLBFGSNewtonCG", "HybridLBFGSNewtonRFO"];
    isDenseP2 = (phase2Opt === "NewtonRFO");
    soCtx = buildSecondOrderContext[ctx, optName, isDenseP2];
    If[soCtx === $Failed, Return[$Failed]];
    scaleVec = soCtx["ScaleVec"];
    invS = soCtx["InvScale"];
    nDim = soCtx["NumDim"];

    {mMem, maxIter, targetFid, gradTol, funTol, stepTol, tLimit, recHist, verbose,
     wolfeOpt, maxLsOpt, monOpt, monIntOpt} =
        {LBFGSMemory, MaxIterations, TargetFidelity, GradientTolerance, FunctionTolerance, StepTolerance,
         TimeConstraint, RecordHistory, Verbose, WolfeParameters, MaxLineSearchEvaluations,
         StepMonitor, MonitorInterval} /. optsList /. Options[GRAPEOptimize];
    fidCrit = Lookup[FilterRules[optsList, {$OCSolverRegistry["FidelityCriterionKey"]}], $OCSolverRegistry["FidelityCriterionKey"], Automatic];
    stageTag = Lookup[FilterRules[optsList, {$OCSolverRegistry["MonitorStageKey"]}], $OCSolverRegistry["MonitorStageKey"], "Optimization"];
    parseMonFn = $OCSolverRegistry["ParseMonitorConfig"];
    emitMonFn = $OCSolverRegistry["EmitMonitorIfDue"];
    evalCritFn = $OCSolverRegistry["EvalCriterionFidelity"];
    resolveLSFn = $OCSolverRegistry["ResolveLineSearchConfig"];
    swLSFn = $OCSolverRegistry["StrongWolfeLineSearch"];

    monCfg = parseMonFn[monOpt, monIntOpt];
    If[monCfg === $Failed, Return[$Failed]];
    lsModeP1 = Lookup[subOpts, "Phase1LineSearch", "StrongWolfe"];
    lsCfg = resolveLSFn[lsModeP1, wolfeOpt, maxLsOpt, "LBFGS", False];
    If[lsCfg === $Failed, Return[$Failed]];

    minP1 = Max[2, Round[Lookup[subOpts, "MinPhase1Iterations", 5]]];
    maxP1 = Max[minP1, Round[Lookup[subOpts, "MaxPhase1Iterations", Max[12, Ceiling[0.4 * maxIter]]]]];
    switchGradTol = Max[1.*^-12, N[Lookup[subOpts, "SwitchGradientNorm", 0.02]]];
    switchRelGrad = Clip[N[Lookup[subOpts, "SwitchRelativeGradient", 0.08]], {1.*^-6, 1.0}];
    carryPrec = TrueQ[Lookup[subOpts, "CarryPreconditioner", True]];
    isFidelityScaleQ = NumericQ[targetFid] && TrueQ[0. < targetFid <= 1.0];

    evalScaledP1[thIn_] := Module[{t0Dyn = AbsoluteTime[], r, gFlat},
        r = evalFlatFn[Developer`ToPackedArray[thIn * invS]];
        dynTime += AbsoluteTime[] - t0Dyn;
        evals++; objEvals++;
        If[r === $Failed || !AssociationQ[r] || !NumericQ[r["Objective"]],
            Return[$Failed]
        ];
        If[KeyExistsQ[r, "Gradient"], gradEvals++];
        gFlat = Developer`ToPackedArray[Flatten[r["Gradient"]] * invS];
        <|
            "Objective" -> N[Re[r["Objective"]]],
            "Fidelity" -> N[Re[r["Fidelity"]]],
            "Penalty" -> N[Re[Lookup[r, "Penalty", 0.]]],
            "Gradient" -> gFlat,
            "Extra" -> Lookup[r, "Extra", <||>]
        |>
    ];
    evalObjScaledP1[thIn_] := If[evalObjOnlyFn === None,
        evalScaledP1[thIn],
        Module[{t0Dyn = AbsoluteTime[], r},
            r = evalObjOnlyFn[Developer`ToPackedArray[thIn * invS]];
            dynTime += AbsoluteTime[] - t0Dyn;
            evals++; objEvals++;
            r
        ]
    ];

    tStart = AbsoluteTime[];
    th = Developer`ToPackedArray[N[x0Vec * scaleVec]];
    res0 = evalScaledP1[th];
    If[res0 === $Failed, Return[$Failed]];
    {fVal, fidVal, penVal, gVal, extra} =
        {res0["Objective"], res0["Fidelity"], res0["Penalty"], res0["Gradient"], res0["Extra"]};
    gNorm = Norm[gVal];
    g0Norm = gNorm;
    bestFid = fidVal;

    objHist = If[recHist, {fVal}, {}];
    fidHist = If[recHist, {fidVal}, {}];
    gHist = If[recHist, {gNorm}, {}];
    recentObjHist = {fVal};

    With[{critFid0 = evalCritFn[fidVal, extra, fidCrit]},
        status = Which[
            NumericQ[targetFid] && TrueQ[critFid0 >= targetFid], "TargetFidelityReached",
            TrueQ[gNorm <= gradTol], "GradientToleranceReached",
            True, "Running"
        ]
    ];
    {lastMonTime, lastMonIter} = emitMonFn[monCfg, 0, evals, fVal, fidVal, bestFid, gNorm, 0., 0., status, extra, stageTag, lastMonTime, lastMonIter];
    If[status =!= "Running",
        With[{elapsed = AbsoluteTime[] - tStart},
            Return[<|
                "ParametersVec" -> Developer`ToPackedArray[th * invS],
                "Objective" -> fVal, "Fidelity" -> fidVal, "Penalty" -> penVal,
                "GradientNorm" -> gNorm, "ScaledGradient" -> gVal, "Optimizer" -> optName,
                "Iterations" -> 0, "Evaluations" -> evals, "ObjectiveEvaluations" -> objEvals,
                "GradientEvaluations" -> gradEvals, "HVPEvaluations" -> 0, "HessianEvaluations" -> 0,
                "CGIterations" -> 0, "NegativeCurvatureCount" -> 0,
                "HybridSwitchIteration" -> None, "HybridSwitchReason" -> "ConvergedAtStart",
                "AcceptedSteps" -> 0, "RejectedSteps" -> 0, "SkippedCurvaturePairs" -> 0,
                "LineSearchMethod" -> lsCfg["Mode"], "StrongWolfeSatisfiedCount" -> 0, "ApproxWolfeCount" -> 0,
                "LastStepLength" -> 0., "Status" -> status,
                "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
                "FidelityHistory" -> Developer`ToPackedArray[fidHist],
                "GradientNormHistory" -> Developer`ToPackedArray[gHist],
                "ElapsedTime" -> elapsed, "DynamicsTime" -> dynTime, "OptimizerTime" -> Max[0., elapsed - dynTime],
                "Extra" -> extra
            |>]
        ]
    ];

    status = "Phase1Active";
    Do[
        dir = If[Length[sList] > 0,
            lbfgsTwoLoop[sList, yList, rhoList, gVal],
            gVal
        ];
        If[!TrueQ[dir . gVal > 1.*^-14 * Norm[dir] * Max[gNorm, 1.*^-30]],
            dir = gVal;
            sList = {}; yList = {}; rhoList = {}
        ];
        initStep = If[Length[sList] === 0, Min[1., 0.35 / Max[Norm[dir], 1.*^-12]], 1.0];
        tolScale = 1.*^-14 * Max[1., Abs[fVal]];

        lsOut = If[lsCfg["Mode"] === "StrongWolfe",
            swLSFn[evalScaledP1, th, fVal, gVal, dir, initStep, lsCfg["C1"], lsCfg["C2"], lsCfg["MaxEvals"], tolScale],
            With[{regLS = Lookup[$OCSolverRegistry, "LineSearches", <||>]},
                If[KeyExistsQ[regLS, lsCfg["Mode"]],
                    (regLS[lsCfg["Mode"]])[<|
                        "EvalScaled" -> evalScaledP1, "EvalObjScaled" -> evalObjScaledP1, "ProjBox" -> Identity,
                        "HasConstraints" -> False, "Theta" -> th, "FVal" -> fVal, "GVal" -> gVal,
                        "Dir" -> dir, "InitStep" -> initStep, "C1" -> lsCfg["C1"], "C2" -> lsCfg["C2"],
                        "MaxEvals" -> lsCfg["MaxEvals"], "TolScale" -> tolScale, "PrevDeltaF" -> prevDeltaF,
                        "NMState" -> nmState, "ExtraOpts" -> lsCfg["ExtraOpts"]
                    |>],
                    $Failed
                ]
            ]
        ];
        accepted = AssociationQ[lsOut] && TrueQ[lsOut["Accepted"]];
        If[AssociationQ[lsOut] && KeyExistsQ[lsOut, "NMState"], nmState = lsOut["NMState"]];
        rejectedSteps += If[AssociationQ[lsOut], Max[0, Lookup[lsOut, "Trials", If[accepted, 1, 0]] - If[accepted, 1, 0]], 1];

        If[!accepted,
            switchReason = "Phase1LineSearchStruggle";
            switchIter = nIter;
            Break[]
        ];

        acceptedSteps++;
        thNew = lsOut["Theta"];
        {fNew, fidNew, penNew, gNew, extraNew} =
            {lsOut["Res"]["Objective"], lsOut["Res"]["Fidelity"], lsOut["Res"]["Penalty"], lsOut["Res"]["Gradient"], lsOut["Res"]["Extra"]};
        lastStep = lsOut["StepLength"];
        If[TrueQ[lsOut["StrongWolfeSatisfied"]], swCount++];
        If[TrueQ[lsOut["ApproxWolfeSatisfied"]], approxWolfeCnt++];
        If[Lookup[lsOut, "Trials", 1] >= 5, struggleStreak++, struggleStreak = 0];

        sVec = Developer`ToPackedArray[thNew - th];
        yVec = Developer`ToPackedArray[gVal - gNew];
        sy = sVec . yVec;
        If[TrueQ[sy > 1.*^-12 * Norm[sVec] * Norm[yVec]],
            If[Length[sList] >= mMem, sList = Rest[sList]; yList = Rest[yList]; rhoList = Rest[rhoList]];
            AppendTo[sList, sVec]; AppendTo[yList, yVec]; AppendTo[rhoList, 1. / sy],
            skippedCurv++; struggleStreak++
        ];

        fOld = fVal;
        prevDeltaF = Max[fNew - fVal, 0.];
        th = thNew; fVal = fNew; fidVal = fidNew; penVal = penNew; gVal = gNew; extra = extraNew;
        gNorm = Norm[gVal];
        bestFid = Max[bestFid, fidVal];
        nIter++;

        AppendTo[recentObjHist, fVal];
        If[TrueQ[recHist],
            AppendTo[objHist, fVal];
            AppendTo[fidHist, fidVal];
            AppendTo[gHist, gNorm]
        ];

        With[{critFidNew = evalCritFn[fidVal, extra, fidCrit]},
            stopNow = Which[
                TrueQ[Abs[fVal - fOld] <= funTol * Max[1., Abs[fOld]]] && TrueQ[gNorm <= Sqrt[gradTol]],
                    If[NumericQ[targetFid] && TrueQ[critFidNew >= targetFid], "TargetFidelityReached", "FunctionToleranceReached"],
                TrueQ[Norm[sVec] <= stepTol * Max[1., Norm[th]]],
                    If[NumericQ[targetFid] && TrueQ[critFidNew >= targetFid], "TargetFidelityReached", "StepToleranceReached"],
                NumericQ[targetFid] && TrueQ[critFidNew >= targetFid],
                    "TargetFidelityReached",
                TrueQ[gNorm <= gradTol],
                    "GradientToleranceReached",
                AbsoluteTime[] - tStart >= tLimit,
                    "TimeConstraintReached",
                nIter >= maxIter,
                    "MaxIterationsReached",
                True,
                    None
            ]
        ];
        If[stopNow =!= None,
            status = stopNow;
            switchIter = None;
            switchReason = "ConvergedInPhase1";
            {lastMonTime, lastMonIter} = emitMonFn[monCfg, nIter, evals, fVal, fidVal, bestFid, gNorm, lastStep, AbsoluteTime[] - tStart, status, extra, stageTag, lastMonTime, lastMonIter];
            Break[]
        ];

        (* Evaluate switching triggers once minimum Phase 1 iterations are complete *)
        If[nIter >= minP1,
            basinReadyQ = TrueQ[gNorm <= switchGradTol || gNorm <= switchRelGrad * g0Norm] &&
                (!isFidelityScaleQ || TrueQ[fidVal >= 0.80] || nIter >= minP1 + 5);
            slowProgQ = (nIter >= minP1 + 5) && If[isFidelityScaleQ,
                TrueQ[(fVal - recentObjHist[[nIter - 4]]) <= 0.08 * Max[1.*^-12, 1. - fVal]],
                TrueQ[Abs[fVal - recentObjHist[[nIter - 4]]] <= 1.*^-3 * Max[1., Abs[fVal]]]
            ];
            Which[
                basinReadyQ,
                    switchReason = "BasinGradientTrigger"; switchIter = nIter; Break[],
                slowProgQ,
                    switchReason = "SlowLinearProgressTrigger"; switchIter = nIter; Break[],
                struggleStreak >= 2,
                    switchReason = "LineSearchStruggleTrigger"; switchIter = nIter; Break[],
                nIter >= maxP1,
                    switchReason = "MaxPhase1BudgetTrigger"; switchIter = nIter; Break[]
            ]
        ];
        {lastMonTime, lastMonIter} = emitMonFn[monCfg, nIter, evals, fVal, fidVal, bestFid, gNorm, lastStep, AbsoluteTime[] - tStart, "Running", extra, stageTag, lastMonTime, lastMonIter],
        {maxIter}
    ];

    If[switchReason === "ConvergedInPhase1" || status =!= "Phase1Active",
        With[{elapsed = AbsoluteTime[] - tStart},
            Return[<|
                "ParametersVec" -> Developer`ToPackedArray[th * invS],
                "Objective" -> fVal, "Fidelity" -> fidVal, "Penalty" -> penVal,
                "GradientNorm" -> gNorm, "ScaledGradient" -> gVal, "Optimizer" -> optName,
                "Iterations" -> nIter, "Evaluations" -> evals, "ObjectiveEvaluations" -> objEvals,
                "GradientEvaluations" -> gradEvals, "HVPEvaluations" -> 0, "HessianEvaluations" -> 0,
                "CGIterations" -> 0, "NegativeCurvatureCount" -> 0,
                "HybridSwitchIteration" -> switchIter, "HybridSwitchReason" -> switchReason,
                "AcceptedSteps" -> acceptedSteps, "RejectedSteps" -> rejectedSteps, "SkippedCurvaturePairs" -> skippedCurv,
                "LineSearchMethod" -> lsCfg["Mode"], "StrongWolfeSatisfiedCount" -> swCount, "ApproxWolfeCount" -> approxWolfeCnt,
                "LastStepLength" -> lastStep, "Status" -> status,
                "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
                "FidelityHistory" -> Developer`ToPackedArray[fidHist],
                "GradientNormHistory" -> Developer`ToPackedArray[gHist],
                "ElapsedTime" -> elapsed, "DynamicsTime" -> dynTime, "OptimizerTime" -> Max[0., elapsed - dynTime],
                "Extra" -> extra
            |>]
        ]
    ];

    If[TrueQ[verbose],
        Print["Hybrid (", optName, ") switching at iter ", switchIter, " (", switchReason, "): Fidelity=",
            NumberForm[fidVal, 6], " |grad|=", ScientificForm[gNorm, 3]]
    ];

    p2RemIter = Max[1, maxIter - nIter];
    p2RemTime = Max[0.01, tLimit - (AbsoluteTime[] - tStart)];
    p2SubOpts = Join[
        subOpts,
        If[phase2Opt === "TrustRegionNewtonCG",
            <|
                "Preconditioner" -> Lookup[subOpts, "Preconditioner", "LBFGS"],
                "InitialPreconditionerHistory" -> If[carryPrec, {sList, yList, rhoList}, None],
                "InitialRadius" -> Lookup[subOpts, "InitialRadius", Clip[2. * If[VectorQ[sVec, NumericQ], Norm[sVec], 0.5], {0.25, 5.0}]]
            |>,
            <||>
        ]
    ];
    p2Opts = Join[
        {
            MaxIterations -> p2RemIter,
            TimeConstraint -> p2RemTime,
            LineSearch -> Lookup[subOpts, "Phase2LineSearch", Automatic]
        },
        FilterRules[optsList, Except[MaxIterations | TimeConstraint | LineSearch | Optimizer]]
    ];
    p2Ctx = Join[ctx, <|
        "X0Vec" -> Developer`ToPackedArray[th * invS],
        "OptsList" -> p2Opts,
        "SubOpts" -> p2SubOpts
    |>];
    p2Res = If[phase2Opt === "TrustRegionNewtonCG",
        RunTrustRegionNewtonCG[p2Ctx],
        runExplicitNewtonFamily[p2Ctx, "NewtonRFO"]
    ];
    If[p2Res === $Failed, Return[$Failed]];

    mergedObj = If[recHist && Length[p2Res["ObjectiveHistory"]] >= 1, Join[objHist, Rest[p2Res["ObjectiveHistory"]]], {}];
    mergedFid = If[recHist && Length[p2Res["FidelityHistory"]] >= 1, Join[fidHist, Rest[p2Res["FidelityHistory"]]], {}];
    mergedG = If[recHist && Length[p2Res["GradientNormHistory"]] >= 1, Join[gHist, Rest[p2Res["GradientNormHistory"]]], {}];
    totalElapsed = AbsoluteTime[] - tStart;

    <|
        "ParametersVec" -> p2Res["ParametersVec"],
        "Objective" -> p2Res["Objective"],
        "Fidelity" -> p2Res["Fidelity"],
        "Penalty" -> p2Res["Penalty"],
        "GradientNorm" -> p2Res["GradientNorm"],
        "ScaledGradient" -> p2Res["ScaledGradient"],
        "Optimizer" -> optName,
        "Iterations" -> nIter + p2Res["Iterations"],
        "Evaluations" -> evals + p2Res["Evaluations"],
        "ObjectiveEvaluations" -> objEvals + p2Res["ObjectiveEvaluations"],
        "GradientEvaluations" -> gradEvals + p2Res["GradientEvaluations"],
        "HVPEvaluations" -> Lookup[p2Res, "HVPEvaluations", 0],
        "HessianEvaluations" -> Lookup[p2Res, "HessianEvaluations", 0],
        "CGIterations" -> Lookup[p2Res, "CGIterations", 0],
        "NegativeCurvatureCount" -> Lookup[p2Res, "NegativeCurvatureCount", 0],
        "HybridSwitchIteration" -> switchIter,
        "HybridSwitchReason" -> switchReason,
        "TrustRegionFinalRadius" -> Lookup[p2Res, "TrustRegionFinalRadius", None],
        "AcceptedSteps" -> acceptedSteps + Lookup[p2Res, "AcceptedSteps", 0],
        "RejectedSteps" -> rejectedSteps + Lookup[p2Res, "RejectedSteps", 0],
        "SkippedCurvaturePairs" -> skippedCurv + Lookup[p2Res, "SkippedCurvaturePairs", 0],
        "LineSearchMethod" -> (lsCfg["Mode"] <> "->" <> Lookup[p2Res, "LineSearchMethod", phase2Opt]),
        "StrongWolfeSatisfiedCount" -> swCount + Lookup[p2Res, "StrongWolfeSatisfiedCount", 0],
        "ApproxWolfeCount" -> approxWolfeCnt + Lookup[p2Res, "ApproxWolfeCount", 0],
        "LastStepLength" -> p2Res["LastStepLength"],
        "Status" -> p2Res["Status"],
        "ObjectiveHistory" -> Developer`ToPackedArray[mergedObj],
        "FidelityHistory" -> Developer`ToPackedArray[mergedFid],
        "GradientNormHistory" -> Developer`ToPackedArray[mergedG],
        "ElapsedTime" -> totalElapsed,
        "DynamicsTime" -> dynTime + Lookup[p2Res, "DynamicsTime", 0.],
        "OptimizerTime" -> Max[0., totalElapsed - (dynTime + Lookup[p2Res, "DynamicsTime", 0.])],
        "Extra" -> p2Res["Extra"]
    |>
];

RunHybridLBFGSNewtonCG[ctx_Association] := runHybridSecondOrder[ctx, "TrustRegionNewtonCG"];
RunHybridLBFGSNewtonRFO[ctx_Association] := runHybridSecondOrder[ctx, "NewtonRFO"];

(* Register all 9 Phase D second-order & hybrid optimizers plus algebraic subproblem primitives *)
If[!KeyExistsQ[$OCSolverRegistry, "Optimizers"] || !AssociationQ[$OCSolverRegistry["Optimizers"]],
    $OCSolverRegistry["Optimizers"] = <||>
];

$OCSolverRegistry["Optimizers"] = Join[$OCSolverRegistry["Optimizers"], <|
    "NewtonRFO" -> Function[c, runExplicitNewtonFamily[c, "NewtonRFO"]],
    "NewtonShift" -> Function[c, runExplicitNewtonFamily[c, "NewtonShift"]],
    "NewtonModifiedCholesky" -> Function[c, runExplicitNewtonFamily[c, "NewtonModifiedCholesky"]],
    "NewtonSpectral" -> Function[c, runExplicitNewtonFamily[c, "NewtonSpectral"]],
    "TrustRegionNewtonCG" -> RunTrustRegionNewtonCG,
    "TrustRegionGLTR" -> RunTrustRegionGLTR,
    "AdaptiveCubicRegularization" -> RunAdaptiveCubicRegularization,
    "HybridLBFGSNewtonCG" -> RunHybridLBFGSNewtonCG,
    "HybridLBFGSNewtonRFO" -> RunHybridLBFGSNewtonRFO
|>];

$OCSolverRegistry["ModifiedCholeskySolve"] = ModifiedCholeskySolve;
$OCSolverRegistry["SolveRFOSecular"] = SolveRFOSecular;
$OCSolverRegistry["SolveSmallTRSubproblem"] = SolveSmallTrustRegionSubproblem;
$OCSolverRegistry["SolveSmallARCSubproblem"] = SolveSmallARCSubproblem;

End[];
EndPackage[];
