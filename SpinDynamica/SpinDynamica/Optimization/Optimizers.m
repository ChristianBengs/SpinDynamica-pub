(* ::Package:: *)

(* SpinDynamica Optimal Control — Advanced Quasi-Newton Optimizers & Line Searches (Phase B) *)

BeginPackage["OptimalControl`Optimizers`", {
    "Global`",
    "SpinDynamica`",
    "OptimalControl`Kernels`",
    "OptimalControl`GRAPE`"
}];

Begin["`Private`"];

(* Cancellation-free Kanzow-Steck (2023) compact regularised L-BFGS linear solver:
   Solves (B_k + mu I) d = -gMin where B_k = delta I - W M^-1 W^T, W = [Y^T, delta S^T]. *)
CompactRegularizedLBFGSStep[sMat_, yMat_, gMin_List, delta_?NumericQ, mu_?NumericQ, gramOpt_: Automatic] := Module[
    {gVec, gNormSq, mCur, denom0, ssMat, syMat, yyMat, dDiag, lMat, uMat,
     tl, bl, tr, br, kMat, rhs, zSol, mIdx, z1, z2, dVec, gd, dd, pred},
    gVec = Developer`ToPackedArray[N[Re[gMin]]];
    gNormSq = gVec . gVec;
    denom0 = N[delta] + N[mu];
    If[!TrueQ[delta > 0. && mu > 0. && denom0 > 0.], Return[$Failed]];
    If[sMat === {} || Length[sMat] === 0,
        dVec = Developer`ToPackedArray[(-1. / denom0) * gVec];
        gd = gVec . dVec;
        dd = dVec . dVec;
        pred = -0.5 * gd + 0.5 * mu * dd;
        Return[<|"Direction" -> dVec, "PredDecrease" -> pred, "Delta" -> N[delta], "Mu" -> N[mu]|>]
    ];

    mCur = Length[sMat];
    If[AssociationQ[gramOpt],
        ssMat = gramOpt["SS"];
        syMat = gramOpt["SY"];
        yyMat = gramOpt["YY"],
        ssMat = sMat . Transpose[sMat];
        syMat = sMat . Transpose[yMat];
        yyMat = yMat . Transpose[yMat]
    ];

    dDiag = DiagonalMatrix[Diagonal[syMat]];
    lMat = LowerTriangularize[syMat, -1];
    uMat = syMat - lMat;
    tl = -denom0 * dDiag - yyMat;
    bl = N[mu] * lMat - N[delta] * uMat;
    tr = Transpose[bl];
    br = (N[delta] * N[mu]) * ssMat;
    kMat = Developer`ToPackedArray[0.5 * (# + Transpose[#])] &@ ArrayFlatten[{{tl, tr}, {bl, br}}];
    rhs = Developer`ToPackedArray[Join[yMat . gVec, N[delta] * (sMat . gVec)]];

    zSol = Quiet[Check[LinearSolve[kMat, rhs], $Failed]];
    If[!VectorQ[zSol, NumericQ] || Length[zSol] =!= 2 * mCur || !AllTrue[zSol, NumberQ[#] && Abs[#] < 1.*^100 &],
        Return[$Failed]
    ];

    mIdx = mCur;
    z1 = zSol[[1 ;; mIdx]];
    z2 = zSol[[mIdx + 1 ;; 2 * mIdx]];
    dVec = Developer`ToPackedArray[(-1. / denom0) * (gVec + (z1 . yMat) + N[delta] * (z2 . sMat))];
    gd = gVec . dVec;
    dd = dVec . dVec;
    If[!TrueQ[gd <= -1.*^-14 * Sqrt[Max[gNormSq, 1.*^-30]] * Sqrt[Max[dd, 1.*^-30]]],
        Return[$Failed]
    ];
    pred = -0.5 * gd + 0.5 * N[mu] * dd;
    <|"Direction" -> dVec, "PredDecrease" -> pred, "Delta" -> N[delta], "Mu" -> N[mu]|>
];

(* Helper for secant interpolation safely inside an interval [a, b] *)
hzSecantPoint[p1_Association, p2_Association] := Module[
    {a1 = p1["Alpha"], g1 = p1["DPhi"], a2 = p2["Alpha"], g2 = p2["DPhi"], lo, hi, span, denom, c},
    lo = Min[a1, a2]; hi = Max[a1, a2]; span = hi - lo;
    If[span <= 1.*^-18, Return[0.5 * (lo + hi)]];
    denom = g2 - g1;
    If[!TrueQ[Abs[denom] > 1.*^-24 * Max[1., Abs[g1], Abs[g2]]], Return[0.5 * (lo + hi)]];
    c = (a1 * g2 - a2 * g1) / denom;
    If[NumericQ[c] && TrueQ[lo + 1.*^-12 * span < c < hi - 1.*^-12 * span],
        c,
        0.5 * (lo + hi)
    ]
];

(* Hager & Zhang (2005 / ACM TOMS Alg 851) line search with safeguarded approximate Wolfe conditions *)
HagerZhangLineSearch[lsCtx_Association] := Module[
    {evalScaled = lsCtx["EvalScaled"], th = lsCtx["Theta"], fVal = lsCtx["FVal"], gVal = lsCtx["GVal"],
     dir = lsCtx["Dir"], initStep = lsCtx["InitStep"], delta = lsCtx["C1"], sigma = lsCtx["C2"],
     maxLsEvals = lsCtx["MaxEvals"], prevDeltaF = Lookup[lsCtx, "PrevDeltaF", 0.],
     extraOpts = Lookup[lsCtx, "ExtraOpts", <||>],
     phi0, dphi0, epsHZ, epsMax, epsK, gammaBis, rhoExp, thetaHZ, alphaMax, lsEvals = 0,
     satisfiesHZQ, approxOnlyQ, bestSafe = None, evalPhiHZ, updateIntervalHZ,
     p0, pA, pB, pPrev, pCur, alpha, chosen = None, awUsed = False,
     pA1, pB1, pC, cTry, c2Try, pC2, pMid, oldWidth},

    phi0 = -fVal;
    dphi0 = -(gVal . dir);
    If[!TrueQ[dphi0 < 0.], Return[$Failed]];
    delta = Clip[N[Lookup[extraOpts, "Delta", delta]], {1.*^-8, 0.49}];
    sigma = Clip[N[Lookup[extraOpts, "Sigma", sigma]], {delta, 0.999999}];
    epsHZ = N[Lookup[extraOpts, "EpsScale", 1.*^-6]];
    epsMax = N[Lookup[extraOpts, "EpsMax", 1.*^-11]];
    gammaBis = N[Lookup[extraOpts, "Gamma", 0.66]];
    rhoExp = N[Lookup[extraOpts, "Rho", 5.0]];
    thetaHZ = N[Lookup[extraOpts, "Theta", 0.5]];
    alphaMax = N[Lookup[extraOpts, "AlphaMax", 10000.]];

    (* Guarded epsilon_k: never enables approximate Wolfe early when J ~ 0.99 while phi0 ~ -1 *)
    epsK = Clip[
        epsHZ * Max[Abs[prevDeltaF], 1.*^-11 * Max[1., Abs[phi0]]],
        {100. * $MachineEpsilon * Max[1., Abs[phi0]], epsMax * Max[1., Abs[phi0]]}
    ];

    satisfiesHZQ[pt_Association] := With[{a = pt["Alpha"], f = pt["Phi"], g = pt["DPhi"]},
        (TrueQ[f <= phi0 + delta * a * dphi0] && TrueQ[g >= sigma * dphi0]) ||
        (TrueQ[f <= phi0 + epsK] && TrueQ[sigma * dphi0 <= g <= (2. * delta - 1.) * dphi0])
    ];
    approxOnlyQ[pt_Association] := With[{a = pt["Alpha"], f = pt["Phi"], g = pt["DPhi"]},
        !(TrueQ[f <= phi0 + delta * a * dphi0] && TrueQ[g >= sigma * dphi0]) &&
        (TrueQ[f <= phi0 + epsK] && TrueQ[sigma * dphi0 <= g <= (2. * delta - 1.) * dphi0])
    ];

    evalPhiHZ[a_?NumericQ] := Module[{res, pt},
        lsEvals++;
        res = evalScaled[Developer`ToPackedArray[th + a * dir]];
        If[res === $Failed || !NumericQ[res["Objective"]] || !VectorQ[res["Gradient"], NumericQ],
            Return[$Failed]
        ];
        pt = <|"Alpha" -> a, "Phi" -> -res["Objective"], "DPhi" -> -(res["Gradient"] . dir), "Res" -> res|>;
        If[TrueQ[pt["Phi"] <= phi0 + epsK] && (bestSafe === None || TrueQ[pt["Phi"] < bestSafe["Phi"]]),
            bestSafe = pt
        ];
        pt
    ];

    (* Algorithm update + U3 opposite-slope bisection from Hager & Zhang (2005) *)
    updateIntervalHZ[aIn_Association, bIn_Association, cPt_Association] := Module[
        {aPt = aIn, bPt = bIn, aBar, bBar, dAlpha, pD},
        If[!TrueQ[aPt["Alpha"] < cPt["Alpha"] < bPt["Alpha"]],
            Return[{aPt, bPt}]
        ];
        If[TrueQ[cPt["DPhi"] >= 0.], Return[{aPt, cPt}]];
        If[TrueQ[cPt["Phi"] <= phi0 + epsK], Return[{cPt, bPt}]];

        aBar = aPt; bBar = cPt;
        While[lsEvals < maxLsEvals,
            If[Abs[bBar["Alpha"] - aBar["Alpha"]] <= 1.*^-14 * Max[1., aBar["Alpha"], bBar["Alpha"]], Break[]];
            dAlpha = (1. - thetaHZ) * aBar["Alpha"] + thetaHZ * bBar["Alpha"];
            pD = evalPhiHZ[dAlpha];
            If[pD === $Failed,
                bBar = <|"Alpha" -> dAlpha, "Phi" -> Indeterminate, "DPhi" -> Indeterminate|>;
                Continue[]
            ];
            If[satisfiesHZQ[pD], chosen = pD; Break[]];
            If[TrueQ[pD["DPhi"] >= 0.],
                Return[{aBar, pD}],
                If[TrueQ[pD["Phi"] <= phi0 + epsK],
                    aBar = pD,
                    bBar = pD
                ]
            ]
        ];
        {aBar, bBar}
    ];

    p0 = <|"Alpha" -> 0., "Phi" -> phi0, "DPhi" -> dphi0|>;
    pPrev = p0;
    alpha = Min[alphaMax, Max[initStep, 1.*^-12]];
    pA = None; pB = None;

    (* Bracketing expansion phase *)
    Do[
        If[lsEvals >= maxLsEvals, Break[]];
        pCur = evalPhiHZ[alpha];
        If[pCur === $Failed,
            alpha = 0.5 * (pPrev["Alpha"] + alpha);
            If[alpha <= 1.*^-14, Break[]];
            Continue[]
        ];
        If[satisfiesHZQ[pCur],
            chosen = pCur; Break[]
        ];
        If[TrueQ[pCur["DPhi"] >= 0.] || TrueQ[pCur["Phi"] > phi0 + epsK],
            pA = pPrev; pB = pCur; Break[]
        ];
        pPrev = pCur;
        If[alpha >= alphaMax, chosen = pCur; Break[]];
        alpha = Min[alphaMax, rhoExp * alpha],
        {lsIter, maxLsEvals}
    ];

    (* secant2 + bisection contraction loop *)
    If[chosen === None && AssociationQ[pA] && AssociationQ[pB],
        While[lsEvals < maxLsEvals && chosen === None,
            oldWidth = pB["Alpha"] - pA["Alpha"];
            If[oldWidth <= 1.*^-14 * Max[1., pA["Alpha"], pB["Alpha"]], Break[]];
            cTry = hzSecantPoint[pA, pB];
            pC = evalPhiHZ[cTry];
            If[pC === $Failed,
                pB = <|"Alpha" -> cTry, "Phi" -> pB["Phi"], "DPhi" -> pB["DPhi"]|>;
                Continue[]
            ];
            If[satisfiesHZQ[pC], chosen = pC; Break[]];
            {pA1, pB1} = updateIntervalHZ[pA, pB, pC];
            If[chosen =!= None || lsEvals >= maxLsEvals, Break[]];

            c2Try = Which[
                TrueQ[Abs[pA1["Alpha"] - pC["Alpha"]] <= 1.*^-18], hzSecantPoint[pA, pA1],
                TrueQ[Abs[pB1["Alpha"] - pC["Alpha"]] <= 1.*^-18], hzSecantPoint[pB, pB1],
                True, Indeterminate
            ];
            If[NumericQ[c2Try] && TrueQ[pA1["Alpha"] < c2Try < pB1["Alpha"]] && lsEvals < maxLsEvals,
                pC2 = evalPhiHZ[c2Try];
                If[AssociationQ[pC2],
                    If[satisfiesHZQ[pC2], chosen = pC2; Break[]];
                    {pA, pB} = updateIntervalHZ[pA1, pB1, pC2],
                    pA = pA1; pB = pB1
                ],
                pA = pA1; pB = pB1
            ];
            If[chosen === None && lsEvals < maxLsEvals && (pB["Alpha"] - pA["Alpha"]) > gammaBis * oldWidth,
                pMid = evalPhiHZ[0.5 * (pA["Alpha"] + pB["Alpha"])];
                If[AssociationQ[pMid],
                    If[satisfiesHZQ[pMid], chosen = pMid; Break[]];
                    {pA, pB} = updateIntervalHZ[pA, pB, pMid]
                ]
            ]
        ]
    ];

    If[chosen === None && AssociationQ[bestSafe] && TrueQ[bestSafe["Phi"] < phi0],
        chosen = bestSafe
    ];
    If[chosen === None,
        $Failed,
        awUsed = approxOnlyQ[chosen];
        <|
            "Accepted" -> True,
            "Theta" -> Developer`ToPackedArray[th + chosen["Alpha"] * dir],
            "Res" -> chosen["Res"],
            "StepLength" -> chosen["Alpha"],
            "StrongWolfeSatisfied" -> (TrueQ[satisfiesHZQ[chosen]] && !awUsed),
            "ApproxWolfeSatisfied" -> awUsed,
            "Trials" -> lsEvals
        |>
    ]
];

(* Shared Grippo-Lampariello-Lucidi (1986) & Zhang-Hager (2004) nonmonotone line search *)
runNonmonotoneLineSearch[lsCtx_Association, mode_String] := Module[
    {evalScaled = lsCtx["EvalScaled"], evalObjScaled = lsCtx["EvalObjScaled"],
     projBox = lsCtx["ProjBox"], hasConstr = TrueQ[lsCtx["HasConstraints"]],
     th = lsCtx["Theta"], fVal = lsCtx["FVal"], gVal = lsCtx["GVal"], dir = lsCtx["Dir"],
     initStep = lsCtx["InitStep"], c1 = lsCtx["C1"], c2 = lsCtx["C2"],
     maxLsEvals = lsCtx["MaxEvals"], tolScale = lsCtx["TolScale"],
     extraOpts = Lookup[lsCtx, "ExtraOpts", <||>], nmState = lsCtx["NMState"],
     etaZH, memGLL, useCurv, phi0, dphi0, phiRef, qPrev, cPrev, bufPrev,
     alphaMax = 10000., lsEvals = 0, step, trialRes, objTrial, gradRes, thTrial, fTrial, phiTrial, gs, tolEff,
     dphiTrial, d0, denom, stepQ, pPrev, pLo, pHi, pJ, alphaJ, chosen = None, swMet = False,
     updateNMState, nmArmijoQ, curvQ, evalObjOrFull},

    phi0 = -fVal;
    dphi0 = -(gVal . dir);
    etaZH = Clip[N[Lookup[extraOpts, "Eta", 0.85]], {0., 0.9999}];
    memGLL = Max[1, Round[Lookup[extraOpts, "Memory", 8]]];
    useCurv = If[hasConstr, False, TrueQ[Lookup[extraOpts, "Curvature", True]]];

    If[mode === "ZhangHagerNonmonotone",
        If[!AssociationQ[nmState],
            qPrev = 1.; cPrev = phi0,
            qPrev = N[Lookup[nmState, "Q", 1.]]; cPrev = N[Lookup[nmState, "C", phi0]]
        ];
        phiRef = Max[cPrev, phi0],
        bufPrev = If[AssociationQ[nmState], N[Lookup[nmState, "Buffer", {phi0}]], {phi0}];
        phiRef = Max[Append[bufPrev, phi0]]
    ];

    updateNMState[phiAcc_?NumericQ] := If[mode === "ZhangHagerNonmonotone",
        With[{qNext = etaZH * qPrev + 1.},
            <|"Q" -> qNext, "C" -> (etaZH * qPrev * cPrev + phiAcc) / qNext|>
        ],
        With[{newBuf = Append[bufPrev, phiAcc]},
            <|"Buffer" -> If[Length[newBuf] > memGLL, Take[newBuf, -memGLL], newBuf]|>
        ]
    ];

    (* Backtracking mode (used when projected OR when Curvature -> False) *)
    If[!useCurv,
        evalObjOrFull[thVec_] := If[evalObjScaled =!= None, evalObjScaled[thVec], evalScaled[thVec]];
        step = Max[initStep, 1.*^-12];
        Do[
            lsEvals++;
            thTrial = Developer`ToPackedArray[projBox[th + step * dir]];
            trialRes = evalObjOrFull[thTrial];
            If[trialRes =!= $Failed && NumericQ[trialRes["Objective"]],
                fTrial = trialRes["Objective"];
                phiTrial = -fTrial;
                gs = gVal . (thTrial - th);
                tolEff = If[TrueQ[gs > 0.], Min[tolScale, 0.25 * c1 * gs], 0.];
                If[TrueQ[gs > 0.] && TrueQ[phiTrial <= phiRef - c1 * gs + tolEff],
                    gradRes = If[KeyExistsQ[trialRes, "Gradient"], trialRes, evalScaled[thTrial]];
                    If[gradRes === $Failed, Return[$Failed, Module]];
                    chosen = <|"Alpha" -> step, "Theta" -> thTrial, "Phi" -> -gradRes["Objective"], "Res" -> gradRes|>;
                    Break[]
                ];
                If[!hasConstr,
                    d0 = -gs / Max[step, 1.*^-20];
                    denom = 2. * (phiTrial - phi0 - d0 * step);
                    stepQ = If[TrueQ[d0 < 0. && denom > 1.*^-20], -d0 * step^2 / denom, 0.5 * step];
                    step = Clip[stepQ, {0.25 * step, 0.5 * step}],
                    step *= 0.5
                ],
                step *= 0.5
            ];
            If[step < 1.*^-13, Break[]],
            {Max[maxLsEvals, 25]}
        ];
        If[chosen === None, Return[$Failed]];
        Return[<|
            "Accepted" -> True,
            "Theta" -> chosen["Theta"],
            "Res" -> chosen["Res"],
            "StepLength" -> chosen["Alpha"],
            "StrongWolfeSatisfied" -> False,
            "Trials" -> lsEvals,
            "NMState" -> updateNMState[chosen["Phi"]]
        |>];
    ];

    (* Unconstrained nonmonotone Wolfe mode: uses phiRef with weak Wolfe curvature dphi >= c2 dphi0 *)
    If[!TrueQ[dphi0 < 0.], Return[$Failed]];
    nmArmijoQ[fT_, a_] := With[{tEff = Min[tolScale, 0.25 * Abs[c1 * a * dphi0]]},
        TrueQ[fT <= phiRef + c1 * a * dphi0 + tEff]
    ];
    curvQ[dfT_] := TrueQ[dfT >= c2 * dphi0];

    pPrev = <|"Alpha" -> 0., "Phi" -> phi0, "DPhi" -> dphi0|>;
    step = Min[alphaMax, Max[initStep, 1.*^-12]];
    Do[
        If[lsEvals >= maxLsEvals, Break[]];
        lsEvals++;
        thTrial = Developer`ToPackedArray[th + step * dir];
        gradRes = evalScaled[thTrial];
        If[gradRes === $Failed || !NumericQ[gradRes["Objective"]] || !VectorQ[gradRes["Gradient"], NumericQ],
            step = 0.5 * (pPrev["Alpha"] + step);
            If[step <= 1.*^-13, Break[]];
            Continue[]
        ];
        objTrial = <|"Alpha" -> step, "Phi" -> -gradRes["Objective"], "DPhi" -> -(gradRes["Gradient"] . dir), "Res" -> gradRes|>;
        If[!nmArmijoQ[objTrial["Phi"], step],
            pLo = pPrev; pHi = objTrial;
            While[lsEvals < maxLsEvals,
                If[Abs[pHi["Alpha"] - pLo["Alpha"]] <= 1.*^-13 * Max[1., pLo["Alpha"], pHi["Alpha"]], Break[]];
                alphaJ = If[NumericQ[pHi["Phi"]] && NumericQ[pHi["DPhi"]],
                    0.5 * (pLo["Alpha"] + pHi["Alpha"]), (* safe fallback; try quadratic first *)
                    0.5 * (pLo["Alpha"] + pHi["Alpha"])
                ];
                If[NumericQ[pHi["Phi"]],
                    denom = 2. * (pHi["Phi"] - pLo["Phi"] - pLo["DPhi"] * (pHi["Alpha"] - pLo["Alpha"]));
                    If[TrueQ[denom > 1.*^-24],
                        stepQ = pLo["Alpha"] - pLo["DPhi"] * (pHi["Alpha"] - pLo["Alpha"])^2 / denom;
                        With[{lo = Min[pLo["Alpha"], pHi["Alpha"]], hi = Max[pLo["Alpha"], pHi["Alpha"]], sp = Abs[pHi["Alpha"] - pLo["Alpha"]]},
                            If[NumericQ[stepQ] && TrueQ[lo + 0.1 * sp <= stepQ <= hi - 0.1 * sp], alphaJ = stepQ]
                        ]
                    ]
                ];
                lsEvals++;
                pJ = With[{r = evalScaled[Developer`ToPackedArray[th + alphaJ * dir]]},
                    If[r === $Failed || !NumericQ[r["Objective"]] || !VectorQ[r["Gradient"], NumericQ],
                        $Failed,
                        <|"Alpha" -> alphaJ, "Phi" -> -r["Objective"], "DPhi" -> -(r["Gradient"] . dir), "Res" -> r|>
                    ]
                ];
                If[pJ === $Failed,
                    pHi = <|"Alpha" -> alphaJ, "Phi" -> Indeterminate, "DPhi" -> Indeterminate|>;
                    Continue[]
                ];
                If[!nmArmijoQ[pJ["Phi"], alphaJ],
                    pHi = pJ,
                    If[curvQ[pJ["DPhi"]],
                        chosen = pJ; swMet = TrueQ[Abs[pJ["DPhi"]] <= -c2 * dphi0]; Break[]
                    ];
                    pLo = pJ;
                    If[TrueQ[pJ["Alpha"] >= pHi["Alpha"]], Break[]]
                ]
            ];
            If[chosen === None && KeyExistsQ[pLo, "Res"] && nmArmijoQ[pLo["Phi"], pLo["Alpha"]], chosen = pLo];
            Break[]
        ];
        If[curvQ[objTrial["DPhi"]],
            chosen = objTrial; swMet = TrueQ[Abs[objTrial["DPhi"]] <= -c2 * dphi0]; Break[]
        ];
        pPrev = objTrial;
        If[step >= alphaMax, chosen = objTrial; Break[]];
        step = Min[alphaMax, 2.0 * step],
        {lsIter, maxLsEvals}
    ];
    If[chosen === None && KeyExistsQ[pPrev, "Res"] && nmArmijoQ[pPrev["Phi"], pPrev["Alpha"]], chosen = pPrev];
    If[chosen === None, Return[$Failed]];
    <|
        "Accepted" -> True,
        "Theta" -> Developer`ToPackedArray[th + chosen["Alpha"] * dir],
        "Res" -> chosen["Res"],
        "StepLength" -> chosen["Alpha"],
        "StrongWolfeSatisfied" -> swMet,
        "Trials" -> lsEvals,
        "NMState" -> updateNMState[chosen["Phi"]]
    |>
];

ZhangHagerNonmonotoneLineSearch[lsCtx_Association] :=
    runNonmonotoneLineSearch[lsCtx, "ZhangHagerNonmonotone"];

GrippoNonmonotoneLineSearch[lsCtx_Association] :=
    runNonmonotoneLineSearch[lsCtx, "GrippoNonmonotone"];

(* Kanzow & Steck (2023) compact regularised L-BFGS driver *)
RunRegularizedLBFGS[ctx_Association] := Module[
    {evalFlatFn = ctx["EvalFn"], evalObjOnlyFn = Lookup[ctx, "EvalObjFn", None],
     x0Vec = ctx["X0Vec"], scaleVec = ctx["ScaleVec"], lbVec = ctx["LBVec"], ubVec = ctx["UBVec"],
     hasBounds = TrueQ[ctx["HasBounds"]], customProjFn = Lookup[ctx, "CustomProjFn", None],
     hasCustomConstr = TrueQ[ctx["HasCustomConstr"]], anyConstr = TrueQ[ctx["AnyConstraints"]],
     optsList = ctx["OptsList"], subOpts = Lookup[ctx, "SubOpts", <||>],
     allowProj, accMode, etaZH, sigmaInit, sigmaMin, sigmaMax, gradPow,
     eta1, eta2, gamma1, gamma2, maxInner, useObjTrial,
     mMem, maxIter, targetFid, fidCrit, gradTol, funTol, stepTol, tLimit, recHist, verbose,
     monOpt, monIntOpt, monCfg, stageTag, lastMonTime = 0., lastMonIter = -1, bestFid,
     parseMonFn, emitMonFn, evalCritFn,
     projBox, tStart, dynTime = 0., lbT, ubT, nDim, evalScaled, evalObjScaled, evalTrialFn,
     th, fVal, phiVal, fidVal, penVal, gVal, gMin, gMinNorm, extra, pg, pgNorm,
     evals = 0, objEvals = 0, gradEvals = 0, acceptedSteps = 0, rejectedSteps = 0, skippedCurv = 0,
     objHist, fidHist, gHist, nIter = 0, status = "MaxIterationsReached", stopNow,
     sMat = {}, yMat = {}, ssMat = {}, syMat = {}, yyMat = {}, appendCurvaturePair,
     qNM = 1., cNM, phiRef, sigma, delta, mu, stepOut, dVec, predDec,
     thTrial, resTrial, fTrial, phiTrial, ared, rhoRatio, accepted, gradRes,
     sVec, yVec, sy, syLast, yyLast, fOld, lastStep = 0., tolScale},

    allowProj = TrueQ[Lookup[subOpts, "AllowProjected", False]];
    If[anyConstr && !allowProj,
        Message[OptimalControlProblem::unconstronly, "RegularizedLBFGS"];
        Return[$Failed]
    ];

    accMode = Lookup[subOpts, "Acceptance", "ZhangHager"];
    etaZH = Clip[N[Lookup[subOpts, "Eta", 0.85]], {0., 0.9999}];
    sigmaInit = Max[N[Lookup[subOpts, "InitialSigma", 0.1]], 1.*^-12];
    sigmaMin = Max[N[Lookup[subOpts, "MinSigma", 1.*^-10]], 1.*^-14];
    sigmaMax = Min[N[Lookup[subOpts, "MaxSigma", 1.*^10]], 1.*^14];
    gradPow = Clip[N[Lookup[subOpts, "GradientPower", 1.0]], {0., 2.}];
    eta1 = Clip[N[Lookup[subOpts, "Eta1", 1.*^-4]], {1.*^-8, 0.25}];
    eta2 = Clip[N[Lookup[subOpts, "Eta2", 0.75]], {eta1 + 0.05, 0.99}];
    gamma1 = Clip[N[Lookup[subOpts, "Gamma1", 0.5]], {0.1, 0.9}];
    gamma2 = Clip[N[Lookup[subOpts, "Gamma2", 4.0]], {1.5, 10.}];
    maxInner = Max[5, Round[Lookup[subOpts, "MaxInnerTrials", 25]]];
    useObjTrial = TrueQ[Lookup[subOpts, "UseObjectiveOnlyTrial", False]] && (evalObjOnlyFn =!= None);

    {mMem, maxIter, targetFid, gradTol, funTol, stepTol, tLimit, recHist, verbose, monOpt, monIntOpt} =
        {LBFGSMemory, MaxIterations, TargetFidelity, GradientTolerance, FunctionTolerance, StepTolerance,
         TimeConstraint, RecordHistory, Verbose, StepMonitor, MonitorInterval} /. optsList /. Options[GRAPEOptimize];
    fidCrit = Lookup[FilterRules[optsList, {$OCSolverRegistry["FidelityCriterionKey"]}], $OCSolverRegistry["FidelityCriterionKey"], Automatic];
    stageTag = Lookup[FilterRules[optsList, {$OCSolverRegistry["MonitorStageKey"]}], $OCSolverRegistry["MonitorStageKey"], "Optimization"];
    parseMonFn = $OCSolverRegistry["ParseMonitorConfig"];
    emitMonFn = $OCSolverRegistry["EmitMonitorIfDue"];
    evalCritFn = $OCSolverRegistry["EvalCriterionFidelity"];
    monCfg = parseMonFn[monOpt, monIntOpt];
    If[monCfg === $Failed, Return[$Failed]];

    nDim = Length[x0Vec];
    tStart = AbsoluteTime[];
    lbT = If[Length[lbVec] === nDim, Developer`ToPackedArray[lbVec * scaleVec], ConstantArray[-Infinity, nDim]];
    ubT = If[Length[ubVec] === nDim, Developer`ToPackedArray[ubVec * scaleVec], ConstantArray[Infinity, nDim]];
    projBox[v_] := Which[
        customProjFn =!= None, Developer`ToPackedArray[customProjFn[v]],
        hasBounds, MapThread[Clip[#1, {#2, #3}] &, {v, lbT, ubT}],
        True, v
    ];

    evalScaled[thVec_] := Module[{tD0, rawEval, gRaw},
        evals++; objEvals++; gradEvals++;
        tD0 = AbsoluteTime[];
        rawEval = evalFlatFn[Developer`ToPackedArray[thVec / scaleVec]];
        dynTime += AbsoluteTime[] - tD0;
        If[rawEval === $Failed, Return[$Failed]];
        gRaw = Developer`ToPackedArray[Flatten[rawEval["Gradient"]] / scaleVec];
        <|
            "Objective" -> rawEval["Objective"],
            "Fidelity" -> rawEval["Fidelity"],
            "Penalty" -> rawEval["Penalty"],
            "Gradient" -> gRaw,
            "Extra" -> KeyDrop[rawEval, {"Objective", "Fidelity", "Penalty", "Gradient"}]
        |>
    ];
    evalObjScaled[thVec_] := Module[{tD0, rawObj},
        evals++; objEvals++;
        tD0 = AbsoluteTime[];
        rawObj = evalObjOnlyFn[Developer`ToPackedArray[thVec / scaleVec]];
        dynTime += AbsoluteTime[] - tD0;
        If[rawObj === $Failed, Return[$Failed]];
        <|
            "Objective" -> rawObj["Objective"],
            "Fidelity" -> rawObj["Fidelity"],
            "Penalty" -> rawObj["Penalty"],
            "Extra" -> KeyDrop[rawObj, {"Objective", "Fidelity", "Penalty"}]
        |>
    ];
    evalTrialFn = If[useObjTrial, evalObjScaled, evalScaled];

    appendCurvaturePair[sNew_List, yNew_List] := Module[{sCol, yCol, syCol, ysNewCol},
        If[Length[sMat] === 0,
            sMat = Developer`ToPackedArray[{sNew}];
            yMat = Developer`ToPackedArray[{yNew}];
            ssMat = Developer`ToPackedArray[{{sNew . sNew}}];
            syMat = Developer`ToPackedArray[{{sNew . yNew}}];
            yyMat = Developer`ToPackedArray[{{yNew . yNew}}];
            Return[]
        ];
        If[Length[sMat] >= mMem,
            sMat = Rest[sMat]; yMat = Rest[yMat];
            ssMat = ssMat[[2 ;;, 2 ;;]]; syMat = syMat[[2 ;;, 2 ;;]]; yyMat = yyMat[[2 ;;, 2 ;;]]
        ];
        sCol = sMat . sNew;
        yCol = yMat . yNew;
        syCol = sMat . yNew;
        ysNewCol = yMat . sNew;
        sMat = Developer`ToPackedArray[Append[sMat, sNew]];
        yMat = Developer`ToPackedArray[Append[yMat, yNew]];
        ssMat = Developer`ToPackedArray[ArrayFlatten[{{ssMat, List /@ sCol}, {{sCol}, {{sNew . sNew}}}}]];
        syMat = Developer`ToPackedArray[ArrayFlatten[{{syMat, List /@ syCol}, {{ysNewCol}, {{sNew . yNew}}}}]];
        yyMat = Developer`ToPackedArray[ArrayFlatten[{{yyMat, List /@ yCol}, {{yCol}, {{yNew . yNew}}}}]];
    ];

    th = Developer`ToPackedArray[projBox[x0Vec * scaleVec]];
    With[{res0 = evalScaled[th]},
        If[res0 === $Failed, Return[$Failed]];
        {fVal, fidVal, penVal, gVal, extra} =
            {res0["Objective"], res0["Fidelity"], res0["Penalty"], res0["Gradient"], res0["Extra"]}
    ];
    phiVal = -fVal;
    cNM = phiVal;
    bestFid = fidVal;
    gMin = Developer`ToPackedArray[-gVal];
    gMinNorm = Norm[gMin];
    pg = Developer`ToPackedArray[projBox[th + gVal] - th];
    pgNorm = Norm[pg];
    sigma = sigmaInit;
    objHist = If[TrueQ[recHist], {fVal}, {}];
    fidHist = If[TrueQ[recHist], {fidVal}, {}];
    gHist = If[TrueQ[recHist], {pgNorm}, {}];

    Do[
        If[NumericQ[targetFid] && TrueQ[evalCritFn[fidVal, extra, fidCrit] >= targetFid], status = "TargetFidelityReached"; Break[]];
        If[TrueQ[pgNorm <= gradTol], status = "GradientToleranceReached"; Break[]];
        If[AbsoluteTime[] - tStart >= tLimit, status = "TimeConstraintReached"; Break[]];

        phiRef = If[accMode === "ZhangHager", Max[cNM, phiVal], phiVal];
        tolScale = 1.*^-14 * Max[1., Abs[phiVal]];
        If[Length[sMat] > 0,
            syLast = syMat[[-1, -1]];
            yyLast = yyMat[[-1, -1]];
            delta = Clip[If[syLast > 1.*^-24, yyLast / syLast, 1.], {1.*^-8, 1.*^8}],
            delta = 1.
        ];

        accepted = False;
        Do[
            mu = Max[1.*^-14, sigma * Max[pgNorm, 1.*^-12]^gradPow];
            stepOut = CompactRegularizedLBFGSStep[sMat, yMat, gMin, delta, mu,
                If[Length[sMat] > 0, <|"SS" -> ssMat, "SY" -> syMat, "YY" -> yyMat|>, Automatic]];
            If[stepOut === $Failed,
                (* Trim corrupted/collinear history once before falling back to damped steepest *)
                If[Length[sMat] > 1,
                    sMat = Rest[sMat]; yMat = Rest[yMat];
                    ssMat = ssMat[[2 ;;, 2 ;;]]; syMat = syMat[[2 ;;, 2 ;;]]; yyMat = yyMat[[2 ;;, 2 ;;]];
                    Continue[]
                ];
                sMat = {}; yMat = {}; ssMat = {}; syMat = {}; yyMat = {}; delta = 1.;
                stepOut = CompactRegularizedLBFGSStep[{}, {}, gMin, delta, mu, Automatic]
            ];
            dVec = stepOut["Direction"];
            predDec = stepOut["PredDecrease"];
            If[anyConstr,
                thTrial = projBox[th + dVec];
                dVec = Developer`ToPackedArray[thTrial - th];
                predDec = -(gMin . dVec);
                If[!TrueQ[predDec > 1.*^-16 * pgNorm * Norm[dVec]],
                    rejectedSteps++;
                    sigma = Min[sigmaMax, gamma2 * sigma];
                    Continue[]
                ],
                thTrial = Developer`ToPackedArray[th + dVec]
            ];

            If[TrueQ[Norm[dVec] <= 1.*^-15 * Max[1., Norm[th]]],
                Break[]
            ];
            resTrial = evalTrialFn[thTrial];
            If[resTrial === $Failed || !NumericQ[resTrial["Objective"]],
                rejectedSteps++;
                sigma = Min[sigmaMax, gamma2 * sigma];
                Continue[]
            ];
            fTrial = resTrial["Objective"];
            phiTrial = -fTrial;
            ared = phiRef - phiTrial;
            rhoRatio = If[TrueQ[predDec > 1.*^-24], ared / predDec, -1.];

            If[(TrueQ[rhoRatio >= eta1] && TrueQ[ared > -Min[tolScale, 0.25 * eta1 * Max[0., predDec]]]) ||
               (TrueQ[phiTrial <= phiRef - eta1 * Max[0., predDec] + Min[tolScale, 0.25 * eta1 * Max[0., predDec]]]),
                gradRes = If[KeyExistsQ[resTrial, "Gradient"], resTrial, evalScaled[thTrial]];
                If[gradRes === $Failed, Return[$Failed, Module]];
                accepted = True;
                Break[],
                rejectedSteps++;
                sigma = Min[sigmaMax, gamma2 * sigma]
            ],
            {inIter, maxInner}
        ];

        If[!accepted,
            status = If[pgNorm <= Sqrt[gradTol], "GradientToleranceReached", "LineSearchTerminated"];
            Break[]
        ];

        acceptedSteps++;
        lastStep = Norm[thTrial - th];
        fOld = fVal;
        {fVal, fidVal, penVal, gVal, extra} =
            {gradRes["Objective"], gradRes["Fidelity"], gradRes["Penalty"], gradRes["Gradient"], gradRes["Extra"]};
        phiVal = -fVal;
        With[{rhoMono = If[TrueQ[predDec > 1.*^-24], (-fOld - phiVal) / predDec, -1.]},
            If[TrueQ[rhoMono >= eta2],
                sigma = Max[sigmaMin, gamma1 * sigma]
            ]
        ];
        If[accMode === "ZhangHager",
            With[{qOld = qNM},
                qNM = etaZH * qOld + 1.;
                cNM = (etaZH * qOld * cNM + phiVal) / qNM
            ],
            cNM = phiVal
        ];

        sVec = Developer`ToPackedArray[thTrial - th];
        yVec = Developer`ToPackedArray[(-gVal) - gMin];
        sy = sVec . yVec;
        If[!anyConstr && TrueQ[sy > 1.*^-12 * Norm[sVec] * Norm[yVec]],
            appendCurvaturePair[sVec, yVec],
            If[anyConstr && TrueQ[sy > 1.*^-12 * Norm[sVec] * Norm[yVec]],
                appendCurvaturePair[sVec, yVec],
                skippedCurv++
            ]
        ];

        th = thTrial;
        gMin = Developer`ToPackedArray[-gVal];
        gMinNorm = Norm[gMin];
        bestFid = Max[bestFid, fidVal];
        pg = Developer`ToPackedArray[projBox[th + gVal] - th];
        pgNorm = Norm[pg];
        nIter++;
        If[TrueQ[verbose] && (Mod[nIter, 10] === 0 || nIter === 1),
            Print["Iter ", nIter, " (RegularizedLBFGS): Fidelity=", NumberForm[fidVal, 6],
                " Objective=", NumberForm[fVal, 6], " sigma=", ScientificForm[sigma, 2],
                " |grad|=", ScientificForm[pgNorm, 3]]
        ];
        If[TrueQ[recHist],
            AppendTo[objHist, fVal];
            AppendTo[fidHist, fidVal];
            AppendTo[gHist, pgNorm]
        ];

        With[{critFidNew = evalCritFn[fidVal, extra, fidCrit]},
            stopNow = Which[
                TrueQ[Abs[fVal - fOld] <= funTol * Max[1., Abs[fOld]]] && TrueQ[pgNorm <= Sqrt[gradTol]],
                    If[NumericQ[targetFid] && TrueQ[critFidNew >= targetFid], "TargetFidelityReached", "FunctionToleranceReached"],
                TrueQ[Norm[sVec] <= stepTol * Max[1., Norm[th]]],
                    If[NumericQ[targetFid] && TrueQ[critFidNew >= targetFid], "TargetFidelityReached", "StepToleranceReached"],
                NumericQ[targetFid] && TrueQ[critFidNew >= targetFid],
                    "TargetFidelityReached",
                TrueQ[pgNorm <= gradTol],
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
            monCfg, nIter, evals, fVal, fidVal, bestFid, pgNorm, lastStep,
            AbsoluteTime[] - tStart, If[stopNow === None, "Running", stopNow],
            extra, stageTag, lastMonTime, lastMonIter
        ];
        If[stopNow =!= None, status = stopNow; Break[]],
        {maxIter}
    ];
    If[lastMonIter < nIter || (lastMonIter < 0 && nIter === 0),
        {lastMonTime, lastMonIter} = emitMonFn[
            monCfg, nIter, evals, fVal, fidVal, bestFid, pgNorm, lastStep,
            AbsoluteTime[] - tStart, status, extra, stageTag, lastMonTime, lastMonIter
        ]
    ];

    With[{elapsed = AbsoluteTime[] - tStart},
        <|
            "ParametersVec" -> Developer`ToPackedArray[th / scaleVec],
            "Objective" -> fVal,
            "Fidelity" -> fidVal,
            "Penalty" -> penVal,
            "GradientNorm" -> pgNorm,
            "ScaledGradient" -> gVal,
            "Optimizer" -> "RegularizedLBFGS",
            "Iterations" -> nIter,
            "Evaluations" -> evals,
            "ObjectiveEvaluations" -> objEvals,
            "GradientEvaluations" -> gradEvals,
            "HVPEvaluations" -> 0,
            "HessianEvaluations" -> 0,
            "AcceptedSteps" -> acceptedSteps,
            "RejectedSteps" -> rejectedSteps,
            "SkippedCurvaturePairs" -> skippedCurv,
            "CompactLBFGSFinalSigma" -> sigma,
            "ExperimentalConstrainedQ" -> anyConstr,
            "LineSearchMethod" -> ("CompactRegularized_" <> accMode),
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

(* Register line searches and optimizers into GRAPE.m $OCSolverRegistry *)
If[!KeyExistsQ[$OCSolverRegistry, "LineSearches"] || !AssociationQ[$OCSolverRegistry["LineSearches"]],
    $OCSolverRegistry["LineSearches"] = <||>
];
If[!KeyExistsQ[$OCSolverRegistry, "Optimizers"] || !AssociationQ[$OCSolverRegistry["Optimizers"]],
    $OCSolverRegistry["Optimizers"] = <||>
];

$OCSolverRegistry["LineSearches"] = Join[$OCSolverRegistry["LineSearches"], <|
    "HagerZhang" -> HagerZhangLineSearch,
    "ZhangHagerNonmonotone" -> ZhangHagerNonmonotoneLineSearch,
    "GrippoNonmonotone" -> GrippoNonmonotoneLineSearch
|>];

$OCSolverRegistry["Optimizers"] = Join[$OCSolverRegistry["Optimizers"], <|
    "RegularizedLBFGS" -> RunRegularizedLBFGS
|>];

$OCSolverRegistry["CompactRegularizedLBFGSStep"] = CompactRegularizedLBFGSStep;

End[];
EndPackage[];
