(* ::Package:: *)

(* SpinDynamica Optimal Control — Pure Numerical GRAPE Kernels & Control Parameterizations *)

BeginPackage["OptimalControl`Kernels`", {"Global`", "SpinDynamica`"}];

ToPackedComplexMatrix::usage = "ToPackedComplexMatrix[mat] converts a numeric matrix or SparseArray into a dense packed Complex128 matrix.";
ToPackedComplexVector::usage = "ToPackedComplexVector[vec] converts a numeric vector or SparseArray into a dense packed Complex128 vector.";
ToPackedRealMatrix::usage = "ToPackedRealMatrix[mat] converts a numeric matrix into a dense packed Real64 matrix.";
HermitianSymmetrize::usage = "HermitianSymmetrize[mat] returns the Hermitian part 0.5*(mat + ConjugateTranspose[mat]) as a packed Complex128 matrix.";
BuildControlBasis::usage = "BuildControlBasis[spec, dtVec, window] constructs a control parameterization Association containing the basis matrix mapping reduced parameters to piecewise-constant slice amplitudes.";
EvaluateHilbertSpectral::usage = "EvaluateHilbertSpectral[h0, hCtrlRows, hFlatT, dtVec, uMat, objFn, returnProps] evaluates closed-system unitary propagation and exact analytical GRAPE gradients via Hermitian eigendecomposition and Fréchet spectral division.";
EvaluateHilbertBlockExp::usage = "EvaluateHilbertBlockExp[h0, hCtrls, dtVec, uMat, objFn] evaluates closed-system unitary propagation and exact analytical gradients via block-matrix exponentials.";
EvaluateHilbertFirstOrder::usage = "EvaluateHilbertFirstOrder[h0, hCtrls, dtVec, uMat, objFn] evaluates closed-system unitary propagation and first-order approximate gradients for benchmarking against exact Fréchet gradients.";
EvaluateHilbertStateVector::usage = "EvaluateHilbertStateVector[h0, hCtrlRows, hFlatT, dtVec, uMat, psiMat0, tarMat, weights, phaseIns] evaluates multi-ket state-vector forward/backward GRAPE propagation without assembling full slice propagators.";
EvaluateLiouvilleExact::usage = "EvaluateLiouvilleExact[l0, lCtrls, dtVec, uMat, rhoMat0, obsMat, weights, objModes, scales] evaluates open-system or superoperator forward/backward GRAPE propagation and exact augmented block-exponential gradients.";
PropagateControls::usage = "PropagateControls[compiled, uMat] performs forward-only numerical propagation under slice control matrix uMat and returns slice and final propagators or states.";

Begin["`Private`"];

ToPackedComplexMatrix[mat_] :=
    Developer`ToPackedArray[N[Normal[mat]] + 0. I];

ToPackedComplexVector[vec_] :=
    Developer`ToPackedArray[N[Normal[vec]] + 0. I];

ToPackedRealMatrix[mat_] :=
    Developer`ToPackedArray[N[Re[Normal[mat]]]];

HermitianSymmetrize[mat_] :=
    Developer`ToPackedArray[0.5 (mat + ConjugateTranspose[mat])];

buildWindowVector[windowSpec_, tMid_List, totalDur_?NumericQ] := Module[{nt, tau, vals},
    nt = Length[tMid];
    If[windowSpec === Automatic || windowSpec === None || windowSpec === False,
        Return[ConstantArray[1., nt]]
    ];
    tau = tMid / totalDur;
    vals = Which[
        windowSpec === "Hann", 0.5 (1. - Cos[2. Pi tau]),
        windowSpec === "Hamming", 0.54 - 0.46 Cos[2. Pi tau],
        windowSpec === "Blackman", 0.42 - 0.5 Cos[2. Pi tau] + 0.08 Cos[4. Pi tau],
        ListQ[windowSpec] && Length[windowSpec] === nt && VectorQ[windowSpec, NumericQ], N[windowSpec],
        True, Quiet[Check[N[Table[windowSpec[u], {u, tau}]], $Failed]]
    ];
    If[!ListQ[vals] || Length[vals] =!= nt || !VectorQ[vals, NumericQ],
        $Failed,
        Developer`ToPackedArray[N[Re[vals]]]
    ]
];

BuildControlBasis[paramSpec_, dtVec_List, windowSpec_: Automatic] := Module[
    {nt, tEdges, totalDur, tMid, winVec, hasWin, rawB, bMat, nH, nOrd, nB, deg, knots, specName},
    nt = Length[dtVec];
    If[nt < 1 || !VectorQ[dtVec, NumericQ], Return[$Failed]];
    tEdges = Prepend[Accumulate[N[dtVec]], 0.];
    totalDur = Last[tEdges];
    If[!TrueQ[totalDur > 0], Return[$Failed]];
    tMid = 0.5 (Most[tEdges] + Rest[tEdges]);
    winVec = buildWindowVector[windowSpec, tMid, totalDur];
    If[winVec === $Failed, Return[$Failed]];
    hasWin = !MatchQ[windowSpec, Automatic | None | False];
    {specName, rawB} = Which[
        paramSpec === Automatic || paramSpec === "PiecewiseConstant",
            {"PiecewiseConstant", N[IdentityMatrix[nt]]},
        MatchQ[paramSpec, {"Fourier", _Integer?Positive}],
            nH = paramSpec[[2]];
            {"Fourier",
                Transpose[Join[
                    {ConstantArray[1., nt]},
                    Flatten[Table[
                        {Cos[2. Pi harmonic tMid / totalDur], Sin[2. Pi harmonic tMid / totalDur]},
                        {harmonic, nH}], 1]
                ]]},
        MatchQ[paramSpec, {"FourierSine", _Integer?Positive}],
            nH = paramSpec[[2]];
            {"FourierSine",
                Transpose[Table[Sin[Pi harmonic tMid / totalDur], {harmonic, nH}]]},
        MatchQ[paramSpec, {"Chebyshev", _Integer?NonNegative}],
            nOrd = paramSpec[[2]];
            {"Chebyshev",
                Transpose[Table[ChebyshevT[ord, 2. tMid / totalDur - 1.], {ord, 0, nOrd}]]},
        MatchQ[paramSpec, {"BSpline", _Integer?Positive} | {"BSpline", _Integer?Positive, _Integer?NonNegative}],
            nB = paramSpec[[2]];
            deg = If[Length[paramSpec] >= 3, Min[paramSpec[[3]], nB - 1], Min[3, nB - 1]];
            If[nB === 1,
                {"BSpline", ConstantArray[1., {nt, 1}]},
                knots = Developer`ToPackedArray[N[Join[
                    ConstantArray[0., deg + 1],
                    If[nB > deg + 1, Range[1, nB - deg - 1] / (nB - deg), {}],
                    ConstantArray[1., deg + 1]
                ]]];
                {"BSpline",
                    Table[N[BSplineBasis[{deg, knots}, basisIdx - 1, tMid[[j]] / totalDur]],
                        {j, nt}, {basisIdx, nB}]}
            ],
        MatchQ[paramSpec, {"BasisMatrix", _?MatrixQ}] && Dimensions[paramSpec[[2]]][[1]] === nt,
            {"BasisMatrix", N[Normal[paramSpec[[2]]]]},
        True,
            {$Failed, $Failed}
    ];
    If[rawB === $Failed || !MatrixQ[rawB, NumericQ], Return[$Failed]];
    bMat = Developer`ToPackedArray[N[Re[If[hasWin, winVec * rawB, rawB]]]];
    <|
        "Type" -> specName,
        "BasisMatrix" -> bMat,
        "BasisTranspose" -> Developer`ToPackedArray[Transpose[bMat]],
        "NumSlices" -> nt,
        "NumParams" -> Dimensions[bMat][[2]],
        "IsIdentityQ" -> (specName === "PiecewiseConstant" && !hasWin),
        "Window" -> winVec
    |>
];

(* Helper routines for intermediate-time running costs and suffix-sum timing gradients *)
prepareRunningRuntime[rcSpec_Association, dtCoarse_List, t0_?NumericQ] := Module[
    {nt = Length[dtCoarse], tEnds, twSpecs = rcSpec["TWSpecs"], wMat, raw},
    tEnds = Developer`ToPackedArray[N[t0] + Accumulate[N[dtCoarse]]];
    wMat = Table[
        If[TrueQ[twSpecs[[i, "StaticQ"]]],
            twSpecs[[i, "Weights"]],
            raw = (twSpecs[[i, "EvalRaw"]])[dtCoarse, t0];
            If[!VectorQ[raw, NumericQ],
                ConstantArray[0., nt],
                Developer`ToPackedArray[If[Total[raw] > 0., raw / Total[raw], ConstantArray[0., nt]]]
            ]
        ],
        {i, Length[twSpecs]}
    ];
    <|"WeightsMat" -> wMat, "TEnds" -> tEnds|>
];

evalRunningBoundary[evalList_List, wMat_, stateJ_, tJ_?NumericQ, j_Integer, zeroCot_] := Module[
    {vSum = 0., cSum = zeroCot, idx, pair},
    Do[
        pair = (evalList[[idx]])[stateJ, tJ, j, wMat[[idx, j]]];
        If[!ListQ[pair], Return[$Failed, Module]];
        vSum += N[Re[pair[[1]]]];
        cSum = cSum + pair[[2]],
        {idx, Length[evalList]}
    ];
    {vSum, ToPackedComplexMatrix[cSum]}
];

evalRunningTimeGradient[rcSpec_Association, evalList_List, statesFine_, nStages_Integer, dtCoarse_List, t0_?NumericQ] := Module[
    {nt = Length[dtCoarse], coarseStates, evalRunTotalDt, dtPlus, dtMinus, hDt, jj},
    If[!TrueQ[rcSpec["HasTimeDependence"]], Return[ConstantArray[0., nt]]];
    coarseStates = statesFine[[nStages + 1 ;; ;; nStages]];
    evalRunTotalDt[dtIn_List] := Module[{rt, wM, tE, sumV = 0., idx, jIdx, pair},
        rt = prepareRunningRuntime[rcSpec, dtIn, t0];
        wM = rt["WeightsMat"];
        tE = rt["TEnds"];
        Do[
            Do[
                pair = (evalList[[idx]])[coarseStates[[jIdx]], tE[[jIdx]], jIdx, wM[[idx, jIdx]]];
                If[ListQ[pair], sumV += N[Re[pair[[1]]]]],
                {idx, Length[evalList]}
            ],
            {jIdx, nt}
        ];
        sumV
    ];
    Developer`ToPackedArray[Table[
        hDt = Min[0.1 * dtCoarse[[jj]], Max[10^-12, 10^-6 * dtCoarse[[jj]], 10^-8 * (Abs[t0] + Total[dtCoarse])]];
        dtPlus = dtCoarse; dtPlus[[jj]] += hDt;
        dtMinus = dtCoarse; dtMinus[[jj]] = Max[10^-14, dtMinus[[jj]] - hDt];
        (evalRunTotalDt[dtPlus] - evalRunTotalDt[dtMinus]) / (dtPlus[[jj]] - dtMinus[[jj]]),
        {jj, nt}
    ]]
];

foldTimingGradient[gExp_List, sigmaSelf_List, sigmaShift_List, sigmaRun_List] := Module[
    {nt = Length[gExp], suf = 0., outDt = ConstantArray[0., Length[gExp]], jj},
    Do[
        outDt[[jj]] = gExp[[jj]] + sigmaSelf[[jj]] + sigmaRun[[jj]] + suf;
        suf += sigmaShift[[jj]],
        {jj, nt, 1, -1}
    ];
    Developer`ToPackedArray[outDt]
];

(* Singularity-free spectral 2nd divided differences G^(2) and tensor contraction helpers *)
compileG2Tensor[] := Quiet[Check[
    Compile[{{lam, _Real, 1}, {dt, _Real}, {phi1, _Complex, 2}},
        Module[{d = Length[lam], z, g2, a, b, c, za, zb, zc, dab, dbc, dac, sMax,
                zbar, da, db, dc, e2, e3, e2sq, val},
            z = (-I * dt) * lam;
            g2 = Table[0. + 0. I, {d}, {d}, {d}];
            Do[
                za = z[[a]]; zb = z[[b]]; zc = z[[c]];
                dab = Abs[za - zb]; dbc = Abs[zb - zc]; dac = Abs[za - zc];
                sMax = Max[dab, dbc, dac];
                val = If[sMax <= 0.15,
                    zbar = (za + zb + zc) / 3.;
                    da = za - zbar; db = zb - zbar; dc = zc - zbar;
                    e2 = da * db + db * dc + dc * da;
                    e3 = da * db * dc;
                    e2sq = e2 * e2;
                    Exp[zbar] * (0.5 - e2 / 24. + e3 / 120. + e2sq / 720. - (2. * e2 * e3) / 5040. +
                        (e3 * e3 - e2 * e2sq) / 40320. + (3. * e2sq * e3) / 362880. +
                        (e2sq * e2sq - 3. * e2 * e3 * e3) / 3628800.),
                    Which[
                        dac >= dab && dac >= dbc, (phi1[[a, b]] - phi1[[b, c]]) / (za - zc),
                        dab >= dbc, (phi1[[a, c]] - phi1[[c, b]]) / (za - zb),
                        True, (phi1[[b, a]] - phi1[[a, c]]) / (zb - zc)
                    ]
                ];
                g2[[a, b, c]] = val;
                If[a != b || b != c,
                    g2[[a, c, b]] = val;
                    g2[[b, a, c]] = val;
                    g2[[b, c, a]] = val;
                    g2[[c, a, b]] = val;
                    g2[[c, b, a]] = val
                ],
                {a, d}, {b, a, d}, {c, b, d}
            ];
            g2
        ],
        CompilationTarget -> "WVM", RuntimeOptions -> "Speed"
    ],
    $Failed
]];

$compiledG2TensorFn = compileG2Tensor[];

buildG2Tensor[lam_List, dt_?NumericQ, phi1Mat_] := Module[
    {d, z, g2, a, b, c, za, zb, zc, dab, dbc, dac, sMax, zbar, da, db, dc, e2, e3, e2sq, val},
    If[$compiledG2TensorFn =!= $Failed,
        Return[ToPackedComplexTensor[$compiledG2TensorFn[Developer`ToPackedArray[N[Re[lam]]], N[dt], ToPackedComplexMatrix[phi1Mat]]]]
    ];
    d = Length[lam];
    z = (-I * N[dt]) * N[Re[lam]];
    g2 = ConstantArray[0. + 0. I, {d, d, d}];
    Do[
        za = z[[a]]; zb = z[[b]]; zc = z[[c]];
        dab = Abs[za - zb]; dbc = Abs[zb - zc]; dac = Abs[za - zc];
        sMax = Max[dab, dbc, dac];
        val = If[sMax <= 0.15,
            zbar = (za + zb + zc) / 3.;
            da = za - zbar; db = zb - zbar; dc = zc - zbar;
            e2 = da * db + db * dc + dc * da;
            e3 = da * db * dc;
            e2sq = e2 * e2;
            Exp[zbar] * (0.5 - e2 / 24. + e3 / 120. + e2sq / 720. - (2. * e2 * e3) / 5040. +
                (e3 * e3 - e2 * e2sq) / 40320. + (3. * e2sq * e3) / 362880. +
                (e2sq * e2sq - 3. * e2 * e3 * e3) / 3628800.),
            Which[
                dac >= dab && dac >= dbc, (phi1Mat[[a, b]] - phi1Mat[[b, c]]) / (za - zc),
                dab >= dbc, (phi1Mat[[a, c]] - phi1Mat[[c, b]]) / (za - zb),
                True, (phi1Mat[[b, a]] - phi1Mat[[a, c]]) / (zb - zc)
            ]
        ];
        g2[[a, b, c]] = val;
        If[a != b || b != c,
            g2[[a, c, b]] = val;
            g2[[b, a, c]] = val;
            g2[[b, c, a]] = val;
            g2[[c, a, b]] = val;
            g2[[c, b, a]] = val
        ],
        {a, d}, {b, a, d}, {c, b, d}
    ];
    Developer`ToPackedArray[g2]
];

ToPackedComplexTensor[t_] := Developer`ToPackedArray[N[t] + 0. I];

compileG2Contract[] := Quiet[Check[
    Compile[{{aMat, _Complex, 2}, {vMat, _Complex, 2}, {g2, _Complex, 3}},
        Module[{d = Length[aMat], res, r, s, p, acc},
            res = Table[0. + 0. I, {d}, {d}];
            Do[
                acc = 0. + 0. I;
                Do[
                    acc += g2[[r, s, p]] * (aMat[[r, p]] * vMat[[p, s]] + vMat[[r, p]] * aMat[[p, s]]),
                    {p, d}
                ];
                res[[r, s]] = acc,
                {r, d}, {s, d}
            ];
            res
        ],
        CompilationTarget -> "WVM", RuntimeOptions -> "Speed"
    ],
    $Failed
]];

$compiledG2ContractFn = compileG2Contract[];

contractG2[aMat_, vMat_, g2_] := If[$compiledG2ContractFn =!= $Failed,
    $compiledG2ContractFn[aMat, vMat, g2],
    Module[{d = Length[aMat]},
        Developer`ToPackedArray[Table[
            ((aMat[[All, All]] * g2[[All, s, All]]) . vMat[[All, s]]) +
            ((vMat[[All, All]] * g2[[All, s, All]]) . aMat[[All, s]]),
            {s, d}
        ]] // Transpose
    ]
];

EnsureSpectralG2Cache[evalState_Association] := Module[{nFine, evals, dtFine, phi1Mats, g2Mats},
    If[KeyExistsQ[evalState, "G2Mats"] && ListQ[evalState["G2Mats"]],
        Return[evalState]
    ];
    dtFine = evalState["DtFine"];
    evals = evalState["Evals"];
    phi1Mats = evalState["Phi1Mats"];
    nFine = Length[dtFine];
    g2Mats = Table[buildG2Tensor[evals[[j]], dtFine[[j]], phi1Mats[[j]]], {j, nFine}];
    Join[evalState, <|"G2Mats" -> g2Mats|>]
];

(* Core closed-system propagator & exact analytical Fréchet gradient using Hermitian Eigensystem *)
EvaluateHilbertSpectral[h0_, hCtrlRows_, hFlatT_, dtVec_List, uMat_, objCotangentFn_, returnProps_: False,
                        nStages_Integer: 1, runningCostSpec_: None, timingGradSpec_: None, returnCache_: False] := Module[
    {nt, nFine, mCtrl, dim, isTD, uFine, hAll, hs, vs, vsH, evals, pHalfs, us, fs, j,
     hj, ev, vec, vecH, pHalf, phases, uj, uFinal, idMat, zeroDim,
     phi, cDag, extra, xMat, aMat, halfDt, lam, lamRow, phiMat, qMat, wMat, gradFine, gradOut,
     coarseUs, coarseFs,
     phi1Mats, aTilMats, xMats, g2Mats, aTil, cacheAssoc = None, needCache = TrueQ[returnCache],
     hasRC, rcEvalList, rcDtCoarse, rcT0, rcRT, rcWMat, rcTEnds, rcTotal = 0., rcPair, coarseJ,
     hasTimingGrad, gExp, sigmaSelf, sigmaShift, sigmaRun, alphas, isTDCont, dSelf, dShift, alphaJ, gradDtOut},
    nt = Length[uMat];
    nFine = Length[dtVec];
    mCtrl = Dimensions[uMat][[2]];
    isTD = ArrayDepth[h0] === 3;
    dim = If[isTD, Dimensions[h0][[2]], Length[h0]];
    uFine = If[nStages === 1, uMat, Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ uMat, 1]]];
    hAll = ArrayReshape[uFine . hCtrlRows, {nFine, dim, dim}];
    hasRC = AssociationQ[runningCostSpec];
    hasTimingGrad = AssociationQ[timingGradSpec];
    If[hasTimingGrad, hs = Table[Null, {nFine}]];
    vs = Table[Null, {nFine}];
    vsH = Table[Null, {nFine}];
    evals = Table[Null, {nFine}];
    pHalfs = Table[Null, {nFine}];
    us = Table[Null, {nFine}];
    fs = Table[Null, {nFine + 1}];
    If[needCache,
        phi1Mats = Table[Null, {nFine}];
        aTilMats = Table[Null, {nFine}];
        xMats = Table[Null, {nFine}];
        g2Mats = Table[Null, {nFine}]
    ];
    idMat = Developer`ToPackedArray[IdentityMatrix[dim, WorkingPrecision -> MachinePrecision] + 0. I];
    fs[[1]] = idMat;
    Do[
        hj = If[isTD, h0[[j]], h0] + hAll[[j]];
        hj = 0.5 (hj + ConjugateTranspose[hj]);
        If[hasTimingGrad, hs[[j]] = hj];
        {ev, vec} = Eigensystem[hj];
        ev = Re[ev];
        vec = Transpose[vec];
        vecH = ConjugateTranspose[vec];
        pHalf = Exp[(-0.5 I dtVec[[j]]) ev];
        phases = pHalf * pHalf;
        uj = vec . (phases * vecH);
        vs[[j]] = vec;
        vsH[[j]] = vecH;
        evals[[j]] = ev;
        pHalfs[[j]] = pHalf;
        us[[j]] = uj;
        fs[[j + 1]] = uj . fs[[j]],
        {j, nFine}
    ];
    uFinal = fs[[nFine + 1]];
    {phi, cDag, extra} = objCotangentFn[uFinal];
    xMat = ToPackedComplexMatrix[cDag];
    If[hasRC,
        zeroDim = ConstantArray[0. + 0. I, {dim, dim}];
        rcEvalList = runningCostSpec["HilbertTermEvaluators"];
        rcDtCoarse = Which[
            hasTimingGrad, timingGradSpec["DtCoarse"],
            nStages === 1, dtVec,
            True, Total[Partition[dtVec, nStages], {2}]
        ];
        rcT0 = Which[
            hasTimingGrad, Lookup[timingGradSpec, "InitialTimePoint", 0.],
            hasRC, Lookup[runningCostSpec, "InitialTimePoint", 0.],
            True, 0.
        ];
        rcRT = prepareRunningRuntime[runningCostSpec, rcDtCoarse, rcT0];
        rcWMat = rcRT["WeightsMat"];
        rcTEnds = rcRT["TEnds"]
    ];
    If[hasTimingGrad,
        gExp = ConstantArray[0., nt];
        sigmaSelf = ConstantArray[0., nt];
        sigmaShift = ConstantArray[0., nt];
        alphas = timingGradSpec["Alphas"];
        isTDCont = TrueQ[timingGradSpec["IsTDContinuous"]];
        dSelf = Lookup[timingGradSpec, "DSelf", None];
        dShift = Lookup[timingGradSpec, "DShift", None];
        sigmaRun = If[hasRC,
            evalRunningTimeGradient[runningCostSpec, rcEvalList, fs, nStages, rcDtCoarse, rcT0],
            ConstantArray[0., nt]
        ]
    ];
    gradFine = ConstantArray[0., {nFine, mCtrl}];
    Do[
        If[hasRC && Mod[j, nStages] === 0,
            coarseJ = Quotient[j, nStages];
            rcPair = evalRunningBoundary[rcEvalList, rcWMat, fs[[j + 1]], rcTEnds[[coarseJ]], coarseJ, zeroDim];
            If[rcPair === $Failed, Return[$Failed, Module]];
            rcTotal += rcPair[[1]];
            xMat = xMat + rcPair[[2]]
        ];
        If[needCache, xMats[[j]] = xMat];
        aMat = fs[[j]] . xMat;
        vec = vs[[j]];
        vecH = vsH[[j]];
        lam = evals[[j]];
        pHalf = pHalfs[[j]];
        halfDt = 0.5 dtVec[[j]];
        lamRow = ConstantArray[lam, dim];
        phiMat = (Transpose[{pHalf}] . {pHalf}) * Sinc[halfDt (Transpose[lamRow] - lamRow)];
        aTil = vecH . aMat . vec;
        qMat = aTil * phiMat;
        wMat = vec . qMat . vecH;
        gradFine[[j]] = dtVec[[j]] * Im[hFlatT . Flatten[wMat]];
        If[needCache,
            phi1Mats[[j]] = phiMat;
            aTilMats[[j]] = aTil;
            g2Mats[[j]] = buildG2Tensor[lam, dtVec[[j]], phiMat]
        ];
        If[hasTimingGrad,
            coarseJ = Ceiling[j / nStages];
            alphaJ = alphas[[j]];
            gExp[[coarseJ]] += alphaJ * Im[Total[Flatten[Transpose[fs[[j + 1]] . xMat] * hs[[j]]]]];
            If[isTDCont,
                sigmaSelf[[coarseJ]] += dtVec[[j]] * Im[Total[Flatten[Transpose[dSelf[[j]]] * wMat]]];
                sigmaShift[[coarseJ]] += dtVec[[j]] * Im[Total[Flatten[Transpose[dShift[[j]]] * wMat]]]
            ]
        ];
        xMat = xMat . us[[j]],
        {j, nFine, 1, -1}
    ];
    gradOut = If[nStages === 1,
        Developer`ToPackedArray[gradFine],
        Developer`ToPackedArray[Total[Partition[gradFine, nStages], {2}]]
    ];
    gradDtOut = If[hasTimingGrad, foldTimingGradient[gExp, sigmaSelf, sigmaShift, sigmaRun], None];
    If[TrueQ[returnProps],
        If[nStages === 1,
            coarseUs = us; coarseFs = fs,
            coarseUs = Table[Fold[Dot, idMat, Reverse[us[[(j - 1) nStages + 1 ;; j nStages]]]], {j, nt}];
            coarseFs = fs[[1 ;; nFine + 1 ;; nStages]]
        ],
        coarseUs = None; coarseFs = None
    ];
    If[needCache,
        cacheAssoc = <|
            "Kernel" -> "Spectral",
            "Dim" -> dim,
            "NumSlices" -> nt,
            "NumFine" -> nFine,
            "NumStages" -> nStages,
            "NumControls" -> mCtrl,
            "DtFine" -> Developer`ToPackedArray[N[dtVec]],
            "Vs" -> vs,
            "VsH" -> vsH,
            "Evals" -> evals,
            "PHalfs" -> pHalfs,
            "Phi1Mats" -> phi1Mats,
            "ATilMats" -> aTilMats,
            "G2Mats" -> g2Mats,
            "Us" -> us,
            "Fs" -> fs,
            "XMats" -> xMats,
            "FinalPropagator" -> uFinal,
            "CDag" -> ToPackedComplexMatrix[cDag],
            "ControlRows" -> hCtrlRows,
            "ControlFlatT" -> hFlatT
        |>
    ];
    <|
        "Fidelity" -> phi,
        "TerminalFidelity" -> phi,
        "RunningCost" -> rcTotal,
        "Gradient" -> gradOut,
        "DurationGradient" -> gradDtOut,
        "FinalPropagator" -> uFinal,
        "Propagators" -> coarseUs,
        "States" -> coarseFs,
        "Extra" -> extra,
        If[needCache, "EvalState" -> cacheAssoc, Nothing]
    |>
];

(* Reference exact closed-system gradient using 2D x 2D block MatrixExp *)
EvaluateHilbertBlockExp[h0_, hCtrls_List, dtVec_List, uMat_, objCotangentFn_, returnProps_: False,
                        nStages_Integer: 1, runningCostSpec_: None, timingGradSpec_: None] := Module[
    {nt, nFine, mCtrl, dim, isTD, uFine, zMat, zeroDim, hs, us, fs, j, k, hj, uj, uFinal, idMat,
     phi, cDag, extra, xMat, aMat, dtj, blk, expBlk, dUj, gradFine, gradOut, coarseUs, coarseFs,
     hasRC, rcEvalList, rcDtCoarse, rcT0, rcRT, rcWMat, rcTEnds, rcTotal = 0., rcPair, coarseJ,
     hasTimingGrad, gExp, sigmaSelf, sigmaShift, sigmaRun, alphas, isTDCont, dSelf, dShift, alphaJ, gradDtOut},
    nt = Length[uMat];
    nFine = Length[dtVec];
    mCtrl = Length[hCtrls];
    isTD = ArrayDepth[h0] === 3;
    dim = If[isTD, Dimensions[h0][[2]], Length[h0]];
    uFine = If[nStages === 1, uMat, Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ uMat, 1]]];
    zMat = ConstantArray[0. + 0. I, {dim, dim}];
    hasRC = AssociationQ[runningCostSpec];
    hasTimingGrad = AssociationQ[timingGradSpec];
    hs = Table[Null, {nFine}];
    us = Table[Null, {nFine}];
    fs = Table[Null, {nFine + 1}];
    idMat = Developer`ToPackedArray[IdentityMatrix[dim, WorkingPrecision -> MachinePrecision] + 0. I];
    fs[[1]] = idMat;
    Do[
        hj = 0.5 ((# + ConjugateTranspose[#])) &[If[isTD, h0[[j]], h0] + Sum[uFine[[j, k]] hCtrls[[k]], {k, mCtrl}]];
        uj = MatrixExp[(-I dtVec[[j]]) hj];
        hs[[j]] = hj;
        us[[j]] = uj;
        fs[[j + 1]] = uj . fs[[j]],
        {j, nFine}
    ];
    uFinal = fs[[nFine + 1]];
    {phi, cDag, extra} = objCotangentFn[uFinal];
    xMat = ToPackedComplexMatrix[cDag];
    If[hasRC,
        zeroDim = zMat;
        rcEvalList = runningCostSpec["HilbertTermEvaluators"];
        rcDtCoarse = Which[
            hasTimingGrad, timingGradSpec["DtCoarse"],
            nStages === 1, dtVec,
            True, Total[Partition[dtVec, nStages], {2}]
        ];
        rcT0 = Which[
            hasTimingGrad, Lookup[timingGradSpec, "InitialTimePoint", 0.],
            hasRC, Lookup[runningCostSpec, "InitialTimePoint", 0.],
            True, 0.
        ];
        rcRT = prepareRunningRuntime[runningCostSpec, rcDtCoarse, rcT0];
        rcWMat = rcRT["WeightsMat"];
        rcTEnds = rcRT["TEnds"]
    ];
    If[hasTimingGrad,
        gExp = ConstantArray[0., nt];
        sigmaSelf = ConstantArray[0., nt];
        sigmaShift = ConstantArray[0., nt];
        alphas = timingGradSpec["Alphas"];
        isTDCont = TrueQ[timingGradSpec["IsTDContinuous"]];
        dSelf = Lookup[timingGradSpec, "DSelf", None];
        dShift = Lookup[timingGradSpec, "DShift", None];
        sigmaRun = If[hasRC,
            evalRunningTimeGradient[runningCostSpec, rcEvalList, fs, nStages, rcDtCoarse, rcT0],
            ConstantArray[0., nt]
        ]
    ];
    gradFine = ConstantArray[0., {nFine, mCtrl}];
    Do[
        If[hasRC && Mod[j, nStages] === 0,
            coarseJ = Quotient[j, nStages];
            rcPair = evalRunningBoundary[rcEvalList, rcWMat, fs[[j + 1]], rcTEnds[[coarseJ]], coarseJ, zeroDim];
            If[rcPair === $Failed, Return[$Failed, Module]];
            rcTotal += rcPair[[1]];
            xMat = xMat + rcPair[[2]]
        ];
        aMat = fs[[j]] . xMat;
        hj = hs[[j]];
        dtj = dtVec[[j]];
        Do[
            blk = Developer`ToPackedArray[ArrayFlatten[{{(-I dtj) hj, (-I dtj) hCtrls[[k]]}, {zMat, (-I dtj) hj}}]];
            expBlk = MatrixExp[blk];
            dUj = expBlk[[1 ;; dim, dim + 1 ;; 2 dim]];
            gradFine[[j, k]] = Re[Tr[aMat . dUj]],
            {k, mCtrl}
        ];
        If[hasTimingGrad,
            coarseJ = Ceiling[j / nStages];
            alphaJ = alphas[[j]];
            gExp[[coarseJ]] += alphaJ * Im[Total[Flatten[Transpose[fs[[j + 1]] . xMat] * hj]]];
            If[isTDCont,
                blk = Developer`ToPackedArray[ArrayFlatten[{{(-I dtj) hj, (-I dtj) dSelf[[j]]}, {zMat, (-I dtj) hj}}]];
                dUj = (MatrixExp[blk])[[1 ;; dim, dim + 1 ;; 2 dim]];
                sigmaSelf[[coarseJ]] += Re[Tr[aMat . dUj]];
                blk = Developer`ToPackedArray[ArrayFlatten[{{(-I dtj) hj, (-I dtj) dShift[[j]]}, {zMat, (-I dtj) hj}}]];
                dUj = (MatrixExp[blk])[[1 ;; dim, dim + 1 ;; 2 dim]];
                sigmaShift[[coarseJ]] += Re[Tr[aMat . dUj]]
            ]
        ];
        xMat = xMat . us[[j]],
        {j, nFine, 1, -1}
    ];
    gradOut = If[nStages === 1,
        Developer`ToPackedArray[gradFine],
        Developer`ToPackedArray[Total[Partition[gradFine, nStages], {2}]]
    ];
    gradDtOut = If[hasTimingGrad, foldTimingGradient[gExp, sigmaSelf, sigmaShift, sigmaRun], None];
    If[TrueQ[returnProps],
        If[nStages === 1,
            coarseUs = us; coarseFs = fs,
            coarseUs = Table[Fold[Dot, idMat, Reverse[us[[(j - 1) nStages + 1 ;; j nStages]]]], {j, nt}];
            coarseFs = fs[[1 ;; nFine + 1 ;; nStages]]
        ],
        coarseUs = None; coarseFs = None
    ];
    <|
        "Fidelity" -> phi,
        "TerminalFidelity" -> phi,
        "RunningCost" -> rcTotal,
        "Gradient" -> gradOut,
        "DurationGradient" -> gradDtOut,
        "FinalPropagator" -> uFinal,
        "Propagators" -> coarseUs,
        "States" -> coarseFs,
        "Extra" -> extra
    |>
];

(* First-order approximate commutator comparator *)
EvaluateHilbertFirstOrder[h0_, hCtrls_List, dtVec_List, uMat_, objCotangentFn_, returnProps_: False,
                          nStages_Integer: 1, runningCostSpec_: None, timingGradSpec_: None] := Module[
    {nt, nFine, mCtrl, dim, isTD, uFine, hs, us, fs, j, k, hj, uj, uFinal, idMat, zeroDim,
     phi, cDag, extra, xMat, bMat, gradFine, gradOut, coarseUs, coarseFs,
     hasRC, rcEvalList, rcDtCoarse, rcT0, rcRT, rcWMat, rcTEnds, rcTotal = 0., rcPair, coarseJ,
     hasTimingGrad, gExp, sigmaSelf, sigmaShift, sigmaRun, alphas, isTDCont, dSelf, dShift, alphaJ, gradDtOut},
    nt = Length[uMat];
    nFine = Length[dtVec];
    mCtrl = Length[hCtrls];
    isTD = ArrayDepth[h0] === 3;
    dim = If[isTD, Dimensions[h0][[2]], Length[h0]];
    uFine = If[nStages === 1, uMat, Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ uMat, 1]]];
    hasRC = AssociationQ[runningCostSpec];
    hasTimingGrad = AssociationQ[timingGradSpec];
    If[hasTimingGrad, hs = Table[Null, {nFine}]];
    us = Table[Null, {nFine}];
    fs = Table[Null, {nFine + 1}];
    idMat = Developer`ToPackedArray[IdentityMatrix[dim, WorkingPrecision -> MachinePrecision] + 0. I];
    fs[[1]] = idMat;
    Do[
        hj = 0.5 ((# + ConjugateTranspose[#])) &[If[isTD, h0[[j]], h0] + Sum[uFine[[j, k]] hCtrls[[k]], {k, mCtrl}]];
        If[hasTimingGrad, hs[[j]] = hj];
        uj = MatrixExp[(-I dtVec[[j]]) hj];
        us[[j]] = uj;
        fs[[j + 1]] = uj . fs[[j]],
        {j, nFine}
    ];
    uFinal = fs[[nFine + 1]];
    {phi, cDag, extra} = objCotangentFn[uFinal];
    xMat = ToPackedComplexMatrix[cDag];
    If[hasRC,
        zeroDim = ConstantArray[0. + 0. I, {dim, dim}];
        rcEvalList = runningCostSpec["HilbertTermEvaluators"];
        rcDtCoarse = Which[
            hasTimingGrad, timingGradSpec["DtCoarse"],
            nStages === 1, dtVec,
            True, Total[Partition[dtVec, nStages], {2}]
        ];
        rcT0 = Which[
            hasTimingGrad, Lookup[timingGradSpec, "InitialTimePoint", 0.],
            hasRC, Lookup[runningCostSpec, "InitialTimePoint", 0.],
            True, 0.
        ];
        rcRT = prepareRunningRuntime[runningCostSpec, rcDtCoarse, rcT0];
        rcWMat = rcRT["WeightsMat"];
        rcTEnds = rcRT["TEnds"]
    ];
    If[hasTimingGrad,
        gExp = ConstantArray[0., nt];
        sigmaSelf = ConstantArray[0., nt];
        sigmaShift = ConstantArray[0., nt];
        alphas = timingGradSpec["Alphas"];
        isTDCont = TrueQ[timingGradSpec["IsTDContinuous"]];
        dSelf = Lookup[timingGradSpec, "DSelf", None];
        dShift = Lookup[timingGradSpec, "DShift", None];
        sigmaRun = If[hasRC,
            evalRunningTimeGradient[runningCostSpec, rcEvalList, fs, nStages, rcDtCoarse, rcT0],
            ConstantArray[0., nt]
        ]
    ];
    gradFine = ConstantArray[0., {nFine, mCtrl}];
    Do[
        If[hasRC && Mod[j, nStages] === 0,
            coarseJ = Quotient[j, nStages];
            rcPair = evalRunningBoundary[rcEvalList, rcWMat, fs[[j + 1]], rcTEnds[[coarseJ]], coarseJ, zeroDim];
            If[rcPair === $Failed, Return[$Failed, Module]];
            rcTotal += rcPair[[1]];
            xMat = xMat + rcPair[[2]]
        ];
        bMat = fs[[j + 1]] . xMat; (* U_j F_{j-1} C^\dagger B_j *)
        Do[
            gradFine[[j, k]] = dtVec[[j]] * Im[Tr[bMat . hCtrls[[k]]]],
            {k, mCtrl}
        ];
        If[hasTimingGrad,
            coarseJ = Ceiling[j / nStages];
            alphaJ = alphas[[j]];
            gExp[[coarseJ]] += alphaJ * Im[Total[Flatten[Transpose[bMat] * hs[[j]]]]];
            If[isTDCont,
                sigmaSelf[[coarseJ]] += dtVec[[j]] * Im[Total[Flatten[Transpose[bMat] * dSelf[[j]]]]];
                sigmaShift[[coarseJ]] += dtVec[[j]] * Im[Total[Flatten[Transpose[bMat] * dShift[[j]]]]]
            ]
        ];
        xMat = xMat . us[[j]],
        {j, nFine, 1, -1}
    ];
    gradOut = If[nStages === 1,
        Developer`ToPackedArray[gradFine],
        Developer`ToPackedArray[Total[Partition[gradFine, nStages], {2}]]
    ];
    gradDtOut = If[hasTimingGrad, foldTimingGradient[gExp, sigmaSelf, sigmaShift, sigmaRun], None];
    If[TrueQ[returnProps],
        If[nStages === 1,
            coarseUs = us; coarseFs = fs,
            coarseUs = Table[Fold[Dot, idMat, Reverse[us[[(j - 1) nStages + 1 ;; j nStages]]]], {j, nt}];
            coarseFs = fs[[1 ;; nFine + 1 ;; nStages]]
        ],
        coarseUs = None; coarseFs = None
    ];
    <|
        "Fidelity" -> phi,
        "TerminalFidelity" -> phi,
        "RunningCost" -> rcTotal,
        "Gradient" -> gradOut,
        "DurationGradient" -> gradDtOut,
        "FinalPropagator" -> uFinal,
        "Propagators" -> coarseUs,
        "States" -> coarseFs,
        "Extra" -> extra
    |>
];

(* Multi-ket pure state-vector kernel: propagates D x R state columns instead of D x D propagators *)
EvaluateHilbertStateVector[h0_, hCtrlRows_, hFlatT_, dtVec_List, uMat_, psiMat0_, tarMat_, weights_List, phaseIns_List,
                           nStages_Integer: 1, runningCostSpec_: None, timingGradSpec_: None, returnCache_: False] := Module[
    {nt, nFine, mCtrl, dim, isTD, uFine, rCount, hAll, hs, vs, vsH, evals, pHalfs, vHpsis, psis,
     j, hj, ev, vec, vecH, pHalf, phases, zVec, subVals, phi, scaleVec, chiMat, zeroState,
     vHchi, aEigMat, lam, halfDt, lamRow, phiMat, qMat, wMat, gradFine, gradOut,
     phi1Mats, aEigMats, vHchis, chis, g2Mats, needCache = TrueQ[returnCache], cacheAssoc = None,
     hasRC, rcEvalList, rcDtCoarse, rcT0, rcRT, rcWMat, rcTEnds, rcTotal = 0., rcPair, coarseJ,
     hasTimingGrad, gExp, sigmaSelf, sigmaShift, sigmaRun, alphas, isTDCont, dSelf, dShift, alphaJ, gradDtOut},
    nt = Length[uMat];
    nFine = Length[dtVec];
    mCtrl = Dimensions[uMat][[2]];
    isTD = ArrayDepth[h0] === 3;
    dim = If[isTD, Dimensions[h0][[2]], Length[h0]];
    rCount = Length[weights];
    uFine = If[nStages === 1, uMat, Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ uMat, 1]]];
    hAll = ArrayReshape[uFine . hCtrlRows, {nFine, dim, dim}];
    hasRC = AssociationQ[runningCostSpec];
    hasTimingGrad = AssociationQ[timingGradSpec];
    If[hasTimingGrad, hs = Table[Null, {nFine}]];
    vs = Table[Null, {nFine}];
    vsH = Table[Null, {nFine}];
    evals = Table[Null, {nFine}];
    pHalfs = Table[Null, {nFine}];
    vHpsis = Table[Null, {nFine}];
    psis = Table[Null, {nFine + 1}];
    If[needCache,
        phi1Mats = Table[Null, {nFine}];
        aEigMats = Table[Null, {nFine}];
        vHchis = Table[Null, {nFine}];
        chis = Table[Null, {nFine}];
        g2Mats = Table[Null, {nFine}]
    ];
    psis[[1]] = psiMat0;
    Do[
        hj = If[isTD, h0[[j]], h0] + hAll[[j]];
        hj = 0.5 (hj + ConjugateTranspose[hj]);
        If[hasTimingGrad, hs[[j]] = hj];
        {ev, vec} = Eigensystem[hj];
        ev = Re[ev];
        vec = Transpose[vec];
        vecH = ConjugateTranspose[vec];
        pHalf = Exp[(-0.5 I dtVec[[j]]) ev];
        phases = pHalf * pHalf;
        vHpsis[[j]] = vecH . psis[[j]];
        psis[[j + 1]] = vec . (phases * vHpsis[[j]]);
        vs[[j]] = vec;
        vsH[[j]] = vecH;
        evals[[j]] = ev;
        pHalfs[[j]] = pHalf,
        {j, nFine}
    ];
    zVec = Total[Conjugate[tarMat] * psis[[nFine + 1]], {1}];
    subVals = Table[If[TrueQ[phaseIns[[idx]]], Abs[zVec[[idx]]]^2, Re[zVec[[idx]]]], {idx, rCount}];
    phi = weights . subVals;
    scaleVec = Table[weights[[idx]] * If[TrueQ[phaseIns[[idx]]], 2. zVec[[idx]], 1.], {idx, rCount}];
    chiMat = Developer`ToPackedArray[Transpose[Transpose[tarMat] * scaleVec]];
    If[hasRC,
        zeroState = ConstantArray[0. + 0. I, {dim, rCount}];
        rcEvalList = runningCostSpec["StateVectorTermEvaluators"];
        rcDtCoarse = Which[
            hasTimingGrad, timingGradSpec["DtCoarse"],
            nStages === 1, dtVec,
            True, Total[Partition[dtVec, nStages], {2}]
        ];
        rcT0 = Which[
            hasTimingGrad, Lookup[timingGradSpec, "InitialTimePoint", 0.],
            hasRC, Lookup[runningCostSpec, "InitialTimePoint", 0.],
            True, 0.
        ];
        rcRT = prepareRunningRuntime[runningCostSpec, rcDtCoarse, rcT0];
        rcWMat = rcRT["WeightsMat"];
        rcTEnds = rcRT["TEnds"]
    ];
    If[hasTimingGrad,
        gExp = ConstantArray[0., nt];
        sigmaSelf = ConstantArray[0., nt];
        sigmaShift = ConstantArray[0., nt];
        alphas = timingGradSpec["Alphas"];
        isTDCont = TrueQ[timingGradSpec["IsTDContinuous"]];
        dSelf = Lookup[timingGradSpec, "DSelf", None];
        dShift = Lookup[timingGradSpec, "DShift", None];
        sigmaRun = If[hasRC,
            evalRunningTimeGradient[runningCostSpec, rcEvalList, psis, nStages, rcDtCoarse, rcT0],
            ConstantArray[0., nt]
        ]
    ];
    gradFine = ConstantArray[0., {nFine, mCtrl}];
    Do[
        If[hasRC && Mod[j, nStages] === 0,
            coarseJ = Quotient[j, nStages];
            rcPair = evalRunningBoundary[rcEvalList, rcWMat, psis[[j + 1]], rcTEnds[[coarseJ]], coarseJ, zeroState];
            If[rcPair === $Failed, Return[$Failed, Module]];
            rcTotal += rcPair[[1]];
            chiMat = chiMat + rcPair[[2]]
        ];
        vec = vs[[j]];
        vecH = vsH[[j]];
        lam = evals[[j]];
        pHalf = pHalfs[[j]];
        phases = pHalf * pHalf;
        halfDt = 0.5 dtVec[[j]];
        vHchi = vecH . chiMat;
        aEigMat = vHpsis[[j]] . ConjugateTranspose[vHchi];
        lamRow = ConstantArray[lam, dim];
        phiMat = (Transpose[{pHalf}] . {pHalf}) * Sinc[halfDt (Transpose[lamRow] - lamRow)];
        qMat = aEigMat * phiMat;
        wMat = vec . qMat . vecH;
        gradFine[[j]] = dtVec[[j]] * Im[hFlatT . Flatten[wMat]];
        If[needCache,
            chis[[j]] = chiMat;
            vHchis[[j]] = vHchi;
            aEigMats[[j]] = aEigMat;
            phi1Mats[[j]] = phiMat;
            g2Mats[[j]] = buildG2Tensor[lam, dtVec[[j]], phiMat]
        ];
        If[hasTimingGrad,
            coarseJ = Ceiling[j / nStages];
            alphaJ = alphas[[j]];
            gExp[[coarseJ]] += alphaJ * Im[Total[Flatten[Conjugate[chiMat] * (hs[[j]] . psis[[j + 1]])]]];
            If[isTDCont,
                sigmaSelf[[coarseJ]] += dtVec[[j]] * Im[Total[Flatten[Transpose[dSelf[[j]]] * wMat]]];
                sigmaShift[[coarseJ]] += dtVec[[j]] * Im[Total[Flatten[Transpose[dShift[[j]]] * wMat]]]
            ]
        ];
        chiMat = vec . (Conjugate[phases] * vHchi),
        {j, nFine, 1, -1}
    ];
    gradOut = If[nStages === 1,
        Developer`ToPackedArray[gradFine],
        Developer`ToPackedArray[Total[Partition[gradFine, nStages], {2}]]
    ];
    gradDtOut = If[hasTimingGrad, foldTimingGradient[gExp, sigmaSelf, sigmaShift, sigmaRun], None];
    If[needCache,
        cacheAssoc = <|
            "Kernel" -> "StateVector",
            "Dim" -> dim,
            "NumSlices" -> nt,
            "NumFine" -> nFine,
            "NumStages" -> nStages,
            "NumControls" -> mCtrl,
            "NumKets" -> rCount,
            "DtFine" -> Developer`ToPackedArray[N[dtVec]],
            "Vs" -> vs,
            "VsH" -> vsH,
            "Evals" -> evals,
            "PHalfs" -> pHalfs,
            "Phi1Mats" -> phi1Mats,
            "AEigMats" -> aEigMats,
            "G2Mats" -> g2Mats,
            "VHpsis" -> vHpsis,
            "Psis" -> psis,
            "VHchis" -> vHchis,
            "Chis" -> chis,
            "TarMat" -> tarMat,
            "Weights" -> Developer`ToPackedArray[N[weights]],
            "PhaseInsensitive" -> phaseIns,
            "ZVec" -> zVec,
            "ControlRows" -> hCtrlRows,
            "ControlFlatT" -> hFlatT
        |>
    ];
    <|
        "Fidelity" -> phi,
        "TerminalFidelity" -> phi,
        "RunningCost" -> rcTotal,
        "Gradient" -> gradOut,
        "DurationGradient" -> gradDtOut,
        "FinalState" -> psis[[nFine + 1]],
        "SubFidelities" -> subVals,
        If[needCache, "EvalState" -> cacheAssoc, Nothing]
    |>
];

(* Open-system Liouville-space GRAPE kernel via MatrixExp and augmented Fréchet derivative *)
EvaluateLiouvilleExact[l0_, lCtrls_List, dtVec_List, uMat_, rhoMat0_, obsMat_, weights_List, objModes_List, scales_List,
                       nStages_Integer: 1, runningCostSpec_: None, timingGradSpec_: None, returnCache_: False] := Module[
    {nt, nFine, mCtrl, dl, isTD, uFine, rCount, zMat, zeroVec, zeroState, ls, ps, rs, j, k, r, lj, pj, rFinal,
     subVals, zVec, phi, csList, csMat, blk, vecIn, dr, expB, dP, gradFine, gradOut, idMat, coarsePs,
     csMats, dPkMats, needCache = TrueQ[returnCache], cacheAssoc = None,
     hasRC, rcEvalList, rcDtCoarse, rcT0, rcRT, rcWMat, rcTEnds, rcTotal = 0., rcPair, coarseJ,
     hasTimingGrad, gExp, sigmaSelf, sigmaShift, sigmaRun, alphas, isTDCont, dSelf, dShift, alphaJ, gradDtOut},
    nt = Length[uMat];
    nFine = Length[dtVec];
    mCtrl = Length[lCtrls];
    isTD = ArrayDepth[l0] === 3;
    dl = If[isTD, Dimensions[l0][[2]], Length[l0]];
    rCount = Length[weights];
    uFine = If[nStages === 1, uMat, Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ uMat, 1]]];
    zMat = ConstantArray[0. + 0. I, {dl, dl}];
    zeroVec = ConstantArray[0. + 0. I, dl];
    hasRC = AssociationQ[runningCostSpec];
    hasTimingGrad = AssociationQ[timingGradSpec];
    ls = Table[Null, {nFine}];
    ps = Table[Null, {nFine}];
    rs = Table[Null, {nFine + 1}];
    If[needCache,
        csMats = Table[Null, {nFine}];
        dPkMats = Table[Null, {nFine}, {mCtrl}]
    ];
    rs[[1]] = rhoMat0;
    Do[
        lj = If[isTD, l0[[j]], l0] + Sum[uFine[[j, k]] lCtrls[[k]], {k, mCtrl}];
        pj = MatrixExp[dtVec[[j]] lj];
        ls[[j]] = lj;
        ps[[j]] = pj;
        rs[[j + 1]] = pj . rs[[j]],
        {j, nFine}
    ];
    rFinal = rs[[nFine + 1]];
    zVec = Total[Conjugate[obsMat] * rFinal, {1}];
    subVals = Table[
        Switch[objModes[[r]],
            "Magnitude", Abs[zVec[[r]] / scales[[r]]]^2,
            "Distance", 1. - Re[Total[Abs[rFinal[[All, r]] - obsMat[[All, r]]]^2]] / scales[[r]],
            _, Re[zVec[[r]]] / scales[[r]]
        ],
        {r, rCount}
    ];
    phi = weights . subVals;
    csList = Table[
        With[{oCol = obsMat[[All, r]], sc = scales[[r]], w = weights[[r]]},
            Switch[objModes[[r]],
                "Magnitude", (2. w / (sc^2)) * zVec[[r]] * oCol,
                "Distance", -(2. w / sc) * (rFinal[[All, r]] - oCol),
                _, (w / sc) * oCol
            ]
        ],
        {r, rCount}
    ];
    csMat = Developer`ToPackedArray[Transpose[csList]];
    If[hasRC,
        zeroState = ConstantArray[0. + 0. I, {dl, rCount}];
        rcEvalList = runningCostSpec["LiouvilleTermEvaluators"];
        rcDtCoarse = Which[
            hasTimingGrad, timingGradSpec["DtCoarse"],
            nStages === 1, dtVec,
            True, Total[Partition[dtVec, nStages], {2}]
        ];
        rcT0 = Which[
            hasTimingGrad, Lookup[timingGradSpec, "InitialTimePoint", 0.],
            hasRC, Lookup[runningCostSpec, "InitialTimePoint", 0.],
            True, 0.
        ];
        rcRT = prepareRunningRuntime[runningCostSpec, rcDtCoarse, rcT0];
        rcWMat = rcRT["WeightsMat"];
        rcTEnds = rcRT["TEnds"]
    ];
    If[hasTimingGrad,
        gExp = ConstantArray[0., nt];
        sigmaSelf = ConstantArray[0., nt];
        sigmaShift = ConstantArray[0., nt];
        alphas = timingGradSpec["Alphas"];
        isTDCont = TrueQ[timingGradSpec["IsTDContinuous"]];
        dSelf = Lookup[timingGradSpec, "DSelf", None];
        dShift = Lookup[timingGradSpec, "DShift", None];
        sigmaRun = If[hasRC,
            evalRunningTimeGradient[runningCostSpec, rcEvalList, rs, nStages, rcDtCoarse, rcT0],
            ConstantArray[0., nt]
        ]
    ];
    gradFine = ConstantArray[0., {nFine, mCtrl}];
    Do[
        If[hasRC && Mod[j, nStages] === 0,
            coarseJ = Quotient[j, nStages];
            rcPair = evalRunningBoundary[rcEvalList, rcWMat, rs[[j + 1]], rcTEnds[[coarseJ]], coarseJ, zeroState];
            If[rcPair === $Failed, Return[$Failed, Module]];
            rcTotal += rcPair[[1]];
            csMat = csMat + rcPair[[2]]
        ];
        If[needCache, csMats[[j]] = csMat];
        lj = ls[[j]];
        Do[
            blk = Developer`ToPackedArray[ArrayFlatten[{{dtVec[[j]] lj, dtVec[[j]] lCtrls[[k]]}, {zMat, dtVec[[j]] lj}}]];
            If[rCount === 1 && !needCache,
                vecIn = Developer`ToPackedArray[Join[zeroVec, rs[[j, All, 1]]]];
                dr = (MatrixExp[blk, vecIn])[[1 ;; dl]];
                gradFine[[j, k]] = Re[Conjugate[csMat[[All, 1]]] . dr],
                expB = MatrixExp[blk];
                dP = expB[[1 ;; dl, dl + 1 ;; 2 dl]];
                If[needCache, dPkMats[[j, k]] = dP];
                gradFine[[j, k]] = Re[Total[Flatten[Conjugate[csMat] * (dP . rs[[j]])]]]
            ],
            {k, mCtrl}
        ];
        If[hasTimingGrad,
            coarseJ = Ceiling[j / nStages];
            alphaJ = alphas[[j]];
            gExp[[coarseJ]] += alphaJ * Re[Total[Flatten[Conjugate[csMat] * (lj . rs[[j + 1]])]]];
            If[isTDCont,
                blk = Developer`ToPackedArray[ArrayFlatten[{{dtVec[[j]] lj, dtVec[[j]] dSelf[[j]]}, {zMat, dtVec[[j]] lj}}]];
                expB = MatrixExp[blk];
                dP = expB[[1 ;; dl, dl + 1 ;; 2 dl]];
                sigmaSelf[[coarseJ]] += Re[Total[Flatten[Conjugate[csMat] * (dP . rs[[j]])]]];
                blk = Developer`ToPackedArray[ArrayFlatten[{{dtVec[[j]] lj, dtVec[[j]] dShift[[j]]}, {zMat, dtVec[[j]] lj}}]];
                expB = MatrixExp[blk];
                dP = expB[[1 ;; dl, dl + 1 ;; 2 dl]];
                sigmaShift[[coarseJ]] += Re[Total[Flatten[Conjugate[csMat] * (dP . rs[[j]])]]]
            ]
        ];
        csMat = ConjugateTranspose[ps[[j]]] . csMat,
        {j, nFine, 1, -1}
    ];
    gradOut = If[nStages === 1,
        Developer`ToPackedArray[gradFine],
        Developer`ToPackedArray[Total[Partition[gradFine, nStages], {2}]]
    ];
    gradDtOut = If[hasTimingGrad, foldTimingGradient[gExp, sigmaSelf, sigmaShift, sigmaRun], None];
    If[nStages === 1,
        coarsePs = ps,
        idMat = Developer`ToPackedArray[IdentityMatrix[dl, WorkingPrecision -> MachinePrecision] + 0. I];
        coarsePs = Table[Fold[Dot, idMat, Reverse[ps[[(j - 1) nStages + 1 ;; j nStages]]]], {j, nt}]
    ];
    If[needCache,
        cacheAssoc = <|
            "Kernel" -> "Liouville",
            "Dim" -> dl,
            "NumSlices" -> nt,
            "NumFine" -> nFine,
            "NumStages" -> nStages,
            "NumControls" -> mCtrl,
            "NumKets" -> rCount,
            "DtFine" -> Developer`ToPackedArray[N[dtVec]],
            "Ls" -> ls,
            "Ps" -> ps,
            "Rs" -> rs,
            "CsMats" -> csMats,
            "DPkMats" -> dPkMats,
            "ControlOps" -> lCtrls,
            "ObsMat" -> obsMat,
            "Weights" -> Developer`ToPackedArray[N[weights]],
            "ObjModes" -> objModes,
            "Scales" -> Developer`ToPackedArray[N[scales]],
            "ZVec" -> zVec
        |>
    ];
    <|
        "Fidelity" -> phi,
        "TerminalFidelity" -> phi,
        "RunningCost" -> rcTotal,
        "Gradient" -> gradOut,
        "DurationGradient" -> gradDtOut,
        "FinalState" -> rFinal,
        "Propagators" -> coarsePs,
        "SubFidelities" -> subVals,
        If[needCache, "EvalState" -> cacheAssoc, Nothing]
    |>
];

(* Forward-only propagation helper *)
PropagateControls[compiled_Association, uMat_] := Module[
    {space = compiled["Space"], dtFine = Lookup[compiled, "FineTimeSteps", compiled["TimeSteps"]],
     nStages = Lookup[compiled, "NumStages", 1], nt = Length[uMat], nFine, uFine,
     d, h0, hCtrls, l0, lCtrls, isTD, us, states, j, k, matj, pj, idMat, coarseUs, coarseStates},
    nFine = Length[dtFine];
    uFine = If[nStages === 1, uMat, Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ uMat, 1]]];
    If[space === "Hilbert",
        h0 = compiled["Drift"];
        hCtrls = compiled["Controls"];
        isTD = ArrayDepth[h0] === 3;
        d = If[isTD, Dimensions[h0][[2]], Length[h0]];
        us = Table[Null, {nFine}];
        states = Table[Null, {nFine + 1}];
        idMat = Developer`ToPackedArray[IdentityMatrix[d, WorkingPrecision -> MachinePrecision] + 0. I];
        states[[1]] = idMat;
        Do[
            matj = 0.5 (# + ConjugateTranspose[#]) &[If[isTD, h0[[j]], h0] + Sum[uFine[[j, k]] hCtrls[[k]], {k, Length[hCtrls]}]];
            pj = MatrixExp[(-I dtFine[[j]]) matj];
            us[[j]] = pj;
            states[[j + 1]] = pj . states[[j]],
            {j, nFine}
        ];
        If[nStages === 1,
            coarseUs = us; coarseStates = states,
            coarseUs = Table[Fold[Dot, idMat, Reverse[us[[(j - 1) nStages + 1 ;; j nStages]]]], {j, nt}];
            coarseStates = states[[1 ;; nFine + 1 ;; nStages]]
        ];
        <|"Space" -> "Hilbert", "FinalPropagator" -> states[[nFine + 1]],
          "Propagators" -> coarseUs, "CumulativePropagators" -> coarseStates|>,
        l0 = compiled["Drift"];
        lCtrls = compiled["Controls"];
        isTD = ArrayDepth[l0] === 3;
        d = If[isTD, Dimensions[l0][[2]], Length[l0]];
        us = Table[Null, {nFine}];
        states = Table[Null, {nFine + 1}];
        states[[1]] = compiled["Objective", "InitialMatrix"];
        Do[
            matj = If[isTD, l0[[j]], l0] + Sum[uFine[[j, k]] lCtrls[[k]], {k, Length[lCtrls]}];
            pj = MatrixExp[dtFine[[j]] matj];
            us[[j]] = pj;
            states[[j + 1]] = pj . states[[j]],
            {j, nFine}
        ];
        If[nStages === 1,
            coarseUs = us; coarseStates = states,
            idMat = Developer`ToPackedArray[IdentityMatrix[d, WorkingPrecision -> MachinePrecision] + 0. I];
            coarseUs = Table[Fold[Dot, idMat, Reverse[us[[(j - 1) nStages + 1 ;; j nStages]]]], {j, nt}];
            coarseStates = states[[1 ;; nFine + 1 ;; nStages]]
        ];
        <|"Space" -> "Liouville", "FinalState" -> states[[nFine + 1]],
          "Propagators" -> coarseUs, "States" -> coarseStates|>
    ]
];

(* Forward-only objective kernels (zero backward sweeps or block-MatrixExp derivatives) *)
EvaluateHilbertSpectralObjectiveOnly[h0_, hCtrlRows_, dtVec_List, uMat_, objCotangentFn_,
                                     nStages_Integer: 1, runningCostSpec_: None] := Module[
    {nFine, dim, isTD, uFine, hAll, hasRC, rcEvalList, rcDtCoarse, rcT0, rcRT, rcWMat, rcTEnds,
     rcTotal = 0., fMat, zeroDim, j, hj, ev, vec, phases, uj, phi, cDag, extra, coarseJ, rcPair},
    nFine = Length[dtVec];
    isTD = ArrayDepth[h0] === 3;
    dim = If[isTD, Dimensions[h0][[2]], Length[h0]];
    uFine = If[nStages === 1, uMat, Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ uMat, 1]]];
    hAll = ArrayReshape[uFine . hCtrlRows, {nFine, dim, dim}];
    hasRC = AssociationQ[runningCostSpec];
    If[hasRC,
        zeroDim = ConstantArray[0. + 0. I, {dim, dim}];
        rcEvalList = runningCostSpec["HilbertTermEvaluators"];
        rcDtCoarse = If[nStages === 1, dtVec, Total[Partition[dtVec, nStages], {2}]];
        rcT0 = Lookup[runningCostSpec, "InitialTimePoint", 0.];
        rcRT = prepareRunningRuntime[runningCostSpec, rcDtCoarse, rcT0];
        rcWMat = rcRT["WeightsMat"];
        rcTEnds = rcRT["TEnds"]
    ];
    fMat = Developer`ToPackedArray[IdentityMatrix[dim, WorkingPrecision -> MachinePrecision] + 0. I];
    Do[
        hj = If[isTD, h0[[j]], h0] + hAll[[j]];
        hj = 0.5 (hj + ConjugateTranspose[hj]);
        {ev, vec} = Eigensystem[hj];
        ev = Re[ev];
        vec = Transpose[vec];
        phases = Exp[(-I dtVec[[j]]) ev];
        uj = vec . (phases * ConjugateTranspose[vec]);
        fMat = uj . fMat;
        If[hasRC && Mod[j, nStages] === 0,
            coarseJ = Quotient[j, nStages];
            rcPair = evalRunningBoundary[rcEvalList, rcWMat, fMat, rcTEnds[[coarseJ]], coarseJ, zeroDim];
            If[rcPair === $Failed, Return[$Failed, Module]];
            rcTotal += rcPair[[1]]
        ],
        {j, nFine}
    ];
    {phi, cDag, extra} = objCotangentFn[fMat];
    <|
        "Fidelity" -> phi,
        "TerminalFidelity" -> phi,
        "RunningCost" -> rcTotal,
        "Extra" -> extra
    |>
];

EvaluateHilbertStateVectorObjectiveOnly[h0_, hCtrlRows_, dtVec_List, uMat_, psiMat0_, tarMat_, weights_List, phaseIns_List,
                                        nStages_Integer: 1, runningCostSpec_: None] := Module[
    {nFine, dim, isTD, rCount, uFine, hAll, hasRC, rcEvalList, rcDtCoarse, rcT0, rcRT, rcWMat, rcTEnds,
     rcTotal = 0., psiMat, zeroState, j, hj, ev, vec, vecH, phases, zVec, subVals, phi, coarseJ, rcPair},
    nFine = Length[dtVec];
    isTD = ArrayDepth[h0] === 3;
    dim = If[isTD, Dimensions[h0][[2]], Length[h0]];
    rCount = Length[weights];
    uFine = If[nStages === 1, uMat, Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ uMat, 1]]];
    hAll = ArrayReshape[uFine . hCtrlRows, {nFine, dim, dim}];
    hasRC = AssociationQ[runningCostSpec];
    If[hasRC,
        zeroState = ConstantArray[0. + 0. I, {dim, rCount}];
        rcEvalList = runningCostSpec["StateVectorTermEvaluators"];
        rcDtCoarse = If[nStages === 1, dtVec, Total[Partition[dtVec, nStages], {2}]];
        rcT0 = Lookup[runningCostSpec, "InitialTimePoint", 0.];
        rcRT = prepareRunningRuntime[runningCostSpec, rcDtCoarse, rcT0];
        rcWMat = rcRT["WeightsMat"];
        rcTEnds = rcRT["TEnds"]
    ];
    psiMat = psiMat0;
    Do[
        hj = If[isTD, h0[[j]], h0] + hAll[[j]];
        hj = 0.5 (hj + ConjugateTranspose[hj]);
        {ev, vec} = Eigensystem[hj];
        ev = Re[ev];
        vec = Transpose[vec];
        vecH = ConjugateTranspose[vec];
        phases = Exp[(-I dtVec[[j]]) ev];
        psiMat = vec . (phases * (vecH . psiMat));
        If[hasRC && Mod[j, nStages] === 0,
            coarseJ = Quotient[j, nStages];
            rcPair = evalRunningBoundary[rcEvalList, rcWMat, psiMat, rcTEnds[[coarseJ]], coarseJ, zeroState];
            If[rcPair === $Failed, Return[$Failed, Module]];
            rcTotal += rcPair[[1]]
        ],
        {j, nFine}
    ];
    zVec = Total[Conjugate[tarMat] * psiMat, {1}];
    subVals = Table[If[TrueQ[phaseIns[[idx]]], Abs[zVec[[idx]]]^2, Re[zVec[[idx]]]], {idx, rCount}];
    phi = weights . subVals;
    <|
        "Fidelity" -> phi,
        "TerminalFidelity" -> phi,
        "RunningCost" -> rcTotal,
        "SubFidelities" -> subVals
    |>
];

EvaluateLiouvilleObjectiveOnly[l0_, lCtrls_List, dtVec_List, uMat_, rhoMat0_, obsMat_, weights_List, objModes_List, scales_List,
                               nStages_Integer: 1, runningCostSpec_: None] := Module[
    {nFine, mCtrl, dl, isTD, rCount, uFine, hasRC, rcEvalList, rcDtCoarse, rcT0, rcRT, rcWMat, rcTEnds,
     rcTotal = 0., rhoMat, zeroState, j, k, lj, zVec, subVals, phi, coarseJ, rcPair},
    nFine = Length[dtVec];
    mCtrl = Length[lCtrls];
    isTD = ArrayDepth[l0] === 3;
    dl = If[isTD, Dimensions[l0][[2]], Length[l0]];
    rCount = Length[weights];
    uFine = If[nStages === 1, uMat, Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ uMat, 1]]];
    hasRC = AssociationQ[runningCostSpec];
    If[hasRC,
        zeroState = ConstantArray[0. + 0. I, {dl, rCount}];
        rcEvalList = runningCostSpec["LiouvilleTermEvaluators"];
        rcDtCoarse = If[nStages === 1, dtVec, Total[Partition[dtVec, nStages], {2}]];
        rcT0 = Lookup[runningCostSpec, "InitialTimePoint", 0.];
        rcRT = prepareRunningRuntime[runningCostSpec, rcDtCoarse, rcT0];
        rcWMat = rcRT["WeightsMat"];
        rcTEnds = rcRT["TEnds"]
    ];
    rhoMat = rhoMat0;
    Do[
        lj = If[isTD, l0[[j]], l0] + Sum[uFine[[j, k]] lCtrls[[k]], {k, mCtrl}];
        rhoMat = If[rCount === 1,
            Developer`ToPackedArray[{MatrixExp[dtVec[[j]] lj, rhoMat[[All, 1]]]}] // Transpose,
            MatrixExp[dtVec[[j]] lj] . rhoMat
        ];
        If[hasRC && Mod[j, nStages] === 0,
            coarseJ = Quotient[j, nStages];
            rcPair = evalRunningBoundary[rcEvalList, rcWMat, rhoMat, rcTEnds[[coarseJ]], coarseJ, zeroState];
            If[rcPair === $Failed, Return[$Failed, Module]];
            rcTotal += rcPair[[1]]
        ],
        {j, nFine}
    ];
    zVec = Total[Conjugate[obsMat] * rhoMat, {1}];
    subVals = Table[
        Switch[objModes[[r]],
            "Magnitude", Abs[zVec[[r]] / scales[[r]]]^2,
            "Distance", 1. - Re[Total[Abs[rhoMat[[All, r]] - obsMat[[All, r]]]^2]] / scales[[r]],
            _, Re[zVec[[r]]] / scales[[r]]
        ],
        {r, rCount}
    ];
    phi = weights . subVals;
    <|
        "Fidelity" -> phi,
        "TerminalFidelity" -> phi,
        "RunningCost" -> rcTotal,
        "SubFidelities" -> subVals
    |>
];

(* Matrix-free O(N_t) analytical HVP & Goodwin-Kuprov interaction-frame Hessian engines *)
EvaluateHilbertSpectralHVP[evalState_Association, vMat_, objDirCotFn_] := Module[
    {state, dim, nt, nFine, nStages, mCtrl, dtFine, vs, vsH, phi1Mats, aTilMats, g2Mats,
     us, fs, xMats, uFinal, hCtrlRows, hFlatT, vFine, vAll, vTilMats, dUs, dFs,
     j, vec, vecH, dU, dXMat, dAMat, wDot, gCross, w2Mat, gDiag, hvpFine},
    state = EnsureSpectralG2Cache[evalState];
    dim = state["Dim"];
    nt = state["NumSlices"];
    nFine = state["NumFine"];
    nStages = state["NumStages"];
    mCtrl = state["NumControls"];
    dtFine = state["DtFine"];
    vs = state["Vs"];
    vsH = state["VsH"];
    phi1Mats = state["Phi1Mats"];
    aTilMats = state["ATilMats"];
    g2Mats = state["G2Mats"];
    us = state["Us"];
    fs = state["Fs"];
    xMats = state["XMats"];
    uFinal = state["FinalPropagator"];
    hCtrlRows = state["ControlRows"];
    hFlatT = state["ControlFlatT"];

    vFine = If[nStages === 1, ToPackedRealMatrix[vMat],
        Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ ToPackedRealMatrix[vMat], 1]]
    ];
    vAll = ArrayReshape[vFine . hCtrlRows, {nFine, dim, dim}];
    vTilMats = Table[Null, {nFine}];
    dUs = Table[Null, {nFine}];
    dFs = Table[Null, {nFine + 1}];
    dFs[[1]] = ToPackedComplexMatrix[ConstantArray[0., {dim, dim}]];
    Do[
        vec = vs[[j]]; vecH = vsH[[j]];
        vTilMats[[j]] = vecH . vAll[[j]] . vec;
        dU = (-I * dtFine[[j]]) * (vec . (vTilMats[[j]] * phi1Mats[[j]]) . vecH);
        dUs[[j]] = dU;
        dFs[[j + 1]] = us[[j]] . dFs[[j]] + dU . fs[[j]],
        {j, nFine}
    ];

    dXMat = ToPackedComplexMatrix[objDirCotFn[uFinal, dFs[[nFine + 1]]]];
    hvpFine = ConstantArray[0., {nFine, mCtrl}];
    Do[
        vec = vs[[j]]; vecH = vsH[[j]];
        dAMat = dFs[[j]] . xMats[[j]] + fs[[j]] . dXMat;
        wDot = vec . ((vecH . dAMat . vec) * phi1Mats[[j]]) . vecH;
        gCross = dtFine[[j]] * Im[hFlatT . Flatten[wDot]];
        w2Mat = vec . contractG2[aTilMats[[j]], vTilMats[[j]], g2Mats[[j]]] . vecH;
        gDiag = -(dtFine[[j]]^2) * Re[hFlatT . Flatten[w2Mat]];
        hvpFine[[j]] = gCross + gDiag;
        dXMat = dXMat . us[[j]] + xMats[[j]] . dUs[[j]],
        {j, nFine, 1, -1}
    ];
    If[nStages === 1,
        Developer`ToPackedArray[hvpFine],
        Developer`ToPackedArray[Total[Partition[hvpFine, nStages], {2}]]
    ]
];

EvaluateHilbertStateVectorHVP[evalState_Association, vMat_] := Module[
    {state, dim, nt, nFine, nStages, mCtrl, rCount, dtFine, vs, vsH, pHalfs, phi1Mats,
     aEigMats, g2Mats, vHpsis, psis, vHchis, chis, tarMat, weights, phaseIns, hCtrlRows, hFlatT,
     vFine, vAll, vTilMats, dUTilMats, dPsis, dPsiFinal, dZVec, dScaleVec, dChiMat,
     j, r, vec, vecH, vTil, dUTil, phases, dVHpsi, dVHchi, dAEig, wDot, gCross, w2Mat, gDiag, hvpFine},
    state = EnsureSpectralG2Cache[evalState];
    dim = state["Dim"];
    nt = state["NumSlices"];
    nFine = state["NumFine"];
    nStages = state["NumStages"];
    mCtrl = state["NumControls"];
    rCount = state["NumKets"];
    dtFine = state["DtFine"];
    vs = state["Vs"];
    vsH = state["VsH"];
    pHalfs = state["PHalfs"];
    phi1Mats = state["Phi1Mats"];
    aEigMats = state["AEigMats"];
    g2Mats = state["G2Mats"];
    vHpsis = state["VHpsis"];
    psis = state["Psis"];
    vHchis = state["VHchis"];
    tarMat = state["TarMat"];
    weights = state["Weights"];
    phaseIns = state["PhaseInsensitive"];
    hCtrlRows = state["ControlRows"];
    hFlatT = state["ControlFlatT"];

    vFine = If[nStages === 1, ToPackedRealMatrix[vMat],
        Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ ToPackedRealMatrix[vMat], 1]]
    ];
    vAll = ArrayReshape[vFine . hCtrlRows, {nFine, dim, dim}];
    vTilMats = Table[Null, {nFine}];
    dUTilMats = Table[Null, {nFine}];
    dPsis = Table[Null, {nFine + 1}];
    dPsis[[1]] = ToPackedComplexMatrix[ConstantArray[0., {dim, rCount}]];
    Do[
        vec = vs[[j]]; vecH = vsH[[j]];
        vTil = vecH . vAll[[j]] . vec;
        dUTil = (-I * dtFine[[j]]) * (vTil * phi1Mats[[j]]);
        phases = pHalfs[[j]] * pHalfs[[j]];
        vTilMats[[j]] = vTil;
        dUTilMats[[j]] = dUTil;
        dPsis[[j + 1]] = vec . (phases * (vecH . dPsis[[j]]) + dUTil . vHpsis[[j]]),
        {j, nFine}
    ];
    dPsiFinal = dPsis[[nFine + 1]];
    dZVec = Total[Conjugate[tarMat] * dPsiFinal, {1}];
    dScaleVec = Table[If[TrueQ[phaseIns[[r]]], 2. * weights[[r]] * dZVec[[r]], 0. + 0. I], {r, rCount}];
    dChiMat = Developer`ToPackedArray[Transpose[Transpose[tarMat] * dScaleVec]];
    hvpFine = ConstantArray[0., {nFine, mCtrl}];
    Do[
        vec = vs[[j]]; vecH = vsH[[j]];
        phases = pHalfs[[j]] * pHalfs[[j]];
        dVHpsi = vecH . dPsis[[j]];
        dVHchi = vecH . dChiMat;
        dAEig = dVHpsi . ConjugateTranspose[vHchis[[j]]] + vHpsis[[j]] . ConjugateTranspose[dVHchi];
        wDot = vec . (dAEig * phi1Mats[[j]]) . vecH;
        gCross = dtFine[[j]] * Im[hFlatT . Flatten[wDot]];
        w2Mat = vec . contractG2[aEigMats[[j]], vTilMats[[j]], g2Mats[[j]]] . vecH;
        gDiag = -(dtFine[[j]]^2) * Re[hFlatT . Flatten[w2Mat]];
        hvpFine[[j]] = gCross + gDiag;
        dChiMat = vec . (Conjugate[phases] * dVHchi + ConjugateTranspose[dUTilMats[[j]]] . vHchis[[j]]),
        {j, nFine, 1, -1}
    ];
    If[nStages === 1,
        Developer`ToPackedArray[hvpFine],
        Developer`ToPackedArray[Total[Partition[hvpFine, nStages], {2}]]
    ]
];

EvaluateLiouvilleHVP[evalState_Association, vMat_] := Module[
    {dl, nt, nFine, nStages, mCtrl, rCount, dtFine, ls, ps, rs, csMats, dPkMats, lCtrls,
     obsMat, weights, objModes, scales, vFine, zMat, zeroVec, dPs, dRs, j, k, r,
     lj, vDir, blk, expB, dPj, dRFinal, dZVec, dCsList, dCsMat, hvpFine,
     d2Rho, blk3A, blk3B, vecIn3, exp3A, exp3B, d2P, dP12A, dP12B},
    dl = evalState["Dim"];
    nt = evalState["NumSlices"];
    nFine = evalState["NumFine"];
    nStages = evalState["NumStages"];
    mCtrl = evalState["NumControls"];
    rCount = evalState["NumKets"];
    dtFine = evalState["DtFine"];
    ls = evalState["Ls"];
    ps = evalState["Ps"];
    rs = evalState["Rs"];
    csMats = evalState["CsMats"];
    dPkMats = evalState["DPkMats"];
    lCtrls = evalState["ControlOps"];
    obsMat = evalState["ObsMat"];
    weights = evalState["Weights"];
    objModes = evalState["ObjModes"];
    scales = evalState["Scales"];

    vFine = If[nStages === 1, ToPackedRealMatrix[vMat],
        Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ ToPackedRealMatrix[vMat], 1]]
    ];
    zMat = ConstantArray[0. + 0. I, {dl, dl}];
    zeroVec = ConstantArray[0. + 0. I, dl];
    dPs = Table[Null, {nFine}];
    dRs = Table[Null, {nFine + 1}];
    dRs[[1]] = ToPackedComplexMatrix[ConstantArray[0., {dl, rCount}]];
    Do[
        lj = dtFine[[j]] * ls[[j]];
        vDir = dtFine[[j]] * Sum[vFine[[j, k]] * lCtrls[[k]], {k, mCtrl}];
        blk = Developer`ToPackedArray[ArrayFlatten[{{lj, vDir}, {zMat, lj}}]];
        expB = MatrixExp[blk];
        dPj = expB[[1 ;; dl, dl + 1 ;; 2 dl]];
        dPs[[j]] = dPj;
        dRs[[j + 1]] = ps[[j]] . dRs[[j]] + dPj . rs[[j]],
        {j, nFine}
    ];
    dRFinal = dRs[[nFine + 1]];
    dZVec = Total[Conjugate[obsMat] * dRFinal, {1}];
    dCsList = Table[
        With[{oCol = obsMat[[All, r]], sc = scales[[r]], w = weights[[r]]},
            Switch[objModes[[r]],
                "Magnitude", (2. * w / (sc^2)) * dZVec[[r]] * oCol,
                "Distance", -(2. * w / sc) * dRFinal[[All, r]],
                _, zeroVec
            ]
        ],
        {r, rCount}
    ];
    dCsMat = Developer`ToPackedArray[Transpose[dCsList]];
    hvpFine = ConstantArray[0., {nFine, mCtrl}];
    Do[
        lj = dtFine[[j]] * ls[[j]];
        vDir = dtFine[[j]] * Sum[vFine[[j, k]] * lCtrls[[k]], {k, mCtrl}];
        Do[
            With[{dPk = dPkMats[[j, k]], lk = dtFine[[j]] * lCtrls[[k]]},
                blk3A = Developer`ToPackedArray[ArrayFlatten[{{lj, vDir, zMat}, {zMat, lj, lk}, {zMat, zMat, lj}}]];
                blk3B = Developer`ToPackedArray[ArrayFlatten[{{lj, lk, zMat}, {zMat, lj, vDir}, {zMat, zMat, lj}}]];
                If[rCount === 1,
                    vecIn3 = Developer`ToPackedArray[Join[zeroVec, zeroVec, rs[[j, All, 1]]]];
                    d2Rho = (MatrixExp[blk3A, vecIn3])[[1 ;; dl]] + (MatrixExp[blk3B, vecIn3])[[1 ;; dl]];
                    hvpFine[[j, k]] = Re[
                        Conjugate[dCsMat[[All, 1]]] . (dPk . rs[[j, All, 1]]) +
                        Conjugate[csMats[[j, All, 1]]] . (dPk . dRs[[j, All, 1]] + d2Rho)
                    ],
                    exp3A = MatrixExp[blk3A];
                    exp3B = MatrixExp[blk3B];
                    d2P = exp3A[[1 ;; dl, 2 dl + 1 ;; 3 dl]] + exp3B[[1 ;; dl, 2 dl + 1 ;; 3 dl]];
                    hvpFine[[j, k]] = Re[Total[Flatten[
                        Conjugate[dCsMat] * (dPk . rs[[j]]) +
                        Conjugate[csMats[[j]]] * (dPk . dRs[[j]] + d2P . rs[[j]])
                    ]]]
                ]
            ],
            {k, mCtrl}
        ];
        dCsMat = ConjugateTranspose[ps[[j]]] . dCsMat + ConjugateTranspose[dPs[[j]]] . csMats[[j]],
        {j, nFine, 1, -1}
    ];
    If[nStages === 1,
        Developer`ToPackedArray[hvpFine],
        Developer`ToPackedArray[Total[Partition[hvpFine, nStages], {2}]]
    ]
];

(* Fast Goodwin-Kuprov interaction-frame exact explicit Hessian for Spectral kernel *)
EvaluateHilbertSpectralHessian[evalState_Association, hCtrls_List, objDirCotFn_] := Module[
    {state, dim, nt, nFine, nStages, mCtrl, dtFine, vs, vsH, phi1Mats, aTilMats, g2Mats,
     fs, uFinal, cDag, hFlatT, zN, nTotFine, pList, dFList, dCDagList, hTilGrid,
     idx, q, k, vec, vecH, hTil, dUqk, pMat, dFMat, matL, matR, dotCDagMat, dFTransMat,
     hProp, diagBlk, w2Mat, i0, i1, hTerm, hFine, hCoarse},
    state = EnsureSpectralG2Cache[evalState];
    dim = state["Dim"];
    nt = state["NumSlices"];
    nFine = state["NumFine"];
    nStages = state["NumStages"];
    mCtrl = state["NumControls"];
    dtFine = state["DtFine"];
    vs = state["Vs"];
    vsH = state["VsH"];
    phi1Mats = state["Phi1Mats"];
    aTilMats = state["ATilMats"];
    g2Mats = state["G2Mats"];
    fs = state["Fs"];
    uFinal = state["FinalPropagator"];
    cDag = state["CDag"];
    hFlatT = state["ControlFlatT"];

    zN = cDag . uFinal;
    nTotFine = nFine * mCtrl;
    pList = Table[Null, {nTotFine}];
    dFList = Table[Null, {nTotFine}];
    dCDagList = Table[Null, {nTotFine}];
    hTilGrid = Table[Null, {nFine}, {mCtrl}];

    idx = 0;
    Do[
        vec = vs[[q]]; vecH = vsH[[q]];
        Do[
            idx++;
            hTil = vecH . hCtrls[[k]] . vec;
            hTilGrid[[q, k]] = hTil;
            dUqk = (-I * dtFine[[q]]) * (vec . (hTil * phi1Mats[[q]]) . vecH);
            pMat = ConjugateTranspose[fs[[q + 1]]] . dUqk . fs[[q]];
            dFMat = uFinal . pMat;
            pList[[idx]] = pMat;
            dFList[[idx]] = dFMat;
            dCDagList[[idx]] = objDirCotFn[uFinal, dFMat],
            {k, mCtrl}
        ],
        {q, nFine}
    ];

    matL = ToPackedComplexMatrix[Flatten[# . zN] & /@ pList];
    matR = ToPackedComplexMatrix[Transpose[Flatten[Transpose[#]] & /@ pList]];
    hProp = UpperTriangularize[Re[matL . matR], 1];
    Do[
        i0 = (q - 1) * mCtrl + 1;
        i1 = q * mCtrl;
        hProp[[i0 ;; i1, i0 ;; i1]] = 0.,
        {q, nFine}
    ];
    hProp = hProp + Transpose[hProp];

    Do[
        vec = vs[[q]]; vecH = vsH[[q]];
        diagBlk = Table[
            w2Mat = vec . contractG2[aTilMats[[q]], hTilGrid[[q, k]], g2Mats[[q]]] . vecH;
            -(dtFine[[q]]^2) * Re[hFlatT . Flatten[w2Mat]],
            {k, mCtrl}
        ];
        i0 = (q - 1) * mCtrl + 1;
        i1 = q * mCtrl;
        hProp[[i0 ;; i1, i0 ;; i1]] = 0.5 * (diagBlk + Transpose[diagBlk]),
        {q, nFine}
    ];

    dotCDagMat = ToPackedComplexMatrix[Flatten /@ dCDagList];
    dFTransMat = ToPackedComplexMatrix[Transpose[Flatten[Transpose[#]] & /@ dFList]];
    hTerm = Re[dotCDagMat . dFTransMat];

    hFine = Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ (hProp + hTerm)];
    If[nStages === 1,
        hFine,
        hCoarse = ArrayReshape[
            Total[Total[ArrayReshape[hFine, {nt, nStages, mCtrl, nt, nStages, mCtrl}], {2}], {4}],
            {nt * mCtrl, nt * mCtrl}
        ];
        Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ hCoarse]
    ]
];

End[];
EndPackage[];
