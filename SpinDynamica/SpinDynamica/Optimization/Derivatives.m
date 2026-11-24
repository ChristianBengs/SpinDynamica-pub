(* ::Package:: *)

BeginPackage["Optimization`Derivatives`",
{"Global`",
 "SpinDynamica`",
 "OptimalControl`Kernels`"
}];

ClearAll[
    EvaluateMemberObjective,
    EvaluateMemberObjectiveAndGrad,
    EvaluateMemberWithCache,
    EvaluateMemberHVP,
    EvaluateProblemOracles,
    BuildFiniteDiffOracles
];

EvaluateMemberObjective::usage = "EvaluateMemberObjective[member, xVec, cacheEnv] evaluates a single compiled member's transformation amplitude with dependency-aware stage caching.";
EvaluateMemberObjectiveAndGrad::usage = "EvaluateMemberObjectiveAndGrad[member, xVec, cacheEnv] evaluates objective and gradient for a compiled member using analytical adjoint (Hilbert, CF4, or Liouville).";
EvaluateMemberWithCache::usage = "EvaluateMemberWithCache[member, xVec, cacheEnv] evaluates objective, gradient, and second-order cached spectral state for a compiled Hilbert member.";
EvaluateMemberHVP::usage = "EvaluateMemberHVP[member, evalState, vVec] computes the exact analytical Hessian-vector product for a compiled Hilbert member in O(K) operations.";
EvaluateProblemOracles::usage = "EvaluateProblemOracles[comp, xVec, mode, opts] evaluates an ensemble or single compiled problem across Strategies A-E (None, FiniteDifference, Analytical, SecondOrder).";
BuildFiniteDiffOracles::usage = "BuildFiniteDiffOracles[comp, evalObjFn, opts] constructs adaptive parameter-scaled finite-difference gradient and optional Gauss-Newton/FD oracles.";

Begin["`Private`"];

(* Pack helpers delegating directly to exported OptimalControl`Kernels` utilities *)
toPackedComplexMat[m_] := ToPackedComplexMatrix[m];
toPackedRealMat[m_] := ToPackedRealMatrix[m];

(* Safe complex scalar trace Tr[ConjugateTranspose[a] . b] *)
innerProductComplex[aMat_, bMat_] := Total[Conjugate[aMat] * bMat, 2];

(* Build 2n x 2n block matrix and extract upper-right Frechet block of MatrixExp[{{A, E}, {0, A}}] *)
frechetExpmBlock[aMat_, eMat_] := Module[{n = Length[aMat], blk, expBlk},
    If[!MatrixQ[aMat, NumericQ] || !MatrixQ[eMat, NumericQ], Return[$Failed]];
    blk = ArrayFlatten[{{aMat, eMat}, {ConstantArray[0. + 0. I, {n, n}], aMat}}];
    expBlk = MatrixExp[blk];
    If[!MatrixQ[expBlk, NumericQ], Return[$Failed]];
    toPackedComplexMat[expBlk[[1 ;; n, n + 1 ;; 2 n]]]
];

(* Hermitize eigendecomposition and spectral phi1 matrix for U = Exp[-I G], G = Sum c_r M_r (Hermitian) *)
spectralHilbertStage[gMat_] := Module[
    {gHerm, vals, vecs, vecH, phi1, uMat},
    gHerm = toPackedComplexMat[0.5 * (gMat + ConjugateTranspose[gMat])];
    If[!MatrixQ[gHerm, NumericQ], Return[$Failed]];
    {vals, vecs} = Eigensystem[gHerm];
    vals = Re[vals];
    vecs = toPackedComplexMat[Transpose[vecs]];
    (* Orthonormalize in case of degenerate eigenvalues *)
    If[Max[Abs[ConjugateTranspose[vecs] . vecs - IdentityMatrix[Length[gHerm]]]] > 1.*^-10,
        vecs = toPackedComplexMat[Transpose[Orthogonalize[Transpose[vecs], Conjugate[#1] . #2 &]]]
    ];
    vecH = toPackedComplexMat[ConjugateTranspose[vecs]];
    phi1 = toPackedComplexMat[
        Table[
            With[{d = vals[[r]] - vals[[c]]},
                If[Abs[d] < 1.*^-7,
                    Exp[-I * vals[[c]]] * (1. - 0.5 I * d - (d^2) / 6. + (I * d^3) / 24.),
                    (Exp[-I * vals[[r]]] - Exp[-I * vals[[c]]]) / (-I * d)
                ]
            ],
            {r, Length[vals]}, {c, Length[vals]}
        ]
    ];
    uMat = toPackedComplexMat[vecs . (Exp[-I * vals] * vecH)];
    <|"U" -> uMat, "V" -> vecs, "VH" -> vecH, "Evals" -> Developer`ToPackedArray[vals], "Phi1" -> phi1|>
];

buildHilbertG2Clean[vals_List, phi1Mat_] := Module[{n = Length[vals], g2, tol = 1.*^-6, a, b, c},
    g2 = ConstantArray[0. + 0. I, {n, n, n}];
    Do[
        a = vals[[i]]; b = vals[[k]]; c = vals[[j]];
        g2[[i, k, j]] = Which[
            Abs[a - c] >= tol,
                (phi1Mat[[i, k]] - phi1Mat[[k, j]]) / (-I * (a - c)),
            Abs[b - c] >= tol,
                (phi1Mat[[k, i]] - phi1Mat[[i, j]]) / (-I * (b - c)),
            Abs[a - b] >= tol,
                (phi1Mat[[i, j]] - phi1Mat[[j, k]]) / (-I * (a - b)),
            True,
                With[{m = (a + b + c) / 3., da = -I * (a - (a + b + c)/3.), db = -I * (b - (a + b + c)/3.), dc = -I * (c - (a + b + c)/3.)},
                    0.5 * Exp[-I * m] * (1. + (da + db + dc)/3. + (da^2 + db^2 + dc^2 + da db + db dc + da dc)/12.)
                ]
        ],
        {i, n}, {k, n}, {j, n}
    ];
    Developer`ToPackedArray[g2]
];

(* Evaluate stage propagator with per-stage dependency-slice cache *)
getStageData[stage_Association, xVec_List, cacheEnv_, needSpectral_?BooleanQ] := Module[
    {kind = stage["Kind"], key, cols, slice, cached, gMats, coeffs, gMat, spec, prop, useHT},
    If[kind === "Static" || kind === "LiouvilleStatic",
        Return[<|"Propagator" -> stage["Propagator"]|>]
    ];
    key = stage["ID"];
    cols = stage["DepCols"];
    slice = If[cols === {}, {}, xVec[[cols]]];
    useHT = Head[cacheEnv] === DataStructure;
    If[useHT && cacheEnv["KeyExistsQ", key],
        cached = cacheEnv["Lookup", key];
        If[cached["Slice"] === slice && (!needSpectral || KeyExistsQ[cached, "Spectral"]),
            If[cacheEnv["KeyExistsQ", "__Hits"], cacheEnv["Insert", "__Hits" -> (cacheEnv["Lookup", "__Hits"] + 1)]];
            Return[cached]
        ]
    ];
    If[useHT && cacheEnv["KeyExistsQ", "__Misses"],
        cacheEnv["Insert", "__Misses" -> (cacheEnv["Lookup", "__Misses"] + 1)]
    ];
    Which[
        kind === "HilbertExpm",
            coeffs = stage["EvalCoeffs"][xVec];
            gMats = stage["Matrices"];
            gMat = toPackedComplexMat[Sum[coeffs[[r]] * gMats[[r]], {r, Length[gMats]}]];
            If[!MatrixQ[gMat, NumericQ], Return[$Failed]];
            If[needSpectral,
                spec = spectralHilbertStage[gMat];
                If[spec === $Failed, Return[$Failed]];
                cached = <|"Slice" -> slice, "Coeffs" -> coeffs, "G" -> gMat,
                           "Propagator" -> spec["U"], "Spectral" -> spec|>,
                prop = toPackedComplexMat[MatrixExp[-I * gMat]];
                cached = <|"Slice" -> slice, "Coeffs" -> coeffs, "G" -> gMat, "Propagator" -> prop|>
            ],
        kind === "HilbertCF4",
            cached = evalCF4Stage[stage, xVec, needSpectral];
            If[cached === $Failed, Return[$Failed]];
            cached = Append[cached, "Slice" -> slice],
        kind === "LiouvilleExpm",
            coeffs = stage["EvalCoeffs"][xVec];
            gMats = stage["Matrices"];
            gMat = toPackedComplexMat[Sum[coeffs[[r]] * gMats[[r]], {r, Length[gMats]}]];
            If[!MatrixQ[gMat, NumericQ], Return[$Failed]];
            prop = toPackedComplexMat[MatrixExp[gMat]];
            cached = <|"Slice" -> slice, "Coeffs" -> coeffs, "A" -> gMat, "Propagator" -> prop|>,
        kind === "LiouvilleCF4",
            cached = evalLiouvilleCF4Stage[stage, xVec];
            If[cached === $Failed, Return[$Failed]];
            cached = Append[cached, "Slice" -> slice]
    ];
    If[cached === $Failed || !AssociationQ[cached], Return[$Failed]];
    If[useHT, cacheEnv["Insert", key -> cached]];
    cached
];

(* CF4 time-dependent stage evaluation: stores all 2*nSub substeps cleanly *)
evalCF4Stage[stage_Association, xVec_List, needSpectral_?BooleanQ] := Module[
    {tStart, dur, nSub, h, c1, c2, a1, a2, subList, k, t1, t2, mats, coef1, coef2, g1, g2,
     s1, s2, uAcc, d = stage["Dim"]},
    tStart = stage["EvalTStart"][xVec];
    dur = stage["EvalDuration"][xVec];
    nSub = stage["NumSubsteps"];
    mats = stage["Matrices"];
    If[!NumericQ[tStart] || !NumericQ[dur], Return[$Failed]];
    If[Abs[N[dur]] < 1.*^-15,
        Return[<|"TStart" -> tStart, "Duration" -> dur, "SubStages" -> {},
                 "Propagator" -> toPackedComplexMat[IdentityMatrix[d]]|>]
    ];
    h = N[dur] / nSub;
    c1 = 0.5 - Sqrt[3.] / 6.; c2 = 0.5 + Sqrt[3.] / 6.;
    a1 = 0.25 - Sqrt[3.] / 6.; a2 = 0.25 + Sqrt[3.] / 6.;
    subList = Table[
        t1 = tStart + (k - 1 + c1) * h;
        t2 = tStart + (k - 1 + c2) * h;
        coef1 = stage["EvalHamCoeffs"][xVec, t1];
        coef2 = stage["EvalHamCoeffs"][xVec, t2];
        If[!VectorQ[coef1, NumericQ] || !VectorQ[coef2, NumericQ], Return[$Failed, Module]];
        (* Earlier factor first (uFirst uses a2*H1 + a1*H2), later factor second (uSecond uses a1*H1 + a2*H2) *)
        With[{cFirst = h * (a2 * coef1 + a1 * coef2), cSecond = h * (a1 * coef1 + a2 * coef2)},
            g1 = toPackedComplexMat[Sum[cFirst[[r]] * mats[[r]], {r, Length[mats]}]];
            g2 = toPackedComplexMat[Sum[cSecond[[r]] * mats[[r]], {r, Length[mats]}]];
            If[!MatrixQ[g1, NumericQ] || !MatrixQ[g2, NumericQ], Return[$Failed, Module]];
            {s1, s2} = If[needSpectral,
                {spectralHilbertStage[g1], spectralHilbertStage[g2]},
                {<|"U" -> toPackedComplexMat[MatrixExp[-I * g1]]|>, <|"U" -> toPackedComplexMat[MatrixExp[-I * g2]]|>}
            ];
            If[s1 === $Failed || s2 === $Failed, Return[$Failed, Module]];
            {{t1, t2, h, cFirst, g1, s1}, {t1, t2, h, cSecond, g2, s2}}
        ],
        {k, nSub}
    ];
    If[subList === $Failed || MemberQ[subList, $Failed], Return[$Failed]];
    subList = Flatten[subList, 1];
    uAcc = toPackedComplexMat[IdentityMatrix[d]];
    Do[uAcc = toPackedComplexMat[subList[[m, 6, "U"]] . uAcc], {m, Length[subList]}];
    <|"TStart" -> tStart, "Duration" -> dur, "SubStages" -> subList, "Propagator" -> uAcc|>
];

evalLiouvilleCF4Stage[stage_Association, xVec_List] := Module[
    {tStart, dur, nSub, h, c1, c2, a1, a2, subList, k, t1, t2, mats, coef1, coef2, aMat1, aMat2,
     p1, p2, pAcc, d = stage["Dim"]},
    tStart = stage["EvalTStart"][xVec];
    dur = stage["EvalDuration"][xVec];
    nSub = stage["NumSubsteps"];
    mats = stage["Matrices"];
    If[Abs[N[dur]] < 1.*^-15,
        Return[<|"TStart" -> tStart, "Duration" -> dur, "SubStages" -> {},
                 "Propagator" -> toPackedComplexMat[IdentityMatrix[d]]|>]
    ];
    h = N[dur] / nSub;
    c1 = 0.5 - Sqrt[3.] / 6.; c2 = 0.5 + Sqrt[3.] / 6.;
    a1 = 0.25 - Sqrt[3.] / 6.; a2 = 0.25 + Sqrt[3.] / 6.;
    subList = Table[
        t1 = tStart + (k - 1 + c1) * h;
        t2 = tStart + (k - 1 + c2) * h;
        coef1 = stage["EvalHamCoeffs"][xVec, t1];
        coef2 = stage["EvalHamCoeffs"][xVec, t2];
        If[!VectorQ[coef1, NumericQ] || !VectorQ[coef2, NumericQ], Return[$Failed, Module]];
        With[{cFirst = h * (a2 * coef1 + a1 * coef2), cSecond = h * (a1 * coef1 + a2 * coef2)},
            aMat1 = toPackedComplexMat[Sum[cFirst[[r]] * mats[[r]], {r, Length[mats]}]];
            aMat2 = toPackedComplexMat[Sum[cSecond[[r]] * mats[[r]], {r, Length[mats]}]];
            If[!MatrixQ[aMat1, NumericQ] || !MatrixQ[aMat2, NumericQ], Return[$Failed, Module]];
            p1 = toPackedComplexMat[MatrixExp[aMat1]];
            p2 = toPackedComplexMat[MatrixExp[aMat2]];
            {{t1, t2, h, cFirst, aMat1, p1}, {t1, t2, h, cSecond, aMat2, p2}}
        ],
        {k, nSub}
    ];
    If[subList === $Failed || MemberQ[subList, $Failed], Return[$Failed]];
    subList = Flatten[subList, 1];
    pAcc = toPackedComplexMat[IdentityMatrix[d]];
    Do[pAcc = toPackedComplexMat[subList[[m, 6]] . pAcc], {m, Length[subList]}];
    <|"TStart" -> tStart, "Duration" -> dur, "SubStages" -> subList, "Propagator" -> pAcc|>
];

(* Assemble operator / vector / scalar from Operator-Atom spec *)
evalMatrixSpec[spec_Association, xVec_List] := If[spec["IsConstant"],
    spec["Constant"],
    With[{c = spec["EvalCoeffs"][xVec], mats = spec["Matrices"]},
        toPackedComplexMat[Sum[c[[r]] * mats[[r]], {r, Length[mats]}]]
    ]
];

evalMatrixSpecDeriv[spec_Association, xVec_List, a_Integer] := Module[{jMat, mats = spec["Matrices"]},
    If[spec["IsConstant"] || !spec["DepMask"][[a]], Return[None]];
    jMat = spec["EvalCoeffJac"][xVec];
    With[{pos = FirstPosition[spec["DepCols"], a, Missing[]]},
        If[MissingQ[pos], None,
            With[{colIdx = First[pos]},
                toPackedComplexMat[Sum[jMat[[colIdx, r]] * mats[[r]], {r, Length[mats]}]]
            ]
        ]
    ]
];

evalVectorSpec[spec_Association, xVec_List] := If[spec["IsConstant"],
    spec["Constant"],
    With[{c = spec["EvalCoeffs"][xVec], vecs = spec["Vectors"]},
        Developer`ToPackedArray[N[Total[c * vecs]] + 0. I]
    ]
];

(* =============================================================== *)
(* Strategy A: Objective Only Evaluation                           *)
(* =============================================================== *)
EvaluateMemberObjective[member_Association, xVec_List, cacheEnv_: None] := Module[
    {space = member["Space"], stages = member["Stages"], rho, obs, normVal, sData, j, nS, denom, rawAmp, v},
    If[space === "BlackBox",
        Return[N[Re[member["EvalBlackBox"][xVec]]]]
    ];
    nS = Length[stages];
    normVal = If[member["NormSpec"]["IsConstant"], member["NormSpec"]["Constant"], member["NormSpec"]["Eval"][xVec]];
    If[!NumericQ[normVal] || Abs[N[normVal]] == 0., Return[$Failed]];

    If[space === "Hilbert",
        rho = evalMatrixSpec[member["RhoSpec"], xVec];
        Do[
            sData = getStageData[stages[[j]], xVec, cacheEnv, False];
            If[sData === $Failed, Return[$Failed, Module]];
            With[{u = sData["Propagator"]},
                rho = toPackedComplexMat[u . rho . ConjugateTranspose[u]]
            ],
            {j, nS}
        ];
        obs = evalMatrixSpec[member["ObsSpec"], xVec];
        denom = normVal * Re[innerProductComplex[obs, obs]];
        If[denom == 0., Return[$Failed]];
        rawAmp = innerProductComplex[obs, rho] / denom;
        Return[N[Re[rawAmp]]]
    ];

    If[space === "Liouville",
        v = If[member["RhoSpec"]["IsConstant"],
            member["RhoSpec"]["Constant"],
            With[{c = member["RhoSpec"]["EvalCoeffs"][xVec], vs = member["RhoSpec"]["Vectors"]},
                Developer`ToPackedArray[Sum[c[[r]] * vs[[r]], {r, Length[vs]}]]
            ]
        ];
        Do[
            sData = getStageData[stages[[j]], xVec, cacheEnv, False];
            If[sData === $Failed, Return[$Failed, Module]];
            v = sData["Propagator"] . v,
            {j, nS}
        ];
        obs = If[member["ObsSpec"]["IsConstant"],
            member["ObsSpec"]["Constant"],
            With[{c = member["ObsSpec"]["EvalCoeffs"][xVec], vs = member["ObsSpec"]["Vectors"]},
                Developer`ToPackedArray[Sum[c[[r]] * vs[[r]], {r, Length[vs]}]]
            ]
        ];
        With[{obsNormSq = If[member["ObsSpec"]["UseOperatorNorm"], member["ObsSpec"]["ConstantNormSq"], Re[Conjugate[obs] . obs]]},
            denom = normVal * obsNormSq;
            If[denom == 0., Return[$Failed]];
            rawAmp = (Conjugate[obs] . v) / denom;
            Return[N[Re[rawAmp]]]
        ]
    ];
    $Failed
];

(* =============================================================== *)
(* Strategy C/D: Generalised Analytical Adjoint (Hilbert/CF4/Liouv)*)
(* =============================================================== *)
EvaluateMemberObjectiveAndGrad[member_Association, xVec_List, cacheEnv_: None] := Module[
    {space = member["Space"], nP = member["NumParams"], stages = member["Stages"], nS,
     rho0, obs, normVal, obsNormSq, denom, num, rawAmp, grad, rhos, sDatas, lambdas, lMat,
     j, r, a, st, sd, bTot, vec, vecH, phi1, w0, jac, mats, colA,
     vList, wList, v0, wVec, dObs, dNorm, dDenom, rhoK, aMat, eMat, dU, dRho0, dObsA,
     subList, m, nSub, subData, dJac, kSub, dH1, dH2, dtStart, dDur, dC1, dC2, t1, t2, dCFirst, dCSecond},
    If[space === "BlackBox" || !TrueQ[member["IsDifferentiable"]],
        Return[$Failed]
    ];
    nS = Length[stages];
    grad = ConstantArray[0., nP];
    normVal = If[member["NormSpec"]["IsConstant"], member["NormSpec"]["Constant"], member["NormSpec"]["Eval"][xVec]];
    If[!NumericQ[normVal] || Abs[N[normVal]] == 0., Return[$Failed]];

    If[space === "Hilbert",
        rho0 = evalMatrixSpec[member["RhoSpec"], xVec];
        rhos = Table[Null, {nS + 1}];
        sDatas = Table[Null, {nS}];
        rhos[[1]] = rho0;
        Do[
            sDatas[[j]] = getStageData[stages[[j]], xVec, cacheEnv, True];
            If[sDatas[[j]] === $Failed, Return[$Failed, Module]];
            With[{u = sDatas[[j]]["Propagator"]},
                rhos[[j + 1]] = toPackedComplexMat[u . rhos[[j]] . ConjugateTranspose[u]]
            ],
            {j, nS}
        ];
        rhoK = rhos[[nS + 1]];
        obs = evalMatrixSpec[member["ObsSpec"], xVec];
        obsNormSq = Re[innerProductComplex[obs, obs]];
        denom = N[Re[normVal]] * obsNormSq;
        If[denom == 0., Return[$Failed]];
        num = innerProductComplex[obs, rhoK];
        rawAmp = N[Re[num / denom]];

        (* Terminal adjoint Lambda_K = O / D (so Re Tr[Lambda_K^H rhoK] = rawAmp) *)
        lMat = toPackedComplexMat[obs / denom];
        lambdas = Table[Null, {nS + 1}];
        lambdas[[nS + 1]] = lMat;

        (* Backward adjoint sweep *)
        For[j = nS, j >= 1, j--,
            st = stages[[j]];
            sd = sDatas[[j]];
            With[{u = sd["Propagator"]},
                lMat = toPackedComplexMat[ConjugateTranspose[u] . lMat . u];
                lambdas[[j]] = lMat
            ];
            If[!AnyTrue[st["DepMask"], TrueQ], Continue[]];

            Which[
                st["Kind"] === "HilbertExpm",
                    vec = sd["Spectral"]["V"];
                    vecH = sd["Spectral"]["VH"];
                    phi1 = sd["Spectral"]["Phi1"];
                    (* General non-Hermitian exact Frechet adjoint contraction:
                       Re Tr[Lambda_{j}^H dU rho_{j-1} U^H + Lambda_j^H U rho_{j-1} dU^H]
                       = Re Tr[ (-I V (Phi1^T * (V^H B_tot V)) V^H) dG ]? CAREFUL WITH TRANSPOSE/CONJUGATE!
                       Let us derive cleanly for dU = -I V (Phi1 * (V^H dG V)) V^H:
                       Note: Tr[A (V (Phi1 * (V^H dG V)) V^H)] = Tr[ (V^H A V) (Phi1 * (V^H dG V)) ]
                       = Sum_{p,q} (V^H A V)_{q,p} Phi1_{p,q} (V^H dG V)_{p,q}
                       = Tr[ (Phi1^T * (V^H A V))^T ??? NO! Notice Transpose indexing:
                       Sum_{p,q} M_{q,p} Phi1_{p,q} X_{p,q} = Sum_{p,q} (Phi1^T * M)_{q,p} X_{p,q} = Tr[(Phi1^T * M) . X]!!
                       YES!! Tr[P . X] = Sum_{p,q} P_{q,p} X_{p,q}, so if M = V^H A V, then P_{q,p} = Phi1_{p,q} M_{q,p} = (Transpose[Phi1] * M)_{q,p}!!
                       Therefore: Tr[A . (-I V (Phi1 * (V^H dG V)) V^H)] = Tr[ W_A . dG ] where W_A = -I V . (Transpose[Phi1] * (V^H A V)) . V^H!
                       Now apply to Re( Tr[L^H dU R U^H] + Tr[L^H U R dU^H] ):
                       Since dG is Hermitian (dG^H = dG), Tr[M dU^H] = Conjugate[Tr[dU M^H]] = Conjugate[Tr[M^H dU]].
                       So Re Tr[L^H dU R U^H + L^H U R dU^H] = Re Tr[ (R U^H L^H + R^H U^H L) . dU ]!
                       Here A = R . U^H . L^H + R^H . U^H . L (notice when R=R^H and L=L^H, A = R U^H L + R U^H L = 2 R U^H L!).
                       And Tr[W_A . dG] for dG Hermitian equals Sum_r (dg_r) Tr[W_A . M_r]!
                       CHECK: when dG is Hermitian and we take Re[Tr[W_A . dG]], we get Re Tr[W_A . M_r]! *)
                    bTot = toPackedComplexMat[
                        rhos[[j]] . ConjugateTranspose[sd["Propagator"]] . ConjugateTranspose[lambdas[[j + 1]]] +
                        ConjugateTranspose[rhos[[j]]] . ConjugateTranspose[sd["Propagator"]] . lambdas[[j + 1]]
                    ];
                    w0 = toPackedComplexMat[-I * vec . (Transpose[phi1] * (vecH . bTot . vec)) . vecH];
                    jac = st["EvalCoeffJac"][xVec];
                    mats = st["Matrices"];
                    With[{traces = Re[Table[Total[Transpose[w0] * mats[[r]], 2], {r, Length[mats]}]]},
                        (* Note Total[Transpose[A] * B, 2] === Tr[A . B]! *)
                        Do[
                            colA = st["DepCols"][[idx]];
                            grad[[colA]] += jac[[idx]] . traces,
                            {idx, Length[st["DepCols"]]}
                        ]
                    ],

                st["Kind"] === "HilbertCF4",
                    (* Propagate substeps forward then adjoint backward through all 2*nSub sub-stages! *)
                    subList = sd["SubStages"];
                    If[Length[subList] > 0,
                        Module[{subRhos, rCur, lSub, mTot = Length[subList], sItem, sp, wSub, jacT,
                                c1 = 0.5 - Sqrt[3.]/6., c2 = 0.5 + Sqrt[3.]/6.,
                                a1 = 0.25 - Sqrt[3.]/6., a2 = 0.25 + Sqrt[3.]/6.,
                                h = sd["Duration"] / st["NumSubsteps"], dur = sd["Duration"], nSub = st["NumSubsteps"]},
                            subRhos = Table[Null, {mTot + 1}];
                            rCur = rhos[[j]];
                            subRhos[[1]] = rCur;
                            Do[
                                rCur = toPackedComplexMat[subList[[m, 6, "U"]] . rCur . ConjugateTranspose[subList[[m, 6, "U"]]]];
                                subRhos[[m + 1]] = rCur,
                                {m, mTot}
                            ];
                            lSub = lambdas[[j + 1]];
                            mats = st["Matrices"];
                            For[m = mTot, m >= 1, m--,
                                sItem = subList[[m]];
                                sp = sItem[[6]];
                                bTot = toPackedComplexMat[
                                    subRhos[[m]] . ConjugateTranspose[sp["U"]] . ConjugateTranspose[lSub] +
                                    ConjugateTranspose[subRhos[[m]]] . ConjugateTranspose[sp["U"]] . lSub
                                ];
                                wSub = toPackedComplexMat[-I * sp["V"] . (Transpose[sp["Phi1"]] * (sp["VH"] . bTot . sp["V"])) . sp["VH"]];
                                With[{traces = Re[Table[Total[Transpose[wSub] * mats[[r]], 2], {r, Length[mats]}]],
                                      kIdx = Ceiling[m / 2.],
                                      isSecond = EvenQ[m],
                                      tNode1 = sItem[[1]], tNode2 = sItem[[2]]},
                                    Module[{wNode1, wNode2, wt1, wt2, dT1, dT2, jH1, jH2, dH1dt, dH2dt, dDurVec, dTStartVec},
                                        {wNode1, wNode2} = If[!isSecond, {a2, a1}, {a1, a2}];
                                        jH1 = st["EvalHamParamJac"][xVec, tNode1]; (* pAct x R *)
                                        jH2 = st["EvalHamParamJac"][xVec, tNode2]; (* pAct x R *)
                                        dH1dt = st["EvalHamTimeDeriv"][xVec, tNode1]; (* length R *)
                                        dH2dt = st["EvalHamTimeDeriv"][xVec, tNode2]; (* length R *)
                                        dDurVec = st["EvalDurationGrad"][xVec]; (* length pAct *)
                                        dTStartVec = st["EvalTStartGrad"][xVec]; (* length pAct *)
                                        (* Chain rule for C_m = (dur/nSub) * (wNode1 * H(tNode1, x) + wNode2 * H(tNode2, x)):
                                           d(C_m)/dx_idx = (dDur_idx / dur) * C_m
                                             + h * wNode1 * (jH1[[idx]] + dH1dt * (dTStart_idx + ((kIdx - 1 + c1)/nSub)*dDur_idx))
                                             + h * wNode2 * (jH2[[idx]] + dH2dt * (dTStart_idx + ((kIdx - 1 + c2)/nSub)*dDur_idx)) *)
                                        With[{cM = sItem[[4]], f1 = (kIdx - 1 + c1) / nSub, f2 = (kIdx - 1 + c2) / nSub},
                                            Do[
                                                colA = st["DepCols"][[idx]];
                                                dT1 = dTStartVec[[idx]] + f1 * dDurVec[[idx]];
                                                dT2 = dTStartVec[[idx]] + f2 * dDurVec[[idx]];
                                                With[{dCm = (dDurVec[[idx]] / dur) * cM +
                                                        h * (wNode1 * (jH1[[idx]] + dT1 * dH1dt) + wNode2 * (jH2[[idx]] + dT2 * dH2dt))},
                                                    grad[[colA]] += dCm . traces
                                                ],
                                                {idx, Length[st["DepCols"]]}
                                            ]
                                        ]
                                    ]
                                ];
                                lSub = toPackedComplexMat[ConjugateTranspose[sp["U"]] . lSub . sp["U"]]
                            ]
                        ]
                    ]
            ]
        ];

        (* Boundary terms: parameter-dependent RhoIni, Obs, and NormFactor *)
        If[!member["RhoSpec"]["IsConstant"],
            Do[
                If[member["RhoSpec"]["DepMask"][[a]],
                    dRho0 = evalMatrixSpecDeriv[member["RhoSpec"], xVec, a];
                    If[dRho0 =!= None,
                        grad[[a]] += Re[innerProductComplex[lambdas[[1]], dRho0]]
                    ]
                ],
                {a, nP}
            ]
        ];
        If[!member["ObsSpec"]["IsConstant"] || !member["NormSpec"]["IsConstant"],
            With[{dNormGrad = If[member["NormSpec"]["IsConstant"], ConstantArray[0., nP], member["NormSpec"]["EvalGrad"][xVec]]},
                Do[
                    dObsA = If[!member["ObsSpec"]["IsConstant"] && member["ObsSpec"]["DepMask"][[a]],
                        evalMatrixSpecDeriv[member["ObsSpec"], xVec, a],
                        None
                    ];
                    If[dObsA =!= None || dNormGrad[[a]] != 0.,
                        dNum = If[dObsA =!= None, Re[innerProductComplex[dObsA, rhoK]], 0.];
                        dDenom = dNormGrad[[a]] * obsNormSq +
                            N[Re[normVal]] * If[dObsA =!= None, 2. * Re[innerProductComplex[dObsA, obs]], 0.];
                        grad[[a]] += dNum / denom - (N[Re[num]] * dDenom) / (denom^2)
                    ],
                    {a, nP}
                ]
            ]
        ];

        Return[<|"Objective" -> rawAmp, "Gradient" -> Developer`ToPackedArray[grad], "Extra" -> <||>|>]
    ];

    If[space === "Liouville",
        v0 = If[member["RhoSpec"]["IsConstant"],
            member["RhoSpec"]["Constant"],
            With[{c = member["RhoSpec"]["EvalCoeffs"][xVec], vs = member["RhoSpec"]["Vectors"]},
                Developer`ToPackedArray[Sum[c[[r]] * vs[[r]], {r, Length[vs]}]]
            ]
        ];
        vList = Table[Null, {nS + 1}];
        sDatas = Table[Null, {nS}];
        vList[[1]] = v0;
        Do[
            sDatas[[j]] = getStageData[stages[[j]], xVec, cacheEnv, False];
            If[sDatas[[j]] === $Failed, Return[$Failed, Module]];
            vList[[j + 1]] = Developer`ToPackedArray[sDatas[[j]]["Propagator"] . vList[[j]]],
            {j, nS}
        ];
        obs = If[member["ObsSpec"]["IsConstant"],
            member["ObsSpec"]["Constant"],
            With[{c = member["ObsSpec"]["EvalCoeffs"][xVec], vs = member["ObsSpec"]["Vectors"]},
                Developer`ToPackedArray[Sum[c[[r]] * vs[[r]], {r, Length[vs]}]]
            ]
        ];
        obsNormSq = If[member["ObsSpec"]["UseOperatorNorm"], member["ObsSpec"]["ConstantNormSq"], Re[Conjugate[obs] . obs]];
        denom = N[Re[normVal]] * obsNormSq;
        If[denom == 0., Return[$Failed]];
        num = Conjugate[obs] . vList[[nS + 1]];
        rawAmp = N[Re[num / denom]];

        wVec = Developer`ToPackedArray[obs / denom];
        wList = Table[Null, {nS + 1}];
        wList[[nS + 1]] = wVec;

        For[j = nS, j >= 1, j--,
            st = stages[[j]];
            sd = sDatas[[j]];
            wVec = Developer`ToPackedArray[ConjugateTranspose[sd["Propagator"]] . wVec];
            wList[[j]] = wVec;
            If[!AnyTrue[st["DepMask"], TrueQ], Continue[]];

            Which[
                st["Kind"] === "LiouvilleExpm",
                    jac = st["EvalCoeffJac"][xVec];
                    mats = st["Matrices"];
                    aMat = sd["A"];
                    Do[
                        colA = st["DepCols"][[idx]];
                        eMat = toPackedComplexMat[Sum[jac[[idx, r]] * mats[[r]], {r, Length[mats]}]];
                        (* Fast shortcut when eMat == alpha * aMat (exact for Delay[tau] with static L!) *)
                        With[{dProp = If[Length[mats] === 1 && AllTrue[Flatten[Abs[aMat]], # == 0 &],
                                ConstantArray[0. + 0. I, Dimensions[aMat]],
                                If[TrueQ[st["IsScalarTimeGenerator"]] && Length[mats] === 1,
                                    (jac[[idx, 1]] * mats[[1]]) . sd["Propagator"],
                                    frechetExpmBlock[aMat, eMat]
                                ]
                            ]},
                            If[dProp === $Failed, Return[$Failed, Module]];
                            grad[[colA]] += Re[Conjugate[wList[[j + 1]]] . (dProp . vList[[j]])]
                        ],
                        {idx, Length[st["DepCols"]]}
                    ],

                st["Kind"] === "LiouvilleCF4",
                    subList = sd["SubStages"];
                    If[Length[subList] > 0,
                        Module[{subVs, vCur, wSub, mTot = Length[subList], sItem,
                                c1 = 0.5 - Sqrt[3.]/6., c2 = 0.5 + Sqrt[3.]/6.,
                                a1 = 0.25 - Sqrt[3.]/6., a2 = 0.25 + Sqrt[3.]/6.,
                                h = sd["Duration"] / st["NumSubsteps"], dur = sd["Duration"], nSub = st["NumSubsteps"]},
                            subVs = Table[Null, {mTot + 1}];
                            vCur = vList[[j]];
                            subVs[[1]] = vCur;
                            Do[vCur = Developer`ToPackedArray[subList[[m, 6]] . vCur]; subVs[[m + 1]] = vCur, {m, mTot}];
                            wSub = wList[[j + 1]];
                            mats = st["Matrices"];
                            For[m = mTot, m >= 1, m--,
                                sItem = subList[[m]];
                                With[{kIdx = Ceiling[m / 2.], isSecond = EvenQ[m], tNode1 = sItem[[1]], tNode2 = sItem[[2]]},
                                    Module[{wNode1, wNode2, jH1, jH2, dH1dt, dH2dt, dDurVec, dTStartVec, cM, f1, f2, dT1, dT2, dCm, eSub, dProp},
                                        {wNode1, wNode2} = If[!isSecond, {a2, a1}, {a1, a2}];
                                        jH1 = st["EvalHamParamJac"][xVec, tNode1];
                                        jH2 = st["EvalHamParamJac"][xVec, tNode2];
                                        dH1dt = st["EvalHamTimeDeriv"][xVec, tNode1];
                                        dH2dt = st["EvalHamTimeDeriv"][xVec, tNode2];
                                        dDurVec = st["EvalDurationGrad"][xVec];
                                        dTStartVec = st["EvalTStartGrad"][xVec];
                                        cM = sItem[[4]]; f1 = (kIdx - 1 + c1) / nSub; f2 = (kIdx - 1 + c2) / nSub;
                                        Do[
                                            colA = st["DepCols"][[idx]];
                                            dT1 = dTStartVec[[idx]] + f1 * dDurVec[[idx]];
                                            dT2 = dTStartVec[[idx]] + f2 * dDurVec[[idx]];
                                            dCm = (dDurVec[[idx]] / dur) * cM +
                                                h * (wNode1 * (jH1[[idx]] + dT1 * dH1dt) + wNode2 * (jH2[[idx]] + dT2 * dH2dt));
                                            eSub = toPackedComplexMat[Sum[dCm[[r]] * mats[[r]], {r, Length[mats]}]];
                                            dProp = frechetExpmBlock[sItem[[5]], eSub];
                                            If[dProp === $Failed, Return[$Failed, Module]];
                                            grad[[colA]] += Re[Conjugate[wSub] . (dProp . subVs[[m]])],
                                            {idx, Length[st["DepCols"]]}
                                        ]
                                    ]
                                ];
                                wSub = Developer`ToPackedArray[ConjugateTranspose[sItem[[6]]] . wSub]
                            ]
                        ]
                    ]
            ]
        ];

        (* Boundary terms in Liouville space *)
        If[!member["RhoSpec"]["IsConstant"],
            With[{jR = member["RhoSpec"]["EvalCoeffJac"][xVec], vecs = member["RhoSpec"]["Vectors"]},
                Do[
                    colA = member["RhoSpec"]["DepCols"][[idx]];
                    With[{dV0 = Sum[jR[[idx, r]] * vecs[[r]], {r, Length[vecs]}]},
                        grad[[colA]] += Re[Conjugate[wList[[1]]] . dV0]
                    ],
                    {idx, Length[member["RhoSpec"]["DepCols"]]}
                ]
            ]
        ];
        If[!member["ObsSpec"]["IsConstant"] || !member["NormSpec"]["IsConstant"],
            With[{dNormGrad = If[member["NormSpec"]["IsConstant"], ConstantArray[0., nP], member["NormSpec"]["EvalGrad"][xVec]],
                  jO = If[member["ObsSpec"]["IsConstant"], None, member["ObsSpec"]["EvalCoeffJac"][xVec]],
                  vecsO = If[member["ObsSpec"]["IsConstant"], {}, member["ObsSpec"]["Vectors"]]},
                Do[
                    dObs = If[!member["ObsSpec"]["IsConstant"] && member["ObsSpec"]["DepMask"][[a]],
                        With[{pos = FirstPosition[member["ObsSpec"]["DepCols"], a, Missing[]]},
                            If[MissingQ[pos], None, Sum[jO[[First[pos], r]] * vecsO[[r]], {r, Length[vecsO]}]]
                        ],
                        None
                    ];
                    If[dObs =!= None || dNormGrad[[a]] != 0.,
                        With[{dNum = If[dObs =!= None, Re[Conjugate[dObs] . vList[[nS + 1]]], 0.],
                              dDen = dNormGrad[[a]] * obsNormSq +
                                  N[Re[normVal]] * If[dObs =!= None && !member["ObsSpec"]["UseOperatorNorm"],
                                      2. * Re[Conjugate[dObs] . obs], 0.]},
                            grad[[a]] += dNum / denom - (N[Re[num]] * dDen) / (denom^2)
                        ]
                    ],
                    {a, nP}
                ]
            ]
        ];

        Return[<|"Objective" -> rawAmp, "Gradient" -> Developer`ToPackedArray[grad], "Extra" -> <||>|>]
    ];

    $Failed
];

(* =============================================================== *)
(* Strategy E: Cached Forward-Adjoint State & Exact O(K) Hilbert HVP*)
(* =============================================================== *)
(* Expand CF4 stage on-the-fly into equivalent HilbertExpm sub-stages so Strategy E covers BOTH *)
flattenMemberStagesForHVP[member_Association, xVec_List, cacheEnv_: None] := Module[
    {stages = member["Stages"], flat = {}, j, st, sd, k, subList, sItem, mTot,
     c1 = 0.5 - Sqrt[3.]/6., c2 = 0.5 + Sqrt[3.]/6.,
     a1 = 0.25 - Sqrt[3.]/6., a2 = 0.25 + Sqrt[3.]/6.,
     nSub, dur, h, mats, nP = member["NumParams"]},
    Do[
        st = stages[[j]];
        sd = getStageData[st, xVec, cacheEnv, True];
        If[sd === $Failed, Return[$Failed, Module]];
        If[st["Kind"] =!= "HilbertCF4",
            AppendTo[flat, {st, sd}],
            subList = sd["SubStages"];
            mTot = Length[subList];
            dur = sd["Duration"]; nSub = st["NumSubsteps"]; h = dur / nSub; mats = st["Matrices"];
            Do[
                sItem = subList[[m]];
                With[{kIdx = Ceiling[m / 2.], isSecond = EvenQ[m], tNode1 = sItem[[1]], tNode2 = sItem[[2]],
                      depCols = st["DepCols"], depMask = st["DepMask"]},
                    Module[{w1 = If[!isSecond, a2, a1], w2 = If[!isSecond, a1, a2],
                            f1, f2, jH1, jH2, dH1dt, dH2dt, dDur, dTStart, jacMat, hessFn, subSt, subSd},
                        f1 = (kIdx - 1 + c1) / nSub; f2 = (kIdx - 1 + c2) / nSub;
                        jH1 = st["EvalHamParamJac"][xVec, tNode1];
                        jH2 = st["EvalHamParamJac"][xVec, tNode2];
                        dH1dt = st["EvalHamTimeDeriv"][xVec, tNode1];
                        dH2dt = st["EvalHamTimeDeriv"][xVec, tNode2];
                        dDur = st["EvalDurationGrad"][xVec];
                        dTStart = st["EvalTStartGrad"][xVec];
                        jacMat = Table[
                            (dDur[[idx]] / dur) * sItem[[4]] +
                            h * (w1 * (jH1[[idx]] + (dTStart[[idx]] + f1 * dDur[[idx]]) * dH1dt) +
                                 w2 * (jH2[[idx]] + (dTStart[[idx]] + f2 * dDur[[idx]]) * dH2dt)),
                            {idx, Length[depCols]}
                        ];
                        hessFn = If[!KeyExistsQ[st, "EvalSubstepHessianFn"] || st["EvalSubstepHessianFn"] === None,
                            None,
                            st["EvalSubstepHessianFn"][m, kIdx, isSecond, w1, w2, f1, f2, xVec]
                        ];
                        subSt = <|
                            "Kind" -> "HilbertExpm",
                            "DepMask" -> depMask,
                            "DepCols" -> depCols,
                            "Matrices" -> mats,
                            "PreJac" -> jacMat,
                            "EvalCoeffHVP" -> hessFn
                        |>;
                        subSd = <|
                            "Propagator" -> sItem[[6, "U"]],
                            "Spectral" -> sItem[[6]]
                        |>;
                        AppendTo[flat, {subSt, subSd}]
                    ]
                ],
                {m, mTot}
            ]
        ],
        {j, Length[stages]}
    ];
    flat
];

EvaluateMemberWithCache[member_Association, xVec_List, cacheEnv_: None] := Module[
    {space = member["Space"], nP = member["NumParams"], flat, nS,
     rho0, rhos, rhoK, obs, normVal, obsNormSq, denom, num, rawAmp, lMat, lambdas,
     grad, j, st, sd, bTot, vec, vecH, phi1, w0, jac, mats, colA, g2s, w0s, bTots, jacs, evalState,
     dRho0, dObsA, dNormGrad, dNum, dDenom},
    If[space =!= "Hilbert" || !TrueQ[member["IsDifferentiable"]], Return[$Failed]];
    normVal = If[member["NormSpec"]["IsConstant"], member["NormSpec"]["Constant"], member["NormSpec"]["Eval"][xVec]];
    If[!NumericQ[normVal] || Abs[N[normVal]] == 0., Return[$Failed]];

    flat = flattenMemberStagesForHVP[member, xVec, cacheEnv];
    If[flat === $Failed, Return[$Failed]];
    nS = Length[flat];
    grad = ConstantArray[0., nP];
    rho0 = evalMatrixSpec[member["RhoSpec"], xVec];
    rhos = Table[Null, {nS + 1}];
    rhos[[1]] = rho0;
    Do[
        With[{u = flat[[j, 2]]["Propagator"]},
            rhos[[j + 1]] = toPackedComplexMat[u . rhos[[j]] . ConjugateTranspose[u]]
        ],
        {j, nS}
    ];
    rhoK = rhos[[nS + 1]];
    obs = evalMatrixSpec[member["ObsSpec"], xVec];
    obsNormSq = Re[innerProductComplex[obs, obs]];
    denom = N[Re[normVal]] * obsNormSq;
    If[denom == 0., Return[$Failed]];
    num = innerProductComplex[obs, rhoK];
    rawAmp = N[Re[num / denom]];

    lMat = toPackedComplexMat[obs / denom];
    lambdas = Table[Null, {nS + 1}];
    lambdas[[nS + 1]] = lMat;
    g2s = Table[Null, {nS}];
    w0s = Table[Null, {nS}];
    bTots = Table[Null, {nS}];
    jacs = Table[Null, {nS}];

    For[j = nS, j >= 1, j--,
        {st, sd} = flat[[j]];
        With[{u = sd["Propagator"]},
            lMat = toPackedComplexMat[ConjugateTranspose[u] . lMat . u];
            lambdas[[j]] = lMat
        ];
        If[!AnyTrue[st["DepMask"], TrueQ], Continue[]];
        vec = sd["Spectral"]["V"];
        vecH = sd["Spectral"]["VH"];
        phi1 = sd["Spectral"]["Phi1"];
        bTot = toPackedComplexMat[
            rhos[[j]] . ConjugateTranspose[sd["Propagator"]] . ConjugateTranspose[lambdas[[j + 1]]] +
            ConjugateTranspose[rhos[[j]]] . ConjugateTranspose[sd["Propagator"]] . lambdas[[j + 1]]
        ];
        bTots[[j]] = bTot;
        w0 = toPackedComplexMat[-I * vec . (Transpose[phi1] * (vecH . bTot . vec)) . vecH];
        w0s[[j]] = w0;
        jac = If[KeyExistsQ[st, "PreJac"], st["PreJac"], st["EvalCoeffJac"][xVec]];
        jacs[[j]] = jac;
        mats = st["Matrices"];
        With[{traces = Re[Table[Total[Transpose[w0] * mats[[r]], 2], {r, Length[mats]}]]},
            Do[
                colA = st["DepCols"][[idx]];
                grad[[colA]] += jac[[idx]] . traces,
                {idx, Length[st["DepCols"]]}
            ]
        ];
        (* Lazily or immediately cache G2 tensor only for active stages *)
        g2s[[j]] = buildHilbertG2Clean[sd["Spectral"]["Evals"], phi1]
    ];

    If[!member["RhoSpec"]["IsConstant"],
        Do[
            If[member["RhoSpec"]["DepMask"][[a]],
                dRho0 = evalMatrixSpecDeriv[member["RhoSpec"], xVec, a];
                If[dRho0 =!= None, grad[[a]] += Re[innerProductComplex[lambdas[[1]], dRho0]]]
            ],
            {a, nP}
        ]
    ];
    If[!member["ObsSpec"]["IsConstant"] || !member["NormSpec"]["IsConstant"],
        dNormGrad = If[member["NormSpec"]["IsConstant"], ConstantArray[0., nP], member["NormSpec"]["EvalGrad"][xVec]];
        Do[
            dObsA = If[!member["ObsSpec"]["IsConstant"] && member["ObsSpec"]["DepMask"][[a]],
                evalMatrixSpecDeriv[member["ObsSpec"], xVec, a],
                None
            ];
            If[dObsA =!= None || dNormGrad[[a]] != 0.,
                dNum = If[dObsA =!= None, Re[innerProductComplex[dObsA, rhoK]], 0.];
                dDenom = dNormGrad[[a]] * obsNormSq +
                    N[Re[normVal]] * If[dObsA =!= None, 2. * Re[innerProductComplex[dObsA, obs]], 0.];
                grad[[a]] += dNum / denom - (N[Re[num]] * dDenom) / (denom^2)
            ],
            {a, nP}
        ]
    ];

    evalState = <|
        "Space" -> "Hilbert",
        "X" -> xVec,
        "FlatStages" -> flat,
        "Rhos" -> rhos,
        "Lambdas" -> lambdas,
        "Obs" -> obs,
        "NormVal" -> normVal,
        "ObsNormSq" -> obsNormSq,
        "Denom" -> denom,
        "Num" -> num,
        "G2s" -> g2s,
        "W0s" -> w0s,
        "BTots" -> bTots,
        "Jacs" -> jacs
    |>;
    <|"Objective" -> rawAmp, "Gradient" -> Developer`ToPackedArray[grad], "EvalState" -> evalState, "Extra" -> <||>|>
];

(* Helper for directional 2nd derivative of f(O, norm) term: D_v(Lambda_K) and boundary Hessians *)
evalBoundaryHVP[member_Association, st_Association, vVec_List, dRhoK_] := Module[
    {xVec = st["X"], nP = member["NumParams"], rho0 = st["Rhos"][[1]], rhoK = Last[st["Rhos"]],
     obs = st["Obs"], normVal = st["NormVal"], obsNormSq = st["ObsNormSq"], denom = st["Denom"], num = st["Num"],
     l0 = st["Lambdas"][[1]], dObsV, dNormV, dDenomV, dNumV, dLK, gradHVP, a, dRhoA, dRho2, dObsA, dObs2, dNormG, dNormH,
     trOO, dDa, d2Dav, trOR, dNa, d2Nav, rhoSpec = member["RhoSpec"], obsSpec = member["ObsSpec"], normSpec = member["NormSpec"]},
    gradHVP = ConstantArray[0., nP];
    dNormG = If[normSpec["IsConstant"], ConstantArray[0., nP], normSpec["EvalGrad"][xVec]];
    dNormH = If[normSpec["IsConstant"], None, normSpec["EvalHess"][xVec]];
    dNormV = dNormG . vVec;

    dObsV = If[obsSpec["IsConstant"], None,
        With[{z = Total[Table[With[{dA = evalMatrixSpecDeriv[obsSpec, xVec, a]}, If[dA === None, 0, vVec[[a]] * dA]], {a, nP}]]},
            If[z === 0, None, toPackedComplexMat[z]]
        ]
    ];

    (* First directional derivative of D = norm * Re Tr[O^H O] *)
    dDenomV = dNormV * obsNormSq + If[dObsV =!= None, N[Re[normVal]] * 2. * Re[innerProductComplex[dObsV, obs]], 0.];
    (* Note: L_K = O / D, so d_v(L_K) = d_v(O)/D - (d_v D / D^2) O.
       And in Re[Num / Denom]:
       Exact derivative wrt theta_a of A = Re[N / D] along variation rhoK' + v:
       d^2 A / (d theta_a d v) = [Re Tr[d_a L_K^H dRhoK] + terminal boundary terms].
       Let's separate cleanly into:
       (1) dLK = directional derivative of L_K with respect to x_v *exclusively* (propagated backward via dLambdas to accumulate with rho_j & dRho_0!)
           Look! At final state rho_K:
           g_term(theta_a) = Re Tr[ (d_a O)^H rho_K ] / D - (Re Num) * d_a D / D^2 + Re Tr[ (O/D)^H d_a rho_K ]!
           Differentiating Re Tr[ (O/D)^H d_a rho_K ] in direction v produces TWO parts:
             (i) Re Tr[ (d_v(O/D))^H d_a rho_K ]  <-- THIS is just dLK propagated backward through adjoint dLambda_{j-1} = U_j^H dLambda_j U_j!
             (ii) Re Tr[ (O/D)^H d_v(d_a rho_K) ] <-- THIS comes from stage second derivatives + tangent dRho_j!
           AND differentiating [Re Tr[(d_a O)^H rho_K] / D - (Re Num) * d_a D / D^2] wrt BOTH v -> dRhoK AND v -> (d_v O, d_v norm) gives exact obs/norm boundary Hessians! *)
    dLK = If[dObsV =!= None, dObsV / denom, ConstantArray[0. + 0. I, Dimensions[obs]]] -
        (dDenomV / denom^2) * obs;
    dLK = toPackedComplexMat[dLK];

    (* Initial state second derivative: D_v D_a rho0 . Lambda_0 + D_a rho0 . dLambda_0 (done in main HVP once dLambda_0 computed) *)
    If[!rhoSpec["IsConstant"] && KeyExistsQ[rhoSpec, "EvalCoeffHVP"] && rhoSpec["EvalCoeffHVP"] =!= None,
        With[{cHvp = rhoSpec["EvalCoeffHVP"][xVec, vVec[[rhoSpec["DepCols"]]]], matsR = rhoSpec["Matrices"]},
            Do[
                With[{colA = rhoSpec["DepCols"][[idx]], d2M = toPackedComplexMat[Sum[cHvp[[idx, r]] * matsR[[r]], {r, Length[matsR]}]]},
                    gradHVP[[colA]] += Re[innerProductComplex[l0, d2M]]
                ],
                {idx, Length[rhoSpec["DepCols"]]}
            ]
        ]
    ];

    (* Obs & Norm second derivatives for each parameter a *)
    If[!obsSpec["IsConstant"] || !normSpec["IsConstant"],
        With[{obsHvpMat = If[!obsSpec["IsConstant"] && KeyExistsQ[obsSpec, "EvalCoeffHVP"] && obsSpec["EvalCoeffHVP"] =!= None,
                obsSpec["EvalCoeffHVP"][xVec, vVec[[obsSpec["DepCols"]]]],
                None
            ]},
            trOO = obsNormSq;
            trOR = Re[num];
            dNumV = Re[innerProductComplex[obs, dRhoK]] + If[dObsV =!= None, Re[innerProductComplex[dObsV, rhoK]], 0.];
            Do[
                dObsA = If[!obsSpec["IsConstant"] && obsSpec["DepMask"][[a]], evalMatrixSpecDeriv[obsSpec, xVec, a], None];
                If[dObsA =!= None || dNormG[[a]] != 0. || (dNormH =!= None && TrueQ[Max[Abs[dNormH[[a]]]] > 0.]),
                    dNa = If[dObsA =!= None, Re[innerProductComplex[dObsA, rhoK]], 0.];
                    dDa = dNormG[[a]] * trOO + N[Re[normVal]] * If[dObsA =!= None, 2. * Re[innerProductComplex[dObsA, obs]], 0.];
                    (* d2(Obs)/dx_a dv *)
                    dObs2 = If[obsHvpMat =!= None && obsSpec["DepMask"][[a]],
                        With[{pos = FirstPosition[obsSpec["DepCols"], a, Missing[]]},
                            If[MissingQ[pos], None,
                                With[{cRow = obsHvpMat[[First[pos]]], matsO = obsSpec["Matrices"]},
                                    toPackedComplexMat[Sum[cRow[[r]] * matsO[[r]], {r, Length[matsO]}]]
                                ]
                            ]
                        ],
                        None
                    ];
                d2Nav = If[dObsA =!= None, Re[innerProductComplex[dObsA, dRhoK]], 0.] +
                        If[dObs2 =!= None, Re[innerProductComplex[dObs2, rhoK]], 0.];
                d2Dav = If[dNormH =!= None, dNormH[[a]] . vVec, 0.] * trOO +
                    dNormG[[a]] * If[dObsV =!= None, 2. * Re[innerProductComplex[dObsV, obs]], 0.] +
                    dNormV * If[dObsA =!= None, 2. * Re[innerProductComplex[dObsA, obs]], 0.] +
                    N[Re[normVal]] * (
                        2. * If[dObs2 =!= None, Re[innerProductComplex[dObs2, obs]], 0.] +
                        2. * If[dObsA =!= None && dObsV =!= None, Re[innerProductComplex[dObsA, dObsV]], 0.]
                    );
                (* Quotient rule 2nd derivative contribution to gradient[a] = Re(d_a N)/D - Re(N) d_a D / D^2:
                   minus (d_v Tr[L_K^H d_a rho_K]) which is handled by dLK backward sweep! *)
                gradHVP[[a]] += d2Nav / denom - (dNa * dDenomV + dNumV * dDa + trOR * d2Dav) / (denom^2) +
                    2. * trOR * dDa * dDenomV / (denom^3)
            ],
            {a, nP}
        ]
        ]
    ];

    {dLK, gradHVP}
];

EvaluateMemberHVP[member_Association, evalState_Association, vIn_List] := Module[
    {nP = member["NumParams"], vVec, flat, nS, rhos, lambdas, g2s, w0s, bTots, jacs,
     xVec, dRhos, dUs, dGs, dRhoCur, dU, dG, st, sd, sp, vec, vecH, phi1, jac, mats,
     j, r, idx, colA, dLK, gradHVP, dLambdas, dLCur, dBTot, bTot, w1, w2, vTil, dGTil, w3,
     dC2Mat, h1, traces1, traces2, traces3, g2, b2, mR},
    If[evalState === $Failed || !AssociationQ[evalState] || evalState["Space"] =!= "Hilbert",
        Return[$Failed]
    ];
    vVec = Developer`ToPackedArray[N[Re[Flatten[vIn]]]];
    xVec = evalState["X"];
    flat = evalState["FlatStages"];
    nS = Length[flat];
    rhos = evalState["Rhos"];
    lambdas = evalState["Lambdas"];
    g2s = evalState["G2s"];
    w0s = evalState["W0s"];
    bTots = evalState["BTots"];
    jacs = evalState["Jacs"];

    dRhos = Table[Null, {nS + 1}];
    dUs = Table[Null, {nS}];
    dGs = Table[Null, {nS}];

    (* Initial tangent rho0' = Sum_a v_a d_a rho0 *)
    dRhoCur = If[member["RhoSpec"]["IsConstant"],
        ConstantArray[0. + 0. I, Dimensions[rhos[[1]]]],
        Sum[
            With[{dA = evalMatrixSpecDeriv[member["RhoSpec"], xVec, a]},
                If[dA === None, 0, vVec[[a]] * dA]
            ],
            {a, nP}
        ]
    ];
    dRhoCur = toPackedComplexMat[dRhoCur];
    dRhos[[1]] = dRhoCur;

    (* Tangent forward pass *)
    Do[
        {st, sd} = flat[[j]];
        With[{u = sd["Propagator"]},
            If[!AnyTrue[st["DepMask"], TrueQ] || AllTrue[vVec[[st["DepCols"]]], # == 0. &],
                dGs[[j]] = None;
                dUs[[j]] = None;
                dRhoCur = toPackedComplexMat[u . dRhoCur . ConjugateTranspose[u]],
                sp = sd["Spectral"];
                vec = sp["V"]; vecH = sp["VH"]; phi1 = sp["Phi1"];
                jac = jacs[[j]];
                mats = st["Matrices"];
                With[{vStage = vVec[[st["DepCols"]]] . jac},
                    dG = toPackedComplexMat[Sum[vStage[[r]] * mats[[r]], {r, Length[mats]}]];
                    dU = toPackedComplexMat[-I * vec . (phi1 * (vecH . dG . vec)) . vecH];
                    dGs[[j]] = dG;
                    dUs[[j]] = dU;
                    dRhoCur = toPackedComplexMat[
                        u . dRhoCur . ConjugateTranspose[u] +
                        dU . rhos[[j]] . ConjugateTranspose[u] +
                        u . rhos[[j]] . ConjugateTranspose[dU]
                    ]
                ]
            ]
        ];
        dRhos[[j + 1]] = dRhoCur,
        {j, nS}
    ];

    {dLK, gradHVP} = evalBoundaryHVP[member, evalState, vVec, dRhos[[nS + 1]]];
    dLCur = dLK;

    (* Tangent backward pass + second-order stage contractions *)
    For[j = nS, j >= 1, j--,
        {st, sd} = flat[[j]];
        With[{u = sd["Propagator"], lNext = lambdas[[j + 1]]},
            (* Note: bTot = R U^H L^H + R^H U^H L. Directional derivative of bTot in direction v involves dR, dU^H, and dL! *)
            If[AnyTrue[st["DepMask"], TrueQ],
                sp = sd["Spectral"];
                vec = sp["V"]; vecH = sp["VH"]; phi1 = sp["Phi1"];
                mats = st["Matrices"];
                jac = jacs[[j]];
                mR = Length[mats];

                (* (A) d_v(bTot) = dR U^H L^H + R (dU)^H L^H + R U^H (dL)^H + conj-transpose swaps *)
                dU = dUs[[j]];
                With[{rMat = rhos[[j]], dRMat = dRhos[[j]], uH = ConjugateTranspose[u], lH = ConjugateTranspose[lNext]},
                    dBTot = dRMat . uH . lH + ConjugateTranspose[dRMat] . uH . lNext +
                            rMat . uH . ConjugateTranspose[dLCur] + ConjugateTranspose[rMat] . uH . dLCur;
                    If[dU =!= None,
                        dBTot += rMat . ConjugateTranspose[dU] . lH + ConjugateTranspose[rMat] . ConjugateTranspose[dU] . lNext
                    ];
                    dBTot = toPackedComplexMat[dBTot];
                    w1 = toPackedComplexMat[-I * vec . (Transpose[phi1] * (vecH . dBTot . vec)) . vecH];
                    traces1 = Re[Table[Total[Transpose[w1] * mats[[r]], 2], {r, mR}]]
                ];

                (* (B) Second Fréchet derivative of propagator D^2 Exp[-I G][-I G_v, -I M_r]:
                   Derived earlier: Tr[ B_tot . (-V [Sum_k G2[p,k,q] (V_til[p,k] M_til[k,q] + M_til[p,k] V_til[k,q])] V^H) ]
                   Let BT = vecH . bTot . vec, VT = vecH . dGs[[j]] . vec.
                   Contraction sum_{p,k,q} BT[q, p] * (-G2[p, k, q]) * (VT[p,k] MT[k,q] + MT[p,k] VT[k,q]):
                   Group terms multiplying MT[alpha, beta]:
                   - From k=alpha, q=beta: sum_p BT[beta, p] * (-G2[p, alpha, beta]) * VT[p, alpha]
                   - From p=alpha, k=beta: sum_q BT[q, alpha] * (-G2[alpha, beta, q]) * VT[beta, q]
                   Look how that equals a single matrix W2Til[beta, alpha] (so Tr[W2Til . MT] = Tr[W2 . M_r] with W2 = vec . W2Til . vecH)!
                   Notice G2[p, k, q] is completely symmetric in {p, k, q}!! Check if it matches contractG2!
                   Wait: let's write out the two contracted matrices C1[beta, alpha] and C2[beta, alpha] directly OR test in 1D!
                   For 1x1: BT[1,1]=b, VT[1,1]=v, G2[1,1,1]=0.5 e^{-i mu}:
                   W2Til[1,1] = -0.5 e^{-i mu} * (b v + b v) = -e^{-i mu} b v.
                   Check: d/dv( d/da Exp[-i G] ) = d/dv( -i e^{-i G} ) = -e^{-i G} * v!!
                   Tr[W2 . M_r] = (-e^{-i mu} b v) * M_r!! EXACT!!
                   Verify matrix order in sum_p: (Term 1)_{beta, alpha} = - Sum_p BT[beta, p] * G2[p, alpha, beta] * VT[p, alpha],
                   (Term 2)_{beta, alpha} = - Sum_q BT[q, alpha] * G2[alpha, beta, q] * VT[beta, q]. *)
                If[dGs[[j]] =!= None,
                    Module[{bTil = vecH . bTots[[j]] . vec, vTil = vecH . dGs[[j]] . vec, g2Mat = g2s[[j]], w2Til},
                        w2Til = -Table[
                            Sum[bTil[[beta, p]] * g2Mat[[p, alpha, beta]] * vTil[[p, alpha]], {p, Length[vec]}] +
                            Sum[bTil[[q, alpha]] * g2Mat[[alpha, beta, q]] * vTil[[beta, q]], {q, Length[vec]}],
                            {beta, Length[vec]}, {alpha, Length[vec]}
                        ];
                        w2 = toPackedComplexMat[vec . w2Til . vecH];
                        traces2 = Re[Table[Total[Transpose[w2] * mats[[r]], 2], {r, mR}]]
                    ],
                    traces2 = ConstantArray[0., mR]
                ];

                (* Combine (A) and (B) with first-order Jac *)
                traces3 = traces1 + traces2;
                Do[
                    colA = st["DepCols"][[idx]];
                    gradHVP[[colA]] += jac[[idx]] . traces3,
                    {idx, Length[st["DepCols"]]}
                ];

                (* (C) Second derivative of stage parameter mapping: Sum_b v_b d^2_{ab} g_{j,r} contracted with W0 *)
                If[KeyExistsQ[st, "EvalCoeffHVP"] && st["EvalCoeffHVP"] =!= None,
                    With[{cHvpMat = st["EvalCoeffHVP"][xVec, vVec[[st["DepCols"]]]],
                          trW0 = Re[Table[Total[Transpose[w0s[[j]]] * mats[[r]], 2], {r, mR}]]},
                        Do[
                            colA = st["DepCols"][[idx]];
                            gradHVP[[colA]] += cHvpMat[[idx]] . trW0,
                            {idx, Length[st["DepCols"]]}
                        ]
                    ]
                ]
            ];

            (* Propagate dLambda backward: dLambda_{j-1} = U^H dLambda_j U + dU^H Lambda_j U + U^H Lambda_j dU *)
            dU = dUs[[j]];
            dLCur = If[dU === None,
                toPackedComplexMat[ConjugateTranspose[u] . dLCur . u],
                toPackedComplexMat[
                    ConjugateTranspose[u] . dLCur . u +
                    ConjugateTranspose[dU] . lNext . u +
                    ConjugateTranspose[u] . lNext . dU
                ]
            ]
        ]
    ];

    (* Final initial-state cross term Re Tr[dLambda_0^H * d_a rho0] *)
    If[!member["RhoSpec"]["IsConstant"],
        Do[
            If[member["RhoSpec"]["DepMask"][[a]],
                With[{dA = evalMatrixSpecDeriv[member["RhoSpec"], xVec, a]},
                    If[dA =!= None, gradHVP[[a]] += Re[innerProductComplex[dLCur, dA]]]
                ]
            ],
            {a, nP}
        ]
    ];

    Developer`ToPackedArray[gradHVP]
];

(* =============================================================== *)
(* Ensemble Dispatcher & Finite Difference Strategy B              *)
(* =============================================================== *)
ensureSubkernelsReady[parallelSpec_, nMembers_Integer] := Module[{targetK, sdEntry},
    If[!((TrueQ[parallelSpec] || (IntegerQ[parallelSpec] && parallelSpec > 1)) && nMembers > 1), Return[False]];
    targetK = If[IntegerQ[parallelSpec], Min[parallelSpec, nMembers], Automatic];
    Quiet[
        If[$KernelCount === 0,
            If[targetK === Automatic, LaunchKernels[], LaunchKernels[targetK]],
            If[IntegerQ[targetK] && $KernelCount < targetK, LaunchKernels[targetK - $KernelCount]]
        ]
    ];
    sdEntry = If[StringQ[SpinDynamica`Private`$SpinDynamicaRootDirectory],
        FileNameJoin[{SpinDynamica`Private`$SpinDynamicaRootDirectory, "SpinDynamica.m"}],
        None
    ];
    If[StringQ[sdEntry] && FileExistsQ[sdEntry],
        With[{entry = sdEntry},
            ParallelEvaluate[
                If[!MemberQ[$Packages, "SpinDynamica`"], Needs["SpinDynamica`", entry]];
                Null
            ]
        ],
        ParallelNeeds["SpinDynamica`"]
    ];
    $KernelCount > 0
];

evalMemberListGrad[members_List, weights_List, xVec_List, cacheList_List, useParQ_] := Module[{res},
    res = If[TrueQ[useParQ],
        ParallelMap[EvaluateMemberObjectiveAndGrad[#, xVec, None] &, members],
        MapThread[EvaluateMemberObjectiveAndGrad[#1, xVec, #2] &, {members, cacheList}]
    ];
    If[MemberQ[res, $Failed] || !AllTrue[res, AssociationQ], Return[$Failed]];
    <|
        "Objective" -> N[weights . (#["Objective"] & /@ res)],
        "Gradient" -> Developer`ToPackedArray[Total[MapThread[#1 * #2["Gradient"] &, {weights, res}]]],
        "Extra" -> <|"MemberObjectives" -> (#["Objective"] & /@ res)|>
    |>
];

evalMemberListObj[members_List, weights_List, xVec_List, cacheList_List, useParQ_] := Module[{res},
    res = If[TrueQ[useParQ],
        ParallelMap[EvaluateMemberObjective[#, xVec, None] &, members],
        MapThread[EvaluateMemberObjective[#1, xVec, #2] &, {members, cacheList}]
    ];
    If[MemberQ[res, $Failed] || !VectorQ[res, NumericQ], Return[$Failed]];
    N[weights . res]
];

evalMemberListWithCache[members_List, weights_List, xVec_List, cacheList_List, useParQ_] := Module[{res},
    res = If[TrueQ[useParQ],
        ParallelMap[EvaluateMemberWithCache[#, xVec, None] &, members],
        MapThread[EvaluateMemberWithCache[#1, xVec, #2] &, {members, cacheList}]
    ];
    If[MemberQ[res, $Failed] || !AllTrue[res, AssociationQ], Return[$Failed]];
    <|
        "Objective" -> N[weights . (#["Objective"] & /@ res)],
        "Gradient" -> Developer`ToPackedArray[Total[MapThread[#1 * #2["Gradient"] &, {weights, res}]]],
        "EvalState" -> <|"MemberStates" -> (#["EvalState"] & /@ res)|>,
        "Extra" -> <|"MemberObjectives" -> (#["Objective"] & /@ res)|>
    |>
];

(* Strategy B: Parameter-scaled central or forward finite differences reusing stage cache *)
BuildFiniteDiffOracles[comp_Association, evalObjFn_, opts_List: {}] := Module[
    {nP = comp["NumParams"], bounds = comp["Bounds"], lb, ub, scaleVec, fdMode, epsRaw, epsBase, invS},
    {lb, ub} = bounds;
    scaleVec = Lookup[comp, "ScaleVec", ConstantArray[1., nP]];
    invS = 1. / scaleVec;
    fdMode = FirstCase[opts,
        (Rule | RuleDelayed)[k_ /; (k === "FiniteDifferenceOrder" || (Head[k] === Symbol && SymbolName[Unevaluated[k]] === "FiniteDifferenceOrder")), v_] :> v,
        "Central"
    ];
    epsRaw = FirstCase[opts,
        (Rule | RuleDelayed)[k_ /; (k === "FiniteDifferenceStep" || (Head[k] === Symbol && SymbolName[Unevaluated[k]] === "FiniteDifferenceStep")), v_] :> v,
        Automatic
    ];
    epsBase = If[!NumericQ[epsRaw] || epsRaw === Automatic,
        If[fdMode === "Forward", 1.5*^-8, 6.0*^-6],
        N[epsRaw]
    ];

    Function[{xIn},
        Module[{x = Developer`ToPackedArray[N[Re[xIn]]], f0, g = ConstantArray[0., nP], a, h,
                xP, xM, fP, fM, hasLo, hasHi},
            f0 = evalObjFn[x];
            If[!NumericQ[f0], Return[$Failed]];
            Do[
                h = epsBase * Max[Abs[x[[a]]], invS[[a]], 1.*^-12];
                hasLo = NumericQ[lb[[a]]] && Abs[N[lb[[a]]]] < Infinity;
                hasHi = NumericQ[ub[[a]]] && Abs[N[ub[[a]]]] < Infinity;
                xP = x; xM = x;
                Which[
                    fdMode === "Forward" || (hasLo && x[[a]] - h < lb[[a]]),
                        If[hasHi && x[[a]] + h > ub[[a]],
                            (* Shrink step to fit inside [lb, ub] *)
                            h = 0.25 * (ub[[a]] - lb[[a]])
                        ];
                        If[hasHi && x[[a]] + h > ub[[a]],
                            xM[[a]] = x[[a]] - h;
                            fM = evalObjFn[xM];
                            If[!NumericQ[fM], Return[$Failed]];
                            g[[a]] = (f0 - fM) / h,
                            xP[[a]] = x[[a]] + h;
                            fP = evalObjFn[xP];
                            If[!NumericQ[fP], Return[$Failed]];
                            g[[a]] = (fP - f0) / h
                        ],
                    hasHi && x[[a]] + h > ub[[a]],
                        xM[[a]] = x[[a]] - h;
                        fM = evalObjFn[xM];
                        If[!NumericQ[fM], Return[$Failed]];
                        g[[a]] = (f0 - fM) / h,
                    True,
                        xP[[a]] = x[[a]] + h;
                        xM[[a]] = x[[a]] - h;
                        fP = evalObjFn[xP];
                        fM = evalObjFn[xM];
                        If[!NumericQ[fP] || !NumericQ[fM], Return[$Failed]];
                        g[[a]] = (fP - fM) / (2. * h)
                ],
                {a, nP}
            ];
            <|"Objective" -> f0, "Gradient" -> Developer`ToPackedArray[g], "Extra" -> <|"DerivativeMethod" -> "FiniteDifference"|>|>
        ]
    ]
];

(* Public Strategy Dispatcher building full oracle set for RunGeneralOptimizer *)
EvaluateProblemOracles[comp_Association, derivMethod_: Automatic, opts_List: {}] := Module[
    {members = comp["Members"], weights = comp["Weights"], par = Lookup[comp, "Parallel", False],
     useParQ, caches, evalObjFn, evalGradFn, evalWithCacheFn, evalHVPFn, evalHessFn, canAnalytic, can2nd},
    useParQ = ensureSubkernelsReady[par, Length[members]];
    caches = If[TrueQ[Lookup[comp, "EnableCaching", True]],
        Table[
            With[{ht = CreateDataStructure["HashTable"]},
                ht["Insert", "__Hits" -> 0];
                ht["Insert", "__Misses" -> 0];
                ht
            ],
            {Length[members]}
        ],
        ConstantArray[None, Length[members]]
    ];

    evalObjFn = Function[{xIn},
        evalMemberListObj[members, weights, Developer`ToPackedArray[N[Re[xIn]]], caches, useParQ]
    ];

    canAnalytic = (comp["CompilationMode"] =!= "BlackBoxFallback") &&
        AllTrue[members, #["Space"] =!= "BlackBox" && TrueQ[#["IsDifferentiable"]] &];
    can2nd = canAnalytic && AllTrue[members, #["Space"] === "Hilbert" && TrueQ[#["SupportsSecondOrder"]] &];

    evalGradFn = Switch[derivMethod,
        "None", None,
        "FiniteDifference", BuildFiniteDiffOracles[comp, evalObjFn, opts],
        "Analytical" | "Adjoint" | "SecondOrder" | Automatic,
            If[canAnalytic,
                Function[{xIn}, evalMemberListGrad[members, weights, Developer`ToPackedArray[N[Re[xIn]]], caches, useParQ]],
                If[derivMethod === Automatic, BuildFiniteDiffOracles[comp, evalObjFn, opts], None]
            ],
        _, None
    ];

    evalWithCacheFn = If[can2nd && MemberQ[{"Analytical", "Adjoint", "SecondOrder", Automatic}, derivMethod],
        Function[{xIn}, evalMemberListWithCache[members, weights, Developer`ToPackedArray[N[Re[xIn]]], caches, useParQ]],
        None
    ];

    evalHVPFn = If[can2nd && MemberQ[{"Analytical", "Adjoint", "SecondOrder", Automatic}, derivMethod],
        Function[{evalState, vIn},
            Module[{mStates = evalState["MemberStates"], mHVPs},
                mHVPs = If[TrueQ[useParQ],
                    ParallelMap[EvaluateMemberHVP[#[[1]], #[[2]], vIn] &, Transpose[{members, mStates}]],
                    MapThread[EvaluateMemberHVP[#1, #2, vIn] &, {members, mStates}]
                ];
                If[MemberQ[mHVPs, $Failed] || !AllTrue[mHVPs, VectorQ[#, NumericQ] &],
                    $Failed,
                    Developer`ToPackedArray[Total[MapThread[#1 * #2 &, {weights, mHVPs}]]]
                ]
            ]
        ],
        None
    ];

    evalHessFn = If[evalHVPFn =!= None,
        Function[{evalState},
            Module[{nP = comp["NumParams"], id = IdentityMatrix[comp["NumParams"]], cols},
                cols = Table[evalHVPFn[evalState, id[[k]]], {k, nP}];
                If[MemberQ[cols, $Failed], $Failed,
                    With[{hRaw = Transpose[cols]},
                        toPackedRealMat[0.5 * (hRaw + Transpose[hRaw])]
                    ]
                ]
            ]
        ],
        None
    ];

    <|
        "EvalObj" -> evalObjFn,
        "EvalGrad" -> evalGradFn,
        "EvalWithCache" -> evalWithCacheFn,
        "EvalHVP" -> evalHVPFn,
        "EvalHessian" -> evalHessFn,
        "Caches" -> caches,
        "EffectiveDerivativeMethod" -> Which[
            derivMethod === "None", "None",
            derivMethod === "FiniteDifference", "FiniteDifference",
            can2nd && derivMethod === "SecondOrder", "SecondOrder",
            canAnalytic && MemberQ[{"Analytical", "Adjoint", "SecondOrder", Automatic}, derivMethod], "Analytical",
            evalGradFn =!= None, "FiniteDifference",
            True, "None"
        ],
        "SupportsSecondOrder" -> (evalHVPFn =!= None)
    |>
];

End[];
EndPackage[];
