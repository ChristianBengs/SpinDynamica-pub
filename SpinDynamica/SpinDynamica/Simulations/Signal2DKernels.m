(* ::Package:: *)

(* :Title: Simulations`Signal2DKernels` *)
(* :Context: Simulations`Signal2DKernels` *)
(* :Summary: Packed complex numerical kernels for high-performance 2D NMR simulation: secular block diagonalization, Hilbert/Liouville double diagonalization, stepped propagator marching, CF4 Magnus periodic MAS COMPUTE2D, and compiled ensemble accumulation. *)

BeginPackage["Simulations`Signal2DKernels`",
{"Global`",
 "SpinDynamica`",
 "General`UserLevel`",
 "General`Miscellaneous`",
 "QM`Hilbert`",
 "QM`SpinOperators`",
 "QM`Liouville`",
 "SD`Propagation`",
 "Signals`SignalExpansionKernels`"
}];

ClearAll[
    DiagonalizeHilbertHermitian,
    DiagonalizeLiouvilleGeneral,
    DiagonalizeLiouvilleBlocked,
    BuildHilbertCoherenceMask,
    ApplyHilbertMixSequence,
    EvaluateHilbertDoubleDiag2D,
    EvaluateLiouvilleDoubleDiag2D,
    EvaluateSteppedPropagator2D,
    CF4PropagatorStep,
    BuildPeriodicPropagatorsCF4,
    EvaluatePeriodicCOMPUTE2D,
    EvaluateCompiledEnsembleChunk2D
];

DiagonalizeHilbertHermitian::usage = "DiagonalizeHilbertHermitian[hMat] eigendecomposes a Hermitian Hilbert-space Hamiltonian matrix into packed Real64 eigenvalues and orthonormal Complex128 eigenvector matrices V and VH, exploiting diagonal fast paths and connected-component block structure.";
DiagonalizeLiouvilleGeneral::usage = "DiagonalizeLiouvilleGeneral[lMat] eigendecomposes a general Liouvillian superoperator matrix into packed Complex128 eigenvalues, X, and XInv, or returns $Failed if defective.";
DiagonalizeLiouvilleBlocked::usage = "DiagonalizeLiouvilleBlocked[lMat, activeIndices] decomposes a Liouvillian matrix into connected-component secular blocks, diagonalizing only blocks intersecting activeIndices (or All).";
BuildHilbertCoherenceMask::usage = "BuildHilbertCoherenceMask[filters, system, basis] constructs a packed coherence-order filtration mask in Hilbert space along with any needed basis transformation to ZeemanBasis.";
ApplyHilbertMixSequence::usage = "ApplyHilbertMixSequence[rhoMat, mixSteps] applies a sequence of unitary propagators and Hilbert coherence masks to a density matrix.";
EvaluateHilbertDoubleDiag2D::usage = "EvaluateHilbertDoubleDiag2D[spec] evaluates one or more 2D time-domain signals in Hilbert space using pre-diagonalized t1/t2 eigenbases and BLAS-3 direct-dimension evaluation.";
EvaluateLiouvilleDoubleDiag2D::usage = "EvaluateLiouvilleDoubleDiag2D[spec] evaluates one or more 2D time-domain signals in Liouville space via pre-diagonalized eigenbases and BLAS-3 contractions.";
EvaluateSteppedPropagator2D::usage = "EvaluateSteppedPropagator2D[spec] evaluates uniform-step 2D time-domain evolution by marching precomputed single-increment propagators across t1 and t2.";
CF4PropagatorStep::usage = "CF4PropagatorStep[evalGenFn, tStart, tEnd, nSub, space, dim] integrates time-dependent evolution from tStart to tEnd using 4th-order Gauss-Legendre commutator-free Magnus (CF4) steps.";
BuildPeriodicPropagatorsCF4::usage = "BuildPeriodicPropagatorsCF4[evalGenFn, t0, period, nDiv, nSubPerDiv, space, dim] precomputes cumulative one-period propagators {A0=I, A1, ..., An} across [t0, t0 + period].";
EvaluatePeriodicCOMPUTE2D::usage = "EvaluatePeriodicCOMPUTE2D[spec] evaluates periodic time-dependent / MAS 2D NMR signals using precomputed one-period Floquet propagators, BLAS-3 t2 contraction, and coherent rotor-phase-locked CarouselAverage.";
EvaluateCompiledEnsembleChunk2D::usage = "EvaluateCompiledEnsembleChunk2D[chunkPayload] evaluates a batch of compiled ensemble members using precomputed operator-atom matrices and accumulates directly into a single packed 2D buffer per component.";

Begin["`Private`"];

toPackedComplexMat[m_] := If[
    Developer`PackedArrayQ[m, Complex],
    m,
    Developer`ToPackedArray[N[Normal[m]] + 0. I]
];

toPackedComplexVec[v_] := If[
    Developer`PackedArrayQ[v, Complex],
    v,
    Developer`ToPackedArray[N[Normal[v]] + 0. I]
];

toPackedRealMat[m_] := If[
    Developer`PackedArrayQ[m, Real],
    m,
    Developer`ToPackedArray[N[Re[Normal[m]]]]
];

toPackedRealVec[v_] := If[
    Developer`PackedArrayQ[v, Real],
    v,
    Developer`ToPackedArray[N[Re[Normal[v]]]]
];

(* Sparsity-graph connected components for basis-agnostic secular block diagonalization *)
findSparsityBlocks[mat_?MatrixQ, tol_?NumericQ] := Module[
    {n = Length[mat], mArr, maxAbs, diagDiff, symMag, sa, adj, visited, blocks, comp, headIdx, cur},
    If[n <= 1, Return[If[n === 0, {}, {{1}}]]];
    mArr = Normal[mat];
    maxAbs = Max[Abs[mArr]];
    If[maxAbs == 0., Return[List /@ Range[n]]];
    diagDiff = Max[Abs[mArr - DiagonalMatrix[Diagonal[mArr]]]];
    If[diagDiff <= tol * maxAbs, Return[List /@ Range[n]]];
    If[n < 4, Return[{Range[n]}]];
    symMag = Abs[mArr] + Abs[Transpose[mArr]];
    sa = SparseArray[UnitStep[symMag - tol * maxAbs], Automatic, 0];
    adj = sa["AdjacencyLists"];
    visited = ConstantArray[False, n];
    blocks = {};
    Do[
        If[!visited[[start]],
            comp = {start};
            visited[[start]] = True;
            headIdx = 1;
            While[headIdx <= Length[comp],
                cur = comp[[headIdx++]];
                Do[
                    If[!visited[[nb]],
                        visited[[nb]] = True;
                        AppendTo[comp, nb]
                    ],
                    {nb, adj[[cur]]}
                ]
            ];
            AppendTo[blocks, Sort[comp]]
        ],
        {start, n}
    ];
    blocks
];

DiagonalizeHilbertHermitian[hMat_?MatrixQ] := Module[
    {gHerm, dim = Length[hMat], blocks, vals, v, vH, bMat, bVals, bVecs, bV, bVH, pos = 1},
    gHerm = toPackedComplexMat[0.5 * (hMat + ConjugateTranspose[hMat])];
    If[!MatrixQ[gHerm, NumericQ], Return[$Failed]];
    blocks = findSparsityBlocks[gHerm, 1.*^-14];
    If[Length[blocks] === dim,
        vals = toPackedRealVec[Re[Diagonal[gHerm]]];
        v = toPackedComplexMat[IdentityMatrix[dim]];
        Return[<|"Evals" -> vals, "V" -> v, "VH" -> v, "Dim" -> dim, "BlockCount" -> dim|>]
    ];
    If[Length[blocks] === 1,
        {bVals, bVecs} = Eigensystem[gHerm];
        vals = toPackedRealVec[Re[bVals]];
        v = toPackedComplexMat[Transpose[bVecs]];
        vH = toPackedComplexMat[ConjugateTranspose[v]];
        If[Max[Abs[vH . v - IdentityMatrix[dim]]] > 1.*^-10,
            v = toPackedComplexMat[Transpose[Orthogonalize[Transpose[v], Conjugate[#1] . #2 &]]];
            vH = toPackedComplexMat[ConjugateTranspose[v]]
        ];
        Return[<|"Evals" -> vals, "V" -> v, "VH" -> vH, "Dim" -> dim, "BlockCount" -> 1|>]
    ];
    vals = ConstantArray[0., dim];
    v = ConstantArray[0. + 0. I, {dim, dim}];
    Do[
        If[Length[blk] === 1,
            With[{idx = First[blk]},
                vals[[pos]] = N[Re[gHerm[[idx, idx]]]];
                v[[idx, pos]] = 1. + 0. I;
                pos++
            ],
            With[{k = Length[blk]},
                bMat = gHerm[[blk, blk]];
                {bVals, bVecs} = Eigensystem[bMat];
                bV = Transpose[bVecs];
                bVH = ConjugateTranspose[bV];
                If[Max[Abs[bVH . bV - IdentityMatrix[k]]] > 1.*^-10,
                    bV = Transpose[Orthogonalize[Transpose[bV], Conjugate[#1] . #2 &]]
                ];
                vals[[pos ;; pos + k - 1]] = N[Re[bVals]];
                v[[blk, pos ;; pos + k - 1]] = bV;
                pos += k
            ]
        ],
        {blk, blocks}
    ];
    vals = toPackedRealVec[vals];
    v = toPackedComplexMat[v];
    vH = toPackedComplexMat[ConjugateTranspose[v]];
    <|"Evals" -> vals, "V" -> v, "VH" -> vH, "Dim" -> dim, "BlockCount" -> Length[blocks]|>
];

DiagonalizeLiouvilleBlocked[lMat_?MatrixQ, activeIndices_] := Module[
    {g, dim = Length[lMat], blocks, neededBlocks, subDim, evals, x, xInv,
     pos = 1, bVals, bVecs, bX, bXInv, failQ = False},
    g = toPackedComplexMat[lMat];
    If[!MatrixQ[g, NumericQ], Return[$Failed]];
    blocks = findSparsityBlocks[g, 1.*^-14];
    neededBlocks = If[activeIndices === All,
        blocks,
        With[{actSet = Union[Select[Flatten[{activeIndices}], IntegerQ]]},
            Select[blocks, IntersectingQ[#, actSet] &]
        ]
    ];
    subDim = Total[Length /@ neededBlocks];
    If[subDim === 0,
        Return[<|
            "Evals" -> toPackedComplexVec[{}],
            "X" -> toPackedComplexMat[ConstantArray[0. + 0. I, {dim, 0}]],
            "XInv" -> toPackedComplexMat[ConstantArray[0. + 0. I, {0, dim}]],
            "Dim" -> dim,
            "SubspaceDim" -> 0,
            "BlockCount" -> Length[blocks],
            "ActiveBlockCount" -> 0
        |>]
    ];
    If[Length[neededBlocks] === 1 && subDim === dim,
        {bVals, bVecs} = Quiet[Check[Eigensystem[g], {$Failed, $Failed}]];
        If[bVals === $Failed, Return[$Failed]];
        x = toPackedComplexMat[Transpose[bVecs]];
        xInv = Quiet[Check[LinearSolve[x, IdentityMatrix[dim]], $Failed]];
        If[!MatrixQ[xInv, NumericQ] || Max[Abs[xInv . x - IdentityMatrix[dim]]] > 1.*^-6, Return[$Failed]];
        Return[<|
            "Evals" -> toPackedComplexVec[bVals],
            "X" -> x,
            "XInv" -> toPackedComplexMat[xInv],
            "Dim" -> dim,
            "SubspaceDim" -> dim,
            "BlockCount" -> Length[blocks],
            "ActiveBlockCount" -> 1
        |>]
    ];
    evals = ConstantArray[0. + 0. I, subDim];
    x = ConstantArray[0. + 0. I, {dim, subDim}];
    xInv = ConstantArray[0. + 0. I, {subDim, dim}];
    Do[
        If[Length[blk] === 1,
            With[{idx = First[blk]},
                evals[[pos]] = g[[idx, idx]];
                x[[idx, pos]] = 1. + 0. I;
                xInv[[pos, idx]] = 1. + 0. I;
                pos++
            ],
            With[{k = Length[blk]},
                {bVals, bVecs} = Quiet[Check[Eigensystem[g[[blk, blk]]], {$Failed, $Failed}]];
                If[bVals === $Failed, failQ = True; Break[]];
                bX = Transpose[bVecs];
                bXInv = Quiet[Check[LinearSolve[bX, IdentityMatrix[k]], $Failed]];
                If[!MatrixQ[bXInv, NumericQ] || Max[Abs[bXInv . bX - IdentityMatrix[k]]] > 1.*^-6,
                    failQ = True; Break[]
                ];
                evals[[pos ;; pos + k - 1]] = bVals;
                x[[blk, pos ;; pos + k - 1]] = bX;
                xInv[[pos ;; pos + k - 1, blk]] = bXInv;
                pos += k
            ]
        ],
        {blk, neededBlocks}
    ];
    If[failQ, Return[$Failed]];
    <|
        "Evals" -> toPackedComplexVec[evals],
        "X" -> toPackedComplexMat[x],
        "XInv" -> toPackedComplexMat[xInv],
        "Dim" -> dim,
        "SubspaceDim" -> subDim,
        "BlockCount" -> Length[blocks],
        "ActiveBlockCount" -> Length[neededBlocks]
    |>
];

DiagonalizeLiouvilleGeneral[lMat_?MatrixQ] := DiagonalizeLiouvilleBlocked[lMat, All];

BuildHilbertCoherenceMask[filters_List, system_List, basis_] := Module[
    {zb = ZeemanBasis[system], dim, mLists, ordLists, mask, tMat, tMatH},
    ordLists = QM`Liouville`Private`canonicalizeCoherenceOrders /@ (Last /@ filters);
    If[MemberQ[ordLists, $Failed], Return[$Failed]];
    dim = Length[BasisKets[zb]];
    mLists = Table[
        Re[Diagonal[Normal[MatrixRepresentation[opI[f[[1]], "z"], zb]]]],
        {f, filters}
    ];
    mask = toPackedRealMat[
        Table[
            Boole[And @@ Table[
                With[{qVal = Round[mLists[[k, r]] - mLists[[k, c]]]},
                    MemberQ[ordLists[[k]], qVal]
                ],
                {k, Length[filters]}
            ]],
            {r, dim}, {c, dim}
        ]
    ];
    If[basis === zb,
        <|"Mask" -> mask, "Transform" -> None, "TransformH" -> None|>,
        tMat = toPackedComplexMat[Transpose[VectorRepresentation[#, zb] & /@ BasisKets[basis]]];
        tMatH = toPackedComplexMat[ConjugateTranspose[tMat]];
        <|"Mask" -> mask, "Transform" -> tMat, "TransformH" -> tMatH|>
    ]
];

ApplyHilbertMixSequence[rhoMat_?MatrixQ, mixSteps_List] := Fold[
    Function[{mat, step},
        Which[
            step[[1]] === "Unitary",
                step[[2]] . mat . step[[3]],
            step[[1]] === "Filter",
                With[{ms = step[[2]]},
                    If[ms["Transform"] === None,
                        ms["Mask"] * mat,
                        ms["TransformH"] . (ms["Mask"] * (ms["Transform"] . mat . ms["TransformH"])) . ms["Transform"]
                    ]
                ],
            True, mat
        ]
    ],
    rhoMat,
    mixSteps
];

applyRowModulators[mat_?MatrixQ, mod_] := If[
    mod === None,
    toPackedComplexMat[mat],
    toPackedComplexMat[toPackedComplexVec[mod] * mat]
];

(* Fast vectorized Hilbert 2D transition-table extractor for Standard single-unitary components *)
extractHilbertComponentPeaks[rMat_, wMat_, dAct_, actPairs_, e1_, w2Act_, weight_, dropTol_: 1.*^-14] := Module[
    {d = Length[e1], uL, uR, maxR, nzPairs, rCut, bagA, bagW1, bagW2, a, b, rab, ampVec, maxAbs, keepIdx, partsA, partsW1, partsW2},
    maxR = Max[Abs[rMat]];
    If[maxR == 0., Return[{{}, {}, {}}]];
    rCut = 1.*^-15 * maxR;
    nzPairs = Position[Abs[rMat], x_ /; x > rCut, {2}, Heads -> False];
    If[nzPairs === {}, Return[{{}, {}, {}}]];
    uL = toPackedComplexMat[wMat[[actPairs[[All, 1]], All]]];
    uR = toPackedComplexMat[Conjugate[wMat[[actPairs[[All, 2]], All]]]];
    bagA = Internal`Bag[]; bagW1 = Internal`Bag[]; bagW2 = Internal`Bag[];
    Do[
        {a, b} = pair;
        rab = rMat[[a, b]] * weight;
        ampVec = toPackedComplexVec[rab * dAct * uL[[All, a]] * uR[[All, b]]];
        maxAbs = Max[Abs[ampVec]];
        If[maxAbs > 0.,
            keepIdx = Pick[Range[Length[dAct]], UnitStep[Abs[ampVec] - dropTol * maxR], 1];
            If[keepIdx =!= {},
                Internal`StuffBag[bagA, toPackedComplexVec[ampVec[[keepIdx]]]];
                Internal`StuffBag[bagW1, toPackedRealVec[ConstantArray[-(e1[[a]] - e1[[b]]), Length[keepIdx]]]];
                Internal`StuffBag[bagW2, toPackedRealVec[w2Act[[keepIdx]]]]
            ]
        ],
        {pair, nzPairs}
    ];
    partsA = Internal`BagPart[bagA, All];
    If[partsA === {}, Return[{{}, {}, {}}]];
    partsW1 = Internal`BagPart[bagW1, All];
    partsW2 = Internal`BagPart[bagW2, All];
    {toPackedComplexVec[If[Length[partsA] === 1, partsA[[1]], Join @@ partsA]],
     toPackedRealVec[If[Length[partsW1] === 1, partsW1[[1]], Join @@ partsW1]],
     toPackedRealVec[If[Length[partsW2] === 1, partsW2[[1]], Join @@ partsW2]]}
];

(* Core Hilbert Double-Diagonalization 2D Engine *)
EvaluateHilbertDoubleDiag2D[spec_Association] := Module[
    {ev1Mode = Lookup[spec, "Ev1Mode", "Standard"],
     h2Diag = spec["H2Diag"], t1Vec = spec["Time1"], t2Vec = spec["Time2"],
     qAdjMat = toPackedComplexMat[spec["QAdjMat"]],
     normDenom = N[Re[spec["NormDenom"]]],
     components = spec["Components"],
     computePeaksQ = TrueQ[Lookup[spec, "ComputePeaks", False]],
     expMeth = Lookup[spec, "ExpansionMethod", "Legacy"],
     expTol = Lookup[spec, "ExpansionTolerance", 1.*^-9],
     compBackend = Lookup[spec, "CompilationBackend", Automatic],
     extractOnlyQ = TrueQ[Lookup[spec, "ExtractPeaksOnly", False]],
     n1, n2, v2, vh2, e2, dMat, dMatT, actPairs, nAct, dAct, w2Act, de2Mat,
     use2DNUFFTQ = False, peaksByKey = <||>, bagsByKey = <||>, modByKey = <||>,
     results = <||>, peaks2D = None},

    n1 = Length[t1Vec];
    n2 = Length[t2Vec];
    v2 = h2Diag["V"];
    vh2 = h2Diag["VH"];
    e2 = h2Diag["Evals"];

    (* Detection matrix in V2 eigenbasis: Tr[D . rho2] = Sum_{c,d} D[[d,c]] * rho2[[c,d]] *)
    dMat = vh2 . qAdjMat . v2;
    dMatT = Transpose[dMat];
    actPairs = Position[Abs[dMatT], x_ /; x > 1.*^-14, {2}, Heads -> False];
    nAct = Length[actPairs];

    If[nAct === 0,
        If[extractOnlyQ,
            Do[peaksByKey[comp["Key"]] = {{}, {}, {}}, {comp, components}];
            Return[<|"PeaksByKey" -> peaksByKey|>]
        ];
        Do[results[comp["Key"]] = toPackedComplexMat[ConstantArray[0. + 0. I, {n1, n2}]], {comp, components}];
        Return[<|"DataMap" -> results, "Peaks2D" -> {}|>]
    ];

    dAct = toPackedComplexVec[Extract[dMatT, actPairs] / normDenom];
    w2Act = toPackedRealVec[-(e2[[actPairs[[All, 1]]]] - e2[[actPairs[[All, 2]]]])];

    use2DNUFFTQ = (!extractOnlyQ) && (ev1Mode === "Standard") && (expMeth =!= "Legacy") &&
        AllTrue[components, TrueQ[Lookup[#, "IsSingleUnitary", False]] && Length[#["RhoStates"]] === 1 &] &&
        (MatchQ[expMeth, "2DNUFFT" | "NUFFT" | "2DSemiSeparable"] ||
         (MatchQ[expMeth, "Fast" | Automatic] && nAct * n1 * n2 >= 2*^6 && n1 >= 32 && n2 >= 128));

    If[extractOnlyQ || use2DNUFFTQ,
        Module[{h1Diag = spec["H1Diag"], v1, vh1, e1, nuMeth2D},
            v1 = h1Diag["V"]; vh1 = h1Diag["VH"]; e1 = h1Diag["Evals"];
            Do[bagsByKey[comp["Key"]] = {Internal`Bag[], Internal`Bag[], Internal`Bag[]}, {comp, components}];
            Do[
                Module[{rEig = vh1 . toPackedComplexMat[comp["RhoStates"][[1]]] . v1,
                        uMix = Lookup[comp, "SingleUnitary", IdentityMatrix[Length[v1]]],
                        wMat, key = comp["Key"], wt = N[Lookup[comp, "Weight", 1.]],
                        extA, extW1, extW2},
                    wMat = vh2 . toPackedComplexMat[uMix] . v1;
                    modByKey[key] = Lookup[comp, "RowModulators", None];
                    {extA, extW1, extW2} = extractHilbertComponentPeaks[rEig, wMat, dAct, actPairs, e1, w2Act, wt, 1.*^-14];
                    If[extA =!= {},
                        Internal`StuffBag[bagsByKey[key][[1]], extA];
                        Internal`StuffBag[bagsByKey[key][[2]], extW1];
                        Internal`StuffBag[bagsByKey[key][[3]], extW2]
                    ]
                ],
                {comp, components}
            ];
            Do[
                With[{pA = Internal`BagPart[bagsByKey[k][[1]], All],
                      pW1 = Internal`BagPart[bagsByKey[k][[2]], All],
                      pW2 = Internal`BagPart[bagsByKey[k][[3]], All]},
                    peaksByKey[k] = If[pA === {},
                        {{}, {}, {}},
                        {toPackedComplexVec[If[Length[pA] === 1, pA[[1]], Join @@ pA]],
                         toPackedRealVec[If[Length[pW1] === 1, pW1[[1]], Join @@ pW1]],
                         toPackedRealVec[If[Length[pW2] === 1, pW2[[1]], Join @@ pW2]]}
                    ]
                ],
                {k, Keys[bagsByKey]}
            ];
            If[extractOnlyQ, Return[<|"PeaksByKey" -> peaksByKey|>]];
            nuMeth2D = If[MatchQ[expMeth, "Fast" | Automatic], "2DNUFFT", expMeth];
            Do[
                With[{pk = peaksByKey[k]},
                    results[k] = If[pk[[1]] === {},
                        toPackedComplexMat[ConstantArray[0. + 0. I, {n1, n2}]],
                        Module[{z0 = toPackedRealVec[ConstantArray[0., Length[pk[[1]]]]], p5, mat},
                            p5 = Transpose[{pk[[1]], pk[[2]], z0, pk[[3]], z0}];
                            mat = ExpandLorentzian2DFast[p5, t1Vec, t2Vec,
                                ExpansionMethod -> nuMeth2D, ExpansionTolerance -> expTol, CompilationBackend -> compBackend];
                            applyRowModulators[mat, Lookup[modByKey, k, None]]
                        ]
                    ]
                ],
                {k, Keys[peaksByKey]}
            ];
            Return[<|"DataMap" -> results, "Peaks2D" -> None|>]
        ]
    ];

    de2Mat = toPackedComplexMat[dAct * Exp[I * Outer[Times, w2Act, t2Vec]]];

    If[ev1Mode === "Standard",
        Module[{h1Diag = spec["H1Diag"], v1, vh1, e1, p1Mat},
            v1 = h1Diag["V"];
            vh1 = h1Diag["VH"];
            e1 = h1Diag["Evals"];
            p1Mat = toPackedComplexMat[Exp[-I * Outer[Times, t1Vec, e1]]];

            Do[
                Module[{rhoList = comp["RhoStates"], rEigList, isSingleU = TrueQ[Lookup[comp, "IsSingleUnitary", False]],
                        mixSteps = Lookup[comp, "MixSteps", {}], nStates, aMat, sMat, key = comp["Key"],
                        wMat, uMix, pVec, wP, rho2},
                    rEigList = (vh1 . toPackedComplexMat[#] . v1 &) /@ rhoList;
                    nStates = Length[rEigList];
                    aMat = ConstantArray[0. + 0. I, {n1, nAct}];

                    If[isSingleU,
                        uMix = Lookup[comp, "SingleUnitary", IdentityMatrix[Length[v1]]];
                        wMat = vh2 . toPackedComplexMat[uMix] . v1;
                        Do[
                            pVec = p1Mat[[k1]];
                            wP = Transpose[pVec * Transpose[wMat]];
                            rho2 = wP . rEigList[[Mod[k1 - 1, nStates] + 1]] . ConjugateTranspose[wP];
                            aMat[[k1]] = Extract[rho2, actPairs],
                            {k1, n1}
                        ];
                        If[computePeaksQ && Length[components] === 1 && nStates === 1 && Length[e1] <= 16,
                            Module[{d = Length[e1], rMat = rEigList[[1]], wConj = Conjugate[wMat], pkList},
                                pkList = Flatten[Table[
                                    With[{amp = (dAct[[idx]] * wMat[[actPairs[[idx, 1]], a]] * rMat[[a, b]] * wConj[[actPairs[[idx, 2]], b]])},
                                        If[Abs[amp] > 1.*^-13,
                                            {N[amp], N[-(e1[[a]] - e1[[b]])], 0., N[w2Act[[idx]]], 0.},
                                            Nothing
                                        ]
                                    ],
                                    {a, d}, {b, d}, {idx, nAct}
                                ], 2];
                                peaks2D = pkList
                            ]
                        ],
                        Do[
                            pVec = p1Mat[[k1]];
                            wP = Transpose[pVec * Transpose[v1]];
                            rho2 = vh2 . ApplyHilbertMixSequence[
                                wP . rEigList[[Mod[k1 - 1, nStates] + 1]] . ConjugateTranspose[wP],
                                mixSteps
                            ] . v2;
                            aMat[[k1]] = Extract[rho2, actPairs],
                            {k1, n1}
                        ]
                    ];
                    sMat = applyRowModulators[toPackedComplexMat[aMat] . de2Mat, Lookup[comp, "RowModulators", None]];
                    If[KeyExistsQ[results, key],
                        results[key] = toPackedComplexMat[results[key] + Lookup[comp, "Weight", 1.] * sMat],
                        results[key] = toPackedComplexMat[Lookup[comp, "Weight", 1.] * sMat]
                    ]
                ],
                {comp, components}
            ]
        ],
        (* "Echo" (t1/2 - Refocus - t1/2) or "ConstantTime" ((Tct + t1)/2 - Refocus - (Tct - t1)/2) *)
        Module[{h1aDiag = spec["H1aDiag"], h1bDiag = spec["H1bDiag"], refSteps = Lookup[spec, "RefocusSteps", {}],
                refSingleU = Lookup[spec, "RefocusSingleUnitary", None],
                v1a, vh1a, e1a, v1b, vh1b, e1b, t1aVec, t1bVec, p1aMat, p1bMat, wRef},
            v1a = h1aDiag["V"]; vh1a = h1aDiag["VH"]; e1a = h1aDiag["Evals"];
            v1b = h1bDiag["V"]; vh1b = h1bDiag["VH"]; e1b = h1bDiag["Evals"];
            If[ev1Mode === "ConstantTime",
                With[{tCT = N[Re[spec["CTDuration"]]]},
                    t1aVec = toPackedRealVec[0.5 * (tCT + t1Vec)];
                    t1bVec = toPackedRealVec[0.5 * (tCT - t1Vec)]
                ],
                t1aVec = toPackedRealVec[0.5 * t1Vec];
                t1bVec = toPackedRealVec[0.5 * t1Vec]
            ];
            p1aMat = toPackedComplexMat[Exp[-I * Outer[Times, t1aVec, e1a]]];
            p1bMat = toPackedComplexMat[Exp[-I * Outer[Times, t1bVec, e1b]]];
            wRef = If[refSingleU =!= None, vh1b . toPackedComplexMat[refSingleU] . v1a, None];

            Do[
                Module[{rhoList = comp["RhoStates"], rEigList, isSingleU = TrueQ[Lookup[comp, "IsSingleUnitary", False]],
                        mixSteps = Lookup[comp, "MixSteps", {}], nStates, aMat, sMat, key = comp["Key"],
                        wMix, uMix, p1a, p1b, rhoHalf1, rhoAfterRef, rho2, wA, wB},
                    rEigList = (vh1a . toPackedComplexMat[#] . v1a &) /@ rhoList;
                    nStates = Length[rEigList];
                    aMat = ConstantArray[0. + 0. I, {n1, nAct}];
                    wMix = If[isSingleU,
                        vh2 . toPackedComplexMat[Lookup[comp, "SingleUnitary", IdentityMatrix[Length[v2]]]] . v1b,
                        None
                    ];
                    Do[
                        p1a = p1aMat[[k1]];
                        p1b = p1bMat[[k1]];
                        rhoAfterRef = If[wRef =!= None,
                            wA = Transpose[p1a * Transpose[wRef]];
                            wA . rEigList[[Mod[k1 - 1, nStates] + 1]] . ConjugateTranspose[wA],
                            wA = Transpose[p1a * Transpose[v1a]];
                            vh1b . ApplyHilbertMixSequence[
                                wA . rEigList[[Mod[k1 - 1, nStates] + 1]] . ConjugateTranspose[wA],
                                refSteps
                            ] . v1b
                        ];
                        rho2 = If[wMix =!= None,
                            wB = Transpose[p1b * Transpose[wMix]];
                            wB . rhoAfterRef . ConjugateTranspose[wB],
                            wB = Transpose[p1b * Transpose[v1b]];
                            vh2 . ApplyHilbertMixSequence[
                                wB . rhoAfterRef . ConjugateTranspose[wB],
                                mixSteps
                            ] . v2
                        ];
                        aMat[[k1]] = Extract[rho2, actPairs],
                        {k1, n1}
                    ];
                    sMat = applyRowModulators[toPackedComplexMat[aMat] . de2Mat, Lookup[comp, "RowModulators", None]];
                    If[KeyExistsQ[results, key],
                        results[key] = toPackedComplexMat[results[key] + Lookup[comp, "Weight", 1.] * sMat],
                        results[key] = toPackedComplexMat[Lookup[comp, "Weight", 1.] * sMat]
                    ]
                ],
                {comp, components}
            ]
        ]
    ];

    <|"DataMap" -> results, "Peaks2D" -> peaks2D|>
];

(* Core Liouville Double-Diagonalization BLAS-3 Engine *)
EvaluateLiouvilleDoubleDiag2D[spec_Association] := Module[
    {ev1Mode = Lookup[spec, "Ev1Mode", "Standard"],
     l2Diag = spec["L2Diag"], t1Vec = spec["Time1"], t2Vec = spec["Time2"],
     qVec = toPackedComplexVec[spec["QVec"]],
     normDenom = N[Re[spec["NormDenom"]]],
     components = spec["Components"],
     n1, n2, x2, x2Inv, e2, d2, act2, nAct2, de2Mat, results = <||>},

    n1 = Length[t1Vec];
    n2 = Length[t2Vec];
    x2 = l2Diag["X"];
    x2Inv = l2Diag["XInv"];
    e2 = l2Diag["Evals"];

    d2 = If[Length[e2] === 0, toPackedComplexVec[{}], toPackedComplexVec[(Conjugate[qVec] . x2) / normDenom]];
    act2 = If[Length[d2] === 0, {}, Pick[Range[Length[d2]], UnitStep[Abs[d2] - 1.*^-14], 1]];
    nAct2 = Length[act2];

    If[nAct2 === 0,
        Do[results[comp["Key"]] = toPackedComplexMat[ConstantArray[0. + 0. I, {n1, n2}]], {comp, components}];
        Return[<|"DataMap" -> results, "Peaks2D" -> None|>]
    ];

    de2Mat = toPackedComplexMat[d2[[act2]] * Exp[Outer[Times, e2[[act2]], t2Vec]]];

    If[ev1Mode === "Standard",
        Module[{l1Diag = spec["L1Diag"], x1, x1Inv, e1},
            x1 = l1Diag["X"];
            x1Inv = l1Diag["XInv"];
            e1 = l1Diag["Evals"];

            Do[
                Module[{rhoList = comp["RhoVecs"], mixMat = toPackedComplexMat[comp["MixMat"]],
                        key = comp["Key"], nStates = Length[comp["RhoVecs"]], sMat,
                        r1Full, act1, mEigFull, mEig, cMat, e1Mat, subIdx, sSub},
                    If[Length[e1] === 0,
                        sMat = toPackedComplexMat[ConstantArray[0. + 0. I, {n1, n2}]],
                    If[nStates === 1,
                        r1Full = x1Inv . toPackedComplexVec[rhoList[[1]]];
                        act1 = Pick[Range[Length[r1Full]], UnitStep[Abs[r1Full] - 1.*^-14], 1];
                        If[act1 === {},
                            sMat = toPackedComplexMat[ConstantArray[0. + 0. I, {n1, n2}]],
                            mEig = x2Inv[[act2]] . mixMat . x1[[All, act1]];
                            cMat = Transpose[r1Full[[act1]] * Transpose[mEig]];
                            e1Mat = toPackedComplexMat[Exp[Outer[Times, t1Vec, e1[[act1]]]]];
                            sMat = toPackedComplexMat[(e1Mat . Transpose[cMat]) . de2Mat]
                        ],
                        (* Periodic phase cycling in t1 (StatesTPPI nStates=2 or TPPI nStates=4): evaluate each subgrid in BLAS-3 *)
                        mEigFull = x2Inv[[act2]] . mixMat . x1;
                        sMat = ConstantArray[0. + 0. I, {n1, n2}];
                        Do[
                            subIdx = Range[stateIdx, n1, nStates];
                            If[subIdx =!= {},
                                r1Full = x1Inv . toPackedComplexVec[rhoList[[stateIdx]]];
                                act1 = Pick[Range[Length[r1Full]], UnitStep[Abs[r1Full] - 1.*^-14], 1];
                                If[act1 =!= {},
                                    cMat = Transpose[r1Full[[act1]] * Transpose[mEigFull[[All, act1]]]];
                                    e1Mat = toPackedComplexMat[Exp[Outer[Times, t1Vec[[subIdx]], e1[[act1]]]]];
                                    sSub = toPackedComplexMat[(e1Mat . Transpose[cMat]) . de2Mat];
                                    sMat[[subIdx]] = sSub
                                ]
                            ],
                            {stateIdx, nStates}
                        ];
                        sMat = toPackedComplexMat[sMat]
                    ]];
                    sMat = applyRowModulators[sMat, Lookup[comp, "RowModulators", None]];
                    If[KeyExistsQ[results, key],
                        results[key] = toPackedComplexMat[results[key] + Lookup[comp, "Weight", 1.] * sMat],
                        results[key] = toPackedComplexMat[Lookup[comp, "Weight", 1.] * sMat]
                    ]
                ],
                {comp, components}
            ]
        ],
        (* Liouville "Echo" or "ConstantTime" *)
        Module[{l1aDiag = spec["L1aDiag"], l1bDiag = spec["L1bDiag"], refMat = toPackedComplexMat[spec["RefocusMat"]],
                x1a, x1aInv, e1a, x1b, x1bInv, e1b, t1aVec, t1bVec, p1aMat, p1bMat, mRef},
            x1a = l1aDiag["X"]; x1aInv = l1aDiag["XInv"]; e1a = l1aDiag["Evals"];
            x1b = l1bDiag["X"]; x1bInv = l1bDiag["XInv"]; e1b = l1bDiag["Evals"];
            If[ev1Mode === "ConstantTime",
                With[{tCT = N[Re[spec["CTDuration"]]]},
                    t1aVec = toPackedRealVec[0.5 * (tCT + t1Vec)];
                    t1bVec = toPackedRealVec[0.5 * (tCT - t1Vec)]
                ],
                t1aVec = toPackedRealVec[0.5 * t1Vec];
                t1bVec = toPackedRealVec[0.5 * t1Vec]
            ];
            p1aMat = toPackedComplexMat[Exp[Outer[Times, t1aVec, e1a]]];
            p1bMat = toPackedComplexMat[Exp[Outer[Times, t1bVec, e1b]]];
            mRef = x1bInv . refMat . x1a;

            Do[
                Module[{rhoList = comp["RhoVecs"], rEigList, nStates, mixMat = toPackedComplexMat[comp["MixMat"]],
                        mMix, aMat, sMat, key = comp["Key"], vMid},
                    rEigList = (x1aInv . toPackedComplexVec[#] &) /@ rhoList;
                    nStates = Length[rEigList];
                    mMix = x2Inv[[act2]] . mixMat . x1b;
                    aMat = ConstantArray[0. + 0. I, {n1, nAct2}];
                    Do[
                        vMid = mRef . (p1aMat[[k1]] * rEigList[[Mod[k1 - 1, nStates] + 1]]);
                        aMat[[k1]] = mMix . (p1bMat[[k1]] * vMid),
                        {k1, n1}
                    ];
                    sMat = applyRowModulators[toPackedComplexMat[aMat] . de2Mat, Lookup[comp, "RowModulators", None]];
                    If[KeyExistsQ[results, key],
                        results[key] = toPackedComplexMat[results[key] + Lookup[comp, "Weight", 1.] * sMat],
                        results[key] = toPackedComplexMat[Lookup[comp, "Weight", 1.] * sMat]
                    ]
                ],
                {comp, components}
            ]
        ]
    ];

    <|"DataMap" -> results, "Peaks2D" -> None|>
];

(* Core Uniform Stepped-Propagator Marching 2D Engine (with optional diagonalized F2 fast finish) *)
EvaluateSteppedPropagator2D[spec_Association] := Module[
    {space = spec["Space"], ev1Mode = Lookup[spec, "Ev1Mode", "Standard"],
     n1 = spec["N1"], n2 = spec["N2"],
     normDenom = N[Re[spec["NormDenom"]]],
     components = spec["Components"],
     initProp1a = Lookup[spec, "InitProp1a", None],
     stepProp1a = Lookup[spec, "StepProp1a", None],
     initProp1b = Lookup[spec, "InitProp1b", None],
     stepProp1b = Lookup[spec, "StepProp1b", None],
     refSteps = Lookup[spec, "RefocusSteps", {}],
     refMatL = Lookup[spec, "RefocusMat", None],
     diag2 = Lookup[spec, "Diag2", None],
     initProp2 = Lookup[spec, "InitProp2", None],
     stepProp2 = Lookup[spec, "StepProp2", None],
     t2Vec = Lookup[spec, "Time2", {}],
     qObj = spec["QObj"],
     results = <||>, useDiag2Q,
     v2, vh2, e2, dMatT, actPairs, nAct, dAct, w2Act, de2Mat,
     x2, x2Inv, d2, act2, nAct2, qAdjFlat, qConjVec},

    useDiag2Q = AssociationQ[diag2];
    If[useDiag2Q,
        If[space === "Hilbert",
            v2 = diag2["V"]; vh2 = diag2["VH"]; e2 = diag2["Evals"];
            dMatT = Transpose[vh2 . toPackedComplexMat[qObj] . v2];
            actPairs = Position[Abs[dMatT], x_ /; x > 1.*^-14, {2}, Heads -> False];
            nAct = Length[actPairs];
            If[nAct === 0,
                Do[results[c["Key"]] = toPackedComplexMat[ConstantArray[0. + 0. I, {n1, n2}]], {c, components}];
                Return[<|"DataMap" -> results, "Peaks2D" -> None|>]
            ];
            dAct = toPackedComplexVec[Extract[dMatT, actPairs] / normDenom];
            w2Act = toPackedRealVec[-(e2[[actPairs[[All, 1]]]] - e2[[actPairs[[All, 2]]]])];
            de2Mat = toPackedComplexMat[dAct * Exp[I * Outer[Times, w2Act, t2Vec]]],
            x2 = diag2["X"]; x2Inv = diag2["XInv"]; e2 = diag2["Evals"];
            d2 = toPackedComplexVec[(Conjugate[toPackedComplexVec[qObj]] . x2) / normDenom];
            act2 = Pick[Range[Length[d2]], UnitStep[Abs[d2] - 1.*^-14], 1];
            nAct2 = Length[act2];
            If[nAct2 === 0,
                Do[results[c["Key"]] = toPackedComplexMat[ConstantArray[0. + 0. I, {n1, n2}]], {c, components}];
                Return[<|"DataMap" -> results, "Peaks2D" -> None|>]
            ];
            de2Mat = toPackedComplexMat[d2[[act2]] * Exp[Outer[Times, e2[[act2]], t2Vec]]]
        ],
        If[space === "Hilbert",
            qAdjFlat = toPackedComplexVec[Flatten[toPackedComplexMat[qObj]] / normDenom],
            qConjVec = toPackedComplexVec[Conjugate[toPackedComplexVec[qObj]] / normDenom]
        ]
    ];

    Do[
        Module[{key = comp["Key"], sMat, aMat, rhoStates, nStates,
                initStatesA, curPropA, curPropB, stepA, stepAH, stepB, stepBH, step2H,
                rhoAfter1, rhoMix, rho2, v1, v2cur, stepProp1aH, stepProp1bH,
                uRefSingle, uRefSingleH, uMixSingle, uMixSingleH, mixSteps, mixMat},
            sMat = ConstantArray[0. + 0. I, {n1, n2}];

            If[space === "Hilbert",
                rhoStates = toPackedComplexMat /@ comp["RhoStates"];
                nStates = Length[rhoStates];
                mixSteps = Lookup[comp, "MixSteps", {}];
                initStatesA = If[initProp1a =!= None,
                    With[{u0 = toPackedComplexMat[initProp1a], u0H = ConjugateTranspose[toPackedComplexMat[initProp1a]]},
                        (u0 . # . u0H &) /@ rhoStates
                    ],
                    rhoStates
                ];
                stepA = toPackedComplexMat[stepProp1a];
                stepAH = ConjugateTranspose[stepA];
                If[ev1Mode =!= "Standard",
                    stepB = toPackedComplexMat[stepProp1b];
                    stepBH = ConjugateTranspose[stepB];
                    uRefSingle = Lookup[spec, "RefocusSingleUnitary", None];
                    If[uRefSingle =!= None,
                        uRefSingle = toPackedComplexMat[uRefSingle];
                        uRefSingleH = ConjugateTranspose[uRefSingle]
                    ]
                ];
                If[TrueQ[Lookup[comp, "IsSingleUnitary", False]],
                    uMixSingle = toPackedComplexMat[comp["SingleUnitary"]];
                    uMixSingleH = ConjugateTranspose[uMixSingle]
                ];
                If[!useDiag2Q, step2H = ConjugateTranspose[toPackedComplexMat[stepProp2]]];
                If[useDiag2Q, aMat = ConstantArray[0. + 0. I, {n1, nAct}]];

                (* Propagate accumulated k1 operator U_1(t1) rather than state when nStates > 1, or march rho when nStates == 1 *)
                If[ev1Mode === "Standard",
                    If[nStates === 1,
                        rhoAfter1 = initStatesA[[1]];
                        Do[
                            rhoMix = If[TrueQ[Lookup[comp, "IsSingleUnitary", False]],
                                uMixSingle . rhoAfter1 . uMixSingleH,
                                ApplyHilbertMixSequence[rhoAfter1, mixSteps]
                            ];
                            If[useDiag2Q,
                                aMat[[k1]] = Extract[vh2 . rhoMix . v2, actPairs],
                                rho2 = If[initProp2 =!= None, initProp2 . rhoMix . ConjugateTranspose[initProp2], rhoMix];
                                Do[
                                    sMat[[k1, k2]] = qAdjFlat . Flatten[Transpose[rho2]];
                                    If[k2 < n2, rho2 = stepProp2 . rho2 . step2H],
                                    {k2, n2}
                                ]
                            ];
                            If[k1 < n1, rhoAfter1 = stepA . rhoAfter1 . stepAH],
                            {k1, n1}
                        ],
                        (* nStates > 1 (StatesTPPI / TPPI): march cumulative propagator curPropA *)
                        curPropA = If[initProp1a =!= None, toPackedComplexMat[initProp1a], toPackedComplexMat[IdentityMatrix[Length[stepA]]]];
                        Do[
                            With[{cAH = ConjugateTranspose[curPropA]},
                                rhoAfter1 = curPropA . rhoStates[[Mod[k1 - 1, nStates] + 1]] . cAH
                            ];
                            rhoMix = If[TrueQ[Lookup[comp, "IsSingleUnitary", False]],
                                uMixSingle . rhoAfter1 . uMixSingleH,
                                ApplyHilbertMixSequence[rhoAfter1, mixSteps]
                            ];
                            If[useDiag2Q,
                                aMat[[k1]] = Extract[vh2 . rhoMix . v2, actPairs],
                                rho2 = If[initProp2 =!= None, initProp2 . rhoMix . ConjugateTranspose[initProp2], rhoMix];
                                Do[
                                    sMat[[k1, k2]] = qAdjFlat . Flatten[Transpose[rho2]];
                                    If[k2 < n2, rho2 = stepProp2 . rho2 . step2H],
                                    {k2, n2}
                                ]
                            ];
                            If[k1 < n1, curPropA = stepA . curPropA],
                            {k1, n1}
                        ]
                    ],
                    (* "Echo" or "ConstantTime": march curPropA and curPropB *)
                    curPropA = If[initProp1a =!= None, toPackedComplexMat[initProp1a], toPackedComplexMat[IdentityMatrix[Length[stepA]]]];
                    curPropB = If[initProp1b =!= None, toPackedComplexMat[initProp1b], toPackedComplexMat[IdentityMatrix[Length[stepB]]]];
                    Do[
                        rhoAfter1 = curPropA . rhoStates[[Mod[k1 - 1, nStates] + 1]] . ConjugateTranspose[curPropA];
                        rhoAfter1 = If[uRefSingle =!= None,
                            uRefSingle . rhoAfter1 . uRefSingleH,
                            ApplyHilbertMixSequence[rhoAfter1, refSteps]
                        ];
                        rhoAfter1 = curPropB . rhoAfter1 . ConjugateTranspose[curPropB];
                        rhoMix = If[TrueQ[Lookup[comp, "IsSingleUnitary", False]],
                            uMixSingle . rhoAfter1 . uMixSingleH,
                            ApplyHilbertMixSequence[rhoAfter1, mixSteps]
                        ];
                        If[useDiag2Q,
                            aMat[[k1]] = Extract[vh2 . rhoMix . v2, actPairs],
                            rho2 = If[initProp2 =!= None, initProp2 . rhoMix . ConjugateTranspose[initProp2], rhoMix];
                            Do[
                                sMat[[k1, k2]] = qAdjFlat . Flatten[Transpose[rho2]];
                                If[k2 < n2, rho2 = stepProp2 . rho2 . step2H],
                                {k2, n2}
                            ]
                        ];
                        If[k1 < n1,
                            curPropA = stepA . curPropA;
                            curPropB = stepB . curPropB
                        ],
                        {k1, n1}
                    ]
                ];
                If[useDiag2Q, sMat = toPackedComplexMat[aMat] . de2Mat],

                (* Liouville SteppedPropagator *)
                rhoStates = toPackedComplexVec /@ comp["RhoVecs"];
                nStates = Length[rhoStates];
                mixMat = toPackedComplexMat[comp["MixMat"]];
                stepA = toPackedComplexMat[stepProp1a];
                If[useDiag2Q, aMat = ConstantArray[0. + 0. I, {n1, nAct2}]];

                If[ev1Mode === "Standard",
                    If[nStates === 1,
                        v1 = If[initProp1a =!= None, initProp1a . rhoStates[[1]], rhoStates[[1]]];
                        Do[
                            v2cur = mixMat . v1;
                            If[useDiag2Q,
                                aMat[[k1]] = x2Inv[[act2]] . v2cur,
                                v2cur = If[initProp2 =!= None, initProp2 . v2cur, v2cur];
                                Do[
                                    sMat[[k1, k2]] = qConjVec . v2cur;
                                    If[k2 < n2, v2cur = stepProp2 . v2cur],
                                    {k2, n2}
                                ]
                            ];
                            If[k1 < n1, v1 = stepA . v1],
                            {k1, n1}
                        ],
                        curPropA = If[initProp1a =!= None, toPackedComplexMat[initProp1a], toPackedComplexMat[IdentityMatrix[Length[stepA]]]];
                        Do[
                            v1 = curPropA . rhoStates[[Mod[k1 - 1, nStates] + 1]];
                            v2cur = mixMat . v1;
                            If[useDiag2Q,
                                aMat[[k1]] = x2Inv[[act2]] . v2cur,
                                v2cur = If[initProp2 =!= None, initProp2 . v2cur, v2cur];
                                Do[
                                    sMat[[k1, k2]] = qConjVec . v2cur;
                                    If[k2 < n2, v2cur = stepProp2 . v2cur],
                                    {k2, n2}
                                ]
                            ];
                            If[k1 < n1, curPropA = stepA . curPropA],
                            {k1, n1}
                        ]
                    ],
                    stepB = toPackedComplexMat[stepProp1b];
                    curPropA = If[initProp1a =!= None, toPackedComplexMat[initProp1a], toPackedComplexMat[IdentityMatrix[Length[stepA]]]];
                    curPropB = If[initProp1b =!= None, toPackedComplexMat[initProp1b], toPackedComplexMat[IdentityMatrix[Length[stepB]]]];
                    Do[
                        v1 = curPropB . (refMatL . (curPropA . rhoStates[[Mod[k1 - 1, nStates] + 1]]));
                        v2cur = mixMat . v1;
                        If[useDiag2Q,
                            aMat[[k1]] = x2Inv[[act2]] . v2cur,
                            v2cur = If[initProp2 =!= None, initProp2 . v2cur, v2cur];
                            Do[
                                sMat[[k1, k2]] = qConjVec . v2cur;
                                If[k2 < n2, v2cur = stepProp2 . v2cur],
                                {k2, n2}
                            ]
                        ];
                        If[k1 < n1,
                            curPropA = stepA . curPropA;
                            curPropB = stepB . curPropB
                        ],
                        {k1, n1}
                    ]
                ];
                If[useDiag2Q, sMat = toPackedComplexMat[aMat] . de2Mat]
            ];

            sMat = applyRowModulators[toPackedComplexMat[sMat], Lookup[comp, "RowModulators", None]];
            If[KeyExistsQ[results, key],
                results[key] = toPackedComplexMat[results[key] + Lookup[comp, "Weight", 1.] * sMat],
                results[key] = toPackedComplexMat[Lookup[comp, "Weight", 1.] * sMat]
            ]
        ],
        {comp, components}
    ];

    <|"DataMap" -> results, "Peaks2D" -> None|>
];

(* ====================================================================== *)
(* Phase 3 Kernels: CF4 Magnus, Periodic COMPUTE2D, and Compiled Ensemble *)
(* Pure System/Developer packed-array code: zero SpinDynamica references.  *)
(* ====================================================================== *)

CF4PropagatorStep[evalGenFn_, tStart_?NumericQ, tEnd_?NumericQ, nSub_Integer, space_String, dim_Integer] := Module[
    {dt = N[Re[tEnd - tStart]], h, pref, c1 = 0.5 - Sqrt[3.] / 6., c2 = 0.5 + Sqrt[3.] / 6.,
     a1 = 0.25 - Sqrt[3.] / 6., a2 = 0.25 + Sqrt[3.] / 6., uAcc, s, ta, g1, g2, uStep, nStep = Max[1, nSub]},
    If[Abs[dt] < 1.*^-15, Return[toPackedComplexMat[IdentityMatrix[dim]]]];
    h = dt / nStep;
    pref = If[space === "Hilbert", -I * h, h];
    uAcc = toPackedComplexMat[IdentityMatrix[dim]];
    Do[
        ta = N[Re[tStart]] + s * h;
        g1 = toPackedComplexMat[evalGenFn[ta + c1 * h]];
        g2 = toPackedComplexMat[evalGenFn[ta + c2 * h]];
        uStep = MatrixExp[pref * (a1 * g1 + a2 * g2)] . MatrixExp[pref * (a2 * g1 + a1 * g2)];
        uAcc = toPackedComplexMat[uStep . uAcc],
        {s, 0, nStep - 1}
    ];
    uAcc
];

BuildPeriodicPropagatorsCF4[evalGenFn_, t0_?NumericQ, period_?NumericQ, nDiv_Integer, nSubPerDiv_Integer, space_String, dim_Integer] := Module[
    {dTau = N[Re[period]] / nDiv, tBase = N[Re[t0]], sSet},
    sSet = Table[
        CF4PropagatorStep[evalGenFn, tBase + (j - 1) * dTau, tBase + j * dTau, nSubPerDiv, space, dim],
        {j, nDiv}
    ];
    FoldList[toPackedComplexMat[#2 . #1] &, toPackedComplexMat[IdentityMatrix[dim]], sSet]
];

diagonalizeUnitaryFloquet[uMat_?MatrixQ] := Module[
    {u = toPackedComplexMat[uMat], dim = Length[uMat], vals, vecs, v, vH, groups, rem, grp,
     gVecs, orthoG},
    {vals, vecs} = Eigensystem[u];
    v = Transpose[vecs];
    (* Orthonormalize columns cleanly across degenerate unitary eigenspaces *)
    vH = ConjugateTranspose[v];
    If[Max[Abs[vH . v - IdentityMatrix[dim]]] > 1.*^-10,
        rem = Range[dim];
        groups = {};
        While[rem =!= {},
            With[{anchor = vals[[First[rem]]]},
                grp = Select[rem, Abs[vals[[#]] - anchor] < 1.*^-8 &];
                rem = Complement[rem, grp];
                AppendTo[groups, grp]
            ]
        ];
        Do[
            gVecs = Transpose[v[[All, g]]];
            orthoG = Orthogonalize[gVecs, Conjugate[#1] . #2 &];
            v[[All, g]] = Transpose[orthoG],
            {g, groups}
        ];
        vH = ConjugateTranspose[v];
        If[Max[Abs[vH . v - IdentityMatrix[dim]]] > 1.*^-10,
            v = Transpose[Orthogonalize[Transpose[v], Conjugate[#1] . #2 &]];
            vH = ConjugateTranspose[v]
        ]
    ];
    <|"Lambda" -> toPackedComplexVec[vals], "V" -> toPackedComplexMat[v], "VH" -> toPackedComplexMat[vH], "Dim" -> dim|>
];

splitFloquetPhaseTime[tRel_?NumericQ, period_?NumericQ, dTau_?NumericQ, nDiv_Integer] := Module[
    {m, rem, jNear, deltaNear, jOut, deltaOut},
    m = Floor[(tRel + 1.*^-12 * period) / period];
    rem = tRel - m * period;
    If[rem < 0., rem = 0.];
    jNear = Round[rem / dTau];
    deltaNear = rem - jNear * dTau;
    If[Abs[deltaNear] <= 1.*^-11 * period,
        If[jNear >= nDiv, {m + 1, 0, 0.}, {m, jNear, 0.}],
        jOut = Floor[rem / dTau];
        deltaOut = rem - jOut * dTau;
        If[jOut >= nDiv, {m + 1, 0, 0.}, {m, jOut, deltaOut}]
    ]
];

applyHilbertTransformMask[mat_?MatrixQ, ms_Association] := If[
    ms["Transform"] === None,
    ms["Mask"] * mat,
    ms["TransformH"] . (ms["Mask"] * (ms["Transform"] . mat . ms["TransformH"])) . ms["Transform"]
];

EvaluatePeriodicCOMPUTE2D[spec_Association] := Module[
    {space = spec["Space"],
     period = N[Re[spec["Period"]]], nDiv = spec["NDiv"], t0 = N[Re[Lookup[spec, "T0", 0.]]],
     dTau,
     aSet1 = toPackedComplexMat /@ spec["ASet1"],
     aSet2 = toPackedComplexMat /@ spec["ASet2"],
     sameQ = TrueQ[Lookup[spec, "SameEv1Ev2", False]],
     evalGen1 = Lookup[spec, "EvalGen1", None],
     evalGen2 = Lookup[spec, "EvalGen2", None],
     microSub = Lookup[spec, "MicroSubsteps", 2],
     tauP = N[Re[Lookup[spec, "PrepDuration", 0.]]],
     tauM = N[Re[Lookup[spec, "MixDuration", 0.]]],
     gammaOffsets = Lookup[spec, "GammaOffsets", {0}],
     t1Vec = toPackedRealVec[spec["Time1"]],
     t2Vec = toPackedRealVec[spec["Time2"]],
     qObj = spec["QObj"],
     normDenom = N[Re[spec["NormDenom"]]],
     components = spec["Components"],
     n1, n2, nG, dim, results = <||>,
     fDiag1, fDiag2, lDiag1 = None, lDiag2 = None, diagMode2Q = True,
     getA1At, getA2At, getA1InvAt, getA2InvAt, uMicroCache1 = <||>, uMicroCache2 = <||>,
     a1MatCache = <||>, a2MatCache = <||>, de2Cache = <||>, powCache1 = <||>, powCache2 = <||>,
     phaseKeyOf, lookupMixStep, buildDe2Hilbert, buildDe2Liouville},

    n1 = Length[t1Vec];
    n2 = Length[t2Vec];
    nG = Length[gammaOffsets];
    dTau = period / nDiv;
    dim = Length[First[aSet1]];
    phaseKeyOf[phi_?NumericQ] := Round[Mod[phi, period], 1.*^-11];

    getA1At[j_Integer, delta_?NumericQ] := If[
        Abs[delta] <= 1.*^-11 * period,
        aSet1[[j + 1]],
        With[{k = {j, Round[delta, 1.*^-13]}},
            Lookup[uMicroCache1, Key[k],
                uMicroCache1[k] = toPackedComplexMat[
                    CF4PropagatorStep[evalGen1, t0 + j * dTau, t0 + j * dTau + delta, microSub, space, dim] . aSet1[[j + 1]]
                ]
            ]
        ]
    ];
    getA2At[j_Integer, delta_?NumericQ] := If[sameQ,
        getA1At[j, delta],
        If[Abs[delta] <= 1.*^-11 * period,
            aSet2[[j + 1]],
            With[{k = {j, Round[delta, 1.*^-13]}},
                Lookup[uMicroCache2, Key[k],
                    uMicroCache2[k] = toPackedComplexMat[
                        CF4PropagatorStep[evalGen2, t0 + j * dTau, t0 + j * dTau + delta, microSub, space, dim] . aSet2[[j + 1]]
                    ]
                ]
            ]
        ]
    ];

    If[space === "Hilbert",
        fDiag1 = diagonalizeUnitaryFloquet[Last[aSet1]];
        fDiag2 = If[sameQ, fDiag1, diagonalizeUnitaryFloquet[Last[aSet2]]];
        getA1InvAt[j_Integer, delta_?NumericQ] := ConjugateTranspose[getA1At[j, delta]];
        getA2InvAt[j_Integer, delta_?NumericQ] := ConjugateTranspose[getA2At[j, delta]],
        lDiag1 = DiagonalizeLiouvilleGeneral[Last[aSet1]];
        lDiag2 = If[sameQ, lDiag1, DiagonalizeLiouvilleGeneral[Last[aSet2]]];
        diagMode2Q = AssociationQ[lDiag2];
        a1MatCache["Inv"] = Table[Inverse[aSet1[[j]]], {j, nDiv}];
        a2MatCache["Inv"] = If[sameQ, a1MatCache["Inv"], Table[Inverse[aSet2[[j]]], {j, nDiv}]];
        getA1InvAt[j_Integer, delta_?NumericQ] := If[
            Abs[delta] <= 1.*^-11 * period,
            a1MatCache["Inv"][[j + 1]],
            Inverse[getA1At[j, delta]]
        ];
        getA2InvAt[j_Integer, delta_?NumericQ] := If[sameQ,
            getA1InvAt[j, delta],
            If[Abs[delta] <= 1.*^-11 * period,
                a2MatCache["Inv"][[j + 1]],
                Inverse[getA2At[j, delta]]
            ]
        ]
    ];

    lookupMixStep[mixByPhase_, endRel_?NumericQ] := If[
        AssociationQ[mixByPhase],
        If[KeyExistsQ[mixByPhase, "Static"],
            mixByPhase["Static"],
            With[{pk = phaseKeyOf[endRel], kBin = Mod[Round[Mod[endRel, period] / dTau], nDiv]},
                Which[
                    KeyExistsQ[mixByPhase, pk], mixByPhase[pk],
                    KeyExistsQ[mixByPhase, kBin], mixByPhase[kBin],
                    True, First[Values[mixByPhase]]
                ]
            ]
        ],
        If[Length[mixByPhase] === 1,
            mixByPhase[[1]],
            mixByPhase[[Mod[Round[Mod[endRel, period] / dTau], nDiv] + 1]]
        ]
    ];

    If[space === "Hilbert",
        Module[{v1 = fDiag1["V"], vh1 = fDiag1["VH"], lam1 = fDiag1["Lambda"], lam1C,
                v2 = fDiag2["V"], vh2 = fDiag2["VH"], lam2 = fDiag2["Lambda"], lam2C,
                qAdjMat = toPackedComplexMat[qObj], dMatCache = <||>, d2PairMat},
            lam1C = Conjugate[lam1];
            lam2C = Conjugate[lam2];
            d2PairMat = Outer[Times, lam2, lam2C];

            buildDe2Hilbert[phiM_?NumericQ] := With[{pk = phaseKeyOf[phiM]},
                Lookup[de2Cache, Key[pk],
                    de2Cache[pk] = Module[{mM, jM, deltaM, aMInv, fMInv, fM, dCols, dFlat, k2,
                                           m2, j2, delta2, a2, dm, powMat, dJ, actIdx},
                        {mM, jM, deltaM} = splitFloquetPhaseTime[phiM, period, dTau, nDiv];
                        aMInv = getA2InvAt[jM, deltaM];
                        fMInv = vh2 . aMInv;
                        fM = ConjugateTranspose[fMInv];
                        dCols = ConstantArray[0. + 0. I, {dim * dim, n2}];
                        Do[
                            {m2, j2, delta2} = splitFloquetPhaseTime[phiM + t2Vec[[k2]], period, dTau, nDiv];
                            dm = m2 - mM;
                            If[Abs[delta2] <= 1.*^-11 * period,
                                dJ = Lookup[dMatCache, Key[j2],
                                    dMatCache[j2] = With[{fj = aSet2[[j2 + 1]] . v2},
                                        Transpose[ConjugateTranspose[fj] . qAdjMat . fj] / normDenom
                                    ]
                                ],
                                a2 = getA2At[j2, delta2] . v2;
                                dJ = Transpose[ConjugateTranspose[a2] . qAdjMat . a2] / normDenom
                            ];
                            powMat = If[dm === 0, 1., Power[d2PairMat, dm]];
                            dCols[[All, k2]] = Flatten[dJ * powMat],
                            {k2, n2}
                        ];
                        actIdx = Pick[Range[dim * dim], UnitStep[(Max /@ Abs[dCols]) - 1.*^-14], 1];
                        <|"FStartInv" -> fMInv, "FStart" -> fM,
                          "ActiveIdx" -> actIdx,
                          "De2" -> If[actIdx === {}, None, toPackedComplexMat[dCols[[actIdx]]]]|>
                    ]
                ]
            ];

            Do[
                Module[{key = comp["Key"], weight = N[Lookup[comp, "Weight", 1.]],
                        mod = Lookup[comp, "RowModulators", None],
                        compTauP = N[Re[Lookup[comp, "PrepDuration", tauP]]],
                        compTauM = N[Re[Lookup[comp, "MixDuration", tauM]]],
                        rhoByG = comp["RhoByGamma"], mixByPhase = comp["MixByPhase"],
                        u1ByK, aMat, sMat, gIdx, tGStart, t1Start, t1End, rhoList, nStates,
                        m1s, j1s, d1s, m1e, j1e, d1e, dm1, u1Prop, p1, p1C, fStartInv, fEnd,
                        rho0, rho1, rhoM, mixStep, pkM, groupsM, deInf, actIdx, k1},
                    sMat = ConstantArray[0. + 0. I, {n1, n2}];
                    Do[
                        tGStart = gammaOffsets[[gIdx]] * dTau;
                        t1Start = tGStart + compTauP;
                        rhoList = toPackedComplexMat /@ If[Length[rhoByG] === 1, rhoByG[[1]], rhoByG[[gIdx]]];
                        nStates = Length[rhoList];
                        {m1s, j1s, d1s} = splitFloquetPhaseTime[t1Start, period, dTau, nDiv];
                        fStartInv = vh1 . getA1InvAt[j1s, d1s];
                        aMat = ConstantArray[0. + 0. I, {n1, dim * dim}];
                        groupsM = <||>;

                        Do[
                            t1End = t1Start + t1Vec[[k1]];
                            {m1e, j1e, d1e} = splitFloquetPhaseTime[t1End, period, dTau, nDiv];
                            dm1 = m1e - m1s;
                            p1 = If[dm1 === 0, ConstantArray[1. + 0. I, dim], Power[lam1, dm1]];
                            fEnd = getA1At[j1e, d1e] . v1;
                            u1Prop = Transpose[p1 * Transpose[fEnd]] . fStartInv;
                            rho0 = rhoList[[Mod[k1 - 1, nStates] + 1]];
                            rho1 = u1Prop . rho0 . ConjugateTranspose[u1Prop];
                            mixStep = lookupMixStep[mixByPhase, t1End];
                            rhoM = Which[
                                mixStep["Kind"] === "Identity", rho1,
                                mixStep["Kind"] === "Unitary", mixStep["U"] . rho1 . mixStep["UH"],
                                True, ApplyHilbertMixSequence[rho1, mixStep["Steps"]]
                            ];
                            With[{phiM = t1End + compTauM},
                                pkM = phaseKeyOf[phiM];
                                deInf = buildDe2Hilbert[phiM];
                                aMat[[k1]] = Flatten[deInf["FStartInv"] . rhoM . deInf["FStart"]];
                                If[KeyExistsQ[groupsM, pkM],
                                    groupsM[pkM] = Append[groupsM[pkM], k1],
                                    groupsM[pkM] = {k1}
                                ]
                            ],
                            {k1, n1}
                        ];

                        Do[
                            deInf = buildDe2Hilbert[pk];
                            actIdx = deInf["ActiveIdx"];
                            If[actIdx =!= {},
                                With[{rows = groupsM[pk]},
                                    sMat[[rows]] = sMat[[rows]] + (aMat[[rows, actIdx]] . deInf["De2"])
                                ]
                            ],
                            {pk, Keys[groupsM]}
                        ],
                        {gIdx, nG}
                    ];
                    sMat = applyRowModulators[toPackedComplexMat[(weight / nG) * sMat], mod];
                    If[KeyExistsQ[results, key],
                        results[key] = toPackedComplexMat[results[key] + sMat],
                        results[key] = sMat
                    ]
                ],
                {comp, components}
            ]
        ],
        (* Liouville space COMPUTE2D *)
        Module[{qConjVec = toPackedComplexVec[Conjugate[qObj] / normDenom],
                x1, x1Inv, lam1, x2, x2Inv, lam2, diag1Q = AssociationQ[lDiag1],
                a1Period = Last[aSet1], a2Period = Last[aSet2], dRowCache = <||>},
            If[diag1Q, x1 = lDiag1["X"]; x1Inv = lDiag1["XInv"]; lam1 = lDiag1["Evals"]];
            If[diagMode2Q, x2 = lDiag2["X"]; x2Inv = lDiag2["XInv"]; lam2 = lDiag2["Evals"]];

            buildDe2Liouville[phiM_?NumericQ] := With[{pk = phaseKeyOf[phiM]},
                Lookup[de2Cache, Key[pk],
                    de2Cache[pk] = Module[{mM, jM, deltaM, aMInv, fMInv, dCols, k2, m2, j2, delta2,
                                           dm, powVec, dRow, actIdx, a2End},
                        {mM, jM, deltaM} = splitFloquetPhaseTime[phiM, period, dTau, nDiv];
                        aMInv = getA2InvAt[jM, deltaM];
                        fMInv = If[diagMode2Q, x2Inv . aMInv, aMInv];
                        dCols = ConstantArray[0. + 0. I, {dim, n2}];
                        Do[
                            {m2, j2, delta2} = splitFloquetPhaseTime[phiM + t2Vec[[k2]], period, dTau, nDiv];
                            dm = m2 - mM;
                            If[diagMode2Q,
                                If[Abs[delta2] <= 1.*^-11 * period,
                                    dRow = Lookup[dRowCache, Key[j2],
                                        dRowCache[j2] = qConjVec . (aSet2[[j2 + 1]] . x2)
                                    ],
                                    dRow = qConjVec . (getA2At[j2, delta2] . x2)
                                ];
                                powVec = If[dm === 0, 1., Power[lam2, dm]];
                                dCols[[All, k2]] = dRow * powVec,
                                a2End = Lookup[powCache2, dm, powCache2[dm] = MatrixPower[a2Period, dm]];
                                dCols[[All, k2]] = qConjVec . (getA2At[j2, delta2] . a2End)
                            ],
                            {k2, n2}
                        ];
                        actIdx = Pick[Range[dim], UnitStep[(Max /@ Abs[dCols]) - 1.*^-14], 1];
                        <|"FStartInv" -> fMInv, "ActiveIdx" -> actIdx,
                          "De2" -> If[actIdx === {}, None, toPackedComplexMat[dCols[[actIdx]]]]|>
                    ]
                ]
            ];

            Do[
                Module[{key = comp["Key"], weight = N[Lookup[comp, "Weight", 1.]],
                        mod = Lookup[comp, "RowModulators", None],
                        compTauP = N[Re[Lookup[comp, "PrepDuration", tauP]]],
                        compTauM = N[Re[Lookup[comp, "MixDuration", tauM]]],
                        rhoByG = comp["RhoByGamma"], mixByPhase = comp["MixByPhase"],
                        aMat, sMat, gIdx, tGStart, t1Start, t1End, rhoList, nStates,
                        m1s, j1s, d1s, m1e, j1e, d1e, dm1, p1, fStartInv, fEnd, u1Prop, a1End,
                        rho0, rho1, rhoM, mixStep, pkM, groupsM, deInf, actIdx, k1},
                    sMat = ConstantArray[0. + 0. I, {n1, n2}];
                    Do[
                        tGStart = gammaOffsets[[gIdx]] * dTau;
                        t1Start = tGStart + compTauP;
                        rhoList = toPackedComplexVec /@ If[Length[rhoByG] === 1, rhoByG[[1]], rhoByG[[gIdx]]];
                        nStates = Length[rhoList];
                        {m1s, j1s, d1s} = splitFloquetPhaseTime[t1Start, period, dTau, nDiv];
                        fStartInv = If[diag1Q, x1Inv . getA1InvAt[j1s, d1s], getA1InvAt[j1s, d1s]];
                        aMat = ConstantArray[0. + 0. I, {n1, dim}];
                        groupsM = <||>;

                        Do[
                            t1End = t1Start + t1Vec[[k1]];
                            {m1e, j1e, d1e} = splitFloquetPhaseTime[t1End, period, dTau, nDiv];
                            dm1 = m1e - m1s;
                            If[diag1Q,
                                p1 = If[dm1 === 0, ConstantArray[1. + 0. I, dim], Power[lam1, dm1]];
                                fEnd = getA1At[j1e, d1e] . x1;
                                u1Prop = Transpose[p1 * Transpose[fEnd]] . fStartInv,
                                a1End = Lookup[powCache1, dm1, powCache1[dm1] = MatrixPower[a1Period, dm1]];
                                u1Prop = getA1At[j1e, d1e] . a1End . fStartInv
                            ];
                            rho0 = rhoList[[Mod[k1 - 1, nStates] + 1]];
                            rho1 = u1Prop . rho0;
                            mixStep = lookupMixStep[mixByPhase, t1End];
                            rhoM = If[mixStep["Kind"] === "Identity", rho1, mixStep["Mat"] . rho1];
                            With[{phiM = t1End + compTauM},
                                pkM = phaseKeyOf[phiM];
                                deInf = buildDe2Liouville[phiM];
                                aMat[[k1]] = deInf["FStartInv"] . rhoM;
                                If[KeyExistsQ[groupsM, pkM],
                                    groupsM[pkM] = Append[groupsM[pkM], k1],
                                    groupsM[pkM] = {k1}
                                ]
                            ],
                            {k1, n1}
                        ];

                        Do[
                            deInf = buildDe2Liouville[pk];
                            actIdx = deInf["ActiveIdx"];
                            If[actIdx =!= {},
                                With[{rows = groupsM[pk]},
                                    sMat[[rows]] = sMat[[rows]] + (aMat[[rows, actIdx]] . deInf["De2"])
                                ]
                            ],
                            {pk, Keys[groupsM]}
                        ],
                        {gIdx, nG}
                    ];
                    sMat = applyRowModulators[toPackedComplexMat[(weight / nG) * sMat], mod];
                    If[KeyExistsQ[results, key],
                        results[key] = toPackedComplexMat[results[key] + sMat],
                        results[key] = sMat
                    ]
                ],
                {comp, components}
            ]
        ]
    ];

    <|"DataMap" -> results, "Peaks2D" -> None|>
];

(* Pure Compiled Ensemble Worker: processes an entire chunk in-place into one packed buffer per key *)
EvaluateCompiledEnsembleChunk2D[chunkPayload_Association] := Module[
    {mode = chunkPayload["Mode"], space = chunkPayload["Space"],
     n1 = chunkPayload["N1"], n2 = chunkPayload["N2"],
     keys = chunkPayload["Keys"], weights = chunkPayload["Weights"],
     nChunk = Length[chunkPayload["Weights"]],
     time1 = chunkPayload["Time1"], time2 = chunkPayload["Time2"],
     qObj = chunkPayload["QObj"], normDenom = chunkPayload["NormDenom"],
     quadMeta = chunkPayload["QuadMeta"],
     atomSpecs = chunkPayload["AtomSpecs"],
     fixedSpec = Lookup[chunkPayload, "FixedSpec", <||>],
     extractPeaksQ = TrueQ[Lookup[Lookup[chunkPayload, "FixedSpec", <||>], "ExtractPeaksOnly", False]],
     acc = <||>, bagA = <||>, bagW1 = <||>, bagW2 = <||>, buildMat, buildCompsForMember, evalMember, memberOut, w, m, k,
     applySpinZPacked, applyQuadTransform, prepBase},

    If[extractPeaksQ,
        Do[bagA[k] = Internal`Bag[]; bagW1[k] = Internal`Bag[]; bagW2[k] = Internal`Bag[], {k, keys}],
        Do[acc[k] = ConstantArray[0. + 0. I, {n1, n2}], {k, keys}]
    ];

    buildMat[atomSpec_Association, coeffsVec_List] := If[
        atomSpec["IsZero"],
        toPackedComplexMat[ConstantArray[0. + 0. I, {atomSpec["Dim"], atomSpec["Dim"]}]],
        Module[{mats = atomSpec["Matrices"], res},
            res = coeffsVec[[1]] * mats[[1]];
            Do[res = res + coeffsVec[[r]] * mats[[r]], {r, 2, Length[mats]}];
            toPackedComplexMat[res]
        ]
    ];

    applySpinZPacked[state_, rotMat_] := If[
        space === "Hilbert",
        toPackedComplexMat[rotMat . state . ConjugateTranspose[rotMat]],
        toPackedComplexVec[rotMat . state]
    ];

    applyQuadTransform[rBase_, tSpec_] := Which[
        tSpec["Kind"] === "Identity", rBase,
        tSpec["Kind"] === "RotZ", applySpinZPacked[rBase, tSpec["Mat"]],
        tSpec["Kind"] === "Filter",
            If[space === "Hilbert",
                applyHilbertTransformMask[rBase, tSpec["MaskSpec"]],
                toPackedComplexVec[tSpec["Mat"] . rBase]
            ],
        tSpec["Kind"] === "RotZThenFilter",
            With[{s1 = applySpinZPacked[rBase, tSpec["RotMat"]]},
                If[space === "Hilbert",
                    applyHilbertTransformMask[s1, tSpec["MaskSpec"]],
                    toPackedComplexVec[tSpec["FilterMat"] . s1]
                ]
            ],
        True, rBase
    ];

    buildCompsForMember[mIdx_Integer] := Module[{rhoBaseList, mSpecs = quadMeta["TemplateComps"], rB},
        (* Build invariant or member-specific RhoBase per base item, then expand quadrature templates *)
        rhoBaseList = If[quadMeta["VarPrepQ"],
            With[{fn = atomSpecs["PrepCoeffFn"], theta = atomSpecs["ThetaChunk"][[mIdx]],
                  uMats = atomSpecs["PrepMats"]},
                (* Only reached if simple variable preparation; default invariant passes BaseStates directly *)
                quadMeta["BaseStates"]
            ],
            quadMeta["BaseStates"]
        ];
        Table[
            With[{tmpl = mSpecs[[cIdx]], idxB = mSpecs[[cIdx, "BaseIdx"]],
                  stateField = If[space === "Hilbert", "RhoStates", "RhoVecs"]},
                rB = rhoBaseList[[idxB]];
                Join[
                    <|"Key" -> tmpl["Key"],
                      "Weight" -> tmpl["Weight"],
                      "RowModulators" -> tmpl["RowModulators"],
                      stateField -> (applyQuadTransform[rB, #] & /@ tmpl["StateOps"])|>,
                    tmpl["MixPart"]
                ]
            ],
            {cIdx, Length[mSpecs]}
        ]
    ];

    evalMember[mIdx_Integer] := Module[{comps = buildCompsForMember[mIdx]},
        Which[
            mode === "DoubleDiag",
                If[space === "Hilbert",
                    Module[{h2Mat, h2Diag, h1Mat, h1Diag, ev1Mode = fixedSpec["Ev1Mode"]},
                        h2Mat = If[KeyExistsQ[atomSpecs, "H2"],
                            buildMat[atomSpecs["H2"], atomSpecs["H2Coeffs"][[mIdx]]],
                            fixedSpec["H2Mat"]
                        ];
                        h2Diag = If[KeyExistsQ[fixedSpec, "H2Diag"] && !KeyExistsQ[atomSpecs, "H2"],
                            fixedSpec["H2Diag"],
                            DiagonalizeHilbertHermitian[h2Mat]
                        ];
                        If[ev1Mode === "Standard",
                            h1Mat = If[KeyExistsQ[atomSpecs, "H1"],
                                buildMat[atomSpecs["H1"], atomSpecs["H1Coeffs"][[mIdx]]],
                                If[TrueQ[fixedSpec["SameH1H2"]], h2Mat, fixedSpec["H1Mat"]]
                            ];
                            h1Diag = Which[
                                TrueQ[fixedSpec["SameH1H2"]], h2Diag,
                                Max[Abs[h1Mat - h2Mat]] == 0., h2Diag,
                                KeyExistsQ[fixedSpec, "H1Diag"] && !KeyExistsQ[atomSpecs, "H1"], fixedSpec["H1Diag"],
                                True, DiagonalizeHilbertHermitian[h1Mat]
                            ];
                            EvaluateHilbertDoubleDiag2D[<|
                                "Ev1Mode" -> "Standard",
                                "H1Diag" -> h1Diag, "H2Diag" -> h2Diag,
                                "Time1" -> time1, "Time2" -> time2,
                                "QAdjMat" -> qObj, "NormDenom" -> normDenom,
                                "Components" -> comps, "ComputePeaks" -> False,
                                "ExtractPeaksOnly" -> extractPeaksQ
                            |>],
                            Module[{h1aMat, h1bMat, h1aDiag, h1bDiag},
                                h1aMat = If[KeyExistsQ[atomSpecs, "H1a"],
                                    buildMat[atomSpecs["H1a"], atomSpecs["H1aCoeffs"][[mIdx]]],
                                    fixedSpec["H1aMat"]
                                ];
                                h1bMat = If[KeyExistsQ[atomSpecs, "H1b"],
                                    buildMat[atomSpecs["H1b"], atomSpecs["H1bCoeffs"][[mIdx]]],
                                    fixedSpec["H1bMat"]
                                ];
                                h1aDiag = If[Max[Abs[h1aMat - h2Mat]] == 0., h2Diag, DiagonalizeHilbertHermitian[h1aMat]];
                                h1bDiag = If[Max[Abs[h1bMat - h1aMat]] == 0., h1aDiag, DiagonalizeHilbertHermitian[h1bMat]];
                                EvaluateHilbertDoubleDiag2D[<|
                                    "Ev1Mode" -> ev1Mode,
                                    "CTDuration" -> Lookup[fixedSpec, "CTDuration", 0.],
                                    "H1aDiag" -> h1aDiag, "H1bDiag" -> h1bDiag,
                                    "RefocusSteps" -> Lookup[fixedSpec, "RefocusSteps", {}],
                                    "RefocusSingleUnitary" -> Lookup[fixedSpec, "RefocusSingleUnitary", None],
                                    "H2Diag" -> h2Diag,
                                    "Time1" -> time1, "Time2" -> time2,
                                    "QAdjMat" -> qObj, "NormDenom" -> normDenom,
                                    "Components" -> comps
                                |>]
                            ]
                        ]
                    ],
                    (* Liouville DoubleDiag member *)
                    Module[{l2Mat, l2Diag, l1Mat, l1Diag, ev1Mode = fixedSpec["Ev1Mode"], actQ = fixedSpec["ActQ"], actRho = fixedSpec["ActRho"], resL},
                        l2Mat = If[KeyExistsQ[atomSpecs, "L2"],
                            buildMat[atomSpecs["L2"], atomSpecs["L2Coeffs"][[mIdx]]],
                            fixedSpec["L2Mat"]
                        ];
                        l2Diag = If[KeyExistsQ[fixedSpec, "L2Diag"] && !KeyExistsQ[atomSpecs, "L2"],
                            fixedSpec["L2Diag"],
                            DiagonalizeLiouvilleBlocked[l2Mat, If[TrueQ[fixedSpec["SameL1L2"]] && ev1Mode === "Standard", Union[actQ, actRho], actQ]]
                        ];
                        If[l2Diag === $Failed,
                            (* Seamless SteppedPropagator fallback for non-diagonalizable member *)
                            With[{dt1 = fixedSpec["Dt1"], dt2 = fixedSpec["Dt2"],
                                  tmin1 = fixedSpec["Tmin1"], tmin2 = fixedSpec["Tmin2"]},
                                l1Mat = If[KeyExistsQ[atomSpecs, "L1"],
                                    buildMat[atomSpecs["L1"], atomSpecs["L1Coeffs"][[mIdx]]],
                                    If[TrueQ[fixedSpec["SameL1L2"]], l2Mat, fixedSpec["L1Mat"]]
                                ];
                                Return[EvaluateSteppedPropagator2D[<|
                                    "Space" -> "Liouville", "Ev1Mode" -> ev1Mode,
                                    "N1" -> n1, "N2" -> n2, "Time2" -> time2, "NormDenom" -> normDenom,
                                    "InitProp1a" -> If[tmin1 != 0., toPackedComplexMat[MatrixExp[l1Mat * tmin1]], None],
                                    "StepProp1a" -> toPackedComplexMat[MatrixExp[l1Mat * dt1]],
                                    "InitProp2" -> If[tmin2 != 0., toPackedComplexMat[MatrixExp[l2Mat * tmin2]], None],
                                    "StepProp2" -> toPackedComplexMat[MatrixExp[l2Mat * dt2]],
                                    "Diag2" -> None, "QObj" -> qObj, "Components" -> comps
                                |>]]
                            ]
                        ];
                        If[ev1Mode === "Standard",
                            l1Mat = If[KeyExistsQ[atomSpecs, "L1"],
                                buildMat[atomSpecs["L1"], atomSpecs["L1Coeffs"][[mIdx]]],
                                If[TrueQ[fixedSpec["SameL1L2"]], l2Mat, fixedSpec["L1Mat"]]
                            ];
                            l1Diag = Which[
                                TrueQ[fixedSpec["SameL1L2"]], l2Diag,
                                Max[Abs[l1Mat - l2Mat]] == 0., DiagonalizeLiouvilleBlocked[l1Mat, actRho],
                                KeyExistsQ[fixedSpec, "L1Diag"] && !KeyExistsQ[atomSpecs, "L1"], fixedSpec["L1Diag"],
                                True, DiagonalizeLiouvilleBlocked[l1Mat, actRho]
                            ];
                            If[l1Diag === $Failed,
                                With[{dt1 = fixedSpec["Dt1"], dt2 = fixedSpec["Dt2"],
                                      tmin1 = fixedSpec["Tmin1"], tmin2 = fixedSpec["Tmin2"]},
                                    Return[EvaluateSteppedPropagator2D[<|
                                        "Space" -> "Liouville", "Ev1Mode" -> ev1Mode,
                                        "N1" -> n1, "N2" -> n2, "Time2" -> time2, "NormDenom" -> normDenom,
                                        "InitProp1a" -> If[tmin1 != 0., toPackedComplexMat[MatrixExp[l1Mat * tmin1]], None],
                                        "StepProp1a" -> toPackedComplexMat[MatrixExp[l1Mat * dt1]],
                                        "InitProp2" -> If[tmin2 != 0., toPackedComplexMat[MatrixExp[l2Mat * tmin2]], None],
                                        "StepProp2" -> toPackedComplexMat[MatrixExp[l2Mat * dt2]],
                                        "Diag2" -> l2Diag, "QObj" -> qObj, "Components" -> comps
                                    |>]]
                                ]
                            ];
                            EvaluateLiouvilleDoubleDiag2D[<|
                                "Ev1Mode" -> "Standard",
                                "L1Diag" -> l1Diag, "L2Diag" -> l2Diag,
                                "Time1" -> time1, "Time2" -> time2,
                                "QVec" -> qObj, "NormDenom" -> normDenom,
                                "Components" -> comps
                            |>],
                            Module[{l1aMat, l1bMat, l1aDiag, l1bDiag},
                                l1aMat = If[KeyExistsQ[atomSpecs, "L1a"],
                                    buildMat[atomSpecs["L1a"], atomSpecs["L1aCoeffs"][[mIdx]]],
                                    fixedSpec["L1aMat"]
                                ];
                                l1bMat = If[KeyExistsQ[atomSpecs, "L1b"],
                                    buildMat[atomSpecs["L1b"], atomSpecs["L1bCoeffs"][[mIdx]]],
                                    fixedSpec["L1bMat"]
                                ];
                                l1aDiag = DiagonalizeLiouvilleBlocked[l1aMat, All];
                                l1bDiag = If[Max[Abs[l1bMat - l1aMat]] == 0., l1aDiag, DiagonalizeLiouvilleBlocked[l1bMat, All]];
                                EvaluateLiouvilleDoubleDiag2D[<|
                                    "Ev1Mode" -> ev1Mode,
                                    "CTDuration" -> Lookup[fixedSpec, "CTDuration", 0.],
                                    "L1aDiag" -> l1aDiag, "L1bDiag" -> l1bDiag,
                                    "RefocusMat" -> Lookup[fixedSpec, "RefocusMat", IdentityMatrix[Length[l2Mat]]],
                                    "L2Diag" -> l2Diag,
                                    "Time1" -> time1, "Time2" -> time2,
                                    "QVec" -> qObj, "NormDenom" -> normDenom,
                                    "Components" -> comps
                                |>]
                            ]
                        ]
                    ]
                ],
            mode === "COMPUTE2D",
                Module[{theta = atomSpecs["ThetaChunk"][[mIdx]], period = fixedSpec["Period"],
                        nDiv = fixedSpec["NDiv"], t0 = fixedSpec["T0"], dTau,
                        nSub = fixedSpec["NSubPerDiv"], microSub = fixedSpec["MicroSubsteps"],
                        dim = fixedSpec["Dim"], sameQ = TrueQ[fixedSpec["SameEv1Ev2"]],
                        td1Spec = atomSpecs["TD1Spec"], td2Spec = atomSpecs["TD2Spec"],
                        prepTDSpec = Lookup[atomSpecs, "PrepTDSpec", None],
                        mixTDSpec = Lookup[atomSpecs, "MixTDSpec", None],
                        prepStageOps = Lookup[fixedSpec, "PrepStageOps", None],
                        mixStageOps = Lookup[fixedSpec, "MixStageOps", None],
                        eval1Fn, eval2Fn, evalPrepFn, evalMixFn, aSet1, aSet2,
                        gammaOffsets = fixedSpec["GammaOffsets"], nG, evalStageSteps,
                        cComps},
                    dTau = period / nDiv;
                    nG = Length[gammaOffsets];
                    evalStageSteps[opsList_List, tStartStage_?NumericQ, bgFn_] := Module[
                        {tCur = N[tStartStage], accU = toPackedComplexMat[IdentityMatrix[dim]], steps = {}, uS},
                        Do[
                            Which[
                                op["Kind"] === "Unitary",
                                    accU = toPackedComplexMat[op["Mat"] . accU],
                                op["Kind"] === "Filter",
                                    If[space === "Hilbert",
                                        If[Max[Abs[accU - IdentityMatrix[dim]]] > 1.*^-13,
                                            AppendTo[steps, {"Unitary", accU, toPackedComplexMat[ConjugateTranspose[accU]]}];
                                            accU = toPackedComplexMat[IdentityMatrix[dim]]
                                        ];
                                        AppendTo[steps, {"Filter", op["MaskSpec"]}],
                                        accU = toPackedComplexMat[op["Mat"] . accU]
                                    ],
                                op["Kind"] === "StaticDelay",
                                    uS = If[space === "Hilbert",
                                        toPackedComplexMat[MatrixExp[-I * op["Mat"] * op["Dur"]]],
                                        toPackedComplexMat[MatrixExp[op["Mat"] * op["Dur"]]]
                                    ];
                                    accU = toPackedComplexMat[uS . accU];
                                    tCur += op["Dur"],
                                op["Kind"] === "DelayCF4",
                                    uS = CF4PropagatorStep[
                                        If[KeyExistsQ[op, "EvMat"],
                                            With[{em = op["EvMat"]}, Function[{tv}, bgFn[tv] + em]],
                                            bgFn
                                        ],
                                        tCur, tCur + op["Dur"], Max[2, microSub], space, dim
                                    ];
                                    accU = toPackedComplexMat[uS . accU];
                                    tCur += op["Dur"]
                            ],
                            {op, opsList}
                        ];
                        If[space === "Hilbert",
                            If[steps === {},
                                If[Max[Abs[accU - IdentityMatrix[dim]]] <= 1.*^-13,
                                    <|"Kind" -> "Identity"|>,
                                    <|"Kind" -> "Unitary", "U" -> accU, "UH" -> toPackedComplexMat[ConjugateTranspose[accU]]|>
                                ],
                                If[Max[Abs[accU - IdentityMatrix[dim]]] > 1.*^-13,
                                    AppendTo[steps, {"Unitary", accU, toPackedComplexMat[ConjugateTranspose[accU]]}]
                                ];
                                <|"Kind" -> "Steps", "Steps" -> steps|>
                            ],
                            If[Max[Abs[accU - IdentityMatrix[dim]]] <= 1.*^-13,
                                <|"Kind" -> "Identity"|>,
                                <|"Kind" -> "Mat", "Mat" -> accU|>
                            ]
                        ]
                    ];

                    eval1Fn = With[{cf = td1Spec["CoeffFn"], ms = td1Spec["Matrices"], th = theta},
                        Function[{tVal},
                            With[{cv = cf[th, tVal]},
                                Sum[cv[[r]] * ms[[r]], {r, Length[ms]}]
                            ]
                        ]
                    ];
                    aSet1 = BuildPeriodicPropagatorsCF4[eval1Fn, t0, period, nDiv, nSub, space, dim];
                    If[sameQ,
                        eval2Fn = eval1Fn;
                        aSet2 = aSet1,
                        eval2Fn = With[{cf = td2Spec["CoeffFn"], ms = td2Spec["Matrices"], th = theta},
                            Function[{tVal},
                                With[{cv = cf[th, tVal]},
                                    Sum[cv[[r]] * ms[[r]], {r, Length[ms]}]
                                ]
                            ]
                        ];
                        aSet2 = BuildPeriodicPropagatorsCF4[eval2Fn, t0, period, nDiv, nSub, space, dim]
                    ];

                    evalPrepFn = If[prepTDSpec === None, eval1Fn,
                        With[{cf = prepTDSpec["CoeffFn"], ms = prepTDSpec["Matrices"], th = theta},
                            Function[{tVal}, With[{cv = cf[th, tVal]}, Sum[cv[[r]] * ms[[r]], {r, Length[ms]}]]]
                        ]
                    ];
                    evalMixFn = If[mixTDSpec === None, eval1Fn,
                        With[{cf = mixTDSpec["CoeffFn"], ms = mixTDSpec["Matrices"], th = theta},
                            Function[{tVal}, With[{cv = cf[th, tVal]}, Sum[cv[[r]] * ms[[r]], {r, Length[ms]}]]]
                        ]
                    ];

                    cComps = Table[
                        With[{tmpl = quadMeta["TemplateComps"][[cIdx]],
                              rStatesFixed = comps[[cIdx, If[space === "Hilbert", "RhoStates", "RhoVecs"]]],
                              mixFixed = fixedSpec["MixCOMPUTE"][[cIdx]]},
                          With[{varRhoByG = If[ListQ[prepStageOps],
                                  Module[{pOps = prepStageOps[[tmpl["BaseIdx"]]], r0Init = quadMeta["Rho0"], stPrep},
                                      Table[
                                          stPrep = evalStageSteps[pOps, t0 + gammaOffsets[[gIdx]] * dTau, evalPrepFn];
                                          With[{rB = Which[
                                                  stPrep["Kind"] === "Identity", r0Init,
                                                  stPrep["Kind"] === "Unitary", toPackedComplexMat[stPrep["U"] . r0Init . stPrep["UH"]],
                                                  stPrep["Kind"] === "Mat", toPackedComplexVec[stPrep["Mat"] . r0Init],
                                                  True, ApplyHilbertMixSequence[r0Init, stPrep["Steps"]]
                                              ]},
                                              applyQuadTransform[rB, #] & /@ tmpl["StateOps"]
                                          ],
                                          {gIdx, nG}
                                      ]
                                  ],
                                  None
                              ],
                              varMixByPhase = If[ListQ[mixStageOps] && KeyExistsQ[fixedSpec, "UniqueMixStartPhases"],
                                  With[{mOps = mixStageOps[[tmpl["BaseIdx"]]], uPhases = fixedSpec["UniqueMixStartPhases"], pkFn = (Round[Mod[#, period], 1.*^-11] &)},
                                      AssociationMap[
                                          Function[{pkVal}, evalStageSteps[mOps, t0 + pkVal, evalMixFn]],
                                          pkFn /@ uPhases
                                      ]
                                  ],
                                  None
                              ]},
                            <|
                                "Key" -> tmpl["Key"],
                                "Weight" -> tmpl["Weight"],
                                "RowModulators" -> tmpl["RowModulators"],
                                "PrepDuration" -> Lookup[tmpl, "PrepDuration", fixedSpec["PrepDuration"]],
                                "MixDuration" -> Lookup[tmpl, "MixDuration", fixedSpec["MixDuration"]],
                                "RhoByGamma" -> If[varRhoByG =!= None, varRhoByG, {rStatesFixed}],
                                "MixByPhase" -> If[varMixByPhase =!= None, varMixByPhase, mixFixed]
                            |>
                          ]
                        ],
                        {cIdx, Length[comps]}
                    ];
                    EvaluatePeriodicCOMPUTE2D[<|
                        "Space" -> space,
                        "Period" -> period, "NDiv" -> nDiv, "T0" -> t0,
                        "ASet1" -> aSet1, "ASet2" -> aSet2, "SameEv1Ev2" -> sameQ,
                        "EvalGen1" -> eval1Fn, "EvalGen2" -> eval2Fn, "MicroSubsteps" -> microSub,
                        "PrepDuration" -> fixedSpec["PrepDuration"], "MixDuration" -> fixedSpec["MixDuration"],
                        "GammaOffsets" -> gammaOffsets,
                        "Time1" -> time1, "Time2" -> time2,
                        "QObj" -> qObj, "NormDenom" -> normDenom,
                        "Components" -> cComps
                    |>]
                ]
        ]
    ];

    If[extractPeaksQ,
        Do[
            w = weights[[m]];
            memberOut = evalMember[m];
            With[{pMap = memberOut["PeaksByKey"]},
                Do[
                    With[{pk = pMap[k]},
                        If[pk[[1]] =!= {},
                            Internal`StuffBag[bagA[k], toPackedComplexVec[w * pk[[1]]]];
                            Internal`StuffBag[bagW1[k], pk[[2]]];
                            Internal`StuffBag[bagW2[k], pk[[3]]]
                        ]
                    ],
                    {k, keys}
                ]
            ],
            {m, nChunk}
        ];
        Return[<|"PeakBags" -> AssociationMap[Function[{k},
            Module[{pA = Internal`BagPart[bagA[k], All], pW1 = Internal`BagPart[bagW1[k], All], pW2 = Internal`BagPart[bagW2[k], All]},
                If[pA === {},
                    {{}, {}, {}},
                    {toPackedComplexVec[If[Length[pA] === 1, pA[[1]], Join @@ pA]],
                     toPackedRealVec[If[Length[pW1] === 1, pW1[[1]], Join @@ pW1]],
                     toPackedRealVec[If[Length[pW2] === 1, pW2[[1]], Join @@ pW2]]}
                ]
            ]], keys]|>],
        Do[
            w = weights[[m]];
            memberOut = evalMember[m];
            With[{dMap = memberOut["DataMap"]},
                Do[
                    acc[k] = acc[k] + w * dMap[k],
                    {k, keys}
                ]
            ],
            {m, nChunk}
        ];
        Do[acc[k] = toPackedComplexMat[acc[k]], {k, keys}];
        <|"DataMap" -> acc|>
    ]
];

End[];
EndPackage[];
