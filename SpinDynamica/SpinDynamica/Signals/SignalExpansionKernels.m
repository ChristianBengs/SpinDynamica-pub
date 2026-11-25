(* ::Package:: *)

(* :Title: Signals`SignalExpansionKernels` *)
(* :Context: Signals`SignalExpansionKernels` *)
(* :Summary: Self-contained high-performance Lorentzian 1D/2D/streaming signal expansion kernels
   (adaptive zero-padded parallel WVM phase recurrence, lazy native LLVM FunctionCompile kernels,
   Type-1 and Type-3 ES-NUFFT with Gauss-Legendre deconvolution, Chebyshev-SVD shared-weight NUFFT
   and dyadic truncated-window hybrid for heterogeneous relaxation rates, 2D separable ES-NUFFT,
   and empirically calibrated automatic method selection). Zero external C dependencies. *)

BeginPackage["Signals`SignalExpansionKernels`",
{"Global`",
 "SpinDynamica`",
 "General`UserLevel`",
 "General`Miscellaneous`"
}];

ClearAll[
    ExpansionMethod,
    ExpansionTolerance,
    OversamplingFactor,
    ReinitStride,
    TimingOrigin,
    DropTolerance,
    PlanCache,
    CompilationBackend,
    ExpandLorentzianUniformFast,
    ExpandLorentzianTimesFast,
    ExpandLorentzian2DFast,
    CreateStreamingSignalAccumulator,
    AccumulateStreamingPeaks,
    FinalizeStreamingSignal,
    StreamingSignalStats,
    StreamingSignalHandle,
    ClearSignalExpansionPlanCache,
    GetSignalExpansionCompileStats,
    WarmupSignalExpansionKernels
];

ExpansionMethod::usage = "ExpansionMethod is an option for Signal, Signal1D, Signal2D, and SignalExpansionKernels selecting \"Legacy\" (default in Signal/Signal1D/Signal2D), \"Fast\"/Automatic (analytical peak-table preserving fast expansion), \"Streaming\" (explicit opt-in bounded-memory numeric FID streaming in Signal1D), \"DirectLLVM\", \"DirectWVM\", \"DirectCompensated\", \"DirectNonUniform\", \"NUFFT\"/\"NUFFTType1\", \"NUFFTType3\", \"HeteroExactGroups\", \"HeteroSeparable\", \"HeteroDyadicHybrid\", \"2DDirect\", \"2DSemiSeparable\", or \"2DNUFFT\".";
ExpansionTolerance::usage = "ExpansionTolerance sets target relative accuracy (default 10^-9) for NUFFT and low-rank signal expansion.";
OversamplingFactor::usage = "OversamplingFactor sets FFT grid oversampling ratio sigma >= 1.1 (default 2.0).";
ReinitStride::usage = "ReinitStride sets block length for exact Exp reinitialization in phase recurrence (default 512).";
TimingOrigin::usage = "TimingOrigin sets decay origin tmin for Exp[-lambda_k * (t - tmin)]. Defaults to First[tGrid].";
DropTolerance::usage = "DropTolerance skips peaks with |A_k| <= tol * max|A|. Default 0. (lossless).";
PlanCache::usage = "PlanCache -> True (default) caches NUFFT plans and separable SVD kernels internally.";
CompilationBackend::usage = "CompilationBackend selects the numerical kernel backend: Automatic (default: zero-JIT WVM/parallel WVM whenever the relevant LLVM kernel has not already been compiled, reusing LLVM automatically once warmed via WarmupSignalExpansionKernels[] or an explicit LLVM call), \"WVM\" (always zero-compile WVM/parallel WVM), or \"LLVM\" (lazy Wolfram FunctionCompile native kernels).";

ExpandLorentzianUniformFast::usage = "ExpandLorentzianUniformFast[peaks, tmin, dt, n, t0, opts] evaluates uniform-grid Lorentzian sums using calibrated A3/A5/A6/A10 kernels.";
ExpandLorentzianTimesFast::usage = "ExpandLorentzianTimesFast[peaks, tvals, tmin, opts] evaluates Lorentzian sums on arbitrary uniform or nonuniform tvals.";
ExpandLorentzian2DFast::usage = "ExpandLorentzian2DFast[peaks2D, t1Grid, t2Grid, opts] evaluates 2D Lorentzian sums {A, w1, l1, w2, l2} via direct BLAS-3, SemiSeparable, or 2D ES-NUFFT.";
CreateStreamingSignalAccumulator::usage = "CreateStreamingSignalAccumulator[tGrid, lam, opts] creates an O(sigma Nt) bounded-memory streaming NUFFT accumulator.";
AccumulateStreamingPeaks::usage = "AccumulateStreamingPeaks[acc, peaks, weight] appends weighted orientation peaks into a streaming NUFFT accumulator.";
FinalizeStreamingSignal::usage = "FinalizeStreamingSignal[acc] flushes remaining staged peaks, runs one FFT, deconvolves, and returns the packed complex FID.";
StreamingSignalStats::usage = "StreamingSignalStats[acc] returns diagnostics for a streaming signal accumulator.";
StreamingSignalHandle::usage = "StreamingSignalHandle[id] is an opaque handle for bounded-memory incremental NUFFT accumulation.";
ClearSignalExpansionPlanCache::usage = "ClearSignalExpansionPlanCache[] clears cached 1D/2D NUFFT plans and separable SVD kernels.";
GetSignalExpansionCompileStats::usage = "GetSignalExpansionCompileStats[] returns lazy LLVM compilation state and timings.";
WarmupSignalExpansionKernels::usage = "WarmupSignalExpansionKernels[] or WarmupSignalExpansionKernels[{keys...}] pre-compiles selected lazy LLVM FunctionCompile kernels so subsequent calls run warm.";

Begin["`Private`"];

kernelDefaults = {
    ExpansionMethod -> Automatic,
    ExpansionTolerance -> 1.*^-9,
    OversamplingFactor -> 2.0,
    ReinitStride -> 512,
    TimingOrigin -> Automatic,
    DropTolerance -> 0.,
    PlanCache -> True,
    CompilationBackend -> Automatic,
    MaxBufferMB -> 32.,
    ParallelExecution -> True,
    NumChunks -> Automatic,
    SlowThresholdD -> 16.,
    BatchPeakCap -> 16384
};

(* ====================================================================== *)
(* Internal Caches & Input Normalization                                  *)
(* ====================================================================== *)
planCache1D = <||>;
planCache2D = <||>;
separableCache = <||>;
nextStreamId = 0;
streamStore = <||>;

ClearSignalExpansionPlanCache[] := (
    planCache1D = <||>;
    planCache2D = <||>;
    separableCache = <||>;
);

normalizePeaks1D[peaks_List, dropTol_?NumericQ, tOrigin_: 0.] := Module[
    {nP = Length[peaks], amps, omegasRaw, omegas, imW, lams, maxA, keep},
    If[nP === 0, Return[{{}, {}, {}}]];
    amps      = Developer`ToPackedArray[N[peaks[[All, 1]]] + 0. I];
    omegasRaw = N[peaks[[All, 2]]];
    omegas    = Developer`ToPackedArray[N[Re[omegasRaw]]];
    imW       = Developer`ToPackedArray[N[Im[omegasRaw]]];
    lams      = If[Length[First[peaks]] >= 3,
        Developer`ToPackedArray[N[Re[peaks[[All, 3]]]]],
        Developer`ToPackedArray[ConstantArray[0., nP]]
    ];
    If[Max[Abs[imW]] > 0.,
        If[N[tOrigin] != 0., amps = Developer`ToPackedArray[amps * Exp[-imW * N[tOrigin]]]];
        lams = Developer`ToPackedArray[lams + imW];
    ];
    maxA = Max[Abs[amps]];
    If[maxA <= 0., Return[{{}, {}, {}}]];
    If[dropTol > 0.,
        keep = Flatten[Position[UnitStep[Abs[amps] - N[dropTol] * maxA], 1, {1}, Heads -> False]];
        If[keep === {}, Return[{{}, {}, {}}]];
        amps = amps[[keep]]; omegas = omegas[[keep]]; lams = lams[[keep]];
    ];
    {amps, omegas, lams}
];

parsePeaks2DLocal[peaks2D_List, dropTol_?NumericQ] := Module[
    {nP = Length[peaks2D], amps, w1, l1, w2, l2, nCol, maxA, keep},
    If[nP === 0, Return[{{}, {}, {}, {}, {}}]];
    nCol = Length[First[peaks2D]];
    amps = Developer`ToPackedArray[N[peaks2D[[All, 1]]] + 0. I];
    w1   = Developer`ToPackedArray[N[peaks2D[[All, 2]]]];
    l1   = Developer`ToPackedArray[N[Re[peaks2D[[All, 3]]]]];
    w2   = Developer`ToPackedArray[N[peaks2D[[All, 4]]]];
    l2   = If[nCol >= 5,
        Developer`ToPackedArray[N[Re[peaks2D[[All, 5]]]]],
        Developer`ToPackedArray[ConstantArray[0., nP]]
    ];
    maxA = Max[Abs[amps]];
    If[maxA <= 0., Return[{{}, {}, {}, {}, {}}]];
    If[dropTol > 0.,
        keep = Flatten[Position[UnitStep[Abs[amps] - N[dropTol] * maxA], 1, {1}, Heads -> False]];
        If[keep === {}, Return[{{}, {}, {}, {}, {}}]];
        amps = amps[[keep]]; w1 = w1[[keep]]; l1 = l1[[keep]]; w2 = w2[[keep]]; l2 = l2[[keep]];
    ];
    {amps, w1, l1, w2, l2}
];

detectUniformGridLocal[tVec_List] := Module[{nT = Length[tVec], dt, span, diffs},
    If[nT < 2, Return[{nT === 1, If[nT >= 1, tVec[[1]], 0.], If[nT === 1, 1., 0.]}]];
    dt = tVec[[2]] - tVec[[1]];
    If[dt == 0., Return[{False, tVec[[1]], 0.}]];
    span = Abs[tVec[[-1]] - tVec[[1]]];
    diffs = Differences[tVec] - dt;
    If[Max[Abs[diffs]] <= 1.*^-11 * Max[Abs[dt], span / nT, 1.*^-15],
        {True, tVec[[1]], dt},
        {False, tVec[[1]], 0.}
    ]
];

padRect[vec_List, cSize_Integer, nChunks_Integer, padVal_: 0.] := Module[{tot = cSize * nChunks, padLen},
    padLen = tot - Length[vec];
    Partition[If[padLen > 0, Join[vec, Developer`ToPackedArray[ConstantArray[padVal, padLen]]], vec], cSize]
];

chooseChunkCountLocal[nP_Integer, nT_Integer, parOpt_, nChunksOpt_] := Module[{work, nCores},
    If[parOpt === False || (IntegerQ[nChunksOpt] && nChunksOpt <= 1) || nP < 512, Return[1]];
    If[IntegerQ[nChunksOpt] && nChunksOpt >= 2, Return[Min[nChunksOpt, nP]]];
    work = nP * nT;
    nCores = Max[1, $ProcessorCount];
    If[parOpt === Automatic && work < 400000, Return[1]];
    Min[nCores, Max[2, Quotient[nP, 256]]]
];

(* ====================================================================== *)
(* WVM Compiled Direct Kernels (Zero MainEvaluate)                        *)
(* ====================================================================== *)
wvmTimeBlockedKernel = Compile[
    {{a0Re, _Real, 1}, {a0Im, _Real, 1}, {wVec, _Real, 1}, {lVec, _Real, 1},
     {zRe, _Real, 1}, {zIm, _Real, 1}, {dt, _Real}, {nT, _Integer}, {stride, _Integer}},
    Module[{outRe, outIm, bufRe, bufIm, nP = Length[a0Re], blkStart, blkEnd, blkLen, dtBlk,
            p, k, dec, ph, cd, sd, aRe, aIm, cRe, cIm, zr, zi, nxtRe},
        outRe = Table[0., {nT}]; outIm = Table[0., {nT}];
        bufRe = Table[0., {stride}]; bufIm = Table[0., {stride}];
        For[blkStart = 1, blkStart <= nT, blkStart += stride,
            blkEnd = Min[blkStart + stride - 1, nT];
            blkLen = blkEnd - blkStart + 1;
            dtBlk = (blkStart - 1) * dt;
            For[k = 1, k <= blkLen, k++, bufRe[[k]] = 0.; bufIm[[k]] = 0.];
            For[p = 1, p <= nP, p++,
                dec = If[blkStart == 1, 1., Exp[-lVec[[p]] * dtBlk]];
                ph  = If[blkStart == 1, 0., wVec[[p]] * dtBlk];
                cd  = dec * Cos[ph]; sd = dec * Sin[ph];
                aRe = a0Re[[p]]; aIm = a0Im[[p]];
                cRe = aRe * cd - aIm * sd; cIm = aRe * sd + aIm * cd;
                zr = zRe[[p]]; zi = zIm[[p]];
                bufRe[[1]] += cRe; bufIm[[1]] += cIm;
                For[k = 2, k <= blkLen, k++,
                    nxtRe = cRe * zr - cIm * zi; cIm = cRe * zi + cIm * zr; cRe = nxtRe;
                    bufRe[[k]] += cRe; bufIm[[k]] += cIm
                ]
            ];
            For[k = 1, k <= blkLen, k++,
                outRe[[blkStart + k - 1]] = bufRe[[k]]; outIm[[blkStart + k - 1]] = bufIm[[k]]
            ]
        ];
        {outRe, outIm}
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
    RuntimeAttributes -> {Listable}, Parallelization -> True
];

wvmCompensatedKernel = Compile[
    {{a0Re, _Real, 1}, {a0Im, _Real, 1}, {wVec, _Real, 1}, {lVec, _Real, 1},
     {zRe, _Real, 1}, {zIm, _Real, 1}, {dt, _Real}, {nT, _Integer}, {stride, _Integer}},
    Module[{outRe, outIm, sumRe, compRe, sumIm, compIm, nP = Length[a0Re],
            blkStart, blkEnd, blkLen, dtBlk, p, k, dec, ph, cd, sd, aRe, aIm,
            cRe, cIm, zr, zi, nxtRe, tVal},
        outRe = Table[0., {nT}]; outIm = Table[0., {nT}];
        sumRe = Table[0., {stride}]; compRe = Table[0., {stride}];
        sumIm = Table[0., {stride}]; compIm = Table[0., {stride}];
        For[blkStart = 1, blkStart <= nT, blkStart += stride,
            blkEnd = Min[blkStart + stride - 1, nT];
            blkLen = blkEnd - blkStart + 1;
            dtBlk = (blkStart - 1) * dt;
            For[k = 1, k <= blkLen, k++,
                sumRe[[k]] = 0.; compRe[[k]] = 0.; sumIm[[k]] = 0.; compIm[[k]] = 0.
            ];
            For[p = 1, p <= nP, p++,
                dec = If[blkStart == 1, 1., Exp[-lVec[[p]] * dtBlk]];
                ph  = If[blkStart == 1, 0., wVec[[p]] * dtBlk];
                cd  = dec * Cos[ph]; sd = dec * Sin[ph];
                aRe = a0Re[[p]]; aIm = a0Im[[p]];
                cRe = aRe * cd - aIm * sd; cIm = aRe * sd + aIm * cd;
                zr = zRe[[p]]; zi = zIm[[p]];
                tVal = sumRe[[1]] + cRe;
                If[Abs[sumRe[[1]]] >= Abs[cRe], compRe[[1]] += (sumRe[[1]] - tVal) + cRe, compRe[[1]] += (cRe - tVal) + sumRe[[1]]];
                sumRe[[1]] = tVal;
                tVal = sumIm[[1]] + cIm;
                If[Abs[sumIm[[1]]] >= Abs[cIm], compIm[[1]] += (sumIm[[1]] - tVal) + cIm, compIm[[1]] += (cIm - tVal) + sumIm[[1]]];
                sumIm[[1]] = tVal;
                For[k = 2, k <= blkLen, k++,
                    nxtRe = cRe * zr - cIm * zi; cIm = cRe * zi + cIm * zr; cRe = nxtRe;
                    tVal = sumRe[[k]] + cRe;
                    If[Abs[sumRe[[k]]] >= Abs[cRe], compRe[[k]] += (sumRe[[k]] - tVal) + cRe, compRe[[k]] += (cRe - tVal) + sumRe[[k]]];
                    sumRe[[k]] = tVal;
                    tVal = sumIm[[k]] + cIm;
                    If[Abs[sumIm[[k]]] >= Abs[cIm], compIm[[k]] += (sumIm[[k]] - tVal) + cIm, compIm[[k]] += (cIm - tVal) + sumIm[[k]]];
                    sumIm[[k]] = tVal
                ]
            ];
            For[k = 1, k <= blkLen, k++,
                outRe[[blkStart + k - 1]] = sumRe[[k]] + compRe[[k]];
                outIm[[blkStart + k - 1]] = sumIm[[k]] + compIm[[k]]
            ]
        ];
        {outRe, outIm}
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
    RuntimeAttributes -> {Listable}, Parallelization -> True
];

wvmNonUniformKernel = Compile[
    {{aRe, _Real, 1}, {aIm, _Real, 1}, {wVec, _Real, 1}, {lVec, _Real, 1},
     {tVec, _Real, 1}, {dtOrigin, _Real, 1}},
    Module[{nP = Length[aRe], nT = Length[tVec], outRe, outIm, p, k, ar, ai, w, l, d, ph, cd, sd},
        outRe = Table[0., {nT}]; outIm = Table[0., {nT}];
        For[p = 1, p <= nP, p++,
            ar = aRe[[p]]; ai = aIm[[p]];
            If[ar == 0. && ai == 0., Continue[]];
            w = wVec[[p]]; l = lVec[[p]];
            For[k = 1, k <= nT, k++,
                ph = w * tVec[[k]];
                d = Exp[-l * dtOrigin[[k]]];
                cd = d * Cos[ph]; sd = d * Sin[ph];
                outRe[[k]] += ar * cd - ai * sd;
                outIm[[k]] += ar * sd + ai * cd
            ]
        ];
        {outRe, outIm}
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
    RuntimeAttributes -> {Listable}, Parallelization -> True
];

wvmCountUniqueUpTo = Compile[
    {{lVec, _Real, 1}, {tol, _Real}, {cap, _Integer}},
    Module[{n = Length[lVec], uniq, count = 0, i, j, val, found},
        If[n == 0, Return[0]];
        uniq = Table[0., {cap + 1}];
        count = 1; uniq[[1]] = lVec[[1]];
        For[i = 2, i <= n, i++,
            val = lVec[[i]]; found = False;
            For[j = 1, j <= count, j++,
                If[Abs[val - uniq[[j]]] <= tol, found = True; Break[]]
            ];
            If[!found,
                If[count >= cap, Return[cap + 1]];
                count++; uniq[[count]] = val
            ]
        ];
        count
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed"
];

sortedLambdaRunsWVM = Compile[
    {{lSorted, _Real, 1}, {lTol, _Real}},
    Module[{n = Length[lSorted], starts, ends, reps, g = 0, anchor = 0., i},
        starts = Table[0, {n}]; ends = Table[0, {n}]; reps = Table[0., {n}];
        If[n == 0, Return[{{0.}, {0.}, {0.}}]];
        g = 1; anchor = lSorted[[1]]; starts[[1]] = 1; reps[[1]] = anchor;
        For[i = 2, i <= n, i++,
            If[Abs[lSorted[[i]] - anchor] > lTol,
                ends[[g]] = i - 1; g++; anchor = lSorted[[i]];
                starts[[g]] = i; reps[[g]] = anchor
            ]
        ];
        ends[[g]] = n;
        {N[Take[starts, g]], N[Take[ends, g]], Take[reps, g]}
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed"
];

(* ====================================================================== *)
(* WVM NUFFT, Hetero, & 2D Kernels (Zero MainEvaluate)                    *)
(* ====================================================================== *)
wvmSpreadESKernel = Compile[
    {{xVec, _Real, 1}, {cRe, _Real, 1}, {cIm, _Real, 1}, {nFFT, _Integer}, {w, _Integer}, {beta, _Real}},
    Module[{nP = Length[xVec], gridRe, gridIm, invHalfW, p, ar, ai, x, j0, r0, s, r, z, arg, wt, idx, eb},
        gridRe = Table[0., {nFFT}]; gridIm = Table[0., {nFFT}];
        invHalfW = 2. / w; eb = Exp[-beta];
        For[p = 1, p <= nP, p++,
            ar = cRe[[p]]; ai = cIm[[p]];
            If[ar == 0. && ai == 0., Continue[]];
            x = xVec[[p]]; j0 = Ceiling[x - 0.5 * w]; r0 = j0 - x;
            If[j0 < 0, j0 += nFFT]; If[j0 >= nFFT, j0 -= nFFT];
            For[s = 0, s < w, s++,
                z = (r0 + s) * invHalfW; arg = 1. - z * z;
                wt = If[arg > 0., Exp[beta * (Sqrt[arg] - 1.)], eb];
                idx = j0 + s + 1; If[idx > nFFT, idx -= nFFT];
                gridRe[[idx]] += ar * wt; gridIm[[idx]] += ai * wt
            ]
        ];
        {gridRe, gridIm}
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
    RuntimeAttributes -> {Listable}, Parallelization -> True
];

wvmSpreadPolyKernel = Compile[
    {{xVec, _Real, 1}, {cRe, _Real, 1}, {cIm, _Real, 1}, {nFFT, _Integer}, {w, _Integer},
     {d1, _Integer}, {polyCoeffs, _Real, 1}},
    Module[{nP = Length[xVec], gridRe, gridIm, p, ar, ai, x, xLeft, j0, u, s, base, val, q, idx},
        gridRe = Table[0., {nFFT}]; gridIm = Table[0., {nFFT}];
        For[p = 1, p <= nP, p++,
            ar = cRe[[p]]; ai = cIm[[p]];
            If[ar == 0. && ai == 0., Continue[]];
            x = xVec[[p]]; xLeft = x - 0.5 * w; j0 = Ceiling[xLeft]; u = j0 - xLeft;
            If[u < 0., u = 0.]; If[u > 1., u = 1.];
            If[j0 < 0, j0 += nFFT]; If[j0 >= nFFT, j0 -= nFFT];
            For[s = 0, s < w, s++,
                base = s * d1; val = polyCoeffs[[base + 1]];
                For[q = 2, q <= d1, q++, val = val * u + polyCoeffs[[base + q]]];
                idx = j0 + s + 1; If[idx > nFFT, idx -= nFFT];
                gridRe[[idx]] += ar * val; gridIm[[idx]] += ai * val
            ]
        ];
        {gridRe, gridIm}
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
    RuntimeAttributes -> {Listable}, Parallelization -> True
];

wvmSpreadPolyMultiKernel = Compile[
    {{xVec, _Real, 1}, {aReFlat, _Real, 1}, {aImFlat, _Real, 1}, {nFFT, _Integer},
     {w, _Integer}, {d1, _Integer}, {polyCoeffs, _Real, 1}, {nChan, _Integer}},
    Module[{nP = Length[xVec], totalGrid, gRe, gIm, p, x, xLeft, j0, u, s, base, val, q, idx,
            ch, pBase, gOff, wtVec},
        totalGrid = nChan * nFFT; gRe = Table[0., {totalGrid}]; gIm = Table[0., {totalGrid}]; wtVec = Table[0., {w}];
        For[p = 1, p <= nP, p++,
            x = xVec[[p]]; xLeft = x - 0.5 * w; j0 = Ceiling[xLeft]; u = j0 - xLeft;
            If[u < 0., u = 0.]; If[u > 1., u = 1.];
            If[j0 < 0, j0 += nFFT]; If[j0 >= nFFT, j0 -= nFFT];
            For[s = 0, s < w, s++,
                base = s * d1; val = polyCoeffs[[base + 1]];
                For[q = 2, q <= d1, q++, val = val * u + polyCoeffs[[base + q]]];
                wtVec[[s + 1]] = val
            ];
            pBase = (p - 1) * nChan;
            For[ch = 0, ch < nChan, ch++,
                gOff = ch * nFFT;
                For[s = 0, s < w, s++,
                    idx = j0 + s + 1; If[idx > nFFT, idx -= nFFT];
                    val = wtVec[[s + 1]];
                    gRe[[gOff + idx]] += aReFlat[[pBase + ch + 1]] * val;
                    gIm[[gOff + idx]] += aImFlat[[pBase + ch + 1]] * val
                ]
            ]
        ];
        {gRe, gIm}
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
    RuntimeAttributes -> {Listable}, Parallelization -> True
];

wvmGatherPolyKernel = Compile[
    {{yVec, _Real, 1}, {gRe, _Real, 1}, {gIm, _Real, 1}, {nFFT, _Integer},
     {w, _Integer}, {d1, _Integer}, {polyCoeffs, _Real, 1}},
    Module[{nY = Length[yVec], outRe, outIm, i, y, yLeft, j0, u, accRe, accIm, s, base, val, q, idx},
        outRe = Table[0., {nY}]; outIm = Table[0., {nY}];
        For[i = 1, i <= nY, i++,
            y = yVec[[i]]; yLeft = y - 0.5 * w; j0 = Ceiling[yLeft]; u = j0 - yLeft;
            If[u < 0., u = 0.]; If[u > 1., u = 1.];
            If[j0 < 0, j0 += nFFT]; If[j0 >= nFFT, j0 -= nFFT];
            accRe = 0.; accIm = 0.;
            For[s = 0, s < w, s++,
                base = s * d1; val = polyCoeffs[[base + 1]];
                For[q = 2, q <= d1, q++, val = val * u + polyCoeffs[[base + q]]];
                idx = j0 + s + 1; If[idx > nFFT, idx -= nFFT];
                accRe += gRe[[idx]] * val; accIm += gIm[[idx]] * val
            ];
            outRe[[i]] = accRe; outIm[[i]] = accIm
        ];
        {outRe, outIm}
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed"
];

wvmEvalChebBasis = Compile[
    {{xVec, _Real, 1}, {chebCoeffs, _Real, 2}},
    Module[{nP = Length[xVec], qEff = Length[chebCoeffs], r = Length[First[chebCoeffs]],
            out, p, k, j, z, twoZ, t0, t1, tNew},
        out = Table[0., {nP}, {r}];
        For[p = 1, p <= nP, p++,
            z = 2. * xVec[[p]] - 1.; twoZ = 2. * z; t0 = 1.;
            For[j = 1, j <= r, j++, out[[p, j]] += chebCoeffs[[1, j]]];
            If[qEff >= 2,
                t1 = z;
                For[j = 1, j <= r, j++, out[[p, j]] += t1 * chebCoeffs[[2, j]]];
                For[k = 3, k <= qEff, k++,
                    tNew = twoZ * t1 - t0; t0 = t1; t1 = tNew;
                    For[j = 1, j <= r, j++, out[[p, j]] += tNew * chebCoeffs[[k, j]]]
                ]
            ]
        ];
        out
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed"
];

wvmTruncatedDirectKernel = Compile[
    {{a0Re, _Real, 1}, {a0Im, _Real, 1}, {zRe, _Real, 1}, {zIm, _Real, 1}, {nActVec, _Integer, 1}, {nT, _Integer}},
    Module[{outRe, outIm, nP = Length[a0Re], p, k, nAct, cRe, cIm, zr, zi, nxtRe},
        outRe = Table[0., {nT}]; outIm = Table[0., {nT}];
        For[p = 1, p <= nP, p++,
            cRe = a0Re[[p]]; cIm = a0Im[[p]]; nAct = nActVec[[p]];
            If[nAct > nT, nAct = nT];
            If[nAct < 1 || (cRe == 0. && cIm == 0.), Continue[]];
            zr = zRe[[p]]; zi = zIm[[p]];
            outRe[[1]] += cRe; outIm[[1]] += cIm;
            For[k = 2, k <= nAct, k++,
                nxtRe = cRe * zr - cIm * zi; cIm = cRe * zi + cIm * zr; cRe = nxtRe;
                outRe[[k]] += cRe; outIm[[k]] += cIm
            ]
        ];
        {outRe, outIm}
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed"
];

wvmSpreadES2DKernel = Compile[
    {{x1Vec, _Real, 1}, {x2Vec, _Real, 1}, {cRe, _Real, 1}, {cIm, _Real, 1},
     {nFFT1, _Integer}, {w1, _Integer}, {beta1, _Real},
     {nFFT2, _Integer}, {w2, _Integer}, {beta2, _Real}},
    Module[{nP = Length[x1Vec], total = nFFT1 * nFFT2, gRe, gIm, wt1, wt2, idx1, idx2,
            invHW1 = 2. / w1, invHW2 = 2. / w2, eb1 = Exp[-beta1], eb2 = Exp[-beta2],
            p, ar, ai, x1, x2, j01, j02, r01, r02, s1, s2, z, arg, rowOff, u1r, u1i, id},
        gRe = Table[0., {total}]; gIm = Table[0., {total}];
        wt1 = Table[0., {w1}]; idx1 = Table[1, {w1}];
        wt2 = Table[0., {w2}]; idx2 = Table[1, {w2}];
        For[p = 1, p <= nP, p++,
            ar = cRe[[p]]; ai = cIm[[p]];
            If[ar == 0. && ai == 0., Continue[]];
            x1 = x1Vec[[p]]; j01 = Ceiling[x1 - 0.5 * w1]; r01 = j01 - x1;
            If[j01 < 0, j01 += nFFT1]; If[j01 >= nFFT1, j01 -= nFFT1];
            For[s1 = 0, s1 < w1, s1++,
                z = (r01 + s1) * invHW1; arg = 1. - z * z;
                wt1[[s1 + 1]] = If[arg > 0., Exp[beta1 * (Sqrt[arg] - 1.)], eb1];
                id = j01 + s1 + 1; If[id > nFFT1, id -= nFFT1]; idx1[[s1 + 1]] = id
            ];
            x2 = x2Vec[[p]]; j02 = Ceiling[x2 - 0.5 * w2]; r02 = j02 - x2;
            If[j02 < 0, j02 += nFFT2]; If[j02 >= nFFT2, j02 -= nFFT2];
            For[s2 = 0, s2 < w2, s2++,
                z = (r02 + s2) * invHW2; arg = 1. - z * z;
                wt2[[s2 + 1]] = If[arg > 0., Exp[beta2 * (Sqrt[arg] - 1.)], eb2];
                id = j02 + s2 + 1; If[id > nFFT2, id -= nFFT2]; idx2[[s2 + 1]] = id
            ];
            For[s1 = 1, s1 <= w1, s1++,
                rowOff = (idx1[[s1]] - 1) * nFFT2; u1r = ar * wt1[[s1]]; u1i = ai * wt1[[s1]];
                For[s2 = 1, s2 <= w2, s2++,
                    id = rowOff + idx2[[s2]];
                    gRe[[id]] += u1r * wt2[[s2]]; gIm[[id]] += u1i * wt2[[s2]]
                ]
            ]
        ];
        {gRe, gIm}
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
    RuntimeAttributes -> {Listable}, Parallelization -> True
];

(* ====================================================================== *)
(* Lazy Native LLVM FunctionCompile Kernels (Zero Startup Overhead)       *)
(* ====================================================================== *)
llvmCache = <|
    "Direct" -> None, "SpreadES" -> None, "SpreadPoly" -> None,
    "SpreadPolyMulti" -> None, "GatherPoly" -> None,
    "EvalCheb" -> None, "TruncatedDirect" -> None, "SpreadES2D" -> None
|>;
llvmCompileTimes = <||>;

shouldUseLLVM[key_String, backendOpt_, workEst_?NumericQ, coldThreshold_?NumericQ] := (
    backendOpt =!= "WVM" && (
        backendOpt === "LLVM" ||
        TrueQ[Head[llvmCache[key]] === CompiledCodeFunction]
    )
);

spreadESWVMLocal[xVec_List, cRe_List, cIm_List, nFFT_Integer, w_Integer, beta_?NumericQ, parOpt_, nChunksOpt_] := Module[
    {nP = Length[xVec], nC, cSize, res},
    nC = If[parOpt === False || nP < 8000 || nFFT * Max[2, $ProcessorCount] > 16000000,
        1,
        chooseChunkCountLocal[nP, Max[16, w], parOpt, nChunksOpt]
    ];
    If[nC > 1,
        cSize = Ceiling[nP / nC];
        res = wvmSpreadESKernel[padRect[xVec, cSize, nC], padRect[cRe, cSize, nC], padRect[cIm, cSize, nC],
            ConstantArray[nFFT, nC], ConstantArray[w, nC], ConstantArray[N[beta], nC]];
        Developer`ToPackedArray[Total[res[[All, 1]]] + I * Total[res[[All, 2]]]],
        res = wvmSpreadESKernel[xVec, cRe, cIm, nFFT, w, N[beta]];
        Developer`ToPackedArray[res[[1]] + I * res[[2]]]
    ]
];

spreadES2DWVMLocal[x1Vec_List, x2Vec_List, cRe_List, cIm_List,
                   nFFT1_Integer, w1_Integer, beta1_?NumericQ,
                   nFFT2_Integer, w2_Integer, beta2_?NumericQ, parOpt_] := Module[
    {nP = Length[x1Vec], nC, cSize, res},
    nC = If[parOpt === False || nP < 8000 || nFFT1 * nFFT2 * Max[2, $ProcessorCount] > 16000000,
        1,
        Min[Max[2, Quotient[$ProcessorCount, 2]], Max[2, Quotient[nP, 2000]]]
    ];
    If[nC > 1,
        cSize = Ceiling[nP / nC];
        res = wvmSpreadES2DKernel[padRect[x1Vec, cSize, nC], padRect[x2Vec, cSize, nC],
            padRect[cRe, cSize, nC], padRect[cIm, cSize, nC],
            ConstantArray[nFFT1, nC], ConstantArray[w1, nC], ConstantArray[N[beta1], nC],
            ConstantArray[nFFT2, nC], ConstantArray[w2, nC], ConstantArray[N[beta2], nC]];
        Developer`ToPackedArray[Total[res[[All, 1]]] + I * Total[res[[All, 2]]]],
        res = wvmSpreadES2DKernel[x1Vec, x2Vec, cRe, cIm, nFFT1, w1, N[beta1], nFFT2, w2, N[beta2]];
        Developer`ToPackedArray[res[[1]] + I * res[[2]]]
    ]
];

compileLazyLLVM[key_String, expr_] := Module[{t0, fn},
    If[llvmCache[key] =!= None, Return[llvmCache[key]]];
    t0 = AbsoluteTime[];
    fn = Quiet[Check[FunctionCompile[expr, ProgressReporting -> False], $Failed]];
    llvmCompileTimes[key] = AbsoluteTime[] - t0;
    llvmCache[key] = If[Head[fn] === CompiledCodeFunction, fn, $Failed];
    llvmCache[key]
];

getLLVMDirect[] := compileLazyLLVM["Direct",
    Function[{
        Typed[a0Re, "PackedArray"::["Real64", 1]], Typed[a0Im, "PackedArray"::["Real64", 1]],
        Typed[wVec, "PackedArray"::["Real64", 1]], Typed[lVec, "PackedArray"::["Real64", 1]],
        Typed[zRe,  "PackedArray"::["Real64", 1]], Typed[zIm,  "PackedArray"::["Real64", 1]],
        Typed[dt, "Real64"], Typed[nT, "MachineInteger"], Typed[stride, "MachineInteger"]
    },
    Module[{nP = Length[a0Re], outRe = Table[0.0, {nT}], outIm = Table[0.0, {nT}],
            bufRe = Table[0.0, {stride}], bufIm = Table[0.0, {stride}],
            blkStart = 1, blkEnd = 1, blkLen = 1, dtBlk = 0.0,
            p = 1, k = 1, dec = 0.0, ph = 0.0, cd = 0.0, sd = 0.0,
            aRe = 0.0, aIm = 0.0, cRe = 0.0, cIm = 0.0, zr = 0.0, zi = 0.0, nxtRe = 0.0},
        While[blkStart <= nT,
            blkEnd = blkStart + stride - 1; If[blkEnd > nT, blkEnd = nT];
            blkLen = blkEnd - blkStart + 1; dtBlk = (blkStart - 1) * dt;
            Do[bufRe[[k]] = 0.0; bufIm[[k]] = 0.0, {k, 1, blkLen}];
            Do[
                If[blkStart == 1, dec = 1.0; cd = 1.0; sd = 0.0,
                    dec = Exp[-lVec[[p]] * dtBlk]; ph = wVec[[p]] * dtBlk;
                    cd = dec * Cos[ph]; sd = dec * Sin[ph]];
                aRe = a0Re[[p]]; aIm = a0Im[[p]];
                cRe = aRe * cd - aIm * sd; cIm = aRe * sd + aIm * cd;
                zr = zRe[[p]]; zi = zIm[[p]];
                bufRe[[1]] += cRe; bufIm[[1]] += cIm;
                Do[
                    nxtRe = cRe * zr - cIm * zi; cIm = cRe * zi + cIm * zr; cRe = nxtRe;
                    bufRe[[k]] += cRe; bufIm[[k]] += cIm,
                    {k, 2, blkLen}
                ],
                {p, 1, nP}
            ];
            Do[outRe[[blkStart + k - 1]] = bufRe[[k]]; outIm[[blkStart + k - 1]] = bufIm[[k]], {k, 1, blkLen}];
            blkStart += stride
        ];
        outRe + I * outIm
    ]]
];

getLLVMSpreadES[] := compileLazyLLVM["SpreadES",
    Function[{
        Typed[xVec, "PackedArray"::["Real64", 1]], Typed[cRe, "PackedArray"::["Real64", 1]],
        Typed[cIm, "PackedArray"::["Real64", 1]], Typed[nFFT, "MachineInteger"],
        Typed[w, "MachineInteger"], Typed[beta, "Real64"]
    },
    Module[{nP = Length[xVec], gridRe = Table[0.0, {nFFT}], gridIm = Table[0.0, {nFFT}],
            invHalfW = 2.0 / w, halfW = 0.5 * w, eb = Exp[-beta],
            p = 1, s = 0, ar = 0.0, ai = 0.0, x = 0.0, j0 = 0, r0 = 0.0,
            z = 0.0, arg = 0.0, wt = 0.0, idx = 1},
        Do[
            ar = cRe[[p]]; ai = cIm[[p]];
            If[ar != 0.0 || ai != 0.0,
                x = xVec[[p]]; j0 = Ceiling[x - halfW]; r0 = j0 - x;
                If[j0 < 0, j0 += nFFT]; If[j0 >= nFFT, j0 -= nFFT];
                If[j0 + w <= nFFT,
                    Do[
                        z = (r0 + s) * invHalfW; arg = 1.0 - z * z;
                        wt = If[arg > 0.0, Exp[beta * (Sqrt[arg] - 1.0)], eb];
                        idx = j0 + s + 1; gridRe[[idx]] += ar * wt; gridIm[[idx]] += ai * wt,
                        {s, 0, w - 1}
                    ],
                    Do[
                        z = (r0 + s) * invHalfW; arg = 1.0 - z * z;
                        wt = If[arg > 0.0, Exp[beta * (Sqrt[arg] - 1.0)], eb];
                        idx = j0 + s + 1; If[idx > nFFT, idx -= nFFT];
                        gridRe[[idx]] += ar * wt; gridIm[[idx]] += ai * wt,
                        {s, 0, w - 1}
                    ]
                ]
            ],
            {p, 1, nP}
        ];
        gridRe + I * gridIm
    ]]
];

getLLVMSpreadPoly[] := compileLazyLLVM["SpreadPoly",
    Function[{
        Typed[xVec, "PackedArray"::["Real64", 1]], Typed[cRe, "PackedArray"::["Real64", 1]],
        Typed[cIm, "PackedArray"::["Real64", 1]], Typed[nFFT, "MachineInteger"],
        Typed[w, "MachineInteger"], Typed[d1, "MachineInteger"], Typed[polyCoeffs, "PackedArray"::["Real64", 1]]
    },
    Module[{nP = Length[xVec], gridRe = Table[0.0, {nFFT}], gridIm = Table[0.0, {nFFT}],
            halfW = 0.5 * w, p = 1, s = 0, q = 2, base = 0,
            ar = 0.0, ai = 0.0, x = 0.0, xLeft = 0.0, j0 = 0, u = 0.0, val = 0.0, idx = 1},
        Do[
            ar = cRe[[p]]; ai = cIm[[p]];
            If[ar != 0.0 || ai != 0.0,
                x = xVec[[p]]; xLeft = x - halfW; j0 = Ceiling[xLeft]; u = j0 - xLeft;
                If[u < 0.0, u = 0.0]; If[u > 1.0, u = 1.0];
                If[j0 < 0, j0 += nFFT]; If[j0 >= nFFT, j0 -= nFFT];
                Do[
                    base = s * d1; val = polyCoeffs[[base + 1]];
                    Do[val = val * u + polyCoeffs[[base + q]], {q, 2, d1}];
                    idx = j0 + s + 1; If[idx > nFFT, idx -= nFFT];
                    gridRe[[idx]] += ar * val; gridIm[[idx]] += ai * val,
                    {s, 0, w - 1}
                ]
            ],
            {p, 1, nP}
        ];
        gridRe + I * gridIm
    ]]
];

getLLVMSpreadPolyMulti[] := compileLazyLLVM["SpreadPolyMulti",
    Function[{
        Typed[xVec, "PackedArray"::["Real64", 1]], Typed[aReFlat, "PackedArray"::["Real64", 1]],
        Typed[aImFlat, "PackedArray"::["Real64", 1]], Typed[nFFT, "MachineInteger"],
        Typed[w, "MachineInteger"], Typed[d1, "MachineInteger"],
        Typed[polyCoeffs, "PackedArray"::["Real64", 1]], Typed[nChan, "MachineInteger"]
    },
    Module[{nP = Length[xVec], gRe = Table[0.0, {nChan * nFFT}], gIm = Table[0.0, {nChan * nFFT}],
            wtVec = Table[0.0, {w}], halfW = 0.5 * w,
            p = 1, s = 0, q = 2, ch = 0, base = 0, pBase = 0, gOff = 0,
            x = 0.0, xLeft = 0.0, j0 = 0, u = 0.0, val = 0.0, idx = 1, ar = 0.0, ai = 0.0},
        Do[
            x = xVec[[p]]; xLeft = x - halfW; j0 = Ceiling[xLeft]; u = j0 - xLeft;
            If[u < 0.0, u = 0.0]; If[u > 1.0, u = 1.0];
            If[j0 < 0, j0 += nFFT]; If[j0 >= nFFT, j0 -= nFFT];
            Do[
                base = s * d1; val = polyCoeffs[[base + 1]];
                Do[val = val * u + polyCoeffs[[base + q]], {q, 2, d1}];
                wtVec[[s + 1]] = val,
                {s, 0, w - 1}
            ];
            pBase = (p - 1) * nChan;
            Do[
                ar = aReFlat[[pBase + ch + 1]]; ai = aImFlat[[pBase + ch + 1]];
                If[ar != 0.0 || ai != 0.0,
                    gOff = ch * nFFT;
                    Do[
                        val = wtVec[[s + 1]];
                        idx = j0 + s + 1; If[idx > nFFT, idx -= nFFT];
                        gRe[[gOff + idx]] += ar * val; gIm[[gOff + idx]] += ai * val,
                        {s, 0, w - 1}
                    ]
                ],
                {ch, 0, nChan - 1}
            ],
            {p, 1, nP}
        ];
        gRe + I * gIm
    ]]
];

getLLVMGatherPoly[] := compileLazyLLVM["GatherPoly",
    Function[{
        Typed[yVec, "PackedArray"::["Real64", 1]], Typed[gRe, "PackedArray"::["Real64", 1]],
        Typed[gIm, "PackedArray"::["Real64", 1]], Typed[nFFT, "MachineInteger"],
        Typed[w, "MachineInteger"], Typed[d1, "MachineInteger"], Typed[polyCoeffs, "PackedArray"::["Real64", 1]]
    },
    Module[{nY = Length[yVec], outRe = Table[0.0, {Length[yVec]}], outIm = Table[0.0, {Length[yVec]}],
            halfW = 0.5 * w, i = 1, s = 0, q = 2, base = 0,
            y = 0.0, yLeft = 0.0, j0 = 0, u = 0.0, accRe = 0.0, accIm = 0.0, val = 0.0, idx = 1},
        Do[
            y = yVec[[i]]; yLeft = y - halfW; j0 = Ceiling[yLeft]; u = j0 - yLeft;
            If[u < 0.0, u = 0.0]; If[u > 1.0, u = 1.0];
            If[j0 < 0, j0 += nFFT]; If[j0 >= nFFT, j0 -= nFFT];
            accRe = 0.0; accIm = 0.0;
            Do[
                base = s * d1; val = polyCoeffs[[base + 1]];
                Do[val = val * u + polyCoeffs[[base + q]], {q, 2, d1}];
                idx = j0 + s + 1; If[idx > nFFT, idx -= nFFT];
                accRe += gRe[[idx]] * val; accIm += gIm[[idx]] * val,
                {s, 0, w - 1}
            ];
            outRe[[i]] = accRe; outIm[[i]] = accIm,
            {i, 1, nY}
        ];
        outRe + I * outIm
    ]]
];

getLLVMEvalCheb[] := compileLazyLLVM["EvalCheb",
    Function[{
        Typed[xVec, "PackedArray"::["Real64", 1]], Typed[chebFlat, "PackedArray"::["Real64", 1]],
        Typed[qEff, "MachineInteger"], Typed[r, "MachineInteger"]
    },
    Module[{nP = Length[xVec], out = Table[0.0, {Length[xVec] * r}],
            p = 1, k = 1, j = 1, oBase = 0, cBase = 0,
            z = 0.0, twoZ = 0.0, t0Val = 1.0, t1Val = 0.0, tNew = 0.0},
        Do[
            z = 2.0 * xVec[[p]] - 1.0; twoZ = 2.0 * z; oBase = (p - 1) * r;
            Do[out[[oBase + j]] = chebFlat[[j]], {j, 1, r}];
            If[qEff >= 2,
                t0Val = 1.0; t1Val = z;
                Do[out[[oBase + j]] += t1Val * chebFlat[[r + j]], {j, 1, r}];
                Do[
                    tNew = twoZ * t1Val - t0Val; t0Val = t1Val; t1Val = tNew; cBase = (k - 1) * r;
                    Do[out[[oBase + j]] += tNew * chebFlat[[cBase + j]], {j, 1, r}],
                    {k, 3, qEff}
                ]
            ],
            {p, 1, nP}
        ];
        out
    ]]
];

getLLVMTruncatedDirect[] := compileLazyLLVM["TruncatedDirect",
    Function[{
        Typed[a0Re, "PackedArray"::["Real64", 1]], Typed[a0Im, "PackedArray"::["Real64", 1]],
        Typed[zRe,  "PackedArray"::["Real64", 1]], Typed[zIm,  "PackedArray"::["Real64", 1]],
        Typed[nActVec, "PackedArray"::["MachineInteger", 1]], Typed[nT, "MachineInteger"]
    },
    Module[{nP = Length[a0Re], outRe = Table[0.0, {nT}], outIm = Table[0.0, {nT}],
            p = 1, k = 1, nAct = 1, cRe = 0.0, cIm = 0.0, zr = 0.0, zi = 0.0, nxtRe = 0.0},
        Do[
            cRe = a0Re[[p]]; cIm = a0Im[[p]]; nAct = nActVec[[p]];
            If[nAct > nT, nAct = nT];
            If[nAct >= 1 && (cRe != 0.0 || cIm != 0.0),
                zr = zRe[[p]]; zi = zIm[[p]];
                outRe[[1]] += cRe; outIm[[1]] += cIm;
                Do[
                    nxtRe = cRe * zr - cIm * zi; cIm = cRe * zi + cIm * zr; cRe = nxtRe;
                    outRe[[k]] += cRe; outIm[[k]] += cIm,
                    {k, 2, nAct}
                ]
            ],
            {p, 1, nP}
        ];
        outRe + I * outIm
    ]]
];

getLLVMSpreadES2D[] := compileLazyLLVM["SpreadES2D",
    Function[{
        Typed[x1Vec, "PackedArray"::["Real64", 1]], Typed[x2Vec, "PackedArray"::["Real64", 1]],
        Typed[cRe,   "PackedArray"::["Real64", 1]], Typed[cIm,   "PackedArray"::["Real64", 1]],
        Typed[nFFT1, "MachineInteger"], Typed[w1, "MachineInteger"], Typed[beta1, "Real64"],
        Typed[nFFT2, "MachineInteger"], Typed[w2, "MachineInteger"], Typed[beta2, "Real64"]
    },
    Module[{nP = Length[x1Vec], gRe = Table[0.0, {nFFT1 * nFFT2}], gIm = Table[0.0, {nFFT1 * nFFT2}],
            wt1 = Table[0.0, {w1}], wt2 = Table[0.0, {w2}], idx1 = Table[1, {w1}], idx2 = Table[1, {w2}],
            halfW1 = 0.5 * w1, halfW2 = 0.5 * w2, invHW1 = 2.0 / w1, invHW2 = 2.0 / w2,
            eb1 = Exp[-beta1], eb2 = Exp[-beta2],
            p = 1, s1 = 0, s2 = 0, ar = 0.0, ai = 0.0, x1 = 0.0, x2 = 0.0,
            j01 = 0, j02 = 0, r01 = 0.0, r02 = 0.0, z = 0.0, arg = 0.0,
            rowOff = 0, u1r = 0.0, u1i = 0.0, id = 1, k1 = 1, k2 = 1},
        Do[
            ar = cRe[[p]]; ai = cIm[[p]];
            If[ar != 0.0 || ai != 0.0,
                x1 = x1Vec[[p]]; j01 = Ceiling[x1 - halfW1]; r01 = j01 - x1;
                If[j01 < 0, j01 += nFFT1]; If[j01 >= nFFT1, j01 -= nFFT1];
                Do[
                    z = (r01 + s1) * invHW1; arg = 1.0 - z * z;
                    wt1[[s1 + 1]] = If[arg > 0.0, Exp[beta1 * (Sqrt[arg] - 1.0)], eb1];
                    id = j01 + s1 + 1; If[id > nFFT1, id -= nFFT1]; idx1[[s1 + 1]] = id,
                    {s1, 0, w1 - 1}
                ];
                x2 = x2Vec[[p]]; j02 = Ceiling[x2 - halfW2]; r02 = j02 - x2;
                If[j02 < 0, j02 += nFFT2]; If[j02 >= nFFT2, j02 -= nFFT2];
                Do[
                    z = (r02 + s2) * invHW2; arg = 1.0 - z * z;
                    wt2[[s2 + 1]] = If[arg > 0.0, Exp[beta2 * (Sqrt[arg] - 1.0)], eb2];
                    id = j02 + s2 + 1; If[id > nFFT2, id -= nFFT2]; idx2[[s2 + 1]] = id,
                    {s2, 0, w2 - 1}
                ];
                Do[
                    rowOff = (idx1[[k1]] - 1) * nFFT2; u1r = ar * wt1[[k1]]; u1i = ai * wt1[[k1]];
                    Do[
                        id = rowOff + idx2[[k2]];
                        gRe[[id]] += u1r * wt2[[k2]]; gIm[[id]] += u1i * wt2[[k2]],
                        {k2, 1, w2}
                    ],
                    {k1, 1, w1}
                ]
            ],
            {p, 1, nP}
        ];
        gRe + I * gIm
    ]]
];

(* ====================================================================== *)
(* NUFFT Plan Construction & Quadrature Helpers                           *)
(* ====================================================================== *)
nextSmoothSize[nMin_Integer] := Module[{candidates = {}, n2, n3, n5, lim},
    lim = Max[8, 2 * nMin + 16];
    n2 = 2;
    While[n2 <= lim,
        n3 = n2;
        While[n3 <= lim,
            n5 = n3;
            While[n5 <= lim,
                If[n5 >= nMin && EvenQ[n5], AppendTo[candidates, n5]];
                n5 *= 5
            ];
            n3 *= 3
        ];
        n2 *= 2
    ];
    If[candidates === {}, Max[nMin, 8], Min[candidates]]
];

(* Golub-Welsch Gauss-Legendre nodes/weights on [0, 1] *)
gaussLegendre01Clean[nPts_Integer] := Module[{k, b, mat, evals, evecs},
    k = Range[1, nPts - 1];
    b = k / Sqrt[4. * k^2 - 1.];
    mat = SparseArray[{Band[{1, 2}] -> b, Band[{2, 1}] -> b}, {nPts, nPts}];
    {evals, evecs} = Eigensystem[N[Normal[mat], 25]];
    {Developer`ToPackedArray[N[0.5 * (1. + evals)]], Developer`ToPackedArray[N[evecs[[All, 1]]^2]]}
];

buildESPolyTableLocal[w_Integer, beta_?NumericQ, d_Integer] := Module[
    {nodes, vals, cheb, uSym, poly, cList, tab, s, k, j, norm},
    nodes = Table[0.5 * (1. + Cos[Pi * j / d]), {j, 0, d}];
    uSym = Symbol["Signals`SignalExpansionKernels`Private`uESPoly"];
    tab = Flatten[Table[
        vals = Table[
            With[{z = (2 * (s - w / 2 + u)) / w},
                Exp[N[beta, 25] * (Sqrt[Max[0, 1 - z^2]] - 1)]
            ],
            {u, nodes}
        ];
        cheb = Table[
            norm = If[k == 0 || k == d, 0.5, 1.];
            norm * (2. / d) * Sum[If[j == 0 || j == d, 0.5, 1.] * vals[[j + 1]] * Cos[k * j * Pi / d], {j, 0, d}],
            {k, 0, d}
        ];
        poly = Expand[Sum[cheb[[k + 1]] * ChebyshevT[k, 2 * uSym - 1], {k, 0, d}]];
        cList = PadRight[CoefficientList[poly, uSym], d + 1, 0];
        N[Reverse[cList]],
        {s, 0, w - 1}
    ]];
    Developer`ToPackedArray[tab]
];

getPlan1DLocal[nT_Integer, dt_?NumericQ, tol_?NumericQ, sigma_?NumericQ, wantPoly_, useCache_] := Module[
    {key, w, beta, gamma = 0.975, nFFT, nc, nGL, glZ, glW, mVec, xiVec, phiZ, deconvVec, d, polyCoeffs, res},
    key = {nT, Round[Log10[Max[1.*^-16, Abs[dt]]], 1.*^-6], Round[Log10[tol], 0.05], Round[sigma, 0.01], TrueQ[wantPoly]};
    If[TrueQ[useCache] && KeyExistsQ[planCache1D, key], Return[planCache1D[key]]];
    If[Length[planCache1D] > 128, planCache1D = <||>];

    w    = Max[3, Min[16, Ceiling[Log[1. / tol] / (gamma * Pi * Sqrt[1. - 1. / sigma]) + 1.5]]];
    beta = N[gamma * Pi * w * (1. - 1. / (2. * sigma))];
    nFFT = nextSmoothSize[Max[2 * w + 4, Ceiling[sigma * nT]]];
    nc   = Floor[nT / 2];

    nGL = Max[32, 2 * w + 12];
    {glZ, glW} = gaussLegendre01Clean[nGL];
    mVec  = Developer`ToPackedArray[N[Range[-nc, nT - 1 - nc]]];
    xiVec = Developer`ToPackedArray[(2. * Pi / nFFT) * mVec];
    phiZ  = Exp[beta * (Sqrt[Max[0., 1. - #^2]] - 1.)] & /@ glZ;
    deconvVec = Developer`ToPackedArray[w * (Cos[Outer[Times, 0.5 * w * xiVec, glZ]] . (glW * phiZ))];

    If[TrueQ[wantPoly],
        d = Min[16, w + 4];
        polyCoeffs = buildESPolyTableLocal[w, beta, d],
        d = 0;
        polyCoeffs = Developer`ToPackedArray[{0.}]
    ];

    res = <|"nT" -> nT, "dt" -> N[dt], "nc" -> nc, "tol" -> tol, "sigma" -> sigma,
            "w" -> w, "beta" -> beta, "nFFT" -> nFFT, "d" -> d,
            "deconvVec" -> deconvVec, "polyCoeffs" -> polyCoeffs|>;
    If[TrueQ[useCache], planCache1D[key] = res];
    res
];

getPlan2DLocal[n1_Integer, dt1_?NumericQ, n2_Integer, dt2_?NumericQ, tol_?NumericQ, sigma_?NumericQ, useCache_] := Module[
    {key, res},
    key = {n1, Round[Log10[Max[1.*^-16, Abs[dt1]]], 1.*^-6], n2, Round[Log10[Max[1.*^-16, Abs[dt2]]], 1.*^-6], Round[Log10[tol], 0.05], Round[sigma, 0.01]};
    If[TrueQ[useCache] && KeyExistsQ[planCache2D, key], Return[planCache2D[key]]];
    If[Length[planCache2D] > 32, planCache2D = <||>];
    res = <|"Plan1" -> getPlan1DLocal[n1, dt1, tol, sigma, False, useCache],
            "Plan2" -> getPlan1DLocal[n2, dt2, tol, sigma, False, useCache]|>;
    If[TrueQ[useCache], planCache2D[key] = res];
    res
];

(* ====================================================================== *)
(* Low-Level Evaluators: Direct, NUFFT Type-1/Batch/Type-3, Hetero        *)
(* ====================================================================== *)
evalDirectBLAS3Small[amps_List, omegas_List, lams_List, t0_?NumericQ, dt_?NumericQ, nT_Integer, tOrigin_?NumericQ] := Module[
    {a0, rates, bSize, vMat, nChunks, chunkTimes, aMat},
    If[nT <= 0, Return[{}]];
    a0 = Developer`ToPackedArray[amps * Exp[I * omegas * t0 - lams * (t0 - tOrigin)]];
    If[nT === 1, Return[Developer`ToPackedArray[{Total[a0]}]]];
    rates = Developer`ToPackedArray[I * omegas - lams];
    bSize = Max[1, Min[nT, Round[Sqrt[N[nT]]]]];
    vMat = Exp[Outer[Times, rates, Developer`ToPackedArray[N[dt] * Range[0, bSize - 1]]]];
    If[nT <= bSize, Return[Developer`ToPackedArray[a0 . vMat]]];
    nChunks = Ceiling[nT / bSize];
    chunkTimes = Developer`ToPackedArray[(bSize * N[dt]) * Range[0, nChunks - 1]];
    aMat = a0 * Exp[Outer[Times, rates, chunkTimes]];
    Developer`ToPackedArray[Take[Flatten[Transpose[Transpose[vMat] . aMat]], nT]]
];

evalDirectUniform1D[amps_List, omegas_List, lams_List, t0_?NumericQ, dt_?NumericQ, nT_Integer, tOrigin_?NumericQ,
                    meth_String, stride_Integer, parOpt_, nChunksOpt_, backendOpt_: Automatic] := Module[
    {nP = Length[amps], a0, zStep, a0Re, a0Im, zRe, zIm, fn, nC, cSize, res, useLLVM},
    a0    = Developer`ToPackedArray[amps * Exp[I * omegas * t0 - lams * (t0 - tOrigin)]];
    zStep = Developer`ToPackedArray[Exp[(I * omegas - lams) * dt]];
    a0Re  = Developer`ToPackedArray[Re[a0]];
    a0Im  = Developer`ToPackedArray[Im[a0]];
    zRe   = Developer`ToPackedArray[Re[zStep]];
    zIm   = Developer`ToPackedArray[Im[zStep]];

    If[meth === "DirectCompensated",
        nC = chooseChunkCountLocal[nP, nT, parOpt, nChunksOpt];
        If[nC > 1,
            cSize = Ceiling[nP / nC];
            res = wvmCompensatedKernel[padRect[a0Re, cSize, nC], padRect[a0Im, cSize, nC],
                padRect[omegas, cSize, nC], padRect[lams, cSize, nC], padRect[zRe, cSize, nC], padRect[zIm, cSize, nC],
                ConstantArray[N[dt], nC], ConstantArray[nT, nC], ConstantArray[Max[64, stride], nC]];
            Return[Developer`ToPackedArray[Total[res[[All, 1]]] + I * Total[res[[All, 2]]]]]
        ];
        res = wvmCompensatedKernel[a0Re, a0Im, omegas, lams, zRe, zIm, N[dt], nT, Max[64, stride]];
        Return[Developer`ToPackedArray[res[[1]] + I * res[[2]]]]
    ];

    useLLVM = meth === "DirectLLVM" || (meth =!= "DirectWVM" && nP <= 6000 && shouldUseLLVM["Direct", backendOpt, N[nP] * N[nT], 1.5*^9]);
    If[useLLVM && Head[fn = getLLVMDirect[]] === CompiledCodeFunction,
        Return[Developer`ToPackedArray[fn[a0Re, a0Im, omegas, lams, zRe, zIm, N[dt], nT, Min[nT, Max[64, stride]]]]]
    ];

    nC = chooseChunkCountLocal[nP, nT, parOpt, nChunksOpt];
    If[nC > 1,
        cSize = Ceiling[nP / nC];
        res = wvmTimeBlockedKernel[padRect[a0Re, cSize, nC], padRect[a0Im, cSize, nC],
            padRect[omegas, cSize, nC], padRect[lams, cSize, nC], padRect[zRe, cSize, nC], padRect[zIm, cSize, nC],
            ConstantArray[N[dt], nC], ConstantArray[nT, nC], ConstantArray[Max[64, stride], nC]];
        Developer`ToPackedArray[Total[res[[All, 1]]] + I * Total[res[[All, 2]]]],
        res = wvmTimeBlockedKernel[a0Re, a0Im, omegas, lams, zRe, zIm, N[dt], nT, Max[64, stride]];
        Developer`ToPackedArray[res[[1]] + I * res[[2]]]
    ]
];

evalDirectNonUniform1D[amps_List, omegas_List, lams_List, tVec_List, tOrigin_?NumericQ, parOpt_, nChunksOpt_] := Module[
    {nP = Length[amps], nT = Length[tVec], dtOrigin, aRe, aIm, nC, cSize, res},
    dtOrigin = Developer`ToPackedArray[tVec - tOrigin];
    aRe = Developer`ToPackedArray[Re[amps]];
    aIm = Developer`ToPackedArray[Im[amps]];
    nC = If[parOpt === False, 1,
        If[IntegerQ[nChunksOpt] && nChunksOpt >= 2, Min[nChunksOpt, nP],
            If[nP * nT < 200000, 1, Min[Max[1, $ProcessorCount], Max[2, Quotient[nP, 64]]]]
        ]
    ];
    If[nC > 1,
        cSize = Ceiling[nP / nC];
        res = wvmNonUniformKernel[padRect[aRe, cSize, nC], padRect[aIm, cSize, nC],
            padRect[omegas, cSize, nC], padRect[lams, cSize, nC],
            ConstantArray[tVec, nC], ConstantArray[dtOrigin, nC]];
        Developer`ToPackedArray[Total[res[[All, 1]]] + I * Total[res[[All, 2]]]],
        res = wvmNonUniformKernel[aRe, aIm, omegas, lams, tVec, dtOrigin];
        Developer`ToPackedArray[res[[1]] + I * res[[2]]]
    ]
];

evalNUFFTType1Common[amps_List, omegas_List, lam_?NumericQ, tVec_List, t0_?NumericQ, dt_?NumericQ,
                     tOrigin_?NumericQ, plan_Association, backendOpt_: Automatic, parOpt_: True, nChunksOpt_: Automatic] := Module[
    {nP = Length[amps], nT = Length[tVec], nc = plan["nc"], nFFT = plan["nFFT"], tMid, cVec, cRe, cIm, xVec,
     fn, uGrid, fAll, fCenter, sig, env},
    If[nT === 1, Return[Developer`ToPackedArray[{Total[amps * Exp[I * omegas * t0 - lam * (t0 - tOrigin)]]}]]];
    tMid = t0 + nc * dt;
    cVec = amps * Exp[I * (omegas * tMid)];
    cRe  = Developer`ToPackedArray[Re[cVec]];
    cIm  = Developer`ToPackedArray[Im[cVec]];
    xVec = Developer`ToPackedArray[(nFFT / (2. * Pi)) * Mod[omegas * dt, 2. * Pi]];

    uGrid = If[shouldUseLLVM["SpreadES", backendOpt, N[nP], 300000.] && Head[fn = getLLVMSpreadES[]] === CompiledCodeFunction,
        Developer`ToPackedArray[fn[xVec, cRe, cIm, nFFT, plan["w"], plan["beta"]]],
        spreadESWVMLocal[xVec, cRe, cIm, nFFT, plan["w"], plan["beta"], parOpt, nChunksOpt]
    ];
    fAll = Fourier[uGrid, FourierParameters -> {1, 1}];
    fCenter = If[nc > 0, Join[fAll[[nFFT - nc + 1 ;; nFFT]], fAll[[1 ;; nT - nc]]], fAll[[1 ;; nT]]];
    sig = fCenter / plan["deconvVec"];
    If[lam != 0.,
        env = Developer`ToPackedArray[Exp[-lam * (tVec - tOrigin)]];
        sig = sig * env
    ];
    Developer`ToPackedArray[sig]
];

evalNUFFTType3Common[amps_List, omegas_List, lam_?NumericQ, tVec_List, tOrigin_?NumericQ,
                     tol_?NumericQ, sigma_?NumericQ, backendOpt_: Automatic] := Module[
    {nP = Length[amps], nT = Length[tVec], wMin, wMax, wMid, W, tMin, tMax, tMid, T,
     gamma = 0.975, w, beta, d, polyCoeffs, nMinReq, nFFT, xShift, yShift, hW, hT,
     nGL, glZ, glW, phiZ, xi1, d1Vec, cVec, cRe, cIm, xGrid, fnS, uGrid, res,
     fAll, kIdx, xi2, d2Vec, gShifted, gRe, gIm, gPad, fnG, gatherRes},
    {wMin, wMax} = MinMax[omegas]; {tMin, tMax} = MinMax[tVec];
    wMid = 0.5 * (wMin + wMax); W = 0.5 * (wMax - wMin);
    tMid = 0.5 * (tMin + tMax); T = 0.5 * (tMax - tMin);
    If[W == 0. || T == 0.,
        Return[Developer`ToPackedArray[Total[amps * Exp[I * omegas * tMid - lam * (tVec - tOrigin)]] * Exp[I * wMid * (tVec - tMid)]]]
    ];
    w    = Max[3, Min[16, Ceiling[Log[1. / tol] / (gamma * Pi * Sqrt[1. - 1. / sigma]) + 1.5]]];
    beta = N[gamma * Pi * w * (1. - 1. / (2. * sigma))];
    d    = Min[16, w + 4];
    polyCoeffs = buildESPolyTableLocal[w, beta, d];
    nMinReq = Max[4 * w + 8, Ceiling[4. * sigma^2 * W * T / Pi + 4. * w]];
    nFFT = nextSmoothSize[nMinReq];

    xShift = omegas - wMid; yShift = tVec - tMid;
    hW = Sqrt[(2. * Pi * W) / (nFFT * T)]; hT = (2. * Pi) / (nFFT * hW);

    nGL = Max[32, 2 * w + 12];
    {glZ, glW} = gaussLegendre01Clean[nGL];
    phiZ = Exp[beta * (Sqrt[Max[0., 1. - #^2]] - 1.)] & /@ glZ;

    xi1   = Developer`ToPackedArray[xShift * hT];
    d1Vec = Developer`ToPackedArray[w * (Cos[Outer[Times, 0.5 * w * xi1, glZ]] . (glW * phiZ))];
    cVec  = (amps * Exp[I * (omegas * tMid)]) / d1Vec;
    cRe   = Developer`ToPackedArray[Re[cVec]];
    cIm   = Developer`ToPackedArray[Im[cVec]];
    xGrid = Developer`ToPackedArray[Mod[xShift / hW, nFFT]];

    uGrid = If[shouldUseLLVM["SpreadPoly", backendOpt, N[nP], 300000.] && Head[fnS = getLLVMSpreadPoly[]] === CompiledCodeFunction,
        Developer`ToPackedArray[fnS[xGrid, cRe, cIm, nFFT, w, d + 1, polyCoeffs]],
        res = wvmSpreadPolyKernel[xGrid, cRe, cIm, nFFT, w, d + 1, polyCoeffs];
        Developer`ToPackedArray[res[[1]] + I * res[[2]]]
    ];
    fAll = Fourier[uGrid, FourierParameters -> {1, 1}];
    kIdx = Developer`ToPackedArray[N[Join[Range[0, nFFT / 2 - 1], Range[-nFFT / 2, -1]]]];
    xi2  = Developer`ToPackedArray[(2. * Pi / nFFT) * kIdx];
    d2Vec = Developer`ToPackedArray[w * (Cos[Outer[Times, 0.5 * w * xi2, glZ]] . (glW * phiZ))];
    With[{dMax = Max[Abs[d2Vec]]},
        d2Vec = Developer`ToPackedArray[Map[If[Abs[#] > 1.*^-6 * dMax, 1. / #, 0.] &, d2Vec]]
    ];
    gShifted = fAll * d2Vec;
    gRe = Developer`ToPackedArray[Re[gShifted]]; gIm = Developer`ToPackedArray[Im[gShifted]];
    gPad = Developer`ToPackedArray[Mod[yShift / hT, nFFT]];
    gatherRes = If[shouldUseLLVM["GatherPoly", backendOpt, N[nT], 300000.] && Head[fnG = getLLVMGatherPoly[]] === CompiledCodeFunction,
        Developer`ToPackedArray[fnG[gPad, gRe, gIm, nFFT, w, d + 1, polyCoeffs]],
        res = wvmGatherPolyKernel[gPad, gRe, gIm, nFFT, w, d + 1, polyCoeffs];
        Developer`ToPackedArray[res[[1]] + I * res[[2]]]
    ];
    Developer`ToPackedArray[gatherRes * Exp[I * (wMid * yShift) - lam * (tVec - tOrigin)]]
];

buildSeparableKernelLocal[dParam_?NumericQ, tol_?NumericQ] := Module[
    {dVal = Max[0., N[dParam]], tVal = Max[1.*^-14, Min[0.1, N[tol]]], cacheKey,
     q, nodes, dctMat, kMat, uMat, sMat, vMat, sVals, sMax, targetErr,
     chkNodes, chkExact, rTry, uScaled, vScaled, cU, cV, uChk, vChk, maxErr,
     rKeep, qEffU, qEffV, qEff, cutTol, res},
    If[dVal <= 1.*^-13,
        Return[<|"D" -> 0., "Rank" -> 1, "ChebQ" -> 1, "QEff" -> 1,
                 "ChebU" -> Developer`ToPackedArray[{{1.}}],
                 "ChebV" -> Developer`ToPackedArray[{{1.}}], "MaxErr" -> 0.|>]
    ];
    cacheKey = {Round[dVal, 1.*^-9], Round[Log10[tVal], 0.05]};
    If[KeyExistsQ[separableCache, cacheKey], Return[separableCache[cacheKey]]];
    If[Length[separableCache] > 64, separableCache = <||>];

    q = Max[24, Min[160, Ceiling[5.5 * Sqrt[dVal] + 2.5 * Log[1. / tVal]]]];
    nodes = Developer`ToPackedArray[Table[0.5 * (1. + Cos[Pi * j / (q - 1)]), {j, 0, q - 1}]];
    dctMat = Developer`ToPackedArray[Table[
        (If[k == 0 || k == q - 1, 0.5, 1.] * (2. / (q - 1))) *
        If[j == 0 || j == q - 1, 0.5, 1.] * Cos[k * j * Pi / (q - 1)],
        {k, 0, q - 1}, {j, 0, q - 1}
    ]];
    kMat = Developer`ToPackedArray[Exp[-dVal * Outer[Times, nodes, nodes]]];
    {uMat, sMat, vMat} = SingularValueDecomposition[kMat];
    sVals = Diagonal[sMat]; sMax = First[sVals];
    targetErr = 0.25 * tVal;
    chkNodes = Developer`ToPackedArray[Range[0.005, 0.995, 1. / 63.]];
    chkExact = Developer`ToPackedArray[Exp[-dVal * Outer[Times, chkNodes, chkNodes]]];

    rKeep = q;
    Do[
        If[rTry < q && sVals[[rTry + 1]] / sMax > 0.5 * tVal, Continue[]];
        uScaled = uMat[[All, 1 ;; rTry]] . DiagonalMatrix[Sqrt[sVals[[1 ;; rTry]]]];
        vScaled = vMat[[All, 1 ;; rTry]] . DiagonalMatrix[Sqrt[sVals[[1 ;; rTry]]]];
        cU = dctMat . uScaled; cV = dctMat . vScaled;
        uChk = wvmEvalChebBasis[chkNodes, cU]; vChk = wvmEvalChebBasis[chkNodes, cV];
        maxErr = Max[Abs[chkExact - uChk . Transpose[vChk]]];
        If[maxErr <= targetErr || rTry === q, rKeep = rTry; Break[]],
        {rTry, 1, q}
    ];
    cU = dctMat . (uMat[[All, 1 ;; rKeep]] . DiagonalMatrix[Sqrt[sVals[[1 ;; rKeep]]]]);
    cV = dctMat . (vMat[[All, 1 ;; rKeep]] . DiagonalMatrix[Sqrt[sVals[[1 ;; rKeep]]]]);
    cutTol = 0.01 * targetErr / q;
    qEffU = Max[1, Max[Join[{1}, Flatten[Position[UnitStep[Max /@ Abs[cU] - cutTol], 1, {1}, Heads -> False]]]]];
    qEffV = Max[1, Max[Join[{1}, Flatten[Position[UnitStep[Max /@ Abs[cV] - cutTol], 1, {1}, Heads -> False]]]]];
    qEff  = Min[q, Max[qEffU, qEffV] + 1];

    res = <|"D" -> dVal, "Rank" -> rKeep, "ChebQ" -> q, "QEff" -> qEff,
            "ChebU" -> Developer`ToPackedArray[cU[[1 ;; qEff, All]]],
            "ChebV" -> Developer`ToPackedArray[cV[[1 ;; qEff, All]]], "MaxErr" -> maxErr|>;
    separableCache[cacheKey] = res;
    res
];

evalChebLocal[xVec_List, chebMat_List, qEff_Integer, r_Integer, backendOpt_: Automatic] := Module[{fn, res},
    If[shouldUseLLVM["EvalCheb", backendOpt, N[Length[xVec]] * N[r], 600000.] && Head[fn = getLLVMEvalCheb[]] === CompiledCodeFunction,
        res = Developer`ToPackedArray[fn[xVec, Developer`ToPackedArray[Flatten[chebMat]], qEff, r]];
        Return[Developer`ToPackedArray[ArrayReshape[res, {Length[xVec], r}]]]
    ];
    Developer`ToPackedArray[wvmEvalChebBasis[xVec, chebMat]]
];

evalHeteroExactGroupsLocal[amps_List, omegas_List, lams_List, tVec_List, t0_?NumericQ, dt_?NumericQ,
                           tOrigin_?NumericQ, tol_?NumericQ, sigma_?NumericQ, useCache_, backendOpt_: Automatic, parOpt_: True, nChunksOpt_: Automatic] := Module[
    {nT = Length[tVec], ord, lS, wS, aS, runs, starts, ends, reps, nGrp, plan, acc, g, s, e, mG},
    ord = Ordering[lams];
    lS = lams[[ord]]; wS = omegas[[ord]]; aS = amps[[ord]];
    runs = sortedLambdaRunsWVM[lS, 0.];
    starts = Round[runs[[1]]]; ends = Round[runs[[2]]]; reps = runs[[3]];
    nGrp = Length[starts];
    plan = getPlan1DLocal[nT, dt, tol, sigma, False, useCache];
    acc  = Developer`ToPackedArray[ConstantArray[0. + 0. I, nT]];
    Do[
        s = starts[[g]]; e = ends[[g]]; mG = e - s + 1;
        If[mG <= 300 || mG * nT <= 250000,
            acc += evalDirectUniform1D[aS[[s ;; e]], wS[[s ;; e]], ConstantArray[reps[[g]], mG],
                t0, dt, nT, tOrigin, "Auto", 512, parOpt, nChunksOpt, backendOpt],
            acc += evalNUFFTType1Common[aS[[s ;; e]], wS[[s ;; e]], reps[[g]], tVec, t0, dt, tOrigin, plan, backendOpt, parOpt, nChunksOpt]
        ],
        {g, nGrp}
    ];
    Developer`ToPackedArray[acc]
];

evalHeteroSeparableLocal[amps_List, omegas_List, lams_List, tVec_List, t0_?NumericQ, dt_?NumericQ,
                         tOrigin_?NumericQ, tol_?NumericQ, sigma_?NumericQ, useCache_, backendOpt_: Automatic, parOpt_: True, nChunksOpt_: Automatic] := Module[
    {nT = Length[tVec], tauVec, tauLo, tauHi, dTau, lamMin, lamMax, dLam, dParam,
     sk, r, qEff, aTilde, xNu, yTau, uMat, vMat, plan, nc, nFFT, tMid, phaseVec, cMat,
     aReFlat, aImFlat, xVec, fn, raw, res, grids, fRows, idxSlice, chanT, sumVec, lamEnv},
    tauVec = Developer`ToPackedArray[tVec - tOrigin];
    {tauLo, tauHi} = MinMax[tauVec]; dTau = tauHi - tauLo;
    {lamMin, lamMax} = MinMax[lams]; dLam = lamMax - lamMin;
    dParam = dLam * dTau;
    plan = getPlan1DLocal[nT, dt, tol, sigma, True, useCache];
    If[dParam <= 1.*^-13, Return[evalNUFFTType1Common[amps, omegas, lamMin, tVec, t0, dt, tOrigin, plan, backendOpt, parOpt, nChunksOpt]]];

    sk = buildSeparableKernelLocal[dParam, tol];
    r  = sk["Rank"]; qEff = sk["QEff"];
    aTilde = If[tauLo == 0., amps, Developer`ToPackedArray[amps * Exp[-lams * tauLo]]];
    xNu  = Developer`ToPackedArray[(lams - lamMin) / dLam];
    yTau = Developer`ToPackedArray[(tauVec - tauLo) / dTau];
    uMat = evalChebLocal[xNu, sk["ChebU"], qEff, r, backendOpt];
    vMat = evalChebLocal[yTau, sk["ChebV"], qEff, r, backendOpt];

    nc = plan["nc"]; nFFT = plan["nFFT"]; tMid = t0 + nc * dt;
    phaseVec = Developer`ToPackedArray[Exp[I * (omegas * tMid)]];
    cMat = Developer`ToPackedArray[Transpose[uMat] * ConstantArray[aTilde * phaseVec, r]];
    aReFlat = Developer`ToPackedArray[Flatten[Re[Transpose[cMat]]]];
    aImFlat = Developer`ToPackedArray[Flatten[Im[Transpose[cMat]]]];
    xVec = Developer`ToPackedArray[(nFFT / (2. * Pi)) * Mod[omegas * dt, 2. * Pi]];

    grids = If[shouldUseLLVM["SpreadPolyMulti", backendOpt, N[Length[amps]] * N[r], 800000.] && Head[fn = getLLVMSpreadPolyMulti[]] === CompiledCodeFunction,
        raw = Developer`ToPackedArray[fn[xVec, aReFlat, aImFlat, nFFT, plan["w"], plan["d"] + 1, plan["polyCoeffs"], r]];
        Developer`ToPackedArray[ArrayReshape[raw, {r, nFFT}]],
        res = wvmSpreadPolyMultiKernel[xVec, aReFlat, aImFlat, nFFT, plan["w"], plan["d"] + 1, plan["polyCoeffs"], r];
        Developer`ToPackedArray[ArrayReshape[res[[1]] + I * res[[2]], {r, nFFT}]]
    ];
    fRows = Fourier[#, FourierParameters -> {1, 1}] & /@ grids;
    idxSlice = Join[If[nc > 0, Range[nFFT - nc + 1, nFFT], {}], Range[1, nT - nc]];
    chanT  = (#[[idxSlice]] / plan["deconvVec"]) & /@ fRows;
    sumVec = Total[Transpose[vMat] * chanT];
    lamEnv = If[lamMin != 0., Developer`ToPackedArray[Exp[-lamMin * (tauVec - tauLo)]], 1.];
    Developer`ToPackedArray[sumVec * lamEnv]
];

evalHeteroDyadicHybridLocal[amps_List, omegas_List, lams_List, tVec_List, t0_?NumericQ, dt_?NumericQ,
                            tOrigin_?NumericQ, tol_?NumericQ, sigma_?NumericQ, dSlow_?NumericQ, useCache_, backendOpt_: Automatic, parOpt_: True, nChunksOpt_: Automatic] := Module[
    {nT = Length[tVec], tauVec, tauLo, tauHi, dTau, lamMin, dProd, slowPos, fastPos, sig,
     fA, fW, fL, a0, zStep, maxA0, tolCut, nActVec, fn, res},
    tauVec = Developer`ToPackedArray[tVec - tOrigin];
    {tauLo, tauHi} = MinMax[tauVec]; dTau = tauHi - tauLo;
    lamMin = Min[lams]; dProd = (lams - lamMin) * dTau;
    slowPos = Flatten[Position[UnitStep[dSlow - dProd], 1, {1}, Heads -> False]];
    fastPos = Flatten[Position[UnitStep[dProd - dSlow] * (1 - UnitStep[dSlow - dProd]), 1, {1}, Heads -> False]];

    sig = Developer`ToPackedArray[ConstantArray[0. + 0. I, nT]];
    If[slowPos =!= {},
        sig += evalHeteroSeparableLocal[amps[[slowPos]], omegas[[slowPos]], lams[[slowPos]],
            tVec, t0, dt, tOrigin, tol, sigma, useCache, backendOpt, parOpt, nChunksOpt]
    ];
    If[fastPos =!= {},
        fA = amps[[fastPos]]; fW = omegas[[fastPos]]; fL = lams[[fastPos]];
        a0 = Developer`ToPackedArray[fA * Exp[I * fW * t0 - fL * (t0 - tOrigin)]];
        zStep = Developer`ToPackedArray[Exp[(I * fW - fL) * dt]];
        maxA0 = Max[Abs[a0]];
        If[maxA0 > 0.,
            tolCut = 0.05 * tol * maxA0 / Sqrt[N[Max[1, Length[fA]]]];
            nActVec = Developer`ToPackedArray[MapThread[
                Function[{aAbs, lam},
                    If[lam * dt <= 1.*^-12, nT, Min[nT, Max[1, Ceiling[Log[Max[1., aAbs / tolCut]] / (lam * dt)] + 1]]]
                ],
                {Abs[a0], fL}
            ]];
            If[shouldUseLLVM["TruncatedDirect", backendOpt, N[Total[nActVec]], 3000000.] && Head[fn = getLLVMTruncatedDirect[]] === CompiledCodeFunction,
                sig += Developer`ToPackedArray[fn[Re[a0], Im[a0], Re[zStep], Im[zStep], nActVec, nT]],
                res = wvmTruncatedDirectKernel[Re[a0], Im[a0], Re[zStep], Im[zStep], nActVec, nT];
                sig += Developer`ToPackedArray[res[[1]] + I * res[[2]]]
            ]
        ]
    ];
    Developer`ToPackedArray[sig]
];

(* ====================================================================== *)
(* Public 1D Fast Expanders (Uniform & General Times)                     *)
(* ====================================================================== *)
ExpandLorentzianUniformFast[peaks_List, tmin_?NumericQ, dtIn_?NumericQ, nT_Integer, t0In_: Automatic, opts___Rule] := Module[
    {t0 = N[If[t0In === Automatic, tmin, t0In]], tOrigin, dt = N[dtIn], tVec,
     methOpt, tol, sigma, stride, dropTol, useCache, parOpt, nChunksOpt, dSlow, backendOpt,
     amps, omegas, lams, nP, lamMin, lamMax, dTau, dParam, directCostS, nufftCostS,
     nUniq, rEst, sepCostS, nFast, slowFrac, avgActFrac, dyadCostS, plan},
    If[nT <= 0, Return[{}]];
    If[peaks === {}, Return[Developer`ToPackedArray[ConstantArray[0. + 0. I, nT]]]];

    methOpt    = ExpansionMethod /. {opts} /. kernelDefaults;
    tol        = Max[1.*^-14, Min[0.1, N[ExpansionTolerance /. {opts} /. kernelDefaults]]];
    sigma      = Max[1.1, N[OversamplingFactor /. {opts} /. kernelDefaults]];
    stride     = Max[64, Round[ReinitStride /. {opts} /. kernelDefaults]];
    dropTol    = N[DropTolerance /. {opts} /. kernelDefaults];
    useCache   = PlanCache /. {opts} /. kernelDefaults;
    parOpt     = ParallelExecution /. {opts} /. kernelDefaults;
    nChunksOpt = NumChunks /. {opts} /. kernelDefaults;
    dSlow      = N[SlowThresholdD /. {opts} /. kernelDefaults];
    backendOpt = CompilationBackend /. {opts} /. kernelDefaults;
    tOrigin    = TimingOrigin /. {opts} /. kernelDefaults;
    If[tOrigin === Automatic, tOrigin = N[tmin], tOrigin = N[tOrigin]];

    {amps, omegas, lams} = normalizePeaks1D[peaks, dropTol, tOrigin];
    If[amps === {}, Return[Developer`ToPackedArray[ConstantArray[0. + 0. I, nT]]]];
    nP = Length[amps];
    tVec = Developer`ToPackedArray[N[t0 + dt * Range[0, nT - 1]]];
    If[nT === 1, Return[Developer`ToPackedArray[{Total[amps * Exp[I * omegas * t0 - lams * (t0 - tOrigin)]]}]]];

    {lamMin, lamMax} = MinMax[lams];
    dTau   = Max[Abs[tVec - tOrigin]];
    dParam = (lamMax - lamMin) * dTau;

    Switch[methOpt,
        "DirectLLVM" | "DirectWVM" | "DirectCompensated",
            Return[evalDirectUniform1D[amps, omegas, lams, t0, dt, nT, tOrigin, methOpt, stride, parOpt, nChunksOpt, backendOpt]],
        "DirectNonUniform",
            Return[evalDirectNonUniform1D[amps, omegas, lams, tVec, tOrigin, parOpt, nChunksOpt]],
        "NUFFT" | "NUFFTType1",
            If[dParam <= 1.*^-12 * Max[1., Abs[lamMin]],
                plan = getPlan1DLocal[nT, dt, tol, sigma, False, useCache];
                Return[evalNUFFTType1Common[amps, omegas, lamMin, tVec, t0, dt, tOrigin, plan, backendOpt, parOpt, nChunksOpt]],
                Return[evalHeteroSeparableLocal[amps, omegas, lams, tVec, t0, dt, tOrigin, tol, sigma, useCache, backendOpt, parOpt, nChunksOpt]]
            ],
        "HeteroExactGroups",
            Return[evalHeteroExactGroupsLocal[amps, omegas, lams, tVec, t0, dt, tOrigin, tol, sigma, useCache, backendOpt, parOpt, nChunksOpt]],
        "HeteroSeparable",
            Return[evalHeteroSeparableLocal[amps, omegas, lams, tVec, t0, dt, tOrigin, tol, sigma, useCache, backendOpt, parOpt, nChunksOpt]],
        "HeteroDyadicHybrid",
            Return[evalHeteroDyadicHybridLocal[amps, omegas, lams, tVec, t0, dt, tOrigin, tol, sigma, dSlow, useCache, backendOpt, parOpt, nChunksOpt]]
    ];

    (* Small-workload bypass to packed baby-step/giant-step BLAS-3 (~1e-15 accuracy, sub-ms latency) *)
    If[MatchQ[methOpt, Automatic | "Fast"] && ((nP <= 64 && nT <= 131072) || nP * nT <= 65536),
        Return[evalDirectBLAS3Small[amps, omegas, lams, t0, dt, nT, tOrigin]]
    ];

    (* Automatic / "Fast" / "Streaming" empirically calibrated dispatch *)
    If[dParam <= 1.*^-12 * Max[1., Abs[lamMin]],
        If[nP <= 300 || nT <= 64 || nP * nT <= 250000,
            Return[evalDirectUniform1D[amps, omegas, lams, t0, dt, nT, tOrigin,
                "Auto", stride, parOpt, nChunksOpt, backendOpt]],
            plan = getPlan1DLocal[nT, dt, tol, sigma, False, useCache];
            Return[evalNUFFTType1Common[amps, omegas, lamMin, tVec, t0, dt, tOrigin, plan, backendOpt, parOpt, nChunksOpt]]
        ]
    ];

    If[nP <= 1500 || nT <= 256 || nP * nT <= 2.5*^6,
        Return[evalDirectUniform1D[amps, omegas, lams, t0, dt, nT, tOrigin,
            "Auto", stride, parOpt, nChunksOpt, backendOpt]]
    ];

    directCostS = 0.002 + 2.15*^-9 * N[nP] * N[nT];
    nUniq = wvmCountUniqueUpTo[lams, 1.*^-12 * Max[1., Abs[lamMin]], 16];
    If[nUniq <= 16,
        With[{grpCost = 0.004 * nUniq + 6.5*^-7 * N[nP] + 1.5*^-6 * nUniq * N[nT]},
            If[grpCost < 0.85 * directCostS,
                Return[evalHeteroExactGroupsLocal[amps, omegas, lams, tVec, t0, dt, tOrigin, tol, sigma, useCache, backendOpt, parOpt, nChunksOpt]]
            ]
        ]
    ];

    rEst = Max[2., 0.30 * Log[1. / tol] + 1.05 * Sqrt[dParam] + 1.5];
    sepCostS = 0.045 + 3.0*^-7 * N[nP] * rEst + 2.5*^-7 * rEst * N[nT];
    If[dParam > dSlow && nP >= 20000,
        nFast = Total[UnitStep[(lams - lamMin) * dTau - dSlow]];
        If[nFast > 0,
            slowFrac = 1. - N[nFast] / N[nP];
            avgActFrac = Min[1., Log[10. / tol] / Max[1., 0.5 * dParam]];
            dyadCostS = 0.04 + 3.0*^-7 * (N[nP] * slowFrac) * Min[rEst, 14.] + 2.5*^-9 * N[nFast] * N[nT] * avgActFrac;
            If[dyadCostS < 0.90 * Min[sepCostS, directCostS] && dParam >= 150.,
                Return[evalHeteroDyadicHybridLocal[amps, omegas, lams, tVec, t0, dt, tOrigin, tol, sigma, dSlow, useCache, backendOpt, parOpt, nChunksOpt]]
            ]
        ]
    ];

    If[sepCostS < 0.85 * directCostS,
        evalHeteroSeparableLocal[amps, omegas, lams, tVec, t0, dt, tOrigin, tol, sigma, useCache, backendOpt, parOpt, nChunksOpt],
        evalDirectUniform1D[amps, omegas, lams, t0, dt, nT, tOrigin, "DirectWVM", stride, parOpt, nChunksOpt, backendOpt]
    ]
];

ExpandLorentzianTimesFast[peaks_List, tvals_List, tmin_?NumericQ, opts___Rule] := Module[
    {tVec = Developer`ToPackedArray[N[tvals]], nT, isUnif, t0, dt, methOpt, tol, sigma,
     dropTol, parOpt, nChunksOpt, backendOpt, amps, omegas, lams, nP, lamMin, lamMax, dParam, wMin, wMax, wtProd},
    nT = Length[tVec];
    If[nT === 0, Return[{}]];
    If[peaks === {}, Return[Developer`ToPackedArray[ConstantArray[0. + 0. I, nT]]]];

    {isUnif, t0, dt} = detectUniformGridLocal[tVec];
    If[isUnif,
        Return[ExpandLorentzianUniformFast[peaks, N[tmin], dt, nT, t0, TimingOrigin -> N[tmin], opts]]
    ];

    methOpt    = ExpansionMethod /. {opts} /. kernelDefaults;
    tol        = Max[1.*^-14, Min[0.1, N[ExpansionTolerance /. {opts} /. kernelDefaults]]];
    sigma      = Max[1.25, N[OversamplingFactor /. {opts} /. kernelDefaults]];
    dropTol    = N[DropTolerance /. {opts} /. kernelDefaults];
    parOpt     = ParallelExecution /. {opts} /. kernelDefaults;
    nChunksOpt = NumChunks /. {opts} /. kernelDefaults;
    backendOpt = CompilationBackend /. {opts} /. kernelDefaults;

    {amps, omegas, lams} = normalizePeaks1D[peaks, dropTol, N[tmin]];
    If[amps === {}, Return[Developer`ToPackedArray[ConstantArray[0. + 0. I, nT]]]];
    nP = Length[amps];
    {lamMin, lamMax} = MinMax[lams];
    dParam = (lamMax - lamMin) * Max[Abs[tVec - N[tmin]]];

    If[methOpt === "DirectNonUniform",
        Return[evalDirectNonUniform1D[amps, omegas, lams, tVec, N[tmin], parOpt, nChunksOpt]]
    ];
    If[methOpt === "NUFFTType3" || methOpt === "NUFFT",
        If[dParam <= 1.*^-12 * Max[1., Abs[lamMin]],
            Return[evalNUFFTType3Common[amps, omegas, lamMin, tVec, N[tmin], tol, sigma, backendOpt]],
            Return[evalDirectNonUniform1D[amps, omegas, lams, tVec, N[tmin], parOpt, nChunksOpt]]
        ]
    ];

    If[dParam <= 1.*^-12 * Max[1., Abs[lamMin]] && nP >= 4000 && nT >= 512,
        {wMin, wMax} = MinMax[omegas];
        wtProd = 0.25 * (wMax - wMin) * (Max[tVec] - Min[tVec]);
        If[(4. * sigma^2 * wtProd / Pi) < 0.08 * nP * nT && (4. * sigma^2 * wtProd / Pi) <= 2^19,
            Return[evalNUFFTType3Common[amps, omegas, lamMin, tVec, N[tmin], tol, sigma, backendOpt]]
        ]
    ];
    evalDirectNonUniform1D[amps, omegas, lams, tVec, N[tmin], parOpt, nChunksOpt]
];

(* ====================================================================== *)
(* Public 2D Fast Expander & Bounded-Memory Streaming Accumulator         *)
(* ====================================================================== *)
ExpandLorentzian2DFast[peaks2D_List, t1Grid_List, t2Grid_List, opts___Rule] := Module[
    {t1 = Developer`ToPackedArray[N[t1Grid]], t2 = Developer`ToPackedArray[N[t2Grid]],
     n1, n2, methOpt, tol, sigma, dropTol, maxMB, useCache, backendOpt, parOpt, amps, w1, l1, w2, l2, nP,
     t1Orig, t2Orig, u1, t10, dt1, u2, t20, dt2, d1, d2, l1Min, l2Min,
     eval2DDir, eval2DSemi, plan2, p1, p2, nc1, nFFT1, tMid1, nc2, nFFT2, tMid2,
     cVec, cRe, cIm, x1Vec, x2Vec, fn, uFlat, uMat, fMat, i1, i2, fCenter, sig, env1, env2},
    n1 = Length[t1]; n2 = Length[t2];
    If[n1 === 0 || n2 === 0, Return[Developer`ToPackedArray[ConstantArray[0. + 0. I, {n1, n2}]]]];
    methOpt    = ExpansionMethod /. {opts} /. kernelDefaults;
    tol        = Max[1.*^-14, Min[0.1, N[ExpansionTolerance /. {opts} /. kernelDefaults]]];
    sigma      = Max[1.1, N[OversamplingFactor /. {opts} /. kernelDefaults]];
    dropTol    = N[DropTolerance /. {opts} /. kernelDefaults];
    maxMB      = N[MaxBufferMB /. {opts} /. kernelDefaults];
    useCache   = PlanCache /. {opts} /. kernelDefaults;
    backendOpt = CompilationBackend /. {opts} /. kernelDefaults;
    parOpt     = ParallelExecution /. {opts} /. kernelDefaults;

    {amps, w1, l1, w2, l2} = parsePeaks2DLocal[peaks2D, dropTol];
    If[amps === {}, Return[Developer`ToPackedArray[ConstantArray[0. + 0. I, {n1, n2}]]]];
    nP = Length[amps];
    t1Orig = t1[[1]]; t2Orig = t2[[1]];

    eval2DDir[] := Module[{dt1V = Developer`ToPackedArray[t1 - t1Orig], dt2V = Developer`ToPackedArray[t2 - t2Orig],
                           tileSz, acc, s, e, e1, e2},
        tileSz = Max[32, Min[nP, Floor[(maxMB * 1024.^2) / (16. * (n1 + n2 + 4.))] ]];
        acc = Developer`ToPackedArray[ConstantArray[0. + 0. I, {n1, n2}]];
        For[s = 1, s <= nP, s += tileSz,
            e = Min[s + tileSz - 1, nP];
            e1 = Developer`ToPackedArray[amps[[s ;; e]] * Exp[I * Outer[Times, w1[[s ;; e]], t1] - Outer[Times, l1[[s ;; e]], dt1V]]];
            e2 = Developer`ToPackedArray[Exp[I * Outer[Times, w2[[s ;; e]], t2] - Outer[Times, l2[[s ;; e]], dt2V]]];
            acc += Transpose[e1] . e2
        ];
        Developer`ToPackedArray[acc]
    ];

    {u1, t10, dt1} = detectUniformGridLocal[t1];
    {u2, t20, dt2} = detectUniformGridLocal[t2];
    {l1Min, d1} = {Min[l1], (Max[l1] - Min[l1]) * Max[1., Abs[t1[[-1]] - t1Orig]]};
    {l2Min, d2} = {Min[l2], (Max[l2] - Min[l2]) * Max[1., Abs[t2[[-1]] - t2Orig]]};

    eval2DSemi[] := Module[{dt1V = Developer`ToPackedArray[t1 - t1Orig], e1Mat, ampsMat, planB, ncB, nFFTB, tMidB,
                            phaseVec, cMat, aReFlat, aImFlat, xVec, fnM, raw, grids, fRows, idxSlice, sigMat, env,
                            rowBlk = Max[1, Quotient[4000000, nP]], s, e, resM},
        e1Mat   = Developer`ToPackedArray[Exp[I * Outer[Times, t1, w1] - Outer[Times, dt1V, l1]]];
        ampsMat = Developer`ToPackedArray[Transpose[amps * Transpose[e1Mat]]];
        planB   = getPlan1DLocal[n2, dt2, tol, sigma, True, useCache];
        ncB = planB["nc"]; nFFTB = planB["nFFT"]; tMidB = t20 + ncB * dt2;
        phaseVec = Developer`ToPackedArray[Exp[I * (w2 * tMidB)]];
        xVec     = Developer`ToPackedArray[(nFFTB / (2. * Pi)) * Mod[w2 * dt2, 2. * Pi]];
        idxSlice = Join[If[ncB > 0, Range[nFFTB - ncB + 1, nFFTB], {}], Range[1, n2 - ncB]];
        env = If[l2Min != 0., Developer`ToPackedArray[Exp[-l2Min * (t2 - t2Orig)]], 1.];
        Developer`ToPackedArray[Join @@ Table[
            e = Min[s + rowBlk - 1, n1];
            With[{nCh = e - s + 1, cSub = Developer`ToPackedArray[Transpose[Transpose[ampsMat[[s ;; e]]] * phaseVec]]},
                aReFlat = Developer`ToPackedArray[Flatten[Re[Transpose[cSub]]]];
                aImFlat = Developer`ToPackedArray[Flatten[Im[Transpose[cSub]]]];
                grids = If[shouldUseLLVM["SpreadPolyMulti", backendOpt, N[nP] * N[nCh], 800000.] && Head[fnM = getLLVMSpreadPolyMulti[]] === CompiledCodeFunction,
                    raw = Developer`ToPackedArray[fnM[xVec, aReFlat, aImFlat, nFFTB, planB["w"], planB["d"] + 1, planB["polyCoeffs"], nCh]];
                    Developer`ToPackedArray[ArrayReshape[raw, {nCh, nFFTB}]],
                    resM = wvmSpreadPolyMultiKernel[xVec, aReFlat, aImFlat, nFFTB, planB["w"], planB["d"] + 1, planB["polyCoeffs"], nCh];
                    Developer`ToPackedArray[ArrayReshape[resM[[1]] + I * resM[[2]], {nCh, nFFTB}]]
                ];
                fRows = Fourier[#, FourierParameters -> {1, 1}] & /@ grids;
                sigMat = (#[[idxSlice]] / planB["deconvVec"]) & /@ fRows;
                If[l2Min != 0., sigMat = (# * env) & /@ sigMat];
                sigMat
            ],
            {s, 1, n1, rowBlk}
        ]]
    ];

    If[methOpt === "2DDirect", Return[eval2DDir[]]];
    If[methOpt === "2DSemiSeparable",
        If[!u2 || n2 < 16 || d2 > 1.*^-12 * Max[1., Abs[l2Min]], Return[eval2DDir[]], Return[eval2DSemi[]]]
    ];

    If[methOpt =!= "2DNUFFT" && methOpt =!= "NUFFT",
        If[!u2 || n2 < 16 || d2 > 1.*^-12 * Max[1., Abs[l2Min]] || nP <= 128 || (nP * n1 * n2 <= 3*^6 && nP <= 1500),
            Return[eval2DDir[]]
        ];
        If[!u1 || d1 > 1.*^-12 * Max[1., Abs[l1Min]] || n1 < 16,
            Return[eval2DSemi[]]
        ];
    ];

    plan2 = getPlan2DLocal[n1, dt1, n2, dt2, tol, sigma, useCache];
    p1 = plan2["Plan1"]; p2 = plan2["Plan2"];
    nc1 = p1["nc"]; nFFT1 = p1["nFFT"]; tMid1 = t10 + nc1 * dt1;
    nc2 = p2["nc"]; nFFT2 = p2["nFFT"]; tMid2 = t20 + nc2 * dt2;
    cVec  = amps * Exp[I * (w1 * tMid1 + w2 * tMid2)];
    cRe   = Developer`ToPackedArray[Re[cVec]];
    cIm   = Developer`ToPackedArray[Im[cVec]];
    x1Vec = Developer`ToPackedArray[(nFFT1 / (2. * Pi)) * Mod[w1 * dt1, 2. * Pi]];
    x2Vec = Developer`ToPackedArray[(nFFT2 / (2. * Pi)) * Mod[w2 * dt2, 2. * Pi]];

    uFlat = If[shouldUseLLVM["SpreadES2D", backendOpt, N[nP], 150000.] && Head[fn = getLLVMSpreadES2D[]] === CompiledCodeFunction,
        Developer`ToPackedArray[fn[x1Vec, x2Vec, cRe, cIm, nFFT1, p1["w"], p1["beta"], nFFT2, p2["w"], p2["beta"]]],
        spreadES2DWVMLocal[x1Vec, x2Vec, cRe, cIm, nFFT1, p1["w"], p1["beta"], nFFT2, p2["w"], p2["beta"], parOpt]
    ];
    uMat = ArrayReshape[uFlat, {nFFT1, nFFT2}];
    fMat = Fourier[uMat, FourierParameters -> {1, 1}];
    i1 = Join[If[nc1 > 0, Range[nFFT1 - nc1 + 1, nFFT1], {}], Range[1, n1 - nc1]];
    i2 = Join[If[nc2 > 0, Range[nFFT2 - nc2 + 1, nFFT2], {}], Range[1, n2 - nc2]];
    sig = fMat[[i1, i2]] / Outer[Times, p1["deconvVec"], p2["deconvVec"]];
    If[l1Min != 0. || l2Min != 0.,
        env1 = Developer`ToPackedArray[Exp[-l1Min * (t1 - t1Orig)]];
        env2 = Developer`ToPackedArray[Exp[-l2Min * (t2 - t2Orig)]];
        sig = sig * Outer[Times, env1, env2]
    ];
    Developer`ToPackedArray[sig]
];

CreateStreamingSignalAccumulator[tGrid_List, lam_?NumericQ, opts___Rule] := Module[
    {tVec = Developer`ToPackedArray[N[tGrid]], nT, dt, t0, tOrigin, tol, sigma, batchCap,
     useCache, backendOpt, parOpt, nChunksOpt, plan, nFFT, nc, tMid, id},
    nT = Length[tVec];
    dt = If[nT >= 2, tVec[[2]] - tVec[[1]], 1.];
    t0 = If[nT >= 1, tVec[[1]], 0.];
    tOrigin = TimingOrigin /. {opts} /. kernelDefaults;
    If[tOrigin === Automatic, tOrigin = t0, tOrigin = N[tOrigin]];
    tol        = Max[1.*^-14, Min[0.1, N[ExpansionTolerance /. {opts} /. kernelDefaults]]];
    sigma      = Max[1.1, N[OversamplingFactor /. {opts} /. kernelDefaults]];
    batchCap   = Max[256, Round[BatchPeakCap /. {opts} /. kernelDefaults]];
    useCache   = PlanCache /. {opts} /. kernelDefaults;
    backendOpt = CompilationBackend /. {opts} /. kernelDefaults;
    parOpt     = ParallelExecution /. {opts} /. kernelDefaults;
    nChunksOpt = NumChunks /. {opts} /. kernelDefaults;
    plan = getPlan1DLocal[nT, dt, tol, sigma, False, useCache];
    nFFT = plan["nFFT"]; nc = plan["nc"]; tMid = t0 + nc * dt;
    id = ++nextStreamId;
    streamStore[id] = <|
        "tVec" -> tVec, "nT" -> nT, "dt" -> dt, "t0" -> t0, "tOrigin" -> tOrigin, "tMid" -> tMid,
        "lam" -> N[Re[lam]], "plan" -> plan, "nFFT" -> nFFT, "nc" -> nc, "batchCap" -> batchCap,
        "backend" -> backendOpt, "parOpt" -> parOpt, "nChunks" -> nChunksOpt,
        "grid" -> Developer`ToPackedArray[ConstantArray[0. + 0. I, nFFT]],
        "bagA" -> Internal`Bag[], "bagW" -> Internal`Bag[],
        "stageCount" -> 0, "totalPeaks" -> 0, "flushes" -> 0
    |>;
    StreamingSignalHandle[id]
];

flushStreamLocal[id_Integer] := Module[
    {st = streamStore[id], pA, pW, amps, omegas, plan, nFFT, dt, tMid, cVec, xVec, fn, delta, nP},
    If[st["stageCount"] === 0, Return[Null]];
    pA = Internal`BagPart[st["bagA"], All]; pW = Internal`BagPart[st["bagW"], All];
    amps   = Developer`ToPackedArray[If[Length[pA] === 1, pA[[1]], Join @@ pA]];
    omegas = Developer`ToPackedArray[If[Length[pW] === 1, pW[[1]], Join @@ pW]];
    nP     = Length[amps];
    plan   = st["plan"]; nFFT = st["nFFT"]; dt = st["dt"]; tMid = st["tMid"];
    cVec   = amps * Exp[I * (omegas * tMid)];
    xVec   = Developer`ToPackedArray[(nFFT / (2. * Pi)) * Mod[omegas * dt, 2. * Pi]];
    delta  = If[shouldUseLLVM["SpreadES", Lookup[st, "backend", Automatic], N[st["totalPeaks"]], 300000.] && Head[fn = getLLVMSpreadES[]] === CompiledCodeFunction,
        Developer`ToPackedArray[fn[xVec, Developer`ToPackedArray[Re[cVec]], Developer`ToPackedArray[Im[cVec]], nFFT, plan["w"], plan["beta"]]],
        spreadESWVMLocal[xVec, Developer`ToPackedArray[Re[cVec]], Developer`ToPackedArray[Im[cVec]], nFFT, plan["w"], plan["beta"], Lookup[st, "parOpt", True], Lookup[st, "nChunks", Automatic]]
    ];
    streamStore[id]["grid"] = Developer`ToPackedArray[st["grid"] + delta];
    streamStore[id]["bagA"] = Internal`Bag[]; streamStore[id]["bagW"] = Internal`Bag[];
    streamStore[id]["stageCount"] = 0; streamStore[id]["flushes"] = st["flushes"] + 1;
];

AccumulateStreamingPeaks[acc:StreamingSignalHandle[id_Integer], peaks_List, weight_: 1.] := Module[
    {st, w = N[weight], nAdd = Length[peaks]},
    If[nAdd === 0 || w == 0. || !KeyExistsQ[streamStore, id], Return[acc]];
    st = streamStore[id];
    Internal`StuffBag[st["bagA"], Developer`ToPackedArray[N[w * peaks[[All, 1]]] + 0. I]];
    Internal`StuffBag[st["bagW"], Developer`ToPackedArray[N[Re[peaks[[All, 2]]]]]];
    streamStore[id]["stageCount"] = st["stageCount"] + nAdd;
    streamStore[id]["totalPeaks"] = st["totalPeaks"] + nAdd;
    If[st["stageCount"] + nAdd >= st["batchCap"], flushStreamLocal[id]];
    acc
];

FinalizeStreamingSignal[acc:StreamingSignalHandle[id_Integer]] := Module[
    {st, nT, nFFT, nc, plan, lam, fAll, fCenter, sig, env},
    If[!KeyExistsQ[streamStore, id], Return[$Failed]];
    flushStreamLocal[id];
    st = streamStore[id];
    nT = st["nT"]; nFFT = st["nFFT"]; nc = st["nc"]; plan = st["plan"]; lam = st["lam"];
    If[nT === 0, Return[{}]];
    If[st["totalPeaks"] === 0, Return[Developer`ToPackedArray[ConstantArray[0. + 0. I, nT]]]];
    fAll = Fourier[st["grid"], FourierParameters -> {1, 1}];
    fCenter = If[nc > 0, Join[fAll[[nFFT - nc + 1 ;; nFFT]], fAll[[1 ;; nT - nc]]], fAll[[1 ;; nT]]];
    sig = fCenter / plan["deconvVec"];
    If[lam != 0.,
        env = Developer`ToPackedArray[Exp[-lam * (st["tVec"] - st["tOrigin"])]];
        sig = sig * env
    ];
    Developer`ToPackedArray[sig]
];

StreamingSignalStats[StreamingSignalHandle[id_Integer]] := With[{st = Lookup[streamStore, id, <||>]},
    If[st === <||>, <||>,
        <|"TotalPeaks" -> st["totalPeaks"], "Flushes" -> st["flushes"],
          "BatchPeakCap" -> st["batchCap"], "nFFT" -> st["nFFT"],
          "GridBufferBytes" -> ByteCount[st["grid"]]|>
    ]
];

GetSignalExpansionCompileStats[] := <|
    "LLVMCompiled" -> AssociationMap[(Head[llvmCache[#]] === CompiledCodeFunction) &, Keys[llvmCache]],
    "LLVMCompileTimes" -> llvmCompileTimes
|>;

WarmupSignalExpansionKernels[keys_: All] := Module[{list},
    list = Which[
        keys === All, {"Direct", "SpreadES", "SpreadPoly", "SpreadPolyMulti", "GatherPoly", "EvalCheb", "TruncatedDirect", "SpreadES2D"},
        StringQ[keys], {keys},
        ListQ[keys], Select[keys, StringQ],
        True, {}
    ];
    Do[
        Switch[k,
            "Direct", getLLVMDirect[],
            "SpreadES", getLLVMSpreadES[],
            "SpreadPoly", getLLVMSpreadPoly[],
            "SpreadPolyMulti", getLLVMSpreadPolyMulti[],
            "GatherPoly", getLLVMGatherPoly[],
            "EvalCheb", getLLVMEvalCheb[],
            "TruncatedDirect", getLLVMTruncatedDirect[],
            "SpreadES2D", getLLVMSpreadES2D[]
        ],
        {k, list}
    ];
    GetSignalExpansionCompileStats[]
];

End[];
EndPackage[];