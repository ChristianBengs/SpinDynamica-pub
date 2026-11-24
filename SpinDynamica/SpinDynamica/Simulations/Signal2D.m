(* ::Package:: *)

(* :Title: Simulations`Signal2D` *)
(* :Context: Simulations`Signal2D` *)
(* :Summary: High-performance 2D NMR simulation framework providing Signal2D with Hilbert/Liouville blocked double-diagonalization, coherence-order masks, split-echo and constant-time indirect evolution, COMPUTE2D periodic MAS with coherent locked-gamma CarouselAverage, compiled operator-atom powder ensembles, hypercomplex/TPPI quadrature, and phase cycling. *)

BeginPackage["Simulations`Signal2D`",
{"Global`",
 "SpinDynamica`",
 "General`UserLevel`",
 "General`Miscellaneous`",
 "QM`Hilbert`",
 "QM`SpinOperators`",
 "QM`Liouville`",
 "SD`Generators`",
 "SD`Propagation`",
 "Events`Events`",
 "Events`PulseSequence`",
 "Ensemble`EnsembleAverage`",
 "Signals`Signal`",
 "Signals`Signal2D`",
 "Simulations`Simulations`",
 "Simulations`Signal2DKernels`"
}];

ClearAll[
    Signal2D,
    Evolution1,
    Evolution2,
    Mixing,
    QuadratureMethod,
    QuadratureOrder,
    QuadratureSpins,
    PhaseCycle,
    CoherencePathway,
    MagnusSubsteps
];

Signal2D::usage = "Signal2D[{timing1, timing2}, opts] simulates a 2D NMR time-domain signal (returned as a SignalData2D object).\n" <>
"timing1 (indirect t1) and timing2 (direct t2) use Signal1D timing syntax: {tmax, dt}, {tmin, tmax, dt}, {tmin, tmax, {n}}, or {{omegaSW, n}} (where omegaSW is spectral width in rad s^-1).\n" <>
"Lifecycle stages obey Chronology (\"LeftToRight\" by default): InitialDensityOperator -> Preparation -> Evolution1 -> Mixing -> Evolution2 (detected with Observable).\n" <>
"Evolution1 may be None (free evolution under BackgroundGenerator), an Operator/Superoperator, {\"Echo\", refocusEvents}, {\"ConstantTime\", Tct, refocusEvents}, or a PeriodicFunction / TimeDependentGeneratorFunction.\n" <>
"QuadratureMethod supports \"States\" (default), \"StatesTPPI\", \"EchoAntiEcho\", \"TPPI\", and \"PhaseModulated\".\n" <>
"SignalCalculationMethod supports Automatic, \"DoubleDiagonalization\", \"SteppedPropagator\", \"COMPUTE2D\", and \"Direct\".";
Evolution1::usage = "Evolution1 is an option for Signal2D specifying evolution during the indirect t1 period: None (background only), a constant generator, {\"Echo\", refocusEvents}, {\"ConstantTime\", Tct, refocusEvents}, {\"Stepped\", dt1Events}, or a time-dependent/periodic generator.";
Evolution2::usage = "Evolution2 is an option for Signal2D specifying evolution during direct detection t2: None (background only), a constant generator, {\"Stepped\", dt2Events}, or a time-dependent/periodic generator.";
Mixing::usage = "Mixing is an option for Signal2D specifying the mixing event sequence between indirect t1 evolution and direct t2 detection (pulses, delays, CoherenceOrderFiltration, etc.).";
QuadratureMethod::usage = "QuadratureMethod is an option for Signal2D selecting indirect-dimension quadrature encoding: \"States\" (default), \"StatesTPPI\", \"EchoAntiEcho\", \"TPPI\", or \"PhaseModulated\".";
QuadratureOrder::usage = "QuadratureOrder is an option for Signal2D specifying the indirect coherence order magnitude p used for quadrature phase shifts and Echo/AntiEcho selection (default 1).";
QuadratureSpins::usage = "QuadratureSpins is an option for Signal2D specifying which spins receive the indirect quadrature phase shift or coherence selection (default All).";
PhaseCycle::usage = "PhaseCycle is an option for Signal2D specifying a list of phase-cycle steps (<|\"Preparation\" -> ..., \"Mixing\" -> ..., \"ReceiverPhase\" -> ..., \"Weight\" -> ...|>).";
CoherencePathway::usage = "CoherencePathway is an option for Signal2D specifying optional coherence-order filtration p1 (scalar/list after Preparation) or {p1, pMix} (after Preparation and Mixing).";
MagnusSubsteps::usage = "MagnusSubsteps is an option for Signal2D specifying the number of 4th-order Gauss-Legendre Magnus (CF4) integration substeps per rotor sub-interval during periodic COMPUTE2D propagation (default Automatic).";

Signal2D::badtimes = "Invalid 2D timing specification `1`.";
Signal2D::invalobs = "Invalid Observable specification `1`.";
Signal2D::invalbg = "BackgroundGenerator must be a generator or a 4-element list {bgPrep, bgEv1, bgMix, bgEv2}.";
Signal2D::badquad = "Invalid QuadratureMethod specification `1`.";
Signal2D::invalchron = "Chronology must be \"LeftToRight\" or \"RightToLeft\".";
Signal2D::badcalc = "Requested SignalCalculationMethod `1` cannot be used for this problem: `2`";
Signal2D::nocarav = "CarouselAverage -> True (or > 1) requires periodic evolution with SignalCalculationMethod -> \"COMPUTE2D\".";
Signal2D::fail = "Signal2D simulation failed: `1`.";
Signal2D::reportscm = "Using SignalCalculationMethod -> `1` (`2` space).";

Options[Signal2D] = {
    OperatorBasis :> ZeemanKetBraOperatorBasis[],
    Basis :> ZeemanBasis[],
    InitialDensityOperator :> opI["z"],
    Preparation :> Pulse[{Pi/2, "x"}],
    Evolution1 -> None,
    Mixing -> None,
    Evolution2 -> None,
    Observable -> All,
    NormalizationFactor -> None,
    BackgroundGenerator -> None,
    QuadratureMethod -> "States",
    QuadratureOrder -> 1,
    QuadratureSpins -> All,
    PhaseCycle -> None,
    CoherencePathway -> Automatic,
    LineBroadening -> None,
    SignalCalculationMethod -> Automatic,
    SignalGeneratorPeriod -> Automatic,
    InitialTimePoint -> 0.,
    CarouselAverage -> False,
    MagnusSubsteps -> Automatic,
    MaxBatchSize -> Infinity,
    UseHilbertSpaceWhenPossible -> True,
    ForceEvenNumberOfPoints -> True,
    EnsembleAverage -> False,
    Parallel -> False,
    Chronology :> $Chronology,
    ReportSignalCalculationMethod -> False
};

SetAttributes[Signal2D, HoldRest];

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

toPackedRealVec[v_] := If[
    Developer`PackedArrayQ[v, Real],
    v,
    Developer`ToPackedArray[N[Re[Normal[v]]]]
];

parsePointCount[k_Integer] /; k >= 2 := k;
parsePointCount[k_String] := Switch[k,
    "1k" | "1K", 1024,
    "2k" | "2K", 2048,
    "4k" | "4K", 4096,
    "8k" | "8K", 8192,
    "16k" | "16K", 16384,
    "32k" | "32K", 32768,
    "64k" | "64K", 65536,
    "128k" | "128K", 131072,
    "256k" | "256K", 262144,
    "512k" | "512K", 524288,
    _, Message[Signal2D::badtimes, k]; $Failed
];
parsePointCount[other_] := (Message[Signal2D::badtimes, other]; $Failed);

parseTimingAxis[spec_, forceEvenQ_?BooleanQ] := Module[
    {norm, tmin, tmax, dt, n, sw, count},
    norm = Which[
        MatchQ[spec, {_?NumericQ, _?NumericQ, _?NumericQ}],
            {tmin, tmax, dt} = N[Re[spec]];
            If[dt <= 0. || tmax < tmin, Return[$Failed]];
            n = Round[(tmax - tmin) / dt] + 1;
            {tmin, tmin + (n - 1) * dt, dt, n},
        MatchQ[spec, {_?NumericQ, _?NumericQ}],
            Return[parseTimingAxis[Prepend[spec, 0], forceEvenQ], Module],
        MatchQ[spec, {_?NumericQ, _?NumericQ, {_}}],
            tmin = N[Re[spec[[1]]]]; tmax = N[Re[spec[[2]]]];
            count = parsePointCount[spec[[3, 1]]];
            If[count === $Failed || tmax <= tmin, Return[$Failed]];
            dt = (tmax - tmin) / N[count - 1];
            {tmin, tmax, dt, count},
        MatchQ[spec, {_?NumericQ, {_}}],
            Return[parseTimingAxis[Prepend[spec, 0], forceEvenQ], Module],
        MatchQ[spec, {_?NumericQ, {_?NumericQ, _}}],
            tmin = N[Re[spec[[1]]]]; sw = N[Re[spec[[2, 1]]]];
            count = parsePointCount[spec[[2, 2]]];
            If[count === $Failed || sw <= 0., Return[$Failed]];
            dt = (2. * Pi) / sw;
            {tmin, tmin + N[count - 1] * dt, dt, count},
        MatchQ[spec, {{_?NumericQ, _}}],
            Return[parseTimingAxis[Prepend[spec, 0], forceEvenQ], Module],
        True,
            Return[$Failed]
    ];
    If[!ListQ[norm], Return[$Failed]];
    {tmin, tmax, dt, n} = norm;
    If[forceEvenQ && OddQ[n] && n > 1,
        n = n - 1;
        tmax = tmin + N[n - 1] * dt
    ];
    If[n < 2 || dt <= 0., Return[$Failed]];
    <|
        "Timing" -> {N[tmin], N[tmax], N[dt]},
        "N" -> n,
        "TimeVec" -> toPackedRealVec[N[tmin] + N[dt] * Range[0, n - 1]]
    |>
];

(* Event normalization respecting Chronology without double-reversing PulseSequence *)
singleTimedEventQ[ev_] := MatchQ[ev,
    {_?OperatorQ | _?SuperoperatorQ | None | NullOperator[] | NullSuperoperator[] | 0 | 0. |
     _Function | _PeriodicFunction | _TimeDependentGeneratorFunction, _?NumericQ, ___Rule}
] || MatchQ[ev, {_?NumericQ, ___Rule}];

normalizeEventsChronological[evSpec_, chron_String] := Which[
    evSpec === None || evSpec === {} || evSpec === False ||
    MatchQ[evSpec, UnitySuperoperator[] | NullOperator[] | NullSuperoperator[]],
        {},
    Head[evSpec] === PulseSequence,
        Module[{std = StandardizeEventSpecification[evSpec, Chronology -> chron], rawList},
            If[std === $Failed, Return[$Failed]];
            rawList = Reverse[EventSpecifications[std][[2]]];
            If[chron === "RightToLeft", rawList = Reverse[rawList]];
            Map[If[MatchQ[#, {_}], First[#], #] &, rawList]
        ],
    Head[evSpec] === EventSpecifications,
        Map[If[MatchQ[#, {_}], First[#], #] &, Reverse[evSpec[[2]]]],
    singleTimedEventQ[evSpec],
        {evSpec},
    ListQ[evSpec],
        If[chron === "RightToLeft", Reverse[evSpec], evSpec],
    ValidEventSpecificationQ[evSpec] || MatchQ[evSpec, _CoherenceOrderFiltration | _CoherenceOrderFiltrationSuperoperator | _RotationOperator | _RotationSuperoperator | _Pulse | _Delay | _ShapedPulse],
        {evSpec},
    True,
        $Failed
];

filterEventQ[ev_] := !FreeQ[ev, _CoherenceOrderFiltration | _CoherenceOrderFiltrationSuperoperator];

extractFiltersFromEvent[ev_] := Module[{std, filts},
    Which[
        MatchQ[ev, _CoherenceOrderFiltrationSuperoperator],
            filts = Cases[ev, {{Except[_List]..}, {Except[_List]..}}, {1}],
        MatchQ[ev, _CoherenceOrderFiltration],
            std = EventSpecifications[ev];
            filts = Cases[std, _CoherenceOrderFiltrationSuperoperator, Infinity];
            If[filts === {}, Return[$Failed]];
            filts = Flatten[Cases[#, {{Except[_List]..}, {Except[_List]..}}, {1}] & /@ filts, 1],
        True,
            Return[$Failed]
    ];
    If[filts === {}, $Failed, filts]
];

nullGenQ[g_] := g === None || g === False || g === NullOperator[] || g === NullSuperoperator[] || MatchQ[g, 0 | 0.];

periodicGenQ[g_] := !nullGenQ[g] && !FreeQ[g, _PeriodicFunction | _TimeDependentGeneratorFunction];

extractGeneratorPeriod[g_] := Which[
    nullGenQ[g], All,
    MatchQ[g, HoldPattern[PeriodicFunction[_, T_?NumericQ, ___]]],
        If[TrueQ[g[[2]] > 0], N[Re[g[[2]]]], None],
    MatchQ[g, HoldPattern[TimeDependentGeneratorFunction[_, _, T_, expr_, ___Rule]]],
        With[{ts = g[[1]], T = g[[3]], ex = g[[4]]},
            Which[
                NumericQ[T] && T > 0, N[Re[T]],
                T === None && FreeQ[ex, ts], All,
                True, None
            ]
        ],
    MatchQ[g, _CommutationSuperoperator | HoldPattern[Superoperator[___]]],
        extractGeneratorPeriod[First[g]],
    MatchQ[g, _Plus | _Times],
        CommonPeriod @@ (extractGeneratorPeriod /@ (List @@ g)),
    OperatorQ[g] || SuperoperatorQ[g],
        If[!periodicGenQ[g], All, None],
    MatchQ[g, _Function],
        With[{std = Quiet[StandardizeGenerator[g]]},
            If[MatchQ[std, _TimeDependentGeneratorFunction | _Operator | _Superoperator | None | NullGeneratorPattern],
                extractGeneratorPeriod[std],
                None
            ]
        ],
    True,
        If[!periodicGenQ[g] && FreeQ[g, _Function], All, None]
];

evalGenAtNumericTime[g_, tVal_?NumericQ] := Quiet[
    g /. {
        PeriodicFunction[ts_, _, b_] :> (b /. ts -> tVal),
        TimeDependentGeneratorFunction[ts_, _, _, b_, ___Rule] :> (b /. ts -> tVal),
        f_Function :> f[tVal]
    }
];

hilbertGenQ[g_] := Which[
    nullGenQ[g] || OperatorQ[g], True,
    MatchQ[g, _CommutationSuperoperator | HoldPattern[Superoperator[_CommutationSuperoperator, ___]]],
        hilbertGenQ[First[g]],
    MatchQ[g, HoldPattern[PeriodicFunction[_, _, _]] | HoldPattern[TimeDependentGeneratorFunction[_, _, _, _, ___Rule]] | _Plus | _Times | _Function],
        With[{sampled = evalGenAtNumericTime[g, 0.1234567]},
            FreeQ[sampled, _PeriodicFunction | _TimeDependentGeneratorFunction | _Superoperator] &&
            (nullGenQ[sampled] || OperatorQ[sampled] ||
             MatchQ[sampled, _CommutationSuperoperator | HoldPattern[Superoperator[_CommutationSuperoperator, ___]]])
        ],
    True, False
];

hilbertEventListQ[evs_List] := AllTrue[evs, Function[ev,
    Which[
        MatchQ[ev, _CoherenceOrderFiltration | _CoherenceOrderFiltrationSuperoperator], True,
        MatchQ[ev, _Pulse | _Delay | _ShapedPulse | _RotationOperator | _RotationSuperoperator | UnitySuperoperator[] | NullSuperoperator[]], True,
        MatchQ[ev, {_?OperatorQ | None | NullOperator[] | 0 | 0. | _PeriodicFunction | _TimeDependentGeneratorFunction, _?NumericQ, ___Rule}],
            hilbertGenQ[First[ev]],
        MatchQ[ev, {_?NumericQ, ___Rule}], True,
        MatchQ[ev, _CommutationSuperoperator | HoldPattern[Superoperator[_CommutationSuperoperator, ___]]], hilbertGenQ[First[ev]],
        ListQ[ev], hilbertEventListQ[ev],
        True, False
    ]
]];

hilbertStageQ[evSpec_, chron_String] := Module[{evs},
    evs = normalizeEventsChronological[evSpec, chron];
    If[evs === $Failed, Return[False]];
    hilbertEventListQ[evs]
];

stageEventDuration[ev_] := Which[
    ev === {} || ev === None || ev === False, 0.,
    MatchQ[ev, _CoherenceOrderFiltration | _CoherenceOrderFiltrationSuperoperator |
        _RotationOperator | _RotationSuperoperator | UnitySuperoperator[] | NullSuperoperator[]], 0.,
    MatchQ[ev, Delay[_?NumericQ]], N[Re[First[ev]]],
    MatchQ[ev, {_, _?NumericQ, ___Rule}], N[Re[ev[[2]]]],
    MatchQ[ev, {_?NumericQ, ___Rule}], N[Re[First[ev]]],
    ListQ[ev], Total[stageEventDuration /@ ev],
    True, Quiet[With[{d = Duration[ev]}, If[NumericQ[d], N[Re[d]], 0.]]]
];

generatorToHilbertMatrix[g_, basis_, dim_Integer] := Module[{std, opMat},
    If[nullGenQ[g], Return[toPackedComplexMat[ConstantArray[0. + 0. I, {dim, dim}]]]];
    std = If[OperatorQ[g] || MatchQ[g, _CommutationSuperoperator | HoldPattern[Superoperator[_CommutationSuperoperator, ___]]],
        g,
        StandardizeGenerator[g]
    ];
    Which[
        nullGenQ[std],
            toPackedComplexMat[ConstantArray[0. + 0. I, {dim, dim}]],
        OperatorQ[std],
            opMat = Quiet[MatrixRepresentation[std, basis]];
            If[!MatrixQ[opMat, NumericQ], $Failed, toPackedComplexMat[opMat]],
        MatchQ[std, HoldPattern[Superoperator[_CommutationSuperoperator, ___]]],
            generatorToHilbertMatrix[First[std], basis, dim],
        MatchQ[std, _CommutationSuperoperator],
            (* CommutationSuperoperator[-I * H] represents evolution U = Exp[-I H t] *)
            opMat = Quiet[MatrixRepresentation[First[std], basis]];
            If[!MatrixQ[opMat, NumericQ], $Failed, toPackedComplexMat[I * opMat]],
        True,
            $Failed
    ]
];

generatorToLiouvilleMatrix[g_, opbasis_, basis_, dimL_Integer] := Module[{std, mat},
    If[nullGenQ[g], Return[toPackedComplexMat[ConstantArray[0. + 0. I, {dimL, dimL}]]]];
    std = Which[
        OperatorQ[g], CommutationSuperoperator[-I * g],
        SuperoperatorQ[g], g,
        True, StandardizeGenerator[g]
    ];
    If[nullGenQ[std], Return[toPackedComplexMat[ConstantArray[0. + 0. I, {dimL, dimL}]]]];
    If[OperatorQ[std], std = CommutationSuperoperator[-I * std]];
    mat = Quiet[SuperoperatorMatrixRepresentation[std, opbasis, basis]];
    If[!MatrixQ[mat, NumericQ], $Failed, toPackedComplexMat[mat]]
];

buildHilbertPropagatorFromSegment[seg_List, bg_, basis_, dim_Integer, Optional[tStart_?NumericQ, 0.]] := Module[
    {cleanSeg, uOp, uMat, totalH, dur, segStaticQ},
    If[seg === {}, Return[toPackedComplexMat[IdentityMatrix[dim]]]];
    segStaticQ = !periodicGenQ[bg] && FreeQ[seg, _PeriodicFunction | _TimeDependentGeneratorFunction | _Function];
    (* Fast analytical MatrixExp when single static {H, dur} or Delay[dur] *)
    If[Length[seg] === 1,
        With[{e = First[seg]},
            Which[
                segStaticQ && MatchQ[e, Delay[_?NumericQ]],
                    dur = N[Re[First[e]]];
                    totalH = generatorToHilbertMatrix[bg, basis, dim];
                    If[totalH === $Failed, Return[$Failed]];
                    Return[toPackedComplexMat[MatrixExp[-I * totalH * dur]]],
                segStaticQ && MatchQ[e, {_?OperatorQ | None | NullOperator[] | 0 | 0., _?NumericQ}],
                    dur = N[Re[e[[2]]]];
                    totalH = generatorToHilbertMatrix[e[[1]], basis, dim] + generatorToHilbertMatrix[bg, basis, dim];
                    If[!MatrixQ[totalH, NumericQ], Return[$Failed]];
                    Return[toPackedComplexMat[MatrixExp[-I * totalH * dur]]],
                MatchQ[e, _RotationOperator],
                    uMat = Quiet[MatrixRepresentation[e, basis]];
                    If[MatrixQ[uMat, NumericQ], Return[toPackedComplexMat[uMat]]],
                MatchQ[e, _RotationSuperoperator],
                    uMat = Quiet[MatrixRepresentation[RotationOperator @@ e, basis]];
                    If[MatrixQ[uMat, NumericQ], Return[toPackedComplexMat[uMat]]]
            ]
        ]
    ];
    cleanSeg = EventSpecifications[seg, Chronology -> "LeftToRight"];
    uOp = Quiet[NPropagationOperator[cleanSeg, Basis -> basis, BackgroundGenerator -> bg, InitialTimePoint -> N[Re[tStart]]]];
    If[uOp === $Failed || !OperatorQ[uOp], Return[$Failed]];
    uMat = Quiet[MatrixRepresentation[uOp, basis]];
    If[!MatrixQ[uMat, NumericQ], $Failed, toPackedComplexMat[uMat]]
];

(* Decompose chronological event list into ordered Unitary and Filter steps in Hilbert space *)
buildHilbertStageSteps[chronEvents_List, bg_, basis_, system_List, Optional[tStart_?NumericQ, 0.]] := Module[
    {dim = Round[Times @@ (2 * (Last /@ system) + 1)], steps = {}, seg = {}, flushSeg, filts, maskSpec,
     tSeg = N[Re[tStart]], dSeg},
    flushSeg[] := If[seg =!= {},
        Module[{u = buildHilbertPropagatorFromSegment[seg, bg, basis, dim, tSeg]},
            If[u === $Failed, Return[$Failed, Module]];
            dSeg = stageEventDuration[seg];
            tSeg += dSeg;
            If[Length[steps] > 0 && Last[steps][[1]] === "Unitary",
                With[{uNew = toPackedComplexMat[u . Last[steps][[2]]]},
                    steps[[-1]] = {"Unitary", uNew, ConjugateTranspose[uNew]}
                ],
                AppendTo[steps, {"Unitary", u, ConjugateTranspose[u]}]
            ];
            seg = {}
        ]
    ];
    Do[
        If[filterEventQ[ev],
            If[flushSeg[] === $Failed, Return[$Failed, Module]];
            filts = extractFiltersFromEvent[ev];
            If[filts === $Failed, Return[$Failed, Module]];
            maskSpec = BuildHilbertCoherenceMask[filts, system, basis];
            AppendTo[steps, {"Filter", maskSpec}],
            AppendTo[seg, ev]
        ],
        {ev, chronEvents}
    ];
    If[flushSeg[] === $Failed, Return[$Failed]];
    steps
];

buildLiouvilleStageMatrix[chronEvents_List, bg_, opbasis_, basis_, dimL_Integer, Optional[tStart_?NumericQ, 0.]] := Module[
    {cleanEvs, pSup, pMat, totalL, dur, staticStageQ},
    If[chronEvents === {}, Return[toPackedComplexMat[IdentityMatrix[dimL]]]];
    staticStageQ = !periodicGenQ[bg] && FreeQ[chronEvents, _PeriodicFunction | _TimeDependentGeneratorFunction | _Function];
    If[Length[chronEvents] === 1,
        With[{e = First[chronEvents]},
            Which[
                staticStageQ && MatchQ[e, Delay[_?NumericQ]],
                    dur = N[Re[First[e]]];
                    totalL = generatorToLiouvilleMatrix[bg, opbasis, basis, dimL];
                    If[totalL === $Failed, Return[$Failed]];
                    Return[toPackedComplexMat[MatrixExp[totalL * dur]]],
                staticStageQ && MatchQ[e, {_, _?NumericQ}] && (OperatorQ[e[[1]]] || SuperoperatorQ[e[[1]]] || nullGenQ[e[[1]]]),
                    dur = N[Re[e[[2]]]];
                    totalL = generatorToLiouvilleMatrix[e[[1]], opbasis, basis, dimL] + generatorToLiouvilleMatrix[bg, opbasis, basis, dimL];
                    If[!MatrixQ[totalL, NumericQ], Return[$Failed]];
                    Return[toPackedComplexMat[MatrixExp[totalL * dur]]],
                MatchQ[e, _CoherenceOrderFiltrationSuperoperator | _CoherenceOrderFiltration],
                    With[{filts = extractFiltersFromEvent[e]},
                        If[filts === $Failed, Return[$Failed]];
                        pSup = Fold[
                            SuperoperatorMatrixRepresentation[CoherenceOrderFiltrationSuperoperator[#2[[1]], #2[[2]]], opbasis, basis] . #1 &,
                            IdentityMatrix[dimL],
                            filts
                        ];
                        Return[toPackedComplexMat[pSup]]
                    ],
                MatchQ[e, _RotationSuperoperator | _RotationOperator],
                    pSup = If[Head[e] === RotationOperator, RotationSuperoperator @@ e, e];
                    pMat = Quiet[SuperoperatorMatrixRepresentation[pSup, opbasis, basis]];
                    If[MatrixQ[pMat, NumericQ], Return[toPackedComplexMat[pMat]]]
            ]
        ]
    ];
    cleanEvs = EventSpecifications[chronEvents, Chronology -> "LeftToRight"];
    pSup = Quiet[NPropagationSuperoperator[cleanEvs, OperatorBasis -> opbasis, Basis -> basis, BackgroundGenerator -> bg, InitialTimePoint -> N[Re[tStart]]]];
    If[pSup === $Failed, Return[$Failed]];
    pMat = Quiet[SuperoperatorMatrixRepresentation[pSup, opbasis, basis]];
    If[!MatrixQ[pMat, NumericQ], $Failed, toPackedComplexMat[pMat]]
];

resolveObservable[obsSpec_, system_List] := Module[{labels = First /@ system},
    Which[
        obsSpec === All, -(I / 2.) * opI[labels, "-"],
        ValidSpinQ[obsSpec, system], -(I / 2.) * opI[{obsSpec}, "-"],
        ListQ[obsSpec] && obsSpec =!= {} && ValidSpinQ[obsSpec, system], -(I / 2.) * opI[obsSpec, "-"],
        OperatorQ[obsSpec], obsSpec,
        True, $Failed
    ]
];

resolveQuadratureSpins[qsSpec_, system_List] := Module[{labels = First /@ system},
    Which[
        qsSpec === All, labels,
        ValidSpinQ[qsSpec, system], {qsSpec},
        ListQ[qsSpec] && qsSpec =!= {} && ValidSpinQ[qsSpec, system], qsSpec,
        True, $Failed
    ]
];

resolveCoherencePathway[cpSpec_, qSpins_List, obsSpins_List] := Which[
    cpSpec === Automatic || cpSpec === None || cpSpec === False,
        {None, None},
    IntegerQ[cpSpec] || MatchQ[cpSpec, {_Integer ..}] && Length[cpSpec] =!= 2,
        (* Wait: ambiguous 2-element flat integer list {-1, -1} means {p1, pMix}; use {{-2, 2}} to pass multiple p1 orders! *)
        {If[ListQ[cpSpec], cpSpec, {cpSpec}], None},
    MatchQ[cpSpec, {_Integer, _Integer}],
        {{cpSpec[[1]]}, {cpSpec[[2]]}},
    MatchQ[cpSpec, {{_Integer ..}}],
        {cpSpec[[1]], None},
    MatchQ[cpSpec, {_Integer | {_Integer ..}, _Integer | {_Integer ..}}],
        {If[ListQ[cpSpec[[1]]], cpSpec[[1]], {cpSpec[[1]]}],
         If[ListQ[cpSpec[[2]]], cpSpec[[2]], {cpSpec[[2]]}]},
    True,
        $Failed
];

staticGenQ[ev_] := (nullGenQ[ev] || OperatorQ[ev] || SuperoperatorQ[ev]) && !periodicGenQ[ev];

(* Parse Evolution1 into static, periodic, echo, constant-time, stepped, or direct *)
parseEvolution1[ev1_, chron_String] := Which[
    staticGenQ[ev1],
        <|"Mode" -> "Standard", "Gen1" -> ev1, "Kind" -> "Static"|>,
    periodicGenQ[ev1],
        <|"Mode" -> "Standard", "Gen1" -> ev1, "Kind" -> "Periodic", "Period" -> extractGeneratorPeriod[ev1]|>,
    MatchQ[ev1, {"Echo", _}],
        With[{ref = normalizeEventsChronological[ev1[[2]], chron]},
            If[ref === $Failed, $Failed,
                <|"Mode" -> "Echo", "Gen1a" -> None, "Gen1b" -> None, "RefocusEvents" -> ref, "Kind" -> "StaticSplit"|>
            ]
        ],
    MatchQ[ev1, {"Echo", _, {_, _}}] && staticGenQ[ev1[[3, 1]]] && staticGenQ[ev1[[3, 2]]],
        With[{ref = normalizeEventsChronological[ev1[[2]], chron]},
            If[ref === $Failed, $Failed,
                <|"Mode" -> "Echo", "Gen1a" -> ev1[[3, 1]], "Gen1b" -> ev1[[3, 2]], "RefocusEvents" -> ref, "Kind" -> "StaticSplit"|>
            ]
        ],
    MatchQ[ev1, {"ConstantTime", _?NumericQ, _}],
        With[{ref = normalizeEventsChronological[ev1[[3]], chron]},
            If[ref === $Failed, $Failed,
                <|"Mode" -> "ConstantTime", "CTDuration" -> N[Re[ev1[[2]]]], "Gen1a" -> None, "Gen1b" -> None, "RefocusEvents" -> ref, "Kind" -> "StaticSplit"|>
            ]
        ],
    MatchQ[ev1, {"ConstantTime", _?NumericQ, _, {_, _}}] && staticGenQ[ev1[[4, 1]]] && staticGenQ[ev1[[4, 2]]],
        With[{ref = normalizeEventsChronological[ev1[[3]], chron]},
            If[ref === $Failed, $Failed,
                <|"Mode" -> "ConstantTime", "CTDuration" -> N[Re[ev1[[2]]]], "Gen1a" -> ev1[[4, 1]], "Gen1b" -> ev1[[4, 2]], "RefocusEvents" -> ref, "Kind" -> "StaticSplit"|>
            ]
        ],
    MatchQ[ev1, {"Stepped", _}],
        With[{sevs = normalizeEventsChronological[ev1[[2]], chron]},
            If[sevs === $Failed, $Failed,
                <|"Mode" -> "Standard", "StepEvents" -> sevs, "Kind" -> "Stepped"|>
            ]
        ],
    Head[ev1] === Function || MatchQ[ev1, _TimeDependentGeneratorFunction | _PeriodicFunction],
        <|"Mode" -> "Standard", "DirectFn" -> ev1, "Kind" -> "Direct"|>,
    True,
        $Failed
];

parseEvolution2[ev2_, chron_String] := Which[
    staticGenQ[ev2],
        <|"Gen2" -> ev2, "Kind" -> "Static"|>,
    periodicGenQ[ev2],
        <|"Gen2" -> ev2, "Kind" -> "Periodic", "Period" -> extractGeneratorPeriod[ev2]|>,
    MatchQ[ev2, {"Stepped", _}],
        With[{sevs = normalizeEventsChronological[ev2[[2]], chron]},
            If[sevs === $Failed, $Failed,
                <|"StepEvents" -> sevs, "Kind" -> "Stepped"|>
            ]
        ],
    Head[ev2] === Function || MatchQ[ev2, _TimeDependentGeneratorFunction | _PeriodicFunction],
        <|"DirectFn" -> ev2, "Kind" -> "Direct"|>,
    True,
        $Failed
];

(* Apply Rz[qSpins, phi] rotation or coherence filter in Hilbert/Liouville *)
applySpinZRotation[state_, space_String, qSpins_List, phi_?NumericQ, basis_, opbasis_] := Module[{m},
    If[Chop[N[phi]] == 0., Return[state]];
    If[space === "Hilbert",
        m = toPackedComplexMat[MatrixRepresentation[RotationOperator[qSpins, {phi, "z"}], basis]];
        toPackedComplexMat[m . state . ConjugateTranspose[m]],
        m = toPackedComplexMat[SuperoperatorMatrixRepresentation[RotationSuperoperator[qSpins, {phi, "z"}], opbasis, basis]];
        toPackedComplexVec[m . state]
    ]
];

applySpinFilter[state_, space_String, spins_List, orders_List, system_List, basis_, opbasis_] := Module[{m, mask},
    If[space === "Hilbert",
        mask = BuildHilbertCoherenceMask[{{spins, orders}}, system, basis];
        ApplyHilbertMixSequence[state, {{"Filter", mask}}],
        m = toPackedComplexMat[SuperoperatorMatrixRepresentation[CoherenceOrderFiltrationSuperoperator[spins, orders], opbasis, basis]];
        toPackedComplexVec[m . state]
    ]
];

(* Direct Reference Kernel via NPropagate + Trajectory with exact global InitialTimePoint tracking *)
evaluateDirect2D[spec_Association] := Module[
    {space = spec["Space"], components = spec["Components"],
     ev1Spec = spec["Ev1Spec"], ev2Spec = spec["Ev2Spec"],
     bg1 = spec["Bg1"], bgMix = spec["BgMix"], bg2 = spec["Bg2"],
     t1Vec = spec["Time1"], t2Vec = spec["Time2"], n1, n2,
     qObs = spec["QObs"], normDenom = spec["NormDenom"], qNormSq = spec["QNormSq"],
     basis = spec["Basis"], opbasis = spec["OperatorBasis"],
     qSpins = spec["QuadratureSpins"], qOrder = spec["QuadratureOrder"], system = spec["System"],
     t0 = N[Re[Lookup[spec, "T0", 0.]]],
     tauP = N[Re[Lookup[spec, "PrepDuration", 0.]]],
     tauM = N[Re[Lookup[spec, "MixDuration", 0.]]],
     gammaShifts = N[Re[Lookup[spec, "GammaTimeShifts", {0.}]]], nG,
     tmin2, tmax2, dt2, useHilbQ, effNorm, results = <||>,
     stateToOp, eval1StateAtT1},

    n1 = Length[t1Vec];
    n2 = Length[t2Vec];
    nG = Length[gammaShifts];
    {tmin2, tmax2, dt2} = spec["Timing2"];
    useHilbQ = (space === "Hilbert");
    effNorm = spec["NormFactor"];

    stateToOp[st_] := If[space === "Hilbert",
        Operator[Normal[st], basis, CheckOperatorType -> False, CheckMatrixDimensions -> False],
        Operator[Normal[st], opbasis, basis, CheckOperatorType -> False, CheckMatrixDimensions -> False]
    ];

    (* Evolves across t1 starting at t1Start and returns {midOp, tEnd1} *)
    eval1StateAtT1[rhoOp_, t1_?NumericQ, t1Start_?NumericQ] := Module[
        {mode = ev1Spec["Mode"], g1, g1a, g1b, halfT, tCT, evsA, evsB, sMid, evs1, durRef, ta, tb},
        Which[
            t1 == 0. && mode === "Standard",
                {rhoOp, t1Start},
            mode === "Standard",
                If[MemberQ[{"Static", "Periodic"}, ev1Spec["Kind"]],
                    g1 = ev1Spec["Gen1"];
                    {NPropagate[{If[nullGenQ[g1], Delay[t1], {g1, t1}]},
                        Basis -> basis, OperatorBasis -> opbasis, BackgroundGenerator -> bg1,
                        InitialTimePoint -> t1Start,
                        UseHilbertSpaceWhenPossible -> useHilbQ][rhoOp], t1Start + t1},
                    (* Function or Stepped during t1 *)
                    evs1 = If[ev1Spec["Kind"] === "Direct",
                        With[{df = ev1Spec["DirectFn"]},
                            If[Head[df] === Function,
                                With[{res = df[t1]}, If[ListQ[res], res, {res}]],
                                {{df, t1}}
                            ]
                        ],
                        Flatten[ConstantArray[ev1Spec["StepEvents"], Max[1, Round[t1 / spec["Timing1"][[3]]]]], 1]
                    ];
                    {NPropagate[EventSpecifications[evs1, Chronology -> "LeftToRight"], Basis -> basis, OperatorBasis -> opbasis,
                        BackgroundGenerator -> bg1, InitialTimePoint -> t1Start,
                        UseHilbertSpaceWhenPossible -> (useHilbQ && !AnyTrue[evs1, filterEventQ])][rhoOp],
                     t1Start + stageEventDuration[evs1]}
                ],
            mode === "Echo",
                halfT = 0.5 * t1;
                g1a = ev1Spec["Gen1a"]; g1b = ev1Spec["Gen1b"];
                durRef = stageEventDuration[ev1Spec["RefocusEvents"]];
                evsA = If[halfT == 0., {}, {If[nullGenQ[g1a], Delay[halfT], {g1a, halfT}]}];
                evsB = If[halfT == 0., {}, {If[nullGenQ[g1b], Delay[halfT], {g1b, halfT}]}];
                sMid = If[evsA === {}, rhoOp,
                    NPropagate[EventSpecifications[evsA, Chronology -> "LeftToRight"], Basis -> basis, OperatorBasis -> opbasis,
                        BackgroundGenerator -> bg1, InitialTimePoint -> t1Start, UseHilbertSpaceWhenPossible -> useHilbQ][rhoOp]];
                sMid = If[ev1Spec["RefocusEvents"] === {}, sMid,
                    NPropagate[EventSpecifications[ev1Spec["RefocusEvents"], Chronology -> "LeftToRight"], Basis -> basis, OperatorBasis -> opbasis,
                        BackgroundGenerator -> bg1, InitialTimePoint -> t1Start + halfT,
                        UseHilbertSpaceWhenPossible -> (useHilbQ && !AnyTrue[ev1Spec["RefocusEvents"], filterEventQ])][sMid]];
                {If[evsB === {}, sMid,
                    NPropagate[EventSpecifications[evsB, Chronology -> "LeftToRight"], Basis -> basis, OperatorBasis -> opbasis,
                        BackgroundGenerator -> bg1, InitialTimePoint -> t1Start + halfT + durRef,
                        UseHilbertSpaceWhenPossible -> useHilbQ][sMid]],
                 t1Start + t1 + durRef},
            mode === "ConstantTime",
                tCT = ev1Spec["CTDuration"];
                g1a = ev1Spec["Gen1a"]; g1b = ev1Spec["Gen1b"];
                durRef = stageEventDuration[ev1Spec["RefocusEvents"]];
                ta = 0.5 * (tCT + t1); tb = 0.5 * (tCT - t1);
                evsA = If[ta == 0., {}, {If[nullGenQ[g1a], Delay[ta], {g1a, ta}]}];
                evsB = If[tb == 0., {}, {If[nullGenQ[g1b], Delay[tb], {g1b, tb}]}];
                sMid = If[evsA === {}, rhoOp,
                    NPropagate[EventSpecifications[evsA, Chronology -> "LeftToRight"], Basis -> basis, OperatorBasis -> opbasis,
                        BackgroundGenerator -> bg1, InitialTimePoint -> t1Start, UseHilbertSpaceWhenPossible -> useHilbQ][rhoOp]];
                sMid = If[ev1Spec["RefocusEvents"] === {}, sMid,
                    NPropagate[EventSpecifications[ev1Spec["RefocusEvents"], Chronology -> "LeftToRight"], Basis -> basis, OperatorBasis -> opbasis,
                        BackgroundGenerator -> bg1, InitialTimePoint -> t1Start + ta,
                        UseHilbertSpaceWhenPossible -> (useHilbQ && !AnyTrue[ev1Spec["RefocusEvents"], filterEventQ])][sMid]];
                {If[evsB === {}, sMid,
                    NPropagate[EventSpecifications[evsB, Chronology -> "LeftToRight"], Basis -> basis, OperatorBasis -> opbasis,
                        BackgroundGenerator -> bg1, InitialTimePoint -> t1Start + ta + durRef,
                        UseHilbertSpaceWhenPossible -> useHilbQ][sMid]],
                 t1Start + tCT + durRef}
        ]
    ];

    Do[
        Module[{key = comp["Key"],
                rhoByG = Lookup[comp, "RhoByGamma", {comp[If[space === "Hilbert", "RhoStates", "RhoVecs"]]}],
                rhoList, nStates, sMat, mixChrons, rSt, rMid, t1End, rMix, tMixEnd, t2Start,
                g2, traj, vals2, mod = Lookup[comp, "RowModulators", None], tGStart, t1Start},
            mixChrons = comp["DirectMixEvents"];
            sMat = ConstantArray[0. + 0. I, {n1, n2}];
            Do[
                tGStart = t0 + gammaShifts[[gIdx]];
                t1Start = tGStart + tauP;
                rhoList = If[Length[rhoByG] === 1, rhoByG[[1]], rhoByG[[gIdx]]];
                nStates = Length[rhoList];
                Do[
                    rSt = stateToOp[rhoList[[Mod[k1 - 1, nStates] + 1]]];
                    {rMid, t1End} = eval1StateAtT1[rSt, t1Vec[[k1]], t1Start];
                    rMix = If[mixChrons === {}, rMid,
                        With[{mixUseHilb = useHilbQ && !AnyTrue[mixChrons, filterEventQ]},
                            NPropagate[EventSpecifications[mixChrons, Chronology -> "LeftToRight"], Basis -> basis, OperatorBasis -> opbasis,
                                BackgroundGenerator -> bgMix, InitialTimePoint -> t1End, UseHilbertSpaceWhenPossible -> mixUseHilb][rMid]
                        ]
                    ];
                    tMixEnd = t1End + tauM;
                    t2Start = tMixEnd + tmin2;
                    If[MemberQ[{"Static", "Periodic"}, ev2Spec["Kind"]],
                        g2 = ev2Spec["Gen2"];
                        traj = Trajectory[rMix -> qObs, {If[nullGenQ[g2], Delay[tmax2 - tmin2], {g2, tmax2 - tmin2}]},
                            Basis -> basis, OperatorBasis -> opbasis, BackgroundGenerator -> bg2,
                            InitialTimePoint -> t2Start, NormalizationFactor -> effNorm,
                            UseHilbertSpaceWhenPossible -> useHilbQ];
                        vals2 = traj /@ (t2Start + (t2Vec - tmin2)),
                        Module[{evs2 = If[ev2Spec["Kind"] === "Stepped", ev2Spec["StepEvents"], {{ev2Spec["DirectFn"], tmax2 - tmin2}}],
                                curState = rMix, stepEffNorm = spec["NormFactor"] * qNormSq, dStep},
                            If[ev2Spec["Kind"] === "Stepped",
                                vals2 = ConstantArray[0. + 0. I, n2];
                                dStep = stageEventDuration[evs2];
                                If[tmin2 != 0., curState = NPropagate[{Delay[tmin2]}, Basis -> basis, OperatorBasis -> opbasis,
                                    BackgroundGenerator -> bg2, InitialTimePoint -> tMixEnd, UseHilbertSpaceWhenPossible -> useHilbQ][curState]];
                                Do[
                                    vals2[[k2]] = TransformationAmplitude[curState -> qObs, {}, Basis -> basis, OperatorBasis -> opbasis,
                                        NormalizationFactor -> stepEffNorm];
                                    If[k2 < n2,
                                        curState = NPropagate[evs2, Basis -> basis, OperatorBasis -> opbasis, BackgroundGenerator -> bg2,
                                            InitialTimePoint -> (t2Start + (k2 - 1) * dStep), UseHilbertSpaceWhenPossible -> useHilbQ][curState]
                                    ],
                                    {k2, n2}
                                ],
                                traj = Trajectory[rMix -> qObs, evs2, Basis -> basis, OperatorBasis -> opbasis,
                                    BackgroundGenerator -> bg2, InitialTimePoint -> t2Start, NormalizationFactor -> effNorm,
                                    UseHilbertSpaceWhenPossible -> useHilbQ];
                                vals2 = traj /@ (t2Start + (t2Vec - tmin2))
                            ]
                        ]
                    ];
                    sMat[[k1]] = sMat[[k1]] + vals2,
                    {k1, n1}
                ],
                {gIdx, nG}
            ];
            sMat = toPackedComplexMat[(1. / nG) * sMat];
            If[mod =!= None, sMat = toPackedComplexMat[toPackedComplexVec[mod] * sMat]];
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
(* Phase 3 Operator-Atom Decomposer & Compiled Ensemble Acceleration      *)
(* ====================================================================== *)

genMatDirect2D[expr_, space_String, basis_, opbasis_, dim_Integer] := If[
    space === "Hilbert",
    generatorToHilbertMatrix[expr, basis, dim],
    generatorToLiouvilleMatrix[expr, opbasis, basis, dim]
];

decomposeSymbolicMatrixEntries2D[matSym_?MatrixQ, dim_Integer] := Module[
    {expMat = Expand[matSym], ht, nz, r, c, v, keys, pairs},
    ht = CreateDataStructure["HashTable"];
    nz = ArrayRules[SparseArray[expMat]];
    Do[
        If[MatchQ[rule, {_Integer, _Integer} -> _],
            {r, c} = rule[[1]];
            v = Expand[rule[[2]]];
            If[v =!= 0 && v =!= 0.,
                If[ht["KeyExistsQ", v],
                    ht["Insert", v -> Append[ht["Lookup", v], {r, c} -> 1.]],
                    ht["Insert", v -> {{r, c} -> 1.}]
                ]
            ]
        ],
        {rule, nz}
    ];
    keys = ht["Keys"];
    Table[{k, toPackedComplexMat[Normal[SparseArray[ht["Lookup", k], {dim, dim}]]]}, {k, keys}]
];

splitGenTerms2D[expr_, symPat_, space_String, basis_, opbasis_, dim_Integer] := Module[
    {facs, constOps, symRest, scExpr, opExpr, m},
    Which[
        nullGenQ[expr], {},
        MatchQ[expr, _Plus],
            Join @@ (splitGenTerms2D[#, symPat, space, basis, opbasis, dim] & /@ (List @@ expr)),
        space === "Hilbert" && MatchQ[expr, HoldPattern[Superoperator[_CommutationSuperoperator, ___]]],
            splitGenTerms2D[First[expr], symPat, space, basis, opbasis, dim],
        space === "Hilbert" && MatchQ[expr, _CommutationSuperoperator],
            ({I * #[[1]], #[[2]]} &) /@ splitGenTerms2D[First[expr], symPat, space, basis, opbasis, dim],
        MatchQ[expr, _Times],
            facs = List @@ expr;
            constOps = Select[facs, (OperatorQ[#] || SuperoperatorQ[#]) && FreeQ[#, symPat] &];
            symRest = Select[facs, !((OperatorQ[#] || SuperoperatorQ[#]) && FreeQ[#, symPat]) &];
            If[Length[constOps] === 1 && FreeQ[Times @@ symRest, _?OperatorQ | _?SuperoperatorQ],
                scExpr = Times @@ symRest;
                m = Quiet[genMatDirect2D[First[constOps], space, basis, opbasis, dim]];
                If[m === $Failed, $Failed, {{scExpr, m}}],
                (* Fallback: compute symbolic matrix representation once *)
                Module[{matSym},
                    matSym = Quiet[
                        If[space === "Hilbert",
                            Normal[MatrixRepresentation[expr, basis]],
                            With[{supExpr = If[OperatorQ[expr] && !SuperoperatorQ[expr], CommutationSuperoperator[-I * expr], expr]},
                                Normal[SuperoperatorMatrixRepresentation[supExpr, opbasis, basis]]
                            ]
                        ]
                    ];
                    If[!MatrixQ[matSym], Return[$Failed]];
                    If[MatrixQ[matSym, NumericQ],
                        {{1., toPackedComplexMat[matSym]}},
                        decomposeSymbolicMatrixEntries2D[matSym, dim]
                    ]
                ]
            ],
        (OperatorQ[expr] || SuperoperatorQ[expr]) && FreeQ[expr, symPat],
            m = Quiet[genMatDirect2D[expr, space, basis, opbasis, dim]];
            If[m === $Failed, $Failed, {{1., m}}],
        True,
            Module[{matSym},
                matSym = Quiet[
                    If[space === "Hilbert",
                        Normal[MatrixRepresentation[expr, basis]],
                        With[{supExpr = If[OperatorQ[expr] && !SuperoperatorQ[expr], CommutationSuperoperator[-I * expr], expr]},
                            Normal[SuperoperatorMatrixRepresentation[supExpr, opbasis, basis]]
                        ]
                    ]
                ];
                If[!MatrixQ[matSym], Return[$Failed]];
                If[MatrixQ[matSym, NumericQ],
                    {{1., toPackedComplexMat[matSym]}},
                    decomposeSymbolicMatrixEntries2D[matSym, dim]
                ]
            ]
    ]
];

(* Decomposes static or time-dependent generator into packed constant matrices + pure Slot Function *)
decomposeGenToAtoms2D[expr_, flatSyms_List, xTest_List, tSym_, space_String, basis_, opbasis_, dim_Integer] := Module[
    {allSyms, symPat, xAllTest, rawExpr, terms, ht, keys, coeffs, mats, keep,
     coeffRules, coeffBody, coeffFn, testVals, mTest, exprTest, mRef, err},
    allSyms = If[tSym === None, flatSyms, Append[flatSyms, tSym]];
    symPat = If[allSyms === {}, Alternatives[], Alternatives @@ allSyms];
    rawExpr = If[tSym === None,
        expr,
        expr /. {
            PeriodicFunction[ts_, _, b_] :> (b /. ts -> tSym),
            TimeDependentGeneratorFunction[ts_, _, _, b_, ___Rule] :> (b /. ts -> tSym),
            f_Function :> f[tSym]
        }
    ];
    If[nullGenQ[rawExpr],
        Return[<|
            "IsZero" -> True, "IsConstant" -> True, "Dim" -> dim,
            "Matrices" -> {}, "CoeffExprs" -> {},
            "CoeffFn" -> If[tSym === None, Function[{}], Function[{}]]
        |>]
    ];
    If[allSyms === {} || FreeQ[rawExpr, symPat],
        With[{m0 = Quiet[genMatDirect2D[rawExpr, space, basis, opbasis, dim]]},
            If[m0 === $Failed || !MatrixQ[m0, NumericQ], Return[$Failed]];
            If[Max[Abs[Flatten[m0]]] == 0.,
                Return[<|"IsZero" -> True, "IsConstant" -> True, "Dim" -> dim,
                         "Matrices" -> {}, "CoeffExprs" -> {}, "CoeffFn" -> Function[{}]|>]
            ];
            Return[<|
                "IsZero" -> False, "IsConstant" -> True, "Dim" -> dim,
                "Matrices" -> {toPackedComplexMat[m0]}, "CoeffExprs" -> {1. + 0. I},
                "CoeffFn" -> Function[{1. + 0. I}]
            |>]
        ]
    ];

    terms = Quiet[splitGenTerms2D[rawExpr, symPat, space, basis, opbasis, dim]];
    If[terms === $Failed || !FreeQ[terms, $Failed], Return[$Failed]];
    If[terms === {},
        Return[<|"IsZero" -> True, "IsConstant" -> True, "Dim" -> dim,
                 "Matrices" -> {}, "CoeffExprs" -> {}, "CoeffFn" -> Function[{}]|>]
    ];

    ht = CreateDataStructure["HashTable"];
    Do[
        With[{sc = Expand[tm[[1]]], mat = toPackedComplexMat[tm[[2]]]},
            If[!FreeQ[sc, _?OperatorQ | _?SuperoperatorQ], Return[$Failed, Module]];
            If[ht["KeyExistsQ", sc],
                ht["Insert", sc -> toPackedComplexMat[ht["Lookup", sc] + mat]],
                ht["Insert", sc -> mat]
            ]
        ],
        {tm, terms}
    ];
    keys = ht["Keys"];
    mats = ht["Lookup", #] & /@ keys;
    keep = Pick[Range[Length[mats]], (Max[Abs[Flatten[#]]] > 0. &) /@ mats];
    If[keep === {},
        Return[<|"IsZero" -> True, "IsConstant" -> True, "Dim" -> dim,
                 "Matrices" -> {}, "CoeffExprs" -> {}, "CoeffFn" -> Function[{}]|>]
    ];
    coeffs = keys[[keep]];
    mats = toPackedComplexMat /@ mats[[keep]];
    If[space === "Hilbert",
        (* Check if rawExpr at xTest is Hermitian? Leave general for safety unless needed; preserve exact matrix! *)
        Null
    ];

    (* Build anonymous pure System Function with Slot[1] (theta vector) and Slot[2] (time) using inert partAt until Function wraps it *)
    Module[{partAt},
        coeffRules = If[flatSyms === {}, {},
            Table[flatSyms[[i]] -> partAt[Slot[1], i], {i, Length[flatSyms]}]
        ];
        If[tSym =!= None, AppendTo[coeffRules, tSym -> Slot[2]]];
        coeffBody = Quiet[Evaluate[N[Expand[coeffs /. coeffRules]] + 0. I]];
        coeffFn = With[{cb = coeffBody}, Function @@ (Hold[cb] /. partAt -> Part)]
    ];

    (* Exact numerical compile-time verification at test point *)
    xAllTest = If[tSym === None, xTest, Append[xTest, 0.173]];
    testVals = Quiet[Check[
        If[tSym === None, coeffFn[xTest], coeffFn[xTest, 0.173]],
        $Failed
    ]];
    If[!VectorQ[testVals, NumericQ] || Length[testVals] =!= Length[mats], Return[$Failed]];
    mTest = toPackedComplexMat[Sum[testVals[[r]] * mats[[r]], {r, Length[mats]}]];

    exprTest = rawExpr;
    If[flatSyms =!= {}, exprTest = exprTest /. Thread[flatSyms -> xTest]];
    If[tSym =!= None, exprTest = exprTest /. tSym -> 0.173];
    mRef = Quiet[genMatDirect2D[exprTest, space, basis, opbasis, dim]];
    If[!MatrixQ[mRef, NumericQ] || Dimensions[mRef] =!= {dim, dim}, Return[$Failed]];
    err = Max[Abs[Flatten[mTest - mRef]]];
    If[err > 1.*^-8 * Max[1., Max[Abs[Flatten[mRef]]]], Return[$Failed]];

    <|
        "IsZero" -> False,
        "IsConstant" -> (flatSyms === {} && tSym === None),
        "Dim" -> dim,
        "Matrices" -> mats,
        "CoeffExprs" -> coeffs,
        "CoeffFn" -> coeffFn
    |>
];

mergeTwoDecomps2D[dA_Association, dB_Association, dim_Integer, hasTimeQ_?BooleanQ] := Module[
    {matsA, matsB, fnA, fnB, allMats, rA, rB},
    Which[
        TrueQ[dA["IsZero"]] && TrueQ[dB["IsZero"]],
            <|"IsZero" -> True, "Dim" -> dim, "Matrices" -> {toPackedComplexMat[ConstantArray[0. + 0. I, {dim, dim}]]},
              "CoeffFn" -> Function[{{0. + 0. I}}]|>,
        TrueQ[dA["IsZero"]], dB,
        TrueQ[dB["IsZero"]], dA,
        True,
            matsA = dA["Matrices"]; matsB = dB["Matrices"];
            rA = Length[matsA]; rB = Length[matsB];
            fnA = dA["CoeffFn"]; fnB = dB["CoeffFn"];
            allMats = Join[matsA, matsB];
            <|
                "IsZero" -> False,
                "Dim" -> dim,
                "Matrices" -> allMats,
                "CoeffFn" -> If[hasTimeQ,
                    With[{f1 = fnA, f2 = fnB}, Function[Join[f1[Slot[1], Slot[2]], f2[Slot[1], Slot[2]]]]],
                    With[{f1 = fnA, f2 = fnB}, Function[Join[f1[Slot[1]], f2[Slot[1]]]]]
                ]
            |>
    ]
];

findRationalDenom2D[ratio_?NumericQ, maxD_Integer: 64] := Module[{d},
    If[Abs[ratio - Round[ratio]] < 1.*^-9, Return[1]];
    For[d = 2, d <= maxD, d++,
        If[Abs[ratio * d - Round[ratio * d]] < 1.*^-9, Return[d]]
    ];
    None
];

resolve2DPeriod[periodOpt_, ev1Spec_Association, ev2Spec_Association, bgList_List] := Module[
    {p1, p2, pb2, pb4, comb},
    If[NumericQ[periodOpt] && periodOpt > 0, Return[N[Re[periodOpt]]]];
    If[!MemberQ[{"Static", "Periodic"}, ev1Spec["Kind"]] || !MemberQ[{"Static", "Periodic"}, ev2Spec["Kind"]],
        Return[None]
    ];
    p1 = If[ev1Spec["Kind"] === "Periodic", ev1Spec["Period"], All];
    p2 = If[ev2Spec["Kind"] === "Periodic", ev2Spec["Period"], All];
    pb2 = extractGeneratorPeriod[bgList[[2]]];
    pb4 = extractGeneratorPeriod[bgList[[4]]];
    comb = CommonPeriod[p1, pb2, p2, pb4];
    Which[
        NumericQ[comb] && comb > 0, N[Re[comb]],
        comb === All, All,
        True, None
    ]
];

buildPackedSpinZMat[qSpins_List, phi_?NumericQ, space_String, basis_, opbasis_, dim_Integer] := If[
    Chop[N[phi]] == 0.,
    toPackedComplexMat[IdentityMatrix[dim]],
    If[space === "Hilbert",
        toPackedComplexMat[MatrixRepresentation[RotationOperator[qSpins, {phi, "z"}], basis]],
        toPackedComplexMat[SuperoperatorMatrixRepresentation[RotationSuperoperator[qSpins, {phi, "z"}], opbasis, basis]]
    ]
];

buildPackedFilterOp[spins_List, orders_List, space_String, system_List, basis_, opbasis_] := If[
    space === "Hilbert",
    <|"Kind" -> "Filter", "MaskSpec" -> BuildHilbertCoherenceMask[{{spins, orders}}, system, basis]|>,
    <|"Kind" -> "Filter", "Mat" -> toPackedComplexMat[SuperoperatorMatrixRepresentation[CoherenceOrderFiltrationSuperoperator[spins, orders], opbasis, basis]]|>
];

(* Compile ensemble problem into pure System`/Developer` chunk payload when possible *)
tryCompileEnsemble2D[heldAll_HoldComplete, stdEns_, rules_List, weights_List, parVal_, maxBatch_] := Module[
    {posArgsHold, optArgsHold, origVars, nEns = Length[rules], flatSyms, varSubRules, thetaList,
     nonStagePat, optRulesRaw, timing1, timing2, ev1Opt, mixOpt, ev2Opt, prepOpt, rho0Opt, obsOpt, bgOpt,
     cpOpt, pcOpt, lbOpt, scmOpt, periodOpt, t0Val, caravOpt, magSubOpt, normVal, useHilbOpt, forceEvenQ, chron,
     system, basis, opbasis, dimH, dimL, dim, qm, qOrder, qSpins, obsRaw, obsOp, obsSpins, cpPair, p1Filt, pMixFilt,
     t1Parsed, t2Parsed, n1Log, t1TimingForKernel, t1TimeVec, n1Kernel, bgList, ev1Spec, ev2Spec, resolvedPeriod,
     isPeriodicProblem, caravActiveQ, scmChosen, space, symPat, prepChron, mixChronsBase, tauPrep, tauMix,
     qAdjMat, qVec, qNormSq, normDenom, qObj, tSym, xTest, mode, keys, rowMod, rotNegHalf, rotPi, rotTppi,
     baseMixPart, templateComps, quadMeta, rho0VecOrMat, baseStates,
     extractStageOps, prepStageOps = None, mixStageOps = None, buildSimpleMixStep,
     prepFixed = True, mixFixed = True, decCache, decompStage, atomSpecs = <||>, fixedSpec = <||>,
     nDiv = 1, nSubPerDiv = 4, microSub = 2, nGamma = 1, gammaOffsets = {0}, gammaShifts = {0.},
     uniqueMixPhases = {0.}, chunkSize, nK, chunks, payloads, resOutList, accFinal, firstAssoc, resSig,
     checkFreeQ, buildCoeffChunk, decomp1, decompBg1, decomp2, decompBg2, decompBgPrep, decompBgMix,
     decomp1a, decomp1b, t1Min, t2Min, dt1, dt2, q1, q2, nBase, evalStageStaticQ},

    posArgsHold = Select[heldAll, !MatchQ[Unevaluated[#], _Rule | _RuleDelayed] &];
    If[Length[posArgsHold] < 1 || Length[posArgsHold] > 2, Return[$Failed]];
    optArgsHold = Select[heldAll, MatchQ[Unevaluated[#], _Rule | _RuleDelayed] &];

    origVars = First[stdEns];
    checkFreeQ[expr_] := FreeQ[expr, Alternatives @@ origVars];

    (* Flatten scalar and multidimensional ensemble variables (e.g. Omega -> {a, b, c}) *)
    With[{blockSyms = Select[origVars, MatchQ[#, _Symbol] &]},
      Apply[Block, Prepend[Hold[
        (* Map each variable rule value dimensions based on first sample *)
        Module[{firstRule = rules[[1]], varNames = Keys[rules[[1]]]},
            varSubRules = {};
            flatSyms = {};
            Do[
                With[{vn = varNames[[k]], val1 = First[Values[firstRule[[k ;; k]]]]},
                    If[ListQ[val1],
                        With[{ss = Table[Unique["sdEns$"], {Length[val1]}]},
                            AppendTo[varSubRules, vn -> ss];
                            flatSyms = Join[flatSyms, ss]
                        ],
                        With[{s = Unique["sdEns$"]},
                            AppendTo[varSubRules, vn -> s];
                            AppendTo[flatSyms, s]
                        ]
                    ]
                ],
                {k, Length[varNames]}
            ];
            thetaList = Developer`ToPackedArray[N[Flatten[Values[#]] & /@ rules]];
            xTest = Normal[thetaList[[1]]];

            (* Extract options while protecting origVars from global OwnValues *)
            optRulesRaw = Apply[List, DeleteCases[optArgsHold, (Rule | RuleDelayed)[EnsembleAverage | Parallel | MaxBatchSize, _], {1}]];
            timing1 = posArgsHold[[1]];
            If[Length[posArgsHold] === 2,
                If[ListQ[posArgsHold[[2]]] && Length[posArgsHold[[2]]] === 3 &&
                   !(NumericQ[posArgsHold[[2, 1]]] && NumericQ[posArgsHold[[2, 2]]]),
                    Return[$Failed],
                    timing2 = posArgsHold[[2]];
                    If[!AllTrue[{timing1, timing2}, checkFreeQ], Return[$Failed]]
                ],
                If[ListQ[timing1] && Length[timing1] === 2,
                    {timing1, timing2} = timing1;
                    If[!AllTrue[{timing1, timing2}, checkFreeQ], Return[$Failed]],
                    Return[$Failed]
                ]
            ];

            chron = Chronology /. optRulesRaw /. Options[Signal2D];
            forceEvenQ = TrueQ[ForceEvenNumberOfPoints /. optRulesRaw /. Options[Signal2D]];
            qm = QuadratureMethod /. optRulesRaw /. Options[Signal2D];
            If[qm === None, qm = "PhaseModulated"];
            qOrder = Abs[N[Re[QuadratureOrder /. optRulesRaw /. Options[Signal2D]]]];
            qSpins = QuadratureSpins /. optRulesRaw /. Options[Signal2D];
            basis = Basis /. optRulesRaw /. Options[Signal2D];
            opbasis = OperatorBasis /. optRulesRaw /. Options[Signal2D];
            obsRaw = Observable /. optRulesRaw /. Options[Signal2D];
            rho0Opt = InitialDensityOperator /. optRulesRaw /. Options[Signal2D];
            prepOpt = Preparation /. optRulesRaw /. Options[Signal2D];
            mixOpt = Mixing /. optRulesRaw /. Options[Signal2D];
            ev1Opt = Evolution1 /. optRulesRaw /. Options[Signal2D];
            ev2Opt = Evolution2 /. optRulesRaw /. Options[Signal2D];
            bgOpt = BackgroundGenerator /. optRulesRaw /. Options[Signal2D];
            cpOpt = CoherencePathway /. optRulesRaw /. Options[Signal2D];
            pcOpt = PhaseCycle /. optRulesRaw /. Options[Signal2D];
            lbOpt = LineBroadening /. optRulesRaw /. Options[Signal2D];
            scmOpt = SignalCalculationMethod /. optRulesRaw /. Options[Signal2D];
            periodOpt = SignalGeneratorPeriod /. optRulesRaw /. Options[Signal2D];
            t0Val = InitialTimePoint /. optRulesRaw /. Options[Signal2D];
            caravOpt = CarouselAverage /. optRulesRaw /. Options[Signal2D];
            magSubOpt = MagnusSubsteps /. optRulesRaw /. Options[Signal2D];
            normVal = NormalizationFactor /. optRulesRaw /. Options[Signal2D];
            useHilbOpt = TrueQ[UseHilbertSpaceWhenPossible /. optRulesRaw /. Options[Signal2D]];

            If[!AllTrue[{chron, qm, qOrder, qSpins, basis, opbasis, obsRaw, rho0Opt, prepOpt, mixOpt, cpOpt, pcOpt, lbOpt, scmOpt, periodOpt, t0Val, caravOpt, normVal}, checkFreeQ],
                Return[$Failed]
            ];
            If[pcOpt =!= None && pcOpt =!= {}, Return[$Failed]];
            If[AssociationQ[prepOpt], Return[$Failed]];
            If[!MemberQ[{"LeftToRight", "RightToLeft"}, chron], Return[$Failed]];
            If[!MemberQ[{"States", "StatesTPPI", "EchoAntiEcho", "TPPI", "PhaseModulated"}, qm], Return[$Failed]];

            t1Parsed = parseTimingAxis[timing1, forceEvenQ];
            t2Parsed = parseTimingAxis[timing2, forceEvenQ];
            If[t1Parsed === $Failed || t2Parsed === $Failed, Return[$Failed]];
            system = SpinSystem[];
            If[!MatrixQ[system], Return[$Failed]];
            dimH = Length[BasisKets[basis]]; dimL = Length[BasisOperators[opbasis]];
            qSpins = resolveQuadratureSpins[qSpins, system];
            obsOp = resolveObservable[obsRaw, system];
            If[qSpins === $Failed || obsOp === $Failed, Return[$Failed]];
            obsSpins = Which[obsRaw === All, First /@ system, ValidSpinQ[obsRaw, system], {obsRaw}, ListQ[obsRaw] && ValidSpinQ[obsRaw, system], obsRaw, True, First /@ system];
            cpPair = resolveCoherencePathway[cpOpt, qSpins, obsSpins];
            If[cpPair === $Failed, Return[$Failed]];
            {p1Filt, pMixFilt} = cpPair;

            normVal = Which[normVal === None || normVal === False, 1., NumericQ[normVal] && Chop[N[Re[normVal]]] != 0., N[Re[normVal]], True, Return[$Failed]];
            t0Val = If[NumericQ[t0Val], N[Re[t0Val]], Return[$Failed]];

            bgList = Which[
                ListQ[bgOpt] && Length[bgOpt] === 4, If[chron === "RightToLeft", Reverse[bgOpt], bgOpt],
                ListQ[bgOpt], Return[$Failed],
                True, {bgOpt, bgOpt, bgOpt, bgOpt}
            ];

            (* Substitute flatSyms BEFORE checking OperatorQ/PeriodicFunction so WignerD etc. expand cleanly *)
            {ev1Opt, ev2Opt, bgList} = {ev1Opt, ev2Opt, bgList} /. varSubRules;
            ev1Spec = parseEvolution1[ev1Opt, chron];
            ev2Spec = parseEvolution2[ev2Opt, chron];
            If[ev1Spec === $Failed || ev2Spec === $Failed, Return[$Failed]];
            If[!MemberQ[{"Static", "StaticSplit", "Periodic"}, ev1Spec["Kind"]] || !MemberQ[{"Static", "Periodic"}, ev2Spec["Kind"]],
                Return[$Failed]
            ];

            (* Determine space (evaluate at test sample so Hilbert compatibility matches non-ensemble runs) *)
            With[{testSub = Thread[flatSyms -> xTest]},
                space = If[useHilbOpt &&
                    OperatorQ[rho0Opt] &&
                    AllTrue[bgList /. testSub, hilbertGenQ] &&
                    Which[
                        ev1Spec["Kind"] === "Static" || ev1Spec["Kind"] === "Periodic", hilbertGenQ[ev1Spec["Gen1"] /. testSub],
                        ev1Spec["Kind"] === "StaticSplit", hilbertGenQ[ev1Spec["Gen1a"] /. testSub] && hilbertGenQ[ev1Spec["Gen1b"] /. testSub] && hilbertEventListQ[ev1Spec["RefocusEvents"]],
                        True, False
                    ] &&
                    hilbertGenQ[ev2Spec["Gen2"] /. testSub] &&
                    hilbertStageQ[prepOpt, chron] && hilbertStageQ[mixOpt, chron],
                    "Hilbert", "Liouville"
                ]
            ];
            dim = If[space === "Hilbert", dimH, dimL];

            (* Check periodic / COMPUTE2D vs DoubleDiag *)
            resolvedPeriod = resolve2DPeriod[periodOpt, ev1Spec /. Thread[flatSyms -> xTest], ev2Spec /. Thread[flatSyms -> xTest], bgList /. Thread[flatSyms -> xTest]];
            isPeriodicProblem = NumericQ[resolvedPeriod] && resolvedPeriod > 0 && MemberQ[{"Static", "Periodic"}, ev1Spec["Kind"]] && MemberQ[{"Static", "Periodic"}, ev2Spec["Kind"]];
            caravActiveQ = (caravOpt === True) || (IntegerQ[caravOpt] && caravOpt > 1);
            If[scmOpt === "Direct" || scmOpt === "SteppedPropagator", Return[$Failed]];
            If[caravActiveQ && (!isPeriodicProblem || scmOpt === "DoubleDiagonalization"), Return[$Failed]];
            scmChosen = Which[
                scmOpt === "COMPUTE2D", If[isPeriodicProblem, "COMPUTE2D", Return[$Failed]],
                scmOpt === "DoubleDiagonalization", If[!isPeriodicProblem && MemberQ[{"Static", "StaticSplit"}, ev1Spec["Kind"]], "DoubleDiagonalization", Return[$Failed]],
                scmOpt === Automatic, If[isPeriodicProblem, "COMPUTE2D", "DoubleDiagonalization"],
                True, Return[$Failed]
            ];
            mode = If[scmChosen === "COMPUTE2D", "COMPUTE2D", "DoubleDiag"];

            n1Log = t1Parsed["N"];
            If[qm === "TPPI",
                With[{t1T = t1Parsed["Timing"]},
                    t1TimingForKernel = {t1T[[1]], t1T[[1]] + (2 * n1Log - 1) * (0.5 * t1T[[3]]), 0.5 * t1T[[3]]};
                    t1TimeVec = toPackedRealVec[t1TimingForKernel[[1]] + t1TimingForKernel[[3]] * Range[0, 2 * n1Log - 1]];
                    n1Kernel = 2 * n1Log
                ],
                t1TimingForKernel = t1Parsed["Timing"];
                t1TimeVec = t1Parsed["TimeVec"];
                n1Kernel = n1Log
            ];

            If[space === "Hilbert",
                With[{qMat = toPackedComplexMat[MatrixRepresentation[obsOp, basis]]},
                    qNormSq = Total[Abs[Flatten[qMat]]^2];
                    qAdjMat = toPackedComplexMat[ConjugateTranspose[qMat]];
                    qObj = qAdjMat;
                    normDenom = normVal * qNormSq
                ],
                qVec = toPackedComplexVec[OperatorVectorRepresentation[obsOp, opbasis, basis]];
                qNormSq = Total[Abs[qVec]^2];
                qObj = qVec;
                normDenom = normVal * qNormSq
            ];
            If[normDenom == 0., Return[$Failed]];

            prepChron = normalizeEventsChronological[prepOpt, chron];
            mixChronsBase = normalizeEventsChronological[mixOpt, chron];
            If[prepChron === $Failed || mixChronsBase === $Failed, Return[$Failed]];
            If[pMixFilt =!= None, AppendTo[mixChronsBase, CoherenceOrderFiltrationSuperoperator[obsSpins, pMixFilt]]];
            tauPrep = stageEventDuration[prepChron];
            tauMix = stageEventDuration[mixChronsBase];

            (* Bail out cleanly if static DoubleDiag has ensemble-dependent Prep/Mix finite propagation *)
            symPat = If[flatSyms === {}, Alternatives[], Alternatives @@ flatSyms];
            If[mode === "DoubleDiag" && ((tauPrep > 0. && !FreeQ[bgList[[1]], symPat]) || (tauMix > 0. && !FreeQ[bgList[[3]], symPat])),
                Return[$Failed]
            ];

            (* Helper: decompose static expression & evaluate all member coefficients in one vectorized call *)
            buildCoeffChunk[dec_Association] := If[
                TrueQ[dec["IsConstant"]],
                ConstantArray[N[dec["CoeffExprs"]] + 0. I, nEns],
                With[{fn = dec["CoeffFn"]},
                    Developer`ToPackedArray[fn /@ thetaList]
                ]
            ];

            tSym = If[mode === "COMPUTE2D", Unique["sdT$"], None];

            If[mode === "DoubleDiag",
                decomp2 = decomposeGenToAtoms2D[ev2Spec["Gen2"], flatSyms, xTest, None, space, basis, opbasis, dim];
                decompBg2 = decomposeGenToAtoms2D[bgList[[4]], flatSyms, xTest, None, space, basis, opbasis, dim];
                If[decomp2 === $Failed || decompBg2 === $Failed, Return[$Failed]];
                With[{m2Dec = mergeTwoDecomps2D[decomp2, decompBg2, dim, False]},
                    If[!FreeQ[{ev2Spec["Gen2"], bgList[[4]]}, symPat],
                        atomSpecs["Gen2Full"] = <|"IsZero" -> m2Dec["IsZero"], "Dim" -> dim, "Matrices" -> m2Dec["Matrices"]|>;
                        atomSpecs["Gen2CoeffsAll"] = buildCoeffChunk[m2Dec],
                        fixedSpec["Gen2Mat"] = If[m2Dec["IsZero"],
                            toPackedComplexMat[ConstantArray[0. + 0. I, {dim, dim}]],
                            toPackedComplexMat[First[m2Dec["CoeffFn"][xTest]] * 0 + Sum[m2Dec["CoeffFn"][xTest][[r]] * m2Dec["Matrices"][[r]], {r, Length[m2Dec["Matrices"]]}]]
                        ]
                    ]
                ];

                fixedSpec["Ev1Mode"] = ev1Spec["Mode"];
                If[ev1Spec["Mode"] === "Standard",
                    If[(ev1Spec["Gen1"] === ev2Spec["Gen2"]) && (bgList[[2]] === bgList[[4]]),
                        fixedSpec[If[space === "Hilbert", "SameH1H2", "SameL1L2"]] = True,
                        fixedSpec[If[space === "Hilbert", "SameH1H2", "SameL1L2"]] = False;
                        decomp1 = decomposeGenToAtoms2D[ev1Spec["Gen1"], flatSyms, xTest, None, space, basis, opbasis, dim];
                        decompBg1 = decomposeGenToAtoms2D[bgList[[2]], flatSyms, xTest, None, space, basis, opbasis, dim];
                        If[decomp1 === $Failed || decompBg1 === $Failed, Return[$Failed]];
                        With[{m1Dec = mergeTwoDecomps2D[decomp1, decompBg1, dim, False]},
                            If[!FreeQ[{ev1Spec["Gen1"], bgList[[2]]}, symPat],
                                atomSpecs["Gen1Full"] = <|"IsZero" -> m1Dec["IsZero"], "Dim" -> dim, "Matrices" -> m1Dec["Matrices"]|>;
                                atomSpecs["Gen1CoeffsAll"] = buildCoeffChunk[m1Dec],
                                fixedSpec["Gen1Mat"] = If[m1Dec["IsZero"],
                                    toPackedComplexMat[ConstantArray[0. + 0. I, {dim, dim}]],
                                    toPackedComplexMat[Sum[m1Dec["CoeffFn"][xTest][[r]] * m1Dec["Matrices"][[r]], {r, Length[m1Dec["Matrices"]]}]]
                                ]
                            ]
                        ]
                    ],
                    (* StaticSplit Echo / ConstantTime *)
                    fixedSpec["CTDuration"] = Lookup[ev1Spec, "CTDuration", 0.];
                    decomp1a = decomposeGenToAtoms2D[ev1Spec["Gen1a"], flatSyms, xTest, None, space, basis, opbasis, dim];
                    decomp1b = decomposeGenToAtoms2D[ev1Spec["Gen1b"], flatSyms, xTest, None, space, basis, opbasis, dim];
                    decompBg1 = decomposeGenToAtoms2D[bgList[[2]], flatSyms, xTest, None, space, basis, opbasis, dim];
                    If[decomp1a === $Failed || decomp1b === $Failed || decompBg1 === $Failed, Return[$Failed]];
                    If[tauPrep == 0. && !FreeQ[bgList[[2]], symPat] && stageEventDuration[ev1Spec["RefocusEvents"]] > 0.,
                        Return[$Failed]
                    ];
                    With[{m1a = mergeTwoDecomps2D[decomp1a, decompBg1, dim, False], m1b = mergeTwoDecomps2D[decomp1b, decompBg1, dim, False]},
                        atomSpecs["Gen1aFull"] = <|"IsZero" -> m1a["IsZero"], "Dim" -> dim, "Matrices" -> m1a["Matrices"]|>;
                        atomSpecs["Gen1aCoeffsAll"] = buildCoeffChunk[m1a];
                        atomSpecs["Gen1bFull"] = <|"IsZero" -> m1b["IsZero"], "Dim" -> dim, "Matrices" -> m1b["Matrices"]|>;
                        atomSpecs["Gen1bCoeffsAll"] = buildCoeffChunk[m1b]
                    ];
                    If[space === "Hilbert",
                        With[{refS = buildHilbertStageSteps[ev1Spec["RefocusEvents"], None, basis, system]},
                            If[refS === $Failed, Return[$Failed]];
                            With[{refFin = If[refS === {}, {{"Unitary", toPackedComplexMat[IdentityMatrix[dim]], toPackedComplexMat[IdentityMatrix[dim]]}}, refS]},
                                fixedSpec["RefocusSteps"] = refFin;
                                fixedSpec["RefocusSingleUnitary"] = If[Length[refFin] === 1 && refFin[[1, 1]] === "Unitary", refFin[[1, 2]], None]
                            ]
                        ],
                        fixedSpec["RefocusMat"] = buildLiouvilleStageMatrix[ev1Spec["RefocusEvents"], None, opbasis, basis, dim]
                    ]
                ];

                (* Align key names with EvaluateCompiledEnsembleChunk2D *)
                If[space === "Hilbert",
                    If[KeyExistsQ[atomSpecs, "Gen2Full"],
                        atomSpecs["H2"] = atomSpecs["Gen2Full"]; atomSpecs["H2CoeffsAll"] = atomSpecs["Gen2CoeffsAll"]
                    ];
                    If[KeyExistsQ[fixedSpec, "Gen2Mat"], fixedSpec["H2Mat"] = fixedSpec["Gen2Mat"]];
                    If[KeyExistsQ[atomSpecs, "Gen1Full"],
                        atomSpecs["H1"] = atomSpecs["Gen1Full"]; atomSpecs["H1CoeffsAll"] = atomSpecs["Gen1CoeffsAll"]
                    ];
                    If[KeyExistsQ[fixedSpec, "Gen1Mat"], fixedSpec["H1Mat"] = fixedSpec["Gen1Mat"]];
                    If[KeyExistsQ[atomSpecs, "Gen1aFull"],
                        atomSpecs["H1a"] = atomSpecs["Gen1aFull"]; atomSpecs["H1aCoeffsAll"] = atomSpecs["Gen1aCoeffsAll"];
                        atomSpecs["H1b"] = atomSpecs["Gen1bFull"]; atomSpecs["H1bCoeffsAll"] = atomSpecs["Gen1bCoeffsAll"]
                    ],
                    If[KeyExistsQ[atomSpecs, "Gen2Full"],
                        atomSpecs["L2"] = atomSpecs["Gen2Full"]; atomSpecs["L2CoeffsAll"] = atomSpecs["Gen2CoeffsAll"]
                    ];
                    If[KeyExistsQ[fixedSpec, "Gen2Mat"], fixedSpec["L2Mat"] = fixedSpec["Gen2Mat"]];
                    If[KeyExistsQ[atomSpecs, "Gen1Full"],
                        atomSpecs["L1"] = atomSpecs["Gen1Full"]; atomSpecs["L1CoeffsAll"] = atomSpecs["Gen1CoeffsAll"]
                    ];
                    If[KeyExistsQ[fixedSpec, "Gen1Mat"], fixedSpec["L1Mat"] = fixedSpec["Gen1Mat"]];
                    If[KeyExistsQ[atomSpecs, "Gen1aFull"],
                        atomSpecs["L1a"] = atomSpecs["Gen1aFull"]; atomSpecs["L1aCoeffsAll"] = atomSpecs["Gen1aCoeffsAll"];
                        atomSpecs["L1b"] = atomSpecs["Gen1bFull"]; atomSpecs["L1bCoeffsAll"] = atomSpecs["Gen1bCoeffsAll"]
                    ];
                    fixedSpec["Dt1"] = t1TimingForKernel[[3]]; fixedSpec["Tmin1"] = t1TimingForKernel[[1]];
                    fixedSpec["Dt2"] = t2Parsed["Timing"][[3]]; fixedSpec["Tmin2"] = t2Parsed["Timing"][[1]]
                ],

                (* Periodic COMPUTE2D compiled ensemble setup *)
                q1 = Quiet[With[{r = findRationalDenom2D[t1TimingForKernel[[3]] / resolvedPeriod, 64]}, If[IntegerQ[r], r, 8]]];
                q2 = Quiet[With[{r = findRationalDenom2D[t2Parsed["Timing"][[3]] / resolvedPeriod, 64]}, If[IntegerQ[r], r, 8]]];
                nGamma = Which[caravOpt === False || caravOpt === None || caravOpt === 1, 1, TrueQ[caravOpt], Max[8, LCM[q1, q2]], True, caravOpt];
                nDiv = LCM[q1, q2, If[TrueQ[caravOpt], 8, nGamma]];
                If[nDiv > 256, nDiv = Min[256, LCM[q1, q2]]; If[nDiv > 256, nDiv = Max[8, Min[64, q1, q2]]]];
                nSubPerDiv = If[IntegerQ[magSubOpt] && magSubOpt >= 1, magSubOpt, Max[4, Ceiling[64 / nDiv]]];
                microSub = Max[2, nSubPerDiv];
                gammaOffsets = If[nGamma === 1, {0}, If[Divisible[nDiv, nGamma], Range[0, nDiv - 1, nDiv / nGamma], Table[N[g * nDiv / nGamma], {g, 0, nGamma - 1}]]];
                gammaShifts = N[gammaOffsets * (resolvedPeriod / nDiv)];

                decomp1 = decomposeGenToAtoms2D[ev1Spec["Gen1"], flatSyms, xTest, tSym, space, basis, opbasis, dim];
                decompBg1 = decomposeGenToAtoms2D[bgList[[2]], flatSyms, xTest, tSym, space, basis, opbasis, dim];
                decomp2 = decomposeGenToAtoms2D[ev2Spec["Gen2"], flatSyms, xTest, tSym, space, basis, opbasis, dim];
                decompBg2 = decomposeGenToAtoms2D[bgList[[4]], flatSyms, xTest, tSym, space, basis, opbasis, dim];
                If[MemberQ[{decomp1, decompBg1, decomp2, decompBg2}, $Failed], Return[$Failed]];

                With[{td1 = mergeTwoDecomps2D[decomp1, decompBg1, dim, True],
                      td2 = mergeTwoDecomps2D[decomp2, decompBg2, dim, True],
                      sameQ = (ev1Spec["Gen1"] === ev2Spec["Gen2"]) && (bgList[[2]] === bgList[[4]])},
                    atomSpecs["TD1Spec"] = <|"Matrices" -> td1["Matrices"], "CoeffFn" -> td1["CoeffFn"]|>;
                    atomSpecs["TD2Spec"] = If[sameQ, atomSpecs["TD1Spec"], <|"Matrices" -> td2["Matrices"], "CoeffFn" -> td2["CoeffFn"]|>];
                    fixedSpec["SameEv1Ev2"] = sameQ
                ];
                If[tauPrep > 0. && periodicGenQ[bgList[[1]] /. Thread[flatSyms -> xTest]],
                    decompBgPrep = decomposeGenToAtoms2D[bgList[[1]], flatSyms, xTest, tSym, space, basis, opbasis, dim];
                    If[decompBgPrep === $Failed, Return[$Failed]];
                    atomSpecs["PrepTDSpec"] = <|"Matrices" -> decompBgPrep["Matrices"], "CoeffFn" -> decompBgPrep["CoeffFn"]|>,
                    atomSpecs["PrepTDSpec"] = None
                ];
                If[tauMix > 0. && periodicGenQ[bgList[[3]] /. Thread[flatSyms -> xTest]],
                    decompBgMix = decomposeGenToAtoms2D[bgList[[3]], flatSyms, xTest, tSym, space, basis, opbasis, dim];
                    If[decompBgMix === $Failed, Return[$Failed]];
                    atomSpecs["MixTDSpec"] = <|"Matrices" -> decompBgMix["Matrices"], "CoeffFn" -> decompBgMix["CoeffFn"]|>,
                    atomSpecs["MixTDSpec"] = None
                ];

                fixedSpec["Period"] = resolvedPeriod;
                fixedSpec["NDiv"] = nDiv;
                fixedSpec["T0"] = t0Val;
                fixedSpec["NSubPerDiv"] = nSubPerDiv;
                fixedSpec["MicroSubsteps"] = microSub;
                fixedSpec["Dim"] = dim;
                fixedSpec["PrepDuration"] = tauPrep;
                fixedSpec["MixDuration"] = tauMix;
                fixedSpec["GammaOffsets"] = gammaOffsets;
                fixedSpec["UniqueMixStartPhases"] = DeleteDuplicates[Round[Mod[Flatten[Table[gammaShifts[[g]] + tauPrep + t1TimeVec, {g, nGamma}]], resolvedPeriod], 1.*^-11]]
            ];

            (* Extract packed Preparation & Mixing stage representations *)
            extractStageOps[chEvs_List, bgStatic_] := Module[{ops = {}, tmp},
                Do[
                    Which[
                        filterEventQ[ev],
                            With[{filts = extractFiltersFromEvent[ev]},
                                If[filts === $Failed, Return[$Failed, Module]];
                                Do[AppendTo[ops, buildPackedFilterOp[f[[1]], f[[2]], space, system, basis, opbasis]], {f, filts}]
                            ],
                        MatchQ[ev, _RotationOperator | _RotationSuperoperator | _Pulse] && stageEventDuration[ev] == 0.,
                            With[{u = If[space === "Hilbert",
                                    buildHilbertPropagatorFromSegment[{ev}, None, basis, dimH, 0.],
                                    buildLiouvilleStageMatrix[{ev}, None, opbasis, basis, dimL, 0.]]},
                                If[u === $Failed, Return[$Failed, Module]];
                                AppendTo[ops, <|"Kind" -> "Unitary", "Mat" -> u|>]
                            ],
                        MatchQ[ev, Delay[_?NumericQ]],
                            AppendTo[ops, <|"Kind" -> "DelayCF4", "Dur" -> N[Re[First[ev]]]|>],
                        MatchQ[ev, {_?OperatorQ | _?SuperoperatorQ | None | NullOperator[] | NullSuperoperator[] | 0 | 0., _?NumericQ}],
                            With[{em = Quiet[genMatDirect2D[ev[[1]], space, basis, opbasis, dim]]},
                                If[em === $Failed, Return[$Failed, Module]];
                                If[Max[Abs[Flatten[em]]] == 0.,
                                    AppendTo[ops, <|"Kind" -> "DelayCF4", "Dur" -> N[Re[ev[[2]]]]|>],
                                    AppendTo[ops, <|"Kind" -> "DelayCF4", "Dur" -> N[Re[ev[[2]]]], "EvMat" -> em|>]
                                ]
                            ],
                        True,
                            Return[$Failed, Module]
                    ],
                    {ev, chEvs}
                ];
                ops
            ];

            (* Build base Rho0 / RhoBase and Mixing MixPart *)
            rho0VecOrMat = If[space === "Hilbert",
                toPackedComplexMat[MatrixRepresentation[rho0Opt, basis]],
                toPackedComplexVec[OperatorVectorRepresentation[rho0Opt, opbasis, basis]]
            ];

            If[mode === "COMPUTE2D" && tauPrep > 0. && periodicGenQ[bgList[[1]] /. Thread[flatSyms -> xTest]],
                prepStageOps = extractStageOps[prepChron, bgList[[1]]];
                If[prepStageOps === $Failed, Return[$Failed]];
                fixedSpec["PrepStageOps"] = {prepStageOps};
                baseStates = {rho0VecOrMat},
                With[{rInit = If[space === "Hilbert",
                        With[{st = buildHilbertStageSteps[prepChron, bgList[[1]], basis, system, t0Val]},
                            If[st === $Failed, Return[$Failed]]; ApplyHilbertMixSequence[rho0VecOrMat, st]],
                        With[{uM = buildLiouvilleStageMatrix[prepChron, bgList[[1]], opbasis, basis, dimL, t0Val]},
                            If[uM === $Failed, Return[$Failed]]; toPackedComplexVec[uM . rho0VecOrMat]]
                    ]},
                    baseStates = {If[p1Filt =!= None, applySpinFilter[rInit, space, qSpins, p1Filt, system, basis, opbasis], rInit]}
                ]
            ];

            If[mode === "COMPUTE2D" && tauMix > 0. && periodicGenQ[bgList[[3]] /. Thread[flatSyms -> xTest]],
                mixStageOps = extractStageOps[mixChronsBase, bgList[[3]]];
                If[mixStageOps === $Failed, Return[$Failed]];
                fixedSpec["MixStageOps"] = {mixStageOps};
                baseMixPart = <||>;
                fixedSpec["MixCOMPUTE"] = {<||>},
                (* Phase-independent static/instantaneous mixing step *)
                If[space === "Hilbert",
                    Module[{steps = buildHilbertStageSteps[mixChronsBase, bgList[[3]], basis, system, t0Val], isU, uS},
                        If[steps === $Failed, Return[$Failed]];
                        If[steps === {}, steps = {{"Unitary", toPackedComplexMat[IdentityMatrix[dimH]], toPackedComplexMat[IdentityMatrix[dimH]]}}];
                        isU = (Length[steps] === 1 && steps[[1, 1]] === "Unitary");
                        uS = If[isU, steps[[1, 2]], None];
                        baseMixPart = <|"MixSteps" -> steps, "IsSingleUnitary" -> isU, "SingleUnitary" -> uS|>;
                        fixedSpec["MixCOMPUTE"] = {
                            <|"Static" -> If[isU,
                                If[Max[Abs[uS - IdentityMatrix[dimH]]] <= 1.*^-13, <|"Kind" -> "Identity"|>, <|"Kind" -> "Unitary", "U" -> uS, "UH" -> toPackedComplexMat[ConjugateTranspose[uS]]|>],
                                <|"Kind" -> "Steps", "Steps" -> steps|>
                            ]|>
                        }
                    ],
                    Module[{mM = buildLiouvilleStageMatrix[mixChronsBase, bgList[[3]], opbasis, basis, dimL, t0Val]},
                        If[mM === $Failed, Return[$Failed]];
                        baseMixPart = <|"MixMat" -> mM|>;
                        fixedSpec["MixCOMPUTE"] = {
                            <|"Static" -> If[Max[Abs[mM - IdentityMatrix[dimL]]] <= 1.*^-13, <|"Kind" -> "Identity"|>, <|"Kind" -> "Mat", "Mat" -> mM|>]|>
                        }
                    ]
                ]
            ];

            If[space === "Liouville" && mode === "DoubleDiag",
                fixedSpec["ActQ"] = Pick[Range[dimL], UnitStep[Abs[qVec] - 1.*^-14], 1]
            ];

            (* Build QuadratureMeta templates once *)
            buildStateOpRotZ[phi_] := <|"Kind" -> "RotZ", "Mat" -> buildPackedSpinZMat[qSpins, phi, space, basis, opbasis, dim]|>;
            buildStateOpFilter[ord_] := buildPackedFilterOp[qSpins, {ord}, space, system, basis, opbasis];
            With[{idOp = <|"Kind" -> "Identity"|>},
                templateComps = Switch[qm,
                    "PhaseModulated",
                        keys = {"Data"};
                        {<|"Key" -> "Data", "BaseIdx" -> 1, "Weight" -> 1., "RowModulators" -> None, "StateOps" -> {idOp}, "MixPart" -> baseMixPart|>},
                    "States",
                        keys = {"Cos", "Sin"};
                        {<|"Key" -> "Cos", "BaseIdx" -> 1, "Weight" -> 1., "RowModulators" -> None, "StateOps" -> {idOp}, "MixPart" -> baseMixPart|>,
                         <|"Key" -> "Sin", "BaseIdx" -> 1, "Weight" -> 1., "RowModulators" -> None, "StateOps" -> {buildStateOpRotZ[-Pi / (2. * qOrder)]}, "MixPart" -> baseMixPart|>},
                    "StatesTPPI",
                        keys = {"Cos", "Sin"};
                        rowMod = toPackedRealVec[Table[(-1.)^(k - 1), {k, n1Kernel}]];
                        {<|"Key" -> "Cos", "BaseIdx" -> 1, "Weight" -> 1., "RowModulators" -> rowMod, "StateOps" -> {idOp, buildStateOpRotZ[-Pi / qOrder]}, "MixPart" -> baseMixPart|>,
                         <|"Key" -> "Sin", "BaseIdx" -> 1, "Weight" -> 1., "RowModulators" -> rowMod, "StateOps" -> {buildStateOpRotZ[-Pi / (2. * qOrder)], buildStateOpRotZ[-3. * Pi / (2. * qOrder)]}, "MixPart" -> baseMixPart|>},
                    "EchoAntiEcho",
                        keys = {"Echo", "AntiEcho"};
                        {<|"Key" -> "Echo", "BaseIdx" -> 1, "Weight" -> 1., "RowModulators" -> None, "StateOps" -> {buildStateOpFilter[-qOrder]}, "MixPart" -> baseMixPart|>,
                         <|"Key" -> "AntiEcho", "BaseIdx" -> 1, "Weight" -> 1., "RowModulators" -> None, "StateOps" -> {buildStateOpFilter[qOrder]}, "MixPart" -> baseMixPart|>},
                    "TPPI",
                        keys = {"Data"};
                        {<|"Key" -> "Data", "BaseIdx" -> 1, "Weight" -> 1., "RowModulators" -> None,
                           "StateOps" -> Prepend[Table[buildStateOpRotZ[(m * Pi) / (2. * qOrder)], {m, 1, 3}], idOp],
                           "MixPart" -> baseMixPart|>}
                ]
            ];
            If[fixedSpec["MixCOMPUTE"] === Missing["KeyAbsent", "MixCOMPUTE"], fixedSpec["MixCOMPUTE"] = ConstantArray[<||>, Length[templateComps]]];
            If[!KeyExistsQ[fixedSpec, "MixCOMPUTE"] || Length[fixedSpec["MixCOMPUTE"]] < Length[templateComps],
                fixedSpec["MixCOMPUTE"] = ConstantArray[First[fixedSpec["MixCOMPUTE"]], Length[templateComps]]
            ];

            If[space === "Liouville" && mode === "DoubleDiag",
                fixedSpec["ActRho"] = Union[Flatten[Table[
                    Flatten[Table[
                        Pick[Range[dimL], UnitStep[Abs[Evaluate[
                            Which[
                                op["Kind"] === "Identity", baseStates[[1]],
                                op["Kind"] === "RotZ", op["Mat"] . baseStates[[1]],
                                op["Kind"] === "Filter", op["Mat"] . baseStates[[1]],
                                True, baseStates[[1]]
                            ]
                        ]] - 1.*^-14], 1],
                        {op, tmpl["StateOps"]}
                    ]],
                    {tmpl, templateComps}
                ]]]
            ];

            quadMeta = <|
                "VarPrepQ" -> (prepStageOps =!= None),
                "Rho0" -> rho0VecOrMat,
                "BaseStates" -> baseStates,
                "TemplateComps" -> templateComps
            |>;

            (* Determine chunking & parallel setup *)
            If[IntegerQ[parVal] && parVal > 1 && $KernelCount < parVal,
                Quiet[LaunchKernels[parVal - $KernelCount]]
            ];
            If[TrueQ[parVal] && $KernelCount === 0,
                Quiet[LaunchKernels[]]
            ];
            nK = Which[
                IntegerQ[parVal] && parVal > 1 && $KernelCount >= 1, Min[parVal, $KernelCount],
                TrueQ[parVal] && $KernelCount >= 1, $KernelCount,
                True, 1
            ];
            chunkSize = Which[
                IntegerQ[maxBatch] && maxBatch >= 1, Min[maxBatch, Ceiling[nEns / Max[1, nK]]],
                nK > 1, Ceiling[nEns / nK],
                True, nEns
            ];
            chunkSize = Max[1, chunkSize];
            chunks = Partition[Range[nEns], UpTo[chunkSize]];

            payloads = Table[
                Module[{chunkAtom = <||>},
                    If[mode === "DoubleDiag",
                        If[KeyExistsQ[atomSpecs, "H2CoeffsAll"],
                            chunkAtom["H2"] = atomSpecs["H2"]; chunkAtom["H2Coeffs"] = atomSpecs["H2CoeffsAll"][[idxs]]
                        ];
                        If[KeyExistsQ[atomSpecs, "H1CoeffsAll"],
                            chunkAtom["H1"] = atomSpecs["H1"]; chunkAtom["H1Coeffs"] = atomSpecs["H1CoeffsAll"][[idxs]]
                        ];
                        If[KeyExistsQ[atomSpecs, "H1aCoeffsAll"],
                            chunkAtom["H1a"] = atomSpecs["H1a"]; chunkAtom["H1aCoeffs"] = atomSpecs["H1aCoeffsAll"][[idxs]];
                            chunkAtom["H1b"] = atomSpecs["H1b"]; chunkAtom["H1bCoeffs"] = atomSpecs["H1bCoeffsAll"][[idxs]]
                        ];
                        If[KeyExistsQ[atomSpecs, "L2CoeffsAll"],
                            chunkAtom["L2"] = atomSpecs["L2"]; chunkAtom["L2Coeffs"] = atomSpecs["L2CoeffsAll"][[idxs]]
                        ];
                        If[KeyExistsQ[atomSpecs, "L1CoeffsAll"],
                            chunkAtom["L1"] = atomSpecs["L1"]; chunkAtom["L1Coeffs"] = atomSpecs["L1CoeffsAll"][[idxs]]
                        ];
                        If[KeyExistsQ[atomSpecs, "L1aCoeffsAll"],
                            chunkAtom["L1a"] = atomSpecs["L1a"]; chunkAtom["L1aCoeffs"] = atomSpecs["L1aCoeffsAll"][[idxs]];
                            chunkAtom["L1b"] = atomSpecs["L1b"]; chunkAtom["L1bCoeffs"] = atomSpecs["L1bCoeffsAll"][[idxs]]
                        ],
                        chunkAtom["ThetaChunk"] = thetaList[[idxs]];
                        chunkAtom["TD1Spec"] = atomSpecs["TD1Spec"];
                        chunkAtom["TD2Spec"] = atomSpecs["TD2Spec"];
                        chunkAtom["PrepTDSpec"] = Lookup[atomSpecs, "PrepTDSpec", None];
                        chunkAtom["MixTDSpec"] = Lookup[atomSpecs, "MixTDSpec", None]
                    ];
                    <|
                        "Mode" -> mode,
                        "Space" -> space,
                        "N1" -> n1Kernel,
                        "N2" -> t2Parsed["N"],
                        "Keys" -> keys,
                        "Weights" -> weights[[idxs]],
                        "Time1" -> t1TimeVec,
                        "Time2" -> t2Parsed["TimeVec"],
                        "QObj" -> qObj,
                        "NormDenom" -> normDenom,
                        "QuadMeta" -> quadMeta,
                        "AtomSpecs" -> chunkAtom,
                        "FixedSpec" -> fixedSpec
                    |>
                ],
                {idxs, chunks}
            ];

            resOutList = If[nK > 1 && Length[payloads] > 1,
                DistributeDefinitions[EvaluateCompiledEnsembleChunk2D];
                ParallelMap[EvaluateCompiledEnsembleChunk2D, payloads, Method -> "CoarsestGrained"],
                EvaluateCompiledEnsembleChunk2D /@ payloads
            ];
            If[MemberQ[resOutList, $Failed] || !AllTrue[resOutList, AssociationQ], Return[$Failed]];

            accFinal = AssociationMap[Function[{k}, toPackedComplexMat[Total[(#["DataMap"][k] &) /@ resOutList]]], keys];
            firstAssoc = <|
                "Data" -> Which[
                    qm === "PhaseModulated" || qm === "TPPI", accFinal["Data"],
                    qm === "States" || qm === "StatesTPPI", <|"Cos" -> accFinal["Cos"], "Sin" -> accFinal["Sin"]|>,
                    qm === "EchoAntiEcho", <|"Echo" -> accFinal["Echo"], "AntiEcho" -> accFinal["AntiEcho"]|>
                ],
                "Timing1" -> t1Parsed["Timing"],
                "Timing2" -> t2Parsed["Timing"],
                "Time1" -> t1Parsed["TimeVec"],
                "Time2" -> t2Parsed["TimeVec"],
                "N1" -> n1Log,
                "N2" -> t2Parsed["N"],
                "QuadratureMethod" -> qm,
                "LineBroadeningApplied" -> None,
                "Peaks2D" -> None,
                "Metadata" -> <|
                    "Space" -> space,
                    "Method" -> scmChosen,
                    "SignalCalculationMethod" -> scmChosen,
                    "Chronology" -> chron,
                    "QuadratureOrder" -> qOrder,
                    "QuadratureSpins" -> qSpins,
                    "NormalizationFactor" -> normVal,
                    "EnsembleMembers" -> nEns,
                    "EnsembleCompiled" -> True,
                    "Period" -> If[mode === "COMPUTE2D", resolvedPeriod, None],
                    "CarouselAverage" -> caravOpt,
                    "InitialTimePoint" -> t0Val
                |>
            |>;
            If[qm === "TPPI",
                firstAssoc["TPPITiming1"] = t1TimingForKernel;
                firstAssoc["TPPITime1"] = t1TimeVec
            ];
            resSig = SignalData2D[firstAssoc];
            If[!MatchQ[lbOpt, None | False],
                resSig = ApplyLineBroadening[resSig, LineBroadening -> lbOpt]
            ];
            resSig
        ]
      ], blockSyms]]
    ]
];

(* HoldRest-safe Ensemble dispatch before evaluating inner problem *)
Signal2D[arg1_, moreArgs___] /;
    Module[{ensHits},
        ensHits = Cases[HoldComplete[moreArgs], (Rule | RuleDelayed)[EnsembleAverage, v_] :> HoldComplete[v], {1}];
        ensHits =!= {} && !MatchQ[ReleaseHold[First[ensHits]], False | None]
    ] := runEnsembleSignal2D[HoldComplete[arg1, moreArgs]];

runEnsembleSignal2D[heldAll_HoldComplete] := Module[
    {posArgsHold, optArgsHold, ensVal, parVal, maxBatch, stdEns, rules, weights, nEns,
     multiVars, multiExpandRules, heldExpanded,
     compiledRes, optsCleanHold, worker, batchIdxs, bRes, accData, keyList, firstSig, resAssoc,
     lbOpt, lbApplied = False, nK},
    posArgsHold = Select[heldAll, !MatchQ[Unevaluated[#], _Rule | _RuleDelayed] &];
    optArgsHold = Select[heldAll, MatchQ[Unevaluated[#], _Rule | _RuleDelayed] &];
    ensVal = ReleaseHold[First[Cases[optArgsHold, (Rule | RuleDelayed)[EnsembleAverage, v_] :> HoldComplete[v], {1}]]];
    parVal = With[{hits = Cases[optArgsHold, (Rule | RuleDelayed)[Parallel, v_] :> HoldComplete[v], {1}]},
        If[hits === {}, False, ReleaseHold[First[hits]]]
    ];
    maxBatch = With[{hits = Cases[optArgsHold, (Rule | RuleDelayed)[MaxBatchSize, v_] :> HoldComplete[v], {1}]},
        If[hits === {}, Infinity, ReleaseHold[First[hits]]]
    ];
    stdEns = StandardizeEnsembleAverage[ensVal];
    If[stdEns === $Failed, Return[$Failed]];
    multiVars = Cases[
        Transpose[{stdEns[[1]], If[ListQ[#], Length[#], 0] & /@ stdEns[[2, 1, 1]]}],
        {_, Except[0]}
    ];
    If[multiVars =!= {},
        (SetAttributes[Evaluate[First[#]], NHoldFirst] &) /@ multiVars;
        multiExpandRules = (Rule[First[#], Array[First[#], Last[#]]] &) /@ multiVars;
        posArgsHold = posArgsHold /. multiExpandRules;
        optArgsHold = Join[
            DeleteCases[optArgsHold, (Rule | RuleDelayed)[EnsembleAverage, _], {1}] /. multiExpandRules,
            HoldComplete @@ Cases[optArgsHold, (Rule | RuleDelayed)[EnsembleAverage, _], {1}]
        ];
        heldExpanded = Join[posArgsHold, optArgsHold],
        heldExpanded = heldAll
    ];
    rules = EnsembleSubstitutionRules[stdEns];
    weights = N[EnsembleWeights[stdEns]];
    nEns = Length[rules];

    (* Fast path: compiled operator-atom engine + in-place packed chunk accumulation *)
    compiledRes = Quiet[tryCompileEnsemble2D[heldExpanded, stdEns, rules, weights, parVal, maxBatch]];
    If[MatchQ[compiledRes, _SignalData2D], Return[compiledRes]];

    (* General symbolic batched fallback with in-place accumulation *)
    optsCleanHold = DeleteCases[optArgsHold, (Rule | RuleDelayed)[EnsembleAverage | Parallel | MaxBatchSize | LineBroadening, _], {1}];
    lbOpt = With[{hits = Cases[optArgsHold, (Rule | RuleDelayed)[LineBroadening, v_] :> HoldComplete[v], {1}]},
        If[hits === {}, None, ReleaseHold[First[hits]]]
    ];

    worker = With[{pH = posArgsHold, oH = optsCleanHold},
        Function[sample,
            Module[{pSub = pH /. sample, oSub = oH /. sample, combined},
                combined = Join[pSub, oSub, HoldComplete[EnsembleAverage -> False, Parallel -> False, LineBroadening -> None, ReportSignalCalculationMethod -> False]];
                Apply[Signal2D, combined]
            ]
        ]
    ];

    If[IntegerQ[parVal] && parVal > 1 && $KernelCount < parVal, Quiet[LaunchKernels[parVal - $KernelCount]]];
    If[TrueQ[parVal] && $KernelCount === 0, Quiet[LaunchKernels[]]];
    nK = Which[IntegerQ[parVal] && parVal > 1 && $KernelCount >= 1, Min[parVal, $KernelCount], TrueQ[parVal] && $KernelCount >= 1, $KernelCount, True, 0];

    batchIdxs = If[IntegerQ[maxBatch] && maxBatch >= 1, Partition[Range[nEns], UpTo[maxBatch]], {Range[nEns]}];
    accData = None; firstSig = None;
    Do[
        bRes = If[nK > 0 && Length[bIdx] > 1,
            ParallelMap[worker, rules[[bIdx]], Method -> "CoarsestGrained"],
            Table[worker[rules[[idx]]], {idx, bIdx}]
        ];
        If[MemberQ[bRes, $Failed] || bRes === {}, Return[$Failed, Module]];
        If[firstSig === None, firstSig = First[bRes]];
        Do[
            With[{w = weights[[bIdx[[j]]]], d = bRes[[j]]["Data"]},
                If[accData === None,
                    accData = If[AssociationQ[d],
                        keyList = Keys[d];
                        AssociationMap[Function[{k}, toPackedComplexMat[w * d[k]]], keyList],
                        toPackedComplexMat[w * d]
                    ],
                    If[AssociationQ[d],
                        Do[accData[k] = accData[k] + w * d[k], {k, keyList}],
                        accData = accData + w * d
                    ]
                ]
            ],
            {j, Length[bIdx]}
        ],
        {bIdx, batchIdxs}
    ];
    If[accData === None || firstSig === None, Return[$Failed]];
    If[AssociationQ[accData],
        accData = AssociationMap[Function[{k}, toPackedComplexMat[accData[k]]], Keys[accData]],
        accData = toPackedComplexMat[accData]
    ];
    resAssoc = firstSig["Association"];
    resAssoc["Data"] = accData;
    resAssoc["Peaks2D"] = None;
    resAssoc["Metadata", "EnsembleMembers"] = nEns;
    resAssoc["Metadata", "EnsembleCompiled"] = False;
    With[{outSig = SignalData2D[resAssoc]},
        If[!MatchQ[lbOpt, None | False], ApplyLineBroadening[outSig, LineBroadening -> lbOpt], outSig]
    ]
];

(* Public user call signatures *)
Signal2D[timing1_List, timing2_List, opts:(_Rule|_RuleDelayed)...] /;
    !(Length[timing2] === 3 && (Head[timing1] === List && Length[timing1] === 2 && ListQ[timing1[[1]]] && ListQ[timing1[[2]]])) :=
    Signal2D[{timing1, timing2}, opts];

Signal2D[{timing1_List, timing2_List}, {ev1_, mix_, ev2_}, opts:(_Rule|_RuleDelayed)...] :=
    Signal2D[{timing1, timing2}, Evolution1 -> ev1, Mixing -> mix, Evolution2 -> ev2, opts];

Signal2D[{timing1_List, timing2_List}, opts:(_Rule|_RuleDelayed)...] :=
    runSingleSignal2D[timing1, timing2, opts];

runSingleSignal2D[timing1_List, timing2_List, opts:(_Rule|_RuleDelayed)...] := Module[
    {optRules, forceEvenQ, chron, t1Parsed, t2Parsed, qm, qOrder, qSpins, system, basis, opbasis,
     rho0Opt, prepOpt, ev1Opt, mixOpt, ev2Opt, obsOpt, normOpt, bgOpt, cpOpt, pcOpt, lbOpt, scmOpt,
     periodOpt, t0Val, caravOpt, magSubOpt,
     useHilbOpt, reportQ, bgList, obsOp, obsSpins, cpPair, p1Filt, pMixFilt,
     ev1Spec, ev2Spec, resolvedPeriod, isPeriodicProblem, caravActiveQ,
     t1TimeVec, t1TimingForKernel, n1Log, n1Kernel,
     pcSteps, space, useDiag2Q = False, scmChosen, dimH, dimL, dim, qAdjMat, qVec, qNormSq, normDenom,
     buildRawMix, sameMixQ, prepChron, buildStateForPrep, prepStateFromSpec,
     addPathwayFilterToMix, baseCompList, finalComps,
     h1Diag, h1aDiag, h1bDiag, h2Diag, l1Diag, l1aDiag, l1bDiag, l2Diag, actQ, actRho,
     refSteps, refSingleU, refMatL, diagFailQ,
     nDiv = 1, nSubPerDiv = 4, microSub = 2, nGamma = 1, gammaOffsets = {0}, gammaShifts = {0.},
     tauPrep = 0., tauMix = 0., tSym, prepPerQ, mixPerQ, prepChronsList, mixChronsList, uPhases,
     kernelOut, dataMap, outData, sigAssoc, resSig},

    optRules = Flatten[{opts}];
    chron = Chronology /. optRules /. Options[Signal2D];
    If[!MemberQ[{"LeftToRight", "RightToLeft"}, chron],
        Message[Signal2D::invalchron]; Return[$Failed]
    ];
    forceEvenQ = TrueQ[ForceEvenNumberOfPoints /. optRules /. Options[Signal2D]];
    qm = QuadratureMethod /. optRules /. Options[Signal2D];
    If[!MemberQ[{"States", "StatesTPPI", "EchoAntiEcho", "TPPI", "PhaseModulated", None}, qm],
        Message[Signal2D::badquad, qm]; Return[$Failed]
    ];
    If[qm === None, qm = "PhaseModulated"];

    t1Parsed = parseTimingAxis[timing1, forceEvenQ];
    t2Parsed = parseTimingAxis[timing2, forceEvenQ];
    If[t1Parsed === $Failed || t2Parsed === $Failed,
        Message[Signal2D::badtimes, {timing1, timing2}]; Return[$Failed]
    ];

    system = SpinSystem[];
    If[!MatrixQ[system], Message[Signal2D::fail, "SpinSystem is not set."]; Return[$Failed]];
    basis = Basis /. optRules /. Options[Signal2D];
    opbasis = OperatorBasis /. optRules /. Options[Signal2D];
    dimH = Length[BasisKets[basis]];
    dimL = Length[BasisOperators[opbasis]];

    qOrder = Abs[N[Re[QuadratureOrder /. optRules /. Options[Signal2D]]]];
    If[!NumericQ[qOrder] || qOrder <= 0., Message[Signal2D::badquad, qOrder]; Return[$Failed]];
    qSpins = resolveQuadratureSpins[QuadratureSpins /. optRules /. Options[Signal2D], system];
    If[qSpins === $Failed, Message[Signal2D::fail, "Invalid QuadratureSpins."]; Return[$Failed]];

    obsRaw = Observable /. optRules /. Options[Signal2D];
    obsOp = resolveObservable[obsRaw, system];
    If[obsOp === $Failed, Message[Signal2D::invalobs, obsRaw]; Return[$Failed]];
    obsSpins = Which[
        obsRaw === All, First /@ system,
        ValidSpinQ[obsRaw, system], {obsRaw},
        ListQ[obsRaw] && ValidSpinQ[obsRaw, system], obsRaw,
        True, First /@ system
    ];

    normVal = NormalizationFactor /. optRules /. Options[Signal2D];
    normVal = Which[
        normVal === None || normVal === False, 1.,
        NumericQ[normVal] && Chop[N[Re[normVal]]] != 0., N[Re[normVal]],
        True, Message[Signal2D::fail, "Invalid NormalizationFactor."]; Return[$Failed]
    ];

    bgOpt = BackgroundGenerator /. optRules /. Options[Signal2D];
    bgList = Which[
        ListQ[bgOpt] && Length[bgOpt] === 4,
            If[chron === "RightToLeft", Reverse[bgOpt], bgOpt],
        ListQ[bgOpt],
            Message[Signal2D::invalbg]; Return[$Failed],
        True,
            {bgOpt, bgOpt, bgOpt, bgOpt}
    ];
    (* {bgPrep, bgEv1, bgMix, bgEv2} in chronological LeftToRight order *)

    ev1Opt = Evolution1 /. optRules /. Options[Signal2D];
    mixOpt = Mixing /. optRules /. Options[Signal2D];
    ev2Opt = Evolution2 /. optRules /. Options[Signal2D];
    prepOpt = Preparation /. optRules /. Options[Signal2D];
    rho0Opt = InitialDensityOperator /. optRules /. Options[Signal2D];
    cpOpt = CoherencePathway /. optRules /. Options[Signal2D];
    pcOpt = PhaseCycle /. optRules /. Options[Signal2D];
    lbOpt = LineBroadening /. optRules /. Options[Signal2D];
    scmOpt = SignalCalculationMethod /. optRules /. Options[Signal2D];
    periodOpt = SignalGeneratorPeriod /. optRules /. Options[Signal2D];
    t0Val = InitialTimePoint /. optRules /. Options[Signal2D];
    If[!NumericQ[t0Val], Message[Signal2D::fail, "InitialTimePoint must be numeric."]; Return[$Failed]];
    t0Val = N[Re[t0Val]];
    caravOpt = CarouselAverage /. optRules /. Options[Signal2D];
    If[!(caravOpt === False || caravOpt === None || caravOpt === True || (IntegerQ[caravOpt] && caravOpt >= 1)),
        Message[Signal2D::nocarav]; Return[$Failed]
    ];
    magSubOpt = MagnusSubsteps /. optRules /. Options[Signal2D];
    useHilbOpt = TrueQ[UseHilbertSpaceWhenPossible /. optRules /. Options[Signal2D]];
    reportQ = TrueQ[ReportSignalCalculationMethod /. optRules /. Options[Signal2D]];

    cpPair = resolveCoherencePathway[cpOpt, qSpins, obsSpins];
    If[cpPair === $Failed, Message[Signal2D::fail, "Invalid CoherencePathway."]; Return[$Failed]];
    {p1Filt, pMixFilt} = cpPair;

    ev1Spec = parseEvolution1[ev1Opt, chron];
    ev2Spec = parseEvolution2[ev2Opt, chron];
    If[ev1Spec === $Failed || ev2Spec === $Failed,
        Message[Signal2D::fail, "Invalid Evolution1 or Evolution2 specification."]; Return[$Failed]
    ];

    (* Normalize PhaseCycle list of Associations *)
    pcSteps = Which[
        pcOpt === None || pcOpt === {}, None,
        ListQ[pcOpt] && AllTrue[pcOpt, AssociationQ[#] || MatchQ[#, {(_Rule|_RuleDelayed)...}] &],
            With[{nPC = Length[pcOpt]},
                Table[
                    With[{a = Association[step]},
                        <|
                            "Preparation" -> Lookup[a, "Preparation", prepOpt],
                            "Mixing" -> Lookup[a, "Mixing", mixOpt],
                            "ReceiverPhase" -> N[Re[Lookup[a, "ReceiverPhase", 0.]]],
                            "Weight" -> N[Lookup[a, "Weight", 1. / nPC]]
                        |>
                    ],
                    {step, pcOpt}
                ]
            ],
        True, Message[Signal2D::fail, "PhaseCycle must be None or a list of Associations/rule lists."]; Return[$Failed]
    ];

    (* Check whether all generators and stages are Hilbert-compatible *)
    space = If[useHilbOpt &&
        OperatorQ[rho0Opt] &&
        AllTrue[bgList, hilbertGenQ] &&
        Which[
            MemberQ[{"Static", "Periodic"}, ev1Spec["Kind"]], hilbertGenQ[ev1Spec["Gen1"]],
            ev1Spec["Kind"] === "StaticSplit", hilbertGenQ[ev1Spec["Gen1a"]] && hilbertGenQ[ev1Spec["Gen1b"]] && hilbertEventListQ[ev1Spec["RefocusEvents"]],
            ev1Spec["Kind"] === "Stepped", hilbertEventListQ[ev1Spec["StepEvents"]],
            True, hilbertGenQ[ev1Spec["DirectFn"]] || MatchQ[ev1Spec["DirectFn"], _Function | _TimeDependentGeneratorFunction | _PeriodicFunction]
        ] &&
        Which[
            MemberQ[{"Static", "Periodic"}, ev2Spec["Kind"]], hilbertGenQ[ev2Spec["Gen2"]],
            ev2Spec["Kind"] === "Stepped", hilbertEventListQ[ev2Spec["StepEvents"]],
            True, hilbertGenQ[ev2Spec["DirectFn"]] || MatchQ[ev2Spec["DirectFn"], _Function | _TimeDependentGeneratorFunction | _PeriodicFunction]
        ] &&
        If[pcSteps === None,
            If[AssociationQ[prepOpt],
                AllTrue[Values[prepOpt], hilbertStageQ[#, chron] &],
                hilbertStageQ[prepOpt, chron]
            ] && hilbertStageQ[mixOpt, chron],
            AllTrue[pcSteps, hilbertStageQ[#["Preparation"], chron] && hilbertStageQ[#["Mixing"], chron] &]
        ],
        "Hilbert",
        "Liouville"
    ];
    dim = If[space === "Hilbert", dimH, dimL];

    (* Resolve periodicity & select SignalCalculationMethod *)
    resolvedPeriod = resolve2DPeriod[periodOpt, ev1Spec, ev2Spec, bgList];
    isPeriodicProblem = NumericQ[resolvedPeriod] && resolvedPeriod > 0 &&
        MemberQ[{"Static", "Periodic"}, ev1Spec["Kind"]] &&
        MemberQ[{"Static", "Periodic"}, ev2Spec["Kind"]];
    caravActiveQ = (caravOpt === True) || (IntegerQ[caravOpt] && caravOpt > 1);
    If[caravActiveQ && (!isPeriodicProblem || MemberQ[{"DoubleDiagonalization", "SteppedPropagator"}, scmOpt]),
        Message[Signal2D::nocarav]; Return[$Failed]
    ];

    scmChosen = Which[
        scmOpt === "Direct", "Direct",
        scmOpt === "COMPUTE2D",
            If[isPeriodicProblem,
                "COMPUTE2D",
                Message[Signal2D::badcalc, "COMPUTE2D", "requires periodic time-dependent Evolution1/Evolution2/BackgroundGenerator with positive SignalGeneratorPeriod."]; Return[$Failed]
            ],
        scmOpt === "SteppedPropagator",
            If[ev1Spec["Kind"] =!= "Direct" && ev2Spec["Kind"] =!= "Direct" && resolvedPeriod === All,
                "SteppedPropagator",
                Message[Signal2D::badcalc, "SteppedPropagator", "time-dependent/periodic generators require COMPUTE2D or Direct."]; Return[$Failed]
            ],
        scmOpt === "DoubleDiagonalization",
            If[MemberQ[{"Static", "StaticSplit"}, ev1Spec["Kind"]] && ev2Spec["Kind"] === "Static" && resolvedPeriod === All,
                "DoubleDiagonalization",
                Message[Signal2D::badcalc, "DoubleDiagonalization", "requires time-independent Evolution1 and Evolution2."]; Return[$Failed]
            ],
        scmOpt === Automatic,
            Which[
                isPeriodicProblem, "COMPUTE2D",
                MemberQ[{"Static", "StaticSplit"}, ev1Spec["Kind"]] && ev2Spec["Kind"] === "Static" && resolvedPeriod === All, "DoubleDiagonalization",
                ev1Spec["Kind"] =!= "Direct" && ev2Spec["Kind"] =!= "Direct" && resolvedPeriod === All, "SteppedPropagator",
                True, "Direct"
            ],
        True,
            Message[Signal2D::badcalc, scmOpt, "unknown SignalCalculationMethod."]; Return[$Failed]
    ];

    (* TPPI doubles t1 grid to {2*N1, dt1/2} while preserving logical Timing1 *)
    n1Log = t1Parsed["N"];
    If[qm === "TPPI",
        With[{t1T = t1Parsed["Timing"]},
            t1TimingForKernel = {t1T[[1]], t1T[[1]] + (2 * n1Log - 1) * (0.5 * t1T[[3]]), 0.5 * t1T[[3]]};
            t1TimeVec = toPackedRealVec[t1TimingForKernel[[1]] + t1TimingForKernel[[3]] * Range[0, 2 * n1Log - 1]];
            n1Kernel = 2 * n1Log
        ],
        t1TimingForKernel = t1Parsed["Timing"];
        t1TimeVec = t1Parsed["TimeVec"];
        n1Kernel = n1Log
    ];

    If[isPeriodicProblem && (scmChosen === "COMPUTE2D" || caravActiveQ),
        Module[{q1, q2},
            q1 = Quiet[With[{r = findRationalDenom2D[t1TimingForKernel[[3]] / resolvedPeriod, 64]}, If[IntegerQ[r], r, 8]]];
            q2 = Quiet[With[{r = findRationalDenom2D[t2Parsed["Timing"][[3]] / resolvedPeriod, 64]}, If[IntegerQ[r], r, 8]]];
            nGamma = Which[caravOpt === False || caravOpt === None || caravOpt === 1, 1, TrueQ[caravOpt], Max[8, LCM[q1, q2]], True, caravOpt];
            nDiv = LCM[q1, q2, If[TrueQ[caravOpt], 8, nGamma]];
            If[nDiv > 256, nDiv = Min[256, LCM[q1, q2]]; If[nDiv > 256, nDiv = Max[8, Min[64, q1, q2]]]];
            nSubPerDiv = If[IntegerQ[magSubOpt] && magSubOpt >= 1, magSubOpt, Max[4, Ceiling[64 / nDiv]]];
            microSub = Max[2, nSubPerDiv];
            gammaOffsets = If[nGamma === 1, {0},
                If[Divisible[nDiv, nGamma], Range[0, nDiv - 1, nDiv / nGamma], Table[N[g * nDiv / nGamma], {g, 0, nGamma - 1}]]
            ];
            gammaShifts = N[gammaOffsets * (resolvedPeriod / nDiv)]
        ]
    ];

    (* Compute stage durations for global time tracking *)
    tauPrep = Which[
        pcSteps =!= None, stageEventDuration[normalizeEventsChronological[pcSteps[[1, "Preparation"]], chron]],
        AssociationQ[prepOpt], stageEventDuration[normalizeEventsChronological[prepOpt["Cos"], chron]],
        True, stageEventDuration[normalizeEventsChronological[prepOpt, chron]]
    ];
    tauMix = If[pcSteps =!= None,
        stageEventDuration[normalizeEventsChronological[pcSteps[[1, "Mixing"]], chron]],
        stageEventDuration[normalizeEventsChronological[mixOpt, chron]]
    ];

    (* Precompute Observable norm & matrix/vector representation *)
    If[space === "Hilbert",
        With[{qMat = toPackedComplexMat[MatrixRepresentation[obsOp, basis]]},
            qNormSq = Total[Abs[Flatten[qMat]]^2];
            qAdjMat = toPackedComplexMat[ConjugateTranspose[qMat]];
            normDenom = normVal * qNormSq
        ],
        qVec = toPackedComplexVec[OperatorVectorRepresentation[obsOp, opbasis, basis]];
        qNormSq = Total[Abs[qVec]^2];
        normDenom = normVal * qNormSq
    ];
    If[normDenom == 0., Message[Signal2D::fail, "Observable has zero norm."]; Return[$Failed]];

    (* Prepare base Preparation state and Mixing operator helpers at explicit start time tStart *)
    prepStateFromSpec[pSpec_, Optional[tStart_?NumericQ, 0.]] := Module[{chEvs, rhoInitMat, steps, uMat},
        chEvs = normalizeEventsChronological[pSpec, chron];
        If[chEvs === $Failed, Return[$Failed]];
        If[space === "Hilbert",
            rhoInitMat = toPackedComplexMat[MatrixRepresentation[rho0Opt, basis]];
            steps = buildHilbertStageSteps[chEvs, bgList[[1]], basis, system, tStart];
            If[steps === $Failed, Return[$Failed]];
            ApplyHilbertMixSequence[rhoInitMat, steps],
            rhoInitMat = toPackedComplexVec[OperatorVectorRepresentation[rho0Opt, opbasis, basis]];
            uMat = buildLiouvilleStageMatrix[chEvs, bgList[[1]], opbasis, basis, dimL, tStart];
            If[uMat === $Failed, Return[$Failed]];
            toPackedComplexVec[uMat . rhoInitMat]
        ]
    ];

    buildMixInfo[mSpec_, Optional[tStart_?NumericQ, 0.]] := Module[{chEvs, steps, isSingleU, uSingle, mMat},
        chEvs = normalizeEventsChronological[mSpec, chron];
        If[chEvs === $Failed, Return[$Failed]];
        If[pMixFilt =!= None,
            AppendTo[chEvs, CoherenceOrderFiltrationSuperoperator[obsSpins, pMixFilt]]
        ];
        If[space === "Hilbert",
            steps = buildHilbertStageSteps[chEvs, bgList[[3]], basis, system, tStart];
            If[steps === $Failed, Return[$Failed]];
            If[steps === {},
                steps = {{"Unitary", toPackedComplexMat[IdentityMatrix[dimH]], toPackedComplexMat[IdentityMatrix[dimH]]}}
            ];
            isSingleU = (Length[steps] === 1 && steps[[1, 1]] === "Unitary");
            uSingle = If[isSingleU, steps[[1, 2]], None];
            <|"MixSteps" -> steps, "IsSingleUnitary" -> isSingleU, "SingleUnitary" -> uSingle, "DirectMixEvents" -> chEvs|>,
            mMat = buildLiouvilleStageMatrix[chEvs, bgList[[3]], opbasis, basis, dimL, tStart];
            If[mMat === $Failed, Return[$Failed]];
            <|"MixMat" -> mMat, "DirectMixEvents" -> chEvs|>
        ]
    ];

    (* Build list of base items before applying QuadratureMethod *)
    prepPerQ = (tauPrep > 0. && periodicGenQ[bgList[[1]]]) || periodicGenQ[prepOpt];
    baseCompList = Module[{makeRhoByG},
        makeRhoByG[pSpec_] := If[prepPerQ && nGamma > 1,
            Table[
                With[{st = prepStateFromSpec[pSpec, t0Val + gs]},
                    If[st === $Failed, Return[$Failed, Module]];
                    If[p1Filt =!= None, applySpinFilter[st, space, qSpins, p1Filt, system, basis, opbasis], st]
                ],
                {gs, gammaShifts}
            ],
            With[{st = prepStateFromSpec[pSpec, t0Val]},
                If[st === $Failed, Return[$Failed, Module]];
                {If[p1Filt =!= None, applySpinFilter[st, space, qSpins, p1Filt, system, basis, opbasis], st]}
            ]
        ];
        If[pcSteps === None,
            If[AssociationQ[prepOpt] && qm === "States",
                With[{rcG = makeRhoByG[prepOpt["Cos"]], rsG = makeRhoByG[prepOpt["Sin"]], mInf = buildMixInfo[mixOpt, t0Val]},
                    If[rcG === $Failed || rsG === $Failed || mInf === $Failed, Return[$Failed]];
                    {<|"ExplicitCosByG" -> rcG, "ExplicitSinByG" -> rsG, "MixInfo" -> mInf, "MixSpec" -> mixOpt, "Weight" -> 1.|>}
                ],
                With[{rG = makeRhoByG[prepOpt], mInf = buildMixInfo[mixOpt, t0Val]},
                    If[rG === $Failed || mInf === $Failed, Return[$Failed]];
                    {<|"RhoByGamma" -> rG, "MixInfo" -> mInf, "MixSpec" -> mixOpt, "Weight" -> 1.|>}
                ]
            ],
            sameMixQ = (Length[Union[#["Mixing"] & /@ pcSteps]] === 1);
            If[sameMixQ,
                Module[{mInf = buildMixInfo[pcSteps[[1, "Mixing"]], t0Val], rGList, stG, nGStates},
                    If[mInf === $Failed, Return[$Failed]];
                    nGStates = If[prepPerQ && nGamma > 1, nGamma, 1];
                    rGList = ConstantArray[0, nGStates];
                    Do[
                        stG = makeRhoByG[step["Preparation"]];
                        If[stG === $Failed, Return[$Failed, Module]];
                        If[Length[stG] === 1 && nGStates > 1, stG = ConstantArray[First[stG], nGStates]];
                        Do[rGList[[g]] = rGList[[g]] + (step["Weight"] * Exp[-I * step["ReceiverPhase"]]) * stG[[g]], {g, nGStates}],
                        {step, pcSteps}
                    ];
                    {<|"RhoByGamma" -> (If[space === "Hilbert", toPackedComplexMat[#], toPackedComplexVec[#]] &) /@ rGList,
                       "MixInfo" -> mInf, "MixSpec" -> pcSteps[[1, "Mixing"]], "Weight" -> 1.|>}
                ],
                Table[
                    With[{stG = makeRhoByG[step["Preparation"]], mInf = buildMixInfo[step["Mixing"], t0Val]},
                        If[stG === $Failed || mInf === $Failed, Return[$Failed, Module]];
                        <|"RhoByGamma" -> stG, "MixInfo" -> mInf, "MixSpec" -> step["Mixing"],
                          "Weight" -> step["Weight"] * Exp[-I * step["ReceiverPhase"]]|>
                    ],
                    {step, pcSteps}
                ]
            ]
        ]
    ];
    If[baseCompList === $Failed, Message[Signal2D::fail, "Failed to construct Preparation or Mixing stage."]; Return[$Failed]];

    (* Expand each base item across QuadratureMethod *)
    finalComps = Flatten[Table[
        Module[{mInf = item["MixInfo"], w = item["Weight"], mSpecItem = item["MixSpec"],
                stateField = If[space === "Hilbert", "RhoStates", "RhoVecs"],
                expandQStates, makeComp, rCosByG, rSinByG, rListByG, tpMod},
            expandQStates[rBase_, modeTag_] := Switch[modeTag,
                "PhaseModulated", {rBase},
                "StatesCos", {rBase},
                "StatesSin", {applySpinZRotation[rBase, space, qSpins, -Pi / (2. * qOrder), basis, opbasis]},
                "StatesTPPICos", {rBase, applySpinZRotation[rBase, space, qSpins, -Pi / qOrder, basis, opbasis]},
                "StatesTPPISin",
                    With[{rs = applySpinZRotation[rBase, space, qSpins, -Pi / (2. * qOrder), basis, opbasis]},
                        {rs, applySpinZRotation[rs, space, qSpins, -Pi / qOrder, basis, opbasis]}
                    ],
                "Echo", {applySpinFilter[rBase, space, qSpins, {-qOrder}, system, basis, opbasis]},
                "AntiEcho", {applySpinFilter[rBase, space, qSpins, {qOrder}, system, basis, opbasis]},
                "TPPI", Table[applySpinZRotation[rBase, space, qSpins, (m * Pi) / (2. * qOrder), basis, opbasis], {m, 0, 3}]
            ];
            makeComp[key_, statesByG_List, mod_] := Join[
                <|"Key" -> key, stateField -> First[statesByG], "RhoByGamma" -> statesByG,
                  "RowModulators" -> mod, "Weight" -> w, "MixSpec" -> mSpecItem|>,
                mInf
            ];
            If[KeyExistsQ[item, "ExplicitCosByG"],
                {makeComp["Cos", ({#} &) /@ item["ExplicitCosByG"], None],
                 makeComp["Sin", ({#} &) /@ item["ExplicitSinByG"], None]},
                With[{rByG = item["RhoByGamma"]},
                    Switch[qm,
                        "PhaseModulated",
                            {makeComp["Data", expandQStates[#, "PhaseModulated"] & /@ rByG, None]},
                        "States",
                            {makeComp["Cos", expandQStates[#, "StatesCos"] & /@ rByG, None],
                             makeComp["Sin", expandQStates[#, "StatesSin"] & /@ rByG, None]},
                        "StatesTPPI",
                            tpMod = toPackedRealVec[Table[(-1.)^(k - 1), {k, n1Kernel}]];
                            {makeComp["Cos", expandQStates[#, "StatesTPPICos"] & /@ rByG, tpMod],
                             makeComp["Sin", expandQStates[#, "StatesTPPISin"] & /@ rByG, tpMod]},
                        "EchoAntiEcho",
                            {makeComp["Echo", expandQStates[#, "Echo"] & /@ rByG, None],
                             makeComp["AntiEcho", expandQStates[#, "AntiEcho"] & /@ rByG, None]},
                        "TPPI",
                            {makeComp["Data", expandQStates[#, "TPPI"] & /@ rByG, None]}
                    ]
                ]
            ]
        ],
        {item, baseCompList}
    ], 1];

    (* Execute chosen simulation kernel with seamless Liouville non-diagonalizable fallback *)
    diagFailQ = False;
    If[scmChosen === "DoubleDiagonalization",
        If[space === "Hilbert",
            Module[{h2Mat = generatorToHilbertMatrix[ev2Spec["Gen2"], basis, dimH] + generatorToHilbertMatrix[bgList[[4]], basis, dimH]},
                If[!MatrixQ[h2Mat, NumericQ], Message[Signal2D::fail, "Non-numeric H2."]; Return[$Failed]];
                h2Diag = DiagonalizeHilbertHermitian[h2Mat];
                If[ev1Spec["Mode"] === "Standard",
                    With[{h1Mat = generatorToHilbertMatrix[ev1Spec["Gen1"], basis, dimH] + generatorToHilbertMatrix[bgList[[2]], basis, dimH]},
                        h1Diag = If[Max[Abs[h1Mat - h2Mat]] == 0., h2Diag, DiagonalizeHilbertHermitian[h1Mat]]
                    ];
                    If[reportQ, Message[Signal2D::reportscm, "DoubleDiagonalization", "Hilbert"]];
                    kernelOut = EvaluateHilbertDoubleDiag2D[<|
                        "Ev1Mode" -> "Standard",
                        "H1Diag" -> h1Diag,
                        "H2Diag" -> h2Diag,
                        "Time1" -> t1TimeVec,
                        "Time2" -> t2Parsed["TimeVec"],
                        "QAdjMat" -> qAdjMat,
                        "NormDenom" -> normDenom,
                        "Components" -> finalComps,
                        "ComputePeaks" -> (qm === "PhaseModulated" && Length[finalComps] === 1 && dimH <= 16)
                    |>],
                    Module[{h1aMat, h1bMat},
                        h1aMat = generatorToHilbertMatrix[ev1Spec["Gen1a"], basis, dimH] + generatorToHilbertMatrix[bgList[[2]], basis, dimH];
                        h1bMat = generatorToHilbertMatrix[ev1Spec["Gen1b"], basis, dimH] + generatorToHilbertMatrix[bgList[[2]], basis, dimH];
                        h1aDiag = If[Max[Abs[h1aMat - h2Mat]] == 0., h2Diag, DiagonalizeHilbertHermitian[h1aMat]];
                        h1bDiag = If[Max[Abs[h1bMat - h1aMat]] == 0., h1aDiag, DiagonalizeHilbertHermitian[h1bMat]];
                        refSteps = buildHilbertStageSteps[ev1Spec["RefocusEvents"], bgList[[2]], basis, system];
                        If[refSteps === $Failed, Return[$Failed, Module]];
                        If[refSteps === {}, refSteps = {{"Unitary", toPackedComplexMat[IdentityMatrix[dimH]], toPackedComplexMat[IdentityMatrix[dimH]]}}];
                        refSingleU = If[Length[refSteps] === 1 && refSteps[[1, 1]] === "Unitary", refSteps[[1, 2]], None];
                        If[reportQ, Message[Signal2D::reportscm, "DoubleDiagonalization", "Hilbert"]];
                        kernelOut = EvaluateHilbertDoubleDiag2D[<|
                            "Ev1Mode" -> ev1Spec["Mode"],
                            "CTDuration" -> Lookup[ev1Spec, "CTDuration", 0.],
                            "H1aDiag" -> h1aDiag,
                            "H1bDiag" -> h1bDiag,
                            "RefocusSteps" -> refSteps,
                            "RefocusSingleUnitary" -> refSingleU,
                            "H2Diag" -> h2Diag,
                            "Time1" -> t1TimeVec,
                            "Time2" -> t2Parsed["TimeVec"],
                            "QAdjMat" -> qAdjMat,
                            "NormDenom" -> normDenom,
                            "Components" -> finalComps
                        |>]
                    ]
                ]
            ],
            (* Liouville DoubleDiagonalization with active connected-component block pruning *)
            Module[{l2Mat = generatorToLiouvilleMatrix[ev2Spec["Gen2"], opbasis, basis, dimL] + generatorToLiouvilleMatrix[bgList[[4]], opbasis, basis, dimL]},
                actQ = Pick[Range[dimL], UnitStep[Abs[qVec] - 1.*^-14], 1];
                actRho = Union[Flatten[Table[
                    Flatten[Table[Pick[Range[dimL], UnitStep[Abs[rv] - 1.*^-14], 1], {rv, comp["RhoVecs"]}]],
                    {comp, finalComps}
                ]]];
                If[ev1Spec["Mode"] === "Standard",
                    With[{l1Mat = generatorToLiouvilleMatrix[ev1Spec["Gen1"], opbasis, basis, dimL] + generatorToLiouvilleMatrix[bgList[[2]], opbasis, basis, dimL]},
                        If[Max[Abs[l1Mat - l2Mat]] == 0.,
                            l2Diag = DiagonalizeLiouvilleBlocked[l2Mat, Union[actQ, actRho]];
                            l1Diag = l2Diag,
                            l2Diag = DiagonalizeLiouvilleBlocked[l2Mat, actQ];
                            l1Diag = DiagonalizeLiouvilleBlocked[l1Mat, actRho]
                        ];
                        If[l2Diag === $Failed || l1Diag === $Failed,
                            diagFailQ = True,
                            If[reportQ, Message[Signal2D::reportscm, "DoubleDiagonalization", "Liouville"]];
                            kernelOut = EvaluateLiouvilleDoubleDiag2D[<|
                                "Ev1Mode" -> "Standard",
                                "L1Diag" -> l1Diag,
                                "L2Diag" -> l2Diag,
                                "Time1" -> t1TimeVec,
                                "Time2" -> t2Parsed["TimeVec"],
                                "QVec" -> qVec,
                                "NormDenom" -> normDenom,
                                "Components" -> finalComps
                            |>]
                        ]
                    ],
                    Module[{l1aMat, l1bMat},
                        l2Diag = DiagonalizeLiouvilleBlocked[l2Mat, actQ];
                        If[l2Diag === $Failed,
                            diagFailQ = True,
                            l1aMat = generatorToLiouvilleMatrix[ev1Spec["Gen1a"], opbasis, basis, dimL] + generatorToLiouvilleMatrix[bgList[[2]], opbasis, basis, dimL];
                            l1bMat = generatorToLiouvilleMatrix[ev1Spec["Gen1b"], opbasis, basis, dimL] + generatorToLiouvilleMatrix[bgList[[2]], opbasis, basis, dimL];
                            l1aDiag = If[Max[Abs[l1aMat - l2Mat]] == 0., DiagonalizeLiouvilleBlocked[l1aMat, All], DiagonalizeLiouvilleBlocked[l1aMat, All]];
                            l1bDiag = If[Max[Abs[l1bMat - l1aMat]] == 0., l1aDiag, DiagonalizeLiouvilleBlocked[l1bMat, All]];
                            If[l1aDiag === $Failed || l1bDiag === $Failed,
                                diagFailQ = True,
                                refMatL = buildLiouvilleStageMatrix[ev1Spec["RefocusEvents"], bgList[[2]], opbasis, basis, dimL];
                                If[refMatL === $Failed, Return[$Failed, Module]];
                                If[reportQ, Message[Signal2D::reportscm, "DoubleDiagonalization", "Liouville"]];
                                kernelOut = EvaluateLiouvilleDoubleDiag2D[<|
                                    "Ev1Mode" -> ev1Spec["Mode"],
                                    "CTDuration" -> Lookup[ev1Spec, "CTDuration", 0.],
                                    "L1aDiag" -> l1aDiag,
                                    "L1bDiag" -> l1bDiag,
                                    "RefocusMat" -> refMatL,
                                    "L2Diag" -> l2Diag,
                                    "Time1" -> t1TimeVec,
                                    "Time2" -> t2Parsed["TimeVec"],
                                    "QVec" -> qVec,
                                    "NormDenom" -> normDenom,
                                    "Components" -> finalComps
                                |>]
                            ]
                        ]
                    ]
                ];
                If[diagFailQ,
                    If[scmOpt === "DoubleDiagonalization",
                        Message[Signal2D::badcalc, "DoubleDiagonalization", "Liouvillian is defective / not diagonalizable."]; Return[$Failed]
                    ];
                    scmChosen = "SteppedPropagator"
                ]
            ]
        ]
    ];

    If[scmChosen === "SteppedPropagator",
        Module[{dt1 = t1TimingForKernel[[3]], tmin1 = t1TimingForKernel[[1]],
                dt2 = t2Parsed["Timing"][[3]], tmin2 = t2Parsed["Timing"][[1]],
                init1a = None, step1a, init1b = None, step1b = None, init2 = None, step2,
                refS = {}, refSingle = None, refML = None, diag2 = None, tCT, genFun, stepFun},
            genFun = If[space === "Hilbert",
                Function[{g}, generatorToHilbertMatrix[g, basis, dimH] + generatorToHilbertMatrix[bgList[[2]], basis, dimH]],
                Function[{g}, generatorToLiouvilleMatrix[g, opbasis, basis, dimL] + generatorToLiouvilleMatrix[bgList[[2]], opbasis, basis, dimL]]
            ];
            stepFun = If[space === "Hilbert",
                Function[{m, tau}, toPackedComplexMat[MatrixExp[-I * m * tau]]],
                Function[{m, tau}, toPackedComplexMat[MatrixExp[m * tau]]]
            ];
            If[ev1Spec["Mode"] === "Standard",
                If[ev1Spec["Kind"] === "Static",
                    With[{m1 = genFun[ev1Spec["Gen1"]]},
                        step1a = stepFun[m1, dt1];
                        If[tmin1 != 0., init1a = stepFun[m1, tmin1]]
                    ],
                    step1a = If[space === "Hilbert",
                        buildHilbertPropagatorFromSegment[ev1Spec["StepEvents"], bgList[[2]], basis, dimH],
                        buildLiouvilleStageMatrix[ev1Spec["StepEvents"], bgList[[2]], opbasis, basis, dimL]
                    ]
                ],
                With[{m1a = genFun[ev1Spec["Gen1a"]], m1b = genFun[ev1Spec["Gen1b"]]},
                    If[ev1Spec["Mode"] === "ConstantTime",
                        tCT = ev1Spec["CTDuration"];
                        init1a = stepFun[m1a, 0.5 * (tCT + tmin1)];
                        init1b = stepFun[m1b, 0.5 * (tCT - tmin1)];
                        step1a = stepFun[m1a, 0.5 * dt1];
                        step1b = stepFun[m1b, -0.5 * dt1],
                        If[tmin1 != 0.,
                            init1a = stepFun[m1a, 0.5 * tmin1];
                            init1b = stepFun[m1b, 0.5 * tmin1]
                        ];
                        step1a = stepFun[m1a, 0.5 * dt1];
                        step1b = stepFun[m1b, 0.5 * dt1]
                    ]
                ];
                If[space === "Hilbert",
                    refS = buildHilbertStageSteps[ev1Spec["RefocusEvents"], bgList[[2]], basis, system];
                    If[refS === {}, refS = {{"Unitary", toPackedComplexMat[IdentityMatrix[dimH]], toPackedComplexMat[IdentityMatrix[dimH]]}}];
                    refSingle = If[Length[refS] === 1 && refS[[1, 1]] === "Unitary", refS[[1, 2]], None],
                    refML = buildLiouvilleStageMatrix[ev1Spec["RefocusEvents"], bgList[[2]], opbasis, basis, dimL]
                ]
            ];
            If[ev2Spec["Kind"] === "Static",
                With[{m2 = If[space === "Hilbert",
                        generatorToHilbertMatrix[ev2Spec["Gen2"], basis, dimH] + generatorToHilbertMatrix[bgList[[4]], basis, dimH],
                        generatorToLiouvilleMatrix[ev2Spec["Gen2"], opbasis, basis, dimL] + generatorToLiouvilleMatrix[bgList[[4]], opbasis, basis, dimL]]},
                    step2 = stepFun[m2, dt2];
                    If[tmin2 != 0., init2 = stepFun[m2, tmin2]];
                    diag2 = If[space === "Hilbert",
                        DiagonalizeHilbertHermitian[m2],
                        With[{d = DiagonalizeLiouvilleGeneral[m2]}, If[d === $Failed, None, d]]
                    ]
                ],
                step2 = If[space === "Hilbert",
                    buildHilbertPropagatorFromSegment[ev2Spec["StepEvents"], bgList[[4]], basis, dimH],
                    buildLiouvilleStageMatrix[ev2Spec["StepEvents"], bgList[[4]], opbasis, basis, dimL]
                ]
            ];
            If[reportQ, Message[Signal2D::reportscm, "SteppedPropagator", space]];
            kernelOut = EvaluateSteppedPropagator2D[<|
                "Space" -> space,
                "Ev1Mode" -> ev1Spec["Mode"],
                "N1" -> n1Kernel,
                "N2" -> t2Parsed["N"],
                "Time2" -> t2Parsed["TimeVec"],
                "NormDenom" -> normDenom,
                "InitProp1a" -> init1a,
                "StepProp1a" -> step1a,
                "InitProp1b" -> init1b,
                "StepProp1b" -> step1b,
                "RefocusSteps" -> refS,
                "RefocusSingleUnitary" -> refSingle,
                "RefocusMat" -> refML,
                "InitProp2" -> init2,
                "StepProp2" -> step2,
                "Diag2" -> diag2,
                "QObj" -> If[space === "Hilbert", qAdjMat, qVec],
                "Components" -> finalComps
            |>]
        ]
    ];

    If[scmChosen === "COMPUTE2D",
        Module[{d1, db1, d2, db2, td1, td2, eval1Fn, eval2Fn, sameQ, aSet1, aSet2,
                mixPerQ, buildFallbackASet, computeMixByPhase, cComps},
            tSym = Unique["sdT$"];
            d1 = decomposeGenToAtoms2D[ev1Spec["Gen1"], {}, {}, tSym, space, basis, opbasis, dim];
            db1 = decomposeGenToAtoms2D[bgList[[2]], {}, {}, tSym, space, basis, opbasis, dim];
            d2 = decomposeGenToAtoms2D[ev2Spec["Gen2"], {}, {}, tSym, space, basis, opbasis, dim];
            db2 = decomposeGenToAtoms2D[bgList[[4]], {}, {}, tSym, space, basis, opbasis, dim];

            buildFallbackASet[evGen_, bgGen_] := Module[{evsDiv, divProps},
                evsDiv = DividedEventPropagators[
                    If[space === "Hilbert", {Operator, resolvedPeriod, nDiv}, {Superoperator, resolvedPeriod, nDiv}],
                    DivideEventSpecifications[
                        EventSpecifications[{If[nullGenQ[evGen], Delay[resolvedPeriod], {evGen, resolvedPeriod}]}, Chronology -> "LeftToRight"],
                        nDiv, IncludeRepetitionNumbers -> True
                    ],
                    opbasis, basis, InitialTimePoint -> t0Val, BackgroundGenerator -> bgGen
                ];
                toPackedComplexMat /@ evsDiv
            ];

            sameQ = (ev1Spec["Gen1"] === ev2Spec["Gen2"]) && (bgList[[2]] === bgList[[4]]);
            If[d1 =!= $Failed && db1 =!= $Failed,
                td1 = mergeTwoDecomps2D[d1, db1, dim, True];
                eval1Fn = With[{cf = td1["CoeffFn"], ms = td1["Matrices"]},
                    Function[{tv}, With[{cv = cf[{}, tv]}, Sum[cv[[r]] * ms[[r]], {r, Length[ms]}]]]
                ];
                aSet1 = BuildPeriodicPropagatorsCF4[eval1Fn, t0Val, resolvedPeriod, nDiv, nSubPerDiv, space, dim],
                eval1Fn = With[{eg = ev1Spec["Gen1"], bg = bgList[[2]]},
                    Function[{tv}, genMatDirect2D[evalGenAtNumericTime[eg, tv], space, basis, opbasis, dim] + genMatDirect2D[evalGenAtNumericTime[bg, tv], space, basis, opbasis, dim]]
                ];
                aSet1 = buildFallbackASet[ev1Spec["Gen1"], bgList[[2]]]
            ];
            If[sameQ,
                eval2Fn = eval1Fn;
                aSet2 = aSet1,
                If[d2 =!= $Failed && db2 =!= $Failed,
                    td2 = mergeTwoDecomps2D[d2, db2, dim, True];
                    eval2Fn = With[{cf = td2["CoeffFn"], ms = td2["Matrices"]},
                        Function[{tv}, With[{cv = cf[{}, tv]}, Sum[cv[[r]] * ms[[r]], {r, Length[ms]}]]]
                    ];
                    aSet2 = BuildPeriodicPropagatorsCF4[eval2Fn, t0Val, resolvedPeriod, nDiv, nSubPerDiv, space, dim],
                    eval2Fn = With[{eg = ev2Spec["Gen2"], bg = bgList[[4]]},
                        Function[{tv}, genMatDirect2D[evalGenAtNumericTime[eg, tv], space, basis, opbasis, dim] + genMatDirect2D[evalGenAtNumericTime[bg, tv], space, basis, opbasis, dim]]
                    ];
                    aSet2 = buildFallbackASet[ev2Spec["Gen2"], bgList[[4]]]
                ]
            ];

            mixPerQ = (tauMix > 0. && periodicGenQ[bgList[[3]]]) || periodicGenQ[mixOpt];
            uPhases = If[mixPerQ,
                DeleteDuplicates[Round[Mod[Flatten[Table[gammaShifts[[g]] + tauPrep + t1TimeVec, {g, Length[gammaShifts]}]], resolvedPeriod], 1.*^-11]],
                {0.}
            ];

            computeMixByPhase[compAssoc_] := Module[{mInf, toStep, pkVal},
                toStep[inf_] := If[space === "Hilbert",
                    If[TrueQ[inf["IsSingleUnitary"]],
                        If[Max[Abs[inf["SingleUnitary"] - IdentityMatrix[dimH]]] <= 1.*^-13,
                            <|"Kind" -> "Identity"|>,
                            <|"Kind" -> "Unitary", "U" -> inf["SingleUnitary"], "UH" -> toPackedComplexMat[ConjugateTranspose[inf["SingleUnitary"]]]|>
                        ],
                        <|"Kind" -> "Steps", "Steps" -> inf["MixSteps"]|>
                    ],
                    If[Max[Abs[inf["MixMat"] - IdentityMatrix[dimL]]] <= 1.*^-13,
                        <|"Kind" -> "Identity"|>,
                        <|"Kind" -> "Mat", "Mat" -> inf["MixMat"]|>
                    ]
                ];
                If[!mixPerQ,
                    <|"Static" -> toStep[compAssoc]|>,
                    AssociationMap[
                        Function[{pk},
                            toStep[buildMixInfo[compAssoc["MixSpec"], t0Val + pk]]
                        ],
                        uPhases
                    ]
                ]
            ];

            cComps = Table[
                <|
                    "Key" -> comp["Key"],
                    "Weight" -> comp["Weight"],
                    "RowModulators" -> comp["RowModulators"],
                    "RhoByGamma" -> comp["RhoByGamma"],
                    "MixByPhase" -> computeMixByPhase[comp]
                |>,
                {comp, finalComps}
            ];

            If[reportQ, Message[Signal2D::reportscm, "COMPUTE2D", space]];
            kernelOut = EvaluatePeriodicCOMPUTE2D[<|
                "Space" -> space,
                "Period" -> resolvedPeriod,
                "NDiv" -> nDiv,
                "T0" -> t0Val,
                "ASet1" -> aSet1,
                "ASet2" -> aSet2,
                "SameEv1Ev2" -> sameQ,
                "EvalGen1" -> eval1Fn,
                "EvalGen2" -> eval2Fn,
                "MicroSubsteps" -> microSub,
                "PrepDuration" -> tauPrep,
                "MixDuration" -> tauMix,
                "GammaOffsets" -> gammaOffsets,
                "Time1" -> t1TimeVec,
                "Time2" -> t2Parsed["TimeVec"],
                "QObj" -> If[space === "Hilbert", qAdjMat, qVec],
                "NormDenom" -> normDenom,
                "Components" -> cComps
            |>]
        ]
    ];

    If[scmChosen === "Direct",
        If[reportQ, Message[Signal2D::reportscm, "Direct", space]];
        kernelOut = evaluateDirect2D[<|
            "Space" -> space,
            "Ev1Spec" -> ev1Spec,
            "Ev2Spec" -> ev2Spec,
            "Bg1" -> bgList[[2]],
            "BgMix" -> bgList[[3]],
            "Bg2" -> bgList[[4]],
            "Time1" -> t1TimeVec,
            "Time2" -> t2Parsed["TimeVec"],
            "Timing1" -> t1TimingForKernel,
            "Timing2" -> t2Parsed["Timing"],
            "T0" -> t0Val,
            "PrepDuration" -> tauPrep,
            "MixDuration" -> tauMix,
            "GammaTimeShifts" -> gammaShifts,
            "QObs" -> obsOp,
            "NormDenom" -> normDenom,
            "QNormSq" -> qNormSq,
            "NormFactor" -> normVal,
            "Basis" -> basis,
            "OperatorBasis" -> opbasis,
            "QuadratureSpins" -> qSpins,
            "QuadratureOrder" -> qOrder,
            "System" -> system,
            "Components" -> finalComps
        |>]
    ];

    dataMap = kernelOut["DataMap"];
    outData = Which[
        qm === "PhaseModulated" || qm === "TPPI",
            dataMap["Data"],
        qm === "States" || qm === "StatesTPPI",
            <|"Cos" -> dataMap["Cos"], "Sin" -> dataMap["Sin"]|>,
        qm === "EchoAntiEcho",
            <|"Echo" -> dataMap["Echo"], "AntiEcho" -> dataMap["AntiEcho"]|>
    ];

    sigAssoc = <|
        "Data" -> outData,
        "Timing1" -> t1Parsed["Timing"],
        "Timing2" -> t2Parsed["Timing"],
        "Time1" -> t1Parsed["TimeVec"],
        "Time2" -> t2Parsed["TimeVec"],
        "N1" -> n1Log,
        "N2" -> t2Parsed["N"],
        "QuadratureMethod" -> qm,
        "LineBroadeningApplied" -> None,
        "Peaks2D" -> Lookup[kernelOut, "Peaks2D", None],
        "Metadata" -> <|
            "Space" -> space,
            "Method" -> scmChosen,
            "SignalCalculationMethod" -> scmChosen,
            "Chronology" -> chron,
            "QuadratureOrder" -> qOrder,
            "QuadratureSpins" -> qSpins,
            "NormalizationFactor" -> normVal,
            "InitialTimePoint" -> t0Val,
            "Period" -> If[NumericQ[resolvedPeriod], resolvedPeriod, None],
            "CarouselAverage" -> caravOpt
        |>
    |>;
    If[qm === "TPPI",
        sigAssoc["TPPITiming1"] = t1TimingForKernel;
        sigAssoc["TPPITime1"] = t1TimeVec
    ];

    resSig = SignalData2D[sigAssoc];
    If[!MatchQ[lbOpt, None | False],
        resSig = ApplyLineBroadening[resSig, LineBroadening -> lbOpt]
    ];
    resSig
];

End[];
EndPackage[];
