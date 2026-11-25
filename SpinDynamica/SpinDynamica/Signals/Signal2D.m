(* ::Package:: *)

(* :Title: Signals`Signal2D` *)
(* :Context: Signals`Signal2D` *)
(* :Summary: Packed 2D NMR time-domain container (SignalData2D), frequency-domain spectrum (Spectrum2D), hypercomplex/TPPI 2D Fourier transforms, apodization, phasing, slicing, projections, and 2D contour plotting. *)

BeginPackage["Signals`Signal2D`",
{"Global`",
 "SpinDynamica`",
 "General`UserLevel`",
 "General`Miscellaneous`",
 "Signals`Signal`",
 "Signals`Bruker`"
}];

ClearAll[
    SignalData2D,
    Spectrum2D,
    ApplyApodization,
    Apodization,
    PhaseSpectrum2D,
    PhaseCorrection,
    ExtractSlice,
    ProjectSpectrum2D,
    PlotSpectrum2D,
    ReverseAxes,
    BaseContourLevel,
    ContourFactor,
    PositiveColor,
    NegativeColor,
    ShowProjections,
    PositionType,
    ProjectToTPPI,
    TPPIProcessing,
    IndirectBandwidth,
    StorageMode,
    ConjugateSpectrum2D,
    QuadratureDiagnostics
];

SignalData2D::usage = "SignalData2D[<|...|>] represents a simulated 2D NMR time-domain dataset with rows indexed by indirect time t1 (length N1, or 2*N1 for \"TPPI\" / \"OversampledStates\") and columns indexed by direct time t2 (length N2). Supports property queries sig[\"prop\"], Re, Im, Abs, Times, Expand, ApplyLineBroadening, ApplyApodization, CombineData, ProjectToTPPI, and FT.";
Spectrum2D::usage = "Spectrum2D[<|...|>] represents a processed 2D NMR frequency-domain spectrum with rows indexed by F1 (indirect) and columns by F2 (direct). Lossless hypercomplex quadrature data are stored via U = S_RR + i1*S_IR and V = S_RI + i1*S_II. Supports property queries spec[\"prop\"] (including \"RR\", \"IR\", \"RI\", \"II\", \"HypercomplexMagnitude\", \"Representation\"), Re, Im, Abs, Times, Expand, ReflectSpectrum, ConjugateSpectrum2D, CombineComplementarySpectra, CombineData, PhaseSpectrum2D, ExtractSlice, ProjectSpectrum2D, PlotSpectrum2D, and ListContourPlot.";
ApplyApodization::usage = "ApplyApodization[sig2D, opts] applies a 2D window function to a SignalData2D object (or 2D matrix) using Apodization -> spec or Apodization -> {directSpec, indirectSpec}.";
Apodization::usage = "Apodization is an option for ApplyApodization and FT specifying a window function or {directSpec, indirectSpec} pair. Supported specs include None, Automatic, fwhh, {\"Lorentzian\", fwhh}, {\"Gaussian\", fwhh}, {\"Gaussian\", lb, gb}, \"SineBell\", {\"SineBell\", shift}, \"SquaredSineBell\", and {\"SquaredSineBell\", shift}.";
PhaseSpectrum2D::usage = "PhaseSpectrum2D[spec2D, phaseSpec] applies zero- and first-order phase corrections (in radians) to a Spectrum2D object. phaseSpec may be phi0, {phi0Direct, phi0Indirect}, or {{phi0Direct, phi1Direct}, {phi0Indirect, phi1Indirect}}.";
PhaseCorrection::usage = "PhaseCorrection is an option for FT on SignalData2D specifying phase corrections in radians as phi0, {phi0Direct, phi0Indirect}, or {{phi0Direct, phi1Direct}, {phi0Indirect, phi1Indirect}}.";
ProjectToTPPI::usage = "ProjectToTPPI[sig2D] projects a SignalData2D acquired with QuadratureMethod -> \"OversampledStates\" onto the equivalent native \"TPPI\" dataset via P_n = C_n Cos[n Pi/2] - S_n Sin[n Pi/2].";
TPPIProcessing::usage = "TPPIProcessing is an option for FT on SignalData2D with QuadratureMethod -> \"TPPI\" selecting Automatic (\"RetainedHalf\"), \"RetainedHalf\", \"RealOnly\", \"Staggered\", or \"Hilbert\".";
IndirectBandwidth::usage = "IndirectBandwidth is an option for FT on SignalData2D with QuadratureMethod -> \"OversampledStates\" selecting \"Logical\" (default, N1*zf1 rows over (-SW1/2, +SW1/2]) or \"Full\" (2*N1*zf1 rows over (-SW1, +SW1]).";
StorageMode::usage = "StorageMode is an option for FT on SignalData2D selecting Automatic/\"Hypercomplex\" (default lossless U+i2*V storage) or \"Compact\"/\"Reduced\" (single-matrix F1-complex storage U only).";
ConjugateSpectrum2D::usage = "ConjugateSpectrum2D[spec2D, axis] conjugates the hypercomplex spectrum along \"F1\"/\"Indirect\" (negating i1) or \"F2\"/\"Direct\" (negating i2).";
QuadratureDiagnostics::usage = "QuadratureDiagnostics[sigOrSpec] returns a structured Association reporting quadrature encoding settings, channel Frobenius norms, phase-shift bookkeeping, and diagnostic warnings for a SignalData2D or Spectrum2D object. QuadratureDiagnostics is also an option for Signal2D (Automatic, True, \"Verbose\", False, or None).";
ExtractSlice::usage = "ExtractSlice[spec2D, axis, pos, opts] extracts a 1D slice along \"F2\"/\"Direct\" (at fixed F1) or \"F1\"/\"Indirect\" (at fixed F2). Returns {{freq, val}, ...}.";
ProjectSpectrum2D::usage = "ProjectSpectrum2D[spec2D, axis, opts] projects a 2D spectrum onto \"F1\" (indirect) or \"F2\" (direct) using Method -> \"Sum\", \"Max\", or \"Skyline\". Returns {{freq, val}, ...}.";
PlotSpectrum2D::usage = "PlotSpectrum2D[spec2D, opts] renders a 2D NMR contour plot with geometric positive/negative contour levels, optional reversed NMR axes (ReverseAxes -> True), and optional 1D projections (ShowProjections -> True).";
ReverseAxes::usage = "ReverseAxes is an option for PlotSpectrum2D specifying False, True, or {revDirect, revIndirect} axis reversal.";
BaseContourLevel::usage = "BaseContourLevel is an option for PlotSpectrum2D specifying the lowest contour level as a fraction of the peak magnitude (default 0.05).";
ContourFactor::usage = "ContourFactor is an option for PlotSpectrum2D specifying the geometric spacing factor between adjacent contours (default 1.4).";
PositiveColor::usage = "PositiveColor is an option for PlotSpectrum2D specifying the color for positive contours.";
NegativeColor::usage = "NegativeColor is an option for PlotSpectrum2D specifying the color for negative contours.";
ShowProjections::usage = "ShowProjections is an option for PlotSpectrum2D controlling whether F2 (top) and F1 (side) 1D projections are displayed alongside the 2D contour plot.";
PositionType::usage = "PositionType is an option for ExtractSlice selecting \"Frequency\" (default) or \"Index\".";

SignalData2D::inval = "Invalid SignalData2D specification: `1`.";
SignalData2D::baddata = "Incompatible dimensions or quadrature modes in SignalData2D operation: `1`.";
SignalData2D::tppiprojequiv = "ProjectToTPPI computed the algebraic projection P_n = C_n Cos[n Pi/2] - S_n Sin[n Pi/2], but equivalence to native TPPI acquisition is not guaranteed because `1`.";
Spectrum2D::inval = "Invalid Spectrum2D specification: `1`.";
Spectrum2D::badaxis = "Unknown axis specification `1`; expected \"F1\", \"F2\", \"Indirect\", or \"Direct\".";
Spectrum2D::badarg = "Invalid argument for 2D spectrum operation: `1`.";
Spectrum2D::nof2phase = "Cannot apply direct-dimension (F2) phase correction to a Spectrum2D with Representation -> `1` because direct-dimension imaginary components (RI, II) are not stored. Use default hypercomplex storage or apply F2 phasing inside FT.";
Spectrum2D::tppif1phase = "Applying indirect-dimension (F1) phase correction to a native real-only TPPI spectrum; finite damping and truncation introduce conjugate-mirror dispersion cross-talk. For exact independent F1 phasing, acquire with QuadratureMethod -> \"States\", \"StatesTPPI\", \"EchoAntiEcho\", or \"OversampledStates\".";
Spectrum2D::pmphase = "Spectrum2D has Representation -> \"Complex\" (PhaseModulated); PhaseSpectrum2D applies a coupled complex phase-twist rotation rather than independent hypercomplex F1/F2 absorption phasing.";
Spectrum2D::nocomp = "Component `1` is unavailable for Spectrum2D with Representation -> `2`.";
Spectrum2D::tppihilbertzf = "TPPIProcessing -> \"Hilbert\" with indirect ZeroFill < 2 suffers from causal wrap-around distortion; use ZeroFill -> {zf2, 2} or TPPIProcessing -> \"RetainedHalf\".";
Spectrum2D::repmerge = "Combining Spectrum2D objects with mixed representations; direct-dimension imaginary components (V) are dropped in the combined result.";
ApplyApodization::badwin = "Unrecognized apodization window specification `1`.";

Options[ApplyApodization] = {Apodization -> None};
Options[ReflectSpectrum] = Join[Options[ReflectSpectrum] /. _Options -> {}, {}];
If[!MemberQ[Keys[Options[ReflectSpectrum]], Axis], Options[ReflectSpectrum] = Append[Options[ReflectSpectrum], Axis -> "F1"]];
Options[ExtractSlice] = {PositionType -> "Frequency", AxisUnits -> "Hz", Component -> "Real"};
Options[ProjectSpectrum2D] = {Method -> "Sum", AxisUnits -> "Hz", Component -> "Real"};
Options[PlotSpectrum2D] = {
    Component -> Automatic,
    AxisUnits -> "Hz",
    ReverseAxes -> {False, False},
    Contours -> 10,
    BaseContourLevel -> 0.05,
    ContourFactor -> 1.4,
    PositiveColor -> RGBColor[0.10, 0.35, 0.75],
    NegativeColor -> RGBColor[0.80, 0.18, 0.18],
    ShowProjections -> False,
    PlotRange -> All
};

Begin["`Private`"];

toPackedComplexMat[m_] := If[
    Developer`PackedArrayQ[m, Complex],
    m,
    Developer`ToPackedArray[N[Normal[m]] + 0. I]
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

(* Option pair resolver preserving SpinDynamica Signal.m {direct, indirect} order *)
singleWindowSpecQ[spec_] := MatchQ[spec,
    None | False | Automatic | _?NumericQ | "SineBell" | "SquaredSineBell" |
    {"Gaussian", _?NumericQ, _?NumericQ}
] || (MatchQ[spec, {"Lorentzian" | "Gaussian" | "SineBell" | "SquaredSineBell", _?NumericQ}]);

resolveAxisPair[spec_, default_: None] := Which[
    spec === Missing["KeyAbsent", "Missing"] || MatchQ[spec, _Missing], {default, default},
    singleWindowSpecQ[spec], {spec, spec},
    MatchQ[spec, {_, _}], spec,
    True, {spec, spec}
];

resolveZeroFillPair[zf_] := Module[{norm1, pair},
    norm1[v_] := Which[
        v === True || v === Automatic, 1 -> 2,
        v === False || v === None, 1 -> 1,
        IntegerQ[v] && v >= 1, 1 -> v,
        True, $Failed
    ];
    pair = If[MatchQ[zf, {_, _}], zf, {zf, zf}];
    Map[
        Which[
            # === True || # === Automatic, 2,
            # === False || # === None, 1,
            IntegerQ[#] && # >= 1, #,
            True, Message[Spectrum2D::badarg, "ZeroFill must be True, False, Automatic, None, or integer >= 1."]; Return[$Failed, Module]
        ] &,
        pair
    ]
];

resolveFpmfPair[f_] := Module[{pair},
    pair = If[MatchQ[f, {_, _}], f, {f, f}];
    If[!VectorQ[pair, NumericQ],
        Message[Spectrum2D::badarg, "FirstPointMultiplicationFactor must be numeric or {dir, ind}."];
        Return[$Failed]
    ];
    N[Re[pair]]
];

resolvePhasePair[ph_] := Which[
    ph === None || ph === False || ph === Automatic, {{0., 0.}, {0., 0.}},
    NumericQ[ph], {{N[Re[ph]], 0.}, {0., 0.}},
    MatchQ[ph, {_?NumericQ, _?NumericQ}], {{N[Re[ph[[1]]]], 0.}, {N[Re[ph[[2]]]], 0.}},
    MatchQ[ph, {{_?NumericQ, _?NumericQ}, {_?NumericQ, _?NumericQ}}], N[Re[ph]],
    True,
        Message[Spectrum2D::badarg, "PhaseCorrection must be phi0, {phi0Dir, phi0Ind}, or {{phi0Dir, phi1Dir}, {phi0Ind, phi1Ind}}."];
        $Failed
];

buildWindow1D[tVec_List, timing_List, spec_] := Module[
    {n = Length[tVec], tMin, tMax, dt, fwhh, lb, gb, shift, span, u},
    If[n === 0, Return[{}]];
    If[spec === None || spec === False, Return[Developer`ToPackedArray[ConstantArray[1., n]]]];
    {tMin, tMax, dt} = N[timing];
    span = Abs[tMax - tMin];
    Which[
        spec === Automatic,
            fwhh = If[span > 0., 4. * Log[10.] / span, 0.];
            toPackedRealVec[Exp[-(fwhh / 2.) * (tVec - tMin)]],
        NumericQ[spec],
            toPackedRealVec[Exp[-(N[Re[spec]] / 2.) * (tVec - tMin)]],
        MatchQ[spec, {"Lorentzian", _?NumericQ}],
            fwhh = N[Re[spec[[2]]]];
            toPackedRealVec[Exp[-(fwhh / 2.) * (tVec - tMin)]],
        MatchQ[spec, {"Gaussian", _?NumericQ}],
            fwhh = N[Re[spec[[2]]]];
            toPackedRealVec[Exp[-((tVec - tMin) * fwhh)^2 / (16. * Log[2.])]],
        MatchQ[spec, {"Gaussian", _?NumericQ, _?NumericQ}],
            lb = N[Re[spec[[2]]]]; gb = N[Re[spec[[3]]]];
            If[span == 0. || gb == 0.,
                toPackedRealVec[Exp[(lb / 2.) * (tVec - tMin)]],
                toPackedRealVec[Exp[(lb / 2.) * (tVec - tMin) - ((tVec - tMin) / (gb * span))^2]]
            ],
        spec === "SineBell" || MatchQ[spec, {"SineBell", _?NumericQ}],
            shift = If[ListQ[spec], N[Re[spec[[2]]]], 0.];
            u = If[n > 1, N[Range[0, n - 1]] / (n - 1.), {0.}];
            toPackedRealVec[Sin[shift + (Pi - shift) * u]],
        spec === "SquaredSineBell" || MatchQ[spec, {"SquaredSineBell", _?NumericQ}],
            shift = If[ListQ[spec], N[Re[spec[[2]]]], 0.];
            u = If[n > 1, N[Range[0, n - 1]] / (n - 1.), {0.}];
            toPackedRealVec[Sin[shift + (Pi - shift) * u]^2],
        True,
            Message[ApplyApodization::badwin, spec];
            $Failed
    ]
];

apply2DWindowToData[data_, w1_List, w2_List] := Module[{wMat},
    wMat = toPackedRealMat[Outer[Times, w1, w2]];
    If[AssociationQ[data],
        AssociationMap[#[[1]] -> toPackedComplexMat[wMat * #[[2]]] &, data],
        toPackedComplexMat[wMat * data]
    ]
];

(* Compact formatting & accessors for SignalData2D *)
SignalData2D /: Format[SignalData2D[assoc_Association]] := RowBox[{
    "SignalData2D[",
    ToBoxes[ToString[Lookup[assoc, "N1", "?"]] <> " x " <> ToString[Lookup[assoc, "N2", "?"]]],
    ", Quadrature -> ",
    ToBoxes[Lookup[assoc, "QuadratureMethod", "Unknown"]],
    ", Method -> ",
    ToBoxes[Lookup[assoc, {"Metadata", "Method"}, Lookup[Lookup[assoc, "Metadata", <||>], "SignalCalculationMethod", "Unknown"]]],
    "]"
}] // DisplayForm;

computeChannelNorms2D[data_] := Which[
    AssociationQ[data],
        AssociationMap[Function[{k}, N[Sqrt[Total[Abs[Flatten[data[k]]]^2]]]], Keys[data]],
    MatrixQ[data],
        <|"Data" -> N[Sqrt[Total[Abs[Flatten[data]]^2]]]|>,
    True,
        <||>
];

computeQuadReport2D[objType_String, assoc_Association] := Module[
    {meta, qm, qOrder, qSpins, autoQuadQ, hasChanPrepQ, hasChanMixQ, norms, vals, maxN, minN,
     zeroChans, chanRatio, warnings, prepShifts, rxShifts, status, baseRep, rep},
    meta = Lookup[assoc, "Metadata", <||>];
    qm = Lookup[assoc, "QuadratureMethod", None];
    qOrder = Lookup[meta, "QuadratureOrder", 1];
    qSpins = Lookup[meta, "QuadratureSpins", All];
    autoQuadQ = Lookup[meta, "AutomaticQuadratureEncoding", True];
    hasChanPrepQ = Lookup[meta, "ChannelDependentPreparation", False];
    hasChanMixQ = Lookup[meta, "ChannelDependentMixing", False];
    norms = If[KeyExistsQ[meta, "ChannelNorms"] && AssociationQ[meta["ChannelNorms"]],
        meta["ChannelNorms"],
        computeChannelNorms2D[Lookup[assoc, "Data", <||>]]
    ];
    vals = Values[norms];
    maxN = If[vals === {}, 0., Max[vals]];
    minN = If[vals === {}, 0., Min[vals]];
    zeroChans = If[Length[vals] > 1,
        Select[Keys[norms], (maxN == 0. || norms[#] <= 1.*^-12 * maxN) &],
        {}
    ];
    chanRatio = Which[
        vals === {}, Missing["Unavailable"],
        Length[vals] == 1, 1.,
        minN == 0., Infinity,
        True, N[maxN / minN]
    ];
    warnings = Lookup[meta, "DiagnosticWarnings", {}];
    prepShifts = Lookup[meta, "PreparationPhaseShifts", <||>];
    rxShifts = Lookup[meta, "ReceiverPhaseShifts", <||>];
    status = If[warnings === {} && zeroChans === {}, "OK", "Warning"];
    baseRep = <|
        "ObjectType" -> objType,
        "QuadratureMethod" -> qm,
        "QuadratureOrder" -> qOrder,
        "QuadratureSpins" -> qSpins,
        "AutomaticQuadratureEncoding" -> autoQuadQ,
        "HeteronuclearDetected" -> Lookup[meta, "HeteronuclearDetected", Missing["Unavailable"]],
        "ChannelDependentPreparation" -> hasChanPrepQ,
        "ChannelDependentMixing" -> hasChanMixQ,
        "ChannelNorms" -> norms,
        "ChannelRatio" -> chanRatio,
        "ZeroChannels" -> zeroChans,
        "PreparationPhaseShifts" -> prepShifts,
        "ReceiverPhaseShifts" -> rxShifts,
        "Warnings" -> warnings,
        "Status" -> status
    |>;
    If[KeyExistsQ[meta, "NativeTPPIEquivalence"],
        baseRep["NativeTPPIEquivalence"] = meta["NativeTPPIEquivalence"]
    ];
    If[objType === "Spectrum2D",
        rep = Lookup[assoc, "Representation", "Complex"];
        Join[baseRep, <|
            "Representation" -> rep,
            "HypercomplexQ" -> (MemberQ[{"Hypercomplex", "TPPIHypercomplex"}, rep] && getSpectrumV[assoc] =!= None),
            "QuadratureCapabilities" -> Lookup[assoc, "QuadratureCapabilities", buildCapabilities[rep, qm, Join[Lookup[assoc, "Processing", <||>], <|"HasTPPIDCRow" -> KeyExistsQ[assoc, "TPPIDCRow"]|>]]],
            "Processing" -> Lookup[assoc, "Processing", <||>]
        |>],
        baseRep
    ]
];

QuadratureDiagnostics[sig:SignalData2D[assoc_Association]] := computeQuadReport2D["SignalData2D", assoc];
QuadratureDiagnostics[spec:Spectrum2D[assoc_Association]] := computeQuadReport2D["Spectrum2D", assoc];

SignalData2D[assoc_Association][prop_String] := Switch[
    prop,
    "Association", assoc,
    "Keys", Keys[assoc],
    "Dimensions", {assoc["N1"], assoc["N2"]},
    "QuadratureReport", QuadratureDiagnostics[SignalData2D[assoc]],
    "Cos" | "Sin" | "Echo" | "AntiEcho",
        If[AssociationQ[assoc["Data"]] && KeyExistsQ[assoc["Data"], prop],
            assoc["Data"][prop],
            Missing["KeyAbsent", prop]
        ],
    _, Lookup[assoc, prop, Lookup[Lookup[assoc, "Metadata", <||>], prop, Missing["KeyAbsent", prop]]]
];

SignalData2D[assoc_Association][k1_Integer, k2_Integer] := If[
    AssociationQ[assoc["Data"]],
    AssociationMap[#[[1]] -> #[[2]][[k1, k2]] &, assoc["Data"]],
    assoc["Data"][[k1, k2]]
];

extractMainSignalMatrix[assoc_Association, compOpt_: Automatic] := Module[{d = assoc["Data"], qm = assoc["QuadratureMethod"]},
    If[!AssociationQ[d], Return[d]];
    Which[
        StringQ[compOpt] && KeyExistsQ[d, compOpt], d[compOpt],
        qm === "EchoAntiEcho", Lookup[d, "Echo", First[Values[d]]],
        True, Lookup[d, "Cos", First[Values[d]]]
    ]
];

SignalData2D /: Re[SignalData2D[assoc_Association]] := With[{d = assoc["Data"]},
    If[AssociationQ[d],
        AssociationMap[#[[1]] -> toPackedRealMat[Re[#[[2]]]] &, d],
        toPackedRealMat[Re[d]]
    ]
];

SignalData2D /: Im[SignalData2D[assoc_Association]] := With[{d = assoc["Data"]},
    If[AssociationQ[d],
        AssociationMap[#[[1]] -> toPackedRealMat[Im[#[[2]]]] &, d],
        toPackedRealMat[Im[d]]
    ]
];

SignalData2D /: Abs[SignalData2D[assoc_Association]] := With[{d = assoc["Data"]},
    If[AssociationQ[d],
        AssociationMap[#[[1]] -> toPackedRealMat[Abs[#[[2]]]] &, d],
        toPackedRealMat[Abs[d]]
    ]
];

SignalData2D /: HoldPattern[Times[a___, sig:SignalData2D[assoc_Association], b___]] /;
    AllTrue[{a, b}, NumericQ] := Module[{c = N[Times[a, b]], d = assoc["Data"], newAssoc},
    newAssoc = assoc;
    newAssoc["Data"] = If[AssociationQ[d],
        AssociationMap[#[[1]] -> toPackedComplexMat[c * #[[2]]] &, d],
        toPackedComplexMat[c * d]
    ];
    If[MatrixQ[Lookup[assoc, "Peaks2D", None]],
        newAssoc["Peaks2D"] = MapAt[c * # &, assoc["Peaks2D"], {All, 1}]
    ];
    SignalData2D[newAssoc]
];

activeIndirectTimeAndTiming[assoc_Association] := Which[
    assoc["QuadratureMethod"] === "TPPI" && KeyExistsQ[assoc, "TPPITime1"],
        {assoc["TPPITime1"], Lookup[assoc, "TPPITiming1", assoc["Timing1"]]},
    assoc["QuadratureMethod"] === "OversampledStates" && KeyExistsQ[assoc, "OversampledTime1"],
        {assoc["OversampledTime1"], Lookup[assoc, "OversampledTiming1", assoc["Timing1"]]},
    True,
        {assoc["Time1"], assoc["Timing1"]}
];

SignalData2D /: Expand[SignalData2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {embedQ, compOpt, mat, t1, t2, n1, n2},
    embedQ = TrueQ[EmbedCoordinates /. {opts} /. {EmbedCoordinates -> True}];
    compOpt = Component /. {opts} /. {Component -> Automatic};
    mat = extractMainSignalMatrix[assoc, compOpt];
    If[!embedQ, Return[mat]];
    t1 = First[activeIndirectTimeAndTiming[assoc]];
    t2 = assoc["Time2"];
    {n1, n2} = Dimensions[mat];
    Flatten[Table[{t1[[i1]], t2[[i2]], mat[[i1, i2]]}, {i1, n1}, {i2, n2}], 1]
];

ApplyApodization[mat_?MatrixQ, opts:(_Rule|_RuleDelayed)...] := Module[
    {spec, dirSpec, indSpec, n1, n2, t1, t2, w1, w2},
    spec = Apodization /. {opts} /. Options[ApplyApodization];
    {dirSpec, indSpec} = resolveAxisPair[spec, None];
    If[dirSpec === None && indSpec === None, Return[mat]];
    {n1, n2} = Dimensions[mat];
    t1 = N[Range[0, n1 - 1]];
    t2 = N[Range[0, n2 - 1]];
    w1 = buildWindow1D[t1, {0., Max[1., n1 - 1.], 1.}, indSpec];
    w2 = buildWindow1D[t2, {0., Max[1., n2 - 1.], 1.}, dirSpec];
    If[w1 === $Failed || w2 === $Failed, Return[$Failed]];
    apply2DWindowToData[mat, w1, w2]
];

SignalData2D /: ApplyApodization[sig:SignalData2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {spec, dirSpec, indSpec, t1Vec, t1Timing, w1, w2, newAssoc},
    spec = Apodization /. {opts} /. Options[ApplyApodization];
    {dirSpec, indSpec} = resolveAxisPair[spec, None];
    If[(dirSpec === None || dirSpec === False) && (indSpec === None || indSpec === False), Return[sig]];
    {t1Vec, t1Timing} = activeIndirectTimeAndTiming[assoc];
    w1 = buildWindow1D[t1Vec, t1Timing, indSpec];
    w2 = buildWindow1D[assoc["Time2"], assoc["Timing2"], dirSpec];
    If[w1 === $Failed || w2 === $Failed, Return[$Failed]];
    newAssoc = assoc;
    newAssoc["Data"] = apply2DWindowToData[assoc["Data"], w1, w2];
    newAssoc["Peaks2D"] = None;
    SignalData2D[newAssoc]
];

ProjectToTPPI[SignalData2D[assoc_Association]] := Module[
    {d, matC, matS, nRows, cVec, sVec, pMat, t1OS, tm1OS, meta, reasons = {}, equivStatus, newAssoc},
    d = Lookup[assoc, "Data", $Failed];
    If[assoc["QuadratureMethod"] =!= "OversampledStates" || !AssociationQ[d] ||
       !KeyExistsQ[d, "Cos"] || !KeyExistsQ[d, "Sin"],
        Message[SignalData2D::baddata, "ProjectToTPPI requires a SignalData2D with QuadratureMethod -> \"OversampledStates\"."];
        Return[$Failed]
    ];
    matC = toPackedComplexMat[d["Cos"]];
    matS = toPackedComplexMat[d["Sin"]];
    nRows = Length[matC];
    cVec = toPackedRealVec[Cos[Range[0., nRows - 1.] * (Pi / 2.)]];
    sVec = toPackedRealVec[Sin[Range[0., nRows - 1.] * (Pi / 2.)]];
    pMat = toPackedComplexMat[cVec * matC - sVec * matS];
    {t1OS, tm1OS} = activeIndirectTimeAndTiming[assoc];
    meta = Lookup[assoc, "Metadata", <||>];
    If[TrueQ[Lookup[meta, "ChannelDependentPreparation", False]],
        AppendTo[reasons, "ChannelDependentPreparation -> True"]
    ];
    If[TrueQ[Lookup[meta, "ChannelDependentMixing", False]],
        AppendTo[reasons, "ChannelDependentMixing -> True"]
    ];
    If[Lookup[meta, "AutomaticQuadratureEncoding", True] === False &&
       !TrueQ[Lookup[meta, "ChannelDependentPreparation", False]] &&
       !TrueQ[Lookup[meta, "ChannelDependentMixing", False]],
        AppendTo[reasons, "AutomaticQuadratureEncoding -> False"]
    ];
    If[reasons =!= {},
        equivStatus = "NotGuaranteed";
        Message[SignalData2D::tppiprojequiv, StringRiffle[reasons, ", "]],
        equivStatus = "ConditionalOnCoherenceOrders"
    ];
    newAssoc = KeyDrop[assoc, {"OversampledTiming1", "OversampledTime1", "OversamplingFactor1"}];
    newAssoc["Data"] = pMat;
    newAssoc["QuadratureMethod"] = "TPPI";
    newAssoc["TPPITiming1"] = tm1OS;
    newAssoc["TPPITime1"] = t1OS;
    newAssoc["Peaks2D"] = None;
    newAssoc["Metadata"] = Join[
        meta,
        <|
            "ProjectedFromOversampledStates" -> True,
            "ProjectionFormula" -> "C_n Cos[n Pi/2] - S_n Sin[n Pi/2]",
            "NativeTPPIEquivalence" -> equivStatus
        |>
    ];
    SignalData2D[newAssoc]
];

SignalData2D /: ApplyLineBroadening[sig:SignalData2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {lbSpec, dirSpec, indSpec, newSig, newAssoc},
    lbSpec = LineBroadening /. {opts} /. Options[ApplyLineBroadening];
    {dirSpec, indSpec} = resolveAxisPair[lbSpec, None];
    If[(dirSpec === None || dirSpec === False) && (indSpec === None || indSpec === False), Return[sig]];
    newSig = ApplyApodization[sig, Apodization -> {dirSpec, indSpec}];
    If[newSig === $Failed, Return[$Failed]];
    newAssoc = First[newSig];
    newAssoc["LineBroadeningApplied"] = {dirSpec, indSpec};
    SignalData2D[newAssoc]
];

CombineData[sigs:{_SignalData2D..}, opts:(_Rule|_RuleDelayed)...] :=
    CombineData[sigs, ConstantArray[1. / Length[sigs], Length[sigs]], opts];

CombineData[sigs:{_SignalData2D..}, weights_List, opts:(_Rule|_RuleDelayed)...] := Module[
    {w = N[weights], firstA, d1, qm, dataList, combData, newAssoc},
    If[Length[sigs] =!= Length[weights],
        Message[CombineData::incompw]; Return[$Failed]
    ];
    firstA = First[sigs]["Association"];
    d1 = firstA["Data"];
    qm = firstA["QuadratureMethod"];
    If[!AllTrue[Rest[sigs], (
        #["Dimensions"] === First[sigs]["Dimensions"] &&
        #["QuadratureMethod"] === qm &&
        Max[Abs[N[#["Timing1"] - firstA["Timing1"]]]] <= 10^-11 * (1. + Max[Abs[N[firstA["Timing1"]]]]) &&
        Max[Abs[N[#["Timing2"] - firstA["Timing2"]]]] <= 10^-11 * (1. + Max[Abs[N[firstA["Timing2"]]]])
    ) &],
        Message[SignalData2D::baddata, "All SignalData2D datasets must match in dimensions, timings, and QuadratureMethod."];
        Return[$Failed]
    ];
    dataList = (#["Data"]) & /@ sigs;
    combData = If[AssociationQ[d1],
        AssociationMap[
            Function[{kv},
                With[{k = kv[[1]]},
                    k -> toPackedComplexMat[Sum[w[[i]] * dataList[[i]][k], {i, Length[sigs]}]]
                ]
            ],
            d1
        ],
        toPackedComplexMat[Sum[w[[i]] * dataList[[i]], {i, Length[sigs]}]]
    ];
    newAssoc = firstA;
    newAssoc["Data"] = combData;
    newAssoc["Peaks2D"] = None;
    SignalData2D[newAssoc]
];

(* Core packed 1D-per-axis and 2D FT helpers matching Signal.m FourierParameters -> {1, -1} *)
apply1DLinearPhase[mat_?MatrixQ, phc0_?NumericQ, phc1_?NumericQ, axis_Integer] := Module[
    {m = N[Re[phc0]], s = N[Re[phc1]], dims = Dimensions[mat], vec},
    If[m == 0. && s == 0., Return[mat]];
    If[axis === 2,
        vec = Developer`ToPackedArray[Exp[I * (m + s * Range[0., dims[[2]] - 1.] / Max[1., dims[[2]] - 1.])]];
        toPackedComplexMat[Transpose[vec * Transpose[mat]]],
        vec = Developer`ToPackedArray[Exp[I * (m + s * Range[0., dims[[1]] - 1.] / Max[1., dims[[1]] - 1.])]];
        toPackedComplexMat[vec * mat]
    ]
];

ftDirectAxis[mat_?MatrixQ, fpmf2_?NumericQ, zf2_Integer, phcDir_List] := Module[
    {n1, n2, n2ZF, work, padded, spec},
    {n1, n2} = Dimensions[mat];
    n2ZF = zf2 * n2;
    work = toPackedComplexMat[mat];
    If[fpmf2 != 1., work[[All, 1]] *= N[fpmf2]];
    padded = If[n2ZF > n2,
        toPackedComplexMat[PadRight[work, {n1, n2ZF}, 0. + 0. I]],
        work
    ];
    spec = Developer`ToPackedArray[Fourier[#, FourierParameters -> {1, -1}] & /@ padded] / N[n2];
    spec = toPackedComplexMat[RotateRight[spec, {0, Ceiling[n2ZF / 2] - 1}]];
    apply1DLinearPhase[spec, phcDir[[1]], phcDir[[2]], 2]
];

ftIndirectAxis[mat_?MatrixQ, fpmf1_?NumericQ, zf1_Integer, phcInd_List] := Module[
    {n1, n2, n1ZF, work, padded, specT},
    {n1, n2} = Dimensions[mat];
    n1ZF = zf1 * n1;
    work = toPackedComplexMat[mat];
    If[fpmf1 != 1., work[[1]] *= N[fpmf1]];
    padded = If[n1ZF > n1,
        toPackedComplexMat[PadRight[work, {n1ZF, n2}, 0. + 0. I]],
        work
    ];
    specT = Developer`ToPackedArray[Fourier[#, FourierParameters -> {1, -1}] & /@ Transpose[padded]] / N[n1];
    specT = toPackedComplexMat[RotateRight[Transpose[specT], {Ceiling[n1ZF / 2] - 1, 0}]];
    apply1DLinearPhase[specT, phcInd[[1]], phcInd[[2]], 1]
];

buildFrequencyAxisHz[nPts_Integer, dt_?NumericQ] := Module[{swHz = 1. / N[dt]},
    If[EvenQ[nPts],
        toPackedRealVec[Range[-(nPts / 2) + 1, nPts / 2] * (swHz / N[nPts])],
        toPackedRealVec[Range[-Floor[nPts / 2], Floor[nPts / 2]] * (swHz / N[nPts])]
    ]
];

resolveTPPIProcessing[opt_] := Which[
    opt === Automatic, "RetainedHalf",
    MemberQ[{"RetainedHalf", "RealOnly", "Staggered", "Hilbert"}, opt], opt,
    True, Message[Spectrum2D::badarg, "TPPIProcessing must be Automatic, \"RetainedHalf\", \"RealOnly\", \"Staggered\", or \"Hilbert\"."]; $Failed
];

resolveIndirectBandwidth[opt_] := Which[
    opt === Automatic || opt === "Logical", "Logical",
    opt === "Full", "Full",
    True, Message[Spectrum2D::badarg, "IndirectBandwidth must be \"Logical\" or \"Full\"."]; $Failed
];

resolveStorageMode[opt_] := Which[
    opt === Automatic || opt === "Hypercomplex" || opt === Full, True,
    opt === "Compact" || opt === "Reduced", False,
    True, Message[Spectrum2D::badarg, "StorageMode must be Automatic, \"Hypercomplex\", or \"Compact\"."]; $Failed
];

buildCapabilities[rep_String, qm_, extra_Association: <||>] := Module[
    {tppiProc, indBand, hasDCRow, halfBinModQ, f1ModQ, apoSpec, zeroingApoQ, presFullHyper, presProc, presRaw},
    tppiProc = Lookup[extra, "TPPIProcessing", Automatic];
    indBand = Lookup[extra, "IndirectBandwidth", "Logical"];
    hasDCRow = TrueQ[Lookup[extra, "HasTPPIDCRow", False]];
    halfBinModQ = TrueQ[Lookup[extra, "TPPIHalfBinModulated", False]];
    f1ModQ = TrueQ[Lookup[extra, "TPPIF1Modified", False]];
    apoSpec = Lookup[extra, "Apodization", {None, None}];
    zeroingApoQ = MemberQ[Flatten[{apoSpec}], "SineBell" | {"SineBell", _} | "SquaredSineBell" | {"SquaredSineBell", _}];
    presFullHyper = (rep === "Hypercomplex") && !(qm === "OversampledStates" && indBand === "Logical");
    presProc = Which[
        rep === "Hypercomplex", presFullHyper,
        rep === "TPPIHypercomplex" && qm === "TPPI",
            !f1ModQ && MemberQ[{"RetainedHalf", "Staggered"}, tppiProc] && (hasDCRow || halfBinModQ),
        True, False
    ];
    presRaw = presProc && !zeroingApoQ;
    <|
        "HasDirectImaginary" -> MemberQ[{"Hypercomplex", "TPPIHypercomplex"}, rep],
        "HasIndirectQuadrature" -> MemberQ[{"Hypercomplex", "ReducedComplex", "TPPIHypercomplex"}, rep],
        "SupportsIndependentF2Phase" -> MemberQ[{"Hypercomplex", "TPPIHypercomplex"}, rep],
        "SupportsIndependentF1Phase" -> (rep === "Hypercomplex"),
        "SupportsMirrorFreeF1Phase" -> (rep === "Hypercomplex"),
        "PreservesFullHypercomplexSpectrum" -> presFullHyper,
        "PreservesProcessedTimeDomain" -> presProc,
        "PreservesRawAcquisition" -> presRaw,
        "InvertibleToTimeDomain" -> presProc,
        "QuadratureMethod" -> qm
    |>
];

ftHypercomplexPair[mCos_?MatrixQ, mSin_?MatrixQ, fpmfPair_List, zfPair_List, phPair_List, storeFullQ_?BooleanQ] := Module[
    {f2c, f2s, uMat, vMat},
    f2c = ftDirectAxis[mCos, fpmfPair[[1]], zfPair[[1]], phPair[[1]]];
    f2s = ftDirectAxis[mSin, fpmfPair[[1]], zfPair[[1]], phPair[[1]]];
    uMat = ftIndirectAxis[toPackedComplexMat[Re[f2c] + I * Re[f2s]], fpmfPair[[2]], zfPair[[2]], phPair[[2]]];
    vMat = If[storeFullQ,
        ftIndirectAxis[toPackedComplexMat[Im[f2c] + I * Im[f2s]], fpmfPair[[2]], zfPair[[2]], phPair[[2]]],
        None
    ];
    {uMat, vMat}
];

ftTPPIRealHalf[realMat_?MatrixQ, fpmf1_?NumericQ, n1Log_Integer, n1ZF_Integer, phInd_List] := Module[
    {n1TPPI, n2ZF, n1TPPIZF, work, padded, raw, dcRow, half, modVec},
    {n1TPPI, n2ZF} = Dimensions[realMat];
    n1TPPIZF = 2 * n1ZF;
    work = toPackedComplexMat[realMat];
    If[fpmf1 != 1., work[[1]] *= N[fpmf1]];
    padded = If[n1TPPIZF > n1TPPI,
        toPackedComplexMat[PadRight[work, {n1TPPIZF, n2ZF}, 0. + 0. I]],
        work
    ];
    If[EvenQ[n1ZF],
        raw = Transpose[Developer`ToPackedArray[Fourier[#, FourierParameters -> {1, -1}] & /@ Transpose[padded]]] / N[n1Log];
        dcRow = toPackedRealVec[Re[raw[[1]]]];
        half = toPackedComplexMat[raw[[2 ;; n1ZF + 1]]];
        {apply1DLinearPhase[half, phInd[[1]], phInd[[2]], 1], dcRow, False},
        modVec = Developer`ToPackedArray[Exp[-I * Pi * Range[0., n1TPPIZF - 1.] / N[n1TPPIZF]]];
        raw = Transpose[Developer`ToPackedArray[Fourier[#, FourierParameters -> {1, -1}] & /@ Transpose[toPackedComplexMat[modVec * padded]]]] / N[n1Log];
        half = toPackedComplexMat[raw[[1 ;; n1ZF]]];
        {apply1DLinearPhase[half, phInd[[1]], phInd[[2]], 1], None, True}
    ]
];

hilbertReconF1Real[reMat_?MatrixQ] := Module[{n1, n2, tHalf, mask, cTmp},
    {n1, n2} = Dimensions[reMat];
    tHalf = toPackedComplexMat[Transpose[Developer`ToPackedArray[InverseFourier[#, FourierParameters -> {1, -1}] & /@ Transpose[toPackedComplexMat[reMat]]]]];
    mask = toPackedRealVec[Table[
        Which[
            k == 1, 1.,
            k <= Ceiling[n1 / 2], 2.,
            EvenQ[n1] && k == n1 / 2 + 1, 1.,
            True, 0.
        ],
        {k, n1}
    ]];
    tHalf = toPackedComplexMat[mask * tHalf];
    cTmp = toPackedComplexMat[Transpose[Developer`ToPackedArray[Fourier[#, FourierParameters -> {1, -1}] & /@ Transpose[tHalf]]]];
    toPackedComplexMat[toPackedRealMat[Re[reMat]] + I * toPackedRealMat[Im[cTmp]]]
];

phaseUVF2[uMat_?MatrixQ, vMat_?MatrixQ, ph02_?NumericQ, ph12_?NumericQ] := Module[
    {n2 = Dimensions[uMat][[2]], ramp, c2, s2},
    If[ph02 == 0. && ph12 == 0., Return[{toPackedComplexMat[uMat], toPackedComplexMat[vMat]}]];
    ramp = N[Re[ph02]] + N[Re[ph12]] * Range[0., n2 - 1.] / Max[1., n2 - 1.];
    c2 = toPackedRealVec[Cos[ramp]];
    s2 = toPackedRealVec[Sin[ramp]];
    {toPackedComplexMat[Transpose[c2 * Transpose[uMat] - s2 * Transpose[vMat]]],
     toPackedComplexMat[Transpose[s2 * Transpose[uMat] + c2 * Transpose[vMat]]]}
];

SignalData2D /: FT[sig:SignalData2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {lbOpt, apoOpt, zfOpt, fpmfOpt, phOpt, tppiProcOpt, indBandOpt, storageOpt,
     zfPair, fpmfPair, phPair, tppiProc, indBand, storeFullQ,
     sigProc, procAssoc, dProc, qm, n1Log, n2, dt1, dt2, n1ZF, n2ZF,
     n1Out, dt1Out, uMat, vMat = None, repType, tppiDC = None, tppiHalfBinModQ = False,
     f2c, f2s, uFull, vFull, idxRange, mCos, mSin, dcU, dcV,
     cEven, sOdd, hdVec, u0, v0, freq1, freq2, specAssoc},

    lbOpt = LineBroadening /. {opts} /. Options[FT];
    apoOpt = Apodization /. {opts} /. Options[ApplyApodization];
    zfOpt = ZeroFill /. {opts} /. Options[FT];
    fpmfOpt = FirstPointMultiplicationFactor /. {opts} /. Options[FT];
    phOpt = PhaseCorrection /. {opts} /. {PhaseCorrection -> {{0., 0.}, {0., 0.}}};
    tppiProcOpt = TPPIProcessing /. {opts} /. {TPPIProcessing -> Automatic};
    indBandOpt = IndirectBandwidth /. {opts} /. {IndirectBandwidth -> "Logical"};
    storageOpt = StorageMode /. {opts} /. {StorageMode -> Automatic};

    zfPair = resolveZeroFillPair[zfOpt];
    fpmfPair = resolveFpmfPair[fpmfOpt];
    phPair = resolvePhasePair[phOpt];
    tppiProc = resolveTPPIProcessing[tppiProcOpt];
    indBand = resolveIndirectBandwidth[indBandOpt];
    storeFullQ = resolveStorageMode[storageOpt];
    If[zfPair === $Failed || fpmfPair === $Failed || phPair === $Failed ||
       tppiProc === $Failed || indBand === $Failed || !BooleanQ[storeFullQ],
        Return[$Failed]
    ];

    sigProc = sig;
    If[!MatchQ[lbOpt, None | False],
        sigProc = ApplyLineBroadening[sigProc, LineBroadening -> lbOpt];
        If[sigProc === $Failed, Return[$Failed]]
    ];
    If[!MatchQ[apoOpt, None | False],
        sigProc = ApplyApodization[sigProc, Apodization -> apoOpt];
        If[sigProc === $Failed, Return[$Failed]]
    ];

    procAssoc = sigProc["Association"];
    dProc = procAssoc["Data"];
    qm = procAssoc["QuadratureMethod"];
    n1Log = procAssoc["N1"];
    n2 = procAssoc["N2"];
    dt1 = N[procAssoc["Timing1"][[3]]];
    dt2 = N[procAssoc["Timing2"][[3]]];
    {n2ZF, n1ZF} = {zfPair[[1]] * n2, zfPair[[2]] * n1Log};
    n1Out = n1ZF;
    dt1Out = dt1;

    Which[
        qm === "PhaseModulated" || qm === None,
            Module[{work = toPackedComplexMat[dProc], padded, s},
                If[fpmfPair[[1]] != 1., work[[All, 1]] *= fpmfPair[[1]]];
                If[fpmfPair[[2]] != 1., work[[1]] *= fpmfPair[[2]]];
                padded = If[n1ZF > n1Log || n2ZF > n2,
                    toPackedComplexMat[PadRight[work, {n1ZF, n2ZF}, 0. + 0. I]],
                    work
                ];
                s = toPackedComplexMat[Fourier[padded, FourierParameters -> {1, -1}] / N[n1Log * n2]];
                s = toPackedComplexMat[RotateRight[s, {Ceiling[n1ZF / 2] - 1, Ceiling[n2ZF / 2] - 1}]];
                s = apply1DLinearPhase[s, phPair[[1, 1]], phPair[[1, 2]], 2];
                uMat = apply1DLinearPhase[s, phPair[[2, 1]], phPair[[2, 2]], 1];
                vMat = None;
                repType = "Complex"
            ],
        qm === "States" || qm === "StatesTPPI",
            {uMat, vMat} = ftHypercomplexPair[dProc["Cos"], dProc["Sin"], fpmfPair, zfPair, phPair, storeFullQ];
            repType = If[storeFullQ, "Hypercomplex", "ReducedComplex"],
        qm === "EchoAntiEcho",
            mCos = toPackedComplexMat[dProc["Echo"] + dProc["AntiEcho"]];
            mSin = toPackedComplexMat[(-I) * (dProc["Echo"] - dProc["AntiEcho"])];
            {uMat, vMat} = ftHypercomplexPair[mCos, mSin, fpmfPair, zfPair, phPair, storeFullQ];
            repType = If[storeFullQ, "Hypercomplex", "ReducedComplex"],
        qm === "OversampledStates",
            {uFull, vFull} = ftHypercomplexPair[dProc["Cos"], dProc["Sin"], fpmfPair, zfPair, {phPair[[1]], {0., 0.}}, storeFullQ];
            If[indBand === "Logical",
                idxRange = If[EvenQ[n1ZF],
                    (n1ZF / 2 + 1) ;; (3 * n1ZF / 2),
                    (n1ZF - Floor[n1ZF / 2]) ;; (n1ZF + Floor[n1ZF / 2])
                ];
                uMat = apply1DLinearPhase[toPackedComplexMat[uFull[[idxRange]]], phPair[[2, 1]], phPair[[2, 2]], 1];
                vMat = If[storeFullQ,
                    apply1DLinearPhase[toPackedComplexMat[vFull[[idxRange]]], phPair[[2, 1]], phPair[[2, 2]], 1],
                    None
                ];
                n1Out = n1ZF;
                dt1Out = dt1,
                uMat = apply1DLinearPhase[uFull, phPair[[2, 1]], phPair[[2, 2]], 1];
                vMat = If[storeFullQ, apply1DLinearPhase[vFull, phPair[[2, 1]], phPair[[2, 2]], 1], None];
                n1Out = Length[uFull];
                dt1Out = 0.5 * dt1
            ];
            repType = If[storeFullQ, "Hypercomplex", "ReducedComplex"],
        qm === "TPPI",
            f2c = ftDirectAxis[dProc, fpmfPair[[1]], zfPair[[1]], phPair[[1]]];
            Switch[tppiProc,
                "RealOnly",
                    {uMat, dcU, tppiHalfBinModQ} = ftTPPIRealHalf[toPackedRealMat[Re[f2c]], fpmfPair[[2]], n1Log, n1ZF, phPair[[2]]];
                    vMat = None;
                    repType = "ReducedComplex";
                    tppiDC = If[dcU =!= None, <|"U" -> dcU|>, None],
                "RetainedHalf",
                    {uMat, dcU, tppiHalfBinModQ} = ftTPPIRealHalf[toPackedRealMat[Re[f2c]], fpmfPair[[2]], n1Log, n1ZF, phPair[[2]]];
                    If[storeFullQ,
                        {vMat, dcV, tppiHalfBinModQ} = ftTPPIRealHalf[toPackedRealMat[Im[f2c]], fpmfPair[[2]], n1Log, n1ZF, phPair[[2]]];
                        repType = "TPPIHypercomplex";
                        tppiDC = If[dcU =!= None, <|"U" -> dcU, "V" -> dcV|>, None],
                        vMat = None;
                        repType = "ReducedComplex";
                        tppiDC = If[dcU =!= None, <|"U" -> dcU|>, None]
                    ],
                "Staggered",
                    cEven = toPackedComplexMat[Table[(-1.)^(m - 1) * f2c[[2 * m - 1]], {m, n1Log}]];
                    sOdd = toPackedComplexMat[Table[(-1.)^m * f2c[[2 * m]], {m, n1Log}]];
                    hdVec = Developer`ToPackedArray[Exp[-I * Pi * buildFrequencyAxisHz[n1ZF, dt1] * dt1]];
                    u0 = toPackedComplexMat[
                        ftIndirectAxis[toPackedRealMat[Re[cEven]], fpmfPair[[2]], zfPair[[2]], {0., 0.}] +
                        I * hdVec * ftIndirectAxis[toPackedRealMat[Re[sOdd]], 1.0, zfPair[[2]], {0., 0.}]
                    ];
                    uMat = apply1DLinearPhase[u0, phPair[[2, 1]], phPair[[2, 2]], 1];
                    tppiHalfBinModQ = OddQ[n1ZF];
                    If[EvenQ[n1ZF],
                        dcU = Module[{wR = toPackedRealMat[Re[f2c]]},
                            If[fpmfPair[[2]] != 1., wR[[1]] *= N[fpmfPair[[2]]]];
                            toPackedRealVec[Total[wR] / N[n1Log]]
                        ],
                        dcU = None
                    ];
                    If[storeFullQ,
                        v0 = toPackedComplexMat[
                            ftIndirectAxis[toPackedRealMat[Im[cEven]], fpmfPair[[2]], zfPair[[2]], {0., 0.}] +
                            I * hdVec * ftIndirectAxis[toPackedRealMat[Im[sOdd]], 1.0, zfPair[[2]], {0., 0.}]
                        ];
                        vMat = apply1DLinearPhase[v0, phPair[[2, 1]], phPair[[2, 2]], 1];
                        repType = "TPPIHypercomplex";
                        If[EvenQ[n1ZF],
                            dcV = Module[{wI = toPackedRealMat[Im[f2c]]},
                                If[fpmfPair[[2]] != 1., wI[[1]] *= N[fpmfPair[[2]]]];
                                toPackedRealVec[Total[wI] / N[n1Log]]
                            ];
                            tppiDC = <|"U" -> dcU, "V" -> dcV|>,
                            tppiDC = None
                        ],
                        vMat = None;
                        repType = "ReducedComplex";
                        tppiDC = If[dcU =!= None, <|"U" -> dcU|>, None]
                    ],
                "Hilbert",
                    If[zfPair[[2]] < 2, Message[Spectrum2D::tppihilbertzf]];
                    {u0, dcU, tppiHalfBinModQ} = ftTPPIRealHalf[toPackedRealMat[Re[f2c]], fpmfPair[[2]], n1Log, n1ZF, {0., 0.}];
                    uMat = apply1DLinearPhase[hilbertReconF1Real[toPackedRealMat[Re[u0]]], phPair[[2, 1]], phPair[[2, 2]], 1];
                    If[storeFullQ,
                        {v0, dcV, tppiHalfBinModQ} = ftTPPIRealHalf[toPackedRealMat[Im[f2c]], fpmfPair[[2]], n1Log, n1ZF, {0., 0.}];
                        vMat = apply1DLinearPhase[hilbertReconF1Real[toPackedRealMat[Re[v0]]], phPair[[2, 1]], phPair[[2, 2]], 1];
                        repType = "TPPIHypercomplex";
                        tppiDC = If[dcU =!= None, <|"U" -> dcU, "V" -> dcV|>, None],
                        vMat = None;
                        repType = "ReducedComplex";
                        tppiDC = If[dcU =!= None, <|"U" -> dcU|>, None]
                    ]
            ],
        True,
            Message[Spectrum2D::badarg, "Unsupported QuadratureMethod in FT: " <> ToString[qm]];
            Return[$Failed]
    ];

    freq1 = buildFrequencyAxisHz[n1Out, dt1Out];
    freq2 = buildFrequencyAxisHz[n2ZF, dt2];

    specAssoc = <|
        "Data" -> uMat,
        "ImaginaryDirectData" -> vMat,
        "Frequency1" -> freq1,
        "Frequency2" -> freq2,
        "SpectralWidth1Hz" -> 1. / dt1Out,
        "SpectralWidth2Hz" -> 1. / dt2,
        "N1" -> n1Out,
        "N2" -> n2ZF,
        "QuadratureMethod" -> qm,
        "Representation" -> repType,
        "QuadratureCapabilities" -> buildCapabilities[repType, qm, <|
            "TPPIProcessing" -> tppiProc,
            "IndirectBandwidth" -> indBand,
            "HasTPPIDCRow" -> AssociationQ[tppiDC],
            "TPPIHalfBinModulated" -> tppiHalfBinModQ,
            "TPPIF1Modified" -> False,
            "Apodization" -> resolveAxisPair[apoOpt, None]
        |>],
        "Processing" -> Join[
            <|
                "ZeroFill" -> zfPair,
                "FirstPointMultiplicationFactor" -> fpmfPair,
                "LineBroadening" -> resolveAxisPair[lbOpt, None],
                "Apodization" -> resolveAxisPair[apoOpt, None],
                "PhaseCorrection" -> phPair,
                "StorageMode" -> If[storeFullQ, "Hypercomplex", "Compact"]
            |>,
            If[qm === "TPPI", <|"TPPIProcessing" -> tppiProc, "TPPIHalfBinModulated" -> tppiHalfBinModQ|>, <||>],
            If[qm === "OversampledStates", <|"IndirectBandwidth" -> indBand|>, <||>]
        ],
        "Metadata" -> Lookup[procAssoc, "Metadata", <||>]
    |>;
    If[qm === "TPPI" && AssociationQ[tppiDC], specAssoc["TPPIDCRow"] = tppiDC];
    If[qm === "OversampledStates",
        specAssoc["LogicalN1"] = n1ZF;
        specAssoc["LogicalSpectralWidth1Hz"] = 1. / dt1
    ];
    Spectrum2D[specAssoc]
];

(* Compact formatting & accessors for Spectrum2D *)
Spectrum2D /: Format[Spectrum2D[assoc_Association]] := RowBox[{
    "Spectrum2D[",
    ToBoxes[ToString[Lookup[assoc, "N1", "?"]] <> " x " <> ToString[Lookup[assoc, "N2", "?"]]],
    ", SW1 -> ",
    ToBoxes[EngineeringForm[N[Lookup[assoc, "SpectralWidth1Hz", 0.]], 4]],
    " Hz, SW2 -> ",
    ToBoxes[EngineeringForm[N[Lookup[assoc, "SpectralWidth2Hz", 0.]], 4]],
    " Hz]"
}] // DisplayForm;

getSpectrumV[assoc_Association] := Lookup[assoc, "ImaginaryDirectData", None];

hypercomplexMagnitudeMat[assoc_Association] := With[
    {u = assoc["Data"], v = getSpectrumV[assoc]},
    If[v === None,
        toPackedRealMat[Abs[u]],
        toPackedRealMat[Sqrt[Abs[u]^2 + Abs[v]^2]]
    ]
];

Spectrum2D[assoc_Association][prop_String] := With[{v = getSpectrumV[assoc], rep = Lookup[assoc, "Representation", "Complex"]},
    Switch[
        prop,
        "Association", assoc,
        "Keys", Keys[assoc],
        "Dimensions", {assoc["N1"], assoc["N2"]},
        "U", assoc["Data"],
        "V", If[v === None, Missing["ComponentUnavailable", "V"], v],
        "Real" | "RR", toPackedRealMat[Re[assoc["Data"]]],
        "Imaginary" | "IR", toPackedRealMat[Im[assoc["Data"]]],
        "RI", If[v === None, Missing["ComponentUnavailable", "RI"], toPackedRealMat[Re[v]]],
        "II", If[v === None, Missing["ComponentUnavailable", "II"], toPackedRealMat[Im[v]]],
        "Magnitude" | "HypercomplexMagnitude", hypercomplexMagnitudeMat[assoc],
        "Components",
            If[v === None,
                <|"RR" -> toPackedRealMat[Re[assoc["Data"]]], "IR" -> toPackedRealMat[Im[assoc["Data"]]]|>,
                <|"RR" -> toPackedRealMat[Re[assoc["Data"]]], "IR" -> toPackedRealMat[Im[assoc["Data"]]],
                  "RI" -> toPackedRealMat[Re[v]], "II" -> toPackedRealMat[Im[v]]|>
            ],
        "HypercomplexData", If[v === None, <|"U" -> assoc["Data"]|>, <|"U" -> assoc["Data"], "V" -> v|>],
        "Representation", rep,
        "HypercomplexQ", MemberQ[{"Hypercomplex", "TPPIHypercomplex"}, rep] && v =!= None,
        "QuadratureCapabilities", Lookup[assoc, "QuadratureCapabilities", buildCapabilities[rep, Lookup[assoc, "QuadratureMethod", None], Join[Lookup[assoc, "Processing", <||>], <|"HasTPPIDCRow" -> KeyExistsQ[assoc, "TPPIDCRow"]|>]]],
        "QuadratureReport", QuadratureDiagnostics[Spectrum2D[assoc]],
        "AngularFrequency1", toPackedRealVec[2. * Pi * assoc["Frequency1"]],
        "AngularFrequency2", toPackedRealVec[2. * Pi * assoc["Frequency2"]],
        _, Lookup[assoc, prop, Lookup[Lookup[assoc, "Metadata", <||>], prop, Missing["KeyAbsent", prop]]]
    ]
];

Spectrum2D /: Re[Spectrum2D[assoc_Association]] := toPackedRealMat[Re[assoc["Data"]]];
Spectrum2D /: Im[Spectrum2D[assoc_Association]] := toPackedRealMat[Im[assoc["Data"]]];
Spectrum2D /: Abs[Spectrum2D[assoc_Association]] := hypercomplexMagnitudeMat[assoc];

Spectrum2D /: HoldPattern[Times[a___, spec:Spectrum2D[assoc_Association], b___]] /;
    AllTrue[{a, b}, NumericQ] := Module[{newAssoc = assoc, c = N[Times[a, b]], v = getSpectrumV[assoc]},
    newAssoc["Data"] = toPackedComplexMat[c * assoc["Data"]];
    If[v =!= None, newAssoc["ImaginaryDirectData"] = toPackedComplexMat[c * v]];
    If[KeyExistsQ[assoc, "TPPIDCRow"] && AssociationQ[assoc["TPPIDCRow"]],
        newAssoc["TPPIDCRow"] = AssociationMap[#[[1]] -> toPackedRealVec[c * #[[2]]] &, assoc["TPPIDCRow"]]
    ];
    Spectrum2D[newAssoc]
];

reflectMatAxis[mat_?MatrixQ, axis_String] := Which[
    axis === "F1" || axis === "Indirect",
        If[EvenQ[Length[mat]],
            toPackedComplexMat[RotateLeft[Reverse[mat], 1]],
            toPackedComplexMat[Reverse[mat]]
        ],
    axis === "F2" || axis === "Direct",
        If[EvenQ[Dimensions[mat][[2]]],
            toPackedComplexMat[RotateLeft[#, 1] & /@ (Reverse /@ mat)],
            toPackedComplexMat[Reverse /@ mat]
        ],
    True,
        Message[Spectrum2D::badaxis, axis]; $Failed
];

reflectVec1D[vec_List] := If[EvenQ[Length[vec]],
    toPackedRealVec[RotateLeft[Reverse[vec], 1]],
    toPackedRealVec[Reverse[vec]]
];

Spectrum2D /: ReflectSpectrum[Spectrum2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {axisOpt, uRef, v = getSpectrumV[assoc], vRef, newAssoc = assoc, rep = Lookup[assoc, "Representation", "Complex"], qm = Lookup[assoc, "QuadratureMethod", None]},
    axisOpt = Axis /. {opts} /. Options[ReflectSpectrum];
    uRef = reflectMatAxis[assoc["Data"], axisOpt];
    If[uRef === $Failed, Return[$Failed]];
    newAssoc["Data"] = uRef;
    If[v =!= None,
        vRef = reflectMatAxis[v, axisOpt];
        newAssoc["ImaginaryDirectData"] = vRef
    ];
    If[(axisOpt === "F2" || axisOpt === "Direct") && KeyExistsQ[assoc, "TPPIDCRow"] && AssociationQ[assoc["TPPIDCRow"]],
        newAssoc["TPPIDCRow"] = AssociationMap[#[[1]] -> reflectVec1D[#[[2]]] &, assoc["TPPIDCRow"]]
    ];
    If[(axisOpt === "F1" || axisOpt === "Indirect") && rep === "TPPIHypercomplex",
        newAssoc = KeyDrop[newAssoc, "TPPIDCRow"];
        newAssoc["Processing"] = Join[Lookup[newAssoc, "Processing", <||>], <|"TPPIF1Modified" -> True|>];
        newAssoc["QuadratureCapabilities"] = buildCapabilities[
            rep, qm, Join[newAssoc["Processing"], <|"HasTPPIDCRow" -> KeyExistsQ[newAssoc, "TPPIDCRow"]|>]
        ]
    ];
    Spectrum2D[newAssoc]
];

ConjugateSpectrum2D[Spectrum2D[assoc_Association], axis_String: "F1"] := Module[
    {v = getSpectrumV[assoc], rep = Lookup[assoc, "Representation", "Complex"], newAssoc = assoc},
    Which[
        axis === "F1" || axis === "Indirect",
            newAssoc["Data"] = toPackedComplexMat[Conjugate[assoc["Data"]]];
            If[v =!= None, newAssoc["ImaginaryDirectData"] = toPackedComplexMat[Conjugate[v]]];
            Spectrum2D[newAssoc],
        axis === "F2" || axis === "Direct",
            Which[
                v =!= None,
                    newAssoc["ImaginaryDirectData"] = toPackedComplexMat[-v];
                    If[KeyExistsQ[assoc, "TPPIDCRow"] && AssociationQ[assoc["TPPIDCRow"]] && KeyExistsQ[assoc["TPPIDCRow"], "V"],
                        newAssoc["TPPIDCRow"] = Join[assoc["TPPIDCRow"], <|"V" -> toPackedRealVec[-assoc["TPPIDCRow"]["V"]]|>]
                    ];
                    Spectrum2D[newAssoc],
                rep === "Complex",
                    newAssoc["Data"] = toPackedComplexMat[Conjugate[assoc["Data"]]];
                    Spectrum2D[newAssoc],
                True,
                    Message[Spectrum2D::nocomp, "V", rep];
                    $Failed
            ],
        True,
            Message[Spectrum2D::badaxis, axis];
            $Failed
    ]
];

CombineComplementarySpectra[specs:{Spectrum2D[_], Spectrum2D[_]}, opts:(_Rule|_RuleDelayed)...] :=
    CombineComplementarySpectra[specs[[1]], specs[[2]], opts];

Spectrum2D /: CombineComplementarySpectra[s1:Spectrum2D[a1_Association], s2:Spectrum2D[a2_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {rep1 = Lookup[a1, "Representation", "Complex"], rep2 = Lookup[a2, "Representation", "Complex"],
     s2Ref, a2Ref, sumEA, diffEA, uMat, vMat, newAssoc},
    If[s1["Dimensions"] =!= s2["Dimensions"] ||
       Max[Abs[N[a1["Frequency1"] - a2["Frequency1"]]]] > 10^-9 * (1. + Max[Abs[a1["Frequency1"]]]) ||
       Max[Abs[N[a1["Frequency2"] - a2["Frequency2"]]]] > 10^-9 * (1. + Max[Abs[a1["Frequency2"]]]),
        Message[SignalData2D::baddata, "Complementary Spectrum2D datasets must share dimensions and frequency axes."];
        Return[$Failed]
    ];
    s2Ref = ReflectSpectrum[s2];
    If[s2Ref === $Failed, Return[$Failed]];
    If[rep1 === "Complex" && rep2 === "Complex",
        a2Ref = First[s2Ref];
        sumEA = 0.5 * (a1["Data"] + a2Ref["Data"]);
        diffEA = 0.5 * (a1["Data"] - a2Ref["Data"]);
        uMat = toPackedComplexMat[toPackedRealMat[Re[sumEA]] + I * toPackedRealMat[Im[diffEA]]];
        vMat = toPackedComplexMat[toPackedRealMat[Im[sumEA]] - I * toPackedRealMat[Re[diffEA]]];
        newAssoc = a1;
        newAssoc["Data"] = uMat;
        newAssoc["ImaginaryDirectData"] = vMat;
        newAssoc["QuadratureMethod"] = "EchoAntiEcho";
        newAssoc["Representation"] = "Hypercomplex";
        newAssoc["QuadratureCapabilities"] = buildCapabilities["Hypercomplex", "EchoAntiEcho"];
        Spectrum2D[newAssoc],
        CombineData[{s1, s2Ref}, {0.5, 0.5}, opts]
    ]
];

CombineData[specs:{_Spectrum2D..}, opts:(_Rule|_RuleDelayed)...] :=
    CombineData[specs, ConstantArray[1. / Length[specs], Length[specs]], opts];

CombineData[specs:{_Spectrum2D..}, weights_List, opts:(_Rule|_RuleDelayed)...] := Module[
    {w = N[weights], firstA = First[specs]["Association"], allV, sameRepQ, reps, qm0, newAssoc},
    If[Length[specs] =!= Length[weights],
        Message[CombineData::incompw]; Return[$Failed]
    ];
    If[!AllTrue[Rest[specs], (
        #["Dimensions"] === First[specs]["Dimensions"] &&
        Max[Abs[N[#["Frequency1"] - firstA["Frequency1"]]]] <= 10^-9 * (1. + Max[Abs[firstA["Frequency1"]]]) &&
        Max[Abs[N[#["Frequency2"] - firstA["Frequency2"]]]] <= 10^-9 * (1. + Max[Abs[firstA["Frequency2"]]])
    ) &],
        Message[SignalData2D::baddata, "All Spectrum2D datasets must share dimensions and frequency axes."];
        Return[$Failed]
    ];
    reps = (Lookup[First[#], "Representation", "Complex"] &) /@ specs;
    allV = (getSpectrumV[First[#]] =!= None &) /@ specs;
    sameRepQ = (Length[DeleteDuplicates[reps]] === 1);
    newAssoc = firstA;
    newAssoc["Data"] = toPackedComplexMat[Sum[w[[i]] * specs[[i]]["Data"], {i, Length[specs]}]];
    If[AllTrue[allV, TrueQ],
        newAssoc["ImaginaryDirectData"] = toPackedComplexMat[Sum[w[[i]] * getSpectrumV[First[specs[[i]]]], {i, Length[specs]}]];
        newAssoc["Representation"] = If[sameRepQ, First[reps], "Hypercomplex"],
        If[MemberQ[allV, True], Message[Spectrum2D::repmerge]];
        newAssoc["ImaginaryDirectData"] = None;
        newAssoc["Representation"] = If[sameRepQ && First[reps] === "Complex", "Complex", "ReducedComplex"]
    ];
    qm0 = Lookup[newAssoc, "QuadratureMethod", None];
    If[AllTrue[specs, KeyExistsQ[First[#], "TPPIDCRow"] && AssociationQ[First[#]["TPPIDCRow"]] &] &&
       AllTrue[allV, TrueQ] && KeyExistsQ[firstA, "TPPIDCRow"],
        newAssoc["TPPIDCRow"] = AssociationMap[
            Function[{kv}, With[{k = kv[[1]]}, k -> toPackedRealVec[Sum[w[[i]] * First[specs[[i]]]["TPPIDCRow"][k], {i, Length[specs]}]]]],
            firstA["TPPIDCRow"]
        ],
        If[KeyExistsQ[newAssoc, "TPPIDCRow"], newAssoc = KeyDrop[newAssoc, "TPPIDCRow"]]
    ];
    newAssoc["QuadratureCapabilities"] = buildCapabilities[
        newAssoc["Representation"], qm0,
        Join[Lookup[newAssoc, "Processing", <||>], <|"HasTPPIDCRow" -> KeyExistsQ[newAssoc, "TPPIDCRow"]|>]
    ];
    Spectrum2D[newAssoc]
];

Spectrum2D /: Expand[Spectrum2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {embedQ, compOpt, uOpt, scale, f1, f2, zMat, n1, n2},
    embedQ = TrueQ[EmbedCoordinates /. {opts} /. {EmbedCoordinates -> True}];
    compOpt = Component /. {opts} /. {Component -> "Complex"};
    uOpt = AxisUnits /. {opts} /. {AxisUnits -> "Hz"};
    zMat = selectSpectrumComponentMat[Spectrum2D[assoc], compOpt];
    If[zMat === $Failed, Return[$Failed]];
    If[!embedQ, Return[zMat]];
    scale = If[uOpt === "rad/s", 2. * Pi, 1.];
    f1 = scale * assoc["Frequency1"];
    f2 = scale * assoc["Frequency2"];
    {n1, n2} = Dimensions[zMat];
    Flatten[Table[{f1[[i1]], f2[[i2]], zMat[[i1, i2]]}, {i1, n1}, {i2, n2}], 1]
];

PhaseSpectrum2D[Spectrum2D[assoc_Association], phaseSpec_] := Module[
    {phPair = resolvePhasePair[phaseSpec], u = assoc["Data"], v = getSpectrumV[assoc],
     rep = Lookup[assoc, "Representation", "Complex"], qm = Lookup[assoc, "QuadratureMethod", None],
     hasF2, hasF1, newAssoc = assoc, prevPh},
    If[phPair === $Failed, Return[$Failed]];
    hasF2 = (phPair[[1]] != {0., 0.});
    hasF1 = (phPair[[2]] != {0., 0.});
    Which[
        rep === "Hypercomplex",
            If[v === None, Message[Spectrum2D::nof2phase, rep]; Return[$Failed]];
            If[hasF2, {u, v} = phaseUVF2[u, v, phPair[[1, 1]], phPair[[1, 2]]]];
            If[hasF1,
                u = apply1DLinearPhase[u, phPair[[2, 1]], phPair[[2, 2]], 1];
                v = apply1DLinearPhase[v, phPair[[2, 1]], phPair[[2, 2]], 1]
            ],
        rep === "TPPIHypercomplex",
            If[v === None, Message[Spectrum2D::nof2phase, rep]; Return[$Failed]];
            If[hasF2,
                {u, v} = phaseUVF2[u, v, phPair[[1, 1]], phPair[[1, 2]]];
                If[KeyExistsQ[assoc, "TPPIDCRow"] && AssociationQ[assoc["TPPIDCRow"]] &&
                   KeyExistsQ[assoc["TPPIDCRow"], "U"] && KeyExistsQ[assoc["TPPIDCRow"], "V"],
                    With[{dcRes = phaseUVF2[{assoc["TPPIDCRow"]["U"]}, {assoc["TPPIDCRow"]["V"]}, phPair[[1, 1]], phPair[[1, 2]]]},
                        newAssoc["TPPIDCRow"] = <|
                            "U" -> toPackedRealVec[Re[dcRes[[1, 1]]]],
                            "V" -> toPackedRealVec[Re[dcRes[[2, 1]]]]
                        |>
                    ]
                ]
            ];
            If[hasF1,
                Message[Spectrum2D::tppif1phase];
                u = apply1DLinearPhase[u, phPair[[2, 1]], phPair[[2, 2]], 1];
                v = apply1DLinearPhase[v, phPair[[2, 1]], phPair[[2, 2]], 1];
                newAssoc = KeyDrop[newAssoc, "TPPIDCRow"];
                newAssoc["Processing"] = Join[Lookup[newAssoc, "Processing", <||>], <|"TPPIF1Modified" -> True|>]
            ],
        rep === "ReducedComplex",
            If[hasF2,
                Message[Spectrum2D::nof2phase, rep];
                Return[$Failed]
            ];
            If[hasF1,
                If[qm === "TPPI", Message[Spectrum2D::tppif1phase]];
                u = apply1DLinearPhase[u, phPair[[2, 1]], phPair[[2, 2]], 1]
            ],
        rep === "Complex",
            If[hasF1 && hasF2, Message[Spectrum2D::pmphase]];
            u = apply1DLinearPhase[u, phPair[[1, 1]], phPair[[1, 2]], 2];
            u = apply1DLinearPhase[u, phPair[[2, 1]], phPair[[2, 2]], 1],
        True,
            Message[Spectrum2D::badarg, "Unknown Spectrum2D Representation: " <> ToString[rep]];
            Return[$Failed]
    ];
    newAssoc["Data"] = toPackedComplexMat[u];
    If[v =!= None, newAssoc["ImaginaryDirectData"] = toPackedComplexMat[v]];
    prevPh = Lookup[Lookup[assoc, "Processing", <||>], "PhaseCorrection", {{0., 0.}, {0., 0.}}];
    newAssoc["Processing"] = Join[
        Lookup[newAssoc, "Processing", <||>],
        <|"PhaseCorrection" -> prevPh + phPair|>
    ];
    newAssoc["QuadratureCapabilities"] = buildCapabilities[
        rep, qm, Join[newAssoc["Processing"], <|"HasTPPIDCRow" -> KeyExistsQ[newAssoc, "TPPIDCRow"]|>]
    ];
    Spectrum2D[newAssoc]
];

selectSpectrumComponentMat[spec:Spectrum2D[assoc_Association], comp_] := With[
    {rep = Lookup[assoc, "Representation", "Complex"], v = getSpectrumV[assoc]},
    Switch[
        comp,
        "Real" | "RR", toPackedRealMat[Re[assoc["Data"]]],
        "Imaginary" | "IR", toPackedRealMat[Im[assoc["Data"]]],
        "RI",
            If[v === None,
                Message[Spectrum2D::nocomp, "RI", rep]; $Failed,
                toPackedRealMat[Re[v]]
            ],
        "II",
            If[v === None,
                Message[Spectrum2D::nocomp, "II", rep]; $Failed,
                toPackedRealMat[Im[v]]
            ],
        "Magnitude" | "HypercomplexMagnitude", hypercomplexMagnitudeMat[assoc],
        "Complex", assoc["Data"],
        _, Message[Spectrum2D::badarg, "Unknown Component -> " <> ToString[comp]]; $Failed
    ]
];

ExtractSlice[spec:Spectrum2D[assoc_Association], axis_String, pos_?NumericQ, opts:(_Rule|_RuleDelayed)...] := Module[
    {posType, units, comp, mat, scale, f1, f2, idx, targetAxis, vals},
    posType = PositionType /. {opts} /. Options[ExtractSlice];
    units = AxisUnits /. {opts} /. Options[ExtractSlice];
    comp = Component /. {opts} /. Options[ExtractSlice];
    mat = selectSpectrumComponentMat[spec, comp];
    If[mat === $Failed, Return[$Failed]];
    scale = If[units === "rad/s", 2. * Pi, 1.];
    f1 = scale * assoc["Frequency1"];
    f2 = scale * assoc["Frequency2"];
    Which[
        axis === "F2" || axis === "Direct",
            (* Slice along F2 at fixed F1 = pos *)
            idx = If[posType === "Index",
                Clip[Round[pos], {1, Length[f1]}],
                First[Ordering[Abs[f1 - N[pos]], 1]]
            ];
            vals = mat[[idx, All]];
            Transpose[{f2, vals}],
        axis === "F1" || axis === "Indirect",
            (* Slice along F1 at fixed F2 = pos *)
            idx = If[posType === "Index",
                Clip[Round[pos], {1, Length[f2]}],
                First[Ordering[Abs[f2 - N[pos]], 1]]
            ];
            vals = mat[[All, idx]];
            Transpose[{f1, vals}],
        True,
            Message[Spectrum2D::badaxis, axis];
            $Failed
    ]
];

ProjectSpectrum2D[spec:Spectrum2D[assoc_Association], axis_String, opts:(_Rule|_RuleDelayed)...] := Module[
    {meth, units, comp, mat, scale, fAxis, reduceAxis, vals},
    meth = Method /. {opts} /. Options[ProjectSpectrum2D];
    units = AxisUnits /. {opts} /. Options[ProjectSpectrum2D];
    comp = Component /. {opts} /. Options[ProjectSpectrum2D];
    mat = selectSpectrumComponentMat[spec, comp];
    If[mat === $Failed, Return[$Failed]];
    scale = If[units === "rad/s", 2. * Pi, 1.];
    Which[
        axis === "F2" || axis === "Direct",
            fAxis = scale * assoc["Frequency2"];
            mat = Transpose[mat],
        axis === "F1" || axis === "Indirect",
            fAxis = scale * assoc["Frequency1"];
            Null,
        True,
            Message[Spectrum2D::badaxis, axis];
            Return[$Failed]
    ];
    vals = Switch[
        meth,
        "Sum", Total /@ mat,
        "Max", Max /@ mat,
        "Skyline", (With[{row = #}, row[[First[Ordering[Abs[row], -1]]]]] &) /@ mat,
        _, Message[Spectrum2D::badarg, "Unknown Method -> " <> ToString[meth]]; Return[$Failed]
    ];
    Transpose[{fAxis, toPackedRealVec[vals]}]
];

PlotSpectrum2D[spec:Spectrum2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {compOpt, unitsOpt, revOpt, nContours, baseFrac, cFactor, posCol, negCol, showProj,
     pRange, passOpts, effComp, zMat, scale, f1, f2, unitLabel, maxAbs, posLevels,
     negLevels, levels, revDir, revInd, sf, mainPlot, topProj, sideProj, pTop, pSide},

    compOpt = Component /. {opts} /. Options[PlotSpectrum2D];
    unitsOpt = AxisUnits /. {opts} /. Options[PlotSpectrum2D];
    revOpt = ReverseAxes /. {opts} /. Options[PlotSpectrum2D];
    nContours = Max[1, Round[Contours /. {opts} /. Options[PlotSpectrum2D]]];
    baseFrac = N[BaseContourLevel /. {opts} /. Options[PlotSpectrum2D]];
    cFactor = N[ContourFactor /. {opts} /. Options[PlotSpectrum2D]];
    posCol = PositiveColor /. {opts} /. Options[PlotSpectrum2D];
    negCol = NegativeColor /. {opts} /. Options[PlotSpectrum2D];
    showProj = TrueQ[ShowProjections /. {opts} /. Options[PlotSpectrum2D]];
    pRange = PlotRange /. {opts} /. Options[PlotSpectrum2D];
    If[MatchQ[pRange, {{_, _}, {_, _}}], pRange = Append[pRange, All]];
    passOpts = FilterRules[{opts}, Options[ListContourPlot]];

    effComp = If[compOpt === Automatic,
        If[assoc["QuadratureMethod"] === "PhaseModulated", "Magnitude", "Real"],
        compOpt
    ];
    zMat = selectSpectrumComponentMat[spec, effComp];
    If[zMat === $Failed, Return[$Failed]];
    If[!Developer`PackedArrayQ[zMat, Real], zMat = toPackedRealMat[zMat]];

    scale = If[unitsOpt === "rad/s", 2. * Pi, 1.];
    unitLabel = If[unitsOpt === "rad/s", "rad s^-1", "Hz"];
    f1 = scale * assoc["Frequency1"];
    f2 = scale * assoc["Frequency2"];

    If[MatchQ[revOpt, {_, _}],
        {revDir, revInd} = TrueQ /@ revOpt,
        {revDir, revInd} = {TrueQ[revOpt], TrueQ[revOpt]}
    ];
    sf = {If[revDir, "Reverse", Identity], If[revInd, "Reverse", Identity]};

    maxAbs = Max[Abs[zMat]];
    mainPlot = If[maxAbs == 0.,
        ListContourPlot[
            {{0., 0.}, {0., 0.}},
            DataRange -> {{First[f2], Last[f2]}, {First[f1], Last[f1]}},
            Contours -> 0,
            ContourShading -> None,
            ScalingFunctions -> sf,
            Frame -> True,
            FrameLabel -> {"F2 (" <> unitLabel <> ")", "F1 (" <> unitLabel <> ")"},
            PlotRange -> If[pRange === All, {{First[f2], Last[f2]}, {First[f1], Last[f1]}, All}, pRange],
            Sequence @@ passOpts
        ],
        posLevels = Table[{baseFrac * maxAbs * (cFactor^k), Directive[posCol]}, {k, 0, nContours - 1}];
        negLevels = If[effComp === "Magnitude",
            {},
            Table[{-baseFrac * maxAbs * (cFactor^k), Directive[negCol]}, {k, 0, nContours - 1}]
        ];
        levels = Join[Reverse[negLevels], posLevels];
        ListContourPlot[
            zMat,
            DataRange -> {{First[f2], Last[f2]}, {First[f1], Last[f1]}},
            Contours -> levels,
            ContourShading -> None,
            ScalingFunctions -> sf,
            Frame -> True,
            FrameLabel -> {"F2 (" <> unitLabel <> ")", "F1 (" <> unitLabel <> ")"},
            PlotRange -> pRange,
            Sequence @@ passOpts
        ]
    ];

    If[!showProj, Return[mainPlot]];

    topProj = ProjectSpectrum2D[spec, "F2", AxisUnits -> unitsOpt, Component -> effComp, Method -> "Skyline"];
    sideProj = ProjectSpectrum2D[spec, "F1", AxisUnits -> unitsOpt, Component -> effComp, Method -> "Skyline"];
    pTop = ListLinePlot[topProj,
        Frame -> True, Axes -> False,
        FrameTicks -> {{None, None}, {None, None}},
        PlotRange -> All,
        ScalingFunctions -> {If[revDir, "Reverse", Identity], Identity},
        PlotStyle -> Directive[posCol, AbsoluteThickness[1.2]],
        ImageSize -> {320, 75}
    ];
    pSide = ListLinePlot[Reverse /@ sideProj,
        Frame -> True, Axes -> False,
        FrameTicks -> {{None, None}, {None, None}},
        PlotRange -> All,
        ScalingFunctions -> {Identity, If[revInd, "Reverse", Identity]},
        PlotStyle -> Directive[posCol, AbsoluteThickness[1.2]],
        ImageSize -> {75, 320}
    ];
    Grid[{
        {pTop, SpanFromLeft},
        {Show[mainPlot, ImageSize -> {320, 320}], pSide}
    }, Alignment -> Center, Spacings -> {0.4, 0.4}]
];

Spectrum2D /: ListContourPlot[spec_Spectrum2D, opts:(_Rule|_RuleDelayed)...] :=
    PlotSpectrum2D[spec, opts];

(* Ensure SignalData2D / Spectrum2D rules precede Signal.m catch-all Except[Signal[___]] rule *)
DownValues[CombineData] = Join[
    Select[DownValues[CombineData], !FreeQ[#, SignalData2D | Spectrum2D] &],
    Select[DownValues[CombineData], FreeQ[#, SignalData2D | Spectrum2D] &]
];
DownValues[CombineComplementarySpectra] = Join[
    Select[DownValues[CombineComplementarySpectra], !FreeQ[#, Spectrum2D] &],
    Select[DownValues[CombineComplementarySpectra], FreeQ[#, Spectrum2D] &]
];

End[];
EndPackage[];
