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
    PositionType
];

SignalData2D::usage = "SignalData2D[<|...|>] represents a simulated 2D NMR time-domain dataset with rows indexed by indirect time t1 (length N1) and columns indexed by direct time t2 (length N2). Supports property queries sig[\"prop\"], Re, Im, Abs, Times, Expand, ApplyLineBroadening, ApplyApodization, CombineData, and FT.";
Spectrum2D::usage = "Spectrum2D[<|...|>] represents a processed 2D NMR frequency-domain spectrum with rows indexed by F1 (indirect) and columns by F2 (direct). Supports property queries spec[\"prop\"], Re, Im, Abs, Times, Expand, ReflectSpectrum, CombineComplementarySpectra, CombineData, PhaseSpectrum2D, ExtractSlice, ProjectSpectrum2D, PlotSpectrum2D, and ListContourPlot.";
ApplyApodization::usage = "ApplyApodization[sig2D, opts] applies a 2D window function to a SignalData2D object (or 2D matrix) using Apodization -> spec or Apodization -> {directSpec, indirectSpec}.";
Apodization::usage = "Apodization is an option for ApplyApodization and FT specifying a window function or {directSpec, indirectSpec} pair. Supported specs include None, Automatic, fwhh, {\"Lorentzian\", fwhh}, {\"Gaussian\", fwhh}, {\"Gaussian\", lb, gb}, \"SineBell\", {\"SineBell\", shift}, \"SquaredSineBell\", and {\"SquaredSineBell\", shift}.";
PhaseSpectrum2D::usage = "PhaseSpectrum2D[spec2D, phaseSpec] applies zero- and first-order phase corrections (in radians) to a Spectrum2D object. phaseSpec may be phi0, {phi0Direct, phi0Indirect}, or {{phi0Direct, phi1Direct}, {phi0Indirect, phi1Indirect}}.";
PhaseCorrection::usage = "PhaseCorrection is an option for FT on SignalData2D specifying phase corrections in radians as phi0, {phi0Direct, phi0Indirect}, or {{phi0Direct, phi1Direct}, {phi0Indirect, phi1Indirect}}.";
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
Spectrum2D::inval = "Invalid Spectrum2D specification: `1`.";
Spectrum2D::badaxis = "Unknown axis specification `1`; expected \"F1\", \"F2\", \"Indirect\", or \"Direct\".";
Spectrum2D::badarg = "Invalid argument for 2D spectrum operation: `1`.";
ApplyApodization::badwin = "Unrecognized apodization window specification `1`.";

Options[ApplyApodization] = {Apodization -> None};
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

SignalData2D[assoc_Association][prop_String] := Switch[
    prop,
    "Association", assoc,
    "Keys", Keys[assoc],
    "Dimensions", {assoc["N1"], assoc["N2"]},
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

SignalData2D /: Expand[SignalData2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {embedQ, compOpt, mat, t1, t2, n1, n2},
    embedQ = TrueQ[EmbedCoordinates /. {opts} /. {EmbedCoordinates -> True}];
    compOpt = Component /. {opts} /. {Component -> Automatic};
    mat = extractMainSignalMatrix[assoc, compOpt];
    If[!embedQ, Return[mat]];
    t1 = If[assoc["QuadratureMethod"] === "TPPI" && KeyExistsQ[assoc, "TPPITime1"],
        assoc["TPPITime1"],
        assoc["Time1"]
    ];
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
    t1Vec = If[assoc["QuadratureMethod"] === "TPPI" && KeyExistsQ[assoc, "TPPITime1"],
        assoc["TPPITime1"], assoc["Time1"]
    ];
    t1Timing = If[assoc["QuadratureMethod"] === "TPPI" && KeyExistsQ[assoc, "TPPITiming1"],
        assoc["TPPITiming1"], assoc["Timing1"]
    ];
    w1 = buildWindow1D[t1Vec, t1Timing, indSpec];
    w2 = buildWindow1D[assoc["Time2"], assoc["Timing2"], dirSpec];
    If[w1 === $Failed || w2 === $Failed, Return[$Failed]];
    newAssoc = assoc;
    newAssoc["Data"] = apply2DWindowToData[assoc["Data"], w1, w2];
    newAssoc["Peaks2D"] = None;
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
    spec = toPackedComplexMat[RotateRight[spec, {0, Floor[n2ZF / 2] - 1}]];
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
    specT = toPackedComplexMat[RotateRight[Transpose[specT], {Floor[n1ZF / 2] - 1, 0}]];
    apply1DLinearPhase[specT, phcInd[[1]], phcInd[[2]], 1]
];

buildFrequencyAxisHz[nPts_Integer, dt_?NumericQ] := Module[{swHz = 1. / N[dt]},
    If[EvenQ[nPts],
        toPackedRealVec[Range[-(nPts / 2) + 1, nPts / 2] * (swHz / N[nPts])],
        toPackedRealVec[Range[-Floor[nPts / 2], Floor[nPts / 2]] * (swHz / N[nPts])]
    ]
];

SignalData2D /: FT[sig:SignalData2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {lbOpt, apoOpt, zfOpt, fpmfOpt, phOpt, zfPair, fpmfPair, phPair,
     sigProc, procAssoc, dProc, qm, n1Log, n2, dt1, dt2, n1ZF, n2ZF,
     specMat, f2c, f2s, hMat, mCos, mSin, f2Real, n1TPPI, n1TPPIZF, rawTPPI,
     freq1, freq2},

    lbOpt = LineBroadening /. {opts} /. Options[FT];
    apoOpt = Apodization /. {opts} /. Options[ApplyApodization];
    zfOpt = ZeroFill /. {opts} /. Options[FT];
    fpmfOpt = FirstPointMultiplicationFactor /. {opts} /. Options[FT];
    phOpt = PhaseCorrection /. {opts} /. {PhaseCorrection -> {{0., 0.}, {0., 0.}}};

    zfPair = resolveZeroFillPair[zfOpt];
    fpmfPair = resolveFpmfPair[fpmfOpt];
    phPair = resolvePhasePair[phOpt];
    If[zfPair === $Failed || fpmfPair === $Failed || phPair === $Failed, Return[$Failed]];

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

    specMat = Which[
        qm === "PhaseModulated" || qm === None,
            Module[{work = toPackedComplexMat[dProc], padded, s},
                If[fpmfPair[[1]] != 1., work[[All, 1]] *= fpmfPair[[1]]];
                If[fpmfPair[[2]] != 1., work[[1]] *= fpmfPair[[2]]];
                padded = If[n1ZF > n1Log || n2ZF > n2,
                    toPackedComplexMat[PadRight[work, {n1ZF, n2ZF}, 0. + 0. I]],
                    work
                ];
                s = toPackedComplexMat[Fourier[padded, FourierParameters -> {1, -1}] / N[n1Log * n2]];
                s = toPackedComplexMat[RotateRight[s, {Floor[n1ZF / 2] - 1, Floor[n2ZF / 2] - 1}]];
                s = apply1DLinearPhase[s, phPair[[1, 1]], phPair[[1, 2]], 2];
                apply1DLinearPhase[s, phPair[[2, 1]], phPair[[2, 2]], 1]
            ],
        qm === "States" || qm === "StatesTPPI",
            f2c = ftDirectAxis[dProc["Cos"], fpmfPair[[1]], zfPair[[1]], phPair[[1]]];
            f2s = ftDirectAxis[dProc["Sin"], fpmfPair[[1]], zfPair[[1]], phPair[[1]]];
            hMat = toPackedComplexMat[Re[f2c] + I * Re[f2s]];
            ftIndirectAxis[hMat, fpmfPair[[2]], zfPair[[2]], phPair[[2]]],
        qm === "EchoAntiEcho",
            mCos = toPackedComplexMat[dProc["Echo"] + dProc["AntiEcho"]];
            mSin = toPackedComplexMat[(-I) * (dProc["Echo"] - dProc["AntiEcho"])];
            f2c = ftDirectAxis[mCos, fpmfPair[[1]], zfPair[[1]], phPair[[1]]];
            f2s = ftDirectAxis[mSin, fpmfPair[[1]], zfPair[[1]], phPair[[1]]];
            hMat = toPackedComplexMat[Re[f2c] + I * Re[f2s]];
            ftIndirectAxis[hMat, fpmfPair[[2]], zfPair[[2]], phPair[[2]]],
        qm === "TPPI",
            n1TPPI = Length[dProc];
            n1TPPIZF = 2 * n1ZF;
            f2c = ftDirectAxis[dProc, fpmfPair[[1]], zfPair[[1]], phPair[[1]]];
            f2Real = toPackedComplexMat[Re[f2c] + 0. I];
            If[fpmfPair[[2]] != 1., f2Real[[1]] *= fpmfPair[[2]]];
            If[n1TPPIZF > n1TPPI,
                f2Real = toPackedComplexMat[PadRight[f2Real, {n1TPPIZF, n2ZF}, 0. + 0. I]]
            ];
            (* Real TPPI splits cosine power across +/- Omega; dividing by n1Log = (n1TPPI/2) matches States amplitude *)
            rawTPPI = Transpose[Developer`ToPackedArray[Fourier[#, FourierParameters -> {1, -1}] & /@ Transpose[f2Real]]] / N[n1Log];
            specMat = toPackedComplexMat[rawTPPI[[2 ;; n1ZF + 1]]];
            apply1DLinearPhase[specMat, phPair[[2, 1]], phPair[[2, 2]], 1],
        True,
            Message[Spectrum2D::badarg, "Unsupported QuadratureMethod in FT: " <> ToString[qm]];
            Return[$Failed]
    ];

    freq1 = buildFrequencyAxisHz[n1ZF, dt1];
    freq2 = buildFrequencyAxisHz[n2ZF, dt2];

    Spectrum2D[<|
        "Data" -> specMat,
        "Frequency1" -> freq1,
        "Frequency2" -> freq2,
        "SpectralWidth1Hz" -> 1. / dt1,
        "SpectralWidth2Hz" -> 1. / dt2,
        "N1" -> n1ZF,
        "N2" -> n2ZF,
        "QuadratureMethod" -> qm,
        "Processing" -> <|
            "ZeroFill" -> zfPair,
            "FirstPointMultiplicationFactor" -> fpmfPair,
            "LineBroadening" -> resolveAxisPair[lbOpt, None],
            "Apodization" -> resolveAxisPair[apoOpt, None],
            "PhaseCorrection" -> phPair
        |>,
        "Metadata" -> Lookup[procAssoc, "Metadata", <||>]
    |>]
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

Spectrum2D[assoc_Association][prop_String] := Switch[
    prop,
    "Association", assoc,
    "Keys", Keys[assoc],
    "Dimensions", {assoc["N1"], assoc["N2"]},
    "Real", toPackedRealMat[Re[assoc["Data"]]],
    "Imaginary", toPackedRealMat[Im[assoc["Data"]]],
    "Magnitude", toPackedRealMat[Abs[assoc["Data"]]],
    "AngularFrequency1", toPackedRealVec[2. * Pi * assoc["Frequency1"]],
    "AngularFrequency2", toPackedRealVec[2. * Pi * assoc["Frequency2"]],
    _, Lookup[assoc, prop, Lookup[Lookup[assoc, "Metadata", <||>], prop, Missing["KeyAbsent", prop]]]
];

Spectrum2D /: Re[Spectrum2D[assoc_Association]] := toPackedRealMat[Re[assoc["Data"]]];
Spectrum2D /: Im[Spectrum2D[assoc_Association]] := toPackedRealMat[Im[assoc["Data"]]];
Spectrum2D /: Abs[Spectrum2D[assoc_Association]] := toPackedRealMat[Abs[assoc["Data"]]];

Spectrum2D /: HoldPattern[Times[a___, spec:Spectrum2D[assoc_Association], b___]] /;
    AllTrue[{a, b}, NumericQ] := Module[{newAssoc = assoc, c = N[Times[a, b]]},
    newAssoc["Data"] = toPackedComplexMat[c * assoc["Data"]];
    Spectrum2D[newAssoc]
];

Spectrum2D /: ReflectSpectrum[Spectrum2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {newAssoc = assoc, mat = assoc["Data"]},
    If[OddQ[Length[mat]],
        Message[Spectrum2D::badarg, "ReflectSpectrum requires an even number of indirect-dimension points."];
        Return[$Failed]
    ];
    newAssoc["Data"] = toPackedComplexMat[RotateLeft[Reverse[mat], 1]];
    Spectrum2D[newAssoc]
];

CombineComplementarySpectra[specs:{Spectrum2D[_], Spectrum2D[_]}, opts:(_Rule|_RuleDelayed)...] :=
    CombineComplementarySpectra[specs[[1]], specs[[2]], opts];

Spectrum2D /: CombineComplementarySpectra[s1:Spectrum2D[_], s2:Spectrum2D[_], opts:(_Rule|_RuleDelayed)...] := Module[
    {s2Ref = ReflectSpectrum[s2]},
    If[s2Ref === $Failed, Return[$Failed]];
    CombineData[{s1, s2Ref}, {0.5, 0.5}, opts]
];

CombineData[specs:{_Spectrum2D..}, opts:(_Rule|_RuleDelayed)...] :=
    CombineData[specs, ConstantArray[1. / Length[specs], Length[specs]], opts];

CombineData[specs:{_Spectrum2D..}, weights_List, opts:(_Rule|_RuleDelayed)...] := Module[
    {w = N[weights], firstA = First[specs]["Association"], newAssoc},
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
    newAssoc = firstA;
    newAssoc["Data"] = toPackedComplexMat[Sum[w[[i]] * specs[[i]]["Data"], {i, Length[specs]}]];
    Spectrum2D[newAssoc]
];

Spectrum2D /: Expand[Spectrum2D[assoc_Association], opts:(_Rule|_RuleDelayed)...] := Module[
    {embedQ, compOpt, uOpt, scale, f1, f2, zMat, n1, n2},
    embedQ = TrueQ[EmbedCoordinates /. {opts} /. {EmbedCoordinates -> True}];
    compOpt = Component /. {opts} /. {Component -> "Complex"};
    uOpt = AxisUnits /. {opts} /. {AxisUnits -> "Hz"};
    zMat = Switch[compOpt,
        "Real", Re[Spectrum2D[assoc]],
        "Imaginary", Im[Spectrum2D[assoc]],
        "Magnitude", Abs[Spectrum2D[assoc]],
        _, assoc["Data"]
    ];
    If[!embedQ, Return[zMat]];
    scale = If[uOpt === "rad/s", 2. * Pi, 1.];
    f1 = scale * assoc["Frequency1"];
    f2 = scale * assoc["Frequency2"];
    {n1, n2} = Dimensions[zMat];
    Flatten[Table[{f1[[i1]], f2[[i2]], zMat[[i1, i2]]}, {i1, n1}, {i2, n2}], 1]
];

PhaseSpectrum2D[Spectrum2D[assoc_Association], phaseSpec_] := Module[
    {phPair = resolvePhasePair[phaseSpec], mat = assoc["Data"], newAssoc = assoc},
    If[phPair === $Failed, Return[$Failed]];
    mat = apply1DLinearPhase[mat, phPair[[1, 1]], phPair[[1, 2]], 2];
    mat = apply1DLinearPhase[mat, phPair[[2, 1]], phPair[[2, 2]], 1];
    newAssoc["Data"] = toPackedComplexMat[mat];
    Spectrum2D[newAssoc]
];

selectSpectrumComponentMat[spec_Spectrum2D, comp_] := Switch[
    comp,
    "Real", spec["Real"],
    "Imaginary", spec["Imaginary"],
    "Magnitude", spec["Magnitude"],
    "Complex", spec["Data"],
    _, Message[Spectrum2D::badarg, "Unknown Component -> " <> ToString[comp]]; $Failed
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
