(* ::Package:: *)

(* SpinDynamica Experimental Bruker 1D Spectrum Data Importer & Signal Integration *)

BeginPackage["Signals`Bruker`", {
    "Global`",
    "SpinDynamica`",
    "General`Miscellaneous`",
    "Signals`Signal`"
}];

ImportBrukerSpectrum::usage = "ImportBrukerSpectrum[path, opts] imports a processed 1D Bruker NMR spectrum from a Bruker experiment directory, a pdata directory, a processed data directory (e.g. pdata/1), or a 1r file, returning a BrukerSpectrum object (or an Association / {x, y} pairs via ReturnFormat).";
BrukerSpectrum::usage = "BrukerSpectrum[<|...|>] represents an imported experimental Bruker NMR spectrum with calibrated chemical-shift (ppm) and frequency (Hz, rad/s) axes, scaled real/imaginary data, scaling metadata, and full processing/acquisition parameters. Supports property queries spec[\"prop\"], Expand, Re, Im, Abs, Times, Length, CombineData, ListPlot, PlotBrukerSpectrum, and ToSignal.";
ReadBrukerParameters::usage = "ReadBrukerParameters[file] parses a Bruker JCAMP-DX parameter file (such as procs, acqus, proc, or acqu) into a Mathematica Association.";
PlotBrukerSpectrum::usage = "PlotBrukerSpectrum[spec, opts] plots a BrukerSpectrum using ListPlot with standard NMR conventions (reversed chemical shift ppm axis by default, or ascending Hz / rad/s frequency axis).";
ToSignal::usage = "ToSignal[spec] converts an imported 1D BrukerSpectrum into a SpinDynamica time-domain Signal object via inverse FFT on the SpinDynamica frequency grid.";

ProcessingNumber::usage = "ProcessingNumber is an option for ImportBrukerSpectrum selecting the pdata subdirectory number (Automatic selects 1 if present or the smallest integer processing directory).";
ReadImaginary::usage = "ReadImaginary is an option for ImportBrukerSpectrum specifying whether to read the imaginary 1i component (Automatic reads if present, True requires 1i, False skips 1i).";
IntensityScaling::usage = "IntensityScaling is an option for ImportBrukerSpectrum controlling spectral intensity scaling: \"Bruker\" (default, applies 2^NC_proc), \"Raw\" (unscaled stored values), \"Normalized\" (applies 2^NC_proc then scales peak magnitude to 1), or a positive numeric multiplier applied after 2^NC_proc.";
IncludeRawData::usage = "IncludeRawData is an option for ImportBrukerSpectrum controlling whether unscaled RawReal/RawImaginary arrays are retained alongside scaled arrays (default True).";
ReadAcquisitionParameters::usage = "ReadAcquisitionParameters is an option for ImportBrukerSpectrum controlling whether acqus parameters from the parent experiment directory are parsed when available (default True).";
ReturnFormat::usage = "ReturnFormat is an option for ImportBrukerSpectrum selecting \"BrukerSpectrum\" (default), \"Association\", or \"Pairs\".";
AxisUnits::usage = "AxisUnits is an option for PlotBrukerSpectrum, Expand, and ImportBrukerSpectrum (with ReturnFormat -> \"Pairs\") selecting \"ppm\" (default), \"Hz\", or \"rad/s\".";
Component::usage = "Component is an option for PlotBrukerSpectrum and Expand selecting \"Real\" (default), \"Imaginary\", \"Both\", \"Complex\", or \"Magnitude\".";
ReverseAxis::usage = "ReverseAxis is an option for PlotBrukerSpectrum controlling whether the horizontal axis decreases left-to-right (Automatic reverses \"ppm\" and keeps \"Hz\"/\"rad/s\" ascending).";
AxisOrder::usage = "AxisOrder is an option for Expand and ImportBrukerSpectrum (ReturnFormat -> \"Pairs\") selecting Automatic (native Bruker order), \"Ascending\", or \"Descending\".";

ReadBrukerParameters::nofile = "Bruker parameter file `1` does not exist or is not a file.";
ReadBrukerParameters::parse = "Failed to parse Bruker parameter file `1`: `2`.";

ImportBrukerSpectrum::nodir = "Bruker dataset path `1` does not exist.";
ImportBrukerSpectrum::badprocno = "Invalid ProcessingNumber `1`; expected Automatic, a positive integer, or a directory-name string.";
ImportBrukerSpectrum::noprocs = "Processed Bruker parameter file 'procs' not found under `1`.";
ImportBrukerSpectrum::no1r = "Processed real spectrum file '1r' not found under `1`.";
ImportBrukerSpectrum::noimag = "Processed imaginary spectrum file '1i' not found under `1` (required by ReadImaginary -> True).";
ImportBrukerSpectrum::badmeta = "Missing or invalid Bruker processing parameter `1` in `2` (got `3`).";
ImportBrukerSpectrum::baddtype = "Unsupported Bruker DTYPP value `1` in `2` (supported: 0 for 32-bit integer, 2 for 64-bit IEEE float).";
ImportBrukerSpectrum::badbyteorder = "Unsupported Bruker BYTORDP value `1` in `2` (supported: 0 for little-endian, 1 for big-endian).";
ImportBrukerSpectrum::truncated = "Bruker binary file `1` is truncated or smaller than expected (`2` bytes found, `3` bytes expected).";
ImportBrukerSpectrum::sizemismatch = "Bruker binary file `1` size (`2` bytes) does not match expected size (`3` bytes).";
ImportBrukerSpectrum::badscaling = "Invalid IntensityScaling `1`; expected \"Bruker\", \"Raw\", \"Normalized\", or a positive number.";
ImportBrukerSpectrum::zeronorm = "Spectrum has zero amplitude; cannot normalize intensities. Returning unnormalized Bruker-scaled data.";
ImportBrukerSpectrum::badret = "Invalid ReturnFormat `1`; expected \"BrukerSpectrum\", \"Association\", or \"Pairs\".";

BrukerSpectrum::noimag = "This BrukerSpectrum does not contain an imaginary component (HasImaginary -> False).";
BrukerSpectrum::badprop = "Unknown BrukerSpectrum property `1`. Use spec[\"Properties\"] to list available keys.";
BrukerSpectrum::incomp = "BrukerSpectrum datasets have incompatible axes or lengths and cannot be combined.";

PlotBrukerSpectrum::badunits = "Invalid AxisUnits `1`; expected \"ppm\", \"Hz\", or \"rad/s\".";
PlotBrukerSpectrum::badcomp = "Invalid Component `1`; expected \"Real\", \"Imaginary\", \"Both\", \"Complex\", or \"Magnitude\".";
ToSignal::noimag = "Converting a BrukerSpectrum to a time-domain Signal uses complex data if available; imaginary part is absent, assuming zero imaginary spectrum.";

Begin["`Private`"];

(* --- Low-level JCAMP-DX parameter parser --- *)

stripJcampComment[line_String] := Module[
    {chars = Characters[line], inAngle = False, n, i, ch},
    n = Length[chars];
    If[!StringContainsQ[line, "$$"], Return[line]];
    For[i = 1, i <= n, i++,
        ch = chars[[i]];
        If[ch === "<", inAngle = True];
        If[ch === ">", inAngle = False];
        If[!inAngle && ch === "$" && i < n && chars[[i + 1]] === "$",
            Return[StringJoin[Take[chars, i - 1]]]
        ]
    ];
    line
];

parseJcampScalar[raw_String] := Module[{s = StringTrim[raw], norm},
    If[s === "", Return[""]];
    If[StringStartsQ[s, "<"] && StringEndsQ[s, ">"],
        Return[StringTrim[StringTake[s, {2, -2}]]]
    ];
    If[s === "yes", Return[True]];
    If[s === "no", Return[False]];
    If[StringMatchQ[s, RegularExpression["[+-]?[0-9]+"]],
        Return[ToExpression[s]]
    ];
    If[StringMatchQ[s, RegularExpression["[+-]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eEdD][+-]?[0-9]+)?"]],
        norm = StringReplace[s, {"e" -> "*^", "E" -> "*^", "d" -> "*^", "D" -> "*^"}];
        Return[N[ToExpression[norm]]]
    ];
    s
];

parseJcampArray[raw_String] := Module[{tokens, parsed},
    tokens = StringCases[
        StringTrim[raw],
        ("<" ~~ Shortest[___] ~~ ">") | (Except[WhitespaceCharacter] ..)
    ];
    parsed = parseJcampScalar /@ tokens;
    If[parsed === {}, Return[{}]];
    If[AllTrue[parsed, IntegerQ], Return[Developer`ToPackedArray[parsed]]];
    If[AllTrue[parsed, NumericQ], Return[Developer`ToPackedArray[N[parsed]]]];
    parsed
];

ReadBrukerParameters[filePath_String] := Module[
    {fullPath, text, lines, n, i, line, clean, body, pos, key, valPart,
     isArr, acc, params = <||>, nextTrim},
    If[!FileExistsQ[filePath] || FileType[filePath] =!= File,
        Message[ReadBrukerParameters::nofile, filePath];
        Return[$Failed]
    ];
    fullPath = ExpandFileName[filePath];
    text = Quiet[ReadString[fullPath]];
    If[!StringQ[text],
        Message[ReadBrukerParameters::parse, fullPath, "unreadable file"];
        Return[$Failed]
    ];
    lines = StringSplit[text, "\r\n" | "\n" | "\r"];
    n = Length[lines];
    i = 1;
    While[i <= n,
        line = lines[[i]];
        If[StringStartsQ[StringTrim[line], "##END="], Break[]];
        If[StringStartsQ[line, "##$"],
            clean = stripJcampComment[line];
            pos = StringPosition[clean, "=", 1];
            If[pos =!= {},
                key = StringTrim[StringTake[clean, {4, pos[[1, 1]] - 1}]];
                valPart = StringTrim[StringDrop[clean, pos[[1, 2]]]];
                isArr = StringMatchQ[valPart, RegularExpression["\\(\\s*[0-9]+\\s*\\.\\.\\s*[0-9]+\\s*\\).*"]];
                If[isArr,
                    body = StringReplace[valPart, RegularExpression["^\\(\\s*[0-9]+\\s*\\.\\.\\s*[0-9]+\\s*\\)\\s*"] -> ""];
                    acc = {body};
                    While[i < n,
                        nextTrim = StringTrim[lines[[i + 1]]];
                        If[StringStartsQ[nextTrim, "##"] || StringStartsQ[nextTrim, "$$"], Break[]];
                        i++;
                        AppendTo[acc, stripJcampComment[lines[[i]]]]
                    ];
                    params[key] = parseJcampArray[StringRiffle[acc, " "]],
                    acc = {valPart};
                    If[StringStartsQ[valPart, "<"] && !StringEndsQ[valPart, ">"],
                        While[i < n,
                            nextTrim = StringTrim[lines[[i + 1]]];
                            If[StringStartsQ[nextTrim, "##"], Break[]];
                            i++;
                            AppendTo[acc, lines[[i]]];
                            If[StringContainsQ[lines[[i]], ">"], Break[]]
                        ]
                    ];
                    params[key] = parseJcampScalar[StringRiffle[acc, "\n"]]
                ]
            ]
        ];
        i++
    ];
    params
];

(* --- Low-level directory/file resolver & binary reader --- *)

selectProcSubdir[pdataDir_String, procNoOpt_] := Module[
    {strOpt, subdirs, validDirs, intMap},
    If[procNoOpt === Automatic,
        If[DirectoryQ[FileNameJoin[{pdataDir, "1"}]] && FileExistsQ[FileNameJoin[{pdataDir, "1", "procs"}]],
            Return["1"]
        ];
        subdirs = Select[FileNames["*", pdataDir], DirectoryQ];
        validDirs = Select[subdirs, FileExistsQ[FileNameJoin[{#, "procs"}]] &];
        If[validDirs === {}, Return[$Failed]];
        intMap = SortBy[
            Select[
                Map[With[{nm = FileNameTake[#]}, {If[StringMatchQ[nm, RegularExpression["[0-9]+"]], FromDigits[nm], Infinity], nm}] &, validDirs],
                NumericQ[#[[1]]] &
            ],
            First
        ];
        If[intMap =!= {}, Return[intMap[[1, 2]]]];
        Return[FileNameTake[First[Sort[validDirs]]]],
        strOpt = Which[
            IntegerQ[procNoOpt] && procNoOpt >= 1, IntegerString[procNoOpt],
            StringQ[procNoOpt] && StringTrim[procNoOpt] =!= "", StringTrim[procNoOpt],
            True, Return[$Failed]
        ];
        If[DirectoryQ[FileNameJoin[{pdataDir, strOpt}]],
            strOpt,
            $Failed
        ]
    ]
];

resolveBrukerPaths[pathIn_String, procNoOpt_] := Module[
    {path, fName, parentDir, pdataDir, procDir, expDir, procNo, procsFile, rFile, iFile, acqusFile, titleFile},
    If[!FileExistsQ[pathIn] && !DirectoryQ[pathIn],
        Message[ImportBrukerSpectrum::nodir, pathIn];
        Return[$Failed]
    ];
    path = ExpandFileName[pathIn];
    Which[
        (* Direct path to 1r or procs *)
        FileType[path] === File,
            fName = FileNameTake[path];
            procDir = DirectoryName[path];
            parentDir = DirectoryName[procDir];
            If[procNoOpt =!= Automatic,
                If[ToLowerCase[FileNameTake[parentDir]] === "pdata",
                    procNo = selectProcSubdir[parentDir, procNoOpt];
                    If[procNo === $Failed,
                        Message[ImportBrukerSpectrum::badprocno, procNoOpt]; Return[$Failed]
                    ];
                    procDir = FileNameJoin[{parentDir, procNo}],
                    Message[ImportBrukerSpectrum::badprocno, procNoOpt]; Return[$Failed]
                ]
            ];
            expDir = If[ToLowerCase[FileNameTake[DirectoryName[procDir]]] === "pdata",
                DirectoryName[DirectoryName[procDir]],
                None
            ];
            procNo = FileNameTake[procDir],

        (* Directly a pdata/<n> directory *)
        FileExistsQ[FileNameJoin[{path, "procs"}]] || FileExistsQ[FileNameJoin[{path, "1r"}]],
            parentDir = DirectoryName[path];
            If[procNoOpt =!= Automatic && ToLowerCase[FileNameTake[parentDir]] === "pdata",
                procNo = selectProcSubdir[parentDir, procNoOpt];
                If[procNo === $Failed,
                    Message[ImportBrukerSpectrum::badprocno, procNoOpt]; Return[$Failed]
                ];
                procDir = FileNameJoin[{parentDir, procNo}],
                procDir = path;
                procNo = FileNameTake[procDir]
            ];
            expDir = If[ToLowerCase[FileNameTake[DirectoryName[procDir]]] === "pdata",
                DirectoryName[DirectoryName[procDir]],
                None
            ],

        (* A pdata directory *)
        ToLowerCase[FileNameTake[path]] === "pdata" || DirectoryQ[FileNameJoin[{path, "1"}]],
            pdataDir = path;
            procNo = selectProcSubdir[pdataDir, procNoOpt];
            If[procNo === $Failed,
                If[procNoOpt === Automatic,
                    Message[ImportBrukerSpectrum::noprocs, pdataDir],
                    Message[ImportBrukerSpectrum::badprocno, procNoOpt]
                ];
                Return[$Failed]
            ];
            procDir = FileNameJoin[{pdataDir, procNo}];
            expDir = DirectoryName[pdataDir],

        (* Experiment directory containing pdata/ *)
        DirectoryQ[FileNameJoin[{path, "pdata"}]],
            expDir = path;
            pdataDir = FileNameJoin[{path, "pdata"}];
            procNo = selectProcSubdir[pdataDir, procNoOpt];
            If[procNo === $Failed,
                If[procNoOpt === Automatic,
                    Message[ImportBrukerSpectrum::noprocs, pdataDir],
                    Message[ImportBrukerSpectrum::badprocno, procNoOpt]
                ];
                Return[$Failed]
            ];
            procDir = FileNameJoin[{pdataDir, procNo}],

        True,
            Message[ImportBrukerSpectrum::noprocs, path];
            Return[$Failed]
    ];

    procsFile = FileNameJoin[{procDir, "procs"}];
    rFile = FileNameJoin[{procDir, "1r"}];
    iFile = FileNameJoin[{procDir, "1i"}];
    If[!FileExistsQ[procsFile],
        Message[ImportBrukerSpectrum::noprocs, procDir];
        Return[$Failed]
    ];
    If[!FileExistsQ[rFile],
        Message[ImportBrukerSpectrum::no1r, procDir];
        Return[$Failed]
    ];
    acqusFile = If[StringQ[expDir] && FileExistsQ[FileNameJoin[{expDir, "acqus"}]],
        FileNameJoin[{expDir, "acqus"}],
        None
    ];
    titleFile = If[FileExistsQ[FileNameJoin[{procDir, "title"}]],
        FileNameJoin[{procDir, "title"}],
        None
    ];

    <|
        "SourceDirectory" -> If[StringQ[expDir], expDir, procDir],
        "ExperimentDirectory" -> expDir,
        "ProcessedDirectory" -> procDir,
        "ProcessingNumber" -> If[StringMatchQ[procNo, RegularExpression["[0-9]+"]], FromDigits[procNo], procNo],
        "ProcsFile" -> procsFile,
        "RealFile" -> rFile,
        "ImaginaryFile" -> If[FileExistsQ[iFile], iFile, None],
        "AcqusFile" -> acqusFile,
        "TitleFile" -> titleFile
    |>
];

readBrukerBinaryVector[filePath_String, dtypp_, bytordp_, expectedPoints_Integer, procsFile_String] := Module[
    {specPair, fmt, bytesPerPoint, bo, actualBytes, expectedBytes, data},
    specPair = Switch[dtypp,
        0, {"Integer32", 4},
        2, {"Real64", 8},
        _, Message[ImportBrukerSpectrum::baddtype, dtypp, procsFile]; $Failed
    ];
    If[specPair === $Failed, Return[$Failed]];
    {fmt, bytesPerPoint} = specPair;
    bo = Switch[bytordp,
        0, -1,
        1, 1,
        _, Message[ImportBrukerSpectrum::badbyteorder, bytordp, procsFile]; Return[$Failed]
    ];
    actualBytes = FileByteCount[filePath];
    expectedBytes = expectedPoints * bytesPerPoint;
    If[!IntegerQ[actualBytes] || actualBytes < expectedBytes,
        Message[ImportBrukerSpectrum::truncated, filePath, actualBytes, expectedBytes];
        Return[$Failed]
    ];
    If[actualBytes =!= expectedBytes,
        Message[ImportBrukerSpectrum::sizemismatch, filePath, actualBytes, expectedBytes];
        Return[$Failed]
    ];
    data = Quiet[BinaryReadList[filePath, fmt, ByteOrdering -> bo]];
    If[!ListQ[data] || Length[data] =!= expectedPoints || !VectorQ[data, NumericQ],
        Message[ImportBrukerSpectrum::truncated, filePath, actualBytes, expectedBytes];
        Return[$Failed]
    ];
    Developer`ToPackedArray[data]
];

(* --- High-level spectrum axis & scaling construction --- *)

buildChemicalShiftAxis[procs_Association, nPts_Integer, procsFile_String] := Module[
    {sf, swP, offset, swPPM, reverseQ, f1p, f2p, leftPPM, rightPPM, ppmAxis},
    sf = Lookup[procs, "SF", Missing["KeyAbsent"]];
    swP = Lookup[procs, "SW_p", Missing["KeyAbsent"]];
    offset = Lookup[procs, "OFFSET", Missing["KeyAbsent"]];
    If[!NumericQ[sf] || !TrueQ[sf > 0.],
        Message[ImportBrukerSpectrum::badmeta, "SF", procsFile, sf]; Return[$Failed]
    ];
    If[!NumericQ[swP] || !TrueQ[swP > 0.],
        Message[ImportBrukerSpectrum::badmeta, "SW_p", procsFile, swP]; Return[$Failed]
    ];
    If[!NumericQ[offset],
        Message[ImportBrukerSpectrum::badmeta, "OFFSET", procsFile, offset]; Return[$Failed]
    ];
    swPPM = N[swP] / N[sf];
    reverseQ = TrueQ[Lookup[procs, "REVERSE", False]];
    f1p = Lookup[procs, "F1P", None];
    f2p = Lookup[procs, "F2P", None];
    If[nPts === 1,
        Return[<|
            "PPMAxis" -> Developer`ToPackedArray[{N[offset]}],
            "SF" -> N[sf],
            "SWHz" -> N[swP],
            "SWPPM" -> swPPM,
            "OffsetPPM" -> N[offset],
            "ReversedQ" -> reverseQ
        |>]
    ];
    (* In Bruker procs, OFFSET is the chemical shift of the first point (left edge)
       and the spectrum steps downward by swPPM across nPts points unless REVERSE=yes
       flipped the stored data with ascending shift. *)
    If[reverseQ,
        If[NumericQ[f1p] && NumericQ[f2p] && TrueQ[f1p < f2p] && Abs[(f2p - f1p) - swPPM] < 1.*^-3 * Max[1., swPPM],
            leftPPM = N[offset]; rightPPM = N[offset + swPPM],
            leftPPM = N[offset - swPPM]; rightPPM = N[offset]
        ],
        leftPPM = N[offset]; rightPPM = N[offset - swPPM]
    ];
    ppmAxis = Developer`ToPackedArray[Subdivide[leftPPM, rightPPM, nPts - 1]];
    <|
        "PPMAxis" -> ppmAxis,
        "SF" -> N[sf],
        "SWHz" -> N[swP],
        "SWPPM" -> swPPM,
        "OffsetPPM" -> N[offset],
        "ReversedQ" -> reverseQ
    |>
];

applyIntensityScaling[rawReal_List, rawImag_, hasImag_?BooleanQ, ncProc_, scalingOpt_, includeRaw_?BooleanQ] := Module[
    {brukerScale, scMode, extraMult = 1., appliedScale, scaledReal, scaledImag, normFactor = 1., normalizedQ = False},
    brukerScale = N[2.^ncProc];
    Which[
        scalingOpt === "Raw",
            scMode = "Raw";
            appliedScale = 1.;
            scaledReal = rawReal;
            scaledImag = rawImag,
        scalingOpt === "Bruker" || scalingOpt === Automatic,
            scMode = "Bruker";
            appliedScale = brukerScale;
            scaledReal = If[TrueQ[brukerScale === 1.] && Developer`PackedArrayQ[rawReal] && Head[First[rawReal]] === Real,
                rawReal,
                Developer`ToPackedArray[N[rawReal] * brukerScale]
            ];
            scaledImag = If[hasImag,
                If[TrueQ[brukerScale === 1.] && Developer`PackedArrayQ[rawImag] && Head[First[rawImag]] === Real,
                    rawImag,
                    Developer`ToPackedArray[N[rawImag] * brukerScale]
                ],
                None
            ],
        scalingOpt === "Normalized",
            scMode = "Normalized";
            scaledReal = Developer`ToPackedArray[N[rawReal] * brukerScale];
            scaledImag = If[hasImag, Developer`ToPackedArray[N[rawImag] * brukerScale], None];
            normFactor = If[hasImag,
                Max[Sqrt[scaledReal^2 + scaledImag^2]],
                Max[Abs[scaledReal]]
            ];
            If[!TrueQ[normFactor > 0.],
                Message[ImportBrukerSpectrum::zeronorm];
                appliedScale = brukerScale;
                normalizedQ = False;
                normFactor = 1.,
                scaledReal = Developer`ToPackedArray[scaledReal / normFactor];
                If[hasImag, scaledImag = Developer`ToPackedArray[scaledImag / normFactor]];
                appliedScale = brukerScale / normFactor;
                normalizedQ = True
            ],
        NumericQ[scalingOpt] && TrueQ[scalingOpt > 0],
            scMode = "Scaled";
            extraMult = N[scalingOpt];
            appliedScale = brukerScale * extraMult;
            scaledReal = Developer`ToPackedArray[N[rawReal] * appliedScale];
            scaledImag = If[hasImag, Developer`ToPackedArray[N[rawImag] * appliedScale], None],
        True,
            Message[ImportBrukerSpectrum::badscaling, scalingOpt];
            Return[$Failed]
    ];
    <|
        "Real" -> scaledReal,
        "Imaginary" -> scaledImag,
        "RawReal" -> If[TrueQ[includeRaw], rawReal, None],
        "RawImaginary" -> If[TrueQ[includeRaw] && hasImag, rawImag, None],
        "ScalingMode" -> scMode,
        "NCProc" -> ncProc,
        "BrukerScaleFactor" -> brukerScale,
        "AppliedScaleFactor" -> appliedScale,
        "NormalizedQ" -> normalizedQ,
        "NormalizationFactor" -> normFactor
    |>
];

Options[ImportBrukerSpectrum] = {
    ProcessingNumber -> Automatic,
    ReadImaginary -> Automatic,
    IntensityScaling -> "Bruker",
    IncludeRawData -> True,
    ReadAcquisitionParameters -> True,
    ReturnFormat -> "BrukerSpectrum",
    AxisUnits -> "ppm",
    Component -> "Real",
    AxisOrder -> Automatic
};

ImportBrukerSpectrum[path_String, opts___Rule] := Module[
    {procNoOpt, readImagOpt, scalingOpt, includeRawOpt, readAcqOpt, retFmt, axisUnitsOpt, compOpt, axisOrderOpt,
     paths, procs, acqus = <||>, titleStr = "", si, stsi, stsr, nPts, dtypp, bytordp, ncProc,
     axisInfo, ppmAxis, freqHz, omegaRad, rawReal, rawImag = None, hasImag = False,
     scaled, nuc, specAssoc, spec},
    {procNoOpt, readImagOpt, scalingOpt, includeRawOpt, readAcqOpt, retFmt, axisUnitsOpt, compOpt, axisOrderOpt} =
        {ProcessingNumber, ReadImaginary, IntensityScaling, IncludeRawData, ReadAcquisitionParameters,
         ReturnFormat, AxisUnits, Component, AxisOrder} /. {opts} /. Options[ImportBrukerSpectrum];

    If[!MemberQ[{"BrukerSpectrum", "Association", "Pairs"}, retFmt],
        Message[ImportBrukerSpectrum::badret, retFmt];
        Return[$Failed]
    ];

    paths = resolveBrukerPaths[path, procNoOpt];
    If[paths === $Failed, Return[$Failed]];

    procs = ReadBrukerParameters[paths["ProcsFile"]];
    If[procs === $Failed, Return[$Failed]];

    If[TrueQ[readAcqOpt] && StringQ[paths["AcqusFile"]],
        acqus = Quiet[ReadBrukerParameters[paths["AcqusFile"]]];
        If[!AssociationQ[acqus], acqus = <||>]
    ];
    If[StringQ[paths["TitleFile"]],
        titleStr = Quiet[StringTrim[ReadString[paths["TitleFile"]]]];
        If[!StringQ[titleStr], titleStr = ""]
    ];

    si = Lookup[procs, "SI", Missing["KeyAbsent"]];
    If[!IntegerQ[si] || si < 1,
        Message[ImportBrukerSpectrum::badmeta, "SI", paths["ProcsFile"], si]; Return[$Failed]
    ];
    stsi = Lookup[procs, "STSI", 0];
    stsr = Lookup[procs, "STSR", 0];
    If[!IntegerQ[stsi] || stsi < 0,
        Message[ImportBrukerSpectrum::badmeta, "STSI", paths["ProcsFile"], stsi]; Return[$Failed]
    ];
    nPts = If[stsi > 0, stsi, si];

    dtypp = Lookup[procs, "DTYPP", 0];
    bytordp = Lookup[procs, "BYTORDP", 0];
    ncProc = Lookup[procs, "NC_proc", 0];
    If[!NumericQ[ncProc],
        Message[ImportBrukerSpectrum::badmeta, "NC_proc", paths["ProcsFile"], ncProc]; Return[$Failed]
    ];

    axisInfo = buildChemicalShiftAxis[procs, nPts, paths["ProcsFile"]];
    If[axisInfo === $Failed, Return[$Failed]];
    ppmAxis = axisInfo["PPMAxis"];
    freqHz = Developer`ToPackedArray[ppmAxis * axisInfo["SF"]];
    omegaRad = Developer`ToPackedArray[(2. * Pi) * freqHz];

    rawReal = readBrukerBinaryVector[paths["RealFile"], dtypp, bytordp, nPts, paths["ProcsFile"]];
    If[rawReal === $Failed, Return[$Failed]];

    Which[
        readImagOpt === True,
            If[!StringQ[paths["ImaginaryFile"]],
                Message[ImportBrukerSpectrum::noimag, paths["ProcessedDirectory"]]; Return[$Failed]
            ];
            rawImag = readBrukerBinaryVector[paths["ImaginaryFile"], dtypp, bytordp, nPts, paths["ProcsFile"]];
            If[rawImag === $Failed, Return[$Failed]];
            hasImag = True,
        readImagOpt === Automatic,
            If[StringQ[paths["ImaginaryFile"]],
                rawImag = readBrukerBinaryVector[paths["ImaginaryFile"], dtypp, bytordp, nPts, paths["ProcsFile"]];
                If[rawImag === $Failed, Return[$Failed]];
                hasImag = True,
                hasImag = False; rawImag = None
            ],
        True,
            hasImag = False; rawImag = None
    ];

    scaled = applyIntensityScaling[rawReal, rawImag, hasImag, ncProc, scalingOpt, TrueQ[includeRawOpt]];
    If[scaled === $Failed, Return[$Failed]];

    nuc = Lookup[procs, "AXNUC", Lookup[acqus, "NUC1", ""]];
    If[nuc === "" || nuc === "off", nuc = Lookup[acqus, "NUC1", ""]];

    specAssoc = <|
        "Real" -> scaled["Real"],
        "Imaginary" -> scaled["Imaginary"],
        "RawReal" -> scaled["RawReal"],
        "RawImaginary" -> scaled["RawImaginary"],
        "HasImaginary" -> hasImag,
        "ChemicalShift" -> ppmAxis,
        "Frequency" -> freqHz,
        "AngularFrequency" -> omegaRad,
        "NumPoints" -> nPts,
        "SI" -> si,
        "STSI" -> stsi,
        "STSR" -> stsr,
        "StripTransformQ" -> (stsi > 0 && stsi < si),
        "ReversedQ" -> axisInfo["ReversedQ"],
        "SpectralWidthHz" -> axisInfo["SWHz"],
        "SpectralWidthPPM" -> axisInfo["SWPPM"],
        "ReferenceFrequencyMHz" -> axisInfo["SF"],
        "OffsetPPM" -> axisInfo["OffsetPPM"],
        "NCProc" -> scaled["NCProc"],
        "BrukerScaleFactor" -> scaled["BrukerScaleFactor"],
        "AppliedScaleFactor" -> scaled["AppliedScaleFactor"],
        "IntensityScaling" -> scaled["ScalingMode"],
        "NormalizedQ" -> scaled["NormalizedQ"],
        "NormalizationFactor" -> scaled["NormalizationFactor"],
        "ByteOrder" -> If[bytordp === 0, "LittleEndian", "BigEndian"],
        "BYTORDP" -> bytordp,
        "DataType" -> If[dtypp === 0, "Integer32", "Real64"],
        "DTYPP" -> dtypp,
        "Nucleus" -> nuc,
        "Title" -> titleStr,
        "SourceDirectory" -> paths["SourceDirectory"],
        "ExperimentDirectory" -> paths["ExperimentDirectory"],
        "ProcessedDirectory" -> paths["ProcessedDirectory"],
        "ProcessingNumber" -> paths["ProcessingNumber"],
        "ProcessingParameters" -> procs,
        "AcquisitionParameters" -> acqus
    |>;

    If[retFmt === "Association", Return[specAssoc]];
    spec = BrukerSpectrum[specAssoc];
    If[retFmt === "Pairs",
        Return[Expand[spec, AxisUnits -> axisUnitsOpt, Component -> compOpt, AxisOrder -> axisOrderOpt]]
    ];
    spec
];

(* --- BrukerSpectrum object methods, properties, formatting, and SpinDynamica UpValues --- *)

getAxisVector[a_Association, units_String] := Switch[units,
    "ppm", a["ChemicalShift"],
    "Hz", a["Frequency"],
    "rad/s", a["AngularFrequency"],
    _, Message[PlotBrukerSpectrum::badunits, units]; $Failed
];

getComponentVector[a_Association, comp_String] := Switch[comp,
    "Real", a["Real"],
    "Imaginary",
        If[!TrueQ[a["HasImaginary"]],
            Message[BrukerSpectrum::noimag]; $Failed,
            a["Imaginary"]
        ],
    "Complex",
        If[TrueQ[a["HasImaginary"]],
            Developer`ToPackedArray[a["Real"] + I * a["Imaginary"]],
            Developer`ToPackedArray[a["Real"] + 0. I]
        ],
    "Magnitude",
        If[TrueQ[a["HasImaginary"]],
            Developer`ToPackedArray[Sqrt[a["Real"]^2 + a["Imaginary"]^2]],
            Developer`ToPackedArray[Abs[a["Real"]]]
        ],
    _, Message[PlotBrukerSpectrum::badcomp, comp]; $Failed
];

buildPairsFromSpec[a_Association, units_String, comp_String, orderOpt_] := Module[
    {x, y, pairs, ord},
    x = getAxisVector[a, units];
    If[x === $Failed, Return[$Failed]];
    y = getComponentVector[a, comp];
    If[y === $Failed, Return[$Failed]];
    ord = Switch[orderOpt,
        Automatic, "Native",
        "Ascending", If[Length[x] > 1 && First[x] > Last[x], "Reverse", "Native"],
        "Descending", If[Length[x] > 1 && First[x] < Last[x], "Reverse", "Native"],
        _, "Native"
    ];
    pairs = If[ord === "Reverse",
        Transpose[{Reverse[x], Reverse[y]}],
        Transpose[{x, y}]
    ];
    Developer`ToPackedArray[pairs]
];

BrukerSpectrum[a_Association]["Properties"] := Sort @ Join[Keys[a], {"Complex", "PPM", "Hz", "Points", "Metadata"}];
BrukerSpectrum[a_Association]["Complex"] := getComponentVector[a, "Complex"];
BrukerSpectrum[a_Association]["PPM"] := a["ChemicalShift"];
BrukerSpectrum[a_Association]["Hz"] := a["Frequency"];
BrukerSpectrum[a_Association]["Points"] := a["NumPoints"];
BrukerSpectrum[a_Association]["Metadata"] := KeyDrop[a, {"Real", "Imaginary", "RawReal", "RawImaginary", "ChemicalShift", "Frequency", "AngularFrequency"}];
BrukerSpectrum[a_Association][key_String] := If[KeyExistsQ[a, key], a[key], Message[BrukerSpectrum::badprop, key]; $Failed];
BrukerSpectrum[a_Association]["Pairs", units_String: "ppm", comp_String: "Real", order_: Automatic] :=
    buildPairsFromSpec[a, units, comp, order];

BrukerSpectrum /: Length[BrukerSpectrum[a_Association]] := a["NumPoints"];

BrukerSpectrum /: Expand[BrukerSpectrum[a_Association], opts___Rule] := Module[{u, c, ord},
    {u, c, ord} = {AxisUnits, Component, AxisOrder} /. {opts} /. {AxisUnits -> "ppm", Component -> "Real", AxisOrder -> Automatic};
    buildPairsFromSpec[a, u, c, ord]
];

BrukerSpectrum /: Re[BrukerSpectrum[a_Association]] :=
    Expand[BrukerSpectrum[a], Component -> "Real"];

BrukerSpectrum /: Im[BrukerSpectrum[a_Association]] :=
    Expand[BrukerSpectrum[a], Component -> "Imaginary"];

BrukerSpectrum /: Abs[BrukerSpectrum[a_Association]] :=
    Expand[BrukerSpectrum[a], Component -> "Magnitude"];

BrukerSpectrum /: HoldPattern[Times[c__?NumericQ, BrukerSpectrum[a_Association]]] := Module[
    {cN = N[Times[c]], z, newReal, newImag, hasImagOut},
    If[TrueQ[Im[cN] == 0],
        newReal = Developer`ToPackedArray[Re[cN] * N[a["Real"]]];
        newImag = If[TrueQ[a["HasImaginary"]], Developer`ToPackedArray[Re[cN] * N[a["Imaginary"]]], None];
        hasImagOut = TrueQ[a["HasImaginary"]],
        z = N[Re[a["Real"]]] + I * If[TrueQ[a["HasImaginary"]], N[a["Imaginary"]], 0.];
        z = cN * z;
        newReal = Developer`ToPackedArray[Re[z]];
        newImag = Developer`ToPackedArray[Im[z]];
        hasImagOut = True
    ];
    BrukerSpectrum[Join[a, <|
        "Real" -> newReal,
        "Imaginary" -> newImag,
        "HasImaginary" -> hasImagOut,
        "AppliedScaleFactor" -> a["AppliedScaleFactor"] * Abs[cN],
        "IntensityScaling" -> "Scaled"
    |>]]
];

CombineData[specs:{_BrukerSpectrum..}, opts___Rule] :=
    CombineData[specs, ConstantArray[1. / Length[specs], Length[specs]], opts];

CombineData[specs:{_BrukerSpectrum..}, weights_List, opts___Rule] := Module[
    {n = Length[specs], a0, pts0, ppm0, wNum, reals, imags, hasAnyImag, combinedComplex, hasImagOut, newReal, newImag},
    If[Length[weights] =!= n,
        Message[CombineData::incompw]; Return[$Failed]
    ];
    If[!AllTrue[weights, NumericQ], Return[$Failed]];
    a0 = specs[[1, 1]];
    pts0 = a0["NumPoints"];
    ppm0 = a0["ChemicalShift"];
    If[!AllTrue[Rest[specs], (#[[1]]["NumPoints"] === pts0 && Max[Abs[#[[1]]["ChemicalShift"] - ppm0]] <= 1.*^-10 * Max[1., Abs[a0["SpectralWidthPPM"]]]) &],
        Message[BrukerSpectrum::incomp]; Return[$Failed]
    ];
    wNum = N[weights];
    hasAnyImag = AnyTrue[specs, TrueQ[#[[1]]["HasImaginary"]] &] || AnyTrue[wNum, !TrueQ[Im[#] == 0] &];
    If[hasAnyImag,
        combinedComplex = Sum[
            wNum[[k]] * (N[specs[[k, 1]]["Real"]] + I * If[TrueQ[specs[[k, 1]]["HasImaginary"]], N[specs[[k, 1]]["Imaginary"]], 0.]),
            {k, n}
        ];
        newReal = Developer`ToPackedArray[Re[combinedComplex]];
        newImag = Developer`ToPackedArray[Im[combinedComplex]];
        hasImagOut = True,
        reals = Sum[wNum[[k]] * N[specs[[k, 1]]["Real"]], {k, n}];
        newReal = Developer`ToPackedArray[reals];
        newImag = None;
        hasImagOut = False
    ];
    BrukerSpectrum[Join[a0, <|
        "Real" -> newReal,
        "Imaginary" -> newImag,
        "RawReal" -> None,
        "RawImaginary" -> None,
        "HasImaginary" -> hasImagOut,
        "IntensityScaling" -> "Combined"
    |>]]
];

(* Ensure BrukerSpectrum CombineData rules precede Signal.m's Except[Signal[___]] catch-all *)
DownValues[CombineData] = Join[
    Select[DownValues[CombineData], !FreeQ[#, BrukerSpectrum] &],
    Select[DownValues[CombineData], FreeQ[#, BrukerSpectrum] &]
];

BrukerSpectrum /: Format[BrukerSpectrum[a_Association]] := RowBox[{
    "BrukerSpectrum[",
    ToBoxes[If[a["Nucleus"] =!= "", a["Nucleus"] <> ", ", ""]],
    ToBoxes[a["NumPoints"]], " pts, SF=", ToBoxes[NumberForm[a["ReferenceFrequencyMHz"], {6, 2}]], " MHz, ",
    ToBoxes[NumberForm[First[a["ChemicalShift"]], 4]], "..", ToBoxes[NumberForm[Last[a["ChemicalShift"]], 4]], " ppm",
    If[TrueQ[a["HasImaginary"]], ", complex", ", real"],
    "]"
}] // DisplayForm;

BrukerSpectrum /: Format[BrukerSpectrum[a_Association], OutputForm] := StringForm[
    "BrukerSpectrum[`1``2` pts, SF=`3` MHz, `4`..`5` ppm, `6`]",
    If[a["Nucleus"] =!= "", a["Nucleus"] <> ", ", ""],
    a["NumPoints"],
    ToString[NumberForm[a["ReferenceFrequencyMHz"], {6, 2}]],
    ToString[NumberForm[First[a["ChemicalShift"]], 4]],
    ToString[NumberForm[Last[a["ChemicalShift"]], 4]],
    If[TrueQ[a["HasImaginary"]], "complex", "real"]
];

(* --- Plotting & Signal Conversion --- *)

Options[PlotBrukerSpectrum] = Join[
    {
        AxisUnits -> "ppm",
        Component -> "Real",
        ReverseAxis -> Automatic
    },
    Options[ListPlot]
];

PlotBrukerSpectrum[BrukerSpectrum[a_Association], opts___Rule] := Module[
    {units, comp, revOpt, revQ, xVec, y1, y2, data, defaultXLabel, defaultYLabel, lpOpts, scOpt},
    {units, comp, revOpt} = {AxisUnits, Component, ReverseAxis} /. {opts} /. Options[PlotBrukerSpectrum];
    xVec = getAxisVector[a, units];
    If[xVec === $Failed, Return[$Failed]];
    If[!MemberQ[{"Real", "Imaginary", "Both", "Complex", "Magnitude"}, comp],
        Message[PlotBrukerSpectrum::badcomp, comp]; Return[$Failed]
    ];
    revQ = Switch[revOpt,
        Automatic, units === "ppm",
        True | False, revOpt,
        _, units === "ppm"
    ];
    data = Switch[comp,
        "Real",
            y1 = getComponentVector[a, "Real"];
            Transpose[{xVec, y1}],
        "Imaginary",
            y1 = getComponentVector[a, "Imaginary"];
            If[y1 === $Failed, Return[$Failed]];
            Transpose[{xVec, y1}],
        "Magnitude",
            y1 = getComponentVector[a, "Magnitude"];
            Transpose[{xVec, y1}],
        "Both" | "Complex",
            y1 = getComponentVector[a, "Real"];
            y2 = getComponentVector[a, "Imaginary"];
            If[y2 === $Failed, Return[$Failed]];
            {Transpose[{xVec, y1}], Transpose[{xVec, y2}]}
    ];
    defaultXLabel = Switch[units,
        "ppm", "Chemical Shift (ppm)",
        "Hz", "Frequency (Hz)",
        "rad/s", "Angular Frequency (rad/s)"
    ];
    defaultYLabel = Switch[a["IntensityScaling"],
        "Normalized", "Normalized Intensity",
        "Raw", "Raw Stored Intensity",
        _, "Intensity (a.u.)"
    ];
    scOpt = If[revQ, ScalingFunctions -> {"Reverse", None}, Sequence @@ {}];
    lpOpts = FilterRules[
        Join[
            {opts},
            {
                FrameLabel -> {defaultXLabel, defaultYLabel},
                Joined -> True,
                PlotRange -> All,
                Frame -> True,
                Axes -> None
            }
        ],
        Options[ListPlot]
    ];
    ListPlot[data, scOpt, Sequence @@ lpOpts]
];

BrukerSpectrum /: ListPlot[spec_BrukerSpectrum, opts___Rule] := PlotBrukerSpectrum[spec, opts];

ToSignal[BrukerSpectrum[a_Association], opts___Rule] := Module[
    {n, swHz, dt, freqs, specC, ordered, unshifted, fid, tmax},
    n = a["NumPoints"];
    swHz = a["SpectralWidthHz"];
    freqs = a["Frequency"];
    If[!TrueQ[a["HasImaginary"]], Message[ToSignal::noimag]];
    specC = getComponentVector[a, "Complex"];
    (* Reorder into ascending frequency order matching Signals`Signal`FT *)
    ordered = If[n > 1 && First[freqs] > Last[freqs], Reverse[specC], specC];
    (* Invert RotateRight[..., Floor[n/2] - 1] and divide-by-n normalization in Signals`Signal`FT *)
    unshifted = RotateLeft[ordered, Floor[n / 2] - 1];
    fid = Developer`ToPackedArray[n * InverseFourier[unshifted, FourierParameters -> {1, -1}]];
    (* Match Signals`Signal`FT grid spacing sw/n: SampleInterval = 1/swHz *)
    dt = 1. / swHz;
    tmax = (n - 1) * dt;
    Signal[{0., tmax, dt}, fid, LineBroadening -> None, ZeroFill -> 1, opts]
];

End[];
EndPackage[];
