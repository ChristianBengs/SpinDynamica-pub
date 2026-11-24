(* ::Package:: *)

(* SpinDynamica Optimal Control — Problem Compilation, Scaled L-BFGS / GRAPE & Converters *)

BeginPackage["OptimalControl`GRAPE`", {
    "Global`",
    "SpinDynamica`",
    "General`Matrices`",
    "QM`Hilbert`",
    "QM`SpinOperators`",
    "QM`Liouville`",
    "SD`Generators`",
    "SD`Propagation`",
    "Events`Events`",
    "Events`PulseSequence`",
    "Simulations`Simulations`",
    "OptimalControl`Kernels`",
    "OptimalControl`Objectives`"
}];

OptimalControlProblem::usage = "OptimalControlProblem[drift, controls, target, timing, opts] represents an optimal-control task on a SpinDynamica spin system with piecewise-constant or parameterized control Hamiltonians.";
CompileOptimalControlProblem::usage = "CompileOptimalControlProblem[problem] precomputes packed numerical matrices, scaling factors, control basis matrices, and objective evaluators from a symbolic OptimalControlProblem.";
EvaluateOptimalControl::usage = "EvaluateOptimalControl[problemOrCompiled, controls] computes the fidelity, regularized objective, and exact analytical gradient for a control or parameter array.";
VerifyOptimalControlGradient::usage = "VerifyOptimalControlGradient[problemOrCompiled, controls, opts] compares analytical GRAPE gradients against central finite differences and returns absolute and relative errors.";
StepSize::usage = "StepSize is an option for VerifyOptimalControlGradient specifying the central finite-difference step size (default 10^-7).";
GRAPEOptimize::usage = "GRAPEOptimize[problemOrEnsemble, opts] or GRAPEOptimize[drift, controls, target, timing, opts] optimizes pulse controls using analytical gradients and scaled quasi-Newton optimization, returning an OptimalControlResult.";
OptimalControlOptimize::usage = "OptimalControlOptimize is an alias for GRAPEOptimize supporting single-system, cooperative, and robust ensemble problems.";
OptimalControlResult::usage = "OptimalControlResult[<|...|>] encapsulates the optimized waveform, convergence diagnostics, fidelity, and problem metadata returned by GRAPEOptimize.";
OptimalControlEventList::usage = "OptimalControlEventList[result, opts] exports an OptimalControlResult as a list of SpinDynamica events {{gen_j, dt_j}, ...} compatible with NPropagate and Signal1D.";
OptimalControlGenerator::usage = "OptimalControlGenerator[result, tSym, opts] exports an OptimalControlResult as an order-0 piecewise-constant TimeDependentGeneratorFunction.";
OptimalControlPulseSequence::usage = "OptimalControlPulseSequence[result, opts] or OptimalControlPulseSequence[result, spins, opts] exports an OptimalControlResult as a SpinDynamica PulseSequence (ShapedPulse for single-channel XY, or simultaneous Delay with TimeDependentGeneratorFunction for multi-channel controls).";
ChannelStructure::usage = "ChannelStructure is an option for OptimalControlPulseSequence, OptimalControlPulseTable, and ExportOptimalControlPulse specifying how control channels pair into spin quadratures: Automatic or {spinSpec -> {chX, chY}, ...}.";
OptimalControlPulseTable::usage = "OptimalControlPulseTable[result, opts] extracts a structured Association and table of per-slice timings, Cartesian amplitudes (Hz), polar amplitudes (Hz and %), and phases (degrees or radians) per RF channel.";
ExportOptimalControlPulse::usage = "ExportOptimalControlPulse[result, format, opts] or ExportOptimalControlPulse[path, result, format, opts] exports optimized waveforms to \"Bruker\" (TopSpin JCAMP-DX), \"SIMPSON\", or \"CSV\" format.";
PhaseUnits::usage = "PhaseUnits is an option for OptimalControlPulseTable selecting \"Degrees\" (default) or \"Radians\".";
AmplitudeScaling::usage = "AmplitudeScaling is an option for OptimalControlPulseTable selecting \"Percent\" (default) or \"Hz\".";
Resample::usage = "Resample is an option for ExportOptimalControlPulse specifying Automatic, None, or an integer target slice count onto which piecewise-constant controls are uniformly resampled.";
Channel::usage = "Channel is an option for ExportOptimalControlPulse selecting which RF channel (1, 2, ..., or All) to export.";
Title::usage = "Title is an option for ExportOptimalControlPulse setting the header title in exported waveform files.";
OptimalControlWaveform::usage = "OptimalControlWaveform[result] returns a reusable numerical Association containing time steps, slice midpoints, control amplitudes, operators, and unit metadata.";
VerifyOptimalControlSimulation::usage = "VerifyOptimalControlSimulation[result, opts] independently simulates the exported OptimalControlEventList using SpinDynamica's NPropagationOperator / NPropagate and compares fidelities and propagators.";

PropagateSpace::usage = "PropagateSpace is an option for OptimalControlProblem and GRAPEOptimize selecting \"Hilbert\", \"Liouville\", or Automatic.";
ControlParameterization::usage = "ControlParameterization is an option selecting \"PiecewiseConstant\", {\"Fourier\", n}, {\"FourierSine\", n}, {\"Chebyshev\", deg}, {\"BSpline\", n, [deg]}, or {\"BasisMatrix\", mat}.";
ControlWindow::usage = "ControlWindow is an option applying an envelope window (None, \"Hann\", \"Hamming\", \"Blackman\", or a function on [0, 1]) to control waveforms.";
ControlBounds::usage = "ControlBounds is an option specifying amplitude bounds {min, max}, per-channel bounds {{min1, max1}, ...}, or None.";
AmplitudeBound::usage = "AmplitudeBound is an option specifying symmetric control bounds [-bound, +bound] across all channels.";
AmplitudeNormBound::usage = "AmplitudeNormBound is an option specifying polar/disk RF amplitude limits (Sqrt[u_x^2 + u_y^2] <= wMax) per slice across channel groups, e.g., wMax, {{1, 2} -> w1, {3, 4} -> w2}, or {{{1, 2}, w1}, ...}, in ControlUnits.";
L1Penalty::usage = "L1Penalty is an option specifying a nonsmooth L1 sparsity regularization weight lam1 * Sum[Abs[u[[j, k]]]] on control amplitudes (default 0).";
TVPenalty::usage = "TVPenalty is an option specifying a nonsmooth total-variation regularization weight lamTV * Sum[Abs[u[[j+1, k]] - u[[j, k]]]] on adjacent slice differences (default 0).";
SlewRateBound::usage = "SlewRateBound is an option specifying a maximum amplitude change rate per second (Abs[u[[j+1, k]] - u[[j, k]]] / (0.5*(dt[[j]] + dt[[j+1]])) <= Smax) in ControlUnits/s, as None, a positive scalar, or a per-channel list.";
FISTARestart::usage = "FISTARestart is an option for Optimizer -> \"FISTA\" selecting adaptive restart mode: \"Gradient\" (default), \"Function\", \"Both\", or None.";
MultiStart::usage = "MultiStart is an option for GRAPEOptimize specifying 1 (single start), an integer n >= 2 of randomized starts, or a list of explicit initial control specifications, returning the highest-objective solution.";
SoftMinSchedule::usage = "SoftMinSchedule is an option for GRAPEOptimize on ensembles specifying a continuation sequence of SoftMin beta parameters (None, Automatic, or a list of non-negative numbers with 0 denoting Mean aggregation).";
TimeStepSchedule::usage = "TimeStepSchedule is an option for GRAPEOptimize specifying a coarse-to-fine list of slice counts {nt1, nt2, ...} preserving TotalDuration with warm-started waveform prolongation.";
MonitorInterval::usage = "MonitorInterval is an option for GRAPEOptimize specifying StepMonitor cadence: a positive integer k (every k completed outer iterations) or {\"Time\", seconds}.";
InitialControls::usage = "InitialControls is an option specifying Automatic, \"Zero\", an explicit array, or a previous OptimalControlResult for warm restarts.";
InitialScale::usage = "InitialScale is an option controlling the dimensionless initial random flip-angle amplitude scale when InitialControls -> Automatic.";
FrequencyUnits::usage = "FrequencyUnits is an option specifying Hamiltonian units: \"rad/s\" (default, angular frequency), \"Hz\", \"kHz\", or \"MHz\".";
ControlUnits::usage = "ControlUnits is an option specifying control-channel units: Automatic (inherits FrequencyUnits), \"rad/s\", \"Hz\", \"kHz\", or \"MHz\".";
Optimizer::usage = "Optimizer is an option selecting the numerical optimizer: Automatic, \"LBFGS\" (default for smooth problems), \"BFGS\", \"ConjugateGradient\", \"SteepestAscent\", \"FISTA\", or an advanced registered optimizer (\"RegularizedLBFGS\", \"NewtonRFO\", \"NewtonShift\", \"NewtonModifiedCholesky\", \"NewtonSpectral\", \"TrustRegionNewtonCG\", \"TrustRegionGLTR\", \"AdaptiveCubicRegularization\", \"HybridLBFGSNewtonCG\", \"HybridLBFGSNewtonRFO\").";
LBFGSMemory::usage = "LBFGSMemory is an option specifying the history depth m used by the \"LBFGS\" optimizer (default 15).";
LineSearch::usage = "LineSearch is an option for GRAPEOptimize selecting Automatic (StrongWolfe when unconstrained, projected Armijo with quadratic interpolation when constrained), \"StrongWolfe\", \"Armijo\", \"HagerZhang\", \"ZhangHagerNonmonotone\", or \"GrippoNonmonotone\".";
WolfeParameters::usage = "WolfeParameters is an option for GRAPEOptimize specifying {c1, c2} for Strong-Wolfe line search (0 < c1 < c2 < 1), or Automatic ({10^-4, 0.9} for LBFGS/BFGS/SteepestAscent, {10^-4, 0.4} for ConjugateGradient).";
MaxLineSearchEvaluations::usage = "MaxLineSearchEvaluations is an option for GRAPEOptimize specifying the maximum number of trial evaluations per line search (default 20).";
TargetFidelity::usage = "TargetFidelity is an option specifying a stopping fidelity threshold (default 0.9999, or None to disable early stopping).";
GradientTolerance::usage = "GradientTolerance is an option specifying the scaled projected gradient-norm stopping tolerance (default 10^-6).";
FunctionTolerance::usage = "FunctionTolerance is an option specifying the relative objective-change stopping tolerance (default 10^-11).";
StepTolerance::usage = "StepTolerance is an option specifying the relative step-norm stopping tolerance (default 10^-10).";
PowerPenalty::usage = "PowerPenalty is an option specifying a quadratic RF power regularization coefficient (default 0).";
SmoothnessPenalty::usage = "SmoothnessPenalty is an option specifying a temporal roughness regularization coefficient on adjacent slices (default 0).";
RecordHistory::usage = "RecordHistory is an option controlling whether iteration-by-iteration objective and gradient histories are stored (default True).";
IncludeDrift::usage = "IncludeDrift is an option for OptimalControlEventList and OptimalControlGenerator controlling whether the drift Hamiltonian is included in exported generators (default True).";
TimeSteps::usage = "TimeSteps is an option specifying the number of piecewise-constant slices when timing is given as a scalar total duration.";
DriftIntegration::usage = "DriftIntegration is an option for OptimalControlProblem and GRAPEOptimize selecting numerical integration for time-dependent drift Hamiltonians: Automatic (\"CF4\" for continuous Hamiltonians, \"Midpoint\" for explicit slice lists), \"Midpoint\", {\"SubSlices\", nSub}, \"CF4\" or \"GaussLegendre2\", or {\"CF4\", nSub}.";
OptimizeDurations::usage = "OptimizeDurations is an option for OptimalControlProblem and GRAPEOptimize enabling joint optimization of durations: False (default), \"Total\" (single scaling variable preserving relative slice proportions), True or \"All\", a list of slice indices {j1, ...}, \"DelaysOnly\", or {\"DelaysOnly\", {j1, ...}}.";
DurationBounds::usage = "DurationBounds is an option for OptimalControlProblem and GRAPEOptimize specifying duration limits when OptimizeDurations is enabled: Automatic, {dtMin, dtMax} (or {tMin, tMax} for \"Total\"), or per-active-slice bounds {{min1, max1}, ...}.";
DurationPenalty::usage = "DurationPenalty is an option for OptimalControlProblem, OptimalControlEnsemble, and GRAPEOptimize specifying a dimensionless linear total-duration penalty weight lam (penalizing lam * T / Tref with Tref = initial TotalDuration) or {lam, Tref} (default 0).";
FixedTotalDuration::usage = "FixedTotalDuration is an option for OptimalControlProblem and GRAPEOptimize constraining Total[TimeSteps] to remain fixed (True or positive target duration) during duration optimization via exact bounded-simplex projection.";

OptimalControlProblem::badtiming = "Invalid timing specification `1`. Expected positive total duration, {duration, nSteps}, or a list of positive slice durations.";
OptimalControlProblem::badctrl = "Control Hamiltonians must be non-empty valid SpinDynamica operators, superoperators, or square numeric matrices.";
OptimalControlProblem::nonherm = "Hilbert-space drift and control Hamiltonians must be Hermitian (`1`); use PropagateSpace -> \"Liouville\" for relaxation/dissipative superoperators.";
OptimalControlProblem::badbasis = "Invalid ControlParameterization `1` or ControlWindow `2`.";
OptimalControlProblem::badunit = "Unknown frequency unit `1`; expected \"rad/s\", \"Hz\", \"kHz\", \"MHz\", or a positive number.";
OptimalControlProblem::badmethod = "Method `1` is incompatible with problem space `2`.";
OptimalControlProblem::badbounds = "Invalid ControlBounds or AmplitudeBound `1`.";
OptimalControlProblem::badnormbound = "Invalid AmplitudeNormBound `1` / SlewRateBound `2`, or infeasible box/disk intersection.";
OptimalControlProblem::badrestart = "Invalid FISTARestart `1`; expected \"Gradient\", \"Function\", \"Both\", or None.";
OptimalControlProblem::nonsmoothopt = "Nonsmooth regularizers, slew bounds, or constrained basis parameterizations require Optimizer -> \"FISTA\"; switching automatically from `1`.";
OptimalControlProblem::badopt = "Unknown Optimizer `1`; expected Automatic, \"LBFGS\", \"BFGS\", \"ConjugateGradient\", \"SteepestAscent\", \"FISTA\", or a registered optimizer specification.";
OptimalControlProblem::unconstronly = "Optimizer `1` only supports unconstrained optimization cleanly and cannot be used with active constraints without explicit experimental opt-in (\"AllowProjected\" -> True).";
OptimalControlProblem::unsupported2nd = "Optimizer `1` requires exact second-order derivatives, which are unavailable for this problem: `2`";
OptimalControlProblem::badls = "Invalid LineSearch `1`; expected Automatic, \"StrongWolfe\", \"Armijo\", \"HagerZhang\", \"ZhangHagerNonmonotone\", \"GrippoNonmonotone\", or a valid {name, opts} tuple.";
OptimalControlProblem::badwolfe = "Invalid WolfeParameters `1`; expected Automatic or {c1, c2} with 0 < c1 < c2 < 1.";
OptimalControlProblem::badmonitor = "Invalid MonitorInterval `1`; expected a positive integer or {\"Time\", positiveSeconds}.";
OptimalControlProblem::badmultistart = "Invalid MultiStart specification `1`; expected 1, None, False, an integer >= 2, or a non-empty list of initial controls.";
OptimalControlProblem::badschedule = "Invalid TimeStepSchedule `1`; expected None, False, or a non-empty list of positive integers.";
OptimalControlProblem::baddriftint = "Invalid DriftIntegration specification `1`; expected Automatic, \"Midpoint\", {\"SubSlices\", n}, \"CF4\", \"GaussLegendre2\", or {\"CF4\", n}.";
OptimalControlProblem::badt0 = "Invalid InitialTimePoint `1`; expected a finite real number.";
OptimalControlProblem::baddrift = "Invalid time-dependent drift specification `1`.";
OptimalControlProblem::nodrifttime = "OptimizeDurations cannot be used with explicit per-slice drift lists because drift values at shifted slice times are unspecified; use a continuous time-dependent drift or static drift.";
OptimalControlProblem::baddur = "Invalid OptimizeDurations (`1`), DurationBounds (`2`), or FixedTotalDuration (`3`) specification.";
OptimalControlProblem::badrun = "Invalid RunningCost specification `1`.";
OptimalControlResult::nopseq = "OptimalControlPulseSequence requires 1 or 2 Cartesian control channels on `1`.";
ExportOptimalControlPulse::nonuniform = "Exporting non-uniform TimeSteps to `1` (`2` slices) without Resample -> nPoints; values assume uniform dwell per line.";
ExportOptimalControlPulse::badfmt = "Unknown ExportOptimalControlPulse format `1`; expected \"Bruker\", \"SIMPSON\", or \"CSV\".";
ExportOptimalControlPulse::badchan = "Invalid Channel or ChannelStructure specification `1`.";

$OCSolverRegistry::usage = "$OCSolverRegistry is an internal hook Association used to dispatch OptimalControlEnsemble evaluations cleanly.";

Begin["`Private`"];

$OCSolverRegistry = <||>;

commonOptions = {
    Basis -> Automatic,
    OperatorBasis -> Automatic,
    PropagateSpace -> Automatic,
    Method -> Automatic,
    Optimizer -> Automatic,
    ControlParameterization -> "PiecewiseConstant",
    ControlWindow -> None,
    ControlBounds -> None,
    AmplitudeBound -> None,
    AmplitudeNormBound -> None,
    SlewRateBound -> None,
    FISTARestart -> "Gradient",
    InitialControls -> Automatic,
    InitialScale -> 0.5,
    FrequencyUnits -> "rad/s",
    ControlUnits -> Automatic,
    TimeSteps -> 50,
    InitialTimePoint -> 0.,
    DriftIntegration -> Automatic,
    RunningCost -> None,
    OptimizeDurations -> False,
    DurationBounds -> Automatic,
    DurationPenalty -> 0.,
    FixedTotalDuration -> False,
    LBFGSMemory -> 15,
    LineSearch -> Automatic,
    WolfeParameters -> Automatic,
    MaxLineSearchEvaluations -> 20,
    MaxIterations -> 200,
    TargetFidelity -> 0.9999,
    GradientTolerance -> 10^-6,
    FunctionTolerance -> 10^-11,
    StepTolerance -> 10^-10,
    PowerPenalty -> 0.,
    SmoothnessPenalty -> 0.,
    L1Penalty -> 0.,
    TVPenalty -> 0.,
    MultiStart -> 1,
    SoftMinSchedule -> None,
    TimeStepSchedule -> None,
    StepMonitor -> None,
    MonitorInterval -> 1,
    RecordHistory -> True,
    RandomSeed -> Automatic,
    TimeConstraint -> Infinity,
    Verbose -> False
};

Options[OptimalControlProblem] = commonOptions;
Options[GRAPEOptimize] = Join[commonOptions, {Parallel -> False}];
Options[OptimalControlOptimize] = Options[GRAPEOptimize];
Options[VerifyOptimalControlGradient] = {StepSize -> 10^-7, Tolerance -> 10^-5, Method -> Automatic};
Options[OptimalControlEventList] = {IncludeDrift -> True, Chronology -> Automatic, InitialTimePoint -> Automatic};
Options[OptimalControlGenerator] = {IncludeDrift -> True, InitialTimePoint -> Automatic};
Options[OptimalControlPulseSequence] = {IncludeDrift -> False, ChannelStructure -> Automatic};
Options[OptimalControlPulseTable] = {ChannelStructure -> Automatic, PhaseUnits -> "Degrees", AmplitudeScaling -> "Percent"};
Options[ExportOptimalControlPulse] = {Channel -> Automatic, ChannelStructure -> Automatic, Resample -> Automatic, Title -> "SpinDynamica Optimal Control Pulse"};
Options[VerifyOptimalControlSimulation] = {Tolerance -> 10^-10, Chronology -> Automatic, InitialTimePoint -> Automatic};

(* Resolve unit conversion factors to SpinDynamica's rad/s baseline *)
resolveUnitFactor[spec_] := Switch[spec,
    "rad/s" | "AngularFrequency" | "RadiansPerSecond", 1.,
    "Hz" | "Frequency" | "Hertz", 2. * Pi,
    "kHz", 2000. * Pi,
    "MHz", 2.0*^6 * Pi,
    _?NumericQ /; TrueQ[spec > 0], N[spec],
    _, Message[OptimalControlProblem::badunit, spec]; $Failed
];

resolveTimingVector[timing_, defaultNt_Integer] := Module[{nt},
    Which[
        NumericQ[timing] && TrueQ[timing > 0],
            If[!IntegerQ[defaultNt] || defaultNt < 1,
                Message[OptimalControlProblem::badtiming, timing]; Return[$Failed]];
            Developer`ToPackedArray[ConstantArray[N[timing / defaultNt], defaultNt]],
        MatchQ[timing, {_?NumericQ, _Integer}] && TrueQ[timing[[1]] > 0] && timing[[2]] >= 1,
            nt = timing[[2]];
            Developer`ToPackedArray[ConstantArray[N[timing[[1]] / nt], nt]],
        VectorQ[timing, NumericQ] && Length[timing] >= 1 && AllTrue[timing, TrueQ[# > 0] &],
            Developer`ToPackedArray[N[timing]],
        True,
            Message[OptimalControlProblem::badtiming, timing]; $Failed
    ]
];

toGenMatrix[expr_, space_String, basis_, opbasis_, dim_Integer] := Module[{mat},
    If[MatchQ[expr, _NullOperator | _NullSuperoperator] || (NumericQ[expr] && TrueQ[expr == 0]),
        Return[ToPackedComplexMatrix[ConstantArray[0., {dim, dim}]]]
    ];
    mat = If[space === "Hilbert",
        Which[
            OperatorQ[expr], MatrixRepresentation[expr, basis],
            MatrixQ[expr, NumericQ], expr,
            True, $Failed
        ],
        Which[
            SuperoperatorQ[expr], SuperoperatorMatrixRepresentation[expr, opbasis, basis],
            OperatorQ[expr], SuperoperatorMatrixRepresentation[CommutationSuperoperator[-I expr], opbasis, basis],
            MatrixQ[expr, NumericQ], expr,
            True, $Failed
        ]
    ];
    If[mat === $Failed || !MatrixQ[mat, NumericQ] || Dimensions[mat] =!= {dim, dim},
        $Failed,
        ToPackedComplexMatrix[mat]
    ]
];

toSpinDynamicaGen[expr_, unitFactor_?NumericQ] := Module[{scaled},
    If[expr === 0 || MatchQ[expr, _NullOperator | _NullSuperoperator], Return[NullOperator[]]];
    scaled = unitFactor * expr;
    scaled
];

(* Time-dependent & MAS drift classification, standardization, hoisting, and quadrature *)
isStaticZeroDriftQ[expr_] :=
    MatchQ[expr, _NullOperator | _NullSuperoperator] || (NumericQ[expr] && TrueQ[expr == 0]);

isStaticKnownDriftQ[expr_] :=
    isStaticZeroDriftQ[expr] || OperatorQ[expr] || SuperoperatorQ[expr] || MatrixQ[expr, NumericQ];

isSliceListDriftQ[expr_] := ListQ[expr] && !MatrixQ[expr, NumericQ];

isContinuousDriftQ[expr_] :=
    !isStaticKnownDriftQ[expr] && !isSliceListDriftQ[expr] && !MatchQ[expr, _String | _Rule];

prepareContinuousDrift[drift_] := Module[{res},
    If[!FreeQ[drift, _PeriodicFunction | _TimeDependentGeneratorFunction],
        res = Quiet[drift /. {
            pf_PeriodicFunction :> StandardizeGenerator[pf, ImplementMatrixRepresentation -> False],
            td_TimeDependentGeneratorFunction :> StandardizeGenerator[td, ImplementMatrixRepresentation -> False]
        }];
        Return[res]
    ];
    If[MatchQ[drift, _Function] && ListQ[First[drift]] && Length[First[drift]] >= 2,
        res = Quiet[Check[StandardizeGenerator[drift, ImplementMatrixRepresentation -> False], drift]];
        Return[res]
    ];
    drift
];

evalContinuousExprAt[prepDrift_, tVal_, t0_: 0., tEnd_: 1.] := Module[{res},
    If[!FreeQ[prepDrift, _TimeDependentGeneratorFunction],
        res = Quiet[
            prepDrift /. td_TimeDependentGeneratorFunction :> Quiet[
                Check[td[tVal], td[tVal, t0, tEnd], {TimeDependentGeneratorFunction::misst}],
                {TimeDependentGeneratorFunction::misst}
            ],
            {TimeDependentGeneratorFunction::misst}
        ];
        Return[res]
    ];
    If[MatchQ[prepDrift, _Function],
        If[ListQ[First[prepDrift]],
            Switch[Length[First[prepDrift]],
                2, Return[prepDrift[tVal, 0.5 * (t0 + tEnd)]],
                3, Return[prepDrift[tVal, t0, tEnd]],
                _, Return[prepDrift[tVal]]
            ],
            Return[prepDrift[tVal]]
        ]
    ];
    prepDrift[tVal]
];

driftContainsSuperoperatorQ[drift_, t0_: 0., totalDur_: 1.] := Module[{prep, probe, tEnd},
    Which[
        isStaticZeroDriftQ[drift], False,
        SuperoperatorQ[drift], True,
        OperatorQ[drift] || MatrixQ[drift, NumericQ], False,
        isSliceListDriftQ[drift], AnyTrue[drift, SuperoperatorQ],
        !FreeQ[drift, _?SuperoperatorQ], True,
        isContinuousDriftQ[drift],
            tEnd = N[t0] + N[totalDur];
            prep = prepareContinuousDrift[drift];
            probe = Quiet[Check[evalContinuousExprAt[prep, N[t0] + 0.5 * N[totalDur], N[t0], tEnd], $Failed]];
            TrueQ[SuperoperatorQ[probe] || !FreeQ[probe, _?SuperoperatorQ]],
        True, False
    ]
];

parseDriftIntegrationSpec[opt_, isListDriftQ_: False] := Module[{nSub},
    If[TrueQ[isListDriftQ],
        Return[Switch[opt,
            Automatic | "Midpoint", <|"Method" -> "Midpoint", "SubSlices" -> 1, "NumStages" -> 1|>,
            {"SubSlices", 1}, <|"Method" -> "Midpoint", "SubSlices" -> 1, "NumStages" -> 1|>,
            _, Message[OptimalControlProblem::baddriftint, opt]; $Failed
        ]]
    ];
    Which[
        opt === Automatic,
            <|"Method" -> "CF4", "SubSlices" -> 1, "NumStages" -> 2|>,
        opt === "Midpoint",
            <|"Method" -> "Midpoint", "SubSlices" -> 1, "NumStages" -> 1|>,
        MatchQ[opt, {"SubSlices", _Integer}] && opt[[2]] >= 1,
            nSub = opt[[2]];
            <|"Method" -> "Midpoint", "SubSlices" -> nSub, "NumStages" -> nSub|>,
        opt === "CF4" || opt === "GaussLegendre2",
            <|"Method" -> "CF4", "SubSlices" -> 1, "NumStages" -> 2|>,
        MatchQ[opt, {"CF4" | "GaussLegendre2", _Integer}] && opt[[2]] >= 1,
            nSub = opt[[2]];
            <|"Method" -> "CF4", "SubSlices" -> nSub, "NumStages" -> 2 * nSub|>,
        True,
            Message[OptimalControlProblem::baddriftint, opt]; $Failed
    ]
];

buildFineDriftSchedule[evalMat_, dtVec_List, t0_?NumericQ, intSpec_Association, space_String] := Module[
    {nt = Length[dtVec], dMethod = intSpec["Method"], nSub = intSpec["SubSlices"],
     nStages = intSpec["NumStages"], nFine, dtFine, hMats, stageMeta, tEdges,
     c1, c2, wPlus, wMinus, sIdx = 0, j, sub, tA, h, t1, t2, h1, h2, hMid, hS1, hS2, hermErr},
    nFine = nt * nStages;
    dtFine = ConstantArray[0., nFine];
    hMats = Table[Null, {nFine}];
    stageMeta = Table[Null, {nFine}];
    tEdges = N[t0] + Prepend[Accumulate[N[dtVec]], 0.];
    c1 = 0.5 - Sqrt[3.] / 6.;
    c2 = 0.5 + Sqrt[3.] / 6.;
    wPlus = 0.5 + Sqrt[3.] / 3.;
    wMinus = 0.5 - Sqrt[3.] / 3.;

    Do[
        h = N[dtVec[[j]]] / nSub;
        Do[
            tA = tEdges[[j]] + sub * h;
            If[dMethod === "Midpoint",
                sIdx++;
                t1 = tA + 0.5 * h;
                hMid = evalMat[t1];
                If[hMid === $Failed || !MatrixQ[hMid, NumericQ],
                    Message[OptimalControlProblem::baddrift, "non-numeric evaluation"]; Return[$Failed, Module]
                ];
                If[space === "Hilbert",
                    hermErr = Max[Abs[hMid - ConjugateTranspose[hMid]]] / Max[1., Max[Abs[hMid]]];
                    If[hermErr > 10^-8, Message[OptimalControlProblem::nonherm, hermErr]; Return[$Failed, Module]];
                    hMid = HermitianSymmetrize[hMid]
                ];
                dtFine[[sIdx]] = h;
                hMats[[sIdx]] = hMid;
                stageMeta[[sIdx]] = <|
                    "SliceIndex" -> j, "Alpha" -> 1. / nSub,
                    "Nodes" -> {t1}, "Weights" -> {1.}, "ShiftPrefixes" -> {(sub + 0.5) / nSub}
                |>,
                (* CF4 exact 2-stage ordering matching Propagation.m#L1493-L1499 *)
                t1 = tA + c1 * h;
                t2 = tA + c2 * h;
                h1 = evalMat[t1];
                h2 = evalMat[t2];
                If[h1 === $Failed || h2 === $Failed || !MatrixQ[h1, NumericQ] || !MatrixQ[h2, NumericQ],
                    Message[OptimalControlProblem::baddrift, "non-numeric evaluation"]; Return[$Failed, Module]
                ];
                If[space === "Hilbert",
                    hermErr = Max[
                        Max[Abs[h1 - ConjugateTranspose[h1]]] / Max[1., Max[Abs[h1]]],
                        Max[Abs[h2 - ConjugateTranspose[h2]]] / Max[1., Max[Abs[h2]]]
                    ];
                    If[hermErr > 10^-8, Message[OptimalControlProblem::nonherm, hermErr]; Return[$Failed, Module]];
                    h1 = HermitianSymmetrize[h1];
                    h2 = HermitianSymmetrize[h2]
                ];
                hS1 = ToPackedComplexMatrix[wPlus * h1 + wMinus * h2];
                hS2 = ToPackedComplexMatrix[wMinus * h1 + wPlus * h2];
                sIdx++;
                dtFine[[sIdx]] = 0.5 * h;
                hMats[[sIdx]] = hS1;
                stageMeta[[sIdx]] = <|
                    "SliceIndex" -> j, "Alpha" -> 0.5 / nSub,
                    "Nodes" -> {t1, t2}, "Weights" -> {wPlus, wMinus},
                    "ShiftPrefixes" -> {(sub + c1) / nSub, (sub + c2) / nSub}
                |>;
                sIdx++;
                dtFine[[sIdx]] = 0.5 * h;
                hMats[[sIdx]] = hS2;
                stageMeta[[sIdx]] = <|
                    "SliceIndex" -> j, "Alpha" -> 0.5 / nSub,
                    "Nodes" -> {t1, t2}, "Weights" -> {wMinus, wPlus},
                    "ShiftPrefixes" -> {(sub + c1) / nSub, (sub + c2) / nSub}
                |>
            ],
            {sub, 0, nSub - 1}
        ],
        {j, nt}
    ];
    <|
        "Drift" -> Developer`ToPackedArray[hMats],
        "FineTimeSteps" -> Developer`ToPackedArray[dtFine],
        "StageMeta" -> stageMeta
    |>
];

compileDriftRepresentation[drift_, space_String, basis_, opbasis_, dim_Integer, hamFactor_?NumericQ,
    dtVec_List, t0_?NumericQ, driftIntOpt_] := Module[
    {nt = Length[dtVec], totalDur = Total[dtVec], t0N = N[t0], tEndN,
     h0Raw, h0Mat, hermErr, isListQ, intSpec, h0List,
     prepDrift, pointFn, tProbe1, tProbe2, probe1, probe2,
     tSym, symExpr, symRaw, symDense, hoistOK = False, mSub1, mSub2, relTol,
     dSym, dSub1, dSub2, fd1, fd2, hFD, derivOK = False, evalMatFn, derivMatFn, sched},
    tEndN = t0N + N[totalDur];
    If[isStaticKnownDriftQ[drift],
        intSpec = parseDriftIntegrationSpec[driftIntOpt, False];
        If[intSpec === $Failed, Return[$Failed]];
        h0Raw = toGenMatrix[drift, space, basis, opbasis, dim];
        If[h0Raw === $Failed, Message[OptimalControlProblem::baddrift, drift]; Return[$Failed]];
        h0Mat = ToPackedComplexMatrix[hamFactor * h0Raw];
        If[space === "Hilbert",
            hermErr = Max[Abs[h0Mat - ConjugateTranspose[h0Mat]]] / Max[1., Max[Abs[h0Mat]]];
            If[hermErr > 10^-8, Message[OptimalControlProblem::nonherm, hermErr]; Return[$Failed]];
            h0Mat = HermitianSymmetrize[h0Mat]
        ];
        Return[<|
            "Drift" -> h0Mat,
            "FineTimeSteps" -> Developer`ToPackedArray[N[dtVec]],
            "NumStages" -> 1,
            "IsTimeDependentDrift" -> False,
            "DriftKind" -> "Static",
            "DriftIntegration" -> "Static",
            "InitialTimePoint" -> t0N,
            "DriftEvaluator" -> None,
            "DriftDerivEvaluator" -> None,
            "StageMeta" -> None
        |>]
    ];

    isListQ = isSliceListDriftQ[drift];
    If[isListQ,
        If[Length[drift] =!= nt, Message[OptimalControlProblem::baddrift, drift]; Return[$Failed]];
        intSpec = parseDriftIntegrationSpec[driftIntOpt, True];
        If[intSpec === $Failed, Return[$Failed]];
        h0List = Table[
            h0Raw = toGenMatrix[drift[[j]], space, basis, opbasis, dim];
            If[h0Raw === $Failed, Message[OptimalControlProblem::baddrift, drift[[j]]]; Return[$Failed, Module]];
            h0Mat = ToPackedComplexMatrix[hamFactor * h0Raw];
            If[space === "Hilbert",
                hermErr = Max[Abs[h0Mat - ConjugateTranspose[h0Mat]]] / Max[1., Max[Abs[h0Mat]]];
                If[hermErr > 10^-8, Message[OptimalControlProblem::nonherm, hermErr]; Return[$Failed, Module]];
                h0Mat = HermitianSymmetrize[h0Mat]
            ];
            h0Mat,
            {j, nt}
        ];
        Return[<|
            "Drift" -> Developer`ToPackedArray[h0List],
            "FineTimeSteps" -> Developer`ToPackedArray[N[dtVec]],
            "NumStages" -> 1,
            "IsTimeDependentDrift" -> True,
            "DriftKind" -> "SliceList",
            "DriftIntegration" -> "Midpoint",
            "InitialTimePoint" -> t0N,
            "DriftEvaluator" -> None,
            "DriftDerivEvaluator" -> None,
            "StageMeta" -> None
        |>]
    ];

    If[!isContinuousDriftQ[drift],
        Message[OptimalControlProblem::baddrift, drift]; Return[$Failed]
    ];
    intSpec = parseDriftIntegrationSpec[driftIntOpt, False];
    If[intSpec === $Failed, Return[$Failed]];

    prepDrift = prepareContinuousDrift[drift];
    pointFn = With[{pr = prepDrift, sp = space, b = basis, ob = opbasis, d = dim,
                    tStart = t0N, tEnd = tEndN, hf = N[hamFactor]},
        Function[tv,
            With[{ev = Quiet[Check[evalContinuousExprAt[pr, N[tv], tStart, tEnd], $Failed]]},
                If[ev === $Failed, $Failed,
                    With[{m = toGenMatrix[ev, sp, b, ob, d]},
                        If[m === $Failed, $Failed, ToPackedComplexMatrix[hf * m]]
                    ]
                ]
            ]
        ]
    ];

    tProbe1 = t0N + 0.17 * N[totalDur];
    tProbe2 = t0N + 0.73 * N[totalDur];
    probe1 = pointFn[tProbe1];
    probe2 = pointFn[tProbe2];
    If[probe1 === $Failed || probe2 === $Failed,
        Message[OptimalControlProblem::baddrift, drift]; Return[$Failed]
    ];

    (* Hoist symbolic SpinDynamica matrix representation once and verify exact match *)
    tSym = Unique["ocDriftT"];
    symExpr = Quiet[Check[evalContinuousExprAt[prepDrift, tSym, t0N, tEndN], $Failed]];
    If[symExpr =!= $Failed,
        symRaw = Quiet[Check[
            If[space === "Hilbert",
                Which[
                    isStaticZeroDriftQ[symExpr], ConstantArray[0., {dim, dim}],
                    OperatorQ[symExpr], MatrixRepresentation[symExpr, basis],
                    MatrixQ[symExpr], symExpr,
                    True, $Failed
                ],
                Which[
                    isStaticZeroDriftQ[symExpr], ConstantArray[0., {dim, dim}],
                    SuperoperatorQ[symExpr], SuperoperatorMatrixRepresentation[symExpr, opbasis, basis],
                    OperatorQ[symExpr], SuperoperatorMatrixRepresentation[CommutationSuperoperator[-I symExpr], opbasis, basis],
                    MatrixQ[symExpr], symExpr,
                    True, $Failed
                ]
            ],
            $Failed
        ]];
        If[symRaw =!= $Failed && MatrixQ[symRaw] && Dimensions[symRaw] === {dim, dim},
            symDense = N[hamFactor] * Normal[symRaw];
            mSub1 = Quiet[Check[ToPackedComplexMatrix[N[symDense /. tSym -> tProbe1]], $Failed]];
            mSub2 = Quiet[Check[ToPackedComplexMatrix[N[symDense /. tSym -> tProbe2]], $Failed]];
            relTol = 10^-11;
            If[MatrixQ[mSub1, NumericQ] && MatrixQ[mSub2, NumericQ] &&
               Max[Abs[mSub1 - probe1]] <= relTol * Max[1., Max[Abs[probe1]]] &&
               Max[Abs[mSub2 - probe2]] <= relTol * Max[1., Max[Abs[probe2]]],
                hoistOK = True
            ]
        ]
    ];

    If[hoistOK,
        evalMatFn = With[{sm = symDense, ts = tSym},
            Function[tv, ToPackedComplexMatrix[N[sm /. ts -> N[tv]]]]
        ];
        dSym = Quiet[Check[D[symDense, tSym], $Failed]];
        hFD = Max[10^-11, 10^-6 * Min[N[totalDur], Mean[N[dtVec]]]];
        If[MatrixQ[dSym] && Dimensions[dSym] === {dim, dim},
            dSub1 = Quiet[Check[ToPackedComplexMatrix[N[dSym /. tSym -> tProbe1]], $Failed]];
            dSub2 = Quiet[Check[ToPackedComplexMatrix[N[dSym /. tSym -> tProbe2]], $Failed]];
            fd1 = (evalMatFn[tProbe1 + hFD] - evalMatFn[tProbe1 - hFD]) / (2. * hFD);
            fd2 = (evalMatFn[tProbe2 + hFD] - evalMatFn[tProbe2 - hFD]) / (2. * hFD);
            If[MatrixQ[dSub1, NumericQ] && MatrixQ[dSub2, NumericQ] &&
               Max[Abs[dSub1 - fd1]] <= 10^-4 * Max[1., Max[Abs[fd1]]] &&
               Max[Abs[dSub2 - fd2]] <= 10^-4 * Max[1., Max[Abs[fd2]]],
                derivOK = True
            ]
        ];
        If[derivOK,
            derivMatFn = With[{dsm = dSym, ts = tSym, spQ = (space === "Hilbert")},
                Function[tv,
                    With[{m = ToPackedComplexMatrix[N[dsm /. ts -> N[tv]]]},
                        If[spQ, HermitianSymmetrize[m], m]
                    ]
                ]
            ],
            derivMatFn = With[{ef = evalMatFn, h = hFD, spQ = (space === "Hilbert")},
                Function[tv,
                    With[{m = ToPackedComplexMatrix[(ef[N[tv] + h] - ef[N[tv] - h]) / (2. * h)]},
                        If[spQ, HermitianSymmetrize[m], m]
                    ]
                ]
            ]
        ],
        evalMatFn = pointFn;
        hFD = Max[10^-11, 10^-6 * Min[N[totalDur], Mean[N[dtVec]]]];
        derivMatFn = With[{ef = evalMatFn, h = hFD, spQ = (space === "Hilbert")},
            Function[tv,
                With[{m = ToPackedComplexMatrix[(ef[N[tv] + h] - ef[N[tv] - h]) / (2. * h)]},
                    If[spQ, HermitianSymmetrize[m], m]
                ]
            ]
        ]
    ];

    sched = buildFineDriftSchedule[evalMatFn, dtVec, t0N, intSpec, space];
    If[sched === $Failed, Return[$Failed]];

    <|
        "Drift" -> sched["Drift"],
        "FineTimeSteps" -> sched["FineTimeSteps"],
        "NumStages" -> intSpec["NumStages"],
        "IsTimeDependentDrift" -> True,
        "DriftKind" -> "Continuous",
        "DriftIntegration" -> Which[
            intSpec["Method"] === "CF4" && intSpec["SubSlices"] > 1, {"CF4", intSpec["SubSlices"]},
            intSpec["Method"] === "CF4", "CF4",
            intSpec["SubSlices"] > 1, {"SubSlices", intSpec["SubSlices"]},
            True, "Midpoint"
        ],
        "DriftIntegrationSpec" -> intSpec,
        "InitialTimePoint" -> t0N,
        "DriftEvaluator" -> evalMatFn,
        "DriftDerivEvaluator" -> derivMatFn,
        "StageMeta" -> sched["StageMeta"]
    |>
];

(* Bounds resolution on parameters of shape {nParams, mCtrl} *)
resolveBoundMatrices[ctrlBounds_, ampBound_, nParams_Integer, mCtrl_Integer] := Module[{spec},
    spec = If[ctrlBounds =!= None && ctrlBounds =!= Automatic, ctrlBounds,
        If[ampBound =!= None && ampBound =!= Automatic,
            If[NumericQ[ampBound] && TrueQ[ampBound > 0], {-N[ampBound], N[ampBound]}, $Failed],
            None
        ]
    ];
    If[spec === None || spec === False,
        Return[{ConstantArray[-Infinity, nParams * mCtrl], ConstantArray[Infinity, nParams * mCtrl], False}]
    ];
    Which[
        MatchQ[spec, {_?NumericQ, _?NumericQ}] && TrueQ[spec[[1]] <= spec[[2]]],
            {
                ConstantArray[N[spec[[1]]], nParams * mCtrl],
                ConstantArray[N[spec[[2]]], nParams * mCtrl],
                True
            },
        MatchQ[spec, {{_?NumericQ, _?NumericQ} ..}] && Length[spec] === mCtrl &&
                AllTrue[spec, TrueQ[#[[1]] <= #[[2]]] &],
            {
                Flatten[ConstantArray[N[spec[[All, 1]]], nParams]],
                Flatten[ConstantArray[N[spec[[All, 2]]], nParams]],
                True
            },
        True,
            Message[OptimalControlProblem::badbounds, spec]; $Failed
    ]
];

(* Polar/disk RF group norm bounds, slew-rate bounds, & composite constraint spec *)
resolveNormBounds[ampNormOpt_, mCtrl_Integer] := Module[
    {groups = {}, lims = {}, parsed, entry, chSpec, wLim, gIdx, flatIdx, disjointQ},
    If[ampNormOpt === None || ampNormOpt === False || ampNormOpt === Automatic,
        Return[<|"Groups" -> {}, "Limits" -> {}, "HasDisk" -> False, "Disjoint" -> True|>]
    ];
    If[NumericQ[ampNormOpt] && TrueQ[ampNormOpt > 0],
        Which[
            mCtrl === 1,
                groups = {{1}}; lims = {N[ampNormOpt]},
            EvenQ[mCtrl],
                groups = Partition[Range[mCtrl], 2];
                lims = ConstantArray[N[ampNormOpt], mCtrl / 2],
            True,
                Message[OptimalControlProblem::badnormbound, ampNormOpt, None]; Return[$Failed]
        ];
        Return[<|"Groups" -> groups, "Limits" -> lims, "HasDisk" -> True, "Disjoint" -> True|>]
    ];
    If[!ListQ[ampNormOpt] || ampNormOpt === {},
        Message[OptimalControlProblem::badnormbound, ampNormOpt, None]; Return[$Failed]
    ];
    Do[
        entry = ampNormOpt[[idx]];
        parsed = Which[
            MatchQ[entry, Rule[_Integer | _List, _?NumericQ]] && TrueQ[entry[[2]] > 0],
                {If[ListQ[entry[[1]]], entry[[1]], {entry[[1]]}], N[entry[[2]]]},
            MatchQ[entry, {_List, _?NumericQ}] && TrueQ[entry[[2]] > 0],
                {entry[[1]], N[entry[[2]]]},
            MatchQ[entry, {_Integer, _?NumericQ}] && TrueQ[entry[[2]] > 0],
                {{entry[[1]]}, N[entry[[2]]]},
            True,
                $Failed
        ];
        If[parsed === $Failed, Message[OptimalControlProblem::badnormbound, ampNormOpt, None]; Return[$Failed, Do]];
        {gIdx, wLim} = parsed;
        If[gIdx === {} || !AllTrue[gIdx, IntegerQ[#] && 1 <= # <= mCtrl &] || Length[DeleteDuplicates[gIdx]] =!= Length[gIdx],
            Message[OptimalControlProblem::badnormbound, ampNormOpt, None]; Return[$Failed, Do]
        ];
        AppendTo[groups, gIdx];
        AppendTo[lims, wLim],
        {idx, Length[ampNormOpt]}
    ];
    If[Length[groups] =!= Length[ampNormOpt], Return[$Failed]];
    flatIdx = Flatten[groups];
    disjointQ = (Length[flatIdx] === Length[DeleteDuplicates[flatIdx]]);
    <|"Groups" -> groups, "Limits" -> lims, "HasDisk" -> True, "Disjoint" -> disjointQ|>
];

resolveSlewBounds[slewOpt_, dtVec_List, mCtrl_Integer] := Module[{nt, rates, dtAdj, deltaMat},
    nt = Length[dtVec];
    If[slewOpt === None || slewOpt === False || slewOpt === Automatic,
        Return[<|"Rates" -> ConstantArray[Infinity, mCtrl], "DeltaMat" -> None, "HasSlew" -> False|>]
    ];
    rates = Which[
        NumericQ[slewOpt] && TrueQ[slewOpt > 0],
            ConstantArray[N[slewOpt], mCtrl],
        VectorQ[slewOpt, NumericQ] && Length[slewOpt] === mCtrl && AllTrue[slewOpt, TrueQ[# > 0] &],
            N[slewOpt],
        True,
            $Failed
    ];
    If[rates === $Failed,
        Message[OptimalControlProblem::badnormbound, None, slewOpt]; Return[$Failed]
    ];
    If[nt <= 1,
        Return[<|"Rates" -> rates, "DeltaMat" -> None, "HasSlew" -> False|>]
    ];
    dtAdj = 0.5 * (Most[dtVec] + Rest[dtVec]);
    deltaMat = Developer`ToPackedArray[Outer[Times, dtAdj, rates]];
    <|"Rates" -> rates, "DeltaMat" -> deltaMat, "HasSlew" -> True|>
];

resolveConstraintSpec[ctrlBounds_, ampBound_, ampNormOpt_, slewOpt_, dtVec_List, nParams_Integer, mCtrl_Integer] := Module[
    {nt = Length[dtVec], boxParam, lbVec, ubVec, hasBoxQ, boxSlice, lbSlice, ubSlice,
     lbMat, ubMat, finiteBoxQ, diskSpec, slewSpec, gIdx, wLim, minSq, okFeas = True},
    boxParam = resolveBoundMatrices[ctrlBounds, ampBound, nParams, mCtrl];
    If[boxParam === $Failed, Return[$Failed]];
    {lbVec, ubVec, hasBoxQ} = boxParam;
    boxSlice = resolveBoundMatrices[ctrlBounds, ampBound, nt, mCtrl];
    If[boxSlice === $Failed, Return[$Failed]];
    {lbSlice, ubSlice, hasBoxQ} = boxSlice;
    lbMat = Partition[N[lbSlice], mCtrl];
    ubMat = Partition[N[ubSlice], mCtrl];
    finiteBoxQ = TrueQ[hasBoxQ] &&
        AllTrue[Flatten[lbMat], NumericQ[#] && # > -Infinity &] &&
        AllTrue[Flatten[ubMat], NumericQ[#] && # < Infinity &];
    If[finiteBoxQ,
        lbMat = ToPackedRealMatrix[lbMat];
        ubMat = ToPackedRealMatrix[ubMat]
    ];

    diskSpec = resolveNormBounds[ampNormOpt, mCtrl];
    If[diskSpec === $Failed, Return[$Failed]];
    slewSpec = resolveSlewBounds[slewOpt, dtVec, mCtrl];
    If[slewSpec === $Failed, Return[$Failed]];

    If[TrueQ[hasBoxQ] && TrueQ[diskSpec["HasDisk"]],
        Do[
            gIdx = (diskSpec["Groups"])[[g]];
            wLim = (diskSpec["Limits"])[[g]];
            minSq = Total[(MapThread[Clip[0., {#1, #2}] &, {lbMat[[1, gIdx]], ubMat[[1, gIdx]]}])^2];
            If[minSq > (wLim + 10^-9 * Max[1., wLim])^2,
                okFeas = False; Break[]
            ],
            {g, Length[diskSpec["Groups"]]}
        ];
        If[!okFeas,
            Message[OptimalControlProblem::badnormbound, ampNormOpt, slewOpt]; Return[$Failed]
        ]
    ];

    <|
        "LowerBounds" -> Developer`ToPackedArray[lbVec],
        "UpperBounds" -> Developer`ToPackedArray[ubVec],
        "LowerBoundsMat" -> lbMat,
        "UpperBoundsMat" -> ubMat,
        "FiniteBox" -> finiteBoxQ,
        "HasBox" -> hasBoxQ,
        "NormGroups" -> diskSpec["Groups"],
        "NormLimits" -> diskSpec["Limits"],
        "DisjointDisks" -> diskSpec["Disjoint"],
        "HasDisk" -> diskSpec["HasDisk"],
        "SlewRates" -> slewSpec["Rates"],
        "SlewDeltaMat" -> slewSpec["DeltaMat"],
        "HasSlew" -> slewSpec["HasSlew"],
        "HasAnyConstraint" -> (hasBoxQ || diskSpec["HasDisk"] || slewSpec["HasSlew"]),
        "NumSlices" -> nt,
        "NumControls" -> mCtrl
    |>
];

(* Exact Condat 1-D Total Variation proximal operator *)
proxTV1D[y_List, lam_?NumericQ] := Module[
    {n = Length[y], x, k = 1, k0 = 1, kMinus = 1, kPlus = 1, vmin, vmax, umin, umax},
    If[n <= 1 || lam <= 0., Return[Developer`ToPackedArray[N[y]]]];
    x = ConstantArray[0., n];
    vmin = y[[1]] - lam; vmax = y[[1]] + lam;
    umin = lam; umax = -lam;
    While[True,
        If[k === n,
            Which[
                umin < 0.,
                    Do[x[[i]] = vmin, {i, k0, kMinus}];
                    k = k0 = kMinus = kPlus = kMinus + 1;
                    If[k <= n, vmin = y[[k]]; vmax = y[[k]] + 2. lam; umin = lam; umax = -lam],
                umax > 0.,
                    Do[x[[i]] = vmax, {i, k0, kPlus}];
                    k = k0 = kMinus = kPlus = kPlus + 1;
                    If[k <= n, vmin = y[[k]] - 2. lam; vmax = y[[k]]; umin = lam; umax = -lam],
                True,
                    With[{val = vmin + umin / (k - k0 + 1)}, Do[x[[i]] = val, {i, k0, n}]];
                    Break[]
            ];
            If[k > n, Break[]];
            Continue[]
        ];
        umin += y[[k + 1]] - vmin;
        umax += y[[k + 1]] - vmax;
        Which[
            umin < -lam,
                Do[x[[i]] = vmin, {i, k0, kMinus}];
                k = k0 = kMinus = kPlus = kMinus + 1;
                vmin = y[[k]]; vmax = y[[k]] + 2. lam; umin = lam; umax = -lam,
            umax > lam,
                Do[x[[i]] = vmax, {i, k0, kPlus}];
                k = k0 = kMinus = kPlus = kPlus + 1;
                vmin = y[[k]] - 2. lam; vmax = y[[k]]; umin = lam; umax = -lam,
            True,
                k++;
                If[umin >= lam,
                    vmin += (umin - lam) / (k - k0 + 1);
                    umin = lam;
                    kMinus = k
                ];
                If[umax <= -lam,
                    vmax += (umax + lam) / (k - k0 + 1);
                    umax = -lam;
                    kPlus = k
                ]
        ]
    ];
    Developer`ToPackedArray[x]
];

softThresholdMat[v_, tau_?NumericQ] := If[tau <= 0., v, Developer`ToPackedArray[Sign[v] * Ramp[Abs[v] - tau]]];

clampBoxMat[v_, lbMat_, ubMat_, finiteQ_?BooleanQ] :=
    If[TrueQ[finiteQ],
        Developer`ToPackedArray[v - Ramp[v - ubMat] + Ramp[lbMat - v]],
        Developer`ToPackedArray[Partition[
            MapThread[Clip[#1, {#2, #3}] &, {Flatten[v], Flatten[lbMat], Flatten[ubMat]}],
            Dimensions[v][[2]]
        ]]
    ];

projDiskGroups[uMat_, groups_List, lims_List, disjoint_?BooleanQ] := Module[
    {out = ToPackedRealMatrix[uMat], mCtrl = Dimensions[uMat][[2]], nG = Length[groups],
     sub, norms, scale, g, w, row, inc, cand, nrm, maxDiff, diff},
    If[nG === 0, Return[out]];
    If[TrueQ[disjoint],
        Do[
            g = groups[[idx]]; w = lims[[idx]];
            sub = out[[All, g]];
            norms = Sqrt[Total[sub^2, {2}]];
            scale = Map[If[# > w, w / #, 1.] &, norms];
            out[[All, g]] = scale * sub,
            {idx, nG}
        ];
        Return[Developer`ToPackedArray[out]]
    ];
    (* Overlapping channel groups: Dykstra cyclic projection per slice *)
    Do[
        row = out[[j]];
        If[AnyTrue[Table[Norm[row[[groups[[idx]]]]] > lims[[idx]] + 10^-13, {idx, nG}], TrueQ],
            inc = ConstantArray[0., {nG, mCtrl}];
            Do[
                maxDiff = 0.;
                Do[
                    g = groups[[idx]]; w = lims[[idx]];
                    cand = row[[g]] + inc[[idx, g]];
                    nrm = Norm[cand];
                    diff = If[nrm > w, (w / nrm) * cand, cand];
                    maxDiff = Max[maxDiff, Max[Abs[diff - row[[g]]]]];
                    inc[[idx, g]] = cand - diff;
                    row[[g]] = diff,
                    {idx, nG}
                ];
                If[maxDiff <= 10^-11 * Max[1., Max[lims]], Break[]],
                {40}
            ];
            out[[j]] = row
        ],
        {j, Length[out]}
    ];
    Developer`ToPackedArray[out]
];

evalNonSmoothPenalty[ctx_Association, uMat_] := Module[{lam1, lamTV, pen = 0.},
    lam1 = ctx["L1Penalty"];
    lamTV = ctx["TVPenalty"];
    If[lam1 > 0., pen += lam1 * Total[Abs[Flatten[uMat]]]];
    If[lamTV > 0. && Length[uMat] > 1, pen += lamTV * Total[Abs[Flatten[Rest[uMat] - Most[uMat]]]]];
    pen
];

computeConstraintViolation[ctx_Association, uMat_] := Module[{viol = 0., lbMat, ubMat, groups, lims, sDelta, sub, nrms},
    If[TrueQ[ctx["HasBox"]],
        lbMat = ctx["LowerBoundsMat"]; ubMat = ctx["UpperBoundsMat"];
        viol = Max[viol, Max[MapThread[Max[#2 - #1, #1 - #3, 0.] &, {Flatten[uMat], Flatten[lbMat], Flatten[ubMat]}]]]
    ];
    If[TrueQ[ctx["HasDisk"]],
        groups = ctx["NormGroups"]; lims = ctx["NormLimits"];
        Do[
            sub = uMat[[All, groups[[idx]]]];
            nrms = Sqrt[Total[sub^2, {2}]];
            viol = Max[viol, Max[nrms - lims[[idx]], 0.]],
            {idx, Length[groups]}
        ]
    ];
    If[TrueQ[ctx["HasSlew"]] && Length[uMat] > 1,
        sDelta = ctx["SlewDeltaMat"];
        viol = Max[viol, Max[Abs[Rest[uMat] - Most[uMat]] - sDelta, 0.]]
    ];
    viol
];

buildOrthonormalBasisReduction[ctrlBasis_Association] := Module[
    {nt, nParams, bMat, uSVD, sSVD, vSVD, sDiag, r, qMat, qT, sKeep, vKeep, rMat, rPinv, qBasis},
    nt = ctrlBasis["NumSlices"];
    nParams = ctrlBasis["NumParams"];
    bMat = ToPackedRealMatrix[ctrlBasis["BasisMatrix"]];
    If[nt < 1 || nParams < 1, Return[$Failed]];
    {uSVD, sSVD, vSVD} = SingularValueDecomposition[bMat, Min[nt, nParams]];
    sDiag = Diagonal[sSVD];
    If[sDiag === {} || !TrueQ[First[sDiag] > 10^-15], Return[$Failed]];
    r = Max[1, Count[sDiag, x_ /; x > 10^-10 * First[sDiag]]];
    qMat = ToPackedRealMatrix[uSVD[[All, 1 ;; r]]];
    qT = Developer`ToPackedArray[Transpose[qMat]];
    sKeep = Take[sDiag, r];
    vKeep = vSVD[[All, 1 ;; r]];
    rMat = ToPackedRealMatrix[DiagonalMatrix[sKeep] . Transpose[vKeep]];
    rPinv = ToPackedRealMatrix[vKeep . DiagonalMatrix[1. / sKeep]];
    qBasis = <|
        "BasisMatrix" -> qMat,
        "BasisTranspose" -> qT,
        "IsIdentityQ" -> False,
        "NumSlices" -> nt,
        "NumParams" -> r
    |>;
    <|"Rank" -> r, "Q" -> qMat, "QT" -> qT, "R" -> rMat, "Rpinv" -> rPinv, "QBasis" -> qBasis|>
];

compileProxContext[ctrlBasis_Association, dtVec_List, cSpec_Association, lam1_?NumericQ, lamTV_?NumericQ, allowOrtho_: True] := Module[
    {nt, nParams, mCtrl, bMat, bT, isIdent, hasBox, hasDisk, hasSlew, hasL1, hasTV,
     mode, mergeL1IntoB, useA, useB, useC, d1Mat, dbMat, dbT, btb, dbtd, rhoA, rhoB, rhoC,
     diagWeight, kMat, kSolver, kScalarInv, lMinP, lMaxP, evD, posD, ortho = None, qProxCtx = None},
    nt = ctrlBasis["NumSlices"];
    nParams = ctrlBasis["NumParams"];
    mCtrl = cSpec["NumControls"];
    bMat = ctrlBasis["BasisMatrix"];
    bT = ctrlBasis["BasisTranspose"];
    isIdent = ctrlBasis["IsIdentityQ"];
    hasBox = TrueQ[cSpec["HasBox"]];
    hasDisk = TrueQ[cSpec["HasDisk"]];
    hasSlew = TrueQ[cSpec["HasSlew"]] && nt > 1;
    hasL1 = TrueQ[lam1 > 0.];
    hasTV = TrueQ[lamTV > 0.] && nt > 1;

    mode = Which[
        !hasBox && !hasDisk && !hasSlew && !hasL1 && !hasTV, "None",
        isIdent && !hasDisk && !hasSlew && !hasTV, "FastBoxL1",
        isIdent && cSpec["DisjointDisks"] && !hasBox && !hasSlew && !hasTV, "FastDisk",
        isIdent && !hasBox && !hasDisk && !hasSlew && hasTV, "FastFusedTV",
        True, "ADMM"
    ];

    mergeL1IntoB = hasDisk && cSpec["DisjointDisks"] && !hasBox && hasL1;
    useA = hasBox || (hasL1 && !mergeL1IntoB);
    useB = hasDisk;
    useC = (hasSlew || hasTV) && nt > 1;

    d1Mat = None; dbMat = None; dbT = None; kSolver = None; kScalarInv = 1.;
    rhoA = 0.; rhoB = 0.; rhoC = 0.;
    If[mode === "ADMM",
        btb = bT . bMat;
        If[useA, rhoA = 2.0 / Max[0.05, Tr[btb] / Max[nParams, 1]]];
        If[useB, rhoB = 2.0 / Max[0.05, Tr[btb] / Max[nParams, 1]]];
        If[useC,
            d1Mat = SparseArray[{Band[{1, 1}] -> -1., Band[{1, 2}] -> 1.}, {nt - 1, nt}];
            dbMat = ToPackedRealMatrix[Normal[d1Mat . bMat]];
            dbT = Developer`ToPackedArray[Transpose[dbMat]];
            dbtd = dbT . dbMat;
            (* Ghadimi minimax-optimal ADMM penalty: (1 + rhoA + rhoB) / Sqrt[lam_min^+ * lam_max] *)
            If[isIdent,
                lMinP = 4. * Sin[Pi / (2. * nt)]^2;
                lMaxP = 4. * Cos[Pi / (2. * nt)]^2,
                evD = Re[Eigenvalues[dbtd]];
                lMaxP = Max[evD];
                posD = Select[evD, # > 10^-10 * Max[lMaxP, 10^-12] &];
                lMinP = If[posD === {}, lMaxP, Min[posD]]
            ];
            rhoC = If[lMaxP > 10^-12,
                Clip[(1. + rhoA + rhoB) / Max[Sqrt[lMinP * lMaxP], 10^-12], {2.0, 50.0}],
                2.0
            ]
        ];
        If[isIdent && !useC,
            kScalarInv = 1. / (1. + rhoA + rhoB),
            diagWeight = rhoA + rhoB;
            kMat = Developer`ToPackedArray[
                N[IdentityMatrix[nParams]] + diagWeight * btb + If[useC, rhoC * dbtd, 0.]
            ];
            kSolver = LinearSolve[kMat]
        ];
        If[TrueQ[allowOrtho] && !isIdent,
            ortho = buildOrthonormalBasisReduction[ctrlBasis];
            If[AssociationQ[ortho],
                qProxCtx = compileProxContext[ortho["QBasis"], dtVec, cSpec, lam1, lamTV, False],
                ortho = None
            ]
        ]
    ];

    <|
        "Mode" -> mode,
        "IsIdentityQ" -> isIdent,
        "BasisMatrix" -> bMat,
        "BasisTranspose" -> bT,
        "NumSlices" -> nt,
        "NumParams" -> nParams,
        "NumControls" -> mCtrl,
        "LowerBoundsMat" -> cSpec["LowerBoundsMat"],
        "UpperBoundsMat" -> cSpec["UpperBoundsMat"],
        "FiniteBox" -> cSpec["FiniteBox"],
        "HasBox" -> hasBox,
        "NormGroups" -> cSpec["NormGroups"],
        "NormLimits" -> cSpec["NormLimits"],
        "DisjointDisks" -> cSpec["DisjointDisks"],
        "HasDisk" -> hasDisk,
        "SlewDeltaMat" -> cSpec["SlewDeltaMat"],
        "HasSlew" -> hasSlew,
        "L1Penalty" -> N[Max[0., lam1]],
        "TVPenalty" -> N[Max[0., lamTV]],
        "HasNonSmoothPenalties" -> (hasL1 || hasTV),
        "HasAnyConstraint" -> (hasBox || hasDisk || hasSlew),
        "UseA" -> useA, "UseB" -> useB, "UseC" -> useC, "MergeL1B" -> mergeL1IntoB,
        "RhoA" -> rhoA, "RhoB" -> rhoB, "RhoC" -> rhoC,
        "DiffBasis" -> dbMat,
        "DiffBasisT" -> dbT,
        "KSolver" -> kSolver,
        "KScalarInv" -> kScalarInv,
        "OrthoReduction" -> ortho,
        "QProxContext" -> qProxCtx
    |>
];

applyProximal[ctx_Association, p0Mat_, eta_?NumericQ] := Module[
    {mode = ctx["Mode"], lam1 = ctx["L1Penalty"], lamTV = ctx["TVPenalty"],
     isIdent = ctx["IsIdentityQ"], bMat = ctx["BasisMatrix"], bT = ctx["BasisTranspose"],
     lbMat = ctx["LowerBoundsMat"], ubMat = ctx["UpperBoundsMat"], finiteBox = ctx["FiniteBox"],
     hasBox = ctx["HasBox"], groups = ctx["NormGroups"], lims = ctx["NormLimits"], disj = ctx["DisjointDisks"],
     hasSlew = ctx["HasSlew"], sDelta = ctx["SlewDeltaMat"], nt = ctx["NumSlices"], mCtrl = ctx["NumControls"],
     useA = ctx["UseA"], useB = ctx["UseB"], useC = ctx["UseC"], mergeL1B = Lookup[ctx, "MergeL1B", False],
     rhoA = ctx["RhoA"], rhoB = ctx["RhoB"], rhoC = ctx["RhoC"],
     dbMat = ctx["DiffBasis"], dbT = ctx["DiffBasisT"], kSolver = ctx["KSolver"], kInv = ctx["KScalarInv"],
     u0, outMat, p, bp, dbp, zA, zB, zC, wA, wB, wC, zAOld, zBOld, zCOld,
     bpR, bpRB, dbpR, rhs, rPri, rDual, tolPri, tolDual, alphaOR = 1.6, tauA, tauB, tauC},
    Switch[mode,
        "None",
            ToPackedRealMatrix[p0Mat],
        "FastBoxL1",
            outMat = If[lam1 > 0. && eta > 0., softThresholdMat[p0Mat, eta * lam1], p0Mat];
            If[hasBox, clampBoxMat[outMat, lbMat, ubMat, finiteBox], ToPackedRealMatrix[outMat]],
        "FastDisk",
            outMat = If[lam1 > 0. && eta > 0., softThresholdMat[p0Mat, eta * lam1], p0Mat];
            projDiskGroups[outMat, groups, lims, disj],
        "FastFusedTV",
            outMat = Transpose[Table[proxTV1D[p0Mat[[All, k]], eta * lamTV], {k, mCtrl}]];
            If[lam1 > 0. && eta > 0., softThresholdMat[outMat, eta * lam1], Developer`ToPackedArray[outMat]],
        "ADMM",
            u0 = If[isIdent, p0Mat, bMat . p0Mat];
            If[eta == 0. && computeConstraintViolation[ctx, u0] <= 10^-13,
                Return[ToPackedRealMatrix[p0Mat]]
            ];
            p = ToPackedRealMatrix[p0Mat];
            bp = If[isIdent, p, Developer`ToPackedArray[bMat . p]];
            dbp = If[useC, Developer`ToPackedArray[dbMat . p], None];
            tauA = If[rhoA > 0., eta * lam1 / rhoA, 0.];
            tauB = If[mergeL1B && rhoB > 0., eta * lam1 / rhoB, 0.];
            tauC = If[rhoC > 0., eta * lamTV / rhoC, 0.];
            If[useA,
                zA = If[hasBox, clampBoxMat[softThresholdMat[bp, tauA], lbMat, ubMat, finiteBox], softThresholdMat[bp, tauA]];
                wA = ConstantArray[0., {nt, mCtrl}],
                zA = 0.; wA = 0.
            ];
            If[useB,
                zB = projDiskGroups[If[tauB > 0., softThresholdMat[bp, tauB], bp], groups, lims, disj];
                wB = ConstantArray[0., {nt, mCtrl}],
                zB = 0.; wB = 0.
            ];
            If[useC,
                zC = softThresholdMat[dbp, tauC];
                If[hasSlew, zC = clampBoxMat[zC, -sDelta, sDelta, True]];
                wC = ConstantArray[0., {nt - 1, mCtrl}],
                zC = 0.; wC = 0.
            ];
            tolPri = 10^-10;
            tolDual = Max[10^-9 * Max[1., Norm[Flatten[p0Mat]]], 10^-10];
            Do[
                rhs = p0Mat +
                    If[useA || useB,
                        With[{comb = If[useA, rhoA * (zA - wA), 0.] + If[useB, rhoB * (zB - wB), 0.]},
                            If[isIdent, comb, bT . comb]
                        ],
                        0.
                    ] +
                    If[useC, dbT . (rhoC * (zC - wC)), 0.];
                p = Developer`ToPackedArray[If[kSolver =!= None, kSolver[rhs], kInv * rhs]];
                bp = If[isIdent, p, Developer`ToPackedArray[bMat . p]];
                If[useC, dbp = Developer`ToPackedArray[dbMat . p]];
                rPri = 0.; rDual = 0.;
                If[useA,
                    zAOld = zA;
                    bpR = alphaOR * bp + (1. - alphaOR) * zAOld;
                    zA = If[hasBox, clampBoxMat[softThresholdMat[bpR + wA, tauA], lbMat, ubMat, finiteBox], softThresholdMat[bpR + wA, tauA]];
                    wA += bpR - zA;
                    rPri = Max[rPri, Max[Abs[bp - zA]]];
                    rDual = Max[rDual, rhoA * Max[Abs[zA - zAOld]]]
                ];
                If[useB,
                    zBOld = zB;
                    bpRB = alphaOR * bp + (1. - alphaOR) * zBOld;
                    zB = projDiskGroups[If[tauB > 0., softThresholdMat[bpRB + wB, tauB], bpRB + wB], groups, lims, disj];
                    wB += bpRB - zB;
                    rPri = Max[rPri, Max[Abs[bp - zB]]];
                    rDual = Max[rDual, rhoB * Max[Abs[zB - zBOld]]]
                ];
                If[useC,
                    zCOld = zC;
                    dbpR = alphaOR * dbp + (1. - alphaOR) * zCOld;
                    zC = softThresholdMat[dbpR + wC, tauC];
                    If[hasSlew, zC = clampBoxMat[zC, -sDelta, sDelta, True]];
                    wC += dbpR - zC;
                    rPri = Max[rPri, Max[Abs[dbp - zC]]];
                    rDual = Max[rDual, rhoC * Max[Abs[zC - zCOld]]]
                ];
                If[rPri <= tolPri && rDual <= tolDual, Break[]],
                {2500}
            ];
            ToPackedRealMatrix[p]
    ]
];

(* Bounded-simplex projection & slice-duration optimization specification helpers *)
projectBoundedSimplex[yVec_List, lbVec_List, ubVec_List, sTarget_?NumericQ] := Module[
    {n = Length[yVec], y = N[yVec], lb = N[lbVec], ub = N[ubVec], sT = N[sTarget],
     sLo, sHi, nuLo, nuHi, nuMid, x, sMid, free, tol},
    If[n === 0, Return[Developer`ToPackedArray[{}]]];
    If[n === 1, Return[Developer`ToPackedArray[{Clip[sT, {First[lb], First[ub]}]}]]];
    sLo = Total[lb];
    sHi = Total[ub];
    If[sT <= sLo + 10^-15 * Max[1., Abs[sLo]], Return[Developer`ToPackedArray[lb]]];
    If[sT >= sHi - 10^-15 * Max[1., Abs[sHi]], Return[Developer`ToPackedArray[ub]]];
    nuLo = Min[y - ub] - 1.;
    nuHi = Max[y - lb] + 1.;
    tol = 10^-15 * Max[1., Abs[sT]];
    Do[
        nuMid = 0.5 * (nuLo + nuHi);
        x = MapThread[Clip[#1 - nuMid, {#2, #3}] &, {y, lb, ub}];
        sMid = Total[x];
        If[Abs[sMid - sT] <= tol || Abs[nuHi - nuLo] <= 10^-15 * Max[1., Abs[nuLo], Abs[nuHi]],
            Break[]
        ];
        If[sMid > sT, nuLo = nuMid, nuHi = nuMid],
        {90}
    ];
    x = MapThread[Clip[#1 - 0.5 * (nuLo + nuHi), {#2, #3}] &, {y, lb, ub}];
    free = Pick[Range[n], MapThread[#1 < #2 && #2 < #3 &, {lb, x, ub}]];
    If[Length[free] >= 1,
        x[[free]] = MapThread[Clip[#1, {#2, #3}] &, {
            x[[free]] + (sT - Total[x]) / Length[free], lb[[free]], ub[[free]]
        }]
    ];
    Developer`ToPackedArray[x]
];

parseDurationPenaltySpec[opt_, total0_?NumericQ] := Module[{tDef = Max[N[total0], 10^-12]},
    Which[
        opt === 0 || opt === 0. || opt === None || opt === False,
            <|"Lambda" -> 0., "Reference" -> tDef|>,
        NumericQ[opt] && TrueQ[Im[N[opt]] == 0 && Re[opt] >= 0.],
            <|"Lambda" -> N[Re[opt]], "Reference" -> tDef|>,
        MatchQ[opt, {_?NumericQ, _?NumericQ}] && TrueQ[Im[N[opt[[1]]]] == 0 && Im[N[opt[[2]]]] == 0 && Re[opt[[1]]] >= 0. && Re[opt[[2]]] > 0.],
            <|"Lambda" -> N[Re[opt[[1]]]], "Reference" -> N[Re[opt[[2]]]]|>,
        True,
            Message[OptimalControlProblem::baddur, "DurationPenalty", opt, tDef];
            $Failed
    ]
];

resolveDurationSpec[optDurOpt_, durBndOpt_, fixedTotOpt_, dtVec_List, initOpt_,
                    ctrlBasis_Association, cSpec_Association, driftKind_String, mCtrl_Integer] := Module[
    {nt = Length[dtVec], isIdent = ctrlBasis["IsIdentityQ"], activeIdx, freezeIdx = {},
     uGuess = None, nAct, dtAct0, dMin, dMax, hasFixed, totTar, actTar,
     dtScaleScalar, cSpecOut, lbMat, ubMat, lbVec, ubVec, minSq, gIdx, wLim, okFeas = True,
     total0 = Total[dtVec], durMode = "Slices", wVec = {}},
    If[optDurOpt === False || optDurOpt === None,
        Return[<|
            "HasDurationOptimization" -> False,
            "DurationMode" -> "None",
            "DurationWeights" -> Developer`ToPackedArray[N[dtVec / Max[total0, 10^-12]]],
            "ActiveIndices" -> {},
            "FreezeControlSlices" -> {},
            "NumActive" -> 0,
            "DurationLowerBounds" -> Developer`ToPackedArray[{}],
            "DurationUpperBounds" -> Developer`ToPackedArray[{}],
            "FixedTotalDuration" -> False,
            "FixedTotalTarget" -> None,
            "ActiveTargetSum" -> None,
            "DurationScale" -> 1.,
            "ModifiedConstraintSpec" -> cSpec
        |>]
    ];
    If[driftKind === "SliceList",
        Message[OptimalControlProblem::nodrifttime];
        Return[$Failed]
    ];
    If[TrueQ[cSpec["HasSlew"]],
        Message[OptimalControlProblem::baddur, optDurOpt, durBndOpt, "SlewRateBound cannot be combined with OptimizeDurations"];
        Return[$Failed]
    ];
    Which[
        optDurOpt === "Total",
            If[fixedTotOpt =!= False && fixedTotOpt =!= None,
                Message[OptimalControlProblem::baddur, optDurOpt, durBndOpt, fixedTotOpt];
                Return[$Failed]
            ];
            durMode = "Total";
            activeIdx = Range[nt];
            wVec = Developer`ToPackedArray[N[dtVec / Max[total0, 10^-12]]];
            Which[
                durBndOpt === Automatic,
                    dMin = {Max[10^-12, 0.01 * total0]};
                    dMax = {Max[10. * total0, 10^-6]},
                MatchQ[durBndOpt, {_?NumericQ, _?NumericQ}] && TrueQ[0. < durBndOpt[[1]] <= durBndOpt[[2]]],
                    dMin = {N[durBndOpt[[1]]]};
                    dMax = {N[durBndOpt[[2]]]},
                MatchQ[durBndOpt, {{_?NumericQ, _?NumericQ}}] && TrueQ[0. < durBndOpt[[1, 1]] <= durBndOpt[[1, 2]]],
                    dMin = {N[durBndOpt[[1, 1]]]};
                    dMax = {N[durBndOpt[[1, 2]]]},
                True,
                    Message[OptimalControlProblem::baddur, optDurOpt, durBndOpt, fixedTotOpt];
                    Return[$Failed]
            ];
            dtScaleScalar = 1. / Max[total0, 10^-12];
            Return[<|
                "HasDurationOptimization" -> True,
                "DurationMode" -> "Total",
                "DurationWeights" -> wVec,
                "ActiveIndices" -> activeIdx,
                "FreezeControlSlices" -> {},
                "NumActive" -> 1,
                "DurationLowerBounds" -> Developer`ToPackedArray[dMin],
                "DurationUpperBounds" -> Developer`ToPackedArray[dMax],
                "FixedTotalDuration" -> False,
                "FixedTotalTarget" -> None,
                "ActiveTargetSum" -> None,
                "DurationScale" -> dtScaleScalar,
                "ModifiedConstraintSpec" -> cSpec
            |>],
        optDurOpt === True || optDurOpt === "All",
            activeIdx = Range[nt],
        ListQ[optDurOpt] && optDurOpt =!= {} &&
                AllTrue[optDurOpt, IntegerQ[#] && 1 <= # <= nt &] &&
                Length[DeleteDuplicates[optDurOpt]] === Length[optDurOpt],
            activeIdx = optDurOpt,
        MatchQ[optDurOpt, {"DelaysOnly", _List}] && isIdent && optDurOpt[[2]] =!= {} &&
                AllTrue[optDurOpt[[2]], IntegerQ[#] && 1 <= # <= nt &] &&
                Length[DeleteDuplicates[optDurOpt[[2]]]] === Length[optDurOpt[[2]]],
            activeIdx = optDurOpt[[2]];
            freezeIdx = activeIdx,
        optDurOpt === "DelaysOnly" && isIdent,
            uGuess = Which[
                MatrixQ[initOpt, NumericQ] && Dimensions[initOpt] === {nt, mCtrl}, N[Re[initOpt]],
                MatchQ[initOpt, _OptimalControlResult] && MatrixQ[initOpt["Controls"], NumericQ] &&
                    Dimensions[initOpt["Controls"]] === {nt, mCtrl}, N[Re[initOpt["Controls"]]],
                VectorQ[initOpt, NumericQ] && Length[initOpt] === nt * mCtrl, Partition[N[Re[initOpt]], mCtrl],
                True, None
            ];
            If[!MatrixQ[uGuess, NumericQ],
                Message[OptimalControlProblem::baddur, optDurOpt, durBndOpt, fixedTotOpt];
                Return[$Failed]
            ];
            freezeIdx = Select[Range[nt], Max[Abs[uGuess[[#]]]] <= 10^-12 &];
            If[freezeIdx === {},
                Message[OptimalControlProblem::baddur, optDurOpt, durBndOpt, fixedTotOpt];
                Return[$Failed]
            ];
            activeIdx = freezeIdx,
        True,
            Message[OptimalControlProblem::baddur, optDurOpt, durBndOpt, fixedTotOpt];
            Return[$Failed]
    ];

    cSpecOut = cSpec;
    If[freezeIdx =!= {},
        lbMat = N[cSpec["LowerBoundsMat"]];
        ubMat = N[cSpec["UpperBoundsMat"]];
        lbMat[[freezeIdx, All]] = 0.;
        ubMat[[freezeIdx, All]] = 0.;
        lbVec = Developer`ToPackedArray[Flatten[lbMat]];
        ubVec = Developer`ToPackedArray[Flatten[ubMat]];
        If[TrueQ[cSpec["HasDisk"]],
            Do[
                gIdx = (cSpec["NormGroups"])[[g]];
                wLim = (cSpec["NormLimits"])[[g]];
                minSq = Total[(MapThread[Clip[0., {#1, #2}] &, {lbMat[[First[freezeIdx], gIdx]], ubMat[[First[freezeIdx], gIdx]]}])^2];
                If[minSq > (wLim + 10^-9 * Max[1., wLim])^2, okFeas = False; Break[]],
                {g, Length[cSpec["NormGroups"]]}
            ];
            If[!okFeas,
                Message[OptimalControlProblem::baddur, optDurOpt, durBndOpt, fixedTotOpt];
                Return[$Failed]
            ]
        ];
        cSpecOut = Join[cSpec, <|
            "LowerBounds" -> lbVec,
            "UpperBounds" -> ubVec,
            "LowerBoundsMat" -> lbMat,
            "UpperBoundsMat" -> ubMat,
            "FiniteBox" -> AllTrue[Flatten[lbMat], NumericQ[#] && # > -Infinity &] &&
                AllTrue[Flatten[ubMat], NumericQ[#] && # < Infinity &],
            "HasBox" -> True,
            "HasAnyConstraint" -> True
        |>]
    ];

    nAct = Length[activeIdx];
    dtAct0 = dtVec[[activeIdx]];
    Which[
        durBndOpt === Automatic,
            dMin = ConstantArray[Max[10^-13, Min[10^-7 * Max[total0, 1.], 0.01 * Min[dtAct0]]], nAct];
            dMax = ConstantArray[Max[10. * Max[dtAct0], total0, 10^-6], nAct],
        MatchQ[durBndOpt, {_?NumericQ, _?NumericQ}] && TrueQ[0. < durBndOpt[[1]] <= durBndOpt[[2]]],
            dMin = ConstantArray[N[durBndOpt[[1]]], nAct];
            dMax = ConstantArray[N[durBndOpt[[2]]], nAct],
        ListQ[durBndOpt] && Length[durBndOpt] === nAct &&
                AllTrue[durBndOpt, MatchQ[#, {_?NumericQ, _?NumericQ}] && TrueQ[0. < #[[1]] <= #[[2]]] &],
            dMin = N[durBndOpt[[All, 1]]];
            dMax = N[durBndOpt[[All, 2]]],
        True,
            Message[OptimalControlProblem::baddur, optDurOpt, durBndOpt, fixedTotOpt];
            Return[$Failed]
    ];

    Which[
        fixedTotOpt === False || fixedTotOpt === None,
            hasFixed = False; totTar = None; actTar = None,
        fixedTotOpt === True || fixedTotOpt === Automatic,
            hasFixed = True; totTar = N[total0];
            actTar = totTar - Total[dtVec[[Complement[Range[nt], activeIdx]]]],
        NumericQ[fixedTotOpt] && TrueQ[Re[fixedTotOpt] > 0] && Im[N[fixedTotOpt]] == 0,
            hasFixed = True; totTar = N[Re[fixedTotOpt]];
            actTar = totTar - Total[dtVec[[Complement[Range[nt], activeIdx]]]],
        True,
            Message[OptimalControlProblem::baddur, optDurOpt, durBndOpt, fixedTotOpt];
            Return[$Failed]
    ];
    If[hasFixed && (!TrueQ[actTar > 0.] || Total[dMin] > actTar + 10^-12 * Max[1., actTar] ||
                    Total[dMax] < actTar - 10^-12 * Max[1., actTar]),
        Message[OptimalControlProblem::baddur, optDurOpt, durBndOpt, fixedTotOpt];
        Return[$Failed]
    ];

    dtScaleScalar = 1. / Max[Mean[dtAct0], 10^-12];
    <|
        "HasDurationOptimization" -> True,
        "DurationMode" -> durMode,
        "DurationWeights" -> Developer`ToPackedArray[N[dtVec / Max[total0, 10^-12]]],
        "ActiveIndices" -> activeIdx,
        "FreezeControlSlices" -> freezeIdx,
        "NumActive" -> nAct,
        "DurationLowerBounds" -> Developer`ToPackedArray[dMin],
        "DurationUpperBounds" -> Developer`ToPackedArray[dMax],
        "FixedTotalDuration" -> hasFixed,
        "FixedTotalTarget" -> totTar,
        "ActiveTargetSum" -> actTar,
        "DurationScale" -> dtScaleScalar,
        "ModifiedConstraintSpec" -> cSpecOut
    |>
];

buildTimingDerivTensors[comp_Association, dtVecCur_List, driftSched_] := Module[
    {nt = Length[dtVecCur], isTDCont = (Lookup[comp, "DriftKind", "Static"] === "Continuous"),
     t0 = Lookup[comp, "InitialTimePoint", 0.], intSpec, dFn, stageMeta, nFine,
     dSelf, dShift, s, meta, nodes, wts, xis, d1, d2, alphas},
    If[!isTDCont,
        Return[<|
            "DtCoarse" -> Developer`ToPackedArray[N[dtVecCur]],
            "InitialTimePoint" -> N[t0],
            "Alphas" -> Developer`ToPackedArray[ConstantArray[1., nt]],
            "IsTDContinuous" -> False,
            "DSelf" -> None,
            "DShift" -> None
        |>]
    ];
    stageMeta = driftSched["StageMeta"];
    nFine = Length[stageMeta];
    intSpec = comp["DriftIntegrationSpec"];
    dFn = comp["DriftDerivEvaluator"];
    dSelf = Table[Null, {nFine}];
    dShift = Table[Null, {nFine}];
    If[intSpec["Method"] === "Midpoint",
        Do[
            meta = stageMeta[[s]];
            d1 = dFn[(meta["Nodes"])[[1]]];
            dShift[[s]] = d1;
            dSelf[[s]] = ToPackedComplexMatrix[((meta["ShiftPrefixes"])[[1]]) * d1],
            {s, nFine}
        ],
        Do[
            meta = stageMeta[[s]];
            {d1, d2} = dFn /@ meta["Nodes"];
            xis = meta["ShiftPrefixes"];
            wts = meta["Weights"];
            dShift[[s]] = ToPackedComplexMatrix[wts[[1]] * d1 + wts[[2]] * d2];
            dSelf[[s]] = ToPackedComplexMatrix[wts[[1]] * xis[[1]] * d1 + wts[[2]] * xis[[2]] * d2];
            wts = stageMeta[[s + 1, "Weights"]];
            dShift[[s + 1]] = ToPackedComplexMatrix[wts[[1]] * d1 + wts[[2]] * d2];
            dSelf[[s + 1]] = ToPackedComplexMatrix[wts[[1]] * xis[[1]] * d1 + wts[[2]] * xis[[2]] * d2],
            {s, 1, nFine, 2}
        ]
    ];
    alphas = Developer`ToPackedArray[#["Alpha"] & /@ stageMeta];
    <|
        "DtCoarse" -> Developer`ToPackedArray[N[dtVecCur]],
        "InitialTimePoint" -> N[t0],
        "Alphas" -> alphas,
        "IsTDContinuous" -> True,
        "DSelf" -> Developer`ToPackedArray[dSelf],
        "DShift" -> Developer`ToPackedArray[dShift]
    |>
];

(* Core compilation of an OptimalControlProblem *)
CompileOptimalControlProblem[prob:OptimalControlProblem[drift_, controls_, target_, timing_, opts___Rule]] := Module[
    {basisOpt, opbasisOpt, spaceOpt, methodOpt, paramOpt, winOpt, freqUnit, ctrlUnit,
     defaultNt, t0Opt, driftIntOpt, runCostOpt, optDurOpt, durBndOpt, durPenOpt, fixedTotOpt, initOpt,
     powerPen, smoothPen, l1Pen, tvPen, ctrlBounds, ampBound, ampNormBound, slewBound,
     basis, opbasis, ctrlList, hasSop, space, dim, dtVec, totalDur, nt, hamFactor, ctrlFactor, t0,
     driftComp, h0Mat, hCtrlsRaw, hCtrlMats, hermErr, mCtrl, compiledObj, runningCostSpec, svAllowed,
     ctrlBasis, nParams, method, hCtrlRows, hFlatT, chanNorms, colScales, ctrlScales, cSpecRaw,
     durSpec, dpSpec, cSpec, proxCtx, hasDur, pfSpec, hasPhi, nPhi, baseScales, lbOut, ubOut, hasBoundsOut},
    {basisOpt, opbasisOpt, spaceOpt, methodOpt, paramOpt, winOpt, freqUnit, ctrlUnit,
     defaultNt, t0Opt, driftIntOpt, runCostOpt, optDurOpt, durBndOpt, durPenOpt, fixedTotOpt, initOpt,
     powerPen, smoothPen, l1Pen, tvPen, ctrlBounds, ampBound, ampNormBound, slewBound} =
        {Basis, OperatorBasis, PropagateSpace, Method, ControlParameterization, ControlWindow,
         FrequencyUnits, ControlUnits, TimeSteps, InitialTimePoint, DriftIntegration,
         RunningCost, OptimizeDurations, DurationBounds, DurationPenalty, FixedTotalDuration, InitialControls,
         PowerPenalty, SmoothnessPenalty, L1Penalty, TVPenalty,
         ControlBounds, AmplitudeBound, AmplitudeNormBound, SlewRateBound} /. {opts} /. Options[OptimalControlProblem];

    basis = If[basisOpt === Automatic, ZeemanBasis[], basisOpt];
    opbasis = If[opbasisOpt === Automatic, ZeemanKetBraOperatorBasis[Sorted -> "Canonical"], opbasisOpt];
    ctrlList = If[ListQ[controls] && !MatrixQ[controls], controls, {controls}];
    If[ctrlList === {}, Message[OptimalControlProblem::badctrl]; Return[$Failed]];

    dtVec = resolveTimingVector[timing, defaultNt];
    If[dtVec === $Failed, Return[$Failed]];
    nt = Length[dtVec];
    totalDur = Total[dtVec];

    If[!NumericQ[t0Opt] || !FreeQ[t0Opt, Complex | DirectedInfinity | Indeterminate] || Im[N[t0Opt]] != 0,
        Message[OptimalControlProblem::badt0, t0Opt]; Return[$Failed]
    ];
    t0 = N[Re[t0Opt]];

    hamFactor = resolveUnitFactor[freqUnit];
    ctrlFactor = If[ctrlUnit === Automatic, hamFactor, resolveUnitFactor[ctrlUnit]];
    If[hamFactor === $Failed || ctrlFactor === $Failed, Return[$Failed]];

    hasSop = driftContainsSuperoperatorQ[drift, t0, totalDur] || AnyTrue[ctrlList, SuperoperatorQ];
    space = If[spaceOpt === Automatic,
        If[hasSop, "Liouville", "Hilbert"],
        If[MemberQ[{"Hilbert", "Liouville"}, spaceOpt], spaceOpt, $Failed]
    ];
    If[space === $Failed, Return[$Failed]];
    If[space === "Hilbert" && hasSop,
        Message[OptimalControlProblem::nonherm, "superoperator generator"]; Return[$Failed]
    ];

    dim = If[space === "Hilbert", BasisDimension[basis], OperatorBasisDimension[opbasis]];
    If[!IntegerQ[dim] || dim < 1, Return[$Failed]];

    driftComp = compileDriftRepresentation[drift, space, basis, opbasis, dim, hamFactor, dtVec, t0, driftIntOpt];
    If[driftComp === $Failed, Return[$Failed]];
    h0Mat = driftComp["Drift"];

    hCtrlMats = toGenMatrix[#, space, basis, opbasis, dim] & /@ ctrlList;
    If[!FreeQ[hCtrlMats, $Failed],
        Message[OptimalControlProblem::badctrl]; Return[$Failed]
    ];
    hCtrlsRaw = hCtrlMats;
    hCtrlMats = ToPackedComplexMatrix[ctrlFactor * #] & /@ hCtrlMats;
    mCtrl = Length[hCtrlMats];

    If[space === "Hilbert",
        hermErr = Max[Map[Max[Abs[# - ConjugateTranspose[#]]] / Max[1., Max[Abs[#]]] &, hCtrlMats]];
        If[hermErr > 10^-8,
            Message[OptimalControlProblem::nonherm, hermErr]; Return[$Failed]
        ];
        hCtrlMats = HermitianSymmetrize /@ hCtrlMats
    ];

    ctrlBasis = BuildControlBasis[paramOpt, dtVec, winOpt];
    If[ctrlBasis === $Failed,
        Message[OptimalControlProblem::badbasis, paramOpt, winOpt]; Return[$Failed]
    ];
    nParams = ctrlBasis["NumParams"];

    compiledObj = CompileObjective[target, space, basis, opbasis];
    If[compiledObj === $Failed, Return[$Failed]];

    runningCostSpec = CompileRunningCost[runCostOpt, compiledObj, space, basis, opbasis, dtVec, t0];
    If[runningCostSpec === $Failed, Return[$Failed]];
    If[AssociationQ[runningCostSpec],
        runningCostSpec = Join[runningCostSpec, <|"InitialTimePoint" -> t0|>]
    ];

    svAllowed = compiledObj["PureKetFastPath"] =!= None &&
        (runningCostSpec === None || TrueQ[runningCostSpec["CompatibleWithStateVectorQ"]]);
    method = Which[
        methodOpt =!= Automatic, methodOpt,
        space === "Liouville", "Liouville",
        svAllowed && dim >= 16, "StateVector",
        True, "Spectral"
    ];
    If[(space === "Hilbert" && !MemberQ[{"Spectral", "BlockMatrixExp", "FirstOrder", "StateVector"}, method]) ||
       (space === "Liouville" && method =!= "Liouville") ||
       (method === "StateVector" && !svAllowed),
        Message[OptimalControlProblem::badmethod, method, space]; Return[$Failed]
    ];

    cSpecRaw = resolveConstraintSpec[ctrlBounds, ampBound, ampNormBound, slewBound, dtVec, nParams, mCtrl];
    If[cSpecRaw === $Failed, Return[$Failed]];
    durSpec = resolveDurationSpec[optDurOpt, durBndOpt, fixedTotOpt, dtVec, initOpt,
        ctrlBasis, cSpecRaw, driftComp["DriftKind"], mCtrl];
    If[durSpec === $Failed, Return[$Failed]];
    dpSpec = parseDurationPenaltySpec[durPenOpt, totalDur];
    If[dpSpec === $Failed, Return[$Failed]];

    cSpec = durSpec["ModifiedConstraintSpec"];
    hasDur = TrueQ[durSpec["HasDurationOptimization"]];
    pfSpec = Lookup[compiledObj, "PhaseFreedomSpec", None];
    hasPhi = AssociationQ[pfSpec] && TrueQ[pfSpec["ParametricQ"]] && TrueQ[pfSpec["NumPhases"] > 0];
    nPhi = If[hasPhi, pfSpec["NumPhases"], 0];
    proxCtx = compileProxContext[ctrlBasis, dtVec, cSpec, N[Max[0., l1Pen]], N[Max[0., tvPen]]];

    chanNorms = Map[Max[Norm[Flatten[#]] / Sqrt[N[dim]], 10^-6] &, hCtrlMats];
    colScales = Map[Max[Norm[dtVec * #], 10^-12] &, Transpose[ctrlBasis["BasisMatrix"]]];
    ctrlScales = Developer`ToPackedArray[Flatten[Outer[Times, colScales, chanNorms]]];
    baseScales = Developer`ToPackedArray[Join[
        ctrlScales,
        If[hasDur, ConstantArray[durSpec["DurationScale"], durSpec["NumActive"]], {}],
        If[hasPhi, ConstantArray[1., nPhi], {}]
    ]];
    lbOut = Developer`ToPackedArray[Join[
        cSpec["LowerBounds"],
        If[hasDur, durSpec["DurationLowerBounds"], {}],
        If[hasPhi, ConstantArray[-Infinity, nPhi], {}]
    ]];
    ubOut = Developer`ToPackedArray[Join[
        cSpec["UpperBounds"],
        If[hasDur, durSpec["DurationUpperBounds"], {}],
        If[hasPhi, ConstantArray[Infinity, nPhi], {}]
    ]];
    hasBoundsOut = TrueQ[cSpec["HasBox"] || hasDur];

    hCtrlRows = ToPackedComplexMatrix[Flatten[#] & /@ hCtrlMats];
    hFlatT = ToPackedComplexMatrix[Flatten[Transpose[#]] & /@ hCtrlMats];

    <|
        "Space" -> space,
        "Method" -> method,
        "Dimension" -> dim,
        "Basis" -> basis,
        "OperatorBasis" -> opbasis,
        "Drift" -> h0Mat,
        "FineTimeSteps" -> driftComp["FineTimeSteps"],
        "NumStages" -> driftComp["NumStages"],
        "IsTimeDependentDrift" -> driftComp["IsTimeDependentDrift"],
        "DriftKind" -> driftComp["DriftKind"],
        "DriftIntegration" -> driftComp["DriftIntegration"],
        "DriftIntegrationSpec" -> Lookup[driftComp, "DriftIntegrationSpec", None],
        "InitialTimePoint" -> driftComp["InitialTimePoint"],
        "DriftEvaluator" -> driftComp["DriftEvaluator"],
        "DriftDerivEvaluator" -> driftComp["DriftDerivEvaluator"],
        "DriftStageMeta" -> driftComp["StageMeta"],
        "Controls" -> hCtrlMats,
        "ControlRows" -> hCtrlRows,
        "ControlFlatT" -> hFlatT,
        "SymbolicDrift" -> drift,
        "SymbolicControls" -> ctrlList,
        "DriftUnitFactor" -> hamFactor,
        "ControlUnitFactor" -> ctrlFactor,
        "FrequencyUnits" -> freqUnit,
        "ControlUnits" -> If[ctrlUnit === Automatic, freqUnit, ctrlUnit],
        "TimeSteps" -> dtVec,
        "TotalDuration" -> totalDur,
        "NumSlices" -> nt,
        "NumControls" -> mCtrl,
        "ControlBasis" -> ctrlBasis,
        "NumParams" -> nParams,
        "Objective" -> compiledObj,
        "PhaseFreedomSpec" -> pfSpec,
        "HasPhaseFreedom" -> hasPhi,
        "NumPhases" -> nPhi,
        "RunningCostSpec" -> runningCostSpec,
        "RawTarget" -> target,
        "PowerPenalty" -> N[Max[0., powerPen]],
        "SmoothnessPenalty" -> N[Max[0., smoothPen]],
        "L1Penalty" -> N[Max[0., l1Pen]],
        "TVPenalty" -> N[Max[0., tvPen]],
        "ControlParameterScales" -> ctrlScales,
        "ParameterScales" -> baseScales,
        "ChannelNorms" -> Developer`ToPackedArray[chanNorms],
        "ControlLowerBounds" -> cSpec["LowerBounds"],
        "ControlUpperBounds" -> cSpec["UpperBounds"],
        "LowerBounds" -> lbOut,
        "UpperBounds" -> ubOut,
        "HasControlBounds" -> cSpec["HasBox"],
        "HasBounds" -> hasBoundsOut,
        "HasDiskBounds" -> cSpec["HasDisk"],
        "HasSlewBounds" -> cSpec["HasSlew"],
        "HasConstraints" -> TrueQ[cSpec["HasAnyConstraint"] || hasDur],
        "ConstraintSpec" -> cSpec,
        "ProxContext" -> proxCtx,
        "HasDurationOptimization" -> hasDur,
        "DurationMode" -> durSpec["DurationMode"],
        "DurationWeights" -> durSpec["DurationWeights"],
        "OptimizeDurations" -> optDurOpt,
        "ActiveDurationIndices" -> durSpec["ActiveIndices"],
        "FreezeControlSlices" -> durSpec["FreezeControlSlices"],
        "NumActiveDurations" -> durSpec["NumActive"],
        "DurationLowerBounds" -> durSpec["DurationLowerBounds"],
        "DurationUpperBounds" -> durSpec["DurationUpperBounds"],
        "FixedTotalDuration" -> durSpec["FixedTotalDuration"],
        "FixedTotalTarget" -> durSpec["FixedTotalTarget"],
        "ActiveDurationTargetSum" -> durSpec["ActiveTargetSum"],
        "DurationScale" -> durSpec["DurationScale"],
        "DurationPenaltyLambda" -> dpSpec["Lambda"],
        "DurationPenaltyReference" -> dpSpec["Reference"],
        "Options" -> {opts},
        "OriginProblem" -> prob
    |>
];

(* Convert parameter matrix / slice matrix / joint (controls, durations, phases) input *)
unpackControlMatrix[compiled_Association, input_] := Module[
    {nt, mCtrl, nParams, basisAssoc, isIdent, pMat, uMat, hasDur, durMode, nAct, activeIdx,
     hasPhi, nPhi, phiVec, nCtrl, dtVec, jointFlag = False, tTot},
    nt = compiled["NumSlices"];
    mCtrl = compiled["NumControls"];
    nParams = compiled["NumParams"];
    basisAssoc = compiled["ControlBasis"];
    isIdent = basisAssoc["IsIdentityQ"];
    hasDur = TrueQ[Lookup[compiled, "HasDurationOptimization", False]];
    durMode = Lookup[compiled, "DurationMode", If[hasDur, "Slices", "None"]];
    nAct = If[hasDur, Lookup[compiled, "NumActiveDurations", 0], 0];
    activeIdx = Lookup[compiled, "ActiveDurationIndices", {}];
    hasPhi = TrueQ[Lookup[compiled, "HasPhaseFreedom", False]];
    nPhi = If[hasPhi, Lookup[compiled, "NumPhases", 0], 0];
    nCtrl = nParams * mCtrl;
    dtVec = Developer`ToPackedArray[N[compiled["TimeSteps"]]];
    phiVec = Developer`ToPackedArray[ConstantArray[0., nPhi]];

    If[(hasDur || hasPhi) && VectorQ[input, NumericQ] && Length[input] === nCtrl + nAct + nPhi,
        jointFlag = True;
        pMat = Partition[Developer`ToPackedArray[N[Re[input[[1 ;; nCtrl]]]]], mCtrl];
        If[hasDur,
            If[durMode === "Total",
                tTot = N[Re[input[[nCtrl + 1]]]];
                If[!TrueQ[tTot > 0.], Return[$Failed]];
                dtVec = Developer`ToPackedArray[tTot * compiled["DurationWeights"]],
                dtVec[[activeIdx]] = Developer`ToPackedArray[N[Re[input[[nCtrl + 1 ;; nCtrl + nAct]]]]]
            ]
        ];
        If[hasPhi, phiVec = Developer`ToPackedArray[N[Re[input[[nCtrl + nAct + 1 ;; -1]]]]]];
        uMat = If[isIdent, pMat, Developer`ToPackedArray[basisAssoc["BasisMatrix"] . pMat]];
        Return[{pMat, uMat, dtVec, phiVec, True}]
    ];
    If[(hasDur || hasPhi) && AssociationQ[input] &&
       (KeyExistsQ[input, "TimeSteps"] || KeyExistsQ[input, "TotalDuration"] || KeyExistsQ[input, "Phases"] || KeyExistsQ[input, "OptimalPhases"]),
        Module[{cIn = Lookup[input, "Parameters", Lookup[input, "Controls", $Failed]],
                dtIn = Lookup[input, "TimeSteps", None],
                tIn = Lookup[input, "TotalDuration", None],
                phIn = Lookup[input, "Phases", If[AssociationQ[Lookup[input, "OptimalPhases", None]], input[["OptimalPhases", "Vector"]], phiVec]],
                parsedCtrl},
            If[cIn === $Failed, Return[$Failed]];
            If[hasDur,
                Which[
                    VectorQ[dtIn, NumericQ] && Length[dtIn] === nt,
                        dtVec = Developer`ToPackedArray[N[Re[dtIn]]],
                    durMode === "Total" && NumericQ[tIn] && TrueQ[Re[tIn] > 0.],
                        dtVec = Developer`ToPackedArray[N[Re[tIn]] * compiled["DurationWeights"]],
                    dtIn =!= None || tIn =!= None,
                        Return[$Failed]
                ]
            ];
            If[hasPhi,
                If[!VectorQ[phIn, NumericQ] || Length[phIn] =!= nPhi, Return[$Failed]];
                phiVec = Developer`ToPackedArray[N[Re[phIn]]]
            ];
            parsedCtrl = unpackControlMatrix[KeyDrop[compiled, {"HasDurationOptimization", "HasPhaseFreedom"}], cIn];
            If[parsedCtrl === $Failed, Return[$Failed]];
            Return[{parsedCtrl[[1]], parsedCtrl[[2]], dtVec, phiVec, True}]
        ]
    ];

    pMat = Which[
        VectorQ[input, NumericQ] && Length[input] === nCtrl,
            Partition[Developer`ToPackedArray[N[Re[input]]], mCtrl],
        MatrixQ[input, NumericQ] && Dimensions[input] === {nParams, mCtrl},
            ToPackedRealMatrix[input],
        !isIdent && MatrixQ[input, NumericQ] && Dimensions[input] === {nt, mCtrl},
            LeastSquares[basisAssoc["BasisMatrix"], ToPackedRealMatrix[input]],
        !isIdent && VectorQ[input, NumericQ] && Length[input] === nt * mCtrl,
            LeastSquares[basisAssoc["BasisMatrix"], Partition[Developer`ToPackedArray[N[Re[input]]], mCtrl]],
        True,
            $Failed
    ];
    If[pMat === $Failed, Return[$Failed]];
    uMat = If[isIdent, pMat, Developer`ToPackedArray[basisAssoc["BasisMatrix"] . pMat]];
    {pMat, uMat, dtVec, phiVec, jointFlag}
];

(* Core objective + analytical gradient evaluation *)
EvaluateOptimalControl[compiled_Association, input_, returnProps_: False] := Module[
    {parsed, pMat, uMat, dtVec, phiVec, jointFlag, basisAssoc, space, method, nt, mCtrl,
     hasDur, durMode, activeIdx, freezeIdx, hasPhi, pfSpec, objCotFn, nPhi, driftCur, dtFine, nStages,
     driftSched, timingSpec, rcSpec, kernRes, fid, runVal, gradU, gradUFid, gDtKern, gPhi = None,
     lamP, lamS, lamT, tRef, penSmooth, penNS, penTotal, gRegU, gDtPen, diffDt, du, duSq,
     gU, gParam, gParamFid, gRegParam, gDt, gDurActive, gDurFidActive, gRegDurActive, gTotDur, gradOut},
    If[KeyExistsQ[$OCSolverRegistry, "IsEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[compiled]],
        Return[($OCSolverRegistry["EvaluateEnsemble"])[compiled, input]]
    ];
    parsed = unpackControlMatrix[compiled, input];
    If[parsed === $Failed, Return[$Failed]];
    {pMat, uMat, dtVec, phiVec, jointFlag} = parsed;
    basisAssoc = compiled["ControlBasis"];
    space = compiled["Space"];
    method = compiled["Method"];
    nt = compiled["NumSlices"];
    mCtrl = compiled["NumControls"];
    hasDur = TrueQ[Lookup[compiled, "HasDurationOptimization", False]];
    durMode = Lookup[compiled, "DurationMode", If[hasDur, "Slices", "None"]];
    activeIdx = Lookup[compiled, "ActiveDurationIndices", {}];
    freezeIdx = Lookup[compiled, "FreezeControlSlices", {}];
    hasPhi = TrueQ[Lookup[compiled, "HasPhaseFreedom", False]];
    pfSpec = Lookup[compiled, "PhaseFreedomSpec", None];
    nPhi = If[hasPhi, Lookup[compiled, "NumPhases", 0], 0];
    If[freezeIdx =!= {} && !basisAssoc["IsIdentityQ"],
        Return[$Failed]
    ];
    If[freezeIdx =!= {},
        uMat = ToPackedRealMatrix[uMat];
        uMat[[freezeIdx, All]] = 0.
    ];

    If[hasDur && !AllTrue[dtVec, TrueQ[# > 0.] &],
        Return[$Failed]
    ];

    nStages = Lookup[compiled, "NumStages", 1];
    If[hasDur,
        If[Lookup[compiled, "DriftKind", "Static"] === "Continuous",
            driftSched = buildFineDriftSchedule[
                compiled["DriftEvaluator"], dtVec, compiled["InitialTimePoint"],
                compiled["DriftIntegrationSpec"], space
            ];
            If[driftSched === $Failed, Return[$Failed]];
            driftCur = driftSched["Drift"];
            dtFine = driftSched["FineTimeSteps"];
            timingSpec = buildTimingDerivTensors[compiled, dtVec, driftSched],
            driftCur = compiled["Drift"];
            dtFine = dtVec;
            timingSpec = buildTimingDerivTensors[compiled, dtVec, None]
        ],
        driftCur = compiled["Drift"];
        dtFine = Lookup[compiled, "FineTimeSteps", dtVec];
        timingSpec = None
    ];
    rcSpec = Lookup[compiled, "RunningCostSpec", None];
    objCotFn = If[hasPhi,
        With[{ewp = pfSpec["EvaluateWithPhases"], pv = phiVec},
            Function[uFinal,
                With[{r = ewp[uFinal, pv]}, {r[[1]], r[[2]], r[[3]]}]
            ]
        ],
        compiled["Objective", "EvaluateCotangent"]
    ];

    kernRes = Switch[method,
        "Spectral",
            EvaluateHilbertSpectral[
                driftCur, compiled["ControlRows"], compiled["ControlFlatT"],
                dtFine, uMat, objCotFn, returnProps, nStages, rcSpec, timingSpec
            ],
        "BlockMatrixExp",
            EvaluateHilbertBlockExp[
                driftCur, compiled["Controls"], dtFine, uMat,
                objCotFn, returnProps, nStages, rcSpec, timingSpec
            ],
        "FirstOrder",
            EvaluateHilbertFirstOrder[
                driftCur, compiled["Controls"], dtFine, uMat,
                objCotFn, returnProps, nStages, rcSpec, timingSpec
            ],
        "StateVector",
            With[{pk = compiled["Objective", "PureKetFastPath"]},
                EvaluateHilbertStateVector[
                    driftCur, compiled["ControlRows"], compiled["ControlFlatT"],
                    dtFine, uMat, pk["PsiMat0"], pk["TarMat"], pk["Weights"], pk["PhaseInsensitive"], nStages, rcSpec, timingSpec
                ]
            ],
        "Liouville",
            With[{ls = compiled["Objective"]},
                EvaluateLiouvilleExact[
                    driftCur, compiled["Controls"], dtFine, uMat,
                    ls["InitialMatrix"], ls["TargetMatrix"], ls["Weights"], ls["Modes"], ls["Scales"], nStages, rcSpec, timingSpec
                ]
            ]
    ];
    If[!AssociationQ[kernRes], Return[$Failed]];

    fid = kernRes["Fidelity"];
    runVal = Lookup[kernRes, "RunningCost", 0.];
    gradU = kernRes["Gradient"];
    gDtKern = Lookup[kernRes, "DurationGradient", None];
    gPhi = If[hasPhi,
        With[{ext = Lookup[kernRes, "Extra", None]},
            If[AssociationQ[ext], Lookup[ext, "PhaseGradient", ConstantArray[0., nPhi]], ConstantArray[0., nPhi]]
        ],
        None
    ];

    lamP = compiled["PowerPenalty"];
    lamS = compiled["SmoothnessPenalty"];
    lamT = Lookup[compiled, "DurationPenaltyLambda", 0.];
    tRef = Max[Lookup[compiled, "DurationPenaltyReference", 1.], 10^-12];
    penSmooth = 0.;
    gRegU = ConstantArray[0., {nt, mCtrl}];
    gDtPen = ConstantArray[0., nt];
    If[lamP > 0.,
        With[{uSq = Total[uMat * uMat, {2}]},
            penSmooth += lamP * Total[dtVec * uSq];
            gRegU += (-2. * lamP) * (dtVec * uMat);
            If[hasDur, gDtPen -= lamP * uSq]
        ]
    ];
    If[lamS > 0. && nt > 1,
        diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]);
        du = Rest[uMat] - Most[uMat];
        duSq = Total[du * du, {2}];
        penSmooth += lamS * Total[duSq / diffDt];
        gRegU[[1 ;; nt - 1]] += (2. * lamS) * (du / diffDt);
        gRegU[[2 ;; nt]] -= (2. * lamS) * (du / diffDt);
        If[hasDur,
            gDtPen[[1 ;; nt - 1]] += (0.5 * lamS) * (duSq / (diffDt^2));
            gDtPen[[2 ;; nt]] += (0.5 * lamS) * (duSq / (diffDt^2))
        ]
    ];
    If[lamT > 0.,
        penSmooth += lamT * (Total[dtVec] / tRef);
        If[hasDur, gDtPen -= ConstantArray[lamT / tRef, nt]]
    ];
    penNS = If[KeyExistsQ[compiled, "ProxContext"], evalNonSmoothPenalty[compiled["ProxContext"], uMat], 0.];
    penTotal = penSmooth + penNS;
    gradUFid = gradU;
    If[freezeIdx =!= {}, gradUFid[[freezeIdx, All]] = 0.];
    gParamFid = If[basisAssoc["IsIdentityQ"],
        gradUFid,
        Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gradUFid]
    ];
    gU = Developer`ToPackedArray[gradU + gRegU];
    If[freezeIdx =!= {}, gU[[freezeIdx, All]] = 0.];
    gParam = If[basisAssoc["IsIdentityQ"],
        gU,
        Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gU]
    ];
    gRegParam = Developer`ToPackedArray[gParam - gParamFid];
    gDt = If[hasDur, Developer`ToPackedArray[gDtKern + gDtPen], None];
    gDurActive = Which[
        !hasDur, {},
        durMode === "Total", Developer`ToPackedArray[{(compiled["DurationWeights"]) . gDt}],
        True, Developer`ToPackedArray[gDt[[activeIdx]]]
    ];
    gDurFidActive = Which[
        !hasDur, {},
        durMode === "Total", Developer`ToPackedArray[{(compiled["DurationWeights"]) . gDtKern}],
        True, Developer`ToPackedArray[gDtKern[[activeIdx]]]
    ];
    gRegDurActive = Developer`ToPackedArray[gDurActive - gDurFidActive];
    gTotDur = If[hasDur && durMode === "Total", First[gDurActive], None];
    gradOut = If[(hasDur || hasPhi) && jointFlag,
        Developer`ToPackedArray[Join[
            Flatten[gParam],
            gDurActive,
            If[hasPhi, gPhi, {}]
        ]],
        gParam
    ];

    Join[kernRes, <|
        "Fidelity" -> fid,
        "TerminalFidelity" -> fid,
        "RunningCost" -> runVal,
        "SmoothObjective" -> fid + runVal - penSmooth,
        "SmoothPenalty" -> penSmooth,
        "NonSmoothPenalty" -> penNS,
        "Objective" -> fid + runVal - penTotal,
        "Penalty" -> penTotal,
        "Gradient" -> gradOut,
        "ControlGradient" -> gParam,
        "SliceGradient" -> gU,
        "FidelitySliceGradient" -> gradUFid,
        "DurationGradient" -> gDt,
        "ActiveDurationGradient" -> gDurActive,
        "TotalDurationGradient" -> gTotDur,
        "FidelityControlGradient" -> gParamFid,
        "FidelityActiveDurationGradient" -> gDurFidActive,
        "SmoothPenaltyControlGradient" -> gRegParam,
        "SmoothPenaltyActiveDurationGradient" -> gRegDurActive,
        "PhaseGradient" -> gPhi,
        "Phases" -> If[hasPhi, phiVec, None],
        "Parameters" -> pMat,
        "Controls" -> uMat,
        "TimeSteps" -> dtVec,
        "TotalDuration" -> Total[dtVec]
    |>]
];

(* Phase C second-order capability check, cached evaluator, smooth-penalty second derivatives, HVP & Hessian *)
checkSingleSecondOrderCapabilities[comp_Association] := Module[
    {hasDur, hasPhi, hasRC, method, hasNS, objOK, badReason, isSpec},
    hasDur = TrueQ[Lookup[comp, "HasDurationOptimization", False]];
    hasPhi = TrueQ[Lookup[comp, "HasPhaseFreedom", False]];
    hasRC = AssociationQ[Lookup[comp, "RunningCostSpec", None]];
    method = Lookup[comp, "Method", "Unknown"];
    hasNS = TrueQ[Lookup[comp, "L1Penalty", 0.] > 0. || Lookup[comp, "TVPenalty", 0.] > 0.];
    objOK = TrueQ[Lookup[comp["Objective"], "ExactHessianSupportedQ", False]];
    isSpec = (method === "Spectral");
    badReason = Which[
        hasDur, "Second-order derivatives are not yet enabled for OptimizeDurations -> True.",
        hasPhi, "Second-order derivatives are not yet enabled for parametric PhaseFreedom.",
        hasRC, "Second-order derivatives are not yet enabled for RunningCost / StatePenalty.",
        hasNS, "Nonsmooth L1Penalty / TVPenalty problems do not have a smooth second-order Hessian.",
        !MemberQ[{"Spectral", "StateVector", "Liouville"}, method],
            "Second-order derivatives require Method -> \"Spectral\", \"StateVector\", or \"Liouville\" (got " <> ToString[method] <> ").",
        !objOK,
            Lookup[comp["Objective"], "UnsupportedSecondOrderReason", "Objective does not support exact second derivatives."],
        True, None
    ];
    <|
        "SupportsSecondOrderQ" -> (badReason === None),
        "SupportsFastExplicitHessianQ" -> (badReason === None && isSpec),
        "Reason" -> badReason,
        "Method" -> method
    |>
];

checkSecondOrderCapabilities[compiled_Association] := If[
    KeyExistsQ[$OCSolverRegistry, "IsEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[compiled]],
    If[KeyExistsQ[$OCSolverRegistry, "CheckEnsembleSecondOrderCapabilities"],
        ($OCSolverRegistry["CheckEnsembleSecondOrderCapabilities"])[compiled],
        <|"SupportsSecondOrderQ" -> False, "SupportsFastExplicitHessianQ" -> False, "Reason" -> "Ensemble second-order handler not registered."|>
    ],
    checkSingleSecondOrderCapabilities[compiled]
];

applySmoothPenaltyHVP[comp_Association, vSliceMat_, dtVec_List] := Module[
    {nt = Length[vSliceMat], mCtrl = Dimensions[vSliceMat][[2]], lamP, lamS, hvp, diffDt, dv},
    lamP = Lookup[comp, "PowerPenalty", 0.];
    lamS = Lookup[comp, "SmoothnessPenalty", 0.];
    hvp = ConstantArray[0., {nt, mCtrl}];
    If[lamP > 0.,
        hvp += (-2. * lamP) * (dtVec * vSliceMat)
    ];
    If[lamS > 0. && nt > 1,
        diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]);
        dv = Rest[vSliceMat] - Most[vSliceMat];
        hvp[[1 ;; nt - 1]] += (2. * lamS) * (dv / diffDt);
        hvp[[2 ;; nt]] -= (2. * lamS) * (dv / diffDt)
    ];
    Developer`ToPackedArray[hvp]
];

buildSmoothPenaltyHessian[comp_Association, dtVec_List] := Module[
    {nt = comp["NumSlices"], mCtrl = comp["NumControls"], nTot, lamP, lamS, hMat, diffDt, j, k, wS, i0, i1},
    nTot = nt * mCtrl;
    lamP = Lookup[comp, "PowerPenalty", 0.];
    lamS = Lookup[comp, "SmoothnessPenalty", 0.];
    hMat = ConstantArray[0., {nTot, nTot}];
    If[lamP > 0.,
        hMat -= 2. * lamP * DiagonalMatrix[Developer`ToPackedArray[Flatten[ConstantArray[#, mCtrl] & /@ dtVec]]]
    ];
    If[lamS > 0. && nt > 1,
        diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]);
        Do[
            wS = (2. * lamS) / diffDt[[j]];
            Do[
                i0 = (j - 1) * mCtrl + k;
                i1 = j * mCtrl + k;
                hMat[[i0, i0]] -= wS;
                hMat[[i1, i1]] -= wS;
                hMat[[i0, i1]] += wS;
                hMat[[i1, i0]] += wS,
                {k, mCtrl}
            ],
            {j, nt - 1}
        ]
    ];
    Developer`ToPackedArray[hMat]
];

pullbackSliceHessianToBasis[basisAssoc_Association, hSliceMat_, mCtrl_Integer] := Module[
    {nt, nP, bMat, bT, h1, h2, h3},
    If[TrueQ[basisAssoc["IsIdentityQ"]], Return[Developer`ToPackedArray[hSliceMat]]];
    bMat = Developer`ToPackedArray[N[basisAssoc["BasisMatrix"]]];
    bT = Developer`ToPackedArray[N[basisAssoc["BasisTranspose"]]];
    {nP, nt} = Dimensions[bT];
    h1 = bT . ArrayReshape[hSliceMat, {nt, mCtrl * nt * mCtrl}];
    h2 = ArrayReshape[Transpose[ArrayReshape[h1, {nP * mCtrl, nt, mCtrl}], {1, 3, 2}], {nP * mCtrl * mCtrl, nt}] . bMat;
    h3 = Transpose[ArrayReshape[h2, {nP * mCtrl, mCtrl, nP}], {1, 3, 2}];
    Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ ArrayReshape[h3, {nP * mCtrl, nP * mCtrl}]]
];

evaluateOptimalControlWithCache[compiled_Association, input_] := Module[
    {parsed, pMat, uMat, dtVec, phiVec, jointFlag, basisAssoc, space, method, nt, mCtrl,
     freezeIdx, caps, needCache, rcSpec, objCotFn, driftCur, dtFine, nStages,
     kernRes, fid, runVal, gradU, gradUFid, lamP, lamS, penSmooth, penNS, penTotal, gRegU,
     gU, gParam, gParamFid, gRegParam, cacheState},
    If[KeyExistsQ[$OCSolverRegistry, "IsEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[compiled]],
        If[KeyExistsQ[$OCSolverRegistry, "EvaluateEnsembleWithCache"],
            Return[($OCSolverRegistry["EvaluateEnsembleWithCache"])[compiled, input]],
            Return[$Failed]
        ]
    ];
    parsed = unpackControlMatrix[compiled, input];
    If[parsed === $Failed, Return[$Failed]];
    {pMat, uMat, dtVec, phiVec, jointFlag} = parsed;
    basisAssoc = compiled["ControlBasis"];
    space = compiled["Space"];
    method = compiled["Method"];
    nt = compiled["NumSlices"];
    mCtrl = compiled["NumControls"];
    freezeIdx = Lookup[compiled, "FreezeControlSlices", {}];
    If[freezeIdx =!= {} && !basisAssoc["IsIdentityQ"], Return[$Failed]];
    If[freezeIdx =!= {},
        uMat = ToPackedRealMatrix[uMat];
        uMat[[freezeIdx, All]] = 0.
    ];

    caps = checkSingleSecondOrderCapabilities[compiled];
    needCache = TrueQ[caps["SupportsSecondOrderQ"]];
    driftCur = compiled["Drift"];
    dtFine = Lookup[compiled, "FineTimeSteps", dtVec];
    nStages = Lookup[compiled, "NumStages", 1];
    rcSpec = Lookup[compiled, "RunningCostSpec", None];
    objCotFn = If[TrueQ[Lookup[compiled, "HasPhaseFreedom", False]],
        With[{ewp = compiled["PhaseFreedomSpec", "EvaluateWithPhases"], pv = phiVec},
            Function[uFinal, With[{r = ewp[uFinal, pv]}, {r[[1]], r[[2]], r[[3]]}]]
        ],
        compiled["Objective", "EvaluateCotangent"]
    ];

    kernRes = Switch[method,
        "Spectral",
            EvaluateHilbertSpectral[
                driftCur, compiled["ControlRows"], compiled["ControlFlatT"],
                dtFine, uMat, objCotFn, False, nStages, rcSpec, None, needCache
            ],
        "BlockMatrixExp",
            EvaluateHilbertBlockExp[
                driftCur, compiled["Controls"], dtFine, uMat,
                objCotFn, False, nStages, rcSpec, None
            ],
        "FirstOrder",
            EvaluateHilbertFirstOrder[
                driftCur, compiled["Controls"], dtFine, uMat,
                objCotFn, False, nStages, rcSpec, None
            ],
        "StateVector",
            With[{pk = compiled["Objective", "PureKetFastPath"]},
                EvaluateHilbertStateVector[
                    driftCur, compiled["ControlRows"], compiled["ControlFlatT"],
                    dtFine, uMat, pk["PsiMat0"], pk["TarMat"], pk["Weights"], pk["PhaseInsensitive"], nStages, rcSpec, None, needCache
                ]
            ],
        "Liouville",
            With[{ls = compiled["Objective"]},
                EvaluateLiouvilleExact[
                    driftCur, compiled["Controls"], dtFine, uMat,
                    ls["InitialMatrix"], ls["TargetMatrix"], ls["Weights"], ls["Modes"], ls["Scales"], nStages, rcSpec, None, needCache
                ]
            ]
    ];
    If[!AssociationQ[kernRes], Return[$Failed]];

    fid = kernRes["Fidelity"];
    runVal = Lookup[kernRes, "RunningCost", 0.];
    gradU = kernRes["Gradient"];
    lamP = compiled["PowerPenalty"];
    lamS = compiled["SmoothnessPenalty"];
    penSmooth = computeSmoothPenaltyValue[compiled, uMat, dtVec];
    gRegU = ConstantArray[0., {nt, mCtrl}];
    If[lamP > 0., gRegU += (-2. * lamP) * (dtVec * uMat)];
    If[lamS > 0. && nt > 1,
        With[{diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]), du = Rest[uMat] - Most[uMat]},
            gRegU[[1 ;; nt - 1]] += (2. * lamS) * (du / diffDt);
            gRegU[[2 ;; nt]] -= (2. * lamS) * (du / diffDt)
        ]
    ];
    penNS = If[KeyExistsQ[compiled, "ProxContext"], evalNonSmoothPenalty[compiled["ProxContext"], uMat], 0.];
    penTotal = penSmooth + penNS;
    gradUFid = gradU;
    If[freezeIdx =!= {}, gradUFid[[freezeIdx, All]] = 0.];
    gParamFid = If[basisAssoc["IsIdentityQ"],
        gradUFid,
        Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gradUFid]
    ];
    gU = Developer`ToPackedArray[gradU + gRegU];
    If[freezeIdx =!= {}, gU[[freezeIdx, All]] = 0.];
    gParam = If[basisAssoc["IsIdentityQ"],
        gU,
        Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gU]
    ];
    gRegParam = Developer`ToPackedArray[gParam - gParamFid];

    cacheState = If[needCache && KeyExistsQ[kernRes, "EvalState"] && AssociationQ[kernRes["EvalState"]],
        Join[kernRes["EvalState"], <|
            "SupportsSecondOrderQ" -> True,
            "BasisAssoc" -> basisAssoc,
            "ControlsMat" -> uMat,
            "ParametersMat" -> pMat,
            "TimeSteps" -> dtVec
        |>],
        <|"SupportsSecondOrderQ" -> False, "Reason" -> caps["Reason"]|>
    ];

    Join[KeyDrop[kernRes, {"EvalState"}], <|
        "Fidelity" -> fid,
        "TerminalFidelity" -> fid,
        "RunningCost" -> runVal,
        "SmoothObjective" -> fid + runVal - penSmooth,
        "SmoothPenalty" -> penSmooth,
        "NonSmoothPenalty" -> penNS,
        "Objective" -> fid + runVal - penTotal,
        "Penalty" -> penTotal,
        "Gradient" -> gParam,
        "ControlGradient" -> gParam,
        "SliceGradient" -> gU,
        "FidelitySliceGradient" -> gradUFid,
        "FidelityControlGradient" -> gParamFid,
        "SmoothPenaltyControlGradient" -> gRegParam,
        "Phases" -> None,
        "Parameters" -> pMat,
        "Controls" -> uMat,
        "TimeSteps" -> dtVec,
        "TotalDuration" -> Total[dtVec],
        "SupportsSecondOrderQ" -> caps["SupportsSecondOrderQ"],
        "UnsupportedSecondOrderReason" -> caps["Reason"],
        "EvalState" -> cacheState
    |>]
];

evaluateOptimalControlHVP[compiled_Association, evalState_Association, dirInput_] := Module[
    {caps, basisAssoc, nParams, mCtrl, nt, isVec, vParamMat, vSliceMat, freezeIdx,
     method, hFidSlice, hPenSlice, hTotSlice, hParamMat},
    If[KeyExistsQ[$OCSolverRegistry, "IsEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[compiled]],
        If[KeyExistsQ[$OCSolverRegistry, "EvaluateEnsembleHVP"],
            Return[($OCSolverRegistry["EvaluateEnsembleHVP"])[compiled, evalState, dirInput]],
            Return[$Failed]
        ]
    ];
    caps = checkSingleSecondOrderCapabilities[compiled];
    If[!TrueQ[caps["SupportsSecondOrderQ"]] || !TrueQ[Lookup[evalState, "SupportsSecondOrderQ", False]],
        Return[$Failed]
    ];
    basisAssoc = compiled["ControlBasis"];
    nParams = compiled["NumParams"];
    mCtrl = compiled["NumControls"];
    nt = compiled["NumSlices"];
    freezeIdx = Lookup[compiled, "FreezeControlSlices", {}];
    isVec = VectorQ[dirInput, NumericQ];
    vParamMat = Which[
        isVec && Length[dirInput] === nParams * mCtrl,
            Partition[Developer`ToPackedArray[N[Re[dirInput]]], mCtrl],
        MatrixQ[dirInput, NumericQ] && Dimensions[dirInput] === {nParams, mCtrl},
            ToPackedRealMatrix[dirInput],
        True,
            Return[$Failed]
    ];
    vSliceMat = If[TrueQ[basisAssoc["IsIdentityQ"]],
        vParamMat,
        Developer`ToPackedArray[basisAssoc["BasisMatrix"] . vParamMat]
    ];
    If[freezeIdx =!= {}, vSliceMat[[freezeIdx, All]] = 0.];

    method = caps["Method"];
    hFidSlice = Switch[method,
        "Spectral",
            OptimalControl`Kernels`Private`EvaluateHilbertSpectralHVP[
                evalState, vSliceMat, compiled["Objective", "EvaluateCotangentDirectional"]
            ],
        "StateVector",
            OptimalControl`Kernels`Private`EvaluateHilbertStateVectorHVP[evalState, vSliceMat],
        "Liouville",
            OptimalControl`Kernels`Private`EvaluateLiouvilleHVP[evalState, vSliceMat],
        _,
            $Failed
    ];
    If[!MatrixQ[hFidSlice, NumericQ], Return[$Failed]];
    hPenSlice = applySmoothPenaltyHVP[compiled, vSliceMat, evalState["TimeSteps"]];
    hTotSlice = Developer`ToPackedArray[hFidSlice + hPenSlice];
    If[freezeIdx =!= {}, hTotSlice[[freezeIdx, All]] = 0.];
    hParamMat = If[TrueQ[basisAssoc["IsIdentityQ"]],
        hTotSlice,
        Developer`ToPackedArray[basisAssoc["BasisTranspose"] . hTotSlice]
    ];
    If[isVec, Developer`ToPackedArray[Flatten[hParamMat]], hParamMat]
];

evaluateOptimalControlHessian[compiled_Association, evalState_Association, mode_: Automatic] := Module[
    {caps, basisAssoc, nParams, mCtrl, nDim, freezeIdx, dtVec, method, hSliceFid, hSlicePen, hSliceTot,
     i, eVec, cols, j, k, idx},
    If[KeyExistsQ[$OCSolverRegistry, "IsEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[compiled]],
        If[KeyExistsQ[$OCSolverRegistry, "EvaluateEnsembleHessian"],
            Return[($OCSolverRegistry["EvaluateEnsembleHessian"])[compiled, evalState, mode]],
            Return[$Failed]
        ]
    ];
    caps = checkSingleSecondOrderCapabilities[compiled];
    If[!TrueQ[caps["SupportsSecondOrderQ"]] || !TrueQ[Lookup[evalState, "SupportsSecondOrderQ", False]],
        Return[$Failed]
    ];
    basisAssoc = compiled["ControlBasis"];
    nParams = compiled["NumParams"];
    mCtrl = compiled["NumControls"];
    nDim = nParams * mCtrl;
    freezeIdx = Lookup[compiled, "FreezeControlSlices", {}];
    dtVec = evalState["TimeSteps"];
    method = caps["Method"];

    If[method === "Spectral" && mode =!= "ColumnHVP",
        hSliceFid = OptimalControl`Kernels`Private`EvaluateHilbertSpectralHessian[
            evalState, compiled["Controls"], compiled["Objective", "EvaluateCotangentDirectional"]
        ];
        If[!MatrixQ[hSliceFid, NumericQ], Return[$Failed]];
        hSlicePen = buildSmoothPenaltyHessian[compiled, dtVec];
        hSliceTot = Developer`ToPackedArray[hSliceFid + hSlicePen];
        If[freezeIdx =!= {},
            Do[
                Do[
                    idx = (j - 1) * mCtrl + k;
                    hSliceTot[[idx, All]] = 0.;
                    hSliceTot[[All, idx]] = 0.,
                    {k, mCtrl}
                ],
                {j, freezeIdx}
            ]
        ];
        Return[pullbackSliceHessianToBasis[basisAssoc, hSliceTot, mCtrl]]
    ];

    (* Column-by-column assembly via cached HVP *)
    eVec = ConstantArray[0., nDim];
    cols = Table[
        eVec[[i]] = 1.;
        With[{c = evaluateOptimalControlHVP[compiled, evalState, eVec]},
            eVec[[i]] = 0.;
            If[!VectorQ[c, NumericQ], Return[$Failed, Module]];
            c
        ],
        {i, nDim}
    ];
    Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ Transpose[cols]]
];

computeSmoothPenaltyValue[comp_Association, uMat_, dtVec_List] := Module[
    {nt = Length[uMat], lamP, lamS, lamT, tRef, penSmooth = 0., diffDt, du, duSq},
    lamP = comp["PowerPenalty"];
    lamS = comp["SmoothnessPenalty"];
    lamT = Lookup[comp, "DurationPenaltyLambda", 0.];
    tRef = Max[Lookup[comp, "DurationPenaltyReference", 1.], 10^-12];
    If[lamP > 0.,
        With[{uSq = Total[uMat * uMat, {2}]},
            penSmooth += lamP * Total[dtVec * uSq]
        ]
    ];
    If[lamS > 0. && nt > 1,
        diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]);
        du = Rest[uMat] - Most[uMat];
        duSq = Total[du * du, {2}];
        penSmooth += lamS * Total[duSq / diffDt]
    ];
    If[lamT > 0.,
        penSmooth += lamT * (Total[dtVec] / tRef)
    ];
    penSmooth
];

(* Forward-only objective evaluation (no gradient calculation) *)
evaluateOptimalControlObjective[compiled_Association, input_] := Module[
    {parsed, pMat, uMat, dtVec, phiVec, jointFlag, basisAssoc, space, method, nt,
     hasDur, freezeIdx, hasPhi, pfSpec, objCotFn, driftCur, dtFine, nStages,
     driftSched, rcSpec, kernRes, fid, runVal, penSmooth, penNS, penTotal},
    If[KeyExistsQ[$OCSolverRegistry, "IsEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[compiled]],
        If[KeyExistsQ[$OCSolverRegistry, "EvaluateEnsembleObjective"],
            Return[($OCSolverRegistry["EvaluateEnsembleObjective"])[compiled, input]],
            Return[($OCSolverRegistry["EvaluateEnsemble"])[compiled, input]]
        ]
    ];
    parsed = unpackControlMatrix[compiled, input];
    If[parsed === $Failed, Return[$Failed]];
    {pMat, uMat, dtVec, phiVec, jointFlag} = parsed;
    basisAssoc = compiled["ControlBasis"];
    space = compiled["Space"];
    method = compiled["Method"];
    nt = compiled["NumSlices"];
    hasDur = TrueQ[Lookup[compiled, "HasDurationOptimization", False]];
    freezeIdx = Lookup[compiled, "FreezeControlSlices", {}];
    hasPhi = TrueQ[Lookup[compiled, "HasPhaseFreedom", False]];
    pfSpec = Lookup[compiled, "PhaseFreedomSpec", None];
    If[freezeIdx =!= {} && !basisAssoc["IsIdentityQ"], Return[$Failed]];
    If[freezeIdx =!= {},
        uMat = ToPackedRealMatrix[uMat];
        uMat[[freezeIdx, All]] = 0.
    ];
    If[hasDur && !AllTrue[dtVec, TrueQ[# > 0.] &], Return[$Failed]];

    nStages = Lookup[compiled, "NumStages", 1];
    If[hasDur,
        If[Lookup[compiled, "DriftKind", "Static"] === "Continuous",
            driftSched = buildFineDriftSchedule[
                compiled["DriftEvaluator"], dtVec, compiled["InitialTimePoint"],
                compiled["DriftIntegrationSpec"], space
            ];
            If[driftSched === $Failed, Return[$Failed]];
            driftCur = driftSched["Drift"];
            dtFine = driftSched["FineTimeSteps"],
            driftCur = compiled["Drift"];
            dtFine = dtVec
        ],
        driftCur = compiled["Drift"];
        dtFine = Lookup[compiled, "FineTimeSteps", dtVec]
    ];
    rcSpec = Lookup[compiled, "RunningCostSpec", None];
    objCotFn = If[hasPhi,
        With[{ewp = pfSpec["EvaluateWithPhases"], pv = phiVec},
            Function[uFinal, With[{r = ewp[uFinal, pv]}, {r[[1]], r[[2]], r[[3]]}]]
        ],
        compiled["Objective", "EvaluateCotangent"]
    ];

    kernRes = Switch[method,
        "Spectral" | "BlockMatrixExp" | "FirstOrder",
            OptimalControl`Kernels`Private`EvaluateHilbertSpectralObjectiveOnly[
                driftCur, compiled["ControlRows"], dtFine, uMat, objCotFn, nStages, rcSpec
            ],
        "StateVector",
            With[{pk = compiled["Objective", "PureKetFastPath"]},
                OptimalControl`Kernels`Private`EvaluateHilbertStateVectorObjectiveOnly[
                    driftCur, compiled["ControlRows"], dtFine, uMat, pk["PsiMat0"], pk["TarMat"],
                    pk["Weights"], pk["PhaseInsensitive"], nStages, rcSpec
                ]
            ],
        "Liouville",
            With[{ls = compiled["Objective"]},
                OptimalControl`Kernels`Private`EvaluateLiouvilleObjectiveOnly[
                    driftCur, compiled["Controls"], dtFine, uMat,
                    ls["InitialMatrix"], ls["TargetMatrix"], ls["Weights"], ls["Modes"], ls["Scales"],
                    nStages, rcSpec
                ]
            ]
    ];
    If[!AssociationQ[kernRes], Return[$Failed]];
    fid = kernRes["Fidelity"];
    runVal = Lookup[kernRes, "RunningCost", 0.];
    penSmooth = computeSmoothPenaltyValue[compiled, uMat, dtVec];
    penNS = If[KeyExistsQ[compiled, "ProxContext"], evalNonSmoothPenalty[compiled["ProxContext"], uMat], 0.];
    penTotal = penSmooth + penNS;
    Join[kernRes, <|
        "Fidelity" -> fid,
        "TerminalFidelity" -> fid,
        "RunningCost" -> runVal,
        "SmoothObjective" -> fid + runVal - penSmooth,
        "SmoothPenalty" -> penSmooth,
        "NonSmoothPenalty" -> penNS,
        "Objective" -> fid + runVal - penTotal,
        "Penalty" -> penTotal,
        "Phases" -> If[hasPhi, phiVec, None],
        "Parameters" -> pMat,
        "Controls" -> uMat,
        "TimeSteps" -> dtVec,
        "TotalDuration" -> Total[dtVec]
    |>]
];

(* Phase C second-order capability check, cached evaluator, smooth-penalty second derivatives, HVP & Hessian *)
checkSingleSecondOrderCapabilities[comp_Association] := Module[
    {hasDur, hasPhi, hasRC, method, hasNS, objOK, badReason, isSpec},
    hasDur = TrueQ[Lookup[comp, "HasDurationOptimization", False]];
    hasPhi = TrueQ[Lookup[comp, "HasPhaseFreedom", False]];
    hasRC = AssociationQ[Lookup[comp, "RunningCostSpec", None]];
    method = Lookup[comp, "Method", "Unknown"];
    hasNS = TrueQ[Lookup[comp, "L1Penalty", 0.] > 0. || Lookup[comp, "TVPenalty", 0.] > 0.];
    objOK = TrueQ[Lookup[comp["Objective"], "ExactHessianSupportedQ", False]];
    isSpec = (method === "Spectral");
    badReason = Which[
        hasDur, "Second-order derivatives are not yet enabled for OptimizeDurations -> True.",
        hasPhi, "Second-order derivatives are not yet enabled for parametric PhaseFreedom.",
        hasRC, "Second-order derivatives are not yet enabled for RunningCost / StatePenalty.",
        hasNS, "Nonsmooth L1Penalty / TVPenalty problems do not have a smooth second-order Hessian.",
        !MemberQ[{"Spectral", "StateVector", "Liouville"}, method],
            "Second-order derivatives require Method -> \"Spectral\", \"StateVector\", or \"Liouville\" (got " <> ToString[method] <> ").",
        !objOK,
            Lookup[comp["Objective"], "UnsupportedSecondOrderReason", "Objective does not support exact second derivatives."],
        True, None
    ];
    <|
        "SupportsSecondOrderQ" -> (badReason === None),
        "SupportsFastExplicitHessianQ" -> (badReason === None && isSpec),
        "Reason" -> badReason,
        "Method" -> method
    |>
];

checkSecondOrderCapabilities[compiled_Association] := If[
    KeyExistsQ[$OCSolverRegistry, "IsEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[compiled]],
    If[KeyExistsQ[$OCSolverRegistry, "CheckEnsembleSecondOrderCapabilities"],
        ($OCSolverRegistry["CheckEnsembleSecondOrderCapabilities"])[compiled],
        <|"SupportsSecondOrderQ" -> False, "SupportsFastExplicitHessianQ" -> False, "Reason" -> "Ensemble second-order handler not registered."|>
    ],
    checkSingleSecondOrderCapabilities[compiled]
];

applySmoothPenaltyHVP[comp_Association, vSliceMat_, dtVec_List] := Module[
    {nt = Length[vSliceMat], mCtrl = Dimensions[vSliceMat][[2]], lamP, lamS, hvp, diffDt, dv},
    lamP = Lookup[comp, "PowerPenalty", 0.];
    lamS = Lookup[comp, "SmoothnessPenalty", 0.];
    hvp = ConstantArray[0., {nt, mCtrl}];
    If[lamP > 0.,
        hvp += (-2. * lamP) * (dtVec * vSliceMat)
    ];
    If[lamS > 0. && nt > 1,
        diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]);
        dv = Rest[vSliceMat] - Most[vSliceMat];
        hvp[[1 ;; nt - 1]] += (2. * lamS) * (dv / diffDt);
        hvp[[2 ;; nt]] -= (2. * lamS) * (dv / diffDt)
    ];
    Developer`ToPackedArray[hvp]
];

buildSmoothPenaltyHessian[comp_Association, dtVec_List] := Module[
    {nt = comp["NumSlices"], mCtrl = comp["NumControls"], nTot, lamP, lamS, hMat, diffDt, j, k, wS, i0, i1},
    nTot = nt * mCtrl;
    lamP = Lookup[comp, "PowerPenalty", 0.];
    lamS = Lookup[comp, "SmoothnessPenalty", 0.];
    hMat = ConstantArray[0., {nTot, nTot}];
    If[lamP > 0.,
        hMat -= 2. * lamP * DiagonalMatrix[Developer`ToPackedArray[Flatten[ConstantArray[#, mCtrl] & /@ dtVec]]]
    ];
    If[lamS > 0. && nt > 1,
        diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]);
        Do[
            wS = (2. * lamS) / diffDt[[j]];
            Do[
                i0 = (j - 1) * mCtrl + k;
                i1 = j * mCtrl + k;
                hMat[[i0, i0]] -= wS;
                hMat[[i1, i1]] -= wS;
                hMat[[i0, i1]] += wS;
                hMat[[i1, i0]] += wS,
                {k, mCtrl}
            ],
            {j, nt - 1}
        ]
    ];
    Developer`ToPackedArray[hMat]
];

pullbackSliceHessianToBasis[basisAssoc_Association, hSliceMat_, mCtrl_Integer] := Module[
    {nt, nP, bMat, bT, h1, h2, h3},
    If[TrueQ[basisAssoc["IsIdentityQ"]], Return[Developer`ToPackedArray[hSliceMat]]];
    bMat = Developer`ToPackedArray[N[basisAssoc["BasisMatrix"]]];
    bT = Developer`ToPackedArray[N[basisAssoc["BasisTranspose"]]];
    {nP, nt} = Dimensions[bT];
    h1 = bT . ArrayReshape[hSliceMat, {nt, mCtrl * nt * mCtrl}];
    h2 = ArrayReshape[Transpose[ArrayReshape[h1, {nP * mCtrl, nt, mCtrl}], {1, 3, 2}], {nP * mCtrl * mCtrl, nt}] . bMat;
    h3 = Transpose[ArrayReshape[h2, {nP * mCtrl, mCtrl, nP}], {1, 3, 2}];
    Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ ArrayReshape[h3, {nP * mCtrl, nP * mCtrl}]]
];

evaluateOptimalControlWithCache[compiled_Association, input_] := Module[
    {parsed, pMat, uMat, dtVec, phiVec, jointFlag, basisAssoc, space, method, nt, mCtrl,
     freezeIdx, caps, needCache, rcSpec, objCotFn, driftCur, dtFine, nStages,
     kernRes, fid, runVal, gradU, gradUFid, lamP, lamS, penSmooth, penNS, penTotal, gRegU,
     gU, gParam, gParamFid, gRegParam, cacheState},
    If[KeyExistsQ[$OCSolverRegistry, "IsEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[compiled]],
        If[KeyExistsQ[$OCSolverRegistry, "EvaluateEnsembleWithCache"],
            Return[($OCSolverRegistry["EvaluateEnsembleWithCache"])[compiled, input]],
            Return[$Failed]
        ]
    ];
    parsed = unpackControlMatrix[compiled, input];
    If[parsed === $Failed, Return[$Failed]];
    {pMat, uMat, dtVec, phiVec, jointFlag} = parsed;
    basisAssoc = compiled["ControlBasis"];
    space = compiled["Space"];
    method = compiled["Method"];
    nt = compiled["NumSlices"];
    mCtrl = compiled["NumControls"];
    freezeIdx = Lookup[compiled, "FreezeControlSlices", {}];
    If[freezeIdx =!= {} && !basisAssoc["IsIdentityQ"], Return[$Failed]];
    If[freezeIdx =!= {},
        uMat = ToPackedRealMatrix[uMat];
        uMat[[freezeIdx, All]] = 0.
    ];

    caps = checkSingleSecondOrderCapabilities[compiled];
    needCache = TrueQ[caps["SupportsSecondOrderQ"]];
    driftCur = compiled["Drift"];
    dtFine = Lookup[compiled, "FineTimeSteps", dtVec];
    nStages = Lookup[compiled, "NumStages", 1];
    rcSpec = Lookup[compiled, "RunningCostSpec", None];
    objCotFn = If[TrueQ[Lookup[compiled, "HasPhaseFreedom", False]],
        With[{ewp = compiled["PhaseFreedomSpec", "EvaluateWithPhases"], pv = phiVec},
            Function[uFinal, With[{r = ewp[uFinal, pv]}, {r[[1]], r[[2]], r[[3]]}]]
        ],
        compiled["Objective", "EvaluateCotangent"]
    ];

    kernRes = Switch[method,
        "Spectral",
            EvaluateHilbertSpectral[
                driftCur, compiled["ControlRows"], compiled["ControlFlatT"],
                dtFine, uMat, objCotFn, False, nStages, rcSpec, None, needCache
            ],
        "BlockMatrixExp",
            EvaluateHilbertBlockExp[
                driftCur, compiled["Controls"], dtFine, uMat,
                objCotFn, False, nStages, rcSpec, None
            ],
        "FirstOrder",
            EvaluateHilbertFirstOrder[
                driftCur, compiled["Controls"], dtFine, uMat,
                objCotFn, False, nStages, rcSpec, None
            ],
        "StateVector",
            With[{pk = compiled["Objective", "PureKetFastPath"]},
                EvaluateHilbertStateVector[
                    driftCur, compiled["ControlRows"], compiled["ControlFlatT"],
                    dtFine, uMat, pk["PsiMat0"], pk["TarMat"], pk["Weights"], pk["PhaseInsensitive"], nStages, rcSpec, None, needCache
                ]
            ],
        "Liouville",
            With[{ls = compiled["Objective"]},
                EvaluateLiouvilleExact[
                    driftCur, compiled["Controls"], dtFine, uMat,
                    ls["InitialMatrix"], ls["TargetMatrix"], ls["Weights"], ls["Modes"], ls["Scales"], nStages, rcSpec, None, needCache
                ]
            ]
    ];
    If[!AssociationQ[kernRes], Return[$Failed]];

    fid = kernRes["Fidelity"];
    runVal = Lookup[kernRes, "RunningCost", 0.];
    gradU = kernRes["Gradient"];
    lamP = compiled["PowerPenalty"];
    lamS = compiled["SmoothnessPenalty"];
    penSmooth = computeSmoothPenaltyValue[compiled, uMat, dtVec];
    gRegU = ConstantArray[0., {nt, mCtrl}];
    If[lamP > 0., gRegU += (-2. * lamP) * (dtVec * uMat)];
    If[lamS > 0. && nt > 1,
        With[{diffDt = 0.5 * (Rest[dtVec] + Most[dtVec]), du = Rest[uMat] - Most[uMat]},
            gRegU[[1 ;; nt - 1]] += (2. * lamS) * (du / diffDt);
            gRegU[[2 ;; nt]] -= (2. * lamS) * (du / diffDt)
        ]
    ];
    penNS = If[KeyExistsQ[compiled, "ProxContext"], evalNonSmoothPenalty[compiled["ProxContext"], uMat], 0.];
    penTotal = penSmooth + penNS;
    gradUFid = gradU;
    If[freezeIdx =!= {}, gradUFid[[freezeIdx, All]] = 0.];
    gParamFid = If[basisAssoc["IsIdentityQ"],
        gradUFid,
        Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gradUFid]
    ];
    gU = Developer`ToPackedArray[gradU + gRegU];
    If[freezeIdx =!= {}, gU[[freezeIdx, All]] = 0.];
    gParam = If[basisAssoc["IsIdentityQ"],
        gU,
        Developer`ToPackedArray[basisAssoc["BasisTranspose"] . gU]
    ];
    gRegParam = Developer`ToPackedArray[gParam - gParamFid];

    cacheState = If[needCache && KeyExistsQ[kernRes, "EvalState"] && AssociationQ[kernRes["EvalState"]],
        Join[kernRes["EvalState"], <|
            "SupportsSecondOrderQ" -> True,
            "BasisAssoc" -> basisAssoc,
            "ControlsMat" -> uMat,
            "ParametersMat" -> pMat,
            "TimeSteps" -> dtVec
        |>],
        <|"SupportsSecondOrderQ" -> False, "Reason" -> caps["Reason"]|>
    ];

    Join[KeyDrop[kernRes, {"EvalState"}], <|
        "Fidelity" -> fid,
        "TerminalFidelity" -> fid,
        "RunningCost" -> runVal,
        "SmoothObjective" -> fid + runVal - penSmooth,
        "SmoothPenalty" -> penSmooth,
        "NonSmoothPenalty" -> penNS,
        "Objective" -> fid + runVal - penTotal,
        "Penalty" -> penTotal,
        "Gradient" -> gParam,
        "ControlGradient" -> gParam,
        "SliceGradient" -> gU,
        "FidelitySliceGradient" -> gradUFid,
        "FidelityControlGradient" -> gParamFid,
        "SmoothPenaltyControlGradient" -> gRegParam,
        "Phases" -> None,
        "Parameters" -> pMat,
        "Controls" -> uMat,
        "TimeSteps" -> dtVec,
        "TotalDuration" -> Total[dtVec],
        "SupportsSecondOrderQ" -> caps["SupportsSecondOrderQ"],
        "UnsupportedSecondOrderReason" -> caps["Reason"],
        "EvalState" -> cacheState
    |>]
];

evaluateOptimalControlHVP[compiled_Association, evalState_Association, dirInput_] := Module[
    {caps, basisAssoc, nParams, mCtrl, nt, isVec, vParamMat, vSliceMat, freezeIdx,
     method, hFidSlice, hPenSlice, hTotSlice, hParamMat},
    If[KeyExistsQ[$OCSolverRegistry, "IsEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[compiled]],
        If[KeyExistsQ[$OCSolverRegistry, "EvaluateEnsembleHVP"],
            Return[($OCSolverRegistry["EvaluateEnsembleHVP"])[compiled, evalState, dirInput]],
            Return[$Failed]
        ]
    ];
    caps = checkSingleSecondOrderCapabilities[compiled];
    If[!TrueQ[caps["SupportsSecondOrderQ"]] || !TrueQ[Lookup[evalState, "SupportsSecondOrderQ", False]],
        Return[$Failed]
    ];
    basisAssoc = compiled["ControlBasis"];
    nParams = compiled["NumParams"];
    mCtrl = compiled["NumControls"];
    nt = compiled["NumSlices"];
    freezeIdx = Lookup[compiled, "FreezeControlSlices", {}];
    isVec = VectorQ[dirInput, NumericQ];
    vParamMat = Which[
        isVec && Length[dirInput] === nParams * mCtrl,
            Partition[Developer`ToPackedArray[N[Re[dirInput]]], mCtrl],
        MatrixQ[dirInput, NumericQ] && Dimensions[dirInput] === {nParams, mCtrl},
            ToPackedRealMatrix[dirInput],
        True,
            Return[$Failed]
    ];
    vSliceMat = If[TrueQ[basisAssoc["IsIdentityQ"]],
        vParamMat,
        Developer`ToPackedArray[basisAssoc["BasisMatrix"] . vParamMat]
    ];
    If[freezeIdx =!= {}, vSliceMat[[freezeIdx, All]] = 0.];

    method = caps["Method"];
    hFidSlice = Switch[method,
        "Spectral",
            OptimalControl`Kernels`Private`EvaluateHilbertSpectralHVP[
                evalState, vSliceMat, compiled["Objective", "EvaluateCotangentDirectional"]
            ],
        "StateVector",
            OptimalControl`Kernels`Private`EvaluateHilbertStateVectorHVP[evalState, vSliceMat],
        "Liouville",
            OptimalControl`Kernels`Private`EvaluateLiouvilleHVP[evalState, vSliceMat],
        _,
            $Failed
    ];
    If[!MatrixQ[hFidSlice, NumericQ], Return[$Failed]];
    hPenSlice = applySmoothPenaltyHVP[compiled, vSliceMat, evalState["TimeSteps"]];
    hTotSlice = Developer`ToPackedArray[hFidSlice + hPenSlice];
    If[freezeIdx =!= {}, hTotSlice[[freezeIdx, All]] = 0.];
    hParamMat = If[TrueQ[basisAssoc["IsIdentityQ"]],
        hTotSlice,
        Developer`ToPackedArray[basisAssoc["BasisTranspose"] . hTotSlice]
    ];
    If[isVec, Developer`ToPackedArray[Flatten[hParamMat]], hParamMat]
];

evaluateOptimalControlHessian[compiled_Association, evalState_Association, mode_: Automatic] := Module[
    {caps, basisAssoc, nParams, mCtrl, nDim, freezeIdx, dtVec, method, hSliceFid, hSlicePen, hSliceTot,
     i, eVec, cols, j, k, idx},
    If[KeyExistsQ[$OCSolverRegistry, "IsEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[compiled]],
        If[KeyExistsQ[$OCSolverRegistry, "EvaluateEnsembleHessian"],
            Return[($OCSolverRegistry["EvaluateEnsembleHessian"])[compiled, evalState, mode]],
            Return[$Failed]
        ]
    ];
    caps = checkSingleSecondOrderCapabilities[compiled];
    If[!TrueQ[caps["SupportsSecondOrderQ"]] || !TrueQ[Lookup[evalState, "SupportsSecondOrderQ", False]],
        Return[$Failed]
    ];
    basisAssoc = compiled["ControlBasis"];
    nParams = compiled["NumParams"];
    mCtrl = compiled["NumControls"];
    nDim = nParams * mCtrl;
    freezeIdx = Lookup[compiled, "FreezeControlSlices", {}];
    dtVec = evalState["TimeSteps"];
    method = caps["Method"];

    If[method === "Spectral" && mode =!= "ColumnHVP",
        hSliceFid = OptimalControl`Kernels`Private`EvaluateHilbertSpectralHessian[
            evalState, compiled["Controls"], compiled["Objective", "EvaluateCotangentDirectional"]
        ];
        If[!MatrixQ[hSliceFid, NumericQ], Return[$Failed]];
        hSlicePen = buildSmoothPenaltyHessian[compiled, dtVec];
        hSliceTot = Developer`ToPackedArray[hSliceFid + hSlicePen];
        If[freezeIdx =!= {},
            Do[
                Do[
                    idx = (j - 1) * mCtrl + k;
                    hSliceTot[[idx, All]] = 0.;
                    hSliceTot[[All, idx]] = 0.,
                    {k, mCtrl}
                ],
                {j, freezeIdx}
            ]
        ];
        Return[pullbackSliceHessianToBasis[basisAssoc, hSliceTot, mCtrl]]
    ];

    (* Column-by-column assembly via cached HVP *)
    eVec = ConstantArray[0., nDim];
    cols = Table[
        eVec[[i]] = 1.;
        With[{c = evaluateOptimalControlHVP[compiled, evalState, eVec]},
            eVec[[i]] = 0.;
            If[!VectorQ[c, NumericQ], Return[$Failed, Module]];
            c
        ],
        {i, nDim}
    ];
    Developer`ToPackedArray[0.5 * (# + Transpose[#]) &@ Transpose[cols]]
];

EvaluateOptimalControl[prob_OptimalControlProblem, input_, returnProps_: False] :=
    With[{comp = CompileOptimalControlProblem[prob]},
        If[comp === $Failed, $Failed, EvaluateOptimalControl[comp, input, returnProps]]
    ];

(* Central finite-difference verification (supports controls, durations, and parametric phases) *)
VerifyOptimalControlGradient[probOrCompiled_, controls_, opts___Rule] := Module[
    {stepSize, tol, methodOpt, comp, hasDur, durMode, hasPhi, evalRes, gAna, nRows, nCols, xMat, gFD,
     xVec, nTotal, nCtrl, nAct, nPhi, activeIdx, stepVec, xPlus, xMinus, fPlus, fMinus, diff, absErr,
     scaleDenom, relErr, ctrlRelErr, durRelErr = 0., phiRelErr = 0., r, col, i, stepC, stepD, stepPhi, outAssoc},
    {stepSize, tol, methodOpt} = {StepSize, Tolerance, Method} /. {opts} /. Options[VerifyOptimalControlGradient];
    comp = Which[
        MatchQ[probOrCompiled, _OptimalControlProblem],
            If[methodOpt === Automatic,
                CompileOptimalControlProblem[probOrCompiled],
                With[{baseOpts = List @@ probOrCompiled[[5 ;;]]},
                    CompileOptimalControlProblem[
                        OptimalControlProblem[probOrCompiled[[1]], probOrCompiled[[2]],
                            probOrCompiled[[3]], probOrCompiled[[4]], Method -> methodOpt, Sequence @@ baseOpts]
                    ]
                ]
            ],
        AssociationQ[probOrCompiled],
            If[methodOpt =!= Automatic && !KeyExistsQ[probOrCompiled, "Members"],
                Join[probOrCompiled, <|"Method" -> methodOpt|>],
                probOrCompiled
            ],
        KeyExistsQ[$OCSolverRegistry, "CompileIfEnsemble"] && TrueQ[($OCSolverRegistry["IsEnsemble"])[probOrCompiled]],
            ($OCSolverRegistry["CompileIfEnsemble"])[probOrCompiled],
        True, $Failed
    ];
    If[comp === $Failed, Return[$Failed]];

    hasDur = TrueQ[Lookup[comp, "HasDurationOptimization", False]];
    durMode = Lookup[comp, "DurationMode", If[hasDur, "Slices", "None"]];
    hasPhi = TrueQ[Lookup[comp, "HasPhaseFreedom", False]];
    If[hasDur || hasPhi,
        activeIdx = Lookup[comp, "ActiveDurationIndices", {}];
        nCtrl = comp["NumParams"] * comp["NumControls"];
        nAct = If[hasDur, Lookup[comp, "NumActiveDurations", Length[activeIdx]], 0];
        nPhi = If[hasPhi, Lookup[comp, "NumPhases", 0], 0];
        xVec = Which[
            VectorQ[controls, NumericQ] && Length[controls] === nCtrl + nAct + nPhi,
                Developer`ToPackedArray[N[Re[controls]]],
            MatchQ[controls, _OptimalControlResult],
                With[{ev0 = EvaluateOptimalControl[comp, controls[[1]]]},
                    If[ev0 === $Failed, Return[$Failed]];
                    Developer`ToPackedArray[Join[
                        Flatten[ev0["Parameters"]],
                        If[hasDur,
                            If[durMode === "Total", {ev0["TotalDuration"]}, (ev0["TimeSteps"])[[activeIdx]]],
                            {}
                        ],
                        If[hasPhi, Lookup[ev0, "Phases", ConstantArray[0., nPhi]], {}]
                    ]]
                ],
            True,
                With[{ev0 = EvaluateOptimalControl[comp, controls]},
                    If[ev0 === $Failed, Return[$Failed]];
                    Developer`ToPackedArray[Join[
                        Flatten[ev0["Parameters"]],
                        If[hasDur,
                            If[durMode === "Total", {ev0["TotalDuration"]}, (ev0["TimeSteps"])[[activeIdx]]],
                            {}
                        ],
                        If[hasPhi, Lookup[ev0, "Phases", ConstantArray[0., nPhi]], {}]
                    ]]
                ]
        ];
        evalRes = EvaluateOptimalControl[comp, xVec];
        If[evalRes === $Failed, Return[$Failed]];
        gAna = evalRes["Gradient"];
        nTotal = Length[xVec];
        {stepC, stepD, stepPhi} = Which[
            ListQ[stepSize] && Length[stepSize] === 3, N[stepSize],
            ListQ[stepSize] && Length[stepSize] === 2, {N[stepSize[[1]]], N[stepSize[[2]]], 10^-6},
            True, {N[stepSize], Automatic, Min[N[stepSize], 10^-6]}
        ];
        stepVec = Table[
            Which[
                i <= nCtrl, stepC,
                i <= nCtrl + nAct,
                    If[NumericQ[stepD], stepD, Max[10^-12, Min[stepC, 10^-5 * Abs[xVec[[i]]]]]],
                True,
                    If[NumericQ[stepPhi], stepPhi, 10^-6]
            ],
            {i, nTotal}
        ];
        gFD = Table[
            xPlus = xVec; xPlus[[i]] += stepVec[[i]];
            xMinus = xVec; xMinus[[i]] -= stepVec[[i]];
            fPlus = EvaluateOptimalControl[comp, xPlus]["Objective"];
            fMinus = EvaluateOptimalControl[comp, xMinus]["Objective"];
            (fPlus - fMinus) / (2. * stepVec[[i]]),
            {i, nTotal}
        ];
        diff = Abs[gAna - gFD];
        absErr = Max[diff];
        ctrlRelErr = Max[diff[[1 ;; nCtrl]]] / Max[Max[Abs[gFD[[1 ;; nCtrl]]]], Max[Abs[gAna[[1 ;; nCtrl]]]], 10^-12];
        If[hasDur,
            durRelErr = Max[diff[[nCtrl + 1 ;; nCtrl + nAct]]] /
                Max[Max[Abs[gFD[[nCtrl + 1 ;; nCtrl + nAct]]]], Max[Abs[gAna[[nCtrl + 1 ;; nCtrl + nAct]]]], 10^-12]
        ];
        If[hasPhi,
            phiRelErr = Max[diff[[nCtrl + nAct + 1 ;; -1]]] /
                Max[Max[Abs[gFD[[nCtrl + nAct + 1 ;; -1]]]], Max[Abs[gAna[[nCtrl + nAct + 1 ;; -1]]]], 10^-12]
        ];
        relErr = Max[ctrlRelErr, durRelErr, phiRelErr];
        outAssoc = <|
            "Objective" -> evalRes["Objective"],
            "Fidelity" -> evalRes["Fidelity"],
            "AnalyticalGradient" -> gAna,
            "FDGradient" -> Developer`ToPackedArray[gFD],
            "MaxAbsoluteError" -> absErr,
            "MaxRelativeError" -> relErr,
            "ControlMaxRelativeError" -> ctrlRelErr,
            "PassedQ" -> TrueQ[relErr <= tol]
        |>;
        If[hasDur, AppendTo[outAssoc, "DurationMaxRelativeError" -> durRelErr]];
        If[hasPhi, AppendTo[outAssoc, "PhaseMaxRelativeError" -> phiRelErr]];
        Return[outAssoc]
    ];

    evalRes = EvaluateOptimalControl[comp, controls];
    If[evalRes === $Failed, Return[$Failed]];
    gAna = evalRes["Gradient"];
    xMat = evalRes["Parameters"];
    {nRows, nCols} = Dimensions[xMat];
    gFD = Table[
        xPlus = xMat; xPlus[[r, col]] += stepSize;
        xMinus = xMat; xMinus[[r, col]] -= stepSize;
        fPlus = EvaluateOptimalControl[comp, xPlus]["Objective"];
        fMinus = EvaluateOptimalControl[comp, xMinus]["Objective"];
        (fPlus - fMinus) / (2. * stepSize),
        {r, nRows}, {col, nCols}
    ];
    diff = Abs[gAna - gFD];
    absErr = Max[diff];
    scaleDenom = Max[Max[Abs[gFD]], Max[Abs[gAna]], 10^-12];
    relErr = absErr / scaleDenom;
    <|
        "Objective" -> evalRes["Objective"],
        "Fidelity" -> evalRes["Fidelity"],
        "AnalyticalGradient" -> gAna,
        "FDGradient" -> Developer`ToPackedArray[gFD],
        "MaxAbsoluteError" -> absErr,
        "MaxRelativeError" -> relErr,
        "PassedQ" -> TrueQ[relErr <= tol]
    |>
];

(* Initial parameter builder *)
resampleControls[oldCtrl_, oldDt_List, newDt_List] := Module[
    {oldEdges, oldMid, tMax, newEdges, newMid, ntOld, ntNew, mCtrl, colVals, interp},
    ntOld = Length[oldDt]; ntNew = Length[newDt]; mCtrl = Dimensions[oldCtrl][[2]];
    If[ntOld === ntNew && Max[Abs[oldDt - newDt]] <= 10^-12 * Max[oldDt], Return[oldCtrl]];
    oldEdges = Prepend[Accumulate[oldDt], 0.];
    oldMid = 0.5 * (Most[oldEdges] + Rest[oldEdges]);
    tMax = Last[oldEdges];
    newEdges = Prepend[Accumulate[newDt], 0.];
    newMid = (0.5 * (Most[newEdges] + Rest[newEdges])) * (tMax / Last[newEdges]);
    Transpose[Table[
        colVals = oldCtrl[[All, k]];
        interp = Quiet[Interpolation[Transpose[{oldMid, colVals}], InterpolationOrder -> 1]];
        Table[interp[Clip[tVal, {First[oldMid], Last[oldMid]}]], {tVal, newMid}],
        {k, mCtrl}
    ]]
];

buildInitialVector[compiled_Association, initSpec_, initScale_?NumericQ, seedSpec_] := Module[
    {nParams, mCtrl, basisAssoc, nt, dtVec, ctrlScales, lbCtrl, ubCtrl, hasCtrlBnd,
     hasDur, durMode, activeIdx, nAct, dLB, dUB, hasFixedTot, actTarSum, dtAct0,
     hasPhi, nPhi, phi0, ctrlInitSpec,
     nCtrl, rawMat, rawVec, ctrlOut, resampled, proxCtx, ortho, qProxCtx, chanNorms, uRand = None, q0Mat, qFeas},
    nParams = compiled["NumParams"];
    mCtrl = compiled["NumControls"];
    basisAssoc = compiled["ControlBasis"];
    nt = compiled["NumSlices"];
    dtVec = compiled["TimeSteps"];
    nCtrl = nParams * mCtrl;
    ctrlScales = Lookup[compiled, "ControlParameterScales", compiled["ParameterScales"]];
    lbCtrl = Lookup[compiled, "ControlLowerBounds", compiled["LowerBounds"]];
    ubCtrl = Lookup[compiled, "ControlUpperBounds", compiled["UpperBounds"]];
    hasCtrlBnd = Lookup[compiled, "HasControlBounds", compiled["HasBounds"]];
    hasDur = TrueQ[Lookup[compiled, "HasDurationOptimization", False]];
    durMode = Lookup[compiled, "DurationMode", If[hasDur, "Slices", "None"]];
    activeIdx = Lookup[compiled, "ActiveDurationIndices", {}];
    nAct = If[hasDur, Lookup[compiled, "NumActiveDurations", 0], 0];
    dLB = Lookup[compiled, "DurationLowerBounds", {}];
    dUB = Lookup[compiled, "DurationUpperBounds", {}];
    hasFixedTot = TrueQ[Lookup[compiled, "FixedTotalDuration", False]];
    actTarSum = Lookup[compiled, "ActiveDurationTargetSum", None];
    hasPhi = TrueQ[Lookup[compiled, "HasPhaseFreedom", False]];
    nPhi = If[hasPhi, Lookup[compiled, "NumPhases", 0], 0];

    proxCtx = Lookup[compiled, "ProxContext", None];
    ortho = If[AssociationQ[proxCtx], Lookup[proxCtx, "OrthoReduction", None], None];
    qProxCtx = If[AssociationQ[proxCtx], Lookup[proxCtx, "QProxContext", None], None];
    chanNorms = Lookup[compiled, "ChannelNorms",
        With[{colMean = Max[Mean[Map[Norm[dtVec * #] &, Transpose[basisAssoc["BasisMatrix"]]]], 10^-12]},
            Table[Max[Mean[ctrlScales[[k ;; ;; mCtrl]]] / colMean, 10^-6], {k, mCtrl}]
        ]
    ];

    {ctrlInitSpec, dtAct0, phi0} = Which[
        (hasDur || hasPhi) && VectorQ[initSpec, NumericQ] && Length[initSpec] === nCtrl + nAct + nPhi,
            {initSpec[[1 ;; nCtrl]],
             If[hasDur, Developer`ToPackedArray[N[Re[initSpec[[nCtrl + 1 ;; nCtrl + nAct]]]]], {}],
             If[hasPhi, Developer`ToPackedArray[N[Re[initSpec[[nCtrl + nAct + 1 ;; -1]]]]], {}]},
        (hasDur || hasPhi) && MatchQ[initSpec, _OptimalControlResult],
            {initSpec,
             If[hasDur,
                 If[durMode === "Total",
                     Developer`ToPackedArray[{N[Re[Lookup[initSpec[[1]], "TotalDuration", Total[initSpec["TimeSteps"]]]]]}],
                     If[VectorQ[initSpec["TimeSteps"], NumericQ] && Length[initSpec["TimeSteps"]] === nt,
                         Developer`ToPackedArray[N[(initSpec["TimeSteps"])[[activeIdx]]]],
                         Developer`ToPackedArray[N[dtVec[[activeIdx]]]]
                     ]
                 ],
                 {}
             ],
             If[hasPhi && AssociationQ[Lookup[initSpec[[1]], "OptimalPhases", None]] &&
                VectorQ[initSpec["OptimalPhases", "Vector"], NumericQ] &&
                Length[initSpec["OptimalPhases", "Vector"]] === nPhi,
                 Developer`ToPackedArray[N[Re[initSpec["OptimalPhases", "Vector"]]]],
                 ConstantArray[0., nPhi]]},
        (hasDur || hasPhi) && AssociationQ[initSpec],
            {Lookup[initSpec, "Parameters", Lookup[initSpec, "Controls", Automatic]],
             If[hasDur,
                 Which[
                     durMode === "Total" && NumericQ[Lookup[initSpec, "TotalDuration", None]],
                         Developer`ToPackedArray[{N[Re[initSpec["TotalDuration"]]]}],
                     VectorQ[Lookup[initSpec, "TimeSteps", None], NumericQ] && Length[initSpec["TimeSteps"]] === nt,
                         If[durMode === "Total",
                             Developer`ToPackedArray[{N[Total[initSpec["TimeSteps"]]]}],
                             Developer`ToPackedArray[N[(initSpec["TimeSteps"])[[activeIdx]]]]
                         ],
                     True,
                         If[durMode === "Total",
                             Developer`ToPackedArray[{N[Total[dtVec]]}],
                             Developer`ToPackedArray[N[dtVec[[activeIdx]]]]
                         ]
                 ],
                 {}
             ],
             If[hasPhi && VectorQ[Lookup[initSpec, "Phases", None], NumericQ] && Length[initSpec["Phases"]] === nPhi,
                 Developer`ToPackedArray[N[Re[initSpec["Phases"]]]],
                 ConstantArray[0., nPhi]]},
        True,
            {initSpec,
             If[hasDur,
                 If[durMode === "Total",
                     Developer`ToPackedArray[{N[Total[dtVec]]}],
                     Developer`ToPackedArray[N[dtVec[[activeIdx]]]]
                 ],
                 {}
             ],
             ConstantArray[0., nPhi]}
    ];

    rawVec = Which[
        ctrlInitSpec === "Zero",
            ConstantArray[0., nCtrl],
        NumericQ[ctrlInitSpec],
            If[basisAssoc["IsIdentityQ"],
                ConstantArray[N[Re[ctrlInitSpec]], nCtrl],
                Flatten[ToPackedRealMatrix[LeastSquares[basisAssoc["BasisMatrix"],
                    ConstantArray[N[Re[ctrlInitSpec]], {nt, mCtrl}]]]]
            ],
        MatchQ[ctrlInitSpec, _OptimalControlResult],
            Which[
                MatrixQ[ctrlInitSpec["Parameters"], NumericQ] && Dimensions[ctrlInitSpec["Parameters"]] === {nParams, mCtrl},
                    Flatten[ToPackedRealMatrix[ctrlInitSpec["Parameters"]]],
                MatrixQ[ctrlInitSpec["Controls"], NumericQ] && Dimensions[ctrlInitSpec["Controls"]][[2]] === mCtrl,
                    resampled = resampleControls[ctrlInitSpec["Controls"], ctrlInitSpec["TimeSteps"], dtVec];
                    rawMat = If[basisAssoc["IsIdentityQ"],
                        resampled,
                        LeastSquares[basisAssoc["BasisMatrix"], resampled]
                    ];
                    Flatten[ToPackedRealMatrix[rawMat]],
                True,
                    $Failed
            ],
        MatrixQ[ctrlInitSpec, NumericQ] && Dimensions[ctrlInitSpec] === {nParams, mCtrl},
            Flatten[ToPackedRealMatrix[ctrlInitSpec]],
        !basisAssoc["IsIdentityQ"] && MatrixQ[ctrlInitSpec, NumericQ] && Dimensions[ctrlInitSpec] === {nt, mCtrl},
            Flatten[ToPackedRealMatrix[LeastSquares[basisAssoc["BasisMatrix"], ToPackedRealMatrix[ctrlInitSpec]]]],
        VectorQ[ctrlInitSpec, NumericQ] && Length[ctrlInitSpec] === nCtrl,
            Developer`ToPackedArray[N[Re[ctrlInitSpec]]],
        ctrlInitSpec === Automatic,
            If[basisAssoc["IsIdentityQ"],
                BlockRandom[
                    If[seedSpec =!= Automatic, SeedRandom[seedSpec]];
                    (initScale / ctrlScales) * RandomReal[{-1., 1.}, nCtrl]
                ],
                BlockRandom[
                    If[seedSpec =!= Automatic, SeedRandom[seedSpec]];
                    uRand = ToPackedRealMatrix[
                        (initScale / (Sqrt[N[nt]] * Outer[Times, dtVec, chanNorms])) * RandomReal[{-1., 1.}, {nt, mCtrl}]
                    ]
                ];
                With[{oUse = If[AssociationQ[ortho], ortho, buildOrthonormalBasisReduction[basisAssoc]]},
                    If[AssociationQ[oUse],
                        Flatten[ToPackedRealMatrix[oUse["Rpinv"] . (oUse["QT"] . uRand)]],
                        Flatten[ToPackedRealMatrix[LeastSquares[basisAssoc["BasisMatrix"], uRand]]]
                    ]
                ]
            ],
        True,
            $Failed
    ];
    If[rawVec === $Failed, Return[$Failed]];
    ctrlOut = If[AssociationQ[proxCtx] && TrueQ[Lookup[proxCtx, "HasAnyConstraint", False] || compiled["HasConstraints"]],
        If[!basisAssoc["IsIdentityQ"] && AssociationQ[qProxCtx] && AssociationQ[ortho],
            q0Mat = If[MatrixQ[uRand, NumericQ],
                ToPackedRealMatrix[ortho["QT"] . uRand],
                ToPackedRealMatrix[ortho["R"] . Partition[rawVec, mCtrl]]
            ];
            qFeas = applyProximal[qProxCtx, q0Mat, 0.];
            Developer`ToPackedArray[Flatten[ToPackedRealMatrix[ortho["Rpinv"] . qFeas]]],
            Developer`ToPackedArray[Flatten[applyProximal[proxCtx, Partition[rawVec, mCtrl], 0.]]]
        ],
        Developer`ToPackedArray[If[hasCtrlBnd, MapThread[Clip[#1, {#2, #3}] &, {rawVec, lbCtrl, ubCtrl}], rawVec]]
    ];
    If[!hasDur && !hasPhi, Return[ctrlOut]];
    If[hasDur,
        dtAct0 = If[hasFixedTot,
            projectBoundedSimplex[dtAct0, dLB, dUB, actTarSum],
            Developer`ToPackedArray[MapThread[Clip[#1, {#2, #3}] &, {dtAct0, dLB, dUB}]]
        ]
    ];
    Developer`ToPackedArray[Join[ctrlOut, dtAct0, phi0]]
];

(* Line search configuration & Strong-Wolfe / Projected Armijo engines *)
resolveLineSearchConfig[lsOpt_, wolfeOpt_, maxLsOpt_, optMethod_String, hasConstraints_?BooleanQ] := Module[
    {lsHead, lsExtra, regLS, isReg, lsMode, wPair, c1, c2, maxLs},
    {lsHead, lsExtra} = Which[
        lsOpt === Automatic || StringQ[lsOpt], {lsOpt, <||>},
        MatchQ[lsOpt, {_String, _Association}] && Length[lsOpt] === 2, {lsOpt[[1]], lsOpt[[2]]},
        MatchQ[lsOpt, {_String, (_Rule | _RuleDelayed) ..}], {lsOpt[[1]], Association[Rest[lsOpt]]},
        True, {$Failed, <||>}
    ];
    regLS = Lookup[$OCSolverRegistry, "LineSearches", <||>];
    isReg = StringQ[lsHead] && KeyExistsQ[regLS, lsHead];
    If[lsHead === $Failed || (!MemberQ[{Automatic, "StrongWolfe", "Armijo"}, lsHead] && !isReg),
        Message[OptimalControlProblem::badls, lsOpt]; Return[$Failed]
    ];
    wPair = Which[
        wolfeOpt === Automatic,
            Which[
                lsHead === "HagerZhang", {0.1, If[optMethod === "ConjugateGradient", 0.4, 0.9]},
                True, {10^-4, If[optMethod === "ConjugateGradient", 0.4, 0.9]}
            ],
        MatchQ[wolfeOpt, {_?NumericQ, _?NumericQ}] && TrueQ[0. < wolfeOpt[[1]] < wolfeOpt[[2]] < 1.] &&
                (lsHead =!= "HagerZhang" || TrueQ[wolfeOpt[[1]] < 0.5]),
            N[wolfeOpt],
        True,
            $Failed
    ];
    If[wPair === $Failed, Message[OptimalControlProblem::badwolfe, wolfeOpt]; Return[$Failed]];
    {c1, c2} = wPair;
    maxLs = If[IntegerQ[maxLsOpt] && maxLsOpt >= 1, maxLsOpt, $Failed];
    If[maxLs === $Failed, Message[OptimalControlProblem::badls, MaxLineSearchEvaluations -> maxLsOpt]; Return[$Failed]];
    lsMode = Which[
        lsHead === "Armijo", "Armijo",
        hasConstraints && MemberQ[{"StrongWolfe", "HagerZhang"}, lsHead], "ProjectedArmijoFallback",
        hasConstraints && lsHead === Automatic, "ProjectedArmijo",
        isReg, lsHead,
        True, "StrongWolfe"
    ];
    <|
        "Mode" -> lsMode,
        "C1" -> c1,
        "C2" -> c2,
        "MaxEvals" -> maxLs,
        "ExtraOpts" -> lsExtra,
        "RequestedLineSearch" -> lsHead
    |>
];

safeguardedZoomStep[pA_Association, pB_Association] := Module[
    {a1 = pA["Alpha"], f1 = pA["Phi"], g1 = pA["DPhi"],
     a2 = pB["Alpha"], f2 = pB["Phi"], g2 = pB["DPhi"],
     lo, hi, span, d1, disc, d2, denomC, aCubic, denomQ, aQuad},
    lo = Min[a1, a2];
    hi = Max[a1, a2];
    span = hi - lo;
    If[Abs[a1 - a2] <= 10^-20, Return[0.5 * (lo + hi)]];
    (* Try cubic Hermite interpolation (Nocedal & Wright eq. 3.59) *)
    d1 = g1 + g2 - 3. * (f1 - f2) / (a1 - a2);
    disc = d1^2 - g1 * g2;
    aCubic = Indeterminate;
    If[TrueQ[disc >= 0.],
        d2 = Sign[a2 - a1] * Sqrt[disc];
        denomC = g2 - g1 + 2. * d2;
        If[TrueQ[Abs[denomC] > 10^-24],
            aCubic = a2 - (a2 - a1) * (g2 + d2 - d1) / denomC
        ]
    ];
    If[NumericQ[aCubic] && TrueQ[lo + 0.1 * span <= aCubic <= hi - 0.1 * span],
        Return[aCubic]
    ];
    (* Quadratic interpolation fallback using (a1, f1, g1) and (a2, f2) *)
    denomQ = 2. * (f2 - f1 - g1 * (a2 - a1));
    aQuad = If[TrueQ[denomQ > 10^-24], a1 - g1 * (a2 - a1)^2 / denomQ, Indeterminate];
    If[NumericQ[aQuad] && TrueQ[lo + 0.1 * span <= aQuad <= hi - 0.1 * span],
        Return[aQuad]
    ];
    0.5 * (lo + hi)
];

strongWolfeZoom[evalPhi_, armijoQ_, swCurvQ_, pLoInit_Association, pHiInit_Association,
                maxLsEvals_Integer, lsEvals_Integer] := Module[
    {pLo = pLoInit, pHi = pHiInit, nEv = lsEvals, alphaJ, aFailUpper = Infinity, pJ, chosen = None, swMet = False},
    While[nEv < maxLsEvals,
        If[Abs[pHi["Alpha"] - pLo["Alpha"]] <= 10^-13 * Max[1., pLo["Alpha"], pHi["Alpha"]], Break[]];
        alphaJ = If[NumericQ[pHi["Phi"]] && NumericQ[pHi["DPhi"]],
            safeguardedZoomStep[pLo, pHi],
            0.5 * (pLo["Alpha"] + Min[pHi["Alpha"], aFailUpper])
        ];
        nEv++;
        pJ = evalPhi[alphaJ];
        If[pJ === $Failed,
            aFailUpper = Min[aFailUpper, alphaJ];
            If[Abs[alphaJ - pLo["Alpha"]] <= 10^-13 * Max[1., pLo["Alpha"]], Break[]];
            alphaJ = 0.5 * (pLo["Alpha"] + alphaJ);
            pHi = <|"Alpha" -> alphaJ, "Phi" -> Indeterminate, "DPhi" -> Indeterminate|>;
            Continue[]
        ];
        If[!armijoQ[pJ["Phi"], alphaJ] || TrueQ[pJ["Phi"] >= pLo["Phi"]],
            pHi = pJ,
            If[swCurvQ[pJ["DPhi"]],
                chosen = pJ; swMet = True; Break[]
            ];
            If[TrueQ[pJ["DPhi"] * (pHi["Alpha"] - pLo["Alpha"]) >= 0.],
                pHi = pLo
            ];
            pLo = pJ
        ]
    ];
    {chosen, swMet, nEv}
];

strongWolfeLineSearch[evalScaled_, th_, fVal_?NumericQ, gVal_, dir_, initStep_?NumericQ,
                      c1_?NumericQ, c2_?NumericQ, maxLsEvals_Integer, tolScale_?NumericQ] := Module[
    {phi0 = -fVal, dphi0, alphaMax = 10000., lsEvals = 0, p0, pPrev, pCur,
     alpha, armijoQ, swCurvQ, bestArmijo = None, chosen = None, swMet = False, evalPhi, zOut},
    dphi0 = -(gVal . dir);
    If[!TrueQ[dphi0 < 0.], Return[$Failed]];
    armijoQ[fTrial_, a_] := With[{tolEff = Min[tolScale, 0.25 * Abs[c1 * a * dphi0]]},
        TrueQ[fTrial <= phi0 + c1 * a * dphi0 + tolEff]
    ];
    swCurvQ[dfTrial_] := TrueQ[Abs[dfTrial] <= -c2 * dphi0];
    evalPhi[a_?NumericQ] := Module[{res},
        res = evalScaled[Developer`ToPackedArray[th + a * dir]];
        If[res === $Failed || !NumericQ[res["Objective"]] || !VectorQ[res["Gradient"], NumericQ],
            Return[$Failed]
        ];
        With[{pt = <|"Alpha" -> a, "Phi" -> -res["Objective"], "DPhi" -> -(res["Gradient"] . dir), "Res" -> res|>},
            If[armijoQ[pt["Phi"], a] && TrueQ[pt["Phi"] < If[bestArmijo === None, phi0, bestArmijo["Phi"]]],
                bestArmijo = pt
            ];
            pt
        ]
    ];

    p0 = <|"Alpha" -> 0., "Phi" -> phi0, "DPhi" -> dphi0|>;
    pPrev = p0;
    alpha = Min[alphaMax, Max[initStep, 10^-12]];

    Do[
        If[lsEvals >= maxLsEvals, Break[]];
        lsEvals++;
        pCur = evalPhi[alpha];
        If[pCur === $Failed,
            alpha = 0.5 * (pPrev["Alpha"] + alpha);
            If[alpha <= 10^-13, Break[]];
            Continue[]
        ];
        If[!armijoQ[pCur["Phi"], alpha] || (lsIter > 1 && TrueQ[pCur["Phi"] >= pPrev["Phi"]]),
            zOut = strongWolfeZoom[evalPhi, armijoQ, swCurvQ, pPrev, pCur, maxLsEvals, lsEvals];
            {chosen, swMet, lsEvals} = zOut;
            Break[]
        ];
        If[swCurvQ[pCur["DPhi"]],
            chosen = pCur; swMet = True; Break[]
        ];
        If[TrueQ[pCur["DPhi"] >= 0.],
            zOut = strongWolfeZoom[evalPhi, armijoQ, swCurvQ, pCur, pPrev, maxLsEvals, lsEvals];
            {chosen, swMet, lsEvals} = zOut;
            Break[]
        ];
        pPrev = pCur;
        If[alpha >= alphaMax, chosen = pCur; Break[]];
        alpha = Min[alphaMax, 2.0 * alpha],
        {lsIter, maxLsEvals}
    ];

    If[chosen === None, chosen = bestArmijo];
    If[chosen === None,
        $Failed,
        <|
            "Accepted" -> True,
            "Theta" -> Developer`ToPackedArray[th + chosen["Alpha"] * dir],
            "Res" -> chosen["Res"],
            "StepLength" -> chosen["Alpha"],
            "StrongWolfeSatisfied" -> swMet,
            "Trials" -> lsEvals
        |>
    ]
];

projectedArmijoLineSearch[evalScaled_, projBox_, th_, fVal_?NumericQ, gVal_, dir_, initStep_?NumericQ,
                          c1_?NumericQ, maxLsEvals_Integer, tolScale_?NumericQ, quadInterp_?BooleanQ] := Module[
    {step = initStep, lsLimit, thTrial, resTrial, fTrial, gs, tolEff, d0, denom, stepQ, trials = 0, out = $Failed},
    lsLimit = Max[maxLsEvals, 36];
    Do[
        trials++;
        thTrial = Developer`ToPackedArray[projBox[th + step * dir]];
        resTrial = evalScaled[thTrial];
        If[resTrial === $Failed || !NumericQ[resTrial["Objective"]] || !VectorQ[resTrial["Gradient"], NumericQ],
            step *= 0.5;
            If[step < 10^-13, Break[]];
            Continue[]
        ];
        fTrial = resTrial["Objective"];
        gs = gVal . (thTrial - th);
        tolEff = If[TrueQ[gs > 0.], Min[tolScale, 0.5 * c1 * gs], 0.];
        If[TrueQ[gs > 0.] && TrueQ[fTrial >= fVal + c1 * gs - tolEff],
            out = <|
                "Accepted" -> True,
                "Theta" -> thTrial,
                "Res" -> resTrial,
                "StepLength" -> step,
                "StrongWolfeSatisfied" -> False,
                "Trials" -> trials
            |>;
            Break[]
        ];
        If[TrueQ[quadInterp],
            d0 = -gs / Max[step, 10^-20];
            denom = 2. * ((-fTrial) - (-fVal) - d0 * step);
            stepQ = If[TrueQ[d0 < 0. && denom > 10^-20], -d0 * step^2 / denom, 0.5 * step];
            step = Clip[stepQ, {0.25 * step, 0.5 * step}],
            step *= 0.5
        ];
        If[step < 10^-13, Break[]],
        {lsLimit}
    ];
    out
];

parseMonitorConfig[monOpt_, monIntOpt_] := Module[{spec, active},
    spec = Which[
        IntegerQ[monIntOpt] && monIntOpt >= 1,
            <|"Kind" -> "Iterations", "Every" -> monIntOpt|>,
        MatchQ[monIntOpt, {"Time", _?NumericQ}] && TrueQ[monIntOpt[[2]] > 0.],
            <|"Kind" -> "Time", "Every" -> N[monIntOpt[[2]]]|>,
        True,
            $Failed
    ];
    If[spec === $Failed,
        Message[OptimalControlProblem::badmonitor, monIntOpt];
        Return[$Failed]
    ];
    active = !MatchQ[monOpt, None | False | Automatic];
    <|"Active" -> active, "Callback" -> monOpt, "Spec" -> spec|>
];

emitMonitorIfDue[monCfg_Association, iter_Integer, evals_Integer, obj_, fid_, bestFid_, gradNorm_,
                 stepLen_, tElapsed_, status_String, extra_Association, stage_,
                 lastMonTime_?NumericQ, lastMonIter_Integer] := Module[{due, payload},
    If[!TrueQ[monCfg["Active"]], Return[{lastMonTime, lastMonIter}]];
    If[iter <= lastMonIter && !(lastMonIter < 0 && iter === 0), Return[{lastMonTime, lastMonIter}]];
    due = (status =!= "Running") ||
        (iter > 0 && monCfg["Spec", "Kind"] === "Iterations" && Mod[iter, monCfg["Spec", "Every"]] === 0) ||
        (iter > 0 && monCfg["Spec", "Kind"] === "Time" && (tElapsed - lastMonTime >= monCfg["Spec", "Every"]));
    If[!due, Return[{lastMonTime, lastMonIter}]];
    payload = <|
        "Iteration" -> iter,
        "Objective" -> obj,
        "Fidelity" -> fid,
        "BestFidelity" -> bestFid,
        "GradientNorm" -> gradNorm,
        "StepLength" -> stepLen,
        "MeanFidelity" -> Lookup[extra, "MeanFidelity", fid],
        "MinFidelity" -> Lookup[extra, "MinFidelity", fid],
        "ElapsedTime" -> tElapsed,
        "Controls" -> Lookup[extra, "Controls", None],
        "Status" -> status,
        "Evaluations" -> evals,
        "Stage" -> stage
    |>;
    Quiet[Check[(monCfg["Callback"])[payload], Null]];
    {tElapsed, iter}
];

evalCriterionFidelity[fid_?NumericQ, ext_, crit_] := Switch[crit,
    "Minimum" | "Min", If[AssociationQ[ext], Lookup[ext, "MinFidelity", fid], fid],
    "Mean", If[AssociationQ[ext], Lookup[ext, "MeanFidelity", fid], fid],
    _, fid
];

(* Dimensionless scaled optimization engine (Projected L-BFGS / BFGS / PR+ CG / Steepest Ascent) *)
runScaledOptimizer[evalFlatFn_, x0Vec_List, scaleVec_List, lbVec_List, ubVec_List, hasBounds_?BooleanQ, optsList_List,
                   customProjFn_: None, hasCustomConstr_: False, evalObjOnlyFn_: None] := Module[
    {method, mMem, maxIter, targetFid, fidCrit, gradTol, funTol, stepTol, tLimit, recHist, verbose,
     lsOpt, wolfeOpt, maxLsOpt, monOpt, monIntOpt, monCfg, stageTag, lastMonTime = 0., lastMonIter = -1, bestFid,
     lsConfig, lsMode, c1Wolfe, c2Wolfe, maxLsEvals, lsExtra, regLS, nmState = Automatic,
     projBox, tStart, dynTime = 0., lbT, ubT, th, fVal, gVal, fidVal, penVal, extra,
     evals = 0, objEvals = 0, gradEvals = 0, swCount = 0, approxWolfeCnt = 0,
     acceptedSteps = 0, rejectedSteps = 0, skippedCurv = 0,
     pg, pgNorm, objHist, fidHist, gHist, sList = {}, yList = {}, rhoList = {}, hInv, nDim,
     nIter = 0, status = "MaxIterationsReached", stopNow, freeMask, gTan, epsTan, qVec, alpha, beta, gamma, dir,
     activePairs, sf, yf, syf, sLast, yLast, syLast, yyLast,
     initStep, lsOut, accepted, thNew, fNew, gNew, fidNew, penNew, extraNew, sVec, yVec, sy,
     oldDir, betaPR, evalScaled, evalObjScaled, tolScale, prevDeltaF = 0., lastStep = 0., fOld, runLineSearchOnce,
     anyConstr = TrueQ[hasBounds || hasCustomConstr]},
    {method, mMem, maxIter, targetFid, gradTol, funTol, stepTol, tLimit, recHist, verbose,
     lsOpt, wolfeOpt, maxLsOpt, monOpt, monIntOpt} =
        {Optimizer, LBFGSMemory, MaxIterations, TargetFidelity, GradientTolerance,
         FunctionTolerance, StepTolerance, TimeConstraint, RecordHistory, Verbose,
         LineSearch, WolfeParameters, MaxLineSearchEvaluations, StepMonitor, MonitorInterval} /. optsList /. Options[GRAPEOptimize];
    fidCrit = Lookup[FilterRules[optsList, {$FidelityCriterionKey}], $FidelityCriterionKey, Automatic];
    stageTag = Lookup[FilterRules[optsList, {$MonitorStage}], $MonitorStage, "Optimization"];
    monCfg = parseMonitorConfig[monOpt, monIntOpt];
    If[monCfg === $Failed, Return[$Failed]];
    If[method === Automatic, method = "LBFGS"];
    If[anyConstr && method === "BFGS", method = "LBFGS"];
    lsConfig = resolveLineSearchConfig[lsOpt, wolfeOpt, maxLsOpt, method, anyConstr];
    If[lsConfig === $Failed, Return[$Failed]];
    {lsMode, c1Wolfe, c2Wolfe, maxLsEvals} = {lsConfig["Mode"], lsConfig["C1"], lsConfig["C2"], lsConfig["MaxEvals"]};
    lsExtra = Lookup[lsConfig, "ExtraOpts", <||>];
    regLS = Lookup[$OCSolverRegistry, "LineSearches", <||>];
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
        evals++;
        objEvals++;
        gradEvals++;
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
    evalObjScaled = If[evalObjOnlyFn === None,
        None,
        Function[thVec,
            Module[{tD0, rawObj},
                evals++;
                objEvals++;
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
            ]
        ]
    ];
    runLineSearchOnce[dVec_, iStep_?NumericQ] := Which[
        KeyExistsQ[regLS, lsMode],
            (regLS[lsMode])[
                <|
                    "EvalScaled" -> evalScaled,
                    "EvalObjScaled" -> evalObjScaled,
                    "ProjBox" -> projBox,
                    "HasConstraints" -> anyConstr,
                    "Theta" -> th,
                    "FVal" -> fVal,
                    "GVal" -> gVal,
                    "Dir" -> dVec,
                    "InitStep" -> iStep,
                    "C1" -> c1Wolfe,
                    "C2" -> c2Wolfe,
                    "MaxEvals" -> maxLsEvals,
                    "TolScale" -> tolScale,
                    "PrevDeltaF" -> prevDeltaF,
                    "ExtraOpts" -> lsExtra,
                    "NMState" -> nmState
                |>
            ],
        lsMode === "StrongWolfe",
            strongWolfeLineSearch[evalScaled, th, fVal, gVal, dVec, iStep, c1Wolfe, c2Wolfe, maxLsEvals, tolScale],
        True,
            projectedArmijoLineSearch[evalScaled, projBox, th, fVal, gVal, dVec, iStep, c1Wolfe, maxLsEvals, tolScale,
                lsMode =!= "Armijo"]
    ];

    th = Developer`ToPackedArray[projBox[x0Vec * scaleVec]];
    With[{res0 = evalScaled[th]},
        If[res0 === $Failed, Return[$Failed]];
        {fVal, fidVal, penVal, gVal, extra} =
            {res0["Objective"], res0["Fidelity"], res0["Penalty"], res0["Gradient"], res0["Extra"]}
    ];
    bestFid = fidVal;
    pg = Developer`ToPackedArray[projBox[th + gVal] - th];
    pgNorm = Norm[pg];
    objHist = If[TrueQ[recHist], {fVal}, {}];
    fidHist = If[TrueQ[recHist], {fidVal}, {}];
    gHist = If[TrueQ[recHist], {pgNorm}, {}];
    If[method === "BFGS", hInv = Developer`ToPackedArray[IdentityMatrix[nDim, WorkingPrecision -> MachinePrecision]]];

    Do[
        If[NumericQ[targetFid] && TrueQ[evalCriterionFidelity[fidVal, extra, fidCrit] >= targetFid], status = "TargetFidelityReached"; Break[]];
        If[TrueQ[pgNorm <= gradTol], status = "GradientToleranceReached"; Break[]];
        If[AbsoluteTime[] - tStart >= tLimit, status = "TimeConstraintReached"; Break[]];

        Which[
            hasBounds && !TrueQ[hasCustomConstr],
                freeMask = Developer`ToPackedArray[N[Boole[MapThread[
                    !((#3 > -Infinity && #1 <= #3 + 10^-11 * Max[1., Abs[#3]] && #2 < 0.) ||
                      (#4 < Infinity && #1 >= #4 - 10^-11 * Max[1., Abs[#4]] && #2 > 0.)) &, {th, gVal, lbT, ubT}]]]];
                gTan = Developer`ToPackedArray[gVal * freeMask];
                qVec = gTan,
            TrueQ[hasCustomConstr],
                freeMask = 1.;
                epsTan = 10^-5 * Max[1., Norm[th]] / Max[Norm[gVal], 10^-12];
                gTan = Developer`ToPackedArray[(projBox[th + epsTan * gVal] - th) / epsTan];
                If[!TrueQ[Norm[gTan] > 10^-14 * Max[1., Norm[gVal]]], gTan = gVal];
                qVec = gTan,
            True,
                freeMask = 1.;
                gTan = gVal;
                qVec = gVal
        ];

        dir = Switch[method,
            "LBFGS",
                If[Length[sList] > 0,
                    If[hasBounds && !TrueQ[hasCustomConstr],
                        activePairs = {};
                        Do[
                            sf = sList[[idx]] * freeMask;
                            yf = yList[[idx]] * freeMask;
                            syf = sf . yf;
                            If[TrueQ[syf > 10^-12 * Norm[sf] * Norm[yf]],
                                AppendTo[activePairs, {sf, yf, 1. / syf}]
                            ],
                            {idx, Length[sList]}
                        ];
                        If[Length[activePairs] > 0,
                            alpha = ConstantArray[0., Length[activePairs]];
                            Do[
                                alpha[[idx]] = activePairs[[idx, 3]] * (activePairs[[idx, 1]] . qVec);
                                qVec = qVec - alpha[[idx]] * activePairs[[idx, 2]],
                                {idx, Length[activePairs], 1, -1}
                            ];
                            sLast = Last[activePairs][[1]];
                            yLast = Last[activePairs][[2]];
                            syLast = sLast . yLast;
                            yyLast = yLast . yLast;
                            gamma = If[yyLast > 10^-24 && syLast > 0., syLast / yyLast, 1.];
                            dir = gamma * qVec;
                            Do[
                                beta = activePairs[[idx, 3]] * (activePairs[[idx, 2]] . dir);
                                dir = dir + activePairs[[idx, 1]] * (alpha[[idx]] - beta),
                                {idx, 1, Length[activePairs]}
                            ];
                            dir * freeMask,
                            gTan
                        ],
                        alpha = ConstantArray[0., Length[sList]];
                        Do[
                            alpha[[idx]] = rhoList[[idx]] * (sList[[idx]] . qVec);
                            qVec = qVec - alpha[[idx]] * yList[[idx]],
                            {idx, Length[sList], 1, -1}
                        ];
                        syLast = Last[sList] . Last[yList];
                        yyLast = Last[yList] . Last[yList];
                        gamma = If[yyLast > 10^-24 && syLast > 0., syLast / yyLast, 1.];
                        dir = gamma * qVec;
                        Do[
                            beta = rhoList[[idx]] * (yList[[idx]] . dir);
                            dir = dir + sList[[idx]] * (alpha[[idx]] - beta),
                            {idx, 1, Length[sList]}
                        ];
                        dir
                    ],
                    gTan
                ],
            "BFGS",
                Developer`ToPackedArray[hInv . gVal],
            "ConjugateGradient",
                If[nIter > 0 && Mod[nIter, Max[nDim, 20]] =!= 0,
                    betaPR = Max[0., (gTan . (gTan - oldDir[[2]])) / Max[oldDir[[2]] . oldDir[[2]], 10^-24]];
                    dir = (gTan + betaPR * oldDir[[1]]) * freeMask,
                    dir = gTan
                ],
            _,
                gTan
        ];
        If[method === "ConjugateGradient", oldDir = {dir, gTan}];
        If[!TrueQ[dir . gVal > 10^-14 * Norm[dir] * Norm[gVal]],
            dir = gTan;
            sList = {}; yList = {}; rhoList = {};
            If[method === "BFGS", hInv = Developer`ToPackedArray[IdentityMatrix[nDim, WorkingPrecision -> MachinePrecision]]]
        ];

        initStep = Which[
            MemberQ[{"StrongWolfe", "HagerZhang", "ZhangHagerNonmonotone", "GrippoNonmonotone"}, lsMode] &&
                MemberQ[{"ConjugateGradient", "SteepestAscent"}, method] && nIter > 0 && prevDeltaF > 0.,
                Clip[2.02 * prevDeltaF / Max[gVal . dir, 10^-12], {10^-4, 2.0}],
            Length[sList] === 0 && ((method === "LBFGS") || (method === "BFGS" && nIter === 0)),
                Min[1., 0.35 / Max[Norm[dir], 10^-12]],
            True,
                1.0
        ];
        tolScale = 10^-14 * Max[1., Abs[fVal]];
        lsOut = runLineSearchOnce[dir, initStep];
        accepted = AssociationQ[lsOut] && TrueQ[lsOut["Accepted"]];
        If[AssociationQ[lsOut] && KeyExistsQ[lsOut, "NMState"], nmState = lsOut["NMState"]];
        rejectedSteps += If[AssociationQ[lsOut], Max[0, Lookup[lsOut, "Trials", If[accepted, 1, 0]] - If[accepted, 1, 0]], 1];

        If[!accepted && Norm[dir - gTan] > 10^-12,
            (* Retry once cleanly along projected steepest ascent before terminating *)
            dir = gTan;
            sList = {}; yList = {}; rhoList = {};
            If[method === "BFGS", hInv = Developer`ToPackedArray[IdentityMatrix[nDim, WorkingPrecision -> MachinePrecision]]];
            initStep = Min[1., 0.25 / Max[Norm[dir], 10^-12]];
            lsOut = runLineSearchOnce[dir, initStep];
            accepted = AssociationQ[lsOut] && TrueQ[lsOut["Accepted"]];
            If[AssociationQ[lsOut] && KeyExistsQ[lsOut, "NMState"], nmState = lsOut["NMState"]];
            rejectedSteps += If[AssociationQ[lsOut], Max[0, Lookup[lsOut, "Trials", If[accepted, 1, 0]] - If[accepted, 1, 0]], 1]
        ];

        If[!accepted,
            status = If[pgNorm <= Sqrt[gradTol], "GradientToleranceReached", "LineSearchTerminated"];
            Break[]
        ];

        acceptedSteps++;
        thNew = lsOut["Theta"];
        With[{resAcc = lsOut["Res"]},
            {fNew, fidNew, penNew, gNew, extraNew} =
                {resAcc["Objective"], resAcc["Fidelity"], resAcc["Penalty"], resAcc["Gradient"], resAcc["Extra"]}
        ];
        lastStep = lsOut["StepLength"];
        If[TrueQ[lsOut["StrongWolfeSatisfied"]], swCount++];
        If[TrueQ[lsOut["ApproxWolfeSatisfied"]], approxWolfeCnt++];

        sVec = Developer`ToPackedArray[thNew - th];
        yVec = Developer`ToPackedArray[gVal - gNew]; (* minimize -f => y = (-gNew) - (-gVal) = gVal - gNew *)
        sy = sVec . yVec;
        If[TrueQ[sy > 10^-12 * Norm[sVec] * Norm[yVec]],
            If[method === "LBFGS",
                If[Length[sList] >= mMem,
                    sList = Rest[sList]; yList = Rest[yList]; rhoList = Rest[rhoList]
                ];
                AppendTo[sList, sVec]; AppendTo[yList, yVec]; AppendTo[rhoList, 1. / sy]
            ];
            If[method === "BFGS",
                With[{rho = 1. / sy, hy = hInv . yVec},
                    With[{yhy = yVec . hy},
                        hInv = Developer`ToPackedArray[
                            hInv + rho * (1. + rho * yhy) * Outer[Times, sVec, sVec] -
                            rho * (Outer[Times, hy, sVec] + Outer[Times, sVec, hy])
                        ]
                    ]
                ]
            ],
            If[MemberQ[{"LBFGS", "BFGS"}, method], skippedCurv++]
        ];

        fOld = fVal;
        prevDeltaF = Max[fNew - fVal, 0.];
        nIter++;
        th = thNew; fVal = fNew; fidVal = fidNew; penVal = penNew; gVal = gNew; extra = extraNew;
        bestFid = Max[bestFid, fidVal];
        pg = Developer`ToPackedArray[projBox[th + gVal] - th];
        pgNorm = Norm[pg];
        If[TrueQ[verbose] && (Mod[nIter, 10] === 0 || nIter === 1),
            Print["Iter ", nIter, ": Fidelity=", NumberForm[fidNew, 6],
                " Objective=", NumberForm[fNew, 6], " step=", NumberForm[lastStep, 4],
                " |grad|=", ScientificForm[pgNorm, 3]]
        ];
        If[TrueQ[recHist],
            AppendTo[objHist, fVal];
            AppendTo[fidHist, fidVal];
            AppendTo[gHist, pgNorm]
        ];

        With[{critFidNew = evalCriterionFidelity[fidVal, extra, fidCrit]},
            stopNow = Which[
                TrueQ[Abs[fNew - fOld] <= funTol * Max[1., Abs[fOld]]] && TrueQ[pgNorm <= Sqrt[gradTol]],
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
        {lastMonTime, lastMonIter} = emitMonitorIfDue[
            monCfg, nIter, evals, fVal, fidVal, bestFid, pgNorm, lastStep,
            AbsoluteTime[] - tStart, If[stopNow === None, "Running", stopNow],
            extra, stageTag, lastMonTime, lastMonIter
        ];
        If[stopNow =!= None, status = stopNow; Break[]],
        {maxIter}
    ];
    If[lastMonIter < nIter || (lastMonIter < 0 && nIter === 0),
        {lastMonTime, lastMonIter} = emitMonitorIfDue[
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
            "Optimizer" -> method,
            "Iterations" -> nIter,
            "Evaluations" -> evals,
            "ObjectiveEvaluations" -> objEvals,
            "GradientEvaluations" -> gradEvals,
            "HVPEvaluations" -> 0,
            "HessianEvaluations" -> 0,
            "AcceptedSteps" -> acceptedSteps,
            "RejectedSteps" -> rejectedSteps,
            "SkippedCurvaturePairs" -> skippedCurv,
            "LineSearchMethod" -> lsMode,
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

(* Proximal Nesterov (FISTA) with Lipschitz backtracking, monotone nonconvex safeguard & adaptive restart *)
runFISTAOptimizer[evalFlatFn_, x0Vec_List, comp_Association, optsList_List] := Module[
    {maxIter, targetFid, fidCrit, gradTol, funTol, stepTol, tLimit, recHist, verbose, restartOpt,
     monOpt, monIntOpt, monCfg, stageTag, lastMonTime = 0., lastMonIter = -1, bestFid, stopNow,
     mCtrl, nParams, nCtrl, proxCtx, hasDur, nAct, dLB, dUB, hasFixedTot, actTarSum,
     hasPhi, nPhi, hasExtra, phiRatio,
     ctrlScales, scaleVec, avgScale, dtRatio, xInitVec, proxState, lMin, lMax, tStart, evals = 0,
     objEvals = 0, gradEvals = 0, restarts = 0, evalUnscaled, xMat, yMat, resX, fX, gX, gradX,
     fidX, penX, objX, extraX, fY, gradY, lVal, lTry, lSafe, dG, epsP, probePt, probeDiff, resProbe, lEst,
     tK = 1., tNext, objHist, fidHist, gHist, nIter = 0, status = "MaxIterationsReached",
     zMat, resZ, fZ, gZ, dz, tolScale, majOK, xNew, resNew, fNewVal, gNewVal, usedSafe,
     vMat, resV, fV, gV, dv, gradRestart, funRestart, doRestart, yNext, resY,
     pgStep, pgNorm, lastStep = 0., smoothPenNew, extraNew, fidNew, objNew, objOld, dStepMat, dStepNorm},

    {maxIter, targetFid, gradTol, funTol, stepTol, tLimit, recHist, verbose, restartOpt, monOpt, monIntOpt} =
        {MaxIterations, TargetFidelity, GradientTolerance, FunctionTolerance, StepTolerance,
         TimeConstraint, RecordHistory, Verbose, FISTARestart, StepMonitor, MonitorInterval} /. optsList /. Options[GRAPEOptimize];
    fidCrit = Lookup[FilterRules[optsList, {$FidelityCriterionKey}], $FidelityCriterionKey, Automatic];
    stageTag = Lookup[FilterRules[optsList, {$MonitorStage}], $MonitorStage, "Optimization"];
    monCfg = parseMonitorConfig[monOpt, monIntOpt];
    If[monCfg === $Failed, Return[$Failed]];
    If[!MemberQ[{"Gradient", "Function", "Both", None, False}, restartOpt],
        Message[OptimalControlProblem::badrestart, restartOpt]; Return[$Failed]
    ];

    nParams = comp["NumParams"];
    mCtrl = comp["NumControls"];
    nCtrl = nParams * mCtrl;
    proxCtx = comp["ProxContext"];
    hasDur = TrueQ[Lookup[comp, "HasDurationOptimization", False]];
    nAct = If[hasDur, Lookup[comp, "NumActiveDurations", 0], 0];
    hasPhi = TrueQ[Lookup[comp, "HasPhaseFreedom", False]];
    nPhi = If[hasPhi, Lookup[comp, "NumPhases", 0], 0];
    hasExtra = hasDur || hasPhi;
    ctrlScales = Lookup[comp, "ControlParameterScales", comp["ParameterScales"]];
    avgScale = Max[Mean[ctrlScales], 10^-12];
    dtRatio = If[hasDur, Max[comp["DurationScale"], 10^-12] / avgScale, 1.];
    phiRatio = 1. / avgScale;
    If[hasDur,
        dLB = comp["DurationLowerBounds"];
        dUB = comp["DurationUpperBounds"];
        hasFixedTot = TrueQ[comp["FixedTotalDuration"]];
        actTarSum = comp["ActiveDurationTargetSum"]
    ];
    scaleVec = Developer`ToPackedArray[Join[
        ctrlScales,
        If[hasDur, ConstantArray[avgScale, nAct], {}],
        If[hasPhi, ConstantArray[avgScale, nPhi], {}]
    ]];
    lMin = 10^-14 * avgScale^2;
    lMax = 10^14 * avgScale^2;
    tStart = AbsoluteTime[];

    proxState[stateIn_, eta_] := If[!hasExtra,
        applyProximal[proxCtx, stateIn, eta],
        Module[{cPart, cOut, dtPart, dtOut, phiOut},
            cPart = Partition[stateIn[[1 ;; nCtrl]], mCtrl];
            cOut = Flatten[applyProximal[proxCtx, cPart, eta]];
            dtOut = If[hasDur,
                dtPart = stateIn[[nCtrl + 1 ;; nCtrl + nAct]] / dtRatio;
                If[hasFixedTot,
                    projectBoundedSimplex[dtPart, dLB, dUB, actTarSum],
                    MapThread[Clip[#1, {#2, #3}] &, {dtPart, dLB, dUB}]
                ] * dtRatio,
                {}
            ];
            phiOut = If[hasPhi, stateIn[[nCtrl + nAct + 1 ;; -1]], {}];
            Developer`ToPackedArray[Join[cOut, dtOut, phiOut]]
        ]
    ];

    evalUnscaled[stateIn_] := Module[{rawEval, vecPhys, gOut},
        evals++;
        objEvals++;
        gradEvals++;
        vecPhys = If[!hasExtra,
            Developer`ToPackedArray[Flatten[stateIn]],
            Developer`ToPackedArray[Join[
                stateIn[[1 ;; nCtrl]],
                If[hasDur, stateIn[[nCtrl + 1 ;; nCtrl + nAct]] / dtRatio, {}],
                If[hasPhi, stateIn[[nCtrl + nAct + 1 ;; -1]] / phiRatio, {}]
            ]]
        ];
        rawEval = evalFlatFn[vecPhys];
        If[rawEval === $Failed || !NumericQ[rawEval["Fidelity"]], Return[$Failed]];
        If[hasExtra,
            gOut = Developer`ToPackedArray[Join[
                rawEval["Gradient"][[1 ;; nCtrl]],
                If[hasDur, rawEval["Gradient"][[nCtrl + 1 ;; nCtrl + nAct]] / dtRatio, {}],
                If[hasPhi, rawEval["Gradient"][[nCtrl + nAct + 1 ;; -1]] / phiRatio, {}]
            ]];
            Join[rawEval, <|"Gradient" -> gOut|>],
            rawEval
        ]
    ];

    xInitVec = If[!hasExtra,
        Partition[Developer`ToPackedArray[N[Re[x0Vec]]], mCtrl],
        Developer`ToPackedArray[Join[
            N[Re[x0Vec[[1 ;; nCtrl]]]],
            If[hasDur, N[Re[x0Vec[[nCtrl + 1 ;; nCtrl + nAct]]]] * dtRatio, {}],
            If[hasPhi, N[Re[x0Vec[[nCtrl + nAct + 1 ;; -1]]]] * phiRatio, {}]
        ]]
    ];
    xMat = proxState[xInitVec, 0.];
    resX = evalUnscaled[xMat];
    If[resX === $Failed, Return[$Failed]];
    fidX = resX["Fidelity"];
    bestFid = fidX;
    penX = resX["Penalty"];
    objX = resX["Objective"];
    fX = -Lookup[resX, "SmoothObjective", fidX - penX];
    gX = Lookup[resX, "NonSmoothPenalty", evalNonSmoothPenalty[proxCtx, resX["Controls"]]];
    gradX = Developer`ToPackedArray[-resX["Gradient"]];
    extraX = KeyDrop[resX, {"Objective", "SmoothObjective", "Fidelity", "Penalty", "SmoothPenalty", "NonSmoothPenalty", "Gradient"}];

    yMat = xMat; fY = fX; gradY = gradX;
    dG = Norm[Flatten[gradX]];
    lVal = Clip[avgScale^2, {lMin, lMax}];
    If[dG > 10^-14,
        epsP = (10^-3 / avgScale) * (gradX / dG);
        probePt = If[hasExtra, proxState[xMat + epsP, 0.], xMat + epsP];
        probeDiff = Norm[Flatten[probePt - xMat]];
        If[TrueQ[probeDiff > 10^-14],
            resProbe = evalUnscaled[probePt];
            If[AssociationQ[resProbe],
                lEst = Norm[Flatten[(-resProbe["Gradient"]) - gradX]] / probeDiff;
                lVal = Clip[Max[lEst, 0.25 * avgScale^2], {lMin, lMax}]
            ]
        ]
    ];

    pgStep = proxState[xMat - (1. / lVal) * gradX, 1. / lVal];
    pgNorm = Norm[Flatten[(xMat - pgStep) * lVal] / scaleVec];
    objHist = If[TrueQ[recHist], {objX}, {}];
    fidHist = If[TrueQ[recHist], {fidX}, {}];
    gHist = If[TrueQ[recHist], {pgNorm}, {}];

    Do[
        If[NumericQ[targetFid] && TrueQ[evalCriterionFidelity[fidX, extraX, fidCrit] >= targetFid], status = "TargetFidelityReached"; Break[]];
        If[TrueQ[pgNorm <= gradTol], status = "GradientToleranceReached"; Break[]];
        If[AbsoluteTime[] - tStart >= tLimit, status = "TimeConstraintReached"; Break[]];

        lTry = If[nIter === 0, lVal, Max[lMin, 0.9 * lVal]];
        tolScale = 10^-14 * Max[1., Abs[fY], Abs[objX]];
        majOK = False;
        resZ = $Failed;
        Do[
            zMat = proxState[yMat - (1. / lTry) * gradY, 1. / lTry];
            resZ = evalUnscaled[zMat];
            If[resZ =!= $Failed,
                fZ = -Lookup[resZ, "SmoothObjective", resZ["Fidelity"] - resZ["Penalty"]];
                gZ = Lookup[resZ, "NonSmoothPenalty", evalNonSmoothPenalty[proxCtx, resZ["Controls"]]];
                dz = zMat - yMat;
                If[TrueQ[fZ <= fY + Total[Flatten[gradY * dz]] + 0.5 * lTry * Total[Flatten[dz * dz]] + tolScale] ||
                   Norm[Flatten[dz]] <= 10^-15 * Max[1., Norm[Flatten[yMat]]],
                    lVal = lTry; majOK = True; Break[]
                ]
            ];
            If[lTry >= lMax, Break[]];
            lTry = Min[lMax, 2.0 * lTry],
            {45}
        ];

        If[!majOK && resZ === $Failed,
            status = If[pgNorm <= Sqrt[gradTol], "GradientToleranceReached", "LineSearchTerminated"];
            Break[]
        ];

        (* Nonconvex monotone safeguard: never let total objective phi(x) + g(x) degrade across iterations *)
        usedSafe = False;
        If[TrueQ[fZ + gZ <= fX + gX + tolScale],
            xNew = zMat; resNew = resZ; fNewVal = fZ; gNewVal = gZ,
            usedSafe = True;
            lSafe = Max[lVal, lTry];
            xNew = xMat; resNew = resX; fNewVal = fX; gNewVal = gX;
            Do[
                vMat = proxState[xMat - (1. / lSafe) * gradX, 1. / lSafe];
                resV = evalUnscaled[vMat];
                If[resV =!= $Failed,
                    fV = -Lookup[resV, "SmoothObjective", resV["Fidelity"] - resV["Penalty"]];
                    gV = Lookup[resV, "NonSmoothPenalty", evalNonSmoothPenalty[proxCtx, resV["Controls"]]];
                    dv = vMat - xMat;
                    If[TrueQ[fV <= fX + Total[Flatten[gradX * dv]] + 0.5 * lSafe * Total[Flatten[dv * dv]] + tolScale] ||
                       Norm[Flatten[dv]] <= 10^-15 * Max[1., Norm[Flatten[xMat]]],
                        If[TrueQ[fV + gV <= fNewVal + gNewVal + tolScale],
                            xNew = vMat; resNew = resV; fNewVal = fV; gNewVal = gV
                        ];
                        lVal = Max[lVal, lSafe];
                        Break[]
                    ]
                ];
                If[lSafe >= lMax, Break[]];
                lSafe = Min[lMax, 2.0 * lSafe],
                {40}
            ];
            If[TrueQ[fZ + gZ < fNewVal + gNewVal],
                xNew = zMat; resNew = resZ; fNewVal = fZ; gNewVal = gZ
            ]
        ];

        lastStep = 1. / lVal;
        fidNew = resNew["Fidelity"];
        objNew = resNew["Objective"];
        smoothPenNew = Lookup[resNew, "SmoothPenalty", resNew["Penalty"]];
        extraNew = KeyDrop[resNew, {"Objective", "SmoothObjective", "Fidelity", "Penalty", "SmoothPenalty", "NonSmoothPenalty", "Gradient"}];

        gradRestart = TrueQ[Total[Flatten[(yMat - xNew) * (xNew - xMat)]] > 0.];
        funRestart = usedSafe || TrueQ[fNewVal + gNewVal > fX + gX + tolScale];
        doRestart = Which[
            usedSafe, True,
            restartOpt === "Gradient", gradRestart,
            restartOpt === "Function", funRestart,
            restartOpt === "Both", gradRestart || funRestart,
            True, False
        ];

        pgStep = proxState[xNew - (1. / lVal) * Developer`ToPackedArray[-resNew["Gradient"]], 1. / lVal];
        pgNorm = Norm[Flatten[(xNew - pgStep) * lVal] / scaleVec];
        objOld = objX;
        dStepMat = xNew - xMat;
        dStepNorm = Norm[Flatten[dStepMat]];
        nIter++;
        If[TrueQ[verbose] && (Mod[nIter, 10] === 0 || nIter === 1),
            Print["Iter ", nIter, " (FISTA): Fidelity=", NumberForm[fidNew, 6],
                " Objective=", NumberForm[objNew, 6], " L=", ScientificForm[lVal, 3],
                " |pg|=", ScientificForm[pgNorm, 3]]
        ];

        xMat = xNew; resX = resNew; fX = fNewVal; gX = gNewVal; objX = objNew; fidX = fidNew;
        bestFid = Max[bestFid, fidX];
        penX = resNew["Penalty"]; gradX = Developer`ToPackedArray[-resNew["Gradient"]]; extraX = extraNew;
        If[TrueQ[recHist],
            AppendTo[objHist, objX];
            AppendTo[fidHist, fidX];
            AppendTo[gHist, pgNorm]
        ];

        With[{critFidX = evalCriterionFidelity[fidX, extraX, fidCrit]},
            stopNow = Which[
                TrueQ[Abs[objNew - objOld] <= funTol * Max[1., Abs[objOld]]] && TrueQ[pgNorm <= Sqrt[gradTol]],
                    If[NumericQ[targetFid] && TrueQ[critFidX >= targetFid], "TargetFidelityReached", "FunctionToleranceReached"],
                TrueQ[dStepNorm <= stepTol * Max[1., Norm[Flatten[xNew]]]],
                    If[NumericQ[targetFid] && TrueQ[critFidX >= targetFid], "TargetFidelityReached", "StepToleranceReached"],
                NumericQ[targetFid] && TrueQ[critFidX >= targetFid],
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
        {lastMonTime, lastMonIter} = emitMonitorIfDue[
            monCfg, nIter, evals, objX, fidX, bestFid, pgNorm, lastStep,
            AbsoluteTime[] - tStart, If[stopNow === None, "Running", stopNow],
            extraX, stageTag, lastMonTime, lastMonIter
        ];
        If[stopNow =!= None, status = stopNow; Break[]];

        If[doRestart,
            tNext = 1.;
            yNext = xMat;
            If[tK > 1., restarts++];
            fY = fX;
            gradY = gradX,
            tNext = 0.5 * (1. + Sqrt[1. + 4. * tK^2]);
            yNext = xMat + ((tK - 1.) / tNext) * dStepMat;
            resY = evalUnscaled[yNext];
            If[resY === $Failed,
                yNext = xMat; tNext = 1.; fY = fX; gradY = gradX,
                fY = -Lookup[resY, "SmoothObjective", resY["Fidelity"] - resY["Penalty"]];
                gradY = Developer`ToPackedArray[-resY["Gradient"]]
            ]
        ];
        yMat = yNext; tK = tNext,
        {maxIter}
    ];
    If[lastMonIter < nIter || (lastMonIter < 0 && nIter === 0),
        {lastMonTime, lastMonIter} = emitMonitorIfDue[
            monCfg, nIter, evals, objX, fidX, bestFid, pgNorm, lastStep,
            AbsoluteTime[] - tStart, status, extraX, stageTag, lastMonTime, lastMonIter
        ]
    ];

    <|
        "ParametersVec" -> If[!hasExtra,
            Developer`ToPackedArray[Flatten[xMat]],
            Developer`ToPackedArray[Join[
                xMat[[1 ;; nCtrl]],
                If[hasDur, xMat[[nCtrl + 1 ;; nCtrl + nAct]] / dtRatio, {}],
                If[hasPhi, xMat[[nCtrl + nAct + 1 ;; -1]] / phiRatio, {}]
            ]]
        ],
        "Objective" -> objX,
        "Fidelity" -> fidX,
        "Penalty" -> penX,
        "GradientNorm" -> pgNorm,
        "ScaledGradient" -> Developer`ToPackedArray[Flatten[gradX] / scaleVec],
        "Optimizer" -> "FISTA",
        "Iterations" -> nIter,
        "Evaluations" -> evals,
        "ObjectiveEvaluations" -> objEvals,
        "GradientEvaluations" -> gradEvals,
        "LineSearchMethod" -> "LipschitzBacktracking",
        "StrongWolfeSatisfiedCount" -> 0,
        "LastStepLength" -> lastStep,
        "FISTARestarts" -> restarts,
        "LipschitzEstimate" -> lVal,
        "Status" -> status,
        "ObjectiveHistory" -> Developer`ToPackedArray[objHist],
        "FidelityHistory" -> Developer`ToPackedArray[fidHist],
        "GradientNormHistory" -> Developer`ToPackedArray[gHist],
        "ElapsedTime" -> AbsoluteTime[] - tStart,
        "Extra" -> extraX
    |>
];

(* Unified optimization dispatcher across smooth LBFGS/BFGS/CG/SteepestAscent, FISTA, and registered optimizers *)
dispatchOptimizationEngine[evalFlatFn_, x0Vec_List, comp_Association, optsList_List, evalObjOnlyFn_: None] := Module[
    {optOpt, optHead, optExtra, regOpts, isRegOpt, hasNonSmoothPen, hasSlew, isIdent, hasBox, hasDisk, hasParamConstr,
     hasDur, hasFixedTot, activeIdx, nAct, dLB, dUB, actTarSum, dtScaleS, nCtrl,
     hasPhi, nPhi, hasExtra,
     ctrlScales, lbCtrl, ubCtrl, mustFISTA, useCtrlProx, proxCtx, ortho, qProxCtx, mCtrl,
     optScales, customProjFn, effOpts, dtVec, chanNorms, qColScales, qScales, qFullScales,
     compQ, nQCtrl, q0Vec, evalQFlatFn, evalQObjOnlyFn, invokeSmoothOrReg, optRes, pOptMat},
    optOpt = Optimizer /. optsList /. Options[GRAPEOptimize];
    If[optOpt === "Automatic", optOpt = Automatic];
    {optHead, optExtra} = Which[
        optOpt === Automatic || StringQ[optOpt], {optOpt, <||>},
        MatchQ[optOpt, {_String, _Association}] && Length[optOpt] === 2, {optOpt[[1]], optOpt[[2]]},
        MatchQ[optOpt, {_String, (_Rule | _RuleDelayed) ..}], {optOpt[[1]], Association[Rest[optOpt]]},
        True, {$Failed, <||>}
    ];
    regOpts = Lookup[$OCSolverRegistry, "Optimizers", <||>];
    isRegOpt = StringQ[optHead] && KeyExistsQ[regOpts, optHead];
    If[optHead === $Failed || (!MemberQ[{Automatic, "LBFGS", "BFGS", "ConjugateGradient", "SteepestAscent", "FISTA"}, optHead] && !isRegOpt),
        Message[OptimalControlProblem::badopt, optOpt]; Return[$Failed]
    ];
    hasNonSmoothPen = TrueQ[comp["L1Penalty"] > 0. || comp["TVPenalty"] > 0.];
    hasSlew = TrueQ[comp["HasSlewBounds"]];
    isIdent = comp["ControlBasis", "IsIdentityQ"];
    hasBox = TrueQ[comp["HasBounds"]];
    hasDisk = TrueQ[comp["HasDiskBounds"]];
    hasDur = TrueQ[Lookup[comp, "HasDurationOptimization", False]];
    hasFixedTot = TrueQ[Lookup[comp, "FixedTotalDuration", False]];
    activeIdx = Lookup[comp, "ActiveDurationIndices", {}];
    nAct = If[hasDur, Lookup[comp, "NumActiveDurations", 0], 0];
    dLB = Lookup[comp, "DurationLowerBounds", {}];
    dUB = Lookup[comp, "DurationUpperBounds", {}];
    actTarSum = Lookup[comp, "ActiveDurationTargetSum", None];
    dtScaleS = Lookup[comp, "DurationScale", 1.];
    hasPhi = TrueQ[Lookup[comp, "HasPhaseFreedom", False]];
    nPhi = If[hasPhi, Lookup[comp, "NumPhases", 0], 0];
    hasExtra = hasDur || hasPhi;
    nCtrl = comp["NumParams"] * comp["NumControls"];
    ctrlScales = Lookup[comp, "ControlParameterScales", comp["ParameterScales"]];
    lbCtrl = Lookup[comp, "ControlLowerBounds", comp["LowerBounds"]];
    ubCtrl = Lookup[comp, "ControlUpperBounds", comp["UpperBounds"]];

    proxCtx = Lookup[comp, "ProxContext", <||>];
    hasParamConstr = (!isIdent) && TrueQ[Lookup[proxCtx, "HasAnyConstraint", False] || (hasDisk || hasSlew || Lookup[comp, "HasControlBounds", False])];
    mustFISTA = hasNonSmoothPen || hasSlew;

    If[mustFISTA && (MemberQ[{"LBFGS", "BFGS", "ConjugateGradient", "SteepestAscent"}, optHead] || isRegOpt),
        Message[OptimalControlProblem::nonsmoothopt, optHead];
        optHead = "FISTA";
        optExtra = <||>;
        isRegOpt = False
    ];
    If[optHead === Automatic,
        optHead = If[mustFISTA || hasParamConstr, "FISTA", "LBFGS"]
    ];
    effOpts = Prepend[FilterRules[optsList, Except[Optimizer]], Optimizer -> optHead];

    invokeSmoothOrReg[evalFnIn_, evalObjIn_, x0In_, scalesIn_, lbIn_, ubIn_, hasBIn_, projIn_, hasCustIn_, compIn_] := If[
        isRegOpt,
        (regOpts[optHead])[<|
            "EvalFn" -> evalFnIn,
            "EvalObjFn" -> evalObjIn,
            "X0Vec" -> x0In,
            "ScaleVec" -> scalesIn,
            "LBVec" -> lbIn,
            "UBVec" -> ubIn,
            "HasBounds" -> hasBIn,
            "CustomProjFn" -> projIn,
            "HasCustomConstr" -> hasCustIn,
            "AnyConstraints" -> TrueQ[hasBIn || hasCustIn],
            "OptsList" -> effOpts,
            "SubOpts" -> optExtra,
            "Comp" -> compIn
        |>],
        runScaledOptimizer[evalFnIn, x0In, scalesIn, lbIn, ubIn, hasBIn, effOpts, projIn, hasCustIn, evalObjIn]
    ];

    ortho = Lookup[proxCtx, "OrthoReduction", None];
    qProxCtx = Lookup[proxCtx, "QProxContext", None];
    mCtrl = comp["NumControls"];

    (* Non-identity parameterized bases with control constraints/penalties run in thin-SVD orthonormal Q-space *)
    If[!isIdent && AssociationQ[qProxCtx] && AssociationQ[ortho],
        dtVec = comp["TimeSteps"];
        chanNorms = Lookup[comp, "ChannelNorms",
            With[{colMean = Max[Mean[Map[Norm[dtVec * #] &, Transpose[comp["ControlBasis", "BasisMatrix"]]]], 10^-12]},
                Table[Max[Mean[ctrlScales[[k ;; ;; mCtrl]]] / colMean, 10^-6], {k, mCtrl}]
            ]
        ];
        qColScales = Map[Max[Norm[dtVec * #], 10^-12] &, Transpose[ortho["Q"]]];
        qScales = Developer`ToPackedArray[Flatten[Outer[Times, qColScales, chanNorms]]];
        nQCtrl = ortho["Rank"] * mCtrl;
        qFullScales = Developer`ToPackedArray[Join[
            qScales,
            If[hasDur, ConstantArray[dtScaleS, nAct], {}],
            If[hasPhi, ConstantArray[1., nPhi], {}]
        ]];
        compQ = Join[comp, <|
            "ControlBasis" -> ortho["QBasis"],
            "NumParams" -> ortho["Rank"],
            "ProxContext" -> qProxCtx,
            "ControlParameterScales" -> qScales,
            "ParameterScales" -> qFullScales
        |>];
        evalQFlatFn = With[{rPinv = ortho["Rpinv"], qT = ortho["QT"], mC = mCtrl, eFn = evalFlatFn,
                            nQC = nQCtrl, nC = nCtrl, hdQ = hasDur, hpQ = hasPhi, hxQ = hasExtra,
                            nA = nAct, aIdx = activeIdx},
            Function[qVecIn,
                Module[{qPart, qProjMat, dtPart, phiPart, vecIn, qProjRes, gQFlat},
                    If[hxQ,
                        qPart = qVecIn[[1 ;; nQC]];
                        dtPart = If[hdQ, qVecIn[[nQC + 1 ;; nQC + nA]], {}];
                        phiPart = If[hpQ, qVecIn[[nQC + nA + 1 ;; -1]], {}];
                        qProjMat = ToPackedRealMatrix[rPinv . Partition[qPart, mC]];
                        vecIn = Developer`ToPackedArray[Join[Flatten[qProjMat], dtPart, phiPart]];
                        qProjRes = eFn[vecIn];
                        If[qProjRes === $Failed, Return[$Failed]];
                        gQFlat = Join[
                            Flatten[qT . qProjRes["SliceGradient"]],
                            If[hdQ, Lookup[qProjRes, "ActiveDurationGradient", (qProjRes["DurationGradient"])[[aIdx]]], {}],
                            If[hpQ, qProjRes["PhaseGradient"], {}]
                        ];
                        Join[qProjRes, <|"Gradient" -> Developer`ToPackedArray[gQFlat]|>],
                        qProjMat = ToPackedRealMatrix[rPinv . Partition[qVecIn, mC]];
                        qProjRes = eFn[Developer`ToPackedArray[Flatten[qProjMat]]];
                        If[qProjRes === $Failed, Return[$Failed]];
                        Join[qProjRes, <|"Gradient" -> ToPackedRealMatrix[qT . qProjRes["SliceGradient"]]|>]
                    ]
                ]
            ]
        ];
        evalQObjOnlyFn = If[evalObjOnlyFn === None,
            None,
            With[{rPinv = ortho["Rpinv"], mC = mCtrl, eObj = evalObjOnlyFn,
                  nQC = nQCtrl, hdQ = hasDur, hpQ = hasPhi, hxQ = hasExtra, nA = nAct},
                Function[qVecIn,
                    If[hxQ,
                        eObj[Developer`ToPackedArray[Join[
                            Flatten[ToPackedRealMatrix[rPinv . Partition[qVecIn[[1 ;; nQC]], mC]]],
                            If[hdQ, qVecIn[[nQC + 1 ;; nQC + nA]], {}],
                            If[hpQ, qVecIn[[nQC + nA + 1 ;; -1]], {}]
                        ]]],
                        eObj[Developer`ToPackedArray[Flatten[ToPackedRealMatrix[rPinv . Partition[qVecIn, mC]]]]]
                    ]
                ]
            ]
        ];
        q0Vec = If[hasExtra,
            Developer`ToPackedArray[Join[
                Flatten[ToPackedRealMatrix[ortho["R"] . Partition[x0Vec[[1 ;; nCtrl]], mCtrl]]],
                x0Vec[[nCtrl + 1 ;; -1]]
            ]],
            Developer`ToPackedArray[Flatten[ToPackedRealMatrix[ortho["R"] . Partition[x0Vec, mCtrl]]]]
        ];
        optRes = If[optHead === "FISTA",
            runFISTAOptimizer[evalQFlatFn, q0Vec, compQ, effOpts],
            With[{uS = Max[Mean[qScales], 10^-12]},
                optScales = Developer`ToPackedArray[Join[
                    ConstantArray[uS, nQCtrl],
                    If[hasDur, ConstantArray[dtScaleS, nAct], {}],
                    If[hasPhi, ConstantArray[1., nPhi], {}]
                ]];
                customProjFn = If[hasExtra,
                    With[{ctx = qProxCtx, mC = mCtrl, sC = uS, sD = dtScaleS, nQC = nQCtrl, nA = nAct,
                          hdQ = hasDur, hpQ = hasPhi,
                          minD = If[hasDur, dLB * dtScaleS, {}], maxD = If[hasDur, dUB * dtScaleS, {}],
                          fxQ = hasFixedTot, tarD = If[hasFixedTot, actTarSum * dtScaleS, 0.]},
                        Function[thIn,
                            Module[{cP, dP, phiP},
                                cP = Flatten[applyProximal[ctx, Partition[thIn[[1 ;; nQC]] / sC, mC], 0.]] * sC;
                                dP = If[hdQ,
                                    If[fxQ,
                                        projectBoundedSimplex[thIn[[nQC + 1 ;; nQC + nA]], minD, maxD, tarD],
                                        MapThread[Clip[#1, {#2, #3}] &, {thIn[[nQC + 1 ;; nQC + nA]], minD, maxD}]
                                    ],
                                    {}
                                ];
                                phiP = If[hpQ, thIn[[nQC + nA + 1 ;; -1]], {}];
                                Developer`ToPackedArray[Join[cP, dP, phiP]]
                            ]
                        ]
                    ],
                    With[{ctx = qProxCtx, mC = mCtrl, sC = uS},
                        (Developer`ToPackedArray[Flatten[applyProximal[ctx, Partition[# / sC, mC], 0.]] * sC] &)
                    ]
                ];
                invokeSmoothOrReg[evalQFlatFn, evalQObjOnlyFn, q0Vec, optScales, {}, {}, True, customProjFn, True, compQ]
            ]
        ];
        If[optRes === $Failed, Return[$Failed]];
        If[hasExtra,
            pOptMat = ToPackedRealMatrix[ortho["Rpinv"] . Partition[(optRes["ParametersVec"])[[1 ;; nQCtrl]], mCtrl]];
            Return[Join[optRes, <|"ParametersVec" -> Developer`ToPackedArray[Join[Flatten[pOptMat], (optRes["ParametersVec"])[[nQCtrl + 1 ;; -1]]]]|>]],
            pOptMat = ToPackedRealMatrix[ortho["Rpinv"] . Partition[optRes["ParametersVec"], mCtrl]];
            Return[Join[optRes, <|"ParametersVec" -> Developer`ToPackedArray[Flatten[pOptMat]]|>]]
        ]
    ];

    If[optHead === "FISTA",
        Return[runFISTAOptimizer[evalFlatFn, x0Vec, comp, effOpts]]
    ];

    useCtrlProx = hasDisk || hasSlew || hasParamConstr;
    If[useCtrlProx || (hasDur && hasFixedTot),
        optScales = If[useCtrlProx,
            With[{uS = Max[Mean[ctrlScales], 10^-12]},
                Developer`ToPackedArray[Join[
                    ConstantArray[uS, nCtrl],
                    If[hasDur, ConstantArray[dtScaleS, nAct], {}],
                    If[hasPhi, ConstantArray[1., nPhi], {}]
                ]]
            ],
            comp["ParameterScales"]
        ];
        customProjFn = If[useCtrlProx,
            If[hasExtra,
                With[{ctx = proxCtx, mC = mCtrl, sC = optScales[[1]], nC = nCtrl, nA = nAct,
                      hdQ = hasDur, hpQ = hasPhi,
                      minD = If[hasDur, dLB * dtScaleS, {}], maxD = If[hasDur, dUB * dtScaleS, {}],
                      fxQ = hasFixedTot, tarD = If[hasFixedTot, actTarSum * dtScaleS, 0.]},
                    Function[thIn,
                        Module[{cP, dP, phiP},
                            cP = Flatten[applyProximal[ctx, Partition[thIn[[1 ;; nC]] / sC, mC], 0.]] * sC;
                            dP = If[hdQ,
                                If[fxQ,
                                    projectBoundedSimplex[thIn[[nC + 1 ;; nC + nA]], minD, maxD, tarD],
                                    MapThread[Clip[#1, {#2, #3}] &, {thIn[[nC + 1 ;; nC + nA]], minD, maxD}]
                                ],
                                {}
                            ];
                            phiP = If[hpQ, thIn[[nC + nA + 1 ;; -1]], {}];
                            Developer`ToPackedArray[Join[cP, dP, phiP]]
                        ]
                    ]
                ],
                With[{ctx = proxCtx, mC = mCtrl, sC = optScales[[1]]},
                    (Developer`ToPackedArray[Flatten[applyProximal[ctx, Partition[# / sC, mC], 0.]] * sC] &)
                ]
            ],
            (* Pure FixedTotalDuration with optional Cartesian control box and optional unconstrained phases *)
            With[{nC = nCtrl, nA = nAct, hpQ = hasPhi, hasB = Lookup[comp, "HasControlBounds", False],
                  minC = lbCtrl * ctrlScales, maxC = ubCtrl * ctrlScales,
                  minD = dLB * dtScaleS, maxD = dUB * dtScaleS, tarD = actTarSum * dtScaleS},
                Function[thIn,
                    Module[{cP, dP, phiP},
                        cP = If[hasB, MapThread[Clip[#1, {#2, #3}] &, {thIn[[1 ;; nC]], minC, maxC}], thIn[[1 ;; nC]]];
                        dP = projectBoundedSimplex[thIn[[nC + 1 ;; nC + nA]], minD, maxD, tarD];
                        phiP = If[hpQ, thIn[[nC + nA + 1 ;; -1]], {}];
                        Developer`ToPackedArray[Join[cP, dP, phiP]]
                    ]
                ]
            ]
        ];
        invokeSmoothOrReg[evalFlatFn, evalObjOnlyFn, x0Vec, optScales, comp["LowerBounds"], comp["UpperBounds"],
            True, customProjFn, True, comp],
        invokeSmoothOrReg[evalFlatFn, evalObjOnlyFn, x0Vec, comp["ParameterScales"], comp["LowerBounds"], comp["UpperBounds"],
            comp["HasBounds"], None, False, comp]
    ]
];

(* P1C: Schedule & MultiStart validators and orchestrators *)
parseTimeStepSchedule[tsOpt_] := Which[
    tsOpt === None || tsOpt === False, {},
    ListQ[tsOpt] && Length[tsOpt] >= 1 && AllTrue[tsOpt, IntegerQ[#] && # >= 1 &], tsOpt,
    True, Message[OptimalControlProblem::badschedule, tsOpt]; $Failed
];

parseMultiStartSpec[msOpt_] := Which[
    msOpt === 1 || msOpt === None || msOpt === False, <|"Mode" -> "Single"|>,
    IntegerQ[msOpt] && msOpt >= 2, <|"Mode" -> "Count", "Count" -> msOpt|>,
    ListQ[msOpt] && Length[msOpt] >= 1, <|"Mode" -> "ExplicitList", "Starts" -> msOpt|>,
    True, Message[OptimalControlProblem::badmultistart, msOpt]; $Failed
];

extractOptimalPhasesAssoc[pfSpec_, phiVec_] := If[AssociationQ[pfSpec] && TrueQ[pfSpec["ParametricQ"]] && VectorQ[phiVec, NumericQ],
    With[{nL = pfSpec["NumLeft"], nR = pfSpec["NumRight"]},
        <|
            "Mode" -> pfSpec["Mode"],
            "Left" -> If[nL > 0, phiVec[[1 ;; nL]], {}],
            "Right" -> If[nR > 0, phiVec[[nL + 1 ;; nL + nR]], {}],
            "Vector" -> Developer`ToPackedArray[N[Re[phiVec]]],
            "LeftLabels" -> pfSpec["LeftLabels"],
            "RightLabels" -> pfSpec["RightLabels"]
        |>
    ],
    None
];

assembleSingleResult[comp_Association, optRes_Association, extraFields_Association: <||>] := Module[
    {finalEval, maxViol, optPhases, leakVal},
    finalEval = EvaluateOptimalControl[comp, optRes["ParametersVec"], True];
    If[finalEval === $Failed, Return[$Failed]];
    maxViol = computeConstraintViolation[comp["ProxContext"], finalEval["Controls"]];
    optPhases = extractOptimalPhasesAssoc[Lookup[comp, "PhaseFreedomSpec", None], Lookup[finalEval, "Phases", None]];
    leakVal = With[{ext = Lookup[finalEval, "Extra", None]},
        If[AssociationQ[ext], Lookup[ext, "Leakage", None], None]
    ];
    OptimalControlResult[Join[<|
        "Kind" -> "Single",
        "Controls" -> finalEval["Controls"],
        "Parameters" -> finalEval["Parameters"],
        "FinalFidelity" -> finalEval["Fidelity"],
        "TerminalFidelity" -> Lookup[finalEval, "TerminalFidelity", finalEval["Fidelity"]],
        "RunningCost" -> Lookup[finalEval, "RunningCost", 0.],
        "Leakage" -> leakVal,
        "OptimalPhases" -> optPhases,
        "FinalObjective" -> finalEval["Objective"],
        "FinalPenalty" -> finalEval["Penalty"],
        "SmoothPenalty" -> finalEval["SmoothPenalty"],
        "NonSmoothPenalty" -> finalEval["NonSmoothPenalty"],
        "MaxConstraintViolation" -> maxViol,
        "Optimizer" -> Lookup[optRes, "Optimizer", "LBFGS"],
        "Gradient" -> finalEval["ControlGradient"],
        "ControlGradient" -> finalEval["ControlGradient"],
        "DurationGradient" -> Lookup[finalEval, "DurationGradient", None],
        "ActiveDurationGradient" -> Lookup[finalEval, "ActiveDurationGradient", None],
        "TotalDurationGradient" -> Lookup[finalEval, "TotalDurationGradient", None],
        "DurationMode" -> Lookup[comp, "DurationMode", "None"],
        "PhaseGradient" -> Lookup[finalEval, "PhaseGradient", None],
        "GradientNorm" -> Norm[Flatten[finalEval["Gradient"]]],
        "ScaledGradientNorm" -> optRes["GradientNorm"],
        "Iterations" -> optRes["Iterations"],
        "Evaluations" -> optRes["Evaluations"],
        "ObjectiveEvaluations" -> optRes["ObjectiveEvaluations"],
        "GradientEvaluations" -> optRes["GradientEvaluations"],
        "LineSearchMethod" -> optRes["LineSearchMethod"],
        "StrongWolfeSatisfiedCount" -> optRes["StrongWolfeSatisfiedCount"],
        "LastStepLength" -> optRes["LastStepLength"],
        "ElapsedTime" -> optRes["ElapsedTime"],
        "Status" -> optRes["Status"],
        "ConvergenceHistory" -> <|
            "Objective" -> optRes["ObjectiveHistory"],
            "Fidelity" -> optRes["FidelityHistory"],
            "GradientNorm" -> optRes["GradientNormHistory"]
        |>,
        "FinalPropagator" -> If[KeyExistsQ[finalEval, "FinalPropagator"], finalEval["FinalPropagator"], None],
        "FinalState" -> If[KeyExistsQ[finalEval, "FinalState"], finalEval["FinalState"], None],
        "Extra" -> Join[Lookup[finalEval, "Extra", <||>], KeyDrop[finalEval, {"Gradient", "ControlGradient", "SliceGradient", "Controls", "Parameters", "FinalPropagator", "Propagators", "States", "FinalState", "SmoothObjective", "SmoothPenalty", "NonSmoothPenalty", "Extra"}]],
        "CompiledProblem" -> comp,
        "TimeSteps" -> Lookup[finalEval, "TimeSteps", comp["TimeSteps"]],
        "TotalDuration" -> Lookup[finalEval, "TotalDuration", comp["TotalDuration"]],
        "InitialTimeSteps" -> comp["TimeSteps"],
        "OptimizeDurations" -> Lookup[comp, "OptimizeDurations", False],
        "FullParameters" -> optRes["ParametersVec"],
        "InitialTimePoint" -> Lookup[comp, "InitialTimePoint", 0.],
        "DriftIntegration" -> Lookup[comp, "DriftIntegration", "Static"],
        "IsTimeDependentDrift" -> Lookup[comp, "IsTimeDependentDrift", False],
        "Space" -> comp["Space"],
        "Method" -> comp["Method"],
        "Basis" -> comp["Basis"],
        "OperatorBasis" -> comp["OperatorBasis"],
        "SymbolicDrift" -> comp["SymbolicDrift"],
        "SymbolicControls" -> comp["SymbolicControls"],
        "DriftUnitFactor" -> comp["DriftUnitFactor"],
        "ControlUnitFactor" -> comp["ControlUnitFactor"],
        "FrequencyUnits" -> comp["FrequencyUnits"],
        "ControlUnits" -> comp["ControlUnits"]
    |>, KeyTake[optRes, {
        "FISTARestarts", "LipschitzEstimate", "DynamicsTime", "OptimizerTime",
        "AcceptedSteps", "RejectedSteps", "HVPEvaluations", "HessianEvaluations",
        "SkippedCurvaturePairs", "ApproxWolfeCount", "ExperimentalConstrainedQ",
        "CompactLBFGSFinalSigma",
        "CGIterations", "NegativeCurvatureCount", "HybridSwitchIteration",
        "HybridSwitchReason", "TrustRegionFinalRadius", "ARCFinalSigma"
    }], extraFields]]
];

optimizeSingleCompiledStart[comp_Association, optsList_List, initOverride_: None, seedOverride_: None] := Module[
    {initOpt, initScale, seedOpt, x0Vec, optRes},
    {initOpt, initScale, seedOpt} = {InitialControls, InitialScale, RandomSeed} /. optsList;
    If[initOverride =!= None, initOpt = initOverride];
    If[seedOverride =!= None, seedOpt = seedOverride];
    x0Vec = buildInitialVector[comp, initOpt, initScale, seedOpt];
    If[x0Vec === $Failed, Return[$Failed]];
    optRes = dispatchOptimizationEngine[
        Function[v, EvaluateOptimalControl[comp, v, False]],
        x0Vec, comp, optsList,
        Function[v, evaluateOptimalControlObjective[comp, v]]
    ];
    If[optRes === $Failed, Return[$Failed]];
    assembleSingleResult[comp, optRes]
];

optimizeSingleWithSchedule[prob:OptimalControlProblem[___], tsList_List, mergedOpts_List] := Module[
    {tAllStart = AbsoluteTime[], defaultNt, dt0, totalDur, maxIterOpt, tLimit,
     s, ntK, stageProb, compK, remTime, stageOpts, stageRes, prevRes = None,
     contHist = {}, totIters = 0, totEvals = 0, totObjEvals = 0, totGradEvals = 0},
    defaultNt = TimeSteps /. mergedOpts;
    dt0 = resolveTimingVector[prob[[4]], defaultNt];
    If[dt0 === $Failed, Return[$Failed]];
    totalDur = Total[dt0];
    maxIterOpt = MaxIterations /. mergedOpts;
    tLimit = TimeConstraint /. mergedOpts;

    Do[
        ntK = tsList[[s]];
        remTime = If[NumericQ[tLimit], Max[0.01, tLimit - (AbsoluteTime[] - tAllStart)], Infinity];
        stageProb = OptimalControlProblem[prob[[1]], prob[[2]], prob[[3]], {totalDur, ntK},
            Sequence @@ FilterRules[Prepend[FilterRules[mergedOpts, Except[TimeSteps]], TimeSteps -> ntK],
                Options[OptimalControlProblem]]
        ];
        compK = CompileOptimalControlProblem[stageProb];
        If[compK === $Failed, Return[$Failed, Module]];
        stageOpts = Join[
            If[s > 1, {InitialControls -> prevRes}, {}],
            If[s < Length[tsList], {MaxIterations -> Max[15, Ceiling[0.5 * maxIterOpt]]}, {}],
            {TimeConstraint -> remTime,
             $MonitorStage -> <|"Kind" -> "TimeStep", "StageIndex" -> s, "TotalStages" -> Length[tsList], "NumSlices" -> ntK|>},
            mergedOpts
        ];
        stageRes = optimizeSingleCompiledStart[compK, stageOpts];
        If[stageRes === $Failed, Return[$Failed, Module]];
        prevRes = stageRes;
        totIters += stageRes["Iterations"];
        totEvals += stageRes["Evaluations"];
        totObjEvals += stageRes["ObjectiveEvaluations"];
        totGradEvals += stageRes["GradientEvaluations"];
        AppendTo[contHist, <|
            "Stage" -> s,
            "Kind" -> "TimeStep",
            "NumSlices" -> ntK,
            "Fidelity" -> stageRes["FinalFidelity"],
            "Objective" -> stageRes["FinalObjective"],
            "Iterations" -> stageRes["Iterations"],
            "Evaluations" -> stageRes["Evaluations"],
            "ElapsedTime" -> stageRes["ElapsedTime"],
            "Status" -> stageRes["Status"]
        |>],
        {s, Length[tsList]}
    ];

    OptimalControlResult[Join[prevRes[[1]], <|
        "Iterations" -> totIters,
        "Evaluations" -> totEvals,
        "ObjectiveEvaluations" -> totObjEvals,
        "GradientEvaluations" -> totGradEvals,
        "ElapsedTime" -> AbsoluteTime[] - tAllStart,
        "ContinuationHistory" -> contHist
    |>]]
];

runMultiStartOrchestrator[singleRunnerFn_, msSpec_Association, optsList_List] := Module[
    {tBatchStart = AbsoluteTime[], initOpt, seedOpt, monOpt, monIntOpt, nStarts,
     startInits, startSeeds, sIdx, resK, summaryList = {}, successList = {},
     failCount = 0, totalEvals = 0, bestIdx, winner, evalAdd},
    {initOpt, seedOpt, monOpt, monIntOpt} =
        {InitialControls, RandomSeed, StepMonitor, MonitorInterval} /. optsList /. Options[GRAPEOptimize];
    If[parseMonitorConfig[monOpt, monIntOpt] === $Failed, Return[$Failed]];

    If[msSpec["Mode"] === "Count",
        nStarts = msSpec["Count"];
        startSeeds = BlockRandom[
            If[seedOpt =!= Automatic, SeedRandom[seedOpt]];
            Table[RandomInteger[{1, 2^30 - 1}], {nStarts}]
        ];
        startInits = ConstantArray[initOpt, nStarts],
        nStarts = Length[msSpec["Starts"]];
        startInits = msSpec["Starts"];
        startSeeds = ConstantArray[seedOpt, nStarts]
    ];

    Do[
        resK = Quiet[Check[
            singleRunnerFn[startInits[[sIdx]], startSeeds[[sIdx]], sIdx, nStarts],
            $Failed
        ], {OptimalControlProblem::badbounds}];
        If[MatchQ[resK, _OptimalControlResult] && NumericQ[resK["FinalObjective"]] && NumericQ[resK["FinalFidelity"]],
            evalAdd = Lookup[resK[[1]], "Evaluations", 0];
            totalEvals += evalAdd;
            AppendTo[successList, {sIdx, resK}];
            AppendTo[summaryList, <|
                "StartIndex" -> sIdx,
                "SuccessQ" -> True,
                "RandomSeed" -> startSeeds[[sIdx]],
                "FinalObjective" -> resK["FinalObjective"],
                "FinalFidelity" -> resK["FinalFidelity"],
                "Iterations" -> resK["Iterations"],
                "Evaluations" -> evalAdd,
                "ElapsedTime" -> resK["ElapsedTime"],
                "Status" -> resK["Status"]
            |>],
            failCount++;
            AppendTo[summaryList, <|
                "StartIndex" -> sIdx,
                "SuccessQ" -> False,
                "RandomSeed" -> startSeeds[[sIdx]],
                "FinalObjective" -> Indeterminate,
                "FinalFidelity" -> Indeterminate,
                "Iterations" -> 0,
                "Evaluations" -> 0,
                "ElapsedTime" -> 0.,
                "Status" -> "Failed"
            |>]
        ],
        {sIdx, nStarts}
    ];

    If[successList === {}, Return[$Failed]];
    bestIdx = First[Ordering[
        Map[{#[[2]]["FinalObjective"], #[[2]]["FinalFidelity"]} &, successList],
        -1
    ]];
    winner = successList[[bestIdx, 2]];
    OptimalControlResult[Join[winner[[1]], <|
        "MultiStartResults" -> summaryList,
        "BestStartIndex" -> successList[[bestIdx, 1]],
        "FailedStarts" -> failCount,
        "TotalEvaluations" -> totalEvals,
        "TotalElapsedTime" -> AbsoluteTime[] - tBatchStart
    |>]]
];

(* Main optimization driver for single OptimalControlProblem *)
GRAPEOptimize[prob:OptimalControlProblem[___], opts___Rule] := Module[
    {mergedOpts, durPenOpt, explicitTFQ, tsOpt, msOpt, monOpt, monIntOpt, tsList, msSpec, comp, singleRunner},
    mergedOpts = Join[{opts}, List @@ (prob[[5 ;;]]), Options[GRAPEOptimize]];
    durPenOpt = DurationPenalty /. mergedOpts;
    explicitTFQ = MemberQ[Keys[Join[{opts}, List @@ (prob[[5 ;;]])]], TargetFidelity];
    If[!explicitTFQ && (NumericQ[durPenOpt] && TrueQ[durPenOpt > 0.] || MatchQ[durPenOpt, {_?NumericQ, _?NumericQ}] && TrueQ[durPenOpt[[1]] > 0.]),
        mergedOpts = Prepend[mergedOpts, TargetFidelity -> None]
    ];
    {tsOpt, msOpt, monOpt, monIntOpt} = {TimeStepSchedule, MultiStart, StepMonitor, MonitorInterval} /. mergedOpts;
    If[parseMonitorConfig[monOpt, monIntOpt] === $Failed, Return[$Failed]];
    tsList = parseTimeStepSchedule[tsOpt];
    msSpec = parseMultiStartSpec[msOpt];
    If[tsList === $Failed || msSpec === $Failed, Return[$Failed]];

    If[msSpec["Mode"] === "Single",
        If[tsList === {},
            comp = CompileOptimalControlProblem[
                OptimalControlProblem[prob[[1]], prob[[2]], prob[[3]], prob[[4]],
                    Sequence @@ FilterRules[mergedOpts, Options[OptimalControlProblem]]]
            ];
            If[comp === $Failed, Return[$Failed]];
            Return[optimizeSingleCompiledStart[comp, mergedOpts]],
            Return[optimizeSingleWithSchedule[prob, tsList, mergedOpts]]
        ]
    ];

    (* MultiStart mode: precompile once when no TimeStepSchedule *)
    If[tsList === {},
        comp = CompileOptimalControlProblem[
            OptimalControlProblem[prob[[1]], prob[[2]], prob[[3]], prob[[4]],
                Sequence @@ FilterRules[mergedOpts, Options[OptimalControlProblem]]]
        ];
        If[comp === $Failed, Return[$Failed]];
        singleRunner = Function[{sInit, sSeed, sIdx, nTot},
            optimizeSingleCompiledStart[comp,
                Prepend[FilterRules[mergedOpts, Except[MultiStart | InitialControls | RandomSeed]],
                    $MonitorStage -> <|"Kind" -> "MultiStart", "StartIndex" -> sIdx, "TotalStarts" -> nTot|>
                ],
                sInit, sSeed
            ]
        ],
        singleRunner = Function[{sInit, sSeed, sIdx, nTot},
            optimizeSingleWithSchedule[prob, tsList,
                Join[{InitialControls -> sInit, RandomSeed -> sSeed,
                      $MonitorStage -> <|"Kind" -> "MultiStart", "StartIndex" -> sIdx, "TotalStarts" -> nTot|>},
                    FilterRules[mergedOpts, Except[MultiStart | InitialControls | RandomSeed]]]
            ]
        ]
    ];
    runMultiStartOrchestrator[singleRunner, msSpec, mergedOpts]
];

GRAPEOptimize[
    drift:Except[_Rule | _OptimalControlProblem | _OptimalControlEnsemble | _CooperativeObjective],
    controls:Except[_Rule], target_, timing:Except[_Rule], opts___Rule] :=
    If[KeyExistsQ[$OCSolverRegistry, "OptimizeCooperativeOrEnsemble"] &&
       MatchQ[target, _CooperativeObjective | _List] && !FreeQ[target, _OptimalControlProblem],
        ($OCSolverRegistry["OptimizeCooperativeOrEnsemble"])[drift, controls, target, timing, {opts}],
        GRAPEOptimize[OptimalControlProblem[drift, controls, target, timing,
            Sequence @@ FilterRules[{opts}, Options[OptimalControlProblem]]], opts]
    ];

OptimalControlOptimize[args___] := GRAPEOptimize[args];

(* Accessor & clean compact SummaryBox formatting for OptimalControlResult *)
OptimalControlResult[assoc_Association][keys__] := assoc[keys];

MakeBoxes[res:OptimalControlResult[assoc_Association], fmt_] :=
    BoxForm`ArrangeSummaryBox[
        OptimalControlResult,
        res,
        None,
        {
            BoxForm`SummaryItem[{"Fidelity: ", NumberForm[assoc["FinalFidelity"], 6]}],
            BoxForm`SummaryItem[{"Status: ", assoc["Status"]}],
            BoxForm`SummaryItem[{"Iterations: ", assoc["Iterations"], " (", assoc["Evaluations"], " evals)"}],
            BoxForm`SummaryItem[{"Duration: ", assoc["TotalDuration"], " s (", Length[assoc["TimeSteps"]], " slices)"}]
        },
        {
            BoxForm`SummaryItem[{"Space / Method: ", assoc["Space"], " / ", assoc["Method"]}],
            BoxForm`SummaryItem[{"Elapsed: ", NumberForm[assoc["ElapsedTime"], 4], " s"}],
            BoxForm`SummaryItem[{"Scaled |grad|: ", ScientificForm[assoc["ScaledGradientNorm"], 3]}]
        },
        fmt
    ] /; AssociationQ[assoc];

(* SpinDynamica Integration Converters *)
OptimalControlWaveform[OptimalControlResult[assoc_Association]] := Module[{dtVec, t0, tEdges},
    dtVec = assoc["TimeSteps"];
    t0 = Lookup[assoc, "InitialTimePoint", 0.];
    tEdges = t0 + Prepend[Accumulate[dtVec], 0.];
    <|
        "Controls" -> assoc["Controls"],
        "Parameters" -> assoc["Parameters"],
        "TimeSteps" -> dtVec,
        "TimeEdges" -> tEdges,
        "TimeMidpoints" -> 0.5 * (Most[tEdges] + Rest[tEdges]),
        "Duration" -> assoc["TotalDuration"],
        "InitialTimePoint" -> t0,
        "DriftIntegration" -> Lookup[assoc, "DriftIntegration", "Static"],
        "Drift" -> assoc["SymbolicDrift"],
        "ChannelOperators" -> assoc["SymbolicControls"],
        "DriftUnitFactor" -> assoc["DriftUnitFactor"],
        "ControlUnitFactor" -> assoc["ControlUnitFactor"],
        "FrequencyUnits" -> assoc["FrequencyUnits"],
        "ControlUnits" -> assoc["ControlUnits"],
        "Space" -> assoc["Space"],
        "Basis" -> assoc["Basis"],
        "OperatorBasis" -> assoc["OperatorBasis"]
    |>
];

toSpaceGenExpr[expr_, "Hilbert", basis_, opbasis_] := Which[
    MatchQ[expr, _NullOperator | _NullSuperoperator] || (NumericQ[expr] && TrueQ[expr == 0]), NullOperator[],
    MatrixQ[expr], Operator[expr, basis, CheckOperatorType -> False],
    True, expr
];

toSpaceGenExpr[expr_, "Liouville", basis_, opbasis_] := Which[
    MatchQ[expr, _NullOperator | _NullSuperoperator] || (NumericQ[expr] && TrueQ[expr == 0]), NullSuperoperator[],
    SuperoperatorQ[expr], expr,
    OperatorQ[expr], CommutationSuperoperator[-I expr],
    MatrixQ[expr], Superoperator[expr, opbasis, basis, CheckSuperoperatorType -> False],
    True, expr
];

buildSliceGeneratorExpr[drift_, ctrls_List, hamFac_, ctrlFac_, vals_List, includeDrift_, space_: "Hilbert"] := Module[{acc},
    acc = If[TrueQ[includeDrift] && !(NumericQ[drift] && TrueQ[drift == 0]) && !MatchQ[drift, _NullOperator | _NullSuperoperator],
        hamFac * drift,
        0
    ];
    Do[
        If[TrueQ[vals[[k]] != 0] && !(NumericQ[ctrls[[k]]] && TrueQ[ctrls[[k]] == 0]) && !MatchQ[ctrls[[k]], _NullOperator | _NullSuperoperator],
            acc += (ctrlFac * vals[[k]]) * ctrls[[k]]
        ],
        {k, Length[ctrls]}
    ];
    If[NumericQ[acc] && TrueQ[acc == 0], If[space === "Liouville", NullSuperoperator[], NullOperator[]], acc]
];

OptimalControlEventList[OptimalControlResult[assoc_Association], opts___Rule] := Module[
    {incDrift, chronOpt, t0Opt, chron, t0, comp, driftKind, rawDrift, drift, ctrls, hamFac, ctrlFac,
     uMat, dtVec, evs, j, sp, basis, opbasis, totalDur, tLab, tLoc, prep, dExpr, dSpace,
     dtFine, nStages, h0Fine, ctrlMats, uFine, matS},
    {incDrift, chronOpt, t0Opt} =
        {IncludeDrift, Chronology, InitialTimePoint} /. {opts} /. Options[OptimalControlEventList];
    chron = If[chronOpt === Automatic, $Chronology, chronOpt];
    t0 = If[t0Opt === Automatic, Lookup[assoc, "InitialTimePoint", 0.], N[Re[t0Opt]]];
    sp = assoc["Space"];
    basis = assoc["Basis"];
    opbasis = assoc["OperatorBasis"];
    comp = Lookup[assoc, "CompiledProblem", <||>];
    driftKind = Lookup[comp, "DriftKind",
        If[TrueQ[Lookup[assoc, "IsTimeDependentDrift", False]], "Continuous", "Static"]
    ];
    rawDrift = assoc["SymbolicDrift"];
    ctrls = toSpaceGenExpr[#, sp, basis, opbasis] & /@ assoc["SymbolicControls"];
    hamFac = assoc["DriftUnitFactor"];
    ctrlFac = assoc["ControlUnitFactor"];
    uMat = assoc["Controls"];
    dtVec = assoc["TimeSteps"];
    totalDur = assoc["TotalDuration"];

    Which[
        !TrueQ[incDrift] || driftKind === "Static",
            drift = If[TrueQ[incDrift], toSpaceGenExpr[rawDrift, sp, basis, opbasis],
                If[sp === "Liouville", NullSuperoperator[], NullOperator[]]];
            evs = Table[
                {buildSliceGeneratorExpr[drift, ctrls, hamFac, ctrlFac, uMat[[j]], incDrift, sp], dtVec[[j]]},
                {j, Length[dtVec]}
            ],
        driftKind === "SliceList",
            evs = Table[
                {buildSliceGeneratorExpr[
                    toSpaceGenExpr[rawDrift[[j]], sp, basis, opbasis],
                    ctrls, hamFac, ctrlFac, uMat[[j]], True, sp], dtVec[[j]]},
                {j, Length[dtVec]}
            ],
        True,
            tLab = Unique["ocLabT"];
            tLoc = Unique["ocLocT"];
            prep = prepareContinuousDrift[rawDrift];
            dExpr = Quiet[Check[evalContinuousExprAt[prep, tLab, t0, t0 + totalDur], $Failed]];
            If[dExpr =!= $Failed && (OperatorQ[dExpr] || SuperoperatorQ[dExpr]),
                dSpace = toSpaceGenExpr[dExpr, sp, basis, opbasis];
                evs = Table[
                    With[{tl = tLab, tloc = tLoc,
                          body = buildSliceGeneratorExpr[dSpace, ctrls, hamFac, ctrlFac, uMat[[j]], True, sp]},
                        {TimeDependentGeneratorFunction[tl, {tloc, 0}, None, body], dtVec[[j]]}
                    ],
                    {j, Length[dtVec]}
                ],
                (* Fallback for non-symbolic black-box matrix callables *)
                nStages = Lookup[comp, "NumStages", 1];
                Module[{dEval = Lookup[comp, "DriftEvaluator", None],
                        dIntSpec = Lookup[comp, "DriftIntegrationSpec", None],
                        t0Comp = Lookup[comp, "InitialTimePoint", 0.],
                        dtComp = Lookup[comp, "TimeSteps", {}], schedFallback},
                    If[dEval =!= None && AssociationQ[dIntSpec] &&
                       (Length[dtVec] =!= Length[dtComp] || Max[Abs[dtVec - dtComp]] > 10^-13 * Max[dtVec] || Abs[t0 - t0Comp] > 10^-13 * Max[1., Abs[t0]]),
                        schedFallback = buildFineDriftSchedule[dEval, dtVec, t0, dIntSpec, sp];
                        h0Fine = schedFallback["Drift"];
                        dtFine = schedFallback["FineTimeSteps"],
                        dtFine = Lookup[comp, "FineTimeSteps", dtVec];
                        h0Fine = comp["Drift"]
                    ]
                ];
                ctrlMats = comp["Controls"];
                uFine = If[nStages === 1, uMat, Developer`ToPackedArray[Flatten[ConstantArray[#, nStages] & /@ uMat, 1]]];
                evs = Table[
                    matS = h0Fine[[s]] + Sum[uFine[[s, k]] * ctrlMats[[k]], {k, Length[ctrlMats]}];
                    {toSpaceGenExpr[matS, sp, basis, opbasis], dtFine[[s]]},
                    {s, Length[dtFine]}
                ]
            ]
    ];
    If[chron === "RightToLeft", Reverse[evs], evs]
];

OptimalControlGenerator[OptimalControlResult[assoc_Association], tSym_Symbol, opts___Rule] := Module[
    {incDrift, t0Opt, t0, sp, basis, opbasis, comp, driftKind, rawDrift, drift, ctrls,
     hamFac, ctrlFac, uMat, dtVec, totalDur, tKnots, interps, expr, k, tLoc,
     prep, dExpr, dSpace, dtFine, tKnotsFine, h0Fine, pwPieces},
    {incDrift, t0Opt} = {IncludeDrift, InitialTimePoint} /. {opts} /. Options[OptimalControlGenerator];
    t0 = If[t0Opt === Automatic, Lookup[assoc, "InitialTimePoint", 0.], N[Re[t0Opt]]];
    sp = assoc["Space"];
    basis = assoc["Basis"];
    opbasis = assoc["OperatorBasis"];
    comp = Lookup[assoc, "CompiledProblem", <||>];
    driftKind = Lookup[comp, "DriftKind",
        If[TrueQ[Lookup[assoc, "IsTimeDependentDrift", False]], "Continuous", "Static"]
    ];
    rawDrift = assoc["SymbolicDrift"];
    ctrls = toSpaceGenExpr[#, sp, basis, opbasis] & /@ assoc["SymbolicControls"];
    hamFac = assoc["DriftUnitFactor"];
    ctrlFac = assoc["ControlUnitFactor"];
    uMat = assoc["Controls"];
    dtVec = assoc["TimeSteps"];
    totalDur = assoc["TotalDuration"];
    tKnots = t0 + Prepend[Accumulate[dtVec], 0.];
    interps = Table[
        With[{vals = uMat[[All, k]]},
            Quiet[Interpolation[Transpose[{tKnots, Prepend[vals, First[vals]]}], InterpolationOrder -> 0]]
        ],
        {k, Length[ctrls]}
    ];

    expr = 0;
    If[TrueQ[incDrift],
        Which[
            driftKind === "Static",
                drift = toSpaceGenExpr[rawDrift, sp, basis, opbasis];
                If[!(NumericQ[drift] && TrueQ[drift == 0]) && !MatchQ[drift, _NullOperator | _NullSuperoperator],
                    expr = hamFac * drift
                ],
            driftKind === "SliceList",
                pwPieces = Table[
                    With[{dJ = toSpaceGenExpr[rawDrift[[j]], sp, basis, opbasis]},
                        {hamFac * dJ, tKnots[[j]] <= tSym < tKnots[[j + 1]]}
                    ],
                    {j, Length[dtVec]}
                ];
                expr = Piecewise[pwPieces, hamFac * toSpaceGenExpr[Last[rawDrift], sp, basis, opbasis]],
            True,
                prep = prepareContinuousDrift[rawDrift];
                dExpr = Quiet[Check[evalContinuousExprAt[prep, tSym, t0, t0 + totalDur], $Failed]];
                If[dExpr =!= $Failed && (OperatorQ[dExpr] || SuperoperatorQ[dExpr]),
                    dSpace = toSpaceGenExpr[dExpr, sp, basis, opbasis];
                    expr = hamFac * dSpace,
                    Module[{dEval = Lookup[comp, "DriftEvaluator", None],
                            dIntSpec = Lookup[comp, "DriftIntegrationSpec", None],
                            t0Comp = Lookup[comp, "InitialTimePoint", 0.],
                            dtComp = Lookup[comp, "TimeSteps", {}], schedFallback},
                        If[dEval =!= None && AssociationQ[dIntSpec] &&
                           (Length[dtVec] =!= Length[dtComp] || Max[Abs[dtVec - dtComp]] > 10^-13 * Max[dtVec] || Abs[t0 - t0Comp] > 10^-13 * Max[1., Abs[t0]]),
                            schedFallback = buildFineDriftSchedule[dEval, dtVec, t0, dIntSpec, sp];
                            h0Fine = schedFallback["Drift"];
                            dtFine = schedFallback["FineTimeSteps"],
                            dtFine = Lookup[comp, "FineTimeSteps", dtVec];
                            h0Fine = comp["Drift"]
                        ]
                    ];
                    tKnotsFine = t0 + Prepend[Accumulate[dtFine], 0.];
                    pwPieces = Table[
                        {toSpaceGenExpr[h0Fine[[s]], sp, basis, opbasis], tKnotsFine[[s]] <= tSym < tKnotsFine[[s + 1]]},
                        {s, Length[dtFine]}
                    ];
                    expr = Piecewise[pwPieces, toSpaceGenExpr[Last[h0Fine], sp, basis, opbasis]]
                ]
        ]
    ];
    Do[
        If[!(NumericQ[ctrls[[k]]] && TrueQ[ctrls[[k]] == 0]) && !MatchQ[ctrls[[k]], _NullOperator | _NullSuperoperator],
            expr += (ctrlFac * interps[[k]][tSym]) * ctrls[[k]]
        ],
        {k, Length[ctrls]}
    ];
    If[NumericQ[expr] && TrueQ[expr == 0], expr = If[sp === "Liouville", NullSuperoperator[], NullOperator[]]];
    With[{tloc = Unique["ocGenLoc"], body = expr},
        TimeDependentGeneratorFunction[tSym, {tloc, 0}, None, body]
    ]
];

resolveChannelStructure[assoc_Association, chOpt_] := Module[
    {ctrls = assoc["SymbolicControls"], mCtrl, activeSpins, matched = {}, used = {},
     findPair, k, expr, base, s, ax, spinList, pairs, explicitList},
    mCtrl = Length[ctrls];
    If[chOpt =!= Automatic,
        explicitList = If[ListQ[chOpt] && chOpt =!= {}, chOpt, Return[$Failed]];
        pairs = Map[
            Function[item,
                Which[
                    MatchQ[item, Rule[_, _List]] && AllTrue[item[[2]], IntegerQ[#] && 1 <= # <= mCtrl &],
                        <|"Spin" -> item[[1]], "Channels" -> item[[2]], "Label" -> item[[1]]|>,
                    MatchQ[item, {_, _List}] && AllTrue[item[[2]], IntegerQ[#] && 1 <= # <= mCtrl &],
                        <|"Spin" -> item[[1]], "Channels" -> item[[2]], "Label" -> item[[1]]|>,
                    ListQ[item] && AllTrue[item, IntegerQ[#] && 1 <= # <= mCtrl &] && Length[item] >= 1,
                        <|"Spin" -> Automatic, "Channels" -> item, "Label" -> "ChGroup"|>,
                    True, $Failed
                ]
            ],
            explicitList
        ];
        If[MemberQ[pairs, $Failed] || pairs === {}, Return[$Failed]];
        Return[pairs]
    ];

    findPair = ConstantArray[None, mCtrl];
    Do[
        expr = ctrls[[k]];
        base = Replace[expr, _?NumericQ * op_ :> op];
        If[MatchQ[base, _opI] && Length[base] >= 2 && MemberQ[{"x", "y"}, base[[2]]],
            findPair[[k]] = {base[[1]], base[[2]]}
        ],
        {k, mCtrl}
    ];
    If[AllTrue[findPair, # =!= None &],
        spinList = DeleteDuplicates[findPair[[All, 1]]];
        matched = Table[
            With[{sp = spinList[[idx]],
                  xs = Flatten[Position[findPair, {spinList[[idx]], "x"}, {1}]],
                  ys = Flatten[Position[findPair, {spinList[[idx]], "y"}, {1}]]},
                Which[
                    Length[xs] >= 1 && Length[ys] >= 1,
                        <|"Spin" -> sp, "Channels" -> {First[xs], First[ys]}, "Label" -> sp|>,
                    Length[xs] >= 1,
                        <|"Spin" -> sp, "Channels" -> {First[xs]}, "Label" -> sp|>,
                    Length[ys] >= 1,
                        <|"Spin" -> sp, "Channels" -> {First[ys]}, "Label" -> sp|>,
                    True, $Failed
                ]
            ],
            {idx, Length[spinList]}
        ];
        If[!MemberQ[matched, $Failed] && Flatten[#["Channels"] & /@ matched] === Range[mCtrl],
            Return[matched]
        ]
    ];

    (* Fallback pairing *)
    activeSpins = Quiet[Check[First /@ SpinSystem[assoc["Basis"]], First /@ SpinSystem[]]];
    If[!ListQ[activeSpins] || activeSpins === {}, activeSpins = {1}];
    Which[
        mCtrl <= 2,
            {<|"Spin" -> If[Length[activeSpins] === 1, First[activeSpins], activeSpins], "Channels" -> Range[mCtrl], "Label" -> 1|>},
        EvenQ[mCtrl],
            Table[<|"Spin" -> If[idx <= Length[activeSpins], activeSpins[[idx]], idx],
                    "Channels" -> {2 idx - 1, 2 idx}, "Label" -> idx|>, {idx, mCtrl / 2}],
        True,
            Table[<|"Spin" -> If[k <= Length[activeSpins], activeSpins[[k]], k],
                    "Channels" -> {k}, "Label" -> k|>, {k, mCtrl}]
    ]
];

OptimalControlPulseSequence[OptimalControlResult[assoc_Association], spins_List, opts___Rule] := Module[
    {incDriftOpt, chOpt, chList, uMat, dtVec, tKnots, totalDur, ctrlFac, hamFac, tSym, tLoc,
     ux, uy, ampVals, phVals, ampInterp, phInterp, ctrls, interps, expr, sp, basis, opbasis,
     comp, driftKind, rawDrift, drift, prep, dExpr},
    {incDriftOpt, chOpt} = {IncludeDrift, ChannelStructure} /. {opts} /. Options[OptimalControlPulseSequence];
    chList = resolveChannelStructure[assoc, chOpt];
    If[chList === $Failed, Message[ExportOptimalControlPulse::badchan, chOpt]; Return[$Failed]];
    uMat = assoc["Controls"];
    dtVec = assoc["TimeSteps"];
    tKnots = Prepend[Accumulate[dtVec], 0.];
    totalDur = assoc["TotalDuration"];
    ctrlFac = assoc["ControlUnitFactor"];
    hamFac = assoc["DriftUnitFactor"];
    sp = assoc["Space"];
    basis = assoc["Basis"];
    opbasis = assoc["OperatorBasis"];

    If[Length[chList] === 1 && MemberQ[{1, 2}, Dimensions[uMat][[2]]] && !TrueQ[incDriftOpt],
        tSym = Unique["ocPulseTime"];
        With[{cIdx = chList[[1, "Channels"]]},
            ux = ctrlFac * uMat[[All, cIdx[[1]]]];
            uy = If[Length[cIdx] === 2, ctrlFac * uMat[[All, cIdx[[2]]]], ConstantArray[0., Length[dtVec]]];
        ];
        ampVals = Sqrt[ux^2 + uy^2];
        phVals = MapThread[If[#1 == 0. && #2 == 0., 0., ArcTan[#1, #2]] &, {ux, uy}];
        ampInterp = Quiet[Interpolation[Transpose[{tKnots, Prepend[ampVals, First[ampVals]]}], InterpolationOrder -> 0]];
        phInterp = Quiet[Interpolation[Transpose[{tKnots, Prepend[phVals, First[phVals]]}], InterpolationOrder -> 0]];
        Return[PulseSequence[spins, {
            ShapedPulse[spins, totalDur, {{tSym, 0}, {ampInterp[tSym], phInterp[tSym], 0}}]
        }]]
    ];

    (* Simultaneous multi-channel (or IncludeDrift -> True) pulse sequence over local time tLoc in [0, totalDur] *)
    tSym = Unique["ocSeqT"];
    tLoc = Unique["ocSeqLoc"];
    ctrls = toSpaceGenExpr[#, sp, basis, opbasis] & /@ assoc["SymbolicControls"];
    interps = Table[
        With[{vals = uMat[[All, k]]},
            Quiet[Interpolation[Transpose[{tKnots, Prepend[vals, First[vals]]}], InterpolationOrder -> 0]]
        ],
        {k, Length[ctrls]}
    ];
    expr = 0;
    If[TrueQ[incDriftOpt],
        comp = Lookup[assoc, "CompiledProblem", <||>];
        driftKind = Lookup[comp, "DriftKind",
            If[TrueQ[Lookup[assoc, "IsTimeDependentDrift", False]], "Continuous", "Static"]
        ];
        rawDrift = assoc["SymbolicDrift"];
        Which[
            driftKind === "Static",
                drift = toSpaceGenExpr[rawDrift, sp, basis, opbasis];
                If[!(NumericQ[drift] && TrueQ[drift == 0]) && !MatchQ[drift, _NullOperator | _NullSuperoperator],
                    expr = hamFac * drift
                ],
            True,
                prep = prepareContinuousDrift[rawDrift];
                dExpr = Quiet[Check[evalContinuousExprAt[prep, tLoc, 0., totalDur], $Failed]];
                If[dExpr =!= $Failed && (OperatorQ[dExpr] || SuperoperatorQ[dExpr]),
                    expr = hamFac * toSpaceGenExpr[dExpr, sp, basis, opbasis]
                ]
        ]
    ];
    Do[
        If[!(NumericQ[ctrls[[k]]] && TrueQ[ctrls[[k]] == 0]) && !MatchQ[ctrls[[k]], _NullOperator | _NullSuperoperator],
            expr += (ctrlFac * interps[[k]][tLoc]) * ctrls[[k]]
        ],
        {k, Length[ctrls]}
    ];
    If[NumericQ[expr] && TrueQ[expr == 0], expr = If[sp === "Liouville", NullSuperoperator[], NullOperator[]]];
    With[{tl = tSym, tloc = tLoc, body = expr},
        PulseSequence[spins, {
            Delay[totalDur, TimeDependentGeneratorFunction[tl, {tloc, 0}, None, body]]
        }]
    ]
];

OptimalControlPulseSequence[res:OptimalControlResult[assoc_Association], opts___Rule] := Module[
    {chOpt = ChannelStructure /. {opts} /. Options[OptimalControlPulseSequence], chList, spins},
    chList = resolveChannelStructure[assoc, chOpt];
    If[chList === $Failed, Message[ExportOptimalControlPulse::badchan, chOpt]; Return[$Failed]];
    spins = DeleteCases[DeleteDuplicates[Flatten[#["Spin"] & /@ chList]], Automatic];
    If[spins === {}, spins = Quiet[Check[First /@ SpinSystem[assoc["Basis"]], First /@ SpinSystem[]]]];
    If[!ListQ[spins] || spins === {}, spins = {1}];
    OptimalControlPulseSequence[res, spins, opts]
];

OptimalControlPulseSequence[res_OptimalControlResult, spin:Except[_List | _Rule], opts___Rule] :=
    OptimalControlPulseSequence[res, {spin}, opts];

(* Structured pulse waveform table & Bruker / SIMPSON / CSV hardware exporters *)
buildPolarChannelEntry[uMat_, dtVec_List, ctrlFac_?NumericQ, chEntry_Association, phUnits_String] := Module[
    {cIdx = chEntry["Channels"], nt = Length[dtVec], ux, uy, ampRad, ampHz, maxRad, maxHz,
     pct, phRad, phDeg, phOut},
    ux = (ctrlFac / (2. * Pi)) * uMat[[All, cIdx[[1]]]];
    uy = If[Length[cIdx] >= 2, (ctrlFac / (2. * Pi)) * uMat[[All, cIdx[[2]]]], ConstantArray[0., nt]];
    ampHz = Developer`ToPackedArray[Sqrt[ux^2 + uy^2]];
    ampRad = Developer`ToPackedArray[(2. * Pi) * ampHz];
    maxHz = Max[ampHz];
    maxRad = Max[ampRad];
    pct = Developer`ToPackedArray[If[maxHz > 0., 100. * ampHz / maxHz, ConstantArray[0., nt]]];
    phRad = Developer`ToPackedArray[MapThread[If[#1 == 0. && #2 == 0., 0., ArcTan[#1, #2]] &, {ux, uy}]];
    phDeg = Developer`ToPackedArray[Mod[phRad * 180. / Pi, 360.]];
    phOut = If[phUnits === "Radians", phRad, phDeg];
    <|
        "Spin" -> chEntry["Spin"],
        "Label" -> chEntry["Label"],
        "ChannelIndices" -> cIdx,
        "UxHz" -> Developer`ToPackedArray[ux],
        "UyHz" -> Developer`ToPackedArray[uy],
        "AmplitudeHz" -> ampHz,
        "AmplitudeRadPerSec" -> ampRad,
        "RelativeAmplitudePercent" -> pct,
        "MaxAmplitudeHz" -> maxHz,
        "MaxAmplitudeRadPerSec" -> maxRad,
        "PhaseRadians" -> phRad,
        "PhaseDegrees" -> phDeg,
        "Phase" -> phOut
    |>
];

OptimalControlPulseTable[OptimalControlResult[assoc_Association], opts___Rule] := Module[
    {chOpt, phUnits, ampScale, chList, uMat, dtVec, ctrlFac, tEdges, tMid, uniformQ,
     chEntries, tableRows, j, c, ampVals},
    {chOpt, phUnits, ampScale} = {ChannelStructure, PhaseUnits, AmplitudeScaling} /. {opts} /. Options[OptimalControlPulseTable];
    chList = resolveChannelStructure[assoc, chOpt];
    If[chList === $Failed, Message[ExportOptimalControlPulse::badchan, chOpt]; Return[$Failed]];
    uMat = assoc["Controls"];
    dtVec = Developer`ToPackedArray[N[assoc["TimeSteps"]]];
    ctrlFac = N[assoc["ControlUnitFactor"]];
    tEdges = Developer`ToPackedArray[Prepend[Accumulate[dtVec], 0.]];
    tMid = Developer`ToPackedArray[0.5 * (Most[tEdges] + Rest[tEdges])];
    uniformQ = TrueQ[Max[dtVec] - Min[dtVec] <= 10^-9 * Mean[dtVec]];
    chEntries = Map[buildPolarChannelEntry[uMat, dtVec, ctrlFac, #, phUnits] &, chList];
    tableRows = Table[
        Flatten[{
            tMid[[j]], dtVec[[j]],
            Table[
                With[{ch = chEntries[[c]]},
                    {If[ampScale === "Hz", (ch["AmplitudeHz"])[[j]], (ch["RelativeAmplitudePercent"])[[j]]],
                     (ch["Phase"])[[j]]}
                ],
                {c, Length[chEntries]}
            ]
        }],
        {j, Length[dtVec]}
    ];
    <|
        "NumSlices" -> Length[dtVec],
        "TotalDuration" -> Total[dtVec],
        "TimeSteps" -> dtVec,
        "TimeEdges" -> tEdges,
        "TimeMidpoints" -> tMid,
        "UniformTimeStepQ" -> uniformQ,
        "PhaseUnits" -> phUnits,
        "AmplitudeScaling" -> ampScale,
        "Channels" -> chEntries,
        "Table" -> Developer`ToPackedArray[tableRows]
    |>
];

fmtCleanReal[x_?NumericQ] := With[{v = Chop[N[Re[x]], 10^-16]},
    If[v == 0., "0.0", ToString[CForm[v]]]
];

fmtFixedDec[x_?NumericQ, digits_Integer: 6] := StringTrim[ToString[
    NumberForm[Chop[N[Re[x]], 10^-14], {14, digits}, ExponentFunction -> (Null &), NumberPadding -> {"", "0"}]
]];

resampleWaveformGrid[uMat_, dtVec_List, nPts_Integer] := Module[
    {totalDur = Total[dtVec], tEdges, tMidsOut, mCtrl = Dimensions[uMat][[2]], dtOut, uOut, idx, s},
    tEdges = Prepend[Accumulate[dtVec], 0.];
    dtOut = Developer`ToPackedArray[ConstantArray[totalDur / nPts, nPts]];
    tMidsOut = (Range[nPts] - 0.5) * (totalDur / nPts);
    uOut = Table[
        s = FirstPosition[Rest[tEdges], e_?NumericQ /; e >= tMidsOut[[k]] - 10^-15 * totalDur, {Length[dtVec]}, {1}, Heads -> False][[1]];
        uMat[[s]],
        {k, nPts}
    ];
    {Developer`ToPackedArray[uOut], dtOut}
];

formatSingleChannelBruker[chEntry_Association, dtVec_List, title_String, chIdx_Integer] := Module[
    {pct = chEntry["RelativeAmplitudePercent"], deg = chEntry["PhaseDegrees"],
     nPts = Length[dtVec], totalDur = Total[dtVec], integFac, lines, header},
    integFac = Mean[pct / 100.];
    header = {
        "##TITLE= " <> title <> " (Ch " <> ToString[chIdx] <> ")",
        "##JCAMP-DX= 5.00 Bruker JCAMP library",
        "##DATA TYPE= Shape Data",
        "##ORIGIN= SpinDynamica OptimalControl",
        "##OWNER= SpinDynamica",
        "##MINX= 0.000000",
        "##MAXX= 100.000000",
        "##MINY= 0.000000",
        "##MAXY= 360.000000",
        "##$SHAPE_EXMODE= None",
        "##$SHAPE_TOTROT= 0.000000",
        "##$SHAPE_BWFAC= 0.000000",
        "##$SHAPE_INTEGFAC= " <> fmtCleanReal[integFac],
        "##$SHAPE_MODE= 1",
        "##$SD_TOTAL_DURATION_S= " <> fmtCleanReal[totalDur],
        "##$SD_MAX_AMPLITUDE_HZ= " <> fmtCleanReal[chEntry["MaxAmplitudeHz"]],
        "##$SD_MAX_AMPLITUDE_RADS= " <> fmtCleanReal[chEntry["MaxAmplitudeRadPerSec"]],
        "##NPOINTS= " <> ToString[nPts],
        "##XYPOINTS= (XY..XY)"
    };
    lines = MapThread[fmtFixedDec[#1, 6] <> ", " <> fmtFixedDec[#2, 6] &, {pct, deg}];
    StringRiffle[Join[header, lines, {"##END="}], "\n"] <> "\n"
];

formatSingleChannelSIMPSON[chEntry_Association] := Module[
    {amp = chEntry["AmplitudeHz"], deg = chEntry["PhaseDegrees"]},
    StringRiffle[MapThread[fmtCleanReal[#1] <> " " <> fmtCleanReal[#2] &, {amp, deg}], "\n"] <> "\n"
];

formatMultiChannelCSV[chEntries_List, dtVec_List, selIndices_List] := Module[
    {tEdges = Prepend[Accumulate[dtVec], 0.], tMid, headerCols, rows},
    tMid = 0.5 * (Most[tEdges] + Rest[tEdges]);
    headerCols = Join[
        {"Slice", "TimeStart_s", "TimeEnd_s", "TimeMid_s", "Duration_s"},
        Flatten[Table[
            With[{p = "Ch" <> ToString[c]},
                {p <> "_Ux_Hz", p <> "_Uy_Hz", p <> "_Amp_Hz", p <> "_Amp_Pct", p <> "_Phase_Deg"}
            ],
            {c, selIndices}
        ]]
    ];
    rows = Table[
        StringRiffle[Flatten[{
            ToString[j], fmtCleanReal[tEdges[[j]]], fmtCleanReal[tEdges[[j + 1]]],
            fmtCleanReal[tMid[[j]]], fmtCleanReal[dtVec[[j]]],
            Table[
                With[{ch = chEntries[[c]]},
                    {fmtCleanReal[(ch["UxHz"])[[j]]], fmtCleanReal[(ch["UyHz"])[[j]]],
                     fmtCleanReal[(ch["AmplitudeHz"])[[j]]], fmtCleanReal[(ch["RelativeAmplitudePercent"])[[j]]],
                     fmtCleanReal[(ch["PhaseDegrees"])[[j]]]}
                ],
                {c, selIndices}
            ]
        }], ","],
        {j, Length[dtVec]}
    ];
    StringRiffle[Prepend[rows, StringRiffle[headerCols, ","]], "\n"] <> "\n"
];

ExportOptimalControlPulse[res:OptimalControlResult[assoc_Association], format_String, opts___Rule] := Module[
    {chOpt, chStructOpt, resampOpt, titleOpt, fmtCanon, chList, nCh, uMat, dtVec,
     uniformQ, resampled, ctrlFac, chEntries, selIdx, outMap},
    If[!MemberQ[{"Bruker", "SIMPSON", "CSV"}, format],
        Message[ExportOptimalControlPulse::badfmt, format]; Return[$Failed]
    ];
    fmtCanon = format;
    {chOpt, chStructOpt, resampOpt, titleOpt} =
        {Channel, ChannelStructure, Resample, Title} /. {opts} /. Options[ExportOptimalControlPulse];
    chList = resolveChannelStructure[assoc, chStructOpt];
    If[chList === $Failed, Message[ExportOptimalControlPulse::badchan, chStructOpt]; Return[$Failed]];
    nCh = Length[chList];
    uMat = assoc["Controls"];
    dtVec = Developer`ToPackedArray[N[assoc["TimeSteps"]]];
    ctrlFac = N[assoc["ControlUnitFactor"]];
    uniformQ = TrueQ[Max[dtVec] - Min[dtVec] <= 10^-9 * Mean[dtVec]];

    Which[
        IntegerQ[resampOpt] && resampOpt >= 1,
            {uMat, dtVec} = resampleWaveformGrid[uMat, dtVec, resampOpt],
        (fmtCanon === "Bruker" || fmtCanon === "SIMPSON") && !uniformQ,
            Message[ExportOptimalControlPulse::nonuniform, fmtCanon, Length[dtVec]]
    ];

    chEntries = Map[buildPolarChannelEntry[uMat, dtVec, ctrlFac, #, "Degrees"] &, chList];

    selIdx = Which[
        chOpt === Automatic, If[fmtCanon === "CSV", Range[nCh], 1],
        chOpt === All, Range[nCh],
        IntegerQ[chOpt] && 1 <= chOpt <= nCh, chOpt,
        ListQ[chOpt] && chOpt =!= {} && AllTrue[chOpt, IntegerQ[#] && 1 <= # <= nCh &], chOpt,
        True, $Failed
    ];
    If[selIdx === $Failed, Message[ExportOptimalControlPulse::badchan, chOpt]; Return[$Failed]];

    If[fmtCanon === "CSV",
        Return[formatMultiChannelCSV[chEntries, dtVec, If[IntegerQ[selIdx], {selIdx}, selIdx]]]
    ];
    If[IntegerQ[selIdx],
        If[fmtCanon === "Bruker",
            formatSingleChannelBruker[chEntries[[selIdx]], dtVec, titleOpt, selIdx],
            formatSingleChannelSIMPSON[chEntries[[selIdx]]]
        ],
        outMap = Association[Map[
            Function[c,
                c -> If[fmtCanon === "Bruker",
                    formatSingleChannelBruker[chEntries[[c]], dtVec, titleOpt, c],
                    formatSingleChannelSIMPSON[chEntries[[c]]]
                ]
            ],
            selIdx
        ]];
        outMap
    ]
];

ExportOptimalControlPulse[path_String, res_OptimalControlResult, format_: Automatic, opts___Rule] := Module[
    {ext = ToLowerCase[FileExtension[path]], effFmt, payload, dir, base, fExt, paths},
    effFmt = If[format =!= Automatic,
        format,
        Switch[ext,
            "csv", "CSV",
            "shp" | "bruker" | "jcamp", "Bruker",
            "rf" | "simpson" | "txt", "SIMPSON",
            _, "CSV"
        ]
    ];
    payload = ExportOptimalControlPulse[res, effFmt, opts];
    If[payload === $Failed, Return[$Failed]];
    dir = DirectoryName[path];
    If[dir =!= "" && !DirectoryQ[dir], CreateDirectory[dir]];
    If[StringQ[payload],
        Export[path, payload, "String"];
        Return[path]
    ];
    (* Association of multiple channel payloads *)
    base = FileNameJoin[{dir, FileBaseName[path]}];
    fExt = If[ext === "", Switch[effFmt, "Bruker", "shp", "SIMPSON", "RF", _, "txt"], ext];
    paths = KeyValueMap[
        With[{p = base <> "_ch" <> ToString[#1] <> "." <> fExt},
            Export[p, #2, "String"];
            p
        ] &,
        payload
    ];
    paths
];

(* Independent verification against SpinDynamica native propagators *)
VerifyOptimalControlSimulation[res:OptimalControlResult[assoc_Association], opts___Rule] := Module[
    {tol, chronOpt, t0Opt, chron, t0, evList, basis, opbasis, space, uSD, rhoFinalMat, sdRes, ocFid, sdFid,
     propErr = 0., fidErr, rawTarget, evalLeaf, optPhases, sdProps, targetNormObj, sdTA},
    If[KeyExistsQ[$OCSolverRegistry, "VerifyEnsemble"] && assoc["Kind"] === "Ensemble",
        Return[($OCSolverRegistry["VerifyEnsemble"])[res, {opts}]]
    ];
    {tol, chronOpt, t0Opt} =
        {Tolerance, Chronology, InitialTimePoint} /. {opts} /. Options[VerifyOptimalControlSimulation];
    chron = If[chronOpt === Automatic, $Chronology, chronOpt];
    t0 = If[t0Opt === Automatic, Lookup[assoc, "InitialTimePoint", 0.], N[Re[t0Opt]]];
    evList = OptimalControlEventList[res, Chronology -> chron, InitialTimePoint -> t0];
    basis = assoc["Basis"];
    opbasis = assoc["OperatorBasis"];
    space = assoc["Space"];
    ocFid = assoc["FinalFidelity"];

    If[space === "Hilbert",
        uSD = ToPackedComplexMatrix[MatrixRepresentation[
            NPropagationOperator[evList, Basis -> basis, Chronology -> chron, InitialTimePoint -> t0], basis]];
        If[assoc["FinalPropagator"] =!= None,
            propErr = Max[Abs[uSD - assoc["FinalPropagator"]]],
            propErr = 0.
        ];
        evalLeaf = CompileObjective[assoc["CompiledProblem", "RawTarget"], "Hilbert", basis, opbasis];
        optPhases = Lookup[assoc, "OptimalPhases", None];
        sdFid = If[AssociationQ[evalLeaf["PhaseFreedomSpec"]] && AssociationQ[optPhases],
            First[(evalLeaf["PhaseFreedomSpec", "EvaluateWithPhases"])[uSD, optPhases["Vector"]]],
            First[(evalLeaf["EvaluateCotangent"])[uSD]]
        ];
        fidErr = Abs[sdFid - ocFid];
        <|
            "FidelityMatchQ" -> TrueQ[fidErr <= tol && propErr <= tol],
            "OptimalControlFidelity" -> ocFid,
            "SpinDynamicaFidelity" -> sdFid,
            "FidelityError" -> fidErr,
            "PropagatorError" -> propErr,
            "SpinDynamicaResult" -> uSD
        |>,
        (* Liouville verification: propagate initial operators and compare with FinalState *)
        Module[{initList, tarList, weights, modes, scales, initMat, tarMat, ocStates, sdStates,
                initOp, finalOp, finalVec, zVec, subFids},
            rawTarget = assoc["CompiledProblem", "RawTarget"];
            evalLeaf = CompileObjective[rawTarget, "Liouville", basis, opbasis];
            initMat = evalLeaf["InitialMatrix"];
            tarMat = evalLeaf["TargetMatrix"];
            weights = evalLeaf["Weights"];
            modes = evalLeaf["Modes"];
            scales = evalLeaf["Scales"];
            sdStates = Transpose[Table[
                initOp = Operator[initMat[[All, r]], opbasis, basis, CheckOperatorType -> False];
                finalOp = NPropagate[evList, OperatorBasis -> opbasis, Basis -> basis,
                    Chronology -> chron, InitialTimePoint -> t0, UseHilbertSpaceWhenPossible -> False][initOp];
                ToPackedComplexVector[Flatten[Normal[OperatorVectorRepresentation[finalOp, opbasis, basis]]]],
                {r, Length[weights]}
            ]];
            If[assoc["FinalState"] =!= None,
                propErr = Max[Abs[sdStates - assoc["FinalState"]]]
            ];
            zVec = Total[Conjugate[tarMat] * sdStates, {1}];
            subFids = Table[
                Switch[modes[[r]],
                    "Magnitude", Abs[zVec[[r]] / scales[[r]]]^2,
                    "Distance", 1. - Re[Total[Abs[sdStates[[All, r]] - tarMat[[All, r]]]^2]] / scales[[r]],
                    _, Re[zVec[[r]]] / scales[[r]]
                ],
                {r, Length[weights]}
            ];
            sdFid = weights . subFids;
            fidErr = Abs[sdFid - ocFid];
            <|
                "FidelityMatchQ" -> TrueQ[fidErr <= tol && propErr <= tol],
                "OptimalControlFidelity" -> ocFid,
                "SpinDynamicaFidelity" -> sdFid,
                "FidelityError" -> fidErr,
                "PropagatorError" -> propErr,
                "SpinDynamicaResult" -> sdStates
            |>
        ]
    ]
];

$OCSolverRegistry["RunScaledOptimizer"] = runScaledOptimizer;
$OCSolverRegistry["RunFISTAOptimizer"] = runFISTAOptimizer;
$OCSolverRegistry["DispatchOptimizer"] = dispatchOptimizationEngine;
$OCSolverRegistry["BuildInitialVector"] = buildInitialVector;
$OCSolverRegistry["ResolveTimingVector"] = resolveTimingVector;
$OCSolverRegistry["ResolveUnitFactor"] = resolveUnitFactor;
$OCSolverRegistry["ResolveBoundMatrices"] = resolveBoundMatrices;
$OCSolverRegistry["ResolveConstraintSpec"] = resolveConstraintSpec;
$OCSolverRegistry["CompileProxContext"] = compileProxContext;
$OCSolverRegistry["ApplyProximal"] = applyProximal;
$OCSolverRegistry["EvalNonSmoothPenalty"] = evalNonSmoothPenalty;
$OCSolverRegistry["ComputeConstraintViolation"] = computeConstraintViolation;
$OCSolverRegistry["ToGenMatrix"] = toGenMatrix;
$OCSolverRegistry["ParseMonitorConfig"] = parseMonitorConfig;
$OCSolverRegistry["ParseTimeStepSchedule"] = parseTimeStepSchedule;
$OCSolverRegistry["ParseMultiStartSpec"] = parseMultiStartSpec;
$OCSolverRegistry["RunMultiStart"] = runMultiStartOrchestrator;
$OCSolverRegistry["MonitorStageKey"] = $MonitorStage;
$OCSolverRegistry["CompileDrift"] = compileDriftRepresentation;
$OCSolverRegistry["DriftContainsSuperoperatorQ"] = driftContainsSuperoperatorQ;
$OCSolverRegistry["BuildFineDriftSchedule"] = buildFineDriftSchedule;
$OCSolverRegistry["ParseDriftIntegrationSpec"] = parseDriftIntegrationSpec;
$OCSolverRegistry["ProjectBoundedSimplex"] = projectBoundedSimplex;
$OCSolverRegistry["ResolveDurationSpec"] = resolveDurationSpec;
$OCSolverRegistry["ParseDurationPenalty"] = parseDurationPenaltySpec;
$OCSolverRegistry["BuildTimingDerivTensors"] = buildTimingDerivTensors;
$OCSolverRegistry["ExtractOptimalPhases"] = extractOptimalPhasesAssoc;
$OCSolverRegistry["ResampleControls"] = resampleControls;
$OCSolverRegistry["AssembleSingleResult"] = assembleSingleResult;
$OCSolverRegistry["FidelityCriterionKey"] = $FidelityCriterionKey;
$OCSolverRegistry["CompileSingle"] = CompileOptimalControlProblem;
$OCSolverRegistry["EvaluateSingleObjective"] = evaluateOptimalControlObjective;
$OCSolverRegistry["ComputeSmoothPenaltyValue"] = computeSmoothPenaltyValue;
$OCSolverRegistry["UnpackControlMatrix"] = unpackControlMatrix;
$OCSolverRegistry["ResolveLineSearchConfig"] = resolveLineSearchConfig;
$OCSolverRegistry["EvalCriterionFidelity"] = evalCriterionFidelity;
$OCSolverRegistry["EmitMonitorIfDue"] = emitMonitorIfDue;
$OCSolverRegistry["ProjectedArmijoLineSearch"] = projectedArmijoLineSearch;
$OCSolverRegistry["StrongWolfeLineSearch"] = strongWolfeLineSearch;
$OCSolverRegistry["SecondOrderCapabilities"] = checkSecondOrderCapabilities;
$OCSolverRegistry["EvaluateWithCache"] = evaluateOptimalControlWithCache;
$OCSolverRegistry["EvaluateHVP"] = evaluateOptimalControlHVP;
$OCSolverRegistry["EvaluateHessian"] = evaluateOptimalControlHessian;
$OCSolverRegistry["ApplySmoothPenaltyHVP"] = applySmoothPenaltyHVP;
$OCSolverRegistry["BuildSmoothPenaltyHessian"] = buildSmoothPenaltyHessian;
$OCSolverRegistry["PullbackSliceHessianToBasis"] = pullbackSliceHessianToBasis;
$OCSolverRegistry["BuildG2Tensor"] = OptimalControl`Kernels`Private`buildG2Tensor;
If[!KeyExistsQ[$OCSolverRegistry, "LineSearches"], $OCSolverRegistry["LineSearches"] = <||>];
If[!KeyExistsQ[$OCSolverRegistry, "Optimizers"], $OCSolverRegistry["Optimizers"] = <||>];

End[];
EndPackage[];
