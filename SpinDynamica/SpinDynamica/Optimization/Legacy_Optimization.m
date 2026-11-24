(* ::Package:: *)

(* :Title: Simulations`Optimization` *)
(* :Context: Simulations`Optimization` *)
(* :Summary: Authoritative .m package providing MaximizeTransformationAmplitude with unified Engine dispatch (Automatic | "Compiled" | "Legacy") and backwards compatibility. *)

BeginPackage["Simulations`Optimization`",
{"Global`",
"SpinDynamica`",
"General`UserLevel`",
"QM`Hilbert`",
"QM`SpinOperators`",
"QM`Liouville`",
"SD`Generators`",
"SD`Approximations`",
"SD`Propagation`",
"Events`Events`",
"Events`PulseSequence`",
"Interactions`Interactions`",
"Interactions`Hamiltonians`",
"Interactions`Relaxation`",
"Simulations`Polarization`",
"Signals`Signal`",
"Ensemble`EnsembleAverage`",
"General`Matrices`",
"General`Miscellaneous`",
"Rotations`Wigner`",
"Units`",
"General`PhysicalConstants`",
"Simulations`Simulations`"
}];

MaximizeTransformationAmplitude::usage="MaximizeTransformationAmplitude[\[Rho]ini\[Rule]obs,events,{{p1,p1ini},{p2,p2ini}..}, opts] numerically optimizes a TransformationAmplitude quantity with respect to parameters {p1,p2..}, starting from initial parameter values {p1ini,p2ini..}. 
The option rules opts may include BackgroundGenerator and any of the other usual option specifications for TransformationAmplitude. The values specified in opts may depend on the parameters {p1,p2..}.\nThe option rules opts may include the directive ConstrainPositive\[Rule]True which constrains all parameters {p1,p2..} to positive values. The syntax ConstrainPositive\[Rule]{p1,p2..} may be used to constrain only a subset of parameters to positive values. This is frequently required for pulse sequence delays, etc., which cannot take negative values. 
More general constraints may be implemented by using MaximizeTransformationAmplitude[\[Rho]ini\[Rule]obs,events,{{p1,p1ini,{p1min,p1max}}..}, opts] subjects the parameter p1 to the constraint p1min\[NonBreakingSpace]\[LessEqual]\[NonBreakingSpace]p1\[NonBreakingSpace]\[LessEqual]\[NonBreakingSpace]p1max. If the maximum value p1max is omitted, there is no upper constraint. Any of the parameters {p1,p2..} may be constrained in this way.  
Constrained optimization is often computationally more intensive than unconstrained optimization.
MaximizeTransformationAmplitude supports Engine -> Automatic (default, uses compiled SpinDynamica optimization kernels with analytical/adjoint derivatives and legacy fallback), Engine -> \"Compiled\", or Engine -> \"Legacy\" (FindMinimum over TransformationAmplitude). By default it returns {maxval, rules}; pass ReturnResultObject -> True for a full TransformationOptimizationResult.";

ConstrainPositive::usage="ConstrainPositive is an option for MaximizeTransformationAmplitude. If ConstrainPositive\[Rule]True or ConstrainPositive\[Rule]All, all variable parameters of  MaximizeTransformationAmplitude are constrained to positive values. This is frequently required for pulse sequence delays, etc., which cannot take negative values. The syntax ConstrainPositive\[Rule]{p1,p2..} may be used to constrain only a subset of parameters to positive values.";

Engine::usage="Engine is an option for MaximizeTransformationAmplitude specifying the optimization backend: Automatic (default, compiled engine with legacy fallback), \"Compiled\" (high-performance compiled kernels and analytical/adjoint derivatives), or \"Legacy\" (FindMinimum over TransformationAmplitude).";

Begin["`Private`"];

DownValues[MaximizeTransformationAmplitude] = {};
SetAttributes[MaximizeTransformationAmplitude, HoldAll];

Options[MaximizeTransformationAmplitude] = {
    Method -> Automatic,
    ConstrainPositive -> Automatic,
    Engine -> Automatic
};

(* ------------------------------------------------------------------ *)
(* Held argument & variable specification helpers                     *)
(* ------------------------------------------------------------------ *)
normalizeBoundSpec[HoldComplete[{v_Symbol, v0_}]] :=
    {Hold[v], N[Re[v0]], {-Infinity, Infinity}, False};
normalizeBoundSpec[HoldComplete[{v_Symbol, v0_, c_List}]] /; Length[c] === 1 :=
    {Hold[v], N[Re[v0]], {N[Re[First[c]]], Infinity}, True};
normalizeBoundSpec[HoldComplete[{v_Symbol, v0_, c_List}]] /; Length[c] === 2 :=
    {
        Hold[v],
        N[Re[v0]],
        {
            If[c[[1]] === -Infinity || c[[1]] === None, -Infinity, N[Re[c[[1]]]]],
            If[c[[2]] === Infinity || c[[2]] === None, Infinity, N[Re[c[[2]]]]]
        },
        True
    };
normalizeBoundSpec[HoldComplete[{v_Symbol, v0_, c:Except[_List]}]] :=
    {Hold[v], N[Re[v0]], {N[Re[c]], Infinity}, True};
normalizeBoundSpec[_] := $Failed;

SetAttributes[extractHeldVarEntries, HoldAllComplete];
extractHeldVarEntries[sym_Symbol] /; OwnValues[sym] =!= {} :=
    Replace[First[OwnValues[sym]], (_ :> rhs_) :> extractHeldVarEntries[rhs]];
extractHeldVarEntries[{v:Except[_List], v0:Except[_List]}] :=
    extractHeldVarEntries[{{v, v0}}];
extractHeldVarEntries[{v:Except[_List], v0:Except[_List], c_}] :=
    extractHeldVarEntries[{{v, v0, c}}];
extractHeldVarEntries[specList_List] := Module[{heldEntries},
    heldEntries = Replace[HoldComplete[specList], HoldComplete[{items___}] :> (HoldComplete /@ HoldComplete[items])];
    List @@ (normalizeBoundSpec /@ heldEntries)
];
extractHeldVarEntries[_] := $Failed;

(* Evaluate held expression inside Block binding all optimization symbols *)
evalHeldAtNumeric[Hold[vs___], bodyHold_Hold, vals_List] := Module[{rules},
    rules = Table[
        With[{pat = Extract[HoldComplete[vs], i, HoldPattern], val = vals[[i]]}, pat -> val],
        {i, Length[vals]}
    ];
    Block[{vs},
        With[{left = HoldComplete[vs]},
            Do[
                With[{s = Extract[left, i, HoldComplete], v = vals[[i]]},
                    Replace[s, HoldComplete[sym_] :> (sym = v)]
                ],
                {i, Length[vals]}
            ]
        ];
        ReleaseHold[bodyHold /. rules]
    ]
];

(* Exact bound reflection preventing out-of-domain probes in FindMinimum *)
reflectIntoBounds[val_, lo_, hi_] := With[{v = N[Re[val]]},
    Which[
        lo =!= -Infinity && hi =!= Infinity,
            With[{a = N[lo], b = N[hi], span = N[hi] - N[lo]},
                If[span <= 0.,
                    a,
                    With[{u = Mod[v - a, 2. * span]},
                        If[u <= span, a + u, b - (u - span)]
                    ]
                ]
            ],
        lo =!= -Infinity && v < lo,
            N[lo] + Abs[v - N[lo]],
        hi =!= Infinity && v > hi,
            N[hi] - Abs[v - N[hi]],
        True,
            v
    ]
];

(* Extract a control option from HoldComplete[opts] without evaluating unrelated option values *)
extractHeldOptionValue[heldOpts_HoldComplete, targetPatt_, default_] := Module[{hits},
    hits = Cases[heldOpts, (Rule|RuleDelayed)[k_ /; MatchQ[Unevaluated[k], targetPatt], v_] :> HoldComplete[v], {1}];
    If[hits === {}, default, ReleaseHold[First[hits]]]
];

hasOptionKeyQ[heldOpts_HoldComplete, keyPatt_] :=
    MemberQ[heldOpts, (Rule|RuleDelayed)[k_ /; MatchQ[Unevaluated[k], keyPatt], _], {1}];

(* Unwrap an expression if passed stored inside a single Symbol *)
SetAttributes[derefHeldOnce, HoldAllComplete];
derefHeldOnce[sym_Symbol] /; OwnValues[sym] =!= {} :=
    Replace[First[OwnValues[sym]], (_ :> rhs_) :> derefHeldOnce[rhs]];
derefHeldOnce[expr_] := HoldComplete[expr];

buildLegacyResultObject[maxval_, varrules_List, vars_List, elapsed_, evals_, iters_] :=
    Optimization`GeneralOptimization`TransformationOptimizationResult[<|
        "Objective" -> N[Re[maxval]],
        "Value" -> N[Re[maxval]],
        "Rules" -> varrules,
        "Parameters" -> Association[varrules],
        "ParametersVec" -> Developer`ToPackedArray[N[Re[vars /. varrules]]],
        "VarNames" -> (SymbolName /@ vars),
        "Optimizer" -> "FindMinimum",
        "DerivativeMethod" -> "FiniteDifference",
        "Space" -> "Legacy",
        "CompilationMode" -> "Legacy",
        "Engine" -> "Legacy",
        "Status" -> "Converged",
        "Iterations" -> iters,
        "Evaluations" -> evals,
        "GradientEvaluations" -> 0,
        "HVPEvaluations" -> 0,
        "ElapsedTime" -> elapsed,
        "CompileTime" -> 0.,
        "LegacyFormat" -> {maxval, varrules}
    |>];

(* ------------------------------------------------------------------ *)
(* Legacy FindMinimum engine (all 6 WP 3 defects resolved)            *)
(* ------------------------------------------------------------------ *)
runLegacyFindMinimum[pairHold_HoldComplete, eventsHold_HoldComplete, varHold_Hold,
        x0Vec_List, lbVec_List, ubVec_List, hasConstrQ_?BooleanQ, heldOpts_HoldComplete, retObjQ_?BooleanQ] :=
    Replace[
        {varHold, pairHold, eventsHold},
        {Hold[vs___], HoldComplete[pairExpr_], HoldComplete[eventsExpr_]} :> Block[{vs},
            Module[{nvar = Length[{vs}], taOptsHeld, taBodyHold,
                    negta, methodOpt, fmMethod, allOptList, fmOptsClean, constraints,
                    evalCount = 0, iterCount = 0, t0, minRes, minVal, varrules, maxval, elapsed,
                    hasStepMon},
                (* Filter TransformationAmplitude options safely inside Block without touching RuleDelayed monitors *)
                taOptsHeld = Cases[
                    heldOpts,
                    r:((Rule|RuleDelayed)[k_, _]) /; MemberQ[First /@ Options[TransformationAmplitude], Unevaluated[k]],
                    {1}
                ];
                taBodyHold = Replace[Fold[Join, HoldComplete[], HoldComplete /@ taOptsHeld],
                    HoldComplete[taO___] :> Hold[TransformationAmplitude[pairExpr, eventsExpr, taO]]
                ];

                negta[vArgs:Repeated[_?NumericQ, {nvar}]] := Module[{probeVals, taVal},
                    evalCount++;
                    probeVals = If[hasConstrQ,
                        MapThread[reflectIntoBounds, {{vArgs}, lbVec, ubVec}],
                        {vArgs}
                    ];
                    taVal = evalHeldAtNumeric[varHold, taBodyHold, probeVals];
                    If[NumericQ[taVal], -N[Re[taVal]], $Failed]
                ];

                constraints = If[hasConstrQ,
                    DeleteCases[
                        MapThread[
                            Which[
                                #1 =!= -Infinity && #3 =!= Infinity, #1 <= #2 <= #3,
                                #1 =!= -Infinity, #1 <= #2,
                                #3 =!= Infinity, #2 <= #3,
                                True, True
                            ] &,
                            {lbVec, {vs}, ubVec}
                        ],
                        True
                    ],
                    {}
                ];

                If[constraints === {},
                    Message[MaximizeTransformationAmplitude::unconstr],
                    Message[MaximizeTransformationAmplitude::constr, If[Length[constraints] === 1, First[constraints], And @@ constraints]]
                ];

                methodOpt = extractHeldOptionValue[heldOpts, Method, Automatic];
                fmMethod = If[constraints === {} && methodOpt === Automatic, "InteriorPoint", methodOpt];
                If[!MatchQ[fmMethod, Automatic],
                    Message[MaximizeTransformationAmplitude::meth, fmMethod]
                ];

                allOptList = List @@ heldOpts;
                fmOptsClean = FilterRules[
                    Join[
                        If[fmMethod === Automatic, {}, {Method -> fmMethod}],
                        DeleteCases[allOptList, (Rule|RuleDelayed)[k_ /; MemberQ[{Method, ConstrainPositive, Engine}, Unevaluated[k]], _]]
                    ],
                    Options[FindMinimum]
                ];
                (* Avoid duplicate StepMonitor injection; increment iterCount cleanly when absent *)
                hasStepMon = hasOptionKeyQ[heldOpts, StepMonitor];
                If[!hasStepMon, fmOptsClean = Append[fmOptsClean, StepMonitor :> (iterCount++)]];

                t0 = AbsoluteTime[];
                minRes = Check[
                    If[constraints === {},
                        FindMinimum[negta[vs], Transpose[{{vs}, x0Vec}], Evaluate[Sequence @@ fmOptsClean]],
                        FindMinimum[{negta[vs], constraints}, Transpose[{{vs}, x0Vec}], Evaluate[Sequence @@ fmOptsClean]]
                    ],
                    $Failed,
                    {Event::negdur, Delay::negdur, Pulse::inval\[Beta], Pulse::inval\[Phi], Pulse::neg\[Tau], ShapedPulse::invaldur, FindMinimum::nrnum}
                ];
                elapsed = AbsoluteTime[] - t0;

                If[!MatchQ[minRes, {_?NumericQ, {_Rule ..}}],
                    Message[MaximizeTransformationAmplitude::fail];
                    Return[$Failed, Module]
                ];
                {minVal, varrules} = minRes;
                maxval = -minVal;
                If[retObjQ,
                    buildLegacyResultObject[maxval, varrules, {vs}, elapsed, evalCount, If[hasStepMon, Missing["UserMonitor"], iterCount]],
                    {maxval, varrules}
                ]
            ]
        ]
    ];

(* ------------------------------------------------------------------ *)
(* Public HoldAll entry points                                        *)
(* ------------------------------------------------------------------ *)
MaximizeTransformationAmplitude[HoldPattern[rhoIni_ -> obs_List], more___] := (
    Message[MaximizeTransformationAmplitude::multobs];
    $Failed
);

MaximizeTransformationAmplitude[pairArg_, eventsArg_, varSpecsArg_, opts:(_Rule|_RuleDelayed)...] := Module[
    {pairHold, eventsHold, heldOpts, entries, nP, varHolds, varHold, x0Vec, lbVec, ubVec, explicitFlags,
     cpHits, cpHeld, cpTargetList, cpIndices, engOpt, retObjQ, knownGenOpts, methodOpt, optOpt,
     hasConstrQ, obsCheckFail = False, filteredOptsHeld, canonicalSpecs, prob, compRes,
     compiledSolverOpts},

    pairHold = derefHeldOnce[pairArg];
    eventsHold = derefHeldOnce[eventsArg];
    If[MatchQ[pairHold, HoldComplete[_ -> _List]],
        Message[MaximizeTransformationAmplitude::multobs];
        Return[$Failed]
    ];

    entries = extractHeldVarEntries[varSpecsArg];
    If[entries === $Failed || MemberQ[entries, $Failed] || !ListQ[entries] || Length[entries] === 0,
        Message[MaximizeTransformationAmplitude::fail];
        Return[$Failed]
    ];
    nP = Length[entries];
    varHolds = entries[[All, 1]];
    varHold = Apply[Hold, Join @@ varHolds];
    x0Vec = entries[[All, 2]];
    lbVec = entries[[All, 3, 1]];
    ubVec = entries[[All, 3, 2]];
    explicitFlags = entries[[All, 4]];

    If[!VectorQ[x0Vec, NumericQ] ||
       AnyTrue[Transpose[{lbVec, ubVec}], !((#[[1]] === -Infinity || NumericQ[#[[1]]]) && (#[[2]] === Infinity || NumericQ[#[[2]]]) && TrueQ[#[[1]] < #[[2]]]) &],
        Message[MaximizeTransformationAmplitude::fail];
        Return[$Failed]
    ];

    (* Verify evaluated observable is not a List inside Block without leaking parameter OwnValues *)
    obsCheckFail = Replace[
        {varHold, pairHold},
        {Hold[vs___], HoldComplete[p_]} :> Block[{vs},
            MatchQ[p, _ -> _List]
        ]
    ];
    If[obsCheckFail,
        Message[MaximizeTransformationAmplitude::multobs];
        Return[$Failed]
    ];

    heldOpts = HoldComplete[opts];
    engOpt = extractHeldOptionValue[heldOpts, Engine, Automatic];
    If[!MemberQ[{Automatic, "Compiled", "Legacy"}, engOpt],
        Message[MaximizeTransformationAmplitude::badeng, engOpt];
        Return[$Failed]
    ];

    retObjQ = TrueQ[extractHeldOptionValue[
        heldOpts,
        s_ /; (s === Optimization`GeneralOptimization`ReturnResultObject || (Head[s] === Symbol && SymbolName[Unevaluated[s]] === "ReturnResultObject")),
        False
    ]];

    (* Handle ConstrainPositive cleanly with held parameter symbols *)
    cpHits = Cases[heldOpts, (Rule|RuleDelayed)[k_ /; MatchQ[Unevaluated[k], ConstrainPositive], v_] :> HoldComplete[v], {1}];
    If[cpHits =!= {},
        cpHeld = First[cpHits];
        Which[
            MatchQ[cpHeld, HoldComplete[Automatic | False]],
                Null,
            MatchQ[cpHeld, HoldComplete[True | All]],
                If[MemberQ[explicitFlags, True],
                    Message[MaximizeTransformationAmplitude::constrposconf, ReleaseHold[cpHeld]],
                    Do[lbVec[[i]] = Max[0., lbVec[[i]]], {i, nP}]
                ],
            True,
                cpTargetList = Replace[
                    cpHeld,
                    {
                        HoldComplete[l_List] :> Replace[HoldComplete[l], HoldComplete[{s___}] :> (Hold /@ HoldComplete[s])],
                        HoldComplete[s_] :> {Hold[s]}
                    }
                ];
                cpTargetList = List @@ cpTargetList;
                If[!AllTrue[cpTargetList, MatchQ[#, Hold[_Symbol]] &] || !SubsetQ[varHolds, cpTargetList],
                    Message[MaximizeTransformationAmplitude::badconstrvar, ReleaseHold[cpHeld]];
                    Return[$Failed]
                ];
                cpIndices = Flatten[FirstPosition[varHolds, #] & /@ cpTargetList];
                If[AnyTrue[explicitFlags[[cpIndices]], TrueQ],
                    Message[MaximizeTransformationAmplitude::constrposconf, ReleaseHold[cpHeld]],
                    Do[lbVec[[i]] = Max[0., lbVec[[i]]], {i, cpIndices}]
                ]
        ]
    ];

    hasConstrQ = AnyTrue[Transpose[{lbVec, ubVec}], !(#[[1]] === -Infinity && #[[2]] === Infinity) &];

    If[engOpt === "Legacy",
        Return[
            runLegacyFindMinimum[pairHold, eventsHold, varHold, x0Vec, lbVec, ubVec, hasConstrQ, heldOpts, retObjQ]
        ]
    ];

    (* Compiled / Automatic engine: delegate to Optimization`GeneralOptimization` *)
    knownGenOpts = {
        "LBFGS", "BFGS", "CG", "SteepestAscent", "RegularizedLBFGS",
        "HybridLBFGSNewtonCG", "HybridLBFGSNewtonRFO", "TrustRegionNewtonCG",
        "TrustRegionGLTR", "AdaptiveCubicRegularization",
        "NewtonRFO", "NewtonShift", "NewtonModifiedCholesky", "NewtonSpectral",
        "NelderMead", "Powell", "FindMinimum"
    };
    optOpt = extractHeldOptionValue[heldOpts, OptimalControl`GRAPE`Optimizer, None];
    methodOpt = extractHeldOptionValue[heldOpts, Method, Automatic];

    compiledSolverOpts = {};
    If[optOpt === None && methodOpt =!= Automatic,
        If[MemberQ[knownGenOpts, If[ListQ[methodOpt], First[methodOpt], methodOpt]],
            AppendTo[compiledSolverOpts, OptimalControl`GRAPE`Optimizer -> methodOpt],
            AppendTo[compiledSolverOpts, OptimalControl`GRAPE`Optimizer -> {"FindMinimum", Method -> methodOpt}]
        ]
    ];

    (* Strip MaximizeTransformationAmplitude-specific options already handled above *)
    filteredOptsHeld = Cases[
        heldOpts,
        r:((Rule|RuleDelayed)[k_, _]) /; !MatchQ[Unevaluated[k], ConstrainPositive | Engine | Method |
            Optimization`GeneralOptimization`ReturnResultObject] &&
            !(Head[Unevaluated[k]] === Symbol && SymbolName[Unevaluated[k]] === "ReturnResultObject"),
        {1}
    ];

    (* Build canonical varspecs with resolved lbVec/ubVec and invoke GeneralOptimization inside Block *)
    canonicalSpecs = MapThread[
        Function[{vh, x0, lo, hi},
            Replace[vh, Hold[s_] :>
                If[lo === -Infinity && hi === Infinity,
                    HoldComplete[{s, x0}],
                    HoldComplete[{s, x0, {lo, hi}}]
                ]
            ]
        ],
        {varHolds, x0Vec, lbVec, ubVec}
    ];
    canonicalSpecs = Apply[List, canonicalSpecs];
    With[{
        specsJoined = Replace[Fold[Join, HoldComplete[], canonicalSpecs], HoldComplete[xs___] :> HoldComplete[{xs}]],
        extraOptsHeld = Join[
            HoldComplete @@ compiledSolverOpts,
            With[{rf = retObjQ}, HoldComplete[Optimization`GeneralOptimization`ReturnResultObject -> rf]]
        ]
    },
        With[{allOptsHeld = Join[Fold[Join, HoldComplete[], HoldComplete /@ filteredOptsHeld], extraOptsHeld]},
            compRes = Replace[
                {varHold, pairHold, eventsHold, specsJoined, allOptsHeld},
                {Hold[vs___], HoldComplete[p_], HoldComplete[evs_], HoldComplete[vsSpecs_], HoldComplete[fOpts___]} :>
                    Block[{vs},
                        Optimization`GeneralOptimization`OptimizeTransformationAmplitude[p, evs, vsSpecs, fOpts]
                    ]
            ]
        ]
    ];

    If[compRes =!= $Failed, Return[compRes]];

    If[engOpt === "Compiled",
        Message[MaximizeTransformationAmplitude::fail];
        Return[$Failed]
    ];

    (* Automatic fallback to Legacy FindMinimum when compiled path cannot handle an edge construct *)
    runLegacyFindMinimum[pairHold, eventsHold, varHold, x0Vec, lbVec, ubVec, hasConstrQ, heldOpts, retObjQ]
];

MaximizeTransformationAmplitude::fail="An event duration, or pulse flip angle, has transgressed into an unphysical region, causing MaximizeTransformationAmplitude to fail. Try again with a different starting point, or incorporate constraints on the parameters (see usage message for MaximizeTransformationAmplitude).";

MaximizeTransformationAmplitude::multobs="Only a single target observable may be specified.";

MaximizeTransformationAmplitude::constr="Deploying FindMinimum with constraints `1`.";

MaximizeTransformationAmplitude::unconstr="Deploying FindMinimum without constraints.";

MaximizeTransformationAmplitude::meth="Executing FindMinimum using Method \[Rule] \"`1`\". See the usage message for FindMinimum for other possible values of Method. ";

MaximizeTransformationAmplitude::constrposconf="Explicit constraints take precedence. The directive ConstrainPositive\[Rule]`1` has been ignored.";

MaximizeTransformationAmplitude::badconstrvar="One or more of the constrained variables `1` are inappropriate.";

MaximizeTransformationAmplitude::badeng="Unsupported Engine -> `1`. Expected Automatic, \"Compiled\", or \"Legacy\".";

End[];
EndPackage[];
