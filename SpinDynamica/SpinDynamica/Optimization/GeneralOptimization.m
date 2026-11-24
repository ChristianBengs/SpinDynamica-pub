(* ::Package:: *)

BeginPackage["Optimization`GeneralOptimization`",
{"Global`",
 "SpinDynamica`",
 "General`UserLevel`",
 "General`Matrices`",
 "General`Miscellaneous`",
 "QM`Hilbert`",
 "QM`SpinOperators`",
 "QM`Liouville`",
 "SD`Generators`",
 "SD`Propagation`",
 "Events`Events`",
 "Events`PulseSequence`",
 "Ensemble`EnsembleAverage`",
 "Simulations`Simulations`",
 "OptimalControl`Kernels`",
 "OptimalControl`Objectives`",
 "OptimalControl`GRAPE`",
 "OptimalControl`Optimizers`",
 "OptimalControl`SecondOrderOptimizers`",
 "Optimization`Solvers`",
 "Optimization`Derivatives`"
}];

ClearAll[
    TransformationOptimizationProblem,
    CompileTransformationOptimization,
    EvaluateTransformationOptimization,
    OptimizeTransformationAmplitude,
    TransformationOptimizationResult,
    VerifyTransformationOptimizationGradient,
    VerifyTransformationOptimizationHVP,
    DerivativeMethod,
    SimulationSpace,
    Substeps,
    EnableCaching,
    FuseStaticEvents,
    FiniteDifferenceOrder,
    FiniteDifferenceStep,
    AllowBlackBoxFallback,
    ReturnResultObject
];

TransformationOptimizationProblem::usage = "TransformationOptimizationProblem[\[Rho]ini -> obs, events, varspecs, opts] constructs a symbolic problem specification for general transformation amplitude optimization.\n" <>
"This constructor holds its arguments to prevent global symbol collisions and accepts scalar, vectorial, bounded, or ConstrainPositive parameter specifications along with BackgroundGenerator, EnsembleAverage, Chronology, NormalizationFactor, Basis, and OperatorBasis options.";

CompileTransformationOptimization::usage = "CompileTransformationOptimization[problem, opts] compiles a TransformationOptimizationProblem into packed operator/superoperator atoms, fused static propagators, dependency masks, symbolic derivatives, and cached evaluators.";

EvaluateTransformationOptimization::usage = "EvaluateTransformationOptimization[compiled, xOrRules, opts] evaluates a compiled transformation amplitude optimization problem at a parameter vector or rule list, returning either an objective value or an Association with \"Objective\" and \"Gradient\".";

OptimizeTransformationAmplitude::usage = "OptimizeTransformationAmplitude[problem, opts] or OptimizeTransformationAmplitude[\[Rho]ini -> obs, events, varspecs, opts] maximizes a SpinDynamica transformation amplitude using compiled simulation kernels, analytical or finite-difference derivatives (or derivative-free search), and unified solver backends.\n" <>
"Returns a TransformationOptimizationResult object by default (or {maxval, rules} when ReturnResultObject -> False).";

TransformationOptimizationResult::usage = "TransformationOptimizationResult[...] represents the output of OptimizeTransformationAmplitude, providing property access such as res[\"Value\"], res[\"Rules\"], res[\"Parameters\"], res[\"Iterations\"], res[\"Evaluations\"], res[\"ElapsedTime\"], res[\"DerivativeMethod\"], and res[\"LegacyFormat\"].";

VerifyTransformationOptimizationGradient::usage = "VerifyTransformationOptimizationGradient[probOrCompiled, xOrRules, opts] compares analytical adjoint gradients against high-accuracy central finite differences.";

VerifyTransformationOptimizationHVP::usage = "VerifyTransformationOptimizationHVP[probOrCompiled, xOrRules, dirVec, opts] compares analytical O(K) Hessian-vector products and explicit Hessians against central finite differences of the analytical gradient.";

DerivativeMethod::usage = "DerivativeMethod is an option for OptimizeTransformationAmplitude specifying derivative evaluation strategy: Automatic, \"Analytical\", \"Adjoint\", \"SecondOrder\", \"FiniteDifference\", or \"None\".";
SimulationSpace::usage = "SimulationSpace is an option for general optimization compilation specifying \"Hilbert\", \"Liouville\", or Automatic.";
Substeps::usage = "Substeps is an option specifying the number of fourth-order Magnus (CF4) integration substeps per time-dependent continuous event.";
EnableCaching::usage = "EnableCaching is an option enabling dependency-aware per-stage propagator caching during evaluation.";
FuseStaticEvents::usage = "FuseStaticEvents is an option controlling whether consecutive parameter-independent events are pre-multiplied once at compile time.";
FiniteDifferenceOrder::usage = "FiniteDifferenceOrder is an option specifying \"Central\" (default) or \"Forward\" finite-difference stencils for DerivativeMethod -> \"FiniteDifference\".";
FiniteDifferenceStep::usage = "FiniteDifferenceStep is an option specifying the relative scaled step size epsilon for finite-difference derivative evaluation.";
AllowBlackBoxFallback::usage = "AllowBlackBoxFallback is an option allowing CompileTransformationOptimization to fall back gracefully to TransformationAmplitude wrapping for opaque or non-decomposable event constructs.";
ReturnResultObject::usage = "ReturnResultObject is an option controlling whether optimization routines return a TransformationOptimizationResult (True) or legacy {value, rules} format (False).";

TransformationOptimizationProblem::multobs = "Only a single target observable may be specified.";
TransformationOptimizationProblem::badvars = "Invalid parameter specification `1`. Expected {{p1, p1ini[, bounds]}..}.";
TransformationOptimizationProblem::badconstrvar = "One or more constrained variables `1` are not among the optimization variables.";
TransformationOptimizationProblem::constrposconf = "Explicit bounds take precedence over ConstrainPositive -> `1`.";
CompileTransformationOptimization::compfail = "Failed to compile optimization problem: `1`.";
OptimizeTransformationAmplitude::badsol = "Optimization failed to produce a valid solution (`1`).";
OptimizeTransformationAmplitude::invalp = "Invalid evaluation point `1`.";

SetAttributes[TransformationOptimizationProblem, HoldAll];
SetAttributes[OptimizeTransformationAmplitude, HoldAll];

Options[CompileTransformationOptimization] = {
    SimulationSpace -> Automatic,
    Substeps -> 16,
    EnableCaching -> True,
    FuseStaticEvents -> True,
    AllowBlackBoxFallback -> True,
    Parallel -> False,
    Verbose -> False
};

Options[OptimizeTransformationAmplitude] = Join[
    {
        DerivativeMethod -> Automatic,
        ReturnResultObject -> True,
        SimulationSpace -> Automatic,
        Substeps -> 16,
        EnableCaching -> True,
        FuseStaticEvents -> True,
        AllowBlackBoxFallback -> True,
        FiniteDifferenceOrder -> "Central",
        FiniteDifferenceStep -> Automatic,
        Parallel -> False
    },
    Options[RunGeneralOptimizer]
];

Begin["`Private`"];

toPackedComplexMat[m_] := ToPackedComplexMatrix[m];

(* Safe isolation of user variable symbols inside Held expressions (including HoldAll TimeDependentGeneratorFunction) *)
evalIsolated[varHold_Hold, phSyms_List, bodyHold_Hold] :=
    Replace[
        {varHold, bodyHold},
        {Hold[vs___], Hold[body_]} :> Module[{rules},
            rules = Table[
                With[{pat = Extract[HoldComplete[vs], i, HoldPattern], val = phSyms[[i]]},
                    pat -> val
                ],
                {i, Length[phSyms]}
            ];
            Block[{vs},
                With[{left = HoldComplete[vs]},
                    Do[With[{s = Extract[left, i, HoldComplete]},
                        Replace[s, HoldComplete[sym_] :> (sym = phSyms[[i]])]
                    ], {i, Length[phSyms]}]
                ];
                ReleaseHold[Hold[body] /. rules]
            ]
        ]
    ];

evalAtNumericValues[varHold_Hold, vals_List, bodyHold_Hold] :=
    Replace[
        {varHold, bodyHold},
        {Hold[vs___], Hold[body_]} :> Module[{rules},
            rules = Table[
                With[{pat = Extract[HoldComplete[vs], i, HoldPattern], val = vals[[i]]},
                    pat -> val
                ],
                {i, Length[vals]}
            ];
            Block[{vs},
                With[{left = HoldComplete[vs]},
                    Do[With[{s = Extract[left, i, HoldComplete]},
                        Replace[s, HoldComplete[sym_] :> (sym = vals[[i]])]
                    ], {i, Length[vals]}]
                ];
                ReleaseHold[Hold[body] /. rules]
            ]
        ]
    ];

normalizeVarEntry[HoldComplete[{v_Symbol, v0_}]] := {Hold[v], N[Re[v0]], {-Infinity, Infinity}, False};
normalizeVarEntry[HoldComplete[{v_Symbol, v0_, c_List}]] /; Length[c] === 1 := {Hold[v], N[Re[v0]], {N[Re[First[c]]], Infinity}, True};
normalizeVarEntry[HoldComplete[{v_Symbol, v0_, c_List}]] /; Length[c] === 2 := {Hold[v], N[Re[v0]], {If[c[[1]] === -Infinity || c[[1]] === None, -Infinity, N[Re[c[[1]]]]], If[c[[2]] === Infinity || c[[2]] === None, Infinity, N[Re[c[[2]]]]]}, True};
normalizeVarEntry[HoldComplete[{v_Symbol, v0_, c:Except[_List]}]] := {Hold[v], N[Re[v0]], {N[Re[c]], Infinity}, True};
normalizeVarEntry[_] := $Failed;

SetAttributes[extractHeldVarSpecs, HoldAllComplete];
extractHeldVarSpecs[sym_Symbol] /; OwnValues[sym] =!= {} :=
    Replace[First[OwnValues[sym]], (_ :> rhs_) :> extractHeldVarSpecs[rhs]];
extractHeldVarSpecs[{v:Except[_List], v0:Except[_List]}] := extractHeldVarSpecs[{{v, v0}}];
extractHeldVarSpecs[{v:Except[_List], v0:Except[_List], c_}] := extractHeldVarSpecs[{{v, v0, c}}];
extractHeldVarSpecs[specList_List] := Module[{heldEntries},
    heldEntries = Replace[HoldComplete[specList], HoldComplete[{items___}] :> (HoldComplete /@ HoldComplete[items])];
    List @@ (normalizeVarEntry /@ heldEntries)
];
extractHeldVarSpecs[_] := $Failed;

(* Constructor with HoldAll protecting all symbols inside rho, obs, events, varspecs, opts *)
TransformationOptimizationProblem[pairSym_Symbol, events:Except[_Rule|_RuleDelayed], varspecs:Except[_Rule|_RuleDelayed], opts:(_Rule|_RuleDelayed)...] /;
    MatchQ[OwnValues[pairSym], {(_ :> (_ -> _)) ..}] :=
    Replace[First[OwnValues[pairSym]],
        (_ :> (r_ -> o_)) :> TransformationOptimizationProblem[r -> o, events, varspecs, opts]
    ];

TransformationOptimizationProblem[HoldPattern[rhoIni_ -> obs_], events:Except[_Rule|_RuleDelayed], varspecs:Except[_Rule|_RuleDelayed], opts:(_Rule|_RuleDelayed)...] := Module[
    {entries, nP, varHolds, varNames, x0Vec, boundsRaw, explicitFlags, lbVec, ubVec,
     cpVal, cpList, symNames, badVars, confQ = False},
    If[MatchQ[Unevaluated[obs], _List],
        Message[TransformationOptimizationProblem::multobs];
        Return[$Failed]
    ];
    entries = extractHeldVarSpecs[varspecs];
    If[entries === $Failed || MemberQ[entries, $Failed] || !ListQ[entries] || Length[entries] === 0,
        Message[TransformationOptimizationProblem::badvars, Short[HoldForm[varspecs], 2]];
        Return[$Failed]
    ];
    nP = Length[entries];
    varHolds = entries[[All, 1]];
    varNames = (Function[vh, Replace[vh, Hold[s_] :> SymbolName[Unevaluated[s]]]]) /@ varHolds;
    If[!DuplicateFreeQ[varNames] && Length[DeleteDuplicates[varHolds]] < nP,
        Message[TransformationOptimizationProblem::badvars, "Duplicate parameter symbols."];
        Return[$Failed]
    ];
    x0Vec = entries[[All, 2]];
    If[!VectorQ[x0Vec, NumericQ],
        Message[TransformationOptimizationProblem::badvars, "Non-numeric initial parameter values."];
        Return[$Failed]
    ];
    boundsRaw = entries[[All, 3]];
    explicitFlags = entries[[All, 4]];
    lbVec = boundsRaw[[All, 1]];
    ubVec = boundsRaw[[All, 2]];

    If[AnyTrue[Transpose[{lbVec, ubVec}], !((#[[1]] === -Infinity || NumericQ[#[[1]]]) && (#[[2]] === Infinity || NumericQ[#[[2]]]) && TrueQ[#[[1]] < #[[2]]]) &],
        Message[TransformationOptimizationProblem::badvars, "Invalid bounds (lower bound must be strictly less than upper bound)."];
        Return[$Failed]
    ];

    (* Process ConstrainPositive cleanly with Hold on symbols *)
    cpVal = If[FreeQ[Hold[opts], HoldPattern[ConstrainPositive -> _]],
        Automatic,
        First[Cases[Hold[opts], HoldPattern[(ConstrainPositive -> v_)] :> HoldComplete[v], {1}]]
    ];
    Which[
        cpVal === Automatic || MatchQ[cpVal, HoldComplete[Automatic | False]],
            Null,
        MatchQ[cpVal, HoldComplete[True | All]],
            If[MemberQ[explicitFlags, True],
                confQ = True;
                Message[TransformationOptimizationProblem::constrposconf, ReleaseHold[cpVal]],
                Do[lbVec[[i]] = Max[0., lbVec[[i]]], {i, nP}]
            ],
        True,
            (* Single symbol or list of symbols *)
            cpList = Replace[cpVal, {HoldComplete[l_List] :> Replace[HoldComplete[l], HoldComplete[{s___}] :> (Hold /@ HoldComplete[s])], HoldComplete[s_] :> {Hold[s]}}];
            cpList = List @@ cpList;
            If[!AllTrue[cpList, MatchQ[#, Hold[_Symbol]] &],
                Message[TransformationOptimizationProblem::badconstrvar, ReleaseHold[cpVal]];
                Return[$Failed]
            ];
            If[!SubsetQ[varHolds, cpList],
                Message[TransformationOptimizationProblem::badconstrvar, ReleaseHold[cpVal]];
                Return[$Failed]
            ];
            If[AnyTrue[PositionIndex[varHolds] /@ cpList, explicitFlags[[First[#]]] &],
                confQ = True;
                Message[TransformationOptimizationProblem::constrposconf, ReleaseHold[cpVal]],
                Do[
                    With[{pos = First[FirstPosition[varHolds, cp]]},
                        lbVec[[pos]] = Max[0., lbVec[[pos]]]
                    ],
                    {cp, cpList}
                ]
            ]
    ];

    (* Clamp initial guess minimally if finite? No, keep x0Vec as given so RunGeneralOptimizer validates feasibility! *)
    With[{probAssoc = <|
        "RhoHold" -> Hold[rhoIni],
        "ObsHold" -> Hold[obs],
        "EventsHold" -> Hold[events],
        "VarHold" -> Apply[Hold, Join @@ varHolds],
        "VarHoldsList" -> varHolds,
        "VarNames" -> varNames,
        "X0Vec" -> Developer`ToPackedArray[N[x0Vec]],
        "LBVec" -> lbVec,
        "UBVec" -> ubVec,
        "OptionsHold" -> Hold[opts],
        "SpinSystem" -> SpinSystem[],
        "Basis" -> Basis[],
        "OperatorBasis" -> Quiet[OperatorBasis[]]
    |>},
        TransformationOptimizationProblem[probAssoc]
    ]
];

TransformationOptimizationProblem[expr_TransformationOptimizationProblem] := expr;

(* =============================================================== *)
(* Symbolic Differentiation & Operator-Atom Matrix Decomposer      *)
(* =============================================================== *)
$NonDiffHeads = _If | _Piecewise | _Which | _Switch | _Max | _Min | _Abs | _Sign | _Floor | _Ceiling | _Round | _Mod | _UnitStep | _Derivative | _D;

differentiableExprQ[expr_, phSyms_List] := Module[{ref, ds},
    If[FreeQ[expr, Alternatives @@ phSyms], Return[True]];
    If[!FreeQ[expr, $NonDiffHeads], Return[False]];
    ref = Quiet[Refine[expr, Assumptions -> (Element[#, Reals] & /@ phSyms)]];
    If[!FreeQ[ref, $NonDiffHeads], Return[False]];
    ds = Table[Quiet[D[ref, ph]], {ph, phSyms}];
    FreeQ[ds, _Derivative | _D | $Failed | Indeterminate | ComplexInfinity]
];

$SpinOperatorHeadsPattern =
    _Operator | _Superoperator | _opI | _opT | _opR |
    _UnityOperator | _NullOperator | _NullSuperoperator | _UnitySuperoperator |
    _PopulationOperator | _CoherenceOperator | _ProjectionOperator |
    _ShiftOperator | _SingleTransitionOperator |
    _ScalarOperatorProduct | _OperatorProduct |
    _CommutationSuperoperator | _DoubleCommutationSuperoperator |
    _RotationSuperoperator | _CoherenceOrderFiltrationSuperoperator |
    _PhenomenologicalRelaxationSuperoperator | _ThermalRelaxationSuperoperator |
    _SuperoperatorProduct;

hasSpinObjectQ[expr_] := OperatorQ[expr] || SuperoperatorQ[expr] || !FreeQ[expr, $SpinOperatorHeadsPattern];

(* Recursively collect additive {scalarExpr, opExpr} pairs without expanding MatrixRepresentation prematurely *)
splitScalarOpTerms[expr_, phSyms_List] := Module[{terms, walk},
    walk[0 | 0. | None | NullOperator[] | NullSuperoperator[]] := {};
    walk[e_Plus] := Flatten[walk /@ (List @@ e), 1];
    walk[e_Times] := Module[{facs = List @@ e, ops, scals, opP, scP},
        ops = Select[facs, hasSpinObjectQ];
        scals = Select[facs, !hasSpinObjectQ[#] &];
        If[ops === {},
            {}, (* pure numeric/scalar doesn't represent an operator unless 0 *)
            opP = If[Length[ops] === 1, First[ops], Times @@ ops];
            scP = If[scals === {}, 1., Times @@ scals];
            If[Head[opP] === Plus,
                ({scP * #[[1]], #[[2]]} &) /@ walk[opP],
                {{scP, opP}}
            ]
        ]
    ];
    walk[e_] := If[hasSpinObjectQ[e], {{1., e}}, {}];
    walk[expr]
];

(* Decompose a single symbolic d x d matrix matSym into {{scalExpr_r, M_r}, ...} *)
decomposeSymbolicMatrix[matSym_, phSyms_List, scalarMult_, xTest_List] := Module[
    {refMat, nzPos, nR, nC, atomMap, pos, entry, terms, term, sc, cPart, facList, numFac, symFac,
     rePart, imPart, rAtoms, keys, mats, keyList},
    {nR, nC} = Dimensions[matSym];
    refMat = Quiet[Refine[matSym, Assumptions -> (Element[#, Reals] & /@ phSyms)]];
    If[!FreeQ[refMat, $Failed | Indeterminate], Return[$Failed]];
    nzPos = Position[refMat, Except[0 | 0. | Complex[0., 0.]], {2}, Heads -> False];
    If[nzPos === {}, Return[{}]];

    atomMap = CreateDataStructure["HashTable"];
    addEntry[gKey_, row_Integer, col_Integer, zVal_] := Module[{lst},
        If[Chop[N[zVal]] == 0., Return[]];
        If[atomMap["KeyExistsQ", gKey],
            lst = atomMap["Lookup", gKey];
            AppendTo[lst, {row, col} -> zVal];
            atomMap["Insert", gKey -> lst],
            atomMap["Insert", gKey -> {{row, col} -> zVal}]
        ]
    ];

    Do[
        entry = Expand[refMat[[pos[[1]], pos[[2]]]]];
        If[entry === 0 || entry === 0., Continue[]];
        terms = If[Head[entry] === Plus, List @@ entry, {entry}];
        Do[
            If[FreeQ[term, Alternatives @@ phSyms],
                addEntry[1, pos[[1]], pos[[2]], N[term]],
                (* Extract parameter-free numeric factor vs parameter-dependent part *)
                facList = If[Head[term] === Times, List @@ term, {term}];
                numFac = Times @@ Select[facList, FreeQ[#, Alternatives @@ phSyms] &];
                symFac = Times @@ Select[facList, !FreeQ[#, Alternatives @@ phSyms] &];
                If[!FreeQ[symFac, Complex] || !FreeQ[symFac, I],
                    {rePart, imPart} = Quiet[ComplexExpand[{Re[symFac], Im[symFac]}, TargetFunctions -> {Re, Im}]];
                    rePart = Refine[rePart, Assumptions -> (Element[#, Reals] & /@ phSyms)];
                    imPart = Refine[imPart, Assumptions -> (Element[#, Reals] & /@ phSyms)];
                    If[rePart =!= 0, addEntry[rePart, pos[[1]], pos[[2]], N[numFac]]];
                    If[imPart =!= 0, addEntry[imPart, pos[[1]], pos[[2]], N[I * numFac]]],
                    addEntry[symFac, pos[[1]], pos[[2]], N[numFac]]
                ]
            ],
            {term, terms}
        ],
        {pos, nzPos}
    ];

    keyList = atomMap["Keys"];
    If[keyList === {}, Return[{}]];
    Table[
        With[{rules = atomMap["Lookup", k]},
            (* Sum duplicate positions *)
            With[{sa = Normal[SparseArray[Merge[rules, Total] // Normal, {nR, nC}]]},
                {Expand[scalarMult * If[k === 1, 1., k]], toPackedComplexMat[sa]}
            ]
        ],
        {k, keyList}
    ]
];

(* Hold-safe pure Function builders for scalar, vector, matrix, and HVP evaluators *)
ClearAll[buildPureScalarFn, buildPureVecFn, buildPureMatFn];

buildPureVecFn[exprList_List, phSyms_List, nP_Integer] := With[
    {e = exprList /. Table[phSyms[[a]] -> Hold[Part][Slot[1], a], {a, nP}]},
    Replace[Hold[e] /. Hold[Part] -> Part, Hold[rhs_] :> Function[Developer`ToPackedArray[N[Re[rhs]]]]]
];

buildPureMatFn[exprMat_List, phSyms_List, nP_Integer] := With[
    {e = exprMat /. Table[phSyms[[a]] -> Hold[Part][Slot[1], a], {a, nP}]},
    Replace[Hold[e] /. Hold[Part] -> Part, Hold[rhs_] :> Function[Developer`ToPackedArray[N[Re[rhs]]]]]
];

buildPureScalarFn[expr_, phSyms_List, nP_Integer] := With[
    {e = expr /. Table[phSyms[[a]] -> Hold[Part][Slot[1], a], {a, nP}]},
    Replace[Hold[e] /. Hold[Part] -> Part, Hold[rhs_] :> Function[N[Re[rhs]]]]
];

(* Decompose general operator or superoperator expression genExpr into MatrixAtom spec *)
decomposeGenerator[expr_, phSyms_List, xTest_List, space_String, basis_, opbasis_, forceHermitianQ_?BooleanQ] := Module[
    {pairs, allAtoms = {}, pair, sc, opP, matSym, matAtoms, spec, dim, mTest, mRef, err},
    If[expr === None || expr === NullOperator[] || expr === NullSuperoperator[] || MatchQ[expr, 0 | 0.],
        Return[<|"IsZero" -> True, "IsConstant" -> True, "DepMask" -> ConstantArray[False, Length[phSyms]], "DepCols" -> {}, "Matrices" -> {}, "IsDifferentiable" -> True|>]
    ];

    pairs = splitScalarOpTerms[expr, phSyms];
    If[pairs === {} && hasSpinObjectQ[expr], pairs = {{1., expr}}];
    If[pairs === {}, Return[$Failed]];

    Do[
        {sc, opP} = pair;
        matSym = Quiet[
            If[space === "Hilbert",
                Normal[MatrixRepresentation[opP, basis]],
                Normal[SuperoperatorMatrixRepresentation[opP, opbasis, basis]]
            ]
        ];
        If[!MatrixQ[matSym], Return[$Failed, Module]];
        If[MatrixQ[matSym, NumericQ],
            AppendTo[allAtoms, {sc, toPackedComplexMat[matSym]}],
            (* Parameterized spin operator/superoperator -> symbolic MatrixRepresentation *)
            matAtoms = decomposeSymbolicMatrix[matSym, phSyms, sc, xTest];
            If[matAtoms === $Failed, Return[$Failed, Module]];
            allAtoms = Join[allAtoms, matAtoms]
        ],
        {pair, pairs}
    ];

    spec = buildScalarAtomSpec[allAtoms, phSyms, Length[phSyms], forceHermitianQ];
    If[spec === $Failed, Return[$Failed]];
    If[TrueQ[spec["IsZero"]], Return[spec]];

    (* Exact compile-time verification at test point xTest *)
    mTest = If[spec["IsConstant"],
        spec["Constant"],
        With[{c = spec["EvalCoeffs"][xTest], ms = spec["Matrices"]},
            toPackedComplexMat[Sum[c[[r]] * ms[[r]], {r, Length[ms]}]]
        ]
    ];
    mRef = Quiet[
        With[{evExpr = expr /. Thread[phSyms -> xTest]},
            If[space === "Hilbert",
                Normal[MatrixRepresentation[evExpr, basis]],
                Normal[SuperoperatorMatrixRepresentation[evExpr, opbasis, basis]]
            ]
        ]
    ];
    If[!MatrixQ[mRef, NumericQ],
        mRef = Quiet[
            N[
                (If[space === "Hilbert",
                    Normal[MatrixRepresentation[expr, basis]],
                    Normal[SuperoperatorMatrixRepresentation[expr, opbasis, basis]]
                ]) /. Thread[phSyms -> xTest]
            ]
        ]
    ];
    If[!MatrixQ[mRef, NumericQ] || Dimensions[mRef] =!= Dimensions[mTest], Return[$Failed]];
    err = Max[Abs[Flatten[mTest - mRef]]];
    If[err > 1.*^-9 * Max[1., Max[Abs[Flatten[mRef]]]], Return[$Failed]];
    spec
];

(* Wait: buildScalarAtomSpec needs vArg placeholder matching! Notice lines above used vArg inside buildScalarAtomSpec & buildPureHVPFn -> make vArg global-safe module constant `$vArgSym`! *)
(* Let's also support extracting Liouville vectors for RhoIni and Obs when space === "Liouville" *)
decomposeLiouvilleState[expr_, phSyms_List, xTest_List, opbasis_, basis_, isObs_?BooleanQ] := Module[
    {pairs, allVecAtoms = {}, pair, sc, opP, vecSym, vecAtoms, vTest, vRef, nP = Length[phSyms],
     keyHT, kList, coeffs, vecs, nR, depMask, depCols, actPh, diffQ, cleanCoeffs, jacExpr, hessAll,
     evalCoeffs, evalJac, evalHVP, normSqConst},
    If[FreeQ[expr, Alternatives @@ phSyms] && FreeQ[expr, _SparseArray],
        vRef = Quiet[Normal[OperatorVectorRepresentation[expr, opbasis, basis]]];
        If[!VectorQ[vRef, NumericQ], Return[$Failed]];
        normSqConst = If[isObs, Quiet[N[Re[OperatorNorm[expr]^2]]], 0.];
        Return[<|
            "IsConstant" -> True,
            "Constant" -> Developer`ToPackedArray[N[vRef] + 0. I],
            "UseOperatorNorm" -> isObs,
            "ConstantNormSq" -> normSqConst,
            "DepMask" -> ConstantArray[False, nP],
            "DepCols" -> {},
            "IsDifferentiable" -> True
        |>]
    ];

    pairs = splitScalarOpTerms[expr, phSyms];
    If[pairs === {} && hasSpinObjectQ[expr], pairs = {{1., expr}}];
    If[pairs === {}, Return[$Failed]];
    Do[
        {sc, opP} = pair;
        vecSym = Quiet[Normal[OperatorVectorRepresentation[opP, opbasis, basis]]];
        If[!VectorQ[vecSym], Return[$Failed, Module]];
        If[VectorQ[vecSym, NumericQ],
            AppendTo[allVecAtoms, {sc, Developer`ToPackedArray[N[vecSym] + 0. I]}],
            (* Symbolic Op in state: decompose symbolic OperatorVectorRepresentation *)
            vecAtoms = decomposeSymbolicMatrix[{vecSym}, phSyms, sc, xTest];
            If[vecAtoms === $Failed, Return[$Failed, Module]];
            allVecAtoms = Join[allVecAtoms, ({#[[1]], First[#[[2]]]} &) /@ vecAtoms]
        ],
        {pair, pairs}
    ];
    keyHT = CreateDataStructure["HashTable"];
    Do[
        With[{s = Expand[Refine[a[[1]], Assumptions -> (Element[#, Reals] & /@ phSyms)]], v = a[[2]]},
            If[keyHT["KeyExistsQ", s],
                keyHT["Insert", s -> Developer`ToPackedArray[keyHT["Lookup", s] + v]],
                keyHT["Insert", s -> v]
            ]
        ],
        {a, allVecAtoms}
    ];
    coeffs = keyHT["Keys"];
    vecs = keyHT["Lookup", #] & /@ coeffs;
    nR = Length[vecs];
    depMask = Table[!FreeQ[coeffs, phSyms[[a]]], {a, nP}];
    depCols = Flatten[Position[depMask, True]];
    actPh = phSyms[[depCols]];
    diffQ = And @@ (differentiableExprQ[#, actPh] & /@ coeffs);
    cleanCoeffs = Refine[coeffs, Assumptions -> (Element[#, Reals] & /@ actPh)];
    evalCoeffs = buildPureVecFn[cleanCoeffs, phSyms, nP];
    evalJac = If[diffQ,
        jacExpr = Table[Refine[D[cleanCoeffs[[r]], ph], Assumptions -> (Element[#, Reals] & /@ actPh)], {ph, actPh}, {r, nR}];
        If[!FreeQ[jacExpr, _Derivative | _D], diffQ = False; None, buildPureMatFn[jacExpr, phSyms, nP]],
        None
    ];
    vTest = Sum[evalCoeffs[xTest][[r]] * vecs[[r]], {r, nR}];
    vRef = Quiet[Normal[OperatorVectorRepresentation[expr /. Thread[phSyms -> xTest], opbasis, basis]]];
    If[!VectorQ[vRef, NumericQ],
        vRef = Quiet[N[Normal[OperatorVectorRepresentation[expr, opbasis, basis]] /. Thread[phSyms -> xTest]]]
    ];
    If[!VectorQ[vRef, NumericQ] || Max[Abs[vTest - vRef]] > 1.*^-9 * Max[1., Max[Abs[vRef]]],
        Return[$Failed]
    ];
    (* Check if Norm[vTest]^2 equals OperatorNorm[obs]^2 for dynamic Obs *)
    With[{opN = If[isObs, Quiet[N[Re[OperatorNorm[expr /. Thread[phSyms -> xTest]]^2]]], 0.]},
        If[isObs && Abs[opN - Re[Conjugate[vTest] . vTest]] > 1.*^-8 * Max[1., opN],
            Return[$Failed]
        ]
    ];
    <|
        "IsConstant" -> False,
        "Vectors" -> vecs,
        "UseOperatorNorm" -> False,
        "DepMask" -> depMask,
        "DepCols" -> depCols,
        "EvalCoeffs" -> evalCoeffs,
        "EvalCoeffJac" -> evalJac,
        "IsDifferentiable" -> diffQ
    |>
];

(* Build scalar spec with explicit Hessian for CF4 substeps and HVP *)
buildScalarAtomSpec[atomPairs_List, phSyms_List, nP_Integer, forceHermitianQ_?BooleanQ] := Module[
    {keyHT, kList, coeffs, mats, nR, depMask, depCols, actPh, diffQ,
     cleanCoeffs, jacExpr, hessAll, evalCoeffs, evalJac, evalHVP, xSym, vSym, hEx},
    If[atomPairs === {},
        Return[<|"IsZero" -> True, "IsConstant" -> True, "DepMask" -> ConstantArray[False, nP], "DepCols" -> {}, "Matrices" -> {}, "IsDifferentiable" -> True|>]
    ];
    keyHT = CreateDataStructure["HashTable"];
    Do[
        With[{sc = Expand[Refine[pair[[1]], Assumptions -> (Element[#, Reals] & /@ phSyms)]], mat = pair[[2]]},
            If[keyHT["KeyExistsQ", sc],
                keyHT["Insert", sc -> toPackedComplexMat[keyHT["Lookup", sc] + mat]],
                keyHT["Insert", sc -> toPackedComplexMat[mat]]
            ]
        ],
        {pair, atomPairs}
    ];
    kList = keyHT["Keys"];
    coeffs = kList;
    mats = keyHT["Lookup", #] & /@ kList;
    If[TrueQ[forceHermitianQ],
        mats = toPackedComplexMat[0.5 * (# + ConjugateTranspose[#])] & /@ mats
    ];
    With[{nzIdx = Pick[Range[Length[mats]], (Max[Abs[Flatten[#]]] > 0. &) /@ mats]},
        coeffs = coeffs[[nzIdx]];
        mats = mats[[nzIdx]]
    ];
    nR = Length[mats];
    If[nR === 0,
        Return[<|"IsZero" -> True, "IsConstant" -> True, "DepMask" -> ConstantArray[False, nP], "DepCols" -> {}, "Matrices" -> {}, "IsDifferentiable" -> True|>]
    ];

    depMask = Table[!FreeQ[coeffs, phSyms[[a]]], {a, nP}];
    depCols = Flatten[Position[depMask, True]];
    If[depCols === {},
        With[{constM = toPackedComplexMat[Sum[N[coeffs[[r]]] * mats[[r]], {r, nR}]]},
            Return[<|
                "IsZero" -> False,
                "IsConstant" -> True,
                "Constant" -> constM,
                "DepMask" -> depMask,
                "DepCols" -> {},
                "Matrices" -> {constM},
                "IsDifferentiable" -> True
            |>]
        ]
    ];

    actPh = phSyms[[depCols]];
    diffQ = And @@ (differentiableExprQ[#, actPh] & /@ coeffs);
    cleanCoeffs = Refine[coeffs, Assumptions -> (Element[#, Reals] & /@ actPh)];
    evalCoeffs = buildPureVecFn[cleanCoeffs, phSyms, nP];
    evalJac = None; evalHVP = None;

    If[TrueQ[diffQ],
        jacExpr = Table[Refine[D[cleanCoeffs[[r]], ph], Assumptions -> (Element[#, Reals] & /@ actPh)], {ph, actPh}, {r, nR}];
        If[!FreeQ[jacExpr, _Derivative | _D | $Failed | Indeterminate],
            diffQ = False,
            evalJac = buildPureMatFn[jacExpr, phSyms, nP];
            hessAll = Table[
                Total[Table[
                    With[{d2 = Refine[D[jacExpr[[idxA, r]], actPh[[idxB]]], Assumptions -> (Element[#, Reals] & /@ actPh)]},
                        If[d2 === 0 || d2 === 0., 0., Hold[Part][Slot[2], idxB] * d2]
                    ],
                    {idxB, Length[actPh]}, {r, nR}
                ]],
                {idxA, Length[actPh]}
            ];
            evalHVP = If[AllTrue[Flatten[hessAll], # === 0 || # === 0. &],
                With[{na = Length[actPh], nr = nR}, Function[{x, v}, ConstantArray[0., {na, nr}]]],
                With[{e = hessAll /. Table[phSyms[[a]] -> Hold[Part][Slot[1], a], {a, nP}]},
                    Replace[Hold[e] /. Hold[Part] -> Part, Hold[rhs_] :> Function[Developer`ToPackedArray[N[Re[rhs]]]]]
                ]
            ]
        ]
    ];

    <|
        "IsZero" -> False,
        "IsConstant" -> False,
        "DepMask" -> depMask,
        "DepCols" -> depCols,
        "Matrices" -> mats,
        "EvalCoeffs" -> evalCoeffs,
        "EvalCoeffJac" -> evalJac,
        "EvalCoeffHVP" -> evalHVP,
        "IsDifferentiable" -> diffQ,
        "NumAtoms" -> nR
    |>
];

(* Expand instant rotation spec into list of {generatorExpr, factor1} *)
rotationToGenerators[spins_, ang_] := Switch[ang,
    {_?(!ListQ[#] &), "x"}, {{ang[[1]] * opI[spins, "x"]}},
    {_?(!ListQ[#] &), "y"}, {{ang[[1]] * opI[spins, "y"]}},
    {_?(!ListQ[#] &), "-x"}, {{-ang[[1]] * opI[spins, "x"]}},
    {_?(!ListQ[#] &), "-y"}, {{-ang[[1]] * opI[spins, "y"]}},
    {_?(!ListQ[#] &), "z"}, {{ang[[1]] * opI[spins, "z"]}},
    {_?(!ListQ[#] &), _?(!ListQ[#] &)},
        {{(ang[[1]] * Cos[ang[[2]]]) * opI[spins, "x"] + (ang[[1]] * Sin[ang[[2]]]) * opI[spins, "y"]}},
    {_, {_, _}},
        With[{xi = ang[[1]], th = ang[[2, 1]], ph = ang[[2, 2]]},
            {{xi * (Sin[th] Cos[ph] opI[spins, "x"] + Sin[th] Sin[ph] opI[spins, "y"] + Cos[th] opI[spins, "z"])}}
        ],
    {_, _, _},
        (* Euler ZYZ: alpha (3rd), beta (2nd), gamma (1st) chronologically *)
        {
            {ang[[3]] * opI[spins, "z"]},
            {ang[[2]] * opI[spins, "y"]},
            {ang[[1]] * opI[spins, "z"]}
        },
    _, $Failed
];

unwrapRotationSuperoperator[sop_, space_String, basis_, opbasis_, phSyms_List, xTest_List] := Module[
    {rotSpecs, genItems, stList = {}},
    (* Extract list of {spins, angleSpec} pairs in chronological application order without re-evaluating RotationSuperoperator *)
    rotSpecs = Which[
        MatchQ[sop, HoldPattern[RotationSuperoperator[spins_, {phi_, "z"}] . RotationSuperoperator[spins_, {b_, "x"}] . RotationSuperoperator[spins_, {negPhi_, "z"}]]] &&
            Quiet[TrueQ[Chop[Simplify[(sop[[1, 2, 1]] + sop[[3, 2, 1]]) /. Thread[phSyms -> xTest]]] == 0]],
            {{sop[[1, 1]], {sop[[2, 2, 1]], sop[[1, 2, 1]]}}},
        Head[sop] === RotationSuperoperator,
            {{sop[[1]], sop[[2]]}},
        Head[sop] === Dot && AllTrue[List @@ sop, Head[#] === RotationSuperoperator &],
            ({#[[1]], #[[2]]} &) /@ Reverse[List @@ sop], (* RightToLeft composition -> Reverse for LeftToRight *)
        True,
            Return[$Failed]
    ];
    Do[
        With[{rSpins = rSpec[[1]], rAng = rSpec[[2]]},
            If[FreeQ[{rSpins, rAng}, Alternatives @@ phSyms],
                With[{m = Quiet[If[space === "Hilbert",
                        toPackedComplexMat[MatrixRepresentation[RotationOperator[rSpins, rAng], basis]],
                        toPackedComplexMat[SuperoperatorMatrixRepresentation[RotationSuperoperator[rSpins, rAng], opbasis, basis]]
                    ]]},
                    If[!MatrixQ[m, NumericQ], Return[$Failed, Module]];
                    AppendTo[stList, <|
                        "Kind" -> If[space === "Hilbert", "Static", "LiouvilleStatic"],
                        "DepMask" -> ConstantArray[False, Length[phSyms]],
                        "DepCols" -> {},
                        "Propagator" -> m,
                        "IsDifferentiable" -> True,
                        "SupportsSecondOrder" -> True
                    |>]
                ],
                genItems = rotationToGenerators[rSpins, rAng];
                If[genItems === $Failed, Return[$Failed, Module]];
                Do[
                    With[{genExpr = If[space === "Hilbert", g[[1]], CommutationSuperoperator[-I * g[[1]]]]},
                        With[{spec = decomposeGenerator[genExpr, phSyms, xTest, space, basis, opbasis, space === "Hilbert"]},
                            If[spec === $Failed, Return[$Failed, Module]];
                            If[!TrueQ[spec["IsZero"]],
                                AppendTo[stList,
                                    If[spec["IsConstant"],
                                        <|"Kind" -> If[space === "Hilbert", "Static", "LiouvilleStatic"],
                                          "DepMask" -> spec["DepMask"],
                                          "DepCols" -> {},
                                          "Propagator" -> toPackedComplexMat[MatrixExp[If[space === "Hilbert", -I, 1.] * spec["Constant"]]],
                                          "IsDifferentiable" -> True,
                                          "SupportsSecondOrder" -> True|>,
                                        <|"Kind" -> If[space === "Hilbert", "HilbertExpm", "LiouvilleExpm"],
                                          "DepMask" -> spec["DepMask"],
                                          "DepCols" -> spec["DepCols"],
                                          "Matrices" -> spec["Matrices"],
                                          "EvalCoeffs" -> spec["EvalCoeffs"],
                                          "EvalCoeffJac" -> spec["EvalCoeffJac"],
                                          "EvalCoeffHVP" -> spec["EvalCoeffHVP"],
                                          "IsScalarTimeGenerator" -> False,
                                          "IsDifferentiable" -> spec["IsDifferentiable"],
                                          "SupportsSecondOrder" -> spec["IsDifferentiable"]|>
                                    ]
                                ]
                            ]
                        ]
                    ],
                    {g, genItems}
                ]
            ]
        ],
        {rSpec, rotSpecs}
    ];
    stList
];

(* Helper to convert SpinDynamica time-dependent generator + background into scalar function of tSym & phSyms *)
extractTDGenBody[gen_, tSym_Symbol, tStartExpr_, durExpr_] := Which[
    gen === None || gen === NullOperator[] || gen === NullSuperoperator[] || MatchQ[gen, 0 | 0.],
        0,
    Head[gen] === TimeDependentGeneratorFunction,
        Module[{tAbsPatt, nPairs, pairRules, exprHold},
            tAbsPatt = Extract[gen, 1, HoldPattern];
            nPairs = Extract[gen, 2, Function[ps, Length[Unevaluated[ps]], HoldAll]];
            pairRules = Table[
                With[{lhsPatt = Extract[gen, {2, k, 1}, HoldPattern], frac = Extract[gen, {2, k, 2}]},
                    lhsPatt -> (tSym - (tStartExpr + frac * durExpr))
                ],
                {k, nPairs}
            ];
            exprHold = Extract[gen, 4, Hold];
            ReleaseHold[exprHold /. Prepend[pairRules, tAbsPatt -> tSym]]
        ],
    Head[gen] === PeriodicFunction,
        With[{tPatt = Extract[gen, 1, HoldPattern], exprHold = Extract[gen, 3, Hold]},
            ReleaseHold[exprHold /. (tPatt -> tSym)]
        ],
    Head[gen] === Function,
        If[Length[Unevaluated[gen]] >= 2 && MatchQ[Extract[gen, 1, Hold], Hold[_Symbol]],
            With[{tPatt = Extract[gen, 1, HoldPattern], exprHold = Extract[gen, 2, Hold]},
                ReleaseHold[exprHold /. (tPatt -> tSym)]
            ],
            gen[tSym]
        ],
    True, gen
];

containsTDGeneratorQ[g_] := !FreeQ[g, _TimeDependentGeneratorFunction | _PeriodicFunction | _Function];

(* Compile a CF4 continuous TD stage cleanly *)
compileTDStage[evGen_, backgGen_, tStartExpr_, durExpr_, phSyms_List, xTest_List, space_String, basis_, opbasis_, nSub_Integer] := Module[
    {tSym = Unique["gt"], tS = Unique["gts"], tD = Unique["gd"], combExpr, augSyms, augTest, pairs, allAtoms = {},
     pair, sc, opP, matSym, matAtoms, keyHT, kList, coeffs, mats, nR, depMask, depCols, actPh, diffQ,
     cleanCoeffs, cTimeExpr, dCdTExpr, jParamExpr, dDurExpr, dTStartExpr, evalHamCoeffs, evalHamTimeDeriv,
     evalHamParamJac, evalDuration, evalDurationGrad, evalTStart, evalTStartGrad, evalSubHess,
     xSym = Unique["gx"], vSym = Unique["gv"], tArg = Unique["gtArg"],
     f1S = Unique["gf1"], f2S = Unique["gf2"], w1S = Unique["gw1"], w2S = Unique["gw2"],
     dim = If[space === "Hilbert", Length[BasisKets[basis]], Length[BasisOperators[opbasis]]]},

    combExpr = If[space === "Hilbert",
        extractTDGenBody[evGen, tSym, tStartExpr, durExpr] + extractTDGenBody[backgGen, tSym, tStartExpr, durExpr],
        (* For Liouville space, note TD generator might be Operator -> wrap CommutationSuperoperator[-I ...] *)
        With[{evPart = extractTDGenBody[evGen, tSym, tStartExpr, durExpr], bgPart = extractTDGenBody[backgGen, tSym, tStartExpr, durExpr]},
            If[OperatorQ[evPart], CommutationSuperoperator[-I * evPart], evPart] +
            If[OperatorQ[bgPart], CommutationSuperoperator[-I * bgPart], bgPart]
        ]
    ];

    augSyms = Append[phSyms, tSym];
    augTest = Append[xTest, 0.5 * (tStartExpr + 0.5 * durExpr /. Thread[phSyms -> xTest])];
    pairs = splitScalarOpTerms[combExpr, augSyms];
    If[pairs === {} && hasSpinObjectQ[combExpr], pairs = {{1., combExpr}}];
    If[pairs === {}, Return[$Failed]];

    Do[
        {sc, opP} = pair;
        matSym = Quiet[
            If[space === "Hilbert",
                Normal[MatrixRepresentation[opP, basis]],
                Normal[SuperoperatorMatrixRepresentation[opP, opbasis, basis]]
            ]
        ];
        If[!MatrixQ[matSym], Return[$Failed, Module]];
        If[MatrixQ[matSym, NumericQ],
            AppendTo[allAtoms, {sc, toPackedComplexMat[matSym]}],
            matAtoms = decomposeSymbolicMatrix[matSym, augSyms, sc, augTest];
            If[matAtoms === $Failed, Return[$Failed, Module]];
            allAtoms = Join[allAtoms, matAtoms]
        ],
        {pair, pairs}
    ];

    keyHT = CreateDataStructure["HashTable"];
    Do[
        With[{s = Expand[Refine[a[[1]], Assumptions -> (Element[#, Reals] & /@ augSyms)]], m = a[[2]]},
            If[keyHT["KeyExistsQ", s],
                keyHT["Insert", s -> toPackedComplexMat[keyHT["Lookup", s] + m]],
                keyHT["Insert", s -> toPackedComplexMat[m]]
            ]
        ],
        {a, allAtoms}
    ];
    kList = keyHT["Keys"];
    coeffs = kList;
    mats = keyHT["Lookup", #] & /@ kList;
    If[space === "Hilbert", mats = toPackedComplexMat[0.5 * (# + ConjugateTranspose[#])] & /@ mats];
    With[{nz = Pick[Range[Length[mats]], (Max[Abs[Flatten[#]]] > 0. &) /@ mats]},
        coeffs = coeffs[[nz]]; mats = mats[[nz]]
    ];
    nR = Length[mats];
    If[nR === 0,
        Return[<|"Kind" -> If[space === "Hilbert", "Static", "LiouvilleStatic"],
                 "DepMask" -> ConstantArray[False, Length[phSyms]], "DepCols" -> {},
                 "Propagator" -> toPackedComplexMat[IdentityMatrix[dim]],
                 "IsDifferentiable" -> True, "SupportsSecondOrder" -> True|>]
    ];

    (* Dependency includes BOTH explicit Hamiltonian parameters AND duration / tStart parameters! *)
    depMask = Table[
        !FreeQ[{coeffs, durExpr, tStartExpr}, phSyms[[a]]],
        {a, Length[phSyms]}
    ];
    depCols = Flatten[Position[depMask, True]];
    actPh = phSyms[[depCols]];

    diffQ = differentiableExprQ[durExpr, actPh] && differentiableExprQ[tStartExpr, actPh] &&
        And @@ (differentiableExprQ[#, Append[actPh, tSym]] & /@ coeffs);

    cleanCoeffs = Refine[coeffs, Assumptions -> (Element[#, Reals] & /@ Append[actPh, tSym])];
    With[{subAll = Append[Table[phSyms[[a]] -> Hold[Part][Slot[1], a], {a, Length[phSyms]}], tSym -> Slot[2]]},
        With[{cEx = cleanCoeffs /. subAll, tsEx = tStartExpr /. Most[subAll], dEx = durExpr /. Most[subAll]},
            evalHamCoeffs = Replace[Hold[cEx] /. Hold[Part] -> Part,
                Hold[rhs_] :> Function[Developer`ToPackedArray[N[Re[rhs]]]]];
            evalTStart = Replace[Hold[tsEx] /. Hold[Part] -> Part,
                Hold[rhs_] :> Function[N[Re[rhs]]]];
            evalDuration = Replace[Hold[dEx] /. Hold[Part] -> Part,
                Hold[rhs_] :> Function[N[Re[rhs]]]]
        ]
    ];

    evalHamTimeDeriv = None; evalHamParamJac = None; evalDurationGrad = None; evalTStartGrad = None; evalSubHess = None;
    If[TrueQ[diffQ],
        dCdTExpr = Table[Refine[D[cleanCoeffs[[r]], tSym], Assumptions -> (Element[#, Reals] & /@ Append[actPh, tSym])], {r, nR}];
        jParamExpr = Table[Refine[D[cleanCoeffs[[r]], ph], Assumptions -> (Element[#, Reals] & /@ Append[actPh, tSym])], {ph, actPh}, {r, nR}];
        dDurExpr = Table[Refine[D[durExpr, ph], Assumptions -> (Element[#, Reals] & /@ actPh)], {ph, actPh}];
        dTStartExpr = Table[Refine[D[tStartExpr, ph], Assumptions -> (Element[#, Reals] & /@ actPh)], {ph, actPh}];
        If[!FreeQ[{dCdTExpr, jParamExpr, dDurExpr, dTStartExpr}, _Derivative | _D],
            diffQ = False,
            With[{subAll = Append[Table[phSyms[[a]] -> Hold[Part][Slot[1], a], {a, Length[phSyms]}], tSym -> Slot[2]]},
                With[{dtEx = dCdTExpr /. subAll, jpEx = jParamExpr /. subAll, ddEx = dDurExpr /. Most[subAll], tsEx = dTStartExpr /. Most[subAll]},
                    evalHamTimeDeriv = Replace[Hold[dtEx] /. Hold[Part] -> Part,
                        Hold[rhs_] :> Function[Developer`ToPackedArray[N[Re[rhs]]]]];
                    evalHamParamJac = Replace[Hold[jpEx] /. Hold[Part] -> Part,
                        Hold[rhs_] :> Function[Developer`ToPackedArray[N[Re[rhs]]]]];
                    evalDurationGrad = Replace[Hold[ddEx] /. Hold[Part] -> Part,
                        Hold[rhs_] :> Function[Developer`ToPackedArray[N[Re[rhs]]]]];
                    evalTStartGrad = Replace[Hold[tsEx] /. Hold[Part] -> Part,
                        Hold[rhs_] :> Function[Developer`ToPackedArray[N[Re[rhs]]]]]
                ]
            ];
            (* Exact 2nd-order HVP generator for CF4 substep (only built for Hilbert space!) *)
            evalSubHess = If[space === "Hilbert" && Length[depCols] > 0,
                Module[{cNode1, cNode2, cmGeneric, h2Tensor, h2Fn, h2Ex, subH2},
                    cNode1 = cleanCoeffs /. tSym -> (tStartExpr + f1S * durExpr);
                    cNode2 = cleanCoeffs /. tSym -> (tStartExpr + f2S * durExpr);
                    cmGeneric = (durExpr / nSub) * (w1S * cNode1 + w2S * cNode2);
                    h2Tensor = Table[
                        Refine[D[cmGeneric[[r]], actPh[[idxA]], actPh[[idxB]]], Assumptions -> (Element[#, Reals] & /@ actPh)],
                        {idxA, Length[actPh]}, {idxB, Length[actPh]}, {r, nR}
                    ];
                    If[!FreeQ[h2Tensor, _Derivative | _D],
                        None,
                        If[!AnyTrue[Flatten[h2Tensor], # =!= 0 && # =!= 0. &],
                            With[{na = Length[actPh], nr = nR},
                                Function[{mIdx, kIdx, isSec, w1Val, w2Val, f1Val, f2Val, xVecIn},
                                    Function[{xDummy, vActIn}, ConstantArray[0., {na, nr}]]
                                ]
                            ],
                            subH2 = Join[{f1S -> Slot[1], f2S -> Slot[2], w1S -> Slot[3], w2S -> Slot[4]}, Table[phSyms[[a]] -> Hold[Part][Slot[5], a], {a, Length[phSyms]}]];
                            With[{e = h2Tensor /. subH2},
                                h2Fn = Replace[Hold[e] /. Hold[Part] -> Part,
                                    Hold[rhs_] :> Function[Developer`ToPackedArray[N[Re[rhs]]]]];
                                With[{hf = h2Fn},
                                    Function[{mIdx, kIdx, isSec, w1Val, w2Val, f1Val, f2Val, xVecIn},
                                        With[{tens = hf[f1Val, f2Val, w1Val, w2Val, xVecIn]},
                                            Function[{xDummy, vActIn}, Developer`ToPackedArray[Table[vActIn . tens[[idxA]], {idxA, Length[vActIn]}]]]
                                        ]
                                    ]
                                ]
                            ]
                        ]
                    ]
                ],
                None
            ]
        ]
    ];

    <|
        "Kind" -> If[space === "Hilbert", "HilbertCF4", "LiouvilleCF4"],
        "Dim" -> dim,
        "DepMask" -> depMask,
        "DepCols" -> depCols,
        "NumSubsteps" -> nSub,
        "Matrices" -> mats,
        "EvalTStart" -> evalTStart,
        "EvalDuration" -> evalDuration,
        "EvalHamCoeffs" -> evalHamCoeffs,
        "EvalHamTimeDeriv" -> evalHamTimeDeriv,
        "EvalHamParamJac" -> evalHamParamJac,
        "EvalDurationGrad" -> evalDurationGrad,
        "EvalTStartGrad" -> evalTStartGrad,
        "EvalSubstepHessianFn" -> evalSubHess,
        "IsDifferentiable" -> diffQ,
        "SupportsSecondOrder" -> (space === "Hilbert" && diffQ && (depCols === {} || evalSubHess =!= None))
    |>
];

(* Compile one ensemble member (or the single standalone member) *)
compileSingleMember[rhoExpr_, obsExpr_, evspecs_, backgGen_, chron_, tIniSpec_, normExpr_, phSyms_List, xTest_List,
        spaceOpt_, basis_, opbasis_, nSub_Integer, fuseQ_?BooleanQ] := Module[
    {evListObj, eventsChrono, space, nP = Length[phSyms], stages = {}, fusedStages = {},
     ev, timing, durExpr, tbExpr, taExpr, genList, genExpr, evOpts, hasOwnBg, bgActive,
     spec, rotStages, cofMat, rhoSpec, obsSpec, normSpec, isDiff = True, supp2nd = True,
     expExpr, gSpec, idCounter = 0, stPrev, propProd, dim, tAccExpr},

    (* Standardize event specifications into SpinDynamica EventList *)
    evListObj = Quiet[
        TimeConstrained[
            EventList[
                EventSpecifications[evspecs, Chronology -> chron],
                If[tIniSpec === Automatic || tIniSpec === None || tIniSpec === Undefined,
                    InitialTimePoint -> 0,
                    InitialTimePoint -> tIniSpec
                ]
            ],
            5, $Failed
        ]
    ];
    If[!MatchQ[evListObj, _EventList], Return[$Failed]];
    eventsChrono = Reverse[evListObj[[3]]]; (* Chronological Left-to-Right execution order *)

    (* Determine target space automatically unless forced *)
    space = Which[
        spaceOpt === "Hilbert" || spaceOpt === "Liouville", spaceOpt,
        True,
            If[!FreeQ[{eventsChrono, backgGen, rhoExpr, obsExpr},
                    _Superoperator | _CommutationSuperoperator | _DoubleCommutationSuperoperator |
                    _PhenomenologicalRelaxationSuperoperator | _ThermalRelaxationSuperoperator |
                    _CoherenceOrderFiltrationSuperoperator | _CoherenceOrderFiltration] ||
               !TrueQ[DefinedBasisQ[basis]] ||
               (opbasis =!= Quiet[ZeemanKetBraOperatorBasis[SpinSystem[basis], Sorted -> "Canonical"]] &&
                !FreeQ[{eventsChrono, backgGen}, _PhenomenologicalRelaxationSuperoperator | _Superoperator]),
                "Liouville",
                "Hilbert"
            ]
    ];
    dim = If[space === "Hilbert", Length[BasisKets[basis]], Length[BasisOperators[opbasis]]];

    tAccExpr = If[NumericQ[tIniSpec], tIniSpec, 0];
    Do[
        If[!MatchQ[ev, _Event], Return[$Failed, Module]];
        timing = ev[[2, 1]];
        durExpr = timing[[1]];
        tbExpr = timing[[2, 1]];
        taExpr = timing[[2, 2]];
        If[taExpr === Undefined, taExpr = tAccExpr];
        genList = ev[[3]];
        evOpts = If[Length[ev] >= 4, List @@ ev[[4 ;;]], {}];
        bgActive = Lookup[evOpts, BackgroundGenerator, backgGen];
        If[bgActive === False, bgActive = None];

        Which[
            (* Instantaneous superoperator / CoherenceOrderFiltration / RotationSuperoperator *)
            MatchQ[genList, {_CoherenceOrderFiltrationSuperoperator, ___}] ||
            MatchQ[genList, {_?SuperoperatorQ}] ||
            (Length[genList] === 1 && Head[genList[[1]]] === Dot && AllTrue[List @@ genList[[1]], Head[#] === RotationSuperoperator &]),
                With[{sop = genList[[1]]},
                    Which[
                        Head[sop] === CoherenceOrderFiltrationSuperoperator,
                            If[space =!= "Liouville" || !FreeQ[sop, Alternatives @@ phSyms], Return[$Failed, Module]];
                            cofMat = Quiet[toPackedComplexMat[SuperoperatorMatrixRepresentation[sop, opbasis, basis]]];
                            If[!MatrixQ[cofMat, NumericQ], Return[$Failed, Module]];
                            AppendTo[stages, <|
                                "Kind" -> "LiouvilleStatic",
                                "DepMask" -> ConstantArray[False, nP],
                                "DepCols" -> {},
                                "Propagator" -> cofMat,
                                "IsDifferentiable" -> True,
                                "SupportsSecondOrder" -> False
                            |>],
                        Head[sop] === RotationSuperoperator || (Head[sop] === Dot && AllTrue[List @@ sop, Head[#] === RotationSuperoperator &]),
                            rotStages = unwrapRotationSuperoperator[sop, space, basis, opbasis, phSyms, xTest];
                            If[rotStages === $Failed, Return[$Failed, Module]];
                            stages = Join[stages, rotStages],
                        MatchQ[sop, _UnitySuperoperator],
                            Null,
                        True,
                            (* Static or general instantaneous superoperator *)
                            If[space =!= "Liouville", Return[$Failed, Module]];
                            cofMat = Quiet[toPackedComplexMat[Normal[SuperoperatorMatrixRepresentation[sop, opbasis, basis]]]];
                            If[!MatrixQ[cofMat, NumericQ], Return[$Failed, Module]];
                            AppendTo[stages, <|
                                "Kind" -> "LiouvilleStatic",
                                "DepMask" -> ConstantArray[False, nP],
                                "DepCols" -> {},
                                "Propagator" -> cofMat,
                                "IsDifferentiable" -> True,
                                "SupportsSecondOrder" -> False
                            |>]
                    ]
                ],

            (* Time-dependent finite-duration event or TD BackgroundGenerator *)
            Length[genList] >= 2 && (containsTDGeneratorQ[genList[[1]]] || containsTDGeneratorQ[bgActive]),
                With[{tdSt = compileTDStage[genList[[1]], bgActive, taExpr, durExpr, phSyms, xTest, space, basis, opbasis, nSub]},
                    If[tdSt === $Failed, Return[$Failed, Module]];
                    AppendTo[stages, tdSt]
                ],

            (* Time-independent finite-duration event {evgen, dur} *)
            Length[genList] >= 2,
                With[{evG = If[genList[[1]] === None, 0, genList[[1]]],
                      bgG = If[bgActive === None, 0, bgActive]},
                    expExpr = If[space === "Hilbert",
                        durExpr * (evG + bgG),
                        durExpr * (
                            If[evG === 0 || SuperoperatorQ[evG], evG, CommutationSuperoperator[-I * evG]] +
                            If[bgG === 0 || SuperoperatorQ[bgG], bgG, CommutationSuperoperator[-I * bgG]]
                        )
                    ];
                    gSpec = decomposeGenerator[expExpr, phSyms, xTest, space, basis, opbasis, space === "Hilbert"];
                    If[gSpec === $Failed, Return[$Failed, Module]];
                    If[!TrueQ[gSpec["IsZero"]],
                        AppendTo[stages,
                            If[gSpec["IsConstant"],
                                <|"Kind" -> If[space === "Hilbert", "Static", "LiouvilleStatic"],
                                  "DepMask" -> gSpec["DepMask"],
                                  "DepCols" -> {},
                                  "Propagator" -> toPackedComplexMat[MatrixExp[If[space === "Hilbert", -I, 1.] * gSpec["Constant"]]],
                                  "IsDifferentiable" -> True,
                                  "SupportsSecondOrder" -> True|>,
                                <|"Kind" -> If[space === "Hilbert", "HilbertExpm", "LiouvilleExpm"],
                                  "DepMask" -> gSpec["DepMask"],
                                  "DepCols" -> gSpec["DepCols"],
                                  "Matrices" -> gSpec["Matrices"],
                                  "EvalCoeffs" -> gSpec["EvalCoeffs"],
                                  "EvalCoeffJac" -> gSpec["EvalCoeffJac"],
                                  "EvalCoeffHVP" -> gSpec["EvalCoeffHVP"],
                                  "IsScalarTimeGenerator" -> (!FreeQ[durExpr, Alternatives @@ phSyms] && FreeQ[evG + bgG, Alternatives @@ phSyms] && gSpec["NumAtoms"] === 1),
                                  "IsDifferentiable" -> gSpec["IsDifferentiable"],
                                  "SupportsSecondOrder" -> (space === "Hilbert" && gSpec["IsDifferentiable"])|>
                            ]
                        ]
                    ]
                ],
            True, Return[$Failed, Module]
        ];
        tAccExpr = If[durExpr === Undefined || durExpr === None, tAccExpr, tAccExpr + durExpr],
        {ev, eventsChrono}
    ];

    If[stages === {},
        stages = {<|"Kind" -> If[space === "Hilbert", "Static", "LiouvilleStatic"],
                    "DepMask" -> ConstantArray[False, nP], "DepCols" -> {},
                    "Propagator" -> toPackedComplexMat[IdentityMatrix[dim]],
                    "IsDifferentiable" -> True, "SupportsSecondOrder" -> True|>}
    ];

    (* Fuse consecutive static propagators *)
    If[TrueQ[fuseQ],
        Do[
            With[{s = stages[[j]]},
                If[Length[fusedStages] > 0 && MemberQ[{"Static", "LiouvilleStatic"}, s["Kind"]] &&
                   Last[fusedStages]["Kind"] === s["Kind"],
                    fusedStages[[-1, "Propagator"]] = toPackedComplexMat[s["Propagator"] . Last[fusedStages]["Propagator"]],
                    AppendTo[fusedStages, s]
                ]
            ],
            {j, Length[stages]}
        ],
        fusedStages = stages
    ];

    (* Assign deterministic integer IDs *)
    fusedStages = Table[Append[fusedStages[[j]], "ID" -> j], {j, Length[fusedStages]}];

    (* Decompose initial density operator, observable, and normalization factor *)
    If[space === "Hilbert",
        rhoSpec = decomposeGenerator[rhoExpr, phSyms, xTest, "Hilbert", basis, opbasis, False];
        obsSpec = decomposeGenerator[obsExpr, phSyms, xTest, "Hilbert", basis, opbasis, False],
        rhoSpec = decomposeLiouvilleState[rhoExpr, phSyms, xTest, opbasis, basis, False];
        obsSpec = decomposeLiouvilleState[obsExpr, phSyms, xTest, opbasis, basis, True]
    ];
    If[rhoSpec === $Failed || obsSpec === $Failed, Return[$Failed]];

    normSpec = If[normExpr === None || normExpr === False,
        <|"IsConstant" -> True, "Constant" -> 1., "DepMask" -> ConstantArray[False, nP], "IsDifferentiable" -> True|>,
        If[FreeQ[normExpr, Alternatives @@ phSyms],
            If[!NumericQ[normExpr] || N[Re[normExpr]] == 0., Return[$Failed]];
            <|"IsConstant" -> True, "Constant" -> N[Re[normExpr]], "DepMask" -> ConstantArray[False, nP], "IsDifferentiable" -> True|>,
            Module[{nDiff = differentiableExprQ[normExpr, phSyms], cNorm, gNormExpr, hNormExpr, xSym = Unique["gnx"]},
                cNorm = Refine[normExpr, Assumptions -> (Element[#, Reals] & /@ phSyms)];
                If[TrueQ[nDiff],
                    gNormExpr = Table[Refine[D[cNorm, ph], Assumptions -> (Element[#, Reals] & /@ phSyms)], {ph, phSyms}];
                    hNormExpr = Table[Refine[D[gNormExpr[[a]], phSyms[[b]]], Assumptions -> (Element[#, Reals] & /@ phSyms)], {a, nP}, {b, nP}];
                    If[!FreeQ[{gNormExpr, hNormExpr}, _Derivative | _D], nDiff = False]
                ];
                <|
                    "IsConstant" -> False,
                    "DepMask" -> Table[!FreeQ[cNorm, phSyms[[a]]], {a, nP}],
                    "Eval" -> buildPureScalarFn[cNorm, phSyms, nP],
                    "EvalGrad" -> If[nDiff, buildPureVecFn[gNormExpr, phSyms, nP], None],
                    "EvalHess" -> If[nDiff, buildPureMatFn[hNormExpr, phSyms, nP], None],
                    "IsDifferentiable" -> nDiff
                |>
            ]
        ]
    ];
    If[normSpec === $Failed, Return[$Failed]];

    isDiff = AllTrue[fusedStages, TrueQ[#["IsDifferentiable"]] &] &&
        TrueQ[rhoSpec["IsDifferentiable"]] && TrueQ[obsSpec["IsDifferentiable"]] && TrueQ[normSpec["IsDifferentiable"]];
    supp2nd = (space === "Hilbert") && isDiff && AllTrue[fusedStages, TrueQ[#["SupportsSecondOrder"]] &];

    <|
        "Space" -> space,
        "Dim" -> dim,
        "NumParams" -> nP,
        "Stages" -> fusedStages,
        "NumRawStages" -> Length[stages],
        "NumFusedStages" -> Length[fusedStages],
        "RhoSpec" -> rhoSpec,
        "ObsSpec" -> obsSpec,
        "NormSpec" -> normSpec,
        "IsDifferentiable" -> isDiff,
        "SupportsSecondOrder" -> supp2nd
    |>
];

(* Global DependencyMask summary across compiled stages and boundary terms *)
computeProblemDepMask[members_List, nP_Integer] := Module[{mask = ConstantArray[False, nP]},
    Do[
        Do[mask = MapThread[Or, {mask, s["DepMask"]}], {s, m["Stages"]}];
        If[KeyExistsQ[m, "RhoSpec"], mask = MapThread[Or, {mask, m["RhoSpec"]["DepMask"]}]];
        If[KeyExistsQ[m, "ObsSpec"], mask = MapThread[Or, {mask, m["ObsSpec"]["DepMask"]}]];
        If[KeyExistsQ[m, "NormSpec"], mask = MapThread[Or, {mask, m["NormSpec"]["DepMask"]}]],
        {m, members}
    ];
    mask
];

(* Safe interior test point strictly inside bounds for compile-time verification *)
safeTestPoint[x0_List, lb_List, ub_List] := Table[
    With[{x = x0[[i]], l = lb[[i]], u = ub[[i]]},
        Which[
            NumericQ[l] && Abs[l] < Infinity && NumericQ[u] && Abs[u] < Infinity,
                If[x > l + 0.05 (u - l) && x < u - 0.05 (u - l), x, 0.5 * (l + u)],
            NumericQ[l] && Abs[l] < Infinity,
                If[x > l + 1.*^-4, x, l + Max[0.01, 0.1 Abs[l]]],
            NumericQ[u] && Abs[u] < Infinity,
                If[x < u - 1.*^-4, x, u - Max[0.01, 0.1 Abs[u]]],
            True,
                If[Abs[x] > 1.*^-4, x, 0.17]
        ]
    ],
    {i, Length[x0]}
];

CompileTransformationOptimization[prob_TransformationOptimizationProblem, opts:(_Rule|_RuleDelayed)...] := Module[
    {tComp0 = AbsoluteTime[], data = First[prob], optsAssoc, phSyms, x0, lb, ub, xTest, nP,
     rhoExpr, obsExpr, evExpr, probOptsEval, allOpts, bgOpt, ensOpt, chronOpt, normOpt, tIniOpt, basisOpt, opBasisOpt,
     spaceOpt, nSub, fuseQ, enableCacheQ, allowBB, parOpt, ensStd, ensVars, ensPairs, weights, members, evalBBFn, compMode = "Compiled"},
    optsAssoc = Join[Association[Options[CompileTransformationOptimization]], Association[FilterRules[{opts}, Except[_RuleDelayed]]]];
    nP = Length[data["X0Vec"]];
    x0 = data["X0Vec"];
    lb = data["LBVec"];
    ub = data["UBVec"];
    xTest = safeTestPoint[x0, lb, ub];
    phSyms = Table[Unique["goptP"], {nP}];

    (* Isolate problem symbols inside Block before evaluating SpinDynamica constructs *)
    {rhoExpr, obsExpr, evExpr, probOptsEval} = evalIsolated[
        data["VarHold"], phSyms,
        Replace[{data["RhoHold"], data["ObsHold"], data["EventsHold"], data["OptionsHold"]},
            {Hold[r_], Hold[o_], Hold[e_], Hold[op___]} :> Hold[{r, o, e, FilterRules[{op}, Except[_RuleDelayed]]}]
        ]
    ];

    allOpts = Join[FilterRules[{opts}, Except[_RuleDelayed]], probOptsEval];
    spaceOpt = SimulationSpace /. allOpts /. {SimulationSpace -> Automatic};
    nSub = Max[2, Round[Substeps /. allOpts /. {Substeps -> 16}]];
    fuseQ = TrueQ[FuseStaticEvents /. allOpts /. {FuseStaticEvents -> True}];
    enableCacheQ = TrueQ[EnableCaching /. allOpts /. {EnableCaching -> True}];
    allowBB = TrueQ[AllowBlackBoxFallback /. allOpts /. {AllowBlackBoxFallback -> True}];
    parOpt = Parallel /. allOpts /. {Parallel -> False};
    chronOpt = Chronology /. allOpts /. {Chronology -> $Chronology};
    normOpt = NormalizationFactor /. allOpts /. {NormalizationFactor -> None};
    tIniOpt = InitialTimePoint /. allOpts /. {InitialTimePoint -> 0};
    basisOpt = Basis /. allOpts /. {Basis -> data["Basis"]};
    opBasisOpt = OperatorBasis /. allOpts /. {OperatorBasis -> data["OperatorBasis"]};
    bgOpt = BackgroundGenerator /. allOpts /. {BackgroundGenerator -> None};
    ensOpt = EnsembleAverage /. allOpts /. {EnsembleAverage -> False};

    evalBBFn = Function[{xIn},
        Module[{taArgs},
            taArgs = evalAtNumericValues[
                data["VarHold"], N[Re[xIn]],
                Replace[{data["RhoHold"], data["ObsHold"], data["EventsHold"], data["OptionsHold"]},
                    {Hold[r_], Hold[o_], Hold[e_], Hold[op___]} :>
                        Hold[Quiet[TransformationAmplitude[r -> o, e,
                            Evaluate[Sequence @@ FilterRules[{op}, Options[TransformationAmplitude]]]]]]
                ]
            ];
            If[NumericQ[taArgs], N[Re[taArgs]], $Failed]
        ]
    ];

    If[evExpr === $Failed || rhoExpr === $Failed || obsExpr === $Failed ||
       !FreeQ[Hold[evExpr, rhoExpr, obsExpr], $Failed],
        If[!allowBB,
            Message[CompileTransformationOptimization::compfail, "Input problem contained $Failed."];
            Return[$Failed]
        ];
        compMode = "BlackBoxFallback";
        members = {<|"Space" -> "BlackBox", "NumParams" -> nP, "Stages" -> {}, "IsDifferentiable" -> False,
                     "SupportsSecondOrder" -> False, "EvalBlackBox" -> evalBBFn|>};
        weights = {1.},
        (* Compile regular or ensemble members *)
        If[ensOpt === False || ensOpt === None,
            With[{m = Quiet[compileSingleMember[rhoExpr, obsExpr, evExpr, bgOpt, chronOpt, tIniOpt, normOpt,
                    phSyms, xTest, spaceOpt, basisOpt, opBasisOpt, nSub, fuseQ]]},
                If[m === $Failed,
                    If[!allowBB,
                        Message[CompileTransformationOptimization::compfail, "Operator/event decomposition failed."];
                        Return[$Failed]
                    ];
                    compMode = "BlackBoxFallback";
                    members = {<|"Space" -> "BlackBox", "NumParams" -> nP, "Stages" -> {}, "IsDifferentiable" -> False,
                                 "SupportsSecondOrder" -> False, "EvalBlackBox" -> evalBBFn|>} ,
                    members = {m}
                ]
            ];
            weights = {1.},
            (* Pre-expand EnsembleAverage grid once at compile time *)
            ensStd = Quiet[Ensemble`EnsembleAverage`StandardizeEnsembleAverage[ensOpt]];
            If[!MatchQ[ensStd, {_List, {{_List, Except[_List]} ..}}],
                compMode = "BlackBoxFallback";
                members = {<|"Space" -> "BlackBox", "NumParams" -> nP, "Stages" -> {}, "IsDifferentiable" -> False,
                             "SupportsSecondOrder" -> False, "EvalBlackBox" -> evalBBFn|>};
                weights = {1.},
                {ensVars, ensPairs} = ensStd;
                weights = Developer`ToPackedArray[N[Re[ensPairs[[All, 2]]]]];
                members = Table[
                    With[{rep = Thread[ensVars -> ensPairs[[k, 1]]]},
                        compileSingleMember[
                            rhoExpr /. rep, obsExpr /. rep, evExpr /. rep,
                            If[bgOpt === None, None, bgOpt /. rep],
                            chronOpt, tIniOpt, If[normOpt === None, None, normOpt /. rep],
                            phSyms, xTest, spaceOpt, basisOpt, opBasisOpt, nSub, fuseQ
                        ]
                    ],
                    {k, Length[ensPairs]}
                ];
                If[MemberQ[members, $Failed],
                    If[!allowBB,
                        Message[CompileTransformationOptimization::compfail, "Ensemble member compilation failed."];
                        Return[$Failed]
                    ];
                    compMode = "BlackBoxFallback";
                    members = {<|"Space" -> "BlackBox", "NumParams" -> nP, "Stages" -> {}, "IsDifferentiable" -> False,
                                 "SupportsSecondOrder" -> False, "EvalBlackBox" -> evalBBFn|>};
                    weights = {1.}
                ]
            ]
        ]
    ];

    <|
        "Problem" -> prob,
        "VarNames" -> data["VarNames"],
        "VarHold" -> data["VarHold"],
        "NumParams" -> nP,
        "X0Vec" -> x0,
        "Bounds" -> {lb, ub},
        "ScaleVec" -> Optimization`Solvers`Private`resolveScaleVector[x0, lb, ub, Automatic],
        "CompilationMode" -> compMode,
        "Space" -> If[compMode === "BlackBoxFallback", "BlackBox", members[[1, "Space"]]],
        "Members" -> members,
        "Weights" -> weights,
        "NumMembers" -> Length[members],
        "NumStagesTotal" -> Total[(Length[#["Stages"]] &) /@ members],
        "NumRawStagesTotal" -> Total[(Lookup[#, "NumRawStages", 0] &) /@ members],
        "IsDifferentiable" -> (compMode === "Compiled" && AllTrue[members, TrueQ[#["IsDifferentiable"]] &]),
        "SupportsSecondOrder" -> (compMode === "Compiled" && AllTrue[members, TrueQ[#["SupportsSecondOrder"]] &]),
        "GlobalDepMask" -> If[compMode === "Compiled", computeProblemDepMask[members, nP], ConstantArray[True, nP]],
        "EnableCaching" -> enableCacheQ,
        "Parallel" -> parOpt,
        "CompileTime" -> AbsoluteTime[] - tComp0
    |>
];

(* =============================================================== *)
(* Evaluation, Verification, Optimization & Result Object          *)
(* =============================================================== *)
extractNumericPoint[comp_Association, xOrRules_] := Module[{names = comp["VarNames"], xVec},
    Which[
        VectorQ[xOrRules, NumericQ] && Length[xOrRules] === comp["NumParams"],
            Developer`ToPackedArray[N[Re[xOrRules]]],
        MatchQ[xOrRules, {_Rule ..}],
            xVec = Table[
                With[{symRule = SelectFirst[xOrRules, ToString[First[#]] === names[[i]] &, Missing[]]},
                    If[MissingQ[symRule], comp["X0Vec"][[i]], N[Re[Last[symRule]]]]
                ],
                {i, comp["NumParams"]}
            ];
            If[VectorQ[xVec, NumericQ], Developer`ToPackedArray[xVec], $Failed],
        True, $Failed
    ]
];

EvaluateTransformationOptimization[comp_Association, opts:(_Rule|_RuleDelayed)...] :=
    EvaluateTransformationOptimization[comp, Automatic, opts];

EvaluateTransformationOptimization[comp_Association, xOrRules:(_List | Automatic): Automatic, opts:(_Rule|_RuleDelayed)...] := Module[
    {xVec, derivOpt, fdOpts, oracles},
    xVec = If[xOrRules === Automatic, comp["X0Vec"], extractNumericPoint[comp, xOrRules]];
    If[xVec === $Failed, Message[OptimizeTransformationAmplitude::invalp, xOrRules]; Return[$Failed]];
    derivOpt = DerivativeMethod /. {opts} /. {DerivativeMethod -> Automatic};
    fdOpts = {
        "FiniteDifferenceOrder" -> (FiniteDifferenceOrder /. {opts} /. {FiniteDifferenceOrder -> "Central"}),
        "FiniteDifferenceStep" -> (FiniteDifferenceStep /. {opts} /. {FiniteDifferenceStep -> Automatic})
    };
    oracles = EvaluateProblemOracles[comp, derivOpt, fdOpts];
    If[derivOpt === "None" || oracles["EvalGrad"] === None,
        oracles["EvalObj"][xVec],
        oracles["EvalGrad"][xVec]
    ]
];

EvaluateTransformationOptimization[prob_TransformationOptimizationProblem, opts:(_Rule|_RuleDelayed)...] :=
    EvaluateTransformationOptimization[prob, Automatic, opts];

EvaluateTransformationOptimization[prob_TransformationOptimizationProblem, xOrRules:(_List | Automatic): Automatic, opts:(_Rule|_RuleDelayed)...] :=
    With[{comp = CompileTransformationOptimization[prob, Sequence @@ FilterRules[{opts}, Options[CompileTransformationOptimization]]]},
        If[comp === $Failed, $Failed, EvaluateTransformationOptimization[comp, xOrRules, opts]]
    ];

VerifyTransformationOptimizationGradient[comp_Association, opts:(_Rule|_RuleDelayed)...] :=
    VerifyTransformationOptimizationGradient[comp, Automatic, opts];

VerifyTransformationOptimizationGradient[comp_Association, xOrRules:(_List | Automatic): Automatic, opts:(_Rule|_RuleDelayed)...] := Module[
    {xVec, eps, oracles, anaRes, gAna, gNum, fdOracles, fdRes, absErr, relErr, tol},
    xVec = If[xOrRules === Automatic, comp["X0Vec"], extractNumericPoint[comp, xOrRules]];
    If[xVec === $Failed, Return[$Failed]];
    eps = FiniteDifferenceStep /. {opts} /. {FiniteDifferenceStep -> 1.*^-6};
    tol = Tolerance /. {opts} /. {Tolerance -> 1.*^-5};
    oracles = EvaluateProblemOracles[comp, "Analytical", {}];
    If[oracles["EvalGrad"] === None, Return[<|"Passed" -> False, "Reason" -> "Analytical gradient unavailable"|>]];
    anaRes = oracles["EvalGrad"][xVec];
    gAna = anaRes["Gradient"];
    fdOracles = BuildFiniteDiffOracles[comp, oracles["EvalObj"], {"FiniteDifferenceOrder" -> "Central", "FiniteDifferenceStep" -> eps}];
    fdRes = fdOracles[xVec];
    gNum = fdRes["Gradient"];
    absErr = Max[Abs[gAna - gNum]];
    relErr = absErr / Max[1., Norm[gNum, Infinity]];
    <|
        "Objective" -> anaRes["Objective"],
        "AnalyticalGradient" -> gAna,
        "NumericalGradient" -> gNum,
        "MaxAbsError" -> absErr,
        "RelativeError" -> relErr,
        "Passed" -> TrueQ[relErr <= tol]
    |>
];

VerifyTransformationOptimizationGradient[prob_TransformationOptimizationProblem, opts:(_Rule|_RuleDelayed)...] :=
    VerifyTransformationOptimizationGradient[prob, Automatic, opts];

VerifyTransformationOptimizationGradient[prob_TransformationOptimizationProblem, xOrRules:(_List | Automatic): Automatic, opts:(_Rule|_RuleDelayed)...] :=
    With[{comp = CompileTransformationOptimization[prob, Sequence @@ FilterRules[{opts}, Options[CompileTransformationOptimization]]]},
        If[comp === $Failed, $Failed, VerifyTransformationOptimizationGradient[comp, xOrRules, opts]]
    ];

VerifyTransformationOptimizationHVP[comp_Association, opts:(_Rule|_RuleDelayed)...] :=
    VerifyTransformationOptimizationHVP[comp, Automatic, Automatic, opts];

VerifyTransformationOptimizationHVP[comp_Association, xOrRules:(_List | Automatic): Automatic, dirIn:(_List | Automatic): Automatic, opts:(_Rule|_RuleDelayed)...] := Module[
    {xVec, dirVec, nP = comp["NumParams"], eps, tol, oracles, cRes, st, hvpAna, gP, gM, hvpFD, hMat, hAsym, absErr, relErr},
    xVec = If[xOrRules === Automatic, comp["X0Vec"], extractNumericPoint[comp, xOrRules]];
    If[xVec === $Failed, Return[$Failed]];
    dirVec = If[dirIn === Automatic,
        Developer`ToPackedArray[Table[If[Mod[k, 2] === 1, 1., -0.7] / comp["ScaleVec"][[k]], {k, nP}]],
        Developer`ToPackedArray[N[Re[Flatten[dirIn]]]]
    ];
    eps = FiniteDifferenceStep /. {opts} /. {FiniteDifferenceStep -> 1.*^-6};
    tol = Tolerance /. {opts} /. {Tolerance -> 1.*^-5};
    oracles = EvaluateProblemOracles[comp, "SecondOrder", {}];
    If[oracles["EvalWithCache"] === None || oracles["EvalHVP"] === None,
        Return[<|"Passed" -> False, "Reason" -> "Second-order HVP unavailable"|>]
    ];
    cRes = oracles["EvalWithCache"][xVec];
    st = cRes["EvalState"];
    hvpAna = oracles["EvalHVP"][st, dirVec];
    With[{hStep = eps * Max[1., Norm[dirVec * comp["ScaleVec"]]] / Max[Norm[dirVec], 1.*^-12]},
        gP = oracles["EvalGrad"][xVec + hStep * dirVec]["Gradient"];
        gM = oracles["EvalGrad"][xVec - hStep * dirVec]["Gradient"];
        hvpFD = Developer`ToPackedArray[(gP - gM) / (2. * hStep)]
    ];
    absErr = Max[Abs[hvpAna - hvpFD]];
    relErr = absErr / Max[1., Norm[hvpFD, Infinity]];
    hMat = oracles["EvalHessian"][st];
    hAsym = Max[Abs[hMat - Transpose[hMat]]];
    <|
        "AnalyticalHVP" -> hvpAna,
        "FiniteDiffHVP" -> hvpFD,
        "MaxAbsError" -> absErr,
        "RelativeError" -> relErr,
        "HessianMaxAsymmetry" -> hAsym,
        "HessianMat" -> hMat,
        "Passed" -> TrueQ[relErr <= tol]
    |>
];

VerifyTransformationOptimizationHVP[prob_TransformationOptimizationProblem, opts:(_Rule|_RuleDelayed)...] :=
    VerifyTransformationOptimizationHVP[prob, Automatic, Automatic, opts];

VerifyTransformationOptimizationHVP[prob_TransformationOptimizationProblem, xOrRules:(_List | Automatic): Automatic, dirIn:(_List | Automatic): Automatic, opts:(_Rule|_RuleDelayed)...] :=
    With[{comp = CompileTransformationOptimization[prob, Sequence @@ FilterRules[{opts}, Options[CompileTransformationOptimization]]]},
        If[comp === $Failed, $Failed, VerifyTransformationOptimizationHVP[comp, xOrRules, dirIn, opts]]
    ];

(* Result wrapper Accessors & Formatting *)
TransformationOptimizationResult[assoc_Association][prop_String] := Switch[prop,
    "Value" | "Objective", assoc["Objective"],
    "Rules" | "Parameters", assoc["Rules"],
    "ParametersVec", assoc["ParametersVec"],
    "LegacyFormat", {assoc["Objective"], assoc["Rules"]},
    "Association", assoc,
    "Keys", Keys[assoc],
    _, Lookup[assoc, prop, Missing["KeyAbsent", prop]]
];

TransformationOptimizationResult[assoc_Association][prop1_String, props__String] :=
    Fold[If[AssociationQ[#1], Lookup[#1, #2, Missing["KeyAbsent", #2]], Missing["KeyAbsent", #2]] &, assoc, {prop1, props}];

Format[TransformationOptimizationResult[assoc_Association]] := StringForm[
    "TransformationOptimizationResult[Value -> ``, Optimizer -> ``, Deriv -> ``, Iters -> ``, Time -> ``s]",
    NumberForm[N[assoc["Objective"]], 6],
    assoc["Optimizer"],
    assoc["DerivativeMethod"],
    assoc["Iterations"],
    NumberForm[N[assoc["ElapsedTime"]], 4]
];

(* Main Optimization Driver accepting TransformationOptimizationProblem *)
optimizeProblemDispatch[prob_TransformationOptimizationProblem, opts:(_Rule|_RuleDelayed)...] := Module[
    {optsNoDelay, comp, derivOpt, optChoice, oracles, fdOpts, solverOpts, solRes, xSol, rules,
     retObjQ, cacheHits, cacheMisses, cacheStats, resAssoc},
    optsNoDelay = FilterRules[{opts}, Except[_RuleDelayed]];
    comp = CompileTransformationOptimization[prob, Sequence @@ FilterRules[optsNoDelay, Options[CompileTransformationOptimization]]];
    If[comp === $Failed, Return[$Failed]];

    derivOpt = DerivativeMethod /. optsNoDelay /. Options[OptimizeTransformationAmplitude];
    optChoice = Optimizer /. optsNoDelay /. Options[OptimizeTransformationAmplitude];

    (* Intelligent Strategy & Solver coordination *)
    If[derivOpt === Automatic,
        derivOpt = Which[
            MemberQ[{"NelderMead", "Powell"}, optChoice] || (ListQ[optChoice] && MemberQ[{"NelderMead", "Powell"}, First[optChoice]]), "None",
            !TrueQ[comp["IsDifferentiable"]] && comp["CompilationMode"] === "Compiled", "None",
            !TrueQ[comp["IsDifferentiable"]], "FiniteDifference",
            True, "Analytical"
        ]
    ];
    If[optChoice === Automatic,
        optChoice = Which[
            derivOpt === "None", "NelderMead",
            derivOpt === "SecondOrder" && TrueQ[comp["SupportsSecondOrder"]], "HybridLBFGSNewtonCG",
            True, "LBFGS"
        ]
    ];

    fdOpts = {
        "FiniteDifferenceOrder" -> (FiniteDifferenceOrder /. optsNoDelay /. Options[OptimizeTransformationAmplitude]),
        "FiniteDifferenceStep" -> (FiniteDifferenceStep /. optsNoDelay /. Options[OptimizeTransformationAmplitude])
    };
    oracles = EvaluateProblemOracles[comp, derivOpt, fdOpts];

    solverOpts = FilterRules[
        Join[FilterRules[{opts}, Except[Optimizer]], {Optimizer -> optChoice}],
        Options[RunGeneralOptimizer]
    ];
    solRes = RunGeneralOptimizer[oracles, comp["X0Vec"], comp["Bounds"], Sequence @@ solverOpts];
    If[solRes === $Failed || !AssociationQ[solRes],
        Message[OptimizeTransformationAmplitude::badsol, "Solver returned $Failed."];
        Return[$Failed]
    ];

    xSol = solRes["ParametersVec"];
    rules = MapThread[
        Function[{vh, val}, Replace[vh, Hold[s_] :> Rule[Unevaluated[s], val]]],
        {First[prob]["VarHoldsList"], xSol}
    ];

    cacheHits = Total[If[Head[#] === DataStructure && #["KeyExistsQ", "__Hits"], #["Lookup", "__Hits"], 0] & /@ oracles["Caches"]];
    cacheMisses = Total[If[Head[#] === DataStructure && #["KeyExistsQ", "__Misses"], #["Lookup", "__Misses"], 0] & /@ oracles["Caches"]];
    cacheStats = <|
        "StagesCompiled" -> comp["NumStagesTotal"],
        "StagesRaw" -> comp["NumRawStagesTotal"],
        "StaticFusionReduction" -> Max[0, comp["NumRawStagesTotal"] - comp["NumStagesTotal"]],
        "CacheHits" -> cacheHits,
        "CacheMisses" -> cacheMisses
    |>;

    resAssoc = Join[solRes, <|
        "Rules" -> rules,
        "VarNames" -> comp["VarNames"],
        "Space" -> comp["Space"],
        "CompilationMode" -> comp["CompilationMode"],
        "DerivativeMethod" -> oracles["EffectiveDerivativeMethod"],
        "CompileTime" -> comp["CompileTime"],
        "ElapsedTime" -> comp["CompileTime"] + solRes["ElapsedTime"],
        "NumMembers" -> comp["NumMembers"],
        "CacheStats" -> cacheStats
    |>];

    retObjQ = TrueQ[ReturnResultObject /. optsNoDelay /. Options[OptimizeTransformationAmplitude]];
    If[retObjQ,
        TransformationOptimizationResult[resAssoc],
        {resAssoc["Objective"], rules}
    ]
];

(* HoldAll dispatchers for pre-constructed problem variable OR direct 3-argument specification *)
OptimizeTransformationAmplitude[pairSym_Symbol, events:Except[_Rule|_RuleDelayed], varspecs:Except[_Rule|_RuleDelayed], opts:(_Rule|_RuleDelayed)...] /;
    MatchQ[OwnValues[pairSym], {(_ :> (_ -> _)) ..}] :=
    Replace[First[OwnValues[pairSym]],
        (_ :> (r_ -> o_)) :> OptimizeTransformationAmplitude[r -> o, events, varspecs, opts]
    ];

OptimizeTransformationAmplitude[HoldPattern[rhoIni_ -> obs_], events:Except[_Rule|_RuleDelayed], varspecs:Except[_Rule|_RuleDelayed], opts:(_Rule|_RuleDelayed)...] := Module[
    {probOptPatt, heldAll, probOptsHeld, optOptsHeld, prob},
    probOptPatt = Alternatives[ConstrainPositive, BackgroundGenerator, EnsembleAverage, Chronology, NormalizationFactor, Basis, OperatorBasis, InitialTimePoint, FinalTimePoint];
    heldAll = HoldComplete[opts];
    probOptsHeld = Cases[heldAll, r:((Rule|RuleDelayed)[p_ /; MatchQ[Unevaluated[p], probOptPatt], _]) :> HoldComplete[r], {1}];
    optOptsHeld = Cases[heldAll, r:((Rule|RuleDelayed)[p_ /; !MatchQ[Unevaluated[p], probOptPatt], _]) :> HoldComplete[r], {1}];
    prob = Replace[
        Fold[Join, HoldComplete[], probOptsHeld],
        HoldComplete[pos___] :> TransformationOptimizationProblem[rhoIni -> obs, events, varspecs, pos]
    ];
    If[prob === $Failed, Return[$Failed]];
    Replace[
        Fold[Join, HoldComplete[], optOptsHeld],
        HoldComplete[oos___] :> optimizeProblemDispatch[prob, oos]
    ]
];

OptimizeTransformationAmplitude[probExpr:Except[_Rule|_RuleDelayed], opts:(_Rule|_RuleDelayed)...] := With[{p = probExpr},
    If[MatchQ[p, _TransformationOptimizationProblem],
        optimizeProblemDispatch[p, opts],
        Message[OptimizeTransformationAmplitude::badsol, "Expected TransformationOptimizationProblem or rhoIni -> obs specification."];
        $Failed
    ]
];

End[];
EndPackage[];
