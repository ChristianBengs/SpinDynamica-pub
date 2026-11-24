(* SpinDynamica 3.10.2 -- products of spin states and operators *)
BeginPackage["QM`SpinDot`"];
With[{contexts = {"QM`Hilbert`", "QM`SpinOperators`", "QM`Liouville`", "General`Matrices`",
        "General`UserLevel`", "SpinDynamica`"}},
    Scan[Needs, contexts];
    $ContextPath = DeleteDuplicates[Join[contexts, $ContextPath]]
];
Begin["`Private`"];

(* Install upvalues on the head of an actual operand. Plus and Times can be
   the heads of superpositions and scaled states. Lists and sparse arrays
   need no rules. Preserve the protection state of every symbol. *)
withHeads[heads_List, install_Function] := Scan[
    Function[tag, With[{protected = MemberQ[Attributes[tag], Protected]},
        Internal`WithLocalSettings[
            If[protected, Unprotect[tag]], install[tag], If[protected, Protect[tag]]]
    ]], heads];

dotFactors[only_] := only;
dotFactors[first_, second_, more___] := Dot[first, second, more];

withHeads[{Ket, Plus, Times}, Function[tag,
    tag /: HoldPattern[Dot[a___, n_?NumericQ, k_tag, b___]] /; KetQ[k] := n dotFactors[a, k, b];
    tag /: HoldPattern[Dot[k1_tag, k2_?KetQ]] /; UserLevel[] === 1 && KetQ[k1] :=
        (Message[Ket::illegketket]; $Failed);
    tag /: HoldPattern[Dot[b_?BraQ, k_tag]] /; KetQ[k] := With[{basis = ZeemanBasis[SpinSystem[Expand[k]]]},
        Dot[VectorRepresentation[b, basis], VectorRepresentation[k, basis]]];
    tag /: HoldPattern[Dot[left___, b_?BraQ, k_tag, right___]] /;
            KetQ[k] && Length[{left, right}] > 0 := Dot[b, k] dotFactors[left, right];
    tag /: HoldPattern[Dot[k1_?KetQ, b_?BraQ][k2_tag]] /; KetQ[k2] := k1 Dot[b, k2];
    tag /: HoldPattern[Dot[op_?OperatorQ, k_tag]] /; KetQ[k] := With[{basis = ZeemanBasis[]},
        Inner[Times, MatrixRepresentation[op, basis] . VectorRepresentation[k, basis], BasisKets[basis], Plus]];
    tag /: HoldPattern[Dot[first_?OperatorQ, rest__?OperatorQ, k_tag]] /; KetQ[k] :=
        Operator[OperatorProduct[first, rest]] . k;
    tag /: HoldPattern[Dot[ops:(_?OperatorQ | _?NumericQ).., k_tag]] /; Length[{ops}] >= 2 && KetQ[k] :=
        OperatorProduct[ops, Sorted -> False] . k;
]];

withHeads[{Bra, Plus, Times}, Function[tag,
    tag /: HoldPattern[Dot[a___, n_?NumericQ, b_tag, z___]] /; BraQ[b] := n dotFactors[a, b, z];
    tag /: HoldPattern[Dot[b1_tag, b2_?BraQ]] /; UserLevel[] === 1 && BraQ[b1] :=
        (Message[Bra::illegbrabra]; $Failed);
    tag /: HoldPattern[Dot[b_tag, op_?OperatorQ]] /; BraQ[b] := Adjoint[Adjoint[op] . Adjoint[b]];
    tag /: HoldPattern[Dot[b_tag, ops:(_?OperatorQ | _?NumericQ)..]] /; Length[{ops}] >= 2 && BraQ[b] :=
        b . OperatorProduct[ops, Sorted -> False];
    tag /: HoldPattern[Dot[b_tag, ops:(_?OperatorQ | _?NumericQ).., k_?KetQ]] /; Length[{ops}] >= 2 && BraQ[b] :=
        Dot[b, OperatorProduct[ops, Sorted -> False], k];
]];

(* Times is Orderless, so one pattern covers a coefficient on either side. *)
withHeads[{Times}, Function[tag,
    tag /: HoldPattern[Dot[Times[coef__, k_?KetQ], b_?BraQ]] := Times[coef, Dot[k, b]];
    tag /: HoldPattern[Dot[k_?KetQ, Times[coef__, b_?BraQ]]] := Times[coef, Dot[k, b]];
    tag /: HoldPattern[Dot[Times[coef__, op_?OperatorQ], k_?KetQ]] := Times[coef, Dot[op, k]];
]];
withHeads[{Plus}, Function[tag,
    tag /: HoldPattern[Dot[Plus[ks__?KetQ], b_?BraQ]] := Total[Dot[#, b] & /@ {ks}];
    tag /: HoldPattern[Dot[k_?KetQ, Plus[bs__?BraQ]]] := Total[Dot[k, #] & /@ {bs}];
]];

(* Anchor the product on an operator. Numeric coefficients may flank it. *)
withHeads[{Operator, opI, opT, opR, ProjectionOperator, ShiftOperator,
    SpinExchangeOperator, SpinPermutationOperator, UnityOperator, NullOperator,
    ScalarOperatorProduct, OperatorProduct, Plus, Times, Power}, Function[tag,
    tag /: HoldPattern[Dot[before:(_?OperatorQ | _?NumericQ)..., op_tag,
        after:(_?OperatorQ | _?NumericQ)...]] /; Length[{before, after}] > 0 && OperatorQ[op] :=
        OperatorProduct[before, op, after, Sorted -> False];
]];

(* An unevaluated interpolation at symbolic t has a compound head. Store
   this rule explicitly on InterpolatingFunction, one level deeper than /: accepts. *)
withHeads[{InterpolatingFunction}, Function[tag,
    With[{rule = HoldPattern[Dot[before:(_?OperatorQ | _?NumericQ)...,
            value:(_InterpolatingFunction[_]), after:(_?OperatorQ | _?NumericQ)...]] /;
            Length[{before, after}] > 0 && OperatorQ[value] :>
                OperatorProduct[before, value, after, Sorted -> False]},
        UpValues[tag] = DeleteDuplicates[Append[UpValues[tag], rule]]
    ];
]];

opI /: HoldPattern[Dot[opI[j_, _?BooleanQ, ___], opI[k_, _?BooleanQ, ___]]] :=
    ScalarOperatorProduct[opI[j], opI[k]];

Superoperator /: HoldPattern[Dot[c:Except[_List].., s_Superoperator, op_?OperatorQ]] := Dot[c, s[op]];

(* These two subvalues implement (s1 . s2)[op]. The kernel only consults
   them when the product itself is applied as a function. *)
withHeads[{Dot}, Function[tag,
    tag /: HoldPattern[Dot[ss:Except[_List]?SuperoperatorQ..][op_?OperatorQ]] :=
        Fold[Function[{v, s}, s[v]], op, Reverse[{ss}]];
    tag /: HoldPattern[Dot[c:Except[_List].., s_Superoperator][op_?OperatorQ]] := Dot[c, s[op]];
]];

End[];
EndPackage[];
