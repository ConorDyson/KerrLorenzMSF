(* ::Package:: *)

(* ::Title:: *)
(*MetricReconstructRadiative Package*)


BeginPackage["MetricReconstructRadiative`",
{
"OrbitalData`",
"KerrOrbitalParameters`",
"SWSHdecomp`",
"SpinWeightedSpheroidalHarmonicsFT`",
"NGrid`"
}];
Needs["SpinWeightedSpheroidalHarmonics`"]
Needs["Teukolsky`"]
Print["Refactor attempt 2..."]


(* ::Chapter:: *)
(*Preamble*)


MetricReconstructRadiative::usage = "Outputs reconstructed metric for given config file"


(* ::Chapter:: *)
(*Private*)


Begin["`Private`"];


(* ::Section:: *)
(*Turn off some messages*)


Off[ClebschGordan::phy];
Off[ClebschGordan::tri];
Off[Infinity::indet];
Off[SpinWeightedSphericalHarmonicY::params];
Off[Power::infy];
Off[FrontEndObject::notavail];


(* ::Section:: *)
(*Helper Functions*)


Get\[CapitalLambda][l_,m_,a\[Omega]_]:=Module[{AA},
SpinWeightedSpheroidalEigenvalue[-2,l,m,a\[Omega]]
]; (* N.B. s = -2. *)


EchoT[lbl_][expr_]:=EchoTiming[expr,lbl]


(* ::Subsection:: *)
(*Teukolsky Equation stuff*)


(* ::Subsubsection::Closed:: *)
(*Teukolsky equation subs*)


(* ::Input::Initialization:: *)
DefTeukSubs[2,{Pm2_Symbol,Pp2_Symbol,Sm2_Symbol,Sp2_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[CapitalLambda]_},{r_Symbol,\[Theta]_Symbol}]:={Derivative[2][Pm2][r]->1/\[CapitalDelta][r]^2 (-K[r]^2 Pm2[r]-2 I K[r] Pm2[r] Derivative[1][\[CapitalDelta]][r]+\[CapitalDelta][r] (Pm2[r] (\[CapitalLambda]+6 I r \[Omega]+I Derivative[1][K][r])+Derivative[1][Pm2][r] Derivative[1][\[CapitalDelta]][r])),\!\(\*SuperscriptBox[\(Pm2\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[r]->1/\[CapitalDelta][r]^3 I (I K[r]^2 (2 (M-r) Pm2[r]+\[CapitalDelta][r] Derivative[1][Pm2][r])+4 K[r] (Pm2[r] (2 (M-r)^2+(-1+I r \[Omega]) \[CapitalDelta][r])+(M-r) \[CapitalDelta][r] Derivative[1][Pm2][r])+\[CapitalDelta][r] (8 \[Omega] Pm2[r] ((M-r) r+\[CapitalDelta][r])-I (2+\[CapitalLambda]+8 I r \[Omega]) \[CapitalDelta][r] Derivative[1][Pm2][r])),\!\(\*SuperscriptBox[\(Pm2\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[r]->1/\[CapitalDelta][r]^4 (8 I (-M+r) K[r]^3 Pm2[r]+K[r]^4 Pm2[r]-2 K[r]^2 (Pm2[r] (14 (M-r)^2+(\[CapitalLambda]+8 I r \[Omega]) \[CapitalDelta][r])+2 (M-r) \[CapitalDelta][r] Derivative[1][Pm2][r])+I \[CapitalDelta][r] (Pm2[r] (48 (M-r)^2 r \[Omega]+(-I \[CapitalLambda]^2+2 \[CapitalLambda] (-I+8 r \[Omega])+24 \[Omega] (M+r (-1+3 I r \[Omega]))) \[CapitalDelta][r])+16 \[Omega] \[CapitalDelta][r] ((M-r) r+\[CapitalDelta][r]) Derivative[1][Pm2][r])+4 I K[r] (Pm2[r] (12 (M-r)^3+2 (M-r) (-3+\[CapitalLambda]+11 I r \[Omega]) \[CapitalDelta][r]+I \[Omega] \[CapitalDelta][r]^2)+2 \[CapitalDelta][r] (2 (M-r)^2+(-1+I r \[Omega]) \[CapitalDelta][r]) Derivative[1][Pm2][r])),Derivative[2][Pp2][r]->1/\[CapitalDelta][r]^2 (-K[r]^2 Pp2[r]+2 I K[r] Pp2[r] Derivative[1][\[CapitalDelta]][r]+\[CapitalDelta][r] (Pp2[r] (\[CapitalLambda]-6 I r \[Omega]-I Derivative[1][K][r])+Derivative[1][Pp2][r] Derivative[1][\[CapitalDelta]][r])),\!\(\*SuperscriptBox[\(Pp2\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[r]->-(1/\[CapitalDelta][r]^3)I (-I K[r]^2 (2 (M-r) Pp2[r]+\[CapitalDelta][r] Derivative[1][Pp2][r])+4 K[r] (Pp2[r] (2 (M-r)^2+(-1-I r \[Omega]) \[CapitalDelta][r])+(M-r) \[CapitalDelta][r] Derivative[1][Pp2][r])+\[CapitalDelta][r] (8 \[Omega] Pp2[r] ((M-r) r+\[CapitalDelta][r])+(2 I+I \[CapitalLambda]+8 r \[Omega]) \[CapitalDelta][r] Derivative[1][Pp2][r])),\!\(\*SuperscriptBox[\(Pp2\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[r]->-(1/\[CapitalDelta][r]^4)I (8 (-M+r) K[r]^3 Pp2[r]+I K[r]^4 Pp2[r]-2 I K[r]^2 (Pp2[r] (14 (M-r)^2+(\[CapitalLambda]-8 I r \[Omega]) \[CapitalDelta][r])+2 (M-r) \[CapitalDelta][r] Derivative[1][Pp2][r])+\[CapitalDelta][r] (Pp2[r] (48 (M-r)^2 r \[Omega]+(I \[CapitalLambda]^2+2 \[CapitalLambda] (I+8 r \[Omega])+24 \[Omega] (M+r (-1-3 I r \[Omega]))) \[CapitalDelta][r])+16 \[Omega] \[CapitalDelta][r] ((M-r) r+\[CapitalDelta][r]) Derivative[1][Pp2][r])+4 K[r] (Pp2[r] (12 (M-r)^3+2 (M-r) (-3+\[CapitalLambda]-11 I r \[Omega]) \[CapitalDelta][r]-I \[Omega] \[CapitalDelta][r]^2)+2 \[CapitalDelta][r] (2 (M-r)^2+(-1-I r \[Omega]) \[CapitalDelta][r]) Derivative[1][Pp2][r])),Derivative[2][Sm2][\[Theta]]->(-1-\[CapitalLambda]-2 a m \[Omega]-4 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2-4 m Cot[\[Theta]] Csc[\[Theta]]+3 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Sm2[\[Theta]]-Cot[\[Theta]] Derivative[1][Sm2][\[Theta]],\!\(\*SuperscriptBox[\(Sm2\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Theta]]->(-Cot[\[Theta]]^3+8 m Cot[\[Theta]]^2 Csc[\[Theta]]+Cot[\[Theta]] (1+\[CapitalLambda]+2 a m \[Omega]+4 a \[Omega] Cos[\[Theta]]-(11+3 m^2) Csc[\[Theta]]^2)+(4 a \[Omega]+a^2 \[Omega]^2 Cos[\[Theta]]+4 m Csc[\[Theta]]^4) Sin[\[Theta]]) Sm2[\[Theta]]+(-2-\[CapitalLambda]-2 a m \[Omega]-4 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2-4 m Cot[\[Theta]] Csc[\[Theta]]+5 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Derivative[1][Sm2][\[Theta]],\!\(\*SuperscriptBox[\(Sm2\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[\[Theta]]->(-a^2 \[Omega]^2-8 m Cot[\[Theta]]^3 Csc[\[Theta]]^3-(1+\[CapitalLambda]+2 a m \[Omega]) Csc[\[Theta]]^4+(11+3 m^2) Csc[\[Theta]]^6+Cot[\[Theta]]^2 (a^2 \[Omega]^2+(25+6 m^2) Csc[\[Theta]]^4)-4 Cot[\[Theta]] (a \[Omega] Csc[\[Theta]]^3+7 m Csc[\[Theta]]^5)) Sin[\[Theta]]^2 Sm2[\[Theta]]+(4 m Cot[\[Theta]]^2 Csc[\[Theta]]-2 (6+m^2) Cot[\[Theta]] Csc[\[Theta]]^2+4 m Csc[\[Theta]]^3+2 a \[Omega] (2+a \[Omega] Cos[\[Theta]]) Sin[\[Theta]]) Derivative[1][Sm2][\[Theta]]+(-Cot[\[Theta]]^3+8 m Cot[\[Theta]]^2 Csc[\[Theta]]+Cot[\[Theta]] (1+\[CapitalLambda]+2 a m \[Omega]+4 a \[Omega] Cos[\[Theta]]-(11+3 m^2) Csc[\[Theta]]^2)+(4 a \[Omega]+a^2 \[Omega]^2 Cos[\[Theta]]+4 m Csc[\[Theta]]^4) Sin[\[Theta]]) Derivative[1][Sm2][\[Theta]]+(-2-\[CapitalLambda]-2 a m \[Omega]-4 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2-4 m Cot[\[Theta]] Csc[\[Theta]]+5 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) ((-1-\[CapitalLambda]-2 a m \[Omega]-4 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2-4 m Cot[\[Theta]] Csc[\[Theta]]+3 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Sm2[\[Theta]]-Cot[\[Theta]] Derivative[1][Sm2][\[Theta]]),Derivative[2][Sp2][\[Theta]]->(-1-\[CapitalLambda]-2 a m \[Omega]+4 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2+4 m Cot[\[Theta]] Csc[\[Theta]]+3 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Sp2[\[Theta]]-Cot[\[Theta]] Derivative[1][Sp2][\[Theta]],\!\(\*SuperscriptBox[\(Sp2\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Theta]]->(-Cot[\[Theta]]^3-8 m Cot[\[Theta]]^2 Csc[\[Theta]]+Cot[\[Theta]] (1+\[CapitalLambda]+2 a m \[Omega]-4 a \[Omega] Cos[\[Theta]]-(11+3 m^2) Csc[\[Theta]]^2)+(a^2 \[Omega]^2 Cos[\[Theta]]-4 (a \[Omega]+m Csc[\[Theta]]^4)) Sin[\[Theta]]) Sp2[\[Theta]]+(-2-\[CapitalLambda]-2 a m \[Omega]+4 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2+4 m Cot[\[Theta]] Csc[\[Theta]]+5 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Derivative[1][Sp2][\[Theta]],\!\(\*SuperscriptBox[\(Sp2\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[\[Theta]]->(-a^2 \[Omega]^2+8 m Cot[\[Theta]]^3 Csc[\[Theta]]^3-(1+\[CapitalLambda]+2 a m \[Omega]) Csc[\[Theta]]^4+(11+3 m^2) Csc[\[Theta]]^6+Cot[\[Theta]]^2 (a^2 \[Omega]^2+(25+6 m^2) Csc[\[Theta]]^4)+4 Cot[\[Theta]] (a \[Omega] Csc[\[Theta]]^3+7 m Csc[\[Theta]]^5)) Sin[\[Theta]]^2 Sp2[\[Theta]]-2 (-a^2 \[Omega]^2 Cos[\[Theta]]+2 m Cot[\[Theta]]^2 Csc[\[Theta]]^2+(6+m^2) Cot[\[Theta]] Csc[\[Theta]]^3+2 (a \[Omega]+m Csc[\[Theta]]^4)) Sin[\[Theta]] Derivative[1][Sp2][\[Theta]]+(-Cot[\[Theta]]^3-8 m Cot[\[Theta]]^2 Csc[\[Theta]]+Cot[\[Theta]] (1+\[CapitalLambda]+2 a m \[Omega]-4 a \[Omega] Cos[\[Theta]]-(11+3 m^2) Csc[\[Theta]]^2)+(a^2 \[Omega]^2 Cos[\[Theta]]-4 (a \[Omega]+m Csc[\[Theta]]^4)) Sin[\[Theta]]) Derivative[1][Sp2][\[Theta]]+(-2-\[CapitalLambda]-2 a m \[Omega]+4 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2+4 m Cot[\[Theta]] Csc[\[Theta]]+5 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) ((-1-\[CapitalLambda]-2 a m \[Omega]+4 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2+4 m Cot[\[Theta]] Csc[\[Theta]]+3 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Sp2[\[Theta]]-Cot[\[Theta]] Derivative[1][Sp2][\[Theta]])}


DefTeukSubs[1,{Pm1_Symbol,Pp1_Symbol,Sm1_Symbol,Sp1_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[Lambda]_},{r_Symbol,\[Theta]_Symbol}]:={Derivative[2][Pm1][r]->(Pm1[r] (-K[r]^2+\[CapitalDelta][r] (\[Lambda]+2 I r \[Omega]+I Derivative[1][K][r])-I K[r] Derivative[1][\[CapitalDelta]][r]))/\[CapitalDelta][r]^2,\!\(\*SuperscriptBox[\(Pm1\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[r]->1/\[CapitalDelta][r]^3 I (I K[r]^2 (4 (M-r) Pm1[r]+\[CapitalDelta][r] Derivative[1][Pm1][r])+2 K[r] (Pm1[r] (4 (M-r)^2+(-1+2 I r \[Omega]) \[CapitalDelta][r])+(M-r) \[CapitalDelta][r] Derivative[1][Pm1][r])+\[CapitalDelta][r] (2 Pm1[r] ((M-r) (-I \[Lambda]+6 r \[Omega])+2 \[Omega] \[CapitalDelta][r])+(-I \[Lambda]+4 r \[Omega]) \[CapitalDelta][r] Derivative[1][Pm1][r])),\!\(\*SuperscriptBox[\(Pm1\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[r]->1/\[CapitalDelta][r]^4 (4 I (-M+r) K[r]^3 Pm1[r]+K[r]^4 Pm1[r]-2 K[r]^2 (Pm1[r] (14 (M-r)^2+(-2+\[Lambda]+4 I r \[Omega]) \[CapitalDelta][r])+4 (M-r) \[CapitalDelta][r] Derivative[1][Pm1][r])+I \[CapitalDelta][r] (Pm1[r] (8 (M-r)^2 (-I \[Lambda]+8 r \[Omega])+(-I \[Lambda]^2+\[Lambda] (2 I+8 r \[Omega])+4 \[Omega] (5 M-9 r+6 I r^2 \[Omega])) \[CapitalDelta][r])+4 \[CapitalDelta][r] ((M-r) (-I \[Lambda]+6 r \[Omega])+2 \[Omega] \[CapitalDelta][r]) Derivative[1][Pm1][r])+4 I K[r] (Pm1[r] (12 (M-r)^3+(M-r) (-6+\[Lambda]+12 I r \[Omega]) \[CapitalDelta][r]+I \[Omega] \[CapitalDelta][r]^2)+\[CapitalDelta][r] (4 (M-r)^2+(-1+2 I r \[Omega]) \[CapitalDelta][r]) Derivative[1][Pm1][r])),Derivative[2][Pp1][r]->-((Pp1[r] (K[r]^2+I \[CapitalDelta][r] (I \[Lambda]+2 r \[Omega]+Derivative[1][K][r])-I K[r] Derivative[1][\[CapitalDelta]][r]))/\[CapitalDelta][r]^2),\!\(\*SuperscriptBox[\(Pp1\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[r]->-(1/\[CapitalDelta][r]^3)I (-I K[r]^2 (4 (M-r) Pp1[r]+\[CapitalDelta][r] Derivative[1][Pp1][r])+2 K[r] (Pp1[r] (4 (M-r)^2+(-1-2 I r \[Omega]) \[CapitalDelta][r])+(M-r) \[CapitalDelta][r] Derivative[1][Pp1][r])+\[CapitalDelta][r] (2 Pp1[r] ((M-r) (I \[Lambda]+6 r \[Omega])+2 \[Omega] \[CapitalDelta][r])+(I \[Lambda]+4 r \[Omega]) \[CapitalDelta][r] Derivative[1][Pp1][r])),\!\(\*SuperscriptBox[\(Pp1\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[r]->-(1/\[CapitalDelta][r]^4)I (4 (-M+r) K[r]^3 Pp1[r]+I K[r]^4 Pp1[r]-2 I K[r]^2 (Pp1[r] (14 (M-r)^2+(-2+\[Lambda]-4 I r \[Omega]) \[CapitalDelta][r])+4 (M-r) \[CapitalDelta][r] Derivative[1][Pp1][r])+\[CapitalDelta][r] (Pp1[r] (8 (M-r)^2 (I \[Lambda]+8 r \[Omega])+(I \[Lambda]^2+\[Lambda] (-2 I+8 r \[Omega])+4 \[Omega] (5 M-9 r-6 I r^2 \[Omega])) \[CapitalDelta][r])+4 \[CapitalDelta][r] ((M-r) (I \[Lambda]+6 r \[Omega])+2 \[Omega] \[CapitalDelta][r]) Derivative[1][Pp1][r])+4 K[r] (Pp1[r] (12 (M-r)^3+(M-r) (-6+\[Lambda]-12 I r \[Omega]) \[CapitalDelta][r]-I \[Omega] \[CapitalDelta][r]^2)+\[CapitalDelta][r] (4 (M-r)^2+(-1-2 I r \[Omega]) \[CapitalDelta][r]) Derivative[1][Pp1][r])),Derivative[2][Sm1][\[Theta]]->(-\[Lambda]-2 a m \[Omega]-2 a \[Omega] Cos[\[Theta]]-2 m Cot[\[Theta]] Csc[\[Theta]]+Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Sm1[\[Theta]]-Cot[\[Theta]] Derivative[1][Sm1][\[Theta]],\!\(\*SuperscriptBox[\(Sm1\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Theta]]->(4 m Cot[\[Theta]]^2 Csc[\[Theta]]+Cot[\[Theta]] (\[Lambda]+2 a m \[Omega]+2 a \[Omega] Cos[\[Theta]]-3 (1+m^2) Csc[\[Theta]]^2)+(2 a \[Omega]+a^2 \[Omega]^2 Cos[\[Theta]]+2 m Csc[\[Theta]]^4) Sin[\[Theta]]) Sm1[\[Theta]]+(-\[Lambda]-2 a m \[Omega]-2 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2-2 m Cot[\[Theta]] Csc[\[Theta]]+2 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Derivative[1][Sm1][\[Theta]],\!\(\*SuperscriptBox[\(Sm1\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[\[Theta]]->(-a^2 \[Omega]^2-4 m Cot[\[Theta]]^3 Csc[\[Theta]]^3-(\[Lambda]+2 a m \[Omega]) Csc[\[Theta]]^4+3 (1+m^2) Csc[\[Theta]]^6+Cot[\[Theta]]^2 (a^2 \[Omega]^2+6 (1+m^2) Csc[\[Theta]]^4)-2 Cot[\[Theta]] (a \[Omega] Csc[\[Theta]]^3+7 m Csc[\[Theta]]^5)) Sin[\[Theta]]^2 Sm1[\[Theta]]+2 (a \[Omega]+a^2 \[Omega]^2 Cos[\[Theta]]+m Cot[\[Theta]]^2 Csc[\[Theta]]^2-(3+m^2) Cot[\[Theta]] Csc[\[Theta]]^3+m Csc[\[Theta]]^4) Sin[\[Theta]] Derivative[1][Sm1][\[Theta]]+(4 m Cot[\[Theta]]^2 Csc[\[Theta]]+Cot[\[Theta]] (\[Lambda]+2 a m \[Omega]+2 a \[Omega] Cos[\[Theta]]-3 (1+m^2) Csc[\[Theta]]^2)+(2 a \[Omega]+a^2 \[Omega]^2 Cos[\[Theta]]+2 m Csc[\[Theta]]^4) Sin[\[Theta]]) Derivative[1][Sm1][\[Theta]]+(-\[Lambda]-2 a m \[Omega]-2 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2-2 m Cot[\[Theta]] Csc[\[Theta]]+2 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) ((-\[Lambda]-2 a m \[Omega]-2 a \[Omega] Cos[\[Theta]]-2 m Cot[\[Theta]] Csc[\[Theta]]+Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Sm1[\[Theta]]-Cot[\[Theta]] Derivative[1][Sm1][\[Theta]]),Derivative[2][Sp1][\[Theta]]->(-\[Lambda]-2 a m \[Omega]+2 a \[Omega] Cos[\[Theta]]+2 m Cot[\[Theta]] Csc[\[Theta]]+Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Sp1[\[Theta]]-Cot[\[Theta]] Derivative[1][Sp1][\[Theta]],\!\(\*SuperscriptBox[\(Sp1\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Theta]]->(-4 m Cot[\[Theta]]^2 Csc[\[Theta]]+Cot[\[Theta]] (\[Lambda]+2 a m \[Omega]-2 a \[Omega] Cos[\[Theta]]-3 (1+m^2) Csc[\[Theta]]^2)+(a^2 \[Omega]^2 Cos[\[Theta]]-2 (a \[Omega]+m Csc[\[Theta]]^4)) Sin[\[Theta]]) Sp1[\[Theta]]+(-\[Lambda]-2 a m \[Omega]+2 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2+2 m Cot[\[Theta]] Csc[\[Theta]]+2 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Derivative[1][Sp1][\[Theta]],\!\(\*SuperscriptBox[\(Sp1\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[\[Theta]]->(-a^2 \[Omega]^2+4 m Cot[\[Theta]]^3 Csc[\[Theta]]^3-(\[Lambda]+2 a m \[Omega]) Csc[\[Theta]]^4+3 (1+m^2) Csc[\[Theta]]^6+Cot[\[Theta]]^2 (a^2 \[Omega]^2+6 (1+m^2) Csc[\[Theta]]^4)+2 Cot[\[Theta]] (a \[Omega] Csc[\[Theta]]^3+7 m Csc[\[Theta]]^5)) Sin[\[Theta]]^2 Sp1[\[Theta]]-2 (a \[Omega]-a^2 \[Omega]^2 Cos[\[Theta]]+m Cot[\[Theta]]^2 Csc[\[Theta]]^2+(3+m^2) Cot[\[Theta]] Csc[\[Theta]]^3+m Csc[\[Theta]]^4) Sin[\[Theta]] Derivative[1][Sp1][\[Theta]]+(-4 m Cot[\[Theta]]^2 Csc[\[Theta]]+Cot[\[Theta]] (\[Lambda]+2 a m \[Omega]-2 a \[Omega] Cos[\[Theta]]-3 (1+m^2) Csc[\[Theta]]^2)+(a^2 \[Omega]^2 Cos[\[Theta]]-2 (a \[Omega]+m Csc[\[Theta]]^4)) Sin[\[Theta]]) Derivative[1][Sp1][\[Theta]]+(-\[Lambda]-2 a m \[Omega]+2 a \[Omega] Cos[\[Theta]]+Cot[\[Theta]]^2+2 m Cot[\[Theta]] Csc[\[Theta]]+2 Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) ((-\[Lambda]-2 a m \[Omega]+2 a \[Omega] Cos[\[Theta]]+2 m Cot[\[Theta]] Csc[\[Theta]]+Csc[\[Theta]]^2+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Sp1[\[Theta]]-Cot[\[Theta]] Derivative[1][Sp1][\[Theta]])}


DefTeukSubs[0,{P0_Symbol,S0_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[Lambda]0_},{r_Symbol,\[Theta]_Symbol}]:={Derivative[2][P0][r]->(-K[r]^2 P0[r]+\[CapitalDelta][r] (\[Lambda]0 P0[r]+2 (M-r) Derivative[1][P0][r]))/\[CapitalDelta][r]^2,\!\(\*SuperscriptBox[\(P0\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[r]->1/\[CapitalDelta][r]^3 (-4 r \[Omega] K[r] P0[r] \[CapitalDelta][r]-K[r]^2 (6 (M-r) P0[r]+\[CapitalDelta][r] Derivative[1][P0][r])+\[CapitalDelta][r] (4 (M-r) \[Lambda]0 P0[r]+(8 (M-r)^2+(-2+\[Lambda]0) \[CapitalDelta][r]) Derivative[1][P0][r])),\!\(\*SuperscriptBox[\(P0\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[r]->1/\[CapitalDelta][r]^4 (K[r]^4 P0[r]-2 K[r]^2 (P0[r] (22 (M-r)^2+(-4+\[Lambda]0) \[CapitalDelta][r])+6 (M-r) \[CapitalDelta][r] Derivative[1][P0][r])-4 \[Omega] K[r] \[CapitalDelta][r] (P0[r] (10 (M-r) r+\[CapitalDelta][r])+2 r \[CapitalDelta][r] Derivative[1][P0][r])+\[CapitalDelta][r] (P0[r] (24 (M-r)^2 \[Lambda]0+(-6 \[Lambda]0+\[Lambda]0^2-8 r^2 \[Omega]^2) \[CapitalDelta][r])+8 (M-r) (6 (M-r)^2+(-3+\[Lambda]0) \[CapitalDelta][r]) Derivative[1][P0][r])),Derivative[2][S0][\[Theta]]->S0[\[Theta]] (-\[Lambda]0-2 a m \[Omega]+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2)-Cot[\[Theta]] Derivative[1][S0][\[Theta]],\!\(\*SuperscriptBox[\(S0\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Theta]]->S0[\[Theta]] (Cot[\[Theta]] (\[Lambda]0+2 a m \[Omega]-3 m^2 Csc[\[Theta]]^2)+a^2 \[Omega]^2 Cos[\[Theta]] Sin[\[Theta]])+(-\[Lambda]0-2 a m \[Omega]+Cot[\[Theta]]^2+(1+m^2) Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) Derivative[1][S0][\[Theta]],\!\(\*SuperscriptBox[\(S0\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[\[Theta]]->S0[\[Theta]] (a^2 \[Omega]^2 Cos[\[Theta]]^2-(\[Lambda]0+2 a m \[Omega]-6 m^2 Cot[\[Theta]]^2) Csc[\[Theta]]^2+3 m^2 Csc[\[Theta]]^4-a^2 \[Omega]^2 Sin[\[Theta]]^2)-2 (-a^2 \[Omega]^2 Cos[\[Theta]]+(2+m^2) Cot[\[Theta]] Csc[\[Theta]]^3) Sin[\[Theta]] Derivative[1][S0][\[Theta]]+(Cot[\[Theta]] (\[Lambda]0+2 a m \[Omega]-3 m^2 Csc[\[Theta]]^2)+a^2 \[Omega]^2 Cos[\[Theta]] Sin[\[Theta]]) Derivative[1][S0][\[Theta]]+(-\[Lambda]0-2 a m \[Omega]+Cot[\[Theta]]^2+(1+m^2) Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2) (S0[\[Theta]] (-\[Lambda]0-2 a m \[Omega]+m^2 Csc[\[Theta]]^2+a^2 \[Omega]^2 Sin[\[Theta]]^2)-Cot[\[Theta]] Derivative[1][S0][\[Theta]])}


DefTeukSubs["\[Kappa]0",{\[Kappa]0_Symbol,h0_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[Lambda]0_},{r_Symbol,\[Theta]_Symbol}]:={Derivative[2][\[Kappa]0][r]->((r^2+a^2 Cos[\[Theta]]^2) h0[r] \[CapitalDelta][r]-2 K[r]^2 \[Kappa]0[r]+2 \[CapitalDelta][r] (\[Lambda]0 \[Kappa]0[r]+2 (M-r) Derivative[1][\[Kappa]0][r]))/(2 \[CapitalDelta][r]^2),\!\(\*SuperscriptBox[\(\[Kappa]0\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[r]->1/(2 \[CapitalDelta][r]^3) (2 h0[r] \[CapitalDelta][r] (2 (M-r) (r^2+a^2 Cos[\[Theta]]^2)+r \[CapitalDelta][r])-8 r \[Omega] K[r] \[CapitalDelta][r] \[Kappa]0[r]-2 K[r]^2 (6 (M-r) \[Kappa]0[r]+\[CapitalDelta][r] Derivative[1][\[Kappa]0][r])+\[CapitalDelta][r] (8 (M-r) \[Lambda]0 \[Kappa]0[r]+16 (M-r)^2 Derivative[1][\[Kappa]0][r]+\[CapitalDelta][r] ((r^2+a^2 Cos[\[Theta]]^2) Derivative[1][h0][r]+2 (-2+\[Lambda]0) Derivative[1][\[Kappa]0][r]))),\!\(\*SuperscriptBox[\(\[Kappa]0\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[r]->1/(2 \[CapitalDelta][r]^4) (h0[r] \[CapitalDelta][r] (12 (M-r)^2 (a^2+2 r^2+a^2 Cos[2 \[Theta]])-(r^2+a^2 Cos[\[Theta]]^2) K[r]^2+(r (12 M+r (-18+\[Lambda]0))+a^2 (-6+\[Lambda]0) Cos[\[Theta]]^2) \[CapitalDelta][r]+2 \[CapitalDelta][r]^2)+2 K[r]^4 \[Kappa]0[r]-4 K[r]^2 ((22 (M-r)^2+(-4+\[Lambda]0) \[CapitalDelta][r]) \[Kappa]0[r]+6 (M-r) \[CapitalDelta][r] Derivative[1][\[Kappa]0][r])-8 \[Omega] K[r] \[CapitalDelta][r] ((10 (M-r) r+\[CapitalDelta][r]) \[Kappa]0[r]+2 r \[CapitalDelta][r] Derivative[1][\[Kappa]0][r])+\[CapitalDelta][r] (2 (24 (M-r)^2 \[Lambda]0+(-6 \[Lambda]0+\[Lambda]0^2-8 r^2 \[Omega]^2) \[CapitalDelta][r]) \[Kappa]0[r]+96 (M-r)^3 Derivative[1][\[Kappa]0][r]+(M-r) \[CapitalDelta][r] (3 (a^2+2 r^2+a^2 Cos[2 \[Theta]]) Derivative[1][h0][r]+16 (-3+\[Lambda]0) Derivative[1][\[Kappa]0][r])+\[CapitalDelta][r]^2 (4 r Derivative[1][h0][r]+(r^2+a^2 Cos[\[Theta]]^2) Derivative[2][h0][r])))}


(* ::Subsubsection::Closed:: *)
(*Teukolksy Starabinskii identities*)


DefTSIP[2,"pm",{Pm2_Symbol,Pp2_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[CapitalLambda]_,CC_},{r_Symbol}]:={Pp2[r]->1/(CC \[CapitalDelta][r]^2) (8 K[r]^4 Pm2[r]+8 K[r]^2 Pm2[r] ((M-r)^2-(1+\[CapitalLambda]+3 I r \[Omega]) \[CapitalDelta][r])+8 I K[r]^3 \[CapitalDelta][r] Derivative[1][Pm2][r]+\[CapitalDelta][r]^2 ((\[CapitalLambda]^2+\[CapitalLambda] (2+4 I r \[Omega])+12 \[Omega] (I M+r^2 \[Omega])) Pm2[r]+8 I \[Omega] ((M-r) r+\[CapitalDelta][r]) Derivative[1][Pm2][r])-4 I K[r] \[CapitalDelta][r] (Pm2[r] ((M-r) (\[CapitalLambda]+8 I r \[Omega])+5 I \[Omega] \[CapitalDelta][r])+(-2 (M-r)^2+(2+\[CapitalLambda]) \[CapitalDelta][r]) Derivative[1][Pm2][r])),Derivative[1][Pp2][r]->-(1/(CC \[CapitalDelta][r]^3))I (8 K[r]^5 Pm2[r]+4 K[r]^3 Pm2[r] (2 (M-r)^2-(2+3 \[CapitalLambda]) \[CapitalDelta][r])+8 I K[r]^4 \[CapitalDelta][r] Derivative[1][Pm2][r]+\[CapitalDelta][r]^3 (-12 \[CapitalLambda] \[Omega] Pm2[r]+(I \[CapitalLambda]^2+\[CapitalLambda] (2 I+4 r \[Omega])+12 \[Omega] (M+I r^2 \[Omega])) Derivative[1][Pm2][r])+4 K[r] \[CapitalDelta][r]^2 ((\[CapitalLambda]+\[CapitalLambda]^2+24 r^2 \[Omega]^2) Pm2[r]+(-((M-r) (\[CapitalLambda]-8 I r \[Omega]))+5 I \[Omega] \[CapitalDelta][r]) Derivative[1][Pm2][r])+8 K[r]^2 \[CapitalDelta][r] (\[Omega] Pm2[r] (7 (M-r) r+4 \[CapitalDelta][r])+I ((M-r)^2-(1+\[CapitalLambda]-3 I r \[Omega]) \[CapitalDelta][r]) Derivative[1][Pm2][r]))}


DefTSIP[2,"mp",{Pm2_Symbol,Pp2_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[CapitalLambda]_,Cstar_},{r_Symbol}]:={Pm2[r]->1/(Cstar \[CapitalDelta][r]^2) (8 K[r]^4 Pp2[r]+8 K[r]^2 Pp2[r] ((M-r)^2-(1+\[CapitalLambda]-3 I r \[Omega]) \[CapitalDelta][r])-8 I K[r]^3 \[CapitalDelta][r] Derivative[1][Pp2][r]-I \[CapitalDelta][r]^2 ((I \[CapitalLambda]^2+\[CapitalLambda] (2 I+4 r \[Omega])+12 \[Omega] (M+I r^2 \[Omega])) Pp2[r]+8 \[Omega] ((M-r) r+\[CapitalDelta][r]) Derivative[1][Pp2][r])-4 I K[r] \[CapitalDelta][r] (Pp2[r] (-((M-r) (\[CapitalLambda]-8 I r \[Omega]))+5 I \[Omega] \[CapitalDelta][r])+(2 (M-r)^2-(2+\[CapitalLambda]) \[CapitalDelta][r]) Derivative[1][Pp2][r])),Derivative[1][Pm2][r]->1/(Cstar \[CapitalDelta][r]^3) I (8 K[r]^5 Pp2[r]+4 K[r]^3 Pp2[r] (2 (M-r)^2-(2+3 \[CapitalLambda]) \[CapitalDelta][r])-8 I K[r]^4 \[CapitalDelta][r] Derivative[1][Pp2][r]+\[CapitalDelta][r]^3 (-12 \[CapitalLambda] \[Omega] Pp2[r]-I (\[CapitalLambda]^2+\[CapitalLambda] (2+4 I r \[Omega])+12 \[Omega] (I M+r^2 \[Omega])) Derivative[1][Pp2][r])+4 K[r] \[CapitalDelta][r]^2 ((\[CapitalLambda]+\[CapitalLambda]^2+24 r^2 \[Omega]^2) Pp2[r]+(-((M-r) (\[CapitalLambda]+8 I r \[Omega]))-5 I \[Omega] \[CapitalDelta][r]) Derivative[1][Pp2][r])+8 K[r]^2 \[CapitalDelta][r] (\[Omega] Pp2[r] (7 (M-r) r+4 \[CapitalDelta][r])-I ((M-r)^2-(1+\[CapitalLambda]+3 I r \[Omega]) \[CapitalDelta][r]) Derivative[1][Pp2][r]))}


DefTSIS[2,"pm",{Sm2_Symbol,Sp2_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[CapitalLambda]_,AA_},{\[Theta]_Symbol}]:={Sp2[\[Theta]]->1/AA ((2+3 \[CapitalLambda]+\[CapitalLambda]^2+6 a m \[Omega]+16 a m \[CapitalLambda] \[Omega]+10 a^2 \[Omega]^2+48 a^2 m^2 \[Omega]^2+2 a^2 \[Omega]^2 Cos[\[Theta]]^2-2 Cot[\[Theta]]^4+8 m Cot[\[Theta]]^3 Csc[\[Theta]]-4 Csc[\[Theta]]^2-9 m^2 Csc[\[Theta]]^2-\[CapitalLambda] Csc[\[Theta]]^2-8 m^2 \[CapitalLambda] Csc[\[Theta]]^2-10 a m \[Omega] Csc[\[Theta]]^2-32 a m^3 \[Omega] Csc[\[Theta]]^2+2 Csc[\[Theta]]^4+m^2 Csc[\[Theta]]^4+8 m^4 Csc[\[Theta]]^4+Cot[\[Theta]]^2 (\[CapitalLambda]-6 a m \[Omega]-9 m^2 Csc[\[Theta]]^2)+Cot[\[Theta]] Csc[\[Theta]] (4 m (2+\[CapitalLambda])-9 a \[Omega]+24 a m^2 \[Omega]-8 m (-1+2 m^2) Csc[\[Theta]]^2)+2 a^2 \[Omega]^2 Sin[\[Theta]]^2-8 a^2 \[CapitalLambda] \[Omega]^2 Sin[\[Theta]]^2-32 a^3 m \[Omega]^3 Sin[\[Theta]]^2+8 a^4 \[Omega]^4 Sin[\[Theta]]^4+a \[Omega] Cos[\[Theta]] (9+9 Cot[\[Theta]]^2-8 a^2 \[Omega]^2 Sin[\[Theta]]^2)) Sm2[\[Theta]]+4 (-a \[Omega] Cos[\[Theta]] Cot[\[Theta]]+(m (2+\[CapitalLambda])+a \[Omega]+6 a m^2 \[Omega]+2 m Cot[\[Theta]]^2) Csc[\[Theta]]-2 m^3 Csc[\[Theta]]^3+a \[Omega] Sin[\[Theta]] (-1-\[CapitalLambda]-6 a m \[Omega]+2 a^2 \[Omega]^2 Sin[\[Theta]]^2)) Derivative[1][Sm2][\[Theta]]),Derivative[1][Sp2][\[Theta]]->1/AA ((-4 (4 a (1+\[CapitalLambda]) \[Omega]+a m^2 (11+9 \[CapitalLambda]) \[Omega]+20 a^2 m^3 \[Omega]^2+m (2+3 \[CapitalLambda]+\[CapitalLambda]^2+20 a^2 \[Omega]^2)) Csc[\[Theta]]+(8 m (2+\[CapitalLambda])+4 m^3 (4+3 \[CapitalLambda])+21 a \[Omega]+52 a m^2 \[Omega]+40 a m^4 \[Omega]) Csc[\[Theta]]^3-8 (m+m^3+m^5) Csc[\[Theta]]^5+Cot[\[Theta]]^2 Csc[\[Theta]] (-8 m (1+\[CapitalLambda])-17 a \[Omega]-20 a m^2 \[Omega]+24 m (-1+2 m^2) Csc[\[Theta]]^2)-2 Cot[\[Theta]]^3 (8 a m \[Omega]+2 a \[Omega] Cos[\[Theta]]+(-4+7 m^2) Csc[\[Theta]]^2)-Cot[\[Theta]] (a \[Omega] (9+8 a m \[Omega]) Cos[\[Theta]]-16 a^2 \[Omega]^2 Cos[\[Theta]]^2+2 (8 a \[Omega] (m+a \[Omega])+(-4+7 m^2-8 a m \[Omega]) Csc[\[Theta]]^2+(4-7 m^2) Csc[\[Theta]]^4))+a \[Omega] ((-5+8 \[CapitalLambda]+4 \[CapitalLambda]^2+40 a m \[Omega]+36 a m \[CapitalLambda] \[Omega]+19 a^2 \[Omega]^2+80 a^2 m^2 \[Omega]^2-3 a^2 \[Omega]^2 Cos[2 \[Theta]]) Sin[\[Theta]]-2 a^2 \[Omega]^2 (-1+6 \[CapitalLambda]+20 a m \[Omega]) Sin[\[Theta]]^3+8 a^4 \[Omega]^4 Sin[\[Theta]]^5+8 a \[Omega] Sin[2 \[Theta]])) Sm2[\[Theta]]+(2+3 \[CapitalLambda]+\[CapitalLambda]^2+6 a m \[Omega]+16 a m \[CapitalLambda] \[Omega]+10 a^2 \[Omega]^2+48 a^2 m^2 \[Omega]^2+2 a^2 \[Omega]^2 Cos[\[Theta]]^2-2 Cot[\[Theta]]^4-8 m Cot[\[Theta]]^3 Csc[\[Theta]]-4 Csc[\[Theta]]^2-9 m^2 Csc[\[Theta]]^2-\[CapitalLambda] Csc[\[Theta]]^2-8 m^2 \[CapitalLambda] Csc[\[Theta]]^2-10 a m \[Omega] Csc[\[Theta]]^2-32 a m^3 \[Omega] Csc[\[Theta]]^2+2 Csc[\[Theta]]^4+m^2 Csc[\[Theta]]^4+8 m^4 Csc[\[Theta]]^4+Cot[\[Theta]]^2 (\[CapitalLambda]-6 a m \[Omega]-9 m^2 Csc[\[Theta]]^2)+Cot[\[Theta]] Csc[\[Theta]] (-4 m (2+\[CapitalLambda])-13 a \[Omega]-24 a m^2 \[Omega]+8 m (-1+2 m^2) Csc[\[Theta]]^2)+2 a^2 \[Omega]^2 Sin[\[Theta]]^2-8 a^2 \[CapitalLambda] \[Omega]^2 Sin[\[Theta]]^2-32 a^3 m \[Omega]^3 Sin[\[Theta]]^2+8 a^4 \[Omega]^4 Sin[\[Theta]]^4+a \[Omega] Cos[\[Theta]] (13+13 Cot[\[Theta]]^2+8 a^2 \[Omega]^2 Sin[\[Theta]]^2)) Derivative[1][Sm2][\[Theta]])}


DefTSIS[2,"mp",{Sm2_Symbol,Sp2_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[CapitalLambda]_,AA_},{\[Theta]_Symbol}]:={Sm2[\[Theta]]->1/AA ((2+3 \[CapitalLambda]+\[CapitalLambda]^2+6 a m \[Omega]+16 a m \[CapitalLambda] \[Omega]+10 a^2 \[Omega]^2+48 a^2 m^2 \[Omega]^2+2 a^2 \[Omega]^2 Cos[\[Theta]]^2-2 Cot[\[Theta]]^4-8 m Cot[\[Theta]]^3 Csc[\[Theta]]-4 Csc[\[Theta]]^2-9 m^2 Csc[\[Theta]]^2-\[CapitalLambda] Csc[\[Theta]]^2-8 m^2 \[CapitalLambda] Csc[\[Theta]]^2-10 a m \[Omega] Csc[\[Theta]]^2-32 a m^3 \[Omega] Csc[\[Theta]]^2+2 Csc[\[Theta]]^4+m^2 Csc[\[Theta]]^4+8 m^4 Csc[\[Theta]]^4+Cot[\[Theta]]^2 (\[CapitalLambda]-6 a m \[Omega]-9 m^2 Csc[\[Theta]]^2)+Cot[\[Theta]] Csc[\[Theta]] (-4 m (2+\[CapitalLambda])+9 a \[Omega]-24 a m^2 \[Omega]+8 m (-1+2 m^2) Csc[\[Theta]]^2)+2 a^2 \[Omega]^2 Sin[\[Theta]]^2-8 a^2 \[CapitalLambda] \[Omega]^2 Sin[\[Theta]]^2-32 a^3 m \[Omega]^3 Sin[\[Theta]]^2+8 a^4 \[Omega]^4 Sin[\[Theta]]^4+a \[Omega] Cos[\[Theta]] (-9-9 Cot[\[Theta]]^2+8 a^2 \[Omega]^2 Sin[\[Theta]]^2)) Sp2[\[Theta]]+4 (a \[Omega] Cos[\[Theta]] Cot[\[Theta]]-(m (2+\[CapitalLambda])+a \[Omega]+6 a m^2 \[Omega]+2 m Cot[\[Theta]]^2) Csc[\[Theta]]+2 m^3 Csc[\[Theta]]^3+a \[Omega] Sin[\[Theta]] (1+\[CapitalLambda]+6 a m \[Omega]-2 a^2 \[Omega]^2 Sin[\[Theta]]^2)) Derivative[1][Sp2][\[Theta]]),Derivative[1][Sm2][\[Theta]]->1/AA ((4 (4 a (1+\[CapitalLambda]) \[Omega]+a m^2 (11+9 \[CapitalLambda]) \[Omega]+20 a^2 m^3 \[Omega]^2+m (2+3 \[CapitalLambda]+\[CapitalLambda]^2+20 a^2 \[Omega]^2)) Csc[\[Theta]]-(8 m (2+\[CapitalLambda])+4 m^3 (4+3 \[CapitalLambda])+21 a \[Omega]+52 a m^2 \[Omega]+40 a m^4 \[Omega]) Csc[\[Theta]]^3+8 (m+m^3+m^5) Csc[\[Theta]]^5+Cot[\[Theta]]^2 Csc[\[Theta]] (8 m (1+\[CapitalLambda])+17 a \[Omega]+20 a m^2 \[Omega]-24 m (-1+2 m^2) Csc[\[Theta]]^2)-2 Cot[\[Theta]]^3 (8 a m \[Omega]-2 a \[Omega] Cos[\[Theta]]+(-4+7 m^2) Csc[\[Theta]]^2)+Cot[\[Theta]] (a \[Omega] (9+8 a m \[Omega]) Cos[\[Theta]]+16 a^2 \[Omega]^2 Cos[\[Theta]]^2-2 (8 a \[Omega] (m+a \[Omega])+(-4+7 m^2-8 a m \[Omega]) Csc[\[Theta]]^2+(4-7 m^2) Csc[\[Theta]]^4))+a \[Omega] (-((-5+8 \[CapitalLambda]+4 \[CapitalLambda]^2+40 a m \[Omega]+36 a m \[CapitalLambda] \[Omega]+22 a^2 \[Omega]^2+80 a^2 m^2 \[Omega]^2-6 a^2 \[Omega]^2 Cos[2 \[Theta]]) Sin[\[Theta]])+4 a^2 \[Omega]^2 (1+3 \[CapitalLambda]+10 a m \[Omega]) Sin[\[Theta]]^3-8 a^4 \[Omega]^4 Sin[\[Theta]]^5+8 a \[Omega] Sin[2 \[Theta]])) Sp2[\[Theta]]+(2+3 \[CapitalLambda]+\[CapitalLambda]^2+6 a m \[Omega]+16 a m \[CapitalLambda] \[Omega]+10 a^2 \[Omega]^2+48 a^2 m^2 \[Omega]^2+2 a^2 \[Omega]^2 Cos[\[Theta]]^2-2 Cot[\[Theta]]^4+8 m Cot[\[Theta]]^3 Csc[\[Theta]]-4 Csc[\[Theta]]^2-9 m^2 Csc[\[Theta]]^2-\[CapitalLambda] Csc[\[Theta]]^2-8 m^2 \[CapitalLambda] Csc[\[Theta]]^2-10 a m \[Omega] Csc[\[Theta]]^2-32 a m^3 \[Omega] Csc[\[Theta]]^2+2 Csc[\[Theta]]^4+m^2 Csc[\[Theta]]^4+8 m^4 Csc[\[Theta]]^4+Cot[\[Theta]]^2 (\[CapitalLambda]-6 a m \[Omega]-9 m^2 Csc[\[Theta]]^2)+Cot[\[Theta]] Csc[\[Theta]] (4 m (2+\[CapitalLambda])+13 a \[Omega]+24 a m^2 \[Omega]-8 m (-1+2 m^2) Csc[\[Theta]]^2)+2 a^2 \[Omega]^2 Sin[\[Theta]]^2-8 a^2 \[CapitalLambda] \[Omega]^2 Sin[\[Theta]]^2-32 a^3 m \[Omega]^3 Sin[\[Theta]]^2+8 a^4 \[Omega]^4 Sin[\[Theta]]^4-a \[Omega] Cos[\[Theta]] (13+13 Cot[\[Theta]]^2+8 a^2 \[Omega]^2 Sin[\[Theta]]^2)) Derivative[1][Sp2][\[Theta]])}


DefTSIP[1,"pm",{Pm1_Symbol,Pp1_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[Lambda]_,BB_},{r_Symbol}]:={Pp1[r]->(-2 K[r]^2 Pm1[r]+(\[Lambda]+2 I r \[Omega]) Pm1[r] \[CapitalDelta][r]-2 I K[r] \[CapitalDelta][r] Derivative[1][Pm1][r])/(BB \[CapitalDelta][r]),Derivative[1][Pp1][r]->1/(BB \[CapitalDelta][r]^2) (2 I K[r]^3 Pm1[r]-2 I \[Lambda] K[r] Pm1[r] \[CapitalDelta][r]-2 K[r]^2 \[CapitalDelta][r] Derivative[1][Pm1][r]+\[CapitalDelta][r]^2 (2 I \[Omega] Pm1[r]+(\[Lambda]-2 I r \[Omega]) Derivative[1][Pm1][r]))}


DefTSIP[1,"mp",{Pm1_Symbol,Pp1_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[Lambda]_,BB_},{r_Symbol}]:={Pm1[r]->(-2 K[r]^2 Pp1[r]+(\[Lambda]-2 I r \[Omega]) Pp1[r] \[CapitalDelta][r]+2 I K[r] \[CapitalDelta][r] Derivative[1][Pp1][r])/(BB \[CapitalDelta][r]),Derivative[1][Pm1][r]->1/(BB \[CapitalDelta][r]^2) (-2 I K[r]^3 Pp1[r]+2 I \[Lambda] K[r] Pp1[r] \[CapitalDelta][r]-2 K[r]^2 \[CapitalDelta][r] Derivative[1][Pp1][r]+\[CapitalDelta][r]^2 (-2 I \[Omega] Pp1[r]+(\[Lambda]+2 I r \[Omega]) Derivative[1][Pp1][r]))}


DefTSIS[1,"pm",{Sm1_Symbol,Sp1_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[Lambda]_,BB_},{\[Theta]_Symbol}]:={Sp1[\[Theta]]->1/BB ((-\[Lambda]-4 a m \[Omega]-2 m Cot[\[Theta]] Csc[\[Theta]]+2 m^2 Csc[\[Theta]]^2+2 a^2 \[Omega]^2 Sin[\[Theta]]^2) Sm1[\[Theta]]+2 (a \[Omega]-m Csc[\[Theta]]^2) Sin[\[Theta]] Derivative[1][Sm1][\[Theta]]),Derivative[1][Sp1][\[Theta]]->1/BB ((2 (m \[Lambda]+a \[Omega]+3 a m^2 \[Omega]+m Cot[\[Theta]]^2) Csc[\[Theta]]-2 m^3 Csc[\[Theta]]^3+2 a \[Omega] Sin[\[Theta]] (-\[Lambda]-3 a m \[Omega]+a^2 \[Omega]^2 Sin[\[Theta]]^2)) Sm1[\[Theta]]+(-\[Lambda]-4 a m \[Omega]+2 m Cot[\[Theta]] Csc[\[Theta]]+2 m^2 Csc[\[Theta]]^2+2 a^2 \[Omega]^2 Sin[\[Theta]]^2) Derivative[1][Sm1][\[Theta]])}


DefTSIS[1,"mp",{Sm1_Symbol,Sp1_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[Lambda]_,BB_},{\[Theta]_Symbol}]:={Sm1[\[Theta]]->1/BB((-\[Lambda]-4 a m \[Omega]+2 m Cot[\[Theta]] Csc[\[Theta]]+2 m^2 Csc[\[Theta]]^2+2 a^2 \[Omega]^2 Sin[\[Theta]]^2) Sp1[\[Theta]]+2 (-a \[Omega]+m Csc[\[Theta]]^2) Sin[\[Theta]] Derivative[1][Sp1][\[Theta]]),Derivative[1][Sm1][\[Theta]]->1/BB((-2 (m \[Lambda]+a \[Omega]+3 a m^2 \[Omega]+m Cot[\[Theta]]^2) Csc[\[Theta]]+2 m^3 Csc[\[Theta]]^3+2 a \[Omega] Sin[\[Theta]] (\[Lambda]+3 a m \[Omega]-a^2 \[Omega]^2 Sin[\[Theta]]^2)) Sp1[\[Theta]]+(-\[Lambda]-4 a m \[Omega]-2 m Cot[\[Theta]] Csc[\[Theta]]+2 m^2 Csc[\[Theta]]^2+2 a^2 \[Omega]^2 Sin[\[Theta]]^2) Derivative[1][Sp1][\[Theta]])}


(* ::Subsection::Closed:: *)
(*Matrices*)


$Chop\[Epsilon]=10^-20


CleanValue[X_?NumericQ]:=Chop[X,$Chop\[Epsilon]];
CleanMatrix[M_]:=Map[CleanValue,M,{2}];


(* ::Input::Initialization:: *)
bmatEV[s_,m_,c_,lmax_,
opts:OptionsPattern[
{
WorkingPrecision:> 32,
AccuracyGoal->16,
"lpadding":> Automatic
}
]]:=
Module[
{
lpad,acc,lmin,btemp,evs
},
acc=OptionValue[AccuracyGoal];
lmin=Max[Abs[m]];
lpad=If[OptionValue["lpadding"]=== Automatic,
Max[If[PossibleZeroQ@c,0,Ceiling[(acc+4(1+Log10[Abs@c]))/(4-2Log10[Abs@c])]],5],
OptionValue["lpadding"]
];

{btemp,evs}=SWSHb[s,m,c,
WorkingPrecision->OptionValue[WorkingPrecision],
"lmax"->lmax+lpad,
"EigenValues"->True
];

{ArrayPad[
With[
{imax=Length@btemp},
Table[
Around[
btemp[[iI,iJ]],
If[iI+iJ<=imax,
10^-OptionValue[WorkingPrecision],
Max[
Abs@btemp[[iI+iJ-imax,imax]],
Abs@btemp[[imax,iJ+iI-imax]]
]
]
]["Value"],{iI,1,lmax+1-Abs[lmin]},{iJ,1,lmax+1-Abs[lmin]}]
]
,
{lmin,0}
],
ArrayPad[-s(s+1)-2m c +c^2-evs[[1;;(lmax+1-Max[Abs[m],Abs[s]])]],{Max[Abs[m],Abs[s]],0}]
}
]



(* ::Subsubsection::Closed:: *)
(*mixing matrices*)


mixsin[1,0,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m},
-Sqrt[2*(2*l2+1)/(2*l1+1)]*ClebschGordan[{l2,mm},{1,0},{l1,mm}]*ClebschGordan[{l2,-1},{1,1},{l1,0}]
]


mixsin[-1,0,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m},
Sqrt[2*(2*l2+1)/(2*l1+1)]*ClebschGordan[{l2,mm},{1,0},{l1,mm}]*ClebschGordan[{l2,1},{1,-1},{l1,0}]
]


mixsin[2,0,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m},
Sqrt[8*(2*l2+1)/(3*(2*l1+1))]*ClebschGordan[{l2,mm},{2,0},{l1,mm}]*ClebschGordan[{l2,-2},{2,2},{l1,0}]
]


mixsin[-2,0,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m},
Sqrt[8*(2*l2+1)/(3*(2*l1+1))]*ClebschGordan[{l2,mm},{2,0},{l1,mm}]*ClebschGordan[{l2,2},{2,-2},{l1,0}]
]


mixsin[2,1,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m},
-Sqrt[2*(2*l2+1)/((2*l1+1))]*ClebschGordan[{l2,mm},{1,0},{l1,mm}]*ClebschGordan[{l2,-2},{1,1},{l1,-1}]
]


mixsin[-2,-1,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m},
Sqrt[2*(2*l2+1)/((2*l1+1))]*ClebschGordan[{l2,mm},{1,0},{l1,mm}]*ClebschGordan[{l2,2},{1,-1},{l1,1}]
]


mixsin[-1,1,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m},
Sqrt[8*(2*l2+1)/(3*(2*l1+1))]*ClebschGordan[{l2,mm},{2,0},{l1,mm}]*ClebschGordan[{l2,1},{2,-2},{l1,-1}]
]


(* Use the symmetry. *)
mixsin[ 0, 1,l_,k_,m_]:=mixsin[ 1, 0,l+k,-k,m]
mixsin[ 0,-1,l_,k_,m_]:=mixsin[-1, 0,l+k,-k,m]
mixsin[ 0,-2,l_,k_,m_]:=mixsin[-2, 0,l+k,-k,m]
mixsin[ 0, 2,l_,k_,m_]:=mixsin[ 2, 0,l+k,-k,m]
mixsin[ 1,-1,l_,k_,m_]:=mixsin[-1, 1,l+k,-k,m]
mixsin[ 1, 2,l_,k_,m_]:=mixsin[ 2, 1,l+k,-k,m]
mixsin[-1,-2,l_,k_,m_]:=mixsin[-2,-1,l+k,-k,m]


mixsin[s1_,s2_,l_,k_,m_]:=0


mixcos[s_,1,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m,ss=s,t0=0},
Sqrt[(2*l2+1)/(2*l1+1)]*ClebschGordan[{l2,mm},{1,0},{l1,mm}]*ClebschGordan[{l2,-ss},{1,0},{l1,-ss}]
]


mixcos[s_,2,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m,ss=s,t0=0},
1/3*KroneckerDelta[l1,l2]+2/3*Sqrt[(2*l2+1)/(2*l1+1)]*ClebschGordan[{l2,mm},{2,0},{l1,mm}]*ClebschGordan[{l2,-ss},{2,0},{l1,-ss}]
]


mixcos[s_,n_,l_,k_,m_]:=0


ClearMixingMatrices[]:=(
Clear[SinMixMatrix,CosMixMatrix];

SinMixMatrix[s1_,s2_,lmax_,m_]:=(SinMixMatrix[s1,s2,lmax,m]=Module[{band,bands,k},
bands={};
For[k=-2,k<=2,k++,
band=Table[mixsin[s1,s2,ll,k,m],{ll,Max[-k,0],Min[lmax,lmax-k]}];
AppendTo[bands,band];
];
SparseArray[{Band[{1,1}]->bands[[3]],Band[{2,1}]->bands[[2]],Band[{3,1}]->bands[[1]],Band[{1,2}]->bands[[4]],Band[{1,3}]->bands[[5]]}, {lmax+1,lmax+1}]
]);

CosMixMatrix[s_,n_,lmax_,m_]:=(CosMixMatrix[s,n,lmax,m]=Module[{band,bands,k},
bands={};
For[k=-2,k<=2,k++,
band=Table[mixcos[s,n,ll,k,m],{ll,Max[-k,0],Min[lmax,lmax-k]}];
AppendTo[bands,band];
];
SparseArray[{Band[{1,1}]->bands[[3]],Band[{2,1}]->bands[[2]],Band[{3,1}]->bands[[1]],Band[{1,2}]->bands[[4]],Band[{1,3}]->bands[[5]]}, {lmax+1,lmax+1}]
]);
)


ClearMixingMatrices[]


(* ::Subsubsection::Closed:: *)
(*Diverse matrices*)


\[Lambda]hat[lmax_,m_]:=\[Lambda]hat[lmax,m]=DiagonalMatrix[Table[If[ll<Abs[m],0,Sqrt[ll*(ll+1)]],{ll,0,lmax}]];
\[Lambda]2hat[lmax_,m_]:=\[Lambda]2hat[lmax,m]=DiagonalMatrix[Table[If[ll<Max[Abs[m],1],0,Sqrt[(ll-1)*(ll+2)]],{ll,0,lmax}]];
\[CapitalLambda]hat[lmax_,m_]:=\[CapitalLambda]hat[lmax,m]=\[Lambda]hat[lmax,m] . \[Lambda]2hat[lmax,m];
signmat[lmax_,m_]:=signmat[lmax,m]=DiagonalMatrix[Table[If[ll<Abs[m],0,(-1)^(ll+m)],{ll,0,lmax}]];
idmat[lmax_,m_]:=idmat[lmax,m]=DiagonalMatrix[Table[If[ll<Abs[m],0,1],{ll,0,lmax}]];


(* ::Subsection::Closed:: *)
(*BuildGrid*)


BuildGrid[a_,{r0_,rmin_,rmax_},opts:OptionsPattern[{
WorkingPrecision->32,
"nterms"-> 8, (* Number of terms in the spherical expansion. *)
"kapord"->6, (* \[Kappa]ord is the maximum order of series expansion of \[Kappa] in spheroidal harmonics. *)
"rinf"->1000., (* rinf is the radius at which the series solutions for \[Kappa]_up should set the initial conditions for the integrator. *)
"rmax"->100.,(* rmax is the maximum value at which the interpolating function can be used, whereas rinf is the starting value for the numerical integration, i.e. the radius at which the initial condition for the UP solutions of kappa is set, using the series expansion. *)
"xhor"->0.001,
"inford"->5, (* The order of the expansion at infinity for the UP function. *)
"horord"->5, (* The order of the expansion at the horizon for the IN function. *)
AccuracyGoal->9,
"rgrid"->1, (* Use a linearly-spaced grid in the variable: 0 = rstar , 1 = r. *)
"rstmin"->Automatic,
"rstmax"->Automatic,
"nres"->4,(* Resolution in the r* direction:  dr* = M / n  (or dr = M / n). *)
"angres"->8 (* Resolution in the \[Theta] direction: number of points = nres * qres. *)
}]]:=Module[{rp,rm,prec,TortoiseCoord,InvTortoise,rmin2,rmax2,rgrid,rstarmin,rstarmax,qres,nres,qpts,dq,qs,r0star,drstar,dr,nleft,nright,rpts,rstarsL,rstarsR,rstars,rsL,rsR,rs},
rp=1+Sqrt[1-a^2];
rm=1+Sqrt[1-a^2];
prec=OptionValue[WorkingPrecision];
(*rmax=SetPrecision[OptionValue["rmax"],prec];*)   (* rmax is the maximum value at which the interpolating function can be used, whereas rinf is the starting value for the numerical integration, i.e. the radius at which the initial condition for the UP solutions of kappa is set, using the series expansion. *)
(*rmin =rhor + SetPrecision[10^-1,prec];*)  (* <--- may need to think again about this choice. *)
rgrid=SetPrecision[OptionValue["rgrid"],prec]; (* Use a linearly-spaced grid in the variable: 0 = rstar , 1 = r. *)
rstarmin=If[OptionValue["rstmin"]===Automatic,SetPrecision[rmin,prec],SetPrecision[OptionValue["rstmin"],prec]]; (* If rgrid=1 then these will be interpreted as rmin and rmax, instead of rstarmin and rstarmax. *)
rstarmax=If[OptionValue["rstmax"]===Automatic,SetPrecision[rmax,prec],SetPrecision[OptionValue["rstmax"],prec]];
nres=OptionValue["nres"];  (* Resolution in the r* direction:  dr* = M / n  (or dr = M / n). *)
qres=OptionValue["angres"]; (* Resolution in the \[Theta] direction: number of points = nres * qres. *)

TortoiseCoord[r_]:=
SetPrecision[r+(rp+rm)/(rp-rm)*(rp*Log[(r-rp)/2]-rm*Log[(r-rm)/2]),prec];
InvTortoise[x_?NumericQ]:=r/.FindRoot[TortoiseCoord[r]-x,{r,SetPrecision[rp+(rp (rm+rp) ProductLog[(2 (E^(((rm-rp) (rp-x))/(rm rp)))^(rm/(rm+rp)) (-rm+rp))/(rp (rm+rp))])/(-rm+rp),prec]},WorkingPrecision->prec];

If[rgrid==0,
rmin2=InvTortoise[rstarmin];
rmax2=InvTortoise[rstarmax];
If[rmin2<rmin,
Print["Warning: rstarmin is outside the allowed range."];
rstarmin=TortoiseCoord[rmin];
];
If[rmax2>rmax,
Print["Warning: rstarmax is outside the allowed range."];
rstarmax=TortoiseCoord[rmax];
];
,
rmin2=rstarmin;
rmax2=rstarmax;
If[rmin2<rmin,
Print["Warning: rstarmin is outside the allowed range."];
rstarmin=rmin;
];
If[rmax2>rmax,
Print["Warning: rstarmax is outside the allowed range."];
rstarmax=rmax;
];
];
qpts=qres*nres+1;
dq=\[Pi]/(qpts-1);
qs=Table[SetPrecision[(qi-1)*dq,prec],{qi,1,qpts}];
r0star=TortoiseCoord[r0];
(* Set up a grid in r* (if rgrid=0) or r (if rgrid=1). *)
If[rgrid==0,
drstar=1/nres;
nleft=Floor[(r0star-rstarmin)/drstar];
nright=Floor[(rstarmax-r0star)/drstar];
rpts=nleft+nright+1;
rstarsL=Table[SetPrecision[r0star-ri*drstar,prec],{ri,0,nleft}];
rstarsR=Table[SetPrecision[r0star+ri*drstar,prec],{ri,0,nright}];
rstars=Table[SetPrecision[r0star+(ri-1-nleft)*drstar,prec],{ri,1,rpts}];
rsL=Map[InvTortoise,rstarsL];
rsR=Map[InvTortoise,rstarsR];
rs=Map[InvTortoise,rstars];
,
(* If rgrid = 1 then we will construct a linear grid in r (rather than rstar). In this case, rstarmin and rstarmax (read from parameter files) are misnomers: they give the grid limits in r. *)
dr=1/nres;
nleft=Floor[(r0-rstarmin)/dr];
nright=Floor[(rstarmax-r0)/dr];
rpts=nleft+nright+1;
rsL=Table[SetPrecision[r0-ri*dr,prec],{ri,0,nleft}];
rsR=Table[SetPrecision[r0+ri*dr,prec],{ri,0,nright}];
rs=Table[SetPrecision[r0+(ri-1-nleft)*dr,prec],{ri,1,rpts}];
rstarsL=Map[TortoiseCoord,rsL];
rstarsR=Map[TortoiseCoord,rsR];
rstars=Map[TortoiseCoord,rs];
];

(*filen=directory<>"data/lm_rs_"<>ToString[iConfig]<>".dat";
Export[filen,Transpose@{rstars,rs}];
filen=directory<>"data/lm_qs_"<>ToString[iConfig]<>".dat";
Export[filen,qs];
filen=directory<>"data/lm_rsL_"<>ToString[iConfig]<>".dat";
Export[filen,Transpose@{rstarsL,rsL}];
filen=directory<>"data/lm_rsR_"<>ToString[iConfig]<>".dat";
Export[filen,Transpose@{rstarsR,rsR}];*)
Return[{rsL,rsR}]
]


(* ::Subsection::Closed:: *)
(*Spin 2 calculation*)


(* ::Subsubsection::Closed:: *)
(*Spin 2 derivatives*)


DdagPp2[a_,m_,\[Omega]_,r_]:=Function[{Pp2,Pp21,Pp22,Pp23,Pp24},Evaluate[
Pp21+(I Pp2 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)
]
]


DdagDdagPp2[a_,m_,\[Omega]_,r_]:=Function[{Pp2,Pp21,Pp22,Pp23,Pp24},Evaluate[
Pp22+(2 I Pp21 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)-(Pp2 (2 I a m+a^2 m^2-2 I a m r-2 I a^2 \[Omega]-2 a^3 m \[Omega]+2 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


DdagDdagDdagPp2[a_,m_,\[Omega]_,r_]:=Function[{Pp2,Pp21,Pp22,Pp23,Pp24},Evaluate[
Pp23+(3 I Pp22 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)-(3 Pp21 (2 I a m+a^2 m^2-2 I a m r-2 I a^2 \[Omega]-2 a^3 m \[Omega]+2 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2+1/(a^2-2 r+r^2)^3 Pp2 (-8 I a m+2 I a^3 m-6 a^2 m^2+I a^3 m^3+12 I a m r+6 a^2 m^2 r-6 I a m r^2+8 I a^2 \[Omega]+12 a^3 m \[Omega]-3 I a^4 m^2 \[Omega]-12 I a^2 r \[Omega]-6 a^3 m r \[Omega]-3 I a^2 m^2 r^2 \[Omega]+4 I r^3 \[Omega]-6 a m r^3 \[Omega]-6 a^4 \[Omega]^2+3 I a^5 m \[Omega]^2+6 I a^3 m r^2 \[Omega]^2+6 r^4 \[Omega]^2+3 I a m r^4 \[Omega]^2-I a^6 \[Omega]^3-3 I a^4 r^2 \[Omega]^3-3 I a^2 r^4 \[Omega]^3-I r^6 \[Omega]^3)
]
]


DPm2[a_,m_,\[Omega]_,r_]:=Function[{Pm2,Pm21,Pm22,Pm23,Pm24},Evaluate[
Pm21-(I Pm2 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)
]
]


DDPm2[a_,m_,\[Omega]_,r_]:=Function[{Pm2,Pm21,Pm22,Pm23,Pm24},Evaluate[
Pm22-(2 I Pm21 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)-(Pm2 (-2 I a m+a^2 m^2+2 I a m r+2 I a^2 \[Omega]-2 a^3 m \[Omega]-2 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


DDDPm2[a_,m_,\[Omega]_,r_]:=Function[{Pm2,Pm21,Pm22,Pm23,Pm24},Evaluate[
Pm23-(3 I Pm22 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)-(3 Pm21 (-2 I a m+a^2 m^2+2 I a m r+2 I a^2 \[Omega]-2 a^3 m \[Omega]-2 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2+1/(a^2-2 r+r^2)^3 Pm2 (8 I a m-2 I a^3 m-6 a^2 m^2-I a^3 m^3-12 I a m r+6 a^2 m^2 r+6 I a m r^2-8 I a^2 \[Omega]+12 a^3 m \[Omega]+3 I a^4 m^2 \[Omega]+12 I a^2 r \[Omega]-6 a^3 m r \[Omega]+3 I a^2 m^2 r^2 \[Omega]-4 I r^3 \[Omega]-6 a m r^3 \[Omega]-6 a^4 \[Omega]^2-3 I a^5 m \[Omega]^2-6 I a^3 m r^2 \[Omega]^2+6 r^4 \[Omega]^2-3 I a m r^4 \[Omega]^2+I a^6 \[Omega]^3+3 I a^4 r^2 \[Omega]^3+3 I a^2 r^4 \[Omega]^3+I r^6 \[Omega]^3)
]
]


DopDdagDdagPp2[a_,m_,\[Omega]_,r_]:=Function[{Pp2,Pp21,Pp22,Pp23,Pp24},Evaluate[
Pp23+(I Pp22 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)+(Pp21 (-6 I a m+a^2 m^2+6 I a m r+6 I a^2 \[Omega]-2 a^3 m \[Omega]-6 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2+1/(a^2-2 r+r^2)^3 Pp2 (-8 I a m+2 I a^3 m-2 a^2 m^2-I a^3 m^3+12 I a m r+2 a^2 m^2 r-6 I a m r^2+8 I a^2 \[Omega]+4 a^3 m \[Omega]+3 I a^4 m^2 \[Omega]-12 I a^2 r \[Omega]-2 a^3 m r \[Omega]+3 I a^2 m^2 r^2 \[Omega]+4 I r^3 \[Omega]-2 a m r^3 \[Omega]-2 a^4 \[Omega]^2-3 I a^5 m \[Omega]^2-6 I a^3 m r^2 \[Omega]^2+2 r^4 \[Omega]^2-3 I a m r^4 \[Omega]^2+I a^6 \[Omega]^3+3 I a^4 r^2 \[Omega]^3+3 I a^2 r^4 \[Omega]^3+I r^6 \[Omega]^3)
]
]


DdagDDPm2[a_,m_,\[Omega]_,r_]:=Function[{Pm2,Pm21,Pm22,Pm23,Pm24},Evaluate[
Pm23-(I Pm22 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)+(Pm21 (6 I a m+a^2 m^2-6 I a m r-6 I a^2 \[Omega]-2 a^3 m \[Omega]+6 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2+1/(a^2-2 r+r^2)^3 Pm2 (8 I a m-2 I a^3 m-2 a^2 m^2+I a^3 m^3-12 I a m r+2 a^2 m^2 r+6 I a m r^2-8 I a^2 \[Omega]+4 a^3 m \[Omega]-3 I a^4 m^2 \[Omega]+12 I a^2 r \[Omega]-2 a^3 m r \[Omega]-3 I a^2 m^2 r^2 \[Omega]-4 I r^3 \[Omega]-2 a m r^3 \[Omega]-2 a^4 \[Omega]^2+3 I a^5 m \[Omega]^2+6 I a^3 m r^2 \[Omega]^2+2 r^4 \[Omega]^2+3 I a m r^4 \[Omega]^2-I a^6 \[Omega]^3-3 I a^4 r^2 \[Omega]^3-3 I a^2 r^4 \[Omega]^3-I r^6 \[Omega]^3)
]
]


(* ::Subsubsection::Closed:: *)
(*Spin 2 MRProjectors*)


(* ::Input::Initialization:: *)
MRProjector[2,"l+l+"][a_,m_,\[Omega]_,lmax_,r_,{Pp2vec_,Pm2vec_,DdagPp2vec_,
DdagDdagPp2vec_,
DdagDdagDdagPp2vec_,
DPm2vec_,
DDPm2vec_,
DDDPm2vec_,
DopDdagDdagPp2vec_,
DdagDDPm2vec_},{ Bmatm2_,Bmatp2_}]:=With[
{
Ss2lplp=((Bmatp2 . (\[CapitalLambda]hat[lmax,m]-2*a \[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[1,0,lmax,m]+(a \[Omega])^2*SinMixMatrix[2,0,lmax,m]))+(signmat[lmax,m] . Bmatm2 . (\[CapitalLambda]hat[lmax,m]-2*a \[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[-1,0,lmax,m]+(a \[Omega])^2*SinMixMatrix[-2,0,lmax,m])))
},
(-1/(6*\[Omega]^2*(a^2-2 r+r^2)^2))*(Pp2vec . Ss2lplp)
]


(* ::Input::Initialization:: *)
MRProjector[2,"l-l-"][a_,m_,\[Omega]_,lmax_,r_,{Pp2vec_,Pm2vec_,DdagPp2vec_,
DdagDdagPp2vec_,
DdagDdagDdagPp2vec_,
DPm2vec_,
DDPm2vec_,
DDDPm2vec_,
DopDdagDdagPp2vec_,
DdagDDPm2vec_},{ Bmatm2_,Bmatp2_}]:=With[
{
Ss2lmlm=((signmat[lmax,m] . Bmatp2 . (\[CapitalLambda]hat[lmax,m]-2*a \[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[1,0,lmax,m]+(a \[Omega])^2*SinMixMatrix[2,0,lmax,m]))+( Bmatm2 . (\[CapitalLambda]hat[lmax,m]-2*a \[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[-1,0,lmax,m]+(a \[Omega])^2*SinMixMatrix[-2,0,lmax,m])))
},
-1/(6*\[Omega]^2*(a^2-2 r+r^2)^2)*Pm2vec . Ss2lmlm
]


(* ::Input::Initialization:: *)
MRProjector[2,"m+m+"][a_,m_,\[Omega]_,lmax_,r_,{Pp2vec_,Pm2vec_,DdagPp2vec_,
DdagDdagPp2vec_,
DdagDdagDdagPp2vec_,
DPm2vec_,
DDPm2vec_,
DDDPm2vec_,
DopDdagDdagPp2vec_,
DdagDDPm2vec_},{ Bmatm2_,Bmatp2_}]:=With[{
tmp=-1/(6*\[Omega]^2)*(DdagDdagPp2vec+signmat[lmax,m] . DDPm2vec),
Ss2mpmp=Bmatp2
},
tmp . Ss2mpmp
]


(* ::Input::Initialization:: *)
MRProjector[2,"m-m-"][a_,m_,\[Omega]_,lmax_,r_,{Pp2vec_,Pm2vec_,DdagPp2vec_,
DdagDdagPp2vec_,
DdagDdagDdagPp2vec_,
DPm2vec_,
DDPm2vec_,
DDDPm2vec_,
DopDdagDdagPp2vec_,
DdagDDPm2vec_},{ Bmatm2_,Bmatp2_}]:=With[{
tmp=-1/(6*\[Omega]^2)*(DdagDdagPp2vec+signmat[lmax,m] . DDPm2vec),
Ss2mmmm=signmat[lmax,m] . Bmatm2
},
tmp . Ss2mmmm
]


(* ::Input::Initialization:: *)
MRProjector[2,"l+m+"][a_,m_,\[Omega]_,lmax_,r_,{Pp2vec_,Pm2vec_,DdagPp2vec_,
DdagDdagPp2vec_,
DdagDdagDdagPp2vec_,
DPm2vec_,
DDPm2vec_,
DDDPm2vec_,
DopDdagDdagPp2vec_,
DdagDDPm2vec_},{ Bmatm2_,Bmatp2_}]:=Module[{
term1,term2,term3,term4,tmp,\[Rho]mat,SSplus,SSminus,MMp2,MMm2,sinmix,sinmix1m1
},
SSplus=(Bmatp2 . (\[CapitalLambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[1,0,lmax,m]+a^2*\[Omega]^2*SinMixMatrix[2,0,lmax,m]));
SSminus=(Bmatm2 . (\[CapitalLambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[-1,0,lmax,m]+a^2*\[Omega]^2*SinMixMatrix[-2,0,lmax,m]));
MMp2=Bmatp2 . (\[Lambda]2hat[lmax,m]-a*\[Omega]*SinMixMatrix[2,1,lmax,m]);
MMm2=Bmatm2 . (\[Lambda]2hat[lmax,m]-a*\[Omega]*SinMixMatrix[-2,-1,lmax,m]);

\[Rho]mat=r*idmat[lmax,m]+I*a*CosMixMatrix[1,1,lmax,m];

sinmix=SinMixMatrix[0,1,lmax,m];
sinmix1m1=SinMixMatrix[-1,1,lmax,m];

term1=1/(12*\[Omega]^2)*(DDPm2vec . signmat[lmax,m] . SSplus-DdagDdagPp2vec . signmat[lmax,m] . SSminus) . (\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix);
term2=-1/(12*\[Omega]^2)*(DDDPm2vec . signmat[lmax,m] . SSplus-DopDdagDdagPp2vec . signmat[lmax,m] . SSminus) . ((\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix) . \[Rho]mat-I*a*sinmix);
term3=-1/(3*I*\[Omega])*(DDDPm2vec . signmat[lmax,m] . MMp2 . \[Rho]mat . \[Rho]mat-2*DDPm2vec . signmat [lmax,m] . MMp2 . \[Rho]mat+2*DPm2vec . signmat[lmax,m] . MMp2);
tmp=(\[Lambda]hat[lmax,m] . \[Lambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]hat[lmax,m] . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]mat . \[Rho]mat-2*I*a*(\[Lambda]hat[lmax,m] . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]mat-2*a^2*sinmix1m1;
term4=-1/(3*I*\[Omega]*(a^2-2 r+r^2))*DdagPp2vec . (signmat[lmax,m] . MMm2 . tmp);
tmp=term1+term2+term3+term4;
Return[tmp/(-6*I*\[Omega])]
]


(* ::Input::Initialization:: *)
MRProjector[2,"l+m-"][a_,m_,\[Omega]_,lmax_,r_,{Pp2vec_,Pm2vec_,DdagPp2vec_,
DdagDdagPp2vec_,
DdagDdagDdagPp2vec_,
DPm2vec_,
DDPm2vec_,
DDDPm2vec_,
DopDdagDdagPp2vec_,
DdagDDPm2vec_},{ Bmatm2_,Bmatp2_}]:=Module[{
term1,term2,term3,term4,tmp,\[Rho]cmat,SSplus,SSminus,MMp2,MMm2,sinmix,sinmix1m1
},
SSplus=(Bmatp2 . (\[CapitalLambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[1,0,lmax,m]+a^2*\[Omega]^2*SinMixMatrix[2,0,lmax,m]));SSminus=(Bmatm2 . (\[CapitalLambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[-1,0,lmax,m]+a^2*\[Omega]^2*SinMixMatrix[-2,0,lmax,m]));
\[Rho]cmat=r*idmat[lmax,m]-I*a*CosMixMatrix[-1,1,lmax,m];
MMp2=Bmatp2 . (\[Lambda]2hat[lmax,m]-a*\[Omega]*SinMixMatrix[2,1,lmax,m]);
MMm2=Bmatm2 . (\[Lambda]2hat[lmax,m]-a*\[Omega]*SinMixMatrix[-2,-1,lmax,m]);
sinmix=SinMixMatrix[0,-1,lmax,m];
sinmix1m1=SinMixMatrix[1,-1,lmax,m];

term1=-1/(12*\[Omega]^2)*(DDPm2vec . SSminus-DdagDdagPp2vec . SSplus) . (\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix);term2=1/(12*\[Omega]^2)*(DDDPm2vec . SSminus-DopDdagDdagPp2vec . SSplus) . ((\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix) . \[Rho]cmat-I*a*sinmix);term3=1/(3*I*\[Omega])*(DDDPm2vec . MMm2 . \[Rho]cmat . \[Rho]cmat-2*DDPm2vec . MMm2 . \[Rho]cmat+2*DPm2vec . MMm2);tmp=((\[Lambda]hat[lmax,m] . \[Lambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]hat[lmax,m] . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]cmat . \[Rho]cmat-2*I*a*(\[Lambda]hat[lmax,m] . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]cmat-2*a^2*sinmix1m1);term4=1/(3*I*\[Omega]*(a^2-2 r+r^2))*DdagPp2vec . MMp2 . tmp;
tmp=term1+term2+term3+term4;
Return[tmp/(-6*I*\[Omega])]
]


(* ::Input::Initialization:: *)
MRProjector[2,"l-m+"][a_,m_,\[Omega]_,lmax_,r_,{Pp2vec_,Pm2vec_,DdagPp2vec_,
DdagDdagPp2vec_,
DdagDdagDdagPp2vec_,
DPm2vec_,
DDPm2vec_,
DDDPm2vec_,
DopDdagDdagPp2vec_,
DdagDDPm2vec_},{ Bmatm2_,Bmatp2_}]:=Module[{
term1,term2,term3,term4,tmp,\[Rho]cmat,SSplus,SSminus,MMp2,MMm2,sinmix,sinmix1m1
},
SSplus=(Bmatp2 . (\[CapitalLambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[1,0,lmax,m]+a^2*\[Omega]^2*SinMixMatrix[2,0,lmax,m]));SSminus=(Bmatm2 . (\[CapitalLambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[-1,0,lmax,m]+a^2*\[Omega]^2*SinMixMatrix[-2,0,lmax,m]));
MMp2=Bmatp2 . (\[Lambda]2hat[lmax,m]-a*\[Omega]*SinMixMatrix[2,1,lmax,m]);
MMm2=Bmatm2 . (\[Lambda]2hat[lmax,m]-a*\[Omega]*SinMixMatrix[-2,-1,lmax,m]);

sinmix=SinMixMatrix[0,1,lmax,m];
sinmix1m1=SinMixMatrix[-1,1,lmax,m];
\[Rho]cmat=r*idmat[lmax,m]-I*a*CosMixMatrix[1,1,lmax,m];

term1=1/(12*\[Omega]^2)*(DDPm2vec . SSminus-DdagDdagPp2vec . SSplus) . (\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix);
term2=-1/(12*\[Omega]^2)*(DdagDDPm2vec . SSminus-DdagDdagDdagPp2vec . SSplus) . ((\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix) . \[Rho]cmat+I*a*sinmix);
term3=-1/(3*I*\[Omega])*(DdagDdagDdagPp2vec . MMp2 . \[Rho]cmat . \[Rho]cmat-2*DdagDdagPp2vec . MMp2 . \[Rho]cmat+2*DdagPp2vec . MMp2);
tmp=((\[Lambda]hat [lmax,m] . \[Lambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]hat[lmax,m] . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]cmat . \[Rho]cmat+2*I*a*(\[Lambda]hat[lmax,m] . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]cmat-2*a^2*sinmix1m1);
term4=-1/(3*I*\[Omega]*(a^2-2 r+r^2))*DPm2vec . MMm2 . tmp;
tmp=term1+term2+term3+term4;
Return[tmp/(-6*I*\[Omega])]
]


(* ::Input::Initialization:: *)
MRProjector[2,"l-m-"][a_,m_,\[Omega]_,lmax_,r_,{Pp2vec_,Pm2vec_,DdagPp2vec_,
DdagDdagPp2vec_,
DdagDdagDdagPp2vec_,
DPm2vec_,
DDPm2vec_,
DDDPm2vec_,
DopDdagDdagPp2vec_,
DdagDDPm2vec_},{ Bmatm2_,Bmatp2_}]:=Module[{
term1,term2,term3,term4,tmp,\[Rho]mat,SSplus,SSminus,MMp2,MMm2,sinmix,sinmix1m1
},
SSplus=(Bmatp2 . (\[CapitalLambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[1,0,lmax,m]+a^2*\[Omega]^2*SinMixMatrix[2,0,lmax,m]));SSminus=(Bmatm2 . (\[CapitalLambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[-1,0,lmax,m]+a^2*\[Omega]^2*SinMixMatrix[-2,0,lmax,m]));
MMp2=Bmatp2 . (\[Lambda]2hat[lmax,m]-a*\[Omega]*SinMixMatrix[2,1,lmax,m]);
MMm2=Bmatm2 . (\[Lambda]2hat[lmax,m]-a*\[Omega]*SinMixMatrix[-2,-1,lmax,m]);

sinmix=SinMixMatrix[0,-1,lmax,m];
sinmix1m1=SinMixMatrix[1,-1,lmax,m];
\[Rho]mat=r*idmat[lmax,m]+I*a*CosMixMatrix[-1,1,lmax,m];

term1=-1/(12*\[Omega]^2)*(DDPm2vec . signmat[lmax,m] . SSplus-DdagDdagPp2vec . signmat[lmax,m] . SSminus) . (\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix);term2=1/(12*\[Omega]^2)*(DdagDDPm2vec . signmat[lmax,m] . SSplus-DdagDdagDdagPp2vec . signmat [lmax,m] . SSminus) . ((\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix) . \[Rho]mat+I*a*sinmix);term3=1/(3*I*\[Omega])*(DdagDdagDdagPp2vec . signmat[lmax,m] . MMm2 . \[Rho]mat . \[Rho]mat-2*DdagDdagPp2vec . signmat[lmax,m] . MMm2 . \[Rho]mat+2*DdagPp2vec . signmat[lmax,m] . MMm2);tmp=((\[Lambda]hat [lmax,m] . \[Lambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]hat[lmax,m] . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]mat . \[Rho]mat+2*I*a*(\[Lambda]hat[lmax,m] . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]mat-2*a^2*sinmix1m1);term4=1/(3*I*\[Omega]*(a^2-2 r+r^2))*DPm2vec . (signmat[lmax,m] . MMp2 . tmp);tmp=term1+term2+term3+term4;Return[tmp/(-6*I*\[Omega])]
]


(* ::Input::Initialization:: *)
MRProjector[2,"l+l-"][a_,m_,\[Omega]_,lmax_,r_,{Pp2vec_,Pm2vec_,DdagPp2vec_,
DdagDdagPp2vec_,
DdagDdagDdagPp2vec_,
DPm2vec_,
DDPm2vec_,
DDDPm2vec_,
DopDdagDdagPp2vec_,
DdagDDPm2vec_},{Pp2vec2_,Pm2vec2_},{ Bmatm2_,Bmatp2_}]:=Module[{
term1,term2,term3,term4,tmp,tmp3,\[Rho]mat,\[Rho]cmat,SSplus,SSminus,MMp2,MMm2,sinmix,sinmix1m1,cosmix01,cosmix02,\[CapitalSigma]mat,sinp10,sinm10,
lplmterm1A ,lplmterm1B ,lplmterm2A ,lplmterm2B 
},
SSplus=(Bmatp2 . (\[CapitalLambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[1,0,lmax,m]+a^2*\[Omega]^2*SinMixMatrix[2,0,lmax,m]));
SSminus=(Bmatm2 . (\[CapitalLambda]hat[lmax,m]-2*a*\[Omega]*\[Lambda]2hat[lmax,m] . SinMixMatrix[-1,0,lmax,m]+a^2*\[Omega]^2*SinMixMatrix[-2,0,lmax,m]));
MMp2=Bmatp2 . (\[Lambda]2hat[lmax,m]-a*\[Omega]*SinMixMatrix[2,1,lmax,m]);
MMm2=Bmatm2 . (\[Lambda]2hat[lmax,m]-a*\[Omega]*SinMixMatrix[-2,-1,lmax,m]);

cosmix01=CosMixMatrix[0,1,lmax,m];
cosmix02=CosMixMatrix[0,2,lmax,m];
\[Rho]mat=r*idmat[lmax,m]+I*a*cosmix01;
\[Rho]cmat=r*idmat[lmax,m]-I*a*cosmix01;
\[CapitalSigma]mat=r^2*idmat[lmax,m]+a^2*cosmix02;
sinp10=SinMixMatrix[1,0,lmax,m];
sinm10=SinMixMatrix[-1,0,lmax,m];


lplmterm1A =alplmterm1A[a,m,\[Omega],r]@@@Pm2vec2;
lplmterm1B =alplmterm1B[a,m,\[Omega],r]@@@Pm2vec2;
lplmterm2A=alplmterm2A[a,m,\[Omega],r]@@@Pp2vec2;
lplmterm2B=alplmterm2B[a,m,\[Omega],r]@@@Pp2vec2;


tmp3=\[Lambda]hat[lmax,m] . (sinm10-sinp10);
term1=1/(24*\[Omega]^2)*(lplmterm1A . SSminus . \[CapitalSigma]mat-4*r*(a^2-2 r+r^2)*lplmterm1B . SSminus-2*a^2*DDPm2vec . SSminus . tmp3 . cosmix01);
term2=-1/(24*\[Omega]^2)*(lplmterm2A . SSplus . \[CapitalSigma]mat-4*r*(a^2-2 r+r^2)*lplmterm2B . SSplus-2*a^2*DdagDdagPp2vec . SSplus . tmp3 . cosmix01);
tmp=(2*a^2*r*cosmix01-I*a*(r^2*idmat[lmax,m]+3*a^2*cosmix02));
term3=1/(3*I*\[Omega])*(DDPm2vec . SSminus . \[CapitalSigma]mat . \[Rho]cmat-DPm2vec . SSminus . \[Rho]cmat . \[Rho]cmat
-DDPm2vec . MMm2 . sinm10 . tmp
-2*I*a*DPm2vec . MMm2 . sinm10 . \[Rho]mat);
term4=1/(3*I*\[Omega])*(DdagDdagPp2vec . SSplus . \[CapitalSigma]mat . \[Rho]cmat-DdagPp2vec . SSplus . \[Rho]cmat . \[Rho]cmat
+DdagDdagPp2vec . MMp2 . sinp10 . tmp
+2*I*a*DdagPp2vec . MMp2 . sinp10 . \[Rho]mat);
tmp=term1+term2+term3+term4;Return[tmp/(-6*I*\[Omega])]
]


(* ::Input::Initialization:: *)
alplmterm1A[a_,m_,\[Omega]_,r_]:=Function[{Pm2,Pm21,Pm22,Pm23,Pm24},Evaluate[
2 Pm24 (a^2-2 r+r^2)+4 Pm23 (-1+I a m+r-I a^2 \[Omega]-I r^2 \[Omega])-(4 I Pm22 (-3 a m+3 a m r+3 a^2 \[Omega]+2 a^2 r \[Omega]-7 r^2 \[Omega]+2 r^3 \[Omega]))/(a^2-2 r+r^2)+1/(a^2-2 r+r^2)^2 4 Pm21 (10 I a m-4 I a^3 m-3 a^2 m^2+I a^3 m^3-12 I a m r+3 a^2 m^2 r+6 I a m r^2-10 I a^2 \[Omega]+6 a^3 m \[Omega]-3 I a^4 m^2 \[Omega]+18 I a^2 r \[Omega]-2 a^3 m r \[Omega]-6 I r^2 \[Omega]-2 a m r^2 \[Omega]-3 I a^2 m^2 r^2 \[Omega]-2 I r^3 \[Omega]-2 a m r^3 \[Omega]-3 a^4 \[Omega]^2+3 I a^5 m \[Omega]^2-a^4 r \[Omega]^2+2 a^2 r^2 \[Omega]^2+6 I a^3 m r^2 \[Omega]^2-2 a^2 r^3 \[Omega]^2+5 r^4 \[Omega]^2+3 I a m r^4 \[Omega]^2-r^5 \[Omega]^2-I a^6 \[Omega]^3-3 I a^4 r^2 \[Omega]^3-3 I a^2 r^4 \[Omega]^3-I r^6 \[Omega]^3)-1/(a^2-2 r+r^2)^3 2 Pm2 (-32 I a m+20 I a^3 m+16 a^2 m^2-4 a^4 m^2-2 I a^3 m^3+a^4 m^4+56 I a m r-20 I a^3 m r-24 a^2 m^2 r+2 I a^3 m^3 r-36 I a m r^2+12 a^2 m^2 r^2+12 I a m r^3+32 I a^2 \[Omega]-12 I a^4 \[Omega]-32 a^3 m \[Omega]+4 a^5 m \[Omega]+6 I a^4 m^2 \[Omega]-4 a^5 m^3 \[Omega]-56 I a^2 r \[Omega]+40 a^3 m r \[Omega]-4 I a^4 m^2 r \[Omega]+48 I a^2 r^2 \[Omega]+2 I a^2 m^2 r^2 \[Omega]-4 a^3 m^3 r^2 \[Omega]-8 I r^3 \[Omega]-8 a m r^3 \[Omega]-4 I a^2 m^2 r^3 \[Omega]-4 I r^4 \[Omega]-4 a m r^4 \[Omega]+16 a^4 \[Omega]^2-6 I a^5 m \[Omega]^2+6 a^6 m^2 \[Omega]^2-16 a^4 r \[Omega]^2+2 I a^5 m r \[Omega]^2-4 I a^3 m r^2 \[Omega]^2+12 a^4 m^2 r^2 \[Omega]^2-16 a^2 r^3 \[Omega]^2+4 I a^3 m r^3 \[Omega]^2+16 r^4 \[Omega]^2+2 I a m r^4 \[Omega]^2+6 a^2 m^2 r^4 \[Omega]^2+2 I a m r^5 \[Omega]^2+2 I a^6 \[Omega]^3-4 a^7 m \[Omega]^3+2 I a^4 r^2 \[Omega]^3-12 a^5 m r^2 \[Omega]^3-2 I a^2 r^4 \[Omega]^3-12 a^3 m r^4 \[Omega]^3-2 I r^6 \[Omega]^3-4 a m r^6 \[Omega]^3+a^8 \[Omega]^4+4 a^6 r^2 \[Omega]^4+6 a^4 r^4 \[Omega]^4+4 a^2 r^6 \[Omega]^4+r^8 \[Omega]^4)
]
]


(* ::Input::Initialization:: *)
alplmterm1B[a_,m_,\[Omega]_,r_]:=Function[{Pm2,Pm21,Pm22,Pm23,Pm24},Evaluate[
Pm23-(2 I Pm22 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)-1/(a^2-2 r+r^2)^3 2 Pm2 (-4 I a m+I a^3 m+2 a^2 m^2+6 I a m r-2 a^2 m^2 r-3 I a m r^2+4 I a^2 \[Omega]-4 a^3 m \[Omega]-6 I a^2 r \[Omega]+2 a^3 m r \[Omega]+2 I r^3 \[Omega]+2 a m r^3 \[Omega]+2 a^4 \[Omega]^2-2 r^4 \[Omega]^2)-(Pm21 (-6 I a m+a^2 m^2+6 I a m r+6 I a^2 \[Omega]-2 a^3 m \[Omega]-6 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


(* ::Input::Initialization:: *)
alplmterm2A[a_,m_,\[Omega]_,r_]:=Function[{Pp2,Pp21,Pp22,Pp23,Pp24},Evaluate[
2 Pp24 (a^2-2 r+r^2)+4 Pp23 (-1-I a m+r+I a^2 \[Omega]+I r^2 \[Omega])+(4 I Pp22 (-3 a m+3 a m r+3 a^2 \[Omega]+2 a^2 r \[Omega]-7 r^2 \[Omega]+2 r^3 \[Omega]))/(a^2-2 r+r^2)+1/(a^2-2 r+r^2)^2 4 Pp21 (-10 I a m+4 I a^3 m-3 a^2 m^2-I a^3 m^3+12 I a m r+3 a^2 m^2 r-6 I a m r^2+10 I a^2 \[Omega]+6 a^3 m \[Omega]+3 I a^4 m^2 \[Omega]-18 I a^2 r \[Omega]-2 a^3 m r \[Omega]+6 I r^2 \[Omega]-2 a m r^2 \[Omega]+3 I a^2 m^2 r^2 \[Omega]+2 I r^3 \[Omega]-2 a m r^3 \[Omega]-3 a^4 \[Omega]^2-3 I a^5 m \[Omega]^2-a^4 r \[Omega]^2+2 a^2 r^2 \[Omega]^2-6 I a^3 m r^2 \[Omega]^2-2 a^2 r^3 \[Omega]^2+5 r^4 \[Omega]^2-3 I a m r^4 \[Omega]^2-r^5 \[Omega]^2+I a^6 \[Omega]^3+3 I a^4 r^2 \[Omega]^3+3 I a^2 r^4 \[Omega]^3+I r^6 \[Omega]^3)-1/(a^2-2 r+r^2)^3 2 Pp2 (32 I a m-20 I a^3 m+16 a^2 m^2-4 a^4 m^2+2 I a^3 m^3+a^4 m^4-56 I a m r+20 I a^3 m r-24 a^2 m^2 r-2 I a^3 m^3 r+36 I a m r^2+12 a^2 m^2 r^2-12 I a m r^3-32 I a^2 \[Omega]+12 I a^4 \[Omega]-32 a^3 m \[Omega]+4 a^5 m \[Omega]-6 I a^4 m^2 \[Omega]-4 a^5 m^3 \[Omega]+56 I a^2 r \[Omega]+40 a^3 m r \[Omega]+4 I a^4 m^2 r \[Omega]-48 I a^2 r^2 \[Omega]-2 I a^2 m^2 r^2 \[Omega]-4 a^3 m^3 r^2 \[Omega]+8 I r^3 \[Omega]-8 a m r^3 \[Omega]+4 I a^2 m^2 r^3 \[Omega]+4 I r^4 \[Omega]-4 a m r^4 \[Omega]+16 a^4 \[Omega]^2+6 I a^5 m \[Omega]^2+6 a^6 m^2 \[Omega]^2-16 a^4 r \[Omega]^2-2 I a^5 m r \[Omega]^2+4 I a^3 m r^2 \[Omega]^2+12 a^4 m^2 r^2 \[Omega]^2-16 a^2 r^3 \[Omega]^2-4 I a^3 m r^3 \[Omega]^2+16 r^4 \[Omega]^2-2 I a m r^4 \[Omega]^2+6 a^2 m^2 r^4 \[Omega]^2-2 I a m r^5 \[Omega]^2-2 I a^6 \[Omega]^3-4 a^7 m \[Omega]^3-2 I a^4 r^2 \[Omega]^3-12 a^5 m r^2 \[Omega]^3+2 I a^2 r^4 \[Omega]^3-12 a^3 m r^4 \[Omega]^3+2 I r^6 \[Omega]^3-4 a m r^6 \[Omega]^3+a^8 \[Omega]^4+4 a^6 r^2 \[Omega]^4+6 a^4 r^4 \[Omega]^4+4 a^2 r^6 \[Omega]^4+r^8 \[Omega]^4)
]
]


(* ::Input::Initialization:: *)
alplmterm2B[a_,m_,\[Omega]_,r_]:=Function[{Pp2,Pp21,Pp22,Pp23,Pp24},Evaluate[
Pp23+Pp22((2 I  (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2))-Pp21 ((6 I a m+a^2 m^2-6 I a m r-6 I a^2 \[Omega]-2 a^3 m \[Omega]+6 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2)/(a^2-2 r+r^2)^2)-Pp2 (1/(a^2-2 r+r^2)^3 2 (4 I a m-I a^3 m+2 a^2 m^2-6 I a m r-2 a^2 m^2 r+3 I a m r^2-4 I a^2 \[Omega]-4 a^3 m \[Omega]+6 I a^2 r \[Omega]+2 a^3 m r \[Omega]-2 I r^3 \[Omega]+2 a m r^3 \[Omega]+2 a^4 \[Omega]^2-2 r^4 \[Omega]^2))
]
]


(* ::Subsubsection::Closed:: *)
(*CalculateSpin2Contributions*)


CalculateSpin2Contributions[a_,m_,\[Omega]_,lmax_,r_,{Pp2vec_,Pm2vec_},{Bmatm2_,Bmatp2_}]:=Module[{
Pvecs},
Pvecs={
Pp2vec[[All,1]],
Pm2vec[[All,1]],
DdagPp2[a,m,\[Omega],r]@@@Pp2vec,
DdagDdagPp2[a,m,\[Omega],r]@@@Pp2vec,
DdagDdagDdagPp2[a,m,\[Omega],r]@@@Pp2vec,
DPm2[a,m,\[Omega],r]@@@Pm2vec,
DDPm2[a,m,\[Omega],r]@@@Pm2vec,
DDDPm2[a,m,\[Omega],r]@@@Pm2vec,
DopDdagDdagPp2[a,m,\[Omega],r]@@@Pp2vec,
DdagDDPm2[a,m,\[Omega],r]@@@Pm2vec
};
Return[
{
MRProjector[2,"l+l+"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm2,Bmatp2}],
MRProjector[2,"l-l-"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm2,Bmatp2}],
MRProjector[2,"m+m+"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm2,Bmatp2}],
MRProjector[2,"m-m-"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm2,Bmatp2}],

MRProjector[2,"l+m+"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm2,Bmatp2}],
MRProjector[2,"l+m-"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm2,Bmatp2}],
MRProjector[2,"l-m+"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm2,Bmatp2}],
MRProjector[2,"l-m-"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm2,Bmatp2}],

MRProjector[2,"l+l-"][a,m,\[Omega],lmax,r,Pvecs,{Pp2vec,Pm2vec},{Bmatm2,Bmatp2}],

0 Pp2vec[[All,1]]
}
]
];


(* ::Subsection::Closed:: *)
(*Spin 0 calculation*)


h\[Kappa]MixingMatrix[a_,m_,\[Omega]_,lmax_,{\[Gamma]mat_,eigs0_}]:=Module[{},
a^2*Table[
If[j!=k,
If[j>=Abs[m]&&k>=Abs[m],
\[Gamma]mat[[j+1,k+1]]/(2*(eigs0[[j+1]]-eigs0[[k+1]]))
,
0
]
,
0
]
,{j,0,lmax},{k,0,lmax}]
];


(* ::Subsubsection::Closed:: *)
(*Spin 0 derivatives*)


DoprDop0[a_,m_,\[Omega]_,r_]:=Function[{h0,h1,h2},Evaluate[
h2 r+(h1 (a^2-2 r+2 I a m r+r^2-2 I a^2 r \[Omega]-2 I r^3 \[Omega]))/(a^2-2 r+r^2)-(h0 (-I a^3 m+a^2 m^2 r+I a m r^2+I a^4 \[Omega]-2 a^3 m r \[Omega]+2 I a^2 r^2 \[Omega]-4 I r^3 \[Omega]-2 a m r^3 \[Omega]+I r^4 \[Omega]+a^4 r \[Omega]^2+2 a^2 r^3 \[Omega]^2+r^5 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


DdagrDdag0[a_,m_,\[Omega]_,r_]:=Function[{h0,h1,h2},Evaluate[
h2 r+(h1 (a^2-2 r-2 I a m r+r^2+2 I a^2 r \[Omega]+2 I r^3 \[Omega]))/(a^2-2 r+r^2)-(h0 (I a^3 m+a^2 m^2 r-I a m r^2-I a^4 \[Omega]-2 a^3 m r \[Omega]-2 I a^2 r^2 \[Omega]+4 I r^3 \[Omega]-2 a m r^3 \[Omega]-I r^4 \[Omega]+a^4 r \[Omega]^2+2 a^2 r^3 \[Omega]^2+r^5 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


DopDop0[a_,m_,\[Omega]_,r_]:=Function[{h0,h1,h2},Evaluate[
h2-(2 I h1 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)-(h0 (-2 I a m+a^2 m^2+2 I a m r+2 I a^2 \[Omega]-2 a^3 m \[Omega]-2 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


DdagDdag0[a_,m_,\[Omega]_,r_]:=Function[{h0,h1,h2},Evaluate[
h2+(2 I h1 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)-(h0 (2 I a m+a^2 m^2-2 I a m r-2 I a^2 \[Omega]-2 a^3 m \[Omega]+2 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


Dop0[a_,m_,\[Omega]_,r_]:=Function[{h0,h1,h2},Evaluate[
h1-(I h0 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)
]
]


Ddag0[a_,m_,\[Omega]_,r_]:=Function[{h0,h1,h2},Evaluate[
h1+(I h0 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)
]
]


d\[CapitalDelta]d0[a_,m_,\[Omega]_,r_]:=Function[{h0,h1,h2},Evaluate[
4 h1 (-1+r)+2 h2 (a^2-2 r+r^2)+(2 h0 (-a m+a^2 \[Omega]+r^2 \[Omega])^2)/(a^2-2 r+r^2)
]
]


(* ::Subsubsection::Closed:: *)
(*Spin 0 MRProjectors*)


(* ::Input::Initialization:: *)
MRProjector[0,"l+l+"][a_,m_,\[Omega]_,lmax_,r_,{hvec_,\[Kappa]mat_,DoprDoph_,DdagrDdagh_,DopDop\[Kappa]_,DdagDdag\[Kappa]_,Dhvec_,D\[Kappa]mat_,Ddaghvec_,Ddag\[Kappa]mat_,dd2\[Kappa]mat_,d\[Kappa]mat_},{ Bmat0_}]:=Module[
{
tmp1,tmp2
},
tmp1=-1/(I*\[Omega])*DoprDoph . Bmat0;
tmp2=-4*DopDop\[Kappa] . Bmat0;
Return[tmp1+tmp2];
]


(* ::Input::Initialization:: *)
MRProjector[0,"l-l-"][a_,m_,\[Omega]_,lmax_,r_,{hvec_,\[Kappa]mat_,DoprDoph_,DdagrDdagh_,DopDop\[Kappa]_,DdagDdag\[Kappa]_,Dhvec_,D\[Kappa]mat_,Ddaghvec_,Ddag\[Kappa]mat_,dd2\[Kappa]mat_,d\[Kappa]mat_},{ Bmat0_}]:=Module[
{
tmp1,tmp2
},
tmp1=1/(I*\[Omega])*DdagrDdagh . Bmat0;
tmp2=-4*DdagDdag\[Kappa] . Bmat0;
Return[tmp1+tmp2];
]


(* ::Input::Initialization:: *)
MRProjector[0,"m+m+"][a_,m_,\[Omega]_,lmax_,r_,{hvec_,\[Kappa]mat_,DoprDoph_,DdagrDdagh_,DopDop\[Kappa]_,DdagDdag\[Kappa]_,Dhvec_,D\[Kappa]mat_,Ddaghvec_,Ddag\[Kappa]mat_,dd2\[Kappa]mat_,d\[Kappa]mat_},{ Bmat0_}]:=Module[
{
 t0,t1,tmpSh1,tmpS\[Kappa]1
},
t0=Bmat0 . (\[CapitalLambda]hat[lmax,m]-2*a \[Omega]*\[Lambda]hat[lmax,m] . SinMixMatrix[1,2,lmax,m]+(a \[Omega])^2*SinMixMatrix[0,2,lmax,m]) . CosMixMatrix[2,1,lmax,m];t1=Bmat0 . (\[Lambda]hat[lmax,m] . SinMixMatrix[1,2,lmax,m]-a \[Omega]*SinMixMatrix[0,2,lmax,m]);tmpSh1=(t0+t1);
tmpS\[Kappa]1=Bmat0 . (\[CapitalLambda]hat[lmax,m]-2*a \[Omega]*\[Lambda]hat[lmax,m] . SinMixMatrix[1,2,lmax,m]+(a \[Omega])^2*SinMixMatrix[0,2,lmax,m]);
Return[(a/\[Omega])*hvec . tmpSh1-4*\[Kappa]mat . tmpS\[Kappa]1]
]


(* ::Input::Initialization:: *)
MRProjector[0,"m-m-"][a_,m_,\[Omega]_,lmax_,r_,{hvec_,\[Kappa]mat_,DoprDoph_,DdagrDdagh_,DopDop\[Kappa]_,DdagDdag\[Kappa]_,Dhvec_,D\[Kappa]mat_,Ddaghvec_,Ddag\[Kappa]mat_,dd2\[Kappa]mat_,d\[Kappa]mat_},{ Bmat0_}]:=Module[
{
 t0,t1,tmpSh2,tmpS\[Kappa]2
},
t0=Bmat0 . (\[CapitalLambda]hat[lmax,m]-2*a \[Omega]*\[Lambda]hat[lmax,m] . SinMixMatrix[-1,-2,lmax,m]+(a \[Omega])^2*SinMixMatrix[0,-2,lmax,m]) . CosMixMatrix[-2,1,lmax,m];t1=-Bmat0 . (\[Lambda]hat[lmax,m] . SinMixMatrix[-1,-2,lmax,m]-a \[Omega]*SinMixMatrix[0,-2,lmax,m]);tmpSh2=(t0+t1);
tmpS\[Kappa]2=Bmat0 . (\[CapitalLambda]hat[lmax,m]-2*a \[Omega]*\[Lambda]hat[ lmax,m] . SinMixMatrix[-1,-2,lmax,m]+(a \[Omega])^2*SinMixMatrix[0,-2,lmax,m]);
Return[-(a/\[Omega])*hvec . tmpSh2-4* \[Kappa]mat . tmpS\[Kappa]2]
]


(* ::Input::Initialization:: *)
MRProjector[0,"l+m+"][a_,m_,\[Omega]_,lmax_,r_,{hvec_,\[Kappa]mat_,DoprDoph_,DdagrDdagh_,DopDop\[Kappa]_,DdagDdag\[Kappa]_,Dhvec_,D\[Kappa]mat_,Ddaghvec_,Ddag\[Kappa]mat_,dd2\[Kappa]mat_,d\[Kappa]mat_},{ Bmat0_}]:=Module[
{
t0,sinmix,cosmix,\[CapitalSigma]mat,MM0,tmp,term1,term2, term3, term4
},
t0=(r^2*idmat[lmax,m]+a^2*CosMixMatrix[1,2,lmax,m]);
sinmix=SinMixMatrix[0,1,lmax,m];
cosmix=CosMixMatrix[1,1,lmax,m];
MM0=(\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix);
term1=1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0;
term2=-a/\[Omega]*(r*Dhvec . Bmat0 . sinmix-hvec . Bmat0 . MM0 . cosmix);
tmp=(r*idmat[lmax,m]+I*a*cosmix);
term3=4*D\[Kappa]mat . Bmat0 . MM0 . tmp;
term4=-4*( \[Kappa]mat . Bmat0 . MM0+I*a*D\[Kappa]mat . Bmat0 . sinmix);
Return[term1+term2+term3+term4];
]


(* ::Input::Initialization:: *)
MRProjector[0,"l+m-"][a_,m_,\[Omega]_,lmax_,r_,{hvec_,\[Kappa]mat_,DoprDoph_,DdagrDdagh_,DopDop\[Kappa]_,DdagDdag\[Kappa]_,Dhvec_,D\[Kappa]mat_,Ddaghvec_,Ddag\[Kappa]mat_,dd2\[Kappa]mat_,d\[Kappa]mat_},{ Bmat0_}]:=Module[
{
t0,sinmix,cosmix,\[CapitalSigma]mat,MM0,tmp,term1,term2, term3, term4
},
t0=(r^2*idmat[lmax,m]+a^2*CosMixMatrix[-1,2,lmax,m]);
sinmix=SinMixMatrix[0,-1,lmax,m];
cosmix=CosMixMatrix[-1,1,lmax,m];
MM0=(\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix);
term1=-1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0;
term2=a/\[Omega]*(r*Dhvec . Bmat0 . sinmix+hvec . Bmat0 . MM0 . cosmix);tmp=(r*idmat[lmax,m]-I*a*cosmix);term3=-4*(D\[Kappa]mat . Bmat0 . MM0 . tmp);term4=4*(\[Kappa]mat . Bmat0 . MM0+I*a*D\[Kappa]mat . Bmat0 . sinmix);Return[term1+term2+term3+term4];
]


(* ::Input::Initialization:: *)
MRProjector[0,"l-m+"][a_,m_,\[Omega]_,lmax_,r_,{hvec_,\[Kappa]mat_,DoprDoph_,DdagrDdagh_,DopDop\[Kappa]_,DdagDdag\[Kappa]_,Dhvec_,D\[Kappa]mat_,Ddaghvec_,Ddag\[Kappa]mat_,dd2\[Kappa]mat_,d\[Kappa]mat_},{ Bmat0_}]:=Module[
{
t0,sinmix,cosmix,\[CapitalSigma]mat,MM0,tmp,term1,term2, term3, term4
},
t0=r^2*idmat[lmax,m]+a^2*CosMixMatrix[1,2,lmax,m];
sinmix=SinMixMatrix[0,1,lmax,m];
cosmix=CosMixMatrix[1,1,lmax,m];
MM0=(\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix);

term1=-1/(2*I*\[Omega])*Ddaghvec . Bmat0 . MM0 . t0;term2=-a/\[Omega]*(r*Ddaghvec . Bmat0 . sinmix-hvec . Bmat0 . MM0 . cosmix);tmp=(r*idmat[lmax,m]-I*a*cosmix);term3=4*( Ddag\[Kappa]mat . Bmat0 . MM0 . tmp);term4=-4*(\[Kappa]mat . Bmat0 . MM0-I*a* Ddag\[Kappa]mat . Bmat0 . sinmix);Return[term1+term2+term3+term4];
]


(* ::Input::Initialization:: *)
MRProjector[0,"l-m-"][a_,m_,\[Omega]_,lmax_,r_,{hvec_,\[Kappa]mat_,DoprDoph_,DdagrDdagh_,DopDop\[Kappa]_,DdagDdag\[Kappa]_,Dhvec_,D\[Kappa]mat_,Ddaghvec_,Ddag\[Kappa]mat_,dd2\[Kappa]mat_,d\[Kappa]mat_},{ Bmat0_}]:=Module[
{
t0,sinmix,cosmix,\[CapitalSigma]mat,MM0,tmp,term1,term2, term3, term4
},
t0=(r^2*idmat[lmax,m]+a^2*CosMixMatrix[-1,2,lmax,m]);
sinmix=SinMixMatrix[0,-1,lmax,m];
cosmix=CosMixMatrix[-1,1,lmax,m];
MM0=(\[Lambda]hat[lmax,m]-a*\[Omega]*sinmix);

term1=1/(2*I*\[Omega])*Ddaghvec . Bmat0 . MM0 . t0;
term2=a/\[Omega]*(r*Ddaghvec . Bmat0 . sinmix+hvec . Bmat0 . MM0 . cosmix);tmp=(r*idmat[lmax,m]+I*a*cosmix);term3=-4*( Ddag\[Kappa]mat . Bmat0 . MM0 . tmp);term4=4*( \[Kappa]mat . Bmat0 . MM0-I*a*Ddag\[Kappa]mat . Bmat0 . sinmix);
Return[term1+term2+term3+term4]
]


(* ::Input::Initialization:: *)
MRProjector[0,"l+l-"][a_,m_,\[Omega]_,lmax_,r_,{hvec_,\[Kappa]mat_,DoprDoph_,DdagrDdagh_,DopDop\[Kappa]_,DdagDdag\[Kappa]_,Dhvec_,D\[Kappa]mat_,Ddaghvec_,Ddag\[Kappa]mat_,dd2\[Kappa]mat_,d\[Kappa]mat_},{ Bmat0_}]:=Module[
{
cosmix,\[CapitalSigma]mat,tmp1,tmp2,tmp3,term1,term2, term3, term4
},
cosmix=CosMixMatrix[0,1,lmax,m];
\[CapitalSigma]mat=(r^2*idmat[lmax,m]+a^2*CosMixMatrix[0,2,lmax,m]);tmp3=\[Lambda]hat[lmax,m] . (SinMixMatrix[-1,0,lmax,m]-SinMixMatrix[1,0,lmax,m]);term1=-hvec . Bmat0 . (((\[Omega]*(r^2+a^2)-a*m)/\[Omega])*idmat[lmax,m]-2*\[CapitalSigma]mat) . \[CapitalSigma]mat;term2=-2*dd2\[Kappa]mat . Bmat0 . \[CapitalSigma]mat;term3=8*r*(a^2-2 r+r^2)*d\[Kappa]mat . Bmat0;term4=4*a^2*\[Kappa]mat . Bmat0 . tmp3 . cosmix;Return[(term1+term2+term3+term4)]
]


(* ::Subsubsection::Closed:: *)
(*CalculateSpin0Contributions*)


CalculateSpin0Contributions[a_,m_,\[Omega]_,lmax_,r_,{hvec_,\[Kappa]vec_},{Bmat0_,\[Gamma]mat_,eigs0_}]:=Module[{
Pvecs,A},
A=h\[Kappa]MixingMatrix[a,m,\[Omega],lmax,{\[Gamma]mat,eigs0}];
Pvecs={
hvec[[All,1]],
\[Kappa]vec[[All,1]]+hvec[[All,1]] . A,

DoprDop0[a,m,\[Omega],r]@@@hvec,
DdagrDdag0[a,m,\[Omega],r]@@@hvec,
DopDop0[a,m,\[Omega],r]@@@\[Kappa]vec +(DopDop0[a,m,\[Omega],r]@@@hvec) . A,
DdagDdag0[a,m,\[Omega],r]@@@\[Kappa]vec +(DdagDdag0[a,m,\[Omega],r]@@@hvec) . A,

Dop0[a,m,\[Omega],r]@@@hvec,
Dop0[a,m,\[Omega],r]@@@\[Kappa]vec+(Dop0[a,m,\[Omega],r]@@@hvec) . A,
Ddag0[a,m,\[Omega],r]@@@hvec,
Ddag0[a,m,\[Omega],r]@@@\[Kappa]vec+(Ddag0[a,m,\[Omega],r]@@@hvec) . A,

d\[CapitalDelta]d0[a,m,\[Omega],r]@@@\[Kappa]vec+(d\[CapitalDelta]d0[a,m,\[Omega],r]@@@hvec) . A,
\[Kappa]vec[[All,2]]+hvec[[All,2]] . A
};
Return[
{
MRProjector[0,"l+l+"][a,m,\[Omega],lmax,r,Pvecs,{Bmat0}],
MRProjector[0,"l-l-"][a,m,\[Omega],lmax,r,Pvecs,{Bmat0}],
MRProjector[0,"m+m+"][a,m,\[Omega],lmax,r,Pvecs,{Bmat0}],
MRProjector[0,"m-m-"][a,m,\[Omega],lmax,r,Pvecs,{Bmat0}],

MRProjector[0,"l+m+"][a,m,\[Omega],lmax,r,Pvecs,{Bmat0}],
MRProjector[0,"l+m-"][a,m,\[Omega],lmax,r,Pvecs,{Bmat0}],
MRProjector[0,"l-m+"][a,m,\[Omega],lmax,r,Pvecs,{Bmat0}],
MRProjector[0,"l-m-"][a,m,\[Omega],lmax,r,Pvecs,{Bmat0}],

MRProjector[0,"l+l-"][a,m,\[Omega],lmax,r,Pvecs,{Bmat0}],
hvec[[All,1]]
}
]
];


(* ::Subsection::Closed:: *)
(*Output Generation*)


(* ::Subsubsection::Closed:: *)
(*GetComp functions*)


GetComp[2][rtry_]:=Module[{comp},
comp=Calcs2components[rtry];
Join[comp,{Table[0,{Length@comp[[1]]}]}]
];

GetComp[1][rtry_]:=Module[{comp},
Join[comp=Calcs1components[rtry],{Table[0,{Length@comp[[1]]}]}]
];

GetComp[0][rtry_]:=Module[{comp},
Calcs0components[rtry]
];


(* ::Subsubsection::Closed:: *)
(*I/O functions*)


GetArray[tbl_]:=Module[{grid,arr},
(*dim=Dimensions[tbl];
arr=Developer`ToPackedArray[Table[SetPrecision[0,prec],{dim[[1]]},{2*dim[[2]]},{dim[[3]]}]];
For[ri=1,ri<=dim[[1]],ri++,
For[qi=1,qi<=10,qi++,
For[ll=0,ll<=lmax,ll++,
arr[[ri,qi,ll+1]]=Re[tbl[[ri,qi,ll+1]]];
arr[[ri,10+qi,ll+1]]=Im[tbl[[ri,qi,ll+1]]];
];
];
];*)
Developer`ToPackedArray@Join[Re[tbl],Im[tbl],2]
];


ExportOutput[{sallL_,sallR_,s2L_,s2R_,s1L_,s1R_,s0L_,s0R_},directory_,iConfig_,dformat_]:=Module[{fn},
fn=directory<>"data/lm_in_"<>ToString[iConfig]<>".bin";
Export[fn,GetArray[sallL],dformat];
fn=directory<>"data/lm_up_"<>ToString[iConfig]<>".bin";
Export[fn,GetArray[sallR],dformat];
(* Save the individual spin contributions as well. *)
fn=directory<>"data/lm_in_s2_"<>ToString[iConfig]<>".bin";
Export[fn,GetArray[s2L],dformat];
fn=directory<>"data/lm_up_s2_"<>ToString[iConfig]<>".bin";
Export[fn,GetArray[s2R],dformat];
fn=directory<>"data/lm_in_s1_"<>ToString[iConfig]<>".bin";
Export[fn,GetArray[s1L],dformat];
fn=directory<>"data/lm_up_s1_"<>ToString[iConfig]<>".bin";
Export[fn,GetArray[s1R],dformat];
fn=directory<>"data/lm_in_s0_"<>ToString[iConfig]<>".bin";
Export[fn,GetArray[s0L],dformat];
fn=directory<>"data/lm_up_s0_"<>ToString[iConfig]<>".bin";
Export[fn,GetArray[s0R],dformat];
]


(* ::Section:: *)
(*Main Function*)


Options[MetricReconstructRadiative]={
WorkingPrecision->32,
"nterms"-> 8, (* Number of terms in the spherical expansion. *)
"kapord"->6, (* \[Kappa]ord is the maximum order of series expansion of \[Kappa] in spheroidal harmonics. *)
"rinf"->1000., (* rinf is the radius at which the series solutions for \[Kappa]_up should set the initial conditions for the integrator. *)
"rmax"->100.,(* rmax is the maximum value at which the interpolating function can be used, whereas rinf is the starting value for the numerical integration, i.e. the radius at which the initial condition for the UP solutions of kappa is set, using the series expansion. *)
"xhor"->0.001,
"inford"->5, (* The order of the expansion at infinity for the UP function. *)
"horord"->5, (* The order of the expansion at the horizon for the IN function. *)
AccuracyGoal->9,
"rgrid"->1, (* Use a linearly-spaced grid in the variable: 0 = rstar , 1 = r. *)
"rstmin"->Automatic,
"rstmax"->Automatic,
"nres"->4,(* Resolution in the r* direction:  dr* = M / n  (or dr = M / n). *)
"angres"->8, (* Resolution in the \[Theta] direction: number of points = nres * qres. *)
"dformat"->"Real64"(* Data format for output files. *)
}


MetricReconstructRadiative[OD_OrbitalData?CircularEquatorialQ, {lmax0_,mm0_},{rgridL_,rgridR_}, iConfig_,primarypath_,opts:OptionsPattern[]]:=Module[
{
r,\[Theta],\[Lambda],l,\[CapitalOmega]0,lmax,mm,
\[CapitalSigma],K,M,\[CapitalDelta],a,m,\[Omega],\[CapitalLambda],\[Lambda]0,
\[CapitalDelta]subs,Ksubs,\[CapitalDelta]Ksubs,\[CapitalDelta]Ksimp,
QQ,CC,Cstar,AA,BB,
Dop,Ddag,Lop,Ldag,
lpvec,lmvec,mpvec,mmvec,gup,
teuks2subs,teuks1subs,teuks0subs,
Pm2,Pp2,Sm2,Sp2,Pm1,Pp1,Sm1,Sp1,P0,S0,h0,\[Kappa]0,
Pm2subs,Pp2subs,Pm1subs,Pp1subs,S0subs,h0subs,\[Kappa]0subs,
Pm2repl,Pp2repl,Sm2repl,Sp2repl,Pm1repl,Pp1repl,Sm1repl,Sp1repl,
sgn,
prec,a0,r0,nterms,\[Omega]0,a\[Omega]0,
\[CapitalDelta]Kr0subs,awsubs,
lmins2,lmins1,lmins0,lmins,
rprmvals,rp,rm,
\[Kappa]ord,rinf,rmax,xhor,inford,horord,accgoal,rstarmin,rstarmax,nres,qres,rhor,rmin,rupmin,rupmax,rgrid,dformat,
Bmats,eigs,Bmatm2,Bmatm1,Bmat0,Bmatp1,Bmatp2,eigs2,eigs1,eigs0,Ys,Ym2s,Ym1s,Y0s,Yp1s,Yp2s,
rsL,rsR,
\[Gamma]mat,\[Lambda]hat,\[Lambda]2hat,signmat,\[CapitalLambda]hat,idmat,
alplmterm1A,alplmterm1B,alplmterm2A,alplmterm2B
},
(*HelperFunctions to be divorced from main loop*)
HeadFn[head_,ll_]:={head[r]->head[ll][r],head'[r]->head[ll]'[r],head''[r]->head[ll]''[r],head'''[r]->head[ll]'''[r],head''''[r]->head[ll]''''[r]};
EvaluateRHS[set_,rval_?NumericQ]:=Module[{kk},Table[Keys[set[[kk]]]->(Values[set[[kk]]]/.{r->rval}),{kk,1,Length[set]}]];

Print["1. Preliminaries."];
\[CapitalSigma]=r^2+a^2 Cos[\[Theta]]^2;
\[CapitalDelta]subs={\[CapitalDelta][r]->r^2-2*M*r+a^2,\[CapitalDelta]'[r]->2*(r-M),\[CapitalDelta]''[r]->2,\[CapitalDelta]'''[r]->0};
Ksubs={K[r]->\[Omega]*(r^2+a^2)-a*m,K'[r]-> 2*\[Omega]*r,K''[r]->2*\[Omega],K'''[r]->0};
\[CapitalDelta]Ksubs=Join[\[CapitalDelta]subs,Ksubs];
\[CapitalDelta]Ksimp={\[CapitalDelta]'[r]->2*(r-M),\[CapitalDelta]''[r]->2,\[CapitalDelta]'''[r]->0,K'[r]-> 2*\[Omega]*r,K''[r]->2*\[Omega],K'''[r]->0};
QQ=m/Sin[\[Theta]]-a*\[Omega]*Sin[\[Theta]];

Dop[n_,R_]:=D[R,r]-I*K[r]/\[CapitalDelta][r]*R+n*\[CapitalDelta]'[r]/\[CapitalDelta][r]*R;
Ddag[n_,R_]:=D[R,r]+I*K[r]/\[CapitalDelta][r]*R+n*\[CapitalDelta]'[r]/\[CapitalDelta][r]*R;
Lop[n_,R_]:=D[R,\[Theta]]+QQ*R+n*Cot[\[Theta]]*R;
Ldag[n_,R_]:=D[R,\[Theta]]-QQ*R+n*Cot[\[Theta]]*R;
Dop[R_]:=Dop[0,R];
Ddag[R_]:=Ddag[0,R];
Lop[R_]:=Lop[0,R];
Ldag[R_]:=Ldag[0,R];

(* Inverse metric *)
lpvec={(r^2+a^2)/\[CapitalDelta][r],1,0,a/\[CapitalDelta][r]};
lmvec={-(r^2+a^2)/\[CapitalDelta][r],1,0,-a/\[CapitalDelta][r]};
mpvec={I*a*Sin[\[Theta]],0,1,I/Sin[\[Theta]]};
mmvec={-I*a*Sin[\[Theta]],0,1,-I/Sin[\[Theta]]};
gup=Table[1/(2*\[CapitalSigma])*(\[CapitalDelta][r]*(lpvec[[i]]*lmvec[[j]]+lmvec[[i]]*lpvec[[j]])+mpvec[[i]]*mmvec[[j]]+mmvec[[i]]*mpvec[[j]]),{i,1,4},{j,1,4}];
Simplify[Det[gup]+1/(\[CapitalSigma]^2*Sin[\[Theta]]^2)];

(* Teuksolky equations for spin 2 *)
teuks2subs=DefTeukSubs[2,{Pm2,Pp2,Sm2,Sp2},{\[CapitalDelta],K},{M,a,m,\[Omega],\[CapitalLambda]},{r,\[Theta]}];
Pm2subs=teuks2subs[[1;;3]];
Pp2subs=teuks2subs[[4;;6]];

(* Teukolsky-Starobinskii identities for spin 2 *)
Pp2repl=DefTSIP[2,"pm",{Pm2,Pp2},{\[CapitalDelta],K},{M,a,m,\[Omega],\[CapitalLambda],CC},{r}];
Pm2repl=DefTSIP[2,"mp",{Pm2,Pp2},{\[CapitalDelta],K},{M,a,m,\[Omega],\[CapitalLambda],Cstar},{r}];
Sp2repl=DefTSIS[2,"pm",{Sm2,Sp2},{\[CapitalDelta],K},{M,a,m,\[Omega],\[CapitalLambda],AA},{\[Theta]}];
Sm2repl=DefTSIS[2,"mp",{Sm2,Sp2},{\[CapitalDelta],K},{M,a,m,\[Omega],\[CapitalLambda],AA},{\[Theta]}];


(* Teukolsky equations for spin 1. *)
teuks1subs=DefTeukSubs[1,{Pm1,Pp1,Sm1,Sp1},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda]},{r,\[Theta]}];
Pm1subs=teuks1subs[[1;;3]];
Pp1subs=teuks1subs[[4;;6]];

(* Teukolsky-Starobinskii identities for spin 1. *)
Pp1repl=DefTSIP[1,"pm",{Pm1,Pp1},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda],BB},{r}];
Pm1repl=DefTSIP[1,"mp",{Pm1,Pp1},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda],BB},{r}];
Sp1repl=DefTSIS[1,"pm",{Sm1,Sp1},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda],BB},{\[Theta]}];
Sm1repl=DefTSIS[1,"mp",{Sm1,Sp1},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda],BB},{\[Theta]}];

(* Teukolsky equations for spin 0. *)
teuks0subs=DefTeukSubs[0,{P0,S0},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda]0},{r,\[Theta]}];
S0subs=teuks0subs[[4;;6]];

(* Teukolsky equations for trace, and for kappa. *)
h0subs=DefTeukSubs[0,{h0,S0},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda]0},{r,\[Theta]}][[1;;3]];
\[Kappa]0subs=DefTeukSubs["\[Kappa]0",{\[Kappa]0,h0},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda]0},{r,\[Theta]}];(* This needs further thought! *)

Print["2. Projection + Jumps."];

(* Read the parameters from a configuration file *)
(* Read the config file *)
extractUpToLastHyphen[str_String] := StringJoin[Riffle[Most[StringSplit[str, "-"]], "-"]];
extractUpToSecondLastHyphen[str_String] := StringJoin[Riffle[Most[StringSplit[extractUpToLastHyphen[str], "-"]], "-"]];
directory=(primarypath<>"/data"<>extractUpToSecondLastHyphen[ToString[iConfig]]<>"/data"<>extractUpToLastHyphen[ToString[iConfig]]<>"/");
If[!DirectoryQ[directory], CreateDirectory[directory, CreateIntermediateDirectories -> True]];

(* Read in the numerical parameters. All parameters should be read in this section only (no longer spread throughout the code). Unlike the above, the following sections are done numerically *)
prec=OptionValue[WorkingPrecision]; (* Number of digits to use where required. This must be read in first. *)
a0=SetPrecision[ODspin[OD],prec];
If[PossibleZeroQ[a0],a0=0]; (* In the Schwarzschild case, a0 should be the integer zero. *)
r0=SetPrecision[ODsemilatusrectum[OD],prec];
lmax=lmax0;
mm=mm0;
nterms=OptionValue["nterms"];(* Number of terms in the spherical expansion. *)
{a0,r0,mm,lmax,nterms};

(* IMPORTANT: Be cautious of using Simplify after introducing the numerical values, as this can lead to huge loss of precision. All Simplify statements must come before replacing symbols numerical values. *)
\[CapitalOmega]0=KerrAzimuthalFrequency[OD];
\[Omega]0=mm \[CapitalOmega]0;
a\[Omega]0=If[a0==0,0,a0*\[Omega]0];
\[CapitalDelta]Kr0subs=Map[#[[1]]->(#[[2]]/.{a->a0,r->r0,M->1,m->mm,\[Omega]->\[Omega]0})&,Join[\[CapitalDelta]subs,Ksubs]];
awsubs=Join[SetPrecision[{a->a0,\[Omega]->\[Omega]0},prec],{m->mm,M->1}];
lmins0=Abs[mm];
lmins1=Max[1,Abs[mm]];
lmins2=Max[2,Abs[mm]];
lmins={lmins2,lmins1,lmins0,lmins1,lmins2};
rprmvals={rp->1+Sqrt[1-a^2],rm->1-Sqrt[1-a^2]}/.awsubs;

(* To handle the \[Kappa] functions, satisfying equations sourced by the radial profile of the trace, I have written some bespoke code (!) *)
\[Kappa]ord=OptionValue["kapord"]; (* \[Kappa]ord is the maximum order of series expansion of \[Kappa] in spheroidal harmonics. *)
rinf=SetPrecision[OptionValue["rinf"],prec]; (* rinf is the radius at which the series solutions for \[Kappa]_up should set the initial conditions for the integrator. *)
rmax=SetPrecision[OptionValue["rmax"],prec];   (* rmax is the maximum value at which the interpolating function can be used, whereas rinf is the starting value for the numerical integration, i.e. the radius at which the initial condition for the UP solutions of kappa is set, using the series expansion. *)
xhor=SetPrecision[OptionValue["xhor"],prec];
rhor=(rp+xhor)/.rprmvals;
rmin =rhor + SetPrecision[10^-1,prec];  (* <--- may need to think again about this choice. *)
inford=OptionValue["inford"]; (* The order of the expansion at infinity for the UP function. *)
horord=OptionValue["horord"]; (* The order of the expansion at the horizon for the IN function. *)
accgoal=OptionValue["AccuracyGoal"];
rupmin=r0/.awsubs;
rupmax=rmax;
rgrid=SetPrecision[OptionValue["rgrid"],prec]; (* Use a linearly-spaced grid in the variable: 0 = rstar , 1 = r. *)
rstarmin=If[OptionValue["rstmin"]===Automatic,SetPrecision[rmin,prec],SetPrecision[OptionValue["rstmin"],prec]]; (* If rgrid=1 then these will be interpreted as rmin and rmax, instead of rstarmin and rstarmax. *)
rstarmax=If[OptionValue["rstmax"]===Automatic,SetPrecision[rmax,prec],SetPrecision[OptionValue["rstmax"],prec]];
nres=OptionValue["nres"];  (* Resolution in the r* direction:  dr* = M / n  (or dr = M / n). *)
qres=OptionValue["angres"]; (* Resolution in the \[Theta] direction: number of points = nres * qres. *)
dformat=ToString@OptionValue["dformat"]; (* Data format for output files. *)
{\[Kappa]ord,rinf,rmin,rmax,xhor,inford,horord,accgoal,rstarmin,rstarmax,nres,qres,\[CapitalOmega]0};

rsL=rgridL;
rsR=rgridR;

(* Set up the b matrices, and the eigenvalue matrices. *)
PrintTemporary["Calculating B-matrices and Eigenvalues"];
{Bmats,eigs}=Transpose@Table[
bmatEV[s,mm,a\[Omega]0,lmax,WorkingPrecision->prec,AccuracyGoal->accgoal],
{s,-2,2}
];
{Bmatm2,Bmatm1,Bmat0,Bmatp1,Bmatp2}=Bmats;
{eigs2,eigs1,eigs0}=eigs[[1;;3]];

(*Ys=Table[Table[If[ll>=lmins[[ss+3]],1/Sqrt[2\[Pi]]SpinWeightedSphericalHarmonic[ss,ll,mm][Cos[\[Theta]]]],{ll,0,lmax}],{ss,-2,2}]/.Null->0;*)
Ys=Table[Table[If[ll>=lmins[[ss+3]],SpinWeightedSphericalHarmonicY[ss,ll,mm,\[Theta],0],0],{ll,0,lmax}],{ss,-2,2}];
{Ym2s,Ym1s,Y0s,Yp1s,Yp2s}=Ys;

(* Now define the S matrices, following the notes *)
\[Lambda]hat=DiagonalMatrix[Table[If[ll<lmins0,0,Sqrt[ll*(ll+1)]],{ll,0,lmax}]];
\[Lambda]2hat=DiagonalMatrix[Table[If[ll<lmins1,0,Sqrt[(ll-1)*(ll+2)]],{ll,0,lmax}]];
\[CapitalLambda]hat=\[Lambda]hat . \[Lambda]2hat;
signmat=DiagonalMatrix[Table[If[ll<lmins0,0,(-1)^(ll+mm)],{ll,0,lmax}]];



(* N.B. I have multiplied calD by -1 because we use the -+++ signature. The results for mode functions and fluxes then agree with the BHPT results using TeukolskyPointParticleMode.
See check_sourceterms_with_BHPtoolkit.nb. *)
PrintTemporary["Jumps"];
ABDsubs={calD->4*\[Pi]/(r^2*ut),calA->1/\[CapitalDelta][r]*(EE*(r^2+a^2)-a*LL),calB->LL-a*EE};
W0=m*(\[CapitalOmega]*(r^2+a^2)-a)/\[CapitalDelta][r];
Q0=m*(1/Sin[\[Theta]]-a*\[CapitalOmega]*Sin[\[Theta]])/.{\[Theta]->\[Pi]/2};
calC2=calB^2*S0;
calC1=2*calB*(I*calA*(S0p+Q0*S0)+calB*(-I*W0+1/r)*S0);
calC0=2*calA*(I*calB*(2/r-I*W0)+calA*(-Q0+I*a/r))*(S0p+Q0*S0)+(calA^2*\[CapitalLambda]+calB^2*(I*D[W0,r]-2*I*W0/r-W0^2))*S0;
V0minus4=\[CapitalDelta][r]*(W0^2+I*D[W0,r])-I*D[\[CapitalDelta][r],r]*W0+6*I*m*\[CapitalOmega]*r-\[CapitalLambda];
plussubs={S0->Sp2[\[Pi]/2],S0p->Sp2'[\[Pi]/2]};
minussubs={m->-mrepl,S0->Sm2[\[Pi]/2],S0p->Sm2'[\[Pi]/2]};
Jump0p=calD*\[CapitalDelta][r]*(calC1-D[\[CapitalDelta][r],r]/\[CapitalDelta][r]*calC2)/.plussubs;
Jump1p=calD*\[CapitalDelta][r]*(calC0-V0minus4/\[CapitalDelta][r]*calC2)/.plussubs;
Jump0m=calD*\[CapitalDelta][r]*((calC1-D[\[CapitalDelta][r],r]/\[CapitalDelta][r]*calC2)/.minussubs)/.{mrepl->m};
Jump1m=calD*\[CapitalDelta][r]*((calC0-V0minus4/\[CapitalDelta][r]*calC2)/.minussubs)/.{mrepl->m};

mysubs={LL->r0*v*(1-2*a*v^3+a^2*v^4)/Sqrt[1-3*v^2+2*a*v^3],EE->(1-2*v^2+a*v^3)/Sqrt[1-3*v^2+2*a*v^3],\[CapitalOmega]->v^3/(1+a*v^3)}/.{v->1/Sqrt[r0]}/.awsubs;
utsubs={ut->-EE*gup[[1,1]]+LL*gup[[1,4]]}/.\[CapitalDelta]Kr0subs/.{\[Theta]->\[Pi]/2,r->r0}/.mysubs/.{a->a0};
orbitsubs=Join[mysubs,utsubs];
\[CapitalOmega]0-(\[CapitalOmega]/.mysubs)/.awsubs;
Jumpnum0p=Jump0p/.ABDsubs/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs/.orbitsubs//Simplify;
Jumpnum1p=Jump1p/.ABDsubs/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs/.orbitsubs//Simplify;
Jumpnum0m=Jump0m/.ABDsubs/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs/.orbitsubs//Simplify;
Jumpnum1m=Jump1m/.ABDsubs/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs/.orbitsubs//Simplify;


If[a\[Omega]0==0,
s2jumps=Table[{Pp2j0[ll]->Jumpnum0p,Pp2j1[ll]->Jumpnum1p,Pm2j0[ll]->Jumpnum0m,Pm2j1[ll]->Jumpnum1m}/.
{Sp2[\[Pi]/2]->SpinWeightedSpheroidalHarmonicS[2,ll,mm,a\[Omega]0][\[Pi]/2,0],
Sm2[\[Pi]/2]->SpinWeightedSpheroidalHarmonicS[-2,ll,mm,a\[Omega]0][\[Pi]/2,0],Sp2'[\[Pi]/2]->D[SpinWeightedSpheroidalHarmonicS[2,ll,mm,a\[Omega]0][\[Theta],0],\[Theta]]/.{\[Theta]->\[Pi]/2},
Sm2'[\[Pi]/2]->D[SpinWeightedSpheroidalHarmonicS[-2,ll,mm,a\[Omega]0][\[Theta],0],\[Theta]]/.{\[Theta]->\[Pi]/2},\[CapitalLambda]->Get\[CapitalLambda][ll,mm,a\[Omega]0]},{ll,lmins2,lmax}]//Join//Flatten
,
s2jumps=Table[{Pp2j0[ll]->Jumpnum0p,Pp2j1[ll]->Jumpnum1p,Pm2j0[ll]->Jumpnum0m,Pm2j1[ll]->Jumpnum1m}/.
{Sp2[\[Pi]/2]->SpinWeightedSpheroidalHarmonicS[2,ll,mm,a\[Omega]0,Method->{"SphericalExpansion","NumTerms"->nterms}][\[Pi]/2,0],
Sm2[\[Pi]/2]->SpinWeightedSpheroidalHarmonicS[-2,ll,mm,a\[Omega]0,Method->{"SphericalExpansion","NumTerms"->nterms}][\[Pi]/2,0],Sp2'[\[Pi]/2]->D[SpinWeightedSpheroidalHarmonicS[2,ll,mm,a\[Omega]0,Method->{"SphericalExpansion","NumTerms"->nterms}][\[Theta],0],\[Theta]]/.{\[Theta]->\[Pi]/2},
Sm2'[\[Pi]/2]->D[SpinWeightedSpheroidalHarmonicS[-2,ll,mm,a\[Omega]0,Method->{"SphericalExpansion","NumTerms"->nterms}][\[Theta],0],\[Theta]]/.{\[Theta]->\[Pi]/2},\[CapitalLambda]->Get\[CapitalLambda][ll,mm,a\[Omega]0]},{ll,lmins2,lmax}]//Join//Flatten
];

(* Now for the jumps in the derivative of the trace. *)
If[a\[Omega]0==0,
s0jumps=Table[hj1[ll]->(-16*\[Pi]/(ut*\[CapitalDelta][r])*S0[\[Pi]/2]/.{S0[\[Pi]/2]->SpinWeightedSpheroidalHarmonicS[0,ll,mm,a\[Omega]0][\[Pi]/2,0]})/.\[CapitalDelta]Kr0subs/.{r->r0}/.orbitsubs,{ll,lmins0,lmax}];
,
s0jumps=Table[hj1[ll]->CleanValue[(-16*\[Pi]/(ut*\[CapitalDelta][r])*S0[\[Pi]/2]/.{S0[\[Pi]/2]->SpinWeightedSpheroidalHarmonicS[0,ll,mm,a\[Omega]0,Method->{"SphericalExpansion","NumTerms"->nterms}][\[Pi]/2,0]})/.\[CapitalDelta]Kr0subs/.{r->r0}/.orbitsubs],{ll,lmins0,lmax}];
];


sourcetermsubs=Table[Simplify[{dhlplp[ll]->-16*\[Pi]/(\[CapitalDelta][r]*ut)*calA^2*SpinWeightedSphericalHarmonicY[0,ll,mm,\[Pi]/2,0],
dhlmlm[ll]->-16*\[Pi]/(\[CapitalDelta][r]*ut)*calA^2*SpinWeightedSphericalHarmonicY[0,ll,mm,\[Pi]/2,0],
dhmpmp[ll]->16*\[Pi]/(\[CapitalDelta][r]*ut)*calB^2*SpinWeightedSphericalHarmonicY[2,ll,mm,\[Pi]/2,0],
dhmmmm[ll]->16*\[Pi]/(\[CapitalDelta][r]*ut)*calB^2*SpinWeightedSphericalHarmonicY[-2,ll,mm,\[Pi]/2,0],
d\[Rho]chlpmm[ll]->-16*\[Pi]/(\[CapitalDelta][r]*ut)*(r*I*calA*calB)*SpinWeightedSphericalHarmonicY[-1,ll,mm,\[Pi]/2,0]}/.ABDsubs]/.\[CapitalDelta]Kr0subs/.{r->r0,a->a0}/.orbitsubs,
{ll,lmins0,lmax}]//Flatten;

(* Radial functions as vectors. *)
ones=Table[1,{ll,0,lmax}];
Pp2vec=Table[If[ll>=lmins2,Pp2[ll][r],0],{ll,0,lmax}];
Pp1vec=Table[If[ll>=lmins1,Pp1[ll][r],0],{ll,0,lmax}];
hvec=Table[If[ll>=lmins0,h[ll][r],0],{ll,0,lmax}];
Pm1vec=Table[If[ll>=lmins1,Pm1[ll][r],0],{ll,0,lmax}];
Pm2vec=Table[If[ll>=lmins2,Pm2[ll][r],0],{ll,0,lmax}];

(* \[Kappa] needs special handling *)
\[Gamma]mat=Bmat0 . CosMixMatrix[0,2,lmax,mm] . (Transpose@Bmat0);
tmp1=DiagonalMatrix[Table[If[ll>=lmins0,\[Kappa][ll][r],0],{ll,0,lmax}]];
tmp2=a^2*Table[If[j!=k,If[j>=lmins0&&k>=lmins0,
\[Gamma]mat[[j+1,k+1]]*h[j][r]/(2*(eigs0[[j+1]]-eigs0[[k+1]]))
,0],0],{j,0,lmax},{k,0,lmax}];
\[Kappa]mat=(tmp1+tmp2)/.awsubs;

(* Awkwardness: to keep the expressions compact, I would like to replace 'r' with the numerical value 'r0', but not in the radial functions themselves, as these I will use as unknowns. So I will replace the radial functions with jump expressions at an early stage. *)
Pp2jump=Table[If[ll>=lmins2,{Pp2[ll][r]->Pp2j0[ll],Pp2[ll]'[r]->Pp2j1[ll]},{}],{ll,0,lmax}]//Flatten;
Pm2jump=Table[If[ll>=lmins2,{Pm2[ll][r]->Pm2j0[ll],Pm2[ll]'[r]->Pm2j1[ll]},{}],{ll,0,lmax}]//Flatten;
Pp1jump=Table[If[ll>=lmins1,{Pp1[ll][r]->Pp1j0[ll],Pp1[ll]'[r]->Pp1j1[ll]},{}],{ll,0,lmax}]//Flatten;
Pm1jump=Table[If[ll>=lmins1,{Pm1[ll][r]->Pm1j0[ll],Pm1[ll]'[r]->Pm1j1[ll]},{}],{ll,0,lmax}]//Flatten;
hjump=Table[If[ll>=lmins0,{h[ll][r]->0,h[ll]'[r]->hj1[ll]},{}],{ll,0,lmax}]//Flatten;
\[Kappa]jump=Table[If[ll>=lmins0,{\[Kappa][ll][r]->\[Kappa]j0[ll],\[Kappa][ll]'[r]->\[Kappa]j1[ll]},{}],{ll,0,lmax}]//Flatten;



teuk\[Kappa]=Dop[\[CapitalDelta][r]*Ddag[\[Kappa][r]]]-2*I*\[Omega]*r*\[Kappa][r]-\[Lambda]0*\[Kappa][r]-1/2*(r^2+a^2*\[Gamma]ll[l])*h[r];  (* This coefficient can be filled in later *)
\[Kappa]subs2=Solve[teuk\[Kappa]==0,\[Kappa]''[r]]//First//Simplify;
\[Kappa]subs3=D[\[Kappa]subs2,r]/.\[Kappa]subs2/.\[CapitalDelta]Ksimp//Simplify;
\[Kappa]subs4=D[\[Kappa]subs3,r]/.\[Kappa]subs2/.\[CapitalDelta]Ksimp//Simplify;
\[Kappa]subs=Join[\[Kappa]subs2,\[Kappa]subs3,\[Kappa]subs4]//Simplify;
(*mat=CosMixMatrix[0,2,lmax,mm];*)
\[Gamma]llsubs=Table[\[Gamma]ll[ll]->\[Gamma]mat[[ll+1,ll+1]],{ll,0,lmax}];

(* 
  Now we need replacement rules from the Teukolsky and Teukolsky-Starobinskii equations, on an industrial scale.
Remember the (-1)^(l+m) factor in the TS identities for spin 2.
*)
(* First the Teukolsky equations *)
Clear[ll];
h0tohrepl={h0[r]->h[r],h0'[r]->h'[r],h0''[r]->h''[r],h0'''[r]->h'''[r],h0''''[r]->h''''[r]};
Pm2lsubs= Table[Pm2subs/.awsubs/.{\[CapitalLambda]->eigs2[[ll+1]]}/.HeadFn[Pm2,ll],{ll,lmins2,lmax}]//Flatten;
Pp2lsubs= Table[Pp2subs/.awsubs/.{\[CapitalLambda]->eigs2[[ll+1]]}/.HeadFn[Pp2,ll],{ll,lmins2,lmax}]//Flatten;
Pm1lsubs= Table[Pm1subs/.awsubs/.{\[Lambda]->eigs1[[ll+1]]}/.HeadFn[Pm1,ll],{ll,lmins1,lmax}]//Flatten;
Pp1lsubs= Table[Pp1subs/.awsubs/.{\[Lambda]->eigs1[[ll+1]]}/.HeadFn[Pp1,ll],{ll,lmins1,lmax}]//Flatten;
hlsubs= Table[h0subs/.awsubs/.{\[Lambda]0->eigs0[[ll+1]]}/.h0tohrepl/.HeadFn[h,ll],{ll,lmins0,lmax}]//Flatten;
\[Kappa]lsubs= Table[\[Kappa]subs/.awsubs/.{\[Lambda]0->eigs0[[ll+1]]}/.{l->ll}/.h0tohrepl/.HeadFn[\[Kappa],ll]/.HeadFn[h,ll],{ll,lmins0,lmax}]//Flatten;

(* Now consider the Teuksolsky-Starobinskii identities *)
AArepl={AA->Sqrt[\[CapitalLambda]^2*(\[CapitalLambda]+2)^2-8*\[Omega]^2*\[CapitalLambda]*(\[Alpha]sq*(5*\[CapitalLambda]+6)-12*a^2)+144*\[Omega]^4*\[Alpha]sq^2]}/.{\[Alpha]sq->a^2-a*m/\[Omega]};
DD2wp=\[Lambda]C^2*(\[Lambda]C-2)^2+8*a*\[Omega]*(m-a*\[Omega])*(\[Lambda]C-2)*(5*\[Lambda]C-4)+48*(a*\[Omega])^2*(2*(\[Lambda]C-2)+3*(m-a*\[Omega])^2); (* Pound and Wardell's version. *)
Simplify[(AA^2-DD2wp)/.AArepl/.{\[Lambda]C->\[CapitalLambda]+2}];(* Check against Pound and Wardell. *)
BBrepl={BB->Sqrt[\[Lambda]^2 +4*a*m*\[Omega]-4*a^2*\[Omega]^2 ]};
Crepl0={CC->AA+12 I M sgn \[Omega],Cstar->AA-12 I M sgn \[Omega]};
Clear[ll];
Pm2lrepl=Table[Pm2repl/.awsubs/.HeadFn[Pp2,ll]/.HeadFn[Pm2,ll]/.Crepl0/.AArepl/.{\[CapitalLambda]->eigs2[[ll+1]]}/.{sgn->(-1)^(ll+mm)}/.awsubs,{ll,lmins2,lmax}]//Flatten;
Pm1lrepl=Table[Pm1repl/.awsubs/.HeadFn[Pp1,ll]/.HeadFn[Pm1,ll]/.BBrepl/.{\[Lambda]->eigs1[[ll+1]]}/.awsubs,{ll,lmins1,lmax}]//Flatten;
Pp2lrepl=Table[Pp2repl/.awsubs/.HeadFn[Pm2,ll]/.HeadFn[Pp2,ll]/.Crepl0/.AArepl/.{\[CapitalLambda]->eigs2[[ll+1]]}/.{sgn->(-1)^(ll+mm)}/.awsubs,{ll,lmins2,lmax}]//Flatten;
Pp1lrepl=Table[Pp1repl/.awsubs/.HeadFn[Pm1,ll]/.HeadFn[Pp1,ll]/.BBrepl/.{\[Lambda]->eigs1[[ll+1]]}/.awsubs,{ll,lmins1,lmax}]//Flatten;

(* Spin 1 projections. *)
(* I have deleted all the projections except lplp, mpmp, lmlm, mmmm and rho h_{lpmp} and rhoc h_{lpmm} *)
Timing[
idmat=DiagonalMatrix[Table[If[ll>=lmins0,1,0],{ll,0,lmax}]];
t0=-(Bmatm1+signmat . Bmatp1) . \[Lambda]hat;
t1=a\[Omega]0*(Bmatm1 . SinMixMatrix[-1,0,lmax,mm]+signmat . Bmatp1 . SinMixMatrix[1,0,lmax,mm]);
SS0=t0+t1;
Ss1lplp=SS0;
Ss1lmlm=SS0;
Ss1mpmp=Bmatp1 . (-\[Lambda]2hat+a\[Omega]0*SinMixMatrix[1,2,lmax,mm]);
Ss1mmmm=Bmatm1 . (\[Lambda]2hat-a\[Omega]0*SinMixMatrix[-1,-2,lmax,mm]);
tmp=2/\[CapitalDelta][r]*Table[Dop[-1,Pp1vec[[ll+1]]],{ll,0,lmax}];
hraws1lplp=tmp . signmat . Ss1lplp;  (* zzz edited line *)
hs1lplp=hraws1lplp/.Pp1lrepl/.awsubs;
tmp=2/\[CapitalDelta][r]*Table[Ddag[-1,Pm1vec[[ll+1]]],{ll,0,lmax}];
hraws1lmlm=tmp . Ss1lmlm;
hs1lmlm=hraws1lmlm/.awsubs;
tmp1=signmat . Table[Ddag[Pp1vec[[ll+1]]],{ll,0,lmax}]/.Pp1lrepl/.awsubs;
tmp2=Table[Dop[Pm1vec[[ll+1]]],{ll,0,lmax}]/.awsubs;
tmp=2*(tmp1+tmp2);
hraws1mpmp=tmp . signmat . Ss1mpmp;  
hraws1mmmm=-tmp . Ss1mmmm;
hs1mpmp=hraws1mpmp/.Pp1lrepl/.awsubs;
hs1mmmm=hraws1mmmm/.Pp1lrepl/.awsubs;
(* Now for the other components *)
tmp1=signmat . Table[Ddag[Pp1vec[[ll+1]]],{ll,0,lmax}]/.Pp1lrepl/.awsubs;
tmp2=Table[Dop[Pm1vec[[ll+1]]],{ll,0,lmax}]/.awsubs;
calPs1=tmp2+tmp1;
DcalPs1=Map[Dop,calPs1]/.Pm1lsubs;
(* Ss1lplp//TableForm *)
(* l+ m+ *)
cosmix=CosMixMatrix[1,1,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
line1=(signmat . (r*DcalPs1-2*calPs1)) . Bmatp1+I*a*(signmat . DcalPs1) . Bmatp1 . cosmix/.Pm1lsubs;
tmp1=1/\[CapitalDelta][r]*signmat . Pp1vec/.Pp1lrepl/.awsubs;
tmp2=-r*\[Lambda]hat-I*a*\[Lambda]hat . cosmix+(a*\[Omega]*r+2*I*a)*sinmix+I*a^2*\[Omega]*sinmix . cosmix;
line2=tmp1 . SS0 . tmp2;
\[Rho]hs1lpmp=(line1+line2)/.awsubs;
(* Now the \[Rho]c h_{l+ m-} term *)
line1=-(r*DcalPs1-2*calPs1) . Bmatm1+I*a*DcalPs1 . Bmatm1 . CosMixMatrix[-1,1,lmax,mm]/.Pm1lsubs;
tmp1=1/\[CapitalDelta][r]*signmat . Pp1vec/.Pp1lrepl/.awsubs;
tmp2=r*\[Lambda]hat-I*a*\[Lambda]hat . CosMixMatrix[-1,1,lmax,mm]-(a*\[Omega]*r+2*I*a)*SinMixMatrix[0,-1,lmax,mm]+I*a^2*\[Omega]*SinMixMatrix[0,-1,lmax,mm] . CosMixMatrix[-1,1,lmax,mm];
line2=tmp1 . SS0 . tmp2;
\[Rho]chs1lpmm=(line1+line2)/.awsubs;
];
(* Spin 2: lplp and mpmp *)
EchoTiming[
idmat=DiagonalMatrix[Table[If[ll>=lmins0,1,0],{ll,0,lmax}]];
t0=Bmatp2 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]2hat . SinMixMatrix[1,0,lmax,mm]+a\[Omega]0^2*SinMixMatrix[2,0,lmax,mm]);
t1=signmat . Bmatm2 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]2hat . SinMixMatrix[-1,0,lmax,mm]+a\[Omega]0^2*SinMixMatrix[-2,0,lmax,mm]);
Ss2lplp=t0+t1;
Ss2lmlm=signmat . Ss2lplp; 
Ss2mpmp=Bmatp2;
Ss2mmmm=signmat . Bmatm2;
hraws2lplp=-1/(6*\[Omega]^2*\[CapitalDelta][r]^2)*Pp2vec . Ss2lplp;
hraws2lmlm=-1/(6*\[Omega]^2*\[CapitalDelta][r]^2)*Pm2vec . Ss2lmlm;
hs2lplp=hraws2lplp/.awsubs/.Pp2lrepl;
hs2lmlm=hraws2lmlm/.awsubs;
tmp1=Table[Ddag[Ddag[Pp2vec[[ll+1]]]],{ll,0,lmax}];
tmp2=signmat . Table[Dop[Dop[Pm2vec[[ll+1]]]],{ll,0,lmax}];
tmp=-1/(6*\[Omega]^2)*(tmp1+tmp2);
hraws2mpmp=tmp . Ss2mpmp;
hraws2mmmm=tmp . Ss2mmmm;
hs2mpmp=hraws2mpmp/.Pp2lsubs/.Pp2lrepl/.Pm2lsubs/.awsubs;
hs2mmmm=hraws2mmmm/.Pp2lsubs/.Pp2lrepl/.Pm2lsubs/.awsubs;
];
EchoTiming[
SSplus=(Bmatp2 . (\[CapitalLambda]hat-2*a*\[Omega]*\[Lambda]2hat . SinMixMatrix[1,0,lmax,mm]+a^2*\[Omega]^2*SinMixMatrix[2,0,lmax,mm]))/.awsubs//Simplify;
SSminus=(Bmatm2 . (\[CapitalLambda]hat-2*a*\[Omega]*\[Lambda]2hat . SinMixMatrix[-1,0,lmax,mm]+a^2*\[Omega]^2*SinMixMatrix[-2,0,lmax,mm]))/.awsubs//Simplify;
DdagPp2vec=Map[Ddag,Pp2vec]/.Pp2lrepl/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DdagDdagPp2vec=Map[Ddag,DdagPp2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DdagDdagDdagPp2vec=Map[Ddag,DdagDdagPp2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DPm2vec=Map[Dop,Pm2vec]/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DDPm2vec=Map[Dop,DPm2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DDDPm2vec=Map[Dop,DDPm2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DopDdagDdagPp2vec=Map[Dop,DdagDdagPp2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DdagDDPm2vec=Map[Ddag,DDPm2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
];
(* lp mp *)
EchoTiming[
cosmix=CosMixMatrix[1,1,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
sinmix1m1=SinMixMatrix[-1,1,lmax,mm];
\[Rho]mat=r*idmat+I*a*cosmix;
MMp2=Bmatp2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[2,1,lmax,mm])/.awsubs;
MMm2=Bmatm2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[-2,-1,lmax,mm])/.awsubs;
term1=1/(12*\[Omega]^2)*(DDPm2vec . signmat . SSplus-DdagDdagPp2vec . signmat . SSminus) . (\[Lambda]hat-a*\[Omega]*sinmix)/.awsubs;
term2=-1/(12*\[Omega]^2)*(DDDPm2vec . signmat . SSplus-DopDdagDdagPp2vec . signmat . SSminus) . ((\[Lambda]hat-a*\[Omega]*sinmix) . \[Rho]mat-I*a*sinmix)/.awsubs;
term3=-1/(3*I*\[Omega])*(DDDPm2vec . signmat . MMp2 . \[Rho]mat . \[Rho]mat-2*DDPm2vec . signmat . MMp2 . \[Rho]mat+2*DPm2vec . signmat . MMp2)/.awsubs;
tmp=(\[Lambda]hat . \[Lambda]hat-2*a*\[Omega]*\[Lambda]hat . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]mat . \[Rho]mat-2*I*a*(\[Lambda]hat . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]mat-2*a^2*sinmix1m1;
term4=-1/(3*I*\[Omega]*\[CapitalDelta][r])*DdagPp2vec . (signmat . MMm2 . tmp)/.awsubs;
(* Don't forget the \[Beta] factor. *)
tmp=term1+term2+term3+term4;
\[Rho]hs2lpmp=(tmp/(-6*I*M*\[Omega]))/.awsubs;
];
(* lp mm *)
EchoTiming[
idmat=DiagonalMatrix[Table[If[ll>=lmins0,1,0],{ll,0,lmax}]];
\[Rho]cmat=r*idmat-I*a*CosMixMatrix[-1,1,lmax,mm];
sinmix=SinMixMatrix[0,-1,lmax,mm];
MMp2=Bmatp2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[2,1,lmax,mm]);
MMm2=Bmatm2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[-2,-1,lmax,mm]);
term1=-1/(12*\[Omega]^2)*(DDPm2vec . SSminus-DdagDdagPp2vec . SSplus) . (\[Lambda]hat-a*\[Omega]*sinmix);
term2=1/(12*\[Omega]^2)*(DDDPm2vec . SSminus-DopDdagDdagPp2vec . SSplus) . ((\[Lambda]hat-a*\[Omega]*sinmix) . \[Rho]cmat-I*a*sinmix);
term3=1/(3*I*\[Omega])*(DDDPm2vec . MMm2 . \[Rho]cmat . \[Rho]cmat-2*DDPm2vec . MMm2 . \[Rho]cmat+2*DPm2vec . MMm2);
sinmix1m1=SinMixMatrix[1,-1,lmax,mm];
tmp=(\[Lambda]hat . \[Lambda]hat-2*a*\[Omega]*\[Lambda]hat . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]cmat . \[Rho]cmat-2*I*a*(\[Lambda]hat . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]cmat-2*a^2*sinmix1m1;
term4=1/(3*I*\[Omega]*\[CapitalDelta][r])*DdagPp2vec . MMp2 . tmp;
tmp=term1+term2+term3+term4;
\[Rho]chs2lpmm=(tmp/(-6*I*M*\[Omega]))/.awsubs;
];

(* Spin 0 *)
(* When we apply \[Kappa]lsubs, we also should fill in the \[Gamma]ll values *)
EchoTiming[
tmp1=-1/(I*\[Omega])*Dop[r*Dop[hvec]] . Bmat0/.hlsubs/.awsubs;
tmp2a=(ones . Dop[Dop[\[Kappa]mat]]) . Bmat0/.\[Kappa]lsubs/.\[Gamma]llsubs/.hlsubs/.awsubs;
tmp2=-4*tmp2a/.awsubs;
hs0lplp=(tmp1+tmp2)/.awsubs;

tmp1=1/(I*\[Omega])*Ddag[r*Ddag[hvec]] . Bmat0/.hlsubs/.awsubs;
tmp2a=(ones . Ddag[Ddag[\[Kappa]mat]]) . Bmat0/.\[Kappa]lsubs/.\[Gamma]llsubs/.hlsubs/.awsubs;
tmp2=-4*tmp2a/.awsubs;
hs0lmlm=(tmp1+tmp2)/.awsubs;

t0=Bmat0 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]hat . SinMixMatrix[1,2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,2,lmax,mm]) . CosMixMatrix[2,1,lmax,mm];
t1=Bmat0 . (\[Lambda]hat . SinMixMatrix[1,2,lmax,mm]-a\[Omega]0*SinMixMatrix[0,2,lmax,mm]);
tmpSh=(t0+t1)/.awsubs;
tmpS\[Kappa]=Bmat0 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]hat . SinMixMatrix[1,2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,2,lmax,mm])/.awsubs;
hs0mpmp=((a/\[Omega])*hvec . tmpSh-4*(ones . \[Kappa]mat) . tmpS\[Kappa])/.awsubs;

t0=Bmat0 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]hat . SinMixMatrix[-1,-2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,-2,lmax,mm]) . CosMixMatrix[-2,1,lmax,mm];
t1=-Bmat0 . (\[Lambda]hat . SinMixMatrix[-1,-2,lmax,mm]-a\[Omega]0*SinMixMatrix[0,-2,lmax,mm]);
tmpSh=(t0+t1)/.awsubs;
tmpS\[Kappa]=Bmat0 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]hat . SinMixMatrix[-1,-2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,-2,lmax,mm])/.awsubs;
hs0mmmm=(-(a/\[Omega])*hvec . tmpSh-4*(ones . \[Kappa]mat) . tmpS\[Kappa])/.awsubs;
];

(* lp mp *)
EchoTiming[
t0=r^2*idmat+a^2*CosMixMatrix[1,2,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
cosmix=CosMixMatrix[1,1,lmax,mm];
MM0=\[Lambda]hat-a*\[Omega]*sinmix;
Dhvec=Map[Dop,hvec];
D\[Kappa]mat=Map[Dop,\[Kappa]mat];
term1=1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0;
term2=-a/\[Omega]*(r*Dhvec . Bmat0 . sinmix-hvec . Bmat0 . MM0 . cosmix);
tmp=r*idmat+I*a*cosmix;
term3=4*(ones . D\[Kappa]mat) . Bmat0 . MM0 . tmp;
term4=-4*((ones . \[Kappa]mat) . Bmat0 . MM0+I*a*(ones . D\[Kappa]mat) . Bmat0 . sinmix);
\[Rho]hs0lpmp=(term1+term2+term3+term4)/.awsubs;
];

(* lp mm *)
EchoTiming[
t0=r^2*idmat+a^2*CosMixMatrix[-1,2,lmax,mm];
sinmix=SinMixMatrix[0,-1,lmax,mm];
MM0=\[Lambda]hat-a*\[Omega]*sinmix;
Dhvec=Map[Dop,hvec];
D\[Kappa]mat=Map[Dop,\[Kappa]mat];
term1=-1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0;
term2=a/\[Omega]*(r*Dhvec . Bmat0 . sinmix+hvec . Bmat0 . MM0 . CosMixMatrix[-1,1,lmax,mm]);
tmp=r*idmat-I*a*CosMixMatrix[-1,1,lmax,mm];
term3=-4*((ones . D\[Kappa]mat) . Bmat0 . MM0 . tmp);
term4=4*((ones . \[Kappa]mat) . Bmat0 . MM0+I*a*(ones . D\[Kappa]mat) . Bmat0 . sinmix);
\[Rho]chs0lpmm=(term1+term2+term3+term4)/.awsubs;
];

(* Spin 0, 1 and spin 2. *)
(* lp-lp and mp-mp components. *)
EchoTiming[
hs0j0lplp=hs0lplp/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
hs0j1lplp=D[hs0lplp,r]/.hlsubs/.\[Kappa]lsubs/.\[Gamma]llsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
hs0j0mpmp=hs0mpmp/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
hs0j1mpmp=D[hs0mpmp,r]/.hlsubs/.\[Kappa]lsubs/.\[Gamma]llsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
hs1j0lplp=hs1lplp/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs1j1lplp=D[hs1lplp,r]/.Pm1lsubs/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs1j0mpmp=hs1mpmp/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs1j1mpmp=D[hs1mpmp,r]/.Pm1lsubs/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs2j0lplp=hs2lplp/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs2j1lplp=D[hs2lplp,r]/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs2j0mpmp=hs2mpmp/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs2j1mpmp=D[hs2mpmp,r]/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
];

(* lmlm and mmmm components. *)
EchoTiming[
hs0j0lmlm=hs0lmlm/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
hs0j1lmlm=D[hs0lmlm,r]/.hlsubs/.\[Kappa]lsubs/.\[Gamma]llsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
hs0j0mmmm=hs0mmmm/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
hs0j1mmmm=D[hs0mmmm,r]/.hlsubs/.\[Kappa]lsubs/.\[Gamma]llsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
hs1j0lmlm=hs1lmlm/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs1j1lmlm=D[hs1lmlm,r]/.Pm1lsubs/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs1j0mmmm=hs1mmmm/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs1j1mmmm=D[hs1mmmm,r]/.Pm1lsubs/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs2j0lmlm=hs2lmlm/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs2j1lmlm=D[hs2lmlm,r]/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs2j0mmmm=hs2mmmm/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs2j1mmmm=D[hs2mmmm,r]/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
];

(* \[Rho]lpmp and \[Rho]clpmm. *)
EchoTiming[
\[Rho]hs0j0lpmp=(\[Rho]hs0lpmp)/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
\[Rho]hs1j0lpmp=(\[Rho]hs1lpmp)/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
\[Rho]hs2j0lpmp=(\[Rho]hs2lpmp)/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
];

EchoTiming[
\[Rho]chs0j0lpmm=(\[Rho]chs0lpmm)/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
\[Rho]chs1j0lpmm=(\[Rho]chs1lpmm)/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
\[Rho]chs2j0lpmm=(\[Rho]chs2lpmm)/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
];

sourcehlplp=Table[dhlplp[ll],{ll,0,lmax}]/.sourcetermsubs;
sourcehlmlm=Table[dhlmlm[ll],{ll,0,lmax}]/.sourcetermsubs;
sourcehmpmp=Table[dhmpmp[ll],{ll,0,lmax}]/.sourcetermsubs;
sourcehmmmm=Table[dhmmmm[ll],{ll,0,lmax}]/.sourcetermsubs;
source\[Rho]chlpmm=Table[d\[Rho]chlpmm[ll],{ll,0,lmax}]/.sourcetermsubs;
(* OK, let's check out whether the linear system of equations has a solution, in principle. *)
(* We are solving for the unknowns: *)

unknowns=Join[Table[Pm1j0[ll],{ll,lmins1,lmax}],Table[Pm1j1[ll],{ll,lmins1,lmax}],Table[\[Kappa]j0[ll],{ll,lmins0,lmax}],Table[\[Kappa]j1[ll],{ll,lmins0,lmax}]];
unknownstozero=Table[unknowns[[kk]]->0,{kk,1,Length[unknowns]}];

(* We have the system of equations: *)
hj0lplp=hs0j0lplp+hs1j0lplp+hs2j0lplp;
hj0mpmp=hs0j0mpmp+hs1j0mpmp+hs2j0mpmp;
hj1lplp=(hs0j1lplp+hs1j1lplp+hs2j1lplp)-sourcehlplp;
hj1mpmp=(hs0j1mpmp+hs1j1mpmp+hs2j1mpmp)-sourcehmpmp;

hj0lmlm=hs0j0lmlm+hs1j0lmlm+hs2j0lmlm;
hj0mmmm=hs0j0mmmm+hs1j0mmmm+hs2j0mmmm;
hj1lmlm=(hs0j1lmlm+hs1j1lmlm+hs2j1lmlm)-sourcehlmlm; 
hj1mmmm=(hs0j1mmmm+hs1j1mmmm+hs2j1mmmm)-sourcehmmmm;

\[Rho]chj0lpmm=\[Rho]chs0j0lpmm+\[Rho]chs1j0lpmm+\[Rho]chs2j0lpmm; (*unnecessary?*)
\[Rho]hj0lpmp=\[Rho]hs0j0lpmp+\[Rho]hs1j0lpmp+\[Rho]hs2j0lpmp;


eqs=Join[Table[hj0lplp[[ll+1]],{ll,lmins0,lmax}],Table[hj1lplp[[ll+1]],{ll,lmins0,lmax}],Table[hj0mpmp[[ll+1]],{ll,lmins2,lmax}],Table[hj1mpmp[[ll+1]],{ll,lmins2,lmax}]]/.s2jumps/.s0jumps;
brhs=Table[CleanValue[-eqs[[kk]]/.unknownstozero],{kk,1,Length[eqs]}];

(* For Abs[mm] != 1, the number of unknowns should equal the number of equations. However, for Abs[mm]=1, there are two more unknowns than equations,
 and this case must be handled separately, using the spin-1 component of the metric perturbation. *)
If[Abs[mm]==1,
ll=1;
(*eqs=Join[eqs,{\[Rho]chj0lpmm[[ll+1]],\[Rho]chj0lpmm[[ll+2]]}]/.s2jumps/.s0jumps;
brhs=Table[CleanValue[-eqs[[kk]]/.unknownstozero],{kk,1,Length[eqs]}]//N;
*)
eqs=Join[eqs,{hj0lmlm[[ll+1]],hj1lmlm[[ll+1]]}]/.s2jumps/.s0jumps;
brhs=Table[CleanValue[-eqs[[kk]]/.unknownstozero],{kk,1,Length[eqs]}]//N;
];

Mmat=CleanMatrix[Table[Coefficient[eqs[[j]],unknowns[[k]]],{j,1,Length[eqs]},{k,1,Length[eqs]}]];
Det[Mmat];
Quiet[sol=LinearSolve[Mmat,brhs]];
knowns=Table[unknowns[[k]]->CleanValue@sol[[k]],{k,1,Length[unknowns]}];


Print["3. Radial functions and lm modes."];
(* First, let's consider the tortoise coordinate. There are two extra phase factors in the definition of the tortoise coordinate to determine here. I need to make sure that rstar to be consistent with the BHP Toolkit. *)
assum=r>rp&&rp>rm&&rm>0;
integrand=(r^2+rp*rm)/((r-rp)*(r-rm));
t0=Integrate[integrand,r];
t1=r+(rp+rm)/(rp-rm)*(rp*Log[(r-rp)/(rp+rm)]-rm*Log[(r-rm)/(rp+rm)]);
Simplify[t0-t1,Assumptions->assum];
t2=Series[t1,{r,\[Infinity],0}]//Normal;
t3=r+(rp+rm)*Log[r];
Simplify[t2-t3,Assumptions->assum];
t2=(rp+rm)/(rp-rm)*rp*Log[r-rp];
(* OK, so the additional phase factor needed at infinity -2M ln(2M). *)
horphase=Simplify[Simplify[(t1-t2),Assumptions->assum]/.{r->rp},Assumptions->assum];
(* Let z = 1/r. Expand at infinity in powers of z. *)
zrepl={r->1/z};
fnrepl={h0''[r]->z^2*D[z^2*Z0'[z],z],h0'[r]->-z^2*Z0'[z],h0[r]->Z0[z]};
eq=Simplify[(Dop[\[CapitalDelta][r]*Ddag[h0[r]]]-2*I*\[Omega]*r*h0[r]-\[Lambda]0*h0[r])/.\[CapitalDelta]subs/.Ksubs];
zeq=\[CapitalDelta][r]*eq/.\[CapitalDelta]subs/.Ksubs/.fnrepl/.zrepl//Simplify;
prefac=Exp[I*\[Omega]/z+\[Sigma]*Log[z]]*2^(-2*I*M*\[Omega]);
t0=zeq/.{Z0''[z]->D[prefac*Z1[z],{z,2}],Z0'[z]->D[prefac*Z1[z],z],Z0[z]->prefac*Z1[z]};
t1=Simplify[t0/prefac];
Zser=Sum[c[kk]*z^kk,{kk,0,inford}]/.{c[0]->1};
t2=t1/.{Z1''[z]->D[Zser,{z,2}],Z1'[z]->D[Zser,z],Z1[z]->Zser};
ser=Series[t2,{z,0,inford-3}];
\[Sigma]subs=Solve[SeriesCoefficient[ser,-3]==0,\[Sigma]]//First;
upsol=Join[{c[0]->1},\[Sigma]subs];
For[kk=1,kk<=inford,kk++,
t0=SeriesCoefficient[ser,-3+kk]/.upsol;
t1=Solve[t0==0,c[kk]]//First;
AppendTo[upsol,t1//First]
];
h0infser=prefac*Zser/.\[Sigma]subs/.upsol;

(* First seek a series solution for the UP mode of the trace h. *)
h0sersol=prefac*Zser;
eqr2=zeq-1/2*\[CapitalDelta][r]*r^2*h0sersol/.\[CapitalDelta]subs/.{r->1/z};
t0=eqr2/.{Z0''[z]->D[prefac*Z1[z],{z,2}],Z0'[z]->D[prefac*Z1[z],z],Z0[z]->prefac*Z1[z]};
t1=Simplify[t0/prefac/.upsol];
Zser2=1/z*Sum[cr2[kk]*z^kk,{kk,0,inford-1}]; (* Note the extra power of z here *)
t2=t1/.{Z1''[z]->D[Zser2,{z,2}],Z1'[z]->D[Zser2,z],Z1[z]->Zser2};
ser=Series[t2,{z,0,inford-4}];
r2sol={cr2[1]->0};
For[kk=0,kk<=inford-1,kk++,
If[kk==1,
Print[""(*SeriesCoefficient[ser,-4+kk]/.r2sol//Simplify*)];
,
t0=SeriesCoefficient[ser,-4+kk]/.r2sol;
t1=Solve[t0==0,cr2[kk]]//First;
AppendTo[r2sol,t1//First];
];
];
\[Kappa]2infser=prefac*Zser2/.\[Sigma]subs/.r2sol;

(* Now seek a series solution for the inhomogeneous mode sourced by the trace radial function.
Here, insert the series solution for the trace used above. *)
eqr0=zeq-1/2*\[CapitalDelta][r]*h0sersol/.\[CapitalDelta]subs/.{r->1/z};
t0=eqr0/.{Z0''[z]->D[prefac*Z1[z],{z,2}],Z0'[z]->D[prefac*Z1[z],z],Z0[z]->prefac*Z1[z]};
t1=Simplify[t0/prefac/.upsol];
Zser0=z*Sum[cr0[kk]*z^kk,{kk,0,inford-1}]; (* Note the extra powers of z here *)
t2=t1/.{Z1''[z]->D[Zser0,{z,2}],Z1'[z]->D[Zser0,z],Z1[z]->Zser0};
ser=Series[t2,{z,0,inford-3}];
r0sol={};
For[kk=0,kk<=inford-1,kk++,
t0=SeriesCoefficient[ser,kk-2]/.r0sol;
t1=Solve[t0==0,cr0[kk]]//First;
AppendTo[r0sol,t1//First];
];
\[Kappa]0infser=prefac*Zser0/.\[Sigma]subs/.r0sol;
ltry=Max[2,Abs[mm]]; (* For testing *)
\[Lambda]0subs={\[Lambda]0->SetPrecision[eigs0[[ltry+1]],prec]};
odes0=Evaluate@{eq==0,(eq/.{h0''[r]->\[Kappa]r0''[r],h0'[r]->\[Kappa]r0'[r],h0[r]->\[Kappa]r0[r]})==1/2*h0[r],(eq/.{h0''[r]->\[Kappa]r2''[r],h0'[r]->\[Kappa]r2'[r],h0[r]->\[Kappa]r2[r]})==1/2*r^2*h0[r]}/.awsubs;
odes=odes0/.\[Lambda]0subs;
icinf={h0[r]==h0infser,h0'[r]==-z^2*D[h0infser,z],\[Kappa]r0[r]==\[Kappa]0infser,\[Kappa]r0'[r]==-z^2*D[\[Kappa]0infser,z],\[Kappa]r2[r]==\[Kappa]2infser,\[Kappa]r2'[r]==-z^2*D[\[Kappa]2infser,z]}/.{r->rinf,z->1/rinf}/.awsubs//Simplify;
ics=icinf/.\[Lambda]0subs;
odeprec=SetPrecision[Join[odes,ics],prec];
dsol=NDSolve[odeprec,{h0[r],\[Kappa]r0[r],\[Kappa]r2[r],h0'[r],\[Kappa]r0'[r],\[Kappa]r2'[r]},{r,r0,rinf}, AccuracyGoal->accgoal, PrecisionGoal->accgoal,WorkingPrecision->prec,MaxSteps->Infinity]//First;

(* Check that my solution for the UP mode of the trace agrees with the BHP Toolkit version. *)
(* tmp=TeukolskyRadial@@({0,ltry,mm,a,\[Omega]}/.awsubs); *)
tmp=TeukolskyRadial[0,ltry,mm,a/.awsubs,\[Omega]/.awsubs,Method->{"NumericalIntegration","Domain"->{"In"->{rmin,r0},"Up"->{r0,rmax}}}];
(*tmp=TeukolskyRadial[0,ltry,mm,a/.awsubs,\[Omega]/.awsubs,Method->"MST"]*)
h0up=tmp["Up"];
(* Let's compare values at r0. *)
rtry=SetPrecision[r0+RandomReal[],prec];
t0={h0[r]/.dsol/.{r->rtry}/.awsubs,h0up[rtry]};
tmp=t0[[1]]-t0[[2]];
{tmp,Abs[tmp]};
(*
LogPlot[Abs[((h0[r]/.dsol)-h0up[rtry])/.{r->rtry}/.awsubs],{rtry,r0,SetPrecision[9,prec]}]
Plot[{Re[h0up[rtry]],Im[h0up[rtry]]},{rtry,r0,SetPrecision[9,prec]}]
*)

(* Now construct all the radial functions that will be needed, for all ell. *)
\[Kappa]upsol={};
vars={h0[r],\[Kappa]r0[r],\[Kappa]r2[r],h0'[r],\[Kappa]r0'[r],\[Kappa]r2'[r]};
Monitor[
For[ll=0,ll<=lmax,ll++,
val={0,0};
If[ll>=lmins0,
odesys=SetPrecision[Join[odes0,icinf]/.{\[Lambda]0->SetPrecision[eigs0[[ll+1]],prec]},prec];
dsol=NDSolve[odesys,vars,{r,rupmin,rupmax},AccuracyGoal->accgoal,PrecisionGoal->accgoal,WorkingPrecision->prec,MaxSteps->Infinity]//First;
val=vars/.dsol;
];
AppendTo[\[Kappa]upsol,val];
];,ll];
Clear[ll];

(* Check that the series solution for the h function is in good agreement with the numerical solution from the ODE, in the appropriate large-r domain. *)

(* Now do all of this again, but for the IN mode. *)
(* First, find a Frobenius series expansion of the trace h *)
eq=Simplify[Dop[\[CapitalDelta][r]*Ddag[h0[r]]]-2*I*\[Omega]*r*h0[r]-\[Lambda]0*h0[r]]/.\[CapitalDelta]subs/.Ksubs;
aMsubs={M->(rp+rm)/2,a->Sqrt[rp*rm]};
horfac=x^\[Nu]*Exp[-I*(\[Omega]-m*Sqrt[rp*rm]/((rp+rm)*rp))*horphase];
horser=Sum[ch[kk]*x^(kk),{kk,0,horord}];
tmp=horfac*horser;
t0=(\[CapitalDelta][r]*eq)/.\[CapitalDelta]subs/.{h0''[r]->D[tmp,{x,2}],h0'[r]->D[tmp,x],h0[r]->tmp}/.aMsubs/.{r->x+rp};
ser=Series[Simplify[t0/horfac],{x,0,horord},Assumptions->assum]//Simplify;
t2=SeriesCoefficient[ser,0];
\[Nu]val=-I*(\[Omega]-m*a/((rp+rm)*rp))*rp*(rp+rm)/(rp-rm)/.aMsubs;
\[Nu]subs={\[Nu]->\[Nu]val};
t3=Simplify[t2/.\[Nu]subs,Assumptions->assum];
horco={ch[0]->1};
rprmvals={rp->1+Sqrt[1-a^2],rm->1-Sqrt[1-a^2]}/.awsubs;
For[kk=1,kk<=horord,kk++,
t0=SeriesCoefficient[ser,kk];
t1=Simplify[First@Solve[t0==0,ch[kk]]/.horco/.\[Nu]subs/.rprmvals/.awsubs,Assumptions->assum];
AppendTo[horco,First@t1];
];
horsol=horfac*horser/.horco/.\[Nu]subs/.rprmvals/.awsubs;
eqr2=eq-1/2*r^2*horfac*horser;
r2ser=x*Sum[c\[Kappa]2[kk]*x^(kk),{kk,0,horord}];
tmp=horfac*r2ser;
t0=(\[CapitalDelta][r]*eqr2)/.\[CapitalDelta]subs/.{h0''[r]->D[tmp,{x,2}],h0'[r]->D[tmp,x],h0[r]->tmp}/.aMsubs/.{r->x+rp};
ser=Series[Simplify[t0/horfac]/.horco/.\[Nu]subs/.rprmvals/.awsubs,{x,0,horord+1},Assumptions->assum]//Simplify;
horr2co={};
For[kk=0,kk<=horord,kk++,
t0=SeriesCoefficient[ser,kk+1];
t1=Simplify[First@Solve[t0==0,c\[Kappa]2[kk]]/.horr2co/.\[Nu]subs/.rprmvals/.awsubs,Assumptions->assum];
AppendTo[horr2co,First@t1];
];
horr2sol=horfac*r2ser/.horr2co/.\[Nu]subs/.rprmvals/.awsubs;
eqr0=eq-1/2*horfac*horser;
r0ser=x*Sum[c\[Kappa]0[kk]*x^(kk),{kk,0,horord}];
tmp=horfac*r0ser;
t0=(\[CapitalDelta][r]*eqr0)/.\[CapitalDelta]subs/.{h0''[r]->D[tmp,{x,2}],h0'[r]->D[tmp,x],h0[r]->tmp}/.aMsubs/.{r->x+rp};
ser=Series[Simplify[t0/horfac]/.horco/.\[Nu]subs/.rprmvals/.awsubs,{x,0,horord+1},Assumptions->assum]//Simplify;
horr0co={};
For[kk=0,kk<=horord,kk++,
t0=SeriesCoefficient[ser,kk+1];
t1=Simplify[First@Solve[t0==0,c\[Kappa]0[kk]]/.horr0co/.\[Nu]subs/.rprmvals/.awsubs,Assumptions->assum];
AppendTo[horr0co,First@t1];
];
horr0sol=horfac*r0ser/.horr0co/.\[Nu]subs/.rprmvals/.awsubs;
ltry=Max[2,Abs[mm]];
rplot=SetPrecision[100,prec];
\[Lambda]0subs={\[Lambda]0->SetPrecision[eigs0[[ltry+1]],prec]};
odes=odes0/.\[Lambda]0subs;
ichor={h0[rhor]==horsol,h0'[rhor]==D[horsol,x],\[Kappa]r0[rhor]==horr0sol,\[Kappa]r0'[rhor]==D[horr0sol,x],\[Kappa]r2[rhor]==horr2sol,\[Kappa]r2'[rhor]==D[horr2sol,x]}/.{x->xhor};
ics=ichor/.\[Lambda]0subs;
odeprec=SetPrecision[Join[odes,ics],prec];
dsol=NDSolve[odeprec,{h0[r],\[Kappa]r0[r],\[Kappa]r2[r]},{r,rhor,r0},AccuracyGoal->accgoal,PrecisionGoal->accgoal,WorkingPrecision->prec]//First;

(* Check that my solution for the IN mode of the trace agrees with the BHP Toolkit version. *)
tmp=TeukolskyRadial@@({0,ltry,mm,a,\[Omega]}/.awsubs);
h0in=tmp["In"];
rtry=SetPrecision[r0-RandomReal[],prec];
t0={h0[r]/.dsol/.{r->rtry}/.awsubs,h0in[rtry/.awsubs]};
tmp=t0[[1]]-t0[[2]];
{t0[[1]],tmp,Abs[tmp],Abs[tmp/t0[[1]]]};

(* Now construct all the radial functions that will be needed. *)
\[Kappa]insol={};
vars={h0[r],\[Kappa]r0[r],\[Kappa]r2[r],h0'[r],\[Kappa]r0'[r],\[Kappa]r2'[r]};
Monitor[
For[ll=0,ll<=lmax,ll++,
val={0,0};
If[ll>=lmins0,
odesys=SetPrecision[Join[odes0,ichor]/.{\[Lambda]0->SetPrecision[eigs0[[ll+1]],prec]},prec];
dsol=NDSolve[odesys,vars,{r,rmin,r0},AccuracyGoal->accgoal,PrecisionGoal->accgoal,WorkingPrecision->prec]//First;
val=vars/.dsol;
];
AppendTo[\[Kappa]insol,val];
];,ll];
Clear[ll];

(* 
  In this section I derive the coefficients \[Alpha] for the radial functions in the form
 R = \[Alpha]up Rup Theta[r-r0] + \[Alpha]in Rin Theta[r0-r]
 from the jumps in R and its derivative.
*)
\[Alpha]s2cos={};
Monitor[
For[ll=0,ll<=lmax,ll++,
val={0,0};
If[ll>=lmins2,
t0=TeukolskyRadial@@({-2,ll,mm,a,\[Omega]}/.awsubs);
Rup=t0["Up"];
Rin=t0["In"];
eq0=\[Alpha]up*Rup[r0]-\[Alpha]in*Rin[r0]==Pm2j0[ll]/.s2jumps;
eq1=\[Alpha]up*Rup'[r0]-\[Alpha]in*Rin'[r0]==Pm2j1[ll]/.s2jumps;
sol=Solve[{eq0,eq1},{\[Alpha]up,\[Alpha]in},WorkingPrecision->prec]//First;
val={\[Alpha]up,\[Alpha]in}/.sol
];
AppendTo[\[Alpha]s2cos,val];
];,ll];
Clear[ll];

\[Alpha]s1cos={};
Monitor[
For[ll=0,ll<=lmax,ll++,
val={0,0};
If[ll>=lmins1,
t0=TeukolskyRadial@@{-1,ll,mm,a,\[Omega]}/.awsubs;
Rup=t0["Up"];
Rin=t0["In"];
eq0=\[Alpha]up*Rup[r0]-\[Alpha]in*Rin[r0]==Pm1j0[ll]/.knowns;
eq1=\[Alpha]up*Rup'[r0]-\[Alpha]in*Rin'[r0]==Pm1j1[ll]/.knowns;
sol=Solve[{eq0,eq1},{\[Alpha]up,\[Alpha]in},WorkingPrecision->prec]//First;
val={\[Alpha]up,\[Alpha]in}/.sol
];
AppendTo[\[Alpha]s1cos,val];
];,ll];
Clear[ll];

\[Alpha]s0hcos={};
Monitor[
For[ll=0,ll<=lmax,ll++,
val={0,0};
If[ll>=lmins0,
t0=TeukolskyRadial@@{0,ll,mm,a,\[Omega]}/.awsubs;
Rup=t0["Up"];
Rin=t0["In"];
eq0=\[Alpha]up*Rup[r0]-\[Alpha]in*Rin[r0]==0/.s0jumps;
eq1=\[Alpha]up*Rup'[r0]-\[Alpha]in*Rin'[r0]==hj1[ll]/.s0jumps;
sol=Solve[{eq0,eq1},{\[Alpha]up,\[Alpha]in},WorkingPrecision->prec]//First;
val={\[Alpha]up,\[Alpha]in}/.sol
];
AppendTo[\[Alpha]s0hcos,val];
];,ll];
Clear[ll];

(* Important: The amplitude of the \[Kappa] part is set by the trace (h) determined above. 
Homogeneous modes are then added to satisfy the jump conditions. *)
\[Alpha]s0\[Kappa]cos={};
Monitor[
For[ll=0,ll<=lmax,ll++,
val={0,0};
If[ll>=lmins0,
Rupfn=\[Alpha]up*\[Kappa]upsol[[ll+1]][[1]]+\[Alpha]s0hcos[[ll+1]][[1]]*(\[Kappa]upsol[[ll+1]][[3]]+a^2*\[Gamma]ll[ll]*\[Kappa]upsol[[ll+1]][[2]])/.\[Gamma]llsubs/.awsubs;
Rup=Rupfn/.{r->r0};
dRupfn=\[Alpha]up*\[Kappa]upsol[[ll+1]][[4]]+\[Alpha]s0hcos[[ll+1]][[1]]*(\[Kappa]upsol[[ll+1]][[6]]+a^2*\[Gamma]ll[ll]*\[Kappa]upsol[[ll+1]][[5]])/.\[Gamma]llsubs/.awsubs;
dRup=dRupfn/.{r->r0};

Rinfn=\[Alpha]in*\[Kappa]insol[[ll+1]][[1]]+\[Alpha]s0hcos[[ll+1]][[2]]*(\[Kappa]insol[[ll+1]][[3]]+a^2*\[Gamma]ll[ll]*\[Kappa]insol[[ll+1]][[2]])/.\[Gamma]llsubs/.awsubs;
Rin=Rinfn/.{r->r0};
dRinfn=\[Alpha]in*\[Kappa]insol[[ll+1]][[4]]+\[Alpha]s0hcos[[ll+1]][[2]]*(\[Kappa]insol[[ll+1]][[6]]+a^2*\[Gamma]ll[ll]*\[Kappa]insol[[ll+1]][[5]])/.\[Gamma]llsubs/.awsubs;
dRin=dRinfn/.{r->r0};

eq0=Rup-Rin==\[Kappa]j0[ll]/.knowns;
eq1=dRup-dRin==\[Kappa]j1[ll]/.knowns;
sol=Solve[{eq0,eq1},{\[Alpha]up,\[Alpha]in},WorkingPrecision->prec]//First;
val={\[Alpha]up,\[Alpha]in}/.sol
];
AppendTo[\[Alpha]s0\[Kappa]cos,val];
];,ll];
Clear[ll];


Print["4. Prepare grid and mode functions."];


TeukMethod={"NumericalIntegration","Domain"->{"In"->{rmin,r0},"Up"->{r0,rmax}}};
(* TeukMethod={"MST","RenormalizedAngularMomentum"\[Rule]"Monodromy"} *)

EchoTiming[
Pm2uprepl={};Pm2inrepl={};
For[ll=lmins2,ll<=lmax,ll++,
t0=TeukolskyRadial@@({-2,ll,mm,a,\[Omega],Method->TeukMethod}/.awsubs);
t1=t0["Up"];
AppendTo[Pm2uprepl,{Pm2[ll][r]->\[Alpha]up*t1[r],Pm2[ll]'[r]->\[Alpha]up*t1'[r]}/.{\[Alpha]up->\[Alpha]s2cos[[ll+1]][[1]]}];
t1=t0["In"];
AppendTo[Pm2inrepl,{Pm2[ll][r]->\[Alpha]in*t1[r],Pm2[ll]'[r]->\[Alpha]in*t1'[r]}/.{\[Alpha]in->\[Alpha]s2cos[[ll+1]][[2]]}];
];
Pm2uprepl=Pm2uprepl//Flatten;
Pm2inrepl=Pm2inrepl//Flatten;
];

EchoTiming[
Pm1uprepl={};Pm1inrepl={};
For[ll=lmins1,ll<=lmax,ll++,
t0=TeukolskyRadial@@({-1,ll,mm,a,\[Omega],Method->TeukMethod}/.awsubs);
t1=t0["Up"];
AppendTo[Pm1uprepl,{Pm1[ll][r]->\[Alpha]up*t1[r],Pm1[ll]'[r]->\[Alpha]up*t1'[r]}/.{\[Alpha]up->\[Alpha]s1cos[[ll+1]][[1]]}];
t1=t0["In"];
AppendTo[Pm1inrepl,{Pm1[ll][r]->\[Alpha]in*t1[r],Pm1[ll]'[r]->\[Alpha]in*t1'[r]}/.{\[Alpha]in->\[Alpha]s1cos[[ll+1]][[2]]}];
];
Pm1uprepl=Pm1uprepl//Flatten;
Pm1inrepl=Pm1inrepl//Flatten;
];

EchoTiming[
h0uprepl={};h0inrepl={};
For[ll=lmins0,ll<=lmax,ll++,
t0=TeukolskyRadial@@({0,ll,mm,a,\[Omega],Method->TeukMethod}/.awsubs);
t1=t0["Up"];
AppendTo[h0uprepl,{h[ll][r]->\[Alpha]up*t1[r],h[ll]'[r]->\[Alpha]up*t1'[r]}/.{\[Alpha]up->\[Alpha]s0hcos[[ll+1]][[1]]}];
t1=t0["In"];
AppendTo[h0inrepl,{h[ll][r]->\[Alpha]in*t1[r],h[ll]'[r]->\[Alpha]in*t1'[r]}/.{\[Alpha]in->\[Alpha]s0hcos[[ll+1]][[2]]}];
];
h0uprepl=h0uprepl//Flatten;
h0inrepl=h0inrepl//Flatten;
];

(* Now for the \[Kappa] functions. *)
(* Important: The amplitude of the \[Kappa] part is set by the trace (h) determined above. 
Homogeneous modes have been added to satisfy the jump conditions. *)
EchoTiming[
\[Kappa]0uprepl={};\[Kappa]0inrepl={};
For[ll=lmins0,ll<=lmax,ll++,
tmp=\[Alpha]up*\[Kappa]upsol[[ll+1]][[1]]+\[Alpha]s0hcos[[ll+1]][[1]]*(\[Kappa]upsol[[ll+1]][[3]]+a^2*\[Gamma]ll[ll]*\[Kappa]upsol[[ll+1]][[2]]);
t1=(tmp/.\[Gamma]llsubs/.{\[Alpha]up->\[Alpha]s0\[Kappa]cos[[ll+1]][[1]]}/.awsubs);
tmp=\[Alpha]up*\[Kappa]upsol[[ll+1]][[4]]+\[Alpha]s0hcos[[ll+1]][[1]]*(\[Kappa]upsol[[ll+1]][[6]]+a^2*\[Gamma]ll[ll]*\[Kappa]upsol[[ll+1]][[5]]);
t2=(tmp/.\[Gamma]llsubs/.{\[Alpha]up->\[Alpha]s0\[Kappa]cos[[ll+1]][[1]]}/.awsubs);
AppendTo[\[Kappa]0uprepl,{\[Kappa][ll][r]->t1,\[Kappa][ll]'[r]->t2}];

tmp=\[Alpha]in*\[Kappa]insol[[ll+1]][[1]]+\[Alpha]s0hcos[[ll+1]][[2]]*(\[Kappa]insol[[ll+1]][[3]]+a^2*\[Gamma]ll[ll]*\[Kappa]insol[[ll+1]][[2]]);
t1=(tmp/.\[Gamma]llsubs/.{\[Alpha]in->\[Alpha]s0\[Kappa]cos[[ll+1]][[2]]}/.awsubs);
tmp=\[Alpha]in*\[Kappa]insol[[ll+1]][[4]]+\[Alpha]s0hcos[[ll+1]][[2]]*(\[Kappa]insol[[ll+1]][[6]]+a^2*\[Gamma]ll[ll]*\[Kappa]insol[[ll+1]][[5]]);
t2=(tmp/.\[Gamma]llsubs/.{\[Alpha]in->\[Alpha]s0\[Kappa]cos[[ll+1]][[2]]}/.awsubs);
AppendTo[\[Kappa]0inrepl,{\[Kappa][ll][r]->t1,\[Kappa][ll]'[r]->t2}];
];
\[Kappa]0uprepl=\[Kappa]0uprepl//Flatten;
\[Kappa]0inrepl=\[Kappa]0inrepl//Flatten;
];

GetPupsubs[ltry_,rtry_?NumericQ]:=
Join[
EvaluateRHS[(HeadFn[Pm2,ltry][[1;;2]]/.Pm2uprepl),rtry],
EvaluateRHS[(HeadFn[Pm1,ltry][[1;;2]]/.Pm1uprepl),rtry]
];

GetPinsubs[ltry_,rtry_?NumericQ]:=
Join[
EvaluateRHS[(HeadFn[Pm2,ltry][[1;;2]]/.Pm2inrepl),rtry],
EvaluateRHS[(HeadFn[Pm1,ltry][[1;;2]]/.Pm1inrepl),rtry]
];

GetPsubs[ltry_,rtry_?NumericQ]:=If[Greater[rtry-r0,0],GetPupsubs[ltry,rtry],GetPinsubs[ltry,rtry]];
EchoTiming[
GetPsubs[4,SetPrecision[5.0,prec]];
GetPsubs[4,SetPrecision[7.0,prec]];
];

(* This is required for filling in the 1D grid. (Bad coding style though). *)
(* HeadFnAlt[head_,ll_]:=Table[D[head[ll][r]->head[ll][r],{r,kk}],{kk,0,3}]; *)

GetPupsubsAlt[spin_,ltry_,rtry_?NumericQ]:=
If[spin==2,EvaluateRHS[Pm2uprepl,rtry],EvaluateRHS[Pm1uprepl,rtry]];

GetPinsubsAlt[spin_,ltry_,rtry_?NumericQ]:=
If[spin==2,EvaluateRHS[Pm2inrepl,rtry],EvaluateRHS[Pm1inrepl,rtry]];

GetPsubsAlt[spin_,ltry_,rtry_?NumericQ]:=
If[Greater[rtry-r0,0],GetPupsubsAlt[spin,ltry,rtry],GetPinsubsAlt[spin,ltry,rtry]];

GetPsubsAll[spin_,rtry_?NumericQ]:=Module[{},
subs1=Table[GetPsubsAlt[spin,ll,rtry],{ll,Abs[mm],lmax}]//Flatten;
\[CapitalDelta]Ktbl=Map[#[[1]]->(#[[2]]/.awsubs/.{r->rtry})&,\[CapitalDelta]Ksubs];
Ppluslrepl=If[spin==2,Pp2lrepl,Pp1lrepl];
Pminuslsubs=If[spin==2,Pm2lsubs,Pm1lsubs];
Ppluslsubs=If[spin==2,Pp2lsubs,Pp1lsubs];
subs2=EvaluateRHS[Ppluslrepl/.awsubs/.\[CapitalDelta]Ktbl/.subs1,rtry];
subs3=EvaluateRHS[Pminuslsubs/.awsubs/.\[CapitalDelta]Ktbl/.subs2/.subs1,rtry];
subs4=EvaluateRHS[Ppluslsubs/.subs2/.awsubs/.\[CapitalDelta]Ktbl,rtry];
Join[subs1,subs2,subs3,subs4]
];

EchoTiming[
GetPsubsAll[1,SetPrecision[7.1,prec]];
];

(* Scalar part needs careful handling. *)
GetSpin0subs[rtry_?NumericQ]:=If[Greater[rtry-r0,0],
Join[EvaluateRHS[h0uprepl,rtry],EvaluateRHS[\[Kappa]0uprepl,rtry]],
Join[EvaluateRHS[h0inrepl,rtry],EvaluateRHS[\[Kappa]0inrepl,rtry]]
];

Print["A. Fill 1D grid."];
(*  Calcs1components ??? *)
idmat=DiagonalMatrix[Table[If[ll>=lmins0,1,0],{ll,0,lmax}]];
t0=-(Bmatm1+signmat . Bmatp1) . \[Lambda]hat;
t1=a\[Omega]0*(Bmatm1 . SinMixMatrix[-1,0,lmax,mm]+signmat . Bmatp1 . SinMixMatrix[1,0,lmax,mm]);
SS0=t0+t1;
Ss1lplp=SS0;
Ss1lmlm=SS0;
Ss1mpmp=Bmatp1 . (-\[Lambda]2hat+a\[Omega]0*SinMixMatrix[1,2,lmax,mm]);
Ss1mmmm=Bmatm1 . (\[Lambda]2hat-a\[Omega]0*SinMixMatrix[-1,-2,lmax,mm]);
(* Spin 1: \[Rho]hlpmp *)
tmp1=signmat . Table[Ddag[Pp1vec[[ll+1]]],{ll,0,lmax}];
tmp2=Table[Dop[Pm1vec[[ll+1]]],{ll,0,lmax}];
calPs1=tmp2+tmp1;
DcalPs1=Map[Dop,calPs1];
(* Spin 1: l+ m+ *)
cosmix=CosMixMatrix[1,1,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
line1=(signmat . (r*DcalPs1-2*calPs1)) . Bmatp1+I*a*(signmat . DcalPs1) . Bmatp1 . cosmix;
tmp1=1/\[CapitalDelta][r]*signmat . Pp1vec;
tmp2=-r*\[Lambda]hat-I*a*\[Lambda]hat . cosmix+(a*\[Omega]*r+2*I*a)*sinmix+I*a^2*\[Omega]*sinmix . cosmix;
line2=tmp1 . SS0 . tmp2;
\[Rho]hraws1lpmp=(line1+line2);
(* Spin 1: l+ m- term *)
line1=-(r*DcalPs1-2*calPs1) . Bmatm1+I*a*DcalPs1 . Bmatm1 . CosMixMatrix[-1,1,lmax,mm];
tmp1=1/\[CapitalDelta][r]*signmat . Pp1vec;
tmp2=r*\[Lambda]hat-I*a*\[Lambda]hat . CosMixMatrix[-1,1,lmax,mm]-(a*\[Omega]*r+2*I*a)*SinMixMatrix[0,-1,lmax,mm]+I*a^2*\[Omega]*SinMixMatrix[0,-1,lmax,mm] . CosMixMatrix[-1,1,lmax,mm];
line2=tmp1 . SS0 . tmp2;
\[Rho]chraws1lpmm=(line1+line2);
(* Spin 1: l- m+ *)
DdagcalPs1=Map[Ddag,calPs1];
cosmix=CosMixMatrix[1,1,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
tmp1=signmat . (r*DdagcalPs1-2*calPs1);
tmp2=signmat . DdagcalPs1;
line1 =tmp1 . Bmatp1-I*a*tmp2 . Bmatp1 . cosmix;
tmp1=1/\[CapitalDelta][r]*Pm1vec;
tmp2=-r*\[Lambda]hat+I*a*\[Lambda]hat . cosmix+(a*\[Omega]*r-2*I*a)*sinmix-I*a^2*\[Omega]*sinmix . cosmix;
line2=tmp1 . SS0 . tmp2;
\[Rho]chraws1lmmp=(line1+line2);
(* Spin 1: l- m- *)
cosmix=CosMixMatrix[-1,1,lmax,mm];
sinmix=SinMixMatrix[0,-1,lmax,mm];
line1=-(r*DdagcalPs1-2*calPs1) . Bmatm1-I*a*DdagcalPs1 . Bmatm1 . cosmix;
tmp1=1/\[CapitalDelta][r]*Pm1vec;
tmp2=r*\[Lambda]hat+I*a*\[Lambda]hat . cosmix-(a*\[Omega]*r-2*I*a)*sinmix-I*a^2*\[Omega]*sinmix . cosmix;
line2=tmp1 . SS0 . tmp2;
\[Rho]hraws1lmmm=(line1+line2);
(* Spin 1: l+ l- *)
cosmix=CosMixMatrix[0,1,lmax,mm];
\[CapitalSigma]mat=r^2*idmat+a^2*CosMixMatrix[0,2,lmax,mm];
term1=calPs1 . SS0 . \[CapitalSigma]mat;
term2=-2*r*(Pm1vec+signmat . Pp1vec) . SS0;
term3=2*a^2*calPs1 . (Bmatm1 . SinMixMatrix[-1,0,lmax,mm]-signmat . Bmatp1 . SinMixMatrix[1,0,lmax,mm]) . cosmix;
\[CapitalSigma]\[CapitalDelta]hraws1lplm=(term1+term2+term3);
idmat=DiagonalMatrix[Table[If[ll>=lmins0,1,0],{ll,0,lmax}]];
t0=Bmatp2 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]2hat . SinMixMatrix[1,0,lmax,mm]+a\[Omega]0^2*SinMixMatrix[2,0,lmax,mm]);
t1=signmat . Bmatm2 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]2hat . SinMixMatrix[-1,0,lmax,mm]+a\[Omega]0^2*SinMixMatrix[-2,0,lmax,mm]);
Ss2lplp=(t0+t1)/.awsubs;
Ss2lmlm=signmat . Ss2lplp/.awsubs; 
Ss2mpmp=Bmatp2/.awsubs;
Ss2mmmm=signmat . Bmatm2/.awsubs;
SSplus=(Bmatp2 . (\[CapitalLambda]hat-2*a*\[Omega]*\[Lambda]2hat . SinMixMatrix[1,0,lmax,mm]+a^2*\[Omega]^2*SinMixMatrix[2,0,lmax,mm]))/.awsubs;
SSminus=(Bmatm2 . (\[CapitalLambda]hat-2*a*\[Omega]*\[Lambda]2hat . SinMixMatrix[-1,0,lmax,mm]+a^2*\[Omega]^2*SinMixMatrix[-2,0,lmax,mm]))/.awsubs;
DdagPp2vec=Map[Ddag,Pp2vec]/.\[CapitalDelta]Ksimp;
DdagDdagPp2vec=Map[Ddag,DdagPp2vec]/.\[CapitalDelta]Ksimp;
DdagDdagDdagPp2vec=Map[Ddag,DdagDdagPp2vec]/.\[CapitalDelta]Ksimp;
DPm2vec=Map[Dop,Pm2vec]/.\[CapitalDelta]Ksimp;
DDPm2vec=Map[Dop,DPm2vec]/.\[CapitalDelta]Ksimp;
DDDPm2vec=Map[Dop,DDPm2vec]/.\[CapitalDelta]Ksimp;
DopDdagDdagPp2vec=Map[Dop,DdagDdagPp2vec]/.\[CapitalDelta]Ksimp;
DdagDDPm2vec=Map[Ddag,DDPm2vec]/.\[CapitalDelta]Ksimp;
(* lp mp *)
cosmix=CosMixMatrix[1,1,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
sinmix1m1=SinMixMatrix[-1,1,lmax,mm];
\[Rho]mat=r*idmat+I*a*cosmix;
MMp2=Bmatp2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[2,1,lmax,mm]);
MMm2=Bmatm2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[-2,-1,lmax,mm]);
term1=1/(12*\[Omega]^2)*(DDPm2vec . signmat . SSplus-DdagDdagPp2vec . signmat . SSminus) . (\[Lambda]hat-a*\[Omega]*sinmix);
term2=-1/(12*\[Omega]^2)*(DDDPm2vec . signmat . SSplus-DopDdagDdagPp2vec . signmat . SSminus) . ((\[Lambda]hat-a*\[Omega]*sinmix) . \[Rho]mat-I*a*sinmix);
term3=-1/(3*I*\[Omega])*(DDDPm2vec . signmat . MMp2 . \[Rho]mat . \[Rho]mat-2*DDPm2vec . signmat . MMp2 . \[Rho]mat+2*DPm2vec . signmat . MMp2);
tmp=(\[Lambda]hat . \[Lambda]hat-2*a*\[Omega]*\[Lambda]hat . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]mat . \[Rho]mat-2*I*a*(\[Lambda]hat . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]mat-2*a^2*sinmix1m1;
term4=-1/(3*I*\[Omega]*\[CapitalDelta][r])*DdagPp2vec . (signmat . MMm2 . tmp);
(* Don't forget the \[Beta] factor. *)
tmp=term1+term2+term3+term4;
\[Rho]hraws2lpmp=(tmp/(-6*I*M*\[Omega]));
(* lp mm *)
\[Rho]cmat=r*idmat-I*a*CosMixMatrix[-1,1,lmax,mm];
sinmix=SinMixMatrix[0,-1,lmax,mm];
MMp2=Bmatp2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[2,1,lmax,mm]);
MMm2=Bmatm2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[-2,-1,lmax,mm]);
term1=-1/(12*\[Omega]^2)*(DDPm2vec . SSminus-DdagDdagPp2vec . SSplus) . (\[Lambda]hat-a*\[Omega]*sinmix);
term2=1/(12*\[Omega]^2)*(DDDPm2vec . SSminus-DopDdagDdagPp2vec . SSplus) . ((\[Lambda]hat-a*\[Omega]*sinmix) . \[Rho]cmat-I*a*sinmix);
term3=1/(3*I*\[Omega])*(DDDPm2vec . MMm2 . \[Rho]cmat . \[Rho]cmat-2*DDPm2vec . MMm2 . \[Rho]cmat+2*DPm2vec . MMm2);
sinmix1m1=SinMixMatrix[1,-1,lmax,mm];
tmp=(\[Lambda]hat . \[Lambda]hat-2*a*\[Omega]*\[Lambda]hat . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]cmat . \[Rho]cmat-2*I*a*(\[Lambda]hat . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]cmat-2*a^2*sinmix1m1;
term4=1/(3*I*\[Omega]*\[CapitalDelta][r])*DdagPp2vec . MMp2 . tmp;
tmp=term1+term2+term3+term4;
\[Rho]chraws2lpmm=(tmp/(-6*I*M*\[Omega]));
(* l- m+ *)
sinmix=SinMixMatrix[0,1,lmax,mm];
cosmix=CosMixMatrix[1,1,lmax,mm];
sinmix1m1=SinMixMatrix[-1,1,lmax,mm];
\[Rho]cmat=r*idmat-I*a*cosmix;
MMp2=Bmatp2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[2,1,lmax,mm]);
MMm2=Bmatm2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[-2,-1,lmax,mm]);
term1=1/(12*\[Omega]^2)*(DDPm2vec . SSminus-DdagDdagPp2vec . SSplus) . (\[Lambda]hat-a*\[Omega]*sinmix);
term2=-1/(12*\[Omega]^2)*(DdagDDPm2vec . SSminus-DdagDdagDdagPp2vec . SSplus) . ((\[Lambda]hat-a*\[Omega]*sinmix) . \[Rho]cmat+I*a*sinmix);
term3=-1/(3*I*\[Omega])*(DdagDdagDdagPp2vec . MMp2 . \[Rho]cmat . \[Rho]cmat-2*DdagDdagPp2vec . MMp2 . \[Rho]cmat+2*DdagPp2vec . MMp2);
tmp=(\[Lambda]hat . \[Lambda]hat-2*a*\[Omega]*\[Lambda]hat . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]cmat . \[Rho]cmat+2*I*a*(\[Lambda]hat . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]cmat-2*a^2*sinmix1m1;
term4=-1/(3*I*\[Omega]*\[CapitalDelta][r])*DPm2vec . MMm2 . tmp;
tmp=term1+term2+term3+term4;
\[Rho]chraws2lmmp=(tmp/(-6*I*M*\[Omega]));
(* l- m- *)
sinmix=SinMixMatrix[0,-1,lmax,mm];
cosmix=CosMixMatrix[-1,1,lmax,mm];
sinmix1m1=SinMixMatrix[1,-1,lmax,mm];
\[Rho]mat=r*idmat+I*a*cosmix;
MMp2=Bmatp2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[2,1,lmax,mm]);
MMm2=Bmatm2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[-2,-1,lmax,mm]);
term1=-1/(12*\[Omega]^2)*(DDPm2vec . signmat . SSplus-DdagDdagPp2vec . signmat . SSminus) . (\[Lambda]hat-a*\[Omega]*sinmix);
term2=1/(12*\[Omega]^2)*(DdagDDPm2vec . signmat . SSplus-DdagDdagDdagPp2vec . signmat . SSminus) . ((\[Lambda]hat-a*\[Omega]*sinmix) . \[Rho]mat+I*a*sinmix);
term3=1/(3*I*\[Omega])*(DdagDdagDdagPp2vec . signmat . MMm2 . \[Rho]mat . \[Rho]mat-2*DdagDdagPp2vec . signmat . MMm2 . \[Rho]mat+2*DdagPp2vec . signmat . MMm2);
tmp=(\[Lambda]hat . \[Lambda]hat-2*a*\[Omega]*\[Lambda]hat . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]mat . \[Rho]mat+2*I*a*(\[Lambda]hat . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]mat-2*a^2*sinmix1m1;
term4=1/(3*I*\[Omega]*\[CapitalDelta][r])*DPm2vec . (signmat . MMp2 . tmp);
(* Don't forget the \[Beta] factor. *)
tmp=term1+term2+term3+term4;
\[Rho]hraws2lmmm=(tmp/(-6*I*M*\[Omega]));
(* l+ l- *)
cosmix=CosMixMatrix[0,1,lmax,mm];
\[CapitalSigma]mat=r^2*idmat+a^2*CosMixMatrix[0,2,lmax,mm];
\[Rho]cmat=r*idmat-I*a*cosmix;
\[Rho]mat=r*idmat+I*a*cosmix;
tmp1=(Map[Dop,\[CapitalDelta][r]*Map[Ddag,DDPm2vec]]+Map[Ddag,\[CapitalDelta][r]*Map[Dop,DDPm2vec]]);
tmp2=Map[D[#,r]&,DDPm2vec];
tmp3=\[Lambda]hat . (SinMixMatrix[-1,0,lmax,mm]-SinMixMatrix[1,0,lmax,mm]);

term1=1/(24*\[Omega]^2)*(tmp1 . SSminus . \[CapitalSigma]mat-4*r*\[CapitalDelta][r]*tmp2 . SSminus-2*a^2*DDPm2vec . SSminus . tmp3 . cosmix);
tmp1=(Map[Dop,\[CapitalDelta][r]*Map[Ddag,DdagDdagPp2vec]]+Map[Ddag,\[CapitalDelta][r]*Map[Dop,DdagDdagPp2vec]]);
tmp2=Map[D[#,r]&,DdagDdagPp2vec];
term2=-1/(24*\[Omega]^2)*(tmp1 . SSplus . \[CapitalSigma]mat-4*r*\[CapitalDelta][r]*tmp2 . SSplus-2*a^2*DdagDdagPp2vec . SSplus . tmp3 . cosmix);
tmp=(2*a^2*r*cosmix-I*a*(r^2*idmat+3*a^2*CosMixMatrix[0,2,lmax,mm]));
term3=1/(3*I*\[Omega])*(DDPm2vec . SSminus . \[CapitalSigma]mat . \[Rho]cmat-DPm2vec . SSminus . \[Rho]cmat . \[Rho]cmat
-DDPm2vec . MMm2 . SinMixMatrix[-1,0,lmax,mm] . tmp
-2*I*a*DPm2vec . MMm2 . SinMixMatrix[-1,0,lmax,mm] . \[Rho]mat);
term4=1/(3*I*\[Omega])*(DdagDdagPp2vec . SSplus . \[CapitalSigma]mat . \[Rho]cmat-DdagPp2vec . SSplus . \[Rho]cmat . \[Rho]cmat
+DdagDdagPp2vec . MMp2 . SinMixMatrix[1,0,lmax,mm] . tmp
+2*I*a*DdagPp2vec . MMp2 . SinMixMatrix[1,0,lmax,mm] . \[Rho]mat);
tmp=term1+term2+term3+term4;
\[CapitalSigma]\[CapitalDelta]hraws2lplm=(tmp/(-6*I*M*\[Omega]));
(* l+ m+ *)
t0=r^2*idmat+a^2*CosMixMatrix[1,2,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
cosmix=CosMixMatrix[1,1,lmax,mm];
MM0=\[Lambda]hat-a*\[Omega]*sinmix;
Dhvec=Map[Dop,hvec];
D\[Kappa]mat=Map[Dop,\[Kappa]mat];
term1=1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0;
term2=-a/\[Omega]*(r*Dhvec . Bmat0 . sinmix-hvec . Bmat0 . MM0 . cosmix);
tmp=r*idmat+I*a*cosmix;
term3=4*(ones . D\[Kappa]mat) . Bmat0 . MM0 . tmp;
term4=-4*((ones . \[Kappa]mat) . Bmat0 . MM0+I*a*(ones . D\[Kappa]mat) . Bmat0 . sinmix);
\[Rho]hraws0lpmp=(term1+term2+term3+term4)/.awsubs;
(* l+ m- *)
t0=r^2*idmat+a^2*CosMixMatrix[-1,2,lmax,mm];
sinmix=SinMixMatrix[0,-1,lmax,mm];
MM0=\[Lambda]hat-a*\[Omega]*sinmix;
Dhvec=Map[Dop,hvec];
D\[Kappa]mat=Map[Dop,\[Kappa]mat];
term1=-1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0;
term2=a/\[Omega]*(r*Dhvec . Bmat0 . sinmix+hvec . Bmat0 . MM0 . CosMixMatrix[-1,1,lmax,mm]);
tmp=r*idmat-I*a*CosMixMatrix[-1,1,lmax,mm];
term3=-4*((ones . D\[Kappa]mat) . Bmat0 . MM0 . tmp);
term4=4*((ones . \[Kappa]mat) . Bmat0 . MM0+I*a*(ones . D\[Kappa]mat) . Bmat0 . sinmix);
\[Rho]chraws0lpmm=(term1+term2+term3+term4)/.awsubs;
(* l- m+ *)
t0=r^2*idmat+a^2*CosMixMatrix[1,2,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
cosmix=CosMixMatrix[1,1,lmax,mm];
MM0=\[Lambda]hat-a*\[Omega]*sinmix;
Dhvec=Map[Ddag,hvec];
D\[Kappa]mat=Map[Ddag,\[Kappa]mat];
term1=-1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0;
term2=-a/\[Omega]*(r*Dhvec . Bmat0 . sinmix-hvec . Bmat0 . MM0 . cosmix);
tmp=r*idmat-I*a*cosmix;
term3=4*((ones . D\[Kappa]mat) . Bmat0 . MM0 . tmp);
term4=-4*((ones . \[Kappa]mat) . Bmat0 . MM0-I*a*(ones . D\[Kappa]mat) . Bmat0 . sinmix);
\[Rho]chraws0lmmp=(term1+term2+term3+term4)/.awsubs;
(* l- m- *)
t0=r^2*idmat+a^2*CosMixMatrix[-1,2,lmax,mm];
sinmix=SinMixMatrix[0,-1,lmax,mm];
cosmix=CosMixMatrix[-1,1,lmax,mm];
MM0=\[Lambda]hat-a*\[Omega]*sinmix;
Dhvec=Map[Ddag,hvec];
D\[Kappa]mat=Map[Ddag,\[Kappa]mat];
term1=1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0;
term2=a/\[Omega]*(r*Dhvec . Bmat0 . sinmix+hvec . Bmat0 . MM0 . cosmix);
tmp=r*idmat+I*a*cosmix;
term3=-4*((ones . D\[Kappa]mat) . Bmat0 . MM0 . tmp);
term4=4*((ones . \[Kappa]mat) . Bmat0 . MM0-I*a*(ones . D\[Kappa]mat) . Bmat0 . sinmix);
\[Rho]hraws0lmmm=(term1+term2+term3+term4)/.awsubs;
(* l+ l- *)
cosmix=CosMixMatrix[0,1,lmax,mm];
\[CapitalSigma]mat=r^2*idmat+a^2*CosMixMatrix[0,2,lmax,mm];
tmp1=(Map[Dop,\[CapitalDelta][r]*Map[Ddag,\[Kappa]mat]]+Map[Ddag,\[CapitalDelta][r]*Map[Dop,\[Kappa]mat]])/.\[Kappa]lsubs/.awsubs;
tmp2=Map[D[#,r]&,\[Kappa]mat]/.awsubs;
tmp3=\[Lambda]hat . (SinMixMatrix[-1,0,lmax,mm]-SinMixMatrix[1,0,lmax,mm]);
term1=-hvec . Bmat0 . ((K[r]/\[Omega])*idmat-2*\[CapitalSigma]mat) . \[CapitalSigma]mat;
term2=-2*(ones . tmp1) . Bmat0 . \[CapitalSigma]mat;
term3=8*r*\[CapitalDelta][r]*(ones . tmp2) . Bmat0;
term4=4*a^2*(ones . \[Kappa]mat) . Bmat0 . tmp3 . cosmix;
\[CapitalSigma]\[CapitalDelta]hraws0lplm=(term1+term2+term3+term4)/.hlsubs/.\[Gamma]llsubs/.awsubs;


t0=-(Bmatm1+signmat . Bmatp1) . \[Lambda]hat;
t1=a\[Omega]0*(Bmatm1 . SinMixMatrix[-1,0,lmax,mm]+signmat . Bmatp1 . SinMixMatrix[1,0,lmax,mm]);
SS0=t0+t1;
Ss1lplp=SS0;
Ss1lmlm=SS0;
Ss1mpmp=Bmatp1 . (-\[Lambda]2hat+a\[Omega]0*SinMixMatrix[1,2,lmax,mm]);
Ss1mmmm=Bmatm1 . (\[Lambda]2hat-a\[Omega]0*SinMixMatrix[-1,-2,lmax,mm]);
aPm1vec=Pm1vec;
aPp1vec=Pp1vec;
aDm1Pp1=Map[Dop[-1,#]&,Pp1vec]; (*Table[Dop[-1,Pp1vec[[ll+1]]],{ll,0,lmax}]; *)
aDdagm1Pm1=Map[Ddag[-1,#]&,Pm1vec]; 
aDdagPp1=Map[Ddag,Pp1vec]/.\[CapitalDelta]Ksimp;
aDPm1=Map[Dop,Pm1vec]/.\[CapitalDelta]Ksimp;
acalPs1=(aDPm1+signmat . aDdagPp1)/.\[CapitalDelta]Ksimp; (* check *)
aDcalPs1=Map[Dop,acalPs1]/.\[CapitalDelta]Ksimp;
aDdagcalPs1=Map[Ddag,acalPs1]/.\[CapitalDelta]Ksimp;
s1cos11=CosMixMatrix[1,1,lmax,mm];
s1cosm11=CosMixMatrix[-1,1,lmax,mm];
s1cos01=CosMixMatrix[0,1,lmax,mm];
s1cos02=CosMixMatrix[0,2,lmax,mm];
s1sin01=SinMixMatrix[0,1,lmax,mm];
s1sin0m1=SinMixMatrix[0,-1,lmax,mm];
s1sinm10=SinMixMatrix[-1,0,lmax,mm];
s1sin10=SinMixMatrix[1,0,lmax,mm];

Calcs1components[rtry_?NumericQ]:=Module[{
line1,line2,tmp1,tmp2,term1,term2,term3,tmp,\[CapitalDelta]Ktbl,myP1subs,
Pm1vec,Pp1vec,Dm1Pp1,Ddagm1Pm1,DdagPp1,DPm1,calPs1,DcalPs1,
DdagcalPs1,cosmix,sinmix,\[CapitalSigma]mat,
hraws1lplp,hraws1lmlm,hraws1mpmp,hraws1mmmm,
\[Rho]hraws1lpmp,\[Rho]chraws1lpmm,\[Rho]chraws1lmmp,\[Rho]hraws1lmmm,\[CapitalSigma]\[CapitalDelta]hraws1lplm},

\[CapitalDelta]Ktbl=EvaluateRHS[\[CapitalDelta]Ksubs,rtry]/.awsubs;
myP1subs=GetPsubsAll[1,rtry];
{Pm1vec,Pp1vec,Dm1Pp1,Ddagm1Pm1,
DdagPp1,DPm1,
calPs1,DcalPs1,DdagcalPs1}={aPm1vec,aPp1vec,aDm1Pp1,aDdagm1Pm1,aDdagPp1,aDPm1,
acalPs1,aDcalPs1,aDdagcalPs1}/.\[CapitalDelta]Ktbl/.awsubs/.myP1subs/.{r->rtry};

(* l+ l+, l- l- *)
tmp=2/\[CapitalDelta][r]*Dm1Pp1;
hraws1lplp=tmp . signmat . Ss1lplp/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};
tmp=2/\[CapitalDelta][r]*Ddagm1Pm1;
hraws1lmlm=tmp . Ss1lmlm/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* m+ m+, m- m- *)
tmp=2*(signmat . DdagPp1+DPm1);
hraws1mpmp=tmp . signmat . Ss1mpmp/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};
hraws1mmmm=-tmp . Ss1mmmm/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* l+ m+ *)
cosmix=s1cos11;
sinmix=s1sin01;
line1=(signmat . (r*DcalPs1-2*calPs1)) . Bmatp1+I*a*(signmat . DcalPs1) . Bmatp1 . cosmix;
tmp1=1/\[CapitalDelta][r]*signmat . Pp1vec;
tmp2=-r*\[Lambda]hat-I*a*\[Lambda]hat . cosmix+(a*\[Omega]*r+2*I*a)*sinmix+I*a^2*\[Omega]*sinmix . cosmix;
line2=tmp1 . SS0 . tmp2;
\[Rho]hraws1lpmp=(line1+line2)/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* l+ m- *)
cosmix=s1cosm11;
sinmix=s1sin0m1;
line1=-(r*DcalPs1-2*calPs1) . Bmatm1+I*a*DcalPs1 . Bmatm1 . CosMixMatrix[-1,1,lmax,mm];
tmp1=1/\[CapitalDelta][r]*signmat . Pp1vec;
tmp2=r*\[Lambda]hat-I*a*\[Lambda]hat . cosmix-(a*\[Omega]*r+2*I*a)*sinmix+I*a^2*\[Omega]*sinmix . cosmix;
line2=tmp1 . SS0 . tmp2;
\[Rho]chraws1lpmm=(line1+line2)/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* l- m+ *)
cosmix=s1cos11;
sinmix=s1sin01;
tmp1=signmat . (r*DdagcalPs1-2*calPs1);
tmp2=signmat . DdagcalPs1;
line1 =tmp1 . Bmatp1-I*a*tmp2 . Bmatp1 . cosmix;
tmp1=1/\[CapitalDelta][r]*Pm1vec;
tmp2=-r*\[Lambda]hat+I*a*\[Lambda]hat . cosmix+(a*\[Omega]*r-2*I*a)*sinmix-I*a^2*\[Omega]*sinmix . cosmix;
line2=tmp1 . SS0 . tmp2;
\[Rho]chraws1lmmp=(line1+line2)/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* l- m- *)
cosmix=s1cosm11;
sinmix=s1sin0m1;
line1=-(r*DdagcalPs1-2*calPs1) . Bmatm1-I*a*DdagcalPs1 . Bmatm1 . cosmix;
tmp1=1/\[CapitalDelta][r]*Pm1vec;
tmp2=r*\[Lambda]hat+I*a*\[Lambda]hat . cosmix-(a*\[Omega]*r-2*I*a)*sinmix-I*a^2*\[Omega]*sinmix . cosmix;
line2=tmp1 . SS0 . tmp2;
\[Rho]hraws1lmmm=(line1+line2)/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* l+ l- *)
cosmix=s1cos01; 
\[CapitalSigma]mat=r^2*idmat+a^2*s1cos02;
term1=calPs1 . SS0 . \[CapitalSigma]mat;
term2=-2*r*(Pm1vec+signmat . Pp1vec) . SS0;
term3=2*a^2*calPs1 . (Bmatm1 . s1sinm10 -signmat . Bmatp1 . s1sin10) . cosmix;
\[CapitalSigma]\[CapitalDelta]hraws1lplm=(term1+term2+term3)/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

{hraws1lplp,hraws1lmlm,hraws1mpmp,hraws1mmmm,\[Rho]hraws1lpmp,\[Rho]chraws1lpmm,\[Rho]chraws1lmmp,\[Rho]hraws1lmmm,\[CapitalSigma]\[CapitalDelta]hraws1lplm}
];


(* This is a bad way to do it, change in future *)
tiny\[Epsilon]=SetPrecision[10^(-prec+3),prec];
rsR[[1]]=r0+tiny\[Epsilon];
If[rsR[[1]]>r0,Print["True"];];

(* Ngrids (currently for testing)*)
lNGrid=NGrid["gridL", Transpose@{rsL,rsL}];
rNGrid=NGrid["gridR", Transpose@{rsR,rsR}];

Module[{buildPp2,buildPm2},
buildPp2=Function[{rtry},
 Table[If[ll>=lmins2,{Pp2[ll][r],Pp2[ll]'[r],Pp2[ll]''[r],Pp2[ll]'''[r],Pp2[ll]''''[r]},{0,0,0,0,0}],{ll,0,lmax}]/.GetPsubsAll[2,rtry]];
 buildPm2=Function[{rtry},
 Table[If[ll>=lmins2,{Pm2[ll][r],Pm2[ll]'[r],Pm2[ll]''[r],Pm2[ll]'''[r],Pm2[ll]''''[r]},{0,0,0,0,0}],{ll,0,lmax}]/.GetPsubsAll[2,rtry]];
 Pp2ngridL= buildPp2/@lNGrid;
 Pp2ngridR= buildPp2/@rNGrid;
 Pm2ngridL= buildPm2/@lNGrid;
 Pm2ngridR= buildPm2/@rNGrid;
]//EchoT["Build Spin-2 solutions Ngrid"];

Module[{buildh,build\[Kappa]},
buildh=Function[{rtry},
	Module[{\[CapitalDelta]Ktbl,tmp,rtest=rtry,s0subs},
	\[CapitalDelta]Ktbl=EvaluateRHS[\[CapitalDelta]Ksubs,rtest]/.awsubs;
	tmp=GetSpin0subs[rtest];
	s0subs=Join[tmp,\[Kappa]lsubs/.tmp];
	Table[If[ll>=lmins2,{h[ll][r],h[ll]'[r],h[ll]''[r]},{0,0,0}],{ll,0,lmax}]/.\[Kappa]lsubs/.hlsubs/.\[Gamma]llsubs/.awsubs/.s0subs/.\[CapitalDelta]Ktbl/.{r->rtest}
	]
 ];
 build\[Kappa]=Function[{rtry},
	Module[{\[CapitalDelta]Ktbl,tmp,rtest=rtry,s0subs},
	\[CapitalDelta]Ktbl=EvaluateRHS[\[CapitalDelta]Ksubs,rtest]/.awsubs;
	tmp=GetSpin0subs[rtest];
	s0subs=Join[tmp,\[Kappa]lsubs/.tmp];
	Table[If[ll>=lmins2,{\[Kappa][ll][r],\[Kappa][ll]'[r],\[Kappa][ll]''[r]},{0,0,0}],{ll,0,lmax}]/.\[Kappa]lsubs/.hlsubs/.\[Gamma]llsubs/.awsubs/.s0subs/.\[CapitalDelta]Ktbl/.{r->rtest}
	]
 ];
 hngridL= buildh/@lNGrid;
 hngridR= buildh/@rNGrid;
 \[Kappa]ngridL= build\[Kappa]/@lNGrid;
 \[Kappa]ngridR= build\[Kappa]/@rNGrid;
]//EchoT["Build Spin-0 solutions Ngrid"];

(*Module[{rtest=7,(*test1,test2,*)Psubs,hvec,\[Kappa]vec,rtry,\[CapitalDelta]Ktbl,myP2subs,tmp,s0subs},
	EchoTiming[test1=Calcs0components[rtest],"old"];
	\[CapitalDelta]Ktbl=EvaluateRHS[\[CapitalDelta]Ksubs,rtest]/.awsubs;
	tmp=GetSpin0subs[rtest];
	s0subs=Join[tmp,\[Kappa]lsubs/.tmp];
	hvec=Table[If[ll>=lmins2,{h[ll][r],h[ll]'[r],h[ll]''[r]},{0,0,0}],{ll,0,lmax}]/.\[Kappa]lsubs/.hlsubs/.\[Gamma]llsubs/.awsubs/.s0subs/.\[CapitalDelta]Ktbl/.{r->rtest};
	\[Kappa]vec=Table[If[ll>=lmins2,{\[Kappa][ll][r],\[Kappa][ll]'[r],\[Kappa][ll]''[r]},{0,0,0}],{ll,0,lmax}]/.\[Kappa]lsubs/.hlsubs/.\[Gamma]llsubs/.awsubs/.s0subs/.\[CapitalDelta]Ktbl/.{r->rtest};
	Print[{hvec,\[Kappa]vec}];
	EchoTiming[test2=CalculateSpin0Contributions[a0,mm,\[Omega]0,lmax,rtest,{hvec,\[Kappa]vec},{Bmat0,\[Gamma]mat,eigs0}],"new"];
];*)


zeros=Table[SetPrecision[0,prec],{ll,0,lmax}];

EchoTiming[s2Rgrid=CalculateSpin2Contributions[a0,mm,\[Omega]0,lmax,rNGrid,{Pp2ngridR,Pm2ngridR},{Bmatm2,Bmatp2}],"s = 2, UP ..."];
s2R=Transpose[s2Rgrid/. NGrid[_,data_]:> data,{2,3,1}];
Print["s = 1, UP ..."];
(s1R=GetComp[1]/@rsR;)//EchoTiming;
(*Print["s = 0, UP ..."];
(s0R=GetComp[0]/@rsR;)//EchoTiming;*)
EchoTiming[s0Rgrid=CalculateSpin0Contributions[a0,mm,\[Omega]0,lmax,rNGrid,{hngridR,\[Kappa]ngridR},{Bmat0,\[Gamma]mat,eigs0}],"s = 0, UP ..."];
s0R=Transpose[s0Rgrid/. NGrid[_,data_]:> data,{2,3,1}];

EchoTiming[s2Lgrid=CalculateSpin2Contributions[a0,mm,\[Omega]0,lmax,lNGrid,{Pp2ngridL,Pm2ngridL},{Bmatm2,Bmatp2}],"s = 2, IN ..."];
s2L=Transpose[s2Lgrid/. NGrid[_,data_]:> data,{2,3,1}];
Print["s = 1, IN ..."];
(s1L=GetComp[1]/@rsL;)//EchoTiming;
(*Print["s = 0, IN ..."];
(s0L=GetComp[0]/@rsL;)//EchoTiming;*)
EchoTiming[s0Lgrid=CalculateSpin0Contributions[a0,mm,\[Omega]0,lmax,lNGrid,{hngridL,\[Kappa]ngridL},{Bmat0,\[Gamma]mat,eigs0}],"s = 0, IN ..."];
s0L=Transpose[s0Lgrid/. NGrid[_,data_]:> data,{2,3,1}];

sallR=s2R+s1R+s0R;
sallL=s2L+s1L+s0L;

(* Save in a binary format. *)
ExportOutput[{sallL,sallR,s2L,s2R,s1L,s1R,s0L,s0R},directory,iConfig,dformat]
]




End[];
EndPackage[];
