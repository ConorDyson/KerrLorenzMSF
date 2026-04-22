(* ::Package:: *)

(* ::Title:: *)
(*MetricReconstructRadiative Package*)


BeginPackage["MetricReconstructRadiative`",
{
"OrbitalData`",
"KerrOrbitalParameters`",
"SWSHdecomp`",
"SpinWeightedSpheroidalHarmonicsFT`",
"NGrid`",
"MSTMode`",
"MST`"
}];
(*Needs["SpinWeightedSpheroidalHarmonics`"]
Needs["Teukolsky`"]*)


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


EchoT[lbl_][expr_]:=EchoTiming[expr,lbl]


(* ::Input::Initialization:: *)
FudgeZero[x_?PossibleZeroQ]:=SetPrecision[10^-Accuracy[x],Max[2,$MinPrecision]];
FudgeZero[0]:= SetPrecision[10^-$TeukolskyWorkingPrecision,Max[2,$MinPrecision]]
FudgeZero[x_]:=x


extractUpToLastHyphen[str_String] := StringJoin[Riffle[Most[StringSplit[str, "-"]], "-"]];
extractUpToSecondLastHyphen[str_String] := StringJoin[Riffle[Most[StringSplit[extractUpToLastHyphen[str], "-"]], "-"]];


ConstructSolution[{insols_,upsols_},{Jump0_,Jump1_}]:=Module[
{Rup0,Rup1,Rin0,Rin1,\[Alpha]in,\[Alpha]up},
Rin0=insols[[1,2,1]];
Rin1=insols[[2,2,1]];
Rup0=upsols[[1,2,1]];
Rup1=upsols[[2,2,1]];
{\[Alpha]up,\[Alpha]in}={(Jump1 Rin0-Jump0 Rin1)/(-Rin1 Rup0+Rin0 Rup1),(Jump1 Rup0-Jump0 Rup1)/(-Rin1 Rup0+Rin0 Rup1)};
Return[{\[Alpha]in insols,\[Alpha]up upsols}]
]


ConstructSolutionh\[Kappa][{{inhs_,in\[Kappa]s_},{uphs_,up\[Kappa]s_}},{hj0_,hj1_},{\[Kappa]j0_,\[Kappa]j1_}]:=Module[
{hup0,hup1,hin0,hin1,\[Kappa]up0,\[Kappa]up1,\[Kappa]in0,\[Kappa]in1,\[Alpha]hin,\[Alpha]hup,\[Alpha]\[Kappa]up,\[Alpha]\[Kappa]in},
hin0=inhs[[1,2,1]];
hin1=inhs[[2,2,1]];
hup0=uphs[[1,2,1]];
hup1=uphs[[2,2,1]];
\[Kappa]in0=in\[Kappa]s[[1,2,1]];
\[Kappa]in1=in\[Kappa]s[[2,2,1]];
\[Kappa]up0=up\[Kappa]s[[1,2,1]];
\[Kappa]up1=up\[Kappa]s[[2,2,1]];
{\[Alpha]hup,\[Alpha]hin}={(hj1 hin0-hj0 hin1)/(-hin1 hup0+hin0 hup1),(hj1 hup0-hj0 hup1)/(-hin1 hup0+hin0 hup1)};
\[Alpha]\[Kappa]up=(-hin1 (\[Alpha]hin \[Kappa]in0+\[Kappa]j0-\[Alpha]hup \[Kappa]up0)+hin0 (\[Alpha]hin \[Kappa]in1+\[Kappa]j1-\[Alpha]hup \[Kappa]up1))/(-hin1 hup0+hin0 hup1);\[Alpha]\[Kappa]in=(-hup1 (\[Alpha]hin \[Kappa]in0+\[Kappa]j0-\[Alpha]hup \[Kappa]up0)+hup0 (\[Alpha]hin \[Kappa]in1+\[Kappa]j1-\[Alpha]hup \[Kappa]up1))/(-hin1 hup0+hin0 hup1);
Return[{\[Alpha]hin inhs,\[Alpha]hup uphs,\[Alpha]\[Kappa]in*inhs+\[Alpha]hin in\[Kappa]s,\[Alpha]\[Kappa]up*uphs+\[Alpha]hup up\[Kappa]s}]
]


(* ::Subsection::Closed:: *)
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


DefTeukSubs["\[Kappa]",{\[Kappa]_Symbol,h_Symbol},{\[CapitalDelta]_Symbol,K_Symbol},{M_,a_,m_,\[Omega]_,\[Lambda]0_,\[Gamma]ll_},{r_Symbol}]:={
Derivative[2][\[Kappa]][r]->((r^2+a^2 \[Gamma]ll) h[r] \[CapitalDelta][r]-2 K[r]^2 \[Kappa][r]+2 \[CapitalDelta][r] (\[Lambda]0 \[Kappa][r]+2 (M-r) Derivative[1][\[Kappa]][r]))/(2 \[CapitalDelta][r]^2),
\!\(\*SuperscriptBox[\(\[Kappa]\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[r]->1/(2 \[CapitalDelta][r]^3) (2 h[r] \[CapitalDelta][r] (2 (M-r) (r^2+a^2 \[Gamma]ll)+r \[CapitalDelta][r])-8 r \[Omega] K[r] \[CapitalDelta][r] \[Kappa][r]-2 K[r]^2 (6 (M-r) \[Kappa][r]+\[CapitalDelta][r] Derivative[1][\[Kappa]][r])+\[CapitalDelta][r] (8 (M-r) \[Lambda]0 \[Kappa][r]+16 (M-r)^2 Derivative[1][\[Kappa]][r]+\[CapitalDelta][r] ((r^2+a^2 \[Gamma]ll) Derivative[1][h][r]+2 (-2+\[Lambda]0) Derivative[1][\[Kappa]][r]))),
\!\(\*SuperscriptBox[\(\[Kappa]\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[r]->1/(2 \[CapitalDelta][r]^4) (h[r] \[CapitalDelta][r] (24 (M-r)^2 (r^2+a^2 \[Gamma]ll)-(r^2+a^2 \[Gamma]ll) K[r]^2+(12 M r+r^2 (-18+\[Lambda]0)+a^2 \[Gamma]ll (-6+\[Lambda]0)) \[CapitalDelta][r]+2 \[CapitalDelta][r]^2)+2 K[r]^4 \[Kappa][r]-4 K[r]^2 ((22 (M-r)^2+(-4+\[Lambda]0) \[CapitalDelta][r]) \[Kappa][r]+6 (M-r) \[CapitalDelta][r] Derivative[1][\[Kappa]][r])-8 \[Omega] K[r] \[CapitalDelta][r] ((10 (M-r) r+\[CapitalDelta][r]) \[Kappa][r]+2 r \[CapitalDelta][r] Derivative[1][\[Kappa]][r])+\[CapitalDelta][r] (2 (24 (M-r)^2 \[Lambda]0+(-6 \[Lambda]0+\[Lambda]0^2-8 r^2 \[Omega]^2) \[CapitalDelta][r]) \[Kappa][r]+96 (M-r)^3 Derivative[1][\[Kappa]][r]+2 (M-r) \[CapitalDelta][r] (3 (r^2+a^2 \[Gamma]ll) Derivative[1][h][r]+8 (-3+\[Lambda]0) Derivative[1][\[Kappa]][r])+\[CapitalDelta][r]^2 (4 r Derivative[1][h][r]+(r^2+a^2 \[Gamma]ll) Derivative[2][h][r])))
}


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


(* ::Input::Initialization:: *)
Pp2fromPm2[{Pm20_,Pm21_,Pm22_,Pm23_,Pm24_},r_NGrid,{l_,m_,a_,\[Omega]_, \[CapitalLambda]_}]:=Module[
{\[CapitalDelta]r,Kr,CC},
\[CapitalDelta]r=r^2-2r+a^2;
Kr=-a m+(a^2+r^2) \[Omega];
CC=(-1)^(l+m) 12 I  \[Omega]+Sqrt[\[CapitalLambda]^2 (2+\[CapitalLambda])^2+96 a^2 \[CapitalLambda] \[Omega]^2+8 a \[CapitalLambda] (6+5 \[CapitalLambda]) \[Omega] (m-a \[Omega])+144 a^2 \[Omega]^2 (m-a \[Omega])^2];
Return[{1/(CC \[CapitalDelta]r^2) (8 Kr^4 Pm20+8 I Kr^3 Pm21 \[CapitalDelta]r+8 Kr^2 Pm20 ((-1+r)^2-\[CapitalDelta]r (1+\[CapitalLambda]+3 I r \[Omega]))-4 I Kr \[CapitalDelta]r (-2 Pm21 (-1+r)^2+Pm21 \[CapitalDelta]r (2+\[CapitalLambda])+5 I Pm20 \[CapitalDelta]r \[Omega]+Pm20 (1-r) (\[CapitalLambda]+8 I r \[Omega]))+\[CapitalDelta]r^2 (8 I Pm21 (r-r^2+\[CapitalDelta]r) \[Omega]+Pm20 (\[CapitalLambda] (2+\[CapitalLambda])+4 I (3+r \[CapitalLambda]) \[Omega]+12 r^2 \[Omega]^2))),1/(CC \[CapitalDelta]r^3) (-8 I Kr^5 Pm20+8 Kr^4 Pm21 \[CapitalDelta]r-4 I Kr^3 Pm20 (2+2 (-2+r) r-\[CapitalDelta]r (2+3 \[CapitalLambda]))+\[CapitalDelta]r^3 (Pm21 \[CapitalLambda] (2+\[CapitalLambda])-4 I (-3 Pm20 \[CapitalLambda]+Pm21 (3+r \[CapitalLambda])) \[Omega]+12 Pm21 r^2 \[Omega]^2)-4 I Kr \[CapitalDelta]r^2 (Pm21 (-1+r) \[CapitalLambda]+I Pm21 (-8 (-1+r) r+5 \[CapitalDelta]r) \[Omega]+Pm20 (\[CapitalLambda]+\[CapitalLambda]^2+24 r^2 \[Omega]^2))+8 Kr^2 \[CapitalDelta]r (I Pm20 (7 (-1+r) r-4 \[CapitalDelta]r) \[Omega]+Pm21 (1-\[CapitalDelta]r (1+\[CapitalLambda])+r (-2+r+3 I \[CapitalDelta]r \[Omega])))),1/(CC \[CapitalDelta]r^4) (-8 Kr^6 Pm20+8 I Kr^5 (2 Pm20 (-1+r)-Pm21 \[CapitalDelta]r)+4 Kr^3 (2 I ((-1+r)^2-\[CapitalDelta]r) (2 Pm20 (-1+r)-Pm21 \[CapitalDelta]r)+3 I \[CapitalDelta]r (Pm20-Pm20 r+Pm21 \[CapitalDelta]r) \[CapitalLambda]+\[CapitalDelta]r (32 Pm20 (-1+r) r-5 Pm20 \[CapitalDelta]r+16 Pm21 r \[CapitalDelta]r) \[Omega])+\[CapitalDelta]r^3 (\[CapitalLambda] (2+\[CapitalLambda]) (2 Pm21 (-1+r)+Pm20 \[CapitalLambda])-4 I (Pm20 (3+r (-2+\[CapitalLambda])) \[CapitalLambda]+Pm21 (-6+6 r-4 r \[CapitalLambda]+4 r^2 \[CapitalLambda]-2 \[CapitalDelta]r \[CapitalLambda])) \[Omega]+4 r (2 Pm21 (-5 (-1+r) r+8 \[CapitalDelta]r)+Pm20 (24+11 r \[CapitalLambda])) \[Omega]^2-96 I Pm20 r^3 \[Omega]^3)+4 Kr \[CapitalDelta]r^2 (-I Pm21 \[CapitalDelta]r \[CapitalLambda] (2+\[CapitalLambda])-2 Pm21 (\[CapitalDelta]r+r (-4-4 (-2+r) r+\[CapitalDelta]r (7+4 \[CapitalLambda]))) \[Omega]+Pm20 \[Omega] (12-12 r+4 r \[CapitalLambda]-4 r^2 \[CapitalLambda]+5 \[CapitalDelta]r \[CapitalLambda]+4 I r (7 (-1+r) r-10 \[CapitalDelta]r) \[Omega]))-8 Kr^4 (Pm20+2 Pm21 (-1+r) \[CapitalDelta]r-Pm20 \[CapitalDelta]r (1+2 \[CapitalLambda])+Pm20 r (-2+r+5 I \[CapitalDelta]r \[Omega]))+Kr^2 \[CapitalDelta]r (16 Pm21 (-1+r) (-(-1+r)^2+\[CapitalDelta]r)+8 I Pm21 (7 (-1+r) r-\[CapitalDelta]r) \[CapitalDelta]r \[Omega]+Pm20 (-\[CapitalLambda] (8+8 (-2+r) r+\[CapitalDelta]r (10+9 \[CapitalLambda]))-4 I (20 r ((-1+r)^2-\[CapitalDelta]r)+7 \[CapitalDelta]r-9 r \[CapitalDelta]r \[CapitalLambda]) \[Omega]-204 r^2 \[CapitalDelta]r \[Omega]^2))),1/(CC \[CapitalDelta]r^5) (8 I Kr^7 Pm20+8 Kr^6 (6 Pm20 (-1+r)-Pm21 \[CapitalDelta]r)+4 Kr \[CapitalDelta]r^2 (-I Pm20 \[CapitalLambda] (2+\[CapitalLambda]) (2 (-1+r)^2+\[CapitalDelta]r \[CapitalLambda])+(8 (-1+r)^2 (3 Pm20-4 Pm21 (-1+r) r)+4 (-3 Pm20+5 Pm21 (-1+r) (-1+2 r)) \[CapitalDelta]r-14 Pm21 \[CapitalDelta]r^2+\[CapitalDelta]r (Pm20 (4-18 r)-4 Pm21 (-1+r) r-3 Pm21 \[CapitalDelta]r) \[CapitalLambda]-9 Pm20 r \[CapitalDelta]r \[CapitalLambda]^2) \[Omega]-4 I (3 Pm21 r \[CapitalDelta]r (-7 (-1+r) r+4 \[CapitalDelta]r)+Pm20 (-44 r^3+22 r^4+29 r \[CapitalDelta]r+10 \[CapitalDelta]r^2+r^2 (22-17 \[CapitalDelta]r+7 \[CapitalDelta]r \[CapitalLambda]))) \[Omega]^2-204 Pm20 r^3 \[CapitalDelta]r \[Omega]^3)+\[CapitalDelta]r^3 (Pm21 \[CapitalDelta]r \[CapitalLambda] (2+\[CapitalLambda])^2-4 I (2+\[CapitalLambda]) (-Pm20 (2 (-1+r) r+\[CapitalDelta]r) \[CapitalLambda]+3 Pm21 (\[CapitalDelta]r+r \[CapitalDelta]r \[CapitalLambda])) \[Omega]+4 (Pm21 (-32 r^3+16 r^4+8 r \[CapitalDelta]r+16 \[CapitalDelta]r^2+r^2 (16-26 \[CapitalDelta]r-5 \[CapitalDelta]r \[CapitalLambda]))+8 Pm20 (3 \[CapitalDelta]r+r (3+r (-3+\[CapitalLambda])-r^2 \[CapitalLambda]+4 \[CapitalDelta]r \[CapitalLambda]))) \[Omega]^2-96 I r^2 (Pm21 r \[CapitalDelta]r+Pm20 (r-r^2+\[CapitalDelta]r)) \[Omega]^3)-4 I Kr^5 (-12 Pm21 (-1+r) \[CapitalDelta]r+Pm20 (14-2 \[CapitalDelta]r+5 \[CapitalDelta]r \[CapitalLambda]+2 r (-14+7 r-12 I \[CapitalDelta]r \[Omega])))+8 Kr^4 (2 Pm20 (-1+r) (3+3 (-2+r) r-3 \[CapitalDelta]r-4 \[CapitalDelta]r \[CapitalLambda])-I Pm20 \[CapitalDelta]r (5 (-1+r) r+4 \[CapitalDelta]r) \[Omega]+Pm21 \[CapitalDelta]r (7-\[CapitalDelta]r+2 \[CapitalDelta]r \[CapitalLambda]+r (-14+7 r-15 I \[CapitalDelta]r \[Omega])))-4 I Kr^3 (Pm21 \[CapitalDelta]r (3 (-1+r) (-4-4 (-2+r) r+4 \[CapitalDelta]r+3 \[CapitalDelta]r \[CapitalLambda])+I \[CapitalDelta]r (-32 (-1+r) r+11 \[CapitalDelta]r) \[Omega])+Pm20 (16+16 r^4+8 r (-8+5 \[CapitalDelta]r+4 \[CapitalDelta]r \[CapitalLambda])-\[CapitalDelta]r (4 (5+4 \[CapitalLambda])+\[CapitalDelta]r (-4+\[CapitalLambda]+4 \[CapitalLambda]^2+42 I \[Omega]))+2 I r \[CapitalDelta]r (-72+45 \[CapitalDelta]r+16 \[CapitalDelta]r \[CapitalLambda]) \[Omega]+16 r^3 (-4-9 I \[CapitalDelta]r \[Omega])-4 r^2 (-24+\[CapitalDelta]r (5+4 \[CapitalLambda]+12 \[Omega] (-6 I+\[CapitalDelta]r \[Omega])))))+Kr^2 \[CapitalDelta]r (2 Pm20 (2 (-1+r) (8 (-1+r)^2+\[CapitalDelta]r) \[CapitalLambda]+9 (-1+r) \[CapitalDelta]r \[CapitalLambda]^2-2 I (4 (-1+r) r-5 \[CapitalDelta]r) \[CapitalDelta]r \[CapitalLambda] \[Omega]+4 \[Omega] (I (40 (-1+r)^3 r-(-1+r) (-31+38 r) \[CapitalDelta]r+10 \[CapitalDelta]r^2)+3 r (49 (-1+r) r-22 \[CapitalDelta]r) \[CapitalDelta]r \[Omega]))+Pm21 (64+64 r^4-\[CapitalDelta]r (8 (10+\[CapitalLambda])+\[CapitalDelta]r (-16+\[CapitalLambda] (10+9 \[CapitalLambda])+84 I \[Omega]))+r^3 (-256-240 I \[CapitalDelta]r \[Omega])+4 r (-64+4 \[CapitalDelta]r (10+\[CapitalLambda])+3 I \[CapitalDelta]r (-20+20 \[CapitalDelta]r+9 \[CapitalDelta]r \[CapitalLambda]) \[Omega])+4 r^2 (96+\[CapitalDelta]r (-20-2 \[CapitalLambda]+15 \[Omega] (8 I+3 \[CapitalDelta]r \[Omega])))))),1/(CC \[CapitalDelta]r^6) (Kr^4 (8 Kr^4 Pm20+8 I Kr^3 Pm21 \[CapitalDelta]r+8 Kr^2 Pm20 ((-1+r)^2-\[CapitalDelta]r (1+\[CapitalLambda]+3 I r \[Omega]))-4 I Kr \[CapitalDelta]r (-2 Pm21 (-1+r)^2+Pm21 \[CapitalDelta]r (2+\[CapitalLambda])+5 I Pm20 \[CapitalDelta]r \[Omega]+Pm20 (1-r) (\[CapitalLambda]+8 I r \[Omega]))+\[CapitalDelta]r^2 (8 I Pm21 (r-r^2+\[CapitalDelta]r) \[Omega]+Pm20 (\[CapitalLambda] (2+\[CapitalLambda])+4 I (3+r \[CapitalLambda]) \[Omega]+12 r^2 \[Omega]^2)))+8 I Kr^3 (1-r) (8 Kr^4 Pm20+8 I Kr^3 Pm21 \[CapitalDelta]r+8 Kr^2 Pm20 ((-1+r)^2-\[CapitalDelta]r (1+\[CapitalLambda]+3 I r \[Omega]))-4 I Kr \[CapitalDelta]r (-2 Pm21 (-1+r)^2+Pm21 \[CapitalDelta]r (2+\[CapitalLambda])+5 I Pm20 \[CapitalDelta]r \[Omega]+Pm20 (1-r) (\[CapitalLambda]+8 I r \[Omega]))+\[CapitalDelta]r^2 (8 I Pm21 (r-r^2+\[CapitalDelta]r) \[Omega]+Pm20 (\[CapitalLambda] (2+\[CapitalLambda])+4 I (3+r \[CapitalLambda]) \[Omega]+12 r^2 \[Omega]^2)))-4 I Kr (-12 (-1+r)^3-I \[CapitalDelta]r^2 \[Omega]+2 (-1+r) \[CapitalDelta]r (3-\[CapitalLambda]+11 I r \[Omega])) (8 Kr^4 Pm20+8 I Kr^3 Pm21 \[CapitalDelta]r+8 Kr^2 Pm20 ((-1+r)^2-\[CapitalDelta]r (1+\[CapitalLambda]+3 I r \[Omega]))-4 I Kr \[CapitalDelta]r (-2 Pm21 (-1+r)^2+Pm21 \[CapitalDelta]r (2+\[CapitalLambda])+5 I Pm20 \[CapitalDelta]r \[Omega]+Pm20 (1-r) (\[CapitalLambda]+8 I r \[Omega]))+\[CapitalDelta]r^2 (8 I Pm21 (r-r^2+\[CapitalDelta]r) \[Omega]+Pm20 (\[CapitalLambda] (2+\[CapitalLambda])+4 I (3+r \[CapitalLambda]) \[Omega]+12 r^2 \[Omega]^2)))-I \[CapitalDelta]r (48 (-1+r)^2 r \[Omega]+\[CapitalDelta]r (I \[CapitalLambda] (2+\[CapitalLambda])+8 (3+r (-3+2 \[CapitalLambda])) \[Omega]-72 I r^2 \[Omega]^2)) (8 Kr^4 Pm20+8 I Kr^3 Pm21 \[CapitalDelta]r+8 Kr^2 Pm20 ((-1+r)^2-\[CapitalDelta]r (1+\[CapitalLambda]+3 I r \[Omega]))-4 I Kr \[CapitalDelta]r (-2 Pm21 (-1+r)^2+Pm21 \[CapitalDelta]r (2+\[CapitalLambda])+5 I Pm20 \[CapitalDelta]r \[Omega]+Pm20 (1-r) (\[CapitalLambda]+8 I r \[Omega]))+\[CapitalDelta]r^2 (8 I Pm21 (r-r^2+\[CapitalDelta]r) \[Omega]+Pm20 (\[CapitalLambda] (2+\[CapitalLambda])+4 I (3+r \[CapitalLambda]) \[Omega]+12 r^2 \[Omega]^2)))-16 \[CapitalDelta]r (r-r^2+\[CapitalDelta]r) \[Omega] (8 Kr^5 Pm20+8 I Kr^4 Pm21 \[CapitalDelta]r+4 Kr^3 Pm20 (2+2 (-2+r) r-\[CapitalDelta]r (2+3 \[CapitalLambda]))+\[CapitalDelta]r^3 (-12 Pm20 \[CapitalLambda] \[Omega]+Pm21 (I \[CapitalLambda] (2+\[CapitalLambda])+4 (3+r \[CapitalLambda]) \[Omega]+12 I r^2 \[Omega]^2))+4 Kr \[CapitalDelta]r^2 (5 I Pm21 \[CapitalDelta]r \[Omega]+Pm21 (-1+r) (\[CapitalLambda]-8 I r \[Omega])+Pm20 (\[CapitalLambda]+\[CapitalLambda]^2+24 r^2 \[Omega]^2))+8 Kr^2 \[CapitalDelta]r (Pm20 (-7 (-1+r) r+4 \[CapitalDelta]r) \[Omega]+I Pm21 (1-\[CapitalDelta]r (1+\[CapitalLambda])+r (-2+r+3 I \[CapitalDelta]r \[Omega]))))-8 Kr (2-\[CapitalDelta]r+r (-4+2 r-I \[CapitalDelta]r \[Omega])) (8 Kr^5 Pm20+8 I Kr^4 Pm21 \[CapitalDelta]r+4 Kr^3 Pm20 (2+2 (-2+r) r-\[CapitalDelta]r (2+3 \[CapitalLambda]))+\[CapitalDelta]r^3 (-12 Pm20 \[CapitalLambda] \[Omega]+Pm21 (I \[CapitalLambda] (2+\[CapitalLambda])+4 (3+r \[CapitalLambda]) \[Omega]+12 I r^2 \[Omega]^2))+4 Kr \[CapitalDelta]r^2 (5 I Pm21 \[CapitalDelta]r \[Omega]+Pm21 (-1+r) (\[CapitalLambda]-8 I r \[Omega])+Pm20 (\[CapitalLambda]+\[CapitalLambda]^2+24 r^2 \[Omega]^2))+8 Kr^2 \[CapitalDelta]r (Pm20 (-7 (-1+r) r+4 \[CapitalDelta]r) \[Omega]+I Pm21 (1-\[CapitalDelta]r (1+\[CapitalLambda])+r (-2+r+3 I \[CapitalDelta]r \[Omega]))))+2 Kr^2 (-((14 (-1+r)^2+\[CapitalDelta]r (\[CapitalLambda]-8 I r \[Omega])) (8 Kr^4 Pm20+8 I Kr^3 Pm21 \[CapitalDelta]r+8 Kr^2 Pm20 ((-1+r)^2-\[CapitalDelta]r (1+\[CapitalLambda]+3 I r \[Omega]))-4 I Kr \[CapitalDelta]r (-2 Pm21 (-1+r)^2+Pm21 \[CapitalDelta]r (2+\[CapitalLambda])+5 I Pm20 \[CapitalDelta]r \[Omega]+Pm20 (1-r) (\[CapitalLambda]+8 I r \[Omega]))+\[CapitalDelta]r^2 (8 I Pm21 (r-r^2+\[CapitalDelta]r) \[Omega]+Pm20 (\[CapitalLambda] (2+\[CapitalLambda])+4 I (3+r \[CapitalLambda]) \[Omega]+12 r^2 \[Omega]^2))))+2 I (1-r) (8 Kr^5 Pm20+8 I Kr^4 Pm21 \[CapitalDelta]r+4 Kr^3 Pm20 (2+2 (-2+r) r-\[CapitalDelta]r (2+3 \[CapitalLambda]))+\[CapitalDelta]r^3 (-12 Pm20 \[CapitalLambda] \[Omega]+Pm21 (I \[CapitalLambda] (2+\[CapitalLambda])+4 (3+r \[CapitalLambda]) \[Omega]+12 I r^2 \[Omega]^2))+4 Kr \[CapitalDelta]r^2 (5 I Pm21 \[CapitalDelta]r \[Omega]+Pm21 (-1+r) (\[CapitalLambda]-8 I r \[Omega])+Pm20 (\[CapitalLambda]+\[CapitalLambda]^2+24 r^2 \[Omega]^2))+8 Kr^2 \[CapitalDelta]r (Pm20 (-7 (-1+r) r+4 \[CapitalDelta]r) \[Omega]+I Pm21 (1-\[CapitalDelta]r (1+\[CapitalLambda])+r (-2+r+3 I \[CapitalDelta]r \[Omega]))))))}
]
]


Pp1fromPm1[{Pm10_,Pm11_,Pm12_},r_NGrid,{l_,m_,a_,\[Omega]_, \[Lambda]_}]:=Module[
{\[CapitalDelta]r,Kr,BB},
\[CapitalDelta]r=r^2-2r+a^2;
Kr=-a m+(a^2+r^2) \[Omega];
BB=Sqrt[\[Lambda]^2 +4*a*m*\[Omega]-4*a^2*\[Omega]^2 ];
Return[{(-2 Kr^2 Pm10-2 I Kr Pm11 \[CapitalDelta]r+Pm10 \[CapitalDelta]r (\[Lambda]+2 I r \[Omega]))/(BB \[CapitalDelta]r),(2 I Kr^3 Pm10-2 Kr^2 Pm11 \[CapitalDelta]r-2 I Kr Pm10 \[CapitalDelta]r \[Lambda]+\[CapitalDelta]r^2 (Pm11 \[Lambda]+2 I (Pm10-Pm11 r) \[Omega]))/(BB \[CapitalDelta]r^2),((Kr^2-2 I Kr (-1+r)-\[CapitalDelta]r \[Lambda]+4 I r \[CapitalDelta]r \[Omega]) (2 Kr^2 Pm10+2 I Kr Pm11 \[CapitalDelta]r-Pm10 \[CapitalDelta]r (\[Lambda]+2 I r \[Omega])))/(BB \[CapitalDelta]r^3)}]
]


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
0(*If[iI+iJ<=imax,
10^-OptionValue[WorkingPrecision],
Max[
Abs@btemp[[iI+iJ-imax,imax]],
Abs@btemp[[imax,iJ+iI-imax]]
]
]*)
](*["Value"]*),{iI,1,lmax+1-Abs[lmin]},{iJ,1,lmax+1-Abs[lmin]}]
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


(* ::Subsubsection:: *)
(*Diverse matrices*)


\[Lambda]hat[lmax_,m_]:=\[Lambda]hat[lmax,m]=DiagonalMatrix[Table[If[ll<Abs[m],0,Sqrt[ll*(ll+1)]],{ll,0,lmax}]];
\[Lambda]2hat[lmax_,m_]:=\[Lambda]2hat[lmax,m]=DiagonalMatrix[Table[If[ll<Max[Abs[m],1],0,Sqrt[(ll-1)*(ll+2)]],{ll,0,lmax}]];
\[CapitalLambda]hat[lmax_,m_]:=\[CapitalLambda]hat[lmax,m]=\[Lambda]hat[lmax,m] . \[Lambda]2hat[lmax,m];
signmat[lmax_,m_]:=signmat[lmax,m]=DiagonalMatrix[Table[If[ll<Abs[m],0,(-1)^(ll+m)],{ll,0,lmax}]];
idmat[lmax_,m_]:=idmat[lmax,m]=DiagonalMatrix[Table[If[ll<Abs[m],0,1],{ll,0,lmax}]];


(* ::Subsection:: *)
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
Return[{
NGrid["gridL", Transpose@{rsL,rsL}],
NGrid["gridR", Transpose@{rsR,rsR}]
}]
]


(* ::Subsection:: *)
(*Jumps*)


(* ::Subsubsection::Closed:: *)
(*sourcesubs*)


Calcsourcesubs[OD_OrbitalData,{dhlplp_,dhlmlm_,dhmpmp_,dhmmmm_,d\[Rho]chlpmm_},{Y0_,Ym1_,Yp2_,Ym2_}]:=
Module[{a,r,ut,calA,calB,EE,LL,M=1},
a=ODspin[OD];
r=ODperiapsis[OD];
EE=KerrEnergy[OD];
LL=KerrAngularMomentum[OD];
ut=(a^3+(-2+r) r^(5/2)+a r (-2 M+r)+a^2 Sqrt[r] (-2+2 M+r))/(Sqrt[1+(2 a)/r^(3/2)-3/r] r^(3/2) (a^2+r (-2 M+r)));
calA=1/(a^2+(-2+r) r)*(EE*(r^2+a^2)-a*LL);
calB=LL-a*EE;
Return[
{
dhlplp->-16*\[Pi]/((a^2+(-2+r) r)*ut)*calA^2*Y0,
dhlmlm->-16*\[Pi]/((a^2+(-2+r) r)*ut)*calA^2*Y0,
dhmpmp->16*\[Pi]/((a^2+(-2+r) r)*ut)*calB^2*Yp2,
dhmmmm->16*\[Pi]/((a^2+(-2+r) r)*ut)*calB^2*Ym2,
d\[Rho]chlpmm->-16*\[Pi]/((a^2+(-2+r) r)*ut)*(r*I*calA*calB)*Ym1
}
]
]


(* ::Subsubsection::Closed:: *)
(*Spin 2*)


CalcJump[{s:(2|-2),m_,\[CapitalLambda]_},OD_?CircularEquatorialQ,{P2j0_,P2j1_},{S2_,S2p_}]:=Module[{a,\[Sigma],\[CapitalOmega],M=1,r,EE,LL,ut,calA,calB,calD, Jump0, Jump1},
\[Sigma]=Sign[s];
a=ODspin[OD];
\[CapitalOmega]=KerrAzimuthalFrequency[OD];
r=ODperiapsis[OD];
EE=KerrEnergy[OD];
LL=KerrAngularMomentum[OD];
ut=(a^3+(-2+r) r^(5/2)+a r (-2 M+r)+a^2 Sqrt[r] (-2+2 M+r))/(Sqrt[1+(2 a)/r^(3/2)-3/r] r^(3/2) (a^2+r (-2 M+r)));
calA=1/(a^2+(-2+r) r)*(EE*(r^2+a^2)-a*LL);
calB=LL-a*EE;
calD=4*\[Pi]/(r^2*ut);

Jump0=-(1/r)2 I calB calD (a^3 calA m r S2 \[Sigma] \[CapitalOmega]+a m r S2 \[Sigma] (-calB+calA (-2+r) r \[CapitalOmega])+a^2 (-calA r (S2p+m S2 \[Sigma])+calB S2 (I+m r \[Sigma] \[CapitalOmega]))+r (-calA (-2+r) r (S2p+m S2 \[Sigma])+calB S2 (-I+m r^2 \[Sigma] \[CapitalOmega])));
Jump1=1/(r (a^2+(-2+r) r)) calD (calB^2 r S2 ((a^2+(-2+r) r) \[CapitalLambda]-6 I m r (a^2+(-2+r) r) \[Sigma] \[CapitalOmega]+2 I m \[Sigma] (a-a r-a^2 \[CapitalOmega]+r^2 \[CapitalOmega])+2 I m (-1+r) \[Sigma] (-a+a^2 \[CapitalOmega]+r^2 \[CapitalOmega])-m^2 \[Sigma]^2 (a-(a^2+r^2) \[CapitalOmega])^2)+2 I calA (a^2+(-2+r) r) (S2p+m S2 \[Sigma] (1-a \[CapitalOmega])) (calA (a^2+(-2+r) r) (a-I m r \[Sigma] (-1+a \[CapitalOmega]))+calB (2 (a^2+(-2+r) r)+I m r \[Sigma] (a-(a^2+r^2) \[CapitalOmega])))+S2 (calA^2 r (a^2+(-2+r) r)^2 \[CapitalLambda]+calB^2 m \[Sigma] (-2 I r (a-a r-a^2 \[CapitalOmega]+r^2 \[CapitalOmega])+2 I (a^2+(-2+r) r) (a-(a^2+r^2) \[CapitalOmega])-m r \[Sigma] (a-(a^2+r^2) \[CapitalOmega])^2)));
Return[{P2j0 -> Jump0,P2j1 -> Jump1}]
]


(* ::Subsubsection::Closed:: *)
(*Spin 0*)


CalcJump[{0,m_},OD_?CircularEquatorialQ,{hj0_,hj1_},{S0_}]:=Module[{a,M=1,r,ut, Jump0, Jump1},
a=ODspin[OD];
r=ODperiapsis[OD];

ut=(a^3+(-2+r) r^(5/2)+a r (-2 M+r)+a^2 Sqrt[r] (-2+2 M+r))/(Sqrt[1+(2 a)/r^(3/2)-3/r] r^(3/2) (a^2+r (-2 M+r)));
Jump0=0;
Jump1=-16*\[Pi]/(ut*(a^2+(-2+r) r))*S0;
Return[{hj0->Jump0,hj1->Jump1}]
]


(* ::Subsubsection:: *)
(*Spin 1 and \[Kappa]*)


CalcJumps\[Kappa]s1[
{\[Kappa]j0_Symbol,\[Kappa]j1_Symbol,Pm1j0_Symbol,Pm1j1_Symbol},
OD_OrbitalData?CircularEquatorialQ, 
{lmax0_,mm0_},
{
{Pm2j0_,Pm2j1_,s2jumps_},
{hj0_,hj1_,s0jumps_},
{dhlplp_Symbol,dhlmlm_Symbol,dhmpmp_Symbol,dhmmmm_Symbol,d\[Rho]chlpmm_Symbol,sourcetermsubs_}
},
{Bmatm2_,Bmatm1_,Bmat0_,Bmatp1_,Bmatp2_},
{eigs0_,eigs1_,eigs2_},
{\[Gamma]mat_},
opts:OptionsPattern[]]:=Module[
{
mm,ll,
r,\[Theta],l,m,a,\[Omega],a\[Omega]0,lmax,M,\[CapitalDelta],K,
a0,\[Omega]0,\[CapitalOmega]0,r0,
\[Lambda],\[Lambda]0,\[CapitalLambda],\[Gamma]ll,sgn,
CC,BB,AA,\[Alpha]sq,
Pp2,Pm2,Pp1,Pm1,h,h0,\[Kappa],Sm2,Sp2,Sm1,Sp1,S0,
h0subs,\[Kappa]subs,h0tohrepl,
teuks1subs,teuks2subs,
Pm1subs,Pp2subs,Pm2subs,
\[CapitalDelta]Ksimp,\[CapitalDelta]subs,Ksubs,
Dop,Ddag,
hvec,\[Kappa]mat,Pp1vec,Pm1vec,Pm2vec,Pp2vec,Pp2lrepl,Pp1lrepl,Pp2repl,Pp1repl,
ones,
knowns,unknowns,unknownstozero,sol,Mmat,brhs,eqs,
hj0lplp,hs0j0lplp,hs1j0lplp,hs2j0lplp,hj0mpmp,hs0j0mpmp,hs1j0mpmp,hs2j0mpmp,hj1lplp,hs0j1lplp,hs1j1lplp,hs2j1lplp,sourcehlplp,hj1mpmp,hs0j1mpmp,hs1j1mpmp,hs2j1mpmp,sourcehmpmp,hj0lmlm,hs0j0lmlm,hs1j0lmlm,hs2j0lmlm,
hj0mmmm,hs0j0mmmm,hs1j0mmmm,hs2j0mmmm,hj1lmlm,hs0j1lmlm,hs1j1lmlm,hs2j1lmlm,sourcehlmlm,hj1mmmm,hs0j1mmmm,hs1j1mmmm,hs2j1mmmm,sourcehmmmm,\[Rho]chj0lpmm,\[Rho]chs0j0lpmm,\[Rho]chs1j0lpmm,\[Rho]chs2j0lpmm, \[Rho]hj0lpmp,\[Rho]hs0j0lpmp,\[Rho]hs1j0lpmp,\[Rho]hs2j0lpmp,
lmins2,lmins1,lmins0,
\[Rho]chs2lpmm,\[Rho]chs1lpmm,\[Rho]chs0lpmm,
\[Rho]hs2lpmp,\[Rho]hs1lpmp,\[Rho]hs0lpmp,
hs2mmmm,hs2lmlm,hs1mmmm,hs1lmlm,hs0mmmm,hs0lmlm,
hs2mpmp,hs2lplp,hs1mpmp,hs1lplp,hs0mpmp,hs0lplp,
hlsubs,\[Kappa]lsubs,hjump,\[Kappa]jump,\[CapitalDelta]Kr0subs,awsubs,\[Gamma]llsubs,
Pm1jump,Pm2jump,Pm1lsubs,Pm2lsubs,Pp2lsubs,
DPm2vec,DDPm2vec,DDDPm2vec,
DdagPp2vec,DdagDdagPp2vec,DopDdagDdagPp2vec,
SSminus,SSplus,
HeadFn,
BBrepl,Crepl0,AArepl
},

HeadFn[head_,ll_]:={head[r]->head[ll][r],head'[r]->head[ll]'[r],head''[r]->head[ll]''[r],head'''[r]->head[ll]'''[r],head''''[r]->head[ll]''''[r]};

\[CapitalDelta]subs={\[CapitalDelta][r]->r^2-2*M*r+a^2,\[CapitalDelta]'[r]->2*(r-M),\[CapitalDelta]''[r]->2,\[CapitalDelta]'''[r]->0};
Ksubs={K[r]->\[Omega]*(r^2+a^2)-a*m,K'[r]-> 2*\[Omega]*r,K''[r]->2*\[Omega],K'''[r]->0};
\[CapitalDelta]Ksimp={\[CapitalDelta]'[r]->2*(r-M),\[CapitalDelta]''[r]->2,\[CapitalDelta]'''[r]->0,K'[r]-> 2*\[Omega]*r,K''[r]->2*\[Omega],K'''[r]->0};

Dop[n_,R_]:=D[R,r]-I*K[r]/\[CapitalDelta][r]*R+n*\[CapitalDelta]'[r]/\[CapitalDelta][r]*R;
Ddag[n_,R_]:=D[R,r]+I*K[r]/\[CapitalDelta][r]*R+n*\[CapitalDelta]'[r]/\[CapitalDelta][r]*R;
Dop[R_]:=Dop[0,R];
Ddag[R_]:=Ddag[0,R];

(* Teuksolky equations for spin 2 *)
teuks2subs=DefTeukSubs[2,{Pm2,Pp2,Sm2,Sp2},{\[CapitalDelta],K},{M,a,m,\[Omega],\[CapitalLambda]},{r,\[Theta]}];
Pm2subs=teuks2subs[[1;;3]];
Pp2subs=teuks2subs[[4;;6]];

(* Teukolsky-Starobinskii identities for spin 2 *)
Pp2repl=DefTSIP[2,"pm",{Pm2,Pp2},{\[CapitalDelta],K},{M,a,m,\[Omega],\[CapitalLambda],CC},{r}];

(* Teukolsky equations for spin 1. *)
teuks1subs=DefTeukSubs[1,{Pm1,Pp1,Sm1,Sp1},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda]},{r,\[Theta]}];
Pm1subs=teuks1subs[[1;;3]];

(* Teukolsky-Starobinskii identities for spin 1. *)
Pp1repl=DefTSIP[1,"pm",{Pm1,Pp1},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda],BB},{r}];

(* Teukolsky equations for trace, and for kappa. *)
h0subs=DefTeukSubs[0,{h0,S0},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda]0},{r,\[Theta]}][[1;;3]];

a0=ODspin[OD];
If[PossibleZeroQ[a0],a0=0]; (* In the Schwarzschild case, a0 should be the integer zero. *)
r0=ODsemilatusrectum[OD];
lmax=lmax0;
mm=mm0;


\[CapitalOmega]0=KerrAzimuthalFrequency[OD];
\[Omega]0=mm \[CapitalOmega]0;
a\[Omega]0=If[a0==0,0,a0*\[Omega]0];
\[CapitalDelta]Kr0subs=Map[#[[1]]->(#[[2]]/.{a->a0,r->r0,M->1,m->mm,\[Omega]->\[Omega]0})&,Join[\[CapitalDelta]subs,Ksubs]];
awsubs=Join[{a->a0,\[Omega]->\[Omega]0},{m->mm,M->1}];
lmins0=Abs[mm];
lmins1=Max[1,Abs[mm]];
lmins2=Max[2,Abs[mm]];

(* Radial functions as vectors. *)
ones=Table[1,{ll,0,lmax}];
Pp2vec=Table[If[ll>=lmins2,Pp2[ll][r],0],{ll,0,lmax}];
Pp1vec=Table[If[ll>=lmins1,Pp1[ll][r],0],{ll,0,lmax}];
hvec=Table[If[ll>=lmins0,h[ll][r],0],{ll,0,lmax}];
Pm1vec=Table[If[ll>=lmins1,Pm1[ll][r],0],{ll,0,lmax}];
Pm2vec=Table[If[ll>=lmins2,Pm2[ll][r],0],{ll,0,lmax}];

(* \[Kappa] needs special handling *)
Module[{tmp1,tmp2},
tmp1=DiagonalMatrix[Table[If[ll>=lmins0,\[Kappa][ll][r],0],{ll,0,lmax}]];
tmp2=a^2*Table[
If[j!=k,
If[j>=lmins0&&k>=lmins0,
\[Gamma]mat[[j+1,k+1]]*h[j][r]/(2*(eigs0[[j+1]]-eigs0[[k+1]]))
,
0
],
0],
{j,0,lmax},{k,0,lmax}];
\[Kappa]mat=(tmp1+tmp2)/.awsubs;
];

(* Awkwardness: to keep the expressions compact, I would like to replace 'r' with the numerical value 'r0', but not in the radial functions themselves, as these I will use as unknowns. So I will replace the radial functions with jump expressions at an early stage. *)

Pm2jump=Table[If[ll>=lmins2,{Pm2[ll][r]->Pm2j0[ll],Pm2[ll]'[r]->Pm2j1[ll]},{}],{ll,0,lmax}]//Flatten;

Pm1jump=Table[If[ll>=lmins1,{Pm1[ll][r]->Pm1j0[ll],Pm1[ll]'[r]->Pm1j1[ll]},{}],{ll,0,lmax}]//Flatten;
hjump=Table[If[ll>=lmins0,{h[ll][r]->0,h[ll]'[r]->hj1[ll]},{}],{ll,0,lmax}]//Flatten;
\[Kappa]jump=Table[If[ll>=lmins0,{\[Kappa][ll][r]->\[Kappa]j0[ll],\[Kappa][ll]'[r]->\[Kappa]j1[ll]},{}],{ll,0,lmax}]//Flatten;



\[Kappa]subs=DefTeukSubs["\[Kappa]",{\[Kappa],h},{\[CapitalDelta],K},{M,a,m,\[Omega],\[Lambda]0,\[Gamma]ll[l]},{r}];
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
(*Pp1lsubs= Table[Pp1subs/.awsubs/.{\[Lambda]->eigs1[[ll+1]]}/.HeadFn[Pp1,ll],{ll,lmins1,lmax}]//Flatten;*)
hlsubs= Table[h0subs/.awsubs/.{\[Lambda]0->eigs0[[ll+1]]}/.h0tohrepl/.HeadFn[h,ll],{ll,lmins0,lmax}]//Flatten;
\[Kappa]lsubs= Table[\[Kappa]subs/.awsubs/.{\[Lambda]0->eigs0[[ll+1]]}/.{l->ll}/.h0tohrepl/.HeadFn[\[Kappa],ll]/.HeadFn[h,ll],{ll,lmins0,lmax}]//Flatten;

(* Now consider the Teuksolsky-Starobinskii identities *)
AArepl={AA->Sqrt[\[CapitalLambda]^2*(\[CapitalLambda]+2)^2-8*\[Omega]^2*\[CapitalLambda]*(\[Alpha]sq*(5*\[CapitalLambda]+6)-12*a^2)+144*\[Omega]^4*\[Alpha]sq^2]}/.{\[Alpha]sq->a^2-a*m/\[Omega]};
(* Check against Pound and Wardell. *)
BBrepl={BB->Sqrt[\[Lambda]^2 +4*a*m*\[Omega]-4*a^2*\[Omega]^2 ]};
Crepl0={CC->AA+12 I M sgn \[Omega],Cstar->AA-12 I M sgn \[Omega]};
Clear[ll];


Pp2lrepl=Table[Pp2repl/.awsubs/.HeadFn[Pm2,ll]/.HeadFn[Pp2,ll]/.Crepl0/.AArepl/.{\[CapitalLambda]->eigs2[[ll+1]]}/.{sgn->(-1)^(ll+mm)}/.awsubs,{ll,lmins2,lmax}]//Flatten;
Pp1lrepl=Table[Pp1repl/.awsubs/.HeadFn[Pm1,ll]/.HeadFn[Pp1,ll]/.BBrepl/.{\[Lambda]->eigs1[[ll+1]]}/.awsubs,{ll,lmins1,lmax}]//Flatten;

(* Spin 1 projections. *)
(* I have deleted all the projections except lplp, mpmp, lmlm, mmmm and rho h_{lpmp} and rhoc h_{lpmm} *)
EchoTiming[
Module[{line1,line2,tmp,tmp1,tmp2,SS0,calPs1,DcalPs1,sinmix,cosmix,hraws1mmmm,hraws1mpmp,hraws1lmlm,hraws1lplp,Ss1mmmm,Ss1mpmp,Ss1lmlm,Ss1lplp,t0,t1},
t0=-(Bmatm1+signmat[lmax,mm] . Bmatp1) . \[Lambda]hat[lmax,mm];
t1=a\[Omega]0*(Bmatm1 . SinMixMatrix[-1,0,lmax,mm]+signmat[lmax,mm] . Bmatp1 . SinMixMatrix[1,0,lmax,mm]);
SS0=t0+t1;
Ss1lplp=SS0;
Ss1lmlm=SS0;
Ss1mpmp=Bmatp1 . (-\[Lambda]2hat[lmax,mm]+a\[Omega]0*SinMixMatrix[1,2,lmax,mm]);
(*Ss1mmmm=Bmatm1 . (\[Lambda]2hat[lmax,mm]-a\[Omega]0*SinMixMatrix[-1,-2,lmax,mm]);*)
tmp=2/\[CapitalDelta][r]*Table[Dop[-1,Pp1vec[[ll+1]]],{ll,0,lmax}];
hraws1lplp=tmp . signmat[lmax,mm] . Ss1lplp;  (* zzz edited line *)
hs1lplp=hraws1lplp/.Pp1lrepl/.awsubs;
tmp=2/\[CapitalDelta][r]*Table[Ddag[-1,Pm1vec[[ll+1]]],{ll,0,lmax}];
hraws1lmlm=tmp . Ss1lmlm;
hs1lmlm=hraws1lmlm/.awsubs;
tmp1=signmat[lmax,mm] . Table[Ddag[Pp1vec[[ll+1]]],{ll,0,lmax}]/.Pp1lrepl/.awsubs;
tmp2=Table[Dop[Pm1vec[[ll+1]]],{ll,0,lmax}]/.awsubs;
tmp=2*(tmp1+tmp2);
hraws1mpmp=tmp . signmat[lmax,mm] . Ss1mpmp;  
(*hraws1mmmm=-tmp . Ss1mmmm;*)
hs1mpmp=hraws1mpmp/.Pp1lrepl/.awsubs;
(*hs1mmmm=hraws1mmmm/.Pp1lrepl/.awsubs;*)
(* Now for the other components *)
(*tmp1=signmat[lmax,mm] . Table[Ddag[Pp1vec[[ll+1]]],{ll,0,lmax}]/.Pp1lrepl/.awsubs;
tmp2=Table[Dop[Pm1vec[[ll+1]]],{ll,0,lmax}]/.awsubs;
calPs1=tmp2+tmp1;
DcalPs1=Map[Dop,calPs1]/.Pm1lsubs;*)
(* Ss1lplp//TableForm *)
(* l+ m+ *)
cosmix=CosMixMatrix[1,1,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
(*line1=(signmat[lmax,mm] . (r*DcalPs1-2*calPs1)) . Bmatp1+I*a*(signmat[lmax,mm] . DcalPs1) . Bmatp1 . cosmix/.Pm1lsubs;
tmp1=1/\[CapitalDelta][r]*signmat[lmax,mm] . Pp1vec/.Pp1lrepl/.awsubs;
tmp2=-r*\[Lambda]hat[lmax,mm]-I*a*\[Lambda]hat[lmax,mm] . cosmix+(a*\[Omega]*r+2*I*a)*sinmix+I*a^2*\[Omega]*sinmix . cosmix;
line2=tmp1 . SS0 . tmp2;
\[Rho]hs1lpmp=(line1+line2)/.awsubs;*)
(* Now the \[Rho]c h_{l+ m-} term *)
(*line1=-(r*DcalPs1-2*calPs1) . Bmatm1+I*a*DcalPs1 . Bmatm1 . CosMixMatrix[-1,1,lmax,mm]/.Pm1lsubs;
tmp1=1/\[CapitalDelta][r]*signmat[lmax,mm] . Pp1vec/.Pp1lrepl/.awsubs;
tmp2=r*\[Lambda]hat[lmax,mm]-I*a*\[Lambda]hat[lmax,mm] . CosMixMatrix[-1,1,lmax,mm]-(a*\[Omega]*r+2*I*a)*SinMixMatrix[0,-1,lmax,mm]+I*a^2*\[Omega]*SinMixMatrix[0,-1,lmax,mm] . CosMixMatrix[-1,1,lmax,mm];
line2=tmp1 . SS0 . tmp2;
\[Rho]chs1lpmm=(line1+line2)/.awsubs;*)
]
];
(* Spin 2: lplp and mpmp *)
EchoTiming[
Module[{hraws2mmmm,hraws2mpmp,hraws2lmlm,hraws2lplp,tmp,tmp1,tmp2,Ss2mmmm,Ss2mpmp,Ss2lmlm,Ss2lplp ,t0,t1},
t0=Bmatp2 . (\[CapitalLambda]hat[lmax,mm]-2*a\[Omega]0*\[Lambda]2hat[lmax,mm] . SinMixMatrix[1,0,lmax,mm]+a\[Omega]0^2*SinMixMatrix[2,0,lmax,mm]);
t1=signmat[lmax,mm] . Bmatm2 . (\[CapitalLambda]hat[lmax,mm]-2*a\[Omega]0*\[Lambda]2hat[lmax,mm] . SinMixMatrix[-1,0,lmax,mm]+a\[Omega]0^2*SinMixMatrix[-2,0,lmax,mm]);
Ss2lplp=t0+t1;
Ss2lmlm=signmat[lmax,mm] . Ss2lplp; 
Ss2mpmp=Bmatp2;
(*Ss2mmmm=signmat[lmax,mm] . Bmatm2;*)
hraws2lplp=-1/(6*\[Omega]^2*\[CapitalDelta][r]^2)*Pp2vec . Ss2lplp;
hraws2lmlm=-1/(6*\[Omega]^2*\[CapitalDelta][r]^2)*Pm2vec . Ss2lmlm;
hs2lplp=hraws2lplp/.awsubs/.Pp2lrepl;
hs2lmlm=hraws2lmlm/.awsubs;
tmp1=Table[Ddag[Ddag[Pp2vec[[ll+1]]]],{ll,0,lmax}];
tmp2=signmat[lmax,mm] . Table[Dop[Dop[Pm2vec[[ll+1]]]],{ll,0,lmax}];
tmp=-1/(6*\[Omega]^2)*(tmp1+tmp2);
hraws2mpmp=tmp . Ss2mpmp;
(*hraws2mmmm=tmp . Ss2mmmm;*)
hs2mpmp=hraws2mpmp/.Pp2lsubs/.Pp2lrepl/.Pm2lsubs/.awsubs;
(*hs2mmmm=hraws2mmmm/.Pp2lsubs/.Pp2lrepl/.Pm2lsubs/.awsubs;*)
]
];
EchoTiming[
SSplus=(Bmatp2 . (\[CapitalLambda]hat[lmax,mm]-2*a*\[Omega]*\[Lambda]2hat[lmax,mm] . SinMixMatrix[1,0,lmax,mm]+a^2*\[Omega]^2*SinMixMatrix[2,0,lmax,mm]))/.awsubs//Simplify;
SSminus=(Bmatm2 . (\[CapitalLambda]hat[lmax,mm]-2*a*\[Omega]*\[Lambda]2hat[lmax,mm] . SinMixMatrix[-1,0,lmax,mm]+a^2*\[Omega]^2*SinMixMatrix[-2,0,lmax,mm]))/.awsubs//Simplify;
DdagPp2vec=Map[Ddag,Pp2vec]/.Pp2lrepl/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DdagDdagPp2vec=Map[Ddag,DdagPp2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
(*DdagDdagDdagPp2vec=Map[Ddag,DdagDdagPp2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;*)
DPm2vec=Map[Dop,Pm2vec]/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DDPm2vec=Map[Dop,DPm2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DDDPm2vec=Map[Dop,DDPm2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
DopDdagDdagPp2vec=Map[Dop,DdagDdagPp2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;
(*DdagDDPm2vec=Map[Ddag,DDPm2vec]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs//Simplify;*)
];
(* lp mp *)
(*EchoTiming[
Module[{term1,term2,term3,term4,tmp,MMp2,MMm2,\[Rho]cmat,sinmix1m1,sinmix,cosmix,\[Rho]mat},
cosmix=CosMixMatrix[1,1,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
sinmix1m1=SinMixMatrix[-1,1,lmax,mm];
\[Rho]mat=r*idmat[lmax,mm]+I*a*cosmix;
MMp2=Bmatp2 . (\[Lambda]2hat[lmax,mm]-a*\[Omega]*SinMixMatrix[2,1,lmax,mm])/.awsubs;
MMm2=Bmatm2 . (\[Lambda]2hat[lmax,mm]-a*\[Omega]*SinMixMatrix[-2,-1,lmax,mm])/.awsubs;
term1=1/(12*\[Omega]^2)*(DDPm2vec . signmat[lmax,mm] . SSplus-DdagDdagPp2vec . signmat[lmax,mm] . SSminus) . (\[Lambda]hat[lmax,mm]-a*\[Omega]*sinmix)/.awsubs;
term2=-1/(12*\[Omega]^2)*(DDDPm2vec . signmat[lmax,mm] . SSplus-DopDdagDdagPp2vec . signmat[lmax,mm] . SSminus) . ((\[Lambda]hat[lmax,mm]-a*\[Omega]*sinmix) . \[Rho]mat-I*a*sinmix)/.awsubs;
term3=-1/(3*I*\[Omega])*(DDDPm2vec . signmat[lmax,mm] . MMp2 . \[Rho]mat . \[Rho]mat-2*DDPm2vec . signmat[lmax,mm] . MMp2 . \[Rho]mat+2*DPm2vec . signmat[lmax,mm] . MMp2)/.awsubs;
tmp=(\[Lambda]hat[lmax,mm] . \[Lambda]hat[lmax,mm]-2*a*\[Omega]*\[Lambda]hat[lmax,mm] . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]mat . \[Rho]mat-2*I*a*(\[Lambda]hat[lmax,mm] . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]mat-2*a^2*sinmix1m1;
term4=-1/(3*I*\[Omega]*\[CapitalDelta][r])*DdagPp2vec . (signmat[lmax,mm] . MMm2 . tmp)/.awsubs;
(* Don't forget the \[Beta] factor. *)
tmp=term1+term2+term3+term4;
\[Rho]hs2lpmp=(tmp/(-6*I*M*\[Omega]))/.awsubs;
]
];*)
(* lp mm *)
(*EchoTiming[
Module[{term1,term2,term3,term4,tmp,MMp2,MMm2,\[Rho]cmat,sinmix1m1,sinmix},
\[Rho]cmat=r*idmat[lmax,mm]-I*a*CosMixMatrix[-1,1,lmax,mm];
sinmix=SinMixMatrix[0,-1,lmax,mm];
MMp2=Bmatp2 . (\[Lambda]2hat[lmax,mm]-a*\[Omega]*SinMixMatrix[2,1,lmax,mm]);
MMm2=Bmatm2 . (\[Lambda]2hat[lmax,mm]-a*\[Omega]*SinMixMatrix[-2,-1,lmax,mm]);
term1=-1/(12*\[Omega]^2)*(DDPm2vec . SSminus-DdagDdagPp2vec . SSplus) . (\[Lambda]hat[lmax,mm]-a*\[Omega]*sinmix);
term2=1/(12*\[Omega]^2)*(DDDPm2vec . SSminus-DopDdagDdagPp2vec . SSplus) . ((\[Lambda]hat[lmax,mm]-a*\[Omega]*sinmix) . \[Rho]cmat-I*a*sinmix);
term3=1/(3*I*\[Omega])*(DDDPm2vec . MMm2 . \[Rho]cmat . \[Rho]cmat-2*DDPm2vec . MMm2 . \[Rho]cmat+2*DPm2vec . MMm2);
sinmix1m1=SinMixMatrix[1,-1,lmax,mm];
tmp=(\[Lambda]hat[lmax,mm] . \[Lambda]hat[lmax,mm]-2*a*\[Omega]*\[Lambda]hat[lmax,mm] . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]cmat . \[Rho]cmat-2*I*a*(\[Lambda]hat[lmax,mm] . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]cmat-2*a^2*sinmix1m1;
term4=1/(3*I*\[Omega]*\[CapitalDelta][r])*DdagPp2vec . MMp2 . tmp;
tmp=term1+term2+term3+term4;
\[Rho]chs2lpmm=(tmp/(-6*I*M*\[Omega]))/.awsubs;
]
];*)

(* Spin 0 *)
(* When we apply \[Kappa]lsubs, we also should fill in the \[Gamma]ll values *)
EchoTiming[
Module[{tmp1,tmp2a,tmp2,Dhvec,D\[Kappa]mat,MM0,sinmix,cosmix,t0,t1,tmpS\[Kappa],tmpSh},
tmp1=-1/(I*\[Omega])*Dop[r*Dop[hvec]] . Bmat0/.hlsubs/.awsubs;
tmp2a=(ones . Dop[Dop[\[Kappa]mat]]) . Bmat0/.\[Kappa]lsubs/.\[Gamma]llsubs/.hlsubs/.awsubs;
tmp2=-4*tmp2a/.awsubs;
hs0lplp=(tmp1+tmp2)/.awsubs;

tmp1=1/(I*\[Omega])*Ddag[r*Ddag[hvec]] . Bmat0/.hlsubs/.awsubs;
tmp2a=(ones . Ddag[Ddag[\[Kappa]mat]]) . Bmat0/.\[Kappa]lsubs/.\[Gamma]llsubs/.hlsubs/.awsubs;
tmp2=-4*tmp2a/.awsubs;
hs0lmlm=(tmp1+tmp2)/.awsubs;

t0=Bmat0 . (\[CapitalLambda]hat[lmax,mm]-2*a\[Omega]0*\[Lambda]hat[lmax,mm] . SinMixMatrix[1,2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,2,lmax,mm]) . CosMixMatrix[2,1,lmax,mm];
t1=Bmat0 . (\[Lambda]hat[lmax,mm] . SinMixMatrix[1,2,lmax,mm]-a\[Omega]0*SinMixMatrix[0,2,lmax,mm]);
tmpSh=(t0+t1)/.awsubs;
tmpS\[Kappa]=Bmat0 . (\[CapitalLambda]hat[lmax,mm]-2*a\[Omega]0*\[Lambda]hat[lmax,mm] . SinMixMatrix[1,2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,2,lmax,mm])/.awsubs;
hs0mpmp=((a/\[Omega])*hvec . tmpSh-4*(ones . \[Kappa]mat) . tmpS\[Kappa])/.awsubs;

(*t0=Bmat0 . (\[CapitalLambda]hat[lmax,mm]-2*a\[Omega]0*\[Lambda]hat[lmax,mm] . SinMixMatrix[-1,-2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,-2,lmax,mm]) . CosMixMatrix[-2,1,lmax,mm];
t1=-Bmat0 . (\[Lambda]hat[lmax,mm] . SinMixMatrix[-1,-2,lmax,mm]-a\[Omega]0*SinMixMatrix[0,-2,lmax,mm]);
tmpSh=(t0+t1)/.awsubs;
tmpS\[Kappa]=Bmat0 . (\[CapitalLambda]hat[lmax,mm]-2*a\[Omega]0*\[Lambda]hat[lmax,mm] . SinMixMatrix[-1,-2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,-2,lmax,mm])/.awsubs;
hs0mmmm=(-(a/\[Omega])*hvec . tmpSh-4*(ones . \[Kappa]mat) . tmpS\[Kappa])/.awsubs;*)
]
];

(* lp mp *)
(*EchoTiming[
Module[{term1,term2,term3,term4,tmp,Dhvec,D\[Kappa]mat,MM0,sinmix,cosmix,t0},
t0=r^2*idmat[lmax,mm]+a^2*CosMixMatrix[1,2,lmax,mm];
sinmix=SinMixMatrix[0,1,lmax,mm];
cosmix=CosMixMatrix[1,1,lmax,mm];
MM0=\[Lambda]hat[lmax,mm]-a*\[Omega]*sinmix;
Dhvec=Map[Dop,hvec];
D\[Kappa]mat=Map[Dop,\[Kappa]mat];
term1=1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0;
term2=-a/\[Omega]*(r*Dhvec . Bmat0 . sinmix-hvec . Bmat0 . MM0 . cosmix);
tmp=r*idmat[lmax,mm]+I*a*cosmix;
term3=4*(ones . D\[Kappa]mat) . Bmat0 . MM0 . tmp;
term4=-4*((ones . \[Kappa]mat) . Bmat0 . MM0+I*a*(ones . D\[Kappa]mat) . Bmat0 . sinmix);
\[Rho]hs0lpmp=(term1+term2+term3+term4)/.awsubs;
]
];*)

(* lp mm *)
(*EchoTiming[
Module[{term1,term2,term3,term4,tmp,Dhvec,D\[Kappa]mat,MM0,sinmix,t0},
t0=r^2*idmat[lmax,mm]+a^2*CosMixMatrix[-1,2,lmax,mm];
sinmix=SinMixMatrix[0,-1,lmax,mm];
MM0=\[Lambda]hat[lmax,mm]-a*\[Omega]*sinmix;
Dhvec=Map[Dop,hvec];
D\[Kappa]mat=Map[Dop,\[Kappa]mat];
term1=-1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0;
term2=a/\[Omega]*(r*Dhvec . Bmat0 . sinmix+hvec . Bmat0 . MM0 . CosMixMatrix[-1,1,lmax,mm]);
tmp=r*idmat[lmax,mm]-I*a*CosMixMatrix[-1,1,lmax,mm];
term3=-4*((ones . D\[Kappa]mat) . Bmat0 . MM0 . tmp);
term4=4*((ones . \[Kappa]mat) . Bmat0 . MM0+I*a*(ones . D\[Kappa]mat) . Bmat0 . sinmix);
\[Rho]chs0lpmm=(term1+term2+term3+term4)/.awsubs;
]
];*)

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
(*hs0j0mmmm=hs0mmmm/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
hs0j1mmmm=D[hs0mmmm,r]/.hlsubs/.\[Kappa]lsubs/.\[Gamma]llsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;*)
hs1j0lmlm=hs1lmlm/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs1j1lmlm=D[hs1lmlm,r]/.Pm1lsubs/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
(*hs1j0mmmm=hs1mmmm/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs1j1mmmm=D[hs1mmmm,r]/.Pm1lsubs/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;*)
hs2j0lmlm=hs2lmlm/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs2j1lmlm=D[hs2lmlm,r]/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
(*hs2j0mmmm=hs2mmmm/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
hs2j1mmmm=D[hs2mmmm,r]/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;*)
];

(* \[Rho]lpmp and \[Rho]clpmm. *)
(*EchoTiming[
(*\[Rho]hs0j0lpmp=(\[Rho]hs0lpmp)/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
\[Rho]hs1j0lpmp=(\[Rho]hs1lpmp)/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
\[Rho]hs2j0lpmp=(\[Rho]hs2lpmp)/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;*)
];
*)
(*EchoTiming[
(*\[Rho]chs0j0lpmm=(\[Rho]chs0lpmm)/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
\[Rho]chs1j0lpmm=(\[Rho]chs1lpmm)/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
\[Rho]chs2j0lpmm=(\[Rho]chs2lpmm)/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;*)
];
*)
sourcehlplp=Table[dhlplp[ll],{ll,0,lmax}]/.sourcetermsubs;
sourcehlmlm=Table[dhlmlm[ll],{ll,0,lmax}]/.sourcetermsubs;
sourcehmpmp=Table[dhmpmp[ll],{ll,0,lmax}]/.sourcetermsubs;
(*sourcehmmmm=Table[dhmmmm[ll],{ll,0,lmax}]/.sourcetermsubs;*)
(*source\[Rho]chlpmm=Table[d\[Rho]chlpmm[ll],{ll,0,lmax}]/.sourcetermsubs;*)

(* OK, let's check out whether the linear system of equations has a solution, in principle. *)
(* We are solving for the unknowns: *)
unknowns=Join[Table[Pm1j0[ll],{ll,lmins1,lmax}],Table[Pm1j1[ll],{ll,lmins1,lmax}],Table[\[Kappa]j0[ll],{ll,lmins0,lmax}],Table[\[Kappa]j1[ll],{ll,lmins0,lmax}]];
unknownstozero=Table[unknowns[[kk]]->0,{kk,1,Length[unknowns]}];

(* We have the system of equations: *)
hj0lplp=hs0j0lplp+hs1j0lplp+hs2j0lplp;
hj1lplp=(hs0j1lplp+hs1j1lplp+hs2j1lplp)-sourcehlplp;

hj0mpmp=hs0j0mpmp+hs1j0mpmp+hs2j0mpmp;
hj1mpmp=(hs0j1mpmp+hs1j1mpmp+hs2j1mpmp)-sourcehmpmp;

hj0lmlm=hs0j0lmlm+hs1j0lmlm+hs2j0lmlm;
hj1lmlm=(hs0j1lmlm+hs1j1lmlm+hs2j1lmlm)-sourcehlmlm;

(*hj0mmmm=hs0j0mmmm+hs1j0mmmm+hs2j0mmmm; (*unnecessary?*)
hj1mmmm=(hs0j1mmmm+hs1j1mmmm+hs2j1mmmm)-sourcehmmmm; (*unnecessary?*)*)

(*\[Rho]chj0lpmm=\[Rho]chs0j0lpmm+\[Rho]chs1j0lpmm+\[Rho]chs2j0lpmm; (*unnecessary?*)
\[Rho]hj0lpmp=\[Rho]hs0j0lpmp+\[Rho]hs1j0lpmp+\[Rho]hs2j0lpmp; (*unnecessary?*)*)


eqs=Join[
	Table[hj0lplp[[ll+1]],{ll,lmins0,lmax}],
	Table[hj1lplp[[ll+1]],{ll,lmins0,lmax}],
	Table[hj0mpmp[[ll+1]],{ll,lmins2,lmax}],
	Table[hj1mpmp[[ll+1]],{ll,lmins2,lmax}]
	]/.s2jumps/.s0jumps;
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
EchoTiming[Quiet[sol=LinearSolve[Mmat,brhs]],"LinearSolve for s1 and \[Kappa] jumps"];
knowns=Table[unknowns[[k]]->CleanValue@sol[[k]],{k,1,Length[unknowns]}];
Return[knowns]
]


(* ::Subsection:: *)
(*h\[Kappa] Solver*)


(* ::Subsubsection::Closed:: *)
(*inf series*)


(* ::Input::Initialization:: *)
h\[Kappa]infSeries[m_,a_,\[Omega]_,\[Lambda]0_,\[Gamma]ll_,r_, nmax_:10]:=
Module[{h,\[Kappa],\[Alpha],\[Beta],\[Gamma],\[Delta],prefac,t1},
h[n_]:=(h[n]=-(1/(\[Beta][2]n))Sum[h[i](i(i-1)\[Alpha][n+3-i]+i \[Beta][n+2-i]+\[Gamma][n+1-i]),{i,Max[0,n-5],n-1}]);
\[Kappa][n_]:=(\[Kappa][n]=-(1/(\[Beta][2]n))(
Sum[\[Kappa][i](i(i-1)\[Alpha][n+3-i]+i \[Beta][n+2-i]+\[Gamma][n+1-i]),{i,Max[-1,n-5],n-1}]+
+Sum[h[n+1-j]\[Delta][j],{j,0,Min[n+1,4]}])
);
h[0]=1;
h[1]= (I (\[Lambda]0+2 (a m-4 \[Omega]) \[Omega]))/(2 \[Omega]);
\[Kappa][-1]= 1/(4I \[Omega]);
\[Kappa][0]=0;
\[Kappa][1]=(2 (a m-4 \[Omega]+a^2 \[Omega])+I (a^2 \[Gamma]ll+h[2]))/(4 \[Omega]);
\[Alpha][0]=0;
\[Alpha][1]=0;
\[Alpha][2]=0;
\[Alpha][3]=0;
\[Alpha][4]=1;
\[Alpha][5]=-4;
\[Alpha][6]=2 (2+a^2);
\[Alpha][7]=-4 a^2;
\[Alpha][8]=a^4;
\[Alpha][n_]/;n>8:=0;
\[Beta][0]=0;
\[Beta][1]=0;
\[Beta][2]=-2 I \[Omega];
\[Beta][3]=2+4 I \[Omega];
\[Beta][4]=-10-4 I (-2+a^2) \[Omega];
\[Beta][5]=2 (6+3 a^2-8 I \[Omega]);
\[Beta][6]=-2 I a^2 (-7 I+(-8+a^2) \[Omega]);
\[Beta][7]=4 a^4 (1-I \[Omega]);
\[Beta][n_]/;n>7:=0;
\[Gamma][-1]=0;
\[Gamma][0]=0;
\[Gamma][1]=0;
\[Gamma][2]=-\[Lambda]0-2 a m \[Omega]+8 \[Omega]^2;
\[Gamma][3]=2 (-1+\[Lambda]0-I (-4+a^2) \[Omega]-2 a^2 \[Omega]^2);
\[Gamma][4]=-2 a^3 m \[Omega]-4 (I+2 \[Omega])^2+a^2 (2+m^2-\[Lambda]0-2 I \[Omega]+8 \[Omega]^2);
\[Gamma][5]=-2 a^2 (I+2 \[Omega]) (-3 I+(-4+a^2) \[Omega]);
\[Gamma][6]=-2 a^4 (-1+3 I \[Omega]+2 \[Omega]^2);
\[Gamma][n_]/;n>6:=0;
\[Delta][0]=-1/2;
\[Delta][1]=1;
\[Delta][2]=-(1/2) a^2 (1+\[Gamma]ll);
\[Delta][3]=a^2 \[Gamma]ll;
\[Delta][4]=-((a^4 \[Gamma]ll)/2);
\[Delta][n_]/;n>4:=0;
prefac=2^(-2 I \[Omega]) E^(I r \[Omega]-(1-2 I \[Omega]) Log[r]);
Return[{
{t1=prefac Sum[h[n]r^-n,{n,0,nmax}], prefac Sum[-n h[n]r^(-n-1),{n,1,nmax}]+t1(-((1-2 I \[Omega])/r)+I \[Omega]), Abs[prefac h[nmax]/r^nmax]},
{t1=prefac Sum[\[Kappa][n]r^-n,{n,-1,nmax}], prefac Sum[-n \[Kappa][n]r^(-n-1),{n,-1,nmax}]+t1(-((1-2 I \[Omega])/r)+I \[Omega]), Abs[prefac \[Kappa][nmax]/r^nmax]}
}]
]


(* ::Subsubsection::Closed:: *)
(*hor series*)


(* ::Input::Initialization:: *)
h\[Kappa]horSeries[m_,a_,\[Omega]_,\[Lambda]0_,\[Gamma]ll_,r_, nmax_:10]:=
Module[{h,\[Kappa],\[Alpha],\[Beta],\[Gamma],\[Delta],prefac,t1,x,\[Nu]},
\[Nu]=(I (a m-2 (1+Sqrt[1-a^2]) \[Omega]))/(2 Sqrt[1-a^2]);
x=r-(1+Sqrt[1-a^2]);h[n_]:=(h[n]=-(1/(\[Alpha][2]n(n-1)+\[Beta][1]n+\[Gamma][0]))Sum[h[i](i(i-1)\[Alpha][n+2-i]+i \[Beta][n+1-i]+\[Gamma][n-i]),{i,Max[0,n-4],n-1}]);
	\[Kappa][n_]:=(\[Kappa][n]=-(1/(\[Alpha][2]n(n-1)+\[Beta][1]n+\[Gamma][0]))(
	Sum[\[Kappa][i](i(i-1)\[Alpha][n+2-i]+i \[Beta][n+1-i]+\[Gamma][n-i]),{i,Max[1,n-4],n-1}]+
	+Sum[h[n-j]\[Delta][j],{j,0,Min[n,4]}])
	);
h[0]=1;
\[Alpha][0]=0;
\[Alpha][1]=0;
\[Alpha][2]=4-4 a^2;
\[Alpha][3]=4 Sqrt[1-a^2];
\[Alpha][4]=1;
\[Alpha][n_]/;n>4:=0;
\[Beta][0]=0;
\[Beta][1]=-4 (-1+a^2) (1+2 \[Nu]);
\[Beta][2]=Sqrt[1-a^2] (6+8 \[Nu]);
\[Beta][3]=2 (1+\[Nu]);
\[Beta][n_]/;n>3:=0;
\[Gamma][-1]=0;
\[Gamma][0]=-4 a m \[Omega]-4 Sqrt[1-a^2] (a m-2 \[Omega]) \[Omega]+4 (\[Nu]^2+2 \[Omega]^2)+a^2 (m^2-4 (\[Nu]^2+\[Omega]^2));
\[Gamma][1]=-4 \[Omega] (a m-4 \[Omega]+2 a^2 \[Omega])+2 Sqrt[1-a^2] (-\[Lambda]0+\[Nu]+2 \[Nu]^2-2 a m \[Omega]+8 \[Omega]^2);
\[Gamma][2]=-\[Lambda]0+\[Nu]+\[Nu]^2+2 \[Omega] (-a m-2 a^2 \[Omega]+6 (1+Sqrt[1-a^2]) \[Omega]);
\[Gamma][3]=4 (1+Sqrt[1-a^2]) \[Omega]^2;
\[Gamma][4]=\[Omega]^2;
\[Gamma][n_]/;n>4:=0;
\[Delta][0]=0;
\[Delta][1]=2 (-1+a^2)+Sqrt[1-a^2] (-2-a^2 (-1+\[Gamma]ll));
\[Delta][2]=-3 (1+Sqrt[1-a^2])-1/2 a^2 (-5+\[Gamma]ll);
\[Delta][3]=-1-2 Sqrt[1-a^2];
\[Delta][4]=-(1/2);
\[Delta][n_]/;n>4:=0;
prefac=E^(-I (-((a m)/(2+2 Sqrt[1-a^2]))+\[Omega]) (1+Sqrt[1-a^2]-Log[4]+1/2 (1-1/Sqrt[1-a^2]) Log[4-4 a^2])) x^\[Nu];
Return[{
{t1=prefac Sum[h[n]x^n,{n,0,nmax}], prefac Sum[n h[n]x^(n-1),{n,1,nmax}]+t1(\[Nu]/x), Abs[prefac h[nmax]x^nmax]},
{t1=prefac Sum[\[Kappa][n]x^n,{n,1,nmax}], prefac Sum[n \[Kappa][n]x^(n-1),{n,1,nmax}]+t1(\[Nu]/x), Abs[prefac \[Kappa][nmax]x^nmax]}
}]
]


(* ::Subsubsection::Closed:: *)
(*Evolution Series*)


(* ::Input::Initialization:: *)
h\[Kappa]Ds[m_,a_,\[Omega]_,\[Lambda]0_, \[Gamma]ll_,r_,{h0_,h1_},{\[Kappa]0_,\[Kappa]1_},nmax_:10]:=
Module[{h,\[Kappa],\[Alpha],\[Beta],\[Gamma],\[Delta],prec=Max[$MinPrecision,Min[$TeukolskyWorkingPrecision,Quiet@Precision[{m,a,\[Omega],\[Lambda]0,r,h0,h1}]]]},
h[n_]:=(h[n]=SetPrecision[-(1/(\[Alpha][0]n(n-1)))(
If[n>2,h[n-1]((n-1)(n-2)\[Alpha][1]+(n-1)\[Beta][0]),0]+
Sum[h[i](i(i-1)\[Alpha][n-i]+i \[Beta][n-i-1]+\[Gamma][n-i-2]),{i,Max[2,n-6],n-2}]+
h[1](\[Beta][n-2]+\[Gamma][n-3])+
h[0]\[Gamma][n-2]),prec]);
\[Kappa][n_]:=(\[Kappa][n]=SetPrecision[-(1/(\[Alpha][0]n(n-1)))(
If[n>2,\[Kappa][n-1]((n-1)(n-2)\[Alpha][1]+(n-1)\[Beta][0]),0]+
Sum[\[Kappa][i](i(i-1)\[Alpha][n-i]+i \[Beta][n-i-1]+\[Gamma][n-i-2]),{i,Max[2,n-6],n-2}]+
\[Kappa][1](\[Beta][n-2]+\[Gamma][n-3])+
\[Kappa][0]\[Gamma][n-2]
+Sum[\[Delta][j]h[n-2-j],{j,0,Min[n-2,4]}]
),prec]);
h[0]=h0;
h[1]= h1;
\[Kappa][0]=\[Kappa]0;
\[Kappa][1]= \[Kappa]1;
\[Alpha][0]=(a^2+(-2+r) r)^2;
\[Alpha][1]=4 (-1+r) (a^2+(-2+r) r);
\[Alpha][2]=2 (2+a^2-6 r+3 r^2);
\[Alpha][3]=4 (-1+r);
\[Alpha][4]=1;
\[Alpha][n_]/;n>4:=0;
\[Beta][0]=2 (-1+r) (a^2+(-2+r) r);
\[Beta][1]=2 (2+a^2-6 r+3 r^2);
\[Beta][2]=6 (-1+r);
\[Beta][3]=2 ;
\[Beta][n_]/;n>3:=0;
\[Gamma][-1]=0;
\[Gamma][0]=-((-2+r) r \[Lambda]0)-2 a^3 m \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+r^4 \[Omega]^2+a^2 (m^2-\[Lambda]0+2 r^2 \[Omega]^2);
\[Gamma][1]=2 (\[Lambda]0-r \[Lambda]0+2 r \[Omega] (-a m+a^2 \[Omega]+r^2 \[Omega]));
\[Gamma][2]=-\[Lambda]0+2 \[Omega] (-a m+a^2 \[Omega]+3 r^2 \[Omega]);
\[Gamma][3]=4 r \[Omega]^2;
\[Gamma][4]=\[Omega]^2;
\[Gamma][n_]/;n>4:=0;
\[Delta][0]=-(1/2) (a^2+(-2+r) r) (r^2+a^2 \[Gamma]ll);
\[Delta][1]=(3-2 r) r^2-a^2 (r-\[Gamma]ll+r \[Gamma]ll);
\[Delta][2]=-3 (-1+r) r-1/2 a^2 (1+\[Gamma]ll);
\[Delta][3]=1-2 r;
\[Delta][4]=-(1/2);
\[Delta][n_]/;n>4:=0;
Return[
Transpose@{
{
{{r},Sequence@@Table[n!h[n],{n,0,nmax}]},
FudgeZero[Abs[h[nmax]]],
(Evaluate[Sum[h[n](#-r)^n,{n,0,nmax}]]&),
Table[n!h[n],{n,0,nmax}]
},
{
{{r},Sequence@@Table[n!\[Kappa][n],{n,0,nmax}]},
FudgeZero[Abs[\[Kappa][nmax]]],
(Evaluate[Sum[\[Kappa][n](#-r)^n,{n,0,nmax}]]&),
Table[n!\[Kappa][n],{n,0,nmax}]
}
}]
]


(* ::Subsubsection::Closed:: *)
(*h\[Kappa]TSolve*)


h\[Kappa]TSolve::icaccuraccy="Warning: asymptotic expansion of order `1` at r=`2` provides insufficiently accurate initial conditions to meet PrecsionGoal, `3`."


(* ::Input::Initialization:: *)
h\[Kappa]TSolve[m_,a_,\[Omega]_,\[Lambda]_,\[Gamma]ll_,{r0_,r1_},seriesF_Symbol,
OptionsPattern[{
PrecisionGoal-> $TeukolskyPrecisionGoal,
WorkingPrecision-> $TeukolskyWorkingPrecision,
MaxSteps-> 10^5,
Order-> 30
}]
]:=
Module[
{
itable,
nmax=OptionValue[Order],
WP=Max[
OptionValue@WorkingPrecision,
OptionValue@PrecisionGoal+OptionValue@Order,
$MinPrecision
],
err=10^-(OptionValue[PrecisionGoal]+2),
rd0,
rd1,
ics
},
Block[
{
$RecursionLimit=$RecursionLimit+OptionValue[MaxSteps],
$MinPrecision=WP
},
ics=seriesF[m,a,\[Omega],\[Lambda],\[Gamma]ll,r0,OptionValue[Order]];
If[ics[[1,3]]>10^-(OptionValue[PrecisionGoal]) Abs@ics[[1,1]]||ics[[2,3]]>10^-(OptionValue[PrecisionGoal]) Abs@ics[[2,1]],
Message[h\[Kappa]TSolve::icaccuraccy,OptionValue[Order],r0,OptionValue[PrecisionGoal]]
];
itable=
Reap[
rd0=h\[Kappa]Ds[m,a,\[Omega],\[Lambda],\[Gamma]ll,r0,ics[[1,1;;2]],ics[[2,1;;2]],nmax];
Sow[First[rd0]];
h\[Kappa]TSolveStep[m,a,\[Omega],\[Lambda],\[Gamma]ll,rd0,r0,r0+ Sign[r1-r0]SetPrecision[
Min[
Max[
Min[(err/rd0[[2]])^(1/nmax)],
(r0^2-2r0+a^2)/(r0^2+a^2) Abs[r1-r0]/OptionValue[MaxSteps]
],
Quiet[Abs[\[Pi]/\[Omega]],{Power::infy}],
Abs[r1-r0]
],WP],
r1,err,nmax,WP,(r0^2-2r0+a^2)/(r0^2+a^2) Abs[r1-r0]/OptionValue[MaxSteps]]
][[2,1]];
{
Interpolation[Sort@itable[[All,1]],InterpolationOrder->2nmax],
Interpolation[Sort@itable[[All,2]],InterpolationOrder->2nmax]
}
]
]


(* ::Subsubsection::Closed:: *)
(*InterpolationStep*)


(* ::Input::Initialization:: *)
h\[Kappa]TSolveStep[m_,a_,\[Omega]_,\[Lambda]_,\[Gamma]ll_,RDprev_,rprev_,rnew_,rmax_,err_,nmax_,WP_,minstep_]:=
With[
{
RDnew=h\[Kappa]Ds[m,a,\[Omega],\[Lambda],\[Gamma]ll,rnew,
SetPrecision[{RDprev[[3,1]][rnew],Derivative[1][RDprev[[3,1]]][rnew]},WP],
SetPrecision[{RDprev[[3,2]][rnew],Derivative[1][RDprev[[3,2]]][rnew]},WP]
,nmax]
},
With[
{
ht0=RDnew[[3,1]][(rnew+rprev)/2],
ht1=RDprev[[3,1]][(rnew+rprev)/2],
\[Kappa]t0=RDnew[[3,2]][(rnew+rprev)/2],
\[Kappa]t1=RDprev[[3,2]][(rnew+rprev)/2]
},
If[
(Abs[rnew-rprev]>=minstep)&&( Abs@((ht0-ht1)/FudgeZero@Max[Abs@ht0+Abs@ht1])>err|| Abs@((\[Kappa]t0-\[Kappa]t1)/FudgeZero@Max[Abs@\[Kappa]t0+Abs@\[Kappa]t1])>err),
h\[Kappa]TSolveStep[m,a,\[Omega],\[Lambda],\[Gamma]ll,RDprev,rprev,(rprev+rnew)/2,rmax,err,nmax,WP,minstep],
Sow[First[RDnew]];
If[
Abs[rnew-rmax]>rmax 10^-WP,
h\[Kappa]TSolveStep[m,a,\[Omega],\[Lambda],\[Gamma]ll,RDnew,rnew,SetPrecision[rnew+Sign[rnew-rprev] 
Min[
Max[
Min[(err/RDnew[[2]])^(1/nmax)],
(rnew^2-2rnew+a^2)/(rnew^2+a^2)/((rprev^2-2rprev+a^2)/(rprev^2+a^2))  minstep
],
Quiet[Abs[\[Pi]/\[Omega]],{Power::infy}],
Abs[rmax-rnew],
2Abs[rnew-rprev]
],
WP],rmax,err,nmax,WP,(rnew^2-2rnew+a^2)/(rnew^2+a^2)/((rprev^2-2rprev+a^2)/(rprev^2+a^2))minstep]
]
]
]
]


(* ::Subsubsection::Closed:: *)
(*h\[Kappa]Mode*)


(* ::Input::Initialization:: *)
h\[Kappa]Mode[{m_,a_,\[Omega]_,\[Lambda]0_,\[Gamma]ll_},{inGrid_NGrid,upGrid_NGrid},opts:OptionsPattern[
{
PrecisionGoal:> $TeukolskyPrecisionGoal,
WorkingPrecision:> $TeukolskyWorkingPrecision,
MaxSteps-> 10^5,
Order-> 30,
"rinf"->Automatic,
"rhor"->Automatic,
"RadialDerivatives"->2
}]
]:=Module[{hins,\[Kappa]ins,hups,\[Kappa]ups,insol,upsol,rhor,rinf,rinmax,rupmin,out},
rhor=If[OptionValue["rhor"]===Automatic,Min@NGridPointList[inGrid],OptionValue["rhor"]];
rinf=If[OptionValue["rinf"]===Automatic,Max@NGridPointList[upGrid],OptionValue["rinf"]];
rinmax=Max@NGridPointList[inGrid];(*Print[{rhor,rinmax}];*)
rupmin=Min@NGridPointList[upGrid];(*Print[{rinf,rupmin}];*)
insol=h\[Kappa]TSolve[m,a,\[Omega], \[Lambda]0,\[Gamma]ll,{rhor,rinmax},h\[Kappa]horSeries,FilterRules[{opts},{PrecisionGoal,WorkingPrecision,MaxSteps,Order}]];
upsol=h\[Kappa]TSolve[m,a,\[Omega], \[Lambda]0,\[Gamma]ll,{rinf,rupmin},h\[Kappa]infSeries,FilterRules[{opts},{PrecisionGoal,WorkingPrecision,MaxSteps,Order}]];
hins={insol[[1]]/@inGrid,Derivative[1][insol[[1]]]/@inGrid};
\[Kappa]ins={insol[[2]]/@inGrid,Derivative[1][insol[[2]]]/@inGrid};
hups={upsol[[1]]/@upGrid,Derivative[1][upsol[[1]]]/@upGrid};
\[Kappa]ups={upsol[[2]]/@upGrid,Derivative[1][upsol[[2]]]/@upGrid};
out={
h\[Kappa]Ds[m,a,\[Omega],\[Lambda]0, \[Gamma]ll,inGrid,hins,\[Kappa]ins,OptionValue["RadialDerivatives"]][[4]],
h\[Kappa]Ds[m,a,\[Omega],\[Lambda]0, \[Gamma]ll,upGrid,hups,\[Kappa]ups,OptionValue["RadialDerivatives"]][[4]]
};
Return[out]
]


(* ::Subsection:: *)
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


(* ::Subsubsection:: *)
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
(*Spin 1 calculation*)


(* ::Subsubsection::Closed:: *)
(*Spin 1 derivatives*)


(* ::Input::Initialization:: *)
Dm1Pp1[a_,m_,\[Omega]_,r_]:=Function[{Pp1,Pp11,Pp12},Evaluate[
Pp11+(Pp1 (2+I a m-2 r-I a^2 \[Omega]-I r^2 \[Omega]))/(a^2-2 r+r^2)
]
]


(* ::Input::Initialization:: *)
Ddagm1Pm1[a_,m_,\[Omega]_,r_]:=Function[{Pm1,Pm11,Pm12},Evaluate[
Pm11+(Pm1 (2-I a m-2 r+I a^2 \[Omega]+I r^2 \[Omega]))/(a^2-2 r+r^2)
]
]


(* ::Input::Initialization:: *)
DdagPp1[a_,m_,\[Omega]_,r_]:=Function[{Pp1,Pp11,Pp12},Evaluate[
Pp11+(I Pp1 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)
]
]


(* ::Input::Initialization:: *)
DPm1[a_,m_,\[Omega]_,r_]:=Function[{Pm1,Pm11,Pm12},Evaluate[
Pm11-(I Pm1 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)
]
]


(* ::Input::Initialization:: *)
DopDPm1[a_,m_,\[Omega]_,r_]:=Function[{Pm1,Pm11,Pm12},Evaluate[
Pm12-(2 I Pm11 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)-(Pm1 (-2 I a m+a^2 m^2+2 I a m r+2 I a^2 \[Omega]-2 a^3 m \[Omega]-2 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


(* ::Input::Initialization:: *)
DopDdagPp1[a_,m_,\[Omega]_,r_]:=Function[{Pp1,Pp11,Pp12},Evaluate[
Pp12+(Pp1 (-2 I a m+a^2 m^2+2 I a m r+2 I a^2 \[Omega]-2 a^3 m \[Omega]-2 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


(* ::Input::Initialization:: *)
DdagDPm1[a_,m_,\[Omega]_,r_]:=Function[{Pm1,Pm11,Pm12},Evaluate[
Pm12+(Pm1 (2 I a m+a^2 m^2-2 I a m r-2 I a^2 \[Omega]-2 a^3 m \[Omega]+2 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


(* ::Input::Initialization:: *)
DdagDdagPp1[a_,m_,\[Omega]_,r_]:=Function[{Pp1,Pp11,Pp12},Evaluate[
Pp12+(2 I Pp11 (-a m+a^2 \[Omega]+r^2 \[Omega]))/(a^2-2 r+r^2)-(Pp1 (2 I a m+a^2 m^2-2 I a m r-2 I a^2 \[Omega]-2 a^3 m \[Omega]+2 I r^2 \[Omega]-2 a m r^2 \[Omega]+a^4 \[Omega]^2+2 a^2 r^2 \[Omega]^2+r^4 \[Omega]^2))/(a^2-2 r+r^2)^2
]
]


(* ::Subsubsection:: *)
(*Spin 1 MRProjectors*)


(* ::Input::Initialization:: *)
MRProjector[1,"l+l+"][a_,m_,\[Omega]_,lmax_,r_,{Pm1vec_,Pp1vec_,Dm1Pp1_,Ddagm1Pm1_,
DdagPp1_,DPm1_,calPs1_,DcalPs1_,DdagcalPs1_},{ Bmatm1_,Bmatp1_}]:=Module[
{
tmp,Ss1lplp,t0,t1
},
t0=-(Bmatm1+signmat[lmax,m] . Bmatp1) . \[Lambda]hat[lmax,m];
t1=a \[Omega]*(Bmatm1 . SinMixMatrix[-1,0,lmax,m]+signmat[lmax,m] . Bmatp1 . SinMixMatrix[1,0,lmax,m]); Ss1lplp=t0+t1;
Ss1lplp=t0+t1;
tmp=2/(a^2-2 r+r^2)*Dm1Pp1;
Return[tmp . signmat[lmax,m] . Ss1lplp]
]


(* ::Input::Initialization:: *)
MRProjector[1,"l-l-"][a_,m_,\[Omega]_,lmax_,r_,{Pm1vec_,Pp1vec_,Dm1Pp1_,Ddagm1Pm1_,
DdagPp1_,DPm1_,calPs1_,DcalPs1_,DdagcalPs1_},{ Bmatm1_,Bmatp1_}]:=Module[
{
tmp,Ss1lmlm,t0,t1
},
t0=-(Bmatm1+signmat[lmax,m] . Bmatp1) . \[Lambda]hat[lmax,m];
t1=a \[Omega]*(Bmatm1 . SinMixMatrix[-1,0,lmax,m]+signmat[lmax,m] . Bmatp1 . SinMixMatrix[1,0,lmax,m]); Ss1lplp=t0+t1;
Ss1lmlm=t0+t1;
tmp=2/(a^2-2 r+r^2)*Ddagm1Pm1;
Return[tmp . Ss1lmlm]
]


(* ::Input::Initialization:: *)



(* ::Input::Initialization:: *)
MRProjector[1,"m+m+"][a_,m_,\[Omega]_,lmax_,r_,{Pm1vec_,Pp1vec_,Dm1Pp1_,Ddagm1Pm1_,
DdagPp1_,DPm1_,calPs1_,DcalPs1_,DdagcalPs1_},{ Bmatm1_,Bmatp1_}]:=Module[
{
tmp,Ss1mpmp
},
Ss1mpmp=Bmatp1 . (-\[Lambda]2hat[lmax,m]+a \[Omega]*SinMixMatrix[1,2,lmax,m]);
tmp=2*(signmat[lmax,m] . DdagPp1+DPm1);
Return[tmp . signmat [lmax,m] . Ss1mpmp]
]


(* ::Input::Initialization:: *)
MRProjector[1,"m-m-"][a_,m_,\[Omega]_,lmax_,r_,{Pm1vec_,Pp1vec_,Dm1Pp1_,Ddagm1Pm1_,
DdagPp1_,DPm1_,calPs1_,DcalPs1_,DdagcalPs1_},{ Bmatm1_,Bmatp1_}]:=Module[
{
tmp,Ss1mmmm
},
Ss1mmmm=Bmatm1 . (\[Lambda]2hat[lmax,m]-a \[Omega]*SinMixMatrix[-1,-2,lmax,m]);
tmp=2*(signmat[lmax,m] . DdagPp1+DPm1);
Return[-tmp . Ss1mmmm]
]


(* ::Input::Initialization:: *)
MRProjector[1,"l+m+"][a_,m_,\[Omega]_,lmax_,r_,{Pm1vec_,Pp1vec_,Dm1Pp1_,Ddagm1Pm1_,
DdagPp1_,DPm1_,calPs1_,DcalPs1_,DdagcalPs1_},{ Bmatm1_,Bmatp1_}]:=Module[
{
sinmix,cosmix,\[CapitalSigma]mat,line1,line2,tmp1,tmp2, SS0,t1,t0
},
t0=-(Bmatm1+signmat[lmax,m] . Bmatp1) . \[Lambda]hat[lmax,m];
t1=a \[Omega]*(Bmatm1 . SinMixMatrix[-1,0,lmax,m]+signmat[lmax,m] . Bmatp1 . SinMixMatrix[1,0,lmax,m]);SS0=t0+t1;
cosmix=CosMixMatrix[1,1,lmax,m];
sinmix=SinMixMatrix[0,1,lmax,m];line1=(signmat[lmax,m] . (r*DcalPs1-2*calPs1)) . Bmatp1+I*a*(signmat[lmax,m] . DcalPs1) . Bmatp1 . cosmix;tmp1=1/(a^2-2 r+r^2)*signmat[lmax,m] . Pp1vec;tmp2=-r*\[Lambda]hat[lmax,m]-I*a*\[Lambda]hat[lmax,m] . cosmix+(a*\[Omega]*r+2*I*a)*sinmix+I*a^2*\[Omega]*Normal[sinmix . cosmix];line2=tmp1 . SS0 . tmp2;Return[line1+line2]
]


(* ::Input::Initialization:: *)
MRProjector[1,"l+m-"][a_,m_,\[Omega]_,lmax_,r_,{Pm1vec_,Pp1vec_,Dm1Pp1_,Ddagm1Pm1_,
DdagPp1_,DPm1_,calPs1_,DcalPs1_,DdagcalPs1_},{ Bmatm1_,Bmatp1_}]:=Module[
{
sinmix,cosmix,\[CapitalSigma]mat,line1,line2,tmp1,tmp2, SS0,t1,t0
},
t0=-(Bmatm1+signmat[lmax,m] . Bmatp1) . \[Lambda]hat[lmax,m];
t1=a \[Omega]*(Bmatm1 . SinMixMatrix[-1,0,lmax,m]+signmat[lmax,m] . Bmatp1 . SinMixMatrix[1,0,lmax,m]);SS0=t0+t1;
cosmix=CosMixMatrix[-1,1,lmax,m];
sinmix=SinMixMatrix[0,-1,lmax,m];
line1=-(r*DcalPs1-2*calPs1) . Bmatm1+I*a*DcalPs1 . Bmatm1 . CosMixMatrix[-1,1,lmax,m];
tmp1=1/(a^2-2 r+r^2)*signmat[lmax,m]  . Pp1vec;
tmp2=r*\[Lambda]hat[lmax,m] -I*a*\[Lambda]hat [lmax,m] . cosmix-(a*\[Omega]*r+2*I*a)*sinmix+I*a^2*\[Omega]*Normal[sinmix . cosmix];
line2=tmp1 . SS0 . tmp2;
Return[line1+line2]
]


(* ::Input::Initialization:: *)
MRProjector[1,"l-m+"][a_,m_,\[Omega]_,lmax_,r_,{Pm1vec_,Pp1vec_,Dm1Pp1_,Ddagm1Pm1_,
DdagPp1_,DPm1_,calPs1_,DcalPs1_,DdagcalPs1_},{ Bmatm1_,Bmatp1_}]:=Module[
{
sinmix,cosmix,\[CapitalSigma]mat,line1,line2,tmp1,tmp2, SS0,t1,t0
},
t0=-(Bmatm1+signmat[lmax,m] . Bmatp1) . \[Lambda]hat[lmax,m];
t1=a \[Omega]*(Bmatm1 . SinMixMatrix[-1,0,lmax,m]+signmat[lmax,m] . Bmatp1 . SinMixMatrix[1,0,lmax,m]);SS0=t0+t1;
cosmix=CosMixMatrix[1,1,lmax,m];
sinmix=SinMixMatrix[0,1,lmax,m];
tmp1=signmat[lmax,m] . (r*DdagcalPs1-2*calPs1);
tmp2=signmat [lmax,m] . DdagcalPs1;
line1 =tmp1 . Bmatp1-I*a*tmp2 . Bmatp1 . cosmix;
tmp1=1/(a^2-2 r+r^2)*Pm1vec;
tmp2=-r*\[Lambda]hat[lmax,m]+I*a*\[Lambda]hat[lmax,m] . cosmix+(a*\[Omega]*r-2*I*a)*sinmix-I*a^2*\[Omega]*Normal[sinmix . cosmix];
line2=tmp1 . SS0 . tmp2;
Return[line1+line2]
]


(* ::Input::Initialization:: *)
MRProjector[1,"l-m-"][a_,m_,\[Omega]_,lmax_,r_,{Pm1vec_,Pp1vec_,Dm1Pp1_,Ddagm1Pm1_,
DdagPp1_,DPm1_,calPs1_,DcalPs1_,DdagcalPs1_},{ Bmatm1_,Bmatp1_}]:=Module[
{
sinmix,cosmix,\[CapitalSigma]mat,line1,line2,tmp1,tmp2, SS0,t1,t0
},
t0=-(Bmatm1+signmat[lmax,m] . Bmatp1) . \[Lambda]hat[lmax,m];
t1=a \[Omega]*(Bmatm1 . SinMixMatrix[-1,0,lmax,m]+signmat[lmax,m] . Bmatp1 . SinMixMatrix[1,0,lmax,m]);SS0=t0+t1;
cosmix=CosMixMatrix[-1,1,lmax,m];
sinmix=SinMixMatrix[0,-1,lmax,m];

line1=-(r*DdagcalPs1-2*calPs1) . Bmatm1-I*a*DdagcalPs1 . Bmatm1 . cosmix;tmp1=1/(a^2-2 r+r^2)*Pm1vec;tmp2=r*\[Lambda]hat[lmax,m]+I*a*\[Lambda]hat[lmax,m] . cosmix-(a*\[Omega]*r-2*I*a)*sinmix-I*a^2*\[Omega]*Normal[sinmix . cosmix];line2=tmp1 . SS0 . tmp2;
Return[line1+line2]
]


(* ::Input::Initialization:: *)
MRProjector[1,"l+l-"][a_,m_,\[Omega]_,lmax_,r_,{Pm1vec_,Pp1vec_,Dm1Pp1_,Ddagm1Pm1_,
DdagPp1_,DPm1_,calPs1_,DcalPs1_,DdagcalPs1_},{ Bmatm1_,Bmatp1_}]:=Module[
{
cosmix,\[CapitalSigma]mat,term1,term2, term3, SS0,t1,t0
},
t0=-(Bmatm1+signmat[lmax,m] . Bmatp1) . \[Lambda]hat[lmax,m];
t1=a \[Omega]*(Bmatm1 . SinMixMatrix[-1,0,lmax,m]+signmat[lmax,m] . Bmatp1 . SinMixMatrix[1,0,lmax,m]);SS0=t0+t1;
cosmix=CosMixMatrix[0,1,lmax,m];
\[CapitalSigma]mat=r^2*idmat[lmax,m]+a^2*CosMixMatrix[0,2,lmax,m];

term1=calPs1 . SS0 . \[CapitalSigma]mat;
term2=-2*r*(Pm1vec+signmat[lmax,m] . Pp1vec) . SS0;
term3=2*a^2*calPs1 . (Bmatm1 . SinMixMatrix[-1,0,lmax,m] -signmat[lmax,m] . Bmatp1 . SinMixMatrix[1,0,lmax,m]) . cosmix;
Return[term1+term2+term3]
]


(* ::Subsubsection:: *)
(*CalculateSpin1Contributions*)


CalculateSpin1Contributions[a_,m_,\[Omega]_,lmax_,r_,{Pp1vec_,Pm1vec_},{Bmatm1_,Bmatp1_}]:=Module[{
Pvecs},
Pvecs={
Pm1vec[[All,1]],
Pp1vec[[All,1]],
Dm1Pp1[a,m,\[Omega],r]@@@Pp1vec,
Ddagm1Pm1[a,m,\[Omega],r]@@@Pm1vec,
DdagPp1[a,m,\[Omega],r]@@@Pp1vec,
DPm1[a,m,\[Omega],r]@@@Pm1vec,
DPm1[a,m,\[Omega],r]@@@Pm1vec+signmat[lmax,m] . DdagPp1[a,m,\[Omega],r]@@@Pp1vec,
DopDPm1[a,m,\[Omega],r]@@@Pm1vec+signmat[lmax,m] . DopDdagPp1[a,m,\[Omega],r]@@@Pp1vec,
DdagDPm1[a,m,\[Omega],r]@@@Pm1vec+signmat[lmax,m] . DdagDdagPp1[a,m,\[Omega],r]@@@Pp1vec
};
Return[
{
MRProjector[1,"l+l+"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm1,Bmatp1}],
MRProjector[1,"l-l-"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm1,Bmatp1}],
MRProjector[1,"m+m+"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm1,Bmatp1}],
MRProjector[1,"m-m-"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm1,Bmatp1}],

MRProjector[1,"l+m+"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm1,Bmatp1}],
MRProjector[1,"l+m-"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm1,Bmatp1}],
MRProjector[1,"l-m+"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm1,Bmatp1}],
MRProjector[1,"l-m-"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm1,Bmatp1}],

MRProjector[1,"l+l-"][a,m,\[Omega],lmax,r,Pvecs,{Bmatm1,Bmatp1}],

0 Pp1vec[[All,1]]
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


(* ::Subsubsection:: *)
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


(* ::Subsection:: *)
(*Output Generation*)


(* ::Subsubsection::Closed:: *)
(*GetComp functions*)


(*GetComp[2][rtry_]:=Module[{comp},
comp=Calcs2components[rtry];
Join[comp,{Table[0,{Length@comp[[1]]}]}]
];

GetComp[1][rtry_]:=Module[{comp},
Join[comp=Calcs1components[rtry],{Table[0,{Length@comp[[1]]}]}]
];

GetComp[0][rtry_]:=Module[{comp},
Calcs0components[rtry]
];*)


(* ::Subsubsection:: *)
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


(*ExportOutput[{sallL_,sallR_,s2L_,s2R_,s1L_,s1R_,s0L_,s0R_},directory_,iConfig_,dformat_]:=Module[{fn},
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
]*)

numStr[x_]:=Module[{s=ToString[x]},
Which[StringMatchQ[s,___~~"."],s<>"0",!StringContainsQ[s,"."],s<>".0",True,s]];

ExportOutput[{sallL_,sallR_},directory_,mm_,a0_,r0_,lmax_,rsL_,rsR_,dformat_]:=Module[
{fn,tolist,comps,datasets},

comps={"h_l+l+","h_l-l-","h_m+m+","h_m-m-","rho_h_l+m+","rhob_h_l+m-","rhob_h_l-m+","rho_h_l-m-","sigma_delta_h_l+l-","h"};

tolist[x_]:=Developer`ToPackedArray[N[x]];

fn=directory<>"data/"<>"h1_a"<>numStr[ Round[a0,0.01] ]<>"_rp"<>numStr[Round[r0,0.01] ]<>"_l"<>ToString[lmax]<>"_m"<>ToString[mm]<>".h5";
datasets=Join[
  Flatten[Table[{
    "/m_"<>ToString[mm]<>"/In/"<>comps[[i]]->tolist[Transpose[sallL[[All,i,All]]]],
    "/m_"<>ToString[mm]<>"/Up/"<>comps[[i]]->tolist[Transpose[sallR[[All,i,All]]]]
  },{i,1,10}],1],
  {"/m_"<>ToString[mm]<>"/r_in"->tolist[rsL],
   "/m_"<>ToString[mm]<>"/r_up"->tolist[rsR]}
];
Export[fn,datasets,"HDF5"];
]


(* ::Section:: *)
(*Main Function*)


Options[MetricReconstructRadiative]={
WorkingPrecision->32,
AccuracyGoal->9,
"rinf"->1000., (* rinf is the radius at which the series solutions for \[Kappa]_up should set the initial conditions for the integrator. *)
"xhor"->0.001,
"dformat"->"Real64"(* Data format for output files. *),
"Interpolator"->Automatic,
MaxSteps-> 10^5,
Order-> 30
}


MetricReconstructRadiative::gridmismatch = "Expecting first element of left and right NGrids (`2` and `3`) to be equal to the particle position r0=`1`."


MetricReconstructRadiative[OD_OrbitalData?CircularEquatorialQ, {lmax0_,mm0_},{lNGrid_NGrid,rNGrid_NGrid}, iConfig_,primarypath_,opts:OptionsPattern[]]:=Module[
{
r,ll,\[CapitalOmega]0,lmax,mm,
prec,a0,r0,\[Omega]0,a\[Omega]0,
lmins2,lmins1,lmins0,
rinf,xhor,accgoal,qres,rhor,dformat,
Bmats,eigs,Bmatm2,Bmatm1,Bmat0,Bmatp1,Bmatp2,
eigs2,eigs1,eigs0,
\[Gamma]mat,
directory,
sallL,sallR,s2L,s2R,s1L,s1R,s0L,s0R,
s0Lgrid,s1Lgrid,s2Lgrid,s0Rgrid,s1Rgrid,s2Rgrid,
MMh0gridL, MMh0gridR,\[Kappa]ngridL, MM\[Kappa]0gridL,MM\[Kappa]0gridR,
MMm1gridL,MMm1gridR,MMp1gridL,MMp1gridR,
MMm2gridL,MMm2gridR,MMp2gridL,MMp2gridR,
lgrid,rgrid,
s1\[Kappa]jumps,s0jumps,s2jumps,
\[Kappa]j1,\[Kappa]j0,hj1,hj0,Pm1j0,Pm1j1,Pm2j0,Pm2j1,Pp2j0,Pp2j1,
MMh\[Kappa]vec,MMm1vec,MMm2vec,
dhlplp,dhlmlm,dhmpmp,dhmmmm,d\[Rho]chlpmm,sourcetermsubs
},
Print["1. Preliminaries."];
(* Read the parameters from a configuration file *)
(* Read the config file *)

directory=(primarypath<>"/data"<>extractUpToSecondLastHyphen[ToString[iConfig]]<>"/data"<>extractUpToLastHyphen[ToString[iConfig]]<>"/");
If[!DirectoryQ[directory], CreateDirectory[directory, CreateIntermediateDirectories -> True]];

(* Read in the numerical parameters. All parameters should be read in this section only (no longer spread throughout the code). Unlike the above, the following sections are done numerically *)
prec=OptionValue[WorkingPrecision]; (* Number of digits to use where required. This must be read in first. *)
a0=SetPrecision[ODspin[OD],prec];
If[PossibleZeroQ[a0],a0=0]; (* In the Schwarzschild case, a0 should be the integer zero. *)
r0=SetPrecision[ODsemilatusrectum[OD],prec];
If[r0!= lNGrid[[2,1]]||r0!= rNGrid[[2,1]], Message[MetricReconstructRadiative::gridmismatch, r0,lNGrid[[2,1]],rNGrid[[2,1]]] ];
lmax=lmax0;
mm=mm0;
\[CapitalOmega]0=KerrAzimuthalFrequency[OD];
\[Omega]0=mm \[CapitalOmega]0;
a\[Omega]0=If[a0==0,0,a0*\[Omega]0];

lmins0=Abs[mm];
lmins1=Max[1,Abs[mm]];
lmins2=Max[2,Abs[mm]];

(* To handle the \[Kappa] functions, satisfying equations sourced by the radial profile of the trace, I have written some bespoke code (!) *)
rinf=SetPrecision[OptionValue["rinf"],prec]; (* rinf is the radius at which the series solutions for \[Kappa]_up should set the initial conditions for the integrator. *)
xhor=SetPrecision[OptionValue["xhor"],prec];
rhor=(1+Sqrt[1-a0^2]+xhor);
accgoal=OptionValue["AccuracyGoal"];
dformat=ToString@OptionValue["dformat"]; (* Data format for output files. *)

(* Set up the b matrices, and the eigenvalue matrices. *)
PrintTemporary["Calculating B-matrices and Eigenvalues"];
{Bmats,eigs}=Transpose@Table[
bmatEV[s,mm,a\[Omega]0,lmax,WorkingPrecision->prec,AccuracyGoal->accgoal],
{s,-2,2}
];
{Bmatm2,Bmatm1,Bmat0,Bmatp1,Bmatp2}=Bmats;
{eigs2,eigs1,eigs0}=eigs[[1;;3]];

\[Gamma]mat=Bmat0 . CosMixMatrix[0,2,lmax,mm] . (Transpose@Bmat0);


Print["2. Generate MSTModes"];
Block[{
	$MST\[Nu]Precision=accgoal,
	$MST\[Nu]WorkingPrecision=prec,
	$CFPrecision=prec-5,
	$TeukolskyWorkingPrecision=prec,
	$TeukolskyPrecisionGoal=accgoal,
	$SWSHaccuracy=accgoal,
	$SWSHEVPrecision=prec
},

EchoTiming[MMm2vec = Table[If[ll>=lmins2,MSTMode[{-2,ll,mm,a0,\[Omega]0},{lNGrid,rNGrid},"RadialDerivatives"-> 4 ,"CleanUp"->False,"Interpolator"->OptionValue["Interpolator"] ],0],{ll,0,lmax}],"Generate s=-2 MSTModes"];
EchoTiming[MMm1vec = Table[If[ll>=lmins2,MSTMode[{-1,ll,mm,a0,\[Omega]0},{lNGrid,rNGrid},"RadialDerivatives"-> 2,"CleanUp"->False,"Interpolator"->OptionValue["Interpolator"] ],0],{ll,0,lmax}],"Generate s=-1 MSTModes"];

];

EchoTiming[Monitor[
MMh\[Kappa]vec=Table[If[ll>=lmins0,h\[Kappa]Mode[{mm,a0,\[Omega]0,eigs0[[ll+1]],\[Gamma]mat[[ll+1,ll+1]]},{lNGrid,rNGrid},
"rinf"->rinf,"rhor"-> rhor,
PrecisionGoal->accgoal,
WorkingPrecision->prec,
MaxSteps-> OptionValue[MaxSteps],
Order-> OptionValue[Order]
],0],{ll,0,lmax}],
{"Generate h\[Kappa]Modes ",ll}
],"Generate h\[Kappa]Modes"];

Print["3. Jumps"];
Block[
{
	$CFPrecision=prec-5,
	$SWSHaccuracy=accgoal,
	$SWSHEVPrecision=prec
},
s2jumps=Flatten@Table[
Join[
CalcJump[{2,mm,eigs2[[ll+1]]},OD,{Pp2j0[ll],Pp2j1[ll]},
 {#[[1]][0]/Sqrt[2\[Pi]],-(#[[2]][0]/Sqrt[2\[Pi]])}&@SpinWeightedSpheroidalHarmonics[2,ll,mm,a\[Omega]0] 
]
,
CalcJump[{-2,mm,eigs2[[ll+1]]},OD,{Pm2j0[ll],Pm2j1[ll]},
 {#[[1]][0]/Sqrt[2\[Pi]],-(#[[2]][0]/Sqrt[2\[Pi]])}&@SpinWeightedSpheroidalHarmonics[-2,ll,mm,a\[Omega]0] 
]
]
,
{ll,lmins2,lmax}
];
 
s0jumps=Flatten@Table[CalcJump[{0,mm},OD,{hj0[ll],hj1[ll]},
{#[[1]][0]/Sqrt[2\[Pi]]}&@SpinWeightedSpheroidalHarmonics[0,ll,mm,a\[Omega]0]
]
,{ll,lmins0,lmax}];

sourcetermsubs=Flatten@Table[
Calcsourcesubs[OD,{dhlplp[ll],dhlmlm[ll],dhmpmp[ll],dhmmmm[ll],d\[Rho]chlpmm[ll]},
{
 SpinWeightedSphericalHarmonic[0,ll,mm][0]/Sqrt[2\[Pi]],
 SpinWeightedSphericalHarmonic[-1,ll,mm][0]/Sqrt[2\[Pi]],
	 SpinWeightedSphericalHarmonic[2,ll,mm][0]/Sqrt[2\[Pi]],
	 SpinWeightedSphericalHarmonic[-2,ll,mm][0]/Sqrt[2\[Pi]]
}],
{ll,lmins0,lmax}
];
];

EchoTiming[
s1\[Kappa]jumps=CalcJumps\[Kappa]s1[
{\[Kappa]j0,\[Kappa]j1,Pm1j0,Pm1j1},
OD, 
{lmax,mm},
{
{Pm2j0,Pm2j1,s2jumps},
{hj0,hj1,s0jumps},
{dhlplp,dhlmlm,dhmpmp,dhmmmm,d\[Rho]chlpmm,sourcetermsubs}
},
{Bmatm2,Bmatm1,Bmat0,Bmatp1,Bmatp2},
{eigs0,eigs1,eigs2},
{\[Gamma]mat}],
"\[Kappa] and s1 jumps"
];

Print["4. Prepare grid and mode functions."];

EchoTiming[
Monitor[{MMm2gridL,MMm2gridR}=
 Transpose@Table[
 If[ll>=lmins2,
 With[
 {MM=MMm2vec[[ll+1]]},
 ConstructSolution[{MSTRhor[MM]/MSTAhor[MM],MSTRinf[MM]/MST\[Beta]inf[MM]}, {Pm2j0[ll],Pm2j1[ll]}/.s2jumps]
 ]
 ,
 {{0,0,0,0,0}lNGrid,{0,0,0,0,0}rNGrid}
 ] 
 ,{ll,0,lmax}]
 ,
 {"s=-2 MSTModes ",ll}
 ],
 "s=-2 MSTModes"
 ];
 
EchoTiming[MMp2gridL=Table[If[ll>=lmins2,Pp2fromPm2[MMm2gridL[[ll+1]],lNGrid,{ll,mm,a0,\[Omega]0,eigs2[[ll+1]]}],{0,0,0,0,0}lNGrid],{ll,0,lmax}],"MMp2gridL"];
EchoTiming[MMp2gridR=Table[If[ll>=lmins2,Pp2fromPm2[MMm2gridR[[ll+1]],rNGrid,{ll,mm,a0,\[Omega]0,eigs2[[ll+1]]}],{0,0,0,0,0}rNGrid],{ll,0,lmax}],"MMp2gridR"];

EchoTiming[
Monitor[{MMm1gridL,MMm1gridR}=
 Transpose@Table[
 If[ll>=lmins2,
 With[
 {MM=MMm1vec[[ll+1]]},
 ConstructSolution[{MSTRhor[MM]/MSTAhor[MM],MSTRinf[MM]/MST\[Beta]inf[MM]}, {Pm1j0[ll],Pm1j1[ll]}/.s1\[Kappa]jumps]
 ]
 ,
 {{0,0,0}lNGrid,{0,0,0}rNGrid}
 ] 
 ,{ll,0,lmax}]
 ,
 {"s=-1 MSTModes ",ll}
 ],
 "s=-1 MSTModes"
 ]; 
   
EchoTiming[MMp1gridL=Table[If[ll>=lmins1,Pp1fromPm1[MMm1gridL[[ll+1]],lNGrid,{ll,mm,a0,\[Omega]0,eigs1[[ll+1]]}],{0,0,0}lNGrid],{ll,0,lmax}],"MMp1gridL"];
EchoTiming[MMp1gridR=Table[If[ll>=lmins1,Pp1fromPm1[MMm1gridR[[ll+1]],rNGrid,{ll,mm,a0,\[Omega]0,eigs1[[ll+1]]}],{0,0,0}rNGrid],{ll,0,lmax}],"MMp1gridR"];

EchoTiming[
Monitor[{MMh0gridL,MMh0gridR,MM\[Kappa]0gridL,MM\[Kappa]0gridR}=
 Transpose@Table[If[ll>=lmins0,
 ConstructSolutionh\[Kappa][MMh\[Kappa]vec[[ll+1]],{hj0[ll],hj1[ll]}/.s0jumps,{\[Kappa]j0[ll],\[Kappa]j1[ll]}/.s1\[Kappa]jumps],
 {{0,0,0}lNGrid,{0,0,0}rNGrid,{0,0,0}lNGrid,{0,0,0}rNGrid}
 ],
 {ll,0,lmax}],
 {"s=0 MSTModes ",ll}],
 "s=0 MSTModes"
 ]; 

Print["A. Fill 1D grid."];
EchoTiming[s2Rgrid=CalculateSpin2Contributions[a0,mm,\[Omega]0,lmax,rNGrid,{MMp2gridR,MMm2gridR},{Bmatm2,Bmatp2}],"s = 2, UP ..."];
s2R=Transpose[s2Rgrid/. NGrid[_,data_]:> data,{2,3,1}];
EchoTiming[s1Rgrid=CalculateSpin1Contributions[a0,mm,\[Omega]0,lmax,rNGrid,{MMp1gridR,MMm1gridR},{Bmatm1,Bmatp1}],"s = 1, UP ..."];
s1R=Transpose[s1Rgrid/. NGrid[_,data_]:> data,{2,3,1}];
EchoTiming[s0Rgrid=CalculateSpin0Contributions[a0,mm,\[Omega]0,lmax,rNGrid,{MMh0gridR,MM\[Kappa]0gridR},{Bmat0,\[Gamma]mat,eigs0}],"s = 0, UP ..."];
s0R=Transpose[s0Rgrid/. NGrid[_,data_]:> data,{2,3,1}];

EchoTiming[s2Lgrid=CalculateSpin2Contributions[a0,mm,\[Omega]0,lmax,lNGrid,{MMp2gridL,MMm2gridL},{Bmatm2,Bmatp2}],"s = 2, IN ..."];
s2L=Transpose[s2Lgrid/. NGrid[_,data_]:> data,{2,3,1}];
EchoTiming[s1Lgrid=CalculateSpin1Contributions[a0,mm,\[Omega]0,lmax,lNGrid,{MMp1gridL,MMm1gridL},{Bmatm1,Bmatp1}],"s = 1, IN ..."];
s1L=Transpose[s1Lgrid/. NGrid[_,data_]:> data,{2,3,1}];
EchoTiming[s0Lgrid=CalculateSpin0Contributions[a0,mm,\[Omega]0,lmax,lNGrid,{MMh0gridL,MM\[Kappa]0gridL},{Bmat0,\[Gamma]mat,eigs0}],"s = 0, IN ..."];
s0L=Transpose[s0Lgrid/. NGrid[_,data_]:> data,{2,3,1}];

sallR=s2R+s1R+s0R;
sallL=s2L+s1L+s0L;

lgrid= lNGrid/. NGrid[_,data_]:> data;
rgrid= rNGrid/. NGrid[_,data_]:> data;

ClearMST[];
ClearSWSH[];

(* Save in a binary format. *)
ExportOutput[{sallL,sallR},directory,mm,a0,r0,lmax,lgrid,rgrid,dformat];
]


End[];
EndPackage[];
