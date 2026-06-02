(* ::Package:: *)

(* ::Title:: *)
(*SWSHdecomp Package*)


(* ::Input::Initialization:: *)
BeginPackage["SWSHdecomp`", 
{
"SpinWeightedSpheroidalHarmonicsFT`"
}]


(* ::Section:: *)
(*Preamble*)


ClearSWSHdecomp::usage="";
SWSHb::usage="";
$SWSHblmax::usage="";


(* ::Section:: *)
(*Private*)


(* ::Input::Initialization:: *)
Begin["`Private`"];


CFm2[s_,m_,lmax_]:=CFm2[s,m,lmax]=SparseArray[
Flatten[Table[
{j+1-Abs[m],l+1-Abs[m]}->1/3 KroneckerDelta[j,l]HeavisideTheta[j-Abs[s]+1/2]+If[j+l>=2,2/3 Sqrt[(2j+1)/(2l+1)]ClebschGordan[{j,m},{2,0},{l,m}]ClebschGordan[{j,-s},{2,0},{l,-s}],0],
{l,Max[Abs[s],Abs[m]],lmax},
{j,Max[Abs[s],Abs[m],l-2],Min[l+2,lmax]}
],1],{lmax+1-Abs[m],lmax+1-Abs[m]}]


CFm1[s_,m_,lmax_]:=CFm1[s,m,lmax]=SparseArray[
Flatten[Table[
{j+1-Abs[m],l+1-Abs[m]}-> If[j+l>=1,Sqrt[(2j+1)/(2l+1)]ClebschGordan[{j,m},{1,0},{l,m}]ClebschGordan[{j,-s},{1,0},{l,-s}],0],
{l,Max[Abs[s],Abs[m]],lmax},
{j,Max[Abs[s],Abs[m],l-1],Min[l+1,lmax]}
],1],{lmax+1-Abs[m],lmax+1-Abs[m]}]


CFm0[s_,m_,lmax_]:=CFm0[s,m,lmax]=DiagonalMatrix[Table[If[l>=Max[Abs[s],Abs[m]], l(l+1),0],{l,Abs[m],lmax}]]


SWSHb[s_,m_,c_,
OptionsPattern[{
WorkingPrecision-> $SWSHdaccuracy,
"lmax"-> $SWSHblmax
}]]:=With[
{$lmax=OptionValue["lmax"]},
 (Sign[Plus@@#]#&)/@#[[2,Join[Range[-1,Abs[m]-Abs[s],-1],Ordering[#[[1,1;;Min[0,Abs[m]-Abs[s]]-1]],All,Greater]]]]&@
Eigensystem[
SetPrecision[
c^2 CFm2[s,m,$lmax]-2c s CFm1[s,m,$lmax]-CFm0[s,m,$lmax],OptionValue[WorkingPrecision]
]
]
]


ClearSWSHdecomp[]:= (
DownValues[CFm0]=DownValues[CFm0][[-1;;-1]];
DownValues[CFm1]=DownValues[CFm1][[-1;;-1]];
DownValues[CFm2]=DownValues[CFm2][[-1;;-1]];
)


(* ::Section:: *)
(*End*)


$SWSHblmax=30;
$SWSHdaccuracy:=$SWSHaccuracy;


(* ::Input::Initialization:: *)
End[];
EndPackage[];
