(* ::Package:: *)

(* ::Title:: *)
(*PTools Package*)


BeginPackage["PTools`"]


(* ::Section:: *)
(*Preamble*)


PTable::usage=""


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


SetAttributes[PTable,HoldAll]


PTable[a_,it:{_Symbol,__}..]:=
WaitAll@With[{vars={it}[[All,1]]},
Table[
ParallelSubmit[vars,a]
,
it
]
]


(* ::Section:: *)
(*End*)


End[];
EndPackage[];
