(* ::Package:: *)

BeginPackage["NContinuedFraction`"];


NContinuedFraction::usage="NContinuedFraction[a[i],b[i],{i,\[Epsilon],\!\(\*SubscriptBox[\(i\), \(min\)]\)}] numericallly evaluates the continued fraction \!\(\*FractionBox[\(a[0]\), \(\(b[0]\)\(+\)\)]\)\!\(\*FractionBox[\(a[1]\), \(\(b[1]\)\(+\)\)]\)... until it converges with relative error \[Epsilon]. \!\(\*SubscriptBox[\(i\), \(min\)]\) is the minimum number of terms to evaluate."


Begin["`Private`"];


NContinuedFraction[a_,b_,{j_Symbol,err_,imin_:3}]:=
Module[{mn=$ModuleNumber,f,g,h,\[CapitalDelta]h,d,i=imin,out},
f[i_]:= a/.j->i;
g[i_]:= b/.j->i;
d[i_]:=(d[i]=1/(g[i]+f[i]d[i-1]));
h[i_]:=(h[i]=(h[i-1]g[i]+h[i-2]d[i-1]f[i])d[i]);

d[0]=1/g[0];h[-1]=0;h[0]=f[0]/g[0];
\[CapitalDelta]h[i_]:=(\[CapitalDelta]h[i]=h[i]-h[i-1]);
While[Abs[\[CapitalDelta]h[i]]>err Abs[h[i]],i++];
out=h[i];
Remove@@{"NContinuedFraction`Private`d$"<>ToString[mn],"NContinuedFraction`Private`f$"<>ToString[mn],"NContinuedFraction`Private`g$"<>ToString[mn],"NContinuedFraction`Private`h$"<>ToString[mn]};
out
];


End[];
EndPackage[];
