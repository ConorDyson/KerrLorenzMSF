(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["SpinWeightedSpheroidalHarmonicsFT`",{"NContinuedFraction`"}];


(* ::Input::Initialization:: *)
$SWSHaccuracy::usage="$SWSHaccuracy set the number of significant digits in which the SpinWeightedSpheroidalHarmonicsFT should do its work. Default value is $MachinePrecision.";
$SWSHEVPrecision::usage="$SWSHaccuracy set the number of significant digits to find the eigenvalues of the Spin Weighted Spheroidal Harmonics. Default value is 2$MachinePrecision.";
SpinWeightedSphericalHarmonic::usage="SpinWeightedSphericalHarmonic[s,l,m,z] gives the (polar part of) the spin-weighted spherical harmonic\!\(\*SubscriptBox[\(\\\ \\\ \), \(s\)]\)\!\(\*SubscriptBox[\(S\), \(lm\)]\)[z], where z=Cos[\[Theta]].";
SpinWeightedSpheroidalHarmonicEV::usage="SpinWeightedSpheriodalHarmonicEV[s,l,m,c] gives the seperation constant of the spin-weighted spheriodal harmonic\!\(\*SubscriptBox[\(\\\ \\\ \), \(s\)]\)\!\(\*SubsuperscriptBox[\(S\), \(lm\), \(c\)]\)[z].";
SpinWeightedSpheroidalHarmonic::usage="SpinWeightedSpheriodalalHarmonic[s,l,m,c][z] gives the (polar part of) the spin-weighted spheriodal harmonic\!\(\*SubscriptBox[\(\\\ \\\ \), \(s\)]\)\!\(\*SubsuperscriptBox[\(S\), \(lm\), \(c\)]\)[z], where z=Cos[\[Theta]]. Normalization is such that \!\(\*SubsuperscriptBox[\(\[Integral]\), \(-1\), \(1\)]\)\!\(\*SubscriptBox[\(\\\ \\\ \), \(s\)]\)\!\(\*SubscriptBox[\(S\), \(lm\)]\)[z\!\(\*SuperscriptBox[\(]\), \(2\)]\)dz == 1. The sign is defined such that\!\(\*SubscriptBox[\(\\\ \), \(s\)]\)\!\(\*SubsuperscriptBox[\(S\), \(lm\), \(c\)]\)[z] is positive in a neighbourhood of z == 1.";
SpinWeightedSpheroidalHarmonics::usage="SpinWeightedSpheriodalalHarmonic[s,l,m,c][z] returns a list of the (polar part of) the spin-weighted spheriodal harmonic\!\(\*SubscriptBox[\(\\\ \\\ \), \(s\)]\)\!\(\*SubsuperscriptBox[\(S\), \(lm\), \(c\)]\)[z], where z=Cos[\[Theta]], and its first two derivatives. Normalization is such that \!\(\*SubsuperscriptBox[\(\[Integral]\), \(-1\), \(1\)]\)\!\(\*SubscriptBox[\(\\\ \\\ \), \(s\)]\)\!\(\*SubscriptBox[\(S\), \(lm\)]\)[z\!\(\*SuperscriptBox[\(]\), \(2\)]\)dz == 1. The sign is defined such that\!\(\*SubscriptBox[\(\\\ \), \(s\)]\)\!\(\*SubsuperscriptBox[\(S\), \(lm\), \(c\)]\)[z] is positive in a neighbourhood of z == 1.";
SpinWeightedSpheroidalHarmonicEVcExpansion::usage="SpinWeightedSpheroidalHarmonicEVcExpansion[s,l,m,M] gives the small c analytic expansion of SpinWeightedSpheriodalHarmonicEV[s,l,m,c] to the Mth order.";
SetSWSHAccuracy::usage="SetSWSHAccuracy[precision] sets the value of $SWSHaccuracy to precision. All stored Spin-Weighted Harmonics and their eigenvalues are reset.";
SWSHevGuess::usage="";
ClearSWSH::usage="";



(* ::Section:: *)
(*Private*)


(* ::Input::Initialization:: *)
Begin["`Private`"];


ModuleReturn[x_,mn_]:=(Remove@@Names["SpinWeightedSpheroidalHarmonicsFT`Private`@$"<>ToString[mn]];x);


SetSWSHAccuracy[precision_]:=($SWSHEVPrecision=2precision;$SWSHaccuracy=precision);



(* ::Subsection:: *)
(*Misc*)


$FRc[{a_-> b_}]:=b
$FRc[_]:=$Failed


(* ::Subsubsection::Closed:: *)
(*SpinWeightedSpheroidalHarmonicEVcExpansion*)


(* ::Input::Initialization:: *)
SpinWeightedSpheroidalHarmonicEVcExpansion[s_,l_,m_,-1]:=(0&);SpinWeightedSpheroidalHarmonicEVcExpansion[s_,l_,m_,M_]:=
SpinWeightedSpheroidalHarmonicEVcExpansion[s,l,m,M]=Module[{mn=$ModuleNumber,km= Abs[m-s],kp= Abs[m+s],pl,A,i,j=Floor[M/2],c,\[Alpha]\[Gamma],\[Beta],R,L},
pl=l-(kp+km)/2;
\[Alpha]\[Gamma][p_]:=-((4 c^2 p (km+p) (kp+p) (km+kp+p) (km+kp+2 p-2 s) (km+kp+2 (p+s)))/((-1+km+kp+2 p) (km+kp+2 p)^2 (1+km+kp+2 p)));
\[Beta][p_]:= A+c^2-(p+(kp+km)/2)(p+(kp+km)/2+1)+(2c s (kp-km)(kp+km))/((2p+kp+km)(2p+kp+km+2));
R[imax_]:=ContinuedFractionK[-\[Alpha]\[Gamma][pl+i+1],\[Beta][pl+i+1],{i,0,imax}];
L[imax_]:=ContinuedFractionK[-\[Alpha]\[Gamma][pl-i],\[Beta][pl-i-1],{i,0,imax}];
Subscript[A, M]=Simplify[Expand[Subscript[A, M]/.Solve[Normal[Series[\[Beta][pl]+R[j]+L[j]==0/.A-> SpinWeightedSpheroidalHarmonicEVcExpansion[s,l,m,M-1][c]+Subscript[A, M] c^M,{c,0,M}]],Subscript[A, M]][[1]]]/.Abs[a_]^x_?EvenQ:> a^x];
ModuleReturn[Function[{c},SpinWeightedSpheroidalHarmonicEVcExpansion[s,l,m,M-1][c]+Subscript[A, M] c^M],mn]
];


(* ::Subsubsection::Closed:: *)
(*guess*)


(* ::Input::Initialization:: *)
SWSHevGuess[s_,l_,m_,c_]:=Module[{mn=$ModuleNumber,f1,f2,A1,A3,l1,l2,q,g,c0},
(*small c*)
f1:=l (1+l)-(2 c m s^2)/(l+l^2)+(c^2 (-l^3 (1+l)^3 (-1+2 l (1+l)-2 m^2)+4 (l+l^2)^2 (l+l^2-3 m^2) s^2+2 (-3 l^2 (1+l)^2+(3+5 l (1+l)) m^2) s^4))/(l^3 (1+l)^3 (-3+4 l (1+l)))+(4 c^3 m s^2 (l^4 (1+l)^4 (-1+3 l (1+l)-5 m^2)-2 l^2 (1+l)^2 (5 l^2 (1+l)^2-(6+7 l (1+l)) m^2) s^2+(-6 m^2+l (1+l) (-19 m^2+l (1+l) (6+7 l (1+l)-9 m^2))) s^4))/((-1+l) l^5 (1+l)^5 (2+l) (-1+2 l) (3+2 l))+(2 c^4 (6 l^7+5 l^8-123 l^9-426 l^10-396 l^11+493 l^12+1352 l^13+722 l^14-824 l^15-1327 l^16-489 l^17+332 l^18+430 l^19+197 l^20+44 l^21+4 l^22+60 l^7 m^2+338 l^8 m^2+678 l^9 m^2+446 l^10 m^2-68 l^11 m^2+656 l^12 m^2+2272 l^13 m^2+2048 l^14 m^2-488 l^15 m^2-2458 l^16 m^2-2214 l^17 m^2-1006 l^18 m^2-240 l^19 m^2-24 l^20 m^2-66 l^7 m^4-469 l^8 m^4-1422 l^9 m^4-2326 l^10 m^4-1982 l^11 m^4-196 l^12 m^4+1666 l^13 m^4+2258 l^14 m^4+1624 l^15 m^4+713 l^16 m^4+180 l^17 m^4+20 l^18 m^4-4 l^6 (1+l)^6 (-2+l+l^2) (l (1+l) (15+l (1+l) (-13+12 l (1+l)))-2 (45+l (1+l) (-87+100 l (1+l))) m^2+9 (-5+28 l (1+l)) m^4) s^2+2 l^4 (1+l)^4 ((-1+l) l^2 (1+l)^2 (2+l) (90+l (1+l) (-141+220 l (1+l)))-2 (-90+l (1+l) (159+l (1+l) (704+3 l (1+l) (-933+476 l (1+l))))) m^2+(3330+l (1+l) (-5883+l (1+l) (-2693+3276 l (1+l)))) m^4) s^4-4 l^2 (1+l)^2 (9 (-1+l) l^4 (1+l)^4 (2+l) (-5+28 l (1+l))-2 l^2 (1+l)^2 (1260+l (1+l) (-2226+l (1+l) (-991+1332 l (1+l)))) m^2+(2430+l (1+l) (567+l (1+l) (-10089+13 l (1+l) (183+220 l (1+l))))) m^4) s^6+(9 (-1+l) l^4 (1+l)^4 (2+l) (-45+l (1+l) (57+68 l (1+l)))-2 l^2 (1+l)^2 (2430+l (1+l) (567+l (1+l) (-10089+13 l (1+l) (183+220 l (1+l))))) m^2+(4050+l (1+l) (9045+l (1+l) (-17955+l (1+l) (-22272+l (1+l) (16489+5876 l (1+l)))))) m^4) s^8))/((-1+l) l^7 (1+l)^7 (2+l) (-3+2 l) (5+2 l) (-3+4 l (1+l))^3)-(4 c^5 m s^2 (138 l^8+151 l^9-3396 l^10-14215 l^11-19842 l^12+3095 l^13+36771 l^14+18086 l^15-60210 l^16-111567 l^17-77978 l^18-8775 l^19+26342 l^20+21945 l^21+8435 l^22+1680 l^23+140 l^24+1920 l^8 m^2+11770 l^9 m^2+25320 l^10 m^2+16140 l^11 m^2-2250 l^12 m^2+73380 l^13 m^2+275520 l^14 m^2+411840 l^15 m^2+297300 l^16 m^2+44010 l^17 m^2-106640 l^18 m^2-99900 l^19 m^2-42330 l^20 m^2-9240 l^21 m^2-840 l^22 m^2-438 l^8 m^4-6341 l^9 m^4-36969 l^10 m^4-118380 l^11 m^4-231708 l^12 m^4-282630 l^13 m^4-195846 l^14 m^4-33156 l^15 m^4+71310 l^16 m^4+73827 l^17 m^4+34423 l^18 m^4+8280 l^19 m^4+828 l^20 m^4-20 l^6 (1+l)^6 (l^2 (1+l)^2 (-192+l (1+l) (359+3 l (1+l) (-129+28 l (1+l))))-6 (-45+l (1+l) (87+l (1+l) (-5+3 l (1+l) (-97+28 l (1+l))))) m^2+(1350+l (1+l) (-2610+l (1+l) (-963+484 l (1+l)))) m^4) s^2+2 l^4 (1+l)^4 (l^2 (1+l)^2 (-1350+l (1+l) (2610+l (1+l) (-588+l (1+l) (-11567+3348 l (1+l)))))-10 (-135+l (1+l) (-9+l (1+l) (4947+l (1+l) (-8696+l (1+l) (-3145+1804 l (1+l)))))) m^2+(45630+l (1+l) (3042+l (1+l) (-153201+l (1+l) (2737+16484 l (1+l))))) m^4) s^4-20 l^2 (1+l)^2 (l^4 (1+l)^4 (1350+l (1+l) (-2610+l (1+l) (-963+484 l (1+l))))-2 l^2 (1+l)^2 (3105+l (1+l) (207+l (1+l) (-10506+l (1+l) (343+1196 l (1+l))))) m^2+3 (1620+l (1+l) (3348+l (1+l) (-5823+l (1+l) (-7581+l (1+l) (1839+700 l (1+l)))))) m^4) s^6+(3 l^4 (1+l)^4 (4860+l (1+l) (324+l (1+l) (-16047+l (1+l) (-231+1508 l (1+l)))))-30 l^2 (1+l)^2 (1620+l (1+l) (3348+l (1+l) (-5823+l (1+l) (-7581+l (1+l) (1839+700 l (1+l)))))) m^2+(34020+l (1+l) (138348+l (1+l) (-7317+l (1+l) (-388737+l (1+l) (-156432+l (1+l) (96187+17884 l (1+l))))))) m^4) s^8))/(If[l==2,1,(-2+l)] (-1+l) l^9 (1+l)^9 (2+l) (3+l) (-3+2 l) (5+2 l) (-3+4 l (1+l))^3);
(*large c*)
l1:=Abs[m+s]+s;
l2:=Abs[m-s]-s;
q:=Which[l>= Max[l1,l2],l+1-(1-(-1)^(l-l1))/2,
l<l1,2l-l1+1,
l<l2,2l-l2+1
];
A1:=-(1/8)(q^3-m^2 q+q-s^2 (q+m));
A3:=1/512 (1/64 ((q-m-1+2s)(q-m-1-2s)(q-m-2s-3)(q-m+2s-3)(q+m-3)^2 (q+m-1)^2-
(q-m-2s+1)(q-m+2s+1)(q-m-2s+3)(q-m+2s+3)(q+m+1)^2 (q+m+3)^2)+
2((q-m+2s-1)(q-m-2s-1)(q-1)^2 (q+m-1)^2-
(q+m+1)^2 (q+1)^2 (q-m-2s+1)(q-m+2s+1))-
2A1((q-m+2s-1)(q-m-2s-1)(q+m-1)^2+(q+m+1)^2 (q-m-2s+1)(q-m+2s+1)));
f2:=-c^2+2q c -1/2 (q^2-m^2+2s+1)+1/c A1-1/(64c^2) (5q^4-(6m^2-10)q^2+m^4-2m^2-4s^2 (q^2-m^2-1)+1)+1/c^3 A3;
(*Combine*)
c0=1/2 (((l^4 (1+l)^4 (-3+4 l (1+l)))/Abs[2 (l+l^2)^3 (-1+l+l^2+m^2)+4 (l+l^2)^2 (l+l^2-3 m^2) s^2+2 (-3 l^2 (1+l)^2+(3+5 l (1+l)) m^2) s^4])^(1/2)+Abs[A3]^(1/3));
g:= (Tanh[(c-c0)/c0]+1)/2 HeavisideTheta[c-3/4 c0];
ModuleReturn[f1 (1-g)+f2 g,mn]
];


(* ::Subsubsection:: *)
(*clear*)


ClearSWSH[]:=(
DownValues[SpinWeightedSpheroidalHarmonicEVcExpansion]=DownValues[SpinWeightedSpheroidalHarmonicEVcExpansion][[-2;;-1]];
DownValues[SpinWeightedSpheroidalHarmonicEV]=
DownValues[SpinWeightedSpheroidalHarmonicEV][[-5;;-1]];
DownValues[SpinWeightedSpheroidalHarmonic]=
DownValues[SpinWeightedSpheroidalHarmonic][[-5;;-1]];
)


(* ::Subsection:: *)
(*Spherical Harmoninic*)


SpinWeightedSphericalHarmonic[s_,l_,m_,diff_:0]/;Abs[m]>l:=(0&)
SpinWeightedSphericalHarmonic[s_,l_,m_,diff_:0]/;Abs[s]>l:=(0&)


SpinWeightedSphericalHarmonic[0,l_,m_]/;Abs[m]>l:=(0&)
SpinWeightedSphericalHarmonic[0,l_,m_]:=Function[{z},Evaluate[Sqrt[(2l+1)/2] Sqrt[(l-m)!/(l+m)!]LegendreP[l,m,z]]]


(* ::Input::Initialization:: *)
SpinWeightedSphericalHarmonic[s_,l_,m_,diff_:0]:=
SpinWeightedSphericalHarmonic[s,l,m,diff]=Function[{z},Evaluate[(-1)^m Simplify[D[Sqrt[((l+m)!(l-m)!(2l+1))/((l+s)!(l-s)!2)] ((1-z)/2)^l \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(r = 0\), \(l - s\)]\((Binomial[l - s, r] Binomial[l + s, r + s - m] 
\*SuperscriptBox[\((\(-1\))\), \(l - r - s\)] 
\*SuperscriptBox[\((
\*FractionBox[\(1 + z\), \(1 - z\)])\), \(r + 
\*FractionBox[\(s - m\), \(2\)]\)])\)\),{z,diff}],Assumptions->{-1<z<1}]]];


(* ::Subsection:: *)
(*Private Functions*)


(* ::Subsubsection:: *)
(*\[Alpha]\[Gamma]*)


\[Alpha]\[Gamma][s_,l_,m_,c_]:=
With[{km= Abs[m-s],kp= Abs[m+s]},
Function[{p},
Evaluate[-((4 c^2 p (km+p) (kp+p) (km+kp+p) (km+kp+2 p-2 s) (km+kp+2 (p+s)))/((-1+km+kp+2 p) (km+kp+2 p)^2 (1+km+kp+2 p)))]
]
 ]


(* ::Subsubsection:: *)
(*\[Alpha]*)


\[Alpha][s_,l_,m_,c_]:=
With[{km= Abs[m-s],kp= Abs[m+s]},
Function[{p},
Evaluate[(4c(p+kp+1)(p+km+1)(p+(kp+km)/2+1-s))/((2p+kp+km+2)(2p+kp+km+3))]
]
 ]


(* ::Subsubsection::Closed:: *)
(*\[Beta]*)


\[Beta][s_,l_,m_,c_]:=
With[{km= Abs[m-s],kp= Abs[m+s]},
Function[{p,x},
Evaluate[x+c^2-(p+(kp+km)/2)(p+(kp+km)/2+1)+(2c s (kp-km)(kp+km))/((2p+kp+km)(2p+kp+km+2))]
]
 ]


\[Beta][s_,l_,m_,c_,A_]:=
With[{km= Abs[m-s],kp= Abs[m+s]},
Function[{p},
Evaluate[A+c^2-(p+(kp+km)/2)(p+(kp+km)/2+1)+(2c s (kp-km)(kp+km))/((2p+kp+km)(2p+kp+km+2))]
]
 ]


(* ::Subsubsection::Closed:: *)
(*\[Gamma]*)


\[Gamma][s_,l_,m_,c_]:=
With[{km= Abs[m-s],kp= Abs[m+s]},
Function[{p},
Evaluate[-((4c p(p+kp+km)(p+(kp+km)/2+s))/((2p+kp+km-1)(2p+kp+km)))]
]
 ]


(* ::Subsubsection::Closed:: *)
(*\[Delta]*)


\[Delta][s_,l_,m_,c_]:=
With[{km= Abs[m-s],kp= Abs[m+s]},
Function[{p},
Evaluate[-((p+(kp+km)/2-s)/(p+(kp+km)/2+s))]
]
 ]


(* ::Subsubsection:: *)
(*Rt*)


Rt[s_,l_,m_,c_,err_]:=
With[{\[Alpha]\[Gamma]= \[Alpha]\[Gamma][s,l,m,c],\[Beta]= \[Beta][s,l,m,c]},
Function[{p,x},
NContinuedFraction[-\[Alpha]\[Gamma][p+CFi+1],\[Beta][p +CFi+1,x],{CFi,err,2l+1}]
]
 ]


(* ::Subsubsection::Closed:: *)
(*Lt*)


Lt[s_,l_,m_,c_,err_]:=
With[{\[Alpha]\[Gamma]= \[Alpha]\[Gamma][s,l,m,c],\[Beta]= \[Beta][s,l,m,c]},
Function[{p,x},
ContinuedFractionK[-\[Alpha]\[Gamma][p-CFi],\[Beta][p-CFi-1,x],{CFi,0,p-1}]
]
 ]


(* ::Subsubsection::Closed:: *)
(*R*)


R[s_,l_,m_,c_,err_,A_]:=
With[{
\[Alpha]\[Gamma]= \[Alpha]\[Gamma][s,l,m,c],
\[Alpha]= \[Alpha][s,l,m,c],
\[Beta]= \[Beta][s,l,m,c,A]},
Function[{p},
1/\[Alpha][p-1] NContinuedFraction[-\[Alpha]\[Gamma][p+i],\[Beta][p+i],{i,err,2l+1}]
]
 ]


(* ::Subsubsection::Closed:: *)
(*L*)


L[s_,l_,m_,c_,err_,A_]:=
With[{
\[Alpha]\[Gamma]= \[Alpha]\[Gamma][s,l,m,c],
\[Gamma]= \[Gamma][s,l,m,c],
\[Beta]= \[Beta][s,l,m,c,A]},
Function[{p},
1/\[Gamma][p+1] ContinuedFractionK[-\[Alpha]\[Gamma][p-i+1],\[Beta][p-i],{i,0,p}]
]
 ]


(* ::Subsubsection::Closed:: *)
(*R2*)


R2[s_,l_,m_,c_,err_,A_]:=
With[{
\[Alpha]\[Gamma]= \[Alpha]\[Gamma][s,l,m,c],
\[Alpha]= \[Alpha][s,l,m,c],
\[Delta]= \[Delta][s,l,m,c],
\[Beta]= \[Beta][s,l,m,c,A]},
Function[{p},
\[Delta][p]/\[Alpha][p-1] NContinuedFraction[-\[Alpha]\[Gamma][p+i],\[Beta][p+i],{i,err,2l+1}]
]
 ]


(* ::Subsubsection::Closed:: *)
(*L2*)


L2[s_,l_,m_,c_,err_,A_]:=
With[{
\[Alpha]\[Gamma]= \[Alpha]\[Gamma][s,l,m,c],
\[Gamma]= \[Gamma][s,l,m,c],
\[Delta]= \[Delta][s,l,m,c],
\[Beta]= \[Beta][s,l,m,c,A]},
Function[{p},
1/(\[Delta][p+1]\[Gamma][p+1]) ContinuedFractionK[-\[Alpha]\[Gamma][p-i+1],\[Beta][p-i],{i,0,p}]
]
 ]


(* ::Subsection:: *)
(*Eigen value*)


Clear[SpinWeightedSpheroidalHarmonicEV];
SpinWeightedSpheroidalHarmonicEV[s_,l_,m_,_?PossibleZeroQ]:=l(l+1);

SpinWeightedSpheroidalHarmonicEV[s_Integer,l_Integer,m_Integer,c_,Ag_:Automatic]/;Abs[s]>l:=0;

SpinWeightedSpheroidalHarmonicEV[s_Integer,l_Integer,m_Integer,c_?Negative,Ag_:Automatic]:=SpinWeightedSpheroidalHarmonicEV[s,l,-m,Abs[c],Ag];

SpinWeightedSpheroidalHarmonicEV[s_Integer,l_Integer,m_Integer,c_?NumericQ]/;Abs[m]<= l:=SpinWeightedSpheroidalHarmonicEV[s,l,m,c,SWSHevGuess[s,l,m,SetPrecision[c,Max[$MinPrecision,$SWSHEVPrecision]]]];


SpinWeightedSpheroidalHarmonicEV[s_Integer,l_Integer,m_Integer,c_?NumericQ,Ag_]/;Abs[m]<= l:= 
With[
{
pl=Max[l-(Abs[m-s]+Abs[m+s])/2,0],
\[Beta]=\[Beta][s,l,m,c],
R=Rt[s,l,m,c,10^-($SWSHEVPrecision)],
L=Lt[s,l,m,c,10^-($SWSHEVPrecision)],
Ag2=If[NumericQ[Ag],
Ag,
SWSHevGuess[s,l,m,N[c,$SWSHEVPrecision]]
]
},
Block[{$MinPrecision=Max[$MinPrecision,$SWSHEVPrecision+1]},
Quiet[
Module[{A,tf},
tf[A_?NumericQ]:=\[Beta][pl,A]+R[pl,A]+L[pl,A];
A=$FRc[FindRoot[tf[A]==0,{A,Ag2},
WorkingPrecision->$SWSHEVPrecision+1,
PrecisionGoal->$SWSHEVPrecision,
Compiled->False
]];
Remove[tf];
A
],
{SetPrecision::precsm}
]
]
];


(* ::Subsection::Closed:: *)
(*SpheroidalHarmonic*)


(* ::Input::Initialization:: *)
Clear[SpinWeightedSpheroidalHarmonic];
SpinWeightedSpheroidalHarmonic[s_Integer,l_Integer,m_Integer,_?PossibleZeroQ,diff_:0]:=SpinWeightedSphericalHarmonic[s,l,m,diff];
SpinWeightedSpheroidalHarmonic[s_Integer,l_Integer,m_Integer,c_,diff_:0]/;c<0:=Function[{z},(-1)^(l+s+diff) SpinWeightedSpheroidalHarmonic[s,l,-m,-c,diff][-z]];


SpinWeightedSpheroidalHarmonic[s_,l_,m_,c_,diff_:0]/;Abs[m]>l:=(0&)
SpinWeightedSpheroidalHarmonic[s_,l_,m_,c_,diff_:0]/;Abs[s]>l:=(0&)


(* ::Input::Initialization:: *)
SpinWeightedSpheroidalHarmonic[s_Integer,l_Integer,m_Integer,c_,diff_:0]:= With[
{
pl=Max[l-(Abs[m-s]+Abs[m+s])/2,0],
km= Abs[m-s],
kp= Abs[m+s],
A=SpinWeightedSpheroidalHarmonicEV[s,l,m,c],
err=10^-($SWSHaccuracy)
},
With[
{
L=L[s,l,m,c,10^-($SWSHaccuracy),A],
R=R[s,l,m,c,10^-($SWSHaccuracy),A],
(*\[Alpha]=\[Alpha][s,l,m,c],
\[Beta]=\[Alpha][s,l,m,c,A],
\[Gamma]=\[Alpha][s,l,m,c],
\[Delta]=\[Alpha][s,l,m,c],
\[Alpha]\[Gamma]=\[Alpha]\[Gamma][s,l,m,c],*)
L2=L2[s,l,m,c,10^-($SWSHaccuracy),A],
R2=R2[s,l,m,c,10^-($SWSHaccuracy),A]
},
Module[{pmax,a,a2,M,sum1,sum2,sum12},
a[pl]=1;
a[p2_]/;p2<pl:=(a[p2]=L[p2]a[p2+1]);
a[p2_]/;p2>pl:=(a[p2]=R[p2]a[p2-1]);
pmax=pl;
a2[pl]=1;
a2[p2_]/;p2<pl:=(a2[p2]=L2[p2]a2[p2+1]);
a2[p2_]/;p2>pl:=(a2[p2]=R2[p2]a2[p2-1]);
sum1=Sum[a[p]Binomial[p+kp,p],{p,0,pmax}];
sum2=Sum[a2[p]Binomial[p+kp,p],{p,0,pmax}];
sum12=Sum[a[p]a2[p] (2Gamma[p+kp+1]Gamma[p+km+1])/((2p+kp+km+1)Gamma[p+1]Gamma[p+kp+km+1]),{p,0,pmax}];
While[
Abs[(a[pmax]Binomial[pmax+kp,pmax])/sum1]>err,
pmax++;
sum1+=a[pmax]Binomial[pmax+kp,pmax];
sum2+=a2[pmax]Binomial[pmax+kp,pmax];
sum12+=a[pmax]a2[pmax] (2(pmax+kp)!(pmax+km)!)/((2pmax+kp+km+1)pmax!(pmax+kp+km)!)
];
M=If[m+s>=0,(-1)^kp,1](-1)^s \[Sqrt](sum2/(Exp[2c]sum1) 1/sum12);
SpinWeightedSpheroidalHarmonic[s,l,m,c]=Function[{z},Evaluate[M Exp[c z]((1-z)/2 )^(kp/2)((1+z)/2 )^(km/2)Sum[a[p]JacobiP[p,kp,km,z],{p,0,pmax}]]];
SpinWeightedSpheroidalHarmonic[s,l,m,c,0]=SpinWeightedSpheroidalHarmonic[s,l,m,c];
SpinWeightedSpheroidalHarmonic[s,l,m,c,1]=Function[{z},Evaluate[M/2 E^(c z) ((1-z)/2 )^(kp/2) ((1+z)/2 )^(km/2) Sum[a[p] ((1+km+kp+p) JacobiP[-1+p,1+kp,1+km,z]-(kp/(1-z)-km/(1+z)-2 c ) JacobiP[p,kp,km,z]),{p,0,pmax}]]];
SpinWeightedSpheroidalHarmonic[s,l,m,c,2]=Function[{z},Evaluate[M/4 E^(c z) ((1-z)/2 )^(kp/2) ((1+z)/2 )^(km/2) Sum[a[p] ((1+km+kp+p) (2+km+kp+p) JacobiP[-2+p,2+kp,2+km,z]+((2 km (1+km+kp+p) )/(1+z)-(2 kp(1+km+kp+p))/(1-z)+ 4 c (1+km+kp+p) )JacobiP[-1+p,1+kp,1+km,z]+(((kp-2) kp )/(1-z)^2-(4c kp )/(1-z)+4 c^2+((km-2) km )/(1+z)^2-(2 (kp+2 c (-1+z))km)/(1-z^2))JacobiP[p,kp,km,z] ),{p,0,pmax}]]];
Switch[diff,
1,SpinWeightedSpheroidalHarmonic[s,l,m,c,1],
2,SpinWeightedSpheroidalHarmonic[s,l,m,c,2],
_,SpinWeightedSpheroidalHarmonic[s,l,m,c]
]
]
]
];



Clear[SpinWeightedSpheroidalHarmonics];
SpinWeightedSpheroidalHarmonics[s_Integer,l_Integer,m_Integer,_?PossibleZeroQ]:={
SpinWeightedSphericalHarmonic[s,l,m,0],
SpinWeightedSphericalHarmonic[s,l,m,1],
SpinWeightedSphericalHarmonic[s,l,m,2]
};

SpinWeightedSpheroidalHarmonics[s_Integer,l_Integer,m_Integer,c_]/;c<0:=
With[{S=SpinWeightedSpheroidalHarmonics[s,l,-m,-c]},
{Function[z,(-1)^(s+l) S[[1]][-z]],
Function[z,(-1)^(s+l+1) S[[2]][-z]],
Function[z,(-1)^(s+l) S[[3]][-z]]
}
]



SpinWeightedSpheroidalHarmonics[s_Integer,l_Integer,m_Integer,c_]/;Abs[s]>l:={0&,0&,0&}
SpinWeightedSpheroidalHarmonics[s_Integer,l_Integer,m_Integer,c_]/;Abs[m]>l:={0&,0&,0&}


SpinWeightedSpheroidalHarmonics[s_Integer,l_Integer,m_Integer,c_]:=  With[
{
pl=Max[l-(Abs[m-s]+Abs[m+s])/2,0],
km= Abs[m-s],
kp= Abs[m+s],
A=SpinWeightedSpheroidalHarmonicEV[s,l,m,c],
err=10^-($SWSHaccuracy)
},
With[
{
L=L[s,l,m,c,10^-($SWSHaccuracy),A],
R=R[s,l,m,c,10^-($SWSHaccuracy),A],
L2=L2[s,l,m,c,10^-($SWSHaccuracy),A],
R2=R2[s,l,m,c,10^-($SWSHaccuracy),A]
},
Module[{pmax,a,a2,M,sum1,sum2,sum12},
a[pl]=1;
a[p2_]/;p2<pl:=(a[p2]=L[p2]a[p2+1]);
a[p2_]/;p2>pl:=(a[p2]=R[p2]a[p2-1]);
pmax=pl;
a2[pl]=1;
a2[p2_]/;p2<pl:=(a2[p2]=L2[p2]a2[p2+1]);
a2[p2_]/;p2>pl:=(a2[p2]=R2[p2]a2[p2-1]);
sum1=Sum[a[p]Binomial[p+kp,p],{p,0,pmax}];
sum2=Sum[a2[p]Binomial[p+kp,p],{p,0,pmax}];
sum12=Sum[a[p]a2[p] (2Gamma[p+kp+1]Gamma[p+km+1])/((2p+kp+km+1)Gamma[p+1]Gamma[p+kp+km+1]),{p,0,pmax}];
While[
Abs[(a[pmax]Binomial[pmax+kp,pmax])/sum1]>err,
pmax++;
sum1+=a[pmax]Binomial[pmax+kp,pmax];
sum2+=a2[pmax]Binomial[pmax+kp,pmax];
sum12+=a[pmax]a2[pmax] (2(pmax+kp)!(pmax+km)!)/((2pmax+kp+km+1)pmax!(pmax+kp+km)!)
];
M=If[m+s>=0,(-1)^kp,1](-1)^s \[Sqrt](sum2/(Exp[2c]sum1) 1/sum12);
{
Function[{z},Evaluate[M Exp[c z]((1-z)/2)^(kp/2) ((1+z)/2)^(km/2) Sum[Expand[a[p]JacP[p,kp,km,z]],{p,0,pmax}]]
],
Function[{z},Evaluate[M/8 E^(c z) ((1-z)/2)^(kp/2-1) ((1+z)/2)^(km/2-1) Sum[Expand[a[p]((1-z^2)(1+km+kp+p) JacP[-1+p,1+kp,1+km,z]-((1+z)kp-(1-z)km-2 c (1-z^2)) JacP[p,kp,km,z])],{p,0,pmax}]]
],Function[{z},Evaluate[M/64 E^(c z) ((1-z)/2)^(kp/2-2) ((1+z)/2)^(km/2-2) Sum[Expand[a[p] ((1-z^2)^2 (1+km+kp+p) (2+km+kp+p) JacP[-2+p,2+kp,2+km,z]+(2 km (1+km+kp+p) (1-z^2)(1-z)-2 kp(1+km+kp+p)(1-z^2)(1+z)+ 4 c (1+km+kp+p)(1-z^2)^2 )JacP[-1+p,1+kp,1+km,z]+((kp-2) kp (1+z)^2-4c kp (1+z)(1-z^2)+4 c^2 (1-z^2)^2+(km-2) km (1-z)^2-(1-z^2)2 (kp+2 c (-1+z))km)JacP[p,kp,km,z] )],{p,0,pmax}]]
]
}
]
]
];


(* ::Subsection:: *)
(*JacP*)


(*JacP[n_Integer,a_Integer,b_Integer,z_]:= JacP[n,a,b,z]=Expand@JacobiP[n,a,b,z]*)


JacP[n_?Negative,a_Integer,b_Integer,z_]:=0


JacP[n:0|1,a_Integer,b_Integer,z_]:= JacP[n,a,b,z]=Expand@JacobiP[n,a,b,z]


JacP[n_?(#>1&),a_Integer,b_Integer,z_]:= JacP[n,a,b,z]=Expand[
1/(2n(n+a+b)(2n+a+b-2)) (
(2n+a+b-1)(2n+a+b)(2n+a+b-2)z JacP[n-1,a,b,z]
+(2n+a+b-1)(a^2-b^2)JacP[n-1,a,b,z]-
2(n+a-1)(n+b-1)(2n+a+b)JacP[n-2,a,b,z]
)
]


(* ::Subsection:: *)
(*End*)


SetSWSHAccuracy[$MachinePrecision]


(* ::Input::Initialization:: *)
End[];


(* ::Input::Initialization:: *)
EndPackage[];
