(* ::Package:: *)

BeginPackage["MetricReconstructRadiative`"];
Needs["SpinWeightedSpheroidalHarmonics`"]
Needs["Teukolsky`"]



MetricReconstructRadiative::usage = "Outputs reconstructed metric for given config file
************************************************************************
This file was generated automatically by the Mathematica front end.  
It contains Initialization cells from a Notebook file, which         
typically will have the same name as this file except ending in      
.nb instead of .m.                                               
                                                                    
This file is intended to be loaded into the Mathematica kernel using 
the package loading commands Get or Needs.  Doing so is equivalent   
to using the Evaluate Initialization Cells menu command in the front 
end.                                                                 
                                                                     
(************************************************************************"





Begin["`Private`"];
Off[ClebschGordan::phy];
Off[ClebschGordan::tri];
Off[Infinity::indet];
Off[SpinWeightedSphericalHarmonicY::params];
Off[Power::infy];
Off[FrontEndObject::notavail];


(*NotebookPut[Notebook[{}]]*)


MetricReconstructRadiative[iConfig_,primarypath_]:=Module[{r,\[Lambda],l,\[CapitalOmega]0},


(*NotebookPut[Notebook[{}]] // UsingFrontEnd;*)
Print["1. Preliminaries."];
\[Rho]=r+I*a*Cos[\[Theta]];
\[Rho]c=r-I*a*Cos[\[Theta]];
\[CapitalSigma]=\[Rho]*\[Rho]c;
\[CapitalDelta]subs={\[CapitalDelta][r]->r^2-2*M*r+a^2,\[CapitalDelta]'[r]->2*(r-M),\[CapitalDelta]''[r]->2,\[CapitalDelta]'''[r]->0};
Ksubs={K[r]->\[Omega]*(r^2+a^2)-a*m,K'[r]-> 2*\[Omega]*r,K''[r]->2*\[Omega],K'''[r]->0};
\[CapitalDelta]Ksubs=Join[\[CapitalDelta]subs,Ksubs];
\[CapitalDelta]Ksimp={\[CapitalDelta]'[r]->2*(r-M),\[CapitalDelta]''[r]->2,\[CapitalDelta]'''[r]->0,K'[r]-> 2*\[Omega]*r,K''[r]->2*\[Omega],K'''[r]->0};
QQ=m/Sin[\[Theta]]-a*\[Omega]*Sin[\[Theta]];
a0subs={a->0};
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
tmp=gup . {-EE,0,0,LL}/.{\[Theta]->\[Pi]/2}//Simplify;
utsubs={ut->tmp[[1]]};

(* Teuksolky equations for spin 2 *)
teukPm=\[CapitalDelta][r]*Ddag[-1,Dop[Pm2[r]]]-6*I*\[Omega]*r*Pm2[r]-\[CapitalLambda]*Pm2[r];
teukPp=\[CapitalDelta][r]*Dop[-1,Ddag[Pp2[r]]]+6*I*\[Omega]*r*Pp2[r]-\[CapitalLambda]*Pp2[r];
teukSm=Lop[-1,Ldag[2,Sm2[\[Theta]]]]+6*a*\[Omega]*Cos[\[Theta]]*Sm2[\[Theta]]+\[CapitalLambda]*Sm2[\[Theta]];
teukSp=Ldag[-1,Lop[2,Sp2[\[Theta]]]]-6*a*\[Omega]*Cos[\[Theta]]*Sp2[\[Theta]]+\[CapitalLambda]*Sp2[\[Theta]];
Pm2subs2=Solve[teukPm==0,Pm2''[r]]//First//Simplify;
Pm2subs3=D[Pm2subs2,r]/.Pm2subs2/.\[CapitalDelta]Ksimp//Simplify;
Pm2subs4=D[Pm2subs3,r]/.Pm2subs2/.\[CapitalDelta]Ksimp//Simplify;
Pm2subs=Join[Pm2subs2,Pm2subs3,Pm2subs4];
Pp2subs2=Solve[teukPp==0,Pp2''[r]]//First//Simplify;
Pp2subs3=D[Pp2subs2,r]/.Pp2subs2/.\[CapitalDelta]Ksimp//Simplify;
Pp2subs4=D[Pp2subs3,r]/.Pp2subs2/.\[CapitalDelta]Ksimp//Simplify;
Pp2subs=Join[Pp2subs2,Pp2subs3,Pp2subs4];
Sm2subs2=Solve[teukSm==0,Sm2''[\[Theta]]]//First//Simplify;
Sm2subs3=D[Sm2subs2,\[Theta]]/.Sm2subs2//Simplify;
Sm2subs4=D[Sm2subs3,\[Theta]]/.Sm2subs2//Simplify;
Sm2subs=Join[Sm2subs2,Sm2subs3,Sm2subs4];
Sp2subs2=Solve[teukSp==0,Sp2''[\[Theta]]]//First//Simplify;
Sp2subs3=D[Sp2subs2,\[Theta]]/.Sp2subs2//Simplify;
Sp2subs4=D[Sp2subs3,\[Theta]]/.Sp2subs2//Simplify;
Sp2subs=Join[Sp2subs2,Sp2subs3,Sp2subs4];
teuks2subs=Join[Pm2subs,Pp2subs,Sm2subs,Sp2subs];

(* Teukolsky-Starobinskii identities for spin 2 *)
tsPp=(\[CapitalDelta][r]^2*Dop[Dop[Dop[Dop[Pm2[r]]]]]-CC*Pp2[r])/.\[CapitalDelta]Ksimp;
tsPm=(\[CapitalDelta][r]^2*Ddag[Ddag[Ddag[Ddag[Pp2[r]]]]]-Cstar*Pm2[r])/.\[CapitalDelta]Ksimp;
tsSp=Ldag[-1,Ldag[0,Ldag[1,Ldag[2,Sm2[\[Theta]]]]]]-AA*Sp2[\[Theta]];
tsSm=Lop[-1,Lop[0,Lop[1,Lop[2,Sp2[\[Theta]]]]]]-AA*Sm2[\[Theta]];
t0=Solve[tsPp==0,Pp2[r]]/.Pm2subs/.\[CapitalDelta]Ksimp//First//Simplify;
Pp2repl=Join[t0,D[t0,r]/.Pm2subs/.\[CapitalDelta]Ksimp//Simplify];
t0=Solve[tsPm==0,Pm2[r]]/.Pp2subs/.\[CapitalDelta]Ksimp//First//Simplify;
Pm2repl=Join[t0,D[t0,r]/.Pp2subs/.\[CapitalDelta]Ksimp//Simplify];
t0=Solve[tsSp==0,Sp2[\[Theta]]]/.Sm2subs//First//Simplify;
Sp2repl=Join[t0,D[t0,\[Theta]]/.Sm2subs//Simplify];
t0=Solve[tsSm==0,Sm2[\[Theta]]]/.Sp2subs//First//Simplify;
Sm2repl=Join[t0,D[t0,\[Theta]]/.Sp2subs//Simplify];


(* Teukolsky-Starobinskii constants for spin 2 *)
Crepl0={CC->AA+sgn*12*I*M*\[Omega],Cstar->AA-sgn*12*I*M*\[Omega]};
Crepl1={Cstar->(AA^2+144*M^2*\[Omega]^2)/CC,CC->(AA^2+144*M^2*\[Omega]^2)/Cstar};
Asqrepl={AA^2->\[CapitalLambda]^2*(\[CapitalLambda]+2)^2-8*\[Omega]^2*\[CapitalLambda]*(\[Alpha]sq*(5*\[CapitalLambda]+6)-12*a^2)+144*\[Omega]^4*\[Alpha]sq^2}/.{\[Alpha]sq->a^2-a*m/\[Omega]};






(* Test *)
t0=Dop[Dop[Dop[Dop[Pm2[r]]]]]/.Pm2subs/.Pm2repl/.\[CapitalDelta]Ksimp//Simplify;
t1=FullSimplify[t0/.Crepl1/.Asqrepl/.\[CapitalDelta]subs/.Ksubs];
Aschwrepl={AA->\[CapitalLambda]*(\[CapitalLambda]+2)};
sgnsq={sgn^2->1};
(* Teukolsky equations for spin 1. *)
teukPm1=\[CapitalDelta][r]*Ddag[Dop[Pm1[r]]]-2*I*\[Omega]*r*Pm1[r]-\[Lambda]*Pm1[r];
teukPp1=\[CapitalDelta][r]*Dop[Ddag[Pp1[r]]]+2*I*\[Omega]*r*Pp1[r]-\[Lambda]*Pp1[r];
teukSm1=Lop[Ldag[1,Sm1[\[Theta]]]]+2*a*\[Omega]*Cos[\[Theta]]*Sm1[\[Theta]]+\[Lambda]*Sm1[\[Theta]];
teukSp1=Ldag[Lop[1,Sp1[\[Theta]]]]-2*a*\[Omega]*Cos[\[Theta]]*Sp1[\[Theta]]+\[Lambda]*Sp1[\[Theta]];
Pm1subs2=Solve[teukPm1==0,Pm1''[r]]//First//Simplify;
Pm1subs3=D[Pm1subs2,r]/.Pm1subs2/.\[CapitalDelta]Ksimp//Simplify;
Pm1subs4=D[Pm1subs3,r]/.Pm1subs2/.\[CapitalDelta]Ksimp//Simplify;
Pm1subs=Join[Pm1subs2,Pm1subs3,Pm1subs4];
Pp1subs2=Solve[teukPp1==0,Pp1''[r]]//First//Simplify;
Pp1subs3=D[Pp1subs2,r]/.Pp1subs2/.\[CapitalDelta]Ksimp//Simplify;
Pp1subs4=D[Pp1subs3,r]/.Pp1subs2/.\[CapitalDelta]Ksimp//Simplify;
Pp1subs=Join[Pp1subs2,Pp1subs3,Pp1subs4];
Sm1subs2=Solve[teukSm1==0,Sm1''[\[Theta]]]//First//Simplify;
Sm1subs3=D[Sm1subs2,\[Theta]]/.Sm1subs2//Simplify;
Sm1subs4=D[Sm1subs3,\[Theta]]/.Sm1subs2//Simplify;
Sm1subs=Join[Sm1subs2,Sm1subs3,Sm1subs4];
Sp1subs2=Solve[teukSp1==0,Sp1''[\[Theta]]]//First//Simplify;
Sp1subs3=D[Sp1subs2,\[Theta]]/.Sp1subs2//Simplify;
Sp1subs4=D[Sp1subs3,\[Theta]]/.Sp1subs2//Simplify;
Sp1subs=Join[Sp1subs2,Sp1subs3,Sp1subs4];
teuks1subs=Join[Pm1subs,Pp1subs,Sm1subs,Sp1subs];
(* Teukolsky-Starobinskii identities for spin 1. *)
tsPp1=\[CapitalDelta][r]*Dop[Dop[Pm1[r]]]-BB*Pp1[r];
tsPm1=\[CapitalDelta][r]*Ddag[Ddag[Pp1[r]]]-BB*Pm1[r];
tsSp1=Ldag[0,Ldag[1,Sm1[\[Theta]]]]-BB*Sp1[\[Theta]];
tsSm1=Lop[0,Lop[1,Sp1[\[Theta]]]]-BB*Sm1[\[Theta]];
t0=Solve[tsPp1==0,Pp1[r]]/.Pm1subs/.\[CapitalDelta]Ksimp//First//Simplify;
Pp1repl=Join[t0,D[t0,r]/.Pm1subs/.\[CapitalDelta]Ksimp//Simplify];
t0=Solve[tsPm1==0,Pm1[r]]/.Pp1subs/.\[CapitalDelta]Ksimp//First//Simplify;
Pm1repl=Join[t0,D[t0,r]/.Pp1subs/.\[CapitalDelta]Ksimp//Simplify];
t0=Solve[tsSp1==0,Sp1[\[Theta]]]/.Sm1subs//First//Simplify;
Sp1repl=Join[t0,D[t0,\[Theta]]/.Sm1subs//Simplify];
t0=Solve[tsSm1==0,Sm1[\[Theta]]]/.Sp1subs//First//Simplify;
Sm1repl=Join[t0,D[t0,\[Theta]]/.Sp1subs//Simplify];
(* Teukolsky equations for spin 0. *)
teukP0=Dop[\[CapitalDelta][r]*Ddag[P0[r]]]-2*I*\[Omega]*r*P0[r]-\[Lambda]0*P0[r];
teukS0=Ldag[1,Lop[S0[\[Theta]]]]+2*a*\[Omega]*Cos[\[Theta]]*S0[\[Theta]]+\[Lambda]0*S0[\[Theta]];
P0subs2=Solve[teukP0==0,P0''[r]]/.\[CapitalDelta]Ksimp//First//Simplify;
P0subs3=D[P0subs2,r]/.P0subs2/.\[CapitalDelta]Ksimp//Simplify;
P0subs4=D[P0subs3,r]/.P0subs2/.\[CapitalDelta]Ksimp//Simplify;
P0subs=Join[P0subs2,P0subs3,P0subs4];
S0subs2=Solve[teukS0==0,S0''[\[Theta]]]//First//Simplify;
S0subs3=D[S0subs2,\[Theta]]/.S0subs2//Simplify;
S0subs4=D[S0subs3,\[Theta]]/.S0subs2//Simplify;
S0subs=Join[S0subs2,S0subs3,S0subs4];
teuks0subs=Join[P0subs,S0subs];
(* Teukolsky equations for trace, and for kappa. *)
teukh0=Dop[\[CapitalDelta][r]*Ddag[h0[r]]]-2*I*\[Omega]*r*h0[r]-\[Lambda]0*h0[r];
h0subs2=Solve[teukh0==0,h0''[r]]/.\[CapitalDelta]Ksimp//First//Simplify;
h0subs3=D[h0subs2,r]/.h0subs2/.\[CapitalDelta]Ksimp//Simplify;
h0subs4=D[h0subs3,r]/.h0subs2/.\[CapitalDelta]Ksimp//Simplify;
h0subs=Join[h0subs2,h0subs3,h0subs4];
teuk\[Kappa]0=Dop[\[CapitalDelta][r]*Ddag[\[Kappa]0[r]]]-2*I*\[Omega]*r*\[Kappa]0[r]-\[Lambda]0*\[Kappa]0[r]-1/2*(r^2+a^2*Cos[\[Theta]]^2)*h0[r];  (* This needs further thought! *)
\[Kappa]0subs2=Solve[teuk\[Kappa]0==0,\[Kappa]0''[r]]/.\[CapitalDelta]Ksimp//First//Simplify;
\[Kappa]0subs3=D[\[Kappa]0subs2,r]/.\[Kappa]0subs2/.\[CapitalDelta]Ksimp//Simplify;
\[Kappa]0subs4=D[\[Kappa]0subs3,r]/.\[Kappa]0subs2/.\[CapitalDelta]Ksimp//Simplify;
\[Kappa]0subs=Join[\[Kappa]0subs2,\[Kappa]0subs3,\[Kappa]0subs4]//Simplify;
Print["2. Projection + Jumps."];
(* Read the parameters from a configuration file *)
(* Read the config file *)



extractUpToLastHyphen[str_String] := StringJoin[Riffle[Most[StringSplit[str, "-"]], "-"]];
extractUpToSecondLastHyphen[str_String] := StringJoin[Riffle[Most[StringSplit[extractUpToLastHyphen[str], "-"]], "-"]];
directory=(primarypath<>"/data"<>extractUpToSecondLastHyphen[ToString[iConfig]]<>"/data"<>extractUpToLastHyphen[ToString[iConfig]]<>"/");
If[!DirectoryQ[directory], CreateDirectory[directory, CreateIntermediateDirectories -> True]];

filepath=directory<>"config/";
filen=filepath<>"config"<>ToString[iConfig]<>".txt";
s=Import[filen,"String"];
s1=StringSplit[s,"\n"];
s2=Map[StringSplit[#,"\t"]&,s1];
configparams=Map[{#[[1]],ToExpression[#[[2]]]}&,s2];
GetParam[key_]:=Module[{ls,val},
ls=Select[configparams,#[[1]]==key&];
If[Length[ls]>0,
val=ls[[1]][[2]];
If[ NumericQ[val]&&(!(IntegerQ[val])),val=SetPrecision[val,prec]];
,
val=Null
];
val
];
GetParam[key_,default_]:=If[GetParam[key]===Null,default,GetParam[key]];
(* Read in the numerical parameters. All parameters should be read in this section only (no longer spread throughout the code). Unlike the above, the following sections are done numerically *)
prec=GetParam["prec",32]; (* Number of digits to use where required. This must be read in first. *)
iround=1000;
a0=SetPrecision[Round[GetParam["a"]*iround]/iround,prec];
If[Abs[a0]<10^(-10),a0=0]; (* In the Schwarzschild case, a0 should be the integer zero. *)
r0=SetPrecision[Round[GetParam["r0"]*iround]/iround,prec];
mm=GetParam["m"];
lmax=GetParam["lmax"];
lplot=GetParam["lplot"];
nterms=GetParam["nterms"];(* Number of terms in the spherical expansion. *)
{a0,r0,mm,lmax,lplot,nterms};
(* IMPORTANT: Be cautious of using Simplify after introducing the numerical values, as this can lead to huge loss of precision. All Simplify statements must come before replacing symbols numerical values. *)
\[CapitalOmega]0=1/(Sqrt[r0^3]+a0);
\[Omega]0=mm*\[CapitalOmega]0;
a\[Omega]0=If[a0==0,0,a0*\[Omega]0];
\[CapitalDelta]Kr0subs=Map[#[[1]]->(#[[2]]/.{a->a0,r->r0,M->1,m->mm,\[Omega]->\[Omega]0})&,Join[\[CapitalDelta]subs,Ksubs]];
awsubs=Join[SetPrecision[{a->a0,\[Omega]->\[Omega]0},prec],{m->mm,M->1}];
lmins0=Abs[mm];
lmins1=Max[1,Abs[mm]];
lmins2=Max[2,Abs[mm]];
lmins={lmins2,lmins1,lmins0,lmins1,lmins2};
rprmvals={rp->1+Sqrt[1-a^2],rm->1-Sqrt[1-a^2]}/.awsubs;
(* To handle the \[Kappa] functions, satisfying equations sourced by the radial profile of the trace, I have written some bespoke code (!) *)
\[Kappa]ord=GetParam["kapord",6]; (* \[Kappa]ord is the maximum order of series expansion of \[Kappa] in spheroidal harmonics. *)
rinf=GetParam["rinf",1000.]; (* rinf is the radius at which the series solutions for \[Kappa]_up should set the initial conditions for the integrator. *)
rmax=GetParam["rmax",100.];   (* rmax is the maximum value at which the interpolating function can be used, whereas rinf is the starting value for the numerical integration, i.e. the radius at which the initial condition for the UP solutions of kappa is set, using the series expansion. *)
xhor=GetParam["xhor",0.001];
rhor=(rp+xhor)/.rprmvals;
rmin =rhor + 0.1;  (* <--- may need to think again about this choice. *)
inford=GetParam["inford",5]; (* The order of the expansion at infinity for the UP function. *)
horord=GetParam["horord",5]; (* The order of the expansion at the horizon for the IN function. *)
accgoal=GetParam["accgoal",9];
rupmin=r0/.awsubs;
rupmax=rmax;
rinfmin=(l+1/2)^2/\[Omega]/.awsubs/.{l->lplot} ; (* Should I insist on a minimum value? *)
rgrid=GetParam["rgrid",1]; (* Use a linearly-spaced grid in the variable: 0 = rstar , 1 = r. *)
rstarmin=GetParam["rstmin",SetPrecision[rmin,prec]]; (* If rgrid=1 then these will be interpreted as rmin and rmax, instead of rstarmin and rstarmax. *)
rstarmax=GetParam["rstmax",SetPrecision[rmax,prec]];
nres=GetParam["n",4];  (* Resolution in the r* direction:  dr* = M / n  (or dr = M / n). *)
qres=GetParam["angres",8]; (* Resolution in the \[Theta] direction: number of points = nres * qres. *)
dformat=ToString@GetParam["dformat","Real64"]; (* Data format for output files. *)
{\[Kappa]ord,rinf,rmax,xhor,inford,horord,accgoal,rstarmin,rstarmax,nres,qres};
(* To remove quantities that are zero to within effective machine precision *)
CleanValue[X_?NumericQ]:=Module[{\[Epsilon]=10^(-20)},If[N[Abs[X]]<\[Epsilon],0,X]];
CleanMatrix[M_]:=Module[{i,j},
M2=M;
For[i=1,i<=Dimensions[M][[1]],i++,
For[j=1,j<=Dimensions[M][[2]],j++,
M2[[i,j]]=CleanValue[M2[[i,j]]];
];
];
M2
];
(* Set up the b matrices, and the eigenvalue matrices. *)
Bmats={};
For[ss=-2,ss<=2,ss++,
barr={};
For[ll=0,ll<lmins[[ss+3]],ll++,
AppendTo[barr,Table[0,{kk,0,lmax}]];
];
For[ll=lmins[[ss+3]],ll<=lmax,ll++,
arr=Table[0,{kk,0,lmax}];
(* NEW version, for BHPToolkit versions >= 1.0 *)
If[a\[Omega]0==0,
arr[[ll+1]]=1;
,
tbl=SpinWeightedSpheroidalHarmonicS[ss,ll,mm,a\[Omega]0,Method->{"SphericalExpansion","NumTerms"->nterms}]["ExpansionCoefficients"];
For[lprime=lmins[[ss+3]],lprime<=lmax,lprime++,
If[MemberQ[Keys[tbl],lprime],arr[[lprime+1]]=(lprime/.tbl)];
];
];
(* * * * *)

(* OLD version, for BHPToolkit versions < 1.0 *)
(*
If[a\[Omega]0==0, (* Handle the Schwarzschild a=0 case separately. *)
tbl={1};
,
tbl=SpinWeightedSpheroidalHarmonicS[ss,ll,mm,a\[Omega]0,Method->{"SphericalExpansion","NumTerms"->nterms}][[5,2]];
];
ind=Ordering[tbl,-1];
For[kk=1,kk<=Length[tbl],kk++,
arr[[ll+kk+1-ind]]=tbl[[kk]];
];
*)
(* * * * *)
AppendTo[barr,arr];
];
AppendTo[Bmats,CleanMatrix@barr];
];
Ys=Table[Table[If[ll>=lmins[[ss+3]],SpinWeightedSphericalHarmonicY[ss,ll,mm,\[Theta],0],0],{ll,0,lmax}],{ss,-2,2}];
(* I want the Chandrasekhar eigenvalues, that is, the ones for -s. *)
eigs=Table[Table[If[ll>=lmins[[ss+3]],SpinWeightedSpheroidalEigenvalue[-Abs[ss],ll,mm,a\[Omega]0],0],{ll,0,lmax}],{ss,-2,2}];
Bmatm2=Bmats[[1]];Bmatm1=Bmats[[2]];Bmat0=Bmats[[3]];Bmatp1=Bmats[[4]]; Bmatp2=Bmats[[5]];
Ym2s=Ys[[1]];Ym1s=Ys[[2]];Y0s=Ys[[3]];Yp1s=Ys[[4]];Yp2s=Ys[[5]];
eigs0=eigs[[3]];
eigs1=eigs[[4]];eigs2=eigs[[5]];
(* Test whether the expansion appears to be working. *)
If[a\[Omega]0!=0,
\[Theta]try=0.9; stry=2;
t0=Bmats[[stry+3]] . Ys[[stry+3]]/.{\[Theta]->\[Theta]try};
t1=Join[Table[0,{ll,0,lmins[[stry+3]]-1}],Table[SpinWeightedSpheroidalHarmonicS[stry,ll,mm,a\[Omega]0,Method->{"SphericalExpansion","NumTerms"->nterms}][\[Theta]try,0],{ll,lmins[[stry+3]],lmax}]];
t0-t1
];
(* Now we need to calculate the arrays of mixing coefficients. *)
mixsin[s1_,s2_,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m,t0=0},
If[(s1==1)&&(s2==0),
t0=-Sqrt[2*(2*l2+1)/(2*l1+1)]*ClebschGordan[{l2,mm},{1,0},{l1,mm}]*ClebschGordan[{l2,-1},{1,1},{l1,0}];
];
If[s1==-1&&s2==0,
t0=Sqrt[2*(2*l2+1)/(2*l1+1)]*ClebschGordan[{l2,mm},{1,0},{l1,mm}]*ClebschGordan[{l2,1},{1,-1},{l1,0}];
];
If[s1==2&&s2==0,
t0=Sqrt[8*(2*l2+1)/(3*(2*l1+1))]*ClebschGordan[{l2,mm},{2,0},{l1,mm}]*ClebschGordan[{l2,-2},{2,2},{l1,0}];
];
If[s1==-2&&s2==0,
t0=Sqrt[8*(2*l2+1)/(3*(2*l1+1))]*ClebschGordan[{l2,mm},{2,0},{l1,mm}]*ClebschGordan[{l2,2},{2,-2},{l1,0}];
];
If[s1==2&&s2==1,
t0=-Sqrt[2*(2*l2+1)/((2*l1+1))]*ClebschGordan[{l2,mm},{1,0},{l1,mm}]*ClebschGordan[{l2,-2},{1,1},{l1,-1}];
];
If[s1==-2&&s2==-1,
t0=Sqrt[2*(2*l2+1)/((2*l1+1))]*ClebschGordan[{l2,mm},{1,0},{l1,mm}]*ClebschGordan[{l2,2},{1,-1},{l1,1}];
];

(* We also need the s1==-1 s2==+1 case here *)
If[s1==-1&&s2==1,
t0=Sqrt[8*(2*l2+1)/(3*(2*l1+1))]*ClebschGordan[{l2,mm},{2,0},{l1,mm}]*ClebschGordan[{l2,1},{2,-2},{l1,-1}];
];
(* Done on 20/11/22 late night, could be checked one more time. *)

(* Use the symmetry. *)
If[s1==0&&s2==1,t0=mixsin[s2,s1,l+k,-k,m];];
If[s1==0&&s2==-1,t0=mixsin[s2,s1,l+k,-k,m];];
If[s1==0&&s2==-2,t0=mixsin[s2,s1,l+k,-k,m];];
If[s1==0&&s2==2,t0=mixsin[s2,s1,l+k,-k,m];];
If[s1==1&&s2==-1,t0=mixsin[s2,s1,l+k,-k,m];];
If[s1==1&&s2==2,t0=mixsin[s2,s1,l+k,-k,m];];
If[s1==-1&&s2==-2,t0=mixsin[s2,s1,l+k,-k,m];];
Simplify[t0]
];
mixcos[s_,n_,l_,k_,m_]:=Module[{l1=l+k,l2=l,mm=m,ss=s,t0=0},
If[n==1,
t0=Sqrt[(2*l2+1)/(2*l1+1)]*ClebschGordan[{l2,mm},{1,0},{l1,mm}]*ClebschGordan[{l2,-ss},{1,0},{l1,-ss}];
];
If[n==2,
t0=1/3*KroneckerDelta[l1,l2]+2/3*Sqrt[(2*l2+1)/(2*l1+1)]*ClebschGordan[{l2,mm},{2,0},{l1,mm}]*ClebschGordan[{l2,-ss},{2,0},{l1,-ss}];
];
Simplify[t0]
];
SinMixMatrix[s1_,s2_,lmax_,m_]:=Module[{band,bands},
bands={};
For[k=-2,k<=2,k++,
band=Table[mixsin[s1,s2,ll,k,mm],{ll,Max[-k,0],Min[lmax,lmax-k]}];
AppendTo[bands,band];
];
SparseArray[{Band[{1,1}]->bands[[3]],Band[{2,1}]->bands[[2]],Band[{3,1}]->bands[[1]],Band[{1,2}]->bands[[4]],Band[{1,3}]->bands[[5]]}, {lmax+1,lmax+1}]
];
CosMixMatrix[s_,n_,lmax_,m_]:=Module[{band,bands},
bands={};
For[k=-2,k<=2,k++,
band=Table[mixcos[s,n,ll,k,mm],{ll,Max[-k,0],Min[lmax,lmax-k]}];
AppendTo[bands,band];
];
SparseArray[{Band[{1,1}]->bands[[3]],Band[{2,1}]->bands[[2]],Band[{3,1}]->bands[[1]],Band[{1,2}]->bands[[4]],Band[{1,3}]->bands[[5]]}, {lmax+1,lmax+1}]
];
(* Now define the S matrices, following the notes *)
\[Lambda]hat=DiagonalMatrix[Table[If[ll<lmins0,0,Sqrt[ll*(ll+1)]],{ll,0,lmax}]];
\[Lambda]2hat=DiagonalMatrix[Table[If[ll<lmins1,0,Sqrt[(ll-1)*(ll+2)]],{ll,0,lmax}]];
\[CapitalLambda]hat=\[Lambda]hat . \[Lambda]2hat;
signmat=DiagonalMatrix[Table[If[ll<lmins0,0,(-1)^(ll+mm)],{ll,0,lmax}]];
(*
ltry=4;
L1dagSm1=Bmatm1.(-\[Lambda]hat+a\[Omega]0*SinMixMatrix[-1,0,lmax,mm]).Y0s;
L1opSp1=Bmatp1.(\[Lambda]hat-a\[Omega]0*SinMixMatrix[1,0,lmax,mm]).Y0s;
L1dagSm1alt=Ldag[1,SpinWeightedSpheroidalHarmonicS[-1,ltry,mm,a\[Omega]0][\[Theta],0]]/.awsubs;
L1opSp1alt=Lop[1,SpinWeightedSpheroidalHarmonicS[1,ltry,mm,a\[Omega]0][\[Theta],0]]/.awsubs;
Plot[{L1dagSm1[[ltry+1]],L1dagSm1alt},{\[Theta],0,\[Pi]}]
Plot[L1dagSm1[[ltry+1]]-L1dagSm1alt,{\[Theta],0,\[Pi]}]
Plot[L1opSp1[[ltry+1]]-L1opSp1alt,{\[Theta],0,\[Pi]}]
*)
Get\[CapitalLambda][l_,m_,a\[Omega]_]:=Module[{AA},
SpinWeightedSpheroidalEigenvalue[-2,l,m,a\[Omega]]
]; (* N.B. s = -2. *)
(* N.B. I have multiplied calD by -1 because we use the -+++ signature. The results for mode functions and fluxes then agree with the BHPT results using TeukolskyPointParticleMode.
See check_sourceterms_with_BHPtoolkit.nb. *)
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
(* Calculate fluxes and make a comparison with BHP Toolkit. *)
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
(* Flux tests *)
ltry=Max[2,Abs[mm]];
orbit=KerrGeoOrbit[SetPrecision[a/.awsubs,prec],SetPrecision[r0/.awsubs,prec],0,1];
{EE-orbit["Energy"],LL-orbit["AngularMomentum"]}/.mysubs;
mode=TeukolskyPointParticleMode[-2,ltry,mm,0,0,orbit];
mode["Fluxes"]["Energy"]//N;
(* Now calculate using my equations. *)
Rfn=TeukolskyRadial[-2,ltry,mm,a/.awsubs,\[Omega]/.awsubs];
eq1=Aup*Rfn["Up"][r0]-Ain*Rfn["In"][r0]==(Pm2j0[ltry]/.s2jumps);
eq2=Aup*Rfn["Up"]'[r0]-Ain*Rfn["In"]'[r0]==(Pm2j1[ltry]/.s2jumps);
AinAupsol=NSolve[{eq1,eq2},{Ain,Aup}, WorkingPrecision->prec]//First;
EFluxInf=1/(4*\[Pi]*\[Omega]^2)*Abs[Aup/4]^2/.AinAupsol/.awsubs;
(EFluxInf-mode["Fluxes"]["Energy"]["\[ScriptCapitalI]"])//N;
(* Now calculate using the +2 Teukolsky functions *)
Rfn=TeukolskyRadial[2,ltry,mm,a0,\[Omega]0];
tmpPpup=Rfn["Up"][r0]*\[CapitalDelta][r]^2/.\[CapitalDelta]Kr0subs/.{r->r0};
tmpPpin=Rfn["In"][r0]*\[CapitalDelta][r]^2/.\[CapitalDelta]Kr0subs/.{r->r0};
tmpPpupd=(Rfn["Up"]'[r0]*\[CapitalDelta][r]^2+D[\[CapitalDelta][r]^2,r]*Rfn["Up"][r0])/.\[CapitalDelta]Kr0subs/.{r->r0};
tmpPpind=(Rfn["In"]'[r0]*\[CapitalDelta][r]^2+D[\[CapitalDelta][r]^2,r]*Rfn["In"][r0])/.\[CapitalDelta]Kr0subs/.{r->r0};
eq1=Aup*tmpPpup-Ain*tmpPpin==(Pp2j0[ltry]/.s2jumps);
eq2=Aup*tmpPpupd-Ain*tmpPpind==(Pp2j1[ltry]/.s2jumps);
AinAupsol=NSolve[{eq1,eq2},{Ain,Aup}, WorkingPrecision->prec]//First;
(* Now consider the Teuksolsky-Starobinskii identities *)
Csqsubs={Csq->AA^2+144*M^2*\[Omega]^2}/.{AA->Sqrt[\[CapitalLambda]^2*(\[CapitalLambda]+2)^2-8*\[Omega]^2*\[CapitalLambda]*(\[Alpha]sq*(5*\[CapitalLambda]+6)-12*a^2)+144*\[Omega]^4*\[Alpha]sq^2]}/.{\[Alpha]sq->a^2-a*m/\[Omega]}/.awsubs;
EFluxInf2=(1/(2*\[Pi])*8*\[Omega]^6*Abs[(Aup/.AinAupsol)]^2/Csq)/.Csqsubs/.{\[CapitalLambda]->Get\[CapitalLambda][ltry,mm,a\[Omega]0]}/.awsubs;
EFluxInf-EFluxInf2;
(EFluxInf-EFluxInf2)/EFluxInf;
(* Now for the jumps in the derivative of the trace. *)
If[a\[Omega]0==0,
s0jumps=Table[hj1[ll]->(-16*\[Pi]/(ut*\[CapitalDelta][r])*S0[\[Pi]/2]/.{S0[\[Pi]/2]->SpinWeightedSpheroidalHarmonicS[0,ll,mm,a\[Omega]0][\[Pi]/2,0]})/.\[CapitalDelta]Kr0subs/.{r->r0}/.orbitsubs,{ll,lmins0,lmax}];
,
s0jumps=Table[hj1[ll]->CleanValue[(-16*\[Pi]/(ut*\[CapitalDelta][r])*S0[\[Pi]/2]/.{S0[\[Pi]/2]->SpinWeightedSpheroidalHarmonicS[0,ll,mm,a\[Omega]0,Method->{"SphericalExpansion","NumTerms"->nterms}][\[Pi]/2,0]})/.\[CapitalDelta]Kr0subs/.{r->r0}/.orbitsubs],{ll,lmins0,lmax}];
];
s0jumps;
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
\[Kappa]mat//TableForm ;
(* Awkwardness: to keep the expressions compact, I would like to replace 'r' with the numerical value 'r0', but not in the radial functions themselves, as these I will use as unknowns. So I will replace the radial functions with jump expressions at an early stage. *)
Pp2jump=Table[If[ll>=lmins2,{Pp2[ll][r]->Pp2j0[ll],Pp2[ll]'[r]->Pp2j1[ll]},{}],{ll,0,lmax}]//Flatten;
Pm2jump=Table[If[ll>=lmins2,{Pm2[ll][r]->Pm2j0[ll],Pm2[ll]'[r]->Pm2j1[ll]},{}],{ll,0,lmax}]//Flatten;
Pp1jump=Table[If[ll>=lmins1,{Pp1[ll][r]->Pp1j0[ll],Pp1[ll]'[r]->Pp1j1[ll]},{}],{ll,0,lmax}]//Flatten;
Pm1jump=Table[If[ll>=lmins1,{Pm1[ll][r]->Pm1j0[ll],Pm1[ll]'[r]->Pm1j1[ll]},{}],{ll,0,lmax}]//Flatten;
hjump=Table[If[ll>=lmins0,{h[ll][r]->0,h[ll]'[r]->hj1[ll]},{}],{ll,0,lmax}]//Flatten;
\[Kappa]jump=Table[If[ll>=lmins0,{\[Kappa][ll][r]->\[Kappa]j0[ll],\[Kappa][ll]'[r]->\[Kappa]j1[ll]},{}],{ll,0,lmax}]//Flatten;
HeadFn[head_,ll_]:={head[r]->head[ll][r],head'[r]->head[ll]'[r],head''[r]->head[ll]''[r],head'''[r]->head[ll]'''[r],head''''[r]->head[ll]''''[r]};
HeadFn[Pm2,3];
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
Timing[
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
Timing[
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
Timing[
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
Timing[
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
(* This is the part that takes a long time. Remove. *)
(*
Timing[
d\[Rho]hs2lpmp=D[\[Rho]hs2lpmp,r]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs;
d\[Rho]chs2lpmm=D[\[Rho]chs2lpmm,r]/.Pm2lsubs/.\[CapitalDelta]Ksimp/.awsubs;
]
*)
(* Spin 0 *)
(* When we apply \[Kappa]lsubs, we also should fill in the \[Gamma]ll values *)
Timing[
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
Timing[
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
Timing[
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
Timing[
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
Timing[
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
Timing[
\[Rho]hs0j0lpmp=(\[Rho]hs0lpmp)/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
\[Rho]hs1j0lpmp=(\[Rho]hs1lpmp)/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
\[Rho]hs2j0lpmp=(\[Rho]hs2lpmp)/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
];
Timing[
\[Rho]chs0j0lpmm=(\[Rho]chs0lpmm)/.hlsubs/.\[Kappa]lsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
\[Rho]chs1j0lpmm=(\[Rho]chs1lpmm)/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
\[Rho]chs2j0lpmm=(\[Rho]chs2lpmm)/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
];
(* These parts take a long time. No longer necessary, even if m = 1. *)
(*
If[mm==1,
Timing[
(* Jumps in the radial derivative. *)
\[Rho]chs0j1lpmm=(D[\[Rho]chs0lpmm,r])/.hlsubs/.\[Kappa]lsubs/.\[Gamma]llsubs/.hjump/.\[Kappa]jump/.\[CapitalDelta]Kr0subs/.{r->r0}/.awsubs//Simplify;
\[Rho]chs1j1lpmm=(D[\[Rho]chs1lpmm,r])/.Pm1lsubs/.Pm1jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
];
Timing[
\[Rho]chs2j1lpmm=(D[\[Rho]chs2lpmm,r])/.Pm2lsubs/.Pm2jump/.\[CapitalDelta]Kr0subs/.{r->r0}//Simplify;
]
]
*)
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

\[Rho]chj0lpmm=\[Rho]chs0j0lpmm+\[Rho]chs1j0lpmm+\[Rho]chs2j0lpmm;
\[Rho]hj0lpmp=\[Rho]hs0j0lpmp+\[Rho]hs1j0lpmp+\[Rho]hs2j0lpmp;

(*
If[mm==1,
\[Rho]chj1lpmm=\[Rho]chs0j1lpmm+\[Rho]chs1j1lpmm+\[Rho]chs2j1lpmm-source\[Rho]chlpmm;
];
*)

eqs=Join[Table[hj0lplp[[ll+1]],{ll,lmins0,lmax}],Table[hj1lplp[[ll+1]],{ll,lmins0,lmax}],Table[hj0mpmp[[ll+1]],{ll,lmins2,lmax}],Table[hj1mpmp[[ll+1]],{ll,lmins2,lmax}]]/.s2jumps/.s0jumps;
Length[unknowns];
Length[eqs];
brhs=Table[CleanValue[-eqs[[kk]]/.unknownstozero],{kk,1,Length[eqs]}];
(* For Abs[mm] != 1, the number of unknowns should equal the number of equations. However, for Abs[mm]=1, there are two more unknowns than equations, and this case must be handled separately, using the spin-1 component of the metric perturbation. *)
If[Abs[mm]==1,
ll=1;
(*eqs=Join[eqs,{\[Rho]chj0lpmm[[ll+1]],\[Rho]chj0lpmm[[ll+2]]}]/.s2jumps/.s0jumps;
brhs=Table[CleanValue[-eqs[[kk]]/.unknownstozero],{kk,1,Length[eqs]}]//N;
*)
eqs=Join[eqs,{hj0lmlm[[ll+1]],hj1lmlm[[ll+1]]}]/.s2jumps/.s0jumps;
brhs=Table[CleanValue[-eqs[[kk]]/.unknownstozero],{kk,1,Length[eqs]}]//N;

];
brhs;
Mmat=CleanMatrix[Table[Coefficient[eqs[[j]],unknowns[[k]]],{j,1,Length[eqs]},{k,1,Length[eqs]}]];
Det[Mmat];
Quiet[sol=LinearSolve[Mmat,brhs]];
Simplify[Mmat . sol-brhs];
knowns=Table[unknowns[[k]]->CleanValue@sol[[k]],{k,1,Length[unknowns]}];
(* Make a table of the derived unknowns for inclusion in a tex file. *)
ldisplay=10;
header={"l","Pm1j0","Pm1j1","\[Kappa]j0","\[Kappa]j1"};
nsf=8;
mytable=Table[Join[{ll},Map[NumberForm[#,nsf]&,{Pm1j0[ll],Pm1j1[ll],Re[\[Kappa]j0[ll]],Re[\[Kappa]j1[ll]]}]],{ll,Max[Abs[mm],1],ldisplay}]/.knowns;
mygrid=Grid[Join[{header},mytable],Frame->All,Alignment->Right];
TeXForm[mygrid];
(* Also make a table of the knowns (from Weyl scalar and trace) for inclusion in a tex file. *)
header={"l","Pm2j0","Pm2j1","h1"};
mytable=Table[Join[{ll},Map[NumberForm[#,nsf]&,{Pm2j0[ll],Pm2j1[ll],Re[hj1[ll]]}]],{ll,Max[Abs[mm],1],ldisplay}]/.s2jumps/.s0jumps;
mygrid=Grid[Join[{header},mytable],Frame->All,Alignment->Right];
TeXForm[mygrid];
(* Check the lplp and mpmp components. *)
Erhceck1 = hj0lplp/.s2jumps/.s0jumps/.knowns//Abs;
Erhceck2 = hj0mpmp/.s2jumps/.s0jumps/.knowns//Abs;
Erhceck3 = hj1lplp/.s2jumps/.s0jumps/.knowns//Abs;
Erhceck4 = hj1mpmp/.s2jumps/.s0jumps/.knowns//Abs;
(* Check the lmlm and mmmm components. *)
Erhceck5 = hj0lmlm/.s2jumps/.s0jumps/.knowns//Abs;
Erhceck6 = hj1lmlm/.s2jumps/.s0jumps/.knowns//Abs;
Erhceck7 = hj0mmmm/.s2jumps/.s0jumps/.knowns//Abs;
Erhceck8 = hj1mmmm/.s2jumps/.s0jumps/.knowns//Abs;
(* Check the lpmp and lpmm components. *)
Erhceck9  = \[Rho]hj0lpmp/.s2jumps/.s0jumps/.knowns//Simplify//Abs;
Erhceck10 = \[Rho]chj0lpmm/.s2jumps/.s0jumps/.knowns//Simplify//Abs;
Erhceck = {Erhceck1,Erhceck2,Erhceck3,Erhceck4,Erhceck5,Erhceck6,Erhceck7,Erhceck8,Erhceck9,Erhceck10};
Export[directory<>"JumpChecks/jumpchecks"<>ToString[iConfig]<>".dat",Erhceck];
(*
If[mm==1,
\[Rho]chj1lpmm/.s2jumps/.s0jumps/.knowns//Simplify//Abs
]
*)
header={"l","Pm2j0","Pm2j1","h1","Pm1j0","Pm1j1","\[Kappa]j0","\[Kappa]j1"};
mytable=Table[Join[{ll},{Pm2j0[ll],Pm2j1[ll],Re[hj1[ll]],Pm1j0[ll],Pm1j1[ll],Re[\[Kappa]j0[ll]],Re[\[Kappa]j1[ll]]}],{ll,Max[Abs[mm],1],lmax}]/.knowns/.s2jumps/.s0jumps;
mygrid=Grid[Join[{header},mytable],Frame->All,Alignment->Right];
Export[directory<>"data/jumps_"<>ToString[iConfig]<>".mx",mytable];
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
(*
Plot[{Re[h0[r]/.dsol/.{r->rr}],Im[h0[r]/.dsol/.{r->rr}],Re[h0infser/.\[Lambda]0subs/.awsubs/.{z->1/rr,r->rr}],Im[h0infser/.\[Lambda]0subs/.awsubs/.{z->1/rr,r->rr}]},{rr,rinf-200,rinf}]
LogPlot[Abs[(h0[r]/.dsol/.{r->rr})-(h0infser/.\[Lambda]0subs/.awsubs/.{z->1/rr,r->rr})],{rr,rinf-200,rinf}]
Plot[{Re[h0[r]/.dsol/.{r->rr}],Re[h0up[rr]]},{rr,rinf-200,rinf},PlotStyle->{Dashed,Dotted}]
Plot[{Re[h0[r]/.dsol/.{r->rr}],Re[h0up[rr]]},{rr,r0,r0+30},PlotStyle->{Dashed,Dotted}]
Plot[Abs[(h0[r]/.dsol/.{r->rr})-(h0up[rr])],{rr,r0,r0+20}]
(* Check that the series solution for the \[Kappa] function is in good agreement with the numerical solution from the ODE, in the appropriate large-r domain. *)
Plot[{Re[\[Kappa]r0[r]/.dsol/.{r->rr}],Im[\[Kappa]r0[r]/.dsol/.{r->rr}]},{rr,rinf-200,rinf}]
Plot[{Re[\[Kappa]r0[r]/.dsol/.{r->rr}],Im[\[Kappa]r0[r]/.dsol/.{r->rr}],Re[\[Kappa]0infser/.\[Lambda]0subs/.awsubs/.{z->1/rr,r->rr}],Im[\[Kappa]0infser/.\[Lambda]0subs/.awsubs/.{z->1/rr,r->rr}]},{rr,rinf-200,rinf}]
Plot[{Re[\[Kappa]r2[r]/.dsol/.{r->rr}],Im[\[Kappa]r2[r]/.dsol/.{r->rr}],Re[\[Kappa]2infser/.\[Lambda]0subs/.awsubs/.{z->1/rr,r->rr}],Im[\[Kappa]2infser/.\[Lambda]0subs/.awsubs/.{z->1/rr,r->rr}]},{rr,rinf-200,rinf}]
*)
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
header={"l","\[Alpha]s2co","\[Alpha]s1co","\[Alpha]s0hco","\[Alpha]s0\[Kappa]co"};
upintext={"up","in"};
For[kk=1,kk<=2,kk++,
mytable=Table[{ll,\[Alpha]s2cos[[ll+1]][[kk]],\[Alpha]s1cos[[ll+1]][[kk]],\[Alpha]s0hcos[[ll+1]][[kk]],\[Alpha]s0\[Kappa]cos[[ll+1]][[kk]]},{ll,0,lmax}];
mygrid=Grid[Join[{header},mytable],Frame->All,Alignment->Right];
Export[directory<>"data/radial_coeffs_"<>upintext[[kk]]<>"_"<>ToString[iConfig]<>".mx",mytable];
];
Print["4. Prepare grid and mode functions."];
TortoiseCoord[r_]:=
SetPrecision[r+(rp+rm)/(rp-rm)*(rp*Log[(r-rp)/2]-rm*Log[(r-rm)/2])/.rprmvals,prec];
InvTortoise[x_?NumericQ]:=r/.FindRoot[TortoiseCoord[r]-x,{r,SetPrecision[(rp+SetPrecision[0.0000001,prec])/.rprmvals,prec]},WorkingPrecision->prec];

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
filen=directory<>"data/lm_rs_"<>ToString[iConfig]<>".dat";

Export[filen,Transpose@{rstars,rs}];
filen=directory<>"data/lm_qs_"<>ToString[iConfig]<>".dat";
Export[filen,qs];
filen=directory<>"data/lm_rsL_"<>ToString[iConfig]<>".dat";
Export[filen,Transpose@{rstarsL,rsL}];
filen=directory<>"data/lm_rsR_"<>ToString[iConfig]<>".dat";
Export[filen,Transpose@{rstarsR,rsR}];
TeukMethod={"NumericalIntegration","Domain"->{"In"->{rmin,r0},"Up"->{r0,rmax}}};
(* TeukMethod={"MST","RenormalizedAngularMomentum"\[Rule]"Monodromy"} *)
Timing[
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
Timing[
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
Timing[
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
Timing[
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
(* a\[Omega]0prec=SetPrecision[a\[Omega]0,50]; *)
sphermethod={"SphericalExpansion","NumTerms"->nterms};
SetOptions[SpinWeightedSpheroidalHarmonicS,Method->sphermethod];
GetSsubs[ltry_,\[Theta]try_?NumericQ]:=Module[{tbl1,tbl2,fns},
a\[Omega]0prec=SetPrecision[a\[Omega]0,2*prec]; 
(* Extended precision appears to be necessary to avoid "failed to converge" messages. *)
tbl1={Sp2[\[Theta]]->SpinWeightedSpheroidalHarmonicS[2,ltry,mm,a\[Omega]0prec][\[Theta]try,0],
Sp2'[\[Theta]]->(D[SpinWeightedSpheroidalHarmonicS[2,ltry,mm,a\[Omega]0prec][\[Theta],0],\[Theta]]/.{\[Theta]->\[Theta]try}),
Sm2[\[Theta]]->SpinWeightedSpheroidalHarmonicS[-2,ltry,mm,a\[Omega]0prec][\[Theta]try,0],
Sm2'[\[Theta]]->(D[SpinWeightedSpheroidalHarmonicS[-2,ltry,mm,a\[Omega]0prec][\[Theta],0],\[Theta]]/.{\[Theta]->\[Theta]try}),
Sp1[\[Theta]]->SpinWeightedSpheroidalHarmonicS[1,ltry,mm,a\[Omega]0prec][\[Theta]try,0],
Sp1'[\[Theta]]->(D[SpinWeightedSpheroidalHarmonicS[1,ltry,mm,a\[Omega]0prec][\[Theta],0],\[Theta]]/.{\[Theta]->\[Theta]try}),
Sm1[\[Theta]]->SpinWeightedSpheroidalHarmonicS[-1,ltry,mm,a\[Omega]0prec][\[Theta]try,0],
Sm1'[\[Theta]]->(D[SpinWeightedSpheroidalHarmonicS[-1,ltry,mm,a\[Omega]0prec][\[Theta],0],\[Theta]]/.{\[Theta]->\[Theta]try})
};
fns={Sin[\[Theta]],Cos[\[Theta]],Cot[\[Theta]],Csc[\[Theta]]};
tbl2=Map[#->(#/.\[Theta]->\[Theta]try)&,fns];
Join[tbl1,tbl2]
];
Timing[GetSsubs[20,SetPrecision[\[Pi]/4,prec]]];
GetS0subs[\[Theta]try_?NumericQ,lval_,\[Kappa]ord_:6]:=Module[{myfn,tbl1,tbl2,fns},
myfn=(Coefficient[#,S0[\[Theta]]]/.{\[Theta]->\[Theta]try})*S0[\[Theta]]+(Coefficient[#,S0'[\[Theta]]]/.{\[Theta]->\[Theta]try})*S0'[\[Theta]]&;
tbl1=Table[S0[ll]''[\[Theta]]->(Map[myfn,S0''[\[Theta]]/.S0subs]/.{S0[\[Theta]]->S0[ll][\[Theta]],S0'[\[Theta]]->S0[ll]'[\[Theta]]}/.{\[Lambda]0->eigs0[[ll+1]]}/.awsubs),
{ll,Max[lmins0,lval-\[Kappa]ord],Min[lmax,lval+\[Kappa]ord],2}];
tbl2=Table[
{S0[ll][\[Theta]]->(SpinWeightedSpheroidalHarmonicS[0,ll,mm,a\[Omega]0prec][\[Theta]try,0]),
S0[ll]'[\[Theta]]->((D[SpinWeightedSpheroidalHarmonicS[0,ll,mm,a\[Omega]0prec][\[Theta],0],\[Theta]]/.{\[Theta]->\[Theta]try}))},
{ll,Max[lmins0,lval-\[Kappa]ord],Min[lmax,lval+\[Kappa]ord],2}]//Flatten;
fns={Sin[\[Theta]],Cos[\[Theta]],Cot[\[Theta]],Csc[\[Theta]]};
tbl3=Map[#->(#/.\[Theta]->\[Theta]try)&,fns];
Join[tbl1/.tbl2,tbl2,tbl3]
];
Timing[GetS0subs[SetPrecision[\[Pi]/4,prec],2]];
EvaluateRHS[set_,rval_?NumericQ]:=Module[{kk},Table[Keys[set[[kk]]]->(Values[set[[kk]]]/.{r->rval}),{kk,1,Length[set]}]];
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
Timing[
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
Timing[
GetPsubsAll[1,SetPrecision[7.1,prec]];
];
(* Scalar part needs careful handling. *)
hppsubs=(h''[r]->(h0''[r]/.h0subs));
Get\[Kappa]subs[ltry_,rtry_?NumericQ]:=Module[{set,repl},
\[Kappa]row=\[Kappa][ltry][r]*S0[ltry][\[Theta]]+Sum[If[ll!=ltry,
a^2*h[ltry][r]*\[Gamma]mat[[ltry+1,ll+1]]*S0[ll][\[Theta]]/(2*(eigs0[[ltry+1]]-eigs0[[ll+1]]))
,0]
,{ll,Max[lmins0,ltry-\[Kappa]ord],Min[lmax,ltry+\[Kappa]ord]}]/.awsubs;
set={\[Kappa][r,\[Theta]]->\[Kappa]row,D[\[Kappa][r,\[Theta]],r]->D[\[Kappa]row,r],D[\[Kappa][r,\[Theta]],{r,2}]->D[\[Kappa]row,{r,2}],
D[\[Kappa][r,\[Theta]],\[Theta]]->D[\[Kappa]row,\[Theta]],D[\[Kappa][r,\[Theta]],{\[Theta],2}]->D[\[Kappa]row,{\[Theta],2}],D[D[\[Kappa][r,\[Theta]],r],\[Theta]]->D[D[\[Kappa]row,r],\[Theta]]}/.\[Kappa]lsubs/.hlsubs/.\[Gamma]llsubs/.\[CapitalDelta]Ksubs;
If[Greater[rtry-r0,0],repl=Join[\[Kappa]0uprepl,h0uprepl],repl=Join[\[Kappa]0inrepl,h0inrepl]];
tmp1=Table[Keys[set[[kk]]]->((Values[set[[kk]]]/.repl)/.awsubs/.{r->rtry}),{kk,1,Length[set]}];
tmp2=If[Greater[rtry-r0,0],EvaluateRHS[(HeadFn[h,ltry][[1;;2]]/.h0uprepl),rtry],EvaluateRHS[(HeadFn[h,ltry][[1;;2]]/.h0inrepl),rtry]];
tmp3={h''[r]->(h0''[r]/.h0subs/.\[CapitalDelta]Ksubs/.awsubs)}/.{h0'[r]->h[ltry]'[r],h0[r]->h[ltry][r]};
tmp4=If[Greater[rtry-r0,0],EvaluateRHS[tmp3/.h0uprepl,rtry],EvaluateRHS[tmp3/.h0inrepl,rtry]]/.{\[Lambda]0->eigs0[[ltry+1]]};
Join[tmp2,tmp1,tmp4]
];
Timing[Get\[Kappa]subs[2,SetPrecision[3.7,prec]];];
GetSpin0subs[rtry_?NumericQ]:=If[Greater[rtry-r0,0],
Join[EvaluateRHS[h0uprepl,rtry],EvaluateRHS[\[Kappa]0uprepl,rtry]],
Join[EvaluateRHS[h0inrepl,rtry],EvaluateRHS[\[Kappa]0inrepl,rtry]]
];
Timing[
GetSpin0subs[SetPrecision[7.1,prec]];
];
GetAll\[Kappa]subs[rtry_?NumericQ]:=Table[Get\[Kappa]subs[ll,rtry],{ll,Abs[mm],lmax}]//Flatten;
Timing[GetAll\[Kappa]subs[SetPrecision[3.7,prec]];];
Print["A. Fill 1D grid."];
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
(* Evals1Fn[X_,rtry1_]:=X/.\[CapitalDelta]Ktbl/.awsubs/.myP1subs/.{r->rtry1}; *)
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
rtry=SetPrecision[7.1,prec];
\[CapitalDelta]Ktbl=EvaluateRHS[\[CapitalDelta]Ksubs,rtry]/.awsubs;
Timing[
tmp1=Calcs1components[rtry];
];
Timing[
myPsubs=GetPsubsAll[1,rtry];
tmp2={hraws1lplp,hraws1lmlm,hraws1mpmp,hraws1mmmm,\[Rho]hraws1lpmp,\[Rho]chraws1lpmm,\[Rho]chraws1lmmp,\[Rho]hraws1lmmm,\[CapitalSigma]\[CapitalDelta]hraws1lplm}/.awsubs/.\[CapitalDelta]Ktbl/.myPsubs/.{r->rtry};
];
(* Factor of 10 speed-up. *)
tmp1-tmp2;
idmat=DiagonalMatrix[Table[If[ll>=lmins0,1,0],{ll,0,lmax}]];
t0=Bmatp2 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]2hat . SinMixMatrix[1,0,lmax,mm]+a\[Omega]0^2*SinMixMatrix[2,0,lmax,mm]);
t1=signmat . Bmatm2 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]2hat . SinMixMatrix[-1,0,lmax,mm]+a\[Omega]0^2*SinMixMatrix[-2,0,lmax,mm]);
Ss2lplp=(t0+t1)/.awsubs;
Ss2lmlm=signmat . Ss2lplp/.awsubs; 
Ss2mpmp=Bmatp2/.awsubs;
Ss2mmmm=signmat . Bmatm2/.awsubs;
SSplus=(Bmatp2 . (\[CapitalLambda]hat-2*a*\[Omega]*\[Lambda]2hat . SinMixMatrix[1,0,lmax,mm]+a^2*\[Omega]^2*SinMixMatrix[2,0,lmax,mm]))/.awsubs;
SSminus=(Bmatm2 . (\[CapitalLambda]hat-2*a*\[Omega]*\[Lambda]2hat . SinMixMatrix[-1,0,lmax,mm]+a^2*\[Omega]^2*SinMixMatrix[-2,0,lmax,mm]))/.awsubs;
(* Attach an a-prefix before they have been evaluated numerically. *)
aPp2vec=Pp2vec/.\[CapitalDelta]Ksimp;
aPm2vec=Pm2vec/.\[CapitalDelta]Ksimp;
aDdagPp2vec=Map[Ddag,aPp2vec]/.\[CapitalDelta]Ksimp;
aDdagDdagPp2vec=Map[Ddag,aDdagPp2vec]/.\[CapitalDelta]Ksimp;
aDdagDdagDdagPp2vec=Map[Ddag,aDdagDdagPp2vec]/.\[CapitalDelta]Ksimp;
aDPm2vec=Map[Dop,aPm2vec]/.\[CapitalDelta]Ksimp;
aDDPm2vec=Map[Dop,aDPm2vec]/.\[CapitalDelta]Ksimp;
aDDDPm2vec=Map[Dop,aDDPm2vec]/.\[CapitalDelta]Ksimp;
aDopDdagDdagPp2vec=Map[Dop,aDdagDdagPp2vec]/.\[CapitalDelta]Ksimp;
aDdagDDPm2vec=Map[Ddag,aDDPm2vec]/.\[CapitalDelta]Ksimp;
alplmterm1A=Map[Dop,\[CapitalDelta][r]*Map[Ddag,aDDPm2vec]]+Map[Ddag,\[CapitalDelta][r]*Map[Dop,aDDPm2vec]];
alplmterm1B=Map[D[#,r]&,aDDPm2vec];
alplmterm2A=Map[Dop,\[CapitalDelta][r]*Map[Ddag,aDdagDdagPp2vec]]+Map[Ddag,\[CapitalDelta][r]*Map[Dop,aDdagDdagPp2vec]];
alplmterm2B=Map[D[#,r]&,aDdagDdagPp2vec];
Clear[cosmix,sinmix];
sinp1=SinMixMatrix[0,1,lmax,mm];
sinm1=SinMixMatrix[0,-1,lmax,mm];
sinmix1m1a=SinMixMatrix[-1,1,lmax,mm];
sinmix1m1b=SinMixMatrix[1,-1,lmax,mm];
\[Rho]mat1=r*idmat+I*a*CosMixMatrix[1,1,lmax,mm];  (* These should be used for the first two components: l+m+ and l+m- *)
\[Rho]cmat1=r*idmat-I*a*CosMixMatrix[-1,1,lmax,mm];
\[Rho]cmat2=r*idmat-I*a*CosMixMatrix[1,1,lmax,mm]; (* These should be used for the second components: l-m+ and l-m- *)
\[Rho]mat2=r*idmat+I*a*CosMixMatrix[-1,1,lmax,mm];
cosmix01=CosMixMatrix[0,1,lmax,mm];
cosmix02=CosMixMatrix[0,2,lmax,mm];
\[Rho]cmat3=r*idmat-I*a*cosmix01;  (* These should be used for l+l- *) 
\[Rho]mat3=r*idmat+I*a*cosmix01;
\[CapitalSigma]mat=r^2*idmat+a^2*cosmix02;
sinp10=SinMixMatrix[1,0,lmax,mm];
sinm10=SinMixMatrix[-1,0,lmax,mm];
MMp2=Bmatp2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[2,1,lmax,mm]);
MMm2=Bmatm2 . (\[Lambda]2hat-a*\[Omega]*SinMixMatrix[-2,-1,lmax,mm]);
(* Pre-define some more matrices for extra speed up. *)
(* This seems to have slowed it down! *)
(*
Timing[
With[{sinmix=sinp1,sinmix1m1=sinmix1m1a,\[Rho]mat=\[Rho]mat1},
s2lpmpM1=signmat.MMm2.((\[Lambda]hat.\[Lambda]hat-2*a*\[Omega]*\[Lambda]hat.sinmix+a^2*\[Omega]^2*sinmix1m1).\[Rho]mat.\[Rho]mat-2*I*a*(\[Lambda]hat.sinmix-a*\[Omega]*sinmix1m1).\[Rho]mat-2*a^2*sinmix1m1)/.awsubs;
];
With[{sinmix=sinm1,sinmix1m1=sinmix1m1b,\[Rho]cmat=\[Rho]cmat1},
s2lpmmM1=MMp2.((\[Lambda]hat.\[Lambda]hat-2*a*\[Omega]*\[Lambda]hat.sinmix+a^2*\[Omega]^2*sinmix1m1).\[Rho]cmat.\[Rho]cmat-2*I*a*(\[Lambda]hat.sinmix-a*\[Omega]*sinmix1m1).\[Rho]cmat-2*a^2*sinmix1m1)/.awsubs;
];
With[{sinmix=sinp1,sinmix1m1=sinmix1m1a,\[Rho]cmat=\[Rho]cmat2},
s2lmmpM1=MMm2.((\[Lambda]hat.\[Lambda]hat-2*a*\[Omega]*\[Lambda]hat.sinmix+a^2*\[Omega]^2*sinmix1m1).\[Rho]cmat.\[Rho]cmat+2*I*a*(\[Lambda]hat.sinmix-a*\[Omega]*sinmix1m1).\[Rho]cmat-2*a^2*sinmix1m1)/.awsubs;
];
With[{sinmix=sinm1,sinmix1m1=sinmix1m1b,\[Rho]mat=\[Rho]mat2},
s2lmmmM1=signmat.MMp2.((\[Lambda]hat.\[Lambda]hat-2*a*\[Omega]*\[Lambda]hat.sinmix+a^2*\[Omega]^2*sinmix1m1).\[Rho]mat.\[Rho]mat+2*I*a*(\[Lambda]hat.sinmix-a*\[Omega]*sinmix1m1).\[Rho]mat-2*a^2*sinmix1m1)/.awsubs;
];
]
*)
Calcs2components[rtry_?NumericQ]:=Module[{
tmp,tmp1,tmp2,tmp3,term1,term2,term3,term4,\[CapitalDelta]Ktbl,myP2subs,
\[Rho]mat,\[Rho]cmat,sinmix,sinmix1m1,
Pp2vec,Pm2vec,DdagPp2vec,
DdagDdagPp2vec,
DdagDdagDdagPp2vec,
DPm2vec,
DDPm2vec,
DDDPm2vec,
DopDdagDdagPp2vec,
DdagDDPm2vec,lplmterm1A,lplmterm1B,lplmterm2A,lplmterm2B,
hraws2lplp,hraws2lmlm,hraws2mpmp,hraws2mmmm,
\[Rho]hraws2lpmp,\[Rho]chraws2lpmm,\[Rho]chraws2lmmp,\[Rho]hraws2lmmm,\[CapitalSigma]\[CapitalDelta]hraws2lplm},
\[CapitalDelta]Ktbl=EvaluateRHS[\[CapitalDelta]Ksubs,rtry]/.awsubs;
myP2subs=GetPsubsAll[2,rtry];

{Pp2vec,Pm2vec,DdagPp2vec,
DdagDdagPp2vec,
DdagDdagDdagPp2vec,
DPm2vec,
DDPm2vec,
DDDPm2vec,
DopDdagDdagPp2vec,
DdagDDPm2vec}={aPp2vec,aPm2vec,aDdagPp2vec,
aDdagDdagPp2vec,
aDdagDdagDdagPp2vec,
aDPm2vec,
aDDPm2vec,
aDDDPm2vec,
aDopDdagDdagPp2vec,
aDdagDDPm2vec}/.\[CapitalDelta]Ktbl/.awsubs/.myP2subs/.{r->rtry};

{lplmterm1A,lplmterm1B,lplmterm2A,lplmterm2B}={alplmterm1A,alplmterm1B,alplmterm2A,alplmterm2B}/.\[CapitalDelta]Ktbl/.awsubs/.myP2subs/.{r->rtry};

(* l+ l+,  l- l- *)
hraws2lplp=-1/(6*\[Omega]^2*\[CapitalDelta][r]^2)*Pp2vec . Ss2lplp/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};
hraws2lmlm=-1/(6*\[Omega]^2*\[CapitalDelta][r]^2)*Pm2vec . Ss2lmlm/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* m+ m+,  m- m- *)
tmp=-1/(6*\[Omega]^2)*(DdagDdagPp2vec+signmat . DDPm2vec);
hraws2mpmp=tmp . Ss2mpmp/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};
hraws2mmmm=tmp . Ss2mmmm/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* l+ m+ *)
sinmix=sinp1;
sinmix1m1=sinmix1m1a;
\[Rho]mat=\[Rho]mat1/.{r->rtry}/.awsubs;
term1=1/(12*\[Omega]^2)*(DDPm2vec . signmat . SSplus-DdagDdagPp2vec . signmat . SSminus) . (\[Lambda]hat-a*\[Omega]*sinmix)/.awsubs;
term2=-1/(12*\[Omega]^2)*(DDDPm2vec . signmat . SSplus-DopDdagDdagPp2vec . signmat . SSminus) . ((\[Lambda]hat-a*\[Omega]*sinmix) . \[Rho]mat-I*a*sinmix)/.awsubs;
term3=-1/(3*I*\[Omega])*(DDDPm2vec . signmat . MMp2 . \[Rho]mat . \[Rho]mat-2*DDPm2vec . signmat . MMp2 . \[Rho]mat+2*DPm2vec . signmat . MMp2)/.awsubs;
tmp=(\[Lambda]hat . \[Lambda]hat-2*a*\[Omega]*\[Lambda]hat . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]mat . \[Rho]mat-2*I*a*(\[Lambda]hat . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]mat-2*a^2*sinmix1m1/.awsubs;
term4=-1/(3*I*\[Omega]*\[CapitalDelta][r])*DdagPp2vec . (signmat . MMm2 . tmp)/.awsubs;
(* term4=-1/(3*I*\[Omega]*\[CapitalDelta][r])*DdagPp2vec.s2lpmpM1;*)
tmp=term1+term2+term3+term4;
\[Rho]hraws2lpmp=(tmp/(-6*I*M*\[Omega]))/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* l+ m- *)
sinmix=sinm1;
sinmix1m1=sinmix1m1b;
\[Rho]cmat=\[Rho]cmat1/.{r->rtry}/.awsubs;
term1=-1/(12*\[Omega]^2)*(DDPm2vec . SSminus-DdagDdagPp2vec . SSplus) . (\[Lambda]hat-a*\[Omega]*sinmix)/.awsubs;
term2=1/(12*\[Omega]^2)*(DDDPm2vec . SSminus-DopDdagDdagPp2vec . SSplus) . ((\[Lambda]hat-a*\[Omega]*sinmix) . \[Rho]cmat-I*a*sinmix)/.awsubs;
term3=1/(3*I*\[Omega])*(DDDPm2vec . MMm2 . \[Rho]cmat . \[Rho]cmat-2*DDPm2vec . MMm2 . \[Rho]cmat+2*DPm2vec . MMm2)/.awsubs;
tmp=((\[Lambda]hat . \[Lambda]hat-2*a*\[Omega]*\[Lambda]hat . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]cmat . \[Rho]cmat-2*I*a*(\[Lambda]hat . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]cmat-2*a^2*sinmix1m1)/.awsubs;
term4=1/(3*I*\[Omega]*\[CapitalDelta][r])*DdagPp2vec . MMp2 . tmp/.awsubs;
(*term4=1/(3*I*\[Omega]*\[CapitalDelta][r])*DdagPp2vec.s2lpmmM1/.awsubs;*)
tmp=term1+term2+term3+term4;
\[Rho]chraws2lpmm=(tmp/(-6*I*M*\[Omega]))/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* l- m+ *)
sinmix=sinp1;
sinmix1m1=sinmix1m1a;
\[Rho]cmat=\[Rho]cmat2/.{r->rtry}/.awsubs;
term1=1/(12*\[Omega]^2)*(DDPm2vec . SSminus-DdagDdagPp2vec . SSplus) . (\[Lambda]hat-a*\[Omega]*sinmix)/.awsubs;
term2=-1/(12*\[Omega]^2)*(DdagDDPm2vec . SSminus-DdagDdagDdagPp2vec . SSplus) . ((\[Lambda]hat-a*\[Omega]*sinmix) . \[Rho]cmat+I*a*sinmix)/.awsubs;
term3=-1/(3*I*\[Omega])*(DdagDdagDdagPp2vec . MMp2 . \[Rho]cmat . \[Rho]cmat-2*DdagDdagPp2vec . MMp2 . \[Rho]cmat+2*DdagPp2vec . MMp2)/.awsubs;
tmp=((\[Lambda]hat . \[Lambda]hat-2*a*\[Omega]*\[Lambda]hat . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]cmat . \[Rho]cmat+2*I*a*(\[Lambda]hat . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]cmat-2*a^2*sinmix1m1)/.awsubs;
term4=-1/(3*I*\[Omega]*\[CapitalDelta][r])*DPm2vec . MMm2 . tmp/.awsubs;
(*term4=-1/(3*I*\[Omega]*\[CapitalDelta][r])*DPm2vec.s2lmmpM1/.awsubs;*)
tmp=term1+term2+term3+term4;
\[Rho]chraws2lmmp=(tmp/(-6*I*M*\[Omega]))/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* l- m- *)
sinmix=sinm1;
sinmix1m1=sinmix1m1b;
\[Rho]mat=\[Rho]mat2/.{r->rtry}/.awsubs;
term1=-1/(12*\[Omega]^2)*(DDPm2vec . signmat . SSplus-DdagDdagPp2vec . signmat . SSminus) . (\[Lambda]hat-a*\[Omega]*sinmix)/.awsubs;
term2=1/(12*\[Omega]^2)*(DdagDDPm2vec . signmat . SSplus-DdagDdagDdagPp2vec . signmat . SSminus) . ((\[Lambda]hat-a*\[Omega]*sinmix) . \[Rho]mat+I*a*sinmix)/.awsubs;
term3=1/(3*I*\[Omega])*(DdagDdagDdagPp2vec . signmat . MMm2 . \[Rho]mat . \[Rho]mat-2*DdagDdagPp2vec . signmat . MMm2 . \[Rho]mat+2*DdagPp2vec . signmat . MMm2)/.awsubs;
tmp=((\[Lambda]hat . \[Lambda]hat-2*a*\[Omega]*\[Lambda]hat . sinmix+a^2*\[Omega]^2*sinmix1m1) . \[Rho]mat . \[Rho]mat+2*I*a*(\[Lambda]hat . sinmix-a*\[Omega]*sinmix1m1) . \[Rho]mat-2*a^2*sinmix1m1)/.awsubs;
term4=1/(3*I*\[Omega]*\[CapitalDelta][r])*DPm2vec . (signmat . MMp2 . tmp)/.awsubs;
(*term4=1/(3*I*\[Omega]*\[CapitalDelta][r])*DPm2vec.s2lmmmM1/.awsubs;*)
tmp=term1+term2+term3+term4;
\[Rho]hraws2lmmm=(tmp/(-6*I*M*\[Omega]))/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

(* l+ l- *)
\[Rho]mat=\[Rho]mat3/.{r->rtry}/.awsubs;
\[Rho]cmat=\[Rho]cmat3/.{r->rtry}/.awsubs;
tmp3=\[Lambda]hat . (sinm10-sinp10);
term1=1/(24*\[Omega]^2)*(lplmterm1A . SSminus . \[CapitalSigma]mat-4*r*\[CapitalDelta][r]*lplmterm1B . SSminus-2*a^2*DDPm2vec . SSminus . tmp3 . cosmix01);
term2=-1/(24*\[Omega]^2)*(lplmterm2A . SSplus . \[CapitalSigma]mat-4*r*\[CapitalDelta][r]*lplmterm2B . SSplus-2*a^2*DdagDdagPp2vec . SSplus . tmp3 . cosmix01);
tmp=(2*a^2*r*cosmix01-I*a*(r^2*idmat+3*a^2*cosmix02));
term3=1/(3*I*\[Omega])*(DDPm2vec . SSminus . \[CapitalSigma]mat . \[Rho]cmat-DPm2vec . SSminus . \[Rho]cmat . \[Rho]cmat
-DDPm2vec . MMm2 . sinm10 . tmp
-2*I*a*DPm2vec . MMm2 . sinm10 . \[Rho]mat);
term4=1/(3*I*\[Omega])*(DdagDdagPp2vec . SSplus . \[CapitalSigma]mat . \[Rho]cmat-DdagPp2vec . SSplus . \[Rho]cmat . \[Rho]cmat
+DdagDdagPp2vec . MMp2 . sinp10 . tmp
+2*I*a*DdagPp2vec . MMp2 . sinp10 . \[Rho]mat);
tmp=term1+term2+term3+term4;
\[CapitalSigma]\[CapitalDelta]hraws2lplm=(tmp/(-6*I*M*\[Omega]))/.\[CapitalDelta]Ktbl/.awsubs/.{r->rtry};

{hraws2lplp,hraws2lmlm,hraws2mpmp,hraws2mmmm,\[Rho]hraws2lpmp,\[Rho]chraws2lpmm,\[Rho]chraws2lmmp,\[Rho]hraws2lmmm,\[CapitalSigma]\[CapitalDelta]hraws2lplm}
];
rtry=SetPrecision[7.1,prec];
Timing[tmp1=Calcs2components[rtry];];
Timing[
myPsubs=GetPsubsAll[2,rtry];
\[CapitalDelta]Ktbl=EvaluateRHS[\[CapitalDelta]Ksubs,rtry]/.awsubs;
tmp2={hraws2lplp,hraws2lmlm,hraws2mpmp,hraws2mmmm}/.\[CapitalDelta]Ktbl/.awsubs/.myPsubs/.{r->rtry};
];
tmp1[[1;;4]]-tmp2;
(* Check that the s2 functions give the same result as the raw components. *)
rtry=SetPrecision[7.1,prec];
\[CapitalDelta]Ktbl=EvaluateRHS[\[CapitalDelta]Ksubs,rtry]/.awsubs;
Timing[
tmp1=Calcs2components[rtry];
];
(*
Timing[
myPsubs=GetPsubsAll[2,rtry];
tmp2={hraws2lplp,hraws2lmlm,hraws2mpmp,hraws2mmmm,\[Rho]hraws2lpmp,\[Rho]chraws2lpmm,\[Rho]chraws2lmmp,\[Rho]hraws2lmmm,\[CapitalSigma]\[CapitalDelta]hraws2lplm}/.awsubs/.\[CapitalDelta]Ktbl/.myPsubs/.{r->rtry};
]
tmp1-tmp2
*)
(* The difference in evaluation times is extreme: the first takes about 4 seconds, the second takes 217 seconds! *)

t0=Bmat0 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]hat . SinMixMatrix[1,2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,2,lmax,mm]) . CosMixMatrix[2,1,lmax,mm];
t1=Bmat0 . (\[Lambda]hat . SinMixMatrix[1,2,lmax,mm]-a\[Omega]0*SinMixMatrix[0,2,lmax,mm]);
tmpSh1=(t0+t1)/.awsubs;
tmpS\[Kappa]1=Bmat0 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]hat . SinMixMatrix[1,2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,2,lmax,mm])/.awsubs;
t0=Bmat0 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]hat . SinMixMatrix[-1,-2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,-2,lmax,mm]) . CosMixMatrix[-2,1,lmax,mm];
t1=-Bmat0 . (\[Lambda]hat . SinMixMatrix[-1,-2,lmax,mm]-a\[Omega]0*SinMixMatrix[0,-2,lmax,mm]);
tmpSh2=(t0+t1)/.awsubs;
tmpS\[Kappa]2=Bmat0 . (\[CapitalLambda]hat-2*a\[Omega]0*\[Lambda]hat . SinMixMatrix[-1,-2,lmax,mm]+a\[Omega]0^2*SinMixMatrix[0,-2,lmax,mm])/.awsubs;
Evals0Fn[X_,mysubs_,rtry1_]:=X/.\[Gamma]llsubs/.awsubs/.mysubs/.\[CapitalDelta]Ktbl/.{r->rtry1};
cosmix11=CosMixMatrix[1,1,lmax,mm];
cosmixm11=CosMixMatrix[-1,1,lmax,mm];
cosmix12=CosMixMatrix[1,2,lmax,mm];
cosmixm12=CosMixMatrix[-1,2,lmax,mm];
ahvector={hvec,\[Kappa]mat,Dop[r*Dop[hvec]],Ddag[r*Ddag[hvec]],Dop[Dop[\[Kappa]mat]],Ddag[Ddag[\[Kappa]mat]],Map[Dop,hvec],Map[Dop,\[Kappa]mat],Map[Ddag,hvec],Map[Ddag,\[Kappa]mat],Map[Dop,\[CapitalDelta][r]*Map[Ddag,\[Kappa]mat]]+Map[Ddag,\[CapitalDelta][r]*Map[Dop,\[Kappa]mat]],Map[D[#,r]&,\[Kappa]mat]}/.\[Kappa]lsubs/.hlsubs;
(* This vector must match with the one in the function below. *)
Calcs0components[rtry_?NumericQ]:=Module[{\[CapitalDelta]Ktbl,tmp,s0subs,\[CapitalSigma]mat,tmp1,tmp2,tmp3,term1,term2,term3,term4,t0,sinmix,cosmix,MM0,
hvec,\[Kappa]mat,DoprDoph,DdagrDdagh,DopDop\[Kappa],DdagDdag\[Kappa],Dhvec,D\[Kappa]mat,Ddaghvec,Ddag\[Kappa]mat,dd2\[Kappa]mat,d\[Kappa]mat,
hraws0lplp,hraws0lmlm,hraws0mpmp,hraws0mmmm,\[Rho]hraws0lpmp,\[Rho]chraws0lpmm,\[Rho]chraws0lmmp,\[Rho]hraws0lmmm,\[CapitalSigma]\[CapitalDelta]hraws0lplm},
\[CapitalDelta]Ktbl=EvaluateRHS[\[CapitalDelta]Ksubs,rtry]/.awsubs;
tmp=GetSpin0subs[rtry];
s0subs=Join[tmp,\[Kappa]lsubs/.tmp];

{hvec,\[Kappa]mat,DoprDoph,DdagrDdagh,DopDop\[Kappa],DdagDdag\[Kappa],Dhvec,D\[Kappa]mat,Ddaghvec,Ddag\[Kappa]mat,dd2\[Kappa]mat,d\[Kappa]mat}=ahvector/.\[Gamma]llsubs/.awsubs/.s0subs/.\[CapitalDelta]Ktbl/.{r->rtry};

tmp1=-1/(I*\[Omega])*DoprDoph . Bmat0;
tmp2=-4*(ones . DopDop\[Kappa]) . Bmat0;
hraws0lplp=(tmp1+tmp2)/.awsubs;

tmp1=1/(I*\[Omega])*DdagrDdagh . Bmat0;
tmp2=-4*(ones . DdagDdag\[Kappa]) . Bmat0;
hraws0lmlm=(tmp1+tmp2)/.awsubs;

hraws0mpmp=((a/\[Omega])*hvec . tmpSh1-4*(ones . \[Kappa]mat) . tmpS\[Kappa]1)/.awsubs;
hraws0mmmm=(-(a/\[Omega])*hvec . tmpSh2-4*(ones . \[Kappa]mat) . tmpS\[Kappa]2)/.awsubs;

(* l+ m+ *)
t0=(r^2*idmat+a^2*cosmix12)/.{r->rtry}/.awsubs;
sinmix=sinp1;
cosmix=cosmix11;
MM0=(\[Lambda]hat-a*\[Omega]*sinmix)/.awsubs;
term1=1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0/.awsubs;
term2=-a/\[Omega]*(r*Dhvec . Bmat0 . sinmix-hvec . Bmat0 . MM0 . cosmix)/.{r->rtry}/.awsubs;
tmp=(r*idmat+I*a*cosmix)/.{r->rtry}/.awsubs;
term3=4*(ones . D\[Kappa]mat) . Bmat0 . MM0 . tmp;
term4=-4*((ones . \[Kappa]mat) . Bmat0 . MM0+I*a*(ones . D\[Kappa]mat) . Bmat0 . sinmix)/.awsubs;
\[Rho]hraws0lpmp=(term1+term2+term3+term4)/.awsubs;

(* l+ m- *)
t0=(r^2*idmat+a^2*cosmixm12)/.{r->rtry}/.awsubs;
sinmix=sinm1;
cosmix=cosmixm11;
MM0=(\[Lambda]hat-a*\[Omega]*sinmix)/.awsubs;
term1=-1/(2*I*\[Omega])*Dhvec . Bmat0 . MM0 . t0/.awsubs;
term2=a/\[Omega]*(r*Dhvec . Bmat0 . sinmix+hvec . Bmat0 . MM0 . cosmix)/.{r->rtry}/.awsubs;
tmp=(r*idmat-I*a*cosmix)/.{r->rtry}/.awsubs;
term3=-4*((ones . D\[Kappa]mat) . Bmat0 . MM0 . tmp);
term4=4*((ones . \[Kappa]mat) . Bmat0 . MM0+I*a*(ones . D\[Kappa]mat) . Bmat0 . sinmix)/.awsubs;
\[Rho]chraws0lpmm=(term1+term2+term3+term4)/.awsubs;

(* l- m+ *)
t0=(r^2*idmat+a^2*cosmix12)/.{r->rtry}/.awsubs;
sinmix=sinp1;
cosmix=cosmix11;
MM0=(\[Lambda]hat-a*\[Omega]*sinmix)/.awsubs;
term1=-1/(2*I*\[Omega])*Ddaghvec . Bmat0 . MM0 . t0/.awsubs;
term2=-a/\[Omega]*(r*Ddaghvec . Bmat0 . sinmix-hvec . Bmat0 . MM0 . cosmix)/.{r->rtry}/.awsubs;
tmp=(r*idmat-I*a*cosmix)/.{r->rtry}/.awsubs;
term3=4*((ones . Ddag\[Kappa]mat) . Bmat0 . MM0 . tmp);
term4=-4*((ones . \[Kappa]mat) . Bmat0 . MM0-I*a*(ones . Ddag\[Kappa]mat) . Bmat0 . sinmix)/.awsubs;
\[Rho]chraws0lmmp=(term1+term2+term3+term4)/.awsubs;

(* l- m- *)
t0=(r^2*idmat+a^2*cosmixm12)/.{r->rtry}/.awsubs;
sinmix=sinm1;
cosmix=cosmixm11;
MM0=(\[Lambda]hat-a*\[Omega]*sinmix)/.awsubs;
term1=1/(2*I*\[Omega])*Ddaghvec . Bmat0 . MM0 . t0/.awsubs;
term2=a/\[Omega]*(r*Ddaghvec . Bmat0 . sinmix+hvec . Bmat0 . MM0 . cosmix)/.{r->rtry}/.awsubs;
tmp=(r*idmat+I*a*cosmix)/.{r->rtry}/.awsubs;
term3=-4*((ones . Ddag\[Kappa]mat) . Bmat0 . MM0 . tmp);
term4=4*((ones . \[Kappa]mat) . Bmat0 . MM0-I*a*(ones . Ddag\[Kappa]mat) . Bmat0 . sinmix)/.awsubs;
\[Rho]hraws0lmmm=(term1+term2+term3+term4)/.awsubs;

(* l+ l- *)
cosmix=cosmix01;
\[CapitalSigma]mat=(r^2*idmat+a^2*cosmix02)/.{r->rtry}/.awsubs;
tmp1=(dd2\[Kappa]mat)/.\[Kappa]lsubs/.awsubs;
tmp2=d\[Kappa]mat/.awsubs;
tmp3=\[Lambda]hat . (sinm10-sinp10);
term1=-hvec . Bmat0 . ((K[r]/\[Omega])*idmat-2*\[CapitalSigma]mat) . \[CapitalSigma]mat;
term2=-2*(ones . tmp1) . Bmat0 . \[CapitalSigma]mat;
term3=8*r*\[CapitalDelta][r]*(ones . tmp2) . Bmat0;
term4=4*a^2*(ones . \[Kappa]mat) . Bmat0 . tmp3 . cosmix;
\[CapitalSigma]\[CapitalDelta]hraws0lplm=(term1+term2+term3+term4)/.hlsubs/.\[Gamma]llsubs/.\[CapitalDelta]Ktbl/.{r->rtry}/.awsubs;

{hraws0lplp,hraws0lmlm,hraws0mpmp,hraws0mmmm,\[Rho]hraws0lpmp,\[Rho]chraws0lpmm,\[Rho]chraws0lmmp,\[Rho]hraws0lmmm,\[CapitalSigma]\[CapitalDelta]hraws0lplm,hvec}
];
rtry=SetPrecision[7.1,prec];
Timing[
tmp1=Calcs0components[rtry];
];
Timing[
\[CapitalDelta]Ktbl=EvaluateRHS[\[CapitalDelta]Ksubs,rtry]/.awsubs;
myPsubs=GetSpin0subs[rtry];
htraces=Table[If[ll<Abs[mm],0,h[ll][r]],{ll,0,lmax}];
tmp2={hs0lplp,hs0lmlm,hs0mpmp,hs0mmmm,\[Rho]hs0lpmp,\[Rho]chs0lpmm,\[Rho]chraws0lmmp,\[Rho]hraws0lmmm,\[CapitalSigma]\[CapitalDelta]hraws0lplm,htraces}/.\[CapitalDelta]Ktbl/.myPsubs/.{r->rtry};
];
tmp1-tmp2;
t0=Calcs0components[2.25][[10,3]];
t1=Calcs0components[SetPrecision[2.25,prec]][[10,3]];
t2=0.038466744124889374` +0.002573275240138173` I;
t0-t1;
t1-t2;



zeros=Table[SetPrecision[0,prec],{ll,0,lmax}];
(*
comps2={hraws2lplp,hraws2lmlm,hraws2mpmp,hraws2mmmm,\[Rho]hraws2lpmp,\[Rho]chraws2lpmm,\[Rho]chraws2lmmp,\[Rho]hraws2lmmm,\[CapitalSigma]\[CapitalDelta]hraws2lplm,zeros}/.awsubs;
comps1={hraws1lplp,hraws1lmlm,hraws1mpmp,hraws1mmmm,\[Rho]hraws1lpmp,\[Rho]chraws1lpmm,\[Rho]chraws1lmmp,\[Rho]hraws1lmmm,\[CapitalSigma]\[CapitalDelta]hraws1lplm,zeros}/.awsubs;
*)
htraces=Table[If[ll<Abs[mm],0,h[ll][r]],{ll,0,lmax}];
comps0={hs0lplp,hs0lmlm,hs0mpmp,hs0mmmm,\[Rho]hraws0lpmp,\[Rho]chraws0lpmm,\[Rho]chraws0lmmp,\[Rho]hraws0lmmm,\[CapitalSigma]\[CapitalDelta]hraws0lplm,htraces};
(* The replacement of Pp2 with Pm2 leads to an horrific slowing-down of evaluation. It is more efficient not to make replacements of P+2 to P-2, or order-reduce, and instead to replace a larger set of functions. Also, calculating the last 5 components is horribly slow if done naively. I have written functions to speed up this evaluation significantly. *)
GetComp[spin_,rtry_]:=Module[{comp},
If[spin==2,
comp=Join[Calcs2components[rtry],{zeros}];
];
If[spin==1,
comp=Join[Calcs1components[rtry],{zeros}];
];
If[spin==0,
comp=Calcs0components[rtry];
];
Evaluate[comp]
];
GetComponents[rtry_]:=Module[{s2comp,s1comp,s0comp},
s2comp=Join[Calcs2components[rtry],{zeros}];
s1comp=Join[Calcs1components[rtry],{zeros}];
s0comp=Calcs0components[rtry];
Evaluate[s2comp+s1comp+s0comp]
];
rtry=SetPrecision[7.1,prec];
Timing[
tmp2=GetComponents[2,rtry];
];
Timing[
tmp1=GetComponents[1,rtry];
];
Timing[
tmp0=GetComponents[0,rtry];
];
(* Check whether the components appear to be continuous at r=r0. *)
tiny\[Epsilon]=SetPrecision[10^(-prec+3),prec];
Timing[
compIN=GetComponents[r0-tiny\[Epsilon]];
compUP=GetComponents[r0+tiny\[Epsilon]];
compIN-compUP;
];
(* Some degradation at large-l, due to the truncation of the sum, as expected. *)
(* This is a bad way to do it, change in future *)
rsR[[1]]=r0+tiny\[Epsilon];
If[rsR[[1]]>r0,Print["True"];];

Print["s = 2, UP ..."];
s2R=Table[GetComp[2,rsR[[ri]]],{ri,1,Length[rsR]}];
Print["s = 1, UP ..."];
s1R=Table[GetComp[1,rsR[[ri]]],{ri,1,Length[rsR]}];
Print["s = 0, UP ..."];
s0R=Table[GetComp[0,rsR[[ri]]],{ri,1,Length[rsR]}];
Print["s = 2, IN ..."];
s2L=Table[GetComp[2,rsL[[ri]]],{ri,1,Length[rsL]}];
Print["s = 1, IN ..."];
s1L=Table[GetComp[1,rsL[[ri]]],{ri,1,Length[rsL]}];
Print["s = 0, IN ..."];
s0L=Table[GetComp[0,rsL[[ri]]],{ri,1,Length[rsL]}];
sallR=s2R+s1R+s0R;
sallL=s2L+s1L+s0L;
sallR[[1]]-sallL[[1]];

(* Save in a binary format. *)
GetArray[tbl_]:=Module[{grid},
dim=Dimensions[tbl];
arr=Developer`ToPackedArray[Table[SetPrecision[0,prec],{dim[[1]]},{2*dim[[2]]},{dim[[3]]}]];
For[ri=1,ri<=dim[[1]],ri++,
For[qi=1,qi<=10,qi++,
For[ll=0,ll<=lmax,ll++,
arr[[ri,qi,ll+1]]=Re[tbl[[ri,qi,ll+1]]];
arr[[ri,10+qi,ll+1]]=Im[tbl[[ri,qi,ll+1]]];
];
];
];
arr
];
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



End[];
EndPackage[];
