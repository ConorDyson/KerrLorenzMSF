(* ::Package:: *)

Needs["SpinWeightedSpheroidalHarmonics`"];
Needs["Teukolsky`"];


Off[ClebschGordan::phy];
Off[ClebschGordan::tri];
Off[Infinity::indet];
Off[SpinWeightedSphericalHarmonicY::params];
Off[Power::infy];


MetricReconstructStatic[iConfig_,primarypath_]:=Block[{$MaxExtraPrecision = Infinity},


Print["1. Preliminaries"];



extractUpToLastHyphen[str_String] := StringJoin[Riffle[Most[StringSplit[str, "-"]], "-"]];
extractUpToSecondLastHyphen[str_String] := StringJoin[Riffle[Most[StringSplit[extractUpToLastHyphen[str], "-"]], "-"]];
datapath=(primarypath<>"/data"<>extractUpToSecondLastHyphen[ToString[iConfig]]<>"/data"<>extractUpToLastHyphen[ToString[iConfig]]<>"/");
If[!DirectoryQ[datapath], CreateDirectory[datapath, CreateIntermediateDirectories -> True]];

mypath="";
filepath=datapath<>"config/";
filen=filepath<>"config"<>ToString[iConfig]<>".txt";

Print[filen];


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
(* Functions for converting r to z and back again. *)
Getz:=(2*#-rp-rm)/(rp-rm)&;
Getr:=1/2 (rm+rp-rm #+rp #)&;
Getz[Getr[2]]//Simplify;
Getr[Getz[6]]//Simplify;
(* Function for initialising the numerical parameters. *)
(* NOTE: It seems necessary to specific the input parameters (r0 and a) with high precision, e.g., 6.0`60. *)
SetParams[rr0_,aa_]:=Module[{},
subs1={M->1,a->aa,rp->1+Sqrt[1-aa^2],rm->1-Sqrt[1-aa^2]};
subs2={r0->rr0,z0->Getz[rr0]}/.subs1;
Join[subs1,subs2]
];
rtozsubs={r->Getr[z]};
SetParams[6.0`60,0.9`60];
iround=1000;
prec=GetParam["prec",50];(* Number of digits to use. This must be read in first. *)
r0val=Round[GetParam["r0"]*iround]/iround;
aval=Round[GetParam["a"]*iround]/iround;
lmax=GetParam["lmax",10];
lplot=GetParam["lplot"];
mypar=SetParams[SetPrecision[r0val,prec],SetPrecision[aval,prec]];
rhor=1+Sqrt[1-aval^2];
rgrid=GetParam["rgrid",1]; (* Use a linearly-spaced grid in the variable: 0 = rstar , 1 = r. *)
rstarmin=GetParam["rstmin",SetPrecision[rhor+0.1,prec]]; (* If rgrid=1 then these will be interpreted as rmin and rmax, instead of rstarmin and rstarmax. *)
rstarmax=GetParam["rstmax",SetPrecision[rmax,prec]];
nres=GetParam["n",4];  (* Resolution in the r* direction:  dr* = M / n  (or dr = M / n). *)
dformat=ToString@GetParam["dformat","Real64"]; 
\[CapitalDelta]=r^2-2*M*r+a^2;
\[CapitalDelta]alt=(r-rp)*(r-rm);
\[CapitalSigma]=r^2+a^2*Cos[\[Theta]]^2;
\[Lambda]2hat=\[CapitalLambda]h/Sqrt[\[Lambda]];
Lop[n_,S_]:=D[S,\[Theta]]+n*Cot[\[Theta]]*S;
Y0eq=Lop[1,Lop[0,Y0[\[Theta]]]]+\[Lambda]*Y0[\[Theta]];
Y1eq=Lop[0,Lop[1,Y1[\[Theta]]]]+\[Lambda]*Y1[\[Theta]];
Y2eq=Lop[-1,Lop[2,Y2[\[Theta]]]]+(\[Lambda]-2)*Y2[\[Theta]];
tmp=Solve[Y0eq==0,Y0''[\[Theta]]]//First;
Y0subs=Join[tmp,{Y0'''[\[Theta]]-> D[Y0''[\[Theta]]/.tmp,\[Theta]]/.tmp}//Simplify];
tmp=Solve[Y1eq==0,Y1''[\[Theta]]]//First;
Y1subs=Join[tmp,{Y1'''[\[Theta]]-> D[Y1''[\[Theta]]/.tmp,\[Theta]]/.tmp}//Simplify];
tmp=Solve[Y2eq==0,Y2''[\[Theta]]]//First;
Y2subs=Join[tmp,{Y2'''[\[Theta]]-> D[Y2''[\[Theta]]/.tmp,\[Theta]]/.tmp}//Simplify];
Y2alt1=(1/\[Lambda]2hat)*Lop[-1,Y1[\[Theta]]]/.Y1subs;
Y2alt0=(1/\[CapitalLambda]h)*Lop[-1,Lop[0,Y0[\[Theta]]]]/.Y0subs;
Y0alt1=-1/Sqrt[\[Lambda]]*Lop[1,Y1[\[Theta]]];
Simplify[(Lop[-1,Lop[2,Y2alt1]]+(\[Lambda]-2)*Y2alt1)/.Y1subs];
Simplify[(Lop[-1,Lop[2,Y2alt0]]+(\[Lambda]-2)*Y2alt0)/.Y0subs];
Simplify[(Lop[1,Lop[0,Y0alt1]]+\[Lambda]*Y0alt1)/.Y1subs];
Y2toY1subs={Y2[\[Theta]]->Y2alt1,Y2'[\[Theta]]->Simplify[ D[Y2alt1,\[Theta]]/.Y1subs]};
Y2toY0subs={Y2[\[Theta]]->Y2alt0,Y2'[\[Theta]]->Simplify[ D[Y2alt0,\[Theta]]/.Y0subs]};
Y0toY1subs={Y0[\[Theta]]->Y0alt1,Y0'[\[Theta]]->Simplify[ D[Y0alt1,\[Theta]]/.Y1subs]};
r0subs={r->r0};
\[CapitalDelta]0subs={\[CapitalDelta]0->r0^2-2*M*r0+a^2};
\[Lambda]subs={\[Lambda]->l*(l+1)};
\[CapitalLambda]hsubs={\[CapitalLambda]h-> Sqrt[(l-1)*l*(l+1)*(l+2)]};
\[Theta]subs={Cos[\[Theta]]->cosq,Sin[\[Theta]]->sinq,Cos[2\[Theta]]-> 2*cosq^2-1,Sin[2\[Theta]]->2*cosq*sinq,Csc[\[Theta]]->1/sinq,Cot[\[Theta]]-> cosq/sinq};

aMsubs={a->Sqrt[rp*rm],M->(rp+rm)/2};
rprmsubs={rp->M+Sqrt[M^2-a^2],rm-> M-Sqrt[M^2-a^2]};
logrepl={Log[(r-rm)/(r-rp)]-> logm-logp,Log[(r-rp)/(r-rm)]-> logp-logm,Log[r-rp]->logp,Log[r-rm]->logm};
reverselogrepl={logp->Log[r-rp],logm-> Log[r-rm]};
(* Metric and principal null tetrad. *)
lvec={(r^2+a^2)/\[CapitalDelta],1,0,a/\[CapitalDelta]};
nvec=1/(2*\[CapitalSigma])*{r^2+a^2,-\[CapitalDelta],0,a};
mvec=1/(Sqrt[2]*(r+I*a*Cos[\[Theta]]))*{I*a*Sin[\[Theta]],0,1,I/Sin[\[Theta]]};
mbarvec=1/(Sqrt[2]*(r-I*a*Cos[\[Theta]]))*{-I*a*Sin[\[Theta]],0,1,-I/Sin[\[Theta]]};
gup=Table[-lvec[[i]]*nvec[[j]]-nvec[[i]]*lvec[[j]]+mvec[[i]]*mbarvec[[j]]+mbarvec[[i]]*mvec[[j]],{i,1,4},{j,1,4}];
gdn=Simplify[Inverse[gup]];
lp=lvec;
lm=-2*\[CapitalSigma]/\[CapitalDelta]*nvec//Simplify;
mp=mvec*Sqrt[2]*(r+I*a*Cos[\[Theta]]);
mm=mbarvec*Sqrt[2]*(r-I*a*Cos[\[Theta]]);
lpdn=gdn . lp//Simplify;
lmdn=gdn . lm//Simplify;
mpdn=gdn . mp//Simplify;
mmdn=gdn . mm//Simplify;

{lpdn . lm-2*\[CapitalSigma]/\[CapitalDelta],mpdn . mm-2*\[CapitalSigma]}//Simplify;
vecs={lpdn,lmdn,mpdn,mmdn};

tensors=Table[1/2*Table[vecs[[i]][[k]]*vecs[[j]][[l]]+vecs[[i]][[l]]*vecs[[j]][[k]],{k,1,4},{l,1,4}],{i,1,4},{j,1,4}]//Simplify;
(* try a check *)

t0=1/2*Table[(lpdn[[k]]*lmdn[[l]]+lmdn[[k]]*lpdn[[l]]),{k,1,4},{l,1,4}]//Simplify;

t1=tensors[[1,2]];
Simplify[t0-t1];
Y0pi2subs={Y0[\[Pi]/2]->SphericalHarmonicY[l,0,\[Pi]/2,0]};
S0subs={S0->Simplify[Y2alt/.{\[Theta]->\[Pi]/2}]}/.Y0pi2subs;

S0psubs={S0p->Simplify[D[Y2alt,\[Theta]]/.{\[Theta]->\[Pi]/2}]}/.{Y0'[\[Pi]/2]->D[SphericalHarmonicY[l,0,\[Theta],0],\[Theta]]/.{\[Theta]->\[Pi]/2}};
numsubs={\[Mu]->1,M->1,r0->r0val,l->lval,rp->2,rm->0};
D0subs={D0-> -4*\[Pi]*\[Mu]/(r0^2*ut)};
fixfac=-2;

gravfac=-2*fixfac/(\[Lambda]*(\[Lambda]-2));
ABkerr={A-> (EE*(r^2+a^2)-a*LL)/\[CapitalDelta],B->LL-a*EE};
Clear[v];
ELkerr=Simplify[{EE-> (1-2*v^2+atil*v^3)/Sqrt[1-3*v^2+2*atil*v^3],LL-> r*v*(1-2*atil*v^3+atil^2*v^4)/Sqrt[1-3*v^2+2*atil*v^3]}/.{atil->a/M,v-> Sqrt[M/r]},Assumptions-> r0>2*M&&a>0&&a<M&&M>0];
utkerr=Simplify[{ut-> Simplify[(gup . {-EE,0,0,LL})[[1]]/.{\[Theta]->\[Pi]/2}]}/.ELkerr];
chk=Simplify[(ut/.utkerr)-(gup[[1,4]]*LL-gup[[1,1]]*EE)/.{\[Theta]->\[Pi]/2}/.ELkerr];
GetABCD[par_]:=Module[{mysubs1,mysubs2,mysubs3},
mysubs1=ABkerr/.ELkerr/.{r->r0}/.par;
mysubs2=Join[ELkerr,utkerr/.ELkerr]/.{r->r0}/.par;
mysubs3=D0subs/.mysubs2/.{r->r0}/.par/.{\[Mu]->1};
Join[mysubs1,mysubs2,mysubs3]
];

(* Various replacement rules for the second derivatives. *)
eq=D[\[CapitalDelta]*D[hl[r],r],r]-\[Lambda]*hl[r];
tmp=Solve[eq==0,hl''[r]]//First//Simplify;
hsubs=Join[tmp,{hl'''[r]-> Simplify[D[hl''[r]/.tmp,r]/.tmp]}];
eq=\[CapitalDelta]*D[D[DeltaB0[r],r],r]-\[Lambda]*DeltaB0[r]+2*r*\[CapitalDelta]*hl[r];
tmp=Solve[eq==0,DeltaB0''[r]]//First//Simplify;
B0subs=Join[tmp,{DeltaB0'''[r]-> Simplify[D[DeltaB0''[r]/.tmp,r]/.tmp]}];
eq=D[\[CapitalDelta]*D[Z[0][r],r],r]-\[Lambda]*Z[0][r]-hl[r];
tmp=Solve[eq==0,Z[0]''[r]]//First//Simplify;
Z0subs=Join[tmp,{Z[0]'''[r]-> Simplify[D[Z[0]''[r]/.tmp,r]/.tmp]}];
eq=D[\[CapitalDelta]*D[R0[r],r],r]-\[Lambda]*R0[r];
tmp=Solve[eq==0,R0''[r]]//First//Simplify;
R0subs=Join[tmp,{R0'''[r]-> Simplify[D[R0''[r]/.tmp,r]/.tmp]}];
Simplify[D[\[CapitalDelta]*D[R0[r],r],r]-\[Lambda]*R0[r]/.R0subs];
Simplify[D[D[\[CapitalDelta]*D[R0[r],r],r]-\[Lambda]*R0[r],r]/.R0subs];
eq=\[CapitalDelta]*P2''[r]-D[\[CapitalDelta],r]*P2'[r]-(\[Lambda]-2)*P2[r];
P2subs=Solve[eq==0,P2''[r]]//First//Simplify;
Simplify[eq/.P2subs];
(* Connection and functions for calculating the Lorenz gauge condition. *)
sqrtg=Sin[\[Theta]]*\[CapitalSigma];
\[Rho]=r+I*a*Cos[\[Theta]];
\[Rho]c=r-I*a*Cos[\[Theta]];
\[Rho]\[Theta]=D[\[Rho],\[Theta]]/\[Rho];
\[Rho]c\[Theta]=D[\[Rho]c,\[Theta]]/\[Rho]c;
xs={t,r,\[Theta],\[Phi]};
\[CapitalGamma]dn=1/2*Table[D[gdn[[k,i]],xs[[j]]]+D[gdn[[i,j]],xs[[k]]]-D[gdn[[j,k]],xs[[i]]],{i,1,4},{j,1,4},{k,1,4}]//Simplify;
\[CapitalGamma]conn=gup . \[CapitalGamma]dn//Simplify;
(* assume that \[Xi] has components downstairs *)
GaugeVecToMP[\[Xi]_]:=Module[{part1,part2},
part1=Table[D[\[Xi][[j]],xs[[k]]]+D[\[Xi][[k]],xs[[j]]],{j,1,4},{k,1,4}];
part2=-2*\[Xi] . \[CapitalGamma]conn;
part1+part2
];
GetTrace[h_]:=Module[{},
Sum[h[[j,k]]*gup[[j,k]],{j,1,4},{k,1,4}]
];
GetZDiv[hb_]:=Module[{t0,t1,t2},
t0=Table[Sum[D[hb[[i,j]],xs[[k]]]*gup[[j,k]],{j,1,4},{k,1,4}],{i,1,4}];
t1=Table[Sum[\[CapitalGamma]conn[[l,i,k]]*hb[[l,j]]*gup[[j,k]],{j,1,4},{k,1,4},{l,1,4}],{i,1,4}];
t2=Table[Sum[\[CapitalGamma]conn[[l,j,k]]*hb[[i,l]]*gup[[j,k]],{j,1,4},{k,1,4},{l,1,4}],{i,1,4}];
(t0-t1-t2)
];
(* Test the connection has been calculated correctly. *)
(* Metric compatibility condition. *)
dergdn=Table[D[gdn[[i,j]],xs[[k]]],{i,1,4},{j,1,4},{k,1,4}];
connterms=Table[Sum[\[CapitalGamma]conn[[l,i,k]]*gdn[[j,l]]+\[CapitalGamma]conn[[l,j,k]]*gdn[[i,l]],{l,1,4}],{i,1,4},{j,1,4},{k,1,4}];
covdn=dergdn-connterms;
Simplify[covdn];
(* Make a metric perturbation from the components. *)
(* The argument should be a list of five components. *)
MakeMetricPerturbation[comp_]:=Module[{hlplp=comp[[1]],hmpmp=comp[[2]],\[Rho]hlpmp=comp[[3]],\[CapitalSigma]hmpmm=comp[[4]],trace=comp[[5]],htmp,\[Rho]chlpmm},
\[Rho]chlpmm=\[Rho]hlpmp/.z_Complex:>Conjugate[z];
htmp=hlplp*(tensors[[1,1]]+tensors[[2,2]])*(\[CapitalDelta]/(2*\[CapitalSigma]))^2
+hmpmp*(tensors[[3,3]]+tensors[[4,4]])*(1/(2*\[CapitalSigma]))^2
+\[Rho]hlpmp/\[Rho]*(2*tensors[[1,3]]+2*tensors[[2,4]])*(\[CapitalDelta]/(2*\[CapitalSigma])^2)
+\[Rho]chlpmm/\[Rho]c*(2*tensors[[1,4]]+2*tensors[[2,3]])*(\[CapitalDelta]/(2*\[CapitalSigma])^2)
+2*(\[CapitalSigma]hmpmm/\[CapitalSigma])*(tensors[[3,4]])*(1/(2*\[CapitalSigma])^2)
-2*\[CapitalDelta]*(\[CapitalSigma]hmpmm/\[CapitalSigma]-\[CapitalSigma]*trace)*(tensors[[1,2]])*(1/(2*\[CapitalSigma])^2);
htry=1/2*Table[htmp[[j,k]]+htmp[[k,j]],{j,1,4},{k,1,4}];
htry
];
htohbarcomponent[comp_]:=Simplify[{comp[[1]],comp[[2]],comp[[3]],comp[[4]]-\[CapitalSigma]^2*comp[[5]],-comp[[5]]},Assumptions->assum];
Makehbar[comp_]:=MakeMetricPerturbation[htohbarcomponent[comp]];
(* Test for consistency. *)
components={H1,H2,H3a+I*H3b,H4,H5};
htmp=MakeMetricPerturbation[components];
Simplify[htmp . lp . lp-H1];
Simplify[htmp . lm . lm-H1];
Simplify[htmp . mp . mp-H2];
Simplify[htmp . mm . mm-H2];
Simplify[\[Rho]*htmp . lp . mp-(H3a+I*H3b)];
Simplify[\[Rho]*htmp . lm . mm-(H3a+I*H3b)];
Simplify[\[Rho]c*htmp . lp . mm-(H3a-I*H3b)];
Simplify[\[Rho]c*htmp . lm . mp-(H3a-I*H3b)];
Simplify[\[CapitalSigma]*htmp . mp . mm-H4];
Simplify[GetTrace[htmp]-H5];
(* Now set up the vat (a matrix), to be filled later *)
GetEmptyVat[lm_]:=Table[0,{qi,1,5},{ll,0,lm}];
nmax=3;
minl={0,2,1,0,0};
spins={0,-2,-1,0,0};
Vat=GetEmptyVat[lmax];
(* vat1 should be smaller than, or the same size as, vat2. *)
AddVat[vat1_, vat2_]:=Module[{},
l1=Length[vat1[[1]]]-1;
l2=Length[vat2[[1]]]-1;
newvat=vat2;
For[qi=1,qi<=5,qi++,
For[ll=0,ll<=l1,ll++,
newvat[[qi]][[ll+1]]+=vat1[[qi]][[ll+1]];
];
];
newvat
];
(* Restore the coupling coefficients. Note that:
Cc[n][l][k] is zero if |k|>n or k+n odd
 and Cs[n][l][k] is zero if |k|>n+1 or k+n even.
*)
GetCsubs[l_,nmax_]:=Module[{},
mysubs={C0c[0]->KroneckerDelta[k,0],C1c[0]->KroneckerDelta[k,0]};
For[nn=1,nn<=nmax,nn++,
tmp=Table[If[Abs[kk]<=nn&&Mod[nn+kk,2]==0,C0c[nn][l][kk]*KroneckerDelta[k,kk],0],{kk,-nn,nn}];
AppendTo[mysubs,C0c[nn]->Sum[tmp[[kk]],{kk,1,Length[tmp]}]];
tmp=Table[If[Abs[kk]<=nn&&Mod[nn+kk,2]==0,C1c[nn][l][kk]*KroneckerDelta[k,kk],0],{kk,-nn,nn}];
AppendTo[mysubs,C1c[nn]->Sum[tmp[[kk]],{kk,1,Length[tmp]}]];
];
mysubs2={};
For[nn=0,nn<=nmax,nn++,
tmp=Table[If[Abs[kk]<=nn+1&&Mod[nn+kk,2]==1,C0s[nn][l][kk]*KroneckerDelta[k,kk],0],{kk,-nn-1,nn+1}];
AppendTo[mysubs2,C0s[nn]->Sum[tmp[[kk]],{kk,1,Length[tmp]}]];
tmp=Table[If[Abs[kk]<=nn+1&&Mod[nn+kk,2]==1,C1s[nn][l][kk]*KroneckerDelta[k,kk],0],{kk,-nn-1,nn+1}];
AppendTo[mysubs2,C1s[nn]->Sum[tmp[[kk]],{kk,1,Length[tmp]}]];
];
Join[mysubs,mysubs2,{C2c[0]-> KroneckerDelta[k,0]}]
];
(* For numerical checks *)
GetRandomIC[fns_]:=Module[{j},
Table[{fns[[j]]'[r]->RandomReal[],fns[[j]][r]->RandomReal[]},{j,1,Length[fns]}]//Flatten
];
(* Load the coefficients that were generated in CalcStaticTensorMode.nb *)
cosubs=Get[mypath<>"static/ProjCo.m"];
tracebccoeffs=Get[mypath<>"static/TraceCo.m"];
SpinWeightedY[ss_,ll_,\[Theta]_]:=If[ll>=Abs[ss],SpinWeightedSphericalHarmonicY[ss,ll,0,\[Theta],0],0];
BrewVat[vat_,lmax_]:=Module[{},
lm=Min[lmax,Length[vat[[1]]]-1];
(* Print[lm]; *)
tmp=Table[Sum[vat[[qi]][[ll+1]]*SpinWeightedY[spins[[qi]],ll,\[Theta]],{ll,0,lm}],{qi,1,5}];
tmp/.cosubs/.tracebccoeffs
];
(* Visualise the mode couplings *)
Conj[x_]:=x/.z_Complex:>Conjugate[z];
visfn=If[ToString[#]== "0",".",If[ToString[Simplify[Conj[#]-#]]=="0","1",If[ToString[Simplify[Conj[#]+#]]== "0", "i","c"]]]&;
visfn[0];

Print["2. Read metric components"];
Csubsfn[dl_]:=Join[Flatten[Table[{C0c[nn]->C0c[dl][nn],C0s[nn]-> C0s[dl][nn],C1c[nn]->C1c[dl][nn], C1s[nn]->C1s[dl][nn],C2c[nn]->C2c[dl][nn]},{nn,0,nmax}]],{Y2-> C2c[dl][0]}];
scalarreplacefn[X_]:={R0''[r]-> X''[r],R0'[r]->X'[r],R0[r]->X[r]};
filen=mypath<>"static/m0components_free_scalar.m";
hsccomponents=Get[filen]/.{Y2->C2c[0]};
hsccomp=hsccomponents;
(* Should we consider the l=0 case separately? The formula above breaks down because of \[Lambda] in denominator. *)
ser=Series[hsccomp,{\[Lambda],0,3}] ;(* The only problematic part is in the third component. *)
myco=Simplify[SpinWeightedSphericalHarmonicY[-1,1,0,\[Theta],0]/Sin[\[Theta]]];
component3=2*I*a/Sqrt[4*\[Pi]]*1/myco*R0'[r]*KroneckerDelta[k,1];
hsc0comp={hsccomp[[1]],hsccomp[[2]],component3,hsccomp[[4]],0};
(* Note that here we have the coupling coefficients C0c[0], C0c[2],C0s[1], C1c[0], C1c[1], and C1s[0]. *)
PourJug[lm_,l1_,jug_,nm_]:=Module[{},
newvat=GetEmptyVat[lm];
For[kk=-nm-1,kk<=nm+1,kk++,
For[qi=1,qi<=5,qi++,
ind=l1+kk+1;
If[(l1+kk>=minl[[qi]])&&(l1+kk<=lm),
tmp=Simplify[jug[[qi]]/.{k->kk}/.{l->l1}];
newvat=ReplacePart[newvat,{qi,ind}->newvat[[qi]][[ind]]+ tmp]
];
];
];
newvat

];


PourFreeScalar[lm_,l_,Aco_,fn_]:=Module[{nm=2,Rrepl},
If[l==0,tmp=hsc0comp/.GetCsubs[l,nm];,tmp=hsccomp/.GetCsubs[l,nm];];
Rrepl=scalarreplacefn[fn];
jug=Aco*tmp/.{\[Lambda]->l*(l+1)}/.{\[CapitalLambda]h->Sqrt[(l-1)*l*(l+1)*(l+2)]}/.Rrepl;

PourJug[lm,l,jug,nm]
];
PourFreeScalar[lm_,l_,Aco_]:=PourFreeScalar[lm,l,Aco,R0[l]];
ltry=4;
PourFreeScalar[lmax,ltry,1];
t0=BrewVat[PourFreeScalar[lmax,ltry,1],10];
(* Only when we call "BrewVat" do the spin-weighted spherical harmonics take their explicit form. *)
(* Now construct the metric perturbation and verify that it is in Lorenz gauge *)
checks[]:=Module[{},
R0eq=(D[\[CapitalDelta]*R0[l]'[r],r]-l*(l+1)*R0[l][r])/.{l->ltry};
tmp=Solve[R0eq==0,R0[ltry]''[r]]//First;
R0subs=Join[{R0[ltry]'''[r]->D[R0[ltry]''[r]/.tmp,r]/.tmp},tmp]//Simplify;
t0=BrewVat[PourFreeScalar[lmax,ltry,1],lmax];
t1=Makehbar[t0];
(* Schwarzschild case *)
t0=Simplify[GetZDiv[t1/.{a->0}]/.R0subs/.{a->0}];
(* Just do a numerical check in the Kerr case *)
numvals={\[Theta]->1.3,a->0.9,M->1, r->6.7};
{t0,N[GetZDiv[t1]/.R0subs/.GetRandomIC[{R0[ltry]}]/.numvals,30]}
];

If[runchecks,checks[]];
nm=4;
Crepl=GetCsubs[l,nm];
tmp=Get[mypath<>"static/m0components_even_noscalar.m"];
hecomp=Simplify[tmp/.Crepl,Assumptions->assum];
GetEvenReplRules[l1_]:=Module[{P2subs,\[Chi]esubs,\[Chi]ek0subs,P3subs,P4subs},
eq=\[CapitalDelta]*P2[l1]''[r]-D[\[CapitalDelta],r]*P2[l1]'[r]-(\[Lambda]-2)*P2[l1][r];
P2subs=Solve[eq==0,P2[l1]''[r]]//First//Simplify;
eq=D[\[CapitalDelta]*\[Chi]e[l1]'[r],r]-\[Lambda]*\[Chi]e[l1][r]-(r*P2[l1]'[r]-2*P2[l1][r]);
\[Chi]esubs=Solve[eq==0,\[Chi]e[l1]''[r]]//First//Simplify;
eq=(D[\[CapitalDelta]*\[Chi]ek0[l1]'[r],r]-\[Lambda]*\[Chi]ek0[l1][r]-P2[l1]''[r])/.P2subs;
\[Chi]ek0subs=Solve[eq==0,\[Chi]ek0[l1]''[r]]//First//Simplify;
P3subs=D[P2subs,r]/.P2subs//Simplify;
P4subs=D[P3subs,r]/.P2subs//Simplify;
Join[P4subs,P3subs,P2subs,\[Chi]esubs,\[Chi]ek0subs]/.\[Lambda]subs/.{l->l1}
];

GetEvenReplRules[3];
PourEvenParity[lm_,l1_,Aco_]:=Module[{nm=4,P2repl,fnrepl,vat0,vat1,vat2,vat3,vat4,newvat},
P2repl={P2'[r]->P2[l]'[r],P2[r]->P2[l][r]};
jug=hecomp/.P2repl/.\[Lambda]subs/.\[CapitalLambda]hsubs/.{l->l1};
(* Now include the scalar part that restores Lorenz gauge. *)
lamsubs={lam[-1]-> (l-2)*(l-1),lam[0]-> l*(l+1),lam[1]->(l+2)*(l+3),\[Lambda]->l*(l+1),\[CapitalLambda]h->Sqrt[(l-1)*l*(l+1)*(l+2)]};
scm1fac=-a^2/(\[Lambda]*\[CapitalLambda]h)*((\[Lambda]+2)*C0s[1][l][-2]-2*\[Lambda]*C0c[2][l][-2])/(\[Lambda]-lam[-1])/.lamsubs/.{l->l1}; (* P2''[r] *)
scp1fac=-a^2/(\[Lambda]*\[CapitalLambda]h)*((\[Lambda]+2)*C0s[1][l][2]-2*\[Lambda]*C0c[2][l][2])/(\[Lambda]-lam[1])/.lamsubs/.{l->l1}; (* P2''[r] *)
sc0fac=-a^2/(\[Lambda]*\[CapitalLambda]h)*((\[Lambda]+2)*C0s[1][l][0]-2*\[Lambda]*C0c[2][l][0]+2*\[Lambda])/.lamsubs/.{l->l1};
\[Chi]efac=-(\[CapitalLambda]h/\[Lambda])/.lamsubs/.{l->l1};
vat0=PourJug[lm,l1,jug,nm];
vat1=PourFreeScalar[lm,l1-2,scm1fac,P2[l1]''];
vat2=PourFreeScalar[lm,l1+2,scp1fac,P2[l1]''];
vat3=PourFreeScalar[lm,l1,sc0fac,\[Chi]ek0[l1]];
vat4=PourFreeScalar[lm,l1,\[Chi]efac,\[Chi]e[l1]];
newvat=vat0+vat1+vat2+vat3+vat4;
Simplify[newvat/.GetEvenReplRules[l1]/.{\[Lambda]->l1*(l1+1)}]
];
checks[]:=Module[{},
ltry=2;
t0=BrewVat[PourEvenParity[lmax,ltry,1],lmax]/.GetEvenReplRules[ltry];
t1=Makehbar[t0];
(* Schwarzschild case *)
t0=Simplify[GetZDiv[t1/.{a->0}]/.GetEvenReplRules[ltry]/.\[Lambda]subs/.{l->ltry}/.{a->0}];
(* Just do a numerical check in the Kerr case *)
fnvals=GetRandomIC[{\[Chi]e[ltry],\[Chi]ek0[ltry],P2[ltry]}];
numvals={\[Theta]->1.3,a->0.9,M->1, r->6.7};
{t0,N[GetZDiv[t1]/.GetEvenReplRules[ltry]/.\[Lambda]subs/.{l->ltry}/.fnvals/.numvals,30]}
];




If[runchecks,checks[]];
nm=4;
Crepl=GetCsubs[l,nm];
tmp=Get[mypath<>"static/m0components_odd_noscalar.m"];
hocomp=Simplify[tmp/.Crepl,Assumptions->assum];
GetOddReplRules[l1_]:=Module[{P2subs,\[Chi]ocsubs,\[Chi]ossubs},
eq=\[CapitalDelta]*P2[l1]''[r]-D[\[CapitalDelta],r]*P2[l1]'[r]-(\[Lambda]-2)*P2[l1][r];
P2subs=Solve[eq==0,P2[l1]''[r]]//First//Simplify;
eqs=Table[D[\[CapitalDelta]*\[Chi]oc[k][l1]'[r],r]-lamo[k]*\[Chi]oc[k][l1][r]-P2[l1]'[r]==0,{k,-1,1,2}];
\[Chi]ocsubs=Solve[eqs,Table[\[Chi]oc[k][l1]''[r],{k,-1,1,2}]]/.P2subs//First//Simplify;
eqs=Table[D[\[CapitalDelta]*\[Chi]os[k][l1]'[r],r]-lamo[k]*\[Chi]os[k][l1][r]-r*P2[l1]''[r]==0,{k,-1,1,2}]/.P2subs;
\[Chi]ossubs=Solve[eqs,Table[\[Chi]os[k][l1]''[r],{k,-1,1,2}]]/.P2subs//First//Simplify;
lamodd={lamo[-1]-> (l-1)*l, lamo[1]-> (l+1)*(l+2),\[Lambda]->l*(l+1),\[CapitalLambda]h->Sqrt[(l-1)*l*(l+1)*(l+2)]};
Flatten[Join[P2subs,\[Chi]ocsubs,\[Chi]ossubs]]/.lamodd/.{l->l1}
];
GetOddReplRules[3];
PourOddParity[lm_,l1_,Aco_]:=Module[{nm=4,P2repl,P2subs,\[Chi]ocsubs,\[Chi]ossubs,vat0,vat1,vat2,vat3,vat4},
P2repl={P2'[r]->P2[l]'[r],P2[r]->P2[l][r]};
jug=hocomp/.P2repl/.\[Lambda]subs/.\[CapitalLambda]hsubs/.{l->l1};
(* Now include the scalar part that restores Lorenz gauge. *)
lamodd={lamo[-1]-> (l-1)*l, lamo[1]-> (l+1)*(l+2),\[Lambda]->l*(l+1),\[CapitalLambda]h->Sqrt[(l-1)*l*(l+1)*(l+2)]};
scocm1fac=(a*\[CapitalLambda]h/\[Lambda])*C0c[1][l][-1]/.lamodd/.{l->l1}; 
scocp1fac=(a*\[CapitalLambda]h/\[Lambda])*C0c[1][l][1]/.lamodd/.{l->l1}; 
scosm1fac=-(a*\[CapitalLambda]h/\[Lambda]^2)*C0s[0][l][-1]/.lamodd/.{l->l1}; 
scosp1fac=-(a*\[CapitalLambda]h/\[Lambda]^2)*C0s[0][l][1]/.lamodd/.{l->l1}; 
vat0=PourJug[lm,l1,jug,nm];
vat1=PourFreeScalar[lm,l1-1,scocm1fac,\[Chi]oc[-1][l1]];
vat2=PourFreeScalar[lm,l1+1,scocp1fac,\[Chi]oc[1][l1]];
vat3=PourFreeScalar[lm,l1-1,scosm1fac,\[Chi]os[-1][l1]];
vat4=PourFreeScalar[lm,l1+1,scosp1fac,\[Chi]os[1][l1]];
(* Apply replacement rules for all the ingredients above *)
Simplify[(vat0+vat1+vat2+vat3+vat4)/.GetOddReplRules[l1]/.lamodd/.{\[Lambda]->l1*(l1+1)}]
];


checks[]:=Module[{},
ltry=3;
t0=BrewVat[PourOddParity[lmax,ltry,1],10]/.GetOddReplRules[ltry];
t1=Makehbar[t0];(* Schwarzschild case *)
t2=Simplify[GetZDiv[t1/.{a->0}]/.\[Lambda]subs/.{l->ltry}/.{a->0}];
(* Numerical check of the Kerr case *)
fnvals = GetRandomIC[{\[Chi]oc[-1][ltry],\[Chi]oc[1][ltry], \[Chi]os[-1][ltry],\[Chi]os[1][ltry],P2[ltry]}];
numvals={\[Theta]->1.3,a->0.9,M->1, r->6.7};
{t2,N[GetZDiv[t1]/.GetOddReplRules[ltry]/.\[Lambda]subs/.{l->ltry}/.fnvals/.numvals,30]}
];
If[runchecks,checks[]];
filen=mypath<>"static/m0components_divB.m";
hdivBcomp=Get[filen];
GetTraceReplRules[l1_]:=Module[{X0subs,Z0subs,DeltaB0subs,hlsubs},
eq=D[\[CapitalDelta]*X0[l1]'[r],r]-\[Lambda]*X0[l1][r]-r^2*hl[l1][r];
tmp=Solve[eq==0,X0[l1]''[r]]//First//Simplify;
X0subs=Join[tmp,D[tmp,r]/.tmp]//Simplify;
eq=D[\[CapitalDelta]*Z0[l1]'[r],r]-\[Lambda]*Z0[l1][r]-hl[l1][r];
tmp=Solve[eq==0,Z0[l1]''[r]]//First//Simplify;
Z0subs=Join[tmp,D[tmp,r]/.tmp]//Simplify;
eq=\[CapitalDelta]*DeltaB0[l1]''[r]-\[Lambda]*DeltaB0[l1][r]+2*r*\[CapitalDelta]*hl[l1][r];
tmp=Solve[eq==0,DeltaB0[l1]''[r]]//First//Simplify;
DeltaB0subs=Join[tmp,D[tmp,r]/.tmp]//Simplify;
eq=D[\[CapitalDelta]*hl[l1]'[r],r]-\[Lambda]*hl[l1][r];
tmp=Solve[eq==0,hl[l1]''[r]]//First//Simplify;
hlsubs=Join[tmp,D[tmp,r]/.tmp]//Simplify;
Join[X0subs,Z0subs,DeltaB0subs,hlsubs]/.{\[Lambda]->l1*(l1+1)}
];
GetTraceReplRules[4];
PourDiv[lm_,l_,Aco_,fn_]:=Module[{nm=2,Rrepl},
tmp=hdivBcomp/.GetCsubs[l,nm];
Brepl={B[r]->fn,B'[r]->D[fn,r],B''[r]->D[fn,{r,2}]};
jug=Aco*tmp/.Brepl/.{\[Lambda]->l*(l+1)}/.{\[CapitalLambda]h->Sqrt[(l-1)*l*(l+1)*(l+2)]};
PourJug[lm,l,jug,nm]
];
PourTraceMode[lm_,l1_,Aco_]:=Module[{nm=4,lamsubs,vat0,vat1,vat2,vat3,vat4,vat5,vat6,vat7,vat8,newvat},
lamsubs={\[Lambda]->l1*(l1+1),lam[0]->l1*(l1+1),lam[-1]->(l1-2)*(l1-1),lam[1]->(l1+2)*(l1+3)};
(* Scalar parts *)
vat0=PourFreeScalar[lm,l1,1/2,X0[l1]];
vat1=PourFreeScalar[lm,l1,(a^2/2)*C0c[2][l1][0],Z0[l1]];
vat2=PourFreeScalar[lm,l1-2,(a^2/2)*C0c[2][l1][-2]/(\[Lambda]-lam[-1])/.lamsubs,hl[l1]];
vat3=PourFreeScalar[lm,l1+2,(a^2/2)*C0c[2][l1][2]/(\[Lambda]-lam[1])/.lamsubs,hl[l1]];
(* Div parts *)




fac=1/(2*l1*(l1+1));
vat4=PourDiv[lm,l1,fac,DeltaB0[l1][r]];
vat5=PourDiv[lm,l1,fac,\[CapitalDelta]*X0[l1]'[r]];
fac=a^2/(2*l1*(l1+1));
vat6=PourDiv[lm,l1,fac*bco[l1][0],\[CapitalDelta]*Z0[l1]'[r]];
vat7=PourDiv[lm,l1-2,fac*bco[l1][-1]/(\[Lambda]-lam[-1])/.lamsubs,\[CapitalDelta]*hl[l1]'[r]];
vat8=PourDiv[lm,l1+2,fac*bco[l1][1]/(\[Lambda]-lam[1])/.lamsubs,\[CapitalDelta]*hl[l1]'[r]];
(* Don't forget to update the trace itself *)
newvat=vat0+vat1+vat2+vat3+vat4+vat5+vat6+vat7+vat8;
If[l1<=lm,
newvat[[5]][[l1+1]]+=hl[l1][r];
];
(* Apply replacement rules for all the ingredients above *)
Simplify[newvat/.GetTraceReplRules[l1]/.lamsubs/.{\[Lambda]->l1*(l1+1)}/.{l->l1}]
];

checks[]:=Module[{},
ltry=4;
t0=BrewVat[PourTraceMode[lmax,ltry,1],10]/.GetTraceReplRules[ltry];
t1=Makehbar[t0];
(* Schwarzschild case *)
t2=Simplify[GetZDiv[t1/.{a->0}]/.GetTraceReplRules[ltry]/.\[Lambda]subs/.{l->ltry}/.{a->0}];
(* Numerical check of the Kerr case *)
fnvals = GetRandomIC[{Z0[ltry],X0[ltry],hl[ltry],DeltaB0[ltry]}];
numvals={\[Theta]->1.3,a->0.9,M->1, r->6.7};
{t2,N[GetZDiv[t1]/.GetTraceReplRules[ltry]/.\[Lambda]subs/.{l->ltry}/.fnvals/.numvals,30]}
];
If[runchecks,checks[]];
(* These modes were previously constructed and saved in Jigsaw.nb *)
path="completion/";
hmodeB=Get[path<>"hBmode.m"];
hmodeC=Get[path<>"hCmode.m"];
hmodeD=Get[path<>"hDmode.m"];
hmodeE=Get[path<>"hEmode.m"];
hmodeF=Get[path<>"hFmode.m"];
hmodeAM=Get[path<>"h_angmom_mode.m"];
AppendTrace[mode_,tr_,linsert_]:=Module[{lmax,ll,newmode},
lmax=Length[mode[[1]]]-1;
newmode=Append[mode,Table[0,{ll,0,lmax}]];
newmode[[5]][[linsert+1]]=tr;
newmode
];
hcompB=AppendTrace[Transpose[hmodeB],6*Sqrt[4*\[Pi]],0];
hcompC=AppendTrace[Transpose[hmodeC],0,0];
hcompD=AppendTrace[Transpose[hmodeD],(2+2*(rp+rm)/(rp-rm)*Log[(r-rp)/(r-rm)])*Sqrt[4*\[Pi]],0];
hcompE=AppendTrace[Transpose[hmodeE],(-4/(rp-rm)*Log[(r-rp)/(r-rm)])*Sqrt[4*\[Pi]],0];
hcompF=AppendTrace[Transpose[hmodeF],4*a*Sqrt[4*\[Pi]],0];
hcompAM=AppendTrace[Transpose[hmodeAM],0,0];
(* Check that these modes are in Lorenz gauge *)
checks[]:=Module[{},
t0=BrewVat[hcompB,Length[hcompB[[1]]]-1]//Simplify;
t1=Makehbar[t0];
t2a=Simplify[GetZDiv[t1],Assumptions->assum];

t0=BrewVat[hcompC,Length[hcompC[[1]]]-1]//Simplify;
t1=Makehbar[t0];
t2b=Simplify[GetZDiv[t1],Assumptions->assum];

t0=BrewVat[hcompAM,Length[hcompAM[[1]]]-1]//Simplify;
t1=Makehbar[t0];
t2c=Simplify[GetZDiv[t1],Assumptions->assum];
{t2a,t2b,t2c}
];
If[runchecks,checks[]];
(* Numerical checks for the rest. *)
checks[]:=Module[{},
numsubs={rp->1.9,rm->0.3,r->10.0,\[Theta]->0.4};

tmp=Simplify[hcompD/.reverselogrepl/.aMsubs,Assumptions->assum2];
t0=Simplify[BrewVat[tmp,Length[hcompD[[1]]]-1],Assumptions->assum2];
t1=Makehbar[t0]/.aMsubs;
t2=GetZDiv[t1]/.aMsubs;
t3a=N[t2/.numsubs];

tmp=Simplify[hcompE/.reverselogrepl/.aMsubs,Assumptions->assum2];
t0=Simplify[BrewVat[tmp,Length[hcompE[[1]]]-1],Assumptions->assum2];
t1=Makehbar[t0]/.aMsubs;
t2=GetZDiv[t1]/.aMsubs;
t3b=N[t2/.numsubs];

tmp=Simplify[hcompF/.reverselogrepl/.aMsubs,Assumptions->assum2];
t0=Simplify[BrewVat[tmp,Length[hcompF[[1]]]-1],Assumptions->assum2];
t1=Makehbar[t0]/.aMsubs;
t2=GetZDiv[t1]/.aMsubs;
t3c=N[t2/.numsubs];
{t3a,t3b,t3c}
];
If[runchecks,checks[]];

Print["3. Fill Matrix"];
(* Here I will add all the pieces together, but with coefficients that remain to be determined in the "Jumps" section. *)
Vat=GetEmptyVat[lmax+6];
Dimensions[Vat];
Vat=AddVat[coB*hcompB,Vat];
Vat=AddVat[coC*hcompC,Vat];
Vat=AddVat[coD*hcompD,Vat];
Vat=AddVat[coE*hcompE,Vat];
Vat=AddVat[coF*hcompF,Vat];
Vat=AddVat[coAM*hcompAM,Vat];
Dimensions[Vat];
For[ll=2,ll<=lmax+6,ll++,
If[Mod[ll,2]==0,
R0eq=(D[\[CapitalDelta]*R0[ll]'[r],r]-ll*(ll+1)*R0[ll][r]);
R0subs=Solve[R0eq==0,R0[ll]''[r]]//First;
Vat=Vat+Ascal[ll]*PourFreeScalar[lmax+6,ll,1]/.R0subs;
Vat=Vat+Atr[ll]*PourTraceMode[lmax+6,ll,1];
Vat=Vat+Aeven[ll]*PourEvenParity[lmax+6,ll,1];
,
Vat=Vat+Aodd[ll]*PourOddParity[lmax+6,ll,1];
];
];
Vat=Vat/.cosubs/.tracebccoeffs;
(* Jumps. *)
(* The only functions with jumps are P2, P2', hl' (not hl) and R0', R0. All other functions are at least C1 at the particle radius. So make replacements here to get the jumps. *)
haystack={"Z0","X0","DeltaB0","\[Chi]ek0","\[Chi]e","\[Chi]oc[-1]","\[Chi]oc[1]","\[Chi]os[-1]","\[Chi]os[1]"};
mysubs1=Flatten[Table[ToExpression[haystack[[jj]]][ll][r]->0,{jj,1,Length[haystack]},{ll,2,lmax+6}]];
mysubs2=Flatten[Table[ToExpression[haystack[[jj]]][ll]'[r]->0,{jj,1,Length[haystack]},{ll,2,lmax+6}]];
VatJumps=Vat/.Table[hl[ll][r]->0,{ll,2,lmax+6}]/.mysubs1/.mysubs2;
(* Takes a while to run *)
(*
Map[visfn,VatJumps,{2}]//TableForm 
*)


Print["4. Jumps"];

(* work out the jumps in the P2 function, in the trace function, and in the monopole. *)
C2=B^2*S0;
C1=2*B*(I*A*S0p+B*S0/r);
C0=2*A*(I*B*2/r+A*I*a/r)*S0p+A^2*(\[Lambda]-2)*S0;
V0minus4=-(\[Lambda]-2);
J0=\[CapitalDelta]*D0*(C1-2*(r-M)/\[CapitalDelta]*C2);
J1=\[CapitalDelta]*D0*(C0-V0minus4*C2/\[CapitalDelta]);
fac=2*I*\[CapitalDelta]*D0*S0p/r;
J0odd=fac*r*A*B;
J1odd=fac*(2*A*B+a*A^2);
fac=D0*S0;
J0even=fac*(-2*(M-a^2/r)*B^2);
J1even=fac*(\[Lambda]-2)*(\[CapitalDelta]*A^2+B^2);
Simplify[(J0-(J0even+J0odd))];
Simplify[J1-(J1even+J1odd)];
(* These are the jumps in the Weyl scalars for each mode. *)
(* don't forget we will still need to multiply by 'gravfac', relating the jump in Weyl scalars to the jump in the metric perturbation. 
Let's do that here. *)
Jumps2={J0o->-I*gravfac* J0odd,J0e->gravfac*J0even,J1o->-I*gravfac*J1odd,J1e->gravfac*J1even}/.{r->r0};
J1h=4*D0*r^2/\[CapitalDelta]*Y0/.{r->r0};
Jumpmono={coE->EE,coAM->(LL-a*EE)/M};

SetNumericalJumps[mypar_]:=Module[{},
ABCDsubs=GetABCD[mypar];
P2evensubs=Table[
{P2[ll][r]-> J0e,P2[ll]'[r]-> J1e}/.Jumps2/.\[Lambda]subs/.{l->ll}/.{S0->SpinWeightedSphericalHarmonicY[-2,ll,0,\[Pi]/2,0]}/.ABCDsubs/.mypar,{ll,2,lmax+6,2}]//Flatten;
P2oddsubs=Table[
{P2[ll][r]-> J0o,P2[ll]'[r]-> J1o}/.Jumps2/.\[Lambda]subs/.{l->ll}/.{S0p->D[SpinWeightedSphericalHarmonicY[-2,ll,0,\[Theta],0],\[Theta]]/.{\[Theta]->\[Pi]/2}}/.ABCDsubs/.mypar,{ll,3,lmax+5,2}]//Flatten;
hl1subs=Table[
(hl[ll]'[r]-> J1h)/.\[Lambda]subs/.{l->ll}/.{Y0->SphericalHarmonicY[ll,0,\[Pi]/2,0]}/.ABCDsubs/.mypar,{ll,2,lmax+6,2}];
Acosubs=Table[{Atr[ll]->1,Aeven[ll]->1,Aodd[ll]->1,R0[ll][r]->Jsc0[ll]/Ascal[ll],R0[ll]'[r]->Jsc1[ll]/Ascal[ll]},{ll,2,lmax+6}]//Flatten;
VatJumps/.P2evensubs/.P2oddsubs/.hl1subs/.Jumpmono/.ABCDsubs/.Acosubs/.reverselogrepl/.{r->r0}/.{z->z0}/.mypar//Simplify
];

SolveJumpEquations[V_]:=Module[{nsol1,checks1},
unknowns0={coB,coC,coD,coF,Jsc0[2],Jsc1[2]};
eqs0={V[[1,1]]==0,V[[4,1]]==0,V[[3,2]]==0,V[[1,3]]==0,V[[2,3]]==0,V[[3,3]]==0};
If[Abs[a/.mypar]<=0.001,
(* In the Schwarzschild case, use the trace equation. *)
eqs0={V[[1,1]]==0,V[[4,1]]==0,V[[3,2]]==0,V[[1,3]]==0,V[[2,3]]==0,V[[5,1]]==0};
];
nsol0=NSolve[eqs0,unknowns0, WorkingPrecision->prec]//First;
(* Now take a larger set of equations, and check a larger set of consistency conditions *)
unknowns1=Join[unknowns0,Table[{Jsc0[ll],Jsc1[ll]},{ll,4,lmax+6,2}]//Flatten,{Jsc0[lmax+4]}];
eqs1=Join[eqs0,Table[{V[[1,ll+1]]==0,V[[2,ll+1]]},{ll,4,lmax+6,2}]//Flatten,{V[[4,lmax+2+1]]==0}];
nsol1=NSolve[eqs1,unknowns1, WorkingPrecision->prec]//First;
checks1=Join[{V[[5,1]],V[[4,3]]},Table[{V[[1,ll+1]],V[[2,ll+1]],V[[3,ll]],V[[3,ll+1]],V[[4,ll+1]]},{ll,4,lmax,2}]//Flatten]/.nsol1;
{nsol1,checks1}
];


VatJumpsNum=SetNumericalJumps[mypar];
{nsol1,checks1}=SolveJumpEquations[VatJumpsNum];

VatJumpsNum=SetNumericalJumps[mypar];
JumpsScalar=SolveJumpEquations[VatJumpsNum][[1]];


(* Make a table of the derived unknowns for inclusion in a tex file. *)
header={"l","Jtensor0","Jtensor1","h1","Jsc0","Jsc1"};
nsf=8;
zerosubs=Flatten@Table[{Jsc0[ll]->0,Jsc1[ll]->0,hl[ll]'[r]->0},{ll,1,lmax,2}];
mytable=Table[Join[{ll},Map[NumberForm[#,nsf]&,{P2[ll][r],P2[ll]'[r],hl[ll]'[r],Jsc0[ll],Jsc1[ll]}/.P2evensubs/.P2oddsubs/.hl1subs/.JumpsScalar/.zerosubs]],{ll,2,lmax}];
mygrid=Grid[Join[{header},mytable],Frame->All,Alignment->Right];

terms={coB,coC,coD,coF,coE,coAM};
mytable=Table[{ToString[terms[[kk]]],NumberForm[terms[[kk]]/.Jumpmono/.JumpsScalar/.GetABCD[mypar]/.mypar,nsf]},{kk,1,Length[terms]}];
mygrid=Grid[mytable,Frame->All,Alignment->Right];

mytable=Table[Join[{ll},{P2[ll][r],P2[ll]'[r],hl[ll]'[r],Jsc0[ll],Jsc1[ll]}/.P2evensubs/.P2oddsubs/.hl1subs/.JumpsScalar/.zerosubs],{ll,0,lmax}];
Export[datapath<>"data/jumps_"<>ToString[iConfig]<>".mx",mytable];
mytable=Table[{ToString[terms[[kk]]],terms[[kk]]/.Jumpmono/.JumpsScalar/.GetABCD[mypar]/.mypar},{kk,1,Length[terms]}];
Export[datapath<>"data/jumps_completion_"<>ToString[iConfig]<>".mx",mytable];

Print["5. Import radial functions"];

scalarfns=Get[mypath<>"static/static_radial_functions.m"];
{Length[scalarfns],Length[scalarfns[[1]]]};


GetSmoothedFunction[fnname_,ll_,Bco_,par_]:=Module[{Zint,Zext,r0,Asol},
ind=Position[haystack,fnname][[1,1]];
fn=scalarfns[[ind]][[ll-1]];
X0jump=Bco[[2]]*fn[[2]]-Bco[[1]]*fn[[1]]/.{z->z0}/.par;
X1jump=D[Bco[[2]]*fn[[2]]-Bco[[1]]*fn[[1]],z]/.{z->z0}/.par;
Asol=Solve[{X0jump==0,X1jump==0},{Ain,Aup},WorkingPrecision->prec]//First;
{Bco[[1]]*fn[[1]],Bco[[2]]*fn[[2]]}/.Asol/.par
];
MakeCompositeFunction[fns_,par_]:=((HeavisideTheta[z0-z]*fns[[1]]+HeavisideTheta[z-z0]*fns[[2]])/.par);
GetSmoothedFunction["X0",3,{1,2},mypar];

ltry=8;
partry=SetParams[SetPrecision[6,prec],SetPrecision[9/10,prec]];
{"Z0","X0","DeltaB0","\[Chi]ek0","\[Chi]e","\[Chi]oc[-1]","\[Chi]oc[1]","\[Chi]os[-1]","\[Chi]os[1]"};
fns=GetSmoothedFunction["DeltaB0",ltry,{1,1},partry];
myf=MakeCompositeFunction[fns,partry];
(*Plot[Re[myf],{z,SetPrecision[1.01,prec],SetPrecision[30,prec]},WorkingPrecision->prec]*)


Gethlfunction[ll_]:=Module[{},
{LegendreP[ll,z],LegendreQ[ll,0,3,z]}
];
GetP2function[ll_]:=Module[{prefac},
prefac=\[CapitalDelta]^2*D[Getz[r],r]^2/.{r->Getr[z]}/.aMsubs//Simplify;
prefac*{D[LegendreP[ll,z],{z,2}], D[LegendreQ[ll,0,3,z],{z,2}]}
];
(* Check these satisfy the ODEs. *)
ltry=4;
Ztry=Gethlfunction[ltry][[2]];
Simplify[D[(z^2-1)*D[Ztry,z],z]-ltry*(ltry+1)*Ztry];
Ztry=GetP2function[ltry][[2]];
Simplify[(z^2-1)*D[D[Ztry,z],z]-2*z*D[Ztry,z]-(ltry-1)*(ltry+2)*Ztry];


Print["6. Numerical calculation and export"];

CalculateLeftRightCoeffs[nsol1_,mypar_]:=Module[{Ainrepl,Aoutrepl,Ain0subs,Aout0subs},
(* Now we need to work out the coefficients from the jumps. *)
dzdr=D[Getz[r],r];
ABCDsubs=GetABCD[mypar];
For[ll=2,ll<=lmax+6,ll++,
Ainsubs[ll]={};
Aoutsubs[ll]={};
(*Print[ll];*)
If[Mod[ll,2]==0,
(* spin 2 even *)
jumps= {J0e,J1e}/.Jumps2/.\[Lambda]subs/.{l->ll}/.{S0->SpinWeightedSphericalHarmonicY[-2,ll,0,\[Pi]/2,0]}/.ABCDsubs/.mypar;
fns=GetP2function[ll];
eqs={Aout*fns[[2]]-Ain*fns[[1]]==jumps[[1]],dzdr*(Aout*D[fns[[2]],z]-Ain*D[fns[[1]],z])==jumps[[2]]}/.{z->z0}/.mypar;
Asol=Solve[eqs,{Ain,Aout},WorkingPrecision->prec]//First;
(* Print["spin 2 even Asol = " <> Asol]; *)
AppendTo[Ainsubs[ll],Aeven[ll]-> Ain/.Asol];
AppendTo[Aoutsubs[ll],Aeven[ll]-> Aout/.Asol];
(* trace *)
fns=Gethlfunction[ll];
jumps={0,J1h}/.\[Lambda]subs/.{l->ll}/.{Y0->SphericalHarmonicY[ll,0,\[Pi]/2,0]}/.ABCDsubs/.mypar;
eqs={Aout*fns[[2]]-Ain*fns[[1]]==jumps[[1]],dzdr*(Aout*D[fns[[2]],z]-Ain*D[fns[[1]],z])==jumps[[2]]}/.{z->z0}/.mypar;
Asol=Solve[eqs,{Ain,Aout},WorkingPrecision->prec]//First;
(* Print["trace Asol = " <> Asol]; *)
AppendTo[Ainsubs[ll],(Atr[ll]->Ain)/.Asol];
AppendTo[Aoutsubs[ll],(Atr[ll]->Aout)/.Asol];
(* scalar *)
fns=Gethlfunction[ll];
jumps={Jsc0[ll],Jsc1[ll]}/.nsol1;
eqs={Aout*fns[[2]]-Ain*fns[[1]]==jumps[[1]],dzdr*(Aout*D[fns[[2]],z]-Ain*D[fns[[1]],z])==jumps[[2]]}/.{z->z0}/.mypar;
Asol=Solve[eqs,{Ain,Aout},WorkingPrecision->prec]//First;
(* Print["scalar Asol = " <> Asol]; *)
AppendTo[Ainsubs[ll],(Ascal[ll]->Ain)/.Asol];
AppendTo[Aoutsubs[ll],(Ascal[ll]->Aout)/.Asol];
,
(* spin 2 odd *)
jumps={J0o,J1o}/.Jumps2/.\[Lambda]subs/.{l->ll}/.{S0p->D[SpinWeightedSphericalHarmonicY[-2,ll,0,\[Theta],0],\[Theta]]/.{\[Theta]->\[Pi]/2}}/.ABCDsubs/.mypar;
fns=GetP2function[ll];
eqs={Aout*fns[[2]]-Ain*fns[[1]]==jumps[[1]],dzdr*(Aout*D[fns[[2]],z]-Ain*D[fns[[1]],z])==jumps[[2]]}/.{z->z0}/.mypar;
Asol=Solve[eqs,{Ain,Aout},WorkingPrecision->prec]//First;
(* Print["spin 2 odd Asol = " <> Asol]; *)
AppendTo[Ainsubs[ll],Aodd[ll]-> Ain/.Asol];
AppendTo[Aoutsubs[ll],Aodd[ll]-> Aout/.Asol];
];
];

Ain0subs={coB->-(coB/.nsol1),coC->((-2*M*rp^2*coD-(a/6)*(rm^3+rm^2*rp+5*rm*rp^2-7*rp^3)*coF)/.nsol1),coD->(-coD/.nsol1),coE->(-M*coD/.nsol1),coAM->((-a*coD-rp*(rp-rm)^2/(rp+rm)*coF)/.nsol1),coF->-(coF/.nsol1)}/.mypar;(* minus signs in the expressions below because now on the inside. *)
Aout0subs={coB->0,coC->((coC-2*M*rp^2*coD-(a/6)*(rm^3+rm^2*rp+5*rm*rp^2-7*rp^3)*coF)/.nsol1),coD->0,coE-> ((coE-M*coD/.nsol1)/.Jumpmono/.ABCDsubs),coAM->((coAM-a*coD-rp*(rp-rm)^2/(rp+rm)*coF)/.nsol1/.Jumpmono/.ABCDsubs/.mypar),coF->0}/.mypar;

Ainrepl=Join[Ain0subs,Table[Ainsubs[ll],{ll,2,lmax+4}]//Flatten];
Aoutrepl=Join[Aout0subs,Table[Aoutsubs[ll],{ll,2,lmax+4}]//Flatten];
{Ainrepl,Aoutrepl,Ain0subs,Aout0subs}
];



GetAllRadialFunctions[{Ainrepl_,Aoutrepl_},mypar_]:=Module[{fnsin,fnsout},
fnsin={};
fnsout={};
For[ll=2,ll<=lmax+4,ll++,
If[Mod[ll,2]==0,
(* Trace functions *)
fnnames={"Z0","X0","DeltaB0"};
ampls={Atr[ll]/.Ainrepl,Atr[ll]/.Aoutrepl};
dzdr=D[Getz[r],r];
For[fnct=1,fnct<=Length[fnnames],fnct++,
fnname=fnnames[[fnct]];
fns=GetSmoothedFunction[fnname,ll,ampls,mypar];
AppendTo[fnsin,ToExpression[fnname][ll][r]->fns[[1]]];
AppendTo[fnsin,ToExpression[fnname][ll]'[r]->dzdr*D[fns[[1]],z]/.mypar];
AppendTo[fnsout,ToExpression[fnname][ll][r]->fns[[2]]];
AppendTo[fnsout,ToExpression[fnname][ll]'[r]->dzdr*D[fns[[2]],z]/.mypar];
];
AppendTo[fnsin,hl[ll][r]->Atr[ll]*Gethlfunction[ll][[1]]/.Ainrepl];
AppendTo[fnsin,hl[ll]'[r]->Atr[ll]*dzdr*D[Gethlfunction[ll][[1]],z]/.Ainrepl];
AppendTo[fnsout,hl[ll][r]->Atr[ll]*Gethlfunction[ll][[2]]/.Aoutrepl];
AppendTo[fnsout,hl[ll]'[r]->Atr[ll]*dzdr*D[Gethlfunction[ll][[2]],z]/.Aoutrepl];
(* Even parity functions *)
fnnames={"\[Chi]ek0","\[Chi]e"};
ampls={Aeven[ll]/.Ainrepl,Aeven[ll]/.Aoutrepl};
dzdr=D[Getz[r],r];
For[fnct=1,fnct<=Length[fnnames],fnct++,
fnname=fnnames[[fnct]];
fns=GetSmoothedFunction[fnname,ll,ampls,mypar];
AppendTo[fnsin,ToExpression[fnname][ll][r]->fns[[1]]];
AppendTo[fnsin,ToExpression[fnname][ll]'[r]->dzdr*D[fns[[1]],z]/.mypar];
AppendTo[fnsout,ToExpression[fnname][ll][r]->fns[[2]]];
AppendTo[fnsout,ToExpression[fnname][ll]'[r]->dzdr*D[fns[[2]],z]/.mypar];
];
AppendTo[fnsin,P2[ll][r]->Aeven[ll]*GetP2function[ll][[1]]/.Ainrepl];
AppendTo[fnsin,P2[ll]'[r]->Aeven[ll]*dzdr*D[GetP2function[ll][[1]],z]/.Ainrepl];
AppendTo[fnsout,P2[ll][r]->Aeven[ll]*GetP2function[ll][[2]]/.Aoutrepl];
AppendTo[fnsout,P2[ll]'[r]->Aeven[ll]*dzdr*D[GetP2function[ll][[2]],z]/.Aoutrepl];
(* Free scalar functions *)
AppendTo[fnsin,R0[ll][r]->Ascal[ll]*Gethlfunction[ll][[1]]/.Ainrepl];
AppendTo[fnsin,R0[ll]'[r]->Ascal[ll]*dzdr*D[Gethlfunction[ll][[1]],z]/.Ainrepl];
AppendTo[fnsout,R0[ll][r]->Ascal[ll]*Gethlfunction[ll][[2]]/.Aoutrepl];
AppendTo[fnsout,R0[ll]'[r]->Ascal[ll]*dzdr*D[Gethlfunction[ll][[2]],z]/.Aoutrepl];
,
(* Odd parity functions *)
fnnames={"\[Chi]oc[-1]","\[Chi]oc[1]","\[Chi]os[-1]","\[Chi]os[1]"};
ampls={Aodd[ll]/.Ainrepl,Aodd[ll]/.Aoutrepl};
dzdr=D[Getz[r],r];
For[fnct=1,fnct<=Length[fnnames],fnct++,
fnname=fnnames[[fnct]];
fns=GetSmoothedFunction[fnname,ll,ampls,mypar];
AppendTo[fnsin,ToExpression[fnname][ll][r]->fns[[1]]];
AppendTo[fnsin,ToExpression[fnname][ll]'[r]->dzdr*D[fns[[1]],z]/.mypar];
AppendTo[fnsout,ToExpression[fnname][ll][r]->fns[[2]]];
AppendTo[fnsout,ToExpression[fnname][ll]'[r]->dzdr*D[fns[[2]],z]/.mypar];
];
AppendTo[fnsin,P2[ll][r]->Aodd[ll]*GetP2function[ll][[1]]/.Ainrepl];
AppendTo[fnsin,P2[ll]'[r]->Aodd[ll]*dzdr*D[GetP2function[ll][[1]],z]/.Ainrepl];
AppendTo[fnsout,P2[ll][r]->Aodd[ll]*GetP2function[ll][[2]]/.Aoutrepl];
AppendTo[fnsout,P2[ll]'[r]->Aodd[ll]*dzdr*D[GetP2function[ll][[2]],z]/.Aoutrepl];
];
];
{fnsin,fnsout}
];


VatJumpsNum=SetNumericalJumps[mypar];
{nsol1,checks1}=SolveJumpEquations[VatJumpsNum];
(* Print[checks1]; *)

(* Quiet[{Ainrepl,Aoutrepl,Ain0subs,Aout0subs}=CalculateLeftRightCoeffs[nsol1,mypar]];
Quiet[{fnsin,fnsout}=GetAllRadialFunctions[{Ainrepl,Aoutrepl},mypar]]; *)

{Ainrepl,Aoutrepl,Ain0subs,Aout0subs}=CalculateLeftRightCoeffs[nsol1,mypar];
{fnsin,fnsout}=GetAllRadialFunctions[{Ainrepl,Aoutrepl},mypar];


ExtendVat[V_]:=Module[{},
dim=Dimensions[V];
row=Table[SetPrecision[0,prec],{ll,0,dim[[2]]-1}];
For[ll=0,ll<=lmax+6,ll++,
If[ll<=lmax&&Mod[ll,2]==0,
row[[ll+1]]=r^4*Vat[[5,ll+1]]-Vat[[4,ll+1]]
+2*a^2*r^2*Sum[C0c[2][ll][kk]*Vat[[5,ll+kk+1]],{kk,-Min[ll,2],2,2}]
+a^4*Sum[C0c[4][ll][kk]*Vat[[5,ll+kk+1]],{kk,-Min[ll,4],4,2}];
;
];
];
V2=Append[Vat,row/.cosubs];
V2
];

ToTenComponents[in_,ll_]:=Module[{fac=-1,out},
(* For m=0, there are five independent degrees of freedom. Above, I have extended the Vat to six rows. 
However, we want output data in the same format as for the m != 0 modes, for which we use a set of 10 metric components. Here I apply the symmetries of the MP to build the set of 10.
*)
(* Here fac = -1 is because, for some reason when I first wrote the code, I projected the \[Rho] h_{l+ m+} component onto the spin = -1 harmonic . This can be easily recitified here, since Y_{-1} = - Y_{+1} for m=0. *)
If[Mod[ll,2]==0,
out={Re[in[[1]]],Re[in[[1]]],Re[in[[2]]],Re[in[[2]]],fac*Re[in[[3]]],-fac*Re[in[[3]]],fac*Re[in[[3]]],-fac*Re[in[[3]]],Re[in[[6]]],Re[in[[5]]]};
,
out={0,0,0,0,fac*I*Im[in[[3]]],fac*I*Im[in[[3]]],-fac*I*Im[in[[3]]],-fac*I*Im[in[[3]]],0,0};
];
out
];


TortoiseCoord[r_]:=
r+(rp+rm)/(rp-rm)*(rp*Log[(r-rp)/2]-rm*Log[(r-rm)/2])/.mypar;
InvTortoise[x_?NumericQ]:=r/.FindRoot[TortoiseCoord[r]-x,{r,(rp+SetPrecision[10^(-5),prec])/.mypar},WorkingPrecision->prec-2];
InvTortoise[TortoiseCoord[6]];
(* Note that the coefficients have been absorbed into the definition of the functions, as is required for the smoothing process. So now we set the coefficients to unity. *)
Getznum[rr_?NumericQ]=(2*rr-rp-rm)/(rp-rm)/.mypar;
Aunitysubs=Table[{Atr[ll]->1,Aeven[ll]->1,Aodd[ll]->1,Ascal[ll]->1},{ll,2,lmax+6}]//Flatten;

r0star=TortoiseCoord[r0val];
If[rgrid==0, (* Is the grid in r* or in r? *)
rmin2=InvTortoise[rstarmin];
rmax2=InvTortoise[rstarmax];
,
rmin2=rstarmin;
rmax2=rstarmax;
];


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
nleft=Floor[(r0val-rstarmin)/dr];
nright=Floor[(rstarmax-r0val)/dr];
rpts=nleft+nright+1;
rsL=Table[SetPrecision[r0val-ri*dr,prec],{ri,0,nleft}];
rsR=Table[SetPrecision[r0val+ri*dr,prec],{ri,0,nright}];
rs=Table[SetPrecision[r0val+(ri-1-nleft)*dr,prec],{ri,1,rpts}];
rstarsL=Map[TortoiseCoord,rsL];
rstarsR=Map[TortoiseCoord,rsR];
rstars=Map[TortoiseCoord,rs];
];


filen=datapath<>"data/lm_rs_"<>ToString[iConfig]<>".dat";
Export[filen,Transpose@{rstars,rs}];

filen=datapath<>"data/lm_qs_"<>ToString[iConfig]<>".dat";
Export[filen,qs];
filen=datapath<>"data/lm_rsL_"<>ToString[iConfig]<>".dat";
Export[filen,Transpose@{rstarsL,rsL}];
filen=datapath<>"data/lm_rsR_"<>ToString[iConfig]<>".dat";
Export[filen,Transpose@{rstarsR,rsR}];


(* Check continuity. Pick arbitrary modes *)
qi=5;ll=2;
fin=N[Vat[[qi]][[ll+1]]/.reverselogrepl/.Ain0subs/.Aunitysubs/.fnsin/.{r->Getr[z]}/.mypar,prec];
fout=N[Vat[[qi]][[ll+1]]/.reverselogrepl/.Aout0subs/.Aunitysubs/.fnsout/.{r->Getr[z]}/.mypar,prec];
{fin,fout,fin-fout}/.{z->z0}/.mypar;
Ain0subs;


Vat2=ExtendVat[Vat];
GetData[ll_]:=Module[{dataL,dataR},
Vatin=Vat2/.reverselogrepl/.Ain0subs/.Aunitysubs/.fnsin;
Vatout=Vat2/.reverselogrepl/.Aout0subs/.Aunitysubs/.fnsout;
dataL=Monitor[Table[ToTenComponents[Table[Vatin[[qi]][[ll+1]]/.{r->rsL[[ri]],z->Getznum[rsL[[ri]]]}/.mypar,{qi,1,6}],ll],{ri,1,Length[rsL]}],ri];
dataR=Monitor[Table[ToTenComponents[Table[Vatout[[qi]][[ll+1]]/.{r->rsR[[ri]],z->Getznum[rsR[[ri]]]}/.mypar,{qi,1,6}],ll],{ri,1,Length[rsR]}],ri];
{dataL,dataR}
];
GetAllData[]:=Module[{allL,allR},
allL=Table[0,{ri,1,Length[rsL]},{qi,1,10},{ll,1,lmax+1}];
allR=Table[0,{ri,1,Length[rsR]},{qi,1,10},{ll,1,lmax+1}];
For[ll=0,ll<=lmax,ll++,
{dataL,dataR}=GetData[ll];
For[ri=1,ri<=Length[rsL],ri++,
For[qi=1,qi<=10,qi++,
allL[[ri,qi,ll+1]]=dataL[[ri,qi]];
];
];
For[ri=1,ri<=Length[rsR],ri++,
For[qi=1,qi<=10,qi++,
allR[[ri,qi,ll+1]]=dataR[[ri,qi]];
];
];
];
{allL,allR}
];
(* Save in a binary format. *)
GetArray[tbl_]:=Module[{grid},
dim=Dimensions[tbl];
arr=Developer`ToPackedArray[Table[SetPrecision[0,prec],{dim[[1]]},{2*dim[[2]]},{dim[[3]]}]];
For[ri=1,ri<=dim[[1]],ri++,
For[qi=1,qi<=dim[[2]],qi++,
For[ll=0,ll<dim[[3]],ll++,
arr[[ri,qi,ll+1]]=Re[tbl[[ri,qi,ll+1]]];
arr[[ri,10+qi,ll+1]]=Im[tbl[[ri,qi,ll+1]]];
];
];
];
arr
];

{dataL,dataR}=GetData[8];
{allL,allR}=GetAllData[];

dataL = GetArray[allL];
dataR = GetArray[allR];

(* Export the data in a binary format. *)

(* clearNumberQ[x_] := NumericQ[x] && FiniteQ[x] *)

(* If[ ! ArrayQ[dataL, clearNumberQ],
  
  (* if any non-numerics are lurking, report their positions *)
  badPos = Position[dataL, x_ /; ! clearNumberQ[x], {3}];
  Print[
    "Non-numeric entries detected in dataL at (ri, qi, l) = ",
    badPos
  ];
  
  (* you can choose to abort or skip this iConfig *)
  Abort[];
  
]; *)

fn=datapath<>"data/lm_in_"<>ToString[iConfig]<>".bin";
Export[fn,dataL,dformat];

(* If[ ! ArrayQ[dataR, clearNumberQ],
  
  (* if any non-numerics are lurking, report their positions *)
  badPos = Position[dataR, x_ /; ! clearNumberQ[x], {3}];
  Print[
    "Non-numeric entries detected in dataR at (ri, qi, l) = ",
    badPos
  ];
  
  (* you can choose to abort or skip this iConfig *)
  Abort[];
  
]; *)

fn=datapath<>"data/lm_up_"<>ToString[iConfig]<>".bin";
Export[fn,dataR,dformat];

];
