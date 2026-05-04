(* ::Package:: *)

Needs["KerrPunctures`"];


BeginPackage["h1Functionality`",{"SpinWeightedSpheroidalHarmonics`", "Developer`", "SimulationTools`", "KerrPunctures`","NContinuedFraction`","SpinWeightedSpheroidalHarmonicsFT`","SWSHdecomp`" }];


"This is a package that containes functionality realting to transformations of handling and regularising radial h1 data contructing the 2D m-mode data and cooridinate transformations on said data "


RetSamDImportRetFunction:"LoadRetSolutions"
RetKevinImportRetFunctionInteral:"Kevins Import function"
PunctureImportPuncFunction:"LoadPuncSolutions"
ResiduleFiedlConstrution:"ResiduleFiedlConstrution[a, rp, RetDataStruct, PuncDataStruct]"
SphericalHamronicsContructor:"Constructor"
MakeResGrid:"ResGridMaker"
h1$tetrad$data:"PuncEval"
h1$tetrad$dataStruct:"Datastruc"
Punc\[Theta]Window:"Theta Pucnture"
h1$Punc$dataStruct:"Punc data struct"
h1$Res$dataStruct:"ResDataStruct"
h1$Window$Punc$dataStruct:"WindowData"
TensorConstruct2DData:"consturct2ddata"
DataStructKinnerslyTetradToBLComp:"Converter"
h1$Ret$Punc$dataStruct:"RetPuncdat"
RetStaticSamDImportRetFunction:""
RetRadiativeKevinImportRetFunction:""
JoininRadiativeandNonePieces:""
h1$Res$mmode$dataStruct:""
h1$Punc$mmode$dataStruct:""
RetRadiativeBenDImportRetFunction:""
JoininRadiativeandStaticPieces:""
RetRadiativeMaartenImportRetFunction:""
RetRadiativeSamDImportRetFunction:""
RetRadiativeMaartenImportRetFunction:""
TensorConstruct2DDataSave:""


Begin["`Private`"];
Get[FindFile["KerrPunctures`"]]
Get[FindFile["SimulationTools`"]]


RetKevinImportRetFunctionInteral[placeholder_]:= Module[{},
Return["placeholder"]
]


RetSamDImportRetFunctionInteral[a_, rp_,m_,filepathin_]:= Module[{filepathMaarten,filen,GetIndexL,GetIndexR,GetParam,a0,iround,prec,configparams,s2,s1,s,subrs,filepath,iConfig,compR,compL,psis,subrstars,configpath,datagpath,r0,mm,lmax,lplot,dformat,rstarsL,rsL,tmp,rstarsR,rsR,rstarmin,rstarmax,fn,dataL,dimsL,dataR,dimsR},

iConfig = ToString[Round[1000*a]]<>"-"<>ToString[Round[10*rp]]<>"-00"<>ToString[m];
filepath=filepathin<>("/a_"<>ToString[NumberForm[a,2]]<>"_r0_"<>ToString[NumberForm[rp,{Infinity,1}]]<>"/");
configpath=filepath<>"config/";
datagpath=filepath<>"data/";
filen=configpath<>"config"<>ToString[iConfig]<>".txt";
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
prec=GetParam["prec",32]; (* Number of digits to use where required. *)
iround=1000;
a0=SetPrecision[Round[GetParam["a"]*iround]/iround,prec];
r0=SetPrecision[Round[GetParam["r0"]*iround]/iround,prec];
mm=GetParam["m"];
lmax=GetParam["lmax"];
lplot=GetParam["lplot"];
dformat=ToString@GetParam["dformat","Real32"];

tmp=Import[datagpath<>"lm_rsL_"<>ToString[iConfig]<>".dat"];
{rstarsL,rsL}=Transpose[tmp];
tmp=Import[datagpath<>"lm_rsR_"<>ToString[iConfig]<>".dat"];

{rstarsR,rsR}=Transpose[tmp];
rstarmin=Last[rstarsL];
rstarmax=Last[rstarsR];

configpath=filepath<>"config/";
datagpath=filepath<>"data/";

fn=datagpath<>"lm_in_"<>ToString[iConfig]<>".bin";
dataL=Import[fn,dformat];
dimsL={Length[rsL],20,lmax+1};
GetIndexL[ri_,qi_,ll_]:=(ri-1)*dimsL[[2]]*dimsL[[3]]+(qi-1)*dimsL[[3]]+ll+1;
compL=Table[dataL[[GetIndexL[ri,qi,ll]]]+I*dataL[[GetIndexL[ri,qi+10,ll]]],{ll,0,lmax},{qi,1,10},{ri,1,Length[rsL]}] ;

fn=datagpath<>"lm_up_"<>ToString[iConfig]<>".bin";
dataR=Import[fn,dformat];
dimsR={Length[rsR],20,lmax+1};

GetIndexR[ri_,qi_,ll_]:=(ri-1)*dimsR[[2]]*dimsR[[3]]+(qi-1)*dimsR[[3]]+ll+1;
compR=Table[dataR[[GetIndexR[ri,qi,ll]]]+I*dataR[[GetIndexR[ri,qi+10,ll]]],{ll,0,lmax},{qi,1,10},{ri,1,Length[rsR]}] ;
subrs = Flatten[Join[{Reverse[rsL],rsR[[2;;]]}]];

psis=Table[Flatten[Join[{Reverse[compL[[l+1,q]]],compR[[l+1,q,2;;]]} ]],{l,0,lmax},{q,1,10}];

Return[Association[ {"h1Dat"-> psis,"rvals"-> subrs,"lmax"->lmax}]]
]


JoininRadiativeandStaticPieces[aVal_, rOrbit_, mmax_, lmax_, RadStruct_, StaticStruct_] := 
Module[{h1Radiative, rvals, rvalsStatic, staticIndices, Joineddata, InheritedFields, h1static,staticIndicesOuter,staticIndicesInner},

rvals       = "rvals" /. RadStruct;
rvalsStatic = "rvals" /. StaticStruct;

(* Automatically find which static grid indices align with the radiative grid *)
staticIndicesInner = Nearest[rvalsStatic -> "Index", rvals[[1]]][[1]];
staticIndicesOuter = Nearest[rvalsStatic -> "Index", rvals[[-1]]][[1]];

h1Radiative = PadLeft[Transpose[("h1Data1D" /. RadStruct), 1 <-> 2][[;; mmax, ;; lmax , ;;]], Dimensions[Transpose[("h1Data1D" /. RadStruct), 1 <-> 2][[;; mmax, ;; lmax , ;;]]]+{0,1,0,0}];
h1static    = ("h1Data1D" /. StaticStruct)[[;; lmax + 1, ;;, staticIndicesInner;;staticIndicesOuter]];
Joineddata     = Join[{h1static}, h1Radiative];
InheritedFields = KeyDrop[RadStruct, {"h1Data1D", "rvals"}];

  Return[Join[
    <|"h1Data1D" -> Transpose[Joineddata, 1 <-> 2], "rvals" -> rvals|>,
    InheritedFields
  ]];
]


RetRadiativeSamDImportRetFunction[a_, rp_,mmax_,lmax_,filepathin_]:= Module[{tmph1rettmp,tmp0,rvals2,h1Radiative,hTraceSphericallm,tmp,tmph1ret,tmprvals ,tmplmax,tmplmaxNew},

tmp = Table[ RetSamDImportRetFunctionInteral[a, rp,m,filepathin],{m,0,mmax}];

tmph1rettmp  = "h1Dat"/.tmp;


tmph1ret = Join[{tmph1rettmp[[1,;;,;;,;;]]},tmph1rettmp[[2;;, ;;, ;;, 2;;]] ];



tmprvals  = "rvals"/.tmp;
tmplmax  = ("lmax"/.tmp)[[3]];

(*Transforms trace from spheroidal to spherical*)
hTraceSphericallm[aa_,r00_,m_]:= Module[{lSumMax,IlmVec,Mixing,Jlms, JlmsBuff},
lSumMax = tmplmax;
IlmVec = Table[ tmph1ret[[m+1,\[Alpha]\[Alpha]+1,10]],{\[Alpha]\[Alpha],Abs[m],lSumMax}] ;
Mixing = SWSHb[0,m,(aa m )/(r00^(3/2)+ a), "lmax"->lSumMax];

(*Note i blelive this shoudl be Conjugate Mix*)
Jlms  = Table[Sum[ Mixing[[\[Beta]\[Beta],lli]] IlmVec[[\[Beta]\[Beta]]] ,{\[Beta]\[Beta],Length[IlmVec]}] ,{lli,Length[IlmVec]}];
Return[Jlms]
];

tmplmaxNew = tmph1ret;

Table[tmplmaxNew[[mi+1, mi +1 ;;, 10]] = hTraceSphericallm[a,rp,mi],{mi,0,mmax}];


h1Radiative =Transpose[ (tmplmaxNew)[[;;, ;;lmax+1, ;;, 10;;180 ]] ,1<->2];

rvals2 = (tmprvals)[[;;,12;;182]];

Return[Association[ {"h1Data1D"-> (h1Radiative)[[;;,2;;]],"rvals"-> (rvals2)[[;;]],"lmax"->tmplmax,"mmax"->mmax,"a0"->a, "rOrbit"->rp,"DataType"->"1RadialData"}]]
]


RetRadiativeMaartenImportRetFunction[a_, rp_,mmax_,lmax_,filepathin_]:= Module[{tmph1rettmp,tmp0,rvals2,h1Radiative,hTraceSphericallm,tmp,tmph1ret,tmprvals ,tmplmax,tmplmaxNew},

tmp = Table[ If[m!=2,0*RetMaartenImportRetFunctionInteral[a, rp,2,filepathin], RetMaartenImportRetFunctionInteral[a, rp,m,filepathin]],{m,0,mmax}];

tmph1rettmp  = "h1Dat"/.tmp;


tmph1ret = Join[{tmph1rettmp[[1,;;,;;,;;]]},tmph1rettmp[[2;;, ;;, ;;, 2;;]] ];



tmprvals  = "rvals"/.tmp;
tmplmax  = ("lmax"/.tmp)[[3]];

(*Transforms trace from spheroidal to spherical*)
hTraceSphericallm[aa_,r00_,m_]:= Module[{lSumMax,IlmVec,Mixing,Jlms, JlmsBuff},
lSumMax = tmplmax;
IlmVec = Table[ tmph1ret[[m+1,\[Alpha]\[Alpha]+1,10]],{\[Alpha]\[Alpha],Abs[m],lSumMax}] ;
Mixing = SWSHb[0,m,(aa m )/(r00^(3/2)+ a), "lmax"->lSumMax];

(*Note i blelive this shoudl be Conjugate Mix*)
Jlms  = Table[Sum[ Mixing[[\[Beta]\[Beta],lli]] IlmVec[[\[Beta]\[Beta]]] ,{\[Beta]\[Beta],Length[IlmVec]}] ,{lli,Length[IlmVec]}];
Return[Jlms]
];

tmplmaxNew = tmph1ret;

Table[tmplmaxNew[[mi+1, mi +1 ;;, 10]] = hTraceSphericallm[a,rp,mi],{mi,0,mmax}];


h1Radiative =Transpose[ (tmplmaxNew)[[;;, ;;lmax+1, ;;, 10;;180 ]] ,1<->2];

rvals2 = (tmprvals)[[;;,12;;182]];

Return[Association[ {"h1Data1D"-> (h1Radiative)[[;;,2;;]],"rvals"-> (rvals2)[[;;]],"lmax"->tmplmax,"mmax"->mmax,"a0"->a, "rOrbit"->rp,"DataType"->"1RadialData"}]]
]


numStr[x_]:=Module[{s=ToString[x]},Which[StringMatchQ[s,___~~"."],s<>"0",!StringContainsQ[s,"."],s<>".0",True,s ]]

LoadMaartenFunction[mi_,mmax_,lmax_,index_,a_,rOrbit_,fullpath_]:=Module[{rvals,loadpath,hIn,hUp,rIn,rUp,rIndata$Test$rIn,data$Test$rUp,datavalsNew,miNew,indexNew,comps,datavals,radialVals,data$Test$In$Re,data$Test$In$Im,data$Test$Up$Re,data$Test$Up$Im,filename},

loadpath = fullpath <> "/h1_a" <> ToString[NumberForm[Round[a, 0.1], {Infinity, 1}]] <> 
           "_rp" <> ToString[NumberForm[Round[rOrbit, 0.1], {Infinity, 1}]] <> 
           "_l" <> ToString[lmax] <> "_m" <> ToString[mi] <> ".h5";

comps={"h_l+l+","h_l-l-","h_m+m+","h_m-m-","rho_h_l+m+","rhob_h_l+m-","rhob_h_l-m+","rho_h_l-m-","sigma_delta_h_l+l-","h"};

rIn=Import[loadpath,{"HDF5","/m_"<>ToString[mi]<>"/r_in"}];
rUp=Import[loadpath,{"HDF5","/m_"<>ToString[mi]<>"/r_up"}];

hIn=Import[loadpath,{"HDF5","/m_"<>ToString[mi]<>"/In/"<>comps[[index]]}];
hUp=Import[loadpath,{"HDF5","/m_"<>ToString[mi]<>"/Up/"<>comps[[index]]}];

rvals = Join[Reverse[rIn][[;;-2]],rUp];

datavalsNew = Table[Join[Reverse[(("Re"/.hIn) + I ("Im"/.hIn ))[[li+1]]][[;;-2]],(("Re"/.hUp) + I ("Im"/.hUp ))[[li+1]]],{li,lmax}];

Return[{rvals,datavalsNew}]

];




RetRadiativeMaartenImportRetFunction[a_, rp_,mmax_,lmax_,filepathin_]:= Module[{fullpath,tmph1rettmp,tmp0,rvals2,h1Radiative,hTraceSphericallm,tmp,tmph1ret,tmprvals ,tmplmax,tmplmaxNew},

h1Radiative = Table[ LoadMaartenFunction[mi,mmax,lmax,index,a,rp,filepathin], {mi,mmax}, {index,10}];

Return[Association[ {"h1Data1D"-> Transpose[Transpose[ (h1Radiative)[[;;,;;,2]],2<->3],1<->2],"rvals"-> (h1Radiative)[[1,1,1]],"lmax"->tmplmax,"mmax"->mmax,"a0"->a, "rOrbit"->rp,"DataType"->"1RadialData"}]]
]


RetStaticSamDImportRetFunction[a_, rp_,mmax_,lmax_,filepathin_]:= Module[{hTraceSphericallm,tmp,tmph1ret,tmprvals ,tmplmax,tmplmaxNew},


tmp = RetSamDImportRetFunctionInteral[a, rp,0,filepathin] ;

tmph1ret  = "h1Dat"/.tmp;
tmprvals  = "rvals"/.tmp;


tmplmaxNew = tmph1ret;


Return[Association[ {"h1Data1D"-> tmplmaxNew,"rvals"-> tmprvals,"lmax"->tmplmax,"mmax"->mmax,"a0"->a, "rOrbit"->rp,"DataType"->"1RadialData"}]]
]


PunctureImportPuncFunction[a_, rp_,mmax_,lmax_,filepath_, rvalsRet_,{rMin_,rMax_}]:=Module[{lhszeros,punctuerdata,rvalsPunc,puncturer1closestelementReduced,punctureRclosestelementReduced,rhszeros,sws,hS1lmModesWindowed,puncturer1,puncturer1closestelementReducedpunc,punctureRclosestelementReducedpunc,puncturer1indexReducedpunc,puncturerR,puncturer1closestelement,punctureRclosestelement,puncturer1index,punctureRindexReducedpunc,punctureRindex},
rvalsPunc = (Get[ (filepath<> "/a_"<>ToString[NumberForm[a,1]]<>"_r0_"<>ToString[NumberForm[rp,{Infinity,1}]]<>"/rvals.dat")]);
If[rvalsPunc[[1]]>rMin,Return["Error: Requested puncutre window larger than generated data. Allowed range, "<>ToString[rvalsPunc[[1]]]<>"<=r<="<>ToString[rvalsPunc[[-1]]]]];
If[rvalsPunc[[-1]]<rMax,Return["Error: Requested puncutre window larger than generated data. Allowed range, "<>ToString[rvalsPunc[[1]]]<>"<=r<="<>ToString[rvalsPunc[[-1]]]]];


puncturer1closestelement = Nearest[rvalsRet->"Element",rMin][[1]];
punctureRclosestelement = Nearest[rvalsRet->"Element",rMax][[1]];
puncturer1index=Position[rvalsRet,puncturer1closestelement][[1,1]];
punctureRindex=Position[rvalsRet,punctureRclosestelement][[1,1]];

puncturer1closestelementReducedpunc= Nearest[rvalsPunc->"Element",rMin][[1]];
punctureRclosestelementReducedpunc  = Nearest[rvalsPunc->"Element",rMax][[1]];
puncturer1indexReducedpunc=Position[rvalsPunc,puncturer1closestelementReducedpunc][[1,1]];
punctureRindexReducedpunc=Position[rvalsPunc,punctureRclosestelementReducedpunc][[1,1]];

sws = {0,0,2,-2,1,-1,1,-1,0,0};

lhszeros = ConstantArray[0,puncturer1index-1];
rhszeros = ConstantArray[0,Length[rvalsRet] - punctureRindex];

punctuerdata=  Table[Table[Join[
lhszeros,
If[l< Abs[m]|| l< Abs[sws[[i]]], 0* rvalsRet[[;;Length[rvalsRet]-(Length[lhszeros]+ Length[rhszeros])]],

hS1lmModesWindowed=(Get[ (filepath<> "/a_"<>ToString[NumberForm[a,1]]<>"_r0_"<>ToString[NumberForm[rp,{Infinity,1}]]<>"/h1Punc_l_"<>ToString[l]<>"_m_"<>ToString[m]<>".dat")]);
hS1lmModesWindowed[[i]][[puncturer1indexReducedpunc;;punctureRindexReducedpunc]]]
,rhszeros],{i,10}],{m,0,mmax},{l,0,lmax}];


Return[Association[ {"h1Data1D"-> Transpose[punctuerdata,1<->2],"rvals"-> rvalsRet ,"lmax"->lmax,"mmax"->mmax,"a0"->a, "rOrbit"->rp, "windowRadialvals"->{rMin,rMax}, "windowRadialindex"->{puncturer1index,punctureRindex},"DataType"->"RadialData"}]]

]



ResiduleFiedlConstrution[a_, rp_,mmax_,lmax_, RetDataStruct_, PuncDataStruct_]:=  Module[{ResiduleDataTest,InheritedFields},

ResiduleDataTest=("h1Data1D"/.RetDataStruct)-("h1Data1D"/.PuncDataStruct);
InheritedFields = KeyDrop[PuncDataStruct,"h1Data1D"];
Return[Join[<|"h1Data1D" -> ResiduleDataTest|>,InheritedFields]]
]


SphericalHamronicsContructor[lmax_,mmax_,\[Theta]vals_]:= Module[{tmp,SpinZeroSpheroidalharmonic},
tmp = Quiet[Table[If[ll<Max[Abs[ss],Abs[mm]],0*\[Theta]vals,SpinWeightedSphericalHarmonicY[ss,ll,mm,\[Theta]vals,0]],{ss,-2,2},{ll,0,lmax},{mm,0,mmax}]];
Return[Association[{"Ys"->tmp,"Slm"->SpinZeroSpheroidalharmonic, "lmax"->lmax, "mmax"->mmax,"\[Theta]vals"->\[Theta]vals}]]
]


MakeResGrid[mmax_,lmaxVal_,ResDataStruct_,YsStruc_]:=Module[{h1TraceSpherical,hTraceSphericallm,\[Theta]vals,Ysdata,Inheritedfields, rvals,dataInput,residlistl,residlist,Yslist,arrL,arrR,arrLsum,arrRsum,arrLR,arr, arrFull,arrlll,sws},
sws = {0,0,2,-2,1,-1,1,-1,0,0};

Ysdata = N["Ys"/.YsStruc];
\[Theta]vals = "\[Theta]vals"/.YsStruc;
dataInput = N["h1Data1D"/.ResDataStruct];

Inheritedfields =  KeyDrop[ResDataStruct,{"h1PuncDat","lmax"}];

arrFull= Table[Sum[Outer[Times,dataInput[[ll+1,mm+1,icomp]],Ysdata[[sws[[icomp]]+3,ll+1, mm+1]]],{ll,Abs[mm],If[mm==0,30,lmaxVal]}],{icomp,1,10},{mm,0,mmax}];
Return[Join[<|"Data2D"->arrFull, "lmaxVal"->lmaxVal,"\[Theta]vals"->\[Theta]vals|>,Inheritedfields]];
];


(*Internal Functions for setting up correct data format*)


h1$tetradSamtoKinnConvert[{ai_,bi_},mi_,Residuleh1DataStruc2D_,{a0_}]:=Module[{\[Theta]vals,Residuleh1$Data$Sam$2D,rvals,r,\[Theta],q,rgrid ,\[Theta]grid,sin\[Theta]grid ,\[CapitalDelta]grid,\[CapitalSigma]grid,\[Rho]grid,coefgrid$vec,coef$Mat},

Residuleh1$Data$Sam$2D = "Data2D"/.Residuleh1DataStruc2D;

rvals = N[("rvals"/.Residuleh1DataStruc2D)];
\[Theta]vals = N[("\[Theta]vals"/.Residuleh1DataStruc2D)];

rgrid =  ToPackedArray[Table[rvals[[i]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];
\[Theta]grid =  ToPackedArray[Table[\[Theta]vals[[j]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];

\[Rho]grid= (r+ I a0 Cos[\[Theta]])/.{r->rgrid, \[Theta]->\[Theta]grid};
\[CapitalSigma]grid= (r^2+a0^2 Cos[\[Theta]]^2)/.{r->rgrid, \[Theta]->\[Theta]grid};
\[CapitalDelta]grid= (r^2-2r + a0^2)/.{r->rgrid, \[Theta]->\[Theta]grid};

coefgrid$vec= {1, ( - \[CapitalDelta]grid/(2 \[CapitalSigma]grid)), 1/(Sqrt[2]\[Rho]grid), 1/( Sqrt[2]Conjugate[\[Rho]grid])} ;
coef$Mat =Table[ coefgrid$vec[[i]]coefgrid$vec[[j]],{i,1,4},{j,1,4}];


If[{ai,bi}=={1,1}, Return[coef$Mat[[1,1]] Residuleh1$Data$Sam$2D[[1,mi+1]]]];

If[{ai,bi}=={1,2}||{ai,bi}=={2,1}, Return[  1/(\[CapitalSigma]grid * \[CapitalDelta]grid) coef$Mat[[1,2]] Residuleh1$Data$Sam$2D[[9,mi+1]]]];

If[{ai,bi}=={1,3}||{ai,bi}=={3,1}, Return[1/\[Rho]grid coef$Mat[[1,3]]  Residuleh1$Data$Sam$2D[[5,mi+1]]]];

If[{ai,bi}=={1,4}||{ai,bi}=={4,1}, Return[ 1/Conjugate[\[Rho]grid] coef$Mat[[1,4]]   Residuleh1$Data$Sam$2D[[6,mi+1]]]];

If[{ai,bi}=={2,2}, Return[coef$Mat[[2,2]]   Residuleh1$Data$Sam$2D[[2,mi+1]]]];

If[{ai,bi}=={2,3}||{ai,bi}=={3,2}, Return[ 1/Conjugate[\[Rho]grid] coef$Mat[[2,3]]  Residuleh1$Data$Sam$2D[[7,mi+1]]]];

If[{ai,bi}=={2,4}||{ai,bi}=={4,2}, Return[1/\[Rho]grid coef$Mat[[2,4]]   Residuleh1$Data$Sam$2D[[8,mi+1]]]];

If[{ai,bi}=={3,3}, Return[ coef$Mat[[3,3]]Residuleh1$Data$Sam$2D[[3,mi+1]]]];

If[{ai,bi}=={3,4}||{ai,bi}=={4,3}, Return[  coef$Mat[[3,4]] (\[CapitalSigma]grid Residuleh1$Data$Sam$2D[[10,mi+1]]-1/\[CapitalSigma]grid Residuleh1$Data$Sam$2D[[9,mi+1]])]];

If[{ai,bi}=={4,4}||{ai,bi}=={4,4}, Return[  coef$Mat[[4,4]] Residuleh1$Data$Sam$2D[[4,mi+1]]]];

]





h1$Punc$Tetrad$dataInteral[{a_,b_},m_,{rGrid_,\[Theta]Grid_},{aVal_,rOrbit_}]:=Module[{M0,q,rgrid ,\[Theta]grid,sin\[Theta]grid ,\[CapitalDelta]grid,\[CapitalSigma]grid,\[Rho]grid,coefgrid$vec,coef$Mat},

M0 = 1;

If[{a,b}=={1,1}, Return[ KerrPunctureHS12D[M0,aVal,rOrbit][1,rGrid,\[Theta]Grid,m] ] ];

If[{a,b}=={1,2}||{a,b}=={2,1}, Return[ KerrPunctureHS12D[M0,aVal,rOrbit][2,rGrid,\[Theta]Grid,m] ]];

If[{a,b}=={1,3}||{a,b}=={3,1}, Return[  KerrPunctureHS12D[M0,aVal,rOrbit][3,rGrid,\[Theta]Grid,m]]];

If[{a,b}=={1,4}||{a,b}=={4,1}, Return[ KerrPunctureHS12D[M0,aVal,rOrbit][4,rGrid,\[Theta]Grid,m]] ];

If[{a,b}=={2,2}, Return[ KerrPunctureHS12D[M0,aVal,rOrbit][5,rGrid,\[Theta]Grid,m] ]];

If[{a,b}=={2,3}||{a,b}=={3,2}, Return[  KerrPunctureHS12D[M0,aVal,rOrbit][6,rGrid,\[Theta]Grid,m] ]];

If[{a,b}=={2,4}||{a,b}=={4,2}, Return[  KerrPunctureHS12D[M0,aVal,rOrbit][7,rGrid,\[Theta]Grid,m] ]];

If[{a,b}=={3,3}, Return[ KerrPunctureHS12D[M0,aVal,rOrbit][8,rGrid,\[Theta]Grid,m]]];

If[{a,b}=={3,4}||{a,b}=={4,3}, Return[  KerrPunctureHS12D[M0,aVal,rOrbit][9,rGrid,\[Theta]Grid,m] ]];

If[{a,b}=={4,4}||{a,b}=={4,4}, Return[ KerrPunctureHS12D[M0,aVal,rOrbit][10,rGrid,\[Theta]Grid,m]]];


]


(*Functions for constructing datastruct for single tetrad component*)


h1$Res$dataStruct[{ai_,bi_},mi_,Residuleh1DataStruc2D_,{a0_,r0_}]:=Module[{\[Theta]grid,rgrid,Inheritedfields,tmp,rvals,\[Theta]vals,dataoutput},

tmp = h1$tetradSamtoKinnConvert[{ai,bi},mi,Residuleh1DataStruc2D,{a0}];

rvals = ("rvals"/.Residuleh1DataStruc2D);
\[Theta]vals = ("\[Theta]vals"/.Residuleh1DataStruc2D);

rgrid =  ToPackedArray[Table[rvals[[i]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];
\[Theta]grid =  ToPackedArray[Table[\[Theta]vals[[j]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];

Inheritedfields =  KeyDrop[Residuleh1DataStruc2D,{"Data2D","h1Data1D","lmax","mmax"}];

dataoutput =  ToDataRegion[   tmp, {rvals[[1]],\[Theta]vals[[1]]},{Abs[rvals[[1]] -rvals[[2]]] ,Abs[\[Theta]vals[[2]] -\[Theta]vals[[1]]] } ];
Return[Join[<|"Data2D"->dataoutput,"rGrid"->rgrid,"\[Theta]Grid"->\[Theta]grid ,"a"->a0,"rOrbit"->r0|>,Inheritedfields]];
];


h1$Punc$dataStruct[{ai_,bi_},mi_,Residuleh1DataStruc2D_,{a0_,r0_}]:=Module[{rpRIndex,rp\[Theta]Index,dataoutputWindow,dataoutputNoWindow,dataoutputNoWindowRegion,dataoutputWindowRegion,dataoutputWindowFull,dataoutputNoWindowFull,tmp1,tmp2,\[Theta]gridWindow,lhsZeroGrid,WindowGrid,rgridWindow,rhsZeroGrid,puncturer1index,punctureRindex,\[Theta]grid,rgrid,Inheritedfields,tmp,rvals,\[Theta]vals,dataoutput,RWindowIndex},

rvals = ("rvals"/.Residuleh1DataStruc2D);
\[Theta]vals = ("\[Theta]vals"/.Residuleh1DataStruc2D);

rpRIndex=First@Ordering[Abs[rvals-r0],1];
rp\[Theta]Index=First@Ordering[Abs[\[Theta]vals-\[Pi]/2],1];

{puncturer1index,punctureRindex} = "windowRadialindex"/.Residuleh1DataStruc2D;

rgrid =  ToPackedArray[Table[rvals[[i]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];
\[Theta]grid =  ToPackedArray[Table[\[Theta]vals[[j]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];

lhsZeroGrid =0* rgrid[[;;puncturer1index-1]];
rhsZeroGrid = 0*rgrid[[punctureRindex+1;;]];
\[Theta]gridWindow =\[Theta]grid[[puncturer1index;;punctureRindex]];
rgridWindow =rgrid[[puncturer1index;;punctureRindex]];


tmp = Quiet[h1$Punc$Tetrad$dataInteral[{ai,bi},mi,{rgrid,\[Theta]grid},{a0,r0}]];
dataoutputWindow =Quiet[h1$Punc$Tetrad$dataInteral[{ai,bi},mi,{rgridWindow,\[Theta]gridWindow},{a0,r0}]];
dataoutputWindowFull = Join[lhsZeroGrid,dataoutputWindow,rhsZeroGrid];
dataoutputWindowFull[[rpRIndex,rp\[Theta]Index]] = 0;
Inheritedfields =  KeyDrop[Residuleh1DataStruc2D,{"Data2D","h1Data1D","lmax","mmax"}];
dataoutputWindowRegion = ToDataRegion[ dataoutputWindowFull, {rvals[[1]],\[Theta]vals[[1]]},{Abs[rvals[[1]] -rvals[[2]]] ,Abs[\[Theta]vals[[2]] -\[Theta]vals[[1]]] } ];

Return[Join[<|"Data2D"->dataoutputWindowRegion,"rGrid"->rgrid,"\[Theta]Grid"->\[Theta]grid,"a"->a0,"rOrbit"->r0 |>,Inheritedfields]];
];


Punc\[Theta]Window[\[Theta]_,x0_,w_,q_,s_]:=Piecewise[{{0,\[Theta]<=x0},{1/2+Tanh[(s*(-(q^2*Cot[(Pi*(-x0+\[Theta]))/(2*w)])+Tan[(Pi*(-x0+\[Theta]))/(2*w)]))/Pi]/2,x0<\[Theta]<w+x0},{1,Pi-w-x0>=\[Theta]>=w+x0},{1/2+Tanh[(s*(-(q^2*Cot[(Pi*(Pi-x0-\[Theta]))/(2*w)])+Tan[(Pi*(Pi-x0-\[Theta]))/(2*w)]))/Pi]/2,Pi-x0>\[Theta]>Pi-w-x0}},0];
WindowFunc = Function[{Global`\[Theta]},Punc\[Theta]Window[Global`\[Theta],0.01`50, 6 \[Pi]/16,1,3],Listable];


h1$Window$Punc$dataStruct[{ai_,bi_},mi_,Residuleh1DataStruc2D_,{a0_,r0_}]:=Module[{rp\[Theta]Index,rpRIndex,dataoutputWindow,dataoutputNoWindow,dataoutputNoWindowRegion,dataoutputWindowRegion,dataoutputWindowFull,dataoutputNoWindowFull,tmp1,tmp2,\[Theta]gridWindow,lhsZeroGrid,WindowGrid,rgridWindow,rhsZeroGrid,puncturer1index,punctureRindex,\[Theta]grid,rgrid,Inheritedfields,tmp,rvals,\[Theta]vals,dataoutput,RWindowIndex},

rvals = ("rvals"/.Residuleh1DataStruc2D);
\[Theta]vals = ("\[Theta]vals"/.Residuleh1DataStruc2D);

rpRIndex=First@Ordering[Abs[rvals-r0],1];
rp\[Theta]Index=First@Ordering[Abs[\[Theta]vals-\[Pi]/2],1];

{puncturer1index,punctureRindex} = "windowRadialindex"/.Residuleh1DataStruc2D;

rgrid =  ToPackedArray[Table[rvals[[i]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];
\[Theta]grid =  ToPackedArray[Table[\[Theta]vals[[j]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];

lhsZeroGrid =0* rgrid[[;;puncturer1index-1]];
rhsZeroGrid = 0*rgrid[[punctureRindex+1;;]];
\[Theta]gridWindow =\[Theta]grid[[puncturer1index;;punctureRindex]];
rgridWindow =rgrid[[puncturer1index;;punctureRindex]];
WindowGrid = Punc\[Theta]Window[\[Theta]grid];

tmp = Quiet[h1$Punc$Tetrad$dataInteral[{ai,bi},mi,{rgrid,\[Theta]grid},{a0,r0}]];
dataoutputWindow =WindowFunc[\[Theta]gridWindow] Quiet[h1$Punc$Tetrad$dataInteral[{ai,bi},mi,{rgridWindow,\[Theta]gridWindow},{a0,r0}]];
dataoutputWindowFull = Join[lhsZeroGrid,dataoutputWindow,rhsZeroGrid];
dataoutputWindowFull[[rpRIndex,rp\[Theta]Index]] = 0;
Inheritedfields =  KeyDrop[Residuleh1DataStruc2D,{"Data2D","h1Data1D","lmax","mmax"}];
dataoutputWindowRegion = ToDataRegion[ dataoutputWindowFull, {rvals[[1]],\[Theta]vals[[1]]},{Abs[rvals[[1]] -rvals[[2]]] ,Abs[\[Theta]vals[[2]] -\[Theta]vals[[1]]] } ];

Return[Join[<|"Data2D"->dataoutputWindowRegion,"rGrid"->rgrid,"\[Theta]Grid"->\[Theta]grid,"a"->a0,"rOrbit"->r0 |>,Inheritedfields]];
];


h1$Ret$Punc$dataStruct[{ai_,bi_},mi_,Residuleh1DataStruc2D_,{a0_,r0_}]:=Module[{Retdata,Inheritedfields,h1windowpuncdat,h1windowresdat},

h1windowresdat = h1$Res$dataStruct[{ai,bi},mi,Residuleh1DataStruc2D,{a0,r0}];
h1windowpuncdat = h1$Window$Punc$dataStruct[{ai,bi},mi,Residuleh1DataStruc2D,{a0,r0}];


Inheritedfields =  KeyDrop[Residuleh1DataStruc2D,{"Data2D"}];
Retdata  =  ("Data2D"/.h1windowpuncdat) + ("Data2D"/.h1windowresdat);
Return[Join[<|"Data2D"->Retdata |>,Inheritedfields]];
];


h1$Res$mmode$dataStruct[{ai_,bi_},mi_,Residuleh1DataStruc2D_,{a0_,r0_}]:=Module[{h1puncdat,Resdata,Inheritedfields,h1windowpuncdat,h1windowresdat},

h1windowresdat = h1$Res$dataStruct[{ai,bi},mi,Residuleh1DataStruc2D,{a0,r0}];
h1windowpuncdat = h1$Window$Punc$dataStruct[{ai,bi},mi,Residuleh1DataStruc2D,{a0,r0}];
h1puncdat = h1$Punc$dataStruct[{ai,bi},mi,Residuleh1DataStruc2D,{a0,r0}];

Inheritedfields =  KeyDrop[Residuleh1DataStruc2D,{"Data2D"}];
Resdata  = ("Data2D"/.h1windowresdat)+ ("Data2D"/.h1windowpuncdat) - ("Data2D"/.h1puncdat);
Return[Join[<|"Data2D"->Resdata |>,Inheritedfields]];
];


h1$Punc$mmode$dataStruct[{ai_,bi_},mi_,Residuleh1DataStruc2D_,{a0_,r0_}]:=Module[{h1puncdat,Retdata,Inheritedfields,h1windowpuncdat,h1windowresdat},

h1puncdat = h1$Punc$dataStruct[{ai,bi},mi,Residuleh1DataStruc2D,{a0,r0}];

Inheritedfields =  KeyDrop[Residuleh1DataStruc2D,{"Data2D"}];
Retdata  =  ("Data2D"/.h1puncdat);
Return[Join[<|"Data2D"->Retdata |>,Inheritedfields]];
];


(*DataStructAllCompoenets*)


TensorConstruct2DData[mi_,SingleCompStructFunc_,Residuleh1DataStruc2D_,{a0_,r0_}]:=Module[{FullTetradGrid,dataoutputNoWindow,dataoutputNoWindowRegion,dataoutputWindowRegion,dataoutputWindowFull,dataoutputNoWindowFull,tmp1,tmp2,\[Theta]gridWindow,lhsZeroGrid,WindowGrid,rgridWindow,rhsZeroGrid,puncturer1index,punctureRindex,\[Theta]grid,rgrid,Inheritedfields,tmp,rvals,\[Theta]vals,dataoutput,RWindowIndex},

rvals = ("rvals"/.Residuleh1DataStruc2D);
\[Theta]vals = ("\[Theta]vals"/.Residuleh1DataStruc2D);

Inheritedfields =  KeyDrop[Residuleh1DataStruc2D,{"Data2D","h1ResData2D","h1ResDat","h1PuncDat","lmax","mmax"}];

rgrid =  ToPackedArray[Table[rvals[[i]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];
\[Theta]grid =  ToPackedArray[Table[\[Theta]vals[[j]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];

FullTetradGrid = Table["Data2D"/.SingleCompStructFunc[{ai,bi},mi,Residuleh1DataStruc2D,{a0,r0}],{ai,4},{bi,4}];

Return[Join[<|"Data2D"->FullTetradGrid,"rGrid"->rgrid,"\[Theta]Grid"->\[Theta]grid,"a"->a0,"rOrbit"->r0,"mVal"->mi |>,Inheritedfields]];
];


TensorConstruct2DDataSave[mi_,SingleCompStructFunc_,Residuleh1DataStruc2D_,{a0_,r0_}]:=Module[{FullTetradGrid,dataoutputNoWindow,dataoutputNoWindowRegion,dataoutputWindowRegion,dataoutputWindowFull,dataoutputNoWindowFull,tmp1,tmp2,\[Theta]gridWindow,lhsZeroGrid,WindowGrid,rgridWindow,rhsZeroGrid,puncturer1index,punctureRindex,\[Theta]grid,rgrid,Inheritedfields,tmp,rvals,\[Theta]vals,dataoutput,RWindowIndex},

rvals = ("rvals"/.Residuleh1DataStruc2D);
\[Theta]vals = ("\[Theta]vals"/.Residuleh1DataStruc2D);

Inheritedfields =  KeyDrop[Residuleh1DataStruc2D,{"Data2D","h1ResData2D","h1ResDat","h1PuncDat","lmax","mmax"}];

rgrid =  ToPackedArray[Table[rvals[[i]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];
\[Theta]grid =  ToPackedArray[Table[\[Theta]vals[[j]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]];

FullTetradGrid = Table["Data2D"/.SingleCompStructFunc[{ai,bi},mi,Residuleh1DataStruc2D,{a0,r0}],{ai,4},{bi,4}];

Return[Join[<|"Data2D"->  Table[ToListOfData[FullTetradGrid[[i,j]]],{i,4},{j,4}],"rGrid"->rgrid,"\[Theta]Grid"->\[Theta]grid,"a"->a0,"rOrbit"->r0,"mVal"->mi |>,Inheritedfields]];
];


(*Convert to BL*)


DataStructKinnerslyTetradToBLComp[DataStruc2D_,{a0_,r0_}]:=Module[{\[CapitalDelta],\[CapitalSigma],Inheritedfields,rgridtmp,\[Theta]gridtmp,rvals,\[Theta]vals,rgrid,\[Theta]grid,ldBLGridK,ndBLGridK,mdBLGridK,mbdBLGridK,ZTildaGrid,D2GBL,Retarded$tetrad$data},
rvals = ("rvals"/.DataStruc2D);
\[Theta]vals = ("\[Theta]vals"/.DataStruc2D);

Retarded$tetrad$data = ("Data2D"/.DataStruc2D);

rgridtmp =  N[ToPackedArray[Table[rvals[[i]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]]];
\[Theta]gridtmp =  N[ToPackedArray[Table[\[Theta]vals[[j]], {i,Length[rvals]},{j,Length[\[Theta]vals]}]]];

rgrid = ToDataRegion[ rgridtmp, {rvals[[1]],\[Theta]vals[[1]]},{Abs[rvals[[1]] -rvals[[2]]] ,Abs[\[Theta]vals[[2]] -\[Theta]vals[[1]]] } ];
\[Theta]grid = ToDataRegion[ \[Theta]gridtmp, {rvals[[1]],\[Theta]vals[[1]]},{Abs[rvals[[1]] -rvals[[2]]] ,Abs[\[Theta]vals[[2]] -\[Theta]vals[[1]]] } ];

\[CapitalDelta] = r^2 - 2 r + a0^2;
\[CapitalSigma] = r^2 + a0^2 Cos[\[Theta]]^2;

ldBLK =  {-Ones,\[CapitalSigma]/\[CapitalDelta],Zero,a0 Sin[\[Theta]]^2};
ndBLK = 1/2 {-\[CapitalDelta]/\[CapitalSigma],-Ones,Zero,a Sin[\[Theta]]^2 \[CapitalDelta]/\[CapitalSigma]};
mdBLK =1/(Sqrt[2] ( r+ I a0 Cos[\[Theta]])) {- I a0 Sin[\[Theta]],Zero,\[CapitalSigma],I (a0^2+r^2) Sin[\[Theta]]};
mbdBLK =1/(Sqrt[2] ( r-I a0 Cos[\[Theta]])) {I a0 Sin[\[Theta]],Zero,\[CapitalSigma],-I (a0^2+r^2) Sin[\[Theta]]};

ldBLGridK=ldBLK/.{M->1,a->a0,\[Theta]->ToPackedArray[\[Theta]grid], r->ToPackedArray[rgrid],Zero->0*ToPackedArray[rgrid],Ones->ToPackedArray[rgrid]/ToPackedArray[rgrid]};
ndBLGridK=ndBLK/.{M->1,a->a0,\[Theta]->ToPackedArray[\[Theta]grid], r->ToPackedArray[rgrid],Zero->0*ToPackedArray[rgrid],Ones->ToPackedArray[rgrid]/ToPackedArray[rgrid]};
mdBLGridK= mdBLK/.{M->1,a->a0,\[Theta]->ToPackedArray[\[Theta]grid], r->ToPackedArray[rgrid],Zero->0*ToPackedArray[rgrid],Ones->ToPackedArray[rgrid]/ToPackedArray[rgrid]};
mbdBLGridK=mbdBLK/.{M->1,a->a0,\[Theta]->ToPackedArray[\[Theta]grid], r->ToPackedArray[rgrid],Zero->0*ToPackedArray[rgrid],Ones->ToPackedArray[rgrid]/ToPackedArray[rgrid]};


ZTildaGrid = {-ndBLGridK,-ldBLGridK,mbdBLGridK,mdBLGridK};

D2GBL = Monitor[Table[ Sum[ZTildaGrid[[k,\[Mu]]]ZTildaGrid[[l,\[Nu]]] Retarded$tetrad$data[[k,l]],{k,1,4},{l,1,4}],{\[Mu],1,4},{\[Nu],1,4}],{\[Mu],\[Nu]}];
Inheritedfields =  KeyDrop[DataStruc2D,{"Data2D","lmax","mmax"}];
Return[Join[<|"Data2D"->D2GBL,"rGrid"->rgrid,"\[Theta]Grid"->\[Theta]grid,"a"->a0,"rOrbit"->r0 |>,Inheritedfields]];
];


End[];


EndPackage[];
