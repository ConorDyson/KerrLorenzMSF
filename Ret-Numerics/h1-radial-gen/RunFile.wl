(* ::Package:: *)

(* ::Chapter:: *)
(*Test*)


(* ::Subsubsection:: *)
(*Prelim*)

(*
SetDirectory[NotebookDirectory[]]; *)
Needs["MetricReconstructRadiative`"]
(* ParallelNeeds["MetricReconstructRadiative`"] *)


(*Save path*)
savepath = "/lustre/hpc/astro/dyson/KerrCircMSFMaarten";
(*Parameters*)
a=0.99;
r0=8;
mmax=30;
tmax=250.0;
n=10;
angres=4;
twr=16;
twq=8;
nu=0.8000;
puncord=4;
gvopt=2;
glgopt=0;
icopt=0;
icampl=0.0;
lmax=50;
lplot=45;
nterms=10;
inford=7;
horord=6;
rinf=4000.0;
xhor=0.0001;
rmax=30.0;
accgoal=20;
kapord=8;
prec=100;
rstmin=-15.0;
rstmax=60.0;
dformat="Real64";


(* ::Input:: *)
(**)


(* ::Subsubsection:: *)
(*Path Directory*)


(*Create the configuration files*)
parituclardirectory=(savepath<>"/data"<>ToString[Round[1000*a]]<>"/data"<>ToString[Round[1000*a]]<>"-"<>ToString[Round[r0*10,1]]);

configformat[m_]:=Module[{text,filename,filetag},text="a\t"<>ToString[a]<>"\nr0\t"<>ToString[r0]<>"\nm\t"<>ToString[m]<>"\ntmax\t"<>ToString[tmax]<>"\nn\t"<>ToString[n]<>"\nangres\t"<>ToString[angres]<>"\ntwr\t"<>ToString[twr]<>"\ntwq\t"<>ToString[twq]<>"\nnu\t"<>ToString[nu]<>"\npuncord\t"<>ToString[puncord]<>"\ngv_opt\t"<>ToString[gvopt]<>"\nglg_opt\t"<>ToString[glgopt]<>"\nic_opt\t"<>ToString[icopt]<>"\nic_ampl\t"<>ToString[icampl]<>"\nlmax\t"<>ToString[lmax]<>"\nlplot\t"<>ToString[lplot]<>"\nnterms\t"<>ToString[nterms]<>"\ninford\t"<>ToString[inford]<>"\nhorord\t"<>ToString[horord]<>"\nrinf\t"<>ToString[rinf]<>"\nxhor\t"<>ToString[xhor]<>"\nrmax\t"<>ToString[rmax]<>"\naccgoal\t"<>ToString[accgoal]<>"\nkapord\t"<>ToString[kapord]<>"\nprec\t"<>ToString[prec]<>"\nrstmin\t"<>ToString[rstmin]<>"\nrstmax\t"<>ToString[rstmax]<>"\ndformat\t"<>ToString[dformat];
filetag=ToString[Round[1000*a]]<>"-"<>ToString[Round[r0*10,1]]<>"-00"<>ToString[m];
filename="config/config"<>filetag<>".txt";
Export[parituclardirectory<>"/"<>filename,text];
{filetag,filename}];


(*Check if the directory exists before creating it*)
If[!DirectoryQ[parituclardirectory],CreateDirectory[parituclardirectory]];
If[!DirectoryQ[parituclardirectory<>"/config"],CreateDirectory[parituclardirectory<>"/config"]];
If[!DirectoryQ[parituclardirectory<>"/data"],CreateDirectory[parituclardirectory<>"/data"]];
If[!DirectoryQ[parituclardirectory<>"/jumpchecks"],CreateDirectory[parituclardirectory<>"/jumpchecks"]];

(*Create the configuration files*)
Files=Table[configformat[i],{i,0,mmax}];

(*Clear all variables except the one you want to keep*)
configTags=Files[[;;,1]];
configFiles=Files[[;;,2]];


(* ::Subsubsection:: *)
(*Run*)


\[Omega]0=2KerrAzimuthalFrequency[a,r0,r0,0];
OD=OrbitalData[{a,r0,0,0,0},
"OrbitPrecision"->prec,
"OrbitWorkingPrecision"->prec
];
mvalslist = Table[mi,{mi,mmax}];


xhor=SetPrecision[xhor,prec];
rprmvals={rp->1+Sqrt[1-a^2],rm->1-Sqrt[1-a^2]};
rmin=3;


{rgridL,rgridR}=MetricReconstructRadiative`Private`BuildGrid[a,{r0,rmin,rmax},{WorkingPrecision->prec,"nterms"-> nterms, (* Number of terms in the spherical expansion. *)"kapord"->kapord, (* \[Kappa]ord is the maximum order of series expansion of \[Kappa] in spheroidal harmonics. *)"rinf"->rinf, (* rinf is the radius at which the series solutions for \[Kappa]_up should set the initial conditions for the integrator. *)"rmax"->rmax,"xhor"->xhor,"inford"->inford, (* The order of the expansion at infinity for the UP function. *)"horord"->horord, (* The order of the expansion at the horizon for the IN function. *)AccuracyGoal->accgoal,"rgrid"->1, (* Use a linearly-spaced grid in the variable: 0 = rstar , 1 = r. *)"rstmin"->rstmin,"rstmax"->rstmax,"nres"->n,(* Resolution in the r* direction:  dr* = M / n  (or dr = M / n). *)"angres"->angres (* Resolution in the \[Theta] direction: number of points = nres * qres. *)}];


Monitor[Table[MetricReconstructRadiative[OD,{lmax,j},{rgridL,rgridR},configTags[[j+1]],savepath,
{
WorkingPrecision->prec,
"rinf"->rinf,
"xhor"->xhor,
AccuracyGoal->accgoal,
"dformat"->dformat,
Order->50
}
],{j,mvalslist}];,{j}]


