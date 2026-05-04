#!/usr/bin/env wolframscript

(* ::Subsubsection:: *)
(*Prelim*)

Needs["MetricReconstructRadiative`"]
ParallelNeeds["MetricReconstructRadiative`"]


(*Save path*)
savepath = "/Users/conordyson/Documents/Research/Open/Gravitational-SF/h1Solvers/KerrLorenzMSF/Radiative-Cluster/Numerics";
(*Parameters*)
a = 0.6;
r0 = 8;
mmax = 10;
tmax = 250.0;
n = 4;
angres = 4;
twr = 16;
twq = 8;
nu = 0.8000;
puncord = 4;
gvopt = 2;
glgopt = 0;
icopt = 0;
icampl = 0.0;
lmax = 5;
lplot = 6;
nterms = 8;
inford = 7;
horord = 6;
rinf = 4000.0;
xhor = 0.0001;
rmax = 60.0;
accgoal = 24;
kapord = 8;
prec = 60;
rstmin = -15.0;
rstmax = 60.0;
dformat = "Real64";


(* ::Subsubsection:: *)
(*Path Directory*)


(*Create the configuration files*)
parituclardirectory = (savepath <> "/data" <> ToString[Round[1000*a]] <> "/data" <> ToString[Round[1000*a]] <> "-" <> ToString[Round[r0*10, 1]]);

configformat[m_] := Module[{text, filename, filetag},
  text = "a\t" <> ToString[a] <> "\nr0\t" <> ToString[r0] <> "\nm\t" <> ToString[m] <> "\ntmax\t" <> ToString[tmax] <> "\nn\t" <> ToString[n] <> "\nangres\t" <> ToString[angres] <> "\ntwr\t" <> ToString[twr] <> "\ntwq\t" <> ToString[twq] <> "\nnu\t" <> ToString[nu] <> "\npuncord\t" <> ToString[puncord] <> "\ngv_opt\t" <> ToString[gvopt] <> "\nglg_opt\t" <> ToString[glgopt] <> "\nic_opt\t" <> ToString[icopt] <> "\nic_ampl\t" <> ToString[icampl] <> "\nlmax\t" <> ToString[lmax] <> "\nlplot\t" <> ToString[lplot] <> "\nnterms\t" <> ToString[nterms] <> "\ninford\t" <> ToString[inford] <> "\nhorord\t" <> ToString[horord] <> "\nrinf\t" <> ToString[rinf] <> "\nxhor\t" <> ToString[xhor] <> "\nrmax\t" <> ToString[rmax] <> "\naccgoal\t" <> ToString[accgoal] <> "\nkapord\t" <> ToString[kapord] <> "\nprec\t" <> ToString[prec] <> "\nrstmin\t" <> ToString[rstmin] <> "\nrstmax\t" <> ToString[rstmax] <> "\ndformat\t" <> ToString[dformat];
  filetag = ToString[Round[1000*a]] <> "-" <> ToString[Round[r0*10, 1]] <> "-00" <> ToString[m];
  filename = "config/config" <> filetag <> ".txt";
  Export[parituclardirectory <> "/" <> filename, text];
  {filetag, filename}];


(*Check if the directory exists before creating it*)
If[!DirectoryQ[parituclardirectory], CreateDirectory[parituclardirectory]];
If[!DirectoryQ[parituclardirectory <> "/config"], CreateDirectory[parituclardirectory <> "/config"]];
If[!DirectoryQ[parituclardirectory <> "/data"], CreateDirectory[parituclardirectory <> "/data"]];
If[!DirectoryQ[parituclardirectory <> "/jumpchecks"], CreateDirectory[parituclardirectory <> "/jumpchecks"]];

(*Create the configuration files*)
Files = Table[configformat[i], {i, 0, mmax}];

(*Clear all variables except the one you want to keep*)
configTags = Files[[;; , 1]];
configFiles = Files[[;; , 2]];


(* ::Subsubsection:: *)
(*Run*)


\[Omega]0 = 2 KerrAzimuthalFrequency[a, r0, r0, 0];
OD = OrbitalData[{a, r0, 0, 0, 0},
  "OrbitPrecision" -> prec,
  "OrbitWorkingPrecision" -> prec
];
mvalslist = Table[mi, {mi, mmax}];


xhor = SetPrecision[xhor, prec];
rprmvals = {rp -> 1 + Sqrt[1 - a^2], rm -> 1 - Sqrt[1 - a^2]};
rhor = (rp + xhor) /. rprmvals;
rmin = rhor + SetPrecision[10^-1, prec];


{rgridL, rgridR} = MetricReconstructRadiative`Private`BuildGrid[a, {r0, rmin, rmax}, {WorkingPrecision -> prec, "nterms" -> nterms, "kapord" -> kapord, "rinf" -> rinf, "rmax" -> rmax, "xhor" -> xhor, "inford" -> inford, "horord" -> horord, AccuracyGoal -> accgoal, "rgrid" -> 1, "rstmin" -> rstmin, "rstmax" -> rstmax, "nres" -> n, "angres" -> angres}];


DistributeDefinitions[OD, lmax, rgridL, rgridR, configTags, savepath, prec, rinf, xhor, accgoal, dformat];

Print["Now we Calculate the Radiative Solutions in Parallel"];
ParallelTable[
  MetricReconstructRadiative[OD, {lmax, j}, {rgridL, rgridR}, configTags[[j + 1]], savepath,
    {
      WorkingPrecision -> prec,
      "rinf" -> rinf,
      "xhor" -> xhor,
      AccuracyGoal -> accgoal,
      "dformat" -> dformat
    }
  ],
  {j, mvalslist}
];
