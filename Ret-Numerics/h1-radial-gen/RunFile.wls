#!/usr/bin/env wolframscript

(* ::Subsubsection:: *)
(*Prelim*)

(*Add package directories to path*)
AppendTo[$Path, "/lustre/hpc/astro/dyson/.Wolfram/Applications/ModGEMS"];
AppendTo[$Path, DirectoryName[$InputFileName]];

Needs["MetricReconstructRadiative`"]


(*Read (a, r0, savepath) from command-line arguments: wolframscript RunFile.wls -- <a> <r0> <savepath>*)
{a, r0} = ToExpression /@ $ScriptCommandLine[[2;;3]];
savepath = $ScriptCommandLine[[4]];
Print["Job started: a = ", a, ", r0 = ", r0, ", savepath = ", savepath];

(*Fixed parameters*)
mmax = 25;
tmax = 250.0;
n = 10;
angres = 4;
twr = 16;
twq = 8;
nu = 0.8000;
puncord = 4;
gvopt = 2;
glgopt = 0;
icopt = 0;
icampl = 0.0;
lmax = 30;
lplot = 25;
nterms = 8;
inford = 7;
horord = 6;
rinf = 4000.0;
xhor = 0.0001;
rmax = 500.0;
accgoal = 30;
kapord = 8;
prec = 120;
rstmin = -15.0;
rstmax = 500.0;
dformat = "Real64";


(* ::Subsubsection:: *)
(*Path Directory*)


parituclardirectory = (savepath <> "/data" <> ToString[Round[1000*a]] <> "/data" <> ToString[Round[1000*a]] <> "-" <> ToString[Round[r0*10, 1]]);

configformat[] := Module[{text, filetag, filename},
  text = "a\t" <> ToString[a] <> "\nr0\t" <> ToString[r0] <> "\nmmax\t" <> ToString[mmax] <> "\ntmax\t" <> ToString[tmax] <> "\nn\t" <> ToString[n] <> "\nangres\t" <> ToString[angres] <> "\ntwr\t" <> ToString[twr] <> "\ntwq\t" <> ToString[twq] <> "\nnu\t" <> ToString[nu] <> "\npuncord\t" <> ToString[puncord] <> "\ngv_opt\t" <> ToString[gvopt] <> "\nglg_opt\t" <> ToString[glgopt] <> "\nic_opt\t" <> ToString[icopt] <> "\nic_ampl\t" <> ToString[icampl] <> "\nlmax\t" <> ToString[lmax] <> "\nlplot\t" <> ToString[lplot] <> "\nnterms\t" <> ToString[nterms] <> "\ninford\t" <> ToString[inford] <> "\nhorord\t" <> ToString[horord] <> "\nrinf\t" <> ToString[rinf] <> "\nxhor\t" <> ToString[xhor] <> "\nrmax\t" <> ToString[rmax] <> "\naccgoal\t" <> ToString[accgoal] <> "\nkapord\t" <> ToString[kapord] <> "\nprec\t" <> ToString[prec] <> "\nrstmin\t" <> ToString[rstmin] <> "\nrstmax\t" <> ToString[rstmax] <> "\ndformat\t" <> ToString[dformat];
  filetag = ToString[Round[1000*a]] <> "-" <> ToString[Round[r0*10, 1]];
  filename = "config/config" <> filetag <> ".txt";
  Export[parituclardirectory <> "/" <> filename, text];
  filetag];


(*Check if the directory exists before creating it*)
If[!DirectoryQ[parituclardirectory], CreateDirectory[parituclardirectory]];
If[!DirectoryQ[parituclardirectory <> "/config"], CreateDirectory[parituclardirectory <> "/config"]];
If[!DirectoryQ[parituclardirectory <> "/data"], CreateDirectory[parituclardirectory <> "/data"]];

(*Create a single meta config file; derive per-m configTags for MetricReconstructRadiative*)
configTag = configformat[];
configTags = Table[configTag <> "-00" <> ToString[i], {i, 0, mmax}];


(* ::Subsubsection:: *)
(*Run*)


\[Omega]0 = 2 KerrAzimuthalFrequency[a, r0, r0, 0];
OD = OrbitalData[{a, r0, 0, 0, 0},
  "OrbitPrecision" -> prec,
  "OrbitWorkingPrecision" -> prec
];
mvalslist = Table[mi, {mi, mmax}];


xhorP = SetPrecision[xhor, prec];
rprmvals = {rp -> 1 + Sqrt[1 - a^2], rm -> 1 - Sqrt[1 - a^2]};
rhor = (rp + xhorP) /. rprmvals;
rmin = rhor + 0.1`100;


{rgridL, rgridR} = MetricReconstructRadiative`Private`BuildGrid[a, {r0, rmin, rmax}, {WorkingPrecision -> prec, "nterms" -> nterms, "kapord" -> kapord, "rinf" -> rinf, "rmax" -> rmax, "xhor" -> xhorP, "inford" -> inford, "horord" -> horord, AccuracyGoal -> accgoal, "rgrid" -> 1, "rstmin" -> rstmin, "rstmax" -> rstmax, "nres" -> n, "angres" -> angres}];


Print["Now we Calculate the Radiative Solutions"];
Table[
  MetricReconstructRadiative[OD, {lmax, j}, {rgridL, rgridR}, configTags[[j + 1]], savepath,
    {
      WorkingPrecision -> prec,
      "rinf" -> rinf,
      "xhor" -> xhorP,
      AccuracyGoal -> accgoal,
      "dformat" -> dformat,
      Order->50
    }
  ],
  {j, mvalslist}
];
