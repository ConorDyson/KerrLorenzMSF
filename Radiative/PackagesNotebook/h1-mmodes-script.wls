#!/usr/bin/env wolframscript
(* ::Package:: *)

Needs["MetricReconstructRadiative`"]
Needs["MetricReconstructRadiative`"]
Needs["Teukolsky`"]
Needs["SpinWeightedSpheroidalHarmonics`"]
Needs["PTools`"];
ParallelNeeds["MetricReconstructRadiative`"];
ParallelNeeds["Teukolsky`"]
ParallelNeeds["SpinWeightedSpheroidalHarmonics`"]
ParallelNeeds["PTools`"];
<<MetricReconstructStatic`

DistributeDefinitions["Global`"];

(* Comment out parallel kernel launch to avoid issues *)
(* Try to read the allocated cores from the environment variable *)
nCoresStr = $SystemEnvironment["SLURM_CPUS_PER_TASK"];
nCores = If[ StringQ[nCoresStr] && nCoresStr != "", ToExpression[nCoresStr], 4 ]; 
(* Fallback to 4 cores if the variable isn’t set *)

(* Ensure you don’t try to launch more kernels than your machine supports *)
LaunchKernels[Min[nCores, $ProcessorCount]];

Print["Command-line arguments: ", $CommandLine[[-1]]];

(* Function to read configuration file and extract parameters *)
readConfig[filename_] := Module[{config, params},
  config = Import[filename, "Table"];
  params = Association[Rule @@@ config];
  params
]

(* Read the configuration file *)
(* configParams = readConfig["ConfigFile.txt"]; *)
configParams = readConfig[$CommandLine[[-1]]];

(* Extract parameters *)
a = configParams["a"];
r0 = configParams["r0"];
mmax = configParams["mmax"];
tmax = configParams["tmax"];
n = configParams["n"];
angres = configParams["angres"];
twr = configParams["twr"];
twq = configParams["twq"];
nu = configParams["nu"];
puncord = configParams["puncord"];
gvopt = configParams["gv_opt"];
glgopt = configParams["glg_opt"];
icopt = configParams["ic_opt"];
icampl = configParams["ic_ampl"];
lmax = configParams["lmax"];
lplot = configParams["lplot"];
nterms = configParams["nterms"];
inford = configParams["inford"];
horord = configParams["horord"];
rinf = configParams["rinf"];
xhor = configParams["xhor"];
rmax = configParams["rmax"];
accgoal = configParams["accgoal"];
kapord = configParams["kapord"];
prec = configParams["prec"];
rstmin = configParams["rstmin"];
rstmax = configParams["rstmax"];
dformat = configParams["dformat"];
  
(* Create the configuration files *)

configformat[m_] := 
 Module[{text, filename, filetag}, 
  text = "a\t" <> ToString[a] <> "\nr0\t" <> ToString[r0] <> "\nm\t" <>
   ToString[m] <> "\ntmax\t" <> ToString[tmax] <> "\nn\t" <> 
  ToString[n] <> "\nangres\t" <> ToString[angres] <> "\ntwr\t" <> 
  ToString[twr] <> "\ntwq\t" <> ToString[twq] <> "\nnu\t" <> 
  ToString[nu] <> "\npuncord\t" <> ToString[puncord] <> 
  "\ngv_opt\t" <> ToString[gvopt] <> "\nglg_opt\t" <> ToString[glgopt] <> 
  "\nic_opt\t" <> ToString[icopt] <> "\nic_ampl\t" <> ToString[icampl] <> 
  "\nlmax\t" <> ToString[lmax] <> "\nlplot\t" <> ToString[lplot] <> 
  "\nnterms\t" <> ToString[nterms] <> "\ninford\t" <> ToString[inford] <> 
  "\nhorord\t" <> ToString[horord] <> "\nrinf\t" <> ToString[rinf] <> 
  "\nxhor\t" <> ToString[xhor] <> "\nrmax\t" <> ToString[rmax] <> 
  "\naccgoal\t" <> ToString[accgoal] <> "\nkapord\t" <> ToString[kapord] <> 
  "\nprec\t" <> ToString[prec] <> "\nrstmin\t" <> ToString[rstmin] <>
   "\nrstmax\t" <> ToString[rstmax] <> "\ndformat\t" <> 
  ToString[dformat];

  filetag =  ToString[Round[1000*a]] <> "-" <>  ToString[Round[r0*10, 1]] <> "-00" <> ToString[m];
  filename = "config/config" <> filetag <> ".txt";
  Export[parituclardirectory<>"/"<>filename, text];
 
  {filetag, filename}
];


primarypath = "../../data/h1_lorenz/july25_r0_6/";
(* Create the directory structure *)

parituclardirectory = (primarypath<>"/data"<>ToString[Round[1000*a]]<>"/data"<>ToString[Round[1000*a]] <> "-" <>  ToString[Round[r0*10, 1]]);

(* Create the directory structure *)


(* Check if the directory exists before creating it *)
If[!DirectoryQ[parituclardirectory],
  CreateDirectory[parituclardirectory]
];


If[!DirectoryQ[parituclardirectory <> "/config"],
  CreateDirectory[parituclardirectory <> "/config"]
];
If[!DirectoryQ[parituclardirectory <> "/data"],
  CreateDirectory[parituclardirectory <> "/data"]
];
If[!DirectoryQ[parituclardirectory <> "/jumpchecks"],
  CreateDirectory[parituclardirectory <> "/jumpchecks"]
];

 
(* Create the configuration files *)
Files = Table[configformat[i], {i, 0, mmax}];

(* Clear all variables except the one you want to keep *)

ClearAll[
  a, r0, tmax, n, angres, twr, twq, nu, puncord, gvopt, glgopt, 
  icopt, icampl, lmax, lplot, nterms, inford, horord, rinf, xhor, rmax, 
  accgoal, kapord, prec, rstmin, rstmax, dformat
];

configTags  = Files[[;;,1]];
configFiles = Files[[;;,2]];

(* Print["First we Calculate the Static Solution"];
Begin["MyPackage`Private`"];
Quiet[MetricReconstructStatic[configTags[[1]] , primarypath ]]; 
End[]; *)

Print["Now we Calculate the Radiative Solutions in Parallel"];
(* Table[ Quiet[MetricReconstructRadiative[configTags[[j+1]]]], {j,1, mmax}];  *)
ParallelTable[MetricReconstructRadiative[configTags[[j+1]] , primarypath], {j, mmax}, Method -> "FinestGrained"];  
 