(* ::Package:: *)

savepath = "/Users/conordyson/Documents/Research/Open/Gravitational-SF/h1Solvers/KerrLorenzMSF/Radiative-Cluster/Numerics";

(* Directory tags — from the cluster run config (a=0.6, r0=8) *)
aDir  = 0.6;
r0Dir = 8.0;

(* Physics parameters in the HDF5 filenames *)
a    = 0.0;
rp   = 8.0;
lmax = 20;
mmax = 4;

component = {"h_l+l+", "h_l-l-", "h_m+m+", "h_m-m-",
             "rho_h_l+m+", "rhob_h_l+m-", "rhob_h_l-m+", "rho_h_l-m-",
             "sigma_delta_h_l+l-", "h"};

numStr[x_] := Module[{s = ToString[x]},
  Which[
    StringMatchQ[s, ___ ~~ "."], s <> "0",
    ! StringContainsQ[s, "."],   s <> ".0",
    True,                         s
  ]
]

dataPath =
  savepath <>
  "/data" <> ToString[Round[1000 aDir]] <>
  "/data" <> ToString[Round[1000 aDir]] <> "-" <> ToString[Round[10 r0Dir]] <>
  "/data/";

h5file[mi_] :=
  dataPath <> "h1_a" <> numStr[a] <> "_rp" <> numStr[rp] <>
  "_l" <> ToString[lmax] <> "_m" <> ToString[mi] <> ".h5";

extractComplex[data_] := Map[("Re" /. #) + I ("Im" /. #) &, data, {2}]

LoadMode[mi_, index_] := Module[
  {miLoad, file, grp, dataIn, dataUp, rIn, rUp, datavals, datavalsNew, rvals},

  miLoad = If[mi == 1 && MemberQ[{8, 9}, index], 2, mi];
  file   = h5file[miLoad];
  grp    = "/m_" <> ToString[miLoad] <> "/";

  dataIn  = extractComplex @ Import[file, {"Datasets", grp <> "In/" <> component[[index]]}];
  dataUp  = extractComplex @ Import[file, {"Datasets", grp <> "Up/" <> component[[index]]}];
  rIn     = Import[file, {"Datasets", grp <> "r_in"}];
  rUp     = Import[file, {"Datasets", grp <> "r_up"}];

  datavals    = Join[dataIn[[;; , ;; -2]], dataUp, 2];
  datavalsNew = PadLeft[datavals, {lmax + 1, Dimensions[datavals][[2]]}];
  rvals       = Join[rIn[[;; -2]], rUp];

  {rvals, datavalsNew}
]

Print["Loading from: ", dataPath];

h1Data = Table[LoadMode[mi, index], {mi, mmax}, {index, 10}];

data = Association[{
  "h1Data1D" -> Transpose[Transpose[h1Data[[;; , ;; , 2]], 2 <-> 3], 1 <-> 2],
  "rvals"    -> h1Data[[1, 1, 1]],
  "lmax"     -> lmax,
  "mmax"     -> mmax,
  "a0"       -> a,
  "rOrbit"   -> rp,
  "DataType" -> "1RadialData"
}];

Print["Done. rvals length: ",       Length[data["rvals"]]];
Print["h1Data1D dimensions: ", Dimensions[data["h1Data1D"]]];
