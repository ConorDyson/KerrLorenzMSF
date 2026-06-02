(* ::Package:: *)

(* h1FunctionalityAsymp.wl

   Package for loading, transforming, and constructing the 2D m-mode h1
   perturbation tensor. Handles radial data in both Sam's binary format and
   Maarten's HDF5 format, assembles the (r, theta) grid, applies the
   Kinnersley-tetrad basis conversion, and optionally converts to
   Boyer-Lindquist components.
*)

Needs["KerrPunctures`"];

BeginPackage["h1FunctionalityAsymp`", {
  "SpinWeightedSpheroidalHarmonics`",
  "Developer`",
  "SimulationTools`",
  "KerrPunctures`",
  "NContinuedFraction`",
  "SpinWeightedSpheroidalHarmonicsFT`",
  "SWSHdecomp`"
}];

(* ------------------------------------------------------------------ *)
(* Public interface                                                     *)
(* ------------------------------------------------------------------ *)

RetSamDImportRetFunction::usage =
  "RetSamDImportRetFunction: alias for the internal single-m loader \
(RetSamDImportRetFunctionInteral).";

RetStaticSamDImportRetFunction::usage =
  "RetStaticSamDImportRetFunction[a, rp, mmax, lmax, filepath] loads the \
m=0 (static) retarded h1 data from Sam's binary format.";

RetRadiativeSamDImportRetFunction::usage =
  "RetRadiativeSamDImportRetFunction[a, rp, mmax, lmax, filepath] loads all \
radiative m-modes from Sam's binary format, applies the spheroidal-to-spherical \
trace transformation, and returns a 1D radial data struct.";

RetRadiativeMaartenImportRetFunction::usage =
  "RetRadiativeMaartenImportRetFunction[a, rp, mmax, lmax, filepath] loads \
radiative h1 data from Maarten's HDF5 format and returns a 1D radial data struct.";

JoininRadiativeandStaticPieces::usage =
  "JoininRadiativeandStaticPieces[a, rOrbit, mmax, lmax, RadStruct, StaticStruct] \
merges the radiative and static 1D data structs onto the radiative radial grid.";

SphericalHamronicsContructor::usage =
  "SphericalHamronicsContructor[lmax, mmax, thetavals] builds a table of \
spin-weighted spherical harmonics Y_s^{lm}(theta, 0) for s = -2,...,2, \
l = 0,...,lmax, m = 0,...,mmax.";

MakeResGrid::usage =
  "MakeResGrid[mmax, lmax, DataStruct1D, YsStruct] projects the 1D radial \
coefficients onto the 2D (r, theta) grid via a sum over l-modes.";

h1$tetrad$data::usage =
  "h1$tetrad$data: Kinnersley-tetrad h1 data on the 2D grid (internal use).";

h1$Ret$dataStruct::usage =
  "h1$Ret$dataStruct[{a,b}, m, DataStruct2D, {a0,r0}] returns the retarded \
h1 Kinnersley-tetrad component (a,b) for m-mode m as a DataRegion on the 2D grid.";

TensorConstruct2DData::usage =
  "TensorConstruct2DData[m, ComponentFunc, DataStruct2D, {a0,r0}] assembles \
all 16 Kinnersley-tetrad components into a 4x4 grid struct for m-mode m.";

TensorConstruct2DDataSave::usage =
  "TensorConstruct2DDataSave[m, ComponentFunc, DataStruct2D, {a0,r0}] is the \
same as TensorConstruct2DData but converts DataRegions to plain arrays for \
HDF5 export.";

DataStructKinnerslyTetradToBLComp::usage =
  "DataStructKinnerslyTetradToBLComp[DataStruct2D, {a0,r0}] contracts the \
Kinnersley tetrad components with the Boyer-Lindquist covectors to produce the \
coordinate-basis h_{mu nu} on the 2D grid.";


Begin["`Private`"];

Get[FindFile["KerrPunctures`"]];
Get[FindFile["SimulationTools`"]];


(* ==================================================================== *)
(* Section 1: Data import — Sam's binary format                         *)
(* ==================================================================== *)

(* Load a single m-mode from Sam's binary .bin format.
   Parses the config file, imports ingoing and outgoing radial data, and
   stitches them across the orbit.
   Returns Association["h1Dat" -> psis, "rvals" -> r, "lmax" -> lmax]. *)
RetSamDImportRetFunctionInteral[a_, rp_, m_, filepathin_] :=
  Module[{
    iConfig, filepath, configpath, datagpath, filen,
    s, s1, s2, configparams, GetParam, prec, iround,
    a0, r0, mm, lmax, lplot, dformat,
    tmp, rstarsL, rsL, rstarsR, rsR,
    fn, dataL, dimsL, GetIndexL, compL,
        dataR, dimsR, GetIndexR, compR,
    subrs, psis
  },

  iConfig    = ToString[Round[1000 a]] <> "-" <> ToString[Round[10 rp]] <> "-00" <> ToString[m];
  filepath   = filepathin <> "/a_" <> ToString[NumberForm[a, 2]] <>
               "_r0_" <> ToString[NumberForm[rp, {Infinity, 1}]] <> "/";
  configpath = filepath <> "config/";
  datagpath  = filepath <> "data/";

  (* Parse tab-delimited config file *)
  filen        = configpath <> "config" <> ToString[iConfig] <> ".txt";
  s            = Import[filen, "String"];
  s1           = StringSplit[s, "\n"];
  s2           = Map[StringSplit[#, "\t"] &, s1];
  configparams = Map[{#[[1]], ToExpression[#[[2]]]} &, s2];

  GetParam[key_] := Module[{ls, val},
    ls = Select[configparams, #[[1]] == key &];
    If[Length[ls] > 0,
      val = ls[[1, 2]];
      If[NumericQ[val] && !IntegerQ[val], val = SetPrecision[val, prec]],
      val = Null
    ];
    val
  ];
  GetParam[key_, default_] := If[GetParam[key] === Null, default, GetParam[key]];

  prec    = GetParam["prec", 32];
  iround  = 1000;
  a0      = SetPrecision[Round[GetParam["a"] * iround] / iround, prec];
  r0      = SetPrecision[Round[GetParam["r0"] * iround] / iround, prec];
  mm      = GetParam["m"];
  lmax    = GetParam["lmax"];
  lplot   = GetParam["lplot"];
  dformat = ToString @ GetParam["dformat", "Real32"];

  (* Radial grids *)
  tmp              = Import[datagpath <> "lm_rsL_" <> ToString[iConfig] <> ".dat"];
  {rstarsL, rsL}   = Transpose[tmp];
  tmp              = Import[datagpath <> "lm_rsR_" <> ToString[iConfig] <> ".dat"];
  {rstarsR, rsR}   = Transpose[tmp];

  (* Ingoing (in) solution *)
  fn    = datagpath <> "lm_in_" <> ToString[iConfig] <> ".bin";
  dataL = Import[fn, dformat];
  dimsL = {Length[rsL], 20, lmax + 1};
  GetIndexL[ri_, qi_, ll_] :=
    (ri - 1) * dimsL[[2]] * dimsL[[3]] + (qi - 1) * dimsL[[3]] + ll + 1;
  compL = Table[
    dataL[[GetIndexL[ri, qi, ll]]] + I * dataL[[GetIndexL[ri, qi + 10, ll]]],
    {ll, 0, lmax}, {qi, 1, 10}, {ri, 1, Length[rsL]}
  ];

  (* Outgoing (up) solution *)
  fn    = datagpath <> "lm_up_" <> ToString[iConfig] <> ".bin";
  dataR = Import[fn, dformat];
  dimsR = {Length[rsR], 20, lmax + 1};
  GetIndexR[ri_, qi_, ll_] :=
    (ri - 1) * dimsR[[2]] * dimsR[[3]] + (qi - 1) * dimsR[[3]] + ll + 1;
  compR = Table[
    dataR[[GetIndexR[ri, qi, ll]]] + I * dataR[[GetIndexR[ri, qi + 10, ll]]],
    {ll, 0, lmax}, {qi, 1, 10}, {ri, 1, Length[rsR]}
  ];

  (* Stitch ingoing (reversed) and outgoing across the orbit *)
  subrs = Flatten[Join[{Reverse[rsL], rsR[[2 ;;]]}]];
  psis  = Table[
    Flatten[Join[{Reverse[compL[[l + 1, q]]], compR[[l + 1, q, 2 ;;]]}]],
    {l, 0, lmax}, {q, 1, 10}
  ];

  Association[{"h1Dat" -> psis, "rvals" -> subrs, "lmax" -> lmax}]
];


(* Load the static (m=0) sector only. *)
RetStaticSamDImportRetFunction[a_, rp_, mmax_, lmax_, filepathin_] :=
  Module[{tmp, tmph1ret, tmprvals, tmplmax},

  tmp      = RetSamDImportRetFunctionInteral[a, rp, 0, filepathin];
  tmph1ret = "h1Dat" /. tmp;
  tmprvals = "rvals" /. tmp;

  Association[{
    "h1Data1D" -> tmph1ret,
    "rvals"    -> tmprvals,
    "lmax"     -> tmplmax,
    "mmax"     -> mmax,
    "a0"       -> a,
    "rOrbit"   -> rp,
    "DataType" -> "1RadialData"
  }]
];


(* Load all m-modes, apply the spheroidal-to-spherical trace transformation
   on component q=10, then trim to the working radial window. *)
RetRadiativeSamDImportRetFunction[a_, rp_, mmax_, lmax_, filepathin_] :=
  Module[{
    tmp, tmph1rettmp, tmph1ret, tmprvals, tmplmax, tmplmaxNew,
    hTraceSphericallm, h1Radiative, rvals2
  },

  tmp         = Table[RetSamDImportRetFunctionInteral[a, rp, m, filepathin], {m, 0, mmax}];
  tmph1rettmp = "h1Dat" /. tmp;

  (* Stitch m=0 full data with m>0 data (dropping the first r-point to avoid duplication) *)
  tmph1ret = Join[{tmph1rettmp[[1, ;;, ;;, ;;]]}, tmph1rettmp[[2 ;;, ;;, ;;, 2 ;;]]];
  tmprvals = "rvals" /. tmp;
  tmplmax  = ("lmax" /. tmp)[[3]];

  (* Transform the trace (q=10) from spheroidal to spherical harmonics.
     Captured variables tmph1ret and tmplmax come from the enclosing Module. *)
  hTraceSphericallm[aa_, r00_, m_] := Module[{lSumMax, IlmVec, Mixing, Jlms},
    lSumMax = tmplmax;
    IlmVec  = Table[tmph1ret[[m + 1, \[Alpha]\[Alpha] + 1, 10]], {\[Alpha]\[Alpha], Abs[m], lSumMax}];
    Mixing  = SWSHb[0, m, (aa m) / (r00^(3/2) + a), "lmax" -> lSumMax];
    Jlms    = Table[
      Sum[Mixing[[\[Beta]\[Beta], lli]] IlmVec[[\[Beta]\[Beta]]], {\[Beta]\[Beta], Length[IlmVec]}],
      {lli, Length[IlmVec]}
    ];
    Jlms
  ];

  tmplmaxNew = tmph1ret;
  Table[
    tmplmaxNew[[mi + 1, mi + 1 ;;, 10]] = hTraceSphericallm[a, rp, mi],
    {mi, 0, mmax}
  ];

  h1Radiative = Transpose[(tmplmaxNew)[[;;, ;; lmax + 1, ;;, 10 ;; 180]], 1 <-> 2];
  rvals2      = (tmprvals)[[;;, 12 ;; 182]];

  Association[{
    "h1Data1D" -> (h1Radiative)[[;;, 2 ;;]],
    "rvals"    -> (rvals2)[[;;]],
    "lmax"     -> tmplmax,
    "mmax"     -> mmax,
    "a0"       -> a,
    "rOrbit"   -> rp,
    "DataType" -> "1RadialData"
  }]
];


(* ==================================================================== *)
(* Section 2: Data import — Maarten's HDF5 format                       *)
(* ==================================================================== *)

(* Load a single m-mode and tensor component from Maarten's HDF5 file.
   index runs 1-10 corresponding to the ten independent tetrad components.
   Returns {rvals, datavalsNew} with datavalsNew indexed as [[l, r]]. *)
LoadMaartenFunction[mi_, mmax_, lmax_, index_, a_, rOrbit_, fullpath_] :=
  Module[{loadpath, comps, rIn, rUp, hIn, hUp, rvals, datavalsNew},

  loadpath = fullpath <>
    "/h1_a"  <> ToString[NumberForm[Round[a,      0.1], {Infinity, 1}]] <>
    "_rp"    <> ToString[NumberForm[Round[rOrbit, 0.1], {Infinity, 1}]] <>
    "_l"     <> ToString[lmax] <>
    "_m"     <> ToString[mi] <> ".h5";

  comps = {
    "h_l+l+", "h_l-l-", "h_m+m+", "h_m-m-",
    "rho_h_l+m+", "rhob_h_l+m-", "rhob_h_l-m+", "rho_h_l-m-",
    "sigma_delta_h_l+l-", "h"
  };

  rIn = Import[loadpath, {"HDF5", "/m_" <> ToString[mi] <> "/r_in"}];
  rUp = Import[loadpath, {"HDF5", "/m_" <> ToString[mi] <> "/r_up"}];
  hIn = Import[loadpath, {"HDF5", "/m_" <> ToString[mi] <> "/In/" <> comps[[index]]}];
  hUp = Import[loadpath, {"HDF5", "/m_" <> ToString[mi] <> "/Up/" <> comps[[index]]}];

  rvals = Join[Reverse[rIn][[;; -2]], rUp];

  datavalsNew = Table[
    Join[
      Reverse[(("Re" /. hIn) + I ("Im" /. hIn))[[li + 1]]][[;; -2]],
      (("Re" /. hUp) + I ("Im" /. hUp))[[li + 1]]
    ],
    {li, lmax}
  ];

  {rvals, datavalsNew}
];


(* Load all m-modes from Maarten's HDF5 format. *)
RetRadiativeMaartenImportRetFunction[a_, rp_, mmax_, lmax_, filepathin_] :=
  Module[{h1Radiative},

  h1Radiative = Table[
    LoadMaartenFunction[mi, mmax, lmax, index, a, rp, filepathin],
    {mi, mmax}, {index, 10}
  ];

  Association[{
    "h1Data1D" -> Transpose[Transpose[(h1Radiative)[[;;, ;;, 2]], 2 <-> 3], 1 <-> 2],
    "rvals"    -> (h1Radiative)[[1, 1, 1]],
    "lmax"     -> lmax,
    "mmax"     -> mmax,
    "a0"       -> a,
    "rOrbit"   -> rp,
    "DataType" -> "1RadialData"
  }]
];


(* ==================================================================== *)
(* Section 3: Joining static and radiative pieces                       *)
(* ==================================================================== *)

(* Merge static (m=0) and radiative 1D data structs onto the radiative
   radial grid. Nearest-neighbour matching aligns the static grid indices. *)
JoininRadiativeandStaticPieces[aVal_, rOrbit_, mmax_, lmax_, RadStruct_, StaticStruct_] :=
  Module[{
    rvals, rvalsStatic,
    staticIndicesInner, staticIndicesOuter,
    h1Radiative, h1static, Joineddata, InheritedFields
  },

  rvals       = "rvals" /. RadStruct;
  rvalsStatic = "rvals" /. StaticStruct;

  staticIndicesInner = Nearest[rvalsStatic -> "Index", rvals[[1]]][[1]];
  staticIndicesOuter = Nearest[rvalsStatic -> "Index", rvals[[-1]]][[1]];

  h1Radiative = PadLeft[
    Transpose[("h1Data1D" /. RadStruct), 1 <-> 2][[;; mmax, ;; lmax, ;;]],
    Dimensions[Transpose[("h1Data1D" /. RadStruct), 1 <-> 2][[;; mmax, ;; lmax, ;;]]] + {0, 1, 0, 0}
  ];
  h1static = ("h1Data1D" /. StaticStruct)[[;; lmax + 1, ;;, staticIndicesInner ;; staticIndicesOuter]];

  Joineddata      = Join[{h1static}, h1Radiative];
  InheritedFields = KeyDrop[RadStruct, {"h1Data1D", "rvals"}];

  Join[
    <|"h1Data1D" -> Transpose[Joineddata, 1 <-> 2], "rvals" -> rvals|>,
    InheritedFields
  ]
];


(* ==================================================================== *)
(* Section 4: Spherical harmonics and 2D grid construction              *)
(* ==================================================================== *)

(* Build the full table of spin-weighted spherical harmonics.
   Note: function name preserves original spelling for notebook compatibility. *)
SphericalHamronicsContructor[lmax_, mmax_, \[Theta]vals_] :=
  Module[{tmp, SpinZeroSpheroidalharmonic},

  tmp = Quiet[Table[
    If[ll < Max[Abs[ss], Abs[mm]],
      0 * \[Theta]vals,
      SpinWeightedSphericalHarmonicY[ss, ll, mm, \[Theta]vals, 0]
    ],
    {ss, -2, 2}, {ll, 0, lmax}, {mm, 0, mmax}
  ]];

  Association[{
    "Ys"           -> tmp,
    "Slm"          -> SpinZeroSpheroidalharmonic,
    "lmax"         -> lmax,
    "mmax"         -> mmax,
    "\[Theta]vals" -> \[Theta]vals
  }]
];


(* Project the 1D radial coefficients h_{lm}^{(comp)}(r) onto the 2D grid by
   summing over l, weighted by Y_s^{lm}(theta).
   sws maps each of the 10 tensor components to its spin weight. *)
MakeResGrid[mmax_, lmaxVal_, ResDataStruct_, YsStruc_] :=
  Module[{sws, Ysdata, \[Theta]vals, dataInput, Inheritedfields, arrFull},

  sws = {0, 0, 2, -2, 1, -1, 1, -1, 0, 0};

  Ysdata    = N["Ys"           /. YsStruc];
  \[Theta]vals   = "\[Theta]vals"      /. YsStruc;
  dataInput = N["h1Data1D"     /. ResDataStruct];

  Inheritedfields = KeyDrop[ResDataStruct, {"h1PuncDat", "lmax"}];

  arrFull = Table[
    Sum[
      Outer[Times, dataInput[[ll + 1, mm + 1, icomp]], Ysdata[[sws[[icomp]] + 3, ll + 1, mm + 1]]],
      {ll, Abs[mm], If[mm == 0, Min[30, lmaxVal], lmaxVal]}
    ],
    {icomp, 1, 10}, {mm, 0, mmax}
  ];

  Join[<|"Data2D" -> arrFull, "lmaxVal" -> lmaxVal, "\[Theta]vals" -> \[Theta]vals|>, Inheritedfields]
];


(* ==================================================================== *)
(* Section 5: Tetrad component conversion                               *)
(* ==================================================================== *)

(* Convert the 2D grid data from Sam's tetrad convention to the Kinnersley
   tetrad for a specific index pair (ai, bi) and m-mode mi.
   Kinnersley indices: 1=l, 2=n, 3=m, 4=mbar.
   Coefficient factors encode the relationship between conventions. *)
h1$tetradSamtoKinnConvert[{ai_, bi_}, mi_, Residuleh1DataStruc2D_, {a0_}] :=
  Module[{
    Residuleh1$Data$Sam$2D, rvals, \[Theta]vals,
    r, \[Theta], rgrid, \[Theta]grid,
    \[Rho]grid, \[CapitalSigma]grid, \[CapitalDelta]grid,
    coefgrid$vec, coef$Mat
  },

  Residuleh1$Data$Sam$2D = "Data2D"      /. Residuleh1DataStruc2D;
  rvals                  = N["rvals"      /. Residuleh1DataStruc2D];
  \[Theta]vals                = N["\[Theta]vals"  /. Residuleh1DataStruc2D];

  rgrid  = ToPackedArray[Table[rvals[[i]],  {i, Length[rvals]}, {j, Length[\[Theta]vals]}]];
  \[Theta]grid = ToPackedArray[Table[\[Theta]vals[[j]], {i, Length[rvals]}, {j, Length[\[Theta]vals]}]];

  \[Rho]grid      = (r + I a0 Cos[\[Theta]])        /. {r -> rgrid, \[Theta] -> \[Theta]grid};
  \[CapitalSigma]grid = (r^2 + a0^2 Cos[\[Theta]]^2) /. {r -> rgrid, \[Theta] -> \[Theta]grid};
  \[CapitalDelta]grid = (r^2 - 2 r + a0^2)           /. {r -> rgrid, \[Theta] -> \[Theta]grid};

  coefgrid$vec = {
    1,
    -\[CapitalDelta]grid / (2 \[CapitalSigma]grid),
    1 / (Sqrt[2] \[Rho]grid),
    1 / (Sqrt[2] Conjugate[\[Rho]grid])
  };
  coef$Mat = Table[coefgrid$vec[[i]] coefgrid$vec[[j]], {i, 4}, {j, 4}];

  Which[
    {ai, bi} == {1, 1},
      coef$Mat[[1, 1]] Residuleh1$Data$Sam$2D[[1, mi + 1]],
    {ai, bi} == {1, 2} || {ai, bi} == {2, 1},
      1 / (\[CapitalSigma]grid * \[CapitalDelta]grid) coef$Mat[[1, 2]] Residuleh1$Data$Sam$2D[[9, mi + 1]],
    {ai, bi} == {1, 3} || {ai, bi} == {3, 1},
      1 / \[Rho]grid coef$Mat[[1, 3]] Residuleh1$Data$Sam$2D[[5, mi + 1]],
    {ai, bi} == {1, 4} || {ai, bi} == {4, 1},
      1 / Conjugate[\[Rho]grid] coef$Mat[[1, 4]] Residuleh1$Data$Sam$2D[[6, mi + 1]],
    {ai, bi} == {2, 2},
      coef$Mat[[2, 2]] Residuleh1$Data$Sam$2D[[2, mi + 1]],
    {ai, bi} == {2, 3} || {ai, bi} == {3, 2},
      1 / Conjugate[\[Rho]grid] coef$Mat[[2, 3]] Residuleh1$Data$Sam$2D[[7, mi + 1]],
    {ai, bi} == {2, 4} || {ai, bi} == {4, 2},
      1 / \[Rho]grid coef$Mat[[2, 4]] Residuleh1$Data$Sam$2D[[8, mi + 1]],
    {ai, bi} == {3, 3},
      coef$Mat[[3, 3]] Residuleh1$Data$Sam$2D[[3, mi + 1]],
    {ai, bi} == {3, 4} || {ai, bi} == {4, 3},
      coef$Mat[[3, 4]] (
        \[CapitalSigma]grid Residuleh1$Data$Sam$2D[[10, mi + 1]] -
        1 / \[CapitalSigma]grid Residuleh1$Data$Sam$2D[[9, mi + 1]]
      ),
    {ai, bi} == {4, 4},
      coef$Mat[[4, 4]] Residuleh1$Data$Sam$2D[[4, mi + 1]]
  ]
];


(* ==================================================================== *)
(* Section 6: Single-component and full-tensor data structs             *)
(* ==================================================================== *)

(* Build a DataRegion-backed struct for one Kinnersley-tetrad component
   (ai, bi) and m-mode mi. *)
h1$Ret$dataStruct[{ai_, bi_}, mi_, Reth1DataStruc2D_, {a0_, r0_}] :=
  Module[{tmp, rvals, \[Theta]vals, rgrid, \[Theta]grid, Inheritedfields, dataoutput},

  tmp    = h1$tetradSamtoKinnConvert[{ai, bi}, mi, Reth1DataStruc2D, {a0}];
  rvals  = "rvals"      /. Reth1DataStruc2D;
  \[Theta]vals = "\[Theta]vals"  /. Reth1DataStruc2D;

  rgrid  = ToPackedArray[Table[rvals[[i]],  {i, Length[rvals]}, {j, Length[\[Theta]vals]}]];
  \[Theta]grid = ToPackedArray[Table[\[Theta]vals[[j]], {i, Length[rvals]}, {j, Length[\[Theta]vals]}]];

  Inheritedfields = KeyDrop[Reth1DataStruc2D, {"Data2D", "h1Data1D", "lmax", "mmax"}];

  dataoutput = ToDataRegion[
    tmp,
    {rvals[[1]], \[Theta]vals[[1]]},
    {Abs[rvals[[1]] - rvals[[2]]], Abs[\[Theta]vals[[2]] - \[Theta]vals[[1]]]}
  ];

  Join[<|"Data2D" -> dataoutput, "rGrid" -> rgrid, "\[Theta]Grid" -> \[Theta]grid, "a" -> a0, "rOrbit" -> r0|>, Inheritedfields]
];


(* Shared helper: build the r and theta grids and inherited fields from a 2D struct. *)
buildGridsAndInheritedFields[Residuleh1DataStruc2D_] :=
  Module[{rvals, \[Theta]vals, rgrid, \[Theta]grid, Inheritedfields},
  rvals  = "rvals"      /. Residuleh1DataStruc2D;
  \[Theta]vals = "\[Theta]vals"  /. Residuleh1DataStruc2D;
  rgrid  = ToPackedArray[Table[rvals[[i]],  {i, Length[rvals]}, {j, Length[\[Theta]vals]}]];
  \[Theta]grid = ToPackedArray[Table[\[Theta]vals[[j]], {i, Length[rvals]}, {j, Length[\[Theta]vals]}]];
  Inheritedfields = KeyDrop[Residuleh1DataStruc2D, {"Data2D", "h1ResData2D", "h1ResDat", "h1PuncDat", "lmax", "mmax"}];
  {rvals, \[Theta]vals, rgrid, \[Theta]grid, Inheritedfields}
];


(* Assemble the full 4x4 tetrad tensor on the 2D grid. *)
TensorConstruct2DData[mi_, SingleCompStructFunc_, Residuleh1DataStruc2D_, {a0_, r0_}] :=
  Module[{rvals, \[Theta]vals, rgrid, \[Theta]grid, Inheritedfields, FullTetradGrid},

  {rvals, \[Theta]vals, rgrid, \[Theta]grid, Inheritedfields} = buildGridsAndInheritedFields[Residuleh1DataStruc2D];

  FullTetradGrid = Table[
    "Data2D" /. SingleCompStructFunc[{ai, bi}, mi, Residuleh1DataStruc2D, {a0, r0}],
    {ai, 4}, {bi, 4}
  ];

  Join[<|
    "Data2D" -> FullTetradGrid,
    "rGrid"  -> rgrid, "\[Theta]Grid" -> \[Theta]grid,
    "a"      -> a0,    "rOrbit"       -> r0, "mVal" -> mi
  |>, Inheritedfields]
];


(* Same as TensorConstruct2DData but converts DataRegions to plain arrays
   for HDF5 export (uses ToListOfData from SimulationTools). *)
TensorConstruct2DDataSave[mi_, SingleCompStructFunc_, Residuleh1DataStruc2D_, {a0_, r0_}] :=
  Module[{rvals, \[Theta]vals, rgrid, \[Theta]grid, Inheritedfields, FullTetradGrid},

  {rvals, \[Theta]vals, rgrid, \[Theta]grid, Inheritedfields} = buildGridsAndInheritedFields[Residuleh1DataStruc2D];

  FullTetradGrid = Table[
    "Data2D" /. SingleCompStructFunc[{ai, bi}, mi, Residuleh1DataStruc2D, {a0, r0}],
    {ai, 4}, {bi, 4}
  ];

  Join[<|
    "Data2D" -> Table[ToListOfData[FullTetradGrid[[i, j]]], {i, 4}, {j, 4}],
    "rGrid"  -> rgrid, "\[Theta]Grid" -> \[Theta]grid,
    "a"      -> a0,    "rOrbit"       -> r0, "mVal" -> mi
  |>, Inheritedfields]
];


(* ==================================================================== *)
(* Section 7: Boyer-Lindquist coordinate conversion                     *)
(* ==================================================================== *)

(* Contract Kinnersley tetrad components with BL covectors to get h_{mu nu}.
   The tetrad covectors in BL coords are l, n, m, mbar (indices 1-4).
   The contraction uses Z = {-n, -l, mbar, m} following the NP convention.

   Note: the 'a' symbol in ndBLK (rather than a0) is intentional — it is
   resolved to a0 via the subsequent replacement rules. *)
DataStructKinnerslyTetradToBLComp[DataStruc2D_, {a0_, r0_}] :=
  Module[{
    rvals, \[Theta]vals, rgridtmp, \[Theta]gridtmp, rgrid, \[Theta]grid,
    r, \[Theta], \[CapitalDelta], \[CapitalSigma],
    ldBLK, ndBLK, mdBLK, mbdBLK,
    ldBLGridK, ndBLGridK, mdBLGridK, mbdBLGridK,
    ZTildaGrid, Retarded$tetrad$data, D2GBL, Inheritedfields
  },

  rvals  = "rvals"      /. DataStruc2D;
  \[Theta]vals = "\[Theta]vals"  /. DataStruc2D;

  Retarded$tetrad$data = "Data2D" /. DataStruc2D;

  rgridtmp  = N[ToPackedArray[Table[rvals[[i]],  {i, Length[rvals]}, {j, Length[\[Theta]vals]}]]];
  \[Theta]gridtmp = N[ToPackedArray[Table[\[Theta]vals[[j]], {i, Length[rvals]}, {j, Length[\[Theta]vals]}]]];

  rgrid  = ToDataRegion[rgridtmp,  {rvals[[1]], \[Theta]vals[[1]]}, {Abs[rvals[[1]] - rvals[[2]]], Abs[\[Theta]vals[[2]] - \[Theta]vals[[1]]]}];
  \[Theta]grid = ToDataRegion[\[Theta]gridtmp, {rvals[[1]], \[Theta]vals[[1]]}, {Abs[rvals[[1]] - rvals[[2]]], Abs[\[Theta]vals[[2]] - \[Theta]vals[[1]]]}];

  \[CapitalDelta] = r^2 - 2 r + a0^2;
  \[CapitalSigma] = r^2 + a0^2 Cos[\[Theta]]^2;

  (* Kinnersley tetrad covectors in BL coordinates *)
  ldBLK  = {-Ones, \[CapitalSigma] / \[CapitalDelta], Zero, a0 Sin[\[Theta]]^2};
  ndBLK  = 1/2 {-\[CapitalDelta] / \[CapitalSigma], -Ones, Zero, a Sin[\[Theta]]^2 \[CapitalDelta] / \[CapitalSigma]};
  mdBLK  = 1 / (Sqrt[2] (r + I a0 Cos[\[Theta]])) {-I a0 Sin[\[Theta]], Zero, \[CapitalSigma], I (a0^2 + r^2) Sin[\[Theta]]};
  mbdBLK = 1 / (Sqrt[2] (r - I a0 Cos[\[Theta]])) { I a0 Sin[\[Theta]], Zero, \[CapitalSigma], -I (a0^2 + r^2) Sin[\[Theta]]};

  (* Evaluate on the numerical grid; Zero/Ones are placeholder grid arrays *)
  ldBLGridK  = ldBLK  /. {M -> 1, a -> a0, \[Theta] -> ToPackedArray[\[Theta]grid], r -> ToPackedArray[rgrid], Zero -> 0 * ToPackedArray[rgrid], Ones -> ToPackedArray[rgrid] / ToPackedArray[rgrid]};
  ndBLGridK  = ndBLK  /. {M -> 1, a -> a0, \[Theta] -> ToPackedArray[\[Theta]grid], r -> ToPackedArray[rgrid], Zero -> 0 * ToPackedArray[rgrid], Ones -> ToPackedArray[rgrid] / ToPackedArray[rgrid]};
  mdBLGridK  = mdBLK  /. {M -> 1, a -> a0, \[Theta] -> ToPackedArray[\[Theta]grid], r -> ToPackedArray[rgrid], Zero -> 0 * ToPackedArray[rgrid], Ones -> ToPackedArray[rgrid] / ToPackedArray[rgrid]};
  mbdBLGridK = mbdBLK /. {M -> 1, a -> a0, \[Theta] -> ToPackedArray[\[Theta]grid], r -> ToPackedArray[rgrid], Zero -> 0 * ToPackedArray[rgrid], Ones -> ToPackedArray[rgrid] / ToPackedArray[rgrid]};

  ZTildaGrid = {-ndBLGridK, -ldBLGridK, mbdBLGridK, mdBLGridK};

  D2GBL = Monitor[
    Table[
      Sum[ZTildaGrid[[k, \[Mu]]] ZTildaGrid[[l, \[Nu]]] Retarded$tetrad$data[[k, l]], {k, 4}, {l, 4}],
      {\[Mu], 4}, {\[Nu], 4}
    ],
    {\[Mu], \[Nu]}
  ];

  Inheritedfields = KeyDrop[DataStruc2D, {"Data2D", "lmax", "mmax"}];

  Join[<|
    "Data2D" -> D2GBL,
    "rGrid"  -> rgrid, "\[Theta]Grid" -> \[Theta]grid,
    "a"      -> a0,    "rOrbit"       -> r0
  |>, Inheritedfields]
];


End[];

EndPackage[];
