#!/usr/bin/env wolframscript


AppendTo[$Path, "/lustre/hpc/astro/dyson/.Wolfram/Applications/KerrPunctures"];
AppendTo[$Path, "/lustre/hpc/astro/dyson/.Wolfram/Applications"];
AppendTo[$Path, DirectoryName[$InputFileName]];

Needs["KerrPunctures`"];
Needs["SpinWeightedSpheroidalHarmonics`"]

RunPuncs[aInput_, rInput_, ll_, mm_, rvalsInput_,
   Prec_ : MachinePrecision, PrecGoal_ : 10] :=
  Module[{a0, r0, M0, func$lplp, func$lmlm, func$mpmp, func$mmmm,
    func$lpmp, func$lpmm, func$lmmp, func$lmmm, func$lplm, func$trace,
     func$list, s$list, method$opt,
    punc$test$worldtube$list$comp$temp2, rvals},
   a0 = SetPrecision[aInput, Prec];
   r0 = SetPrecision[rInput, Prec];
   rvals = SetPrecision[rvalsInput, Prec];
   M0 = 1;

   punc$hs1 =
    SetPrecision[LoadKerrPuncture["hS1", M0, a0, r0], Prec];

   Clear[integrand$NB];

   window[\[Theta]_, x0_, w_, q_, s_] :=
    Piecewise[{{0, \[Theta] <= x0}, {1/2 +
        Tanh[(s*(-(q^2*Cot[(Pi*(-x0 + \[Theta]))/(2*w)]) +
               Tan[(Pi*(-x0 + \[Theta]))/(2*w)]))/Pi]/2,
       x0 < \[Theta] < w + x0}, {1,
       Pi - w - x0 >= \[Theta] >= w + x0}, {1/2 +
        Tanh[(s*(-(q^2*Cot[(Pi*(Pi - x0 - \[Theta]))/(2*w)]) +
               Tan[(Pi*(Pi - x0 - \[Theta]))/(2*w)]))/Pi]/2,
       Pi - x0 > \[Theta] > Pi - w - x0}}, 0];
   WindowFunc =
    Function[{Global`\[Theta]},
     window[Global`\[Theta], 1/100, 6 \[Pi]/16, 1, 3], Listable];


   integrand$NB[s_, l_, m_, r_, \[Theta]_, func_] :=
    Module[{punc, harm$l},
     If[l < Abs[s] || Abs[m] > l, Return[0*r] ];
     punc = WindowFunc[\[Theta]] func[r, \[Theta], m];
     harm$l =
      Sin[\[Theta]] SpinWeightedSphericalHarmonicY[s, l, m, \[Theta],
        0];
     Return[punc harm$l]
     ];

   func$lplp[r_, \[Theta]_, m_] :=
    (punc$hs1["ll", r, \[Theta], m]);
   func$lmlm[r_, \[Theta]_, m_] :=
    4 (r^2 + a0^2 Cos[\[Theta]]^2)^2/(r^2 + a0^2 - 2 r)^2 (punc$hs1[
        "nn", r, \[Theta], m]);
   func$mpmp[r_, \[Theta]_, m_] :=
    2 (r + I a0 Cos[\[Theta]])^2 (punc$hs1["mm", r, \[Theta], m]);
   func$mmmm[r_, \[Theta]_, m_] :=
    2 (r - I a0 Cos[\[Theta]])^2 (punc$hs1["mbmb", r, \[Theta],
        m]);
   func$lpmp[r_, \[Theta]_, m_] :=
    Sqrt[2] (r + I a0 Cos[\[Theta]])^2 (punc$hs1["lm", r, \[Theta],
        m]);
   func$lpmm[r_, \[Theta]_, m_] :=
    Sqrt[2] (r - I a0 Cos[\[Theta]])^2 (punc$hs1["lmb", r, \[Theta],
        m]);
   func$lmmp[r_, \[Theta]_, m_] :=
    -2 Sqrt[2] (r^2 + a0^2 Cos[\[Theta]]^2)^2/(r^2 + a0^2 - 2 r) (punc$hs1[
        "nm", r, \[Theta], m]);
   func$lmmm[r_, \[Theta]_, m_] :=
    -2 Sqrt[2] (r^2 + a0^2 Cos[\[Theta]]^2)^2/(r^2 + a0^2 - 2 r) (punc$hs1[
        "nmb", r, \[Theta], m]);
   func$lplm[r_, \[Theta]_, m_] :=
    -2 (r^2 + a0^2 Cos[\[Theta]]^2)^2 (punc$hs1["ln",
        r, \[Theta], m]);
   func$trace[r_, \[Theta]_, m_] :=
    -2 (punc$hs1["ln", r, \[Theta], m] -
        punc$hs1["mmb", r, \[Theta], m]);

   func$list = {func$lplp, func$lmlm, func$mpmp, func$mmmm, func$lpmp,
      func$lpmm, func$lmmp, func$lmmm, func$lplm, func$trace};
   s$list = {0, 0, 2, -2, 1, -1, 1, -1, 0, 0};

   method$opt = {Automatic, "SymbolicProcessing" -> 0};
   punc$test$worldtube$list$comp$temp2 =
    Table[If[(s$list[[i]] == 0 && OddQ[ll + mm]) ||
       ll < Abs[s$list[[i]]], 0*rvals,

      NIntegrate[
       SetPrecision[
        integrand$NB[s$list[[i]], ll, mm, rvals, \[Theta],
         func$list[[i]]], Prec], {\[Theta], 0, Pi/2, Pi},
       Exclusions -> "Singularities", "Method" -> method$opt,
       PrecisionGoal -> PrecGoal, WorkingPrecision -> Prec - 4]]

     , {i, 1, 10}];
   Return[punc$test$worldtube$list$comp$temp2];
   ];


(* ── Parameters — edit these before submitting ── *)
a0       = ToExpression[$ScriptCommandLine[[2]]];
r0       = ToExpression[$ScriptCommandLine[[3]]];
rmin     = ToExpression[$ScriptCommandLine[[4]]];
rmax     = ToExpression[$ScriptCommandLine[[5]]];
rstep    = ToExpression[$ScriptCommandLine[[6]]];
lmax     = ToExpression[$ScriptCommandLine[[7]]];
mmax     = ToExpression[$ScriptCommandLine[[8]]];
savepath = $ScriptCommandLine[[9]];

workingprec = 100;
precgoal    = 16;

rvals = Table[SetPrecision[i, workingprec], {i, rmin, rmax, rstep}];

Print["Job started: a=", a0, " r0=", r0, " lmax=", lmax, " mmax=", mmax];
Print["rvals = ", rvals];

(* ── Directory setup ── *)
atag    = IntegerString[Round[1000*a0], 10, 3];
datatag = atag <> "-" <> ToString[Round[r0*10, 1]];
datadir = savepath <> "/data" <> atag <> "/data" <> datatag <> "/data";
If[!DirectoryQ[datadir], CreateDirectory[datadir, CreateIntermediateDirectories -> True]];

(* ── Parallel setup ── *)
jobStart = AbsoluteTime[];
Print["License server port open: ", Length[Select[ReadList["!ss -tlnp 2>/dev/null", String], StringContainsQ[#, "16286"]&]] > 0];
LaunchKernels[];
Print["Kernels launched: ", Length[Kernels[]]];
ParallelEvaluate[AppendTo[$Path, "/lustre/hpc/astro/dyson/.Wolfram/Applications/KerrPunctures"]];
ParallelEvaluate[AppendTo[$Path, "/lustre/hpc/astro/dyson/.Wolfram/Applications"]];
ParallelEvaluate[AppendTo[$Path, DirectoryName[$InputFileName]]];
ParallelNeeds["KerrPunctures`"];
ParallelNeeds["SpinWeightedSpheroidalHarmonics`"];
DistributeDefinitions[a0, r0, rvals, workingprec, precgoal, window, WindowFunc, RunPuncs, integrand$NB];

(* ── Compute (parallel) ── *)
totalPairs = Sum[lmax - mm + 1, {mm, 1, mmax}];
Print["Total (mm,ll) pairs to compute: ", totalPairs];
allResults = ParallelTable[
  Module[{res, t},
    {t, res} = AbsoluteTiming[RunPuncs[a0, r0, ll, mm, rvals, workingprec, precgoal]];
    Print["Done: mm=", mm, " ll=", ll, " in ", Round[t, 0.1], "s"];
    res
  ],
  {mm, 1, mmax}, {ll, mm, lmax}
];

(* ── Save HDF5 ── *)
numStr[x_] := Module[{s = ToString[x]},
  Which[StringMatchQ[s, ___ ~~ "."], s <> "0",
        !StringContainsQ[s, "."],    s <> ".0",
        True,                        s]];

comps = {"h_l+l+", "h_l-l-", "h_m+m+", "h_m-m-", "rho_h_l+m+",
         "rhob_h_l+m-", "rhob_h_l-m+", "rho_h_l-m-", "sigma_delta_h_l+l-", "h"};

tolist[x_] := Developer`ToPackedArray[N[x]];

nr    = Length[rvals];
zeros = Table[0. + 0. I, {nr}];

normaliseComp[x_] := Which[
  ListQ[x] && Length[x] == nr, N[x],
  ListQ[x],                    N[Flatten[x, 1]],
  True,                        Table[N[x], {nr}]
];

Table[
  Module[{fn, datasets},
    fn = datadir <> "/h1_a" <> numStr[Round[a0, 0.01]] <>
         "_rp" <> numStr[Round[r0, 0.01]] <>
         "_l" <> ToString[lmax] <> "_m" <> ToString[mm] <> ".h5";
    datasets = Flatten[Table[
      With[{result = If[ll < mm,
              Table[zeros, {10}],
              Table[normaliseComp[allResults[[mm, ll - mm + 1, i]]], {i, 1, 10}]]},
        Table["/m_" <> ToString[mm] <> "/l_" <> ToString[ll] <> "/" <> comps[[i]] -> tolist[result[[i]]], {i, 1, 10}]
      ],
      {ll, 0, lmax}], 1];
    datasets = Append[datasets, "/m_" <> ToString[mm] <> "/rvals" -> tolist[rvals]];
    Export[fn, datasets, "HDF5"];
    Print["Saved m=", mm, " lmax=", lmax];
  ],
  {mm, 1, mmax}
];

Print["Total job time: ", Round[(AbsoluteTime[] - jobStart)/60, 0.1], " minutes"];
