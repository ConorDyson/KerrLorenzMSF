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
     punc = WindowFunc[\[Theta]] func[r, \[Theta], m];
     harm$l = 
      Sin[\[Theta]] SpinWeightedSphericalHarmonicY[s, l, m, \[Theta], 
        0];
     Return[punc harm$l]
     ];
   
   func$lplp[r_, \[Theta]_, m_] := 
    Evaluate[(punc$hs1["ll", r, \[Theta], m])];
   func$lmlm[r_, \[Theta]_, m_] := 
    Evaluate[
     4 (r^2 + a0^2 Cos[\[Theta]]^2)^2/(r^2 + a0^2 - 2 r)^2 (punc$hs1[
        "nn", r, \[Theta], m])];
   func$mpmp[r_, \[Theta]_, m_] := 
    Evaluate[
     2 (r + I a0 Cos[\[Theta]])^2 (punc$hs1["mm", r, \[Theta], m])];
   func$mmmm[r_, \[Theta]_, m_] := 
    Evaluate[
     2 (r - I a0 Cos[\[Theta]])^2 (punc$hs1["mbmb", r, \[Theta], 
        m])];
   func$lpmp[r_, \[Theta]_, m_] := 
    Evaluate[
     Sqrt[2] (r + I a0 Cos[\[Theta]])^2 (punc$hs1["lm", r, \[Theta], 
        m])];
   func$lpmm[r_, \[Theta]_, m_] := 
    Evaluate[
     Sqrt[2] (r - I a0 Cos[\[Theta]])^2 (punc$hs1["lmb", r, \[Theta], 
        m])];
   func$lmmp[r_, \[Theta]_, m_] := 
    Evaluate[-2 Sqrt[
       2] (r^2 + a0^2 Cos[\[Theta]]^2)^2/(r^2 + a0^2 - 2 r) (punc$hs1[
        "nm", r, \[Theta], m])];
   func$lmmm[r_, \[Theta]_, m_] := 
    Evaluate[-2 Sqrt[
       2] (r^2 + a0^2 Cos[\[Theta]]^2)^2/(r^2 + a0^2 - 2 r) (punc$hs1[
        "nmb", r, \[Theta], m])];
   func$lplm[r_, \[Theta]_, m_] := 
    Evaluate[-2 (r^2 + a0^2 Cos[\[Theta]]^2)^2 (punc$hs1["ln", 
        r, \[Theta], m])];
   func$trace[r_, \[Theta]_, m_] := 
    Evaluate[-2 (punc$hs1["ln", r, \[Theta], m] - 
        punc$hs1["mmb", r, \[Theta], m])];
   
   Print[func$lplp[10,0.3,2] ];

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


workingprec = 100;
rvals = Table[SetPrecision[i, workingprec], {i, 79/10, 81/10, 1/10}];
a0 = 6/10;
r0 = 8;
test1 = RunPuncs[a0, r0, 2, 2, rvals, workingprec, 10]//AbsoluteTiming;
Print[test1]   