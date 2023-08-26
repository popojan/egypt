#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*sqrtt[x_,n_]:=1+\!\( *)
(*\*UnderoverscriptBox[\(\[Sum]\), \(j = 1\), \(n\)]*)
(*\*FractionBox[\(1\), \(1 + *)
(*\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(j\)]*)
(*\*SuperscriptBox[\(2\), \(i - 1\)] *)
(*\*FractionBox[\(Factorial[j + i]\), \(Factorial[j - i] Factorial[2  i]\)] *)
(*\*SuperscriptBox[\(x\), \(i\)]\)]\)*)
(*sqrtl[x_,n_]:=Join[{1},Table[1/(1+\!\( *)
(*\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(j\)]\( *)
(*\*SuperscriptBox[\(2\), \(i - 1\)] *)
(*\*FractionBox[\(Factorial[j + i]\), \(Factorial[j - i] Factorial[2  i]\)] *)
(*\*SuperscriptBox[\(x\), \(i\)]\)\)),{j,1,n}]]*)
(*sqrtpell[n_,acc_]:=(x-1)/y sqrtt[x-1,acc]/.Normal@First@Solve[x^2-n y^2==1,{x,y},PositiveIntegers]/.C[1]->1 (*Pell Equation*)*)
(*sqrtpelll[n_,acc_]:={(x-1)/y ,sqrtl[x-1,acc]}/.Normal@First@Solve[x^2-n y^2==1,{x,y},PositiveIntegers]/.C[1]->1 (*Pell Equation*)*)


(* ::Input:: *)
(*fmt[x_]:={Round[First@#,0.001],Last@#}&@MantissaExponent[x]*)
(*cmp[n_,k_]:=Module[{a,b},Block[{$MaxExtraPrecision=10000},*)
(*a=Timing[{HoldForm@n HoldForm@#&@(1/n #),fmt@N[#-Sqrt[n],10000]}&@sqrtpell[n,k]];*)
(*b=Timing[{#,fmt@N[#-Sqrt[n],10000]}&@Rationalize[Sqrt[n],10^(Last@Last@Last@a)]]];*)
(*{a,b}]*)


(* ::Input:: *)
(*cmp[61,7]*)


(* ::Input:: *)
(*First/@cmp[2^MersennePrimeExponent@7-1,10]*)
