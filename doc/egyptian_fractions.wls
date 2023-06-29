#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*idx[Rational[x_,y_]]:=SortBy[1/x {ChineseRemainder[{0,-1},{x,y}],ChineseRemainder[{0,1},{x,y}]},N[Round[x/y #]/#-x/y]&]*)
(*best[Rational[x_,y_],i_,d_]:=Module[{b=(i-1) (2Floor[y/2]+1)+idx[x/y][[If[d<0,1,2]]]},Round[x/y b]/b]*)
(*lo[Rational[x_,y_],n_]:=Table[best[x/y,i,-1],{i,1,n}]*)
(*hi[Rational[x_,y_],n_]:=Table[best[x/y,i,1],{i,1,n}]*)


(* ::Input:: *)
(*lo1[Rational[x_,y_]]:=1/# Round[x/y #]&@If[y-x==1,x,(1/x ChineseRemainder[{0,1},{x,y}])]*)
(*lo1[Rational[x_,y_]]:=1/# Round[x/y #]&@If[y-x==1,x,(1/x Mod[x First@Last@ExtendedGCD[x,y]+x y,x y])]*)
(*lo1[Rational[x_,y_]]:=1/# Round[x/y #]&@If[y-x==1,x,Mod[First@Last@ExtendedGCD[x,y]+ y,y]]*)
(*EgyptPrime[Rational[x_,y_],ret_]:=Module[{z0,z,retz=ret,whole={},steps=0},*)
(*z0=If[x/y>=1,AppendTo[whole,Floor[x/y]];x/y-Floor[x/y],x/y];*)
(*While[Numerator@z0>1,z=lo1[z0];PrependTo[retz,z0-z];z0=z];*)
(*Join[whole,{z0},retz]]*)
(*EgyptPrime[x_]:={x}*)
(*EgyptMerge[eg_]:=Module[{i,j,ret={}},*)
(*For[i=1,i<=Length@eg,If[Denominator@eg[[i]]>1,*)
(*j=Flatten@Position[Numerator[Accumulate@eg[[i;;]]],1];*)
(*j=If[j!={},i+Last@j-1,i];*)
(*AppendTo[ret,Total@eg[[i;;j]] ];i=j+1,*)
(*AppendTo[ret,eg[[i]]];i+=1]];*)
(*ret*)
(*]*)
(*EgyptPrime[Rational[x_,y_]]:=If[x==1,{x/y},*)
(*Module[{ret=EgyptPrime[x/y,{}],mret},*)
(*While[mret=EgyptMerge@ret;ret!=mret,ret=mret];*)
(*mret]*)
(*]*)


(* ::Input:: *)
(*EgyptPrime[lst_List]:=EgyptPrime/@lst*)


(* ::Input:: *)
(*Timing@EgyptPrime[2023/2024]*)


(* ::Input:: *)
(*coprimes[m_,n_]:={*)
(*{2m-n,m},*)
(*{2m+n,m},*)
(*{m+2n,n}*)
(*}*)


(* ::Input:: *)
(*Coprimes[n_]:=Divide@@#->EgyptPrime[Divide@@#]&/@Sort[Sort/@Select[Join@@NestList[Join@@((coprimes@@#)&/@#)&,{{2,1},{3,1}},n],Last@#>0&&First@#>0&]]*)


(* ::Input:: *)
(*Coprimes[5]*)


(* ::Input:: *)
(*Timing@EgyptPrime[Prime@3999/Prime@4000]*)
