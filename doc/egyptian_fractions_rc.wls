#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*p[{u_,v_,i_,j_},lim_]:=If[j-i+1>lim,Module[{ret={},a,b},*)
(*{a,b}=NumeratorDenominator@First@t@{{u,v,i,j}};*)
(*If[OddQ@a,AppendTo[ret,Floor[a/2]/b];AppendTo[ret,Ceiling[a/2]/b],*)
(*AppendTo[ret,(a/2-1)/b];AppendTo[ret,(a/2+1)/b]];Flatten[f/@DeleteCases[ret,0],1]*)
(*],{{u,v,i,j}}]*)
(*f[q_Rational]:=Module[{e={},x,v,a,b,t,r=q},*)
(*While[{a,b}=NumeratorDenominator@r;r>0&&b>1,*)
(*v=b-Denominator[1/# Floor[a/b #]&@If[b-a==1,a,Mod[ExtendedGCD[a,b][[2,1]]+b,b]]];*)
(*t=Floor[a b / (1+a v)];*)
(*PrependTo[e,{b-t v,v,1,t}];r-=t/(b(b-t v))];*)
(*If[r>0,Prepend[e,{r,0,1,1}],e]](*raw quadruplets*)*)
(*f0[q_List]:=If[#3==#4,1/#1/(#1-#2),HoldForm[Sum[#1,{k,#2,#3}]]&@@{1/(#1-#2 k)/(#1-#2(k-1)),#3,#4}]&@@@q(*symbolic sums of egyptian terms*)*)
(*v[q_List]:=If[#3==#4,1/(#1+#2 #3)/(#1+#2(#3-1)),HoldForm@Sum[#1,{k,#2,#3}]&@@{1/((#1+k #2)(#1 +(k-1)#2)),#3,#4}]&@@@q(*symbolic sums of egyptian terms*)*)
(*e[q_List]:=ReverseSort@Flatten[Table[1/(#1+#2 k)/(#1+#2(k-1)),{k,#3,#4}]&@@@q] (*evaluate all egyptian terms - worst case slow*)*)
(*h[q_List]:=HoldForm/@Table[1/(#1+#2(HoldForm[#]&@k))/(#1+#2(HoldForm@#&@(k-1))),{k,#3,#4}]&@@@q //Total/@#&(*evaluate all egyptian terms - worst case slow*)*)
(*t[x_List]:=(1-#3+#4)/((#1-#2+#2 #3) (#1+#2 #4))&@@@x (*evaluate symbolic sums totals*)*)
(*m[eg_]:=Module[{i,j,ret={}},*)
(*	For[i=1,i<=Length@eg,If[Denominator@eg[[i]]>1,*)
(*		j=Select[Flatten@Position[Numerator[Accumulate@eg[[i;;]]],1],Position[eg,Total@eg[[i;;i+#-1]]]=={}&&Position[ret,Total@eg[[i;;i+#-1]]]=={}&];*)
(*		j=If[j!={},i+Last@j-1,i];*)
(*		AppendTo[ret,Total@eg[[i;;j]] ];i=j+1,*)
(*		AppendTo[ret,eg[[i]]];i+=1]*)
(*];*)
(*ReverseSort@ret*)
(*]*)
(*split[eg_List,lim_]:=FixedPoint[Join@@(p[#,lim]&/@#)&,eg]*)
(*correct[eg_List]:=Module[{x=Split@eg,i},*)
(*i=FirstPosition[x,_List?(Length@#>1&)];*)
(*ReverseSort@Flatten@If[i=={},x,i=First@i;*)
(*Join[x[[;;i-1]],{e@f@Total@x[[i]]},x[[i+1;;]]]*)
(*]]*)
(*F[q_Rational]:=correct@e@split[f@q,8]*)


(* ::Input:: *)
(*mr[exp_]:=Divide@@(Power[2,MersennePrimeExponent@{exp,exp+1}]-1)*)


(* ::Input:: *)
(*Timing[Total@F@#-#]&@mr@9*)


(* ::Input:: *)
(*Timing@{F@(7/11),m@F@(7/11)}*)
