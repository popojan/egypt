#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*lo1[Rational[x_,y_]]:=1/# MixedFractionParts[x/y #]&@If[y-x==1,x,Mod[First@Last@ExtendedGCD[x,y]+ y,y]]*)
(*ToPow[x0_]:=Module[{x=x0,s=0},While[Mod[x,2]==0,s+=1;x/=2];If[s>1,2^HoldForm@#&@s x,x 2^s]]*)
(*sum[b_,v_,to_]:=If[to==1,1/((b-v)b ),HoldForm[Sum[#1,{k,0,#2}]]&@@{1/((b-(k+1) v)(b-k v) ),to-1}]*)
(*EgyptExpandAll[eg_]:=ReleaseHold[#/.Sum->Table]&/@eg//Flatten//Sort//Reverse*)
(*EgyptTermCount[eg_]:=If[Depth@#>1,1+Part[#/.Sum->Table,1,2,3],1]&/@eg//Total*)
(*EgyptExpandUpTo[eg_,lim_]:=If[1+Part[#/.Sum->Table,1,2,3]<=lim,ReleaseHold[#/.Sum->Table],#]&/@eg//Flatten//SortBy[#,ReleaseHold]&//Reverse*)
(*Etc[a_,0]:=a*)
(*Etc[Rational[a_,b_],v_]:=sum[b,v,Floor@((a b)/(1+a v))]*)
(*EgyptNew[Rational[a_,b_]]:=Module[{ret={},x,rest=a/b},*)
(*While[rest>0,*)
(*x=Etc[rest,Denominator@rest-Denominator@First@lo1[rest]];PrependTo[ret,x];rest=rest-ReleaseHold@ReleaseHold@First@ret];ret]*)


(* ::Input:: *)
(*mr[exp_]:=(Divide@@(Power[2,MersennePrimeExponent@{exp,exp+1}]-1))*)


(* ::Input:: *)
(*EgyptNew@(2023/2024)*)


(* ::Input:: *)
(*{Total@ReleaseHold@#&,Identity,EgyptExpandAll,EgyptExpandUpTo[#,2]&}[EgyptNew@(29/67)]//Through*)


(* ::Input:: *)
(*EgyptNew@mr@5*)


(* ::Input:: *)
(*EgyptNew[Prime@1000/Prime@2000]*)


(* ::Input:: *)
(*allMethods={"New","Greedy","Harmonic","OddGreedy","BinaryRemainder","ContinuedFraction",*)
(*"Splitting","GroupedContinuedFraction" ,"Pairing", "Hybrid","Binary"};*)
(*fasterMethods = {"New","Greedy","Harmonic","OddGreedy","BinaryRemainder","ContinuedFraction"};*)
(*fastestMethods={"New","Greedy","Harmonic","OddGreedy","BinaryRemainder"};*)


(* ::Input:: *)
(*compare[Rational[x_,y_],methods_]:=Table[Flatten[{m,{AccountingForm[First@#,6],If[m=="New",Row[{Length@Last@#,"/",EgyptTermCount@Last@#}],If[m=="GroupedContinuedFraction"||m=="Hybrid",Min[Length/@Last@#],Length@Last@#]],StringLength@Compress@If[m=="BinaryRemainder",1/(ToPow/@Denominator@Last@#),Last@#]}&@If[m=="New",Timing@EgyptNew@(x/y),Timing@ResourceFunction["EgyptianFraction"][x/y,Method->m]]},1],{m,methods}]*)
(**)


(* ::Input:: *)
(*timeit[Rational[x_,y_],methods_]:=compare[x/y,methods]//SortBy[#,#[[2]]&]&//TableForm[#[[;;,2;;]],TableHeadings->{First/@#, {"time", "terms","compressed size"}}]&*)


(* ::Input:: *)
(*timeit[16/77,DeleteCases[allMethods,"Splitting"]]*)


(* ::Input:: *)
(*timeit[mr@8,fasterMethods]*)


(* ::Input:: *)
(*timeit[999999/1000000,fasterMethods]*)


(* ::Input:: *)
(*timeit[Prime@1000/Prime@2000,fasterMethods]*)


(* ::Input:: *)
(*timeit[Prime@123456789/Prime@987654321,fasterMethods]*)


(* ::Input:: *)
(*timeit[mr@15,{"New","BinaryRemainder"}]*)


(* ::Input:: *)
(*timeit[mr@19,{"New","BinaryRemainder"}]*)


(* ::Input:: *)
(*timeit[mr@33,{"New"}]*)


(* ::Input:: *)
(*EgyptNew@mr@9*)


(* ::Input:: *)
(*ContinuedFraction@mr@9*)


(* ::Input:: *)
(*Equal@Through[{*)
(*ResourceFunction["EgyptianFraction"][#,Method->"ContinuedFraction"]&,*)
(*EgyptExpandAll@EgyptNew@#&}[Divide@@#]*)
(*]&/@Partition[RandomPrime[{2^10,2^30},2000],2]//Counts*)
