#!/usr/bin/env wolframscript
If[Length[$ScriptCommandLine] < 3, Print["Usage: wl/test.wls <num_tests> <max_bits>"]; Exit[0]];
testnum=ToExpression[$ScriptCommandLine[[2]]];
testmag=ToExpression[$ScriptCommandLine[[3]]];
file=CreateFile[];
For[i=1,i<=testnum,
  original=Divide@@RandomInteger[{1,2^testmag},2];
  WriteString["stdout", StringJoin[{"\rTesting ", ToString@i}]];(*, ": ", Sequence@@ToString@NumeratorDenominator@original}]]*);
  output = RunProcess[Join[{"target/release/egypt"}, NumeratorDenominator@original], "StandardOutput"];
  fractions = ReadList[StringToStream@output, {Number, Number}];
  error = If[Not@AllTrue[Rest[First/@ fractions], # == 1 &], Print[{original, "bad num\n"}]; True, False];
  error = error || If[Not@DuplicateFreeQ[Last /@ fractions], Print[{original, Commonest[Last/@fractions], "dupes\n"}]; True, False];
  recovered = Divide @@@ fractions // Total;
  error = error || If[recovered != original, Print[{original, "wrong\n"}]; True, False];
  If[error, Exit[1]];
  i+=1;
]
Print["\nDone."]
DeleteFile[file]

