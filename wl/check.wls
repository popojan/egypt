#!/usr/bin/env wolframscript

fractions = ReadList[$ScriptCommandLine[[4]], {Number,Number}];
If[Not@AllTrue[Rest[First/@ fractions], # == 1 &], Print["bad num"]; Exit[1]];
If[Not@DuplicateFreeQ[Last /@ fractions], Print["dupes"]; Exit[2]];
original = Divide @@ (ToExpression /@ $ScriptCommandLine[[{2,3}]]);
recovered = Divide @@@ fractions // Total;
If[recovered != original, Print[{recovered - original,"wrong"}]; Exit[3] ];
Exit[0];
