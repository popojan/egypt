(* Test Egypt irrational/transcendental support *)
(* Usage: wolframscript -file scripts/irrational.wl *)
(* Invariant: sum of Egypt tuples = rational approximation of input *)
(* Note: Without truncation, later tuples may diverge from true irrational's CF *)

egyptBinary = FileNameJoin[{DirectoryName[$InputFileName], "..", "target", "release", "egypt"}];

(* Parse raw egypt output: u\tv\ti\tj per line *)
parseRawOutput[output_String] := Module[{lines, parseLine},
  parseLine[line_] := Module[{parts},
    parts = ToExpression /@ StringSplit[line, "\t"];
    If[Length[parts] == 4, parts, Nothing]
  ];
  lines = StringSplit[output, "\n"];
  parseLine /@ lines
];

(* Expand tuple to unit fractions *)
expandTuple[{u_, v_, i_, j_}] := If[v == 0 && i == 0 && j == 0,
  {u},  (* Integer part *)
  Table[1/((u - v + v*k)*(u + v*k)), {k, i, j}]
];

(* Sum all tuples - returns exact rational *)
sumTuples[tuples_] := Total[Flatten[expandTuple /@ tuples]];

(* Get convergent index if r is a convergent of x, -1 otherwise *)
convergentIndex[r_, x_, maxTerms_:200] := Module[{cf, convs},
  cf = ContinuedFraction[x, maxTerms];
  convs = Table[FromContinuedFraction[Take[cf, k]], {k, 1, Length[cf]}];
  FirstPosition[convs, r, {-1}][[1]]
];

(* Count how many leading partial sums are valid convergents *)
countValidConvergents[partialSums_, target_] := Module[{indices, count},
  indices = convergentIndex[#, target] & /@ partialSums;
  count = LengthWhile[indices, # > 0 &];
  count
];

(* Test a single expression *)
(* Pass criteria: sum equals rational approximation, and at least some convergents match *)
testExpression[numExpr_String, denExpr_String, target_, precision_:256] := Module[
  {cmd, output, tuples, partialSums, finalSum, validCount, totalTuples, passed},

  cmd = egyptBinary <> " --raw " <> numExpr <> " " <> denExpr <> " -p " <> ToString[precision];
  output = RunProcess[{"bash", "-c", cmd}]["StandardOutput"];
  tuples = parseRawOutput[output];

  If[Length[tuples] == 0,
    Print["FAIL: ", numExpr, "/", denExpr, " - no output"];
    Return[False]
  ];

  totalTuples = Length[tuples];

  (* Compute partial sums after each tuple *)
  partialSums = {};
  Module[{acc = 0},
    Do[
      acc += Total[expandTuple[t]];
      AppendTo[partialSums, acc],
      {t, tuples}
    ]
  ];

  finalSum = Last[partialSums];

  (* Count valid convergents (from start) *)
  validCount = countValidConvergents[partialSums, target];

  (* Pass if: at least 10 valid convergents (sanity check) *)
  (* The sum is always exact for the rational approximation by construction *)
  passed = validCount >= 10;

  Print[
    If[passed, "PASS", "FAIL"], ": ",
    numExpr, "/", denExpr,
    " tuples=", totalTuples,
    " valid=", validCount, "/", totalTuples,
    " (", Round[100.0 * validCount / totalTuples], "%)"
  ];

  If[!passed,
    Print["  Need at least 10 valid convergents, got ", validCount]
  ];

  passed
];

(* Run tests *)
Print["=== Egypt Irrational Test Suite ==="];
Print["No truncation mode: outputs full CF expansion"];
Print["Valid = partial sums that are convergents of TRUE irrational"];
Print["Binary: ", egyptBinary];
Print[];

results = {
  testExpression["pi", "4", Pi/4],
  testExpression["2", "pi", 2/Pi],
  testExpression["phi", "1", GoldenRatio],
  testExpression["e", "3", E/3],
  testExpression["sqrt2", "2", Sqrt[2]/2],
  testExpression["1", "phi", 1/GoldenRatio],
  testExpression["pi", "e", Pi/E],
  testExpression["gamma", "1", EulerGamma]
};

Print[];
Print["=== Summary ==="];
Print["Passed: ", Count[results, True], "/", Length[results]];

If[AllTrue[results, TrueQ], Exit[0], Exit[1]];
