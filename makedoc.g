#
# SL2Reps: Constructing symmetric representations of SL(2,Z).
#
# This file is a script which compiles the package manual.
#


if fail = LoadPackage("AutoDoc", "2018.02.14") then
    Error("AutoDoc version 2018.02.14 or newer is required.");
fi;

AutoDoc( rec(
    scaffold := true,
    autodoc := true,
    gapdoc := rec(
        LaTeXOptions := rec(
            LateExtraPreamble := "\\usepackage{amsmath}"
        )
    )
) );
