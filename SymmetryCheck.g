SymmetryCheck := function(m)
    local i,j,sym;

    sym := true;

    for i in [1 .. Length(m[1])] do
        for j in [i+1 .. Length(m[1])] do
            if m[1][i][j] <> m[1][j][i] then
                Print(i, ", ", j, ", ", m[1][i][j] / m[1][j][i], "\n");
                sym := false;
            fi;
        od;
    od;

    if sym then
        Print("Symmetric!\n");
    fi;

    if not IsDiagonalMat(m[2]) then
        Print("t is not diagonal!\n");
    fi;
end;

WhichNonSym := function(list)
    local rep;

    for rep in list do
        if rep[1] <> TransposedMat(rep[1]) or not IsDiagonalMat(r[2]) then
            Print(rep[5], "\n");
        fi;
    od;
end;
