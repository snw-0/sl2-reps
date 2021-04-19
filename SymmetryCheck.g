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
    local rho, ok;

    ok := true;

    for rho in list do
        if rho[1] <> TransposedMat(rho[1]) or not IsDiagonalMat(rho[2]) then
            Print(rho[5], "\n");
            ok := false;
        fi;
    od;
    if ok then
        Print("All symmetric.\n");
    fi;
end;
