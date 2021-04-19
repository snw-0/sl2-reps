#-------------------------------------------------------------
# The unary type R, i.e. R_{lambda}^{sigma} in the case sigma = lambda.
# This is defined only for odd p (because with p = 2 we have ld-si-1 = -1).
#
# Unlike other cases, here we do not decompose by character, instead considering
# the whole rep. R_ld(r); this contains two irreducible subreps. of the same level.
#
# Note that R_{lambda}^{sigma} has an additional parameter, t, but this is irrelevant here:
# t affects only the second component of M, which is trivial in this case.
#-------------------------------------------------------------

RepRUnary := function(p, ld, r, silent)
    local l, M, B_f, B_g, x, y, k, c, s, S_p, T_p, deg_p, S_m, T_m, deg_m, xi, U, V, i;

    l := p^ld;

    M := [0 .. l-1];

    # There are two types of basis elements; see notes S. 1.6.1.
    if ld = 1 then
        # f_{x, +-}
        B_f := [];
        for x in [1 .. (l-1) / 2] do
            Add(B_f, x);
        od;

        # For +, there is one additional basis element: delta_0.
        deg_p := Length(B_f) + 1;
        # For -, there are no additional elements.
        deg_m := Length(B_f);
    else
        # f_{x, +-}
        B_f := [];
        for x in PrimeResidues(l) do
            if x > (l-1) / 2 then
                break;
            else
                Add(B_f, x);
            fi;
        od;

        # g_{y, k, +-}
        B_g := [];
        for k in [1 .. (p-1) / 2] do
            Add(B_g, [0, k]);
            for y in [1 .. ((p^(ld-2)) - 1)/2] do
                Append(B_g, [[y, k], [p^(ld-2)-y, k]]);
            od;
        od;

        deg_p := Length(B_f) + Length(B_g);
        deg_m := Length(B_f) + Length(B_g);
    fi;

    # Construct T matrix.
    s := List(B_f, x -> E(l)^(r * x^2));

    if ld = 1 then
        T_m := DiagonalMat(s);
        Add(s, 1);
        T_p := DiagonalMat(s);
    else
        for x in B_g do
            Add(s, E(l)^(r * p^2 * x[1]^2));
        od;
        T_m := DiagonalMat(s);
        T_p := DiagonalMat(s);
    fi;

    # Scale factor, used for the S matrix:
    # S_Q(-1) / Sqrt(|M|) = (1/|M|) * Sum(e(Q(x))).
    c := Sum(M, x -> E(l)^(r * x^2)) / l;

    # Construct S matrix.
    S_p := NullMat(deg_p, deg_p);
    S_m := NullMat(deg_m, deg_m);

    for x in [1 .. Length(B_f)] do
        for y in [1 .. Length(B_f)] do
            # s_{f,f}
            s := 2 * r * B_f[x] * B_f[y];

            S_p[x][y] := c * (E(l)^(s) + E(l)^(-s));
            S_m[x][y] := c * (E(l)^(s) - E(l)^(-s));
        od;
    od;

    if ld = 1 then
        # Row and column corr. to delta_0 (in R_+).
        for x in [1 .. Length(B_f)] do
            S_p[x][Length(B_f) + 1] := c * Sqrt(2);
            S_p[Length(B_f) + 1][x] := c * Sqrt(2);
        od;

        S_p[Length(B_f) + 1][Length(B_f) + 1] := c;
    else
        for x in [1 .. Length(B_f)] do
            for y in [1 .. Length(B_g)] do
                # s_{f,g} and s_{g,f}
                # Note that, in general, s_{f,g} != s_{g,f}.
                S_p[x][y + Length(B_f)] := (1/Sqrt(p)) * c * Sum([0 .. (p-1)], a ->
                        E(p)^(a * B_g[y][2])
                        * (E(l)^(2*r*B_f[x]*(p*B_g[y][1] + a*p^(ld-1))) + E(l)^(-2*r*B_f[x]*(p*B_g[y][1] + a*p^(ld-1))))
                        );
                S_p[y + Length(B_f)][x] := (1/Sqrt(p)) * c * Sum([0 .. (p-1)], a ->
                        E(p)^(-a * B_g[y][2])
                        * (E(l)^(2*r*B_f[x]*(p*B_g[y][1] + a*p^(ld-1))) + E(l)^(-2*r*B_f[x]*(p*B_g[y][1] + a*p^(ld-1))))
                        );;

                S_m[x][y + Length(B_f)] := (1/Sqrt(p)) * c * Sum([0 .. (p-1)], a ->
                        E(p)^(a * B_g[y][2])
                        * (E(l)^(2*r*B_f[x]*(p*B_g[y][1] + a*p^(ld-1))) - E(l)^(-2*r*B_f[x]*(p*B_g[y][1] + a*p^(ld-1))))
                        );
                S_m[y + Length(B_f)][x] := (1/Sqrt(p)) * c * Sum([0 .. (p-1)], a ->
                        E(p)^(-a * B_g[y][2])
                        * (E(l)^(2*r*B_f[x]*(p*B_g[y][1] + a*p^(ld-1))) - E(l)^(-2*r*B_f[x]*(p*B_g[y][1] + a*p^(ld-1))))
                        );
            od;
        od;

        for x in [1 .. Length(B_g)] do
            for y in [1 .. Length(B_g)] do
                # s_{g,g}
                S_p[x + Length(B_f)][y + Length(B_f)] := (1/p) * c
                        * Sum([0 .. (p-1)], a ->
                        Sum([0 .. (p-1)], b ->
                        E(p)^(a * B_g[x][2] - b * B_g[y][2])
                        * (E(l)^(2*r*(p^2)*B_g[x][1]*B_g[y][1]) + E(l)^(-2*r*(p^2)*B_g[x][1]*B_g[y][1]))
                        ));

                S_m[x + Length(B_f)][y + Length(B_f)] := (1/p) * c
                        * Sum([0 .. (p-1)], a ->
                        Sum([0 .. (p-1)], b ->
                        E(p)^(a * B_g[x][2] - b * B_g[y][2])
                        * (E(l)^(2*r*(p^2)*B_g[x][1]*B_g[y][1]) - E(l)^(-2*r*(p^2)*B_g[x][1]*B_g[y][1]))
                        ));
            od;
        od;
    fi;

    if ld = 1 then
        return [[S_p, T_p, deg_p], [S_m, T_m, deg_m]];
    else
        U:=IdentityMat(Length(B_f)+Length(B_g));
        V:=IdentityMat(Length(B_f)+Length(B_g));
        i:=1;
        while i <= Length(B_g) do
            if B_g[i][1]>0 then
                xi := E(2*p)^(B_g[i][2])/Sqrt(2);
                U[Length(B_f)+i][Length(B_f)+i] := xi;
                U[Length(B_f)+i][Length(B_f)+i+1] := E(4)*xi;
                U[Length(B_f)+i+1][Length(B_f)+i] := xi;
                U[Length(B_f)+i+1][Length(B_f)+i+1] := -E(4)*xi;

                xi := E(2*p)^(B_g[i][2])/Sqrt(2);
                V[Length(B_f)+i][Length(B_f)+i] := E(4)*xi;
                V[Length(B_f)+i][Length(B_f)+i+1] := xi;
                V[Length(B_f)+i+1][Length(B_f)+i] := E(4)*xi;
                V[Length(B_f)+i+1][Length(B_f)+i+1] := -xi;

                i := i+2;
            else
                V[Length(B_f)+i][Length(B_f)+i]:=E(4);

                i := i+1;
            fi;
        od;
        return [[S_p^U, T_p^U, deg_p], [S_m^V, T_m^V, deg_m]];
    fi;
end;

