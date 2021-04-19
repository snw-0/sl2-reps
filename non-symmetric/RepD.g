#-------------------------------------------------------------
# MD gives the module and character information of type D.
# The basis for a type D representation for a primitive chi with chi^2 != 1 depends only
# on p and ld, so it can also be produced here.
#
# For p = 2 and ld = 2 or ld = 3, the rep. is reducible (for all primitive chi);
# the basis will be converted in RepD.
#-------------------------------------------------------------

MD := function(p, ld)
    local l, M, alpha, ord, omicron, Aind, Agrp, Bp, B1, Char, a;

    l := p^ld;

    M := Tuples([0 .. (l - 1)], 2);

    # For type D, A = A_ld^times.
    if p > 2 or ld < 3 then
        # A is cyclic.
        alpha := GeneratorsPrimeResidues(l).generators[1];
        ord := OrderMod(alpha, l);
        Agrp := List([0 .. ord - 1], x ->
                [x, , (alpha^x) mod l, (alpha^(ord-x)) mod l]);
        if ld = 1 then
            omicron := Agrp[1];
        else
            for a in Agrp do
                if a[3] = 1 + p then
                    omicron := a;
                    break;
                fi;
            od;
        fi;

        Char := function(i, j)
            local Chi;
            Chi := function(x)
                return E(OrderMod(alpha, l))^(x[1]*i);
            end;
            return Chi;
        end;
    else
        # A = <-1> * <5>.
        ord := OrderMod(5, l);
        Aind := Cartesian([0 .. 1], [0 .. ord - 1]);
        Agrp := List(Aind, x ->
                [x[1], x[2], (((-1)^x[1]) * (5^x[2])) mod l,
                (((-1)^x[1]) * (5^(ord-x[2]))) mod l]);
        for a in Agrp do
            if a[3] = 5 then
                omicron := a;
                break;
            fi;
        od;
        Char := function(i, j)
            local Chi;
            Chi := function(x)
                return (-1)^(x[1]*i) * E(ord)^(x[2]*j);
            end;
            return Chi;
        end;
    fi;

    # Find basis.
    # Note: for p^ld = 4 or 8, the resulting rep. is reducible,
    # so we must perform a change of basis later.  This basis covers
    # the resulting subspaces.
    Bp := [];
    B1 := [];
    for a in [0 .. l-1] do
        if Gcd(a,l) = 1 then
            Append(Bp, [[a,1]]);
        else
            Append(B1, [[a,1], [1,a]]);
        fi;
    od;
    Bp := Concatenation(Bp, B1);

    return [Agrp, omicron, Char, Bp];
end;

#-------------------------------------------------------------
# Representation of type D. Input: p, ld and the label of the character.
# Index is of the form [i,j], but note that j is only relevant for p = 2, ld > 2.
#-------------------------------------------------------------

RepD := function(p, ld, chi_index, silent)
    local i, j, l, M, Agrp, omicron, Chi, Bp, sxy, S, T, deg, w, U;

    if p <= 3 and ld = 1 then
        if not silent then
            Print("Invalid level.");
        fi;
        return fail;
    fi;

    l := p^ld;
    M := MD(p, ld);
    Agrp := M[1];
    omicron := M[2];

    Chi := M[3](chi_index[1], chi_index[2]);

    Bp := M[4];

    # Check for primitivity.  Primitive if chi is injective on <omicron>.
    if OrderMod(omicron[3], l) =
            Length(AsSet(List([0..Length(Agrp)-1], x -> (Chi(omicron))^x))) then
        if not silent then
            Print("chi is primitive.\n");
        fi;

        S := List(Bp, x -> List(Bp, y ->
                (1/l) * Sum(Agrp, a ->
                Chi(a) * E(l)^(a[3] * x[2] * y[1] + a[4] * x[1] * y[2]))));
        T := DiagonalMat(List(Bp, x -> E(l)^(x[1] * x[2])));

        if p = 2 and (ld = 2 or ld = 3) then
            # Convert to block diagonal.
            w := 1 / Sqrt(2);
            if ld = 2 then
                deg := 3;
                U := [
                    [ 1, 0, 0, 0, 0, 0],
                    [ 0, 0, 0, 1, 0, 0],
                    [ 0, w, 0, 0, w, 0],
                    [ 0, w, 0, 0,-w, 0],
                    [ 0, 0, w, 0, 0, w],
                    [ 0, 0, w, 0, 0,-w],
                ];
            else
                deg := 6;
                if Chi([1,0]) = 1 then
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, w, 0, 0, 0, 0, 0, w, 0, 0, 0],
                        [ 0, 0, w, 0, 0, 0, 0, 0,-w, 0, 0, 0],
                        [ 0, 0, 0, w, 0, 0, 0, 0, 0, w, 0, 0],
                        [ 0, 0, 0, w, 0, 0, 0, 0, 0,-w, 0, 0],
                        [ 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w, 0],
                        [ 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                else
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, w, 0, 0, 0, 0, 0, w, 0, 0, 0],
                        [ 0, 0, w, 0, 0, 0, 0, 0,-w, 0, 0, 0],
                        [ 0, 0, 0, w, 0, 0, 0, 0, 0, w, 0, 0],
                        [ 0, 0, 0, w, 0, 0, 0, 0, 0,-w, 0, 0],
                        [ 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w, 0],
                        [ 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                fi;
            fi;
            S := Inverse(U) * S * U;
            T := Inverse(U) * T * U;

            if not silent then
                Print("The character is primitive of order 2 and the representation is reducible. It decomposes into two irreducible components. The output is of the form [D_", ld, "(chi)_+, D_", ld, "(chi)_-].\n\n");
            fi;

            return [
                [S{[1..deg]}{[1..deg]}, T{[1..deg]}{[1..deg]}, deg],
                [S{[deg+1..deg*2]}{[deg+1..deg*2]}, T{[deg+1..deg*2]}{[deg+1..deg*2]}, deg]
            ];
        else
            deg := Length(Bp);

            if not silent then
                Print("chi^2 != 1.\n");
                Print("The dimension of the representation is ", deg, ".\n");
                Print("\n");
            fi;

            return [S, T, deg];
        fi;
    else
        if not silent then
            Print("chi is NOT primitive.\n");
        fi;
        return fail;
    fi;
end;
