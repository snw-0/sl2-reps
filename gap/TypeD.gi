#
# SL2Reps: Constructs representations of SL2(Z).
#
# Representations of type D.
#
# Implementations
#

InstallGlobalFunction( SL2Reps_ModuleD,
function(p, ld)
    local l, M, alpha, ord, omicron, Agrp, Bp, B1, Char, IsPrim, a;

    if not IsPrime(p) then
        Error("p must be prime.");
    elif not ld in PositiveIntegers then
        Error("ld must be a positive integer.");
        # TODO: technically ld = 0 is fine, it just gives us the trivial rep.
    elif p <= 3 and ld = 1 then
        Error("Type D is not defined for p <= 3 and ld = 1.");
    fi;

    l := p^ld;

    M := Tuples([0 .. (l - 1)], 2);

    # For type D, A = A_ld^times.
    if p > 2 or ld < 3 then
        # A is cyclic.
        alpha := GeneratorsPrimeResidues(l).generators[1];
        ord := OrderMod(alpha, l);
        Agrp := List([0 .. ord - 1], x ->
                [
                    [x], # Power of alpha.
                    (alpha^x) mod l, # a.
                    (alpha^(ord-x)) mod l # a^(-1), used for the S matrix.
                ]);
        if ld = 1 then
            omicron := Agrp[1];
        else
            for a in Agrp do
                if a[2] = 1 + p then
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
        Agrp := List(Cartesian([0 .. 1], [0 .. ord - 1]), x ->
                [
                    [x[1], x[2]], # Powers of -1 and 5.
                    (((-1)^x[1]) * (5^x[2])) mod l, # a.
                    (((-1)^x[1]) * (5^(ord-x[2]))) mod l # a^(-1), used for the S matrix.
                ]);
        for a in Agrp do
            if a[2] = 5 then
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

    # A character is primitive iff injective on omicron.
    IsPrim := function(chi)
        return OrderMod(omicron[2], l) = Order(chi(omicron[1]));
    end;

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

    return rec(
        Agrp := Agrp,
        Char := Char,
        IsPrim := IsPrim,
        Bp := Bp,
    );
end );

InstallGlobalFunction( SL2Reps_RepD,
function(p, ld, chi_index)
    local i, j, l, M_rec, Agrp, Chi, IsPrim, Bp, sxy, S, T, deg, w, U, a, b, k;

    if (not chi_index[1] in Integers) or (not chi_index[2] in Integers) then
        Error("chi must have integer indices.");
        # n.b.: it's fine if they're negative or 0
    fi;

    # this will check if p,ld are valid
    M_rec := SL2Reps_ModuleD(p, ld);

    l := p^ld;

    Agrp := M_rec.Agrp;
    Chi := M_rec.Char(chi_index[1], chi_index[2]);
    IsPrim := M_rec.IsPrim;
    Bp := M_rec.Bp;

    # Check for primitivity.  Primitive if chi is injective on <omicron>.
    if not IsPrim(Chi) then
        Error("chi is NOT primitive.");
        return fail;
        # TODO: is returning fail the correct action here?
    else
        Info(InfoSL2Reps, 2, "SL2Reps : chi is primitive.");

        S := List(Bp, x -> List(Bp, y ->
                (1/l) * Sum(Agrp, a ->
                Chi(a[1]) * E(l)^(a[2] * x[1] * y[2] + a[3] * x[2] * y[1]))));
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

            Info(InfoSL2Reps, 2, "SL2Reps : chi has order 2, so the representation is reducible. It decomposes into two irreducible components of degree ", deg, ". The output is of the form [D_", ld, "(chi)_+, D_", ld, "(chi)_-].");

            return [
                [S{[1..deg]}{[1..deg]}, T{[1..deg]}{[1..deg]}],
                [S{[deg+1..deg*2]}{[deg+1..deg*2]}, T{[deg+1..deg*2]}{[deg+1..deg*2]}]
            ];
        else
            deg := Length(Bp);

            Info(InfoSL2Reps, 2, "SL2Reps : chi has order greater than 2, so the representation is irreducible. The dimension of the representation is ", deg, ".");

            # Construct change of basis to make S symmetric.
            U := IdentityMat(Length(Bp));
            k := 1;
            repeat
                b := Bp[k];
                a := First(Agrp, x -> x[2] = b[1]);
                if a = fail then
                    # b[1] is not invertible, which means we have two basis elements,
                    # [b[1],1] and [1,b[1]]; they should be paired up.
                    U{[k,k+1]}{[k,k+1]} := (1 / Sqrt(2)) * [
                        [1, E(4)],
                        [1, -E(4)],
                    ];
                    k := k + 2;
                else
                    # b[1] is invertible, so we scale by Sqrt(Chi(a)).
                    U[k][k] := 1 / _SL2Reps_SqrtOfRootOfUnity(Chi(a[1]));
                    k := k + 1;
                fi;
            until k > Length(Bp);

            S := S ^ U;
            T := T ^ U;

            return [[S, T]];
        fi;
    fi;
end );
