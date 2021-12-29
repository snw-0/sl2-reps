#
# SL2Reps: Constructing symmetric representations of SL(2,Z).
#
# Representations of type D.
#
# Implementations
#


InstallGlobalFunction( SL2ModuleD,
function(p, ld)
    local l, M, alpha, beta, ord, Agrp, Bp, B1, Char, IsPrim, a;

    # TODO: add Nm and Prod functions, even though we don't need them

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
    if p > 2 then
        # p is odd.

        # A is cyclic of order p^(ld-1)(p-1). We use generators alpha = p+1 (of order p^(ld-1))
        # and beta (of order p-1; we search for this, as there seems to be no known closed form).
        # Note that if ld = 1 then alpha = 1 and beta is the only non-trivial generator.
        alpha := (p + 1) mod l;
        ord := p^(ld-1);
        beta := First([1 .. l], x -> OrderMod(x, l) = p-1);
        Agrp := List(Cartesian([0 .. ord - 1], [0 .. p-2]), x ->
                [
                    [x[1], x[2]], # Powers of alpha and beta.
                    ((alpha^x[1]) * (beta^x[2])) mod l, # a.
                    ((alpha^x[1]) * (beta^x[2]))^(-1) mod l # a^(-1), used for the S matrix.
                ]);

        Char := function(i, j)
            local Chi;
            Chi := function(x)
                return E(ord)^(x[1]*i) * E(p-1)^(x[2]*j);
            end;
            return Chi;
        end;

        if ld = 1 then
            # alpha is trivial. Hence, a character is primitive as long as it's not trivial.
            IsPrim := function(chi)
                return chi([0,1]) <> 1;
            end;
        else
            # A character is primitive iff injective on omega = alpha.
            IsPrim := function(chi)
                return OrderMod(alpha, l) = Order(chi([1,0]));
            end;
        fi;
    else
        # p = 2. Note: not defined when ld=1.
        if ld = 2 then
            # A = <-1>.
            Agrp := [
                [ [0,0], 1, 1 ], [ [1,0], 3, 3 ]
            ];
            Char := function(i, j)
                # j is irrelevant
                local Chi;
                Chi := function(x)
                    return (-1)^(x[1]*i);
                end;
                return Chi;
            end;

            # A character is primitive as long as it's non-trivial <=> chi(-1) != 1.
            IsPrim := function(chi)
                return chi([1,0]) <> 1;
            end;
        else
            # A = <5> * <-1>.
            ord := OrderMod(5, l);
            Agrp := List(Cartesian([0 .. ord - 1], [0 .. 1]), x ->
                    [
                        [x[1], x[2]], # Powers of 5 and -1.
                        (((5)^x[1]) * ((-1)^x[2])) mod l, # a.
                        (((5)^x[1]) * ((-1)^x[2]))^(-1) mod l # a^(-1), used for the S matrix.
                    ]);
            # for a in Agrp do
            #     if a[2] = 5 then
            #         omega := a;
            #         break;
            #     fi;
            # od;
            Char := function(i, j)
                local Chi;
                Chi := function(x)
                    return E(ord)^(x[1]*i) * (-1)^(x[2]*j);
                end;
                return Chi;
            end;

            # A character is primitive iff injective on omega = alpha = 5.
            IsPrim := function(chi)
                return OrderMod(5, l) = Order(chi([1,0]));
            end;
        fi;
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

    return rec(
        Agrp := Agrp,
        Char := Char,
        IsPrim := IsPrim,
        Bp := Bp,
    );
end );

InstallGlobalFunction( SL2IrrepD,
function(p, ld, chi_index)
    local i, j, l, M_rec, Agrp, Chi, IsPrim, Bp, sxy, S, T, deg, w, U, a, b, k;

    if (not chi_index[1] in Integers) or (not chi_index[2] in Integers) then
        Error("chi must have integer indices.");
        # n.b.: it's fine if they're negative or 0
    fi;

    # this will check if p,ld are valid
    M_rec := SL2ModuleD(p, ld);

    l := p^ld;

    Agrp := M_rec.Agrp;
    Chi := M_rec.Char(chi_index[1], chi_index[2]);
    IsPrim := M_rec.IsPrim;
    Bp := M_rec.Bp;

    # Check for primitivity.  Primitive if chi is injective on <omega>.
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
                # check if Chi is injective on <-1>
                if Chi([0,1]) = 1 then
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
                    U[k][k] := 1 / _SL2SqrtOfRootOfUnity(Chi(a[1]));
                    k := k + 1;
                fi;
            until k > Length(Bp);

            S := S ^ U;
            T := T ^ U;

            return [[S, T]];
        fi;
    fi;
end );
