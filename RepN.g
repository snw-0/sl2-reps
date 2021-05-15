#-------------------------------------------------------------
# MN gives the module and character information of type N.
#
# The basis for a type N representation for a primitive chi depends only on p and ld,
# so it can also be produced here.
#-------------------------------------------------------------

MN := function(p, ld)
    local l, t, M, pM, tM, A, Nm, Prod, Pow, Ord, alpha, zeta, Agrp, Char, i, j, u, theta, Bp;

    l := p^ld;

    # Find t.
    if p = 2 then t := 3;
    else
        t := First([0..4*l], i -> (i mod 4 = 3) and (Jacobi(-i, p) = -1));
    fi;

    # Construct the module.
    M := Tuples([0..l-1], 2);
    tM := Tuples([0..l-1], 2);
    pM := Tuples(List([0..l/p-1], x -> p*x mod l), 2);
    SubtractSet(tM, pM);

    # Product, power and order.
    Prod := function(a, x)
        return [(a[1]*x[1] - (1 + t)*a[2]*x[2]/4) mod l, (a[1]*x[2]+a[2]*x[1]+a[2]*x[2]) mod l];
    end;

    Pow := function(a, n)
        local y, i;
        y := [a[1] mod l, a[2] mod l];
        if n > 1 then
            for i in [1..n-1] do
                y := Prod(y, a);
            od;
        elif n = 0 then
            y := [1,0];
        fi;
        return y;
    end;

    Ord := function(a)
        local x, b;
        x := 1;
        b := a;
        while b <> [1,0] do
            b := Prod(a, b);
            x := x + 1;
        od;
        return x;
    end;

    # Norm and the group A of norm 1.
    Nm := function(a)
        return (a[1]^2 + a[1]*a[2] + (1 + t)*a[2]^2/4) mod l;
    end;

    A := Filtered(M, a -> Nm(a) = 1);

    # Use Ord(zeta) = p+1 or 6 to find zeta.
    if ld = 1 or ((p > 2) and (ld > 1)) then
        zeta := First(A, x -> Ord(x) = p+1);
    else
        zeta := First(A, x -> Ord(x) = 6);
    fi;

    # Use Ord(alpha) = 2^(ld - 2) or p^(ld - 1) to find alpha.
    if ld = 1 then
        alpha := [1,0];
    elif p = 2 then
        alpha := First(A, x -> (Ord(x) = 2^(ld-2)) and (x[1] mod 4 = 1) and (x[2] mod 4 = 0));
    else
        alpha := First(A, x -> Ord(x) = p^(ld-1) and (x[1] mod p = 1) and (x[2] mod p = 0));
    fi;

    # Use powers of alpha and zeta to index elements in A.
    Agrp := List(Cartesian([0..Ord(alpha)-1], [0..Ord(zeta)-1]), x ->
            [
                [x[1], x[2]],
                Prod(Pow(alpha,x[1]), Pow(zeta, x[2]))
            ]);

    # Use (i, j) to index chars. Each chi(i, j) is a function Agrp -> C^*.
    Char := function(i, j)
        local Chi;
        Chi := function(x)
            return E(Ord(alpha))^(x[1]*i) * E(Ord(zeta))^(x[2]*j);
        end;
        return Chi;
    end;

    # Find the bases for primitive chars, which depend only on p and ld.
    # For odd p, pick u to be the smallest quad nonres mod p (any quad nonres works).
    if p > 2 then
        u := First([0..p-1], x -> Jacobi(x, p) = -1);
    fi;

    # Find the sets theta. Only works for p = odd or p = 2 and ld >= 2.
    theta := function(a)
        local i, eta, t1;
        eta := [1, 0];
        if a > 1 then
            eta := First(tM, x -> Nm(x) = a);
        fi;
        if p > 2 then
            t1 := Filtered([1..(l-1)/2], n -> Gcd(n, p) = 1);
        elif p = 2 and ld > 2 then
            t1 := Filtered([1..2^(ld-2)-1], n -> n mod 2 = 1);
        elif p = 2 and ld = 2 then
            t1 := [1];
        fi;
        return List(t1, x -> Prod([x, 0], eta));
    end;

    # The basis for primitive chars.
    if p = 2 and ld = 1 then
        Bp := [[1, 0]];
    elif p = 2 and ld = 2 then
        Bp := [[1, 0], theta(3)[1]];
    elif p = 2 and ld > 2 then
        Bp := Concatenation(theta(1), theta(3), theta(5), theta(7));
    else
        Bp := Concatenation(theta(1), theta(u));
    fi;

    # Return.
    return rec(
        Agrp := Agrp,
        Char := Char,
        Bp := Bp,
        Nm := Nm,
        Prod := Prod
    );
end;

SqrtOfRootOfUnity := function(b)
    local desc;

    if b = 1 then
        return 1;
    else
        # this returns [q,p] such that b = E(q)^p
        desc := DescriptionOfRootOfUnity(b);

        return E(desc[1]*2)^desc[2];
    fi;
end;

#-------------------------------------------------------------
# Representation of type N. Input: p, ld and [i, j], where [i, j] labels the character.
#-------------------------------------------------------------

RepN := function(p, ld, chi_index, silent)
    local l, M_rec, Agrp, Chi, Bp, beta, Nm, Prod, Tr, sxy, S, T, deg,
            N, B, O, tO, BQ, a, j, k, VInd, Prim1, Prim2, Prim3, U;

    l := p^ld;

    M_rec := MN(p, ld);

    Agrp := M_rec.Agrp;
    Chi := M_rec.Char(chi_index[1], chi_index[2]);
    Bp := M_rec.Bp;
    Nm := M_rec.Nm;
    Prod := M_rec.Prod;

    Tr := function(x)
        return (2*x[1] + x[2]) mod l;
    end;

    BQ := function(x, y)
        return (Nm([x[1]+y[1], x[2]+y[2]]) - Nm(x) - Nm(y)) mod l;
    end;

    Prim1 := (((p = 2) and (ld > 2)) or ((p > 2) and (ld > 1))) and (Gcd(chi_index[1], p) = 1);

    Prim2 := (p = 2) and (ld = 2) and (Chi([ 0, 3 ]) = -1);

	Prim3 := (p>2) and (ld = 1) and ((chi_index[2] mod (p+1)/2) <> 0);

    if Prim1 or Prim2 or Prim3 then
        if not silent then
            Print("Chi is primitive.\n");

            if Length(AsSet(List(Agrp, x -> Chi(x[1])^2))) > 1 then
                Print("Chi^2 != 1.\n");
            fi;
        fi;

        sxy := function(x, y)
            local z;
            z := Prod(x, [y[1] + y[2], -y[2]]);
            return ((-1)^ld / l) * Sum(Agrp, a -> Chi(a[1]) * E(l)^(Tr(Prod(a[2], z))));
        end;

        S := List(Bp, x -> List(Bp, y -> sxy(x, y)));
        T := DiagonalMat(List(Bp, x -> E(l)^(Nm(x))));

        if not (p = 2 and ld = 1) then
            # We need to perform a change of basis to make S symmetric.

            # beta is the element of Agrp equal to x / x-bar. This is used to make the s-matrix symmetric below.
            beta := function(x)
                local b, n;

                n := 1 / Nm(x) mod l;
                b := Prod(Prod(x,x), [n,0]); # x / x-bar
                b := First(Agrp, y -> y[2] = b);

                return b;
            end;

            U := [];

            if p = 2 and ld > 2 then
                # Bp consists of theta_1, theta_3, theta_5, theta_7.
                for j in [1 .. Length(Bp) / 4] do
                    Add(U, 1);
                od;
                for j in [1 .. Length(Bp) / 4] do
                    Add(U, SqrtOfRootOfUnity(1 / Chi(beta(Bp[1 + Length(Bp)/4])[1])));
                od;
                for j in [1 .. Length(Bp) / 4] do
                    Add(U, SqrtOfRootOfUnity(1 / Chi(beta(Bp[1 + 2 * Length(Bp)/4])[1])));
                od;
                for j in [1 .. Length(Bp) / 4] do
                    Add(U, SqrtOfRootOfUnity(1 / Chi(beta(Bp[1 + 3 * Length(Bp)/4])[1])));
                od;
            else
                # For p = 2, ld = 2, Bp consists of theta_1 and theta_3.
                # Otherwise, Bp consists of theta_1 and theta_u, where u is a quad. non-residue.
                # The process we need here is the same, so we just handle both cases together.

                for j in [1 .. Length(Bp) / 2] do
                    Add(U, 1);
                od;
                for j in [1 .. Length(Bp) / 2] do
                    Add(U, SqrtOfRootOfUnity(1 / Chi(beta(Bp[1 + Length(Bp)/2])[1])));
                od;
            fi;

            U := DiagonalMat(U);
            S := S ^ U;
            T := T ^ U;
        fi;

        deg := Length(Bp);
    elif (ld = 1 and Length(AsSet(List(Agrp, x -> Chi(x[1])))) = 1) then
        if not silent then
            Print("ld = 1, and Chi is the trivial character. This representation is also called the Steinberg representation.\n");
        fi;

        Bp := Concatenation([[0,0]], Bp);

        sxy := function(x, y)
            local z;
            if x = [0,0] and y = [0,0] then
                return -1 / p;
            elif x = [0,0] or y = [0,0] then
                return -Sqrt(p+1) / p;
            else
                z := Prod(x, [y[1] + y[2], -y[2]]);
                return ((-1)^ld / l) * Sum(Agrp, a -> Chi(a[1])*E(l)^(Tr(Prod(a[2], z))));
            fi;
        end;

        S := List(Bp, x -> List(Bp, y -> sxy(x, y)));
        T := DiagonalMat(List(Bp, x -> E(l)^(Nm(x))));
        deg := Length(Bp);
    else
        if not silent then
            Print("Chi is not primitive or Ord(Chi) <= 2. It is also not the Steinberg representation. The first decomposition method gives the following representation corresponding to Chi that is REDUCIBLE.\n");
        fi;

        N := Tuples([0..l-1], 2);;
        B := [];
        O := [];
        while Length(N) > 0 do
            Add(B,N[1]);
            tO := Set(Agrp, a -> Prod(a[2], N[1]));
            Add(O, tO);
            SubtractSet(N, tO);
        od;

        VInd := [];
        for k in [1..Length(B)] do
            for a in Agrp do
                if a[1] <> [0, 0] and Chi(a[1]) <> 1 and Prod(a[2], B[k]) = B[k] then
                    Add(VInd, k); break;
                fi;
            od;
        od;

        SubtractSet(B, List(VInd, k -> B[k]));
        SubtractSet(O, List(VInd, k -> O[k]));

        deg := Length(B);
        sxy := function(x, y)
            return ((-1)^ld / l) * Sum(Agrp, a ->
                    Sum(Agrp, b ->
                    Chi(a[1]) * ComplexConjugate(Chi(b[1])) * E(l)^(BQ(Prod(a[2], x), Prod(b[2], y)))));
        end;

        S := List([1..deg], x ->
                List([1..deg], y ->
                Sqrt(Length(O[x]) * Length(O[y])) * sxy(B[x], B[y]) / (Length(Agrp))^2));

        T := DiagonalMat(List(B, x -> E(l)^(Nm(x))));
    fi;

    if [p, ld, chi_index[1] mod 2, chi_index[2] mod 6] = [2, 3, 1, 0] then
        if not silent then
            Print("Special case [2, 3, 1, 0]. The character is primitive of order 2 and the representation is reducible. It decomposes into two irreducible components. The output is of the form [N_3(chi)_+, N_3(chi)_-].\n");
        fi;

        return [
            [S{[1,2]}{[1,2]}, T{[1,2]}{[1,2]}],
            [S{[3,4]}{[3,4]}, T{[3,4]}{[3,4]}]
        ];
    elif [p, ld, chi_index[1] mod 2, chi_index[2] mod 6] = [2, 3, 1, 3] then
        if not silent then
            Print("Special case [2, 3, 1, 3]. The character is primitive of order 2 and the representation is reducible. It decomposes into two irreducible components. The output is of the form [N_3(chi)_+, N_3(chi)_-].\n");
        fi;

        U := [
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [0, 1, 0, 0],
        ];
        S := S ^ U;
        T := T ^ U;

        return [
            [S{[1,2]}{[1,2]}, T{[1,2]}{[1,2]}],
            [S{[3,4]}{[3,4]}, T{[3,4]}{[3,4]}]
        ];
    else
        return [[S, T]];
    fi;
end;
