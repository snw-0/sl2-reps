#-------------------------------------------------------------
# MR gives the module and character information of type R_{lambda}^{sigma} for most cases.
#-------------------------------------------------------------

MR := function(p, ld, si, r, t)
    local l, ls, m1, m2, M, tM, pM,
            Prod, Pow, Ord, Nm,
            t_rep, r_rep,
            A, A01, Aind, Agrp, alpha, zeta, zeta_coords, AOrbit,
            Char, omicron, IsPrim, c,
            a, b;

    l := p^ld;
    ls := p^si;

    if p = 2 then
        m1 := p^(ld-1);
        m2 := p^(ld-si-1);
    else
        m1 := p^(ld);
        m2 := p^(ld-si);
    fi;

    # Module.
    M := Cartesian([0 .. (m1 - 1)], [0 .. (m2 - 1)]);
    tM := Filtered(M, x -> x[1] mod p <> 0 or x[2] mod p <> 0);
    pM := Filtered(M, x -> x[1] mod p = 0 and x[2] mod p = 0);

    Prod := function(x, y)
        return [((x[1] * y[1]) - (x[2] * y[2] * ls * t)) mod m1,
                ((x[1] * y[2]) + (y[1] * x[2])) mod m2];
    end;

    Pow := function(x, n)
        local y, i;
        y := [x[1] mod m1, x[2] mod m2];
        if n > 1 then
            for i in [1..n-1] do
                y := Prod(y, x);
            od;
        elif n = 0 then
            y := [1,0];
        fi;
        return y;
    end;

    Ord := function(x)
        local i, y;
        i := 1;
        y := x;
        while y <> [1,0] do
            y := Prod(x, y);
            i := i + 1;
        od;
        return i;
    end;

    Nm := function(x)
        return ((x[1]^2) + ls * t * (x[2]^2)) mod l;
    end;

    # Scale factor, used for the S matrix:
    # S_Q(-1) / Sqrt(|M|) = (1/|M|) * Sum(e(Q(x))).
    c := Sum(M, x -> E(l)^(r * Nm(x))) / Length(M);

    # Construct A.
    A := Filtered(M, x -> Nm(x) = 1);

    # Find generators for A and construct Agrp.
    if p = 2 then
        if si = ld - 2 then
            # NW S. 2.5 and 6.3.
            if ld = 3 then
                Print("This type is the same as TypeN.\n");
                return;
            elif ld = 4 then
                alpha := [1,0];
                zeta := [7,0];
                Agrp := [[0,0, [1, 0]], [0,1, [7,0]]];
                omicron := [[1,0], alpha];

                Char := function(i, j)
                        local Chi;
                        Chi := function(x)
                                return (-1)^(x[2]*j);
                        end;
                        return Chi;
                end;

                # XXXXX refactor away this separate return point
                IsPrim := function(chi)
                    return Ord(omicron[2]) = Order(chi(omicron[1]));
                end;
                return rec(
                    Agrp := Agrp,
                    Char := Char,
                    IsPrim := IsPrim,
                    tM := tM,
                    Nm := Nm,
                    Prod := Prod,
                    Ord := Ord,
                    c := c
                );
            else
                # A = <alpha> x <-1> with ord(alpha) = 2.
                #
                # This case is an exception to the usual definition of a primitive char.:
                # here there are NO non-trivial elements of A that fix pM pointwise, so
                # we instead call a character primitive if it is injective on alpha
                # (i.e. if chi(alpha) = -1).

                if ld = 5 then
                    alpha := [(1 + 4 * t) mod m1, 1];
                else
                    alpha := [(1 - (2^(ld-3) * t)) mod m1, 1];
                fi;
                zeta := [-1 mod m1, 0];
                omicron := [[1,0], alpha];
            fi;

            Aind := Cartesian([0 .. 1], [0 .. 1]);
            Agrp := List(Aind, x -> [
                    [x[1], x[2]],
                    Prod(Pow(alpha, x[1]), Pow(zeta, x[2]))
                    ]);

            Char := function(i, j)
                local Chi;
                Chi := function(x)
                        return (-1)^(x[1]*i) * (-1)^(x[2]*j);
                end;
                return Chi;
            end;
        else
            # NW S. 2.4.
            #
            # First, we find the standard representative for (r,t),
            # as per table on p. 475, NW part I.
            #
            # Then, we apply table on p. 496, NW part II to find generators of A.
            # In all cases except ld=3, si=0, t=5 (for which see below)
            # we have A = <alpha> x <zeta> with zeta = [-1,0] or [0,1].
            #
            # A character is primitive if injective on:
            # <-1>, when ld = 3 and si = 0,
            # <-alpha^2>, when ld = 4, si = 0, t = 5,
            # <alpha>, otherwise.
            t_rep := t mod Minimum(8, 2^(ld-si));
            if si = 0 then
                r_rep := Minimum(r, (r*t_rep) mod 4);
                if r_rep in [1,3] and t_rep = 1 then
                    if ld = 3 then
                        # Unique case: ord(alpha) = 1.
                        alpha := [1,0];
                        zeta := [0,1];
                        omicron := [[0,2], [-1 mod m1, 0]];
                    else
                        # Note: NW give alpha = (4, 1 mod 4) here.
                        #
                        # This is incorrect for ld = 4,5, where we instead need (1 mod 4, 4);
                        # for ld >= 6 it is irrelevant which of the two we use, as they both
                        # have the same order and both generate the primitive element,
                        # [1, 2^(ld-2)], in the same way.
                        #
                        # We therefore use the latter throughout.
                        alpha := First(A, x -> (x[1] mod 4 = 1) and x[2] = 4);
                        zeta := [0,1];
                        omicron := [[1,0], alpha];
                    fi;
                elif r_rep in [1,3] and t_rep = 5 then
                    if ld = 3 then
                        # Unique case: see footnote NW p. 496.
                        alpha := [1,0];
                        zeta := [2,1];
                        omicron := [[0,2], [-1 mod m1, 0]];
                    else
                        alpha := First(A, x -> x[1] = 2 and (x[2] mod 4) = 3);
                        zeta := [-1 mod m1, 0];
                        if ld = 4 then
                            omicron := [[2,1],
                                    Prod([-1 mod m1, 0], Pow(alpha, 2))];
                        else
                            omicron := [[1,0], alpha];
                        fi;
                    fi;
                elif r_rep = 1 and t_rep in [3,7] then
                    if ld = 3 then
                        alpha := [1,0];
                        zeta := [-1 mod m1, 0];
                        omicron := [[0,1], [-1 mod m1, 0]];
                    else
                        alpha := First(A, x -> (x[1] mod 4 = 1) and x[2] = 4);
                        zeta := [-1 mod m1, 0];
                        omicron := [[1,0], alpha];
                    fi;
                fi;
            elif si = 1 then
                r_rep := Minimum(r, (r + 2*r*t_rep) mod 8);
                # There are two cases here (namely [[1,5],[1,5]] and [[1,3],[3,7]])
                # but for this question they differ only in which alpha is selected,
                # not the criterion therefore.
                alpha := First(A, x -> (x[1] mod 4 = 1) and x[2] = 2);
                zeta := [-1 mod m1, 0];
                omicron := [[1,0], alpha];
            elif si = 2 then
                r_rep := r mod 4;
                alpha := First(A, x -> (x[1] mod 4 = 1) and x[2] = 2);
                zeta := [-1 mod m1, 0];
                omicron := [[1,0], alpha];
            else
                r_rep := r mod 8;
                alpha := First(A, x -> (x[1] mod 4 = 1) and x[2] = 1);
                zeta := [-1 mod m1, 0];
                omicron := [[1,0], alpha];
            fi;

            Aind := Cartesian([0 .. Ord(alpha) - 1], [0 .. Ord(zeta) - 1]);
            Agrp := List(Aind, x -> [
                    [x[1], x[2]],
                    Prod(Pow(alpha, x[1]), Pow(zeta, x[2]))
                    ]);

            Char := function(i, j)
                local Chi;
                Chi := function(x)
                    return E(Ord(alpha))^(x[1]*i) * E(Ord(zeta))^(x[2]*j);
                end;
                return Chi;
            end;
        fi;
    else
        if p = 3 and ld >= 3 and si = 1 and t = 1 then
            # Special case:
            # A is not cyclic; instead, A = <alpha> x <zeta>
            # where Ord(alpha) = 3^(ld-2) and Ord(zeta) = 6.
            #
            # alpha and -zeta are both found in the group A_0 - A_1.
            #
            # pM is fixed pointwise by the subgroup generated by
            # alpha^(3^(ld-3)) = [1, 3^(ld-2)].
            # Hence, a character is primitive iff it is injective on <alpha>.

            if ld = 3 then
                # Unique case; alpha = alpha^(3^(ld-3)).
                alpha := [1,3];
                zeta := [23,7];
            else
                zeta_coords := List([1,3,5],
                        x -> (1 + 3 * (1/2) * ((x * 3^(ld-2)) - 1)) mod m1);
                A01 := Filtered(A, x -> (x[1] mod (ls*t) = 1) and (x[2] mod p <> 0));
                alpha := First(A01, x -> not x[1] in zeta_coords);
                zeta := Prod([-1 mod m1, 0], First(A01, x -> x[1] in zeta_coords));
            fi;
            omicron := [[1,0], alpha];

            Aind := Cartesian([0 .. (3^(ld-2)) - 1], [0 .. 5]);
            Agrp := List(Aind, x -> [
                    [x[1], x[2]],
                    Prod(Pow(alpha, x[1]), Pow(zeta, x[2]))
                    ]);

            Char := function(i, j)
                local Chi;
                Chi := function(x)
                    return E(3^(ld-2))^(x[1]*i) * E(6)^(x[2]*j);
                end;
                return Chi;
            end;
        else
            # A is cyclic; A = <alpha> x <-1> where Ord(alpha) = p^(ld-si).
            # A character is primitive iff it is injective on <alpha>.

            alpha := First(A, x -> Ord(x) = p^(ld-si));
            omicron := [[1,0], alpha];

            Aind := Cartesian([0 .. (p^(ld-si)) - 1], [0 .. 1]);
            Agrp := List(Aind, x -> [
                    [x[1], x[2]],
                    Prod(Pow(alpha, x[1]), [(-1)^(x[2]) mod m1, 0])
                    ]);

            Char := function(i, j)
                local Chi;
                Chi := function(x)
                    return E(p^(ld-si))^(x[1]*i) * (-1)^(x[2]*j);
                end;
                return Chi;
            end;
        fi;
    fi;

    # A character is primitive iff injective on omicron.
    IsPrim := function(chi)
        return Ord(omicron[2]) = Order(chi(omicron[1]));
    end;

    # Return.
    return rec(
        Agrp := Agrp,
        Char := Char,
        IsPrim := IsPrim,
        tM := tM,
        Nm := Nm,
        Prod := Prod,
        Ord := Ord,
        c := c
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
# Representation of type R_{lambda}^{sigma}. Input: p, ld, sigma, r, t, and [i,j] the label of the character.
#-------------------------------------------------------------

RepR := function(p, ld, si, r, t, chi_index, silent)
    local l, M_rec, Agrp, Chi, IsPrim, tM, Nm, Prod, Ord, c,
            Tr, AOrbit, Epsilon,
            tM1, theta, b, B1, Bp, BaseChangeMat, w, U, U_index, B_Q, sxy, S, T, deg,
            N, B, O, VInd, tO, a, k;

    l := p^ld;

    M_rec := MR(p, ld, si, r, t);

    Agrp := M_rec.Agrp;
    Chi := M_rec.Char(chi_index[1], chi_index[2]);
    IsPrim := M_rec.IsPrim;
    tM := M_rec.tM;
    Nm := M_rec.Nm;
    Prod := M_rec.Prod;
    Ord := M_rec.Ord;
    c := M_rec.c;

    Tr := function(x)
        return (2 * x[1]) mod l;
    end;

    AOrbit := function(x)
        return List(Agrp, a -> Prod(a[2],x));
    end;

    # A character is primitive if it is injective on <omicron>;
    # usually omicron = alpha, but not for all cases.
    if IsPrim(Chi) then
        if not silent then
            Print("chi is primitive.\n");
        fi;

        # Find basis for primitive characters; this can depend on the character,
        # so we have to do it here instead of in MR.
        if p = 2 then
            if ld >= 3 and si <= ld - 3 then
                theta := function(a)
                    # find theta_a, with norm a
                    local eta, x;
                    if a = 1 then
                        eta := [1, 0];
                    else
                        eta := First(tM, x -> Nm(x) mod l = a);
                        if eta = fail then
                            Error("No eta found with norm ", a);
                        fi;
                    fi;
                    return List(List([0 .. (2^(ld-3) - 1)], x -> (2*x + 1)),
                            x -> Prod([x, 0], eta));
                end;

                # see table in NW sec. 5.
                if si = 0 then
                    if t mod 4 = 1 then
                        Bp := Concatenation(theta(1), theta(5));
                    else
                        Bp := Concatenation(theta(1), theta(3), theta(5), theta(7));
                    fi;
                elif si = 1 then
                    if t mod 4 = 1 then
                        Bp := Concatenation(theta(1), theta(3));
                    else
                        Bp := Concatenation(theta(1), theta(7));
                    fi;
                elif si = 2 then
                    Bp := Concatenation(theta(1), theta(5));
                    U := IdentityMat(3*2^(ld-3));
                else
                    Bp := theta(1);
                fi;

                # now find representatives for the remaining A orbits
                tM1 := ShallowCopy(tM);
                B1 := [];
                for b in Bp do
                    SubtractSet(tM1, AOrbit(b));
                od;
                while Length(tM1) > 0 do
                    b := tM1[1];
                    # Does kappa(b) lie in the same A-orbit as b?
                    a := First(Agrp, x -> Prod(x[2],b) = [b[1], -b[2] mod p^(ld-si-1)]);
                    if a = fail then
                        # No, it doesn't. We should therefore pair them up at the end.
                        Add(B1,b);
                        SubtractSet(tM1, AOrbit(b));
                        Add(B1,[b[1], -b[2] mod p^(ld-si-1)]);
                        SubtractSet(tM1, AOrbit([b[1], -b[2] mod p^(ld-si-1)]));
                    else
                        # Yes, it does. We will need to scale by 1 / Sqrt(Chi(a)) later.
                        Add(Bp,b);
                        SubtractSet(tM1, AOrbit(b));
                    fi;
                od;

                Bp := Concatenation(Bp, B1);
            elif ld >= 5 and si = ld - 2 then
                # Depends on the character; see NW p. 511.
                Bp := [];

                # theta_1
                b := 1;
                while b <= 2^(ld-2) - 1 do
                    Add(Bp, [b, 0]);
                    b := b + 2;
                od;

                # theta_2
                b := 2;
                while b <= 2^(ld-3) - 2 do
                    Add(Bp, [b, 0]);
                    b := b + 4;
                od;

                # theta_3
                if ld > 5 then
                    b := 4;
                    while b <= 2^(ld-3) - 4 do
                        Add(Bp, [b, 1]);
                        b := b + 4;
                    od;
                fi;

                # The (single) basis element in theta_4 or theta_5 corresponds
                # to an A-orbit of size 2 (instead of |A| = 4), so we must normalize it
                # below by dividing the corresponding entries of the S matrix by sqrt(2).
                # Note that this isn't done the same way as other base changes - here the
                # problem is that our basis is not orthonormal, so we're correcting for that.
                # We could also use a double sum as with the non-primitive cases, but this
                # is simpler.
                if chi_index[2] mod 2 = 0 then
                    # chi_{-1}
                    # theta_4
                    Add(Bp, [0,1]);
                else
                    # chi_{-alpha}
                    # theta_5
                    Add(Bp, [2^(ld-3), 1]);
                fi;

                # Separate return point in order to normalize as described above.
                deg := Length(Bp);

                sxy := function(x, y)
                    return c * Sum(Agrp, a ->
                            Chi(a[1]) * E(l)^(r * Tr(Prod(a[2], Prod(x, [y[1], -y[2]])))));
                end;

                S := List(Bp, x -> List(Bp, y -> sxy(x, y)));
                for b in [1 .. deg] do
                    S[deg][b] := S[deg][b] / Sqrt(2);
                    S[b][deg] := S[b][deg] / Sqrt(2);
                od;
                T := DiagonalMat(List(Bp, x -> E(l)^(r * Nm(x))));

                return [[S, T]];
            elif ld = 4 and si = 2 then
                N := Cartesian([0 .. 7], [0 .. 1]);
                B := [];
                O := [];
                while Length(N) > 0 do
                    Add(B, N[1]);
                    tO := Set(Agrp, a -> Prod(a[3], N[1]));
                    Add(O, tO);
                    SubtractSet(N, tO);
                od;

                VInd := [];
                for k in [1 .. Length(B)] do
                    for a in Agrp do
                        if a <> [0, 0, [1, 0]] and Chi(a) <> 1 and Prod(a[3], B[k]) = B[k] then
                            Add(VInd, k);
                            break;
                        fi;
                    od;
                od;

                SubtractSet(B, List(VInd, k -> B[k]));
                SubtractSet(O, List(VInd, k -> O[k]));

                deg := Length(B);

                B_Q := function(a,b)
                    # B((x,y),(w,z)) = (2r/p^ld) * (xw + (p^si)tyz)
                    return 2 * r * (a[1] * b[1] + (2^si) * t * a[2] * b[2]);
                end;

                sxy := function(x, y)
                    return c * Sum(Agrp, a ->
                            Sum(Agrp, b ->
                            Chi(a) * ComplexConjugate(Chi(b)) * E(l)^(B_Q(Prod(a[3], x), Prod(b[3], y)))));
                end;

                S := List([1..deg], x ->
                        List([1..deg], y ->
                        Sqrt(Length(O[x]) * Length(O[y])) * sxy(B[x], B[y]) / (Length(Agrp))^2));

                T := DiagonalMat(List(B, x -> E(l)^(r*Nm(x))));
                return [[S, T]];
            elif ld = 3 and si = 1 then
                if not silent then
                    Print("This type is the same as TypeN, use TypeN code.\n");
                fi;
                return fail;
            else
                if not silent then
                    Print("Not regular!");
                fi;
                return fail;
            fi;
        else
            # theta_1 = theta intersect M^times
            Bp := List(Filtered(PrimeResidues(l),
                    a -> a >= 1 and a <= (l-1)/2),
                    a -> [a,0]);

            # now find representatives for the remaining A orbits
            tM1 := ShallowCopy(tM);
            B1 := [];
            for b in Bp do
                SubtractSet(tM1, AOrbit(b));
            od;
            while Length(tM1) > 0 do
                b := tM1[1];
                # Does kappa(b) lie in the same A-orbit as b?
                a := First(Agrp, x -> Prod(x[2],b) = [b[1], -b[2] mod p^(ld-si)]);
                if a = fail then
                    # No, it doesn't. We should therefore pair them up.
                    Add(B1,b);
                    SubtractSet(tM1, AOrbit(b));
                    Add(B1,[b[1], -b[2] mod p^(ld-si)]);
                    SubtractSet(tM1, AOrbit([b[1], -b[2] mod p^(ld-si)]));
                else
                    # Yes, it does. So we need to scale by Chi(a).
                    Add(Bp,b);
                    SubtractSet(tM1, AOrbit(b));
                fi;
            od;

            Bp := Concatenation(Bp, B1);
        fi;

        sxy := function(x, y)
            return c * Sum(Agrp, a ->
                    Chi(a[1]) * E(l)^(r * Tr(Prod(a[2], Prod(y, [x[1], -x[2]])))));
        end;

        S := List(Bp, x -> List(Bp, y -> sxy(x, y)));
        T := DiagonalMat(List(Bp, x -> E(l)^(r * Nm(x))));
        deg := Length(Bp);

        # Handle +- cases, where chi is primitive but squares to 1.
        if p=2 and ld=3 and si=0 and r=1 and t in [3,7] and chi_index=[0,1] then
            w := 1 / Sqrt(2);
            if t = 3 then
                U := [
                    [ 1,  0,  0,  0,  0,  0],
                    [ 0,  0,  0,  1,  0,  0],
                    [ 0,  1,  0,  0,  0,  0],
                    [ 0,  0,  0,  0,  1,  0],
                    [ 0,  0,  w,  0,  0,  w],
                    [ 0,  0,  w,  0,  0, -w],
                ];
            else
                # basis is in a different order
                U := [
                    [ 1,  0,  0,  0,  0,  0],
                    [ 0,  0,  0,  0,  1,  0],
                    [ 0,  1,  0,  0,  0,  0],
                    [ 0,  0,  0,  1,  0,  0],
                    [ 0,  0,  w,  0,  0,  w],
                    [ 0,  0,  w,  0,  0, -w],
                ];
            fi;
            deg := 3;

            S := Inverse(U) * S * U;
            T := Inverse(U) * T * U;

            if not silent then
                Print("The character is primitive of order 2 and the representation is reducible. It decomposes into two irreducible components.\nThe output is of the form [R_3^0(1,", t, ",chi)_+, R_3^0(1,", t, ",chi)_-].\n\n");
            fi;

            return [
                [S{[1..deg]}{[1..deg]}, T{[1..deg]}{[1..deg]}],
                [S{[deg+1..deg*2]}{[deg+1..deg*2]}, T{[deg+1..deg*2]}{[deg+1..deg*2]}]
            ];
        elif p=2 and ld=4 and si=0 and r in [1,3] and ((t=1 and chi_index in [[1,0],[1,2]]) or (t=5 and chi_index in [[0,1],[2,1]])) then
            deg := 3;
            if not silent then
                Print("The character is primitive of order 2 and the representation is reducible. It decomposes into two irreducible components.\nThe output is of the form [R_4^0(", r, ",", t, ",chi)_+, R_4^0(", r, ",", t, ",chi)_-].\n\n");
            fi;

            # Here the subreps. are conveniently generated by sub-bases
            if (t = 1 and chi_index = [1,0]) or (t = 5 and chi_index = [0,1]) then
                return [
                    [S{[1,2,5]}{[1,2,5]}, T{[1,2,5]}{[1,2,5]}],
                    [S{[3,4,6]}{[3,4,6]}, T{[3,4,6]}{[3,4,6]}]
                ];
            else
                return [
                    [S{[1,2,6]}{[1,2,6]}, T{[1,2,6]}{[1,2,6]}],
                    [S{[3,4,5]}{[3,4,5]}, T{[3,4,5]}{[3,4,5]}]
                ];
            fi;
        elif p=2 and ld=4 and si=0 and r=1 and t in [3,7] and chi_index in [[1,0],[1,1]] then
            # See notes for bases. Note that they are in different orders for t = 3,7.
            w := 1 / Sqrt(2);
            if t = 3 then
                if chi_index = [1,0] then
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
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
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w, 0],
                        [ 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                fi;
            else
                if chi_index = [1,0] then
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w, 0],
                        [ 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                else
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w, 0],
                        [ 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                fi;
            fi;
            deg := 6;

            S := Inverse(U) * S * U;
            T := Inverse(U) * T * U;

            if not silent then
                Print("The character is primitive of order 2 and the representation is reducible. It decomposes into two irreducible components.\nThe output is of the form [R_4^0(1,", t, ",chi)_+, R_4^0(1,", t, ",chi)_-].\n\n");
            fi;

            return [
                [S{[1..deg]}{[1..deg]}, T{[1..deg]}{[1..deg]}],
                [S{[deg+1..deg*2]}{[deg+1..deg*2]}, T{[deg+1..deg*2]}{[deg+1..deg*2]}]
            ];
        elif p=2 and ld=5 and si=2 and r in [1,3] and t in [1,3,5,7] and chi_index in [[1,0],[1,1]] then
            # See notes for bases. Note that they are in different orders for t = 1,3,5,7.
            w := 1 / Sqrt(2);
            if t = 1 then
                if chi_index = [1,0] then
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                else
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                fi;
            elif t = 3 then
                if chi_index = [1,0] then
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                else
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                fi;
            elif t = 5 then
                if chi_index = [1,0] then
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                else
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                fi;
            else
                if chi_index = [1,0] then
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                else
                    U := [
                        [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0, w],
                        [ 0, 0, 0, 0, 0, w, 0, 0, 0, 0, 0,-w],
                    ];
                fi;
            fi;
            deg := 6;

            S := Inverse(U) * S * U;
            T := Inverse(U) * T * U;

            if not silent then
                Print("The character is primitive of order 2 and the representation is reducible. It decomposes into two irreducible components.\nThe output is of the form [R_5,2(", r, ",", t, ",chi)_+, R_5,2(", r, ",", t, ",chi)_-].\n\n");
            fi;

            return [
                [S{[1..deg]}{[1..deg]}, T{[1..deg]}{[1..deg]}],
                [S{[deg+1..deg*2]}{[deg+1..deg*2]}, T{[deg+1..deg*2]}{[deg+1..deg*2]}]
            ];
        else
            # Apply a change of basis to make S symmetric.
            U := IdentityMat(Length(Bp));
            k := 1;
            repeat
                b := Bp[k];
                # Is kappa(b) in the same A-orbit as b?
                if p = 2 then
                    a := First(Agrp, x -> Prod(x[2],b) = [b[1], -b[2] mod p^(ld-si-1)]);
                else
                    a := First(Agrp, x -> Prod(x[2],b) = [b[1], -b[2] mod p^(ld-si)]);
                fi;
                if a = fail then
                    # No, it isn't. This and the next basis element should be paired.
                    U{[k, k+1]}{[k, k+1]} := (1 / Sqrt(2)) * [
                        [1, E(4)],
                        [1, -E(4)]
                    ];
                    k := k + 2;
                else
                    # Yes, it is.  We should scale by 1 / Sqrt(Chi(a)).
                    U[k][k] := 1 / SqrtOfRootOfUnity(Chi(a[1]));
                    k := k + 1;
                fi;
            until k > Length(Bp);

            S := S ^ U;
            T := T ^ U;

            return [[S, T]];
        fi;
    else
        if not silent then
            Print("chi is NOT primitive.\n");
        fi;

        if p = 2 and ((ld >= 5 and si = ld - 2)
                or (ld = 5 and si = 2 and r in [1,3] and t = 1)
                or (ld >= 6 and si = ld - 3)) then
            # Handle cases of non-primitive characters chi such that V(chi) still has some
            # irred. subrep. of the same level as the main rep.
            if not silent then
                Print("There is a single irreducible subrepresentation of the same level as V.\n");
            fi;

            # Function to construct a particular type of base change matrix used several times below
            BaseChangeMat := function(len, i, sign)
                local U, b;
                U := NullMat(len, len);
                for b in [1..i-1] do
                    U[b][b] := 1;
                od;
                b := 0;
                while i + 2*b < len do
                    U[i + 2*b][i + b] := 1 / Sqrt(2);
                    U[i + 1 + 2*b][i + b] := sign / Sqrt(2);

                    U[i + 2*b][i + (len - i + 1)/2 + b] := 1 / Sqrt(2);
                    U[i + 1 + 2*b][i + (len - i + 1)/2 + b] := -1 * sign / Sqrt(2);

                    b := b + 1;
                od;
                return U;
            end;

            # We need to construct a basis for a larger representation, and then cut down
            # to the (unique) irreducible subrep. by applying a base-change matrix.
            if ld >= 5 and si = ld - 2 then
                # See NW S. 6.3.

                # Depends on the character; see NW p. 511. Note that |A| = 4.
                Bp := [];

                # Construct a basis that spans the subrep. we need.
                if chi_index[2] mod 2 = 0 then
                    # chi_1
                    deg := 0;

                    # theta_1. These elements are of the usual form, with A-orbit size |A| = 4.
                    b := 1;
                    while b <= 2^(ld-2) - 1 do
                        Add(Bp, [[b, 0], 4]);
                        deg := deg + 1;
                        b := b + 2;
                    od;

                    # Remaining basis elements are of the form f_{xi} - f_{2^(ld-2) - xi};
                    # we therefore add both summands to Bp and use a base-change matrix later.
                    U_index := Length(Bp) + 1; # Index where these compound elements begin.
                    # These have A-orbits of length 1.
                    Add(Bp, [[0, 0], 1]);
                    Add(Bp, [[2^(ld-2), 0], 1]);
                    deg := deg + 1;
                    # The remainder have A-orbits of length 2.
                    b := 4;
                    while b <= 2^(ld-3) - 4 do
                        Add(Bp, [[b, 0], 2]);
                        Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 0], 2]);
                        deg := deg + 1;
                        b := b + 4;
                    od;
                    b := 2;
                    while b <= 2^(ld-3) - 2 do
                        Add(Bp, [[b, 1], 2]);
                        Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 1], 2]);
                        deg := deg + 1;
                        b := b + 4;
                    od;

                    # Now construct the base-change matrix.
                    U := BaseChangeMat(Length(Bp), U_index, -1);
                else
                    # chi_{alpha}
                    deg := 0;

                    # theta_1. These elements are of the usual form, with A-orbit size |A| = 4.
                    b := 1;
                    while b <= 2^(ld-2) - 1 do
                        Add(Bp, [[b, 0], 4]);
                        deg := deg + 1;
                        b := b + 2;
                    od;

                    # This element corresponds to an A-orbit of length 2.
                    Add(Bp, [[2^(ld-3), 0], 2]);
                    deg := deg + 1;

                    # Remaining basis elements are of the form f_{xi} + f_{2^(ld-2) - xi};
                    # we therefore add both summands to Bp and use a base-change matrix later.
                    # All these basis elements have A-orbits of length 2.
                    U_index := Length(Bp) + 1; # Index where these compound elements begin.
                    b := 2;
                    while b <= 2^(ld-3) - 2 do
                        Add(Bp, [[b, 1], 2]);
                        Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 1], 2]);
                        deg := deg + 1;
                        b := b + 4;
                    od;
                    if ld > 5 then
                        b := 4;
                        while b <= 2^(ld-3) - 4 do
                            Add(Bp, [[b, 0], 2]);
                            Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 0], 2]);
                            deg := deg + 1;
                            b := b + 4;
                        od;
                    fi;

                    # Now construct the base-change matrix.
                    U := BaseChangeMat(Length(Bp), U_index, 1);
                fi;
            elif ld = 5 and si = 2 and r in [1,3] and t = 1 then
                # Depends on the character. See NW. p. 524 for basis.
                Bp := [];

                # A = <alpha> x <-1> with ord(alpha) = 2.
                # Here we are handling the non-primitive characters, namely [0,0] and [0,1].

                deg := 12;

                # These basis elements are of the usual form, with A-orbit size |A| = 4.
                for b in [1,3,5,7] do
                    Add(Bp, [[b,0], 4]);
                od;
                for b in [1,3,5,7] do
                    Add(Bp, [[b,1], 4]);
                od;

                if chi_index = [0,0] then
                    # chi_1 = nu.

                    # The remaining basis elements are of the form f_{xi} - f_{8 - xi}.
                    # These have A-orbit length 2.
                    Add(Bp, [[2,0], 2]);
                    Add(Bp, [[6,0], 2]);

                    Add(Bp, [[2,2], 2]);
                    Add(Bp, [[6,2], 2]);

                    # These have A-orbit length 1 (i.e. are fixed by A).
                    Add(Bp, [[0,0], 1]);
                    Add(Bp, [[8,0], 1]);

                    Add(Bp, [[0,2], 1]);
                    Add(Bp, [[8,2], 1]);

                    # Construct a base-change matrix to extract the appropriate subrep. (dim 12)
                    # from the dim. 16 space spanned by the above basis.
                    w := 1 / Sqrt(2);
                    U := BaseChangeMat(16,9,-1);
                else
                    # chi_alpha, indexed by [0,1].

                    # All of these basis elements have A-orbit length 2.
                    # These two are still of the form f_{xi}.
                    Add(Bp, [[4,0], 2]);

                    Add(Bp, [[4,2], 2]);

                    # The remainder are of the form f_{xi} + f_{8 - xi}.
                    Add(Bp, [[2,0], 2]);
                    Add(Bp, [[6,0], 2]);

                    Add(Bp, [[2,2], 2]);
                    Add(Bp, [[6,2], 2]);

                    # Construct a base-change matrix to extract the appropriate subrep. (dim 12)
                    # from the dim. 14 space spanned by the above basis.
                    w := 1 / Sqrt(2);
                    U := BaseChangeMat(14,11,1);
                fi;
            elif ld >= 6 and si = ld - 3 then
                # See NW S. 6.4.

                # Depends on the character; see NW p. 512. Note that |A| = 8.
                Bp := [];

                # Construct a basis that spans the subrep. we need.
                # Ord(alpha) = 4, and chars are primitive if Chi(alpha) = E(4) or E(4)^3.
                # So the cases we are currently handling have Chi(alpha) = +-1.
                if chi_index[1] mod 4 = 0 then
                    # chi(alpha) = 1
                    if chi_index[2] mod 2 = 0 then
                        # chi_1 = nu.
                        deg := 0;

                        # theta_1. These elements are of the usual form, with A-orbit size |A| = 8.
                        b := 1;
                        while b <= 2^(ld-2) - 1 do
                            Add(Bp, [[b, 0], 8]);
                            deg := deg + 1;
                            b := b + 2;
                        od;

                        # Remaining basis elements are of the form f_{xi} - f_{2^(ld-2) - xi};
                        # we therefore add both summands to Bp and use a base-change matrix later.
                        U_index := Length(Bp) + 1; # Index where these compound elements begin.
                        # These have A-orbits of length 4.
                        b := 2;
                        while b <= 2^(ld-3) - 2 do
                            Add(Bp, [[b, 0], 4]);
                            Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 0], 4]);
                            deg := deg + 1;
                            b := b + 4;
                        od;
                        # These have A-orbits of length 2.
                        b := 4;
                        while b <= 2^(ld-3) - 4 do
                            Add(Bp, [[b, 2], 2]);
                            Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 2], 2]);
                            deg := deg + 1;
                            b := b + 8;
                        od;
                        # This pair has A-orbits of length 1.
                        Add(Bp, [[0, 0], 1]);
                        Add(Bp, [[2^(ld-2), 0], 1]);
                        deg := deg + 1;
                        # These have A-orbits of length 2.
                        b := 8;
                        while b <= 2^(ld-3) - 8 do
                            Add(Bp, [[b, 0], 2]);
                            Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 0], 2]);
                            deg := deg + 1;
                            b := b + 8;
                        od;

                        # Now construct the base-change matrix.
                        U := BaseChangeMat(Length(Bp), U_index, -1);
                    else
                        # chi_{alpha}
                        deg := 0;

                        # theta_1. These elements are of the usual form, with A-orbit size |A| = 8.
                        b := 1;
                        while b <= 2^(ld-2) - 1 do
                            Add(Bp, [[b, 0], 8]);
                            deg := deg + 1;
                            b := b + 2;
                        od;

                        # SEE REMARK BELOW. This has A-orbit size 2.
                        Add(Bp, [[2^(ld-3), 0], 2]);
                        deg := deg + 1;

                        # Remaining basis elements are of the form f_{xi} + f_{2^(ld-2) - xi};
                        # we therefore add both summands to Bp and use a base-change matrix later.
                        U_index := Length(Bp) + 1; # Index where these compound elements begin.
                        # These have A-orbits of length 4.
                        b := 2;
                        while b <= 2^(ld-3) - 2 do
                            Add(Bp, [[b, 0], 4]);
                            Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 0], 4]);
                            deg := deg + 1;
                            b := b + 4;
                        od;
                        # These have A-orbits of length 2.
                        b := 4;
                        while b <= 2^(ld-3) - 4 do
                            Add(Bp, [[b, 2], 2]);
                            Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 2], 2]);
                            deg := deg + 1;
                            b := b + 8;
                        od;
                        # These have A-orbits of length 2.
                        # REMARK: NW include b = 2^(ld-3) in this sequence, but then we get
                        # 2f_{2^(ld-3)} (instead of sum of two distinct fs).
                        # Moving it up as a single basis element for now.
                        b := 8;
                        while b <= 2^(ld-3) - 8 do
                            Add(Bp, [[b, 0], 2]);
                            Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 0], 2]);
                            deg := deg + 1;
                            b := b + 8;
                        od;

                        # Now construct the base-change matrix.
                        U := BaseChangeMat(Length(Bp), U_index, 1);
                    fi;
                else
                    # chi(alpha) = -1
                    if chi_index[2] mod 2 = 0 then
                        # chi_{-1}
                        deg := 0;

                        # theta_1. These elements are of the usual form, with A-orbit size |A| = 8.
                        b := 1;
                        while b <= 2^(ld-2) - 1 do
                            Add(Bp, [[b, 0], 8]);
                            deg := deg + 1;
                            b := b + 2;
                        od;

                        # These elements have A-orbits of length 4.
                        b := 4;
                        while b <= 2^(ld-3) - 4 do
                            Add(Bp, [[b, 0], 4]);
                            deg := deg + 1;
                            b := b + 8;
                        od;

                        # This element has A-orbit of length 2.
                        Add(Bp, [[0, 2], 2]);
                        deg := deg + 1;

                        # These elements have A-orbits of length 4.
                        b := 8;
                        while b <= 2^(ld-3) - 8 do
                            Add(Bp, [[b, 2], 4]);
                            deg := deg + 1;
                            b := b + 8;
                        od;

                        # Remaining basis elements are of the form f_{xi} - f_{2^(ld-2) - xi};
                        # we therefore add both summands to Bp and use a base-change matrix later.
                        U_index := Length(Bp) + 1; # Index where these compound elements begin.
                        # These have A-orbits of length 4.
                        b := 2;
                        while b <= 2^(ld-3) - 2 do
                            Add(Bp, [[b, 0], 4]);
                            Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 0], 4]);
                            deg := deg + 1;
                            b := b + 4;
                        od;

                        # Now construct the base-change matrix.
                        U := BaseChangeMat(Length(Bp), U_index, -1);
                    else
                        # chi_{-alpha}
                        deg := 0;

                        # theta_1. These elements are of the usual form, with A-orbit size |A| = 8.
                        b := 1;
                        while b <= 2^(ld-2) - 1 do
                            Add(Bp, [[b, 0], 8]);
                            deg := deg + 1;
                            b := b + 2;
                        od;

                        # These elements have A-orbits of length 4.
                        b := 4;
                        while b <= 2^(ld-3) - 4 do
                            Add(Bp, [[b, 0], 4]);
                            deg := deg + 1;
                            b := b + 8;
                        od;

                        # These elements have A-orbits of length 4.
                        b := 8;
                        while b <= 2^(ld-3) - 8 do
                            Add(Bp, [[b, 2], 4]);
                            deg := deg + 1;
                            b := b + 8;
                        od;

                        # This element has A-orbit of length 2.
                        Add(Bp, [[2^(ld-3), 2], 2]);
                        deg := deg + 1;

                        # Remaining basis elements are of the form f_{xi} + f_{2^(ld-2) - xi};
                        # we therefore add both summands to Bp and use a base-change matrix later.
                        U_index := Length(Bp) + 1; # Index where these compound elements begin.
                        # The remainder have A-orbits of length 4.
                        b := 2;
                        while b <= 2^(ld-3) - 2 do
                            Add(Bp, [[b, 0], 4]);
                            Add(Bp, [[(2^(ld-2) - b) mod 2^(ld-1), 0], 4]);
                            deg := deg + 1;
                            b := b + 4;
                        od;

                        # Now construct the base-change matrix.
                        U := BaseChangeMat(Length(Bp), U_index, 1);
                    fi;
                fi;
            fi;

            # Here we must use the double sum to find S_{x,y}: for the standard formula
            # used earlier, we assumed that both basis elements correspond to A-orbits
            # with length |A|; that is guaranteed for prim. characters, but not here.
            B_Q := function(a,b)
                # B((x,y),(w,z)) = (2r/p^ld) * (xw + (p^si)tyz)
                return 2 * r * (a[1] * b[1] + (2^si) * t * a[2] * b[2]);
            end;
            sxy := function(x, y)
                return (Sqrt(x[2] * y[2]) / (Length(Agrp)^2)) * c * Sum(Agrp, a -> Sum(Agrp, b ->
                        Chi(a[1]) * ComplexConjugate(Chi(b[1]))
                        * (E(l)^B_Q(Prod(a[2], x[1]), Prod(b[2], y[1])))));
            end;

            S := List(Bp, x -> List(Bp, y -> sxy(x, y)));
            T := DiagonalMat(List(Bp, x -> E(l)^(r * Nm(x[1]))));

            # Now we change basis and take the irreducible subrep.
            S := S ^ U;
            T := T ^ U;

            return [[S{[1..deg]}{[1..deg]}, T{[1..deg]}{[1..deg]}]];
        else
            if not silent then
                Print("V_chi is reducible and has no subreps of full level.\n");
                return fail;
            fi;
        fi;
    fi;
end;
