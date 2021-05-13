#----------------------------------------------------------------------------------------
# Output: the module M and basic functions on it. 
#----------------------------------------------------------------------------------------
MRBasic := function(p, ld, si, r, t)
    local l, ls, m1, m2, M, tM,
        Prod, Pow, Ord, Nm, c, A;

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

    # Scale factor, used for the S matrix: S_Q(-1) / Sqrt(|M|) = (1/|M|) * Sum(e(Q(x))).
    c := Sum(M, x -> E(l)^(r * Nm(x))) / Length(M);

    # The abelian subgroup A of Aut(M, Q).
    A := Filtered(M, x -> Nm(x) = 1);

    return rec(
        M := M, 
        tM := tM, 
        Prod := Prod, 
        Pow := Pow, 
        Ord := Ord,
        Nm := Nm, 
        c := c,
        A := A);
end;

#----------------------------------------------------------------------------------------
# Output: when p > 2, the group A, characters on A and primitivity criterion of the characters.
# First find two generators alpha and zeta of the group A (which can be 1 in some cases). Then write every element a in M as a = alpha^i * zeta^j, and record the element as [[i, j], a]. Note that a is also a list with two entries. 
#----------------------------------------------------------------------------------------
AgrpOdd := function(p, ld, si, r, t)
    local basic, l, ls, m1, m2, M, tM, 
        Prod, Pow, Ord, Nm,
        A, A01, Aind, Agrp, alpha, zeta, zeta_coords,
        Char, omicron, IsPrim,
        a, b;
    
    basic := MRBasic(p, ld, si, r, t);
    M := basic.M;
    tM := basic.tM;
    Pow := basic.Pow;
    Prod := basic.Prod;
    Ord := basic.Ord;
    Nm := basic.Nm;
    A := basic.A;

    l := p^ld;
    ls := p^si;
    m1 := p^(ld);
    m2 := p^(ld-si);

    # Find generators for A and construct Agrp.
    if p = 3 and ld >= 3 and si = 1 and t = 1 then
        # Special case:
        # A is not cyclic; instead, A = <alpha> x <zeta>
        # where Ord(alpha) = 3^(ld-2) and Ord(zeta) = 6.
        #
        # alpha and -zeta are both found in the group A_0 - A_1.
        #
        # A character is primitive iff it is injective on <alpha>.

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


    # A character is primitive iff injective on omicron.
    IsPrim := function(chi)
        return Ord(omicron[2]) = Order(chi(omicron[1]));
    end;

    # Return.
    return rec(
        Agrp := Agrp,
        Char := Char,
        IsPrim := IsPrim
    );
end;

#----------------------------------------------------------------------------------------
# Output: when p = 2, the group A, characters on A and primitivity criterion of the characters.
# First find two generators alpha and zeta of the group A (which can be 1 in some cases). Then write every element a in M as a = alpha^i * zeta^j, and record the element as [i, j, [a]]. Note that [a] is a list with two entries. 
#----------------------------------------------------------------------------------------
AgrpEven := function(p, ld, si, r, t)
    local basic, l, ls, m1, m2, M, tM, pM,
            Prod, Pow, Ord, Nm,
            t_rep, r_rep,
            A, A01, Aind, Agrp, alpha, zeta,
            Char, omicron, IsPrim, c,
            a, b;

    basic := MRBasic(p, ld, si, r, t);
    M := basic.M;
    tM := basic.tM;
    Pow := basic.Pow;
    Prod := basic.Prod;
    Nm := basic.Nm;
    Ord := basic.Ord;
    A := basic.A;

    l := p^ld;
    ls := p^si;
    m1 := p^(ld-1);
    m2 := p^(ld-si-1);


    # Find generators for A and construct Agrp.    
    if si = ld - 2 then
        if ld = 3 then
            Print("This type is the same as TypeN.\n");
            return;
        elif ld = 4 then
            alpha := [1,0];
            zeta := [7,0];
            Agrp := [[[0,0], [1, 0]], [[0,1], [7,0]]];
            omicron := [[1,0], alpha];

            Char := function(i, j)
                    local Chi;
                    Chi := function(x)
                            return (-1)^(x[2]*j);
                    end;
                    return Chi;
            end;

            IsPrim := function(chi)
                return Ord(omicron[2]) = Order(chi(omicron[1]));
            end;
            return rec(
                Agrp := Agrp,
                Char := Char,
                IsPrim := IsPrim,
            );
        else
            # A = <alpha> x <-1> with ord(alpha) = 2.
            #
            # This case is an exception to the usual definition of a primitive char:
            # here there are NO non-trivial elements of A that fix pM pointwise, so
            # we instead call a character primitive if it is injective on alpha
            # (i.e., if chi(alpha) = -1).

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
        # Use the standard representative for (r,t) as per table on p. 475, NW part I.
        # Then apply table on p. 496, NW part II to find generators of A.
        t_rep := t mod Minimum(8, 2^(ld-si));
        if si = 0 then
            r_rep := Minimum(r, r*t_rep mod 4);
            if r in [1,3] and t = 1 then
                if ld = 3 then
                    alpha := [1,0];
                    zeta := [0,1];
                    omicron := [[0,2], [-1 mod m1, 0]];
                else
                    alpha := First(A, x -> (x[1] mod 4 = 1) and x[2] = 4);
                    zeta := [0,1];
                    omicron := [[1,0], alpha];
                fi;
            elif r in [1,3] and t = 5 then
                if ld = 3 then
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
            elif r = 1 and t in [3,7] then
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

    # A character is primitive iff injective on omicron.
    IsPrim := function(chi)
        return Ord(omicron[2]) = Order(chi(omicron[1]));
    end;

    # Return.
    return rec(
        Agrp := Agrp,
        Char := Char,
        IsPrim := IsPrim,
    );
end;

#-----------------------------------------------------------------------------
# Output: module and Agrp data that will be used in RepRs.
#-----------------------------------------------------------------------------
MR := function(p, ld, si, r, t)
    local basic, c, Nm, Ord, Pow, Prod, tM,
        fA, Agrp, Char, IsPrim;


    basic := MRBasic(p, ld, si, r, t);
    c := basic.c;
    Nm := basic.Nm;
    Ord := basic.Ord;
    Pow := basic.Pow;
    Prod := basic.Prod;
    tM := basic.tM;

    if p = 2 then
        fA := AgrpEven(p, ld, si, r, t);
    else 
        fA := AgrpOdd(p, ld, si, r, t);
    fi;

    Agrp := fA.Agrp;
    Char := fA.Char;
    IsPrim := fA.IsPrim;

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

#-----------------------------------------------------------------------------
# Output: square root of roots of unity.
#-----------------------------------------------------------------------------

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