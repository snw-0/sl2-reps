Read("./RepN.g");
Read("./RepD.g");
Read("./RepRs_regular.g");
Read("./RepRu.g");
Read("./Conj_Class.g");

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Construct SL_2(Z/p^ld) with appropriate conjugacy classes for our use.
SL2Conj := function(p, ld)
    local l, CC, o, s, t, G, C1, c, ccl;

    l := p^ld;
    CC := CClasses(p, ld);
    o := ZmodnZObj(1, l);
    s := [[0,1],[-1,0]] * o;
    t := [[1,1],[0,1]] * o;
    G := Group([s,t]);

    C1 := [];
    if p > 2 then
        for c in CC[1] do
            if Length(c[1]) = 1 then
                Add(C1, t^(c[1][1])*s*(-1)^(c[2]));
            elif (Length(c[1])=4 and c[1][1]<>1) or Length(c[1]) = 2 then
                Add(C1, Product(c[1], m -> t^m*s)*(-1)^(c[2]));
            fi;
        od;
    else
        for c in CC[1] do
            if Length(c[1])=4 and c[1][1]<>1 then
                Add(C1, Product(c[1], m -> t^m*s)*(-1)^(c[2]));
            elif Length(c[1])=2 then
                Add(C1, Product(c[1], m -> t^m*s)*(-1)^(c[2]));
            elif Length(c[1])=3 then
                Add(C1, Product(c[1], m -> t^m*s)*(-1)^(c[2]));
            elif Length(c[1])=5 then
                Add(C1, Product(c[1], m -> t^m*s)*(-1)^(c[2]));
            elif Length(c[1])=1 then
                Add(C1, Product(c[1], m -> t^m*s)*(-1)^(c[2]));
            fi;
        od;
    fi;

    ccl := List(C1, x -> ConjugacyClass(G, x));
    if Sum(ccl, x-> Size(x)) <> Size(G) then
        Print("Wrong Conjugacy Classes\n");
    fi;

    SetConjugacyClasses( G, ccl );

    return G;
end;

# Construct chi according to the above conjugacy classes.
ChiST := function(S, T)
    local l, primepower, p, ld, CC, s2, TS, i, ProdTrace,
            C1, Du, c,
            Dh, h, e1, e1inv, e2, e2inv, j, jinv, Dhmap;

    # S,T are normalized S and T matrices
    l := Order(T);
    primepower:=Factors(l);
    if Length(AsSet(primepower))>1 then
        Print("The level is not a prime power. We cannot proceed!", "\n");
        return fail;
    fi;
    p := primepower[1]; ld := Length(primepower);
    CC := CClasses(p, ld);

    s2 := S^2;
    s2 := s2[1][1]; # assuming S^2 = +-1
    TS := [S];
    for i in [1 .. l-1] do
        Add(TS, T * TS[i]);
    od;

    ProdTrace := function(A, B)
        # assumes they are both square matrices of same size
        return Sum([1..Length(A)], x -> Sum([1..Length(A)], y -> A[x][y]*B[y][x]));
    end;

    C1:=[];
    if p > 2 then
        Du := TS[CC[2]+1] * TS[CC[3]+1] * TS[CC[2]+1];
        for c in CC[1] do
            if Length(c[1]) = 1 then
                Add(C1, Trace(TS[(c[1][1] mod l) +1] * s2^(c[2])));
            elif Length(c[1]) = 4 and c[1][1]<>1 then
                Add(C1, ProdTrace(Du,TS[(c[1][4] mod l) +1]) * s2^(c[2]));
            elif Length(c[1]) = 2 then
                Add(C1, ProdTrace(TS[(c[1][1] mod l)+1], TS[(c[1][2] mod l)+1]) * s2^(c[2]));
            fi;
        od;
    else
        Dh:=[];

        for h in [3..ld] do
            e1 := 1+2^(h-1);
            e1inv := e1^-1 mod l;
            Add(Dh, [e1inv,TS[e1inv+1] * TS[e1+1] * TS[e1inv+1]]);
            e2 := (-e1) mod l;
            e2inv := (-e1inv) mod l;
            Add(Dh, [e2inv, TS[e2inv+1] * TS[e2+1] * TS[e2inv+1]]);
        od;

        for i in [3,5,7] do
            j := i mod l;
            jinv:=j^-1 mod l;
            Add(Dh, [jinv,TS[jinv+1] * TS[j+1] * TS[jinv+1]]);
        od;
        Dh := AsSet(Dh);

        Dhmap:=function(uinv)
            local x;
            x := First(Dh, y -> y[1] = uinv);
            return x[2];
        end;

        for c in CC[1] do
            if Length(c[1])=4 and c[1][1]<>1 then
                Add(C1, ProdTrace(Dhmap(c[1][1]), TS[(c[1][4] mod l) +1]) * s2^(c[2]));
            elif Length(c[1])=2 then
                Add(C1, ProdTrace(TS[(c[1][1] mod l)+1], TS[(c[1][2] mod l)+1]) * s2^(c[2]));
            elif Length(c[1])=3 then
                Add(C1, Trace(Dhmap(c[1][1] mod l)) * s2^(c[2]));
            elif Length(c[1])=5 then
                Add(C1, ProdTrace((TS[(c[1][1] mod l)+1] * Dhmap(c[1][2] mod l)), TS[(c[1][5] mod l)+1]) * s2^(c[2]));
            elif Length(c[1])=1 then
                Add(C1, Trace(TS[(c[1][1] mod l)+1]) * s2^(c[2]));
            fi;
        od;
    fi;

    return C1;
end;

SomeQuadNonRes := function(p)
    local j;

    j := 2;
    while Jacobi(j, p) = 1 do
        j := j + 1;
    od;
    return j;
end;

RecordIrrep := function(irrep_list, name, rho, l)
    Print(name, " [level ", l, "]");
    # rho = [S, T, degree, level, name]
    Add(irrep_list, [rho[1], rho[2], rho[3], l, name]);
    Print("\n");
end;

_FactorizationsStep := function(factors)
    local c, m, o, rest, r;

    if Length(factors) <= 1 then
        return [factors];
    fi;

    o := [];

    for c in Combinations([1 .. Length(factors)]) do
        if Length(c) = 0 then
            continue;
        fi;

        m := Product(factors{c});

        rest := [1 .. Length(factors)];
        SubtractSet(rest, c);
        rest := factors{rest};

        for r in _FactorizationsStep(rest) do
            Add(r, m);
            Add(o, r);
        od;
    od;

    return o;
end;

Factorizations := function(n)
    local o, f;

    if n = 1 then
        return [[1]];
    fi;

    o := [];

    for f in _FactorizationsStep(Factors(n)) do
        Sort(f);
        f := Reversed(f);
        Add(f, 1);
        if not f in o then
            Add(o,f);
        fi;
    od;

    Sort(o);
    return Reversed(o);
end;

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Input: degree. Output: the number of all the irreps of dimension = degree.

PrimePowerIrrepsOfDegree := function(deg)
    local irrep_list, p, ld, pmax, ldmax, pset, ldset, name, rho, l,
            i, j, u, w, si, r, t, x;

    irrep_list := [];

    if 1 = deg then
        # The trivial irrep Xi_0 = C_1; level 1.
        name := "Xi_0";
        RecordIrrep(irrep_list, name, [[[1]], [[1]], 1], 1);
    fi;

    # p odd, ld = 1.
    pmax := 2*deg + 1;
    pset := Filtered([3..pmax], x -> IsPrime(x));
    for p in pset do
        l := p^1;
        if p + 1 = deg and p > 3 then
            # D_1(chi), for chi primitive, chi^2 != 1.
            # Defined for p >= 3.
            # A = <alpha>, ord(alpha) = p-1.
            # Relevant character indices: [(1..(p-3)/2), 0].
            for i in [1 .. (p-3)/2] do
                rho := RepD(p, 1, [i, 0], true);
                name := Concatenation("D_1([", String(i), ",0])");
                RecordIrrep(irrep_list, name, rho, l);
            od;
        fi;

        if p - 1 = deg then
            # N_1(chi), for chi primitive.
            # A = <zeta>, ord(zeta) = p+1.
            # Relevant character indices: [0, (1..(p-1)/2)].
            for j in [1 .. (p-1)/2] do
                rho := RepN(p, 1, [0, j], true);
                name := Concatenation("N_1([0,", String(j), "])");
                RecordIrrep(irrep_list, name, rho, l);
            od;
        fi;

        if (p+1)/2 = deg then
            # R_1(1)+ and R_1(u)+, for u a non-residue.
            # RepRUnary gives a list containing two reps, + and - .
            rho := RepRUnary(p, 1, 1, true);
            name := "R_1(1)+";
            RecordIrrep(irrep_list, name, rho[1], l);

            u := SomeQuadNonRes(p);
            rho := RepRUnary(p, 1, u, true);
            name := Concatenation("R_1(", String(u), ")+");
            RecordIrrep(irrep_list, name, rho[1], l);
        fi;

        if (p-1)/2 = deg then
            # R_1(1)- and R_1(u)-, for u a non-residue.
            # RepRUnary gives a list containing two reps, + and - .
            # For p=3, R_1(1)- = Xi_4 and R_1(2)- = Xi_8, the two linear reps. of level 3.
            rho := RepRUnary(p, 1, 1, true);
            name := "R_1(1)-";
            RecordIrrep(irrep_list, name, rho[2], l);

            u := SomeQuadNonRes(p);
            rho := RepRUnary(p, 1, u, true);
            name := Concatenation("R_1(", String(u), ")-");
            RecordIrrep(irrep_list, name, rho[2], l);
        fi;

        if p = deg then
            # N_1(nu), the Steinberg representation.
            rho := RepN(p, 1, [0, 0], true);
            name := "N_1(nu)";
            RecordIrrep(irrep_list, name, rho, l);
        fi;
    od;

    # p odd, ld >= 2.
    pmax := Maximum([RootInt(2*deg+1, 2), 3]);
    ldmax := Maximum([LogInt(deg, 3) + 2, 2]);
    pset := Filtered([3..pmax], x -> IsPrime(x));
    ldset := [2..ldmax];
    for p in pset do
        for ld in ldset do
            l := p^ld;
            if (p+1)*p^(ld-1) = deg then
                # D_ld(chi), for chi primitive, chi^2 != 1.
                # A = <alpha>, with ord(alpha) = p^ld - p^(ld-1).
                # chi is primitive iff it is injective on <1+p = omicron>; this occurs
                # when the index of chi is coprime to p.
                for i in [1 .. (l - l/p)/2 - 1] do
                    if Gcd(i, p) = 1 then
                        rho := RepD(p, ld, [i, 0], true);
                        name := Concatenation("D_", String(ld), "([", String(i), ",0])");
                        RecordIrrep(irrep_list, name, rho, l);
                    fi;
                od;
            fi;

            if (p-1)*p^(ld-1) = deg then
                # N_ld(chi), for chi primitive.
                # A = <alpha> x <zeta>; ord(alpha) = p^(ld-1), ord(zeta) = p+1.
                # chi is primitive if injective on <alpha>.
                # Remember chi and bar(chi) give same rep; for j = 0 or (p+1)/2,
                # bar(chi(i,j)) = chi(-i, j).
                for i in PrimeResidues(p^(ld-1)) do
                    if i > p^(ld-1) / 2 then
                        break;
                    fi;
                    rho := RepN(p, ld, [i, 0], true);
                    name := Concatenation("N_", String(ld), "([", String(i), ",0])");
                    RecordIrrep(irrep_list, name, rho, l);

                    rho := RepN(p, ld, [i, (p+1)/2], true);
                    name := Concatenation("N_", String(ld), "([", String(i), ",", String((p+1)/2), "])");
                    RecordIrrep(irrep_list, name, rho, l);
                od;
                for i in PrimeResidues(p^(ld-1)) do
                    for j in [1 .. (p-1)/2] do
                        rho := RepN(p, ld, [i, j], true);
                        name := Concatenation("N_", String(ld), "([", String(i), ",", String(j), "])");
                        RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            fi;

            if p^(ld-2)*(p^2-1)/2 = deg then
                # R_ld^si(r,t,chi), for 1 <= si <= ld-1, r,t either 1 or a non-residue, and chi primitive
                for si in [1 .. ld-1] do
                    u := SomeQuadNonRes(p);
                    for r in [1,u] do
                        for t in [1,u] do
                            if p = 3 and ld >= 3 and si = 1 and t = 1 then
                                # Special case, see RepRs.
                                # A = <alpha> x <zeta> where ord(alpha) = 3^(ld-2) and ord(zeta) = 6.
                                # chi is primitive if injective on <alpha>.
                                # Remember chi and bar(chi) give same rep; for j = 0 or 3,
                                # bar(chi(i,j)) = chi(-i, j).
                                for i in PrimeResidues(3^(ld-2)) do
                                    if i > 3^(ld-2) / 2 then
                                        break;
                                    fi;
                                    for j in [0,3] do
                                        rho := RepR(p, ld, si, r, t, [i,j], true);
                                        name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                        RecordIrrep(irrep_list, name, rho, l);
                                    od;
                                od;
                                for i in PrimeResidues(3^(ld-2)) do
                                    for j in [1,2] do
                                        rho := RepR(p, ld, si, r, t, [i,j], true);
                                        name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                        RecordIrrep(irrep_list, name, rho, l);
                                    od;
                                od;
                            else
                                # A is cyclic; A = <alpha> x <-1> where Ord(alpha) = p^(ld-si).
                                # A character is primitive iff it is injective on <alpha>.
                                for i in PrimeResidues(p^(ld-si)) do
                                    if i > p^(ld-si) / 2 then
                                        break;
                                    fi;
                                    for j in [0,1] do
                                        rho := RepR(p, ld, si, r, t, [i,j], true);
                                        name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                        RecordIrrep(irrep_list, name, rho, l);
                                    od;
                                od;
                            fi;
                        od;
                    od;
                od;

                # (R_ld(1)+-)_1 and (R_ld(u)+-)_1, for u a non-residue
                rho := RepRUnary(p, ld, 1, true);
                name := Concatenation("(R_", String(ld), "(1)+)_1");
                RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("(R_", String(ld), "(1)-)_1");
                RecordIrrep(irrep_list, name, rho[2], l);

                u := SomeQuadNonRes(p);
                rho := RepRUnary(p, ld, u, true);
                name := Concatenation("(R_", String(ld), "(", String(u), ")+)_1");
                RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("(R_", String(ld), "(", String(u), ")-)_1");
                RecordIrrep(irrep_list, name, rho[2], l);
            fi;
        od;
    od;

    # p = 2, ld = 1.
    l := 2;
    if 1 = deg then
        # Xi_6 = C_2, with s = t = -1.
        rho := [
            [[-1]],
            [[-1]],
            1
        ];
        name := "Xi_6";
        RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 2 = deg then
        # N_1(nu), the Steinberg representation.
        rho := RepN(2, 1, [0, 0], true);
        name := "N_1(nu)";
        RecordIrrep(irrep_list, name, rho, l);
    fi;

    # p = 2, ld = 2.
    l := 4;
    if 3 = deg then
        # D_2(chi)+-, for chi != 1. Note that there is only one non-trivial character,
        # indexed by [1,0]. RepD returns a list of the two subreps, + and - .
        rho := RepD(2, 2, [1, 0], true);
        name := "D_2([1,0])+";
        RecordIrrep(irrep_list, name, rho[1], l);
        name := "D_2([1,0])-";
        RecordIrrep(irrep_list, name, rho[2], l);

        # R_2^0(1,3)_1 and C_2 tensor R_2^0(1,3)_1 (exceptional).
        # Matrix given in notes sec. 1.6.2.
        w := Sqrt(2);
        rho := [
            1/2 * [
                [-1,  1,  w],
                [ 1, -1,  w],
                [ w,  w,  0]
            ],
            DiagonalMat([E(4), -E(4), 1]),
            3
        ];
        name := "R_2^0(1,3)_1";
        RecordIrrep(irrep_list, name, rho, l);
        rho := [-1*rho[1], -1*rho[2], 3];
        name := "Xi_6 tensor R_2^0(1,3)_1";
        RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 2 = deg then
        # N_2(chi), for chi primitive, chi^2 != 1.
        # A = <zeta> with ord(zeta) = 6.
        # For this specific case, chi is primitive iff chi(-1) = -1, which leaves only
        # a single relevant character, indexed by [0,1].
        rho := RepN(2, 2, [0, 1], true);
        name := "N_2([0,1])";
        RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 1 = deg then
        # Xi_3 = C_4, with s = -1, t = i.
        rho := [
            [[-E(4)]],
            [[E(4)]],
            1
        ];
        name := "Xi_3";
        RecordIrrep(irrep_list, name, rho, l);

        # Xi_9 = C_3, with s = i, t = -i.
        rho := [
            [[E(4)]],
            [[-E(4)]],
            1
        ];
        name := "Xi_9";
        RecordIrrep(irrep_list, name, rho, l);
    fi;

    # p = 2, ld = 3.
    l := 8;
    if 6 = deg then
        # D_3(chi)+-, for chi primitive.
        # A = <-1> x <5>, and a character is primitive if injective on <5>.
        # This gives two characters, [0,1] and [1,1] (note that both square to 1).
        # RepD returns a list of the two subreps, + and - .
        rho := RepD(2, 3, [0, 1], true);
        name := "D_3([0,1])+";
        RecordIrrep(irrep_list, name, rho[1], l);
        name := "D_3([0,1])-";
        RecordIrrep(irrep_list, name, rho[2], l);

        rho := RepD(2, 3, [1, 1], true);
        name := "D_3([1,1])+";
        RecordIrrep(irrep_list, name, rho[1], l);
        name := "D_3([1,1])-";
        RecordIrrep(irrep_list, name, rho[2], l);

        # R_3^0(1,3,nu)_1 and C_3 tensor R_3^0(1,3,nu)_1.
        # Matrix given in notes S. 1.6.3.
        rho := [
            -1/2 * [
                [ 0,  0,  1,  1,  1,  1],
                [ 0,  0, -1, -1,  1,  1],
                [ 1, -1,  0,  0,  1, -1],
                [ 1, -1,  0,  0, -1,  1],
                [ 1,  1,  1, -1,  0,  0],
                [ 1,  1, -1,  1,  0,  0],
            ],

            DiagonalMat([1, -1, E(8), E(8)^5, E(8)^3, E(8)^7]),

            6
        ];
        name := "R_3^0(1,3,nu)_1";
        RecordIrrep(irrep_list, name, rho, l);

        rho := [E(4) * rho[1], -E(4) * rho[2], 6];
        name := "Xi_9 tensor R_3^0(1,3,nu)_1";
        RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 4 = deg then
        # N_3(chi), for chi primitive, chi^2 != 1.
        # A = <alpha> x <zeta> where ord(alpha) = 2 and ord(zeta) = 6.
        # chi is primitive if injective on <alpha>.
        # This gives two relevant characters, [1,1] and [1,2].
        rho := RepN(2, 3, [1, 1], true);
        name := "N_3([1,1])";
        RecordIrrep(irrep_list, name, rho, l);

        rho := RepN(2, 3, [1, 2], true);
        name := "N_3([1,2])";
        RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 2 = deg then
        # N_3(chi)+-, for chi primitive, chi^2 = 1.
        # As above; relevant characters are those which square to 1: [1,0] and [1,3].
        # RepN returns a list of the two subreps, + and - .
        rho := RepN(2, 3, [1, 0], true);
        name := "N_3([1,0])+";
        RecordIrrep(irrep_list, name, rho[1], l);
        name := "N_3([1,0])-";
        RecordIrrep(irrep_list, name, rho[2], l);

        rho := RepN(2, 3, [1, 3], true);
        name := "N_3([1,3])+";
        RecordIrrep(irrep_list, name, rho[1], l);
        name := "N_3([1,3])-";
        RecordIrrep(irrep_list, name, rho[2], l);
    fi;
    if 3 = deg then
        # R_3^0(r,t,chi), for chi(zeta) = i and r in {1,3} and t in {1,5}.
        # A = <zeta> with ord(zeta) = 4.
        # Specifically, zeta = (0,1) when t = 1, and zeta = (2,1) when t = 5.
        # The given character is indexed by [0,1].
        for r in [1,3] do
            for t in [1,5] do
                rho := RepR(2, 3, 0, r, t, [0,1], true);
                name := Concatenation("R_3^0(", String(r), ",", String(t), ",[0,1])");
                RecordIrrep(irrep_list, name, rho, l);
            od;
        od;

        # R_3^0(1, t, chi)_+-, for chi != 1 and t in {3,7}.
        # |A| = 2, so there's only one non-trivial chi, indexed by [0,1].
        # chi squares to 1, giving two subreps.
        # RepR returns a list with the two supreps, + and -.
        for t in [3,7] do
            rho := RepR(2,3,0,1,t,[0,1],true);
            name := Concatenation("R_3^0(1,", String(t), ",[0,1])+");
            RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_3^0(1,", String(t), ",[0,1])-");
            RecordIrrep(irrep_list, name, rho[2], l);
        od;
    fi;

    # p = 2, ld = 4.
    l := 16;
    if 24 = deg then
        # D_4(chi), for chi primitive, chi^2 != 1.
        # A = <-1> x <5>, and a character is primitive if injective on <5>.
        # Note that 5 has order 4.
        # This gives two relevant characters, [0,1] and [1,1].
        rho := RepD(2, 4, [0, 1], true);
        name := "D_4([0,1])";
        RecordIrrep(irrep_list, name, rho, l);

        rho := RepD(2, 4, [1, 1], true);
        name := "D_4([1,1])";
        RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 8 = deg then
        # N_4(chi), for chi primitive, chi^2 != 1.
        # A = <alpha> x <zeta> where ord(alpha) = 4 and ord(zeta) = 6.
        # chi is primitive if injective on <alpha>.
        # This gives six relevant characters: [1,0],[1,3],[1,1],[1,2],[3,1],[3,2].
        for x in [[1,0],[1,3],[1,1],[1,2],[3,1],[3,2]] do
            rho := RepN(2, 4, x, true);
            name := Concatenation("N_4([", String(x[1]), ",", String(x[2]), "])");
            RecordIrrep(irrep_list, name, rho, l);
        od;
    fi;
    if 6 = deg then
        # R_4^0(r,t,chi), for chi primitive with chi^2 != 1, r in {1,3}, t in {1,5}.
        # For t = 1, A = <alpha> x <zeta>, with ord(alpha) = 2 and ord(zeta) = 4,
        # and a character is primitive iff injective on alpha, so the only relevant
        # character is [1,1].
        # For t = 5, A = <alpha> x <-1>, with ord(alpha) = 4; in this case a character
        # is primitive iff injective on <-alpha^2> (see NW p. 496), which means [1,0].
        for r in [1,3] do
            rho := RepR(2, 4, 0, r, 1, [1, 1], true);
            name := Concatenation("R_4^0(", String(r), ",1,[1,1])");
            RecordIrrep(irrep_list, name, rho, l);
        od;
        for r in [1,3] do
            rho := RepR(2, 4, 0, r, 5, [1, 0], true);
            name := Concatenation("R_4^0(", String(r), ",5,[1,0])");
            RecordIrrep(irrep_list, name, rho, l);
        od;

        # R_4^0(1,t,chi)+-, for chi primitive, t in {3,7}.
        # A = <alpha> x <-1>, with ord(alpha) = 2. There are therefore two primitive chars,
        # indexed by [1,0] and [1,1]. Both square to 1, giving two subreps.
        # RepR returns a list of the two subreps, + and - .
        for t in [3,7] do
            rho := RepR(2, 4, 0, 1, t, [1,0], true);
            name := Concatenation("R_4^0(1,", String(t), ",[1,0])+");
            RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(1,", String(t), ",[1,0])-");
            RecordIrrep(irrep_list, name, rho[2], l);

            rho := RepR(2, 4, 0, 1, t, [1,1], true);
            name := Concatenation("R_4^0(1,", String(t), ",[1,1])+");
            RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(1,", String(t), ",[1,1])-");
            RecordIrrep(irrep_list, name, rho[2], l);
        od;

        # R_4^2(r,t,chi) and C_2 tensor R_4^2(r,t,chi), for chi != 1, r,t in {1,3}.
        # A = {+-1}, so there is a single non-trivial character, indexed by [0,1].
        # NOTE: C_2 tensor R_4^2(1,1,[0,1]) is iso to R_4^2(1,3,nu)_1, and
        # C_2 tensor R_4^2(3,1,[0,1]) is iso to R_4^2(3,3,nu)_1.
        # NW lists them under the latter name, so this is what we do here; see below.
        for r in [1,3] do
            rho := RepR(2,4,2,r,1,[0,1],true);
            name := Concatenation("R_4^2(", String(r), ",1,[0,1])");
            RecordIrrep(irrep_list, name, rho, l);

            rho := RepR(2,4,2,r,3,[0,1],true);
            name := Concatenation("R_4^2(", String(r), ",3,[0,1])");
            RecordIrrep(irrep_list, name, rho, l);

            rho := [-1 * rho[1], -1 * rho[2], 6];
            name := Concatenation("Xi_6 tensor R_4^2(", String(r), ",3,[0,1])");
            RecordIrrep(irrep_list, name, rho, l);
        od;

        # R_4^2(r,3,nu)_1, for r in {1,3}.
        # Basis is given on NW p. 524.
        # See above.
        w := Sqrt(2);
        rho := [
            (1 / (2*w)) * [
                [ 1, -1,  1, -1,  w,  w],
                [-1,  1, -1,  1,  w,  w],
                [ 1, -1, -1,  1,  w, -w],
                [-1,  1,  1, -1,  w, -w],
                [ w,  w,  w,  w,  0,  0],
                [ w,  w, -w, -w,  0,  0],
            ],

            DiagonalMat([
                E(16),
                E(16)^9,
                E(16)^13,
                E(16)^5,
                1,
                E(16)^12
            ]),

            6
        ];
        name := "R_4^2(1,3,nu)_1";
        RecordIrrep(irrep_list, name, rho, l);

        rho := [
            (1 / (2*w)) * [
                [-1,  1, -1,  1,  w,  w],
                [ 1, -1,  1, -1,  w,  w],
                [-1,  1,  1, -1,  w, -w],
                [ 1, -1, -1,  1,  w, -w],
                [ w,  w,  w,  w,  0,  0],
                [ w,  w, -w, -w,  0,  0],
            ],

            DiagonalMat([
                E(16)^3,
                E(16)^11,
                E(16)^7,
                E(16)^15,
                1,
                E(16)^4
            ]),

            6
        ];
        name := "R_4^2(3,3,nu)_1";
        RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 3 = deg then
        # R_4^0(r,t,chi)+-, for chi primitive, chi^2 = 1, r in {1,3}, t in {1,5}.
        # See previous. For t = 1, the relevant chars. are [1,0] and [1,2].
        # For t = 5, the relevant chars. are [0,1] and [2,1].
        # RepR returns a list with the two supreps, + and -.
        for r in [1,3] do
            rho := RepR(2, 4, 0, r, 1, [1,0], true);
            name := Concatenation("R_4^0(", String(r), ",1,[1,0])+");
            RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(", String(r), ",1,[1,0])-");
            RecordIrrep(irrep_list, name, rho[2], l);

            rho := RepR(2, 4, 0, r, 1, [1,2], true);
            name := Concatenation("R_4^0(", String(r), ",1,[1,2])+");
            RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(", String(r), ",1,[1,2])-");
            RecordIrrep(irrep_list, name, rho[2], l);
        od;
        for r in [1,3] do
            rho := RepR(2, 4, 0, r, 5, [0,1], true);
            name := Concatenation("R_4^0(", String(r), ",5,[0,1])+");
            RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(", String(r), ",5,[0,1])-");
            RecordIrrep(irrep_list, name, rho[2], l);

            rho := RepR(2, 4, 0, r, 5, [2,1], true);
            name := Concatenation("R_4^0(", String(r), ",5,[2,1])+");
            RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(", String(r), ",5,[2,1])-");
            RecordIrrep(irrep_list, name, rho[2], l);
        od;
    fi;
    if 12 = deg then
        # N_3(chi)+ tensor R_4^0(1,7,psi)+, for psi(alpha) = -1, psi(-1) = 1,
        # chi primitive and chi^2 = 1.
        # Calculated via Kronecker product.
        # For the N part, chars. are [1,0] and [1,3] (see the ld=3 section earlier).
        #
        # First construct R_4^0(1,7,psi)+. Char is [1,1]. RepR gives a list with + and then -.
        x := RepR(2, 4, 0, 1, 7, [1,1], true);
        for j in [0,3] do
            # Construct N_3(chi)+. RepR again returns a list.
            rho := RepN(2, 3, [1, j], true);
            rho := [
                KroneckerProduct(rho[1][1], x[1][1]),
                KroneckerProduct(rho[1][2], x[1][2]),
                12
            ];
            name := Concatenation("N_3([1,", String(j), "])+ tensor R_4^0(1,7,[1,1])+");
            RecordIrrep(irrep_list, name, rho, l);
        od;
    fi;

    # p = 2, ld = 5.
    l := 32;
    if 48 = deg then
        # D_5(chi), for chi primitive.
        # A = <-1> x <5>, and a character is primitive if injective on <5>.
        # Note that 5 has order 8.
        # This gives four relevant characters, [0,1], [1,1], [0,3], [1,3].
        for i in [0,1] do
            for j in [1,3] do
                rho := RepD(2, 5, [i, j], true);
                name := Concatenation("D_5([", String(i), ",", String(j), "])");
                RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
    fi;
    if 16 = deg then
        # N_5(chi), for chi primitive.
        # A = <alpha> x <zeta> where ord(alpha) = 8 and ord(zeta) = 6.
        # chi is primitive if injective on <alpha>.
        for i in [1,3] do
            for j in [0,3] do
                rho := RepN(2, 5, [i,j], true);
                name := Concatenation("N_5([", String(i), ",", String(j), "])");
                RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
        for i in [1,3,5,7] do
            for j in [1,2] do
                rho := RepN(2, 5, [i,j], true);
                name := Concatenation("N_5([", String(i), ",", String(j), "])");
                RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
    fi;
    if 12 = deg then
        # R_5^0(r,t,chi), for chi primitive, r in {1,3}, t in {1,5}.
        # For t = 1, A = <alpha> x <zeta>, with ord(alpha) = 4 and ord(zeta) = 4.
        # For t = 5, A = <alpha> x <-1>, with ord(alpha) = 8.
        # In both cases, a character is primitive iff injective on alpha.
        for r in [1,3] do
            for j in [0,2] do
                rho := RepR(2, 5, 0, r, 1, [1,j], true);
                name := Concatenation("R_5^0(", String(r), ",1,[1,", String(j), "])");
                RecordIrrep(irrep_list, name, rho, l);
            od;
            for i in [1,3] do
                rho := RepR(2, 5, 0, r, 1, [i,1], true);
                name := Concatenation("R_5^0(", String(r), ",1,[", String(i), ",1])");
                RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
        for r in [1,3] do
            for i in [1,3] do
                for j in [0,1] do
                    rho := RepR(2, 5, 0, r, 5, [i,j], true);
                    name := Concatenation("R_5^0(", String(r), ",5,[", String(i), ",", String(j), "])");
                    RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
        od;

        # R_5^1(r,t,chi), for chi primitive and:
        # r in {1,5}, t in {1,5}, or
        # r in {1,3}, t in {3,7}.
        # A = <alpha> x <-1>; ord(alpha) = 4. Prim. iff injective on alpha.
        # Chars are therefore [0,1] and [1,1].
        for r in [1,5] do
            for t in [1,5] do
                for j in [0,1] do
                    rho := RepR(2, 5, 1, r, t, [1,j], true);
                    name := Concatenation("R_5^1(", String(r), ",", String(t), ",[1,", String(j), "])");
                    RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
        od;
        for r in [1,3] do
            for t in [3,7] do
                for j in [0,1] do
                    rho := RepR(2, 5, 1, r, t, [1,j], true);
                    name := Concatenation("R_5^1(", String(r), ",", String(t), ",[1,", String(j), "])");
                    RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
        od;

        # R_5^2(r,1,chi)_1 and C_3 tensor R_5^2(r,1,chi)_1,
        # for chi NOT primitive and r in {1,3}.
        # A = <alpha> x <-1>, so non-primitive characters are indexed by [0,0] and [0,1].
        for r in [1,3] do
            for j in [0,1] do
                rho := RepR(2, 5, 2, r, 1, [0,j], true);
                name := Concatenation("R_5^2(", String(r), ",1,[0,", String(j), "])_1");
                RecordIrrep(irrep_list, name, rho, l);

                rho := [E(4) * rho[1], -E(4) * rho[2], 12];
                name := Concatenation("Xi_9 tensor R_5^2(", String(r), ",1,[0,", String(j), "])_1");
                RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
    fi;
    if 24 = deg then
        # R_5^0(1,t,chi), for chi primitive and t in {3,7}.
        # A = <alpha> x <-1> with ord(alpha) = 4.
        # A character is primitive iff injective on <alpha>.
        # Relevant chars. are therefore [1,0] and [1,1].
        for t in [3,7] do
            for j in [0,1] do
                rho := RepR(2, 5, 0, 1, t, [1,j], true);
                name := Concatenation("R_5^0(1,", String(t), ",[1,", String(j), "])");
                RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
    fi;
    if 6 = deg then
        # R_5^2(r,t,chi)+-, for chi primitive, r in {1,3}, and t in {1,3,5,7}.
        # A = <alpha> x <-1>, with ord(alpha) = 2.
        # Relevant chars are therefore [1,0] and [1,1], both of which square to 1.
        # RepR returns a list of the two resulting subreps, + and - .
        for r in [1,3] do
            for t in [1,3,5,7] do
                for j in [0,1] do
                    rho := RepR(2, 5, 2, r, t, [1,j], true);
                    name := Concatenation("R_5^2(", String(r), ",", String(t), ",[1,", String(j), "])+");
                    RecordIrrep(irrep_list, name, rho[1], l);
                    name := Concatenation("R_5^2(", String(r), ",", String(t), ",[1,", String(j), "])-");
                    RecordIrrep(irrep_list, name, rho[2], l);
                od;
            od;
        od;
    fi;

    # p = 2, ld > 5.
    ldmax := Maximum([LogInt(deg, 2) + 4, 6]);
    ldset := [6..ldmax];
    for ld in ldset do
        l := 2^ld;
        if 3*2^(ld-1) = deg then
            # D_ld(chi), for chi primitive.
            # A = <-1> x <5>, and a character is primitive if injective on <5>.
            # Note that 5 has order 2^(ld-2).
            for i in [0,1] do
                for j in PrimeResidues(2^(ld-2) / 2) do
                    rho := RepD(2, ld, [i, j], true);
                    name := Concatenation("D_", String(ld), "([", String(i), ",", String(j), "])");
                    RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
        fi;
        if 2^(ld-1) = deg then
            # N_ld(chi), for chi primitive.
            # A = <alpha> x <zeta> where ord(alpha) = 2^(ld-2) and ord(zeta) = 6.
            # chi is primitive if injective on <alpha>.
            for i in PrimeResidues(2^(ld-2) / 2) do
                for j in [0,3] do
                    rho := RepN(2, ld, [i,j], true);
                    name := Concatenation("N_", String(ld), "([", String(i), ",", String(j), "])");
                    RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
            for i in PrimeResidues(2^(ld-2)) do
                for j in [1,2] do
                    rho := RepN(2, ld, [i,j], true);
                    name := Concatenation("N_", String(ld), "([", String(i), ",", String(j), "])");
                    RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
        fi;
        if 3*2^(ld-2) = deg then
            # R_ld^0(1,t,chi), for chi primitive and t in {3,7}.
            # A = <alpha> x <-1> with ord(alpha) = 2^(ld-3).
            # A character is primitive iff injective on <alpha>.
            for t in [3,7] do
                for i in PrimeResidues(2^(ld-3) / 2) do
                    for j in [0,1] do
                        rho := RepR(2, ld, 0, 1, t, [i,j], true);
                        name := Concatenation("R_", String(ld), "^0(1,", String(t), ",[", String(i), ",", String(j), "])");
                        RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            od;
        fi;
        if 3*2^(ld-3) = deg then
            # R_ld^si(r,t,chi), for chi primitive and:
            # si = 0, r in {1,3}, t = 1 (ord(alpha) = 2^(ld-3), ord(zeta) = 4),
            # si = 0, r in {1,3}, t = 5 (ord(alpha) = 2^(ld-2), ord(zeta) = 2),
            # si = 1, r in {1,5}, t in {1,5} or r in {1,3}, t in {3,7} (ord(alpha) = 2^(ld-3), ord(zeta) = 2),
            # si = 2, r in {1,3}, t in {1,3,5,7} (ord(alpha) = 2^(ld-4), ord(zeta) = 2).
            # In all cases, a character is primitive iff injective on <alpha>.
            for r in [1,3] do
                for i in PrimeResidues(2^(ld-3) / 2) do
                    for j in [0,2] do
                        rho := RepR(2, ld, 0, r, 1, [i,j], true);
                        name := Concatenation("R_", String(ld), "^0(", String(r), ",1,[", String(i), ",", String(j), "])");
                        RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
                for i in PrimeResidues(2^(ld-3)) do
                    rho := RepR(2, ld, 0, r, 1, [i,1], true);
                    name := Concatenation("R_", String(ld), "^0(", String(r), ",1,[", String(i), ",1])");
                    RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
            for r in [1,3] do
                for i in PrimeResidues(2^(ld-2) / 2) do
                    for j in [0,1] do
                        rho := RepR(2, ld, 0, r, 5, [i,j], true);
                        name := Concatenation("R_", String(ld), "^0(", String(r), ",5,[", String(i), ",", String(j), "])");
                        RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            od;
            for r in [1,5] do
                for t in [1,5] do
                    for i in PrimeResidues(2^(ld-3) / 2) do
                        for j in [0,1] do
                            rho := RepR(2, ld, 1, r, t, [i,j], true);
                            name := Concatenation("R_", String(ld), "^1(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                            RecordIrrep(irrep_list, name, rho, l);
                        od;
                    od;
                od;
            od;
            for r in [1,3] do
                for t in [3,7] do
                    for i in PrimeResidues(2^(ld-3) / 2) do
                        for j in [0,1] do
                            rho := RepR(2, ld, 1, r, t, [i,j], true);
                            name := Concatenation("R_", String(ld), "^1(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                            RecordIrrep(irrep_list, name, rho, l);
                        od;
                    od;
                od;
            od;
            for r in [1,3] do
                for t in [1,3,5,7] do
                    for i in PrimeResidues(2^(ld-4) / 2) do
                        for j in [0,1] do
                            rho := RepR(2, ld, 2, r, t, [i,j], true);
                            name := Concatenation("R_", String(ld), "^2(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                            RecordIrrep(irrep_list, name, rho, l);
                        od;
                    od;
                od;
            od;
        fi;
        if 3*2^(ld-4) = deg then
            # R_ld^si(r,t,chi), for chi primitive, si in {3, ..., ld-3}, and r,t in {1,3,5,7}.
            # ord(alpha) = 2^(ld-si-1), ord(zeta) = 2.
            # A character is primitive iff injective on <alpha>.
            for si in [3 .. (ld-3)] do
                for r in [1,3,5,7] do
                    for t in [1,3,5,7] do
                        for i in PrimeResidues(2^(ld-si-1) / 2) do
                            for j in [0,1] do
                                rho := RepR(2, ld, si, r, t, [i,j], true);
                                name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                RecordIrrep(irrep_list, name, rho, l);
                            od;
                        od;
                    od;
                od;
            od;

            # R_ld^(ld-2)(r,t,chi), for chi primitive, r in {1,3,5,7}, t in {1,3}.
            # A = <alpha> x <-1>, with ord(alpha) = 2.
            # A character is primitive iff injective on alpha (but note alternate def.
            # of primitive; see NW p. 496). Two such exist, indexed by [1,0] and [1,1].
            for r in [1,3,5,7] do
                for t in [1,3] do
                    rho := RepR(2, ld, ld-2, r, t, [1,0], true);
                    name := Concatenation("R_", String(ld), "^", String(ld-2), "(", String(r), ",", String(t), ",[1,0])");
                    RecordIrrep(irrep_list, name, rho, l);

                    rho := RepR(2, ld, ld-2, r, t, [1,1], true);
                    name := Concatenation("R_", String(ld), "^", String(ld-2), "(", String(r), ",", String(t), ",[1,1])");
                    RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            if ld = 6 then
                # R_6^4(r,t,nu)_1 and C_2 tensor R_6^4(r,t,nu)_1, for r in {1,3,5,7}, t in {1,3}.
                for r in [1,3,5,7] do
                    for t in [1,3] do
                        rho := RepR(2,6,4,r,t,[0,0],true);
                        name := Concatenation("R_6^4(", String(r), ",", String(t), ",nu)_1");
                        RecordIrrep(irrep_list, name, rho, l);

                        rho[1] := -1 * rho[1];
                        rho[2] := -1 * rho[2];
                        name := Concatenation("Xi_6 tensor R_6^4(", String(r), ",", String(t), ",nu)_1");
                        RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            else
                # R_ld^(ld-3)(r,t,chi)_1, for chi in {[0,0],[2,0]} (which both square to 1),
                # r in {1,3,5,7}, t in {1,3}.
                for r in [1,3,5,7] do
                    for t in [1,3] do
                        rho := RepR(2,ld,ld-3,r,t,[0,0],true);
                        name := Concatenation("R_", String(ld), "^", String(ld-3), "(", String(r), ",", String(t), ",nu)_1");
                        RecordIrrep(irrep_list, name, rho, l);

                        rho := RepR(2,ld,ld-3,r,t,[2,0],true);
                        name := Concatenation("R_", String(ld), "^", String(ld-3), "(", String(r), ",", String(t), ",[2,0])_1");
                        RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            fi;
        fi;
    od;

    Print("\n", Length(irrep_list), " irreps of degree ", deg, " found.\n");
    return irrep_list;
end;

PrimePowerIrrepsOfDegreeAtMost := function(max_deg)
    local output, deg, count;
    output := [];
    count := 0;

    for deg in [1..max_deg] do
        Print("Degree ", deg, ":\n");
        output[deg] := PrimePowerIrrepsOfDegree(deg);
        count := count + Length(output[deg]);
        Print("\n");
    od;

    Print(count, " total irreps found.\n");
    return output;
end;

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

IrrepsOfDegree := function(deg)
    local linears, prime_power_reps, i, ConstructIrreps, triv, factorizations, output, f;

    prime_power_reps := [];

    # collect prime-power-level irreps of deg dividing the given degree
    Print("Constructing irreps of prime-power level.\n");
    # The linear reps are denoted Xi_n, n in Z/12Z, with T = [zeta_12^n] and S = [i^n].
    # We handle these separately just for brevity in the output.
    prime_power_reps[1] := [
        [[[1]], [[1]], 1, 1, "Xi_0"], # Xi_0 = C_1

        [[[-1]], [[-1]], 1, 2, "Xi_6"], # Xi_6 = C_2

        [[[1]], [[E(3)]], 1, 3, "Xi_4"], # Xi_4 = R_1(1)- with p=3
        [[[1]], [[E(3)^2]], 1, 3, "Xi_8"], # Xi_8 = R_1(2)- with p=3

        [[[-E(4)]], [[E(4)]], 1, 4, "Xi_3"], # Xi_3 = C_4
        [[[E(4)]], [[-E(4)]], 1, 4, "Xi_9"], # Xi_9 = C_3

        [[[-1]], [[-1]], 1, 6, "Xi_2"], # = Xi_6 * Xi_8
        [[[-1]], [[-1]], 1, 6, "Xi_10"], # = Xi_6 * Xi_4

        [[[E(4)]], [[E(12)]], 1, 12, "Xi_1"], # = Xi_9 * Xi_4
        [[[E(4)]], [[E(12)^5]], 1, 12, "Xi_5"], # = Xi_9 * Xi_8
        [[[-E(4)]], [[E(12)^7]], 1, 12, "Xi_7"], # = Xi_3 * Xi_4
        [[[-E(4)]], [[E(12)^11]], 1, 12, "Xi_11"] # = Xi_3 * Xi_8
    ];
    for i in DivisorsInt(deg) do
        if i = 1 then
            continue;
        fi;
        Print("\nDegree ", i, ":\n");
        prime_power_reps[i] := PrimePowerIrrepsOfDegree(i);
    od;

    Print("\nConstructing tensor products.\n");

    ConstructIrreps := function(rho, factors, start)
        local i, eta, output, name, new_start;

        if Length(factors) = 0 then
            # Nothing left to do.
            return [rho];
        fi;

        output := [];

        # Tensor the given rep., rho, with all coprime irreps of degree (factors[1])
        # and recurse.
        for i in [start .. Length(prime_power_reps[factors[1]])] do
            eta := prime_power_reps[factors[1]][i];

            if Gcd(eta[4], rho[4]) = 1 then
                # No point listing Xi_0 = 1.
                if rho[5] = "Xi_0" then
                    name := Concatenation(eta[5], " [d:", String(eta[3]), ", l:", String(eta[4]), "]");
                elif eta[5] = "Xi_0" then
                    name := rho[5];
                else
                    name := Concatenation(rho[5], " (tensor) ", eta[5], " [d:", String(eta[3]), ", l:", String(eta[4]), "]");
                fi;

                if Length(factors) > 1 and factors[2] = factors[1] then
                    # Next degree is same as current. Need to keep irreps in order.
                    new_start := i+1;
                else
                    new_start := 1;
                fi;

                output := Concatenation(
                    output,
                    ConstructIrreps(
                        [
                            KroneckerProduct(rho[1], eta[1]),
                            KroneckerProduct(rho[2], eta[2]),
                            rho[3] * eta[3],
                            rho[4] * eta[4], # level is Lcm(level(rho), level(eta)), but they're coprime.
                            name
                        ],
                        factors{[2..Length(factors)]},
                        new_start
                    )
                );
            fi;
        od;

        return output;
    end;

    triv := [[[1]], [[1]], 1, 1, "Xi_0"];

    factorizations := Factorizations(deg);
    output := [];
    for f in factorizations do
        Append(output, ConstructIrreps(triv, f, 1));
    od;

    # Sort by level.
    SortBy(output, x -> x[4]);

    for i in [1 .. Length(output)] do
        Print(i, ": ( ", output[i][5], " ) [d: ", output[i][3], ", l: ", output[i][4] ,"]\n");
    od;
    Print("\nTotal count: ", Length(output), "\n");
    return output;
end;
