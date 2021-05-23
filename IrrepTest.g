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

PositionTest := function(irreps, rho, pos_list)
    local pos;

    pos := Position(irreps, ChiST(rho[1], rho[2]));
    if pos = fail then
        Print("FAILED!\n");
    else
        Print("pos ", pos, ", deg ", Length(rho[1]), "\n");
        Add(pos_list, pos);
    fi;
end;

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

IrrepTestOdd := function(p, ld)
    local l, G, irreps, pos_list, i, j, rho, u, si, r, t;

    l := p^ld;

    Print("Constructing irreps...\n");
    G := SL2Conj(p, ld);
    irreps := Irr(G);
    Print("# of irreps: ", Length(irreps), "\n\n");

    pos_list := [];

    if ld = 1 then
        # D_1(chi), for chi primitive, chi^2 != 1.
        # Defined for p >= 3.
        # A = <alpha>, ord(alpha) = p-1.
        # Relevant character indices: [(1..(p-3)/2), 0].
        if p > 3 then
            for i in [1 .. (p-3)/2] do
                rho := RepD(p, 1, [i, 0], true);
                Print("D_1([", i, ",0]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        fi;

        # N_1(chi), for chi primitive.
        # A = <zeta>, ord(zeta) = p+1.
        # Relevant character indices: [0, (1..(p-1)/2)].
        for j in [1 .. (p-1)/2] do
            rho := RepN(p, 1, [0, j], true);
            Print("N_1([0,", j, "]) : ");
            PositionTest(irreps, rho[1], pos_list);
        od;

        # R_1(1)+- and R_1(u)+-, for u a non-residue.
        # RepRUnary gives a list containing two reps, + and - .
        rho := RepRUnary(p, 1, 1, true);
        Print("R_1(1)+ : ");
        PositionTest(irreps, rho[1], pos_list);
        Print("R_1(1)- : ");
        PositionTest(irreps, rho[2], pos_list);

        u := SomeQuadNonRes(p);
        rho := RepRUnary(p, 1, u, true);
        Print("R_1(", u, ")+ : ");
        PositionTest(irreps, rho[1], pos_list);
        Print("R_1(", u, ")- : ");
        PositionTest(irreps, rho[2], pos_list);

        # N_1(nu), the Steinberg representation.
        rho := RepN(p, 1, [0, 0], true);
        Print("N_1(nu) : ");
        PositionTest(irreps, rho[1], pos_list);
    else
        # ld >= 2.

        # D_ld(chi), for chi primitive, chi^2 != 1.
        # A = <alpha>, with ord(alpha) = p^ld - p^(ld-1).
        # chi is primitive iff it is injective on <1+p = omicron>; this occurs
        # when the index of chi is coprime to p.
        for i in [1 .. (l - l/p)/2 - 1] do
            if Gcd(i, p) = 1 then
                rho := RepD(p, ld, [i, 0], true);
                Print("D_", ld, "([", i, ",0]) : ");
                PositionTest(irreps, rho[1], pos_list);
            fi;
        od;

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
            Print("N_", ld, "([", i, ",0]) : ");
            PositionTest(irreps, rho[1], pos_list);

            rho := RepN(p, ld, [i, (p+1)/2], true);
            Print("N_", ld, "([", i, ",", (p+1)/2, "]) : ");
            PositionTest(irreps, rho[1], pos_list);
        od;
        for i in PrimeResidues(p^(ld-1)) do
            for j in [1 .. (p-1)/2] do
                rho := RepN(p, ld, [i, j], true);
                Print("N_", ld, "([", i, ",", j, "]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;

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
                                Print("R_", ld, "^", si, "(", r, ",", t, ",[", i, ",", j, "]) : ");
                                PositionTest(irreps, rho[1], pos_list);
                            od;
                        od;
                        for i in PrimeResidues(3^(ld-2)) do
                            for j in [1,2] do
                                rho := RepR(p, ld, si, r, t, [i,j], true);
                                Print("R_", ld, "^", si, "(", r, ",", t, ",[", i, ",", j, "]) : ");
                                PositionTest(irreps, rho[1], pos_list);
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
                                Print("R_", ld, "^", si, "(", r, ",", t, ",[", i, ",", j, "]) : ");
                                PositionTest(irreps, rho[1], pos_list);
                            od;
                        od;
                    fi;
                od;
            od;
        od;

        # (R_ld(1)+-)_1 and (R_ld(u)+-)_1, for u a non-residue
        rho := RepRUnary(p, ld, 1, true);
        Print("(R_", ld, "(1)+)_1 : ");
        PositionTest(irreps, rho[1], pos_list);
        Print("(R_", ld, "(1)-)_1 : ");
        PositionTest(irreps, rho[2], pos_list);

        u := SomeQuadNonRes(p);
        rho := RepRUnary(p, ld, u, true);
        Print("(R_", ld, "(", u, ")+)_1 : ");
        PositionTest(irreps, rho[1], pos_list);
        Print("(R_", ld, "(", u, ")-)_1 : ");
        PositionTest(irreps, rho[2], pos_list);
    fi;

    if Length(pos_list) = Length(AsSet(pos_list)) then
        Print("\n", Length(pos_list), " full-level irreps found:\n", AsSet(pos_list), "\n");
    else
        Print("\nWARNING: duplicates found:\n", pos_list, "\n");
    fi;
end;

IrrepTestEven := function(ld)
    local l, G, irreps, pos_list, i, j, rho, si, w, r, t, x;

    l := 2^ld;

    Print("Constructing irreps...\n");
    G := SL2Conj(2, ld);
    irreps := Irr(G);
    Print("# of irreps: ", Length(irreps), "\n\n");

    pos_list := [];

    if ld = 1 then
        # C_2, with s = t = -1.
        rho := [
            [[-1]],
            [[-1]]
        ];
        Print("C_2 : ");
        PositionTest(irreps, rho, pos_list);

        # N_1(nu), the Steinberg representation.
        rho := RepN(2, 1, [0, 0], true);
        Print("N_1(nu) : ");
        PositionTest(irreps, rho[1], pos_list);
    elif ld = 2 then
        # D_2(chi)+-, for chi != 1. Note that there is only one non-trivial character,
        # indexed by [1,0]. RepD returns a list of the two subreps, + and - .
        rho := RepD(2, 2, [1, 0], true);
        Print("D_2([1,0])+ : ");
        PositionTest(irreps, rho[1], pos_list);
        Print("D_2([1,0])- : ");
        PositionTest(irreps, rho[2], pos_list);

        # R_2^0(1,3)_1 and C_2 tensor R_2^0(1,3)_1 (exceptional).
        # Matrix given in notes sec. 1.6.2.
        w := Sqrt(2);
        rho := [
            1/2 * [
                [-1,  1,  w],
                [ 1, -1,  w],
                [ w,  w,  0]
            ],
            DiagonalMat([E(4), -E(4), 1])
        ];
        Print("R_2^0(1,3)_1 : ");
        PositionTest(irreps, rho, pos_list);
        rho := [-1*rho[1], -1*rho[2]];
        Print("C_2 tensor R_2^0(1,3)_1 : ");
        PositionTest(irreps, rho, pos_list);

        # N_2(chi), for chi primitive, chi^2 != 1.
        # A = <zeta> with ord(zeta) = 6.
        # For this specific case, chi is primitive iff chi(-1) = -1, which leaves only
        # a single relevant character, indexed by [0,1].
        rho := RepN(2, 2, [0, 1], true);
        Print("N_2([0,1]) : ");
        PositionTest(irreps, rho[1], pos_list);

        # C_3, with s = i, t = -i.
        rho := [
            [[E(4)]],
            [[-E(4)]]
        ];
        Print("C_3 : ");
        PositionTest(irreps, rho, pos_list);

        # C_4, with s = -1, t = i.
        rho := [
            [[-E(4)]],
            [[E(4)]]
        ];
        Print("C_4 : ");
        PositionTest(irreps, rho, pos_list);
    elif ld = 3 then
        # D_3(chi)+-, for chi primitive.
        # A = <-1> x <5>, and a character is primitive if injective on <5>.
        # This gives two characters, [0,1] and [1,1] (note that both square to 1).
        # RepD returns a list of the two subreps, + and - .
        rho := RepD(2, 3, [0, 1], true);
        Print("D_3([0,1])+ : ");
        PositionTest(irreps, rho[1], pos_list);
        Print("D_3([0,1])- : ");
        PositionTest(irreps, rho[2], pos_list);

        rho := RepD(2, 3, [1, 1], true);
        Print("D_3([1,1])+ : ");
        PositionTest(irreps, rho[1], pos_list);
        Print("D_3([1,1])- : ");
        PositionTest(irreps, rho[2], pos_list);

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

            DiagonalMat([1, -1, E(8), E(8)^5, E(8)^3, E(8)^7])
        ];
        Print("R_3^0(1,3,nu)_1 : ");
        PositionTest(irreps, rho, pos_list);

        rho := [E(4) * rho[1], -E(4) * rho[2]];
        Print("C_3 tensor R_3^0(1,3,nu)_1 : ");
        PositionTest(irreps, rho, pos_list);

        # N_3(chi), for chi primitive, chi^2 != 1.
        # A = <alpha> x <zeta> where ord(alpha) = 2 and ord(zeta) = 6.
        # chi is primitive if injective on <alpha>.
        # This gives two relevant characters, [1,1] and [1,2].
        rho := RepN(2, 3, [1, 1], true);
        Print("N_3([1,1]) : ");
        PositionTest(irreps, rho[1], pos_list);

        rho := RepN(2, 3, [1, 2], true);
        Print("N_3([1,2]) : ");
        PositionTest(irreps, rho[1], pos_list);

        # N_3(chi)+-, for chi primitive, chi^2 = 1.
        # As above; relevant characters are those which square to 1: [1,0] and [1,3].
        # RepN returns a list of the two subreps, + and - .
        rho := RepN(2, 3, [1, 0], true);
        Print("N_3([1,0])+ : ");
        PositionTest(irreps, rho[1], pos_list);
        Print("N_3([1,0])- : ");
        PositionTest(irreps, rho[2], pos_list);

        rho := RepN(2, 3, [1, 3], true);
        Print("N_3([1,3])+ : ");
        PositionTest(irreps, rho[1], pos_list);
        Print("N_3([1,3])- : ");
        PositionTest(irreps, rho[2], pos_list);

        # R_3^0(r,t,chi), for chi(zeta) = i and r in {1,3} and t in {1,5}.
        # A = <zeta> with ord(zeta) = 4.
        # Specifically, zeta = (0,1) when t = 1, and zeta = (2,1) when t = 5.
        # The given character is indexed by [0,1].
        for r in [1,3] do
            for t in [1,5] do
                rho := RepR(2, 3, 0, r, t, [0,1], true);
                Print("R_3^0(", r, ",", t, ",[0,1]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;

        # R_3^0(1, t, chi)_+-, for chi != 1 and t in {3,7}.
        # |A| = 2, so there's only one non-trivial chi, indexed by [0,1].
        # chi squares to 1, giving two subreps.
        # RepR returns a list with the two supreps, + and -.
        for t in [3,7] do
            rho := RepR(2,3,0,1,t,[0,1],true);
            Print("R_3^0(1,", t, ",[0,1])+ : ");
            PositionTest(irreps, rho[1], pos_list);
            Print("R_3^0(1,", t, ",[0,1])- : ");
            PositionTest(irreps, rho[2], pos_list);
        od;
    elif ld = 4 then
        # D_4(chi), for chi primitive, chi^2 != 1.
        # A = <-1> x <5>, and a character is primitive if injective on <5>.
        # Note that 5 has order 4.
        # This gives two relevant characters, [0,1] and [1,1].
        rho := RepD(2, 4, [0, 1], true);
        Print("D_4([0,1]) : ");
        PositionTest(irreps, rho[1], pos_list);

        rho := RepD(2, 4, [1, 1], true);
        Print("D_4([1,1]) : ");
        PositionTest(irreps, rho[1], pos_list);

        # N_4(chi), for chi primitive, chi^2 != 1.
        # A = <alpha> x <zeta> where ord(alpha) = 4 and ord(zeta) = 6.
        # chi is primitive if injective on <alpha>.
        # This gives six relevant characters: [1,0],[1,3],[1,1],[1,2],[3,1],[3,2].
        for x in [[1,0],[1,3],[1,1],[1,2],[3,1],[3,2]] do
            rho := RepN(2, 4, x, true);
            Print("N_4([", x[1], ",", x[2], "]) : ");
            PositionTest(irreps, rho[1], pos_list);
        od;

        # R_4^0(r,t,chi), for chi primitive with chi^2 != 1, r in {1,3}, t in {1,5}.
        # For t = 1, A = <alpha> x <zeta>, with ord(alpha) = 2 and ord(zeta) = 4,
        # and a character is primitive iff injective on alpha, so the only relevant
        # character is [1,1].
        # For t = 5, A = <alpha> x <-1>, with ord(alpha) = 4; in this case a character
        # is primitive iff injective on <-alpha^2> (see NW p. 496), which means [1,0].
        for r in [1,3] do
            rho := RepR(2, 4, 0, r, 1, [1, 1], true);
            Print("R_4^0(", r, ",1,[1,1]) : ");
            PositionTest(irreps, rho[1], pos_list);
        od;
        for r in [1,3] do
            rho := RepR(2, 4, 0, r, 5, [1, 0], true);
            Print("R_4^0(", r, ",5,[1,0]) : ");
            PositionTest(irreps, rho[1], pos_list);
        od;

        # R_4^0(r,t,chi)+-, for chi primitive, chi^2 = 1, r in {1,3}, t in {1,5}.
        # See previous. For t = 1, the relevant chars. are [1,0] and [1,2].
        # For t = 5, the relevant chars. are [0,1] and [2,1].
        # RepR returns a list with the two supreps, + and -.
        for r in [1,3] do
            rho := RepR(2, 4, 0, r, 1, [1,0], true);
            Print("R_4^0(", r, ",1,[1,0])+ : ");
            PositionTest(irreps, rho[1], pos_list);
            Print("R_4^0(", r, ",1,[1,0])- : ");
            PositionTest(irreps, rho[2], pos_list);

            rho := RepR(2, 4, 0, r, 1, [1,2], true);
            Print("R_4^0(", r, ",1,[1,2])+ : ");
            PositionTest(irreps, rho[1], pos_list);
            Print("R_4^0(", r, ",1,[1,2])- : ");
            PositionTest(irreps, rho[2], pos_list);
        od;
        for r in [1,3] do
            rho := RepR(2, 4, 0, r, 5, [0,1], true);
            Print("R_4^0(", r, ",5,[0,1])+ : ");
            PositionTest(irreps, rho[1], pos_list);
            Print("R_4^0(", r, ",5,[0,1])- : ");
            PositionTest(irreps, rho[2], pos_list);

            rho := RepR(2, 4, 0, r, 5, [2,1], true);
            Print("R_4^0(", r, ",5,[2,1])+ : ");
            PositionTest(irreps, rho[1], pos_list);
            Print("R_4^0(", r, ",5,[2,1])- : ");
            PositionTest(irreps, rho[2], pos_list);
        od;

        # R_4^0(1,t,chi)+-, for chi primitive, t in {3,7}.
        # A = <alpha> x <-1>, with ord(alpha) = 2. There are therefore two primitive chars,
        # indexed by [1,0] and [1,1]. Both square to 1, giving two subreps.
        # RepR returns a list of the two subreps, + and - .
        for t in [3,7] do
            rho := RepR(2, 4, 0, 1, t, [1,0], true);
            Print("R_4^0(1,", t, ",[1,0])+ : ");
            PositionTest(irreps, rho[1], pos_list);
            Print("R_4^0(1,", t, ",[1,0])- : ");
            PositionTest(irreps, rho[2], pos_list);

            rho := RepR(2, 4, 0, 1, t, [1,1], true);
            Print("R_4^0(1,", t, ",[1,1])+ : ");
            PositionTest(irreps, rho[1], pos_list);
            Print("R_4^0(1,", t, ",[1,1])- : ");
            PositionTest(irreps, rho[2], pos_list);
        od;

        # R_4^2(r,t,chi) and C_2 tensor R_4^2(r,t,chi), for chi != 1, r,t in {1,3}.
        # A = {+-1}, so there is a single non-trivial character, indexed by [0,1].
        # NOTE: C_2 tensor R_4^2(1,1,[0,1]) is iso to R_4^2(1,3,nu)_1, and
        # C_2 tensor R_4^2(3,1,[0,1]) is iso to R_4^2(3,3,nu)_1.
        # NW lists them under the latter name, so this is what we do here; see below.
        for r in [1,3] do
            rho := RepR(2,4,2,r,1,[0,1],true);
            Print("R_4^2(", r, ",1,[0,1]) : ");
            PositionTest(irreps, rho[1], pos_list);

            rho := RepR(2,4,2,r,3,[0,1],true);
            Print("R_4^2(", r, ",3,[0,1]) : ");
            PositionTest(irreps, rho[1], pos_list);

            rho := [-1 * rho[1][1], -1 * rho[1][2]];
            Print("C_2 tensor R_4^2(", r, ",3,[0,1]) : ");
            PositionTest(irreps, rho, pos_list);
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
            ])
        ];
        Print("R_4^2(1,3,nu)_1 : ");
        PositionTest(irreps, rho, pos_list);

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
            ])
        ];
        Print("R_4^2(3,3,nu)_1 : ");
        PositionTest(irreps, rho, pos_list);

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
                KroneckerProduct(rho[1][2], x[1][2])
            ];
            Print("N_3([1,", j, "])+ tensor R_4^0(1,7,[1,1])+ : ");
            PositionTest(irreps, rho, pos_list);
        od;
    elif ld = 5 then
        # D_5(chi), for chi primitive.
        # A = <-1> x <5>, and a character is primitive if injective on <5>.
        # Note that 5 has order 8.
        # This gives four relevant characters, [0,1], [1,1], [0,3], [1,3].
        for i in [0,1] do
            for j in [1,3] do
                rho := RepD(2, 5, [i, j], true);
                Print("D_5([", i, ",", j, "]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;

        # N_5(chi), for chi primitive.
        # A = <alpha> x <zeta> where ord(alpha) = 8 and ord(zeta) = 6.
        # chi is primitive if injective on <alpha>.
        for i in [1,3] do
            for j in [0,3] do
                rho := RepN(2, 5, [i,j], true);
                Print("N_5([", i, ",", j, "]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;
        for i in [1,3,5,7] do
            for j in [1,2] do
                rho := RepN(2, 5, [i,j], true);
                Print("N_5([", i, ",", j, "]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;

        # R_5^0(r,t,chi), for chi primitive, r in {1,3}, t in {1,5}.
        # For t = 1, A = <alpha> x <zeta>, with ord(alpha) = 4 and ord(zeta) = 4.
        # For t = 5, A = <alpha> x <-1>, with ord(alpha) = 8.
        # In both cases, a character is primitive iff injective on alpha.
        for r in [1,3] do
            for j in [0,2] do
                rho := RepR(2, ld, 0, r, 1, [1,j], true);
                Print("R_", ld, "^0(", r, ",1,[1,", j, "]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
            for i in [1,3] do
                rho := RepR(2, ld, 0, r, 1, [i,1], true);
                Print("R_", ld, "^0(", r, ",1,[", i, ",1]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;
        for r in [1,3] do
            for i in [1,3] do
                for j in [0,1] do
                    rho := RepR(2, ld, 0, r, 5, [i,j], true);
                    Print("R_", ld, "^0(", r, ",5,[", i, ",", j, "]) : ");
                    PositionTest(irreps, rho[1], pos_list);
                od;
            od;
        od;

        # R_5^0(1,t,chi), for chi primitive and t in {3,7}.
        # A = <alpha> x <-1> with ord(alpha) = 4.
        # A character is primitive iff injective on <alpha>.
        # Relevant chars. are therefore [1,0] and [1,1].
        for t in [3,7] do
            for j in [0,1] do
                rho := RepR(2, ld, 0, 1, t, [1,j], true);
                Print("R_", ld, "^0(1,", t, ",[1,", j, "]) : ");
                PositionTest(irreps, rho[1], pos_list);
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
                    rho := RepR(2, ld, 1, r, t, [1,j], true);
                    Print("R_", ld, "^1(", r, ",", t, ",[1,", j, "]) : ");
                    PositionTest(irreps, rho[1], pos_list);
                od;
            od;
        od;
        for r in [1,3] do
            for t in [3,7] do
                for j in [0,1] do
                    rho := RepR(2, ld, 1, r, t, [1,j], true);
                    Print("R_", ld, "^1(", r, ",", t, ",[1,", j, "]) : ");
                    PositionTest(irreps, rho[1], pos_list);
                od;
            od;
        od;

        # R_5^2(r,t,chi)+-, for chi primitive, r in {1,3}, and t in {1,3,5,7}.
        # A = <alpha> x <-1>, with ord(alpha) = 2.
        # Relevant chars are therefore [1,0] and [1,1], both of which square to 1.
        # RepR returns a list of the two resulting subreps, + and - .
        for r in [1,3] do
            for t in [1,3,5,7] do
                for j in [0,1] do
                    rho := RepR(2, 5, 2, r, t, [1,j], true);
                    Print("R_5^2(", r, ",", t, ",[1,", j, "])+ : ");
                    PositionTest(irreps, rho[1], pos_list);
                    Print("R_5^2(", r, ",", t, ",[1,", j, "])- : ");
                    PositionTest(irreps, rho[2], pos_list);
                od;
            od;
        od;

        # R_5^2(r,1,chi)_1 and C_3 tensor R_5^2(r,1,chi)_1,
        # for chi NOT primitive and r in {1,3}.
        # A = <alpha> x <-1>, so non-primitive characters are indexed by [0,0] and [0,1].
        for r in [1,3] do
            for j in [0,1] do
                rho := RepR(2, 5, 2, r, 1, [0,j], true);
                Print("R_5^2(", r, ",1,[0,", j, "])_1 : ");
                PositionTest(irreps, rho[1], pos_list);

                rho := [E(4) * rho[1][1], -E(4) * rho[1][2]];
                Print("C_3 tensor R_5^2(", r, ",1,[0,", j, "])_1 : ");
                PositionTest(irreps, rho, pos_list);
            od;
        od;
    else
        # ld >= 6.

        # D_ld(chi), for chi primitive.
        # A = <-1> x <5>, and a character is primitive if injective on <5>.
        # Note that 5 has order 2^(ld-2).
        for i in [0,1] do
            for j in PrimeResidues(2^(ld-2) / 2) do
                rho := RepD(2, ld, [i, j], true);
                Print("D_", ld, "([", i, ",", j, "]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;

        # N_ld(chi), for chi primitive.
        # A = <alpha> x <zeta> where ord(alpha) = 2^(ld-2) and ord(zeta) = 6.
        # chi is primitive if injective on <alpha>.
        for i in PrimeResidues(2^(ld-2) / 2) do
            for j in [0,3] do
                rho := RepN(2, ld, [i,j], true);
                Print("N_", ld, "([", i, ",", j, "]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;
        for i in PrimeResidues(2^(ld-2)) do
            for j in [1,2] do
                rho := RepN(2, ld, [i,j], true);
                Print("N_", ld, "([", i, ",", j, "]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;

        # R_ld^0(1,t,chi), for chi primitive and t in {3,7}.
        # A = <alpha> x <-1> with ord(alpha) = 2^(ld-3).
        # A character is primitive iff injective on <alpha>.
        for t in [3,7] do
            for i in PrimeResidues(2^(ld-3) / 2) do
                for j in [0,1] do
                    rho := RepR(2, ld, 0, 1, t, [i,j], true);
                    Print("R_", ld, "^0(1,", t, ",[", i, ",", j, "]) : ");
                    PositionTest(irreps, rho[1], pos_list);
                od;
            od;
        od;

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
                    Print("R_", ld, "^0(", r, ",1,[", i, ",", j, "]) : ");
                    PositionTest(irreps, rho[1], pos_list);
                od;
            od;
            for i in PrimeResidues(2^(ld-3)) do
                rho := RepR(2, ld, 0, r, 1, [i,1], true);
                Print("R_", ld, "^0(", r, ",1,[", i, ",1]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;
        for r in [1,3] do
            for i in PrimeResidues(2^(ld-2) / 2) do
                for j in [0,1] do
                    rho := RepR(2, ld, 0, r, 5, [i,j], true);
                    Print("R_", ld, "^0(", r, ",5,[", i, ",", j, "]) : ");
                    PositionTest(irreps, rho[1], pos_list);
                od;
            od;
        od;
        for r in [1,5] do
            for t in [1,5] do
                for i in PrimeResidues(2^(ld-3) / 2) do
                    for j in [0,1] do
                        rho := RepR(2, ld, 1, r, t, [i,j], true);
                        Print("R_", ld, "^1(", r, ",", t, ",[", i, ",", j, "]) : ");
                        PositionTest(irreps, rho[1], pos_list);
                    od;
                od;
            od;
        od;
        for r in [1,3] do
            for t in [3,7] do
                for i in PrimeResidues(2^(ld-3) / 2) do
                    for j in [0,1] do
                        rho := RepR(2, ld, 1, r, t, [i,j], true);
                        Print("R_", ld, "^1(", r, ",", t, ",[", i, ",", j, "]) : ");
                        PositionTest(irreps, rho[1], pos_list);
                    od;
                od;
            od;
        od;
        for r in [1,3] do
            for t in [1,3,5,7] do
                for i in PrimeResidues(2^(ld-4) / 2) do
                    for j in [0,1] do
                        rho := RepR(2, ld, 2, r, t, [i,j], true);
                        Print("R_", ld, "^2(", r, ",", t, ",[", i, ",", j, "]) : ");
                        PositionTest(irreps, rho[1], pos_list);
                    od;
                od;
            od;
        od;

        # R_ld^si(r,t,chi), for chi primitive, si in {3, ..., ld-3}, and r,t in {1,3,5,7}.
        # ord(alpha) = 2^(ld-si-1), ord(zeta) = 2.
        # A character is primitive iff injective on <alpha>.
        for si in [3 .. (ld-3)] do
            for r in [1,3,5,7] do
                for t in [1,3,5,7] do
                    for i in PrimeResidues(2^(ld-si-1) / 2) do
                        for j in [0,1] do
                            rho := RepR(2, ld, si, r, t, [i,j], true);
                            Print("R_", ld, "^", si, "(", r, ",", t, ",[", i, ",", j, "]) : ");
                            PositionTest(irreps, rho[1], pos_list);
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
                Print("R_", ld, "^", ld-2, "(", r, ",", t, ",[1,0]) : ");
                PositionTest(irreps, rho[1], pos_list);

                rho := RepR(2, ld, ld-2, r, t, [1,1], true);
                Print("R_", ld, "^", ld-2, "(", r, ",", t, ",[1,1]) : ");
                PositionTest(irreps, rho[1], pos_list);
            od;
        od;

        if ld = 6 then
            # R_6^4(r,t,nu)_1 and C_2 tensor R_6^4(r,t,nu)_1, for r in {1,3,5,7}, t in {1,3}.
            for r in [1,3,5,7] do
                for t in [1,3] do
                    rho := RepR(2,6,4,r,t,[0,0],true);
                    Print("R_6^4(", r, ",", t, ",nu)_1 : ");
                    PositionTest(irreps, rho[1], pos_list);

                    rho := [-1 * rho[1][1], -1 * rho[1][2]];
                    Print("C_2 tensor R_6^4(", r, ",", t, ",nu)_1 : ");
                    PositionTest(irreps, rho, pos_list);
                od;
            od;
        else
            # R_ld^(ld-3)(r,t,chi)_1, for chi in {[0,0],[2,0]} (which both square to 1),
            # r in {1,3,5,7}, t in {1,3}.
            for r in [1,3,5,7] do
                for t in [1,3] do
                    rho := RepR(2,ld,ld-3,r,t,[0,0],true);
                    Print("R_", ld, "^", ld-3, "(", r, ",", t, ",nu)_1 : ");
                    PositionTest(irreps, rho[1], pos_list);

                    rho := RepR(2,ld,ld-3,r,t,[2,0],true);
                    Print("R_", ld, "^", ld-3, "(", r, ",", t, ",[2,0])_1 : ");
                    PositionTest(irreps, rho[1], pos_list);
                od;
            od;
        fi;
    fi;

    if Length(pos_list) = Length(AsSet(pos_list)) then
        Print("\n", Length(pos_list), " full-level irreps found:\n", AsSet(pos_list), "\n");
    else
        Print("\nWARNING: duplicates found:\n", pos_list, "\n");
    fi;
end;

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# IrrepTest := function(max_level)
#     local p, ld;

#     LogTo("./output/irrep_test.txt");

#     # Note: Primes is a built-in, global list of the 168 primes less than 1000.
#     # If we need higher primes we can calculate them.
#     for p in Primes do
#         if p > max_level then
#             break;
#         fi;

#         if p = 2 then
#             ld := 1;
#             while 2^ld < max_level do
#                 IrrepTestEven(ld);
#                 ld := ld + 1;
#             od;
#         else
#             ld := 1;
#             while p^ld < max_level do
#                 IrrepTestOdd(p, ld);
#                 ld := ld + 1;
#             od;
#         fi;
#     od;

#     LogTo();
# end;
