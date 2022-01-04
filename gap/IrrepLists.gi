#
# SL2Reps: Constructing symmetric representations of SL(2,Z).
#
# Lists of Representations.
#
# Implementations
#


InstallGlobalFunction( _SL2IrrepsPPLOfDegree,
function(degree)
    local irrep_list, p, ld, pmax, ldmax, pset, ldset, name, rho, l,
            i, j, u, w, si, r, t, x;

    if not degree in PositiveIntegers then
        Error("degree must be a positive integer.");
    fi;

    irrep_list := [];

    if 1 = degree then
        # The trivial irrep Xi_0 = C_1; level 1.
        name := "Xi_0";
        _SL2RecordIrrep(irrep_list, name, [[[1]], [[1]]], 1);
    fi;

    # p odd, ld = 1.
    pmax := 2*degree + 1;
    pset := Filtered([3..pmax], x -> IsPrime(x));
    for p in pset do
        l := p^1;
        if p + 1 = degree and p > 3 then
            # D_1(chi), for chi primitive, chi^2 != 1.
            # Defined for p >= 3.
            # A = <beta> (alpha is trivial) with ord(beta) = p-1.
            # All non-trivial characters are primitive.
            # Relevant character indices: [0, (1..(p-3)/2)].
            for i in [1 .. (p-3)/2] do
                rho := SL2IrrepD(p, 1, [0, i])[1];
                name := Concatenation("D_1([0,", String(i), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
        fi;

        if p - 1 = degree then
            # N_1(chi), for chi primitive.
            # A = <zeta>, ord(zeta) = p+1.
            # Relevant character indices: [0, (1..(p-1)/2)].
            for j in [1 .. (p-1)/2] do
                rho := SL2IrrepN(p, 1, [0, j])[1];
                name := Concatenation("N_1([0,", String(j), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
        fi;

        if (p+1)/2 = degree then
            # R_1(1)+ and R_1(u)+, for u a non-residue.
            # SL2IrrepRUnary gives a list containing two reps, + and - .
            rho := SL2IrrepRUnary(p, 1, 1)[1];
            name := "R_1(1)+";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            u := _SL2QuadNonRes(p);
            rho := SL2IrrepRUnary(p, 1, u)[1];
            name := Concatenation("R_1(", String(u), ")+");
            _SL2RecordIrrep(irrep_list, name, rho, l);
        fi;

        if (p-1)/2 = degree then
            # R_1(1)- and R_1(u)-, for u a non-residue.
            # SL2IrrepRUnary gives a list containing two reps, + and - .
            # For p=3, R_1(1)- = Xi_4 and R_1(2)- = Xi_8, the two linear reps. of level 3.
            rho := SL2IrrepRUnary(p, 1, 1)[2];
            name := "R_1(1)-";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            u := _SL2QuadNonRes(p);
            rho := SL2IrrepRUnary(p, 1, u)[2];
            name := Concatenation("R_1(", String(u), ")-");
            _SL2RecordIrrep(irrep_list, name, rho, l);
        fi;

        if p = degree then
            # N_1(nu), the Steinberg representation.
            rho := SL2IrrepN(p, 1, [0, 0])[1];
            name := "N_1(nu)";
            _SL2RecordIrrep(irrep_list, name, rho, l);
        fi;
    od;

    # p odd, ld >= 2.
    pmax := Maximum([RootInt(2*degree+1, 2), 3]);
    ldmax := Maximum([LogInt(degree, 3) + 2, 2]);
    pset := Filtered([3..pmax], x -> IsPrime(x));
    ldset := [2..ldmax];
    for p in pset do
        for ld in ldset do
            l := p^ld;
            if (p+1)*p^(ld-1) = degree then
                # D_ld(chi), for chi primitive, chi^2 != 1.
                # A = <alpha> x <beta>, with ord(alpha) = p^(ld-1) and ord(beta) = p-1.
                # chi is primitive iff it is injective on <alpha>.
                # Remember chi and bar(chi) give same rep; for j = 0 or (p-1)/2,
                # bar(chi(i,j)) = chi(-i, j).
                # We handle that case first:
                for i in PrimeResidues(p^(ld-1)) do
                    if i > p^(ld-1) / 2 then
                        break;
                    fi;
                    rho := SL2IrrepD(p, ld, [i, 0])[1];
                    name := Concatenation("D_", String(ld), "([", String(i), ",0])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);

                    rho := SL2IrrepD(p, ld, [i, (p-1)/2])[1];
                    name := Concatenation("D_", String(ld), "([", String(i), ",", String((p-1)/2), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
                # Then the rest:
                for i in PrimeResidues(p^(ld-1)) do
                    for j in [1 .. (p-3)/2] do
                        if Gcd(i, p) = 1 then
                            rho := SL2IrrepD(p, ld, [i, j])[1];
                            name := Concatenation("D_", String(ld), "([", String(i), ",", String(j), "])");
                            _SL2RecordIrrep(irrep_list, name, rho, l);
                        fi;
                    od;
                od;
            fi;

            if (p-1)*p^(ld-1) = degree then
                # N_ld(chi), for chi primitive.
                # A = <alpha> x <zeta>; ord(alpha) = p^(ld-1), ord(zeta) = p+1.
                # chi is primitive if injective on <alpha>.
                # Remember chi and bar(chi) give same rep; for j = 0 or (p+1)/2,
                # bar(chi(i,j)) = chi(-i, j).
                # We handle that case first:
                for i in PrimeResidues(p^(ld-1)) do
                    if i > p^(ld-1) / 2 then
                        break;
                    fi;
                    rho := SL2IrrepN(p, ld, [i, 0])[1];
                    name := Concatenation("N_", String(ld), "([", String(i), ",0])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);

                    rho := SL2IrrepN(p, ld, [i, (p+1)/2])[1];
                    name := Concatenation("N_", String(ld), "([", String(i), ",", String((p+1)/2), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
                # Then the rest:
                for i in PrimeResidues(p^(ld-1)) do
                    for j in [1 .. (p-1)/2] do
                        rho := SL2IrrepN(p, ld, [i, j])[1];
                        name := Concatenation("N_", String(ld), "([", String(i), ",", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            fi;

            if p^(ld-2)*(p^2-1)/2 = degree then
                # R_ld^si(r,t,chi), for 1 <= si <= ld-1, r,t either 1 or a non-residue, and chi primitive
                for si in [1 .. ld-1] do
                    u := _SL2QuadNonRes(p);
                    for r in [1,u] do
                        for t in [1,u] do
                            if p = 3 and ld >= 3 and si = 1 and t = 1 then
                                # Special case, see RepR.
                                # A = <alpha> x <zeta> where ord(alpha) = 3^(ld-2) and ord(zeta) = 6.
                                # chi is primitive if injective on <alpha>.
                                # Remember chi and bar(chi) give same rep; for j = 0 or 3,
                                # bar(chi(i,j)) = chi(-i, j).
                                for i in PrimeResidues(3^(ld-2)) do
                                    if i > 3^(ld-2) / 2 then
                                        break;
                                    fi;
                                    for j in [0,3] do
                                        rho := SL2IrrepR(p, ld, si, r, t, [i,j])[1];
                                        name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                        _SL2RecordIrrep(irrep_list, name, rho, l);
                                    od;
                                od;
                                for i in PrimeResidues(3^(ld-2)) do
                                    for j in [1,2] do
                                        rho := SL2IrrepR(p, ld, si, r, t, [i,j])[1];
                                        name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                        _SL2RecordIrrep(irrep_list, name, rho, l);
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
                                        rho := SL2IrrepR(p, ld, si, r, t, [i,j])[1];
                                        name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                        _SL2RecordIrrep(irrep_list, name, rho, l);
                                    od;
                                od;
                            fi;
                        od;
                    od;
                od;

                # (R_ld(1)+-)_1 and (R_ld(u)+-)_1, for u a non-residue
                rho := SL2IrrepRUnary(p, ld, 1);
                name := Concatenation("(R_", String(ld), "(1)+)_1");
                _SL2RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("(R_", String(ld), "(1)-)_1");
                _SL2RecordIrrep(irrep_list, name, rho[2], l);

                u := _SL2QuadNonRes(p);
                rho := SL2IrrepRUnary(p, ld, u);
                name := Concatenation("(R_", String(ld), "(", String(u), ")+)_1");
                _SL2RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("(R_", String(ld), "(", String(u), ")-)_1");
                _SL2RecordIrrep(irrep_list, name, rho[2], l);
            fi;
        od;
    od;

    # p = 2, ld = 1.
    l := 2;
    if 1 = degree then
        # Xi_6 = C_2, with s = t = -1.
        rho := [
            [[-1]],
            [[-1]]
        ];
        name := "Xi_6";
        _SL2RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 2 = degree then
        # N_1(nu), the Steinberg representation.
        rho := SL2IrrepN(2, 1, [0, 0])[1];
        name := "N_1(nu)";
        _SL2RecordIrrep(irrep_list, name, rho, l);
    fi;

    # p = 2, ld = 2.
    l := 4;
    if 3 = degree then
        # D_2(chi)+-, for chi != 1. Note that there is only one non-trivial character,
        # indexed by [1,0]. SL2IrrepD returns a list of the two subreps, + and - .
        rho := SL2IrrepD(2, 2, [1, 0]);
        name := "D_2([1,0])+";
        _SL2RecordIrrep(irrep_list, name, rho[1], l);
        name := "D_2([1,0])-";
        _SL2RecordIrrep(irrep_list, name, rho[2], l);

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
        name := "R_2^0(1,3)_1";
        _SL2RecordIrrep(irrep_list, name, rho, l);
        rho := [-1*rho[1], -1*rho[2]];
        name := "Xi_6 tensor R_2^0(1,3)_1";
        _SL2RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 2 = degree then
        # N_2(chi), for chi primitive, chi^2 != 1.
        # A = <zeta> with ord(zeta) = 6.
        # For this specific case, chi is primitive iff chi(-1) = -1, which leaves only
        # a single relevant character, indexed by [0,1].
        rho := SL2IrrepN(2, 2, [0, 1])[1];
        name := "N_2([0,1])";
        _SL2RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 1 = degree then
        # Xi_3 = C_4, with s = -1, t = i.
        rho := [
            [[-E(4)]],
            [[E(4)]]
        ];
        name := "Xi_3";
        _SL2RecordIrrep(irrep_list, name, rho, l);

        # Xi_9 = C_3, with s = i, t = -i.
        rho := [
            [[E(4)]],
            [[-E(4)]]
        ];
        name := "Xi_9";
        _SL2RecordIrrep(irrep_list, name, rho, l);
    fi;

    # p = 2, ld = 3.
    l := 8;
    if 6 = degree then
        # D_3(chi)+-, for chi primitive.
        # A = <5> x <-1>, and a character is primitive if injective on <5>.
        # This gives two characters, [1,0] and [1,1] (note that both square to 1).
        # SL2IrrepD returns a list of the two subreps, + and - .
        rho := SL2IrrepD(2, 3, [1, 0]);
        name := "D_3([1,0])+";
        _SL2RecordIrrep(irrep_list, name, rho[1], l);
        name := "D_3([1,0])-";
        _SL2RecordIrrep(irrep_list, name, rho[2], l);

        rho := SL2IrrepD(2, 3, [1, 1]);
        name := "D_3([1,1])+";
        _SL2RecordIrrep(irrep_list, name, rho[1], l);
        name := "D_3([1,1])-";
        _SL2RecordIrrep(irrep_list, name, rho[2], l);

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
        name := "R_3^0(1,3,nu)_1";
        _SL2RecordIrrep(irrep_list, name, rho, l);

        rho := [E(4) * rho[1], -E(4) * rho[2]];
        name := "Xi_9 tensor R_3^0(1,3,nu)_1";
        _SL2RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 4 = degree then
        # N_3(chi), for chi primitive, chi^2 != 1.
        # A = <alpha> x <zeta> where ord(alpha) = 2 and ord(zeta) = 6.
        # chi is primitive if injective on <alpha>.
        # This gives two relevant characters, [1,1] and [1,2].
        rho := SL2IrrepN(2, 3, [1, 1])[1];
        name := "N_3([1,1])";
        _SL2RecordIrrep(irrep_list, name, rho, l);

        rho := SL2IrrepN(2, 3, [1, 2])[1];
        name := "N_3([1,2])";
        _SL2RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 2 = degree then
        # N_3(chi)+-, for chi primitive, chi^2 = 1.
        # As above; relevant characters are those which square to 1: [1,0] and [1,3].
        # SL2IrrepN returns a list of the two subreps, + and - .
        rho := SL2IrrepN(2, 3, [1, 0]);
        name := "N_3([1,0])+";
        _SL2RecordIrrep(irrep_list, name, rho[1], l);
        name := "N_3([1,0])-";
        _SL2RecordIrrep(irrep_list, name, rho[2], l);

        rho := SL2IrrepN(2, 3, [1, 3]);
        name := "N_3([1,3])+";
        _SL2RecordIrrep(irrep_list, name, rho[1], l);
        name := "N_3([1,3])-";
        _SL2RecordIrrep(irrep_list, name, rho[2], l);
    fi;
    if 3 = degree then
        # R_3^0(r,t,chi), for chi(zeta) = i and r in {1,3} and t in {1,5}.
        # A = <zeta> with ord(zeta) = 4.
        # Specifically, zeta = (0,1) when t = 1, and zeta = (2,1) when t = 5.
        # The given character is indexed by [0,1].
        for r in [1,3] do
            for t in [1,5] do
                rho := SL2IrrepR(2, 3, 0, r, t, [0,1])[1];
                name := Concatenation("R_3^0(", String(r), ",", String(t), ",[0,1])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
        od;

        # R_3^0(1, t, chi)_+-, for chi != 1 and t in {3,7}.
        # |A| = 2, so there's only one non-trivial chi, indexed by [0,1].
        # chi squares to 1, giving two subreps.
        # SL2IrrepR returns a list with the two supreps, + and -.
        for t in [3,7] do
            rho := SL2IrrepR(2,3,0,1,t,[0,1]);
            name := Concatenation("R_3^0(1,", String(t), ",[0,1])+");
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_3^0(1,", String(t), ",[0,1])-");
            _SL2RecordIrrep(irrep_list, name, rho[2], l);
        od;
    fi;

    # p = 2, ld = 4.
    l := 16;
    if 24 = degree then
        # D_4(chi), for chi primitive, chi^2 != 1.
        # A = <5> x <-1>, and a character is primitive if injective on <5>.
        # Note that 5 has order 4.
        # This gives two relevant characters, [1,0] and [1,1].
        rho := SL2IrrepD(2, 4, [1, 0])[1];
        name := "D_4([1,0])";
        _SL2RecordIrrep(irrep_list, name, rho, l);

        rho := SL2IrrepD(2, 4, [1, 1])[1];
        name := "D_4([1,1])";
        _SL2RecordIrrep(irrep_list, name, rho, l);
    fi;
    if 8 = degree then
        # N_4(chi), for chi primitive, chi^2 != 1.
        # A = <alpha> x <zeta> where ord(alpha) = 4 and ord(zeta) = 6.
        # chi is primitive if injective on <alpha>.
        # This gives six relevant characters: [1,0],[1,3],[1,1],[1,2],[3,1],[3,2].
        for x in [[1,0],[1,3],[1,1],[1,2],[3,1],[3,2]] do
            rho := SL2IrrepN(2, 4, x)[1];
            name := Concatenation("N_4([", String(x[1]), ",", String(x[2]), "])");
            _SL2RecordIrrep(irrep_list, name, rho, l);
        od;
    fi;
    if 6 = degree then
        # R_4^0(r,t,chi), for chi primitive with chi^2 != 1, r in {1,3}, t in {1,5}.
        # For t = 1, A = <alpha> x <zeta>, with ord(alpha) = 2 and ord(zeta) = 4,
        # and a character is primitive iff injective on alpha, so the only relevant
        # character is [1,1].
        # For t = 5, A = <alpha> x <-1>, with ord(alpha) = 4; in this case a character
        # is primitive iff injective on <-alpha^2> (see NW p. 496), which means [1,0].
        for r in [1,3] do
            rho := SL2IrrepR(2, 4, 0, r, 1, [1, 1])[1];
            name := Concatenation("R_4^0(", String(r), ",1,[1,1])");
            _SL2RecordIrrep(irrep_list, name, rho, l);
        od;
        for r in [1,3] do
            rho := SL2IrrepR(2, 4, 0, r, 5, [1, 0])[1];
            name := Concatenation("R_4^0(", String(r), ",5,[1,0])");
            _SL2RecordIrrep(irrep_list, name, rho, l);
        od;

        # R_4^0(1,t,chi)+-, for chi primitive, t in {3,7}.
        # A = <alpha> x <-1>, with ord(alpha) = 2. There are therefore two primitive chars,
        # indexed by [1,0] and [1,1]. Both square to 1, giving two subreps.
        # SL2IrrepR returns a list of the two subreps, + and - .
        for t in [3,7] do
            rho := SL2IrrepR(2, 4, 0, 1, t, [1,0]);
            name := Concatenation("R_4^0(1,", String(t), ",[1,0])+");
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(1,", String(t), ",[1,0])-");
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

            rho := SL2IrrepR(2, 4, 0, 1, t, [1,1]);
            name := Concatenation("R_4^0(1,", String(t), ",[1,1])+");
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(1,", String(t), ",[1,1])-");
            _SL2RecordIrrep(irrep_list, name, rho[2], l);
        od;

        # R_4^2(r,t,chi) and C_2 tensor R_4^2(r,t,chi), for chi != 1, r,t in {1,3}.
        # A = {+-1}, so there is a single non-trivial character, indexed by [0,1].
        # NOTE: C_2 tensor R_4^2(1,1,[0,1]) is iso to R_4^2(1,3,nu)_1, and
        # C_2 tensor R_4^2(3,1,[0,1]) is iso to R_4^2(3,3,nu)_1.
        # NW lists them under the latter name, so this is what we do here; see below.
        for r in [1,3] do
            rho := SL2IrrepR(2,4,2,r,1,[0,1])[1];
            name := Concatenation("R_4^2(", String(r), ",1,[0,1])");
            _SL2RecordIrrep(irrep_list, name, rho, l);

            rho := SL2IrrepR(2,4,2,r,3,[0,1])[1];
            name := Concatenation("R_4^2(", String(r), ",3,[0,1])");
            _SL2RecordIrrep(irrep_list, name, rho, l);

            rho := [-1 * rho[1], -1 * rho[2]];
            name := Concatenation("Xi_6 tensor R_4^2(", String(r), ",3,[0,1])");
            _SL2RecordIrrep(irrep_list, name, rho, l);
        od;

        # R_4^2(r,3,nu)_1, for r in {1,3}.
        # Basis is given on NW p. 524.
        # See above.
        for r in [1,3] do
            rho := SL2IrrepR(2,4,2,r,3,[0,0])[1];
            name := Concatenation("R_4^2(", String(r), ",3,nu)_1");
            _SL2RecordIrrep(irrep_list, name, rho, l);
        od;
    fi;
    if 3 = degree then
        # R_4^0(r,t,chi)+-, for chi primitive, chi^2 = 1, r in {1,3}, t in {1,5}.
        # See previous. For t = 1, the relevant chars. are [1,0] and [1,2].
        # For t = 5, the relevant chars. are [0,1] and [2,1].
        # SL2IrrepR returns a list with the two supreps, + and -.
        for r in [1,3] do
            rho := SL2IrrepR(2, 4, 0, r, 1, [1,0]);
            name := Concatenation("R_4^0(", String(r), ",1,[1,0])+");
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(", String(r), ",1,[1,0])-");
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

            rho := SL2IrrepR(2, 4, 0, r, 1, [1,2]);
            name := Concatenation("R_4^0(", String(r), ",1,[1,2])+");
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(", String(r), ",1,[1,2])-");
            _SL2RecordIrrep(irrep_list, name, rho[2], l);
        od;
        for r in [1,3] do
            rho := SL2IrrepR(2, 4, 0, r, 5, [0,1]);
            name := Concatenation("R_4^0(", String(r), ",5,[0,1])+");
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(", String(r), ",5,[0,1])-");
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

            rho := SL2IrrepR(2, 4, 0, r, 5, [2,1]);
            name := Concatenation("R_4^0(", String(r), ",5,[2,1])+");
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_4^0(", String(r), ",5,[2,1])-");
            _SL2RecordIrrep(irrep_list, name, rho[2], l);
        od;
    fi;
    if 12 = degree then
        # N_3(chi)+ tensor R_4^0(1,7,psi)+, for psi(alpha) = -1, psi(-1) = 1,
        # chi primitive and chi^2 = 1.
        # Calculated via Kronecker product.
        # For the N part, chars. are [1,0] and [1,3] (see the ld=3 section earlier).
        #
        # First construct R_4^0(1,7,psi)+. Char is [1,1]. SL2IrrepR gives a list with + and then -.
        x := SL2IrrepR(2, 4, 0, 1, 7, [1, 1]);
        for j in [0,3] do
            # Construct N_3(chi)+. SL2IrrepR again returns a list.
            rho := SL2IrrepN(2, 3, [1, j]);
            rho := [
                KroneckerProduct(rho[1][1], x[1][1]),
                KroneckerProduct(rho[1][2], x[1][2]),
                12
            ];
            name := Concatenation("N_3([1,", String(j), "])+ tensor R_4^0(1,7,[1,1])+");
            _SL2RecordIrrep(irrep_list, name, rho, l);
        od;
    fi;

    # p = 2, ld = 5.
    l := 32;
    if 48 = degree then
        # D_5(chi), for chi primitive.
        # A = <5> x <-1>, and a character is primitive if injective on <5>.
        # Note that 5 has order 8.
        # This gives four relevant characters, [1,0], [1,1], [3,0], [3,1].
        for i in [1,3] do
            for j in [0,1] do
                rho := SL2IrrepD(2, 5, [i, j])[1];
                name := Concatenation("D_5([", String(i), ",", String(j), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
    fi;
    if 16 = degree then
        # N_5(chi), for chi primitive.
        # A = <alpha> x <zeta> where ord(alpha) = 8 and ord(zeta) = 6.
        # chi is primitive if injective on <alpha>.
        for i in [1,3] do
            for j in [0,3] do
                rho := SL2IrrepN(2, 5, [i,j])[1];
                name := Concatenation("N_5([", String(i), ",", String(j), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
        for i in [1,3,5,7] do
            for j in [1,2] do
                rho := SL2IrrepN(2, 5, [i,j])[1];
                name := Concatenation("N_5([", String(i), ",", String(j), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
    fi;
    if 12 = degree then
        # R_5^0(r,t,chi), for chi primitive, r in {1,3}, t in {1,5}.
        # For t = 1, A = <alpha> x <zeta>, with ord(alpha) = 4 and ord(zeta) = 4.
        # For t = 5, A = <alpha> x <-1>, with ord(alpha) = 8.
        # In both cases, a character is primitive iff injective on alpha.
        for r in [1,3] do
            for j in [0,2] do
                rho := SL2IrrepR(2, 5, 0, r, 1, [1,j])[1];
                name := Concatenation("R_5^0(", String(r), ",1,[1,", String(j), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
            for i in [1,3] do
                rho := SL2IrrepR(2, 5, 0, r, 1, [i,1])[1];
                name := Concatenation("R_5^0(", String(r), ",1,[", String(i), ",1])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
        for r in [1,3] do
            for i in [1,3] do
                for j in [0,1] do
                    rho := SL2IrrepR(2, 5, 0, r, 5, [i,j])[1];
                    name := Concatenation("R_5^0(", String(r), ",5,[", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
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
                    rho := SL2IrrepR(2, 5, 1, r, t, [1,j])[1];
                    name := Concatenation("R_5^1(", String(r), ",", String(t), ",[1,", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
        od;
        for r in [1,3] do
            for t in [3,7] do
                for j in [0,1] do
                    rho := SL2IrrepR(2, 5, 1, r, t, [1,j])[1];
                    name := Concatenation("R_5^1(", String(r), ",", String(t), ",[1,", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
        od;

        # R_5^2(r,1,chi)_1 and C_3 tensor R_5^2(r,1,chi)_1,
        # for chi NOT primitive and r in {1,3}.
        # A = <alpha> x <-1>, so non-primitive characters are indexed by [0,0] and [0,1].
        for r in [1,3] do
            for j in [0,1] do
                rho := SL2IrrepR(2, 5, 2, r, 1, [0,j])[1];
                name := Concatenation("R_5^2(", String(r), ",1,[0,", String(j), "])_1");
                _SL2RecordIrrep(irrep_list, name, rho, l);

                rho := [E(4) * rho[1], -E(4) * rho[2], 12];
                name := Concatenation("Xi_9 tensor R_5^2(", String(r), ",1,[0,", String(j), "])_1");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
    fi;
    if 24 = degree then
        # R_5^0(1,t,chi), for chi primitive and t in {3,7}.
        # A = <alpha> x <-1> with ord(alpha) = 4.
        # A character is primitive iff injective on <alpha>.
        # Relevant chars. are therefore [1,0] and [1,1].
        for t in [3,7] do
            for j in [0,1] do
                rho := SL2IrrepR(2, 5, 0, 1, t, [1,j])[1];
                name := Concatenation("R_5^0(1,", String(t), ",[1,", String(j), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
        od;
    fi;
    if 6 = degree then
        # R_5^2(r,t,chi)+-, for chi primitive, r in {1,3}, and t in {1,3,5,7}.
        # A = <alpha> x <-1>, with ord(alpha) = 2.
        # Relevant chars are therefore [1,0] and [1,1], both of which square to 1.
        # SL2IrrepR returns a list of the two resulting subreps, + and - .
        for r in [1,3] do
            for t in [1,3,5,7] do
                for j in [0,1] do
                    rho := SL2IrrepR(2, 5, 2, r, t, [1,j]);
                    name := Concatenation("R_5^2(", String(r), ",", String(t), ",[1,", String(j), "])+");
                    _SL2RecordIrrep(irrep_list, name, rho[1], l);
                    name := Concatenation("R_5^2(", String(r), ",", String(t), ",[1,", String(j), "])-");
                    _SL2RecordIrrep(irrep_list, name, rho[2], l);
                od;
            od;
        od;
    fi;

    # p = 2, ld > 5.
    ldmax := Maximum([LogInt(degree, 2) + 4, 6]);
    ldset := [6..ldmax];
    for ld in ldset do
        l := 2^ld;
        if 3*2^(ld-1) = degree then
            # D_ld(chi), for chi primitive.
            # A = <5> x <-1>, and a character is primitive if injective on <5>.
            # Note that 5 has order 2^(ld-2).
            for i in PrimeResidues(2^(ld-2) / 2) do
                for j in [0,1] do
                    rho := SL2IrrepD(2, ld, [i, j])[1];
                    name := Concatenation("D_", String(ld), "([", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
        fi;
        if 2^(ld-1) = degree then
            # N_ld(chi), for chi primitive.
            # A = <alpha> x <zeta> where ord(alpha) = 2^(ld-2) and ord(zeta) = 6.
            # chi is primitive if injective on <alpha>.
            for i in PrimeResidues(2^(ld-2) / 2) do
                for j in [0,3] do
                    rho := SL2IrrepN(2, ld, [i,j])[1];
                    name := Concatenation("N_", String(ld), "([", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
            for i in PrimeResidues(2^(ld-2)) do
                for j in [1,2] do
                    rho := SL2IrrepN(2, ld, [i,j])[1];
                    name := Concatenation("N_", String(ld), "([", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
        fi;
        if 3*2^(ld-2) = degree then
            # R_ld^0(1,t,chi), for chi primitive and t in {3,7}.
            # A = <alpha> x <-1> with ord(alpha) = 2^(ld-3).
            # A character is primitive iff injective on <alpha>.
            for t in [3,7] do
                for i in PrimeResidues(2^(ld-3) / 2) do
                    for j in [0,1] do
                        rho := SL2IrrepR(2, ld, 0, 1, t, [i,j])[1];
                        name := Concatenation("R_", String(ld), "^0(1,", String(t), ",[", String(i), ",", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            od;
        fi;
        if 3*2^(ld-3) = degree then
            # R_ld^si(r,t,chi), for chi primitive and:
            # si = 0, r in {1,3}, t = 1 (ord(alpha) = 2^(ld-3), ord(zeta) = 4),
            # si = 0, r in {1,3}, t = 5 (ord(alpha) = 2^(ld-2), ord(zeta) = 2),
            # si = 1, r in {1,5}, t in {1,5} or r in {1,3}, t in {3,7} (ord(alpha) = 2^(ld-3), ord(zeta) = 2),
            # si = 2, r in {1,3}, t in {1,3,5,7} (ord(alpha) = 2^(ld-4), ord(zeta) = 2).
            # In all cases, a character is primitive iff injective on <alpha>.
            for r in [1,3] do
                for i in PrimeResidues(2^(ld-3) / 2) do
                    for j in [0,2] do
                        rho := SL2IrrepR(2, ld, 0, r, 1, [i,j])[1];
                        name := Concatenation("R_", String(ld), "^0(", String(r), ",1,[", String(i), ",", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
                for i in PrimeResidues(2^(ld-3)) do
                    rho := SL2IrrepR(2, ld, 0, r, 1, [i,1])[1];
                    name := Concatenation("R_", String(ld), "^0(", String(r), ",1,[", String(i), ",1])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
            for r in [1,3] do
                for i in PrimeResidues(2^(ld-2) / 2) do
                    for j in [0,1] do
                        rho := SL2IrrepR(2, ld, 0, r, 5, [i,j])[1];
                        name := Concatenation("R_", String(ld), "^0(", String(r), ",5,[", String(i), ",", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            od;
            for r in [1,5] do
                for t in [1,5] do
                    for i in PrimeResidues(2^(ld-3) / 2) do
                        for j in [0,1] do
                            rho := SL2IrrepR(2, ld, 1, r, t, [i,j])[1];
                            name := Concatenation("R_", String(ld), "^1(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                            _SL2RecordIrrep(irrep_list, name, rho, l);
                        od;
                    od;
                od;
            od;
            for r in [1,3] do
                for t in [3,7] do
                    for i in PrimeResidues(2^(ld-3) / 2) do
                        for j in [0,1] do
                            rho := SL2IrrepR(2, ld, 1, r, t, [i,j])[1];
                            name := Concatenation("R_", String(ld), "^1(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                            _SL2RecordIrrep(irrep_list, name, rho, l);
                        od;
                    od;
                od;
            od;
            for r in [1,3] do
                for t in [1,3,5,7] do
                    for i in PrimeResidues(2^(ld-4) / 2) do
                        for j in [0,1] do
                            rho := SL2IrrepR(2, ld, 2, r, t, [i,j])[1];
                            name := Concatenation("R_", String(ld), "^2(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                            _SL2RecordIrrep(irrep_list, name, rho, l);
                        od;
                    od;
                od;
            od;
        fi;
        if 3*2^(ld-4) = degree then
            # R_ld^si(r,t,chi), for chi primitive, si in {3, ..., ld-3}, and r,t in {1,3,5,7}.
            # ord(alpha) = 2^(ld-si-1), ord(zeta) = 2.
            # A character is primitive iff injective on <alpha>.
            for si in [3 .. (ld-3)] do
                for r in [1,3,5,7] do
                    for t in [1,3,5,7] do
                        for i in PrimeResidues(2^(ld-si-1) / 2) do
                            for j in [0,1] do
                                rho := SL2IrrepR(2, ld, si, r, t, [i,j])[1];
                                name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                _SL2RecordIrrep(irrep_list, name, rho, l);
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
                    rho := SL2IrrepR(2, ld, ld-2, r, t, [1,0])[1];
                    name := Concatenation("R_", String(ld), "^", String(ld-2), "(", String(r), ",", String(t), ",[1,0])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);

                    rho := SL2IrrepR(2, ld, ld-2, r, t, [1,1])[1];
                    name := Concatenation("R_", String(ld), "^", String(ld-2), "(", String(r), ",", String(t), ",[1,1])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            if ld = 6 then
                # R_6^4(r,t,nu)_1 and C_2 tensor R_6^4(r,t,nu)_1, for r in {1,3,5,7}, t in {1,3}.
                for r in [1,3,5,7] do
                    for t in [1,3] do
                        rho := SL2IrrepR(2,6,4,r,t,[0,0])[1];
                        name := Concatenation("R_6^4(", String(r), ",", String(t), ",nu)_1");
                        _SL2RecordIrrep(irrep_list, name, rho, l);

                        rho := [-1 * rho[1], -1 * rho[2]];
                        name := Concatenation("Xi_6 tensor R_6^4(", String(r), ",", String(t), ",nu)_1");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            else
                # R_ld^(ld-3)(r,t,chi)_1, for chi in {[0,0],[2,0]} (which both square to 1),
                # r in {1,3,5,7}, t in {1,3}.
                for r in [1,3,5,7] do
                    for t in [1,3] do
                        rho := SL2IrrepR(2,ld,ld-3,r,t,[0,0])[1];
                        name := Concatenation("R_", String(ld), "^", String(ld-3), "(", String(r), ",", String(t), ",nu)_1");
                        _SL2RecordIrrep(irrep_list, name, rho, l);

                        rho := SL2IrrepR(2,ld,ld-3,r,t,[2,0])[1];
                        name := Concatenation("R_", String(ld), "^", String(ld-3), "(", String(r), ",", String(t), ",[2,0])_1");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            fi;
        fi;
    od;

    Info(InfoSL2Reps, 1, "SL2Reps : ", Length(irrep_list), " irreps of degree ", degree, " found.");
    return irrep_list;
end );

InstallGlobalFunction( _SL2IrrepsPPLOfMaxDegree,
function(max_degree)
    local output, degree, count;

    if not max_degree in PositiveIntegers then
        Error("max_degree must be a positive integer.");
    fi;

    output := [];
    count := 0;

    for degree in [1..max_degree] do
        Info(InfoSL2Reps, 1, "SL2Reps : Degree ", degree, ":");
        output[degree] := _SL2IrrepsPPLOfDegree(degree);
        count := count + Length(output[degree]);
    od;

    Info(InfoSL2Reps, 1, "SL2Reps : ", count, " total irreps found.");
    return output;
end );

InstallGlobalFunction( _SL2IrrepsPPLOfLevel,
function(p, ld)
    local irrep_list, l, si, r, t, u, w, x, i, j, name, rho;

    if not IsPrime(p) then
        Error("p must be a prime.");
    elif (not ld in Integers) or (ld < 0) then
        Error("ld must be a non-negative integer.");
    fi;

    l := p^ld;

    irrep_list := [];

    if ld = 0 then
        # The trivial irrep Xi_0 = C_1.
        _SL2RecordIrrep(irrep_list, "Xi_0", [[[1]], [[1]]], 1);
    elif p > 2 then
        if ld = 1 then
            # D_1(chi), for chi primitive, chi^2 != 1.
            # Defined for p >= 3.
            # A = <beta> (alpha is trivial) with ord(beta) = p-1.
            # All non-trivial characters are primitive.
            # Relevant character indices: [0, (1..(p-3)/2)].
            for i in [1 .. (p-3)/2] do
                rho := SL2IrrepD(p, 1, [0, i])[1];
                name := Concatenation("D_1([0,", String(i), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;

            # N_1(chi), for chi primitive.
            # A = <zeta>, ord(zeta) = p+1.
            # Relevant character indices: [0, (1..(p-1)/2)].
            for j in [1 .. (p-1)/2] do
                rho := SL2IrrepN(p, 1, [0, j])[1];
                name := Concatenation("N_1([0,", String(j), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;

            # R_1(1)+- and R_1(u)+-, for u a non-residue.
            # SL2IrrepRUnary gives a list containing two reps, + and - .
            # For p=3, R_1(1)- = Xi_4 and R_1(2)- = Xi_8, the two linear reps. of level 3.
            rho := SL2IrrepRUnary(p, 1, 1);
            name := "R_1(1)+";
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := "R_1(1)-";
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

            u := _SL2QuadNonRes(p);
            rho := SL2IrrepRUnary(p, 1, u);
            name := Concatenation("R_1(", String(u), ")+");
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("R_1(", String(u), ")-");
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

            # N_1(nu), the Steinberg representation.
            rho := SL2IrrepN(p, 1, [0, 0])[1];
            name := "N_1(nu)";
            _SL2RecordIrrep(irrep_list, name, rho, l);
        else
            # ld >= 2.

            # D_ld(chi), for chi primitive, chi^2 != 1.
            # A = <alpha> x <beta>, with ord(alpha) = p^(ld-1) and ord(beta) = p-1.
            # chi is primitive iff it is injective on <alpha>.
            # Remember chi and bar(chi) give same rep; for j = 0 or (p-1)/2,
            # bar(chi(i,j)) = chi(-i, j).
            # We handle that case first:
            for i in PrimeResidues(p^(ld-1)) do
                if i > p^(ld-1) / 2 then
                    break;
                fi;
                rho := SL2IrrepD(p, ld, [i, 0])[1];
                name := Concatenation("D_", String(ld), "([", String(i), ",0])");
                _SL2RecordIrrep(irrep_list, name, rho, l);

                rho := SL2IrrepD(p, ld, [i, (p-1)/2])[1];
                name := Concatenation("D_", String(ld), "([", String(i), ",", String((p-1)/2), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
            # Then the rest:
            for i in PrimeResidues(p^(ld-1)) do
                for j in [1 .. (p-3)/2] do
                    if Gcd(i, p) = 1 then
                        rho := SL2IrrepD(p, ld, [i, j])[1];
                        name := Concatenation("D_", String(ld), "([", String(i), ",", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    fi;
                od;
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
                rho := SL2IrrepN(p, ld, [i, 0])[1];
                name := Concatenation("N_", String(ld), "([", String(i), ",0])");
                _SL2RecordIrrep(irrep_list, name, rho, l);

                rho := SL2IrrepN(p, ld, [i, (p+1)/2])[1];
                name := Concatenation("N_", String(ld), "([", String(i), ",", String((p+1)/2), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
            for i in PrimeResidues(p^(ld-1)) do
                for j in [1 .. (p-1)/2] do
                    rho := SL2IrrepN(p, ld, [i, j])[1];
                    name := Concatenation("N_", String(ld), "([", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            # R_ld^si(r,t,chi), for 1 <= si <= ld-1, r,t either 1 or a non-residue, and chi primitive
            for si in [1 .. ld-1] do
                u := _SL2QuadNonRes(p);
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
                                    rho := SL2IrrepR(p, ld, si, r, t, [i,j])[1];
                                    name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                    _SL2RecordIrrep(irrep_list, name, rho, l);
                                od;
                            od;
                            for i in PrimeResidues(3^(ld-2)) do
                                for j in [1,2] do
                                    rho := SL2IrrepR(p, ld, si, r, t, [i,j])[1];
                                    name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                    _SL2RecordIrrep(irrep_list, name, rho, l);
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
                                    rho := SL2IrrepR(p, ld, si, r, t, [i,j])[1];
                                    name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                    _SL2RecordIrrep(irrep_list, name, rho, l);
                                od;
                            od;
                        fi;
                    od;
                od;
            od;

            # (R_ld(1)+-)_1 and (R_ld(u)+-)_1, for u a non-residue
            rho := SL2IrrepRUnary(p, ld, 1);
            name := Concatenation("(R_", String(ld), "(1)+)_1");
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("(R_", String(ld), "(1)-)_1");
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

            u := _SL2QuadNonRes(p);
            rho := SL2IrrepRUnary(p, ld, u);
            name := Concatenation("(R_", String(ld), "(", String(u), ")+)_1");
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := Concatenation("(R_", String(ld), "(", String(u), ")-)_1");
            _SL2RecordIrrep(irrep_list, name, rho[2], l);
        fi;
    else
        # p = 2
        if ld = 1 then
            # Xi_6 = C_2, with s = t = -1.
            rho := [
                [[-1]],
                [[-1]]
            ];
            name := "Xi_6";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            # N_1(nu), the Steinberg representation.
            rho := SL2IrrepN(2, 1, [0, 0])[1];
            name := "N_1(nu)";
            _SL2RecordIrrep(irrep_list, name, rho, l);
        elif ld = 2 then
            # D_2(chi)+-, for chi != 1. Note that there is only one non-trivial character,
            # indexed by [1,0]. SL2IrrepD returns a list of the two subreps, + and - .
            rho := SL2IrrepD(2, 2, [1, 0]);
            name := "D_2([1,0])+";
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := "D_2([1,0])-";
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

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
            name := "R_2^0(1,3)_1";
            _SL2RecordIrrep(irrep_list, name, rho, l);
            rho := [-1*rho[1], -1*rho[2]];
            name := "Xi_6 tensor R_2^0(1,3)_1";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            # N_2(chi), for chi primitive, chi^2 != 1.
            # A = <zeta> with ord(zeta) = 6.
            # For this specific case, chi is primitive iff chi(-1) = -1, which leaves only
            # a single relevant character, indexed by [0,1].
            rho := SL2IrrepN(2, 2, [0, 1])[1];
            name := "N_2([0,1])";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            # Xi_3 = C_4, with s = -1, t = i.
            rho := [
                [[-E(4)]],
                [[E(4)]]
            ];
            name := "Xi_3";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            # Xi_9 = C_3, with s = i, t = -i.
            rho := [
                [[E(4)]],
                [[-E(4)]]
            ];
            name := "Xi_9";
            _SL2RecordIrrep(irrep_list, name, rho, l);
        elif ld = 3 then
            # D_3(chi)+-, for chi primitive.
            # A = <5> x <-1>, and a character is primitive if injective on <5>.
            # This gives two characters, [1,0] and [1,1] (note that both square to 1).
            # SL2IrrepD returns a list of the two subreps, + and - .
            rho := SL2IrrepD(2, 3, [1, 0]);
            name := "D_3([1,0])+";
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := "D_3([1,0])-";
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

            rho := SL2IrrepD(2, 3, [1, 1]);
            name := "D_3([1,1])+";
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := "D_3([1,1])-";
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

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
            name := "R_3^0(1,3,nu)_1";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            rho := [E(4) * rho[1], -E(4) * rho[2]];
            name := "Xi_9 tensor R_3^0(1,3,nu)_1";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            # N_3(chi), for chi primitive, chi^2 != 1.
            # A = <alpha> x <zeta> where ord(alpha) = 2 and ord(zeta) = 6.
            # chi is primitive if injective on <alpha>.
            # This gives two relevant characters, [1,1] and [1,2].
            rho := SL2IrrepN(2, 3, [1, 1])[1];
            name := "N_3([1,1])";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            rho := SL2IrrepN(2, 3, [1, 2])[1];
            name := "N_3([1,2])";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            # N_3(chi)+-, for chi primitive, chi^2 = 1.
            # As above; relevant characters are those which square to 1: [1,0] and [1,3].
            # SL2IrrepN returns a list of the two subreps, + and - .
            rho := SL2IrrepN(2, 3, [1, 0]);
            name := "N_3([1,0])+";
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := "N_3([1,0])-";
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

            rho := SL2IrrepN(2, 3, [1, 3]);
            name := "N_3([1,3])+";
            _SL2RecordIrrep(irrep_list, name, rho[1], l);
            name := "N_3([1,3])-";
            _SL2RecordIrrep(irrep_list, name, rho[2], l);

            # R_3^0(r,t,chi), for chi(zeta) = i and r in {1,3} and t in {1,5}.
            # A = <zeta> with ord(zeta) = 4.
            # Specifically, zeta = (0,1) when t = 1, and zeta = (2,1) when t = 5.
            # The given character is indexed by [0,1].
            for r in [1,3] do
                for t in [1,5] do
                    rho := SL2IrrepR(2, 3, 0, r, t, [0,1])[1];
                    name := Concatenation("R_3^0(", String(r), ",", String(t), ",[0,1])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            # R_3^0(1, t, chi)_+-, for chi != 1 and t in {3,7}.
            # |A| = 2, so there's only one non-trivial chi, indexed by [0,1].
            # chi squares to 1, giving two subreps.
            # SL2IrrepR returns a list with the two subreps, + and -.
            for t in [3,7] do
                rho := SL2IrrepR(2,3,0,1,t,[0,1]);
                name := Concatenation("R_3^0(1,", String(t), ",[0,1])+");
                _SL2RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("R_3^0(1,", String(t), ",[0,1])-");
                _SL2RecordIrrep(irrep_list, name, rho[2], l);
            od;
        elif ld = 4 then
            # D_4(chi), for chi primitive, chi^2 != 1.
            # A = <5> x <-1>, and a character is primitive if injective on <5>.
            # Note that 5 has order 4.
            # This gives two relevant characters, [1,0] and [1,1].
            rho := SL2IrrepD(2, 4, [1, 0])[1];
            name := "D_4([1,0])";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            rho := SL2IrrepD(2, 4, [1, 1])[1];
            name := "D_4([1,1])";
            _SL2RecordIrrep(irrep_list, name, rho, l);

            # N_4(chi), for chi primitive, chi^2 != 1.
            # A = <alpha> x <zeta> where ord(alpha) = 4 and ord(zeta) = 6.
            # chi is primitive if injective on <alpha>.
            # This gives six relevant characters: [1,0],[1,3],[1,1],[1,2],[3,1],[3,2].
            for x in [[1,0],[1,3],[1,1],[1,2],[3,1],[3,2]] do
                rho := SL2IrrepN(2, 4, x)[1];
                name := Concatenation("N_4([", String(x[1]), ",", String(x[2]), "])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;

            # R_4^0(r,t,chi), for chi primitive with chi^2 != 1, r in {1,3}, t in {1,5}.
            # For t = 1, A = <alpha> x <zeta>, with ord(alpha) = 2 and ord(zeta) = 4,
            # and a character is primitive iff injective on alpha, so the only relevant
            # character is [1,1].
            # For t = 5, A = <alpha> x <-1>, with ord(alpha) = 4; in this case a character
            # is primitive iff injective on <-alpha^2> (see NW p. 496), which means [1,0].
            for r in [1,3] do
                rho := SL2IrrepR(2, 4, 0, r, 1, [1, 1])[1];
                name := Concatenation("R_4^0(", String(r), ",1,[1,1])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
            for r in [1,3] do
                rho := SL2IrrepR(2, 4, 0, r, 5, [1, 0])[1];
                name := Concatenation("R_4^0(", String(r), ",5,[1,0])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;

            # R_4^0(1,t,chi)+-, for chi primitive, t in {3,7}.
            # A = <alpha> x <-1>, with ord(alpha) = 2. There are therefore two primitive chars,
            # indexed by [1,0] and [1,1]. Both square to 1, giving two subreps.
            # SL2IrrepR returns a list of the two subreps, + and - .
            for t in [3,7] do
                rho := SL2IrrepR(2, 4, 0, 1, t, [1,0]);
                name := Concatenation("R_4^0(1,", String(t), ",[1,0])+");
                _SL2RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("R_4^0(1,", String(t), ",[1,0])-");
                _SL2RecordIrrep(irrep_list, name, rho[2], l);

                rho := SL2IrrepR(2, 4, 0, 1, t, [1,1]);
                name := Concatenation("R_4^0(1,", String(t), ",[1,1])+");
                _SL2RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("R_4^0(1,", String(t), ",[1,1])-");
                _SL2RecordIrrep(irrep_list, name, rho[2], l);
            od;

            # R_4^2(r,t,chi) and C_2 tensor R_4^2(r,t,chi), for chi != 1, r,t in {1,3}.
            # A = {+-1}, so there is a single non-trivial character, indexed by [0,1].
            # NOTE: C_2 tensor R_4^2(1,1,[0,1]) is iso to R_4^2(1,3,nu)_1, and
            # C_2 tensor R_4^2(3,1,[0,1]) is iso to R_4^2(3,3,nu)_1.
            # NW lists them under the latter name, so this is what we do here; see below.
            for r in [1,3] do
                rho := SL2IrrepR(2,4,2,r,1,[0,1])[1];
                name := Concatenation("R_4^2(", String(r), ",1,[0,1])");
                _SL2RecordIrrep(irrep_list, name, rho, l);

                rho := SL2IrrepR(2,4,2,r,3,[0,1])[1];
                name := Concatenation("R_4^2(", String(r), ",3,[0,1])");
                _SL2RecordIrrep(irrep_list, name, rho, l);

                rho := [-1 * rho[1], -1 * rho[2]];
                name := Concatenation("Xi_6 tensor R_4^2(", String(r), ",3,[0,1])");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;

            # R_4^2(r,3,nu)_1, for r in {1,3}.
            # Basis is given on NW p. 524.
            # See above.
            for r in [1,3] do
                rho := SL2IrrepR(2,4,2,r,3,[0,0])[1];
                name := Concatenation("R_4^2(", String(r), ",3,nu)_1");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;

            # R_4^0(r,t,chi)+-, for chi primitive, chi^2 = 1, r in {1,3}, t in {1,5}.
            # See previous. For t = 1, the relevant chars. are [1,0] and [1,2].
            # For t = 5, the relevant chars. are [0,1] and [2,1].
            # SL2IrrepR returns a list with the two supreps, + and -.
            for r in [1,3] do
                rho := SL2IrrepR(2, 4, 0, r, 1, [1,0]);
                name := Concatenation("R_4^0(", String(r), ",1,[1,0])+");
                _SL2RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("R_4^0(", String(r), ",1,[1,0])-");
                _SL2RecordIrrep(irrep_list, name, rho[2], l);

                rho := SL2IrrepR(2, 4, 0, r, 1, [1,2]);
                name := Concatenation("R_4^0(", String(r), ",1,[1,2])+");
                _SL2RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("R_4^0(", String(r), ",1,[1,2])-");
                _SL2RecordIrrep(irrep_list, name, rho[2], l);
            od;
            for r in [1,3] do
                rho := SL2IrrepR(2, 4, 0, r, 5, [0,1]);
                name := Concatenation("R_4^0(", String(r), ",5,[0,1])+");
                _SL2RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("R_4^0(", String(r), ",5,[0,1])-");
                _SL2RecordIrrep(irrep_list, name, rho[2], l);

                rho := SL2IrrepR(2, 4, 0, r, 5, [2,1]);
                name := Concatenation("R_4^0(", String(r), ",5,[2,1])+");
                _SL2RecordIrrep(irrep_list, name, rho[1], l);
                name := Concatenation("R_4^0(", String(r), ",5,[2,1])-");
                _SL2RecordIrrep(irrep_list, name, rho[2], l);
            od;

            # N_3(chi)+ tensor R_4^0(1,7,psi)+, for psi(alpha) = -1, psi(-1) = 1,
            # chi primitive and chi^2 = 1.
            # Calculated via Kronecker product.
            # For the N part, chars. are [1,0] and [1,3] (see the ld=3 section earlier).
            #
            # First construct R_4^0(1,7,psi)+. Char is [1,1]. SL2IrrepR gives a list with + and then -.
            x := SL2IrrepR(2, 4, 0, 1, 7, [1, 1]);
            for j in [0,3] do
                # Construct N_3(chi)+. SL2IrrepR again returns a list.
                rho := SL2IrrepN(2, 3, [1, j]);
                rho := [
                    KroneckerProduct(rho[1][1], x[1][1]),
                    KroneckerProduct(rho[1][2], x[1][2]),
                    12
                ];
                name := Concatenation("N_3([1,", String(j), "])+ tensor R_4^0(1,7,[1,1])+");
                _SL2RecordIrrep(irrep_list, name, rho, l);
            od;
        elif ld = 5 then
            # D_5(chi), for chi primitive.
            # A = <5> x <-1>, and a character is primitive if injective on <5>.
            # Note that 5 has order 8.
            # This gives four relevant characters, [1,0], [1,1], [3,0], [3,1].
            for i in [1,3] do
                for j in [0,1] do
                    rho := SL2IrrepD(2, 5, [i, j])[1];
                    name := Concatenation("D_5([", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            # N_5(chi), for chi primitive.
            # A = <alpha> x <zeta> where ord(alpha) = 8 and ord(zeta) = 6.
            # chi is primitive if injective on <alpha>.
            for i in [1,3] do
                for j in [0,3] do
                    rho := SL2IrrepN(2, 5, [i,j])[1];
                    name := Concatenation("N_5([", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
            for i in [1,3,5,7] do
                for j in [1,2] do
                    rho := SL2IrrepN(2, 5, [i,j])[1];
                    name := Concatenation("N_5([", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            # R_5^0(r,t,chi), for chi primitive, r in {1,3}, t in {1,5}.
            # For t = 1, A = <alpha> x <zeta>, with ord(alpha) = 4 and ord(zeta) = 4.
            # For t = 5, A = <alpha> x <-1>, with ord(alpha) = 8.
            # In both cases, a character is primitive iff injective on alpha.
            for r in [1,3] do
                for j in [0,2] do
                    rho := SL2IrrepR(2, 5, 0, r, 1, [1,j])[1];
                    name := Concatenation("R_5^0(", String(r), ",1,[1,", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
                for i in [1,3] do
                    rho := SL2IrrepR(2, 5, 0, r, 1, [i,1])[1];
                    name := Concatenation("R_5^0(", String(r), ",1,[", String(i), ",1])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
            for r in [1,3] do
                for i in [1,3] do
                    for j in [0,1] do
                        rho := SL2IrrepR(2, 5, 0, r, 5, [i,j])[1];
                        name := Concatenation("R_5^0(", String(r), ",5,[", String(i), ",", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
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
                        rho := SL2IrrepR(2, 5, 1, r, t, [1,j])[1];
                        name := Concatenation("R_5^1(", String(r), ",", String(t), ",[1,", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            od;
            for r in [1,3] do
                for t in [3,7] do
                    for j in [0,1] do
                        rho := SL2IrrepR(2, 5, 1, r, t, [1,j])[1];
                        name := Concatenation("R_5^1(", String(r), ",", String(t), ",[1,", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            od;

            # R_5^2(r,1,chi)_1 and C_3 tensor R_5^2(r,1,chi)_1,
            # for chi NOT primitive and r in {1,3}.
            # A = <alpha> x <-1>, so non-primitive characters are indexed by [0,0] and [0,1].
            for r in [1,3] do
                for j in [0,1] do
                    rho := SL2IrrepR(2, 5, 2, r, 1, [0,j])[1];
                    name := Concatenation("R_5^2(", String(r), ",1,[0,", String(j), "])_1");
                    _SL2RecordIrrep(irrep_list, name, rho, l);

                    rho := [E(4) * rho[1], -E(4) * rho[2], 12];
                    name := Concatenation("Xi_9 tensor R_5^2(", String(r), ",1,[0,", String(j), "])_1");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            # R_5^0(1,t,chi), for chi primitive and t in {3,7}.
            # A = <alpha> x <-1> with ord(alpha) = 4.
            # A character is primitive iff injective on <alpha>.
            # Relevant chars. are therefore [1,0] and [1,1].
            for t in [3,7] do
                for j in [0,1] do
                    rho := SL2IrrepR(2, 5, 0, 1, t, [1,j])[1];
                    name := Concatenation("R_5^0(1,", String(t), ",[1,", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            # R_5^2(r,t,chi)+-, for chi primitive, r in {1,3}, and t in {1,3,5,7}.
            # A = <alpha> x <-1>, with ord(alpha) = 2.
            # Relevant chars are therefore [1,0] and [1,1], both of which square to 1.
            # SL2IrrepR returns a list of the two resulting subreps, + and - .
            for r in [1,3] do
                for t in [1,3,5,7] do
                    for j in [0,1] do
                        rho := SL2IrrepR(2, 5, 2, r, t, [1,j]);
                        name := Concatenation("R_5^2(", String(r), ",", String(t), ",[1,", String(j), "])+");
                        _SL2RecordIrrep(irrep_list, name, rho[1], l);
                        name := Concatenation("R_5^2(", String(r), ",", String(t), ",[1,", String(j), "])-");
                        _SL2RecordIrrep(irrep_list, name, rho[2], l);
                    od;
                od;
            od;
        else
            # ld > 5

            # D_ld(chi), for chi primitive.
            # A = <5> x <-1>, and a character is primitive if injective on <5>.
            # Note that 5 has order 2^(ld-2).
            for i in PrimeResidues(2^(ld-2) / 2) do
                for j in [0,1] do
                    rho := SL2IrrepD(2, ld, [i, j])[1];
                    name := Concatenation("D_", String(ld), "([", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            # N_ld(chi), for chi primitive.
            # A = <alpha> x <zeta> where ord(alpha) = 2^(ld-2) and ord(zeta) = 6.
            # chi is primitive if injective on <alpha>.
            for i in PrimeResidues(2^(ld-2) / 2) do
                for j in [0,3] do
                    rho := SL2IrrepN(2, ld, [i,j])[1];
                    name := Concatenation("N_", String(ld), "([", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
            for i in PrimeResidues(2^(ld-2)) do
                for j in [1,2] do
                    rho := SL2IrrepN(2, ld, [i,j])[1];
                    name := Concatenation("N_", String(ld), "([", String(i), ",", String(j), "])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            # R_ld^0(1,t,chi), for chi primitive and t in {3,7}.
            # A = <alpha> x <-1> with ord(alpha) = 2^(ld-3).
            # A character is primitive iff injective on <alpha>.
            for t in [3,7] do
                for i in PrimeResidues(2^(ld-3) / 2) do
                    for j in [0,1] do
                        rho := SL2IrrepR(2, ld, 0, 1, t, [i,j])[1];
                        name := Concatenation("R_", String(ld), "^0(1,", String(t), ",[", String(i), ",", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
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
                        rho := SL2IrrepR(2, ld, 0, r, 1, [i,j])[1];
                        name := Concatenation("R_", String(ld), "^0(", String(r), ",1,[", String(i), ",", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
                for i in PrimeResidues(2^(ld-3)) do
                    rho := SL2IrrepR(2, ld, 0, r, 1, [i,1])[1];
                    name := Concatenation("R_", String(ld), "^0(", String(r), ",1,[", String(i), ",1])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;
            for r in [1,3] do
                for i in PrimeResidues(2^(ld-2) / 2) do
                    for j in [0,1] do
                        rho := SL2IrrepR(2, ld, 0, r, 5, [i,j])[1];
                        name := Concatenation("R_", String(ld), "^0(", String(r), ",5,[", String(i), ",", String(j), "])");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            od;
            for r in [1,5] do
                for t in [1,5] do
                    for i in PrimeResidues(2^(ld-3) / 2) do
                        for j in [0,1] do
                            rho := SL2IrrepR(2, ld, 1, r, t, [i,j])[1];
                            name := Concatenation("R_", String(ld), "^1(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                            _SL2RecordIrrep(irrep_list, name, rho, l);
                        od;
                    od;
                od;
            od;
            for r in [1,3] do
                for t in [3,7] do
                    for i in PrimeResidues(2^(ld-3) / 2) do
                        for j in [0,1] do
                            rho := SL2IrrepR(2, ld, 1, r, t, [i,j])[1];
                            name := Concatenation("R_", String(ld), "^1(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                            _SL2RecordIrrep(irrep_list, name, rho, l);
                        od;
                    od;
                od;
            od;
            for r in [1,3] do
                for t in [1,3,5,7] do
                    for i in PrimeResidues(2^(ld-4) / 2) do
                        for j in [0,1] do
                            rho := SL2IrrepR(2, ld, 2, r, t, [i,j])[1];
                            name := Concatenation("R_", String(ld), "^2(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                            _SL2RecordIrrep(irrep_list, name, rho, l);
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
                                rho := SL2IrrepR(2, ld, si, r, t, [i,j])[1];
                                name := Concatenation("R_", String(ld), "^", String(si), "(", String(r), ",", String(t), ",[", String(i), ",", String(j), "])");
                                _SL2RecordIrrep(irrep_list, name, rho, l);
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
                    rho := SL2IrrepR(2, ld, ld-2, r, t, [1,0])[1];
                    name := Concatenation("R_", String(ld), "^", String(ld-2), "(", String(r), ",", String(t), ",[1,0])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);

                    rho := SL2IrrepR(2, ld, ld-2, r, t, [1,1])[1];
                    name := Concatenation("R_", String(ld), "^", String(ld-2), "(", String(r), ",", String(t), ",[1,1])");
                    _SL2RecordIrrep(irrep_list, name, rho, l);
                od;
            od;

            if ld = 6 then
                # R_6^4(r,t,nu)_1 and C_2 tensor R_6^4(r,t,nu)_1, for r in {1,3,5,7}, t in {1,3}.
                for r in [1,3,5,7] do
                    for t in [1,3] do
                        rho := SL2IrrepR(2,6,4,r,t,[0,0])[1];
                        name := Concatenation("R_6^4(", String(r), ",", String(t), ",nu)_1");
                        _SL2RecordIrrep(irrep_list, name, rho, l);

                        rho := [-1 * rho[1], -1 * rho[2]];
                        name := Concatenation("Xi_6 tensor R_6^4(", String(r), ",", String(t), ",nu)_1");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            else
                # R_ld^(ld-3)(r,t,chi)_1, for chi in {[0,0],[2,0]} (which both square to 1),
                # r in {1,3,5,7}, t in {1,3}.
                for r in [1,3,5,7] do
                    for t in [1,3] do
                        rho := SL2IrrepR(2,ld,ld-3,r,t,[0,0])[1];
                        name := Concatenation("R_", String(ld), "^", String(ld-3), "(", String(r), ",", String(t), ",nu)_1");
                        _SL2RecordIrrep(irrep_list, name, rho, l);

                        rho := SL2IrrepR(2,ld,ld-3,r,t,[2,0])[1];
                        name := Concatenation("R_", String(ld), "^", String(ld-3), "(", String(r), ",", String(t), ",[2,0])_1");
                        _SL2RecordIrrep(irrep_list, name, rho, l);
                    od;
                od;
            fi;
        fi;
    fi;

    Info(InfoSL2Reps, 1, "SL2Reps : ", Length(irrep_list), " irreps of level ", l, " found.");
    return irrep_list;
end );

#

InstallGlobalFunction( SL2IrrepsOfDegree,
function(degree)
    local linears, prime_power_reps, i, ConstructIrreps, triv, factorizations, output, f;

    if not degree in PositiveIntegers then
        Error("degree must be a positive integer.");
    fi;

    prime_power_reps := [];

    # collect prime-power-level irreps of degree dividing the given degree
    Info(InfoSL2Reps, 1, "SL2Reps : Constructing irreps of prime-power level.");

    Info(InfoSL2Reps, 1, "SL2Reps : Degree 1:");

    # The linear reps are denoted Xi_n, n in Z/12Z, with T = [zeta_12^n] and S = [i^n].
    # We handle these separately just for brevity in the output.
    prime_power_reps[1] := [];

    _SL2RecordIrrep(prime_power_reps[1], "Xi_0", [[[1]], [[1]]], 1); # Xi_0 = C_1

    _SL2RecordIrrep(prime_power_reps[1], "Xi_6", [[[-1]], [[-1]]], 2); # Xi_6 = C_2

    _SL2RecordIrrep(prime_power_reps[1], "Xi_4", [[[1]], [[E(3)]]], 3); # Xi_4 = R_1(1)- with p=3
    _SL2RecordIrrep(prime_power_reps[1], "Xi_8", [[[1]], [[E(3)^2]]], 3); # Xi_8 = R_1(2)- with p=3

    _SL2RecordIrrep(prime_power_reps[1], "Xi_3", [[[-E(4)]], [[E(4)]]], 4); # Xi_3 = C_4
    _SL2RecordIrrep(prime_power_reps[1], "Xi_9", [[[E(4)]], [[-E(4)]]], 4); # Xi_9 = C_3

    _SL2RecordIrrep(prime_power_reps[1], "Xi_2", [[[-1]], [[-E(3)^2]]], 6); # = Xi_6 * Xi_8
    _SL2RecordIrrep(prime_power_reps[1], "Xi_10", [[[-1]], [[-E(3)]]], 6); # = Xi_6 * Xi_4

    _SL2RecordIrrep(prime_power_reps[1], "Xi_1", [[[E(4)]], [[E(12)]]], 12); # = Xi_9 * Xi_4
    _SL2RecordIrrep(prime_power_reps[1], "Xi_5", [[[E(4)]], [[E(12)^5]]], 12); # = Xi_9 * Xi_8
    _SL2RecordIrrep(prime_power_reps[1], "Xi_7", [[[-E(4)]], [[E(12)^7]]], 12); # = Xi_3 * Xi_4
    _SL2RecordIrrep(prime_power_reps[1], "Xi_11", [[[-E(4)]], [[E(12)^11]]], 12); # = Xi_3 * Xi_8

    Info(InfoSL2Reps, 1, "SL2Reps : 12 irreps of degree 1 found.");

    for i in DivisorsInt(degree) do
        if i = 1 then
            continue;
        fi;
        Info(InfoSL2Reps, 1, "SL2Reps : Degree ", i, ":");
        prime_power_reps[i] := _SL2IrrepsPPLOfDegree(i);
    od;

    Info(InfoSL2Reps, 1, "SL2Reps : Constructing tensor products.");

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

            if Gcd(eta.level, rho.level) = 1 then
                if Length(factors) > 1 and factors[2] = factors[1] then
                    # Next degree is same as current. Need to keep irreps in order.
                    new_start := i+1;
                else
                    new_start := 1;
                fi;

                output := Concatenation(
                    output,
                    ConstructIrreps(
                        rec(
                            S := KroneckerProduct(rho.S, eta.S),
                            T := KroneckerProduct(rho.T, eta.T),
                            degree := rho.degree * eta.degree,
                            level := rho.level * eta.level, # level is Lcm(level(rho), level(eta)), but they're coprime.
                            name := _SL2ConcatNames(rho.name, eta.name)
                        ),
                        factors{[2..Length(factors)]},
                        new_start
                    )
                );
            fi;
        od;

        return output;
    end;

    triv := rec(
        S := [[1]],
        T := [[1]],
        degree := 1,
        level := 1,
        name := "Xi_0"
    );

    factorizations := _SL2Factorizations(degree);
    output := [];
    for f in factorizations do
        Append(output, ConstructIrreps(triv, f, 1));
    od;

    # Sort by level.
    SortBy(output, x -> x.level);

    for i in [1 .. Length(output)] do
        Info(InfoSL2Reps, 1, "SL2Reps : ", i, ": ( ", output[i].name, " ) [d: ", output[i].degree, ", l: ", output[i].level ,"]");
    od;
    Info(InfoSL2Reps, 1, "SL2Reps : Total count: ", Length(output));
    return output;
end );

InstallGlobalFunction( SL2IrrepsOfMaxDegree,
function(max_degree)
    local output, degree;

    if not max_degree in PositiveIntegers then
        Error("max_degree must be a positive integer.");
    fi;

    output := [];

    for degree in [1..max_degree] do
        Info(InfoSL2Reps, 1, "SL2Reps : Degree ", degree, ":");
        Append(output, SL2IrrepsOfDegree(degree));
    od;

    Info(InfoSL2Reps, 1, "SL2Reps : ", Length(output), " total irreps found.");
    return output;
end );

InstallGlobalFunction( SL2IrrepsOfLevel,
function(l)
    local output, pp_lists, factors, i, rho, eta;

    if (not l in Integers) or (l <= 0) then
        Error("level must be a positive integer.");
    fi;

    if l = 1 then
        output := [];

        Info(InfoSL2Reps, 1, "SL2Reps : At level 1, there is only the trivial irrep.");
        return [
            rec(
                S := [[1]],
                T := [[1]],
                degree := 1,
                level := 1,
                name := "Xi_0"
            )
        ];
    else
        factors := PrimePowersInt(l);

        if Length(factors) = 2 then
            # level is a prime power
            Info(InfoSL2Reps, 1, "SL2Reps : Constructing irreps of level ", l, ", which is a prime power.");
            return _SL2IrrepsPPLOfLevel(factors[1], factors[2]);
        else
            # tensor irreps of the prime power factors
            Info(InfoSL2Reps, 1, "SL2Reps : Constructing irreps of prime-power level.");
            pp_lists := [];
            for i in [1..(Length(factors) / 2)] do
                Info(InfoSL2Reps, 1, "SL2Reps : Level ", (factors[2*i-1])^(factors[2*i]), ":");
                pp_lists[i] := _SL2IrrepsPPLOfLevel(factors[2*i-1], factors[2*i]);
            od;

            # Possible alternate approach: make lists of the indices and take the Cartesian
            # thereof instead.
            Info(InfoSL2Reps, 1, "SL2Reps : Constructing tensor products.");
            pp_lists := Cartesian(pp_lists);

            output := [];

            for i in pp_lists do
                rho := rec(
                    S := [[1]],
                    T := [[1]],
                    degree := 1,
                    level := 1,
                    name := "Xi_0"
                );

                for eta in i do
                    rho := rec(
                        S := KroneckerProduct(rho.S, eta.S),
                        T := KroneckerProduct(rho.T, eta.T),
                        degree := rho.degree * eta.degree,
                        level := rho.level * eta.level, # level is Lcm(level(rho), level(eta)), but they're coprime.
                        name := _SL2ConcatNames(rho.name,eta.name)
                    );
                od;

                Add(output, rho);
            od;

            # Sort by degree.
            SortBy(output, x -> x.degree);

            for i in [1 .. Length(output)] do
                Info(InfoSL2Reps, 1, "SL2Reps : ", i, ": ( ", output[i].name, " ) [d: ", output[i].degree, ", l: ", output[i].level ,"]");
            od;
            Info(InfoSL2Reps, 1, "SL2Reps : Total count: ", Length(output));
            return output;
        fi;
    fi;
end );

InstallGlobalFunction( SL2IrrepsExceptional,
function()
    local irrep_list, rho, name, w, x, j, r, t;

    irrep_list := [];

    # C_2 tensor R_2^0(1,3)_1.
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
    rho := [-1*rho[1], -1*rho[2]];
    name := "Xi_6 tensor R_2^0(1,3)_1";
    _SL2RecordIrrep(irrep_list, name, rho, 2^2);

    # C_3 tensor R_3^0(1,3,nu)_1.
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
    rho := [E(4) * rho[1], -E(4) * rho[2]];
    name := "Xi_9 tensor R_3^0(1,3,nu)_1";
    _SL2RecordIrrep(irrep_list, name, rho, 2^3);

    # C_2 tensor R_4^2(r,3,chi), for chi != 1, r in {1,3}.
    # A = {+-1}, so there is a single non-trivial character, indexed by [0,1].
    # NOTE: C_2 tensor R_4^2(1,1,[0,1]) is iso to R_4^2(1,3,nu)_1, and
    # C_2 tensor R_4^2(3,1,[0,1]) is iso to R_4^2(3,3,nu)_1.
    # Hence they are not exceptional.
    for r in [1,3] do
        rho := SL2IrrepR(2,4,2,r,3,[0,1])[1];
        rho := [-1 * rho[1], -1 * rho[2]];
        name := Concatenation("Xi_6 tensor R_4^2(", String(r), ",3,[0,1])");
        _SL2RecordIrrep(irrep_list, name, rho, 2^4);
    od;

    # N_3(chi)+ tensor R_4^0(1,7,psi)+, for psi(alpha) = -1, psi(-1) = 1,
    # chi primitive and chi^2 = 1.
    # Calculated via Kronecker product.
    # For the N part, chars. are [1,0] and [1,3] (see the ld=3 section earlier).
    #
    # First construct R_4^0(1,7,psi)+. Char is [1,1]. SL2IrrepR gives a list with + and then -.
    x := SL2IrrepR(2, 4, 0, 1, 7, [1, 1]);
    for j in [0,3] do
        # Construct N_3(chi)+. SL2IrrepR again returns a list.
        rho := SL2IrrepN(2, 3, [1, j]);
        rho := [
            KroneckerProduct(rho[1][1], x[1][1]),
            KroneckerProduct(rho[1][2], x[1][2]),
            12
        ];
        name := Concatenation("N_3([1,", String(j), "])+ tensor R_4^0(1,7,[1,1])+");
        _SL2RecordIrrep(irrep_list, name, rho, 2^4);
    od;

    # C_3 tensor R_5^2(r,1,chi)_1, for chi NOT primitive and r in {1,3}.
    # A = <alpha> x <-1>, so non-primitive characters are indexed by [0,0] and [0,1].
    for r in [1,3] do
        for j in [0,1] do
            rho := SL2IrrepR(2, 5, 2, r, 1, [0,j])[1];
            rho := [E(4) * rho[1], -E(4) * rho[2], 12];
            name := Concatenation("Xi_9 tensor R_5^2(", String(r), ",1,[0,", String(j), "])_1");
            _SL2RecordIrrep(irrep_list, name, rho, 2^5);
        od;
    od;

    # C_2 tensor R_6^4(r,t,nu)_1, for r in {1,3,5,7}, t in {1,3}.
    for r in [1,3,5,7] do
        for t in [1,3] do
            rho := SL2IrrepR(2,6,4,r,t,[0,0])[1];
            rho := [-1 * rho[1], -1 * rho[2]];
            name := Concatenation("Xi_6 tensor R_6^4(", String(r), ",", String(t), ",nu)_1");
            _SL2RecordIrrep(irrep_list, name, rho, 2^6);
        od;
    od;

    Info(InfoSL2Reps, 1, "SL2Reps : ", Length(irrep_list), " exceptional irreps found.");
    return irrep_list;
end );
