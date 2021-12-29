#
# SL2Reps: Constructing symmetric representations of SL(2,Z).
#
# Miscellaneous and helper functions.
#
# Implementations
#


#----------------------------------------
# Returns a square root of the given root of unity.
#---------------------------------------
InstallGlobalFunction( _SL2SqrtOfRootOfUnity,
function(b)
    local desc;

    if b = 1 then
        return 1;
    else
        # this returns [q,p] such that b = E(q)^p
        desc := DescriptionOfRootOfUnity(b);

        return E(desc[1]*2)^desc[2];
    fi;
end );

#----------------------------------------
# Returns a quadratic non-residue mod p.
#---------------------------------------
InstallGlobalFunction( _SL2QuadNonRes,
function(p)
    local j;

    j := 2;
    while Jacobi(j, p) = 1 do
        j := j + 1;
    od;
    return j;
end );

#----------------------------------------
# Returns a list of all factorizations of n.
#---------------------------------------
InstallGlobalFunction( _SL2Factorizations,
function(n)
    local o, f, Step;

    # recursion step
    Step := function(factors)
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

            for r in Step(rest) do
                Add(r, m);
                Add(o, r);
            od;
        od;

        return o;
    end;

    if n = 1 then
        return [[1]];
    fi;

    o := [];

    for f in Step(Factors(n)) do
        Sort(f);
        f := Reversed(f);
        Add(f, 1);
        if not f in o then
            Add(o,f);
        fi;
    od;

    Sort(o);
    return Reversed(o);
end );

#----------------------------------------
# Records an irrep into irrep_list.
#---------------------------------------
InstallGlobalFunction( _SL2RecordIrrep,
function(irrep_list, name, rho, l)
    Info(InfoSL2Reps, 1, "SL2Reps : ", name, " [level ", l, "]");
    # record rho in the form [S, T, degree, level, name]
    Add(irrep_list, rec(
        S := rho[1],
        T := rho[2],
        degree := Length(rho[1]),
        level := l,
        name := name
    ));
end );

#----------------------------------------
# Constructs a name for a tensor product.
#---------------------------------------
InstallGlobalFunction( _SL2ConcatNames,
function(a,b)
    if a = "Xi_0" then
        return b;
    elif b = "Xi_0" then
        return a;
    else
        return Concatenation(a, " tensor ", b);
    fi;
end );

#----------------------------------------
# Returns conj class representative information of SL(2, Z/p^ld Z) for p (odd prime) and ld (>= 1).
# Each conj class representative is represented by a list of length 2. To rewrite the conj class rep as products of s and t (NW notation) corresponding to a, we can compute
# Prod(a[1], n -> t^n*s) * s^(2*a[2]).
# Note that s^2 is central and can be computed easily.
# Update: 1. Merged cases, deleted size of the conj classes; 2. Simplified the expression of the word presentation of the representatives; 3. Sorted the output so that id and -id are the first two conj classes.
#---------------------------------------
InstallGlobalFunction( _SL2ConjClassesOdd,
function(p, ld)
    local l, u, uinv, Ald, aset1, type11, type1u, type12, type13, Type1, type2h, Type2, i;

    l := p^ld;

    # Find a quad non-residue u and its mult inverse mod l.
    u:=First([1..p], i->Jacobi(i,p)=-1);
    uinv := u^-1 mod l;

    # There are 2 types of conj classes.

    # Type 1.   There are 2 cases depending on whether the upper left entry is 1 or u.
    Ald := [0..l-1];
    aset1 := Filtered(Ald, a -> (a mod p = 2) or (a mod p = p-2));

    type1u := List(aset1, a -> [[uinv, u, uinv, a*u], 1]);
    type11 := List(Ald, a -> [[a], 1]);

    Type1 := Concatenation(type1u, type11);

    # Type 2. For there is a type of conj class for each 1 <= h <= ld.
    type2h := function(h)
        local ans, dset, e;
        if h = ld then
            return [ [[0, 0], 1], [[0, 0], 0] ];
        else
            ans := [];
            for e in [1, -1] do
                dset := Filtered([0..p^(ld-h)-1], x -> x mod p = 0);
                ans := Concatenation(ans, List(dset, d -> [[-e*p^h*d*uinv, -e*p^h*u],  (1+e)/2]));
                ans := Concatenation(ans, List([0..p^(ld-h)-1], d -> [[-e*p^h*d, -e*p^h], (1+e)/2]));
            od;
            return ans;
        fi;
    end;

    # Type2 is the collection of all type2h's. Among them, type2h(ld) contains id and -id.
    Type2 := type2h(ld);
    for i in [1..ld-1] do
        Type2 := Concatenation(Type2, type2h(i));
    od;

    # Finally, combine Type1 and Type2 and return. Put Type2 in the front to make sure id and -id are the first two conj classes.
    return [Concatenation(Type2, Type1), uinv, u];
end );

#----------------------------------------
# Returns conj class representative information of SL(2, Z/2^ld Z) for ld (>= 1).
# Each conj class representative is represented by a list of length 2. To rewrite the conj class rep as products of s and t (NW notation) corresponding to a, we can compute
# Prod(a[1], n -> t^n*s) * s^(2*a[2]).
# Sorted the output so that id and -id are the first two conj classes.
#---------------------------------------
InstallGlobalFunction( _SL2ConjClassesEven,
function(ld)
    local l, 3inv, 5inv, 7inv, Type1, Type2, d, h, e1, e1inv, e2, e2inv, a;

    l := 2^ld;

    # Find the multiplicative inverse of 3, 5 and 7.
    3inv := 3^-1 mod l;
    5inv := 5^-1 mod l;
    7inv := 7^-1 mod l;

    # There are 2 types of conj classes.
    # Type 1.
    Type1 := List([0..l-1], a -> [[a], 1]);

    for a in List([0..l/2-1], x -> 2*x) do
        if ld > 1 then
            Add(Type1, [[3inv, 3, 3inv, 3*a], 1]);
        fi;

        if (ld > 2) and (a mod 4 = 2) then
            Add(Type1, [[5inv, 5, 5inv, 5*a], 1]);
            Add(Type1, [[7inv, 7, 7inv, 7*a], 1]);
        fi;
    od;

    # Type 2 are labelled by 1 <= h <= ld.
    # When h = ld, the corresponding conj classes are elements in the center. id and -id are put in the front.
    if ld = 1 then
        Type2 := [ [[0, 0], 1] ];
    elif ld = 2 then
        Type2 := [ [[0, 0], 1], [[0, 0], 0] ];
    else
        Type2 := [ [[0,0], 1], [[0, 0], 0], [[1+2^(ld-1), 1+2^(ld-1), 1+2^(ld-1)], 0], [[-(1+2^(ld-1)), -(1+2^(ld-1)), -(1+2^(ld-1))], 0] ];
    fi;

    # h = 1.
    if ld = 2 then
        Type2 := Concatenation(Type2, List([0..1], d -> [[-2*d mod l, -2 mod l], 1]));
    elif ld = 3 then
        Type2 := Concatenation(Type2, List([0..3], d -> [[-2*d mod l, -2 mod l], 1]));
        Type2 := Concatenation(Type2, List([0..1], d -> [[-2*d*3inv mod l, -2*3 mod l], 1]));
    elif ld > 3 then
        for d in [0..2^(ld-1)-1] do
            Add(Type2, [[-2*d mod l, -2 mod l], 1]);
            if (d mod 8 = 0) or (d mod 8 = 1) then
                Add(Type2, [[-2*d*3inv mod l, -2*3 mod l], 1]);
                Add(Type2, [[-2*d*5inv mod l, -2*5 mod l], 1]);
                Add(Type2, [[-2*d*7inv mod l, -2*7 mod l], 1]);
            elif (d mod 8 = 4) or (d mod 8 = 5) then
                Add(Type2, [[-2*d*3inv mod l, -2*3 mod l], 1]);
            else
                Add(Type2, [[-2*d*5inv mod l, -2*5 mod l], 1]);
            fi;
        od;
    fi;

    # ld >= 3 and h = ld - 1.
    if ld > 2 then
        for d in [0..1] do
            Add(Type2, [[-2^(ld-1)*d mod l, -2^(ld-1) mod l], 1]);
            Add(Type2, [[2^(ld-1)*d mod l, 2^(ld-1) mod l], 0]);
        od;
    fi;

    # ld >= 4 and h = ld - 2.
    if ld > 3 then
        for d in [0..3] do
            Add(Type2, [[-2^(ld-2)*d mod l, -2^(ld-2) mod l], 1]);
            Add(Type2, [[2^(ld-2)*d mod l, 2^(ld-2) mod l], 0]);
            if d = 0 or d = 1 then
                Add(Type2, [[-2^(ld-2)*d*3inv mod l, -2^(ld-2)*3 mod l], 1]);
                Add(Type2, [[2^(ld-2)*d*3inv mod l, 2^(ld-2)*3 mod l], 0]);
            fi;
        od;
    fi;

    # ld >= 5, 2 <= h <= ld - 3.
    if ld > 4 then
        for h in [2..ld-3] do
            for d in [0..2^(ld-h)-1] do
                Add(Type2, [[-2^h*d mod l, -2^h mod l], 1]);
                Add(Type2, [[2^h*d mod l, 2^h mod l], 0]);
                if (d mod 4 = 0) or (d mod 4 = 1) then
                    Add(Type2, [[-2^h*d*3inv mod l, -2^h*3 mod l], 1]);
                    Add(Type2, [[2^h*d*3inv mod l, 2^h*3 mod l], 0]);
                fi;
                if (d mod 8 = 0) or (d mod 8 = 2) or (d mod 8 = 6) then
                    Add(Type2, [[-2^h*d*5inv mod l, -2^h*5 mod l], 1]);
                    Add(Type2, [[2^h*d*5inv mod l, 2^h*5 mod l], 0]);
                fi;
                if (d mod 8 = 0) then
                    Add(Type2, [[-2^h*d*7inv mod l, -2^h*7 mod l], 1]);
                    Add(Type2, [[2^h*d*7inv mod l, 2^h*7 mod l], 0]);
                fi;
            od;
        od;
    fi;

    # ld > 3, 3 <= h <= ld - 1.
    if ld > 3 then
        for h in [3..ld-1] do
            e1 := 1 + 2^(h-1);
            e1inv := e1^(-1) mod l;
            e2 := -e1 mod l;
            e2inv := -e1inv mod l;
            for d in [0..2^(ld-h)-1] do
                Add(Type2, [[-2^h*d*e1inv mod l, e1inv, e1, e1inv, -2^h*e1inv mod l], 1]);
                Add(Type2, [[-2^h*d*e2inv mod l, e2inv, e2, e2inv, -2^h*e2inv mod l], 1]);
            od;
        od;
    fi;

    # Finally, combine Type1 and Type2 and return. Put Type2 in the front to make sure id and -id are the first two conj classes.
    return [Concatenation(Type2, Type1), 3inv, 5inv, 7inv];
end );

#----------------------------------------
# Returns conj class representative information of SL(2, Z/p^ld Z) for ld (>= 1) using ConjClassesOdd or ConjClassesEven.
#---------------------------------------
f := MemoizePosIntFunction(function(l)
    local v;
    v := PrimePowersInt(l);
    if not Length(v) = 2 then
        Error("level must be a prime power.");
    else
        if v[1] = 2 then
            return _SL2ConjClassesEven(v[2]);
        elif v[1] > 2 then
            return _SL2ConjClassesOdd(v[1], v[2]);
        fi;
    fi;
end);

InstallGlobalFunction( _SL2ConjClasses, f );

#-----------------------------------------------
# Irreducibility check
# Input: p, ld, S, and T for a rep of SL(2,Z/p^ldZ), given by RepN, RepD and RepRs. Returns the inner product of the character of this representation with itself. If needed, we can also list the character values on conjugacy classes.
#-----------------------------------------------
# InstallGlobalFunction( _SL2Reps_CharNorm,
# function(S, T)
#     local s2, chi, Csize, G, n, p, ld, CC, TS, Du, c, i, primepower;
#     # Size of SL(2, Z/p^ldZ).
#     n := Order(T);
#     primepower := Factors(n);
#     if Length(AsSet(primepower)) > 1 then
#         Error("The level is not a prime power. We cannot proceed!");
#     fi;
#     p := primepower[1];
#     ld := Length(primepower);
#     G := p^(3*ld) - p^(3*ld-2);

#     CC := _SL2Reps_ConjClassesOdd(p, ld);
#     s2 := S^2;
#          s2 := s2[1][1];
#          TS:=[S];
#     for i in [1..n-1] do
#         Add(TS, T*TS[i]);
#     od;
#          Du := TS[CC[2]+1]*TS[CC[3]+1]*TS[CC[2]+1];

#     chi := []; Csize := [];
#     for c in CC[1] do
#         Add(Csize, c[3]);
#         if Length(c[1])=4 and c[1][1]=1 then
#             Add(chi, Trace(TS[(c[1][4] mod n) + 1] * s2^(c[2])));
#         elif Length(c[1])=4 and c[1][1]<>1 then
#             Add(chi, Trace(Du*TS[(c[1][4] mod n) + 1] * s2^(c[2])));
#         elif Length(c[1])=2 then
#             Add(chi, Trace(Product(c[1], m -> TS[(m mod n) + 1]) * s2^(c[2])));
#                 fi;
#     od;

#     return ComplexConjugate(chi)*ListN(chi, Csize, \*)/G;
# end );

#-----------------------------------------------
# TODO
#-----------------------------------------------
# InstallGlobalFunction( _SL2Reps_RepChi,
# function(S, T)
#     local id, TS, Du, Dh, h, Dhmap, DS, DR, e1, e1inv, e2, e2inv, CC, ccl, s2, n, G, s, t,  C1, o, chi, i, j, jinv, c, pos, primepower, p, ld, ProdTrace;

#     # S,T are normalized S and T matrices
#     n := Order(T);
#     primepower := Factors(n);
#     if Length(AsSet(primepower))>1 then
#         Error("The level is not a prime power. We cannot proceed!");
#     fi;
#     p := primepower[1]; ld := Length(primepower);
#     CC := _SL2Reps_ConjClasses(p, ld);
#     o := ZmodnZObj(1,n);
#     s := [[0,1],[-1,0]] * o;
#     t := [[1,1],[0,1]] * o;
#     id := s^0;
#     G := Group([s,t]);

#     s2 := S^2;
#     s2 := s2[1][1];
#     TS := [S];
#     for i in [1..n-1] do
#         Add(TS, T*TS[i]);
#     od;

#     ProdTrace := function(A, B)
#         # assumes they are both square matrices of same size
#         return Sum([1..Length(A)], x -> Sum([1..Length(A)], y -> A[x][y]*B[y][x]));
#     end;

#     C1 := [];
#     if p > 2 then
#         Du := TS[CC[2]+1]*TS[CC[3]+1]*TS[CC[2]+1];
#         for c in CC[1] do
#             if Length(c[1])=1  then
#                 Add(C1, [t^(c[1][1]) * s * (-1)^(c[2]), Trace(TS[(c[1][1] mod n) + 1] * s2^(c[2]))]);
#             elif Length(c[1])=4 and c[1][1]<>1 then
#                 Add(C1, [Product(c[1], m -> t^m * s) * (-1)^(c[2]), ProdTrace(Du,TS[(c[1][4] mod n) + 1]) * s2^(c[2])]);
#             elif Length(c[1])=2 then
#                 Add(C1, [Product(c[1], m -> t^m * s) * (-1)^(c[2]), ProdTrace(TS[(c[1][1] mod n)+1], TS[(c[1][2] mod n) + 1]) * s2^(c[2])]);
#             fi;
#         od;
#     fi;

#     if p = 2 then
#         Dh := [];
#         for h in [3..ld] do
#             e1 := 1 + 2^(h-1);
#             e1inv := e1^-1 mod n;
#             Add(Dh, [e1inv,TS[e1inv + 1] * TS[e1 + 1] * TS[e1inv + 1]]);
#                     e2 := (-e1) mod n;
#             e2inv := (-e1inv) mod n;
#             Add(Dh, [e2inv,TS[e2inv + 1] * TS[e2 + 1] * TS[e2inv + 1]]);
#         od;
#         for i in [3,5,7] do
#             j := i mod n;
#             jinv := j^-1 mod n;
#             Add(Dh, [jinv,TS[jinv + 1] * TS[j + 1] * TS[jinv + 1]]);
#         od;
#         Dh := AsSet(Dh);

#         Dhmap := function(uinv)
#             local x;
#             x := First(Dh, y -> y[1] = uinv);
#             return x[2];
#         end;

#         for c in CC[1] do
#             if Length(c[1]) = 4 and c[1][1] <> 1 then
#                 Add(C1, [Product(c[1], m -> t^m * s) * (-1)^(c[2]), ProdTrace(Dhmap(c[1][1]), TS[(c[1][4] mod n) + 1]) * s2^(c[2])]);
#             elif Length(c[1]) = 2 then
#                 Add(C1, [Product(c[1], m -> t^m*s)*(-1)^(c[2]), ProdTrace(TS[(c[1][1] mod n)+1], TS[(c[1][2] mod n)+1])*s2^(c[2])]);
#             elif Length(c[1])=3 then
#                 Add(C1, [Product(c[1], m -> t^m*s)*(-1)^(c[2]), Trace(Dhmap(c[1][1] mod n))*s2^(c[2])]);
#             elif Length(c[1])=5 then
#                 Add(C1, [Product(c[1], m -> t^m*s)*(-1)^(c[2]), ProdTrace((TS[(c[1][1] mod n)+1]*Dhmap(c[1][2] mod n)), TS[(c[1][5] mod n)+1])*s2^(c[2])]);
#             elif Length(c[1])=1 then
#                 Add(C1, [Product(c[1], m -> t^m*s)*(-1)^(c[2]), Trace(TS[(c[1][1] mod n)+1])*s2^(c[2])]);
#             fi;
#         od;
#     fi;

#     chi := List(C1, x -> x[2]);
#     ccl := List(C1, x -> ConjugacyClass( G, x[1] ) );
#     if Sum(ccl, x -> Size(x)) <> Size(G) then
#         Error("Wrong conjugacy classes.");
#     fi;

#     SetConjugacyClasses( G, ccl );

#     return Character(G, chi);
# end );

#-----------------------------------------------
# Let G be SL(2, Z/nZ),  C a list of Conjugacy Classes of G,
# and CC a list of representatives of Conjugacy classes of G
# The function return a list pos such that C[i]=CC[pos[i]]^G.
#-----------------------------------------------
# InstallGlobalFunction( _SL2Reps_ClassMap,
# function(C, CC)

#     local  rep, c, C1, u, v, traces, tracelist, traceu, vlist, i, pos, PowTrace, urep, perm;


#     PowTrace:=function(x)
#         local j, temp, power;
#         temp:=[]; power:=x;
#         for j in [1..Order(x)-1] do
#             Add(temp, Trace(power));
#             power:=power*x;
#         od;
#         return temp;
#     end;

#     C1:=List(CC, c->[c, PowTrace(c)]);
#     traces:=Set(C1, c->c[2]);
#     tracelist:=Set(traces, x->[x, Filtered(C1, c->c[2]=x)]);
#     perm:=[];

#     for u in C do
#         urep:=Representative(u);
#         traceu:=PowTrace(urep);
#         vlist:=First(tracelist, x->traceu=x[1]);
#         if Length(vlist[2])=1 then
#             Add(perm, Position(CC, vlist[2][1][1]));
#             RemoveSet(tracelist, vlist);
#         elif Length(vlist[2])>1 then
#             v:=First(vlist[2], x->x[1] in u);
#             Add(perm, Position(CC, v[1]));
#             pos:=Position(vlist[2], v);
#             Remove(vlist[2], pos);
#         fi;
#     od;

#     return perm;
# end );

#-----------------------------------------------
# Find the position of the given representation among the irreducibles of its underlying group.
#-----------------------------------------------
# InstallGlobalFunction( _SL2Reps_ChiTest,
# function(chi)
#     local G;
#     G := UnderlyingGroup(chi);
#     return Position(Irr(G), chi);
# end );
