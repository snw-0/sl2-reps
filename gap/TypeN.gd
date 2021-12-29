#
# SL2Reps: Constructing symmetric representations of SL(2,Z).
#
# Representations of type N.
#


#! @Chapter Irreps

#! @Section Representations of type N
#!
#! See Section <Ref Sect="Chapter_Description_Section_Weil_Subsection_Type_N"/>.


#! @Arguments p,lambda
#! @Returns a record `rec(Agrp, Bp, Char, IsPrim, Nm, Prod)` describing $(M,Q)$.
#! @Description
#!  Constructs information about the underlying quadratic module $(M,Q)$ of type $N$, for
#!  $p$ a prime and $\lambda \geq 1$.
#!
#!  `Agrp` is a list describing the elements of $\mathfrak{A}$.
#!  Each element $a \in \mathfrak{A}$ is represented in `Agrp` by a list `[v, a]`,
#!  where `v` is a list defined by $a = \alpha^{\mathtt{v[1]}} \beta^{\mathtt{v[2]}}$.
#!  Note that $\alpha$ is trivial, and hence `v[1]` is irrelevant, when $\mathfrak{A}$ is cyclic.
#!
#!  `Bp` is a list of representatives for the $\mathfrak{A}$-orbits on $M^\times$, which
#!  correspond to a basis for the $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$-invariant subspace
#!  associated to any primitive character $\chi \in \widehat{\mathfrak{A}}$ with $\chi^2 \not\equiv 1$.
#!  This is the basis given by <Cite Key="NW76"/>, which may result in a non-symmetric representation;
#!  if this occurs, we perform a change of basis in <Ref Func="SL2IrrepD"/> to obtain a symmetric
#!  representation.
#!  For non-primitive characters, we must use different bases which are particular to each case.
#!
#!  `Char(i,j)` converts two integers $i$, $j$ to a function representing the character $\chi_{i,j} \in \widehat{\mathfrak{A}}$.
#!
#!  `IsPrim(chi)` tests whether the output of `Char(i,j)` represents a primitive character.
#!
#!  `Nm(a)` and `Prod(a,b)` are the norm and product functions on $M$, respectively.
DeclareGlobalFunction( "SL2ModuleN" );

#! @Arguments p,lambda,chi_index
#! @Returns a list of lists of the form $[S,T]$.
#! @Description
#!  Constructs the modular data for the irreducible representation(s) of type $N$ with
#!  level $p^\lambda$, for $p$ a prime and $\lambda \geq 1$, corresponding to the
#!  character $\chi$ indexed by `chi_index = [i,j]`
#!  (see the discussion of `Char(i,j)` in <Ref Func="SL2ModuleN"/>).
#!
#!  Here $S$ is symmetric and unitary and $T$ is diagonal.
#!
#!  Depending on the parameters, $W(M,Q)$ will contain either 1 or 2 such irreps.
DeclareGlobalFunction( "SL2IrrepN" );
