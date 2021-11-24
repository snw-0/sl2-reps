#
# SL2Reps: Constructs representations of SL2(Z).
#
# Representations of type D.


#! @Chapter Irreps
#! @ChapterTitle Irreducible representations of prime-power level
#!  Methods for generating individual irreducible representations of
#!  $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$ for a given level $p^\lambda$.
#!
#!  In each case (except the unary type $R$, for which see <Ref Func="SL2IrrepRUnary"/>),
#!  the underlying module $M$ is of rank $2$, so its elements have the form $(x,y)$
#!  and are thus represented by lists $[x,y]$.
#!
#!  Characters of the abelian group $\mathfrak{A} = \langle\alpha\rangle \times \langle\beta\rangle$,
#!  have the form $\chi_{i,j}$, given by
#!  \[\chi_{i,j}(\alpha^{v}\zeta^{w}) \mapsto \mathbf{e}\left(\frac{vi}{|\alpha|}\right) \mathbf{e}\left(\frac{wj}{|\zeta|}\right)~,\]
#!  where $i$ and $j$ are integers.  We therefore represent each character by a list $[i,j]$.
#!  Note that $\mathfrak{A}$ may be cyclic, in which case $\alpha$ or $\beta$ will be trivial.

#! @Section Representations of type D
#!
#! See Section <Ref Sect="Chapter_Description_Section_Weil_Subsection_Type_D"/>.


#! @Arguments p,ld
#! @Returns a record `rec(Agrp, Bp, Char, IsPrim)` describing $(M,Q)$
#! @Description
#!  Constructs information about the underlying quadratic module $(M,Q)$ of type $D$, for
#!  $p$ a prime and $\lambda \geq 1$.
#!
#!  `Agrp` is a list describing the elements of $\mathfrak{A}$.
#!  Each element $a \in \mathfrak{A}$ is represented in `Agrp` by a list `[v, a, a_inv]`,
#!  where `v` is a list defined by $a = \alpha^{\mathtt{v[1]}} \zeta^{\mathtt{v[2]}}$.
#!  Note that $\zeta$ is trivial, and hence `v[2]` is irrelevant, when $\mathfrak{A}$ is cyclic.
#!
#!  `Bp` is a list of representatives for the $\mathfrak{A}$-orbits on $M^\times$, which
#!  correspond to a basis for the $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$-invariant subspace
#!  associated to any primitive character $\chi \in \widehat{\mathfrak{A}}$ with $\chi^2 \not\equiv 1$.
#!  For other characters, we must use different bases which are particular to each case.
#!
#!  `Char(i,j)` converts two integers $i$, $j$ to a function representing a character of $\mathfrak{A}$.
#!  Each character in $\hat{\mathfrak{A}}$ is of the form $\chi_{i,j}$, given by
#!  \[\chi_{i,j}(\alpha^{v}\zeta^{w}) \mapsto \mathbf{e}\left(\frac{vi}{|\alpha|}\right) \mathbf{e}\left(\frac{wj}{|\zeta|}\right)~.\]
#!  Note that $j$ is irrelevant when $\mathfrak{A}$ is cyclic.
#!
#!  `IsPrim(chi)` tests whether the output of `Char(i,j)` represents a primitive character.
DeclareGlobalFunction( "SL2ModuleD" );

#! @Arguments p,ld,chi_index
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the modular data for the irreducible representation(s) of type $D$ with
#!  level $p^\lambda$, for $p$ a prime and $\lambda \geq 1$, corresponding to the
#!  character $\chi$ indexed by `chi_index = [i,j]`
#!  (see the discussion of `Char(i,j)` in <Ref Func="SL2ModuleD"/>).
#!
#!  Depending on the parameters, $W(M,Q)$ will contain either 1 or 2 such irreps.
DeclareGlobalFunction( "SL2IrrepD" );
