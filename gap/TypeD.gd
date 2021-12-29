#
# SL2Reps: Constructing symmetric representations of SL(2,Z).
#
# Representations of type D.
#


#! @Chapter Irreps
#! @ChapterTitle Irreducible representations of prime-power level
#!  Methods for generating individual irreducible representations of
#!  $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$ for a given level $p^\lambda$.
#!
#!  After generating a representation $\rho$ by means of the bases in <Cite Key="NW76"/>, we perform
#!  a change of basis that results in a symmetric representation equivalent to $\rho$.
#!
#!  In each case (except the unary type $R$, for which see <Ref Func="SL2IrrepRUnary"/>),
#!  the underlying module $M$ is of rank $2$, so its elements have the form $(x,y)$
#!  and are thus represented by lists `[x,y]`.
#!
#!  Characters of the abelian group $\mathfrak{A} = \langle\alpha\rangle \times \langle\beta\rangle$
#!  have the form $\chi_{i,j}$, given by
#!  $$\chi_{i,j}(\alpha^{v}\beta^{w}) \mapsto \mathbf{e}\left(\frac{vi}{|\alpha|}\right) \mathbf{e}\left(\frac{wj}{|\beta|}\right)~,$$
#!  where $i$ and $j$ are integers.  We therefore represent each character by a list `[i,j]`.
#!  Note that in some cases $\alpha$ or $\beta$ is trivial, and the corresponding index
#!  $i$ or $j$ is therefore irrelevant.
#!
#!  We write `p=`$p$, `lambda=`$\lambda$, `sigma=`$\sigma$, and `chi=`$\chi$.

#! @Section Representations of type D
#!
#! See Section <Ref Sect="Chapter_Description_Section_Weil_Subsection_Type_D"/>.


#! @Arguments p,lambda
#! @Returns a record `rec(Agrp, Bp, Char, IsPrim)` describing $(M,Q)$.
#! @Description
#!  Constructs information about the underlying quadratic module $(M,Q)$ of type $D$, for
#!  $p$ a prime and $\lambda \geq 1$.
#!
#!  `Agrp` is a list describing the elements of $\mathfrak{A}$.
#!  Each element $a \in \mathfrak{A}$ is represented in `Agrp` by a list `[v, a, a_inv]`,
#!  where `v` is a list defined by $a = \alpha^{\mathtt{v[1]}} \beta^{\mathtt{v[2]}}$.
#!  Note that $\beta$ is trivial, and hence `v[2]` is irrelevant, when $\mathfrak{A}$ is cyclic.
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
DeclareGlobalFunction( "SL2ModuleD" );

#! @Arguments p,lambda,chi_index
#! @Returns a list of lists of the form $[S,T]$.
#! @Description
#!  Constructs the modular data for the irreducible representation(s) of type $D$ with
#!  level $p^\lambda$, for $p$ a prime and $\lambda \geq 1$, corresponding to the
#!  character $\chi$ indexed by `chi_index = [i,j]`
#!  (see the discussion of `Char(i,j)` in <Ref Func="SL2ModuleD"/>).
#!
#!  Here $S$ is symmetric and unitary and $T$ is diagonal.
#!
#!  Depending on the parameters, $W(M,Q)$ will contain either 1 or 2 such irreps.
DeclareGlobalFunction( "SL2IrrepD" );
