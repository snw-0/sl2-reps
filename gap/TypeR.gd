#
# SL2Reps: Constructing symmetric representations of SL(2,Z).
#
# Representations of type R.
#


#! @Chapter Irreps

#! @Section Representations of type R
#!
#! See Section <Ref Sect="Chapter_Description_Section_Weil_Subsection_Type_R"/>.


#! @Arguments p,lambda,sigma,r,t
#! @Returns a record `rec(Agrp, Bp, Char, IsPrim, Nm, Ord, Prod, c, tM)` describing $(M,Q)$.
#! @Description
#!  Constructs information about the underlying quadratic module $(M,Q)$ of type $R$, for
#!  $p$ a prime. The additional parameters $\lambda$, $\sigma$, $r$, and $t$ should
#!  be integers chosen as follows.
#!
#!  If $p$ is an odd prime, let $\lambda \geq 2$, $\sigma \in \{1, \dots, \lambda - 1\}$,
#!  and $r,t \in \{1,u\}$ with $u$ a quadratic non-residue mod $p$.  Note that $\sigma = \lambda$
#!  is a valid choice for type $R$, however, this gives the unary case, and so is not handled by this
#!  function, as it is decomposed in a different way; for this case, use <Ref Func="SL2IrrepRUnary"/> instead.
#!
#!  If $p=2$, let $\lambda \geq 2$, $\sigma \in \{0, \dots, \lambda-2\}$ and $r,t \in \{1,3,5,7\}$.
#!
#!  `Agrp` is a list describing the elements of $\mathfrak{A}$. Each element $a$ of
#!  $\mathfrak{A}$ is represented in `Agrp` by a list `[v, a]`,
#!  where `v` is a list defined by $a = \alpha^{\mathtt{v[1]}} \beta^{\mathtt{v[2]}}$.
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
#!  `Nm(a)`, `Ord(a)`, and `Prod(a,b)` are the norm, order, and product functions on $M$, respectively.
#!
#!  `c` is a scalar used in calculating the $S$-matrix; namely
#!  $c = \frac{1}{|M|} \sum_{x \in M} \mathbf{e}(Q(x))$.
#!  Note that this is equal to $S_Q(-1) / \sqrt{|M|}$, where
#!  $S_Q(-1)$ is the central charge (see Section <Ref Sect="Chapter_Description_Section_Construction_Subsection_Weil_representations"/>).
#!
#!  `tM` is a list describing the elements of the group $M - pM$.
DeclareGlobalFunction( "SL2ModuleR" );

#! @Arguments p,lambda,sigma,r,t,chi_index
#! @Returns a list of lists of the form $[S,T]$.
#! @Description
#!  Constructs the modular data for the irreducible representation(s) of type $R$ with
#!  parameters $p$, $\lambda$, $\sigma$, $r$, and $t$,
#!  corresponding to the character $\chi$ indexed by `chi_index = [i,j]`
#!  (see the discussions of $\sigma$, $r$, $t$, and `Char(i,j)` in <Ref Func="SL2ModuleR"/>).
#!
#!  Here $S$ is symmetric and unitary and $T$ is diagonal.
#!
#!  Depending on the parameters, $W(M,Q)$ will contain either 1 or 2 such irreps.
#!
#!  If $\sigma = \lambda$ for $p \neq 2$, then the second factor of $M$ is trivial
#!  (and hence $t$ is irrelevant), so this falls through to <Ref Func="SL2IrrepRUnary"/>.
DeclareGlobalFunction( "SL2IrrepR" );

#! @Arguments p,lambda,r
#! @Returns a list of lists of the form $[S,T]$.
#! @Description
#!  Constructs the modular data for the irreducible representation(s) of unary type $R$
#!  (that is, the special case where $\sigma = \lambda$) with $p$ an odd prime,
#!  $\lambda$ a positive integer, and $r \in \{1,u\}$ with $u$ a quadratic non-residue mod $p$.
#!
#!  Here $S$ is symmetric and unitary and $T$ is diagonal.
#!
#!  In this case, $W(M,Q)$ always contains exactly 2 such irreps.
DeclareGlobalFunction( "SL2IrrepRUnary" );
