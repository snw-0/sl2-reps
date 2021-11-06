#
# SL2Reps: Constructs representations of SL2(Z).
#
# Representations of type R.

#! @Chapter Irreps

#! @Section Representations of type R
#!
#! See <Ref Sect="Chapter_Description_Section_Weil_Subsection_Type_R"/>.

#! @Arguments p,ld,sigma,r,t
#! @Returns a record `rec(Agrp, Char, IsPrim, Nm, Ord, Prod, c, tM)` describing $(M,Q)$
#! @Description
#!  Constructs information about the underlying quadratic module $(M,Q)$ of type $R$.
#!
#!  `Agrp` describes the elements of $\mathfrak{A} = \{\varepsilon \in M^\times \mid \operatorname{Nm}(\varepsilon) = 1 \}$
#!  (see <Cite Key="NW76" Where="Section 2.3 - 2.5"/>).
#!
#!  Representatives for the $\mathfrak{A}$-orbits on $M^\times$ can depend on the choice of
#!  character, even for primitive characters $\chi$ with $\chi^2 \not\equiv 1$.  Thus, we
#!  cannot provide them here, and they are instead calculated by <Ref Func="SL2IrrepR"/>.
#!
#!  `Char(i,j)` converts the `chi_index` used in <Ref Func="SL2IrrepR"/> to a function.
#!
#!  `IsPrim(chi)` tests whether a given character (e.g. from `Char`) is primitive.
#!
#!  `Nm(a)`, `Ord(a)`, and `Prod(a,b)` are the norm, order, and product functions on $M$, respectively.
#!
#!  `c` is a scalar used in calculating the $S$-matrix; namely $c = \frac{1}{|M|} \sum_{x \in M} \mathbf{e}(Q(x))$.
#!
#!  `tM` is the group $M - pM$.
DeclareGlobalFunction( "SL2ModuleR" );

#! @Arguments p,ld,sigma,r,t,chi_index
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the modular data for the irreducible representation(s) of type $R$ with level $p^\lambda$
#!  corresponding to the character $\chi$ indexed by `chi_index`.
#!
#!  When $\sigma = \lambda$, this falls through to <Ref Func="SL2IrrepRUnary"/>.
DeclareGlobalFunction( "SL2IrrepR" );

#! @Arguments p,ld,r
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the modular data for the irreducible representation(s) of unary type $R$ (that is, with $\sigma = \lambda$)
#!  with level $p^\lambda$.
DeclareGlobalFunction( "SL2IrrepRUnary" );
