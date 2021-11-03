#
# SL2Reps: Constructs representations of SL2(Z).
#
# Representations of type D.

#! @Chapter Irreps
#! @ChapterTitle Irreducible representations of prime-power level

#! @Section Representations of type D
#!
#! See <Ref Sect="Chapter_Description_Section_Weyl_Subsection_Type_D"/>.

#! @Arguments p,ld
#! @Returns a record `rec(Agrp, Bp, Char, IsPrim)` describing $(M,Q)$
#! @Description
#!  Constructs information about the underlying quadratic module $(M,Q)$ of type $D$.
#!
#!  `Agrp` describes the elements of $\mathfrak{A} = (\mathbb{Z}/p^\lambda\mathbb{Z})^\times$
#!  (see <Cite Key="NW76" Where="Section 2.1"/>).
#!
#!  `Bp` describes a set of representatives for the $\mathfrak{A}$-orbits on $M^\times$, which
#!  correspond to a basis the $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$-invariant subspace
#!  associated to any primitive character $\chi \in \hat{\mathfrak{A}}$ with $\chi^2 \not\equiv 1$.
#!  For other characters, we must use different bases which are particular to each case.
#!
#!  `Char(i,j)` converts the `chi_index` used in <Ref Func="SL2Reps_RepD"/> to a function.
#!
#!  `IsPrim(chi)` tests whether a given character (e.g. from `Char(i,j)`) is primitive.
DeclareGlobalFunction( "SL2Reps_ModuleD" );

#! @Arguments p,ld,chi_index
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the modular data for the irreducible representation(s) of type $D$ with
#!  level $p^\lambda$ corresponding to the character $\chi$ indexed by `chi_index`.
DeclareGlobalFunction( "SL2Reps_RepD" );
