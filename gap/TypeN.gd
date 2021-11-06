#
# SL2Reps: Constructs representations of SL2(Z).
#
# Representations of type N.

#! @Chapter Irreps

#! @Section Representations of type N
#!
#! See <Ref Sect="Chapter_Description_Section_Weyl_Subsection_Type_N"/>.

#! @Arguments p,ld
#! @Returns a record `rec(Agrp, Bp, Char, Nm, Prod)` describing $(M,Q)$
#! @Description
#!  Constructs information about the underlying quadratic module $(M,Q)$ of type $N$.
#!
#!  `Agrp` describes the elements of $\mathfrak{A} = \{\varepsilon \in M^\times \mid \operatorname{Nm}(\varepsilon) = 1 \}$
#!  (see <Cite Key="NW76" Where="Section 2.2"/>).
#!
#!  `Bp` describes a set of representatives for the $\mathfrak{A}$-orbits on $M^\times$, which
#!  correspond to a basis the $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$-invariant subspace
#!  associated to any primitive character $\chi \in \hat{\mathfrak{A}}$ with $\chi^2 \not\equiv 1$.
#!  For other characters, we must use different bases which are particular to each case.
#!
#!  `Char(i,j)` converts the `chi_index` used in <Ref Func="SL2IrrepN"/> to a function.
#!
#!  `Nm(a)` and `Prod(a,b)` are the norm and product functions on $M$, respectively.
DeclareGlobalFunction( "SL2ModuleN" );

#! @Arguments p,ld,chi_index
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the modular data for the irreducible representation(s) of type $N$ with level $p^\lambda$
#!  corresponding to the character $\chi$ indexed by `chi_index`.
DeclareGlobalFunction( "SL2IrrepN" );
