#
# SL2Reps: Constructs representations of SL2(Z).
#
# Representations of type R.

#! @Chapter Methods for constructing representations

#! @Section Representations of type R
#!
#! TODO: Description of type $R$.

#! @Arguments p,ld,sigma,r,t
#! @Returns a record describing $(M,Q)$
#! @Description
#!  Constructs information about the underlying quadratic module $(M,Q)$ of type $R$.
DeclareGlobalFunction( "SL2Reps_ModuleR" );

#! @Arguments p,ld,sigma,r,t,chi_index
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the irreducible representation(s) of type $R$ with level $p^\ld$
#!  corresponding to the character $\chi$ indexed by `chi_index`.
DeclareGlobalFunction( "SL2Reps_RepR" );

#! @Arguments p,ld,r
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the irreducible representation(s) of unary type $R$ (that is, with $\sigma = \ld$)
#!  with level $p^\ld$.
DeclareGlobalFunction( "SL2Reps_RepRUnary" );
