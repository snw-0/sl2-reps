#
# SL2Reps: Constructs representations of SL2(Z).
#
# Representations of type R.

#! @Chapter Methods for Constructing Representations

#! @Section Representations of Type R
#!
#! TODO: Description of type $R$.

#! @Arguments p,lambda,sigma,r,t
#! @Returns a record describing $M$
#! @Description
#!  Constructs information about the underlying module $M$ of type $R$.
DeclareGlobalFunction( "SL2Reps_ModuleR" );

#! @Arguments p,lambda,sigma,r,t,chi_index
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the irreducible representation(s) of type $R$ with level $p^\lambda$
#!  corresponding the character $\chi$ indexed by `chi_index`.
DeclareGlobalFunction( "SL2Reps_RepR" );

#! @Arguments p,lambda,r
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the irreducible representation(s) of unary type $R$ (that is, with $\sigma = \lambda$)
#!  with level $p^\lambda$.
DeclareGlobalFunction( "SL2Reps_RepRUnary" );
