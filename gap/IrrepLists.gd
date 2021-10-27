#
# SL2Reps: Constructs representations of SL2(Z).
#

# Lists of Representations.

#! @Chapter Methods for constructing lists of representations

#! @Section Irreducible representations
#!
#! TODO: Description of this.

#! @Arguments degree
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  that are exactly the given degree and have prime power level.
DeclareGlobalFunction( "SL2Reps_PrimePowerIrrepsOfDegree" );

#! @Arguments max_degree
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  that are at most the given degree and have prime power level.
DeclareGlobalFunction( "SL2Reps_PrimePowerIrrepsOfDegreeAtMost" );

#! @Arguments degree
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  that are exactly the given degree.
DeclareGlobalFunction( "SL2Reps_IrrepsOfDegree" );

#! @Arguments degree
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  that are at most the given degree.
DeclareGlobalFunction( "SL2Reps_IrrepsOfDegreeAtMost" );

#! @Arguments p, lambda
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  with level exactly $p^\lambda$, or equivalently, those which factor through
#!  $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$ but not through
#!  $\mathrm{SL}_2(\mathbb{Z}/m\mathbb{Z})$ for any $m < p^\lambda$.
DeclareGlobalFunction( "SL2Reps_PrimePowerIrrepsOfLevel" );

#! @Returns a list of records of the form `rec(S, T, degree, level, name)`
#! @Description
#!  Constructs a list of the 18 exceptional irreps of $\mathrm{SL}_2(\mathbb{Z})$.
#!  All but 18 of the irreps of $\mathrm{SL}_2(\mathbb{Z})$ with prime power level
#!  may be found as subrepresentations of quadratic modules; this function constructs the
#!  18 exceptions, all of which may be found as the tensor product of two
#!  subrepresentations of quadratic modules.
DeclareGlobalFunction( "SL2Reps_ExceptionalPrimePowerIrreps" );
