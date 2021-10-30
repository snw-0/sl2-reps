#
# SL2Reps: Constructs representations of SL2(Z).
#

# Lists of Representations.

#! @Chapter Lists
#! @ChapterTitle Lists of representations

#! @Section Degree
#! @SectionTitle Lists by degree

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

#! @Section Level
#! @SectionTitle Lists by level

#! @Arguments p, lambda
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  with level exactly $p^\lambda$.
DeclareGlobalFunction( "SL2Reps_PrimePowerIrrepsOfLevel" );

#! @Section Exceptions
#! @SectionTitle Lists of exceptional representations

#! @Returns a list of records of the form `rec(S, T, degree, level, name)`
#! @Description
#!  Constructs a list of the 18 exceptional irreps of $\mathrm{SL}_2(\mathbb{Z})$.
DeclareGlobalFunction( "SL2Reps_Exceptions" );
