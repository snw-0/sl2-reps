#
# SL2Reps: Constructs representations of SL2(Z).
#

# Lists of Representations.

#! @Chapter Lists
#! @ChapterTitle Lists of representations

#! The **degree** of a representation is also known as the **dimension**.
#! The **level** of the congruent representation determined by the pair $(S,T)$ is equal to the order of $T$.
#!
#! We assign to each representation a **name** according to the conventions of <Cite Key="NW76"/>.

#! @Section Degree
#! @SectionTitle Lists by degree

#! @Arguments degree
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`.
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  that have the given degree and prime power level.
DeclareGlobalFunction( "SL2PrimePowerIrrepsOfDegree" );

#! @Arguments maximum_degree
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`.
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  that have at most the given maximum degree and prime power level.
DeclareGlobalFunction( "SL2PrimePowerIrrepsOfDegreeAtMost" );

#! @Arguments degree
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`.
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  that have the given degree.
DeclareGlobalFunction( "SL2IrrepsOfDegree" );

#! @Arguments maximum_degree
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`.
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  that have at most the given maximum degree.
DeclareGlobalFunction( "SL2IrrepsOfDegreeAtMost" );

#! @Section Level
#! @SectionTitle Lists by level

#! @Arguments p, lambda
#! @Returns a list of records of the form `rec(S, T, degree, level, name)`.
#! @Description
#!  Constructs a list of all irreps of $\mathrm{SL}_2(\mathbb{Z})$
#!  with level exactly $p^\lambda$.
DeclareGlobalFunction( "SL2PrimePowerIrrepsOfLevel" );

#! @Section Exceptions
#! @SectionTitle Lists of exceptional representations

#! @Returns a list of records of the form `rec(S, T, degree, level, name)`.
#! @Description
#!  Constructs a list of the 18 exceptional irreps of $\mathrm{SL}_2(\mathbb{Z})$.
DeclareGlobalFunction( "SL2ExceptionalIrreps" );
