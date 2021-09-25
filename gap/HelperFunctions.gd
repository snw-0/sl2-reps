#
# SL2Reps: Constructs representations of SL2(Z).
#

#! @Chapter Introduction
#!
#! This package constructs matrix presentations of the representations of $\mathrm{SL}_2(\mathbb{Z})$.

#! @Chapter Usage
#!
#! @Section InfoClass
#!  This package uses an `InfoClass`, `InfoSL2Reps`. It may be set to `0` (silent), `1` (info), or `2` (verbose).
#!  To change it, use `SetInfoLevel(InfoSL2Reps, k)`.

DeclareInfoClass( "InfoSL2Reps" );

DeclareGlobalFunction( "_SL2Reps_SqrtOfRootOfUnity" );
DeclareGlobalFunction( "_SL2Reps_SomeQuadraticNonResidue" );
DeclareGlobalFunction( "_SL2Reps_Factorizations" );

DeclareGlobalFunction( "_SL2Reps_RecordIrrep" );

DeclareGlobalFunction( "_SL2Reps_ConjClassesOdd" );
DeclareGlobalFunction( "_SL2Reps_ConjClassesEven" );
DeclareGlobalFunction( "_SL2Reps_ConjClasses" );
# DeclareGlobalFunction( "_SL2Reps_SL2Conj" );

DeclareGlobalFunction( "_SL2Reps_CharNorm" );
DeclareGlobalFunction( "_SL2Reps_RepChi" );
DeclareGlobalFunction( "_SL2Reps_ClassMap" );
DeclareGlobalFunction( "_SL2Reps_ChiTest" );
