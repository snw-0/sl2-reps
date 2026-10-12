gap> START_TEST( "Irrep position test, table classes reordered" );

# a character table of G that does not list the classes as ConjugacyClasses(G) does
gap> G := SL2WithConjClasses(5,1);;
gap> ccl := ConjugacyClasses(G);;
gap> perm := [ 1, 2, 9, 5, 3, 7, 6, 4, 8 ];;
gap> tbl := Objectify( NewType( NearlyCharacterTablesFamily,
>                 IsOrdinaryTable and IsAttributeStoringRep ), rec() );;
gap> SetUnderlyingGroup( tbl, G );
gap> SetUnderlyingCharacteristic( tbl, 0 );
gap> SetConjugacyClasses( tbl, ccl{ perm } );
gap> SetIdentificationOfConjugacyClasses( tbl, perm );
gap> SetOrdinaryCharacterTable( G, tbl );
gap> IsIdenticalObj( UnderlyingCharacterTable( Irr(G)[1] ), tbl );
true
gap> _SL2TestPositionsInGroup( G, 5, 1 );
true
gap> STOP_TEST( "classorder.tst" );
