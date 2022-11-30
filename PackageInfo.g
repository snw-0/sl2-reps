#
# SL2Reps: Constructing symmetric representations of SL(2,Z).
#


SetPackageInfo( rec(

PackageName := "SL2Reps",
Subtitle := "Constructing symmetric representations of SL(2,Z).",
Version := "1.1",
Date := "26/11/2022", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    FirstNames := "Siu-Hung",
    LastName := "Ng",
    WWWHome := "https://www.math.lsu.edu/~rng/",
    Email := "rng@math.lsu.edu",
    IsAuthor := true,
    IsMaintainer := false,
    PostalAddress := "Louisiana State University, Baton Rouge, LA, 70803, USA",
    Place := "Baton Rouge, LA, USA",
    Institution := "Louisiana State University",
  ),
  rec(
    FirstNames := "Yilong",
    LastName := "Wang",
    WWWHome := "https://yilongwang11.github.io",
    Email := "wyl@bimsa.cn",
    IsAuthor := true,
    IsMaintainer := false,
    PostalAddress := "Louisiana State University, Baton Rouge, LA, 70803, USA <Br/> Current: Beijing Institute of Mathematical Sciences and Applications (BIMSA), Huairou, Beijing, China",
    #Place := TODO,
    Institution := "Louisiana State University; currently Beijing Institute of Mathematical Sciences and Applications (BIMSA)",
  ),
  rec(
    FirstNames := "Samuel",
    LastName := "Wilson",
    WWWHome := "https://snw-0.github.io",
    Email := "swil311@lsu.edu",
    IsAuthor := true,
    IsMaintainer := true,
    PostalAddress := "Louisiana State University, Baton Rouge, LA, 70803, USA",
    Place := "Baton Rouge, LA, USA",
    Institution := "Louisiana State University",
  ),
],

Status := "deposited",

SourceRepository := rec( Type := "git", URL := "https://github.com/snw-0/sl2-reps" ),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome := "https://snw-0.github.io/sl2-reps",
PackageInfoURL := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
README_URL     := Concatenation( ~.PackageWWWHome, "/README" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                  "/releases/download/v", ~.Version,
                                  "/sl2-reps-", ~.Version ),
ArchiveFormats := ".tar.gz",

AbstractHTML :=
  "The <span class=\"pkgname\">SL2Reps</span> package provides provides methods for constructing and testing matrix presentations of the representations of SL(2,Z) whose kernels are congruence subgroups of SL(2,Z).",

PackageDoc := rec(
  BookName  := "SL2Reps",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Constructing symmetric representations of SL(2,Z).",
),

Dependencies := rec(
  GAP := ">= 4.10",
  NeededOtherPackages := [ ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

Keywords := [ "representations" ],

AutoDoc := rec(
  TitlePage := rec(
    Copyright := """
      <Index>License</Index>
      &copyright; 2022 by Siu-Hung Ng, Yilong Wang, and Samuel Wilson<P/>
      This package is free software;
      you can redistribute it and/or modify it under the terms of the
      <URL Text="GNU General Public License">https://www.fsf.org/licenses/gpl.html</URL>
      as published by the Free Software Foundation; either version 2 of the License,
      or (at your option) any later version.
      """,
    Acknowledgements := """
      This project is partially supported by NSF grant DMS 1664418.
      """,
  ),
),

));
