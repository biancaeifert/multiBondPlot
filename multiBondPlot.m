(* ::Package:: *)

(* :Title: Multiple Bond Plot *)

(* :Context: multiBondPlot` *)

(* :Author: Bianca Eifert *)

(* :Summary: 
	This package provides the function multiBondPlot, which generates plots of molecules with double and triple bonds.
    It is essentially a packaged version of this Demonstration:
    http://demonstrations.wolfram.com/DisplayingMoleculesWithMultipleBonds/
*)

(* :Package Version: 1.0 *)

(* :Mathematica Version: 9.0+ *)

(* :Copyright: Bianca Eifert, under the MIT license *)

(* :History:  *)

(* :Keywords:  *)

(* :Limitations: 
    Input (member of ChemicalData[] or Chemcial Table File) must contain all bonds and information on their multiplicity.
*)

(* :Discussion:  *)


BeginPackage["multiBondPlot`"];

multiBondPlot::usage="multiBondPlot[molecule] returns a plot of a molecule with single, double and/or triple bonds. \
The argument is a substance from ChemicalData or the name of a MOL file. Options of Graphics3D can be given.";

Begin["`Private`"];


(*the underlying plotting engine: *)
mbp[atomPositions_,vertexTypes_,edgeRules_,edgeTypes_,options___]:=
Module[{atomplot,bondlist,shift,bond,bondplot},

(*Graphics3D directives for atoms; colours are chosen based on atom type (which will always be available for the types of input used here), 
radii are all the same except for hydrogen because smaller H atoms tend to look better in the overall plot: *)
atomplot={
    ColorData["Atoms",vertexTypes[[#]]],
    Sphere[atomPositions[[#]],If[vertexTypes[[#]]=="H",25,35]]
}&/@Range[Length[vertexTypes]];

(*types of bonds, and atom pairs connected by them: *)
bondlist=Flatten/@(Transpose[{edgeRules,edgeTypes}]/.Rule->List/.{"Single"->1,"Double"->2,"Triple"->3});

(*construct the bonds: *)
(*shift bonds away from the vector connecting the two atoms to make room for several bonds: *)
shift[bondindex_]:=Module[{thisbond,dir,spacing},
    (*get the bond from the list of bonds: *)
    thisbond=bondlist[[bondindex]];
    (*the shift vector must be perpendicular to the bond itself and to the local molecular plane;
    first, identify neighbours of the bond to establish the local molecular plane: *)
    dir=Select[
        bondlist[[Delete[Range[Length[bondlist]],bondindex],;;2]],
        !Intersection[#,thisbond[[;;2]]]==={}&
    ];
    (*... now construct a vector that is perpendicular to that plane and the bond;
    if there aren't enough neighbouring atoms to establish a plane, just assume a fixed vector of {0,0,1}: *)
    dir=If[Length[dir]<2,
        {0,0,1},
        Normalize[Cross[Subtract@@atomPositions[[dir[[1]]]],Subtract@@atomPositions[[dir[[2]]]]]]
    ];
    (*fixed parameter for the spacing of the bonds, i.e. how far they are shifted out from the vector connecting the two atoms;
    a value of 7 seems to look good with the chosen atom radii (see "atomplot" above), so that's what we'll use: *)
    spacing=7;
    (*apply the shift: *)
    thisbond[[3]]*spacing*Normalize[Cross[Subtract@@atomPositions[[thisbond[[;;2]]]],dir]]
];

(*Graphics3D directives for one bond of given multiplicity; this is either a single bond or one of the tubes of a multiple bond: *)
bond[bondindex_,partials_]:=Module[{boundatoms,totalshift,bondrad},
    (*coordinates of the two atoms involved in the bond: *)
    boundatoms=atomPositions[[bondlist[[bondindex,1;;2]]]];
    (*calculate the overall shift for this bond: *)
    totalshift=partials*shift[bondindex];
    (*bond radius is hard-coded to a value that looks good with the chosen atom radii; scales with multiplicity: *)
    bondrad=12-2*bondlist[[bondindex,3]];
    {
        ColorData["Atoms",vertexTypes[[bondlist[[bondindex,#]]]]],
        Tube[{boundatoms[[#]]+totalshift,.5*Total[boundatoms]+totalshift},bondrad]
    }&/@{1,2}
];

(*all bonds: *)
bondplot=Table[
    bond[#,partials],
    {partials,Switch[bondlist[[#,3]],3,{-1,0,1},2,{-1,1},1,{0}]}
]&/@Range[Length[bondlist]];

(*return the plot: *)
Graphics3D[{atomplot,bondplot},
    FilterRules[Flatten[{options}],Options[Graphics3D]],
    Boxed->False,
    SphericalRegion->True,
    Lighting->"Neutral",
    ViewPoint->{1.3,-2.4,-2.},
    ImageSize->400,
    BaseStyle->{Specularity[White,100]}
]
];


(*from ChemicalData; 'options' can contain any options of Graphics3D: *)
multiBondPlot[molecule_String,options___]:=
Catch[Module[{datasets},

(*return Missing if any of the required data is not available: *)
If[
    Or[
        Head[ChemicalData[molecule]]===ChemicalData,
        MemberQ[Head[ChemicalData[molecule,#]]&/@{"VertexTypes","AtomPositions","EdgeRules","EdgeTypes"},Missing|ChemicalData],
        Head[ChemicalData[molecule,"AtomPositions"][[1]]]===Missing
    ],
    Throw[Missing["NotAvailable"]]
];

(*data can have varying depths depending on the molecule in question, correct for that: *)
datasets=If[MatchQ[ChemicalData[molecule,"VertexTypes"],{{_String..}..}],1,All];

(*pass data to the plotting engine: *)
mbp[
    If[$VersionNumber<=9.,
        ChemicalData[molecule,"AtomPositions"][[datasets]],
        QuantityMagnitude[ChemicalData[molecule,"AtomPositions"][[datasets]]]
    ],
    ChemicalData[molecule,"VertexTypes"][[datasets]],
    ChemicalData[molecule,"EdgeRules"][[datasets]],
    ChemicalData[molecule,"EdgeTypes"][[datasets]],
    options
]
]];


(*from file; 'options' can contain any options of Graphics3D: *)
multiBondPlot[molecule_String/;FileExistsQ[molecule],options___]:=
Catch[Module[{datasets},

(*return Missing if any of the required data is not available: *)
If[
    Or[
        MemberQ[Import[molecule,#]&/@{"VertexTypes","VertexCoordinates","EdgeRules","EdgeTypes"},{}],
        Length[Import[molecule,"VertexCoordinates"][[1]]]===2
    ],
    Throw[Missing["NotAvailable"]]
];

(*data can have varying depths depending on the molecule in question, correct for that: *)
datasets=If[MatchQ[Import[molecule,"VertexTypes"],{{_String..}..}],1,All];

(*pass data to the plotting engine: *)
mbp[
    Import[molecule,"VertexCoordinates"][[datasets]],
    Import[molecule,"VertexTypes"][[datasets]],
    Import[molecule,"EdgeRules"][[datasets]],
    Import[molecule,"EdgeTypes"][[datasets]],
    options
]
]];


End[];
EndPackage[];
