//--------------------------------*- C++ -*----------------------------------
// blockMesh :  Block mesh description file
//
// adapted from:
// http://www.cfd-online.com/Forums/openfoam-meshing-blockmesh/61796-help-could-anyone-post-simple-cylinder-mesh.html
//
// JJS, 1/8/16
//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
FoamFile
{
  version  2.0;
  format   ascii;
  class dictionary;
  object blockMeshDict;
}
// ************************************
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
define(calcint, [esyscmd(perl -e 'printf int($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

   convertToMeters 0.001;

   define(D2, 30.0)   // outer cylinder diameter
   define(L1, 300.0)  // outer cylinder length
   define(D1, calc(0.707*D2))   //  inner cylinder diameter
   
   define(PI, 3.14159265)
   
   define(R1, calc(D1/2))
   define(R2, calc(D2/2))
   define(CW, calc(D1/4)) //Width of middle square section
   
   define(CY1, calc(R1*cos((PI/180)*45)))
   define(CZ1, calc(R1*sin((PI/180)*45)))

   define(CY2, calc(R2*cos((PI/180)*45)))
   define(CZ2, calc(R2*sin((PI/180)*45)))

   define(NPS, 8) //how many cells in the square section
   define(NPD, 8) //how many cells from square section to perimeter
   define(NPX, 200) // how many cells from top to bottom

   vertices
   (
    //Left surface
    //==========

    //square
    (0.0  CW  CW) vlabel(one_l)
    (0.0 -CW  CW) vlabel(two_l)
    (0.0 -CW -CW) vlabel(three_l)
    (0.0  CW -CW) vlabel(four_l)

    //inner circle
    (0.0  CY1  CZ1) vlabel(five_l)
    (0.0 -CY1  CZ1) vlabel(six_l)
    (0.0 -CY1 -CZ1) vlabel(seven_l)
    (0.0  CY1 -CZ1) vlabel(eight_l)

    //outer circle
    (0.0  CY2  CZ2) vlabel(nine_l)
    (0.0 -CY2  CZ2) vlabel(ten_l)
    (0.0 -CY2 -CZ2) vlabel(eleven_l)
    (0.0  CY2 -CZ2) vlabel(twelve_l)
   

    //Right surface
    //==========
    
    //square
    (L1  CW  CW) vlabel(one_r)
    (L1 -CW  CW) vlabel(two_r)
    (L1 -CW -CW) vlabel(three_r)
    (L1  CW -CW) vlabel(four_r)

    //inner circle
    (L1  CY1  CZ1) vlabel(five_r)
    (L1 -CY1  CZ1) vlabel(six_r)
    (L1 -CY1 -CZ1) vlabel(seven_r)
    (L1  CY1 -CZ1) vlabel(eight_r)

    //outer circle
    (L1  CY2  CZ2) vlabel(nine_r)
    (L1 -CY2  CZ2) vlabel(ten_r)
    (L1 -CY2 -CZ2) vlabel(eleven_r)
    (L1  CY2 -CZ2) vlabel(twelve_r)

   );				

   blocks
   (
    //square block
    //===========

    hex (
       one_l two_l three_l four_l
       one_r two_r three_r four_r
       )
    (NPS NPS NPX)
    simpleGrading (1 1 1)

    //inner circle
    //===========

    //slice1
    hex (
	five_l six_l two_l one_l
	five_r six_r two_r one_r
       )
    (NPS NPD NPX)
    simpleGrading (1 1 1)
    
    //slice2
    hex (
	six_l seven_l three_l two_l
	six_r seven_r three_r two_r
       )
    (NPS NPD NPX)
    simpleGrading (1 1 1)
    
    //slice3
    hex (
	seven_l eight_l four_l three_l
	seven_r eight_r four_r three_r
       )
    (NPS NPD NPX)
    simpleGrading (1 1 1)
    
    //slice4
    hex (
	eight_l five_l one_l four_l
	eight_r five_r one_r four_r
       )
    (NPS NPD NPX)
    simpleGrading (1 1 1)
   
 
    //outer circle
    //===========

   //slice1
    hex (
	nine_l ten_l six_l five_l
	nine_r ten_r six_r five_r
       )
    (NPS NPD NPX)
    simpleGrading (1 1 1)
    
    //slice2
    hex (
	ten_l eleven_l seven_l six_l
	ten_r eleven_r seven_r six_r
       )
    (NPS NPD NPX)
    simpleGrading (1 1 1)
    
    //slice3
    hex (
	eleven_l twelve_l eight_l seven_l
	eleven_r twelve_r eight_r seven_r
       )
    (NPS NPD NPX)
    simpleGrading (1 1 1)
    
    //slice4
    hex (
	twelve_l nine_l five_l eight_l
	twelve_r nine_r five_r eight_r
       )
    (NPS NPD NPX)
    simpleGrading (1 1 1)
   );


   //create the quarter circles
   edges
   (
    arc five_l six_l       (0.0 0.0 R1)
    arc six_l seven_l      (0.0 -R1 0.0)
    arc seven_l eight_l    (0.0 0.0 -R1)
    arc eight_l five_l     (0.0 R1 0.0)
    
    arc nine_l ten_l       (0.0 0.0 R2)
    arc ten_l eleven_l     (0.0 -R2 0.0)
    arc eleven_l twelve_l  (0.0 0.0 -R2)
    arc twelve_l nine_l    (0.0 R2 0.0)

    arc five_r six_r       (L1 0.0 R1)
    arc six_r seven_r      (L1 -R1 0.0)
    arc seven_r eight_r    (L1 0.0 -R1)
    arc eight_r five_r     (L1 R1 0.0)
    
    arc nine_r ten_r       (L1 0.0 R2)
    arc ten_r eleven_r     (L1 -R2 0.0)
    arc eleven_r twelve_r  (L1 0.0 -R2)
    arc twelve_r nine_r    (L1 R2 0.0)
   );

   patches
   (
    patch inlet
    (
	//outer circle patches
	(nine_l ten_l six_l five_l)
	(ten_l eleven_l seven_l six_l)
	(eleven_l twelve_l eight_l seven_l)
	(twelve_l nine_l five_l eight_l)
    )

    patch outlet
    (
       	//square patch - inlet
	(one_l two_l three_l four_l)

	//outer circle patches - outlet
	(nine_r ten_r six_r five_r)
	(ten_r eleven_r seven_r six_r)
	(eleven_r twelve_r eight_r seven_r)
	(twelve_r nine_r five_r eight_r)

    )

    wall wallbdry
    (

        //tube walls
	(nine_l ten_l ten_r nine_r)
	(ten_l eleven_l eleven_r ten_r)
	(eleven_l twelve_l twelve_r eleven_r)
	(twelve_l nine_l nine_r twelve_r)

	//square patch - outlet
	(one_r two_r three_r four_r)
	
	//inner circle patches - outlet
	(five_r six_r two_r one_r)
	(six_r seven_r three_r two_r)
	(seven_r eight_r four_r three_r)
	(eight_r five_r one_r four_r)

	//inner circle patches - inlet region
	(five_l six_l two_l one_l)
	(six_l seven_l three_l two_l)
	(seven_l eight_l four_l three_l)
	(eight_l five_l one_l four_l)

    )

);

mergePatchPairs
(
);
