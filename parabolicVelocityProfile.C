/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    setHotRoom

Description
    Set the initial field of T for the hot room problem.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OSspecific.H"
#include "fixedValueFvPatchFields.H"
#include "fvPatchFieldFields.H"

#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

  argList::validArgs.append("boundaryName");
  argList::validArgs.append("maximum Velocity");
  argList::validOptions.insert("z_noncenter","");
  argList::validOptions.insert("y_noncenter","");

#   include "setRootCase.H"

  scalar maxVel(readScalar(IStringStream(args.args()[4])()));
  word bcName(args.args()[3]);

  bool z_noncenter=false;
  bool y_noncenter=false;

  if(args.options().found("z_noncenter")) {
    z_noncenter=true;
  }
  if(args.options().found("y_noncenter")) {
    y_noncenter=true;
  }

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  fvPatchVectorFieldField& Upatches = U.boundaryField();

  forAll(Upatches, patchI)
   {
     if
       (
        (typeid(Upatches[patchI]) == typeid(fixedValueFvPatchVectorField))
        &&
        (mesh.boundaryMesh()[patchI].name() == bcName) 
	)
       {
	 Info << " Patching inlet\n";
      
	 const vectorField& faceCentres = 
	   mesh.Cf().boundaryField()[patchI];
	 
	 fixedValueFvPatchVectorField& Upatch =
            refCast<fixedValueFvPatchVectorField>(Upatches[patchI]);

	 scalar maxY=max(faceCentres.component(1));
	 scalar minY=min(faceCentres.component(1));
	 scalar maxZ=max(faceCentres.component(2));
	 scalar minZ=min(faceCentres.component(2));

	 Info << " y [ " << minY << " , " << maxY << " ]  z [ " << minZ << " , " << maxZ << " ] \n";

	 scalar lenY=maxY-minY;
	 scalar lenZ=maxZ-minZ;
	 scalar offY=minY;
	 scalar offZ=minZ;
	 if(y_noncenter) {
	   offY -= lenY;
	   lenY *= 2;
	 }
	 if(z_noncenter) {
	   offZ -= lenZ;
	   lenZ *= 2;
	 }

	 Info << " => y [ " << offY << " , " << lenY << " ]  z [ " << offZ << " , " << lenZ << " ]\n";
 
	 forAll(faceCentres, facei)
	   {
	     scalar y,z;
	     if(lenY>1e-8)
	       y=(faceCentres[facei].y()-offY)/(lenY)*2-1;
	     else 
	       y=0;

	     if(lenZ>1e-8)
	       z=(faceCentres[facei].z()-offZ)/(lenZ)*2-1;
	     else
	       z=0;

	     scalar vel=maxVel*(1-y*y)*(1-z*z);

	     //	     Info << y << " \t " << z << " \t " << vel << endl;

	     Upatch[facei] = vector(vel,0,0);
	   }
	 
       }
   };

 Info<< "Writing modified field U\n" << endl;
 U.write();

 Info<< "End\n" << endl;

 return(0);
}


// ************************************************************************* //
