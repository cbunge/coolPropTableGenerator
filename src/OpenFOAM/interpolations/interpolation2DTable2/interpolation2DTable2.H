/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::interpolation2DTable2

Description
    2D table interpolation. The data must be in ascending order in both
    dimensions x and y.

SourceFiles
    interpolation2DTable2.C

\*---------------------------------------------------------------------------*/

#ifndef interpolation2DTable2_H
#define interpolation2DTable2_H

#include "List.H"
#include "Tuple2.H"
#include "tableReader.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
/*
template<class Type> class interpolation2DTable2;

template<class Type>
inline interpolation2DTable2<Type> operator+
(
    const interpolation2DTable2<Type>&,
    const interpolation2DTable2<Type>&
);


template<class Type>
inline interpolation2DTable2<Type> operator-
(
    const interpolation2DTable2<Type>&,
    const interpolation2DTable2<Type>&
);


template<class Type>
inline interpolation2DTable2<Type> operator*
(
    const scalar,
    const interpolation2DTable2<Type>&
);

*/
/*---------------------------------------------------------------------------*\
                    Class interpolation2DTable2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolation2DTable2
:
    public List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound values
        enum boundsHandling
        {
            ERROR,          //!< Exit with a FatalError
            WARN,           //!< Issue warning and clamp value (default)
            CLAMP           //!< Clamp value to the start/end value
        };

        //- Cconvenience typedef
        typedef List<Tuple2<scalar, List<Tuple2<scalar, Type>>>> table;


private:

    // Private data

        //- Enumeration for handling out-of-bound values
        boundsHandling boundsHandling_;

        //- File name
        fileName fileName_;

        //- The actual reader
        autoPtr<tableReader<Type>> reader_;

	//- true or yes if all the values are 0
	//Switch isNull_;

    // Private Member Functions

        //- Read the table of data from file
        void readTable();

        //- Return interpolated value in List
        Type interpolateValue
        (
            const List<Tuple2<scalar, Type>>& data,
            const scalar
        ) const;
/*
        //- Return interpolated value in List
        Type inverseInterpolate
        (
            const List<Tuple2<scalar, Type> >& data,
            const scalar
        ) const;
*/

        //- Return an X index from the matrix
        template<class BinaryOp>
        label Xi
        (
            const BinaryOp& bop,
            const scalar valueX,
            const bool reverse
        ) const;


public:

    // Constructors

        //- Construct null
        interpolation2DTable2();

        //- Construct from components
        interpolation2DTable2
        (
            const List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& values,
            const boundsHandling bounds,
            const fileName& fName
	    //const Switch isNull
        );

        //- Construct given the name of the file containing the table of data
        interpolation2DTable2(const fileName& fName);

        //- Construct by reading the fileName and boundsHandling from dictionary
        interpolation2DTable2(const dictionary& dict);

        //- Construct copy
        interpolation2DTable2(const interpolation2DTable2& interpTable);


    // Member Functions

        //- Return the out-of-bounds handling as a word
        word boundsHandlingToWord(const boundsHandling& bound) const;

        //- Return the out-of-bounds handling as an enumeration
        boundsHandling wordToBoundsHandling(const word& bound) const;

        //- Set the out-of-bounds handling from enum, return previous setting
        boundsHandling outOfBounds(const boundsHandling& bound);

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void checkOrder() const;

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        //- Return an element of constant Tuple2<scalar, Type>
        const List<Tuple2<scalar, Type>>& operator[](const label) const;

        //- Return an interpolated value
        Type operator()(const scalar, const scalar) const;

	//inline interpolation2DTable2& operator=(const interpolation2DTable2&);
/*
    // Friend operators

        friend interpolation2DTable2 operator+ <Type>
        (
            const interpolation2DTable2&,
            const interpolation2DTable2&
        );

        friend interpolation2DTable2 operator- <Type>
        (
            const interpolation2DTable2&,
            const interpolation2DTable2&
        );

        friend interpolation2DTable2 operator* <Type>
	(
            const scalar,
            const interpolation2DTable2&
        );*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation2DTable2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
