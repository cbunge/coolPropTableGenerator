/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "openFoamTableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::interpolation2DTable2<Type>::readTable()
{
    fileName fName(fileName_);
    fName.expand();

    // Read data from file
    reader_()(fName, *this);

    if (this->empty())
    {
        FatalErrorInFunction
            << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are in ascending order
    checkOrder();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolation2DTable2<Type>::interpolation2DTable2()
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(),
    boundsHandling_(interpolation2DTable2::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(NULL)
    //isNull_(true)
{}


template<class Type>
Foam::interpolation2DTable2<Type>::interpolation2DTable2
(
    const List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& values,
    const boundsHandling bounds,
    const fileName& fName
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(values),
    boundsHandling_(bounds),
    fileName_(fName),
    reader_(NULL)
    //isNull_(isNull)
{}


template<class Type>
Foam::interpolation2DTable2<Type>::interpolation2DTable2(const fileName& fName)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(),
    boundsHandling_(interpolation2DTable2::WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>(dictionary()))
    //isNull_(false)
{
    readTable();
}


template<class Type>
Foam::interpolation2DTable2<Type>::interpolation2DTable2(const dictionary& dict)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(),
    boundsHandling_(wordToBoundsHandling(dict.lookup("outOfBounds"))),
    fileName_(dict.lookup("file")),
    reader_(tableReader<Type>::New(dict))
    //isNull_(false)
{
    readTable();
}


template<class Type>
Foam::interpolation2DTable2<Type>::interpolation2DTable2
(
     const interpolation2DTable2& interpTable
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(interpTable),
    boundsHandling_(interpTable.boundsHandling_),
    fileName_(interpTable.fileName_),
    reader_(interpTable.reader_)    // note: steals reader. Used in write().
    //isNull_(interpTable.isNull_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolation2DTable2<Type>::interpolateValue
(
    const List<Tuple2<scalar, Type>>& data,
    const scalar lookupValue
) const
{
    label n = data.size();

    scalar minLimit = data.first().first();
    scalar maxLimit = data.last().first();

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable2::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable2::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable2::CLAMP:
            {
                return data.first().second();
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable2::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable2::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable2::CLAMP:
            {
                return data.last().second();
                break;
            }
        }
    }

    // look for the correct range in X
    label lo = 0;
    label hi = 0;

    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= data[i].first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        return data[lo].second();
    }
    else
    {
        Type m =
            (data[hi].second() - data[lo].second())
           /(data[hi].first() - data[lo].first());

        // normal interpolation
        return data[lo].second() + m*(lookupValue - data[lo].first());
    }
}
/*
template<class Type>
Type Foam::interpolation2DTable2<Type>::inverseInterpolate
(
    const List<Tuple2<scalar, Type> >& data,
    const scalar lookupValue
) const
{
    label n = data.size();

    scalar minLimit = data.first().second();
    scalar maxLimit = data.last().second();

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable2::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable2::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable2::CLAMP:
            {
                return data.first().second();
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable2::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable2::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable2::CLAMP:
            {
                return data.last().second();
                break;
            }
        }
    }


    // look for the correct range in X
    label lo = 0;
    label hi = 0;

    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= data[i].second())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        return data[lo].second();
    }
    else
    {
        Type m =
            (data[hi].first() - data[lo].first())
           /(data[hi].second() - data[lo].second());

        // inverse interpolation
        return data[lo].first() + m*(lookupValue - data[lo].second());
    }
}
*/

template<class Type>
template<class BinaryOp>
Foam::label Foam::interpolation2DTable2<Type>::Xi
(
    const BinaryOp& bop,
    const scalar valueX,
    const bool reverse
) const
{
    const table& t = *this;

    label limitI = 0;
    if (reverse)
    {
        limitI = t.size() - 1;
    }

    if (bop(valueX, t[limitI].first()))
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable2::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << valueX << ") out of bounds"
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable2::WARN:
            {
                WarningInFunction
                    << "value (" << valueX << ") out of bounds"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable2::CLAMP:
            {
                return limitI;
            }
            default:
            {
                FatalErrorInFunction
                    << "Un-handled enumeration " << boundsHandling_
                    << abort(FatalError);
            }
        }
    }

    label i = 0;
    if (reverse)
    {
        label nX = t.size();
        i = 0;
        while ((i < nX) && (valueX > t[i].first()))
        {
            i++;
        }
    }
    else
    {
        i = t.size() - 1;
        while ((i > 0) && (valueX < t[i].first()))
        {
            i--;
        }
    }

    return i;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolation2DTable2<Type>::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{
    // Considers all of the list in Y being equal
    label nX = this->size();

    const table& t = *this;

    if (nX == 0)
    {
        WarningInFunction
            << "cannot interpolate a zero-sized table - returning zero" << endl;

        return Zero;
    }
    else if (nX == 1)
    {
        // only 1 column (in X) - interpolate to find Y value
        return interpolateValue(t.first().second(), valueY);
 
    }
    else
    {
        // have 2-D data, interpolate

        // find low and high indices in the X range that bound valueX
        label x0i = Xi(lessOp<scalar>(), valueX, false);
        label x1i = Xi(greaterOp<scalar>(), valueX, true);

        if (x0i == x1i) 
        {
            return interpolateValue(t[x0i].second(), valueY);
        }
        else
        {
            Type y0(interpolateValue(t[x0i].second(), valueY));
            Type y1(interpolateValue(t[x1i].second(), valueY));

            // gradient in X
            scalar x0 = t[x0i].first();
            scalar x1 = t[x1i].first();
            Type mX = (y1 - y0)/(x1 - x0);

            // interpolate
            return y0 + mX*(valueX - x0);
        }
    }
}


template<class Type>
Foam::word Foam::interpolation2DTable2<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case interpolation2DTable2::ERROR:
        {
            enumName = "error";
            break;
        }
        case interpolation2DTable2::WARN:
        {
            enumName = "warn";
            break;
        }
        case interpolation2DTable2::CLAMP:
        {
            enumName = "clamp";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::interpolation2DTable2<Type>::boundsHandling
Foam::interpolation2DTable2<Type>::wordToBoundsHandling
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return interpolation2DTable2::ERROR;
    }
    else if (bound == "warn")
    {
        return interpolation2DTable2::WARN;
    }
    else if (bound == "clamp")
    {
        return interpolation2DTable2::CLAMP;
    }
    else
    {
        WarningInFunction
            << "bad outOfBounds specifier " << bound << " using 'warn'" << endl;

        return interpolation2DTable2::WARN;
    }
}

/*
template<class Type>
inline Foam::interpolation2DTable2<Type>&
Foam::interpolation2DTable2<Type>::operator=
(
    const interpolation2DTable2<Type>& et
)
{
    boundsHandling_ = et.boundsHandling_;
    fileName_ = et.fileName_;
    reader_ = et.reader_;
    isNull_ = et.isNull_;

    return *this;
}


template<class Type>
inline Foam::interpolation2DTable2<Type> Foam::operator+
(
    const interpolation2DTable2<Type>& et1,
    const interpolation2DTable2<Type>& et2
)
{
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > etn = et1, ett = et2;

    if (et1.size() != et2.size())
    {
	FatalErrorInFunction
            << "attempt to sum list with different sizes."
            << abort(FatalError);
    }

    if (et1.isNull_)
    {
    	return interpolation2DTable2<Type>
        (
            ett,
    	    et2.boundsHandling_,
    	    et2.fileName_,
    	    et2.isNull_
        );
    }
    else if (et2.isNull_)
    {
    	return interpolation2DTable2<Type>
        (
            etn,
    	    et1.boundsHandling_,
    	    et1.fileName_,
    	    et1.isNull_
        );
    }

    for (int i = 0 ; i < etn.size() ; ++i)
    {
	for (int j = 0 ; j < etn[i].second().size() ; ++j)
	{
	    etn[i].second()[j].second() += ett[i].second()[j].second();
	}
    }

    return interpolation2DTable2<Type>
    (
	etn,
	et1.boundsHandling_,
	et1.fileName_,
	false
    );
}


template<class Type>
inline Foam::interpolation2DTable2<Type> Foam::operator-
(
    const interpolation2DTable2<Type>& et1,
    const interpolation2DTable2<Type>& et2
)
{
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > etn = et1, ett = et2;

    if (et1.size() != et2.size())
    {
	FatalErrorInFunction
            << "attempt to sum list with different sizes."
            << abort(FatalError);
    }

    for (int i = 0 ; i < etn.size() ; ++i)
    {
	for (int j = 0 ; j < etn[i].second().size() ; ++j)
	{
	    etn[i].second()[j].second() -= ett[i].second()[j].second();
	}
    }

    return interpolation2DTable2<Type>
    (
	 etn,
	 et1.boundsHandling_,
	 et1.fileName_,
	 et1.isNull_
    );
}


template<class Type>
inline Foam::interpolation2DTable2<Type> Foam::operator*
(
    const scalar s,
    const interpolation2DTable2<Type>& et
)
{
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > etn = et;

    // I'm not sure that it saves time but I try it.
    if (s == 1 || et.isNull_)
    {
    	return interpolation2DTable2<Type>
        (
    	    etn,
    	    et.boundsHandling_,
    	    et.fileName_,
	    et.isNull_
        );
    }
    else
    {
    	for (int i = 0 ; i < etn.size() ; ++i)
    	{
    	    for (int j = 0 ; j < etn[i].second().size() ; ++j)
    	    {
    		etn[i].second()[j].second() *= s;
    	    }
    	}
    }

    if (s == 0)
    {
    	return interpolation2DTable2<Type>
    	(
    	    etn,
    	    et.boundsHandling_,
    	    et.fileName_,
    	    true
    	);
    }
    else
    {
	return interpolation2DTable2<Type>
	(
	    etn,
	    et.boundsHandling_,
	    et.fileName_,
	    et.isNull_
	);
    }
}

*/
template<class Type>
typename Foam::interpolation2DTable2<Type>::boundsHandling
Foam::interpolation2DTable2<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;
    return prev;
}


template<class Type>
void Foam::interpolation2DTable2<Type>::checkOrder() const
{
    label n = this->size();
    const table& t = *this;

    scalar prevValue = t[0].first();

    for (label i=1; i<n; ++i)
    {
        const scalar currValue = t[i].first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class Type>
void Foam::interpolation2DTable2<Type>::write(Ostream& os) const
{
    os.writeKeyword("file")
        << fileName_ << token::END_STATEMENT << nl;
    os.writeKeyword("outOfBounds")
        << boundsHandlingToWord(boundsHandling_) << token::END_STATEMENT << nl;

    *this >> os;
}


// ************************************************************************* //
