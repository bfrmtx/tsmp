/***************************************************************************
 *   Copyright (C) 2003 by Bernhard Friedrichs                             *
 *   bfriedr@isi2.metronix.de                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef ADU06CAL_H
#define ADU06CAL_H
#include "allinclude.h"
#include "icalibration.h"
#include "my_valarray.h"
#include "my_valarray_double.h"
#include <string>
#include "string_utils.h"
#include <algorithm>
using namespace std;

/**
@author Bernhard Friedrichs
*/
class ADU06cal : public ICalibration {
public:
    ADU06cal();

    ~ADU06cal();

    virtual size_t read(const valarray<double>& newfreq, string chopper, const double& samplefreq, const double& gain,
                        const double& factor, const double& offset, string adc, string tags);
    virtual size_t dump(const string& fname);
    /*!
    return frequencies
    */
    virtual valarray<double>& freq();

    virtual void set_filename(const string& fname) {
        this->filename = fname;
    }

    /*!
    return complex transfer function
    */
    virtual valarray<complex<double> >& trfkt();



    /*!
    interpolate to a given list
    */
    virtual size_t interpolate(const valarray<double>& newfreq);

    /*!
    interpolate and extrapolate
    */
    virtual size_t interpolate_extend_theoretical(const valarray<double>& newfreq, string chopper, const double& samplefreq, const double& gain,
            const double& factor, const double& offset, string adc, string tags);

    /*!
    use a build in transfer function
    */
    virtual size_t hardwired(const valarray<double>& newfreq, string chopper, const double& samplefreq, const double& gain, const double& factor,
                             const double& offset, string adc, string tags);

    /*!
    use a theoretical function
    */
    virtual size_t theoretical(const valarray<double>& newfreq, string chopper, const double& samplefreq, const double& gain, const double& factor,
                               const double& offset, string adc, string tags);




    virtual size_t set_calib(const valarray<double>& f, const valarray<complex<double> >& trfn);
    virtual size_t set_calib(const valarray<double>& f, const valarray<double>& ampl, const valarray<double>& ph_deg);

		virtual string my_name() {return name;}


private:
    valarray<double> freqs;
    valarray<complex<double> > trf;
    string filename;
		string name;

};

#endif
