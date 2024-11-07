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
#include "mtxcal.h"

MTXCal::MTXCal()
    : ICalibration() {

    name = "mtxcal";
}


MTXCal::~MTXCal() {}


/*!
    \fn MTXCal::read(string filename, string chopper)
     reads from a file; chopper: on, off, auto, none
     simple check for double entry,
     makes sure that f[0] < f[1]
 */
size_t MTXCal::read(const valarray<double>& newfreq, string chopper = "on", const double& samplefreq = 0, const double& gain = 0,
                    const double& factor = 0,  const double& offset = 0, string adc = "", string tags = "") {
    if (!filename.size())
        return 0;
    ifstream ifs;


    ifs.open(filename.c_str());
    string share_dir;
    if (!ifs) {
        ifs.close();
        ifs.clear();
        share_dir =  "/usr/local/share/" + filename;
        ifs.open(share_dir.c_str());
        if (!ifs) {
            ifs.close();
            ifs.clear();
            cerr << "MTXCal::read -> can nor open file " << filename << " neither " << share_dir << endl;
            cerr << "try uppper case" << endl;
            upper(filename);

            ifs.open(filename.c_str());
            if (!ifs) {
                ifs.close();
                ifs.clear();
            }
            share_dir =  "/usr/local/share/" + filename;

            ifs.open(share_dir.c_str());
            if (!ifs) {
                cerr << "MTXCal::read -> can nor open file " << filename << " neither " << share_dir << " exit" << endl;
                cerr << "giving up" << endl; exit(0);
            }
            else {
                cerr << "MTXCal::read -> " << share_dir << endl;

            }
        } else
            cerr << "MTXCal::read -> local file not found, taking " << share_dir << " instead" << endl;
    }

    vector<double> ix;
    vector<complex<double> > iz;
    double dx, dy, dz;
    complex<double> z;
    string line;

    lower(chopper);

    if ( ( chopper == "on") || (chopper == "off") || (chopper == "none")) {
        // metronix uses chopper on or chopper off or chopper none
        chopper = "chopper " + chopper;    // is "chopper on" now for example
    }

    if ( (chopper == "auto") && (samplefreq > 512) )
        chopper = "chopper off";
    else if ( (chopper == "auto") && (samplefreq <= 512) )
        chopper = "chopper on";

    do {
        getline(ifs, line, '\n');
        lower(line);
        if ( (line.find(chopper) != line.npos))
            break;
    } while (ifs);

    if (ifs.eof()) {
        cerr << "MTXCal::read -> can not find key words chopper on or chopper off; I try a file with numbers only now.... check!!!" << endl;

        ifs.clear();
        ifs.seekg(0);
        chopper = "no chopper info, blank file, freq[hz], ampl [V/(nT*Hz)], phase[deg]";
    }


    cerr << "MTXCal::read -> using " << chopper << endl;
    while(ifs) {
        ifs >> dx >> dy >> dz;              // read temp variables
        dy *=  (dx * 1000.);                // mtx data is normalized by frequency and we do NOT keep that!!!!!!!!!

        // file should have V/(nT *Hz) as units; so multiply
        // with 1000 to get mV as (stored in ats files)
        // hence that mapros might uses nT instead of pT


        z = polar(dy ,(M_PI/180.) * dz);
        // if last operation was successful
        if(ifs) {
            ix.push_back(dx);
            iz.push_back(z);

        }

    }

    freqs.resize(ix.size());
    trf.resize(ix.size());

    // copy data into valarrays real and complex
    copy(&ix[0], &ix[ix.size()], &freqs[0]);
    copy(&iz[0], &iz[iz.size()], &trf[0]);

    // some checks
    size_t i = size_t ( adjacent_find(&freqs[0], &freqs[freqs.size()]) - &freqs[0]);

    if (i != freqs.size()) {
        cerr << "MTXCal::read -> at least on frequency with double entry!!! " << endl;
        exit(0);
    }

    // sort frequency if neccessaray
    if (freqs[1] > freqs[2]) {
        reverse(&freqs[0], &freqs[freqs.size()]);
        reverse(&trf[0], &trf[trf.size()]);
    }


    cerr << "MTXCal:: read trf: " << freqs[0] << " " << trf[0] << " " << freqs[freqs.size()-1] << " "  << trf[freqs.size()-1] << endl;
    return freqs.size();



}

/*!
    \fn MTXCal::dump(string filename)
     debug function; just dump values to file
 */

size_t MTXCal::dump(const string& filename) {

    ofstream ofs;
    size_t i;
    if (!filename.size()) {
        cerr << "MTXCal::dump -> no filename given" << endl;
        return 0;
    }
    ofs.open(filename.c_str());
    if (!ofs) {
        cerr << "MTXCal::dump -> not open output stream" << endl;
        return 0;

    }

    if (!freqs.size() || !trf.size()) {
        cerr << "MTXCal::dump -> no frequencies / complex data to dump" << endl;
        return 0;

    }

    for (i = 1; i < freqs.size(); i++) {

        ofs << freqs[i] << "  " << trf[i].real() << "  " << trf[i].imag() << "                     "
            << abs(trf[i]) << "  " << (180. / M_PI) * arg(trf[i]) << "               "
            << (abs(trf[i])  / freqs[i]) / 1000. << "  " << (180. / M_PI) * arg(trf[i]) << endl;

    }



    return freqs.size();
}




/*!
interpolate to a given list
f(wl/2 +1)
for ( i = 0; i < f.size(); i++) f[i] = i * samplefreq / wl;

 since induction coils having a frequency dependend output interpolation is bad for low frequencies
 where interpolation especially for tha phase is very bad; solution: normalize again by f, interpolate, and
 multiply by f again
*/
size_t MTXCal::interpolate(const valarray<double>& newfreq) {

    size_t i;
    valarray<double> b, c, d;

    // transfer function must be loaded here!!
    if (trf.size()) {


        //normalize for better interpolation

        for (size_t ii = 1; ii < trf.size(); ii++)
            trf[ii] /= freqs[ii];
        // can be DC, check
        if (freqs[0])
            trf[0] /= freqs[0];

        //copy original data into temp
        valarray<double>temp_re(trf.size()), temp_im(trf.size());
        my_valarray::cplx2ri(trf, temp_re, temp_im);



        // generate new size - empty
        trf.resize(newfreq.size());

        cerr << "MTXCal::interpolate -> interpolating sensor trf (" << newfreq[1] << " - "
             << newfreq[newfreq.size()-1] << " Hz), " << newfreq.size() << "pts" << flush << endl;

        //
        // use the old data to prepare
        my_valarray_double::aspline(freqs, temp_re, b, c, d);
        for (i = 0; i < trf.size(); i++)
            trf[i] = complex<double> (my_valarray_double::seval(newfreq[i], freqs, temp_re, b, c, d), 0.);

        // use the old data to prepare
        my_valarray_double::aspline(freqs, temp_im, b, c, d);
        for (i = 0; i < trf.size(); i++)
            trf[i] = complex<double> (trf[i].real(),
                                      my_valarray_double::seval(newfreq[i], freqs, temp_im, b, c, d));


    } else
        cerr << "MTXCal::interpolate -> interpolating ERROR" << endl;



    freqs.resize(newfreq.size());
    freqs = newfreq;

    // undo normalization
    for (size_t ii = 1; ii < trf.size(); ii++)
        trf[ii] *= freqs[ii];
    // can be DC, check
    if (freqs[0])
        trf[0] *= freqs[0];
    // else correct DC part
    else
        trf[0] = complex<double> (1.0, 0.);

    cerr << "MTXCal::interpolate -> " << freqs[0] << " " << trf[0] << ", " << freqs[1] << " " << trf[1] << " " << freqs[freqs.size()-1] << " "  << trf[freqs.size()-1] << flush << endl;


    return freqs.size();
}



/*!
interpolate and extrapolate
*/
size_t MTXCal::interpolate_extend_theoretical(const valarray<double>& newfreq, string chopper = "on", const double& samplefreq = 0,
                                              const double& gain = 0, const double& factor = 0,  const double& offset = 0, string adc = "", string tags = "") {

    valarray<double> old_f(freqs);  // save the frequencies
    MTXCal::interpolate(newfreq);



    // interpolate with theoretical transfer function instead of spline
    // outside give values from read in transfer function
    // phase wil jump about 0.1 deg @ 0.1 Hz ... where the measured trf function ends
    // there is no theoretical function for chopper off!
    // interpolation gives bad results for imaginary part!
    complex<double> p1, p2, im(0,1.);

    if (freqs[freqs.size()-1] < 513) {
        cerr << "MTXCal::interpolate_extend_theoretical -> extrapolating trf function for mfs coils 05 06 08U KIM" << endl;
        for (size_t i = 0; i < freqs.size(); i++) {
            if ( (freqs[i] < old_f.min()) || (freqs[i] > old_f.max()) ) {

                p1 = (freqs[i]/4.) * im;
                p2 = (freqs[i]/8192.) * im;
                trf[i] = 1000. * 0.8 * (  (p1 / ( 1. + p1)) * (1. / (1. + p2))  );
                // formula returns V - > scale with 1000 to obtain mV; non normalized equaton

            }
        }

    }


    if (!freqs[0]) trf[0] = complex<double> (1.0, 0.);

    return freqs.size();
}

/*!
use a build in transfer function
*/
size_t MTXCal::hardwired(const valarray<double>& newfreq, string chopper, const double& samplefreq, const double& gain, const double& factor,
                         const double& offset, string adc, string tags) {


    return freqs.size();
}

/*!
use a theoretical function
*/
size_t MTXCal::theoretical(const valarray<double>& newfreq, string chopper = "on", const double& samplefreq = 0, const double& gain = 0,
                           const double& factor = 0,  const double& offset = 0, string adc = "", string tags = "") {



    //    valarray<double> old_f(freqs);  // save the frequencies
    //     MTXCal::interpolate(newfreq);


    freqs.resize(newfreq.size());
    trf.resize(newfreq.size());
    freqs = newfreq;


    // interpolate with theoretical transfer function instead of spline
    // outside give values from read in transfer function
    // phase wil jump about 0.1 deg @ 0.1 Hz ... where the measured trf function ends
    // there is no theoretical function for chopper off!
    // interpolation gives bad results for imaginary part!
    complex<double> p1, p2, im(0,1.);

    if (freqs[freqs.size()-1] < 513) {
        cerr << "MTXCal::theoretical -> trf function for mfs coils 05 06 08U KIM" << " " << name << endl;
        for (size_t i = 0; i < freqs.size(); i++) {
            p1 = (freqs[i]/4.) * im;
            p2 = (freqs[i]/8192.) * im;
            trf[i] = 1000. * 0.8 * (  (p1 / ( 1. + p1)) * (1. / (1. + p2))  );
            // formula returns V - > scale with 1000 to obtain mV; non normalized equaton
        }

    } else
        cerr << " no frequency range found" << endl;

    if (!freqs[0])
        trf[0] = complex<double>(1.,0.);

    return freqs.size();
}




/*!
    \fn MTXCal::set_calib(const double& f, const complex<double>& trfn)
 */
size_t MTXCal::set_calib(const valarray<double>& f, const valarray<complex<double> >& trfn) {
    if (f.size() != trfn.size()) {
        cerr << "MTXCal::set_calib -> f and trf are not of the same size, exit!" << endl;
        exit(0);
    }

    freqs.resize(f.size());
    trf.resize(f.size());
    freqs = f;
    trf = trfn;

    return freqs.size();
}





/*!
    \fn ICalbiration::set_calib(const double& f, const double& ampl, const double& ph_deg)
 */
size_t MTXCal::set_calib(const valarray<double>& f, const valarray<double>& ampl, const valarray<double>& ph_deg) {
    /// @todo implement me
    return 0;
}

