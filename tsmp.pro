include (../../../../compile.pri)
include (../../../../bins.pri)

CONFIG +=console

INCLUDEPATH += ../fftReal


#TARGET = tsmp
QT -= gui
QT += core

SOURCES += main.cpp \
 adu06cal.cpp \
 atmfile.cpp \
 atsheader72.cpp \
 atsheader75.cpp \
 coherency.cpp \
 efp05cal.cpp \
 emi_bfcal.cpp \
 emimt24header.cpp \
 emisyscal.cpp \
 fgs01cal.cpp \
 fgs02cal.cpp \
 fubheader.cpp \
 icalibration.cpp \
 itsheader.cpp \
 ltscal.cpp \
 mtstd.cpp \
 mtxcal.cpp \
 spectra.cpp \
 tsdata.cpp \
 tsheader.cpp \
 date_utils.cpp \
 koordi.cpp \
 my_valarray_double.cpp \
 string_utils.cpp \
 atsfilename.cpp \
 mfs07cal.cpp \
 adu07cal.cpp \
 eDateTime.cpp \
 plaincal.cpp \
 atsheader80.cpp

HEADERS += adu06cal.h \
 allinclude.h \
 atmfile.h \
 ats_consts.h \
 atsheader72.h \
 atsheader75.h \
 coherency.h \
 efp05cal.h \
 emi_bfcal.h \
 emimt24header.h \
 emisyscal.h \
 fgs01cal.h \
 fgs02cal.h \
 fubheader.h \
 icalibration.h \
 itsheader.h \
 ltscal.h \
 mtstd.h \
 mtxcal.h \
 spectra.h \
 tsdata.h \
 tsheader.h \
 bfr_iter.h \
 bt_utils.h \
 date_utils.h \
 koordi.h \
 mysql_simple_vector.h \
 my_valarray_double.h \
 my_valarray.h \
 myvector.h \
 string_utils.h \
 atsfilename.h \
 mfs07cal.h \
 adu07cal.h \
 mtx_iter.h \
 eDateTime.h \
 c_slice_iter.h \
 slice_iter.h \
 plaincal.h \
 atsheader80.h
