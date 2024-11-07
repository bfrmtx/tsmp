/***************************************************************************
                          ats_consts.h  -  description
                             -------------------
    begin                : Mon Aug 20 2001
    copyright            : (C) 2001 by B. Friedrichs
    email                : bfr@metronix.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef ATS_CONSTS_H
#define ATS_CONSTS_H
using namespace std;


/***************************************************************************
*                                                                         *
*                                                                         *
* filter coefficients                                                     *
*                                                                         *
*                                                                         *
*                                                                         *
***************************************************************************/


// filter values for 32x decimation
// 0...470 = 471 coeff
const double coeff32[] = {
    7.636836537E-6,
    -4.584081140E-7,
    -6.420621773E-7,
    -9.585163552E-7,
    -1.406648604E-6,
    -1.998638145E-6,
    -2.738873598E-6,
    -3.644118238E-6,
    -4.722818807E-6,
    -5.995041318E-6,
    -7.472232203E-6,
    -9.176794514E-6,
    -1.112234516E-5,
    -1.333254225E-5,
    -1.582208887E-5,
    -1.861485657E-5,
    -2.173198454E-5,
    -2.516514786E-5,
    -2.898004296E-5,
    -3.317231449E-5,
    -3.775996159E-5,
    -4.276666661E-5,
    -4.820699838E-5,
    -5.410173248E-5,
    -6.046320172E-5,
    -6.730941714E-5,
    -7.465003370E-5,
    -8.249954407E-5,
    -9.086434157E-5,
    -9.975482472E-5,
    -1.091733763E-4,
    -1.191253534E-4,
    -1.296074562E-4,
    -1.406239993E-4,
    -1.521536162E-4,
    -1.641994695E-4,
    -1.767433788E-4,
    -1.897693608E-4,
    -2.032547074E-4,
    -2.171759590E-4,
    -2.315025723E-4,
    -2.462023283E-4,
    -2.612362165E-4,
    -2.765628091E-4,
    -2.921342447E-4,
    -3.078995543E-4,
    -3.238019452E-4,
    -3.397810303E-4,
    -3.557713869E-4,
    -3.717022806E-4,
    -3.875041326E-4,
    -4.030925790E-4,
    -4.183881624E-4,
    -4.333045055E-4,
    -4.477517782E-4,
    -4.616362436E-4,
    -4.748620065E-4,
    -4.873296221E-4,
    -4.989382902E-4,
    -5.095842466E-4,
    -5.191635197E-4,
    -5.275698387E-4,
    -5.346977357E-4,
    -5.404403178E-4,
    -5.446927481E-4,
    -5.473501339E-4,
    -5.483090804E-4,
    -5.474720373E-4,
    -5.447404020E-4,
    -5.400204594E-4,
    -5.332247422E-4,
    -5.242688413E-4,
    -5.130760420E-4,
    -4.995743979E-4,
    -4.837010952E-4,
    -4.653995222E-4,
    -4.446235738E-4,
    -4.213348542E-4,
    -3.955068689E-4,
    -3.671220976E-4,
    -3.361761343E-4,
    -3.026747054E-4,
    -2.666383911E-4,
    -2.280967805E-4,
    -1.870984044E-4,
    -1.437014503E-4,
    -9.798027317E-5,
    -5.002295223E-5,
    6.650805861E-8,
    5.216999048E-5,
    1.061534694E-4,
    1.618682996E-4,
    2.191501347E-4,
    2.778199969E-4,
    3.376837961E-4,
    3.985329065E-4,
    4.601444760E-4,
    5.222817499E-4,
    5.846951505E-4,
    6.471211557E-4,
    7.092871295E-4,
    7.709064997E-4,
    8.316839987E-4,
    8.913148587E-4,
    9.494867270E-4,
    1.005879618E-3,
    1.060168711E-3,
    1.112024067E-3,
    1.161113718E-3,
    1.207103285E-3,
    1.249659522E-3,
    1.288450003E-3,
    1.323146621E-3,
    1.353425440E-3,
    1.378970268E-3,
    1.399473034E-3,
    1.414635734E-3,
    1.424173640E-3,
    1.427815957E-3,
    1.425306899E-3,
    1.416409765E-3,
    1.400906398E-3,
    1.378601223E-3,
    1.349320687E-3,
    1.312917368E-3,
    1.269269289E-3,
    1.218283797E-3,
    1.159896754E-3,
    1.094076440E-3,
    1.020822503E-3,
    9.401694729E-4,
    8.521856064E-4,
    7.569767290E-4,
    6.546832637E-4,
    5.454852368E-4,
    4.295998542E-4,
    3.072833495E-4,
    1.788304330E-4,
    4.457541058E-5,
    -9.510877341E-5,
    -2.398100744E-4,
    -3.890782669E-4,
    -5.424250214E-4,
    -6.993249817E-4,
    -8.592167665E-4,
    -1.021503690E-3,
    -1.185555638E-3,
    -1.350709832E-3,
    -1.516273716E-3,
    -1.681524623E-3,
    -1.845714934E-3,
    -2.008071371E-3,
    -2.167799162E-3,
    -2.324083450E-3,
    -2.476093196E-3,
    -2.622982610E-3,
    -2.763895500E-3,
    -2.897967000E-3,
    -3.024328216E-3,
    -3.142107998E-3,
    -3.250437738E-3,
    -3.348453605E-3,
    -3.435301325E-3,
    -3.510138398E-3,
    -3.572139058E-3,
    -3.620496926E-3,
    -3.654428959E-3,
    -3.673178894E-3,
    -3.676021889E-3,
    -3.662266157E-3,
    -3.631258445E-3,
    -3.582385690E-3,
    -3.515080065E-3,
    -3.428820433E-3,
    -3.323137079E-3,
    -3.197613170E-3,
    -3.051889060E-3,
    -2.885663125E-3,
    -2.698696085E-3,
    -2.490811485E-3,
    -2.261899363E-3,
    -2.011916230E-3,
    -1.740888906E-3,
    -1.448913226E-3,
    -1.136157251E-3,
    -8.028608895E-4,
    -4.493370356E-4,
    -7.597124000E-5,
    3.167778112E-4,
    7.283788660E-4,
    1.158228692E-3,
    1.605652920E-3,
    2.069907434E-3,
    2.550179259E-3,
    3.045589330E-3,
    3.555193420E-3,
    4.077985536E-3,
    4.612899392E-3,
    5.158812763E-3,
    5.714548608E-3,
    6.278880471E-3,
    6.850534598E-3,
    7.428194815E-3,
    8.010505412E-3,
    8.596076712E-3,
    9.183488074E-3,
    9.771293872E-3,
    1.035802663E-2,
    1.094220357E-2,
    1.152233001E-2,
    1.209690556E-2,
    1.266442813E-2,
    1.322340035E-2,
    1.377233339E-2,
    1.430975326E-2,
    1.483420529E-2,
    1.534425967E-2,
    1.583851574E-2,
    1.631560838E-2,
    1.677421125E-2,
    1.721304312E-2,
    1.763087094E-2,
    1.802651596E-2,
    1.839885639E-2,
    1.874683295E-2,
    1.906945131E-2,
    1.936578751E-2,
    1.963498933E-2,
    1.987628129E-2,
    2.008896602E-2,
    2.027242813E-2,
    2.042613519E-2,
    2.054964111E-2,
    2.064258657E-2,
    2.070470105E-2,
    2.073580358E-2,
    2.073580358E-2,
    2.070470105E-2,
    2.064258657E-2,
    2.054964111E-2,
    2.042613519E-2,
    2.027242813E-2,
    2.008896602E-2,
    1.987628129E-2,
    1.963498933E-2,
    1.936578751E-2,
    1.906945131E-2,
    1.874683295E-2,
    1.839885639E-2,
    1.802651596E-2,
    1.763087094E-2,
    1.721304312E-2,
    1.677421125E-2,
    1.631560838E-2,
    1.583851574E-2,
    1.534425967E-2,
    1.483420529E-2,
    1.430975326E-2,
    1.377233339E-2,
    1.322340035E-2,
    1.266442813E-2,
    1.209690556E-2,
    1.152233001E-2,
    1.094220357E-2,
    1.035802663E-2,
    9.771293872E-3,
    9.183488074E-3,
    8.596076712E-3,
    8.010505412E-3,
    7.428194815E-3,
    6.850534598E-3,
    6.278880471E-3,
    5.714548608E-3,
    5.158812763E-3,
    4.612899392E-3,
    4.077985536E-3,
    3.555193420E-3,
    3.045589330E-3,
    2.550179259E-3,
    2.069907434E-3,
    1.605652920E-3,
    1.158228692E-3,
    7.283788660E-4,
    3.167778112E-4,
    -7.597124000E-5,
    -4.493370356E-4,
    -8.028608895E-4,
    -1.136157251E-3,
    -1.448913226E-3,
    -1.740888906E-3,
    -2.011916230E-3,
    -2.261899363E-3,
    -2.490811485E-3,
    -2.698696085E-3,
    -2.885663125E-3,
    -3.051889060E-3,
    -3.197613170E-3,
    -3.323137079E-3,
    -3.428820433E-3,
    -3.515080065E-3,
    -3.582385690E-3,
    -3.631258445E-3,
    -3.662266157E-3,
    -3.676021889E-3,
    -3.673178894E-3,
    -3.654428959E-3,
    -3.620496926E-3,
    -3.572139058E-3,
    -3.510138398E-3,
    -3.435301325E-3,
    -3.348453605E-3,
    -3.250437738E-3,
    -3.142107998E-3,
    -3.024328216E-3,
    -2.897967000E-3,
    -2.763895500E-3,
    -2.622982610E-3,
    -2.476093196E-3,
    -2.324083450E-3,
    -2.167799162E-3,
    -2.008071371E-3,
    -1.845714934E-3,
    -1.681524623E-3,
    -1.516273716E-3,
    -1.350709832E-3,
    -1.185555638E-3,
    -1.021503690E-3,
    -8.592167665E-4,
    -6.993249817E-4,
    -5.424250214E-4,
    -3.890782669E-4,
    -2.398100744E-4,
    -9.510877341E-5,
    4.457541058E-5,
    1.788304330E-4,
    3.072833495E-4,
    4.295998542E-4,
    5.454852368E-4,
    6.546832637E-4,
    7.569767290E-4,
    8.521856064E-4,
    9.401694729E-4,
    1.020822503E-3,
    1.094076440E-3,
    1.159896754E-3,
    1.218283797E-3,
    1.269269289E-3,
    1.312917368E-3,
    1.349320687E-3,
    1.378601223E-3,
    1.400906398E-3,
    1.416409765E-3,
    1.425306899E-3,
    1.427815957E-3,
    1.424173640E-3,
    1.414635734E-3,
    1.399473034E-3,
    1.378970268E-3,
    1.353425440E-3,
    1.323146621E-3,
    1.288450003E-3,
    1.249659522E-3,
    1.207103285E-3,
    1.161113718E-3,
    1.112024067E-3,
    1.060168711E-3,
    1.005879618E-3,
    9.494867270E-4,
    8.913148587E-4,
    8.316839987E-4,
    7.709064997E-4,
    7.092871295E-4,
    6.471211557E-4,
    5.846951505E-4,
    5.222817499E-4,
    4.601444760E-4,
    3.985329065E-4,
    3.376837961E-4,
    2.778199969E-4,
    2.191501347E-4,
    1.618682996E-4,
    1.061534694E-4,
    5.216999048E-5,
    6.650805861E-8,
    -5.002295223E-5,
    -9.798027317E-5,
    -1.437014503E-4,
    -1.870984044E-4,
    -2.280967805E-4,
    -2.666383911E-4,
    -3.026747054E-4,
    -3.361761343E-4,
    -3.671220976E-4,
    -3.955068689E-4,
    -4.213348542E-4,
    -4.446235738E-4,
    -4.653995222E-4,
    -4.837010952E-4,
    -4.995743979E-4,
    -5.130760420E-4,
    -5.242688413E-4,
    -5.332247422E-4,
    -5.400204594E-4,
    -5.447404020E-4,
    -5.474720373E-4,
    -5.483090804E-4,
    -5.473501339E-4,
    -5.446927481E-4,
    -5.404403178E-4,
    -5.346977357E-4,
    -5.275698387E-4,
    -5.191635197E-4,
    -5.095842466E-4,
    -4.989382902E-4,
    -4.873296221E-4,
    -4.748620065E-4,
    -4.616362436E-4,
    -4.477517782E-4,
    -4.333045055E-4,
    -4.183881624E-4,
    -4.030925790E-4,
    -3.875041326E-4,
    -3.717022806E-4,
    -3.557713869E-4,
    -3.397810303E-4,
    -3.238019452E-4,
    -3.078995543E-4,
    -2.921342447E-4,
    -2.765628091E-4,
    -2.612362165E-4,
    -2.462023283E-4,
    -2.315025723E-4,
    -2.171759590E-4,
    -2.032547074E-4,
    -1.897693608E-4,
    -1.767433788E-4,
    -1.641994695E-4,
    -1.521536162E-4,
    -1.406239993E-4,
    -1.296074562E-4,
    -1.191253534E-4,
    -1.091733763E-4,
    -9.975482472E-5,
    -9.086434157E-5,
    -8.249954407E-5,
    -7.465003370E-5,
    -6.730941714E-5,
    -6.046320172E-5,
    -5.410173248E-5,
    -4.820699838E-5,
    -4.276666661E-5,
    -3.775996159E-5,
    -3.317231449E-5,
    -2.898004296E-5,
    -2.516514786E-5,
    -2.173198454E-5,
    -1.861485657E-5,
    -1.582208887E-5,
    -1.333254225E-5,
    -1.112234516E-5,
    -9.176794514E-6,
    -7.472232203E-6,
    -5.995041318E-6,
    -4.722818807E-6,
    -3.644118238E-6,
    -2.738873598E-6,
    -1.998638145E-6,
    -1.406648604E-6,
    -9.585163552E-7,
    -6.420621773E-7,
    -4.584081140E-7
};


//  filter values for 4x decimation       0 .. 70 = 71 **************************************
const double coeff4[] = {
    6.854913864E-6,
    2.933716833E-5,
    8.267682092E-5,
    1.822009396E-4,
    3.346982688E-4,
    5.253549867E-4,
    7.050447753E-4,
    7.848495715E-4,
    6.456473839E-4,
    1.680573886E-4,
    -7.184734519E-4,
    -1.978042517E-3,
    -3.421991913E-3,
    -4.698704960E-3,
    -5.340810348E-3,
    -4.877338474E-3,
    -2.996021341E-3,
    2.840936053E-4,
    4.488289189E-3,
    8.671796886E-3,
    1.155862537E-2,
    1.183726918E-2,
    8.564797604E-3,
    1.581379992E-3,
    -8.191622918E-3,
    -1.869064668E-2,
    -2.696017386E-2,
    -2.970906692E-2,
    -2.407566432E-2,
    -8.426506187E-3,
    1.702308969E-2,
    4.995732286E-2,
    8.618845263E-2,
    1.203573989E-1,
    1.469705098E-1,
    1.615376929E-1,
    1.615376929E-1,
    1.469705098E-1,
    1.203573989E-1,
    8.618845263E-2,
    4.995732286E-2,
    1.702308969E-2,
    -8.426506187E-3,
    -2.407566432E-2,
    -2.970906692E-2,
    -2.696017386E-2,
    -1.869064668E-2,
    -8.191622918E-3,
    1.581379992E-3,
    8.564797604E-3,
    1.183726918E-2,
    1.155862537E-2,
    8.671796886E-3,
    4.488289189E-3,
    2.840936053E-4,
    -2.996021341E-3,
    -4.877338474E-3,
    -5.340810348E-3,
    -4.698704960E-3,
    -3.421991913E-3,
    -1.978042517E-3,
    -7.184734519E-4,
    1.680573886E-4,
    6.456473839E-4,
    7.848495715E-4,
    7.050447753E-4,
    5.253549867E-4,
    3.346982688E-4,
    1.822009396E-4,
    8.267682092E-5,
    2.933716833E-5
};

// 501 (0-500) coeff for filtering 25x
// not checked!
const double coeff25[] = {

    6.033421866E-6,
    2.676325810E-6,
    3.222534925E-6,
    3.827166873E-6,
    4.453934884E-6,
    5.124141231E-6,
    5.800775889E-6,
    6.505914843E-6,
    7.197343949E-6,
    7.884980426E-6,
    8.516606921E-6,
    9.114891746E-6,
    9.625235634E-6,
    1.003648676E-5,
    1.031824274E-5,
    1.044943685E-5,
    1.038617157E-5,
    1.011188154E-5,
    9.572239634E-6,
    8.752664909E-6,
    7.595561728E-6,
    6.081250392E-6,
    4.150454235E-6,
    1.789225480E-6,
    -1.065808633E-6,
    -4.430573845E-6,
    -8.359587903E-6,
    -1.287124785E-5,
    -1.801261028E-5,
    -2.379862795E-5,
    -3.027002047E-5,
    -3.743172185E-5,
    -4.531537798E-5,
    -5.391767489E-5,
    -6.325699327E-5,
    -7.331555427E-5,
    -8.410042835E-5,
    -9.557776617E-5,
    -1.077378938E-4,
    -1.205313165E-4,
    -1.339280209E-4,
    -1.478624139E-4,
    -1.622839197E-4,
    -1.771065961E-4,
    -1.922603218E-4,
    -2.076404218E-4,
    -2.231557200E-4,
    -2.386825730E-4,
    -2.541121312E-4,
    -2.693021876E-4,
    -2.841284830E-4,
    -2.984328933E-4,
    -3.120773912E-4,
    -3.248937697E-4,
    -3.367334631E-4,
    -3.474213494E-4,
    -3.568048631E-4,
    -3.647079355E-4,
    -3.709788142E-4,
    -3.754470898E-4,
    -3.779688976E-4,
    -3.783859922E-4,
    -3.765705658E-4,
    -3.723827379E-4,
    -3.657183535E-4,
    -3.564650030E-4,
    -3.445484294E-4,
    -3.298909654E-4,
    -3.124562161E-4,
    -2.922080968E-4,
    -2.691542843E-4,
    -2.433066376E-4,
    -2.147221058E-4,
    -1.834659442E-4,
    -1.496490141E-4,
    -1.133921552E-4,
    -7.486367851E-5,
    -3.424196288E-5,
    8.247412141E-6,
    5.236884961E-5,
    9.784101759E-5,
    1.443740374E-4,
    1.916347771E-4,
    2.392845909E-4,
    2.869453901E-4,
    3.342366967E-4,
    3.807448877E-4,
    4.260584841E-4,
    4.697388935E-4,
    5.113564642E-4,
    5.504613603E-4,
    5.866202820E-4,
    6.193877286E-4,
    6.483427819E-4,
    6.730610289E-4,
    6.931517916E-4,
    7.082294153E-4,
    7.179508008E-4,
    7.219877825E-4,
    7.200626478E-4,
    7.119220394E-4,
    6.973709513E-4,
    6.762470510E-4,
    6.484537892E-4,
    6.139343821E-4,
    5.727036117E-4,
    5.248221025E-4,
    4.704266009E-4,
    4.097031683E-4,
    3.429168574E-4,
    2.703839955E-4,
    1.924998527E-4,
    1.097109919E-4,
    2.254020501E-5,
    -6.844100913E-5,
    -1.625896783E-4,
    -2.592198443E-4,
    -3.575815770E-4,
    -4.568898156E-4,
    -5.563055942E-4,
    -6.549657827E-4,
    -7.519653048E-4,
    -8.463883441E-4,
    -9.372913391E-4,
    -1.023735145E-3,
    -1.104769398E-3,
    -1.179464486E-3,
    -1.246897229E-3,
    -1.306183083E-3,
    -1.356462282E-3,
    -1.396931670E-3,
    -1.426831634E-3,
    -1.445477047E-3,
    -1.452244350E-3,
    -1.446601729E-3,
    -1.428095363E-3,
    -1.396378986E-3,
    -1.351199147E-3,
    -1.292422892E-3,
    -1.220022383E-3,
    -1.134100498E-3,
    -1.034874080E-3,
    -9.226971497E-4,
    -7.980428757E-4,
    -6.615237753E-4,
    -5.138725193E-4,
    -3.559587663E-4,
    -1.887682440E-4,
    -1.341699361E-5,
    1.688719162E-4,
    3.567569569E-4,
    5.488033375E-4,
    7.434744001E-4,
    9.391588386E-4,
    1.134164614E-3,
    1.326748102E-3,
    1.515110304E-3,
    1.697427268E-3,
    1.871848875E-3,
    2.036529976E-3,
    2.189630828E-3,
    2.329349171E-3,
    2.453921567E-3,
    2.561655869E-3,
    2.650932956E-3,
    2.720239007E-3,
    2.768167089E-3,
    2.793448797E-3,
    2.794954825E-3,
    2.771725254E-3,
    2.722969136E-3,
    2.648092202E-3,
    2.546695009E-3,
    2.418597710E-3,
    2.263835998E-3,
    2.082682612E-3,
    1.875640462E-3,
    1.643460332E-3,
    1.387131025E-3,
    1.107892930E-3,
    8.072245142E-4,
    4.868517704E-4,
    1.487310630E-4,
    -2.049460614E-4,
    -5.717734569E-4,
    -9.491303179E-4,
    -1.334205093E-3,
    -1.723999814E-3,
    -2.115357462E-3,
    -2.504970297E-3,
    -2.889410113E-3,
    -3.265140259E-3,
    -3.628548543E-3,
    -3.975962028E-3,
    -4.303682165E-3,
    -4.608001668E-3,
    -4.885241117E-3,
    -5.131767357E-3,
    -5.344030530E-3,
    -5.518583172E-3,
    -5.652117097E-3,
    -5.741482169E-3,
    -5.783721971E-3,
    -5.776091629E-3,
    -5.716091302E-3,
    -5.601482126E-3,
    -5.430316754E-3,
    -5.200952413E-3,
    -4.912077878E-3,
    -4.562722931E-3,
    -4.152280771E-3,
    -3.680513387E-3,
    -3.147568876E-3,
    -2.553982101E-3,
    -1.900686349E-3,
    -1.189009041E-3,
    -4.206772836E-4,
    4.021914379E-4,
    1.277090101E-3,
    2.201136628E-3,
    3.171080776E-3,
    4.183323892E-3,
    5.233931659E-3,
    6.318658779E-3,
    7.432967122E-3,
    8.572054943E-3,
    9.730879853E-3,
    1.090419199E-2,
    1.208656106E-2,
    1.327241254E-2,
    1.445605819E-2,
    1.563173433E-2,
    1.679363471E-2,
    1.793595004E-2,
    1.905290209E-2,
    2.013878354E-2,
    2.118799213E-2,
    2.219506970E-2,
    2.315473541E-2,
    2.406192299E-2,
    2.491181170E-2,
    2.569986079E-2,
    2.642183753E-2,
    2.707384757E-2,
    2.765235921E-2,
    2.815422883E-2,
    2.857672063E-2,
    2.891752642E-2,
    2.917478012E-2,
    2.934707129E-2,
    2.943345408E-2,
    2.943345408E-2,
    2.934707129E-2,
    2.917478012E-2,
    2.891752642E-2,
    2.857672063E-2,
    2.815422883E-2,
    2.765235921E-2,
    2.707384757E-2,
    2.642183753E-2,
    2.569986079E-2,
    2.491181170E-2,
    2.406192299E-2,
    2.315473541E-2,
    2.219506970E-2,
    2.118799213E-2,
    2.013878354E-2,
    1.905290209E-2,
    1.793595004E-2,
    1.679363471E-2,
    1.563173433E-2,
    1.445605819E-2,
    1.327241254E-2,
    1.208656106E-2,
    1.090419199E-2,
    9.730879853E-3,
    8.572054943E-3,
    7.432967122E-3,
    6.318658779E-3,
    5.233931659E-3,
    4.183323892E-3,
    3.171080776E-3,
    2.201136628E-3,
    1.277090101E-3,
    4.021914379E-4,
    -4.206772836E-4,
    -1.189009041E-3,
    -1.900686349E-3,
    -2.553982101E-3,
    -3.147568876E-3,
    -3.680513387E-3,
    -4.152280771E-3,
    -4.562722931E-3,
    -4.912077878E-3,
    -5.200952413E-3,
    -5.430316754E-3,
    -5.601482126E-3,
    -5.716091302E-3,
    -5.776091629E-3,
    -5.783721971E-3,
    -5.741482169E-3,
    -5.652117097E-3,
    -5.518583172E-3,
    -5.344030530E-3,
    -5.131767357E-3,
    -4.885241117E-3,
    -4.608001668E-3,
    -4.303682165E-3,
    -3.975962028E-3,
    -3.628548543E-3,
    -3.265140259E-3,
    -2.889410113E-3,
    -2.504970297E-3,
    -2.115357462E-3,
    -1.723999814E-3,
    -1.334205093E-3,
    -9.491303179E-4,
    -5.717734569E-4,
    -2.049460614E-4,
    1.487310630E-4,
    4.868517704E-4,
    8.072245142E-4,
    1.107892930E-3,
    1.387131025E-3,
    1.643460332E-3,
    1.875640462E-3,
    2.082682612E-3,
    2.263835998E-3,
    2.418597710E-3,
    2.546695009E-3,
    2.648092202E-3,
    2.722969136E-3,
    2.771725254E-3,
    2.794954825E-3,
    2.793448797E-3,
    2.768167089E-3,
    2.720239007E-3,
    2.650932956E-3,
    2.561655869E-3,
    2.453921567E-3,
    2.329349171E-3,
    2.189630828E-3,
    2.036529976E-3,
    1.871848875E-3,
    1.697427268E-3,
    1.515110304E-3,
    1.326748102E-3,
    1.134164614E-3,
    9.391588386E-4,
    7.434744001E-4,
    5.488033375E-4,
    3.567569569E-4,
    1.688719162E-4,
    -1.341699361E-5,
    -1.887682440E-4,
    -3.559587663E-4,
    -5.138725193E-4,
    -6.615237753E-4,
    -7.980428757E-4,
    -9.226971497E-4,
    -1.034874080E-3,
    -1.134100498E-3,
    -1.220022383E-3,
    -1.292422892E-3,
    -1.351199147E-3,
    -1.396378986E-3,
    -1.428095363E-3,
    -1.446601729E-3,
    -1.452244350E-3,
    -1.445477047E-3,
    -1.426831634E-3,
    -1.396931670E-3,
    -1.356462282E-3,
    -1.306183083E-3,
    -1.246897229E-3,
    -1.179464486E-3,
    -1.104769398E-3,
    -1.023735145E-3,
    -9.372913391E-4,
    -8.463883441E-4,
    -7.519653048E-4,
    -6.549657827E-4,
    -5.563055942E-4,
    -4.568898156E-4,
    -3.575815770E-4,
    -2.592198443E-4,
    -1.625896783E-4,
    -6.844100913E-5,
    2.254020501E-5,
    1.097109919E-4,
    1.924998527E-4,
    2.703839955E-4,
    3.429168574E-4,
    4.097031683E-4,
    4.704266009E-4,
    5.248221025E-4,
    5.727036117E-4,
    6.139343821E-4,
    6.484537892E-4,
    6.762470510E-4,
    6.973709513E-4,
    7.119220394E-4,
    7.200626478E-4,
    7.219877825E-4,
    7.179508008E-4,
    7.082294153E-4,
    6.931517916E-4,
    6.730610289E-4,
    6.483427819E-4,
    6.193877286E-4,
    5.866202820E-4,
    5.504613603E-4,
    5.113564642E-4,
    4.697388935E-4,
    4.260584841E-4,
    3.807448877E-4,
    3.342366967E-4,
    2.869453901E-4,
    2.392845909E-4,
    1.916347771E-4,
    1.443740374E-4,
    9.784101759E-5,
    5.236884961E-5,
    8.247412141E-6,
    -3.424196288E-5,
    -7.486367851E-5,
    -1.133921552E-4,
    -1.496490141E-4,
    -1.834659442E-4,
    -2.147221058E-4,
    -2.433066376E-4,
    -2.691542843E-4,
    -2.922080968E-4,
    -3.124562161E-4,
    -3.298909654E-4,
    -3.445484294E-4,
    -3.564650030E-4,
    -3.657183535E-4,
    -3.723827379E-4,
    -3.765705658E-4,
    -3.783859922E-4,
    -3.779688976E-4,
    -3.754470898E-4,
    -3.709788142E-4,
    -3.647079355E-4,
    -3.568048631E-4,
    -3.474213494E-4,
    -3.367334631E-4,
    -3.248937697E-4,
    -3.120773912E-4,
    -2.984328933E-4,
    -2.841284830E-4,
    -2.693021876E-4,
    -2.541121312E-4,
    -2.386825730E-4,
    -2.231557200E-4,
    -2.076404218E-4,
    -1.922603218E-4,
    -1.771065961E-4,
    -1.622839197E-4,
    -1.478624139E-4,
    -1.339280209E-4,
    -1.205313165E-4,
    -1.077378938E-4,
    -9.557776617E-5,
    -8.410042835E-5,
    -7.331555427E-5,
    -6.325699327E-5,
    -5.391767489E-5,
    -4.531537798E-5,
    -3.743172185E-5,
    -3.027002047E-5,
    -2.379862795E-5,
    -1.801261028E-5,
    -1.287124785E-5,
    -8.359587903E-6,
    -4.430573845E-6,
    -1.065808633E-6,
    1.789225480E-6,
    4.150454235E-6,
    6.081250392E-6,
    7.595561728E-6,
    8.752664909E-6,
    9.572239634E-6,
    1.011188154E-5,
    1.038617157E-5,
    1.044943685E-5,
    1.031824274E-5,
    1.003648676E-5,
    9.625235634E-6,
    9.114891746E-6,
    8.516606921E-6,
    7.884980426E-6,
    7.197343949E-6,
    6.505914843E-6,
    5.800775889E-6,
    5.124141231E-6,
    4.453934884E-6,
    3.827166873E-6,
    3.222534925E-6,
    2.676325810E-6


};

// 35 coeff
const double coeff2[] = {
  -1.372521126E-4,
  -6.667804822E-4,
  -1.186883027E-3,
  -3.696243078E-4,
  2.433044326E-3,
  4.527364971E-3,
  1.082752294E-3,
  -7.812171581E-3,
  -1.213328444E-2,
  -4.843352987E-4,
  2.144382608E-2,
  2.653128881E-2,
  -6.473890927E-3,
  -5.546122574E-2,
  -5.601629060E-2,
  4.114742442E-2,
  2.065781657E-1,
  3.365394460E-1,
  3.365394460E-1,
  2.065781657E-1,
  4.114742442E-2,
  -5.601629060E-2,
  -5.546122574E-2,
  -6.473890927E-3,
  2.653128881E-2,
  2.144382608E-2,
  -4.843352987E-4,
  -1.213328444E-2,
  -7.812171581E-3,
  1.082752294E-3,
  4.527364971E-3,
  2.433044326E-3,
  -3.696243078E-4,
  -1.186883027E-3,
  -6.667804822E-4
};

// 159
const double coeff8[] = {
  1.274638071E-5,
  1.871833207E-5,
  3.051448966E-5,
  4.508428102E-5,
  6.158830100E-5,
  7.858234792E-5,
  9.395673097E-5,
  1.050224364E-4,
  1.087241754E-4,
  1.018846125E-4,
  8.161889406E-5,
  4.577110201E-5,
  -6.565664053E-6,
  -7.462226237E-5,
  -1.555564299E-4,
  -2.442676894E-4,
  -3.333774120E-4,
  -4.135467518E-4,
  -4.740400810E-4,
  -5.036386902E-4,
  -4.917498406E-4,
  -4.297519242E-4,
  -3.123499630E-4,
  -1.389044065E-4,
  8.554531293E-5,
  3.496159862E-4,
  6.355282534E-4,
  9.196709247E-4,
  1.173951543E-3,
  1.367818261E-3,
  1.470963341E-3,
  1.456425398E-3,
  1.303939654E-3,
  1.003123800E-3,
  5.562490955E-4,
  -1.986065669E-5,
  -6.930342056E-4,
  -1.416414222E-3,
  -2.130629748E-3,
  -2.767585135E-3,
  -3.255632794E-3,
  -3.525914544E-3,
  -3.519325575E-3,
  -3.193601619E-3,
  -2.529775410E-3,
  -1.537385494E-3,
  -2.576911868E-4,
  1.235589439E-3,
  2.838487276E-3,
  4.422025666E-3,
  5.840473903E-3,
  6.942157624E-3,
  7.582235943E-3,
  7.636471028E-3,
  7.014961223E-3,
  5.674575302E-3,
  3.628945214E-3,
  9.548541383E-4,
  -2.205820901E-3,
  -5.649263137E-3,
  -9.118130628E-3,
  -1.231507218E-2,
  -1.492064668E-2,
  -1.661462901E-2,
  -1.709924220E-2,
  -1.612268602E-2,
  -1.350111910E-2,
  -9.137342056E-3,
  -3.034517438E-3,
  4.696367013E-3,
  1.383622388E-2,
  2.406513571E-2,
  3.497700849E-2,
  4.610082290E-2,
  5.692720913E-2,
  6.693858403E-2,
  7.564080472E-2,
  8.259409732E-2,
  8.744105514E-2,
  8.992966248E-2,
  8.992966248E-2,
  8.744105514E-2,
  8.259409732E-2,
  7.564080472E-2,
  6.693858403E-2,
  5.692720913E-2,
  4.610082290E-2,
  3.497700849E-2,
  2.406513571E-2,
  1.383622388E-2,
  4.696367013E-3,
  -3.034517438E-3,
  -9.137342056E-3,
  -1.350111910E-2,
  -1.612268602E-2,
  -1.709924220E-2,
  -1.661462901E-2,
  -1.492064668E-2,
  -1.231507218E-2,
  -9.118130628E-3,
  -5.649263137E-3,
  -2.205820901E-3,
  9.548541383E-4,
  3.628945214E-3,
  5.674575302E-3,
  7.014961223E-3,
  7.636471028E-3,
  7.582235943E-3,
  6.942157624E-3,
  5.840473903E-3,
  4.422025666E-3,
  2.838487276E-3,
  1.235589439E-3,
  -2.576911868E-4,
  -1.537385494E-3,
  -2.529775410E-3,
  -3.193601619E-3,
  -3.519325575E-3,
  -3.525914544E-3,
  -3.255632794E-3,
  -2.767585135E-3,
  -2.130629748E-3,
  -1.416414222E-3,
  -6.930342056E-4,
  -1.986065669E-5,
  5.562490955E-4,
  1.003123800E-3,
  1.303939654E-3,
  1.456425398E-3,
  1.470963341E-3,
  1.367818261E-3,
  1.173951543E-3,
  9.196709247E-4,
  6.355282534E-4,
  3.496159862E-4,
  8.554531293E-5,
  -1.389044065E-4,
  -3.123499630E-4,
  -4.297519242E-4,
  -4.917498406E-4,
  -5.036386902E-4,
  -4.740400810E-4,
  -4.135467518E-4,
  -3.333774120E-4,
  -2.442676894E-4,
  -1.555564299E-4,
  -7.462226237E-5,
  -6.565664053E-6,
  4.577110201E-5,
  8.161889406E-5,
  1.018846125E-4,
  1.087241754E-4,
  1.050224364E-4,
  9.395673097E-5,
  7.858234792E-5,
  6.158830100E-5,
  4.508428102E-5,
  3.051448966E-5,
  1.871833207E-5

};


#endif