// Tell emacs that this is -*-c++-*- mode

/*
 *  SBN.h
 *
 *
 *  Created by Jason Immerman on 5/29/09.
 *  Altered significantly on 4/7/10.
 *  Copyright 2009 Bowdoin College. All rights reserved.
 *
 */

#ifndef _SBN_H
#define _SBN_H
#include <iostream>
#include <iomanip>
#include "nr3.h"
#include <math.h>
 // using namespace std;


class SBN {
public:

    //===================================================================================================
    //                 Constructor (initialize values in all arrays)
    //
    //		NOTES ON ARRAY LENGTH:	CONFORMAL FACTOR (for maxAreal = 25)
    //					  1e5 for resolution of 3e-2 minimum
    //					  1e6 for resolution of 6e-4 minimum
    //
    //					LAPSE (for maxAreal = 25)
    //					  1e5 for resolution of 5.3e-3 minimum
    //					  1e6 for resolution of 2e-3 minimum
    //					  1e7 for resolution of 3e-4 minimum
    //
    //===================================================================================================
    SBN(Doub maxAreal, int arrayLength) {
        m = 1.0;
        Doub minAreal = 3.0 * m / 2.0 + 1e-15;
        length = arrayLength;
        arealMax = maxAreal;

        arealArray_p = new VecDoub(arrayLength);
        isotropicArray_p = new VecDoub(arrayLength);
        conFactorArray_p = new VecDoub(arrayLength);
        shiftArray_p = new VecDoub(arrayLength);
        lapseArray_p = new VecDoub(arrayLength);
        dConFactorArray_p = new VecDoub(arrayLength);
        dLapseArray_p = new VecDoub(arrayLength);
        ddConFactorArray_p = new VecDoub(arrayLength);
        ddLapseArray_p = new VecDoub(arrayLength);
        dshiftArray_p = new VecDoub(arrayLength);
        ddshiftArray_p = new VecDoub(arrayLength);

        VecDoub& arealArray = *arealArray_p, & isotropicArray = *isotropicArray_p,
            & conFactorArray = *conFactorArray_p, & shiftArray = *shiftArray_p,
            & lapseArray = *lapseArray_p;

        for (int i = 0; i < length; i++) {
            arealArray[i] = minAreal + i * ((arealMax - minAreal) / (length - 1.0));
            isotropicArray[i] = (2.0 * arealArray[i] + m + sqrt(4.0 * arealArray[i] * arealArray[i] + 4.0 * m * arealArray[i] + 3.0 * m * m)) / 4.0 * pow(((4.0 + 3.0 * sqrt(2.0)) * (2.0 * arealArray[i] - 3.0 * m)) / (8.0 * arealArray[i] + 6.0 * m + 3.0 * sqrt(8.0 * arealArray[i] * arealArray[i] + 8.0 * m * arealArray[i] + 6.0 * m * m)), 1.0 / sqrt(2.0));
            conFactorArray[i] = sqrt(arealArray[i] / isotropicArray[i]);
            lapseArray[i] = sqrt(1.0 - 2.0 * m / arealArray[i] + 27.0 * pow(m, 4.0) / (16.0 * pow(arealArray[i], 4.0)));
            shiftArray[i] = 3.0 * sqrt(3.0) * pow(m, 2.0) / 4.0 * isotropicArray[i] / pow(arealArray[i], 3.0);

            //if(lapseArray[i] =! lapseArray[i]) cout << lapseArray[i] << endl;

            if (isotropicArray[i] < .02) conFactorArray[i] = sqrt(3.0 * m / (2.0 * isotropicArray[i]));

            //cout << isotropicArray[i]  << endl;
        }
        firstDeriv(conFactorArray_p, dConFactorArray_p);
        firstDeriv(lapseArray_p, dLapseArray_p);
        firstDeriv(shiftArray_p, dshiftArray_p);
        secondDeriv(conFactorArray_p, ddConFactorArray_p);
        secondDeriv(lapseArray_p, ddLapseArray_p);
        secondDeriv(shiftArray_p, ddshiftArray_p);
    }




    //===================================================================================================
    //                  Destructor
    //===================================================================================================
    ~SBN() {
        delete arealArray_p;
        delete isotropicArray_p;
        delete conFactorArray_p;
        delete shiftArray_p;
        delete lapseArray_p;
        delete dConFactorArray_p;
        delete dshiftArray_p;
        delete dLapseArray_p;
        delete ddConFactorArray_p;
        delete ddshiftArray_p;
        delete ddLapseArray_p;
    }




    //===================================================================================================
    //                  Reset array values if mass has changed
    //===================================================================================================
    void setMass() {
        Doub minAreal = 3.0 * m / 2.0 + 1e-15;

        VecDoub& arealArray = *arealArray_p, & isotropicArray = *isotropicArray_p, & conFactorArray = *conFactorArray_p, & shiftArray = *shiftArray_p, & lapseArray = *lapseArray_p;

        for (int i = 0; i < length; i++) {
            arealArray[i] = minAreal + i * ((arealMax - minAreal) / (length - 1.0));
            isotropicArray[i] = (2.0 * arealArray[i] + m + sqrt(4.0 * arealArray[i] * arealArray[i] + 4.0 * m * arealArray[i] + 3.0 * m * m)) / 4.0 * pow(((4.0 + 3.0 * sqrt(2.0)) * (2.0 * arealArray[i] - 3.0 * m)) / (8.0 * arealArray[i] + 6.0 * m + 3.0 * sqrt(8.0 * arealArray[i] * arealArray[i] + 8.0 * m * arealArray[i] + 6.0 * m * m)), 1.0 / sqrt(2.0));
            conFactorArray[i] = sqrt(arealArray[i] / isotropicArray[i]);
            lapseArray[i] = sqrt(1.0 - 2.0 * m / arealArray[i] + 27.0 * pow(m, 4.0) / (16.0 * pow(arealArray[i], 4.0)));
            shiftArray[i] = 3.0 * sqrt(3.0) * pow(m, 2.0) / 4.0 * isotropicArray[i] / pow(arealArray[i], 3.0);


            if (isotropicArray[i] < .02) conFactorArray[i] = sqrt(3.0 * m / (2.0 * isotropicArray[i]));

            //cout << isotropicArray[i] << "   " << conFactorArray[i] << endl;
        }

        firstDeriv(conFactorArray_p, dConFactorArray_p);
        firstDeriv(lapseArray_p, dLapseArray_p);
        firstDeriv(shiftArray_p, dshiftArray_p);
        secondDeriv(conFactorArray_p, ddConFactorArray_p);
        secondDeriv(lapseArray_p, ddLapseArray_p);
        secondDeriv(shiftArray_p, ddshiftArray_p);
    }




    //===================================================================================================
    //                  Externally Called routines to find value in arrays
    //===================================================================================================
    Doub areal(Doub inputIso, Doub mass) {
        Doub value = searchPolint(inputIso, mass, arealArray_p);
        return value;
    }

    Doub conFactor(Doub inputIso, Doub mass) {
        Doub value = searchRatint(inputIso, mass, conFactorArray_p, 1);
        return value;
    }

    Doub shift(Doub inputIso, Doub mass) {
        Doub value = searchPolint(inputIso, mass, shiftArray_p);
        return value;
    }

    Doub lapse(Doub inputIso, Doub mass) {
        Doub value = searchPolint(inputIso, mass, lapseArray_p);
        return value;
    }

    Doub dConFactor(Doub inputIso, Doub mass) {
        Doub value = searchRatint(inputIso, mass, dConFactorArray_p, 2);
        return value;
    }

    Doub dLapse(Doub inputIso, Doub mass) {
        Doub value = searchPolint(inputIso, mass, dLapseArray_p);
        return value;
    }

    Doub ddLapse(Doub inputIso, Doub mass) {
        Doub value = searchRatint(inputIso, mass, ddLapseArray_p, 0);
        return value;
    }
    Doub dShift(Doub inputIso, Doub mass) {
        Doub value = searchPolint(inputIso, mass, dshiftArray_p);
        return value;
    }

    Doub ddShift(Doub inputIso, Doub mass) {
        Doub value = searchPolint(inputIso, mass, ddshiftArray_p);
        return value;
    }

    Doub ddConFactor(Doub inputIso, Doub mass) {
        Doub value = searchRatint(inputIso, mass, ddConFactorArray_p, 3);
        return value;
    }




    //===================================================================================================
    // Hunt for value in Array (rational versus polynomial interpolation...)
    //===================================================================================================
    Doub searchPolint(Doub inputIso, Doub mass, VecDoub* array_p) {
        if (mass != m) {
            m = mass;
            setMass();
        }

        VecDoub& array = *array_p;
        VecDoub& isotropicArray = *isotropicArray_p;

        int index = 0;
        hunt(isotropicArray, inputIso, index);
        if (index == length - 1) { cout << "Value sent to lapse above given range." << endl; return -1; }
        if (index == -1) { cout << "Value sent to lapse: " << inputIso << ", below minimum radius." << endl; return -1; }

        VecDoub truncF(4), truncIso(4);
        int i;

        if (index > 0 && index < length - 2) {
            for (i = -1; i < 3; i++) {
                truncF[i + 1] = array[index + i];
                truncIso[i + 1] = isotropicArray[index + i];
            }
        }

        if (index == 0) {
            for (i = 0; i < 4; i++) {
                truncF[i] = array[index + i];
                truncIso[i] = isotropicArray[index + i];
            }
        }

        if (index == length - 2) {
            for (i = -2; i < 2; i++) {
                truncF[i + 2] = array[index + i];
                truncIso[i + 2] = isotropicArray[index + i];
            }
        }

        Doub err;
        Doub val;
        polint(truncIso, truncF, inputIso, val, err);
        return val;
    }

    Doub searchRatint(Doub inputIso, Doub mass, VecDoub* array_p, int num) {
        // num is for exact values of limits:
        //	1: conformal Factor
        //	2: dConformalFactor
        //	3: ddConformalFactor
        //  0: All Others...

        if (mass != m) {
            m = mass;
            setMass();
        }

        VecDoub& array = *array_p;
        VecDoub& isotropicArray = *isotropicArray_p;

        int index = 0;
        hunt(isotropicArray, inputIso, index);
        if (index == length - 1) { cout << "Value sent to conFactor above given range." << endl; return -1; }
        if (index == -1) { cout << "Value sent to conFactor below minimum radius." << endl; return -1; }
        /*if(inputIso < .001)
          {
          if(num == 1) return sqrt(3.0*m/(2.0*inputIso));
          if(num == 2) return -.5*sqrt(3.0*m/(2.0*inputIso))/inputIso;
          if(num == 3) return .75*sqrt(3.0*m/(2.0*inputIso))/(inputIso*inputIso);
          }*/

        VecDoub truncF(4), truncIso(4);
        int i;

        if (index > 0 && index < length - 2) {
            for (i = -1; i < 3; i++) {
                truncF[i + 1] = array[index + i];
                truncIso[i + 1] = isotropicArray[index + i];
            }
        }

        if (index == 0) {
            for (i = 0; i < 4; i++) {
                //return sqrt(3.0*m/(2.0*inputIso));
                truncF[i] = array[index + i];
                truncIso[i] = isotropicArray[index + i];
            }
        }

        if (index == length - 2) {
            for (i = -2; i < 2; i++) {
                truncF[i + 2] = array[index + i];
                truncIso[i + 2] = isotropicArray[index + i];
            }
        }

        Doub val;
        Doub err;

        ratint(truncIso, truncF, inputIso, val, err);
        return val;
    }



    //===================================================================================================
    //                  Take Derivatives (DOES NOT FILL FIRST INDEX..but that is vere close to 0...)
    //===================================================================================================
    void firstDeriv(VecDoub* f_p, VecDoub* dF_p) {
        VecDoub& f = *f_p, & dF = *dF_p;
        VecDoub& isotropicArray = *isotropicArray_p;
        Doub dr;

        for (int i = 1; i < length - 1; i++) {
            dr = isotropicArray[i + 1] - isotropicArray[i - 1];
            dF[i] = (f[i + 1] - f[i - 1]) / dr;
        }
    }

    void secondDeriv(VecDoub* f_p, VecDoub* ddF_p) {
        VecDoub& f = *f_p, & ddF = *ddF_p;
        VecDoub& isotropicArray = *isotropicArray_p;
        // Doub dr;
        // Doub dr2;

        for (int i = 1; i < length - 1; i++) {
            Doub rp1 = isotropicArray[i + 1];
            Doub rm1 = isotropicArray[i - 1];
            Doub r = isotropicArray[i];
            Doub dfp1o2 = (f[i + 1] - f[i]) / (rp1 - r);
            Doub dfm1o2 = (f[i] - f[i - 1]) / (r - rm1);
            ddF[i] = (dfp1o2 - dfm1o2) / ((rp1 - rm1) / 2.0);
            //		    (f[i+1] - 2.0*f[i] + f[i-1])/dr2;
        }
    }

private:
    VecDoub* arealArray_p, * isotropicArray_p, * conFactorArray_p, * shiftArray_p,
        * lapseArray_p, * dConFactorArray_p, * dLapseArray_p, * ddConFactorArray_p,
        * ddLapseArray_p, * dshiftArray_p, * ddshiftArray_p;
    Doub length, m, arealMax;

    void ratint(VecDoub& xa, VecDoub& ya, const Doub x, Doub& y, Doub& dy) {
        const Doub TINY = 1.0e-25;
        int m, i, ns = 0;
        Doub w, t, hh, h, dd;

        int n = xa.size();
        VecDoub c(n), d(n);
        hh = fabs(x - xa[0]);
        for (i = 0;i < n;i++) {
            h = fabs(x - xa[i]);
            if (h == 0.0) {
                y = ya[i];
                dy = 0.0;
                return;
            } else if (h < hh) {
                ns = i;
                hh = h;
            }
            c[i] = ya[i];
            d[i] = ya[i] + TINY;
        }
        y = ya[ns--];
        for (m = 1;m < n;m++) {
            for (i = 0;i < n - m;i++) {
                w = c[i + 1] - d[i];
                h = xa[i + m] - x;
                t = (xa[i] - x) * d[i] / h;
                dd = t - c[i + 1];
                if (dd == 0.0) toss("Error in routine ratint");
                dd = w / dd;
                d[i] = c[i + 1] * dd;
                c[i] = t * dd;
            }
            y += (dy = (2 * (ns + 1) < (n - m) ? c[ns + 1] : d[ns--]));
        }
    }

    void polint(VecDoub& xa, VecDoub& ya, const Doub x, Doub& y, Doub& dy) {
        int i, m, ns = 0;
        Doub den, dif, dift, ho, hp, w;

        int n = xa.size();
        VecDoub c(n), d(n);
        dif = fabs(x - xa[0]);
        for (i = 0;i < n;i++) {
            if ((dift = fabs(x - xa[i])) < dif) {
                ns = i;
                dif = dift;
            }
            c[i] = ya[i];
            d[i] = ya[i];
        }
        y = ya[ns--];
        for (m = 1;m < n;m++) {
            for (i = 0;i < n - m;i++) {
                ho = xa[i] - x;
                hp = xa[i + m] - x;
                w = c[i + 1] - d[i];
                if ((den = ho - hp) == 0.0) toss("Error in routine polint");
                den = w / den;
                d[i] = hp * den;
                c[i] = ho * den;
            }
            y += (dy = (2 * (ns + 1) < (n - m) ? c[ns + 1] : d[ns--]));
        }
    }

    void hunt(VecDoub& xx, const Doub x, int& jlo) {
        int jm, jhi, inc;
        bool ascnd;

        int n = xx.size();
        ascnd = (xx[n - 1] >= xx[0]);
        if (jlo < 0 || jlo > n - 1) {
            jlo = -1;
            jhi = n;
        } else {
            inc = 1;
            if (x >= xx[jlo] == ascnd) {
                if (jlo == n - 1) return;
                jhi = jlo + 1;
                while (x >= xx[jhi] == ascnd) {
                    jlo = jhi;
                    inc += inc;
                    jhi = jlo + inc;
                    if (jhi > n - 1) {
                        jhi = n;
                        break;
                    }
                }
            } else {
                if (jlo == 0) {
                    jlo = -1;
                    return;
                }
                jhi = jlo--;
                while (x < xx[jlo] == ascnd) {
                    jhi = jlo;
                    inc <<= 1;
                    if (inc >= jhi) {
                        jlo = -1;
                        break;
                    } else jlo = jhi - inc;
                }
            }
        }
        while (jhi - jlo != 1) {
            jm = (jhi + jlo) >> 1;
            if (x >= xx[jm] == ascnd)
                jlo = jm;
            else
                jhi = jm;
        }
        if (x == xx[n - 1]) jlo = n - 2;
        if (x == xx[0]) jlo = 0;
    }
};

#endif
