#ifndef FREEZEIN_H
#define FREEZEIN_H

/********************************/
/* Standard and Boost Libraries */
/********************************/

//Standard libraries
#include <iostream>/*standard i/o library*/
#include <fstream>/*to read/write file*/
#include <cstdio>/*provides printf*/
#include <string>/*to use string data type*/
#include <cmath>/*provides pow, sqrt, ...*/
#include <vector>/*provides std::vector*/
#include <functional> /*std::function*/
#include <limits> /* std::numeric_limits */

//Boost C++ library
#include <boost/math/special_functions/bessel.hpp>/*provides bessel-K function*/
#include <boost/math/quadrature/gauss.hpp>/*provides Gauss-Legendre quadrature*/
#include <boost/math/quadrature/exp_sinh.hpp>/*provides exp_sinh quadrature*/
#include <boost/math/tools/roots.hpp> /*roots of a function*/

//Namespaces
using namespace std;
using namespace boost::math;
using namespace boost::math::quadrature;
using namespace boost::math::tools;

/*********************************************************/
/* Masses, widths, coupling constants, and mixing angles */
/*********************************************************/

//Masses in GeV
#define Me 0.51099895e-3L
#define Mmu 105.6583755e-3L
#define Mta 1.77686L
#define Mu 2.16e-3L
#define Mc 1.27L
#define Mt 172.69L
#define Md 4.67e-3L
#define Ms 93.4e-3L
#define Mb 4.18L
#define Mpip 139.57039e-3L
#define MKp 493.677e-3L
#define MZ 91.1876L
#define MW 80.379L
#define MPl 2.435323077e18L

//Widths in GeV
#define WZ 2.4952L

//Alpha strong at mZ
#define alphaS 0.1179L

//Fine structure constant
#define alphaEM 7.2973525664e-3L
const long double eEM = sqrt(4.0L*M_PI*alphaEM); /*EM coupling*/

//Weinberg angle
#define sW2 0.23121L /*Sin-squared Weinberg angle*/
const long double sW = sqrt(sW2); /*Sin(ThetaW)*/
const long double cW = sqrt(1.0L - sW2); /*Cos(ThetaW)*/
const long double tW = sqrt(sW2/(1.0L - sW2)); /*Tan(ThetaW)*/
const long double s2W = 2.0L*sqrt(sW2 - sW2*sW2); /*Sin(2 ThetaW)*/

//Conversions
#define GeVinvtocm 1.97326937e-14 /*Inverse GeV in cm*/

/************************************************************************/
/* Useful fuctions: linear interpolation, finite-difference derivatives */
/************************************************************************/

//Linear interpolation function
long double interp(long double x, const vector<long double> &xData,
              const vector<long double> &yData, bool extrapolate) {

    //Check if x is increasing/decreasing
    bool increasing = xData[1] > xData[0];

    //Perform binary search to find the interval for interpolation
    int mid, low = 0, high = xData.size() - 1;
    while (high - low > 1) {
        
        mid = (low + high)/2;

        if (
                (increasing && x >= xData[mid]) ||
                (!increasing && x <= xData[mid])
           ) {
            low = mid;
        }
        else {
            high = mid;
        }
    }

    //Points on either side (unless beyond ends)
    long double xL = xData[low], yL = yData[low];
    long double xR = xData[low+1], yR = yData[low+1];

    if ( !extrapolate ) {//if beyond ends of array and not extrapolating
       if ( (increasing && x < xL) || (!increasing && x > xL) ) yR = yL;
       if ( (increasing && x > xR) || (!increasing && x < xR) ) yL = yR;
    }

    long double dydx = (yR - yL)/(xR - xL);//gradient
    
    return yL + dydx * ( x - xL );//linear interpolation
}

//Array of slopes using Finite-difference formulas
vector<long double> slopearray(const vector<long double>& xData,
                          const vector<long double>& yData) {

    long double x, y, x1, y1, x2, y2;
    vector<long double> dydxData;

    for (size_t i = 0; i < xData.size(); i++) {

        x = xData[i]; y = yData[i];
        if ( i == 0 ) {
            x1 = xData[i+1]; y1 = yData[i+1]; x2 = xData[i+2]; y2 = yData[i+2];
        }
        else if ( i == xData.size() - 1 ) {
            x1 = xData[i-1]; y1 = yData[i-1]; x2 = xData[i-2]; y2 = yData[i-2];
        }
        else {
            x1 = xData[i-1]; y1 = yData[i-1]; x2 = xData[i+1]; y2 = yData[i+1];
        }

        dydxData.push_back(y*(2.0L*x - (x1 + x2))/(x - x1)/(x - x2) +
                           y1*(x - x2)/(x1 - x)/(x1 - x2) +
                           y2*(x - x1)/(x2 - x)/(x2 - x1));
    }

    return dydxData;
}

//Root-solver for minimizing an increasing/decreasing function f(x)
long double MinimizeCpp(std::function<long double(long double)> f,
                        long double xguess, bool increasing,
                        long double factor, int niter) {

    /* Max no. of iterations while finding root */
    const boost::uintmax_t maxit = niter;

    /* Actual no. of iterations, updated during root finding */
    boost::uintmax_t it = maxit;

    bool is_rising = increasing;

    /* Max. possible binary digits accuracy for long double type */
    int digits = numeric_limits<long double>::digits;

    /* Solve for root that is accurate to chosen number of digits */
    int get_digits = (digits * 1) /2;

    eps_tolerance<long double> tol(get_digits); /* Tolerance */

    pair<long double, long double> r = bracket_and_solve_root(f, xguess,
                                                              factor, is_rising,
                                                              tol, it);

    /* If root solver is unable to find root with niter iterations */
    if(it >= maxit){
        if ( 2.0L*abs(r.second - r.first)/(r.second + r.first) >= 0.1 ) {
            cout << "guess: x = " << xguess << endl;
            cout << "Unable to locate solution with " << maxit <<
                    " iterations, and current best guess is between " <<
                    r.first << " and " << r.second << endl;
        }
    }

    /* Root is midway between brackets */
    long double result = r.first + (r.second - r.first)/2.0L;

    return result; 
}

/***********************************************************/
/* g*(S): Effective number of degrees of freedom in the SM */
/***********************************************************/

//Define arrays for Temperature T and gstar(S) in the SM
vector<long double> Tvec;

vector<long double> gstarvec;
vector<long double> dlngstardlnTvec;

vector<long double> gstarSvec;
vector<long double> dlngstarSdlnTvec;

//Read gstar(S) data from various .tab files in the gstar folder, and compute
//their finite-difference derivatives as a function of T
void Read_gstar(const string& choice, const string& gstarpath) {
    
    long double r1, r2, r3;
    vector <long double> tempvec;

    //Clear global vectors before filling them
    Tvec.clear(), Tvec.shrink_to_fit();
    gstarvec.clear(), gstarvec.shrink_to_fit();
    gstarSvec.clear(), gstarSvec.shrink_to_fit();
    dlngstardlnTvec.clear(), dlngstardlnTvec.shrink_to_fit();
    dlngstarSdlnTvec.clear(), dlngstarSdlnTvec.shrink_to_fit();

    //Open a file corresponding to the input choice
    string filename;
    if (choice == "standard") { filename = "gstar/std.tab"; }
    else if (choice == "HP_A") { filename = "gstar/HP_A.tab"; }
    else if (choice == "HP_B") { filename = "gstar/HP_B.tab"; }
    else if (choice == "HP_B2") { filename = "gstar/HP_B2.tab"; }
    else if (choice == "HP_B3") { filename = "gstar/HP_B3.tab"; }
    else if (choice == "HP_C") { filename = "gstar/HP_C.tab"; }

    //prepend the path to the gstar folder to the filename
    filename = gstarpath + "/" + filename;

    ifstream file(filename);
    //If a file is open, read line-by-line to extract three values (r1, r2, r3)
    //from each line and store them in the corresponding global vectors
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;//Skip lines with '#'
            istringstream iss(line);
            iss >> r1 >> r2 >> r3;
            Tvec.push_back(r1);
            gstarvec.push_back(r3);
            gstarSvec.push_back(r2);
        }
        file.close();
    }
    else {
        cout << "Unable to open the file:" << filename << endl;
    }

    tempvec = slopearray(Tvec, gstarvec);
    for (size_t i = 0; i < Tvec.size(); i++) {
        dlngstardlnTvec.push_back((Tvec[i]/gstarvec[i])*tempvec[i]);
    }
    tempvec = slopearray(Tvec, gstarSvec);
    for (size_t i = 0; i < Tvec.size(); i++) {
        dlngstarSdlnTvec.push_back((Tvec[i]/gstarSvec[i])*tempvec[i]);
    }

    return;
}

//g*
long double gstar(long double T) {
    return interp(T, Tvec, gstarvec, false);
}

//g*S
long double gstarS(long double T) {
    return interp(T, Tvec, gstarSvec, false);
}

//dlng*S/dlnT
long double dlngstarSdlnT(long double T) {
    return interp(T, Tvec, dlngstarSdlnTvec, false);
}

//dlng*/dlnT
long double dlngstardlnT(long double T) {
    return interp(T, Tvec, dlngstardlnTvec, false);
}

/***************************************************************************/
/* Energy density, comoving entropy, and Hubble rate in the Visible sector */
/***************************************************************************/

//Rho Visible
long double RhoVisible(long double T) {
    return (M_PI*M_PI/30.0L)*gstar(T)*pow(T, 4.0L);
}

//Comoving entropy of the Visible sector
long double EntropyVisible(long double T) {
    return (2.0L*M_PI*M_PI/45.0L)*gstarS(T)*T*T*T;
}

//Hubble rate
long double Hubble(long double T) {
    return sqrt(M_PI*M_PI*gstar(T)/90.0L)*T*T/MPl;
}

//(H / Hbar) to account for varying gstarS only in the Visible sector
long double HoverHbarVisible(long double T) {
    return (1.0L + (1.0L/3.0L)*dlngstarSdlnT(T));
}

/**************************************/
/* Equilibrium number density for Chi */
/**************************************/

long double NumEq(long double T, long double m, int dof) {

    return (dof/(2.0L*M_PI*M_PI))*T*m*m*boost::math::cyl_bessel_k(2, m/T);

}

/**********************/
/* Plasmon properties */
/**********************/

//Define arrays for temperature, Omega1, Omegap, and vstar
vector<long double> Tvector;

vector<long double> Omega1vec;
vector<long double> Omegapvec;
vector<long double> vstarvec;

//Read Omega1, Omegap, and vstar data as a function of T
//from OmegaInterpolation.tab file in the gstar folder
void Read_OmegasFile(const string& gstarpath) {
    
    long double r1, r2, r3, r4;

    //Clear global vectors before filling them
    Tvector.clear(), Tvector.shrink_to_fit();
    Omega1vec.clear(), Omega1vec.shrink_to_fit();
    Omegapvec.clear(), Omegapvec.shrink_to_fit();
    vstarvec.clear(), vstarvec.shrink_to_fit();

    //Open the file 
    string filename = "gstar/OmegaInterpolation.tab"; //Computed in Jupyterlab
    
    //prepend the path to the gstar folder to the filename
    filename = gstarpath + "/" + filename;
    
    ifstream file(filename);
    //If a file is open, read line-by-line to extract 4 values
    //(r1, r2, r3, r4) from each line and store them in
    //the corresponding global vectors
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            if (line.empty() || line[0] == '#') continue; //Skip lines with '#'
            istringstream iss(line);
            iss >> r1 >> r2 >> r3 >> r4;
            Tvector.push_back(r1);
            Omega1vec.push_back(r2);
            Omegapvec.push_back(r3);
            vstarvec.push_back(r4);
        }
        file.close();
    }
    else {
        cout << "Unable to open the file:" << filename << endl;
    }
    
    return;
}

//Plasma frequency (Eq B2 in DLS 2019)
long double Omegap(long double T, bool UseInterpolation) {

    if ( UseInterpolation ) { return interp(T, Tvector, Omegapvec, false); }

    if ( T >= 50.0L*Me ) {//high T
        return eEM*T/3.0L;
    }
    else if ( T <= Me/500.0L ) {//low T
        return eEM*sqrt( (NumEq(T, Me, 2)/Me)*(1.0L - 2.5L*(T/Me)) );
    }
    else {
        auto integrand_p = [=] (long double p) {

            long double E = sqrt(p*p + Me*Me);
            long double v = p/E;

            //Using Fermi-Dirac statistics
            return (p*p/E)*(1.0L - v*v/3.0L)*(1.0L/(1.0L + exp(E/T)));
        };
        
        return sqrt(2.0L * (4.0L*alphaEM/M_PI) *
                exp_sinh<long double>().integrate(integrand_p, 0.0L, INFINITY));
    }
}

//First-mode frequency (Eq B1 in DLS 2019)
long double Omega1(long double T,  bool UseInterpolation) {

    if ( UseInterpolation ) { return interp(T, Tvector, Omega1vec, false); }

    if ( T >= 50.0L*Me ) {//high T
        long double pF = pow(3.0L*M_PI*M_PI*NumEq(T, Me, 2), 1.0L/3.0L);

        return (eEM*T/3.0L)*(pF/sqrt(pF*pF + Me*Me));
    }
    else if ( T <= Me/500.0L ) {//low T
        return (eEM/Me)*sqrt( 5.0L*T*NumEq(T, Me, 2)*(1.0L - 2.5L*(T/Me)) );
    }
    else {
        auto integrand_p = [=] (long double p) {

            long double E = sqrt(p*p + Me*Me);
            long double v = p/E;

            //Using Fermi-Dirac statistics
            return (p*p/E)*((5.0L/3.0L)*v*v - v*v*v*v) *
                   (1.0L/(1.0L + exp(E/T)));
        };
        
        return sqrt(2.0L * (4.0L*alphaEM/M_PI) * 
                exp_sinh<long double>().integrate(integrand_p, 0.0L, INFINITY));
    }
}

//Typical electron velocity (Above Eq B1 in DLS 2019)
long double vstar(long double T, bool UseInterpolation) {

    if ( UseInterpolation ) { return interp(T, Tvector, vstarvec, false); }

    if ( T >= 50.0L*Me ) {//high T
        long double pF = pow(3.0L*M_PI*M_PI*NumEq(T, Me, 2), 1.0L/3.0L);

        return pF/sqrt(pF*pF + Me*Me);
    }
    else if ( T <= Me/500.0L ) {//low T
        return sqrt(5.0L*T/Me);
    }
    else {
        return Omega1(T, false)/Omegap(T, false);
    }
}

/**********************/
/* Transverse plasmon */
/**********************/

//Transverse plasmon frequency (Eq B8 in DLS 2019)
long double Omegat(long double k, long double T,
                   long double factor, int niter,
                   bool UseInterpolation) {

    long double omp = Omegap(T, UseInterpolation);

    if ( T <= Me/500.0L ) {
        if ( k >= 1e3L*omp ) {//low T, high k
            long double mt = omp*sqrt(1.0L + T/Me);
            
            return sqrt(mt*mt + k*k);
        }
        else {//low T, low/mid k
            long double temp = k*k + omp*omp;

            return sqrt(0.5L*temp*(1.0L +
                               sqrt(1.0L + 4.0L*T*k*k*omp*omp/(Me*temp*temp))));
        }
    }
    else if ( (T >= 50.0L*Me) && (k >= 1e3L*omp) ) {//high T, high k
        return sqrt(1.5L*omp*omp + k*k);
    }
    else {
        long double om1 = Omega1(T, UseInterpolation);
        long double v = om1/omp;
        long double vk = v*k;

        if ( k <= 1e-3L*omp ) {//mid/high T, low k
            return omp;
        }
        else if ( k >= 1e3L*omp ) {//mid T, high k
            long double mt = omp*sqrt((3.0L/(2.0L*v*v)) *
                             (1.0L - ((1.0L - v*v)/(2.0L*v)) *
                                     log((1.0L + v)/(1.0L - v)) ));

            return sqrt(mt*mt + k*k);
        }
        else {//mid/high T, mid k
            
            //Transverse plasmon polarization function (Eq B5 in DLS 2019)
            auto solve_omega = [=] (long double omega) {
                long double x, PlasmonPIt;

                if ( omega < vk ) {
                    PlasmonPIt = std::numeric_limits<long double>::max();
                }
                else if ( vk <= 1e-9L*omega ) {
                    PlasmonPIt = omp*omp;
                }
                else if ( abs(omega - vk) <= 1e-9L*vk ) {
                    PlasmonPIt = (3.0L/2.0L)*omp*omp;
                }
                else {
                    x = omega/vk;
                    PlasmonPIt = (3.0L/2.0L)*omp*omp*(x*x -
                                0.5L*x*(x*x - 1.0L)*log((x + 1.0L)/(x - 1.0L)));
                }

                return PlasmonPIt + k*k - omega*omega;
            };

            long double omegaguess = sqrt(omp*omp + k*k);
            bool increasing = false;
            
            return MinimizeCpp(solve_omega, omegaguess, increasing,
                               factor, niter);  
        }
    }
}

//Transverse plasmon mass
long double PlasmonMasst(long double k, long double T,
                         long double factor, int niter,
                         bool UseInterpolation) {

    long double omt = Omegat(k, T, factor, niter, UseInterpolation);

    return ( omt >= k ) ? sqrt(omt*omt - k*k) : 0.0L;

}

//Transverse renormalization factor (Eq B11 in DLS 2019)
long double PlasmonZt(long double k, long double T,
                      long double factor, int niter,
                      bool UseInterpolation) {

    long double omp = Omegap(T, UseInterpolation);
    long double v = vstar(T, UseInterpolation);
    long double mt = PlasmonMasst(k, T, factor, niter, UseInterpolation);

    long double result = 0.0L;

    if ( mt >= 0.0L ) {
        if ( k >= 1e3L*omp ) {
            result = 1.0L;
        }
        else {
            long double temp = mt*mt + k*k*(1.0L - v*v);
            long double omt = sqrt(mt*mt + k*k);

            result = 2.0L*omt*omt*temp/( 3.0L*omp*omp*omt*omt +
                                         temp*(omt*omt + k*k) -
                                         2.0L*omt*omt*(omt*omt - k*k) );
        }
    }

    return result;
}

/************************/
/* Longitudinal plasmon */
/************************/

//kmax
long double kMax(long double T, bool UseInterpolation) {

    long double omp = Omegap(T, UseInterpolation);
    long double om1 = Omega1(T, UseInterpolation);
    long double v = om1/omp;

    return omp*sqrt((3.0L/(v*v)) * 
           ((1.0L/(2.0L*v))*log((1.0L + v)/(1.0L - v)) - 1.0L));
} 

//Longitudinal plasmon frequency (Eq B8 in DLS 2019)
long double Omegal(long double k, long double T,
                   long double factor, int niter,
                   bool UseInterpolation) {

    long double omp = Omegap(T, UseInterpolation);

    if ( T <= Me/500.0L ) {
        long double kmax = omp*sqrt(1 + 3.0L*T/Me);
        
        if ( k >= kmax) {//low T, high k
            return k;
        }
        else {//low T, low/mid k
            return omp *
                      sqrt(0.5L*(1.0L + sqrt(1.0L + 12.0L*T*k*k/(Me*omp*omp))));
        }
    }
    else {
        long double om1 = Omega1(T, UseInterpolation);
        long double v = om1/omp;
        long double vk = v*k;

        long double kmax = omp*sqrt((3.0L/(v*v)) * 
                           ((1.0L/(2.0L*v))*log((1.0L + v)/(1.0L - v)) - 1.0L));

        if ( vk <= 1e-3L*omp ) {//mid/high T, low k
            return omp;
        }
        else if ( k >= kmax ) {//mid/high T, high k 
            return k;
        }
        else {//mid/high T, mid k

            //Longitudinal plasmon polarization function (Eq B4 in DLS 2019)
            auto solve_omega = [=] (long double omega) {

                long double x, PlasmonPIl;

                if ( omega < vk ) {
                    PlasmonPIl = std::numeric_limits<long double>::max();
                }
                else if ( vk <= 1e-9L*omega ) {
                    PlasmonPIl = omp*omp*k*k/(omega*omega);
                }
                else if ( abs(omega - vk) <= 1e-9L*vk ) {
                    PlasmonPIl = std::numeric_limits<long double>::max();
                }
                else {
                    x = omega/vk;
                    PlasmonPIl = 3.0L*(omp*omp/(v*v)) * 
                                 ( (x/2.0L)*log((x + 1.0L)/(x - 1.0L)) - 1.0L );
                }

                return (omega*omega/(k*k))*PlasmonPIl - omega*omega;
            };
            
            long double mlguess = omp*pow(1.0L - k/kmax, 6.0L);
            long double omegaguess = sqrt(mlguess*mlguess + k*k);
            bool increasing = false;
            
            return MinimizeCpp(solve_omega, omegaguess, increasing,
                               factor, niter);  
        }
    }
}

//Longitudinal plasmon mass
long double PlasmonMassl(long double k, long double T,
                         long double factor, int niter,
                         bool UseInterpolation) {

    long double oml = Omegal(k, T, factor, niter, UseInterpolation);

    return ( oml >= k ) ? sqrt(oml*oml - k*k) : 0.0L;
}

//Longitudinal renormalization factor (Eq B10 in DLS 2019)
long double PlasmonZl(long double k, long double T,
                      long double factor, int niter,
                      bool UseInterpolation) {

    long double omp = Omegap(T, UseInterpolation);
    long double v = vstar(T, UseInterpolation);
    long double ml = PlasmonMassl(k, T, factor, niter, UseInterpolation);

    long double result = 0.0L;

    if ( ml > 0.0L ) {
        long double temp = ml*ml + k*k*(1.0L - v*v);

        result = 2.0L*temp/(3.0L*omp*omp - temp);
    }

    return result;
}

/******************************************/
/* Fully averaged squared Matrix elements */
/******************************************/

//Fully averaged matrix element squared for f f -> Aprime/Z -> chi chi
long double M2_ffchichi(long double s, long double mchi, long double mf,
                        long double kappa, long double Nf, long double Qf,
                        long double T3f) {

    //Vector (Vf) and Axial (Af) pieces of Z-f-f couplings in the SM
    long double Vf = (T3f - 2.0L*sW2*Qf)/s2W;
    long double Af = T3f/s2W;

    return (32.0L/3.0L)*Nf*M_PI*M_PI*alphaEM*alphaEM*kappa*kappa*(
            //Aprime
            Qf*Qf*(1.0L + 2.0L*mf*mf/s)*(1.0L + 2.0L*mchi*mchi/s)
            +
            //Aprime/Z interference
            tW*tW*(Vf*Vf*(s + 2.0L*mf*mf) + Af*Af*(s - 4.0L*mf*mf))*(
                (s + 2.0L*mchi*mchi)/(pow(s - MZ*MZ, 2.0L) + MZ*MZ*WZ*WZ)
                )
            -
            //Z
            2.0L*Qf*Vf*tW*(s + 2.0L*mf*mf)*(s + 2.0L*mchi*mchi)*(
                (1.0L - MZ*MZ/s)/(pow(s - MZ*MZ, 2.0L) + MZ*MZ*WZ*WZ)
                )
            );
}

//Fully averaged matrix element squared for phi+ phi- -> Aprime -> chi chi
long double M2_sschichi(long double s, long double mchi, long double mcs,
                        long double kappa) {

    return (32.0L/3.0L)*M_PI*M_PI*alphaEM*alphaEM*kappa*kappa*(
            (1.0L - 4.0L*mcs*mcs/s)*(1.0L + 2.0L*mchi*mchi/s)
            );
}

//Fully averaged matrix element squared for W+ W- -> Aprime/Z -> chi chi
long double M2_WWchichi(long double s, long double mchi, long double kappa) {

    return (8.0L/27.0L)*M_PI*M_PI*alphaEM*alphaEM*kappa*kappa*(
            (MZ*MZ*MZ*MZ)/(MW*MW*MW*MW)
            )*(1.0L + 2.0L*mchi*mchi/s)*(1.0L - 4.0L*MW*MW/s)*(
                (s*s + 20.0L*MW*MW*s + 12.0L*MW*MW*MW*MW)/(
                    pow(s - MZ*MZ, 2.0L) + MZ*MZ*WZ*WZ
                    )
            );
}

/****************************************/
/* Collision terms for number densities */
/****************************************/

//Number-density collision term for f f -> Aprime/Z -> Chi Chi
long double CollisionNum_ffchichi(long double T, long double mchi,
                                  long double mf, long double kappa,
                                  long double Nf, long double Qf,
                                  long double T3f, long double LambdaQCD) {

    if ( ( Nf == 1.0L ) || ( (Nf == 3.0L) && (T > LambdaQCD) ) ) {

        auto integrand_s = [=] (long double s) {
            return M2_ffchichi(s, mchi, mf, kappa, Nf, Qf, T3f) *
                   sqrt(1.0L - 4.0L*mchi*mchi/s) *
                   sqrt(1.0L - 4.0L*mf*mf/s) *
                   sqrt(s) *
                   boost::math::cyl_bessel_k(1, sqrt(s)/T);
        };
        
        return (16.0L*T/pow(4.0L*M_PI, 5.0L)) *
               exp_sinh<long double>().integrate(integrand_s,
                                            max(4.0L*mf*mf, 4.0L*mchi*mchi),
                                            INFINITY);
    }
    else { return 0.0L; }
}

//Number-density collision term for phi+ phi- -> Aprime -> Chi Chi
long double CollisionNum_sschichi(long double T, long double mchi,
                                  long double mcs, long double kappa,
                                  long double LambdaQCD) {

    if ( T <= LambdaQCD ) {

        auto integrand_s = [=] (long double s) {
            return M2_sschichi(s, mchi, mcs, kappa) *
                   sqrt(1.0L - 4.0L*mchi*mchi/s) *
                   sqrt(1.0L - 4.0L*mcs*mcs/s) *
                   sqrt(s) *
                   boost::math::cyl_bessel_k(1, sqrt(s)/T);
        };
        
        return (4.0L*T/pow(4.0L*M_PI, 5.0L)) *
               exp_sinh<long double>().integrate(integrand_s,
                                            max(4.0L*mcs*mcs, 4.0L*mchi*mchi),
                                            INFINITY);
    }
    else { return 0.0L; }
}

//Number-density collision term for W+ W- -> Aprime -> Chi Chi
long double CollisionNum_WWchichi(long double T, long double mchi,
                                  long double kappa) {

    auto integrand_s = [=] (long double s) {
        return M2_WWchichi(s, mchi, kappa) *
               sqrt(1.0L - 4.0L*mchi*mchi/s) *
               sqrt(1.0L - 4.0L*MW*MW/s) *
               sqrt(s) *
               boost::math::cyl_bessel_k(1, sqrt(s)/T);
    };
    
    return (36.0L*T/pow(4.0L*M_PI, 5.0L)) *
           exp_sinh<long double>().integrate(integrand_s,
                                        max(4.0L*MW*MW, 4.0L*mchi*mchi),
                                        INFINITY);
}

//Number-density collision term for plasmon_t -> chi chi (Eq 22 in DLS 2019)
long double CollisionNum_Plasmontchichi(long double T, long double mchi,
                                        long double kappa,
                                        bool UseInterpolation) {
    
    long double omp = Omegap(T, UseInterpolation);
    long double v = vstar(T, UseInterpolation);

    auto integrand_kt = [=] (long double k) {

        long double mt = PlasmonMasst(k, T, 2.0L, 25, UseInterpolation);
        long double result = 0.0L;

        if ( mt >= 2.0L*mchi ) {
            long double omt = sqrt(mt*mt + k*k);
            long double temp = mt*mt + k*k*(1.0L - v*v);
            long double Zt;

            if ( k >= 1e3L*omp ) {
                Zt = 1.0L;
            }
            else {
                Zt = 2.0L*omt*omt*temp/( 3.0L*omp*omp*omt*omt +          
                                         temp*(omt*omt + k*k) -
                                         2.0L*omt*omt*(omt*omt - k*k) );
            }

            result = (k*k/3.0L)*Zt*(1.0L/(omt*(exp(omt/T) - 1.0L))) *
                     (mt*mt + 2.0L*mchi*mchi) *
                     sqrt(1.0L - 4.0L*mchi*mchi/(mt*mt));

        }

        return result;
        
    };

    return ((alphaEM*kappa*kappa)/(M_PI*M_PI)) *
           gauss<long double, 1501>().integrate(integrand_kt, 0.0L, INFINITY);
}

//Number-density collision term for plasmon_l -> chi chi (Eq 21 in DLS 2019)
long double CollisionNum_Plasmonlchichi(long double T, long double mchi,
                                        long double kappa,
                                        bool UseInterpolation) {

    long double omp = Omegap(T, UseInterpolation);
    long double v = vstar(T, UseInterpolation);

    auto integrand_kl = [=] (long double k) {

        long double ml = PlasmonMassl(k, T, 2.0L, 25, UseInterpolation);
        long double result = 0.0L;

        if ( ml >= 2.0L*mchi ) {
            long double temp = ml*ml + k*k*(1.0L - v*v);
            long double Zl = 2.0L*temp/(3.0L*omp*omp - temp);
            long double oml = sqrt(ml*ml + k*k);

            result = (k*k/3.0L)*Zl*(oml/(exp(oml/T) - 1.0L)) *
                     (1.0L + 2.0L*mchi*mchi/(ml*ml)) *
                     sqrt(1.0L - 4.0L*mchi*mchi/(ml*ml));
        }

        return result;
    };

    return ((alphaEM*kappa*kappa)/(2.0L*M_PI*M_PI)) *
           exp_sinh<long double>().integrate(integrand_kl, 0.0L, INFINITY);
}

//Sum of all number-density collision terms for portal freeze-in
long double CollisionNum_chi(long double T, long double mchi,
                             long double kappa, long double LambdaQCD,
                             bool IncludePlasmons, bool UseInterpolation) {

    long double result = 0.0L;

    //SM SM -> Chi Chi
    result += CollisionNum_ffchichi(T, mchi, 0.0L, kappa,
                                    1.0L, 0.0L, 0.5L,
                                    LambdaQCD)*3.0L + /*neutrinos*/
              CollisionNum_ffchichi(T, mchi, Me, kappa,
                                    1.0L, -1.0L, -0.5L,
                                    LambdaQCD) + /*e*/
              CollisionNum_ffchichi(T, mchi, Mmu, kappa,
                                    1.0L, -1.0L, -0.5L,
                                    LambdaQCD) + /*mu*/
              CollisionNum_ffchichi(T, mchi, Mta, kappa,
                                    1.0L, -1.0L, -0.5L,
                                    LambdaQCD) + /*ta*/
              CollisionNum_ffchichi(T, mchi, Mu, kappa,
                                    3.0L, 2.0L/3.0L, 0.5L,
                                    LambdaQCD) + /*u*/
              CollisionNum_ffchichi(T, mchi, Mc, kappa,
                                    3.0L, 2.0L/3.0L, 0.5L,
                                    LambdaQCD) + /*c*/
              CollisionNum_ffchichi(T, mchi, Mt, kappa,
                                    3.0L, 2.0L/3.0L, 0.5L,
                                    LambdaQCD) + /*t*/
              CollisionNum_ffchichi(T, mchi, Md, kappa,
                                    3.0L, -1.0L/3.0L, -0.5L,
                                    LambdaQCD) + /*d*/
              CollisionNum_ffchichi(T, mchi, Ms, kappa,
                                    3.0L, -1.0L/3.0L, -0.5L,
                                    LambdaQCD) + /*s*/
              CollisionNum_ffchichi(T, mchi, Mb, kappa,
                                    3.0L, -1.0L/3.0L, -0.5L,
                                    LambdaQCD) + /*b*/
              CollisionNum_sschichi(T, mchi, Mpip, kappa,
                                    LambdaQCD) + /*pi+*/
              CollisionNum_sschichi(T, mchi, MKp, kappa,
                                    LambdaQCD) + /*K+*/
              CollisionNum_WWchichi(T, mchi, kappa); /*W*/

    //Plasmon -> Chi Chi
    if (
            ( IncludePlasmons ) &&
            (mchi >= 9.9e-10L) && (mchi <= 1.1L) &&
            (T >= Me/20.0L)  && (T <= 1.0L)
       ) {

        result += CollisionNum_Plasmontchichi(T, mchi, kappa,
                                              UseInterpolation); /*gammastar_t*/
        result += CollisionNum_Plasmonlchichi(T, mchi, kappa,
                                              UseInterpolation); /*gammastar_l*/

    }

    return result;
}

/*********************************************************/
/* Thermally-averaged cross-section for portal freeze-in */
/*********************************************************/

//Thermally-averaged cross section for SM SM -> chi chi
long double SigmaV_chi(long double T, long double mchi, long double kappa,
                       long double LambdaQCD, bool IncludePlasmons,
                       bool UseInterpolation) {
    
    return CollisionNum_chi(T, mchi, kappa, LambdaQCD,
                            IncludePlasmons, UseInterpolation) /
           pow(NumEq(T, mchi, 2), 2.0L);
}

/*****************************/
/* Freeze-in portal coupling */
/*****************************/

//Portal Yield for chi
long double Yield_FreezeIn(long double mchi, long double kappa,
                           long double LambdaQCD, bool IncludePlasmons,
                           bool UseInterpolation) {

    long double result;
    long double CollisionTerm;

    auto integrand_T = [&] (long double T) {

        CollisionTerm = CollisionNum_chi(T, mchi, kappa, LambdaQCD,
                                         IncludePlasmons, UseInterpolation);

        return HoverHbarVisible(T) *
               CollisionTerm /
               (gstarS(T)*sqrt(gstar(T))*pow(T, 6.0L));
    };

    result = gauss<long double, 3001>().integrate(integrand_T, 0.0L, INFINITY);

    return (135.0L*sqrt(10.0L)*MPl/(2.0L*pow(M_PI, 3.0L)))*result;
}

//Portal coupling, kappa, for freezing-in the required relic abundance
long double kappa_FreezeIn(long double mchi, long double LambdaQCD,
                           bool IncludePlasmons, bool UseInterpolation) {
    return sqrt(
                4.37e-10L /
                (2.0L * mchi * Yield_FreezeIn(mchi, 1.0L, LambdaQCD,
                                              IncludePlasmons,
                                              UseInterpolation))
               );
}

/**********************************/
/* Direct detection cross section */
/**********************************/

//Reduced mass of Chi and electron
long double Muchie(long double mchi) {
    return mchi*Me/(mchi + Me);
}

//Direct detection cross section in squared-centimeter: \overline{\sigma}_e
long double SigmaDDe(long double mchi, long double kappa) {
    return 16.0L*M_PI*alphaEM*alphaEM*kappa*kappa *
           pow(Muchie(mchi), 2.0L) * pow(GeVinvtocm, 2.0L) /
           pow((alphaEM * Me), 4.0L);
}

#endif
