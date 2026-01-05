#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

int main() {
    // ---- Parameters (from your Fortran) ----
    const double ETA  = 0.4;     // not used in this snippet, but kept
    const int    NTOT = 64;

    const double AL   = 6.2;     // box side length in original code
    const double ALH  = 0.5 * AL;

    const int ND = static_cast<int>(AL) + 1; // Fortran: ND = INT(AL) + 1
    const double DR = 1.0;

    std::cout << std::fixed << std::setprecision(16);
    std::cout << "ETA, AL, DR = " << ETA << " " << AL << " " << DR << "\n";

    // ---- Hard-sphere diameter implied by overlap test R2 < 1.0 ----
    const double d = 1.0;
    const double a = 0.5 * d;

    // ---- Spherical container: radius R = ALH, but keep centers inside R-a ----
    const double Redge  = ALH - a;
    const double Redge2 = Redge * Redge;

    std::vector<double> X, Y, Z;
    X.reserve(NTOT);
    Y.reserve(NTOT);
    Z.reserve(NTOT);

    int IC = 0;

    for (int I = 1; I <= ND - 1; ++I) {
        for (int J = 1; J <= ND - 1; ++J) {
            for (int K = 1; K <= ND - 1; ++K) {

                const double Xtry = -ALH + (I - 0.5) * DR;
                const double Ytry = -ALH + (J - 0.5) * DR;
                const double Ztry = -ALH + (K - 0.5) * DR;

                // ---- Spherical inside check ----
               // const double r2 = Xtry*Xtry + Ytry*Ytry + Ztry*Ztry;
                //if (r2 > Redge2) continue;

                // (A) SPHERE inside check  (COMMENTED OUT)
                // ============================================================
                // const double r2 = Xtry*Xtry + Ytry*Ytry + Ztry*Ztry;
                // if (r2 > Redge2) continue;

                // ============================================================
                // (B) CUBE inside check  (ACTIVE)
                // Keep center inside cube walls by radius a:
                // |x| <= AL/2 - a, |y| <= AL/2 - a, |z| <= AL/2 - a
                // ============================================================
                  if (std::abs(Xtry) > ALH - a) continue;
                  if (std::abs(Ytry) > ALH - a) continue;
                  if (std::abs(Ztry) > ALH - a) continue;






                // ---- Overlap check with previously accepted particles ----
                bool overlap = false;
                for (int m = 0; m < IC; ++m) {
                    const double tx = X[m] - Xtry;
                    const double ty = Y[m] - Ytry;
                    const double tz = Z[m] - Ztry;
                    const double rij2 = tx*tx + ty*ty + tz*tz;

                    // Fortran: IF(R2+1e-6 .LT. 1.0) reject
                    if (rij2 + 1e-6 < d*d) { // d=1.0
                        overlap = true;
                        break;
                    }
                }
                if (overlap) continue;

                // ---- Accept ----
                X.push_back(Xtry);
                Y.push_back(Ytry);
                Z.push_back(Ztry);
                ++IC;

                std::cout << "IC " << IC << "  " << Xtry << " " << Ytry << " " << Ztry << "\n";

                if (IC == NTOT) goto done;
            }
        }
    }

done:
    if (IC != NTOT) {
        std::cerr << "NOT ENOUGH PARTICLES\n";
        return 1;
    }

// added code for INPUT file

 const double AREA = 6.0 * (AL + 1.0) * (AL + 1.0);

    std::ofstream input("INPUT");
    if (!input) {
        std::cerr << "ERROR: cannot create INPUT file\n";
        return 1;
    }

    input << 1000000 << "\n";

    input << std::fixed << std::setprecision(16);
    
    input <<std::fixed << std::setprecision(16)<< AL << " " << AL << " " << AL << "\n";
    input << NTOT << "\n";
    input << AREA << std::fixed << std::setprecision(14)<<"\n";

    input.close();



    // ---- Write POS.TXT (like Fortran) ----
    std::ofstream pos("POS.TXT");
    if (!pos) {
        std::cerr << "ERROR: cannot create POS.TXT\n";
        return 1;
    }

   pos << std::fixed << std::setprecision(16);

    for (int i = 0; i < NTOT; ++i) {

        
        pos << std::setw(22)<< X[i] <<std::setw(22) << Y[i] <<std::setw(22) << Z[i] << "\n";
    }
    pos.close();

    return 0;
}
