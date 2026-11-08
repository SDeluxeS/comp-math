#include <vector>

struct ShGrav {
    std::vector<double> cilm;
    int cilm_dim;
    int lmax;
    double gm;
    double r0;
};

inline int shGravIndex(int i, int l, int m, int cilm_dim) {
    return i + 2*l + 2*cilm_dim*m;
}

ShGrav shGravCreate(int lmax, double gm, double r0);
int shGravLoadCoefficients(ShGrav &sh, const char *filename);

void shGravCalculate(const ShGrav &sh,
                     double x, double y, double z,
                     double result[3]);
void shGravGetPartialsMatrix(const ShGrav &sh,
                             double x, double y, double z,
                             double M[9]);

int shGravGetCoefficientC(const ShGrav &sh, int l, int m, double *c);
int shGravGetCoefficientS(const ShGrav &sh, int l, int m, double *s);
int shGravSetCoefficientC(ShGrav &sh, int l, int m, double c);
int shGravSetCoefficientS(ShGrav &sh, int l, int m, double s);
