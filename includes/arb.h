#ifndef ARB_H
#define ARB_H




#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#ifndef _CRT_NONSTDC_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#endif


#ifdef __GNUC__
#include   <dirent.h>
#include   <unistd.h>
#elif defined(_WIN32)
#include   <dirent_win.h>
#elif defined(_WIN64)
#include   <dirent_win.h>
#endif

#include <inttypes.h>
#include <stdio.h>

#if defined(_WIN32) || defined(_WIN64)
#define ARB_EXPORTS 1
#ifdef ARB_EXPORTS
#define ARB_API __declspec(dllexport)
#else
#define ARB_API __declspec(dllimport)
#endif
#else
#define ARB_API
#endif

#ifndef NULL
#define NULL   ((void *) 0)
#endif

#define IM1     2147483563
#define IM2     2147483399
#define AM      (1.0/IM1)
#define IMM1    (IM1-1.0)
#define IA1     40014
#define IA2     40692
#define INFD    HUGE_VAL//1.0E308
#define INFL    HUGE_VALL//2147483647
#define IQ1     53668
#define IQ2     52774
#define IR1     12211
#define IR2     3791
#define LN2     0.693147180559945309417232121458176
#define LN10    2.302585092994045684017991454684364
#define LNPI    1.144729885849400174143427351353058
#define NTAB    32
#define NDIV    (1.0+IMM1/NTAB)
#define EPS     1.2e-7
#define RNMX    (1.0-EPS)
#define PI      3.141592653589793238462643383279502
#define EULER   2.718281828459045235360287471352662
#define M_EULER 0.57721566490153286060651209008
#define ROOT2   1.414213562373095048801688724209698
#define ROOT3   1.732050807568877293527446341505872
#define ROOT5   2.236067977499789696409173668731276
#define STREX0  "\?\'\"!#$%&()*<=>@[]^_`{}~ABCFGHIJKLMNOPQRSTUVWXYZabcfghijklmn\
opqrstuvwxyz"
#define STREX00 "\?\'\"!#$%&()*<=>@[]^_`{}~ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghij\
klmnopqrstuvwxyz"
#define STREXN  "\?\'\"!#$%&()*<=>@[]^_`{}~:+-,.;\\/ |ABCDEFGHIJKLMNOPQRSTUVWXY\
Zabcdefghijklmnopqrstuvwxyz"
#define STREX1  "\n\t\v\b\r\f\a\\\?\'\" !#$%&()+-*/,:;<=>@[]^_`{|}~."
#define STREX2  "\n\t\v\b\r\f\a "
#define STREX3  "\n\t\v\b\r\f\a\\ ,/:;|"
#define STREX4  "0123456789"
#define STREX5  "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
#define MSUB    7
#define NSTACK  500
#define MAXDEG  50
#define DIM_MAX  20
#define NBITS   31
#define MAX_DIGITS 15
#define EPSILON 1.0E-15

typedef enum { memory_stream, file_stream }stream_type;
typedef struct cmplx cmplx;
typedef struct MEM_FILE MEM_FILE;
typedef struct bit_array bit_array;

typedef struct cmplx {

  double r;
  double i;
}cmplx;

ARB_API double          AbsC(cmplx);
ARB_API double        **AbsCM(cmplx **);
ARB_API double         *AbsCV(cmplx *);
ARB_API cmplx           AddC(cmplx,cmplx);
ARB_API cmplx           AddNC(cmplx *);
ARB_API double         *AddV(double *,double *);
ARB_API cmplx          *AddVC(cmplx *,cmplx *);
ARB_API void            AddIPM(double **,double **,double **);
ARB_API double        **AddM(double **,double **);
ARB_API cmplx         **AddMC(cmplx **,cmplx **);
ARB_API void            AddPrimeFactors(double,double ,bit_array *,long *,
                                double **);
ARB_API cmplx           Add3C(cmplx,cmplx,cmplx);
ARB_API void            AdjacentAveraging(double *,long);
ARB_API void            AdjacentAveragingC(cmplx *,long);
ARB_API void            AdjacentAveragingML(long **,long,long);
ARB_API void            AdjacentAveragingInterpolation(double *,long);
ARB_API void            AdjacentAveragingInterpolationC(cmplx *,long);
ARB_API void            AdjacentInterpolation(double *,long);
ARB_API void            AdjacentInterpolationC(cmplx *,long);
ARB_API cmplx         **Adjoint(cmplx **);
ARB_API void            AdjointIP(cmplx **,cmplx **);
ARB_API void            AdjointSquareIP(cmplx **);
ARB_API void            AiryFuncSub(double,double *,double *,double *,double *);
ARB_API long            AlphabeticStringComparator(char *,char *);
ARB_API double         *AllocDouble(double);
ARB_API long           *AllocLong(long);
ARB_API void            Amoeba(double (*)(double *),double **,double *,double,long *);
ARB_API double          AmoebaTry(double **,double *,double *,double (*)(double *),
                          long,double,double *);
ARB_API cmplx         **AnalyticSignal(double **);
ARB_API void            AppendFilenameR(char **,char *);
ARB_API void            AppendMatrix(double ****,double **);
ARB_API void            AppendMatrixVector(double ****,double ***);
ARB_API char           *AppendS(char **, char *);
ARB_API char          **AppendString(char ***, char *);
ARB_API cmplx           ArccosC(cmplx);
ARB_API cmplx           ArccoshC(cmplx);
ARB_API cmplx           ArcsinC(cmplx);
ARB_API cmplx           ArcsinhC(cmplx);
ARB_API cmplx           ArctanC(cmplx);
ARB_API cmplx           ArctanhC(cmplx);
ARB_API void            AugmentMatrixHC(cmplx ***,cmplx **);
ARB_API void            AugmentMatrixHL(long ***,long **);
ARB_API void            AugmentMatrixV(double ***,double **);
ARB_API void            AugmentMatrixVC(cmplx ***,cmplx **);
ARB_API void            AugmentMatrixVL(long ***,long **);
ARB_API char           *AugmentString(char **,long);
ARB_API void            AugmentVector(double **,double *);
ARB_API void            AugmentVectorC(cmplx **,cmplx *);
ARB_API void            AugmentVectorL(long **,long *);
ARB_API double        **Autocorrelation(double **);
ARB_API double        **AutoSample1D(double (*)(double),double,double,double);
ARB_API cmplx         **AutoSample1DC(cmplx (*)(double),double,double,double);
ARB_API long           *BanDec(double **,long,long,double **,double *);
ARB_API long           *BanDecC(cmplx **,long,long,cmplx **,double *);
ARB_API void            BanSolve(double **,long,long,double *);
ARB_API void            BanSolveC(cmplx **,long,long,cmplx *);
ARB_API unsigned char  *Base64Decoding(char *,long *);
ARB_API char           *Base64Encoding(unsigned char *, long);
ARB_API void            Base64EncodingTriplet(uint8_t,uint8_t,uint8_t,char *,char *);
ARB_API double          BernoulliNumber(double, char *abort);
ARB_API double         *BesselFuncRoots(double,double,long);
ARB_API double          BesselFuncSeries(double,double);
ARB_API long            BesselFuncSub(double,double,double *,double *,double *,
                              double *);
ARB_API void            BesselGammas(double,double *,double *,double *,double *);
ARB_API cmplx           BesselIFunc(double,double);
ARB_API cmplx           BesselJFunc(double,double);
ARB_API cmplx           BesselKFunc(double,double);
ARB_API long            BesselModFuncSub(double,double,double *,double *,double *,
                                 double *);
ARB_API cmplx           BesselYFunc(double,double);
ARB_API double          BetaFuncCF(double,double,double);
ARB_API double          BetaFuncI(double,double,double);
ARB_API void            BicubicCoefficients(double *,double *,double *,double *,double,
                                    double,double **,double **);
ARB_API double       ***BicubicInterpolation(double ***,long,long);
ARB_API double        **BicubicInterpolationM(double **,long,long);
ARB_API void            BicubicInterpolationSquare(double *,double *,double *,double *,
                                           double,double,double,double,double,
                                           double,double *,double *,double *,
                                           double **);
ARB_API long            BinarySearchPos(double,double *);
ARB_API long            BinarySearchPosC(double,cmplx *);
ARB_API long            BinarySearchPosS(char *, char **);
ARB_API void            BinToAscii(char *);
ARB_API void            BinToAsciiB(char *);
ARB_API void            BinToAsciiI(char *);
ARB_API void            BinToAsciiL(char *);
ARB_API long            Bisection(double (*)(double),double *,double,double,double,
                          long);
ARB_API double          BisectionN(double (*)(double *),double *,long,double,double,
                           double,long);
ARB_API double        **BlockDiagonalMatrix(double ***);
ARB_API double         *BoxMuller(long, long *);
ARB_API double          BoxMuller1(long *,long *,long *,long *);
ARB_API double         *BoxMuller2(long,double *);
ARB_API void            BracketMinimum(double (*)(double),double,double,double *,
                               double *,double *,double *,double *,double *);
ARB_API double          BrentMinimization(double (*)(double),double,double,double);
ARB_API char           *ByteArrayToString(unsigned char *,long);
ARB_API long           *ByteArrayToVectorL(unsigned char *,long);
ARB_API uint32_t       *ByteArrayToVectorUint32(unsigned char *, long);
ARB_API double        **CartesianProduct(double *, double *,char *);
ARB_API double          Ceiling(double);
ARB_API long            CenterDistr(double *);
ARB_API long            CenterDistrC(cmplx *);
ARB_API double          CenterDistribution(double **);
ARB_API double          CenterDistributionC(cmplx **);
ARB_API double          CenterDistributionWindowC(cmplx **,long,long,double);
ARB_API long            CenterDistrWindowC(cmplx *,long,long);
ARB_API long            CenterDistr0C(cmplx *);
ARB_API void            CenterDomain(double *,double);
ARB_API void            CenterDomainC(cmplx *,double);
ARB_API double          CenterDomainFFT(double *);
ARB_API double          CenterDomainFFTC(cmplx *);
ARB_API void            ChangeExtension(char **,char *);
ARB_API char           *CharString(char);
ARB_API char           *CharToHex(char);
ARB_API double          ChebychevExpansion(double,double,double *,long,double);
ARB_API double          ChiSquare(double **,double (*)(double,double *),double *);
ARB_API cmplx           CI(void);
ARB_API double         *Col(double **,long);
ARB_API double          Comb(long,long);
ARB_API double          CombFunc(double,double);
ARB_API double        **CombineMatrices(double **,double **);
ARB_API long          **CombineMatricesL(long **,long **);
ARB_API long            CompareFloats(double,double);
ARB_API double        **ComplementaryM(double **);
ARB_API cmplx           ComplexP(double,double);
ARB_API cmplx           ComplexR(double,double);
ARB_API cmplx          *ComplexRoots(cmplx,long);
ARB_API char          **ComplexString(cmplx,long);
ARB_API char           *Concat(char *, char *);
ARB_API char           *ConcatFilenameL(char *,char *);
ARB_API char           *ConcatFilenameR(char *,char *);
ARB_API char           *ConcatN(char **);
ARB_API char          **ConcatStringVector(char **,char **);
ARB_API char           *Concat3(char *,char *,char *);
ARB_API char           *Concat4(char *,char *,char *,char *);
ARB_API char           *Concat5(char *,char *,char *,char *,char *);
ARB_API cmplx           ConjC(cmplx);
ARB_API double        **ComplementM(double **,double **);
ARB_API double         *ComplementV(double *,double *);
ARB_API double       ***ConvertMVVM(double ***,long);
ARB_API cmplx        ***ConvertMVVMC(cmplx ***,long);
ARB_API double       ***ConvertVMMV(double ***);
ARB_API cmplx        ***ConvertVMMVC(cmplx ***);
ARB_API double        **ConvolutionC(double (*)(double),double (*)(double),double,
                             double,long,double,double,long);
ARB_API double        **ConvolutionD(double **,double **);
ARB_API cmplx         **ConvolutionFFT(cmplx **,cmplx **);
ARB_API cmplx          *ConvolutionFFTV(cmplx *,cmplx *);
ARB_API double        **Convolution1D(double *,double *,double *);
ARB_API cmplx           CopyC(cmplx);
ARB_API double        **CopyMatrix(double **);
ARB_API cmplx         **CopyMatrixC(cmplx **);
ARB_API void            CopyMatrixIPC(cmplx **,cmplx **);
ARB_API int           **CopyMatrixI(int **);
ARB_API void            CopyMatrixIP(double **,double **);
ARB_API void            CopyMatrixIPL(long **,long **);
ARB_API float         **CopyMatrixF(float **);
ARB_API long          **CopyMatrixL(long **);
ARB_API uint32_t      **CopyMatrixUInt32(uint32_t **);
ARB_API double       ***CopyMatrixVector(double ***);
ARB_API char           *CopyS(char *);
ARB_API char         ***CopyStringMatrix(char ***);
ARB_API char        ****CopyStringMatrixVector(char ****);
ARB_API char          **CopyStringVector(char **);
ARB_API cmplx        ***CopyMatrixVectorC(cmplx ***);
ARB_API void           *CopyMemFile(void *v);
ARB_API double         *CopyVector(double *);
ARB_API cmplx          *CopyVectorC(cmplx *);
ARB_API int            *CopyVectorI(int *);
ARB_API void            CopyVectorIP(double *,double *);
ARB_API void            CopyVectorIPC(cmplx *,cmplx *);
ARB_API void            CopyVectorIPL(long *,long *);
ARB_API float          *CopyVectorF(float *);
ARB_API long           *CopyVectorL(long *);
ARB_API double          CorrelationCoefficient(double *,double *);
ARB_API cmplx           CosC(cmplx);
ARB_API cmplx           CoshC(cmplx);
ARB_API double          Cos2PulseFunc(double);
ARB_API double          Cos4PulseFunc(double);
ARB_API double          Cos8PulseFunc(double);
ARB_API long            CoulombWaveFunc(double, long, double *, double *, double *,
                                long, double, long);
ARB_API double          CoulombWaveFuncBess(double,double,long,long);
ARB_API void            CoulombWaveFuncBesselSph(double,double,long,double *,long);
ARB_API void            CoulombWaveFuncCF(double,double,long,double *,double *,
                                  double *,double *,long,double);
ARB_API long            CoulombWaveFuncF(double,long,double *,double *,long,double,
                                 long);
ARB_API long            CoulombWaveFuncFPF(double,double,long,double *,long);
ARB_API long            CoulombWaveFuncG(double,long,double *,double *,long,double,
                                 double,long);
ARB_API double          CoulombWaveFuncFPrimeSeries(double,double,long,long);
ARB_API double          CoulombWaveFuncFSeries(double,double,long,long);
ARB_API void            CoulombWaveFuncHPH(double,double,long,double *,double *,double);
ARB_API void            CoulombWaveFuncIntegrateFABM(double,long,double *,double *,
                                             long,long,double,double,long);
ARB_API void            CoulombWaveFuncIntegrateF0ABM(double,double *,double *,long,
                                              long,double,double,long);
ARB_API void            CoulombWaveFuncIntegrateFPrimeABM(double *,double *,double *,
                                                  long,long);
ARB_API void            CoulombWaveFuncIntegrateGABM(double,long,double *,double *,
                                             long,long,double,double,long);
ARB_API void            CoulombWaveFuncIntegrateG0ABM(double,double *,double *,long,
                                              long,double,double,long);
ARB_API void            CoulombWaveFuncIntegrateGPrimeABM(double *,double *,double *,
                                                  long,long);
ARB_API void            CoulombWaveFuncNear(double *,double *,double *,double,double);
ARB_API double          Covariance(double *,double *);
ARB_API double        **CropMatrix(double **,long,long,long,long);
ARB_API cmplx         **CropMatrixC(cmplx **,long,long,long,long);
ARB_API double        **CropMatrixCentered(double **,long,long);
ARB_API cmplx         **CropMatrixCenteredC(cmplx **,long,long);
ARB_API char           *CropString(char *, long, long);
ARB_API char          **CropStringVector(char **,long,long);
ARB_API double         *CropVector(double *,long,long);
ARB_API cmplx          *CropVectorC(cmplx *,long,long);
ARB_API double         *CropVectorCentered(double *,long);
ARB_API cmplx          *CropVectorCenteredC(cmplx *,long);
ARB_API double        **Crop1D(double **,double,double);
ARB_API double       ***Crop2D(double ***,double,double,double,double);
ARB_API double         *CrossP(double *,double *);
ARB_API double          CubicInterpolation1(double *,double *,double);
ARB_API double         *Curl(double *(*)(double *),double *,double);
ARB_API double         *CurlCyl(double *(*)(double *),double *,double);
ARB_API double         *CurlSph(double *(*)(double *),double *,double);
ARB_API cmplx           CZ(void);
ARB_API cmplx           C1(void);
ARB_API void            DampenLeft(double *,double *,double,double);
ARB_API void            DampenLeftFit(double *,long,long);
ARB_API void            DampenLeftGauss(double *,long,double);
ARB_API void            DampenLeftLinear(double *,long);
ARB_API void            DampenLeftZero(double *,long,double);
ARB_API void            DampenRight(double *,double *,double,double);
ARB_API void            DampenRightExp(double *,long,double);
ARB_API void            DampenRightFit(double *,long,long);
ARB_API void            DampenRightGauss(double *,long,double);
ARB_API void            DampenRightLinear(double *,long);
ARB_API double          DeltaFunc(double);
ARB_API double          Derivative(double (*)(double),double,double);
ARB_API double        **DerivativeDCyclic(double **);
ARB_API double          DerivativeN(double (*)(double),double,long,double);
ARB_API double          DerivativeP(double (*)(double *),double *,long,double);
ARB_API double         *DerivativePV(double *(*)(double *),double *,long,double);
ARB_API double        **DerivativeSpectral(double *,double *,double,double);
ARB_API double         *DerivativeV(double *(*)(double),double,double);
ARB_API double          DerivativeWeights(double *,double *,double);
ARB_API double        **Derivative1D(double *,double *);
ARB_API cmplx         **Derivative1DC(cmplx **);
ARB_API double          Derivative2(double (*)(double),double,double);
ARB_API double        **Derivative2D(double **);
ARB_API double        **Derivative2DCyclic(double **);
ARB_API double          Derivative2P(double (*)(double *),double *,long,double);
ARB_API double         *Derivative2PV(double *(*)(double *),double *,long,double);
ARB_API double         *Derivative2V(double *(*)(double),double,double);
ARB_API double          Det(double **);
ARB_API cmplx           DetC(cmplx **);
ARB_API double        **DFT(double **,double,double,long,int);
ARB_API double        **DiagonalMatrix(double *);
ARB_API cmplx         **DiagonalMatrixC(cmplx *);
ARB_API cmplx           DigammaFunc(cmplx, char *);
ARB_API cmplx         **Discretize1DC(cmplx (*)(double),double,double,long);
ARB_API double        **DiscretizeAdaptive(double (*)(double),double,double,double);
ARB_API void            DiscretizeAdaptiveSub(double (*)(double),double ,double,double,
                                      double,double,double,double,double ***,
                                      long *);
ARB_API cmplx         **DiscretizeC(cmplx (*)(cmplx),cmplx,cmplx,long);
ARB_API cmplx         **DiscretizeCC(cmplx (*)(cmplx),cmplx,cmplx,long,long);
ARB_API double        **DiscretizeDFT(cmplx (*)(double),double,double,long);
ARB_API double         *DiscretizeDomain(double,double,long);
ARB_API cmplx          *DiscretizeDomainC(double,double,long);
ARB_API double         *DiscretizeDomainFFT(double,double,long);
ARB_API cmplx          *DiscretizeDomainFFTC(double,double,long);
ARB_API void            DiscretizeDomainIP(double *,double,double);
ARB_API double        **DiscretizeDomainN(double *,double *,long *);
ARB_API cmplx         **DiscretizeFFT(cmplx (*)(double),double,double,long);
ARB_API double        **DiscretizeFFTN(cmplx (*)(double *),double *,double *,long *);
ARB_API double        **DiscretizeN(double (*)(double *),double *,double *,long *);
ARB_API double        **DiscretizeNU(double (*)(double),double,double,double *);
ARB_API cmplx         **DiscretizeNUC(cmplx (*)(double),double,double,double *);
ARB_API double        **DiscretizeN1(double (*)(double *),double,double,long,long,
                             double *);
ARB_API double        **DiscretizeP(double (*)(double,double *),double,double,long,
                            double *);
ARB_API cmplx         **DiscretizePC(cmplx (*)(double,double *),double,double,long,
                             double *);
ARB_API double        **DiscretizePN(double (*)(double *,double *),double *,double *,
                             long *,double *);
ARB_API cmplx         **DiscretizePNC(cmplx (*)(double *,double *),double *,double *,
                              long *,double *);
ARB_API double        **DiscretizeRungeKutta1D(void (*)(double,double *,double *),
                                       double,double);
ARB_API cmplx         **DiscretizeRungeKutta1DC(void (*)(double,cmplx *,cmplx *),
                                        double,double);
ARB_API double       ***DiscretizeSphericalCS3D(double (*)(double,double),double *,
                                        double *);
ARB_API double       ***DiscretizeSpherical3D(double (*)(double,double,double),
                                      double *,double *,double *);
ARB_API double        **DiscretizeV(double (*)(double *),double,double,long,long,
                            double *);
ARB_API double        **Discretize1D(double (*)(double),double,double,long);
ARB_API double       ***Discretize2D(double (*)(double,double),double,double,double,
                             double,long,long);
ARB_API double          Distance(double *,double *);
ARB_API double          Div(double *(*)(double *),double *,double);
ARB_API cmplx           DivC(cmplx,cmplx);
ARB_API cmplx           DivCR(cmplx,double);
ARB_API double          DivCyl(double *(*)(double *),double *,double);
ARB_API long            Divisible(double,double);
ARB_API double         *Divisors(double,char *);
ARB_API double          DivSph(double *(*)(double *),double *,double);
ARB_API double          DotP(double *,double *);
ARB_API char           *DoubleString(double,long);
ARB_API char           *DoubleStringCleanZeros(char **);
ARB_API char           *DoubleStringDecimal(double,long);
ARB_API char           *DoubleStringExponential(double,long);
ARB_API char           *DoubleStringRound(char **,long);
ARB_API double        **DownSampleCols(double **,long);
ARB_API cmplx         **DownSampleColsC(cmplx **,long);
ARB_API double         *DownSampleV(double *,long);
ARB_API cmplx          *DownSampleVC(cmplx *,long);
ARB_API void            EigenDecomposition(double **,double *,double **);
ARB_API cmplx          *EigenVector(double **mat,cmplx,long,long *);
ARB_API long          **EnumerateCombinations(long, long, char *);
ARB_API void            EnumerateCombinationsSub(long **, long, long, long, long, long,
                                         char *);
ARB_API long          **EnumeratePermutations(long, char *);
ARB_API void            EnumeratePermutationsSub(long **, long, long, int *, char *);
ARB_API long          **EnumerateTuples(long *);
ARB_API void            EnumerateTuplesSub(long **,long,long,long *);
ARB_API double         *Envelope(double *,long);
ARB_API int             EqualMatrices(double **,double **);
ARB_API int             EqualMatricesC(cmplx **,cmplx **);
ARB_API int             EqualMatricesF(float **,float **);
ARB_API int             EqualMatricesI(int **,int **);
ARB_API int             EqualMatricesL(long **,long **);
ARB_API int             EqualMatricesS(char ***,char ***);
ARB_API long            EqualRows(double **,long,double **,long);
ARB_API int             EqualVectors(double *,double *);
ARB_API int             EqualVectorsC(cmplx *,cmplx *);
ARB_API int             EqualVectorsS(char **,char **);
ARB_API int             EqualVectorsF(float *,float *);
ARB_API int             EqualVectorsI(int *,int *);
ARB_API int             EqualVectorsL(long *,long *);
ARB_API double          EuclideanDistM(double **,double **);
ARB_API double          EuclideanDistV(double *,double *);
ARB_API double          EuclideanNormM(double **);
ARB_API double          EuclideanNormV(double *);
ARB_API double          EuclideanNorm2(double,double);
ARB_API double          EuclideanNorm3(double,double,double);
ARB_API int             Even(double);
ARB_API cmplx           ExpC(double);
ARB_API cmplx           ExpCC(cmplx);
ARB_API void            ExtendVector(double **,long);
ARB_API void            ExtendVectorC(cmplx **,long);
ARB_API void            ExtendVectorF(float **,long);
ARB_API void            ExtendVectorI(int **,long);
ARB_API void            ExtendVectorL(long **,long);
ARB_API void            ExtendMatrix(double ***,long,long);
ARB_API void            ExtendMatrixC(cmplx ***,long,long);
ARB_API void            ExtendMatrixF(float ***,long,long);
ARB_API void            ExtendMatrixI(int ***,long,long);
ARB_API void            ExtendMatrixL(long ***,long,long);
ARB_API double          Factorial(long);
ARB_API double          FactorialLn(long);
ARB_API double          FactorialProduct(long *);
ARB_API double          FactorialProductLn(long *);
ARB_API double          FactorialProductRatio(long *,long *);
ARB_API double          FactorialProductRatioLn(long *,long *);
ARB_API double          FactorialProduct3(long,long,long);
ARB_API double          FactorialRatio(long,long);
ARB_API double          FactorialRatioLn(long,long);
ARB_API double          FactorialZeros(double,char *);
ARB_API void            FFTCols(cmplx **,long);
ARB_API void            FFTDomainChange(double *);
ARB_API void            FFTDomainChangeC(cmplx *);
ARB_API void            FFTM(cmplx **,long);
ARB_API void            FFTND(double **,long);
ARB_API void            FFTNDomainChange(double **);
ARB_API long           *FFTNGetNPts(double **);
ARB_API void            FFTNShift(double **);
ARB_API void            FFTRows(cmplx **,long);
ARB_API void            FFTV(cmplx *,long);
ARB_API void            FFTVErr(cmplx *,cmplx *,long);
ARB_API void            FFT1D(cmplx **,long);
ARB_API char           *FGetS(char *,long *,FILE *);
ARB_API char           *Filename(char *);
ARB_API double         *FindPeaks1D(double *,double *);
ARB_API double         *FindPeaks1DC(cmplx *,double *);
ARB_API double         *FindRoots1D(double *,double *,double);
ARB_API double         *FindRoots1DC(cmplx *,double *,double);
ARB_API long           *FindS(char *, char *);
ARB_API double         *FindSimpleRoots1D(double *,double *,double);
ARB_API double         *FindSimpleRoots1DC(cmplx *,double *,double);
ARB_API double          Floor(double);
ARB_API cmplx         **FormatFFT(cmplx **);
ARB_API cmplx         **FormatZeroPadFFT(cmplx **,long);
ARB_API void            FourierFilter(double **,double,double,double);
ARB_API void            FourierFilterC(cmplx **,double,double,double);
ARB_API double          FractionalPart(double);
ARB_API double          FredholmInterpolation1(double **,double,double (*)(double),
                                       double (*)(double,double));
ARB_API cmplx           FredholmInterpolation1C(cmplx **,double,cmplx (*)(double),
                                        cmplx (*)(double,double));
ARB_API double        **Fredholm2(double (*)(double),double (*)(double,double),long,
                          long,long,double);
ARB_API cmplx         **Fredholm21DC(cmplx (*)(double),cmplx (*)(double,double),
                             double,double,long,double);
ARB_API void            FreeMatrixUInt32(uint32_t **t);
ARB_API void            Free(void *);
ARB_API void            FreeMatrix(double **);
ARB_API void            FreeMatrixC(cmplx **);
ARB_API void            FreeMatrixF(float **);
ARB_API void            FreeMatrixI(int **);
ARB_API void            FreeMatrixL(long **);
ARB_API void            FreeMatrixVector(double ***);
ARB_API void            FreeMatrixVectorC(cmplx ***);
ARB_API void            FreeMatrixVectorL(long ***);
ARB_API void            FreeString(char *s);
ARB_API void            FreeStringMatrix(char ***);
ARB_API void            FreeStringMatrixVector(char ****);
ARB_API void            FreeStringVector(char **);
ARB_API void            FreeVector(double *);
ARB_API void            FreeVectorC(cmplx *);
ARB_API void            FreeVectorL(long *);
ARB_API void            FreeVectorMatrix(double ***);
ARB_API void            FreeVectorMatrixC(cmplx ***);
ARB_API double        **FSeries(double (*)(double),double,double,long,long);
ARB_API double          FWHM(double **,long,long);
ARB_API double          FWHMC(cmplx **,long,long);
ARB_API double          FWHMIntensityC(cmplx **,long,long);
ARB_API double          FWHMWindow(double **, long, long, long);
ARB_API double          GammaLn(double);
ARB_API cmplx           GammaLnC(cmplx);
ARB_API void            GaussianFit(char *,double **);
ARB_API double          GaussianFitFunction(double,double *);
ARB_API double        **GaussLaguerre(long,double);
ARB_API double        **GaussLegendre(double,double,long,double);
ARB_API long            GreatestCommonDivisor(long,long);
ARB_API double        **GenerateComb(double,double,double,long);
ARB_API double        **GenerateCombinations(double *, long, long,char *);
ARB_API long          **GenerateCombinationsL(long *, long, long,char *);
ARB_API double        **GenerateRowCombination(double **,long **,long);
ARB_API double        **GenerateRowPermutation(double **,long **,long);
ARB_API double        **GeneratePermutations(double *,char *);
ARB_API long          **GeneratePermutationsL(long *,char *);
ARB_API char           *GetDirectory(char *);
ARB_API double          GetMin(double *,double *);
ARB_API long           *GetNPts(double **);
ARB_API double         *Grad(double (*)(double *),double *,double);
ARB_API double         *GradCyl(double (*)(double *),double *,double);
ARB_API double         *Gradient1D(double *,double);
ARB_API cmplx          *Gradient1DC(cmplx *,double);
ARB_API double       ***Gradient2D(double **,double,double);
ARB_API double         *GradSph(double (*)(double *),double *,double);
ARB_API void            GramSchmidt(double **);
ARB_API cmplx         **GramSchmidtC(cmplx **);
ARB_API double          HalfCircleFunc(double);
ARB_API double          HannFunc(double,double,double);
ARB_API char           *HashFromByteArray(MEM_FILE *,long,long);
ARB_API char           *HashFromString(char *,long);
ARB_API long            HashFromStringL(char *);
ARB_API char           *HashWithSalt(char *,char *);
ARB_API double          HeavisideFunc(double);
ARB_API uint8_t         HexDigitToDecimal(char);
ARB_API uint8_t         HexByteToDecimal(char *);
ARB_API double        **Histogram(double *,long);
ARB_API long          **HistogramD(long *);
ARB_API double         *HMSTime(double,double,double);
ARB_API void            HouseholderReduction(double **,double *,double *);
ARB_API uint32_t        HSVToLong(double,double,double);
ARB_API void            HSVToRGB(double *,double *,double *,double,double,double);
ARB_API cmplx           HydrogenWaveFunc(double,double,double,long,long,long);
ARB_API double        **IdentityMatrix(long);
ARB_API double          Im(cmplx);
ARB_API double        **ImM(cmplx **);
ARB_API double         *ImV(cmplx *);
ARB_API cmplx           ImagC(double);
ARB_API void            IncM(double **,double **);
ARB_API void            IncV(double *,double *);
ARB_API double          InnerP(double *,double *);
ARB_API cmplx           InnerPC(cmplx *,cmplx *);
ARB_API char           *InsertChar(char,char **,long);
ARB_API void            InsertCol(double *,double ***,long);
ARB_API void            InsertColC(cmplx *,cmplx ***,long);
ARB_API void            InsertElement(double,double **,long);
ARB_API void            InsertElementC(cmplx,cmplx **,long);
ARB_API void            InsertElementL(long,long **,long);
ARB_API char          **InsertElementS(char *,char ***,long);
ARB_API void            InsertElementSM(char ***,char *****,long);
ARB_API void            InsertElementSV(char **,char ****,long);
ARB_API void            InsertElementUint32(uint32_t,uint32_t **,long);
ARB_API void            InsertMatrix(double **,double ****,long);
ARB_API void            InsertRow(double *,double ***,long);
ARB_API void            InsertRowC(cmplx *,cmplx ***,long);
ARB_API void            InsertRowS(char **,char ****,long);
ARB_API char           *InsertStr(char *, char **, long);
ARB_API int             InsideInterval(double,double,double);
ARB_API double          Integral(double (*)(double),double,double,long);
ARB_API double          Integral2(double (*)(double,double),double,double,double,
                          double,long,long);
ARB_API double          Integral2D(double **);
ARB_API double          IntegralD(double **,long,long);
ARB_API double          IntegralP(double (*)(double,double *),double,double,long,
                          double *);
ARB_API double          IntegralPath(double *(*)(double *),double *(*)(double),double,
                             double,long);
ARB_API cmplx           IntegralPC(cmplx (*)(double,double *),double,double,long,
                           double *);
ARB_API double          IntegralSurfaceS(double (*)(double *),double *(*)(double *),
                                 double *,double *,long *);
ARB_API double          IntegralSurfaceV(double *(*)(double *),double *(*)(double *),
                                 double *,double *,long *);
ARB_API double         *IntegralV(double *(*)(double),double,double,long);
ARB_API void            Integrator(double *,double *,double *,long);
ARB_API void            IntegratorC(cmplx *,double *,cmplx *,long);
ARB_API double         *IntersectionV(double *,double *);
ARB_API double        **IntersectionM(double **,double **);
ARB_API char          **IntersectionSV(char ***);
ARB_API void            InvertCols(double **);
ARB_API void            InvertColsC(cmplx **);
ARB_API void            InvertM(double **);
ARB_API void            InvertMC(cmplx **);
ARB_API void            InvertRows(double **);
ARB_API void            InvertRowsC(cmplx **);
ARB_API void            Invert2M(double **);
ARB_API char           *InvertCase(char *in);
ARB_API cmplx           InvC(cmplx);
ARB_API long            IsNAN(double);
ARB_API long            IsINF(double);
ARB_API double        **JacobianMatrix(double *(*)(double *),double *,double);
ARB_API double        **JoinMatrixH(double **,double **);
ARB_API long          **JoinMatrixHL(long **,long **);
ARB_API double       ***JoinMatrixVector(double ***,double ***);
ARB_API double        **JoinMatrixV(double **,double **);
ARB_API long          **JoinMatrixVL(long **,long **);
ARB_API double         *JoinVector(double *,double *);
ARB_API long           *JoinVectorL(long *,long *);
ARB_API double        **JoinVectorM(double *,double *);
ARB_API long            Kronecker(long,long);
ARB_API double        **Lagrange(double **,double,double,long);
ARB_API double          LaguerreFunc(double,long,double);
ARB_API double          Lapl(double (*)(double *),double *,double);
ARB_API double         *Laplacian1D(double *,double);
ARB_API cmplx          *Laplacian1DC(cmplx *,double);
ARB_API double        **Laplacian2D(double **,double,double);
ARB_API double          LaplCyl(double (*)(double *),double *,double);
ARB_API double          LaplSph(double (*)(double *),double *,double);
ARB_API char           *LeadingCharacters(char **,char,char,long);
ARB_API uint32_t        LeftRotate32(uint32_t,int32_t);
ARB_API uint64_t        LeftRotate64(uint64_t,int64_t);
ARB_API double          LegendreFunc(double,long,long);
ARB_API double          LegendrePolynomial(double,long);
ARB_API double        **LinearInterpolation1D(double **,double,double,long);
ARB_API cmplx         **LinearInterpolation1DC(cmplx **,double,double,long);
ARB_API double        **LinearInterpolationCols(double **,long);
ARB_API cmplx           LinearInterpolation1C(cmplx *,double *,double);
ARB_API double        **LinearInterpolationM(double **,long,long);
ARB_API long          **LinearInterpolationML(long **,long,long);
ARB_API void            LinearInterpolationPts1D(double *,double *,double *,double *);
ARB_API void            LinearInterpolationPts1DC(cmplx *,double *,cmplx *,double *);
ARB_API void            LinearInterpolationPts2D(double **,double *,double *,double **,
                                         double *,double *);
ARB_API double        **LinearInterpolationRows(double **,long);
ARB_API double         *LinearInterpolationV(double *,long);
ARB_API cmplx          *LinearInterpolationVC(cmplx *,long);
ARB_API double          LinearInterpolation1(double *,double *,double);
ARB_API double       ***LinearInterpolation2D(double ***,long,long);
ARB_API double          LinearInterpolation2D1(double **,double *,double *,double,
                                       double);
ARB_API double          LinearPhase(cmplx *);
ARB_API void            LinSolve(double **,double **);
ARB_API unsigned char  *LoadFileToByteArray(char *,long *);
ARB_API char           *LoadFileToString(char *);
ARB_API double        **LoadMatrix(char *);
ARB_API double        **LoadMatrixBin(char *);
ARB_API cmplx         **LoadMatrixC(char *);
ARB_API cmplx         **LoadMatrixCBin(char *);
ARB_API long          **LoadMatrixL(char *);
ARB_API long          **LoadMatrixLBin(char *);
ARB_API double       ***LoadMatrixVector(char *);
ARB_API double       ***LoadMatrixVectorBin(char *);
ARB_API long            LoadMatrixVectorSize(char *);
ARB_API cmplx         **LoadMatrix1DC(char *);
ARB_API double       ***LoadMV(char *);
ARB_API long            LoadMVSize(char *);
ARB_API long          **LoadMVSpecs(char *);
ARB_API double        **LoadOneMatrix(FILE *);
ARB_API char         ***LoadStringMatrix(char *);
ARB_API char          **LoadStringVector(char *);
ARB_API char          **LoadStringVectorN(char *,long);
ARB_API double         *LoadVector(char *);
ARB_API double         *LoadVectorBin(char *);
ARB_API cmplx          *LoadVectorC(char *);
ARB_API cmplx          *LoadVectorCBin(char *);
ARB_API long           *LoadVectorL(char *);
ARB_API long           *LoadVectorLBin(char *);
ARB_API double       ***LoadVectorMatrix(char *);
ARB_API double       ***LoadVectorMatrixBin(char *);
ARB_API cmplx        ***LoadVectorMatrixBinC(char *);
ARB_API cmplx        ***LoadVectorMatrixC(char *);
ARB_API void            LoadVectorMatrix1D(char *,double ****,double **);
ARB_API void            LoadVectorMatrix1DC(char *,cmplx ****,double **);
ARB_API double       ***LoadVM(char *);
ARB_API long            LODE1(double *,double *,double *,double,double *);
ARB_API void            LODE1C(cmplx *,cmplx *,double *,cmplx,cmplx *,long);
ARB_API cmplx           LogC(cmplx);
ARB_API cmplx           Log10C(cmplx);
ARB_API char           *LongString(long);
ARB_API char           *LongStringLZ(long,long);
ARB_API void            LongToHSV(long,double *,double *,double *);
ARB_API void            LongToRGB(long,unsigned char *,unsigned char *,unsigned char *);
ARB_API void            LorentzianFit(char *,double **);
ARB_API double          LorentzianFitFunction(double,double *);
ARB_API char           *LowerCase(char *);
ARB_API long           *LUDecomposition(double **);
ARB_API long           *LUDecompositionC(cmplx **);
ARB_API void            LUSolve(double **,double **);
ARB_API void            LUSolveC(cmplx **,cmplx **);
ARB_API double        **Matrix(long,long);
ARB_API cmplx         **MatrixC(long,long);
ARB_API void            MatrixCopy(double **,double **);
ARB_API void            MatrixCopyC(cmplx **,cmplx **);
ARB_API void            MatrixCopyF(float **,float **);
ARB_API void            MatrixCopyI(int **,int **);
ARB_API void            MatrixCopyL(long **,long **);
ARB_API cmplx         **Matrix0C(long,long);
ARB_API uint32_t      **MatrixUInt32(uint32_t,uint32_t);
ARB_API float         **MatrixF(long,long);
ARB_API double        **MatrixFunction(double (*)(double),double **);
ARB_API double        **MatrixFunctionTridag(double (*)(double),double *,double *);
ARB_API int           **MatrixI(long,long);
ARB_API long          **MatrixL(long,long);
ARB_API double        **MatrixV(double *);
ARB_API cmplx         **MatrixVC(cmplx *);
ARB_API double       ***MatrixVector(long);
ARB_API long         ***MatrixVectorL(long);
ARB_API cmplx        ***MatrixVectorC(long);
ARB_API double       ***MatrixVector2D(long,long);
ARB_API cmplx        ***MatrixVector2DC(long,long);
ARB_API double        **MatrixS(double);
ARB_API double        **Matrix0(long,long);
ARB_API long          **Matrix0L(long,long);
ARB_API double          Max(double *);
ARB_API double          MaxAbs(double *);
ARB_API cmplx           MaxC(cmplx *);
ARB_API float           MaxF(float *);
ARB_API int             MaxI(int *);
ARB_API long            MaxL(long *);
ARB_API double          MaxM(double **);
ARB_API double          MaxMC(cmplx **);
ARB_API long            MaxML(long **);
ARB_API long            MaxPos(double *,long,long);
ARB_API long            MaxPosAbs(double *,long,long);
ARB_API long            MaxPosAbsC(cmplx *,long,long);
ARB_API double          MaxWindow(double *,long,long);
ARB_API double          Max2(double,double);
ARB_API long            Max2L(long,long);
ARB_API double          Max3(double,double,double);
ARB_API long            Max3L(long,long,long);
ARB_API double          Mean(double *,double *,double);
ARB_API double          MeanV(double *);
ARB_API double          MeanAbs(double *);
ARB_API double          MeanAbsC(cmplx *);
ARB_API double          MeanM(double **);
ARB_API double          MeanWindow(double *,long,long);
ARB_API double          Median(double *);
ARB_API MEM_FILE       *MemOpen(long);
ARB_API void            MemClose(MEM_FILE *);
ARB_API MEM_FILE       *MemFRead(FILE *f1); 
ARB_API void            MemWrapperClose(MEM_FILE *f);
ARB_API MEM_FILE       *MemWrapperOpen(unsigned char *bytes, unsigned long size);
ARB_API void            MemFWrite(MEM_FILE *m1, FILE *f1);
ARB_API unsigned char  *MemDump(MEM_FILE *f, unsigned long startpos, unsigned long endpos);
ARB_API unsigned char  *MemGetBufHandle(MEM_FILE *f);;
ARB_API unsigned long   MemGetBufPos(MEM_FILE *f);
ARB_API unsigned long   MemGetBufSize(MEM_FILE *f);
ARB_API size_t          MemRead(void *,size_t,size_t,MEM_FILE *); 
ARB_API void            MemSetRelativeBufPos(MEM_FILE *f, long rpos);
ARB_API void            MemSetBufPos(MEM_FILE *f, unsigned long pos);
ARB_API MEM_FILE       *MemWrapperOpen(unsigned char *bytes, unsigned long size);
ARB_API size_t          MemWrite(void *,size_t,size_t,MEM_FILE *);
ARB_API uint32_t        MergeBytesToLong(unsigned char *);
ARB_API double          Metropolis1(double (*)(double *,double *),
                            double (*)(double,double *),long,double *,long,
                            long);
ARB_API double          Min(double *);
ARB_API double          MinAbs(double *);
ARB_API double          MinC(cmplx *);
ARB_API float           MinF(float *);
ARB_API int             MinI(int *);
ARB_API double          MinimizeChiSquare(double **,double (*)(double,double *),
                                  double *,double,double);
ARB_API long            MinL(long *);
ARB_API double          MinM(double **);
ARB_API double          MinMC(cmplx **);
ARB_API long            MinML(long **);
ARB_API long            MinPos(double *);
ARB_API long            MinPosAbs(double *);
ARB_API long            MinPosL(long *);
ARB_API double          MinWindow(double *,long,long);
ARB_API double          Min2(double,double);
ARB_API long            Min2L(long,long);
ARB_API double          Min3(double,double,double);
ARB_API long            Min3L(long,long,long);
ARB_API double          ModularExponentiation(double, double, double);
ARB_API double          Moment(double *,double *,long);
ARB_API double          MomentAbsC(cmplx *,double *,long);
ARB_API double          MomentAbsWindowC(cmplx *,double *,long,long,long,double);
ARB_API long            MomentD(double *,long);
ARB_API long            MomentDC(cmplx *,long);
ARB_API double          MomentWindow(double *,double *,long,long,long,double);
ARB_API double          MomentWindowD(double *,long,long,long,double);
ARB_API double         *MonteCarloIntegral(double (*)(double *), double *, double *,
                                   long, long *);
ARB_API double          MostFrequent(double *);
ARB_API cmplx           MulC(cmplx,cmplx);
ARB_API cmplx           MulCR(cmplx,double);
ARB_API void            MulIPM(double **,double **,double **);
ARB_API void            MulIPMC(cmplx **,cmplx **,cmplx **);
ARB_API double        **MulM(double **,double **);
ARB_API double         *MulMV(double **,double *);
ARB_API double          MulVMV(double *,double **,double *);
ARB_API double         *MulVM(double *,double **);
ARB_API cmplx         **MulMC(cmplx **,cmplx **);
ARB_API cmplx           Mul3C(cmplx,cmplx,cmplx);
ARB_API int             MutuallyExclusive(double,double,double,double);
ARB_API long            NCols(char *);
ARB_API long            NextPow(long,long);
ARB_API double          NewtonRaphson(double (*)(double),double,double,long);
ARB_API double          NewtonRaphsonN(double (*)(double *),double *,long,double,long);
ARB_API double        **Niederreiter(long,long,long *);
ARB_API void            NiederreiterCalcC2(long,long **);
ARB_API void            NiederreiterCalcV2(long,long,long *,long **,long **,long **,
                                   long *,long *,long *);
ARB_API void            NiederreiterPlyMul2(long **,long **,long,long *,long,long *,
                                    long *,long *pc);
ARB_API void            NiederreiterSetfld2(long **add,long **mul,long **);
ARB_API void            Niederreiter1(long,long *,double *);
ARB_API long            NOcc(double,double *);
ARB_API long            NOccL(long,long *);
ARB_API long            NOccS(char,char *);
ARB_API long            NOccWindowS(char,char *,long,long);
ARB_API double          Norm(double *v);
ARB_API cmplx           NormalizeC(cmplx *);
ARB_API double          NormalizeD(double **);
ARB_API double          NormalizeM(double **);
ARB_API double          NormalizeMC(cmplx **);
ARB_API double          NormalizeV(double *);
ARB_API double          NormalizeVC(cmplx *);
ARB_API double          NormalizeWindowV(double *,long,long);
ARB_API cmplx           NormalizeWindowVC(cmplx *,long,long);
ARB_API double          NormalizeVMax(double *);
ARB_API double          NormalizeVMaxC(cmplx *);
ARB_API double          NormC(cmplx *);
ARB_API double          NormWindow(double *,long,long);
ARB_API long            NRows(char *);
ARB_API void            NumerovIntegrator(double *,double *,double *,long,long);
ARB_API int             Odd(double);
ARB_API long            OrderOfMagnitude10(double);
ARB_API long            OrderOfMagnitude1000(double);
ARB_API double        **OuterP(double *,double *);
ARB_API double         *OverlapInterval(double,double,double,double);
ARB_API double          ParabolicExtremum(double,double,double,double,double,double);
ARB_API double          ParabolicInterpolation1(double *,double *,double);
ARB_API double          ParabolicInterpolation11(double *,double *,long,double);
ARB_API double        **ParabolicInterpolation1D(double **,double,double,long);
ARB_API char           *ParentFolder(char *,long);
ARB_API double        **PartialDerivativeMixed2D(double **,double,double);
ARB_API double        **PartialDerivativeX2D(double **, double);
ARB_API double        **PartialDerivativeY2D(double **, double);
ARB_API double        **PartialDerivative2X2D(double **, double);
ARB_API double        **PartialDerivative2Y2D(double **, double);
ARB_API double         *Pascal(long);
ARB_API double          Perm(long,long);
ARB_API double          PhaC(cmplx);
ARB_API double         *Phase(double *,long);
ARB_API long            PoissonRandomNumber(double, long *);
ARB_API long           *PoissonRandomNumbers(double, long, long *);
ARB_API double          PolyEval(double,double *);
ARB_API double       ***PolyFit(double **,long,double,double,long);
ARB_API double          PolyFitFunction(double,double *);
ARB_API void            PolynomialFit(double *,double *,double *,double *,long,long);
ARB_API double        **PolynomialFit2(double *,double *,double *,long);
ARB_API void            PolynomialFitPhase(cmplx **,double *,long,long);
ARB_API double          PolynomialFitFunc(double,double *);
ARB_API long            Pos(double,double *);
ARB_API long            PosC(cmplx,cmplx *);
ARB_API long            PosS(char *,char **);
ARB_API double        **PositiveDomain(double **);
ARB_API cmplx         **PositiveDomainC(cmplx **);
ARB_API long            PosL(long,long *);
ARB_API cmplx           PowC(cmplx,cmplx);
ARB_API cmplx           PowCR(cmplx,double);
ARB_API cmplx           PowIR(double,double);
ARB_API double          PowerMethod(double **,double *,long);
ARB_API cmplx           PowerMethodC(cmplx **,cmplx *,long);
ARB_API double        **PowM(double **,long);
ARB_API cmplx         **PowMC(cmplx **,long);
ARB_API cmplx           PowRC(double,cmplx);
ARB_API cmplx           PowRR(double,double);
ARB_API double          PrimeCount(double,long);
ARB_API double        **PrimeFactorsSmall(double,double *);
ARB_API double        **PrimeFactors2(double, char *);
ARB_API long            PrimeNumber(long);
ARB_API double         *PrimeNumbers(double,char *);
ARB_API double         *PrimeNumbersRange(double, double, char *);
ARB_API void            PrimeNumbersSieve(double,uint64_t,uint64_t ,bit_array *,
                                  double **);
ARB_API double         *PrimeNumbersSmall(double);
ARB_API void            PrimeNumbersSub(double,uint64_t,uint64_t,bit_array *,uint64_t *,
                                double **);
ARB_API double          PrimitiveRootPrime(double, char *);
ARB_API void            Print(double);
ARB_API void            PrintC(cmplx);
ARB_API void            PrintF(float);
ARB_API void            PrintI(int);
ARB_API void            PrintL(long);
ARB_API void            PrintMatrix(double **);
ARB_API void            PrintMatrixC(cmplx **);
ARB_API void            PrintMatrixF(float **);
ARB_API void            PrintMatrixI(int **);
ARB_API void            PrintMatrixL(long **);
ARB_API void            PrintMatrixVector(double ***);
ARB_API void            PrintMatrixVectorC(cmplx ***);
ARB_API void            PrintNL(void);
ARB_API void            PrintS(char *);
ARB_API void            PrintStringMatrixVector(char ****);
ARB_API void            PrintStringVector(char **);
ARB_API void            PrintStringVectorArray(char ***);
ARB_API void            PrintVector(double *);
ARB_API void            PrintVectorC(cmplx *);
ARB_API void            PrintVectorF(float *);
ARB_API void            PrintVectorI(int *);
ARB_API void            PrintVectorL(long *);
ARB_API double          Product(double (*)(double),long,long,long);
ARB_API cmplx           ProductC(cmplx (*)(double),long,long,long);
ARB_API double          ProductP(double (*)(double *),double *,long,long,long);
ARB_API double        **ProximalInterpolationM(double **,long,long);
ARB_API cmplx          *QRAlgorithm(double **);
ARB_API cmplx        ***QRAlgorithm2(double **);
ARB_API void            QRBalance(double **);
ARB_API double         *QRBalance2(double **);
ARB_API void            QRHessenberg(double **);
ARB_API long           *QRHessenberg2(double **,double **);
ARB_API void            QRSortVecs(double **,cmplx *);
ARB_API double          QuarterCircleFunc(double);
ARB_API double          Ran2(long *);
ARB_API double          Ran2P(int32_t *);
ARB_API double          Ran2P2(int32_t *, int32_t, int32_t, int32_t *);
ARB_API double          Ran22(long *,long *,long *,long *);
ARB_API long            RandomInteger(long,long,long *);
ARB_API int32_t         RandomIntegerP(int32_t,int32_t,int32_t *);
ARB_API long            RandomInteger2(long,long,long *,long *,long *,long *);
ARB_API void            RandomizeMatrix(double **,double,double,long *);
ARB_API void            RandomizeMatrixC(cmplx **,double,double,long *);
ARB_API void            RandomizeVector(double *,double,double,long *);
ARB_API void            RandomizeVectorC(cmplx *,double,double,long *);
ARB_API double          RandomReal(double,double,long *);
ARB_API double          RandomRealNormal(double,double,long *);
ARB_API char           *RandomString(long,long *);
ARB_API char           *RandomStringAlphaNum(long,long *);
ARB_API void            RandomStringAlphaNumIP(long ,long *,char *);
ARB_API void            RandomStringAlphaNumIP2(long,long *,long *,long *,long *,
                                        char *);
ARB_API char           *RandomStringAlphaNumP(int32_t,int32_t *);
ARB_API char           *RandomStringBracketed(char *, char *, long *);
ARB_API char           *RandomStringBracketedExtASCII(char *,char *,long *);
ARB_API char           *RandomStringExtASCII(long,long *);
ARB_API double          Range(double *);
ARB_API double          Re(cmplx);
ARB_API void           *ReadMemFile(void *f, stream_type st);
ARB_API double        **ReM(cmplx **);
ARB_API double         *ReV(cmplx *);
ARB_API cmplx           RealC(double);
ARB_API cmplx         **RealMC(double **);
ARB_API cmplx          *RealVC(double *);
ARB_API cmplx        ***RealVMC(double ***);
ARB_API double        **ReferenceDomain(double *,double *);
ARB_API char           *RemoveCharacters(char *, char **);
ARB_API char           *RemoveCharacter(char,char **);
ARB_API void            RemoveCol(long,double ***);
ARB_API void            RemoveColS(long,char  ****);
ARB_API double          RemoveDC(double *);
ARB_API cmplx           RemoveDCC(cmplx *);
ARB_API double          RemoveDCWindow(double *,long,long);
ARB_API void            RemoveDuplicateCols(double ***,char *);
ARB_API void            RemoveDuplicateColsS(char ****,char *);
ARB_API void            RemoveDuplicateRows(double ***,char *);
ARB_API void            RemoveDuplicates(double **,char *);
ARB_API void            RemoveDuplicatesSorted(double **);
ARB_API void            RemoveDuplicateStrings(char ***,char *);
ARB_API void            RemoveElement(double,double **);
ARB_API void            RemoveElements(double *,double **);
ARB_API void            RemoveElementsL(long *,long **);
ARB_API double          RemoveLinearPhase(cmplx *);
ARB_API double          RemovePosition(long,double **);
ARB_API char           *RemovePositionCh(long,char **);
ARB_API long            RemovePositionL(long,long **);
ARB_API void            RemovePositions(long *,double **);
ARB_API char           *RemovePositionS(long,char ***);
ARB_API char           *RemovePositionsC(long *, char **);
ARB_API char         ***RemovePositionSM(long,char *****);
ARB_API void            RemovePositionsS(long *,char ***);
ARB_API long            RemovePositionUint32(long, uint32_t **);
ARB_API void            RemoveElementL(long,long **);
ARB_API double         *RemoveRow(long,double ***);
ARB_API void            RemoveRows(long *,double ***);
ARB_API char          **RemoveRowS(long,char ****);
ARB_API void            ReorderMV(double ****,long *);
ARB_API void            ReorderV(double **,long *);
ARB_API void            ReverseS(char *);
ARB_API void            ReverseV(double *);
ARB_API void            ReverseVC(cmplx *);
ARB_API uint32_t        RGBToLong(unsigned char,unsigned char,unsigned char);
ARB_API void            RGBToHSV(double,double,double,double *,double *,double *);
ARB_API uint32_t        RightRotate32(uint32_t,int32_t);
ARB_API uint64_t        RightRotate64(uint64_t,int64_t);
ARB_API double        **RKN(void (*)(double,double *,double *,double *),double,
                    double *,double *,double,long);
ARB_API double        **RK1(double (*)(double,double),double,double,double,long);
ARB_API double        **RK2(double (*)(double,double,double),double,double,double,
                    double,long);
ARB_API double          Round(double);
ARB_API char           *RoundDoubleString(char **,long);
ARB_API double         *Row(double **,long);
ARB_API void            RungeKutta(double (*)(double,double,double),double,double,
                           double,double,double *,double *,long);
ARB_API void            RungeKuttaN(void (*)(double,double *,double *,double *),double,
                            double *,double *,double,double *,double *,long);
ARB_API double          RungeKuttaIntegral1D(void (*)(double,double *,double *),double,
                                     double);
ARB_API double        **RungeKuttaIntegrator(double *,double,double,long,double,double,
                                     double,double,long *,long *,
                                     void (*)(double,double *,double *));
ARB_API void            RungeKuttaCashKarp(double *,double *,double,double,double *,
                                   double *,void (*)(double,double *,double *));
ARB_API void            RungeKuttaStep(double *,double *,double *,double,double,double,
                               double *,double *,double *,
                               void (*)(double,double *,double *));
ARB_API cmplx           RungeKuttaIntegral1DC(void (*)(double,cmplx *,cmplx *),double,
                                      double);
ARB_API cmplx         **RungeKuttaIntegratorC(cmplx *,double,double,long,double,double,
                                      double,double,long *,long *,
                                      void (*)(double,cmplx *,cmplx *));
ARB_API void            RungeKuttaCashKarpC(cmplx *,cmplx *,double,double,cmplx *,
                                    cmplx *,void (*)(double,cmplx *,cmplx *));
ARB_API void            RungeKuttaStepC(cmplx *,cmplx *,double *,double,double,double,
                                double *,double *,double *,
                                void (*)(double,cmplx *,cmplx *));
ARB_API double        **RungeKuttaSampler(double *,double,double,double,double,double,
                                  double,long *,long *,
                                  void (*)(double,double *,double *));
ARB_API cmplx         **RungeKuttaSamplerC(cmplx *,double,double,double,double,double,
                                   double,long *,long *,
                                   void (*)(double,cmplx *,cmplx *));
ARB_API char           *SaltedHash(char *,long *);
ARB_API char           *SaltFromHash(char *);
ARB_API double        **SavitzkyGolay(double *,double *,long,long,long,long);
ARB_API double         *SavitzkyGolayCoefficients(long,long,long,long);
ARB_API double       ***ScalarField2D(double **,double *,double *);
ARB_API long            SeedFrom4Bytes(MEM_FILE *,unsigned char *,long);
ARB_API long            SeedFrom4String(char *,long);
ARB_API long            SeedFromByteArray(MEM_FILE *,long,long,long);
ARB_API long            SeedFromString(char *);
ARB_API int32_t         SeedFromStringP(char *);
ARB_API char           *SHA256(char *);
ARB_API char           *SHA512(char *);
ARB_API void            Shift(double *,long);
ARB_API void            ShiftC(cmplx *,long);
ARB_API void            ShiftCol(double **,long);
ARB_API void            ShiftColC(cmplx **,long);
ARB_API void            ShiftCont(double **,double);
ARB_API void            ShiftCont1DC(cmplx **,double);
ARB_API void            ShiftRow(double **,long);
ARB_API void            ShiftRowC(cmplx **,long);
ARB_API void            Shift0C(cmplx *,long);
ARB_API double          Sign(double);
ARB_API double          SimpsonIntegralAdaptive(double (*)(double),double,double,
                                        double);
ARB_API double          SimpsonIntegralAdaptiveSub(double (*)(double),double,double,
                                           double,double,double,double,double);
ARB_API double          SimpsonIntegalSphericalCS3D(double ***,double *,double *,double,
                                            double,double,double);
ARB_API double          SimpsonIntegalSpherical3D(double ***,double *,double *,
                                          double *,double,double,double,double,
                                          double,double);
ARB_API double          SimpsonIntegral1(double, double, double, double, double);
ARB_API double          SimpsonIntegral1D(double *,double *,double,double);
ARB_API cmplx           SimpsonIntegral1DC(cmplx *,double *,double,double);
ARB_API void            SimulatedAnnealing(double (*)(double *),double *,double **,
                                   double *,long *);
ARB_API cmplx           SinC(cmplx);
ARB_API double          SincFunc(double);
ARB_API cmplx           SinhC(cmplx);
ARB_API char           *SmileyHash(char **, long);
ARB_API char           *SmileyHexHash(char **,long);
ARB_API void            SmoothFFT(cmplx **,long,long);
ARB_API void            SmoothFFTV(cmplx *,long,long);
ARB_API void            SmoothSpline(double *,double *,long);
ARB_API void            SmoothSplinePolar1DC(cmplx *,double *,long);
ARB_API void            SmoothSplineRectangular1DC(cmplx *,double *,long);
ARB_API double          SmoothSquareFunc(double,double,double,double);
ARB_API cmplx           SmoothSquareFuncC(double,double,double,double);
ARB_API double          SmoothSquareFunc2D(double,double,double,double,double,double,
                                   double,double);
ARB_API double          SmoothStepFunc(double);
ARB_API void            Sort(double *);
ARB_API void            SortL(long *);
ARB_API void            SortMatrix(double **,long,long,long);
ARB_API void            SortMatrixL(long **,long,long,long);
ARB_API void            SortMatrixCols(double **,double *);
ARB_API void            SortMatrixColsC(cmplx **,double *);
ARB_API void            SortMatrixRows(double **,double *);
ARB_API void            SortMatrixUInt32(uint32_t **, long, long, long);
ARB_API void            SortMatrixVector(double ***,double *);
ARB_API void            SortMatrix1DC(cmplx **);
ARB_API void            SortStringMatrix(char ****,long,long,long);
ARB_API void            SortStringVector(char ***);
ARB_API void            SortStringVectorRef(char ***,double *);
ARB_API char          **SortStringVectorSizes(char ***);
ARB_API double        **Spectrum(cmplx **);
ARB_API cmplx           SphericalHarmonic(double,double,long,long);
ARB_API double        **SplineInterpolation(double **,double,double,long);
ARB_API double        **SplineInterpolationCols(double **,long);
ARB_API void            SplineInterpolationPts(double *,double *,double *,double *);
ARB_API double         *SplineInterpolationV(double *,long);
ARB_API cmplx         **SplineInterpolation1DC(cmplx **,double,double,long);
ARB_API double       ***SplineInterpolation2D(double ***,long,long);
ARB_API void            SplitLongToBytes(unsigned long,unsigned char *);
ARB_API void            SplitShortToBytes(short,unsigned char *);
ARB_API double          Sqr(double);
ARB_API cmplx           SqrC(cmplx);
ARB_API double          Sqrt(double);
ARB_API cmplx           SqrtC(cmplx);
ARB_API double          SquareFunc(double);
ARB_API double          StandardDeviation(double *,double *,double);
ARB_API double          StandardDeviationV(double *);
ARB_API cmplx        ***STFT(cmplx **,double **,double,double,long);
ARB_API void           *StrDup(void *);
ARB_API char           *StringChar(char);
ARB_API long            StringCols(char *);
ARB_API char         ***StringCombinations(char **, long, char *);
ARB_API char           *StringEncoding(long,char *);
ARB_API char           *StringInit(long);
ARB_API char         ***StringMatrix(long,long);
ARB_API long            StringMatrixCols(char ***);
ARB_API long            StringMatrixRows(char ***);
ARB_API char        ****StringMatrixVector(long);
ARB_API long            StringMatrixVectorLength(char ****);
ARB_API char         ***StringPermutations(char **, long, char *);
ARB_API long            StringRows(char *);
ARB_API char           *StringSV(char **);
ARB_API unsigned char  *StringToByteArray(char *,long *);
ARB_API char          **StringTokens(char *, char *);
ARB_API char          **StringTokensS(char *, char *);
ARB_API long           *StringToVectorL(char *);
ARB_API char          **StringVector(long);
ARB_API char          **StringVectorNULL(long);
ARB_API char          **StringVector2(char *, char *);
ARB_API char          **StringVector3(char *,char *,char *);
ARB_API long            StringVectorLength(char **);
ARB_API char           *StrStr(char *, char *);
ARB_API char           *StrTokR(char *, const char *, char **);
ARB_API double          StudentCriticalValue(double,double);
ARB_API double          StudentPValue(double,double);
ARB_API cmplx           SubC(cmplx,cmplx);
ARB_API double        **SubM(double **,double **);
ARB_API cmplx         **SubMC(cmplx **,cmplx **);
ARB_API char           *SubstituteString(char **,char *,char *);
ARB_API char           *SubstituteStringEqualSizes(char **,char *,char *);
ARB_API double         *SubV(double *,double *);
ARB_API cmplx          *SubVC(cmplx *,cmplx *);
ARB_API double          Summation(double (*)(double),long,long,long);
ARB_API double          SummationP(double (*)(double *),double *,long,long,long);
ARB_API cmplx           SummationC(cmplx (*)(double),long,long,long);
ARB_API double          SuperGaussianFilter(double,double,double,double);
ARB_API double       ***SVD(double **);
ARB_API double          SVDPythag(double,double);
ARB_API void            Swap(double *,double *);
ARB_API void            SwapC(cmplx *,cmplx *);
ARB_API void            SwapCols(double **,long,long);
ARB_API void            SwapColsC(cmplx **,long,long);
ARB_API void            SwapL(long *,long *);
ARB_API void            SwapRows(double **,long,long);
ARB_API void            SwapRowsC(cmplx **,long,long);
ARB_API void            SwapRowsL(long **,long,long);
ARB_API double        **SymmetricDifferenceM(double **,double **);
ARB_API double         *SymmetricDifferenceV(double *,double *);
ARB_API double         *SymmetrizeTridiagonalMatrix(double *,double *,double *);
ARB_API cmplx           TanC(cmplx);
ARB_API cmplx           TanhC(cmplx);
ARB_API long            TestAgainstSaltedHash(char *,char *);
ARB_API long            TestElementS(char *,char **);
ARB_API long            TestPrimeFactors(double, double **,long);
ARB_API double          TimeBandwidthProduct(cmplx **);
ARB_API double          TimeBandwidthProductSpectrum(double **,double);
ARB_API double          ToothFunc(double);
ARB_API char           *TrailingCharacters(char **,char,char,long);
ARB_API double        **Transpose(double **);
ARB_API cmplx         **TransposeC(cmplx **);
ARB_API void            TransposeIPC(cmplx **,cmplx **);
ARB_API void            TransposeIP(double **,double **);
ARB_API void            TransposeSquareIPC(cmplx **);
ARB_API void            TransposeSquareIP(double **);
ARB_API double          TrapezoidIntegral(double *,double *,long,long);
ARB_API cmplx           TrapezoidIntegralC(cmplx *,cmplx *,long,long);
ARB_API double          TrapezoidIntegral1D(double *,double *,double,double);
ARB_API double          TriangleFunc(double);
ARB_API void            TriagonalEigenDecomposition(double *,double *,double **);
ARB_API void            TriagonalEigenValues(double *,double *);
ARB_API double        **TriDiagonalMatrix(double *,double *,double *);
ARB_API cmplx         **TriDiagonalMatrixC(cmplx *,cmplx *,cmplx *);
ARB_API void            TriDiagSolve(double *,double *,double *,double *,double *);
ARB_API void            TriDiagSolveC(cmplx *,cmplx *,cmplx *,cmplx *,cmplx *);
ARB_API uint64_t        UInt64FromString(char *);
ARB_API char           *UInt64String(uint64_t);
ARB_API uint32_t        UnicodeBytesToInteger(unsigned char *);
ARB_API void            UniformScale(double **);
ARB_API double        **UnionM(double **,double **,char *);
ARB_API double         *UnionV(double *,double *,char *);
ARB_API double         *UnitV(double *);
ARB_API cmplx          *UnitVC(cmplx *);
ARB_API double          UnityFunc(double);
ARB_API cmplx           UnityFunc1DC(double);
ARB_API void            UnwrapM(double **,double);
ARB_API void            UnwrapV(double *,double);
ARB_API void            UnwrapVWindow(double *,double,long,long);
ARB_API char           *UpperCase(char *);
ARB_API double         *Vector(long);
ARB_API cmplx          *VectorC(long);
ARB_API void            VectorCopy(double *,double *);
ARB_API void            VectorCopyC(cmplx *,cmplx *);
ARB_API void            VectorCopyF(float *,float *);
ARB_API void            VectorCopyI(int *,int *);
ARB_API void            VectorCopyL(long *,long *);
ARB_API float          *VectorF(long);
ARB_API int            *VectorI(long);
ARB_API long           *VectorL(long);
ARB_API double       ***VectorMatrix(long,long,long);
ARB_API cmplx        ***VectorMatrixC(long,long,long);
ARB_API double       ***VectorMatrix0(long,long,long);
ARB_API double         *VectorS(double);
ARB_API void            VectorToByteArrayL(long *,unsigned char *,long);
ARB_API void            VectorToByteArrayUint32(uint32_t *, unsigned char *, long);
ARB_API char           *VectorToStringL(long *);
ARB_API uint32_t       *VectorUint32(long );
ARB_API double         *Vector0(long);
ARB_API cmplx          *Vector0C(long);
ARB_API long           *Vector0L(long);
ARB_API long           *Vector1L(long);
ARB_API double         *Vector2(double,double);
ARB_API double         *Vector3(double,double,double);
ARB_API long           *Vector3L(long,long,long);
ARB_API double         *Vector6(double,double,double,double,double,double);
ARB_API long           *Vector7L(long,long,long,long,long,long,long);
ARB_API double         *VonNeumannRejection(double (*)(double), double, double, double,
                                    long, long *);
ARB_API double 		      WeightedMean(double *,double *);
ARB_API double 		      WeightedMeanWindow(double *,double *,long,long);
ARB_API double       ***WignerDistribution(cmplx **,long);
ARB_API double          Wigner3J(double,double,double,double,double,double);
ARB_API double        **Wigner3JSymbols(double,double,double,double);
ARB_API double          Wigner3JSymbolSpecial1(double,double,double);
ARB_API double          Wigner3JSymbolSpecial2(double,double,double);
ARB_API void            Wigner3JSymbolsSpecial(double,double,double *,double *);
ARB_API double          Wrap(double,double);
ARB_API void            Write(char *,double);
ARB_API void            WriteByteArray(char *,unsigned char *,long);
ARB_API void            WriteMemFile(void *val, void *f, stream_type st);
ARB_API void            WriteMV(char *,double ***);
ARB_API void            WriteMatrix(char *, double **);
ARB_API void            WriteMatrixUInt32(char *,uint32_t **);
ARB_API void            WriteMatrixAbsC(char *,cmplx **);
ARB_API void            WriteMatrixBin(char *,double **);
ARB_API void            WriteMatrixC(char *,cmplx **);
ARB_API void            WriteMatrixCBin(char *,cmplx **);
ARB_API void            WriteMatrixF(char *,float **);
ARB_API void            WriteMatrixI(char *,int **);
ARB_API void            WriteMatrixL(char *,long **);
ARB_API void            WriteMatrixLBin(char *,long **);
ARB_API void            WriteMatrixVector(char *,double ***);
ARB_API void            WriteMatrixVectorBin(char *,double ***);
ARB_API void            WriteMatrixVectorC(char *,cmplx ***);
ARB_API void            WriteMatrixVectorFiles(char *,double ***);
ARB_API void            WriteMatrixVectorL(char *,long ***);
ARB_API void            WriteMatrix1DC(char *,cmplx **);
ARB_API void            WriteMatrix1DCErr(char *,cmplx **,cmplx *);
ARB_API void            WriteOneMatrix(FILE *,double **);
ARB_API void            WriteS(FILE *,char *);
ARB_API void            WriteString(char *,char *);
ARB_API void            WriteStringMatrix(char *,char ***);
ARB_API void            WriteStringVector(char *,char **);
ARB_API void            WriteVector(char *,double *);
ARB_API void            WriteVectorAbsC(char *,cmplx *);
ARB_API void            WriteVectorBin(char *,double *);
ARB_API void            WriteVectorC(char *, cmplx *);
ARB_API void            WriteVectorCBin(char *,cmplx *);
ARB_API void            WriteVectorF(char *,float *);
ARB_API void            WriteVectorI(char *,int *);
ARB_API void            WriteVectorL(char *,long *);
ARB_API void            WriteVectorLBin(char *,long *);
ARB_API void            WriteVectorMatrix(char *,double ***);
ARB_API void            WriteVectorMatrixC(char *,cmplx ***);
ARB_API void            WriteVectorMatrixBin(char *,double ***);
ARB_API void            WriteVectorMatrixBinC(char *,cmplx ***);
ARB_API void            WriteVectorMatrix1D(char *,double ***,double *);
ARB_API void            WriteVectorMatrix1DC(char *,cmplx ***,double *);
ARB_API void            WriteVM(char *,double ***);
ARB_API double          ZeroEndPoints(double **);
ARB_API double          ZeroFunc(double);
ARB_API cmplx           ZeroFunc1DC(double);
ARB_API cmplx         **ZeroPadFFT(cmplx **,long);
ARB_API cmplx          *ZeroPadVC(cmplx *,long);
ARB_API double        **ZeroPad1D(double **,long);
ARB_API cmplx         **ZeroPad1DC(cmplx **,long);

ARB_API bit_array *InitializeBitArray(uintmax_t size);
ARB_API void FreeBitArray(void *bav);
ARB_API void SetBit(bit_array *ba, uintmax_t pos);
ARB_API void ClearBit(bit_array *ba, uintmax_t pos);
ARB_API uint64_t GetBit(bit_array *ba, uintmax_t pos);
ARB_API char *PrintBitArray(bit_array *ba);
ARB_API void SetAllBits(bit_array *ba);
ARB_API void SetBitsFromInt64(bit_array *ba, uintmax_t pos, uint64_t num);
ARB_API uint32_t GetInt32(bit_array *ba, uintmax_t pos);
ARB_API uint64_t GetInt64(bit_array *ba, uintmax_t pos);
ARB_API void SetBitsFromChar(bit_array *ba, uintmax_t pos, char c);
ARB_API void ClearAllBits(bit_array *ba);

extern long arb_counter;
extern int prime_numbers[2048];
extern size_t max_size_t;
extern uint64_t one_bit_array64_h;
extern uint32_t one_bit_array32_h;


#endif // ARB_H
