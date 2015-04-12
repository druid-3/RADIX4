#include "trigArray.h"
#include "aT_FFTK2.h"
#include "aT_FFTK4.h"
#include "FFTwrap.h"

#define STEP(SMPL, TAB)  ((TAB)/(SMPL))

#define PN4 16
#define PN2 8
#define N_FFTL 1024

int main(int argc, char *argv[])
{
#if 1
    /* ***************************** 2 point kernel ***************************** */
    float a_ReTr[N_FFTL], a_ImTr[N_FFTL];

    fn_Test_DDS_C(a_ReTr, N_FFTL, 128., true);
    fn_Test_DDS_C(a_ImTr, N_FFTL, 63.,  true);
    fn_Test_DDS_S(a_ReTr, N_FFTL, 128., false);
    fn_Test_DDS_S(a_ImTr, N_FFTL, 48.,  false);

    printf("input\n");
    fn_aPrint(a_ReTr, N_FFTL);
    fn_aPrint(a_ImTr, N_FFTL);

    printf(" \\\// radix4 start \\\// \n");
    printf("\n");

    ifn_dirFFT4
    (
        a_ReTr,
        a_ImTr,
        N_FFTL,
        &acc_TT4[N_FFTL>>2],
        acc_TT4,
        STEP(N_FFTL, N_FFTL)
    );

    printf("output\n");
    fn_aPrint(a_ReTr, N_FFTL);
    fn_aPrint(a_ImTr, N_FFTL);

    ifn_invFFT4
    (
        a_ReTr,
        a_ImTr,
        N_FFTL,
        &acc_TT4[N_FFTL>>2],
        acc_TT4,
        STEP(N_FFTL, N_FFTL)
    );

    printf("input again\n");
    fn_aPrint(a_ReTr, N_FFTL);
    fn_aPrint(a_ImTr, N_FFTL);

    /*----------------------------------------------------------------------------*/
    printf(" \\\// radix2 start \\\// \n");
    printf("\n");

    ifn_dirFFT2
    (
        a_ReTr,
        a_ImTr,
        N_FFTL,
        acc_TT2,
        STEP(N_FFTL, N_FFTL)
    );

    printf("output\n");
    fn_aPrint(a_ReTr, N_FFTL);
    fn_aPrint(a_ImTr, N_FFTL);

    ifn_invFFT2
    (
        a_ReTr,
        a_ImTr,
        N_FFTL,
        acc_TT2,
        STEP(N_FFTL, N_FFTL)
    );

    printf("input again\n");
    fn_aPrint(a_ReTr, N_FFTL);
    fn_aPrint(a_ImTr, N_FFTL);
    /* ***************************** 4 point kernel ***************************** */
#else
    /* ***************************** 2 point kernel ***************************** */
    /* ***************************** 4 point kernel ***************************** */

    /* 4 point FFT

    float i[4] = {0.5,  0.,-0.5, 0.};
    float q[4] = {  0., 0., 0.,  0.};

    printf("input\n");
    fn_aPrint(i, 4);
    fn_aPrint(q, 4);

    fn_a_FFT4(i, q);

    printf("output\n");
    fn_aPrint(i, 4);
    fn_aPrint(q, 4);

    fn_a_invFFTdrSciling(i, q, 4);
    fn_a_FFT4(q, i);

    printf("input again\n");
    fn_aPrint(i, 4);
    fn_aPrint(q, 4);

    */
    float i2[PN2] = {4., 2., 1., 4., 6., 3., 5., 2.};
    float q2[PN2] = {0., 0., 0., 0., 0., 0., 0., 0.};

    float i4[PN4] = {1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13., 14., 15., 16.};
    float q4[PN4] = {0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.};

    /*----------------------------------------------------------------------------*/
    printf("input\n");
    fn_aPrint(i4, PN4);
    fn_aPrint(q4, PN4);

    printf(" \\\// radix4 start \\\// \n");
    printf("\n");

    DIRFFT4
    (
        i4,
        q4,
        PN4,
        &acc_TT4[N_FFTL>>2],
        acc_TT4,
        STEP(PN4, N_FFTL)
    );

    printf("output\n");
    fn_aPrint(i4, PN4);
    fn_aPrint(q4, PN4);

    INVFFT4
    (
        i4,
        q4,
        PN4,
        &acc_TT4[N_FFTL>>2],
        acc_TT4,
        STEP(PN4, N_FFTL)
    );

    printf("input again\n");
    fn_aPrint(i4, PN4);
    fn_aPrint(q4, PN4);

    /*----------------------------------------------------------------------------*/
    printf(" \\\// radix2 start \\\// \n");

    printf("input\n");
    fn_aPrint(i2, PN2);
    fn_aPrint(q2, PN2);

    printf("\n");

    DIRFFT2
    (
        i2,
        q2,
        PN2,
        acc_TT2,
        STEP(PN2, N_FFTL)
    );

    printf("output\n");
    fn_aPrint(i2, PN2);
    fn_aPrint(q2, PN2);

    INVFFT2
    (
        i2,
        q2,
        PN2,
        acc_TT2,
        STEP(PN2, N_FFTL)
    );

    printf("input again\n");
    fn_aPrint(i2, PN2);
    fn_aPrint(q2, PN2);
#endif

    return 0;
}
/*------------------------------------------------------------------------------
 I[time]{4., 2., 1., 4., 6., 3., 5., 2.} --->
--------------------------------------------------------------------------------
 I[freq]
  {
   27.,
   -4.121319,
    4.,
    0.121321,
    5.,
    0.121330,
    4.,
   -4.121333
  }
--------------------------------------------------------------------------------
 Q[time]{0., 0., 0., 0., 0., 0., 0., 0.} --->
--------------------------------------------------------------------------------
 Q[freq]
  {
    0.,
    3.292894,
    1.,
   -4.707105,
    0.,
    4.707109,
   -1.0,
   -3.292882

  }
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 4 point test

 I[time]{0.5, 0.,-0.5, 0.}  --->  I[freq]{ 0., 1., 0., 1.}
 Q[time]{ 0., 0.,  0., 0.}  --->  Q[freq]{ 0., 0., 0., 0.}

------------------------------------------------------------------------------*/
