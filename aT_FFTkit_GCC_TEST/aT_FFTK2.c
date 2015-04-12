/*////////////////////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////////////////////*/
/*  developed by DRUID's labs.                                                */
/*  vitaliy.druid3@yandex.ru    druidthree@gmail.com                          */
/*  http://dsp.la.net.ua                                                      */
/*                                                                            */
/*                 00  000   0  0  0  000   0   000                           */
/*                000  0  0  0  0  0  0  0     00        *******              */
/*               0100  000   0  0  0  0  0       00     * ****  *             */
/*              01100  0  0   000  0  000      000     *  *   *  *            */
/*             011 00                                  *  ****   *            */
/*            011  00     11   1      111               * *   * *             */
/*           00000000    1  1  111   11                  *******              */
/*          00000000011  1111  1  1    11                                     */
/*          11111111111  1  1  111   111                                      */
/*                                                                            */
/*////////////////////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////////////////////*/
/* aT_ditRadix2Radix4FFTkit.c                                                 */
/* distributed under GNU BSD                                                  */
/*////////////////////////////////////////////////////////////////////////////*/

#include "aT_FFTK2.h"

#define NCMUL(Iout, Qout, I0, Q0, I1, Q1)                                      \
        Iout = ((I0)*(I1) + (Q0)*(Q1));                                        \
        Qout = ((I1)*(Q0) - (I0)*(Q1));

#define CMUL(Iout, Qout, I0, Q0, I1, Q1)                                       \
        Iout = ((I0)*(I1) - (Q0)*(Q1));                                        \
        Qout = ((I0)*(Q1) + (I1)*(Q0));

/*============================================================================*/
static inline void ifn_NCMUL
(
    float *Iout,
    float *Qout,

    float  I0,
    float  Q0,
    float  I1,
    float  Q1
)
/*----------------------------------------------------------------------------*/
{
    Iout[0] = I0*I1 + Q0*Q1;
    Qout[0] = I1*Q0 - I0*Q1;

    return;
}
/*============================================================================*/
static inline void ifn_CMUL
(
    float *Iout,
    float *Qout,

    float  I0,
    float  Q0,
    float  I1,
    float  Q1
)
/*----------------------------------------------------------------------------*/
{
    Iout[0] = I0*I1 - Q0*Q1;
    Qout[0] = I0*Q1 + I1*Q0;

    return;
}
/*============================================================================*/
void fn_bRev2
(
    float                    a_Re[],
    float                    a_Im[],

    const unsigned int const cuic_N
)
/*----------------------------------------------------------------------------*/
{
    register float buff;

    unsigned int
    ui_temp,
    ui_bNum,
    ui_bRev,
    ui_k;

    for (ui_bNum=1; ui_bNum<cuic_N; ui_bNum++)
        {
            ui_bRev = 0;
            ui_temp = ui_bNum;                                                  /* buff */

            for (ui_k=1; ui_k<cuic_N; ui_k<<=1)                                 /* cascades */
                {
                    ui_bRev <<= 1;                                              /* "1" with in "bit stack", first time - idle */
                    ui_bRev  += (ui_temp&1);                                    /* "1" in bit stack 1 if "1" in bit stack 2 */
                    ui_temp >>= 1;                                              /* low digit is gone forever, last time - idle */
                }
            if (ui_bRev>ui_bNum)
                {
                    buff          = a_Re[ui_bRev];
                    a_Re[ui_bRev] = a_Re[ui_bNum];
                    a_Re[ui_bNum] = buff;

                    buff          = a_Im[ui_bRev];
                    a_Im[ui_bRev] = a_Im[ui_bNum];
                    a_Im[ui_bNum] = buff;
                }
        }
    return;
}
/**
--------------------------------------------------------------------------------
int fn_aT_TR2FFTkern
(
  const unsigned int const  ,                             // - number of samples
  float                    *,                                  // - real samples
  float                    *,                             // - imaginary samples
  const float const        *,                           // - trigonometry tables
  const unsigned int const                      // - step on trigonometry tables
);
--------------------------------------------------------------------------------
*/
/*============================================================================*/
void fn_aT_TR2FFTkern
(
    float                    a_Re[],
    float                    a_Im[],

    const unsigned int const cuic_N,

    const float const        acc_tt[],
    const unsigned int const cuic_jmp
)
/*----------------------------------------------------------------------------*/
{
    unsigned int
    ui_S,                                             /* sinus' table counter */
    ui_C,                                           /* cosinus' table counter */
    ui_kk,                                                /* cascaded counter */
    ui_wc,                                          /* constant roots counter */
    ui_wd;                                        /* differents roots counter */

    float
    tmpRe,                                                /* buffer variables */
    tmpIm;                                                /* buffer variables */

    for (ui_kk=1; ui_kk<cuic_N; ui_kk<<=1)                                      /* cascades */
        {
            for (ui_wc=0; ui_wc<ui_kk; ui_wc++)                                 /* constant roots */
                {
                    ui_S = ((cuic_jmp*ui_wc*cuic_N)>>1)/ui_kk ;
                    ui_C =   cuic_jmp*(cuic_N>>2) + ui_S;                       /* cos(x) = sin(x+90deg) */

                    for (ui_wd=ui_wc; ui_wd<cuic_N; ui_wd+=(ui_kk<<1))          /* differents roots */
                        {
                            ifn_NCMUL
                            (
                                &tmpRe,
                                &tmpIm,
                                a_Re[ui_wd+ui_kk],
                                a_Im[ui_wd+ui_kk],
                                acc_tt[ui_C],                                   /*     cos((3.14*ui_wc)/ui_kk) */
                                acc_tt[ui_S]                                   /* -1.*sin((3.14*ui_wc)/ui_kk) */
                            );

                            a_Re[ui_wd+ui_kk] = a_Re[ui_wd] - tmpRe;
                            a_Im[ui_wd+ui_kk] = a_Im[ui_wd] - tmpIm;

                            a_Re[ui_wd]       = a_Re[ui_wd] + tmpRe;
                            a_Im[ui_wd]       = a_Im[ui_wd] + tmpIm;
                        }
                }
        }
    return;
}

