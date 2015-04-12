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

#include "aT_FFTK4.h"

#define NCMUL(I0, Q0, I1, Q1, Iout, Qout)                                      \
        Iout = ((I0)*(I1) + (Q0)*(Q1));                                        \
        Qout = ((I1)*(Q0) - (I0)*(Q1));

#define CMUL(I0, Q0, I1, Q1, Iout, Qout)                                       \
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
void fn_bRev4
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
    ui_qRev,
    ui_k;

    for (ui_bNum=1; ui_bNum<cuic_N; ui_bNum++)
        {
            ui_qRev = 0;
            ui_temp = ui_bNum;                                                  /* buff */

            for (ui_k=1; ui_k<cuic_N; ui_k<<=2)                                 /* cascades */
                {
                    ui_qRev <<= 2;                                              /* "11" with in "bit stack", first time - idle */
                    ui_qRev  += (ui_temp&3);                                    /* "11" in bit stack 1 and put result bit stack 2 */
                    ui_temp >>= 2;                                              /* two low digits is gone forever, last time - idle */
                }
            if (ui_qRev>ui_bNum)
                {
                    buff          = a_Re[ui_qRev];
                    a_Re[ui_qRev] = a_Re[ui_bNum];
                    a_Re[ui_bNum] = buff;

                    buff          = a_Im[ui_qRev];
                    a_Im[ui_qRev] = a_Im[ui_bNum];
                    a_Im[ui_bNum] = buff;
                }
        }
    return;
}
/**
--------------------------------------------------------------------------------
int fn_aT_TR4FFTkern
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
void fn_aT_TR4FFTkern
(
    float                    a_Re[],
    float                    a_Im[],

    const unsigned int const cuic_N,

    const float const        acc_ct[],
    const float const        acc_st[],
    const unsigned int const cuic_jmp
)
/*----------------------------------------------------------------------------*/
{
    unsigned int
    ui_i,                                      /* trigonometry tables counter */
    ui_kk,                                                /* cascaded counter */
    ui_wc,                                          /* constant roots counter */
    ui_wd;                                        /* differents roots counter */

    float
    tmpRe0,                                               /* buffer variables */
    tmpIm0,                                               /* buffer variables */
    tmpRe1,                                               /* buffer variables */
    tmpIm1,                                               /* buffer variables */
    tmpRe2,                                               /* buffer variables */
    tmpIm2,                                               /* buffer variables */
    tmpRe3,                                               /* buffer variables */
    tmpIm3;                                               /* buffer variables */

    float
    re_buff0,                                             /* buffer variables */
    im_buff0,                                             /* buffer variables */
    re_buff1,                                             /* buffer variables */
    im_buff1,                                             /* buffer variables */
    re_buff2,                                             /* buffer variables */
    im_buff2,                                             /* buffer variables */
    re_buff3,                                             /* buffer variables */
    im_buff3;                                             /* buffer variables */

    unsigned int
    ui_stp0,                                                /* buffer counter */
    ui_stp1,                                                /* buffer counter */
    ui_stp2,                                                /* buffer counter */
    ui_stp3;                                                /* buffer counter */

    for (ui_kk=1; ui_kk<cuic_N; ui_kk<<=2)                                      /* cascades */
        {
            for (ui_wc=0; ui_wc<ui_kk; ui_wc++)                                 /* constant roots */
                {
                    ui_i  = (cuic_jmp*cuic_N*ui_wc)/(ui_kk<<2);                 /* first step - max (ui_kk<<2 because ui_kk<<=2) */
                    for (ui_wd=ui_wc; ui_wd<cuic_N; ui_wd+=(ui_kk<<2))          /* differents roots */
                        {
                            ui_stp0 = ui_wd;
                            ui_stp1 = ui_wd + ui_kk;
                            ui_stp2 = ui_wd + (ui_kk<<1);
                            ui_stp3 = ui_wd + (ui_kk<<1) + ui_kk;

                            tmpRe0 = a_Re[ui_stp0];
                            tmpIm0 = a_Im[ui_stp0];

                            /* NCMUL(a_Re[ui_stp1],a_Im[ui_stp1],cos(   (3.14*ui_wc)/(ui_kk<<1)),sin(   (3.14*ui_wc)/(ui_kk<<1)),tmpRe1,tmpIm1); */
                            /* NCMUL(a_Re[ui_stp2],a_Im[ui_stp2],cos(2.*(3.14*ui_wc)/(ui_kk<<1)),sin(2.*(3.14*ui_wc)/(ui_kk<<1)),tmpRe2,tmpIm2); */
                            /* NCMUL(a_Re[ui_stp3],a_Im[ui_stp3],cos(3.*(3.14*ui_wc)/(ui_kk<<1)),sin(3.*(3.14*ui_wc)/(ui_kk<<1)),tmpRe3,tmpIm3); */

                            ifn_NCMUL
                            (
                                &tmpRe1,
                                &tmpIm1,

                                a_Re[ui_stp1],
                                a_Im[ui_stp1],
                                //     cos((3.14*ui_wc)/(ui_kk<<1)),
                                // -1.*sin((3.14*ui_wc)/(ui_kk<<1))
                                acc_ct[ui_i],
                                acc_st[ui_i]
                            );

                            ifn_NCMUL
                            (
                                &tmpRe2,
                                &tmpIm2,

                                a_Re[ui_stp2],
                                a_Im[ui_stp2],
                                //     cos(2.*(3.14*ui_wc)/(ui_kk<<1)),
                                // -1.*sin(2.*(3.14*ui_wc)/(ui_kk<<1))
                                acc_ct[2*ui_i],
                                acc_st[2*ui_i]
                            );

                            ifn_NCMUL
                            (
                                &tmpRe3,
                                &tmpIm3,

                                a_Re[ui_stp3],
                                a_Im[ui_stp3],
                                //     cos(3.*(3.14*ui_wc)/(ui_kk<<1)),
                                // -1.*sin(3.*(3.14*ui_wc)/(ui_kk<<1))
                                acc_ct[3*ui_i],
                                acc_st[3*ui_i]
                            );

                            re_buff0      = tmpRe0 + tmpRe2;
                            im_buff0      = tmpIm0 + tmpIm2;

                            re_buff1      = tmpRe1 + tmpRe3;
                            im_buff1      = tmpIm1 + tmpIm3;

                            re_buff2      = tmpRe0 - tmpRe2;
                            im_buff2      = tmpIm0 - tmpIm2;

                            re_buff3      = tmpIm1 - tmpIm3;
                            im_buff3      = tmpRe3 - tmpRe1;


                            a_Re[ui_stp0] = re_buff0 + re_buff1;
                            a_Im[ui_stp0] = im_buff0 + im_buff1;

                            a_Re[ui_stp1] = re_buff2 + re_buff3;
                            a_Im[ui_stp1] = im_buff2 + im_buff3;

                            a_Re[ui_stp2] = re_buff0 - re_buff1;
                            a_Im[ui_stp2] = im_buff0 - im_buff1;

                            a_Re[ui_stp3] = re_buff2 - re_buff3;
                            a_Im[ui_stp3] = im_buff2 - im_buff3;
                        }
                }
        }
    return;
}
/*============================================================================*/
void fn_a_FFT4
(
    float *a_Re,
    float *a_Im
)
/*----------------------------------------------------------------------------*/
{
    float
    re_temp0,
    im_temp0,
    re_temp1,
    im_temp1,
    re_temp2,
    im_temp2,
    re_temp3,
    im_temp3;

    re_temp0 = a_Re[0] + a_Re[2];
    im_temp0 = a_Im[0] + a_Im[2];

    re_temp1 = a_Re[1] + a_Re[3];
    im_temp1 = a_Im[1] + a_Im[3];

    re_temp2 = a_Re[0] - a_Re[2];
    im_temp2 = a_Im[0] - a_Im[2];

    re_temp3 = a_Im[1] - a_Im[3];
    im_temp3 = a_Re[3] - a_Re[1];

    a_Re[0] = re_temp0 + re_temp1;
    a_Im[0] = im_temp0 + im_temp1;

    a_Re[1] = re_temp2 + re_temp3;
    a_Im[1] = im_temp2 + im_temp3;

    a_Re[2] = re_temp0 - re_temp1;
    a_Im[2] = im_temp0 - im_temp1;

    a_Re[3] = re_temp2 - re_temp3;
    a_Im[3] = im_temp2 - im_temp3;

    return;
}

