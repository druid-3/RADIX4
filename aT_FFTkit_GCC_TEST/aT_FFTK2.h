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

#ifndef _C_FFTK2
#define _C_FFTK2

#include "testF.h"

/*============================================================================*/
static inline void ifn_NCMUL
(
    float *,
    float *,

    float  ,
    float  ,
    float  ,
    float
);
/*============================================================================*/
static inline void ifn_CMUL
(
    float *,
    float *,

    float  ,
    float  ,
    float  ,
    float
);
/*============================================================================*/
void fn_bRev2
(
    float                    *,
    float                    *,
    const unsigned int const
);
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
    float                    *,
    float                    *,
    const unsigned int const ,
    const float const        *,
    const unsigned int const
);

#endif /* __FFTK2_H__ */
