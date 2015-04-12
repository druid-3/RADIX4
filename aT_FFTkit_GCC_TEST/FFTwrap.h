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

#ifndef _C_FFTWRAP
#define _C_FFTWRAP

#define DIRFFT4(i, q, N, CT, ST, S)                                            \
        fn_bRev4((i), (q), (N));                                               \
        fn_aT_TR4FFTkern((i), (q), (N), (CT), (ST), (S));

#define INVFFT4(i, q, N, CT, ST, S)                                            \
        fn_bRev4((i), (q), (N));                                               \
        fn_aT_TR4FFTkern((q), (i), (N), (CT), (ST), (S));                      \
        fn_a_invFFTdrSciling((q), (i), (N));

#define FASTCONV4(i2, q2, i0, q0, i1, q1, N, CT, ST, S)                        \
        fn_aT_TR4FFTkern((i0), (q0), (N), (CT), (ST), (S));                    \
        fn_aT_TR4FFTkern((i1), (q1), (N), (CT), (ST), (S));                    \
        fn_a_complexF0xF1((i2), (q2), (i0), (q0), (i1), (q1), (N));            \
        fn_aT_TR4FFTkern((q2), (i2), (N), (CT), (ST), (S));                    \
        fn_a_invFFTdrSciling((i2), (q2), (N));

#define DIRFFT2(i, q, N, T, S)                                                 \
        fn_bRev2((i), (q), (N));                                               \
        fn_aT_TR2FFTkern((i), (q), (N), (T), (S));

#define INVFFT2(i, q, N, T, S)                                                 \
        fn_bRev2((i), (q), (N));                                               \
        fn_aT_TR2FFTkern((q), (i), (N), (T), (S));                             \
        fn_a_invFFTdrSciling((q), (i), (N));

#define FASTCONV2(i2, q2, i0, q0, i1, q1, N, T, S)                             \
        fn_aT_TR2FFTkern((i0), (q0), (N), (T), (S));                           \
        fn_aT_TR2FFTkern((i1), (q1), (N), (T), (S));                           \
        fn_a_complexF0xF1((i2), (q2), (i0), (q0), (i1), (q1), (N));            \
        fn_aT_TR2FFTkern((q2), (i2), (N), (T), (S));                           \
        fn_a_invFFTdrSciling((i2), (q2), (N));

/*============================================================================*/
void fn_a_invFFTdrSciling
(
    float                    *,
    float                    *,
    const unsigned int const
);
/*============================================================================*/
void fn_a_complexF0xF1
(
    float                    *,
    float                    *,

    float                    *,
    float                    *,
    float                    *,
    float                    *,

    const unsigned int const
);
/*============================================================================*/
void fn_a_F0xF1
(
    float                    *,

    float                    *,
    float                    *,
    const unsigned int const
);
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_dirFFT4
(
    float        *,
    float        *,
    unsigned int  ,
    float        *,
    float        *,
    unsigned int
);
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_invFFT4
(
    float        *,
    float        *,
    unsigned int  ,
    float        *,
    float        *,
    unsigned int
);
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_FASTCONV4
(
    float        *,
    float        *,

    float        *,
    float        *,
    float        *,
    float        *,

    unsigned int  ,

    float        *,
    float        *,
    unsigned int
);
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_dirFFT2
(
    float        *,
    float        *,

    unsigned int  ,

    float        *,
    unsigned int
);
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_invFFT2
(
    float        *,
    float        *,

    unsigned int  ,
    float        *,
    unsigned int
);
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_FASTCONV2
(
    float        *,
    float        *,

    float        *,
    float        *,
    float        *,
    float        *,

    unsigned int  ,
    float        *,
    unsigned int
);
/* -------------------------------------------------------------------------- */

#endif
