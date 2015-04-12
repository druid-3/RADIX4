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

#include "FFTwrap.h"

/*============================================================================*/
void fn_a_invFFTdrSciling
(
    float                    *a_Re,
    float                    *a_Im,

    const unsigned int const  cuic_N
)
/*----------------------------------------------------------------------------*/
{
    unsigned int ui_o;

    for (ui_o=0; ui_o<cuic_N; ui_o++)
        {
            a_Re[ui_o] = a_Re[ui_o]/(float)cuic_N;
            a_Im[ui_o] = a_Im[ui_o]/(float)cuic_N;
        }
    return;
}
/*============================================================================*/
void fn_a_complexF0xF1
(
    float                    *a_ReOut,
    float                    *a_ImOut,

    float                    *a_Re0,
    float                    *a_Im0,
    float                    *a_Re1,
    float                    *a_Im1,

    const unsigned int const  cuic_N
)
/*----------------------------------------------------------------------------*/
{
    unsigned int ui_o;

    for (ui_o=0; ui_o<cuic_N; ui_o++)
        {
            a_ReOut[ui_o] = a_Re0[ui_o]*a_Re1[ui_o] - a_Im0[ui_o]*a_Im1[ui_o];
            a_ImOut[ui_o] = a_Im0[ui_o]*a_Re1[ui_o] + a_Im1[ui_o]*a_Re0[ui_o];
        }
    return;
}
/*============================================================================*/
void fn_a_F0xF1
(
    float                    *a_s2,

    float                    *a_s0,
    float                    *a_s1,

    const unsigned int const  cuic_N
)
/*----------------------------------------------------------------------------*/
{
    unsigned int ui_o;

    for (ui_o=0; ui_o<cuic_N; ui_o++)
        {
            a_s2[ui_o] = a_s0[ui_o]*a_s1[ui_o];
        }
    return;
}
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_dirFFT4
(
    float        *p_i,
    float        *p_q,

    unsigned int  ui_N,

    float        *p_CT,
    float        *p_ST,
    unsigned int  ui_S
)
/* -------------------------------------------------------------------------- */
{
    fn_bRev4(p_i, p_q, ui_N);
    fn_aT_TR4FFTkern(p_i, p_q, ui_N, p_CT, p_ST, ui_S);
}
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_invFFT4
(
    float        *p_i,
    float        *p_q,

    unsigned int  ui_N,

    float        *p_CT,
    float        *p_ST,
    unsigned int  ui_S
)
/* -------------------------------------------------------------------------- */
{
    fn_bRev4(p_i, p_q, ui_N);
    fn_aT_TR4FFTkern(p_q, p_i, ui_N, p_CT, p_ST, ui_S);
    fn_a_invFFTdrSciling(p_q, p_i, ui_N);
}
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_FASTCONV4
(
    float        *p_iOut,
    float        *p_qOut,

    float        *p_i0,
    float        *p_q0,
    float        *p_i1,
    float        *p_q1,

    unsigned int  ui_N,

    float        *p_CT,
    float        *p_ST,
    unsigned int  ui_S
)
/* -------------------------------------------------------------------------- */
{
    fn_aT_TR4FFTkern(p_i0, p_q0, ui_N, p_CT, p_ST, ui_S);
    fn_aT_TR4FFTkern(p_i1, p_q1, ui_N, p_CT, p_ST, ui_S);
    fn_a_complexF0xF1(p_iOut, p_qOut, p_i0, p_q0, p_i1, p_q1, ui_N);
    fn_aT_TR4FFTkern(p_qOut, p_iOut, ui_N, p_CT, p_ST, ui_S);
    fn_a_invFFTdrSciling(p_iOut, p_qOut, ui_N);
}
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_dirFFT2
(
    float        *p_i,
    float        *p_q,

    unsigned int  ui_N,

    float        *p_T,
    unsigned int  ui_S
)
/* -------------------------------------------------------------------------- */
{
    fn_bRev2(p_i, p_q, ui_N);
    fn_aT_TR2FFTkern(p_i, p_q, ui_N, p_T, ui_S);
}
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_invFFT2
(
    float        *p_i,
    float        *p_q,

    unsigned int  ui_N,

    float        *p_T,
    unsigned int  ui_S
)
/* -------------------------------------------------------------------------- */
{
    fn_bRev2(p_i, p_q, ui_N);
    fn_aT_TR2FFTkern(p_q, p_i, ui_N, p_T, ui_S);
    fn_a_invFFTdrSciling(p_q, p_i, ui_N);
}
/* iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii */
inline void ifn_FASTCONV2
(
    float        *p_iOut,
    float        *p_qOut,

    float        *p_i0,
    float        *p_q0,
    float        *p_i1,
    float        *p_q1,

    unsigned int  ui_N,

    float        *p_T,
    unsigned int  ui_S
)
/* -------------------------------------------------------------------------- */
{
    fn_aT_TR2FFTkern(p_i0, p_q0, ui_N, p_T, ui_S);
    fn_aT_TR2FFTkern(p_i1, p_q1, ui_N, p_T, ui_S);
    fn_a_complexF0xF1(p_iOut, p_qOut, p_i0, p_q0, p_i1, p_q1, ui_N);
    fn_aT_TR2FFTkern(p_qOut, p_iOut, ui_N, p_T, ui_S);
    fn_a_invFFTdrSciling(p_iOut, p_qOut, ui_N);
}

