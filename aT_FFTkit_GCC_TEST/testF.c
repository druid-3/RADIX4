#include "testF.h"

/*============================================================================*/
void fn_aPrint
(
    float                    a_a[],
    const unsigned int const ui_N
)
/*----------------------------------------------------------------------------*/
{
    unsigned int ui_ii;

    printf(" DATA OK ---------------------------------------------------------------------- \n");

    for(ui_ii=0; ui_ii<ui_N; ui_ii++)
        {
            printf("[%4i] -> %f\n", ui_ii, a_a[ui_ii]);
        }

    printf("\n");
}
/*============================================================================*/
void fn_Test_DDS_C
(
    float                    *a_s,
    const unsigned int const  cuic_N,
    float                     garm,
    bool                      b_reset
)
/*----------------------------------------------------------------------------*/
{
    unsigned int ui_jj;
    float phase = 0.;
    static float oldPhase = 0.;

    printf(" GEN OK ---------------------------------------------------------------------- \n");

    if(b_reset)
        {
            memset(a_s, 0, cuic_N*sizeof(float));
        }

    for(ui_jj=0; ui_jj<cuic_N; ui_jj++)
        {
            phase = (((float)ui_jj)*garm*df_2PI + oldPhase)/((float)cuic_N);
            a_s[ui_jj] += sin(phase);
        }

    oldPhase =  remainder(phase, (df_2PI/garm));

    return;
}
/*============================================================================*/
void fn_Test_DDS_S
(
    float                    *a_s,
    const unsigned int const  cuic_N,
    float                     garm,
    bool                      b_reset
)
/*----------------------------------------------------------------------------*/
{
    unsigned int ui_jj;
    float phase = 0.;
    static float oldPhase = 0.;

    printf(" GEN OK ---------------------------------------------------------------------- \n");

    if(b_reset)
        {
            memset(a_s, 0, cuic_N*sizeof(float));
        }

    for(ui_jj=0; ui_jj<cuic_N; ui_jj++)
        {
            phase = (((float)ui_jj)*garm*df_2PI + oldPhase)/((float)cuic_N);
            a_s[ui_jj] += cos(phase);
        }

    oldPhase =  remainder(phase, (df_2PI/garm));

    return;
}
