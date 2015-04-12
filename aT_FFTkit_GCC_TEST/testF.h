#ifndef _C_INOUTDATAPRINT
#define _C_INOUTDATAPRINT

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define df_2PI 6.283185
#define df_PI  3.141593

/*============================================================================*/
void fn_aPrint
(
    float                    *,
    const unsigned int const
);
/*============================================================================*/
void fn_Test_DDS_C
(
    float                    *,
    const unsigned int const  ,
    float                     ,
    bool
);
/*============================================================================*/
void fn_Test_DDS_S
(
    float                    *,
    const unsigned int const  ,
    float                     ,
    bool
);
/*----------------------------------------------------------------------------*/

#endif
