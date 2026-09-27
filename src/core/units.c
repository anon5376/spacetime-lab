#include "stc/units.h"

double stc_schwarzschild_radius(double m_kg)
{
    double c = STC_C;
    return 2.0 * STC_G * m_kg / (c * c);
}

double stc_geometric_M(double m_kg)
{
    double c = STC_C;
    return STC_G * m_kg / (c * c);
}

double stc_arcsec(void)
{
    return STC_PI / (180.0 * 3600.0);
}

double stc_microarcsec(void)
{
    return stc_arcsec() / 1000000.0;
}
