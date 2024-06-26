#include <cmath>

#define m_PI       3.14159265358979323846

template<typename m_T>
inline m_T const m_max(const m_T& _l, const m_T& _r){ return (_l>_r)?_l:_r;}

template<typename m_T>
inline m_T const m_min(const m_T& _l, const m_T& _r){ return (_l<_r)?_l:_r;}