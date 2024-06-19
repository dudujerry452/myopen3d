#ifndef _M_UTIL_H
#define _M_UTIL_H

#include <cmath>

#define m_PI       3.14159265358979323846
#define m_E        2.71828182845904523536028747135266250

template<typename m_T>
inline m_T const m_max(const m_T& _l, const m_T& _r){ return (_l>_r)?_l:_r;}

template<typename m_T>
inline m_T const m_min(const m_T& _l, const m_T& _r){ return (_l<_r)?_l:_r;}

#define chnow() std::chrono::steady_clock::now()
#define mildiff(x) std::chrono::duration_cast<std::chrono::milliseconds>(x).count()
#endif