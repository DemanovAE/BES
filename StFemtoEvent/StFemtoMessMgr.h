#ifndef StFemtoMessMgr_h
#define StFemtoMessMgr_h

#if defined(_VANILLA_ROOT_)
#include "Riostream.h"
#define LOG_INFO  std::cout
#define LOG_DEBUG std::cout
#define LOG_WARN  std::cout
#define LOG_ERROR std::cout
#define endm      std::endl
#else
#include "St_base/StMessMgr.h"
#endif

#endif
