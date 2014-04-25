#ifndef FILTER_GLOBAL_H
#define FILTER_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(FILTER_LIBRARY)
#  define FILTERSHARED_EXPORT Q_DECL_EXPORT
#else
#  define FILTERSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // FILTER_GLOBAL_H
