/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_MSVCCOMPATIBILITY_H
#define READUCT_MSVCCOMPATIBILITY_H

#if defined(_MSC_VER)
#  define SCINE_DLLEXPORT __declspec(dllexport)
#else
#  define SCINE_DLLEXPORT
#endif

#endif // READUCT_MSVCCOMPATIBILITY_H
