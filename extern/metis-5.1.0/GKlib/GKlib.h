/*
 * GKlib.h
 * 
 * George's library of most frequently used routines
 *
 * $Id: GKlib.h 13005 2012-10-23 22:34:36Z karypis $
 *
 */

#ifndef _GKLIB_H_
#define _GKLIB_H_ 1

#define GKMSPACE

#if defined(_MSC_VER)
#define __MSC__
#endif
#if defined(__ICC)
#define __ICC__
#endif


#include "gk_arch.h" /*!< This should be here, prior to the includes */


/*************************************************************************
* Header file inclusion section
**************************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <signal.h>
#include <setjmp.h>
#include <assert.h>
#include <sys/stat.h>

#if defined(__WITHPCRE__)
  #include <pcreposix.h>
#else
  #if defined(USE_GKREGEX)
    #include "gkregex.h"
  #else
    #include <regex.h>
  #endif /* defined(USE_GKREGEX) */
#endif /* defined(__WITHPCRE__) */



#if defined(__OPENMP__) 
#include <omp.h>
#endif




#include <gk_types.h>
#include <gk_struct.h>
#include <gk_externs.h>
#include <gk_defs.h>
#include <gk_macros.h>
#include <gk_getopt.h>

#include <gk_mksort.h>
#include <gk_mkblas.h>
#include <gk_mkmemory.h>
#include <gk_mkpqueue.h>
#include <gk_mkpqueue2.h>
#include <gk_mkrandom.h>
#include <gk_mkutils.h>

/* Private namespacing, matching libmetis/rename.h. Everything else
 * GKlib exports already carries its own ``gk_`` prefix; these are the
 * names that do not, so without the renames below this private copy
 * would define them as bare globals and could bind against another
 * GKlib/METIS in the same link. Placed ahead of gk_proto.h so both the
 * declarations and the definitions in GKlib/*.c see them. */
#define ComputeAccuracy                  gk_MUMPS_ComputeAccuracy
#define ComputeMean                      gk_MUMPS_ComputeMean
#define ComputeMedianRFP                 gk_MUMPS_ComputeMedianRFP
#define ComputeROCn                      gk_MUMPS_ComputeROCn
#define ComputeStdDev                    gk_MUMPS_ComputeStdDev
#define decodeblock                      gk_MUMPS_decodeblock
#define encodeblock                      gk_MUMPS_encodeblock
#define errexit                          gk_MUMPS_errexit
#define getpathname                      gk_MUMPS_getpathname
#define GKDecodeBase64                   gk_MUMPS_GKDecodeBase64
#define GKEncodeBase64                   gk_MUMPS_GKEncodeBase64
#define gkfooo                           gk_MUMPS_gkfooo
#define HTable_Create                    gk_MUMPS_HTable_Create
#define HTable_Delete                    gk_MUMPS_HTable_Delete
#define HTable_Destroy                   gk_MUMPS_HTable_Destroy
#define HTable_GetNext                   gk_MUMPS_HTable_GetNext
#define HTable_HFunction                 gk_MUMPS_HTable_HFunction
#define HTable_Insert                    gk_MUMPS_HTable_Insert
#define HTable_Reset                     gk_MUMPS_HTable_Reset
#define HTable_Resize                    gk_MUMPS_HTable_Resize
#define HTable_Search                    gk_MUMPS_HTable_Search
#define HTable_SearchAndDelete           gk_MUMPS_HTable_SearchAndDelete
#define itemsets_find_frequent_itemsets  gk_MUMPS_itemsets_find_frequent_itemsets
#define itemsets_project_matrix          gk_MUMPS_itemsets_project_matrix
#define PrintBackTrace                   gk_MUMPS_PrintBackTrace

#include <gk_proto.h>


#endif  /* GKlib.h */


