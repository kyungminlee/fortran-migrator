/*
 * GKlib.h
 * 
 * George's library of most frequently used routines
 *
 * $Id: GKlib.h 14866 2013-08-03 16:40:04Z karypis $
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
#include <memory.h>
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

#if defined(USE_PCRE) && defined(HAVE_PCREPOSIX_H)
  #include <pcreposix.h>
#elif defined(HAVE_REGEX_H)
  #include <regex.h>
#else
  #include "gkregex.h"
#endif


#if defined(__OPENMP__) 
#include <omp.h>
#endif





/* Private namespacing, matching libmetis/rename.h: every GKlib
   global gets MUMPS_ prepended, so this copy cannot clash with
   another GKlib in the link — gk_ is GKlib's prefix, not ours.
   Generated — see scripts/vendor_metis.sh. */
#define ComputeAccuracy                          MUMPS_ComputeAccuracy
#define ComputeMean                              MUMPS_ComputeMean
#define ComputeMedianRFP                         MUMPS_ComputeMedianRFP
#define ComputeROCn                              MUMPS_ComputeROCn
#define ComputeStdDev                            MUMPS_ComputeStdDev
#define decodeblock                              MUMPS_decodeblock
#define encodeblock                              MUMPS_encodeblock
#define errexit                                  MUMPS_errexit
#define getpathname                              MUMPS_getpathname
#define gk_AllocMatrix                           MUMPS_gk_AllocMatrix
#define gk_array2csr                             MUMPS_gk_array2csr
#define gk_cacheCreate                           MUMPS_gk_cacheCreate
#define gk_cacheDestroy                          MUMPS_gk_cacheDestroy
#define gk_cacheGetHitRate                       MUMPS_gk_cacheGetHitRate
#define gk_cacheLoad                             MUMPS_gk_cacheLoad
#define gk_cacheReset                            MUMPS_gk_cacheReset
#define gk_cAllocMatrix                          MUMPS_gk_cAllocMatrix
#define gk_cargmax                               MUMPS_gk_cargmax
#define gk_cargmax_n                             MUMPS_gk_cargmax_n
#define gk_cargmin                               MUMPS_gk_cargmin
#define gk_caxpy                                 MUMPS_gk_caxpy
#define gk_ccopy                                 MUMPS_gk_ccopy
#define gk_cdot                                  MUMPS_gk_cdot
#define gk_cFreeMatrix                           MUMPS_gk_cFreeMatrix
#define gk_cincset                               MUMPS_gk_cincset
#define gk_ckvAllocMatrix                        MUMPS_gk_ckvAllocMatrix
#define gk_ckvcopy                               MUMPS_gk_ckvcopy
#define gk_ckvFreeMatrix                         MUMPS_gk_ckvFreeMatrix
#define gk_ckvmalloc                             MUMPS_gk_ckvmalloc
#define gk_ckvrealloc                            MUMPS_gk_ckvrealloc
#define gk_ckvset                                MUMPS_gk_ckvset
#define gk_ckvSetMatrix                          MUMPS_gk_ckvSetMatrix
#define gk_ckvsmalloc                            MUMPS_gk_ckvsmalloc
#define gk_ckvsortd                              MUMPS_gk_ckvsortd
#define gk_ckvsorti                              MUMPS_gk_ckvsorti
#define gk_cmalloc                               MUMPS_gk_cmalloc
#define gk_cmax                                  MUMPS_gk_cmax
#define gk_cmin                                  MUMPS_gk_cmin
#define gk_cnorm2                                MUMPS_gk_cnorm2
#define gk_CPUSeconds                            MUMPS_gk_CPUSeconds
#define gk_crand                                 MUMPS_gk_crand
#define gk_crandArrayPermute                     MUMPS_gk_crandArrayPermute
#define gk_crandArrayPermuteFine                 MUMPS_gk_crandArrayPermuteFine
#define gk_crandInRange                          MUMPS_gk_crandInRange
#define gk_creadfilebin                          MUMPS_gk_creadfilebin
#define gk_crealloc                              MUMPS_gk_crealloc
#define gk_cscale                                MUMPS_gk_cscale
#define gk_cset                                  MUMPS_gk_cset
#define gk_cSetMatrix                            MUMPS_gk_cSetMatrix
#define gk_csmalloc                              MUMPS_gk_csmalloc
#define gk_csortd                                MUMPS_gk_csortd
#define gk_csorti                                MUMPS_gk_csorti
#define gk_csrand                                MUMPS_gk_csrand
#define gk_csr_CompactColumns                    MUMPS_gk_csr_CompactColumns
#define gk_csr_ComputeBestFOrderingSymmetric     MUMPS_gk_csr_ComputeBestFOrderingSymmetric
#define gk_csr_ComputeBFSOrderingSymmetric       MUMPS_gk_csr_ComputeBFSOrderingSymmetric
#define gk_csr_ComputeNorms                      MUMPS_gk_csr_ComputeNorms
#define gk_csr_ComputePairSimilarity             MUMPS_gk_csr_ComputePairSimilarity
#define gk_csr_ComputeSimilarity                 MUMPS_gk_csr_ComputeSimilarity
#define gk_csr_ComputeSquaredNorms               MUMPS_gk_csr_ComputeSquaredNorms
#define gk_csr_ComputeSums                       MUMPS_gk_csr_ComputeSums
#define gk_csr_Create                            MUMPS_gk_csr_Create
#define gk_csr_CreateIndex                       MUMPS_gk_csr_CreateIndex
#define gk_csr_DetermineFormat                   MUMPS_gk_csr_DetermineFormat
#define gk_csr_Dup                               MUMPS_gk_csr_Dup
#define gk_csr_ExtractPartition                  MUMPS_gk_csr_ExtractPartition
#define gk_csr_ExtractRows                       MUMPS_gk_csr_ExtractRows
#define gk_csr_ExtractSubmatrix                  MUMPS_gk_csr_ExtractSubmatrix
#define gk_csr_FindConnectedComponents           MUMPS_gk_csr_FindConnectedComponents
#define gk_csr_Free                              MUMPS_gk_csr_Free
#define gk_csr_FreeContents                      MUMPS_gk_csr_FreeContents
#define gk_csr_GetSimilarRows                    MUMPS_gk_csr_GetSimilarRows
#define gk_csr_Init                              MUMPS_gk_csr_Init
#define gk_csr_LowFilter                         MUMPS_gk_csr_LowFilter
#define gk_csr_MakeSymmetric                     MUMPS_gk_csr_MakeSymmetric
#define gk_csr_Normalize                         MUMPS_gk_csr_Normalize
#define gk_csr_Prune                             MUMPS_gk_csr_Prune
#define gk_csr_Read                              MUMPS_gk_csr_Read
#define gk_csr_ReorderSymmetric                  MUMPS_gk_csr_ReorderSymmetric
#define gk_csr_Scale                             MUMPS_gk_csr_Scale
#define gk_csr_Shuffle                           MUMPS_gk_csr_Shuffle
#define gk_csr_SortIndices                       MUMPS_gk_csr_SortIndices
#define gk_csr_Split                             MUMPS_gk_csr_Split
#define gk_csr_TopKPlusFilter                    MUMPS_gk_csr_TopKPlusFilter
#define gk_csr_Transpose                         MUMPS_gk_csr_Transpose
#define gk_csr_Write                             MUMPS_gk_csr_Write
#define gk_csr_ZScoreFilter                      MUMPS_gk_csr_ZScoreFilter
#define gk_csum                                  MUMPS_gk_csum
#define gk_cur_jbufs                             MUMPS_gk_cur_jbufs
#define gk_cwritefilebin                         MUMPS_gk_cwritefilebin
#define gk_dAllocMatrix                          MUMPS_gk_dAllocMatrix
#define gk_dargmax                               MUMPS_gk_dargmax
#define gk_dargmax_n                             MUMPS_gk_dargmax_n
#define gk_dargmin                               MUMPS_gk_dargmin
#define gk_daxpy                                 MUMPS_gk_daxpy
#define gk_dcopy                                 MUMPS_gk_dcopy
#define gk_ddot                                  MUMPS_gk_ddot
#define GKDecodeBase64                           MUMPS_GKDecodeBase64
#define gk_dexists                               MUMPS_gk_dexists
#define gk_dfkvkselect                           MUMPS_gk_dfkvkselect
#define gk_dFreeMatrix                           MUMPS_gk_dFreeMatrix
#define gk_dincset                               MUMPS_gk_dincset
#define gk_dkvAllocMatrix                        MUMPS_gk_dkvAllocMatrix
#define gk_dkvcopy                               MUMPS_gk_dkvcopy
#define gk_dkvFreeMatrix                         MUMPS_gk_dkvFreeMatrix
#define gk_dkvmalloc                             MUMPS_gk_dkvmalloc
#define gk_dkvrealloc                            MUMPS_gk_dkvrealloc
#define gk_dkvset                                MUMPS_gk_dkvset
#define gk_dkvSetMatrix                          MUMPS_gk_dkvSetMatrix
#define gk_dkvsmalloc                            MUMPS_gk_dkvsmalloc
#define gk_dkvsortd                              MUMPS_gk_dkvsortd
#define gk_dkvsorti                              MUMPS_gk_dkvsorti
#define gk_dmalloc                               MUMPS_gk_dmalloc
#define gk_dmax                                  MUMPS_gk_dmax
#define gk_dmin                                  MUMPS_gk_dmin
#define gk_dnorm2                                MUMPS_gk_dnorm2
#define gk_dpqCheckHeap                          MUMPS_gk_dpqCheckHeap
#define gk_dpqCreate                             MUMPS_gk_dpqCreate
#define gk_dpqDelete                             MUMPS_gk_dpqDelete
#define gk_dpqDestroy                            MUMPS_gk_dpqDestroy
#define gk_dpqFree                               MUMPS_gk_dpqFree
#define gk_dpqGetTop                             MUMPS_gk_dpqGetTop
#define gk_dpqInit                               MUMPS_gk_dpqInit
#define gk_dpqInsert                             MUMPS_gk_dpqInsert
#define gk_dpqLength                             MUMPS_gk_dpqLength
#define gk_dpqReset                              MUMPS_gk_dpqReset
#define gk_dpqSeeKey                             MUMPS_gk_dpqSeeKey
#define gk_dpqSeeTopKey                          MUMPS_gk_dpqSeeTopKey
#define gk_dpqSeeTopVal                          MUMPS_gk_dpqSeeTopVal
#define gk_dpqUpdate                             MUMPS_gk_dpqUpdate
#define gk_drand                                 MUMPS_gk_drand
#define gk_drandArrayPermute                     MUMPS_gk_drandArrayPermute
#define gk_drandArrayPermuteFine                 MUMPS_gk_drandArrayPermuteFine
#define gk_drandInRange                          MUMPS_gk_drandInRange
#define gk_dreadfilebin                          MUMPS_gk_dreadfilebin
#define gk_drealloc                              MUMPS_gk_drealloc
#define gk_dscale                                MUMPS_gk_dscale
#define gk_dset                                  MUMPS_gk_dset
#define gk_dSetMatrix                            MUMPS_gk_dSetMatrix
#define gk_dsmalloc                              MUMPS_gk_dsmalloc
#define gk_dsortd                                MUMPS_gk_dsortd
#define gk_dsorti                                MUMPS_gk_dsorti
#define gk_dsrand                                MUMPS_gk_dsrand
#define gk_dsum                                  MUMPS_gk_dsum
#define gk_dwritefilebin                         MUMPS_gk_dwritefilebin
#define GKEncodeBase64                           MUMPS_GKEncodeBase64
#define gk_errexit                               MUMPS_gk_errexit
#define gk_fAllocMatrix                          MUMPS_gk_fAllocMatrix
#define gk_fargmax                               MUMPS_gk_fargmax
#define gk_fargmax_n                             MUMPS_gk_fargmax_n
#define gk_fargmin                               MUMPS_gk_fargmin
#define gk_faxpy                                 MUMPS_gk_faxpy
#define gk_fclose                                MUMPS_gk_fclose
#define gk_fcopy                                 MUMPS_gk_fcopy
#define gk_fdot                                  MUMPS_gk_fdot
#define gk_fexists                               MUMPS_gk_fexists
#define gk_fFreeMatrix                           MUMPS_gk_fFreeMatrix
#define gk_fincset                               MUMPS_gk_fincset
#define gk_find_frequent_itemsets                MUMPS_gk_find_frequent_itemsets
#define gk_fkvAllocMatrix                        MUMPS_gk_fkvAllocMatrix
#define gk_fkvcopy                               MUMPS_gk_fkvcopy
#define gk_fkvFreeMatrix                         MUMPS_gk_fkvFreeMatrix
#define gk_fkvmalloc                             MUMPS_gk_fkvmalloc
#define gk_fkvrealloc                            MUMPS_gk_fkvrealloc
#define gk_fkvset                                MUMPS_gk_fkvset
#define gk_fkvSetMatrix                          MUMPS_gk_fkvSetMatrix
#define gk_fkvsmalloc                            MUMPS_gk_fkvsmalloc
#define gk_fkvsortd                              MUMPS_gk_fkvsortd
#define gk_fkvsorti                              MUMPS_gk_fkvsorti
#define gk_flog2                                 MUMPS_gk_flog2
#define gk_fmalloc                               MUMPS_gk_fmalloc
#define gk_fmax                                  MUMPS_gk_fmax
#define gk_fmin                                  MUMPS_gk_fmin
#define gk_fnorm2                                MUMPS_gk_fnorm2
#define gkfooo                                   MUMPS_gkfooo
#define gk_fopen                                 MUMPS_gk_fopen
#define gk_fpqCheckHeap                          MUMPS_gk_fpqCheckHeap
#define gk_fpqCreate                             MUMPS_gk_fpqCreate
#define gk_fpqDelete                             MUMPS_gk_fpqDelete
#define gk_fpqDestroy                            MUMPS_gk_fpqDestroy
#define gk_fpqFree                               MUMPS_gk_fpqFree
#define gk_fpqGetTop                             MUMPS_gk_fpqGetTop
#define gk_fpqInit                               MUMPS_gk_fpqInit
#define gk_fpqInsert                             MUMPS_gk_fpqInsert
#define gk_fpqLength                             MUMPS_gk_fpqLength
#define gk_fpqReset                              MUMPS_gk_fpqReset
#define gk_fpqSeeKey                             MUMPS_gk_fpqSeeKey
#define gk_fpqSeeTopKey                          MUMPS_gk_fpqSeeTopKey
#define gk_fpqSeeTopVal                          MUMPS_gk_fpqSeeTopVal
#define gk_fpqUpdate                             MUMPS_gk_fpqUpdate
#define gk_frand                                 MUMPS_gk_frand
#define gk_frandArrayPermute                     MUMPS_gk_frandArrayPermute
#define gk_frandArrayPermuteFine                 MUMPS_gk_frandArrayPermuteFine
#define gk_frandInRange                          MUMPS_gk_frandInRange
#define gk_freadfilebin                          MUMPS_gk_freadfilebin
#define gk_frealloc                              MUMPS_gk_frealloc
#define gk_free                                  MUMPS_gk_free
#define gk_FreeMatrix                            MUMPS_gk_FreeMatrix
#define gk_freetokenslist                        MUMPS_gk_freetokenslist
#define gk_fscale                                MUMPS_gk_fscale
#define gk_fset                                  MUMPS_gk_fset
#define gk_fSetMatrix                            MUMPS_gk_fSetMatrix
#define gk_fsmalloc                              MUMPS_gk_fsmalloc
#define gk_fsortd                                MUMPS_gk_fsortd
#define gk_fsorti                                MUMPS_gk_fsorti
#define gk_fsrand                                MUMPS_gk_fsrand
#define gk_fsum                                  MUMPS_gk_fsum
#define gk_fwritefilebin                         MUMPS_gk_fwritefilebin
#define gk_getbasename                           MUMPS_gk_getbasename
#define gk_GetCurMemoryUsed                      MUMPS_gk_GetCurMemoryUsed
#define gk_getextname                            MUMPS_gk_getextname
#define gk_getfilename                           MUMPS_gk_getfilename
#define gk_getfilestats                          MUMPS_gk_getfilestats
#define gk_getfsize                              MUMPS_gk_getfsize
#define gk_getline                               MUMPS_gk_getline
#define gk_GetMaxMemoryUsed                      MUMPS_gk_GetMaxMemoryUsed
#define gk_getopt                                MUMPS_gk_getopt
#define gk_getopt_initialized                    MUMPS_gk_getopt_initialized
#define gk_getopt_long                           MUMPS_gk_getopt_long
#define gk_getopt_long_only                      MUMPS_gk_getopt_long_only
#define gk_GetProcVmPeak                         MUMPS_gk_GetProcVmPeak
#define gk_GetStringID                           MUMPS_gk_GetStringID
#define gk_GetVMInfo                             MUMPS_gk_GetVMInfo
#define gk_gkmcoreAdd                            MUMPS_gk_gkmcoreAdd
#define gk_gkmcoreCreate                         MUMPS_gk_gkmcoreCreate
#define gk_gkmcoreDel                            MUMPS_gk_gkmcoreDel
#define gk_gkmcoreDestroy                        MUMPS_gk_gkmcoreDestroy
#define gk_gkmcorePop                            MUMPS_gk_gkmcorePop
#define gk_gkmcorePush                           MUMPS_gk_gkmcorePush
#define gk_graph_ComputeBestFOrdering            MUMPS_gk_graph_ComputeBestFOrdering
#define gk_graph_ComputeBestFOrdering0           MUMPS_gk_graph_ComputeBestFOrdering0
#define gk_graph_ComputeBFSOrdering              MUMPS_gk_graph_ComputeBFSOrdering
#define gk_graph_Create                          MUMPS_gk_graph_Create
#define gk_graph_Dup                             MUMPS_gk_graph_Dup
#define gk_graph_ExtractSubgraph                 MUMPS_gk_graph_ExtractSubgraph
#define gk_graph_FindComponents                  MUMPS_gk_graph_FindComponents
#define gk_graph_Free                            MUMPS_gk_graph_Free
#define gk_graph_FreeContents                    MUMPS_gk_graph_FreeContents
#define gk_graph_Init                            MUMPS_gk_graph_Init
#define gk_graph_MakeSymmetric                   MUMPS_gk_graph_MakeSymmetric
#define gk_graph_Read                            MUMPS_gk_graph_Read
#define gk_graph_Reorder                         MUMPS_gk_graph_Reorder
#define gk_graph_SingleSourceShortestPaths       MUMPS_gk_graph_SingleSourceShortestPaths
#define gk_graph_SortAdjacencies                 MUMPS_gk_graph_SortAdjacencies
#define gk_graph_Transpose                       MUMPS_gk_graph_Transpose
#define gk_graph_Write                           MUMPS_gk_graph_Write
#define gk_i16AllocMatrix                        MUMPS_gk_i16AllocMatrix
#define gk_i16copy                               MUMPS_gk_i16copy
#define gk_i16FreeMatrix                         MUMPS_gk_i16FreeMatrix
#define gk_i16kvAllocMatrix                      MUMPS_gk_i16kvAllocMatrix
#define gk_i16kvcopy                             MUMPS_gk_i16kvcopy
#define gk_i16kvFreeMatrix                       MUMPS_gk_i16kvFreeMatrix
#define gk_i16kvmalloc                           MUMPS_gk_i16kvmalloc
#define gk_i16kvrealloc                          MUMPS_gk_i16kvrealloc
#define gk_i16kvset                              MUMPS_gk_i16kvset
#define gk_i16kvSetMatrix                        MUMPS_gk_i16kvSetMatrix
#define gk_i16kvsmalloc                          MUMPS_gk_i16kvsmalloc
#define gk_i16malloc                             MUMPS_gk_i16malloc
#define gk_i16realloc                            MUMPS_gk_i16realloc
#define gk_i16set                                MUMPS_gk_i16set
#define gk_i16SetMatrix                          MUMPS_gk_i16SetMatrix
#define gk_i16smalloc                            MUMPS_gk_i16smalloc
#define gk_i2cc2i_create_common                  MUMPS_gk_i2cc2i_create_common
#define gk_i32AllocMatrix                        MUMPS_gk_i32AllocMatrix
#define gk_i32argmax                             MUMPS_gk_i32argmax
#define gk_i32argmax_n                           MUMPS_gk_i32argmax_n
#define gk_i32argmin                             MUMPS_gk_i32argmin
#define gk_i32axpy                               MUMPS_gk_i32axpy
#define gk_i32copy                               MUMPS_gk_i32copy
#define gk_i32dot                                MUMPS_gk_i32dot
#define gk_i32FreeMatrix                         MUMPS_gk_i32FreeMatrix
#define gk_i32incset                             MUMPS_gk_i32incset
#define gk_i32kvAllocMatrix                      MUMPS_gk_i32kvAllocMatrix
#define gk_i32kvcopy                             MUMPS_gk_i32kvcopy
#define gk_i32kvFreeMatrix                       MUMPS_gk_i32kvFreeMatrix
#define gk_i32kvmalloc                           MUMPS_gk_i32kvmalloc
#define gk_i32kvrealloc                          MUMPS_gk_i32kvrealloc
#define gk_i32kvset                              MUMPS_gk_i32kvset
#define gk_i32kvSetMatrix                        MUMPS_gk_i32kvSetMatrix
#define gk_i32kvsmalloc                          MUMPS_gk_i32kvsmalloc
#define gk_i32kvsortd                            MUMPS_gk_i32kvsortd
#define gk_i32kvsorti                            MUMPS_gk_i32kvsorti
#define gk_i32malloc                             MUMPS_gk_i32malloc
#define gk_i32max                                MUMPS_gk_i32max
#define gk_i32min                                MUMPS_gk_i32min
#define gk_i32norm2                              MUMPS_gk_i32norm2
#define gk_i32pqCheckHeap                        MUMPS_gk_i32pqCheckHeap
#define gk_i32pqCreate                           MUMPS_gk_i32pqCreate
#define gk_i32pqDelete                           MUMPS_gk_i32pqDelete
#define gk_i32pqDestroy                          MUMPS_gk_i32pqDestroy
#define gk_i32pqFree                             MUMPS_gk_i32pqFree
#define gk_i32pqGetTop                           MUMPS_gk_i32pqGetTop
#define gk_i32pqInit                             MUMPS_gk_i32pqInit
#define gk_i32pqInsert                           MUMPS_gk_i32pqInsert
#define gk_i32pqLength                           MUMPS_gk_i32pqLength
#define gk_i32pqReset                            MUMPS_gk_i32pqReset
#define gk_i32pqSeeKey                           MUMPS_gk_i32pqSeeKey
#define gk_i32pqSeeTopKey                        MUMPS_gk_i32pqSeeTopKey
#define gk_i32pqSeeTopVal                        MUMPS_gk_i32pqSeeTopVal
#define gk_i32pqUpdate                           MUMPS_gk_i32pqUpdate
#define gk_i32rand                               MUMPS_gk_i32rand
#define gk_i32randArrayPermute                   MUMPS_gk_i32randArrayPermute
#define gk_i32randArrayPermuteFine               MUMPS_gk_i32randArrayPermuteFine
#define gk_i32randInRange                        MUMPS_gk_i32randInRange
#define gk_i32readfile                           MUMPS_gk_i32readfile
#define gk_i32readfilebin                        MUMPS_gk_i32readfilebin
#define gk_i32realloc                            MUMPS_gk_i32realloc
#define gk_i32scale                              MUMPS_gk_i32scale
#define gk_i32set                                MUMPS_gk_i32set
#define gk_i32SetMatrix                          MUMPS_gk_i32SetMatrix
#define gk_i32smalloc                            MUMPS_gk_i32smalloc
#define gk_i32sortd                              MUMPS_gk_i32sortd
#define gk_i32sorti                              MUMPS_gk_i32sorti
#define gk_i32srand                              MUMPS_gk_i32srand
#define gk_i32sum                                MUMPS_gk_i32sum
#define gk_i32writefilebin                       MUMPS_gk_i32writefilebin
#define gk_i64AllocMatrix                        MUMPS_gk_i64AllocMatrix
#define gk_i64argmax                             MUMPS_gk_i64argmax
#define gk_i64argmax_n                           MUMPS_gk_i64argmax_n
#define gk_i64argmin                             MUMPS_gk_i64argmin
#define gk_i64axpy                               MUMPS_gk_i64axpy
#define gk_i64copy                               MUMPS_gk_i64copy
#define gk_i64dot                                MUMPS_gk_i64dot
#define gk_i64FreeMatrix                         MUMPS_gk_i64FreeMatrix
#define gk_i64incset                             MUMPS_gk_i64incset
#define gk_i64kvAllocMatrix                      MUMPS_gk_i64kvAllocMatrix
#define gk_i64kvcopy                             MUMPS_gk_i64kvcopy
#define gk_i64kvFreeMatrix                       MUMPS_gk_i64kvFreeMatrix
#define gk_i64kvmalloc                           MUMPS_gk_i64kvmalloc
#define gk_i64kvrealloc                          MUMPS_gk_i64kvrealloc
#define gk_i64kvset                              MUMPS_gk_i64kvset
#define gk_i64kvSetMatrix                        MUMPS_gk_i64kvSetMatrix
#define gk_i64kvsmalloc                          MUMPS_gk_i64kvsmalloc
#define gk_i64kvsortd                            MUMPS_gk_i64kvsortd
#define gk_i64kvsorti                            MUMPS_gk_i64kvsorti
#define gk_i64malloc                             MUMPS_gk_i64malloc
#define gk_i64max                                MUMPS_gk_i64max
#define gk_i64min                                MUMPS_gk_i64min
#define gk_i64norm2                              MUMPS_gk_i64norm2
#define gk_i64pqCheckHeap                        MUMPS_gk_i64pqCheckHeap
#define gk_i64pqCreate                           MUMPS_gk_i64pqCreate
#define gk_i64pqDelete                           MUMPS_gk_i64pqDelete
#define gk_i64pqDestroy                          MUMPS_gk_i64pqDestroy
#define gk_i64pqFree                             MUMPS_gk_i64pqFree
#define gk_i64pqGetTop                           MUMPS_gk_i64pqGetTop
#define gk_i64pqInit                             MUMPS_gk_i64pqInit
#define gk_i64pqInsert                           MUMPS_gk_i64pqInsert
#define gk_i64pqLength                           MUMPS_gk_i64pqLength
#define gk_i64pqReset                            MUMPS_gk_i64pqReset
#define gk_i64pqSeeKey                           MUMPS_gk_i64pqSeeKey
#define gk_i64pqSeeTopKey                        MUMPS_gk_i64pqSeeTopKey
#define gk_i64pqSeeTopVal                        MUMPS_gk_i64pqSeeTopVal
#define gk_i64pqUpdate                           MUMPS_gk_i64pqUpdate
#define gk_i64readfile                           MUMPS_gk_i64readfile
#define gk_i64readfilebin                        MUMPS_gk_i64readfilebin
#define gk_i64realloc                            MUMPS_gk_i64realloc
#define gk_i64scale                              MUMPS_gk_i64scale
#define gk_i64set                                MUMPS_gk_i64set
#define gk_i64SetMatrix                          MUMPS_gk_i64SetMatrix
#define gk_i64smalloc                            MUMPS_gk_i64smalloc
#define gk_i64sortd                              MUMPS_gk_i64sortd
#define gk_i64sorti                              MUMPS_gk_i64sorti
#define gk_i64sum                                MUMPS_gk_i64sum
#define gk_i64writefilebin                       MUMPS_gk_i64writefilebin
#define gk_i8AllocMatrix                         MUMPS_gk_i8AllocMatrix
#define gk_i8copy                                MUMPS_gk_i8copy
#define gk_i8FreeMatrix                          MUMPS_gk_i8FreeMatrix
#define gk_i8kvAllocMatrix                       MUMPS_gk_i8kvAllocMatrix
#define gk_i8kvcopy                              MUMPS_gk_i8kvcopy
#define gk_i8kvFreeMatrix                        MUMPS_gk_i8kvFreeMatrix
#define gk_i8kvmalloc                            MUMPS_gk_i8kvmalloc
#define gk_i8kvrealloc                           MUMPS_gk_i8kvrealloc
#define gk_i8kvset                               MUMPS_gk_i8kvset
#define gk_i8kvSetMatrix                         MUMPS_gk_i8kvSetMatrix
#define gk_i8kvsmalloc                           MUMPS_gk_i8kvsmalloc
#define gk_i8malloc                              MUMPS_gk_i8malloc
#define gk_i8realloc                             MUMPS_gk_i8realloc
#define gk_i8set                                 MUMPS_gk_i8set
#define gk_i8SetMatrix                           MUMPS_gk_i8SetMatrix
#define gk_i8smalloc                             MUMPS_gk_i8smalloc
#define gk_iAllocMatrix                          MUMPS_gk_iAllocMatrix
#define gk_iargmax                               MUMPS_gk_iargmax
#define gk_iargmax_n                             MUMPS_gk_iargmax_n
#define gk_iargmin                               MUMPS_gk_iargmin
#define gk_iaxpy                                 MUMPS_gk_iaxpy
#define gk_icopy                                 MUMPS_gk_icopy
#define gk_idot                                  MUMPS_gk_idot
#define gk_idxAllocMatrix                        MUMPS_gk_idxAllocMatrix
#define gk_idxargmax                             MUMPS_gk_idxargmax
#define gk_idxargmax_n                           MUMPS_gk_idxargmax_n
#define gk_idxargmin                             MUMPS_gk_idxargmin
#define gk_idxaxpy                               MUMPS_gk_idxaxpy
#define gk_idxcopy                               MUMPS_gk_idxcopy
#define gk_idxdot                                MUMPS_gk_idxdot
#define gk_idxFreeMatrix                         MUMPS_gk_idxFreeMatrix
#define gk_idxincset                             MUMPS_gk_idxincset
#define gk_idxkvAllocMatrix                      MUMPS_gk_idxkvAllocMatrix
#define gk_idxkvcopy                             MUMPS_gk_idxkvcopy
#define gk_idxkvFreeMatrix                       MUMPS_gk_idxkvFreeMatrix
#define gk_idxkvmalloc                           MUMPS_gk_idxkvmalloc
#define gk_idxkvrealloc                          MUMPS_gk_idxkvrealloc
#define gk_idxkvset                              MUMPS_gk_idxkvset
#define gk_idxkvSetMatrix                        MUMPS_gk_idxkvSetMatrix
#define gk_idxkvsmalloc                          MUMPS_gk_idxkvsmalloc
#define gk_idxkvsortd                            MUMPS_gk_idxkvsortd
#define gk_idxkvsorti                            MUMPS_gk_idxkvsorti
#define gk_idxmalloc                             MUMPS_gk_idxmalloc
#define gk_idxmax                                MUMPS_gk_idxmax
#define gk_idxmin                                MUMPS_gk_idxmin
#define gk_idxnorm2                              MUMPS_gk_idxnorm2
#define gk_idxpqCheckHeap                        MUMPS_gk_idxpqCheckHeap
#define gk_idxpqCreate                           MUMPS_gk_idxpqCreate
#define gk_idxpqDelete                           MUMPS_gk_idxpqDelete
#define gk_idxpqDestroy                          MUMPS_gk_idxpqDestroy
#define gk_idxpqFree                             MUMPS_gk_idxpqFree
#define gk_idxpqGetTop                           MUMPS_gk_idxpqGetTop
#define gk_idxpqInit                             MUMPS_gk_idxpqInit
#define gk_idxpqInsert                           MUMPS_gk_idxpqInsert
#define gk_idxpqLength                           MUMPS_gk_idxpqLength
#define gk_idxpqReset                            MUMPS_gk_idxpqReset
#define gk_idxpqSeeKey                           MUMPS_gk_idxpqSeeKey
#define gk_idxpqSeeTopKey                        MUMPS_gk_idxpqSeeTopKey
#define gk_idxpqSeeTopVal                        MUMPS_gk_idxpqSeeTopVal
#define gk_idxpqUpdate                           MUMPS_gk_idxpqUpdate
#define gk_idxrand                               MUMPS_gk_idxrand
#define gk_idxrandArrayPermute                   MUMPS_gk_idxrandArrayPermute
#define gk_idxrandArrayPermuteFine               MUMPS_gk_idxrandArrayPermuteFine
#define gk_idxrandInRange                        MUMPS_gk_idxrandInRange
#define gk_idxrealloc                            MUMPS_gk_idxrealloc
#define gk_idxscale                              MUMPS_gk_idxscale
#define gk_idxset                                MUMPS_gk_idxset
#define gk_idxSetMatrix                          MUMPS_gk_idxSetMatrix
#define gk_idxsmalloc                            MUMPS_gk_idxsmalloc
#define gk_idxsortd                              MUMPS_gk_idxsortd
#define gk_idxsorti                              MUMPS_gk_idxsorti
#define gk_idxsrand                              MUMPS_gk_idxsrand
#define gk_idxsum                                MUMPS_gk_idxsum
#define gk_ifkvkselect                           MUMPS_gk_ifkvkselect
#define gk_iFreeMatrix                           MUMPS_gk_iFreeMatrix
#define gk_iincset                               MUMPS_gk_iincset
#define gk_ikvAllocMatrix                        MUMPS_gk_ikvAllocMatrix
#define gk_ikvcopy                               MUMPS_gk_ikvcopy
#define gk_ikvFreeMatrix                         MUMPS_gk_ikvFreeMatrix
#define gk_ikvmalloc                             MUMPS_gk_ikvmalloc
#define gk_ikvrealloc                            MUMPS_gk_ikvrealloc
#define gk_ikvset                                MUMPS_gk_ikvset
#define gk_ikvSetMatrix                          MUMPS_gk_ikvSetMatrix
#define gk_ikvsmalloc                            MUMPS_gk_ikvsmalloc
#define gk_ikvsortd                              MUMPS_gk_ikvsortd
#define gk_ikvsorti                              MUMPS_gk_ikvsorti
#define gk_imalloc                               MUMPS_gk_imalloc
#define gk_imax                                  MUMPS_gk_imax
#define gk_imin                                  MUMPS_gk_imin
#define gk_inorm2                                MUMPS_gk_inorm2
#define gk_ipqCheckHeap                          MUMPS_gk_ipqCheckHeap
#define gk_ipqCreate                             MUMPS_gk_ipqCreate
#define gk_ipqDelete                             MUMPS_gk_ipqDelete
#define gk_ipqDestroy                            MUMPS_gk_ipqDestroy
#define gk_ipqFree                               MUMPS_gk_ipqFree
#define gk_ipqGetTop                             MUMPS_gk_ipqGetTop
#define gk_ipqInit                               MUMPS_gk_ipqInit
#define gk_ipqInsert                             MUMPS_gk_ipqInsert
#define gk_ipqLength                             MUMPS_gk_ipqLength
#define gk_ipqReset                              MUMPS_gk_ipqReset
#define gk_ipqSeeKey                             MUMPS_gk_ipqSeeKey
#define gk_ipqSeeTopKey                          MUMPS_gk_ipqSeeTopKey
#define gk_ipqSeeTopVal                          MUMPS_gk_ipqSeeTopVal
#define gk_ipqUpdate                             MUMPS_gk_ipqUpdate
#define gk_irand                                 MUMPS_gk_irand
#define gk_irandArrayPermute                     MUMPS_gk_irandArrayPermute
#define gk_irandArrayPermuteFine                 MUMPS_gk_irandArrayPermuteFine
#define gk_irandInRange                          MUMPS_gk_irandInRange
#define gk_irealloc                              MUMPS_gk_irealloc
#define gk_iscale                                MUMPS_gk_iscale
#define gk_iset                                  MUMPS_gk_iset
#define gk_iSetMatrix                            MUMPS_gk_iSetMatrix
#define gk_ismalloc                              MUMPS_gk_ismalloc
#define gk_isortd                                MUMPS_gk_isortd
#define gk_isorti                                MUMPS_gk_isorti
#define gk_ispow2                                MUMPS_gk_ispow2
#define gk_isrand                                MUMPS_gk_isrand
#define gk_isum                                  MUMPS_gk_isum
#define gk_jbuf                                  MUMPS_gk_jbuf
#define gk_jbufs                                 MUMPS_gk_jbufs
#define gk_log2                                  MUMPS_gk_log2
#define gk_malloc                                MUMPS_gk_malloc
#define gk_malloc_cleanup                        MUMPS_gk_malloc_cleanup
#define gk_malloc_init                           MUMPS_gk_malloc_init
#define gk_mcoreAdd                              MUMPS_gk_mcoreAdd
#define gk_mcoreCreate                           MUMPS_gk_mcoreCreate
#define gk_mcoreDel                              MUMPS_gk_mcoreDel
#define gk_mcoreDestroy                          MUMPS_gk_mcoreDestroy
#define gk_mcoreMalloc                           MUMPS_gk_mcoreMalloc
#define gk_mcorePop                              MUMPS_gk_mcorePop
#define gk_mcorePush                             MUMPS_gk_mcorePush
#define gk_mkpath                                MUMPS_gk_mkpath
#define gk_NonLocalExit_Handler                  MUMPS_gk_NonLocalExit_Handler
#define gk_optarg                                MUMPS_gk_optarg
#define gk_opterr                                MUMPS_gk_opterr
#define gk_optind                                MUMPS_gk_optind
#define gk_optopt                                MUMPS_gk_optopt
#define gk_randinit                              MUMPS_gk_randinit
#define gk_randint32                             MUMPS_gk_randint32
#define gk_randint64                             MUMPS_gk_randint64
#define gk_RandomPermute                         MUMPS_gk_RandomPermute
#define gk_read                                  MUMPS_gk_read
#define gk_readfile                              MUMPS_gk_readfile
#define gk_realloc                               MUMPS_gk_realloc
#define gk_rmpath                                MUMPS_gk_rmpath
#define gk_rw_PageRank                           MUMPS_gk_rw_PageRank
#define gk_seq_free                              MUMPS_gk_seq_free
#define gk_seq_init                              MUMPS_gk_seq_init
#define gk_seq_ReadGKMODPSSM                     MUMPS_gk_seq_ReadGKMODPSSM
#define gk_set_exit_on_error                     MUMPS_gk_set_exit_on_error
#define gk_SetSignalHandlers                     MUMPS_gk_SetSignalHandlers
#define gk_sigthrow                              MUMPS_gk_sigthrow
#define gk_sigtrap                               MUMPS_gk_sigtrap
#define gk_siguntrap                             MUMPS_gk_siguntrap
#define gk_skvAllocMatrix                        MUMPS_gk_skvAllocMatrix
#define gk_skvcopy                               MUMPS_gk_skvcopy
#define gk_skvFreeMatrix                         MUMPS_gk_skvFreeMatrix
#define gk_skvmalloc                             MUMPS_gk_skvmalloc
#define gk_skvrealloc                            MUMPS_gk_skvrealloc
#define gk_skvset                                MUMPS_gk_skvset
#define gk_skvSetMatrix                          MUMPS_gk_skvSetMatrix
#define gk_skvsmalloc                            MUMPS_gk_skvsmalloc
#define gk_skvsortd                              MUMPS_gk_skvsortd
#define gk_skvsorti                              MUMPS_gk_skvsorti
#define gk_str2time                              MUMPS_gk_str2time
#define gk_strcasecmp                            MUMPS_gk_strcasecmp
#define gk_strchr_replace                        MUMPS_gk_strchr_replace
#define gk_strdup                                MUMPS_gk_strdup
#define gk_strerror                              MUMPS_gk_strerror
#define gk_strhprune                             MUMPS_gk_strhprune
#define gk_strrcmp                               MUMPS_gk_strrcmp
#define gk_strstr_replace                        MUMPS_gk_strstr_replace
#define gk_strtokenize                           MUMPS_gk_strtokenize
#define gk_strtolower                            MUMPS_gk_strtolower
#define gk_strtoupper                            MUMPS_gk_strtoupper
#define gk_strtprune                             MUMPS_gk_strtprune
#define gk_time2str                              MUMPS_gk_time2str
#define gk_ui16AllocMatrix                       MUMPS_gk_ui16AllocMatrix
#define gk_ui16copy                              MUMPS_gk_ui16copy
#define gk_ui16FreeMatrix                        MUMPS_gk_ui16FreeMatrix
#define gk_ui16malloc                            MUMPS_gk_ui16malloc
#define gk_ui16realloc                           MUMPS_gk_ui16realloc
#define gk_ui16set                               MUMPS_gk_ui16set
#define gk_ui16SetMatrix                         MUMPS_gk_ui16SetMatrix
#define gk_ui16smalloc                           MUMPS_gk_ui16smalloc
#define gk_ui32AllocMatrix                       MUMPS_gk_ui32AllocMatrix
#define gk_ui32copy                              MUMPS_gk_ui32copy
#define gk_ui32FreeMatrix                        MUMPS_gk_ui32FreeMatrix
#define gk_ui32malloc                            MUMPS_gk_ui32malloc
#define gk_ui32realloc                           MUMPS_gk_ui32realloc
#define gk_ui32set                               MUMPS_gk_ui32set
#define gk_ui32SetMatrix                         MUMPS_gk_ui32SetMatrix
#define gk_ui32smalloc                           MUMPS_gk_ui32smalloc
#define gk_ui32sortd                             MUMPS_gk_ui32sortd
#define gk_ui32sorti                             MUMPS_gk_ui32sorti
#define gk_ui64AllocMatrix                       MUMPS_gk_ui64AllocMatrix
#define gk_ui64copy                              MUMPS_gk_ui64copy
#define gk_ui64FreeMatrix                        MUMPS_gk_ui64FreeMatrix
#define gk_ui64malloc                            MUMPS_gk_ui64malloc
#define gk_ui64realloc                           MUMPS_gk_ui64realloc
#define gk_ui64set                               MUMPS_gk_ui64set
#define gk_ui64SetMatrix                         MUMPS_gk_ui64SetMatrix
#define gk_ui64smalloc                           MUMPS_gk_ui64smalloc
#define gk_ui64sortd                             MUMPS_gk_ui64sortd
#define gk_ui64sorti                             MUMPS_gk_ui64sorti
#define gk_ui8AllocMatrix                        MUMPS_gk_ui8AllocMatrix
#define gk_ui8copy                               MUMPS_gk_ui8copy
#define gk_ui8FreeMatrix                         MUMPS_gk_ui8FreeMatrix
#define gk_ui8malloc                             MUMPS_gk_ui8malloc
#define gk_ui8realloc                            MUMPS_gk_ui8realloc
#define gk_ui8set                                MUMPS_gk_ui8set
#define gk_ui8SetMatrix                          MUMPS_gk_ui8SetMatrix
#define gk_ui8smalloc                            MUMPS_gk_ui8smalloc
#define gk_UnsetSignalHandlers                   MUMPS_gk_UnsetSignalHandlers
#define gk_WClockSeconds                         MUMPS_gk_WClockSeconds
#define gk_write                                 MUMPS_gk_write
#define gk_zAllocMatrix                          MUMPS_gk_zAllocMatrix
#define gk_zargmax                               MUMPS_gk_zargmax
#define gk_zargmax_n                             MUMPS_gk_zargmax_n
#define gk_zargmin                               MUMPS_gk_zargmin
#define gk_zaxpy                                 MUMPS_gk_zaxpy
#define gk_zcopy                                 MUMPS_gk_zcopy
#define gk_zdot                                  MUMPS_gk_zdot
#define gk_zFreeMatrix                           MUMPS_gk_zFreeMatrix
#define gk_zincset                               MUMPS_gk_zincset
#define gk_zkvAllocMatrix                        MUMPS_gk_zkvAllocMatrix
#define gk_zkvcopy                               MUMPS_gk_zkvcopy
#define gk_zkvFreeMatrix                         MUMPS_gk_zkvFreeMatrix
#define gk_zkvmalloc                             MUMPS_gk_zkvmalloc
#define gk_zkvrealloc                            MUMPS_gk_zkvrealloc
#define gk_zkvset                                MUMPS_gk_zkvset
#define gk_zkvSetMatrix                          MUMPS_gk_zkvSetMatrix
#define gk_zkvsmalloc                            MUMPS_gk_zkvsmalloc
#define gk_zkvsortd                              MUMPS_gk_zkvsortd
#define gk_zkvsorti                              MUMPS_gk_zkvsorti
#define gk_zmalloc                               MUMPS_gk_zmalloc
#define gk_zmax                                  MUMPS_gk_zmax
#define gk_zmin                                  MUMPS_gk_zmin
#define gk_znorm2                                MUMPS_gk_znorm2
#define gk_zrand                                 MUMPS_gk_zrand
#define gk_zrandArrayPermute                     MUMPS_gk_zrandArrayPermute
#define gk_zrandArrayPermuteFine                 MUMPS_gk_zrandArrayPermuteFine
#define gk_zrandInRange                          MUMPS_gk_zrandInRange
#define gk_zreadfile                             MUMPS_gk_zreadfile
#define gk_zreadfilebin                          MUMPS_gk_zreadfilebin
#define gk_zrealloc                              MUMPS_gk_zrealloc
#define gk_zscale                                MUMPS_gk_zscale
#define gk_zset                                  MUMPS_gk_zset
#define gk_zSetMatrix                            MUMPS_gk_zSetMatrix
#define gk_zsmalloc                              MUMPS_gk_zsmalloc
#define gk_zsrand                                MUMPS_gk_zsrand
#define gk_zsum                                  MUMPS_gk_zsum
#define gk_zuAllocMatrix                         MUMPS_gk_zuAllocMatrix
#define gk_zuargmax                              MUMPS_gk_zuargmax
#define gk_zuargmax_n                            MUMPS_gk_zuargmax_n
#define gk_zuargmin                              MUMPS_gk_zuargmin
#define gk_zuaxpy                                MUMPS_gk_zuaxpy
#define gk_zucopy                                MUMPS_gk_zucopy
#define gk_zudot                                 MUMPS_gk_zudot
#define gk_zuFreeMatrix                          MUMPS_gk_zuFreeMatrix
#define gk_zuincset                              MUMPS_gk_zuincset
#define gk_zukvAllocMatrix                       MUMPS_gk_zukvAllocMatrix
#define gk_zukvcopy                              MUMPS_gk_zukvcopy
#define gk_zukvFreeMatrix                        MUMPS_gk_zukvFreeMatrix
#define gk_zukvmalloc                            MUMPS_gk_zukvmalloc
#define gk_zukvrealloc                           MUMPS_gk_zukvrealloc
#define gk_zukvset                               MUMPS_gk_zukvset
#define gk_zukvSetMatrix                         MUMPS_gk_zukvSetMatrix
#define gk_zukvsmalloc                           MUMPS_gk_zukvsmalloc
#define gk_zukvsortd                             MUMPS_gk_zukvsortd
#define gk_zukvsorti                             MUMPS_gk_zukvsorti
#define gk_zumalloc                              MUMPS_gk_zumalloc
#define gk_zumax                                 MUMPS_gk_zumax
#define gk_zumin                                 MUMPS_gk_zumin
#define gk_zunorm2                               MUMPS_gk_zunorm2
#define gk_zurand                                MUMPS_gk_zurand
#define gk_zurandArrayPermute                    MUMPS_gk_zurandArrayPermute
#define gk_zurandArrayPermuteFine                MUMPS_gk_zurandArrayPermuteFine
#define gk_zurandInRange                         MUMPS_gk_zurandInRange
#define gk_zurealloc                             MUMPS_gk_zurealloc
#define gk_zuscale                               MUMPS_gk_zuscale
#define gk_zuset                                 MUMPS_gk_zuset
#define gk_zuSetMatrix                           MUMPS_gk_zuSetMatrix
#define gk_zusmalloc                             MUMPS_gk_zusmalloc
#define gk_zusrand                               MUMPS_gk_zusrand
#define gk_zusum                                 MUMPS_gk_zusum
#define gk_zwritefilebin                         MUMPS_gk_zwritefilebin
#define HTable_Create                            MUMPS_HTable_Create
#define HTable_Delete                            MUMPS_HTable_Delete
#define HTable_Destroy                           MUMPS_HTable_Destroy
#define HTable_GetNext                           MUMPS_HTable_GetNext
#define HTable_HFunction                         MUMPS_HTable_HFunction
#define HTable_Insert                            MUMPS_HTable_Insert
#define HTable_Reset                             MUMPS_HTable_Reset
#define HTable_Resize                            MUMPS_HTable_Resize
#define HTable_Search                            MUMPS_HTable_Search
#define HTable_SearchAndDelete                   MUMPS_HTable_SearchAndDelete
#define itemsets_find_frequent_itemsets          MUMPS_itemsets_find_frequent_itemsets
#define itemsets_project_matrix                  MUMPS_itemsets_project_matrix
#define PrintBackTrace                           MUMPS_PrintBackTrace
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

#include <gk_proto.h>


#endif  /* GKlib.h */


