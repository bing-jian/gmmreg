/*=========================================================================
$Author: bingjian $
$Date: 2013-01-04 01:39:25 -0500 (Fri, 04 Jan 2013) $
$Revision: 145 $
=========================================================================*/

/**
 * \file gmmreg_api.h
 * \brief  The interface of calling gmmreg_api
 */

#ifndef GMMREG_API_H_
#define GMMREG_API_H_

#ifdef __cplusplus
extern "C"
#endif
int gmmreg_api(const char* f_config, const char* method);

#ifdef __cplusplus
extern "C"
#endif
void print_usage();

#endif  // GMMREG_API_H_


