//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

#include "PrismHybrid.h"
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <util.h>
#include <cudd.h>
#include <dd.h>
#include "hybrid.h"
#include "PrismHybridGlob.h"
#include "jnipointer.h"

#define MAX_LOG_STRING_LEN 1024
#define MAX_ERR_STRING_LEN 1024

//------------------------------------------------------------------------------
// hybrid engine global variables
//------------------------------------------------------------------------------

// cudd manager
DdManager *ddman;

// logs
// global refs to log classes
static jclass main_log_cls = NULL;
static jclass tech_log_cls = NULL;
// global refs to log objects
static jobject main_log_obj = NULL;
static jobject tech_log_obj = NULL;
// method ids for print method in logs
static jmethodID main_log_mid = NULL;
static jmethodID main_log_warn = NULL;
static jmethodID tech_log_mid = NULL;

// error message
char error_message[MAX_ERR_STRING_LEN];

// details from numerical computation which may be queried
double last_unif;

//------------------------------------------------------------------------------
// cudd manager
//------------------------------------------------------------------------------

JNIEXPORT void JNICALL Java_hybrid_PrismHybrid_PH_1SetCUDDManager(JNIEnv *env, jclass cls, jlong __jlongpointer ddm)
{
	ddman = jlong_to_DdManager(ddm);
}

//------------------------------------------------------------------------------
// logs
//------------------------------------------------------------------------------

JNIEXPORT void JNICALL Java_hybrid_PrismHybrid_PH_1SetMainLog(JNIEnv *env, jclass cls, jobject log)
{
	// if main log has been set previously, we need to delete existing global refs first
	if (main_log_obj != NULL) {
		env->DeleteGlobalRef(main_log_cls);
		main_log_cls = NULL;
		env->DeleteGlobalRef(main_log_obj);
		main_log_obj = NULL;
	}
	
	// make a global reference to the log object
	main_log_obj = env->NewGlobalRef(log);
	// get the log class and make a global reference to it
	main_log_cls = (jclass)env->NewGlobalRef(env->GetObjectClass(main_log_obj));
	// get the method id for the print method
	main_log_mid = env->GetMethodID(main_log_cls, "print", "(Ljava/lang/String;)V");
    main_log_warn = env->GetMethodID(main_log_cls, "printWarning", "(Ljava/lang/String;)V");
}

//------------------------------------------------------------------------------

JNIEXPORT void JNICALL Java_hybrid_PrismHybrid_PH_1SetTechLog(JNIEnv *env, jclass cls, jobject log)
{
	// if tech log has been set previously, we need to delete existing global refs first
	if (tech_log_obj != NULL) {
		env->DeleteGlobalRef(tech_log_cls);
		tech_log_cls = NULL;
		env->DeleteGlobalRef(tech_log_obj);
		tech_log_obj = NULL;
	}
	
	// make a global reference to the log object
	tech_log_obj = env->NewGlobalRef(log);
	// get the log class and make a global reference to it
	tech_log_cls = (jclass)env->NewGlobalRef(env->GetObjectClass(tech_log_obj));
	// get the method id for the print method
	tech_log_mid = env->GetMethodID(tech_log_cls, "print", "(Ljava/lang/String;)V");
}

//------------------------------------------------------------------------------

void PH_PrintToMainLog(JNIEnv *env, const char *str, ...)
{
	va_list argptr;
	char full_string[MAX_LOG_STRING_LEN];
	
	va_start(argptr, str);
	vsnprintf(full_string, MAX_LOG_STRING_LEN, str, argptr);
	va_end(argptr);
	
	if (env)
		env->CallVoidMethod(main_log_obj, main_log_mid, env->NewStringUTF(full_string));
	else
		printf("%s", full_string);
}

//------------------------------------------------------------------------------

void PH_PrintWarningToMainLog(JNIEnv *env, const char *str, ...)
{
	va_list argptr;
	char full_string[MAX_LOG_STRING_LEN];
	
	va_start(argptr, str);
	vsnprintf(full_string, MAX_LOG_STRING_LEN, str, argptr);
	va_end(argptr);
	
	if (env)
		env->CallVoidMethod(main_log_obj, main_log_warn, env->NewStringUTF(full_string));
	else
		printf("\nWarning: %s\n", full_string);
}

//------------------------------------------------------------------------------

void PH_PrintToTechLog(JNIEnv *env, const char *str, ...)
{
	va_list argptr;
	char full_string[MAX_LOG_STRING_LEN];
	
	va_start(argptr, str);
	vsnprintf(full_string, MAX_LOG_STRING_LEN, str, argptr);
	va_end(argptr);
	
	if (env)
		env->CallVoidMethod(tech_log_obj, tech_log_mid, env->NewStringUTF(full_string));
	else
		printf("%s", full_string);
}

//------------------------------------------------------------------------------

// Print formatted memory info to main log

void PH_PrintMemoryToMainLog(JNIEnv *env, const char *before, double mem, const char *after)
{
	char full_string[MAX_LOG_STRING_LEN];
	
	if (mem > 1048576)
		snprintf(full_string, MAX_LOG_STRING_LEN, "%s%.1f GB%s", before, mem/1048576.0, after);
	else if (mem > 1024)
		snprintf(full_string, MAX_LOG_STRING_LEN, "%s%.1f MB%s", before, mem/1024.0, after);
	else
		snprintf(full_string, MAX_LOG_STRING_LEN, "%s%.1f KB%s", before, mem, after);
	
	if (env) {
		env->CallVoidMethod(main_log_obj, main_log_mid, env->NewStringUTF(full_string));
	}
	else {
		printf("%s", full_string);
	}
}

//------------------------------------------------------------------------------
// error message handling
//------------------------------------------------------------------------------

void PH_SetErrorMessage(const char *str, ...)
{
	va_list argptr;
	
	va_start(argptr, str);
	vsnprintf(error_message, MAX_LOG_STRING_LEN, str, argptr);
	va_end(argptr);
}

char *PH_GetErrorMessage()
{
	return error_message;
}

JNIEXPORT jstring JNICALL Java_hybrid_PrismHybrid_PH_1GetErrorMessage(JNIEnv *env, jclass cls)
{
	return env->NewStringUTF(error_message);
}

//------------------------------------------------------------------------------
// numerical computation detail queries
//------------------------------------------------------------------------------

JNIEXPORT jdouble JNICALL Java_hybrid_PrismHybrid_PH_1GetLastUnif(JNIEnv *env, jclass cls)
{
	return last_unif;
}

//------------------------------------------------------------------------------
// tidy up
//------------------------------------------------------------------------------

JNIEXPORT void JNICALL Java_hybrid_PrismHybrid_PH_1FreeGlobalRefs(JNIEnv *env, jclass cls)
{
	// delete all global references
	env->DeleteGlobalRef(main_log_cls);
	main_log_cls = NULL;
	env->DeleteGlobalRef(tech_log_cls);
	tech_log_cls = NULL;
	env->DeleteGlobalRef(main_log_obj);
	main_log_obj = NULL;
	env->DeleteGlobalRef(tech_log_obj);
	tech_log_obj = NULL;
}

//------------------------------------------------------------------------------
