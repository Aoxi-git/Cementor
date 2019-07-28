// 2006-2008 © Václav Šmilauer <eudoxos@arcig.cz>
// 2019        Janek Kozicki
#pragma once
/*
 * This file defines various useful logging-related macros - userspace stuff is
 * - LOG_* for actual logging,
 * - DECLARE_LOGGER; that should be used in class header to create separate logger for that class,
 * - CREATE_LOGGER(className); that must be used in class implementation file to create the static variable.
 *
 * Note that the latter 2 may change their name to something like LOG_DECLARE and LOG_CREATE, to be consistent.
 * Some other macros will be very likely added, to allow for easy variable tracing etc. Suggestions welcome.
 *
 *
 * Yade has the logging config file by default in ~/.yade-$VERSION/logging.conf.
 *
 */


/* BOOST_LOGGER 014b11496
=======
#ifdef YADE_LOG4CXX

#	include<log4cxx/logger.h>
#	include<log4cxx/basicconfigurator.h>
#	include<log4cxx/propertyconfigurator.h>
#	include<log4cxx/helpers/exception.h>

#	define _LOG_HEAD __FILE__":"<<__LINE__<<" "<<__FUNCTION__<<": "
	// logger is local for every class, but if it is missing, we will use the parent's class logger automagically.
// TRACE doesn't really exist in log4cxx 0.9.7 (does in 0.10), otput through DEBUG then
#ifdef NDEBUG
#	define LOG_TRACE(msg){}
#	define LOG_DEBUG(msg){}
#else
#	ifdef LOG4CXX_TRACE
#		define LOG_TRACE(msg) {LOG4CXX_TRACE(logger, _LOG_HEAD<<msg);}
#	else
#		define LOG_TRACE(msg) {LOG4CXX_DEBUG(logger, _LOG_HEAD<<msg);}
#	endif
#	define LOG_DEBUG(msg) {LOG4CXX_DEBUG(logger, _LOG_HEAD<<msg);}
#endif
#	define LOG_INFO(msg)  {LOG4CXX_INFO(logger,  _LOG_HEAD<<msg);}
#	define LOG_WARN(msg)  {LOG4CXX_WARN(logger,  _LOG_HEAD<<msg);}
#	define LOG_ERROR(msg) {LOG4CXX_ERROR(logger, _LOG_HEAD<<msg);}
#	define LOG_FATAL(msg) {LOG4CXX_FATAL(logger, _LOG_HEAD<<msg);}

#	define DECLARE_LOGGER public: static log4cxx::LoggerPtr logger
#	define CREATE_LOGGER(classname) log4cxx::LoggerPtr classname::logger = log4cxx::Logger::getLogger("yade." #classname)

#else

#include<iostream>

>>>>>>> parent of 014b11496... Remove log4cxx support and yade.log module.

*/

#include <iostream>

#ifdef YADE_BOOST_LOG
#	include <boost/log/trivial.hpp>
#	define _POOR_MANS_LOG(level,msg) {BOOST_LOG_TRIVIAL(trace)<<" "<<_LOG_HEAD<<msg<<std::endl;}
#	define _LOG_HEAD __FILE__ ":"<<__LINE__<<" "<<__FUNCTION__<<": "
#else
#	define _POOR_MANS_LOG(level,msg) {std::cerr<<level " "<<_LOG_HEAD<<msg<<std::endl;}
#	define _LOG_HEAD __FILE__ ":"<<__LINE__<<" "<<__FUNCTION__<<": "
#endif

#if defined(YADE_DEBUG) ||  defined(YADE_TRACEONLY)
	# define LOG_TRACE(msg) _POOR_MANS_LOG("TRACE",msg)
#else
	# define LOG_TRACE(msg)
#endif

#if defined(YADE_DEBUG) && !defined(YADE_TRACEONLY)
	# define LOG_INFO(msg)  _POOR_MANS_LOG("INFO ",msg)
	# define LOG_DEBUG(msg) _POOR_MANS_LOG("DEBUG",msg)
#else
	# define LOG_INFO(msg)
	# define LOG_DEBUG(msg)
#endif


#define LOG_WARN(msg)  _POOR_MANS_LOG("WARN ",msg)
#define LOG_ERROR(msg) _POOR_MANS_LOG("ERROR",msg)
#define LOG_FATAL(msg) _POOR_MANS_LOG("FATAL",msg)

#define DECLARE_LOGGER
#define CREATE_LOGGER(classname)


// macros for quick debugging
#define TRACE LOG_TRACE("Been here")
#define _TRVHEAD cerr<<__FILE__<<":"<<__LINE__<<":"<<__FUNCTION__<<": "
#define _TRV(x) #x"="<<x<<"; "
#define TRVAR1(a) LOG_TRACE( _TRV(a) );
#define TRVAR2(a,b) LOG_TRACE( _TRV(a) << _TRV(b) );
#define TRVAR3(a,b,c) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) );
#define TRVAR4(a,b,c,d) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) << _TRV(d) );
#define TRVAR5(a,b,c,d,e) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) << _TRV(d) << _TRV(e) );
#define TRVAR6(a,b,c,d,e,f) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) << _TRV(d) << _TRV(e) << _TRV(f) );

