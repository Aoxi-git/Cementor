// 2006-2008 © Václav Šmilauer <eudoxos@arcig.cz>
// 2019        Janek Kozicki
#pragma once
/*
 * This file defines various useful logging-related macros - userspace stuff is
 * - LOG_* for actual logging,
 * - DECLARE_LOGGER; that should be used in class header to create separate logger for that class,
 * - CREATE_LOGGER(className); that must be used in class implementation file to create the static variable.
 * - CREATE_LOCAL_LOGGER(name); use this inside a *.cpp file which has code that does not belong to any class, and needs logging. The name will be used for filtering logging.
 *
 * Note that the latter 2 may change their name to something like LOG_DECLARE and LOG_CREATE, to be consistent.
 * Some other macros will be very likely added, to allow for easy variable tracing etc. Suggestions welcome.
 *
 *
 * Yade has the logging config file by default in ~/.yade-$VERSION/logging.conf.
 *
 */

// boost::log inspired by git 014b11496
#ifdef YADE_BOOST_LOG
	#include <lib/base/Singleton.hpp>
	#include <boost/log/expressions.hpp>
	#include <boost/log/trivial.hpp>
	#include <boost/log/utility/setup.hpp>

	#define _LOG_HEAD ":"<<__LINE__<<" "<<__PRETTY_FUNCTION__<<": "
	// If you get "error: ‘logger’ was not declared in this scope" then you have to declare logger.
	// Use DECLARE_LOGGER; inside class and CREATE_LOGGER(className); inside .cpp file
	// or use CREATE_LOCAL_LOGGER("SomeName") if you need logging outside some class.
	#define LOG_TRACE(msg)  { BOOST_LOG_SEV(logger, severity_level::eTRACE) << _LOG_HEAD << msg; }
	#define LOG_MORE(msg)   { BOOST_LOG_SEV(logger, severity_level::eMORE)  << _LOG_HEAD << msg; }
	#define LOG_DEBUG(msg)  { BOOST_LOG_SEV(logger, severity_level::eDEBUG) << _LOG_HEAD << msg; }
	#define LOG_INFO(msg)   { BOOST_LOG_SEV(logger, severity_level::eINFO)  << _LOG_HEAD << msg; }
	#define LOG_WARN(msg)   { BOOST_LOG_SEV(logger, severity_level::eWARN)  << _LOG_HEAD << msg; }
	#define LOG_ERROR(msg)  { BOOST_LOG_SEV(logger, severity_level::eERROR) << _LOG_HEAD << msg; }
	#define LOG_THROW(msg)  { BOOST_LOG_SEV(logger, severity_level::eTHROW) << _LOG_HEAD << msg; }
	#define LOG_FATAL(msg)  { BOOST_LOG_SEV(logger, severity_level::eFATAL) << _LOG_HEAD << msg; }

	enum severity_level { eFATAL=1, eTHROW=2, eERROR=3, eWARN=4, eINFO=5, eDEBUG=6, eMORE=7, eTRACE=8 };
	inline std::ostream& operator<< (std::ostream& strm, severity_level level) // necessary for formatting output.
	{
		static const char* strings[] = { "UNKNOWN" , "FATAL_1" , "THROW_2" , "ERROR_3" , "WARN__4" , "INFO__5" , "DEBUG_6" , "MORE__7" , "TRACE_8" };
		if (static_cast< std::size_t >(level) < sizeof(strings) / sizeof(*strings))
			strm << strings[level];
		else
			strm << static_cast< int >(level);
		return strm;
	}

	BOOST_LOG_ATTRIBUTE_KEYWORD(severity      , "Severity"    , severity_level )
	BOOST_LOG_ATTRIBUTE_KEYWORD(class_name_tag, "ClassNameTag", std::string    )

	struct Logging : public Singleton<Logging> {
		std::map<std::string,int> classLogLevels;
		FRIEND_SINGLETON(Logging);
	};

	inline boost::log::sources::severity_logger< severity_level > createNamedLogger(std::string name) {
		boost::log::sources::severity_logger< severity_level > l;
		l.add_attribute("ClassNameTag", boost::log::attributes::constant< std::string >(name));
		Logging::instance().classLogLevels[name] = -1;
		return l;
	};
	// logger is local for every class, but if it is missing, we will use the parent's class logger automagically.
	#define DECLARE_LOGGER public: static boost::log::sources::severity_logger< severity_level > logger;
	#define CREATE_LOGGER(classname) boost::log::sources::severity_logger< severity_level > classname::logger=createNamedLogger(#classname);
	#define CREATE_LOCAL_LOGGER(filtername) namespace{ boost::log::sources::severity_logger< severity_level > logger=createNamedLogger(filtername); };
#else
	#include <iostream>
	#define _POOR_MANS_LOG(level,msg) {std::cerr<<level " "<<_LOG_HEAD<<msg<<std::endl;}
	#define _LOG_HEAD __FILE__ ":"<<__LINE__<<" "<<__PRETTY_FUNCTION__<<": "

	#ifdef YADE_DEBUG
	// when compiling with debug symbols and without boost::log it will print everything.
		#define LOG_TRACE(msg) _POOR_MANS_LOG("TRACE_8",msg)
		#define LOG_MORE(msg)  _POOR_MANS_LOG("MORE__7",msg)
		#define LOG_DEBUG(msg) _POOR_MANS_LOG("DEBUG_6",msg)
		#define LOG_INFO(msg)  _POOR_MANS_LOG("INFO__5",msg)
	#else
		#define LOG_TRACE(msg)
		#define LOG_MORE(msg)
		#define LOG_INFO(msg)
		#define LOG_DEBUG(msg)
	#endif

	#define LOG_WARN(msg)  _POOR_MANS_LOG("WARN__4",msg)
	#define LOG_ERROR(msg) _POOR_MANS_LOG("ERROR_3",msg)
	#define LOG_THROW(msg) _POOR_MANS_LOG("THROW_2",msg)
	#define LOG_FATAL(msg) _POOR_MANS_LOG("FATAL_1",msg)

	#define DECLARE_LOGGER
	#define CREATE_LOGGER(classname)
	#define CREATE_LOCAL_LOGGER(name)
#endif

// macros for quick debugging
#define TRACE LOG_TRACE("Been here")
#define _TRV(x) #x"="<<x<<"; "
#define TRVAR1(a) LOG_TRACE( _TRV(a) );
#define TRVAR2(a,b) LOG_TRACE( _TRV(a) << _TRV(b) );
#define TRVAR3(a,b,c) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) );
#define TRVAR4(a,b,c,d) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) << _TRV(d) );
#define TRVAR5(a,b,c,d,e) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) << _TRV(d) << _TRV(e) );
#define TRVAR6(a,b,c,d,e,f) LOG_TRACE( _TRV(a) << _TRV(b) << _TRV(c) << _TRV(d) << _TRV(e) << _TRV(f) );

// Logger aliases:
#define LOG_8_TRACE(msg) LOG_TRACE(msg)
#define LOG_7_MORE(msg)  LOG_MORE(msg)
#define LOG_6_DEBUG(msg) LOG_DEBUG(msg)
#define LOG_5_INFO(msg)  LOG_INFO(msg)
#define LOG_4_WARN(msg)  LOG_WARN(msg)
#define LOG_3_ERROR(msg) LOG_ERROR(msg)
#define LOG_2_THROW(msg) LOG_THROW(msg)
#define LOG_1_FATAL(msg) LOG_FATAL(msg)

#define LOG_8(msg) LOG_TRACE(msg)
#define LOG_7(msg) LOG_MORE(msg)
#define LOG_6(msg) LOG_DEBUG(msg)
#define LOG_5(msg) LOG_INFO(msg)
#define LOG_4(msg) LOG_WARN(msg)
#define LOG_3(msg) LOG_ERROR(msg)
#define LOG_2(msg) LOG_THROW(msg)
#define LOG_1(msg) LOG_FATAL(msg)

