// 2006-2008 © Václav Šmilauer
// 2019 Janek Kozicki
// hint: follow changes in d067b0696a8 to add new modules.

#include<core/Omega.hpp>
#include<lib/pyutil/doc_opts.hpp>
#include<lib/base/Logging.hpp>
#include<string>

CREATE_CPP_LOCAL_LOGGER("_log.cpp");

namespace py = boost::python;

void printNoBoostLogWarning() {
	std::cerr << "\nWarning: yade was compiled with cmake option -DBOOST_LOGGER=OFF, any attempts to manipulate log filter levels will not have effect.\n\n";
}

int getDefaultLogLevel() {
#ifdef YADE_BOOST_LOG
	return Logging::instance().getDefaultLogLevel();
#else
	return std::min(MAX_LOG_LEVEL,MAX_HARDCODED_LOG_LEVEL);
#endif
}

void setDefaultLogLevel(int level) {
#ifdef YADE_BOOST_LOG
	return Logging::instance().setDefaultLogLevel(level);
#else
	return std::min(MAX_LOG_LEVEL,MAX_HARDCODED_LOG_LEVEL);
#endif
}

void testAllLevels() {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
// ignore unused variable warning because when log level is low, they are not used (not printed).
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
	int testInt     = 0;
	std::string testStr = "test_string";
	Real testReal(11);
	Vector3r testVec(1,2,3);
	Matrix3r testMat = (Matrix3r() << 1, 2, 3, 4, 5, 6, 7, 8, 9).finished();
	std::complex<Real> testComplex(-1,1);

	LOG_0_NOFILTER("Current \"Default\" filter level is " << getDefaultLogLevel());

	LOG_6_TRACE   ("Test log level: LOG_6_TRACE   , test int: " << testInt++ << ", test string: "<< testStr);
	LOG_5_DEBUG   ("Test log level: LOG_5_DEBUG   , test int: " << testInt++ << ", test string: "<< testStr);
	LOG_4_INFO    ("Test log level: LOG_4_INFO    , test int: " << testInt++ << ", test string: "<< testStr);
	LOG_3_WARN    ("Test log level: LOG_3_WARN    , test int: " << testInt++ << ", test string: "<< testStr);
	LOG_2_ERROR   ("Test log level: LOG_2_ERROR   , test int: " << testInt++ << ", test string: "<< testStr);
	LOG_1_FATAL   ("Test log level: LOG_1_FATAL   , test int: " << testInt++ << ", test string: "<< testStr);
	LOG_0_NOFILTER("Test log level: LOG_0_NOFILTER, test int: " << testInt++ << ", test string: "<< testStr);

	LOG_0_NOFILTER("Below 6 variables are printed at filter level TRACE");
	TRVAR1(testInt);
	TRVAR2(testInt,testStr);
	TRVAR3(testInt,testStr,testReal);
	TRVAR4(testInt,testStr,testReal,testVec);
	TRVAR5(testInt,testStr,testReal,testVec,testMat);
	TRVAR6(testInt,testStr,testReal,testVec,testMat,testComplex);

	LOG_0_NOFILTER("Below is macro TRACE;");
	TRACE;
#pragma GCC diagnostic pop
}

// accepted streams: "clog", "cerr", "cout", "filename"
// It is possible to set different levels per log file, see notes about that in Logging::setOutputStream(…)
void setOutputStream(std::string streamName, bool reset /*, int level */ ) {
#ifdef YADE_BOOST_LOG
	Logging::instance().setOutputStream( streamName , reset );
	if(reset) {
		LOG_INFO("Log output stream has been set to "<< streamName <<". Other output streams were removed.");
	} else {
		LOG_INFO("Additional output stream has been set to "<< streamName <<".");
	}
#else
	printNoBoostLogWarning();
#endif
}

// resets all sinks to default values: 'cerr' for those below logLevel, 'none' for those above logLevel
void resetOutputStream() {
#ifdef YADE_BOOST_LOG
	Logging::instance().setOutputStream("clog" , true);
	LOG_INFO("Log output stream has been reset to std::clog. File sinks are not removed.");
#else
	printNoBoostLogWarning();
#endif
}

void setLevel(std::string className, int level) {
#ifdef YADE_BOOST_LOG
	Logging::instance().setNamedLogLevel(className , level);
	LOG_INFO("filter log level for " << className << " has been set to " << Logging::instance().getNamedLogLevel(className));
#else
	printNoBoostLogWarning();
#endif
}

void unsetLevel(std::string className) {
#ifdef YADE_BOOST_LOG
	Logging::instance().unsetNamedLogLevel(className);
	LOG_INFO("filter log level for " << className << " has been unset to " << Logging::instance().getNamedLogLevel(className));
#else
	printNoBoostLogWarning();
#endif
}

py::dict getAllLevels() {
	py::dict ret{};
#ifdef YADE_BOOST_LOG
	for(const auto& a : Logging::instance().getClassLogLevels()) {
		ret[a.first]=a.second;
	}
#else
	printNoBoostLogWarning();
#endif
	return ret;
}

py::dict getUsedLevels() {
	py::dict ret{};
#ifdef YADE_BOOST_LOG
	for(const auto& a : Logging::instance().getClassLogLevels()) {
		if(a.second != -1) {
			ret[a.first]=a.second;
		}
	}
#else
	printNoBoostLogWarning();
#endif
	return ret;
}

BOOST_PYTHON_MODULE(_log){
	YADE_SET_DOCSTRING_OPTS;
// We can use C++ string literal just like """ """ in python to write docstrings (see. https://en.cppreference.com/w/cpp/language/string_literal )
// The """ is a custom delimeter, we could use    R"RAW( instead, or any other delimeter. This decides what will be the termination delimeter.
// The docstrings can use syntax :param ……: ……… :return: ……. For details see https://thomas-cokelaer.info/tutorials/sphinx/docstring_python.html
	py::def("testAllLevels", testAllLevels, R"""(
This function prints test messages on all log levels. Can be used to see how filtering works and to what streams the logs are written.
	)""");
	py::def("getDefaultLogLevel", getDefaultLogLevel, R"""(
:return: The current ``Default`` filter log level.
	)""");
	py::def("setDefaultLogLevel", setDefaultLogLevel, R"""(
:param int level: Sets the ``Default`` filter log level, same as calling ``log.setLevel("Default",level)``.
	)""");
	py::def("setOutputStream", setOutputStream, R"""(
:param str streamName: sets the output stream, special names ``cout``, ``cerr``, and default ``clog`` use the ``std::cout``, ``std::cerr``, ``std::clog`` counterpart. Every other name means that log will be written to a filename.
:param bool reset: dictates whether all previously set output streams are to be removed. When set to false: the new output stream is set additionally to the current one.
	)""");
	py::def("resetOutputStream", resetOutputStream, R"""(
Resets log output stream to default state: all logs are printed on ``std::clog`` channel, which usually redirects to ``std::cerr``.
	)""");
	py::def("setLevel", setLevel , R"""(
:param str className: The logger name for which the filter level is to be set. Use name ``Default`` to change the default filter level.
:param int level: The filter level to be set.
.. warning:: setting ``Default`` log level higher than ``MAX_LOG_LEVEL`` provided during compilation will have no effect. Logs will not be printed because they are removed during compilation.
	)""");
	py::def("unsetLevel", unsetLevel , R"""(
:param str className: The logger name for which the filter level is to be unset, so that a ``Default`` will be used instead. Unsetting the ``Default`` level will change it to max level and print everything.
	)""");
	py::def("getAllLevels", getAllLevels , R"""(
:return: A python dictionary with all known loggers in yade. Those without a debug level set will have value -1 to indicate that ``Default`` filter log level is to be used for them.
	)""");
	py::def("getUsedLevels", getUsedLevels , R"""(
:return: A python dictionary with all used log levels in yade. Those without a debug level (value -1) are omitted.
	)""");

	py::scope().attr("TRACE")=int(6);
	py::scope().attr("DEBUG")=int(5);
	py::scope().attr("INFO") =int(4);
	py::scope().attr("WARN") =int(3);
	py::scope().attr("ERROR")=int(2);
	py::scope().attr("FATAL")=int(1);
	py::scope().attr("NOFILTER")=int(0);
}

/* this was in git revision 014b11496

#include<boost/python.hpp>
#include<string>
#include<lib/base/Logging.hpp>
#include<lib/pyutil/doc_opts.hpp>
using namespace boost;
enum{ll_TRACE,ll_DEBUG,ll_INFO,ll_WARN,ll_ERROR,ll_FATAL};

#ifdef YADE_LOG4CXX

	log4cxx::LoggerPtr logger=log4cxx::Logger::getLogger("yade.log");

	#include<log4cxx/logmanager.h>

	void logSetLevel(std::string loggerName,int level){
		std::string fullName(loggerName.empty()?"yade":("yade."+loggerName));
		if(!log4cxx::LogManager::exists(fullName)){
			LOG_WARN("No logger named "<<loggerName<<", ignoring level setting.");
			// throw std::invalid_argument("No logger named `"+fullName+"'");
		}
		log4cxx::LevelPtr l;
		switch(level){
			#ifdef LOG4CXX_TRACE
				case ll_TRACE: l=log4cxx::Level::getTrace(); break;
				case ll_DEBUG: l=log4cxx::Level::getDebug(); break;
				case ll_INFO:  l=log4cxx::Level::getInfo(); break;
				case ll_WARN:  l=log4cxx::Level::getWarn(); break;
				case ll_ERROR: l=log4cxx::Level::getError(); break;
				case ll_FATAL: l=log4cxx::Level::getFatal(); break;
			#else
				case ll_TRACE: l=log4cxx::Level::DEBUG; break;
				case ll_DEBUG: l=log4cxx::Level::DEBUG; break;
				case ll_INFO:  l=log4cxx::Level::INFO; break;
				case ll_WARN:  l=log4cxx::Level::WARN; break;
				case ll_ERROR: l=log4cxx::Level::ERROR; break;
				case ll_FATAL: l=log4cxx::Level::FATAL; break;
			#endif
			default: throw std::invalid_argument("Unrecognized logging level "+lexical_cast<std::string>(level));
		}
		log4cxx::LogManager::getLogger("yade."+loggerName)->setLevel(l);
	}
	void logLoadConfig(std::string f){ log4cxx::PropertyConfigurator::configure(f); }
#else
	bool warnedOnce=false;
	void logSetLevel(std::string loggerName, int level){
		// better somehow python's raise RuntimeWarning, but not sure how to do that from c++
		// it shouldn't be trapped by boost::python's exception translator, just print warning
		// Do like this for now.
		if(warnedOnce) return;
		LOG_WARN("Yade was compiled without boost::log support. Setting log levels from python will have no effect (warn once).");
		warnedOnce=true;
	}
	void logLoadConfig(std::string f){
		if(warnedOnce) return;
		LOG_WARN("Yade was compiled without boost::log support. Loading log file will have no effect (warn once).");
		warnedOnce=true;
	}
#endif

BOOST_PYTHON_MODULE(log){
	python::scope().attr("__doc__") = "Access and manipulation of log4cxx loggers.";

	YADE_SET_DOCSTRING_OPTS;

	python::def("setLevel",logSetLevel,(python::arg("logger"),python::arg("level")),"Set minimum severity *level* (constants ``TRACE``, ``DEBUG``, ``INFO``, ``WARN``, ``ERROR``, ``FATAL``) for given logger. \nLeading 'yade.' will be appended automatically to the logger name; if logger is '', the root logger 'yade' will be operated on.");
	python::def("loadConfig",logLoadConfig,(python::arg("fileName")),"Load configuration from file (log4cxx::PropertyConfigurator::configure)");
	python::scope().attr("TRACE")=(int)ll_TRACE;
	python::scope().attr("DEBUG")=(int)ll_DEBUG;
	python::scope().attr("INFO")= (int)ll_INFO;
	python::scope().attr("WARN")= (int)ll_WARN;
	python::scope().attr("ERROR")=(int)ll_ERROR;
	python::scope().attr("FATAL")=(int)ll_FATAL;
}

*/
