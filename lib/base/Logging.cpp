/*************************************************************************
*  2019      © Janek Kozicki                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef YADE_BOOST_LOG

#include <ostream>
#include <fstream>
#include <lib/base/Logging.hpp>
#include <boost/core/null_deleter.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/phoenix/bind/bind_function.hpp>

SINGLETON_SELF(Logging);

bool logFilterLevels(  boost::log::value_ref< Logging::SeverityLevel , tag::severity      > const& level
		     , boost::log::value_ref< std::string            , tag::class_name_tag> const& name )
{
	auto itEnd          = Logging::instance().getClassLogLevels().end();
	short int itDefault = Logging::instance().getDefaultLogLevel();
	if(not level) { // something is seriously broken, we resort to std::cerr to report that.
		std::cerr << "LOGGER Warning: Logging::SeverityLevel is missing. Expect problems with logging.\n";
		return true;
	}
	if(name) {
		auto it = Logging::instance().getClassLogLevels().find(name.get());
		if( ( it != itEnd ) and ( it->second >= 0 ) ) {
			return level <= it->second;
		}
	}
	// this is triggered for LOG_NOFILTER macros. I comment this out, because they are legal now. Although I could 
	// slg.add_attribute("NameTag", boost::log::attributes::constant< std::string >("NoFilter"));
	// inside the LOG_NOFILTER macro. And then they would become filterable. Then they would need a different name. Because eNOFILTER is level 0 now.
	// I don't expect that a logger inside .hpp without a class logger available will be needed, so I don't do this now.
	//else {
	//	std::cerr << "LOGGER Warning: class_name_tag needed for filtering is missing. Expect problems with logging.\n";
	//}
	return level <= itDefault;
}

// Setup the common formatter for all sinks
Logging::Logging()
	: defaultLogLevel{(short int)(SeverityLevel::eWARN)}
	, classLogLevels{{"Default",defaultLogLevel}}
	, sink{boost::make_shared< TextSink >()}
	, streamClog(&std::clog, boost::null_deleter())
	, streamCerr(&std::cerr, boost::null_deleter())
	, streamCout(&std::cout, boost::null_deleter())
	, colors{true}
{
	boost::log::formatter fmt = boost::log::expressions::stream
		<< "<" << severity << "> "
		<< boost::log::expressions::if_(boost::log::expressions::has_attr(class_name_tag))
		[
			boost::log::expressions::stream << class_name_tag
		]
		<< boost::log::expressions::smessage;

	sink->locked_backend()->add_stream(streamClog);
	sink->set_formatter(fmt);
	sink->set_filter( boost::phoenix::bind(&logFilterLevels, severity.or_none(), class_name_tag.or_none() ));
	boost::log::core::get()->add_sink(sink);
}

// It is possible in boost::log to have different output sinks (e.g. different log files), each with a different filtering level.
// For now it is not necessary so I don't use it. If it becomes necessary in the future, see following examples:
//   https://www.boost.org/doc/libs/1_70_0/libs/log/doc/html/log/detailed/expressions.html#log.detailed.expressions.formatters.conditional
//   https://www.boost.org/doc/libs/1_70_0/libs/log/example/doc/tutorial_filtering.cpp
// To do this multiple variables like the variable 'sink' below would have to be defined. Each with a different filter level.
// Below all add_stream(…), remove_stream(…) calls are to the same sink. So they all have the same filtering level.
void Logging::setOutputStream(const std::string& name , bool reset) {
	if(reset) {
		sink->locked_backend()->remove_stream(streamClog);
		sink->locked_backend()->remove_stream(streamCerr);
		sink->locked_backend()->remove_stream(streamCout);
		sink->locked_backend()->remove_stream(streamFile);
		streamFile=boost::shared_ptr< std::ostream >{};
		for(const auto& oldFile : streamOld) {
			sink->locked_backend()->remove_stream(oldFile);
		}
		streamOld.clear();
	}
	switch(hash(name.c_str())) {
		case hash("clog") : sink->locked_backend()->add_stream(streamClog); break;
		case hash("cerr") : sink->locked_backend()->add_stream(streamCerr); break;
		case hash("cout") : sink->locked_backend()->add_stream(streamCout); break;
		default           : {
			if(not reset and streamFile) {
				std::cerr << "LOGGER Warning: adding a new log file without resetting the old one means that the logs will go to both files.\n";
				streamOld.push_back(streamFile);
			}
			streamFile = boost::make_shared< std::ofstream >(name.c_str());
			sink->locked_backend()->add_stream(streamFile);
			sink->locked_backend()->auto_flush(true); // make sure that log is written immediately to file without long time buffering.
		}
	}
}

void Logging::readConfigFile(const std::string& fname) {
	// TODO ...
}

void Logging::setDefaultLogLevel(short int level) {
	if(level < (short int)(SeverityLevel::eNOFILTER) or level > (short int)(SeverityLevel::eTRACE)) {
		throw std::runtime_error("The level must be >= NOFILTER (0) and <= TRACE (6), other values are not meaningful. To unset level to \"Default\" level use function unsetLevel(…).");
	}
	classLogLevels["Default"] = level;
	defaultLogLevel           = level;
}

short int Logging::getNamedLogLevel  (const std::string& name) {
	return findFilterName(name)->second;
}

void Logging::setNamedLogLevel  (const std::string& name , short int level) {
	if(level < (short int)(SeverityLevel::eNOFILTER) or level > (short int)(SeverityLevel::eTRACE)) {
		throw std::runtime_error("The level must be >= NOFILTER (0) and <= TRACE (6), other values are not meaningful. To unset level to \"Default\" level use function unsetLevel(…).");
	}
	if(level > maxLogLevel) {
		std::cerr << "LOGGER Warning: setting \""<<name<<"\" log level higher than MAX_LOG_LEVEL="<<maxLogLevel<<" will have no effect. Logs will not be printed, they were removed during compilation.\n";
		std::cerr << "LOGGER Warning: to be able to use \""<<name<<"\"="<<level<<" you have to recompile yade with cmake option MAX_LOG_LEVEL="<<level<<" or higher.\n";
	}
	if( name == "Default" ) {
		setDefaultLogLevel(level);
	} else {
		findFilterName(name)->second = level;
	}
}

void Logging::unsetNamedLogLevel(const std::string& name) {
	if(name == "Default") {
		// unsetting Default will result in printing everything.
		classLogLevels["Default"] = (short int)(SeverityLevel::eTRACE);
	} else {
		// unsetting anything else will result in printing it at Default level.
		findFilterName(name)->second = -1;
	}
}

boost::log::sources::severity_logger< Logging::SeverityLevel > Logging::createNamedLogger(std::string name) {
	boost::log::sources::severity_logger< Logging::SeverityLevel > l;
	l.add_attribute("NameTag", boost::log::attributes::constant< std::string >(name));
	classLogLevels[name] = -1;
	return l;
};

std::map<std::string,short int>::iterator Logging::findFilterName(const std::string& name) {
	auto it = classLogLevels.find(name);
	if( it == classLogLevels.end() ) {
		throw std::runtime_error(name+" is not recognized. Did you forget CREATE_LOGGER; and DECLARE_LOGGER(Classname); macros? Or maybe CREATE_CPP_LOCAL_LOGGER(\"filename.cpp\"); macro?\n");
	}
	return it;
}

#endif

