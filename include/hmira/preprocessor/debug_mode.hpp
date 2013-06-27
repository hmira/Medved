#ifndef __HMIRA_DEBUG_MODE__
#define __HMIRA_DEBUG_MODE__



	/**
	 * preprocessor requires flag 'HMIRA_DEBUG' and Debug mode enabled
	 * when the algorithms are running, it prints the debug messages on the err output
	 */
	#ifdef HMIRA_DEBUG

	#include <boost/preprocessor/seq/for_each.hpp>
	#include <boost/preprocessor/variadic/to_seq.hpp>

		#warning debug mode "HMIRA" is enabled

		#define H_DEBUG_STDERR(...) H_DEBUG_STD(std::cerr, __VA_ARGS__)
		#define H_DEBUG_STDOUT(...) H_DEBUG_STD(std::cout, __VA_ARGS__)

		#define H_DEBUG_STD(stream_name,...) \
			BOOST_PP_SEQ_FOR_EACH(H_DBG_PRINT, stream_name, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
			stream_name << std::endl;
			
		#define H_DBG_PRINT(a, stream_name, data)  stream_name << data;

	#else

		#define H_DEBUG_STDERR(...)
		#define H_DEBUG_STDOUT(...)
		#define H_DEBUG_STD(...)

	#endif //HMIRA_DEBUG


	/**
	 * preprocessor requires flag 'HMIRA_DEBUG_CUSTOM' and Debug mode enabled
	 * another level of debugging reserved for the needs of the user
	 * separatedly from the 'HMIRA_DEBUG' purposes
	 */
	#ifdef HMIRA_DEBUG_CUSTOM

	#include <boost/preprocessor/seq/for_each.hpp>
	#include <boost/preprocessor/variadic/to_seq.hpp>

		#warning debug mode "HMIRA" is enabled

		#define HC_DEBUG_STDERR(...) HC_DEBUG_STD(std::cerr, __VA_ARGS__)
		#define HC_DEBUG_STDOUT(...) HC_DEBUG_STD(std::cout, __VA_ARGS__)

		#define HC_DEBUG_STD(stream_name,...) \
			BOOST_PP_SEQ_FOR_EACH(HC_DBG_PRINT, stream_name, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
			stream_name << std::endl;
			
		#define HC_DBG_PRINT(a, stream_name, data)  stream_name << data;

	#else

		#define HC_DEBUG_STDERR(...)
		#define HC_DEBUG_STDOUT(...)
		#define HC_DEBUG_STD(...)

	#endif //HMIRA_DEBUG
		
#endif // __HMIRA_DEBUG_MODE__