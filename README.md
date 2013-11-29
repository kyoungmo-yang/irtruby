Item Response Theory library for Ruby (IRTRuby)
http://code.google.com/p/irtruby/

IRTRuby is an ruby library to estimates items' parameters and examinees' ability 
based on Item Response Theory (IRT) ,and that is available over Unix(OSX)/Linux Systems.

IRTRuby consists of external c++ classes/libraries(boost1.45.0, ETIRM, SCPPNT, UNCMIN) 
and SWIG Ruby wrappers(irt_ruby.i, irt_ruby.cpp, irt_ruby.h).

The files and directories in this project are:

	* RELEASE - File containing release notes for IRTRuby.

	* LICENSE - File containing the license under which IRTRuby is distributed.

	* doc - This directory contains some documentation for IRTRuby.
	
	* test - This directory contains some unit-test/example programs for IRTRuby.
	
	* lib : 
		* lib/include - This directory contains external classes/libraries.
		* lib/extconf.rb - File containing script code to makes suitable "Makefile" 
						   for OSX or Linux Systems.
		* lib/include/README - File containing notes to makes "Makefile" 
							   for connecting C++ classes with IRT language.

