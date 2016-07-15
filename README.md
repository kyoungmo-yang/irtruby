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

** Install
1. Building the IRTRuby Library
 1) Install SWIG (http://www.swig.org/)

    `sudo apt-get install swig`

 2) Installs external classes/libraries
    * Boost Library : 
		- Download boost 1.45.0 library at http://www.boost.org/
		  (Downlaod: https://sourceforge.net/projects/boost/files/boost/1.45.0/)
		- Unzip the downloaded "boost_1_45_0.zip" 
					   into "{IRT_RUBY_HOME}/lib/include/bosst_1_45_0" directory

    * ETIRM Classes : Unzip "{IRT_RUBY_HOME}/lib/include/etirm/ETIRM.zip" into "{IRT_RUBY_HOME}/lib/include/etirm" directory
	              Copy/Change files in "{IRT_RUBY_HOME}/lib/etrim/modify" directory into "{IRT_RUBY_HOME}/lib/include/etirm/src" directory
	  (Download : http://www.smallwaters.com/index.html)
    
	* SCPPNT Classes : Unzip "{IRT_RUBY_HOME}/lib/include/scppnt/SCPPNT.zip"  into "{IRT_RUBY_HOME}/lib/include/scppnt" directory (Download : http://www.smallwaters.com/index.html)
	
	* UNCMIN Classes : Unzip "{IRT_RUBY_HOME}/lib/include/uncmin/uncmin.zip" into "{IRT_RUBY_HOME}/lib/include/uncmin" directory (Download : http://www.smallwaters.com/index.html)
	
	
 3) Generates "Makefile"

    cd {IRT_RUBY_HOME}/lib
    ruby extconf.rb --with-etirm-dir={IRT_RUBY_HOME}/lib/include/etirm/src \\
                    --with-scppnt-dir={IRT_RUBY_HOME}/lib/include/scppnt/src/include \\
                    --with-uncmin-dir={IRT_RUBY_HOME}/lib/include/uncmin/src/include \\
                    --with-boost-dir={IRT_RUBY_HOME}/lib/include/boost_1_45_0
    
 4) Compile
    
    make
