1. Building the IRTRuby Library
 1) Installs external classes/libraries
    * Boost Library : 
		- Download boost 1.45.0 library at http://www.boost.org/
		- Unzip the downloaded "boost_1_45_0.zip" 
					   into "{IRT_RUBY_HOME}/lib/include/bosst_1_45_0" directory

    * ETIRM Classes : Unzip "{IRT_RUBY_HOME}/lib/include/etirm/ETIRM.zip" 
					   into "{IRT_RUBY_HOME}/lib/include/etirm" directory
	    - Copy/Change files in "{IRT_RUBY_HOME}/lib/etrim/modify" directory 
		    			  into "{IRT_RUBY_HOME}/lib/include/etirm/src" directory
	  (Download : http://www.smallwaters.com/index.html)
    
	* SCPPNT Classes : Unzip "{IRT_RUBY_HOME}/lib/include/scppnt/SCPPNT.zip" 
						into "{IRT_RUBY_HOME}/lib/include/scppnt" directory
	  (Download : http://www.smallwaters.com/index.html)
	
	* UNCMIN Classes : Unzip "{IRT_RUBY_HOME}/lib/include/uncmin/uncmin.zip" 
						into "{IRT_RUBY_HOME}/lib/include/uncmin" directory
	  (Download : http://www.smallwaters.com/index.html)
	
	
 2) Generates "Makefile"
    - cd {IRT_RUBY_HOME}/lib
    - ruby extconf.rb --with-etirm-dir={IRT_RUBY_HOME}/lib/include/etirm/src \\
                      --with-scppnt-dir={IRT_RUBY_HOME}/lib/include/scppnt/src/include \\
					  --with-uncmin-dir={IRT_RUBY_HOME}/lib/include/uncmin/src/include \\
					  --with-boost-dir={IRT_RUBY_HOME}/lib/include/boost_1_45_0