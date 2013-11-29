/*
	SWIG input to generate Ruby wrapper for Item Response Theory library for Ruby (IRTRuby)
	
	Item Response Theory library for Ruby (IRTRuby)
	http://code.google.com/p/irtruby/

	Author(s): 
	Embian, Inc., maintenance (Site : http://www.embian.com, Email : mo@embian.com, yjj0309@gmail.com)

	Copyright (C) 2010, Embian Inc.
*/

%module irt_ruby


%{
	#include "irt_ruby.h"
	using namespace etirm;
%}

/* Typemap for an argument that is a string giving an
   IRT model: 3PL, 2PL, 1Pl, GPCM, PCM.
   The string is mapped to an IRTModel enumeration */
%typemap(in) IRTModel irtModel
{
	//char *cstr = RSTRING($input) -> ptr;
	char *cstr = RSTRING_PTR($input);
	
	std::string modelStr(cstr);
	if (modelStr == "3PL") $1 = ThreePL;
	else if (modelStr == "GPCM") $1 = GPCM;
	else if (modelStr == "2PL") $1 = TwoPL;
	else if (modelStr == "1PL") $1 = OnePL;
	else if (modelStr == "PCM") $1 = PCM;
	else {rb_raise(rb_eArgError, "invalid IRT Model(model=%s), that should be one of the dichotomous model names(3PL, 2PL, 1PL) or the polytomous model names(PCM, GPCM)", cstr); SWIG_fail;}
}

/////////////////typemap : array - int_vector/////////////////////

%typemap(typecheck) int_vector *
{
	$1 = TYPE($input) == T_ARRAY ? 1 : 0;
}

%typemap(in) int_vector *
{
	Check_Type($input, T_ARRAY);
	//int len = RARRAY($input)->len;
	int len = RARRAY_LEN($input);
	int_vector *target = new int_vector(len);
	int_vector::iterator iv = target->begin();
	for (int i=0; i<len; ++i, ++iv){
		VALUE inst = rb_ary_entry($input, i);
		if(TYPE(inst) != T_FIXNUM && TYPE(inst) != T_BIGNUM){
			//rb_warn("wrong argument type's input value  (expected Fixnum or Bignum), so the input value will be cast integer type");
		}
		*iv = NUM2INT(inst);
	}
	$1 = target;
}

%typemap(freearg) int_vector *
{
 	if($1) delete $1;
}

%typemap(out) int_vector *
{
	VALUE arr = rb_ary_new2($1->size());
	int_vector::iterator iv = $1->begin(), iend = $1->end();
	for ( ; iv < iend; ++iv ) {		
		rb_ary_push(arr, INT2NUM(*iv));
	}
	delete $1;
	$result = arr;
}

/////////////////typemap : array - double_vector/////////////////////

%typemap(typecheck) double_vector *
{
	$1 = TYPE($input) == T_ARRAY ? 1 : 0;
}

%typemap(in) double_vector *
{
	Check_Type($input, T_ARRAY);
	//int len = RARRAY($input)->len;
	int len = RARRAY_LEN($input);	
	double_vector *target = new double_vector(len);
	double_vector::iterator iv = target->begin();
	for (int i=0; i<len; ++i, ++iv){
		VALUE inst = rb_ary_entry($input, i);
		if(TYPE(inst) != T_FLOAT){
			//rb_warn("wrong argument type's input value  (expected Float), so the input value will be cast float type");
		}
		*iv = NUM2DBL(inst);
	}
	$1 = target;
}

%typemap(freearg) double_vector *
{
 	if($1) delete $1;
}

%typemap(out) double_vector *
{
	VALUE arr = rb_ary_new2($1->size());
	double_vector::iterator iv = $1->begin(), iend = $1->end();
	for ( ; iv < iend; ++iv ) rb_ary_push(arr, rb_float_new(*iv));
	delete $1;
	$result = arr;
}

/////////////////typemap : string - char */////////////////////

%typemap(typecheck) char *
{
	$1 = TYPE($input) == T_STRING ? 1 : 0;
}

%typemap(in) char *
{
	Check_Type($input, T_STRING);
	//$1 = STR2CSTR($input);
	$1 = StringValuePtr($input);	
}

%typemap(freearg) char *
{
 	//if($1) free($1);
}

%typemap(out) char *
{
	VALUE str = rb_str_new2($1);
	delete $1;
	$result = str;
}


/////////////////exception/////////////////////
%exception {

	try
	{
		$action
	}
	catch (SCPPNT::Exception &e)
	{
		//rb_raise(rb_eRuntimeError, (char *) e.what()); 
		rb_raise(rb_eRuntimeError, "SCPPNT exception"); 
		SWIG_fail;
	}
	catch (std::exception &e)
	{
		//rb_raise(rb_eRuntimeError, (char *) e.what()); 
		rb_raise(rb_eRuntimeError, "Undefined exception"); 
		SWIG_fail;
	}
}


/* Declarations of commands to be wrapped by SWIG */
%include "swig_etirm.h"
%include "irt_ruby.h"