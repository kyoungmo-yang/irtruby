/*
	Functions and classes for use with SWIG to generate scripting
	Item Response Theory library for Ruby (IRTRuby)
	
	Item Response Theory library for Ruby (IRTRuby)
	http://code.google.com/p/irtruby/

	Author(s): 
	Embian, Inc., maintenance (Site : www.embian.com, Email : mo@embian.com, yjj0309@gmail.com)
	Brad Hanson (http://www.b-a-h.com/)
	See the file LICENSE for information on usage and redistribution.

	Copyright (c) 2000-2002, Bradley A. Hanson
	Copyright (C) 2010, Embian, Inc.
	
*/


#include "irt_ruby.h"

#ifdef ETIRM_NO_DIR_PREFIX
#include "DiscreteNormalDist.h"
#include "Start3PL.h"
#include "ICCLogistic.h"
#include "ICRF_GPCM.h"
#include "ICRF_PCM.h"
#include "ItemDichotomous.h"
#include "ItemPolytomous.h"
#include "ItemParamPriorBeta4.h"
#include "ItemParamPriorLogNormal.h"
#include "ItemParamPriorNormal.h"

#include "MStepIRT.h"
#include "BootstrapSample.h"
#include "SimulateResponses.h"
#include "ExamineeThetaMLE.h"

#else
#include "etirm/DiscreteNormalDist.h"
#include "etirm/Start3PL.h"
#include "etirm/ICCLogistic.h"
#include "etirm/ICRF_GPCM.h"
#include "etirm/ICRF_PCM.h"
#include "etirm/ItemDichotomous.h"
#include "etirm/ItemPolytomous.h"
#include "etirm/ItemParamPriorBeta4.h"
#include "etirm/ItemParamPriorLogNormal.h"
#include "etirm/ItemParamPriorNormal.h"

#include "etirm/MStepIRT.h"
#include "etirm/BootstrapSample.h"
#include "etirm/SimulateResponses.h"
#include "etirm/ExamineeThetaMLE.h"

#endif

#include <boost/random/uniform_int.hpp>
#include <cstdio>
#include <string>

// Use normal random number generator from Boost library
// (www.boost.org).
// #include <boost/random/normal_distribution.hpp>


namespace etirm {

	// Prototype of local function
	int ParamIndex(const char *paramname, const char *funcname);

	/* Definitions of static data members */
	template <class L>
	Response Item<L>::notPresentedResponse = '.';
	template <class L>
	Response Item<L>::defaultFirstResponse = '0';

	// Item types for 1, 2, and 3 parameter logistic models
	typedef ItemDichotomous<DiscreteLatentDist<double>, ICCLogistic<3> > item3PL_type;
	typedef ItemDichotomous<DiscreteLatentDist<double>, ICCLogistic<2> > item2PL_type;
	typedef ItemDichotomous<DiscreteLatentDist<double>, ICCLogistic<1> > item1PL_type;

	// Item types for polytomous items
	typedef ItemPolytomous<DiscreteLatentDist<double>, ICRF_GPCM > itemGPCM_type;
	typedef ItemPolytomous<DiscreteLatentDist<double>, ICRF_PCM > itemPCM_type;
	
	//############################################################################################################################################
	// global variable holding information about a run defined in
	// swig_etirm.cpp
	// extern SwigEtirmRun *gEtirmRun;

	// global variable holding default values used in creating item objects
	// DefaultValues gIclDefaults;
	//############################################################################################################################################
	
	
	//############################################################################################################################################
	/*	
		Return a pointer to a item object.
	
		model		Model to use for item (ThreePL, TwoPL, OnePL, GPCM, PCM)
		index		Index of item in examinee response string (zero-offset)
		nrespcat	Number of response categories for polytomous items (ignored for dichtomous items)
		allPriors	Priors for a, b, and c parameters
		D			Logistic scaling constant (D=1.7 make logistic curve close to normal)
		latentDist	Latent variable distribution (discrete points used in SetLatentVarPoints call)
		funcname	Name of calling function (used in error messages)
	*/
	static item_type *CreateItem(IRTModel model, int index, int nrespcat, PriorVector &allPriors, double D, DiscreteLatentDist<Real> &latentDist, const char *funcname)
	{

		item_type *item;
	
		// Vector to hold default parameters assigned to item
		item_type::param_vector default_param(3);
	
		// Vector to hold prior assigned to item parameters
		PriorVector priors(3);
		switch (model)
		{
		case ThreePL:
		{
			ICCLogistic<3> icc(D);
			item = new item3PL_type(index, icc, &latentDist);
		
			// Set default parameters
			default_param[0] = 1.0; // default a parameter
			default_param[1] = 0.0; // default b parameter
			default_param[2] = 0.2; // default c parameter
			item->SetParameters(default_param);

			// Set default parameter priors
			for (int i=0; i<3; ++i)
			{
				if (allPriors[i])
				{
					priors[i] = CreatePrior(allPriors[i]->DistributionName(), 
						allPriors[i]->GetParameters(), funcname);
				}
				else
				{
					priors[i] = 0;
				}
			}
			item->SetPriors(priors);
		}
			break;
		
		case TwoPL:
		{
			ICCLogistic<2> icc(D, 0.0);
			item = new item2PL_type(index, icc, &latentDist);
		
			// Set default parameters
			default_param.newsize(2);
			default_param[0] = 1.0; // default a parameter
			default_param[1] = 0.0; // default b parameter
			item->SetParameters(default_param);

			// Set default priors for a and b parameters
			priors.resize(2);
			for (int i=0; i<2; ++i)
			{
				if (allPriors[i])
				{
					priors[i] = CreatePrior(allPriors[i]->DistributionName(), 
						allPriors[i]->GetParameters(), funcname);
				}
				else
				{
					priors[i] = 0;
				}
			}
			item->SetPriors(priors);
		}
			break;
		
		case OnePL:
		{
			ICCLogistic<1> icc(D, 1.0, 0.0);
			item = new item1PL_type(index, icc, &latentDist);
		
			// Set default parameter
			default_param.newsize(1);
			default_param[0] = 0.0; // default b parameter

			// Set default prior for b parameter
			priors.resize(1);
			if (allPriors[1])
			{
				priors[0] = CreatePrior(allPriors[1]->DistributionName(), 
					allPriors[1]->GetParameters(), funcname);
			}
			else
			{
				priors[0] = 0;
			}

			item->SetPriors(priors);
			item->SetParameters(default_param);
		}
			break;
		
		case GPCM:
		{
			ICRF_GPCM icrf(nrespcat, item_type::DefaultFirstResponse());
			item = new itemGPCM_type(index, icrf, &latentDist);
		
			default_param.newsize(nrespcat);
			priors.resize(nrespcat);
		
			// Set default parameter and prior for a
			default_param[0] = 1.0;
			if (allPriors[0])
			{
				priors[0] = CreatePrior(allPriors[0]->DistributionName(), 
						allPriors[0]->GetParameters(), funcname);
			}
			else priors[0] = 0;
		
			// Set default parameters and priors for b
			for (int i=1; i<nrespcat; ++i)
			{
				default_param[i] = 0.0;
				if (allPriors[1])
				{
					priors[i] = CreatePrior(allPriors[1]->DistributionName(), 
						allPriors[1]->GetParameters(), funcname);
				}
				else priors[i] = 0;
			}

			item->SetPriors(priors);
			item->SetParameters(default_param);
		}
			break;

		case PCM:
		{
			ICRF_PCM icrf(nrespcat, item_type::DefaultFirstResponse(), 1.0);
			item = new itemPCM_type(index, icrf, &latentDist);
		
			default_param.newsize(nrespcat-1);
			priors.resize(nrespcat-1);
		
			// Set default parameters and priors for b
			for (int i=0; i<nrespcat-1; ++i)
			{
				default_param[i] = 0.0;
				if (allPriors[1])
				{
					priors[i] = CreatePrior(allPriors[1]->DistributionName(), 
						allPriors[1]->GetParameters(), funcname);
				}
				else priors[i] = 0;
			}

			item->SetPriors(priors);
			item->SetParameters(default_param);
		}
			break;
		default:
			throw RuntimeError("Invalid model for item", funcname);
		}
	
		return item;

	}
	//############################################################################################################################################
	
	// Set initial default values
	DefaultValues::DefaultValues()
		: D(1.7), priors(3), dichModel(ThreePL), polyModel(GPCM), baseGroup(1), base_rand_normal(0), rand_normal(0)
	{

		priors[0] = new ItemParamPriorBeta4(1.75, 3.0, 0.0, 3.0);
		priors[1]= new ItemParamPriorBeta4(1.01, 1.01, -6.0, 6.0);
		priors[2] = new ItemParamPriorBeta4(3.5, 4.0, 0.0, 0.5);

	}

	// Release memory allocated by constructor
	DefaultValues::~DefaultValues()
	{
	
		PriorVector::iterator ip = priors.begin();
		for (int n = priors.size(); n--; ++ip)
		{
			if (*ip) delete *ip;
		}
	
		if (base_rand_normal) delete base_rand_normal;
		if (rand_normal) delete rand_normal;
	}
	
	//############################################################################################################################################
	
	//############################################################################################################################################
	IRTModule::IRTModule()
	{
		gIclDefaults = new DefaultValues();
		gEtirmRun = 0;
	}

	// Release memory allocated by constructor
	IRTModule::~IRTModule()
	{
	
		if (gIclDefaults) delete gIclDefaults;
		if (gEtirmRun) delete gEtirmRun;
	}
	//############################################################################################################################################
	
	//############################################################################################################################################
	// SwigIclRun constructor
	// nitems		Total number of items
	// ncat			Vector giving number of response categories for each item (1 = dichotomous item,
	//				>1 = polytomous item).
	// nlatentcat	Number of categories of discrete latent variable
	// ngroups		Number of examinee groups
	// minTheta		Minimum point of discrete latent variable distribution
	// maxTheta		Maximum point of discrete latent variable distribution
	// uniquePoints If true use different latent distribution points for different examinee groups
	IRTModule::SwigIclRun::SwigIclRun(int nitems, int_vector *ncat, int nlatentcat, int ngroups, double minTheta, double maxTheta, bool uniquePoints, DefaultValues *gIclDefaults)
			: SwigEtirmRun(nitems, nlatentcat, ngroups, minTheta, maxTheta, uniquePoints)
	{

		int i;

		// Create items
		for (i=0; i<numItems; i++) 
		{
			IRTModel model = gIclDefaults->dichModel;
			int n = 1;
			if (ncat)
			{
				n = (*ncat)[i];
				if (n > 1) model = gIclDefaults->polyModel;
			}
			items[i] = CreateItem(model, i, n, gIclDefaults->priors , gIclDefaults->D, latentDist, 
				"SwigIclRun::SwigIclRun");

			/* Set up Uncmin objects to use for optimization in M-step */
			minProc[i] = new Uncmin<RealVector, RealMatrix, item_type>(items[i]->NumParameters());
		
			itemStats[i] = new ItemRespCounts(items[i]->NumRespCat(),ngroups);
		}

	}

	// Release memory allocated by constructor
	IRTModule::SwigIclRun::~SwigIclRun()
	{
	}
	//############################################################################################################################################

	/* translate string containing item parameter name to zero-offset index.
		Prototype is in swig_etirm.h */
	int ParamIndex(const char *paramname, const char *funcname)
	{
		switch (*paramname)
		{
			case 'a' : return 0;
			case 'b' : return 1;
			case 'c' : return 2;
		}
	
		char errstr[100];
		std::sprintf(errstr, "Invalid item parameter name: %.6s", paramname);
		throw RuntimeError(errstr, funcname);

		return -1;
	}
	

	/** Functions wrapped by SWIG **/
	
	// Check that new_etirm has been called
	void IRTModule::CheckRunInit()
	{
		if (!gEtirmRun) throw RuntimeError("The new_items_dist command has not been executed", 0);

	}
	
	SwigEtirmRun * IRTModule::GetEtirmObj()
	{
		CheckRunInit();
		return gEtirmRun;
	}
	
	// Initialize memory to store information for items and the latent variable
	// distributions. Create new SwigIclRun object.
	// nitems		Total number of test items.
	// ntheta		Number of discrete categories in latent variable distribution
	// ngroup		Number of examinee groups. If > 1 than multiple group estimation is used
	// ilist4		Vector containing number of response categories for each item (1 = dichotomous
	//				item modeled by dichotomous model, >2 = polytomous item with that number of
	//				response categories)
	// dlist5		Vector containing minimum and maximum points to use for discrete 
	//				latent variable distribution.
	// uniquePoints If nonzero then unique latent distribution points are used for each examinee group
	//YKM : chage int uniquePoints to bool uniquePoints
	void IRTModule::new_items_dist(int nitems, int ntheta, int ngroups, int_vector *ilist4 , double_vector *dlist5, bool uniquePoints)
	{
		const char *funcname = "new_items_dist";	
		double minTheta,maxTheta; // theta minimum and maximum
	
		if (gEtirmRun) delete gEtirmRun;
	
		if (nitems < 0 || ntheta < 0 || ngroups < 0)
		{
			throw InvalidArgument("Invalid negative argument", funcname);
		}
	
		if (gIclDefaults->baseGroup > ngroups)
		{
			throw InvalidArgument("Number of groups specified makes default group invalid",
				funcname);
		}
	
		int_vector *models;
		if (ilist4 != 0)
		{
			if (ilist4->size() != nitems)
			{
				throw InvalidArgument("Number of element in response categories list does not equal number of items",
					funcname);
			}
			models = ilist4;
		}
		else
		{
			// If list of models is not specified use dichotomous model for all items
			models = new int_vector(nitems, 1);
		}
	
		if (dlist5 == 0)
		{
			minTheta = -4.0;
			maxTheta = 4.0;
		}
		else if (dlist5->size() != 2)
		{
			if (!ilist4) delete models;
			throw InvalidArgument("Fifth argument must be a list containing two elements",
				funcname);
		}
		else
		{
			minTheta = (*dlist5)[0];
			maxTheta = (*dlist5)[1];
		}
	
		// If there is only one group uniquePoints should be zero
		if (uniquePoints && ngroups == 1) uniquePoints = 0;

		gEtirmRun = new SwigIclRun(nitems, models, ntheta, ngroups, minTheta, maxTheta, uniquePoints, gIclDefaults);
	
		if (!ilist4) delete models;
	}

	// End run and delete gEtirmRun object
	void IRTModule::delete_items_dist()
	{
		if (gEtirmRun) delete gEtirmRun;
		gEtirmRun = 0;
	}

	// Return base examinee group
	int IRTModule::get_base_group()
	{
		return gIclDefaults->baseGroup;
	}

	// Change base examinee group
	void IRTModule::set_base_group(int group)
	{
		if (group < 0 || (gEtirmRun != 0 && gEtirmRun->numGroups < group))
		{
			throw InvalidArgument("Invalid group", "set_base_group");
		}
	
		gIclDefaults->baseGroup = group;
	}

	// Change value used for D in creating item objects
	void IRTModule::set_default_D(double D)
	{
		gIclDefaults->D = D;
	}

	// Change default item parameter priors used in creating item objects
	//
	// param		Parameter to set prior for (a, b, c)
	// priortype	Prior distribution to use as default (beta, normal, lognormal, none)
	// priorparam	Parameters of default prior distribution.
	void IRTModule::set_default_prior(const char *param, const char *priortype, double_vector *priorparam)
	{
		const char *funcname = "default_prior";
		int i = ParamIndex(param, funcname);
	
		if ((gIclDefaults->priors)[i]) delete (gIclDefaults->priors)[i];
		std::string pstr(priortype);
		(gIclDefaults->priors)[i] = CreatePrior(pstr, *priorparam, funcname);
	}

	// Return type of default prior for one item parameter (normal, lognormal, beta, none)
	//	paramname	Name of parameter (a, b, or c)
	const char * IRTModule::get_default_prior_type(const char *paramname)
	{
		const char *fname = "get_default_prior_type";
		CheckRunInit();
		std::string &priorName = gEtirmRun->returnString;
		int index = ParamIndex(paramname, fname);
	
		PriorVector::iterator pv = gIclDefaults->priors.begin() + index;

		if (*pv != 0)
		{
			priorName = (*pv)->DistributionName();
		}
		else
		{
			priorName = "none";
		}
	
		return priorName.c_str();
	
	
	}

	// Return vector of parameters of default prior for one item parameter
	//	paramname	Name of parameter (a, b, or c)
	double_vector * IRTModule::get_default_prior_param(const char *paramname)
	{
		const char *fname = "get_default_prior_param";
		CheckRunInit();
		int index = ParamIndex(paramname, fname);
	
		PriorVector::iterator pv = gIclDefaults->priors.begin() + index;
	
		double_vector *v;
		if (*pv != 0)
		{
			int n = (*pv)->NumParameters();
	
			v = new double_vector(n);
	
			*v = (*pv)->GetParameters();
		}
		else
		{
			v = new double_vector();
		}
	
		return v;
	}


	// Change default model used in creating item objects
	//
	// model		Model to set as default (ThreePL, TwoPL, OnePL)
	void IRTModule::set_default_model_dichotomous(IRTModel irtModel)
	{	
		gIclDefaults->dichModel = irtModel;
	}

	// Change default model used in creating item objects
	//
	// model		Model to set as default (ThreePL, TwoPL, OnePL)
	void IRTModule::set_default_model_polytomous(IRTModel irtModel)
	{	
		gIclDefaults->polyModel = irtModel;
	}


	// Set model to use for an item
	//	itemno	Number of item to set model for (1-offset)
	//	model	Model to use for item (3PL, 2PL, 1PL)
	void IRTModule::item_set_model(int itemno, IRTModel irtModel)
	{
		const char *fname = "item_set_model";
		CheckRunInit();
		gEtirmRun->CheckItemNumber(itemno, fname);

		item_type *item = (gEtirmRun->items)[itemno-1];
	
		int ncat = item->NumRespCat();
	
		// Don't let polytomous item be assigned a dichotomous model
		if (ncat > 2)
		{
			if (irtModel == ThreePL || irtModel == TwoPL || irtModel == OnePL)
			{
				throw InvalidArgument("Cannot assign a dichotomous model to a polytomous item", fname);
			}
		}
	
		// erase current information about item
		int index = item->Index();
		item->DeletePriors();
		delete item;
	
		// Recreate item information with new model
		item = CreateItem(irtModel, index, ncat, gIclDefaults->priors, gIclDefaults->D, 
			gEtirmRun->latentDist, fname);
		(gEtirmRun->items)[itemno-1] = item;
	
		// Recreate Uncmin object associated with item
		Uncmin<RealVector, RealMatrix, item_type> *uncmin = (gEtirmRun->minProc)[itemno-1];
		delete uncmin;
		uncmin = new Uncmin<RealVector, RealMatrix, item_type>(item->NumParameters());
		(gEtirmRun->minProc)[itemno-1] = uncmin;

	}

	/*
		item_3PL_starting_values
	
		Calculate item parameter starting values for 3PL, 2PL, and 1PL items.
		Returns number of items for which minimization procedure
		used to compute starting values failed.
	
		use_all
			If non-zero calculate PROX estimates for all examinees and items, even examinees
			who get all items right or all items wrong, and even for items answered all correctly 
			or all incorrectly. (default = 0)
			 
		item_all
			If non-zero all items are used to compute PROX estimates of examinee
			ability, otherwise only items for which starting values are computed
			are used. (default = 0)
	
		itemno3
			Vector of integers giving items numbers (1-offset)
			for which starting values are calculated. All these items must
			be modeled by the 3PL, 2PL, or 1Pl model. Note all 3PL items are
			still used to compute PROX estimates of examinee abilities which
			are used to compute the starting values for the items given by
			itemno3 if item_all is non-zero.
			If null compute starting values for all 3PL items (default = null).
		
	*/
	// YKM : change int use_all and int item_all to bool use_all and bool item_all, respectivelly
	int IRTModule::item_3PL_starting_values(bool use_all, bool item_all, int_vector *itemno3)
	{
		const char *fname = "item_starting_values";
		CheckRunInit();
		gEtirmRun->CheckExaminees(fname);

		int_vector *itemnos, *itemall3pl = 0;

		if (itemno3)
		{
			itemnos = itemno3;

			// Check that all items in list are 3PL, 2PL, or 1PL items
			ItemVector &all_items = gEtirmRun->items;
			int_vector::iterator ii = itemnos->begin();
			for (int i = itemnos->size(); i--; ++ii)
			{
				IRTModel m = all_items[*ii-1]->Model();
				if (m != ThreePL && m != TwoPL && m != OnePL)
				{
					throw InvalidArgument("Item number specified for which starting values can not be computed",
						"item_3PL_starting_values");
				}
			}
		}

		// Find item numbers for all 3PL, 2PL, and 1Pl items.
		if (item_all || !itemno3)
		{
			// Use std::vector to accumulate item numbers
			// since int_vector does not have a push_back function.
			std::vector<int> num;

			ItemVector::iterator ii = (gEtirmRun->items).begin();
			for (int i = (gEtirmRun->items).size(); i--; ++ii)
			{
				IRTModel m = (*ii)->Model();
				if (m == ThreePL || m == TwoPL || m == OnePL)
				{
					num.push_back((*ii)->Index()+1);
				}
			}
		
			itemall3pl = new int_vector(num.begin(), num.end());
		
			if (!itemno3) itemnos = itemall3pl;
		}
		
		UncminVector minsub = ItemSubset(itemnos,gEtirmRun->minProc, fname);
		ItemVector itemsub = ItemSubset(itemnos,gEtirmRun->items, fname);
	
		int n;
		try
		{
			if (item_all)
			{
				/* Create set containing pointers to items for which starting values are computed */
				std::set<item_type *> start_items;
				ItemVector::iterator ii = itemsub.begin();
				for (int i = itemsub.size(); i--; ++ii)
				{
					start_items.insert(*ii);
				}
	
				/* Vector of all 3PL items */
				ItemVector itemall = ItemSubset(itemall3pl, gEtirmRun->items, fname);

				n = StartingValues3PL<examinee_type::response_iterator, UncminVector::iterator,
					std::vector<examinee_type *>::iterator, item_type, std::vector<item_type *>::iterator>
					(minsub.begin(), (gEtirmRun->examinees).begin(), (gEtirmRun->examinees).end(), 
					itemall.begin(), itemall.end(), 
					item_type::NotPresentedResponse(), use_all,
					gIclDefaults->baseGroup, &start_items);
			}
			else
			{
			
				n = StartingValues3PL<examinee_type::response_iterator, UncminVector::iterator,
					std::vector<examinee_type *>::iterator, item_type, std::vector<item_type *>::iterator>
					(minsub.begin(), (gEtirmRun->examinees).begin(), (gEtirmRun->examinees).end(), 
					itemsub.begin(), itemsub.end(),  
					item_type::NotPresentedResponse(), use_all, gIclDefaults->baseGroup);
			}
		}
		catch (...)
		{
			if (itemall3pl) delete itemall3pl;
			throw;
		}
		if (itemall3pl) delete itemall3pl;
		return n;

	}

	// Set seed of random number generator used for simulating from a normal distribution
	void IRTModule::normal_seed(unsigned long seed)
	{
		/*YKM commented
		// convert to type used for seed in boost random number library
		boost::uint32_t boost_seed = seed;
	
		if (gIclDefaults->base_rand_normal)
		{
			gIclDefaults->base_rand_normal->seed(boost_seed);
		}
		else
		{
			gIclDefaults->base_rand_normal = new random_type(boost_seed);
			gIclDefaults->rand_normal = new boost::normal_distribution<random_type>(*(gIclDefaults->base_rand_normal));
		}*/
	}

	// Return a generated number from a normal distribution with the specified mean and s.d.
	double IRTModule::rand_normal(double mean, double sd)
	{
		/* YKM commented
		if (!gIclDefaults->base_rand_normal)
		{
			gIclDefaults->base_rand_normal = new random_type();
			gIclDefaults->rand_normal = new boost::normal_distribution<random_type>(*(gIclDefaults->base_rand_normal));
		}
	
		return (*gIclDefaults->rand_normal)() * sd + mean;
		*/
		return 0.0;

	}
	
	//##################################################################################################################################################################
	
	 /*!
	  \brief
	  Assigns missing response code to identify examinees who did not respond to an item.
	 */
	void IRTModule::set_missing_resp(char nr)
	{
	  if (gEtirmRun)
	  {
	    throw RuntimeError("The not presented response can only be set before any items are initialized", 0);
	  }
	  Item<Real>::SetNotPresentedResponse(nr);
	}

	/*!
	  \brief
	  Returns the number of items.
	 */
	int IRTModule::num_items()
	{
	  CheckRunInit();

	  return gEtirmRun->numItems;
	}

	/*!
	  \brief
	  Returns number of categories of the discrete theta distribution.
	 */
	int IRTModule::num_latent_dist_points()
	{
	  CheckRunInit();

	  return gEtirmRun->numLatentCat;
	}

	/*!
	  \brief
	  Returns number of examinee groups.
	 */
	int IRTModule::num_groups()
	{
	  CheckRunInit();

	  return gEtirmRun->numGroups;
	}

	/*!
	  \brief
	  Returns number of examinees.
	 */
	int IRTModule::num_examinees()
	{
	  CheckRunInit();

	  return (gEtirmRun->examinees).size();
	}

	/*!
	  \brief
	  Returns the name of model used for an item.

	  \section function_args Function Parameters

	  \param[in]  itemno Item number (1-offset).
	 */
	const char * IRTModule::item_get_model(int itemno)
	{

	  CheckRunInit();
	  std::string &modelName = gEtirmRun->returnString;
	  gEtirmRun->CheckItemNumber(itemno, "item_get_model");

	  item_type *item = (gEtirmRun->items)[itemno-1];
	  modelName = item->ModelName();

	  return modelName.c_str();
	}

	/*!
	  \brief
	  Assigns a value to one item parameter of an item.

	  \section function_args Function Parameters

	  \param[in] paramno  1-offset index of parameter in parameter vector for item.
	  \param[in] itemno   Item number (1-offset).
	  \param[in] paramvalue	Value parameter is set to.
	 */
	void IRTModule::item_set_param(int paramno, int itemno, double paramvalue)
	{
	  const char *fname = "item_set_param";
	  CheckRunInit();
	  int index = paramno - 1;

	  gEtirmRun->CheckItemNumber(itemno, fname);
	  item_type *item = (gEtirmRun->items)[itemno-1];
	  CheckItemParam(item, index, fname, gEtirmRun);

	  /* Check that value of parameter has nonzero prior density */
	  item_type::prior_iterator pri = item->PriorsIterator() + index;
	  if (*pri && (*pri)->ZeroDensity(paramvalue))
	  {
	    char merr[100];
	    std::sprintf(merr, "Parameter specified for item %d has zero prior density", itemno);
	    gEtirmRun->returnString = merr;
	    throw RuntimeError((gEtirmRun->returnString).c_str(), fname);
	  }

	  item_type::param_iterator pv = item->ParametersIterator();

	  pv[index] = paramvalue;
	}

	/*!
	  \brief
	  Assigns values to all item parameters of an item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Item number (1-offset).
	  \param[in]  *params  Pointer to params vector containing the values to set the 
	      parameters to.

	  Note: The order of the parameters in the params vector is: 
	  a (if there is an a parameter), one or more b's, and
	  c (if there is a c parameter).
	 */
	void IRTModule::item_set_params(int itemno, double_vector *params)
	{
	  const char *fname = "item_set_params";
	  CheckRunInit();
	  gEtirmRun->CheckItemNumber(itemno, fname);

	  item_type *item = (gEtirmRun->items)[itemno-1];

	  // Assign new parameters
	  item_type::param_iterator pv = item->ParametersIterator();
	  double_vector::iterator iv = params->begin();
	  int n;
	  for (n = item->NumParameters(); n--; ++pv, ++iv)
	  {
	    *pv = *iv;
	  }

	  /* Check that values of all parameters have nonzero prior density */
	  iv = params->begin();
	  item_type::prior_iterator pri = item->PriorsIterator();
	  for (n = item->NumParameters(); n--; ++pri, ++iv)
	  {
	    if (*pri && (*pri)->ZeroDensity(*iv))
	    {
	      char merr[100];
	      std::sprintf(merr, "Parameter specified for item %d has zero prior density", itemno);
	      throw InvalidArgument(merr, fname);
	    }
	  }
	}

	/*!
	  \brief
	  Assigns values to all fixed and estimated item parameters of an item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Item number (1-offset).
	  \param[in]  *params  Pointer to params vector containing the values to set the 
	      parameters to.

	  Note: The order of the parameters in the params vector is: 
	  a (if there is an a parameter), one or more b's, and
	  c (if there is a c parameter).
	 */
	void IRTModule::item_set_all_params(int itemno, double_vector *params)
	{
	  const char *fname = "item_set_all_params";
	  CheckRunInit();
	  gEtirmRun->CheckItemNumber(itemno, fname);

	  item_type *item = (gEtirmRun->items)[itemno-1];

	  // Assign new fixed and estimated parameters
	  item->SetAllParameters(*params);

	  /* Check that values of estimated parameters all have nonzero prior density */
	  item_type::param_iterator iv = item->ParametersIterator();
	  item_type::prior_iterator pri = item->PriorsIterator();
	  for (int n = item->NumParameters(); n--; ++pri, ++iv)
	  {
	    if (*pri && (*pri)->ZeroDensity(*iv))
	    {
	      char merr[100];
	      std::sprintf(merr, "Parameter specified for item %d has zero prior density", itemno);
	      gEtirmRun->returnString = merr;
	      throw InvalidArgument((gEtirmRun->returnString).c_str(), fname);
	    }
	  }
	}

	/*!
	  \brief
	  Returns value of one item parameter of an item.

	  \section function_args Function Parameters

	  \param[in]  paramno  1-offset index of the parameter in the item's parameter vector.
	  \param[in]  itemno  Item number (1-offset).

	  Note: The order of the parameters in the params vector is: 
	  a (if there is an a parameter), one or more b's, and
	  c (if there is a c parameter).
	 */
	double IRTModule::item_get_param(int paramno, int itemno)
	{
	  const char *fname = "item_get_param";
	  CheckRunInit();
	  int index = paramno - 1;

	  gEtirmRun->CheckItemNumber(itemno, fname);
	  item_type *item = (gEtirmRun->items)[itemno-1];
	  CheckItemParam(item, index, fname, gEtirmRun);

	  item_type::param_iterator pv = item->ParametersIterator();

	  return pv[index];
	}

	/*!
	  \brief
	  Returns vector of all estimated item parameters values of an item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Item number (1-offset).

	  Note: The order of the parameters in the params vector is: 
	  a (if there is an a parameter), one or more b's, and
	  c (if there is a c parameter).
	 */
	double_vector * IRTModule::item_get_params(int itemno)
	{
	  CheckRunInit();
	  gEtirmRun->CheckItemNumber(itemno, "item_get_params");

	  item_type *item = (gEtirmRun->items)[itemno-1];

	  int n = item->NumParameters();
	  double_vector *params = new double_vector(n);

	  item_type::param_iterator pv = item->ParametersIterator();
	  double_vector::iterator iv = params->begin();
	  for (; n--; ++pv, ++iv)
	  {
	    *iv = *pv;
	  }

	  return params;
	}

	/*!
	  \brief
	  Returns vector of all fixed and estimated item parameters values of an item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Item number (1-offset).

	  Note: The order of the parameters in the params vector is: 
	  a (if there is an a parameter), one or more b's, and
	  c (if there is a c parameter).
	 */
	double_vector * IRTModule::item_get_all_params(int itemno)
	{
	  CheckRunInit();
	  gEtirmRun->CheckItemNumber(itemno, "item_get_all_params");

	  item_type *item = (gEtirmRun->items)[itemno-1];

	  RealVector allparam = item->GetAllParameters();
	  int n = allparam.size();
	  double_vector *params = new double_vector(n);

	  RealVector::iterator pv = allparam.begin();
	  double_vector::iterator iv = params->begin();
	  for (; n--; ++pv, ++iv)
	  {
	    *iv = *pv;
	  }

	  return params;
	}

	/*!
	  \brief
	  Returns the number of parameters of an item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Item number (1-offset).
	 */
	int IRTModule::item_num_params(int itemno)
	{
	  CheckRunInit();
	  gEtirmRun->CheckItemNumber(itemno, "item_num_params");

	  item_type *item = (gEtirmRun->items)[itemno-1];

	  return item->NumParameters();
	}

	/*!
	  \brief
	  Returns the number of response categories of an item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Item number (1-offset).
	 */
	int IRTModule::item_num_resp_cat(int itemno)
	{
	  CheckRunInit();
	  gEtirmRun->CheckItemNumber(itemno, "item_num_resp_cat");

	  item_type *item = (gEtirmRun->items)[itemno-1];

	  return item->NumRespCat();
	}

	/*!
	  \brief
	  Assigns prior distribution parameters for one item parameter.

	  \section function_args Function Parameters

	  \param[in]  paramno 1-offset index of parameter in parameter vector for item.
	  \param[in]  itemno  Number of item for which prior is set (1-offset).
	  \param[in]  *priortype Pointer to type of prior distribution ("normal", "lognormal", "beta", "none").
	  \param[in]  *dlist Pointer to vector of prior distribution parameters.
	 */
	void IRTModule::item_set_prior(int paramno, int itemno, char *priortype, double_vector *dlist)
	{
	  const char *funcname = "item_set_prior";
	  CheckRunInit();
	  int index = paramno - 1;

	  gEtirmRun->CheckItemNumber(itemno, funcname);
	  item_type *item = (gEtirmRun->items)[itemno-1];
	  CheckItemParam(item, index, funcname, gEtirmRun);

	  item_type::prior_iterator pv = item->PriorsIterator() + index;

	  if (*pv)
	    delete *pv;

	  std::string pstr(priortype);
	  *pv = CreatePrior(pstr, *dlist, funcname);
	}

	/*!
	  \brief
	  Returns type of prior distribution ("normal", "lognormal", "beta", "none")
	  for one item parameter.

	  \section function_args Function Parameters

	  \param[in]  paramno 1-offset index of parameter in parameter vector for item.
	  \param[in]  itemno  Number of item for which prior is returned (1-offset).
	 */
	const char * IRTModule::item_get_prior_type(int paramno, int itemno)
	{
	  const char *fname = "item_get_prior_type";
	  CheckRunInit();
	  std::string &priorName = gEtirmRun->returnString;
	  int index = paramno - 1;

	  gEtirmRun->CheckItemNumber(itemno, fname);
	  item_type *item = (gEtirmRun->items)[itemno-1];
	  CheckItemParam(item, index, fname, gEtirmRun);

	  item_type::prior_iterator pv = item->PriorsIterator() + index;

	  if (*pv != 0)
	  {
	    priorName = (*pv)->DistributionName();
	  }
	  else
	  {
	    priorName = "none";
	  }

	  return priorName.c_str();
	}

	/*!
	  \brief
	  Returns vector of prior distribution parameters for one item parameter.

	  \section function_args Function Parameters

	  \param[in]  paramno 1-offset index of parameter in parameter vector for item.
	  \param[in]  itemno  Number of item for which prior is returned (1-offset).
	 */
	double_vector * IRTModule::item_get_prior_param(int paramno, int itemno)
	{
	  const char *fname = "item_get_prior_param";
	  CheckRunInit();
	  int index = paramno - 1;

	  gEtirmRun->CheckItemNumber(itemno, fname);
	  item_type *item = (gEtirmRun->items)[itemno-1];
	  CheckItemParam(item, index, fname, gEtirmRun);

	  item_type::prior_iterator pv = item->PriorsIterator() + index;

	  double_vector *v;
	  if (*pv != 0)
	  {
	    int n = (*pv)->NumParameters();

	    v = new double_vector(n);

	    *v = (*pv)->GetParameters();
	  }
	  else
	  {
	    v = new double_vector();
	  }

	  return v;
	}

	/*!
	  \brief
	  Returns vector of response counts in each response category of an item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Number of item for which response counts are returned (1-offset).
	  \param[in]  group Group to return counts for (1-offset). Specify group = 0 to return 
	      counts across all groups.
	 */
	double_vector * IRTModule::item_cat_counts(int itemno, int group)
	{
	  const char *fname = "item_cat_counts";
	  CheckRunInit();

	  gEtirmRun->CheckItemNumber(itemno, fname);
	  if (group != 0)
	    gEtirmRun->CheckGroup(group, "item_resp_count");

	  ItemRespCounts *stats = (gEtirmRun->itemStats)[itemno-1];

	  double_vector *v = new double_vector(stats->CategoryCounts(group));

	  return v;
	}

	/*!
	  \brief
	  Returns the number of examinees responding to an item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Number of item for which response counts are returned (1-offset).
	  \param[in]  group Group to return counts for (1-offset). Specify group = 0 to return 
	      count across all groups.
	 */
	double IRTModule::item_resp_count(int itemno, int group)
	{
	  CheckRunInit();

	  gEtirmRun->CheckItemNumber(itemno, "item_resp_count");
	  if (group != 0)
	    gEtirmRun->CheckGroup(group, "item_resp_count");

	  return (gEtirmRun->itemStats)[itemno-1]->RespCount(group);
	}
	/*!
	  \brief
	  Transforms the parameter estimates of an item to a different latent variable scale.

	  Returns zero if parameters are successfully scaled, or returns nonzero if
	  scaling would result in an invalid parameter value.   

	  \section function_args Function Parameters

	  \param[in]  itemno  Number of item for which to transform parameters (1-offset).
	  \param[in]	slope  Slope of scale transformation to apply to item parameters.
	  \param[in]	intercept	Intercept of scale transformation to apply to item parameters.
	  \param[in]  ignorePriorError  If true then do not report an error if transformed parameter
	      has zero density in prior used for that parameter.
	 */
	int IRTModule::item_scale_params(int itemno, double slope, double intercept, bool ignorePriorError) 
	{ // Retyped ignorePriorErro from "int" to "bool", ww, 2-24-2008.
	  CheckRunInit();

	  gEtirmRun->CheckItemNumber(itemno, "item_scale_params");
	  item_type *item = (gEtirmRun->items)[itemno-1];

	  return item->ScaleParameters(slope, intercept, ignorePriorError);
	}

	/*!
	  \brief
	  Returns the probability that an examinee with a particular theta value will give a 
	  particular response to a particular item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Number (1-based) of item for which response probability is returned.
	  \param[in]  response  Integer representing response, where a response in the first
	      response category is 0, a response in the second response category is 1, etc.
	  \param[in]  theta	Value of latent variable for which response probability is returned.
	 */
	double IRTModule::item_prob_resp(int itemno, int response, double theta)
	{
	  CheckRunInit();
	  gEtirmRun->CheckItemNumber(itemno, "item_prob_resp");

	  item_type *item = (gEtirmRun->items)[itemno-1];

	  return item->ProbResp(item->IndexResponse(response), theta);
	}

	/*!
	  \brief
	  Returns a vector of values of the test characteristic curve, given a vector of theta values and
	  a selection of items. 

	  It is assumed for each item that the lowest response category corresponds to a score of 0, the
	  second response category corresponds to a score of 1, etc.

	  \section function_args Function Parameters

	  \param[in]  *thetas	Pointer to vector of theta values.
	  \param[in]  *ilist2	Pointer to vector of item sequence numbers (1-based).
	 */
	double_vector * IRTModule::test_characteristic_curve(double_vector *thetas, int_vector *ilist2)
	{
	  CheckRunInit();

	  double_vector tcc(thetas->size(), 0.0);

	  // Vector to hold response probabilities for an item.
	  double_vector probs;

	  ItemVector *sitems;
	  if (ilist2)
	  {
	    sitems = new ItemVector(ItemSubset(ilist2, gEtirmRun->items, "test_characteristic_curve"));
	  }
	  else
	  {
	    sitems = &(gEtirmRun->items);
	  }

	  // loop over items
	  ItemVector::iterator ii = sitems->begin();
	  for (int i = sitems->size(); i--; ++ii)
	  {
	    int ncat = (*ii)->NumRespCat();
	    probs.newsize(ncat);

	    double_vector::iterator icc = tcc.begin();
	    double_vector::iterator it = thetas->begin();
	    for (int j = tcc.size(); j--; ++icc, ++it) // loop over thetas
	    {
	      (*ii)->ProbRespAll(*it, probs.begin());

	      // first response category corresponds to 0, second response category
	      // corresponds to 1, etc.
	      double dresp = 1.0; // Start with response 1, 0 is skipped
	      double sum = 0.0;
	      double_vector::iterator ip = probs.begin()+1;
	      for (int k = ncat-1; k--; ++ip, ++dresp)
	      {
	        sum += dresp * *ip;
	      }

	      // Add expected score for item to TCC
	      *icc += sum;
	    }
	  }

	  if (ilist2)
	    delete sitems;

	  return new double_vector(tcc);
	}

	/*!
	  \brief
	  Assigns quadrature point values to each discrete category of the discrete latent variable distribution.

	  \section function_args Function Parameters

	  \param[in]  *dlist Pointer to vector of quadrature points of the latent variable distribution.
	  \param[in]  group Group to use (1, 2, ..., number of groups), default = 1.
	 */
	void IRTModule::dist_set_points(double_vector *dlist, int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "dist_set_points");

	  int n = dlist->size();
	  if (n != (gEtirmRun->latentDist).size())
	  {
	    throw RuntimeError("Number of points does not match number of categories in latent distribution",
	        "dist_set_points");
	  }

	  DiscreteLatentDist<Real>::point_iterator ip = (gEtirmRun->latentDist).begin_points(group);

	  double_vector::iterator id = dlist->begin();
	  for (; n--; ++ip, ++id)
	  {
	    *ip = *id;
	  }
	}

	/*!
	  \brief
	  Assigns a quadrature point value to one discrete category of the discrete latent variable distribution.

	  \section function_args Function Parameters

	  \param[in]  index Index of latent variable category (1, 2, ..., number of categories).
	  \param[in]  p Value to assign to point.
	  \param[in]  group Group to use (1, 2, ..., number of groups), default = 1.
	 */
	void IRTModule::dist_set_point(int index, double p, int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "dist_set_point");

	  if (index < 1 || index> (gEtirmRun->latentDist).size())
	  {
	    throw RuntimeError("Invalid category of latent variable distribution","dist_set_point");
	  }

	  DiscreteLatentDist<Real>::point_iterator ip = (gEtirmRun->latentDist).begin_points(group);

	  ip[index-1] = p;
	}

	/*!
	  \brief
	  Returns vector of quadrature point values of the discrete latent variable distribution.

	  \section function_args Function Parameters

	  \param[in]  group Group to use (1, 2, ..., number of groups), default = 1.
	 */
	double_vector * IRTModule::dist_get_points(int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "dist_get_points");

	  DiscreteLatentDist<Real>::point_iterator ip = (gEtirmRun->latentDist).begin_points(group);

	  int n = (gEtirmRun->latentDist).size();
	  double_vector *outv = new double_vector(n);
	  double_vector::iterator id = outv->begin();
	  for (; n--; ++ip, ++id)
	  {
	    *id = *ip;
	  }

	  return outv;
	}

	/*!
	  \brief
	  Returns the quadrature point value of one discrete category of the discrete latent variable distribution.

	  \section function_args Function Parameters

	  \param[in]  index Index of latent variable category (1, 2, ..., number of categories).
	  \param[in]  group Group to use (1, 2, ..., number of groups), default = 1.
	 */
	double IRTModule::dist_get_point(int index, int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "dist_get_point");

	  if (index < 1 || index> (gEtirmRun->latentDist).size())
	  {
	    throw RuntimeError("Invalid category of latent variable distribution","dist_get_point");
	  }

	  DiscreteLatentDist<Real>::point_iterator ip = (gEtirmRun->latentDist).begin_points(group);

	  return ip[index-1];
	}

	/*!
	  \brief
	  Assigns quadrature weights to the quadrature points of the latent variable distribution.

	  \section function_args Function Parameters

	  \param[in]  *dlist Pointer to vector of quadrature weights of the latent variable distribution.
	  \param[in]  group Group to use (1, 2, ..., number of groups), default = 1.
	 */
	void IRTModule::dist_set_probs(double_vector *dlist, int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "dist_set_probs");

	  int n = dlist->size();
	  if (n != (gEtirmRun->latentDist).size())
	  {
	    throw RuntimeError("Number of probabilities does not match number of categories in latent distribution",
	        "dist_set_probs");
	  }

	  // assign probabilities
	  DiscreteLatentDist<Real>::weight_iterator ip = (gEtirmRun->latentDist).begin_weights(group);
	  double_vector::iterator id = dlist->begin();
	  Real sum = 0.0;
	  for (; n--; ++ip, ++id)
	  {
	    *ip = *id;
	    sum += *ip;
	  }

	  // standardize probabilities to they sum to 1
	  ip = (gEtirmRun->latentDist).begin_weights(group);
	  for (n = dlist->size(); n--; ++ip)
	  {
	    *ip /= sum;
	  }
	}

	/*!
	  \brief
	  Assigns the quadrature weight to one discrete category of the discrete latent variable distribution.

	  \section function_args Function Parameters

	  \param[in]  index Index of latent variable category (1, 2, ..., number of categories).
	  \param[in]  w Value to assign to weight.
	  \param[in]  group Group to use (1, 2, ..., number of groups), default = 1.
	 */
	void IRTModule::dist_set_prob(int index, double w, int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "dist_set_prob");

	  if (index < 1 || index> (gEtirmRun->latentDist).size())
	  {
	    throw RuntimeError("Invalid category of latent variable distribution","dist_set_prob");
	  }

	  DiscreteLatentDist<Real>::weight_iterator ip = (gEtirmRun->latentDist).begin_weights(group);

	  ip[index-1] = w;
	}

	/*!
	  \brief
	  Returns vector of quadrature weights of the latent variable distribution for one group of
	  examinees.

	  \section function_args Function Parameters

	  \param[in]  group Group to use (1, 2, ..., number of groups), default = 1.
	 */
	double_vector *IRTModule::dist_get_probs(int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "dist_get_probs");

	  DiscreteLatentDist<Real>::weight_iterator ip = (gEtirmRun->latentDist).begin_weights(group);

	  int n = (gEtirmRun->latentDist).size();
	  double_vector *outv = new double_vector(n);
	  double_vector::iterator id = outv->begin();
	  for (; n--; ++ip, ++id)
	  {
	    *id = *ip;
	  }

	  return outv;

	}

	/*!
	  \brief
	  Returns the quadrature weight of one discrete category of the latent variable 
	  distribution for one group of examinees.

	  \section function_args Function Parameters

	  \param[in]  index Index of latent variable category (1, 2, ..., number of categories).
	  \param[in]  group Group to use (1, 2, ..., number of groups), default = 1.
	 */
	double IRTModule::dist_get_prob(int index, int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "dist_get_prob");

	  if (index < 1 || index> (gEtirmRun->latentDist).size())
	  {
	    throw RuntimeError("Invalid category of latent variable distribution","dist_get_prob");
	  }

	  DiscreteLatentDist<Real>::weight_iterator ip = (gEtirmRun->latentDist).begin_weights(group);

	  return ip[index-1];
	}

	/*!
	  \brief
	  Transforms the quadrature points of latent variable distribution to a new scale.

	  \section function_args Function Parameters

	  \param[in]  slope Slope parameter of linear scale transformation.
	  \param[in]  intercept Intercept parameter of linear scale transformation.
	 */
	void IRTModule::dist_transform(double slope, double intercept)
	{
	  CheckRunInit();

	  (gEtirmRun->latentDist).Transform(slope, intercept);
	}

	/*!
	  \brief
	  Scales to the quadrature points of latent variable distribution to yield a 
	  specfic mean and standard deviation in one group.

	  \section function_args Function Parameters

	  \param[in]  mean  Mean of reference group after scaling the latent distribution.
	  \param[in]  sd  Standard deviation of reference group after scaling the latent distribution.
	  \param[in]  group Group to use as reference group (1, 2, ..., number of groups).
	 */
	double_vector *IRTModule::dist_scale(double mean, double sd, int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "dist_scale");

	  double slope, intercept;
	  (gEtirmRun->latentDist).Scale(mean, sd, group, slope, intercept);

	  double_vector *outv = new double_vector(2);

	  (*outv)[0] = slope;
	  (*outv)[1] = intercept;

	  return outv;
	}

	/*!
	  \brief
	  Returns a vector with the mean and standard deviation of the latent variable distribution 
	  for a selected group of examinees.

	  \section function_args Function Parameters

	  \param[in]  group Group to use as reference group (1, 2, ..., number of groups).
	 */
	double_vector *IRTModule::dist_mean_sd(int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "dist_mean_sd");

	  double mean, sd;
	  (gEtirmRun->latentDist).MeanSD(group, mean, sd);

	  double_vector *outv = new double_vector(2);

	  (*outv)[0] = mean;
	  (*outv)[1] = sd;

	  return outv;
	}

	/*!
	  \brief
	  Returns 1 if unique sets of quadrature points are used for two or more groups of 
	  examinees, or returns 0 otherwise.
	 */
	int IRTModule::dist_unique_points()
	{
	  CheckRunInit();

	  return ((gEtirmRun->latentDist).NumGroupsUnique() > 1 && (gEtirmRun->latentDist).NumGroups() > 1);
	}

	/*!
	  \brief
	  Returns a vector of probabilities for a discrete distribution to approximate a normal distribution 
	  over a set of equally spaced points.

	  Note: The mean and s.d. parameters only effect the quadrature points themselves, so it does not
	  matter which values are passed for the mean and s.d.

	  \section function_args Function Parameters

	  \param[in]  npoints		Number of discrete points in the latent distribution.
	  \param[in]  minPoint	Minimum point of the discrete distribution (must be smaller than the mean).
	  \param[in]  maxPoint	Maximum point of the discrete distribution (must be larger than the mean).
	  \param[in]  mean  Mean of the normal distribution.
	  \param[in]  sd  Standard deviation of the normal distribution.
	 */
	double_vector * IRTModule::normal_dist_prob(int npoints, double minPoint, double maxPoint, double mean, double sd)
	{
	  RealVector points(npoints);
	  RealVector *prob = new RealVector(npoints);

	  DiscreteNormalDist(npoints, minPoint, maxPoint, points.begin(), prob->begin(), mean, sd);

	  return prob;
	}

	/*!
	  \brief
	  Returns a vector of quadrature points for a discrete distribution over a set of equally-spaced 
	  points.
	  \section function_args Function Parameters

	  \param[in]  npoints   Number of discrete points in the latent distribution.
	  \param[in]  minPoint  Minimum point of the discrete distribution (must be smaller than the mean).
	  \param[in]  maxPoint  Maximum point of the discrete distribution (must be larger than the mean).
	  \param[in]  mean  Mean of the normal distribution.
	  \param[in]  sd  Standard deviation of the normal distribution.
	 */
	double_vector * IRTModule::normal_dist_points(int npoints, double minPoint, double maxPoint, double mean, double sd)
	{
	  RealVector prob(npoints);
	  RealVector *points = new RealVector(npoints);

	  // Probabilities do not matter so make up mean and sd
	  DiscreteNormalDist(npoints, minPoint, maxPoint, points->begin(), prob.begin(), mean, sd);

	  return points;
	}

	/*!
	  \brief
	  CalculateS M-step for a set of items.

	  \section function_args Function Parameters

	  \param[in]  ignore_max_iter Flag: If ignore_max_iter == TRUE, then exceeding the maximum number of 
	      iterations in optimization procedure is NOT a fatal error; if ignore_max_iter == FALSE, then the
	      exceeding the max_iter is treated as a fatal error (default = FALSE).
	  \param[in]  *itemno  Pointer to vector of integers giving item numbers (1-offset)
	     for which the M-step is calculated. If itemno is the null pointer, then all items are used
	     (default = null).
	 */
	int IRTModule::mstep_items(bool ignore_max_iter, int_vector *itemno)
	{  // Retyped ignore_max_iter from "int" to "bool", ww, 2-25-2008.
	  CheckRunInit();

	  if (itemno)
	  {
	    const char *fname = "mstep_items";
	    ItemVector itemsub = ItemSubset(itemno, gEtirmRun->items, fname);
	    UncminVector minsub = ItemSubset(itemno, gEtirmRun->minProc, fname);
	    return MStepItems(minsub.begin(), itemsub.begin(), itemsub.end(), gEtirmRun->mstepMaxDiff, ignore_max_iter);
	  }
	  else
	  {
	    return MStepItems((gEtirmRun->minProc).begin(), (gEtirmRun->items).begin(), (gEtirmRun->items).end(), gEtirmRun->mstepMaxDiff, ignore_max_iter);
	  }
	}

	/*!
	  \brief
	  Returns message from M-step minimization for one item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Integer giving item number (1-offset).
	 */
	int IRTModule::mstep_message(int itemno)
	{
	  CheckRunInit();
	  gEtirmRun->CheckItemNumber(itemno, "mstep_message");

	  return ((gEtirmRun->minProc)[itemno-1])->GetMessage();
	}

	/*!
	  \brief
	  Assigns the maximum number of iterations permitted with optimization procedure
	  in mstep_items for one item.

	  \section function_args Function Parameters

	  \param[in]  itemno  Integer giving item number (1-offset).
	  \param[in]  maxiter  Maximum number of iterations for stepwise optimization.
	 */
	void IRTModule::mstep_max_iter(int itemno, int maxiter)
	{
	  CheckRunInit();
	  gEtirmRun->CheckItemNumber(itemno, "mstep_max_iter");

	  return ((gEtirmRun->minProc)[itemno-1])->SetMaxIter(maxiter);
	}

	/*!
	  \brief
	  Returns the maximum relative difference between parameter estimates in two
	  successive EM iterations, as computed in last call to mstep_items.
	 */
	double IRTModule::mstep_max_diff()
	{
	  CheckRunInit();

	  return gEtirmRun->mstepMaxDiff;
	}

	/*!
	  \brief
	  Executes M-step for the latent distribution in one group.

	  Returns maximum relative difference between new and old probabilities.

	  \section function_args Function Parameters

	  \param[in]  *e  Pointer to not sure what(?!). To-Do: Find out what *e relates to! 
	  \param[in]  group Number of examinee group (1, 2, ..., number of groups).
	 */
	double IRTModule::mstep_dist(Estep *e, int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, "mstep_dist");

	  estep_type::ngroup_iterator i = (e->GetEStep())->GetNGroup(group);
	  RealVector::size_type n = (e->GetEStep())->size();
	  RealVector prob(i, i+n);
	  return (gEtirmRun->latentDist).MStep(prob, group);
	}

	/*!
	  \brief
	  Add an examinee objec to the end of the examinee vector.

	  Returns the examinee number corresponding to this examinee.

	  \section function_args Function Parameters

	  \param[in]  *responses	Pointer to vector of integer-formatted item responses. 
	      A negative integer indicates the examinee did not respond to the item.
	      The integers representing responses are zero-offset (A response in the 
	      first response category is represented as zero).
	  \param[in]  group		Integer from 1 to maximum number of examinee groups
	      giving the group the examinee belongs to.
	  \param[in]  count		Frequency given to response pattern (usually, count = 1).
	 */
	int IRTModule::add_examinee(int_vector *responses, int group, double count)
	{
	  const char *fname = "add_examinee";
	  CheckRunInit();
	  gEtirmRun->CheckGroup(group, fname);

	  int len = responses->size();

	  if (len != gEtirmRun->numItems)
	  {
	    throw RuntimeError("Invalid number of item responses", fname);
	  }

	  examinee_type *e = new examinee_type(len, group);

	  ResponseVector r(len);
	  ItemVector::iterator iitem = (gEtirmRun->items).begin();
	  int_vector::iterator iir = responses->begin();
	  ResponseVector::iterator ir = r.begin();
	  ItemStatsVector::iterator is = (gEtirmRun->itemStats).begin();
	  Response np = Item<Real>::NotPresentedResponse();
	  for (; len--; ++iitem, ++iir, ++ir, ++is)
	  {
	    if (*iir < 0)
	      *ir = np;
	    else
	    {
	      *ir = (*iitem)->IndexResponse(*iir);
	      if (!((*iitem)->ValidResponse(*ir)))
	      {
	        throw RuntimeError("Invalid item response", fname);
	      }
	      (*is)->AddResponse(*iir+1, count, group);
	    }
	  }

	  e->SetResponses(r);
	  e->SetCount(count);

	  (gEtirmRun->examinees).push_back(e);

	  (gEtirmRun->examineeCounts)[0] += count; // total count
	  (gEtirmRun->examineeCounts)[group] += count; // group count
	  return (gEtirmRun->examinees).size();
	}

	/*!
	  \brief
	  Returns pointer to an integer vector of examinee item responses to all items.

	  A non-negative integer indicates an examinee response in the category 
	  represented by the integer (0 = first response category, 1 = second response 
	  category, etc.). A negative integer indicates the examinee did not respond to 
	  the item.

	  \section function_args Function Parameters

	  \param[in]  examineeno	Number of examinee for whom item responses are returned.
	 */
	int_vector * IRTModule::examinee_responses(int examineeno)
	{
	  const char *fname = "examinee_responses";
	  CheckRunInit();
	  gEtirmRun->CheckExamineeNumber(examineeno, fname);

	  examinee_type *examinee = (gEtirmRun->examinees)[examineeno-1];

	  if ((gEtirmRun->items).size() != examinee->NumItems())
	  {
	    throw RuntimeError("Number of examinee responses does not match number of items", fname);
	  }

	  // iterators to first and one past last response
	  examinee_type::response_iterator first = examinee->responses_begin();
	  examinee_type::response_iterator last = examinee->responses_end();

	  int_vector *responses = new int_vector(examinee->NumItems());
	  int_vector::iterator ir = responses->begin();
	  ItemVector::iterator ii = (gEtirmRun->items).begin();
	  Response np = Item<Real>::NotPresentedResponse();
	  while (first != last)
	  {
	    if (*first == np)
	      *ir = -1;
	    else
	      *ir = (*ii)->ResponseIndex(*first);
	    ++ir;
	    ++ii;
	    ++first;
	  }

	  return responses;
	}

	/*!
	  \brief
	  Returns string containing (character) examinee responses to all items.

	  A response in the first response category is '0', a response in the second category is 
	  '1', etc. Because each response is represented by a single character care must be 
	  taken when some items have more than 10 response categories. For example, a response 
	  in the 11th response category would be represented by the ascii character ':', which 
	  is the 10th character greater than '0' in the ascii character sequence.

	  \section function_args Function Parameters

	  \param[in]  examineeno  Number of examinee for whom item responses are returned.
	 */

	const char * IRTModule::examinee_response_str(int examineeno)
	{
	  const char *fname = "examinee_response_str";
	  CheckRunInit();
	  gEtirmRun->CheckExamineeNumber(examineeno, fname);

	  examinee_type *examinee = (gEtirmRun->examinees)[examineeno-1];

	  if ((gEtirmRun->items).size() != examinee->NumItems())
	  {
	    throw RuntimeError("Number of examinee responses does not match number of items", fname);
	  }

	  // iterators to first and one past last response
	  examinee_type::response_iterator first = examinee->responses_begin();
	  examinee_type::response_iterator last = examinee->responses_end();

	  std::string &responseStr = gEtirmRun->returnString;
	  responseStr.clear();
	  ItemVector::iterator ii = (gEtirmRun->items).begin();
	  while (first != last)
	  {
	    responseStr.push_back(Resp2Char(*first, *ii));
	    ++ii;
	    ++first;
	  }

	  return responseStr.c_str();
	}

	/*!
	  \brief
	  Return the number of the group the examinee belongs to.

	  \section function_args Function Parameters

	  \param[in]  examineeno  Number of examinee for whom group membership is to be returned.
	 */

	int IRTModule::examinee_get_group(int examineeno)
	{
	  CheckRunInit();
	  gEtirmRun->CheckExamineeNumber(examineeno, "examinee_group");

	  examinee_type *examinee = (gEtirmRun->examinees)[examineeno-1];

	  return examinee->Group();
	}

	/*!
	  \brief
	  Assigns which group the examinee belongs to.

	  \section function_args Function Parameters

	  \param[in]  examineeno  Number of examinee for whom group membership is to be assigned.
	  \param[in]  group   Group number to be assigned to examinee.
	 */
	void IRTModule::examinee_set_group(int examineeno, int group)
	{
	  CheckRunInit();
	  gEtirmRun->CheckExamineeNumber(examineeno, "examinee_set_group");

	  examinee_type *examinee = (gEtirmRun->examinees)[examineeno-1];

	  return examinee->SetGroup(group);
	}

	/*!
	  \brief
	  Assigns the count (or weight) of the examinee.

	  \section function_args Function Parameters

	  \param[in]  examineeno  Number of examinee whose count (or weight) is to be assigned.
	  \param[in]  count   Count to be assigned to examinee.
	 */
	void IRTModule::examinee_set_count(int examineeno, double count)
	{
	  CheckRunInit();
	  gEtirmRun->CheckExamineeNumber(examineeno, "examinee_set_count");

	  examinee_type *examinee = (gEtirmRun->examinees)[examineeno-1];

	  examinee->SetCount(count);
	}

	/*!
	  \brief
	  Returns the count (or weight) of the examinee.

	  \section function_args Function Parameters

	  \param[in]  examineeno  Number of examinee whose count (or weight) is to be returned.
	 */
	double IRTModule::examinee_get_count(int examineeno)
	{
	  CheckRunInit();
	  gEtirmRun->CheckExamineeNumber(examineeno, "examinee_get_count");

	  examinee_type *examinee = (gEtirmRun->examinees)[examineeno-1];

	  return examinee->Count();
	}

	/*!
	  \brief
	  Assigns posterior distribution of an examinee.

	  \section function_args Function Parameters

	  \param[in]  examineeno  Number of examinee whose posterior is to be assigned.
	  \param[in]  *posterior  Pointer to vector of posterior probabilities of the examinee.
	 */
	void IRTModule::examinee_set_posterior(int examineeno, double_vector *posterior)
	{
	  CheckRunInit();
	  gEtirmRun->CheckExamineeNumber(examineeno, "examinee_set_posterior");

	  if (posterior->size() != (gEtirmRun->latentDist).size())
	  {
	    throw RuntimeError("Invalid number of posterior probabilities",
	        "examinee_set_posterior");
	  }

	  examinee_type *examinee = (gEtirmRun->examinees)[examineeno-1];

	  examinee_type::posterior_vector epost(gEtirmRun->numLatentCat);
	  examinee_type::posterior_vector::iterator iep = epost.begin();
	  double_vector::iterator ip = posterior->begin();

	  // assign probabilities
	  int i;
	  Real sum = 0.0;
	  for (i = gEtirmRun->numLatentCat; i--; ++iep, ++ip)
	  {
	    *iep = *ip;
	    sum += *iep;
	  }

	  // standardize probabilities to sum to 1.0
	  iep = epost.begin();
	  for (i = gEtirmRun->numLatentCat; i--; ++iep)
	  {
	    *iep /= sum;
	  }

	  examinee->SetPosterior(epost);
	}

	/*!
	  \brief
	  Returns pointer to vector of the posterior distribution of an examinee.

	  \section function_args Function Parameters

	  \param[in]  examineeno  Number of examinee whose posterior is to be returned.
	 */ 
	double_vector* IRTModule::examinee_get_posterior(int examineeno)
	{
	  //Embian Inc. char *fname = "examinee_get_posterior";
	  char *fname = (char *)("examinee_get_posterior");
	  CheckRunInit();
	  gEtirmRun->CheckExamineeNumber(examineeno, fname);

	  examinee_type *examinee = (gEtirmRun->examinees)[examineeno-1];
	  int n = examinee->NumLatentVarCat();

	  if (n == 0)
	  {
	    char merr[100];
	    std::sprintf(merr, "Posterior has not been computed for examinee %d", examineeno);
	    gEtirmRun->returnString = merr;
	    throw RuntimeError((gEtirmRun->returnString).c_str(), fname);
	  }

	  double_vector *post = new double_vector(n);

	  examinee_type::posterior_vector::iterator ie = examinee->posterior_begin();
	  double_vector::iterator ip = post->begin();
	  for (; n--; ++ie, ++ip)
	  {
	    *ip = *ie;
	  }

	  return post;
	}

	/*!
	  \brief
	  Returns means of an examinee's posterior distribution.

	  Posterior must have already been computed, for example, by estep_compute.

	  \section function_args Function Parameters

	  \param[in]  examineeno  Number of examinee whose posterior mean is to be returned.
	 */ 
	double IRTModule::examinee_posterior_mean(int examineeno)
	{
	  //Embian Inc. char *fname = "examinee_posterior_mean";
	  char *fname = (char*)("examinee_posterior_mean");
	  CheckRunInit();
	  gEtirmRun->CheckExamineeNumber(examineeno, fname);

	  examinee_type *examinee = (gEtirmRun->examinees)[examineeno-1];
	  int n = examinee->NumLatentVarCat();

	  if (n == 0)
	  {
	    char merr[100];
	    std::sprintf(merr, "Posterior distribution has not been computed for examinee %d", examineeno);
	    gEtirmRun->returnString = merr;
	    throw RuntimeError((gEtirmRun->returnString).c_str(), fname);
	  }

	  examinee_type::posterior_vector::iterator iw = examinee->posterior_begin();
	  lvdist_type::point_iterator ip = (gEtirmRun->latentDist).begin_points();
	  double mean = 0.0;
	  for (; n--; ++iw, ++ip)
	  {
	    mean += *ip * *iw;
	  }

	  return mean;
	}

	/*!
	  \brief
	  Returns an examinee's maximum likelihood estimate of theta.

	  Posterior must have already been computed, for example, by estep_compute.

	  \section function_args Function Parameters

	  \param[in]  examineeno  Number of examinee whose theta MLE is to be returned (1 = first examinee).
	  \param[in]  minTheta	Minimum value of theta estimate.
	  \param[in]  maxTheta	Maximum value of theta estimate.	
	  \param[in]  precision	Length of interval in which MLE is determined to lie. This should be greater 
	      than, roughly, 3.0e-8.
	  \param[in]  *itemno5		Pointer to a list of item numbers of the items to use in computing examinee MLE. 
	      If itemno5 = NULL, then all items will be used.
	 */
	double IRTModule::examinee_theta_MLE(int examineeno, double minTheta, double maxTheta, double precision,
	    int_vector *itemno5)
	{
	  CheckRunInit();
	  gEtirmRun->CheckExamineeNumber(examineeno, "examinee_theta_MLE");

	  if (!gEtirmRun->base_rand_simulate)
	  {
	    gEtirmRun->base_rand_simulate = new random_type();
	    gEtirmRun->rand_simulate = new boost::uniform_01<random_type>(*(gEtirmRun->base_rand_simulate));
	  }

	  ItemVector *sitems;
	  if (itemno5)
	  {
	    sitems = new ItemVector(ItemSubset(itemno5, gEtirmRun->items, "examinee_theta_MLE"));
	  }
	  else
	  {
	    sitems = &(gEtirmRun->items);
	  }

	  double theta =
	      ExamineeThetaMLE<ItemVector::iterator, examinee_type::response_vector::iterator>(minTheta,
	          maxTheta, precision, sitems->begin(), sitems->end(), (gEtirmRun->examinees[examineeno-1])->responses_begin());

	  if (itemno5)
	    delete sitems;

	  return theta;
	}

	/*!
	  \brief
	  Returns total examinee count in an examinee group (1, 2, ...), or across all groups.

	  \section function_args Function Parameters

	  \param[in]  group Group to return count for (1-offset). Specify group = 0 to return 
	      counts across all groups.
	 */
	double IRTModule::examinees_count(int group)
	{
	  CheckRunInit();
	  if (group != 0)
	    gEtirmRun->CheckGroup(group, "examinee_count");

	  return (gEtirmRun->examineeCounts)[group];
	}

	/*!
	  \brief
	  Assigns seed of random number generator used for bootstrap samples.

	  \section function_args Function Parameters

	  \param[in]  seed  Seed value, 0 <= seed <= 4,294,967,295.
	 */
	void IRTModule::bootstrap_seed(unsigned long seed)
	{
	  CheckRunInit();

	  // convert to type used for seed in boost random number library
	  boost::uint32_t boost_seed = seed;

	  if (gEtirmRun->rand_boot)
	  {
	    gEtirmRun->rand_boot->seed(boost_seed);
	  }
	  else
	  {
	    gEtirmRun->rand_boot = new random_type(boost_seed);
	  }
	}

	/*!
	  \brief
	  Generate bootstrap sample of examinees.
	 */
	void IRTModule::bootstrap_sample()
	{
	  /* Embian Inc. commented
	  CheckRunInit();
	  gEtirmRun->CheckExaminees("bootstrap_sample");

	  if (!gEtirmRun->rand_boot)
	    gEtirmRun->rand_boot = new random_type();

	  boost::uniform_int<random_type> urand(*(gEtirmRun->rand_boot), 1, (gEtirmRun->examinees).size());

	  BootstrapSample((gEtirmRun->examinees).begin(), (gEtirmRun->examinees).end(), urand);

	  // Set vector of examinee counts
	  gEtirmRun->examineeCounts = 0.0; // initialize to zero
	  ExamineeVector::iterator ie = (gEtirmRun->examinees).begin();
	  for (int i = (gEtirmRun->examinees).size(); i--; ++ie)
	  {
	    double count = (*ie)->Count();
	    (gEtirmRun->examineeCounts)[0] += count; // total count
	    int group = (*ie)->Group();
	    (gEtirmRun->examineeCounts)[group] += count; // group count
	  }*/
	}

	/*!
	  \brief
	  Assigns seed of random number generator used for simulating item responses.

	  \section function_args Function Parameters

	  \param[in]  seed  Seed value, 0 <= seed <= 4,294,967,295.
	 */  
	void IRTModule::simulate_seed(unsigned long seed)
	{
	  CheckRunInit();

	  // convert to type used for seed in boost random number library
	  boost::uint32_t boost_seed = seed;

	  if (gEtirmRun->base_rand_simulate)
	  {
	    gEtirmRun->base_rand_simulate->seed(boost_seed);
	  }
	  else
	  {
	    gEtirmRun->base_rand_simulate = new random_type(boost_seed);
	    gEtirmRun->rand_simulate = new boost::uniform_01<random_type>(*(gEtirmRun->base_rand_simulate));
	  }
	}

	/*!
	  \brief
	  Simulates item responses (integers) for a specific value of the latent variable.  

	  Returns pointer to integer vector containing item responses, where
	  item responses are integers from 0 to one minus the maximum
	  number of response categories for the item.

	  \section function_args Function Parameters

	  \param[in]  theta Value of latent variable for which item responses are generated.
	  \param[in]  itemno2	List of item numbers of the items for which responses will be
	      simulated.
	 */
	int_vector * IRTModule::simulate_responses(double theta, int_vector *itemno2)
	{
	  CheckRunInit();

	  if (!gEtirmRun->base_rand_simulate)
	  {
	    gEtirmRun->base_rand_simulate = new random_type();
	    gEtirmRun->rand_simulate = new boost::uniform_01<random_type>(*(gEtirmRun->base_rand_simulate));
	  }

	  int_vector *resp;
	  if (itemno2)
	  {
	    resp = new int_vector(itemno2->size());
	    ItemVector sitems = ItemSubset(itemno2, gEtirmRun->items, "simulate_responses");
	    resp->newsize(sitems.size());
	    SimulateResponses(sitems.begin(), sitems.end(), theta, *(gEtirmRun->rand_simulate),
	        resp->begin(), true);
	  }
	  else
	  {
	    resp = new int_vector((gEtirmRun->items).size());
	    SimulateResponses((gEtirmRun->items).begin(), (gEtirmRun->items).end(), theta, *(gEtirmRun->rand_simulate), resp->begin(), true);
	  }

	  return resp;
	}

	/*!
	  \brief
	  Simulates item responses (characters) for a specific value of the latent variable.  

	  Returns simulated responses in a string, where each character of the string gives a 
	  response to one item. The first response category of an item is represented by a 
	  '0', the second response category by a '1', etc. 

	  Because each response is represented by a single character care must be taken when 
	  some items have more than 10 response categories. For example, a response in the 
	  11th response category would be represented by the ascii character ':', which is 
	  the 10th character greater than '0' in the ascii character sequence.

	  \section function_args Function Parameters

	  \param[in]  theta Value of latent variable for which item responses are generated.
	  \param[in]  itemno2 List of item numbers of the items for which responses will be
	      simulated.
	 */
	const char * IRTModule::simulate_response_str(double theta, int_vector *itemno2)
	{
	  CheckRunInit();

	  // Vector of integers representing simulated responses
	  int_vector *iresp = simulate_responses(theta, itemno2);

	  // Resize string to contain item responses
	  std::string &responseStr = gEtirmRun->returnString;
	  responseStr.resize(iresp->size());

	  // Convert integers to characters
	  int_vector::iterator ii = iresp->begin();
	  std::string::iterator ic = responseStr.begin();
	  for (int i = iresp->size(); i--; ++ii, ++ic)
	  {
	    if (*ii < 0)
	    {
	      *ic = item_type::NotPresentedResponse();
	    }
	    else
	    {
	      *ic = *ii + '0';
	    }
	  }

	  delete iresp;

	  return responseStr.c_str();
	}

	/*!
	  \brief
	  Reads item responses from an input record.  

	  The item responses are assumed to be integers from 0 to one minus the maximum
	  number of response categories for the item. A response representing a missing 
	  response is assigned -1.

	  Returns poiner to an integer vector containing thew item responses.

	  \section function_args Function Parameters

	  \param[in]  *line  Pointer to character string to read the item responses from.
	  \param[in]  *offset  Pointer to vector of zero-based offsets of each item response in 'line'.
	  \param[in]  *len   Pointer to vector of field widths for each item response.  
	 */
	int_vector * IRTModule::get_responses(char *line, int_vector *offset, int_vector *len)
	{
	  const char *fname = "get_responses";
	  CheckRunInit();

	  int n = offset->size();

	  if (n != len->size())
	  {
	    throw InvalidArgument("Lengths of offset and length vectors do not match", fname);
	  }

	  int_vector *resp = new int_vector(n);

	  int_vector::iterator ir = resp->begin();
	  int_vector::iterator ioff = offset->begin();
	  int_vector::iterator ilen = len->begin();
	  Response np = Item<Real>::NotPresentedResponse();
	  for (; n--; ++ioff, ++ilen, ++ir)
	  {
	    // Convert response to integer
	    char *pos = line + *ioff + *ilen - 1;
	    *ir = 0;
	    int power = 1;
	    for (int i=*ilen; i--; power *= 10, --pos)
	    {
	      if (*pos == np) // check for missing item response
	      {
	        *ir = -1;
	        break;
	      }
	      else if ( (*pos < '0') && (*pos > '9')) // check that item response is valid // Sytax edited, ww, 2-24-2008.
	      {
	        char errstr[50];
	        std::sprintf(errstr, "Invalid item response: %1c", *pos);
	        gEtirmRun->returnString = errstr;
	        throw RuntimeError((gEtirmRun->returnString).c_str(), fname);
	      }
	      *ir += power * (*pos - '0');
	    }
	  }

	  return resp;
	}

	/*!
	  \brief
	  Read item responses from a string where responses to only some items are
	  present. The responses to the remaining items are assumed to be missing.

	  The item responses are assumed to be integers from 0 to one minus the maximum
	  number of response categories for the item. A response representing a missing 
	  response is assigned -1.

	  Returns poiner to an integer vector containing thew item responses.

	  \section function_args Function Parameters

	  \param[in]  *line  Pointer to character string to read the item responses from.
	  \param[in]  *offset  Pointer to vector of zero-based offsets of each item response in 'line'.
	  \param[in]  *len   Pointer to vector of field widths for each item response.  
	  \param[in]  *items  Pointer to vector 1-offset indices of items for which responses are read.
	 */
	int_vector * IRTModule::get_responses_missing(char *line, int_vector *offset, int_vector *len,
	    int_vector *items)
	{
	  const char *fname = "get_responses_missing";
	  CheckRunInit();

	  int n = offset->size();

	  if (n != len->size() || n != items->size())
	  {
	    throw InvalidArgument("Lengths of vector arguments do not match", fname);
	  }

	  int_vector *resp = get_responses(line, offset, len);

	  // initially assign all missing responses
	  int nall = gEtirmRun->numItems;
	  int_vector *allresp = new int_vector(nall, -1);

	  // Assign responses read
	  int_vector::iterator ir = resp->begin();
	  int_vector::iterator ii = items->begin();
	  for (; n--; ++ii, ++ir)
	  {
	    if ( (*ii < 1) || (*ii > nall)) // Syntax edited, ww, 2-24-2008.
	    {
	      delete resp;
	      throw InvalidArgument("Invalid item number", fname);
	    }
	    (*allresp)[*ii-1] = *ir;
	  }

	  delete resp;

	  return allresp;
	}

} // namespace etirm
