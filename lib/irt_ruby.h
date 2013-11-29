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
	Copyright (C) 2010, Embian Inc.
	
*/

#ifndef SWIGEPDIRM_H_
#define SWIGEPDIRM_H_

#ifdef ETIRM_NO_DIR_PREFIX
#include "swig_etirm.h"
#else
#include "etirm/swig_etirm.h"
#endif

#include <boost/random/normal_distribution.hpp>

#ifndef SWIG
namespace etirm {
#endif

//############################################################################################################################################
// Holds information needed in creating item objects.
// These default values are used when new_etirm is called.
// Also holds random number objects used in generating normal random variables
class DefaultValues
{

public:

	DefaultValues();
	~DefaultValues();

	double D;
		// Logistic function scaling constant. D=1.7 makes logistic similar to
		// normal.

	PriorVector priors;
		// Vector of pointers to priors for a, b, and c parameters.

	IRTModel dichModel;
		// Default model to use for dichotomous items (ThreePL, TwoPL, OnePL)

	IRTModel polyModel;
		// Default model to use for polytomous items (GPCM, PCM)

	int baseGroup;
		// Examinee group to use as base group for which
		// mean and s.d. of latent variable distribution are
		// set at zero and one.

	random_type *base_rand_normal;
		// Base random number generator object for normal random generator

	boost::normal_distribution<random_type> *rand_normal;
		// Normal random number generator object 


};
//############################################################################################################################################

class IRTModule
{
	//############################################################################################################################################
	// Subclass of SwigEtirmRun to implement Icl
	class SwigIclRun : public SwigEtirmRun
	{

	public:

		SwigIclRun(int nitems, int_vector *ncat, int nlatentcat, int ngroups, double minTheta, double maxTheta,
			bool uniquePoints, DefaultValues *gIclDefaults);
		virtual ~SwigIclRun();

	};
	//############################################################################################################################################
	
private:
	SwigEtirmRun *gEtirmRun;
	DefaultValues *gIclDefaults;
public:
	IRTModule();
	virtual ~IRTModule();

	void CheckRunInit();
	
	SwigEtirmRun * GetEtirmObj();
	
	// Return base examinee group
	int get_base_group();

	// Change base examinee group
	void set_base_group(int group);

	// Set logistic scaling constant (D=1.7 make logistic curve close to normal)
	void set_default_D(double D);

	// Change default model used in creating item objects for dichotomous items
	void set_default_model_dichotomous(IRTModel irtModel);

	// Change default model used in creating item objects for dichotomous items
	void set_default_model_polytomous(IRTModel irtModel);

	// Set model used for item
	void item_set_model(int itemno, IRTModel irtModel);

	// Compute starting values for item parameters
	// YKM : change int use_all and int item_all to bool use_all and bool item_all, respectivelly
	int item_3PL_starting_values(bool use_all = 0, bool item_all = 0, int_vector *ilist3 = 0);

	// Initialize run
	// nitems		Total number of test items.
	// ntheta		Number of discrete categories in latent variable distribution
	// ngroup		Number of examinee groups. If > 1 than multiple group estimation is used
	// ilist4		Vector containing number of response categories for each items (1 = dichotomous
	//				item modeled by dichotomous model, >2 = polytomous item with that number of
	//				response categories)
	// dlist5		Vector containing minimum and maximum points to use for discrete 
	//				latent variable distribution.
	// uniquePoints If nonzero then unique latent distribution points are used for each examinee group
	//YKM : chage int uniquePoints to bool uniquePoints
	void new_items_dist(int nitems, int ntheta = 40, int ngroups = 1, int_vector *ilist4 = 0, double_vector *dlist5 = 0, bool uniquePoints = 0);

	// release memory allocated for run in new_epdirm
	void delete_items_dist();

	// Change default item parameter priors used in creating item objects
	void set_default_prior(const char *param, const char *priortype, double_vector *priorparam);

	// Return type of default prior for one item parameter (normal, lognormal, beta, none)
	const char *get_default_prior_type(const char *paramname);

	// Return vector of parameters of default prior for one item parameter
	double_vector *get_default_prior_param(const char *paramname);

	// Set seed of random number generator used for simulating from a normal distribution
	void normal_seed(unsigned long seed);

	// Return a generated number from a normal distribution with the specified mean and s.d.
	double rand_normal(double mean = 0.0, double sd = 1.0);
	
	
	//##################################################################################################################
	
	
	// Set response which indicates an examinee did
	// not respond to an item
	void set_missing_resp(char nr);

	// Return number of items, number of theta points, and number of groups
	// specified in new_etirm
	int num_items();
	int num_latent_dist_points();
	int num_groups();

	// Return number of examinees added by calling add_examinee
	int num_examinees();

	/*** item functions ***/

	// Return name of model used for item
	const char * item_get_model(int itemno);

	// Assign a value to one parameter for one item
	void item_set_param(int paramno, int itemno, double paramvalue);

	// Assign values to all estimated parameters for one item
	void item_set_params(int itemno, double_vector *params);

	// Assign values to all fixed and estimated parameters for one item
	void item_set_all_params(int itemno, double_vector *params);

	// Get one parameter for one item
	double item_get_param(int paramno, int itemno);

	// Get all estimated parameters for one item
	double_vector *item_get_params(int itemno);

	// Get all estimated and fixed parameters for one item
	double_vector *item_get_all_params(int itemno);

	// Return number of parameters for one item
	int item_num_params(int itemno);

	// Return number of response categories for an item.
	int item_num_resp_cat(int itemno);

	// Set prior for one item
	void item_set_prior(int paramno, int itemno, char *priortype, double_vector *dlist4 = 0);

	// Get prior distribution for one parameter and one item
	const char *item_get_prior_type(int paramno, int itemno);

	// Get parameters of prior
	double_vector *item_get_prior_param(int paramno, int item);

	// Get probability of an item response
	double item_prob_resp(int itemno, int response, double theta);

	// Transform item parameters to different IRT scale
	int item_scale_params(int itemno, double slope, double intercept, bool ignorePriorError = 0); 
	// Retyped ignorePriorError from "int" to "bool", ww, 2-26-2008.

	// Return vector of counts of responses in each response category for an item
	// in an examinee group (0=all groups)
	double_vector *item_cat_counts(int itemno, int group = 0);

	// Return counts of number of examinees responding to an item
	// in an examinee group (0=all groups)
	double item_resp_count(int itemno, int group = 0);

	// Compute M-step for indicated items
	double_vector *test_characteristic_curve(double_vector *thetas, int_vector *ilist2 = 0);

	/** Latent distribution functions ***/

	// Set point value for one category of discrete latent variable distribution
	void dist_set_point(int index, double p, int group = 1);

	// Set points of latent discrete distribution
	void dist_set_points(double_vector *dlist, int group = 1);

	// Return value of point for one category of latent variable distribution
	double dist_get_point(int index, int group = 1);

	// Get points of latent discrete distribution
	double_vector *dist_get_points(int group = 1);

	// Set value of the probability for one category of discrete latent variable distribution for one group
	void dist_set_prob(int index, double w, int group = 1);

	// Set probabilities for latent discrete distribution for one group
	void dist_set_probs(double_vector *dlist, int group = 1);

	// Return value of probability for one category of latent variable distribution for one group
	double dist_get_prob(int index, int group = 1);

	// Get probabilities for latent discrete distribution for one group
	double_vector *dist_get_probs(int group = 1);

	// Transform points of latent variable distribution to a new scale
	void dist_transform(double slope, double intercept);

	// Scale points to give specific mean and sd in one group
	// Returns slope and intercept of linear scale transformation
	double_vector *dist_scale(double mean, double sd, int group = 1);

	// Return mean and sd of distribution in one group
	double_vector *dist_mean_sd(int group = 1);

	// Returns 1 if unique discrete latent distribution points are used for
	// different examinee groups and number of examinee groups is greater than 1, otherwise returns 0
	int dist_unique_points();

	// Return vector of probabilities for discrete distribution to approximate a normal distribution
	double_vector *normal_dist_prob(int npoints, double minPoint, double maxPoint,
	double mean = 0.0, double sd = 1.0);

	// Return vector of points for discrete distribution to approximate a normal distribution
	double_vector *normal_dist_points(int npoints, double minPoint, double maxPoint, double mean = 0.0,
	double sd = 1.0);

	/*** M-step functions ***/

	// Compute M-step for indicated items
	//Embian Inc. int mstep_items(bool ignore_max_iter = FALSE, int_vector *ilist2 = 0);
	int mstep_items(bool ignore_max_iter = false, int_vector *ilist2 = 0);
	// Retyped ignore_max_iter from "int" to "bool", ww, 2-25-2008.

	// Return maximum relative difference between parameters in current and previous
	// EM iteration computed in last call to mstep_items
	double mstep_max_diff();

	// Return message from optimization procedure generated in
	// last call to mstep_items for an item
	int mstep_message(int itemno);

	// Set maximum number of iterations used by optimization procedure
	// in mstep_items for one item.
	void mstep_max_iter(int itemno, int maxiter);

	// M-step for latent distribution in one group
	// Returns maximum relative difference between new and old probabilities
	double mstep_dist(Estep *e, int group);

	/*** Examinee functions ***/

	// add one examinee
	int add_examinee(int_vector *responses, int group = 1, double count = 1.0);

	// Return vector of all item responses for an examinee
	int_vector *examinee_responses(int examineeno);

	// Returns string containing examinee responses to all items.
	// A response in the first response category is '0'.
	const char *examinee_response_str(int examineeno);

	// Return group examinee belongs to
	int examinee_get_group(int examineeno);

	// Sets group examinee belongs to
	void examinee_set_group(int examineeno, int group);

	// Set count for an examinee
	void examinee_set_count(int examineeno, double count);

	// Return count for an examinee
	double examinee_get_count(int examineeno);

	// Return total examinee counts in an examinee group (1, 2, ...).
	// If group==0 then return total count across all groups.
	double examinees_count(int group = 0);

	// Set posterior distribution for an examinee
	void examinee_set_posterior(int examineeno, double_vector *posterior);

	// Return copy of posterior distribution for an examinee
	double_vector* examinee_get_posterior(int examineeno);

	// Return mean of posterior distribution for an examinee
	double examinee_posterior_mean(int examineeno);

	// Return MLE theta estimate for an examinee
	double examinee_theta_MLE(int examineeno, double minTheta, double maxTheta, double precision = 0.001, int_vector *ilist5 = 0);

	// Set seed of random number generator used for bootstrap samples
	void bootstrap_seed(unsigned long seed);

	// Generate bootstrap sample of examinees
	void bootstrap_sample();

	// Set seed of random number generator used for simulating item responses
	void simulate_seed(unsigned long seed);

	// Simulate item responses corresponding to latent variable value theta
	int_vector *simulate_responses(double theta, int_vector *ilist2 = 0);

	// Return a string containing simulated item responses
	const char *simulate_response_str(double theta, int_vector *ilist2 = 0);

	// Read item responses from a string
	int_vector *get_responses(char *line, int_vector *offset, int_vector *len);

	// Read item responses for a subset of items from a string
	int_vector *get_responses_missing(char *line, int_vector *offset, int_vector *len, int_vector *items);
};

#ifndef SWIG
} // namespace etirm
#endif

#endif
