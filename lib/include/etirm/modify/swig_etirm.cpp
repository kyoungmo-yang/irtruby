/*! \file swig_etirm.cpp
 
  \brief
  Definitions of functions and classes to generate scripting language wrapper 
  functions for ETIRM using SWIG.

  The file swig_etirm_types.h (which is included in swig_etirm.h) must be 
  provided which should contain typedef's for the following types:
 
  ResponseVector	Vector containing item responses.
  item_type		Item class.
  examinee_type	Examinee class.
  lvdist_type		Latent variable distribution class.
  random_type		Type of random number generator used for bootstrap.

  For example:

  // Do not process with SWIG, SWIG cannot handle the template arguments
  xx_ifndef SWIG
 
  // Vector of item responses
  typedef std::vector<etirm::Response> ResponseVector;
 
  // Class to hold information about discrete latent variable distribution
  typedef etirm::DiscreteLatentDist<etirm::Real> lvdist_type;
 
  // Class to hold information about examinees
  typedef etirm::ExamineeGrpCov<ResponseVector, etirm::RealVector> examinee_type;
 
  // Type used for array of items. Using ItemNR allows different items to be modeled
  // by different classes descendent from ItemNR (e.g., the 3PL model could be used for
  // some items, and the 2PL model used for other items).
  typedef etirm::ItemNR<lvdist_type> item_type;
 
  // Use the Mersenne Twister (http://www.math.keio.ac.jp/matumoto/emt.html)
  // from the Boost random number library
  // (http://www.boost.org/libs/random/index.html) 
  // as random number generator for bootstrap samples.
  typedef boost::mt19937 random_type;

  xx_endif // SWIG
 
  An application using these wrapper functions should declare a subclass
  of SwigEtirmRun whose constructor initializes the item object pointers, and 
  assigns them to the 'items' data member of SwigEtirmRun. 
  An application should provide functions for creating and deleting an object
  that is a subclass of SwigEtirmRun which is stored in the global variable
  gEtirmRun defined in this file.

  Definitions of the function CheckRunInit must be provided. This
  function is declared in swig_etirm.h, but each application must define it.
  It is not defined in swig_etirm.cpp.
 
  The function CheckRunInit checks whether gEtirmRun has been
  initialized, and if not throws an exception. The message contained
  in the exception will differ for different applications. 

  Estimation Toolkit for Item Response Models (ETIRM)
  http://www.smallwaters.com/software/cpp/etirm.html

  Author(s): 
  Embian, Inc., maintenance (Site : http://www.embian.com, Email : mo@embian.com, yjj0309@gmail.com)
  Werner Wothke, maintenance (http://www.smallwaters.com)
  Brad Hanson (http://www.b-a-h.com/)
  See the file LICENSE for information on usage and redistribution.

  Copyright (C) 2008, Werner Wothke
  Copyright (c) 2000-2001, Bradley A. Hanson
  Copyright (C) 2010, Embian, Inc.
 */

#ifdef ETIRM_NO_DIR_PREFIX
#include "swig_etirm.h"
#include "MStepIRT.h"
#include "DiscreteNormalDist.h"
#include "ItemParamPriorBeta4.h"
#include "ItemParamPriorLogNormal.h"
#include "ItemParamPriorNormal.h"
#include "BootstrapSample.h"
#include "SimulateResponses.h"
#include "ExamineeThetaMLE.h"
#else
#include "etirm/swig_etirm.h"
#include "etirm/MStepIRT.h"
#include "etirm/DiscreteNormalDist.h"
#include "etirm/ItemParamPriorBeta4.h"
#include "etirm/ItemParamPriorLogNormal.h"
#include "etirm/ItemParamPriorNormal.h"
#include "etirm/BootstrapSample.h"
#include "etirm/SimulateResponses.h"
#include "etirm/ExamineeThetaMLE.h"
#endif

// Use integer random number generator from boost library
// (http://www.boost.org).
#include <boost/random/uniform_int.hpp>

#include <cstdio>  // sprintf prototype
// for compilers which do not put C library functions in std namespace
#ifdef BOOST_NO_STDC_NAMESPACE
namespace std
{ using ::sprintf;}
#endif

// Definitions not wrapped by SWIG
namespace etirm
{

  /*!
    \brief
    Increments counts based on item response.
   
    \section function_args Function Parameters
   
    \param[in]  respIndex "One" offset index corresponding to response to increment 
        count for (1 = first response category, 2 = second response category, etc.).
    \param[in]  count		Amount to increment counts.
    \param[in]  group		Group to increment counts for (1, 2, ...).
   */
  void ItemRespCounts::AddResponse(int respIndex, Real count, int group)
  {

    if (respIndex > catCounts.size() || respIndex < 1)
    {
      throw InvalidArgument("Invalid response index", "ItemRespCounts::AddResponse");
    }

    totalCount[0] += count; // total count
    totalCount[group] += count; // count for group

    catCounts(1, respIndex) += count; // total count
    catCounts(group+1, respIndex) += count; // count for group
  }

  /*!
    \brief
    Class constructor.
   
    \section function_args Function Parameters
   
    \param[in]  nitems		Total number of items.
    \param[in]  nlatentcat	Number of categories of discrete latent variable.
    \param[in]  ngroups		Number of examinee groups.
    \param[in]  minTheta  Minimum value of theta points for discrete latent 
        variable distribution.
    \param[in]  maxTheta	Maximum value of theta points for discrete latent 
        variable distribution.
    \param[in]  uniquePoints	If true then unique latent distribution points 
        are used for each examinee group.
   */
  SwigEtirmRun::SwigEtirmRun(int nitems, int nlatentcat, int ngroups, Real minTheta, Real maxTheta,
      bool uniquePoints) :
    numItems(nitems), numLatentCat(nlatentcat), numGroups(ngroups), examineeCounts(ngroups+1, 0.0),
        items(nitems), itemStats(nitems), minProc(nitems), latentDist(nlatentcat, ngroups,
            uniquePoints), rand_boot(0), base_rand_simulate(0), rand_simulate(0)
  {
    int i;

    // Set default points and weights for latent variable distribution to standard normal
    // for first group
    DiscreteNormalDist(nlatentcat, minTheta, maxTheta, latentDist.begin_points(1),
        latentDist.begin_weights(1));

    // Set points and weights for other groups equal to the points and weights for group 1
    for (i=2; i<=ngroups; ++i)
    {
      DiscreteLatentDist<Real>::weight_iterator w1 = latentDist.begin_weights(1);
      DiscreteLatentDist<Real>::weight_iterator w2 = latentDist.begin_weights(i);
      for (int j = nlatentcat; j--; ++w1, ++w2)
      {
        *w2 = *w1;
      }

      if (uniquePoints)
      {
        // initialize unique points for examinee group
        DiscreteLatentDist<Real>::point_iterator w1p = latentDist.begin_points(1);
        DiscreteLatentDist<Real>::point_iterator w2p = latentDist.begin_points(i);
        for (int j = nlatentcat; j--; ++w1p, ++w2p)
        {
          *w2p = *w1p;
        }
      }
    }

    /* initialize vectors of pointers to null */
    for (i=0; i<nitems; ++i)
    {
      items[i] = 0;
      itemStats[i] = 0;
      minProc[i] = 0;
    }

  }

  //! Releases memory allocated by constructor.
  SwigEtirmRun::~SwigEtirmRun()
  {
    int i;

    for (i=0; i<numItems; i++)
    {
      if (items[i])
      {
        items[i]->DeletePriors();
        delete items[i];
      }

      if (itemStats[i])
        delete itemStats[i];
      if (minProc[i])
        delete minProc[i];
    }

    int n = examinees.size();
    for (i=0; i<n; ++i)
    {
      delete examinees[i];
    }

    delete rand_boot;
    delete base_rand_simulate;
    delete rand_simulate;
  }

  // Check that data for examinees exist
  void SwigEtirmRun::CheckExaminees(const char *funcname)
  {
    if ((this->examinees).size() == 0)
    {
      throw RuntimeError("Error: No examinee data", funcname);
    }
  }

  // Check that examinee number is valid
  void SwigEtirmRun::CheckExamineeNumber(int examno, const char *funcname)
  {
    if ( (examno < 1) || ( examno > (this->examinees).size() ) )
    {
      char errstr[50];
      std::sprintf(errstr, "Invalid examinee number: %d", examno);
      this->returnString = errstr;
      throw RuntimeError((this->returnString).c_str(), funcname);
    }
  }

  // Check that item number is valid
  void SwigEtirmRun::CheckItemNumber(int itemno, const char *funcname)
  {
    if (itemno < 1 || itemno >(this->numItems) ) // Syntax edited, ww, 2-24-2008.
    {
      char errstr[50];
      std::sprintf(errstr, "Invalid item number: %d", itemno);
      this->returnString = errstr;
      throw RuntimeError((this->returnString).c_str(), funcname);
    }
  }

  // Check that group number is valid
  // 	group		Value to check.
  // 	funcname	Name of calling function (used in error message)
  void SwigEtirmRun::CheckGroup(int group, const char *funcname)
  {
    if (group < 1 || group >(this->numGroups) ) // Syntax edited, ww, 2-24-2008
    {
      char errstr[50];
      std::sprintf(errstr, "Invalid examinee group: %d", group);
      this->returnString = errstr;
      throw RuntimeError((this->returnString).c_str(), funcname);
    }
  }

  /*! 
    \brief
    Tests whether item parameter index is valid for the item.
   
    \section function_args Function Parameters
    
    \param[in]  *item		Pointer to item object for which index is checked.
    \param[in]  index		Zero-offset index of item parameter.
    \param[in]  *funcname	pointer to name of calling function (used in error message).
   */
  void CheckItemParam(item_type *item, int index, const char *funcname, SwigEtirmRun *gEtirmRun)
  {
    if (item->NumParameters() <= index || index < 0)
    {
      char errstr[100];
      std::sprintf(errstr, "Invalid item parameter index for item %d: %d", item->Index()+1, index);
      gEtirmRun->returnString = errstr;
      throw RuntimeError((gEtirmRun->returnString).c_str(), funcname);
    }
  }

  /*!
    \brief
    Create a prior distribution object and return a pointer to it.
   
    \section function_args Function Parameters

    \param[in]  pstr  Type of prior distribution, valued as ("none", "normal", "lognormal", "beta").
    \param[in]  priorparam	Parameters of prior distribution.
    \param[in]  funcname	Name of calling function (used in error messages)
   */
  ItemParamPrior *CreatePrior(const std::string &pstr, const double_vector &priorparam,
      const char *funcname)
  {
    const char *err_message = "Invalid number of prior parameters";

    if (pstr.compare("none") == 0)
    {
      return 0;
    }
    else if (pstr.compare("normal") == 0)
    {
      if (priorparam.size() != 2)
        throw RuntimeError(err_message, funcname);
      return new ItemParamPriorNormal(priorparam[0], priorparam[1]);
    }
    else if (pstr.compare("lognormal") == 0)
    {
      if (priorparam.size() != 2)
        throw RuntimeError(err_message, funcname);
      return new ItemParamPriorLogNormal(priorparam[0], priorparam[1]);
    }
    else if (pstr.compare("beta") == 0)
    {
      if (priorparam.size() != 4)
        throw RuntimeError(err_message, funcname);
      return new ItemParamPriorBeta4(priorparam[0], priorparam[1], priorparam[2], priorparam[3]);
    }
    else
    {
      throw RuntimeError("Invalid prior type", funcname);
    }
    return 0;
  }

  /*!
    \brief
    Converts a response for an item into a character, where '0' represents the 
    first response, '1' represents the second response, etc.
   
    If the number of response categories for the item is greater than 10, then the 
    characters returned will be ascii characters greater than '9'. For example, a 
    ':' is return for a response in response category 11, since ':' is one greater
    than '9' in the ascii sequence.

    \section function_args Function Parameters
    
    \param[in]  r Item response to be converted.
    \param[in]  *item Pointer to item object.
   */
  char Resp2Char(Response r, const item_type *item)
  {
    const char zero = '0';
    Response np = item_type::NotPresentedResponse();
    char cr;
    if (r == np)
      cr = np - item->FirstResponse() + zero;
    else
      cr = item->ResponseIndex(r) + zero;

    return cr;
  }

  /****** Functions for which SWIG wrappers are generated ******/

  /* Member functions for estep class */

  /*!
    \brief
    Initialize estep object.
   
    \section function_args Function Parameters
   
    \param[in]  *itemno	Pointer to list of item numbers to use for computing 
    examinee posteriors distribution in E-Step. If itemno = NUL, use all items (default).
   */
  Estep::Estep(SwigEtirmRun *runner, int_vector *itemno)
  {
    // CheckRunInit();
	gEtirmRun = runner;
	if (!gEtirmRun) throw RuntimeError("The new_items_dist command has not been executed", 0);
    if (itemno)
    {
      mItems = new ItemVector(ItemSubset(itemno,gEtirmRun->items,"new_estep"));

      mEStep = new estep_type(mItems->begin(), mItems->end(), gEtirmRun->latentDist);
    }
    else
    {
      mItems = 0;
      mEStep = new estep_type((gEtirmRun->items).begin(), (gEtirmRun->items).end(), gEtirmRun->latentDist);
    }
  }

  /*!
    \brief
    Release memory allocated in constructor.
   */
  Estep::~Estep()
  {
    if (mItems) delete mItems;
    if (mEStep) delete mEStep;
	gEtirmRun = 0;
  }

  /*!
    \brief
    Calls mEStep->DoEStep to perform E-step.
   
    Returns marginal loglikelihood of examinees' responses (sum over examinees of the 
    marginal loglikelihood of an examinee's responses) plus sum of prior likelihoods 
    over all item parameters.  This is the value of the marginal posterior density that 
    the EM algorithm is maximizing at the values of the item parameters computed in the 
    last M-step.
   
    \section function_args Function Parameters
   
    \param[in]  compute_post  Flag: If compute_post == TRUE, then compute posterior 
    distributions for all examinees. If compute_post == FALSE, then use previously 
    stored posteriors for examinees. (Default: TRUE).
   
    \param[in]  store_post  Flag: If store_post == TRUE, then store the posterior 
    distribution computed for each examinee. These posterior distributions are 
    stored as part of the examinee objects.
    (Default: FALSE).
   
    \param[in]  *estep_items  Pointer to list of items for which n and r are updated. If a 
    null pointer is passed then n arnd r are updated for all items used in the E-step 
    to compute examinee posterior distributions.
   */
  double Estep::compute(bool compute_post, bool store_post, int_vector *estep_items) // compute_post and store_post retyped as "bool", ww, 2-24-2008.
  {
    const char *fname = "estep_compute";
    // CheckRunInit();
	if (!gEtirmRun) throw RuntimeError("The new_items_dist command has not been executed", 0);
    gEtirmRun->CheckExaminees(fname);

    if (estep_items)
    {
      if (estep_items->size() > 0)
      {
        ItemVector itemsub = ItemSubset(estep_items, gEtirmRun->items, fname);
        return mEStep->DoEStep((gEtirmRun->examinees).begin(), (gEtirmRun->examinees).end(), itemsub.begin(), itemsub.end(), compute_post, store_post);
      }
      else // do not update n and r for any items
      {
        return mEStep->DoEStep((gEtirmRun->examinees).begin(), (gEtirmRun->examinees).end(), (gEtirmRun->items).begin(), (gEtirmRun->items).begin(), compute_post, store_post);
      }
    }
    else
    {
      return mEStep->DoEStep((gEtirmRun->examinees).begin(), (gEtirmRun->examinees).end(), compute_post, store_post);
    }
  }

} // namespace etirm

