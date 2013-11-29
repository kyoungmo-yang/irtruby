/*! \file swig_etirm.h
  
  \brief
  Declarations of functions and classes to generate scripting
  language wrapper functions for ETIRM using SWIG.

  See swig_etirm.cpp for documentation.

  Estimation Toolkit for Item Response Models (ETIRM)
  http://www.smallwaters.com/software/cpp/etirm.html

  Author(s): 
  Embian, Inc., maintenance (Site : http://www.embian.com, Email : mo@embian.com, yjj0309@gmail.com)
  Werner Wothke, maintenance (http://www.smallwaters.com)
  Brad Hanson (http://www.b-a-h.com/)
  See the file LICENSE for information on usage and redistribution.

  Copyright (C) 2008, Werner Wothke
  Copyright (c) 2000-2001, Bradley A. Hanson
  Copyright (C) 2010, Embian Inc.
 */

#ifndef SWIGETIRM_H_
#define SWIGETIRM_H_

// This must be provided - declares types unique to application 
// (see comments at the top of swig_etirm.cpp).
// It is not provided in ETIRM.
#include "swig_etirm_types.h"

#ifdef ETIRM_NO_DIR_PREFIX
#include "EStepDiscrete.h"
#else
#include "etirm/EStepDiscrete.h"
#endif

#include "Uncmin.h"

#include <vector>
#include <string>

// Use uniform random number generator from boost library
// (http://www.boost.org).
#include <boost/random/uniform_01.hpp>

#ifndef SWIG
namespace etirm
{
#endif

  // Declarations not wrapped with SWIG
#ifndef SWIG
  //! Vector type for double precision values.
  typedef RealVector double_vector;

  //! Vector type for integer values.
  typedef SCPPNT::Vector<int> int_vector;

  //! Vector type for item pointers
  typedef std::vector<item_type *> ItemVector;

  //! Class used in E-step calculation
  typedef EStepDiscrete<examinee_type, item_type, ItemVector::iterator, lvdist_type> estep_type;

  //! Container to hold Uncmin optimization objects
  typedef std::vector<Uncmin<RealVector, RealMatrix, item_type> *> UncminVector;

  //! Container holding examinee object pointers
  typedef std::vector<examinee_type *> ExamineeVector;

  /*!
    \brief
    Class to hold counts of number of examinees who responded in each response 
    category of an item.
   */
  class ItemRespCounts
  {
public:

    /*! 
      \brief
      Class constructor.
      
      \section function_args Function Parameters
   
      \param[in]  ncat Number of response categories for the item.
      \param[in]  ngroup  Number of examinee groups.
     */
    ItemRespCounts(int ncat, int ngroup) :
      catCounts(ngroup+1, ncat, 0.0), totalCount(ngroup+1, 0.0)
    {
    }

    //! Deconstructor
    virtual ~ItemRespCounts()
    {
    }

    // Increment counts for a particular item response (1-offset)
    void AddResponse(int respIndex, Real count, int group);


    /*!
      \brief
      Returns total number of responses to item.
      
      \section function_args Function Parameters
   
      \param[in]  group Group(1, 2, ...) for which counts are to be returned.
      
      Note: group = 0 returns counts across all groups.
     */
    Real RespCount(int group)
    {
      return totalCount[group];
    }

    /*!
      \brief
      Returns number of responses in each response category.
      
      \section function_args Function Parameters
   
      \param[in]  group Group(1, 2, ...) for which counts are to be returned.
      
      Note: group = 0 returns counts across all groups.
     */      
    RealVector CategoryCounts(int group)
    {
      return RealVector(catCounts.begin_row(group+1), catCounts.begin_row(group+1)
          +catCounts.num_columns());
    } // ww: Edit "num_cols" -> "num_columns", 1/8/2008.


protected:
    /*!
      \brief
      Total number of examinees who have responses for the item.
      
      totalCount[0] is total count across all groups
      totalCount[i], i>0, is count for group i
     */
    RealVector totalCount;

    /*!
      \brief
      Number of responses in each response category.
      
      Columns give counts for response categories
      Row 1 gives counts across all groups.
      Row i, i>1,  gives counts for group i-1.
     */
    RealMatrix catCounts;
  };

  //! Container holding pointers to objects containing item statistics for each item
  typedef std::vector<ItemRespCounts *> ItemStatsVector;

  //! Class holding information used for one modeling problem.
  class SwigEtirmRun
  {

public:

    SwigEtirmRun(int nitems, int nlatentcat, int ngroups, Real minTheta = -4.0,
        Real maxTheta = 4.0, bool uniquePoints = false);

    virtual ~SwigEtirmRun();

    ExamineeVector examinees;
    //!< Vector holding pointers to examinee objects.
    
    ItemVector items;
    //!< Vector holding pointers to item objects.
    
    ItemStatsVector itemStats;
    //!< Vector holding pointers to item count objects.
    
    UncminVector minProc;
    //!< Vector holding pointers to minimization objects.
    
    lvdist_type latentDist;
    //!< Latent distribution object.

    random_type *rand_boot;
    //!< Base random number generator object for bootstrap.

    random_type *base_rand_simulate;
    //!< Base random number generator object for simulating item responses.

    boost::uniform_01<random_type> *rand_simulate;
    //!< Uniform random number generator object for simulating item responses.

    std::string returnString;
    //!< Temporary space to hold strings returned by functions.

    int numItems;
    //!< Number of operational/field test items to model.
    
    int numLatentCat;
    //!< Number of discrete categories to approximate the latent ability distribution.
    
    int numGroups;
    //!< Number of non-equivalents groups of examinees.
    
    RealVector examineeCounts;
    //!< Examinee counts in each group: examineeCounts[0] contains count across all groups; examineeCounts[i], i>0, contains count for group i.

    double mstepMaxDiff;
    //!< Maximum relative difference between item parameters computed in last M-Step call.

    /*!
      \brief
      Checks that data for examinees exist.
      
      \section function_args Function Parameters
   
      \param[in] *funcname  Pointer to name of calling routine.
     */
    void CheckExaminees(const char *funcname);

    /*!
      \brief
      Check that the examinee number is within its valid range.
      
      \section function_args Function Parameters
   
      \param[in] examno Examinee number.
      \param[in] *funcname  Pointer to name of calling routine.
     */
    void CheckExamineeNumber(int examno, const char *funcname);

    /*!
      \brief
      Checks for valid item number.
      
      \section function_args Function Parameters
   
      \param[in] itemno Item number.
      \param[in] *funcname  Pointer to name of calling routine.
     */
    void CheckItemNumber(int itemno, const char *funcname);

    /*!
      \brief
      Checks for valid group.
      
      \section function_args Function Parameters
      
      \param[in]  group Group number.
      \param[in] *funcname  Pointer to name of calling routine.
     */
    void CheckGroup(int group, const char *funcname);
  };





  /*!
    \brief
    Check that new_etirm has been called
  
    Note: This function is not defined in swig_etirm.cpp, it must be defined 
    somewhere in the application.
   */
  void CheckRunInit();

  /*! 
    \brief
    Tests whether item parameter index is valid for the item.
   
    \section function_args Function Parameters
    
    \param[in]  *item   Pointer to item object for which index is checked.
    \param[in]  index   Zero-offset index of item parameter.
    \param[in]  *funcname pointer to name of calling function (used in error message).
   */
  void CheckItemParam(item_type *item, int index, const char *funcname, SwigEtirmRun *gEtirmRun);

  //Create a prior distribution object and return a pointer to it.
  ItemParamPrior *CreatePrior(const std::string &pstr, const double_vector &priorparam,
      const char *funcname);

  /*!
    \brief 
    Creates a container holding the elements corresponding to a subset of items.

    ItemSubset is based on a container holding elements corresponding to all items. 
    This container is typically used to create subsets of etirmrun->items (which 
    is an ItemVector) and etirmrun->minProc (which is a UncminVector). 
    
    \section template_args Template Parameters
   
    \param V ItemVector type.
   
    \section function_args Function Parameters
   
    \param[in]  *itemnums		Pointer to vector containing numbers of items to include in subset.
    \param[in]  &all  Address of container containing all items.
    \param[in]  *funcname		Pointer to name of calling function (used in error messages)
   */
  template <class V> V ItemSubset(int_vector *itemnums, V &all, const char *funcname)
  {
    int i;
    int nsub = itemnums->size();
    int nall = all.size();

    if (nsub < 1 || nsub> nall)
    {
      throw RuntimeError("Invalid number of items", funcname);
    }

    V sub(nsub);
    int *flag = new int[nall]; // flags indicating which items are in subset
    for (i = 0; i < nall; ++i) flag[i] = 0;

    typename V::iterator iall = all.begin();
    typename V::iterator isub = sub.begin();
    int_vector::iterator ino = itemnums->begin();
    for(i = 0; i < nsub; ++i, ++isub, ++ino)
    {
      int index = *ino - 1; // elements of itemnums are 1-offset indices
      if (index < 0 || index >= nall)
      {
        delete [] flag;
        throw RuntimeError("Invalid item number", funcname);
      }
      if (flag[index] == 1)
      {
        delete [] flag;
        throw RuntimeError("Duplicate item number", funcname);
      }
      flag[index] = 1;
      *isub = iall[index];
    }
    delete [] flag;

    return sub;
  }

  // Convert an item response to a character, where '0'
  // corresponds to the first response category
  char Resp2Char(Response r, const item_type *item);

#endif // SWIG
  /*** Classes and functions to be wrapped by SWIG ***/

  //! SWIG wrapper for EStepDiscrete object
  class Estep
  {

  public:

    Estep(SwigEtirmRun *runner, int_vector *ilist1 = 0);
    ~Estep();

    //Embian Inc. double compute(bool compute_prior = TRUE, bool store_prior = FALSE, int_vector *ilist4 = 0);
    double compute(bool compute_prior = true, bool store_prior = false, int_vector *ilist4 = 0);
    // calls mEStep->DoEStep  // Retyped compute_prior and store_prior as "bool", ww, 2-24-2008.

#ifndef SWIG
    /*! 
      \brief
      Returns EStepDiscrete object. 
      
      Do not create SWIG wrapper for this member function!
     */
    estep_type *GetEStep()
    { return mEStep;}
#endif

  private:

    estep_type *mEStep; //!< EStepDiscrete object.
    ItemVector *mItems; //!< List of items used to initialize mEStep
	SwigEtirmRun *gEtirmRun;
  };

#ifndef SWIG
} // namespace etirm
#endif // SWIG

#endif
