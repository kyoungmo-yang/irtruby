/*! \file EStepDiscrete.h
 
  \brief 
  Class definitions of EStepDiscrete and NMatrixVec.
  
  Used to perform E-step calculation and store results for IRT model with a 
  discrete latent variable distribution. For use with dichotomous items where 
  n and r are stored for each latent variable category.
  
  Note that the main definitions of the ExamineePosterior and DoEStep member 
  templates of the EStepDiscrete class are outside of the class definition, 
  although there are duplicate definitions given inside the class definition 
  for use by compilers which do not allow member templates to defined outside
  the class definition. 
  
  If the symbol BOOST_MSVC6_MEMBER_TEMPLATES is defined the member template 
  definitions inside the class definition are used, otherwise the member 
  template definitions outside the class definition are used. If the 
  ExamineePosterior and DoEStep member templates are modified, the 
  modifications must be made in both definitions.
 
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

#ifndef ETIRM_ESTEPDISCRETE_H_
#define ETIRM_ESTEPDISCRETE_H_

#ifdef ETIRM_NO_DIR_PREFIX
#include "etirmtypes.h"
#include "ItemParamPrior.h"
#else
#include "etirm/etirmtypes.h"
#include "etirm/ItemParamPrior.h"
#endif

#include <cmath> // for exp and log
#include <vector>

// for compilers which do not put C library functions in std namespace
#ifdef BOOST_NO_STDC_NAMESPACE
namespace std
{ using ::exp; using ::log;}
#endif

namespace etirm
{
  //! Value used to represent the log of zero probabilities
  const double logZero = -1021.0;

  /*!
    \brief
    Helper class template used by EStepDiscrete for managing a vector 
    containing pointers to RealMatrix objects.
 
    \section template_args Template Parameters
   
    \param II  Iterator over item objects.
   */
  template <class II> class NMatrixVec
  {

  public:

    //! Type of iterator over item matrices.
    typedef std::vector<RealMatrix *>::iterator iterator;

    NMatrixVec(II bitem, II eitem, int nlatentcat);
    ~NMatrixVec();

    iterator begin()
    {
      return mVector.begin();
    }
    //!< Iterator to matrices over items.

    RealMatrix *operator[](int Index)
    {
      return mVector[Index];
    }
    //!< Return pointer to matrix for item Index+1 (Index is zero-offset).

  private:

    std::vector<RealMatrix *> mVector;
    /*!< Each matrix corresponds to an item.
      The rows of each matrix contain a discrete latent variable distribution
      corresponding to one response.
     */
  };

  /*!
    \brief
    Constructor 
    
    \section template_args Template Parameters
    
    \param II  Iterator over item objects.
    
    \section function_args Function Parameters
    
    \param[in]  bitem Iterator to first item.
    \param[in]  eitem Iterator to last item.
    \param[in]  nlatentcat  Number of discrete categories of latent variable distribution.
   
   */
  template <class II> NMatrixVec<II>::NMatrixVec(II bitem, II eitem, int nlatentcat) :
  mVector(eitem-bitem)
  {
    iterator ri = mVector.begin();
    for (II i = bitem; i != eitem; ++ri, ++i)
    {
      *ri = new RealMatrix((*i)->NumRespCat(), nlatentcat);
    }
  }

  /*! Destructor */
  template <class II> NMatrixVec<II>::~NMatrixVec()
  {
    iterator ri = mVector.begin();
    for (int i = mVector.size(); i--; ++ri)
    {
      if (*ri)
      delete *ri;
    }
  }

  /*!
    \brief
    Class template to perform E-step calculation and store results
    for IRT model with a discrete latent variable distribution.
 
    \section template_args Template Parameters
    
    \param E Examinee type.
    \param I Item type. 
    \param II  Iterator over item objects.
    \param D  Class for discrete latent variable distribution.
   */
  template <class E, class I, class II, class D> class EStepDiscrete
  {

  public:

    typedef RealMatrix::row_iterator ngroup_iterator;
    //!< Type of iterator over marginal group probabilities.

    typedef typename D::point_iterator point_iterator;
    //!< Type of iterator over points of latent variable distribution.

    EStepDiscrete(II bitem, II eitem, D &dist);
    //!< Constructor.

    ~EStepDiscrete();
    // Destructor.

#ifndef BOOST_MSVC6_MEMBER_TEMPLATES
    /*!
      \brief
      Computes the E-step of the EM algorithm for IRT models with a
      discrete latent variable distribution.
      
      Returns marginal loglikelihood of examinees' responses
      (sum over examinees of the marginal loglikelihood of an examinee's responses)
      plus sum of prior likelihoods over all item parameters.
      This is the value of the marginal posterior density that the
      EM algorithm is maximizing at the values of the item parameters
      computed in the last M-step. The log of the priors for the item
      parameters are added to this value for the items for which
      n and r are calculated.
      
      Results of the E-step are stored in data member nGroups. The posterior
      distributions for examinees are stored in the examinee objects
      if storeExamineePosterior is true, and updated n's and r's are computed
      the items given by itemsNR_begin and itemsNR_end.
      
      A duplicate of this definition is given inside the class definition
      for use when BOOST_MSVC6_MEMBER_TEMPLATES is defined.
 
      \section template_args Template Parameters
      
      \param D  Class for discrete latent variable distribution.
      \param E  Examinee type.
      \param EI Iterator over pointers to examinee objects.
      \param I  Item type. 
      \param II Iterator over item objects.
 
      \section function_args Function Parameters
      
      \param[in]  examinees_begin Iterator to pointer to first examinee
      \param[in]  examinees_end Iterator to pointer to one past last examinee
      \param[in]  itemsNR_begin Iterator to first item pointer for which n and r will be updated.
          The items for which n and r are updated can be different from the items
          used to compute the posterior distribution of the latent variable for
          each examinee. If itemsNR_end - itemsNR_begin == 0 then n and r are not
          updated for any items.
      \param[in]  itemsNR_end Iterator to one past last item pointer for which n and r will be 
          updated.
      \param[in]  computeExamineePosterior If TRUE posterior latent variable distribution
          for each examinee is computed. If FALSE previously stored posterior
          latent variable distribution for each examinee is used.
      \param[in]  storeExamineePosterior If TRUE posterior latent variable distribution
          is stored for each examinee. If 'computeExamineePosterior' is
          FALSE the value of the argument is not used (in that case a
          previously stored posterior distribution for each examinee is
          being used in this function).
     */
    template <class EI>
    Real DoEStep(EI examinees_begin, EI examinees_end, II itemsNR_begin, II itemsNR_end,
        bool computeExamineePosterior, bool storeExamineePosterior);

    /*!
      \brief
      Version of DoEStep in which the items that are used to compute the
      posterior distributions for examinees are also the items for which
      n and r are updated, i.e., itemsNR_begin == items_begin and
      itemsNR_end == items_end.
 
      A duplicate of this definition is given inside the class definition
      for use when BOOST_MSVC6_MEMBER_TEMPLATES is defined.
 
      \section template_args Template Parameters
      
      \param D  Class for discrete latent variable distribution.
      \param E  Examinee type.
      \param EI Iterator over pointers to examinee objects.
      \param I  Item type. 
      \param II Iterator over item objects.
 
      \section function_args Function Parameters
      
      \param[in]  examinees_begin Iterator to pointer to first examinee
      \param[in]  examinees_end Iterator to pointer to one past last examinee
      \param[in]  computeExamineePosterior If TRUE posterior latent variable distribution
          for each examinee is computed. If FALSE previously stored posterior
          latent variable distribution for each examinee is used.
      \param[in]  storeExamineePosterior If TRUE posterior latent variable distribution
          is stored for each examinee. If 'computeExamineePosterior' is
          FALSE the value of the argument is not used (in that case a
          previously stored posterior distribution for each examinee is
          being used in this function).
     */
    template <class EI>
    Real DoEStep(EI examinees_begin, EI examinees_end, bool computeExamineePosterior,
        bool storeExamineePosterior);

    /*!
     \brief
     Computes posterior distribution of discrete latent variable for an examinee.
     Returns marginal likelihood of the examinee's responses.
     
     A duplicate of this definition is given inside the class definition
     for use when BOOST_MSVC6_MEMBER_TEMPLATES is defined.

     \section template_args Template Parameters
     
     \param PI Iterator over posterior probabilities.

     \section function_args Function Parameters
     
     \param[in]  examinee Object holding information about an examinee's item responses
     \param[in]  begin_posterior Iterator pointing to first element of container holding posterior
         probabilities.
     \param[in]  end_posterior Iterator pointing to one past last element of container holding
         posterior probabilities. This argument is only used to make sure there
         is enough space in the container which will hold the posterior probabilities.
     */
    template <class PI>
    Real ExamineePosterior(E &examinee, PI begin_posterior, PI end_posterior);
    // Calculates posterior latent variable distribution for an examinee.
#endif

    void CalcResponseProb();
    //!< Calculates matrices of response probabilities used in DoEStep.

    ngroup_iterator GetNGroup(int group)
    {
      return (*nGroups).begin_row(group);
    }
    //!< Returns iterator to latent variable distribution for a group.

    int size()
    {
      return numLatentVarCat;
    }
    //!< Returns number of categories in discrete latent variable distribution.

    point_iterator GetPoints()
    {
      return latentvar_dist.begin_points();
    }
    //!< Returns Iterator to points of latent variable distribution.

  protected:

    RealMatrix *nGroups;
    /*!< 
      Row i gives the expected number of examinees at each category of the latent variable
      for group i computed in DoEStep.
     */

    std::vector<NMatrixVec<II> *> mRespProb;
    /*!<	
      Element i is a pointer to a NMatrixVec object for examinee group i.
      Different NMatrixVec objects are needed for different examinee
      groups if unique latent distribution points are used for each group.
      The NMatrixVec object for each examinee group holds the log 
      of the response probabilities for each response of each item in each
      latent variable category as computed in CalcResponseProb.
      The k-th element is a pointer to a matrix containing log response probabilities
      for item k. Row i of the k-th matrix contains the log response probabilities
      for response i of item k over the latent variable categories.
     */

    II items_begin;
    //!< Iterator to first item.

    II items_end;
    //!< Iterator to one past last item.

    D &latentvar_dist;
    //!< Object containing discrete latent variable distribution.

    RealMatrix logLatentProb;
    /*!< 
      Stores logarithms of probilities in latentvar_dist for
      for each examinee group for use in ExamineePosterior().
      Row i contains log of latent probabilities for examinee
      group i.
     */

    int numItems; //!< number of items.

    int numLatentVarCat; //!< number of categories in discrete latent variable distribution.

    int numGroupUnique; //!< Number of groups with unique latent distribution points.

    std::vector<int> *itemIndices;
    /*!< This vector contains zero-offset indices in the examinee response
      vector corresponding to the sequence of items given by items_begin and items_end.
     */

    Response notPresentedResponse;
    /*!< Response indicating an item was not presented to an examinee.
      Assumed to be the same for all items.
     */

    /***************
     Define member templates in class declaration for Visual C++ 6
     ****************/
#ifdef BOOST_MSVC6_MEMBER_TEMPLATES
  public:

    /*! 
     \brief
      Returns the posterior for an examinee.
     
      Also update n and r for items the examinee responded to.
     
     \section template_args Template Parameters
     
     \param E  Examinee type.
     \param PI Iterator over posterior distribution vector.
     
     \section function_args Function Parameters
     
     \param[in]  &examinee Address of examinee object.
     \param[in]  begin_posterior Iterator pointing to first element of posterior vector.
     \param[in]  end_posterior Iterator pointing to last element of posterior vector.
     */
    template <class PI> Real ExamineePosterior(E &examinee, PI begin_posterior, PI end_posterior)
    {

      int i, il;

      int group = examinee.Group();

      if ((end_posterior - begin_posterior) != numLatentVarCat)
      {
        throw InvalidArgument("Incorrect size of vector to hold posterior probabilities",
            "EStepDiscrete::ExamineePosterior");
      }

      /* Constants used for loop unrolling */
      int Ndiv4 = numLatentVarCat / 4;
      int Nmod4 = numLatentVarCat - Ndiv4*4;

      /* initialize posterior probabilities */
      RealVector::iterator ipost = begin_posterior;
      RealMatrix::row_iterator iwt = logLatentProb.begin_row(group);
      for (i = Ndiv4; i--; ipost+=4, iwt+=4)
      {
        *ipost = *iwt;
        ipost[1] = iwt[1];
        ipost[2] = iwt[2];
        ipost[3] = iwt[3];
      }
      for (i = Nmod4; i--; ++ipost, ++iwt)
      {
        *ipost = *iwt;
      }
	
      // Embian Inc. added "typename" - 2010
      typename NMatrixVec<II>::iterator item = (numGroupUnique == 1) ? mRespProb[0]->begin()
      : mRespProb[group-1]->begin();
      II iitem = items_begin;
      typename E::response_iterator presp = examinee.responses_begin();
      std::vector<int>::iterator ii = itemIndices->begin();
      for (i = numItems; i--; ++item, ++iitem, ++ii)
      {
        Response resp = presp[*ii];
        if (resp != notPresentedResponse)
        {
          ipost = begin_posterior;
          int index = (*iitem)->ResponseIndex(resp);
          RealMatrix::row_iterator ir = (*item)->begin_row(index+1);
          for (il=Ndiv4; il--; ipost+=4, ir+=4)
          {
            *ipost += *ir;
            ipost[1] += ir[1];
            ipost[2] += ir[2];
            ipost[3] += ir[3];
          }
          for (il=Nmod4; il--; ++ipost, ++ir)
          {
            *ipost += *ir;
          }
        }
      }

      /* find sum in order to standardize posterior */
      ipost = begin_posterior;
      iwt = logLatentProb.begin_row(group);
      Real sum = 0.0;
      for (i = numLatentVarCat; i--; ++ipost, ++iwt)
      {
        if (*iwt != logZero)
        sum += std::exp(*ipost);
        else
        *ipost = logZero;
      }

      /* standardize */
      ipost = begin_posterior;
      Real logsum = std::log(sum);
      for (i = numLatentVarCat; i--; ++ipost)
      {
        if (*ipost != logZero)
        {
          *ipost -= logsum;
          *ipost = std::exp(*ipost);
        }
        else
        *ipost = 0.0;
      }

      return sum;
    }

    /*! 
      \brief
      Computes the E-step of the EM algorithm for IRT models with a
      discrete latent variable distribution.
      
      Returns marginal loglikelihood of examinees' responses
      (sum over examinees of the marginal loglikelihood of an examinee's responses)
      plus sum of prior likelihoods over all item parameters.
      This is the value of the marginal posterior density that the
      EM algorithm is maximizing at the values of the item parameters
      computed in the last M-step. The log of the priors for the item
      parameters are added to this value for the items for which
      n and r are calculated.
      
      Results of the E-step are stored in data member nGroups. The posterior
      distributions for examinees are stored in the examinee objects
      if storeExamineePosterior is true, and updated n's and r's are computed
      the items given by itemsNR_begin and itemsNR_end.
      
      A duplicate of this definition is given outside the class definition
      for use when BOOST_MSVC6_MEMBER_TEMPLATES is not defined.

      \section template_args Template Parameters
      
      \param EI Iterator over pointers to examinee objects.
      \param II Iterator over item objects.
 
      \section function_args Function Parameters
      
      \param[in]  examinees_begin Iterator to pointer to first examinee
      \param[in]  examinees_end Iterator to pointer to one past last examinee
      \param[in]  itemsNR_begin Iterator to first item pointer for which n and r will be updated.
           The items for which n and r are updated can be different from the items
           used to compute the posterior distribution of the latent variable for
           each examinee. If itemsNR_end - itemsNR_begin == 0 then n and r are not
           updated for any items.
      \param[in]  itemsNR_end Iterator to one past last item pointer for which n and r will be 
           updated.
      \param[in]  computeExamineePosterior If TRUE posterior latent variable distribution
           for each examinee is computed. If FALSE previously stored posterior
           latent variable distribution for each examinee is used.
      \param[in]  storeExamineePosterior If TRUE posterior latent variable distribution
           is stored for each examinee. If 'computeExamineePosterior' is
           FALSE the value of the argument is not used (in that case a
           previously stored posterior distribution for each examinee is
           being used in this function).     
     */
    template <class EI> Real DoEStep(EI examinees_begin, EI examinees_end, II itemsNR_begin,
        II itemsNR_end, bool computeExamineePosterior, bool storeExamineePosterior)
    {

      int i, j;
      Real loglikelihood = 0.0;

      RealVector posterior(numLatentVarCat);

      *nGroups = 0.0;
      int numItemsNR = itemsNR_end - itemsNR_begin;
      II iitem = itemsNR_begin;
      for (i=numItemsNR; i--; ++iitem)
      {
        /*
          Check that n and r for each item use the same number of latent variable 
          categories as in latent variable distribution used to compute
          examinee posterior distributions (does not check that points match,
          just that number of points match, it is assumed that points also match).
         */
        j = (*iitem)->NumLatentVarCat();
        if (j != numLatentVarCat)
        {
          throw RuntimeError("Mismatch in number of latent variable categories",
              "EStepDiscrete::DoEStep");
        }

        /* Initialize n and r for item to zero */
        (*iitem)->InitializeNR();
      }

      if (computeExamineePosterior)
      {
        // Check that number of latent variable categories for latentvar_dist
        // has not changed.
        if (numLatentVarCat != latentvar_dist.size())
        {
          throw RuntimeError("Number of latent variable categories has changed",
              "EStepDiscrete::DoEStep");
        }

        // Compute log probabilities of each response to each item using
        // current item parameter estimates
        CalcResponseProb();

        /* Set logLatentProb to log of current probabilities of latentvar_dist */
        for (i = 1; i <= latentvar_dist.NumGroups(); ++i)
        {
          RealMatrix::row_iterator ip = logLatentProb.begin_row(i);
          typename D::weight_iterator iwt = latentvar_dist.begin_weights(i);
          for (j = numLatentVarCat; j--; ++ip, ++iwt)
          {
            *ip = (*iwt != 0.0) ? std::log(*iwt) : logZero;
          }
        }
      }

      /* For each examinee compute posterior distribution and
         update n and r for items the examinee responded to.
       */
      for (EI examinee_i = examinees_begin; examinee_i != examinees_end; ++examinee_i)
      {

        Real marginalLikelihood;
        if (computeExamineePosterior) // Compute posterior distribution for examinee.

        {
          marginalLikelihood = ExamineePosterior(**examinee_i, posterior.begin(), posterior.end());

          if (storeExamineePosterior)
          {
            typename E::posterior_vector epost(numLatentVarCat);
            typename E::posterior_vector::iterator iep = epost.begin();
            RealVector::iterator ip = posterior.begin();
            for (i = numLatentVarCat; i--; ++iep, ++ip)
            *iep = *ip;
            (*examinee_i)->SetPosterior(epost);

            (*examinee_i)->SetMarginalRespLikelihood(marginalLikelihood);
          }
        }
        else // use examinee posterior distribution already computed

        {
          typename E::posterior_vector::iterator iep = (*examinee_i)->posterior_begin();
          RealVector::iterator ip = posterior.begin();
          for (i = numLatentVarCat; i--; ++iep, ++ip)
          *ip = *iep;

          marginalLikelihood = (*examinee_i)->GetMarginalRespLikelihood();
        }

        /* update marginal loglikelihood */
        loglikelihood += std::log(marginalLikelihood);

        typename E::response_iterator iresp = (*examinee_i)->responses_begin();
        Real casewt = (*examinee_i)->Count();
        int group = (*examinee_i)->Group();
        iitem = itemsNR_begin;
        for (i = numItemsNR; i--; ++iitem)
        {
          /* Update n and r for each item */
          Response resp = iresp[(*iitem)->Index()];
          if (resp != notPresentedResponse)
          {
            typename I::r_iterator ir = (*iitem)->RVector(resp, group);
            typename I::n_iterator in = (*iitem)->NVector(group);
            RealVector::iterator ipost = posterior.begin();
            for (j = numLatentVarCat; j--; ++ir, ++in, ++ipost)
            {
              *ir += *ipost * casewt;
              *in += *ipost * casewt;
            }
          }

        }

        /* Update marginal distribution for group examinee belongs to */
        RealVector::iterator ipost = posterior.begin();
        RealMatrix::row_iterator igroup = nGroups->begin_row(group);
        for (j = numLatentVarCat; j--; ++ipost, ++igroup)
        {
          *igroup += *ipost * casewt;
        }

      }

      /* Add log of prior densities of item parameter estimates to loglikelihood */
      for (II ii = itemsNR_begin; ii != itemsNR_end; ++ii)
      {
        PriorVector::iterator iprior = (*ii)->PriorsIterator();
        RealVector::iterator iparam = (*ii)->ParametersIterator();
        for (i = (*ii)->NumParameters(); i--; ++iprior, ++iparam)
        {
          if (*iprior)
          loglikelihood += (*iprior)->LogDensity(*iparam);
        }
      }

      return loglikelihood;

    }

    /*!
      \brief
      Version of DoEStep in which the items that are used to compute the
      posterior distributions for examinees are also the items for which
      n and r are updated, i.e., itemsNR_begin == items_begin and
      itemsNR_end == items_end.
 
      A duplicate of this definition is given outside the class definition
      for use when BOOST_MSVC6_MEMBER_TEMPLATES is not defined.
 
      \section template_args Template Parameters
      
      \param EI Iterator over pointers to examinee objects.
 
      \section function_args Function Parameters
      
      \param[in]  examinees_begin Iterator to pointer to first examinee
      \param[in]  examinees_end Iterator to pointer to one past last examinee
      \param[in]  computeExamineePosterior If TRUE posterior latent variable distribution
          for each examinee is computed. If FALSE previously stored posterior
          latent variable distribution for each examinee is used.
      \param[in]  storeExamineePosterior If TRUE posterior latent variable distribution
          is stored for each examinee. If 'computeExamineePosterior' is
          FALSE the value of the argument is not used (in that case a
          previously stored posterior distribution for each examinee is
          being used in this function).
     */
    template <class EI> Real DoEStep(EI examinees_begin, EI examinees_end,
        bool computeExamineePosterior, bool storeExamineePosterior)
    {
      return DoEStep(examinees_begin, examinees_end, items_begin, items_end,
          computeExamineePosterior, storeExamineePosterior);
    }
#endif // BOOST_MSVC6_MEMBER_TEMPLATES
  };

  /*!
    \brief 
    Class template to perform E-step calculation and store results
    for IRT model with a discrete latent variable distribution.
 
    \section template_args Template Parameters
    
    \param E Examinee type.
    \param I Item type. 
    \param II  Iterator over item objects.
    \param D  Class for discrete latent variable distribution.
    
    \section function_args Function Parameters
    
    \param[in] bitem Iterator pointing to first item.
    \param[in] eitem Iterator pointing to one past last item.
    \param[in]  &dist Address of latent variable distribution object.
   */
  template <class E, class I, class II, class D>
  EStepDiscrete<E, I, II, D>::EStepDiscrete( II bitem, II eitem, D &dist) :
  items_begin(bitem), items_end(eitem), latentvar_dist(dist), nGroups(0),
  logLatentProb( dist.NumGroups(), dist.size()),
  mRespProb(dist.NumGroupsUnique())
  {
    int i;
    numItems = items_end - items_begin;
    numLatentVarCat = dist.size();
    numGroupUnique = dist.NumGroupsUnique();

    notPresentedResponse = I::NotPresentedResponse();

    nGroups = new RealMatrix(dist.NumGroups(), numLatentVarCat);

    for (i=0; i<numGroupUnique; ++i)
    {
      mRespProb[i] = new NMatrixVec<II>(bitem, eitem, dist.size());
    }

    /* Store indices of items in examinee response vector */
    II iitem = items_begin;
    itemIndices = new std::vector<int>(numItems);
    std::vector<int>::iterator ii = itemIndices->begin();
    for (i=numItems; i--; ++ii, ++iitem)
    {
      *ii = (*iitem)->Index();
    }

    /* Initialize logLatentProb to log of probabilities in latentvar_dist */
    for (i = 1; i <= latentvar_dist.NumGroups(); ++i)
    {
      RealMatrix::row_iterator ip = logLatentProb.begin_row(i);
      typename D::weight_iterator iwt = latentvar_dist.begin_weights(i);
      for (int j = numLatentVarCat; j--; ++ip, ++iwt)
      {
        *ip = (*iwt != 0.0) ? std::log(*iwt) : logZero;
      }
    }
  }

  /*! Destructor */
  template <class E, class I, class II, class D> EStepDiscrete<E, I, II, D>::~EStepDiscrete()
  {
    delete nGroups;

    delete itemIndices;

    for (int i = 0; i<numGroupUnique; ++i)
    delete mRespProb[i];
  }

  /*!
    \brief
    Fills NMatrixVec objects with log of probabilities of each response to each item
    for each level of the discrete latent variable.
 
    \section template_args Template Parameters
    
    \param E Examinee type.
    \param I Item type. 
    \param II  Iterator over item objects.
    \param D  Class for discrete latent variable distribution.
   */
  template <class E, class I, class II, class D> void EStepDiscrete<E, I, II, D>::CalcResponseProb()
  {
    for (int g=0; g<numGroupUnique; ++g)
    {
      typename NMatrixVec<II>::iterator vi = mRespProb[g]->begin(); // "typename" keyword added. ww, 1/10/2008
      for (II iitem = items_begin; iitem != items_end; ++iitem, ++vi) // loop over items

      {
        Response response;
        int n = (*iitem)->NumRespCat();
        for (int j = 1; j <= n; ++j) // loop over responses to items

        {
          typename D::point_iterator ipoint = latentvar_dist.begin_points(g+1);
          response = (*iitem)->IndexResponse(j-1);
          RealMatrix::row_iterator il = (**vi).begin_row(j);
          for (int k = numLatentVarCat; k--; ++il, ++ipoint) // loop over latent variable categories

          {
            *il = std::log((*iitem)->ProbResp(response, *ipoint));
          }
        }
      }
    }
  }

  /* Definitions of member templates for compilers which can handle member template definitions
     outside the class declaration.
   */
#ifndef BOOST_MSVC6_MEMBER_TEMPLATES
  /*!
    \brief
    Computes posterior distribution of discrete latent variable for an examinee.
    Returns marginal likelihood of the examinee's responses.
    
    A duplicate of this definition is given inside the class definition
    for use when BOOST_MSVC6_MEMBER_TEMPLATES is defined.
 
    \section template_args Template Parameters
    
    \param PI Iterator over posterior probabilities.
 
    \section function_args Function Parameters
    
    \param[in]	examinee Object holding information about an examinee's item responses
    \param[in]	begin_posterior Iterator pointing to first element of container holding posterior
		    probabilities.
    \param[in]	end_posterior Iterator pointing to one past last element of container holding
		    posterior probabilities. This argument is only used to make sure there
		    is enough space in the container which will hold the posterior probabilities.
   */
  template <class E, class I, class II, class D> template <class PI> Real
  EStepDiscrete<E, I, II, D>::ExamineePosterior( E &examinee,
      PI begin_posterior, PI end_posterior)
  {
    int i, il;

    int group = examinee.Group();

    if ((end_posterior - begin_posterior) != numLatentVarCat)
    {
      throw InvalidArgument("Incorrect size of vector to hold posterior probabilities",
          "EStepDiscrete::ExamineePosterior");
    }

    /* Constants used for loop unrolling */
    int Ndiv4 = numLatentVarCat / 4;
    int Nmod4 = numLatentVarCat - Ndiv4*4;

    /* initialize posterior probabilities */
    RealVector::iterator ipost = begin_posterior;
    RealMatrix::row_iterator iwt = logLatentProb.begin_row(group);
    for (i = Ndiv4; i--; ipost+=4, iwt+=4)
    {
      *ipost = *iwt;
      ipost[1] = iwt[1];
      ipost[2] = iwt[2];
      ipost[3] = iwt[3];
    }
    for (i = Nmod4; i--; ++ipost, ++iwt)
    {
      *ipost = *iwt;
    }

    typename NMatrixVec<II>::iterator item = // "typename" keyword added. ww, 1/10/2008.
    (numGroupUnique == 1) ? mRespProb[0]->begin() : mRespProb[group-1]->begin();
    II iitem = items_begin;
    typename E::response_iterator presp = examinee.responses_begin();
    std::vector<int>::iterator ii = itemIndices->begin();
    for (i = numItems; i--; ++item, ++iitem, ++ii)
    {
      Response resp = presp[*ii];
      if (resp != notPresentedResponse)
      {
        ipost = begin_posterior;
        int index = (*iitem)->ResponseIndex(resp);
        RealMatrix::row_iterator ir = (*item)->begin_row(index+1);
        for (il=Ndiv4; il--; ipost+=4, ir+=4)
        {
          *ipost += *ir;
          ipost[1] += ir[1];
          ipost[2] += ir[2];
          ipost[3] += ir[3];
        }
        for (il=Nmod4; il--; ++ipost, ++ir)
        {
          *ipost += *ir;
        }
      }
    }

    /* find sum in order to standardize posterior */
    ipost = begin_posterior;
    iwt = logLatentProb.begin_row(group);
    Real sum = 0.0;
    for (i = numLatentVarCat; i--; ++ipost, ++iwt)
    {
      if (*iwt != logZero)
      sum += std::exp(*ipost);
      else
      *ipost = logZero;
    }

    /* standardize */
    ipost = begin_posterior;
    Real logsum = std::log(sum);
    for (i = numLatentVarCat; i--; ++ipost)
    {
      if (*ipost != logZero)
      {
        *ipost -= logsum;
        *ipost = std::exp(*ipost);
      }
      else
      *ipost = 0.0;
    }

    return sum;
  }

  /*!
    \brief
    Computes the E-step of the EM algorithm for IRT models with a
    discrete latent variable distribution.
    
    Returns marginal loglikelihood of examinees' responses
    (sum over examinees of the marginal loglikelihood of an examinee's responses)
    plus sum of prior likelihoods over all item parameters.
    This is the value of the marginal posterior density that the
    EM algorithm is maximizing at the values of the item parameters
    computed in the last M-step. The log of the priors for the item
    parameters are added to this value for the items for which
    n and r are calculated.
    
    Results of the E-step are stored in data member nGroups. The posterior
    distributions for examinees are stored in the examinee objects
    if storeExamineePosterior is true, and updated n's and r's are computed
    the items given by itemsNR_begin and itemsNR_end.
    
    A duplicate of this definition is given inside the class definition
    for use when BOOST_MSVC6_MEMBER_TEMPLATES is defined.
 
    \section template_args Template Parameters
    
    \param D  Class for discrete latent variable distribution.
    \param E  Examinee type.
    \param EI Iterator over pointers to examinee objects.
    \param I  Item type. 
    \param II Iterator over item objects.
 
    \section function_args Function Parameters
    
    \param[in]	examinees_begin Iterator to pointer to first examinee
    \param[in]	examinees_end Iterator to pointer to one past last examinee
    \param[in]	itemsNR_begin Iterator to first item pointer for which n and r will be updated.
        The items for which n and r are updated can be different from the items
        used to compute the posterior distribution of the latent variable for
        each examinee. If itemsNR_end - itemsNR_begin == 0 then n and r are not
        updated for any items.
    \param[in]	itemsNR_end Iterator to one past last item pointer for which n and r will be 
        updated.
    \param[in]  computeExamineePosterior If TRUE posterior latent variable distribution
        for each examinee is computed. If FALSE previously stored posterior
        latent variable distribution for each examinee is used.
    \param[in]  storeExamineePosterior If TRUE posterior latent variable distribution
        is stored for each examinee. If 'computeExamineePosterior' is
        FALSE the value of the argument is not used (in that case a
        previously stored posterior distribution for each examinee is
        being used in this function).
   */
  template <class E, class I, class II, class D> template <class EI> Real EStepDiscrete<E, I, II, D>::DoEStep(
      EI examinees_begin, EI examinees_end, II itemsNR_begin, II itemsNR_end,
      bool computeExamineePosterior, bool storeExamineePosterior)
  {
    int i, j;
    Real loglikelihood = 0.0;

    RealVector posterior(numLatentVarCat);

    *nGroups = 0.0;
    int numItemsNR = itemsNR_end - itemsNR_begin;
    II iitem = itemsNR_begin;
    for (i=numItemsNR; i--; ++iitem)
    {
      /* Check that n and r for each item use the same number of latent variable 
         categories as in latent variable distribution used to compute
         examinee posterior distributions (does not check that points match,
         just that number of points match, it is assumed 
         that points also match).
       */
      j = (*iitem)->NumLatentVarCat();
      if (j != numLatentVarCat)
      {
        throw RuntimeError("Mismatch in number of latent variable categories",
            "EStepDiscrete::DoEStep");
      }

      /* Initialize n and r for item to zero */
      (*iitem)->InitializeNR();
    }

    if (computeExamineePosterior)
    {
      // Check that number of latent variable categories for latentvar_dist
      // has not changed.
      if (numLatentVarCat != latentvar_dist.size())
      {
        throw RuntimeError("Number of latent variable categories has changed",
            "EStepDiscrete::DoEStep");
      }

      // Compute log probabilities of each response to each item using
      // current item parameter estimates
      CalcResponseProb();

      /* Set logLatentProb to log of current probabilities of latentvar_dist */
      for (i = 1; i <= latentvar_dist.NumGroups(); ++i)
      {
        RealMatrix::row_iterator ip = logLatentProb.begin_row(i);
        typename D::weight_iterator iwt = latentvar_dist.begin_weights(i);
        for (j = numLatentVarCat; j--; ++ip, ++iwt)
        {
          *ip = (*iwt != 0.0) ? std::log(*iwt) : logZero;
        }
      }
    }

    /* For each examinee compute posterior distribution and
       update n and r for items the examinee responded to.
     */
    for (EI examinee_i = examinees_begin; examinee_i != examinees_end; ++examinee_i)
    {

      Real marginalLikelihood;
      if (computeExamineePosterior) // Compute posterior distribution for examinee

      {
        marginalLikelihood = ExamineePosterior(**examinee_i, posterior.begin(), posterior.end());

        if (storeExamineePosterior)
        {
          typename E::posterior_vector epost(numLatentVarCat);
          typename E::posterior_vector::iterator iep = epost.begin();
          RealVector::iterator ip = posterior.begin();
          for (i = numLatentVarCat; i--; ++iep, ++ip)
          *iep = *ip;
          (*examinee_i)->SetPosterior(epost);

          (*examinee_i)->SetMarginalRespLikelihood(marginalLikelihood);
        }
      }
      else // use examinee posterior distribution already computed

      {
        typename E::posterior_vector::iterator iep = (*examinee_i)->posterior_begin();
        RealVector::iterator ip = posterior.begin();
        for (i = numLatentVarCat; i--; ++iep, ++ip)
        *ip = *iep;

        marginalLikelihood = (*examinee_i)->GetMarginalRespLikelihood();
      }

      /* update marginal loglikelihood */
      loglikelihood += std::log(marginalLikelihood);

      typename E::response_iterator iresp = (*examinee_i)->responses_begin();
      Real casewt = (*examinee_i)->Count();
      int group = (*examinee_i)->Group();
      iitem = itemsNR_begin;
      for (i = numItemsNR; i--; ++iitem)
      {
        /* Update n and r for each item */
        Response resp = iresp[(*iitem)->Index()];
        if (resp != notPresentedResponse)
        {
          typename I::r_iterator ir = (*iitem)->RVector(resp, group);
          typename I::n_iterator in = (*iitem)->NVector(group);
          RealVector::iterator ipost = posterior.begin();
          for (j = numLatentVarCat; j--; ++ir, ++in, ++ipost)
          {
            *ir += *ipost * casewt;
            *in += *ipost * casewt;
          }
        }

      }

      /* Update marginal distribution for group examinee belongs to */
      RealVector::iterator ipost = posterior.begin();
      RealMatrix::row_iterator igroup = nGroups->begin_row(group);
      for (j = numLatentVarCat; j--; ++ipost, ++igroup)
      {
        *igroup += *ipost * casewt;
      }

    }

    /* Add log of prior densities of item parameter estimates to loglikelihood */
    for (II ii = itemsNR_begin; ii != itemsNR_end; ++ii)
    {
      PriorVector::iterator iprior = (*ii)->PriorsIterator();
      RealVector::iterator iparam = (*ii)->ParametersIterator();
      for (i = (*ii)->NumParameters(); i--; ++iprior, ++iparam)
      {
        if (*iprior)
        loglikelihood += (*iprior)->LogDensity(*iparam);
      }
    }

    return loglikelihood;

  }

  /*!
    \brief
    Version of DoEStep in which the items that are used to compute the
    posterior distributions for examinees are also the items for which
    n and r are updated, i.e., itemsNR_begin == items_begin and
    itemsNR_end == items_end.
 
    A duplicate of this definition is given inside the class definition
    for use when BOOST_MSVC6_MEMBER_TEMPLATES is defined.
 
    \section template_args Template Parameters
    
    \param D  Class for discrete latent variable distribution.
    \param E  Examinee type.
    \param EI Iterator over pointers to examinee objects.
    \param I  Item type. 
    \param II Iterator over item objects.
 
    \section function_args Function Parameters
    
    \param[in]  examinees_begin Iterator to pointer to first examinee
    \param[in]  examinees_end Iterator to pointer to one past last examinee
    \param[in]  computeExamineePosterior If TRUE posterior latent variable distribution
        for each examinee is computed. If FALSE previously stored posterior
        latent variable distribution for each examinee is used.
    \param[in]  storeExamineePosterior If TRUE posterior latent variable distribution
        is stored for each examinee. If 'computeExamineePosterior' is
        FALSE the value of the argument is not used (in that case a
        previously stored posterior distribution for each examinee is
        being used in this function).
   */
  template <class E, class I, class II, class D> template <class EI> Real EStepDiscrete<E, I, II, D>::DoEStep(
      EI examinees_begin, EI examinees_end, bool computeExamineePosterior,
      bool storeExamineePosterior)
  {
    return DoEStep(examinees_begin, examinees_end, items_begin, items_end,
        computeExamineePosterior, storeExamineePosterior);
  }

#endif // BOOST_MSVC6_MEMBER_TEMPLATES
} // namespace etirm

#endif // ETIRM_ESTEPDISCRETE_H_
