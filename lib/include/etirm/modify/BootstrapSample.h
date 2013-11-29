/*! \file BootstrapSample.h

  \brief Defines function for creating a bootstrap sample.
 
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

#ifndef ETIRM_BOOTSTRAPSAMPLE_H_
#define ETIRM_BOOTSTRAPSAMPLE_H_

#ifdef ETIRM_NO_DIR_PREFIX
#include "etirmtypes.h"
#else
#include "etirm/etirmtypes.h"
#endif

#include <iterator>

namespace etirm
{

  /*!

   \brief Generate a bootstrap sample for a group of examinees.
   
   Template parameters
   
   \param URNG	Type of function object returning uniform random numbers over the integers 1 to 
   the number of examines. An example of a type that meets the requirements for this function
   is the uniform_int class from the boost library (http://www.boost.org).
   
   \param EI	Type of iterator over pointers to examinee objects.

   Arguments
   
   \param  examinees_begin Iterator to pointer to first examinee.
   \param	examinees_end - Iterator to one past pointer to last examinee.
   */
  template <class URNG, class EI> void BootstrapSample(EI examinees_begin, EI examinees_end,
      URNG &rand)
  {
    // Number of examinees
	// Embian Inc. added "typename" - 2010
    typename std::iterator_traits<EI>::difference_type n = std::distance(examinees_begin, examinees_end);

    if (rand.min() != 1 || rand.max() != n)
    {
      throw InvalidArgument("Invalid range for random number generator",
          "SCPPNT::BootstrapSample");
    }

    // Create vector to hold counts for each examinee from bootstrap sample - initially zero
    std::vector<int> counts(n, 0);

    // generate bootstrap counts
    for (int i=0; i<n; ++i)
    {
      ++counts[rand()-1];
    }

    // Assign counts to examinees
    for (std::vector<int>::iterator ic = counts.begin(); examinees_begin != examinees_end; ++ic,
        ++examinees_begin)
    {
      (*examinees_begin)->SetCount(*ic);
    }
  }

} // namespace etirm

#endif // ETIRM_BOOTSTRAPSAMPLE_H_
