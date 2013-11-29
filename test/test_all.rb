require 'test/unit'
require './irt_item_test'
require './irt_examinee_test'
require './irt_exam_test'
require './irt_algorithm_test'
require './irt_model_test'
require './irt_test'

require '../irt'
require '../irt_item'
require '../irt_examinee'
require '../irt_exam'
require '../irt_model'
require '../irt_algorithm'

#* test_all.rb
#* Copyright (C) 2010  Embian Inc.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or 
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with This program.  If not, see <http://www.gnu.org/licenses/>.

class TestAll < Test::Unit::TestCase
  def setup
    #Do nothing
  end
  
  def teardown
    #Do nothing
  end
  
  def test_all
    #Call Other Unit Tests
  end
  
  def test_with_realistic_data
=begin
    Example data from Chapters 4 and 6 Kolen and Brennan (1995).
    
    Estimate item parameters for Form Y and Form X items while at the same time estimating 
    the latent variable distribution of the groups that took Form X and Form Y. 
    
    After each M-step the points of the discrete latent variable distribution are linearly 
    transformed so the mean and standard deviation of the Form X distribution are 0 and 1.
    
    The item parameters are also tranformed using the same scale transformation.
=end

    # Read examinee item responses from file mondat.dat.
    # Each record contains the responses to 
    # 60 items for an examinee in columns 2-61.
    items = []
    60.times.each do |no|
      items << IRTItem.new("id#{no+1}")
    end
    
    # The first 24 items are the unique items on
    # Form Y, the second 12 items are common items,
    # and the last 24 items are unique items on
    # Form X. An integer in column 1 gives
    # the examinee group: 1 for examinees
    # who took Form Y, and 2 for examinees
    # who took Form X
    
    examinees = []
    no = 1
    File.open('./mondat.dat', 'r') do |file|
      while line = file.gets
        group = line[0].chr.to_i
        responses = _get_responses(line[1..(line.size-1)])
        examinees << IRTExaminee.new("id##{no}", responses,"00000-0000#{no}", group)
        no = no.next
      end
    end

    # 24 unique items on each of two forms and
    # 12 common items for a total of 60
    # items. Two groups specified
    # for multiple group estimation
    exam = IRTExam.new(examinees, items, 2)
    
    # Perform EM iterations for computing item parameter estimates
    # and probabilities of latent variable distributions for
    # groups 1 and 2. Scale points of latent variable distribution
    # after each M-step so the mean and s.d. in group 1 are
    # 0 and 1. Allow a maximum of 200 EM iterations.
    
    model = IRTModel.new(IRTModel::DICHOTOMOUS_2PL)
    algorithm = IRTAlgorithm.new(Irt_ruby::IRTModule.new)
    algorithm.estimates_dist = true
    algorithm.scale_points = true
    algorithm.max_iter = 200
    
    irt = IRT.new(exam, model, algorithm)
    
    #estimates item paramters
    irt.estimates_item_params
    irt.estimates_examinee_theta
    # @irt.estimates_examinee_score

    # Write item parameter estimates, discrete latent
    # variable distributions, and mean and s.d. of
    # latent variable distributions.
    File.open('./mondat.log', 'w') do |file|
      file.puts irt.to_s
    end

    # end of run
    irt.release
  end
  
  private
  def _get_responses(str_responses)
    responses = []
    str_responses.each_char do |res|
      case res
      when '.'
        responses << IRTExam::MISSING_RESPONSE
      when '1'
        responses << IRTExam::CORRECT_RESPONSE
      when '0'
        responses << IRTExam::INCORRECT_RESPONSE
      end
    end
    return responses
  end

end