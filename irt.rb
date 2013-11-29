require File.dirname(__FILE__) + '/lib/irt_ruby'
require File.dirname(__FILE__) + '/irt_model'
require File.dirname(__FILE__) + '/irt_algorithm'

#* irt.rb
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
#  
#* IRT class includes methods to estimates item parameters and thetas of examinees based on Item Response Theory.
#* Author: KyoungMo Yang(E-mail : mo@embian.com)
#* Modified: 2010-11-05
#* References: Irt_ruby module, IRTModel class and IRTAlgorithm class
class IRT
  
  attr_accessor :test
  
  #* Functionailities
  # * Constructs an instance of IRT class.
  #* 
  # * [test] An instance of IRTExam class
  # * [model] An instance of IRTModel class. if do not passing this parameter, use default model.
  # * [algorithm] An instance of IRTAlgorithm class. if do not passing this parameter, use default algorithm.
  # * [inc] increment between adjacent scores
  #* Return
  # * An instance of IRT class
  #* Note
  # * N/A
  def initialize(test, model=nil, algorithm=nil)
    # init test
    @test = test
    
    # init default Expectation-Maximization and MLE, EAP algorithm
    _init_algorithm(algorithm)
    
    # init default model
    _init_model(model)
    
    # init  information of items, examinees and model
    _init_irt
  end
    
  #* Functionailities
  # * Estimates item parameters(i.e., difficulty, discriment, guessing)
  #* Parameters 
  # * N/A
  #* Return
  # * N/A
  #* Note
  # * N/A
  def estimates_item_params
    @algorithm.estimates_item_params(@test.items, @model.use_dichotomous_model?)
  end
  
  def em_iter_logs
    return @algorithm.em_iter_logs
  end
  
  #* Functionailities
  # * Estimates examinees' theta
  #* Parameters 
  # * [examinees] : if do not passing this parameter, estimates thetas of all examinees
  #* Return
  # * N/A
  #* Note
  # * N/A
  def estimates_examinee_theta(examinees = nil)
    if examinees
      @examinee_abilities = @algorithm.estimates_examinee_theta(examinees, @model.range_of_ability)
    else
      @examinee_abilities = @algorithm.estimates_examinee_theta(@test.examinees, @model.range_of_ability)  
    end
  end
  
  #* Functionailities
  # * Load estimated item parameters
  #* Parameters 
  # * N/A
  #* Return
  # * N/A
  #* Note
  # * N/A
  def load_estimated_item_params
    @test.items.each do |item|
      @irt_ruby_module.item_set_all_params(item.no, [item.irt_a, item.irt_b, item.irt_c])
    end
  end
  
  #* Functionailities
  # * Estimates examinees' score
  #* Parameters 
  # * [score_per_items] : if this paramter is true, calculates score of each items for examinees. Otherwise, calculates total score for each examinees.
  # * [item_nos] : if this parameter is specified, calculates score for the items. Otherwise, calculates score for all items
  # * [examinees] : if this parameter is specified, calculates socre for the examinees. Otherwise, calculates score for all examinees.
  #* Return
  # * N/A
  #* Note
  # * N/A
  def estimates_examinee_score(score_per_items = false, item_nos = nil, examinees = nil)
    if  examinees
      s_examinees = examinees
    else
      s_examinees = @test.examinees
    end
    
    s_examinees.each do |examinee|
      theta = examinee.theta
      if(theta == nil)
        raise("Should be estimate score afeter estimates examinees' thetas")
      end
      
      if item_nos
        score = @irt_ruby_module.test_characteristic_curve([theta], item_nos)
        responses = examinee.responses[(item_nos[0]-1)..(item_nos[item_nos.size-1]-1)]
      else
        score = @irt_ruby_module.test_characteristic_curve([theta])
      end
      examinee.score = score[0]
      
      if(score_per_items)
        #scores of each items for an examinee
        unless item_nos
          item_nos = (1..(@test.number_of_items)).to_a
        end
        item_nos.each do |no|
          examinee.scores[@test.get_item_no(no-1)] = @irt_ruby_module.test_characteristic_curve([theta], [no])[0]
        end
      end
    end
  end
  
  #* Functionailities
  # * Returns result of EM iteration
  #* Parameters 
  # * N/A
  #* Return
  # * N/A
  #* Note
  # * N/A
  def em_logs
    return @algorithm.em_iter_logs
  end
  
  #* Functionailities
  # * Release memory allocated in the _init_irt method
  #* Parameters 
  # * N/A
  #* Return
  # * N/A
  #* Note
  # * N/A
  def release
    @irt_ruby_module.delete_items_dist()
  end
  
  #* Functionailities
  # * Returns a string representing instance
  #* Parameters 
  # * N/A
  #* Return
  # * [result] : a string representing instance
  #* Note
  # * N/A
  def to_s
    s = "Number of items: #{@test.number_of_items}\n"
    s += "Number of latent variable points: #{@model.points_of_ability}\n"
    s += "Number of examinee groups: #{@test.number_of_groups}\n"
    s += "Number of examinees: #{@test.number_of_examinees}\n\n"
    
    s += "Default prior for a-parameters:\n"
    params = @irt_ruby_module.get_default_prior_param('a')
    s += "beta a: %.3f b: %.3f lower limite: %.3f upper limit: %.3f\n" % [params[0], params[1], params[2], params[3]]
    s += "Default prior for b-parameters:\n"
    params = @irt_ruby_module.get_default_prior_param('b')
    s += "beta a: %.3f b: %.3f lower limite: %.3f upper limit: %.3f\n" % [params[0], params[1], params[2], params[3]]
    s += "Default prior for c-parameters:\n"
    params = @irt_ruby_module.get_default_prior_param('c')
    s += "beta a: %.3f b: %.3f lower limite: %.3f upper limit: %.3f\n" % [params[0], params[1], params[2], params[3]]
    s += "EM iterations\n"
    logs = em_iter_logs
    len = logs.size
    1.upto(len).each do |i|
      s += "\t#{i}: %.6f %.6f\n" % logs[i]
    end
    
    s += "\nItem Parameter Estimates (a, b, c)\n"
    @test.items.each do |item|
      s += "#{item.no}\t%.6f\t%.6f\t%.6f\n" % [item.irt_a, item.irt_b, item.irt_c]
    end
    
    s += "\nDiscrete Latent Variable Distributions (point, prob)\n"
    1.upto(@test.number_of_groups).each do |i|
      s += "Group #{i}\n"
      points = @irt_ruby_module.dist_get_points(i)
      probs = @irt_ruby_module.dist_get_probs(i)
      idx = 0
      points.each do |point|
        s += "%.6f\t%.6f\n" % [point, probs[idx]]
        idx = idx.next
      end
    end
    
    s += "\nMoments of Latent Variable Distributions\n"
    1.upto(@test.number_of_groups).each do |i|
      s += "Group #{i}\n"
      mean_sd = @irt_ruby_module.dist_mean_sd(i)
      s += "Mean: %.6f\ts.d.: %.6f\n" % [mean_sd[0], mean_sd[1]]      
    end
    return s
  end
  
  ################################
  # private methods
  ################################
  private
  #* Functionailities
  # * Adds examinee and the examinee's responses
  #* Parameters 
  # * [examinees] : an instance of IRTExaminee
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _add_examinee(an_examinee)
    @irt_ruby_module.add_examinee(an_examinee.responses, an_examinee.group)
  end
  
  #* Functionailities
  # * Adds examinees and the responses of each examinee
  #* Parameters 
  # * [examinees] : array of examinees
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _add_examinees(examinees)
    no = (@irt_ruby_module.examinees_count + 1).to_i
    examinees.each do |examinee|
      examinee.no = no
      _add_examinee(examinee)
      no = no.next
    end
  end
  
  #* Functionailities
  # * Initializes model and options.
  #* Parameters 
  # * [examinees] : array of examinees
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _init_model(a_model = nil)
    if(a_model)
      @model = a_model
    else
      @model = IRTModel.new
    end
    if (@model.use_dichotomous_model?)
      @irt_ruby_module.set_default_model_dichotomous(@model.model)
    else
      @irt_ruby_module.set_default_model_polytomous(@model.model)
    end
    # _init_prior_params(@model)
    @irt_ruby_module.set_default_D(@model.D)
  end
  
  #* Functionailities
  # * Initializes prior item parameters.
  #* Parameters 
  # * [model] : instance of IRTModel
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _init_prior_params(model)
   param = model.prior_a
   @irt_ruby_module.set_default_prior("a", param[:type], param[:value])
   
   param = model.prior_b
   @irt_ruby_module.set_default_prior("b", param[:type], param[:value])
   
   param = model.prior_c
   @irt_ruby_module.set_default_prior("c", param[:type], param[:value])
  end
  
  #* Functionailities
  # * Initialize estimation algorithm options.
  #* Parameters 
  # * [algorithm] : instance of IRTAlgorithm
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _init_algorithm(algorithm = nil)
    if(algorithm)
      @algorithm = algorithm
    else
      @algorithm = IRTAlgorithm.new(Irt_ruby::IRTModule.new)
    end
    @irt_ruby_module = @algorithm.irt_ruby_module
  end
  
  #* Functionailities
  # * Initializes information of items, examinees and model
  #* Parameters 
  # * N/A
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _init_irt
    number_of_items = @test.number_of_items
    number_of_groups = @test.number_of_groups
    points_of_ability = @model.points_of_ability
    range_of_ability = @model.range_of_ability
    use_unique_points = @model.use_unique_points
    
    #array of integers that indicate models for each item,
    #where 1 is dichotomous model, >1 polytomous model. default is dichotomous (SEE : irt_ruby.h)
    if (@model.use_dichotomous_model?)
      #setting defaul dichotomous model
      #@models = [1, 1, 1, 1, 1, ... , number_of_items]
      models = Array.new(number_of_items, 1)
    else
      #if you handles polytomous moudle, you should modify below code
      # the second parameter '2' means category of item 
      models = Array.new(number_of_items, 2)
    end
    #init information for items and latent variable(ability) distribution.
    @irt_ruby_module.new_items_dist(number_of_items, points_of_ability, number_of_groups, models, range_of_ability, use_unique_points)
    
    #init information for examinees
    _add_examinees(@test.examinees)
  end
  
  #* Functionailities
  # * number formating
  #* Parameters 
  # * [value] : number
  #* Return
  # * [value] : formated float number
  #* Note
  # * N/A
  def _to_format(value)
    if value
      return ("%.5f" % [value.to_f]).to_f      
    end
    return value
  end
end