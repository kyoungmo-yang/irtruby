#* irt_algorithm.rb
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
#* The IRTEmAlgorithm class provides some options to estimate item paramters(i.g., difficulty, discrimet and guessing) and examinee's theta
#* Author: KyoungMo Yang(E-mail : mo@embian.com)
#* Modified: 2010-11-05
#* References: Irt_ruby module
class IRTAlgorithm
  EAP = "EAP"
  MLE = "MLE"
  
  attr_accessor :max_iter, :crit, :estimates_dist, :scale_points, :no_initial_estep
  attr_accessor :no_mstep_iter_error, :estimates_mean_sd_of_dist, :estimates_mean_of_dist
  attr_accessor :max_optimize_iter, :algorithm_to_estimates_theta, :em_iter_logs
  attr_accessor :use_all, :use_all_items, :ignore_error
  attr_accessor :error_items
  
  #* Functionailities
  # * Constructs an instance of IRTAlgorithm class and initialzes options.
  #* Parameters 
  # * [irt_ruby_module] : an instance of Irt_ruby::IRTModule
  #* Return
  # * An instance of IRTAlgorithm class
  #* Note
  # * N/A
  def initialize(irt_ruby_module)
    #init opitons for algorithm
    _init(irt_ruby_module)
    #init options for dichotomous model
    _init_dichotomous_options
  end
  
  
  #* Functionailities
  # * Returns Irt_ruby::IRTModule instance
  #* Parameters 
  # * N/A
  #* Return
  # * [irt_ruby_module] : an instance of Irt_ruby::IRTModel
  #* Note
  # * N/A
  def irt_ruby_module
    return @irt_ruby_module
  end
  
  #* Functionailities
  # * Estimates items' params(i.e., difficulty, discriment, guessing)
  #* Parameters 
  # * [items] : array of IRTItem's instances
  # * [use_dichotomous_model] : if true, estimates starting values. Default value is true.
  #* Return
  # * N/A
  #* Note
  # * N/A
  def estimates_item_params(items, use_dichotomous_model=true)
    @em_iter_logs.clear
    
    number_of_items = items.size
    #if > 1, the latent distributions for each group can be estimated with estimates_dist option
    number_of_groups = @irt_ruby_module.num_groups
    
    #Ruby wrapper for item_3PL_starting_values that handles errors
    if(use_dichotomous_model)
      _estimates_initial_values(number_of_items)
    end
    
    # create estep object using all items
    @estep_obj = Irt_ruby::Estep.new(@irt_ruby_module.GetEtirmObj())
    # first E-step
    if(@no_initial_estep == false)
      @estep_obj.compute
    end
    
    # EM iterations
    rel_diff = 0
    1.upto(@max_iter){|iter|
      # M-step for item parameters
      rel_diff = _mstep_item_param([], @no_mstep_iter_error)
      # M-step for discrete latent variable distribution
      if(@estimates_mean_of_dist || @estimates_mean_sd_of_dist)
        # Estimate only mean and s.d. of distribution
        reld_diff = _mstep_latent_dist_moments(@estep_obj, @estimates_mean_of_dist, @estimates_dist)
      else
        # Estimate latent distribution probabilities
        reld_diff = _mstep_latent_dist(@estep_obj, number_of_items, number_of_groups, @estimates_dist, @scale_points)
      end
      # next E-step
      log_like = @estep_obj.compute
      
      #logging results of EM iterations
      @em_iter_logs[iter] = [rel_diff, log_like]
      
      if(rel_diff < @crit)
        break
      end
    }
    
    #Report warning if convergence criterion not met
    if(rel_diff >= @crit)
      #handling error or warning
      if(! @no_mstep_iter_error)
        raise("\nConvergence criterion not met after #{@max_iter} EM iterations\n")
      end
    end
    
    #Setting item paramters
    _setting_item_params(items)
    @estep_obj = nil
  end
  
  #* Functionailities
  # * Estimates examinees' theta. This methods can used after estimates item parameters.
  #* Parameters 
  # * [examinees] : array of IRTExaminees's instances
  # * [range_of_ability] : range of ability such as -4..4
  #* Return
  # * N/A
  #* Note
  # * N/A
  def estimates_examinee_theta(examinees, range_of_ability)
    @estep_obj = nil
    @estep_obj = Irt_ruby::Estep.new(@irt_ruby_module.GetEtirmObj())
    @estep_obj.compute(true, true)
    @estep_obj = nil
    if(@algorithm_to_estimates_theta == EAP)
      _estimates_EAP_theta(examinees)
    else
      _estimates_MLE_theta(examinees, range_of_ability)
    end
  end
  
  
  #############################
  # private methods
  ############################# 
  private
  #* Functionailities
  # * Setting item paramters for each item
  #* Parameters 
  # * [items] : array of IRTItem's instances
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _setting_item_params(items)
    items.each do |item|
      params = @irt_ruby_module.item_get_all_params(item.no)
      #discriment
      item.irt_a = params[0]
      #difficulty
      item.irt_b = params[1]
      #guessing
      item.irt_c = params[2]
    end
  end
  
  #* Functionailities
  # * Initializes default otpions
  #* Parameters 
  # * [irt_ruby_module] : an instance of Irt_ruby::IRTModule
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _init(irt_ruby_module)
    #IRTRuby Module
    @irt_ruby_module = irt_ruby_module
    
    @error_items = Hash.new
    #logging result of em_steps
    @em_iter_logs = Hash.new
    
    #algorithm for estimates examinee theta, possible values = ["MLE", "EAP"]
    #default is MLE
    @algorithm_to_estimates_theta = EAP
    # maximum number of EM iterations, default is 20
    @max_iter = 20
    
    #criterion, default = 0.001
    @crit = 0.001
    
    #if this value is true and number of groups are more that one, 
    #estimates latent variable distribution(of ability) in base group, default is false
    # @estimates_dist = false
    @estimates_dist = true
    
    #if true, points of ability(distribution) for the base group are scaled,
    #so the mean is 0 and s.d. is 1 after the weights(prob) are estimatesd in multiple group estimation.
    #The item parameter estimatess are also rescaled using the same scale transformation.
    #The option only has an effect if the -estimates_dist option is used.
    # @scale_points = false
    @scale_points = true
    
    #if true, do not compute one E-step before loop over EM interations
    @no_initial_estep = false
    
    #if true, an error of the maximum of iterations being exceeded in the optimization performed, 
    #so the M-step is not considered an error and the calculation will continue.
    @no_mstep_iter_error = true
    
    #if true, the mean and s.d. of the latent distributions(of abilities) for all groups(except the base group)
    #are estimatesd rather than the discrete probabilities.
    #To use the option, the unique_points option must be used(in IRTModel class)
    @estimates_mean_sd_of_dist = false
    
    #if ture, the mean of the latent distributions(of abilities) for all groups(except base group)
    #are estimatesd rather than the discrete probabilities.
    #The s.d. is not estimatesd.
    #To use the option, the unique_points option must be used(in IRTModel class)
    @estimates_mean_of_dist = false
    
    # Set maximum number of iterations used by optimization procedure
    # in mstep for each item. default is 150 (SEE : Uncmin.h)
    @max_optimize_iter = 150
  end
  
  #* Functionailities
  # * Initializes default options for dichotomous model
  #* Parameters 
  # * N/A
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _init_dichotomous_options
    #if ture, all examinees and items are used to compute initial thetas used for starting values,
    #even examinees who get all items correct or all items incorrect. default is false
    @use_all = false
    
    #if ture, using all items
    @use_all_items = false
    
    #if false, program will stop if an error occurred in computing start values,
    #otherwise the error is just reported and the program continues. default is false
    @ignore_error = false
  end
  
  def _add_error_item_no(no)
    @error_items[no] = no
  end
  
  #* Functionailities
  # * Estimates initial values for dichotnomous model
  #* Parameters 
  # * [number_of_items] : the number of items
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _estimates_initial_values(number_of_items)
    item_nos = (1..number_of_items).to_a
    begin
      @irt_ruby_module.item_3PL_starting_values(@use_all, @use_all_items, item_nos)
    rescue => err
      puts "#{err} In item_3PL_starting_values"
    end
    err_maxiter = Array.new
    err_other = Array.new
    item_nos.each{ |item_no|
      err = @irt_ruby_module.mstep_message(item_no)
      #what means "err < 0"? and what means "err > 4"?
      if(err < 0 || err >= 4)
        _add_error_item_no(item_no)
        # assumes err == 4 means maximum number of iterations exceeded
        # SEE : MStemIRT.h
        if(err == 4)
          err_maxiter << item_no
        elsif
          err_other << item_no
        end
      end
    }
    
    if(! err_maxiter.empty?)
      #error handling
      puts "Maximum number of iterations exceeded in which compute starting values for items : #{err_maxiter.join(", ")}"
    end
    
    if(! err_other.empty?)
      #error handling
      puts "Minimization defedure used to compute starting values failed for items : #{err_other.join(", ")}"
      if(@ignore_error == false)
        raise("Starting values not computed for #{e} items in command starting_values")
      else
        puts "Starting values not computed for #{e} items in command starting_values"
      end
    end
  end
  
  #* Functionailities
  # * Computes M-step for items.
  #* Parameters 
  # * [item_nos] : array of item numbers(question_no or index) for items to perform M-steps. Default is empty that indicates using all items.
  # * [no_max_iter_error] : if true, an error in the optimization performed in the M-step is not considered an error and the calculation will continue.
  #* Return
  # * [mstep_max_diff] : maximum relative difference in paramter from lat iteration to current iteration.
  #* Note
  # * N/A
  def _mstep_item_param(item_nos=[], no_max_iter_error=false)
    if(item_nos.empty?)
      error_item = @irt_ruby_module.mstep_items(no_max_iter_error)
    else
      error_item = @irt_ruby_module.mstep_items(no_max_iter_error, item_nos)
    end
    if(error_item > 0)
      message = @irt_ruby_module.mstep_message(error_item)
      raise("M-step failed for item #{error_item} with error number #{message}")
    end
    
    #Display item number for which maximum number of iterations was exceeded
    if(error_item < 0)
      if(item_nos.empty?)
        number_of_items = @irt_ruby_module.num_items
        1.upto(number_of_items){|item_no| item_nos << item_no}
      end
      err_maxiter = Array.new
      item_nos.each do |item_no|
        err = @irt_ruby_module.mstep_message(item_no)
        if(err == 4)
          _add_error_item_no(item_no)
          err_maxiter << item_no
        end
      end
      
      if(err_maxiter.size > 0)
        puts "Maximum nuber of iterations exceeded in M-step calculation for items: #{err_maxiter.join(", ")}"
      end      
    end
    return @irt_ruby_module.mstep_max_diff()
  end
  
  #* Functionailities
  # * Computes M-step to estimates mean and standard deviation of latent distributions.
  #* Parameters 
  # * [estep_obj] : an instance of Irt_ruby::Estep
  # * [mean_only] : if true, only estimates the means of latent distribution and do not estimates the standard deviations.
  # * [estim_base_group] :  - if ture, estimates latent distribution for the base group.
  #* Return
  # * [[mean_diff, sd_diff]] : maximum relative difference in the latent distribution means from last iteration to current iteration.
  #* Note
  # * N/A
  def _mstep_latent_dist_moments(estep_obj, mean_only, estim_base_group)
    # If only one group then there is nothing to estimates
    num_groups = @number_of_groups
    if(num_groups == 1)
      return []
    end
    
    #initialize differences
    mean_diff = 0
    sd_diff = 0
    
    #list of indices for latent distributions points
    num_latent_dist_points = @points_of_ability
    point_list = (1..num_latent_dist_points).to_a
    
    #Compute M-step for discrete latent variable distribution in each group
    base_group = @irt_ruby_module.get_base_group()
    1.upto(num_groups){|group|
      #skip base group when estim_base_group is false
      if(estim_base_group || group != base_group)
        #Compute moments based on original points
        moments = @irt_ruby_module.dist_mean_sd(group)
        old_mean = moments[0]
        old_sd = moments[1]
        
        #store original probabilities for group
        old_prob = @irt_ruby_module.dist_get_probs(group)
        
        #M-step for probabilities
        @irt_ruby_module.mstep_dist(estep_obj, group)
        
        #Compute moments based on new probabilities
        moments = @irt_ruby_module.dist_mean_sd(group)
        mean = moments[0]
        sd = moments[1]
        
        #Compute difference in old and new mean
        diff = (mean-old_mean).abs
        if(diff > mean_diff)
          mean_diff = diff
        end
        
        #Compute difference in old and new s.d.
        if(!mean_only)
          diff = (sd-old_sd).abs
          if(diff > sd_diff)
            sd_diff = diff
          end
        end
        
        #Compute slope and intercept to transform points
        #so distribution has new mean and s.d.
        if(mean_only)
          slop = 1.0
        else
          slop = sd / old_sd
        end
        intercept = mean - slope * old_mean
        
        #Transform points for group distribution
        point_list.each do |point|
          old_point = @irt_ruby_module.dist_get_point(point, group)
          @irt_ruby_module.dist_set_point(point, old_point*slope+intercept, group)
        end
        
        #Restore original probabilities for group
        @irt_ruby_module.dist_set_probs(old_prob, group)
      end
    }
    return [mean_diff, sd_diff]
  end
  
  #* Functionailities
  # * Computes M-step for latent distributions from last iteration to current iteration.
  #* Parameters 
  # * [estep_obj] : an instance of Irt_ruby::Estep
  # * [number_of_items] : the number of items
  # * [number_of_groups] : the number of groups. default value is one.
  # * [estim_base_group] :  If true, estimates latent distribution for the base group.
  # * [scale_points] :  If true then scale points of latent distribution so mean and s.d. in base group are 0 and 1, and correspondingly scale item parameter estimatess in multiple group estimation.  This option only has an effect when the "estim_base_group" option is used.
  #* Return
  # * [reld_diff] : maximum relative difference in latent distribution probabilities.
  #* Note
  # * N/A
  def _mstep_latent_dist(estep_obj, number_of_items, number_of_groups=1, estim_base_group=false, scale_points=false)
    reld_diff = Array.new
    base_group = @irt_ruby_module.get_base_group()
    
    #Compute M-step for discrete latent variable distribution in base group
    if(estim_base_group)
      reld_diff << @irt_ruby_module.mstep_dist(estep_obj, base_group)
      
      #Transform latent variable scale so that mean and s.d. of latent variable for base group are 0 and 1.
      #The points of the latent variable distribution are changed to reflect the new scale as are the item
      #paramter estimatess.
      if(scale_points)
        _standardize_scale(0.0, 1.0, number_of_items, base_group)
      end
    else
      reld_diff << -1.0
    end
    
    #Estimate distributions for group 2 through last group
    if(number_of_groups > 1)
      1.upto(number_of_groups){|group|
        if(base_group != group)
          t_diff = @irt_ruby_module.mstep_dist(estep_obj, group)
          if(t_diff > reld_diff[0])
            reld_diff[0] = t_diff
          end
        end
      }
    end    
    return reld_diff        
  end
  
  #* Functionailities
  # * Standardize scale of latent variable so that the mean and standard deviation are equal to specific values in one examinee group. Item parameters for all items are correspondingly transformed. If different latent variable points are used for different groups,  the points for all groups are transformed to be on the new scale.
  #* Parameters 
  # * [mean] : Target mean of latent variable distribution.
  # * [sd] : 
  # * [number_of_items] : the number of items
  # * [group] :  Number of group in which mean and s.d. of latent variable distribution should be equal to target values. (1, 2, ...)
  # * [ignore_prior] :  If true then do not report an error if a transformed parameter has zero density in the prior used for that parameter.
  #* Return
  # * [[slope, intercept]] : slope and intercept of scale transformation.
  #* Note
  # * N/A
  def _standardize_scale(mean, sd, number_of_items, group=1, ignore_prior=false)
    #check for valid s.d.
    if(sd <= 0.0)
      raise("Invalid s.d. : #{sd}")
    end
    
    #Find transformation of points of latent distribution that will give the specified mean and s.d.
    trans = @irt_ruby_module.dist_scale(mean, sd, group)
    slope = trans[0]
    intercept = trans[1]
    
    #Transform item parameters using scale transformation computed for points of latent variable distribution
    1.upto(number_of_items){|item_no|
      if(@irt_ruby_module.item_scale_params(item_no, slope, intercept, ignore_prior) != 0)
#        raise("Transformation would result in an invalid paramter for item #{item_no}")
        puts ("Transformation would result in an invalid paramter for item #{item_no}")
      end
    }
    
    #Return slope and intercept of scale transformation
    return [slope, intercept]
  end
  
  #* Functionailities
  # * Estimates the MLE theta for each examinee.
  #* Parameters 
  # * [examinees] : array of IRTExaminee
  # * [range_of_ability] : range of ability such as -4..4
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _estimates_MLE_theta(examinees, range_of_ability)
    #compute MLE estimates and number correct for each examinee
    examinees.each do |examinee|
      #Get examinee MLE estimate
      min_theta = range_of_ability[0]
      max_theta = range_of_ability[1]
      mle = @irt_ruby_module.examinee_theta_MLE(examinee.no, min_theta, max_theta)
      examinee.theta = mle
    end
  end
  
  #* Functionailities
  # * Estimates the EAP theta for each examinee.
  #* Parameters 
  # * [examinees] : array of IRTExaminee
  # * [estep_items] : array of items' numbers for which n and r are updated. If nil is passed then n arnd r are updated for all items used in the E-step to compute examinee posterior distributions.
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _estimates_EAP_theta(examinees, estep_items=[])
    examinees.each do |examinee|
      #Get examinee posterior mean (EAP estimate)
      eap = @irt_ruby_module.examinee_posterior_mean(examinee.no)
      examinee.theta = eap
    end
  end
end