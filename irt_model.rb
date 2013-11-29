#* irt_model.rb
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
#* IRTModel is a class to handle options of irt model 
#* Author: KyoungMo Yang(E-mail : mo@embian.com)
#* Modified: 2010-11-05
#* References: N/A
class IRTModel
  DICHOTOMOUS_1PL = "1PL"
  DICHOTOMOUS_2PL = "2PL"
  DICHOTOMOUS_3PL = "3PL"
  
  POLYTOMOUS_PCM  = "PCM"
  POLYTOMOUS_GPCM = "GPCM"
  
  attr_accessor :model, :range_of_ability, :points_of_ability, :D
  attr_accessor :prior_a, :prior_b, :prior_c, :use_unique_points
  
  #* Functionailities
  # * Constructs an instance of IRTModel class and initializes options.
  #* Parameters 
  # * [new_model] : model name. default value is IRTModel::DICHOTOMOUS_2PL
  #* Return
  # * An instance of IRTModel class
  #* Note
  # * N/A
  def initialize(new_model=IRTModel::DICHOTOMOUS_2PL)
    _init_model(new_model)
    _init_prior_for_parameters
  end
  
  #* Functionailities
  # * return value that indicates whether using dichotomous model.
  #* Parameters 
  # * N/A
  #* Return
  # * if true, using dichotmomous model otherwise, return false
  #* Note
  # * N/A
  def use_dichotomous_model?
    if (@model == IRTModel::DICHOTOMOUS_1PL ||
      @model == IRTModel::DICHOTOMOUS_2PL ||
      @model == IRTModel::DICHOTOMOUS_3PL)
      return true
    end
    return false
  end
  
  ################################
  # private methods
  ################################  
  private
  #* Functionailities
  # * Initializes default model and options
  #* Parameters 
  # * [new_model] : model name
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _init_model(new_model)
    #mathmatical model, default is two-parameter logistic model
    @model = new_model
    
    #range of ability(theta), default is from -4 to 4
    # @range_of_ability = [-4, 4]
    @range_of_ability = [-3, 3]
    
    #number of ability points, default is 40 points 
    @points_of_ability = 40
    
    #D is logistic scaling constant for 3PL
    #(the default value 1.7 makes logistic courve, close to normal) 
    @D = 1.7
    
    #if true, different points of ability(latent distribution points) are used for
    #different examinee groups, default is false
    #SEE : num_groups option of IRT class
    @use_unique_points = false
  end
  
  #* Functionailities
  # * Initializes default prior paramters
  #* Parameters 
  # * N/A
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _init_prior_for_parameters
   #default prior used for dicrimination of item
   @prior_a = {:type => 'beta', :value => [1.75, 3.0, 0.0, 3.0]}
   
   #default prior used for difficulty of item
   @prior_b = {:type => 'beta', :value => [1.01, 1.01, -6.0, 6.0]}
   
   #default prior used for guessing of item
   @prior_c = {:type => 'beta', :value => [3.5, 4.0, 0.0, 0.5]}

    #default prior used for dicrimination of item
    # @prior_a = {:type => 'beta', :value => [5.0, 17.0, 0.0, 1.0]}
    # 
    # #default prior used for difficulty of item
    # @prior_b = {:type => 'lognormal', :value => [0.0, 0.5]}
    # 
    # #default prior used for guessing of item
    # @prior_c = {:type => 'beta', :value => [1.01, 1.01, -6.0, 6.0]}
  end
end