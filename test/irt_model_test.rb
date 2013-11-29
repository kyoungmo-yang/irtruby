require 'test/unit'
require '../irt_model'

#* irt_model_test.rb
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
class IRTModelTest < Test::Unit::TestCase
  def setup
    @range_of_ability = [-4, 4]
    @points_of_ability = 30
    @D = 1.7
    @prior_a = {:type => 'beta', :value => [1.75, 3.0, 0.0, 3.0]}
    @prior_b = {:type => 'lognormal', :value => [0.0, 0.5]}
    @prior_c = {:type => 'normal', :value => [0.5, 0.5]}
    @use_unique_points = false
  end
  
  def teardown
    #Do nothing
  end
  
  def test_instance_methods_of_item_class
    irt_model = IRTModel.new(IRTModel::DICHOTOMOUS_3PL)
    assert_not_nil(irt_model)
    assert_equal(IRTModel::DICHOTOMOUS_3PL, irt_model.model)
    
    irt_model.range_of_ability = @range_of_ability
    assert_equal(@range_of_ability, irt_model.range_of_ability)
    
    irt_model.points_of_ability = @points_of_ability
    assert_equal(@points_of_ability, irt_model.points_of_ability)
    
    irt_model.D = @D
    assert_equal(@D, irt_model.D)
    
    irt_model.prior_a = @prior_a
    assert_equal(@prior_a, irt_model.prior_a)
    
    irt_model.prior_b = @prior_b
    assert_equal(@prior_b, irt_model.prior_b)
    
    irt_model.prior_c = @prior_c
    assert_equal(@prior_c, irt_model.prior_c)
    
    irt_model.use_unique_points = @use_unique_points
    assert_equal(@use_unique_points, irt_model.use_unique_points)
    
    assert(irt_model.use_dichotomous_model?)
  end
end