require 'test/unit'
require '../lib/irt_ruby'
require '../irt_algorithm'
require '../irt_item'
require '../irt_examinee'

#* irt_algorithm_test.rb
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
class IRTAlgorithmTest < Test::Unit::TestCase
  def setup
    @expected_item_params = [
      {:irt_a=>1.1589553583, :irt_b=>-0.9637488867, :irt_c=>0.2208181695},
      {:irt_a=>0.5245147546, :irt_b=>-0.5362510839, :irt_c=>0.2307967442},
      {:irt_a=>1.0403761972, :irt_b=>0.2636368398, :irt_c=>0.2253940896},
      {:irt_a=>0.9649196229, :irt_b=>-0.2708693551, :irt_c=>0.2267874693},
      {:irt_a=>1.2080574422, :irt_b=>-0.7281504136, :irt_c=>0.2158702943},
    ]

    @expected_thetas = 
    [
      0.8266582137,
      -0.2927585773,
      0.5439140479,
      -0.9818123114,
      0.8266582137,
      -0.0466947116,
      -0.5437991172,
      0.8266582137,
      0.4098901131,
      -1.5359306825,
    ]

    answers = ['A', 'B', 'C', 'D']
    @responses = [
      [1, 1, 1, 1, 1],
      [1, 1, 0, 0, 1],
      [1, 0, 1, 1, 1],
      [0, 1, 1, 1, 0],
      [1, 1, 1, 1, 1],
      [1, 0, 0, 1, 1],
      [1, 1, 0, 0, -1],
      [1, 1, 1, 1, 1],
      [1, 1, -1, -1, -1],
      [0, 0, 0, 0, 0]
    ]
    
    @items = []
    5.times.each do |no|
      item = IRTItem.new("id##{no}")
      item.answer = answers[rand(answers.size)]
      item.no = no+1
      @items << item
    end

    @examinees = []
    10.times.each do |no|
      @examinees << IRTExaminee.new("id##{no}", @responses[no], "00000-0000#{no}")
    end
    
    #initialize
    @irt_ruby = Irt_ruby::IRTModule.new
    @irt_ruby.new_items_dist(@expected_item_params.size, 40, 1, [1, 1, 1, 1, 1], [-3, 3], false)
    @examinees.each do |examinee|
      no = (@irt_ruby.examinees_count + 1).to_i
      @irt_ruby.add_examinee(examinee.responses, examinee.group)
      examinee.no = no
    end
    
    @algorithm = IRTAlgorithm.new(@irt_ruby)
  end
  
  def teardown
    @irt_ruby.delete_items_dist
  end
  
  def test_estimates_item_params
    @items.each do |item|
      assert_equal(0.0, item.irt_a)
      assert_equal(0.0, item.irt_b)
      assert_equal(0.0, item.irt_c)
    end
    
    use_dichotomous_model=true
    @algorithm.estimates_item_params(@items, use_dichotomous_model)
    
    @items.each_with_index do |item, i|
      params = @expected_item_params[i]
      assert_equal(params[:irt_a], ("%.10f" % item.irt_a).to_f)
      assert_equal(params[:irt_b], ("%.10f" % item.irt_b).to_f)
      assert_equal(params[:irt_c], ("%.10f" % item.irt_c).to_f)
    end
  end
  
  def test_estimates_examinee_theta
    use_dichotomous_model=true
    @algorithm.estimates_item_params(@items, use_dichotomous_model)
    
    @examinees.each do |examinee|
      assert_equal(0.0, examinee.theta)
    end
    
    @algorithm.estimates_examinee_theta(@examinees, [-4, 4])
    
    @examinees.each_with_index do |examinee, i|
      assert_equal(@expected_thetas[i], ("%.10f" % examinee.theta).to_f)      
    end
  end
end