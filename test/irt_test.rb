require 'test/unit'
require '../irt'
require '../irt_item'
require '../irt_examinee'
require '../irt_exam'
require '../irt_model'
require '../irt_algorithm'

#* irt_test.rb
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

class IRTTest < Test::Unit::TestCase
  def setup
    @examinees = []
    @items = []

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

    @expected_score = 
    [
      4.4522617514,
      3.3079929514,
      4.2417963495,
      2.3644541704,
      4.4522617514,
      3.6256603344,
      2.9618310621,
      4.4522617514,
      4.1226182523,
      1.7764694147
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
      @items << item
    end

    @examinees = []
    10.times.each do |no|
      @examinees << IRTExaminee.new("id##{no}", @responses[no], "00000-0000#{no}")
    end

    @exam = IRTExam.new(@examinees, @items)
    @irt = IRT.new(@exam, IRTModel.new(IRTModel::DICHOTOMOUS_3PL))
  end
  
  def teardown
    @irt.release
  end

  def test_estimates_item_params
    @items.each do |item|
      assert_equal(0.0, item.irt_a)
      assert_equal(0.0, item.irt_b)
      assert_equal(0.0, item.irt_c)
    end
    @irt.estimates_item_params

    @items.each_with_index do |item, i|
      params = @expected_item_params[i]
      assert_equal(params[:irt_a], ("%.10f" % item.irt_a).to_f)
      assert_equal(params[:irt_b], ("%.10f" % item.irt_b).to_f)
      assert_equal(params[:irt_c], ("%.10f" % item.irt_c).to_f)
    end
  end

  def test_estimates_examinee_theta
    @irt.estimates_item_params
    @examinees.each do |examinee|
      assert_equal(0.0, examinee.theta)
    end

    @irt.estimates_examinee_theta

    @examinees.each_with_index do |examinee, i|
      assert_equal(@expected_thetas[i], ("%.10f" % examinee.theta).to_f)
    end
  end

  def test_estimates_examinee_score
    @irt.estimates_item_params
    @irt.estimates_examinee_theta
    @examinees.each do |examinee|
      assert_equal(0.0, examinee.score)
      assert(examinee.scores.empty?)
    end

    @irt.estimates_examinee_score(true)

    @examinees.each_with_index do |examinee, i|
      assert_equal(@expected_score[i], ("%.10f" % examinee.score).to_f)
      @items.each do |item|
        assert_not_nil(examinee.scores[item.no])
      end
    end
  end

  def test_load_estimated_item_params
    @irt.estimates_item_params
    @irt.estimates_examinee_theta
    @irt.estimates_examinee_score(true)
    
    #----------------------------------
    
    estimated_items = []
    5.times.each do |no|
      item = IRTItem.new("id##{no}")
      item.answer = @items[no].answer
  
      item.irt_a = @items[no].irt_a
      item.irt_b = @items[no].irt_b
      item.irt_c = @items[no].irt_c
  
      estimated_items << item
    end
    
    new_examinees = []
    10.times.each do |no|
      new_examinees << IRTExaminee.new("id##{no}", @responses[no], "00000-0000#{no}")
    end
    
    test = IRTExam.new(new_examinees, estimated_items)
    irt = IRT.new(test, @model_3PL)
    irt.load_estimated_item_params
    
    estimated_items.each_with_index do |item, i|
      params = @expected_item_params[i]
      assert_equal(params[:irt_a], ("%.10f" % item.irt_a).to_f)
      assert_equal(params[:irt_b], ("%.10f" % item.irt_b).to_f)
      assert_equal(params[:irt_c], ("%.10f" % item.irt_c).to_f)
    end
  
    irt.estimates_examinee_theta
  
    #not equals pretested result 
    expected_thetas = 
    [    
      0.8754095948,
      -0.3143090248,
      0.4527003898,
      -0.8396082799,
      0.8754095948,
      -0.1153879772,
      -0.4960945621,
      0.8754095948,
      0.3704384356,
      -1.549457165
    ]
  
    new_examinees.each_with_index do |examinee, i|
      assert_equal(expected_thetas[i], ("%.10f" % examinee.theta).to_f)
    end
  
    #not equals pretested result 
    expected_score  = 
    [    
      4.4831201496,
      3.2789321125,
      4.1620558427,
      2.5516283623,
      4.4831201496,
      3.5399176403,
      3.0285407549,
      4.4831201496,
      4.0851422551,
      1.7654646958
    ]
     
    irt.estimates_examinee_score(true)
    new_examinees.each_with_index do |examinee, i|
      assert_equal(expected_score[i], ("%.10f" % examinee.score).to_f)
      @items.each do |item|
        assert_not_nil(examinee.scores[item.no])
      end
    end
    
    irt.release
  end
end