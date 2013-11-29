require 'test/unit'
require '../irt_item'
require '../irt_examinee'
require '../irt_exam'

#* irt_exam_test.rb
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
class IRTExamTest < Test::Unit::TestCase
  def setup
    answers = ['A', 'B', 'C', 'D']
    @items = []
    5.times.each do |no|
      item = IRTItem.new("id##{no}")
      item.answer = answers[rand(answers.size)]
      @items << item
    end
    
    @examinees = []
    10.times.each do |no|
      responses = [rand(2), rand(2), rand(2), rand(2), rand(2)]
      responses[rand(responses.size)] = -1 if no % 2 == 0
      @examinees << IRTExaminee.new("id##{no}", responses, "00000-0000#{no}")
    end
  end
  
  def teardown
    #Do nothing
  end
  
  def test_instance_methods_of_item_class
    exam = IRTExam.new(@examinees, @items)
    assert_not_nil(exam)
    
    @items.each_with_index do |item, i|
      assert_equal(i+1, item.no)
    end
    
    assert_equal(@examinees, exam.examinees)
    assert_equal(@items, exam.items)
    
    assert_equal(1, exam.number_of_groups)
    assert_equal(10, exam.number_of_examinees)
    assert_equal(5, exam.number_of_items)
    
    assert_equal(@items[2].no, exam.get_item_no(2))
  end
end