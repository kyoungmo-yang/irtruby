require 'test/unit'
require '../irt_examinee'

#* irt_examinee_test.rb
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
class IRTItemTest < Test::Unit::TestCase
  def setup
    @responses = [1, 1, 1, 0, -1, 0, 1, 1]
    @no = 1
    @group = 1
  end
  
  def teardown
    #Do nothing
  end
  
  def test_instance_methods_of_examinee_class
    examinee = IRTExaminee.new("id##{@no}", @responses, "810000-00000#{@no}", @group)
    assert_not_nil(examinee)
    assert_equal("id#1", examinee.eid)
    assert_equal(@responses, examinee.responses)
    assert_equal("810000-00000#{@no}", examinee.social_no)
    assert_equal(@group, examinee.group)
    assert_equal(@no, examinee.no)
    
    examinee.theta = 1.15
    assert_equal(1.15, examinee.theta)
    
    examinee.no = 2
    assert_equal(5, examinee.ctt_score)
  end
end