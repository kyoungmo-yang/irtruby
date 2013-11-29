require 'test/unit'
require '../irt_item'

#* irt_item_test.rb
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
    #Do nothing
  end
  
  def teardown
    #Do nothing
  end
  
  def test_instance_methods_of_item_class
    item = IRTItem.new("id##{1}")
    assert_not_nil(item)
    assert_equal("id#1", item.iid)
    assert_equal(0, item.no)
    
    item.irt_a = 0.159
    assert_equal(0.159, item.irt_a)
    
    item.irt_b = 1.222
    assert_equal(1.222, item.irt_b)
    
    item.irt_c = 0.000
    assert_equal(0.000, item.irt_c)
    
    item.answer = 'A'
    assert_equal('A', item.answer)
  end
end