#* irt_exam.rb
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
#* IRTExam is a class to handle data of test form (e.g., test's score, items, examinees, examinees' responses, etc.).
#* Author: KyoungMo Yang(E-mail : mo@embian.com)
#* Modified: 2010-11-05
#* References: N/A
class IRTExam
  attr_accessor :items, :examinees, :number_of_groups
  #constants for examinees' missing response  
  MISSING_RESPONSE    = -1
  #constants for examinees' correct response
  CORRECT_RESPONSE    = 1
  #constants for examinees' incorrect response
  INCORRECT_RESPONSE  = 0
  
  #instance fields
  @items
  @number_of_items
  
  @examinees
  @number_of_examinees
  @number_of_groups
  
  #* Functionailities
  # * Constructs an instance of IRTTset class.
  #* Parameters 
  # * [new_examinees] : array of IRTExaminee
  # * [new_items] : array of IRTItem
  # * [new_number_of_groups] : the number of groups. default value is one.
  #* Return
  # * An instance of IRTExam class
  #* Note
  # * N/A
  def initialize(new_examinees, new_items, new_number_of_groups = 1)
    @examinees = new_examinees
    @items = new_items
    @number_of_groups = new_number_of_groups
    _init_question_no_of_items
  end
  
  #* Functionailities
  # * return the number of examinees.
  #* Parameters 
  # * N/A
  #* Return
  # * [result] : if true, using dichotmomous model otherwise, return false
  #* Note
  # * N/A
  def number_of_examinees
    return @examinees.size
  end

  #* Functionailities
  # * return the number of items.
  #* Parameters 
  # * N/A
  #* Return
  # * [num_items] : the number of items.
  #* Note
  # * N/A
  def number_of_items
    return @items.size
  end
  
  #* Functionailities
  # * return question number at given index
  #* Parameters 
  # * N/A
  #* Return
  # * [no] : question number
  #* Note
  # * N/A
  def get_item_no(index)
    return @items[index].no
  end
  
  private
  #* Functionailities
  # * initialize question No. of items by array order
  #* Parameters 
  # * N/A
  #* Return
  # * N/A
  #* Note
  # * N/A
  def _init_question_no_of_items
    @items.each_with_index do |item, i|
      item.no = i+1
    end
  end
end