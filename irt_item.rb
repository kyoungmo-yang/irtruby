#* irt_item.rb
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
#* IRTItem is class to handle item's data. (i.g., difficulty - irt_b, descriment - irt_a, guessing - irt_c)
#* Author: KyoungMo Yang(E-mail : mo@embian.com)
#* Modified: 2010-11-05
#* References: N/A
class IRTItem
  #item paramters
  attr_accessor :iid, :no, :irt_a, :irt_b, :irt_c, :answer
  
  #instance fields
  @iid #Item ID
  @no #question No.
  @irt_a  #discriment
  @irt_b  #difficulty
  @irt_c  #guessing
  
  #* Functionailities
  # * Constructs an instance of IRT class.
  #* Parameters 
  # * [new_id] : item id 
  # * [q_no] : question No.
  #* Return
  # * An instance of IRTItem class
  #* Note
  # * N/A
  def initialize(new_id, irt_a=0.0, irt_b=0.0, irt_c=0.0)
    @iid = new_id
    @no = 0
    @irt_a = irt_a
    @irt_b = irt_b
    @irt_c = irt_c
  end
end