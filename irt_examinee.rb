#* irt_examinee.rb
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
#* The IRTExaminee is an class to handle examinee's data.
#* Author: KyoungMo Yang(E-mail : mo@embian.com)
#* Modified: 2010-11-05
#* References: N/A
class IRTExaminee
  attr_accessor :eid, :no, :theta, :scores, :score, :group, :responses, :social_no
  
  #instance fields
  @eid
  @no
  @theta
  @score
  @scores
  @ctt_score
  @responses
  @social_no
  @group
  
  #* Functionailities
  # * Constructs an instance of IRTExaminee.
  #* Parameters 
  # * [new_id] : id of examinee.
  # * [new_responses] : responses of examinee.
  # * [social_no] : social number of examinee. default value is NULL.
  # * [new_group] : group number(index). default value is one.
  #* Return
  # * An instance of IRTExaminee class
  #* Note
  # * N/A
  def initialize(new_id, new_responses, social_no = nil, new_group = 1)
    @eid = new_id
    @responses = new_responses
    @scores = []
    @group = new_group
    @social_no = social_no
    @ctt_score = 0
    @no = 1
    @theta = 0.0
    @score = 0
  end
  
  #* Functionailities
  # * Calculates CTT Score for an examinee
  #* Parameters 
  # * N/A
  #* Return
  # * CTT Score
  #* Note
  # * This method is used in Dichotomous Model
  def ctt_score
    ctt_score = @responses.select{|v| v.to_i > 0 }.inject{|sum, v2| sum+=v2.to_i}
    @ctt_score = ctt_score
    return @ctt_score
  end
end