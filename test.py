# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 23:33:26 2014

@author: pascal
"""

#nested_list = [['x', 'y', 'z'], ['a', 'b', 'c']].
nested_list = [['x', 'y', 'z'], ['a', 'b', 'c']]
for line in nested_list:
	print '\t'.join(line)



nested_list = [[0.03434355, 1.03434355, 0.13434355], [0.03434355, 0.03434355, 34.03434355]]
for line in nested_list:
	print '\t'.join("%.3f".format(line))
