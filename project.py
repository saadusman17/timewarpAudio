import time
import random
import matplotlib.pylab as plt
from math import sin, cos
from pysoundfile import SoundFile
from numpy import *
import scipy.io.wavfile as wavreader
from copy import deepcopy
import os
from xml.dom import minidom


def print_matrix(sky):
    if(len(sky) > 0):
        number_of_characters = 1
        print "=" * len(sky[0]) * (number_of_characters * 6 + 1)
        for i in range(0, len(sky)):
            for j in range(0, len(sky[i])):
            	
 	            print "{:6.2f}".format(sky[i][j]), #.zfill(number_of_characters),
            print
        print "=" * len(sky[0]) * (number_of_characters* 6  + 1)


def time_warp(a, b): #standard time warp with no space or time efficiency
#	print a, b
	len_a = len(a)
	len_b = len(b)
	
	if(len_a < len_b):
		temp = a
		temp_l = len_a
		
		a = b
		len_a = len_b
		
		b = temp
		len_b = temp_l
	
	matrix = [ [0]* (len(a) + 2) for i in xrange(len(b) + 2)]

	for i in range(2, len(a) + 2):
		matrix[0][i] = a[i - 2]

	for i in range(2, len(b) + 2):
		matrix[i][0] = b[i - 2]

	for i in range(1, len(matrix)):
		for j in range(1, len(matrix[0])):
			if(i == 1):
				matrix[i][j] = matrix[i][j-1] + abs(matrix[i][0] - matrix[0][j]) / float(2)
			elif(j == 1):
				matrix[i][j] = matrix[i-1][j] + abs(matrix[i][0] - matrix[0][j]) / float(2)	
			else:
				matrix[i][j] = min(matrix[i-1][j-1] + abs(matrix[i][0] - matrix[0][j]),
							   matrix[i-1][j] + abs(matrix[i][0] - matrix[0][j]) / float(2),
							   matrix[i][j-1] + abs(matrix[i][0] - matrix[0][j]) / float(2))
	print_matrix(matrix)
	return matrix[-1][-1]

def time_warp2(a, b): #Time warp which is time efficient, takes only some cells in the middle
	len_a = len(a)
	len_b = len(b)
	
	if(len_a < len_b):
		temp = a
		temp_l = len_a
		
		a = b
		len_a = len_b
		
		b = temp
		len_b = temp_l
	

	matrix = [ [0]* (len(a) + 2) for i in xrange(len(b) + 2)]

	for i in range(2, len(a) + 2):
		matrix[0][i] = a[i - 2]

	for i in range(2, len(b) + 2):
		matrix[i][0] = b[i - 2]
		
	global rows_to_consider
	stop_coming_down = False
	
	for j in range(1, len(matrix[0])):
		start = (len(b) + 2)*(j - 1) /float(len(a)) - rows_to_consider
		start = int(start)
		
		end = start + 2*rows_to_consider
		
		if(start < 1):
			start = 1
			
		if(end > len(matrix) - 1):
			end = len(matrix)
		
		if(start - 1 >= 1):
			matrix[start -1][j] = float("inf")
			
		if(end < len(matrix)):
			matrix[end][j] = float("inf")
		
		for i in range(start,end):
			if(i == 1):
				matrix[i][j] = matrix[i][j-1] + abs(matrix[i][0] - matrix[0][j]) / float(2)
			elif(j == 1):
				matrix[i][j] = matrix[i-1][j] + abs(matrix[i][0] - matrix[0][j]) / float(2)	
			else:
				matrix[i][j] = min(matrix[i-1][j-1] + abs(matrix[i][0] - matrix[0][j]),
							   matrix[i-1][j] + abs(matrix[i][0] - matrix[0][j]) / float(2),
							   matrix[i][j-1] + abs(matrix[i][0] - matrix[0][j]) / float(2))
	print_matrix(matrix)
	return matrix[-1][-1]


def time_warp_space_efficient(a, b, rows_to_consider): #implements space efficiency so rather than O(n*n), it uses just O(n) space
	len_a = len(a)
	len_b = len(b)
	
	if(len_a < len_b):
		temp = a
		temp_l = len_a
		
		a = b
		len_a = len_b
		
		b = temp
		len_b = temp_l
	
	C = 0
	C_new = 0
	m = zeros((len(b) + 1, 1), dtype=float128)
	
	stop_coming_down = False	
	
	#full_matrix = zeros((len(b) + 1, 1))
	
	for j in range(0, len_a + 1):
		start = (len(b) + 2)*(j - 1) /float(len(a)) - rows_to_consider
		start = int(start)
		
		end = start + 2*rows_to_consider
		
		if(start < 0):
			start = 0
			
		if(end > len_b + 1):
			end = len_b + 1
		
		if(start - 1 >= 0):
			m[start -1, 0] = float("inf")
			
		if(end < len_b + 1):
			m[end, 0] = float("inf")

		for i in range(start,end):			
			if(i ==0 and j ==0):
				pass
			elif(i == 0):
				C_new = m[i, 0]
				m[i, 0] = m[i, 0] + abs(a[j - 1]) / float(2)
			elif(j == 0):
				C_new = m[i, 0]
				m[i, 0] = m[i-1, 0] + abs(b[i - 1]) / float(2)	
			else:
				C_new = m[i, 0]
				m[i, 0] = min(C + abs(b[i - 1] - a[j - 1]),
							   m[i-1, 0] + abs(b[i - 1] - a[j - 1]) / float(2),
							   m[i, 0] + abs(b[i - 1] - a[j - 1]) / float(2))				
			C = C_new	
		#print i, j
	print C, m	
		#making whole matrix for debugging purposes
		#full_matrix = append(full_matrix,m,1)
			
	#print_matrix(full_matrix)
	return m[-1, 0]

def clip_empty_area(a, threshold):
	starting_index = 0
	for i in range(0, len(a)):
		#print a[i]
		if(logical_and(a[i] < threshold, a[i] > -threshold)):
			continue
		else:
			starting_index = i
			break;
		
	ending_index = len(a) - 1
	for i in range(len(a) - 1, -1 , -1):
		if(logical_and(a[i] < threshold, a[i] > -threshold)):
			continue
		else:
			ending_index = i
			break;
	return a[starting_index:ending_index]

def plot_audio_signals(data1, data2):
	axes1 = plt.subplot(1,2,1)	
	axes1.plot(data1)
	
	axes2 = plt.subplot(1,2,2)
	axes2.plot(data2)
	plt.show()	

def time_warp_sound_files(file1, file2, word1, word2, threshold):
	start_time = time.time()
	#opening data
	data1 = wavreader.read("eng-wims-mary_flac/flac/" + file1)
	data1=  data1[1].astype('float64')

	data2 = wavreader.read("eng-wims-mary_flac/flac/" + file2)
	data2 = data2[1].astype('float64')

	#axes1 = plt.subplot(2,2,1)	
	#axes1.set_title(word1)
	#axes1.plot(data1)
	
	#axes2 = plt.subplot(2,2,2)
	#axes2.set_title(word2)
	#axes2.plot(data2)

	#Clipping data
	data3 = clip_empty_area (data1, 1000)
	data4 = clip_empty_area(data2, 1000)

	#normalizing
	data3 /= max(abs(data3))
	data4 /= max(abs(data4))
	
	axes1 = plt.subplot(1,2,1)	
	axes1.set_title(word1)
	axes1.plot(data3)
	
	axes2 = plt.subplot(1,2,2)
	axes2.set_title(word2)
	axes2.plot(data4)
	plt.show()

	#Test the situation where table has lot more rows than columns
	#a = array([1,2,-5,2,3,-1,2,1,-1,2,-5,2,3,-1])
	#b = array([2,-3,-9,2,1,2,-2,-3,-9,2,1,2,-2,-3,-9,2,1,2,-2])
	
	#print time_warp(data1, data2)
	#print "first time warp", time_warp(deepcopy(a),deepcopy(b))
	#print "second time warp", time_warp2(a,b)
	#print "data type, shape", type(data1), data1.shape
	#print " a type", type(a), a.shape
	#result = time_warp_space_efficient(data1, data2, threshold)
	result = 0
	##print "final version: ", time_warp_space_efficient(a,b)

	
	#some code
	if(len(data1) > len(data2)):
		result = sum(abs(data1))/float(2)
	else:
		result = sum(abs(data2))/float(2)
	print word1, word2, result		
	return result


def main():
	
	#Debugging Code
	#Test the situation where table has lot more rows than columns
	#a = array([1,2,-5,2,3,-1,2,1,-1,2,-5,2,3,-1])
	#b = array([2,-3,-9,2,1,2,-2,-3,-9,2,1,2,-2,-3,-9,2,1,2,-2,2,55,-2])
	#print "sum_of_b ", sum(abs(b))/float(2)
	
	#print time_warp(a,b)
	#print "first time warp", time_warp(deepcopy(a),deepcopy(b))
	#print "second time warp", time_warp2(a,b)
	#print "final version: ", time_warp_space_efficient(a,b, rows_to_consider)
	#exit()
	
	#End of Debugging Code

	word_to_file_map = []

	xmldoc = minidom.parse('eng-wims-mary_flac/flac/index.xml')
	itemlist = xmldoc.getElementsByTagName('file') 
	print len(itemlist)
	#print itemlist[0].attributes['name'].value
	for s in itemlist :
		filename = s.attributes['path'].value
		
		#changing extension to wav
		root,ext = os.path.splitext(filename)
		filename = root + ".wav"
		
		tag = s.getElementsByTagName('tag')
		word = tag[0].attributes['swac_text'].value
		word_to_file_map.append([filename, word])
	
	#print word_to_file_map
	
	#restricting for scalability
	word_to_file_map = word_to_file_map[10:20]
	
	output_fh = open("output", "w")
	output_fh.write("File 1, File 2, Word 1, Word 2, Time Warp Value\n")
	time_warp_sound_files("eng-043252a4.wav" , "eng-05cda17c.wav",  "bizonal","distinctness", 500)
	
	for i in range(0, len(word_to_file_map)):
		for j in range(i + 1, len(word_to_file_map)):
			time_warp_result = time_warp_sound_files(word_to_file_map[i][0] , word_to_file_map[j][0],  word_to_file_map[i][1], word_to_file_map[j][1], 500)
			output_fh.write(word_to_file_map[i][0] + "," + word_to_file_map[j][0] + "," + word_to_file_map[i][1] + "," + word_to_file_map[j][1] + "," + str(time_warp_result) + "\n")
		
	output_fh.close()

	

rows_to_consider = 5

if __name__ == '__main__':
	main()
