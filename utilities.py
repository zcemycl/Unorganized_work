import numpy as np
import glob
import cv2
class UtilityFunctions:
	@staticmethod
	def getImageFileList(file_name_string):
        	result = sorted(glob.glob(file_name_string))
        	return result
	@staticmethod
	def getImages(file_name_string):
		list=UtilityFunctions.getImageFileList(file_name_string)
		result=[]
		for fname in list:
			img=cv2.imread(fname)
  			gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
			result.append(gray)
		return result

	@staticmethod
	def writePlyFile(file_name, vertices, colors):
		ply_header = '''ply
				format ascii 1.0
				element vertex %(vert_num)d
				property float x
				property float y
				property float z
				property uchar red
				property uchar green
				property uchar blue
				end_header
		     	     '''
		vertices = vertices.reshape(-1, 3)
		colors = colors.reshape(-1, 3)
		vertices = np.hstack([vertices, colors])
		with open(file_name, 'w') as f:
			f.write(ply_header % dict(vert_num=len(vertices)))
			np.savetxt(f, vertices, '%f %f %f %d %d %d')

