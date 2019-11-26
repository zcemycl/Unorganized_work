# The python file for Exercise 1 - Calibration

# The goal of this exercise is to build a framework for camera intrinsic and extrinsic parameter calibration from a sequence of calibration images. 
# This code then is used to display a virtual object (3D cube) on a calibrated chess-board pattern. 
#
# Note that in order to complete an exercise several lines need to be inserted into this code. Places where such lines need to be added are marked as follows:
#
# Task: Read image file 'images-ex01/calibration/calibration-01.png'
# [--------LINE-TO-BE-INSERTED
# 
# --------]
#
# A successful insertion would look like this:

#
import numpy as np
import cv2
import glob
import os
from PIL import Image
from utilities import UtilityFunctions #importing functions from utilities.py file provided with the exercise, make sure to copy it in the same location together with cv_ex01_calibration.py


# Function performing drawing of a projected 3D cube on the provided image
def draw(image, points):
    points = np.int32(points).reshape(-1,2)

    # draw ground floor in green
    image = cv2.drawContours(image, [points[:4]],-1,(0,255,0),-3)

    # draw pillars in blue color
    for i,j in zip(range(4),range(4,8)):
        image = cv2.line(image, tuple(points[i]), tuple(points[j]),(255),3)

    # draw top layer in red color
    image = cv2.drawContours(image, [points[4:]],-1,(0,0,255),3)

    return image

#Function finding chess-board calibration pattern in a provided image
def findChessBoardCorners(image, corner_size):

    # Termination criteria
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)

    # Prepare calibration pattern points (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
    object_points = np.zeros((corner_size[0]*corner_size[1],3),np.float32)
    object_points[:,:2] = np.mgrid[0:corner_size[0],0:corner_size[1]].T.reshape(-1,2)
# ------------------------------------------------------------------------------------- #
    # [--------LINE-TO-BE-INSERTED
    # Task: find the chess board corners using OpenCV function findChessboardCorners
    retval, corners = cv2.findChessboardCorners(image, corner_size, None)    
    
    # --------]
    
# ------------------------------------------------------------------------------------- #

    # If pattern detected, refine corners
    if retval == True:
 
        refined_corners = cv2.cornerSubPix(image,corners, corner_size, (-1,-1), criteria)##
	#print(refined_corners)
	return True, object_points, refined_corners
    else:

        return False, None, None

# Function handling finding a calibration pattern in an image (currently - chessboard)
def processImage(image, corner_size, object_point_list, image_point_list):

    # Can replace this in order to use a different calibration pattern
    retval, object_points, image_points = findChessBoardCorners(image, corner_size)

    if (retval):

	object_point_list.append(object_points)
        image_point_list.append(image_points)

    return retval

# Function drawing chess board pattern found onto an image
def drawCorners(image, pattern_size, corners, retval):
        # Draw and display the corners
        img = cv2.drawChessboardCorners(image, pattern_size, corners, retval)

# Function finding the corners of chess board pattern in a list of images
def findCorners(images, pattern_size):
	# Arrays to store object points and image points from all the images.
	object_points = [] # 3d point in real world space
	image_points = [] # 2d points in image plane.

	for id,gray in enumerate(images):

		retval=processImage(gray, pattern_size, object_points, image_points)

		if(retval): #for visualisation purposes only
			drawCorners(gray, pattern_size, image_points[-1], retval)

	        cv2.imwrite('images-ex01/results/calibration_'+str(id).zfill(5)+'_LINES.png',gray)
	return object_points, image_points

# Function calculating reprojection error of a calibrated chess-board pattern
def getReprojectionError(object_points, image_points, camera_matrix, distortion_parameters, rotation_vectors, translation_vectors):
        mean_error = 0
	for i in range(len(object_points)):

# ------------------------------------------------------------------------------------- #

		# [--------LINE-TO-BE-INSERTED
		# Task: Find projection of calibration pattern OpenCV function projectPoints
                projected_points, _ = cv2.projectPoints(object_points[i], rotation_vectors[i], translation_vectors[i], camera_matrix, distortion_parameters)
		# ---------]

# ------------------------------------------------------------------------------------- #


		error = cv2.norm(image_points[i], projected_points, cv2.NORM_L2) / len(projected_points)
		mean_error += error
	mean_error/=len(object_points)
	return mean_error

# Function performing camera calibration from a list of images
def calibrateCamera(images, pattern_size, image_size):
	print('Calibrating camera ...')
	object_points,image_points=findCorners(images, pattern_size)

	# [--------LINE-TO-BE-INSERTED
	# Task: perform camera calibration with OpenCV function calibrateCamera 
        retval, camera_matrix, distortion_parameters, rotation_vectors, translation_vectors = cv2.calibrateCamera(object_points, image_points, image_size[::-1],None,None)
        
	# ============================================ 
	#print(rotation_vectors) #32
    	#print(type(rotation_vectors))
	#print(len(rotation_vectors))
	#rotation_matrix = cv2.Rodrigues(rotation_vectors[0].T)
	#print(rotation_matrix)
	# ============================================ #

	# Comparing reprojection error returned by camera calibration code and by your evaluation
	reprojection_error = getReprojectionError(object_points, image_points, camera_matrix, distortion_parameters, rotation_vectors, translation_vectors)

	print('Reprojection error  = '+str(reprojection_error))

	print('Calibrating camera - finished')
	return reprojection_error, camera_matrix, distortion_parameters, rotation_vectors, translation_vectors

# Function displaying a virtual object on an image list from given intrinsic camera parameters
def displayPattern(image_list, result_dir, corner_size, camera_matrix, distortion_parameters):
	print('Calculating virtual reality patterns ...')

	cube_coordinates = np.float32([[0,0,0], [0,3,0], [3,3,0], [3,0,0],
                   [0,0,-3],[0,3,-3],[3,3,-3],[3,0,-3] ])

	for fname in image_list:
		print('Processing image: '+fname)
		image = cv2.imread(fname)
		gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

		retval, object_points, image_points = findChessBoardCorners(gray_image, corner_size)
		# ============================================ #
		#print(object_points) #56
		#print(image_points)
		# ============================================ #
    		if retval == True:

			object_points_Nx1x3 = np.zeros((object_points.shape[0], 1, 3), np.float32)
			object_points_Nx1x3[:, 0, :] = object_points[:, :]
			
# ------------------------------------------------------------------------------------- #
			# [--------LINE-TO-BE-INSERTED
			# Task: find the rotation and translation vectors using OpenCV function solvePnPRansac. Inputs include object points, image points, camera matrix and distortion parameters
                        retval, rvecs, tvecs = cv2.solvePnP(object_points_Nx1x3, image_points, camera_matrix, distortion_parameters)
			# ============================================ #
			#print(camera_matrix)
			#print(distortion_parameters)
			print(rvecs)
			
			# ============================================ #
			# ---------]
			
# ------------------------------------------------------------------------------------- #


			# Projecting 3D points to image plane
        		projected_cube_points, jac = cv2.projectPoints(cube_coordinates, rvecs, tvecs, camera_matrix, distortion_parameters)

			image = draw(image, projected_cube_points)
			#print(projected_cube_points[0,0])
			#print(projected_cube_points[4,0])
			#print(type(projected_cube_points))
			#print(projected_cube_points[0,0][0])	
		# ============================================ #
			#project_specific, jj = cv2.projectPoints(np.float32([[0,0,-3]]),rvecs,tvecs, camera_matrix, distortion_parameters)
			#image = cv2.circle(image, (projected_cube_points[0,0][0], projected_cube_points[0,0][1]), 10, (255,0,0), -1)
			#image = cv2.circle(image, (projected_cube_points[4,0][0], projected_cube_points[4,0][1]), 10, (255,0,0), -1)
		# ============================================ #
			cv2.imwrite(result_dir+fname.split('/')[-1].split('.')[0]+'.png', image)

	print('Calculating virtual reality patterns - finished')

def main():
	calibration_image_list = UtilityFunctions.getImages('images-ex01/calibration/*.png')

	reprojection_error, camera_matrix, distortion_parameters, rotation_vectors, translation_vectors = calibrateCamera(calibration_image_list,(9,7),(1920, 1080))

	vrpattern_image_list = UtilityFunctions.getImageFileList('images-ex01/ar-images/*.png')

	displayPattern(vrpattern_image_list,'images-ex01/results/',(9,7),camera_matrix, distortion_parameters)

	print("Finished: cv-ex01-calibration.py")

if __name__ == "__main__":
   main()

