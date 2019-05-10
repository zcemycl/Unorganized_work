#Title: The python file for Exercise 3 - 3D Reconstruction
#
#Goal: The goal of this exercise is to build a framework for 3D reconstruction from disparity map with and witout explicit stereo rectification.
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
# [--------LINE-TO-BE-INSERTED
# image = cv2.imread('images-ex01/calibration/calibration-01.png')
# --------]

import numpy as np
import os
import pickle
import cv2
from matplotlib import pyplot as plt
from cv_ex01_calibration import findChessBoardCorners, calibrateCamera,displayPattern 
# importing functions from cv_ex01_calibration.py file provided with the exercise, make sure to copy a working version of it in the
										      # same location together with cv_ex03_reconstruction.py
import cv2
from utilities import UtilityFunctions # importing functions from utilities.py file provided with the exercise, make sure to copy it in the same location together with cv_ex03_reconstruction.py
import math

# Function performing stereo calibration using a list of image pairs
def stereoCalibrate(image_pairs,w,h):
	pattern=(9,7)
	assert(len(image_pairs[0])==len(image_pairs[1])), "Unequal number of images"
	stereocalib_criteria = (cv2.TERM_CRITERIA_MAX_ITER + cv2.TERM_CRITERIA_EPS, 1000, 1e-5)
        stereocalib_flags = cv2.CALIB_FIX_INTRINSIC | cv2.CALIB_USE_INTRINSIC_GUESS | cv2.CALIB_FIX_FOCAL_LENGTH | cv2.CALIB_ZERO_TANGENT_DIST | cv2.CALIB_FIX_K3 | cv2.CALIB_FIX_K4 | cv2.CALIB_FIX_K5

	object_points = [[],[]]
	image_points = [[],[]]
	images=[[],[]]

	# Extracting chessboard pattern from left and right camera images
	for image_id in range(len(image_pairs[0])):
		print('Extracting corners for image id '+str(image_id))
        	ret_L,obj_points_L,img_points_L=findChessBoardCorners(image_pairs[0][image_id],pattern)
		ret_R,obj_points_R,img_points_R=findChessBoardCorners(image_pairs[1][image_id],pattern)

	        if(ret_L==False) or (ret_R==False):
			continue
		assert(np.array_equal(obj_points_L,obj_points_R)),"Wrong object points"

		images[0].append(image_pairs[0][image_id])
		images[1].append(image_pairs[1][image_id])

	        object_points[0].append(obj_points_L)
        	image_points[0].append(img_points_L)
		object_points[1].append(obj_points_R)
        	image_points[1].append(img_points_R)


	assert(len(object_points[0])==len(object_points[1])), "Wrong number of points"
	assert(len(image_points[0])==len(image_points[1])), "Wrong number of points"
	assert(len(image_points[0])==len(object_points[0])), "Wrong number of points"

	# Calibrating left and right camera images
	retL, mtxL, distL, rvecsL, tvecsL = cv2.calibrateCamera(object_points[0], image_points[0], (w,h),None,None)
	retR, mtxR, distR, rvecsR, tvecsR = cv2.calibrateCamera(object_points[1], image_points[1], (w,h),None,None)

	# Overlaying virtual cube images onto individually calibrated left and right  camera images
        displayPattern(UtilityFunctions.getImageFileList('images-ex03/stereo-calibration/*l.png'),'images-ex03/results/',pattern,mtxL,distL)
	displayPattern(UtilityFunctions.getImageFileList('images-ex03/stereo-calibration/*r.png'),'images-ex03/results/',pattern,mtxL,distL)

	print('Calibration finished')
	cameraMatrix1 = mtxL.copy()
	cameraMatrix2 = mtxR.copy()
	distCoeffs1 = distL.copy()
	distCoeffs2 = distR.copy()
	R =None
	T = None
	E = None
	F = None
# -------------------------------------------------------------------------------------------- #
	# Task(1): Perform stereo calibration using OpenCV function stereoCalibrate.
	# [--------LINE-TO-BE-INSERTED
	## obj,img?
	retval,cameraMatrix1,distCoeffs1,cameraMatrix2,distCoeffs2,R,T,E,F=cv2.stereoCalibrate(object_points[0],image_points[0],image_points[1],cameraMatrix1,distCoeffs1,cameraMatrix2,distCoeffs2,R,T,E,F,pattern,criteria = stereocalib_criteria, flags = stereocalib_flags)
	# --------]
# -------------------------------------------------------------------------------------------- #
	print("Stereo calibration finished = "+str(retval))

	return {'mtx1':cameraMatrix1,'dist1':distCoeffs1,'mtx2':cameraMatrix2,'dist2':distCoeffs2, 'R':R, 'T':T,'E': E,'F': F}

# Function performing stereo rectification on a left and right image given camera parameters
def rectifyStereoPair(image_prefix, camera, left_image, right_image):

	w = left_image.shape[1]
	h = right_image.shape[0]
# -------------------------------------------------------------------------------------------- #
        # Task(2): Obtain undistortion and rectification maps for left and right images using OpenCV function stereoRectify
        # [--------LINE-TO-BE-INSERTED
	## camera['image']? alpha?
	R1,R2,P1,P2,Q,roi1,roi2 = cv2.stereoRectify(camera['mtx1'],camera['dist1'],camera['mtx2'],camera['dist2'],(w,h),camera['R'],camera['T'],alpha=-1)
	# --------]
# -------------------------------------------------------------------------------------------- #
	left_maps = cv2.initUndistortRectifyMap(camera['mtx1'], camera['dist1'], R1, P1, (w,h), cv2.CV_16SC2)
	right_maps = cv2.initUndistortRectifyMap(camera['mtx2'],camera['dist2'], R2, P2, (w,h), cv2.CV_16SC2)

	left_image_remap = cv2.remap(left_image, left_maps[0], left_maps[1], cv2.INTER_LANCZOS4)
	right_image_remap = cv2.remap(right_image, right_maps[0], right_maps[1], cv2.INTER_LANCZOS4)

	cv2.imwrite('images-ex03/results/'+image_prefix+'_l_rectified.png',left_image_remap)
	cv2.imwrite('images-ex03/results/'+image_prefix+'_r_rectified.png',right_image_remap)

	return left_image_remap, right_image_remap

# Function calculating disparity map for a left and right image pair
def calculateDisparity(left_image,right_image):
        disparity_size= 32
        block_size= 15

        stereo = cv2.StereoBM_create(numDisparities = disparity_size, blockSize = block_size)
# -------------------------------------------------------------------------------------------- #
        # Task(3): Compute disparity map from left and right image
        # [--------LINE-TO-BE-INSERTED
       	disparity = stereo.compute(left_image,right_image)
	# --------]
# -------------------------------------------------------------------------------------------- #
        disparity = disparity.astype(np.float32)/(disparity_size*1.0)

	return disparity

# Function creating a perspective transformation matrix for the visualistion of the disparity map as a point cloud. Note that estimated Q could be used for images which are stereo-rectified within this code
def obtainPerspectiveTransformationMatrix(f,w,h):
	f = 0.8*w
        Q = np.float32([[1, 0, 0, -0.5*w],
                        [0,-1, 0,  0.5*h],
                        [0, 0, 0,     -f],
                        [0, 0, 1,      0]])
	return Q

# Function generating a 3D point cloud from a disparity map
def generate3DPointsFromDisparity(disparity, left_image):

	Q = obtainPerspectiveTransformationMatrix(0.8,left_image.shape[1],left_image.shape[0])
# -------------------------------------------------------------------------------------------- #
	# Task(4): Calculate 3D points from disparity and projection matrix using OpenCV function reprojectImageTo3D
        # [--------LINE-TO-BE-INSERTED
	points = cv2.reprojectImageTo3D(disparity,Q)
	# --------]
# -------------------------------------------------------------------------------------------- #
        colors = cv2.cvtColor(left_image,cv2.COLOR_GRAY2RGB)

	# Filtering out  low disparity values for cleaner visualisation purposes
        mask = (disparity > 2)

        out_points = points[mask]
        out_colors = colors[mask]

        sum = abs(out_points[:,0])+ abs(out_points[:,1])+abs(out_points[:,2])
        mask2 = ~((sum == float('Inf'))|(np.isnan(sum)))

	# Filtering out points with infinte or NaN coordinates
        points = out_points[mask2,:]
        colors = out_colors[mask2,:]

	return points, colors

# Function performing 3D point cloud reconstruction from rectified image pair
def reconstruct3DPointCloudFromRectifiedImages(image_prefix,left_image,right_image):

	disparity = calculateDisparity(left_image, right_image)

	cv2.imwrite('images-ex03/results/'+image_prefix+'_disparity.png',(disparity/16.0)*255.0)

	points, colors = generate3DPointsFromDisparity(disparity, left_image)

	UtilityFunctions.writePlyFile('images-ex03/results/'+image_prefix+'_3dpoints.ply',points, colors)

# Function loading (or recalculating) stereo calibration parameters
def obtainStereoCameraInformation():
	if os.path.isfile('images-ex03/results/stereo_camera.pic'):
        	camera = pickle.load( open( "images-ex03/results/stereo_camera.pic", "rb" ) )
        else:
        	camera  = stereoCalibrate([UtilityFunctions.getImages('images-ex03/stereo-calibration/*_l.png'), UtilityFunctions.getImages('images-ex03/stereo-calibration/*_r.png')],1280,720)
                pickle.dump( camera, open( "images-ex03/results/stereo_camera.pic", "wb" ) )
	return camera

# Function performing  3D point cloud reconstruction with or without stereo-rectification
def performReconstruction(image_prefix, rectify, left_image, right_image):
	Q = None
	if(rectify == True):
	     	camera=obtainStereoCameraInformation()
                left_image, right_image = rectifyStereoPair(image_prefix, camera, left_image, right_image)
                left_image, right_image = cv2.pyrDown(cv2.pyrDown(left_image)),cv2.pyrDown(cv2.pyrDown(right_image))

        reconstruct3DPointCloudFromRectifiedImages(image_prefix, left_image,right_image)

def main():
	experiment_info=[{'name':'internal', 'rectify':True, 'left':UtilityFunctions.getImages('images-ex03/reconstruction/internal/*_l.png'), 'right':UtilityFunctions.getImages('images-ex03/reconstruction/internal/*_r.png')},
        	        {'name':'middlebury','rectify':False,'left':UtilityFunctions.getImages('images-ex03/reconstruction/middlebury/*1.ppm'), 'right':UtilityFunctions.getImages('images-ex03/reconstruction/middlebury/*2.ppm')}]

	for image_info in experiment_info:
		for i in range(len(image_info['left'])):
			performReconstruction(image_info['name']+'_'+str(i).zfill(3),image_info['rectify'], image_info['left'][i],image_info['right'][i])

if __name__ == "__main__":
	main()

