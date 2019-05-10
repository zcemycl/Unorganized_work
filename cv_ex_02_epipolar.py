#Title: The python file for Exercise 2 - Epipolar Geometry
#
#Goal: The goal of this exercise is to build a framework for estimating epipolar lines on a pair of images.
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

import cv2
import cv2.xfeatures2d
import numpy as np
from matplotlib import pyplot as plt
import random
from utilities import UtilityFunctions #importing functions from utilities.py file provided with the exercise, make sure to copy it in the same location together with cv_ex02_epipolar.py

#Function drawing detected keypoints and epipolar lines on an image
def drawLines(image, lines, points, colors):
    r,c = image.shape
    image = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)

    for r, point, color in zip(lines, points, colors):

        x0, y0 = map(int, [0, -r[2] / r[1]])
        x1, y1 = map(int, [c, -(r[2] + r[0] * c) / r[1]])

        image = cv2.line(image, (x0, y0), (x1, y1), color, 1)
        image = cv2.circle(image, tuple(point), 5, color, -1)

    return image

#Function performing image keypoint and corresponding descriptor extraction for a list of images
def extractKeyPointsAndDescriptors(images):
	print('Extracting keypoints and descriptors...')
	sift = cv2.xfeatures2d.SIFT_create()
	keypoints=[]
	descriptors= []
	for image in images:

# ------------------------------------------------------------------------------------------------------ #
		# Obtain image keypoints and corresponding descriptors from an image using SIFT's detectAndCompute function
		# [--------LINE-TO-BE-INSERTED
		image_keypoints, image_descriptors = sift.detectAndCompute(image, None)

		# --------]

# ------------------------------------------------------------------------------------------------------ #
		keypoints.append(image_keypoints)
		descriptors.append(image_descriptors)

	return keypoints, descriptors

# Function finding epipolar lines corresponding given image point list and fundamental matrix
def findEpiLines(image_points, F):
	result=[]
	for i in range(2):
# ------------------------------------------------------------------------------------------------------ #
		# Compute epipolar line for image points using OpenCV function computeCorrespondEpilines
		# [--------LINE-TO-BE-INSERTED
                lines = cv2.computeCorrespondEpilines(image_points[i].reshape(-1,1,2), 2-i,  F)
		# --------]

# ------------------------------------------------------------------------------------------------------ #
		result.append(lines)

	return result

#Function finding matches of keypoint descriptors using Fann based matcher
def obtainDescriptorMatches(descriptors):

	# FLANN parameters
	FLANN_INDEX_KDTREE = 0
	index_params = dict(algorithm = FLANN_INDEX_KDTREE, trees = 5)
	search_params = dict(checks = 50)

	flann = cv2.FlannBasedMatcher(index_params, search_params)

# ------------------------------------------------------------------------------------------------------ #
	# Obtain matches between descriptors using Flann based nearest neighbour matcher
	# [--------LINE-TO-BE-INSERTED
        matches = flann.knnMatch(descriptors[0],descriptors[1], k = len(descriptors))
	# --------]

# ------------------------------------------------------------------------------------------------------ #
	return matches

#Function filtering strongly matching points
def filterStronglyMatchingPoints(keypoints, matches):
	points = [[],[]]

	# Ratio test
	for i,(m,n) in enumerate(matches):
		if m.distance < 0.6 * n.distance:
			points[1].append(keypoints[1][m.trainIdx].pt)
			points[0].append(keypoints[0][m.queryIdx].pt)

	points[1]=np.int32(points[1])
	points[0]=np.int32(points[0])
	return points

#Function finding strong keypoint matches and estimating fundamental matrix
def estimateFundamentalMatrixAndMatchingPoints(keypoints, descriptors):

	print('Estimating fundamental matrix...')

	matches= obtainDescriptorMatches(descriptors)

	points=filterStronglyMatchingPoints(keypoints, matches)

# ------------------------------------------------------------------------------------------------------ #
	# Find fundamental matrix and valid point mask using OpenCV function findFundamentalMat
	# [--------LINE-TO-BE-INSERTED
        F, mask = cv2.findFundamentalMat(points[0], points[1], cv2.FM_LMEDS)
	# --------]

# ------------------------------------------------------------------------------------------------------ #
	masked_points= [points[0][mask.ravel() == 1], points[1][mask.ravel() == 1]]

	return F, masked_points

#Function performing epipolar line visualisation on a list of images
def drawEpiLines(images,points, epi_lines, pair_id):

	print('Finding epipolar lines...')

	colors = [tuple(np.random.randint(0,255,3).tolist()) for i in range(len(points[0]))]

	result= [drawLines(images[i],epi_lines[i].reshape(-1,3),points[i],colors) for i in range(2)]

	cv2.imwrite('images-ex02/results/'+pair_id+'.png',np.concatenate(tuple(result),axis=1))

def main():

	for pair_id in ['pair-01','pair-02']:
		print('------------\nAnalysing image pair: '+pair_id)

		images = UtilityFunctions.getImages('images-ex02/'+pair_id+'/*')

		keypoints, descriptors = extractKeyPointsAndDescriptors(images)

		F, points = estimateFundamentalMatrixAndMatchingPoints(keypoints,descriptors)

		epi_lines= findEpiLines(points,F)

		drawEpiLines(images,points,epi_lines,pair_id)


	print('Finished: cv-ex02-epipolar.py')

if __name__ == "__main__":
   main()
