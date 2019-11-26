# The python file for Exercise 4 - Neural Networks and 3D Shape Estimation
#
# In this exercise you will learn to build, train and test various simple neural networks. You will also apply a pre-trained model performing human body silhouette segmentation
# and use it silhouette to obtain a 3D human body model which you will render on a set of augmented reality images.
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
import cv2
import glob
import os
from utilities import UtilityFunctions # Importing class of utility functions from utilities.py file provided with the exercise, make sure to copy it in the same location together
				       # with cv_ex04_nn_and_shape.py
from human_model import HumanModel # Importing class containing functionality of fitting 3D statistical human body model parameters into a silhoutte. Make sure to copy this file in the same
				   # location as cv_ex04_nn_and_shape.py.
from cv_ex01_calibration import calibrateCamera, findChessBoardCorners, draw # Importing functions from cv_ex01_calibration.py file provided with the exercise, make sure to copy a working
									     # version of it in the same location together with cv_ex04_nn_and_shape.py
from keras.models import Sequential,load_model
from keras.layers import Dense, Dropout, Activation, Flatten,Input, Convolution2D, MaxPooling2D
from keras.engine.topology import Layer
from keras.utils import np_utils
from keras.datasets import mnist
from keras import backend as K
import tensorflow as tf
np.random.seed(123)  # for reproducibility

# Custom implementation of a dense/inner product layer in Keras
class MyInnerProductLayer(Layer):

    def __init__(self, output_dim, **kwargs):
        self.output_dim = output_dim
        super(MyInnerProductLayer, self).__init__(**kwargs)

    def build(self, input_shape):

        # Creating trainable weight variables for this layer.
        self.kernel = self.add_weight(name='kernel', shape=(input_shape[1], self.output_dim),initializer='uniform', trainable=True)
	self.bias = self.add_weight(shape=(self.output_dim,),initializer='uniform', name='bias',trainable=True)

        super(MyInnerProductLayer, self).build(input_shape)

    def call(self, x):
# ------------------------------------------------------------------------------------------- #
	# Task(1): Calculate a dot product between variable x and kernel
	# [--------LINE-TO-BE-INSERTED
	output = K.dot(x,self.kernel)
	# -------]
# ------------------------------------------------------------------------------------------- #
	output = K.bias_add(output, self.bias)
        return output

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_dim)

# Function creating a simple neural network model with two dense/inner product layers
def buildSimpleDigitRecognitionModel():

	model = Sequential()
	model.add(Flatten(input_shape=(28,28,1)))
# ------------------------------------------------------------------------------------------- #
	# Task(2): Add a 512 unit dense layer with reLU activation
        # [--------LINE-TO-BE-INSERTED
	model.add(Dense(512, activation='relu'))
	# ---------]
# ------------------------------------------------------------------------------------------- #
	model.add(Dropout(0.2))

	model.add(Dense(10, activation='softmax'))

	return model

# Function creating a simple neural network model with two custom dense/inner product layers
def buildCustomDigitRecognitionModel():
        model = Sequential()
        model.add(Flatten(input_shape=(28,28,1)))
# ------------------------------------------------------------------------------------------- #
	# Task(3): Add a 512 unit custom inner product layer MyInnerProductLayer
        # [--------LINE-TO-BE-INSERTED
	model.add(Dense(512, activation='tanh'))
	# -------]
# ------------------------------------------------------------------------------------------- #
	model.add(Activation('relu'))
        model.add(Dropout(0.2))

        model.add(Dense(10, activation='softmax'))
        return model

# Function creating a simple convolutional neural network model with two convolutional layers
def buildCNNBasedDigitRecognitionModel():
	model = Sequential ()

	model.add(Convolution2D(32, 3, 3, activation='relu',input_shape=(28,28,1)))
# ------------------------------------------------------------------------------------------- #
	# Task(4): Apply a 3x3 convolution with 32 output filters and relu activation function
        # [--------LINE-TO-BE-INSERTED
	model.add(Convolution2D(32, 3, 3, activation='relu'))
	# ---------]
# ------------------------------------------------------------------------------------------- #
	model.add(MaxPooling2D(pool_size=(2,2)))
	model.add(Dropout(0.25))
	model.add(Flatten())
	model.add(Dense(128, activation='relu'))
	model.add(Dropout(0.5))
	model.add(Dense(10, activation='softmax'))

	return model

# Function loading and preparing MNIST data for training and testing
def loadAndPrepareMNISTData():

        # Loading pre-shuffled MNIST data
        (train_images, train_labels), (test_images, test_labels) = mnist.load_data()

        # Reshaping images represented as 2 dimensional (N, w*h) into 4 dimensional arrays (N, w, h, 1) and normalising
	# gray scale intensity values to reside in an interval [0,1]
        train_images = train_images.reshape (train_images.shape[0], 28, 28,1).astype('float32')/255.0
        test_images = test_images.reshape (test_images.shape[0],  28, 28,1).astype('float32')/255.0

        # Converting digit class id's into a one-hot encoding
        train_labels = np_utils.to_categorical(train_labels, 10)
        test_labels = np_utils.to_categorical(test_labels, 10)

	return (train_images, train_labels), (test_images, test_labels)

# Function performing training of a Keras neural network model
def trainModels(train_images,train_labels,models):

	for model_name in models:

		model=models[model_name]

		print('Training model: ' +str(model_name))

		# Instructing Keras to print model summary
		model.summary()

		# Compiling a model selecting a categorical cross-entropy loss, 'adam' optimiser and accuracy metric
        	model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

		# Fitting model on training data
        	model.fit(train_images, train_labels, batch_size=32, nb_epoch=1, verbose=1)

# Function evaluating model acccuracy on provided MNIST data
def evaluateModels(test_images, test_labels, models):

	for model_name in models:

		model = models[model_name]

		print('-------------------------\nEvaluating model: '+ str(model_name))

		# Evaluating model on test data
	        score = model.evaluate(test_images, test_labels, verbose=0)

        	print('Test loss: ' + str(score[0]))
	        print('Test accuracy: ' + str(score[1]))

# Function building three different Keras models for digit recognition
def defineDigitRecognitionModels():
        models={}

        # Building a simple two layer neural network architecture based on fully connected layers
        models['simple_fully_connected'] = buildSimpleDigitRecognitionModel()

        # Building a simple two layer neural network architecture based on custom implementation of fully connected layers
        models['custom_fully_connected'] = buildCustomDigitRecognitionModel()

        # Building a neural network making use of convolutional layers
        models['convolution_based'] = buildCNNBasedDigitRecognitionModel()

	return models

#Function defining, training and evaluating various digit recognition models
def performDigitRecognitionExperiments():

	#Loading MNIST data using Keras functions
	(train_images, train_labels), (test_images, test_labels) = loadAndPrepareMNISTData()

	#Creating model definitions of various digit recognition models
	models = defineDigitRecognitionModels()

	#Training models on MNIST data
	trainModels(train_images, train_labels, models)

	#Evaluating models
	evaluateModels(test_images, test_labels, models)


# Function displaying a virtual cube and 3D human body model a list of images with chessboard pattern
def displayAugmentedObjects(image_list,result_dir,corner_size,camera_matrix,distortion_parameters, human_model):
	print('Calculating augmented reality patterns ...')

	# Initializing coordinates of a virtual cube
	cube_coordinates = np.float32([[0,0,0], [0,3,0], [3,3,0], [3,0,0], [0,0,-3],[0,3,-3],[3,3,-3],[3,0,-3] ])

        count=0

	for fname in image_list:

		print('Processing image: '+fname)
		image = cv2.imread(fname)
		gray_image = cv2.cvtColor(image,cv2.COLOR_BGR2GRAY)

		# Finding chessboard corners in an image
		ret, object_points, image_points = findChessBoardCorners(gray_image,corner_size)

    		if ret == True:

			# Converting chessboard corner points into a format required by solvePnPRansac
			object_points_Nx1x3 = np.zeros((object_points.shape[0], 1, 3), np.float32)
			object_points_Nx1x3[:, 0, :] = object_points[:, :]

			# Finding camera rotation and translation parameters for an augmented reality pattern
			_,rvecs, tvecs, inliers = cv2.solvePnPRansac(object_points_Nx1x3, image_points, camera_matrix, distortion_parameters)

			# Projecting virtual cube points using estimated camera parameters
        		projected_image_points, jac = cv2.projectPoints(cube_coordinates, rvecs, tvecs, camera_matrix,distortion_parameters)

			# Drawing projected virtual cube points into image coordinates
			image = draw(image,projected_image_points)

			# Obtaining 3D mesh coordinates of human body model. Note that counting variable is used as a state variable in order to create a "live" 3D human body model
			human_coordinates=human_model.generate3DMeshCoordinates(count)

			# Projecting obtained 3D human body model mesh coordinates into image coordinates
			projected_human_points, jac = cv2.projectPoints(human_coordinates, rvecs, tvecs, camera_matrix,distortion_parameters)

			# Drawing projected 3D human body model mesh into an image
			image=human_model.drawHumanModelOnImage(image, projected_human_points)

			count=count+1

			cv2.imwrite(result_dir+fname.split('/')[-1].split('.')[0]+'.png', image)

	print('Calculating virtual reality patterns - finished')

# Function performing estimated human body model visualisation
def visualiseHumanBodyModel():
	print('Visualising 3D human body model on a set of augmented reality images ...')

	# Loading SMPL human body model from file
	human_model = HumanModel.loadFromFile('m',1024,1024,'./images-ex04/results/mymodel_param.pic')

	calibration_image_list = UtilityFunctions.getImages('images-ex01/calibration/*.png')

	# Calibrating camera. Note that this function is borrowed from Exercise 1 code (cv_ex01_calibration.py) and assumes calibration images to be placed in images-ex01/calibration/
        reprojection_error, camera_matrix, distortion_parameters, rotation_vectors, translation_vectors = calibrateCamera(calibration_image_list,(9,7),(1920, 1080))

        vrpattern_image_list = UtilityFunctions.getImageFileList('images-ex04/ar-images-02/*.png')

	# Visualising virtual cube and 3D human body model avatar on a set of images
        displayAugmentedObjects(vrpattern_image_list,'images-ex04/results/ar-images/',(9,7),camera_matrix, distortion_parameters,human_model)

# Function loading a tensor consisting of images from provided list of directories
def loadImageTensor(image_dir_list):
	image_list=[]
        for image_dir in image_dir_list:
                for file_name in UtilityFunctions.getImageFileList(image_dir+'*.png'):
			assert (os.path.isfile(file_name)), "File: "+file_name + " does not exist!"
                        image = cv2.cvtColor(cv2.imread(file_name), cv2.COLOR_BGR2RGB)
                        image_list.append(image/255.0)
	return np.array(image_list)

def performHumanSilhouetteSegmentation(train_data,model,img_wh,image_dir_list):
	# Loading images into a tensor
        image_tensor=loadImageTensor(image_dir_list)
# ------------------------------------------------------------------------------------------- #
	# Task(5): Obtain segmentation network predition on image_tensor
        # [------LINE-TO-BE-INSERTED
 	output = model.predict(image_tensor, batch_size=2)
	# ----------]
# ------------------------------------------------------------------------------------------- #
	output=np.reshape(output,(image_tensor.shape[0],img_wh,img_wh,2))

        results_dir='images-ex04/results/silhouettes/'

	# Storing predicted silhouettes as png images in a results directory along with query images (for convenience)
        for id in range(0,image_tensor.shape[0]):
                cv2.imwrite(results_dir+'img_'+str(id).zfill(5)+'_lab.png',(255*output[id,:,:,1]))
                cv2.imwrite(results_dir+'img_'+str(id).zfill(5)+'.png',cv2.cvtColor((255*image_tensor[id,:,:,:]).astype(np.uint8),cv2.COLOR_BGR2RGB))


# Function loading a pre-trained silhouette segmentation network and applying to so segment on a set of images from a provided directory
def predictSilhouettes():
	print('Performing human silhouette segmentation from images ...')

	# Loading a pre-trained Keras model for human body model segmenantation
	silhouette_segmentation_model=load_model('images-ex04/segmentation-network/segmentation_model_weights.hdf5')

	# Segmenting human body silhoutttes in images from a given directory
	performHumanSilhouetteSegmentation(None,silhouette_segmentation_model,256,['images-ex04/person-01/'])

# Function extracting 3D human body model parameters from a human body silhoutte image and storing them on disk
def fit3DShape():
	print('Performing 3D model fitting into a given silhouette image ...')
	print('Note that fitting 3D model to a silhouette may take a few minutes. Once calculated result is stored on disk and loaded from for the purposed of further computation.')

	human_model_parameter_file_path='./images-ex04/results/mymodel_param.pic'
	silhouette_image_file_path='./images-ex04/results/silhouettes/img_00000_lab.png'

	if not (os.path.isfile(human_model_parameter_file_path)):
		assert (os.path.isfile(silhouette_image_file_path)), "Could not find silhouette image: "+silhouette_image_file_path
		human_model=HumanModel.fitHumanModelToSilhouette(silhouette_image_file_path)
		HumanModel.saveToFile(human_model,human_model_parameter_file_path)

# Function performing silhouette extraction from images, 3D body shape and pose fitting and its visualisation in an augmented reality demonstration
def performShapeEstimationExperiment():
	predictSilhouettes()
        fit3DShape()
        visualiseHumanBodyModel()

def main():

	# Commands used to force Keras and Tensorflow to use CPU as GPU resources are very limited. Do not run your model training  and testing on GPUs on MSALT machines
	os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
	os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
	tf.logging.set_verbosity(tf.logging.ERROR)

        with tf.device('/cpu:0'):
		performDigitRecognitionExperiments()
		performShapeEstimationExperiment()

	print("Finished: cv-ex04-nn-and-shape.py")

if __name__ == "__main__":
	main()

