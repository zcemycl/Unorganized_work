# This file performs various functions to do with SMPL statistical human body model loading, saving and fitting it to an image silhouette. Note that students are not required to understand
# the contents of this file. If interested to learn more, visit the tutorials and libraries referred in this file as well as the main exercise.

import os
import numpy as np
from opendr.renderer import ColoredRenderer  # OpenDR is a differenciable renderer. See GitHub on:  https://github.com/mattloper/opendr/wiki
from opendr.lighting import LambertianPointLight, SphericalHarmonics
from opendr.camera import ProjectPoints
from opendr.filters import gaussian_pyramid
import chumpy as ch
import cPickle as pickle
from opendr.geometry import VertNormals
import sys
sys.path.append('/usr/local/teach/smpl/')  # Appending a system path variable with folder to SMPL python libraries
from smpl_webuser.serialization import load_model  # You can find SMPL library code here: http://smpl.is.tue.mpg.de/downloads
import cv2

# This class encapsulates function performing saving, loading and fitting of SMPL model
class HumanModel:
	# Initialising class elements
	def __init__(self,gender,render_w,render_h):

		# Both male and female body models are available
		assert (gender in['m','f']), "Invalid gender token"
		if(gender == 'm'):
			self.m = load_model('./images-ex04/smpl-models/basicModel_f_lbs_10_207_0_v1.0.0.pkl')
		else:
			self.m = load_model('./images-ex04/smpl-models/basicModel_m_lbs_10_207_0_v1.0.0.pkl')

		self.rn=None
		self.rn_color=None
		self.initializeParameters()
		self.loadParametersToModel(self.parameters)
		self.loadCameraToRender(render_w,render_h)
	        self.N_period=21
	        self.steps=[i for i in range(self.N_period)]
	        for i in reversed(range(1,self.N_period-1)):
        	        self.steps.append(i)

	# Setting up camera parameters for SMPL model to be rendereed by OpenDR
	def loadCameraToRender(self,w,h):
		# Creating OpenDR renderer for a silhouette of a 3D model
		self.rn = ColoredRenderer()
		self.rn.camera = ProjectPoints(v=self.m, rt=np.zeros(3), t=np.array([0, 0, 2.]), f=np.array([w,w])/2., c=np.array([w,h])/2., k=np.zeros(5))
		self.rn.frustum = {'near': 1., 'far': 10., 'width': w, 'height': h}
		self.rn.set(v=self.rn.v, f=self.m.f, vc=np.ones_like(self.m), bgcolor=np.array([0.0,0.0,0.0]))

		# Creating OpenDR renderer for a full 3D mesh of a 3D model
                self.rn_color = ColoredRenderer()
                self.rn_color.camera = ProjectPoints(v=self.m, rt=np.zeros(3), t=np.array([0, 0, 2.]), f=np.array([w,w])/2., c=np.array([w,h])/2., k=np.zeros(5))
                self.rn_color.frustum = {'near': 1., 'far': 10., 'width': w, 'height': h}
                self.rn_color.set(v=self.rn.v, f=self.m.f, vc=np.ones_like(self.m), bgcolor=np.array([0.0,0.0,0.0]))
                self.rn_color.vc = LambertianPointLight(f=self.m.f, v=self.rn_color.v, num_verts=len(self.m), light_pos=np.array([0,0,-20]),  vc=np.ones_like(self.m)*.9,   light_color=np.array([1., 1., 1.]))

	# Function saving 82 human body model parameters into a file : pose (72)  and shape (10)
	@staticmethod
	def saveToFile(human_model,file_name):
		pickle.dump(human_model.parameters,open(file_name,'wb'))

	# Function loading 82 human body model parameters from files
	@staticmethod
	def loadFromFile(gender,w,h,file_name):
		assert(os.path.isfile(file_name)),"File "+file_name + " does not exist."
		hm = HumanModel(gender,w,h)
		hm.parameters = pickle.load(open(file_name,'rb'))
		hm.loadParametersToModel(hm.parameters)
		return hm

	# Function transferring internal human body model parameters of the parent class into SMPL human body model object
	def loadParametersToModel(self,parameters):
		self.m.pose[:]=parameters['pose'][:].copy()
		self.m.betas[:]=parameters['betas'][:].copy()

	# Function instantiating standard pose and shape parameters
	def initializeParameters(self):
		self.parameters={'pose': np.zeros(self.m.pose.size),
				 'betas':np.zeros(self.m.betas.size)}
		self.parameters['pose'][0] = np.pi

	# Function obtaining a perturbed set of parameters, sampled from chosen intervals uniformly. Note that only hip (5,8) and shoulder (41, 44) parameters are
	# changed, along with several shape (betas) parameters.
	def obtainRandomParameters(self,pose_ids,betas_ids):
		assert( len((set([5,8,41,44]))&set(pose_ids)) == 4), "Pose parameter list is too limited"
		result={'pose': self.parameters['pose'],'betas': self.parameters['betas']}
		result['pose'][5]= np.pi / 32 + (1-2*np.random.rand())*np.pi/(32*2)
		result['pose'][8]=- (np.pi / 32 + (1-2*np.random.rand())*np.pi/(32*2))
		result['pose'][41]=-(10*np.pi/32) +(1-2*np.random.rand())*(np.pi/16)
                result['pose'][44]=-(-(10*np.pi/32) +(1-2*np.random.rand())*(np.pi/16))
		result['betas'][betas_ids]=(np.ones(self.m.betas.size)[betas_ids]-2*np.random.rand(self.m.betas.size)[betas_ids])*10.0

		return result

	# Function updating parameters stored by HumanModel class
	def updateParameters(self,parameters):
		self.parameters=parameters.copy()

	# Function extracting boolean silhouette image from RGB image
	@staticmethod
	def extractSilhouetteFromRGB(image):
		result=(image[:,:,0]+image[:,:,1]+image[:,:,2]>0)
		return result

	# Function calculating total number of pixels with correctly predicted foreground (person) and background values
	@staticmethod
	def calculateOverlap(prediction_sil,gt_sil):
		return np.count_nonzero((prediction_sil==gt_sil))

	# Function creating RGB image from a boolean image
	@staticmethod
	def rgbFromBoolean(bool_img):
		result=np.zeros((bool_img.shape[0],bool_img.shape[1],3))
		result[bool_img==True,:]=(255,255,255)
		return result

	# Function generating 3D coordinates of a human body model mesh. Count variable is used to encode the state of the mesh for the purpose of the visualisation of the waving hand
	def generate3DMeshCoordinates(self,count):
                self.m.pose[0]= -np.pi/2
                self.m.pose[41] = 0
                self.m.pose[56] =(1-2*((self.steps[count%(2*self.N_period-2)])/(self.N_period*1.0)))*np.pi/2+np.pi/2
                self.m.pose[60] =-np.pi/2
                coordinates=np.array(self.rn.v*10)
                return coordinates

	# Function, drawing human body visualisation onto an image. Note that projected 2D coordinates of 3D points are used. Visualisation consists of edges of the triangles of the mesh points.
	def drawHumanModelOnImage(self, image, image_points):
		faces=np.array(self.rn.f)
		image_points = np.int32(image_points).reshape(-1,2)
		for id in range(faces.shape[0]):
			fid=np.int32(faces[id])
		        pt0=(int(image_points[fid[0]][0]),int(image_points[fid[0]][1]))
		        pt1=(int(image_points[fid[1]][0]),int(image_points[fid[1]][1]))
		        pt2=(int(image_points[fid[2]][0]),int(image_points[fid[2]][1]))
		        points = np.array([[pt0[0],pt0[1]],[pt1[0],pt1[1]],[pt2[0],pt2[1]]])
		        image=cv2.polylines(image,np.int32([points]),True,(100,250,250))

		return image
	# Function cropping, centering and resising a silhouette from an image
	@staticmethod
	def crop_silhouette(silhouette_image):
		silhouette_coordinates = np.where(silhouette_image != 0)
		min_0=np.min(silhouette_coordinates[0])
		max_0=np.max(silhouette_coordinates[0])
		min_1=np.min(silhouette_coordinates[1])
		max_1=np.max(silhouette_coordinates[1])
		len_0=max_0-min_0+1
		len_1=max_1-min_1+1
		if(len_0 % 2 ==0):
			if(min_0>0):
				min_0-=1
				len_0+=1
			else:
				min_0+=1
				len_0-=1

		if(len_1 % 2 ==0):
	        	if(min_1>0):
	                        min_1-=1
	                        len_1+=1
	                else:
	                        min_1+=1
	                        len_1-=1


		max_len=max(len_1,len_0)

		bounding_box = silhouette_image[min_0:(max_0+1), min_1:(max_1+1),:]
		squared_bounding_box=np.zeros((max_len,max_len,3))
		squared_bounding_box[(0+(max_len-len_0)/2):(max_len-(max_len-len_0)/2),(0+(max_len-len_1)/2):(max_len-(max_len-len_1)/2),:]=bounding_box
		resized_sqb_box=cv2.resize(squared_bounding_box,(224,224),interpolation=cv2.INTER_NEAREST)
		result=np.zeros((256,256,3))

		result[((256-224)/2):(256-((256-224)/2)),((256-224)/2):(256-((256-224)/2)),:] = resized_sqb_box

		return result

	# Function implementing a single step of a very simple optimisation procedure. It takes a set of shape and pose ids and uniformly samples values from selected intervals. If the silhoutte of the 
	# resulting body model has a higher overlap count, than the best found so far, the former value is updated.
	def updateStep(self,gt_silhouette,total_iter,iter,best_overlap, best_param):
		if(iter% 100 ==0):
        	        print('Iteration: ' +str(iter) + ' of '+str(total_iter))
                pose_ids = [5,8,41,44]
                shape_ids = [0,1]
                param = self.obtainRandomParameters(pose_ids,shape_ids)
                self.loadParametersToModel(param)
                predicted_cropped_silhouette = HumanModel.crop_silhouette(self.rn.r)
                boolean_pc_silhouette = HumanModel.extractSilhouetteFromRGB(predicted_cropped_silhouette)
                overlap = HumanModel.calculateOverlap(boolean_pc_silhouette,gt_silhouette)
                if(overlap > best_overlap):
                        best_overlap = overlap
                        self.updateParameters(param)
                        result_image=np.concatenate((cv2.resize((self.rn_color.r*255),(256,256),interpolation=cv2.INTER_NEAREST),
						HumanModel.rgbFromBoolean(boolean_pc_silhouette),
						HumanModel.rgbFromBoolean(gt_silhouette),
						HumanModel.rgbFromBoolean(gt_silhouette==boolean_pc_silhouette)), axis=1)
                        cv2.imwrite('images-ex04/results/res_img_'+str(iter).zfill(6)+'.png',result_image)
			best_param={'pose':param['pose'][:].copy(),'betas':param['betas'][:].copy()}

		return best_overlap, best_param

	# Function performing a simple optimisation routine to fit a 3D model into a silhouette image
	@staticmethod
	def fitHumanModelToSilhouette(image_name):
		gt_silhouette = cv2.imread(image_name)
		gt_silhouette [gt_silhouette[:,:,0] < 127, :] = [0, 0, 0]

		cropped_gt_silhouette = HumanModel.crop_silhouette(gt_silhouette)
		gt_silhouette = HumanModel.extractSilhouetteFromRGB(cropped_gt_silhouette)

		human_model = HumanModel('m', 1024, 1024)
		os.system("rm images-ex04/results/res_img*.png")

		best_overlap=0
		best_parameters=None
		total_iter = 1000

		for iter in range(0,total_iter):
			best_overlap, best_parameters = human_model.updateStep(gt_silhouette, total_iter, iter, best_overlap, best_parameters)

		human_model.loadParametersToModel(best_parameters)
		human_model.updateParameters(best_parameters)
		return human_model

