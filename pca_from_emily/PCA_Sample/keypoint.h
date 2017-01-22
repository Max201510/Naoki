/* -*-  Mode:C++; c-basic-offset:8; tab-width:8; indent-tabs-mode:t -*- */
/*
Author: Yan Ke
Dec 2003
*/

#ifndef KEYPOINT_H
#define KEYPOINT_H

#include <vector>
using namespace std;

#include "image.h"

/**
   Number of components in PCA basis vector (in input file).  This is
   always 36, or else needs retraining.
*/
#define PCALEN 36

#define EPCALEN 20 //< Effective number of pca components (calculation and file output).

/**
   Window size of image patch for keypoint.
   PatchSize must be odd
   Must be 41 because the PCA basis vectors were trained on
   this size.  Changing this requires retraining.  If the window
   size of image is larger, than that window is resized to this size.
*/
#define PatchSize 41

#define PatchLength (PatchSize * PatchSize)

/**
   Length one element gradient basis vector
*/
#define GPLEN ((PatchSize - 2) * (PatchSize - 2) * 2)


/**
 * Normalizes a vector to unit length.
 * Used for normalizing gradients of an image patch.
 *
 * @param v Vector to normalize.
 */
void Normvec(float * v, unsigned int len);

/**
 * Class for one keypoint.
 * This class contains a keypoint's location, orientation, scale,
 * and PCA local descriptor.
 */

class Keypoint {
public:


	float x; ///< x location of keypoint in base image
	float y; ///< y location of keypoint in base image
	
	float sx; ///< x location of keypoint in scaled image
	float sy; ///< y location of keypoint in scaled image
	
	int octave; ///< octave index of keypoint
	int scale; ///< scale index of keypoint
	float fscale; ///< interpolated scale index

	float ori; ///< orientation of keypoint.  Ranges from (-PI, PI)

	float gscale; ///< scale of keypoint in pixels in base image.

        vector<float> ld; ///< PCA local descriptor

	Image * patch; ///< image patch
                                                                                                                       
        Keypoint(int ldlen) {
		ld.resize(ldlen);
                patch = NULL;
        }
                                                                                                                       
        ~Keypoint() {
                if (patch != NULL)
                        delete (patch);
                patch = NULL;
        }


};

/**
 * Keypoint detector class.
 * Finds keypoints using David Lowe's SIFT algorithm.
 * Computes local descriptors using PCASIFT.
 *
 * Limiations: Does not fully implement Lowe's SIFT algorithm.
 * Does not interpolate in location and scale.
 */
class KeypointDetector {
private:

	float avgs[GPLEN]; ///< average local descriptor value.  Used in PCA.
	float eigs[EPCALEN][GPLEN]; ///< PCA basis vectors

	/**
	 * Scales and blurs input image to make base image.
	 * Possibly doubles image size, based on DOUBLE_BASE_IMAGE_SIZE.
	 * Blur image based on difference of SIGMA and INITSIGMA.
	 *
	 * @param image Raw image to find keypoints on.
	 * @return Base image for the bottom of the Gaussian pyramid.
	 */
	static Image * ScaleInitImage(Image * image);

	/**
	 * Given a base image, computes a gaussian pyramid.
	 * Base image must be blurred and scaled appropriately.
	 *
	 * @param image Base image for the bottom of the Gaussian pyramid.
	 * @return Vector of vector of Image *.  Outside vector are the octaves.
	 *   Inside vector are the scales of images.
	 * @see BuildGaussianScales()
	 * @see BuildDOGOctaves()
	 * @see ScaleInitImage()
	 */
	static vector<vector<Image *> > BuildGaussianOctaves(Image * image);

	/**
	 * Given a base image, computes a gaussian pyramid for a particular octave.
	 * Number of scales in an octave based on SCALESPEROCTAVE.
	 *
	 * @param image Base image for the bottom of this octave.
	 * @return Vector of Image *, the scales of this octave.
	 * @see BuildGaussianOctaves()
	 */
	static vector<Image *> BuildGaussianScales(Image * image);

	/**
	 * Given a Gaussian pyramid, builds the difference of Gaussians (DoG) pyramid.
	 * Simply takes the pixel difference between adjacent scales in an octave.
	 * 
	 * @param GOctaves Gaussian pyramid.
	 * @return DoG pyramid.  Outside vector are the octaves.  Inside vector are
	 *   the scales.
	 *   Number of octaves returned is the same as the input
	 *   Gaussian pyramid.  Number of scales per octave returned will be 1 less
	 *   than the input Gaussian pyramid.
	 * @see BuildDOGScales()
	 * @see BuildGaussianOctaves()
	 */
	vector<vector<Image *> > BuildDOGOctaves(vector<vector<Image *> > & GOctaves);

	/**
	 * Given a Gaussian octave, builds the difference of Gaussians (DoG) octave.
	 * Simply takes the pixel difference between adjacent scales.
	 * 
	 * @param LImages One Gaussian octave.
	 * @return OneDoG octave.\ Number of scales returned will be 1 less
	 *   than the number of scales in the input Gaussian octave.
	 * @see BuildDOGOctave()
	 */
	vector<Image *> BuildDOGScales(vector<Image *> & LImages);

	/**
	 * Given a set of DoG octaves, returns a list of Keypoints.
	 * Localizes keypoints in space and scale.
	 * Only finds x, y, sx, sy, ocative, scale,  and gscale of a keypoint.
	 * Does NOT compute orientation or local descriptor of keypoint.
	 *
	 * @param DOctaves DoG pyramid.
	 * @return Vector of Keypoint *.
	 * @see FindPeaksScales()
	 * @see BuildDOGOctaves()
	 * @see ComputeLocalDecr()
	 * @see FindOrientation()
	 */
	vector<Keypoint *> FindPeaksOctaves(vector<vector<Image *> > & DOctaves);

	/**
	 * Given a set of DoG scales, returns a list of Keypoints.
	 * Localizes keypoints in space and scale.
	 * Only finds x, y, sx, sy, ocative, scale, and gscale of a keypoint.
	 * Does NOT compute orientation or local descriptor  of keypoint.
	 * 
	 * @param octave octave number of this set of scales.
	 * @param DImages One DoG octave.
	 * @return Vector of Keypoint * found in this octave.
	 * @see FindPeaksOctaves()
	 */
	vector<Keypoint *> FindPeaksScales(int octave, vector<Image *> & DImages);
	/**
	 * Interpolates keypoints to subpixel and subscale accuracy.
	 *
	 * @param x x pixel location
	 * @param y y pixel location
	 * @param s s Scale
	 * @param DImages DoG pyramid.
	 * @param fx interpolated x position
	 * @param fy interpolated y position
	 * @param fs interpolated scale
	 * @return true if key should be added
	 */
	bool InterpKey(int x, int y, int s, vector<Image *> & DImages,
		       float * fx, float * fy, float * fs);

	/**
	 * Interpolates keypoints to subpixel and subscale accuracy.
	 *
	 * @param x x pixel location
	 * @param y y pixel location
	 * @param s s Scale
	 * @param DImages DoG pyramid.
	 * @param dx delta x
	 * @param dy delta y
	 * @param ds delta s
	 * @return function value at interpolated position
	 */
	float InterpKeyStep(int x, int y, int s, vector<Image *> & DImages,
		       float * dx, float * dy, float * ds);

	/**
	 * Checks to see if the peak lies along an edge.
	 * Uses the ratio of the eigenvalues of the Hessian matrix to determine
	 * the ratio of the curvature of the two principal axis.
	 * If the ratio is above EDGE_PEAK_THRESH, then the peak is on an edge.
	 *
	 * @param x x pixel location of peak
	 * @param y y pixel location of peak
	 * @param dimage DoG image for a particular scale.
	 * @return true if peak lies along an edge
	 */	
	bool isEdgePeak(int x, int y, Image * dimage);

	/**
	 * Determins if a pixel is a local peak.
	 * Determines if a pixel is at a local maximum or a local minimum.
	 * Checks the 26 pixels (9 above, 9 below, 8 on the same plane) around
	 * the center pixel to see if it's a local peak.
	 *
	 * @param aimage DoG image (scale below, smaller scale).
	 * @param bimage DoG image to detect peak medium scale).
	 * @param cimage DoG image (scale above, larger scale).
	 * @param x x pixel location.
	 * @param y y pixel location.
	 * @return true if pixel at (x, y) is a local peak.
	 */
	bool isPeak(Image * aimage, Image * bimage,
		   Image * cimage, int x, int y);

	/** 
	 * Calculates the gradient centered at a pixel.
	 * Gradient is calculated using a 3x3 square around center.
	 * 
	 * @param x x pixel location
	 * @param y y pixel location
	 * @param image Image to find gradient.
	 * @param m returns magnitude of gradient in m.
	 * @param theta returns the orientation of gradient in theta (-PI, PI).
	 * @return true if pixel is a valid location and gradient is found.
	    Also writes *m and *theta.
	 */
	bool GetPixOrientation(int x, int y, Image * image, float * m,  float * theta);
	
	/**
	 * Finds the dominant orientation(s) of a keypoint.
	 * If more than one dominant orientation is found, then generate additional
	 * keypoints.  New keypoints are generated with threshold NEW_KP_OR_THRESH.
	 *
	 * @param kps Keypoints localized in scale and space, but without orientations.
	 * @param LOctaves Gaussian pyramid.
	 * @return vector of new Keypoint *.  This should be merged with kps.
	 * @see BuildGaussianOctaves().
	 * @see FindPeaksOctaves().
	 */
	vector<Keypoint *> FindOrientation(vector<Keypoint *> & kps,
					   vector<vector<Image *> > & LOctaves);
	
	/**
	 * Computes the PCA local descriptor of a keypoint
	 * Stores local descriptor in keypoint ld.
	 *
	 * @param key A keypoint localized in space, scale, and orientation.
	 * @param blur The Gaussian scale closest to the scale of the keypoint.
	 */
	void MakeKeypointPCA(Keypoint * key, Image * blur);

	/**
	 * Extracts the image gradient patch of a keypoint.
	 * Size of image patch determined by PatchSize.  Gradient patch is 2
	 * pixels smaller in each dimension.  We compute the gradient in
	 * both the x and y direction.
	 *
	 * @param key Keypoint localized in space, scale, and orientation.
	 * @param blur Gaussian image closest in scale to the keypoint.
	 * @param v Return.  Gradient patch flattened into one long vector.  We store
	 *   gx and gy interleaved.
	 */
	void MakeKeypointPatch(Keypoint * key, Image * blur, float v[GPLEN]);

	/**
	 * Computes the PCA local descriptor for keypoints.
	 * Stores local desciptor in the ld field of the Keypoint.
	 *
	 * @param kps Keypoints localized in space, scale, and orientation.
	 * @param LOctaves Gaussian pyramid.
	 * @see BuildGaussianOctaves()
	 */
	void ComputeLocalDescr(vector<Keypoint *> & kps,
			       vector<vector<Image *> > & LOctaves);


	static void RecalcKeyIndices(vector <Keypoint *> & keys);

	static void MakeLocalPatch(Keypoint * key, Image * blur, int windowsize);

	static void ComputeLocalPatches(vector<Keypoint *> & kps,
					vector<vector<Image *> > & LOctaves, int patchsize);

public:	
	/**
	 * Constructs a new KeypointDetector object.
	 *
	 * @param basisfile File name of the PCA basis vector file and averages.
	 */
	KeypointDetector(char * basisfile);

	/**
	 * Finds all of the keypoints in an Image.
	 *
	 * @param image Input image.
	 * @return Vector of Keypoint *
	 */
	vector <Keypoint *> FindKeypoints(Image * image);


	/**
	 * Recalculates local descriptor of keys on the given image.
	 * Given localized keys, calculate local descriptor.
	 *
	 * @param image Image to calculate.
	 * @param keys Vector of keys with position already calculated.
	 */
	void RecalcKeys(Image * image, vector <Keypoint *> & keys);

	/**
	 * Writes keypoints to a file
	 *
	 * @param keys Keypoints to be written.
	 * @param fn File name of output.
	 */
	static void writeKeysToFile(vector<Keypoint *> & keys, char * fn);

	/**
	 * Read keys from a file.
	 * 
	 * @param fn File name of input.
	 * @return Vector of keys.
	 */
	static vector<Keypoint *> readKeysFromFile(char * fn);

	/**
	 * Read patches from a file.
	 * 
	 * @param fn File name of input.
	 * @return Vector of keys.
	 */
	static vector<Keypoint *> readPatchesFromFile(char * fn);


       /** Get image patch.
         * @param image Image to calculate.
         * @param keys Vector of keys with position already calculated.
         * @param patchsize size of square patch in pixels.
         */
        static void getPatches(Image * im, vector <Keypoint *> & keys, int patchsize);
 
        /** Write patches to file.
         * @param keys Keypoints to be written.
         * @param fn File name of output.
         * @param patchsize size of square patch in pixels.
         */
        static void writePatchesToFile(vector<Keypoint *> & keys, char * fn, int patchsize);

};

#endif
