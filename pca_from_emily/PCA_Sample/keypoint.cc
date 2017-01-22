/* -*-  Mode:C++; c-basic-offset:8; tab-width:8; indent-tabs-mode:t -*- */
/*
Author: Yan Ke
Dec 2003
*/

#include <math.h>
#include <assert.h>

#include "keypoint.h"

#define PI 3.14159265358979323846


/**
   Keypoint window sample size is PatchMag * scale
   Smallest scale is typically around 2.0, so window size
   is at least 40.
*/
#define PatchMag 20


/**
   Sigma of base image -- See D.L.'s paper.
*/
#define INITSIGMA 0.5

/**
   Sigma of each octave -- See D.L.'s paper.
*/
#define SIGMA 1.6

/**
   Number of scales per octave.  See D.L.'s paper.
*/
#define SCALESPEROCTAVE 3
#define MAXOCTAVES 14

/**
   Number of bins in histogram bin
*/
#define NUMORIENTATIONS 36
#define DEGPERBIN (360 / NUMORIENTATIONS)

/**
   If peak value in histogram bin is at least this multiple
   of the maximum, then create new keypoint. See
   D.L.'s paper.
*/
#define NEW_KP_OR_THRESH 0.8

/**
   Intensity of DoG for consideration for keypoint.
   See D.L.'s paper.  Intensity value must be above this,
   divided by SCALESPEROCTAVE.
*/
#define LOW_CONTRAST_THRESH 0.04

/**
   Double image size before looking searching for keypoints
   Doubling finds more keypoints, but takes longer. See
   D.L.'s paper.
*/
#define DOUBLE_BASE_IMAGE_SIZE 1

/**
   Rejection threshold the edginess of a peak.  We can't
   accurately localize peaks that lie on edges.  Se D.L.'s paper.
*/
#define EDGE_PEAK_THRESH 10.0

/**
   Turn on interpolation in space and scale for keypoints.
   See D.L.'s paper.
*/
#define INTERP_KEYS 1

/**
   Don't look for keypoints that are within this distance
   from the edge of the image.
*/
#define BORDER 5

int scalestat[100] = {0};

void printstat() {
	printf("Scale stat:\n");
	for (int i = 0; i < 20; i++) {
		printf("%d\n", scalestat[i]);
	}
}

/**
 * Smoothes the values in a histogram
 * 
 * @param x Input histogram array
 * @param n size of histogram array
 */

void SmoothHistogram(float * x, int n) {
	assert(x);

	float temp = 0, result;
	int prev, next;

	for ( int i = 0; i < n; i ++){
		prev = (i-1) < 0? n-1: i-1;
		next =  (i+1)%n;
		result = (x[prev]+ x[i]  + x[next])/3.0;
		x[prev] = temp;
		temp = result;
	}
	
	x[n-1] = temp;
}

/**
 * Finds the maximum value in a array.
 *
 * @param v Input array
 * @param len Number of elements in array.
 * @return maximum value found
 */
float maxvec(float * v, int len) {
	assert(v);

	int maxbin = 0;
	float maxval = v[0];

        // find maximum
	for (int j = 1; j < NUMORIENTATIONS; j++) {
		if (v[j] > maxval) {
			maxval = v[j];
			maxbin = j;
		}
	}

        return maxval;
}

/**
 * Mod function that returns all positive values.
 * 
 * @param x value to take mod
 * @param b modulo
 * @return positive value of x%b.
 */
int mod(int x, int b) {
        x %= b;
        if (x < 0)
                x += b;

        return x;

}


/**
 * Calculates the inverse of a 3x3 matrix.
 *
 * @param a Input matrix
 * @param b Output inverse matrix.
 */
void mInv33(float a[3][3], float b[3][3]) {
	assert(a);
	assert(b);
	
	float det = a[0][0]*(a[2][2]*a[1][1] - a[2][1]*a[1][2])
		- a[1][0]*(a[2][2]*a[0][1] - a[2][1]*a[0][2])
		+ a[2][0]*(a[1][2]*a[0][1] - a[1][1] * a[0][2]);

	//printf("det: %f\n", det);

	b[0][0] = (a[2][2]*a[1][1] - a[2][1]*a[1][2]) / det;
	b[0][1] = (-a[2][2]*a[0][1] + a[2][1]*a[0][2]) / det;
	b[0][2] = (a[1][2]*a[0][1] - a[1][1]*a[0][2]) / det;

	b[1][0] = (-a[2][2]*a[1][0] + a[2][0]*a[1][2]) / det;
	b[1][1] = (a[2][2]*a[0][0] - a[2][0]*a[0][2]) / det;
	b[1][2] = (-a[1][2]*a[0][0] + a[1][0]*a[0][2]) / det;

	b[2][0] = (a[2][1]*a[1][0] - a[2][0]*a[1][1]) / det;
	b[2][1] = (-a[2][1]*a[0][0] + a[2][0]*a[0][1]) / det;
	b[2][2] = (a[1][1]*a[0][0] - a[1][0]*a[0][1]) / det;
}

/**
 * Finds the x coordinate of the maximum (or minimum) of the parabola described by 3 points.
 *
 * @param x1 x coordinate of point 1
 * @param y1 y coordinate of point 1
 * @param x2 x coordinate of point 2
 * @param y2 y coordinate of point 2
 * @param x3 x coordinate of point 3
 * @param y3 y coordinate of point 3
 * @return x coordinate of the maximum (or minimum) of the parabola.
 */
float getMaxParabola(float x1, float y1,
                     float x2, float y2,
                     float x3, float y3) {

	float mb[3] = {y1, y2, y3};

	float mA[3][3] = {{x1*x1, x1, 1},
			  {x2*x2, x2, 1},
			  {x3*x3, x3, 1}};
	float mAInv[3][3];
	
	mInv33(mA, mAInv);
	
	float a = 0;

	// mx = maInv * mb
	for (int i = 0; i < 3; i++)
		a += mAInv[0][i] * mb[i];

	float b = 0;
	for (int i = 0; i < 3; i++)
		b += mAInv[1][i] * mb[i];

        return -b / (2*a);
}


void Normvec(float * v, unsigned int len) {
	float total = 0;
	

	for (unsigned int i = 0; i < len; i++) {
		total += fabs(v[i]);
	}
	

	total /= len;

	
	for (unsigned int i = 0; i < len; i++) {
		v[i] = v[i] / total * 100.0;
	}	
}

KeypointDetector::KeypointDetector(char * fn) {
	assert(fn);

	FILE * pcaf = fopen(fn, "rb");

	printf("Reading averages\n");
 
        for (int i = 0; i < GPLEN; i++) {
                float val;
                if (fscanf(pcaf, "%f", &val) != 1) {
                       printf("Invalid pca vector file (avg).\n");
		       exit(1);

		}
                avgs[i] = val;
        }
 
        printf("Reading pca vector %dx%d\n", GPLEN, PCALEN);
 
        // read in vector, transposed
 
        for (int i = 0; i < GPLEN; i++) {
                for (int j = 0; j < PCALEN; j++) {
 
                        float val;
                        if (fscanf(pcaf, "%f", &val) != 1) {
                                printf("Invalid pca vector file (eig).\n");
				exit(1);
			}
                        
			if (j < EPCALEN)
				eigs[j][i] = val;
                }
        }
 
        fclose(pcaf);


}


bool KeypointDetector::isPeak(Image * aimage, Image * bimage,
			      Image * cimage, int x, int y) {

	assert(aimage);
	assert(bimage);
	assert(cimage);


	Image * ims[3];

	ims[0] = aimage;
	ims[1] = bimage;
	ims[2] = cimage;

	float center = bimage->getPixel(x, y);

	if (center > 0.0) {

		for (int si = 0; si < 3; si++) {
			for (int yi = y - 1; yi <= y + 1; yi++) {
				for (int xi = x - 1; xi <= x + 1; xi++) {
					float p2 = ims[si]->getPixel(xi, yi);
					
					if (center < p2)
						return false;
				}
			}
		}
	}

	else {
		for (int si = 0; si < 3; si++) {
			for (int yi = y - 1; yi <= y + 1; yi++) {
				for (int xi = x - 1; xi <= x + 1; xi++) {
					float p2 = ims[si]->getPixel(xi, yi);
					
					if (center > p2)
						return false;
				}
			}
		}
	}

	return true;
}
    

vector<Keypoint *> KeypointDetector::FindPeaksOctaves(vector<vector<Image *> > & DOctaves) {
	vector<Keypoint *> peaks;

	printf("FindPeaks(): Finding peaks...\n");                                       

	for (unsigned int i = 0; i < DOctaves.size(); i++) {
		vector<Keypoint *> ps;
	        ps = FindPeaksScales(i, DOctaves[i]);
		peaks.insert(peaks.begin(), ps.begin(), ps.end());
	}

	return peaks;
}

bool KeypointDetector::isEdgePeak(int x, int y, Image * image) {
	assert(image);

	float r = EDGE_PEAK_THRESH;
	float b = pow(r + 1, 2) / r;

	float D = image->getPixel(x, y);
	float Dxx = image->getPixel(x + 1, y) + image->getPixel(x - 1, y) - 2*D;
	float Dyy = image->getPixel(x, y + 1) + image->getPixel(x, y - 1) - 2*D;
	float Dxy = ((image->getPixel(x + 1, y + 1) - image->getPixel(x - 1, y + 1))
		     - (image->getPixel(x + 1, y - 1) - image->getPixel(x - 1, y - 1))) / 4.0;

	float TrH = Dxx + Dyy;
	float DetH = Dxx*Dyy - Dxy * Dxy;

	float a = TrH * TrH / DetH;

	
	if (a < b) {
		//printf("Not edge %f %f\n", a, b);
		return false;
	}
	else {
		//printf("Edge\n");
		return true;
	}
}

float KeypointDetector::InterpKeyStep(int x, int y, int s, vector<Image *> & DI,
				 float * dx, float * dy, float * ds) {
	
       
	float Dp[3] = {0}; // first derivative of D with respect to x, y, s

	float Dpp[3][3] = {{0}}; // Hessian of D

	Dp[0] = (DI[s]->getPixel(x+1, y) - DI[s]->getPixel(x-1, y)) / 2.0; // Dx
	Dp[1] = (DI[s]->getPixel(x, y+1) - DI[s]->getPixel(x, y-1)) / 2.0; // Dy
	Dp[2] = (DI[s+1]->getPixel(x, y) - DI[s-1]->getPixel(x, y)) / 2.0; // Ds

	
	/*
	  Hessian (3x3) matrix is defined as follows:
	  Symmetric matrix
 
	  Dxx Dxy Dxs
	  Dyx Dyy Dys
	  Dsx Dsy Dss
	
	*/

        // Dxx
	Dpp[0][0] = (DI[s]->getPixel(x+1, y) + DI[s]->getPixel(x-1, y)
		     - 2.0 * DI[s]->getPixel(x, y));

        // Dyy
	Dpp[1][1] = (DI[s]->getPixel(x, y+1) + DI[s]->getPixel(x, y-1)
		     - 2.0 * DI[s]->getPixel(x, y));

	// Dzz
	Dpp[2][2] = (DI[s+1]->getPixel(x, y) + DI[s-1]->getPixel(x, y)
		     - 2.0 * DI[s]->getPixel(x, y));


	// Dxy = Dyx
	Dpp[0][1] = Dpp[1][0] = (DI[s]->getPixel(x+1, y+1) - DI[s]->getPixel(x-1, y+1)
				 - DI[s]->getPixel(x+1, y-1) + DI[s]->getPixel(x-1, y-1)) / 4.0;
		     
        // Dxs = Dsx
	Dpp[0][2] = Dpp[2][0] = (DI[s+1]->getPixel(x+1, y) - DI[s+1]->getPixel(x-1, y)
				 - DI[s-1]->getPixel(x+1, y) + DI[s-1]->getPixel(x-1, y)) / 4.0;

	// Dys = Dsy
	Dpp[1][2] = Dpp[2][1] = (DI[s+1]->getPixel(x, y+1) - DI[s+1]->getPixel(x, y-1)
				 - DI[s-1]->getPixel(x, y+1) + DI[s-1]->getPixel(x, y-1)) / 4.0;


	float invDpp[3][3];

	mInv33(Dpp, invDpp);

	// Solve for delta positions
	*dx = 0;
	for (int i = 0; i < 3; i++)
		*dx -= invDpp[0][i] * Dp[i];

	*dy = 0;
	for (int i = 0; i < 3; i++)
		*dy -= invDpp[1][i] * Dp[i];

	*ds = 0;
	for (int i = 0; i < 3; i++)
		*ds -= invDpp[2][i] * Dp[i];

	//printf("Interp: %f %f %f\n", *dx, *dy, *ds);

	float val = DI[s]->getPixel(x, y);
	val += 0.5 * (Dp[0] * *ds + Dp[1] * *dy + Dp[2] * *ds);

	return fabs(val);
}

bool KeypointDetector::InterpKey(int x, int y, int s, vector<Image *> & DImages,
	       float * fx, float * fy, float * fs) {

	bool addkey = true;

	
	int moves_left = 5;
	int tx = x;
	int ty = y;
	int ts = s;
		
	float dx, dy, ds;

	bool updated;
	
	do {
		moves_left--;
		updated = false;

		float val = InterpKeyStep(tx, ty, ts, DImages, &dx, &dy, &ds);
		
		if (val < (float) LOW_CONTRAST_THRESH / (float) SCALESPEROCTAVE) {
			addkey = false;
			continue;
		}
		
	
		if (dx > 0.6 && tx < DImages[0]->width - 3) {
			tx++;
			updated = true;
		}
		else if (dx < -0.6 && tx > 3) {
			tx--;
			updated = true;
		}


		if (dy > 0.6 && ty < DImages[0]->height - 3) {
			ty++;
			updated = true;
		}
		else if (dy < -0.6 && ty > 3) {
			ty--;
			updated = true;
		}

	} while (moves_left > 0 && updated);

	if (addkey && fabs(dx) < 1.5
	    && fabs(dy) < 1.5 && fabs(ds) < 1.5) {
		*fx = tx + dx;
		*fy = ty + dy;
		*fs = ts + ds;

		//printf ("Interp: %f %f %f\n", *fx - x, *fy - y, *fs - s);
		return true;
	}

	return false;
}

vector<Keypoint *> KeypointDetector::FindPeaksScales(int octave, vector<Image *> & DImages) {
	vector<Keypoint *> peaks;


        //printf("DImages.size %d\n", DImages.size());

        Image * kpfound = new Image(DImages[0]->width, DImages[0]->height);
	
	float contrast_thresh = (float) LOW_CONTRAST_THRESH / (float) SCALESPEROCTAVE;
	
	if (INTERP_KEYS)
		contrast_thresh *= 0.8;

	for (unsigned int s = 1; s < DImages.size() - 1; s++) {

		for (int y = BORDER; y < DImages[0]->height - BORDER; y++) {
	
			for (int x = BORDER; x < DImages[0]->width - BORDER; x++) {
			
		
				//printf("FindPeaks(): %d %d %d\n", x, y, s); 

				// check for low contrast
				if (fabs(DImages[s]->getPixel(x, y))
				    <= contrast_thresh)
					continue;
				
				//printf("\nK\t%d\t%d\t%d\t", x, y, s);

				//printf("C ");


                                if (kpfound->getPixel(x, y) == 1)
					continue;

				if (!isPeak(DImages[s-1], DImages[s], DImages[s+1], x, y))
				       	continue;

				//printf("P ");

				//printf("(%d, %d, %d) = %f\n", x, y, s, DImages[s]->getPixel(x, y));

				if (isEdgePeak(x, y, DImages[s]))
					continue;

				//printf("E ");

				float fx = x;
				float fy = y;
				float fs = s;

				if (INTERP_KEYS) {
					if (!InterpKey(x, y, s, DImages, &fx, &fy, &fs))
					    continue;
				}

				Keypoint * peak = new Keypoint(EPCALEN);
				peak->x = fx * pow(2.0, octave);
				peak->y = fy * pow(2.0, octave);
				peak->octave = octave;
				peak->gscale = SIGMA * pow(2.0, octave + fs / (float) SCALESPEROCTAVE);
				if (DOUBLE_BASE_IMAGE_SIZE)
					peak->gscale /= 2.0;

				peak->scale = s;
				peak->fscale = fs;
				peak->sx = fx;
				peak->sy = fy;
				
				
				peaks.push_back(peak);
                                kpfound->setPixel((int)(fx + 0.5), (int)(fy + 0.5), 1);

				scalestat[octave * SCALESPEROCTAVE + s]++;

				//printf("Found Key at (%d, %d, %d) -> (%f, %f, %f)\n",
				//       x, y, s, fx - x, fy - y, fs - s);

			}
		}
	}

        delete kpfound;

	return peaks;
}


bool KeypointDetector::GetPixOrientation(int x, int y, Image * image, float * m,  float * theta) {
	assert(image);

	if (x < 1 || y < 1 || x > image->width - 2 || y > image->height - 2)
		return false;

	float a = image->getPixel(x + 1, y) - image->getPixel(x - 1, y);
	float b = image->getPixel(x, y - 1) - image->getPixel(x, y + 1);
	
	*m = sqrt(a*a +  b*b);
        
	*theta = atan2(b, a);

	return true;
}

vector<Keypoint *> KeypointDetector::FindOrientation(vector<Keypoint *> & kps, vector<vector<Image *> > & LOctaves) {
	vector<Keypoint * > newkps;

        printf("FindOrientation(): Finding keypoint orientations\n");
        
	for (unsigned int i = 0; i < kps.size(); i++) {
		float sigma = 1.5 * pow(2.0, (kps[i]->fscale)/(float) SCALESPEROCTAVE) * SIGMA;
		vector<vector<float> > gmat = GaussianKernel2D(sigma);
		int c = gmat.size()/2;

		float thetas[NUMORIENTATIONS] = {0};

		for (int y = -c; y <= c; y++) {
			for (int x = -c; x <= c; x++) {
				float m, theta;
                                if (sqrt((float) x*x + y*y) > 3.0 * sigma)
                                        continue;
                                
				bool valid = GetPixOrientation((int) (kps[i]->sx + x + 0.5), (int) (kps[i]->sy + y + 0.5),
							       LOctaves[kps[i]->octave][kps[i]->scale],
							       &m, &theta);
				if (valid) {
					float degree = theta / PI * 180.0 + 180.0;
					int index = ((int) (degree / DEGPERBIN));
					thetas[index] += m * gmat[y + c][x + c];
				}
					
			}
		}

                for (int j = 0; j < 6; j++)
                        SmoothHistogram(thetas, NUMORIENTATIONS);
                        
                float maxval = maxvec(thetas, NUMORIENTATIONS);

		
                for (int j = 0; j < NUMORIENTATIONS; j++) {
                        if (thetas[j] < maxval * NEW_KP_OR_THRESH)
                                continue;

                        int indexa = mod(j - 1, NUMORIENTATIONS);
                        int indexb = j;
                        int indexc = mod(j + 1, NUMORIENTATIONS);

                        float thetaa = thetas[indexa];
                        float thetab = thetas[indexb];
                        float thetac = thetas[indexc];

                        if (!(thetab > thetaa && thetab > thetac))
                                continue;
                        
                        float maxp = getMaxParabola(-1, thetaa, 0, thetab, 1, thetac);
			//printf("maxp %f : %f %f %f\n", maxp, thetaa, thetab, thetac);

                        if (thetas[j] == maxval)
                                kps[i]->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NUMORIENTATIONS - PI;
                        else {
                                Keypoint * kp = new Keypoint(EPCALEN);
				*kp = *kps[i];
                                kp->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NUMORIENTATIONS - PI;
                                newkps.push_back(kp);
                        }
                                
		 }
	}

	return newkps;
}

void KeypointDetector::MakeKeypointPatch(Keypoint * key, Image * blur, float v[GPLEN])
{
	assert(key);
	assert(blur);

	int patchsize;
	int iradius;
	float  sine, cosine;
	float sizeratio;
	
	float scale = SIGMA * pow(2.0, key->fscale / (float) SCALESPEROCTAVE);


	// sampling window size
	patchsize = (int) (PatchMag * scale);
	
	// make odd
	patchsize /= 2;
	patchsize = patchsize * 2 + 1;
	
	if (patchsize < PatchSize)
		patchsize = PatchSize;
	
	sizeratio = (float) patchsize / (float) PatchSize;
	
	Image * win = new Image(patchsize, patchsize);
	
	sine = sin(key->ori);
	cosine = cos(key->ori);
	
	iradius = patchsize / 2;
	
	/* Examine all points from the gradient image that could lie within the
	   index square.
	*/
	
	//fprintf(stderr, "Scale %f  %d\n", scale, patchsize);
	
	//fprintf(stderr, "Key Patch of orientation %f\n", key->ori);
	for (int y = -iradius; y <= iradius; y++)
		for (int x = -iradius; x <= iradius; x++) {
			
			// calculate sample window coordinates (rotated along keypoint)
			float cpos = ((cosine * x + sine * (float) y) + key->sx);
			float rpos = ((-sine * x + cosine * (float) y) + key->sy);
			
			win->setPixel(x + iradius, y + iradius, blur->getPixelBI(cpos, rpos));
			
			//fprintf(stderr, "  (%d, %d) -> (%f, %f)\n", j, i, cpos, rpos);
		}
	
	int count = 0;

	for (int y = 1; y < PatchSize - 1; y++) {
		for (int x = 1; x < PatchSize - 1; x++) {
			float x1 = win->getPixelBI((float) (x+1) * sizeratio, (float) y * sizeratio);
			float x2 = win->getPixelBI((float) (x-1) * sizeratio, (float) y * sizeratio);
			
			float y1 = win->getPixelBI((float) x * sizeratio, (float) (y + 1) * sizeratio);
			float y2 = win->getPixelBI((float) x * sizeratio, (float) (y - 1) * sizeratio);
			
			// would need to divide by 2 (span 2 pixels), but we normalize anyway
			// so it's not necessary
			float gx = x1 - x2;
			float gy = y1 - y2;


			v[count] = gx;
			v[count + 1] = gy;

			count += 2;
		}
		//fprintf(stderr, "\n");
	}
	
}



void KeypointDetector::MakeKeypointPCA(Keypoint * key, Image * blur) {
	assert(key);
	assert(blur);

	float v[GPLEN];
	MakeKeypointPatch(key, blur, v);

	Normvec(v, GPLEN);

	
	for (unsigned int i = 0; i < GPLEN; i++) {
		v[i] -= avgs[i];
	}

	for (int ldi = 0; ldi < EPCALEN; ldi++) {
		key->ld[ldi] = 0;

		for (int x = 0; x < GPLEN; x++)
			key->ld[ldi] += eigs[ldi][x] * v[x];
	}
}

void KeypointDetector::ComputeLocalDescr(vector<Keypoint *> & kps, vector<vector<Image *> > & LOctaves) {
	printf("ComputeLocalDescr(): Computing local descriptors.\n");

	for (unsigned int i = 0; i < kps.size(); i++) {
                MakeKeypointPCA(kps[i], LOctaves[kps[i]->octave][kps[i]->scale]);
	}
}


vector<vector<Image *> > KeypointDetector::BuildGaussianOctaves(Image * image) {
	assert(image);

	vector<vector<Image *> > octaves;
	
	printf("BuildGaussianOctaves(): Base image dimension is %dx%d\n", image->width, image->height);
	int dim = min(image->height, image->width);
	int numoctaves = int (log((double) dim) / log(2.0)) - 2;
	
	numoctaves = min(numoctaves, MAXOCTAVES);
	
	printf("BuildGaussianOctaves(): Building %d octaves\n", numoctaves);

        // start with initial source image
	Image * timage = image->clone();

	for (int i = 0; i < numoctaves; i++) {
		printf("Building octave %d of dimesion (%d, %d)\n", i, timage->width, timage->height);
		vector<Image *> scales = BuildGaussianScales(timage);
		octaves.push_back(scales);

		// halve the image size for next iteration
		Image * simage  = scales[SCALESPEROCTAVE]->halfSizeImage();
		delete timage;
		timage = simage;

	}

	delete timage;

	return octaves;
}

vector<Image *> KeypointDetector::BuildGaussianScales(Image * image) {
	assert(image);

	vector<Image *> GScales;

	//printf("Building scales of dimension (%d, %d)\n", image->width, image->height);

	double k = pow(2, 1.0/(float)SCALESPEROCTAVE);
	
	GScales.push_back(image->clone());

	for (int i =  1; i < SCALESPEROCTAVE + 3; i++) {
		Image * dst = new Image(image->width, image->height);
		//Image * dst = GScales[GScales.size() - 1]->clone();

		// 2 passes of 1D on original
		float sigma1 = pow(k, i - 1) * SIGMA;
		float sigma2 = pow(k, i) * SIGMA;
		float sigma = sqrt(sigma2*sigma2 - sigma1*sigma1);
		
		//printf("Blur %f\n", sigma);
		BlurImage(GScales[GScales.size() - 1], dst, sigma);
		

		GScales.push_back(dst);
	}

	return GScales;
}


vector<vector<Image *> > KeypointDetector::BuildDOGOctaves(vector<vector<Image *> > & GOctaves) {
	vector<vector<Image *> > DOctaves;

	printf("BuildDOGs(): Building %d DOG Octaves\n", GOctaves.size());

	for (unsigned int i = 0; i < GOctaves.size(); i++) {
		DOctaves.push_back(BuildDOGScales(GOctaves[i]));
	}

	return DOctaves;
}

vector<Image *> KeypointDetector::BuildDOGScales(vector<Image *> & LImages) {
	vector<Image *> DImages;

	for (unsigned int i = 1; i < LImages.size(); i++) {
		Image * dog = new Image(LImages[i]->width, LImages[i]->height);
		LImages[i]->sub(LImages[i-1], dog);
		DImages.push_back(dog);
	}

	return DImages;
}

Image * KeypointDetector::ScaleInitImage(Image * image) {
	assert(image);

	Image * dst;

	if (DOUBLE_BASE_IMAGE_SIZE) {
		Image * im = image->doubleSizeImage2();
		dst = new Image(im->width, im->height);
		//dst = im->clone();
		double sigma = sqrt(SIGMA * SIGMA - INITSIGMA * INITSIGMA * 4);
		//printf("Init Sigma: %f\n", sigma);
		BlurImage(im, dst, sigma);
		
		delete im;
	} else {
		dst = new Image(image->width, image->height);
		double sigma = sqrt(SIGMA * SIGMA - INITSIGMA * INITSIGMA);
		//printf("Init Sigma: %f\n", sigma);
		BlurImage(image, dst, sigma);
	}

	return dst;
}

vector <Keypoint *> KeypointDetector::FindKeypoints(Image * im) {
	assert(im);
	Image * image = ScaleInitImage(im);

	vector<vector<Image *> > GOctaves = BuildGaussianOctaves(image);
	vector<vector<Image *> > DOctaves = BuildDOGOctaves(GOctaves);

	vector<Keypoint *> peaks = FindPeaksOctaves(DOctaves);
	vector<Keypoint *> newkps = FindOrientation(peaks, GOctaves);
	peaks.insert(peaks.begin(), newkps.begin(), newkps.end());

	/*
	for (unsigned int i = 0; i < peaks.size(); i++) {
		Keypoint * k = peaks[i];
		printf("Found key at (%d, %d, o %f, s %f)\n", (int) k->x, (int) k->y,
			 k->ori, k->gscale);
	}
	*/

        ComputeLocalDescr(peaks, GOctaves);
	printf("Done\n");

	delete image;


	//printstat();

	return peaks;
}


void KeypointDetector::RecalcKeyIndices(vector <Keypoint *> & keys) {
	float log2 = log(2);

	for (unsigned int i = 0; i < keys.size(); i++) {
		
		Keypoint * k = keys[i];

		/* need to recalc:
		   sx, sy, octave, scale, fscale
		*/
		//printf("gscale %f\n", k->gscale);
		
		float tmp = log(k->gscale/SIGMA) / log2 + 1.0;
		//printf("  tmp %f\n", tmp);

		k->octave = (int) tmp;

		//printf("  octave %d\n", k->octave);
		k->fscale = (tmp - k->octave) * (float) SCALESPEROCTAVE;
		//printf("  fscale %f\n", k->fscale);

		k->scale = (int) round(k->fscale);

		if (k->scale == 0 && k->octave > 0) {
			k->scale = SCALESPEROCTAVE;
			k->octave--;
			k->fscale += SCALESPEROCTAVE;
		}
			
		//printf("  scale %d\n", k->scale);

		k->sx = k->x / pow(2.0, k->octave);
		k->sy = k->y / pow(2.0, k->octave);

		if (DOUBLE_BASE_IMAGE_SIZE) {
			k->sx *= 2.0;
			k->sy *= 2.0;
			k->x *= 2.0;
			k->y *= 2.0;
		}
		
		//printf("(%f, %f) -> (%f, %f)\n", k->x, k->y, k->sx, k->sy);
	}
}

void KeypointDetector::RecalcKeys(Image * im, vector <Keypoint *> & keys) {
	assert(im);
	Image * image = ScaleInitImage(im);

	vector<vector<Image *> > GOctaves = BuildGaussianOctaves(image);
	//vector<vector<Image *> > DOctaves = BuildDOGOctaves(GOctaves);

	RecalcKeyIndices(keys);

        ComputeLocalDescr(keys, GOctaves);

	delete image;

}

void KeypointDetector::MakeLocalPatch(Keypoint * key, Image * blur, int windowsize)
{
        assert(key);
        assert(blur);
                                                                                                                       
        int patchsize;
        int iradius;
        float sine, cosine;
        float sizeratio;
                                                                                                                       
        float scale = SIGMA * pow(2.0, key->fscale / (float) SCALESPEROCTAVE);
                                                                                                                       
                                                                                                                       
        // sampling window size
        patchsize = (int) (PatchMag * scale);
                                                                                                                       
        // make odd
        patchsize /= 2;
        patchsize = patchsize * 2 + 1;
                                                                                                                       
        /* The following two lines are not entirely correct.
           In theory, patchsize should never be smaller than PatchSize.
           If it is, then PatchMag isn't set large enough, or PatchSize is set too small.
           It's kept here for compatibility, and is a bug.
           It should eventually be removed.
        */
        if (patchsize < PatchSize)
                patchsize = PatchSize;
                                                                                                                       
        sizeratio = (float) patchsize / (float) PatchSize;
                                                                                                                       
        key->patch = new Image(windowsize, windowsize);
                                                                                                                       
        sine = sin(key->ori);
        cosine = cos(key->ori);
                                                                                                                       
        iradius = windowsize / 2;
                                                                                                                       
                                                                                                                       
        //fprintf(stderr, "Scale %f  %d\n", scale, patchsize);
                                                                                                                       
        //fprintf(stderr, "Key Patch of orientation %f\n", key->ori);
        for (int y = -iradius; y <= iradius; y++) {
                for (int x = -iradius; x <= iradius; x++) {
                                                                                                                       
                        // calculate sample window coordinates (rotated along keypoint)
                        float cpos = ((cosine * x * sizeratio
                                       + sine * (float) y * sizeratio) + key->sx);
                        float rpos = ((-sine * x * sizeratio
                                       + cosine * (float) y * sizeratio) + key->sy);
                                                                                                                       
                        key->patch->setPixel(x + iradius, y + iradius, blur->getPixelBI(cpos, rpos));
                                                                                                                       
                        //fprintf(stderr, "  (%d, %d) -> (%f, %f)\n", j, i, cpos, rpos);
                }
        }
 
}
 


void KeypointDetector::ComputeLocalPatches(vector<Keypoint *> & kps, vector<vector<Image *> > & LOctaves, int patchsize) {
        printf("ComputeLocalDescr(): Computing local patches.\n");
 
        for (unsigned int i = 0; i < kps.size(); i++) {
                assert(kps[i]->octave >= 0 && kps[i]->octave < (int) LOctaves.size());
                assert(kps[i]->scale >= 0 && kps[i]->scale < (int) LOctaves[kps[i]->octave].size());
 
 
                MakeLocalPatch(kps[i], LOctaves[kps[i]->octave][kps[i]->scale], patchsize);
        }
}


void KeypointDetector::writeKeysToFile(vector<Keypoint *> & keys, char * fn) {
	assert(fn);

	printf("Writing to %s\n", fn);
                                                                                
        FILE * fp = fopen(fn, "wb");

	/* Output total number of keypoints and VecLength. */
        fprintf(fp, "%d %d\n", keys.size(), EPCALEN);
                                                                                
        for (unsigned int i = 0; i < keys.size(); i++) {

		
		if (DOUBLE_BASE_IMAGE_SIZE) {
			// don't divide gscale by 2 here because we already did.
			fprintf(fp, "%4.2f %4.2f %4.3f %4.4f", keys[i]->y/2.0, keys[i]->x/2.0, keys[i]->gscale,
				keys[i]->ori);
		} else
			fprintf(fp, "%4.2f %4.2f %4.3f %4.4f", keys[i]->y, keys[i]->x,
				keys[i]->gscale, keys[i]->ori);
                                                                                
                for (int j = 0; j < EPCALEN; j++) {
                        if (j % 12 == 0)
                                fprintf(fp, "\n");
                        fprintf(fp, " %d", (int) keys[i]->ld[j]);
                }
                                                                                
                                                                                
                fprintf(fp, "\n");
        }
                                                                                
                                                                                
        fflush(fp);
        fclose(fp);

}


vector<Keypoint *> KeypointDetector::readKeysFromFile(char * filename)
{
	assert(filename);
	
	int num, pcalen;
	
	FILE * fp;
	vector<Keypoint *> keys;
	
	fp = fopen (filename, "rb");
	
	if (! fp) {
		fprintf(stderr, "Could not open file: %s", filename);
		exit(1);
	}         
	
	fscanf(fp, "%d %d", &num, &pcalen);
	

	for (int i = 0; i < num; i++) {
		/* Allocate memory for the keypoint. */
		Keypoint * key = new Keypoint(pcalen);
		
		if (fscanf(fp, "%f %f %f %f", &(key->y), &(key->x), &(key->gscale),
			   &(key->ori)) != 4) {
			fprintf(stderr, "Invalid keypoint file format.");
			exit(1);
		}
                                                                                
		for (int j = 0; j < pcalen; j++) {
			float fval;
			fscanf(fp, "%f", &fval);
			key->ld[j] = fval;
		}
		
		keys.push_back(key);
		
	}
	
	fclose(fp);
	
	return keys;
}



vector<Keypoint *> KeypointDetector::readPatchesFromFile(char * filename)
{
	assert(filename);
	
	int num, pcalen;
	
	FILE * fp;
	vector<Keypoint *> keys;
	
	fp = fopen (filename, "rb");
	
	if (! fp) {
		fprintf(stderr, "Could not open file: %s", filename);
		exit(1);
	}         
	
	fscanf(fp, "%d %d", &num, &pcalen);
	
	int sqrtlen = (int) sqrt(pcalen);

	if (sqrtlen * sqrtlen != pcalen) {
		fprintf(stderr, "Invalid patch file - dimensions incorrect: %d\n", pcalen);
		exit(1);
	}

	for (int i = 0; i < num; i++) {
		/* Allocate memory for the keypoint. */
		Keypoint * key = new Keypoint(0);
		
		if (fscanf(fp, "%f %f %f %f", &(key->y), &(key->x), &(key->gscale),
			   &(key->ori)) != 4) {
			fprintf(stderr, "Invalid keypoint file format.");
			exit(1);
		}

		key->patch = new Image(sqrtlen, sqrtlen);

		for (int y = 0; y < sqrtlen; y++) {
			for (int x = 0; x < sqrtlen; x++) {
				float fval;
				fscanf(fp, "%f", &fval);
				key->patch->setPixel(x, y, fval);
			}
		}
		
		keys.push_back(key);
		
	}
	
	fclose(fp);
	
	return keys;
}


void KeypointDetector::getPatches(Image * im, vector <Keypoint *> & keys, int patchsize) {
        assert(im);
        Image * image = ScaleInitImage(im);
 
        vector<vector<Image *> > GOctaves = BuildGaussianOctaves(image);
 
        RecalcKeyIndices(keys);
 
        ComputeLocalPatches(keys, GOctaves, patchsize);
 
        delete image;
 
}

 
void KeypointDetector::writePatchesToFile(vector<Keypoint *> & keys, char * fn, int patchsize) {
        assert(fn);
 
        printf("Writing to %s\n", fn);
                                                                                 
        FILE * fp = fopen(fn, "wb");
 
        /* Output total number of keypoints and VecLength. */
        fprintf(fp, "%d %d\n", keys.size(), patchsize * patchsize);
                                                                                 
        for (unsigned int i = 0; i < keys.size(); i++) {
 
                 
                if (DOUBLE_BASE_IMAGE_SIZE) {
                        // don't divide gscale by 2 here because we already did.
                        fprintf(fp, "%4.2f %4.2f %4.2f %4.3f\n",
                                keys[i]->y/2.0, keys[i]->x/2.0, keys[i]->gscale,
                                keys[i]->ori);
                } else
                        fprintf(fp, "%4.2f %4.2f %4.2f %4.3f\n",
                                keys[i]->y, keys[i]->x,
                                keys[i]->gscale, keys[i]->ori);
                           
                for (int y = 0; y < patchsize; y++) {
                        for (int x = 0; x < patchsize; x++)
                                fprintf(fp, "%0.4f ", keys[i]->patch->getPixel(x, y));
                        fprintf(fp, "\n");
 
                }
                 
        }
                                                                                 
                                                                                 
        fflush(fp);
        fclose(fp);
 
}
