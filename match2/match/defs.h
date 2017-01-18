/************************************************************************
Demo software: Invariant keypoint matching.
Author: David Lowe

defs.h:
This file contains the headers for a sample program to read images and
  keypoints, then perform simple keypoint matching.

Modified by: Yan Ke
Dec 2003

*************************************************************************/

/* From the standard C libaray: */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

/*---------------------------- Constants ---------------------------------*/
/* by Yan Ke.
   1 if we check the absolute distance to the nearest neighbor only.
   This uses SIFT_THRESH and PCA_THRESH.

   0 if we check the relative distance between the closest
   and second closest neighbors.
   This uses SIFT_MULT and PCA_MULT.

   In practice, SIFT works better if we check for relative distances.
   (NEAREST_ONLY = 0). PCA-SIFT works better if we check for absolute
   distances (NEAREST_ONLY = 1).

   Further, NEAREST_ONLY = 0 works better if we're matching exactly two images.

   However, NEAREST_ONLY = 0 works only if we're searching for one instance
   of a keypoint.  If we're looking for multiple instances (for example from
   two of the same objects in the scene), then we must use NEAREST_ONLY = 1.
   Similarly, if we're looking for keypoints in a large database and we don't
   don't know if the second closest neighbor is supposed to be from the same
   or a different keypoint, then we must use NEAREST_ONLY = 1.
*/
#define NEAREST_ONLY 1

#define SIFT_MULT 0.6
#define PCA_MULT 0.6

#define SIFT_THRESH 0.1
#define PCA_THRESH 10000

#define TRUE 1
#define FALSE 0

/*------------------------------ Macros  ---------------------------------*/

#define ABS(x)    (((x) > 0) ? (x) : (-(x)))
#define MAX(x,y)  (((x) > (y)) ? (x) : (y))
#define MIN(x,y)  (((x) < (y)) ? (x) : (y))


/*---------------------------- Structures --------------------------------*/

/* Data structure for a float image.
*/
typedef struct ImageSt {
  int rows, cols;          /* Dimensions of image. */
  float **pixels;          /* 2D array of image pixels. */
  struct ImageSt *next;    /* Pointer to next image in sequence. */
} *Image;


/* Data structure for a keypoint.  Lists of keypoints are linked
   by the "next" field.
*/
typedef struct KeypointSt {
  float row, col;             /* Subpixel location of keypoint. */
  float scale, ori;           /* Scale and orientation (range [-PI,PI]) */
  int *descrip;     /* Vector of descriptor values */
  struct KeypointSt *next;    /* Pointer to next keypoint in list. */
} *Keypoint;

struct keyDist {
	float dist;
	int LIPID;
	Keypoint key2;
};

struct pair {
	int pairNo;
  	Keypoint key1;
	Keypoint key2;
};
struct pair2 {
	int pairNo;
	int k1No;
	float dist[100];
  	Keypoint key1[100];
	Keypoint key2;
};





/* whether we use PCA-SIFT (=1) or SIFT keys(=0) */
extern int PCAKEYS;

/* length of descriptor vector */
extern int DLEN;

/*-------------------------- Function prototypes -------------------------*/
/* These are prototypes for the external functions that are shared
   between files.
*/

/* From util.c */
void FatalError(char *fmt, ...);
Image CreateImage(int rows, int cols);
Image ReadPGMFile(char *filename);
Image ReadPGM(FILE *fp);
void WritePGM(FILE *fp, Image image);
void DrawLine(Image image, int r1, int c1, int r2, int c2);
Keypoint ReadKeyFile(char *filename);
Keypoint ReadKeys(FILE *fp);
