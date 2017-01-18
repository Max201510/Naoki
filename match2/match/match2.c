/************************************************************************
Demo software: Invariant keypoint matching.
Author: David Lowe

match.c:
This file contains a sample program to read images and keypoints, then
   draw lines connecting matched keypoints.


Modified by: Yan Ke
December 2003
Adds matching of PCA-SIFT local descriptors.
Adds option of matching nearest neighbor with and without looking
at second nearest neighbor.

*************************************************************************/

#include <time.h>
#include <sys/time.h>
#include "defs.h"

int find=0;
int dsq;
int min=0, max=0;
int p1Max=0, p2Max=0;
int HashSim = 0, cc =0; 
int LIP = 0, LIP2 = 0;
int *temp=0,*temp2=0;
int i, j,f;
int number = 0,k2Number=0;
int count = 0;
int hash;
int hash1[3000][20];
int hash2[3000][20];
struct keyDist KDistance[3000];
struct pair NNpair[3000];
struct pair2 FinalPair[3000];
Keypoint k,k2, match;
float normalize, normalize2, NV[5000][20];
Image result;
char image1[20], image2[20], IK1[20],IK2[20];


struct timeval startTime, endTime; 
/* -------------------- Local function prototypes ------------------------ */

int Hash(float value);
void FindMatches(Image im1, Keypoint keys1, Image im2, Keypoint keys2);
Image CombineImagesVertically(Image im1, Image im2);
void SimCalculation(Keypoint keys1, Keypoint keys2);
void SecondKeyCheck(Keypoint k, Keypoint keys2, int LIP);
void HashCalculation(Keypoint keys2);
void UPFile(int count,float Ksim);
/*----------------------------- Routines ----------------------------------*/

/* Top level routine.  Read PGM images and keypoints from files given
   in command line arguments, then call FindMatches.
*/
int main (int argc, char **argv)
{


    int arg = 0;
    Image im1 = NULL, im2 = NULL;
    Keypoint k1 = NULL, k2 = NULL;


    while (++arg < argc) {
      if (! strcmp(argv[arg], "-im1")) 
	{im1 = ReadPGMFile(argv[++arg]);
	strcpy(image1, argv[arg]);
	}
      else if (! strcmp(argv[arg], "-im2")) {
	im2 = ReadPGMFile(argv[++arg]);
	strcpy(image2, argv[arg]);
	}
      else if (! strcmp(argv[arg], "-k1")){
	k1 = ReadKeyFile(argv[++arg]);
	strcpy(IK1, argv[arg]);
	}
      else if (! strcmp(argv[arg], "-k2")){
	k2 = ReadKeyFile(argv[++arg]);
	strcpy(IK2, argv[arg]);
	}
      else
	FatalError("Invalid command line argument: %s", argv[arg]);
    }
    if (im1 == NULL || im2 == NULL || k1 == NULL || k2 == NULL)
      FatalError("Command line does not specify all images and keys.");

	  

    FindMatches(im1, k1, im2, k2);

    exit(0);
}


void FindMatches(Image im1, Keypoint keys1, Image im2, Keypoint keys2)
{
	float sumOfSim;
    /* Create a new image that joins the two images vertically. */
    result = CombineImagesVertically(im1, im2);
    /*calculate the hash value for Q points and store in to the array*/
    HashCalculation(keys2);


    /* Match the keys in list keys1 to their best matches in keys2.*/
	SimCalculation(keys1,keys2);
	
	
	for(i=0;i<LIP2;i++)
	{
		float preDist=0;
		int KeyNo=0;

		if(FinalPair[i].pairNo == 1)
		{
			//fprintf(stderr,"Checking P = %i ,  Q = %i,  %f\n" ,FinalPair[i].k1No ,i,FinalPair[i].dist[0] );
			count++;
			sumOfSim+= FinalPair[i].dist[0];
			DrawLine(result, (int)FinalPair[i].key1[0]->row, (int) FinalPair[i].key1[0]->col,
				 (int) (FinalPair[i].key2->row + im1->rows), (int) FinalPair[i].key2->col);
		}
		else if(FinalPair[i].pairNo > 1)
		{
			for(j = 0; j< FinalPair[i].pairNo; j++){
				if(preDist > FinalPair[i].dist[j]){
				preDist = FinalPair[i].dist[j];
				KeyNo = j;
				}
			}

			count++;
			sumOfSim+= FinalPair[i].dist[KeyNo];
			//fprintf(stderr,"Checking2 P = %i , Q = %i,  %f\n" , FinalPair[i].k1No,i, FinalPaair[i].dist[KeyNo] );
			DrawLine(result, (int)FinalPair[i].key1[KeyNo]->row, (int) FinalPair[i].key1[KeyNo]->col,
				 (int) (FinalPair[i].key2->row + im1->rows), (int) FinalPair[i].key2->col);
		}

	}

//fprintf(stderr,"\n\nChecking dupilicate for k2: %f\n", sumOfSim);
float Ksim = 0;
Ksim = sumOfSim/count;

    /* Write result image to standard output. */
    WritePGM(stdout, result);
   // fprintf(stderr,"Found %d matches. and similarity between two images is = %f\n\n", count, Ksim);
	
	
    UPFile(count,Ksim);
}

/*make a hash function for all points Q in C2*/
void HashCalculation(Keypoint keys2)
{
	for (k= keys2; k != NULL; k = k->next) {	
		temp = k->descrip;
		p2Max=0;

		for (i = 0; i < DLEN; i++) {
		 p2Max += *temp * *temp;
		 *temp++;
		// fprintf(stderr,"Checking %i %i\n" ,i, p2Max);
		}
	
		p2Max = sqrt(p2Max);
	 	//fprintf(stderr,"Checking %i %i\n" ,i, p2Max);
		temp = k->descrip;
		for (i = 0; i < DLEN; i++) {
			normalize = (float)*temp++/(float)p2Max;
			NV[LIP][i] = normalize;
			hash1[LIP][i] = Hash(normalize);
		}
	
		LIP++;
	}
	LIP=0;
}

/*This is a first loop to check the similarities between the points*/
void SimCalculation(Keypoint keys1, Keypoint keys2)
{
	for (k= keys1; k != NULL; k = k->next) {
		
		temp = k->descrip;
		p1Max=0;
		for (i = 0; i < DLEN; i++) {
		 p1Max += *temp * *temp;
		 *temp++;
		}
		p1Max = sqrt(p1Max);
		
		cc=0;
		LIP2=0;

		SecondKeyCheck(k, keys2,LIP);

	/*something wrong here*/
		//find the nn
		/*NN pair contain the LIP2 ID
		 setcc is the NN to LIP2*/
		if (cc!=0) {
			float geto = 0;
			int setcc=0;
			if(cc >1){
				for (i = 0; i < cc; i++) {
					if(KDistance[i].dist  > geto){
					geto=KDistance[i].dist;
					setcc = i;
					NNpair[LIP].pairNo = KDistance[i].LIPID;
					}
			//fprintf(stderr,"Checking %i, %f\n" ,i, KDistance[i].dist);
				}
			}
			else
			{
				geto=KDistance[0].dist ;
				setcc = 0;
				NNpair[LIP].pairNo = KDistance[0].LIPID;
				//fprintf(stderr,"Checking  %f\n" , KDistance[0].dist);
			}
	
		//pair the NN
		/*k2Number is the key 2 points number*/
		k2Number = NNpair[LIP].pairNo;
		NNpair[LIP].key1 = k;
		NNpair[LIP].key2 = KDistance[setcc].key2;

		//setting the NN to pair 2
		  //putting the LIP2 to the array no so FinalPair store k to k2, k2 = key
			number = FinalPair[k2Number].pairNo;

			FinalPair[k2Number].dist[number] = geto;
			FinalPair[k2Number].key1[number] = k;
			FinalPair[k2Number].key2 = KDistance[setcc].key2;
			FinalPair[k2Number].k1No = LIP;
			FinalPair[k2Number].pairNo++;
		//fprintf(stderr,"Checking the pair no  %i, %i\n" , k2Number, FinalPair[k2Number].pairNo);
		
	      }

	LIP++;
	}//first loop end
		//   realsec= diffsec+diffsub*1e-6; 
	//  fprintf(stderr,"The time for second run is: %f\n",realsec );	

}

/*This is the second loop where it calculate
  the similarities of hash and check the cosine angle similarity*/
void SecondKeyCheck(Keypoint k, Keypoint keys2, int LIP)
{
  float sim;
  float SimDistance = 0;
  int distsq=0;
  int dif;

for (k2= keys2; k2 != NULL; k2 = k2->next) 
	{
	distsq=0;
	temp = k->descrip;
	temp2 = k2->descrip;
   	HashSim = 0;
	p2Max=0;
	SimDistance = 0;
	
		/*calculating the similarity between the H(p) and H(q)*/
		/*
		*/
		for (j = 0; j < DLEN; j++) 
		{
			p2Max += *temp2 * *temp2;
			*temp2++;
			/************0.1****************/
			/*find cosine sim*/
			normalize2 = (float)*temp++/(float)p1Max;
    			sim =  NV[LIP2][j] * normalize2;
			SimDistance += sim;
			/***************0.3**************/
			/*find hash distance*/
			hash = Hash(normalize2);
			
			//if(abs(hash - hash1[LIP2][j]) <=1 )
			if(((hash - hash1[LIP2][j]) <=1 )||((hash - hash1[LIP2][j])>=(-1)))
			  HashSim+=  1;
		}
		distsq = (int) sqrt(distsq);
   
		//p2Max = sqrt(p2Max);
//	fprintf(stderr,"Checking %i, %i, %i %f\n" ,LIP, LIP2, HashSim, SimDistance);	fprintf(stderr,"Checking %i, %i, %i %f\n" ,LIP, LIP2, HashSim, SimDistance);
		/*"HashSim" is the similarity of hashes between LIPs where M is set as 0*/
		if(HashSim >= 20)
		{    
			if (SimDistance > 0.87)		
			{
			KDistance[cc].LIPID = LIP2;
			KDistance[cc].dist = SimDistance;
			KDistance[cc].key2 = k2;
			cc++;
			//fprintf(stderr,"Checking %i, %i, %f\n" ,LIP, LIP2, SimDistance);
			}
	
		}
				

	LIP2++;
	}	

}

/*Hash Function */
/*int Hash(float value)
{	
	int HV = floor((value+1)/0.25);
	//float HV = fmod((value+1), 0.25);
	//fprintf(stderr,"find fmod %f, %f, %f\n" ,value, 0.25, HV);
	return HV;
}*/
int Hash(float x)
{
  int i = (int)x; 
  return i - ( i > x ); 
}

/* Return a new image that contains the two images with im1 above im2.
*/
Image CombineImagesVertically(Image im1, Image im2)
{
    int rows, cols, r, c;
    Image result;

    rows = im1->rows + im2->rows;
    cols = MAX(im1->cols, im2->cols);
    result = CreateImage(rows, cols);

    /* Set all pixels to 0,5, so that blank regions are grey. */
    for (r = 0; r < rows; r++)
      for (c = 0; c < cols; c++)
	result->pixels[r][c] = 0.5;

    /* Copy images into result. */
    for (r = 0; r < im1->rows; r++)
      for (c = 0; c < im1->cols; c++)
	result->pixels[r][c] = im1->pixels[r][c];
    for (r = 0; r < im2->rows; r++)
      for (c = 0; c < im2->cols; c++)
	result->pixels[r + im1->rows][c] = im2->pixels[r][c];
    
    return result;
}

void UPFile(int count,float Ksim){
	FILE *fp;	
	
	if ((fp = fopen("list.txt", "a")) == NULL) {
		printf("file open error!!\n");
		exit(EXIT_FAILURE);
	}
	
	fprintf(fp, " File No 1: %s, File No 2: %s, Found %d matches. similarity = %f\n\n", image1, image2, count, Ksim);
	fclose(fp);
}


