
/* -*-  Mode:C++; c-basic-offset:8; tab-width:8; indent-tabs-mode:t -*- */

/*

Author: Yan Ke
Dec 2003

*/

#include <stdio.h>

#include "config.h"
#include "image.h"
#include "keypoint.h"

int main(int argc, char* argv[])
{
	if (argc != 4 ) {
                printf("Usage: %s gpcavects.txt image.pgm image.keys\n", argv[0]);
                return -1;
        }

	char * gpcafn = argv[1];
	char * imagefn = argv[2];
	char * keysfn = argv[3];

	Image * im = new Image(imagefn);

	printf("Image size (%d, %d)\n", im->width, im->height);
	
	KeypointDetector kpd(gpcafn);

	vector<Keypoint *> keys = kpd.FindKeypoints(im);
	
	printf("Found %d keys\n", keys.size());
	
	kpd.writeKeysToFile(keys, keysfn);

	delete(im);
	
	return 0;
}
