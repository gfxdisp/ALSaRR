/*
 * MATLAB MEX function for writing EXR images.
 * 
 * Only supports EXR images with half-precision floating-point data.
 *
 */


#include "ImathBox.h"
#include "ImfRgba.h"
#include "ImfRgbaFile.h"
#include "ImfArray.h"
#include "ImfChannelList.h"
#include "ImfPixelType.h"
#include "Iex.h"

#include "mex.h" 
#include "matrix.h"

using namespace Imf;
using namespace Imath;
using namespace Iex;

/*
 * Check inputs
 * Two or three input arguments, last is a string
 * no output arguments
 * If a mask is specified, it must be 2D and the same size as the image
 * 
 */
void checkInputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if (nlhs != 0)
		mexErrMsgTxt("Too many output arguments");

	if (nrhs < 2 || nrhs > 3)
		mexErrMsgTxt("Incorrect number of input arguments.");

	// Check that last argument is a string
	if (mxIsChar(prhs[nrhs-1]) != 1)
		mexErrMsgTxt("Last argument must be a string.");

	if (mxGetM(prhs[nrhs-1]) != 1)
		mexErrMsgTxt("Filename must be a row vector.");

	// Check the number of dimensions
	int ndims = mxGetNumberOfDimensions(prhs[0]);
	if (ndims < 2 || ndims > 3)
		mexErrMsgTxt("image must be 2D or 3D.");

	// Make sure that image is double precision
	if (mxGetElementSize(prhs[0]) != 8)
		mexErrMsgTxt("image must be double precision.");

	// If a mask is specified, it must be the same size as the image
	if (nrhs == 3) {

		int ndmask = mxGetNumberOfDimensions(prhs[1]);
		if (ndmask != 2) 
			mexErrMsgTxt("mask must be 2D");

		int ydim = mxGetM(prhs[0]);
		int xdim = mxGetN(prhs[0]);
		if (ndims == 3) {
			xdim = mxGetN(prhs[0])/3;
		}

		int ymask = mxGetM(prhs[1]);
		int xmask = mxGetN(prhs[1]);

		if (ydim != ymask || xdim != xmask)
			mexErrMsgTxt("mask and image must be the same size.");

		if (mxGetElementSize(prhs[1]) != 8)
			mexErrMsgTxt("mask must be double precision.");
	}

	return;
}

/*
 * Write an EXR file.
 * Code follows examples from ReadingAndWritingImageFiles.pdf, found
 * here:
 * http://www.openexr.com/ReadingAndWritingImageFiles.pdf
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	checkInputs(nlhs, plhs, nrhs, prhs);
	char *filename = mxArrayToString(prhs[nrhs-1]);

	try {
		// Get image size
		int ndims = mxGetNumberOfDimensions(prhs[0]);
		int ydim = mxGetM(prhs[0]);
		int xdim = mxGetN(prhs[0]);
		if (ndims == 3) 
			xdim = mxGetN(prhs[0])/3;
		
		double *img = mxGetPr(prhs[0]);
		struct Rgba *px = new struct Rgba[xdim*ydim];
		
		int step = 0;
		if (ndims == 3)
			step = xdim*ydim;

		RgbaChannels mode = WRITE_RGBA;

		/*
		 * alpha specified as second argument
		 */
		if (nrhs == 3) {
			double *msk = mxGetPr(prhs[1]);

			for (int i = 0; i < ydim; ++i) {
				for (int j = 0; j < xdim; ++j) {
					int k = j*ydim+i;
					px[i*xdim+j].r = img[k];
					px[i*xdim+j].g = img[k+step];
					px[i*xdim+j].b = img[k+2*step];
					px[i*xdim+j].a = msk[k];
				}
			}

			if (ndims == 2)
				mode = WRITE_YA;

		/*
		 * No alpha specified
		 */
		} else {

			for (int i = 0; i < ydim; ++i) {
				for (int j = 0; j < xdim; ++j) {
					int k = j*ydim+i;
					px[i*xdim+j].r = img[k];
					px[i*xdim+j].g = img[k+step];
					px[i*xdim+j].b = img[k+2*step];
					px[i*xdim+j].a = 1.0f;
				}
			}

			mode = WRITE_RGB;
			if (ndims == 2)
				mode = WRITE_Y;

		}

		RgbaOutputFile file(filename, xdim, ydim, mode);
		file.setFrameBuffer(px, 1, xdim);
		file.writePixels(ydim);
		mxFree(filename);
	
		delete[] px; 
		px = NULL;

	} catch (const std::exception &exc) {
		mxFree(filename);
		mexErrMsgTxt(exc.what());
	}

	return;

} 

