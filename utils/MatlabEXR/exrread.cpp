/*
 * MATLAB MEX function for reading EXR images.
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

using namespace Imf;
using namespace Imath;
using namespace Iex;

/*
 * Check inputs
 * Only one input argument that is a string (row vector of chars)
 * one or two output arguments
 * 
 * These checks were copied from the MATLAB example file revord.c
 */
void checkInputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if (nrhs != 1)
		mexErrMsgTxt("Incorrect number of input arguments.");

	if (nlhs < 1 || nlhs > 2)
		mexErrMsgTxt("Incorrect number of output arguments.");

	if (mxIsChar(prhs[0]) != 1)
		mexErrMsgTxt("Input must be a string.");

	if (mxGetM(prhs[0]) != 1)
		mexErrMsgTxt("Input must be a row vector.");

	return;
}

/*
 * Check that the image is one of the supported formats (1-4 channels of
 * half-precision floating-point data).
 */
int numChannels(const RgbaInputFile &file) {

	const ChannelList &ch = file.header().channels();
	int nchannels = 0;
	for (ChannelList::ConstIterator i = ch.begin(); i != ch.end(); ++i) {
		const Channel &channel = i.channel(); 
		PixelType type = channel.type;

		if (type == HALF) 
			++nchannels;
	}

	if (nchannels > 4) {
		mexWarnMsgTxt("Image has more than 4 channels.");
		nchannels = 4;
	}

	return nchannels;
}


/*
 * Read an EXR file.
 * Code follows examples from ReadingAndWritingImageFiles.pdf, found
 * here:
 * http://www.openexr.com/ReadingAndWritingImageFiles.pdf
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 

	checkInputs(nlhs, plhs, nrhs, prhs);
	char *filename = mxArrayToString(prhs[0]);

	try {
		RgbaInputFile file(filename);
		mxFree(filename);

		// Get the image type
		int nchannels = numChannels(file);
		if (nchannels == 0) {
			mexErrMsgTxt("Unsupported image type.");
		}
		
		Box2i dw = file.dataWindow();

		int xdim  = dw.max.x - dw.min.x + 1;
		int ydim = dw.max.y - dw.min.y + 1;

		Array2D<Rgba> px(ydim,xdim);

		file.setFrameBuffer(&px[0][0]-dw.min.x-dw.min.y*xdim, 1, xdim);
		file.readPixels(dw.min.y, dw.max.y);
		
		int dims[3];
		int sz = xdim*ydim;
		dims[0] = ydim; 
		dims[1] = xdim; 
		dims[2] = 1;

		// Initialize mask if call has 2 arguments on left-hand side,
		// but the image does not have a mask
		if (nlhs == 2 && (nchannels % 2) == 1) {
			plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); 
			double *msk = mxGetPr(plhs[1]);

			for (int i = 0; i < sz; ++i) {
				msk[i] = 1.0;
			}
		}

		dims[2] = 3;
		if (nchannels < 3) {
			sz      = 0;
			dims[2] = 1;
		}
	
		if (nlhs == 1) {

			plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); 
			double *img = mxGetPr(plhs[0]);

			for (int i = 0; i < ydim; ++i) {
				for (int j = 0; j < xdim; ++j) {
					int k = j*ydim + i;

					img[k]      = px[i][j].r;
					img[sz+k]   = px[i][j].g;
					img[2*sz+k] = px[i][j].b;
				}
			}

		/*
		 * Floating point image and alpha channel
		 */
		} else {
			plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); 
			double *img = mxGetPr(plhs[0]);

			dims[2] = 1;
			plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); 
			double *msk = mxGetPr(plhs[1]);

			for (int i = 0; i < ydim; ++i) {
				for (int j = 0; j < xdim; ++j) {
					int k = j*ydim + i;

					img[k]      = px[i][j].r;
					img[sz+k]   = px[i][j].g;
					img[2*sz+k] = px[i][j].b;
 					msk[k]      = px[i][j].a;
				}
			}
		}

	} catch (const std::exception &exc) {
		mxFree(filename);
		mexErrMsgTxt(exc.what());
	}


	return;
} 


