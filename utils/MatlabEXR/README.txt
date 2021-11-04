
Install ilmbase and openexr to default location /usr/local

Compile from within MATLAB:

mex exrread.cpp  -lIlmImf -I/usr/local/include/OpenEXR -L/usr/local/lib
mex exrwrite.cpp -lIlmImf -I/usr/local/include/OpenEXR -L/usr/local/lib
mex exrinfo.cpp  -lIlmImf -I/usr/local/include/OpenEXR -L/usr/local/lib

If OpenEXR is installed in another location (such as /sw), try these commands:

mex exrread.cpp  -lIlmImf -lIex -lImath -lHalf -I/sw/include/OpenEXR -L/sw/lib
mex exrwrite.cpp -lIlmImf -lIex -lImath -lHalf -I/sw/include/OpenEXR -L/sw/lib
mex exrinfo.cpp  -lIlmImf -lIex -lImath -lHalf -I/sw/include/OpenEXR -L/sw/lib


-- Usage --

Read image
>> im = exrread(filename);

Read image and alpha channel
>> [im,mask] = exrread(filename);

Image can be 1 or 3 channels of floating-point data
Mask will be 1.0 if there is no alpha channel in the file

Write image
>> exrwrite(im,filename)

Write image with mask in alpha channel
>> exrwrite(im,mask,filename)

Image can be 1 or 3 channels
Mask must be 1 channel the same size as the image


-- Details --

Does not support EXR images with uint16 data or float data
All data is returned as MATLAB type double
Negative, NaN, or Inf values in Y are set to 0 by exrwrite
Negative, NaN, or Inf in RGB or alpha are preserved by exrwrite

