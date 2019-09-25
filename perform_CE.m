function [normalCE,antiCE] = perform_CE(image,gamma) 
%function [normalCE,antiCE] = perform_CE(image,gamma) 
%This function performs Contrast Enhancement on an input image 'image'
%while making sure the first order (Gray level histogram) and second order
%statistics (GLCM) of the enhanced image resembles the original. This is an
%implementation of the technique described in

%Ravi, H., Subramanyam, A. V., Emmanuel, S., "ACE - An Effective
%Anti-forensic Contrast Enhancement Technique", accepted in IEEE Signal
%Processing Letters (IEEE SPL), 2015. Please cite the paper if you use this
%code.

%Input
% image - any input image that has to be enhanced
% gamma - value of gamma as in Gamma correction. If value is 100, then
% histogram stretching enhancement is performed.
%Output
% normalCE - Normally enhanced image
% antiCE - Image enhanced using our anti forensic method

%This code might not be optimized and completely free of errors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Code written by, @ Hareesh Ravi (Research Associate at IIITD)    %%%%  
%%%%  (haree.24@gmail.com)                                             %%%% 
%%%%  code can be used and modified for research purposes.             %%%% 
%%%%  Kindly let me know of mistakes by mail                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert image to grayscale
if ndims(image) == 3
    img = uint8(rgb2gray(image));
else
    img = uint8(image);
end

%Perform normal contrast enhancement
if gamma == 100
    normalCE = uint8(imadjust(img));
else
    normalCE = uint8(255.*((double(img)./255).^(gamma)));
end

%Perfom anti-forensic contrast enhancement using TV norm optimizaiton
antiCE = uint8(reshape(Anti_TV(img,3,gamma),[size(img,1) size(img,2)]));


end

