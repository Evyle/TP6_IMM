/*!
  \file usImageGradient.cpp
  \brief Process images
*/

#include "usImageGradient.h"

#define TO_COMPLETE

/*!
  \brief Compute the gradient images along x, y and z directions from a thin volume composed of 3 2D slices
  \param imga, img0, imgb: parallel images captured
  \param ip0,ipf: upper left and lower right corners of the ROI (in pixels) inside the gradient is computed
  \param dIdx,dIdy,dIdz: gradient images along x,y and z directions to compute
*/
void usImageGradient::Grad3DF3x3x3 (vpImage<unsigned char> &imga, vpImage<unsigned char> &img0, vpImage<unsigned char> &imgb, vpImage<double> &dIdx, vpImage<double> &dIdy, vpImage<double> &dIdz, vpImagePoint ip0, vpImagePoint ipf )
{

    int Wmin = (int)ip0.get_u(); // pixel u coordinate of the ROI left-top corner
    int Wmax = (int)ipf.get_u(); // pixel u coordinate of the ROI right-bottom corner
    int Hmin = (int)ip0.get_v(); // pixel v coordinate of the ROI left-top corner
    int Hmax = (int)ipf.get_v(); // pixel v coordinate of the ROI right-bottom corner

#ifdef TO_COMPLETE
    // Question 3
    // Compute the gradient images along x, y and z directions from a thin volume composed of three 2D slices
    // corresponding to images imga, img0, imgb
    // Note that dIdx, dIdy and dIdz have the same size that the whole 2D ultrasound image but the gradient could
    // be computed only for the pixels contained in the ROI

    vpMatrix patchX(3,3);
    patchX[0][0] = -1;  patchX[1][0] = -2;  patchX[2][0] = -1;
    patchX[0][1] = 0;   patchX[1][1] = 0;   patchX[2][1] = 0;
    patchX[0][2] = 1;   patchX[1][2] = 2;   patchX[2][2] = 1;

    vpMatrix patchY = patchX.transpose();

    vpMatrix patchZ(3,3);
    patchZ[0][0] = 1;  patchZ[1][0] = 2;  patchZ[2][0] = 1;
    patchZ[0][1] = 2;   patchZ[1][1] = 4;   patchZ[2][1] = 2;
    patchZ[0][2] = 1;   patchZ[1][2] = 2;   patchZ[2][2] = 1;

    for(int i=Hmin; i<Hmax; i++) {
        for(int j=Wmin; j<Wmax; j++) {
            double img0Fx = ApplyPatch(2*patchX, img0, i, j);
            double img0Fy = ApplyPatch(2*patchY, img0, i, j);
            double img0Fz = 0;

            double imgaFx = ApplyPatch(patchX, imga, i, j);
            double imgaFy = ApplyPatch(patchY, imga, i, j);
            double imgaFz = ApplyPatch(-patchZ, imga, i, j);

            double imgbFx = ApplyPatch(patchX, imgb, i, j);
            double imgbFy = ApplyPatch(patchY, imgb, i, j);
            double imgbFz = ApplyPatch(patchZ, imgb, i, j);

            dIdx[i][j] = img0Fx + imgaFx +  imgbFx;
            dIdy[i][j] = img0Fy + imgaFy +  imgbFy;
            dIdz[i][j] = img0Fz + imgaFz +  imgbFz;
        }
    }


#endif

}

/**
 * @brief usImageGradient::ApplyPatch This method apply a patch to a point (i, j) in a image img
 * @param patch : the patch to apply (ex : a 3x3 image)
 * @param img : the original image
 * @param i : coordonate i of the point
 * @param j : coordonate j of the point
 * @return the result of the patch applying to this point
 */
double usImageGradient::ApplyPatch(const vpMatrix &patch, vpImage<unsigned char> &img, int & i, int & j)
{
    double result = 0;

    for(int iPatch = -(int)(patch.getRows()/2); iPatch <= (int)patch.getRows()/2; iPatch++) {
        for(int jPatch = -(int)(patch.getCols()/2); jPatch <= (int)patch.getCols()/2; jPatch++) {
            result += patch[iPatch+patch.getRows()/2][jPatch+patch.getCols()/2] * img[i+iPatch][j+jPatch];
        }
    }
    return result;
}


