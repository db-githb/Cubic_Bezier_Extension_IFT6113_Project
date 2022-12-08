/*
* An implementation of Cem Yuksel's paper, "A Class of C2 Interpolating Spline" (2020)
* with an extension to use Cubic Bezier as an interpolating function
* by Damian Bowness (20257155)
* for Prof. Mikhail Bessmeltsev (U de M)
* IFT 6113, Final Project, Fall 2022
*/


#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void showBezierPoints(int** image, vector<vector<double>> bVec)
{
  // show bezier control points
  for (int i = 0; i < bVec.size(); i++) {
      int row = (int)bVec[i][0];
      int col = (int)bVec[i][1];
      if (row < 0) row = 0;
      else if (row > 511) row = 511;
      if (col < 0) row = 0;
      else if (col > 511) row = 511;
      image[row][col] = 3;
  }
}

// debug function
void printBvec(vector<vector<double>> bVec) {
  for (int i = 0; i < bVec.size(); i += 1) {
    std::cout << "(" << bVec[i][0] << ", " << bVec[i][1] << ")  ";
  }
  std::cout << endl;
}

// read in curve control points from file
vector<vector<double>> readCtrlPoints(bool loopOption, double scale, ifstream &cp_file) {
  vector<vector<double>> ctrlP;
  int numPoints;
  cp_file >> numPoints;
  double col;
  double row;
  for (int i = 0; i < numPoints; i++) {
    cp_file >> col;
    cp_file >> row;
    ctrlP.push_back({ row * scale, col * scale });
  }
  if (loopOption) {
    // Need to blend first, last, and gap curves
    ctrlP.push_back(ctrlP[0]);
    ctrlP.push_back(ctrlP[1]);
    ctrlP.push_back(ctrlP[2]);
  }
  return ctrlP;
}

// need cubic curve to find t that forces max curvature on p1
double cubicRoot(double d, double c, double b, double a) {
    double value = (d + 3.0 * c + 3.0 * b + a) / 8.0;
    if (value >= .000001) return cubicRoot(d, (d + c) / 2.0, (d + 2.0 * c + b) / 4.0, value) / 2.0;
    if (value <= -.000001) return 0.5 + cubicRoot(value, (c + 2.0 * b + a) / 4.0, (b + a) / 2.0, a) / 2.0;
    return 0.5;
}

// application supports option to use Zhipei Yan's method for finding t that forces max curvature to be on p1
// https://github.com/zhipeiyan/kappa-Curves
double kCurvesT(vector<double> p0, vector<double> p1, vector<double> p2)
{
  vector<double> p2Mp0 = { p2[0] - p0[0], p2[1] - p0[1] };
  vector<double> p0Mp1 = { p0[0] - p1[0], p0[1] - p1[1] };

  double a = p2Mp0[0] * p2Mp0[0] + p2Mp0[1] * p2Mp0[1]; //(p2 - p0).dot(p2 - p0);
  double b = 3 * (p2Mp0[0] * p0Mp1[0] + p2Mp0[1] * p0Mp1[1]); //3 * (p2 - p0).dot(p0 - p1);

  double c = 0; // = (3 * p0 - 2 * p1 - p2).dot(p0 - p1);
  for (int i = 0; i < 2; i += 1) {
    c += (3 * p0[i] - 2 * p1[i] - p2[i]) * p0Mp1[i];
  }

  double d = -(p0Mp1[0] * p0Mp1[0] + p0Mp1[1] * p0Mp1[1]);// -(p0 - p1).dot(p0 - p1);

  double p = (3 * a * c - b * b) / 3 / a / a;
  double q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / 27 / a / a / a;
  // solve u^3 + p u + q = 0
  if (4 * p * p * p + 27 * q * q >= 0)
  {
    // single real root
    return cbrt(-q / 2 + sqrt(q * q / 4 + p * p * p / 27)) + cbrt(-q / 2 - sqrt(q * q / 4 + p * p * p / 27)) - b / 3 / a;
  }
  else
  {
    // three real roots
    for (int k = 0; k < 3; ++k)
    {
      double t = 2 * sqrt(-p / 3) * cos(1. / 3 * acos(3 * q / 2. / p * sqrt(-3 / p)) - 2 * M_PI * k / 3.) - b / 3 / a;
      if (0 <= t && t <= 1)
        return t;
    }
  }
  // ignore Yan's error and use default t=0.5
  return 0.5;
}

// use Yuksel's method for finding t that forces max curvature to be on p1 for quadratic Bezier 
// (implemented from paper)
double qBezierT(vector<double> p0, vector<double> p1, vector<double> p2) {

    vector<double> v0 = { p0[0] - p1[0], p0[1] - p1[1] };
    vector<double> v2 = { p2[0] - p1[0], p2[1] - p1[1] };

    double d = -(v0[0]*v0[0] + v0[1] * v0[1]); // dot product with itself is the square of the magnitude
    double c = (v0[0] * v2[0] + v0[1] * v2[1]) / 3.0;
    double b = -c;
    double a = v2[0] * v2[0] + v2[1] * v2[1];

    double t = cubicRoot(d, c, b, a);

    return t;
}

// general method to select between Yan's k-curves or Yuksel's quadratic bezier
double getReferenceT(bool kCurves, vector<double> p0, vector<double> p1, vector<double> p2) {
  return kCurves ? kCurvesT(p0, p1, p2) : qBezierT(p0, p1, p2);
}

// compute Bezier control points for quadratic Bezier
vector<vector<double>> qBezierParams(vector<double> p0, vector<double> p1, vector<double> p2, double t) {
    
    vector<vector<double>> bVec(3);
    
    // first and last bezier control points bound to first and last curve control points
    bVec[0] = p0;
    bVec[2] = p2;

    // calculating b1: equation 6 from paper
    double x2 = (p1[0] - (1.0 - t) * (1.0 - t) * p0[0] - t * t * p2[0]) / (2.0 * (1.0 - t) * t);
    double y2 = (p1[1] - (1.0 - t) * (1.0 - t) * p0[1] - t * t * p2[1]) / (2.0 * (1.0 - t) * t);
    bVec[1] = { x2, y2 };

    return bVec;
}

// Interpolate curve given quadratic Bezier control points at position t
vector<double> qBezierInterpolation(vector<vector<double>> bVec, double t) {
       
       double temp = 1 - t;
       double tempSq = temp * temp;

       // Implemented from https://en.wikipedia.org/wiki/B%C3%A9zier_curve (or rearrange equation 6...)
       double x = bVec[0][0] * tempSq + 2 * bVec[1][0] * temp*t + bVec[2][0] * t*t;
       double y = bVec[0][1] * tempSq + 2 * bVec[1][1] * temp*t + bVec[2][1] * t*t;

       return { x,y };
}

// Implementation of one part of extension to paper
// Define vector that connects B1 to B2 based on relationship between curve control points
vector<double> getB2Offset(int option, vector<double> p0, vector<double> p1, vector<double> p2) {

    switch (option){
    default:
        // B2-B1 = p1-p0
        return { p1[0] - p0[0], p1[1] - p0[1] };

    case 1:
        // B2-B1 = p2-p0
        return { p2[0] - p1[0], p2[1] - p1[1] };

    case 2:
        // B2-B1 = p2-p1
        return { p2[0] - p0[0], p2[1] - p0[1] };

    case 3:
        // B2-B1 = avg(p1-p0, p2-p0, p2-p0) :: Produced best looking results
        vector<double> retVec = { 0,0 };
        for (int i = 0; i < 2; i++) {
            retVec[i] = (p1[i] - p0[i] + p2[i] - p1[i] + p2[i] - p0[i]) / 3.0;
        }
        return retVec;
    }
}

// Implementation of one part of extension to paper
// compute Bezier control points for quadratic Bezier
vector<vector<double>> cBezierParams(int option, vector<double> p0, vector<double> p1, vector<double> p2, double t) {

    vector<vector<double>> bVec(4);

    bVec[0] = p0;
    bVec[3] = p2;

    double temp = 1 - t;
    double tempSq = temp * temp;
    double tempCb = tempSq * temp;
    double temptsq = temp * t * t;
    double tcb = t * t * t;
    double tempSqT = tempSq * t;

    bVec[1] = { 0,0 };
    bVec[2] = { 0,0 };

    // Alternative relationships between B1 and B2
    if (option == 4 || (option == 6 && t >= 0.5)) {

      // define B1 and B2 by aligning B2 vertically above p2
      bVec[2] = { p2[0], 2. * (p2[1] - p1[1]) };
      for (int i = 0; i < 2; i++) {
        bVec[1][i] = (p1[i] - tempCb * p0[i] - tcb * p2[i]) / (3 * tempSq * t) - bVec[2][i] * (t / temp);
      }
    }
    else if (option == 5 || option == 6) {

        // define B1 and B2 by aligning B1 vertically above p0
        bVec[1] = { p0[0], 2. * (p1[1] - p0[1]) };
        for (int i = 0; i < 2; i++) {
          bVec[2][i] = (p1[i] - tempCb * p0[i] - tcb * p2[i]) / (3 * temptsq) - bVec[1][i] * (temp / t);
        }
    }
    // Alternative relationship produce horrible results so define B1 to B2 based on relationship between curve control points
    else {
      vector<double> b2_offset = getB2Offset(option, p0, p1, p2);
      for (int i = 0; i < 2; i++) {
        bVec[1][i] = (p1[i] - tempCb * p0[i] - tcb * p2[i]) / (3 * temp * t) - t * b2_offset[i];
        bVec[2][i] = { bVec[1][i] + b2_offset[i] };
      }
    }

    return bVec;
}

// Interpolate curve given cubic Bezier control points at position t
vector<double> cBezierInterpolation(vector<vector<double>> bVec, double t) {
    double temp = 1 - t;
    double tempSq = temp * temp;
    double tempCb = tempSq * temp;
    double temptsq = temp * t * t;
    double tcb = t * t * t;
    double tempSqT = tempSq * t;

    // Implemented from // https://en.wikipedia.org/wiki/B%C3%A9zier_curve
    double x = bVec[0][0] * tempCb + 3 * bVec[1][0] * tempSqT + 3 * bVec[2][0] * temptsq + bVec[3][0] * tcb;
    double y = bVec[0][1] * tempCb + 3 * bVec[1][1] * tempSqT + 3 * bVec[2][1] * temptsq + bVec[3][1] * tcb;

    return { x,y };

}

// Implementation of Yuksel's blend function from paper
vector<int> blend(vector<vector<double>> bVec1, vector<vector<double>> bVec2, double t, double t1, double t2) {
    
    vector<double> pt_from_C1;
    vector<double> pt_from_C2;

    if (bVec1.size() == 4) {
        pt_from_C1 = cBezierInterpolation(bVec1, (1 - t1) * t + t1);
        pt_from_C2 = cBezierInterpolation(bVec2, t * t2);
    }
    else {
        pt_from_C1 = qBezierInterpolation(bVec1, (1 - t1) * t + t1);
        pt_from_C2 = qBezierInterpolation(bVec2, t * t2);
    }

    // now blend point from curve 1 with point from curve 2
    double x = cos(M_PI_2 * t) * cos(M_PI_2 * t) * pt_from_C1[0] + sin(M_PI_2 * t) * sin(M_PI_2 * t) * pt_from_C2[0];
    double y = cos(M_PI_2 * t) * cos(M_PI_2 * t) * pt_from_C1[1] + sin(M_PI_2 * t) * sin(M_PI_2 * t) * pt_from_C2[1];

    return { (int) floor(x), (int) floor(y) };
}

// write point to bitmap for first and last Bezier curve;s (i.e. no loop option selected)
void drawBezierCurve(int** qImage, vector<vector<double>> qbVec, int** cImage, vector<vector<double>> cbVec, double startT, double stopT) {
  if (startT >= stopT) {
    return;
  }
  for (double t = startT; t <= stopT; t += .001) {
    vector<double> cp1 = cBezierInterpolation(cbVec, t);
    cImage[(int)cp1[0]][(int)cp1[1]] = 4;

    vector<double> qp1 = qBezierInterpolation(qbVec, t);
    qImage[(int)qp1[0]][(int)qp1[1]] = 4;
  }
}

// ==========================
// || MAIN DRIVER FUNCTION || 
// ==========================
int main(int argc, char* argv[])
{
  if (argc != 2) {
    std::cout << "Usage: IFT6113_Project.exe curveControlPointsFileName" << endl;
    return 1;
  }

  std::cout << "Reading...";

  // get curve data from file
  string path = "Data/";
  string fileName = path.append(argv[1]);

  ifstream cp_file(fileName);
  bool isOpen = cp_file.is_open();

  int imgHeight;
  int imgWidth;
  cp_file >> imgHeight;
  cp_file >> imgWidth;

  // rescale all images to 255 x 255 for faster write to file
  int size = (imgHeight > imgWidth) ? imgHeight : imgWidth;
  double scale = 512.0 / size;
  size = 512;

  int numCurves;
  cp_file >> numCurves;

  int option;
  cp_file >> option;

  // bitmap to decode application options
  int b2OffsetOption = option & 3;
  bool kCurvesToption = option & 4;
  bool loopOption = option & 8;
  bool showBvecOption = option & 16;

  vector<vector<double>> ctrlP = readCtrlPoints(loopOption, scale, cp_file);

  std::cout << endl << "Processing...";

  // intialize cubic bitmap and quadratic bitmap
  // application outputs two image files
  int** cImage = new int* [size];
  int** qImage = new int* [size];

  for (int i = 0; i < size; i++) {
    cImage[i] = new int[size];
    qImage[i] = new int[size];
    for (int j = 0; j < size; j++) {
      cImage[i][j] = 2;
      qImage[i][j] = 2;
    }
  }

  // Initialize bezier control points. Each iteration needs 2 sets of control
  // points, but 1st set is the 2nd set of previous iteration
  double t1 = 0.;
  vector<vector<double>> cbVec1 = { {0,0},{0,0},{0,0},{0,0} };
  double t2 = getReferenceT(kCurvesToption, ctrlP[0], ctrlP[1], ctrlP[2]);
  vector<vector<double>> cbVec2 = cBezierParams(b2OffsetOption, ctrlP[0], ctrlP[1], ctrlP[2], t2);

  vector<vector<double>> qbVec1 = { {0,0},{0,0},{0,0} };
  vector<vector<double>> qbVec2 = qBezierParams(ctrlP[0], ctrlP[1], ctrlP[2], t2);

  for (size_t i = 1; i < ctrlP.size() - 2; i++) {
    
     // t1 of current curve is t2 of previous curve
    t1 = t2;
    t2 = getReferenceT(kCurvesToption, ctrlP[i], ctrlP[i + 1], ctrlP[i + 2]);

    // cubic - compute bezier control points
    cbVec1 = cbVec2;
    cbVec2 = cBezierParams(b2OffsetOption, ctrlP[i], ctrlP[i + 1], ctrlP[i + 2], t2);

    // quadratic - compute bezier control points
    qbVec1 = qbVec2;
    qbVec2 = qBezierParams(ctrlP[i], ctrlP[i + 1], ctrlP[i + 2], t2);

    // curve between first two curve control points (not blended)
    if (i == 1 && !loopOption) {
      drawBezierCurve(qImage, qbVec1, cImage, cbVec1, 0, t1);
    }

    // loop for all blended curves (i.e. everything but first and last curve)
    for (double t = 0; t <= 1.0; t += .001) {
      vector<int> cPoint = blend(cbVec1, cbVec2, t, t1, t2);
      cImage[cPoint[0]][cPoint[1]] = 1;

      vector<int> qPoint = blend(qbVec1, qbVec2, t, t1, t2);
      qImage[qPoint[0]][qPoint[1]] = 1;
    }

    // curve between last two curve control points (not blended)
    if (i == ctrlP.size() - 3 && !loopOption) {
      drawBezierCurve(qImage, qbVec2, cImage, cbVec2, t2, 1);
    }
    if (showBvecOption) {
      showBezierPoints(cImage, cbVec1);
      showBezierPoints(qImage, qbVec1);
    }
  }
  if (showBvecOption) {
    showBezierPoints(cImage, cbVec2);
    showBezierPoints(qImage, qbVec2);
  }

  // color curve control points NOT CHECKING FOR IMAGE EDGE
  for (size_t i = 0; i < ctrlP.size(); i++) {
    int x = (int)ctrlP[i][0];
    int y = (int)ctrlP[i][1];

    for (int xOff = -1; xOff < 2; xOff++) {
      for (int yOff = -1; yOff < 2; yOff++) {
        cImage[x + xOff][y + yOff] = 0; // BE CAREFUL HERE
        qImage[x + xOff][y + yOff] = 0; // BE CAREFUL HERE
      }
    }
  }

  // color map
  vector<vector<int>> color_map = { { 0,0,0 }, { 255,0,0 }, {255, 255, 255}, {127, 127, 127}, {0, 0, 255} };

  std::cout << endl << "Saving...";
  ofstream cOutImg("interpolation_cb.ppm");
  ofstream qOutImg("interpolation_qb.ppm");

  cOutImg << "P3\n" << size << " " << size << endl << 255 << endl;
  qOutImg << "P3\n" << size << " " << size << endl << 255 << endl;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      for (int c = 0; c < 3; c++) {
        cOutImg << color_map[cImage[i][j]][c] << " ";
        qOutImg << color_map[qImage[i][j]][c] << " ";
      }
    }
    cOutImg << endl;
    qOutImg << endl;
  }

  cOutImg.close();
  qOutImg.close();
  std::cout << "Done" << endl;
}