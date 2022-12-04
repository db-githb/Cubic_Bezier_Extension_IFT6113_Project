#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

double cubicRoot(double d, double c, double b, double a) {
    double value = (d + 3.0 * c + 3.0 * b + a) / 8.0;
    if (value >= .000001) return cubicRoot(d, (d + c) / 2.0, (d + 2.0 * c + b) / 4.0, value) / 2.0;
    if (value <= -.000001) return 0.5 + cubicRoot(value, (c + 2.0 * b + a) / 4.0, (b + a) / 2.0, a) / 2.0;
    return 0.5;
}

// forces max curvature at Pi
double getReferenceT(std::vector<double> p0, std::vector<double> p1, std::vector<double> p2) {

    std::vector<double> v0 = { p0[0] - p1[0], p0[1] - p1[1] };
    std::vector<double> v2 = { p2[0] - p1[0], p2[1] - p1[1] };

    double d = -(v0[0]*v0[0] + v0[1] * v0[1]); // dot product with itself is the square of the magnitude
    double c = (v0[0] * v2[0] + v0[1] * v2[1]) / 3.0;
    double b = -c;
    double a = v2[0] * v2[0] + v2[1] * v2[1];

    double t = cubicRoot(d, c, b, a);

    return t;
}

double  getPt(double n1, double n2, double t)
{
    int diff = n2 - n1;

    return n1 + (diff * t);
}

std::vector<std::vector<double>> qBezierParams(std::vector<double> p0, std::vector<double> p1, std::vector<double> p2, double t) {
    std::vector<std::vector<double>> bVec(3);

    bVec[0] = p0;
    bVec[2] = p2;

    // calculating b
    // double t = maxCurvature(p0, p1, p2);
    double x2 = (p1[0] - (1.0 - t) * (1.0 - t) * p0[0] - t * t * p2[0]) / (2.0 * (1.0 - t) * t);
    double y2 = (p1[1] - (1.0 - t) * (1.0 - t) * p0[1] - t * t * p2[1]) / (2.0 * (1.0 - t) * t);
    bVec[1] = { x2, y2 };

    return bVec;
}

// from https://stackoverflow.com/questions/785097/how-do-i-implement-a-b%C3%A9zier-curve-in-c
std::vector<double> qBezierInterpolation(std::vector<std::vector<double>> bVec, double t) {
        // The Green Line
       //double xa = getPt(bVec[0][0], bVec[1][0], t);
       //double ya = getPt(bVec[0][1], bVec[1][1], t);
       //double xb = getPt(bVec[1][0], bVec[2][0], t);
       //double yb = getPt(bVec[1][1], bVec[2][1], t);

       // // The Black Dot
       //double x = getPt(xa, xb, t);
       //double y = getPt(ya, yb, t);
       
       double temp = 1 - t;
       double tempSq = temp * temp;

       double x = bVec[0][0] * tempSq + 2 * bVec[1][0] * temp*t + bVec[2][0] * t*t;
       double y = bVec[0][1] * tempSq + 2 * bVec[1][1] * temp*t + bVec[2][1] * t*t;

       return { x,y };
}

std::vector<std::vector<double>> cBezierParams(std::vector<double> p0, std::vector<double> p1, std::vector<double> p2, double t) {

    std::vector<std::vector<double>> bVec(4);

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
    for (int i = 0; i < 2; i++) {
        bVec[1][i] = (p1[i] - tempCb * p0[i] - tcb * p2[i])/(3*temp*t) - t*p1[i] + t*p0[i];
        bVec[2][i] = { bVec[1][i] + (p1[i] - p0[i])};
    }

    // less rounded option
    /*bVec[1][i] = (p1[i] - tempCb * p0[i] - tcb * p2[i])/(3*temp*t) - t*p2[i] + t*p0[i];
     bVec[2][i] = { bVec[1][i] + (p2[i] - p0[i])};*/

    return bVec;
}


std::vector<double> cBezierInterpolation(std::vector<std::vector<double>> bVec, double t) {
    double temp = 1 - t;
    double tempSq = temp * temp;
    double tempCb = tempSq * temp;
    double temptsq = temp * t * t;
    double tcb = t * t * t;
    double tempSqT = tempSq * t;

    double x = bVec[0][0] * tempCb + 3 * bVec[1][0] * tempSqT + 3 * bVec[2][0] * temptsq + bVec[3][0] * tcb;
    double y = bVec[0][1] * tempCb + 3 * bVec[1][1] * tempSqT + 3 * bVec[2][1] * temptsq + bVec[3][1] * tcb;

    return { x,y };

}

std::vector<int> blend(std::vector<std::vector<double>> bVec1, std::vector<std::vector<double>> bVec2, double t, double t1, double t2) {
    
    std::vector<double> p1;
    std::vector<double> p2;

    if (bVec1.size() == 4) {
        p1 = cBezierInterpolation(bVec1, (1 - t1) * t + t1);
        p2 = cBezierInterpolation(bVec2, t * t2);
    }
    else {
        p1 = qBezierInterpolation(bVec1, (1 - t1) * t + t1);
        p2 = qBezierInterpolation(bVec2, t * t2);
    }
    

    double x = cos(M_PI_2 * t) * cos(M_PI_2 * t) * p1[0] + sin(M_PI_2 * t) * sin(M_PI_2 * t) * p2[0];
    double y = cos(M_PI_2 * t) * cos(M_PI_2 * t) * p1[1] + sin(M_PI_2 * t) * sin(M_PI_2 * t) * p2[1];

    return { (int) floor(x), (int) floor(y) };
}

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cout << "Usage: IFT6113_Project.exe curveControlPointsFileName" << std::endl;
        return 1;
    }

    // get curve data from file
    std::string path = "Data/";
    std::string fileName = path.append(argv[1]);

    std::ifstream cp_file(fileName);
    bool isOpen = cp_file.is_open(); 

    int imgHeight;
    int imgWidth;
    int size;

    cp_file >> imgHeight;
    cp_file >> imgWidth;

    // rescale all images to 255 x 255 for faster write to file
    size = (imgHeight > imgWidth) ? imgHeight : imgWidth;
    double scale = 256.0 / size;
    size = 256;

    int numCurves;
    cp_file >> numCurves;

    int option;
    cp_file >> option;

    int numPoints;
    cp_file >> numPoints;

    std::vector<double> point;
    double x;
    double y;
    std::vector<std::vector<double>> ctrlP;

    for(int i = 0; i < numPoints; i++){
        cp_file >> x;
        cp_file >> y;
        point = { x*scale,y*scale };
        ctrlP.push_back(point);
    }

    // std::vector<std::vector<double>> cp = { {42, 42}, { 84,84 }, { 64, 42 }, { 168, 126 }, {84, 235} };
    /*const int height = 256;
    const int width = 256;*/

    // intialize cubic bitmap and quadratic bitmap
	int** cImage = new int*[size];
    int** qImage = new int* [size];

	for (int i = 0; i < size; i++) {
        cImage[i] = new int[size];
        qImage[i] = new int[size];
		for (int j = 0; j < size; j++) {
			cImage[i][j] = 2;
            qImage[i][j] = 2;
		}
	}
    
    for (int i = 1; i < ctrlP.size() - 2; i++) {
        
        // force position of max curvature at middle curve control points (only for quadratic)
        double t1 = getReferenceT(ctrlP[i-1], ctrlP[i], ctrlP[i+1]);
        double t2 = getReferenceT(ctrlP[i], ctrlP[i+1], ctrlP[i+2]); // .431378 seems ideal!!!

        // cubic - compute bezier control points
        std::vector<std::vector<double>> cbVec1 = cBezierParams(ctrlP[i-1], ctrlP[i], ctrlP[i+1], t1);
        std::vector<std::vector<double>> cbVec2 = cBezierParams(ctrlP[i], ctrlP[i+1], ctrlP[i+2], t2);

        // quadratic - compute bezier control points
        std::vector<std::vector<double>> qbVec1 = qBezierParams(ctrlP[i - 1], ctrlP[i], ctrlP[i + 1], t1);
        std::vector<std::vector<double>> qbVec2 = qBezierParams(ctrlP[i], ctrlP[i + 1], ctrlP[i + 2], t2);

        // show bezier control points
       /* for (int i = 0; i < bVec1.size(); i++) {
            int x1 = (int)bVec1[i][0];
            int y1 = (int)bVec1[i][1];
            image[x1][y1] = 0;

            int x2 = (int)bVec2[i][0];
            int y2 = (int)bVec2[i][1];
            image[x2][y2] = 0;
        }*/

        // loop for curve between first two curve control points (not blended)
        if (i == 1) {
            for (double t = 0; t <= t1; t += .001) {
                std::vector<double> cp1 = cBezierInterpolation(cbVec1, t);
                std::vector<double> qp1 = qBezierInterpolation(qbVec1, t);

                cImage[(int)cp1[0]][(int)cp1[1]] = 1;
                qImage[(int)qp1[0]][(int)qp1[1]] = 1;
            }
        }

        // loop for all blended curves (i.e. everything but first and last curve)
        for (double t = 0; t <= 1.0; t += .001) {
            std::vector<int> cPoint = blend(cbVec1, cbVec2, t, t1, t2);
            cImage[cPoint[0]][cPoint[1]] = 1;

            std::vector<int> qPoint = blend(qbVec1, qbVec2, t, t1, t2);
            qImage[qPoint[0]][qPoint[1]] = 1;
        }

        // loop for curve between last two curve control points (not blended)
        if (i == ctrlP.size() - 3) {
            for (double t = t2; t <= 1; t += .001) {
                std::vector<double> cp2 = cBezierInterpolation(cbVec2, t);
                cImage[(int)cp2[0]][(int)cp2[1]] = 1;

                std::vector<double> qp2 = cBezierInterpolation(cbVec2, t);
                qImage[(int)qp2[0]][(int)qp2[1]] = 1;
            }
        }

        
    }

    // color curve control points NOT CHECKING FOR IMAGE EDGE
    for (int i = 0; i < ctrlP.size(); i++) {
        int x = (int) ctrlP[i][0];
        int y = (int) ctrlP[i][1];

        for (int xOff = -1; xOff < 2; xOff++) {
            for (int yOff = -1; yOff < 2; yOff++) {
                cImage[x + xOff][y+ yOff] = 0; // BE CAREFUL HERE
                qImage[x + xOff][y + yOff] = 0; // BE CAREFUL HERE
            }
        }
    }

    //std::ofstream outImg1("control_points.pbm");
    //outImg1 << "P1\n" << size << " " << size << std::endl;
    //for (int i = 0; i < size; i++) {
    //    for (int j = 0; j < size; j++) {
    //        int outpixel = 0;
    //        for (int k = 0; k < cp.size(); k++) {
    //            if (i == cp[k][0] && j == cp[k][1]) {
    //                outpixel = 1;
    //                break;
    //            }
    //        }
    //        outImg1 << outpixel << " ";

    //    }
    //    outImg1 << std::endl;
    //}

    //outImg1.close();

    // color map
    std::vector<std::vector<int>> color_map = { { 0,0,0 }, { 255,0,0 }, {255, 255, 255} };

    std::cout << "saving pgm file";
    std::ofstream cOutImg("interpolation_cb.ppm");
    std::ofstream qOutImg("interpolation_qb.ppm");

    cOutImg << "P3\n" << size << " " << size << std::endl << 255 << std::endl;
    qOutImg << "P3\n" << size << " " << size << std::endl << 255 << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int c = 0; c < 3; c++) {
                cOutImg << color_map[cImage[i][j]][c] << " ";
                qOutImg << color_map[qImage[i][j]][c] << " ";
            }
       }
        cOutImg << std::endl;
        qOutImg << std::endl;
    }

    cOutImg.close();
    qOutImg.close();
}