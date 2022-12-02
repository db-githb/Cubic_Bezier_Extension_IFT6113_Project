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
double maxCurvature(std::vector<double> p0, std::vector<double> p1, std::vector<double> p2) {

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
       double xa = getPt(bVec[0][0], bVec[1][0], t);
       double ya = getPt(bVec[0][1], bVec[1][1], t);
       double xb = getPt(bVec[1][0], bVec[2][0], t);
       double yb = getPt(bVec[1][1], bVec[2][1], t);

        // The Black Dot
       double x = getPt(xa, xb, t);
       double y = getPt(ya, yb, t);

       return { x,y };
}

//std::vector<double> ci(std::vector<std::vector<double>> bVec, double t) {
//        // The Green Lines
//        double xa = getPt(bVec[0][0], bVec[1][0], t);
//        double ya = getPt(bVec[0][1], bVec[1][1], t);
//        double xb = getPt(bVec[1][0], bVec[2][0], t);
//        double yb = getPt(bVec[1][1], bVec[2][1], t);
//        double xc = getPt(bVec[2][0], bVec[3][0], t);
//        double yc = getPt(bVec[2][1], bVec[3][1], t);
//        
//        // The Blue Line
//        double xm = getPt(xa, xb, t);
//        double ym = getPt(ya, yb, t);
//        double xn = getPt(xb, xc, t);
//        double yn = getPt(yb, yc, t);
//
//        // The Black Dot
//        double x = getPt(xm, xn, t);
//        double y = getPt(ym, yn, t);
//
//        return { x,y };
//}

std::vector<int> blend(std::vector<std::vector<double>> bVec1, std::vector<std::vector<double>> bVec2, double t, double t1, double t2) {
    std::vector<double> p1 = qBezierInterpolation(bVec1, (1-t1)*t + t1);
    std::vector<double> p2 = qBezierInterpolation(bVec2, t * t2);

    double x = cos(M_PI_2 * t) * cos(M_PI_2 * t) * p1[0] + sin(M_PI_2 * t) * sin(M_PI_2 * t) * p2[0];
    double y = cos(M_PI_2 * t) * cos(M_PI_2 * t) * p1[1] + sin(M_PI_2 * t) * sin(M_PI_2 * t) * p2[1];

    return { (int) floor(x), (int) floor(y) };
}

std::vector<std::vector<double>> cBezierParams(std::vector<double> p0, std::vector<double> p1, std::vector<double> p2, double t) {
    
    std::vector<std::vector<double>> bVec(4);

    bVec[0] = p0;
    // bVec[2] = p2;
    bVec[3] = p2;

    double temp = 1 - t;
    double tempSq = temp * temp;
    double tempCb = tempSq * temp;
    double temptsq = temp * t * t;
    double tcb = t * t * t;
    double tempSqT = tempSq * t;
    
    /*bVec[1] = { (p1[0] - p0[0] * tempCb - 3 * p2[0] * temptsq - p3[0] * tcb) / (3 * tempSqT),
              (p1[1] - p0[1] * tempCb - 3 * p2[1] * temptsq - p3[1] * tcb) / (3 * tempSqT) };*/

    bVec[1] = { (1.0 / 2.0) * ((5.0/3.0) * p1[0] + (2.0/3.0)*p0[0] - (1.0/3.0)* p2[0]), (1.0 / 2.0) * ((5.0 / 3.0) * p1[1] + (2.0 / 3.0) * p0[1] - (1.0 / 3.0) * p2[1]) };
    bVec[2] = { bVec[1][0] + (p1[0]-p0[0]), bVec[1][1] + (p1[1] - p0[1]) };

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

int main() {
    std::vector<std::vector<double>> cp = { {42, 42}, { 84,84 }, { 126, 42 }, { 168, 126 } };
	const int size = 256;
	int** image = new int*[size];
	for (int i = 0; i < size; i++) {
        image[i] = new int[size];
		for (int j = 0; j < size; j++) {
			image[i][j] = 2;
		}
	}
    double t1 = .5; // maxCurvature(cp[0], cp[1], cp[2]);
    double t2 = .5; //maxCurvature(cp[1], cp[2], cp[3]); // .431378 seems ideal!!!
    std::vector<std::vector<double>> bVec1 = cBezierParams(cp[0], cp[1], cp[2], t1);
    std::vector<std::vector<double>> bVec2 = cBezierParams(cp[1], cp[2], cp[3], t2);
    
    for (double t = 0; t <= 1; t += .001) {
        std::vector<double> p1 = cBezierInterpolation(bVec1, t);
        image[(int)p1[0]][(int)p1[1]] = 1;
    }

    /*for (double t = 0; t <=1.0; t += .001) {
         std::vector<int> point = blend(bVec1, bVec2, t, t1, t2);
         image[point[0]][point[1]] = 1;
    }*/

    for (double t = 0; t <= 1; t += .001) {
        std::vector<double> p2 = cBezierInterpolation(bVec2, t);
        image[(int)p2[0]][(int)p2[1]] = 1;
    }

    // color curve control points
    for (int i = 0; i < cp.size(); i++) {
        int x = (int) cp[i][0];
        int y = (int) cp[i][1];
        image[x][y] = 0;
    }

    // bezier control points
    for (int i = 0; i < bVec1.size(); i++) {
        int x1 = (int)bVec1[i][0];
        int y1 = (int)bVec1[i][1];
        image[x1][y1] = 0;

        int x2 = (int)bVec2[i][0];
        int y2 = (int)bVec2[i][1];
        image[x2][y2] = 0;
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

    std::ofstream outImg("cb_interpolation.pgm");
    outImg << "P2\n" << size << " " << size << std::endl << 2 << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            outImg << image[i][j]<< " ";
       }
        outImg << std::endl;
    }

    outImg.close();
}