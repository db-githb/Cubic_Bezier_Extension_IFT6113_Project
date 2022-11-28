#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

double  getPt(double n1, double n2, double perc)
{
    int diff = n2 - n1;

    return n1 + (diff * perc);
}

std::vector<std::vector<double>> qBezierParams(std::vector<double> p0, std::vector<double> p1, std::vector<double> p2) {
    std::vector<std::vector<double>> bVec(3);

    bVec[0] = p0;
    bVec[2] = p2;

    // calculating b
    double t = .5;
    double x2 = (p1[0] - (1 - t) * (1 - t) * p0[0] - t * t * p2[0]) / (2 * (1 - t) * t);
    double y2 = (p1[1] - (1 - t) * (1 - t) * p0[1] - t * t * p2[1]) / (2 * (1 - t) * t);
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

std::vector<int> blend(std::vector<std::vector<double>> bVec1, std::vector<std::vector<double>> bVec2, double t) {
    std::vector<double> p1 = qBezierInterpolation(bVec1, t);
    std::vector<double> p2 = qBezierInterpolation(bVec2, t);

    double x = cos(M_PI_2 * t) * cos(M_PI_2 * t) * p1[0] + sin(M_PI_2 * t) * sin(M_PI_2 * t) * p2[0];
    double y = cos(M_PI_2 * t) * cos(M_PI_2 * t) * p1[1] + sin(M_PI_2 * t) * sin(M_PI_2 * t) * p2[1];

    return { (int) floor(x), (int) floor(y) };
}

int main() {
    std::vector<std::vector<double>> cp = { {42, 42}, { 84,84 }, { 126, 42 }, { 168, 126 }};
	const int size = 256;
	int** image = new int*[size];
	for (int i = 0; i < size; i++) {
        image[i] = new int[size];
		for (int j = 0; j < size; j++) {
			image[i][j] = 2;
		}
	}
    std::vector<std::vector<double>> bVec1 = qBezierParams(cp[0], cp[1], cp[2]);
    std::vector<std::vector<double>> bVec2 = qBezierParams(cp[1], cp[2], cp[3]);

    for (double t = 0; t <= 1; t += .001) {
        //std::vector<int> point = blend(bVec1, bVec2, t);
        std::vector<double> p1 = qBezierInterpolation(bVec1, t);
        std::vector<double> p2 = qBezierInterpolation(bVec2, t);

        image[(int) p1[0]][(int) p1[1]] = 1;
        image[(int) p2[0]][(int) p2[1]] = 1;

    }

    // color control points
    for (int i = 0; i < cp.size(); i++) {
        int x = (int) cp[i][0];
        int y = (int) cp[i][1];
        image[x][y] = 0;
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

    std::ofstream outImg("qb_interpolation.pgm");
    outImg << "P2\n" << size << " " << size << std::endl << 2 << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            outImg << image[i][j]<< " ";
       }
        outImg << std::endl;
    }

    outImg.close();
}