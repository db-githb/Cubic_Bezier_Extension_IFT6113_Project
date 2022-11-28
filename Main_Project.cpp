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

// from https://stackoverflow.com/questions/785097/how-do-i-implement-a-b%C3%A9zier-curve-in-c
void qBezierInterpolation(int** image,  std::vector<double> p0, std::vector<double> p1, std::vector<double> p2) {
    for (double i = 0; i < 1; i += 0.001)
    {
        double x1 = p0[0];
        double y1 = p0[1];
        double x2 = p1[0];
        double y2 = p1[1];
        double x3 = p2[0];
        double y3 = p2[1];

        // The Green Line
       double xa = getPt(x1, x2, i);
       double ya = getPt(y1, y2, i);
       double xb = getPt(x2, x3, i);
       double yb = getPt(y2, y3, i);

        // The Black Dot
       int x = floor(getPt(xa, xb, i));
       int y = floor(getPt(ya, yb, i));

       image[x][y] = 1;
    }
}


int main() {
    std::vector<std::vector<double>> cp = { {42, 42}, { 84,84 }, { 126, 42 }, { 168, 126 }};
	const int size = 256;
	int** image = new int*[size];
	for (int i = 0; i < size; i++) {
        image[i] = new int[size];
		for (int j = 0; j < size; j++) {
			image[i][j] = 0;
		}
	}
    
    qBezierInterpolation(image, cp[0], cp[1], cp[2]);

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

    std::ofstream outImg("qb_interpolation.pbm");
    outImg << "P1\n" << size << " " << size << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            outImg << image[i][j]<< " ";
       }
        outImg << std::endl;
    }

    outImg.close();
}