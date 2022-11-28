#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

// function from https://www.geeksforgeeks.org/equation-of-circle-when-three-points-on-the-circle-are-given/
std::vector<double> getCircleParams(std::vector<double> p0, std::vector<double> p1, std::vector<double> p2)
{

    int x1 = p0[0];
    int y1 = p0[1];
    int x2 = p1[0];
    int y2 = p1[1];
    int x3 = p2[0];
    int y3 = p2[1];

    int x12 = x1 - x2;
    int x13 = x1 - x3;

    int y12 = y1 - y2;
    int y13 = y1 - y3;

    int y31 = y3 - y1;
    int y21 = y2 - y1;

    int x31 = x3 - x1;
    int x21 = x2 - x1;

    // x1^2 - x3^2
    int sx13 = pow(x1, 2) - pow(x3, 2);

    // y1^2 - y3^2
    int sy13 = pow(y1, 2) - pow(y3, 2);

    int sx21 = pow(x2, 2) - pow(x1, 2);
    int sy21 = pow(y2, 2) - pow(y1, 2);

    int f = ((sx13) * (x12)
        +(sy13) * (x12)
        +(sx21) * (x13)
        +(sy21) * (x13))
        / (2 * ((y31) * (x12)-(y21) * (x13)));
    int g = ((sx13) * (y12)
        +(sy13) * (y12)
        +(sx21) * (y13)
        +(sy21) * (y13))
        / (2 * ((x31) * (y12)-(x21) * (y13)));

    int c = -pow(x1, 2) - pow(y1, 2) - 2 * g * x1 - 2 * f * y1;

    // eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0
    // where centre is (h = -g, k = -f) and radius r
    // as r^2 = h^2 + k^2 - c
    int h = -g;
    int k = -f;
    int sqr_of_r = h * h + k * k - c;

    // r is the radius
    double  r = sqrt(sqr_of_r);

    std::vector<double> returnVec(3);
    returnVec[0] = h;
    returnVec[1] = k;
    returnVec[2] = r;
    return returnVec;
}

std::vector<double> circleInterpolation(double theta, std::vector<double> params) {
    double h = params[0];
    double k = params[1];
    double r = params[2];

    double x = r * cos(theta) + h; // add h to offset because circle origin is not (0,0)
    double y = r * sin(theta) + k; // add h to offset because circle origin is not (0,0)

    std::vector<double> point = { x,y }; // interpolated point on circumference

    return point;
}

int main() {
    std::vector<std::vector<double>> cp = { {42, 42}, { 84,84 }, { 126, 42 }, { 168, 126 }};
	const int size = 256;
	int image[size][size];
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			image[i][j] = 0;
		}
	}

	for (int i = 1; i < cp.size() - 2; i++) {
		std::vector<double> p0 = cp[i - 1];
		std::vector<double> p1 = cp[i];
		std::vector<double> p2 = cp[i + 1];
		std::vector<double> p3 = cp[i + 2];

        // getParams returns vector <h,k,r>
        std::vector<double> f1_params = getCircleParams(p0, p1, p2);
        std::vector<double> f2_params = getCircleParams(p1, p2, p3);
        double theta_p1 = atan2(p1[1] - f1_params[1], p1[0] - f1_params[0]);
        double theta_p2 = atan2(p2[1] - f2_params[1], p2[0] - f2_params[0]);

        for (double theta = 0; theta <= M_PI_2; theta += M_PI_2 / 1000.0) {
            double cos2 = cos(theta) * cos(theta);
            double sin2 = sin(theta) * sin(theta);

            std::vector<double> f1 = circleInterpolation(theta+theta_p1, f1_params);
            std::vector<double> f2 = circleInterpolation(theta+theta_p2, f2_params);

            int cx = floor(cos2 * f1[0] + sin2 * f2[0]);
            int cy = floor(cos2 * f1[1] + sin2 * f2[1]);

            image[cx][cy] = 1;
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

    std::ofstream outImg("interpolation.pbm");
    outImg << "P1\n" << size << " " << size << std::endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            outImg << image[i][j]<< " ";
       }
        outImg << std::endl;
    }

    outImg.close();
}