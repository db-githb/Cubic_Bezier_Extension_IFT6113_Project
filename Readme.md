# IFT 6113 - Final Project

An implementation of Cem Yuksel's paper, "A Class of C2 Interpolating Spline" (2020)  
with an extension to use Cubic Bezier as an interpolating function  
By Damian Bowness (20257155)  
for Prof. Mikhail Bessmeltsev  
U de M, Fall 2022  

# Contents
- Data:  
  - input data files  
- Debug:  
  - executable: IFT6113_Project.exe
- Output:
  - generated image files
- Presentation
  - Powerpoint slides
  - presentation notes    
- Main_Project.cpp:
  - source code

# Usage
Application runs on Windows.  

From main folder (IFT6113_Project) run executable as follows:  

> Debug\IFT6113_Project.exe InputFileName

where InputFileName is the name of a file in the Data folder which has a .txt extension.  
This generates files InputFileName_cb.ppm and InputFileName_qb.ppm in the Output folder.  
These are image files in the portable pixmap format.  
_cb and _qb are the blended interpolation using cubic Bezier and quadratic Bezier respectively.

## Input File Format
1. Image size: rows x columns
2. Option: [see details below](#options)
3. Number of Points
4. Curve Control Points: Column, Row (x,y) for each curve control point
    1. Note: values must be within image size.  Not validated by code. 

# Options
Options number is a bit encoded set of options  

Bits | Description
-----|----
0-1 | B1 -> B2 mode: 0: P1 - P0, 1: P2 - P1, 2: P2 - P0, 3: avg(P1-P0,P2-P1,P)
 2 | Use _kappa_-curves _t_ parameter, otherwise use _cy_ quadratic Bezer _t_
 3 | "Loop" effect, otherwise draw Bezier for first and last segments
 4 | Show computed Bezier control points in image

 Examples:
 - 11 = avg + loop options
 - 7 = avg + show Bezier points

 # References
 1. C. Yuksel, “A class of *C*<sup>2</sup> interpolating splines,” ACM Transactions on Graphics, vol. 39, no. 5, pp. 1–14, 2020.  
 2. D. J. Walton and D. S. Meek, “Curvature extrema of planar parametric polynomial cubic curves,” Journal of Computational and Applied Mathematics, vol. 134, no. 1-2, pp. 69–83, 2001.   
 3. Z. Yan, S. Schiller, G. Wilensky, N. Carr, and S. Schaefer, “&kappa;-curves,” ACM Transactions on Graphics, vol. 36, no. 4, pp. 1–7, 2017. 
