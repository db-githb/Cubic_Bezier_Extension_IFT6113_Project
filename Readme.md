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

where InputFileName is the name of a file in the data file.  

## Input File Format
1. Image size: rows x columns
2. 

# Options

Bits | Description
-----|----
0-1 | B1 -> B2 mode: 0: P1 - P0, 1: P2 - P1, 2: P2 - P0, 3: avg(P1-P0,P2-P1,P)
 2 | Use _kappa_-curves _t_ parameter, otherwise use _cy_ quadratic Bezer _t_
 3 | "Loop" effect, otherwise draw Bezier for first and last segments
 4 | Show computed Bezier control points in image

 Examples:
 - 11 = avg + loop options
 - 7 = avg + show Bezier points

 # References: