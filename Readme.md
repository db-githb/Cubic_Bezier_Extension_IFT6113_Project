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