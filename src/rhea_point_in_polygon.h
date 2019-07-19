/* RHEA_POINT_IN_POLYGON  Tests wheter a 2D point is inside a polygon. */

#ifndef RHEA_POINT_IN_POLYGON_H
#define RHEA_POINT_IN_POLYGON_H

#include <sc.h>

/**
 * Tests wheter a point with `(testx,testy)` coordinates is inside (or possibly
 * on the boarder of) a polygon with vertices `[(vertx,verty)]_nvert`.
 *
 * \return            True (1) if the test point is inside polygon, false (0)
 *                    otherwise.
 * \param[in] nvert   Number of vertices in the polygon. (Whether to repeat the
 *                    first vertex at the end is discussed below.)
 * \param[in] vertx   Array containing the x-coord's of the polygon's vertices.
 * \param[in] verty   Array containing the y-coord's of the polygon's vertices.
 * \param[in] testx   X-coordinate of the test point.
 * \param[in] testy   Y-coordinate of the test point.
 *
 * THE METHOD:
 * Run a semi-infinite ray horizontally (increasing x, fixed y) out from the
 * test point, and count how many edges it crosses.  At each crossing, the ray
 * switches between inside and outside.  This is according to the Jordan curve
 * theorem.
 *
 * CONCAVE COMPONENTS, MULTIPLE COMPONENTS, AND HOLES:
 * - The polygon may be concave.  However, if a vertex is very close to an edge
 *   (that the vertex is not an end of) then beware of roundoff errors.
 * - The direction that you list the vertices (clockwise or counterclockwise)
 *   does not matter.
 * - The polygon may contain multiple separate components, and/or holes, which
 *   may be concave, provided that you separate the components and holes with a
 *   (0,0) vertex, as follows.
 *     1. First, include a (0,0) vertex.
 *     2. Then include the first component' vertices, repeating its first
 *        vertex after the last vertex.
 *     3. Include another (0,0) vertex.
 *     4. Include another component or hole, repeating its first vertex after
 *        the last vertex.
 *     5. Repeat the above two steps for each component and hole.
 *     6. Include a final (0,0) vertex.
 *   For example, let three components' vertices be
 *     A1, A2, A3, B1, B2, B3, and C1, C2, C3.
 *   Let two holes be
 *     H1, H2, H3, and I1, I2, I3.
 *   Let O be the point (0,0). List the vertices thus:
 *     O, A1, A2, A3, A1, O, B1, B2, B3, B1, O, C1, C2, C3, C1, O, ...
 *     H1, H2, H3, H1, O, I1, I2, I3, I1, O.
 * - Each component or hole's vertices may be listed either clockwise or
 *   counter-clockwise.
 * - If there is only one connected component, then it is optional to repeat
 *   the first vertex at the end. It's also optional to surround the component
 *   with zero vertices.
 *
 ******************************************************************************
 *
 * Algorithm from `pnpoly()` function by Franklin with license as follows.
 *
 * Copyright (c) 1970-2003, Wm. Randolph Franklin
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *   1. Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimers.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice in the documentation and/or other materials provided with the
 *      distribution.
 *   3. The name of W. Randolph Franklin may not be used to endorse or promote
 *      products derived from this Software without specific prior written
 *      permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 * URL: https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
 */
int                 rhea_point_in_polygon_is_inside (
                                              const float testx,
                                              const float testy,
                                              const float *_sc_restrict vertx,
                                              const float *_sc_restrict verty,
                                              const size_t nvert);

#endif /* RHEA_POINT_IN_POLYGON_H */
