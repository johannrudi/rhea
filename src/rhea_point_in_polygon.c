/*
 */

#include <rhea_point_in_polygon.h>

int
rhea_point_in_polygon_is_inside (const float testx,
                                 const float testy,
                                 const float *_sc_restrict vertx,
                                 const float *_sc_restrict verty,
                                 const size_t nvert)
{
  size_t              i, j;
  int                 c = 0;

  for (i = 0, j = nvert-1; i < nvert; j = i++) { /* loop over all vertices */
    /* check if ray intersects an edge of the polygon */
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
         (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) ) {
       /* flip even/odd #intersections flag */
       c = !c;
    }
  }

  /* return even/odd #intersections flag (odd means pt. is inside polygon) */
  return c;
}
