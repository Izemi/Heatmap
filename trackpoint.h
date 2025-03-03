#ifndef __TRACKPOINT_H__
#define __TRACKPOINT_H__

#include "location.h"

typedef struct trackpoint trackpoint;

/**
 * Creates a trackpoint with the given values.  The caller takes
 * ownership of the new trackpoint and so is responsible for ensuring
 * that it is eventually destryed.
 *
 * @param lat a value between -90 and 90 (inclusive)
 * @param lon a value between -180 (inclusive) and 180 (exclusive)
 * @param time a long
 * @return a pointer to the trackpoint created
 */
trackpoint *trackpoint_create(double lat, double lon, long time);

/**
 * Returns a copy of the given trackpoint.  The caller takes ownership
 * of the returned copy, and so is responsible for ensuring that it
 * is eventually destroyed. 
 */
trackpoint *trackpoint_copy(const trackpoint *pt);

/**
 * Destroys the given trackpoint, rendering it invalid.
 *
 * @param pt a pointer to a valid trackpoint
 */
void trackpoint_destroy(trackpoint *pt);

/**
 * Returns the location of the given trackpoint.
 *
 * @param pt a pointer to a valid track point
 * @return the location of that trackpoint
 */
location trackpoint_location(const trackpoint *pt);

/**
 * Returns the timestamp of the given trackpoint.
 *
 * @param pt a pointer to a valid track point
 * @return the timestamp of that trackpoint
 */
long trackpoint_time(const trackpoint *pt);

#endif
