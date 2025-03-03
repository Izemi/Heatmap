#ifndef __TRACK_UNIT_H__
#define __TRACK_UNIT_H__

#include <stddef.h>

#include "track.h"
#include "trackpoint.h"
#include "location.h"

/**
 * Makes a track from the given 2-D array of locations.  Each row will
 * be a different segment in the track.  The times of the corresponding
 * track points start with the given value and increase by 1 each track
 * point.
 *
 * @param loc a pointer to a 2-D array of locations, non-NULL
 * @param num_segs the number of rows in loc
 * @param num_pts a pointer to an array containing the number of columns
 * in each row of locs, non-NULL
 * @param t an integer
 */
track *make_track(const location **loc,
                  size_t num_segs,
                  size_t *num_pts,
                  long t);

/**
 * Tests creating an destroying a track.  Outputs an error message to
 * standard output if there were problems.
 */
void create_destroy_track();

/**
 * Tests segment and point counts for an empty track.  Outputs an error
 * message to standard output if there were problems.
 */
void empty_track_counts();

/**
 * Adds the given points to a track and checks that they are retrieved
 * correctly using track_get_point.  Outputs an error message to standard
 * output if there were problems.
 *
 * @param pts a pointer to a 2-D array of locations, non-NULL
 * @param num_segs the number of rows in loc
 * @param num_pts a pointer to an array containing the number of columns
 * in each row of locs, non-NULL
 * @param t an integer
 */
void add_and_get(const location **pts,
                 size_t num_segs,
                 size_t *num_pts,
                 long t);

/**
 * Adds the given points to a track in a single segment and checks that
 * they are retrieved correctly using track_get_point.  Outputs an error
 * message to standard output if there were problems.
 *
 * @param pts a pointer to an array of locations, non-NULL
 * @param n the number of elements in pts
 */
void add_and_get_single(const location *pts, size_t n);

/**
 * Adds the given points to a track and verifies that the counts returned
 * by track_count_segments and track_count_points are correct.  Outputs
 * an error message to standard output if there were problems.
 *
 * @param pts a pointer to a 2-D array of locations, non-NULL
 * @param num_segs the number of rows in loc
 * @param num_pts a pointer to an array containing the number of columns
 * in each row of locs, non-NULL
 * @param t an integer
 */
void check_counts(const location **pts,
                  size_t num_segs,
                  size_t *num_pts,
                  long t);

/**
 * Checks that track_merge_segments works correctly.  Creates a track
 * with the given number of segments and then merges the indicated
 * segments.  Outputs an error message to standard output if there were
 * problems.
 *
 * @param start_segments a nonnegative integer
 * @param start an integer greater than or equal to 0 and strictly less than
 * start_segments
 * @param end an integer greater than or equal to start and less than or
 * equal to start_segments
 */
void merge(size_t start_segments, size_t merge_start, size_t merge_end);

/**
 * Verifies that track_add_point creates a copy of the trackpoint passed
 * to it.  Outputs an error message if there were problems.
 */
void copy_in_add();

/**
 * Verifies that track_add_point and track_get_point work for tracks
 * of the given size.  Builds a track with the given number of points
 * on each of the given number of segments.  Locations are computed
 * arbirarily.  Outputs a message to standard output if there were problems.
 *
 * @param num_segs a positive integer
 * @param num_ptss a positive integer
 */
void many_segments_and_points(size_t num_segs, size_t num_pts);

/**
 * Verifies that track_heatmap works correctly.  Builds a track so
 * that the size and number of points in each cell of a heatmap will
 * be as given by the counts argument.  Outputs an error message to
 * standard output if the heatmap returned by track_heatmap is not
 * correct.
 *
 * @param rows a positive integer
 * @param cols a positive integer
 * @param counts a 2-D array of nonnegative integers
 * @param north a valid latitude
 * @param west a valid longitude
 */
void heatmap(size_t rows,
             size_t cols,
             size_t counts[][cols],
             double north,
             double west);

/**
 * Frees a heatmap returned by track_heatmap.
 *
 * @param map a pointer to a 2-D array of counts, non-NULL
 * @param rows the number of rows in map
 */
void free_heatmap(size_t **map, size_t rows);

#endif
