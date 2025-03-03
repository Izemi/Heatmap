/**
 * Unit tests for the track ADT.
 *
 * @author CPSC 223 Staff
 */

#include "track_unit.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "track.h"
#include "trackpoint.h"
#include "location.h"

// a 3-point segment for testing
location short_segment[] = {{41.3078680, -72.9342120},
                          {41.3078780, -72.9342340},
                          {41.3078810, -72.9342590}};

// a 4-point segment for testing
location parallel_short_segment[] = {{41.3078680, -72.8342120},
                                     {41.3078780, -72.8342340},
                                     {41.3078790, -72.8342260},
                                     {41.3078810, -72.8342590}};

// a 2-segment track for testing
const location *two_segment[] = {short_segment, parallel_short_segment};

// the lengths of the segments in the above track
size_t two_segment_lengths[] = {3, 4};

// expected counts in heatmap cells (will build a track to match)
size_t small_map_counts[][3] = {{1, 2, 3}, {4, 5, 6}};
size_t small_map_rows = 2;
size_t small_map_cols = 3;

int main(int argc, char **argv)
{
  if (argc < 2)
    {
      fprintf(stderr, "USAGE: %s test-number\n", argv[0]);
      return 1;
    }

  int test = atoi(argv[1]);
  switch (test)
    {
    case 1:
      create_destroy_track();
      break;

    case 2:
      empty_track_counts();
      break;

    case 3:
      copy_in_add();
      break;

    case 4:
      add_and_get_single(short_segment,
                         sizeof(short_segment) / sizeof(location));
      break;

    case 5:
      add_and_get(two_segment, 2, two_segment_lengths, 2000);
      break;

    case 6:
      check_counts(two_segment, 2, two_segment_lengths, 2000);
      break;

    case 7:
      // generic merge
      merge(6, 2, 4);
      break;

    case 8:
      // merge including 1st segment
      merge(6, 0, 3);
      break;

    case 9:
      // merge including only 1 segment
      merge(6, 2, 3);
      break;

    case 10:
      // the resulting heatmap doesn't cross 180/-180
      heatmap(small_map_rows, small_map_cols, small_map_counts, 45.0, 10.0);
      break;

    case 11:
      // the resulting heatmap will cross 180/-180
      heatmap(small_map_rows, small_map_cols, small_map_counts, 45.0, 179.0);
      break;

    case 12:
      {
        // read extra command line arguments for number of segments
        // and number of points on each segment
        int num_segs = -1;
        int num_pts = -1;
        if (argc >= 4)
          {
            num_segs = atoi(argv[2]);
            num_pts = atoi(argv[3]);
          }
        if (argc < 4 || num_segs <= 0 || num_pts <= 0)
          {
            fprintf(stderr,
                    "USAGE: %s %s num-segments num-points\n\twhere num-segments and num-points are positive\n",
                    argv[0],
                    argv[1]);
            return 1;
          }
        many_segments_and_points((size_t)num_segs, (size_t)num_pts);
      }
      break;

    default:
      fprintf(stderr, "%s: invalid test number %s\n", argv[0], argv[1]);
      return 1;
    }

  return 0;
}


track *make_track(const location **loc,
                  size_t num_seg,
                  size_t *num_pts,
                  long t)
{
  track *trk = track_create();
  if (trk == NULL)
    {
      return NULL;
    }

  long time = t;
  for (size_t seg = 0; seg < num_seg; seg++)
    {
      if (seg > 0)
        {
          track_start_segment(trk);
        }
      
      for (size_t i = 0; i < num_pts[seg]; i++)
        {
          trackpoint *pt = trackpoint_create(loc[seg][i].lat,
                                             loc[seg][i].lon,
                                             time++);
          if (pt != NULL)
            {
              track_add_point(trk, pt);
              trackpoint_destroy(pt);
            }
          else
            {
              fprintf(stdout, "ERROR: failed to create point\n");
              track_destroy(trk);
              return NULL;
            }
        }

    }
  return trk;
}


void create_destroy_track()
{
  track *trk = track_create();
  track_destroy(trk);
  fprintf(stdout, "PASSED\n");
}


void empty_track_counts()
{
  track *trk = track_create();
  if (trk == NULL)
    {
      fprintf(stdout, "ERROR: trackpoint_create returned NULL\n");
      return;
    }

  size_t segment_count;
  if ((segment_count = track_count_segments(trk)) != 1)
    {
      fprintf(stdout, "ERROR: empty track has %zu segments\n", segment_count);
      return;
    }

  size_t segment_size;
  if ((segment_size = track_count_points(trk, 0)) != 0)
    {
      fprintf(stdout, "ERROR: empty segment has %zu points\n", segment_size);
      return;
    }
  
  track_destroy(trk);
  fprintf(stdout, "PASSED\n");
}


void add_and_get_single(const location *pts, size_t n)
{
  long start_time = 1000;

  size_t segment_lengths[] = {n};
  add_and_get(&pts, 1, segment_lengths, start_time);
}


void add_and_get(const location **pts,
                 size_t num_segs,
                 size_t *num_pts,
                 long t)
{
  // build a track with all the points
  track *trk = make_track(pts, num_segs, num_pts, t);
  if (trk == NULL)
    {
      fprintf(stdout, "ERROR: couldn't make track\n");
      return;
    }

  // retrieve and compare all the points
  long time = t;
  for (size_t seg = 0; seg < num_segs; seg++)
    {
      for (size_t i = 0; i < num_pts[seg]; i++)
        {
          trackpoint *pt = track_get_point(trk, seg, i);
          if (pt == NULL)
            {
              fprintf(stdout, "ERROR: couldn't get point\n");
              track_destroy(trk);
              return;
            }
      
          location gotten = trackpoint_location(pt);
          if (gotten.lat != pts[seg][i].lat)
            {
              fprintf(stdout,
                      "ERROR: saved then retrieved latitude doesn't match -- %f %f\n",
                      gotten.lat,
                      pts[seg][i].lat);
              trackpoint_destroy(pt);
              track_destroy(trk);
              return;
            }
      
          if (gotten.lon != pts[seg][i].lon)
            {
              fprintf(stdout,
                      "ERROR: saved then retrieved longitude doesn't match -- %f %f\n",
                      gotten.lon,
                      pts[seg][i].lon);
              trackpoint_destroy(pt);
              track_destroy(trk);
              return;
            }
          
          if (trackpoint_time(pt) != time)
            {
              fprintf(stdout,
                      "ERROR: saved then retrieved time doesn't match -- %ld %ld\n",
                      trackpoint_time(pt),
                      time);
              trackpoint_destroy(pt);
              track_destroy(trk);
              return;
            }
          trackpoint_destroy(pt);
          
          time++;
        }
    }
  
  track_destroy(trk);
  fprintf(stdout, "PASSED\n");
}


void check_counts(const location **pts,
                  size_t num_segs,
                  size_t *num_pts,
                  long t)
{
  // make the track
  track *trk = make_track(pts, num_segs, num_pts, t);
  if (trk == NULL)
    {
      fprintf(stdout, "ERROR: couldn't make track\n");
      return;
    }

  // check number of segments
  if (track_count_segments(trk) != num_segs)
    {
      fprintf(stdout,
              "ERROR: number of segments %zu is incorrect\n",
              track_count_segments(trk));
      track_destroy(trk);
      return;
    }
  else
    {
      // try to find a row where the point count doesn't match
      size_t i = 0;
      while (i < num_segs && track_count_points(trk, i) == num_pts[i])
        {
          i++;
        }
      if (i < num_segs)
        {
          // found one -- display error
          fprintf(stdout,
                  "ERROR: number of points %zu on segment %zu is incorrect\n",
                  track_count_points(trk, i),
                  i);
          track_destroy(trk);
          return;
        }
    }
  
  track_destroy(trk);
  fprintf(stdout, "PASSED\n");
}


void merge(size_t start_segments, size_t merge_start, size_t merge_end)
{
  long start_time = 20000000L;

  // start with an empty track
  track *trk = track_create();
  if (trk == NULL)
    {
      fprintf(stdout, "ERROR: could not create track\n");
      return;
    }

  // build a track with known segments and points
  double lat = 0.0;
  double lon = 0.0;
  long time = start_time;
  for (size_t i = 0; i < start_segments; i++)
    {
      if (i > 0)
        {
          track_start_segment(trk);
        }

      // add i+1 points to segment i
      for (size_t j = 0; j < i + 1; j++)
        {
          trackpoint *pt = trackpoint_create(lat, lon, time);
          
          // arbitrarily update coordinates and timestamp of next point
          lat += 0.1;
          lon += 0.1;
          time += 10;
          if (pt == NULL)
            {
              fprintf(stdout, "ERROR: could not create trackpoint\n");
              track_destroy(trk);
              return;
            }
          track_add_point(trk, pt);

          trackpoint_destroy(pt);
        }
    }

  // do the requested merge
  track_merge_segments(trk, merge_start, merge_end);

  // compute what should have happened
  size_t deleted_segments;
  if (merge_start >= 0
      && merge_start < start_segments
      && merge_end > merge_start
      && merge_end <= start_segments)
    {
      deleted_segments = merge_end - merge_start - 1;
    }
  else
    {
      deleted_segments = 0;
    }
  
  size_t expected_segments = start_segments - deleted_segments;

  // verify number of segments
  if (track_count_segments(trk) != expected_segments)
    {
      fprintf(stdout,
              "ERROR: incorrect number of segments %zu after merge\n",
              track_count_segments(trk));
      track_destroy(trk);
      return;
    }

  // go through each segment and verify count and points on segment
  lat = 0.0;
  lon = 0.0;
  time = start_time;
  
  for (size_t i = 0; i < expected_segments; i++)
    {
      size_t expected_length = 0;
      if (i < merge_start || deleted_segments == 0)
        {
          // segment is pre-merge, so has original number of points
          expected_length = i + 1;
        }
      else if (i == merge_start)
        {
          // segment is the merged segment, so add points on merged segments
          for (size_t s = merge_start; s < merge_end; s++)
            {
              expected_length += (s + 1);
            }
        }
      else
        {
          // segment is post-merge, so has original number of points, but
          // its index changed
          expected_length = i + deleted_segments + 1;
        }

      // verify point count
      if (track_count_points(trk, i) != expected_length)
        {
          fprintf(stdout,
                  "ERROR: segment %zu has size %zu after merge\n",
                  i,
                  track_count_points(trk, i));
          track_destroy(trk);
          return;
        }

      // go through points and make sure coordinates and timestamp match
      for (size_t j = 0; j < expected_length; j++)
        {
          trackpoint *pt = track_get_point(trk, i, j);
          if (pt == NULL)
            {
              fprintf(stdout,
                      "ERROR: segment %zu point %zu is NULL after merge\n",
                      i,
                      j);
              track_destroy(trk);
              return;
            }

          location loc = trackpoint_location(pt);
          if (loc.lat != lat || loc.lon != lon)
            {
              fprintf(stdout,
                      "ERROR: got point %f %f from segment %zu point %zu after merge\n",
                      loc.lat,
                      loc.lon,
                      i,
                      j);
              track_destroy(trk);
              trackpoint_destroy(pt);
              return;
            }

          if (trackpoint_time(pt) != time)
            {
              fprintf(stdout,
                      "ERROR: got time %ld from segment %zu point %zu after merge\n",
                      trackpoint_time(pt),
                      i,
                      j);
              track_destroy(trk);
              trackpoint_destroy(pt);
              return;
            }

          trackpoint_destroy(pt);

          // update expected values for next point
          lat += 0.1;
          lon += 0.1;
          time += 10;
        }
    }

  track_destroy(trk);
  fprintf(stdout, "PASSED\n");
}

void copy_in_add()
{
  track *trk = track_create();

  if (trk == NULL)
    {
      fprintf(stdout, "ERROR: couldn't make track\n");
      return;
    }

  // add an arbitrary point
  trackpoint *pt = trackpoint_create(47.5, -125.0, 1000);
  track_add_point(trk, pt);

  trackpoint_destroy(pt);

  // make a new point to increase probablity we overwrite the old point
  pt = trackpoint_create(0.0, 0.0, 0);
  trackpoint_destroy(pt);

  pt = track_get_point(trk, 0, 0);
  if (pt == NULL)
    {
      fprintf(stdout, "ERROR: could not retrieve point\n");
      track_destroy(trk);
      return;
    }

  // check that we get the original values back
  location loc = trackpoint_location(pt);
  if (loc.lat != 47.5 || loc.lon != -125.0 || trackpoint_time(pt) != 1000)
    {
      fprintf(stdout,
              "ERROR: trackpoint in track affected by destruction of original\n");
      track_destroy(trk);
      trackpoint_destroy(pt);
      return;
    }

  track_destroy(trk);
  trackpoint_destroy(pt);
  fprintf(stdout, "PASSED\n");
}

void many_segments_and_points(size_t num_segs, size_t num_pts)
{
  // this is similar to add_and_get but there's no array passed in giving
  // the locations; instead we compute locations in the test function to
  // facilitate creating LOTS of segments and points
  track *trk = track_create();

  if (trk == NULL)
    {
      fprintf(stdout, "ERROR: couldn't make track\n");
      return;
    }

  // make a track with a bunch of points and segments
  long time = 0;
  for (size_t s = 0; s < num_segs; s++)
    {
      double lat = (double)s / (1 << 20); // each segment follows lat s/2^20
      if (s > 0)
        {
          track_start_segment(trk);
        }

      for (size_t p = 0; p < num_pts; p++)
        {
          double lon = (double)p / (1 << 20); // longitude is p/2^20 each time
          trackpoint *pt =  trackpoint_create(lat, lon, time);
          time++;
          track_add_point(trk, pt);
          trackpoint_destroy(pt);
        }
    }

  // retrieve all the points to make sure they're where we put them
  long expected_time = 0;
  for (size_t s = 0; s < num_segs; s++)
    {
      double expected_lat = (double)s / (1 << 20);
      for (size_t p = 0; p < num_pts; p++)
        {
          double expected_lon = (double)p / (1 << 20);
          trackpoint *pt =  track_get_point(trk, s, p);
          if (trackpoint_time(pt) != expected_time)
            {
              fprintf(stdout,
                      "ERROR: got time %ld when expecting %ld\n",
                      trackpoint_time(pt),
                      expected_time);
              track_destroy(trk);
              return;
            }
          location loc = trackpoint_location(pt);
          if (loc.lat != expected_lat || loc.lon != expected_lon)
            {
              fprintf(stdout,
                      "ERROR: got location %lf %lf when expecting %lf %lf\n",
                      loc.lat,
                      expected_lat,
                      loc.lon,
                      expected_lon);
              track_destroy(trk);
              return;
            }
          expected_time++;
          trackpoint_destroy(pt);
        }
    }
  
  track_destroy(trk);
  fprintf(stdout, "PASSED\n");
}

void heatmap(size_t rows,
             size_t cols,
             size_t counts[][cols],
             double north,
             double west)
{
  double cell_width = 1.0;
  double cell_height = 1.0;
  long time = 1000;
  
  track *trk = track_create();
  
  if (trk == NULL)
    {
      fprintf(stdout, "ERROR: couldn't make track\n");
      return;
    }

  bool made_north = false;
  bool made_west = false;
  
  for (size_t r = 0; r < rows; r++)
    {
      for (size_t c = 0; c < cols; c++)
        {
          for (size_t p = 0; p < counts[r][c]; p++)
            {
              // compute lat and lon for point to add so all are distinct
              // and there is always a point on the north edge and
              // a point on the west edge, otherwise all points are
              // not on the borders
              double lat;
              double lon;
              if (p == 0 && r == 0 && c == 0)
                {
                  // first point in NW cell goes on NW corner
                  lat = north;
                  lon = west;
                  made_north = true;
                  made_west = true;
                }
              else if (p == 0 && r == 0 && !made_north)
                {
                  // first point on top row goes on top if there wasn't
                  // already something in the NW corner of NW cell
                  lat = north;
                  lon = west
                    + cell_width * c
                    + (p + 1) * (cell_width / (counts[r][c] + 1));
                  made_north = true;
                }
              else if (p == 0 && c == 0 && !made_west)
                {
                  // first point in first column goes on west edge if
                  // there wasn't already something in NW corner of NW cell
                  lat = north
                    - cell_height * r
                    - (p + 1) * (cell_height / (counts[r][c] + 1));
                  lon = west;
                  made_west = true;
                }
              else
                {
                  // other points go diagonally through cell
                  lat = north
                    - cell_height * r
                    - (p + 1) * (cell_height / (counts[r][c] + 1));
                  lon = west
                    + cell_width * c
                    + (p + 1) * (cell_width / (counts[r][c] + 1));
                }

              if (lon >= 180.0)
                {
                  // normalize longitude
                  lon -= 360.0;
                }
              
              trackpoint *pt = trackpoint_create(lat, lon, time);
              track_add_point(trk, pt);
              trackpoint_destroy(pt);
              time++;
            }
        }
    }

  // get the heatmap
  size_t **map;
  size_t map_rows;
  size_t map_cols;
  track_heatmap(trk, 1.0, 1.0, &map, &map_rows, &map_cols);

  if (map == NULL)
    {
      fprintf(stdout, "ERROR: couldn't make heatmap\n");
      track_destroy(trk);
      return;
    }

  // verify that the dimensions are what we expected
  if (map_rows != rows || map_cols != cols)
    {
      fprintf(stdout,
              "ERROR: heatmap dimensions %zu %zu incorrect\n",
              map_rows,
              map_cols);
      free_heatmap(map, map_rows);
      track_destroy(trk);
      return;
    }

  // verify that each cell in the heatmap has the expected count
  for (size_t r = 0; r < map_rows; r++)
    {
      if (map[r] == NULL)
        {
          fprintf(stdout,
                  "ERROR: heatmap row %zu is NULL\n",
                  r);
          free_heatmap(map, map_rows);
          track_destroy(trk);
          return;
        }
      
      for (size_t c = 0; c < map_cols; c++)
        {
          if (map[r][c] != counts[r][c])
            {
              fprintf(stdout,
                      "ERROR: heatmap entry %zu %zu is incorrect %zu\n",
                      r,
                      c,
                      map[r][c]);
              free_heatmap(map, map_rows);
              track_destroy(trk);
              return;
            }
        }
    }

  free_heatmap(map, map_rows);
  track_destroy(trk);
  fprintf(stdout, "PASSED\n");
}

void free_heatmap(size_t **map, size_t rows)
{
  for (size_t r = 0; r < rows; r++)
    {
      free(map[r]);
    }
  free(map);
}
