#ifndef __HEATMAP_H__
#define __HEATMAP_H__

#include <stdio.h>
#include <stdbool.h>

#include "track.h"

/**
 * Parameters read from the command line.
 */

typedef struct {
  double cell_width;
  double cell_height;
  const char *colors; // or "colors" -- no scare quotes in C identifiers!
  size_t num_colors;
  size_t range_width;
} command_line;

/**
 * Parses the given array of strings and parses them according to the
 * rules for valid command-line arguments.  If there are invalid arguments
 * then the return value is false and the contents of the struct args
 * points to are undefined.
 *
 * @param argc a positive integer
 * @param argv an array of argc strings
 * @param a pointer to a command_line struct, non-NULL
 * @return true if and only if the arguments were parsed without error
 */
bool validate_command_line(int argc, char **argv, command_line *args);

/**
 * Reads track points and segments from the given stream and adds them
 * to the given track.  The input should be in the format "lat lon time"
 * with blank lines between segments.  The first point is added to the
 * current segment of the track.
 *
 * @param trk a pointer to a track, non-NULL
 * @return true if and only if the input was read with no errors
 */
bool read_input(FILE *in, track *trk);

/**
 * Prints a heatmap for the given track to the given output stream according
 * to the given parameters.
 *
 * @param out an output stream, non-NULL
 * @param track a pointer to a track, non-NULL
 * @param args a pointer to valid parsed command-line arguments, non-NULL
 * @return true if and only if the heatmap was created and displayed w/o error
 */
bool show_heatmap(FILE *out, const track *trk, const command_line *args);

#endif
