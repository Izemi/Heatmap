#include "heatmap.h"
#include "track.h"
#include "trackpoint.h"
#include "string_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* Forward decls. */
static bool parse_line(const char *line, double *lat, double *lon, long *t);

/*  validate_command_line */
bool validate_command_line(int argc, char **argv, command_line *args)
{
    if (argc < 5) {
        return false;
    }
    char *endptr;

    double w = strtod(argv[1], &endptr);
    if (!isfinite(w) || w <= 0.0) return false;

    double h = strtod(argv[2], &endptr);
    if (!isfinite(h) || h <= 0.0) return false;

    const char *color_str = argv[3];
    if (strlen(color_str) == 0) return false;

    int rw = atoi(argv[4]);
    if (rw <= 0) return false;

    args->cell_width = w;
    args->cell_height = h;
    args->colors = color_str;
    args->num_colors = strlen(color_str);
    args->range_width = (size_t)rw;
    return true;
}

/* read_input, read lines until EOF, parse "lat lon time" or blank line => new segment
   returns false if parse error
 */
bool read_input(FILE *in, track *trk)
{
    if (!in || !trk) return false;
    char *line;
    while ((line = read_line(in)) != NULL) {
        // trim leading spaces
        char *p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0') {
            // blank => new segment
            track_start_segment(trk);
            free(line);
            continue;
        }
        // parse lat lon time
        double lat, lon;
        long t;
        if (!parse_line(p, &lat, &lon, &t)) {
            free(line);
            return false;
        }
        // create trackpoint
        trackpoint *pt = trackpoint_create(lat, lon, t);
        if (!pt) {
            free(line);
            return false;
        }
        track_add_point(trk, pt);
        trackpoint_destroy(pt);

        free(line);
    }
    // treat EOF as normal end
    return true;
}

static bool parse_line(const char *line, double *lat, double *lon, long *t)
{
    if (!line || !lat || !lon || !t) return false;
    double la, lo;
    long tt;
    int rc = sscanf(line, "%lf %lf %ld", &la, &lo, &tt);
    if (rc != 3) return false;
    *lat = la;
    *lon = lo;
    *t = tt;
    return true;
}

/* show_heatmap calls track_heatmap to get the 2D array */
bool show_heatmap(FILE *out, const track *trk, const command_line *args)
{
    if (!out || !trk || !args) return false;
    size_t **map = NULL;
    size_t rows = 0, cols = 0;
    track_heatmap(trk, args->cell_width, args->cell_height,
                  &map, &rows, &cols);
    if (!map) return false;
    // rows x cols map
    for (size_t r = 0; r < rows; r++) {
        for (size_t c = 0; c < cols; c++) {
            size_t val = map[r][c];
            size_t idx;
            // if val >= (num_colors-1)*range_width => last color
            size_t high_cutoff = (args->num_colors - 1) * args->range_width;
            if (val >= high_cutoff) {
                idx = args->num_colors - 1;
            } else {
                idx = val / args->range_width;
            }
            fputc(args->colors[idx], out);
        }
        fputc('\n', out);
    }
    // free map
    for (size_t r = 0; r < rows; r++) {
        free(map[r]);
    }
    free(map);
    return true;
}

/* main for Heatmap */
int main(int argc, char **argv)
{
    command_line args;
    if (!validate_command_line(argc, argv, &args)) {
        fprintf(stderr, "%s: invalid arguments\n", argv[0]);
        return 1;
    }
    track *trk = track_create();
    if (!trk) {
        fprintf(stderr, "%s: could not create track\n", argv[0]);
        return 1;
    }
    if (!read_input(stdin, trk)) {
        track_destroy(trk);
        fprintf(stderr, "%s: invalid input\n", argv[0]);
        return 1;
    }
    if (!show_heatmap(stdout, trk, &args)) {
        track_destroy(trk);
        fprintf(stderr, "%s: could not generate heatmap\n", argv[0]);
        return 1;
    }
    track_destroy(trk);
    return 0;
}
