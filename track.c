#include "track.h"
#include "trackpoint.h"
#include "location.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

/* Segment data structure */
typedef struct segment {
    trackpoint **pts;
    size_t size;       // how many points
    size_t capacity;   // capacity
    double length;     // sum of distances
} segment;

struct track {
    segment **segments;
    size_t seg_size;      // used
    size_t seg_capacity;  // capacity
};

#define INITIAL_SEG_CAPACITY 4
#define INITIAL_PT_CAPACITY  10

/* segment_create/destroy */
static segment *segment_create(void)
{
    segment *seg = malloc(sizeof(*seg));
    if (!seg) return NULL;
    seg->pts = malloc(INITIAL_PT_CAPACITY * sizeof(*seg->pts));
    if (!seg->pts) {
        free(seg);
        return NULL;
    }
    seg->size = 0;
    seg->capacity = INITIAL_PT_CAPACITY;
    seg->length = 0.0;
    return seg;
}

static void segment_destroy(segment *seg)
{
    if (!seg) return;
    for (size_t i = 0; i < seg->size; i++) {
        trackpoint_destroy(seg->pts[i]);
    }
    free(seg->pts);
    free(seg);
}

/* distance helper */
static double distance_between(const trackpoint *a, const trackpoint *b)
{
    if (!a || !b) return 0.0;
    location la = trackpoint_location(a);
    location lb = trackpoint_location(b);
    return location_distance(&la, &lb);
}

/* add point to segment, amortized O(1) */
static void segment_add_point(segment *seg, const trackpoint *pt)
{
    if (seg->size == seg->capacity) {
        size_t new_cap = seg->capacity * 2;
        trackpoint **new_pts = realloc(seg->pts, new_cap * sizeof(*new_pts));
        seg->pts = new_pts; 
        seg->capacity = new_cap;
    }
    trackpoint *copy = trackpoint_copy(pt);
    if (seg->size > 0) {
        trackpoint *prev = seg->pts[seg->size - 1];
        seg->length += distance_between(prev, copy);
    }
    seg->pts[seg->size++] = copy;
}

/* Track ADT */
track *track_create(void)
{
    track *trk = malloc(sizeof(*trk));
    if (!trk) return NULL;

    trk->seg_capacity = INITIAL_SEG_CAPACITY;
    trk->segments = malloc(trk->seg_capacity * sizeof(*trk->segments));
    if (!trk->segments) {
        free(trk);
        return NULL;
    }
    trk->seg_size = 1;
    trk->segments[0] = segment_create();
    if (!trk->segments[0]) {
        free(trk->segments);
        free(trk);
        return NULL;
    }
    return trk;
}

void track_destroy(track *trk)
{
    if (!trk) return;
    for (size_t i = 0; i < trk->seg_size; i++) {
        segment_destroy(trk->segments[i]);
    }
    free(trk->segments);
    free(trk);
}

size_t track_count_segments(const track *trk)
{
    if (!trk) return 0;
    return trk->seg_size;
}

size_t track_count_points(const track *trk, size_t seg_index)
{
    if (!trk) return 0;
    if (seg_index >= trk->seg_size) return 0;
    return trk->segments[seg_index]->size;
}

trackpoint *track_get_point(const track *trk, size_t seg_index, size_t pt_index)
{
    if (!trk) return NULL;
    if (seg_index >= trk->seg_size) return NULL;
    segment *seg = trk->segments[seg_index];
    if (pt_index >= seg->size) return NULL;
    return trackpoint_copy(seg->pts[pt_index]);
}

/* track_get_lengths */
double *track_get_lengths(const track *trk)
{
    if (!trk) return NULL;
    if (trk->seg_size == 0) return NULL;

    double *arr = malloc(trk->seg_size * sizeof(*arr));
    if (!arr) return NULL; // assume success
    for (size_t i = 0; i < trk->seg_size; i++) {
        arr[i] = trk->segments[i]->length;
    }
    return arr;
}

void track_add_point(track *trk, const trackpoint *pt)
{
    if (!trk || !pt) return;
    segment_add_point(trk->segments[trk->seg_size - 1], pt);
}

void track_start_segment(track *trk)
{
    if (!trk) return;
    if (trk->seg_size == trk->seg_capacity) {
        size_t new_cap = trk->seg_capacity * 2;
        segment **new_arr = realloc(trk->segments, new_cap * sizeof(*new_arr));
        trk->segments = new_arr; 
        trk->seg_capacity = new_cap;
    }
    trk->segments[trk->seg_size] = segment_create();
    trk->seg_size++;
}

/* merges [start..end-1], O(n) worst-case */
void track_merge_segments(track *trk, size_t start, size_t end)
{
    if (!trk) return;
    if (start >= end - 1) return;
    if (start >= trk->seg_size) return;
    if (end > trk->seg_size) end = trk->seg_size;

    segment *base = trk->segments[start];
    for (size_t i = start + 1; i < end; i++) {
        segment *curr = trk->segments[i];
        if (curr->size == 0) continue;
        // distance from last of base to first of curr
        if (base->size > 0 && curr->size > 0) {
            trackpoint *bLast = base->pts[base->size - 1];
            trackpoint *cFirst = curr->pts[0];
            base->length += distance_between(bLast, cFirst);
        }
        // append
        for (size_t j = 0; j < curr->size; j++) {
            trackpoint *p = curr->pts[j];
            if (j == 0 && base->size > 0) {
                if (base->size == base->capacity) {
                    size_t new_cap = base->capacity * 2;
                    base->pts = realloc(base->pts, new_cap*sizeof(*base->pts));
                    base->capacity = new_cap;
                }
                base->pts[base->size++] = trackpoint_copy(p);
            } else {
                if (base->size == base->capacity) {
                    size_t new_cap = base->capacity * 2;
                    base->pts = realloc(base->pts, new_cap*sizeof(*base->pts));
                    base->capacity = new_cap;
                }
                trackpoint *cpy = trackpoint_copy(p);
                trackpoint *prev = base->pts[base->size - 1];
                base->length += distance_between(prev, cpy);
                base->pts[base->size++] = cpy;
            }
        }
    }
    // free merged
    size_t removed = end - 1 - start;
    for (size_t i = start+1; i < end; i++) {
        segment_destroy(trk->segments[i]);
    }
    // shift
    for (size_t i = end; i < trk->seg_size; i++) {
        trk->segments[i-removed] = trk->segments[i];
    }
    trk->seg_size -= removed;
}

/* track_heatmap */
static size_t track_total_points(const track *trk)
{
    size_t total = 0;
    for (size_t i = 0; i < trk->seg_size; i++) {
        total += trk->segments[i]->size;
    }
    return total;
}

/* gather all points in an array for convenience */
static trackpoint **gather_points(const track *trk, size_t total)
{
    trackpoint **arr = malloc(total*sizeof(*arr));
    size_t idx = 0;
    for (size_t s = 0; s < trk->seg_size; s++) {
        segment *seg = trk->segments[s];
        for (size_t p = 0; p < seg->size; p++) {
            arr[idx++] = seg->pts[p];
        }
    }
    return arr;
}

/* normalize lon => [-180,180). */
static double norm_lon(double lon)
{
    while (lon >= 180.0) lon -= 360.0;
    while (lon < -180.0) lon += 360.0;
    return lon;
}

/* measure eastward difference from west->x in [0..360]. */
static double eastward_diff(double west, double x)
{
    double w = norm_lon(west);
    double d = norm_lon(x);
    if (d >= w) {
        return d - w;
    } else {
        return 360.0 - (w - d);
    }
}

/* qsort comparator for doubles */
static int cmp_double(const void *a, const void *b)
{
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/* find minimal wedge for all n longitudes. */
static void find_min_longitude_wedge(trackpoint **pts, size_t n,
                                     double *west_out, double *east_out)
{
    if (n == 0) {
        *west_out = 0.0;
        *east_out = 0.0;
        return;
    }
    if (n == 1) {
        location loc = trackpoint_location(pts[0]);
        double w = norm_lon(loc.lon);
        *west_out = w;
        *east_out = w;
        return;
    }
    double *L = malloc(n*sizeof(*L));
    for (size_t i = 0; i < n; i++) {
        location loc = trackpoint_location(pts[i]);
        L[i] = norm_lon(loc.lon);
    }
    qsort(L,n,sizeof(double),cmp_double);

    double *Lext = malloc(2*n*sizeof(*Lext));
    for (size_t i = 0; i < n; i++) {
        Lext[i] = L[i];
        Lext[i+n] = L[i] + 360.0;
    }

    double best_size = 1e9;
    double best_west = L[0];
    for (size_t i = 0; i < n; i++) {
        double w = Lext[i];
        double e = Lext[i + n - 1];
        double size = e - w;
        if (size < best_size) {
            best_size = size;
            best_west = w;
        } else if (fabs(size - best_size) < 1e-14) {
            if (w < best_west) {
                best_west = w;
            }
        }
    }
    free(Lext);
    free(L);

    double wnorm = norm_lon(best_west);
    double enorm = wnorm + best_size;
    *west_out = wnorm;
    *east_out = enorm;
}


/* picks row for lat */
static size_t pick_row(double lat, double north, double cell_height, size_t rows)
{
    if (rows == 1) return 0; // trivial
    double exact = (north - lat) / cell_height; // how many cells down
    // clamp
    if (exact < 0.0) exact = 0.0;
    double max_r = (double)(rows - 1);
    if (exact > max_r) exact = max_r;

    // integer part
    double f = floor(exact);
    double frac = exact - f;
    size_t base = (size_t)f;

    
    const double EPS = 1e-14;
    if (frac < EPS && base < rows-1) {
        // lat is exactly top boundary 
        return base;
    } else if ((1.0 - frac) < EPS) {
        // lat is exactly the bottom boundary 
        // unless base=rows-1 
        if (base < rows-1) {
            return base + 1;
        } else {
            return base;
        }
    } else {
        // normal case
        return base;
    }
}

/* picks col for lon */
static size_t pick_col(double lon, double west, double cell_width, size_t cols)
{
    if (cols == 1) return 0;
    double diff = eastward_diff(west, lon);
    double exact = diff / cell_width;
    if (exact < 0.0) exact = 0.0;
    double max_c = (double)(cols - 1);
    if (exact > max_c) exact = max_c;

    double f = floor(exact);
    double frac = exact - f;
    size_t base = (size_t)f;

    const double EPS = 1e-14;
    if (frac < EPS && base < cols-1) {
        // boundary belongs to col=base 
        return base;
    } else if ((1.0 - frac) < EPS) {
        // boundary belongs to next col if not last
        if (base < cols-1) {
            return base + 1;
        } else {
            return base;
        }
    } else {
        return base;
    }
}

void track_heatmap(const track *trk,
                   double cell_width, double cell_height,
                   size_t ***map_out, size_t *rows_out, size_t *cols_out)
{
    if (!trk || !map_out || !rows_out || !cols_out) return;

    size_t total = track_total_points(trk);
    if (total == 0) {
        // empty 
        size_t **one = malloc(sizeof(*one));
        one[0] = calloc(1, sizeof(**one));
        *map_out = one;
        *rows_out = 1;
        *cols_out = 1;
        return;
    }

    trackpoint **pts = gather_points(trk, total);

    // bounding lat
    double north = -90.0, south = 90.0;
    for (size_t i = 0; i < total; i++) {
        location loc = trackpoint_location(pts[i]);
        if (loc.lat > north) north = loc.lat;
        if (loc.lat < south) south = loc.lat;
    }

    // minimal wedge for lon
    double west, east;
    find_min_longitude_wedge(pts, total, &west, &east);

    double lat_span = north - south;
    if (lat_span < 0) lat_span = 0;
    // rows = ceil
    size_t rows = (size_t)ceil(lat_span / cell_height);
    if (rows < 1) rows = 1;

    double lon_span = east - west;
    if (lon_span < 0) lon_span += 360.0;
    if (lon_span > 360.0) lon_span = 360.0;
    size_t cols = (size_t)ceil(lon_span / cell_width);
    if (cols < 1) cols = 1;

    // allocate map
    size_t **map = malloc(rows*sizeof(*map));
    for (size_t r = 0; r < rows; r++) {
        map[r] = calloc(cols, sizeof(**map));
    }

    // place each point
    for (size_t i = 0; i < total; i++) {
        location loc = trackpoint_location(pts[i]);
        double lat = loc.lat;
        double lon = norm_lon(loc.lon);

        size_t rr = pick_row(lat, north, cell_height, rows);
        if (rr >= rows) rr = rows-1; 
        size_t cc = pick_col(lon, west, cell_width, cols);
        if (cc >= cols) cc = cols-1; 

        map[rr][cc]++;
    }

    *map_out = map;
    *rows_out = rows;
    *cols_out = cols;

    free(pts);
}
