#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.141592653
#endif

#define USE_KMEANS   // use k-means or k-means++
#define PTS 100000   // number of points
#define K 6          // number of clusters

#define W 400       // width and height of output eps file
#define H 400

// https://rosettacode.org/wiki/K-means%2B%2B_clustering

typedef struct{ double x, y; int group; } point_t, *point;

double randf(double m)
{
	return m * rand() / (RAND_MAX - 1.);
}

point gen_xy(int count, double radius)
{
	double ang, r;
	point p = new point_t[count];
	point pt = new point_t[count];

	/* note: this is not a uniform 2-d distribution */
	for (p = pt + count; p-- > pt;)
	{
		ang = randf(2 * M_PI);
		r = randf(radius);
		p->x = r * cos(ang);
		p->y = r * sin(ang);
	}

	return pt;
}

inline double dist2(point a, point b)
{
	double x = a->x - b->x;
	double y = a->y - b->y;
	return x*x + y*y;
}

inline int
nearest(point pt, point cent, int n_cluster, double *d2)
{
	int i, min_i;
	point c;
	double d, min_d;

	for (c = cent, i = 0; i < n_cluster; i++, c++)
	{
		min_d = HUGE_VAL;
		min_i = pt->group;
		for (c = cent, i = 0; i < n_cluster; i++, c++)
		{
			if (min_d > (d = dist2(c, pt)))
			{
				min_d = d; min_i = i;
			}
		}
	}

	if (d2) *d2 = min_d;
	return min_i;
}

void kpp(point pts, int len, point cent, int n_cent)
{
	int i, j;
	int n_cluster;
	double sum;
	double *d = new double[len];

	point p, c;
	cent[0] = pts[rand() % len];
	for (n_cluster = 1; n_cluster < n_cent; n_cluster++) 
	{
		sum = 0;
		for (j = 0, p = pts; j < len; j++, p++) 
		{
			nearest(p, cent, n_cluster, d + j);
			sum += d[j];
		}
		sum = randf(sum);
		for (j = 0, p = pts; j < len; j++, p++) 
		{
			if ((sum -= d[j]) > 0) 
				continue;
			cent[n_cluster] = pts[j];
			break;
		}
	}
	for (j = 0, p = pts; j < len; j++, p++)
		p->group = nearest(p, cent, n_cluster, 0);
	delete[] d;
}

point lloyd(point pts, int len, int n_cluster)
{
	int i, j, min_i;
	int changed;
	point cent = new point_t[n_cluster];
	point p, c;


#ifdef USE_KMEANS
	/* assign init grouping randomly */
	for (j = 0, p = pts; j < len; j++, p++) 
		//p->group = j % n_cluster;
		p->group = rand() % n_cluster;
	
#else
	/* or call k++ init */
	kpp(pts, len, cent, n_cluster);
#endif // USE_KMEANSPP


	do 
	{
		/* group element for centroids are used as counters */
		for (c = cent, i = 0; i < n_cluster; i++, c++) 
		{ 
			c->group = 0;
			c->x = c->y = 0; 
		}
		for (j = 0, p = pts; j < len; j++, p++) 
		{
			c = cent + p->group;
			c->group++;
			c->x += p->x; c->y += p->y;
		}
		for (c = cent, i = 0; i < n_cluster; i++, c++) 
		{ 
			c->x /= c->group; 
			c->y /= c->group;
		}

		changed = 0;
		/* find closest centroid of each point */
		for (j = 0, p = pts; j < len; j++, p++) 
		{
			min_i = nearest(p, cent, n_cluster, 0);
			if (min_i != p->group) 
			{
				changed++;
				p->group = min_i;
			}
		}
	} while (changed > (len >> 10)); /* stop when 99.9% of points are good */

	for (c = cent, i = 0; i < n_cluster; i++, c++) 
	{ 
		c->group = i;
	}

	return cent;
}

void print_eps(point pts, int len, point cent, int n_cluster)
{
	int i, j;
	point p, c;
	double min_x, max_x, min_y, max_y, scale, cx, cy;
	double *colors = new double[n_cluster * 3]; 

	for (c = cent, i = 0; i < n_cluster; i++, c++) 
	{
		colors[3 * i + 0] = (3 * (i + 1) % 11) / 11.;
		colors[3 * i + 1] = (7 * i % 11) / 11.;
		colors[3 * i + 2] = (9 * i % 11) / 11.;
	}

	max_x = max_y = -(min_x = min_y = HUGE_VAL);
	for (j = 0, p = pts; j < len; j++, p++) 
	{
		if (max_x < p->x) max_x = p->x;
		if (min_x > p->x) min_x = p->x;
		if (max_y < p->y) max_y = p->y;
		if (min_y > p->y) min_y = p->y;
	}
	scale = W / (max_x - min_x);
	if (scale > H / (max_y - min_y)) 
		scale = H / (max_y - min_y);
	cx = (max_x + min_x) / 2;
	cy = (max_y + min_y) / 2;

#ifdef USE_KMEANS
	FILE *stream = fopen("kmeans.eps", "w");
#else
	FILE *stream = fopen("kmeanspp.eps", "w");
#endif // USE_KMEANS
	
	fprintf(stream, "%%!PS-Adobe-3.0\n%%%%BoundingBox: -5 -5 %d %d\n", W + 10, H + 10);
	fprintf(stream, "/l {rlineto} def /m {rmoveto} def\n"
		"/c { .25 sub exch .25 sub exch .5 0 360 arc fill } def\n"
		"/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
		"	gsave 1 setgray fill grestore gsave 3 setlinewidth"
		" 1 setgray stroke grestore 0 setgray stroke }def\n"
	);
	for (c = cent, i = 0; i < n_cluster; i++, c++)
	{
		fprintf(stream, "%g %g %g setrgbcolor\n", colors[3 * i], colors[3 * i + 1], colors[3 * i + 2]);
		for (j = 0, p = pts; j < len; j++, p++) 
		{
			if (p->group != i) 
				continue;
			fprintf(stream, "%.3f %.3f c\n", (p->x - cx) * scale + W / 2, (p->y - cy) * scale + H / 2);
		}
		fprintf(stream, "\n0 setgray %g %g s\n", (c->x - cx) * scale + W / 2, (c->y - cy) * scale + H / 2);
	}
	fprintf(stream, "\n%%%%EOF");
	fclose(stream);
	delete[] colors;
}


int main()
{
	int i;
	point v = gen_xy(PTS, 10);
	point c = lloyd(v, PTS, K);
	print_eps(v, PTS, c, K);
	delete v;
	delete c;
	return 0;
}