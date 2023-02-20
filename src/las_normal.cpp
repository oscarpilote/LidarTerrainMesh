#include <assert.h>

#include "las_point_cloud.h"
#include "las_source.h"
#include "las_tree.h"
#include "plane_fitting.h"
#include "vec3.h"

#include "las_normal.h"

void estimate_unoriented_normals(const struct LasCoords &coords,
				 const LasTree &tree, TArray<Vec3> &nml,
				 TArray<float> &qual, int probes)
{
	TArray<unsigned> knn_idx(probes);
	TArray<Vec3> neigh(probes);
	TArray<float> sqr_dist(probes);
	printf(" [");
	fflush(stdout);
	size_t step = coords.size / 20;
	size_t done = 0;
	for (size_t i = 0; i < coords.size; ++i) {
		const float *query = &coords[i].x;
		tree.knnSearch(query, probes, &knn_idx[0], &sqr_dist[0]);
		for (int j = 0; j < probes; ++j) {
			neigh[j] = coords[knn_idx[j]];
		}
		qual[i] = fit_plane_vcg<float>(&neigh[0], probes, nml[i]);
		done += 1;
		if (done == step) {
			printf("-");
			fflush(stdout);
			done = 0;
		}
	}
	printf("]\n");
}

size_t orient_normals_with_scanlines(const TArray<struct LasPoint> &points,
				     const TArray<struct SourceFlightLine> &fls,
				     TArray<float> &qual, TArray<Vec3> &nml,
				     TArray<EOrient> &oriented, float tol)
{
	size_t unsettled = 0;
	for (size_t i = 0; i < points.size; ++i) {
		const SourceFlightLine &fl = fls[points[i].source_idx];
		if (!fl.is_valid) {
			unsettled++;
			continue;
		}
		float scan_angle = points[i].scan_angle * M_PI / 180;
		float theta = fl.theta_across;
		Vec3 beam;
		beam.x = cos(theta) * sin(scan_angle);
		beam.y = sin(theta) * sin(scan_angle);
		beam.z = -cos(scan_angle);
		float test = dot(nml[i], beam);
		if (fabs(test) > tol + 2 * (1 - qual[i])) {
			oriented[i] = EScanline;
			if (test > 0)
				nml[i] *= -1;
		} else {
			unsettled++;
		}
	}
	return (unsettled);
}

size_t propagate_normals_once(const struct LasCoords &pos, const LasTree &tree,
			      TArray<float> &qual, TArray<Vec3> &nml,
			      TArray<EOrient> &oriented, size_t probes,
			      float tol)
{
	size_t newly_settled = 0;
	TArray<unsigned> knn_idx(probes);
	TArray<float> sqr_dist(probes);
	for (size_t i = 0; i < nml.size; ++i) {
		if (oriented[i] >= EPlague)
			continue;
		const float *query = &pos[i].x;
		tree.knnSearch(query, probes, &knn_idx[0], &sqr_dist[0]);
		float majority = 0;
		for (size_t j = 0; j < probes; ++j) {
			unsigned k = knn_idx[j];
			if (oriented[k] >= EPlague) {
				float test = dot(nml[k], nml[i]);
				if (fabs(test) > (1 + tol - qual[i])) {
					majority += qual[k] * test;
				}
			}
		}
		if (majority != 0) {
			oriented[i] = ETmpPlague;
			newly_settled++;
			if (majority < 0)
				nml[i] *= -1;
		}
	}
	for (size_t i = 0; i < nml.size; ++i) {
		if (oriented[i] == ETmpPlague)
			oriented[i] = EPlague;
	}
	return (newly_settled);
};
