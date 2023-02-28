#include "plane_fitting.h"
#include "vec3.h"

#include "las_point_cloud.h"
#include "las_source.h"
#include "las_tree.h"

#include "las_normal.h"

void estim_unoriented_nml(const Vec3 *pos, size_t point_num, Vec3 *nml,
			  float *qual, const LasTree &tree,
			  void (*cb)(float progress), int probes)
{
	TArray<unsigned> knn_idx(probes);
	TArray<Vec3> neigh(probes);
	TArray<float> sqr_dist(probes);
	size_t tick = point_num / 100;
	size_t step = 0;
	for (size_t i = 0; i < point_num; ++i) {
		const float *query = &pos[i].x;
		tree.knnSearch(query, probes, &knn_idx[0], &sqr_dist[0]);
		for (int j = 0; j < probes; ++j) {
			neigh[j] = pos[knn_idx[j]];
		}
		qual[i] = fit_plane_vcg<float>(&neigh[0], probes, nml[i]);
		if (cb) {
			step += 1;
			if (step == tick) {
				cb((float)i * 100 / point_num);
				step = 0;
			}
		}
	}
}

size_t orient_nml_with_z(Vec3 *nml, EOrient *oriented, size_t point_num)
{
	size_t unsettled = 0;
	for (size_t i = 0; i < point_num; ++i) {
		if (fabs(nml[i].z) > NML_Z_THRESH) {
			oriented[i] = EPositiveZ;
			if (nml[i].z < 0)
				nml[i] *= -1;
		} else {
			unsettled++;
		}
	}
	return (unsettled);
}

size_t orient_nml_with_scan(const LasPoint *points, size_t point_num,
			    const SourceFlightLine *fls, size_t source_num,
			    const float *qual, Vec3 *nml, EOrient *oriented,
			    float tol)
{
	size_t unsettled = 0;
	for (size_t i = 0; i < point_num; ++i) {
		if (oriented[i] >= EScanline)
			continue;
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

size_t propagate_nml_once(const Vec3 *pos, size_t point_num,
			  const LasTree &tree, const float *qual, Vec3 *nml,
			  EOrient *oriented, size_t probes, float tol)
{
	size_t newly_settled = 0;
	TArray<unsigned> knn_idx(probes);
	TArray<float> sqr_dist(probes);
	for (size_t i = 0; i < point_num; ++i) {
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
	for (size_t i = 0; i < point_num; ++i) {
		if (oriented[i] == ETmpPlague)
			oriented[i] = EPlague;
	}
	return (newly_settled);
};
