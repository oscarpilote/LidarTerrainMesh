#include <pthread.h>

#include "kd_tree.h"
#include "plane_fitting.h"
#include "vec3.h"

#include "las_point_cloud.h"
#include "las_source.h"

#include "las_normal.h"

#define NML_NUMTHREAD 8

#if 0
static void estim_unoriented_nml_on_thread(const Vec3 *pos, size_t point_num,
					   Vec3 *nml, float *qual,
					   const KdTree<3> &tree,
					   void (*cb)(float progress),
					   int probes)
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
#endif

struct NMLArg1 {
	const Vec3 *pos;
	size_t point_start;
	size_t point_stop;
	Vec3 *nml;
	float *qual;
	const KdTree<3> *tree;
	int probes;
	void (*cb)(float);
};

static void *estim_unoriented_nml_threaded(void *arg)
{
	const Vec3 *pos = ((NMLArg1 *)arg)->pos;
	size_t point_start = ((NMLArg1 *)arg)->point_start;
	size_t point_stop = ((NMLArg1 *)arg)->point_stop;
	Vec3 *nml = ((NMLArg1 *)arg)->nml;
	float *qual = ((NMLArg1 *)arg)->qual;
	const KdTree<3> *tree = ((NMLArg1 *)arg)->tree;
	int probes = ((NMLArg1 *)arg)->probes;
	void (*cb)(float) = ((NMLArg1 *)arg)->cb;

	size_t payload(point_stop - point_start);
	size_t tick = (point_stop - point_start) * 0.01f;
	size_t step = 0;

	TArray<unsigned> knn_idx(probes);
	TArray<Vec3> neigh(probes);
	TArray<float> sqr_dist(probes);
	for (size_t i = point_start; i < point_stop; ++i) {
		const float *query = &pos[i].x;
		tree->knnSearch(query, probes, &knn_idx[0], &sqr_dist[0]);
		for (int j = 0; j < probes; ++j) {
			neigh[j] = pos[knn_idx[j]];
		}
		qual[i] = fit_plane_vcg<float>(&neigh[0], probes, nml[i]);
		if (cb) {
			step += 1;
			if (step == tick) {
				cb((float)i * 100 / payload);
				step = 0;
			}
		}
	}
	pthread_exit(NULL);
}

void estim_unoriented_nml(const Vec3 *pos, size_t point_num, Vec3 *nml,
			  float *qual, const KdTree<3> &tree,
			  void (*cb)(float progress), int probes)
{
	NMLArg1 args[NML_NUMTHREAD];
	pthread_t threads[NML_NUMTHREAD];
	size_t payload = point_num / NML_NUMTHREAD;
	for (int i = 0; i < NML_NUMTHREAD; ++i) {
		args[i] = NMLArg1{pos,	  i * payload, (i + 1) * payload,
				  nml,	  qual,	       &tree,
				  probes, NULL};
		if (i == 0)
			args[i].cb = cb;
		if (i == NML_NUMTHREAD - 1)
			args[i].point_stop = point_num;
		pthread_create(&threads[i], NULL, estim_unoriented_nml_threaded,
			       &args[i]);
	}
	for (int i = 0; i < NML_NUMTHREAD; ++i) {
		pthread_join(threads[i], NULL);
	}
}

size_t orient_nml_with_z(Vec3 *nml, EOrient *oriented, size_t point_num,
			 const float *qual, float tol)
{
	size_t unsettled = 0;
	for (size_t i = 0; i < point_num; ++i) {
		if (oriented[i] >= EPositiveZ)
			continue;
		if (fabs(nml[i].z) > tol + 2 * (1 - qual[i])) {
			oriented[i] = EPositiveZ;
			if (nml[i].z < 0)
				nml[i] *= -1;
		} else if (oriented[i] < EOriented) {
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
		} else if (oriented[i] < EOriented) {
			unsettled++;
		}
	}
	return (unsettled);
}

struct NmlOrientArg {
	const Vec3 *pos;
	size_t point_start;
	size_t point_stop;
	const KdTree<3> *tree;
	const float *qual;
	Vec3 *nml;
	EOrient *oriented;
	bool *to_visit;
	int probes;
	float tol;
};

void *propagate_nml_once_threaded(void *arg)
{
	const Vec3 *pos = ((NmlOrientArg *)arg)->pos;
	size_t point_start = ((NmlOrientArg *)arg)->point_start;
	size_t point_stop = ((NmlOrientArg *)arg)->point_stop;
	const KdTree<3> *tree = ((NmlOrientArg *)arg)->tree;
	const float *qual = ((NmlOrientArg *)arg)->qual;
	Vec3 *nml = ((NmlOrientArg *)arg)->nml;
	EOrient *oriented = ((NmlOrientArg *)arg)->oriented;
	bool *to_visit = ((NmlOrientArg *)arg)->to_visit;
	int probes = ((NmlOrientArg *)arg)->probes;
	float tol = ((NmlOrientArg *)arg)->tol;

	TArray<unsigned> knn_idx(probes);
	TArray<float> sqr_dist(probes);
	for (size_t i = point_start; i < point_stop; ++i) {
		if (oriented[i] >= EOriented || oriented[i] == EDead)
			continue;
		const float *query = &pos[i].x;
		tree->knnSearch(query, probes, &knn_idx[0], &sqr_dist[0]);
		int majority = 0;
		int votes = 0;
		for (int j = 0; j < probes; ++j) {
			unsigned k = knn_idx[j];
			if (oriented[k] >= EOriented) {
				float test = dot(nml[k], nml[i]);
				if (fabs(test) >
				    (tol +
				     (1 - tol) * (1 - qual[i] * qual[k]))) {
					majority += test > 0 ? 1 : -1;
					votes++;
				}
			}
		}
		if (majority != 0 && abs(majority) >= votes / 2) {
			oriented[i] = ETmpPlague;
			if (majority < 0)
				nml[i] *= -1;
			/* Mark neigbours visited */
			/* Note : possible race write but should be harmless
			 * because all threads which to do the same flase ->
			 * true update */
			for (int j = 0; j < probes; ++j) {
				unsigned k = knn_idx[j];
				to_visit[k] = true;
			}
		}
	}
	pthread_exit(NULL);
}

size_t propagate_nml_once(const Vec3 *pos, size_t point_num,
			  const KdTree<3> &tree, const float *qual, Vec3 *nml,
			  EOrient *oriented, int probes, float tol)
{
	TArray<EOrient> next_status(point_num, EDead);
	TArray<bool> to_visit(point_num, false);

	NmlOrientArg args[NML_NUMTHREAD];
	pthread_t threads[NML_NUMTHREAD];
	size_t payload = point_num / NML_NUMTHREAD;
	for (int i = 0; i < NML_NUMTHREAD; ++i) {
		args[i] = NmlOrientArg{
		    pos, i * payload, (i + 1) * payload, &tree,	 qual,
		    nml, oriented,    &to_visit[0],	 probes, tol};
		if (i == NML_NUMTHREAD - 1)
			args[i].point_stop = point_num;
		pthread_create(&threads[i], NULL, propagate_nml_once_threaded,
			       &args[i]);
	}
	for (int i = 0; i < NML_NUMTHREAD; ++i) {
		pthread_join(threads[i], NULL);
	}
	size_t newly_settled = 0;
	for (size_t i = 0; i < point_num; ++i) {
		if (oriented[i] < ETmpPlague) {
			oriented[i] = to_visit[i] ? ENone : EDead;
		} else if (oriented[i] == ETmpPlague) {
			oriented[i] = EPlague;
			newly_settled++;
		}
	}
	return (newly_settled);
}

#undef NML_NUMTHREAD
