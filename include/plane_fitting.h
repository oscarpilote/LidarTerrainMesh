#include "vec3.h"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>

Vec3 fit_plane(const Vec3 *points, size_t point_num)
{
	TVec3<double> bary = TVec3<double>::Zero;
	for (size_t i = 0; i < point_num; ++i) {
		bary.x += points[i].x;
		bary.y += points[i].y;
		bary.z += points[i].z;
	}
	bary /= point_num;

	double xx = 0;
	double xy = 0;
	double xz = 0;
	double yx = 0;
	double yy = 0;
	double yz = 0;
	double zx = 0;
	double zy = 0;
	double zz = 0;

	for (size_t i = 0; i < point_num; ++i) {
		TVec3<double> v;
		v.x = points[i].x - bary.x;
		v.y = points[i].y - bary.y;
		v.z = points[i].z - bary.z;
		xx += v.x * v.x;
		xy += v.x * v.y;
		xz += v.x * v.z;
		yx += v.y * v.x;
		yy += v.y * v.y;
		yz += v.y * v.z;
		zx += v.z * v.x;
		zy += v.z * v.y;
		zz += v.z * v.z;
	}

	TVec3<double> weighted_dir = TVec3<double>::Zero;

	{
		double det_x = yy * zz - yz * yz;
		TVec3<double> axis_dir{det_x, xz * yz - xy * zz,
				       xy * yz - xz * yy};
		double weight = det_x * det_x;
		weighted_dir += axis_dir * weight;
	}

	{
		double det_y = xx * zz - xz * xz;
		TVec3<double> axis_dir{xz * yz - xy * zz, det_y,
				       xy * xz - yz * xx};
		double weight = det_y * det_y;
		if (dot(weighted_dir, axis_dir) < 0) {
			weight = -weight;
		}
		weighted_dir += axis_dir * weight;
	}

	{
		double det_z = xx * yy - xy * xy;
		TVec3<double> axis_dir{xy * yz - xz * yy, xy * xz - yz * xx,
				       det_z};
		double weight = det_z * det_z;
		if (dot(weighted_dir, axis_dir) < 0) {
			weight = -weight;
		}
		weighted_dir += axis_dir * weight;
	}

	TVec3<double> n = normalized(weighted_dir);

	return Vec3{(float)n.x, (float)n.y, (float)n.z};
}

template <class T>
void ComputeCovarianceMatrix(const TVec3<T> *points, size_t point_num,
			     TVec3<T> &barycenter, Eigen::Matrix<T, 3, 3> &m)
{
	barycenter = TVec3<T>::Zero;
	for (size_t i = 0; i < point_num; ++i) {
		barycenter += points[i];
	}
	barycenter /= point_num;

	m.setZero();
	Eigen::Matrix<T, 3, 1> p;
	for (size_t i = 0; i < point_num; ++i) {
		p(0) = points[i].x - barycenter.x;
		p(1) = points[i].y - barycenter.y;
		p(2) = points[i].z - barycenter.z;
		m += p * p.transpose();
	}
	m /= point_num;
}

template <class T>
static inline T quality(T eig_min, T eig_mid, T eig_max, T perp_dist2)
{
	/* Note : both qual1 and qual2 fall into [0,1] */

	/* qual1 : flatness of the ellipsoid */
	T qual1 =
	    (eig_mid == 0)
		? 0
		: std::max<T>(1 - eig_min * eig_max / (eig_mid * eig_mid), 0);

	/* qual2 : target point not an outlier wrt the ellipsoid */
	T qual2 =
	    (eig_max == 0) ? 0 : std::max<T>(1 - 3 * perp_dist2 / eig_max, 0);

	/* TODO Profile */
	/* qual : heuristic based on qual1 and qual2 */
	return qual1 * qual2;
	// return std::max<T>(2 * qual1 * qual2 - 1, 0);
}

template <class T>
T fit_plane_vcg(const TVec3<T> *points, size_t point_num, TVec3<T> &normal)
{
	Eigen::Matrix<T, 3, 3> covMat = Eigen::Matrix<T, 3, 3>::Zero();
	TVec3<T> bary;
	ComputeCovarianceMatrix(points, point_num, bary, covMat);

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, 3, 3>> eig(covMat);
	Eigen::Matrix<T, 3, 1> eval = eig.eigenvalues();
	Eigen::Matrix<T, 3, 3> evec = eig.eigenvectors();
	int idx = 0;
	idx = std::abs(eval[1]) < std::abs(eval[0]) ? 1 : 0;
	idx = std::abs(eval[2]) < std::abs(eval[idx]) ? 2 : idx;
	T eig_min = std::abs(eval[idx]);
	T eig_mid = std::min(std::abs(eval[(idx + 1) % 3]),
			     std::abs(eval[(idx + 2) % 3]));
	T eig_max = std::max(std::abs(eval[(idx + 1) % 3]),
			     std::abs(eval[(idx + 2) % 3]));
	normal[0] = evec(0, idx);
	normal[1] = evec(1, idx);
	normal[2] = evec(2, idx);
	normal = normalized(normal);
	T perp_dist = std::abs<T>(dot(normal, points[0] - bary));
	return quality(eig_min, eig_mid, eig_max, perp_dist * perp_dist);
}

template <class T>
void fit_plane_qual(const TVec3<T> *points, size_t point_num, TVec3<T> &normal,
		    T &qual1, T &qual2)
{
	Eigen::Matrix<T, 3, 3> covMat = Eigen::Matrix<T, 3, 3>::Zero();
	TVec3<T> bary;
	ComputeCovarianceMatrix(points, point_num, bary, covMat);

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, 3, 3>> eig(covMat);
	Eigen::Matrix<T, 3, 1> eval = eig.eigenvalues();
	Eigen::Matrix<T, 3, 3> evec = eig.eigenvectors();
	int idx = 0;
	idx = std::abs(eval[1]) < std::abs(eval[0]) ? 1 : 0;
	idx = std::abs(eval[2]) < std::abs(eval[idx]) ? 2 : idx;
	T eig_min = std::abs(eval[idx]);
	T eig_mid = std::min(std::abs(eval[(idx + 1) % 3]),
			     std::abs(eval[(idx + 2) % 3]));
	T eig_max = std::max(std::abs(eval[(idx + 1) % 3]),
			     std::abs(eval[(idx + 2) % 3]));
	normal[0] = evec(0, idx);
	normal[1] = evec(1, idx);
	normal[2] = evec(2, idx);
	normal = normalized(normal);
	T perp_dist = std::abs(dot(normal, points[0] - bary));
	qual1 = (eig_mid == 0) ? 0 : 1 - eig_min / eig_mid;
	qual2 = eig_max == 0
		    ? 0
		    : std::max<T>(1 - 3 * (perp_dist * perp_dist) / eig_max, 0);
}
