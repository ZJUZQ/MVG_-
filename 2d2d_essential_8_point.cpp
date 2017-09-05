
namespace openMVG {
namespace fundamental {
namespace kernel {

template<typename TMatX, typename TMatA>
inline void EncodeEpipolarEquation(const TMatX &x1, const TMatX &x2, TMatA *A) {
  for (typename TMatX::Index i = 0; i < x1.cols(); ++i) {
    const Vec2 xx1 = x1.col(i);
    const Vec2 xx2 = x2.col(i);
    A->row(i) <<
      xx2(0) * xx1(0),  // 0 represents x coords,
      xx2(0) * xx1(1),  // 1 represents y coords.
      xx2(0),
      xx2(1) * xx1(0),
      xx2(1) * xx1(1),
      xx2(1),
      xx1(0),
      xx1(1),
      1.0;
  }
}

}  // namespace kernel
}  // namespace fundamental
}  // namespace openMVG



namespace openMVG {
namespace essential {
namespace kernel {

/**
 * Eight-point algorithm for solving for the essential matrix from normalized
 * image coordinates of point correspondences.
 * See page 294 in HZ Result 11.1.
 *
 */
	
void EightPointRelativePoseSolver::Solve(const Mat &x1, const Mat &x2, std::vector<Mat3> *Es) {
  assert(2 == x1.rows());
  assert(8 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  MatX9 A(x1.cols(), 9);
  fundamental::kernel::EncodeEpipolarEquation(x1, x2, &A);

  Vec9 e;
  Nullspace(A, e);
  Mat3 E = Map<RMat3>(e.data());

  // Find the closest essential matrix to E in frobenius norm
  // E = UD'VT
  if (x1.cols() > 8) {
    Eigen::JacobiSVD<Mat3> USV(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Vec3 d = USV.singularValues();
    const double a = d[0];
    const double b = d[1];
    d << (a+b)/2., (a+b)/2., 0.0;
    E = USV.matrixU() * d.asDiagonal() * USV.matrixV().transpose();
  }
  Es->push_back(E);
}

}  // namespace kernel
}  // namespace essential
}  // namespace openMVG
