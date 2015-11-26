#include "MSCKF.h"
#include "RK.h"
#include "MVG.h"

using namespace std;
using namespace Eigen;

// 计算 Givens 旋转用到的系数
// 其中 a 为用来置零的元素， b 为需要置零的元素
// 返回 c 为对应的 cos 分量，s 为 sin 分量
inline void Givens(double a, double b, double &c, double &s) {
    if (abs(b)<1.0e-15) {
        c = copysign(1.0, a);
        s = 0;
    }
    else if (abs(a)<1.0e-15) {
        c = 0;
        s = -copysign(1.0, b);
    }
    else if (abs(b) > abs(a)) {
        double t = a / b;
        double u = copysign(sqrt(1 + t*t), b);
        s = -1.0 / u;
        c = -s*t;
    }
    else {
        double t = b / a;
        double u = copysign(sqrt(1 + t*t), a);
        c = 1.0 / u;
        s = -c*t;
    }
}

MSCKF::MSCKF() {
    // 各种参数的初始化
    m_state_limit = 10;
    Matrix3d R_imu_to_cam;
    R_imu_to_cam << -Vector3d::UnitY(), -Vector3d::UnitX(), -Vector3d::UnitZ();
    m_q_imu_to_cam = HamiltonToJPL(Quaterniond(R_imu_to_cam));
    m_p_cam_in_imu << 0.0065, 0.0638, 0.0;

    m_cov_ng.setIdentity();
    m_cov_nwg.setIdentity();
    m_cov_na.setIdentity();
    m_cov_nwa.setIdentity();
    m_sigma_im_squared = 1.0;
}

void MSCKF::setNoiseCov(const Matrix3d &cov_ng, const Matrix3d &cov_nwg, const Matrix3d &cov_na, const Matrix3d &cov_nwa, double sigma_im_squared) {
    m_cov_ng = cov_ng;
    m_cov_nwg = cov_nwg;
    m_cov_na = cov_na;
    m_cov_nwa = cov_nwa;
    m_sigma_im_squared = sigma_im_squared;
}

void MSCKF::initialize(const JPL_Quaternion &q, const Vector3d &bg, const Vector3d &v, const Vector3d &ba, const Vector3d &p, double g) {
    m_q = q;
    m_bg = bg;
    m_v = v;
    m_ba = ba;
    m_p = p;
    m_g = g;
}

// [1] 中的 (9)(12)和(13)，bg 和 ba 是固定的，所以没有包含在这里，而PII和Phi都需要在相同的时间内进行积分所以包含在了这里
// StateVector := q, v, p, PII, Phi
typedef tuple<Vector4d, Vector3d, Vector3d, Matrix15d, Matrix15d> StateVector;

// StateVector 的数乘，用于 RK 进行积分
inline StateVector operator*(double s, const StateVector &v2) {
    return make_tuple(s*get<0>(v2), s*get<1>(v2), s*get<2>(v2), s*get<3>(v2), s*get<4>(v2));
}

// StateVector 的加法，用于 RK 进行积分
inline void operator+=(StateVector &v, const StateVector &v1) {
    get<0>(v) += get<0>(v1);
    get<1>(v) += get<1>(v1);
    get<2>(v) += get<2>(v1);
    get<3>(v) += get<3>(v1);
    get<4>(v) += get<4>(v1);
}

class MotionSystem {
public:
    MotionSystem(const MSCKF &filter, double t, const Vector3d &w, const Vector3d &a) :
        dt(t - filter.m_t_old),
        w0(filter.m_w_old - filter.m_bg),
        dw(w - filter.m_w_old),
        a0(filter.m_a_old - filter.m_ba),
        da(a - filter.m_a_old),
        filter(filter)
    {}

    StateVector operator() (double t, StateVector s) {
        StateVector r;
        double tt = t - filter.m_t_old;
        Vector3d wt = w0 + dw*tt / dt;
        Vector3d at = a0 + da*tt / dt;

        const JPL_Quaternion& sq = get<0>(s);
        const Vector3d& sv = get<1>(s);
        const Vector3d& sp = get<2>(s);
        const Matrix15d& sPII = get<3>(s);
        const Matrix15d& sPhi = get<4>(s);

        JPL_Quaternion& rq = get<0>(r);
        Vector3d& rv = get<1>(r);
        Vector3d& rp = get<2>(r);
        Matrix15d& rPII = get<3>(r);
        Matrix15d& rPhi = get<4>(r);

        Matrix3d cqt = JPL_CT(sq);
        Matrix3d neg_cross_w = -JPL_Cross(wt);
        Matrix3d cross_a = JPL_Cross(at);

        // (9)
        rq = 0.5*(JPL_Omega(wt)*sq);
        rv = cqt*at; rv.z() += filter.m_g;
        rp = sv;

        // (12)
        Matrix15d FxPII;
        multiplyFM(FxPII, neg_cross_w, -cqt, cross_a, sPII);
        rPII = FxPII + FxPII.transpose();
        rPII.block<3, 3>(0, 0) += filter.m_cov_ng;
        rPII.block<3, 3>(3, 3) += filter.m_cov_nwg;
        rPII.block<3, 3>(6, 6) += filter.m_cov_na;
        rPII.block<3, 3>(9, 9) += filter.m_cov_nwa;

        // (13)
        multiplyFM(rPhi, neg_cross_w, -cqt, cross_a, sPhi);

        return r;
    }
private:
    void multiplyFM(Matrix15d &Result, const Matrix3d &neg_cross_w, const Matrix3d &neg_cqt, const Matrix3d &cross_a, const Matrix15d &M) {
        Result.setZero();
        for (int i = 0; i < 15; i += 3) {
            Result.block<3, 3>(0, i) = neg_cross_w*M.block<3, 3>(0, i) - M.block<3, 3>(3, i);
            Result.block<3, 3>(6, i) = neg_cqt*(cross_a*M.block<3, 3>(0, i) + M.block<3, 3>(9, i));
            Result.block<3, 3>(12, i) = M.block<3, 3>(6, i);
        }
    }

    double dt;
    Vector3d w0, dw;
    Vector3d a0, da;
    const MSCKF &filter;
};

RungeKutta<double, StateVector> rk;

void MSCKF::propagate(double t, const Vector3d &w, const Vector3d &a) {
    static const Matrix15d eye = Matrix15d::Identity();
    Matrix15d Phi;

    if (m_has_old) {
        MotionSystem system(*this, t, w, a);
        tie(m_q, m_v, m_p, m_PII, Phi) = rk.integrate(system, tie(m_q, m_v, m_p, m_PII, eye), m_t_old, t);
        m_q = JPL_Normalize(m_q);
    }

    // eq. after (11)
    for (size_t i = 0; i < m_PIC.size(); ++i) {
        m_PIC[i] = Phi*m_PIC[i];
    }

    m_t_old = t;
    m_w_old = w;
    m_a_old = a;
    m_has_old = true;
}

void MSCKF::track(double t, const unordered_map<size_t, pair<size_t, Vector2d>> &matches) {
    //
    // 首先进行跟踪，找出当前丢失的 track
    //
    unordered_map<size_t, vector<Vector2d>> continued_tracks; // 其中包含了所有跟踪下来的 track
    size_t useful_state_length = 0; // 记录跟踪下来的 track 的最大长度
    for (const auto & f : matches) {
        if (m_tracks.count(f.second.first)>0) { // 这个特征跟踪到了已有的特征
            continued_tracks[f.first].swap(m_tracks.at(f.second.first)); // 将跟踪到的旧 track 放到 continued_tracks 中
        }
        useful_state_length = max(useful_state_length, continued_tracks[f.first].size());
        continued_tracks[f.first].push_back(f.second.second); // 将跟踪到的特征的位置加入到 track 中
    }

    vector<vector<Vector2d>> lost_tracks; // 其中包含了丢失的 track
    for (auto &t : m_tracks) {
        lost_tracks.emplace_back(vector<Vector2d>());
        lost_tracks.back().swap(t.second);
        useful_state_length = max(useful_state_length, lost_tracks.back().size());
    }

    //
    // 去除老旧无用的 states
    //
    if (useful_state_length < m_states.size()) {
        size_t num_to_erase = m_states.size() - useful_state_length;
        m_states.erase(m_states.begin(), m_states.begin() + num_to_erase);
        m_PIC.erase(m_PIC.begin(), m_PIC.begin() + num_to_erase);
        m_PCC.erase(m_PCC.begin(), m_PCC.begin() + num_to_erase);
        for (auto &p : m_PCC) {
            p.erase(p.begin(), p.begin() + num_to_erase);
        }
    }

    //
    // 如果状态数依旧达到上限，在加入新的 state 前我们需要删除一些
    //

    // 将删除过的 state 和 track 记录在这里
    vector<size_t> remaining_states_id;
    unordered_map<size_t, vector<Vector2d>> new_tracks;

    // 将被删除的 track 记录在这里，它还要结合 m_state 进行 update
    unordered_map<size_t, vector<pair<Vector2d, size_t>>> jumping_tracks;
    if (m_states.size() == m_state_limit) {
        for (size_t i = 0; i < m_states.size(); ++i) {
            if (i % 3 != 1) {
                remaining_states_id.push_back(i);
            }
        }
        for (auto& t : continued_tracks) { // 检查每个跟踪上的 track，从中去除相关联的特征组成独立的 track
            size_t tbegin = m_states.size() - (t.second.size() - 1);
            vector<Vector2d>& nt = new_tracks[t.first];
            vector<pair<Vector2d, size_t>>& jt = jumping_tracks[t.first];
            for (size_t i = tbegin; i < m_states.size(); ++i) {
                if (i % 3 != 1) {
                    nt.emplace_back(t.second[i - tbegin]);
                }
                else {
                    jt.emplace_back(t.second[i - tbegin], i);
                }
            }
            nt.emplace_back(t.second.back());
        }
    }

    // 把两类 track 统一，方便用同一段逻辑处理
    vector<vector<pair<Vector2d, size_t>>> track_for_update; // 用于 update 的所有 track
    vector<Vector3d> point_for_update; // 所有用于 update 的三维点
    for (size_t i = 0; i < lost_tracks.size(); ++i) {
        if (lost_tracks[i].size() > 1) {
            Vector3d p = LinearLSTriangulation(lost_tracks[i], m_states);
            Vector3d q = m_states.back().first*p + m_states.back().second;
            if (q.z()>0.1 && q.z() < 100.0) {
                point_for_update.push_back(p);
                track_for_update.emplace_back(vector<pair<Vector2d, size_t>>());
                size_t jstart = m_states.size() - lost_tracks[i].size();
                for (size_t j = 0; j < lost_tracks[i].size(); ++j) {
                    track_for_update.back().emplace_back(lost_tracks[i][j], jstart + j);
                }
            }
        }
    }

    for (auto &t : jumping_tracks) {
        if (t.second.size() > 1) {
            Vector3d p = LinearLSTriangulation(t.second, m_states);
            size_t tid = t.second.back().second;
            Vector3d q = m_states[tid].first*p + m_states[tid].second;
            if (q.z()>0.1 && q.z() < 100.0) {
                point_for_update.push_back(p);
                track_for_update.emplace_back(move(t.second));
            }
        }
    }

    for (size_t i = 0; i < track_for_update.size(); ++i) {
        point_for_update[i] = RefineTriangulation(point_for_update[i], track_for_update[i], m_states);
    }

    //
    // 进行 EKF 更新
    //

    int Hrows = 0;
    int Hcols = (int)m_states.size() * 6;

    // 预先计算好 HX 的大小
    for (size_t i = 0; i < track_for_update.size(); ++i) {
        Hrows += (int)track_for_update[i].size() * 2 - 3;
    }

    if (Hrows > 0) {
        MatrixXd HX(Hrows, Hcols);
        VectorXd ro(Hrows, 1);

        int row_start = 0;
        for (size_t j = 0; j < track_for_update.size(); ++j) {
            int nrow = (int)track_for_update[j].size() * 2;
            MatrixXd HXj(nrow, Hcols);
            MatrixXd Hfj(nrow, 3);
            VectorXd rj(nrow);
            HXj.setZero();
            Hfj.setZero();
            const Vector3d &pj = point_for_update[j];
            for (size_t i = 0; i < track_for_update[j].size(); ++i) {
                const size_t ii = track_for_update[j][i].second;
                const Vector2d &zij = track_for_update[j][i].first;
                const Matrix3d &Ri = m_states[ii].first;
                const Vector3d &Ti = m_states[ii].second;
                Vector3d Xji = Ri*pj + Ti;                                       // eq. after (20)
                Vector2d zij_triangulated(Xji.x() / Xji.z(), Xji.y() / Xji.z()); // eq. after (20)
                MatrixXd Jij(2, 3);                                              // eq. after (23)
                Jij.setIdentity();                                               // eq. after (23)
                Jij.col(2) = -zij_triangulated;                                  // eq. after (23)
                Jij /= Xji.z();                                                  // eq. after (23)
                MatrixXd Hfij = Jij*Ri;                                          // (23)
                HXj.block<2, 3>(i * 2, ii * 6) = Jij*JPL_Cross(Xji);             // (22)
                HXj.block<2, 3>(i * 2, ii * 6 + 3) = -Hfij;                      // (22)
                Hfj.block<2, 3>(i * 2, 0) = Hfij;                                // (23)
                rj.block<2, 1>(i * 2, 0) = zij - zij_triangulated;               // (20)
            }

            // (25)(26)：将 HXj 投影到 Hfj 的零空间，通过 Givens 旋转对 Hfj 进行 QR 分解完成
            for (int col = 0; col < 3; ++col) {
                for (int row = (int)Hfj.rows() - 1; row > col; --row) {
                    if (abs(Hfj(row, col)) > 1e-15) {
                        double c, s;
                        Givens(Hfj(row - 1, col), Hfj(row, col), c, s);
                        for (int k = 0; k < 3; ++k) { // 这里我们只在乎右上三角，左下角可以直接置零或者忽略
                            double a = c*Hfj(row - 1, k) - s*Hfj(row, k);
                            double b = s*Hfj(row - 1, k) + c*Hfj(row, k);
                            Hfj(row - 1, k) = a;
                            Hfj(row, k) = b;
                        }
                        for (int k = 0; k < HXj.cols(); ++k) {
                            double a = c*HXj(row - 1, k) - s*HXj(row, k);
                            double b = s*HXj(row - 1, k) + c*HXj(row, k);
                            HXj(row - 1, k) = a;
                            HXj(row, k) = b;
                        }
                        double a = c*rj(row - 1) - s*rj(row);
                        double b = s*rj(row - 1) + c*rj(row);
                        rj(row - 1) = a;
                        rj(row) = b;
                    }
                }
            }
            HX.block(row_start, 0, nrow - 3, Hcols) = HXj.block(3, 0, nrow - 3, Hcols);
            ro.block(row_start, 0, nrow - 3, 1) = rj.block(3, 0, nrow - 3, 1);
            row_start += nrow - 3;
        }

        // (28)(29)：将 HX 投影到它的 range，同样使用 Givens 旋转进行 QR 分解
        for (int col = 0; col < HX.cols(); ++col) {
            for (int row = (int)HX.rows() - 1; row > col; --row) {
                if (abs(HX(row, col)) > 1e-15) {
                    double c, s;
                    Givens(HX(row - 1, col), HX(row, col), c, s);
                    for (int k = 0; k < HX.cols(); ++k) { // 这里我们同样只关心右上三角，但其余部分需要置零
                        double a = c*HX(row - 1, k) - s*HX(row, k);
                        double b = s*HX(row - 1, k) + c*HX(row, k);
                        HX(row - 1, k) = a;
                        HX(row, k) = b;
                    }
                    double a = c*ro(row - 1) - s*ro(row);
                    double b = s*ro(row - 1) + c*ro(row);
                    ro(row - 1) = a;
                    ro(row) = b;
                }
            }
        }

        // 准备完整的协方差矩阵
        MatrixXd P(m_states.size() * 6 + 15, m_states.size() * 6 + 15);
        P.block<15, 15>(0, 0) = m_PII;
        for (int i = 0; i < m_PIC.size(); ++i) {
            P.block<15, 6>(0, 15 + 6 * i) = m_PIC[i];
            P.block<6, 15>(15 + 6 * i, 0) = m_PIC[i].transpose();
            for (int j = 0; j <= i; ++j) {
                P.block<6, 6>(15 + 6 * j, 15 + 6 * i) = m_PCC[i][j];
                if (i != j) {
                    P.block<6, 6>(15 + 6 * i, 15 + 6 * j) = m_PCC[i][j].transpose();
                }
            }
        }

        int Trows = min(Hrows, Hcols);

        // 真实的 TH 左边 15 行为 0
        MatrixXd TH(Trows, 15 + Hcols);
        TH.block(0, 0, Trows, 15).setZero();
        TH.block(0, 15, Trows, Hcols) = HX.block(0, 0, Trows, Hcols);

        // (31)：计算卡尔曼增益
        MatrixXd PTHT = P*TH.transpose();
        MatrixXd PR = TH*PTHT;
        for (int i = 0; i < PR.rows(); ++i) {
            PR(i, i) += m_sigma_im_squared;
        }
        MatrixXd K = PTHT*PR.inverse();
        // (32)：EKF状态更新
        VectorXd dX = K*ro.block(0, 0, Trows, 1);

        m_q = JPL_Correct(m_q, dX.block<3, 1>(0, 0));
        m_bg += dX.block<3, 1>(3, 0);
        m_v += dX.block<3, 1>(6, 0);
        m_ba += dX.block<3, 1>(9, 0);
        m_p += dX.block<3, 1>(12, 0);

        for (size_t i = 0; i < m_states.size(); ++i) {
            JPL_Quaternion q = HamiltonToJPL(Quaterniond(m_states[i].first));
            q = JPL_Correct(q, dX.block<3, 1>(15 + i * 6, 0));
            Vector3d p = -(m_states[i].first.transpose()*m_states[i].second);
            p += dX.block<3, 1>(18 + i * 6, 0);
            m_states[i].first = JPL_toHamilton(q).toRotationMatrix();
            m_states[i].second = -m_states[i].first*p;
        }

        // (33)：EKF方差更新
        MatrixXd KTH = K*TH;
        for (int i = 0; i < KTH.rows(); ++i) {
            KTH(i, i) -= 1;
        }
        MatrixXd Pnew = KTH*P*KTH.transpose() + (m_sigma_im_squared*K)*K.transpose();

        // 将 P 拆分为内部的表达
        m_PII = Pnew.block<15, 15>(0, 0);
        for (int i = 0; i < m_PIC.size(); ++i) {
            m_PIC[i] = Pnew.block<15, 6>(0, 15 + 6 * i);
            for (int j = 0; j <= i; ++j) {
                m_PCC[i][j] = Pnew.block<6, 6>(15 + 6 * j, 15 + 6 * i);
            }
        }
    }

    //
    // 进行状态扩充
    //

    if (m_states.size() == m_state_limit) {
        m_tracks.swap(new_tracks);
        vector<pair<Matrix3d, Vector3d>> new_states(remaining_states_id.size());
        vector<MatrixXd> new_PIC(remaining_states_id.size());
        vector<vector<MatrixXd>> new_PCC(remaining_states_id.size());
        for (size_t i = 0; i < remaining_states_id.size(); ++i) {
            new_states[i] = m_states[remaining_states_id[i]];
            new_PIC[i] = m_PIC[remaining_states_id[i]];
            new_PCC[i].resize(i+1);
            for (size_t j = 0; j < new_PCC[i].size(); ++j) {
                new_PCC[i][j] = m_PCC[remaining_states_id[i]][remaining_states_id[j]];
            }
        }
        m_states.swap(new_states);
        m_PIC.swap(new_PIC);
        m_PCC.swap(new_PCC);
    }
    else { // 如果我们没从中间删除任何状态，所有track都在continued_tracks中
        m_tracks.swap(continued_tracks);
    }

    // 计算当前相机姿态
    // (14)
    JPL_Quaternion q_world_to_cam = JPL_Multiply(m_q_imu_to_cam, m_q);
    Vector3d imu_to_cam_shift_in_world = JPL_CT(m_q)*m_p_cam_in_imu;
    Vector3d p_cam_in_world = m_p + imu_to_cam_shift_in_world;

    // (16)
    Matrix<double, 6, 15> Jc; 
    Jc.setZero();
    Jc.block<3, 3>(0, 0) = JPL_C(m_q_imu_to_cam);
    Jc.block<3, 3>(3, 0) = JPL_Cross(imu_to_cam_shift_in_world);
    Jc.block<3, 3>(3, 12).setIdentity();
    Matrix<double, 15, 6> JcT = Jc.transpose();
    
    // (15)
    m_PCC.push_back(vector<MatrixXd>());
    for (size_t i = 0; i < m_PIC.size(); ++i) {
        m_PCC.back().push_back((Jc*m_PIC[i]).transpose());
    }
    m_PCC.back().push_back(Jc*m_PII*JcT);
    m_PIC.push_back(m_PII*JcT);

    // 增加新的状态
    Matrix3d R_world_to_cam = JPL_C(q_world_to_cam);
    Vector3d T_world_to_cam = -R_world_to_cam*p_cam_in_world;
    m_states.emplace_back(R_world_to_cam, T_world_to_cam);
}
